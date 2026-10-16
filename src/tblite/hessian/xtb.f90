! This file is part of tblite.
! SPDX-Identifier: LGPL-3.0-or-later
!
! tblite is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! tblite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with tblite.  If not, see <https://www.gnu.org/licenses/>.

!> @file tblite/hessian/xtb.f90
!> Provides the analytic molecular Hessian

!> Analytic second derivatives of the energy of a self-consistent xTB
!> calculation w.r.t. the atomic positions
module tblite_hessian_xtb
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_blas, only : gemv, gemm
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_hessian_coulomb, only : get_coulomb_hessian
   use tblite_hessian_dispersion, only : get_dispersion_hessian
   use tblite_hessian_halogen, only : get_halogen_hessian
   use tblite_hessian_hamiltonian, only : get_hamiltonian_hessian
   use tblite_hessian_multipole, only : get_multipole_hessian, get_multipole_band_hessian
   use tblite_hessian_reference, only : reference_data, new_reference_data, &
      & get_nuclear_perturbation, check_hessian_support
   use tblite_hessian_repulsion, only : get_repulsion_hessian
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_polarizability, only : check_response_support
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_hamiltonian_gradient
   implicit none
   private

   public :: xtb_hessian

contains

subroutine xtb_hessian(ctx, mol, calc, wfn, hessian, dipgrad)
   !***********************************
   !* Analytic Hessian, ordered as hessian(3*(jat-1)+jc, 3*(iat-1)+ic).
   !* Each column is the derivative of the gradient w.r.t. one coordinate
   !* and consists of two parts. The response part contracts the first-order
   !* wavefunction of the coupled-perturbed equations like a gradient, the
   !* explicit part is the second derivative of all terms at fixed density
   !* and is added for the complete matrix at once.
   !* The derivative of the dipole moment w.r.t. the atomic positions is
   !* available from the same first-order wavefunctions.
   !* Requires the converged wavefunction of a preceding single point
   !* calculation with the same calculator. Failures are reported in ctx.
   !***********************************
   !> Calculation context
   type(context_type), intent(inout) :: ctx
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(out) :: hessian(:, :)
   !> Derivative of the dipole moment, dipgrad(:, ic, iat) for the
   !> displacement of atom iat in direction ic
   real(wp), intent(out), optional :: dipgrad(:, :, :)

   integer :: nat, nao, nsh, iat
   logical :: multipoles, dipole
   real(wp), allocatable :: wmat(:, :), dqshdr(:, :), dmu(:, :, :)
   type(reference_data) :: ref
   type(error_type), allocatable :: error, local_error

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)
   dipole = present(dipgrad)
   hessian(:, :) = 0.0_wp
   if (dipole) dipgrad(:, :, :) = 0.0_wp

   call check_response_support(mol, calc, wfn, error)
   if (.not.allocated(error)) call check_hessian_support(calc, error)
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if

   call new_reference_data(ref, mol, calc, wfn, error)
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if
   allocate(dqshdr(nsh, 3*nat))

   !>-- response part, one coupled-perturbed solution per coordinate
   allocate(dmu(3, 3, nat), source=0.0_wp)
   !$omp parallel do schedule(dynamic) default(none) &
   !$omp shared(mol, calc, wfn, ref, nat, dipole, hessian, dqshdr, dmu, error) &
   !$omp private(iat, local_error)
   do iat = 1, nat
      call get_response_columns(mol, calc, wfn, ref, iat, dipole, &
         & hessian(:, 3*iat-2:3*iat), dqshdr(:, 3*iat-2:3*iat), dmu(:, :, iat), local_error)
      if (allocated(local_error)) then
         !$omp critical (tblite_hessian_xtb_error)
         if (.not.allocated(error)) call move_alloc(local_error, error)
         !$omp end critical (tblite_hessian_xtb_error)
      end if
   end do
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if
   if (dipole) dipgrad(:, :, :) = dmu

   hessian(:, :) = 0.5_wp * (hessian + transpose(hessian))

   !>-- explicit part, second derivatives at fixed density
   allocate(wmat(nao, nao))
   call get_wdensity(wfn, wmat)
   call get_hamiltonian_hessian(mol, calc%bas, calc%h0, calc%ncoord, ref%selfenergy, &
      & ref%dsedcn, wfn%density(:, :, 1), wmat, ref%pot%vao(:, 1), ref%dcndr, ref%lattr, &
      & ref%cutoff, hessian)

   !>-- the potential depends on the geometry and acts on the charges, whose
   !>   derivative at fixed density stems from the overlap
   call gemm(dqshdr, ref%dvshdr, hessian, transa="t", alpha=0.5_wp, beta=1.0_wp)
   call gemm(ref%dvshdr, dqshdr, hessian, transa="t", alpha=0.5_wp, beta=1.0_wp)

   if (multipoles) then
      call get_multipole_band_hessian(mol, calc%bas, wfn%density(:, :, 1), &
         & ref%pot%vdp(:, :, 1), ref%pot%vqp(:, :, 1), ref%dvdp, ref%dvqp, ref%lattr, &
         & ref%cutoff, hessian)
      call get_multipole_hessian(calc%coulomb%aes2, mol, wfn%qat(:, 1), &
         & wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), hessian)
   end if

   call get_coulomb_hessian(calc, mol, wfn, hessian, error)
   if (allocated(calc%repulsion)) then
      call get_repulsion_hessian(calc%repulsion, mol, hessian)
   end if
   if (allocated(calc%halogen) .and. .not.allocated(error)) then
      call get_halogen_hessian(calc%halogen, mol, hessian, error)
   end if
   if (allocated(calc%dispersion) .and. .not.allocated(error)) then
      call get_dispersion_hessian(calc%dispersion, mol, wfn%qat(:, 1), hessian, error)
   end if
   if (allocated(error)) call ctx%set_error(error)

end subroutine xtb_hessian

subroutine get_response_columns(mol, calc, wfn, ref, iat, dipole, columns, dqshdr, dmu, error)
   !***********************************
   !* Response part of the three Hessian columns belonging to the
   !* displacements of one atom, together with the derivative of the shell
   !* charges at fixed density and, on request, of the dipole moment.
   !* All work arrays are local to allow evaluating atoms in parallel.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Quantities of the unperturbed system
   type(reference_data), intent(in) :: ref
   !> Displaced atom
   integer, intent(in) :: iat
   !> Evaluate the derivative of the dipole moment
   logical, intent(in) :: dipole
   !> Response part of the Hessian columns
   real(wp), intent(out) :: columns(:, :)
   !> Derivative of the shell charges at fixed density
   real(wp), intent(out) :: dqshdr(:, :)
   !> Derivative of the dipole moment
   real(wp), intent(inout) :: dmu(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh, ic
   logical :: multipoles
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: ds(:, :, :), fock0(:, :, :), dp(:, :, :), dw(:, :, :)
   real(wp), allocatable :: zero(:, :, :)
   real(wp), allocatable :: ddpat(:, :, :), dqpat(:, :, :), dpat(:, :, :)
   real(wp), allocatable :: grad(:, :), dEdcn(:)
   type(container_cache) :: ccache, dcache
   type(potential_type) :: pot_r
   type(wavefunction_type) :: wfn_r, wfn_m

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)

   !>-- the caches serve as scratch space of the containers
   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)
   call new_potential(pot_r, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   allocate(ds(nao, nao, 3), fock0(nao, nao, 3), dp(nao, nao, 1), dw(nao, nao, 1))
   allocate(grad(3, nat), dEdcn(nat), dpat(3, nat, 1), ddpat(3, nat, 3), dqpat(6, nat, 3))
   allocate(zero(nao, nao, 1), source=0.0_wp)

   call get_nuclear_perturbation(mol, calc, wfn, ref, iat, dipole, ds, fock0, dqshdr, &
      & ddpat, dqpat)

   do ic = 1, 3
      !>-- derivative of the charges and moments at fixed density
      wfn_m%qsh(:, 1) = dqshdr(:, ic)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if

      call xtb_nuclear_response(mol, calc, wfn, ref%ints, ccache, dcache, ds(:, :, ic), &
         & fock0(:, :, ic), wfn_m, dp, dw(:, :, 1), pot_r, wfn_r, error)
      if (allocated(error)) return

      !>-- first-order density with the unperturbed potential and
      !>   unperturbed density with the first-order potential
      grad(:, :) = -ref%gradient
      dEdcn(:) = -ref%dEdcn
      call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
         & ref%selfenergy, ref%dsedcn, ref%pot, dp, dw, dEdcn, grad, sigma)
      call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
         & ref%selfenergy, ref%dsedcn, pot_r, wfn%density, zero, dEdcn, grad, sigma)
      call gemv(ref%dcndr, dEdcn, grad, beta=1.0_wp)

      !>-- first-order charges and moments with the derivative of the potential
      if (multipoles) then
         call add_moment_gradient(wfn_r, ref%dvdp, ref%dvqp, grad)
      end if
      call gemv(ref%dvshdr, wfn_r%qsh(:, 1), columns(:, ic), trans="t")
      columns(:, ic) = columns(:, ic) + reshape(grad, [3*nat])

      !>-- mu = sum_A (R_A q_A + d_A), the atomic dipoles enter for all methods
      if (dipole) then
         call get_mulliken_atomic_multipoles(calc%bas, ref%ints%dipole, dp, dpat)
         dmu(:, ic) = matmul(mol%xyz, wfn_r%qat(:, 1)) &
            & + sum(dpat(:, :, 1) + ddpat(:, :, ic), dim=2)
         dmu(ic, ic) = dmu(ic, ic) + wfn%qat(iat, 1)
      end if
   end do
end subroutine get_response_columns

subroutine add_moment_gradient(wfn_r, dvdp, dvqp, gradient)
   !***********************************
   !* Add the contraction of first-order atomic multipole moments with the
   !* derivative of the multipole potentials to the gradient.
   !***********************************
   !> First-order charges and multipole moments
   type(wavefunction_type), intent(in) :: wfn_r
   !> Derivative of the potential acting on the dipole moments
   real(wp), intent(in) :: dvdp(:, :, :, :)
   !> Derivative of the potential acting on the quadrupole moments
   real(wp), intent(in) :: dvqp(:, :, :, :)
   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: kat, ic

   do kat = 1, size(gradient, 2)
      do ic = 1, 3
         gradient(ic, kat) = gradient(ic, kat) &
            & + sum(wfn_r%dpat(:, :, 1) * dvdp(:, :, ic, kat)) &
            & + sum(wfn_r%qpat(:, :, 1) * dvqp(:, :, ic, kat))
      end do
   end do
end subroutine add_moment_gradient

subroutine get_wdensity(wfn, wmat)
   !***********************************
   !* Energy-weighted density matrix of a restricted wavefunction,
   !* W = sum_p f_p e_p C_p C_p^T with the occupations of both spin channels.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Energy-weighted density matrix
   real(wp), intent(out) :: wmat(:, :)

   integer :: iao
   real(wp), allocatable :: scratch(:, :)

   allocate(scratch(size(wmat, 1), size(wmat, 2)))
   do iao = 1, size(wmat, 2)
      scratch(:, iao) = wfn%coeff(:, iao, 1) &
         & * (wfn%focc(iao, 1) + wfn%focc(iao, 2)) * wfn%emo(iao, 1)
   end do
   call gemm(scratch, wfn%coeff(:, :, 1), wmat, transb="t")
end subroutine get_wdensity

end module tblite_hessian_xtb
