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
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_blas, only : gemv, gemm
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_coulomb, only : get_coulomb_hessian, get_potential_derivative
   use tblite_hessian_dispersion, only : get_dispersion_hessian, &
      & get_dispersion_potential_derivative
   use tblite_hessian_halogen, only : get_halogen_hessian
   use tblite_hessian_hamiltonian, only : get_scaled_overlap, &
      & get_hamiltonian_derivative, get_hamiltonian_hessian
   use tblite_hessian_multipole, only : get_multipole_hessian, &
      & get_multipole_potential_derivative, get_multipole_band_hessian
   use tblite_hessian_repulsion, only : get_repulsion_hessian
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_response_integral, only : get_overlap_derivative, get_multipole_derivative
   use tblite_response_kernel, only : has_dispersion_kernel
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_polarizability, only : check_response_support
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_scf_potential, only : add_pot_to_h1
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_hamiltonian_gradient
   implicit none
   private

   public :: xtb_hessian

contains

subroutine xtb_hessian(ctx, mol, calc, wfn, hessian)
   !***********************************
   !* Analytic Hessian, ordered as hessian(3*(jat-1)+jc, 3*(iat-1)+ic).
   !* Each column is the derivative of the gradient w.r.t. one coordinate
   !* and consists of two parts. The response part contracts the first-order
   !* wavefunction of the coupled-perturbed equations like a gradient, the
   !* explicit part is the second derivative of all terms at fixed density
   !* and is added for the complete matrix at once.
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

   integer :: nat, nao, nsh, ndim, iat, ic, idx, iao, jao
   logical :: multipoles
   real(wp) :: cutoff, sigma(3, 3)
   real(wp), allocatable :: lattr(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: selfenergy(:), dsedcn(:), bmat(:, :), wmat(:, :), h1(:, :, :)
   real(wp), allocatable :: dvshdr(:, :), dqshdr(:, :), n0sh(:), dvao(:), fock0(:, :)
   real(wp), allocatable :: dvat(:, :, :), dvdp(:, :, :, :), dvqp(:, :, :, :)
   real(wp), allocatable :: ddpat(:, :, :), dqpat(:, :, :)
   real(wp), allocatable :: grad0(:, :), dEdcn0(:), grad(:, :), dEdcn(:)
   real(wp), allocatable :: ds(:, :, :), dh0(:, :, :), dp(:, :, :), dw(:, :, :)
   real(wp), allocatable :: zero(:, :, :)
   type(integral_type) :: ints
   type(container_cache) :: ccache, dcache
   type(adjacency_list) :: list
   type(potential_type) :: pot0, pot_zero, pot_r
   type(wavefunction_type) :: wfn_r, wfn_m
   type(error_type), allocatable :: error

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   ndim = 3*nat
   multipoles = allocated(calc%coulomb%aes2)
   hessian(:, :) = 0.0_wp

   call check_response_support(mol, calc, wfn, error)
   if (.not.allocated(error)) call check_hessian_support(calc, error)
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if

   !>-- integrals and potential of the unperturbed system
   cutoff = get_cutoff(calc%bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)

   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call calc%ncoord%get_cn(mol, cn, dcndr, dcndL)
   allocate(selfenergy(nsh), dsedcn(nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=selfenergy, dsedcn=dsedcn)

   call new_integral(ints, nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)
   allocate(bmat(nao, nao))
   call get_scaled_overlap(mol, calc%bas, calc%h0, ints%overlap, cutoff, bmat)

   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)

   !>-- folding the potential is only required to expand it to the orbitals
   call new_potential(pot0, mol, calc%bas, 1)
   call pot0%reset
   call calc%coulomb%get_potential(mol, ccache, wfn, pot0)
   if (has_dispersion_kernel(calc)) then
      call calc%dispersion%get_potential(mol, dcache, wfn, pot0)
   end if
   allocate(h1(nao, nao, 1))
   call add_pot_to_h1(calc%bas, ints, pot0, h1)
   call new_potential(pot_zero, mol, calc%bas, 1)
   call pot_zero%reset

   !>-- derivative of the potential at fixed charges and multipole moments,
   !>   the atom-resolved parts are added to the shell-resolved potential
   allocate(dvshdr(nsh, ndim), dqshdr(nsh, ndim), dvat(nat, 3, nat))
   call get_potential_derivative(mol, calc, ccache, wfn, dvshdr)
   allocate(dvdp(3, nat, 3, nat), dvqp(6, nat, 3, nat), source=0.0_wp)
   if (multipoles) then
      call get_multipole_potential_derivative(calc%coulomb%aes2, mol, wfn%qat(:, 1), &
         & wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), dvat, dvdp, dvqp)
      call add_atomic_potential_derivative(calc, dvat, dvshdr)
   end if
   if (has_dispersion_kernel(calc)) then
      call get_dispersion_potential_derivative(calc%dispersion, mol, wfn%qat(:, 1), &
         & dvat, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call add_atomic_potential_derivative(calc, dvat, dvshdr)
   end if

   allocate(wmat(nao, nao))
   call get_wdensity(wfn, wmat)

   !>-- contribution of the core Hamiltonian to the gradient, which is
   !>   contained in the contraction with the first-order potential
   allocate(grad0(3, nat), dEdcn0(nat), zero(nao, nao, 1), source=0.0_wp)
   call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & dsedcn, pot_zero, wfn%density, zero, dEdcn0, grad0, sigma)

   call new_potential(pot_r, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   allocate(ds(nao, nao, 3), dh0(nao, nao, 3), dp(nao, nao, 1), dw(nao, nao, 1))
   allocate(grad(3, nat), dEdcn(nat), dvao(nao), fock0(nao, nao))
   allocate(n0sh(nsh), ddpat(3, nat, 3), dqpat(6, nat, 3), source=0.0_wp)

   !>-- response part, one coupled-perturbed solution per coordinate
   do iat = 1, nat
      call get_overlap_derivative(mol, calc%bas, lattr, cutoff, iat, ds)
      call get_hamiltonian_derivative(mol, calc%bas, calc%h0, bmat, selfenergy, dsedcn, &
         & dcndr, lattr, cutoff, iat, dh0)
      if (multipoles) then
         call get_multipole_derivative(mol, calc%bas, lattr, cutoff, iat, &
            & wfn%density(:, :, 1), pot0%vdp(:, :, 1), pot0%vqp(:, :, 1), dh0, ddpat, dqpat)
      end if
      do ic = 1, 3
         idx = 3*(iat - 1) + ic

         !>-- derivative of the Hamiltonian at fixed charges and moments, the
         !>   potential enters the Hamiltonian as -1/2 (v_m + v_n) S_mn
         do iao = 1, nao
            dvao(iao) = dvshdr(calc%bas%ao2sh(iao), idx)
         end do
         do iao = 1, nao
            do jao = 1, nao
               fock0(jao, iao) = dh0(jao, iao, ic) &
                  & - 0.5_wp*(pot0%vao(jao, 1) + pot0%vao(iao, 1))*ds(jao, iao, ic) &
                  & - 0.5_wp*(dvao(jao) + dvao(iao))*ints%overlap(jao, iao)
            end do
         end do
         if (multipoles) then
            call add_multipole_fock(calc, ints, dvdp(:, :, ic, iat), dvqp(:, :, ic, iat), &
               & fock0)
         end if

         !>-- derivative of the charges and moments at fixed density
         call get_mulliken_shell_charges(calc%bas, ds(:, :, ic), wfn%density, n0sh, &
            & wfn_m%qsh)
         call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
         dqshdr(:, idx) = wfn_m%qsh(:, 1)

         call xtb_nuclear_response(mol, calc, wfn, ints, ccache, dcache, ds(:, :, ic), &
            & fock0, wfn_m, dp, dw(:, :, 1), pot_r, wfn_r, error)
         if (allocated(error)) then
            call ctx%set_error(error)
            return
         end if

         !>-- first-order density with the unperturbed potential and
         !>   unperturbed density with the first-order potential
         grad(:, :) = -grad0
         dEdcn(:) = -dEdcn0
         call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
            & dsedcn, pot0, dp, dw, dEdcn, grad, sigma)
         call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
            & dsedcn, pot_r, wfn%density, zero, dEdcn, grad, sigma)
         call gemv(dcndr, dEdcn, grad, beta=1.0_wp)

         !>-- first-order charges and moments with the derivative of the potential
         if (multipoles) then
            call add_moment_gradient(wfn_r, dvdp, dvqp, grad)
         end if
         call gemv(dvshdr, wfn_r%qsh(:, 1), hessian(:, idx), trans="t")
         hessian(:, idx) = hessian(:, idx) + reshape(grad, [ndim])
      end do
   end do

   hessian(:, :) = 0.5_wp * (hessian + transpose(hessian))

   !>-- explicit part, second derivatives at fixed density
   call get_hamiltonian_hessian(mol, calc%bas, calc%h0, calc%ncoord, selfenergy, dsedcn, &
      & wfn%density(:, :, 1), wmat, pot0%vao(:, 1), dcndr, lattr, cutoff, hessian)

   !>-- the potential depends on the geometry and acts on the charges, whose
   !>   derivative at fixed density stems from the overlap
   call gemm(dqshdr, dvshdr, hessian, transa="t", alpha=0.5_wp, beta=1.0_wp)
   call gemm(dvshdr, dqshdr, hessian, transa="t", alpha=0.5_wp, beta=1.0_wp)

   if (multipoles) then
      call get_multipole_band_hessian(mol, calc%bas, wfn%density(:, :, 1), &
         & pot0%vdp(:, :, 1), pot0%vqp(:, :, 1), dvdp, dvqp, lattr, cutoff, hessian)
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

subroutine check_hessian_support(calc, error)
   !***********************************
   !* Reject Hamiltonians for which second derivatives are not implemented.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (calc%h0%do_diat_scale) then
      call fatal_error(error, "Analytic Hessian is not available with diatomic scaling")
      return
   end if
   if (.not.allocated(calc%ncoord)) then
      call fatal_error(error, "Analytic Hessian requires a coordination number")
      return
   end if
end subroutine check_hessian_support

subroutine add_atomic_potential_derivative(calc, dvat, dvshdr)
   !***********************************
   !* Add the derivative of an atom-resolved potential to the one of the
   !* shell-resolved potential.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Derivative of the atom-resolved potential
   real(wp), intent(in) :: dvat(:, :, :)
   !> Derivative of the shell-resolved potential
   real(wp), intent(inout) :: dvshdr(:, :)

   integer :: iat, ish, kat, ic

   do kat = 1, size(dvat, 3)
      do ic = 1, 3
         do iat = 1, size(dvat, 1)
            do ish = 1, calc%bas%nsh_at(iat)
               dvshdr(calc%bas%ish_at(iat)+ish, 3*(kat-1)+ic) = &
                  & dvshdr(calc%bas%ish_at(iat)+ish, 3*(kat-1)+ic) + dvat(iat, ic, kat)
            end do
         end do
      end do
   end do
end subroutine add_atomic_potential_derivative

subroutine add_multipole_fock(calc, ints, dvdp, dvqp, fock)
   !***********************************
   !* Add the derivative of the Hamiltonian due to the derivative of the
   !* multipole potentials at unchanged integrals.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Integrals at the current geometry
   type(integral_type), intent(in) :: ints
   !> Derivative of the potential acting on the dipole moments
   real(wp), intent(in) :: dvdp(:, :)
   !> Derivative of the potential acting on the quadrupole moments
   real(wp), intent(in) :: dvqp(:, :)
   !> Derivative of the Hamiltonian
   real(wp), intent(inout) :: fock(:, :)

   integer :: iao, jao, iat, jat
   real(wp) :: vij

   do iao = 1, calc%bas%nao
      iat = calc%bas%ao2at(iao)
      do jao = 1, iao
         jat = calc%bas%ao2at(jao)
         vij = -0.5_wp*( &
            & dot_product(ints%dipole(:, jao, iao), dvdp(:, iat)) &
            & + dot_product(ints%dipole(:, iao, jao), dvdp(:, jat)) &
            & + dot_product(ints%quadrupole(:, jao, iao), dvqp(:, iat)) &
            & + dot_product(ints%quadrupole(:, iao, jao), dvqp(:, jat)))
         fock(jao, iao) = fock(jao, iao) + vij
         if (jao /= iao) fock(iao, jao) = fock(iao, jao) + vij
      end do
   end do
end subroutine add_multipole_fock

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
