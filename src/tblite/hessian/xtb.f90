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
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_hamiltonian_gradient
   implicit none
   private

   public :: xtb_hessian

   !> Quantities of the unperturbed system shared by all displacements
   type :: reference_data
      !> Realspace cutoff of the integrals
      real(wp) :: cutoff
      !> Lattice points within the cutoff
      real(wp), allocatable :: lattr(:, :)
      !> Neighbor list
      type(adjacency_list) :: list
      !> Derivative of the coordination numbers w.r.t. the atomic positions
      real(wp), allocatable :: dcndr(:, :, :)
      !> Self-energies and their derivative w.r.t. the coordination number
      real(wp), allocatable :: selfenergy(:), dsedcn(:)
      !> Integrals
      type(integral_type) :: ints
      !> Scaled overlap matrix of the core Hamiltonian
      real(wp), allocatable :: bmat(:, :)
      !> Converged potential and a vanishing one
      type(potential_type) :: pot, pot_zero
      !> Derivative of the shell-resolved potential at fixed charges and moments
      real(wp), allocatable :: dvshdr(:, :)
      !> Derivative of the multipole potentials at fixed charges and moments
      real(wp), allocatable :: dvdp(:, :, :, :), dvqp(:, :, :, :)
      !> Gradient of the core Hamiltonian contracted with the converged density
      real(wp), allocatable :: gradient(:, :), dEdcn(:)
   end type reference_data

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
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: cn(:), dcndL(:, :, :), wmat(:, :), h1(:, :, :)
   real(wp), allocatable :: dqshdr(:, :), dvat(:, :, :), zero(:, :, :), dmu(:, :, :)
   type(reference_data) :: ref
   type(container_cache) :: ccache, dcache
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

   !>-- integrals and potential of the unperturbed system
   ref%cutoff = get_cutoff(calc%bas)
   call get_lattice_points(mol%periodic, mol%lattice, ref%cutoff, ref%lattr)
   call new_adjacency_list(ref%list, mol, ref%lattr, ref%cutoff)

   allocate(cn(nat), ref%dcndr(3, nat, nat), dcndL(3, 3, nat))
   call calc%ncoord%get_cn(mol, cn, ref%dcndr, dcndL)
   allocate(ref%selfenergy(nsh), ref%dsedcn(nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=ref%selfenergy, dsedcn=ref%dsedcn)

   call new_integral(ref%ints, nao)
   call get_hamiltonian(mol, ref%lattr, ref%list, calc%bas, calc%h0, ref%selfenergy, &
      & ref%ints%overlap, ref%ints%dipole, ref%ints%quadrupole, ref%ints%hamiltonian)
   allocate(ref%bmat(nao, nao))
   call get_scaled_overlap(mol, calc%bas, calc%h0, ref%ints%overlap, ref%cutoff, ref%bmat)

   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)

   !>-- folding the potential is only required to expand it to the orbitals
   call new_potential(ref%pot, mol, calc%bas, 1)
   call ref%pot%reset
   call calc%coulomb%get_potential(mol, ccache, wfn, ref%pot)
   if (has_dispersion_kernel(calc)) then
      call calc%dispersion%get_potential(mol, dcache, wfn, ref%pot)
   end if
   allocate(h1(nao, nao, 1))
   call add_pot_to_h1(calc%bas, ref%ints, ref%pot, h1)
   call new_potential(ref%pot_zero, mol, calc%bas, 1)
   call ref%pot_zero%reset

   !>-- derivative of the potential at fixed charges and multipole moments,
   !>   the atom-resolved parts are added to the shell-resolved potential
   allocate(ref%dvshdr(nsh, 3*nat), dqshdr(nsh, 3*nat), dvat(nat, 3, nat))
   call get_potential_derivative(mol, calc, ccache, wfn, ref%dvshdr)
   allocate(ref%dvdp(3, nat, 3, nat), ref%dvqp(6, nat, 3, nat), source=0.0_wp)
   if (multipoles) then
      call get_multipole_potential_derivative(calc%coulomb%aes2, mol, wfn%qat(:, 1), &
         & wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), dvat, ref%dvdp, ref%dvqp)
      call add_atomic_potential_derivative(calc, dvat, ref%dvshdr)
   end if
   if (has_dispersion_kernel(calc)) then
      call get_dispersion_potential_derivative(calc%dispersion, mol, wfn%qat(:, 1), &
         & dvat, error)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if
      call add_atomic_potential_derivative(calc, dvat, ref%dvshdr)
   end if

   !>-- contribution of the core Hamiltonian to the gradient, which is
   !>   contained in the contraction with the first-order potential
   allocate(ref%gradient(3, nat), ref%dEdcn(nat), zero(nao, nao, 1), source=0.0_wp)
   call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
      & ref%selfenergy, ref%dsedcn, ref%pot_zero, wfn%density, zero, ref%dEdcn, &
      & ref%gradient, sigma)

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

   integer :: nat, nao, nsh, ic, idx, iao, jao
   logical :: multipoles
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: ds(:, :, :), dh0(:, :, :), dp(:, :, :), dw(:, :, :)
   real(wp), allocatable :: zero(:, :, :), fock0(:, :), dvao(:), n0sh(:)
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
   allocate(ds(nao, nao, 3), dh0(nao, nao, 3), dp(nao, nao, 1), dw(nao, nao, 1))
   allocate(grad(3, nat), dEdcn(nat), dvao(nao), fock0(nao, nao), dpat(3, nat, 1))
   allocate(zero(nao, nao, 1), n0sh(nsh), ddpat(3, nat, 3), dqpat(6, nat, 3), &
      & source=0.0_wp)

   call get_overlap_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, ds)
   call get_hamiltonian_derivative(mol, calc%bas, calc%h0, ref%bmat, ref%selfenergy, &
      & ref%dsedcn, ref%dcndr, ref%lattr, ref%cutoff, iat, dh0)
   !>-- the multipole potentials vanish if the moments do not enter the
   !>   Hamiltonian, then only the derivative of the atomic dipoles is used
   if (multipoles .or. dipole) then
      call get_multipole_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, &
         & wfn%density(:, :, 1), ref%pot%vdp(:, :, 1), ref%pot%vqp(:, :, 1), dh0, &
         & ddpat, dqpat)
   end if

   do ic = 1, 3
      idx = 3*(iat - 1) + ic

      !>-- derivative of the Hamiltonian at fixed charges and moments, the
      !>   potential enters the Hamiltonian as -1/2 (v_m + v_n) S_mn
      do iao = 1, nao
         dvao(iao) = ref%dvshdr(calc%bas%ao2sh(iao), idx)
      end do
      do iao = 1, nao
         do jao = 1, nao
            fock0(jao, iao) = dh0(jao, iao, ic) &
               & - 0.5_wp*(ref%pot%vao(jao, 1) + ref%pot%vao(iao, 1))*ds(jao, iao, ic) &
               & - 0.5_wp*(dvao(jao) + dvao(iao))*ref%ints%overlap(jao, iao)
         end do
      end do
      if (multipoles) then
         call add_multipole_fock(calc, ref%ints, ref%dvdp(:, :, ic, iat), &
            & ref%dvqp(:, :, ic, iat), fock0)
      end if

      !>-- derivative of the charges and moments at fixed density
      call get_mulliken_shell_charges(calc%bas, ds(:, :, ic), wfn%density, n0sh, &
         & wfn_m%qsh)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if
      dqshdr(:, ic) = wfn_m%qsh(:, 1)

      call xtb_nuclear_response(mol, calc, wfn, ref%ints, ccache, dcache, ds(:, :, ic), &
         & fock0, wfn_m, dp, dw(:, :, 1), pot_r, wfn_r, error)
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
