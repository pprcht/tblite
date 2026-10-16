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

!> @file tblite/hessian/reference.f90
!> Provides the perturbation of the Hamiltonian by a nuclear displacement

!> Derivative of the Hamiltonian, the overlap and the atomic moments w.r.t.
!> the atomic positions at fixed density, charges and multipole moments
module tblite_hessian_reference
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_coulomb, only : get_potential_derivative
   use tblite_hessian_dispersion, only : get_dispersion_potential_derivative
   use tblite_hessian_hamiltonian, only : get_scaled_overlap, get_hamiltonian_derivative
   use tblite_hessian_multipole, only : get_multipole_potential_derivative
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_response_integral, only : get_overlap_derivative, get_multipole_derivative
   use tblite_response_kernel, only : has_dispersion_kernel
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_potential, only : add_pot_to_h1
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian, get_hamiltonian_gradient
   implicit none
   private

   public :: reference_data, new_reference_data, get_nuclear_perturbation
   public :: add_atomic_potential_derivative, add_multipole_fock
   public :: check_hessian_support

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

subroutine new_reference_data(ref, mol, calc, wfn, error)
   !***********************************
   !* Collect the integrals, the converged potential and its derivative at
   !* fixed charges and multipole moments for the unperturbed system.
   !***********************************
   !> Quantities of the unperturbed system
   type(reference_data), intent(out) :: ref
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh
   logical :: multipoles
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: cn(:), dcndL(:, :, :), h1(:, :, :), dvat(:, :, :)
   real(wp), allocatable :: zero(:, :, :)
   type(container_cache) :: ccache, dcache

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)

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
   allocate(ref%dvshdr(nsh, 3*nat), dvat(nat, 3, nat))
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

end subroutine new_reference_data

subroutine get_nuclear_perturbation(mol, calc, wfn, ref, iat, dipole, ds, fock0, dqsh, &
      & ddpat, dqpat)
   !***********************************
   !* Perturbation by the displacement of one atom in the three Cartesian
   !* directions, which is the last dimension of all results: derivative of
   !* the overlap, derivative of the Hamiltonian at fixed charges and
   !* multipole moments, and derivative of the shell charges and atomic
   !* multipole moments at fixed density.
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
   !> Evaluate the derivative of the atomic dipole moments for all methods
   logical, intent(in) :: dipole
   !> Derivative of the overlap matrix
   real(wp), intent(out) :: ds(:, :, :)
   !> Derivative of the Hamiltonian at fixed charges and multipole moments
   real(wp), intent(out) :: fock0(:, :, :)
   !> Derivative of the shell charges at fixed density
   real(wp), intent(out) :: dqsh(:, :)
   !> Derivative of the atomic dipole moments at fixed density
   real(wp), intent(out) :: ddpat(:, :, :)
   !> Derivative of the atomic quadrupole moments at fixed density
   real(wp), intent(out) :: dqpat(:, :, :)

   integer :: ic, idx, iao, jao
   logical :: multipoles
   real(wp), allocatable :: dvao(:), n0sh(:), qsh(:, :)

   multipoles = allocated(calc%coulomb%aes2)
   allocate(dvao(calc%bas%nao), qsh(calc%bas%nsh, 1))
   allocate(n0sh(calc%bas%nsh), source=0.0_wp)
   ddpat(:, :, :) = 0.0_wp
   dqpat(:, :, :) = 0.0_wp

   call get_overlap_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, ds)
   call get_hamiltonian_derivative(mol, calc%bas, calc%h0, ref%bmat, ref%selfenergy, &
      & ref%dsedcn, ref%dcndr, ref%lattr, ref%cutoff, iat, fock0)
   !>-- the multipole potentials vanish if the moments do not enter the
   !>   Hamiltonian, then only the derivative of the atomic dipoles is used
   if (multipoles .or. dipole) then
      call get_multipole_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, &
         & wfn%density(:, :, 1), ref%pot%vdp(:, :, 1), ref%pot%vqp(:, :, 1), fock0, &
         & ddpat, dqpat)
   end if

   do ic = 1, 3
      idx = 3*(iat - 1) + ic

      !>-- the potential enters the Hamiltonian as -1/2 (v_m + v_n) S_mn
      do iao = 1, calc%bas%nao
         dvao(iao) = ref%dvshdr(calc%bas%ao2sh(iao), idx)
      end do
      do iao = 1, calc%bas%nao
         do jao = 1, calc%bas%nao
            fock0(jao, iao, ic) = fock0(jao, iao, ic) &
               & - 0.5_wp*(ref%pot%vao(jao, 1) + ref%pot%vao(iao, 1))*ds(jao, iao, ic) &
               & - 0.5_wp*(dvao(jao) + dvao(iao))*ref%ints%overlap(jao, iao)
         end do
      end do
      if (multipoles) then
         call add_multipole_fock(calc, ref%ints, ref%dvdp(:, :, ic, iat), &
            & ref%dvqp(:, :, ic, iat), fock0(:, :, ic))
      end if

      call get_mulliken_shell_charges(calc%bas, ds(:, :, ic), wfn%density, n0sh, qsh)
      dqsh(:, ic) = qsh(:, 1)
   end do
end subroutine get_nuclear_perturbation

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

end module tblite_hessian_reference
