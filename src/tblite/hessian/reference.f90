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
   use tblite_adjlist, only : adjacency_list
   use tblite_hessian_hamiltonian, only : get_scaled_overlap, get_hamiltonian_derivative
   use tblite_integral_type, only : integral_type
   use tblite_response_cpscc, only : get_fock_response
   use tblite_response_integral, only : get_overlap_derivative, get_multipole_derivative
   use tblite_response_kernel, only : interaction_cache, get_response_derivative
   use tblite_response_type, only : response_type, new_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_hamiltonian_gradient
   implicit none
   private

   public :: reference_data, new_reference_data, get_nuclear_perturbation
   public :: add_potential_fock, check_hessian_support

   !> Quantities of the unperturbed system shared by all displacements.
   !> The integrals, the neighbor list and the converged potential belong
   !> to the single point calculation and are only referenced here.
   type :: reference_data
      !> Realspace cutoff of the integrals
      real(wp) :: cutoff
      !> Lattice points within the cutoff
      real(wp), pointer :: lattr(:, :) => null()
      !> Neighbor list
      type(adjacency_list), pointer :: list => null()
      !> Derivative of the coordination numbers w.r.t. the atomic positions
      real(wp), pointer :: dcndr(:, :, :) => null()
      !> Self-energies
      real(wp), pointer :: selfenergy(:) => null()
      !> Derivative of the self-energies w.r.t. the coordination number
      real(wp), pointer :: dsedcn(:) => null()
      !> Integrals
      type(integral_type), pointer :: ints => null()
      !> Converged potential expanded to the orbitals
      type(potential_type), pointer :: pot => null()
      !> Scaled overlap matrix of the core Hamiltonian
      real(wp), allocatable :: bmat(:, :)
      !> Vanishing potential
      type(potential_type) :: pot_zero
      !> Gradient of the core Hamiltonian contracted with the converged density
      real(wp), allocatable :: gradient(:, :), dEdcn(:)
   end type reference_data

contains

subroutine new_reference_data(ref, mol, calc, wfn, ints, lattr, list, cutoff, selfenergy, &
      & dsedcn, dcndr, pot)
   !***********************************
   !* Reference the integrals and the converged potential of a single point
   !* calculation and add the quantities of the unperturbed system that
   !* only the nuclear derivatives require. All targets must outlive ref.
   !***********************************
   !> Quantities of the unperturbed system
   type(reference_data), intent(out) :: ref
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Integrals of the single point calculation
   type(integral_type), target, intent(in) :: ints
   !> Lattice points within the cutoff
   real(wp), target, intent(in) :: lattr(:, :)
   !> Neighbor list
   type(adjacency_list), target, intent(in) :: list
   !> Realspace cutoff of the integrals
   real(wp), intent(in) :: cutoff
   !> Self-energies
   real(wp), target, intent(in) :: selfenergy(:)
   !> Derivative of the self-energies w.r.t. the coordination number
   real(wp), target, intent(in) :: dsedcn(:)
   !> Derivative of the coordination numbers w.r.t. the atomic positions
   real(wp), target, intent(in) :: dcndr(:, :, :)
   !> Converged potential expanded to the orbitals
   type(potential_type), target, intent(in) :: pot

   integer :: nat, nao
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: zero(:, :, :)

   nat = mol%nat
   nao = calc%bas%nao

   ref%cutoff = cutoff
   ref%lattr => lattr
   ref%list => list
   ref%selfenergy => selfenergy
   ref%dsedcn => dsedcn
   ref%dcndr => dcndr
   ref%ints => ints
   ref%pot => pot

   allocate(ref%bmat(nao, nao))
   call get_scaled_overlap(mol, calc%bas, calc%h0, ints%overlap, cutoff, ref%bmat)
   call new_potential(ref%pot_zero, mol, calc%bas, 1)
   call ref%pot_zero%reset

   !>-- contribution of the core Hamiltonian to the gradient, which is
   !>   contained in the contraction with the first-order potential
   allocate(ref%gradient(3, nat), ref%dEdcn(nat), zero(nao, nao, 1), source=0.0_wp)
   call get_hamiltonian_gradient(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & dsedcn, ref%pot_zero, wfn%density, zero, ref%dEdcn, ref%gradient, sigma)

end subroutine new_reference_data

subroutine get_nuclear_perturbation(mol, calc, wfn, ref, cache, iat, dipole, ds, fock0, &
      & dqsh, ddpat, dqpat, error)
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
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
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
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, iao, jao
   logical :: multipoles
   real(wp), allocatable :: n0sh(:), qsh(:, :)
   type(response_type) :: response

   multipoles = allocated(calc%coulomb%aes2)
   allocate(qsh(calc%bas%nsh, 1))
   allocate(n0sh(calc%bas%nsh), source=0.0_wp)
   ddpat(:, :, :) = 0.0_wp
   dqpat(:, :, :) = 0.0_wp

   !>-- derivative of the potential at fixed charges and multipole moments
   call new_response(response, mol, calc%bas, [3*iat-2, 3*iat-1, 3*iat])
   call get_response_derivative(mol, calc, cache, wfn, response, error)
   if (allocated(error)) return

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
      !>-- the potential enters the Hamiltonian as -1/2 (v_m + v_n) S_mn
      do iao = 1, calc%bas%nao
         do jao = 1, calc%bas%nao
            fock0(jao, iao, ic) = fock0(jao, iao, ic) &
               & - 0.5_wp*(ref%pot%vao(jao, 1) + ref%pot%vao(iao, 1))*ds(jao, iao, ic)
         end do
      end do
      call add_potential_fock(mol, calc, ref%ints, response, ic, fock0(:, :, ic))

      call get_mulliken_shell_charges(calc%bas, ds(:, :, ic), wfn%density, n0sh, qsh)
      dqsh(:, ic) = qsh(:, 1)
   end do
end subroutine get_nuclear_perturbation

subroutine add_potential_fock(mol, calc, ints, response, ip, fock)
   !***********************************
   !* Add the change of the Hamiltonian due to the derivative of the
   !* potential w.r.t. the perturbed coordinate ip at unchanged integrals.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Integrals at the current geometry
   type(integral_type), intent(in) :: ints
   !> Derivative of the potential w.r.t. the perturbed coordinates
   type(response_type), intent(in) :: response
   !> Perturbed coordinate
   integer, intent(in) :: ip
   !> Perturbation of the Hamiltonian
   real(wp), intent(inout) :: fock(:, :)

   real(wp), allocatable :: h1(:, :, :)
   type(potential_type) :: dpot

   call new_potential(dpot, mol, calc%bas, 1)
   call dpot%reset
   dpot%vat(:, 1) = response%vat(:, ip)
   dpot%vsh(:, 1) = response%vsh(:, ip)
   dpot%vdp(:, :, 1) = response%vdp(:, :, ip)
   dpot%vqp(:, :, 1) = response%vqp(:, :, ip)
   allocate(h1(calc%bas%nao, calc%bas%nao, 1))
   call get_fock_response(calc, ints, dpot, h1)
   fock(:, :) = fock + h1(:, :, 1)
end subroutine add_potential_fock

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
