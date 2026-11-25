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

!> @file tblite/xtb/response.f90
!> Provides the analytic response properties of a single point calculation

!> Analytic second derivatives of a converged xTB wavefunction: molecular
!> Hessian, dipole derivatives, static polarizability and its derivative.
!> All properties are evaluated from the integrals, the potential and the
!> interaction caches of the single point calculation.
module tblite_xtb_response
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list
   use tblite_hessian_polarizability, only : xtb_polarizability_gradient
   use tblite_hessian_reference, only : reference_data, new_reference_data, &
      & check_hessian_support
   use tblite_hessian_xtb, only : xtb_hessian
   use tblite_integral_type, only : integral_type
   use tblite_response_density, only : has_integer_occupation
   use tblite_response_kernel, only : interaction_cache, get_response_kernel
   use tblite_response_polarizability, only : xtb_polarizability, check_response_support, &
      & field_data
   use tblite_response_type, only : kernel_type
   use tblite_scf, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_response, check_response, response_accuracy

   !> Accuracy of the wavefunction required for response properties, the
   !> analytic response is not variational in the density
   real(wp), parameter :: response_accuracy = 1.0e-3_wp

contains

subroutine check_response(mol, calc, wfn, nuclear, error)
   !***********************************
   !* Reject calculations for which the requested response properties are
   !* not available before the wavefunction is converged.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Derivatives w.r.t. the atomic positions are requested
   logical, intent(in) :: nuclear
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call check_response_support(mol, calc, wfn, error)
   if (allocated(error)) return
   if (nuclear) call check_hessian_support(calc, error)
end subroutine check_response

subroutine xtb_response(mol, calc, wfn, ints, lattr, list, cutoff, selfenergy, dsedcn, &
      & dcndr, pot, cache, error, hessian, wdensity, dipgrad, alpha, dalpha)
   !***********************************
   !* Evaluate the requested response properties of a converged single
   !* point calculation. One response kernel serves all properties, the
   !* field responses of the polarizability enter its derivative.
   !* Input:  integrals, lattice points, neighbor list, cutoff, self-energies,
   !*         CN derivative, converged potential and caches of the single point
   !* Output: hessian(3*nat, 3*nat) and dipgrad(3, 3, nat), alpha(3, 3) and
   !*         dalpha(3, 3, 3, nat), each one evaluated if present.
   !*         wdensity is required for the Hessian.
   !***********************************
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
   !> Derivative of the coordination numbers w.r.t. the atomic positions,
   !> only referenced for derivatives w.r.t. the atomic positions
   real(wp), allocatable, target, intent(in) :: dcndr(:, :, :)
   !> Converged potential expanded to the orbitals
   type(potential_type), target, intent(in) :: pot
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular Hessian
   real(wp), contiguous, intent(out), optional :: hessian(:, :)
   !> Energy-weighted density matrix
   real(wp), intent(in), optional :: wdensity(:, :, :)
   !> Derivative of the dipole moment w.r.t. the atomic positions
   real(wp), intent(out), optional :: dipgrad(:, :, :)
   !> Static dipole polarizability
   real(wp), intent(out), optional :: alpha(:, :)
   !> Derivative of the polarizability w.r.t. the atomic positions
   real(wp), intent(out), optional :: dalpha(:, :, :, :)

   real(wp) :: alpha_field(3, 3)
   type(kernel_type) :: kernel
   type(field_data) :: field(3)
   type(reference_data) :: ref

   if (present(dalpha) .and. .not.has_integer_occupation(wfn)) then
      call fatal_error(error, "Analytic polarizability derivative requires integer "// &
         & "occupation numbers")
      return
   end if
   if (present(hessian) .and. .not.present(wdensity)) then
      call fatal_error(error, "Analytic Hessian requires the energy-weighted density matrix")
      return
   end if
   if ((present(hessian) .or. present(dalpha)) .and. .not.allocated(dcndr)) then
      call fatal_error(error, "Analytic Hessian requires a coordination number")
      return
   end if

   call get_response_kernel(mol, calc, cache, wfn, kernel, error)
   if (allocated(error)) return

   if (present(alpha) .or. present(dalpha)) then
      call xtb_polarizability(mol, calc, wfn, ints, kernel, field, alpha_field, error)
      if (allocated(error)) return
      if (present(alpha)) alpha(:, :) = alpha_field
   end if

   if (present(hessian) .or. present(dalpha)) then
      call new_reference_data(ref, mol, calc, wfn, ints, lattr, list, cutoff, selfenergy, &
         & dsedcn, dcndr, pot)
   end if

   if (present(dalpha)) then
      call xtb_polarizability_gradient(mol, calc, wfn, ref, cache, kernel, field, dalpha, &
         & error)
      if (allocated(error)) return
   end if

   if (present(hessian)) then
      call xtb_hessian(mol, calc, wfn, ref, cache, kernel, wdensity, hessian, error, &
         & dipgrad)
      if (allocated(error)) return
   end if
end subroutine xtb_response

end module tblite_xtb_response
