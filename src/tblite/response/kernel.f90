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

!> @file tblite/response/kernel.f90
!> Provides the second derivatives of all interactions of a calculator

!> Second derivatives of the interactions of an xTB calculator w.r.t. the
!> charges, multipole moments and atomic positions, collected from the
!> electrostatic, dispersion and additional interaction containers
module tblite_response_kernel
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_response_type, only : kernel_type, new_kernel, response_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: get_response_kernel, get_response_derivative, get_interaction_hessian

   !> Caches of the interaction containers of a calculator, each one is
   !> allocated if the respective container is present
   type, public :: interaction_cache
      !> Cache of the Coulomb interactions
      type(container_cache), allocatable :: coulomb
      !> Cache of the dispersion interaction
      type(container_cache), allocatable :: dispersion
      !> Cache of the additional interaction containers
      type(container_cache), allocatable :: interactions
      !> Cache of the repulsion interaction
      type(container_cache), allocatable :: repulsion
      !> Cache of the halogen bonding correction
      type(container_cache), allocatable :: halogen
   end type interaction_cache

contains

subroutine get_response_kernel(mol, calc, cache, wfn, kernel, error, dwfn)
   !***********************************
   !* Set up the response kernel, the second derivative of the energy
   !* w.r.t. the charges and multipole moments at the state given in wfn.
   !* With dwfn, the derivative of the kernel along this change of the
   !* state is set up instead.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Response kernel
   type(kernel_type), intent(out) :: kernel
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order change of the state
   type(wavefunction_type), intent(in), optional :: dwfn

   call new_kernel(kernel, mol, calc%bas, calc%variable_info())
   if (allocated(calc%coulomb)) then
      call calc%coulomb%get_response_kernel(mol, cache%coulomb, wfn, kernel, error, dwfn)
      if (allocated(error)) return
   end if
   if (allocated(calc%dispersion)) then
      call calc%dispersion%get_response_kernel(mol, cache%dispersion, wfn, kernel, error, &
         & dwfn)
      if (allocated(error)) return
   end if
   if (allocated(calc%interactions)) then
      call calc%interactions%get_response_kernel(mol, cache%interactions, wfn, kernel, &
         & error, dwfn)
      if (allocated(error)) return
   end if
end subroutine get_response_kernel

subroutine get_response_derivative(mol, calc, cache, wfn, response, error, dwfn)
   !***********************************
   !* Evaluate the derivative of the potential w.r.t. the perturbed
   !* coordinates at fixed moments and its contraction with the first-order
   !* moments in response. Previous results in response are overwritten.
   !* With dwfn, the derivative of the kernel applied to dwfn is evaluated.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> First-order quantities of the perturbations
   type(response_type), intent(inout) :: response
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order change of the state
   type(wavefunction_type), intent(in), optional :: dwfn

   call response%reset
   if (allocated(calc%coulomb)) then
      call calc%coulomb%get_response_derivative(mol, cache%coulomb, wfn, response, error, &
         & dwfn)
      if (allocated(error)) return
   end if
   if (allocated(calc%dispersion)) then
      call calc%dispersion%get_response_derivative(mol, cache%dispersion, wfn, response, &
         & error, dwfn)
      if (allocated(error)) return
   end if
   if (allocated(calc%interactions)) then
      call calc%interactions%get_response_derivative(mol, cache%interactions, wfn, &
         & response, error, dwfn)
      if (allocated(error)) return
   end if
end subroutine get_response_derivative

subroutine get_interaction_hessian(mol, calc, cache, wfn, hessian, error)
   !***********************************
   !* Add the second derivative of all interactions w.r.t. the atomic
   !* positions at fixed charges and multipole moments.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), contiguous, intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (allocated(calc%repulsion)) then
      call calc%repulsion%get_hessian(mol, cache%repulsion, wfn, hessian, error)
      if (allocated(error)) return
   end if
   if (allocated(calc%halogen)) then
      call calc%halogen%get_hessian(mol, cache%halogen, wfn, hessian, error)
      if (allocated(error)) return
   end if
   if (allocated(calc%coulomb)) then
      call calc%coulomb%get_hessian(mol, cache%coulomb, wfn, hessian, error)
      if (allocated(error)) return
   end if
   if (allocated(calc%dispersion)) then
      call calc%dispersion%get_hessian(mol, cache%dispersion, wfn, hessian, error)
      if (allocated(error)) return
   end if
   if (allocated(calc%interactions)) then
      call calc%interactions%get_hessian(mol, cache%interactions, wfn, hessian, error)
      if (allocated(error)) return
   end if
end subroutine get_interaction_hessian

end module tblite_response_kernel
