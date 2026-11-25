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

!> @file tblite/response/polarizability.f90
!> Provides the analytic static dipole polarizability

!> Static dipole polarizability of a self-consistent xTB wavefunction
module tblite_response_polarizability
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_partition, only : serial_work_partition, operator(==)
   use tblite_response_field, only : xtb_field_response
   use tblite_response_type, only : kernel_type
   use tblite_scf, only : potential_type, new_potential
   use tblite_wavefunction_mulliken, only : get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_polarizability, check_response_support

   !> First-order response to a homogeneous electric field in one direction
   type, public :: field_data
      !> First-order density matrix and Hamiltonian
      real(wp), allocatable :: dp(:, :, :), dh(:, :, :)
      !> First-order potential including the field
      type(potential_type) :: pot
      !> First-order charges and multipole moments
      type(wavefunction_type) :: wfn
      !> Derivative of the response kernel along the first-order moments
      type(kernel_type) :: kernel
   end type field_data

contains

subroutine xtb_polarizability(mol, calc, wfn, ints, kernel, field, alpha, error)
   !***********************************
   !* Static dipole polarizability alpha(a, b) = dmu_a/dF_b from one
   !* coupled-perturbed field response per Cartesian direction.
   !* The field responses are returned for the derivative of the
   !* polarizability, their kernel derivatives are not set up here.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Integrals of the single point calculation
   type(integral_type), intent(in) :: ints
   !> Response kernel at the converged state
   type(kernel_type), intent(in) :: kernel
   !> Field response for the three directions
   type(field_data), intent(out) :: field(3)
   !> Polarizability tensor
   real(wp), intent(out) :: alpha(3, 3)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: idir, nat, nao, nsh
   real(wp), allocatable :: ddpat(:, :, :)

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   alpha(:, :) = 0.0_wp
   allocate(ddpat(3, nat, 1))

   do idir = 1, 3
      associate(fld => field(idir))
         allocate(fld%dp(nao, nao, 1), fld%dh(nao, nao, 1))
         call new_potential(fld%pot, mol, calc%bas, 1)
         call new_wavefunction(fld%wfn, nat, nsh, nao, 1, wfn%kt)
         call xtb_field_response(mol, calc, wfn, ints, kernel, idir, fld%dp, fld%dh, &
            & fld%pot, fld%wfn, error)
         if (allocated(error)) return

         !>-- mu = sum_A (R_A q_A + d_A), the atomic dipoles enter for all methods
         call get_mulliken_atomic_multipoles(calc%bas, ints%dipole, fld%dp, ddpat)
         alpha(:, idir) = matmul(mol%xyz, fld%wfn%qat(:, 1)) + sum(ddpat(:, :, 1), dim=2)
      end associate
   end do

end subroutine xtb_polarizability

subroutine check_response_support(mol, calc, wfn, error)
   !***********************************
   !* Reject calculations outside the validity of the analytic response:
   !* open-shell, non-self-consistent, periodic and distributed. Interaction
   !* containers without analytic second derivatives report an error when
   !* setting up the response kernel.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (wfn%nspin /= 1) then
      call fatal_error(error, "Analytic response requires a restricted wavefunction")
      return
   end if
   if (.not.allocated(calc%coulomb)) then
      call fatal_error(error, "Analytic response requires a self-consistent method")
      return
   end if
   if (any(mol%periodic)) then
      call fatal_error(error, "Analytic response is not available for periodic systems")
      return
   end if
   if (.not.(calc%partition == serial_work_partition)) then
      call fatal_error(error, "Analytic response is not available with a work partition")
      return
   end if
end subroutine check_response_support

end module tblite_response_polarizability
