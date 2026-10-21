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

!> @file tblite/hessian/solvation.f90
!> Provides the implicit solvation models entering the analytic second derivatives

!> Implicit solvation in the analytic second derivatives
module tblite_hessian_solvation
   use mctc_env, only : error_type, fatal_error
   use tblite_container, only : container_type
   use tblite_partition, only : serial_work_partition, operator(==)
   use tblite_solvation, only : alpb_solvation, cds_solvation, shift_solvation, born_kernel
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: check_solvation_support

contains

subroutine check_solvation_support(calc, error)
   !***********************************
   !* Reject interaction containers without analytic second derivatives.
   !* Supported are the generalized Born models (ALPB and GBSA with the P16
   !* or Still kernel), the CDS surface term and the free energy shift, all
   !* evaluated over the complete work of a serial calculation.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic
   logical :: has_cm5, use_cm5
   class(container_type), allocatable :: cont

   if (.not.allocated(calc%interactions)) return
   if (.not.(calc%interactions%partition == serial_work_partition)) then
      call fatal_error(error, "Analytic second derivatives with solvation require "// &
         & "a serial calculation")
      return
   end if

   !>-- the Born and surface terms act on the same CM5-shifted charges
   has_cm5 = .false.
   use_cm5 = .false.
   do ic = 1, calc%interactions%get_size()
      call calc%interactions%get_container(ic, cont)
      if (.not.(cont%partition == serial_work_partition)) then
         call fatal_error(error, "Analytic second derivatives with solvation require "// &
            & "a serial calculation")
         return
      end if
      select type(cont)
      type is(alpb_solvation)
         if (cont%kernel /= born_kernel%p16 .and. cont%kernel /= born_kernel%still) then
            call fatal_error(error, "Analytic second derivatives are not available "// &
               & "for this generalized Born kernel")
            return
         end if
         call check_cm5(cont%useCM5)
      type is(cds_solvation)
         call check_cm5(cont%useCM5)
      type is(shift_solvation)
         continue
      class default
         call fatal_error(error, "Analytic second derivatives are not available for "// &
            & "the interaction '"//container_label(cont)//"'")
         return
      end select
      if (allocated(error)) return
   end do

contains

   subroutine check_cm5(cm5)
      !> The container acts on CM5-shifted charges
      logical, intent(in) :: cm5

      if (has_cm5 .and. (cm5 .neqv. use_cm5)) then
         call fatal_error(error, "Analytic second derivatives require the Born and "// &
            & "surface terms to agree on the CM5 charges")
      end if
      has_cm5 = .true.
      use_cm5 = cm5
   end subroutine check_cm5

end subroutine check_solvation_support

function container_label(cont) result(label)
   !> Interaction container
   class(container_type), intent(in) :: cont
   !> Label of the container, or a placeholder
   character(len=:), allocatable :: label

   if (allocated(cont%label)) then
      label = cont%label
   else
      label = "unknown"
   end if
end function container_label

end module tblite_hessian_solvation
