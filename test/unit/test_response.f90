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

module test_response
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_type
   use tblite_context_type, only : context_type
   use tblite_external_field, only : electric_field
   use tblite_response_polarizability, only : xtb_polarizability
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_response

   real(wp), parameter :: acc = 1.0e-6_wp
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

contains


!> Collect all exported unit tests
subroutine collect_response(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1-polarizability", test_alpha_gfn1), &
      new_unittest("gfn2-polarizability", test_alpha_gfn2), &
      new_unittest("gfn2-polarizability-smearing", test_alpha_gfn2_smearing), &
      new_unittest("unsupported-interaction", test_unsupported, should_fail=.true.) &
      ]

end subroutine collect_response


subroutine test_polarizability_gen(error, mol, method, etemp)
   !***********************************
   !* Compare the analytic polarizability against central finite differences
   !* of the dipole moment in a homogeneous electric field.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the method, gfn1 or gfn2
   character(len=*), intent(in) :: method
   !> Electronic temperature
   real(wp), intent(in) :: etemp

   type(context_type) :: ctx
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn, wfn0
   class(container_type), allocatable :: cont
   real(wp), parameter :: step = 2.0e-4_wp, thr = 1.0e-5_wp
   real(wp), parameter :: disp(4) = [-2.0_wp, -1.0_wp, 1.0_wp, 2.0_wp]
   real(wp), parameter :: weight(4) = [1.0_wp, -8.0_wp, 8.0_wp, -1.0_wp] / 12.0_wp
   real(wp) :: energy, efield(3), dipole(3), alpha(3, 3), numalpha(3, 3)
   integer :: i, k

   call new_calculator(calc, mol, method, error)
   if (allocated(error)) return
   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, 1, etemp)

   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)
   call xtb_polarizability(ctx, mol, calc, wfn0, alpha)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   !>-- five-point stencil, the leading error of the central difference is
   !>   the second hyperpolarizability and too large for a tight threshold
   efield(:) = 0.0_wp
   cont = electric_field(efield)
   call calc%push_back(cont)
   do i = 1, 3
      numalpha(:, i) = 0.0_wp
      do k = 1, 4
         wfn = wfn0
         efield(:) = 0.0_wp
         efield(i) = disp(k) * step
         call calc%pop(cont)
         cont = electric_field(efield)
         call calc%push_back(cont)
         call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
         dipole(:) = matmul(mol%xyz, wfn%qat(:, 1)) + sum(wfn%dpat(:, :, 1), 2)
         numalpha(:, i) = numalpha(:, i) + weight(k) * dipole / step
      end do
   end do

   if (any(abs(alpha - numalpha) > thr)) then
      call test_failed(error, "Polarizability does not match finite differences")
      print "(3es21.14)", alpha
      print "(a)", "---"
      print "(3es21.14)", numalpha
      print "(a)", "---"
      print "(3es21.14)", alpha - numalpha
   end if

end subroutine test_polarizability_gen

subroutine new_calculator(calc, mol, method, error)
   !> Single-point calculator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the method, gfn1 or gfn2
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(method)
   case("gfn1")
      call new_gfn1_calculator(calc, mol, error)
   case default
      call new_gfn2_calculator(calc, mol, error)
   end select
end subroutine new_calculator

subroutine test_alpha_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "03")
   call test_polarizability_gen(error, mol, "gfn1", kt)

end subroutine test_alpha_gfn1

subroutine test_alpha_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_polarizability_gen(error, mol, "gfn2", kt)

end subroutine test_alpha_gfn2

subroutine test_alpha_gfn2_smearing(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_polarizability_gen(error, mol, "gfn2", 0.02_wp)

end subroutine test_alpha_gfn2_smearing

subroutine test_unsupported(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   class(container_type), allocatable :: cont
   real(wp) :: energy, alpha(3, 3)

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)

   cont = electric_field([0.0_wp, 0.0_wp, 1.0e-3_wp])
   call calc%push_back(cont)
   call xtb_singlepoint(ctx, mol, calc, wfn, acc, energy, verbosity=0)
   call xtb_polarizability(ctx, mol, calc, wfn, alpha)
   if (ctx%failed()) call ctx%get_error(error)

end subroutine test_unsupported

end module test_response
