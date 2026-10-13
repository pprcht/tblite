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

module test_hessian
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io, only : structure_type
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_hessian_dispersion, only : get_dispersion_hessian
   use tblite_hessian_repulsion, only : get_repulsion_hessian
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   implicit none
   private

   public :: collect_hessian

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: kt = 300.0_wp * 3.166808578545117e-06_wp

   abstract interface
      !> Evaluate the gradient of an energy contribution
      subroutine gradient_interface(mol, gradient)
         import :: structure_type, wp
         !> Molecular structure data
         type(structure_type), intent(in) :: mol
         !> Molecular gradient
         real(wp), intent(out) :: gradient(:, :)
      end subroutine gradient_interface
   end interface

contains


!> Collect all exported unit tests
subroutine collect_hessian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("repulsion-gfn1", test_repulsion_gfn1), &
      new_unittest("repulsion-gfn2", test_repulsion_gfn2), &
      new_unittest("dispersion-d3", test_dispersion_d3), &
      new_unittest("dispersion-d4", test_dispersion_d4) &
      ]

end subroutine collect_hessian


subroutine get_numerical_hessian(mol, get_gradient, hessian)
   !***********************************
   !* Central finite differences of an analytic gradient.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Gradient evaluator
   procedure(gradient_interface) :: get_gradient
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(out) :: hessian(:, :)

   integer :: iat, ic
   type(structure_type) :: displ
   real(wp), allocatable :: gr(:, :), gl(:, :)
   real(wp), parameter :: step = 1.0e-4_wp

   allocate(gr(3, mol%nat), gl(3, mol%nat))
   displ = mol
   do iat = 1, mol%nat
      do ic = 1, 3
         displ%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_gradient(displ, gr)
         displ%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call get_gradient(displ, gl)
         displ%xyz(ic, iat) = mol%xyz(ic, iat)
         hessian(:, 3*(iat - 1) + ic) = reshape(gr - gl, [3*mol%nat]) / (2*step)
      end do
   end do
end subroutine get_numerical_hessian

subroutine compare_hessian(error, hessian, numhess)
   !***********************************
   !* Check an analytic Hessian for symmetry and against a numerical one.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Analytic Hessian
   real(wp), intent(in) :: hessian(:, :)
   !> Numerical Hessian
   real(wp), intent(in) :: numhess(:, :)

   if (any(abs(hessian - transpose(hessian)) > thr)) then
      call test_failed(error, "Hessian is not symmetric")
      print "(a, es12.4)", "max. asymmetry:", maxval(abs(hessian - transpose(hessian)))
      return
   end if
   if (any(abs(hessian - numhess) > thr)) then
      call test_failed(error, "Hessian does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(hessian - numhess))
      print "(a, es12.4)", "max. element:  ", maxval(abs(numhess))
   end if
end subroutine compare_hessian

subroutine test_repulsion_gen(error, mol, calc)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc

   real(wp), allocatable :: hessian(:, :), numhess(:, :)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat))
   hessian(:, :) = 0.0_wp
   call get_repulsion_hessian(calc%repulsion, mol, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(container_cache) :: cache
      real(wp) :: energies(displ%nat), sigma(3, 3)

      energies(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      call calc%repulsion%update(displ, cache)
      call calc%repulsion%get_engrad(displ, cache, energies, gradient, sigma)
   end subroutine get_gradient

end subroutine test_repulsion_gen

subroutine test_repulsion_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "01")
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_repulsion_gen(error, mol, calc)

end subroutine test_repulsion_gfn1

subroutine test_repulsion_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "02")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_repulsion_gen(error, mol, calc)

end subroutine test_repulsion_gfn2

subroutine test_dispersion_gen(error, mol, calc)
   !***********************************
   !* The self-consistent part of the dispersion energy is evaluated for a
   !* fixed set of atomic charges.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc

   integer :: iat
   type(wavefunction_type) :: wfn
   real(wp), allocatable :: hessian(:, :), numhess(:, :)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   do iat = 1, mol%nat
      wfn%qat(iat, 1) = 0.3_wp * sin(real(iat, wp))
   end do

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat))
   hessian(:, :) = 0.0_wp
   call get_dispersion_hessian(calc%dispersion, mol, wfn%qat(:, 1), hessian, error)
   if (allocated(error)) return
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(container_cache) :: cache
      real(wp) :: energies(displ%nat), sigma(3, 3)

      energies(:) = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      call calc%dispersion%update(displ, cache)
      call calc%dispersion%get_engrad(displ, cache, energies, gradient, sigma)
      call calc%dispersion%get_gradient(displ, cache, wfn, gradient, sigma)
   end subroutine get_gradient

end subroutine test_dispersion_gen

subroutine test_dispersion_d3(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "03")
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_dispersion_gen(error, mol, calc)

end subroutine test_dispersion_d3

subroutine test_dispersion_d4(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "04")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_dispersion_gen(error, mol, calc)

end subroutine test_dispersion_d4

end module test_hessian
