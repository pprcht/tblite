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
   use mctc_io, only : structure_type, new
   use mstore, only : get_structure
   use tblite_container, only : container_cache
   use tblite_context_type, only : context_type
   use tblite_hessian_dispersion, only : get_dispersion_hessian, &
      & get_dispersion_potential_derivative
   use tblite_hessian_multipole, only : get_multipole_hessian, &
      & get_multipole_potential_derivative
   use tblite_hessian_repulsion, only : get_repulsion_hessian
   use tblite_hessian_xtb, only : xtb_hessian
   use tblite_scf, only : potential_type, new_potential
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_singlepoint, only : xtb_singlepoint
   implicit none
   private

   public :: collect_hessian

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: acc = 1.0e-6_wp
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
      new_unittest("dispersion-d4", test_dispersion_d4), &
      new_unittest("dispersion-potential", test_dispersion_potential), &
      new_unittest("multipole-hessian", test_multipole_hessian), &
      new_unittest("multipole-potential", test_multipole_potential), &
      new_unittest("gfn1-hessian", test_hessian_gfn1), &
      new_unittest("gfn2-hessian", test_hessian_gfn2), &
      new_unittest("gfn1-hessian-smearing", test_hessian_gfn1_smearing), &
      new_unittest("gfn1-hessian-halogen", test_hessian_gfn1_halogen), &
      new_unittest("gfn1-dipgrad", test_dipgrad_gfn1), &
      new_unittest("gfn2-dipgrad", test_dipgrad_gfn2), &
      new_unittest("unsupported-spin", test_unsupported_spin, should_fail=.true.) &
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

subroutine test_dispersion_potential(error)
   !***********************************
   !* Compare the derivative of the dispersion potential at fixed charges
   !* against finite differences of the potential.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   type(structure_type) :: mol, displ
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(potential_type) :: potr, potl
   real(wp), allocatable :: dvat(:, :, :), nvat(:, :, :)
   real(wp), parameter :: step = 1.0e-4_wp

   call get_structure(mol, "MB16-43", "04")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call set_moments(wfn)
   call new_potential(potr, mol, calc%bas, 1)
   call new_potential(potl, mol, calc%bas, 1)

   allocate(dvat(mol%nat, 3, mol%nat), nvat(mol%nat, 3, mol%nat))
   call get_dispersion_potential_derivative(calc%dispersion, mol, wfn%qat(:, 1), dvat, error)
   if (allocated(error)) return

   displ = mol
   do iat = 1, mol%nat
      do ic = 1, 3
         displ%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_potential(displ, potr)
         displ%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call get_potential(displ, potl)
         displ%xyz(ic, iat) = mol%xyz(ic, iat)
         nvat(:, ic, iat) = (potr%vat(:, 1) - potl%vat(:, 1)) / (2*step)
      end do
   end do

   if (any(abs(dvat - nvat) > thr)) then
      call test_failed(error, "Potential derivative does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(dvat - nvat))
      print "(a, es12.4)", "max. element:  ", maxval(abs(nvat))
   end if

contains

   subroutine get_potential(geo, pot)
      type(structure_type), intent(in) :: geo
      type(potential_type), intent(inout) :: pot

      type(container_cache) :: cache

      call pot%reset
      call calc%dispersion%update(geo, cache)
      call calc%dispersion%get_potential(geo, cache, wfn, pot)
   end subroutine get_potential

end subroutine test_dispersion_potential

subroutine set_moments(wfn)
   !***********************************
   !* Fill a wavefunction with an arbitrary fixed set of atomic charges,
   !* dipole and traceless quadrupole moments.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(inout) :: wfn

   integer :: iat, k

   do iat = 1, size(wfn%qat, 1)
      wfn%qat(iat, 1) = 0.3_wp * sin(real(iat, wp))
      do k = 1, 3
         wfn%dpat(k, iat, 1) = 0.2_wp * cos(real(3*iat + k, wp))
      end do
      do k = 1, 6
         wfn%qpat(k, iat, 1) = 0.4_wp * sin(real(7*iat + 2*k, wp))
      end do
      wfn%qpat(6, iat, 1) = -wfn%qpat(1, iat, 1) - wfn%qpat(3, iat, 1)
   end do
end subroutine set_moments

subroutine test_multipole_hessian(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp), allocatable :: hessian(:, :), numhess(:, :)

   call get_structure(mol, "MB16-43", "06")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call set_moments(wfn)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat))
   hessian(:, :) = 0.0_wp
   call get_multipole_hessian(calc%coulomb%aes2, mol, wfn%qat(:, 1), wfn%dpat(:, :, 1), &
      & wfn%qpat(:, :, 1), hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(container_cache) :: cache
      real(wp) :: sigma(3, 3)

      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      call calc%coulomb%aes2%update(displ, cache)
      call calc%coulomb%aes2%get_gradient(displ, cache, wfn, gradient, sigma)
   end subroutine get_gradient

end subroutine test_multipole_hessian

subroutine test_multipole_potential(error)
   !***********************************
   !* Compare the derivative of the multipole potentials at fixed moments
   !* against finite differences of the potentials.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   type(structure_type) :: mol, displ
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   type(potential_type) :: potr, potl
   real(wp), allocatable :: dvat(:, :, :), dvdp(:, :, :, :), dvqp(:, :, :, :)
   real(wp), allocatable :: nvat(:, :, :), nvdp(:, :, :, :), nvqp(:, :, :, :)
   real(wp), parameter :: step = 1.0e-4_wp

   call get_structure(mol, "MB16-43", "06")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call set_moments(wfn)
   call new_potential(potr, mol, calc%bas, 1)
   call new_potential(potl, mol, calc%bas, 1)

   allocate(dvat(mol%nat, 3, mol%nat), dvdp(3, mol%nat, 3, mol%nat), &
      & dvqp(6, mol%nat, 3, mol%nat), nvat(mol%nat, 3, mol%nat), &
      & nvdp(3, mol%nat, 3, mol%nat), nvqp(6, mol%nat, 3, mol%nat))
   call get_multipole_potential_derivative(calc%coulomb%aes2, mol, wfn%qat(:, 1), &
      & wfn%dpat(:, :, 1), wfn%qpat(:, :, 1), dvat, dvdp, dvqp)

   displ = mol
   do iat = 1, mol%nat
      do ic = 1, 3
         displ%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_potential(displ, potr)
         displ%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call get_potential(displ, potl)
         displ%xyz(ic, iat) = mol%xyz(ic, iat)
         nvat(:, ic, iat) = (potr%vat(:, 1) - potl%vat(:, 1)) / (2*step)
         nvdp(:, :, ic, iat) = (potr%vdp(:, :, 1) - potl%vdp(:, :, 1)) / (2*step)
         nvqp(:, :, ic, iat) = (potr%vqp(:, :, 1) - potl%vqp(:, :, 1)) / (2*step)
      end do
   end do

   if (any(abs(dvat - nvat) > thr) .or. any(abs(dvdp - nvdp) > thr) &
      & .or. any(abs(dvqp - nvqp) > thr)) then
      call test_failed(error, "Potential derivative does not match finite differences")
      print "(a, 3es12.4)", "max. deviation:", maxval(abs(dvat - nvat)), &
         & maxval(abs(dvdp - nvdp)), maxval(abs(dvqp - nvqp))
      print "(a, 3es12.4)", "max. element:  ", maxval(abs(nvat)), &
         & maxval(abs(nvdp)), maxval(abs(nvqp))
   end if

contains

   subroutine get_potential(geo, pot)
      type(structure_type), intent(in) :: geo
      type(potential_type), intent(inout) :: pot

      type(container_cache) :: cache

      call pot%reset
      call calc%coulomb%aes2%update(geo, cache)
      call calc%coulomb%aes2%get_potential(geo, cache, wfn, pot)
   end subroutine get_potential

end subroutine test_multipole_potential

subroutine test_hessian_gen(error, mol, calc, etemp)
   !***********************************
   !* Compare the analytic Hessian against finite differences of the
   !* gradient of fully converged single point calculations.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Electronic temperature
   real(wp), intent(in) :: etemp

   type(context_type) :: ctx
   type(wavefunction_type) :: wfn0
   real(wp) :: energy
   real(wp), allocatable :: hessian(:, :), numhess(:, :)
   real(wp), parameter :: thr_scf = 1.0e-6_wp

   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, 1, etemp)
   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat))
   call xtb_hessian(ctx, mol, calc, wfn0, hessian)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   call get_numerical_hessian(mol, get_gradient, numhess)

   if (any(abs(hessian - numhess) > thr_scf)) then
      call test_failed(error, "Hessian does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(hessian - numhess))
      print "(a, es12.4)", "max. element:  ", maxval(abs(numhess))
   end if

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(wavefunction_type) :: wfn
      real(wp) :: sigma(3, 3)

      wfn = wfn0
      call xtb_singlepoint(ctx, displ, calc, wfn, acc, energy, gradient, sigma, &
         & verbosity=0)
   end subroutine get_gradient

end subroutine test_hessian_gen

subroutine test_hessian_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "01")
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_hessian_gen(error, mol, calc, kt)

end subroutine test_hessian_gfn1

subroutine test_hessian_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "02")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_hessian_gen(error, mol, calc, kt)

end subroutine test_hessian_gfn2

subroutine test_hessian_gfn1_smearing(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "03")
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_hessian_gen(error, mol, calc, 0.02_wp)

end subroutine test_hessian_gfn1_smearing

subroutine test_hessian_gfn1_halogen(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   character(len=*), parameter :: sym(6) = [character(len=4) :: &
      & "Br", "Br", "O", "C", "H", "H"]
   real(wp), parameter :: xyz(3, 6) = reshape([&
      & -1.78533374700000_wp, -3.12608299900000_wp,  0.00000000000000_wp, &
      &  0.00000000000000_wp,  0.81604226400000_wp,  0.00000000000000_wp, &
      &  2.65828699900000_wp,  5.29707580600000_wp,  0.00000000000000_wp, &
      &  4.88597158600000_wp,  4.86116137300000_wp,  0.00000000000000_wp, &
      &  5.61550975300000_wp,  2.90822215900000_wp,  0.00000000000000_wp, &
      &  6.28907612600000_wp,  6.39963643500000_wp,  0.00000000000000_wp], &
      & shape(xyz))

   call new(mol, sym, xyz)
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_hessian_gen(error, mol, calc, kt)

end subroutine test_hessian_gfn1_halogen

subroutine test_dipgrad_gen(error, mol, calc)
   !***********************************
   !* Compare the analytic derivative of the dipole moment against finite
   !* differences of the dipole moment of converged single point calculations.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc

   integer :: iat, ic
   type(context_type) :: ctx
   type(structure_type) :: displ
   type(wavefunction_type) :: wfn0, wfn
   real(wp) :: energy, dipr(3), dipl(3)
   real(wp), allocatable :: hessian(:, :), dipgrad(:, :, :), numdip(:, :, :)
   real(wp), parameter :: step = 1.0e-4_wp, thr_scf = 1.0e-5_wp

   call new_wavefunction(wfn0, mol%nat, calc%bas%nsh, calc%bas%nao, 1, kt)
   call xtb_singlepoint(ctx, mol, calc, wfn0, acc, energy, verbosity=0)

   allocate(hessian(3*mol%nat, 3*mol%nat), dipgrad(3, 3, mol%nat), numdip(3, 3, mol%nat))
   call xtb_hessian(ctx, mol, calc, wfn0, hessian, dipgrad)
   if (ctx%failed()) then
      call ctx%get_error(error)
      return
   end if

   displ = mol
   do iat = 1, mol%nat
      do ic = 1, 3
         wfn = wfn0
         displ%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call xtb_singlepoint(ctx, displ, calc, wfn, acc, energy, verbosity=0)
         dipr(:) = matmul(displ%xyz, wfn%qat(:, 1)) + sum(wfn%dpat(:, :, 1), 2)
         wfn = wfn0
         displ%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call xtb_singlepoint(ctx, displ, calc, wfn, acc, energy, verbosity=0)
         dipl(:) = matmul(displ%xyz, wfn%qat(:, 1)) + sum(wfn%dpat(:, :, 1), 2)
         displ%xyz(ic, iat) = mol%xyz(ic, iat)
         numdip(:, ic, iat) = (dipr - dipl) / (2*step)
      end do
   end do

   if (any(abs(dipgrad - numdip) > thr_scf)) then
      call test_failed(error, "Dipole derivative does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(dipgrad - numdip))
      print "(a, es12.4)", "max. element:  ", maxval(abs(numdip))
   end if

end subroutine test_dipgrad_gen

subroutine test_dipgrad_gfn1(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "01")
   call new_gfn1_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_dipgrad_gen(error, mol, calc)

end subroutine test_dipgrad_gfn1

subroutine test_dipgrad_gfn2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(xtb_calculator) :: calc

   call get_structure(mol, "MB16-43", "02")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call test_dipgrad_gen(error, mol, calc)

end subroutine test_dipgrad_gfn2

subroutine test_unsupported_spin(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(context_type) :: ctx
   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(wavefunction_type) :: wfn
   real(wp), allocatable :: hessian(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_gfn2_calculator(calc, mol, error)
   if (allocated(error)) return
   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 2, kt)

   allocate(hessian(3*mol%nat, 3*mol%nat))
   call xtb_hessian(ctx, mol, calc, wfn, hessian)
   if (ctx%failed()) call ctx%get_error(error)

end subroutine test_unsupported_spin

end module test_hessian
