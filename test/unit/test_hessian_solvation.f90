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

module test_hessian_solvation
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mstore, only : get_structure
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_container, only : container_cache, container_type
   use tblite_hessian_solvation, only : solvation_kernel, new_solvation_kernel, &
      & get_solvation_potential, get_solvation_potential_derivative, get_solvation_hessian
   use tblite_hessian_solvation_born, only : add_psi_hessian, add_born_radii_hessian
   use tblite_hessian_solvation_cm5, only : add_cm5_hessian
   use tblite_hessian_solvation_coupling, only : add_born_coupling_hessian
   use tblite_hessian_solvation_kernel, only : add_born_kernel_hessian
   use tblite_hessian_solvation_shape, only : add_shape_hessian
   use tblite_hessian_solvation_surface, only : add_surface_hessian
   use tblite_solvation_alpb, only : alpb_solvation, alpb_input, new_alpb, born_kernel, &
      & get_adet_deriv
   use tblite_solvation_born, only : born_integrator, new_born_integrator, compute_psi
   use tblite_solvation_cm5, only : get_cm5_charges
   use tblite_scf_potential, only : potential_type, new_potential
   use tblite_solvation, only : solvation_input, cds_input, shift_input, solvent_data, &
      & get_solvent_data, new_solvation, new_solvation_cds, new_solvation_shift, &
      & solvation_type
   use tblite_solvation_data, only : get_vdw_rad_d3
   use tblite_solvation_surface, only : surface_integrator, new_surface_integrator
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   implicit none
   private

   public :: collect_hessian_solvation

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: step = 1.0e-4_wp
   !> Dielectric constant of water
   real(wp), parameter :: eps_water = 80.2_wp

   abstract interface
      !> Evaluate the gradient of a weighted sum of atomic quantities
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
subroutine collect_hessian_solvation(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("descreening", test_psi), &
      new_unittest("born-radii", test_born_radii), &
      new_unittest("born-p16", test_born_p16), &
      new_unittest("born-still", test_born_still), &
      new_unittest("alpb-charged", test_alpb_charged), &
      new_unittest("shape", test_shape), &
      new_unittest("surface", test_surface), &
      new_unittest("cm5", test_cm5), &
      new_unittest("kernel-gfn2-alpb", test_kernel_gfn2_alpb), &
      new_unittest("kernel-gfn1-gbsa", test_kernel_gfn1_gbsa), &
      new_unittest("potential-gfn2-alpb", test_potential_gfn2_alpb), &
      new_unittest("potential-gfn1-alpb", test_potential_gfn1_alpb), &
      new_unittest("potential-gfn2-gbsa", test_potential_gfn2_gbsa), &
      new_unittest("hessian-gfn2-alpb", test_hessian_gfn2_alpb), &
      new_unittest("hessian-gfn2-gbsa", test_hessian_gfn2_gbsa), &
      new_unittest("hessian-gfn1-alpb", test_hessian_gfn1_alpb) &
      ]

end subroutine collect_hessian_solvation


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

subroutine compare_hessian(error, hessian, numhess, tol)
   !***********************************
   !* Check an analytic Hessian for symmetry and against a numerical one.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Analytic Hessian
   real(wp), intent(in) :: hessian(:, :)
   !> Numerical Hessian
   real(wp), intent(in) :: numhess(:, :)
   !> Tolerance for the deviation from the numerical Hessian
   real(wp), intent(in) :: tol

   if (any(abs(hessian - transpose(hessian)) > thr)) then
      call test_failed(error, "Hessian is not symmetric")
      print "(a, es12.4)", "max. asymmetry:", maxval(abs(hessian - transpose(hessian)))
      return
   end if
   if (any(abs(hessian - numhess) > tol)) then
      call test_failed(error, "Hessian does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(hessian - numhess))
      print "(a, es12.4)", "max. element:  ", maxval(abs(numhess))
   end if
end subroutine compare_hessian

pure function get_weights(nat) result(weight)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Weights of the atomic quantities, without any symmetry
   real(wp) :: weight(nat)

   integer :: iat

   do iat = 1, nat
      weight(iat) = cos(1.3_wp*iat) + 0.1_wp*iat
   end do
end function get_weights

pure function get_charges(nat) result(qat)
   !> Number of atoms
   integer, intent(in) :: nat
   !> Atomic charges of a neutral molecule
   real(wp) :: qat(nat)

   integer :: iat

   do iat = 1, nat
      qat(iat) = 0.4_wp*sin(0.7_wp*iat)
   end do
   qat(:) = qat - sum(qat)/nat
end function get_charges


subroutine test_psi(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: weight(:), hessian(:, :), numhess(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_born_integrator(gbobc, mol, get_vdw_rad_d3(mol%num))
   weight = get_weights(mol%nat)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call new_adjacency_list(list, mol, trans, gbobc%lrcut)
   call add_psi_hessian(mol, list, gbobc%vdwr, gbobc%rho, weight, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess, thr)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(adjacency_list) :: dlist
      real(wp), allocatable :: psi(:), dpsidr(:, :, :)
      integer :: iat

      allocate(psi(displ%nat), dpsidr(3, displ%nat, displ%nat))
      call new_adjacency_list(dlist, displ, trans, gbobc%lrcut)
      call compute_psi(displ%nat, displ%xyz, dlist, gbobc%vdwr, gbobc%rho, psi, dpsidr)
      gradient(:, :) = 0.0_wp
      do iat = 1, displ%nat
         gradient(:, :) = gradient + weight(iat)*dpsidr(:, :, iat)
      end do
   end subroutine get_gradient

end subroutine test_psi


subroutine test_born_radii(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(born_integrator) :: gbobc
   real(wp), allocatable :: weight(:), hessian(:, :), numhess(:, :)

   call get_structure(mol, "MB16-43", "01")
   call new_born_integrator(gbobc, mol, get_vdw_rad_d3(mol%num))
   weight = get_weights(mol%nat)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call add_born_radii_hessian(gbobc, mol, weight, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess, thr)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      real(wp), allocatable :: rad(:), draddr(:, :, :)
      integer :: iat

      allocate(rad(displ%nat), draddr(3, displ%nat, displ%nat))
      call gbobc%get_rad(displ, rad, draddr)
      gradient(:, :) = 0.0_wp
      do iat = 1, displ%nat
         gradient(:, :) = gradient + weight(iat)*draddr(:, :, iat)
      end do
   end subroutine get_gradient

end subroutine test_born_radii


subroutine test_born_gen(error, mol, input, qat)
   !***********************************
   !* Compare the Hessian of the generalized Born energy at fixed charges
   !* against finite differences of the gradient of the container.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Input for the generalized Born model
   type(alpb_input), intent(in) :: input
   !> Atomic charges
   real(wp), intent(in) :: qat(:)

   type(alpb_solvation) :: solv
   real(wp), allocatable :: rad(:), draddr(:, :, :), dEda(:), hessian(:, :), numhess(:, :)

   call new_alpb(solv, mol, input)

   allocate(rad(mol%nat), draddr(3, mol%nat, mol%nat), dEda(mol%nat))
   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call solv%gbobc%get_rad(mol, rad, draddr)
   call add_born_kernel_hessian(mol, solv%kernel, solv%keps, qat, rad, hessian)
   call add_born_coupling_hessian(mol, solv%kernel, solv%keps, qat, rad, draddr, &
      & hessian, dEda)
   call add_born_radii_hessian(solv%gbobc, mol, dEda, hessian)
   call add_shape_hessian(mol, solv%gbobc%vdwr, solv%keps*solv%alpbet, sum(qat), hessian)

   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess, thr)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      type(container_cache) :: cache
      type(wavefunction_type) :: wfn
      real(wp) :: sigma(3, 3)

      call new_wavefunction(wfn, displ%nat, 0, 0, 1, 0.0_wp)
      wfn%qat(:, 1) = qat
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      call solv%update(displ, cache)
      call solv%get_gradient(displ, cache, wfn, gradient, sigma)
   end subroutine get_gradient

end subroutine test_born_gen


subroutine test_born_p16(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_born_gen(error, mol, &
      & alpb_input(eps_water, kernel=born_kernel%p16, alpb=.false.), get_charges(mol%nat))

end subroutine test_born_p16


subroutine test_born_still(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call test_born_gen(error, mol, &
      & alpb_input(eps_water, kernel=born_kernel%still, alpb=.false.), get_charges(mol%nat))

end subroutine test_born_still


subroutine test_alpb_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: qat(:)

   call get_structure(mol, "MB16-43", "01")
   qat = get_charges(mol%nat)
   qat(1) = qat(1) + 1.0_wp
   call test_born_gen(error, mol, &
      & alpb_input(eps_water, kernel=born_kernel%p16, alpb=.true.), qat)

end subroutine test_alpb_charged


subroutine test_shape(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: rad(:), qat(:), hessian(:, :), numhess(:, :)
   real(wp), parameter :: kepsalpbet = -0.3_wp

   call get_structure(mol, "MB16-43", "02")
   rad = get_vdw_rad_d3(mol%num)
   qat = get_charges(mol%nat)
   qat(2) = qat(2) - 1.0_wp

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call add_shape_hessian(mol, rad, kepsalpbet, sum(qat), hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   !>-- the shape term is small, an absolute threshold would not resolve it
   call compare_hessian(error, hessian, numhess, thr*maxval(abs(numhess)))

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      gradient(:, :) = 0.0_wp
      call get_adet_deriv(displ%nat, displ%xyz, rad, kepsalpbet, qat, gradient)
   end subroutine get_gradient

end subroutine test_shape


subroutine test_surface(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(surface_integrator) :: sasa
   real(wp), allocatable :: weight(:), hessian(:, :), numhess(:, :)
   real(wp), parameter :: probe = 1.4_wp * aatoau
   integer, parameter :: nang = 110

   call get_structure(mol, "MB16-43", "01")
   call new_surface_integrator(sasa, mol%id, get_vdw_rad_d3(mol%num), probe, nang)
   weight = get_weights(mol%nat)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call add_surface_hessian(sasa, mol, weight, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   !>-- the large third derivatives of the switching function dominate the
   !>   finite difference error, which decreases with the square of the step
   call compare_hessian(error, hessian, numhess, 1.0e-7_wp*maxval(abs(numhess)))

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      real(wp), allocatable :: surface(:), dsdr(:, :, :)
      integer :: iat

      allocate(surface(displ%nat), dsdr(3, displ%nat, displ%nat))
      call sasa%get_surface(displ, surface, dsdr)
      gradient(:, :) = 0.0_wp
      do iat = 1, displ%nat
         gradient(:, :) = gradient + weight(iat)*dsdr(:, :, iat)
      end do
   end subroutine get_gradient

end subroutine test_surface


subroutine test_cm5(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: weight(:), hessian(:, :), numhess(:, :)

   call get_structure(mol, "MB16-43", "01")
   weight = get_weights(mol%nat)

   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call add_cm5_hessian(mol, weight, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   call compare_hessian(error, hessian, numhess, thr)

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      real(wp), allocatable :: cm5(:), dcm5dr(:, :, :)
      integer :: iat

      allocate(cm5(displ%nat), dcm5dr(3, displ%nat, displ%nat))
      call get_cm5_charges(displ, cm5, dcm5dr)
      gradient(:, :) = 0.0_wp
      do iat = 1, displ%nat
         gradient(:, :) = gradient + weight(iat)*dcm5dr(:, :, iat)
      end do
   end subroutine get_gradient

end subroutine test_cm5


subroutine get_solvated_calculator(calc, mol, method, alpb, kernel, error)
   !***********************************
   !* Calculator with the parametrized solvation models for water attached
   !* in the same way as by the command line driver.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(out) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Method for parameter selection
   character(len=*), intent(in) :: method
   !> Use the analytical linearized Poisson-Boltzmann model
   logical, intent(in) :: alpb
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(solvation_input) :: input
   type(solvent_data) :: solvent
   class(container_type), allocatable :: cont
   class(solvation_type), allocatable :: solv

   select case(method)
   case("gfn1")
      call new_gfn1_calculator(calc, mol, error)
   case default
      call new_gfn2_calculator(calc, mol, error)
   end select
   if (allocated(error)) return

   solvent = get_solvent_data("water")
   input%alpb = alpb_input(solvent%eps, solvent=solvent%solvent, kernel=kernel, alpb=alpb)
   input%cds = cds_input(alpb=alpb, solvent=solvent%solvent)
   input%shift = shift_input(alpb=alpb, solvent=solvent%solvent)

   call new_solvation(solv, mol, input, error, method)
   if (allocated(error)) return
   call move_alloc(solv, cont)
   call calc%push_back(cont)
   call new_solvation_cds(solv, mol, input, error, method)
   if (allocated(error)) return
   call move_alloc(solv, cont)
   call calc%push_back(cont)
   call new_solvation_shift(solv, input, error, method)
   if (allocated(error)) return
   call move_alloc(solv, cont)
   call calc%push_back(cont)
end subroutine get_solvated_calculator

subroutine get_container_potential(calc, mol, qat, vat)
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Atom-resolved potential of the interaction containers
   real(wp), intent(out) :: vat(:)

   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   type(potential_type) :: pot

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 0.0_wp)
   call new_potential(pot, mol, calc%bas, 1)
   wfn%qat(:, 1) = qat
   call pot%reset
   call calc%interactions%update(mol, cache)
   call calc%interactions%get_potential(mol, cache, wfn, pot)
   vat(:) = pot%vat(:, 1)
end subroutine get_container_potential

subroutine get_container_gradient(calc, mol, qat, gradient)
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Gradient of the interaction containers at fixed charges
   real(wp), intent(out) :: gradient(:, :)

   type(container_cache) :: cache
   type(wavefunction_type) :: wfn
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: energies(:)

   call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, 1, 0.0_wp)
   wfn%qat(:, 1) = qat
   allocate(energies(mol%nat), source=0.0_wp)
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call calc%interactions%update(mol, cache)
   call calc%interactions%get_engrad(mol, cache, energies, gradient, sigma)
   call calc%interactions%get_gradient(mol, cache, wfn, gradient, sigma)
end subroutine get_container_gradient


subroutine test_kernel_gen(error, method, alpb, kernel)
   !***********************************
   !* The potential of the containers is affine in the charges, which makes
   !* v(e_A) - v(0) the column A of the response kernel up to rounding.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), intent(in) :: method
   !> Use the analytical linearized Poisson-Boltzmann model
   logical, intent(in) :: alpb
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(solvation_kernel) :: solv
   integer :: iat
   real(wp), allocatable :: qat(:), vzero(:), vat(:), kmat(:, :)

   call get_structure(mol, "MB16-43", "01")
   call get_solvated_calculator(calc, mol, method, alpb, kernel, error)
   if (allocated(error)) return
   call new_solvation_kernel(solv, mol, calc, error)
   if (allocated(error)) return

   allocate(qat(mol%nat), vzero(mol%nat), vat(mol%nat), kmat(mol%nat, mol%nat))
   qat(:) = 0.0_wp
   call get_container_potential(calc, mol, qat, vzero)
   do iat = 1, mol%nat
      qat(:) = 0.0_wp
      qat(iat) = 1.0_wp
      call get_container_potential(calc, mol, qat, vat)
      kmat(:, iat) = vat - vzero
   end do

   if (any(abs(solv%kmat - kmat) > thr)) then
      call test_failed(error, "Response kernel does not match the containers")
      print "(a, es12.4)", "max. deviation:", maxval(abs(solv%kmat - kmat))
   end if
end subroutine test_kernel_gen

subroutine test_kernel_gfn2_alpb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_kernel_gen(error, "gfn2", .true., born_kernel%p16)
end subroutine test_kernel_gfn2_alpb

subroutine test_kernel_gfn1_gbsa(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_kernel_gen(error, "gfn1", .false., born_kernel%still)
end subroutine test_kernel_gfn1_gbsa


subroutine test_potential_gen(error, method, alpb, kernel)
   !***********************************
   !* Compare the potential at fixed charges against the containers and its
   !* derivative w.r.t. the positions against finite differences.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), intent(in) :: method
   !> Use the analytical linearized Poisson-Boltzmann model
   logical, intent(in) :: alpb
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel

   type(structure_type) :: mol, displ
   type(xtb_calculator) :: calc
   type(solvation_kernel) :: solv
   integer :: iat, ic
   real(wp), allocatable :: qat(:), vat(:), vref(:), vr(:), vl(:)
   real(wp), allocatable :: dvat(:, :, :), numdv(:, :, :)

   call get_structure(mol, "MB16-43", "01")
   call get_solvated_calculator(calc, mol, method, alpb, kernel, error)
   if (allocated(error)) return
   call new_solvation_kernel(solv, mol, calc, error)
   if (allocated(error)) return

   qat = get_charges(mol%nat)
   qat(3) = qat(3) - 1.0_wp
   allocate(vat(mol%nat), vref(mol%nat), vr(mol%nat), vl(mol%nat))
   allocate(dvat(mol%nat, 3, mol%nat), numdv(mol%nat, 3, mol%nat))

   call get_solvation_potential(solv, qat, vat)
   call get_container_potential(calc, mol, qat, vref)
   if (any(abs(vat - vref) > thr)) then
      call test_failed(error, "Potential does not match the containers")
      print "(a, es12.4)", "max. deviation:", maxval(abs(vat - vref))
      return
   end if

   call get_solvation_potential_derivative(solv, mol, qat, dvat)
   displ = mol
   do iat = 1, mol%nat
      do ic = 1, 3
         displ%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_container_potential(calc, displ, qat, vr)
         displ%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call get_container_potential(calc, displ, qat, vl)
         displ%xyz(ic, iat) = mol%xyz(ic, iat)
         numdv(:, ic, iat) = (vr - vl) / (2*step)
      end do
   end do

   !>-- the surface area limits the accuracy of the finite differences
   if (any(abs(dvat - numdv) > 1.0e-7_wp*maxval(abs(numdv)))) then
      call test_failed(error, "Potential derivative does not match finite differences")
      print "(a, es12.4)", "max. deviation:", maxval(abs(dvat - numdv))
      print "(a, es12.4)", "max. element:  ", maxval(abs(numdv))
   end if
end subroutine test_potential_gen

subroutine test_potential_gfn2_alpb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_potential_gen(error, "gfn2", .true., born_kernel%p16)
end subroutine test_potential_gfn2_alpb

subroutine test_potential_gfn1_alpb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_potential_gen(error, "gfn1", .true., born_kernel%p16)
end subroutine test_potential_gfn1_alpb

subroutine test_potential_gfn2_gbsa(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_potential_gen(error, "gfn2", .false., born_kernel%still)
end subroutine test_potential_gfn2_gbsa


subroutine test_hessian_gen(error, method, alpb, kernel, charge)
   !***********************************
   !* Compare the Hessian of all solvation models at fixed charges against
   !* finite differences of the gradient of the containers.
   !***********************************
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Method for parameter selection
   character(len=*), intent(in) :: method
   !> Use the analytical linearized Poisson-Boltzmann model
   logical, intent(in) :: alpb
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel
   !> Total charge
   real(wp), intent(in) :: charge

   type(structure_type) :: mol
   type(xtb_calculator) :: calc
   type(solvation_kernel) :: solv
   real(wp), allocatable :: qat(:), hessian(:, :), numhess(:, :)

   call get_structure(mol, "MB16-43", "01")
   call get_solvated_calculator(calc, mol, method, alpb, kernel, error)
   if (allocated(error)) return
   call new_solvation_kernel(solv, mol, calc, error)
   if (allocated(error)) return

   qat = get_charges(mol%nat)
   qat(5) = qat(5) + charge
   allocate(hessian(3*mol%nat, 3*mol%nat), numhess(3*mol%nat, 3*mol%nat), source=0.0_wp)
   call get_solvation_hessian(solv, mol, qat, hessian)
   call get_numerical_hessian(mol, get_gradient, numhess)
   !>-- the surface area limits the accuracy of the finite differences
   call compare_hessian(error, hessian, numhess, 1.0e-7_wp*maxval(abs(numhess)))

contains

   subroutine get_gradient(displ, gradient)
      type(structure_type), intent(in) :: displ
      real(wp), intent(out) :: gradient(:, :)

      call get_container_gradient(calc, displ, qat, gradient)
   end subroutine get_gradient

end subroutine test_hessian_gen

subroutine test_hessian_gfn2_alpb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_hessian_gen(error, "gfn2", .true., born_kernel%p16, 1.0_wp)
end subroutine test_hessian_gfn2_alpb

subroutine test_hessian_gfn2_gbsa(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_hessian_gen(error, "gfn2", .false., born_kernel%still, 0.0_wp)
end subroutine test_hessian_gfn2_gbsa

subroutine test_hessian_gfn1_alpb(error)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call test_hessian_gen(error, "gfn1", .true., born_kernel%p16, -1.0_wp)
end subroutine test_hessian_gfn1_alpb

end module test_hessian_solvation
