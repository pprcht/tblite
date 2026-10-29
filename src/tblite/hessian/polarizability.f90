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

!> @file tblite/hessian/polarizability.f90
!> Provides the analytic derivative of the polarizability

!> Derivative of the static dipole polarizability of a self-consistent xTB
!> calculation w.r.t. the atomic positions
module tblite_hessian_polarizability
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_hessian_coulomb, only : get_potential_derivative
   use tblite_hessian_dispersion, only : get_dispersion_potential_derivative
   use tblite_hessian_solvation, only : get_solvation_kernel_derivative
   use tblite_hessian_multipole, only : get_multipole_potential_derivative
   use tblite_hessian_reference, only : reference_data, new_reference_data, &
      & get_nuclear_perturbation, add_atomic_potential_derivative, add_multipole_fock, &
      & check_hessian_support
   use tblite_response_cpscc, only : solve_cpscc, get_fock_response
   use tblite_response_density, only : density_response, density_response_derivative, &
      & has_integer_occupation
   use tblite_response_field, only : xtb_field_response
   use tblite_response_integral, only : get_multipole_derivative
   use tblite_response_kernel, only : get_kernel_potential, &
      & get_kernel_potential_derivative, has_dispersion_kernel
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_polarizability, only : check_response_support
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_polarizability_gradient

   !> First-order response to a homogeneous electric field in one direction
   type :: field_data
      !> First-order density matrix and Hamiltonian
      real(wp), allocatable :: dp(:, :, :), dh(:, :, :)
      !> First-order potential including the field
      type(potential_type) :: pot
      !> First-order charges and multipole moments
      type(wavefunction_type) :: wfn
      !> Derivative of the first-order shell-resolved potential w.r.t. the
      !> atomic positions at fixed first-order moments
      real(wp), allocatable :: dvshdr(:, :)
      !> Derivative of the first-order multipole potentials
      real(wp), allocatable :: dvdp(:, :, :, :), dvqp(:, :, :, :)
   end type field_data

contains

subroutine xtb_polarizability_gradient(ctx, mol, calc, wfn, dalpha)
   !***********************************
   !* Derivative of the polarizability w.r.t. the atomic positions,
   !* dalpha(a, b, ic, iat) for the displacement of atom iat in direction ic.
   !* The coupled-perturbed equations of a nuclear displacement are
   !* differentiated w.r.t. the field. The resulting equations for the mixed
   !* second-order density matrix have the same response operator, their
   !* right-hand side contains first-order quantities only.
   !* Requires integer occupation numbers and the converged wavefunction of
   !* a preceding single point calculation with the same calculator.
   !* Failures are reported in ctx.
   !***********************************
   !> Calculation context
   type(context_type), intent(inout) :: ctx
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Derivative of the polarizability tensor
   real(wp), intent(out) :: dalpha(:, :, :, :)

   integer :: nat, nao, nsh, iat, idir, ish
   real(wp), allocatable :: dvat(:, :, :)
   type(reference_data) :: ref
   type(field_data) :: field(3)
   type(container_cache) :: ccache, dcache
   type(error_type), allocatable :: error, local_error

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   dalpha(:, :, :, :) = 0.0_wp

   call check_response_support(mol, calc, wfn, error)
   if (.not.allocated(error)) call check_hessian_support(calc, error)
   if (.not.allocated(error) .and. .not.has_integer_occupation(wfn)) then
      call fatal_error(error, "Analytic polarizability derivative requires integer "// &
         & "occupation numbers")
   end if
   if (.not.allocated(error)) call new_reference_data(ref, mol, calc, wfn, error)
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if

   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)

   !>-- field response and the derivative of its potential at fixed moments
   allocate(dvat(nat, 3, nat))
   do idir = 1, 3
      associate(fld => field(idir))
         allocate(fld%dp(nao, nao, 1), fld%dh(nao, nao, 1), fld%dvshdr(nsh, 3*nat))
         allocate(fld%dvdp(3, nat, 3, nat), fld%dvqp(6, nat, 3, nat), source=0.0_wp)
         call new_potential(fld%pot, mol, calc%bas, 1)
         call new_wavefunction(fld%wfn, nat, nsh, nao, 1, wfn%kt)
         call xtb_field_response(mol, calc, wfn, ref%ints, ccache, dcache, idir, fld%dp, &
            & fld%dh, fld%pot, fld%wfn, error, kat=ref%solv%kmat)
         if (allocated(error)) exit

         call get_potential_derivative(mol, calc, ccache, fld%wfn, fld%dvshdr)
         if (allocated(calc%coulomb%aes2)) then
            call get_multipole_potential_derivative(calc%coulomb%aes2, mol, &
               & fld%wfn%qat(:, 1), fld%wfn%dpat(:, :, 1), fld%wfn%qpat(:, :, 1), dvat, &
               & fld%dvdp, fld%dvqp)
            call add_atomic_potential_derivative(calc, dvat, fld%dvshdr)
         end if
         if (has_dispersion_kernel(calc)) then
            call get_dispersion_potential_derivative(calc%dispersion, mol, wfn%qat(:, 1), &
               & dvat, error, dqat=fld%wfn%qat(:, 1))
            if (allocated(error)) exit
            call add_atomic_potential_derivative(calc, dvat, fld%dvshdr)
         end if
         if (allocated(calc%interactions)) then
            call get_solvation_kernel_derivative(ref%solv, mol, fld%wfn%qat(:, 1), dvat)
            call add_atomic_potential_derivative(calc, dvat, fld%dvshdr)
         end if
         !>-- the potential of a unit field on an atom is minus its position
         do iat = 1, nat
            do ish = 1, calc%bas%nsh_at(iat)
               fld%dvshdr(calc%bas%ish_at(iat)+ish, 3*(iat-1)+idir) = &
                  & fld%dvshdr(calc%bas%ish_at(iat)+ish, 3*(iat-1)+idir) - 1.0_wp
            end do
         end do
      end associate
   end do
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if

   !$omp parallel do schedule(dynamic) default(none) &
   !$omp shared(mol, calc, wfn, ref, field, nat, dalpha, error) &
   !$omp private(iat, local_error)
   do iat = 1, nat
      call get_mixed_response(mol, calc, wfn, ref, field, iat, dalpha(:, :, :, iat), &
         & local_error)
      if (allocated(local_error)) then
         !$omp critical (tblite_hessian_polarizability_error)
         if (.not.allocated(error)) call move_alloc(local_error, error)
         !$omp end critical (tblite_hessian_polarizability_error)
      end if
   end do
   if (allocated(error)) call ctx%set_error(error)

end subroutine xtb_polarizability_gradient

subroutine get_mixed_response(mol, calc, wfn, ref, field, iat, dalpha, error)
   !***********************************
   !* Derivative of the polarizability w.r.t. the three displacements of one
   !* atom. For the displacement c and the field b the right-hand side of the
   !* coupled-perturbed equations consists of the derivative of the nuclear
   !* density response w.r.t. the state changed by the field, and of the
   !* uncoupled response to the mixed derivative of the Hamiltonian at fixed
   !* second-order moments. The latter contains the derivative of the
   !* first-order field Hamiltonian w.r.t. the position, the potential of
   !* the change of the field moments due to the derivative of the integrals
   !* and the derivative of the kernel acting on the nuclear moments.
   !* All work arrays are local to allow evaluating atoms in parallel.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Quantities of the unperturbed system
   type(reference_data), intent(in) :: ref
   !> Field response for the three directions
   type(field_data), intent(in) :: field(3)
   !> Displaced atom
   integer, intent(in) :: iat
   !> Derivative of the polarizability, the last dimension is the direction
   real(wp), intent(out) :: dalpha(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh, ic, idir, iao, jao
   logical :: multipoles
   real(wp), allocatable :: ds(:, :, :), fock0(:, :, :), dqsh(:, :), ddpat(:, :, :)
   real(wp), allocatable :: dqpat(:, :, :), dp(:, :, :), dw(:, :), fock_r(:, :)
   real(wp), allocatable :: dfock(:, :, :, :), ddpat_x(:, :, :, :), dqpat_x(:, :, :, :)
   real(wp), allocatable :: fock(:, :, :), rhs(:, :, :), tmp(:, :), dpx(:, :, :)
   real(wp), allocatable :: dvao(:), n0sh(:), qsh(:, :), qat(:, :), dpat(:, :, :)
   type(container_cache) :: ccache, dcache
   type(potential_type) :: pot_r, pot_x
   type(wavefunction_type) :: wfn_r, wfn_m, wfn_x

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)
   dalpha(:, :, :) = 0.0_wp

   !>-- the caches serve as scratch space of the containers
   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)
   call new_potential(pot_r, mol, calc%bas, 1)
   call new_potential(pot_x, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_x, nat, nsh, nao, 1, wfn%kt)
   allocate(ds(nao, nao, 3), fock0(nao, nao, 3), dqsh(nsh, 3), ddpat(3, nat, 3), &
      & dqpat(6, nat, 3), dp(nao, nao, 1), dw(nao, nao), fock_r(nao, nao), &
      & ddpat_x(3, nat, 3, 3), dqpat_x(6, nat, 3, 3), fock(nao, nao, 1), rhs(nao, nao, 1), &
      & tmp(nao, nao), dpx(nao, nao, 1), dvao(nao), qsh(nsh, 1), qat(nat, 1), &
      & dpat(3, nat, 1))
   allocate(dfock(nao, nao, 3, 3), n0sh(nsh), source=0.0_wp)

   call get_nuclear_perturbation(mol, calc, wfn, ref, iat, .true., ds, fock0, dqsh, &
      & ddpat, dqpat)

   !>-- derivative of the multipole integrals with the field response
   do idir = 1, 3
      call get_multipole_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, &
         & field(idir)%dp(:, :, 1), field(idir)%pot%vdp(:, :, 1), &
         & field(idir)%pot%vqp(:, :, 1), dfock(:, :, :, idir), ddpat_x(:, :, :, idir), &
         & dqpat_x(:, :, :, idir))
   end do

   do ic = 1, 3
      wfn_m%qsh(:, 1) = dqsh(:, ic)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if
      call xtb_nuclear_response(mol, calc, wfn, ref%ints, ccache, dcache, ds(:, :, ic), &
         & fock0(:, :, ic), wfn_m, dp, dw, pot_r, wfn_r, error, fock_r, kat=ref%solv%kmat)
      if (allocated(error)) return

      do idir = 1, 3
         associate(fld => field(idir))
            !>-- change of the field moments due to the derivative of the integrals
            call get_mulliken_shell_charges(calc%bas, ds(:, :, ic), fld%dp, n0sh, wfn_x%qsh)
            call get_qat_from_qsh(calc%bas, wfn_x%qsh, wfn_x%qat)
            if (multipoles) then
               wfn_x%dpat(:, :, 1) = ddpat_x(:, :, ic, idir)
               wfn_x%qpat(:, :, 1) = dqpat_x(:, :, ic, idir)
            end if

            call pot_x%reset
            call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_x, pot_x, error, &
               & kat=ref%solv%kmat)
            if (allocated(error)) return
            call get_kernel_potential_derivative(mol, calc, dcache, wfn, fld%wfn, wfn_r, &
               & pot_x, error)
            if (allocated(error)) return
            call get_fock_response(calc, ref%ints, pot_x, fock)

            !>-- derivative of the field Hamiltonian at fixed first-order moments
            do iao = 1, nao
               dvao(iao) = fld%dvshdr(calc%bas%ao2sh(iao), 3*(iat-1)+ic)
            end do
            do iao = 1, nao
               do jao = 1, nao
                  fock(jao, iao, 1) = fock(jao, iao, 1) + dfock(jao, iao, ic, idir) &
                     & - 0.5_wp*(fld%pot%vao(jao, 1) + fld%pot%vao(iao, 1))*ds(jao, iao, ic) &
                     & - 0.5_wp*(dvao(jao) + dvao(iao))*ref%ints%overlap(jao, iao)
               end do
            end do
            if (multipoles) then
               call add_multipole_fock(calc, ref%ints, fld%dvdp(:, :, ic, iat), &
                  & fld%dvqp(:, :, ic, iat), fock(:, :, 1))
            end if

            call density_response(wfn, fock(:, :, 1), rhs(:, :, 1))
            call density_response_derivative(wfn, fld%dh(:, :, 1), fock_r, ds(:, :, ic), tmp)
            rhs(:, :, 1) = rhs(:, :, 1) + tmp

            call solve_cpscc(mol, calc, wfn, ref%ints, ccache, dcache, rhs, dpx, error, &
               & kat=ref%solv%kmat)
            if (allocated(error)) return

            !>-- mu = sum_A (R_A q_A + d_A), the position enters explicitly
            call get_mulliken_shell_charges(calc%bas, ref%ints%overlap, dpx, n0sh, qsh)
            qsh(:, :) = qsh + wfn_x%qsh
            call get_qat_from_qsh(calc%bas, qsh, qat)
            call get_mulliken_atomic_multipoles(calc%bas, ref%ints%dipole, dpx, dpat)
            dalpha(:, idir, ic) = matmul(mol%xyz, qat(:, 1)) &
               & + sum(dpat(:, :, 1) + ddpat_x(:, :, ic, idir), dim=2)
            dalpha(ic, idir, ic) = dalpha(ic, idir, ic) + fld%wfn%qat(iat, 1)
         end associate
      end do

      dalpha(:, :, ic) = 0.5_wp * (dalpha(:, :, ic) + transpose(dalpha(:, :, ic)))
   end do
end subroutine get_mixed_response

end module tblite_hessian_polarizability
