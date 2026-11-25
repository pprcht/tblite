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
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_hessian_reference, only : reference_data, get_nuclear_perturbation, &
      & add_potential_fock
   use tblite_response_cpscc, only : solve_cpscc, get_fock_response
   use tblite_response_density, only : density_response, density_response_derivative
   use tblite_response_integral, only : get_multipole_derivative
   use tblite_response_kernel, only : interaction_cache, get_response_kernel, &
      & get_response_derivative
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_polarizability, only : field_data
   use tblite_response_type, only : kernel_type, response_type, new_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_polarizability_gradient

contains

subroutine xtb_polarizability_gradient(mol, calc, wfn, ref, cache, kernel, field, dalpha, &
      & error)
   !***********************************
   !* Derivative of the polarizability w.r.t. the atomic positions,
   !* dalpha(a, b, ic, iat) for the displacement of atom iat in direction ic.
   !* The coupled-perturbed equations of a nuclear displacement are
   !* differentiated w.r.t. the field. The resulting equations for the mixed
   !* second-order density matrix have the same response operator, their
   !* right-hand side contains first-order quantities only.
   !* Requires integer occupation numbers and the field responses of the
   !* polarizability, whose kernel derivatives are set up here.
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
   !> Response kernel at the converged state
   type(kernel_type), intent(in) :: kernel
   !> Field response for the three directions
   type(field_data), intent(inout) :: field(3)
   !> Derivative of the polarizability tensor
   real(wp), intent(out) :: dalpha(:, :, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, iat, idir
   type(error_type), allocatable :: local_error

   nat = mol%nat
   dalpha(:, :, :, :) = 0.0_wp

   !>-- derivative of the kernel along the field response
   do idir = 1, 3
      call get_response_kernel(mol, calc, cache, wfn, field(idir)%kernel, error, &
         & dwfn=field(idir)%wfn)
      if (allocated(error)) return
   end do

   !$omp parallel do schedule(dynamic) default(none) &
   !$omp shared(mol, calc, wfn, ref, cache, kernel, field, nat, dalpha, error) &
   !$omp private(iat, local_error)
   do iat = 1, nat
      call get_mixed_response(mol, calc, wfn, ref, cache, kernel, field, iat, &
         & dalpha(:, :, :, iat), local_error)
      if (allocated(local_error)) then
         !$omp critical (tblite_hessian_polarizability_error)
         if (.not.allocated(error)) call move_alloc(local_error, error)
         !$omp end critical (tblite_hessian_polarizability_error)
      end if
   end do

end subroutine xtb_polarizability_gradient

subroutine get_mixed_response(mol, calc, wfn, ref, cache, kernel, field, iat, dalpha, error)
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
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Response kernel at the converged state
   type(kernel_type), intent(in) :: kernel
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
   real(wp), allocatable :: n0sh(:), qsh(:, :), qat(:, :), dpat(:, :, :)
   type(potential_type) :: pot_r, pot_x
   type(wavefunction_type) :: wfn_r, wfn_m, wfn_x
   type(response_type) :: dfield(3)

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)
   dalpha(:, :, :) = 0.0_wp

   call new_potential(pot_r, mol, calc%bas, 1)
   call new_potential(pot_x, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_x, nat, nsh, nao, 1, wfn%kt)
   allocate(ds(nao, nao, 3), fock0(nao, nao, 3), dqsh(nsh, 3), ddpat(3, nat, 3), &
      & dqpat(6, nat, 3), dp(nao, nao, 1), dw(nao, nao), fock_r(nao, nao), &
      & ddpat_x(3, nat, 3, 3), dqpat_x(6, nat, 3, 3), fock(nao, nao, 1), rhs(nao, nao, 1), &
      & tmp(nao, nao), dpx(nao, nao, 1), qsh(nsh, 1), qat(nat, 1), dpat(3, nat, 1))
   allocate(dfock(nao, nao, 3, 3), n0sh(nsh), source=0.0_wp)

   call get_nuclear_perturbation(mol, calc, wfn, ref, cache, iat, .true., ds, fock0, dqsh, &
      & ddpat, dqpat, error)
   if (allocated(error)) return

   do idir = 1, 3
      !>-- derivative of the multipole integrals with the field response
      call get_multipole_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, &
         & field(idir)%dp(:, :, 1), field(idir)%pot%vdp(:, :, 1), &
         & field(idir)%pot%vqp(:, :, 1), dfock(:, :, :, idir), ddpat_x(:, :, :, idir), &
         & dqpat_x(:, :, :, idir))

      !>-- derivative of the first-order field potential at fixed first-order
      !>   moments, the potential of a unit field on an atom is minus its position
      call new_response(dfield(idir), mol, calc%bas, [3*iat-2, 3*iat-1, 3*iat])
      call get_response_derivative(mol, calc, cache, wfn, dfield(idir), error, &
         & dwfn=field(idir)%wfn)
      if (allocated(error)) return
      dfield(idir)%vat(iat, idir) = dfield(idir)%vat(iat, idir) - 1.0_wp
   end do

   do ic = 1, 3
      wfn_m%qsh(:, 1) = dqsh(:, ic)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if
      call xtb_nuclear_response(mol, calc, wfn, ref%ints, kernel, ds(:, :, ic), &
         & fock0(:, :, ic), wfn_m, dp, dw, pot_r, wfn_r, error, fock_r)
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
            call kernel%get_potential(wfn_x, pot_x)
            call fld%kernel%get_potential(wfn_r, pot_x)
            call get_fock_response(calc, ref%ints, pot_x, fock)

            !>-- derivative of the field Hamiltonian at fixed first-order moments
            do iao = 1, nao
               do jao = 1, nao
                  fock(jao, iao, 1) = fock(jao, iao, 1) + dfock(jao, iao, ic, idir) &
                     & - 0.5_wp*(fld%pot%vao(jao, 1) + fld%pot%vao(iao, 1))*ds(jao, iao, ic)
               end do
            end do
            call add_potential_fock(mol, calc, ref%ints, dfield(idir), ic, fock(:, :, 1))

            call density_response(wfn, fock(:, :, 1), rhs(:, :, 1))
            call density_response_derivative(wfn, fld%dh(:, :, 1), fock_r, ds(:, :, ic), tmp)
            rhs(:, :, 1) = rhs(:, :, 1) + tmp

            call solve_cpscc(mol, calc, wfn, ref%ints, kernel, rhs, dpx, error)
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
