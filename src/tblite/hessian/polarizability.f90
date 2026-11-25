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
!> calculation w.r.t. the atomic positions from first-order responses only.
!> The mixed second-order density matrix of a field and a displacement only
!> enters via the dipole moment. By the interchange theorem its contraction
!> with the dipole operator equals the contraction of the right-hand side of
!> its coupled-perturbed equations with the first-order field response,
!> mu_a(dp_bx) = -Tr(dp_a F_bx) - Tr(dh_a T_bx), where F_bx is the mixed
!> perturbation of the Hamiltonian and T_bx the derivative of the uncoupled
!> response. Only the field responses and the nuclear responses of the
!> Hessian are therefore required.
module tblite_hessian_polarizability
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_hessian_reference, only : reference_data, add_potential_fock
   use tblite_response_cpscc, only : get_fock_response
   use tblite_response_density, only : density_response_derivative
   use tblite_response_integral, only : get_multipole_derivative
   use tblite_response_kernel, only : interaction_cache, get_response_kernel, &
      & get_response_derivative
   use tblite_response_polarizability, only : field_data
   use tblite_response_type, only : kernel_type, response_type, new_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: mixed_data, get_field_kernels, new_mixed_data, get_polarizability_derivative

   !> Mixed perturbation by the fields and the displacements of one atom at
   !> fixed first-order moments, the last dimension is the field direction
   type :: mixed_data
      !> Derivative of the multipole part of the first-order field Hamiltonian
      !> w.r.t. the position, dfock(:, :, ic, idir)
      real(wp), allocatable :: dfock(:, :, :, :)
      !> Derivative of the atomic dipole moments of the field response
      real(wp), allocatable :: ddpat(:, :, :, :)
      !> Derivative of the atomic quadrupole moments of the field response
      real(wp), allocatable :: dqpat(:, :, :, :)
      !> Derivative of the first-order field potential at fixed moments
      type(response_type) :: dfield(3)
   end type mixed_data

contains

subroutine get_field_kernels(mol, calc, cache, wfn, field, error)
   !***********************************
   !* Derivative of the response kernel along each field response, which
   !* couples the field and the nuclear first-order moments.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Caches of the interaction containers
   type(interaction_cache), intent(inout) :: cache
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Field response for the three directions
   type(field_data), intent(inout) :: field(3)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: idir

   do idir = 1, 3
      call get_response_kernel(mol, calc, cache, wfn, field(idir)%kernel, error, &
         & dwfn=field(idir)%wfn)
      if (allocated(error)) return
   end do
end subroutine get_field_kernels

subroutine new_mixed_data(mol, calc, wfn, ref, cache, field, iat, mixed, error)
   !***********************************
   !* Derivative of the first-order field Hamiltonian and moments w.r.t. the
   !* displacements of one atom at fixed first-order density and moments.
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
   !> Field response for the three directions
   type(field_data), intent(in) :: field(3)
   !> Displaced atom
   integer, intent(in) :: iat
   !> Mixed perturbation
   type(mixed_data), intent(out) :: mixed
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: idir

   allocate(mixed%dfock(calc%bas%nao, calc%bas%nao, 3, 3), source=0.0_wp)
   allocate(mixed%ddpat(3, mol%nat, 3, 3), mixed%dqpat(6, mol%nat, 3, 3))
   do idir = 1, 3
      call get_multipole_derivative(mol, calc%bas, ref%lattr, ref%cutoff, iat, &
         & field(idir)%dp(:, :, 1), field(idir)%pot%vdp(:, :, 1), &
         & field(idir)%pot%vqp(:, :, 1), mixed%dfock(:, :, :, idir), &
         & mixed%ddpat(:, :, :, idir), mixed%dqpat(:, :, :, idir))

      !>-- the potential of a unit field on an atom is minus its position
      call new_response(mixed%dfield(idir), mol, calc%bas, [3*iat-2, 3*iat-1, 3*iat])
      call get_response_derivative(mol, calc, cache, wfn, mixed%dfield(idir), error, &
         & dwfn=field(idir)%wfn)
      if (allocated(error)) return
      mixed%dfield(idir)%vat(iat, idir) = mixed%dfield(idir)%vat(iat, idir) - 1.0_wp
   end do
end subroutine new_mixed_data

subroutine get_polarizability_derivative(mol, calc, wfn, ref, kernel, field, mixed, iat, &
      & ic, ds, fock_r, wfn_r, dalpha)
   !***********************************
   !* Derivative of the polarizability w.r.t. the displacement of atom iat
   !* in direction ic, dalpha(a, b) for the fields a and b.
   !* For the field b the mixed perturbation of the Hamiltonian contains the
   !* derivative of the first-order field Hamiltonian at fixed moments, the
   !* potential of the change of the field moments due to the derivative of
   !* the integrals and the derivative of the kernel acting on the nuclear
   !* moments. It is contracted with the field response a, the derivative of
   !* the uncoupled response with the first-order field Hamiltonian a.
   !* Input:  first-order nuclear overlap ds, Hamiltonian fock_r and moments
   !*         wfn_r of the displacement, field responses and mixed data
   !* Output: dalpha(3, 3), symmetric
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Quantities of the unperturbed system
   type(reference_data), intent(in) :: ref
   !> Response kernel at the converged state
   type(kernel_type), intent(in) :: kernel
   !> Field response for the three directions
   type(field_data), intent(in) :: field(3)
   !> Mixed perturbation of the displaced atom
   type(mixed_data), intent(in) :: mixed
   !> Displaced atom
   integer, intent(in) :: iat
   !> Displacement direction
   integer, intent(in) :: ic
   !> Derivative of the overlap matrix
   real(wp), intent(in) :: ds(:, :)
   !> First-order Hamiltonian of the displacement
   real(wp), intent(in) :: fock_r(:, :)
   !> Total first-order charges and moments of the displacement
   type(wavefunction_type), intent(in) :: wfn_r
   !> Derivative of the polarizability
   real(wp), intent(out) :: dalpha(:, :)

   integer :: nat, nao, nsh, idir, ia, iao, jao
   logical :: multipoles
   real(wp), allocatable :: fock(:, :, :), tmp(:, :), n0sh(:)
   type(potential_type) :: pot_x
   type(wavefunction_type) :: wfn_x

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)

   call new_potential(pot_x, mol, calc%bas, 1)
   call new_wavefunction(wfn_x, nat, nsh, nao, 1, wfn%kt)
   allocate(fock(nao, nao, 1), tmp(nao, nao))
   allocate(n0sh(nsh), source=0.0_wp)

   do idir = 1, 3
      associate(fld => field(idir))
         !>-- change of the field moments due to the derivative of the integrals
         call get_mulliken_shell_charges(calc%bas, ds, fld%dp, n0sh, wfn_x%qsh)
         call get_qat_from_qsh(calc%bas, wfn_x%qsh, wfn_x%qat)
         if (multipoles) then
            wfn_x%dpat(:, :, 1) = mixed%ddpat(:, :, ic, idir)
            wfn_x%qpat(:, :, 1) = mixed%dqpat(:, :, ic, idir)
         end if

         call pot_x%reset
         call kernel%get_potential(wfn_x, pot_x)
         call fld%kernel%get_potential(wfn_r, pot_x)
         call get_fock_response(calc, ref%ints, pot_x, fock)

         !>-- derivative of the field Hamiltonian at fixed first-order moments
         do iao = 1, nao
            do jao = 1, nao
               fock(jao, iao, 1) = fock(jao, iao, 1) + mixed%dfock(jao, iao, ic, idir) &
                  & - 0.5_wp*(fld%pot%vao(jao, 1) + fld%pot%vao(iao, 1))*ds(jao, iao)
            end do
         end do
         call add_potential_fock(mol, calc, ref%ints, mixed%dfield(idir), ic, fock(:, :, 1))

         call density_response_derivative(wfn, fld%dh(:, :, 1), fock_r, ds, tmp)

         !>-- mu = sum_A (R_A q_A + d_A), the second-order density by the
         !>   interchange theorem, the position enters explicitly
         do ia = 1, 3
            dalpha(ia, idir) = -sum(field(ia)%dp(:, :, 1)*fock(:, :, 1)) &
               & - sum(field(ia)%dh(:, :, 1)*tmp)
         end do
         dalpha(:, idir) = dalpha(:, idir) + matmul(mol%xyz, wfn_x%qat(:, 1)) &
            & + sum(mixed%ddpat(:, :, ic, idir), dim=2)
         dalpha(ic, idir) = dalpha(ic, idir) + fld%wfn%qat(iat, 1)
      end associate
   end do

   dalpha(:, :) = 0.5_wp * (dalpha + transpose(dalpha))
end subroutine get_polarizability_derivative

end module tblite_hessian_polarizability
