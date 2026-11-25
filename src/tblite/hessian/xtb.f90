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

!> @file tblite/hessian/xtb.f90
!> Provides the analytic molecular Hessian

!> Analytic second derivatives of the energy of a self-consistent xTB
!> calculation w.r.t. the atomic positions
module tblite_hessian_xtb
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_blas, only : gemv
   use tblite_hessian_hamiltonian, only : get_hamiltonian_hessian
   use tblite_hessian_multipole, only : get_multipole_band_hessian
   use tblite_hessian_reference, only : reference_data, get_nuclear_perturbation
   use tblite_response_kernel, only : interaction_cache, get_response_kernel, &
      & get_response_derivative, get_interaction_hessian
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_type, only : kernel_type, response_type, new_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_hamiltonian_gradient
   implicit none
   private

   public :: xtb_hessian

contains

subroutine xtb_hessian(mol, calc, wfn, ref, cache, kernel, wdensity, hessian, error, dipgrad)
   !***********************************
   !* Analytic Hessian, ordered as hessian(3*(jat-1)+jc, 3*(iat-1)+ic).
   !* Each column is the derivative of the gradient w.r.t. one coordinate
   !* and consists of two parts. The response part contracts the first-order
   !* wavefunction of the coupled-perturbed equations like a gradient, the
   !* explicit part is the second derivative of all terms at fixed density
   !* and is added for the complete matrix at once.
   !* The derivative of the dipole moment w.r.t. the atomic positions is
   !* available from the same first-order wavefunctions.
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
   !> Energy-weighted density matrix
   real(wp), intent(in) :: wdensity(:, :, :)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), contiguous, intent(out) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Derivative of the dipole moment, dipgrad(:, ic, iat) for the
   !> displacement of atom iat in direction ic
   real(wp), intent(out), optional :: dipgrad(:, :, :)

   integer :: nat, iat
   logical :: dipole
   real(wp), allocatable :: dmu(:, :, :)
   type(error_type), allocatable :: local_error

   nat = mol%nat
   dipole = present(dipgrad)
   hessian(:, :) = 0.0_wp
   if (dipole) dipgrad(:, :, :) = 0.0_wp

   !>-- response part, one coupled-perturbed solution per coordinate,
   !>   the caches are only read from here on
   allocate(dmu(3, 3, nat), source=0.0_wp)
   !$omp parallel do schedule(dynamic) default(none) &
   !$omp shared(mol, calc, wfn, ref, cache, kernel, nat, dipole, hessian, dmu, error) &
   !$omp private(iat, local_error)
   do iat = 1, nat
      call get_response_columns(mol, calc, wfn, ref, cache, kernel, iat, dipole, &
         & hessian(:, 3*iat-2:3*iat), dmu(:, :, iat), local_error)
      if (allocated(local_error)) then
         !$omp critical (tblite_hessian_xtb_error)
         if (.not.allocated(error)) call move_alloc(local_error, error)
         !$omp end critical (tblite_hessian_xtb_error)
      end if
   end do
   if (allocated(error)) return
   if (dipole) dipgrad(:, :, :) = dmu

   hessian(:, :) = 0.5_wp * (hessian + transpose(hessian))

   !>-- explicit part, second derivatives at fixed density
   call get_hamiltonian_hessian(mol, calc%bas, calc%h0, calc%ncoord, ref%selfenergy, &
      & ref%dsedcn, wfn%density(:, :, 1), wdensity(:, :, 1), ref%pot%vao(:, 1), &
      & ref%dcndr, ref%lattr, ref%cutoff, hessian)
   if (allocated(calc%coulomb%aes2)) then
      call get_multipole_band_hessian(mol, calc%bas, wfn%density(:, :, 1), &
         & ref%pot%vdp(:, :, 1), ref%pot%vqp(:, :, 1), ref%lattr, ref%cutoff, hessian)
   end if
   call get_interaction_hessian(mol, calc, cache, wfn, hessian, error)

end subroutine xtb_hessian

subroutine get_response_columns(mol, calc, wfn, ref, cache, kernel, iat, dipole, columns, &
      & dmu, error)
   !***********************************
   !* Response part of the three Hessian columns belonging to the
   !* displacements of one atom and, on request, the derivative of the
   !* dipole moment. The derivative of the potential at fixed moments is
   !* contracted with the first-order moments of the density response and
   !* twice with the ones at fixed density. The symmetrization of the
   !* Hessian turns the latter into the contraction of both orders.
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
   !> Displaced atom
   integer, intent(in) :: iat
   !> Evaluate the derivative of the dipole moment
   logical, intent(in) :: dipole
   !> Response part of the Hessian columns
   real(wp), intent(out) :: columns(:, :)
   !> Derivative of the dipole moment
   real(wp), intent(inout) :: dmu(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh, ic
   logical :: multipoles
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: ds(:, :, :), fock0(:, :, :), dp(:, :, :), dw(:, :, :)
   real(wp), allocatable :: zero(:, :, :), dqshdr(:, :)
   real(wp), allocatable :: ddpat(:, :, :), dqpat(:, :, :), dpat(:, :, :)
   real(wp), allocatable :: grad(:, :), dEdcn(:)
   type(potential_type) :: pot_r
   type(wavefunction_type) :: wfn_r, wfn_m
   type(response_type) :: response

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)

   call new_potential(pot_r, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   call new_response(response, mol, calc%bas, [3*iat-2, 3*iat-1, 3*iat])
   allocate(ds(nao, nao, 3), fock0(nao, nao, 3), dp(nao, nao, 1), dw(nao, nao, 1))
   allocate(grad(3, nat), dEdcn(nat), dpat(3, nat, 1), ddpat(3, nat, 3), dqpat(6, nat, 3), &
      & dqshdr(nsh, 3))
   allocate(zero(nao, nao, 1), source=0.0_wp)

   call get_nuclear_perturbation(mol, calc, wfn, ref, cache, iat, dipole, ds, fock0, &
      & dqshdr, ddpat, dqpat, error)
   if (allocated(error)) return

   do ic = 1, 3
      !>-- derivative of the charges and moments at fixed density
      wfn_m%qsh(:, 1) = dqshdr(:, ic)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if

      call xtb_nuclear_response(mol, calc, wfn, ref%ints, kernel, ds(:, :, ic), &
         & fock0(:, :, ic), wfn_m, dp, dw(:, :, 1), pot_r, wfn_r, error)
      if (allocated(error)) return

      !>-- first-order density with the unperturbed potential and
      !>   unperturbed density with the first-order potential
      grad(:, :) = -ref%gradient
      dEdcn(:) = -ref%dEdcn
      call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
         & ref%selfenergy, ref%dsedcn, ref%pot, dp, dw, dEdcn, grad, sigma)
      call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
         & ref%selfenergy, ref%dsedcn, pot_r, wfn%density, zero, dEdcn, grad, sigma)
      call gemv(ref%dcndr, dEdcn, grad, beta=1.0_wp)
      columns(:, ic) = reshape(grad, [3*nat])

      !>-- moments contracted with the derivative of the potential, wfn_r
      !>   holds the total first-order moments
      response%qsh(:, ic) = wfn_r%qsh(:, 1) + wfn_m%qsh(:, 1)
      response%qat(:, ic) = wfn_r%qat(:, 1) + wfn_m%qat(:, 1)
      response%dpat(:, :, ic) = wfn_r%dpat(:, :, 1) + wfn_m%dpat(:, :, 1)
      response%qpat(:, :, ic) = wfn_r%qpat(:, :, 1) + wfn_m%qpat(:, :, 1)

      !>-- mu = sum_A (R_A q_A + d_A), the atomic dipoles enter for all methods
      if (dipole) then
         call get_mulliken_atomic_multipoles(calc%bas, ref%ints%dipole, dp, dpat)
         dmu(:, ic) = matmul(mol%xyz, wfn_r%qat(:, 1)) &
            & + sum(dpat(:, :, 1) + ddpat(:, :, ic), dim=2)
         dmu(ic, ic) = dmu(ic, ic) + wfn%qat(iat, 1)
      end if
   end do

   call get_response_derivative(mol, calc, cache, wfn, response, error)
   if (allocated(error)) return
   columns(:, :) = columns + reshape(response%hessian, [3*nat, 3])
end subroutine get_response_columns

end module tblite_hessian_xtb
