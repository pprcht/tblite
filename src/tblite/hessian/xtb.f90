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
!> Provides the analytic derivatives w.r.t. the atomic positions

!> Analytic second derivatives of the energy of a self-consistent xTB
!> calculation w.r.t. the atomic positions, and the derivatives of the
!> dipole moment and of the polarizability, from one coupled-perturbed
!> solution per Cartesian displacement
module tblite_hessian_xtb
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_blas, only : gemv
   use tblite_hessian_hamiltonian, only : get_hamiltonian_hessian
   use tblite_hessian_multipole, only : get_multipole_band_hessian
   use tblite_hessian_polarizability, only : mixed_data, new_mixed_data, &
      & get_polarizability_derivative
   use tblite_hessian_reference, only : reference_data, get_nuclear_perturbation
   use tblite_response_kernel, only : interaction_cache, get_response_derivative, &
      & get_interaction_hessian
   use tblite_response_nuclear, only : xtb_nuclear_response
   use tblite_response_polarizability, only : field_data
   use tblite_response_type, only : kernel_type, response_type, new_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_wavefunction_mulliken, only : get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_hamiltonian_gradient
   implicit none
   private

   public :: xtb_nuclear_derivatives

contains

subroutine xtb_nuclear_derivatives(mol, calc, wfn, ref, cache, kernel, field, error, &
      & hessian, wdensity, dipgrad, dalpha)
   !***********************************
   !* Derivatives w.r.t. the atomic positions, each one evaluated if present.
   !* The Hessian is ordered as hessian(3*(jat-1)+jc, 3*(iat-1)+ic). Each
   !* column is the derivative of the gradient w.r.t. one coordinate and
   !* consists of two parts. The response part contracts the first-order
   !* wavefunction of the coupled-perturbed equations like a gradient, the
   !* explicit part is the second derivative of all terms at fixed density
   !* and is added for the complete matrix at once.
   !* The derivatives of the dipole moment, dipgrad(:, ic, iat), and of the
   !* polarizability, dalpha(:, :, ic, iat), are available from the same
   !* first-order wavefunctions, the latter with the field responses.
   !* Input:  wdensity for the Hessian, field responses with their kernel
   !*         derivatives for dalpha, field is not referenced otherwise
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
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), contiguous, intent(out), optional :: hessian(:, :)
   !> Energy-weighted density matrix, required for the Hessian
   real(wp), intent(in), optional :: wdensity(:, :, :)
   !> Derivative of the dipole moment
   real(wp), intent(out), optional :: dipgrad(:, :, :)
   !> Derivative of the polarizability
   real(wp), intent(out), optional :: dalpha(:, :, :, :)

   integer :: nat, iat, nhess, ndip, npol
   logical :: do_hess, do_dip, do_pol
   real(wp), allocatable :: hess(:, :), dmu(:, :, :), dpol(:, :, :, :)
   type(error_type), allocatable :: local_error

   nat = mol%nat
   do_hess = present(hessian)
   do_dip = present(dipgrad)
   do_pol = present(dalpha)
   nhess = merge(3*nat, 0, do_hess)
   ndip = merge(nat, 0, do_dip)
   npol = merge(nat, 0, do_pol)
   allocate(hess(3*nat, nhess), dmu(3, 3, ndip), dpol(3, 3, 3, npol), source=0.0_wp)

   !>-- one coupled-perturbed solution per coordinate, the caches are only
   !>   read from here on
   !$omp parallel do schedule(dynamic) default(none) &
   !$omp shared(mol, calc, wfn, ref, cache, kernel, field, nat, do_hess, do_dip) &
   !$omp shared(do_pol, hess, dmu, dpol, error) private(iat, local_error)
   do iat = 1, nat
      call get_atom_derivatives(mol, calc, wfn, ref, cache, kernel, field, iat, &
         & do_hess, do_dip, do_pol, hess, dmu, dpol, local_error)
      if (allocated(local_error)) then
         !$omp critical (tblite_hessian_xtb_error)
         if (.not.allocated(error)) call move_alloc(local_error, error)
         !$omp end critical (tblite_hessian_xtb_error)
      end if
   end do
   if (allocated(error)) return
   if (do_dip) dipgrad(:, :, :) = dmu
   if (do_pol) dalpha(:, :, :, :) = dpol
   if (.not.do_hess) return

   hessian(:, :) = 0.5_wp * (hess + transpose(hess))

   !>-- explicit part, second derivatives at fixed density
   call get_hamiltonian_hessian(mol, calc%bas, calc%h0, calc%ncoord, ref%selfenergy, &
      & ref%dsedcn, wfn%density(:, :, 1), wdensity(:, :, 1), ref%pot%vao(:, 1), &
      & ref%dcndr, ref%lattr, ref%cutoff, hessian)
   if (allocated(calc%coulomb%aes2)) then
      call get_multipole_band_hessian(mol, calc%bas, wfn%density(:, :, 1), &
         & ref%pot%vdp(:, :, 1), ref%pot%vqp(:, :, 1), ref%lattr, ref%cutoff, hessian)
   end if
   call get_interaction_hessian(mol, calc, cache, wfn, hessian, error)

end subroutine xtb_nuclear_derivatives

subroutine get_atom_derivatives(mol, calc, wfn, ref, cache, kernel, field, iat, do_hess, &
      & do_dip, do_pol, hess, dmu, dpol, error)
   !***********************************
   !* Derivatives w.r.t. the three displacements of one atom from one
   !* coupled-perturbed solution each. The response part of the Hessian
   !* contracts the derivative of the potential at fixed moments with the
   !* first-order moments of the density response and twice with the ones
   !* at fixed density. The symmetrization of the Hessian turns the latter
   !* into the contraction of both orders.
   !* Only the columns of atom iat in hess, dmu and dpol are written, all
   !* work arrays are local to allow evaluating atoms in parallel.
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
   !> Field response for the three directions, only referenced for do_pol
   type(field_data), intent(in) :: field(3)
   !> Displaced atom
   integer, intent(in) :: iat
   !> Evaluate the response part of the Hessian
   logical, intent(in) :: do_hess
   !> Evaluate the derivative of the dipole moment
   logical, intent(in) :: do_dip
   !> Evaluate the derivative of the polarizability
   logical, intent(in) :: do_pol
   !> Response part of the Hessian
   real(wp), intent(inout) :: hess(:, :)
   !> Derivative of the dipole moment
   real(wp), intent(inout) :: dmu(:, :, :)
   !> Derivative of the polarizability
   real(wp), intent(inout) :: dpol(:, :, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, nao, nsh, ic
   logical :: multipoles
   real(wp) :: sigma(3, 3)
   real(wp), allocatable :: ds(:, :, :), fock0(:, :, :), dp(:, :, :), dw(:, :, :)
   real(wp), allocatable :: zero(:, :, :), dqshdr(:, :), fock_r(:, :)
   real(wp), allocatable :: ddpat(:, :, :), dqpat(:, :, :), dpat(:, :, :)
   real(wp), allocatable :: grad(:, :), dEdcn(:)
   type(potential_type) :: pot_r
   type(wavefunction_type) :: wfn_r, wfn_m
   type(response_type) :: response
   type(mixed_data) :: mixed

   nat = mol%nat
   nao = calc%bas%nao
   nsh = calc%bas%nsh
   multipoles = allocated(calc%coulomb%aes2)

   call new_potential(pot_r, mol, calc%bas, 1)
   call new_wavefunction(wfn_r, nat, nsh, nao, 1, wfn%kt)
   call new_wavefunction(wfn_m, nat, nsh, nao, 1, wfn%kt)
   allocate(ds(nao, nao, 3), fock0(nao, nao, 3), dp(nao, nao, 1), dw(nao, nao, 1), &
      & fock_r(nao, nao))
   allocate(ddpat(3, nat, 3), dqpat(6, nat, 3), dqshdr(nsh, 3), dpat(3, nat, 1))
   if (do_hess) then
      call new_response(response, mol, calc%bas, [3*iat-2, 3*iat-1, 3*iat])
      allocate(grad(3, nat), dEdcn(nat))
      allocate(zero(nao, nao, 1), source=0.0_wp)
   end if

   call get_nuclear_perturbation(mol, calc, wfn, ref, cache, iat, do_dip .or. do_pol, &
      & ds, fock0, dqshdr, ddpat, dqpat, error)
   if (allocated(error)) return
   if (do_pol) then
      call new_mixed_data(mol, calc, wfn, ref, cache, field, iat, mixed, error)
      if (allocated(error)) return
   end if

   do ic = 1, 3
      !>-- derivative of the charges and moments at fixed density
      wfn_m%qsh(:, 1) = dqshdr(:, ic)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         wfn_m%dpat(:, :, 1) = ddpat(:, :, ic)
         wfn_m%qpat(:, :, 1) = dqpat(:, :, ic)
      end if

      call xtb_nuclear_response(mol, calc, wfn, ref%ints, kernel, ds(:, :, ic), &
         & fock0(:, :, ic), wfn_m, dp, dw(:, :, 1), pot_r, wfn_r, error, fock_r)
      if (allocated(error)) return

      if (do_hess) then
         !>-- first-order density with the unperturbed potential and
         !>   unperturbed density with the first-order potential
         grad(:, :) = -ref%gradient
         dEdcn(:) = -ref%dEdcn
         call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
            & ref%selfenergy, ref%dsedcn, ref%pot, dp, dw, dEdcn, grad, sigma)
         call get_hamiltonian_gradient(mol, ref%lattr, ref%list, calc%bas, calc%h0, &
            & ref%selfenergy, ref%dsedcn, pot_r, wfn%density, zero, dEdcn, grad, sigma)
         call gemv(ref%dcndr, dEdcn, grad, beta=1.0_wp)
         hess(:, 3*(iat-1)+ic) = reshape(grad, [3*nat])

         !>-- moments contracted with the derivative of the potential, wfn_r
         !>   holds the total first-order moments
         response%qsh(:, ic) = wfn_r%qsh(:, 1) + wfn_m%qsh(:, 1)
         response%qat(:, ic) = wfn_r%qat(:, 1) + wfn_m%qat(:, 1)
         response%dpat(:, :, ic) = wfn_r%dpat(:, :, 1) + wfn_m%dpat(:, :, 1)
         response%qpat(:, :, ic) = wfn_r%qpat(:, :, 1) + wfn_m%qpat(:, :, 1)
      end if

      !>-- mu = sum_A (R_A q_A + d_A), the atomic dipoles enter for all methods
      if (do_dip) then
         call get_mulliken_atomic_multipoles(calc%bas, ref%ints%dipole, dp, dpat)
         dmu(:, ic, iat) = matmul(mol%xyz, wfn_r%qat(:, 1)) &
            & + sum(dpat(:, :, 1) + ddpat(:, :, ic), dim=2)
         dmu(ic, ic, iat) = dmu(ic, ic, iat) + wfn%qat(iat, 1)
      end if

      if (do_pol) then
         call get_polarizability_derivative(mol, calc, wfn, ref, kernel, field, mixed, &
            & iat, ic, ds(:, :, ic), fock_r, wfn_r, dpol(:, :, ic, iat))
      end if
   end do

   if (do_hess) then
      call get_response_derivative(mol, calc, cache, wfn, response, error)
      if (allocated(error)) return
      hess(:, 3*iat-2:3*iat) = hess(:, 3*iat-2:3*iat) + reshape(response%hessian, [3*nat, 3])
   end if
end subroutine get_atom_derivatives

end module tblite_hessian_xtb
