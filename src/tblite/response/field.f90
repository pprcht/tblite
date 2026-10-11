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

!> @file tblite/response/field.f90
!> Provides the coupled-perturbed response to a homogeneous electric field

!> First-order response of a self-consistent xTB wavefunction to a
!> homogeneous electric field
module tblite_response_field
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_external_field, only : electric_field, new_electric_field
   use tblite_integral_type, only : integral_type
   use tblite_response_bicgstab, only : bicgstab_solve
   use tblite_response_density, only : density_response
   use tblite_response_kernel, only : get_kernel_potential
   use tblite_scf, only : potential_type, new_potential
   use tblite_scf_iterator, only : get_qat_from_qsh
   use tblite_scf_potential, only : add_pot_to_h1
   use tblite_wavefunction_mulliken, only : get_mulliken_shell_charges, &
      & get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_field_response

   !> Convergence threshold for the relative residual of the response equations
   real(wp), parameter :: rtol = 1.0e-9_wp

contains

subroutine xtb_field_response(mol, calc, wfn, ints, ccache, dcache, idir, dp, dh, pot, &
      & wfn_r, error)
   !***********************************
   !* Solve the coupled-perturbed equations (1 - A) dp = b for a unit field
   !* along the Cartesian direction idir. Here b is the uncoupled density
   !* response to the field and A maps a density to the uncoupled response
   !* to the potential it generates.
   !* On exit wfn_r holds the first-order density and moments, pot the
   !* first-order potential including the field and dh the corresponding
   !* perturbation of the Hamiltonian.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Integrals at the current geometry
   type(integral_type), intent(in) :: ints
   !> Cache of the electrostatic container
   type(container_cache), intent(inout) :: ccache
   !> Cache of the dispersion container
   type(container_cache), intent(inout) :: dcache
   !> Cartesian direction of the electric field
   integer, intent(in) :: idir
   !> First-order density matrix
   real(wp), contiguous, intent(out) :: dp(:, :, :)
   !> First-order Hamiltonian
   real(wp), contiguous, intent(out) :: dh(:, :, :)
   !> First-order potential
   type(potential_type), intent(inout) :: pot
   !> First-order charges and multipole moments
   type(wavefunction_type), intent(inout) :: wfn_r
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   logical :: multipoles
   real(wp) :: resnorm, efield(3)
   real(wp), allocatable :: rhs(:, :, :), dhk(:, :, :), n0sh(:)
   type(potential_type) :: potk
   type(container_cache) :: fcache
   type(electric_field) :: field
   type(wavefunction_type) :: wfn_k
   type(error_type), allocatable :: kernel_error

   multipoles = allocated(calc%coulomb%aes2)

   call new_potential(potk, mol, calc%bas, wfn%nspin)
   call new_wavefunction(wfn_k, mol%nat, calc%bas%nsh, calc%bas%nao, wfn%nspin, wfn%kt)
   allocate(rhs(calc%bas%nao, calc%bas%nao, 1), dhk(calc%bas%nao, calc%bas%nao, 1))
   allocate(n0sh(calc%bas%nsh), source=0.0_wp)

   efield(:) = 0.0_wp
   efield(idir) = 1.0_wp
   call new_electric_field(field, efield)

   call potk%reset
   call field%get_potential(mol, fcache, wfn, potk)
   call get_fock_response(potk, dhk)
   call density_response(wfn, dhk(:, :, 1), rhs(:, :, 1))

   call bicgstab_solve(rhs, apply_operator, rtol, dp, resnorm)
   if (allocated(kernel_error)) then
      call move_alloc(kernel_error, error)
      return
   end if
   !>-- written to fail for a NaN residual as well
   if (.not.(resnorm < rtol)) then
      call fatal_error(error, "Coupled-perturbed field response did not converge")
      return
   end if

   call get_moments(dp, wfn_r)
   wfn_r%density(:, :, 1) = dp(:, :, 1)
   call pot%reset
   call field%get_potential(mol, fcache, wfn, pot)
   call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_r, pot, error)
   if (allocated(error)) return
   call get_fock_response(pot, dh)

contains

   subroutine get_fock_response(vpot, fock)
      !***********************************
      !* Perturbation of the Hamiltonian generated by a potential.
      !***********************************
      !> Potential
      type(potential_type), intent(inout) :: vpot
      !> Perturbation of the Hamiltonian
      real(wp), intent(out) :: fock(:, :, :)

      call add_pot_to_h1(calc%bas, ints, vpot, fock)
      fock(:, :, 1) = fock(:, :, 1) - ints%hamiltonian
   end subroutine get_fock_response

   subroutine get_moments(pmat, wfn_m)
      !***********************************
      !* Mulliken charges and atomic multipole moments of a first-order
      !* density matrix.
      !***********************************
      !> First-order density matrix
      real(wp), intent(in) :: pmat(:, :, :)
      !> Wavefunction to store the moments
      type(wavefunction_type), intent(inout) :: wfn_m

      call get_mulliken_shell_charges(calc%bas, ints%overlap, pmat, n0sh, wfn_m%qsh)
      call get_qat_from_qsh(calc%bas, wfn_m%qsh, wfn_m%qat)
      if (multipoles) then
         call get_mulliken_atomic_multipoles(calc%bas, ints%dipole, pmat, wfn_m%dpat)
         call get_mulliken_atomic_multipoles(calc%bas, ints%quadrupole, pmat, wfn_m%qpat)
      else
         wfn_m%dpat(:, :, :) = 0.0_wp
         wfn_m%qpat(:, :, :) = 0.0_wp
      end if
   end subroutine get_moments

   subroutine apply_operator(x, y)
      !***********************************
      !* Action of the response operator, y = x - A x.
      !***********************************
      !> Trial vector
      real(wp), intent(in) :: x(:, :, :)
      !> Image of the trial vector
      real(wp), intent(out) :: y(:, :, :)

      type(error_type), allocatable :: local_error

      call get_moments(x, wfn_k)
      call potk%reset
      call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_k, potk, local_error)
      if (allocated(local_error) .and. .not.allocated(kernel_error)) then
         call move_alloc(local_error, kernel_error)
      end if
      call get_fock_response(potk, dhk)
      call density_response(wfn, dhk(:, :, 1), y(:, :, 1))
      y(:, :, :) = x - y
   end subroutine apply_operator

end subroutine xtb_field_response

end module tblite_response_field
