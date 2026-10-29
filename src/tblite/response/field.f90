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
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_external_field, only : electric_field, new_electric_field
   use tblite_integral_type, only : integral_type
   use tblite_response_cpscc, only : solve_cpscc, get_response_moments, &
      & get_fock_response
   use tblite_response_density, only : density_response
   use tblite_response_kernel, only : get_kernel_potential
   use tblite_scf, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_field_response

contains

subroutine xtb_field_response(mol, calc, wfn, ints, ccache, dcache, idir, dp, dh, pot, &
      & wfn_r, error, kat)
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
   !> Atom-resolved charge kernel of additional interactions
   real(wp), intent(in), optional :: kat(:, :)

   real(wp) :: efield(3)
   real(wp), allocatable :: rhs(:, :, :)
   type(container_cache) :: fcache
   type(electric_field) :: field

   allocate(rhs(calc%bas%nao, calc%bas%nao, 1))

   efield(:) = 0.0_wp
   efield(idir) = 1.0_wp
   call new_electric_field(field, efield)

   call pot%reset
   call field%get_potential(mol, fcache, wfn, pot)
   call get_fock_response(calc, ints, pot, dh)
   call density_response(wfn, dh(:, :, 1), rhs(:, :, 1))

   call solve_cpscc(mol, calc, wfn, ints, ccache, dcache, rhs, dp, error, kat)
   if (allocated(error)) return

   call get_response_moments(calc, ints, dp, wfn_r)
   wfn_r%density(:, :, 1) = dp(:, :, 1)
   !>-- folding a potential accumulates, hence it is built from scratch
   call pot%reset
   call field%get_potential(mol, fcache, wfn, pot)
   call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_r, pot, error, kat)
   if (allocated(error)) return
   call get_fock_response(calc, ints, pot, dh)

end subroutine xtb_field_response

end module tblite_response_field
