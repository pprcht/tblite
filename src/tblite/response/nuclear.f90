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

!> @file tblite/response/nuclear.f90
!> Provides the coupled-perturbed response to a nuclear displacement

!> First-order response of a self-consistent xTB wavefunction to the
!> displacement of an atom
module tblite_response_nuclear
   use mctc_env, only : wp, error_type
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_integral_type, only : integral_type
   use tblite_response_cpscc, only : solve_cpscc, get_response_moments, &
      & get_fock_response
   use tblite_response_density, only : density_response, wdensity_response
   use tblite_response_kernel, only : get_kernel_potential
   use tblite_scf, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: xtb_nuclear_response

contains

subroutine xtb_nuclear_response(mol, calc, wfn, ints, ccache, dcache, ds, fock0, wfn_m, &
      & dp, dw, pot, wfn_r, error, fock_r)
   !***********************************
   !* Solve the coupled-perturbed equations for one Cartesian displacement.
   !* The perturbation enters via the derivative of the overlap ds and of
   !* the Hamiltonian at fixed charges and multipole moments fock0.
   !* The moments change already at fixed density due to the derivative of
   !* the integrals, given in wfn_m. Their potential is added to the
   !* right-hand side such that the response operator is the one of a
   !* perturbation without overlap dependence.
   !* On exit dp and dw are the first-order density and energy-weighted
   !* density matrix, wfn_r holds the total first-order charges and moments
   !* and pot the potential they generate.
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
   !> Derivative of the overlap matrix
   real(wp), intent(in) :: ds(:, :)
   !> Derivative of the Hamiltonian at fixed charges and multipole moments
   real(wp), intent(in) :: fock0(:, :)
   !> Derivative of the charges and multipole moments at fixed density
   type(wavefunction_type), intent(in) :: wfn_m
   !> First-order density matrix
   real(wp), contiguous, intent(out) :: dp(:, :, :)
   !> First-order energy-weighted density matrix
   real(wp), intent(out) :: dw(:, :)
   !> First-order potential
   type(potential_type), intent(inout) :: pot
   !> First-order charges and multipole moments
   type(wavefunction_type), intent(inout) :: wfn_r
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order Hamiltonian including the response of the potential
   real(wp), intent(out), optional :: fock_r(:, :)

   real(wp), allocatable :: fock(:, :, :), rhs(:, :, :)

   allocate(fock(calc%bas%nao, calc%bas%nao, 1), rhs(calc%bas%nao, calc%bas%nao, 1))

   call pot%reset
   call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_m, pot, error)
   if (allocated(error)) return
   call get_fock_response(calc, ints, pot, fock)
   fock(:, :, 1) = fock(:, :, 1) + fock0
   call density_response(wfn, fock(:, :, 1), rhs(:, :, 1), ds)

   call solve_cpscc(mol, calc, wfn, ints, ccache, dcache, rhs, dp, error)
   if (allocated(error)) return

   call get_response_moments(calc, ints, dp, wfn_r)
   wfn_r%qsh(:, :) = wfn_r%qsh + wfn_m%qsh
   wfn_r%qat(:, :) = wfn_r%qat + wfn_m%qat
   wfn_r%dpat(:, :, :) = wfn_r%dpat + wfn_m%dpat
   wfn_r%qpat(:, :, :) = wfn_r%qpat + wfn_m%qpat

   call pot%reset
   call get_kernel_potential(mol, calc, ccache, dcache, wfn, wfn_r, pot, error)
   if (allocated(error)) return
   call get_fock_response(calc, ints, pot, fock)
   fock(:, :, 1) = fock(:, :, 1) + fock0
   call wdensity_response(wfn, fock(:, :, 1), dw, ds)
   if (present(fock_r)) fock_r(:, :) = fock(:, :, 1)
end subroutine xtb_nuclear_response

end module tblite_response_nuclear
