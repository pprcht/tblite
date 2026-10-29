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

!> @file tblite/response/polarizability.f90
!> Provides the analytic static dipole polarizability

!> Static dipole polarizability of a self-consistent xTB wavefunction
module tblite_response_polarizability
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_basis_type, only : get_cutoff
   use tblite_container, only : container_cache
   use tblite_context, only : context_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_solvation, only : solvation_kernel, new_solvation_kernel
   use tblite_integral_type, only : integral_type, new_integral
   use tblite_partition, only : serial_work_partition, operator(==)
   use tblite_response_field, only : xtb_field_response
   use tblite_scf, only : potential_type, new_potential
   use tblite_wavefunction_mulliken, only : get_mulliken_atomic_multipoles
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_h0, only : get_selfenergy, get_hamiltonian
   implicit none
   private

   public :: xtb_polarizability, check_response_support

contains

subroutine xtb_polarizability(ctx, mol, calc, wfn, alpha)
   !***********************************
   !* Static dipole polarizability alpha(a, b) = dmu_a/dF_b from one
   !* coupled-perturbed field response per Cartesian direction.
   !* Requires the converged wavefunction of a preceding single point
   !* calculation with the same calculator. Failures are reported in ctx.
   !***********************************
   !> Calculation context
   type(context_type), intent(inout) :: ctx
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Polarizability tensor
   real(wp), intent(out) :: alpha(3, 3)

   integer :: idir
   real(wp) :: cutoff
   real(wp), allocatable :: cn(:), lattr(:, :), selfenergy(:)
   real(wp), allocatable :: dh(:, :, :), dp(:, :, :), ddpat(:, :, :)
   type(integral_type) :: ints
   type(container_cache) :: ccache, dcache
   type(adjacency_list) :: list
   type(potential_type) :: pot
   type(wavefunction_type) :: wfn_r
   type(solvation_kernel) :: solv
   type(error_type), allocatable :: error

   alpha(:, :) = 0.0_wp

   call check_response_support(mol, calc, wfn, error)
   !>-- implicit solvation enters the field response only via its charge kernel
   if (.not.allocated(error) .and. allocated(calc%interactions)) then
      call new_solvation_kernel(solv, mol, calc, error)
   end if
   if (allocated(error)) then
      call ctx%set_error(error)
      return
   end if

   if (allocated(calc%ncoord)) then
      allocate(cn(mol%nat))
      call calc%ncoord%get_cn(mol, cn)
   end if
   allocate(selfenergy(calc%bas%nsh))
   call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, cn=cn, &
      & selfenergy=selfenergy)

   cutoff = get_cutoff(calc%bas)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call new_adjacency_list(list, mol, lattr, cutoff)
   call new_integral(ints, calc%bas%nao)
   call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
      & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian)

   call calc%coulomb%update(mol, ccache)
   if (allocated(calc%dispersion)) call calc%dispersion%update(mol, dcache)

   call new_potential(pot, mol, calc%bas, wfn%nspin)
   call new_wavefunction(wfn_r, mol%nat, calc%bas%nsh, calc%bas%nao, wfn%nspin, wfn%kt)
   allocate(dh(calc%bas%nao, calc%bas%nao, 1), dp(calc%bas%nao, calc%bas%nao, 1), &
      & ddpat(3, mol%nat, 1))

   do idir = 1, 3
      call xtb_field_response(mol, calc, wfn, ints, ccache, dcache, idir, dp, dh, pot, &
         & wfn_r, error, kat=solv%kmat)
      if (allocated(error)) then
         call ctx%set_error(error)
         return
      end if

      !>-- mu = sum_A (R_A q_A + d_A), the atomic dipoles enter for all methods
      call get_mulliken_atomic_multipoles(calc%bas, ints%dipole, dp, ddpat)
      alpha(:, idir) = matmul(mol%xyz, wfn_r%qat(:, 1)) + sum(ddpat(:, :, 1), dim=2)
   end do

end subroutine xtb_polarizability

subroutine check_response_support(mol, calc, wfn, error)
   !***********************************
   !* Reject calculations outside the validity of the analytic response:
   !* open-shell, non-self-consistent, periodic and distributed. Interaction
   !* containers are checked when collecting the implicit solvation models.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (wfn%nspin /= 1) then
      call fatal_error(error, "Analytic response requires a restricted wavefunction")
      return
   end if
   if (.not.allocated(calc%coulomb)) then
      call fatal_error(error, "Analytic response requires a self-consistent method")
      return
   end if
   if (any(mol%periodic)) then
      call fatal_error(error, "Analytic response is not available for periodic systems")
      return
   end if
   if (.not.(calc%partition == serial_work_partition)) then
      call fatal_error(error, "Analytic response is not available with a work partition")
      return
   end if
end subroutine check_response_support

end module tblite_response_polarizability
