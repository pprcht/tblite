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

!> @file tblite/hessian/coulomb.f90
!> Provides derivatives of the electrostatic energy required for the Hessian

!> Second derivatives of the isotropic second-order electrostatic energy
module tblite_hessian_coulomb
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container, only : container_cache
   use tblite_coulomb_charge, only : effective_coulomb
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_scf, only : potential_type, new_potential
   use tblite_wavefunction_type, only : wavefunction_type, new_wavefunction
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: get_coulomb_hessian, get_potential_derivative

contains

subroutine get_coulomb_hessian(calc, mol, wfn, hessian, error)
   !***********************************
   !* Add the Hessian of the second-order electrostatic energy at fixed shell
   !* charges, E = 1/2 sum q_i gamma_ij(r) q_j with the Klopman-Ohno kernel
   !* gamma = (r**g + eta**(-g))**(-1/g). The onsite third-order energy does
   !* not depend on the geometry. Only available for molecular systems.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
   real(wp) :: g, vec(3), r1, u, dgam, d2gam, qq, block(3, 3)

   if (.not.allocated(calc%coulomb%es2)) return

   select type(es2 => calc%coulomb%es2)
   type is (effective_coulomb)
      g = es2%gexp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = es2%offset(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jj = es2%offset(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r1 = norm2(vec)
            do ish = 1, es2%nshell(iat)
               do jsh = 1, es2%nshell(jat)
                  u = r1**g + es2%hubbard(jsh, ish, jzp, izp)**(-g)
                  dgam = -r1**(g - 1.0_wp) * u**(-1.0_wp/g - 1.0_wp)
                  d2gam = -(g - 1.0_wp) * r1**(g - 2.0_wp) * u**(-1.0_wp/g - 1.0_wp) &
                     & + (g + 1.0_wp) * r1**(2*g - 2.0_wp) * u**(-1.0_wp/g - 2.0_wp)
                  qq = wfn%qsh(ii+ish, 1) * wfn%qsh(jj+jsh, 1)
                  call get_pair_block(qq*dgam, qq*d2gam, r1, vec, block)
                  call add_pair_hessian(hessian, iat, jat, block)
               end do
            end do
         end do
      end do
   class default
      call fatal_error(error, "Hessian is not available for this electrostatic model")
   end select
end subroutine get_coulomb_hessian

subroutine get_potential_derivative(mol, calc, ccache, wfn, dvshdr)
   !***********************************
   !* Derivative of the shell-resolved electrostatic potential w.r.t. the
   !* atomic positions at fixed charges, dvshdr(ish, 3*(iat-1)+ic).
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Cache of the electrostatic container
   type(container_cache), intent(inout) :: ccache
   !> Converged wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Derivative of the shell-resolved potential
   real(wp), intent(out) :: dvshdr(:, :)

   integer :: iat, ic, ish
   type(wavefunction_type) :: wfn_q
   type(potential_type) :: pot

   !>-- without derivatives of the charges only the explicit dependence remains
   call new_wavefunction(wfn_q, mol%nat, calc%bas%nsh, calc%bas%nao, 1, wfn%kt, .true.)
   wfn_q%qat(:, 1) = wfn%qat(:, 1)
   wfn_q%qsh(:, 1) = wfn%qsh(:, 1)
   wfn_q%dqatdr(:, :, :, :) = 0.0_wp
   wfn_q%dqshdr(:, :, :, :) = 0.0_wp
   wfn_q%dqatdL(:, :, :, :) = 0.0_wp
   wfn_q%dqshdL(:, :, :, :) = 0.0_wp

   call new_potential(pot, mol, calc%bas, 1, .true.)
   call pot%reset
   call calc%coulomb%get_potential_gradient(mol, ccache, wfn_q, pot)

   do iat = 1, mol%nat
      do ic = 1, 3
         do ish = 1, calc%bas%nsh
            dvshdr(ish, 3*(iat-1)+ic) = pot%dvshdr(ic, iat, ish, 1)
         end do
      end do
   end do
end subroutine get_potential_derivative

end module tblite_hessian_coulomb
