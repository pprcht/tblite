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

!> @file tblite/hessian/solvation/cm5.f90
!> Provides the second derivatives of the CM5 charges

!> Second derivatives of the CM5 charge corrections
module tblite_hessian_solvation_cm5
   use mctc_data_atomicrad, only : get_atomic_rad
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_solvation_cm5, only : cm5_alpha, get_cm5_pairpar
   implicit none
   private

   public :: add_cm5_hessian

contains

subroutine add_cm5_hessian(mol, weight, hessian)
   !***********************************
   !* Add sum_A weight_A d2c_A/dR dR for the CM5 corrections
   !* c_A = sum_B p_AB exp(-alpha (r_AB - R_A - R_B)). The pair parameter is
   !* antisymmetric, so both directions of a pair combine into one radial
   !* term weighted by (weight_A - weight_B) p_AB.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Weight of each CM5 correction
   real(wp), intent(in) :: weight(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat, nati, natj
   real(wp) :: vec(3), r1, radi, radj, bkk, wpair, block(3, 3)

   do iat = 1, mol%nat
      nati = mol%num(mol%id(iat))
      radi = get_atomic_rad(nati)
      do jat = 1, iat - 1
         natj = mol%num(mol%id(jat))
         if (nati == natj) cycle
         radj = get_atomic_rad(natj)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         bkk = exp(-cm5_alpha*(r1 - radi - radj))
         wpair = (weight(iat) - weight(jat))*get_cm5_pairpar(nati, natj)
         call get_pair_block(-wpair*cm5_alpha*bkk, wpair*cm5_alpha**2*bkk, r1, vec, block)
         call add_pair_hessian(hessian, iat, jat, block)
      end do
   end do
end subroutine add_cm5_hessian

end module tblite_hessian_solvation_cm5
