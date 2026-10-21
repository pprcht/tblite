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

!> @file tblite/hessian/solvation/coupling.f90
!> Provides the Born radii dependent second derivatives of the Born energy

!> Second derivatives of the generalized Born energy which involve the
!> geometry dependence of the Born radii
module tblite_hessian_solvation_coupling
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_blas, only : gemm
   use tblite_hessian_solvation_kernel, only : get_kernel_derivs, get_radii_combination
   implicit none
   private

   public :: add_born_coupling_hessian

contains

subroutine add_born_coupling_hessian(mol, kernel, keps, qat, rad, draddr, hessian, dEda)
   !***********************************
   !* Add the parts of the Hessian of E = 1/2 q^T J(R, a(R)) q containing the
   !* first derivatives of the Born radii a, with G the Jacobian of the radii:
   !*   G Theta G^T  with Theta = d2E/da da, including the self terms,
   !*   d2E/dr da (dr/dR G^T + G dr/dR^T) for every pair distance r.
   !* The remaining radii term dE/da d2a/dR dR needs dEda, returned here.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel
   !> Dielectric screening
   real(wp), intent(in) :: keps
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Born radii
   real(wp), intent(in) :: rad(:)
   !> Derivative of the Born radii w.r.t. the atomic positions
   real(wp), intent(in) :: draddr(:, :, :)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Derivative of the energy w.r.t. the Born radii
   real(wp), intent(out) :: dEda(:)

   integer :: nat, iat, jat, ic, ii, jj
   real(wp) :: vec(3), ehat(3), r1, wq, rho, da, db, daa, dab, dbb
   real(wp) :: dr, drho, drr, drrho, drhorho
   real(wp), allocatable :: gmat(:, :), theta(:, :), tmp(:, :), hvec(:)

   nat = mol%nat
   allocate(gmat(3*nat, nat), theta(nat, nat), source=0.0_wp)
   allocate(tmp(3*nat, nat), hvec(3*nat))
   gmat(:, :) = reshape(draddr, [3*nat, nat])
   dEda(:nat) = 0.0_wp

   do iat = 1, nat
      ii = 3*(iat - 1)
      do jat = 1, iat - 1
         jj = 3*(jat - 1)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         ehat(:) = vec / r1
         wq = keps*qat(iat)*qat(jat)

         call get_radii_combination(kernel, rad(iat), rad(jat), rho, da, db, daa, dab, dbb)
         call get_kernel_derivs(kernel, r1, rho, dr, drho, drr, drrho, drhorho)

         dEda(iat) = dEda(iat) + wq*drho*da
         dEda(jat) = dEda(jat) + wq*drho*db

         theta(iat, iat) = theta(iat, iat) + wq*(drhorho*da*da + drho*daa)
         theta(jat, jat) = theta(jat, jat) + wq*(drhorho*db*db + drho*dbb)
         theta(iat, jat) = theta(iat, jat) + wq*(drhorho*da*db + drho*dab)
         theta(jat, iat) = theta(jat, iat) + wq*(drhorho*da*db + drho*dab)

         !>-- the distance moves with +ehat on iat and -ehat on jat
         hvec(:) = wq*drrho*(da*gmat(:, iat) + db*gmat(:, jat))
         do ic = 1, 3
            hessian(ii+ic, :) = hessian(ii+ic, :) + ehat(ic)*hvec
            hessian(jj+ic, :) = hessian(jj+ic, :) - ehat(ic)*hvec
            hessian(:, ii+ic) = hessian(:, ii+ic) + hvec*ehat(ic)
            hessian(:, jj+ic) = hessian(:, jj+ic) - hvec*ehat(ic)
         end do
      end do
   end do

   !>-- self terms keps q_A**2 / (2 a_A)
   do iat = 1, nat
      theta(iat, iat) = theta(iat, iat) + keps*qat(iat)**2/rad(iat)**3
      dEda(iat) = dEda(iat) - 0.5_wp*keps*qat(iat)**2/rad(iat)**2
   end do

   call gemm(gmat, theta, tmp)
   call gemm(tmp, gmat, hessian, transb="t", beta=1.0_wp)
end subroutine add_born_coupling_hessian

end module tblite_hessian_solvation_coupling
