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

!> @file tblite/hessian/solvation/shape.f90
!> Provides the second derivative of the ALPB shape correction

!> Second derivative of the shape correction of the analytical linearized
!> Poisson-Boltzmann model
module tblite_hessian_solvation_shape
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: add_shape_hessian

contains

subroutine add_shape_hessian(mol, rad, kepsalpbet, qtotal, hessian)
   !***********************************
   !* Add the Hessian of E = keps alpbet qtot**2 / (2 A_det), with
   !* 1/A_det = sqrt(2/5 V) det(I)**(-1/6) from the inertia tensor I of the
   !* atomic spheres. The total charge does not change with the geometry.
   !* With C = det(I) I**-1 the derivatives of the determinant are
   !*   d det/dx_a       = tr(C dI_a),
   !*   d2 det/dx_a dx_b = tr(C d2I_ab) + det [tr(I**-1 dI_a) tr(I**-1 dI_b)
   !*                                          - tr(I**-1 dI_a I**-1 dI_b)],
   !* where tr(C d2I_ab) = 2 g_ab (tr(C) delta_mn - C_mn) for the Cartesian
   !* directions m, n of a, b and g_ab accounts for the moving center.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic radii of the spheres
   real(wp), intent(in) :: rad(:)
   !> Product of the dielectric screening and the ALPB constant
   real(wp), intent(in) :: kepsalpbet
   !> Total charge
   real(wp), intent(in) :: qtotal
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   real(wp), parameter :: tof = 2.0_wp/5.0_wp
   integer :: nat, iat, jat, ic, jc, kc, lc, ii, jj
   real(wp) :: vol, center(3), r2, inertia(3, 3), iinv(3, 3), det, cmat(3, 3)
   real(wp) :: mmat(3, 3), trc, cvol, kbig, gab, d2det
   real(wp), allocatable :: r3(:), vec(:, :), xd(:, :, :, :), ixd(:, :, :, :)
   real(wp), allocatable :: trcd(:, :), trid(:, :)

   if (kepsalpbet == 0.0_wp .or. qtotal == 0.0_wp) return

   nat = mol%nat
   allocate(r3(nat), vec(3, nat))
   vol = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nat
      r3(iat) = rad(iat)**3
      vol = vol + r3(iat)
      center(:) = center + mol%xyz(:, iat)*r3(iat)
   end do
   center(:) = center / vol

   inertia(:, :) = 0.0_wp
   do iat = 1, nat
      vec(:, iat) = mol%xyz(:, iat) - center
      r2 = sum(vec(:, iat)**2)
      do lc = 1, 3
         do kc = 1, 3
            inertia(kc, lc) = inertia(kc, lc) - r3(iat)*vec(kc, iat)*vec(lc, iat)
         end do
         inertia(lc, lc) = inertia(lc, lc) + r3(iat)*(r2 + tof*rad(iat)**2)
      end do
   end do

   det = matdet_3x3(inertia)
   iinv(:, :) = matinv_3x3(inertia)
   cmat(:, :) = det*iinv
   trc = cmat(1, 1) + cmat(2, 2) + cmat(3, 3)
   !>-- tr(C) delta_mn - C_mn, contracted with the second derivative of I
   mmat(:, :) = -cmat
   do kc = 1, 3
      mmat(kc, kc) = mmat(kc, kc) + trc
   end do

   cvol = sqrt(tof*vol)
   kbig = 0.5_wp*kepsalpbet*qtotal**2*cvol

   !>-- first derivative dI_a of the inertia tensor for every coordinate
   allocate(xd(3, 3, 3, nat), ixd(3, 3, 3, nat), trcd(3, nat), trid(3, nat))
   do iat = 1, nat
      do ic = 1, 3
         xd(:, :, ic, iat) = 0.0_wp
         do kc = 1, 3
            xd(kc, kc, ic, iat) = 2*r3(iat)*vec(ic, iat)
            xd(ic, kc, ic, iat) = xd(ic, kc, ic, iat) - r3(iat)*vec(kc, iat)
            xd(kc, ic, ic, iat) = xd(kc, ic, ic, iat) - r3(iat)*vec(kc, iat)
         end do
         ixd(:, :, ic, iat) = matmul(iinv, xd(:, :, ic, iat))
         !>-- C and dI_a are symmetric, the trace of their product is their dot product
         trcd(ic, iat) = sum(cmat*xd(:, :, ic, iat))
         trid(ic, iat) = ixd(1, 1, ic, iat) + ixd(2, 2, ic, iat) + ixd(3, 3, ic, iat)
      end do
   end do

   do jat = 1, nat
      do iat = 1, nat
         !>-- the center of the spheres moves with every atom
         gab = r3(iat)*(merge(1.0_wp, 0.0_wp, iat == jat) - r3(jat)/vol)
         do jc = 1, 3
            jj = 3*(jat - 1) + jc
            do ic = 1, 3
               ii = 3*(iat - 1) + ic
               d2det = 2*gab*mmat(ic, jc) + det*(trid(ic, iat)*trid(jc, jat) &
                  & - sum(ixd(:, :, ic, iat)*transpose(ixd(:, :, jc, jat))))
               hessian(ii, jj) = hessian(ii, jj) &
                  & + (7*kbig/36.0_wp)*det**(-13.0_wp/6.0_wp)*trcd(ic, iat)*trcd(jc, jat) &
                  & - (kbig/6.0_wp)*det**(-7.0_wp/6.0_wp)*d2det
            end do
         end do
      end do
   end do
end subroutine add_shape_hessian

end module tblite_hessian_solvation_shape
