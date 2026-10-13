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

!> @file tblite/hessian/pair.f90
!> Provides helpers to assemble Hessians of pairwise interactions

!> Cartesian second derivatives of pairwise radial energy contributions
module tblite_hessian_pair
   use mctc_env, only : wp
   implicit none
   private

   public :: get_pair_block, add_pair_hessian, add_hessian_block

contains

pure subroutine get_pair_block(df, d2f, r1, vec, block)
   !***********************************
   !* Second derivative of a radial function f(r) w.r.t. the position of the
   !* first atom, d2f e e^T + df/r (1 - e e^T) with e the unit vector along
   !* the distance vector.
   !***********************************
   !> First derivative of the function w.r.t. the distance
   real(wp), intent(in) :: df
   !> Second derivative of the function w.r.t. the distance
   real(wp), intent(in) :: d2f
   !> Distance between the two atoms
   real(wp), intent(in) :: r1
   !> Distance vector between the two atoms
   real(wp), intent(in) :: vec(3)
   !> Cartesian second derivative
   real(wp), intent(out) :: block(3, 3)

   integer :: ic
   real(wp) :: ehat(3)

   ehat(:) = vec / r1
   block(:, :) = (d2f - df/r1) * spread(ehat, 1, 3) * spread(ehat, 2, 3)
   do ic = 1, 3
      block(ic, ic) = block(ic, ic) + df/r1
   end do
end subroutine get_pair_block

pure subroutine add_pair_hessian(hessian, iat, jat, block)
   !***********************************
   !* Add the second derivative of a pair energy to the Hessian, with the
   !* block entering the diagonal positively and the off-diagonal negatively.
   !***********************************
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Atoms forming the pair
   integer, intent(in) :: iat, jat
   !> Cartesian second derivative w.r.t. the position of the first atom
   real(wp), intent(in) :: block(3, 3)

   integer :: ii, jj

   ii = 3*(iat - 1)
   jj = 3*(jat - 1)
   hessian(ii+1:ii+3, ii+1:ii+3) = hessian(ii+1:ii+3, ii+1:ii+3) + block
   hessian(jj+1:jj+3, jj+1:jj+3) = hessian(jj+1:jj+3, jj+1:jj+3) + block
   hessian(ii+1:ii+3, jj+1:jj+3) = hessian(ii+1:ii+3, jj+1:jj+3) - block
   hessian(jj+1:jj+3, ii+1:ii+3) = hessian(jj+1:jj+3, ii+1:ii+3) - block
end subroutine add_pair_hessian

pure subroutine add_hessian_block(hessian, iat, jat, block)
   !***********************************
   !* Add a Cartesian block to the entry of the Hessian belonging to the
   !* atoms iat and jat.
   !***********************************
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Atoms addressing the block
   integer, intent(in) :: iat, jat
   !> Cartesian block
   real(wp), intent(in) :: block(3, 3)

   integer :: ii, jj

   ii = 3*(iat - 1)
   jj = 3*(jat - 1)
   hessian(ii+1:ii+3, jj+1:jj+3) = hessian(ii+1:ii+3, jj+1:jj+3) + block
end subroutine add_hessian_block

end module tblite_hessian_pair
