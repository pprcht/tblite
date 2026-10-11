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

!> @file tblite/response/bicgstab.f90
!> Provides an iterative solver for the coupled-perturbed equations

!> Matrix-free stabilized biconjugate gradient solver for linear response
module tblite_response_bicgstab
   use mctc_env, only : wp
   implicit none
   private

   public :: bicgstab_operator, bicgstab_solve

   abstract interface
      !> Application of the linear operator, y = L x
      subroutine bicgstab_operator(x, y)
         import :: wp
         !> Trial vector
         real(wp), intent(in) :: x(:, :, :)
         !> Image of the trial vector
         real(wp), intent(out) :: y(:, :, :)
      end subroutine bicgstab_operator
   end interface

   !> Maximum number of iterations
   integer, parameter :: max_iter = 500
   !> Threshold for a breakdown of the recurrence
   real(wp), parameter :: eps = 1.0e-30_wp

contains

subroutine bicgstab_solve(rhs, lmat, tol, x, resnorm)
   !***********************************
   !* Solve L x = rhs by BiCGSTAB, starting from x = 0. The operator is only
   !* available through its action on a vector.
   !* The relative residual is returned without being checked against the
   !* threshold, the caller has to decide on convergence.
   !***********************************
   !> Right-hand side of the linear system
   real(wp), intent(in), contiguous :: rhs(:, :, :)
   !> Linear operator
   procedure(bicgstab_operator) :: lmat
   !> Convergence threshold for the relative residual
   real(wp), intent(in) :: tol
   !> Solution vector
   real(wp), intent(out), contiguous :: x(:, :, :)
   !> Relative residual norm of the solution
   real(wp), intent(out) :: resnorm

   integer :: iter
   real(wp) :: rho, rho_new, alpha, omega, beta, bnorm, rtv, tt, ts
   real(wp), allocatable :: r(:, :, :), rhat(:, :, :), p(:, :, :), v(:, :, :)
   real(wp), allocatable :: s(:, :, :), t(:, :, :)

   x(:, :, :) = 0.0_wp
   resnorm = 0.0_wp
   bnorm = sqrt(sum(rhs*rhs))
   if (bnorm < eps) return

   allocate(r, rhat, p, v, s, t, mold=x)
   r(:, :, :) = rhs
   rhat(:, :, :) = r
   v(:, :, :) = 0.0_wp
   p(:, :, :) = 0.0_wp
   rho = 1.0_wp
   alpha = 1.0_wp
   omega = 1.0_wp
   resnorm = 1.0_wp

   do iter = 1, max_iter
      rho_new = sum(rhat*r)
      if (abs(rho_new) < eps) exit
      beta = (rho_new/rho) * (alpha/omega)
      p(:, :, :) = r + beta*(p - omega*v)
      call lmat(p, v)
      rtv = sum(rhat*v)
      if (abs(rtv) < eps) exit
      alpha = rho_new/rtv
      s(:, :, :) = r - alpha*v
      resnorm = sqrt(sum(s*s))/bnorm
      if (resnorm < tol) then
         x(:, :, :) = x + alpha*p
         exit
      end if
      call lmat(s, t)
      tt = sum(t*t)
      ts = sum(t*s)
      if (tt < eps) then
         x(:, :, :) = x + alpha*p
         exit
      end if
      omega = ts/tt
      x(:, :, :) = x + alpha*p + omega*s
      r(:, :, :) = s - omega*t
      resnorm = sqrt(sum(r*r))/bnorm
      if (resnorm < tol) exit
      rho = rho_new
   end do
end subroutine bicgstab_solve

end module tblite_response_bicgstab
