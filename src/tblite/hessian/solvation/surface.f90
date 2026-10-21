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

!> @file tblite/hessian/solvation/surface.f90
!> Provides the second derivatives of the solvent accessible surface area

!> Second derivatives of the atomic solvent accessible surface areas
module tblite_hessian_solvation_surface
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_solvation_surface, only : surface_integrator, compute_w_spg
   implicit none
   private

   public :: add_surface_hessian

   !> Threshold for grid points considered as buried, identical to the one of the area
   real(wp), parameter :: tolsesp = 1.0e-6_wp

contains

subroutine add_surface_hessian(self, mol, weight, hessian)
   !***********************************
   !* Add sum_A weight_A d2S_A/dR dR. Each grid point p of atom A contributes
   !* w_p L_p with the product L_p = prod_k sigma(t_k) over the neighbours k
   !* at distance t_k from the grid point, which moves with atom A, hence
   !*   d2L_p = L_p [dlnL dlnL^T + sum_k d2 ln sigma(t_k)].
   !* The grid is fixed in the laboratory frame like for the area.
   !***********************************
   !> Instance of the surface integrator
   class(surface_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Weight of each atomic surface area
   real(wp), intent(in) :: weight(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, ip, ino, nno, nni, kk, ii, jj, jdx
   integer, allocatable :: grdi(:)
   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp) :: xyzp(3), sasap, wsa, vec(3), tj, ehat(3), dlns, uj, sigma, d2sigma
   real(wp) :: block(3, 3)
   real(wp), allocatable :: grds(:, :), dlnl(:)

   nat = mol%nat
   call new_adjacency_list(list, mol, trans, self%srcut, .true.)
   allocate(grds(3, maxval(list%nnl)), grdi(maxval(list%nnl)), dlnl(3*nat))

   do iat = 1, nat
      if (weight(iat) == 0.0_wp) cycle
      ino = list%inl(iat)
      nno = list%nnl(iat)
      ii = 3*(iat - 1)

      do ip = 1, size(self%ang_grid, 2)
         xyzp(:) = mol%xyz(:, iat) + self%vdwsa(iat)*self%ang_grid(:, ip)
         call compute_w_spg(nat, list%nlat(ino+1:ino+nno), self%trj2, self%vdwsa, &
            & mol%xyz, xyzp, self%ah0, self%ah1, self%ah3, sasap, grds, nni, grdi)
         if (sasap <= tolsesp) cycle

         wsa = weight(iat)*self%ang_weight(ip)*self%wrp(iat)*sasap

         !>-- gradient of ln L_p, grds holds sigma'/sigma times the unit vector
         dlnl(:) = 0.0_wp
         do kk = 1, nni
            jj = 3*(grdi(kk) - 1)
            dlnl(ii+1:ii+3) = dlnl(ii+1:ii+3) + grds(:, kk)
            dlnl(jj+1:jj+3) = dlnl(jj+1:jj+3) - grds(:, kk)
         end do
         do jdx = 1, 3*nat
            hessian(:, jdx) = hessian(:, jdx) + wsa*dlnl*dlnl(jdx)
         end do

         !>-- curvature of ln sigma along each distance, sigma = ah0 + ah1 u + ah3 u**3
         do kk = 1, nni
            jat = grdi(kk)
            vec(:) = xyzp - mol%xyz(:, jat)
            tj = norm2(vec)
            ehat(:) = vec / tj
            dlns = dot_product(grds(:, kk), ehat)
            uj = tj - self%vdwsa(jat)
            sigma = self%ah0 + (self%ah1 + self%ah3*uj*uj)*uj
            d2sigma = 6*self%ah3*uj
            call get_pair_block(wsa*dlns, wsa*(d2sigma/sigma - dlns*dlns), tj, vec, block)
            call add_pair_hessian(hessian, iat, jat, block)
         end do
      end do
   end do
end subroutine add_surface_hessian

end module tblite_hessian_solvation_surface
