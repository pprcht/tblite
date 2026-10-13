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

!> @file tblite/hessian/repulsion.f90
!> Provides the second derivative of the repulsion energy

!> Analytic Hessian of the screened Coulomb repulsion
module tblite_hessian_repulsion
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_repulsion_effective, only : tb_repulsion
   implicit none
   private

   public :: get_repulsion_hessian

contains

subroutine get_repulsion_hessian(self, mol, hessian)
   !***********************************
   !* Add the Hessian of the repulsion energy, a sum of pair terms
   !* f(r) = zeff exp(-alpha r**kexp) / r**rexp. With g = f'/f the
   !* derivatives are f' = g f and f'' = (g' + g**2) f.
   !* Interactions of an atom with its own periodic images are not included.
   !***********************************
   !> Instance of the repulsion container
   class(tb_repulsion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r1, r2, vec(3), cutoff2, r1k, fval, g, dg, block(3, 3)
   real(wp), allocatable :: trans(:, :)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff, trans)
   cutoff2 = self%cutoff**2

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = sum(vec**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)
            r1k = r1**self%kexp(jzp, izp)
            fval = self%zeff(jzp, izp) * exp(-self%alpha(jzp, izp)*r1k) &
               & / r1**self%rexp(jzp, izp)
            g = -(self%alpha(jzp, izp)*self%kexp(jzp, izp)*r1k + self%rexp(jzp, izp)) / r1
            dg = (self%rexp(jzp, izp) - self%alpha(jzp, izp)*self%kexp(jzp, izp) &
               & * (self%kexp(jzp, izp) - 1.0_wp)*r1k) / r2
            call get_pair_block(g*fval, (dg + g*g)*fval, r1, vec, block)
            call add_pair_hessian(hessian, iat, jat, block)
         end do
      end do
   end do
end subroutine get_repulsion_hessian

end module tblite_hessian_repulsion
