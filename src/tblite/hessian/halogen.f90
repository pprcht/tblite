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

!> @file tblite/hessian/halogen.f90
!> Provides the second derivative of the halogen bonding correction

!> Analytic Hessian of the classical halogen bonding correction
module tblite_hessian_halogen
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_classical_halogen, only : halogen_correction
   use tblite_hessian_pair, only : add_hessian_block
   implicit none
   private

   public :: get_halogen_hessian

   !> Exponents of the angular damping and of the Lennard-Jones like potential,
   !> identical to the ones used for the energy
   real(wp), parameter :: alp = 6.0_wp, lj2 = 6.0_wp

   real(wp), parameter :: unity(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])

contains

subroutine get_halogen_hessian(self, mol, hessian, error)
   !***********************************
   !* Add the Hessian of the halogen bonding correction. For a halogen X,
   !* an acceptor J and the nearest neighbor K of the halogen the energy is
   !* E = c A(t) L(|a|) with a = R_J - R_X, b = R_K - R_X,
   !* A = (1/2 - t/4)**alp, t = 2 a.b/(|a||b|) and
   !* L = (u**2 - d u)/(1 + u**2), u = (r0/|a|)**lj2.
   !* Only available for molecular systems.
   !***********************************
   !> Instance of the halogen bonding correction
   type(halogen_correction), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: xat, jat, kat, iat, xzp, jzp, i, j
   real(wp) :: cc, r0, dist, rk, a(3), b(3), ra, rb, ra2, rb2, ab
   real(wp) :: u, ll, dldu, d2ldu2, dudr, d2udr2, dl, d2l, s, aa, daa, d2aa
   real(wp) :: ta(3), tb(3), taa(3, 3), tab(3, 3), tbb(3, 3)
   real(wp) :: la(3), laa(3, 3), ga(3), gb(3), gaa(3, 3), gab(3, 3), gbb(3, 3)
   real(wp) :: haa(3, 3), hab(3, 3), hbb(3, 3), hxj(3, 3), hxk(3, 3)

   if (any(mol%periodic)) then
      call fatal_error(error, "Halogen bond Hessian is not available for periodic systems")
      return
   end if
   if (count(self%halogen) * count(self%acceptor) == 0) return

   do xat = 1, mol%nat
      xzp = mol%id(xat)
      if (.not.self%halogen(xzp)) cycle

      !>-- nearest neighbor of the halogen, shared by all its acceptors
      dist = huge(1.0_wp)
      kat = 0
      do iat = 1, mol%nat
         rk = norm2(mol%xyz(:, xat) - mol%xyz(:, iat))
         if (rk < dist .and. rk > 0.0_wp) then
            dist = rk
            kat = iat
         end if
      end do
      if (kat == 0) cycle

      do jat = 1, mol%nat
         jzp = mol%id(jat)
         if (.not.self%acceptor(jzp)) cycle
         a(:) = mol%xyz(:, jat) - mol%xyz(:, xat)
         if (norm2(a) > self%cutoff) cycle
         b(:) = mol%xyz(:, kat) - mol%xyz(:, xat)
         cc = self%bond_strength(xzp)
         r0 = self%rad(xzp) + self%rad(jzp)

         ra2 = sum(a*a)
         rb2 = sum(b*b)
         ra = sqrt(ra2)
         rb = sqrt(rb2)
         ab = sum(a*b)

         !>-- radial factor and its derivatives w.r.t. |a|
         u = (r0/ra)**lj2
         ll = (u*u - self%damping*u) / (1.0_wp + u*u)
         dldu = (2*u - self%damping + self%damping*u*u) / (1.0_wp + u*u)**2
         d2ldu2 = (2.0_wp - 6*u*u + 6*self%damping*u - 2*self%damping*u**3) &
            & / (1.0_wp + u*u)**3
         dudr = -lj2*u/ra
         d2udr2 = lj2*(lj2 + 1.0_wp)*u/ra2
         dl = dldu*dudr
         d2l = d2ldu2*dudr*dudr + dldu*d2udr2

         !>-- angular factor and its derivatives w.r.t. t
         s = 0.5_wp - 0.5_wp*ab/(ra*rb)
         aa = s**alp
         daa = -0.25_wp*alp*s**(alp - 1.0_wp)
         d2aa = 0.0625_wp*alp*(alp - 1.0_wp)*s**(alp - 2.0_wp)

         !>-- derivatives of t w.r.t. the vectors a and b
         ta(:) = 2.0_wp/(ra*rb) * (b - ab/ra2*a)
         tb(:) = 2.0_wp/(ra*rb) * (a - ab/rb2*b)
         do j = 1, 3
            do i = 1, 3
               taa(i, j) = 2.0_wp/(rb*ra**3) * (3*ab/ra2*a(i)*a(j) &
                  & - (a(i)*b(j) + b(i)*a(j)) - ab*unity(i, j))
               tbb(i, j) = 2.0_wp/(ra*rb**3) * (3*ab/rb2*b(i)*b(j) &
                  & - (a(i)*b(j) + b(i)*a(j)) - ab*unity(i, j))
               tab(i, j) = 2.0_wp/(ra*rb) * (unity(i, j) - b(i)*b(j)/rb2 &
                  & - a(i)*a(j)/ra2 + ab/(ra2*rb2)*a(i)*b(j))
               laa(i, j) = d2l*a(i)*a(j)/ra2 + dl/ra*(unity(i, j) - a(i)*a(j)/ra2)
            end do
         end do
         la(:) = dl*a/ra
         ga(:) = daa*ta
         gb(:) = daa*tb
         do j = 1, 3
            do i = 1, 3
               gaa(i, j) = d2aa*ta(i)*ta(j) + daa*taa(i, j)
               gab(i, j) = d2aa*ta(i)*tb(j) + daa*tab(i, j)
               gbb(i, j) = d2aa*tb(i)*tb(j) + daa*tbb(i, j)
            end do
         end do

         do j = 1, 3
            do i = 1, 3
               haa(i, j) = cc*(gaa(i, j)*ll + ga(i)*la(j) + la(i)*ga(j) + aa*laa(i, j))
               hab(i, j) = cc*(gab(i, j)*ll + la(i)*gb(j))
               hbb(i, j) = cc*gbb(i, j)*ll
            end do
         end do

         !>-- the halogen enters both vectors with negative sign
         hxj(:, :) = -(haa + transpose(hab))
         hxk(:, :) = -(hab + hbb)
         call add_hessian_block(hessian, xat, xat, haa + hab + transpose(hab) + hbb)
         call add_hessian_block(hessian, jat, jat, haa)
         call add_hessian_block(hessian, kat, kat, hbb)
         call add_hessian_block(hessian, xat, jat, hxj)
         call add_hessian_block(hessian, jat, xat, transpose(hxj))
         call add_hessian_block(hessian, xat, kat, hxk)
         call add_hessian_block(hessian, kat, xat, transpose(hxk))
         call add_hessian_block(hessian, jat, kat, hab)
         call add_hessian_block(hessian, kat, jat, transpose(hab))
      end do
   end do
end subroutine get_halogen_hessian

end module tblite_hessian_halogen
