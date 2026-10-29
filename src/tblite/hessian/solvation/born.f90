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

!> @file tblite/hessian/solvation/born.f90
!> Provides the second derivatives of the Born radii

!> Second derivatives of the GBOBC Born radii a_A = A(psi_A), with the
!> descreening integral psi_A a sum of pair contributions
module tblite_hessian_solvation_born
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_adjlist, only : adjacency_list, new_adjacency_list
   use tblite_blas, only : gemm
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_solvation_born, only : born_integrator, compute_psi
   implicit none
   private

   public :: get_descreening_derivs, get_obc_derivs
   public :: add_psi_hessian, add_born_radii_hessian

contains

pure subroutine get_descreening_derivs(r, rho, rvdw, g, dg, d2g)
   !***********************************
   !* Contribution g(r) of a descreening sphere with radius rho at distance r
   !* to the integral psi of an atom with van der Waals radius rvdw, and its
   !* first and second derivative w.r.t. r. The branches follow compute_psi:
   !* separated spheres, overlapping spheres, and a buried descreening sphere
   !* that does not contribute.
   !***********************************
   !> Distance between the two atoms
   real(wp), intent(in) :: r
   !> Descreening radius of the other atom
   real(wp), intent(in) :: rho
   !> Van der Waals radius of the atom owning the integral
   real(wp), intent(in) :: rvdw
   !> Contribution and its derivatives w.r.t. the distance
   real(wp), intent(out) :: g, dg, d2g

   real(wp) :: ap, am, ab, rhab, lnab, wov, dwov, d2wov

   if (r >= rvdw + rho) then
      ap = r + rho
      am = r - rho
      ab = ap*am
      rhab = rho/ab
      lnab = 0.5_wp*log(am/ap)/r
      g = rhab + lnab
      dg = -2*r*rho/ab**2 + (rhab - lnab)/r
      d2g = -4*rhab/ab + 8*r*r*rhab/ab**2 + 2*(lnab - rhab)/r**2
   else if (r + rho > rvdw) then
      ap = r + rho
      am = r - rho
      wov = 0.5_wp*am*(1.0_wp/ap - ap/rvdw**2) - log(ap/rvdw)
      dwov = 0.5_wp*((1.0_wp/ap - ap/rvdw**2) + am*(-1.0_wp/ap**2 - 1.0_wp/rvdw**2)) &
         & - 1.0_wp/ap
      d2wov = am/ap**3 - 1.0_wp/rvdw**2
      g = 1.0_wp/rvdw - 1.0_wp/ap + wov/(2*r)
      dg = 1.0_wp/ap**2 + dwov/(2*r) - wov/(2*r**2)
      d2g = -2.0_wp/ap**3 + d2wov/(2*r) - dwov/r**2 + wov/r**3
   else
      g = 0.0_wp
      dg = 0.0_wp
      d2g = 0.0_wp
   end if
end subroutine get_descreening_derivs

pure subroutine get_obc_derivs(psi, svdw, vdwr, born_scale, obc, da, d2a)
   !***********************************
   !* First and second derivative of the Born radius
   !* a = s / (1/svdw - tanh(P(x))/vdwr) w.r.t. the descreening integral,
   !* with x = svdw psi / 2 and the volume polynomial P(x) = x(c1 - c2 x + c3 x**2).
   !***********************************
   !> Descreening integral
   real(wp), intent(in) :: psi
   !> Van der Waals radius with offset
   real(wp), intent(in) :: svdw
   !> Van der Waals radius
   real(wp), intent(in) :: vdwr
   !> Scaling factor of the Born radii
   real(wp), intent(in) :: born_scale
   !> Volume polynomial correction
   real(wp), intent(in) :: obc(3)
   !> First and second derivative of the Born radius w.r.t. psi
   real(wp), intent(out) :: da, d2a

   real(wp) :: lam, x, p, dp, d2p, th, sech2, denom

   lam = 0.5_wp*svdw
   x = lam*psi
   p = x*(obc(1) + x*(obc(3)*x - obc(2)))
   dp = obc(1) - 2*obc(2)*x + 3*obc(3)*x*x
   d2p = -2*obc(2) + 6*obc(3)*x
   th = tanh(p)
   sech2 = 1.0_wp - th*th
   denom = 1.0_wp/svdw - th/vdwr
   da = born_scale*(lam/vdwr)*sech2*dp/denom**2
   d2a = born_scale*(lam*lam/vdwr)*sech2*((d2p - 2*th*dp*dp)/denom**2 &
      & + 2*sech2*dp*dp/(vdwr*denom**3))
end subroutine get_obc_derivs

subroutine add_psi_hessian(mol, list, vdwr, rho, weight, hessian)
   !***********************************
   !* Add sum_A weight_A d2psi_A/dR dR. Both directional contributions of a
   !* pair share the distance and combine into one radial block. The pairs
   !* are taken from the same neighbour list as the descreening integral.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Neighbour list of the descreening integral
   type(adjacency_list), intent(in) :: list
   !> Van der Waals radii
   real(wp), intent(in) :: vdwr(:)
   !> Descreening radii
   real(wp), intent(in) :: rho(:)
   !> Weight of each descreening integral
   real(wp), intent(in) :: weight(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat, img, inl
   real(wp) :: vec(3), r1, g, dgij, d2gij, dgji, d2gji, block(3, 3)

   do iat = 1, mol%nat
      inl = list%inl(iat)
      do img = 1, list%nnl(iat)
         jat = list%nlat(inl+img)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         call get_descreening_derivs(r1, rho(jat), vdwr(iat), g, dgij, d2gij)
         call get_descreening_derivs(r1, rho(iat), vdwr(jat), g, dgji, d2gji)
         call get_pair_block(weight(iat)*dgij + weight(jat)*dgji, &
            & weight(iat)*d2gij + weight(jat)*d2gji, r1, vec, block)
         call add_pair_hessian(hessian, iat, jat, block)
      end do
   end do
end subroutine add_psi_hessian

subroutine add_born_radii_hessian(self, mol, weight, hessian)
   !***********************************
   !* Add sum_A weight_A d2a_A/dR dR
   !*   = sum_A weight_A [A'_A d2psi_A/dR dR + A''_A dpsi_A/dR dpsi_A/dR].
   !***********************************
   !> Instance of the Born integrator
   class(born_integrator), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Weight of each Born radius
   real(wp), intent(in) :: weight(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, nat
   type(adjacency_list) :: list
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: psi(:), dpsidr(:, :, :), wda(:), gpsi(:, :), wgpsi(:, :)
   real(wp) :: da, d2a

   nat = mol%nat
   call new_adjacency_list(list, mol, trans, self%lrcut)
   allocate(psi(nat), dpsidr(3, nat, nat), wda(nat))
   call compute_psi(nat, mol%xyz, list, self%vdwr, self%rho, psi, dpsidr)

   allocate(gpsi(3*nat, nat), wgpsi(3*nat, nat))
   gpsi(:, :) = reshape(dpsidr, [3*nat, nat])
   do iat = 1, nat
      call get_obc_derivs(psi(iat), self%svdw(iat), self%vdwr(iat), self%born_scale, &
         & self%obc, da, d2a)
      wda(iat) = weight(iat)*da
      wgpsi(:, iat) = weight(iat)*d2a*gpsi(:, iat)
   end do

   call add_psi_hessian(mol, list, self%vdwr, self%rho, wda, hessian)
   call gemm(wgpsi, gpsi, hessian, transb="t", beta=1.0_wp)
end subroutine add_born_radii_hessian

end module tblite_hessian_solvation_born
