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

!> @file tblite/hessian/solvation/kernel.f90
!> Provides the second derivatives of the generalized Born interaction kernels

!> Generalized Born kernels J_AB = keps Phi(r_AB, rho_AB) and their second
!> derivatives w.r.t. the distance and the combined Born radius rho_AB
module tblite_hessian_solvation_kernel
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   use tblite_solvation_alpb, only : born_kernel
   implicit none
   private

   public :: get_kernel_derivs, get_radii_combination, add_born_kernel_hessian

   !> Parameter of the P16 kernel, identical to the one of the energy
   real(wp), parameter :: zeta = 1.028_wp
   real(wp), parameter :: zeta16 = zeta / 16.0_wp

contains

pure subroutine get_kernel_derivs(kernel, r, rho, dr, drho, drr, drrho, drhorho)
   !***********************************
   !* First and second derivatives of the kernel Phi(r, rho) w.r.t. the
   !* distance r and the combined Born radius rho.
   !* P16:   Phi = 1/f,       f = r + rho u**16, u = rho / (rho + zeta r/16)
   !* Still: Phi = F**(-1/2), F = r**2 + rho exp(-r**2/(4 rho))
   !***********************************
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel
   !> Distance between the two atoms
   real(wp), intent(in) :: r
   !> Combined Born radius of the pair
   real(wp), intent(in) :: rho
   !> First derivatives w.r.t. r and rho
   real(wp), intent(out) :: dr, drho
   !> Second derivatives w.r.t. r and rho
   real(wp), intent(out) :: drr, drrho, drhorho

   real(wp) :: f, fr, frho, frr, frrho, frhorho, u, lam, t, d, e, f32, f52

   select case(kernel)
   case default
      dr = 0.0_wp
      drho = 0.0_wp
      drr = 0.0_wp
      drrho = 0.0_wp
      drhorho = 0.0_wp
   case(born_kernel%p16)
      u = rho / (rho + zeta16*r)
      f = r + rho*u**16
      fr = 1.0_wp - zeta*u**17
      frho = u**16 + zeta*r*u**17/rho
      !>-- the second derivatives of f form a rank-one matrix
      lam = (17.0_wp*zeta**2/16.0_wp) * u**18/rho
      t = r/rho
      frr = lam
      frrho = -lam*t
      frhorho = lam*t*t

      dr = -fr/f**2
      drho = -frho/f**2
      drr = -frr/f**2 + 2*fr*fr/f**3
      drrho = -frrho/f**2 + 2*fr*frho/f**3
      drhorho = -frhorho/f**2 + 2*frho*frho/f**3
   case(born_kernel%still)
      d = 0.25_wp*r*r/rho
      e = exp(-d)
      f = r*r + rho*e
      fr = 2*r*(1.0_wp - 0.25_wp*e)
      frho = e*(1.0_wp + d)
      frr = 2.0_wp - 0.5_wp*e + d*e
      frrho = -r*e*d/(2*rho)
      frhorho = e*d*d/rho

      f32 = f*sqrt(f)
      f52 = f32*f
      dr = -fr/(2*f32)
      drho = -frho/(2*f32)
      drr = -frr/(2*f32) + 3*fr*fr/(4*f52)
      drrho = -frrho/(2*f32) + 3*fr*frho/(4*f52)
      drhorho = -frhorho/(2*f32) + 3*frho*frho/(4*f52)
   end select
end subroutine get_kernel_derivs

pure subroutine get_radii_combination(kernel, rada, radb, rho, da, db, daa, dab, dbb)
   !***********************************
   !* Combined Born radius of a pair entering the kernel and its derivatives
   !* w.r.t. both Born radii, sqrt(a b) for P16 and a b for the Still kernel.
   !***********************************
   !> Generalized Born interaction kernel
   integer, intent(in) :: kernel
   !> Born radii of the two atoms
   real(wp), intent(in) :: rada, radb
   !> Combined Born radius
   real(wp), intent(out) :: rho
   !> First derivatives w.r.t. both Born radii
   real(wp), intent(out) :: da, db
   !> Second derivatives w.r.t. both Born radii
   real(wp), intent(out) :: daa, dab, dbb

   select case(kernel)
   case default
      rho = 0.0_wp
      da = 0.0_wp
      db = 0.0_wp
      daa = 0.0_wp
      dab = 0.0_wp
      dbb = 0.0_wp
   case(born_kernel%p16)
      rho = sqrt(rada*radb)
      da = rho/(2*rada)
      db = rho/(2*radb)
      daa = -rho/(4*rada*rada)
      dbb = -rho/(4*radb*radb)
      dab = rho/(4*rada*radb)
   case(born_kernel%still)
      rho = rada*radb
      da = radb
      db = rada
      daa = 0.0_wp
      dbb = 0.0_wp
      dab = 1.0_wp
   end select
end subroutine get_radii_combination

subroutine add_born_kernel_hessian(mol, kernel, keps, qat, rad, hessian)
   !***********************************
   !* Add the Hessian of 1/2 q^T J q at fixed Born radii, the pair terms
   !* keps q_A q_B Phi(r_AB, rho_AB). The self terms keps q_A**2 / (2 a_A)
   !* depend on the geometry only via the Born radii.
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
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, jat
   real(wp) :: vec(3), r1, wq, rho, da, db, daa, dab, dbb
   real(wp) :: dr, drho, drr, drrho, drhorho, block(3, 3)

   do iat = 1, mol%nat
      do jat = 1, iat - 1
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         wq = keps*qat(iat)*qat(jat)
         call get_radii_combination(kernel, rad(iat), rad(jat), rho, da, db, daa, dab, dbb)
         call get_kernel_derivs(kernel, r1, rho, dr, drho, drr, drrho, drhorho)
         call get_pair_block(wq*dr, wq*drr, r1, vec, block)
         call add_pair_hessian(hessian, iat, jat, block)
      end do
   end do
end subroutine add_born_kernel_hessian

end module tblite_hessian_solvation_kernel
