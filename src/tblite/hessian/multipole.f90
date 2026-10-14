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

!> @file tblite/hessian/multipole.f90
!> Provides derivatives of the anisotropic electrostatics required for the Hessian

!> Second derivatives of the damped multipole electrostatic energy. With the
!> distance vector u between two atoms, the damped kernels k_n = f_n(r, R0)/r**n
!> and the multipole radius R0 depending on the coordination numbers, the
!> energy consists of the charge-dipole, the dipole-dipole and the
!> charge-quadrupole interaction.
module tblite_hessian_multipole
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_coulomb_multipole, only : damped_multipole
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_pair, only : add_pair_hessian, add_ncoord_hessian
   use tblite_integral_native_integrals, only : multipole_grad_cgto, &
      & multipole_hess_cgto, msao
   implicit none
   private

   public :: get_multipole_hessian, get_multipole_potential_derivative
   public :: get_multipole_band_hessian

   real(wp), parameter :: unity(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])

   !> Inverse power of the distance and index of the damping function for the
   !> charge-dipole, dipole-dipole and charge-quadrupole interaction
   integer, parameter :: rpow(3) = [3, 5, 5], idmp(3) = [1, 2, 2]

contains

pure subroutine get_damping(kexp, r, rad, f, dfdr, dfdm, d2fdr2, d2fdrdm, d2fdm2)
   !***********************************
   !* Damping function f = 1/(1 + 6 (rad/r)**kexp) with its first and second
   !* derivatives w.r.t. the distance r and the multipole radius rad.
   !***********************************
   !> Exponent of the damping function
   real(wp), intent(in) :: kexp
   !> Distance between the atoms
   real(wp), intent(in) :: r
   !> Multipole radius of the pair
   real(wp), intent(in) :: rad
   !> Damping function and its derivatives
   real(wp), intent(out) :: f, dfdr, dfdm, d2fdr2, d2fdrdm, d2fdm2

   real(wp) :: t, den, dfdt, d2fdt2, dpdt, d2pdt2

   t = rad/r
   den = 1.0_wp + 6.0_wp*t**kexp
   dpdt = 6.0_wp*kexp*t**(kexp - 1.0_wp)
   d2pdt2 = 6.0_wp*kexp*(kexp - 1.0_wp)*t**(kexp - 2.0_wp)
   f = 1.0_wp/den
   dfdt = -dpdt/den**2
   d2fdt2 = 2*dpdt**2/den**3 - d2pdt2/den**2
   dfdr = -dfdt*t/r
   dfdm = dfdt/r
   d2fdr2 = d2fdt2*(t/r)**2 + 2*dfdt*t/r**2
   d2fdrdm = -d2fdt2*t/r**2 - dfdt/r**2
   d2fdm2 = d2fdt2/r**2
end subroutine get_damping

pure subroutine get_multipole_radii(self, mol, cn, rad, drad, d2rad)
   !***********************************
   !* Multipole radius of each atom with its first and second derivative
   !* w.r.t. the coordination number.
   !***********************************
   !> Instance of the multipole container
   type(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Coordination numbers
   real(wp), intent(in) :: cn(:)
   !> Multipole radii and their derivatives
   real(wp), intent(out) :: rad(:), drad(:), d2rad(:)

   integer :: iat, izp
   real(wp) :: expterm, delta

   do iat = 1, mol%nat
      izp = mol%id(iat)
      expterm = exp(-self%kexp*(cn(iat) - self%valence_cn(izp) - self%shift))
      delta = self%rmax - self%rad(izp)
      rad(iat) = self%rad(izp) + delta/(1.0_wp + expterm)
      drad(iat) = delta*self%kexp*expterm/(1.0_wp + expterm)**2
      d2rad(iat) = delta*self%kexp**2*expterm*(expterm - 1.0_wp)/(1.0_wp + expterm)**3
   end do
end subroutine get_multipole_radii

pure function unpack(qp) result(theta)
   !> Quadrupole moment packed as xx, xy, yy, xz, yz, zz
   real(wp), intent(in) :: qp(6)
   !> Quadrupole moment tensor
   real(wp) :: theta(3, 3)

   theta(:, :) = reshape([qp(1), qp(2), qp(4), qp(2), qp(3), qp(5), &
      & qp(4), qp(5), qp(6)], [3, 3])
end function unpack

pure function pack(theta) result(qp)
   !> Symmetric tensor
   real(wp), intent(in) :: theta(3, 3)
   !> Packed tensor as xx, xy, yy, xz, yz, zz with summed off-diagonal elements
   real(wp) :: qp(6)

   qp(:) = [theta(1, 1), theta(1, 2) + theta(2, 1), theta(2, 2), &
      & theta(1, 3) + theta(3, 1), theta(2, 3) + theta(3, 2), theta(3, 3)]
end function pack

pure function outer(a, b) result(ab)
   real(wp), intent(in) :: a(3), b(3)
   real(wp) :: ab(3, 3)

   ab(:, :) = spread(a, 2, 3) * spread(b, 1, 3)
end function outer

subroutine get_multipole_hessian(self, mol, qat, dpat, qpat, hessian)
   !***********************************
   !* Add the Hessian of the multipole electrostatic energy at fixed atomic
   !* charges, dipole and quadrupole moments. For a pair with u = R_j - R_i
   !* the energy is the sum of k_n(r, R0) P_n(u) over the three interactions
   !* with the polynomials P_n in u. Besides the direct second derivative the
   !* dependence of the multipole radius R0 on the coordination numbers is
   !* accounted for. Only available for molecular systems.
   !***********************************
   !> Instance of the multipole container
   type(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, n, ii, jj
   real(wp) :: vec(3), r, r2, eh(3), mrad, kdmp(2)
   real(wp) :: f(2), dfdr(2), dfdm(2), d2fdr2(2), d2fdrdm(2), d2fdm2(2)
   real(wp) :: pol(3), dpol(3, 3), d2pol(3, 3, 3), ti(3, 3), tj(3, 3), tiv(3), tjv(3)
   real(wp) :: g, dg, d2g, k, dk, d2k, dkdm, d2kdrdm, block(3, 3)
   real(wp) :: dedm, d2edm2, d2edrdm(3)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), lattr(:, :)
   real(wp), allocatable :: rad(:), drad(:), d2rad(:)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :), dEdcndcn(:, :)

   nat = mol%nat
   kdmp(:) = [self%kdmp3, self%kdmp5]
   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat), rad(nat), drad(nat), d2rad(nat))
   call self%ncoord%get_cn(mol, cn, dcndr, dcndL)
   call get_multipole_radii(self, mol, cn, rad, drad, d2rad)
   allocate(dEdcn(nat), dEdcndr(3*nat, nat), dEdcndcn(nat, nat), source=0.0_wp)

   do iat = 1, nat
      ii = 3*(iat - 1)
      do jat = 1, iat - 1
         jj = 3*(jat - 1)
         vec(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = sum(vec**2)
         r = sqrt(r2)
         eh(:) = vec/r
         mrad = 0.5_wp*(rad(iat) + rad(jat))
         do n = 1, 2
            call get_damping(kdmp(n), r, mrad, f(n), dfdr(n), dfdm(n), d2fdr2(n), &
               & d2fdrdm(n), d2fdm2(n))
         end do

         ti(:, :) = unpack(qpat(:, iat))
         tj(:, :) = unpack(qpat(:, jat))
         tiv(:) = matmul(ti, vec)
         tjv(:) = matmul(tj, vec)

         !>-- charge-dipole
         pol(1) = qat(jat)*dot_product(dpat(:, iat), vec) &
            & - qat(iat)*dot_product(dpat(:, jat), vec)
         dpol(:, 1) = qat(jat)*dpat(:, iat) - qat(iat)*dpat(:, jat)
         d2pol(:, :, 1) = 0.0_wp
         !>-- dipole-dipole
         pol(2) = dot_product(dpat(:, iat), dpat(:, jat))*r2 &
            & - 3*dot_product(dpat(:, jat), vec)*dot_product(dpat(:, iat), vec)
         dpol(:, 2) = 2*dot_product(dpat(:, iat), dpat(:, jat))*vec &
            & - 3*(dot_product(dpat(:, jat), vec)*dpat(:, iat) &
            & + dot_product(dpat(:, iat), vec)*dpat(:, jat))
         d2pol(:, :, 2) = 2*dot_product(dpat(:, iat), dpat(:, jat))*unity &
            & - 3*(outer(dpat(:, iat), dpat(:, jat)) + outer(dpat(:, jat), dpat(:, iat)))
         !>-- charge-quadrupole
         pol(3) = qat(jat)*dot_product(vec, tiv) + qat(iat)*dot_product(vec, tjv)
         dpol(:, 3) = 2*(qat(jat)*tiv + qat(iat)*tjv)
         d2pol(:, :, 3) = 2*(qat(jat)*ti + qat(iat)*tj)

         block(:, :) = 0.0_wp
         dedm = 0.0_wp
         d2edm2 = 0.0_wp
         d2edrdm(:) = 0.0_wp
         do n = 1, 3
            g = r**(-rpow(n))
            dg = -rpow(n)*g/r
            d2g = rpow(n)*(rpow(n) + 1)*g/r2
            k = f(idmp(n))*g
            dk = dfdr(idmp(n))*g + f(idmp(n))*dg
            d2k = d2fdr2(idmp(n))*g + 2*dfdr(idmp(n))*dg + f(idmp(n))*d2g
            dkdm = dfdm(idmp(n))*g
            d2kdrdm = d2fdrdm(idmp(n))*g + dfdm(idmp(n))*dg

            block(:, :) = block + k*d2pol(:, :, n) &
               & + dk*(outer(eh, dpol(:, n)) + outer(dpol(:, n), eh)) &
               & + pol(n)*(d2k*outer(eh, eh) + dk/r*(unity - outer(eh, eh)))

            dedm = dedm + pol(n)*dkdm
            d2edm2 = d2edm2 + pol(n)*d2fdm2(idmp(n))*g
            d2edrdm(:) = d2edrdm + dkdm*dpol(:, n) + pol(n)*d2kdrdm*eh
         end do
         call add_pair_hessian(hessian, iat, jat, block)

         !>-- the multipole radius of the pair is the average of the atomic ones
         dEdcn(iat) = dEdcn(iat) + 0.5_wp*dedm*drad(iat)
         dEdcn(jat) = dEdcn(jat) + 0.5_wp*dedm*drad(jat)
         dEdcndr(ii+1:ii+3, iat) = dEdcndr(ii+1:ii+3, iat) - 0.5_wp*drad(iat)*d2edrdm
         dEdcndr(jj+1:jj+3, iat) = dEdcndr(jj+1:jj+3, iat) + 0.5_wp*drad(iat)*d2edrdm
         dEdcndr(ii+1:ii+3, jat) = dEdcndr(ii+1:ii+3, jat) - 0.5_wp*drad(jat)*d2edrdm
         dEdcndr(jj+1:jj+3, jat) = dEdcndr(jj+1:jj+3, jat) + 0.5_wp*drad(jat)*d2edrdm
         dEdcndcn(iat, iat) = dEdcndcn(iat, iat) + 0.25_wp*d2edm2*drad(iat)**2 &
            & + 0.5_wp*dedm*d2rad(iat)
         dEdcndcn(jat, jat) = dEdcndcn(jat, jat) + 0.25_wp*d2edm2*drad(jat)**2 &
            & + 0.5_wp*dedm*d2rad(jat)
         dEdcndcn(iat, jat) = dEdcndcn(iat, jat) + 0.25_wp*d2edm2*drad(iat)*drad(jat)
         dEdcndcn(jat, iat) = dEdcndcn(jat, iat) + 0.25_wp*d2edm2*drad(iat)*drad(jat)
      end do
   end do

   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, lattr)
   call self%ncoord%add_coordination_number_hessian(mol, lattr, dEdcn, hessian)
   call add_ncoord_hessian(hessian, dcndr, dEdcndr, dEdcndcn)
end subroutine get_multipole_hessian

subroutine get_multipole_potential_derivative(self, mol, qat, dpat, qpat, dvat, dvdp, dvqp)
   !***********************************
   !* Derivative of the multipole electrostatic potentials w.r.t. the atomic
   !* positions at fixed charges and multipole moments. The derivatives are
   !* indexed as (potential component, atom of the potential, direction,
   !* displaced atom), the atomic potential has no component.
   !* The quadrupole potential acts on the packed quadrupole moments, hence
   !* its off-diagonal components are summed.
   !***********************************
   !> Instance of the multipole container
   type(damped_multipole), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Atomic dipole moments
   real(wp), intent(in) :: dpat(:, :)
   !> Atomic quadrupole moments
   real(wp), intent(in) :: qpat(:, :)
   !> Derivative of the potential acting on the charges
   real(wp), intent(out) :: dvat(:, :, :)
   !> Derivative of the potential acting on the dipole moments
   real(wp), intent(out) :: dvdp(:, :, :, :)
   !> Derivative of the potential acting on the quadrupole moments
   real(wp), intent(out) :: dvqp(:, :, :, :)

   integer :: nat, iat, jat, kat, ic
   real(wp) :: vec(3), r, r2, eh(3), mrad, f3, f5, df3dr, df5dr, df3dm, df5dm, dum(3)
   real(wp) :: k3, k5, dk3, dk5, dk3dm, dk5dm, sd, sq, tj(3, 3), tjv(3), dj(3)
   real(wp) :: dv(3), dvd(3, 3), dvq(6, 3), dvdm, dvddm(3), dvqdm(6), dpm(3, 3)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: rad(:), drad(:), d2rad(:), dmrad(:, :)

   nat = mol%nat
   dvat(:, :, :) = 0.0_wp
   dvdp(:, :, :, :) = 0.0_wp
   dvqp(:, :, :, :) = 0.0_wp

   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat), rad(nat), drad(nat), &
      & d2rad(nat), dmrad(3, nat))
   call self%ncoord%get_cn(mol, cn, dcndr, dcndL)
   call get_multipole_radii(self, mol, cn, rad, drad, d2rad)

   do iat = 1, nat
      do jat = 1, nat
         if (iat == jat) cycle
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         r = sqrt(r2)
         eh(:) = vec/r
         mrad = 0.5_wp*(rad(iat) + rad(jat))
         call get_damping(self%kdmp3, r, mrad, f3, df3dr, df3dm, dum(1), dum(2), dum(3))
         call get_damping(self%kdmp5, r, mrad, f5, df5dr, df5dm, dum(1), dum(2), dum(3))
         k3 = f3/r**3
         dk3 = df3dr/r**3 - 3*k3/r
         dk3dm = df3dm/r**3
         k5 = f5/r**5
         dk5 = df5dr/r**5 - 5*k5/r
         dk5dm = df5dm/r**5

         dj(:) = dpat(:, jat)
         tj(:, :) = unpack(qpat(:, jat))
         tjv(:) = matmul(tj, vec)
         sd = dot_product(vec, dj)
         sq = dot_product(vec, tjv)

         !>-- potential on the charge, k3 (u.d_j) + k5 (u.T_j.u)
         dv(:) = (dk3*sd + dk5*sq)*eh + k3*dj + 2*k5*tjv
         dvdm = dk3dm*sd + dk5dm*sq

         !>-- potential on the dipole, -q_j k3 u + k5 r**2 d_j - 3 k5 (u.d_j) u,
         !>   the second term is written with the damped kernel of third order
         dvddm(:) = -qat(jat)*dk3dm*vec + df5dm/r**3*dj - 3*dk5dm*sd*vec
         do ic = 1, 3
            dvd(:, ic) = -qat(jat)*(dk3*eh(ic)*vec + k3*unity(:, ic)) &
               & + (df5dr/r**3 - 3*f5/r**4)*eh(ic)*dj &
               & - 3*(dk5*eh(ic)*sd*vec + k5*(dj(ic)*vec + sd*unity(:, ic)))
         end do

         !>-- potential on the quadrupole, q_j k5 u u^T
         dvqdm(:) = qat(jat)*dk5dm*pack(outer(vec, vec))
         do ic = 1, 3
            dpm(:, :) = outer(unity(:, ic), vec) + outer(vec, unity(:, ic))
            dvq(:, ic) = qat(jat)*pack(dk5*eh(ic)*outer(vec, vec) + k5*dpm)
         end do

         dvat(iat, :, iat) = dvat(iat, :, iat) + dv
         dvat(iat, :, jat) = dvat(iat, :, jat) - dv
         dvdp(:, iat, :, iat) = dvdp(:, iat, :, iat) + dvd
         dvdp(:, iat, :, jat) = dvdp(:, iat, :, jat) - dvd
         dvqp(:, iat, :, iat) = dvqp(:, iat, :, iat) + dvq
         dvqp(:, iat, :, jat) = dvqp(:, iat, :, jat) - dvq

         !>-- dependence of the multipole radius on the coordination numbers
         dmrad(:, :) = 0.5_wp*(drad(iat)*dcndr(:, :, iat) + drad(jat)*dcndr(:, :, jat))
         dvat(iat, :, :) = dvat(iat, :, :) + dvdm*dmrad
         do kat = 1, nat
            do ic = 1, 3
               dvdp(:, iat, ic, kat) = dvdp(:, iat, ic, kat) + dvddm*dmrad(ic, kat)
               dvqp(:, iat, ic, kat) = dvqp(:, iat, ic, kat) + dvqdm*dmrad(ic, kat)
            end do
         end do
      end do
   end do
end subroutine get_multipole_potential_derivative

subroutine get_multipole_band_hessian(mol, bas, pmat, vdp, vqp, dvdp, dvqp, trans, &
      & cutoff, hessian)
   !***********************************
   !* Add the second derivative of the interaction of the density with the
   !* multipole potentials at fixed density and moments. The potentials act
   !* via the dipole and quadrupole integrals, hence the second derivative of
   !* the integrals with the potential and the first derivative of the
   !* integrals with the derivative of the potential contribute.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :)
   !> Potential acting on the dipole moments
   real(wp), intent(in) :: vdp(:, :)
   !> Potential acting on the quadrupole moments
   real(wp), intent(in) :: vqp(:, :)
   !> Derivative of the potential acting on the dipole moments
   real(wp), intent(in) :: dvdp(:, :, :, :)
   !> Derivative of the potential acting on the quadrupole moments
   real(wp), intent(in) :: dvqp(:, :, :, :)
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, kat, izp, jzp, itr, is, js, ish, jsh, ii, jj
   integer :: iao, jao, nao, ij, ic, jc, kc, ia, ja, ka
   real(wp) :: r2, vec(3), cutoff2, pij, block(3, 3), gdi(3, 3), gdj(3, 3)
   real(wp) :: gqi(3, 6), gqj(3, 6), val(3)
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), dstmp(:, :)
   real(wp), allocatable :: ddj(:, :, :), dqj(:, :, :), ddi(:, :, :), dqi(:, :, :)
   real(wp), allocatable :: d2s(:, :, :), d2dj(:, :, :, :), d2qj(:, :, :, :)
   real(wp), allocatable :: d2di(:, :, :, :), d2qi(:, :, :, :), cross(:, :)

   nat = mol%nat
   cutoff2 = cutoff**2
   nao = msao(bas%maxl)**2
   allocate(stmp(nao), dtmp(3, nao), qtmp(6, nao), dstmp(3, nao), ddj(3, 3, nao), &
      & dqj(3, 6, nao), ddi(3, 3, nao), dqi(3, 6, nao), d2s(3, 3, nao), &
      & d2dj(3, 3, 3, nao), d2qj(3, 3, 6, nao), d2di(3, 3, 3, nao), d2qi(3, 3, 6, nao))
   allocate(cross(3*nat, 3*nat), source=0.0_wp)

   do iat = 1, nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      ia = 3*(iat - 1)
      do jat = 1, nat
         if (iat == jat) cycle
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         ja = 3*(jat - 1)
         block(:, :) = 0.0_wp
         gdi(:, :) = 0.0_wp
         gdj(:, :) = 0.0_wp
         gqi(:, :) = 0.0_wp
         gqj(:, :) = 0.0_wp
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_hess_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddj, dqj, ddi, dqi, &
                     & d2s, d2dj, d2qj, d2di, d2qi)
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        pij = pmat(jj+jao, ii+iao)
                        do jc = 1, 3
                           do ic = 1, 3
                              block(ic, jc) = block(ic, jc) - pij*( &
                                 & dot_product(d2di(ic, jc, :, ij), vdp(:, iat)) &
                                 & + dot_product(d2dj(ic, jc, :, ij), vdp(:, jat)) &
                                 & + dot_product(d2qi(ic, jc, :, ij), vqp(:, iat)) &
                                 & + dot_product(d2qj(ic, jc, :, ij), vqp(:, jat)))
                           end do
                        end do
                        gdi(:, :) = gdi + pij*ddi(:, :, ij)
                        gdj(:, :) = gdj + pij*ddj(:, :, ij)
                        gqi(:, :) = gqi + pij*dqi(:, :, ij)
                        gqj(:, :) = gqj + pij*dqj(:, :, ij)
                     end do
                  end do
               end do
            end do
         end do

         !>-- both orders of each pair are visited
         call add_pair_hessian(hessian, iat, jat, 0.5_wp*block)

         do kat = 1, nat
            ka = 3*(kat - 1)
            do kc = 1, 3
               val(:) = matmul(gdi, dvdp(:, iat, kc, kat)) &
                  & + matmul(gdj, dvdp(:, jat, kc, kat)) &
                  & + matmul(gqi, dvqp(:, iat, kc, kat)) &
                  & + matmul(gqj, dvqp(:, jat, kc, kat))
               cross(ia+1:ia+3, ka+kc) = cross(ia+1:ia+3, ka+kc) - val
               cross(ja+1:ja+3, ka+kc) = cross(ja+1:ja+3, ka+kc) + val
            end do
         end do
      end do
   end do

   hessian(:, :) = hessian + 0.25_wp*(cross + transpose(cross))
end subroutine get_multipole_band_hessian

end module tblite_hessian_multipole
