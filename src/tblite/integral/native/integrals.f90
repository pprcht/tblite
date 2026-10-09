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

!> @file tblite/integral/native/integrals.f90
!> Provides native evaluation of Gaussian integrals

!> Implementation of native overlap, dipole, and quadrupole integrals
module tblite_integral_native_integrals
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_integral_diat_trafo, only: diat_trafo_cache, setup_diat_trafo, &
      & diat_trafo
   use tblite_integral_trafo, only : transform0, transform1, transform2
   implicit none
   private

   public :: overlap_cgto, overlap_grad_cgto, get_overlap
   public :: dipole_cgto, dipole_grad_cgto, get_dipole_integrals
   public :: multipole_cgto, multipole_grad_cgto, get_multipole_integrals
   public :: overlap_hess_cgto, multipole_hess_cgto
   public :: maxl, msao, smap, sdim

   interface get_overlap
      module procedure :: get_overlap_lat
      module procedure :: get_overlap_diat_lat
   end interface get_overlap

   interface get_dipole_integrals
      module procedure :: get_dipole_integrals_lat
      module procedure :: get_dipole_integrals_diat_lat
   end interface get_dipole_integrals

   interface get_multipole_integrals
      module procedure :: get_multipole_integrals_lat
      module procedure :: get_multipole_integrals_diat_lat
   end interface get_multipole_integrals

   integer, parameter :: maxl = 6
   integer, parameter :: maxl2 = maxl*2
   integer, parameter :: msao(0:maxl) = [1, 3, 5, 7, 9, 11, 13]
   integer, parameter :: mlao(0:maxl) = [1, 3, 6, 10, 15, 21, 28]
   integer, parameter :: smap(0:maxl) = [0, 1, 4, 9, 16, 25, 36]
   integer, parameter :: lmap(0:maxl) = [0, 1, 4, 10, 20, 35, 56]
   integer, parameter :: sdim(0:maxl) = [1, 4, 9, 16, 25, 36, 49]
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrtpi3 = sqrtpi**3

   !> Cartesian exponents in CCA ordering. For angular momentum l, components
   !> are ordered by decreasing lx and, for equal lx, by decreasing ly, with
   !> lz = l - lx - ly.
   integer, parameter :: lx(3, 84) = reshape([&
      & 0, &
      & 1,0,0, &
      & 2,1,1,0,0,0, &
      & 3,2,2,1,1,1,0,0,0,0, &
      & 4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, &
      & 5,4,4,3,3,3,2,2,2,2,1,1,1,1,1,0,0,0,0,0,0, &
      & 6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1,0,0,0,0,0,0,0, &
      & 0, &
      & 0,1,0, &
      & 0,1,0,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0, &
      & 0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,5,4,3,2,1,0,6,5,4,3,2,1,0, &
      & 0, &
      & 0,0,1, &
      & 0,0,1,0,1,2, &
      & 0,0,1,0,1,2,0,1,2,3, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5, &
      & 0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,0,1,2,3,4,5,0,1,2,3,4,5,6], &
      & shape(lx), order=[2, 1])

contains

elemental function overlap_1d(moment, alpha) result(overlap)
   integer, intent(in) :: moment
   real(wp), intent(in) :: alpha
   real(wp) :: overlap
   real(wp), parameter :: dfactorial(0:8) = & ! see OEIS A001147
      & [1.0_wp,1.0_wp,3.0_wp,15.0_wp,105.0_wp,945.0_wp,10395.0_wp,135135.0_wp, &
      & 2027025.0_wp]

   if (modulo(moment, 2) == 0) then
      overlap = (0.5_wp/alpha)**(moment/2) * dfactorial(moment/2)
   else
      overlap = 0.0_wp
   end if
end function overlap_1d

pure subroutine horizontal_shift(ae, l, cfs)
   integer, intent(in) :: l
   real(wp), intent(in) :: ae
   ! allow(C071): assumed-size retained from the reference integral kernels
   real(wp), intent(inout) :: cfs(*)
   select case(l)
   case(0) ! s
      continue
   case(1) ! p
      cfs(1)=cfs(1)+ae*cfs(2)
   case(2) ! d
      cfs(1)=cfs(1)+ae*ae*cfs(3)
      cfs(2)=cfs(2)+ 2*ae*cfs(3)
   case(3) ! f
      cfs(1)=cfs(1)+ae*ae*ae*cfs(4)
      cfs(2)=cfs(2)+ 3*ae*ae*cfs(4)
      cfs(3)=cfs(3)+ 3*ae*cfs(4)
   case(4) ! g
      cfs(1)=cfs(1)+ae*ae*ae*ae*cfs(5)
      cfs(2)=cfs(2)+ 4*ae*ae*ae*cfs(5)
      cfs(3)=cfs(3)+ 6*ae*ae*cfs(5)
      cfs(4)=cfs(4)+ 4*ae*cfs(5)
   case default
      continue
   end select
end subroutine horizontal_shift

pure subroutine form_product(a, b, la, lb, d)
   integer, intent(in) :: la, lb
   ! allow(C071): assumed-size retained from the reference integral kernels
   real(wp), intent(in) :: a(*), b(*)
   ! allow(C071): assumed-size retained from the reference integral kernels
   real(wp), intent(inout) :: d(*)
   if(la>=4.or.lb>=4) goto 40
   if(la>=3.or.lb>=3) goto 30
   if(la>=2.or.lb>=2) goto 20
   ! <s|s> = <s>
   d(1)=a(1)*b(1)
   if(la==0.and.lb==0) return
   ! <s|p> = <s|*(|s>+|p>)
   !       = <s> + <p>
   d(2)=a(1)*b(2)+a(2)*b(1)
   if(la==0.or.lb==0) return
   ! <p|p> = (<s|+<p|)*(|s>+|p>)
   !       = <s> + <p> + <d>
   d(3)=a(2)*b(2)
   return
20 continue
   ! <s|d> = <s|*(|s>+|p>+|d>)
   !       = <s> + <p> + <d>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   if(la==0.or.lb==0) return
   ! <p|d> = (<s|+<p|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f>
   d(3)=d(3)+a(2)*b(2)
   d(4)=a(2)*b(3)+a(3)*b(2)
   if(la<=1.or.lb<=1) return
   ! <d|d> = (<s|+<p|+<d|)*(|s>+|p>+|d>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(5)=a(3)*b(3)
   return
30 continue
   ! <s|f> = <s|*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   if(la==0.or.lb==0) return
   ! <p|f> = (<s|+<p|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=a(2)*b(4)+a(4)*b(2)
   if(la<=1.or.lb<=1) return
   ! <d|f> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(5)=d(5)+a(3)*b(3)
   d(6)=a(3)*b(4)+a(4)*b(3)
   if(la<=2.or.lb<=2) return
   ! <f|f> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(7)=a(4)*b(4)
   return
40 continue
   ! <s|g> = <s|*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g>
   d(1)=a(1)*b(1)
   d(2)=a(1)*b(2)+a(2)*b(1)
   d(3)=a(1)*b(3)+a(3)*b(1)
   d(4)=a(1)*b(4)+a(4)*b(1)
   d(5)=a(1)*b(5)+a(5)*b(1)
   if(la==0.or.lb==0) return
   ! <p|g> = (<s|+<p|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h>
   d(3)=d(3)+a(2)*b(2)
   d(4)=d(4)+a(2)*b(3)+a(3)*b(2)
   d(5)=d(5)+a(2)*b(4)+a(4)*b(2)
   d(6)=a(2)*b(5)+a(5)*b(2)
   if(la<=1.or.lb<=1) return
   ! <d|g> = (<s|+<p|+<d|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i>
   d(5)=d(5)+a(3)*b(3)
   d(6)=d(6)+a(3)*b(4)+a(4)*b(3)
   d(7)=a(3)*b(5)+a(5)*b(3)
   if(la<=2.or.lb<=2) return
   ! <f|g> = (<s|+<p|+<d|+<f|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k>
   d(7)=d(7)+a(4)*b(4)
   d(8)=a(4)*b(5)+a(5)*b(4)
   if(la<=3.or.lb<=3) return
   ! <g|g> = (<s|+<p|+<d|+<f|+<g|)*(|s>+|p>+|d>+|f>+|g>)
   !       = <s> + <p> + <d> + <f> + <g> + <h> + <i> + <k> + <l>
   d(9)=a(5)*b(5)

end subroutine form_product


pure subroutine overlap_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3)

   v1d(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp

      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         v1d(k) = v1d(k) + s1d(l) * vv(l)
      end do
   end do

   s3d = v1d(1) * v1d(2) * v1d(3)

end subroutine overlap_3d

pure subroutine overlap_grad_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, ds3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: ds3d(3)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3)
   real(wp) :: gi(0:maxl), gg(0:maxl2), g1d(3)

   v1d(:) = 0.0_wp
   g1d(:) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rpi(k), li(k)-1, gi)
      call horizontal_shift(rpi(k), li(k)+1, gi)
      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         v1d(k) = v1d(k) + s1d(l) * vv(l)
         g1d(k) = g1d(k) + s1d(l) * gg(l)
      end do
   end do

   s3d = v1d(1) * v1d(2) * v1d(3)
   ds3d(1) = g1d(1) * v1d(2) * v1d(3)
   ds3d(2) = v1d(1) * g1d(2) * v1d(3)
   ds3d(3) = v1d(1) * v1d(2) * g1d(3)

end subroutine overlap_grad_3d

pure subroutine overlap_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call overlap_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)

end subroutine overlap_cgto

pure subroutine overlap_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, doverlap)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, grad(3), pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 1
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call overlap_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, grad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)

end subroutine overlap_grad_cgto

!> Evaluate overlap for a molecular structure
subroutine get_overlap_lat(mol, trans, cutoff, bas, overlap)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:)

   overlap(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh) &
   !$omp private(ii, jj, iao, jao, nao, r2, vec, stmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_overlap_lat

!> Evaluate overlap integrals and diatomic frame scaled overlap
subroutine get_overlap_diat_lat(mol, trans, cutoff, bas, ksig, kpi, kdel, &
   & overlap, overlap_diat)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling factors for the diatomic frame for sigma orbitals
   real(wp), intent(in) :: ksig(:, :)
   !> Scaling factors for the diatomic frame for pi orbitals
   real(wp), intent(in) :: kpi(:, :)
   !> Scaling factors for the diatomic frame for delta orbitals
   real(wp), intent(in) :: kdel(:, :)
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Overlap matrix with diatomic frame scaled elements in the diatomic frame
   real(wp), intent(out) :: overlap_diat(:, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, nsi, nsj, ii, jj, ij, iao, jao, iaosh, jaosh, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), block_overlap(:, :)
   type(diat_trafo_cache) :: dt_cache

   overlap(:, :) = 0.0_wp
   overlap_diat(:, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), block_overlap(sdim(bas%maxl), sdim(bas%maxl)))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, ksig, kpi, kdel, overlap, overlap_diat) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, nsi, nsj, ii, jj, ij) &
   !$omp private(iao, jao, iaosh, jaosh, nao, r2, vec, stmp, block_overlap, dt_cache)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle

            ! Calculate pairwise overlap and multipole integrals
            block_overlap = 0.0_wp
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)
                  call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)

                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(ij)
                     end do
                  end do
               end do
            end do

            ! Skip diatomic frame scaling for the same atom
            if (r2 > epsilon(1.0_wp)) then
               ! Perform diatomic frame scaling of current block
               call setup_diat_trafo(dt_cache, vec, nsj-1, nsi-1)
               call diat_trafo(dt_cache, ksig(izp, jzp), kpi(izp, jzp), &
                  & kdel(izp, jzp), block_overlap)
            end if

            ! Distribute the diatomic frame scaled overlap elements
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao)
                     end do
                  end do
               end do
            end do

         end do
      end do
   end do

end subroutine get_overlap_diat_lat

pure subroutine dipole_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: d3d(3)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 2)
   real(wp), parameter :: s3 = sqrt(3.0_wp), s3_4 = s3 * 0.5_wp

   v1d(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp

      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
         v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpi(k)*s1d(l)) * vv(l)
      end do
   end do

   s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)

end subroutine dipole_3d

pure subroutine dipole_grad_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d, &
      & ds3d, dd3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: d3d(3)
   real(wp), intent(out) :: ds3d(3)
   real(wp), intent(out) :: dd3d(3, 3)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 2)
   real(wp) :: gi(0:maxl), gg(0:maxl2), g1d(3, 2)

   v1d(:, :) = 0.0_wp
   g1d(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rpi(k), li(k)-1, gi)
      call horizontal_shift(rpi(k), li(k)+1, gi)
      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
         v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpi(k)*s1d(l)) * vv(l)
         g1d(k, 1) = g1d(k, 1) + s1d(l) * gg(l)
         g1d(k, 2) = g1d(k, 2) + (s1d(l+1) + rpi(k)*s1d(l)) * gg(l)
      end do
   end do

   s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)

   ds3d(1) = g1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   ds3d(2) = v1d(1, 1) * g1d(2, 1) * v1d(3, 1)
   ds3d(3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 1)
   dd3d(1, 1) = g1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   dd3d(2, 1) = v1d(1, 2) * g1d(2, 1) * v1d(3, 1)
   dd3d(3, 1) = v1d(1, 2) * v1d(2, 1) * g1d(3, 1)
   dd3d(1, 2) = g1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   dd3d(2, 2) = v1d(1, 1) * g1d(2, 2) * v1d(3, 1)
   dd3d(3, 2) = v1d(1, 1) * v1d(2, 2) * g1d(3, 1)
   dd3d(1, 3) = g1d(1, 1) * v1d(2, 1) * v1d(3, 2)
   dd3d(2, 3) = v1d(1, 1) * g1d(2, 1) * v1d(3, 2)
   dd3d(3, 3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 2)

end subroutine dipole_grad_3d

pure subroutine dipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, &
      & doverlap, ddpint)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpint(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3)
   real(wp) :: pre, grad(3), ddip(3, 3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3d(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   dd3d(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 2
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call dipole_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, grad, ddip)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
               dd3d(:, :, mlj, mli) = dd3d(:, :, mlj, mli) + cc*ddip
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3d, ddpint, .true., .true.)

end subroutine dipole_grad_cgto

pure subroutine multipole_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d, q3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: d3d(3)
   real(wp), intent(out) :: q3d(6)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 3)

   v1d(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp

      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      do l = 0, li(k) + lj(k)
         v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
         v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpi(k)*s1d(l)) * vv(l)
         v1d(k, 3) = v1d(k, 3) + (s1d(l+2) + 2*rpi(k)*s1d(l+1) + rpi(k)*rpi(k)*s1d(l)) * vv(l)
      end do
   end do

   s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)
   q3d(1) = v1d(1, 3) * v1d(2, 1) * v1d(3, 1)
   q3d(2) = v1d(1, 2) * v1d(2, 2) * v1d(3, 1)
   q3d(3) = v1d(1, 1) * v1d(2, 3) * v1d(3, 1)
   q3d(4) = v1d(1, 2) * v1d(2, 1) * v1d(3, 2)
   q3d(5) = v1d(1, 1) * v1d(2, 2) * v1d(3, 2)
   q3d(6) = v1d(1, 1) * v1d(2, 1) * v1d(3, 3)

end subroutine multipole_3d

pure subroutine multipole_grad_3d(rpj, rpi, aj, ai, lj, li, s1d, s3d, d3d, q3d, &
      & ds3d, dd3d, dq3d)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: ai
   real(wp), intent(in) :: aj
   integer, intent(in) :: li(3)
   integer, intent(in) :: lj(3)
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: s3d
   real(wp), intent(out) :: d3d(3)
   real(wp), intent(out) :: q3d(6)
   real(wp), intent(out) :: ds3d(3)
   real(wp), intent(out) :: dd3d(3, 3)
   real(wp), intent(out) :: dq3d(3, 6)

   integer :: k, l
   real(wp) :: vi(0:maxl), vj(0:maxl), vv(0:maxl2), v1d(3, 3)
   real(wp) :: gi(0:maxl), gg(0:maxl2), g1d(3, 3), rpc

   v1d(:, :) = 0.0_wp
   g1d(:, :) = 0.0_wp

   do k = 1, 3
      vv(:) = 0.0_wp
      gg(:) = 0.0_wp
      vi(:) = 0.0_wp
      vj(:) = 0.0_wp
      gi(:) = 0.0_wp
      rpc = rpj(k)

      vi(li(k)) = 1.0_wp
      vj(lj(k)) = 1.0_wp
      gi(li(k)+1) = 2*ai
      if (li(k) > 0) gi(li(k)-1) = -li(k)

      call horizontal_shift(rpi(k), li(k)-1, gi)
      call horizontal_shift(rpi(k), li(k)+1, gi)
      call horizontal_shift(rpi(k), li(k), vi)
      call horizontal_shift(rpj(k), lj(k), vj)
      call form_product(vi, vj, li(k), lj(k), vv)
      call form_product(gi, vj, li(k)+1, lj(k), gg)
      do l = 0, li(k) + lj(k) + 1
         v1d(k, 1) = v1d(k, 1) + s1d(l) * vv(l)
         v1d(k, 2) = v1d(k, 2) + (s1d(l+1) + rpc*s1d(l)) * vv(l)
         v1d(k, 3) = v1d(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * vv(l)
         g1d(k, 1) = g1d(k, 1) + s1d(l) * gg(l)
         g1d(k, 2) = g1d(k, 2) + (s1d(l+1) + rpc*s1d(l)) * gg(l)
         g1d(k, 3) = g1d(k, 3) + (s1d(l+2) + 2*rpc*s1d(l+1) + rpc*rpc*s1d(l)) * gg(l)
      end do
   end do

   s3d = v1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   d3d(1) = v1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   d3d(2) = v1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   d3d(3) = v1d(1, 1) * v1d(2, 1) * v1d(3, 2)
   q3d(1) = v1d(1, 3) * v1d(2, 1) * v1d(3, 1)
   q3d(2) = v1d(1, 2) * v1d(2, 2) * v1d(3, 1)
   q3d(3) = v1d(1, 1) * v1d(2, 3) * v1d(3, 1)
   q3d(4) = v1d(1, 2) * v1d(2, 1) * v1d(3, 2)
   q3d(5) = v1d(1, 1) * v1d(2, 2) * v1d(3, 2)
   q3d(6) = v1d(1, 1) * v1d(2, 1) * v1d(3, 3)

   ds3d(1) = g1d(1, 1) * v1d(2, 1) * v1d(3, 1)
   ds3d(2) = v1d(1, 1) * g1d(2, 1) * v1d(3, 1)
   ds3d(3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 1)
   dd3d(1, 1) = g1d(1, 2) * v1d(2, 1) * v1d(3, 1)
   dd3d(2, 1) = v1d(1, 2) * g1d(2, 1) * v1d(3, 1)
   dd3d(3, 1) = v1d(1, 2) * v1d(2, 1) * g1d(3, 1)
   dd3d(1, 2) = g1d(1, 1) * v1d(2, 2) * v1d(3, 1)
   dd3d(2, 2) = v1d(1, 1) * g1d(2, 2) * v1d(3, 1)
   dd3d(3, 2) = v1d(1, 1) * v1d(2, 2) * g1d(3, 1)
   dd3d(1, 3) = g1d(1, 1) * v1d(2, 1) * v1d(3, 2)
   dd3d(2, 3) = v1d(1, 1) * g1d(2, 1) * v1d(3, 2)
   dd3d(3, 3) = v1d(1, 1) * v1d(2, 1) * g1d(3, 2)
   dq3d(1, 1) = g1d(1, 3) * v1d(2, 1) * v1d(3, 1)
   dq3d(2, 1) = v1d(1, 3) * g1d(2, 1) * v1d(3, 1)
   dq3d(3, 1) = v1d(1, 3) * v1d(2, 1) * g1d(3, 1)
   dq3d(1, 2) = g1d(1, 2) * v1d(2, 2) * v1d(3, 1)
   dq3d(2, 2) = v1d(1, 2) * g1d(2, 2) * v1d(3, 1)
   dq3d(3, 2) = v1d(1, 2) * v1d(2, 2) * g1d(3, 1)
   dq3d(1, 3) = g1d(1, 1) * v1d(2, 3) * v1d(3, 1)
   dq3d(2, 3) = v1d(1, 1) * g1d(2, 3) * v1d(3, 1)
   dq3d(3, 3) = v1d(1, 1) * v1d(2, 3) * g1d(3, 1)
   dq3d(1, 4) = g1d(1, 2) * v1d(2, 1) * v1d(3, 2)
   dq3d(2, 4) = v1d(1, 2) * g1d(2, 1) * v1d(3, 2)
   dq3d(3, 4) = v1d(1, 2) * v1d(2, 1) * g1d(3, 2)
   dq3d(1, 5) = g1d(1, 1) * v1d(2, 2) * v1d(3, 2)
   dq3d(2, 5) = v1d(1, 1) * g1d(2, 2) * v1d(3, 2)
   dq3d(3, 5) = v1d(1, 1) * v1d(2, 2) * g1d(3, 2)
   dq3d(1, 6) = g1d(1, 1) * v1d(2, 1) * v1d(3, 3)
   dq3d(2, 6) = v1d(1, 1) * g1d(2, 1) * v1d(3, 3)
   dq3d(3, 6) = v1d(1, 1) * v1d(2, 1) * g1d(3, 3)

end subroutine multipole_grad_3d

pure subroutine shift_operator(vec, s, di, qi, ds, ddi, dqi, ddj, dqj)
   real(wp),intent(in) :: vec(:)
   real(wp),intent(in) :: s
   real(wp),intent(in) :: di(:)
   real(wp),intent(in) :: qi(:)
   real(wp),intent(in) :: ds(:)
   real(wp),intent(in) :: ddi(:, :)
   real(wp),intent(in) :: dqi(:, :)
   real(wp),intent(out) :: ddj(:, :)
   real(wp),intent(out) :: dqj(:, :)

   ddj(:, 1) = ddi(:, 1) - vec(1)*ds
   ddj(:, 2) = ddi(:, 2) - vec(2)*ds
   ddj(:, 3) = ddi(:, 3) - vec(3)*ds
   ddj(1, 1) = ddj(1, 1) - s
   ddj(2, 2) = ddj(2, 2) - s
   ddj(3, 3) = ddj(3, 3) - s

   dqj(:, 1) = dqi(:, 1) - 2*vec(1)*ddi(:, 1) + vec(1)**2*ds
   dqj(:, 3) = dqi(:, 3) - 2*vec(2)*ddi(:, 2) + vec(2)**2*ds
   dqj(:, 6) = dqi(:, 6) - 2*vec(3)*ddi(:, 3) + vec(3)**2*ds
   dqj(:, 2) = dqi(:, 2) - vec(1)*ddi(:, 2) - vec(2)*ddi(:, 1) + vec(1)*vec(2)*ds
   dqj(:, 4) = dqi(:, 4) - vec(1)*ddi(:, 3) - vec(3)*ddi(:, 1) + vec(1)*vec(3)*ds
   dqj(:, 5) = dqi(:, 5) - vec(2)*ddi(:, 3) - vec(3)*ddi(:, 2) + vec(2)*vec(3)*ds
   dqj(1, 1) = dqj(1, 1) - 2*di(1) + 2*vec(1)*s
   dqj(2, 3) = dqj(2, 3) - 2*di(2) + 2*vec(2)*s
   dqj(3, 6) = dqj(3, 6) - 2*di(3) + 2*vec(3)*s
   dqj(1, 2) = dqj(1, 2) - di(2) + vec(2)*s
   dqj(2, 2) = dqj(2, 2) - di(1) + vec(1)*s
   dqj(1, 4) = dqj(1, 4) - di(3) + vec(3)*s
   dqj(3, 4) = dqj(3, 4) - di(1) + vec(1)*s
   dqj(2, 5) = dqj(2, 5) - di(3) + vec(3)*s
   dqj(3, 5) = dqj(3, 5) - di(2) + vec(2)*s

end subroutine shift_operator


pure subroutine dipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), pre
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 1
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call dipole_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint, .true., .true.)

end subroutine dipole_cgto


subroutine get_dipole_integrals_lat(mol, trans, cutoff, bas, overlap, dpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, dpint) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj) &
   !$omp private(iao, jao, nao, r2, vec, stmp, dtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call dipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_dipole_integrals_lat


subroutine get_dipole_integrals_diat_lat(mol, trans, cutoff, bas, &
   & ksig, kpi, kdel, overlap, overlap_diat, dpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling factors for the diatomic frame for sigma orbitals
   real(wp), intent(in) :: ksig(:, :)
   !> Scaling factors for the diatomic frame for pi orbitals
   real(wp), intent(in) :: kpi(:, :)
   !> Scaling factors for the diatomic frame for delta orbitals
   real(wp), intent(in) :: kdel(:, :)
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Overlap matrix with diatomic frame scaled elements
   real(wp), intent(out) :: overlap_diat(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, nsi, nsj, ii, jj, ij, iao, jao, iaosh, jaosh, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), block_overlap(:, :)
   type(diat_trafo_cache) :: dt_cache

   overlap(:, :) = 0.0_wp
   overlap_diat(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), &
      & block_overlap(sdim(bas%maxl), sdim(bas%maxl)))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) shared(mol, bas, trans)&
   !$omp shared(cutoff2, ksig, kpi, kdel, overlap, overlap_diat, dpint) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, nsi, nsj) &
   !$omp private(ii, jj, ij, iao, jao, iaosh, jaosh, nao, r2, vec) &
   !$omp private(stmp, dtmp, block_overlap, dt_cache)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle

            ! Calculate pairwise overlap and dipole integrals
            block_overlap = 0.0_wp
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)
                  call dipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)

                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(ij)

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, ij)
                     end do
                  end do
               end do
            end do

            ! Skip diatomic frame scaling for the same atom
            if (r2 > epsilon(1.0_wp)) then
               ! Perform diatomic frame scaling of current block
               call setup_diat_trafo(dt_cache, vec, nsj-1, nsi-1)
               call diat_trafo(dt_cache, ksig(izp, jzp), kpi(izp, jzp), &
                  & kdel(izp, jzp), block_overlap)
            end if

            ! Distribute the diatomic frame scaled overlap elements
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao)
                     end do
                  end do
               end do
            end do

         end do
      end do
   end do

end subroutine get_dipole_integrals_diat_lat


pure subroutine multipole_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6), pre, tr
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3d(6, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3d(:, :, :) = 0.0_wp
   q3d(:, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 2
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3d(:, mlj, mli) = d3d(:, mlj, mli) + cc*dip
               q3d(:, mlj, mli) = q3d(:, mlj, mli) + cc*quad
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3d, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, q3d, qpint, .true., .true.)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

end subroutine multipole_cgto


pure subroutine multipole_grad_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint, &
      & doverlap, ddpintj, dqpintj, ddpinti, dqpinti)
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpinti(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpinti(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpintj(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpintj(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2), rpi(3), rpj(3), cc, val, dip(3), quad(6)
   real(wp) :: pre, grad(3), ddip(3, 3), dquad(3, 6), tr, dtr(3)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3dj(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3dj(6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3di(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3di(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3dj(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3dj(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   d3dj(:, :, :) = 0.0_wp
   q3dj(:, :, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   dd3di(:, :, :, :) = 0.0_wp
   dq3di(:, :, :, :) = 0.0_wp
   dd3dj(:, :, :, :) = 0.0_wp
   dq3dj(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 3
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_grad_3d(rpj, rpi, cgtoj%alpha(jp), cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, val, dip, quad, grad, ddip, dquad)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val
               d3dj(:, mlj, mli) = d3dj(:, mlj, mli) + cc*dip
               q3dj(:, mlj, mli) = q3dj(:, mlj, mli) + cc*quad
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad
               dd3dj(:, :, mlj, mli) = dd3dj(:, :, mlj, mli) + cc*ddip
               dq3dj(:, :, mlj, mli) = dq3dj(:, :, mlj, mli) + cc*dquad
            end do
         end do
      end do
   end do

   do mli = 1, mlao(cgtoi%ang)
      do mlj = 1, mlao(cgtoj%ang)
         call shift_operator(vec, s3d(mlj, mli), d3dj(:, mlj, mli), q3dj(:, mlj, mli), &
            & ds3d(:, mlj, mli), dd3dj(:, :, mlj, mli), dq3dj(:, :, mlj, mli), &
            & dd3di(:, :, mlj, mli), dq3di(:, :, mlj, mli))
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3dj, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, q3dj, qpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3dj, ddpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3dj, dqpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3di, ddpinti, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3di, dqpinti, .true., .true.)

   ! remove trace from quadrupole integrals (transfrom to spherical harmonics and back)
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         tr = 0.5_wp * (qpint(1, mlj, mli) + qpint(3, mlj, mli) + qpint(6, mlj, mli))
         qpint(1, mlj, mli) = 1.5_wp * qpint(1, mlj, mli) - tr
         qpint(2, mlj, mli) = 1.5_wp * qpint(2, mlj, mli)
         qpint(3, mlj, mli) = 1.5_wp * qpint(3, mlj, mli) - tr
         qpint(4, mlj, mli) = 1.5_wp * qpint(4, mlj, mli)
         qpint(5, mlj, mli) = 1.5_wp * qpint(5, mlj, mli)
         qpint(6, mlj, mli) = 1.5_wp * qpint(6, mlj, mli) - tr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpinti(:, 1, mlj, mli) + dqpinti(:, 3, mlj, mli) + dqpinti(:, 6, mlj, mli)
         dqpinti(:, 1, mlj, mli) = 1.5_wp * dqpinti(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 2, mlj, mli) = 1.5_wp * dqpinti(:, 2, mlj, mli)
         dqpinti(:, 3, mlj, mli) = 1.5_wp * dqpinti(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpinti(:, 4, mlj, mli) = 1.5_wp * dqpinti(:, 4, mlj, mli)
         dqpinti(:, 5, mlj, mli) = 1.5_wp * dqpinti(:, 5, mlj, mli)
         dqpinti(:, 6, mlj, mli) = 1.5_wp * dqpinti(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do

   ! remove trace from quadrupole integral derivatives
   do mli = 1, msao(cgtoi%ang)
      do mlj = 1, msao(cgtoj%ang)
         dtr = dqpintj(:, 1, mlj, mli) + dqpintj(:, 3, mlj, mli) + dqpintj(:, 6, mlj, mli)
         dqpintj(:, 1, mlj, mli) = 1.5_wp * dqpintj(:, 1, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 2, mlj, mli) = 1.5_wp * dqpintj(:, 2, mlj, mli)
         dqpintj(:, 3, mlj, mli) = 1.5_wp * dqpintj(:, 3, mlj, mli) - 0.5_wp * dtr
         dqpintj(:, 4, mlj, mli) = 1.5_wp * dqpintj(:, 4, mlj, mli)
         dqpintj(:, 5, mlj, mli) = 1.5_wp * dqpintj(:, 5, mlj, mli)
         dqpintj(:, 6, mlj, mli) = 1.5_wp * dqpintj(:, 6, mlj, mli) - 0.5_wp * dtr
      end do
   end do

end subroutine multipole_grad_cgto


subroutine get_multipole_integrals_lat(mol, trans, cutoff, bas, overlap, dpint, qpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, ii, jj, iao, jao, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :)

   overlap(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, bas, trans, cutoff2, overlap, dpint, qpint) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj) &
   !$omp private(iao, jao, nao, r2, vec, stmp, dtmp, qtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(jao + nao*(iao-1))

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, jao + nao*(iao-1))

                        qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                           & + qtmp(:, jao + nao*(iao-1))
                     end do
                  end do

               end do
            end do

         end do
      end do
   end do

end subroutine get_multipole_integrals_lat


subroutine get_multipole_integrals_diat_lat(mol, trans, cutoff, bas, &
   & ksig, kpi, kdel, overlap, overlap_diat, dpint, qpint)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Scaling factors for the diatomic frame for sigma orbitals
   real(wp), intent(in) :: ksig(:, :)
   !> Scaling factors for the diatomic frame for pi orbitals
   real(wp), intent(in) :: kpi(:, :)
   !> Scaling factors for the diatomic frame for delta orbitals
   real(wp), intent(in) :: kdel(:, :)
   !> Overlap matrix
   real(wp), intent(out) :: overlap(:, :)
   !> Overlap matrix with diatomic frame scaled elements
   real(wp), intent(out) :: overlap_diat(:, :)
   !> Dipole moment integral matrix
   real(wp), intent(out) :: dpint(:, :, :)
   !> Quadrupole moment integral matrix
   real(wp), intent(out) :: qpint(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js
   integer :: ish, jsh, nsi, nsj, ii, jj, ij, iao, jao, iaosh, jaosh, nao
   real(wp) :: r2, vec(3), cutoff2
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), block_overlap(:, :)
   type(diat_trafo_cache) :: dt_cache

   overlap(:, :) = 0.0_wp
   overlap_diat(:, :) = 0.0_wp
   dpint(:, :, :) = 0.0_wp
   qpint(:, :, :) = 0.0_wp

   allocate(stmp(msao(bas%maxl)**2), dtmp(3, msao(bas%maxl)**2), qtmp(6, msao(bas%maxl)**2), &
      & block_overlap(sdim(bas%maxl), sdim(bas%maxl)))
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) shared(mol, bas, trans) &
   !$omp shared(cutoff2, ksig, kpi, kdel, overlap, overlap_diat, dpint, qpint) &
   !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, nsi, nsj) &
   !$omp private(ii, jj, ij, iao, jao, iaosh, jaosh, nao, r2, vec) &
   !$omp private(stmp, dtmp, qtmp, block_overlap, dt_cache)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      nsi = bas%nsh_id(izp)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         nsj = bas%nsh_id(jzp)
         do itr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle

            ! Calculate pairwise overlap and multipole integrals
            block_overlap = 0.0_wp
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)
                  call multipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        block_overlap(jaosh+jao, iaosh+iao) = stmp(ij)

                        overlap(jj+jao, ii+iao) = overlap(jj+jao, ii+iao) &
                           & + stmp(ij)

                        dpint(:, jj+jao, ii+iao) = dpint(:, jj+jao, ii+iao) &
                           & + dtmp(:, ij)

                        qpint(:, jj+jao, ii+iao) = qpint(:, jj+jao, ii+iao) &
                           & + qtmp(:, ij)
                     end do
                  end do
               end do
            end do

            ! Skip diatomic frame scaling for the same atom
            if (r2 > epsilon(1.0_wp)) then
               ! Perform diatomic frame scaling of current block
               call setup_diat_trafo(dt_cache, vec, nsj-1, nsi-1)
               call diat_trafo(dt_cache, ksig(izp, jzp), kpi(izp, jzp), &
                  & kdel(izp, jzp), block_overlap)
            end if

            ! Distribute the diatomic frame scaled overlap elements
            do ish = 1, nsi
               ii = bas%iao_sh(is+ish)
               iaosh = smap(ish-1)
               do jsh = 1, nsj
                  jj = bas%iao_sh(js+jsh)
                  jaosh = smap(jsh-1)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  !$omp simd collapse(2)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)

                        overlap_diat(jj+jao, ii+iao) = overlap_diat(jj+jao, ii+iao) &
                           + block_overlap(jaosh+jao, iaosh+iao)
                     end do
                  end do
               end do
            end do

         end do
      end do
   end do

end subroutine get_multipole_integrals_diat_lat

pure subroutine add_shifted_monomial(ae, l, scale, cfs)
   !***********************************
   !* Add the coefficients of scale*(x + ae)**l to the polynomial cfs,
   !* for any l within the bounds of cfs. Negative l is a no-op.
   !***********************************
   real(wp), intent(in) :: ae
   integer, intent(in) :: l
   real(wp), intent(in) :: scale
   real(wp), intent(inout) :: cfs(0:)

   integer :: m
   real(wp) :: bin, pw

   if (l < 0) return
   bin = 1.0_wp
   pw = scale
   do m = l, 0, -1
      cfs(m) = cfs(m) + bin * pw
      bin = bin * real(m, wp) / real(l - m + 1, wp)
      pw = pw * ae
   end do
end subroutine add_shifted_monomial

pure subroutine multipole_hess_1d(rpj, rpi, ai, lj, li, s1d, v1d, g1d, h1d)
   !***********************************
   !* One-dimensional factor of a primitive pair with its first and second
   !* derivative w.r.t. center i. The three entries are the moment orders
   !* 0, 1 and 2 of the operator about center j.
   !***********************************
   real(wp), intent(in) :: rpj
   real(wp), intent(in) :: rpi
   real(wp), intent(in) :: ai
   integer, intent(in) :: lj
   integer, intent(in) :: li
   real(wp), intent(in) :: s1d(0:)
   real(wp), intent(out) :: v1d(3)
   real(wp), intent(out) :: g1d(3)
   real(wp), intent(out) :: h1d(3)

   integer :: l, m
   real(wp) :: vi(0:maxl+2), gi(0:maxl+2), hi(0:maxl+2), vj(0:maxl)
   real(wp) :: vv(0:maxl2+2), gg(0:maxl2+2), hh(0:maxl2+2), mom(3)

   vi(:) = 0.0_wp
   gi(:) = 0.0_wp
   hi(:) = 0.0_wp
   vj(:) = 0.0_wp
   vv(:) = 0.0_wp
   gg(:) = 0.0_wp
   hh(:) = 0.0_wp

   call add_shifted_monomial(rpj, lj, 1.0_wp, vj)
   call add_shifted_monomial(rpi, li, 1.0_wp, vi)
   !>-- d/dRi x**l = 2a x**(l+1) - l x**(l-1)
   call add_shifted_monomial(rpi, li+1, 2*ai, gi)
   call add_shifted_monomial(rpi, li-1, -real(li, wp), gi)
   !>-- d2/dRi2 x**l = 4a**2 x**(l+2) - 2a(2l+1) x**l + l(l-1) x**(l-2)
   call add_shifted_monomial(rpi, li+2, 4*ai*ai, hi)
   call add_shifted_monomial(rpi, li, -2*ai*real(2*li+1, wp), hi)
   call add_shifted_monomial(rpi, li-2, real(li*(li-1), wp), hi)

   do m = 0, lj
      do l = 0, li + 2
         vv(l+m) = vv(l+m) + vi(l) * vj(m)
         gg(l+m) = gg(l+m) + gi(l) * vj(m)
         hh(l+m) = hh(l+m) + hi(l) * vj(m)
      end do
   end do

   v1d(:) = 0.0_wp
   g1d(:) = 0.0_wp
   h1d(:) = 0.0_wp
   do l = 0, li + lj + 2
      mom(1) = s1d(l)
      mom(2) = s1d(l+1) + rpj*s1d(l)
      mom(3) = s1d(l+2) + 2*rpj*s1d(l+1) + rpj*rpj*s1d(l)
      v1d(:) = v1d + mom * vv(l)
      g1d(:) = g1d + mom * gg(l)
      h1d(:) = h1d + mom * hh(l)
   end do
end subroutine multipole_hess_1d

pure subroutine multipole_hess_3d(rpj, rpi, ai, lj, li, s1d, ncomp, val, grad, hess)
   !***********************************
   !* Value, gradient and Hessian w.r.t. center i of a primitive pair for the
   !* first ncomp operator components: overlap (1), dipole about center j
   !* (2-4) and raw quadrupole about center j (5-10 as xx, xy, yy, xz, yz, zz).
   !***********************************
   real(wp), intent(in) :: rpj(3)
   real(wp), intent(in) :: rpi(3)
   real(wp), intent(in) :: ai
   integer, intent(in) :: lj(3)
   integer, intent(in) :: li(3)
   real(wp), intent(in) :: s1d(0:)
   integer, intent(in) :: ncomp
   real(wp), intent(out) :: val(:)
   real(wp), intent(out) :: grad(:, :)
   real(wp), intent(out) :: hess(:, :, :)

   !> Moment order plus one in x, y and z for each operator component
   integer, parameter :: mo(3, 10) = reshape([&
      & 1, 1, 1, &
      & 2, 1, 1,  1, 2, 1,  1, 1, 2, &
      & 3, 1, 1,  2, 2, 1,  1, 3, 1,  2, 1, 2,  1, 2, 2,  1, 1, 3], shape(mo))
   integer :: k, ic, mx, my, mz
   real(wp) :: v(3, 3), g(3, 3), h(3, 3)

   do k = 1, 3
      call multipole_hess_1d(rpj(k), rpi(k), ai, lj(k), li(k), s1d, &
         & v(:, k), g(:, k), h(:, k))
   end do

   do ic = 1, ncomp
      mx = mo(1, ic)
      my = mo(2, ic)
      mz = mo(3, ic)
      val(ic) = v(mx, 1) * v(my, 2) * v(mz, 3)
      grad(1, ic) = g(mx, 1) * v(my, 2) * v(mz, 3)
      grad(2, ic) = v(mx, 1) * g(my, 2) * v(mz, 3)
      grad(3, ic) = v(mx, 1) * v(my, 2) * g(mz, 3)
      hess(1, 1, ic) = h(mx, 1) * v(my, 2) * v(mz, 3)
      hess(2, 2, ic) = v(mx, 1) * h(my, 2) * v(mz, 3)
      hess(3, 3, ic) = v(mx, 1) * v(my, 2) * h(mz, 3)
      hess(2, 1, ic) = g(mx, 1) * g(my, 2) * v(mz, 3)
      hess(3, 1, ic) = g(mx, 1) * v(my, 2) * g(mz, 3)
      hess(3, 2, ic) = v(mx, 1) * g(my, 2) * g(mz, 3)
      hess(1, 2, ic) = hess(2, 1, ic)
      hess(1, 3, ic) = hess(3, 1, ic)
      hess(2, 3, ic) = hess(3, 2, ic)
   end do
end subroutine multipole_hess_3d

pure subroutine shift_operator_hess(vec, s, ds, d2s, ddj, d2dj, d2qj, d2di, d2qi)
   !***********************************
   !* Move the second derivatives of the dipole and raw quadrupole integrals
   !* from the origin at center j to the origin at center i, which moves
   !* with the differentiated center. Second-order analogue of shift_operator.
   !***********************************
   real(wp), intent(in) :: vec(3)
   real(wp), intent(in) :: s
   real(wp), intent(in) :: ds(3)
   real(wp), intent(in) :: d2s(3, 3)
   real(wp), intent(in) :: ddj(3, 3)
   real(wp), intent(in) :: d2dj(3, 3, 3)
   real(wp), intent(in) :: d2qj(3, 3, 6)
   real(wp), intent(out) :: d2di(3, 3, 3)
   real(wp), intent(out) :: d2qi(3, 3, 6)

   integer, parameter :: cmap(6) = [1, 1, 2, 1, 2, 3], dmap(6) = [1, 2, 2, 3, 3, 3]
   real(wp), parameter :: unity(3, 3) = reshape([1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])
   integer :: a, b, c, d, ic

   do c = 1, 3
      do b = 1, 3
         do a = 1, 3
            d2di(a, b, c) = d2dj(a, b, c) - vec(c)*d2s(a, b) &
               & - unity(a, c)*ds(b) - unity(b, c)*ds(a)
         end do
      end do
   end do

   do ic = 1, 6
      c = cmap(ic)
      d = dmap(ic)
      do b = 1, 3
         do a = 1, 3
            d2qi(a, b, ic) = d2qj(a, b, ic) &
               & - vec(d)*d2dj(a, b, c) - vec(c)*d2dj(a, b, d) &
               & - unity(a, d)*ddj(b, c) - unity(b, d)*ddj(a, c) &
               & - unity(a, c)*ddj(b, d) - unity(b, c)*ddj(a, d) &
               & + (unity(a, c)*unity(b, d) + unity(a, d)*unity(b, c))*s &
               & + (unity(a, c)*vec(d) + unity(a, d)*vec(c))*ds(b) &
               & + (unity(b, c)*vec(d) + unity(b, d)*vec(c))*ds(a) &
               & + vec(c)*vec(d)*d2s(a, b)
         end do
      end do
   end do
end subroutine shift_operator_hess

pure subroutine remove_trace(qp)
   !***********************************
   !* Project the trace out of one set of Cartesian quadrupole components
   !* ordered as xx, xy, yy, xz, yz, zz.
   !***********************************
   real(wp), intent(inout) :: qp(:)

   real(wp) :: tr

   tr = 0.5_wp * (qp(1) + qp(3) + qp(6))
   qp(:) = 1.5_wp * qp
   qp(1) = qp(1) - tr
   qp(3) = qp(3) - tr
   qp(6) = qp(6) - tr
end subroutine remove_trace

pure subroutine overlap_hess_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, doverlap, &
      & d2overlap)
   !***********************************
   !* Overlap integrals of a shell pair with their first and second
   !* derivative w.r.t. the position of center i.
   !***********************************
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2overlap(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l
   real(wp) :: eab, oab, est, s1d(0:maxl2+4), rpi(3), rpj(3), cc, pre
   real(wp) :: val(1), grad(3, 1), hess(3, 3, 1)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2s3d(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))

   s3d(:, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   d2s3d(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 4
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, mlao(cgtoi%ang)
            do mlj = 1, mlao(cgtoj%ang)
               call multipole_hess_3d(rpj, rpi, cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, 1, val, grad, hess)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val(1)
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad(:, 1)
               d2s3d(:, :, mlj, mli) = d2s3d(:, :, mlj, mli) + cc*hess(:, :, 1)
            end do
         end do
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, d2s3d, d2overlap, .true., .true.)

end subroutine overlap_hess_cgto

pure subroutine multipole_hess_cgto(cgtoj, cgtoi, r2, vec, intcut, overlap, dpint, qpint, &
      & doverlap, ddpintj, dqpintj, ddpinti, dqpinti, d2overlap, d2dpintj, d2qpintj, &
      & d2dpinti, d2qpinti)
   !***********************************
   !* Overlap, dipole and traceless quadrupole integrals of a shell pair with
   !* their first and second derivative w.r.t. the position of center i.
   !* Values, first derivatives and the j/i naming of the derivative arrays
   !* follow multipole_grad_cgto. The second derivatives are indexed as
   !* (direction, direction, operator component, ao j, ao i).
   !***********************************
   !> Description of contracted Gaussian function on center j
   type(cgto_type), intent(in) :: cgtoj
   !> Description of contracted Gaussian function on center i
   type(cgto_type), intent(in) :: cgtoi
   !> Square distance between center i and j
   real(wp), intent(in) :: r2
   !> Distance vector between center i and j, ri - rj
   real(wp), intent(in) :: vec(3)
   !> Maximum value of integral prefactor to consider
   real(wp), intent(in) :: intcut
   !> Overlap integrals for the given pair i  and j
   real(wp), intent(out) :: overlap(msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integrals for the given pair i  and j
   real(wp), intent(out) :: dpint(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integrals for the given pair i  and j
   real(wp), intent(out) :: qpint(6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral gradient for the given pair i  and j
   real(wp), intent(out) :: doverlap(3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpintj(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpintj(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: ddpinti(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral gradient for the given pair i  and j
   real(wp), intent(out) :: dqpinti(3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Overlap integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2overlap(3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2dpintj(3, 3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2qpintj(3, 3, 6, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Dipole moment integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2dpinti(3, 3, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   !> Quadrupole moment integral Hessian for the given pair i  and j
   real(wp), intent(out) :: d2qpinti(3, 3, 6, msao(cgtoj%ang), msao(cgtoi%ang))

   integer :: ip, jp, mli, mlj, l, a, b, nci, ncj, nsi, nsj
   real(wp) :: eab, oab, est, s1d(0:maxl2+4), rpi(3), rpj(3), cc, pre
   real(wp) :: val(10), grad(3, 10), hess(3, 3, 10), d2di(3, 3, 3), d2qi(3, 3, 6)
   real(wp) :: s3d(mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d3dj(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: q3dj(6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: ds3d(3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3di(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3di(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dd3dj(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: dq3dj(3, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2s3d(3, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   !>-- the two directions share the leading dimension, as transform2 expects
   real(wp) :: d2d3dj(9, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2q3dj(9, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2d3di(9, 3, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2q3di(9, 6, mlao(cgtoj%ang), mlao(cgtoi%ang))
   real(wp) :: d2dsph(9, 3, msao(cgtoj%ang), msao(cgtoi%ang))
   real(wp) :: d2qsph(9, 6, msao(cgtoj%ang), msao(cgtoi%ang))

   ncj = mlao(cgtoj%ang)
   nci = mlao(cgtoi%ang)
   nsj = msao(cgtoj%ang)
   nsi = msao(cgtoi%ang)

   s3d(:, :) = 0.0_wp
   d3dj(:, :, :) = 0.0_wp
   q3dj(:, :, :) = 0.0_wp
   ds3d(:, :, :) = 0.0_wp
   dd3dj(:, :, :, :) = 0.0_wp
   dq3dj(:, :, :, :) = 0.0_wp
   d2s3d(:, :, :, :) = 0.0_wp
   d2d3dj(:, :, :, :) = 0.0_wp
   d2q3dj(:, :, :, :) = 0.0_wp

   do ip = 1, cgtoi%nprim
      do jp = 1, cgtoj%nprim
         eab = cgtoi%alpha(ip) + cgtoj%alpha(jp)
         oab = 1.0_wp/eab
         est = cgtoi%alpha(ip) * cgtoj%alpha(jp) * r2 * oab
         if (est > intcut) cycle
         pre = exp(-est) * sqrtpi3*sqrt(oab)**3
         rpi = -vec * cgtoj%alpha(jp) * oab
         rpj = +vec * cgtoi%alpha(ip) * oab
         do l = 0, cgtoi%ang + cgtoj%ang + 4
            s1d(l) = overlap_1d(l, eab)
         end do
         cc = cgtoi%coeff(ip) * cgtoj%coeff(jp) * pre
         do mli = 1, nci
            do mlj = 1, ncj
               call multipole_hess_3d(rpj, rpi, cgtoi%alpha(ip), &
                  & lx(:, mlj+lmap(cgtoj%ang)), lx(:, mli+lmap(cgtoi%ang)), &
                  & s1d, 10, val, grad, hess)
               s3d(mlj, mli) = s3d(mlj, mli) + cc*val(1)
               d3dj(:, mlj, mli) = d3dj(:, mlj, mli) + cc*val(2:4)
               q3dj(:, mlj, mli) = q3dj(:, mlj, mli) + cc*val(5:10)
               ds3d(:, mlj, mli) = ds3d(:, mlj, mli) + cc*grad(:, 1)
               dd3dj(:, :, mlj, mli) = dd3dj(:, :, mlj, mli) + cc*grad(:, 2:4)
               dq3dj(:, :, mlj, mli) = dq3dj(:, :, mlj, mli) + cc*grad(:, 5:10)
               d2s3d(:, :, mlj, mli) = d2s3d(:, :, mlj, mli) + cc*hess(:, :, 1)
               d2d3dj(:, :, mlj, mli) = d2d3dj(:, :, mlj, mli) &
                  & + cc*reshape(hess(:, :, 2:4), [9, 3])
               d2q3dj(:, :, mlj, mli) = d2q3dj(:, :, mlj, mli) &
                  & + cc*reshape(hess(:, :, 5:10), [9, 6])
            end do
         end do
      end do
   end do

   do mli = 1, nci
      do mlj = 1, ncj
         call shift_operator(vec, s3d(mlj, mli), d3dj(:, mlj, mli), q3dj(:, mlj, mli), &
            & ds3d(:, mlj, mli), dd3dj(:, :, mlj, mli), dq3dj(:, :, mlj, mli), &
            & dd3di(:, :, mlj, mli), dq3di(:, :, mlj, mli))
         call shift_operator_hess(vec, s3d(mlj, mli), ds3d(:, mlj, mli), &
            & d2s3d(:, :, mlj, mli), dd3dj(:, :, mlj, mli), &
            & reshape(d2d3dj(:, :, mlj, mli), [3, 3, 3]), &
            & reshape(d2q3dj(:, :, mlj, mli), [3, 3, 6]), d2di, d2qi)
         d2d3di(:, :, mlj, mli) = reshape(d2di, [9, 3])
         d2q3di(:, :, mlj, mli) = reshape(d2qi, [9, 6])
      end do
   end do

   call transform0(cgtoj%ang, cgtoi%ang, s3d, overlap, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, d3dj, dpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, q3dj, qpint, .true., .true.)
   call transform1(cgtoj%ang, cgtoi%ang, ds3d, doverlap, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3dj, ddpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3dj, dqpintj, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dd3di, ddpinti, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, dq3di, dqpinti, .true., .true.)
   call transform2(cgtoj%ang, cgtoi%ang, d2s3d, d2overlap, .true., .true.)

   call transform2(cgtoj%ang, cgtoi%ang, d2d3dj, d2dsph, .true., .true.)
   d2dpintj(:, :, :, :, :) = reshape(d2dsph, [3, 3, 3, nsj, nsi])
   call transform2(cgtoj%ang, cgtoi%ang, d2d3di, d2dsph, .true., .true.)
   d2dpinti(:, :, :, :, :) = reshape(d2dsph, [3, 3, 3, nsj, nsi])
   call transform2(cgtoj%ang, cgtoi%ang, d2q3dj, d2qsph, .true., .true.)
   d2qpintj(:, :, :, :, :) = reshape(d2qsph, [3, 3, 6, nsj, nsi])
   call transform2(cgtoj%ang, cgtoi%ang, d2q3di, d2qsph, .true., .true.)
   d2qpinti(:, :, :, :, :) = reshape(d2qsph, [3, 3, 6, nsj, nsi])

   !>-- remove trace from quadrupole integrals and their derivatives
   do mli = 1, nsi
      do mlj = 1, nsj
         call remove_trace(qpint(:, mlj, mli))
         do a = 1, 3
            call remove_trace(dqpintj(a, :, mlj, mli))
            call remove_trace(dqpinti(a, :, mlj, mli))
            do b = 1, 3
               call remove_trace(d2qpintj(b, a, :, mlj, mli))
               call remove_trace(d2qpinti(b, a, :, mlj, mli))
            end do
         end do
      end do
   end do

end subroutine multipole_hess_cgto

end module tblite_integral_native_integrals
