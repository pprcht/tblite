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

!> @file tblite/hessian/dispersion.f90
!> Provides the second derivative of the dispersion energy

!> Analytic Hessian of the D3 and D4 dispersion corrections
module tblite_hessian_dispersion
   use dftd3, only : get_d3_dispersion => get_dispersion
   use dftd4_cutoff, only : smooth_cutoff_r2
   use dftd4_disp, only : get_dispersion3_hessian
   use dftd4_model_d4, only : d4_model
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_blas, only : gemm
   use tblite_cutoff, only : get_lattice_points
   use tblite_disp, only : dispersion_type
   use tblite_disp_d3, only : d3_dispersion
   use tblite_disp_d4, only : d4_dispersion
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian
   implicit none
   private

   public :: get_dispersion_hessian

contains

subroutine get_dispersion_hessian(disp, mol, qat, hessian, error)
   !***********************************
   !* Add the Hessian of the dispersion energy. For the self-consistent D4
   !* model the atomic charges are kept fixed, their response has to be
   !* accounted for by the caller.
   !***********************************
   !> Instance of the dispersion container
   class(dispersion_type), intent(in) :: disp
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: energies(:), tmp(:, :)

   allocate(tmp(3*mol%nat, 3*mol%nat))

   select type(disp)
   type is (d3_dispersion)
      allocate(energies(mol%nat))
      call get_d3_dispersion(error, mol, disp%model, disp%param, disp%cutoff, &
         & energies, hessian=tmp)
      if (allocated(error)) return
      hessian(:, :) = hessian + tmp
   type is (d4_dispersion)
      !>-- the three-body term is evaluated for neutral atoms
      call get_dispersion3_hessian(error, mol, disp%model, disp%param, disp%cutoff, tmp)
      if (allocated(error)) return
      hessian(:, :) = hessian + tmp
      call get_d4_twobody_hessian(disp, mol, qat, hessian, error)
   class default
      call fatal_error(error, "Hessian is not available for this dispersion correction")
   end select
end subroutine get_dispersion_hessian

subroutine get_d4_twobody_hessian(disp, mol, qat, hessian, error)
   !***********************************
   !* Add the Hessian of the two-body D4 dispersion energy at fixed charges,
   !* E = -sum_ij C6_ij(CN_i, CN_j) F(r_ij). The direct second derivative of
   !* the damped radial function F is added per pair, the dependence on the
   !* coordination numbers is collected and contracted with their derivatives.
   !* Interactions of an atom with its own periodic images are not included.
   !***********************************
   !> Instance of the dispersion container
   type(d4_dispersion), intent(in) :: disp
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, ndim, mref, iat, jat, kat, izp, jzp, itr, ii, jj
   real(wp) :: vec(3), r1, r2, cutoff2, rrij, r0ij, t6, t8, f, dfdr2, d2fdr22
   real(wp) :: sw, dswdr2, d2swdr22, fsw, df, d2f, block(3, 3), dfvec(3)
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), lattr(:, :)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwd2cn(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :), dEdcndcn(:, :), dr(:, :), tmp(:, :)

   nat = mol%nat
   ndim = 3*nat

   allocate(cn(nat), dcndr(3, nat, nat), dcndL(3, 3, nat))
   call get_lattice_points(mol%periodic, mol%lattice, disp%cutoff%cn, lattr)
   call disp%ncoord%get_coordination_number(mol, lattr, cn, dcndr, dcndL)

   select type(model => disp%model)
   type is (d4_model)
      mref = maxval(model%ref)
      allocate(gwvec(mref, nat, model%ncoup), gwdcn(mref, nat, model%ncoup), &
         & gwd2cn(mref, nat, model%ncoup))
      allocate(c6(nat, nat), dc6dcn(nat, nat), d2c6dcn2(nat, nat), d2c6dcnij(nat, nat))
      call model%weight_references_hessian(mol, cn, qat, gwvec, gwdcn, gwd2cn)
      call model%get_atomic_c6_hessian(mol, gwvec, gwdcn, gwd2cn, c6, dc6dcn, &
         & d2c6dcn2, d2c6dcnij)
   class default
      call fatal_error(error, "Dispersion Hessian is only available for the D4 model")
      return
   end select

   allocate(dEdcn(nat), dEdcndr(ndim, nat), dEdcndcn(nat, nat), source=0.0_wp)

   cutoff2 = disp%cutoff%disp2**2
   call get_lattice_points(mol%periodic, mol%lattice, disp%cutoff%disp2, lattr)
   do iat = 1, nat
      izp = mol%id(iat)
      ii = 3*(iat - 1)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         jj = 3*(jat - 1)
         rrij = 3*disp%model%r4r2(izp)*disp%model%r4r2(jzp)
         r0ij = disp%param%a1*sqrt(rrij) + disp%param%a2
         do itr = 1, size(lattr, 2)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - lattr(:, itr)
            r2 = sum(vec**2)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r1 = sqrt(r2)

            !>-- damped radial function and its derivatives w.r.t. r**2
            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)
            f = disp%param%s6*t6 + disp%param%s8*rrij*t8
            dfdr2 = -3*disp%param%s6*r2**2*t6**2 - 4*disp%param%s8*rrij*r2**3*t8**2
            d2fdr22 = disp%param%s6*(18*r2**4*t6**3 - 6*r2*t6**2) &
               & + disp%param%s8*rrij*(32*r2**6*t8**3 - 12*r2**2*t8**2)
            call smooth_cutoff_r2(r2, disp%cutoff%disp2, disp%cutoff%width2, sw, &
               & dswdr2, d2swdr22)
            fsw = sw*f
            !>-- convert to derivatives w.r.t. the distance
            df = 2*r1*(dswdr2*f + sw*dfdr2)
            d2f = 2*(dswdr2*f + sw*dfdr2) &
               & + 4*r2*(d2swdr22*f + 2*dswdr2*dfdr2 + sw*d2fdr22)

            call get_pair_block(-c6(iat, jat)*df, -c6(iat, jat)*d2f, r1, vec, block)
            call add_pair_hessian(hessian, iat, jat, block)

            dEdcn(iat) = dEdcn(iat) - fsw*dc6dcn(iat, jat)
            dEdcn(jat) = dEdcn(jat) - fsw*dc6dcn(jat, iat)

            dfvec(:) = df*vec/r1
            dEdcndr(ii+1:ii+3, iat) = dEdcndr(ii+1:ii+3, iat) - dc6dcn(iat, jat)*dfvec
            dEdcndr(jj+1:jj+3, iat) = dEdcndr(jj+1:jj+3, iat) + dc6dcn(iat, jat)*dfvec
            dEdcndr(ii+1:ii+3, jat) = dEdcndr(ii+1:ii+3, jat) - dc6dcn(jat, iat)*dfvec
            dEdcndr(jj+1:jj+3, jat) = dEdcndr(jj+1:jj+3, jat) + dc6dcn(jat, iat)*dfvec

            dEdcndcn(iat, iat) = dEdcndcn(iat, iat) - fsw*d2c6dcn2(iat, jat)
            dEdcndcn(jat, jat) = dEdcndcn(jat, jat) - fsw*d2c6dcn2(jat, iat)
            dEdcndcn(iat, jat) = dEdcndcn(iat, jat) - fsw*d2c6dcnij(iat, jat)
            dEdcndcn(jat, iat) = dEdcndcn(jat, iat) - fsw*d2c6dcnij(jat, iat)
         end do
      end do
   end do

   call get_lattice_points(mol%periodic, mol%lattice, disp%cutoff%cn, lattr)
   call disp%ncoord%add_coordination_number_hessian(mol, lattr, dEdcn, hessian)

   !>-- derivative of the coordination number of kat w.r.t. the position of iat
   allocate(dr(ndim, nat), tmp(ndim, nat))
   do kat = 1, nat
      do iat = 1, nat
         dr(3*(iat - 1)+1:3*iat, kat) = dcndr(:, iat, kat)
      end do
   end do

   call gemm(dEdcndr, dr, hessian, transb="t", beta=1.0_wp)
   call gemm(dr, dEdcndr, hessian, transb="t", beta=1.0_wp)
   call gemm(dr, dEdcndcn, tmp)
   call gemm(tmp, dr, hessian, transb="t", beta=1.0_wp)
end subroutine get_d4_twobody_hessian

end module tblite_hessian_dispersion
