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
!> Provides the second derivative of the multipole integrals with the potential

!> Second derivative of the interaction of the density with the multipole
!> potentials of the anisotropic electrostatics via the multipole integrals
module tblite_hessian_multipole
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_hessian_pair, only : add_pair_hessian
   use tblite_integral_native_integrals, only : multipole_hess_cgto, msao
   implicit none
   private

   public :: get_multipole_band_hessian

contains

subroutine get_multipole_band_hessian(mol, bas, pmat, vdp, vqp, trans, cutoff, hessian)
   !***********************************
   !* Add the second derivative of the interaction of the density with the
   !* multipole potentials at fixed density and potentials. The potentials
   !* act via the dipole and quadrupole integrals, whose second derivative is
   !* contracted with the density and the potentials. The first derivative
   !* of the integrals with the derivative of the potentials is part of the
   !* contraction of the moments at fixed density in the response part.
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
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj
   integer :: iao, jao, nao, ij, ic, jc
   real(wp) :: r2, vec(3), cutoff2, pij, block(3, 3)
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), dstmp(:, :)
   real(wp), allocatable :: ddj(:, :, :), dqj(:, :, :), ddi(:, :, :), dqi(:, :, :)
   real(wp), allocatable :: d2s(:, :, :), d2dj(:, :, :, :), d2qj(:, :, :, :)
   real(wp), allocatable :: d2di(:, :, :, :), d2qi(:, :, :, :)

   nat = mol%nat
   cutoff2 = cutoff**2
   nao = msao(bas%maxl)**2
   allocate(stmp(nao), dtmp(3, nao), qtmp(6, nao), dstmp(3, nao), ddj(3, 3, nao), &
      & dqj(3, 6, nao), ddi(3, 3, nao), dqi(3, 6, nao), d2s(3, 3, nao), &
      & d2dj(3, 3, 3, nao), d2qj(3, 3, 6, nao), d2di(3, 3, 3, nao), d2qi(3, 3, 6, nao))

   do iat = 1, nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, nat
         if (iat == jat) cycle
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         block(:, :) = 0.0_wp
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
                     end do
                  end do
               end do
            end do
         end do

         !>-- both orders of each pair are visited
         call add_pair_hessian(hessian, iat, jat, 0.5_wp*block)
      end do
   end do
end subroutine get_multipole_band_hessian

end module tblite_hessian_multipole
