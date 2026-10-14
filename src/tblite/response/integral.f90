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

!> @file tblite/response/integral.f90
!> Provides derivatives of the integral matrices w.r.t. a single atom

!> Derivatives of the integral matrices w.r.t. the position of one atom
module tblite_response_integral
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use tblite_basis_type, only : basis_type
   use tblite_integral_native_integrals, only : overlap_grad_cgto, &
      & multipole_grad_cgto, msao
   implicit none
   private

   public :: get_overlap_derivative, get_multipole_derivative

contains

subroutine get_overlap_derivative(mol, bas, trans, cutoff, kat, doverlap)
   !***********************************
   !* Derivative of the overlap matrix w.r.t. the Cartesian position of the
   !* atom kat. Only blocks between kat and another atom are nonzero.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Atom to differentiate for
   integer, intent(in) :: kat
   !> Derivative of the overlap matrix, the last dimension is the direction
   real(wp), intent(out) :: doverlap(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij
   real(wp) :: r2, vec(3), cutoff2, sgn
   real(wp), allocatable :: stmp(:), dstmp(:, :)

   doverlap(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2
   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         if (iat == jat) cycle
         if (iat /= kat .and. jat /= kat) cycle
         !>-- the integral derivatives refer to the position of atom iat
         sgn = merge(1.0_wp, -1.0_wp, iat == kat)
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
                  call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dstmp)
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        doverlap(jj+jao, ii+iao, :) = doverlap(jj+jao, ii+iao, :) &
                           & + sgn*dstmp(:, ij)
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine get_overlap_derivative

subroutine get_multipole_derivative(mol, bas, trans, cutoff, kat, pmat, vdp, vqp, &
      & dfock, ddpat, dqpat)
   !***********************************
   !* Contributions of the derivative of the dipole and quadrupole integrals
   !* w.r.t. the Cartesian position of the atom kat. Returned are the change
   !* of the atomic multipole moments at fixed density and the derivative of
   !* the Hamiltonian at fixed multipole potentials, which is added to dfock.
   !* The last dimension of all results is the direction.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Atom to differentiate for
   integer, intent(in) :: kat
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :)
   !> Potential acting on the dipole moments
   real(wp), intent(in) :: vdp(:, :)
   !> Potential acting on the quadrupole moments
   real(wp), intent(in) :: vqp(:, :)
   !> Derivative of the Hamiltonian
   real(wp), intent(inout) :: dfock(:, :, :)
   !> Derivative of the atomic dipole moments at fixed density
   real(wp), intent(out) :: ddpat(:, :, :)
   !> Derivative of the atomic quadrupole moments at fixed density
   real(wp), intent(out) :: dqpat(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, ic
   real(wp) :: r2, vec(3), cutoff2, sgn, vij, pij
   real(wp), allocatable :: stmp(:), dtmp(:, :), qtmp(:, :), dstmp(:, :)
   real(wp), allocatable :: ddj(:, :, :), dqj(:, :, :), ddi(:, :, :), dqi(:, :, :)

   ddpat(:, :, :) = 0.0_wp
   dqpat(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2
   nao = msao(bas%maxl)**2
   allocate(stmp(nao), dtmp(3, nao), qtmp(6, nao), dstmp(3, nao), ddj(3, 3, nao), &
      & dqj(3, 6, nao), ddi(3, 3, nao), dqi(3, 6, nao))

   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         if (iat == jat) cycle
         if (iat /= kat .and. jat /= kat) cycle
         sgn = merge(1.0_wp, -1.0_wp, iat == kat)
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
                  call multipole_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmp, qtmp, dstmp, ddj, dqj, ddi, dqi)
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        pij = pmat(jj+jao, ii+iao)
                        !>-- the multipole operator is centered on the atom of
                        !>   the second orbital, here atom iat
                        do ic = 1, 3
                           vij = -0.5_wp*sgn*(dot_product(ddi(ic, :, ij), vdp(:, iat)) &
                              & + dot_product(dqi(ic, :, ij), vqp(:, iat)))
                           dfock(jj+jao, ii+iao, ic) = dfock(jj+jao, ii+iao, ic) + vij
                           dfock(ii+iao, jj+jao, ic) = dfock(ii+iao, jj+jao, ic) + vij
                           ddpat(:, iat, ic) = ddpat(:, iat, ic) - sgn*pij*ddi(ic, :, ij)
                           dqpat(:, iat, ic) = dqpat(:, iat, ic) - sgn*pij*dqi(ic, :, ij)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine get_multipole_derivative

end module tblite_response_integral
