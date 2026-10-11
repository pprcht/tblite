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

!> @file tblite/response/density.f90
!> Provides the uncoupled response of the density matrix

!> First-order response of the density matrix to a perturbation of the
!> Hamiltonian at fixed potential
module tblite_response_density
   use mctc_env, only : wp
   use tblite_blas, only : gemm
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: density_response, smear_occupations, fermi_shift, occtol

   !> Occupation below which an orbital is treated as empty
   real(wp), parameter :: occtol = 1.0e-12_wp
   !> Orbital energy difference below which two orbitals count as degenerate
   real(wp), parameter :: degtol = 1.0e-8_wp

contains

subroutine density_response(wfn, dh, dp)
   !***********************************
   !* Uncoupled first-order density matrix for a symmetric perturbation dh
   !* of the Hamiltonian with unchanged overlap, summed over both spin
   !* channels of a restricted wavefunction.
   !* With M = C^T dh C the response in the orbital basis is
   !* (f_p - f_q)/(e_p - e_q) M_pq, which turns into f'_p M_pq for
   !* degenerate pairs. At finite electronic temperature the shift of the
   !* Fermi level is subtracted on the diagonal to conserve the electron
   !* number.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Perturbation of the Hamiltonian in the atomic orbital basis
   real(wp), intent(in) :: dh(:, :)
   !> Response of the density matrix
   real(wp), intent(out) :: dp(:, :)

   integer :: nao, nact, i, q, spin
   real(wp) :: de, def
   real(wp), allocatable :: hmo(:, :), wmo(:, :), tmp(:, :), z(:, :), w2(:, :)
   real(wp), allocatable :: f(:, :), fp(:, :), dep(:)

   !>-- only rows of orbitals with nonzero occupation carry weight
   nao = size(dh, 1)
   call smear_occupations(wfn%focc(:, 1), wfn%focc(:, 2), wfn%kt, f, fp)
   nact = count(max(f(:, 1), f(:, 2)) > occtol)
   if (nact == 0) then
      dp(:, :) = 0.0_wp
      return
   end if

   allocate(tmp(nact, nao), hmo(nact, nao), dep(nact))
   allocate(wmo(nact, nao), source=0.0_wp)
   call gemm(wfn%coeff(:, 1:nact, 1), dh, tmp, transa="t")
   call gemm(tmp, wfn%coeff(:, :, 1), hmo)

   do i = 1, nact
      dep(i) = hmo(i, i)
   end do

   do spin = 1, 2
      def = fermi_shift(fp(1:nact, spin), dep)
      do q = 1, nao
         do i = 1, nact
            de = wfn%emo(i, 1) - wfn%emo(q, 1)
            if (i == q .or. abs(de) < degtol) then
               wmo(i, q) = wmo(i, q) + fp(i, spin) * hmo(i, q)
            else
               wmo(i, q) = wmo(i, q) + (f(i, spin) - f(q, spin)) / de * hmo(i, q)
            end if
         end do
      end do
      do i = 1, nact
         wmo(i, i) = wmo(i, i) - def * fp(i, spin)
      end do
   end do

   !>-- dp = C_act wmo C^T + transpose, less the double counted active block
   allocate(z(nact, nao), w2(nact, nao))
   call gemm(wmo, wfn%coeff(:, :, 1), z, transb="t")
   call gemm(wmo(:, 1:nact), wfn%coeff(:, 1:nact, 1), w2, transb="t")
   call gemm(wfn%coeff(:, 1:nact, 1), z, dp)
   call gemm(z, wfn%coeff(:, 1:nact, 1), dp, transa="t", transb="t", beta=1.0_wp)
   call gemm(wfn%coeff(:, 1:nact, 1), w2, dp, alpha=-1.0_wp, beta=1.0_wp)

end subroutine density_response

subroutine smear_occupations(focca, foccb, kt, f, fp)
   !***********************************
   !* Occupations of both spin channels and their derivative w.r.t. the
   !* orbital energy, f' = -f(1-f)/kT for Fermi smearing. The derivative
   !* vanishes without electronic temperature.
   !***********************************
   !> Occupation numbers of the alpha channel
   real(wp), intent(in) :: focca(:)
   !> Occupation numbers of the beta channel
   real(wp), intent(in) :: foccb(:)
   !> Electronic temperature
   real(wp), intent(in) :: kt
   !> Occupation numbers per spin channel
   real(wp), allocatable, intent(out) :: f(:, :)
   !> Derivative of the occupation numbers w.r.t. the orbital energy
   real(wp), allocatable, intent(out) :: fp(:, :)

   allocate(f(size(focca), 2), fp(size(focca), 2))
   f(:, 1) = focca
   f(:, 2) = foccb

   if (kt > 0.0_wp) then
      fp(:, :) = -f * (1.0_wp - f) / kt
   else
      fp(:, :) = 0.0_wp
   end if
end subroutine smear_occupations

pure function fermi_shift(fp, dep) result(def)
   !***********************************
   !* Response of the Fermi level which conserves the number of electrons
   !* in one spin channel, sum_p f'_p de_p / sum_p f'_p. Vanishes for
   !* integer occupations.
   !***********************************
   !> Derivative of the occupation numbers w.r.t. the orbital energy
   real(wp), intent(in) :: fp(:)
   !> Response of the orbital energies
   real(wp), intent(in) :: dep(:)
   !> Response of the Fermi level
   real(wp) :: def

   real(wp) :: sfp

   def = 0.0_wp
   sfp = sum(fp)
   if (abs(sfp) < 1.0e-14_wp) return
   def = sum(fp*dep)/sfp
end function fermi_shift

end module tblite_response_density
