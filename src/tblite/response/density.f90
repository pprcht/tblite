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

!> First-order response of the density matrix and the energy-weighted density
!> matrix to a perturbation of the Hamiltonian and the overlap at fixed potential
module tblite_response_density
   use mctc_env, only : wp
   use tblite_blas, only : gemm
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: density_response, wdensity_response, density_response_derivative
   public :: has_integer_occupation

   !> Occupation below which an orbital is treated as empty
   real(wp), parameter :: occtol = 1.0e-12_wp
   !> Orbital energy difference below which two orbitals count as degenerate
   real(wp), parameter :: degtol = 1.0e-8_wp

contains

subroutine density_response(wfn, dh, dp, ds)
   !***********************************
   !* First-order density matrix for a symmetric perturbation dh of the
   !* Hamiltonian and, optionally, ds of the overlap, at fixed potential and
   !* summed over both spin channels of a restricted wavefunction.
   !* With H = C^T dh C and S = C^T ds C the response in the orbital basis is
   !* ((f_p - f_q) H_pq - (f_p e_p - f_q e_q) S_pq)/(e_p - e_q), which turns
   !* into f'_p (H_pq - e_p S_pq) - f_p S_pq for degenerate pairs. At finite
   !* electronic temperature the shift of the Fermi level is subtracted on
   !* the diagonal to conserve the electron number.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Perturbation of the Hamiltonian in the atomic orbital basis
   real(wp), intent(in) :: dh(:, :)
   !> Response of the density matrix
   real(wp), intent(out) :: dp(:, :)
   !> Perturbation of the overlap in the atomic orbital basis
   real(wp), intent(in), optional :: ds(:, :)

   integer :: nact, i, q, spin
   real(wp) :: de, def
   real(wp), allocatable :: hmo(:, :), smo(:, :), wmo(:, :)
   real(wp), allocatable :: f(:, :), fp(:, :), dep(:)

   call smear_occupations(wfn%focc(:, 1), wfn%focc(:, 2), wfn%kt, f, fp)
   call transform_perturbation(wfn, f, dh, ds, nact, hmo, smo, dep)
   if (nact == 0) then
      dp(:, :) = 0.0_wp
      return
   end if

   allocate(wmo(nact, size(dh, 1)), source=0.0_wp)
   do spin = 1, 2
      def = fermi_shift(fp(1:nact, spin), dep)
      do q = 1, size(dh, 1)
         do i = 1, nact
            de = wfn%emo(i, 1) - wfn%emo(q, 1)
            if (i == q .or. abs(de) < degtol) then
               wmo(i, q) = wmo(i, q) - f(i, spin) * smo(i, q) &
                  & + fp(i, spin) * (hmo(i, q) - wfn%emo(i, 1) * smo(i, q))
            else
               wmo(i, q) = wmo(i, q) + ((f(i, spin) - f(q, spin)) * hmo(i, q) &
                  & - (f(i, spin) * wfn%emo(i, 1) - f(q, spin) * wfn%emo(q, 1)) &
                  & * smo(i, q)) / de
            end if
         end do
      end do
      do i = 1, nact
         wmo(i, i) = wmo(i, i) - def * fp(i, spin)
      end do
   end do

   call backtransform_response(wfn, nact, wmo, dp)

end subroutine density_response

subroutine wdensity_response(wfn, dh, dw, ds)
   !***********************************
   !* First-order energy-weighted density matrix, W = sum_p f_p e_p C_p C_p^T,
   !* for the perturbations dh and ds, with the conventions of
   !* density_response. With g = f e the response in the orbital basis is
   !* ((g_p - g_q) H_pq - (g_p e_p - g_q e_q) S_pq)/(e_p - e_q) and
   !* (f_p + e_p f'_p) (H_pq - e_p S_pq) - g_p S_pq for degenerate pairs.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Perturbation of the Hamiltonian in the atomic orbital basis
   real(wp), intent(in) :: dh(:, :)
   !> Response of the energy-weighted density matrix
   real(wp), intent(out) :: dw(:, :)
   !> Perturbation of the overlap in the atomic orbital basis
   real(wp), intent(in), optional :: ds(:, :)

   integer :: nact, i, q, spin
   real(wp) :: de, def, gi, gq
   real(wp), allocatable :: hmo(:, :), smo(:, :), wmo(:, :)
   real(wp), allocatable :: f(:, :), fp(:, :), dep(:)

   call smear_occupations(wfn%focc(:, 1), wfn%focc(:, 2), wfn%kt, f, fp)
   call transform_perturbation(wfn, f, dh, ds, nact, hmo, smo, dep)
   if (nact == 0) then
      dw(:, :) = 0.0_wp
      return
   end if

   allocate(wmo(nact, size(dh, 1)), source=0.0_wp)
   do spin = 1, 2
      def = fermi_shift(fp(1:nact, spin), dep)
      do q = 1, size(dh, 1)
         do i = 1, nact
            de = wfn%emo(i, 1) - wfn%emo(q, 1)
            gi = f(i, spin) * wfn%emo(i, 1)
            gq = f(q, spin) * wfn%emo(q, 1)
            if (i == q .or. abs(de) < degtol) then
               wmo(i, q) = wmo(i, q) - gi * smo(i, q) &
                  & + (f(i, spin) + wfn%emo(i, 1) * fp(i, spin)) &
                  & * (hmo(i, q) - wfn%emo(i, 1) * smo(i, q))
            else
               wmo(i, q) = wmo(i, q) + ((gi - gq) * hmo(i, q) &
                  & - (gi * wfn%emo(i, 1) - gq * wfn%emo(q, 1)) * smo(i, q)) / de
            end if
         end do
      end do
      do i = 1, nact
         wmo(i, i) = wmo(i, i) - def * wfn%emo(i, 1) * fp(i, spin)
      end do
   end do

   call backtransform_response(wfn, nact, wmo, dw)

end subroutine wdensity_response

subroutine density_response_derivative(wfn, dhs, dh, ds, dp)
   !***********************************
   !* Derivative of the density response to the perturbations dh and ds
   !* w.r.t. the converged state, which changes due to a perturbation dhs of
   !* the Hamiltonian with unchanged overlap. With Hs = C^T dhs C the
   !* orbitals rotate by U_pq = Hs_pq/(e_q - e_p) and the orbital energies
   !* shift by Hs_pp. Only valid for integer occupation numbers, rotations
   !* between degenerate orbitals are not included.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Perturbation of the Hamiltonian changing the state
   real(wp), intent(in) :: dhs(:, :)
   !> Perturbation of the Hamiltonian in the atomic orbital basis
   real(wp), intent(in) :: dh(:, :)
   !> Perturbation of the overlap in the atomic orbital basis
   real(wp), intent(in) :: ds(:, :)
   !> Derivative of the response of the density matrix
   real(wp), intent(out) :: dp(:, :)

   integer :: nao, p, q
   real(wp) :: de
   real(wp), allocatable :: hmo(:, :), smo(:, :), hsmo(:, :), tmp(:, :), f(:), fe(:)
   real(wp), allocatable :: umat(:, :), w0(:, :), dw0(:, :), hcom(:, :), scom(:, :)

   nao = size(dh, 1)
   allocate(hmo(nao, nao), smo(nao, nao), hsmo(nao, nao), tmp(nao, nao), umat(nao, nao), &
      & w0(nao, nao), dw0(nao, nao), hcom(nao, nao), scom(nao, nao))
   call gemm(wfn%coeff(:, :, 1), dh, tmp, transa="t")
   call gemm(tmp, wfn%coeff(:, :, 1), hmo)
   call gemm(wfn%coeff(:, :, 1), ds, tmp, transa="t")
   call gemm(tmp, wfn%coeff(:, :, 1), smo)
   call gemm(wfn%coeff(:, :, 1), dhs, tmp, transa="t")
   call gemm(tmp, wfn%coeff(:, :, 1), hsmo)

   !>-- half of the total occupation, the factor two is restored at the end
   f = 0.5_wp * (wfn%focc(:, 1) + wfn%focc(:, 2))
   fe = f * wfn%emo(:, 1)

   do q = 1, nao
      do p = 1, nao
         de = wfn%emo(p, 1) - wfn%emo(q, 1)
         if (p == q) then
            umat(p, q) = 0.0_wp
            w0(p, q) = -f(p) * smo(p, q)
         else if (abs(de) < degtol) then
            umat(p, q) = 0.0_wp
            w0(p, q) = 0.0_wp
         else
            umat(p, q) = -hsmo(p, q) / de
            w0(p, q) = ((f(p) - f(q)) * hmo(p, q) - (fe(p) - fe(q)) * smo(p, q)) / de
         end if
      end do
   end do

   !>-- commutators of the perturbations with the orbital rotation
   call gemm(hmo, umat, hcom)
   call gemm(umat, hmo, hcom, alpha=-1.0_wp, beta=1.0_wp)
   call gemm(smo, umat, scom)
   call gemm(umat, smo, scom, alpha=-1.0_wp, beta=1.0_wp)

   do q = 1, nao
      do p = 1, nao
         de = wfn%emo(p, 1) - wfn%emo(q, 1)
         if (p == q) then
            dw0(p, q) = -f(p) * scom(p, q)
         else if (abs(de) < degtol) then
            dw0(p, q) = 0.0_wp
         else
            dw0(p, q) = ((f(p) - f(q)) * hcom(p, q) &
               & - (f(p) * hsmo(p, p) - f(q) * hsmo(q, q)) * smo(p, q) &
               & - (fe(p) - fe(q)) * scom(p, q) &
               & - w0(p, q) * (hsmo(p, p) - hsmo(q, q))) / de
         end if
      end do
   end do

   call gemm(umat, w0, dw0, beta=1.0_wp)
   call gemm(w0, umat, dw0, alpha=-1.0_wp, beta=1.0_wp)
   call gemm(wfn%coeff(:, :, 1), dw0, tmp)
   call gemm(tmp, wfn%coeff(:, :, 1), dp, transb="t", alpha=2.0_wp)
end subroutine density_response_derivative

pure function has_integer_occupation(wfn) result(integer_occupation)
   !***********************************
   !* Whether all orbitals are either occupied or empty in both spin channels.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> All occupation numbers are zero or one
   logical :: integer_occupation

   integer_occupation = all(abs(wfn%focc(:, 1:2) * (1.0_wp - wfn%focc(:, 1:2))) < 1.0e-8_wp)
end function has_integer_occupation

subroutine transform_perturbation(wfn, f, dh, ds, nact, hmo, smo, dep)
   !***********************************
   !* Transform the perturbations to the orbital basis, restricted to rows
   !* of orbitals with nonzero occupation. Occupations decrease with the
   !* orbital energy, hence those orbitals are the first nact ones.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Occupation numbers per spin channel
   real(wp), intent(in) :: f(:, :)
   !> Perturbation of the Hamiltonian in the atomic orbital basis
   real(wp), intent(in) :: dh(:, :)
   !> Perturbation of the overlap in the atomic orbital basis
   real(wp), intent(in), optional :: ds(:, :)
   !> Number of orbitals with nonzero occupation
   integer, intent(out) :: nact
   !> Perturbation of the Hamiltonian in the orbital basis
   real(wp), allocatable, intent(out) :: hmo(:, :)
   !> Perturbation of the overlap in the orbital basis
   real(wp), allocatable, intent(out) :: smo(:, :)
   !> First-order orbital energies
   real(wp), allocatable, intent(out) :: dep(:)

   integer :: nao, i
   real(wp), allocatable :: tmp(:, :)

   nao = size(dh, 1)
   nact = count(max(f(:, 1), f(:, 2)) > occtol)
   if (nact == 0) return

   allocate(tmp(nact, nao), hmo(nact, nao), dep(nact))
   allocate(smo(nact, nao), source=0.0_wp)
   call gemm(wfn%coeff(:, 1:nact, 1), dh, tmp, transa="t")
   call gemm(tmp, wfn%coeff(:, :, 1), hmo)
   if (present(ds)) then
      call gemm(wfn%coeff(:, 1:nact, 1), ds, tmp, transa="t")
      call gemm(tmp, wfn%coeff(:, :, 1), smo)
   end if

   do i = 1, nact
      dep(i) = hmo(i, i) - wfn%emo(i, 1) * smo(i, i)
   end do
end subroutine transform_perturbation

subroutine backtransform_response(wfn, nact, wmo, dp)
   !***********************************
   !* Transform a response given for the occupied rows of the orbital basis
   !* to the atomic orbital basis, dp = C_act wmo C^T + transpose, less the
   !* double counted block of the occupied orbitals.
   !***********************************
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Number of orbitals with nonzero occupation
   integer, intent(in) :: nact
   !> Response in the orbital basis
   real(wp), intent(in) :: wmo(:, :)
   !> Response in the atomic orbital basis
   real(wp), intent(out) :: dp(:, :)

   real(wp), allocatable :: z(:, :), w2(:, :)

   allocate(z(nact, size(dp, 1)), w2(nact, size(dp, 1)))
   call gemm(wmo, wfn%coeff(:, :, 1), z, transb="t")
   call gemm(wmo(:, 1:nact), wfn%coeff(:, 1:nact, 1), w2, transb="t")
   call gemm(wfn%coeff(:, 1:nact, 1), z, dp)
   call gemm(z, wfn%coeff(:, 1:nact, 1), dp, transa="t", transb="t", beta=1.0_wp)
   call gemm(wfn%coeff(:, 1:nact, 1), w2, dp, alpha=-1.0_wp, beta=1.0_wp)
end subroutine backtransform_response

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
