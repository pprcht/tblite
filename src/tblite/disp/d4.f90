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

!> @file tblite/disp/d4.f90
!> Provides a proxy for the [DFT-D4 dispersion correction](https://dftd4.readthedocs.io)

!> Generally applicable charge-dependent London-dispersion correction, DFT-D4.
module tblite_disp_d4
   use dftd4, only : dispersion_model, d4_model, d4s_model, &
      & rational_damping_param, realspace_cutoff, &
      & new_d4_model, new_d4s_model
   use dftd4_cutoff, only : smooth_cutoff, smooth_cutoff_r2
   use dftd4_disp, only : get_dispersion3_hessian
   use dftd4_model, only : d4_qmod
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : new_ncoord, ncoord_type, cn_count
   use tblite_blas, only : gemv
   use tblite_container_cache, only : container_cache
   use tblite_cutoff, only : get_lattice_points
   use tblite_disp_cache, only : dispersion_cache
   use tblite_disp_type, only : dispersion_type
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian, add_ncoord_hessian
   use tblite_partition, only : work_partition, owns_pair
   use tblite_response_type, only : kernel_type, response_type
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   public :: new_d4_dispersion, new_d4s_dispersion


   !> Container for self-consistent D4 dispersion interactions
   type, public, extends(dispersion_type) :: d4_dispersion
      !> Instance of the actual D4 dispersion model
      class(dispersion_model), allocatable :: model
      !> Rational damping parameters
      type(rational_damping_param) :: param
      !> Selected real space cutoffs for this instance
      type(realspace_cutoff) :: cutoff
      !> Coordination number instance
      class(ncoord_type), allocatable :: ncoord
   contains
      !> Update dispersion cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the dispersion correction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the dispersion correction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the dispersion correction
      procedure :: get_potential
      !> Evaluate gradient contributions from the selfconsistent dispersion correction
      procedure :: get_gradient
      !> Evaluate second derivative w.r.t. the charges
      procedure :: get_response_kernel
      !> Evaluate mixed second derivative w.r.t. the charges and positions
      procedure :: get_response_derivative
      !> Evaluate second derivative w.r.t. the positions at fixed charges
      procedure :: get_hessian
   end type d4_dispersion

   character(len=*), parameter :: label_d4 = "self-consistent DFT-D4 dispersion"
   character(len=*), parameter :: label_d4s = "self-consistent DFT-D4S dispersion"
   real(wp), parameter :: default_disp2_width = 0.0_wp
   real(wp), parameter :: default_disp3_width = 0.0_wp


contains


!> Create a new instance of a self-consistent D4 dispersion correction
subroutine new_d4_dispersion(self, mol, s6, s8, a1, a2, s9, error, disp2_width, disp3_width)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   real(wp), intent(in) :: s6, s8, a1, a2, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Width of smooth two-body interaction cutoff
   real(wp), intent(in), optional :: disp2_width
   !> Width of smooth three-body interaction cutoff
   real(wp), intent(in), optional :: disp3_width

   type(d4_model), allocatable :: tmp
   real(wp) :: width2, width3

   width2 = default_disp2_width
   width3 = default_disp3_width
   if (present(disp2_width)) width2 = disp2_width
   if (present(disp3_width)) width3 = disp3_width

   self%label = label_d4

   ! Create a new instance of the D4 model
   allocate(tmp)
   call new_d4_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)
   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp, &
      & width2=width2, width3=width3)

   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4_dispersion


!> Create a new instance of a self-consistent D4S dispersion correction
subroutine new_d4s_dispersion(self, mol, s6, s8, a1, a2, s9, error, disp2_width, disp3_width)
   !> Instance of the dispersion correction
   type(d4_dispersion), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   real(wp), intent(in) :: s6, s8, a1, a2, s9
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Width of smooth two-body interaction cutoff
   real(wp), intent(in), optional :: disp2_width
   !> Width of smooth three-body interaction cutoff
   real(wp), intent(in), optional :: disp3_width

   type(d4s_model), allocatable :: tmp
   real(wp) :: width2, width3

   width2 = default_disp2_width
   width3 = default_disp3_width
   if (present(disp2_width)) width2 = disp2_width
   if (present(disp3_width)) width3 = disp3_width

   self%label = label_d4s

   ! Create a new instance of the D4S model
   allocate(tmp)
   call new_d4s_model(error, tmp, mol, qmod=d4_qmod%gfn2)
   if(allocated(error)) return
   call move_alloc(tmp, self%model)

   self%param = rational_damping_param(s6=s6, s8=s8, s9=s9, a1=a1, a2=a2)
   self%cutoff = realspace_cutoff(disp3=25.0_wp, disp2=50.0_wp, &
      & width2=width2, width3=width3)

   call new_ncoord(self%ncoord, mol, cn_count%dftd4, error, &
      & cutoff=self%cutoff%cn, rcov=self%model%rcov, en=self%model%en)
end subroutine new_d4s_dispersion


!> Update dispersion cache
subroutine update(self, mol, cache)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache

   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr
   integer :: mref

   call taint(cache, ptr)
   mref = maxval(self%model%ref)

   if (.not.allocated(ptr%cn)) allocate(ptr%cn(mol%nat))
   if (.not.allocated(ptr%dcndr)) allocate(ptr%dcndr(3, mol%nat, mol%nat))
   if (.not.allocated(ptr%dcndL)) allocate(ptr%dcndL(3, 3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
   call self%ncoord%get_coordination_number(mol, lattr, ptr%cn, ptr%dcndr, ptr%dcndL)

   if (.not.allocated(ptr%dispmat)) allocate(ptr%dispmat(mref, mol%nat, mref, mol%nat))
   if (.not.allocated(ptr%vvec)) allocate(ptr%vvec(mref, mol%nat, self%model%ncoup))
   if (.not.allocated(ptr%gwvec)) allocate(ptr%gwvec(mref, mol%nat, self%model%ncoup))
   if (.not.allocated(ptr%dgwdq)) allocate(ptr%dgwdq(mref, mol%nat, self%model%ncoup))

   ! Keep both symmetric blocks of each owned pair. Energy and potential then
   ! follow the same partition without reducing the matrix itself.
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call get_dispersion_matrix(mol, self%model, self%param, lattr, self%cutoff%disp2, &
      & self%cutoff%width2, self%model%r4r2, ptr%dispmat, self%partition)
end subroutine update


!> Evaluate non-selfconsistent part of the dispersion correction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)
   !> Dispersion gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)

   type(dispersion_cache), pointer :: ptr

   logical :: grad
   integer :: mref
   real(wp), allocatable :: qat(:)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:)
   real(wp), allocatable :: lattr(:, :)

   call view(cache, ptr)

   mref = maxval(self%model%ref)
   grad = present(gradient).and.present(sigma)

   allocate(gwvec(mref, mol%nat, self%model%ncoup), qat(mol%nat), c6(mol%nat, mol%nat))
   if (grad) then
      allocate(gwdcn(mref, mol%nat, self%model%ncoup), gwdq(mref, mol%nat, self%model%ncoup), &
         & dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   end if
   qat(:) = 0.0_wp
   if (grad) then
      allocate(dEdcn(mol%nat), dEdq(mol%nat))
      dEdcn(:) = 0.0_wp
      dEdq(:) = 0.0_wp
   end if

   call self%model%weight_references(mol, ptr%cn, qat, gwvec, gwdcn, gwdq)
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp3, lattr)
   call self%param%get_dispersion3(mol, lattr, self%cutoff%disp3, self%cutoff%width3, &
      & self%model%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma, &
      & partition=self%partition%get_d4())
   if (grad) then
      call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
      call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
   end if

end subroutine get_engrad


!> Evaluate selfconsistent energy of the dispersion correction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion energy
   real(wp), intent(inout) :: energies(:)

   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), ptr%gwvec)

   if (self%model%ncoup > 1) then
      call contract_dispersion(self, mol, ptr%dispmat, ptr%gwvec, ptr%gwvec, &
         & energies, 0.5_wp)
   else
      ! Dispersion energy with atom-wise weighting
      call gemv(ptr%dispmat, ptr%gwvec(:, :, 1), ptr%vvec(:, :, 1), alpha=0.5_wp)
      ptr%vvec(:, :, 1) = ptr%vvec(:, :, 1) * ptr%gwvec(:, :, 1)
      energies(:) = energies + sum(ptr%vvec(:, :, 1), 1)
   end if

end subroutine get_energy


!> Evaluate charge dependent potential shift from the dispersion correction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot

   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)

   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), ptr%gwvec, ptr%vvec, &
      & ptr%dgwdq)

   if (self%model%ncoup > 1) then
      call contract_dispersion(self, mol, ptr%dispmat, ptr%dgwdq, ptr%gwvec, &
         & pot%vat(:, 1), 1.0_wp)
   else
      ! Dispersion energy with atom-wise weighting
      call gemv(ptr%dispmat, ptr%gwvec(:, :, 1), ptr%vvec(:, :, 1))
      ptr%vvec(:, :, 1) = ptr%vvec(:, :, 1) * ptr%dgwdq(:, :, 1)
      pot%vat(:, 1) = pot%vat(:, 1) + sum(ptr%vvec(:, :, 1), 1)
   end if

end subroutine get_potential


!> Contract the pairwise D4S weights for energy or charge potential.
subroutine contract_dispersion(self, mol, dispmat, left, right, values, scale)
   class(d4_dispersion), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: dispmat(:, :, :, :), left(:, :, :), right(:, :, :)
   real(wp), intent(inout) :: values(:)
   real(wp), intent(in) :: scale

   integer :: iat, jat, izp, jzp, iref, jref
   real(wp) :: val

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(self, mol, dispmat, left, right, values, scale) &
   !$omp private(iat, jat, izp, jzp, iref, jref, val)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      val = 0.0_wp
      do jat = 1, mol%nat
         if (.not.owns_pair(self%partition, max(iat, jat), min(iat, jat))) cycle
         jzp = mol%id(jat)
         do iref = 1, self%model%ref(izp)
            do jref = 1, self%model%ref(jzp)
               val = val + scale * dispmat(iref, iat, jref, jat) * &
                  & left(iref, iat, jat) * right(jref, jat, iat)
            end do
         end do
      end do
      values(iat) = values(iat) + val
   end do
end subroutine contract_dispersion


!> Evaluate gradient contributions from the selfconsistent dispersion correction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Dispersion gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Dispersion virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   integer :: mref
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), dc6dq(:, :)
   real(wp), allocatable :: dEdcn(:), dEdq(:), energies(:)
   real(wp), allocatable :: lattr(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   mref = maxval(self%model%ref)

   allocate(gwvec(mref, mol%nat, self%model%ncoup), gwdcn(mref, mol%nat, self%model%ncoup), &
      &  gwdq(mref, mol%nat, self%model%ncoup))
   call self%model%weight_references(mol, ptr%cn, wfn%qat(:, 1), gwvec, gwdcn, gwdq)

   allocate(c6(mol%nat, mol%nat), dc6dcn(mol%nat, mol%nat), dc6dq(mol%nat, mol%nat))
   call self%model%get_atomic_c6(mol, gwvec, gwdcn, gwdq, c6, dc6dcn, dc6dq)

   allocate(energies(mol%nat), dEdcn(mol%nat), dEdq(mol%nat))
   energies(:) = 0.0_wp
   dEdcn(:) = 0.0_wp
   dEdq(:) = 0.0_wp
   call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
   call self%param%get_dispersion2(mol, lattr, self%cutoff%disp2, self%cutoff%width2, &
      & self%model%r4r2, c6, dc6dcn, dc6dq, energies, dEdcn, dEdq, gradient, sigma, &
      & partition=self%partition%get_d4())
   call gemv(ptr%dcndr, dEdcn, gradient, beta=1.0_wp)
   call gemv(ptr%dcndL, dEdcn, sigma, beta=1.0_wp)
end subroutine get_gradient


subroutine get_dispersion_matrix(mol, disp, param, trans, cutoff, width, r4r2, dispmat, partition)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Damping parameters
   type(rational_damping_param), intent(in) :: param
   !> Instance of the dispersion model
   class(dispersion_model), intent(in) :: disp
   !> Lattice points
   real(wp), intent(in) :: trans(:, :)
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Width of smooth cutoff
   real(wp), intent(in) :: width
   !> Expectation values for r4 over r2 operator
   real(wp), intent(in) :: r4r2(:)
   !> Dispersion matrix
   real(wp), intent(out) :: dispmat(:, :, :, :)
   !> Work partition of the atom pairs
   type(work_partition), intent(in) :: partition

   integer :: iat, jat, izp, jzp, jtr, iref, jref
   real(wp) :: vec(3), r2, r, cutoff2, r0ij, rrij, t6, t8
   real(wp) :: edisp, dE, sw, dswdr

   dispmat(:, :, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do schedule(runtime) default(none) &
   !$omp shared(mol, param, disp, trans, cutoff, width, cutoff2, r4r2, dispmat, partition) &
   !$omp private(iat, jat, izp, jzp, jtr, vec, r2, r0ij, rrij, &
   !$omp& t6, t8, edisp, dE, r, sw, dswdr)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         if (.not.owns_pair(partition, iat, jat)) cycle
         jzp = mol%id(jat)
         rrij = 3*r4r2(izp)*r4r2(jzp)
         r0ij = param%a1 * sqrt(rrij) + param%a2
         dE = 0.0_wp
         do jtr = 1, size(trans, 2)
            vec(:) = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, jtr))
            r2 = vec(1)*vec(1) + vec(2)*vec(2) + vec(3)*vec(3)
            if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
            r = sqrt(r2)
            call smooth_cutoff(r, cutoff, width, sw, dswdr)
            if (sw <= 0.0_wp) cycle

            t6 = 1.0_wp/(r2**3 + r0ij**6)
            t8 = 1.0_wp/(r2**4 + r0ij**8)

            edisp = sw * (param%s6*t6 + param%s8*rrij*t8)

            dE = dE - edisp
         end do

         do iref = 1, disp%ref(izp)
            do jref = 1, disp%ref(jzp)
               dispmat(iref, iat, jref, jat) = dE * disp%c6(iref, jref, izp, jzp)
               dispmat(jref, jat, iref, iat) = dE * disp%c6(jref, iref, jzp, izp)
            end do
         end do
      end do
   end do
end subroutine get_dispersion_matrix


subroutine get_response_kernel(self, mol, cache, wfn, kernel, error, dwfn)
   !***********************************
   !* Add the second charge derivative of the self-consistent D4 energy to
   !* the atomic kernel. The energy is bilinear in the reference weights
   !* gw = w(CN) zeta(q) with the damped reference C6 matrix of the cache,
   !* hence the kernel of a pair consists of the second derivative of the
   !* weights of one atom and the first derivatives of the weights of both.
   !* With dwfn, the third charge derivative contracted with the first-order
   !* charges is added instead.
   !***********************************
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Response kernel
   type(kernel_type), intent(inout) :: kernel
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order change of the state
   type(wavefunction_type), intent(in), optional :: dwfn

   integer :: iat, jat, iref, jref
   real(wp) :: ha, hb, w6
   real(wp), allocatable :: gw(:, :, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   call get_charge_weights(self, mol, ptr, wfn%qat(:, 1), gw, error)
   if (allocated(error)) return

   do iat = 1, mol%nat
      do jat = 1, mol%nat
         ha = 0.0_wp
         hb = 0.0_wp
         do iref = 1, size(gw, 1)
            do jref = 1, size(gw, 1)
               w6 = ptr%dispmat(iref, iat, jref, jat)
               if (present(dwfn)) then
                  ha = ha + w6 * (gw(iref, iat, 3)*gw(jref, jat, 0)*dwfn%qat(iat, 1) &
                     & + gw(iref, iat, 2)*gw(jref, jat, 1)*dwfn%qat(jat, 1))
                  hb = hb + w6 * (gw(iref, iat, 2)*gw(jref, jat, 1)*dwfn%qat(iat, 1) &
                     & + gw(iref, iat, 1)*gw(jref, jat, 2)*dwfn%qat(jat, 1))
               else
                  ha = ha + gw(iref, iat, 2) * gw(jref, jat, 0) * w6
                  hb = hb + gw(iref, iat, 1) * gw(jref, jat, 1) * w6
               end if
            end do
         end do
         kernel%kat(iat, iat) = kernel%kat(iat, iat) + ha
         kernel%kat(iat, jat) = kernel%kat(iat, jat) + hb
      end do
   end do
end subroutine get_response_kernel


subroutine get_response_derivative(self, mol, cache, wfn, response, error, dwfn)
   !***********************************
   !* Mixed second derivative of the self-consistent D4 energy w.r.t. the
   !* charges and the atomic positions. The atomic potential is
   !* v_i = -sum_j dC6_ij/dq_i F(r_ij) and depends on the positions via the
   !* damped radial function and via the coordination numbers in the
   !* reference weights. The latter part is collected per pair of atoms and
   !* contracted with the derivative of the coordination numbers. With dwfn,
   !* the charge derivative of v contracted with the first-order charges is
   !* differentiated instead. Interactions of an atom with its own periodic
   !* images are not included.
   !***********************************
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> First-order quantities of the perturbations
   type(response_type), intent(inout) :: response
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order change of the state
   type(wavefunction_type), intent(in), optional :: dwfn

   integer :: nat, np, mref, iat, jat, izp, jzp, itr, iref, jref, ip, kat, kc
   real(wp) :: vec(3), r2, cutoff2, fsw, df, dfvec(3), refc6
   real(wp) :: dc6dq, d2c6dqdcni, d2c6dqdcnj, c, u, g, qi, qj
   real(wp), allocatable :: lattr(:, :), qat(:)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), allocatable :: gwdcndq(:, :), gwd2q(:, :), gwdcnd2q(:, :), dvdcn(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   nat = mol%nat
   np = size(response%coord)
   qat = wfn%qat(:, 1)

   select type(model => self%model)
   type is (d4_model)
      mref = maxval(model%ref)
      allocate(gwvec(mref, nat, model%ncoup), gwdcn(mref, nat, model%ncoup), &
         & gwdq(mref, nat, model%ncoup))
      allocate(gwdcndq(mref, nat), gwd2q(mref, nat), gwdcnd2q(mref, nat), &
         & dvdcn(nat, nat), source=0.0_wp)
      call model%weight_references(mol, ptr%cn, qat, gwvec, gwdcn, gwdq)

      !>-- the weights factorize, gw = w(CN) zeta(q). With the logarithmic
      !>   derivative g = zeta'/zeta the mixed derivative is g dw/dCN and
      !>   the second charge derivative follows from g' = g (c*qref/u**2 - 2/u)
      do iat = 1, nat
         izp = mol%id(iat)
         c = model%eta(izp) * model%gc
         u = qat(iat) + model%zeff(izp)
         do iref = 1, model%ref(izp)
            if (abs(gwvec(iref, iat, 1)) < tiny(1.0_wp)) cycle
            g = gwdq(iref, iat, 1) / gwvec(iref, iat, 1)
            gwdcndq(iref, iat) = gwdcn(iref, iat, 1) * g
            if (abs(u) < 1.0e-12_wp) cycle
            g = g * (g + c*(model%q(iref, izp) + model%zeff(izp))/u**2 - 2.0_wp/u)
            gwd2q(iref, iat) = gwvec(iref, iat, 1) * g
            gwdcnd2q(iref, iat) = gwdcn(iref, iat, 1) * g
         end do
      end do

      cutoff2 = self%cutoff%disp2**2
      call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
      do iat = 1, nat
         izp = mol%id(iat)
         do jat = 1, nat
            if (iat == jat) cycle
            jzp = mol%id(jat)
            dc6dq = 0.0_wp
            d2c6dqdcni = 0.0_wp
            d2c6dqdcnj = 0.0_wp
            do iref = 1, model%ref(izp)
               do jref = 1, model%ref(jzp)
                  refc6 = model%c6(iref, jref, izp, jzp)
                  if (present(dwfn)) then
                     qi = dwfn%qat(iat, 1)
                     qj = dwfn%qat(jat, 1)
                     dc6dq = dc6dq + refc6 &
                        & * (gwd2q(iref, iat)*gwvec(jref, jat, 1)*qi &
                        & + gwdq(iref, iat, 1)*gwdq(jref, jat, 1)*qj)
                     d2c6dqdcni = d2c6dqdcni + refc6 &
                        & * (gwdcnd2q(iref, iat)*gwvec(jref, jat, 1)*qi &
                        & + gwdcndq(iref, iat)*gwdq(jref, jat, 1)*qj)
                     d2c6dqdcnj = d2c6dqdcnj + refc6 &
                        & * (gwd2q(iref, iat)*gwdcn(jref, jat, 1)*qi &
                        & + gwdq(iref, iat, 1)*gwdcndq(jref, jat)*qj)
                  else
                     dc6dq = dc6dq + gwdq(iref, iat, 1)*gwvec(jref, jat, 1)*refc6
                     d2c6dqdcni = d2c6dqdcni + gwdcndq(iref, iat)*gwvec(jref, jat, 1)*refc6
                     d2c6dqdcnj = d2c6dqdcnj + gwdq(iref, iat, 1)*gwdcn(jref, jat, 1)*refc6
                  end if
               end do
            end do
            do itr = 1, size(lattr, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - lattr(:, itr)
               r2 = sum(vec**2)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
               call get_damped_radial(self, model, izp, jzp, r2, fsw, df)
               dfvec(:) = dc6dq*df*vec/sqrt(r2)
               do ip = 1, np
                  response%hessian(:, iat, ip) = response%hessian(:, iat, ip) &
                     & - response%qat(iat, ip)*dfvec
                  response%hessian(:, jat, ip) = response%hessian(:, jat, ip) &
                     & + response%qat(iat, ip)*dfvec
                  kat = (response%coord(ip) - 1)/3 + 1
                  kc = response%coord(ip) - 3*(kat - 1)
                  if (kat == iat) response%vat(iat, ip) = response%vat(iat, ip) - dfvec(kc)
                  if (kat == jat) response%vat(iat, ip) = response%vat(iat, ip) + dfvec(kc)
               end do
               dvdcn(iat, iat) = dvdcn(iat, iat) - fsw*d2c6dqdcni
               dvdcn(iat, jat) = dvdcn(iat, jat) - fsw*d2c6dqdcnj
            end do
         end do
      end do

      !>-- dependence of the reference weights on the coordination numbers
      do ip = 1, np
         call gemv(ptr%dcndr, matmul(response%qat(:, ip), dvdcn), &
            & response%hessian(:, :, ip), beta=1.0_wp)
         kat = (response%coord(ip) - 1)/3 + 1
         kc = response%coord(ip) - 3*(kat - 1)
         response%vat(:, ip) = response%vat(:, ip) + matmul(dvdcn, ptr%dcndr(kc, kat, :))
      end do
   class default
      call fatal_error(error, "Analytic second derivatives are only available for "// &
         & "the D4 model")
   end select
end subroutine get_response_derivative


subroutine get_hessian(self, mol, cache, wfn, hessian, error)
   !***********************************
   !* Add the Hessian of the D4 dispersion energy at fixed charges. The
   !* three-body term is evaluated for neutral atoms by dftd4. For the
   !* two-body term, E = -sum_ij C6_ij(CN_i, CN_j) F(r_ij), the direct second
   !* derivative of the damped radial function F is added per pair, the
   !* dependence on the coordination numbers is collected and contracted
   !* with their derivatives. Interactions of an atom with its own periodic
   !* images are not included.
   !***********************************
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different dispersion runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), contiguous, intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, mref, iat, jat, izp, jzp, itr, ii, jj
   real(wp) :: vec(3), r1, r2, cutoff2, rrij, r0ij, t6, t8, f, dfdr2, d2fdr22
   real(wp) :: sw, dswdr2, d2swdr22, fsw, df, d2f, block(3, 3), dfvec(3)
   real(wp), allocatable :: lattr(:, :), tmp(:, :)
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwd2cn(:, :, :)
   real(wp), allocatable :: c6(:, :), dc6dcn(:, :), d2c6dcn2(:, :), d2c6dcnij(:, :)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :), dEdcndcn(:, :)
   type(dispersion_cache), pointer :: ptr

   call view(cache, ptr)
   nat = mol%nat

   select type(model => self%model)
   type is (d4_model)
      allocate(tmp(3*nat, 3*nat))
      call get_dispersion3_hessian(error, mol, model, self%param, self%cutoff, tmp)
      if (allocated(error)) return
      hessian(:, :) = hessian + tmp

      mref = maxval(model%ref)
      allocate(gwvec(mref, nat, model%ncoup), gwdcn(mref, nat, model%ncoup), &
         & gwd2cn(mref, nat, model%ncoup))
      allocate(c6(nat, nat), dc6dcn(nat, nat), d2c6dcn2(nat, nat), d2c6dcnij(nat, nat))
      call model%weight_references_hessian(mol, ptr%cn, wfn%qat(:, 1), gwvec, gwdcn, &
         & gwd2cn)
      call model%get_atomic_c6_hessian(mol, gwvec, gwdcn, gwd2cn, c6, dc6dcn, &
         & d2c6dcn2, d2c6dcnij)
      allocate(dEdcn(nat), dEdcndr(3*nat, nat), dEdcndcn(nat, nat), source=0.0_wp)

      cutoff2 = self%cutoff%disp2**2
      call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%disp2, lattr)
      do iat = 1, nat
         izp = mol%id(iat)
         ii = 3*(iat - 1)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jj = 3*(jat - 1)
            rrij = 3*model%r4r2(izp)*model%r4r2(jzp)
            r0ij = self%param%a1*sqrt(rrij) + self%param%a2
            do itr = 1, size(lattr, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - lattr(:, itr)
               r2 = sum(vec**2)
               if (r2 > cutoff2 .or. r2 < epsilon(1.0_wp)) cycle
               r1 = sqrt(r2)

               !>-- damped radial function and its derivatives w.r.t. r**2
               t6 = 1.0_wp/(r2**3 + r0ij**6)
               t8 = 1.0_wp/(r2**4 + r0ij**8)
               f = self%param%s6*t6 + self%param%s8*rrij*t8
               dfdr2 = -3*self%param%s6*r2**2*t6**2 - 4*self%param%s8*rrij*r2**3*t8**2
               d2fdr22 = self%param%s6*(18*r2**4*t6**3 - 6*r2*t6**2) &
                  & + self%param%s8*rrij*(32*r2**6*t8**3 - 12*r2**2*t8**2)
               call smooth_cutoff_r2(r2, self%cutoff%disp2, self%cutoff%width2, sw, &
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

      call get_lattice_points(mol%periodic, mol%lattice, self%cutoff%cn, lattr)
      call self%ncoord%add_coordination_number_hessian(mol, lattr, dEdcn, hessian)
      call add_ncoord_hessian(hessian, ptr%dcndr, dEdcndr, dEdcndcn)
   class default
      call fatal_error(error, "Analytic second derivatives are only available for "// &
         & "the D4 model")
   end select
end subroutine get_hessian


subroutine get_charge_weights(self, mol, ptr, qat, gw, error)
   !***********************************
   !* Reference weights of the D4 model and their first three derivatives
   !* w.r.t. the atomic charge, gw(:, :, 0:3). The weights are gw = w zeta(q)
   !* with the logarithmic derivative g = zeta'/zeta, for which g' = g s and
   !* g'' = g (s**2 + s') hold with s = c*qref/u**2 - 2/u, u = q + zeff and
   !* c = eta*gc.
   !***********************************
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Dispersion cache
   type(dispersion_cache), intent(in) :: ptr
   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)
   !> Reference weights and their charge derivatives
   real(wp), allocatable, intent(out) :: gw(:, :, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: mref, iat, izp, iref
   real(wp) :: c, u, qref, g, s, ds
   real(wp), allocatable :: gwvec(:, :, :), gwdcn(:, :, :), gwdq(:, :, :)
   real(wp), parameter :: small = 1.0e-12_wp

   select type(model => self%model)
   type is (d4_model)
      mref = maxval(model%ref)
      allocate(gwvec(mref, mol%nat, model%ncoup), gwdcn(mref, mol%nat, model%ncoup), &
         & gwdq(mref, mol%nat, model%ncoup))
      allocate(gw(mref, mol%nat, 0:3), source=0.0_wp)
      call model%weight_references(mol, ptr%cn, qat, gwvec, gwdcn, gwdq)

      gw(:, :, 0) = gwvec(:, :, 1)
      gw(:, :, 1) = gwdq(:, :, 1)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         c = model%eta(izp) * model%gc
         u = qat(iat) + model%zeff(izp)
         if (abs(u) < small) cycle
         do iref = 1, model%ref(izp)
            if (abs(gw(iref, iat, 0)) < tiny(1.0_wp)) cycle
            qref = model%q(iref, izp) + model%zeff(izp)
            g = gw(iref, iat, 1) / gw(iref, iat, 0)
            s = c*qref/u**2 - 2.0_wp/u
            ds = -2*c*qref/u**3 + 2.0_wp/u**2
            gw(iref, iat, 2) = gw(iref, iat, 0) * g * (g + s)
            gw(iref, iat, 3) = gw(iref, iat, 0) * g * (g*g + 3*g*s + s*s + ds)
         end do
      end do
   class default
      call fatal_error(error, "Analytic second derivatives are only available for "// &
         & "the D4 model")
   end select
end subroutine get_charge_weights


pure subroutine get_damped_radial(self, model, izp, jzp, r2, fsw, df)
   !***********************************
   !* Damped radial function of the two-body dispersion energy including
   !* the smooth cutoff and its derivative w.r.t. the distance.
   !***********************************
   !> Instance of the dispersion correction
   class(d4_dispersion), intent(in) :: self
   !> D4 dispersion model
   type(d4_model), intent(in) :: model
   !> Species of the two atoms
   integer, intent(in) :: izp, jzp
   !> Squared distance
   real(wp), intent(in) :: r2
   !> Radial function
   real(wp), intent(out) :: fsw
   !> Derivative of the radial function w.r.t. the distance
   real(wp), intent(out) :: df

   real(wp) :: rrij, r0ij, t6, t8, f, dfdr2, sw, dswdr2, d2swdr22

   rrij = 3*model%r4r2(izp)*model%r4r2(jzp)
   r0ij = self%param%a1*sqrt(rrij) + self%param%a2
   t6 = 1.0_wp/(r2**3 + r0ij**6)
   t8 = 1.0_wp/(r2**4 + r0ij**8)
   f = self%param%s6*t6 + self%param%s8*rrij*t8
   dfdr2 = -3*self%param%s6*r2**2*t6**2 - 4*self%param%s8*rrij*r2**3*t8**2
   call smooth_cutoff_r2(r2, self%cutoff%disp2, self%cutoff%width2, sw, dswdr2, d2swdr22)
   fsw = sw*f
   df = 2*sqrt(r2)*(dswdr2*f + sw*dfdr2)
end subroutine get_damped_radial


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   use tblite_scf_info, only : scf_info, atom_resolved
   !> Instance of the electrostatic container
   class(d4_dispersion), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info(charge=atom_resolved)
end function variable_info


!> Inspect container cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(dispersion_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not.allocated(cache%raw)) then
      block
         type(dispersion_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to container cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the container cache
   type(container_cache), target, intent(inout) :: cache
   !> Reference to the container cache
   type(dispersion_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(dispersion_cache)
      ptr => target
   end select
end subroutine view

end module tblite_disp_d4
