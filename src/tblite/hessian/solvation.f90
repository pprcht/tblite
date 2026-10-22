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

!> @file tblite/hessian/solvation.f90
!> Provides the implicit solvation models entering the analytic second derivatives

!> Implicit solvation in the analytic second derivatives.
!>
!> The supported models give an energy which is exactly quadratic in the
!> atomic charges q, acting on Q = q + c with the CM5 corrections c(R),
!>
!>   E = 1/2 Q^T J(R) Q + sum_A [gamma_A + h_A Q_A**2] S_A(R),
!>
!> with the generalized Born matrix J including the ALPB shape correction and
!> the surface areas S. The response kernel K = J + diag(2 h S) therefore does
!> not depend on the charges and the solvation energy has no third derivative.
module tblite_hessian_solvation
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_blas, only : gemm, symv
   use tblite_container, only : container_type, container_list
   use tblite_hessian_solvation_born, only : add_born_radii_hessian
   use tblite_hessian_solvation_cm5, only : add_cm5_hessian
   use tblite_hessian_solvation_coupling, only : add_born_coupling_hessian
   use tblite_hessian_solvation_kernel, only : get_kernel_derivs, get_radii_combination, &
      & add_born_kernel_hessian
   use tblite_hessian_solvation_shape, only : add_shape_hessian
   use tblite_hessian_solvation_surface, only : add_surface_hessian
   use tblite_partition, only : serial_work_partition, operator(==)
   use tblite_solvation, only : alpb_solvation, cds_solvation, shift_solvation, born_kernel
   use tblite_solvation_alpb, only : add_born_mat_p16, add_born_mat_still, get_adet, &
      & get_adet_deriv
   use tblite_solvation_cm5, only : get_cm5_charges
   use tblite_xtb_calculator, only : xtb_calculator
   implicit none
   private

   public :: check_solvation_support
   public :: solvation_kernel, new_solvation_kernel
   public :: get_solvation_potential, get_solvation_potential_derivative
   public :: get_solvation_kernel_derivative
   public :: get_solvation_hessian

   !> Copy of an interaction container
   type :: container_leaf
      class(container_type), allocatable :: raw
   end type container_leaf

   !> Geometry dependent quantities of the solvation models for one structure
   type :: solvation_kernel
      !> Generalized Born model, if present
      type(alpb_solvation), allocatable :: born
      !> Surface area model, if present
      type(cds_solvation), allocatable :: cds
      !> The models act on the charges shifted by the CM5 corrections
      logical :: use_cm5 = .false.
      !> Born radii and their derivative w.r.t. the atomic positions
      real(wp), allocatable :: rad(:), draddr(:, :, :)
      !> Shape correction keps alpbet / A_det, added to every element of J
      real(wp) :: kshape = 0.0_wp
      !> Derivative of the shape correction w.r.t. the atomic positions
      real(wp), allocatable :: dkshape(:, :)
      !> Surface areas and their derivative w.r.t. the atomic positions
      real(wp), allocatable :: surface(:), dsdr(:, :, :)
      !> Surface tension and hydrogen bonding strength of each atom
      real(wp), allocatable :: tension(:), hbond(:)
      !> CM5 corrections and their derivative w.r.t. the atomic positions
      real(wp), allocatable :: cm5(:), dcm5dr(:, :, :)
      !> Response kernel, second derivative of the energy w.r.t. the charges
      real(wp), allocatable :: kmat(:, :)
   end type solvation_kernel

contains

subroutine check_solvation_support(calc, error)
   !***********************************
   !* Reject interaction containers without analytic second derivatives.
   !* Supported are the generalized Born models (ALPB and GBSA with the P16
   !* or Still kernel), the CDS surface term and the free energy shift, all
   !* evaluated over the complete work of a serial calculation.
   !***********************************
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, nborn, ncds
   logical :: has_cm5, use_cm5
   type(container_leaf), allocatable :: leaves(:)

   if (.not.allocated(calc%interactions)) return
   if (.not.(calc%interactions%partition == serial_work_partition)) then
      call fatal_error(error, "Analytic second derivatives with solvation require "// &
         & "a serial calculation")
      return
   end if

   !>-- the Born and surface terms act on the same CM5-shifted charges
   has_cm5 = .false.
   use_cm5 = .false.
   nborn = 0
   ncds = 0
   call get_leaf_containers(calc%interactions, leaves)
   do ic = 1, size(leaves)
      associate(cont => leaves(ic)%raw)
      if (.not.(cont%partition == serial_work_partition)) then
         call fatal_error(error, "Analytic second derivatives with solvation require "// &
            & "a serial calculation")
         return
      end if
      select type(cont)
      type is(alpb_solvation)
         if (cont%kernel /= born_kernel%p16 .and. cont%kernel /= born_kernel%still) then
            call fatal_error(error, "Analytic second derivatives are not available "// &
               & "for this generalized Born kernel")
            return
         end if
         call check_cm5(cont%useCM5)
         nborn = nborn + 1
      type is(cds_solvation)
         call check_cm5(cont%useCM5)
         ncds = ncds + 1
      type is(shift_solvation)
         continue
      class default
         call fatal_error(error, "Analytic second derivatives are not available for "// &
            & "the interaction '"//container_label(cont)//"'")
         return
      end select
      end associate
      if (allocated(error)) return
   end do
   if (nborn > 1 .or. ncds > 1) then
      call fatal_error(error, "Analytic second derivatives support only one generalized "// &
         & "Born and one surface area model")
   end if

contains

   subroutine check_cm5(cm5)
      !> The container acts on CM5-shifted charges
      logical, intent(in) :: cm5

      if (has_cm5 .and. (cm5 .neqv. use_cm5)) then
         call fatal_error(error, "Analytic second derivatives require the Born and "// &
            & "surface terms to agree on the CM5 charges")
      end if
      has_cm5 = .true.
      use_cm5 = cm5
   end subroutine check_cm5

end subroutine check_solvation_support

subroutine get_leaf_containers(list, leaves)
   !***********************************
   !* Copies of all containers in a list, with nested lists, like the ones
   !* created by the C API for a solvation model, replaced by their content.
   !***********************************
   !> List of interaction containers
   class(container_list), intent(in) :: list
   !> Copies of the containers which are not lists themselves
   type(container_leaf), allocatable, intent(out) :: leaves(:)

   integer :: nleaf

   allocate(leaves(list%get_size()))
   nleaf = 0
   call add_leaf_containers(list, leaves, nleaf)
   leaves = leaves(:nleaf)
end subroutine get_leaf_containers

recursive subroutine add_leaf_containers(list, leaves, nleaf)
   !> List of interaction containers
   class(container_list), intent(in) :: list
   !> Copies of the containers which are not lists themselves
   type(container_leaf), allocatable, intent(inout) :: leaves(:)
   !> Number of collected containers
   integer, intent(inout) :: nleaf

   integer :: ic, il
   logical :: nested
   class(container_type), allocatable :: cont
   type(container_leaf), allocatable :: tmp(:)

   do ic = 1, list%get_size()
      call list%get_container(ic, cont)
      nested = .false.
      select type(cont)
      class is(container_list)
         nested = .true.
         call add_leaf_containers(cont, leaves, nleaf)
      end select
      if (nested) cycle

      if (nleaf >= size(leaves)) then
         allocate(tmp(2*size(leaves) + 1))
         do il = 1, nleaf
            call move_alloc(leaves(il)%raw, tmp(il)%raw)
         end do
         call move_alloc(tmp, leaves)
      end if
      nleaf = nleaf + 1
      call move_alloc(cont, leaves(nleaf)%raw)
   end do
end subroutine add_leaf_containers

function container_label(cont) result(label)
   !> Interaction container
   class(container_type), intent(in) :: cont
   !> Label of the container, or a placeholder
   character(len=:), allocatable :: label

   if (allocated(cont%label)) then
      label = cont%label
   else
      label = "unknown"
   end if
end function container_label

subroutine new_solvation_kernel(self, mol, calc, error)
   !***********************************
   !* Collect the geometry dependent quantities of the solvation models of the
   !* calculator, which are shared by all derivatives at this structure.
   !***********************************
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Single-point calculator
   type(xtb_calculator), intent(in) :: calc
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: nat, ic, iat
   real(wp) :: adet
   real(wp), allocatable :: qunit(:)
   type(container_leaf), allocatable :: leaves(:)

   call check_solvation_support(calc, error)
   if (allocated(error)) return
   if (.not.allocated(calc%interactions)) return

   call get_leaf_containers(calc%interactions, leaves)
   do ic = 1, size(leaves)
      select type(cont => leaves(ic)%raw)
      type is(alpb_solvation)
         allocate(self%born, source=cont)
         self%use_cm5 = cont%useCM5
      type is(cds_solvation)
         allocate(self%cds, source=cont)
         self%use_cm5 = cont%useCM5
      end select
   end do

   nat = mol%nat
   allocate(self%kmat(nat, nat), source=0.0_wp)

   if (self%use_cm5) then
      allocate(self%cm5(nat), self%dcm5dr(3, nat, nat))
      call get_cm5_charges(mol, self%cm5, self%dcm5dr)
   end if

   if (allocated(self%born)) then
      associate(born => self%born)
         allocate(self%rad(nat), self%draddr(3, nat, nat))
         call born%gbobc%get_rad(mol, self%rad, self%draddr)
         select case(born%kernel)
         case(born_kernel%p16)
            call add_born_mat_p16(nat, mol%xyz, born%keps, self%rad, self%kmat)
         case(born_kernel%still)
            call add_born_mat_still(nat, mol%xyz, born%keps, self%rad, self%kmat)
         end select

         !>-- the shape correction depends on the total charge only
         allocate(self%dkshape(3, nat), source=0.0_wp)
         if (born%alpbet > 0.0_wp) then
            call get_adet(nat, mol%xyz, born%gbobc%vdwr, adet)
            self%kshape = born%keps*born%alpbet/adet
            self%kmat(:, :) = self%kmat + self%kshape
            !>-- the derivative of kshape/2 is the gradient at unit total charge
            allocate(qunit(nat), source=0.0_wp)
            qunit(1) = 1.0_wp
            call get_adet_deriv(nat, mol%xyz, born%gbobc%vdwr, born%keps*born%alpbet, &
               & qunit, self%dkshape)
            self%dkshape(:, :) = 2*self%dkshape
         end if
      end associate
   end if

   !>-- a geometry dependent surface tension is not provided by the CDS model
   if (allocated(self%cds)) then
      associate(cds => self%cds)
         allocate(self%surface(nat), self%dsdr(3, nat, nat))
         call cds%sasa%get_surface(mol, self%surface, self%dsdr)
         self%tension = cds%tension(mol%id)
         if (allocated(cds%hbond)) then
            self%hbond = cds%hbond(mol%id)
            do iat = 1, nat
               self%kmat(iat, iat) = self%kmat(iat, iat) + 2*self%hbond(iat)*self%surface(iat)
            end do
         end if
      end associate
   end if
end subroutine new_solvation_kernel

pure function get_shifted_charges(self, qat) result(qeff)
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(in) :: self
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Charges the solvation models act on
   real(wp) :: qeff(size(qat))

   qeff(:) = qat
   if (self%use_cm5) qeff(:) = qeff + self%cm5
end function get_shifted_charges

subroutine get_solvation_potential(self, qat, vat)
   !***********************************
   !* Atom-resolved solvation potential v = K Q at the atomic charges q.
   !***********************************
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(in) :: self
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Atom-resolved potential
   real(wp), intent(out) :: vat(:)

   real(wp), allocatable :: qeff(:)

   allocate(qeff(size(qat)))
   qeff(:) = get_shifted_charges(self, qat)
   call symv(self%kmat, qeff, vat)
end subroutine get_solvation_potential

subroutine get_solvation_potential_derivative(self, mol, qat, dvat)
   !***********************************
   !* Derivative of the atom-resolved solvation potential w.r.t. the atomic
   !* positions at fixed charges q, dvat(A, c, B) = dv_A/dR_cB,
   !*   dv/dR = dK/dR Q + K dc/dR.
   !***********************************
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Derivative of the atom-resolved potential
   real(wp), intent(out) :: dvat(:, :, :)

   integer :: nat
   real(wp), allocatable :: qeff(:), dvflat(:, :)

   nat = mol%nat
   allocate(qeff(size(qat)))
   qeff(:) = get_shifted_charges(self, qat)
   call get_solvation_kernel_derivative(self, mol, qeff, dvat)

   if (self%use_cm5) then
      allocate(dvflat(nat, 3*nat))
      call gemm(self%kmat, reshape(self%dcm5dr, [3*nat, nat]), dvflat, transb="t")
      dvat(:, :, :) = dvat + reshape(dvflat, [nat, 3, nat])
   end if
end subroutine get_solvation_potential_derivative

subroutine get_solvation_kernel_derivative(self, mol, xat, dvat)
   !***********************************
   !* Derivative of the response kernel w.r.t. the atomic positions applied
   !* to atomic charges x, dvat(A, c, B) = sum_C dK_AC/dR_cB x_C.
   !* The Born matrix depends on the positions directly and via the Born
   !* radii, whose contribution is collected in W_AX = sum_C dJ_AC/da_X x_C
   !* and contracted once with the Jacobian of the radii.
   !***********************************
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic charges the kernel acts on
   real(wp), intent(in) :: xat(:)
   !> Derivative of the kernel applied to the charges
   real(wp), intent(out) :: dvat(:, :, :)

   integer :: nat, iat, jat
   real(wp) :: vec(3), r1, ehat(3), rho, da, db, daa, dab, dbb
   real(wp) :: dr, drho, drr, drrho, drhorho, xtot
   real(wp), allocatable :: wmat(:, :), dvflat(:, :)

   nat = mol%nat
   dvat(:, :, :) = 0.0_wp

   if (allocated(self%born)) then
      associate(born => self%born)
         allocate(wmat(nat, nat), source=0.0_wp)
         allocate(dvflat(nat, 3*nat))
         do iat = 1, nat
            do jat = 1, iat - 1
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
               r1 = norm2(vec)
               ehat(:) = vec / r1
               call get_radii_combination(born%kernel, self%rad(iat), self%rad(jat), &
                  & rho, da, db, daa, dab, dbb)
               call get_kernel_derivs(born%kernel, r1, rho, dr, drho, drr, drrho, drhorho)
               dr = born%keps*dr
               dvat(iat, :, iat) = dvat(iat, :, iat) + dr*xat(jat)*ehat
               dvat(iat, :, jat) = dvat(iat, :, jat) - dr*xat(jat)*ehat
               dvat(jat, :, iat) = dvat(jat, :, iat) + dr*xat(iat)*ehat
               dvat(jat, :, jat) = dvat(jat, :, jat) - dr*xat(iat)*ehat
               wmat(iat, iat) = wmat(iat, iat) + born%keps*drho*da*xat(jat)
               wmat(iat, jat) = wmat(iat, jat) + born%keps*drho*db*xat(jat)
               wmat(jat, iat) = wmat(jat, iat) + born%keps*drho*da*xat(iat)
               wmat(jat, jat) = wmat(jat, jat) + born%keps*drho*db*xat(iat)
            end do
            wmat(iat, iat) = wmat(iat, iat) - born%keps*xat(iat)/self%rad(iat)**2
         end do
         call gemm(wmat, reshape(self%draddr, [3*nat, nat]), dvflat, transb="t")
         dvat(:, :, :) = dvat + reshape(dvflat, [nat, 3, nat])

         xtot = sum(xat)
         do iat = 1, nat
            dvat(iat, :, :) = dvat(iat, :, :) + xtot*self%dkshape
         end do
      end associate
   end if

   if (allocated(self%hbond)) then
      do iat = 1, nat
         dvat(iat, :, :) = dvat(iat, :, :) &
            & + 2*self%hbond(iat)*xat(iat)*self%dsdr(:, :, iat)
      end do
   end if
end subroutine get_solvation_kernel_derivative

subroutine get_solvation_hessian(self, mol, qat, hessian)
   !***********************************
   !* Add the second derivative of the solvation energy w.r.t. the atomic
   !* positions at fixed charges q. The models are evaluated at Q = q + c, the
   !* CM5 corrections add with the potential v = K Q and its derivative dv/dR
   !*   v c'' + dc^T dv + dv^T dc - dc^T K dc.
   !***********************************
   !> Geometry dependent quantities of the solvation models
   type(solvation_kernel), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic charges
   real(wp), intent(in) :: qat(:)
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat
   real(wp), allocatable :: qeff(:), dEda(:), vat(:), dvat(:, :, :), dcdr(:, :), tmp(:, :)

   nat = mol%nat
   allocate(qeff(size(qat)))
   qeff(:) = get_shifted_charges(self, qat)

   if (allocated(self%born)) then
      associate(born => self%born)
         allocate(dEda(nat))
         call add_born_kernel_hessian(mol, born%kernel, born%keps, qeff, self%rad, hessian)
         call add_born_coupling_hessian(mol, born%kernel, born%keps, qeff, self%rad, &
            & self%draddr, hessian, dEda)
         call add_born_radii_hessian(born%gbobc, mol, dEda, hessian)
         call add_shape_hessian(mol, born%gbobc%vdwr, born%keps*born%alpbet, sum(qeff), &
            & hessian)
      end associate
   end if

   if (allocated(self%cds)) then
      if (allocated(self%hbond)) then
         call add_surface_hessian(self%cds%sasa, mol, self%tension + self%hbond*qeff**2, &
            & hessian)
      else
         call add_surface_hessian(self%cds%sasa, mol, self%tension, hessian)
      end if
   end if

   if (self%use_cm5) then
      allocate(vat(nat), dvat(nat, 3, nat), tmp(3*nat, nat))
      call get_solvation_potential(self, qat, vat)
      call get_solvation_potential_derivative(self, mol, qat, dvat)
      call add_cm5_hessian(mol, vat, hessian)

      dcdr = reshape(self%dcm5dr, [3*nat, nat])
      call gemm(dcdr, reshape(dvat, [nat, 3*nat]), hessian, beta=1.0_wp)
      call gemm(reshape(dvat, [nat, 3*nat]), dcdr, hessian, transa="t", transb="t", &
         & beta=1.0_wp)
      call gemm(dcdr, self%kmat, tmp)
      call gemm(tmp, dcdr, hessian, transb="t", alpha=-1.0_wp, beta=1.0_wp)
   end if
end subroutine get_solvation_hessian

end module tblite_hessian_solvation
