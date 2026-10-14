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

!> @file tblite/hessian/hamiltonian.f90
!> Provides derivatives of the core Hamiltonian required for the Hessian

!> First and second derivatives of the extended Hückel core Hamiltonian,
!> H_mn = 1/2 K_mn (e_m + e_n) p_mn(r) S_mn, where the self-energies e depend
!> on the coordination numbers and p is the distance dependent scaling
module tblite_hessian_hamiltonian
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_ncoord, only : ncoord_type
   use tblite_basis_type, only : basis_type, cgto_type
   use tblite_cutoff, only : get_lattice_points
   use tblite_hessian_pair, only : get_pair_block, add_pair_hessian, add_ncoord_hessian
   use tblite_integral_native_integrals, only : overlap_grad_cgto, overlap_hess_cgto, msao
   use tblite_xtb_h0, only : tb_hamiltonian
   implicit none
   private

   public :: get_scaled_overlap, get_hamiltonian_derivative, get_hamiltonian_hessian

contains

pure subroutine get_shpoly(h0, izp, jzp, ish, jsh, vec, r2, sh, dsh, d2sh)
   !***********************************
   !* Distance dependent scaling p = (1 + k_i sqrt(r/R))(1 + k_j sqrt(r/R))
   !* of a shell pair and its Cartesian derivatives w.r.t. the distance vector.
   !***********************************
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Species of the two atoms
   integer, intent(in) :: izp, jzp
   !> Shells of the two atoms
   integer, intent(in) :: ish, jsh
   !> Distance vector
   real(wp), intent(in) :: vec(3)
   !> Squared distance
   real(wp), intent(in) :: r2
   !> Scaling function
   real(wp), intent(out) :: sh
   !> First derivative of the scaling function
   real(wp), intent(out) :: dsh(3)
   !> Second derivative of the scaling function
   real(wp), intent(out) :: d2sh(3, 3)

   real(wp) :: r1, rad, a, b

   r1 = sqrt(r2)
   rad = h0%rad(izp) + h0%rad(jzp)
   a = (h0%shpoly(ish, izp) + h0%shpoly(jsh, jzp)) / sqrt(rad)
   b = h0%shpoly(ish, izp) * h0%shpoly(jsh, jzp) / rad
   sh = 1.0_wp + a*sqrt(r1) + b*r1
   dsh(:) = (0.5_wp*a/sqrt(r1) + b) * vec/r1
   call get_pair_block(0.5_wp*a/sqrt(r1) + b, -0.25_wp*a/(r1*sqrt(r1)), r1, vec, d2sh)
end subroutine get_shpoly

subroutine get_scaled_overlap(mol, bas, h0, overlap, cutoff, bmat)
   !***********************************
   !* Overlap matrix scaled such that H_mn = (e_m + e_n) B_mn, i.e.
   !* B = 1/2 K p S between different atoms and B = 1/2 S within an atom.
   !* Only available for molecular systems.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Overlap matrix
   real(wp), intent(in) :: overlap(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Scaled overlap matrix
   real(wp), intent(out) :: bmat(:, :)

   integer :: iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, ni, nj
   real(wp) :: r2, vec(3), sh, dsh(3), d2sh(3, 3), fac

   bmat(:, :) = 0.0_wp
   do iat = 1, mol%nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, mol%nat
         jzp = mol%id(jat)
         js = bas%ish_at(jat)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(vec**2)
         if (r2 > cutoff**2) cycle
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is+ish)
            ni = msao(bas%cgto(ish, izp)%ang)
            do jsh = 1, bas%nsh_id(jzp)
               jj = bas%iao_sh(js+jsh)
               nj = msao(bas%cgto(jsh, jzp)%ang)
               fac = 0.5_wp
               if (iat /= jat) then
                  call get_shpoly(h0, izp, jzp, ish, jsh, vec, r2, sh, dsh, d2sh)
                  fac = 0.5_wp * h0%hscale(jsh, ish, jzp, izp) * sh
               end if
               bmat(jj+1:jj+nj, ii+1:ii+ni) = fac * overlap(jj+1:jj+nj, ii+1:ii+ni)
            end do
         end do
      end do
   end do
end subroutine get_scaled_overlap

subroutine get_hamiltonian_derivative(mol, bas, h0, bmat, selfenergy, dsedcn, dcndr, &
      & trans, cutoff, kat, dh0)
   !***********************************
   !* Derivative of the core Hamiltonian w.r.t. the Cartesian position of the
   !* atom kat. The derivative of the self-energies via the coordination
   !* numbers enters all blocks, the derivative of the scaled overlap only the
   !* blocks between kat and another atom.
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Scaled overlap matrix
   real(wp), intent(in) :: bmat(:, :)
   !> Self-energy of each shell
   real(wp), intent(in) :: selfenergy(:)
   !> Derivative of the self-energy w.r.t. the coordination number
   real(wp), intent(in) :: dsedcn(:)
   !> Derivative of the coordination numbers w.r.t. the atomic positions
   real(wp), intent(in) :: dcndr(:, :, :)
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Atom to differentiate for
   integer, intent(in) :: kat
   !> Derivative of the Hamiltonian, the last dimension is the direction
   real(wp), intent(out) :: dh0(:, :, :)

   integer :: iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, ic
   real(wp) :: r2, vec(3), cutoff2, sgn, sh, dsh(3), d2sh(3, 3), hij
   real(wp), allocatable :: dse(:, :), stmp(:), dstmp(:, :)

   !>-- derivative of the self-energy of each orbital
   allocate(dse(bas%nao, 3))
   do iao = 1, bas%nao
      dse(iao, :) = dsedcn(bas%ao2sh(iao)) * dcndr(:, kat, bas%ao2at(iao))
   end do
   do ic = 1, 3
      do iao = 1, bas%nao
         do jao = 1, bas%nao
            dh0(jao, iao, ic) = bmat(jao, iao) * (dse(jao, ic) + dse(iao, ic))
         end do
      end do
   end do

   cutoff2 = cutoff**2
   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2))
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
                  call get_shpoly(h0, izp, jzp, ish, jsh, vec, r2, sh, dsh, d2sh)
                  hij = 0.5_wp * h0%hscale(jsh, ish, jzp, izp) &
                     & * (selfenergy(is+ish) + selfenergy(js+jsh))
                  call overlap_grad_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dstmp)
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        dh0(jj+jao, ii+iao, :) = dh0(jj+jao, ii+iao, :) &
                           & + sgn*hij*(sh*dstmp(:, ij) + dsh*stmp(ij))
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do
end subroutine get_hamiltonian_derivative

subroutine get_hamiltonian_hessian(mol, bas, h0, ncoord, selfenergy, dsedcn, pmat, wmat, &
      & vao, dcndr, trans, cutoff, hessian)
   !***********************************
   !* Add the second derivative of the band energy at fixed density,
   !* Tr[P d2H0] - Tr[W d2S]. The self-energies depend on the positions via
   !* the coordination numbers only, the dependence on the charges is part
   !* of the response. The fixed potential of the self-consistent terms acts
   !* via the overlap and is absorbed in the energy-weighted density matrix,
   !* W_mn + 1/2 P_mn (v_m + v_n).
   !***********************************
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Basis set information
   type(basis_type), intent(in) :: bas
   !> Hamiltonian interaction data
   type(tb_hamiltonian), intent(in) :: h0
   !> Coordination number container
   class(ncoord_type), intent(in) :: ncoord
   !> Self-energy of each shell
   real(wp), intent(in) :: selfenergy(:)
   !> Derivative of the self-energy w.r.t. the coordination number
   real(wp), intent(in) :: dsedcn(:)
   !> Density matrix
   real(wp), intent(in) :: pmat(:, :)
   !> Energy-weighted density matrix
   real(wp), intent(in) :: wmat(:, :)
   !> Potential of the self-consistent terms for each orbital
   real(wp), intent(in) :: vao(:)
   !> Derivative of the coordination numbers w.r.t. the atomic positions
   real(wp), intent(in) :: dcndr(:, :, :)
   !> Lattice points within a given realspace cutoff
   real(wp), intent(in) :: trans(:, :)
   !> Realspace cutoff
   real(wp), intent(in) :: cutoff
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj
   integer :: iao, jao, nao, ij, ic, jc
   real(wp) :: r2, vec(3), cutoff2, sh, dsh(3), d2sh(3, 3), hscale, hij, pij, wij
   real(wp) :: block(3, 3), dfs(3)
   real(wp), allocatable :: stmp(:), dstmp(:, :), d2stmp(:, :, :)
   real(wp), allocatable :: dEdcn(:), dEdcndr(:, :, :), lattr(:, :)

   nat = mol%nat
   cutoff2 = cutoff**2
   allocate(stmp(msao(bas%maxl)**2), dstmp(3, msao(bas%maxl)**2), &
      & d2stmp(3, 3, msao(bas%maxl)**2))
   allocate(dEdcn(nat), dEdcndr(3, nat, nat), source=0.0_wp)

   do iat = 1, nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do jat = 1, iat - 1
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
                  call get_shpoly(h0, izp, jzp, ish, jsh, vec, r2, sh, dsh, d2sh)
                  hscale = 0.5_wp * h0%hscale(jsh, ish, jzp, izp)
                  hij = hscale * (selfenergy(is+ish) + selfenergy(js+jsh))
                  call overlap_hess_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dstmp, d2stmp)
                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao-1)
                        pij = pmat(jj+jao, ii+iao)
                        wij = wmat(jj+jao, ii+iao) &
                           & + 0.5_wp*pij*(vao(jj+jao) + vao(ii+iao))

                        !>-- direct second derivative of the scaled overlap
                        do jc = 1, 3
                           do ic = 1, 3
                              block(ic, jc) = block(ic, jc) - wij*d2stmp(ic, jc, ij) &
                                 & + pij*hij*(sh*d2stmp(ic, jc, ij) + d2sh(ic, jc)*stmp(ij) &
                                 & + dsh(ic)*dstmp(jc, ij) + dsh(jc)*dstmp(ic, ij))
                           end do
                        end do

                        !>-- dependence on the coordination numbers, the factor
                        !>   two accounts for both triangles of the matrices
                        pij = 2*pij*hscale
                        dfs(:) = sh*dstmp(:, ij) + dsh*stmp(ij)
                        dEdcn(iat) = dEdcn(iat) + pij*dsedcn(is+ish)*sh*stmp(ij)
                        dEdcn(jat) = dEdcn(jat) + pij*dsedcn(js+jsh)*sh*stmp(ij)
                        dEdcndr(:, iat, iat) = dEdcndr(:, iat, iat) + pij*dsedcn(is+ish)*dfs
                        dEdcndr(:, jat, iat) = dEdcndr(:, jat, iat) - pij*dsedcn(is+ish)*dfs
                        dEdcndr(:, iat, jat) = dEdcndr(:, iat, jat) + pij*dsedcn(js+jsh)*dfs
                        dEdcndr(:, jat, jat) = dEdcndr(:, jat, jat) - pij*dsedcn(js+jsh)*dfs
                     end do
                  end do
               end do
            end do
         end do
         call add_pair_hessian(hessian, iat, jat, 2*block)
      end do
   end do

   do iat = 1, nat
      izp = mol%id(iat)
      is = bas%ish_at(iat)
      do ish = 1, bas%nsh_id(izp)
         ii = bas%iao_sh(is+ish)
         do iao = 1, msao(bas%cgto(ish, izp)%ang)
            dEdcn(iat) = dEdcn(iat) + dsedcn(is+ish)*pmat(ii+iao, ii+iao)
         end do
      end do
   end do

   call get_lattice_points(mol%periodic, mol%lattice, ncoord%cutoff, lattr)
   call ncoord%add_coordination_number_hessian(mol, lattr, dEdcn, hessian)

   call add_ncoord_hessian(hessian, dcndr, reshape(dEdcndr, [3*nat, nat]))
end subroutine get_hamiltonian_hessian

end module tblite_hessian_hamiltonian
