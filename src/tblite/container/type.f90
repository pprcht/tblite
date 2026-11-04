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

!> @file tblite/container/type.f90
!> Provides an abstract base class for interaction containers

!> Definition of abstract base class for general interactions
module tblite_container_type
   use mctc_env, only : wp, error_type, fatal_error
   use mctc_io, only : structure_type
   use tblite_container_cache, only : container_cache
   use tblite_partition, only : work_partition, serial_work_partition
   use tblite_response_type, only : kernel_type, response_type
   use tblite_scf_info, only : scf_info, not_used
   use tblite_scf_potential, only : potential_type
   use tblite_wavefunction_type, only : wavefunction_type
   implicit none
   private

   !> Abstract base class for interactions containers
   type, public, abstract :: container_type
      !> Label identifying this contribution
      character(len=:), allocatable :: label
      !> Share of the interaction loops evaluated by this instance
      type(work_partition) :: partition = serial_work_partition
   contains
      !> Assign an externally managed share of the interaction loops
      procedure :: set_partition
      !> Update container cache
      procedure :: update
      !> Get information about density dependent quantities used in the energy
      procedure :: variable_info
      !> Evaluate non-selfconsistent part of the interaction
      procedure :: get_engrad
      !> Evaluate selfconsistent energy of the interaction
      procedure :: get_energy
      !> Evaluate charge dependent potential shift from the interaction
      procedure :: get_potential
      !> Evaluate gradient of charge dependent potential shift from the interaction
      procedure :: get_potential_gradient
      !> Evaluate gradient contributions from the selfconsistent interaction
      procedure :: get_gradient
      !> Evaluate second derivative w.r.t. the density variables
      procedure :: get_response_kernel
      !> Evaluate mixed second derivative w.r.t. the density variables and positions
      procedure :: get_response_derivative
      !> Evaluate second derivative w.r.t. the positions at fixed density variables
      procedure :: get_hessian
      !> Information on container
      procedure :: info
   end type container_type


contains


!> Assign an externally managed share of the interaction loops
subroutine set_partition(self, partition)
   !> Instance of the interaction container
   class(container_type), intent(inout) :: self
   !> Share of the interaction loops evaluated by this instance
   type(work_partition), intent(in) :: partition

   self%partition = partition
end subroutine set_partition


!> Update container cache
subroutine update(self, mol, cache)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
end subroutine update


!> Evaluate non-selfconsistent part of the interaction
subroutine get_engrad(self, mol, cache, energies, gradient, sigma)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
   !> Interaction gradient
   real(wp), contiguous, intent(inout), optional :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout), optional :: sigma(:, :)
end subroutine get_engrad


!> Evaluate selfconsistent energy of the interaction
subroutine get_energy(self, mol, cache, wfn, energies)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Interaction energy
   real(wp), intent(inout) :: energies(:)
end subroutine get_energy


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential


!> Evaluate charge dependent potential shift from the interaction
subroutine get_potential_gradient(self, mol, cache, wfn, pot)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Density dependent potential
   type(potential_type), intent(inout) :: pot
end subroutine get_potential_gradient


!> Evaluate gradient contributions from the selfconsistent interaction
subroutine get_gradient(self, mol, cache, wfn, gradient, sigma)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Interaction gradient
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   !> Interaction virial
   real(wp), contiguous, intent(inout) :: sigma(:, :)
end subroutine get_gradient


subroutine get_response_kernel(self, mol, cache, wfn, kernel, error, dwfn)
   !***********************************
   !* Add the second derivative of the energy w.r.t. the charges and
   !* multipole moments at the state in wfn to the kernel. With dwfn, the
   !* derivative of the kernel along the first-order state dwfn is added
   !* instead. The default covers interactions without density dependence
   !* and reports an error for all others.
   !***********************************
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Response kernel
   type(kernel_type), intent(inout) :: kernel
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> First-order change of the state
   type(wavefunction_type), intent(in), optional :: dwfn

   if (uses_density(self%variable_info())) call missing_derivative(self, error)
end subroutine get_response_kernel


subroutine get_response_derivative(self, mol, cache, wfn, response, error)
   !***********************************
   !* Mixed second derivative of the energy w.r.t. the charges and multipole
   !* moments and the atomic positions at the state in wfn. Added to response
   !* are its contraction with the perturbed coordinates (derivative of the
   !* potential at fixed moments) and with the first-order moments (response
   !* part of the Hessian columns). The default covers interactions without
   !* density dependence and reports an error for all others.
   !***********************************
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> First-order quantities of the perturbations
   type(response_type), intent(inout) :: response
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   if (uses_density(self%variable_info())) call missing_derivative(self, error)
end subroutine get_response_derivative


subroutine get_hessian(self, mol, cache, wfn, hessian, error)
   !***********************************
   !* Add the second derivative of the energy w.r.t. the atomic positions at
   !* fixed charges and multipole moments to the Hessian, ordered as
   !* hessian(3*(jat-1)+jc, 3*(iat-1)+ic). The default reports an error.
   !***********************************
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Cached data between different runs
   type(container_cache), intent(inout) :: cache
   !> Wavefunction data
   type(wavefunction_type), intent(in) :: wfn
   !> Molecular Hessian, dimension 3*nat by 3*nat
   real(wp), contiguous, intent(inout) :: hessian(:, :)
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call missing_derivative(self, error)
end subroutine get_hessian


pure function uses_density(info) result(used)
   !***********************************
   !* Whether the energy depends on any of the density variables.
   !***********************************
   !> Information on the density dependence
   type(scf_info), intent(in) :: info
   !> At least one density variable is used
   logical :: used

   used = info%charge /= not_used .or. info%dipole /= not_used &
      & .or. info%quadrupole /= not_used
end function uses_density


subroutine missing_derivative(self, error)
   !***********************************
   !* Report that a container provides no analytic second derivatives.
   !***********************************
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call fatal_error(error, "Analytic second derivatives are not available for '"// &
      & self%info(0, "")//"'")
end subroutine missing_derivative


!> Get information about density dependent quantities used in the energy
pure function variable_info(self) result(info)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Information on the required potential data
   type(scf_info) :: info

   info = scf_info()
end function variable_info


!> Information on container
pure function info(self, verbosity, indent) result(str)
   !> Instance of the interaction container
   class(container_type), intent(in) :: self
   !> Verbosity level
   integer, intent(in) :: verbosity
   !> Indentation level
   character(len=*), intent(in) :: indent
   !> Information on the container
   character(len=:), allocatable :: str

   if (allocated(self%label)) then
      str = self%label
   else
      str = "Unknown"
   end if
end function info


end module tblite_container_type
