# This file is part of tblite.
# SPDX-Identifier: LGPL-3.0-or-later
#
# tblite is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# tblite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License

set(dir "${CMAKE_CURRENT_SOURCE_DIR}")

list(
  APPEND srcs
  "${dir}/coulomb.f90"
  "${dir}/dispersion.f90"
  "${dir}/halogen.f90"
  "${dir}/hamiltonian.f90"
  "${dir}/multipole.f90"
  "${dir}/pair.f90"
  "${dir}/polarizability.f90"
  "${dir}/reference.f90"
  "${dir}/repulsion.f90"
  "${dir}/xtb.f90"
)

set(srcs "${srcs}" PARENT_SCOPE)
