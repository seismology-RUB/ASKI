!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief auxiliary fortran module conaining all derived-types stuff (to be used by module shore_f2py in file create_shore_lines_f2py.f90)
!!
!! \author Florian Schumacher
!! \date May 2016
!! 
!! \copyright &copy; GNU Public License
!
module create_shore_lines_f2pyVar
  use inversionGrid
  use errorMessage
  implicit none
  save

  type (inversion_grid) :: invgrid
  type (error_message) :: errmsg

  logical :: invgrid_created = .false.

end module create_shore_lines_f2pyVar
