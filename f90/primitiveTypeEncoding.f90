!--------------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 0.3.
!
!   Gemini II and ASKI version 0.3 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 0.3 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 0.3.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------------
!  conventions for identifying pritive data types
!----------------------------------------------------------------------
 module primitiveTypeEncoding
	implicit none
	integer, parameter :: T_INTEGER = 1
	integer, parameter :: T_REAL = 2
	integer, parameter :: T_DOUBLE = 3
	integer, parameter :: T_COMPLEX = 4
	integer, parameter :: T_DOUBLE_COMPLEX = 5
	integer, parameter :: T_CHAR = 6
	integer, parameter :: T_FLEXIBLE = 7
	integer, parameter :: T_LOGICAL = 8
end module

