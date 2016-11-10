!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!   This file is part of ASKI version 1.2.
!
!   GEMINI_UNIFIED version 1.0 and ASKI version 1.2 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 and ASKI 1.2 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!  module for a string
!----------------------------------------------------------------------------
 module simpleString
    implicit none
    interface new; module procedure createSimpleString; end interface
    interface dealloc; module procedure deallocSimpleString; end interface
    interface operator (.length.); module procedure lengthSimpleString; end interface
    type simple_string
        private
        character (len=1), dimension(:), pointer :: zeichen => null()
    end type simple_string
!
 contains
!----------------------------------------------------------------------------
!  create a string object from a given character sequence
!  strip trailing blanks
!
    subroutine createSimpleString(this,string)
    type (simple_string) :: this
    character (len=*) :: string
    integer :: i
!
!  copy string onto zeichen
!
    allocate(this%zeichen(len_trim(string)))
    do i=1,len_trim(string)
        this%zeichen(i) = string(i:i)
    enddo
    end subroutine createSimpleString
!----------------------------------------------------------------------------
!  deallocate
!
    subroutine deallocSimpleString(this)
    type (simple_string) :: this
    if (associated(this%zeichen)) deallocate(this%zeichen)
    end subroutine deallocSimpleString
!----------------------------------------------------------------------------
!  get length of string
!
    integer function lengthSimpleString(this)
    type (simple_string), intent(in) :: this
    if(associated(this%zeichen)) then
        lengthSimpleString = size(this%zeichen)
    else
        lengthSimpleString = 0
    end if
    end function lengthSimpleString
!----------------------------------------------------------------------------
!  convert string back to character sequence
!
    subroutine convertToCharSimpleString(this,string)
    type (simple_string) :: this
    character (len=*) :: string
    integer :: n,i,size_this_zeichen
    if(associated(this%zeichen)) then
        size_this_zeichen = size(this%zeichen)
    else
        size_this_zeichen = 0
    end if
    n=min(size_this_zeichen,len(string))
    do i=1,n
        string(i:i) = this%zeichen(i)
    enddo
    do i=n+1,len(string)
        string(i:i) = ' '
    enddo
    end subroutine convertToCharSimpleString
!----------------------------------------------------------------------------    
 end module simpleString
