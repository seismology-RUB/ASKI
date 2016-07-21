!--------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 1.1.
!
!   Gemini II and ASKI version 1.1 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 1.1 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 1.1.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!-----------------------------------------------------------------------
!> \brief Module to handle file units
!! Avoids opening of different files with same unit
!! Starts with a filled stack of file unit numbers.
!! A file unit number can be taken from the stack and
!! file units of closed file can be put back onto the stack
!! thus becoming available for the next request.
!-----------------------------------------------------------------------
 module fileUnitHandler
	use realloc
	implicit none
	interface new; module procedure createFileUnitHandler; end interface
	interface dealloc; module procedure deallocFileUnitHandler; end interface
	interface get; module procedure getFileUnitHandler; end interface
	interface add; module procedure addFileUnitHandler; end interface
	interface undo; module procedure undoLastGetFileUnitHandler; end interface
	type file_unit_handler
		private
		integer :: top                                                 ! number of top item on file unit stack
		integer :: maxitems                                            ! max number of registered file units
		integer :: nchunk                                              ! stack extension if empty
		integer, dimension(:), pointer :: funit => null()              ! array of units closed after use
	end type
!
 contains
!-----------------------------------------------------------------------
!> \brief Create a file unit handler object
!
	subroutine createFileUnitHandler(this,maxitems,nchunk)
	type (file_unit_handler) :: this
	integer :: maxitems
	integer, optional :: nchunk
	integer :: j
	allocate(this%funit(maxitems))
	this%funit = (/ (j, j = 11,maxitems+10) /)
	this%top = 1
	this%maxitems = maxitems
	if (present(nchunk)) then; this%nchunk = nchunk; else; this%nchunk = maxitems; endif
	end subroutine createFileUnitHandler
!-----------------------------------------------------------------------
!> \brief Deallocate file unit handler
!
	subroutine deallocFileUnitHandler(this)
	type (file_unit_handler) :: this
	if (associated(this%funit)) deallocate(this%funit)
	end subroutine deallocFileUnitHandler
!-----------------------------------------------------------------------
!> \brief Add a unit to the stack, ignore if stack is full
!
	subroutine addFileUnitHandler(this,lu)
	type (file_unit_handler) :: this
	integer :: lu
	if (this%top > 1) then
		this%funit(this%top-1) = lu
		this%top = this%top-1
	endif
	end subroutine addFileUnitHandler
!-----------------------------------------------------------------------
!> \brief Get a new file unit, increases stack automatically if empty
!
	function getFileUnitHandler(this) result(lu)
	type (file_unit_handler) :: this
	integer :: lu,j
!
	lu = this%funit(this%top)
	if (this%top == this%maxitems) then
		this%funit => reallocate(this%funit,this%top+this%nchunk)
		this%maxitems = this%maxitems+this%nchunk
		this%funit(this%top+1:this%maxitems) = (/ (j, j = this%top+1+10,this%maxitems+10) /)
	endif
	this%top = this%top+1
	end function getFileUnitHandler
!-------------------------------------------------------------------
!> \brief Undo the last get
!
	subroutine undoLastGetFileUnitHandler(this)
	type (file_unit_handler) :: this
	this%top = this%top-1
	end subroutine undoLastGetFileUnitHandler
!
 end module
