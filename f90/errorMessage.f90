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
!> \brief  Module to handle errors occuring in functions or subroutines
!> \par Description
!> Functions producing exceptions may return an error_message object to
!! specify the name of the function where the exception happened, to
!! provide a description of the error and to rate its severity
!! ranging from OK (no error), warning to error.
!<--------------------------------------------------------------------- 
 module errorMessage
       use iso_fortran_env
	use realloc
	implicit none
	interface new
		module procedure createErrorMessage
		module procedure createOKErrorMessage
	end interface
	interface add; module procedure updateErrorMessage; end interface
	interface dealloc; module procedure deallocErrorMessage; end interface
	interface print; module procedure printErrorMessage; end interface
	interface addTrace; module procedure addTraceErrorMessage; end interface
	interface operator (.level.); module procedure getLevelErrorMessage; end interface
	type error_message
		private
		integer :: level = 0                 !< error level: 0 = success, 1 = warning, 2 = error
		character (len=400), dimension(:), pointer :: messages => null()    !< error/warning messages
		character (len=132) :: fctname = ''    !< function name where error message was created
		character (len=132), dimension(:), pointer :: trace => null()  !< function names through which error was propagated
	end type
!
 contains
!-----------------------------------------------------------------
!> \brief Create an error message
!
	subroutine createErrorMessage(this,level,message,fctname)
	type (error_message) :: this
	integer :: level
	character (len=*) :: message,fctname
	call deallocErrorMessage(this)
	this%level = level
	call addMessageErrorMessage(this,level,message,fctname)
	this%fctname = fctname
	end subroutine createErrorMessage
!------------------------------------------------------------------
!> \brief Create a default OK message
!
	subroutine createOKErrorMessage(this,fctname)
	type (error_message) :: this
	character (len=*) :: fctname
	call deallocErrorMessage(this)
	this%level = 0
	this%fctname = fctname
	end subroutine createOKErrorMessage
!-----------------------------------------------------------------
!> \brief Update an error message
!
	subroutine updateErrorMessage(this,level,message,fctname)
	type (error_message) :: this
	integer :: level
	character (len=*) :: message,fctname
	! error levels can only increase! (from success to warning or error, from warning to error)
	if(level>this%level) this%level = level
	if(this%fctname == '') then
		! error message was actually not created yet
		this%fctname = fctname
	endif
	call addMessageErrorMessage(this,level,message,fctname)
	end subroutine updateErrorMessage
!-----------------------------------------------------------------
!> \brief Deallocate error message
!
	subroutine deallocErrorMessage(this)
	type (error_message) :: this
	if (associated(this%trace)) deallocate(this%trace)
	if (associated(this%messages)) deallocate(this%messages)
	this%level = 0
	this%fctname = ''
	end subroutine deallocErrorMessage
!------------------------------------------------------------------
!> \brief Print error message
!
	subroutine printErrorMessage(this)
	type (error_message), intent(in) :: this
	integer :: j
	write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	select case (this%level)
	case (0); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>   SUCCESS   >>>>>>>>>>>>>>>'
	case (1); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>   WARNING   >>>>>>>>>>>>>>>'
	case (2); write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>   ERROR   >>>>>>>>>>>>>>>>'
	end select
	if (associated(this%messages)) then
		write(OUTPUT_UNIT,*) '>>>> '
		do j=1,size(this%messages)
			write(OUTPUT_UNIT,*) '>>>> ',trim(this%messages(j))
		enddo
	endif
	write(OUTPUT_UNIT,*) '>>>> '
	write(OUTPUT_UNIT,*) '>>>> message created in --> ',trim(this%fctname)
	if (associated(this%trace)) then
		do j=1,size(this%trace)
			write(OUTPUT_UNIT,*) '>>>>     passed through --> ',trim(this%trace(j))
		enddo
	endif
	write(OUTPUT_UNIT,*) '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	write(OUTPUT_UNIT,*) ''
	end subroutine printErrorMessage
!---------------------------------------------------------------
!> \brief Add a trace to error message
!
	subroutine addTraceErrorMessage(this,fctname)
	type (error_message) :: this
	character (len=*) :: fctname
	integer :: n
	if (associated(this%trace)) then
		n = size(this%trace)
		this%trace => reallocate(this%trace,n+1)
		this%trace(n+1) = trim(fctname)
	else
		allocate(this%trace(1))
		this%trace(1) = trim(fctname)
	endif
	end subroutine addTraceErrorMessage
!---------------------------------------------------------------
!> \brief Add a message to error message object
!
	subroutine addMessageErrorMessage(this,level,message,fctname)
	type (error_message) :: this
	integer :: level
	character (len=*) :: message,fctname
	integer :: n
	character (len=400) :: string
	select case (level)
	case(0); string = 'SUCCESS'
	case(1); string = 'WARNING'
	case(2); string = 'ERROR'
	end select
	string = trim(string)//' in '//trim(fctname)//' --> '//trim(message)
	if (associated(this%messages)) then
		n = size(this%messages)
		this%messages => reallocate(this%messages,n+1)
		this%messages(n+1) = trim(string)
	else
		allocate(this%messages(1))
		this%messages(1) = trim(string)
	endif
	end subroutine addMessageErrorMessage
!---------------------------------------------------------------
!> \brief Get level of error message
!
	integer function getLevelErrorMessage(this)
	type (error_message), intent(in) :: this
	getLevelErrorMessage = this%level
	end function getLevelErrorMessage
!
 end module errorMessage
