!--------------------------------------------------------------------------
!   Copyright 2013 Thomas Forbriger () and Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------------------
!  Command line module
!
!  maxopt:  total number of possible options
!  maxmanarg:  total number of mandatory arguments
!  optid:   array of option names
!  optarg:  array of arguments to options
!  hasarg:  array of 1 and 0 indicating whether an option has an argument
!  optset:  array of 1 and 0 indicating which option was set on the command line
!  lastarg: index of last argument that belongs to the options
!           next argument is a mandatory one
!
!  for initialization use strings "opts" and "mask" to enter
!  a list of options and hasarg-flags:
!  opts:    string filled with blank separated options without minus signs
!  mask:    string with blank separated hasarg-flags
!--------------------------------------------------------------------
 module commandLine
	private :: clInit
	interface new
		module procedure clInit
	end interface
	interface dealloc
		module procedure clDealloc
	end interface
	interface clOptset
		module procedure clOptsetIndex
		module procedure clOptsetString
	end interface
	interface clOptarg
		module procedure clOptargIndex
		module procedure clOptargString
	end interface
	interface operator (.optset.); module procedure clOptsetIndex; end interface
	interface operator (.optarg.); module procedure clOptargIndex; end interface
	interface operator (.manarg.); module procedure clManarg; end interface
	type cmdLine
		private
		integer :: maxopt,maxmanarg,lastarg
		integer, dimension(:), pointer :: hasarg => null()
		integer, dimension(:), pointer :: optset => null()
		character (len=6), dimension(:), pointer :: optid => null()
		character (len=232), dimension(:), pointer :: optarg => null()
		character (len=232), dimension(:), pointer :: manarg => null()
	end type cmdLine
!
	contains
!-------------------------------------------------------------------------------
	subroutine clInit(this,maxopt,maxmanarg,opts,mask,printhelp,verbose)
	implicit none
	type (cmdLine) :: this
	integer :: maxopt, maxmanarg
	character (len=*) :: opts, mask
	external printhelp
	integer :: i,jarg,lastarg,jopt
	integer :: iargc
	character (len=7) argopt
	logical, optional :: verbose
	logical :: printflag
!
	if (present(verbose)) then; printflag = verbose; else; printflag = .false.; endif
!
	this%maxopt = maxopt
	this%maxmanarg = maxmanarg
!
	allocate(this%optid(maxopt),this%hasarg(maxopt),this%optarg(maxopt),&
	       & this%optset(maxopt))
	if (maxmanarg > 0) allocate(this%manarg(maxmanarg))
!
	read(opts,*,end=11) (this%optid(i),i=1,maxopt)
	read(mask,*,end=12) (this%hasarg(i),i=1,maxopt)
!
!  process command line
!
	do i=1,maxopt
		this%optset(i)=0
	enddo
	jarg=1
	lastarg=0
	do while (jarg <= iargc())           ! loop over optional command line arguments
		call getarg(jarg,argopt)
		jopt=1
		do while (jopt <= maxopt)     ! check which option has been set
			if(argopt == '-'//this%optid(jopt)) then
				this%optset(jopt)=1
				if(this%hasarg(jopt) == 1) then
					if(jarg+1 > iargc()) then
						print *,'commandLIne: Argument for option ',this%optid(jopt),' is missing!'
						call printhelp
						stop
					endif
					call getarg(jarg+1,this%optarg(jopt))
					jarg=jarg+1
				endif
				lastarg=jarg
				exit
			endif
			jopt=jopt+1
		enddo
		jarg=jarg+1
	enddo
	this%lastarg=lastarg
!
!  get mandatory arguments only if -h is not set and iargc() > 0
!
	if(this%optset(1) == 0 .and. iargc() > 0) then
		do i=1,maxmanarg
			if(i+lastarg > iargc()) then
				print *,'commandLine: not enough mandatory arguments'
				call printhelp
				stop
			endif
			call getarg(i+lastarg,this%manarg(i))
		enddo
	endif
!
!  optional output of options
!
	if (printflag) call clPrint(this)  
!
!  if -h is set, call printhelp and stop
!
	if(this%optset(1) == 1) then
		call printhelp
		stop
	endif
!
!  if iargc() == 0, trigger printhelp and stop
!
	if(iargc() == 0 .and. maxmanarg > 0) then
		call printhelp
		stop
	endif
	return
!
 11	print *,'commandLine: maxopt and option string are inconsistent'
 	stop
 12	print *,'commandLine: maxopt and hasarg string are inconsistent'
	stop
	end subroutine clInit
!-------------------------------------------------------------------------------
!  delete allocated pointer
!
	subroutine clDealloc(this)
	type (cmdLine) :: this
!
	if(associated(this%optid)) deallocate(this%optid)
	if(associated(this%hasarg)) deallocate(this%hasarg)
	if(associated(this%optarg)) deallocate(this%optarg)
	if(associated(this%optset)) deallocate(this%optset)
	if(associated(this%manarg)) deallocate(this%manarg)
	end subroutine clDealloc
!-------------------------------------------------------------------------------
!  return optset(i) by index
!
	logical function clOptsetIndex(this,n)
	type (cmdLine), intent(in) :: this
	integer, intent(in) :: n
	if(this%optset(n) == 1) then
		clOptsetIndex = .true.
	else
		clOptsetIndex = .false.
	endif
	end function clOptsetIndex
!-------------------------------------------------------------------------------
!  return optset(i) by option string
!
	logical function clOptsetString(this,string)
	type (cmdLine), intent(in) :: this
	character (len=*), intent(in) :: string
	integer :: i
	clOptsetString = .false.
	do i = 1,this%maxopt
		if (trim(this%optid(i)) == trim(string)) then
			clOptsetString = (this%optset(i) == 1)
			exit
		endif
	enddo
	end function clOptsetString
!-------------------------------------------------------------------------------
!  return optarg(i) by index
!
	character (len=232) function clOptargIndex(this,n)
	type (cmdLine), intent(in) :: this
	integer, intent(in) :: n
	clOptargIndex = this%optarg(n)
	end function clOptargIndex
!-------------------------------------------------------------------------------
!  return optarg(i) by string
!
	character (len=232) function clOptargString(this,string)
	type (cmdLine), intent(in) :: this
	character (len=*), intent(in) :: string
	integer :: i
	do i = 1,this%maxopt
		if (trim(this%optid(i)) == trim(string)) then
			clOptargString = this%optarg(i)
			exit
		endif
	enddo
	end function clOptargString
!-------------------------------------------------------------------------------
!  return manarg(i)
!
	character (len=232) function clManarg(this,n)
	type (cmdLine), intent(in) :: this
	integer, intent(in) :: n
	clManarg = this%manarg(n)
	end function clManarg
!---------------------------------------------------------------------------
!> \brief Print options with arguments to screen
!
	subroutine clPrint(this)
	type (cmdLine) :: this
	integer :: j
!
	print *,'-----------------------------------------------------------------------'
	print *,'                     Options and arguments'
	print *,''
	do j = 1,this%maxopt
		if (this%optset(j) == 1) print *,'Option: ',this%optid(j),' Argument: ',trim(this%optarg(j))
	enddo
	print *,'-----------------------------------------------------------------------'
	end subroutine
!-------------------------------------------------------------------------------
  end module commandLine
