!----------------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 0.3.
!
!   ASKI version 0.3 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 0.3 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 0.3.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!> \brief General module to read in parameters from a parameter file
!-------------------------------------------------------------------------------
 module inputParameter
	use errorMessage
	implicit none
	interface readInputParameter
		module procedure readFunctionInputParameter
	end interface
	interface dealloc; module procedure deallocInputParameter; end interface
	interface print; module procedure printInputParameter; end interface
	! keep .value. operator to make this module compatible with old programs
	interface operator (.value.); module procedure getRealValueFromKeywordInputParameter; end interface
	interface operator (.ival.); module procedure getIntegerValueFromKeywordInputParameter; end interface
	interface ival; module procedure getIostatIntegerValueFromKeywordInputParameter; end interface
	interface operator (.rval.); module procedure getRealValueFromKeywordInputParameter; end interface
	interface rval; module procedure getIostatRealValueFromKeywordInputParameter; end interface
	interface operator (.dval.); module procedure getDoubleValueFromKeywordInputParameter; end interface
	interface dval; module procedure getIostatDoubleValueFromKeywordInputParameter; end interface
	interface operator (.cval.); module procedure getComplexValueFromKeywordInputParameter; end interface
	interface cval; module procedure getIostatComplexValueFromKeywordInputParameter; end interface
	interface operator (.dcval.); module procedure getDoubleComplexValueFromKeywordInputParameter; end interface
	interface dcval; module procedure getIostatDoubleComplexValueFromKeywordInputParameter; end interface
	interface operator (.sval.); module procedure getStringValueFromKeywordInputParameter; end interface
	interface sval; module procedure getIostatStringValueFromKeywordInputParameter; end interface
	interface operator (.lval.); module procedure getLogicalValueFromKeywordInputParameter; end interface
	interface lval; module procedure getIostatLogicalValueFromKeywordInputParameter; end interface
	interface ivec; module procedure getIostatIntegerVectorFromKeywordInputParameter; end interface
	interface ivecp; module procedure getIostatPointerIntegerVectorFromKeywordInputParameter; end interface
	interface rvec; module procedure getIostatRealVectorFromKeywordInputParameter; end interface
	interface rvecp; module procedure getIostatPointerRealVectorFromKeywordInputParameter; end interface
	interface cvec; module procedure getIostatComplexVectorFromKeywordInputParameter; end interface
	interface cvecp; module procedure getIostatPointerComplexVectorFromKeywordInputParameter; end interface
	interface svec; module procedure getIostatStringVectorFromKeywordInputParameter; end interface
	interface svecp; module procedure getIostatPointerStringVectorFromKeywordInputParameter; end interface
	type input_parameter
		private
		integer :: npar
		character (len=80), dimension(:), pointer :: keyword => null()
		character (len=300), dimension(:), pointer :: value => null()
	end type input_parameter
!
 contains
!-------------------------------------------------------------------------------
!> \brief Create keywords in input_parameter object
!
	subroutine createKeywordsInputParameter(this,word)
	type (input_parameter) :: this
	character (len=*), dimension(:) :: word
	integer :: i
	this%npar = size(word)
	allocate(this%keyword(this%npar))
	allocate(this%value(this%npar))
	do i = 1,this%npar
		this%keyword(i) = trim(adjustl(word(i)))
	enddo
	end subroutine createKeywordsInputParameter
!-------------------------------------------------------------------------------
!> \brief Read input parameters from parameter file
!
	function readFunctionInputParameter(this,lu,parfile) result(errmsg)
	type (input_parameter) :: this
	integer :: lu
	character (len=*) :: parfile
	type (error_message) :: errmsg
	character (len=300) :: line,line_adjl
	character (len=400) :: errstr
	logical, dimension(:), allocatable :: foundkey
	character (len=18) :: myname = 'readInputParameter'
	integer :: ios,i,keyidx,eqindx,linelen,iline
!
	call new(errmsg,myname)
!
!  open file
!
	open(lu,file = parfile,status = 'old',iostat = ios)
	if (ios /= 0) then
		call add(errmsg,2,"Parameter file '"//trim(parfile)//"' could not be opened",myname)
		return
	endif
!
!  read in file line by line
!
	allocate(foundkey(this%npar))
	ios = 0
	foundkey = .false.
	iline = 0
	do while (ios == 0)
		read(lu,'(a)',iostat = ios) line
		iline = iline + 1
		if (ios < 0) then
			cycle
		else if (ios > 0) then
			write(errstr,*) "could not read line ",iline," of parameter file '"//trim(parfile)//"'"
			call add(errmsg,2,trim(errstr),myname)
			deallocate(foundkey); close(lu)
			return
		endif
		line_adjl = adjustl(line)
		if (line_adjl(1:1) == '#') cycle        ! ignore comment lines
		linelen = len_trim(line)
		if(linelen == 0) cycle                  ! ignore empty lines (or with spaces only)
		do i = 1,this%npar
			keyidx = index(line,trim(this%keyword(i)))
			if (keyidx > 0) then
				eqindx = index(line,'=')
				if ( trim(line(keyidx:eqindx-1)) == trim(this%keyword(i)) .and. (.not.foundkey(i)) ) then
					this%value(i) = trim(adjustl(line(eqindx+1:linelen)))
					foundkey(i) = .true.
					exit
				endif
			endif ! keyidx > 0
		enddo ! i
	enddo ! while (ios == 0)
!
!  check if all requested keywords were found
!
	do i = 1,this%npar
		if(.not.foundkey(i)) then
			write(errstr,*) "requested keyword '"//trim(this%keyword(i))//&
			  "' was not found in parameter file '"//trim(parfile)//"'"
			call add(errmsg,2,trim(errstr),myname)
		endif
	enddo ! i
	close(lu)
	deallocate(foundkey)
	end function readFunctionInputParameter
!-------------------------------------------------------------------------------
!> \brief Read input parameters from parameter file
!
	subroutine readSubroutineInputParameter(this,lu,parfile,errmsg)
	type (input_parameter) :: this
	integer :: lu
	character (len=*) :: parfile
	type (error_message) :: errmsg
	character (len=300) :: line,line_adjl
	character (len=400) :: errstr
	logical, dimension(:), allocatable :: foundkey
	character (len=18) :: myname = 'readInputParameter'
	integer :: ios,i,keyidx,eqindx,linelen,iline
!
	call addTrace(errmsg,myname)
!
!  open file
!
	open(lu,file = parfile,status = 'old',iostat = ios)
	if (ios /= 0) then
		call add(errmsg,2,"Parameter file '"//trim(parfile)//"' could not be opened",myname)
		return
	endif
!
!  read in file line by line
!
	allocate(foundkey(this%npar))
	ios = 0
	foundkey = .false.
	iline = 0
	do while (ios == 0)
		read(lu,'(a)',iostat = ios) line
		iline = iline + 1
		if (ios < 0) then
			cycle
		else if (ios > 0) then
			write(errstr,*) "could not read line ",iline," of parameter file '"//trim(parfile)//"'"
			call add(errmsg,2,trim(errstr),myname)
			deallocate(foundkey); close(lu)
			return
		endif
		line_adjl = adjustl(line)
		if (line_adjl(1:1) == '#') cycle        ! ignore comment lines
		linelen = len_trim(line)
		if(linelen == 0) cycle                  ! ignore empty lines (or with spaces only)
		do i = 1,this%npar
			keyidx = index(line,trim(this%keyword(i)))
			if (keyidx > 0) then
				eqindx = index(line,'=')
				if ( trim(line(keyidx:eqindx-1)) == trim(this%keyword(i)) .and. (.not.foundkey(i)) ) then
					this%value(i) = trim(adjustl(line(eqindx+1:linelen)))
					foundkey(i) = .true.
					exit
				endif
			endif ! keyidx > 0
		enddo ! i
	enddo ! while (ios == 0)
!
!  check if all requested keywords were found
!
	do i = 1,this%npar
		if(.not.foundkey(i)) then
			write(errstr,*) "requested keyword '"//trim(this%keyword(i))//&
			  "' was not found in parameter file '"//trim(parfile)//"'"
			call add(errmsg,2,trim(errstr),myname)
		endif
	enddo ! i
	close(lu)
	deallocate(foundkey)
	end subroutine readSubroutineInputParameter
!-------------------------------------------------------------------------------
!> \brief Deallocate inputParameter object
!
	subroutine deallocInputParameter(this)
	type (input_parameter) :: this
	if (associated(this%keyword)) deallocate(this%keyword)
	if (associated(this%value)) deallocate(this%value)
	end subroutine deallocInputParameter
!-------------------------------------------------------------------------------
!> \brief get integer value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res integer value for keyword word
!! \return integer value read from value string
!
	function getIntegerValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer :: res
	integer :: i,ios
	res = 0
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getIntegerValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get integer value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res integer value for keyword word
!! \return integer value read from value string
!
	function getIostatIntegerValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = 0
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatIntegerValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get real value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res real value for keyword word
!! \return real value read from value string
!
	function getRealValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	real :: res
	integer :: i,ios
	res = 0.
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getRealValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get real value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res real value for keyword word
!! \return real value read from value string
!
	function getIostatRealValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	real :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = 0.
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatRealValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get double precision value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res double precision value for keyword word
!! \return double precision value read from value string
!
	function getDoubleValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	double precision :: res
	integer :: i,ios
	res = 0.d0
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getDoubleValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get double precision value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res double precision value for keyword word
!! \return double precision value read from value string
!
	function getIostatDoubleValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	double precision :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = 0.d0
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatDoubleValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get complex value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res complex value for keyword word
!! \return complex value read from value string
!
	function getComplexValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	complex :: res
	integer :: i,ios
	res = (0.,0.)
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getComplexValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get complex value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res complex value for keyword word
!! \return complex value read from value string
!
	function getIostatComplexValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	complex :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = (0.,0.)
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatComplexValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get double complex value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res double complex value for keyword word
!! \return double complex value read from value string
!
	function getDoubleComplexValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	double complex :: res
	integer :: i,ios
	res = (0.d0,0.d0)
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getDoubleComplexValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get double complex value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res double complex value for keyword word
!! \return double complex value read from value string
!
	function getIostatDoubleComplexValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	double complex :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = (0.d0,0.d0)
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatDoubleComplexValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get string value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res string value for keyword word
!! \return real string read from value string
!
	function getStringValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	character (len=300) :: res
	integer :: i
	res = ''
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			res = trim(this%value(i))
			return
		endif
	enddo
	end function getStringValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get string value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res string value for keyword word
!! \return real string read from value string
!
	function getIostatStringValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	character (len=300) :: res
	integer :: i
	if(present(iostat)) iostat = -1
	res = ''
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			res = trim(this%value(i))
			if(present(iostat)) iostat = 0 ! no sensible iostat handling here, just indicate if the keyword is present
			return
		endif
	enddo
	end function getIostatStringValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get logical value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param res logical value for keyword word
!! \return logical value read from value string
!
	function getLogicalValueFromKeywordInputParameter(this,word) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	logical :: res
	integer :: i,ios
	res = .false.
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			return
		endif
	enddo
	end function getLogicalValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get logical value for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param res logical value for keyword word
!! \return logical value read from value string
!
	function getIostatLogicalValueFromKeywordInputParameter(this,word,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	logical :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = .false.
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatLogicalValueFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get integer 1D array for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res integer 1D array for keyword word
!! \return integer 1D array read from value string
!
	function getIostatIntegerVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	integer, dimension(n) :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = 0
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatIntegerVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get pointer to integer 1D array for given keyword
!! \details There is memory allocated for the return pointer within this function.
!! Hence, this function should not be used implicitely ( as in "write(*,*) ivecp(inpar,'keyword',5)" ),
!! because then there is no way of deallocating the memory. Only use as in "p => ivecp(inpar,'keyword',5)".
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res pointer to integer 1D array for keyword word
!! \return pointer to integer 1D array read from value string
!
	function getIostatPointerIntegerVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	integer, dimension(:), pointer :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res => null()
	if(n.le.0) return
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			allocate(res(n))
			read(this%value(i),*,iostat=ios) res
			if(ios/=0) deallocate(res)
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatPointerIntegerVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get real 1D array for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res real 1D array for keyword word
!! \return real 1D array read from value string
!
	function getIostatRealVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	real, dimension(n) :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = 0.
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatRealVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get pointer to real 1D array for given keyword
!! \details There is memory allocated for the return pointer within this function.
!! Hence, this function should not be used implicitely ( as in "write(*,*) rvecp(inpar,'keyword',5)" ),
!! because then there is no way of deallocating the memory. Only use as in "p => rvecp(inpar,'keyword',5)".
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res pointer to real 1D array for keyword word
!! \return pointer to real 1D array read from value string
!
	function getIostatPointerRealVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	real, dimension(:), pointer :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res => null()
	if(n.le.0) return
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			allocate(res(n))
			read(this%value(i),*,iostat=ios) res
			if(ios/=0) deallocate(res)
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatPointerRealVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get complex 1D array for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res complex 1D array for keyword word
!! \return complex 1D array read from value string
!
	function getIostatComplexVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	complex, dimension(n) :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = (0.,0.)
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatComplexVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get pointer to complex 1D array for given keyword
!! \details There is memory allocated for the return pointer within this function.
!! Hence, this function should not be used implicitely ( as in "write(*,*) cvecp(inpar,'keyword',5)" ),
!! because then there is no way of deallocating the memory. Only use as in "p => cvecp(inpar,'keyword',5)".
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res pointer to complex 1D array for keyword word
!! \return pointer to complex 1D array read from value string
!
	function getIostatPointerComplexVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	complex, dimension(:), pointer :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res => null()
	if(n.le.0) return
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			allocate(res(n))
			read(this%value(i),*,iostat=ios) res
			if(ios/=0) deallocate(res)
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatPointerComplexVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get 1D array of character for given keyword
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res 1D array of character for keyword word
!! \return 1D array of character read from value string
!
!####################################################################################
!# ATTENTION:
!# this routine works fine, if all n substrings (members of 1D character array of size n)
!# are separated by blanks (no blanks allowed within a substring) and contain only letters
!# and digits.
!# THIS ROUTINE FAILS, due to command "read(this%value(i),*) res", if some substring contains
!# a "/"-character (others not tested /noticed yet)
!####################################################################################
	function getIostatStringVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	character (len=300), dimension(n) :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res = ''
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			read(this%value(i),*,iostat=ios) res
			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatStringVectorFromKeywordInputParameter
!-------------------------------------------------------------------------------
!> \brief get pointer to 1D array of character for given keyword
!! \details There is memory allocated for the return pointer within this function.
!! Hence, this function should not be used implicitely ( as in "write(*,*) cvecp(inpar,'keyword',5)" ),
!! because then there is no way of deallocating the memory. Only use as in "p => cvecp(inpar,'keyword',5)".
!! \param[in] this input_parameter object
!! \param word keyword
!! \param iostat optional integer, contain iostat returned by read() on exit
!! \param n size of 1D array
!! \param res pointer to 1D array of character for keyword word
!! \return pointer to 1D array of character read from value string
!
!####################################################################################
!# ATTENTION:
!# this routine works fine, if all n substrings (members of 1D character array of size n)
!# are separated by blanks (no blanks allowed within a substring) and contain only letters
!# and digits.
!# THIS ROUTINE FAILS, due to command "read(this%value(i),*) res", if some substring contains
!# a "/"-character (others not tested /noticed yet)
!####################################################################################
	function getIostatPointerStringVectorFromKeywordInputParameter(this,word,n,iostat) result(res)
	type (input_parameter), intent(in) :: this
	character (len=*), intent(in) :: word
	integer, optional :: iostat
	integer :: n
	character (len=300), dimension(:), pointer :: res
	integer :: i,ios
	if(present(iostat)) iostat = -1
	res => null()
	if(n.le.0) return
	do i = 1,this%npar
		if (trim(this%keyword(i)) == trim(adjustl(word))) then
			allocate(res(n))
			read(this%value(i),*,iostat=ios) res
			if(ios/=0) deallocate(res)
 			if(present(iostat)) iostat = ios
			return
		endif
	enddo
	end function getIostatPointerStringVectorFromKeywordInputParameter
!----------------------------------------------------------------------------------
! \brief Print overview of input parameters
!
	subroutine printInputParameter(this)
	type (input_parameter) :: this
	integer :: i
	print *,'The following values for the input parameters are used: '
	do i = 1,this%npar
		write(6,*) "'"//trim(this%keyword(i))//"' = '"//trim(this%value(i)),"'"
	enddo
	end subroutine printInputParameter
!
 end module
