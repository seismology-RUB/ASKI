!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.0.
!
!   ASKI version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief module containing routines to read/write data from/to ascii files
!!
!! \details This module contains routines which open formatted
!!  files and read specific kinds of (numerical) data from file and return the data
!!  or write some given data to file.
!!  E.g. read in one (integer, real, complex) value per line of an ascii file and
!!  return that data vector to main program.
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module asciiDataIO
!
	use realloc
	use errorMessage
!
implicit none
!
!< overload all functions reading data from ascii files (select by type and rank of pointer p)
	interface readAsciiData
		!module procedure readIntegerVectorAsciiDataIO
		!module procedure readRealVectorAsciiDataIO
		!module procedure readDoublePrecisionVectorAsciiDataIO
		module procedure readComplexVectorAsciiDataIO
		!module procedure readIntegerMatrixAsciiDataIO
		module procedure readRealMatrixAsciiDataIO
		!module procedure readDoublePrecisionMatrixAsciiDataIO
		module procedure readComplexMatrixAsciiDataIO
	end interface
!< overload all functions writing data to ascii files (select by type and rank of array of incoming data)
	interface writeAsciiData
		module procedure writeComplexVectorAsciiDataIO
		module procedure writeComplexMatrixAsciiDataIO
		module procedure writeRealVectorAsciiDataIO
		module procedure writeRealMatrixAsciiDataIO
	end interface
!
contains
!------------------------------------------------------------------------
!> \brief read one complex value from each line of ascii file
!! \param filename ascii file name
!! \param lu file unit
!! \param ndata optional number of data values to be read from file
!! \param nskip optional number of lines to be skipped before reading the data
!! \param errmsg error message
!! \return error message
!
	function readComplexVectorAsciiDataIO(filename,lu,p,ndata,nskip) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	type (error_message) :: errmsg
	integer, optional :: ndata
	integer, optional :: nskip
	! outgoing
	complex, dimension(:), pointer :: p
	! local
	integer, parameter :: nchunk = 100 !< size by which data array is dynamically reallocated (if read all)
	complex :: z
	integer :: ios,ios_line,idata,iline
	logical :: dont_read_all
	character(len=400) :: errstr,line
	character(len=28) :: myname = 'readComplexVectorAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
	p => null()
!
	! define the required number of lines to read
	if(present(ndata)) then
		if(ndata<1) then
			write(errstr,*) 'desired number of data ( =',ndata,') is less than 1'
			call add(errmsg,2,trim(errstr),myname)
			return
		else
			dont_read_all = .true.
			write(errstr,*) 'trying to read ',ndata,' data values from file, one per line'
			call add(errmsg,0,trim(errstr),myname)
		endif
	else
		dont_read_all = .false.
		call add(errmsg,0,'trying to read all data values from the lines of the file',myname)
	endif
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	! if required, skip the first nskip lines of the file
	if(present(nskip)) then
		write(errstr,*) 'the first ',nskip,' lines will be skipped'
		call add(errmsg,0,trim(errstr),myname)
		do iline = 1,nskip
			! read new line from file
			read(lu,"(a400)",iostat=ios_line) line
			if(ios_line/=0) then
				! write maximum number of lines read from file to error message and return
				write(errstr,*) "Only ",iline-1," lines could be read from file. No data was read."
				call add(errmsg,1,trim(errstr),myname)
				return
			endif ! ios_line/=0
				
		enddo ! iline
	endif ! present(nskip)
!
	if(dont_read_all) then
		allocate(p(ndata))
	else
		allocate(p(nchunk))
	endif
	idata = 0
	iline = 0
	ios_line=0
!
	do while(ios_line==0)
		! counter on lines (for error messages)
		iline = iline+1
		! read new line from file
		read(lu,"(a400)",iostat=ios_line) line
		if(ios_line/=0) then
			! write maximum number of lines read from file to error message and terminate loop
			write(errstr,*) iline-1," lines could be read from file"
			call add(errmsg,0,trim(errstr),myname)
			exit
		endif ! ios_line/=0
!
		! try to read complex number from line string
		read(line,*,iostat=ios) z
		if(ios==0) then
			! reallocate if necessary
			if(idata==size(p)) p => reallocate(p,size(p)+nchunk)
			! add data value to array
			idata = idata+1
			p(idata) = z
		else
			! give warning and ignore this line
			write(errstr,*) "could not read complex value from ",iline,"-th line '",trim(line),"', iostat = ",ios
			call add(errmsg,1,trim(errstr),myname)
		endif
!
		if(dont_read_all) then
			! terminate loop if required number of data values is reached
			if(idata==ndata) exit
		endif
	enddo ! while(ios_line==0)
!
	close(lu)
!
	if(dont_read_all) then
		if(idata/=ndata) then
			! raise error
			write(errstr,*) "could not read the required number of data values, instead ",idata," values were read"
			call add(errmsg,2,trim(errstr),myname)
			! reallocate
			p => reallocate(p,idata)
		else
			call add(errmsg,0,'all required data values were found',myname)
		endif
	else ! dont_read_all
		write(errstr,*) "there were ",idata," data values found"
		call add(errmsg,0,trim(errstr),myname)
		! reallocate
		p => reallocate(p,idata)
	endif ! dont_read_all
!
	if(idata==0) then
		! fix problem of relallocating size 0 arrays (they have size 0, but pointers are associated)
		if(associated(p)) deallocate(p)
		p => null()
	endif
	end function readComplexVectorAsciiDataIO
!------------------------------------------------------------------------
!> \brief read the same number (>=0) of real values from each line of ascii file
!! \param filename ascii file name
!! \param lu file unit
!! \param nrow optional number of rows to be read from file
!! \param nskip optional number of lines to be skipped before reading the data
!! \param errmsg error message
!! \return error message
!
	function readRealMatrixAsciiDataIO(filename,lu,p,ncol,nrow,nskip) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	type (error_message) :: errmsg
	integer :: ncol
	integer, optional :: nrow,nskip
	! outgoing
	real, dimension(:,:), pointer :: p
	! local
	integer, parameter :: nchunk = 100 !< size by which data array is dynamically reallocated (if read all)
	real, dimension(:), allocatable :: r
	integer :: ios,ios_line,irow,iline
	logical :: dont_read_all
	character(len=400) :: errstr,line
	character(len=25) :: myname = 'readRealMatrixAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
	p => null()
!
	! check on ncol
	if(ncol.le.0) then
		write(errstr,*) 'desired number of columns ( =',ncol,') is less than 1, must be at least 1'
		call add(errmsg,2,trim(errstr),myname)
		return
	endif
!
	! define the required number of lines=rows to read
	if(present(nrow)) then
		if(nrow<1) then
			write(errstr,*) 'desired number of rows ( =',nrow,') is less than 1'
			call add(errmsg,2,trim(errstr),myname)
			return
		else
			dont_read_all = .true.
			write(errstr,*) 'trying to read ',nrow,' rows from file'
			call add(errmsg,0,trim(errstr),myname)
		endif
	else
		dont_read_all = .false.
		call add(errmsg,0,'trying to read all rows from the file',myname)
	endif
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	! if required, skip the first nskip lines of the file
	if(present(nskip)) then
		write(errstr,*) 'the first ',nskip,' lines of the file will be skipped'
		call add(errmsg,0,trim(errstr),myname)
		do iline = 1,nskip
			! read new line from file
			read(lu,"(a400)",iostat=ios_line) line
			if(ios_line/=0) then
				! write maximum number of lines read from file to error message and return
				write(errstr,*) "Only ",iline-1," lines could be read from file. No data was read."
				call add(errmsg,1,trim(errstr),myname)
				return
			endif ! ios_line/=0
				
		enddo ! iline
	endif ! present(nskip)
!
	allocate(r(ncol))
	if(dont_read_all) then
		allocate(p(nrow,ncol))
	else
		allocate(p(nchunk,ncol))
	endif
	irow = 0
	iline = 0
	ios_line=0
!
	do while(ios_line==0)
		! counter on lines (for error messages)
		iline = iline+1
		! read new line from file
		read(lu,"(a400)",iostat=ios_line) line
		if(ios_line/=0) then
			! write maximum number of lines read from file to error message and terminate loop
			write(errstr,*) iline-1," lines could be read from file"
			call add(errmsg,0,trim(errstr),myname)
			exit
		endif ! ios_line/=0
!
		! try to read complex number from line string
		read(line,*,iostat=ios) r
		if(ios==0) then
			! reallocate if necessary
			if(irow==size(p,1)) p => reallocate(p,size(p,1)+nchunk,ncol)
			! add data value to array
			irow = irow+1
			p(irow,:) = r
		else
			! give warning and ignore this line
			write(errstr,*) "could not read ",ncol," complex values from ",iline,"-th line '",trim(line),"', iostat = ",ios
			call add(errmsg,1,trim(errstr),myname)
		endif
!
		if(dont_read_all) then
			! terminate loop if required number of data values is reached
			if(irow==nrow) exit
		endif
	enddo ! while(ios_line==0)
!
	close(lu)
!
	if(dont_read_all) then
		if(irow/=nrow) then
			! raise error
			write(errstr,*) "could not read the required number of rows, instead ",irow," rows were read"
			call add(errmsg,2,trim(errstr),myname)
			! reallocate
			p => reallocate(p,irow,ncol)
		else
			call add(errmsg,0,'all required data values were found',myname)
		endif
	else ! dont_read_all
		write(errstr,*) "there were ",irow," rows of values found"
		call add(errmsg,0,trim(errstr),myname)
		! reallocate
		p => reallocate(p,irow,ncol)
	endif ! dont_read_all
!
	if(irow==0) then
		! fix problem of relallocating size 0 arrays (they have size 0, but pointers are associated)
		if(associated(p)) deallocate(p)
		p => null()
	endif
	if(allocated(r)) deallocate(r)
	end function readRealMatrixAsciiDataIO
!------------------------------------------------------------------------
!> \brief read the same number (>=0) of complex values from each line of ascii file
!! \param filename ascii file name
!! \param lu file unit
!! \param nrow optional number of rows to be read from file
!! \param nskip optional number of lines to be skipped before reading the data
!! \param errmsg error message
!! \return error message
!
	function readComplexMatrixAsciiDataIO(filename,lu,p,ncol,nrow,nskip) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	type (error_message) :: errmsg
	integer :: ncol
	integer, optional :: nrow,nskip
	! outgoing
	complex, dimension(:,:), pointer :: p
	! local
	integer, parameter :: nchunk = 100 !< size by which data array is dynamically reallocated (if read all)
	complex, dimension(:), allocatable :: z
	integer :: ios,ios_line,irow,iline
	logical :: dont_read_all
	character(len=400) :: errstr,line
	character(len=28) :: myname = 'readComplexMatrixAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
	p => null()
!
	! check on ncol
	if(ncol.le.0) then
		write(errstr,*) 'desired number of columns ( =',ncol,') is less than 1, must be at least 1'
		call add(errmsg,2,trim(errstr),myname)
		return
	endif
!
	! define the required number of lines=rows to read
	if(present(nrow)) then
		if(nrow<1) then
			write(errstr,*) 'desired number of rows ( =',nrow,') is less than 1'
			call add(errmsg,2,trim(errstr),myname)
			return
		else
			dont_read_all = .true.
			write(errstr,*) 'trying to read ',nrow,' rows from file'
			call add(errmsg,0,trim(errstr),myname)
		endif
	else
		dont_read_all = .false.
		call add(errmsg,0,'trying to read all rows from the file',myname)
	endif
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	! if required, skip the first nskip lines of the file
	if(present(nskip)) then
		write(errstr,*) 'the first ',nskip,' lines of the file will be skipped'
		call add(errmsg,0,trim(errstr),myname)
		do iline = 1,nskip
			! read new line from file
			read(lu,"(a400)",iostat=ios_line) line
			if(ios_line/=0) then
				! write maximum number of lines read from file to error message and return
				write(errstr,*) "Only ",iline-1," lines could be read from file. No data was read."
				call add(errmsg,1,trim(errstr),myname)
				return
			endif ! ios_line/=0
				
		enddo ! iline
	endif ! present(nskip)
!
	allocate(z(ncol))
	if(dont_read_all) then
		allocate(p(nrow,ncol))
	else
		allocate(p(nchunk,ncol))
	endif
	irow = 0
	iline = 0
	ios_line=0
!
	do while(ios_line==0)
		! counter on lines (for error messages)
		iline = iline+1
		! read new line from file
		read(lu,"(a400)",iostat=ios_line) line
		if(ios_line/=0) then
			! write maximum number of lines read from file to error message and terminate loop
			write(errstr,*) iline-1," lines could be read from file"
			call add(errmsg,0,trim(errstr),myname)
			exit
		endif ! ios_line/=0
!
		! try to read complex number from line string
		read(line,*,iostat=ios) z
		if(ios==0) then
			! reallocate if necessary
			if(irow==size(p,1)) p => reallocate(p,size(p,1)+nchunk,ncol)
			! add data value to array
			irow = irow+1
			p(irow,:) = z
		else
			! give warning and ignore this line
			write(errstr,*) "could not read ",ncol," complex values from ",iline,"-th line '",trim(line),"', iostat = ",ios
			call add(errmsg,1,trim(errstr),myname)
		endif
!
		if(dont_read_all) then
			! terminate loop if required number of data values is reached
			if(irow==nrow) exit
		endif
	enddo ! while(ios_line==0)
!
	close(lu)
!
	if(dont_read_all) then
		if(irow/=nrow) then
			! raise error
			write(errstr,*) "could not read the required number of rows, instead ",irow," rows were read"
			call add(errmsg,2,trim(errstr),myname)
			! reallocate
			p => reallocate(p,irow,ncol)
		else
			call add(errmsg,0,'all required data values were found',myname)
		endif
	else ! dont_read_all
		write(errstr,*) "there were ",irow," rows of values found"
		call add(errmsg,0,trim(errstr),myname)
		! reallocate
		p => reallocate(p,irow,ncol)
	endif ! dont_read_all
!
	if(irow==0) then
		! fix problem of relallocating size 0 arrays (they have size 0, but pointers are associated)
		if(associated(p)) deallocate(p)
		p => null()
	endif
	if(allocated(z)) deallocate(z)
	end function readComplexMatrixAsciiDataIO
!------------------------------------------------------------------------
!> \brief write complex values to ascii file (one value per line)
!! \param filename ascii file name
!! \param lu file unit
!! \param val complex data vector
!! \return error message
!
! TODO: * add optional array of character strings, which contains arbitrary header lines
!         which are written first before the data values are written
!       * add sophisticated overwrite handling
!
	function writeComplexVectorAsciiDataIO(filename,lu,val) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	complex, dimension(:) :: val
	! outgoing
	type (error_message) :: errmsg
	! local
	!complex :: z
	integer :: ios,idata,ndata!,ios_line,idata,iline
	!logical :: dont_read_all
	character(len=400) :: errstr!,line
	character(len=29) :: myname = 'writeComplexVectorAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='UNKNOWN',action='WRITE',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	ndata = size(val)
	write(errstr,*) "incoming data values = ",ndata
	call new(errmsg,0,trim(errstr),myname)
!
	do idata = 1,size(val)
		write(lu,*,iostat=ios) val(idata)
		if(ios/=0) then
			write(errstr,*) "could not write ",idata,"'th value; iostat = ",ios
			call add(errmsg,2,trim(errstr),myname)
			exit
		endif
	end do ! idata
!
	close(lu)
	end function writeComplexVectorAsciiDataIO
!------------------------------------------------------------------------
!> \brief write complex values to ascii file (n values per line)
!! \param filename ascii file name
!! \param lu file unit
!! \param val complex data matrix
!! \return error message
!
! TODO: * add optional array of character strings, which contains arbitrary header lines
!         which are written first before the data values are written
!       * add sophisticated overwrite handling
!
	function writeComplexMatrixAsciiDataIO(filename,lu,val) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	complex, dimension(:,:) :: val
	! outgoing
	type (error_message) :: errmsg
	! local
	integer :: ios,idata,ndata,icol,ncol!,ios_line,idata,iline
	!logical :: dont_read_all
	character(len=400) :: errstr!,line
	character(len=29) :: myname = 'writeComplexMatrixAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='UNKNOWN',action='WRITE',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	ndata = size(val,1)
	ncol = size(val,2)
	write(errstr,*) "incoming data values = ",ndata
	call new(errmsg,0,trim(errstr),myname)
	write(errstr,*) "incoming number of traces = ",ncol
	call new(errmsg,0,trim(errstr),myname)
!
	do idata = 1,ndata
		write(lu,*,iostat=ios) (val(idata,icol),icol = 1,ncol)
		if(ios /= 0) then
			write(errstr,*) "could not write ",idata,"'th values; iostat = ",ios
			call add(errmsg,2,trim(errstr),myname)
			exit
		endif
	end do ! idata
!
	close(lu)
	end function writeComplexMatrixAsciiDataIO
!------------------------------------------------------------------------
!> \brief write real values to ascii file (one value per line)
!! \param filename ascii file name
!! \param lu file unit
!! \param val real data vector
!! \return error message
!
! TODO: * add optional array of character strings, which contains arbitrary header lines
!         which are written first before the data values are written
!       * add sophisticated overwrite handling
!
	function writeRealVectorAsciiDataIO(filename,lu,val) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	real, dimension(:) :: val
	! outgoing
	type (error_message) :: errmsg
	! local
	integer :: ios,idata,ndata
	character(len=400) :: errstr
	character(len=26) :: myname = 'writeRealVectorAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='UNKNOWN',action='WRITE',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	ndata = size(val)
	write(errstr,*) "incoming data values = ",ndata
	call new(errmsg,0,trim(errstr),myname)
!
	do idata = 1,size(val)
		write(lu,*,iostat=ios) val(idata)
		if(ios/=0) then
			write(errstr,*) "could not write ",idata,"'th value; iostat = ",ios
			call add(errmsg,2,trim(errstr),myname)
			exit
		endif
	end do ! idata
!
	close(lu)
	end function writeRealVectorAsciiDataIO
!------------------------------------------------------------------------
!> \brief write real values to ascii file (n values per line)
!! \param filename ascii file name
!! \param lu file unit
!! \param val real data matrix
!! \return error message
!
! TODO: * add optional array of character strings, which contains arbitrary header lines
!         which are written first before the data values are written
!       * add sophisticated overwrite handling
!
	function writeRealMatrixAsciiDataIO(filename,lu,val) result(errmsg)
	! incoming
	character(len=*) :: filename
	integer :: lu
	real, dimension(:,:) :: val
	! outgoing
	type (error_message) :: errmsg
	! local
	integer :: ios,idata,ndata,icol,ncol!,ios_line,idata,iline
	!logical :: dont_read_all
	character(len=400) :: errstr!,line
	character(len=29) :: myname = 'writeRealMatrixAsciiDataIO'
!
	call new(errmsg,0,"incoming filename: '"//trim(filename)//"'",myname)
!
	! open file
	open(unit=lu,file=trim(filename),form='FORMATTED',status='UNKNOWN',action='WRITE',iostat=ios)
	if(ios/=0) then
		write(errstr,*) "could not open file, iostat = ",ios
		call add(errmsg,2,trim(errstr),myname)
		close(lu)
		return
	endif
!
	ndata = size(val,1)
	ncol = size(val,2)
	write(errstr,*) "incoming data values = ",ndata
	call new(errmsg,0,trim(errstr),myname)
	write(errstr,*) "incoming number of traces = ",ncol
	call new(errmsg,0,trim(errstr),myname)
!
	do idata = 1,ndata
		write(lu,*,iostat=ios) (val(idata,icol),icol = 1,ncol)
		if(ios /= 0) then
			write(errstr,*) "could not write ",idata,"'th values; iostat = ",ios
			call add(errmsg,2,trim(errstr),myname)
			exit
		endif
	end do ! idata
!
	close(lu)
	end function writeRealMatrixAsciiDataIO
!
end module asciiDataIO
