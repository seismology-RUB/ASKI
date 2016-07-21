!--------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 1.0.
!
!   Gemini II and ASKI version 1.0 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 1.0 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!-------------------------------------------------------------------------
! module for a variable of flexible type
! can take integer, real, double precision, complex, double complex, string
!-------------------------------------------------------------------------
 module flexibleType
	use kindDefinitions
	use primitiveTypeEncoding
	use simpleString
	implicit none
	interface assignment (=);
		module procedure createIntegerFlexibleType
		module procedure createRealFlexibleType
		module procedure createDoubleFlexibleType
		module procedure createComplexFlexibleType
		module procedure createDoubleComplexFlexibleType
		module procedure createCharFlexibleType
		module procedure valueIntegerFlexibleType
		module procedure valueRealFlexibleType
		module procedure valueDoubleFlexibleType
		module procedure valueComplexFlexibleType
		module procedure valueDoubleComplexFlexibleType
		module procedure valueCharFlexibleType
	end interface	
	interface print; module procedure printFlexibleType; end interface
	type flexible
		private
		integer :: datatype
		integer :: val_int
		real :: val_real
		double precision :: val_double
		complex :: val_complex
		double complex :: val_double_complex
		type (simple_string):: val_string
	end type
!
 contains
!----------------------------------------------------------------------
!  create a flexibleType object with integer content
!
	subroutine createIntegerFlexibleType(this,value)
	type (flexible), intent(out) :: this
	integer, intent(in) :: value
	this%val_int = value
	this%datatype = T_INTEGER
	end subroutine createIntegerFlexibleType
!----------------------------------------------------------------------
!  create a flexibleType object with real content
!
	subroutine createRealFlexibleType(this,value)
	type (flexible), intent(out) :: this
	real, intent(in) :: value
	this%val_real = value
	this%datatype = T_REAL
	end subroutine createRealFlexibleType
!----------------------------------------------------------------------
!  create a flexibleType object with double content
!
	subroutine createDoubleFlexibleType(this,value)
	type (flexible), intent(out) :: this
	double precision, intent(in) :: value
	this%val_double = value
	this%datatype = T_DOUBLE
	end subroutine createDoubleFlexibleType
!----------------------------------------------------------------------
!  create a flexibleType object with complex content
!
	subroutine createComplexFlexibleType(this,value)
	type (flexible), intent(out) :: this
	complex, intent(in) :: value
	this%val_complex = value
	this%datatype = T_COMPLEX
	end subroutine createComplexFlexibleType
!----------------------------------------------------------------------
!  create a flexibleType object with double complex content
!
	subroutine createDoubleComplexFlexibleType(this,value)
	type (flexible), intent(out) :: this
	double complex, intent(in) :: value
	this%val_double_complex = value
	this%datatype = T_DOUBLE_COMPLEX
	end subroutine createDoubleComplexFlexibleType
!----------------------------------------------------------------------
!  create a flexibleType object with string content
!
	subroutine createCharFlexibleType(this,value)
	type (flexible), intent(out) :: this
	character (len=*), intent(in) :: value
	call createSimpleString(this%val_string,value)
	this%datatype = T_CHAR
	end subroutine createCharFlexibleType
!----------------------------------------------------------------------
!  deallocate flexible type object
!
	subroutine deallocFlexibleType(this)
	type (flexible) :: this
	if (this%datatype == T_CHAR) then
		call dealloc(this%val_string)
	endif
	end subroutine deallocFlexibleType
!--------------------------------------------------------------------
!  get type inside flexible object
!
	function typeFlexibleType(this) result(typ)
	type (flexible) :: this
	integer :: typ
	typ = this%datatype
	end function typeFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueIntegerFlexibleType(val,this)
	type (flexible), intent(in) :: this
	integer, intent(out) :: val
	val = this%val_int
	end subroutine valueIntegerFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueRealFlexibleType(val,this)
	type (flexible), intent(in) :: this
	real, intent(out) :: val
	val = this%val_real
	end subroutine valueRealFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueDoubleFlexibleType(val,this)
	type (flexible), intent(in) :: this
	double precision, intent(out) :: val
	val = this%val_double
	end subroutine valueDoubleFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueComplexFlexibleType(val,this)
	type (flexible), intent(in) :: this
	complex, intent(out) :: val
	val = this%val_complex
	end subroutine valueComplexFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueDoubleComplexFlexibleType(val,this)
	type (flexible), intent(in) :: this
	double complex, intent(out) :: val
	val = this%val_double_complex
	end subroutine valueDoubleComplexFlexibleType
!---------------------------------------------------------------------
!  get back value from flexible type object
!
	subroutine valueCharFlexibleType(val,this)
	type (flexible), intent(in) :: this
	character (len=*), intent(out) :: val
	call convertToCharSimpleString(this%val_string,val)
	end subroutine valueCharFlexibleType
!---------------------------------------------------------------------
!  write a flexible type object to a stream access file
!  returns number of bytes written
!
	subroutine writeSAFlexibleType(this,lu,filepos,nbytes)
	type (flexible), intent(in) :: this
	integer :: lu,nbytes
	integer :: clen
	integer (longint) :: filepos
	character (len=132) :: string
	select case (this%datatype)
	case (T_INTEGER); write(lu,pos=filepos) this%datatype,this%val_int; nbytes = 2*kind(1)
	case (T_REAL);    write(lu,pos=filepos) this%datatype,this%val_real; nbytes = kind(1)+kind(1.0)
	case (T_DOUBLE);  write(lu,pos=filepos) this%datatype,this%val_double; nbytes = kind(1)+kind(1.d0)
	case (T_COMPLEX); write(lu,pos=filepos) this%datatype,this%val_complex; nbytes = kind(1)+2*kind(1.0)
	case (T_DOUBLE_COMPLEX)
		write(lu,pos=filepos) this%datatype,this%val_double_complex; nbytes = kind(1)+2*kind(1.d0)
	case (T_CHAR)
		call valueCharFlexibleType(string,this)
		clen = len_trim(string)
		write(lu,pos=filepos) this%datatype,clen,trim(string)
		nbytes = 2*kind(1)+clen
	end select
	end subroutine writeSAFlexibleType
!---------------------------------------------------------------------
!  read a flexible type object from stream access file
!  returns number of bytes read in
!
	subroutine readSAFlexibleType(this,lu,filepos,nbytes)
	type (flexible) :: this
	integer :: lu,nbytes
	integer (longint) :: filepos
	integer :: clen
	character (len=132) :: string
	read(lu,pos = filepos) this%datatype
	select case (this%datatype)
	case (T_INTEGER); read(lu,pos=filepos+kind(1)) this%val_int; nbytes = kind(1)
	case (T_REAL);    read(lu,pos=filepos+kind(1)) this%val_real; nbytes = kind(1.0)
	case (T_DOUBLE);  read(lu,pos=filepos+kind(1)) this%val_double; nbytes = kind(1.d0)
	case (T_COMPLEX); read(lu,pos=filepos+kind(1)) this%val_complex; nbytes = 2*kind(1.0)
	case (T_DOUBLE_COMPLEX); read(lu,pos=filepos+kind(1)) this%val_double_complex; nbytes = 2*kind(1.d0)
	case (T_CHAR)
		read(lu,pos=filepos+kind(1)) clen,string(1:clen)
		nbytes = kind(1)+clen
		call createSimpleString(this%val_string,string(1:clen))
	end select
	nbytes = nbytes+kind(1)
	end subroutine readSAFlexibleType
!---------------------------------------------------------------------
!  print value of flexible object
!
	subroutine printFlexibleType(this)
	type (flexible), intent(in) :: this
	character (len=132) :: string
	select case (this%datatype)
	case (T_INTEGER); print *,this%val_int
	case (T_REAL);    print *,this%val_real
	case (T_DOUBLE);  print *,this%val_double
	case (T_COMPLEX); print *,this%val_complex
	case (T_DOUBLE_COMPLEX); print *,this%val_double_complex
	case (T_CHAR);    call convertToCharSimpleString(this%val_string,string); print *,trim(string)
	end select
	end subroutine printFlexibleType
!
 end module flexibleType
