!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!------------------------------------------------------------------------
!> \brief Module defining vector pointers
!
 module vectorPointer
	implicit none
	interface associateVectorPointer
		module procedure associateIntegerVectorPointer
		module procedure associateRealVectorPointer
		module procedure associateDoubleVectorPointer
		module procedure associateComplexVectorPointer
		module procedure associateDoubleComplexVectorPointer
	end interface
	interface fillVectorPointer
		module procedure fillIntegerVectorPointer
		module procedure fillRealVectorPointer
		module procedure fillDoubleVectorPointer
		module procedure fillComplexVectorPointer
		module procedure fillDoubleComplexVectorPointer
	end interface
	interface getVectorPointer
		module procedure getIntegerVectorPointer
		module procedure getRealVectorPointer
		module procedure getDoubleVectorPointer
		module procedure getComplexVectorPointer
		module procedure getDoubleComplexVectorPointer
	end interface
	interface allocateVectorPointer
		module procedure allocateIntegerVectorPointer
		module procedure allocateRealVectorPointer
		module procedure allocateDoubleVectorPointer
		module procedure allocateComplexVectorPointer
		module procedure allocateDoubleComplexVectorPointer
	end interface
	interface dealloc
		module procedure deallocIntegerVectorPointer
		module procedure deallocRealVectorPointer
		module procedure deallocDoubleVectorPointer
		module procedure deallocComplexVectorPointer
		module procedure deallocDoubleComplexVectorPointer
	end interface
	interface nullifyVectorPointer
		module procedure nullifyIntegerVectorPointer
	end interface
	interface extendArrayVectorPointer
		module procedure extendArrayIntegerVectorPointer
	end interface
!!
	type integer_vector_pointer
		private
		integer, dimension(:), pointer :: p => null()
	end type
	type real_vector_pointer
		private
		real, dimension(:), pointer :: p => null()
	end type
	type double_vector_pointer
		private
		double precision, dimension(:), pointer :: p => null()
	end type
	type complex_vector_pointer
		private
		complex, dimension(:), pointer :: p => null()
	end type
	type double_complex_vector_pointer
		private
		double complex, dimension(:), pointer :: p => null()
	end type
!
 contains
!------------------------------------------------------------------------
!> \brief Allocate space for a integer vector pointer
!
	subroutine allocateIntegerVectorPointer(this,n)
	type (integer_vector_pointer) :: this
	integer :: n
	allocate(this%p(n))
	end subroutine allocateIntegerVectorPointer
!------------------------------------------------------------------------
!> \brief Allocate space for a real vector pointer
!
	subroutine allocateRealVectorPointer(this,n)
	type (real_vector_pointer) :: this
	integer :: n
	allocate(this%p(n))
	end subroutine allocateRealVectorPointer
!------------------------------------------------------------------------
!> \brief Allocate space for a double vector pointer
!
	subroutine allocateDoubleVectorPointer(this,n)
	type (double_vector_pointer) :: this
	integer :: n
	allocate(this%p(n))
	end subroutine allocateDoubleVectorPointer
!------------------------------------------------------------------------
!> \brief Allocate space for a complex vector pointer
!
	subroutine allocateComplexVectorPointer(this,n)
	type (complex_vector_pointer) :: this
	integer :: n
	allocate(this%p(n))
	end subroutine allocateComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Allocate space for a double complex vector pointer
!
	subroutine allocateDoubleComplexVectorPointer(this,n)
	type (double_complex_vector_pointer) :: this
	integer :: n
	allocate(this%p(n))
	end subroutine allocateDoubleComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Fill a integer vector pointer
!
	subroutine fillIntegerVectorPointer(this,d,ns)
	type (integer_vector_pointer) :: this
	integer, dimension(:) :: d
	integer :: ns
	this%p(ns:ns+size(d)-1) = d
	end subroutine fillIntegerVectorPointer
!------------------------------------------------------------------------
!> \brief Fill a real vector pointer
!
	subroutine fillRealVectorPointer(this,d,ns)
	type (real_vector_pointer) :: this
	real, dimension(:) :: d
	integer :: ns
	this%p(ns:ns+size(d)-1) = d
	end subroutine fillRealVectorPointer
!------------------------------------------------------------------------
!> \brief Fill a double vector pointer
!
	subroutine fillDoubleVectorPointer(this,d,ns)
	type (double_vector_pointer) :: this
	double precision, dimension(:) :: d
	integer :: ns
	this%p(ns:ns+size(d)-1) = d
	end subroutine fillDoubleVectorPointer
!------------------------------------------------------------------------
!> \brief Fill a complex vector pointer
!
	subroutine fillComplexVectorPointer(this,d,ns)
	type (complex_vector_pointer) :: this
	complex, dimension(:) :: d
	integer :: ns
	this%p(ns:ns+size(d)-1) = d
	end subroutine fillComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Fill a double complex vector pointer
!
	subroutine fillDoubleComplexVectorPointer(this,d,ns)
	type (double_complex_vector_pointer) :: this
	double complex, dimension(:) :: d
	integer :: ns
	this%p(ns:ns+size(d)-1) = d
	end subroutine fillDoubleComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Associate a integer vector pointer
!
	subroutine associateIntegerVectorPointer(this,id)
	type (integer_vector_pointer) :: this
	integer, dimension(:), target :: id
	this%p => id
	end subroutine associateIntegerVectorPointer
!------------------------------------------------------------------------
!> \brief Associate a real vector pointer
!
	subroutine associateRealVectorPointer(this,rd)
	type (real_vector_pointer) :: this
	real, dimension(:), target :: rd
	this%p => rd
	end subroutine associateRealVectorPointer
!------------------------------------------------------------------------
!> \brief Associate a double vector pointer
!
	subroutine associateDoubleVectorPointer(this,dd)
	type (double_vector_pointer) :: this
	double precision, dimension(:), target :: dd
	this%p => dd
	end subroutine associateDoubleVectorPointer
!------------------------------------------------------------------------
!> \brief Associate a complex vector pointer
!
	subroutine associateComplexVectorPointer(this,cd)
	type (complex_vector_pointer) :: this
	complex, dimension(:), target :: cd
	this%p => cd
	end subroutine associateComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Associate a dopuble complex vector pointer
!
	subroutine associateDoubleComplexVectorPointer(this,dcd)
	type (double_complex_vector_pointer) :: this
	double complex, dimension(:), target :: dcd
	this%p => dcd
	end subroutine associateDoubleComplexVectorPointer
!---------------------------------------------------------------------
!> \brief Deallocate a integer vector pointer
!
	subroutine deallocIntegerVectorPointer(this)
	type (integer_vector_pointer) :: this
	if (associated(this%p)) deallocate(this%p)
	end subroutine deallocIntegerVectorPointer
!---------------------------------------------------------------------
!> \brief Deallocate a real vector pointer
!
	subroutine deallocRealVectorPointer(this)
	type (real_vector_pointer) :: this
	if (associated(this%p)) deallocate(this%p)
	end subroutine deallocRealVectorPointer
!---------------------------------------------------------------------
!> \brief Deallocate a double vector pointer
!
	subroutine deallocDoubleVectorPointer(this)
	type (double_vector_pointer) :: this
	if (associated(this%p)) deallocate(this%p)
	end subroutine deallocDoubleVectorPointer
!---------------------------------------------------------------------
!> \brief Deallocate a complex vector pointer
!
	subroutine deallocComplexVectorPointer(this)
	type (complex_vector_pointer) :: this
	if (associated(this%p)) deallocate(this%p)
	end subroutine deallocComplexVectorPointer
!---------------------------------------------------------------------
!> \brief Deallocate a double complex vector pointer
!
	subroutine deallocDoubleComplexVectorPointer(this)
	type (double_complex_vector_pointer) :: this
	if (associated(this%p)) deallocate(this%p)
	end subroutine deallocDoubleComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Get a integer vector pointer
!
	function getIntegerVectorPointer(this) result(p)
	type (integer_vector_pointer) :: this
	integer, dimension(:), pointer :: p
	p => this%p
	end function getIntegerVectorPointer
!------------------------------------------------------------------------
!> \brief Get a real vector pointer
!
	function getRealVectorPointer(this) result(p)
	type (real_vector_pointer) :: this
	real, dimension(:), pointer :: p
	p => this%p
	end function getRealVectorPointer
!------------------------------------------------------------------------
!> \brief Get a double vector pointer
!
	function getDoubleVectorPointer(this) result(p)
	type (double_vector_pointer) :: this
	double precision, dimension(:), pointer :: p
	p => this%p
	end function getDoubleVectorPointer
!------------------------------------------------------------------------
!> \brief Get a complex vector pointer
!
	function getComplexVectorPointer(this) result(p)
	type (complex_vector_pointer) :: this
	complex, dimension(:), pointer :: p
	p => this%p
	end function getComplexVectorPointer
!------------------------------------------------------------------------
!> \brief Get a double complex vector pointer
!
	function getDoubleComplexVectorPointer(this) result(p)
	type (double_complex_vector_pointer) :: this
	double complex, dimension(:), pointer :: p
	p => this%p
	end function getDoubleComplexVectorPointer
!---------------------------------------------------------------------
!> \brief Nullify an integer vector pointer
!
	subroutine nullifyIntegerVectorPointer(this)
	type (integer_vector_pointer) :: this
	if (associated(this%p)) nullify(this%p)
	end subroutine nullifyIntegerVectorPointer
!-----------------------------------------------------------------------
!> \brief Extend an array of vector pointers
!
	function extendArrayIntegerVectorPointer(array,n) result(newarray)
	type (integer_vector_pointer), dimension(:), pointer :: array
	type (integer_vector_pointer), dimension(:), pointer :: newarray
	integer :: n,nold,i
!
	allocate(newarray(n))
	if (.not. associated(array)) return
	nold = min(size(array),n)
	newarray(1:nold) = array(1:nold)
	do i = 1,nold
		call nullifyIntegerVectorPointer(array(i))
	enddo
	do i = nold+1,size(array)
		call dealloc(array(i))
	enddo
	deallocate(array)
	end function extendArrayIntegerVectorPointer
!
 end module vectorPointer
