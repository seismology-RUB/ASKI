!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------
!> \brief handle linear systems and solve them with LAPACK
!! 
!! \details Module that handles linear systems of equations which
!!  can be solved with serial LAPACK algorithms. 
!!  Several solvers for different types of matrices
!!  are supported.
!! 
!! \author Florian Schumacher
!! \date 2015
!---------------------------------------------------------------
module serialLinearSystem
!
	use errorMessage
!
	implicit none
!
!!! THIS DOES NOT WORK, BECAUSE THEN, MODULE linearSystem WHICH USES THIS MODULE serialLinearSystem CANNOT USE NAMES "solveLeastSquaresLinearSystem" etc
!!$!> \brief overloading with name consistent with module linearSystem
!!$	interface createLinearSystem
!!$		module procedure createSerialLinearSystem
!!$	end interface createLinearSystem
!!$!> \brief overloading with name consistent with module linearSystem
!!$	interface solveLeastSquaresLinearSystem
!!$		module procedure solveLeastSquaresSerialLinearSystem
!!$	end interface solveLeastSquaresLinearSystem
	interface createSerialLinearSystem
		module procedure createUnpackedSerialLinearSystem
		module procedure createPackedSerialLinearSystem
	end interface createSerialLinearSystem
!> \brief overloading deallocate routine with shorter name
	interface dealloc
		module procedure deallocSerialLinearSystem
	end interface dealloc
	interface operator (.A.); module procedure getMatrixSerialLinearSystem; end interface
	interface operator (.AP.); module procedure getPackedMatrixSerialLinearSystem; end interface
!
!> \brief all specifications of the linear system
	type serial_linear_system
		private
		integer :: nrow = 0 !< number of rows of linear system
		integer :: ncol = 0 !< number of columns of linear system
		integer :: nrhs = 0 !< number of right-hand-sides of linear system
		real, dimension(:,:), pointer :: A => null() !< nrow-by-ncol system matrix
		character(len=1) :: transpose = '' !< either 'N' (default) or 'T', indicating whether A or A^T is involved in solving the linear system
		real, dimension(:,:), pointer :: b => null() !< array of size(nrow,nrhs) ('N') or size(ncol,nrhs) ('T'), right-hand-side(s) of system
		! FOR PACKED STORAGE OF SYMMETRIC LINEAR SYSTEMS:
		character(len=1) :: UPLO = '' !< indicating if A contains columnwisely the 'U'pper or 'L'ower triangular matrix
		! if UPLO = 'U', A_packed(i + (j-1)*j/2) = A(i,j) for 1<=i<=j
		! if UPLO = 'L', A_packed(i + (j-1)*(2nrow-j)/2) = A(i,j) for j<=i<=nrow
		real, dimension(:), pointer :: A_packed => null() !< nrow*(nrow+1)/2-sized array for packed storage of symmetric matrices
	end type serial_linear_system
!
contains
!
!------------------------------------------------------------------------
!> \brief create serial linear system
!! \details Set the fields in this serial_linear_system object. System matrix A
!!  and right-hand-side(s) b are just pointers, i.e. there is no space allocated here.
!!  Take care of the memory for A and b outside this module.
!! \param this serial_linear_system object
!! \param nrow number of rows of linear system
!! \param ncol number of columns of linear system
!! \param nrhs number of right-hand-sides of linear system
!! \param A pointer to system matrix
!! \param b pointer to matrix of right-hand-sides
!
	subroutine createUnpackedSerialLinearSystem(this,nrow,ncol,nrhs,A,b,transpose)
	type (serial_linear_system) :: this
	integer :: nrow,ncol,nrhs
	character(len=*), optional :: transpose
	real, dimension(:,:), pointer :: A
	real, dimension(:,:), pointer :: b
	call deallocSerialLinearSystem(this)
	this%nrow=nrow; this%ncol=ncol; this%nrhs=nrhs
	this%A => A
	this%A_packed => null()
	this%b => b
	if(present(transpose)) then
		this%transpose = transpose
	else
		this%transpose = 'N'
	endif
	end subroutine createUnpackedSerialLinearSystem
!------------------------------------------------------------------------
!> \brief create packed symmetric serial linear system
!! \details Set the fields in this serial_linear_system object. Packed system matrix A_packed
!!  and right-hand-side(s) b are just pointers, i.e. there is no space allocated here.
!!  Take care of the memory for A_packed and b outside this module.
!! \param this serial_linear_system object
!! \param nrow number of rows of linear system
!! \param ncol number of columns of linear system
!! \param nrhs number of right-hand-sides of linear system
!! \param A pointer to system matrix
!! \param b pointer to matrix of right-hand-sides
!
	subroutine createPackedSerialLinearSystem(this,nrow,ncol,nrhs,A_packed,UPLO,b)
	type (serial_linear_system) :: this
	integer :: nrow,ncol,nrhs
	real, dimension(:), pointer :: A_packed
	character(len=1) :: UPLO
	real, dimension(:,:), pointer :: b
	call deallocSerialLinearSystem(this)
	this%nrow=nrow; this%ncol=ncol; this%nrhs=nrhs
	this%A => null()
	this%A_packed => A_packed
	this%UPLO = UPLO
	this%b => b
	end subroutine createPackedSerialLinearSystem
!-----------------------------------------------------------------
!> \brief nullify pointers of system matrix (A) and right-hand-side (b)
!! \param this serial_linear_system object
!
	subroutine deallocSerialLinearSystem(this)
	type (serial_linear_system) :: this
	if(associated(this%A)) this%A => null()
	if(associated(this%A_packed)) this%A_packed => null()
	this%UPLO = ''
	if(associated(this%b)) this%b => null()
	this%nrow = 0; this%ncol = 0; this%nrhs = 0
	end subroutine deallocSerialLinearSystem
!-----------------------------------------------------------------
!> \brief solve linear system by least squares
!! \details Calling LAPACK routine "SGELS" which does a QR decomposition of the system matrix, i.e. approximates
!!  the solution by least squares (if nrow >= ncol), or computes the minimum norm solution (if nrow < ncol).
!! \param this serial_linear_system object
!! \param x ncol-by-nrhs solution array (first ncol columns of argument "B" of LAPACK routine "SGELS")
!! \param overwrite_A logical to indicate whether system matrix A should be overwritten or not
!! \param errmsg error_message object
!! \return error message
!! \pre You need LAPACK libraries. Routine "SGELS" is called here.
!
	function solveLeastSquaresSerialLinearSystem(this,x,res,overwrite_A) result(errmsg)
	type (serial_linear_system) :: this
	real, dimension(:,:), pointer :: x,res
	logical, optional :: overwrite_A
	type (error_message) :: errmsg
	character (len=35) :: myname = 'solveLeastSquaresSerialLinearSystem'
	character (len=400) :: errstr
	real, dimension(:,:), pointer :: A_tmp,x_tmp
	integer :: LDB,LWORK,INFO
	real, dimension(:), allocatable :: WORK
	logical :: overwrite
	call new(errmsg,myname)
	nullify(x,res)
	if(present(overwrite_A)) then
		overwrite = overwrite_A
	else
		overwrite = .true.
	endif
	if(overwrite) then
		call add(errmsg,0,"this routine will overwrite the system matrix in oder to "//&
		   "reduce the memory requirements",myname)
	else
		call add(errmsg,0,"this routine will NOT overwrite the system matrix, "//&
		   "hence twice as much memory for the matrix required",myname)		
	endif
!
	if(.not.associated(this%A) .or. .not.associated(this%b)) then
		call add(errmsg,2,'the unpacked system matrix, or right hand side(s) are not defined',myname)
		return
	endif
	select case(this%transpose)
	case('N')
		if(this%nrow<this%ncol) call add(errmsg,1,'calculating minimum norm solution of untransposed system, as nrow < ncol',myname)
		if(this%nrow==this%ncol) call add(errmsg,1,'calculating least square solution of untransposed system with nrow == ncol',myname)
	case('T')
		if(this%nrow>this%ncol) call add(errmsg,1,'calculating minimum norm solution of transposed system, as nrow > ncol',myname)
		if(this%nrow==this%ncol) call add(errmsg,1,'calculating least square solution of transposed system with nrow == ncol',myname)
	case default
		call add(errmsg,2,"transpose flag '"//trim(this%transpose)//"' must be either 'N' or 'T'",myname)
		return
	end select
	LDB = max(this%nrow,this%ncol)
	allocate(x_tmp(LDB,this%nrhs))
	x_tmp(:,:) = 0.
	select case(this%transpose)
	case('N'); x_tmp(1:this%nrow,:) = this%b(:,:)
	case('T'); x_tmp(1:this%ncol,:) = this%b(:,:)
	end select
	if(overwrite) then
		A_tmp => this%A
	else
		allocate(A_tmp(this%nrow,this%ncol))
		A_tmp = this%A
	endif
	! first carry out workspace query in order to optimally allocate space for WORK array
	allocate(WORK(1)); LWORK = -1
	call SGELS( this%transpose, this%nrow, this%ncol, this%nrhs, A_tmp, this%nrow, x_tmp, LDB, WORK, LWORK,INFO )
	if(INFO/=0) then
		write(errstr,"('LAPACK routine SGELS (called for workspace query) exited with INFO = ',i7)") INFO
		call add(errmsg,2,trim(errstr),myname)
 		deallocate(WORK,x_tmp)
		if(.not.overwrite) deallocate(A_tmp)
		return
	endif
	LWORK = WORK(1)
	write(errstr,"('using optimal size of WORK array = ',i10)") LWORK
	call add(errmsg,0,trim(errstr),myname)
	deallocate(WORK); allocate(WORK(LWORK))
	! now solve the system
	call SGELS( this%transpose, this%nrow, this%ncol, this%nrhs, A_tmp, this%nrow, x_tmp, LDB, WORK, LWORK,INFO )
	if(INFO/=0) then
		write(errstr,"('LAPACK routine SGELS exited with INFO = ',i7)") INFO
		call new(errmsg,2,trim(errstr),myname)
	endif
	select case(this%transpose)
	case('N') !x_tmp(1:this%nrow,:) = this%b(:,:)
		allocate(x(this%ncol,this%nrhs))
		x(:,:) = x_tmp(1:this%ncol,:)
		if(this%nrow>this%ncol) then
			allocate(res(this%nrow-this%ncol,this%nrhs))
			res(:,:) = x_tmp(this%ncol+1:this%nrow,:)
		else
			res => null()
		endif	
	case('T') !x_tmp(1:this%ncol,:) = this%b(:,:)
		allocate(x(this%nrow,this%nrhs))
		x(:,:) = x_tmp(1:this%nrow,:)
		if(this%ncol>this%nrow) then
			allocate(res(this%ncol-this%nrow,this%nrhs))
			res(:,:) = x_tmp(this%nrow+1:this%ncol,:)
		else
			res => null()
		endif	
	end select
	deallocate(WORK,x_tmp)
	if(.not.overwrite) deallocate(A_tmp)
	end function solveLeastSquaresSerialLinearSystem
!-----------------------------------------------------------------
!> \brief solve general linear system
!! \details Calling LAPACK routine "SGESV" which solves a general quadratic linear system.
!!  In this case, nrow must equal ncol. 
!! \param this serial_linear_system object
!! \param x solution array (argument "B" of LAPACK routine "SGESV")
!! \param overwrite_A logical to indicate whether system matrix A should be overwritten or not
!! \param errmsg error_message object
!! \return error message
!! \pre You need LAPACK libraries. Routine "SGESV" is called here.
!
	function solveGeneralSerialLinearSystem(this,x,overwrite_A) result(errmsg)
	type (serial_linear_system) :: this
	real, dimension(:,:), pointer :: x
	logical, optional :: overwrite_A
	type (error_message) :: errmsg
	character (len=30) :: myname = 'solveGeneralSerialLinearSystem'
	character (len=400) :: errstr
	real, dimension(:,:), pointer :: A_tmp
	integer :: INFO
	integer, dimension(:), allocatable :: IPIV
	logical :: overwrite
	call new(errmsg,myname)
	nullify(x)
	if(present(overwrite_A)) then
		overwrite = overwrite_A
	else
		overwrite = .true.
	endif
	if(overwrite) then
		call add(errmsg,0,"this routine will overwrite the system matrix in oder to "//&
		   "reduce the memory requirements",myname)
	else
		call add(errmsg,0,"this routine will NOT overwrite the system matrix, "//&
		   "hence twice as much memory for the matrix required",myname)		
	endif
!
	if(this%nrow/=this%ncol) then
		call add(errmsg,2,'nrow must equal ncol!',myname)
		return
	endif
	if(.not.associated(this%A) .or. .not.associated(this%b)) then
		call add(errmsg,2,'the unpacked system matrix, or right hand side(s) are not defined',myname)
		return
	endif
	allocate(x(this%nrow,this%nrhs))
	x = this%b
	if(overwrite) then
		A_tmp => this%A
	else
		allocate(A_tmp(this%nrow,this%ncol))
		A_tmp = this%A
	endif
	allocate(IPIV(this%nrow))
	call SGESV(this%nrow, this%nrhs, this%A, this%nrow, IPIV, x, this%nrow, INFO)
	if(INFO/=0) then
		write(errstr,"('SGESV exited with INFO = ',i7)") INFO
		call add(errmsg,2,trim(errstr),myname)
	endif
	deallocate(IPIV)
	if(.not.overwrite) deallocate(A_tmp)
	end function solveGeneralSerialLinearSystem
!-----------------------------------------------------------------
!> \brief solve packed symmetric linear system serially
!! \details Calling LAPACK routine "SSPSV" which solves an indefinite symmetric 
!!  linear system where the triangular system matrix is stored in a 1d array. 
!! \param this serial_linear_system object
!! \param A columnwise packed triangular matrix of linear system
!! \param uplo one charcter indicating if A contains columnwisely the 'U'pper or 'L'ower triangular matrix
!! \param x solution array (argument "B" of LAPACK routine "SSPSV")
!! \param overwrite logical to indicate whether system matrix A should be overwritten or not
!! \param errmsg error_message object
!! \return error message
!! \pre You need LAPACK libraries. Routine "SSPSV" is called here.
!
	function solvePackedSymmetricSerialLinearSystem(this,x,overwrite_A) result(errmsg)
	type (serial_linear_system) :: this
	real, dimension(:,:), pointer :: x
	logical, optional :: overwrite_A
	type (error_message) :: errmsg
	character (len=38) :: myname = 'solvePackedSymmetricSerialLinearSystem'
	character (len=400) :: errstr
	real, dimension(:), pointer :: A_tmp
	integer, dimension(:), allocatable :: IPIV
	integer :: INFO
	logical :: overwrite
	call new(errmsg,myname)
	nullify(x)
	if(present(overwrite_A)) then
		overwrite = overwrite_A
	else
		overwrite = .true.
	endif
	if(overwrite) then
		call add(errmsg,0,"this routine will overwrite the system matrix in oder to "//&
		   "reduce the memory requirements",myname)
	else
		call add(errmsg,0,"this routine will NOT overwrite the system matrix, "//&
		   "hence twice as much memory for the matrix required",myname)		
	endif
!
	if(this%nrow/=this%ncol) then
		call add(errmsg,2,'nrow must equal ncol!',myname)
		return
	endif
	if(.not.associated(this%A_packed) .or. .not.associated(this%b)) then
		call add(errmsg,2,'the packed system matrix, or right hand side(s) are not defined',myname)
		return
	endif
	if(size(this%A_packed) /= (this%nrow*(this%nrow+1))/2) then
		write(errstr,"(a,i15,a,i7,a)") "the size of packed system matrix ( = ",size(this%A_packed), &
		  " ) does not equal N*(N+1)/2, where N is the number of rows ( = number of columns = ", &
		  this%nrow," )"
		call add(errmsg,2,trim(errstr),myname)
		return
	endif
	allocate(x(this%nrow,this%nrhs))
	x = this%b
	if(overwrite) then
		A_tmp => this%A_packed
	else
		allocate(A_tmp(size(this%A_packed)))
		A_tmp = this%A_packed
	endif
	! solve the system
	allocate(IPIV(this%nrow))
	!subroutine SSPSV (UPLO, N, NRHS, AP, IPIV, B, LDB, INFO)
	call SSPSV(this%UPLO, this%nrow, this%nrhs, A_tmp, IPIV, x, this%nrow, INFO )
	if(INFO/=0) then
		write(errstr,"('SSPSV exited with INFO = ',i7)") INFO
		call add(errmsg,2,trim(errstr),myname)
	endif
	deallocate(IPIV)
	if(.not.overwrite) deallocate(A_tmp)
	end function solvePackedSymmetricSerialLinearSystem
!-----------------------------------------------------------------
!> \brief get pointer to system matrix
!! \param this linear_system object
!! \param p pointer to system matrix
!! \return pointer to system matrix
!
	function getMatrixSerialLinearSystem(this) result(p)
	type (serial_linear_system), intent(in) :: this
	real, dimension(:,:), pointer :: p
	p => this%A
	end function getMatrixSerialLinearSystem
!-----------------------------------------------------------------
!> \brief get pointer to packed triangular system matrix
!! \param this linear_system object
!! \param p pointer to packed triangular system matrix
!! \return pointer to system matrix
!
	function getPackedMatrixSerialLinearSystem(this) result(p)
	type (serial_linear_system), intent(in) :: this
	real, dimension(:), pointer :: p
	p => this%A_packed
	end function getPackedMatrixSerialLinearSystem
!
end module serialLinearSystem
