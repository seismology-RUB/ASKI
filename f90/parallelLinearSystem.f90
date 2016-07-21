!----------------------------------------------------------------------------
!   Copyright 2013 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------
! ########################################################
! ### THIS MODULE IS WORK IN PROGRESS (obviously)
! ########################################################
!----------------------------------------------------------------
!> \brief handle linear systems and solve them with SCALAPACK
!!
!! \details Module that handles linear systems of equations which
!!  can be solved with parallel SCALAPACK algorithms.
!!  Several solvers for different types of matrices
!!  are supported.
!!
!! \author Florian Schumacher
!! \date 2012
!---------------------------------------------------------------
module parallelLinearSystem
!
	use errorMessage
!
	implicit none
	!include 'mpif.h'
!
!!$  propably does not work in connection with module linearSystem (??) (compare module sequentialLinearSystem)
!!$!> \brief overloading with name consistent with module linearSystem
!!$	interface createLinearSystem
!!$		module procedure createParallelLinearSystem
!!$	end interface createLinearSystem
!> \brief overloading deallocate routine with shorter name
	interface dealloc
		module procedure deallocParallelLinearSystem
	end interface dealloc
	!interface operator (.A.); module procedure getMatrixLinearSystem; end interface
!
!> \brief all specifications of the linear system
	type parallel_linear_system
		private
		integer :: ICTXT !< index of BLACS context
		integer :: NPROW !< ??
		integer :: NPCOL !< ??
		integer :: NBROW !< ??
		integer :: NBCOL !< ??
		integer :: MYROW !< ??
		integer :: MYCOL !< ??
		!integer :: nrow,ncol,nrhs                  ! number of rows, number of columns, number of right-hand-sides  of linear system
		!real, dimension(:,:), pointer :: A         ! nrow-by-ncol system matrix
		!real, dimension(:,:), pointer :: b         ! right-hand-side(s) of system. nrow-by-nrhs array
	end type parallel_linear_system
!
contains
!
!------------------------------------------------------------------------
!> \brief create parallel linear system by allocating and setting LSE_seq
!! \details Set parallel_linear_system object and initiate parallel process...
!! \param this parallel_linear_system object
!! \param NPROW ??
!! \param NPCOL ??
!! \param NBROW ??
!! \param NBCOL ??
!! <!-- \param A system matrix -->
!! <!-- \param b matrix of right-hand-sides -->
!! \pre You need SCALAPACK libraries (including BLACS etc.).
!
	subroutine createParallelLinearSystem(this,NPROW,NPCOL,NBROW,NBCOL)!,A,b)
	type (parallel_linear_system) :: this
	integer :: NPROW,NPCOL,NBROW,NBCOL
	integer :: ICTXT,MYROW,MYCOL
	!if(associated(this%A) .or. associated(this%b)) call deallocParallelLinearSystem(this)
	!! if(this%ICTXT=??? -> check if parallel environment is running!, if so: deallocate (in order to call BLACKS_EXIT)
	this%NPROW=NPROW; this%NPCOL=NPCOL
	this%NBROW=NBROW; this%NBCOL=NBCOL
	print *,"vor SL_INIT, ICTXT, NPROW, NPCOL=",ICTXT, NPROW, NPCOL
	call SL_INIT(ICTXT,NPROW,NPCOL)
	print *,"nach SL_INIT, ICTXT, NPROW, NPCOL=",ICTXT, NPROW, NPCOL
	this%ICTXT=ICTXT
	print *,"createParallelLinearSystem,ICTXT=",this%ICTXT
	if(ICTXT<0) return
	call BLACS_GRIDINFO(ICTXT,NPROW,NPCOL,MYROW,MYCOL)
	this%MYROW=MYROW; this%MYCOL=MYCOL
	print *,"MYROW,MYCOL=",MYROW,MYCOL,";  ICTXT=",ICTXT,";  NPROW,NPCOL,NBROW,NBCOL=",NPROW,NPCOL,NBROW,NBCOL
!Initialize the array descriptors for the matrices A and b
	!call DESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO )
	!call DESCINIT( DESCB, N, NRHS, NB, NBRHS, RSRC, CSRC, ICTXT, MXLLDB, INFO )
	end subroutine createParallelLinearSystem
!-----------------------------------------------------------------
!> \brief deallocate parallel linear system
!! \details Close parallel process. Nullify pointers to...
!! \param this sequential_linear_system object
!! \pre You need SCALAPACK libraries (including BLACS etc.).
!
	subroutine deallocParallelLinearSystem(this)
	type (parallel_linear_system) :: this
	!if(associated(this%A)) this%A => null()
	!if(associated(this%b)) this%b => null()
	if(this%ICTXT.ge.0) call BLACS_GRIDEXIT(this%ICTXT)
	!Argument "0" of BLACS_EXIT: Flag indicating whether message passing continues after the BLACS are done. 
	!If continue is non-zero, the user is assumed to continue using the machine after completing the BLACS. 
	!Otherwise, no message passing is assumed after calling this routine.
	call BLACS_EXIT( 0 ) 
	end subroutine deallocParallelLinearSystem
!-----------------------------------------------------------------
!> \brief solve linear system by least squares
!! \details Calling SCALAPACK routine (??) which ...
!! \param this parallel_linear_system object
!! \param errmsg error_message object
!! \return error message
!! \pre You need SCALAPACK libraries (including BLACS etc.).
!
	function solveLeastSquaresParallelLinearSystem(this) result(errmsg)
	type (parallel_linear_system) :: this
	type (error_message) :: errmsg
	character (len=37) :: myname = 'solveLeastSquaresParallelLinearSystem'
	call new(errmsg,myname)
	call new(errmsg,1,'this is a warning message from this fake routine!!',myname)
	end function solveLeastSquaresParallelLinearSystem
!-----------------------------------------------------------------
!> \brief solve general linear system in parallel
!! \details Calling SCALAPACK routine (??) which ...
!! \param this parallel_linear_system object
!! \param errmsg error_message object
!! \return error message
!! \pre You need SCALAPACK libraries (including BLACS etc.).
!
	function solveGeneralParallelLinearSystem(this) result(errmsg)
	type (parallel_linear_system) :: this
	type (error_message) :: errmsg
	character (len=32) :: myname = 'solveGeneralParallelLinearSystem'
	call new(errmsg,myname)
	call new(errmsg,1,'this is a test warning message from this fake routine',myname)
	end function solveGeneralParallelLinearSystem
!------------------------------------------------------------------
!!$ THERE IS NO PARALLEL SOLVER FOR INDEFINITE SYMMETRIC LINEAR SYSTEMS
!!$!-----------------------------------------------------------------
!!$!> \brief solve symmetric linear system in parallel
!!$!! \details Calling SCALAPACK routine (??) which ...
!!$!! \param this parallel_linear_system object
!!$!! \param errmsg error_message object
!!$!! \return error message
!!$!! \pre You need SCALAPACK libraries (including BLACS etc.).
!!$	function solveSymmetricParallelLinearSystem(this) result(errmsg)
!!$	type (parallel_linear_system) :: this
!!$	type (error_message) :: errmsg
!!$	character (len=34) :: myname = 'solveSymmetricParallelLinearSystem'
!!$	call new(errmsg,myname)
!!$	call new(errmsg,1,'this is a test warning message from this fake routine',myname)
!!$	end function solveSymmetricParallelLinearSystem
!-----------------------------------------------------------------
!> \brief get pointer to system matrix
!! \param this linear_system object
!! \param p pointer to (local?) system matrix (?)
!! \return pointer to (local?) system matrix(?)
!
	function getMatrixParallelLinearSystem(this) result(p)
	type (parallel_linear_system), intent(in) :: this
	real, dimension(:,:), pointer :: p
	p => null()
	end function getMatrixParallelLinearSystem
!
end module parallelLinearSystem
