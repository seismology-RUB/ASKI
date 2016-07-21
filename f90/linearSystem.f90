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
!> \brief handle linear systems and solve them with LAPACK or SCALAPACK
!!
!! \details Generic module that handles linear systems of equations by using
!!  modules serialLinearSystem and parallelLinearSystem. \n
!!  By allocating either one of the fields LSE_ser or LSE_par of an object of type linear_system, 
!!  this module allows to define a linear system which 
!!  is either solved serially or in parallel, but using the same generic commands (which 
!!  then call the respective serial or parallel routines from the respective module).
!!
!! \author Florian Schumacher
!! \date Mar 2012
!---------------------------------------------------------------
module linearSystem
!
	use serialLinearSystem
	use parallelLinearSystem
	use errorMessage
!
	implicit none
	!include 'mpif.h'
!
!> \brief overloading create routines
	interface createLinearSystem
		module procedure createSeqLinearSystem
		module procedure createParLinearSystem
	end interface createLinearSystem
!> \brief overloading deallocate routine with shorter name
	interface dealloc
		module procedure deallocLinearSystem
	end interface dealloc
	interface operator (.A.); module procedure getMatrixLinearSystem; end interface
!
!> \brief generic type which branches either to type serial_linear_system or 
!! parallel_linear_system
	type linear_system
		private
		type (serial_linear_system), pointer :: LSE_ser => null() !< pointer to serial linear system
		type (parallel_linear_system), pointer :: LSE_par => null() !< pointer to parallel linear system
	end type linear_system
!
contains
!
!------------------------------------------------------------------------
!> \brief create serial linear system
!! \details Allocate the pointer to field LSE_ser and then call createSerialLinearSystem
!!  from module serialLinearSystem. This way, this generic linear_system object will
!!  be handled as a serial_linear_system object.
!! \param this linear_system object
!! \param nrow number of rows of linear system
!! \param ncol number of columns of linear system
!! \param nrhs number of right-hand-sides of linear system
!! \param A pointer to system matrix
!! \param b pointer to matrix of right-hand-sides
	subroutine createSeqLinearSystem(this,nrow,ncol,nrhs,A,b)
	type (linear_system) :: this
	integer :: nrow,ncol,nrhs
	real, dimension(:,:), pointer :: A
	real, dimension(:,:), pointer :: b
	if(associated(this%LSE_ser) .or. associated(this%LSE_par)) call deallocLinearSystem(this)
	allocate(this%LSE_ser)
	call createSerialLinearSystem(this%LSE_ser,nrow,ncol,nrhs,A,b)
	end subroutine createSeqLinearSystem
!------------------------------------------------------------------------
! THIS ROUTINE IS WORK IN PROGRESS (obviously)
!> \brief create parallel linear system
!! \details Allocate the pointer to field LSE_par and then call createParallelLinearSystem
!!  from module parallelLinearSystem. This way, this generic linear_system object will
!!  be handled as a parallel_linear_system object.
!! \param this linear_system object
!! \param NPROW number of rows of linear system
!! \param NPCOL number of columns of linear system
!! \param NBROW block size (?)
!! \param NBCOL block size (?)
!$!! \param A pointer to local system matrix block
!$!! \param b pointer to local matrix block of right-hand-sides
	subroutine createParLinearSystem(this,NPROW,NPCOL,NBROW,NBCOL)!,A,b)
	type (linear_system) :: this
	integer :: NPROW,NPCOL,NBROW,NBCOL
	integer :: ICTXT,MYROW,MYCOL
	if(associated(this%LSE_ser) .or. associated(this%LSE_par)) call deallocLinearSystem(this)
	allocate(this%LSE_par)
	call createParallelLinearSystem(this%LSE_par,NPROW,NPCOL,NBROW,NBCOL)!,A,b)
	end subroutine createParLinearSystem
!-----------------------------------------------------------------
!> \brief deallocate linear system
!! \details Generic rountine that deallocates fields LSE_ser and LSE_par contained
!!  in linear_system object after calling routines deallocSerialLinearSystem and
!!  deallocParallelLinearSystem from modules serialLinearSystem and parallelLinearSystem,
!!  respectively.
!! \param this linear_system object
	subroutine deallocLinearSystem(this)
	type (linear_system) :: this
		if(associated(this%LSE_ser)) then
			call deallocSerialLinearSystem(this%LSE_ser)
			deallocate(this%LSE_ser)
		endif
		if(associated(this%LSE_par)) then
			call deallocParallelLinearSystem(this%LSE_par)
			deallocate(this%LSE_par)
		endif
	end subroutine deallocLinearSystem
!-----------------------------------------------------------------
!> \brief solve linear system by least squares
!! \details Generic routine that calls either solveLeastSquaresSerialLinearSystem from 
!!  module serialLinearSystem or solveLeastSquaresParallelLinearSystem from module
!!  parallelLinearSystem, dependet on which pointer is associated in linear_system object.
!! \param this linear_system object
!! \param x pointer to solution array returned by the routine which was called
!! \param override_A optional logical to indicate whether system matrix A should be overwritten or not (default: .true.)
!! \param errmsg error_message object
!! \return error message
	function solveLeastSquaresLinearSystem(this,x,res,override_A) result(errmsg)
	type (linear_system) :: this
	real, dimension(:,:), pointer :: x,res
	logical, optional :: override_A
	type (error_message) :: errmsg
	character (len=29) :: myname = 'solveLeastSquaresLinearSystem'
	logical :: override
	if(present(override_A)) then 
		override = override_A
	else
		override = .true.
	endif
	if (associated(this%LSE_ser)) then
		errmsg = solveLeastSquaresSerialLinearSystem(this%LSE_ser,x,res,override)
	else if (associated(this%LSE_par)) then
		errmsg = solveLeastSquaresParallelLinearSystem(this%LSE_par)
	endif
	call addTraceErrorMessage(errmsg,myname)
	if (.level.errmsg == 2) return
	end function solveLeastSquaresLinearSystem
!-----------------------------------------------------------------
!> \brief solve general linear system
!! \details Generic routine that calls either solveGeneralSerialLinearSystem from 
!!  module serialLinearSystem or solveGeneralParallelLinearSystem from module
!!  parallelLinearSystem, dependet on which pointer is associated in linear_system object.
!! \param this linear_system object
!! \param x pointer to solution array returned by the routine which was called
!! \param override_A optional logical to indicate whether system matrix A should be overwritten or not (default: .true.)
!! \return error message
	function solveGeneralLinearSystem(this,x,override_A) result(errmsg)
	type (linear_system) :: this
	real, dimension(:,:), pointer :: x
	logical, optional :: override_A
	type (error_message) :: errmsg
	character (len=24) :: myname = 'solveGeneralLinearSystem'
	logical :: override
	if(present(override_A)) then 
		override = override_A
	else
		override = .true.
	endif
	if (associated(this%LSE_ser)) then
		errmsg = solveGeneralSerialLinearSystem(this%LSE_ser,x,override)
	else if (associated(this%LSE_par)) then
		errmsg = solveGeneralParallelLinearSystem(this%LSE_par)
	endif
	call addTraceErrorMessage(errmsg,myname)
	if (.level.errmsg == 2) return
	end function solveGeneralLinearSystem
!------------------------------------------------------------------------
!!$ THERE IS NO PARALLEL SOLVER FOR INDEFINITE SYMMETRIC LINEAR SYSTEMS, SO GENERIC ROUTINE IS USELESS
!!$!------------------------------------------------------------------------
!!$!> \brief solve symmetric linear system
!!$!! \details Generic routine that calls either solveSymmetricSerialLinearSystem from 
!!$!!  module serialLinearSystem or solveSymmetricParallelLinearSystem from module
!!$!!  parallelLinearSystem, dependet on which pointer is associated in linear_system object.
!!$!! \param this linear_system object
!!$!! \return 
!!$!
!!$	function solveSymmetricLinearSystem(this,x,override_A) result(errmsg)
!!$	type (linear_system) :: this
!!$	real, dimension(:,:), pointer :: x
!!$	logical, optional :: override_A
!!$	type (error_message) :: errmsg
!!$	character (len=26) :: myname = 'solveSymmetricLinearSystem'
!!$	logical :: override
!!$	if(present(override_A)) then 
!!$		override = override_A
!!$	else
!!$		override = .true.
!!$	endif
!!$	if (associated(this%LSE_ser)) then
!!$		errmsg = solveSymmetricSerialLinearSystem(this%LSE_ser,x,override)
!!$	else if (associated(this%LSE_par)) then
!!$		errmsg = solveSymmetricParallelLinearSystem(this%LSE_par)
!!$	endif
!!$	call addTraceErrorMessage(errmsg,myname)
!!$	if (.level.errmsg == 2) return
!!$	end function solveSymmetricLinearSystem
!-----------------------------------------------------------------
!> \brief get pointer to system matrix
!! \param this linear_system object
!! \param p pointer to system matrix
!! \return pointer to system matrix
	function getMatrixLinearSystem(this) result(p)
	type (linear_system), intent(in) :: this
	real, dimension(:,:), pointer :: p
	if(associated(this%LSE_ser)) then
		p => getMatrixSerialLinearSystem(this%LSE_ser)
	elseif(associated(this%LSE_par)) then
		p => null()
	else
		p => null()
	endif
	end function getMatrixLinearSystem
!
end module linearSystem
