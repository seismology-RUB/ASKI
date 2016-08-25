!----------------------------------------------------------------------------
!	Copyright 2016 Wolfgang Friederich
!
!	This file is part of Gemini II.
!
!	Gemini II is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 2 of the License, or
!	any later version.
!
!	Gemini II is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------
!> \brief Support tools for using MPI

!> \author Wolfgang Friederich

!> \par Description
!>  Module to ease MPI parallelization of code 
!!  Creates an object that contains information about the parallel
!!  environment.
!<---------------------------------------------------------------
 module mpiSupport
       use mpi
	implicit none
	interface new; module procedure createMpiSupport; end interface
	interface dealloc; module procedure deallocMpiSupport; end interface
	interface abort; module procedure abortRunMpiSupport; end interface
	interface barrier; module procedure barrierMpiSupport; end interface
	interface operator (.myrank.); module procedure myrankMpiSupport; end interface
	interface operator (.numtasks.); module procedure numtasksMpiSupport; end interface
	type mpi_support
		private
		integer myrank                           ! Rank of this process
		integer numtasks                         ! Number of parallel tasks
	end type
!
 contains
!-----------------------------------------------------------------
!> \brief  Initialize MPI
!
	subroutine createMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr,rc
!
	call mpi_init(ierr)
	if (ierr .ne. MPI_SUCCESS) then
		print *,'Error starting MPI program. Terminating.'
		call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
	endif
!
	call MPI_COMM_RANK(MPI_COMM_WORLD, this%myrank, ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, this%numtasks, ierr)
	end subroutine createMpiSupport
!-----------------------------------------------------------------
!> \brief Dealloc support object and finalize MPI
!
	subroutine deallocMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr
	call MPI_FINALIZE(ierr)
	end subroutine deallocMpiSupport
!-----------------------------------------------------------------
!> \brief Abort mpi run
!
	subroutine abortRunMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr,rc
	call MPI_ABORT(MPI_COMM_WORLD,rc,ierr)
	end subroutine abortRunMpiSupport
!-----------------------------------------------------------------
!> \brief Mpi Barrier
!
	subroutine barrierMpiSupport(this)
	type (mpi_support) :: this
	integer :: ierr
	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
	end subroutine barrierMpiSupport
!-----------------------------------------------------------------
!> \brief get my rank
!
	integer function myrankMpiSupport(this)
	type (mpi_support), intent(in) :: this
	myrankMpiSupport = this%myrank
	end function myrankMpiSupport
!-----------------------------------------------------------------
!> \brief get number of parallel taks
!
	integer function numtasksMpiSupport(this)
	type (mpi_support), intent(in) :: this
	numtasksMpiSupport = this%numtasks
	end function numtasksMpiSupport
!
 end module
