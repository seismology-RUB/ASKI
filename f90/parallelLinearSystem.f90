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
     module procedure deallocateParallelLinearSystem
  end interface dealloc
  interface im_in
     module procedure getImInTheGridParallelLinearSystem
  end interface im_in

  interface operator (.myrow.); module procedure getMyrowParallelLinearSystem; end interface
  interface operator (.mycol.); module procedure getMycolParallelLinearSystem; end interface
  !interface operator (.A.); module procedure getMatrixLinearSystem; end interface
  !
  !> \brief all specifications of the linear system
  type parallel_linear_system
     private
     ! IT IS ASSUMED THAT IRSRC == ICSRC == 0, i.e. the distribution of the blocks of the global matrix starts at process {0,0} in the process grid

     ! PARALLEL ENVIRONMENT
     integer :: blacs_context = -1
     integer :: nprow = -1 !< total number of rows in process grid
     integer :: npcol = -1 !< total number of columns in process grid
     integer :: myrow = -1 !< row coordinate of this process in process grid
     integer :: mycol = -1 !< column coordinate of this process in process grid
     integer :: nbrow = -1 !< number of rows in a block (blacs blocking factor)
     integer :: nbcol = -1 !< number of columns in a block (blacs blocking factor)
     logical :: im_in_the_grid = .false. !< indicating whether this process is in the grid or not
     logical :: initiated = .false. !< .true. if initParallelLinearSystem exited successfully
     logical :: blacs_initiated = .false. !< solely tells if the parallel environment was invoked yet by SL_INIT or not

     ! LINEAR SYSTEM

     ! GLOBAL INFORMATION
     integer :: nrow = -1 !< number of rows in global system
     integer :: ncol = -1 !< number of columns in global system
     integer :: nrhs = -1 !< number of right-hand-sides in global system
     ! global rhs-array is assumed to have max(nrow,ncol) rows in order to store the solution of the system

     ! LOCAL INFORMATION
     integer :: lrowA = -1 !< local number of rows (i.e. leading dimension) of (local) submatrix LA
     integer :: lcolA = -1 !< local number of columns of (local) submatrix LA
     integer :: lrowb = -1 !< local number of rows of rhs / solution array Lb (computed from max(nrow,ncol))
     integer :: lcolb = -1 !< local number of (columns of) right-hand-sides of (local) subrhs Lb
     real, dimension(:,:), pointer :: LA => null() ! local submatrix (size lrowA-times-lcolA)
     integer, dimension(9) :: descA !< SCALAPACL descriptor for matrix A
     real, dimension(:,:), pointer :: Lb => null()  ! local right-hand-side(s) of system (size lrowb-times-lrhs)
     integer, dimension(9) :: descb !< SCALAPACL descriptor for rhs b
  end type parallel_linear_system
  !
contains
  !
  !------------------------------------------------------------------------
  !> \brief initiate parallel linear system by allocating the matrices and initiating the parallel environment
  !! \details Set parallel_linear_system object and initiate parallel process...
  !! \param this parallel_linear_system object
  !! \param NPROW requested number of processes for the rows
  !! \param NPCOL requested number of processes for the columns
  !! \pre You need SCALAPACK libraries (including BLACS etc.).
  !
  subroutine initParallelLinearSystem(this,nrow,ncol,nrhs,nprow,npcol,nbrow,nbcol,errmsg)
    ! incoming
    type (parallel_linear_system) :: this
    integer :: nrow,ncol,nrhs
    integer :: nprow,npcol,nbrow,nbcol
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=24) :: myname = 'initParallelLinearSystem'
    character(len=400) :: errstr
    integer :: INFO

    integer NUMROC
    external NUMROC
!
    call addTrace(errmsg,myname)
!
    if( any((/nrow,ncol,nrhs/) <= 0) ) then ! ,nprow,npcol,nbrow,nbcol
       write(errstr,*) "some of the incoming global system specifications are invalid, must be > 0: number of rows = ",&
            nrow,", number of columns = ",ncol,", number of right-hand-sides of the system = ",nrhs
       call add(errmsg,2,errstr,myname)
       return
    end if
    if( any((/nprow,npcol,nbrow,nbcol/) <= 0) ) then
       write(errstr,*) "some of the incoming parallelization specifications are invalid, must be > 0: ",&
            "number of rows,columns in ScaLAPACK submatrix block = ",nbrow,", ",nbcol,&
            "; number of processes in rows,columns of BLACS process grid = ",npcol,", ",npcol
       call add(errmsg,2,errstr,myname)
       return
    end if
    this%nrow = nrow; this%ncol = ncol; this%nrhs = nrhs
    this%nbrow = nbrow; this%nbcol = nbcol


!!$!  BLACS_PINFO tells this process its number, and the total number of
!!$!  processes.
!!$    call blacs_pinfo(this%blacks_iam,nproc)
!!$
!!$!  If machine needs additional set up, do it now:
!!$!  we only want one process to do the setup, choose process 0 
!!$!  to do it
!!$  if(nproc < 1) then
!!$    if(this%blacks_iam == 0) nproc = nprow * npcol
!!$    call blacs_setup(this%blacks_iam,nprocs)
!!$  end if
!!$    
!!$!  ask BLACS_GET to return the default system context
!!$  call blacs_get(-1,0,context)
!!$
!!$!  BLACS_GRIDINIT sets up the process grid of NPROW by NPCOL
!!$!  processes, in row-major order.  
!!$!  BLACS_GRIDINIT replaces the default context with the BLACS context for this grid.
!!$  call blacs_gridinit (context,'Row-major',nprow,npcol)

    ! the above commands are called by SL_INIT in exactly the same way
    call SL_INIT(this%blacs_context,nprow,npcol)
    this%blacs_initiated = .true.

!  BLACS_GRIDINFO tells each process the "position" it has been
!  assigned in the process grid.
!
!  CONTEXT is input, the context handle.
!  NPROW and NPCOL are output, the number of processor rows and columns. ( (HOW) ARE THEY MODIFIED?!)
!  MYROW and MYCOL are output, this process's row and column.
!
    call blacs_gridinfo ( this%blacs_context, this%nprow, this%npcol, this%myrow, this%mycol )
    if(all((/this%nprow, this%npcol, this%myrow, this%mycol/) ==-1 )) then
       this%im_in_the_grid = .false.
       return ! if I am not in the grid, return, there is nothing to do. wait for blacs_exit()     
    end if

    if(this%nprow /= nprow) then
       write(errstr,*) "nprow before calling blacs_gridinfo = ", nprow,"; after = ",this%nprow,&
            "; myrow,mycol = ",this%myrow, this%mycol
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(this%npcol /= npcol) then
       write(errstr,*) "npcol before calling blacs_gridinfo = ", npcol,"; after = ",this%npcol,&
            "; myrow,mycol = ",this%myrow, this%mycol
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(this%nprow <= 0 .or. this%npcol <= 0) then
       call add(errmsg,2,"after blacs_gridinfo, nprow or npcol are <= 0, there is a problem with the parallel environment",myname)
       return
    end if

    ! determine whether this process is in the BLACS process grid or not
    if(this%myrow >= this%nprow .or. this%mycol >= this%npcol) then
       this%im_in_the_grid = .false.
       return ! if I am not in the grid, return, there is nothing to do. wait for blacs_exit()
    else
       this%im_in_the_grid = .true.
    end if

    ! FROM HERE ON, I KNOW THAT I AM IN THE GRID, SO PREPARE LOCAL ARRAY INFORMATION

    ! determine the local sizes of the local submatrix and the local rhs vector(s) by function
    ! NUMROC (NUMber of Rows Or Columns)
    ! INTEGER FUNCTION NUMROC( N, NB, IPROC, ISRCPROC, NPROCS )
    this%lrowA = NUMROC( this%nrow, this%nbrow, this%myrow, 0, this%nprow )
    this%lcolA = NUMROC( this%ncol, this%nbcol, this%mycol, 0, this%npcol )
    !write(*,*) "myrow,mycol,size(A) = ",this%myrow,this%mycol,this%lrowA,this%lcolA
    ! global rhs-array is assumed to have max(nrow,ncol) rows in order to store the solution of the system:
    this%lrowb = NUMROC( max(this%nrow,this%ncol), this%nbrow, this%myrow, 0, this%nprow )
    this%lcolb = NUMROC( this%nrhs, this%nbcol, this%mycol, 0, this%npcol )
    !write(*,*) "myrow,mycol,size(b) = ",this%myrow,this%mycol,this%lrowb,this%lcolb

    ! Initialize the array descriptors for the local sub matrices A and b
    ! SUBROUTINE DESCINIT( DESC, M, N, MB, NB, IRSRC, ICSRC, ICTXT, LLD, INFO )
    call DESCINIT( this%descA, this%nrow, this%ncol, this%nbrow, this%nbcol, 0, 0, this%blacs_context, &
         max(1,this%lrowA), INFO )
    if(INFO/=0) then
       write(errstr,*) "descriptor initiation of A failed, returned INFO = ",INFO
       call add(errmsg,2,errstr,myname)
       return
    end if
    !write(*,*) "myrow,mycol,descA = ",this%myrow,this%mycol,this%descA
    call DESCINIT( this%descb, max(this%nrow,this%ncol), this%nrhs, this%nbrow, this%nbcol, 0, 0, &
         this%blacs_context, max(1,this%lrowb), INFO )
    if(INFO/=0) then
       call add(errmsg,2,"descriptor initiation of b failed",myname)
       return
    end if
    !write(*,*) "myrow,mycol,descb = ",this%myrow,this%mycol,this%descb
    this%initiated = .true.
  end subroutine initParallelLinearSystem
  !------------------------------------------------------------------------
!> \brief allocate local storage for arrays LA and Lb
  subroutine allocateLocalStorageParallelLinearSystem(this,A,b,errmsg)
    type (parallel_linear_system) :: this
    real, dimension(:,:), pointer :: A,b
    type (error_message) :: errmsg
    call addTrace(errmsg,'allocateLocalStorageParallelLinearSystem')
    nullify(A,b)
    if(.not.this%im_in_the_grid) return
    if(.not.this%initiated) then
       call add(errmsg,2,"parallel linear system not yet initiated",'allocateLocalStorageParallelLinearSystem')
       return
    end if
    if(associated(this%LA)) deallocate(this%LA)
    if(associated(this%Lb)) deallocate(this%Lb)
    allocate(this%LA(this%lrowA,this%lcolA),this%Lb(this%lrowb,this%lcolb))
    A => this%LA
    b => this%lb
  end subroutine allocateLocalStorageParallelLinearSystem
  !-----------------------------------------------------------------
  !> \brief deallocate parallel linear system
  !! \details Close parallel process. Nullify pointers to...
  !! \param this sequential_linear_system object
  !! \pre You need SCALAPACK libraries (including BLACS etc.).
  !
  subroutine deallocateParallelLinearSystem(this)
    type (parallel_linear_system) :: this
    !if(associated(this%A)) this%A => null()
    !if(associated(this%b)) this%b => null()
!
    ! write(*,*) "###",this%blacs_context,this%myrow,this%mycol,this%nprow,this%npcol,this%im_in_the_grid,&
    !      this%initiated,this%blacs_initiated
    ! free the process grid, if I am in the process grid
    if(this%blacs_context.ge.0 .and. this%blacs_initiated) then
       !call BLACS_BARRIER(this%blacs_context,'A') ! a barrier does not seem to be necessary
       call BLACS_GRIDEXIT(this%blacs_context)
    end if
!
    !  For processes which are not in the process grid, this%blacs_context has value -1 (so no BLACS_GRIDEXIT requirerd),
    !  but still BLACS_EXIT is required to be called! So additionally check whether parallel environment was initiated or not before calling BLACS_EXIT. Otherwise mpi will complain that it finalizes before it was initiated
    !Argument "0" of BLACS_EXIT: Flag indicating whether message passing continues after the BLACS are done. 
    !If argument is non-zero, the user is assumed to continue using the machine after completing the BLACS. 
    !Otherwise, no message passing is assumed after calling this routine.
    if(this%blacs_initiated) call BLACS_EXIT( 0 )
!
    this%blacs_context = -1
    this%nprow = -1
    this%npcol = -1
    this%myrow = -1
    this%mycol = -1
    this%nbrow = -1
    this%nbcol = -1
    this%im_in_the_grid = .false.
    this%initiated = .false.
    this%blacs_initiated = .false.
!
    this%nrow = -1
    this%ncol = -1
    this%nrhs = -1
!
    this%lrowA = -1
    this%lcolA = -1
    this%lrowb = -1
    this%lcolb = -1
    if(associated(this%LA)) deallocate(this%LA)
    this%descA(:) = -1
    if(associated(this%Lb)) deallocate(this%Lb)
    this%descb(:) = -1
  end subroutine deallocateParallelLinearSystem
  !------------------------------------------------------------------------
  subroutine distributeGlobalSubmatrixParallelLinearSystem(this,which_array,mat_send,nrow_send,ncol_send,&
       ig_start,jg_start,sendprow,sendpcol,errmsg)
    ! incoming
    type (parallel_linear_system) :: this
    character(len=*) :: which_array
    real, dimension(:,:) :: mat_send
    integer :: nrow_send,ncol_send,ig_start,jg_start
    integer :: sendprow,sendpcol
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=45) :: myname = 'distributeGlobalSubmatrixParallelLinearSystem'
    character(len=400) :: errstr
    ! local array to add values to, fork to either A or b
    real, dimension(:,:), pointer :: A_or_b
    integer, dimension(9) :: desc_A_or_b
    integer :: nrow_A_or_b,ncol_A_or_b,lrow_A_or_b,lcol_A_or_b
    ! regarding rows
    integer :: nrow_send_tmp,i_block,i_block_start,i_block_end,m_block_send
    integer :: ig,il,nr,i_mat_send,recvprow
    integer, dimension(:), allocatable :: nrow_of_block,ig_start_of_block
    ! regarding columns
    integer :: ncol_send_tmp,j_block,j_block_start,j_block_end,n_block_send
    integer :: jg,jl,nc,j_mat_send,recvpcol
    integer, dimension(:), allocatable :: ncol_of_block,jg_start_of_block
!
    ! IT IS ASSUMED THAT ALL PROCESSES OF THE PROCESS GRID GET THE SAME VALUES OF 
    !   which_array,nrow_send,ncol_send,ig_start,jg_start,sendprow,sendpcol !!
!
    call addTrace(errmsg,myname)
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
    if(.not.this%initiated) then
       call add(errmsg,2,"parallel linear system not yet initiated",myname)
       return
    end if
    if(.not.associated(this%LA)) then
       call add(errmsg,2,"parallel linear system not yet allocated",myname)
       return
    end if
    if(sendprow < 0 .or. sendprow > this%nprow .or. sendpcol < 0 .or. sendpcol > this%npcol) then
       write(errstr,*) "requested process of distribution is (row,col) = ",sendprow,sendpcol,&
            "; however, maximum dimension of process grid is (nrow,ncol) = ",this%nprow,this%npcol
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! check incoming type of values (matrix or rhs) to distribute
    select case(which_array)
    case ('matrix')
       A_or_b => this%LA
       nrow_A_or_b = this%nrow
       ncol_A_or_b = this%ncol
       lrow_A_or_b = this%lrowA
       lcol_A_or_b = this%lcolA
       desc_A_or_b = this%descA
       write(errstr,*) "will distribute incoming values to system matrix, which has global dimensions ",&
            nrow_A_or_b," - by - ",ncol_A_or_b
       call add(errmsg,0,errstr,myname)
    case ('rhs')
       A_or_b => this%Lb
       nrow_A_or_b = max(this%nrow,this%ncol)
       ncol_A_or_b = this%nrhs
       lrow_A_or_b = this%lrowb
       lcol_A_or_b = this%lcolb
       desc_A_or_b = this%descb
       write(errstr,*) "will distribute incoming values to right-hand-side array of system, which has global dimensions ",&
            nrow_A_or_b," - by - ",ncol_A_or_b
       call add(errmsg,0,errstr,myname)
    case default
       call add(errmsg,2,"incoming type of array '"//trim(which_array)//&
            "' not supported, must be 'matrix' or 'rhs'",myname)
       return
    end select
!
    ! the distributing process checks incoming dimensions of distribution array depending on type
    if(this%myrow == sendprow .and. this%mycol == sendpcol) then
       write(errstr,*) "I am the process which distributes the values: (row,col) = ",sendprow,sendpcol
       call add(errmsg,0,errstr,myname)
       if(nrow_send < 1 .or. ncol_send < 1) then
          write(errstr,*) "incoming pretended size of distribution matrix = (",nrow_send,",",ncol_send,&
               ") must be valid (positive dimensions)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(mat_send,1) /= nrow_send .or. size(mat_send,2) /= ncol_send) then
          write(errstr,*) "pretended size of distribution matrix (",nrow_send,",",ncol_send,&
               ") does not match actual size of incoming array (",size(mat_send,1),",",size(mat_send,2),")"
          call add(errmsg,2,errstr,myname)
          return
       end if
       write(errstr,*) "dimensions of incoming subarray to distribute: ",&
            nrow_send," - by - ",ncol_send
       call add(errmsg,0,errstr,myname)
       if(ig_start < 1 .or. ig_start+nrow_send-1 > nrow_A_or_b) then
          write(errstr,*) "requested row-position of submatrix = ",ig_start," is out of permitted range: must be between 1 and ",&
               nrow_A_or_b-nrow_send," in order to fit in range of global dimensions"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(jg_start < 1 .or. jg_start+ncol_send-1 > ncol_A_or_b) then
          write(errstr,*) "requested column-position of submatrix = ",jg_start,&
               " is out of permitted range: must be between 1 and ",ncol_A_or_b-ncol_send,&
               " in order to fit in range of global dimensions"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if ! this%myrow == sendprow .and. this%mycol == sendpcol
!
!
    !! START TO DISTRIBUTE mat_send (send by distproc) TO ARRAY A_or_b (received by everyone)

    ! BLACS send subroutine:
    ! subroutine SGESD2D( ICONTXT,  M, N, A, LDA, RDEST, CDEST )

    ! BLACS receive subroutine:
    ! subroutine SGERV2D( ICONTXT,  M, N, A, LDA, RSRC, CSRC )
!
    ! figure out in which blocks (row, column blocks) the incoming global submatrix is divided by the SCALAPACK block-cyclic distribution techinque
    ! be aware that the first block and the last block can be a partial one!
    i_block_start = (ig_start-1) / this%nbrow + 1
    i_block_end = (ig_start+nrow_send-1 -1) / this%nbrow + 1
    m_block_send = i_block_end - i_block_start + 1  ! must be at least 1, since it was checked that nrow_send > 0
    if(m_block_send < 1) then
       call add(errmsg,2,"There are no row blocks to distribute. This error should actually not occur!! "//&
            "Please contact the developer",myname)
       return
    end if

    j_block_start = (jg_start-1) / this%nbcol + 1
    j_block_end = (jg_start+ncol_send-1 -1) / this%nbcol + 1
    n_block_send = j_block_end - j_block_start + 1 ! must be at least 1, since it was checked that ncol_send > 0
    if(n_block_send < 1) then
       call add(errmsg,2,&
            "There are no column blocks to distribute. This error should actually not occur!! "//&
            "Please contact the developer",myname)
       return
    end if


    ! THIS PIECE OF CODE NO LONGER NEEDED, since explicitely calling INFOG2L in nested loops below
    ! ! INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    ! ! curprow = INDXG2P(ig_start,this%nbrow,0,this%nprow) ! this intrinsic does the same computation as is done in the following line:
    ! curprow = mod(i_block_start-1,this%nprow) ! row in process grid, to which the first values will be sent
    ! curpcol = mod(j_block_start-1,this%nprow) ! column in process grid, to which the first values will be sent


    ! define arrays for each block to distribute, containing
    !  - number of rows and columns in submatrix to be sent/received per block
    !  - starting indices ig_start_block,jg_start_block in global array where block starts
    allocate(nrow_of_block(m_block_send),ncol_of_block(n_block_send),&
         ig_start_of_block(m_block_send),jg_start_of_block(n_block_send))
!
    nrow_send_tmp = nrow_send ! nrow_send_tmp contains the number of rows that are left to be distributed (decreases in following loop)
    ! EXCLUDE THE FIRST BLOCK FROM THE LOOP, DO IT MANUALLY
    ig_start_of_block(1) = ig_start
                           ! all rows until the next full block starts:
                           ! mod(ig_start-1,this%nbrow) is the number of rows PRECEEDING row ig_start inside the current block
                           ! simply substract this number from the total number of rows in a block to get the number of rows up to the next full block (including the first row ig_start
    nrow_of_block(1) = min(this%nbrow - mod(ig_start_of_block(1)-1,this%nbrow) , nrow_send_tmp)
    nrow_send_tmp = nrow_send_tmp - nrow_of_block(1)
    ! now loop on the rest of the row blocks
    do i_block = 2,m_block_send
       ig_start_of_block(i_block) = ig_start_of_block(i_block-1) + nrow_of_block(i_block-1)
       nrow_of_block(i_block) = min(this%nbrow , nrow_send_tmp)
       nrow_send_tmp = nrow_send_tmp - nrow_of_block(i_block)
    end do ! i_block
!
    ! now do the same for the colums
    ncol_send_tmp = ncol_send ! ncol_send_tmp contains the number of columns that are left to be distributed (decreases in following loop)
    ! EXCLUDE THE FIRST BLOCK FROM THE LOOP, DO IT MANUALLY
    jg_start_of_block(1) = jg_start
                           ! all cols until the next full block starts:
                           ! mod(jg_start-1,this%nbcol) is the number of cols PRECEEDING col jg_start inside the current block
                           ! simply substract this number from the total number of cols in a block to get the number of cols up to the next full block (including the first col jg_start
    ncol_of_block(1) = min(this%nbcol - mod(jg_start_of_block(1)-1,this%nbcol) , ncol_send_tmp)
    ncol_send_tmp = ncol_send_tmp - ncol_of_block(1)
    ! now loop on the rest of the col blocks
    do j_block = 2,n_block_send
       jg_start_of_block(j_block) = jg_start_of_block(j_block-1) + ncol_of_block(j_block-1)
       ncol_of_block(j_block) = min(this%nbcol , ncol_send_tmp)
       ncol_send_tmp = ncol_send_tmp - ncol_of_block(j_block)
    end do ! j_block
    

    ! loop on all the blocks (two nested loops on block rows and block columns) and set local arrays:
    !   if I'm sending proc: simiply set values in local array if I'm also receiving
    !                    else send values to correct receiving block
    !   else receive from sending proc if I'm the receiver
!
    do i_block = 1,m_block_send
       nr = nrow_of_block(i_block)
       ig = ig_start_of_block(i_block)
       i_mat_send = ig-ig_start+1 ! row index of array mat_send at which the currently sent rows start

       do j_block = 1,n_block_send
          nc = ncol_of_block(j_block)
          jg = jg_start_of_block(j_block)
          j_mat_send = jg-jg_start+1 ! column index of array mat_send at which the currently sent columns start

          ! ScaLAPAPACK tools routine to get local info of global matrix entry:
          !SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW, MYCOL, LRINDX, LCINDX, RSRC, CSRC )
          call INFOG2L(ig, jg, desc_A_or_b, this%nprow, this%npcol, this%myrow, this%mycol, il, jl, recvprow, recvpcol)

          if(this%myrow == sendprow .and. this%mycol == sendpcol) then
             if(this%myrow == recvprow .and. this%mycol == recvpcol) then
                ! if the sending proc (sendprow,sendpcol) also is the one holding the local subblock, 
                ! simply assign the values and do not send via BLACS (obviously)
                A_or_b(il:il+nr-1,jl:jl+nc-1) = mat_send(i_mat_send:i_mat_send+nr-1,j_mat_send:j_mat_send+nc-1)
             else
                ! otherwise send to receiving proc
                ! BLACS send subroutine:
                ! subroutine SGESD2D( ICONTXT,  M, N, A, LDA, RDEST, CDEST )
                call SGESD2D( this%blacs_context,  nr, nc, mat_send(i_mat_send:i_mat_send+nr-1,j_mat_send:j_mat_send+nc-1), &
                     nr, recvprow, recvpcol )
             end if ! I am receiver
          else ! I am sender
             if(this%myrow == recvprow .and. this%mycol == recvpcol) then
                ! receive block from sending proc (sendprow,sendpcol)
                ! BLACS receive subroutine:
                ! subroutine SGERV2D( ICONTXT,  M, N, A, LDA, RSRC, CSRC )
                call SGERV2D( this%blacs_context,  nr, nc, A_or_b(il:il+nr-1,jl:jl+nc-1), nr, sendprow, sendpcol )
             end if ! I am receiver
          end if ! I am sender
       end do ! j_block
    end do ! i_block
    
    deallocate(nrow_of_block,ig_start_of_block,ncol_of_block,jg_start_of_block)

  end subroutine distributeGlobalSubmatrixParallelLinearSystem
  !-----------------------------------------------------------------
  subroutine collectGlobalSubmatrixParallelLinearSystem(this,which_array,mat_recv,nrow_recv,ncol_recv,&
       ig_start,jg_start,recvprow,recvpcol,errmsg)
    ! incoming
    type (parallel_linear_system) :: this
    character(len=*) :: which_array
    real, dimension(:,:) :: mat_recv
    integer :: nrow_recv,ncol_recv,ig_start,jg_start
    integer :: recvprow,recvpcol
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=42) :: myname = 'collectGlobalSubmatrixParallelLinearSystem'
    character(len=400) :: errstr
    ! local array to add values to, fork to either A or b
    real, dimension(:,:), pointer :: A_or_b
    integer, dimension(9) :: desc_A_or_b
    integer :: nrow_A_or_b,ncol_A_or_b,lrow_A_or_b,lcol_A_or_b
    ! regarding rows
    integer :: nrow_recv_tmp,i_block,i_block_start,i_block_end,m_block_recv
    integer :: ig,il,nr,i_mat_recv,sendprow
    integer, dimension(:), allocatable :: nrow_of_block,ig_start_of_block
    ! regarding columns
    integer :: ncol_recv_tmp,j_block,j_block_start,j_block_end,n_block_recv
    integer :: jg,jl,nc,j_mat_recv,sendpcol
    integer, dimension(:), allocatable :: ncol_of_block,jg_start_of_block
!
    ! IT IS ASSUMED THAT ALL PROCESSES OF THE PROCESS GRID GET THE SAME VALUES OF 
    !   which_array,nrow_recv,ncol_recv,ig_start,jg_start,recvprow,recvpcol !!
!
    call addTrace(errmsg,myname)
!
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
    if(.not.this%initiated) then
       call add(errmsg,2,"parallel linear system not yet initiated",myname)
       return
    end if
    if(.not.associated(this%LA)) then
       call add(errmsg,2,"parallel linear system not yet allocated",myname)
       return
    end if
    if(recvprow < 0 .or. recvprow > this%nprow .or. recvpcol < 0 .or. recvpcol > this%npcol) then
       write(errstr,*) "requested process of collection is (row,col) = ",recvprow,recvpcol,&
            "; however, maximum dimension of process grid is (nrow,ncol) = ",this%nprow,this%npcol
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! check incoming type of values (matrix or rhs) to collect
    select case(which_array)
    case ('matrix')
       A_or_b => this%LA
       nrow_A_or_b = this%nrow
       ncol_A_or_b = this%ncol
       lrow_A_or_b = this%lrowA
       lcol_A_or_b = this%lcolA
       desc_A_or_b = this%descA
       write(errstr,*) "will collect system submatrix, which has global dimensions ",&
            nrow_A_or_b," - by - ",ncol_A_or_b
       call add(errmsg,0,errstr,myname)
    case ('rhs')
       A_or_b => this%Lb
       nrow_A_or_b = max(this%nrow,this%ncol)
       ncol_A_or_b = this%nrhs
       lrow_A_or_b = this%lrowb
       lcol_A_or_b = this%lcolb
       desc_A_or_b = this%descb
       write(errstr,*) "will collect right-hand-side subarray of system, which has global dimensions ",&
            nrow_A_or_b," - by - ",ncol_A_or_b
       call add(errmsg,0,errstr,myname)
    case default
       call add(errmsg,2,"incoming type of array '"//trim(which_array)//&
            "' not supported, must be 'matrix' or 'rhs'",myname)
       return
    end select
!
    ! the collecting process checks incoming dimensions of the collection array depending on type
    if(this%myrow == recvprow .and. this%mycol == recvpcol) then
       write(errstr,*) "I am the collecting process: (row,col) = ",recvprow,recvpcol
       call add(errmsg,0,errstr,myname)
       if(nrow_recv < 1 .or. ncol_recv < 1) then
          write(errstr,*) "incoming pretended size of collection array = (",nrow_recv,",",ncol_recv,&
               ") must be valid (positive dimensions)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(mat_recv,1) /= nrow_recv .or. size(mat_recv,2) /= ncol_recv) then
          write(errstr,*) "pretended size of collection array (",nrow_recv,",",ncol_recv,&
               ") does not match actual size of incoming array (",size(mat_recv,1),",",size(mat_recv,2),")"
          call add(errmsg,2,errstr,myname)
          return
       end if
       write(errstr,*) "dimensions of incoming collection array: ",&
            nrow_recv," - by - ",ncol_recv
       call add(errmsg,0,errstr,myname)
       if(ig_start < 1 .or. ig_start+nrow_recv-1 > nrow_A_or_b) then
          write(errstr,*) "requested row-position of submatrix to collect = ",ig_start,&
               " is out of permitted range: must be between 1 and ",&
               nrow_A_or_b-nrow_recv," in order to fit in range of global dimensions"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(jg_start < 1 .or. jg_start+ncol_recv-1 > ncol_A_or_b) then
          write(errstr,*) "requested column-position of submatrix to collect = ",jg_start,&
               " is out of permitted range: must be between 1 and ",ncol_A_or_b-ncol_recv,&
               " in order to fit in range of global dimensions"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if ! this%myrow == recvprow .and. this%mycol == recvpcol
!
!
    !! START TO COLLECT mat_recv (send to recvproc) FROM ARRAY A_or_b (sent by everyone)

    ! BLACS send subroutine:
    ! subroutine SGESD2D( ICONTXT,  M, N, A, LDA, RDEST, CDEST )

    ! BLACS receive subroutine:
    ! subroutine SGERV2D( ICONTXT,  M, N, A, LDA, RSRC, CSRC )
!
    ! figure out in which blocks (row, column blocks) the incoming global submatrix is divided by the SCALAPACK block-cyclic distribution techinque
    ! be aware that the first block and the last block can be a partial one!
    i_block_start = (ig_start-1) / this%nbrow + 1
    i_block_end = (ig_start+nrow_recv-1 -1) / this%nbrow + 1
    m_block_recv = i_block_end - i_block_start + 1  ! must be at least 1, since it was checked that nrow_recv > 0
    if(m_block_recv < 1) then
       call add(errmsg,2,"There are no row blocks to collect. This error should actually not occur!! "//&
            "Please contact the developer",myname)
       return
    end if

    j_block_start = (jg_start-1) / this%nbcol + 1
    j_block_end = (jg_start+ncol_recv-1 -1) / this%nbcol + 1
    n_block_recv = j_block_end - j_block_start + 1 ! must be at least 1, since it was checked that ncol_recv > 0
    if(n_block_recv < 1) then
       call add(errmsg,2,&
            "There are no column blocks to collect. This error should actually not occur!! "//&
            "Please contact the developer",myname)
       return
    end if


    ! THIS PIECE OF CODE NO LONGER NEEDED, since explicitely calling INFOG2L in nested loops below
    ! ! INTEGER FUNCTION INDXG2P( INDXGLOB, NB, IPROC, ISRCPROC, NPROCS )
    ! ! curprow = INDXG2P(ig_start,this%nbrow,0,this%nprow) ! this intrinsic does the same computation as is done in the following line:
    ! curprow = mod(i_block_start-1,this%nprow) ! row in process grid, to which the first values will be sent
    ! curpcol = mod(j_block_start-1,this%nprow) ! column in process grid, to which the first values will be sent


    ! define arrays for each block to collect, containing
    !  - number of rows and columns in submatrix to be sent/received per block
    !  - starting indices ig_start_block,jg_start_block in global array where block starts
    allocate(nrow_of_block(m_block_recv),ncol_of_block(n_block_recv),&
         ig_start_of_block(m_block_recv),jg_start_of_block(n_block_recv))
!
    nrow_recv_tmp = nrow_recv ! nrow_recv_tmp contains the number of rows that are left to be collected (decreases in following loop)
    ! EXCLUDE THE FIRST BLOCK FROM THE LOOP, DO IT MANUALLY
    ig_start_of_block(1) = ig_start
                           ! all rows until the next full block starts:
                           ! mod(ig_start-1,this%nbrow) is the number of rows PRECEEDING row ig_start inside the current block
                           ! simply substract this number from the total number of rows in a block to get the number of rows up to the next full block (including the first row ig_start
    nrow_of_block(1) = min(this%nbrow - mod(ig_start_of_block(1)-1,this%nbrow) , nrow_recv_tmp)
    nrow_recv_tmp = nrow_recv_tmp - nrow_of_block(1)
    ! now loop on the rest of the row blocks
    do i_block = 2,m_block_recv
       ig_start_of_block(i_block) = ig_start_of_block(i_block-1) + nrow_of_block(i_block-1)
       nrow_of_block(i_block) = min(this%nbrow , nrow_recv_tmp)
       nrow_recv_tmp = nrow_recv_tmp - nrow_of_block(i_block)
    end do ! i_block
!
    ! now do the same for the colums
    ncol_recv_tmp = ncol_recv ! ncol_recv_tmp contains the number of columns that are left to be collected (decreases in following loop)
    ! EXCLUDE THE FIRST BLOCK FROM THE LOOP, DO IT MANUALLY
    jg_start_of_block(1) = jg_start
                           ! all cols until the next full block starts:
                           ! mod(jg_start-1,this%nbcol) is the number of cols PRECEEDING col jg_start inside the current block
                           ! simply substract this number from the total number of cols in a block to get the number of cols up to the next full block (including the first col jg_start
    ncol_of_block(1) = min(this%nbcol - mod(jg_start_of_block(1)-1,this%nbcol) , ncol_recv_tmp)
    ncol_recv_tmp = ncol_recv_tmp - ncol_of_block(1)
    ! now loop on the rest of the col blocks
    do j_block = 2,n_block_recv
       jg_start_of_block(j_block) = jg_start_of_block(j_block-1) + ncol_of_block(j_block-1)
       ncol_of_block(j_block) = min(this%nbcol , ncol_recv_tmp)
       ncol_recv_tmp = ncol_recv_tmp - ncol_of_block(j_block)
    end do ! j_block
    

    ! loop on all the blocks (two nested loops on block rows and block columns) and set local arrays:
    !   if I'm receiving proc: simiply set values in return array mat_recv if I'm also sending
    !                    else receive values from correct sending block
    !   else send to receiving proc if I'm the sender
!
    do i_block = 1,m_block_recv
       nr = nrow_of_block(i_block)
       ig = ig_start_of_block(i_block)
       i_mat_recv = ig-ig_start+1 ! row index of array mat_recv at which the currently sent rows start

       do j_block = 1,n_block_recv
          nc = ncol_of_block(j_block)
          jg = jg_start_of_block(j_block)
          j_mat_recv = jg-jg_start+1 ! column index of array mat_recv at which the currently sent columns start

          ! ScaLAPAPACK tools routine to get local info of global matrix entry:
          !SUBROUTINE INFOG2L( GRINDX, GCINDX, DESC, NPROW, NPCOL, MYROW, MYCOL, LRINDX, LCINDX, RSRC, CSRC )
          call INFOG2L(ig, jg, desc_A_or_b, this%nprow, this%npcol, this%myrow, this%mycol, il, jl, sendprow, sendpcol)

          if(this%myrow == recvprow .and. this%mycol == recvpcol) then
             if(this%myrow == sendprow .and. this%mycol == sendpcol) then
                ! if the receiving proc (recvprow,recvpcol) also is the one holding the local subblock, 
                ! simply assign the values and do not receive via BLACS (obviously)
                mat_recv(i_mat_recv:i_mat_recv+nr-1,j_mat_recv:j_mat_recv+nc-1) = A_or_b(il:il+nr-1,jl:jl+nc-1)
             else
                ! otherwise receive from sending proc
                ! BLACS receive subroutine:
                ! subroutine SGERV2D( ICONTXT,  M, N, A, LDA, RSRC, CSRC )
                call SGERV2D( this%blacs_context,  nr, nc, mat_recv(i_mat_recv:i_mat_recv+nr-1,j_mat_recv:j_mat_recv+nc-1), &
                     nr, sendprow, sendpcol )
             end if ! I am sender
          else ! I am receiver
             if(this%myrow == sendprow .and. this%mycol == sendpcol) then
                ! send block to receiving proc (recvprow,recvpcol)
                ! BLACS send subroutine:
                ! subroutine SGESD2D( ICONTXT,  M, N, A, LDA, RDEST, CDEST )
                call SGESD2D( this%blacs_context,  nr, nc, A_or_b(il:il+nr-1,jl:jl+nc-1), &
                     nr, recvprow, recvpcol )
             end if ! I am sender
          end if ! I am receiver
       end do ! j_block
    end do ! i_block
    
    deallocate(nrow_of_block,ig_start_of_block,ncol_of_block,jg_start_of_block)

  end subroutine collectGlobalSubmatrixParallelLinearSystem
  !-----------------------------------------------------------------
  subroutine solveLeastSquaresParallelLinearSystem(this,errmsg)
    type (parallel_linear_system) :: this
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=37) :: myname = 'solveLeastSquaresParallelLinearSystem'
    integer :: LWORK,INFO
    real, dimension(:), allocatable :: WORK
!
    call addTrace(errmsg,myname)
!
    if(.not.this%im_in_the_grid) return
    if(.not.this%initiated) then
       call add(errmsg,2,"parallel linear system not yet initiated",myname)
       return
    end if
    if(.not.associated(this%LA)) then
       call add(errmsg,2,"parallel linear system not yet allocated",myname)
       return
    end if
!
! NOW SOLVE THE SYSTEM IN A LEAST-SQUARES-SENSE BY SCALAPACK
! SUBROUTINE PSGELS( TRANS, M, N, NRHS, A, IA, JA, DESCA, B, IB, JB, DESCB, WORK, LWORK, INFO )
!
    ! do a workspace query first
    allocate(WORK(1)); LWORK = -1
    call PSGELS( 'N', this%nrow, this%ncol, this%nrhs, this%LA, 1, 1, this%descA, this%Lb, 1, 1, this%descb, &
         WORK, LWORK, INFO )
    if(INFO /= 0 .or. WORK(1) <= 0) then
       write(errstr,*) "workspace query was not successful, ScaLAPACK routine PSGELS returned INFO = ",&
            INFO,"; optimal size WORK(1) = ",WORK(1)
       call add(errmsg,2,errstr,myname)
       goto 1
    else
       write(errstr,*) "workspace query was successful, ScaLAPACK routine PSGELS returned INFO = ",&
            INFO,"; optimal size WORK(1) = ",WORK(1)
       call add(errmsg,0,errstr,myname)
    end if
    LWORK = WORK(1)
    deallocate(WORK)
    allocate(WORK(LWORK))
!
    ! now actually solve the system
    call PSGELS( 'N', this%nrow, this%ncol, this%nrhs, this%LA, 1, 1, this%descA, this%Lb, 1, 1, this%descb, &
         WORK, LWORK, INFO )
    if(INFO /= 0) then
       write(errstr,*) "solving the linear system was not successful, ScaLAPACK routine PSGELS returned INFO = ",&
            INFO
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
1   if(allocated(WORK)) deallocate(WORK)
  end subroutine solveLeastSquaresParallelLinearSystem
  !-----------------------------------------------------------------
  subroutine solveQRParallelLinearSystem(this,errmsg)
    type (parallel_linear_system) :: this
    type (error_message) :: errmsg
    ! local
    !character(len=400) :: errstr
    character (len=37) :: myname = 'solveLeastSquaresParallelLinearSystem'
    !integer :: LWORK,INFO
    !real, dimension(:), allocatable :: WORK
!
    call addTrace(errmsg,myname)
!
    call add(errmsg,2,"this routine is not yet fully implemented",myname)
!
!            CALL PSORMQR( 'Left', 'Transpose', M, NRHS, N, A, IA, JA,
!     $                    DESCA, WORK, B, IB, JB, DESCB, WORK( IPW ),
!     $                    LWORK-LTAU, INFO )
!*
!*           workspace at least NRHS, optimally NRHS*NB
!*
!*           B(IB:IB+N-1,JB:JB+NRHS-1) := inv(R) *
!*                                        B(IB:IB+N-1,JB:JB+NRHS-1)
!*
!            CALL PSTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
!     $                   NRHS, ONE, A, IA, JA, DESCA, B, IB, JB, DESCB )

  end subroutine solveQRParallelLinearSystem
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  function getMyrowParallelLinearSystem(this) result(myrow)
    type (parallel_linear_system), intent(in) :: this
    integer :: myrow
    myrow = this%myrow
  end function getMyrowParallelLinearSystem
  !-----------------------------------------------------------------
  function getMycolParallelLinearSystem(this) result(mycol)
    type (parallel_linear_system), intent(in) :: this
    integer :: mycol
    mycol = this%mycol
  end function getMycolParallelLinearSystem
  !-----------------------------------------------------------------
  function getImInTheGridParallelLinearSystem(this) result(l)
    type (parallel_linear_system), intent(in) :: this
    logical :: l
    l = this%im_in_the_grid
  end function getImInTheGridParallelLinearSystem




!###############################################################################################
!## THE FOLLOWING ROUTINES ARE STILL REQUIRED TO EXIST (could also be empty, though), SINCE 
!## THEY ARE FORKED TO BY GENERIC MODULE linearSystem, WHICH HOWEVER IS DEPRECIATED!!
!## LEAVE THE ROUTINES BE FOR COMPATIBILITY WITH OLD PROGRAMS USING linearSystem
!###############################################################################################

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
!###################
!## THIS IS AN OLD PRELIMINARY VERSION TESTING THE MODULE
!###################
    this%NPROW=NPROW; this%NPCOL=NPCOL
    this%NBROW=NBROW; this%NBCOL=NBCOL
    print *,"vor SL_INIT, ICTXT, NPROW, NPCOL=",ICTXT, NPROW, NPCOL
    call SL_INIT(ICTXT,NPROW,NPCOL)
    print *,"nach SL_INIT, ICTXT, NPROW, NPCOL=",ICTXT, NPROW, NPCOL
    this%blacs_context=ICTXT
    print *,"createParallelLinearSystem,ICTXT=",this%blacs_context
    if(ICTXT<0) return
    call BLACS_GRIDINFO(ICTXT,NPROW,NPCOL,MYROW,MYCOL)
    this%MYROW=MYROW; this%MYCOL=MYCOL
    print *,"MYROW,MYCOL=",MYROW,MYCOL,";  ICTXT=",ICTXT,";  NPROW,NPCOL,NBROW,NBCOL=",NPROW,NPCOL,NBROW,NBCOL
    !Initialize the array descriptors for the matrices A and b
    !call DESCINIT( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT, MXLLDA, INFO )
    !call DESCINIT( DESCB, N, NRHS, NB, NBRHS, RSRC, CSRC, ICTXT, MXLLDB, INFO )
  end subroutine createParallelLinearSystem
  !-----------------------------------------------------------------
  !> \brief solve linear system by least squares
  !! \details Calling SCALAPACK routine (??) which ...
  !! \param this parallel_linear_system object
  !! \param errmsg error_message object
  !! \return error message
  !! \pre You need SCALAPACK libraries (including BLACS etc.).
  !
  function solveLeastSquaresFunctionParallelLinearSystem(this) result(errmsg)
    type (parallel_linear_system) :: this
    type (error_message) :: errmsg
    character (len=45) :: myname = 'solveLeastSquaresFunctionParallelLinearSystem'
    call new(errmsg,myname)
    call new(errmsg,1,'this is a warning message from this fake routine!!',myname)
  end function solveLeastSquaresFunctionParallelLinearSystem
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
!###############################################################################################
!###############################################################################################


end module parallelLinearSystem
