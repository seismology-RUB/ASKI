!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.2 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief set up and handle parallelized kernel matrix system
!!
!! \details Given a data_sample_info object and general prequesits
!!  of the inversion and iteration step, the kernel matrix is set up
!!  in a parallel environment using BLACS libraries
!!  and solved using module parallelLinearSystem. The solution is collected by
!!  the master process from the process grid. 
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module parKernelLinearSystem
!
  use kernelLinearSystem
  use parallelLinearSystem
  use dataModelSpaceInfo
  use seismicEvent
  use seismicStation
  use linearModelRegularization
  use modelParametrization
  use vectorPointer
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateParKernelLinearSystem; end interface
  interface operator (.nrow.); module procedure getNrowParKernelLinearSystem; end interface
  interface operator (.ncol.); module procedure getNcolParKernelLinearSystem; end interface
!
  type par_kernel_linear_system
     private
     ! KERNEL SYSTEM SPECIFICATIONS
     type (data_model_space_info) :: dmspace !< data and model space of GLOBAL kernel system
     integer :: nrow = 0 !< global number of rows of kernel matrix, nrow=ndata+nsmooth
     integer :: ndata = 0 !< global number of rows reserved for kernels or synthetic data (first rows of the system)
     integer :: nrowreg = 0 !< global number rows for regularization conditions (rows) that were added to the system (last rows of the system)
     integer :: ncol = 0 !< global number of columns of kernel matrix, ncol = nmval+ncolreg
     integer :: nmval = 0 !< global number of columns reserved for kernel values on inversion grid cells
     integer :: ncolreg = 0 !< global number of columns for column regularization (e.g. in kernel focussing, last columns of the matrix)
     integer :: nrhs = 0 !< global number of right-hand-side vectors

     ! BLACS GRID, PROCS AND BLOCKS
     integer :: nprow = 0
     integer :: npcol = 0
     integer :: nbrow = 0
     integer :: nbcol = 0
     integer :: nrow_handle = 0
     ! PARALLEL LINEAR SYSTEM OBJECT
     type (parallel_linear_system) :: PLSE
     integer :: myrow = -1
     integer :: mycol = -1
     logical :: im_in_the_grid = .false.
     ! INFORMATION OF GLOBAL KERNEL MATRIX
     real, dimension(:,:), pointer :: minmaxval_K_per_param => null() !< array of size 2-by-nparam_pmtrz, only allocated, defined and used by master process
     logical :: kernel_matrix_set = .false.
     logical :: row_regularization_set = .false.
     logical :: column_regularization_set = .false.
     logical :: rhs_set = .false.
     ! LOCAL KERNEL MATRIX
     real, dimension(:,:), pointer :: LK => null() !< local kernel matrix (only pointer to PLSE%LA !! do not deallocate)
     ! LOCAL RHS ARRAY
     real, dimension(:,:), pointer :: Lrhs => null() !< local rhs array (only pointer to PLSE%Lb !! do not deallocate)
  end type par_kernel_linear_system
!
contains
!
  subroutine initiateParKernelLinearSystem(this,dmspace,nrowreg,ncolreg,nrhs,&
       nbrow,nbcol,nprow,npcol,nrow_handle,im_in_the_grid,errmsg)
    ! incoming
    type (par_kernel_linear_system) :: this
    type (data_model_space_info) :: dmspace
    integer :: nrowreg,ncolreg,nrhs
    integer :: nbrow,nbcol,nprow,npcol,nrow_handle
    ! returning
    logical :: im_in_the_grid
    type (error_message) :: errmsg
    ! local
    character (len=29) :: myname = 'initiateParKernelLinearSystem'
    character (len=400) :: errstr
    integer :: ndata,nmval
!
    call addTrace(errmsg,myname)
!
    ndata = .ndata.dmspace
    nmval = .nmval.dmspace
    if(any((/ndata,nmval,nrhs/) <= 0)) then
       write(errstr,*) "some of the incoming values are invalid, must be > 0: number of data = ",&
            ndata,", number of model values = ",nmval,", number of right-hand-sides of the system = ",&
            nrhs
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any((/nbrow,nbcol,nprow,npcol,nrow_handle/) <= 0)) then
       write(errstr,*) "some of the incoming parallelization parameters are invalid, must be > 0: ",&
            "number of rows,columns in ScaLAPACK submatrix block = ",nbrow,", ",nbcol,&
            "; number of processes in rows,columns of BLACS process grid = ",npcol,", ",npcol,&
            "; number of rows to handle by memory of master process = ",nrow_handle
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    if(any((/nrowreg,ncolreg/) < 0)) then
       write(errstr,*) "some of the incoming regularization parameters are invalid, must be >= 0: ",&
            "number of smoothing/damping conditions = ",nrowreg,", number of column regularization = ",ncolreg
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    this%ndata = ndata
    this%nrowreg = nrowreg
    if(this%nrowreg == 0) then
       this%row_regularization_set = .true.
    else
       this%row_regularization_set = .false.
    end if
    this%nrow = this%ndata + this%nrowreg

    this%nmval = nmval
    this%ncolreg = ncolreg
    if(this%ncolreg == 0) then
       this%column_regularization_set = .true.
    else
       this%column_regularization_set = .false.
    end if
    this%ncol = this%nmval  + this%ncolreg

    this%nrhs = nrhs

    this%nbrow = nbrow
    this%nbcol = nbcol

    this%nprow = nprow
    this%npcol = npcol

    this%nrow_handle = nrow_handle
!
    this%kernel_matrix_set = .false.
    this%rhs_set = .false.
!
    ! subroutine initParallelLinearSystem(this,nrow,ncol,nrhs,nprow,npcol,nbrow,nbcol,errmsg)
    call initParallelLinearSystem(this%PLSE,this%nrow,this%ncol,this%nrhs,&
         this%nprow,this%npcol,this%nbrow,this%nbcol,errmsg)
    if(.level.errmsg == 2) then
       call add(errmsg,2,"initiation of parllel linear system failed",myname)
       return
    end if
!
    this%im_in_the_grid = im_in(this%PLSE)
    this%myrow = .myrow.(this%PLSE)
    this%mycol = .mycol.(this%PLSE)
!
    call copyDataModelSpaceInfo(this%dmspace,dmspace)
!
    im_in_the_grid = this%im_in_the_grid
!
    return
!
1   continue
    !call deallocateParKernelLinearSystem(this)
  end subroutine initiateParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine allocateParKernelLinearSystem(this,errmsg)
    type (par_kernel_linear_system) :: this
    type (error_message) :: errmsg
    ! subroutine allocateLocalStorageParallelLinearSystem(this,A,b,errmsg)
    call allocateLocalStorageParallelLinearSystem(this%PLSE,this%LK,this%Lrhs,errmsg)
  end subroutine allocateParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine readMatrixParKernelLinearSystem(this,df_measured_data,&
       nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
       lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
       ignore_data_weights,apply_kernel_normalization)
    ! incoming
    type (par_kernel_linear_system) :: this
    character(len=*) :: path_sensitivity_kernels
    real :: df_measured_data
    integer, dimension(:) ::  ifreq_measured_data
    integer :: nfreq_measured_data,ntot_invgrid,lu1,lu2
    type (parameter_correlation) :: pcorr
    character(len=*), optional :: path_event_filter,path_station_filter
    logical, optional :: apply_event_filter,apply_station_filter,ignore_data_weights,apply_kernel_normalization
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=400)  :: errstr
    character(len=31) :: myname = 'readMatrixParKernelLinearSystem'
    type (error_message) :: errmsg2
    ! serial kernel linear system
    type (kernel_linear_system) :: KLSE
    real, dimension(:,:), pointer :: K_KLSE
    type (data_model_space_info) :: dmspace2
    ! scaling for model smoothing
    character(len=character_length_param) :: param_name
    integer :: iparam_pmtrz,nparam_pmtrz
    character(len=character_length_param), dimension(:), pointer :: pparam
    type (integer_vector_pointer), dimension(:), allocatable :: idx_mspace
    integer, dimension(:), pointer :: idx
    real :: minval_K_rowblock,maxval_K_rowblock
    ! BLACS grid stuff
    integer :: ig_start
    integer :: isection,nrow_section
!! FS FS DEBUGGING START
!!$integer :: i,j
!!$character(len=100) :: filename
!!$real, dimension(:,:), pointer :: mat
!! FS FS DEBUGGING END
!
    nullify(K_KLSE,pparam,idx)
!
    call addTrace(errmsg,myname)
    ! if there was no call to initiateParKernelLinearSystem before, then always this%im_in_the_grid == .false.
    ! so the following check should be sufficient to test whether this par_kernel_linear_system was initiated (e.g. this%dmspace is set etc.)
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
!
    if(.not.associated(this%LK)) then
       call add(errmsg,2,"parallel kernel linear system not allocated yet",myname)
       goto 1
    end if
!
    K_KLSE => null()
!
    if(this%myrow == 0 .and. this%mycol == 0) then
       nparam_pmtrz = numberOfParamModelParametrization(.pmtrz.this%dmspace)
       if(nparam_pmtrz <= 0) then
          call add(errmsg,2,"invalid model parametrization contained in data model space info object",myname)
          goto 1
       end if
       allocate(this%minmaxval_K_per_param(2,nparam_pmtrz),idx_mspace(nparam_pmtrz))
       this%minmaxval_K_per_param = 0.
       ! prepare for calculating global max and min values of kernel matrix 
       nullify(pparam)
       do while (nextParamModelParametrization(.pmtrz.this%dmspace,param_name,iparam_pmtrz))
          if(associated(pparam)) deallocate(pparam)
          allocate(pparam(1)); pparam(1) = param_name
          idx => getIndxModelValues(this%dmspace,param=pparam)
          call associateVectorPointer(idx_mspace(iparam_pmtrz),idx)
          nullify(idx)
       end do ! while nextParamModelPmtrz
    end if ! this%myrow == this%mycol == 0
!
    isection = 0
!
    write(errstr,*) "distributing global matrix in sections of rows containing ",this%nrow_handle," rows each"
    call add(errmsg2,0,errstr,myname)
!
    ! loop over subset of rows (always use a section of this%nrow_handle rows)
    do ig_start = 1,this%ndata,this%nrow_handle ! ig_start is the global index of the first row of the current submatrix
       isection = isection +1
       call new(errmsg2,myname)

       nrow_section = min(this%nrow_handle,this%ndata-ig_start+1) ! number of rows in the current section of rows

       write(errstr,*) "current section of rows: ",isection,"; number of rows in it: ",nrow_section
       call add(errmsg2,0,errstr,myname)

       ! extract the current section of rows from data model space info object:
       call copyDataModelSpaceInfo(dmspace2,this%dmspace,idata1=ig_start,idata2=ig_start+nrow_section-1)
       if(.level.errmsg2 == 2) then
          call print(errmsg2)
          call add(errmsg,2,"could not extract data samples for current section of rows from "//&
               "data model space info object",myname)
          goto 1
       end if

       if(this%myrow == 0 .and. this%mycol == 0) then
          ! read in complete section of rows of global kernel matrix:
          write(*,*) "read and send parallel kernel matrix, ",isection,"'th section of rows out of approx. ",&
               (this%ndata/this%nrow_handle)+1,"; first row of this section = ",ig_start,&
               " ( handling ",this%nrow_handle," rows at once ) "

          ! subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
          call initiateSerialKernelLinearSystem(KLSE,dmspace2,0,0,errmsg2)
          if(.level.errmsg2 ==2) then
             call print(errmsg2)
             call add(errmsg,2,"could not initiate serial kernel system",myname)
             goto 1
          end if

          ! subroutine readMatrixSerialKernelLinearSystem(this,df_measured_data,&
          !      nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
          !      lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
          !      ignore_data_weights,apply_kernel_normalization)
          call readMatrixSerialKernelLinearSystem(KLSE,df_measured_data,&
               nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
               lu1,lu2,errmsg2,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
               ignore_data_weights,apply_kernel_normalization)
          if(.level.errmsg2 == 2) goto 1

          ! get serial kernel matrix
          K_KLSE => .KM.KLSE
          
          ! memorize (update) the max and min values of the total kernel matrix per model parameter
          ! for this reason loop on all model parameters contained in the model space and separately
          ! update min and max values
          do iparam_pmtrz = 1,nparam_pmtrz
             idx => getVectorPointer(idx_mspace(iparam_pmtrz))
             if(.not.associated(idx)) cycle

             minval_K_rowblock = minval(K_KLSE(:,idx))
             maxval_K_rowblock = maxval(K_KLSE(:,idx))
             if(ig_start == 1) then
                this%minmaxval_K_per_param(1,iparam_pmtrz) = minval_K_rowblock
                this%minmaxval_K_per_param(2,iparam_pmtrz) = maxval_K_rowblock
             else
                this%minmaxval_K_per_param(1,iparam_pmtrz) = min(this%minmaxval_K_per_param(1,iparam_pmtrz),minval_K_rowblock)
                this%minmaxval_K_per_param(2,iparam_pmtrz) = max(this%minmaxval_K_per_param(2,iparam_pmtrz),maxval_K_rowblock)
             end if
          end do ! iparam_pmtrz

!! FS FS DEBUGGING START
!!$do i = ig_start,ig_start+nrow_section-1
!!$   do j = 1,size(K_KLSE,2)
!!$      K_KLSE(i-ig_start+1,j) = i*1000 + j
!!$   end do ! j
!!$end do ! i
!! FS FS DEBUGGING END
       end if ! this%myrow == this%mycol == 0

       ! distribute section of rows on the blacs process grid
!
       !subroutine distributeGlobalSubmatrixParallelLinearSystem(this,which_array,mat_send,nrow_send,ncol_send,&
       !   ig_start,jg_start,sendprow,sendpcol,errmsg)
       call distributeGlobalSubmatrixParallelLinearSystem(this%PLSE,'matrix',K_KLSE,nrow_section,.nmval.dmspace2,&
            ig_start,1,0,0,errmsg2)
       if(.level.errmsg2 == 2) then
          call print(errmsg2)
          write(errstr,*) isection,"'th section of rows of global kernel matrix could not be distributed on the process grid"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       call dealloc(dmspace2)
       call dealloc(errmsg2)

       if(this%myrow == 0 .and. this%mycol == 0) call dealloc(KLSE)
    end  do ! ig_start

!! FS FS DEBUGGING START
!!$write(filename,"('DEBUG_local_A_',I2.2,'_'I2.2,'.dat')") this%myrow,this%mycol
!!$open(unit=20,file=filename,form='formatted',status='unknown',action='write')
!!$write(20,*) size(this%PLSE%LA,1),size(this%PLSE%LA,2)
!!$if(size(this%PLSE%LA,1) > 0 .and. size(this%PLSE%LA,2)>0 ) then
!!$   do i = 1,size(this%PLSE%LA,1)
!!$      write(20,*) this%PLSE%LA(i,:)
!!$   end do
!!$end if
!!$close(20)
!!$call sleep(10)
!!$stop
!! FS FS DEBUGGING END

!! FS FS DEBUGGING START
!!$if(this%myrow == 0 .and. this%mycol == 0) then
!!$   allocate(mat(this%nrow,this%ncol))
!!$else
!!$   nullify(mat)
!!$end if
!!$call collectGlobalSubmatrixParallelLinearSystem(this%PLSE,'matrix',mat,this%nrow,this%ncol,&
!!$     1,1,0,0,errmsg2)
!!$if(.level.errmsg2 == 2)then
!!$   call print(errmsg2)
!!$   stop
!!$end if
!!$if(this%myrow == 0 .and. this%mycol == 0) then
!!$   write(filename,"('DEBUG_global_A_',I2.2,'_'I2.2,'.dat')") this%myrow,this%mycol
!!$   open(unit=20,file=filename,form='formatted',status='unknown',action='write')
!!$   write(20,*) size(mat,1),size(mat,2)
!!$   if(size(mat,1) > 0 .and. size(mat,2)>0 ) then
!!$      do i = 1,size(mat,1)
!!$         write(20,*) mat(i,:)
!!$      end do
!!$   end if
!!$   close(20)   
!!$   deallocate(mat)
!!$end if
!!$call sleep(10)
!!$stop
!! FS FS DEBUGGING END

    ! if everything went OK, indicate by flag
    this%kernel_matrix_set = .true.

    ! CLEAN UP
1   if(allocated(idx_mspace)) then
       do iparam_pmtrz = 1,size(idx_mspace)
          call dealloc(idx_mspace(iparam_pmtrz))
       end do
       deallocate(idx_mspace)
    end if
    call dealloc(errmsg2)
    call dealloc(dmspace2)
    if(this%myrow == 0 .and. this%mycol == 0) call dealloc(KLSE)
  end subroutine readMatrixParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine addLinearModelRegularizationParKernelLinearSystem(this,lmreg,errmsg)
    type (par_kernel_linear_system) :: this
    type (linear_model_regularization) :: lmreg
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    type (error_message) :: errmsg2
    character(len=49) :: myname = 'addLinearModelRegularizationParKernelLinearSystem'
    ! row regularization conditions
    character(len=character_length_regscal_type) :: scltyp
    logical :: require_minmaxval_K_per_param
    real, dimension(:,:), pointer :: SA
    real, dimension(:), pointer :: Sb
    real, dimension(:,:), allocatable :: Srhs
    ! other
    integer :: i_rowreg,ig_start,isection,nrow_section,irhs
!! FS FS DEBUGGING START
!!$integer :: i
!!$character(len=100) :: filename
!!$real, dimension(:,:), pointer :: mat
!! FS FS DEBUGGING END
!
    nullify(SA,Sb)
!
    call addTrace(errmsg,myname)
    ! if there was no call to initiateParKernelLinearSystem before, then always this%im_in_the_grid == .false.
    ! so the following check should be sufficient to test whether this par_kernel_linear_system was initiated (e.g. this%dmspace is set etc.)
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
!
    if(.not.associated(this%LK)) then
       call add(errmsg,2,"parallel kernel linear system not allocated yet",myname)
       goto 1
    end if
!
    if(.neq.lmreg /= this%nrowreg) then
       write(errstr,*) "number of incoming row regularization constraints (smoothing/damping) = ",.neq.lmreg,&
            " does not equal number of smoothinc conditions for which parallel kernel system was allocated = ",&
            this%nrowreg
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    require_minmaxval_K_per_param = lmreg.requires.'minmaxval_K_per_param'
    if(require_minmaxval_K_per_param .and. this%myrow == 0 .and. this%mycol == 0) then
       if(.not.(associated(this%minmaxval_K_per_param).and.this%kernel_matrix_set)) then
          scltyp = getScalingTypeLinearModelRegularization(lmreg)
          call add(errmsg,2,"it seems that the kernel matrix was not yet read in, however this is required "//&
               "for model regularization of scaling type '"//trim(scltyp)//"'",myname)
          goto 1
       end if
    end if
!
    nullify(SA,Sb)
!
    ! for the same number of rows per section for which the kernel matrix was read, get the equations from the linear smoothing object
    isection = 0
!
    write(errstr,*) "distributing global matrix in sections of rows containing ",this%nrow_handle," rows each"
    call add(errmsg2,0,errstr,myname)
!
    ! loop over subset of rows (always use a section of this%nrow_handle rows)
    do i_rowreg = 1,this%nrowreg,this%nrow_handle ! ig_start is the global index of the first row of the current submatrix
       ig_start = this%ndata + i_rowreg ! global row index of system matrix
       isection = isection +1
       call new(errmsg2,myname)

       nrow_section = min(this%nrow_handle,this%nrowreg-i_rowreg+1) ! number of rows in the current section of rows

       write(errstr,*) "current section of rows: ",isection,"; number of rows in it: ",nrow_section
       call add(errmsg2,0,errstr,myname)

       if(this%myrow == 0 .and. this%mycol == 0) then
          ! get smoothing equations for this section of rows
          write(*,*) "get and send regularization conditions, ",isection,"'th section of rows, ig_start",ig_start,&
               "; i_rowreg ",i_rowreg," ( handling ",this%nrow_handle," rows at once ) "

          if(associated(SA)) deallocate(SA)
          if(associated(Sb)) deallocate(Sb)

          ! get linear model regularization equations for indices i_rowreg:i_rowreg+nrow_section-1
          if(require_minmaxval_K_per_param) then
             call getEquations(lmreg,SA,Sb,errmsg2,ieq_start=i_rowreg,ieq_end=i_rowreg+nrow_section-1,&
                  minmaxval_K_per_param=this%minmaxval_K_per_param)
          else
             call getEquations(lmreg,SA,Sb,errmsg2,ieq_start=i_rowreg,ieq_end=i_rowreg+nrow_section-1)
          end if
          if(.level.errmsg2 == 2) then
             call print(errmsg2)
             call add(errmsg,2,"could not get regularization equations",myname)
             goto 1
          end if
          if(.not.(associated(SA).and.associated(Sb))) then
             call print(errmsg2)
             call add(errmsg,2,"could not get regularization equations, arrays not associated",myname)
             goto 1
          end if
          if(size(SA,2) /= this%nmval) then
             write(errstr,*) "number of columns of returned regularization equations (smoothing/damping) = ",size(SA,2),&
                  "; expected ",this%nmval,"; this suggests that regularization conditions were produced w.r.t. ",&
                  "a different model space than used for initiation of parallell kernel system"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          if(allocated(Srhs)) deallocate(Srhs)
          allocate(Srhs(nrow_section,this%nrhs))
          do irhs = 1,this%nrhs
             Srhs(:,irhs) = Sb
          end do
          deallocate(Sb)
       end if ! this%myrow == this%mycol == 0

       ! distribute section of rows on the blacs process grid
!
       !subroutine distributeGlobalSubmatrixParallelLinearSystem(this,which_array,mat_send,nrow_send,ncol_send,&
       !   ig_start,jg_start,sendprow,sendpcol,errmsg)
       call distributeGlobalSubmatrixParallelLinearSystem(this%PLSE,'matrix',SA,nrow_section,this%nmval,&
            ig_start,1,0,0,errmsg2)
       if(.level.errmsg2 == 2) then
          call print(errmsg2)
          write(errstr,*) isection,"'th section of rows of global regulariztaion matrix could not be distributed ",&
               "on the process grid"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       call distributeGlobalSubmatrixParallelLinearSystem(this%PLSE,'rhs',Srhs,nrow_section,this%nrhs,&
            ig_start,1,0,0,errmsg2)
       if(.level.errmsg2 == 2) then
          call print(errmsg2)
          write(errstr,*) isection,"'th section of rows of global regularization rhs vectors could not be distributed ",&
               "on the process grid"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       call dealloc(errmsg2)
!
    end do ! ig_start
!
!! FS FS DEBUGGING START
!!$if(this%myrow == 0 .and. this%mycol == 0) then
!!$   allocate(mat(this%nrow,this%ncol))
!!$else
!!$   nullify(mat)
!!$end if
!!$call collectGlobalSubmatrixParallelLinearSystem(this%PLSE,'matrix',mat,this%nrow,this%ncol,&
!!$     1,1,0,0,errmsg2)
!!$if(.level.errmsg2 == 2)then
!!$   call print(errmsg2)
!!$   stop
!!$end if
!!$call dealloc(errmsg2)
!!$if(this%myrow == 0 .and. this%mycol == 0) then
!!$   filename = 'DEBUG_matrix_par.dat'
!!$   open(unit=20,file=filename,form='formatted',status='unknown',action='write')
!!$   write(20,*) size(mat,1),size(mat,2)
!!$   if(size(mat,1) > 0 .and. size(mat,2)>0 ) then
!!$      do i = 1,size(mat,1)
!!$         write(20,*) mat(i,:)
!!$      end do
!!$   end if
!!$   close(20)
!!$   deallocate(mat)
!!$   allocate(mat(max(this%nrow,this%ncol),this%nrhs))
!!$end if
!!$call collectGlobalSubmatrixParallelLinearSystem(this%PLSE,'rhs',mat,max(this%nrow,this%ncol),this%nrhs,&
!!$     1,1,0,0,errmsg2)
!!$if(.level.errmsg2 == 2)then
!!$   call print(errmsg2)
!!$   stop
!!$end if
!!$call dealloc(errmsg2)
!!$if(this%myrow == 0 .and. this%mycol == 0) then
!!$   filename = 'DEBUG_rhs_par.dat'
!!$   open(unit=20,file=filename,form='formatted',status='unknown',action='write')
!!$   write(20,*) size(mat,1),size(mat,2)
!!$   if(size(mat,1) > 0 .and. size(mat,2)>0 ) then
!!$      do i = 1,size(mat,1)
!!$         write(20,*) mat(i,:)
!!$      end do
!!$   end if
!!$   close(20)
!!$   deallocate(mat)
!!$end if
!!$call sleep(10)
!!$stop
!! FS FS DEBUGGING END
!
    ! if everything went OK, indicate by flag
    this%row_regularization_set = .true.

    ! CLEAN UP
1   if(associated(SA)) deallocate(SA)
    if(associated(Sb)) deallocate(Sb)
    if(allocated(Srhs)) deallocate(Srhs)
  end subroutine addLinearModelRegularizationParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine setRhsAsDataResidualParKernelLinearSystem(this,&
       nfreq_measured_data,ifreq_measured_data,nfreq_synthetic_data,ifreq_synthetic_data,&
       path_synthetic_data,path_measured_data,lu,misfit,errmsg,&
       apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
       ignore_data_weights,apply_sdata_normalization,apply_mdata_normalization,set_corrected_residual)
    ! incoming
    type (par_kernel_linear_system) :: this
    character(len=*) :: path_synthetic_data,path_measured_data
    integer, dimension(:) ::  ifreq_measured_data,ifreq_synthetic_data
    integer :: nfreq_measured_data,nfreq_synthetic_data,lu
    character(len=*), optional :: path_event_filter,path_station_filter
    logical, optional :: apply_event_filter,apply_station_filter,ignore_data_weights,apply_sdata_normalization,&
         apply_mdata_normalization,set_corrected_residual
    ! returning
    real :: misfit
    type (error_message) :: errmsg
    ! serial kernel linear system
    type (kernel_linear_system) :: KLSE
    real, dimension(:), pointer :: mdata,sdata,scdata
    real, dimension(:,:), allocatable :: rhs
    ! local
    character(len=400)  :: errstr
    character(len=41) :: myname = 'setRhsAsDataResidualParKernelLinearSystem'
    logical :: corrected_residual
!
    nullify(mdata,sdata,scdata)
!
    call addTrace(errmsg,myname)
    ! if there was no call to initiateParKernelLinearSystem before, then always this%im_in_the_grid == .false.
    ! so the following check should be sufficient to test whether this par_kernel_linear_system was initiated (e.g. this%dmspace is set etc.)
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
!
    corrected_residual = .false.
    if(present(set_corrected_residual)) corrected_residual = set_corrected_residual
    if(corrected_residual) then
       write(errstr,*) "Setting right-hand-side as corrected data residuals"
       call add(errmsg,0,errstr,myname)
    end if
!
    if(.not.associated(this%LK)) then
       call add(errmsg,2,"parallel kernel linear system not allocated yet",myname)
       goto 1
    end if
!
    if(this%nrhs /= 1) then
       write(errstr,*) "number of right-hand-sides initiated = ",this%nrhs,"; must be 1 for data residuals as rhs"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    if(this%myrow == 0 .and. this%mycol == 0) then
       call initiateSerialKernelLinearSystem(KLSE,this%dmspace,0,0,errmsg)
       if(.level.errmsg == 2) goto 1
!
       call readMeasuredDataSerialKernelLinearSystem(KLSE,nfreq_measured_data,ifreq_measured_data,&
            path_measured_data,lu,errmsg,ignore_data_weights,apply_mdata_normalization)
       if(.level.errmsg == 2) goto 1
       mdata =>.md.KLSE
       if(.not.associated(mdata)) then
          call add(errmsg,2,"measured data vector returned from serial Kernel system not associated",myname)
          goto 1
       end if
       if(size(mdata) /= this%ndata) then
          write(errstr,*) "there are ",size(mdata)," measured data values returned from serial kernel system, but expected ",&
               this%ndata
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       call readSyntheticDataSerialKernelLinearSystem(KLSE,nfreq_measured_data,ifreq_measured_data,&
            nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
            apply_event_filter=apply_event_filter,path_event_filter=path_event_filter,&
            apply_station_filter=apply_station_filter,path_station_filter=path_station_filter,&
            ignore_data_weights=ignore_data_weights,apply_sdata_normalization=apply_sdata_normalization)
       if(.level.errmsg == 2) goto 1
       sdata =>.sd.KLSE
       if(.not.associated(sdata)) then
          call add(errmsg,2,"synthetic data vector returned from serial Kernel system not associated",myname)
          goto 1
       end if
       if(size(sdata) /= this%ndata) then
          write(errstr,*) "there are ",size(sdata)," synthetic data values returned from serial kernel system, but expected ",&
               this%ndata
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       if(corrected_residual) then
          call readSyntheticDataSerialKernelLinearSystem(KLSE,nfreq_measured_data,ifreq_measured_data,&
               nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
               apply_event_filter=apply_event_filter,path_event_filter=path_event_filter,&
               apply_station_filter=apply_station_filter,path_station_filter=path_station_filter,&
               ignore_data_weights=ignore_data_weights,apply_sdata_normalization=apply_sdata_normalization,&
               read_synthetic_corrections=.true.)
          if(.level.errmsg == 2) goto 1
          scdata =>.scd.KLSE
          if(.not.associated(scdata)) then
             call add(errmsg,2,"synthetic correction data vector returned from serial Kernel system not associated",myname)
             goto 1
          end if
          if(size(scdata) /= this%ndata) then
             write(errstr,*) "there are ",size(scdata),&
                  " synthetic correction data values returned from serial kernel system, but expected ",&
                  this%ndata
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          allocate(rhs(this%ndata,1))       
          rhs(:,1) = mdata - sdata -scdata
       else ! corrected_residual
          allocate(rhs(this%ndata,1))       
          rhs(:,1) = mdata - sdata
       end if ! corrected_residual
!
       misfit = sum(rhs(:,1)**2) ! return the sum of squares of the differences data - synthetics
!
       call dealloc(KLSE)
    end if ! myrow == mycol == 0
!
    ! subroutine distributeGlobalSubmatrixParallelLinearSystem(this,which_array,mat_send,nrow_send,ncol_send,&
    !      ig_start,jg_start,sendprow,sendpcol,errmsg)
    call distributeGlobalSubmatrixParallelLinearSystem(this%PLSE,'rhs',rhs,this%ndata,1,1,1,0,0,errmsg)
!
    ! if everything went OK, indicate by flag
    this%rhs_set = .true.
!
    ! CLEAN UP
1   if(this%myrow == 0 .and. this%mycol == 0) then
       call dealloc(KLSE)
       if(allocated(rhs)) deallocate(rhs)
    end if
  end subroutine setRhsAsDataResidualParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine solveParKernelLinearSystem(this,solution,errmsg)
    type (par_kernel_linear_system) :: this
    ! returning
    real, dimension(:,:), pointer :: solution
    type (error_message) :: errmsg
    ! local
    character(len=400)  :: errstr
    character(len=26) :: myname = 'solveParKernelLinearSystem'
!
    call addTrace(errmsg,myname)
    nullify(solution)
    if(.not.this%im_in_the_grid) return
!
    write(errstr,*) "I am in the process grid at myrow,mycol = ",this%myrow,this%mycol
    call add(errmsg,0,errstr,myname)
!
    if(.not.associated(this%LK)) then
       call add(errmsg,2,"parallel kernel linear system not allocated yet",myname)
       goto 1
    end if
    if(.not.(this%kernel_matrix_set .and. this%rhs_set .and. &
         this%row_regularization_set .and. this%column_regularization_set)) then
       call add(errmsg,2,"either kernel values, or right-hand-sides, or smoothing equations (if initiated) "//&
            "or column regularization (if initiated) was not yet set. Should not solve kernel system in this case.",&
            myname)
       goto 1
    end if
    if(.not.associated(this%Lrhs)) then
       call add(errmsg,2,"right-hand-sides of kernel linear system not yet defined",myname)
       goto 1
    end if
!
    call solveLeastSquaresParallelLinearSystem(this%PLSE,errmsg)
!
    if(this%myrow == 0 .and. this%mycol == 0) then
       allocate(solution(this%nmval,this%nrhs))
    else
       solution => null()
    end if
    call collectGlobalSubmatrixParallelLinearSystem(this%PLSE,'rhs',solution,this%nmval,this%nrhs,1,1,0,0,errmsg)
 !
1   continue
  end subroutine solveParKernelLinearSystem
!------------------------------------------------------------------------
  subroutine deallocateParKernelLinearSystem(this)
    type (par_kernel_linear_system) :: this
    call dealloc(this%PLSE)

    call dealloc(this%dmspace)
    this%nprow = 0
    this%npcol = 0
    this%nbrow = 0
    this%nbcol = 0
    this%myrow = -1
    this%mycol = -1
    this%im_in_the_grid = .false.

    this%nrow = 0
    this%ndata = 0
    this%nrowreg = 0
    this%ncol = 0
    this%nmval = 0
    this%ncolreg = 0
    this%nrhs = 0

    this%kernel_matrix_set = .false.
    this%row_regularization_set = .false.
    this%column_regularization_set = .false.
    this%rhs_set = .false.

    if(associated(this%minmaxval_K_per_param)) deallocate(this%minmaxval_K_per_param)
  end subroutine deallocateParKernelLinearSystem
!------------------------------------------------------------------------
  function getNrowParKernelLinearSystem(this) result(nrow)
    type (par_kernel_linear_system), intent(in) :: this
    integer :: nrow
    nrow = this%nrow
  end function getNrowParKernelLinearSystem
!------------------------------------------------------------------------
  function getNcolParKernelLinearSystem(this) result(ncol)
    type (par_kernel_linear_system), intent(in) :: this
    integer :: ncol
    ncol = this%ncol
  end function getNcolParKernelLinearSystem
  !
end module parKernelLinearSystem
