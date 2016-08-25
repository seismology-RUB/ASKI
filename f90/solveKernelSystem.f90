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
program solveKernelSystem
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use linearModelRegularization
  use kernelLinearSystem
  use modelParametrization
  use kernelInvertedModel
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  ! basic
  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,outfile,str

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=17) :: myname = 'solveKernelSystem'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  type (kernel_linear_system) :: KLSE

  ! linear model regularization
  type (linear_model_regularization) :: lmreg
  logical :: add_smoothing,add_damping,outfile_is_relative
  character(len=character_length_regscal_type) :: regularization_scaling_type
  integer :: nsmooth_added,ndamp_added
  real, dimension(:), pointer :: smoothing_scaling_values,damping_scaling_values
  character(len=100) :: sm_boundary_condition

  ! data normalization
  logical :: normalize_data
  character(len=100) :: data_normalization_type
  type (kernel_linear_system) :: KLSE_norm
  real, dimension(:), pointer :: mdata_norm,sdata_norm
  
  ! results
  real :: misfit
  real, dimension(:,:), pointer :: delta_mval

  type (kernel_inverted_model) :: kim_new,kim_up_abs,kim_ref
  character(len=character_length_param) :: param_name
  character(len=character_length_param), dimension(:), pointer :: pparam
  integer, dimension(:), pointer :: idx_dmspace,pcell

  ! other
  integer :: ios,lu_out,lu1,lu2
  logical :: outfile_exists


! FS FS TEST
!!$real, dimension(:,:), pointer :: KM
!!$real, dimension(:), pointer :: data
!!$real, dimension(:,:), pointer :: data2
!!$integer :: lu,i
! FS FS TEST

  nullify(smoothing_scaling_values,damping_scaling_values,mdata_norm,sdata_norm,delta_mval,pparam,idx_dmspace,pcell)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Do inversion step by solving the kernel linear system defined by data model space info "//&
       "and regularization constraints. Serial programm using Lapack library.")
  call addPosarg(ap,"dmspace_file","sval","File of data model space info")
  call addPosarg(ap,"outfile_base","sval","Base name of output files (will be used for all files, with "//&
       "suitable extensions)")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-regscal",.true.,"type of scaling of regularization constraints, at the moment only "//&
       "'absmax_per_param,overall_factor', 'absmax_per_param,param_factors' are supported","sval","none")
  call addOption(ap,"-smoothing",.true.,"If set, smoothing conditions are applied. Give a vector of scaling "//&
       "values here consistent with -regscal. 'absmax_per_param,overall_factor': one single factor; "//&
       "'absmax_per_param,param_factors': one factor per parameter name of current parametrization "//&
       "(in conventional order)","rvec","")
  call addOption(ap,"-smoothbnd",.true.,"Define the way (non-existing) neighbours are treated in smoothing "//&
       "conditions at outer/inner boundaries of the inversion grid. Supported types: 'zero_all_outer_bnd',"//&
       "'zero_burried_outer_bnd,cont_free_surface'. If not set, standard average is used everywhere.","sval","")
  call addOption(ap,"-damping",.true.,"If set, damping conditions are applied. Give a vector of scaling "//&
       "values here consistent with -regscal. 'absmax_per_param,overall_factor': one single factor; "//&
       "'absmax_per_param,param_factors': one factor per parameter name of current parametrization "//&
       "(in conventional order); 'none' : values given here are ignored","rvec","")
  call addOption(ap,"-odir",.false.,"If set, outfile_base will be assumed relatively to iteration step output "//&
       "files directory. If not set, outfile_base will be assumed relatively to path from where program is run.")
  call addOption(ap,'-normalize',.true.,&
       "type of data normalization; supported types: 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
       "'scale_maxamp_mdata_by_paths'",'sval','')
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  dmspace_file = ap.sval.'dmspace_file'
  outfile = ap.sval.'outfile_base'
  main_parfile = ap.sval.'main_parfile'
!
  add_smoothing = ap.optset.'-smoothing'
  add_damping = ap.optset.'-damping'
!
  outfile_is_relative = ap.optset.'-odir'
!
  ! -regscal
  str = ap.sval.'-regscal'
  regularization_scaling_type = str
!
  ! -smoothbnd
  if(ap.optset.'-smoothbnd') then
     str = ap.sval.'-smoothbnd'
     sm_boundary_condition = str
  else
     sm_boundary_condition = 'standard'
  end if
!
  normalize_data = ap.optset.'-normalize'
  if(normalize_data) then
     str = ap.sval.'-normalize'
     data_normalization_type = str
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency','scale_maxamp_mdata_by_paths')
     case default
        write(*,*) "ERROR: data normalization type '"//trim(data_normalization_type)//&
            "' not supported; supported types are 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
            "'scale_maxamp_mdata_by_paths'"
        call usage(ap)
        goto 1
     end select
  end if
!
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  if((.not.add_smoothing).and.(ap.optset.'-smoothbnd')) then
     write(*,*) "ERROR: -smoothbnd can only be set when adding smoothing equations by option -smoothing."
     call usage(ap)
     goto 1
  end if
!
  select case(regularization_scaling_type)
  case('absmax_per_param,overall_factor', 'absmax_per_param,param_factors')
     if(.not.(add_smoothing.or.add_damping)) then
        write(*,*) "ERROR: -regscal can only be requested when adding any smoothing or damping equations by ",&
             "options -smoothing , -damping."
        call usage(ap)
        goto 1
     end if
  end select
!
  ! -smoothing , scaling values
  nullify(smoothing_scaling_values)
  if(add_smoothing) then
     smoothing_scaling_values => ap.rvec.'-smoothing'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(smoothing_scaling_values)) then
        write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
             "even though there was no error parsing argument -smoothing. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
  end if
!
  ! -damping , scaling values
  nullify(damping_scaling_values)
  if(add_damping) then
     damping_scaling_values => ap.rvec.'-damping'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(damping_scaling_values)) then
        write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
             "even though there was no error parsing argument -damping. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
  end if
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(main_parfile),get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  if(outfile_is_relative) outfile = trim(.iterpath.iterbasics)//&
       trim((.inpar.iterbasics).sval.'PATH_OUTPUT_FILES')//trim(outfile)
!------------------------------------------------------------------------
!  check if output files already exist
!
  write(*,*) "base name of output files will be '"//trim(outfile)//"'"
!
  ! check if output text file exists
  inquire(file=trim(outfile)//'_out.txt',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "output text file '"//trim(outfile)//"_out.txt' already exists. Please (re)move it."
     goto 1
  end if
  ! check if kernel inverted model update file exists
  inquire(file=trim(outfile)//'_up_abs.kim',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "model uptdate output file '"//trim(outfile)//"_up_abs.kim' already exists. Please (re)move it."
     goto 1
  end if
  ! check if kernel inverted model new file exists
  inquire(file=trim(outfile)//'_new.kim',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "new inverted model output file '"//trim(outfile)//"_new.kim' already exists. Please (re)move it."
     goto 1
  end if
  ! also check if vtk output files already exist ?! (-> very extensive, as we would have to assume here naming of vtk files as done by kernelInvertedModel module)
!------------------------------------------------------------------------
!  open output textfile to write
!
  lu_out = get(fuh)
  open(unit=lu_out,file=trim(outfile)//"_out.txt",status='unknown',form='formatted',action='write',iostat=ios)
  if(ios/=0) then
     write(*,*) "could not open output text file '"//trim(outfile)//"_out.txt' to write. Raised iostat = ",ios
     close(lu_out)
     goto 1
  end if
  write(lu_out,*) ""; write(lu_out,*) ""
  write(lu_out,*) "welcome to ASKI (put more general information about everything here)"; write(lu_out,*) ""
  write(lu_out,*) "inverting data now"; write(lu_out,*) ""
  write(lu_out,*) "base name of output files will be '"//trim(outfile)//"'"; write(lu_out,*) ""
!------------------------------------------------------------------------
!  setup data model space info object
!
  write(*,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
  write(lu_out,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
!
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,sval(.inpar.invbasics,'MODEL_PARAMETRIZATION'),&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  call dealloc(errmsg)
!
  write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model parameters"
  write(lu_out,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model parameters"; write(lu_out,*) ""
!
  ! in case of data normalization, read in the unweighted data and synthetics vectors (if required)
  if(normalize_data) then
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency')
        write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
             "', therefore reading in unweighted complete vectors of measured data and synthetic data"
        write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
             "', therefore reading in unweighted complete vectors of measured data and synthetic data"; write(lu_out,*) ""
!
        ! in those cases, BOTH mdata and sdata (unweighted) are needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace,0,0,errmsg)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        mdata_norm => .md.KLSE_norm
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readSyntheticDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
             ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
             trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
             apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
             path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
             apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
             path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
             ignore_data_weights=.true.,apply_sdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        call dealloc(errmsg)
        sdata_norm => .sd.KLSE_norm
!
        write(*,*) "finished reading measured and synthetic data vectors, now computing normalization factors"
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace,data_normalization_type,errmsg,mdata=mdata_norm,&
             sdata=sdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        call dealloc(errmsg)
!
        nullify(mdata_norm,sdata_norm)
        call dealloc(KLSE_norm)
     case('scale_maxamp_mdata_by_paths')
        write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
             "', therefore reading in unweighted complete vector of measured data"
        write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
             "', therefore reading in unweighted complete vector of measured data"; write(lu_out,*) ""
!
        ! in those cases, only mdata (unweighted) is needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace,0,0,errmsg)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        call dealloc(errmsg)
        mdata_norm => .md.KLSE_norm
!
        write(*,*) "finished reading measured data vector, now computing normalization factors"
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace,data_normalization_type,errmsg,mdata=mdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        call dealloc(errmsg)
!
        nullify(mdata_norm)
        call dealloc(KLSE_norm)
     end select ! data_normalization_type
  end if ! normalize_data
!------------------------------------------------------------------------
!  create smoothing conditions
!
  if(add_smoothing.or.add_damping) then
     select case(regularization_scaling_type)
     case('','none')
        write(*,*) "initiating linear regularization constraints, no scaling"
        write(lu_out,*) "initiating linear regularization constraints, no scaling"
     case default
        write(*,*) "initiating linear regularization constraints, with scaling '"//&
             trim(regularization_scaling_type)//"'"
        write(lu_out,*) "initiating linear regularization constraints, with scaling '"//&
             trim(regularization_scaling_type)//"'"
     end select
!
     call new(errmsg,myname)
     call init(lmreg,dmspace,errmsg,regularization_scaling_type)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  end if

  if(add_smoothing) then
     if(associated(smoothing_scaling_values)) then
        write(*,*) "creating smoothing constraints with boundary conditions '",trim(sm_boundary_condition),&
             "', scaling values given: ",smoothing_scaling_values
        write(lu_out,*) "creating smoothing constraints with boundary conditions '",trim(sm_boundary_condition),&
             "', scaling values given: ",smoothing_scaling_values
     else
        write(*,*) "creating smoothing constraints with boundary conditions '",trim(sm_boundary_condition),&
             "'; there are no scaling values given"
        write(lu_out,*) "creating smoothing constraints with boundary conditions '",trim(sm_boundary_condition),&
             "'; there are no scaling values given"
     end if
!
     ! use the same error message as for initiating lmreg (errmsg was not allocated after init(lmreg,...) )
     call addSmoothing(lmreg,.invgrid.iterbasics,errmsg,scaling_values=smoothing_scaling_values,&
          boundary_conditions=sm_boundary_condition,neq_added=nsmooth_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
!
     write(*,*) "there were ",nsmooth_added," linear smoothing constraints created"
     write(lu_out,*) "there were ",nsmooth_added," linear smoothing constraints created"
  end if ! add_smoothing

  if(add_damping) then
     if(associated(damping_scaling_values)) then
        write(*,*) "creating damping constraints, scaling values given: ",damping_scaling_values
        write(lu_out,*) "creating damping constraints, scaling values given: ",damping_scaling_values
     else
        write(*,*) "creating damping constraints; there are no scaling values given"
        write(lu_out,*) "creating damping constraints; there are no scaling values given"
     end if
!
     ! use the same error message as for initiating lmreg and adding smoothing (errmsg was not allocated after init(lmreg,...) )
     call addDamping(lmreg,errmsg,scaling_values=damping_scaling_values,neq_added=ndamp_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
!
     write(*,*) "there were ",ndamp_added," linear damping constraints created"
     write(lu_out,*) "there were ",ndamp_added," linear damping constraints created"
  end if ! add_damping
!
  if(add_smoothing.or.add_damping) then
     write(*,*) "in total, there were ",.neq.lmreg," linear regularization constraints created"
     write(lu_out,*) "in total, there were ",.neq.lmreg," linear regularization constraints created"; write(lu_out,*) ""
     call dealloc(errmsg)
  end if
!
!------------------------------------------------------------------------
!  read in kernel matrix
!
  write(*,*) "initiating kernel linear system now"
  !   subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
  call new(errmsg,myname)
  call initiateSerialKernelLinearSystem(KLSE,dmspace,.neq.lmreg,0,errmsg)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  call dealloc(errmsg)
!
  write(*,*) "reading in kernel matrix now"
  !subroutine readMatrixSerialKernelLinearSystem(this,df_measured_data,&
  !     nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
  !     lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !     ignore_data_weights,apply_kernel_normalization)
  call new(errmsg,myname)
  lu1 = get(fuh)
  lu2 = get(fuh)
  call readMatrixSerialKernelLinearSystem(KLSE,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
       ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
       .ncell.(.invgrid.iterbasics),.pcorr.invbasics,lu1,lu2,errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
       apply_kernel_normalization=normalize_data)
  call add(fuh,lu1); call add(fuh,lu2)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  call dealloc(errmsg)
!!$!! FS FS TEST
!!$write(*,*) "store linear system for debugging reasons in text file 'test-ASKI_solveKernelSystem_Matrix.dat'"
!!$KM => .KM.KLSE
!!$lu = get(fuh)
!!$open(unit=lu,file='test-ASKI_solveKernelSystem_Matrix.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$do i = 1,size(KM,1)
!!$   !write(lu,*) residuals(i),KM(i,:)
!!$   write(lu,*) KM(i,:)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$stop
!!$!! FS FS TEST
!------------------------------------------------------------------------
!  read in measured and synthetic data, compute difference residual and misfit
!
  write(*,*) "reading in measured data now"
  !subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !     path_measured_data,lu,errmsg,ignore_data_weights,apply_mdata_normalization)
  call new(errmsg,myname)
  call readMeasuredDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,apply_mdata_normalization=normalize_data)
  call undo(fuh)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  call dealloc(errmsg)
!!$!! FS FS TEST
!!$write(*,*) "store measured data for debugging reasons in text file 'test-ASKI_solveKernelSystem_measured.dat'"
!!$data => .md.KLSE
!!$lu = get(fuh)
!!$open(unit=lu,file='test-ASKI_solveKernelSystem_measured.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$do i = 1,size(data)
!!$   write(lu,*) data(i)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$stop
!!$!! FS FS TEST
!
  write(*,*) "reading in synthetic data now"
  call new(errmsg,myname)
  ! subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !     nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
  !     apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !     ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
  call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
       ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
       apply_sdata_normalization=normalize_data,&
       read_synthetic_corrections=.false.)
  call undo(fuh)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
  call dealloc(errmsg)
!!$!! FS FS TEST
!!$write(*,*) "store synthetics for debugging reasons in text file 'test-ASKI_solveKernelSystem_syntetics.dat'"
!!$data => .sd.KLSE
!!$lu = get(fuh)
!!$open(unit=lu,file='test-ASKI_solveKernelSystem_synthetics.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$do i = 1,size(data)
!!$   write(lu,*) data(i)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$stop
!!$!! FS FS TEST
!
! IN CASE OF PATH SPECIFIC MODELS, ACCOUNT FOR SYNTHETIC CORRECTIONS
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then

     write(*,*) "reading in synthetics correction data now"
     ! subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
     !  nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
     !  apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
     !  ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
     call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
          ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
          trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
          apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
          path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
          apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
          path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
          apply_sdata_normalization=normalize_data,&
          read_synthetic_corrections=.true.)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
        call dealloc(errmsg)
!!$!! FS FS TEST
!!$write(*,*) "store synthetics correction for debugging reasons in text file 'test-ASKI_solveKernelSystem_corr.dat'"
!!$data => .scd.KLSE
!!$lu = get(fuh)
!!$open(unit=lu,file='test-ASKI_solveKernelSystem_corr.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$do i = 1,size(data)
!!$   write(lu,*) data(i)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$stop
!!$!! FS FS TEST
!
     write(*,*) "setting right-hand-side to corrected data residuals now"
     !subroutine setRhsAsCorrectedDataResidualKernelLinearSystem(this,errmsg)
     call new(errmsg,myname)
     call setRhsAsCorrectedDataResidualKernelLinearSystem(KLSE,errmsg)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
     call dealloc(errmsg)
!
     misfit = getCorrectedMisfitKernelLinearSystem(KLSE,iostat=ios)
     if(ios/=0) then
        write(*,*) "there was an error computing the corrected misfit, raised iostat = ",ios
        write(lu_out,*) "there was an error computing the misfit, raised iostat = ",ios
        close(lu_out); goto 1
     end if
     write(*,*) "the corrected misfit (i.e. sum of squares of components of corrected residual vector) is ",misfit
     write(lu_out,*) "the corrected misfit (i.e. sum of squares of components of corrected residual vector) is ",&
          misfit; write(lu_out,*) ""

  else ! lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')

     write(*,*) "setting right-hand-side to data residuals now"
     !subroutine setRhsAsDataResidualKernelLinearSystem(this,errmsg)
     call new(errmsg,myname)
     call setRhsAsDataResidualKernelLinearSystem(KLSE,errmsg)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
     call dealloc(errmsg)
!
     misfit = getMisfitKernelLinearSystem(KLSE,iostat=ios)
     if(ios/=0) then
        write(*,*) "there was an error computing the misfit, raised iostat = ",ios
        write(lu_out,*) "there was an error computing the misfit, raised iostat = ",ios
        close(lu_out); goto 1
     end if
     write(*,*) "the misfit (i.e. sum of squares of components of residual vector) is ",misfit
     write(lu_out,*) "the misfit (i.e. sum of squares of components of residual vector) is ",misfit; write(lu_out,*) ""
  end if ! lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')

!------------------------------------------------------------------------
!  add smoothing constraints to kernel linear system
!
   if(add_smoothing.or.add_damping) then
      write(*,*) "adding linear regularization constrains to kernel linear system now"
      !subroutine addToKernelLinearSystemLinearModelSmoothing(this,KLSE,errmsg)
      call new(errmsg,myname)
      call addToKernelLinearSystemLinearModelRegularization(lmreg,KLSE,errmsg)
      !if (.level.errmsg /= 0) call print(errmsg)
      call print(errmsg)
      if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
         call dealloc(errmsg)
   end if
!! FS FS TEST
!!$write(*,*) "store linear system and rhs for debugging reasons in text file 'DEBUG_matrix_serial.dat'"
!!$KM => .KM.KLSE
!!$data2 => .rhs.KLSE
!!$lu = get(fuh)
!!$!!open(unit=lu,file='test-ASKI_solveKernelSystem_system.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$open(unit=lu,file='DEBUG_matrix_serial.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$write(lu,*) size(KM,1),size(KM,2)
!!$do i = 1,size(KM,1)
!!$   write(lu,*) KM(i,:) !data2(i,1),KM(i,:)
!!$end do ! i , rows of system
!!$close(lu)
!!$open(unit=lu,file='DEBUG_rhs_serial.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$write(lu,*) size(data2,1),size(data2,2)
!!$do i = 1,size(data2,1)
!!$   write(lu,*) data2(i,1) !data2(i,1),KM(i,:)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$stop
!! FS FS TEST
!------------------------------------------------------------------------
!  solve kernel linear system now
!
   write(*,*) "solve kernel linear system now, having ",.nrow.KLSE," rows and ",.ncol.KLSE," columns"
   write(lu_out,*) "solve kernel linear system now, having ",.nrow.KLSE," rows and ",.ncol.KLSE," columns"
   !subroutine solveSerialKernelLinearSystem(this,errmsg)
   call new(errmsg,myname)
   call solveSerialKernelLinearSystem(KLSE,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
   delta_mval => .sol.KLSE
   if(associated(delta_mval)) then
      write(*,*) "there are a total of ",size(delta_mval,1)," model update values in the solution vector"
      write(lu_out,*) "there are a total of ",size(delta_mval,1)," model update values in the solution vector"
   else
      write(*,*) "after solving linear system: solution vector not associated! There must have gone sth. wrong"
      write(lu_out,*) "after solving linear system: solution vector not associated! There must have gone sth. wrong"
      close(lu_out)
      goto 1 
   end if
!
   write(lu_out,*) "successfully solved linear system"; write(lu_out,*) ""
!! FS FS TEST
!!$write(*,*) "store solution for debugging reasons in text file 'DEBUG_solution_serial.dat'"
!!$lu = get(fuh)
!!$!!open(unit=lu,file='test-ASKI_solveKernelSystem_system.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$open(unit=lu,file='DEBUG_solution_serial.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$write(lu,*) size(delta_mval,1),size(delta_mval,2)
!!$if(size(delta_mval,1) > 0 .and. size(delta_mval,2)>0 ) then
!!$   do i = 1,size(delta_mval,1)
!!$      write(lu,*) delta_mval(i,:)
!!$   end do
!!$end if
!!$close(lu)
!!$call undo(fuh)
!!$stop
!! FS FS TEST
!------------------------------------------------------------------------
!  retrieve model update from solution vector
!
   write(*,*) "creating model objects from solution vector and reference model and computing new model"
!
! create kernel_inverted_model object model update
!
   call new(errmsg,myname)
   call packVectorToKernelInvertedModel(kim_up_abs,delta_mval(:,1),dmspace,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)

!
! make a copy of kim_up_abs in order to produce both, absolute new model and relative model update
   call copyKernelInvertedModel(kim_new,kim_up_abs)
!
! write absolute update to file(s)
!
   !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim_up_abs,trim(outfile)//'_up_abs.kim',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
   !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim_up_abs,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile)//'_up_abs',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
! create kernel_inverted_model object of reference model 
!
   !subroutine interpolateKernelReferenceToKernelInvertedModel(this,krm,parametrization,invgrid,intw,errmsg)
   call new(errmsg,myname)
   call interpolateKernelReferenceToKernelInvertedModel(kim_ref,.krm.iterbasics,.pmtrz.dmspace,&
        .invgrid.iterbasics,.intw.iterbasics,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
! create kernel_inverted_model object new abs model = abs update + reference by modifying object kim_new, which at this point contains values of kim_up_abs
! print new max/min model values
! write new model to files
!
   write(*,*) "updating model now, computing new absolute model values"
   write(lu_out,*) "updating model now, computing new absolute model values"
   !subroutine summateInstancesKernelInvertedModel(this,that,errmsg,c1,c2,relative)
   call new(errmsg,myname)
   call summateInstancesKernelInvertedModel(kim_new,kim_ref,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
   do while(nextParamModelParametrization(.pmtrz.kim_new,param_name))
         write(*,*) "  new model  max('"//trim(param_name)//"') = ",maxValue(kim_new,param_name),&
              "min('"//trim(param_name)//"') = ",minValue(kim_new,param_name)
         write(lu_out,*) "  new model  max('"//trim(param_name)//"') = ",maxValue(kim_new,param_name),&
              "min('"//trim(param_name)//"') = ",minValue(kim_new,param_name)
   end do
   write(lu_out,*) ""
!
   !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim_new,trim(outfile)//'_new.kim',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
   !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim_new,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile)//'_new',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
! create kernel_inverted_model object relative update [percent] = (0.0 * reference + 100.0 * abs_ update) / reference
! by modifying object kim_ref (which at this point contains values of kim_up_abs)
   write(*,*) "computing relative model update in [percent] now"
   write(lu_out,*) "computing relative model update in [percent] now"
   call new(errmsg,myname)
   call summateInstancesKernelInvertedModel(kim_ref,kim_up_abs,errmsg,0.0,100.0,relative=.true.)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
   do while(nextParamModelParametrization(.pmtrz.kim_ref,param_name))
         write(*,*) "  relative model update [percent]  max('"//trim(param_name)//"') = ",maxValue(kim_ref,param_name),&
              "  min('"//trim(param_name)//"') = ",minValue(kim_ref,param_name)
         write(lu_out,*) "  relative model update [percent]  max('"//trim(param_name)//"') = ",maxValue(kim_ref,param_name),&
              "  min('"//trim(param_name)//"') = ",minValue(kim_ref,param_name)
   end do
   write(lu_out,*) ""
!
   !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim_ref,trim(outfile)//'_up_rel.kim',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
   !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim_ref,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile)//'_up_rel',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); goto 1; endif
   call dealloc(errmsg)
!
   write(lu_out,*) "successfully written all output files"; write(lu_out,*) ""
!
!------------------------------------------------------------------------
!  clean up
!
   write(lu_out,*) "good bye"; write(lu_out,*) ""
   close(lu_out); call add(fuh,lu_out)
!
1  if(associated(smoothing_scaling_values)) deallocate(smoothing_scaling_values)
   if(associated(pparam)) deallocate(pparam)
   if(associated(idx_dmspace)) deallocate(idx_dmspace)
   if(associated(pcell)) deallocate(pcell)
!
   nullify(mdata_norm,sdata_norm)
   call dealloc(KLSE_norm)
!
   call dealloc(lmreg)
   call dealloc(KLSE)
   call dealloc(dmspace)
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
!   
   call dealloc(kim_up_abs); call dealloc(kim_new); call dealloc(kim_ref)
!
end program solveKernelSystem
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'    solveKernelSystem [-h] [-smooth] [-scltyp smoothing_scaling_type] '
!   print *,'        [-sclval smoothing_scaling_values] [-smbnd] [-odir] dmspace_file outfile_base parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    dmspace_file: data model space input file which defines data and model space"
!   print *,"    outfile_base: base name of output files (will be used for all files, with suitable extensions)"
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print *,'-smooth  : indicates if linear smoothing constraints (average neighbour) are added to the kernel linear system'
!   print *,''
!   print *,"-scltyp  : type of scaling of smoothing constraints, at the moment only 'absmax_per_param,overall_factor' "
!   print *,"           is supported: requires TWO value of -sclval:"
!   print *,"           - the first value gives the smoothing intensity factor (set 0.0 if no smoothing should be done)"
!   print *,"           - the second value gives the damping intensity (set to 0.0 if no damping should be applied)"
!   print *,''
!   print *,'-sclval  : dependent on -scltyp, smoothing_scaling_values is of form "n val1 .. valn" giving the number of '//&
!        'values n followed by n real numbers'
!   print *,''
!   print *,'-smbnd  : is ignored at the moment, will in the future define the way (non-existing) neighbours are treated '//&
!        'at the boundaries of the inversion grid'
!   print *,''
!   print *,'-odir  : if set, outfile_base will be assumed relatively to iteration step output files directory'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
