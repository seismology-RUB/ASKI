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
program solveParKernelSystem
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use linearModelRegularization
  use parKernelLinearSystem
  use modelParametrization
  use kernelInvertedModel
  use inputParameter
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,outfile_base,mpi_parfile,str

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=20) :: myname = 'solveParKernelSystem'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  type (par_kernel_linear_system) :: PKLSE

  type (linear_model_regularization) :: lmreg
  logical :: add_smoothing,add_damping
  character(len=character_length_regscal_type) :: regularization_scaling_type
  real, dimension(:), pointer :: smoothing_scaling_values,damping_scaling_values
  integer :: neq_lmreg,nsmooth_added,ndamp_added
  character(len=100) :: sm_boundary_condition

  ! data normalization
  logical :: normalize_data
  character(len=max_length_string) :: data_normalization_type
  type (kernel_linear_system) :: KLSE_norm
  real, dimension(:), pointer :: mdata_norm,sdata_norm

  real, dimension(:,:), pointer :: solution

  ! paralellization of linear system
  integer :: mypnum,nproc
  integer :: nbrow,nbcol,nprow,npcol,nrows_read_at_once
  logical :: im_in_the_grid
  ! input parameter defining parallelization of linear system
  type (input_parameter) :: mpi_inpar
  character (len=80), dimension(5) :: mpi_inpar_keys
  data mpi_inpar_keys/'NPROC_ROWS', 'NPROC_COLUMNS', &
       'NROW_PER_BLOCK', 'NCOL_PER_BLOCK', 'NROW_TO_PROCESS_AT_ONCE'/

  type (kernel_inverted_model) :: kim_up,kim_ref
  character(len=character_length_param) :: param_name
  character(len=character_length_param), dimension(:), pointer :: pparam
  integer, dimension(:), pointer :: idx_dmspace,pcell

  integer :: ios,lu_out,lu1,lu2
  logical :: outfile_exists,outfile_is_relative
  real :: misfit

!------------------------------------------------------------------------
!  get info from parallel environment
!
  call blacs_pinfo(mypnum,nproc)
  if(nproc < 1) then
     write(*,*) "ERROR: THERE SEEMS TO BE NO PARALLEL ENVIRONMENT, blacs_pinfo returns mypnum,nproc = ",mypnum,nproc
     goto 20
  end if
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Do inversion step by solving the kernel linear system defined by data model space info "//&
       "and regularization constraints. Parallel programm using ScaLAPACK library.")
  ! define positional arguments
  call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
  call addPosarg(ap,'outfile_base','sval','base name of output files (will be used for all files, with suitable extensions')
  call addPosarg(ap,'mpi_parfile','sval','parameter file defining everything related to the ScaLAPACK '//&
       'parallelization of the linear system')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  ! define optional arguments
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
  call addOption(ap,'-odir',.false.,"if set, outfile_base will be assumed relatively to iteration step output files directory. "//&
       "If not set, outfile_base will be assumed relatively to path from where program is run.")
  call addOption(ap,'-normalize',.true.,&
       "type of data normalization; supported types: 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
       "'scale_maxamp_mdata_by_paths'",'sval','')
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     if(mypnum==0) then
        call print(.errmsg.ap)
        call usage(ap)
     end if
     goto 10 !call blacs_abort(0,1)
  end if
  ! get values of positional arguments
  main_parfile = ap.sval.'main_parfile'
  dmspace_file = ap.sval.'dmsi_file'
  outfile_base = ap.sval.'outfile_base'
  mpi_parfile = ap.sval.'mpi_parfile'
  if (.level.(.errmsg.ap) == 2) goto 10

  ! get values of optional arguments, check if all mandatory options are set
  add_smoothing = ap.optset.'-smoothing'
  add_damping = ap.optset.'-damping'

  outfile_is_relative = ap.optset.'-odir'

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
     data_normalization_type = ap.sval.'-normalize'
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency','scale_maxamp_mdata_by_paths')
     case default
        write(*,*) "ERROR: data normalization type '"//trim(data_normalization_type)//&
            "' not supported; supported types are 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
            "'scale_maxamp_mdata_by_paths'"
        goto 10
     end select
  end if
!
   if (.level.(.errmsg.ap) == 2) then
      if(mypnum==0) then
         call print(.errmsg.ap)
         call usage(ap)
      end if
      goto 10
   end if
!
   if((.not.add_smoothing).and.(ap.optset.'-smoothbnd')) then
      if(mypnum==0) then
         write(*,*) "ERROR: -smoothbnd can only be set when adding smoothing equations by option -smoothing. In ",&
              "any case, -smoothbnd is currently ignored!"
         call usage(ap)
      end if
      goto 10
   end if
!
   select case(regularization_scaling_type)
   case('absmax_per_param,overall_factor', 'absmax_per_param,param_factors')
      if(.not.(add_smoothing.or.add_damping)) then
         if(mypnum==0) then
            write(*,*) "ERROR: -regscal can only be requested when adding any smoothing or damping equations by ",&
                 "options -smoothing , -damping."
            call usage(ap)
         end if
         goto 10
      end if
   end select
!
  ! -smoothing , scaling values
  nullify(smoothing_scaling_values)
  if(add_smoothing) then
     smoothing_scaling_values => ap.rvec.'-smoothing'
     if (.level.(.errmsg.ap) == 2) then
        if(mypnum==0) then
           call print(.errmsg.ap)
           call usage(ap)
        end if
        goto 10
     end if
     if(.not.associated(smoothing_scaling_values)) then
        if(mypnum==0) then
           write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
             "even though there was no error parsing argument -smoothing. This is strange..."
           write(*,*) ""
           call usage(ap)
        end if
        goto 10
     end if
  end if
!
  ! -damping , scaling values
  nullify(damping_scaling_values)
  if(add_damping) then
     damping_scaling_values => ap.rvec.'-damping'
     if (.level.(.errmsg.ap) == 2) then
        if(mypnum==0) then
           call print(.errmsg.ap)
           call usage(ap)
        end if
        goto 10
     end if
     if(.not.associated(damping_scaling_values)) then
        if(mypnum==0) then
           write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
                "even though there was no error parsing argument -damping. This is strange..."
           write(*,*) ""
           call usage(ap)
        end if
        goto 10
     end if
  end if
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup mpi parallelization parameter file
  call createKeywordsInputParameter(mpi_inpar,mpi_inpar_keys)
  call new(errmsg,myname)
  call readSubroutineInputParameter(mpi_inpar,get(fuh),mpi_parfile,errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2)  goto 10 !blacs_abort(0,1)
  call dealloc(errmsg)
!
  nprow = ival(mpi_inpar,'NPROC_ROWS',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'NPROC_ROWS' in parallelization parameter file '"//&
          trim(mpi_parfile)//"' is no valid integer"
     goto 10
  end if
  npcol = ival(mpi_inpar,'NPROC_COLUMNS',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'NPROC_COLUMNS' in parallelization parameter file '"//&
          trim(mpi_parfile)//"' is no valid integer"
     goto 10
  end if
  nbrow = ival(mpi_inpar,'NROW_PER_BLOCK',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'NROW_PER_BLOCK' in parallelization parameter file '"//&
          trim(mpi_parfile)//"' is no valid integer"
     goto 10
  end if
  nbcol = ival(mpi_inpar,'NCOL_PER_BLOCK',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'NCOL_PER_BLOCK' in parallelization parameter file '"//&
          trim(mpi_parfile)//"' is no valid integer"
     goto 10
  end if
  nrows_read_at_once = ival(mpi_inpar,'NROW_TO_PROCESS_AT_ONCE',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'NROW_TO_PROCESS_AT_ONCE' in parallelization parameter file '"//&
          trim(mpi_parfile)//"' is no valid integer"
     goto 10
  end if
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2)  goto 10 !blacs_abort(0,1)
  call dealloc(errmsg)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 10 !blacs_abort(0,1)
  call dealloc(errmsg)
!
  if(outfile_is_relative) outfile_base = trim(.iterpath.iterbasics)//&
       trim((.inpar.iterbasics).sval.'PATH_OUTPUT_FILES')//trim(outfile_base)
!
!------------------------------------------------------------------------
!  check if output files already exist,
!  then open output textfile to write
!
  if(mypnum == 0) then
     call document(ap)

     ! check if output text file exists
     inquire(file=trim(outfile_base)//'_out.txt',exist=outfile_exists)
     if(outfile_exists) then
        write(*,*) "ERROR: output text file '"//trim(outfile_base)//"_out.txt' already exists. Please (re)move it."
        goto 20
     end if
     ! check if kernel inverted model update file exists
     inquire(file=trim(outfile_base)//'_up.kim',exist=outfile_exists)
     if(outfile_exists) then
        write(*,*) "ERROR: model uptdate output file '"//trim(outfile_base)//&
             "_up.kim' already exists. Please (re)move it."
        goto 20
     end if
     ! check if kernel inverted model new file exists
     inquire(file=trim(outfile_base)//'_new.kim',exist=outfile_exists)
     if(outfile_exists) then
        write(*,*) "ERROR: new inverted model output file '"//trim(outfile_base)//&
             "_new.kim' already exists. Please (re)move it."
        goto 20
     end if
     ! also check if vtk output files already exist ?! (-> very extensive, as we would have to assume here naming of vtk files as done by kernelInvertedModel module)
!
     lu_out = get(fuh)
     open(unit=lu_out,file=trim(outfile_base)//"_out.txt",status='unknown',form='formatted',action='write',iostat=ios)
     if(ios/=0) then
        write(*,*) "could not open output text file '"//trim(outfile_base)//"_out.txt' to write. Raised iostat = ",ios
        goto 20 !call blacs_abort(0,1)
     end if
     write(lu_out,*) ""; write(lu_out,*) ""
     write(lu_out,*) "welcome to ASKI (put more general information about everything here)"; write(lu_out,*) ""
     write(lu_out,*) "solving the kernel system in parallel, the parallel environment offers ",nproc," processes"
     write(lu_out,*) ""
     write(lu_out,*) "parallelization of the matrix defined in parameter file '"//trim(mpi_parfile)//"': "
     write(lu_out,*) "   number of processes in the rows of the process grid = ",nprow
     write(lu_out,*) "   number of processes in the columns of the process grid = ",npcol
     write(lu_out,*) "   -> total number of required processes for the process grid = ",nprow*npcol
     write(lu_out,*) "   number of rows of a submatrix block used in ScaLAPACK block decomposition = ",nbrow
     write(lu_out,*) "   number of columns of a submatrix block used in ScaLAPACK block decomposition = ",nbcol
     write(lu_out,*) "   number of rows of the kernel system which will be processed in memory of master process = ",&
          nrows_read_at_once
     write(lu_out,*) ""
     write(lu_out,*) "base name of output files will be '"//trim(outfile_base)//"'"
     write(lu_out,*) ""
!
  end if ! mypnum == 0
!------------------------------------------------------------------------
!  setup data model space info object
!
  if(mypnum == 0) then
     write(*,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
     write(lu_out,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
  end if
!
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,sval(.inpar.invbasics,'MODEL_PARAMETRIZATION'),&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 10
  call dealloc(errmsg)
!
  if(mypnum == 0) then
     write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"
     write(lu_out,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"; write(lu_out,*) ""
  end if
!
  ! in case of data normalization, read in the unweighted data and synthetics vectors (if required) 
  if(normalize_data) then
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency')
        if(mypnum == 0) then
           write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vectors of measured data and synthetic data"
           write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vectors of measured data and synthetic data"; write(lu_out,*) ""
        end if
!
        ! in those cases, BOTH mdata and sdata (unweighted) are needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace,0,0,errmsg)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) goto 10
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) goto 10
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
        if (.level.errmsg == 2) goto 10
        call dealloc(errmsg)
        sdata_norm => .sd.KLSE_norm
!
        if(mypnum == 0) then
           write(*,*) "finished reading measured and synthetic data vectors, now computing normalization factors"
        end if
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace,data_normalization_type,errmsg,mdata=mdata_norm,&
             sdata=sdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 10
        call dealloc(errmsg)
!
        nullify(mdata_norm,sdata_norm)
        call dealloc(KLSE_norm)
     case('scale_maxamp_mdata_by_paths')
        if(mypnum == 0) then
           write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vector of measured data"
           write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vector of measured data"; write(lu_out,*) ""
        end if
!
        ! in those cases, only mdata (unweighted) is needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace,0,0,errmsg)
        if (.level.errmsg == 2) goto 10
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 10
        call dealloc(errmsg)
        mdata_norm => .md.KLSE_norm
!
        if(mypnum == 0) then
           write(*,*) "finished reading measured data vector, now computing normalization factors"
        end if
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace,data_normalization_type,errmsg,mdata=mdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 10
        call dealloc(errmsg)
!
        nullify(mdata_norm)
        call dealloc(KLSE_norm)
     end select ! data_normalization_type
  end if ! normalize_data
!------------------------------------------------------------------------
!  create regularization conditions (smoothing/damping)
!
  if(add_smoothing.or.add_damping) then
     select case(regularization_scaling_type)
     case('','none')
        if(mypnum == 0) then
           write(*,*) "initiating linear regularization constraints, no scaling"
           write(lu_out,*) "initiating linear regularization constraints, no scaling"
        end if
     case default
        if(mypnum == 0) then
           write(*,*) "initiating linear regularization constraints, with scaling '"//&
                trim(regularization_scaling_type)//"'"
           write(lu_out,*) "initiating linear regularization constraints, with scaling '"//&
                trim(regularization_scaling_type)//"'"
        end if
     end select
!
     call new(errmsg,myname)
     call init(lmreg,dmspace,errmsg,regularization_scaling_type)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 10
  end if

  if(add_smoothing) then
     if(mypnum == 0) then
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
     end if ! mypnum == 0
!
     ! use the same error message as for initiating lmreg (errmsg was not allocated after init(lmreg,...) )
     call addSmoothing(lmreg,.invgrid.iterbasics,errmsg,scaling_values=smoothing_scaling_values,&
          boundary_conditions=sm_boundary_condition,neq_added=nsmooth_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 10
!
     if(mypnum == 0) then
        write(*,*) "there were ",nsmooth_added," linear smoothing constraints created"
        write(lu_out,*) "there were ",nsmooth_added," linear smoothing constraints created"
     end if ! mypnum == 0
  end if ! add_smoothing
!
  if(add_damping) then
     if(mypnum == 0) then
        if(associated(damping_scaling_values)) then
           write(*,*) "creating damping constraints, scaling values given: ",damping_scaling_values
           write(lu_out,*) "creating damping constraints, scaling values given: ",damping_scaling_values
        else
           write(*,*) "creating damping constraints; there are no scaling values given"
           write(lu_out,*) "creating damping constraints; there are no scaling values given"
        end if
     end if ! mypnum == 0
!
     ! use the same error message as for initiating lmreg and adding smoothing (errmsg was not allocated after init(lmreg,...) )
     call addDamping(lmreg,errmsg,scaling_values=damping_scaling_values,neq_added=ndamp_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 10
!
     if(mypnum == 0) then
        write(*,*) "there were ",ndamp_added," linear damping constraints created"
        write(lu_out,*) "there were ",ndamp_added," linear damping constraints created"
     end if ! mypnum == 0
  end if ! add_damping
!
  if(add_smoothing.or.add_damping) then
     neq_lmreg = .neq.lmreg
     if(mypnum == 0) then
        write(*,*) "in total, there were ",neq_lmreg," linear regularization constraints created"
        write(lu_out,*) "in total, there were ",neq_lmreg," linear regularization constraints created"; write(lu_out,*) ""
     end if ! mypnum == 0
     call dealloc(errmsg)
  else
     neq_lmreg = 0
  end if
!
!------------------------------------------------------------------------
!  read in kernel matrix
!
  if(mypnum == 0) then
     write(*,*) "initiating parallel kernel linear system now"
  end if
  ! subroutine initiateParKernelLinearSystem(this,dmspace,nrowreg,ncolreg,nrhs,&
  !      nbrow,nbcol,nprow,npcol,nrow_handle,im_in_the_grid,errmsg)
  call new(errmsg,myname)
  call initiateParKernelLinearSystem(PKLSE,dmspace,neq_lmreg,0,1,&
       nbrow,nbcol,nprow,npcol,nrows_read_at_once,im_in_the_grid,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if(.level.errmsg == 2) goto 20
  call dealloc(errmsg)

  if(.not.im_in_the_grid) goto 10

  if(mypnum == 0) then
     write(*,*) "allocating parallel kernel linear system now"
  end if
  call new(errmsg,myname)
  !subroutine allocateParKernelLinearSystem(this,errmsg)
  call allocateParKernelLinearSystem(PKLSE,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20 
  call dealloc(errmsg)
!
  if(mypnum == 0) then
     write(*,*) "reading in kernel matrix now and distributing to parallel process grid"
  end if
  ! subroutine readMatrixParKernelLinearSystem(this,df_measured_data,&
  !      nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
  !      lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !      ignore_data_weights,apply_kernel_normalization)
  call new(errmsg,myname)
  lu1 = get(fuh)
  lu2 = get(fuh)
  call readMatrixParKernelLinearSystem(PKLSE,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
       ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
       .ncell.(.invgrid.iterbasics),.pcorr.invbasics,lu1,lu2,errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'))
  call add(fuh,lu1); call add(fuh,lu2)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20
  call dealloc(errmsg)
!!$!! FS FS TEST
!!$write(*,*) "store linear system for debugging reasons in text file 'test-ASKI_solveParKernelSystem_Matrix.dat'"
!!$KM => .KM.KLSE
!!$lu = get(fuh)
!!$open(unit=lu,file='test-ASKI_solveParKernelSystem_Matrix.dat',status='unknown',form='formatted',action='write',iostat=ios)
!!$do lu1 = 1,size(KM,1)
!!$   !write(lu,*) residuals(i),KM(lu1,:)
!!$   write(lu,*) KM(lu1,:)
!!$end do ! i , rows of system
!!$close(lu)
!!$call undo(fuh)
!!$call blacs_abort(0,1)
!!$!! FS FS TEST
!------------------------------------------------------------------------
!  read in measured and synthetic data, compute difference residual and misfit
!
  ! IN CASE OF PATH SPECIFIC MODELS, ACCOUNT FOR SYNTHETIC CORRECTIONS
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
     if(mypnum == 0) then
        write(*,*) "reading in right-hand-side as CORRECTED data residuals and distributing to parallel process grid"
     end if
  else
     if(mypnum == 0) then
        write(*,*) "reading in right-hand-side as regular data residuals and distributing to parallel process grid"
     end if
  end if
!
  ! subroutine setRhsAsDataResidualParKernelLinearSystem(this,&
  !      nfreq_measured_data,ifreq_measured_data,nfreq_synthetic_data,ifreq_synthetic_data,&
  !      path_synthetic_data,path_measured_data,lu,misfit,errmsg,&
  !      apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !      ignore_data_weights,apply_sdata_normalization,apply_mdata_normalization,set_corrected_residual)
  call new(errmsg,myname)
  call setRhsAsDataResidualParKernelLinearSystem(PKLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
       ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),misfit,errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
       set_corrected_residual=lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS'))
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20 
  call dealloc(errmsg)

  if(mypnum == 0) then
     write(*,*) "the misfit (i.e. sum of squares of components of residual vector) is ",misfit
     write(lu_out,*) "the misfit (i.e. sum of squares of components of residual vector) is ",&
          misfit; write(lu_out,*) ""
  end if
!------------------------------------------------------------------------
!  add regularization constraints to kernel linear system
!
  !subroutine addLinearModelSmoothingParKernelLinearSystem(this,lmreg,errmsg)
  call new(errmsg,myname)
  call addLinearModelRegularizationParKernelLinearSystem(PKLSE,lmreg,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20 
  call dealloc(errmsg)
!------------------------------------------------------------------------
!  solve kernel linear system now
!
  if(mypnum == 0) then
     write(*,*) "solve kernel linear system now, having ",.nrow.PKLSE," rows and ",.ncol.PKLSE," columns"
     write(lu_out,*) "solve kernel linear system now, having ",.nrow.PKLSE," rows and ",.ncol.PKLSE," columns"
  end if
  call new(errmsg,myname)
  call solveParKernelLinearSystem(PKLSE,solution,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20 
  call dealloc(errmsg)
!
  if(mypnum == 0) then
     if(.not.associated(solution)) then
        write(*,*) "ERROR: although I am rank 0, I did not get the solution vector by routine '"//&
             "'solveParKernelLinearSystem"
        write(lu_out,*) "ERROR: although I am rank 0, I did not get the solution vector by routine '"//&
             "'solveParKernelLinearSystem"
        goto 20
     end if
     if(size(solution,1) /= .nmval.dmspace .or. size(solution,2) /= 1) then
        write(*,*) "ERROR: returned solution vector has size ",size(solution,1),"-by-",size(solution,2),&
             "; however, expected size is nmval-by-1, where nmval = ",.nmval.dmspace
        write(lu_out,*) "ERROR: returned solution vector has size ",size(solution,1),"-by-",size(solution,2),&
             "; however, expected size is nmval-by-1, where nmval = ",.nmval.dmspace
        goto 20
     end if
!
!! FS FS DEBUG START
!!$     open(unit=20,file='DEBUG_solution_par.dat',form='formatted',status='unknown',action='write')
!!$     write(20,*) size(solution,1),size(solution,2)
!!$     if(size(solution,1) > 0 .and. size(solution,2)>0 ) then
!!$        do lu1 = 1,size(solution,1)
!!$           write(20,*) solution(lu1,:)
!!$        end do
!!$     end if
!!$     close(20)
!!$call sleep(5)
!!$stop
!! FS FS DEBUT END
!------------------------------------------------------------------------
!  create kernel_inverted_model object model update
!
     write(*,*) "creating model objects from solution vector and reference model and computing new model"
!
! create kernel_inverted_model object model update
!
     call new(errmsg,myname)
     call init(kim_up,.pmtrz.dmspace,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     do while(nextParamModelParametrization(.pmtrz.dmspace,param_name))
        if(associated(pparam)) deallocate(pparam)
        if(associated(pcell)) deallocate(pcell); nullify(pcell)
        if(associated(idx_dmspace)) deallocate(idx_dmspace)
        allocate(pparam(1)); pparam(1) = param_name
        idx_dmspace => getIndxModelValues(dmspace,param=pparam,cell=pcell)
        if(associated(idx_dmspace)) then
           write(*,*) "there are ",size(idx_dmspace)," '"//trim(param_name)//"' model update values in the solution vector"
           write(lu_out,*) "there are ",size(idx_dmspace)," '"//trim(param_name)//"' model update values in the solution vector"
           call new(errmsg,myname)
           call addValuesKernelInvertedModel(kim_up,param_name,pcell,solution(idx_dmspace,1),errmsg)
           if (.level.errmsg /= 0) call print(errmsg)
           !call print(errmsg)
           if (.level.errmsg == 2) goto 20
           call dealloc(errmsg)
           deallocate(idx_dmspace)
        else
           write(*,*) "there are no '"//trim(param_name)//"' model update values in the solution vector"
           write(lu_out,*) "there are no '"//trim(param_name)//"' model update values in the solution vector"
        end if
     end do ! while next param
     write(lu_out,*) ""
!
! write update to file(s) now, as model object kim_up will be modified later on when adding kim_ref to kim_up (to get new model)
!
     !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
     call new(errmsg,myname)
     call writeFileKernelInvertedModel(kim_up,trim(outfile_base)//'_up.kim',get(fuh),errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim_up,.invgrid.iterbasics,&
          trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile_base)//'_up',get(fuh),errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
! create kernel_inverted_model object reference model 
!
     !subroutine interpolateKernelReferenceToKernelInvertedModel(this,krm,parametrization,invgrid,intw,errmsg)
     call new(errmsg,myname)
     call interpolateKernelReferenceToKernelInvertedModel(kim_ref,.krm.iterbasics,.pmtrz.dmspace,&
          .invgrid.iterbasics,.intw.iterbasics,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
! create kernel_inverted_model object new model = reference + update by modifying object kim_up
! print new max/min model values
! write new model to files
!
     write(*,*) "updating model now, computing new model values"
     write(lu_out,*) "updating model now, computing new model values"
     !subroutine summateInstancesKernelInvertedModel(this,that,errmsg,c1,c2,relative)
     call new(errmsg,myname)
     call summateInstancesKernelInvertedModel(kim_up,kim_ref,errmsg)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
     do while(nextParamModelParametrization(.pmtrz.kim_up,param_name))
        write(*,*) "max('"//trim(param_name)//"') = ",maxValue(kim_up,param_name),&
             "min('"//trim(param_name)//"') = ",minValue(kim_up,param_name)
        write(lu_out,*) "max('"//trim(param_name)//"') = ",maxValue(kim_up,param_name),&
             "min('"//trim(param_name)//"') = ",minValue(kim_up,param_name)
     end do
     write(lu_out,*) ""
!
     !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
     call new(errmsg,myname)
     call writeFileKernelInvertedModel(kim_up,trim(outfile_base)//'_new.kim',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim_up,.invgrid.iterbasics,&
          trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile_base)//'_new',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     write(lu_out,*) "successfully written all output files"; write(lu_out,*) ""
!
  end if ! mypnum == 0


  ! if program reaches regular end, indicate so in the outputt file
  if(mypnum == 0) then
     write(lu_out,*) "everything went OK, good bye"; write(lu_out,*) ""
  end if


!------------------------------------------------------------------------
!  clean up
!
  ! exit the parallel environment, which is initiated by constructing object PKLSE and exited by destructing it
  ! BY CALLING dealloc(PKLSE), INTERNALLY BLACS_EXIT(0) IS CALLED IN MODULE parallelLinearSystem
  ! THIS CAUSES AN MPI_EXIT AND, HENCE, SHOULD BE CALLED IN ORDER TO SAVELY EXIT THE MPI ENVIRONMENT
10 call dealloc(PKLSE)
!
  if(mypnum == 0) then
     close(lu_out); 
     call add(fuh,lu_out)
  end if
!
  nullify(mdata_norm,sdata_norm)
  call dealloc(KLSE_norm)
!
  call dealloc(dmspace)
!
  call dealloc(invbasics); call dealloc(iterbasics)
  call dealloc(fuh)
  call dealloc(ap)
!
  ! clean stop, everythin went OK
  stop
!
  ! in case an error has occurred, better call blacs_abort (killing all processes), since the error may have occurred only 
  ! on a few but not all processes. Just 'stop' could cause some mpi processes to continue as phantom processes
20 if(mypnum == 0) then
     close(lu_out)
  end if
  call blacs_abort(0,1)
  stop
end program solveParKernelSystem
