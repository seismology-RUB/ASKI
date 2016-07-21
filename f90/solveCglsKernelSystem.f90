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

!  Program to perform the solution the least-squares
!  problem |(d|s)-(K|S)x|^2 -> min for kernel matrix 
!  K and regularization matrix S using conjugate gradients. 
!  The algorithm applied here is the "CGLS1" from paper:
!  
!     Ake Bjoerck, Tommy Elfving and Zdenek Strakos
!     "Stability of conjugate gradient and Lanczos methods
!     for linear least squares problems"
!     SIAM J. Matrix Anal. APPL., 19, 720-736, 1998.'
!
!  By request (controlled by parameter niter_recompute_r), 
!  this algorithm is modified in the way that circularly
!  after a certain number of iterations always the true residual vector 
!  is recomputed instead of always updating it.
!  SEEMS NOT TO WORK (like that), THE ALGORITHM SLOWS DOWN
!  SIGNIFICANTLY! Hence, this feature is switched off in the code
!  (in a hard-coded way, ignoring the value of niter_recompute_r
!  which is given in the parameter file).

program solveCglsKernelSystem

  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use linearModelRegularization
  use kernelLinearSystem
  use modelParametrization
  use kernelInvertedModel
  use vectorPointer
  use inputParameter
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage
  use mpiSupport

  implicit none

  integer, parameter :: SIZE_REAL = 4
  integer, parameter :: SIZE_DOUBLE = 8
  integer :: CUSTOM_MPI_TYPE

!#########################################################################################
!####   SWITCH HERE BETWEEN SINGLE AND DOUBLE PRECISION FOR CG ALGORITHM
!
! in order to use single precision for the CG algorithm variables, set to SIZE_REAL
! in order to use single precision for the CG algorithm variables, set to SIZE_DOUBLE
  integer, parameter :: CUSTOM_REAL = SIZE_REAL
!
!####   DO NOT MODIFY ANYTHING ELSE (unless you know what you're doing)
!#########################################################################################


  ! GLOBAL VARIABLES (on all procs the same content)

  real, parameter :: HUGEVAL = huge(1.0)

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,outfile_base,CG_parfile,startsol_file,str
  logical :: use_nontrivial_starting_solution

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace_glob
  character(len=character_length_pmtrz) :: parametrization
  integer :: npath_glob,ndata_glob,ncol_glob,nparam_pmtrz

  type (linear_model_regularization) :: lmreg
  logical :: add_smoothing,add_damping,add_regularization
  character(len=character_length_regscal_type) :: regularization_scaling_type
  real, dimension(:), pointer :: smoothing_scaling_values,damping_scaling_values
  integer :: neq_lmreg_glob,nsmooth_added,ndamp_added
  character(len=100) :: sm_boundary_condition

  real, dimension(:,:), allocatable :: minmaxval_K_glob_per_param


  ! LOCAL VARIABLES (locally different values, dependent on proc)

  ! data normalization
  logical :: normalize_data
  character(len=max_length_string) :: data_normalization_type
  type (kernel_linear_system) :: KLSE_norm
  real, dimension(:), pointer :: mdata_norm,sdata_norm

  type (kernel_linear_system) :: KLSE_loc
  type (data_model_space_info) :: dmspace_loc
  real, dimension(:,:), pointer :: K_loc
  double precision, dimension(:,:), allocatable :: K_loc_dble
  real, dimension(:,:), pointer :: rhs_loc
  double precision, dimension(:,:), allocatable :: rhs_loc_dble
  real, dimension(:,:), allocatable :: minmaxval_K_loc_per_param
  integer :: ndata_loc,npath_loc,ipath1_loc,ipath2_loc

  integer :: nregul_loc,iregul1_loc,iregul2_loc,ieq,imval
  type (real_vector_pointer), dimension(:), pointer :: regeq_coef,regeq_T_coef
  type (integer_vector_pointer), dimension(:), pointer :: regeq_indx,regeq_T_indx
  real, dimension(:), pointer :: regeq_rhs,coef,coef_T
  double precision, dimension(:), allocatable :: regeq_rhs_dble
  integer, dimension(:), pointer :: indx,indx_T

  logical :: there_is_local_data,there_is_only_local_data,&
       there_is_local_regularization,there_is_only_local_regularization,&
       there_is_local_data_and_regularization


  ! cgls1 system matrix and rhs vector
  ! cgls1 vectors for kernel/data part of the system
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: qd,rd,sd_loc
  ! cgls1 vectors for regularization part of the system
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: qs,rs
  ! other cgls1 variables
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: p,s_loc,s,solution
  real(kind=CUSTOM_REAL) :: alpha,gamma,gamma_old
  real(kind=CUSTOM_REAL) :: norm_qd_squared,norm_qs_squared,norm_q_squared,&
       norm_rd_squared_loc,norm_rs_squared_loc,norm_r,&
       norm_rd_squared,norm_rs_squared,&
       norm_solution,atol,atol_norm_A,btol,btol_norm_b
  integer :: max_niter,i_iter,niter_recompute_r
  logical :: terminate,recompute_r_at_all,recompute_r
  double precision, dimension(:), allocatable :: norm_r_array_sta,norm_r_array_lta
  integer :: nsta,nlta,max_nlta_nsta,ista,ilta
  double precision :: sta,lta,lta_old!,err_sta,err_lta !! uncomment, when using err_sta,err_lta below in convergence_criterion_is_met
  real :: eps_single_precision!,r_tmp uncomment, when using r_tmp below in convergence_criterion_is_met
  real(kind=CUSTOM_REAL) :: cr_tmp


  ! OTHER STUFF

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=21) :: myname = 'solveCglsKernelSystem'

  ! parallellization of linear system
  type (mpi_support) :: mpisup
  integer, dimension(:), allocatable :: ipath1_loc_procs,ipath2_loc_procs,npath_loc_procs,ndata_loc_procs,&
       iregul1_loc_procs,iregul2_loc_procs,nregul_loc_procs
  ! input parameter defining details of CG algorithm
  type (input_parameter) :: CG_inpar
  character (len=80), dimension(6) :: CG_inpar_keys
  data CG_inpar_keys/'MAX_NUM_CG_ITERATIONS', 'NITER_WINDOW_STA', 'NITER_WINDOW_LTA', &
       'USE_BLAS_LEVEL_1', 'USE_BLAS_LEVEL_2', 'NITER_RECOMPUTE_RESIDUAL'/

  type (kernel_inverted_model) :: kim_up_abs,kim_new,kim_ref
  character(len=character_length_param) :: param_name
  character(len=character_length_param), dimension(:), pointer :: pparam
  integer, dimension(:), pointer :: idx_dmspace,pcell,idx_kim
  real, dimension(:), pointer :: model_vector
  integer :: iparam_pmtrz

  integer :: n_remain,n,j
  integer :: ios,lu_out,lu_plot,lu1,lu2
  logical :: file_exists,outfile_is_relative

  logical :: use_dble_cgls,use_blas1,use_blas2

  ! EXTERNAL BLAS FUNCTIONS

  real :: SDOT
  external SDOT

  double precision :: DDOT
  external DDOT

  real :: SNRM2
  external SNRM2

  double precision :: DNRM2
  external DNRM2

  nullify(smoothing_scaling_values,damping_scaling_values,mdata_norm,sdata_norm,K_loc,rhs_loc,&
       regeq_coef,regeq_T_coef,regeq_indx,regeq_T_indx,regeq_rhs,coef,coef_T,indx,indx_T,&
       pparam,idx_dmspace,pcell,idx_kim,model_vector)

!------------------------------------------------------------------------
!  single or double precision for CG algorithm?
!
  if(CUSTOM_REAL == SIZE_REAL) then
     use_dble_cgls = .false.
  elseif(CUSTOM_REAL == SIZE_DOUBLE) then
     use_dble_cgls = .true.
  else
     write(*,*) "ERROR: CUSTOM_REAL = ",CUSTOM_REAL,"; CUSTOM_REAL MUST EQUAL EITHER SIZE_REAL = ",SIZE_REAL,&
          " OR SIZE_DOUBLE = ",SIZE_DOUBLE,"; NOTHING ELSE SUPPORTED!"
     goto 20
  end if
!
!------------------------------------------------------------------------
!  initiate MPI environment
!
  call new(mpisup)
  if(.numtasks.mpisup <= 0) then
     write(*,*) "ERROR: THERE SEEMS TO BE NO PARALLEL ENVIRONMENT, as number of mpi ranks equals ",.numtasks.mpisup
     goto 20
  end if
!
  ! set mpi type of float data which is communicated via mpi
  if(use_dble_cgls) then
     CUSTOM_MPI_TYPE = MPI_DOUBLE_PRECISION
  else
     CUSTOM_MPI_TYPE = MPI_REAL
  end if
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Do inversion step by solving the kernel linear system defined by data model space info "//&
       "and regularization constraints in least squares optimization by parallelized conjugate-gradient method")
  ! define positional arguments
  call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
  call addPosarg(ap,'outfile_base','sval','base name of output files (will be used for all files, with suitable extensions)')
  call addPosarg(ap,'CG_parfile','sval','parameter file defining details related to the conjugate gradient algorithm '//&
       'used to solve the linear system')
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
  call addOption(ap,"-odir",.false.,"If set, outfile_base will be assumed relatively to iteration step output "//&
       "files directory. If not set, outfile_base will be assumed relatively to path from where program is run.")
  call addOption(ap,'-startsol',.true.,&
       "defines starting solution of linear system, which is optionally to be used (kim file containing model update values",&
       'sval','up.kim')
  call addOption(ap,'-normalize',.true.,&
       "type of data normalization; supported types: 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
       "'scale_maxamp_mdata_by_paths'",'sval','')

  call parse(ap)
  if (.level.(.errmsg.ap) == 2) goto 30

  ! get values of positional arguments
  main_parfile = ap.sval.'main_parfile'
  dmspace_file = ap.sval.'dmsi_file'
  CG_parfile = ap.sval.'CG_parfile'
  outfile_base = ap.sval.'outfile_base'
  if (.level.(.errmsg.ap) == 2) goto 30

  add_smoothing = ap.optset.'-smoothing'
  add_damping = ap.optset.'-damping'
  add_regularization = add_smoothing.or.add_damping

  outfile_is_relative = ap.optset.'-odir'

  ! -regscal
  str = ap.sval.'-regscal'
  regularization_scaling_type = str

  ! -smoothbnd
  if(ap.optset.'-smoothbnd') then
     str = ap.sval.'-smoothbnd'
     sm_boundary_condition = str
  else
     sm_boundary_condition = 'standard'
  end if

  normalize_data = ap.optset.'-normalize'
  if(normalize_data) then
     data_normalization_type = ap.sval.'-normalize'
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency','scale_maxamp_mdata_by_paths')
     case default
        write(*,*) "ERROR: data normalization type '"//trim(data_normalization_type)//&
            "' not supported; supported types are 'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', "//&
            "'scale_maxamp_mdata_by_paths'"
        goto 30
     end select
  end if
!
  ! read in starting solution (if requested) from kim file, interpret it using current inversion grid and data model space info file.
  use_nontrivial_starting_solution = ap.optset.'-startsol'
  if(use_nontrivial_starting_solution) startsol_file = ap.sval.'-startsol'
!
  if (.level.(.errmsg.ap) == 2) goto 30
!
   if((.not.add_smoothing).and.(ap.optset.'-smoothbnd')) then
      write(*,*) "ERROR: -smoothbnd can only be set when adding smoothing equations by option -smoothing."
      goto 30
   end if
!
   select case(regularization_scaling_type)
   case('absmax_per_param,overall_factor', 'absmax_per_param,param_factors')
      if(.not.(add_regularization)) then
         write(*,*) "ERROR: -regscal can only be requested when adding any smoothing or damping equations by ",&
              "options -smoothing , -damping."
         goto 30
      end if
   end select
!
  ! -smoothing , scaling values
  nullify(smoothing_scaling_values)
  if(add_smoothing) then
     smoothing_scaling_values => ap.rvec.'-smoothing'
     if (.level.(.errmsg.ap) == 2) goto 30
     if(.not.associated(smoothing_scaling_values)) then
        write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
             "even though there was no error parsing argument -smoothing. This is strange..."
        write(*,*) ""
        goto 30
     end if
  end if
!
  ! -damping , scaling values
  nullify(damping_scaling_values)
  if(add_damping) then
     damping_scaling_values => ap.rvec.'-damping'
     if (.level.(.errmsg.ap) == 2) goto 30
     if(.not.associated(damping_scaling_values)) then
        write(*,*) "ERROR: for some reason, there is no vector of scaling values returned by argument parser, "//&
             "even though there was no error parsing argument -damping. This is strange..."
        write(*,*) ""
        goto 30
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
  call createKeywordsInputParameter(CG_inpar,CG_inpar_keys)
  call new(errmsg,myname)
  call readSubroutineInputParameter(CG_inpar,get(fuh),CG_parfile,errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2)  goto 20
  call dealloc(errmsg)

  ! 'MAX_NUM_CG_ITERATIONS'
  max_niter = ival(CG_inpar,'MAX_NUM_CG_ITERATIONS',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value for key 'MAX_NUM_CG_ITERATIONS' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid 4-byte signed integer (maximum value is 2.147.483.647)"
     goto 20
  end if
  if(max_niter < 1) then
     write(*,*) "ERROR: value ",max_niter," for key 'MAX_NUM_CG_ITERATIONS' in CG parameter file '"//&
          trim(CG_parfile)//"' is < 1"
     goto 20
  end if

  ! 'NITER_WINDOW_STA'
  nsta = ival(CG_inpar,'NITER_WINDOW_STA',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value '"//trim(CG_inpar.sval.'NITER_WINDOW_STA')//&
          "' for key 'NITER_WINDOW_STA' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid integer"
     goto 20
  end if
  if(nsta < 1) then
     write(*,*) "ERROR: value ",nsta," for key 'NITER_WINDOW_STA' in CG parameter file '"//&
          trim(CG_parfile)//"' is < 1"
     goto 20
  end if
  ! 'NITER_WINDOW_LTA'
  nlta = ival(CG_inpar,'NITER_WINDOW_LTA',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value '"//trim(CG_inpar.sval.'NITER_WINDOW_LTA')//&
          "' for key 'NITER_WINDOW_LTA' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid integer"
     goto 20
  end if
  if(nlta < 1) then
     write(*,*) "ERROR: value ",nlta," for key 'NITER_WINDOW_LTA' in CG parameter file '"//&
          trim(CG_parfile)//"' is < 1"
     goto 20
  end if
  if(nsta>nlta) then
     write(*,*) "ERROR: value ",nsta," for key 'NITER_WINDOW_STA' in CG parameter file '"//&
          trim(CG_parfile)//"' must not be larger than value ",nlta," for key 'NITER_WINDOW_LTA'"
     goto 20
  end if

  ! 'USE_BLAS_LEVEL_1'
  use_blas1 = lval(CG_inpar,'USE_BLAS_LEVEL_1',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value '"//trim(CG_inpar.sval.'USE_BLAS_LEVEL_1')//&
          "' for key 'USE_BLAS_LEVEL_1' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid logical"
     goto 20
  end if
  ! 'USE_BLAS_LEVEL_2'
  use_blas2 = lval(CG_inpar,'USE_BLAS_LEVEL_2',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: value '"//trim(CG_inpar.sval.'USE_BLAS_LEVEL_2')//&
          "' for key 'USE_BLAS_LEVEL_2' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid logical"
     goto 20
  end if

  ! 'NITER_RECOMPUTE_RESIDUAL'
  !niter_recompute_r = ival(CG_inpar,'NITER_RECOMPUTE_RESIDUAL',iostat=ios) ! HARDCODED SWITCHED OFF: uncomment this line and delete (comment) the following, if you want to use this feature
  niter_recompute_r = 0; ios = 0 ! this featuer is hard-coded switched off!, see documentation in solveCglsKernelSystem_parfile_template
  if(ios/=0) then
     write(*,*) "ERROR: value '"//trim(CG_inpar.sval.'NITER_RECOMPUTE_RESIDUAL')//&
          "' for key 'NITER_RECOMPUTE_RESIDUAL' in CG parameter file '"//&
          trim(CG_parfile)//"' is no valid integer"
     goto 20
  end if
  recompute_r_at_all = niter_recompute_r > 0

  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20
  call dealloc(errmsg)
!
  if(outfile_is_relative) outfile_base = trim(.iterpath.iterbasics)//&
       trim((.inpar.iterbasics).sval.'PATH_OUTPUT_FILES')//trim(outfile_base)
!
!------------------------------------------------------------------------
!  check if output files already exist,
!  then open output textfile to write
!
  if(.myrank.mpisup == 0) then
     call document(ap)
     write(*,*) ""

     ! check if output text file exists
     inquire(file=trim(outfile_base)//'_out.txt',exist=file_exists)
     if(file_exists) then
        write(*,*) "ERROR: output text file '"//trim(outfile_base)//"_out.txt' already exists. Please (re)move it."
        goto 20
     end if
     ! check if plotting text file exists
     inquire(file=trim(outfile_base)//'_plot.txt',exist=file_exists)
     if(file_exists) then
        write(*,*) "ERROR: plotting text file '"//trim(outfile_base)//"_plot.txt' already exists. Please (re)move it."
        goto 20
     end if
     ! check if kernel inverted model update file exists
     inquire(file=trim(outfile_base)//'_up_abs.kim',exist=file_exists)
     if(file_exists) then
        write(*,*) "ERROR: model uptdate output file '"//trim(outfile_base)//&
             "_up_abs.kim' already exists. Please (re)move it."
        goto 20
     end if
     ! check if kernel inverted model new file exists
     inquire(file=trim(outfile_base)//'_new.kim',exist=file_exists)
     if(file_exists) then
        write(*,*) "ERROR: new inverted model output file '"//trim(outfile_base)//&
             "_new.kim' already exists. Please (re)move it."
        goto 20
     end if
     ! also check if vtk output files already exist ?! (-> very extensive, as we would have to assume here naming of vtk files as done by kernelInvertedModel module)
!
     lu_plot = get(fuh)
     open(unit=lu_plot,file=trim(outfile_base)//"_plot.txt",status='unknown',form='formatted',action='write',iostat=ios)
     if(ios/=0) then
        write(*,*) "could not open plotting text file '"//trim(outfile_base)//"_plot.txt' to write. Raised iostat = ",ios
        goto 20 !call blacs_abort(0,1)
     end if
     lu_out = get(fuh)
     open(unit=lu_out,file=trim(outfile_base)//"_out.txt",status='unknown',form='formatted',action='write',iostat=ios)
     if(ios/=0) then
        write(*,*) "could not open output text file '"//trim(outfile_base)//"_out.txt' to write. Raised iostat = ",ios
        goto 20 !call blacs_abort(0,1)
     end if
     write(lu_out,*) ""; write(lu_out,*) ""
     write(lu_out,*) "welcome to ASKI (put more general information about everything here)"; write(lu_out,*) ""
     write(lu_out,*) "solving the kernel system in parallel by conjugate gradient method, the parallel environment offers ",&
          .numtasks.mpisup," processes"
     write(lu_out,*) ""
     write(lu_out,*) "CG parameter file '"//trim(CG_parfile)//"' tells: "
     write(lu_out,*) "   upper limit of number of iterations will be ",max_niter
     write(lu_out,*) "   number of iterations over which the short-term average is build will be ",nsta
     write(lu_out,*) "   number of iterations over which the long-term average is build will be ",nlta
     if(use_dble_cgls) then
        write(lu_out,*) "   will use double precision routines for CG algorithm"
     else
        write(lu_out,*) "   will use single precision routines for CG algorithm: IF THERE ARE ERRORS DUE TO NaN ",&
             "VALUES, PLEASE CONSIDER TO TRY DOUBLE PRECISION (change the appropriate line at the ",&
             "beginning of code file solveCglsKernelSystem.f90 and recompile program solveCglsKernelSystem "
        write(*,*) "ATTENTION: will use single precision routines for CG algorithm: IF THERE ARE ERRORS DUE TO NaN ",&
             "VALUES, PLEASE CONSIDER TO TRY DOUBLE PRECISION (change the appropriate line at the ",&
             "beginning of code file solveCglsKernelSystem.f90 and recompile program solveCglsKernelSystem "
     end if
     if(use_blas1) then
        write(lu_out,*) "   WILL use BLAS level 1 routines"
     else
        write(lu_out,*) "   will NOT use BLAS level 1 routines"
     end if
     if(use_blas2) then
        write(lu_out,*) "   WILL use BLAS level 2 routines"
     else
        write(lu_out,*) "   will NOT use BLAS level 2 routines"
     end if
     if(recompute_r_at_all) then
        write(lu_out,*) "   residual vector will be recomputed every ",niter_recompute_r," iterations"
     else
        write(lu_out,*) "   residual vector will NOT be recomputed, i.e. using approximations only as in original CG algorithm"
     end if
     write(lu_out,*) ""
     if(use_nontrivial_starting_solution) then
        write(lu_out,*) "will use starting solution for conjugate gradient method taken from 'kim' file '"//&
             trim(startsol_file)//"'"
     else
        write(lu_out,*) "will use trivial solution x = 0 for conjugate gradient method"
     end if
     write(lu_out,*) ""
     write(lu_out,*) "base name of output files will be '"//trim(outfile_base)//"'"
     write(lu_out,*) ""
     write(lu_out,*) "for plotting, file '"//trim(outfile_base)//"_plot.txt' will be written"
     write(lu_out,*) "plot e.g. by gnuplot using the following commands:"
     write(lu_out,*) "   plot residual norm and sta, lta of it:"
     write(lu_out,*) '      plot ',&
          '"',trim(outfile_base)//'_plot.txt" using 1:2 w l title column, ',&
          '"',trim(outfile_base)//'_plot.txt" using 1:3 w l title column, ',&
          '"',trim(outfile_base)//'_plot.txt" using 1:4 w l title column'
     write(lu_out,*) "   plot norm of solution vector x:"
     write(lu_out,*) '      plot ',&
          '"',trim(outfile_base)//'_plot.txt" using 1:5 w l title column'
     write(lu_out,*) ""
!
  end if ! .myrank.mpisup == 0
  call barrier(mpisup)
!------------------------------------------------------------------------
!  setup data model space info object
!
  if(.myrank.mpisup == 0) then
     write(*,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
     write(lu_out,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
  end if
!
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace_glob,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,sval(.inpar.invbasics,'MODEL_PARAMETRIZATION'),&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 20
  call dealloc(errmsg)
!
  ndata_glob = .ndata.dmspace_glob
  npath_glob = .npath.dmspace_glob
  ncol_glob = .nmval.dmspace_glob
  parametrization = .pmtrz.dmspace_glob
  nparam_pmtrz = numberOfParamModelParametrization(parametrization)
!
  if(.myrank.mpisup == 0) then
     write(*,*) "there are ",.ndata.dmspace_glob," data samples and ",.nmval.dmspace_glob," model values"
     write(lu_out,*) "there are ",.ndata.dmspace_glob," data samples and ",.nmval.dmspace_glob," model values"; write(lu_out,*) ""
  end if
!
  ! in case of data normalization, read in the unweighted data and synthetics vectors (if required) 
  if(normalize_data) then
     select case(data_normalization_type)
     case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency')
        if(.myrank.mpisup == 0) then
           write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vectors of measured data and synthetic data"
           write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vectors of measured data and synthetic data"; write(lu_out,*) ""
        end if
!
        ! in those cases, BOTH mdata and sdata (unweighted) are needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace_glob,0,0,errmsg)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) goto 20
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg == 2) call print(errmsg)
        if (.level.errmsg == 2) goto 20
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
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
        sdata_norm => .sd.KLSE_norm
!
        if(.myrank.mpisup == 0) then
           write(*,*) "finished reading measured and synthetic data vectors, now computing normalization factors"
        end if
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace_glob,data_normalization_type,errmsg,mdata=mdata_norm,&
             sdata=sdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
!
        nullify(mdata_norm,sdata_norm)
        call dealloc(KLSE_norm)
     case('scale_maxamp_mdata_by_paths')
        if(.myrank.mpisup == 0) then
           write(*,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vector of measured data"
           write(lu_out,*) "computing data normalization of type '",trim(data_normalization_type),&
                "', therefore reading in unweighted complete vector of measured data"; write(lu_out,*) ""
        end if
!
        ! in those cases, only mdata (unweighted) is needed to compute the normalization
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE_norm,dmspace_glob,0,0,errmsg)
        if (.level.errmsg == 2) goto 20
!
        ! keep using the above error message, for better debugging in case there is an error reading in data
        call readMeasuredDataSerialKernelLinearSystem(KLSE_norm,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,ignore_data_weights=.true.,&
             apply_mdata_normalization=.false.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
        mdata_norm => .md.KLSE_norm
!
        if(.myrank.mpisup == 0) then
           write(*,*) "finished reading measured data vector, now computing normalization factors"
        end if
        call new(errmsg,myname)
        !subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
        call createDataNormalizationDataModelSpaceInfo(dmspace_glob,data_normalization_type,errmsg,mdata=mdata_norm)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
!
        nullify(mdata_norm)
        call dealloc(KLSE_norm)
     end select ! data_normalization_type
  end if ! normalize_data
!------------------------------------------------------------------------
!  read in non-trivial starting solution, if requested, otherwise initiate trivial solution
!
  allocate(solution(ncol_glob))
  if(use_nontrivial_starting_solution) then
     call new(errmsg,myname)
     call readFileKernelInvertedModel(kim_up_abs,startsol_file,get(fuh),errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     call new(errmsg,myname)
     call unpackToVectorKernelInvertedModel(model_vector,kim_up_abs,dmspace_glob,errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     if(.not.associated(model_vector)) then
        write(*,*) "ERROR: model vector returned from routine 'unpackToVectorKernelInvertedModel' is empty"
        goto 20
     end if
     if(size(model_vector) /= size(solution)) then
        write(*,*) "ERROR: model vector returned from routine 'unpackToVectorKernelInvertedModel' has size ",size(model_vector),&
             ", but requested values (number of model parameters in model space) is ",size(solution)
        goto 20
     end if
     solution = model_vector ! here possible conversion from single to double precision if use_dble_cgls
     deallocate(model_vector)
!
  else ! use_nontrivial_starting_solution
!
     solution = 0._CUSTOM_REAL
!
  end if ! use_nontrivial_starting_solution
!------------------------------------------------------------------------
!  create regularization conditions
!
  if(add_regularization) then
     if(.myrank.mpisup == 0) then
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
     end if ! myrank == 0
!
     call new(errmsg,myname)
     call init(lmreg,dmspace_glob,errmsg,regularization_scaling_type)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
  end if

  if(add_smoothing) then
     if(.myrank.mpisup == 0) then
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
     end if ! myrank == 0
!
     ! use the same error message as for initiating lmreg (errmsg was not allocated after init(lmreg,...) )
     call addSmoothing(lmreg,.invgrid.iterbasics,errmsg,scaling_values=smoothing_scaling_values,&
          boundary_conditions=sm_boundary_condition,neq_added=nsmooth_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
!
     if(.myrank.mpisup == 0) then
        write(*,*) "there were ",nsmooth_added," linear smoothing constraints created"
        write(lu_out,*) "there were ",nsmooth_added," linear smoothing constraints created"
     end if ! myrank == 0
  end if ! add_smoothing

  if(add_damping) then
     if(.myrank.mpisup == 0) then
        if(associated(damping_scaling_values)) then
           write(*,*) "creating damping constraints, scaling values given: ",damping_scaling_values
           write(lu_out,*) "creating damping constraints, scaling values given: ",damping_scaling_values
        else
           write(*,*) "creating damping constraints; there are no scaling values given"
           write(lu_out,*) "creating damping constraints; there are no scaling values given"
        end if
     end if ! myrank == 0
!
     ! use the same error message as for initiating lmreg and adding smoothing (errmsg was not allocated after init(lmreg,...) )
     call addDamping(lmreg,errmsg,scaling_values=damping_scaling_values,neq_added=ndamp_added)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
!
     if(.myrank.mpisup == 0) then
        write(*,*) "there were ",ndamp_added," linear damping constraints created"
        write(lu_out,*) "there were ",ndamp_added," linear damping constraints created"
     end if ! myrank == 0
  end if ! add_damping

  if(add_regularization) then
     neq_lmreg_glob = .neq.lmreg
     if(.myrank.mpisup == 0) then
        write(*,*) "in total, there were ",neq_lmreg_glob," linear regularization constraints created"
        write(lu_out,*) "in total, there were ",neq_lmreg_glob," linear regularization constraints created"; write(lu_out,*) ""
     end if ! myrank == 0
     call dealloc(errmsg)
  else
     neq_lmreg_glob = 0
  end if
!
!------------------------------------------------------------------------
!  each local process computes which portion of kernel system and regularization equations to manage, proc 0 documents this subdivision in output file
!
  ! distribute the paths contained in dmspace evenly (as even as possible) among the procs
  npath_loc = npath_glob / .numtasks.mpisup ! every proc gets at least npath_glob/.numtasks.mpisup paths to handle
  n_remain = mod(npath_glob,.numtasks.mpisup) ! the remaining paths should be distributed among the first procs
  if(.myrank.mpisup+1 <= n_remain) then
     ! if I'm among the first n_remain procs, I'll take one more path
     npath_loc = npath_loc+1
     ipath1_loc = (.myrank.mpisup)*npath_loc + 1 ! all procs smaller than me must also have npath_loc=npath_loc+1 paths (i.e. one more); since .myrank.mpisup has offset 0 (i.e. 0,1,2,...), (.myrank.mpisup)*npath_loc is the number of paths that all procs smaller than me were assigned
     ipath2_loc = ipath1_loc + npath_loc -1
  else
     ipath1_loc = (.myrank.mpisup)*npath_loc + n_remain + 1 ! all the remaining paths were assigned to procs smaller than me, so shift my first index by the complete number of remaining paths. note, .myrank.mpisup has offset 0 (i.e. tasks are enumerated 0,1,2,...)
     ipath2_loc = ipath1_loc + npath_loc -1
  end if
  if(npath_loc == 0) then
     ndata_loc = 0 ! there are no data samples to be handled
     ipath1_loc = -1 ! artificially set first and last path index to -1 indicating that those values are not used
     ipath2_loc = -1
  else
     call copyDataModelSpaceInfo(dmspace_loc,dmspace_glob,ipath1=ipath1_loc,ipath2=ipath2_loc)
     ndata_loc = .ndata.dmspace_loc
  end if
  there_is_local_data = ndata_loc > 0
!
  ! distribute the regularization equations as evenly as possible among the procs
  if(add_regularization) then
     nregul_loc = neq_lmreg_glob / .numtasks.mpisup ! every proc gets at least neq_lmreg_glob/.numtasks.mpisup regularization equations to handle
     n_remain = mod(neq_lmreg_glob,.numtasks.mpisup) ! the remaining equations should be distributed among the first procs
     if(.myrank.mpisup+1 <= n_remain) then
        ! if I'm among the first n_remain procs, I'll take one more equation
        nregul_loc = nregul_loc+1
        iregul1_loc = (.myrank.mpisup)*nregul_loc + 1 ! all procs smaller than me must also have nregul_loc=nregul_loc+1 eqs (i.e. one more); since .myrank.mpisup has ooffset 0 (i.e. 0,1,2,...), (.myrank.mpisup)*nregul_loc is the number of eqs that all procs smaller than me were assigned
        iregul2_loc = iregul1_loc + nregul_loc -1
     else
        iregul1_loc = (.myrank.mpisup)*nregul_loc + n_remain + 1 !  all the remaining eqs were assigned to procs smaller than me, so shift my first index by the complete number of remaining eqs. note, .myrank.mpisup has offset 0 (i.e. tasks are enumerated 0,1,2,...)
        iregul2_loc = iregul1_loc + nregul_loc -1
     end if
     if(nregul_loc == 0) then
        ! in this case, iregul1_loc,iregul2_loc will have strange values
        ! actually it's not important to reset them, since they will not be used below when nregul_loc == 0
        ! but do it anyway
        iregul1_loc = 0
        iregul2_loc = 0
     end if
  else
     nregul_loc = 0
  end if ! add_regularization
  there_is_local_regularization = nregul_loc > 0
!
  there_is_only_local_data = there_is_local_data .and. (.not.there_is_local_regularization)
  there_is_only_local_regularization = there_is_local_regularization .and. (.not.there_is_local_data)
  there_is_local_data_and_regularization = there_is_local_data .and. there_is_local_regularization
!
  ! on root process, gather specifications on all local portions of data and regularization part of the system
  if(.myrank.mpisup == 0) allocate(ipath1_loc_procs(.numtasks.mpisup),ipath2_loc_procs(.numtasks.mpisup),&
       npath_loc_procs(.numtasks.mpisup),ndata_loc_procs(.numtasks.mpisup),&
       iregul1_loc_procs(.numtasks.mpisup),iregul2_loc_procs(.numtasks.mpisup),&
       nregul_loc_procs(.numtasks.mpisup))
  call barrier(mpisup)
  call MPI_GATHER(ipath1_loc,1,MPI_INTEGER,ipath1_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER ipath1_loc"
     goto 20
  end if
  call MPI_GATHER(ipath2_loc,1,MPI_INTEGER,ipath2_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER ipath2_loc"
     goto 20
  end if
  call MPI_GATHER(npath_loc,1,MPI_INTEGER,npath_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER npath_loc"
     goto 20
  end if
  call MPI_GATHER(ndata_loc,1,MPI_INTEGER,ndata_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER ndata_loc"
     goto 20
  end if
  call MPI_GATHER(iregul1_loc,1,MPI_INTEGER,iregul1_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER iregul1_loc"
     goto 20
  end if
  call MPI_GATHER(iregul2_loc,1,MPI_INTEGER,iregul2_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER iregul2_loc"
     goto 20
  end if
  call MPI_GATHER(nregul_loc,1,MPI_INTEGER,nregul_loc_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT GATHER nregul_loc"
     goto 20
  end if
  if(.myrank.mpisup == 0) then
     write(lu_out,*) "FOR ALL MPI PROCESSES DISPLAY NOW THE FOLLOWING LOCAL VALUES (rank, path range, range of "//&
          "regularization equations):"
     write(lu_out,*) "          myrank     ipath1     ipath2     npath     ndata     iregul1     iregul2     nregul"
     do j = 1,.numtasks.mpisup
        write(lu_out,*) j-1,ipath1_loc_procs(j),ipath2_loc_procs(j),npath_loc_procs(j),&
             ndata_loc_procs(j),iregul1_loc_procs(j),iregul2_loc_procs(j),nregul_loc_procs(j)
     end do
     write(lu_out,*) ""
     deallocate(ipath1_loc_procs,ipath2_loc_procs,npath_loc_procs,ndata_loc_procs,&
       iregul1_loc_procs,iregul2_loc_procs,nregul_loc_procs)
  end if
!
!------------------------------------------------------------------------
!  read in local portion of kernel matrix (and right-hand-side)
!
  if(.myrank.mpisup == 0) then
     write(*,*) "allocating local kernel linear system now"
  end if
  if(there_is_local_data) then
     !subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
     call new(errmsg,myname)
     call initiateSerialKernelLinearSystem(KLSE_loc,dmspace_loc,0,0,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
  end if
!
  if(.myrank.mpisup == 0) then
     write(*,*) "reading in kernel matrix now"
  end if
  allocate(minmaxval_K_loc_per_param(2,nparam_pmtrz))
  minmaxval_K_loc_per_param(1,:) = HUGEVAL
  minmaxval_K_loc_per_param(2,:) = -HUGEVAL
  if(there_is_local_data) then
     call new(errmsg,myname)
     lu1 = get(fuh)
     lu2 = get(fuh)
     call readMatrixSerialKernelLinearSystem(KLSE_loc,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
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
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     K_loc => .KM.KLSE_loc
     if(.not.associated(K_loc)) then
        write(*,*) "ERROR: K_loc not associated (rank ",.myrank.mpisup,")"
        goto 20
     end if
     if(.not.(size(K_loc,1)==ndata_loc .and. size(K_loc,2)==ncol_glob)) then
        write(*,*) "ERROR: K_loc has size ",size(K_loc,1),size(K_loc,2),"; expected size ",&
             ndata_loc,ncol_glob," (rank ",.myrank.mpisup,")"
        goto 20
     end if
!
     ! compute min and max values of local portion of kernel matrix per parameter of model parametrization
     nullify(pparam,idx_dmspace)
     do while (nextParamModelParametrization(parametrization,param_name,iparam_pmtrz))
        if(associated(pparam)) deallocate(pparam)
        if(associated(idx_dmspace)) deallocate(idx_dmspace)
        allocate(pparam(1)); pparam(1) = param_name
        idx_dmspace => getIndxModelValues(dmspace_loc,param=pparam)
        if(.not.associated(idx_dmspace)) cycle

        minmaxval_K_loc_per_param(1,iparam_pmtrz) = minval(K_loc(:,idx_dmspace))
        minmaxval_K_loc_per_param(2,iparam_pmtrz) = maxval(K_loc(:,idx_dmspace))
     end do ! nextParamModelParametrization
     if(associated(pparam)) deallocate(pparam)
     if(associated(idx_dmspace)) deallocate(idx_dmspace)
     nullify(pparam,idx_dmspace)
  else ! there_is_local_data
     nullify(K_loc)
  end if ! there_is_local_data
!
  if(.myrank.mpisup == 0) then
     write(*,*) "reading in measured data now"
  end if
  if(there_is_local_data) then
     !subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
     !  path_measured_data,lu,errmsg,ignore_data_weights,apply_mdata_normalization)
     call new(errmsg,myname)
     call readMeasuredDataSerialKernelLinearSystem(KLSE_loc,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg,&
          apply_mdata_normalization=normalize_data)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
  end if ! there_is_local_data
!
  if(.myrank.mpisup == 0) then
     write(*,*) "reading in synthetic data now"
  end if
  if(there_is_local_data) then
     ! subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
     !      nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
     !      apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
     !      ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
     call new(errmsg,myname)
     call readSyntheticDataSerialKernelLinearSystem(KLSE_loc,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
          ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
          trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
          apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
          path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
          apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
          path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
          apply_sdata_normalization=normalize_data)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
  end if ! there_is_local_data
!
! IN CASE OF PATH SPECIFIC MODELS, READ IN SYNTHETIC CORRECTIONS
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
     if(.myrank.mpisup == 0) then
        write(*,*) "reading in synthetics correction data now"
     end if
     if(there_is_local_data) then
        ! subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
        !      nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
        !      apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
        !      ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
        call new(errmsg,myname)
        call readSyntheticDataSerialKernelLinearSystem(KLSE_loc,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
             ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
             trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
             apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
             path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
             apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
             path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
             apply_sdata_normalization=normalize_data,read_synthetic_corrections=.true.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
     end if ! there_is_local_data
  end if ! USE_PATH_SPECIFIC_MODELS
!
!------------------------------------------------------------------------
!  now set the residuals vector
!  accounting for the cases of using path specific corrections or not
!  
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
     if(.myrank.mpisup == 0) then
        write(*,*) "setting right-hand-side to corrected data residuals now"
     end if
     if(there_is_local_data) then
        !subroutine setRhsAsCorrectedDataResidualKernelLinearSystem(this,errmsg)
        call new(errmsg,myname)
        call setRhsAsCorrectedDataResidualKernelLinearSystem(KLSE_loc,errmsg)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
     end if ! there_is_local_data
!
  else ! USE_PATH_SPECIFIC_MODELS
!
     if(.myrank.mpisup == 0) then
        write(*,*) "setting right-hand-side to data residuals now"
     end if
     if(there_is_local_data) then
        !subroutine setRhsAsDataResidualKernelLinearSystem(this,errmsg)
        call new(errmsg,myname)
        call setRhsAsDataResidualKernelLinearSystem(KLSE_loc,errmsg)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 20
        call dealloc(errmsg)
     end if ! there_is_local_data
  end if ! USE_PATH_SPECIFIC_MODELS
!
  if(there_is_local_data) then
     rhs_loc => .rhs.KLSE_loc
     if(.not.associated(rhs_loc)) then
        write(*,*) "ERROR: rhs_loc not associated (rank ",.myrank.mpisup,")"
        goto 20
     end if
     if(.not.(size(rhs_loc,1)==ndata_loc .and. size(rhs_loc,2)==1)) then
        write(*,*) "ERROR: rhs_loc has size ",size(rhs_loc,1),size(rhs_loc,2),"; expected size ",&
             ndata_loc,"      1  (rank ",.myrank.mpisup,")"
        goto 20
     end if
  else
     nullify(rhs_loc)
  end if ! there_is_local_data
!
!------------------------------------------------------------------------
!  get local portion of regularization constraints
!  first comunicate global values minmax_K_per_param
!
  allocate(minmaxval_K_glob_per_param(2,nparam_pmtrz))
  call barrier(mpisup)
!
  call MPI_ALLREDUCE(minmaxval_K_loc_per_param(1,:),minmaxval_K_glob_per_param(1,:),nparam_pmtrz,&
       MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT ALLREDUCE MIN K loc per param"
     goto 20
  end if
  call MPI_ALLREDUCE(minmaxval_K_loc_per_param(2,:),minmaxval_K_glob_per_param(2,:),nparam_pmtrz,&
       MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ios)
  if(ios.ne.MPI_SUCCESS) then
     write(*,*) "rank ",.myrank.mpisup," COULD NOT ALLREDUCE MAX K loc per param"
     goto 20
  end if
!
  if(.myrank.mpisup == 0) then
     write(*,*) "defining explicit local regularization equations now"
  end if

  if(there_is_local_regularization) then
     call new(errmsg,myname)
     call getEquations(lmreg,regeq_indx,regeq_coef,regeq_rhs,errmsg,&
          ieq_start=iregul1_loc,ieq_end=iregul2_loc,minmaxval_K_per_param=minmaxval_K_glob_per_param)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20 
     call dealloc(errmsg)
     if(.not.(associated(regeq_indx).and.associated(regeq_coef).and.associated(regeq_rhs))) then
        write(*,*) "rank ",.myrank.mpisup," ERROR: NO REGULARIZATION EQUATIONS RETURNED, EVEN IF THERE SHOULD BE"
        goto 20
     end if
     if(use_dble_cgls) then
        ! convert vector regeq_rhs to double precision allocating an extra variable regeq_rhs_dble
        allocate(regeq_rhs_dble(size(regeq_rhs)))
        regeq_rhs_dble = dble(regeq_rhs)
        deallocate(regeq_rhs)
     end if
  end if
  ! free the memory of all regularization constraints, the local constraints remain in memory in regeq_indx,regeq_coef,regeq_rhs
  call dealloc(lmreg)
!------------------------------------------------------------------------
!  solve kernel linear system now
!
  if(.myrank.mpisup == 0) then
     write(*,*) "solving linear system now by CG-least-squares algorithm, having ",ndata_glob+neq_lmreg_glob,&
          " rows and ",ncol_glob,"columns"
  end if

  ! Solve the linear system by calling subroutines which are ALL contained in this program (no modules used).
  ! On return, the above declared variable solution is allocated and contains the solution of the CG algorithm.
  ! Since the variable solution has kind CUSTOM_REAL, we need to convert to single precision when writing
  ! the inverted model output files after solving the system. The output model values in ASKI, namely 
  ! are handled in single precision only.

  if(there_is_local_data) then
     if(use_dble_cgls) then
        ! Create exact double precision copies of the kernel matrix K_loc and the residual vector rhs_loc (which 
        ! are just pointers to the arrays in object KLSE_loc). 
        allocate(K_loc_dble(size(K_loc,1),size(K_loc,2)),rhs_loc_dble(size(rhs_loc,1),size(rhs_loc,2)))
        K_loc_dble = dble(K_loc)
        rhs_loc_dble = dble(rhs_loc)

        ! Deallocate object KLSE_loc now, since the kernel matrix is probably 
        ! the biggest variable in memory.
        call dealloc(KLSE_loc)
        nullify(K_loc,rhs_loc)
     end if ! use_dble_cgls

  else ! there_is_local_data

     ! there is no data to handle on this proc, hence object KLSE should not have been created at all
     nullify(K_loc,rhs_loc)
  end if ! there_is_local_data


  ! SOLVE THE LINEAR SYSTEM
  call barrier(mpisup)
  call solveCgls()

  if(there_is_local_data) then
     if(use_dble_cgls) then
        ! deallocate the double precision copies as soon as possible
        deallocate(K_loc_dble,rhs_loc_dble)
     else
        ! deallocate kernel matrix (and rhs vector) as soon as possible by deallocating object KLSE_loc
        nullify(K_loc,rhs_loc)
        call dealloc(KLSE_loc)
     end if
  end if
!
!------------------------------------------------------------------------
!  create kernel_inverted_model object model update
!
  if(.myrank.mpisup == 0) then
!
     write(*,*) "creating model objects from solution vector and reference model and computing new model"
!
! create kernel_inverted_model object model update
!
     call new(errmsg,myname)
     call packVectorToKernelInvertedModel(kim_up_abs,real(solution),dmspace_glob,errmsg) ! no need to fork "if(use_dble_cgls)..." for both cases this is compilable and gives the correct result
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
! make a copy of kim_up_abs in order to produce both, absolute new model and relative model update
     call copyKernelInvertedModel(kim_new,kim_up_abs)
!
! write update to file(s) now
!
     !subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
     call new(errmsg,myname)
     call writeFileKernelInvertedModel(kim_up_abs,trim(outfile_base)//'_up_abs.kim',get(fuh),errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim_up_abs,.invgrid.iterbasics,&
          trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile_base)//'_up_abs',get(fuh),errmsg)
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
     call interpolateKernelReferenceToKernelInvertedModel(kim_ref,.krm.iterbasics,parametrization,&
          .invgrid.iterbasics,.intw.iterbasics,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
! create kernel_inverted_model object new abs model = abs update + reference by modifying object kim_new, which at this point contains values of kim_up_abs
! print new max/min model values
! write new model to files
!
     write(*,*) "updating model now, computing new absolute model values"
     write(lu_out,*) "updating model now, computing absolute new model values"
     !subroutine summateInstancesKernelInvertedModel(this,that,errmsg,c1,c2,relative)
     call new(errmsg,myname)
     call summateInstancesKernelInvertedModel(kim_new,kim_ref,errmsg)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
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
     call writeFileKernelInvertedModel(kim_new,trim(outfile_base)//'_new.kim',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim_new,.invgrid.iterbasics,&
          trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile_base)//'_new',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
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
     if (.level.errmsg == 2) goto 20
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
     call writeFileKernelInvertedModel(kim_ref,trim(outfile_base)//'_up_rel.kim',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     !subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim_ref,.invgrid.iterbasics,&
          trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile_base)//'_up_rel',get(fuh),errmsg)
     call undo(fuh)
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) goto 20
     call dealloc(errmsg)
!
     write(lu_out,*) "successfully written all output files"; write(lu_out,*) ""
!
  end if ! .myrank.mpisup == 0


  ! if program reaches regular end, indicate so in the output file
  if(.myrank.mpisup == 0) then
     write(lu_out,*) "everything went OK, good bye"; write(lu_out,*) ""
  end if


!------------------------------------------------------------------------
!  clean up
!
  if(.myrank.mpisup == 0) then
     close(lu_out)
     call add(fuh,lu_out)
     close(lu_plot)
     call add(fuh,lu_plot)
  end if
  call dealloc(mpisup)
!
!!$  if(allocated(pd)) deallocate(pd) ! etc
!!$  if(allocated()) deallocate()
!!$  if(allocated()) deallocate()
!!$  if(allocated()) deallocate()
!!$  if(allocated()) deallocate()
!
  call dealloc(dmspace_loc)
!
  if(allocated(ipath1_loc_procs)) deallocate(ipath1_loc_procs)
  if(allocated(ipath2_loc_procs)) deallocate(ipath2_loc_procs)
  if(allocated(npath_loc_procs)) deallocate(npath_loc_procs)
  if(allocated(ndata_loc_procs)) deallocate(ndata_loc_procs)
  if(allocated(iregul1_loc_procs)) deallocate(iregul1_loc_procs)
  if(allocated(iregul2_loc_procs)) deallocate(iregul2_loc_procs)
  if(allocated(nregul_loc_procs)) deallocate(nregul_loc_procs)

  if(allocated(minmaxval_K_loc_per_param)) deallocate(minmaxval_K_loc_per_param)
  if(allocated(minmaxval_K_glob_per_param)) deallocate(minmaxval_K_glob_per_param)
  if(associated(smoothing_scaling_values)) deallocate(smoothing_scaling_values)
  if(associated(damping_scaling_values)) deallocate(damping_scaling_values)
  if(associated(regeq_indx)) then
     do n = 1,size(regeq_indx)
        call dealloc(regeq_indx(n))
     end do
     deallocate(regeq_indx)
  end if
  if(associated(regeq_coef)) then
     do n = 1,size(regeq_coef)
        call dealloc(regeq_coef(n))
     end do
     deallocate(regeq_coef)
  end if
  if(associated(regeq_rhs)) deallocate(regeq_rhs)
  if(allocated(regeq_rhs_dble)) deallocate(regeq_rhs_dble)
!
  if(associated(pparam)) deallocate(pparam)
  if(associated(idx_dmspace)) deallocate(idx_dmspace)
  if(associated(idx_kim)) deallocate(idx_kim)
  if(associated(pcell)) deallocate(pcell)
!
  call dealloc(kim_up_abs); call dealloc(kim_new); call dealloc(kim_ref)
!
  call dealloc(dmspace_glob)
  call dealloc(invbasics); call dealloc(iterbasics)
  call dealloc(fuh)
  call dealloc(ap)
!
  ! clean stop, everythin went OK
  stop
!
  ! in case an error has occurred, better call blacs_abort (killing all processes), since the error may have occurred only 
  ! on a few but not all processes. Just 'stop' could cause some mpi processes to continue as phantom processes
20 if(.myrank.mpisup == 0) then
     close(lu_out)
  end if
  call abort(mpisup)
  stop
!
  ! in case an erro has occurred processing the command line arguments, only rank 0 prints errors and usage on screen
  ! and all others wait until all rank 0 has written all output to screen before program terminates
30 if(.myrank.mpisup == 0) then
     if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
     call usage(ap)
  end if
  call barrier(mpisup)
  call dealloc(mpisup)
  stop
!
!###########################################################################################
!###########################################################################################
!
contains

  subroutine solveCgls()

    call setup_transposed_regularization_equations()

    ! ALLOCATE ALL CGLS VECTORS

    if(there_is_local_data) allocate(rd(ndata_loc),qd(ndata_loc))
    if(there_is_local_regularization) allocate(rs(nregul_loc),qs(nregul_loc))
    if(there_is_local_data_and_regularization) allocate(sd_loc(ncol_glob))

    allocate(s_loc(ncol_glob),s(ncol_glob),p(ncol_glob))


    ! INITIATION AND FIRST ITERATION

    ! initiate convergence criterion
    call initate_convergence_criterion()

    if(use_nontrivial_starting_solution) then
       ! starting solution vector was read in above

       call r_is_b_minus_A_solution()
!
    else ! use_nontrivial_starting_solution
       ! in this case solution=x was initiated to 0.0, hence the vector r = b - Ax is trivially assigned r = b

       call r_is_b()

    end if ! .not.use_nontrivial_starting_solution

    call compute_norms_r_solution()
    if(.myrank.mpisup == 0) then
       write(*,*) "after iteration XXX  ;  residual norm ||r||  ;  sta  ;  lta   ;  solution norm ||x||"
       write(lu_plot,*) "after_iteration     residual_norm_||r||     sta_norm_||r||      lta_norm_||r||    solution_norm_||x||"

       ! print the norms of the starting solutions as "after iteration 0"
       i_iter = 0
       write(*,*) i_iter,norm_r,sta,lta,norm_solution
       write(lu_plot,*) i_iter,norm_r,sta,lta,norm_solution
    end if

! GATHER ALL PROCS BEFORE ITERATING
    call barrier(mpisup)

!
! FIRST ITERATION STARTS HERE
!
    i_iter = 1

    call s_is_A_transposed_r()
    call gamma_is_norm_s_squared()

    ! in the first iteration, there is no gamma/gamma_old yet to update p, instead p simply is s
    if(use_blas1) then
       if(use_dble_cgls) then
          call DCOPY(ncol_glob, s, 1, p, 1) 
       else
          call SCOPY(ncol_glob, s, 1, p, 1) 
       end if
    else ! use_blas1
       p = s
    end if ! use_blas1

    call q_is_A_p()
    call compute_norm_q_squared()

    call update_r_solution()
    call compute_norms_r_solution()

    ! if the iteration process already terminates here, skip the loop over the other iterations
    if(convergence_criterion_is_met()) goto 1

    if(.myrank.mpisup == 0) then
       write(*,*) i_iter,norm_r,sta,lta,norm_solution
       write(lu_plot,*) i_iter,norm_r,sta,lta,norm_solution
    end if

!
! ITERATE INVERSION SCHEME, SECOND ITERATION AND SO FORTH
!

    do i_iter = 2,max_niter

       call s_is_A_transposed_r()
       call gamma_is_norm_s_squared()

       call update_p()

       call q_is_A_p()
       call compute_norm_q_squared()

       call update_r_solution()
       call compute_norms_r_solution()

       if(convergence_criterion_is_met()) exit

       ! print this after calling convergence_criterion_is_met (since in that function, sta,lta are updated)
       if(.myrank.mpisup == 0) then
          write(*,*) i_iter,norm_r,sta,lta,norm_solution
          write(lu_plot,*) i_iter,norm_r,sta,lta,norm_solution
       end if
    end do ! i_iter


1   if(.myrank.mpisup == 0) then
       ! print this after the loop has terminated since in function convergence_criterion_is_met sta,lta are updated (so print here the last values)
       write(*,*) i_iter,norm_r,sta,lta,norm_solution
       write(lu_plot,*) i_iter,norm_r,sta,lta,norm_solution

       write(*,*) ""
       write(*,*) "conjugate gradient scheme terminated after ",i_iter," iterations"
       write(*,*) ""
       write(lu_out,*) ""
       write(lu_out,*) "conjugate gradient scheme terminated after ",i_iter," iterations"
       write(lu_out,*) ""
    end if

    ! in the very end, recompute the exact residuals (since in algorithm above, the residuals are only updated!)
    ! also recompute the norms of the residuals (data residuals, regularization residuals)
    call r_is_b_minus_A_solution()
    call compute_norms_r_solution()
    if(.myrank.mpisup == 0) then
       write(*,*) "since in the CG algorithm the residual vector is only updated for efficiency reasons, ",&
            "the exact residuals were recomputed here:"
       write(*,*) "squared residual of data part of the system ||(d-s) - Kx||^2 = ",norm_rd_squared
       write(*,*) "squared residual of regularization part of the system ||Sx||^2 = ",norm_rs_squared
       write(*,*) ""
       write(lu_out,*) "since in the algorithm, the residual vector is only updated for efficiency reasons, ",&
            "the exact residuals were recomputed here:"
       write(lu_out,*) "squared residual norm of data part of the system ||(d-s) - Kx||^2 = ",norm_rd_squared
       write(lu_out,*) "squared residual norm of regularization part of the system ||Sx||^2 = ",norm_rs_squared
       write(lu_out,*) ""
    end if


    ! DEALLOCATE EVERYTHING WHICH WAS ALLOCATED WITHIN THIS SUBROUTINE (except solution, of course)


    if(associated(regeq_T_indx)) then
       do n = 1,size(regeq_T_indx)
          call dealloc(regeq_T_indx(n))
       end do
       deallocate(regeq_T_indx)
    end if
    if(associated(regeq_T_coef)) then
       do n = 1,size(regeq_T_coef)
          call dealloc(regeq_T_coef(n))
       end do
       deallocate(regeq_T_coef)
    end if
    if(allocated(rd)) deallocate(rd)
    if(allocated(qd)) deallocate(qd)
    if(allocated(rs)) deallocate(rs)
    if(allocated(qs)) deallocate(qs)
    if(allocated(s_loc)) deallocate(s_loc)
    if(allocated(s)) deallocate(s)
    if(allocated(p)) deallocate(p)
    if(allocated(norm_r_array_sta)) deallocate(norm_r_array_sta)
    if(allocated(norm_r_array_lta)) deallocate(norm_r_array_lta)

  end subroutine solveCgls

!------------------------------------------------------------------------

  subroutine setup_transposed_regularization_equations()
    if(.not.there_is_local_regularization) return

    ! the transposed regularization matrix is ALWAYS single precision, 
    ! (just as regeq_coef is always single precision)
    allocate(regeq_T_indx(ncol_glob),regeq_T_coef(ncol_glob))

    do ieq = 1,nregul_loc
       ! get this regularization equation
       indx => getVectorPointer(regeq_indx(ieq))
       coef => getVectorPointer(regeq_coef(ieq))

       ! check if there is any inconsistency with this regularization equation
       if(.not.(associated(indx).and.associated(coef))) then
          write(*,*) "ERROR: local regularization equation ",ieq," is not associated, this should not happen!!, myrank ",&
               .myrank.mpisup
          call abort(mpisup)
       end if
       if(size(indx)/=size(coef) .or. size(indx) == 0) then
          write(*,*) "ERROR: in local regularization equation ",ieq,", sizes of indx and coef are not the same or 0, ",&
               "this should not happen!!, myrank ",.myrank.mpisup
          call abort(mpisup)
       end if

       ! for each coefficient in this regularization equation, define the corresponding entry in the transposed equations
       do j = 1,size(indx)
          ! indx(j) contains the model value index, the transposed equation of which needs to be appended
          ! by coefficient coef(j)
          imval = indx(j)

          ! get transposed equation corresponding to the model parameter index "indx(j)"
          indx_T => getVectorPointer(regeq_T_indx(imval))
          coef_T => getVectorPointer(regeq_T_coef(imval))
!
          if(associated(indx_T)) then
             n = size(indx_T)
          else
             n = 0
          end if
!
          ! extend this transposed equation by the pair of current equation index ieq and coefficient coef(j)
          indx_T => reallocate(indx_T,n+1)
          coef_T => reallocate(coef_T,n+1)

          indx_T(n+1) = ieq
          coef_T(n+1) = coef(j)

          ! re-associate reallocated pointers
          call associateVectorPointer(regeq_T_indx(imval),indx_T)
          call associateVectorPointer(regeq_T_coef(imval),coef_T) 
      end do ! j
    end do ! ieq

    ! NOTE THAT FOR THIS LOCAL PORTION OF THE REGULARIZATION MATRIX, SOME OF THE TRANSPOSED EQUATIONS
    ! MAY HAVE REMAINED UNDEFINED! IN THIS CASE THE MATRIX-VECTOR PRODUCT OF THAT ROW IS 0
    ! (however, other procs should have contribution to that entry, since overall the regularization matrix should have full rank)

  end subroutine setup_transposed_regularization_equations

!------------------------------------------------------------------------

  subroutine r_is_b()

    if(use_blas1) then

       if(there_is_local_data) then
          if(use_dble_cgls) then
             call DCOPY(ndata_loc, rhs_loc_dble(:,1), 1, rd, 1)
          else
             call SCOPY(ndata_loc, rhs_loc(:,1), 1, rd, 1)
          end if
       end if ! there_is_local_data

       if(there_is_local_regularization) then
          if(use_dble_cgls) then
             call DCOPY(nregul_loc, regeq_rhs_dble, 1, rs, 1)
          else
             call SCOPY(nregul_loc, regeq_rhs, 1, rs, 1)
          end if
       end if ! there_is_local_regularization

    else ! use_blas1

       if(there_is_local_data) then
          if(use_dble_cgls) then
             rd = rhs_loc_dble(:,1)
          else
             rd = rhs_loc(:,1)
          end if
       end if ! there_is_local_data

       if(there_is_local_regularization) then
          if(use_dble_cgls) then
             rs = regeq_rhs_dble
          else
             rs = regeq_rhs
          end if
       end if

    end if ! use_blas1

  end subroutine r_is_b

!------------------------------------------------------------------------

  subroutine s_is_A_transposed_r()

    if(there_is_only_local_data) then

       ! compute data Part of A^T r, i.e.  K_loc^T rd
       ! use variable s_loc for the result (this is already the final result, as there is no regularization part to add)
       if(use_blas2) then
          if(use_dble_cgls) then
             call DGEMV('Transpose', ndata_loc,ncol_glob, 1._CUSTOM_REAL , K_loc_dble, ndata_loc, rd, 1, 0._CUSTOM_REAL, s_loc, 1)
          else
             call SGEMV('Transpose', ndata_loc,ncol_glob, 1._CUSTOM_REAL , K_loc, ndata_loc, rd, 1, 0._CUSTOM_REAL, s_loc, 1)
          end if
       else ! use_blas2
          if(use_dble_cgls) then
             s_loc = matmul(transpose(K_loc_dble),rd)
          else
             s_loc = matmul(transpose(K_loc),rd)
          end if
       end if ! use_blas2

    elseif(there_is_only_local_regularization) then
       
       ! compute regularization Part of A^T r, i.e. S^T rs
       call s_loc_is_regeq_transpose_rs()

    elseif(there_is_local_data_and_regularization) then

       ! compute data Part of A^T r, i.e.  K_loc^T rd
       ! use variable sd_loc for the result, add to regularization part below
       if(use_blas2) then
          if(use_dble_cgls) then
             call DGEMV('Transpose', ndata_loc,ncol_glob, 1._CUSTOM_REAL , K_loc_dble, ndata_loc, rd, 1, 0._CUSTOM_REAL, sd_loc, 1)
          else
             call SGEMV('Transpose', ndata_loc,ncol_glob, 1._CUSTOM_REAL , K_loc, ndata_loc, rd, 1, 0._CUSTOM_REAL, sd_loc, 1)
          end if
       else ! use_blas2
          if(use_dble_cgls) then
             sd_loc = matmul(transpose(K_loc_dble),rd)
          else
             sd_loc = matmul(transpose(K_loc),rd)
          end if
       end if ! use_blas2
       
       ! compute regularization Part of A^T r, i.e. S^T rs
       call s_loc_is_regeq_transpose_rs()

       ! add to s_loc the result of the data part sd_loc
       if(use_blas1) then
          ! BLAS1 routine for operation s_loc = 1.0*sd_loc + s_loc
          if(use_dble_cgls) then
             call DAXPY(ncol_glob, 1._CUSTOM_REAL, sd_loc, 1, s_loc, 1)
          else
             call SAXPY(ncol_glob, 1._CUSTOM_REAL, sd_loc, 1, s_loc, 1)
          end if
       else ! use_blas1
          s_loc = s_loc + sd_loc
       end if ! use_blas1

    else

       ! there is neither data nor regularization to be handled on this local proc, 
       ! hence set s_loc = 0. in order to account for proper ALLREDUCE operation below
       s_loc = 0._CUSTOM_REAL

    end if

    !! FS FS
    ! The following used to be a check when testing the code. For some unrealistic test scenarios there
    ! happened to be NaNs in the vectors. 
    ! If you notice that this slows down your calculations significantly, uncomment this check.
    if(any(isnan(s_loc))) then
       write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  I have NaNs in s_loc"
       call abort(mpisup)
    end if

    ! ALLREDUCE s_loc to s
    call MPI_ALLREDUCE(s_loc,s,ncol_glob,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE s IN ITERATION ",i_iter
       call abort(mpisup)
    end if

    !! FS FS
    ! The following used to be a check when testing the code. For some unrealistic test scenarios there
    ! happened to be NaNs in the vectors. 
    ! If you notice that this slows down your calculations significantly, uncomment this check.
    if(any(isnan(s))) then
       write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  I have NaNs in s"
       call abort(mpisup)
    end if
  end subroutine s_is_A_transposed_r

!------------------------------------------------------------------------

  subroutine s_loc_is_regeq_transpose_rs()
    if(.not.there_is_local_regularization) return

    do imval = 1,ncol_glob

       indx_T => getVectorPointer(regeq_T_indx(imval))
       coef_T => getVectorPointer(regeq_T_coef(imval))

       if(associated(indx_T)) then
          ! DO NOT USE BLAS HERE, VERY SMALL ARRAY SIZE TO SUM AND NOT EVENLY DISTRIBUTED IN MEMORY
          ! I GUESS BLAS WOULD NOT BE EFFICIENT HERE (Florian Schumacher, May 2014)
          s_loc(imval) = sum( coef_T * rs(indx_T) )
       else
          s_loc(imval) = 0.
       end if
    end do ! imval
  end subroutine s_loc_is_regeq_transpose_rs

!------------------------------------------------------------------------

  subroutine gamma_is_norm_s_squared()
    ! backup current value of gamma for computing beta = gamma/gamma_old
    gamma_old = gamma

    if(use_blas1) then
       ! BLAS1 function for single precision dot product
       if(use_dble_cgls) then
          gamma = DDOT(ncol_glob, s, 1, s, 1)
       else
          gamma = SDOT(ncol_glob, s, 1, s, 1)
       end if
    else ! use_blas1
       gamma = sum(s*s)
    end if ! use_blas1

    !! FS FS
    ! The following used to be a check when testing the code. For some unrealistic test scenarios there
    ! happened to be NaNs in the vectors. 
    ! If you notice that this slows down your calculations significantly, uncomment this check.
    if(isnan(gamma)) then
       write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  gamma is NaN"
       call abort(mpisup)
    end if
  end subroutine gamma_is_norm_s_squared

!------------------------------------------------------------------------

  subroutine update_p()
    if(use_blas1) then
       ! SINCE THERE IS NO BLAS ROUTINE like "xAXPBY", need to scale p first by beta=gamma/gamma_old, before adding to s
       if(use_dble_cgls) then
          call DSCAL(ncol_glob, gamma/gamma_old, p, 1) ! p = (gamma/gamma_old)*p
          call DAXPY(ncol_glob, 1._CUSTOM_REAL, s, 1, p, 1) ! p = 1.0*s + p
       else
          call SSCAL(ncol_glob, gamma/gamma_old, p, 1) ! p = (gamma/gamma_old)*p
          call SAXPY(ncol_glob, 1._CUSTOM_REAL, s, 1, p, 1) ! p = 1.0*s + p
       end if
    else ! use_blas1
       p = s + (gamma/gamma_old)*p
    end if ! use_blas1
  end subroutine update_p

!------------------------------------------------------------------------

  subroutine q_is_A_p()
    ! if there are data equations on this proc, compute qd = K_loc p 
    if(there_is_local_data) then
       if(use_blas2) then
          ! BLAS2 routine for qd = K_loc p
          if(use_dble_cgls) then
             call DGEMV('No transpose', ndata_loc, ncol_glob, 1._CUSTOM_REAL, K_loc_dble, ndata_loc, p, 1 , 0._CUSTOM_REAL, qd, 1)
          else
             call SGEMV('No transpose', ndata_loc, ncol_glob, 1._CUSTOM_REAL, K_loc, ndata_loc, p, 1 , 0._CUSTOM_REAL, qd, 1)
          end if
       else ! use_blas2
          if(use_dble_cgls) then
             qd = matmul(K_loc_dble,p)
          else
             qd = matmul(K_loc,p)
          end if
       end if ! use_blas2

       !! FS FS
       ! The following used to be a check when testing the code. For some unrealistic test scenarios there
       ! happened to be NaNs in the vectors. 
       ! If you notice that this slows down your calculations significantly, uncomment this check.
       if(any(isnan(qd))) then
          write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  I have NaNs in vector qd"
          call abort(mpisup)   
       end if

    end if

    ! if there are regularization equations on this proc, compute qs = regeq p 
    if(there_is_local_regularization) then
       do ieq = 1,nregul_loc

          indx => getVectorPointer(regeq_indx(ieq))
          coef => getVectorPointer(regeq_coef(ieq))

if(.not.(associated(indx).and.associated(coef))) then
   write(*,*) "########## rank ",.myrank.mpisup,", iteration,ieq ",i_iter,ieq," :  indx or coef not associated"
   call abort(mpisup)   
end if

          ! in subroutine setup_transposed_regularization_equations it was already checked that
          ! both, indx and coef are associated and are of same length /= 0
          ! so no need to check this again here inside CG loop

          ! DO NOT USE BLAS HERE, VERY SMALL ARRAY SIZE TO SUM AND NOT EVENLY DISTRIBUTED IN MEMORY
          ! I GUESS BLAS WOULD NOT BE EFFICIENT HERE (but guess only!!) (Florian Schumacher, May 2014)
          qs(ieq) = sum( coef * p(indx) )
       end do ! ieq

       !! FS FS
       ! The following used to be a check when testing the code. For some unrealistic test scenarios there
       ! happened to be NaNs in the vectors. 
       ! If you notice that this slows down your calculations significantly, uncomment this check.
       if(any(isnan(qs))) then
          write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  I have NaNs in vector qs"
          call abort(mpisup)   
       end if

    end if
  end subroutine q_is_A_p

!------------------------------------------------------------------------

  subroutine compute_norm_q_squared()

    if(there_is_local_data) then
       if(use_blas1) then
          ! BLAS1 function for single precision dot product
          if(use_dble_cgls) then
             norm_qd_squared = DDOT(ndata_loc, qd, 1, qd, 1)
          else
             norm_qd_squared = SDOT(ndata_loc, qd, 1, qd, 1)
          end if
       else ! use_blas1
          norm_qd_squared = sum( qd * qd )
       end if ! use_blas1
    else
       norm_qd_squared = 0.
    end if

    if(there_is_local_regularization) then
       if(use_blas1) then
          ! BLAS1 function for single precision dot product
          if(use_dble_cgls) then
             norm_qs_squared = DDOT(nregul_loc, qs, 1, qs, 1)
          else
             norm_qs_squared = SDOT(nregul_loc, qs, 1, qs, 1)
          end if
       else ! use_blas1
          norm_qs_squared = sum( qs * qs )
       end if ! use_blas1
    else
       norm_qs_squared = 0.
    end if

    ! total sum of squares of all local entries of vector q
    cr_tmp = norm_qd_squared + norm_qs_squared

    !! FS FS
    ! The following used to be a check when testing the code. For some unrealistic test scenarios there
    ! happened to be NaNs in the vectors. 
    ! If you notice that this slows down your calculations significantly, uncomment this check.
    if(isnan(cr_tmp)) then
       write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  cr_tmp is NaN in compute_norm_q_squared()"
       call abort(mpisup)   
    end if

    ! ALLREDUCE cr_tmp to norm_q_squared in order to get the sum of squares of entries in vector q
    call MPI_ALLREDUCE(cr_tmp,norm_q_squared,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE norm_q_squared IN ITERATION ",i_iter
       call abort(mpisup)
    end if

    !! FS FS
    ! The following used to be a check when testing the code. For some unrealistic test scenarios there
    ! happened to be NaNs in the vectors. 
    ! If you notice that this slows down your calculations significantly, uncomment this check.
    if(isnan(norm_q_squared)) then
       write(*,*) "########## rank ",.myrank.mpisup,", iteration, ",i_iter," :  norm_q_squared is NaN in compute_norm_q_squared()"
       call abort(mpisup)   
    end if

  end subroutine compute_norm_q_squared

!------------------------------------------------------------------------

  subroutine update_r_solution()

    alpha = gamma/norm_q_squared
    if(alpha <= 0. .or. isnan(alpha)) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," THERE IS NO UPDATE OF THE SOLUTION IN ITERATION ",&
           i_iter,": (alpha = ",alpha,"), probably dividing through by 0 (or by Infinity)?! (norm_q_squared = ",&
           norm_q_squared,"), "
       if(use_dble_cgls) then
          write(*,*) "          you are already using double precision, so I don't know what to do. You could ",&
               "try to scale up (or down) all values in your linear system"
       else
          write(*,*) "          this could be caused by numbers being too small (or too large) for single ",&
               "precision: I strongly recommend to use double precision (change the appropriate line at the ",&
               "beginning of code file solveCglsKernelSystem.f90 and recompile program solveCglsKernelSystem)"
       end if
       call abort(mpisup)
    end if

    if(recompute_r_at_all) then
       recompute_r = mod(i_iter,niter_recompute_r) == 0
    else
       recompute_r = .false.
    end if

    if(recompute_r) then
       ! RECOMPUTE RESIDUAL VECTOR INSTEAD OF UPDATING IT
       call r_is_b_minus_A_solution()

    else ! recompute_r

       ! update local residual vectors rd = rd - alpha*qd , rs = -alpha*qs
       if(there_is_local_data) then
          if(use_blas1) then
             if(use_dble_cgls) then
                call DAXPY(ndata_loc, -alpha, qd, 1, rd, 1) ! rd = -alpha*qd + rd
             else
                call SAXPY(ndata_loc, -alpha, qd, 1, rd, 1) ! rd = -alpha*qd + rd
             end if
          else ! use_blas1
             rd = rd - alpha*qd
          end if ! use_blas1
       end if
       if(there_is_local_regularization) then
          if(use_blas1) then
             if(use_dble_cgls) then
                call DAXPY(nregul_loc, -alpha, qs, 1, rs, 1) ! rs = -alpha*qs + rs
             else
                call SAXPY(nregul_loc, -alpha, qs, 1, rs, 1) ! rs = -alpha*qs + rs
             end if
          else ! use_blas1
             rs = rs - alpha*qs
          end if ! use_blas1
       end if
    end if ! recompute_r

    ! update solution vector
    if(use_blas1) then
       if(use_dble_cgls) then
          call DAXPY(ncol_glob, alpha, p, 1, solution, 1) ! solution = alpha*p + solution
       else
          call SAXPY(ncol_glob, alpha, p, 1, solution, 1) ! solution = alpha*p + solution
       end if
    else ! use_blas1
       solution = solution + alpha*p
    end if ! use_blas1

  end subroutine update_r_solution

!------------------------------------------------------------------------

  subroutine r_is_b_minus_A_solution()
    ! compute rd = rhs_loc(:,1) - K_loc*solution
    !         rs = regeq_rhs - S_loc*solution , where S_loc symbolizes the regularization equation matrix (do it by explicit summation)
!
    if(there_is_local_data) then
       if(use_blas2) then
          if(use_blas1) then
             if(use_dble_cgls) then
                call DCOPY(ndata_loc, rhs_loc_dble(:,1), 1, rd, 1)
             else
                call SCOPY(ndata_loc, rhs_loc(:,1), 1, rd, 1)
             end if
          else ! use_blas1
             if(use_dble_cgls) then
                rd = rhs_loc_dble(:,1)
             else
                rd = rhs_loc(:,1)
             end if
          end if ! use_blas1
          if(use_dble_cgls) then
             call DGEMV('No-Transpose', ndata_loc,ncol_glob, -1._CUSTOM_REAL , K_loc_dble, ndata_loc, &
                  solution, 1, 1._CUSTOM_REAL, rd,1)
          else
             call SGEMV('No-Transpose', ndata_loc,ncol_glob, -1._CUSTOM_REAL , K_loc, ndata_loc, &
                  solution, 1, 1._CUSTOM_REAL, rd, 1)
          end if
       else ! use_blas2
          if(use_dble_cgls) then
             rd = rhs_loc_dble(:,1) - matmul(K_loc_dble,solution)
          else
             rd = rhs_loc(:,1) - matmul(K_loc,solution)
          end if
       end if ! use_blas2
    end if ! there_is_local_data
!
    if(there_is_local_regularization) then
       do ieq = 1,nregul_loc
          indx => getVectorPointer(regeq_indx(ieq))
          coef => getVectorPointer(regeq_coef(ieq))
          if(.not.(associated(indx).and.associated(coef))) then
             write(*,*) "########## rank ",.myrank.mpisup,", iteration,ieq ",i_iter,ieq," :  indx or coef not associated"
             call abort(mpisup)   
          end if
          ! in subroutine setup_transposed_regularization_equations it was already checked that
          ! both, indx and coef are associated and are of same length /= 0
          ! so actually no need to check this again here!
          ! DO NOT USE BLAS HERE, VERY SMALL ARRAY SIZE TO SUM AND NOT EVENLY DISTRIBUTED IN MEMORY
          ! I GUESS BLAS WOULD NOT BE EFFICIENT HERE (Florian Schumacher, May 2014)
          if(use_dble_cgls) then
             rs(ieq) = regeq_rhs_dble(ieq) - sum( dble(coef) * solution(indx) )
          else
             rs(ieq) = regeq_rhs(ieq) - sum( coef * solution(indx) )
          end if
       end do ! ieq
    end if ! there_is_local_regularization
  end subroutine r_is_b_minus_A_solution

!------------------------------------------------------------------------

  subroutine compute_norms_r_solution()

    ! COMPUTE THE NORM OF r
    ! FIRST SUM UP SCALAR PRODUCTS OF ALL PARTS OF VECTOR r, AT LAST TAKE SQUARE-ROOT
    if(there_is_local_data) then
       if(use_blas1) then
          ! BLAS1 function for single precision dot product
          if(use_dble_cgls) then
             norm_rd_squared_loc = DDOT(ndata_loc, rd, 1, rd, 1)
          else
             norm_rd_squared_loc = SDOT(ndata_loc, rd, 1, rd, 1)
          end if
       else ! use_blas1
          norm_rd_squared_loc = sum(rd*rd)
       end if ! use_blas1
    else
       norm_rd_squared_loc = 0.
    end if

    if(there_is_local_regularization) then
       if(use_blas1) then
          ! BLAS1 function for single precision dot product
          if(use_dble_cgls) then
             norm_rs_squared_loc = DDOT(nregul_loc, rs, 1, rs, 1)
          else
             norm_rs_squared_loc = SDOT(nregul_loc, rs, 1, rs, 1)
          end if
       else ! use_blas1
          norm_rs_squared_loc = sum(rs*rs)
       end if ! use_blas1
    else
       norm_rs_squared_loc = 0.
    end if

    ! ALLREDUCE norm_rd_squared_loc to norm_rd_squared in order to get the sum of squares of entries in vector r
    ! we also then have independent information on the residuals of the data equations, knowing norm_rd_squared
    call MPI_ALLREDUCE(norm_rd_squared_loc,norm_rd_squared,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE norm_rd_squared IN ITERATION ",i_iter
       call abort(mpisup)
    end if

    ! ALLREDUCE norm_rs_squared_loc to norm_rs_squared in order to get the sum of squares of entries in vector r
    ! we also then have independent information on the residuals of the regularization equations, knowing norm_rs_squared
    call MPI_ALLREDUCE(norm_rs_squared_loc,norm_rs_squared,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE norm_rs_squared IN ITERATION ",i_iter
       call abort(mpisup)
    end if

    ! this is used for sta,lta check if the solution of the linear system has converged.
    norm_r = sqrt(norm_rd_squared + norm_rs_squared)


    ! COMPUTE THE NORM OF solution
    if(use_blas1) then
       ! BLAS1 function for single precision 2-norm
       if(use_dble_cgls) then
          norm_solution = DNRM2(ncol_glob, solution, 1)
       else
          norm_solution = SNRM2(ncol_glob, solution, 1)
       end if
    else ! use_blas1
       norm_solution = sqrt( sum(solution*solution) )
    end if ! use_blas1

  end subroutine compute_norms_r_solution

!------------------------------------------------------------------------

  subroutine initate_convergence_criterion()

    ! DEFINE TOLERANCE VALUES atol, btol FOR NORM CRITERION IN FUNCTION convergence_criterion_is_met
    atol = 1.0
    btol = 1.0

    ! COMPUTE NORM OF COMPLETE rhs-VECTOR b OF REGULARIZED LINEAR SYSTEM

    cr_tmp = 0._CUSTOM_REAL

    ! sum of squares of local data part of rhs
    if(there_is_local_data) then
       if(use_blas1) then
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + DDOT(ndata_loc, rhs_loc_dble(:,1), 1, rhs_loc_dble(:,1), 1)
          else
             cr_tmp = cr_tmp + SDOT(ndata_loc, rhs_loc(:,1), 1, rhs_loc(:,1), 1)
          end if
       else ! use_blas1
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + sum(rhs_loc_dble(:,1)*rhs_loc_dble(:,1))
          else
             cr_tmp = cr_tmp + sum(rhs_loc(:,1)*rhs_loc(:,1))
          end if
       end if ! use_blas1
    end if
    if(there_is_local_regularization) then
       if(use_blas1) then
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + DDOT(nregul_loc, regeq_rhs_dble, 1, regeq_rhs_dble, 1)
          else
             cr_tmp = cr_tmp + SDOT(nregul_loc, regeq_rhs, 1, regeq_rhs, 1)
          end if
       else ! use_blas1
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + sum(regeq_rhs_dble*regeq_rhs_dble)
          else
             cr_tmp = cr_tmp + sum(regeq_rhs*regeq_rhs)
          end if
       end if ! use_blas1
    end if ! there_is_local_regularization

    ! ALLREDUCE SUM OF SQUARES OF GLOBAL rhs-VECTOR
    call MPI_ALLREDUCE(cr_tmp,btol_norm_b,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE norm_b"
       call abort(mpisup)
    end if
    btol_norm_b = btol*sqrt(btol_norm_b)

    ! COMPUTE FROBENIUS NORM OF COMPLETE REGULARIZED SYSTEM MATRIX

    cr_tmp = 0._CUSTOM_REAL

    ! sum of squares of local data part of rhs
    if(there_is_local_data) then
       if(use_blas1) then
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + DDOT(ndata_loc*ncol_glob, K_loc_dble, 1, K_loc_dble, 1) ! interpret K_loc as a column-wise stored vector
          else
             cr_tmp = cr_tmp + SDOT(ndata_loc*ncol_glob, K_loc, 1, K_loc, 1) ! interpret K_loc as a column-wise stored vector
          end if
       else ! use_blas1
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + sum(K_loc_dble*K_loc_dble)
          else
             cr_tmp = cr_tmp + sum(K_loc*K_loc)
          end if
          ! ! build the sum of squares over the columns separaterly
          ! ! does this perform better than sum(K_loc*K_loc) ??
          ! ! (usually the size of the columns will be smaller than the size of the rows, i.e. the number of columns)
          ! ! woud additionally need to fork if(use_dble_cgls)... here!
          ! do j = 1,ncol_glob
          !    cr_tmp + cr_tmp + sum(K_loc(:,j)*K_loc(:,j))
          ! end do ! j
       end if ! use_blas1
    end if
    if(there_is_local_regularization) then
       do ieq = 1,nregul_loc
          coef => getVectorPointer(regeq_coef(ieq))
          ! in subroutine setup_transposed_regularization_equations it was already checked that
          ! all coef are associated so no need to check this again here
          if(use_dble_cgls) then
             cr_tmp = cr_tmp + sum(dble(coef)*dble(coef))
          else
             cr_tmp = cr_tmp + sum(coef*coef)
          end if
       end do ! ieq
    end if
    ! ALLREDUCE SUM OF SQUARES OF GLOBAL MATRIX ENTRIES
    call MPI_ALLREDUCE(cr_tmp,atol_norm_A,1,CUSTOM_MPI_TYPE,MPI_SUM,MPI_COMM_WORLD,ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT ALLREDUCE norm_A"
       call abort(mpisup)
    end if
    atol_norm_A = atol*sqrt(atol_norm_A)

    if(.myrank.mpisup == 0) write(*,*) "BEFORE ITERATIONS: atol*norm_A,btol*norm_b = ",atol_norm_A ,btol_norm_b
    if(.myrank.mpisup == 0) write(lu_out,*) "BEFORE ITERATIONS: atol*norm_A,btol*norm_b = ",atol_norm_A ,btol_norm_b


    ! INITIATE sta/lta TRIGGER

    allocate(norm_r_array_sta(nsta),norm_r_array_lta(nlta))
    norm_r_array_sta = 0. ; norm_r_array_lta = 0.
    sta = 0.; lta = 0.
    max_nlta_nsta = max(nlta,nsta)
    !eps = epsilon(sta+lta) ! sta+lta gives the largest kind of sta and lta
    eps_single_precision = epsilon(1.0)

    if(.myrank.mpisup == 0) write(*,*) "BEFORE ITERATIONS: single-precision epsilon for termination check is ",eps_single_precision
    if(.myrank.mpisup == 0) write(lu_out,*) "BEFORE ITERATIONS: single-precision epsilon for termination check is ",&
         eps_single_precision
    if(.myrank.mpisup == 0) write(lu_out,*) ""
  end subroutine initate_convergence_criterion

!------------------------------------------------------------------------

  logical function convergence_criterion_is_met()
    ! ignore this routine for now, introduce a proper convergence criterion in due time

    convergence_criterion_is_met = .false.

    terminate =  (i_iter == max_niter)
    if(terminate) then
       if(.myrank.mpisup == 0) write(*,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
            ": maximum number of iterations reached"
       if(.myrank.mpisup == 0) write(lu_out,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
            ": maximum number of iterations reached"
    end if
    ! update function output 
    convergence_criterion_is_met = convergence_criterion_is_met .or. terminate

    !convergence_criterion_is_met = convergence_criterion_is_met .or. (norm_r <= btol_norm_b + atol_norm_A * norm_solution) ! DO NOT YET USE THIS CRITERION UNTIL IT IS FULY UNDERSTOOD
    !! FOR BETTER UNDERSTANDING, OUTPUT THE RESPECTIVE NUMBERS OF THE CRITERION
    !if(.myrank.mpisup == 0) write(*,*) "    btol*||b|| + atol*||A||*||x|| = ",btol_norm_b+atol_norm_a*norm_solution
    !if(.myrank.mpisup == 0) write(lu_out,*) "    btol*||b|| + atol*||A||*||x|| = ",btol_norm_b+atol_norm_a*norm_solution 

    ! ADDITIONALLY CHECK WHETHER THERE IS A TERMINATION REQUEST BY FILE outfile_base_TERMINATE.txt
    if(.myrank.mpisup == 0) then
       terminate = .false.
       ! master process checks whether the file exists 
       inquire(file=trim(outfile_base)//'_TERMINATE.txt',exist=terminate)
    end if
    ! master process broadcasts the result of the TERMINATE request query to the other processes
    call MPI_BCAST(terminate, 1, MPI_LOGICAL, 0 , MPI_COMM_WORLD , ios)
    if(ios.ne.MPI_SUCCESS) then
       write(*,*) "ERROR: rank ",.myrank.mpisup," COULD NOT BROADCAST external terminate request"
       call abort(mpisup)
    end if
    ! check if criterion is met
    if(terminate) then
       if(.myrank.mpisup == 0) write(*,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
            ": external termination was requested by existence of file '",trim(outfile_base)//"_TERMINATE.txt'"
       if(.myrank.mpisup == 0) write(lu_out,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
            ": external termination was requested by existence of file '",trim(outfile_base)//"_TERMINATE.txt'"
    end if
    ! update function output 
    convergence_criterion_is_met = convergence_criterion_is_met .or. terminate


    ! AT LAST CHECK IF THE RESIDUAL NORM DOES NOT CHANGE ANYMORE

    ! circularly loop over the entries in arrays norm_r_array_sta,norm_r_array_lta
    ! replacing the oldest entry with the current entry of norm_r
    ! define index "pointers" ista (ilta) pointing to the oldest entry in the arrays
    ista = mod(i_iter,nsta) + 1 ! ista starts with 2, but that doesn't matter... (circular loop)
    norm_r_array_sta(ista) = norm_r ! replace oldest value by current value or norm_r
    sta = sum(norm_r_array_sta) / dble(nsta) ! compute current value short-term average

    ! r_tmp = norm_r_array_sta( mod(max(i_iter-nsta,1),nsta) + 1 ) ! this is the oldest value in array norm_array_sta (by assumption the largest one)
    ! err_sta = 1.e8*norm_r/r_tmp ! norm_r is the newest value in array norm_r_array_sta, by assumpition the smallest one
    ! err_sta = (err_sta - aint(err_sta))*1.e-8 ! represents the digits of norm_r/r_tmp smaller than 10^-8
    ! err_sta = err_sta * r_tmp!*nsta/nsta  !! upper estimate absolute error done by computing the average value sta

    ! now do the same as for sta for the long-term average
    ilta = mod(i_iter,nlta) + 1
    norm_r_array_lta(ilta) = norm_r
    lta_old = lta ! backup the lta value of the last iteration, used for safty criterion
    lta = sum(norm_r_array_lta) / dble(nlta)

    ! r_tmp = norm_r_array_lta( mod(max(i_iter-nlta,1),nlta) + 1 ) ! this is the oldes value in array norm_array_lta (by assumption the largest one)
    ! err_lta = 1.e8*norm_r/r_tmp ! norm_r is the newest value in array norm_r_array_lta, by assumpition the smallest one
    ! err_lta = (err_lta - aint(err_lta))*1.e-8 ! represents the digits of norm_r/r_tmp smaller than 10^-8
    ! err_lta = err_lta * r_tmp!*nlta/nlta  !! upper estimate absolute error done by computing the average value lta

    ! do not use values sta,lta before the arrays norm_r_array_sta, norm_r_array_lta are
    ! completely filled (one complete circular loop for both, achieved after max(nlta,sta) iterations)
    if(i_iter > max_nlta_nsta) then
       ! always expect sta < lta, assuming monotonically decreasing values of norm_r
       ! terminate if sta is very close to lta in terms of single precision, i.e. sta + a_bit >= lta <=> sta/lta + eps_single_precision >= 1; note that eps_single_precision is the smallest real such that 1 + eps > 1
       terminate = (sta/lta + eps_single_precision > 1.d0)
       if(terminate) then
          if(.myrank.mpisup == 0) write(*,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
               ": short-term average of residual norm has reached the value of its long-term average (in single precision)"
          if(.myrank.mpisup == 0) write(lu_out,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
               ": short-term average of residual norm has reached the value of its long-term average (in single precision)"
       end if
       ! update function output 
       convergence_criterion_is_met = convergence_criterion_is_met .or. terminate

       ! ! terminate if sta is very close to lta, accounting for errors in the summation
       ! ! assume that for two summations the precision is eps, i.e. | a + b - true_sum| <= eps
       ! ! Hence, | sum(norm_r_array_sta) - true_sum | <= nsta*eps
       ! ! Hence, | sum(norm_r_array_sta)/nsta - true_sum/nsta | <= eps
       ! ! Use + eps for the smaller value to compare (sta) and -eps for the larger value to compare
       ! terminate = sta + err_sta + eps >= lta - err_lta !! DOES THIS CRITERION MAKE SENSE (do the error bounds computed above even make sense?)

       ! as a sort of safety criterion (in case sta/lta trigger is not configured 
       ! correctly), always terminate if lta is not monotonically decreasing anymore
       ! in order to prevent norm_r to leave minimum and explode
       terminate = lta > lta_old 
       if(terminate) then
          if(.myrank.mpisup == 0) write(*,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
               ": long-term average of residual norm is not monotonically decreasing ",&
               "anymore (violates assumption of generally monotonically decreasing residual norm)"
          if(.myrank.mpisup == 0) write(lu_out,*) "CONVERGENCE CRITERION IS MET IN ITERATION ",i_iter,&
               ": long-term average of residual norm is not monotonically decreasing ",&
               "anymore (violates assumption of generally monotonically decreasing residual norm)"
       end if
       convergence_criterion_is_met = convergence_criterion_is_met .or. terminate

    end if ! i_iter > max_nlta_nsta

  end function convergence_criterion_is_met


end program solveCglsKernelSystem
