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
!> \brief set up and handle kernel matrix system by serial processing (not parallel)
!!
!! \details Given an object of type data_model_space_info, the kernel matrix as well as
!!  the right-hand-side of the linear system can be set up, possibly adding regularization
!!  equations or regularizing columns to the kernel matrix.
!!  The set-up linear system can be solved serially by module serialLinearSystem. The 
!!  computed solution can be related to the original entries in the data_model_space_info
!!  object. An overall weighted sensitivity may be computed. 
!!
!! \author Florian Schumacher
!! \date July 2015
!
module kernelLinearSystem
!
  use serialLinearSystem
  use dataModelSpaceInfo
  use modelParametrization
  use parameterCorrelation
  use realloc
  use seismicEvent
  use seismicStation
  use asciiDataIO
  use spectralWaveformKernel
  use vectorPointer
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateKernelLinearSystem; end interface
  interface getSumColumns; module procedure getSumColumnsKernelLinearSystem; end interface
  interface setColumnRegularization; module procedure setColumnRegKernelLinearSystem; end interface
  interface operator (.KM.); module procedure getMatrixKernelLinearSystem; end interface
  interface operator (.sol.); module procedure getSolutionKernelLinearSystem; end interface
  interface operator (.rhs.); module procedure getRhsKernelLinearSystem; end interface
  interface operator (.nrow.); module procedure getNrowKernelLinearSystem; end interface
  interface operator (.ndata.); module procedure getNdataKernelLinearSystem; end interface
  interface operator (.nmval.); module procedure getNmvalKernelLinearSystem; end interface
  interface operator (.ncol.); module procedure getNcolKernelLinearSystem; end interface
  interface operator (.nrowreg.); module procedure getNrowregKernelLinearSystem; end interface
  interface operator (.md.); module procedure getMeasuredDataKernelLinearSystem; end interface
  interface operator (.sd.); module procedure getSyntheticDataKernelLinearSystem; end interface
  interface operator (.scd.); module procedure getSyntheticCorrectionDataKernelLinearSystem; end interface
  interface operator (.dmspace.); module procedure getDmspaceKernelLinearSystem; end interface
  interface isInitiated; module procedure isInitiatedKernelLinearSystem; end interface
!
!> \brief kernel linear system object
  type kernel_linear_system
     private
     !
     logical :: initiated = .false. !< flag to indicate whether the object was initiated by calling initiateSerialKernelLinearSystem
     !
     real, dimension(:,:), pointer :: K => null() !< kernel system matrix, including regularization conditions (if any)
     !
     type (data_model_space_info) :: dmspace !< complete data and model space corresponding to the kernel system
     !
     integer :: nrow = 0 !< (local) number of rows of kernel matrix, nrow=ndata+nrowreg
     integer :: ndata = 0 !< (local) number of rows reserved for kernels or synthetic data (first rows of the system)
     integer :: nrowreg = 0 !< (local) number rows for data regularization conditions (smoothing, damping of rows) that were added to the system (last rows of the system)
     integer :: ncol = 0 !< (local) number of columns of kernel matrix, ncol = npar+ncolreg
     integer :: nmval = 0 !< (local) number of columns reserved for kernel values on inversion grid cells (corresponding to model values)
     integer :: ncolreg = 0 !< (local) number of columns for column regularization (e.g. in kernel focussing, last columns of the matrix)
     !
     character(len=1) :: transpose = '' !< either 'N' or 'T': indicates whether K*sol=rhs is solved (N) or (K^T)*sol=rhs is solved (T). THE MEANING OF rhs AND sol IS DEPENDENT ON transpose!!
     !
     integer :: nrhs = 0 !< number of right-hand-side vectors (columns of array rhs)
     real, dimension(:,:), pointer :: rhs => null() !< right-hand-side vector(s) of kernel system (has size(nrow,nrhs) if 'N', or size(ncol,nrhs) if 'T')
     real, dimension(:,:), pointer :: sol => null() !< solution vector array (has size(ncol,nrhs) if 'N', or size(nrow,nrhs) if 'T')
     !
     real, dimension(:), pointer :: mdata => null() !< vector of measured data according to data space
     real, dimension(:), pointer :: sdata => null() !< vector of synthetic data according to data space
     real, dimension(:), pointer :: scdata => null() !< vector of synthetic corrections according to data space
  end type kernel_linear_system
!
contains
!------------------------------------------------------------------------
!> \brief allocate system matrix, according to number of data, model values and regularization conditions
!! \param this kernel linear system object
!! \param ndata number of data in data space, which will be used to fill kernel matrix
!! \param nmval number of model values in model space, which will be used to fill kernel matrix
!! \param nrowreg number of regularization conditions that will be added to the rows of the system
!! \param errmsg error message
!
  subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
    type (kernel_linear_system) :: this
    type (data_model_space_info) :: dmspace
    integer :: nrowreg,ncolreg
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=32) :: myname = 'initiateSerialKernelLinearSystem'
    integer :: ndata,nmval
!
    call addTrace(errmsg,myname)
!
    if(this%initiated) then
       call add(errmsg,2,"this object is already initiated. deallocated first, before initiating anew",myname)
       return
    end if
!
    ndata = .ndata.dmspace
    nmval = .nmval.dmspace
!
    if(ndata .le. 0) then
       write(errstr,*) "requested number of data (rows) = ",ndata,&
            ", is invalid: must be positive numbers. Incoming data space is empty."
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(nmval .le. 0) then
       write(errstr,*) "requested number of model values (columns of Kernel matrix) = ",nmval,&
            "; There are no columns of the kernel matrix, hence this object will only be usable ",&
            "to read in measured or synthetic data"
       call add(errmsg,1,errstr,myname)
    end if
!
    if(nrowreg < 0) then
       write(errstr,*) "requested number of row regularization conditions = ",nrowreg," must not be negative"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(ncolreg < 0) then
       write(errstr,*) "requested number of column regularization conditions = ",ncolreg," must not be negative"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    write(errstr,*) "incoming: ndata in dmspace = ",ndata,"; nmval in dmspace = ",nmval,"; nrowreg = ",nrowreg,&
         "; ncolreg = ",ncolreg,"; i.e. nrow of kernel system = ",ndata + nrowreg,"; ncol of kernel system = ",nmval + ncolreg
    call add(errmsg,0,errstr,myname)
!
    this%nrow = ndata + nrowreg
    this%ndata = ndata
    this%nrowreg = nrowreg
    this%ncol = nmval + ncolreg
    this%nmval = nmval
    this%ncolreg = ncolreg
    call copyDataModelSpaceInfo(this%dmspace,dmspace)
!
    this%initiated = .true.
!
  end subroutine initiateSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief setup kernel matrix (not parallelized)
!! \details loop over paths and read in only kernel frequencies which are needed for current path, 
!!  then place the kernel values for those data samples in the respective rows of the kernel matrix,
!!  which are in same order as the data samples. Account for possible correlation of model parameters
!!  and apply event and station filters, and normalize if requested.
!! \param this kernel system
!! \param df_measured_data frequency step of measured data
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param path_sensitivity_kernels path where sensitivity kernel files are located
!! \param ntot_invgrid total number of inversion grid cells of the inversion grid in use (used to check kernel files)
!! \param pcorr parameter correlation object which contains information about parameter correlation (it also knows whether a correlation should be applied at all)
!! \param lu1 first file unit
!! \param lu2 second file unit
!! \param errmsg error message
!! \param apply_event_filter optional logical flag whether event filters are to be accounted for
!! \param path_event_filter optional path were event filter files are (must be present if apply_event_filter == .true., otherwise ignored)
!! \param apply_station_filter optional logical flag whether station filters are to be accounted for
!! \param path_station_filter optional path were station filter files are (must be present if apply_station_filter == .true., otherwise ignored)
!! \param ignore_data_weights optional logical flag whether to ignore the data weights contained in data space
!! \param apply_kernel_normalization optional logical flag whether to normalize the kernel values by respective normalization factors contained in data space
!
  subroutine readMatrixSerialKernelLinearSystem(this,df_measured_data,&
       nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
       lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
       ignore_data_weights,apply_kernel_normalization)
    ! incoming
    type (kernel_linear_system) :: this
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
    character(len=35) :: myname = 'readMatrixSerialKernellLinearSystem'
    type (error_message) :: errmsg2,errmsg3
    integer :: j,n,ndata_count,ndata_div_50,status
    logical, dimension(:), allocatable :: idata_added,imval_added
    ! data space
    integer :: ipath,idata,icomp
    character(len=character_length_evid) :: evid,evid_old
    character(len=character_length_staname) :: staname,staname_old
    character(len=character_length_component), dimension(:), pointer :: comp,pcomp
    character(len=2), dimension(:), pointer :: pimre
    integer, dimension(:), pointer :: idata_path,idata_path_ifreq,idata_path_ifreq_comp,pifreq,ifreq
    real, dimension(:), pointer :: pwdata
    real :: wdata
    logical :: ignore_wdata,normalize_kernel,ldummy
    real, dimension(:), pointer :: normalization_factors
    ! model space
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name,param_name_corr
    integer :: nparam,iparam,kparam
    character(len=character_length_param), dimension(:), pointer :: pparam,param_dmspace,param_corr
    character(len=character_length_param), dimension(:), allocatable :: param
    integer, dimension(:), pointer :: idx_mval,idx_cell,pcell
    type (integer_vector_pointer), dimension(:,:), allocatable :: idx_mval_cell
    ! filtering
    complex, dimension(:), pointer :: event_filter,station_comp_filter
    logical :: do_event_filtering,do_station_filtering,do_filtering
    integer :: ifilter
    complex :: filter_value
    ! kernels
    character(len=500) :: file_kernel
    type (spectral_waveform_kernel) :: kernel
    integer :: jfreq,ifreq_current
    complex, dimension(:,:), pointer :: kernel_values
    complex, dimension(:), allocatable :: filtered_kernel_values
    logical :: apply_parameter_correlation
    real :: c_corr
!
    nullify(comp,pcomp,pimre,idata_path,idata_path_ifreq,idata_path_ifreq_comp,pifreq,ifreq,pwdata,normalization_factors,&
         pparam,param_dmspace,param_corr,idx_mval,idx_cell,pcell,event_filter,station_comp_filter,kernel_values)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before reading in kernel matrix",myname)
       return
    end if
    if(this%nmval == 0) then
       call add(errmsg,2,"kernel linear system was initiated without any model values (only to read in data). "//&
            "cannot read kernel matrix",myname)
       return
    end if
    if(associated(this%K)) then
       call add(errmsg,2,"kernel matrix seems to be read in already. please deallocate linear system first",myname)
       return
    else
       allocate(this%K(this%nrow,this%ncol),stat=status)
       if(status/=0) then
          write(errstr,*) "could not allocate kernel matrix for ",this%nrow," rows and ",this%ncol,&
               " columns, allocate status = ",status
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    do_event_filtering = .false.
    if(present(apply_event_filter)) then
       if(apply_event_filter) then
          if(present(path_event_filter)) then
             do_event_filtering = .true.
          else
             call add(errmsg,2,&
                  "event filtering is requested by flag apply_event_filter, but no path of event filters is given ",myname)
             return
          end if
       end if
    end if
    if(do_event_filtering) then
       call add(errmsg,0,"incoming path of event filters : '"//trim(path_event_filter)//"'",myname)
    else
       call add(errmsg,0,"no incoming path of event filters, hence will NOT apply any event filters",myname)
    end if
!
    do_station_filtering = .false.
    if(present(apply_station_filter)) then
       if(apply_station_filter) then
          if(present(path_station_filter)) then
             do_station_filtering = .true.
          else
             call add(errmsg,2,&
                  "station filtering is requested by flag apply_station_filter, but no path of station filters is given ",myname)
             return
          end if
       end if
    end if
    if(do_station_filtering) then
       call add(errmsg,0,"incoming path of station filters: '"//trim(path_station_filter)//"'",myname)
    else
       call add(errmsg,0,"no incoming path of station filters, hence will NOT apply any station filters",myname)
    end if
!
    do_filtering = do_event_filtering .or. do_station_filtering
!
    call add(errmsg,0,"incoming path of spectral sensitivity kernels: '"//trim(path_sensitivity_kernels)//"'",myname)
!
    if(present(ignore_data_weights)) then
       ignore_wdata = ignore_data_weights
    else
       ! by default do not ignore, but apply the data weights!
       ignore_wdata = .false.
    end if
    if(ignore_wdata) then
       call add(errmsg,0,"WILL IGNORE data weights defined in data space",myname)
    else
       call add(errmsg,0,"will apply data weights defined in data space normally (will NOT ignore them)",myname)
    end if
!
    if(present(apply_kernel_normalization)) then
       normalize_kernel = apply_kernel_normalization
    else
       ! by default do not normalize the kernel values
       normalize_kernel = .false.
    end if
    if(normalize_kernel) then
       call add(errmsg,0,"WILL NORMALIZE kernel values as defined in data space",myname)
    else
       call add(errmsg,0,"will NOT normalize kernel values",myname)
    end if
!
    if(nfreq_measured_data.le.0) then
       call add(errmsg,2,"number frequencies of the measured data is less or equal to zero",myname)
       return
    end if
    if(nfreq_measured_data /= size(ifreq_measured_data)) then
       write(errstr,*) "number frequencies of the measured data ( = ",nfreq_measured_data,&
            ") does not match the size of the vector of frequency indices ( = ",size(ifreq_measured_data),")"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    parametrization = .pmtrz.this%dmspace
!
    ! get all different model parameters contained in model space
    param_dmspace => allParam(this%dmspace)
    if(.not.associated(param_dmspace)) then
       write(errstr,*) "no parameters returned from model space object: this error should not occur!"
       call add(errmsg,2,errstr,myname)
       return
    end if
    nparam = size(param_dmspace)
    apply_parameter_correlation = correlateAnyParameters(pcorr)
    if(apply_parameter_correlation) then
       ! define all ADDITIONAL parameters which are correlated to some parameter in the model space, 
       ! but which are NOT in the model space (also need kernel values for those)
       j = 0
       n = 0
       param_corr => null()
       do iparam = 1,nparam
          param_name = param_dmspace(iparam)
          if(correlateAnyParameters(pcorr,param_name)) then
             j = j+1
          else
             cycle
          end if
          do while(nextCorrelationParameter(pcorr,param_name,param_name_corr))
             ! only memorize the parameter param_name_corr, if it is not already in the model space
             if(.not.any(param_dmspace == param_name_corr)) then
                param_corr => reallocate(param_corr,n+1)
                param_corr(n+1) = param_name_corr
                n = n+1
             end if
          end do ! while nextCorrleation
       end do ! iparam
       if(associated(param_corr)) then
          ! increase the total number of parameters which must be accounted for in this routine
          nparam = nparam + size(param_corr)
       end if
       if(j>0) then
          write(errstr,*) "there are ",j," parameters in the model space which will be correlated to others "//&
               "(which should not be in the model space). a respective parameter correlation will be applied"
          call add(errmsg,0,errstr,myname)
       else
          call add(errmsg,1,"even though parameter correlation should be used in general, there is no "//&
               "parameter in the specific current model space for which any correlation applies.",myname)
       end if
    else
       call add(errmsg,0,"there is no parameter correlation to be applied",myname)
    end if
!
    ! finally define information arrays for all required parameters 
    ! (the order of the parameter names corresponds to the general order of parameters, possibly with gaps)
    allocate(param(nparam),idx_mval_cell(2,nparam))
    iparam = 0
    do while(nextParamModelParametrization(parametrization,param_name))
       if(any(param_dmspace == param_name)) then
          iparam = iparam + 1
          param(iparam) = param_name
          allocate(pparam(1)); pparam(1) = param_name
          idx_mval => getIndxModelValues(this%dmspace,param=pparam,cell=pcell)
          ! since param_name is in model space, idx_mval should be associated here (don't check again)
          call associateVectorPointer(idx_mval_cell(1,iparam),idx_mval)
          call associateVectorPointer(idx_mval_cell(2,iparam),pcell)
          nullify(idx_mval,pcell)
          deallocate(pparam)
       else if(apply_parameter_correlation .and. associated(param_corr)) then
          if(any(param_corr == param_name)) then
             ! for those parameters, the pointers in idx_mval_cell stay unassociated by intention
             iparam = iparam + 1
             param(iparam) = param_name
          end if
       end if
    end do ! while(nextParam)
    deallocate(param_dmspace)
    if(associated(param_corr)) deallocate(param_corr)
!
    ! allocate space for filtered kernel values. this is also used for correlated (and filtered) kernel values
    allocate(filtered_kernel_values(ntot_invgrid))
!
    ! as a security measure, check whether the model values defined by the model space are consistent (no double 
    ! indices, i.e. every idx_mval occurs exactly once)
    allocate(imval_added(this%nmval))
    imval_added(:) = .false.
    do iparam = 1,nparam
       idx_mval => getVectorPointer(idx_mval_cell(1,iparam))
       ! if there are no model values in model space for this parameter (i.e. parameter is only used for correlation) just cycle
       if(.not.associated(idx_mval)) cycle
!
       ! check if all returned model value indices are valid
       if(any(idx_mval>this%nmval .or. idx_mval<1)) then
          write(errstr,*) "there are invalid model value indices for ",iparam,"'th parameter '",trim(param(iparam)),&
               "' in model space: smaller than 1 or larger than nmval = ",&
               this%nmval,". data_model_space_info object may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 1          
       end if
       ! check if those model value indices idx_mval already occurred before (for different param_name)
       if(any(imval_added(idx_mval))) then
          call add(errmsg,2,"model space contains same model value indices for different model parameters. "//&
               "data_model_space_info object may be corrupt",myname)
          goto 1
       end if
       ! memorize these model value indices
       imval_added(idx_mval) = .true.
    end do ! iparam
    ! check if there are some model value indices which did not occur in one of the idx_mval above
    if(any(.not.imval_added)) then
       write(errstr,*) "model space did not correctly return all model values in expected index range from 1 to ",&
            this%nmval,". data_model_space_info object may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    deallocate(imval_added)
!
    ! get the kernel normalization factors, if requested
    if(normalize_kernel) then
       call getNormalization(this%dmspace,(/ (idata, idata=1,this%ndata) /),errmsg,norm_kernel=normalization_factors)
       if(.level.errmsg == 2) goto 1
    end if ! normalize_kernel
!
    ! keep track on which rows were added to the kernel matrix (to detect rows which are not initiated on exit (those rows might equal zero on exit, which is bad))
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
    ndata_count = 0
    ndata_div_50 = this%ndata / 50 + 1 ! give at most 50 messages of progress
!
    nullify(event_filter,station_comp_filter)
    evid_old = ''; staname_old = ''
    nullify(pifreq,pcomp,pimre,pwdata)
    nullify(idata_path,comp,ifreq,idata_path_ifreq,idata_path_ifreq_comp)
!
    ! loop over paths
    ipath = 0
    do while(nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp,all_ifreq=ifreq))
       ipath = ipath + 1
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       write(errstr,*) "adding kernel values to system for ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
       call add(errmsg2,0,errstr,myname)
!
       ! read in event_filter for this path, if necessary
       if(do_event_filtering .and. evid_old /= evid) then
          if(associated(event_filter)) deallocate(event_filter)
          errmsg3 = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu1,event_filter,ndata=nfreq_measured_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) call print(errmsg3)
          if(.level.errmsg3 == 2) then
             call add(errmsg2,2,"error reading event filter from formatted file '"//trim(path_event_filter)//&
                  "filter_"//trim(evid)//"'",myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") reading event filter"
             call add(errmsg,2,errstr,myname)
             goto 2
          endif
          call dealloc(errmsg3)
          ! remember from which event the filter values are
          evid_old = evid
       end if
!
       ! initiate kernel object
       call initiateSpectralWaveformKernel(kernel,parametrization,param,ntot_invgrid,comp,errmsg2,kernel_on_wp=.false.)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") initiating kernel object"
          call add(errmsg,2,errstr,myname)
          call print(errmsg2)
          goto 2
       end if
!
       ! open kernel file to read
       file_kernel = trim(path_sensitivity_kernels)//"spectral_kernel_"//trim(parametrization)//&
            "_"//trim(evid)//"_"//trim(staname)
       call initialReadSpectralWaveformKernel(kernel,file_kernel,lu1,errmsg2)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") initially reading kernel file '"//trim(file_kernel)//"'"
          call add(errmsg,2,errstr,myname)
          call print(errmsg2)
          goto 2
       end if
       call add(errmsg2,0,"successfully openend kernel file '"//trim(file_kernel)//"' to read",myname)
!
       ! check df of kernel
       if( (df_measured_data-.df.kernel)/df_measured_data > 1.e-4) then
          write(errstr,*) "the frequency step of the frequencies in the kernel file (",.df.kernel,&
               ") differs significantly from the frequency step of the measured data (",df_measured_data,&
               "), which suggests that the kernels were computed w.r.t. a different frequency discretization"
          call add(errmsg2,2,errstr,myname)
          call print(errmsg2)
          write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
!
       ! iterate over all frequencies contained in the data space for this path
       do jfreq = 1,size(ifreq)
          ifreq_current = ifreq(jfreq)

          if(associated(pifreq)) deallocate(pifreq)
          if(associated(idata_path_ifreq)) deallocate(idata_path_ifreq)
          allocate(pifreq(1)); pifreq(1) = ifreq_current
          idata_path_ifreq => getIndxDataSamples(this%dmspace,ifreq=pifreq,indx_in=idata_path)
          ! idata_path_ifreq should be associated here, don't check. 
          ! otherwise function nextPathDataModelSpaceInfo is corrupt
          write(errstr,*) "there are ",size(idata_path_ifreq),&
               " data samples in data space for this path and frequency index ",ifreq_current
          call add(errmsg2,0,errstr,myname)
!
          ! find correct index of filter value for this frequency
          if(do_filtering) then
             ifilter = -1
             do j = 1,nfreq_measured_data
                if(ifreq_measured_data(j) == ifreq_current) then
                   ifilter = j
                   exit
                end if
             end do ! j
             if(ifilter .le. 0) then
                write(errstr,*) "frequency index ",ifreq_current,", of this path is not contained in the ",&
                     "vector of frequency indices of measured data = ",ifreq_measured_data
                call add(errmsg2,2,errstr,myname)
                call print(errmsg2)
                write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                     ") mapping current frequency index to measured data indices"
                call add(errmsg,2,errstr,myname)
                goto 2
             end if
          end if ! do_filtering
!
          ! read kernel for this frequency
          call readSpectralWaveformKernel(kernel,ifreq_current,errmsg2)
          if(.level.errmsg2 ==  2) then
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") reading frequency ",ifreq_current," of kernel"
             call add(errmsg,2,errstr,myname)
             call print(errmsg2)
             goto 2
          end if
          write(errstr,*) "successfully read frequency ",ifreq_current," of kernel"
          call add(errmsg2,0,errstr,myname)
!
          do icomp = 1,size(comp)
             kernel_values => getValuesByCompSpectralWaveformKernel(kernel,comp(icomp))
!
             if(associated(pcomp)) deallocate(pcomp)
             if(associated(pimre)) deallocate(pimre); nullify(pimre)
             if(associated(pwdata)) deallocate(pwdata); nullify(pwdata)
             if(associated(idata_path_ifreq_comp)) deallocate(idata_path_ifreq_comp)
             
             allocate(pcomp(1)); pcomp(1) = comp(icomp)
             idata_path_ifreq_comp => getIndxDataSamples(this%dmspace,comp=pcomp,imre=pimre,wdata=pwdata,indx_in=idata_path_ifreq)
!print *, "sample:: evid = '"//trim(evid)//"', staname = '"//trim(staname)//"', jf = ",jf,&
!", comp = '"//trim(comp(icomp))//"', idata_path_ifreq = ",idata_path_ifreq,", idata_path_ifreq_comp = ",idata_path_ifreq_comp
             ! idata_path_ifreq_comp should be associated here, don't check. 
             ! otherwise function nextPathDataModelSpaceInfo is corrupt
             write(errstr,*) "there are ",size(idata_path_ifreq_comp),&
                  " data samples in data space for this path and frequency and component '"//trim(comp(icomp))//"'"
             call add(errmsg2,0,errstr,myname)
!
             ! in a sensible setup, size(idata_path_ifreq_comp) == 2, one 'im' sample and one 're' sample
             do idata=1,size(idata_path_ifreq_comp)
                ! check, if for any reasons these rows were already filled
                if(idata_added(idata_path_ifreq_comp(idata))) then
                   write(errstr,*) idata_path_ifreq_comp(idata),"'th row of kernel matrix was already filled with ",&
                        "kernel values. data_model_space_info object may be corrupt"
                   call add(errmsg2,2,errstr,myname)
                   call print(errmsg2)
                   write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
                   call add(errmsg,2,errstr,myname)
                   goto 2
                end if
             end do ! idata
!
!
             ! neutralize the data weights, if they are to be ignored
             if(ignore_wdata) pwdata = 1.0
!
             ! read in station_comp_filter for this path and component
             if(do_station_filtering) then
                if(associated(station_comp_filter)) deallocate(station_comp_filter)
                errmsg3 = readAsciiData(trim(path_station_filter)//"filter_"//trim(staname)//"_"//trim(comp(icomp)),&
                     lu2,station_comp_filter,ndata=nfreq_measured_data)
                call addTrace(errmsg3,myname)
                if(.level.errmsg3 /= 0) call print(errmsg3)
                if(.level.errmsg3 == 2) then
                   call add(errmsg2,2,"error reading station filter from formatted file '"//trim(path_station_filter)//&
                        "filter_"//trim(staname)//"_"//trim(comp(icomp))//"'",myname)
                   call print(errmsg2)
                   write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                        ") reading station filter for component '"//trim(comp(icomp))//"'"
                   call add(errmsg,2,errstr,myname)
                   goto 2
                endif
                call dealloc(errmsg3)
             end if ! do_station_filtering
!
             ! define filter value for this path,frequency and component
             if(do_filtering) then
                filter_value = (1.,0.)
                if(do_event_filtering) filter_value = filter_value * event_filter(ifilter)
                if(do_station_filtering) filter_value = filter_value * station_comp_filter(ifilter)
             end if
!
             ! now, for these particular rows of the kernel matrix which all depend on the same complex kernel 
             ! values (usually two rows, one for real and one for imaginary part), 
             ! loop on all parameters and see
             ! if they occur in the model space. If yes, the corresponding model value indices and 
             ! inversion grid cell indices are contained in vector pointer array idx_mval_cell, that 
             ! was initiated above (If no, then the parameter is used in correlation below).
             ! inside the parameter loop, apply the complex filtering before distinguishing between the real-valued rows
             do iparam = 1,nparam
                idx_mval => getVectorPointer(idx_mval_cell(1,iparam))
                idx_cell => getVectorPointer(idx_mval_cell(2,iparam))
                ! if this parameter is only correlated to another one but is actually not in the model space, cycle
                if(.not.associated(idx_mval)) cycle

                ! start off with the plain kernel values for this parameter
                filtered_kernel_values = kernel_values(:,iparam)

                ! then, secondly, correlate this parameter to any other parameters (if any correlation)
                if(apply_parameter_correlation .and. correlateAnyParameters(pcorr,param(iparam))) then
                   ! loop on all other parameters and check whether the current parameter should be correlated to 
                   ! those. If so, add the kernel values of those parameters to the temporary correlated kernel values
                   do kparam = 1,nparam
                      if(correlateParameters(pcorr,param(iparam),param(kparam),c_corr)) &
                           filtered_kernel_values = filtered_kernel_values + c_corr*kernel_values(:,kparam)
                   end do ! kparam
                end if ! any correlation

                ! thirdly, apply complex-valued filters (if any filtering should be applied) to bring syntetics to the form of the data
                if(do_filtering) filtered_kernel_values = filter_value * filtered_kernel_values
!
                ! at last, loop on the specific real-valued rows of the kernel matrix.
                ! in a sensible setup, size(idata_path_ifreq_comp) == 2, one 'im' sample and one 're' sample
                do idata=1,size(idata_path_ifreq_comp)

                   ! weight for this datum, including possible normalization
                   if(normalize_kernel) then
                      wdata = pwdata(idata)*normalization_factors(idata_path_ifreq_comp(idata))
                   else
                      wdata = pwdata(idata)
                   end if

                   ! apply weighting of this particular datum
                   select case (pimre(idata))
                   case('im')
                      this%K(idata_path_ifreq_comp(idata),idx_mval) = &
                           wdata * imag( filtered_kernel_values(idx_cell) )
                   case('re')
                      this%K(idata_path_ifreq_comp(idata),idx_mval) = &
                           wdata * real( filtered_kernel_values(idx_cell) )
                   end select

                end do ! idata
!
             end do ! iparam
!
             !ndata_count = ndata_count + size(idata_path_ifreq_comp)
             !if(mod(ndata_count,ndata_div_50)==0) write(*,*) "read serial kernel matrix row ", &
             !     ndata_count," out of",this%ndata,"; ",real(ndata_count)/real(this%ndata)*100.,"%"

             ! remember these rows were filled
             do idata=1,size(idata_path_ifreq_comp)
                idata_added(idata_path_ifreq_comp(idata)) = .true.
             end do
!
          end do ! icomp
!
       end do ! jfreq
!
       call finalReadSpectralWaveformKernel(kernel)
       call dealloc(kernel)
!
       if(.level.errmsg2 /= 0) call print(errmsg2)
       call dealloc(errmsg2)
    end do ! while(nextPath)
!
    ! check if there are any rows which were not filled
    if(.not.all(idata_added)) then
       write(errstr,*) "maybe some data samples occur more than once in the data space (not allowed!), as ",&
            count(.not.idata_added)," rows of the kernel matrix were not filled with values"
            !,", namely these ones (by index) :",pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(pimre)) deallocate(pimre)
    if(associated(pwdata)) deallocate(pwdata)
    if(associated(idata_path_ifreq)) deallocate(idata_path_ifreq)
    if(associated(idata_path_ifreq_comp)) deallocate(idata_path_ifreq_comp)
    if(associated(normalization_factors)) deallocate(normalization_factors)
!
    if(allocated(idx_mval_cell)) then
       do iparam=1,nparam
          call dealloc(idx_mval_cell(1,iparam))
          call dealloc(idx_mval_cell(2,iparam))
       end do
       deallocate(idx_mval_cell)
    end if
!
    if(associated(event_filter)) deallocate(event_filter)
    if(associated(station_comp_filter)) deallocate(station_comp_filter)
    if(allocated(idata_added)) deallocate(idata_added)
    if(allocated(imval_added)) deallocate(imval_added)
    if(allocated(param)) deallocate(param)
    if(allocated(filtered_kernel_values)) deallocate(filtered_kernel_values)
!
    ! if code comes here, return normally
    return
!
    ! if there is an error inside the loop on nextPath, need to reset the iterator (and deallocate all the pointers, etc.)
2   ldummy = nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp,all_ifreq=ifreq,reset=.true.)
    goto 1
  end subroutine readMatrixSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief copy kernel matrix from another kernelLinearSystem object
!! \details 
!!  
!! \param this kernel linear system for which the kernel matrix should be filled by copying from that%K
!! \param that kernel linear system whose kernel matrix is copied to this%K
!! \param errmsg error message
!
  subroutine copyMatrixKernelLinearSystem(this,that,errmsg)
    type (kernel_linear_system) :: this,that
    type (error_message) :: errmsg
    ! local
    character(len=28) :: myname = 'copyMatrixKernelLinearSystem'
    integer, dimension(:), pointer :: dindx,mindx
!
    nullify(dindx,mindx)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the kernel system for which the matrix is to be filled here, is not initiated yet. "//&
            "call initiateSerialKernelLinearSystem first using a data and model subspace (or the same dmspace)",myname)
       return
    end if
    if(.not.that%initiated) then
       call add(errmsg,2,"the kernel system from which the kernel matrix is to be copied, is not initiated yet. "//&
            "call initiateSerialKernelLinearSystem first and read in a matrix.",myname)
       return
    end if
    if(this%nmval == 0) then
       call add(errmsg,2,"kernel linear system for which the matrix is to be filled here, was initiated without "//&
            "any model values (only to read in data). cannot copy any kernel matrix",myname)
       return
    end if
    if(.not.associated(that%K)) then
       call add(errmsg,2,"it seems that the kernel matrix which is to be copied is not yet read in",myname)
       return
    end if
!
    dindx => getIndxDataSamples(this%dmspace,that%dmspace)
    if(.not.associated(dindx)) then
       call add(errmsg,2,"the data space of the new kernel system is not a subspace of that of the original kernel "//&
            "system",myname)
       goto 1
    end if
    if(size(dindx) /= this%ndata .or. any(dindx == -1)) then
       call add(errmsg,2,"the data space of the new kernel system is not completely contained in the data space "//&
            "of the original kernel system",myname)
       goto 1
    end if
    mindx => getIndxModelValues(this%dmspace,that%dmspace)
    if(.not.associated(mindx)) then
       call add(errmsg,2,"the model space of the new kernel system is not a subspace of that of the original kernel "//&
            "system",myname)
       goto 1
    end if
    if(size(mindx) /= this%nmval .or. any(mindx == -1)) then
       call add(errmsg,2,"the model space of the new kernel system is not completely contained in the model space "//&
            "of the original kernel system",myname)
       goto 1
    end if
!
    this%K(1:this%ndata,1:this%nmval) = that%K(dindx,mindx)
!
1   if(associated(dindx)) deallocate(dindx)
    if(associated(mindx)) deallocate(mindx)
  end subroutine copyMatrixKernelLinearSystem
!------------------------------------------------------------------------
!> \brief read in synthetic data files for system right-hand-side (not parallelized)
!! \details Path-wise construction. Apply event and station filters, and normalize if requested.
!! \param this kernel system
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param nfreq_synthetic_data number of frequencies of synthetic data (must equal length of vector ifreq_synthetic_data)
!! \param ifreq_synthetic_data vector containing all indices of frequencies of synthetic data for current iteration step (used to interpret filter files)
!! \param path_synthetic_data path where synthetic data files are
!! \param lu file unit
!! \param errmsg error message
!! \param apply_event_filter optional logical flag whether event filters are to be accounted for
!! \param path_event_filter optional path were event filter files are (must be present if apply_event_filter == .true., otherwise ignored)
!! \param apply_station_filter optional logical flag whether station filters are to be accounted for
!! \param path_station_filter optional path were station filter files are (must be present if apply_station_filter == .true., otherwise ignored)
!! \param ignore_data_weights optional logical flag whether to ignore the data weights contained in data space
!! \param apply_sdata_normalization optional logical flag whether to normalize the synthetic (correction) values by respective normalization factors contained in data space
!! \param read_synthetic_corrections optional logical flag whether to read in synthetic correction data from "corr_*" files instead of regular synthetic data from "synthetics_*" files
!
  subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
       nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
       apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
       ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
    ! incoming
    type (kernel_linear_system) :: this
    integer :: nfreq_measured_data,nfreq_synthetic_data,lu
    integer, dimension(:) ::  ifreq_measured_data,ifreq_synthetic_data
    character(len=*) :: path_synthetic_data
    character(len=*), optional :: path_event_filter,path_station_filter
    logical, optional :: apply_event_filter,apply_station_filter,ignore_data_weights,apply_sdata_normalization,&
         read_synthetic_corrections
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=41) :: myname='readSyntheticDataSerialKernelLinearSystem'
    character(len=400)  :: errstr
    type (error_message) :: errmsg2,errmsg3
    integer :: j
    logical, dimension(:), allocatable :: idata_added
    ! data space
    integer :: ipath,icomp
    character(len=character_length_evid) :: evid,evid_old
    character(len=character_length_staname) :: staname,staname_old
    character(len=character_length_component), dimension(:), pointer :: comp,pcomp
    character(len=2), dimension(:), pointer :: pimre
    real, dimension(:), pointer :: pwdata
    integer, dimension(:), pointer :: idata_path,idata_path_comp,pifreq
    logical :: ignore_wdata,normalize_sdata,ldummy
    real, dimension(:), pointer :: normalization_factors
    ! filtering
    complex, dimension(:), pointer :: event_filter,station_comp_filter
    logical :: do_event_filtering,do_station_filtering,do_filtering
    integer, dimension(:), allocatable :: map_ifreq_to_synthetic_index,map_ifreq_to_measured_index
    ! synthetic data
    character(len=11) :: sdata_file_base
    logical :: read_corrections
    complex, dimension(:), pointer :: sdata_cmplx
    real, dimension(:), pointer :: s_or_sc_data
!
    nullify(comp,pcomp,pimre,pwdata,idata_path,idata_path_comp,pifreq,normalization_factors,&
         event_filter,station_comp_filter,sdata_cmplx,s_or_sc_data)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before reading in synthetic data",myname)
       return
    end if
!
    if(present(read_synthetic_corrections)) then
       read_corrections = read_synthetic_corrections
    else
       read_corrections = .false.
    end if
   if(read_corrections) then
      call add(errmsg,0,"BY INCOMING FLAG: WILL READ IN SYNTHETIC CORRECTION DATA",myname)
      if(associated(this%scdata)) then
         call add(errmsg,2,"synthetic correction data seems to be read in already. please deallocate linear "//&
              "system first",myname)
         return
      end if
   else ! read_corrections
      if(associated(this%sdata)) then
         call add(errmsg,2,"synthetic data seems to be read in already. please deallocate linear "//&
              "system first",myname)
         return
      end if
   end if ! read_corrections
!
    do_event_filtering = .false.
    if(present(apply_event_filter)) then
       if(apply_event_filter) then
          if(present(path_event_filter)) then
             do_event_filtering = .true.
          else
             call add(errmsg,2,&
                  "event filtering is requested by flag apply_event_filter, but no path of event filters is given ",myname)
             return
          end if
       end if
    end if
    if(do_event_filtering) then
       call add(errmsg,0,"incoming path of event filters : '"//trim(path_event_filter)//"'",myname)
    else
       call add(errmsg,0,"no incoming path of event filters, hence will NOT apply any event filters",myname)
    end if
!
    do_station_filtering = .false.
    if(present(apply_station_filter)) then
       if(apply_station_filter) then
          if(present(path_station_filter)) then
             do_station_filtering = .true.
          else
             call add(errmsg,2,&
                  "station filtering is requested by flag apply_station_filter, but no path of station filters is given ",myname)
             return
          end if
       end if
    end if
    if(do_station_filtering) then
       call add(errmsg,0,"incoming path of station filters: '"//trim(path_station_filter)//"'",myname)
    else
       call add(errmsg,0,"no incoming path of station filters, hence will NOT apply any station filters",myname)
    end if
!
    do_filtering = do_event_filtering .or. do_station_filtering
!
    call add(errmsg,0,"incoming path of synthetic data: '"//trim(path_synthetic_data)//"'",myname)
!
    if(present(ignore_data_weights)) then
       ignore_wdata = ignore_data_weights
    else
       ! by default do not ignore, but apply the data weights!
       ignore_wdata = .false.
    end if
    if(ignore_wdata) then
       call add(errmsg,0,"WILL IGNORE data weights defined in data space",myname)
    else
       call add(errmsg,0,"will apply data weights defined in data space normally (will not ignore them)",myname)
    end if
!
    if(present(apply_sdata_normalization)) then
       normalize_sdata = apply_sdata_normalization
    else
       ! by default do not normalize the sdata values
       normalize_sdata = .false.
    end if
    if(normalize_sdata) then
       call add(errmsg,0,"WILL NORMALIZE synthetic data values as defined in data space",myname)
    else
       call add(errmsg,0,"will NOT normalize synthetic data values",myname)
    end if
!
    if(nfreq_measured_data.le.0) then
       call add(errmsg,2,"number frequencies of the measured data is less or equal to zero",myname)
       return
    end if
    if(nfreq_measured_data /= size(ifreq_measured_data)) then
       write(errstr,*) "number frequencies of the measured data ( = ",nfreq_measured_data,&
            ") does not match the size of the corresponding vector of frequency indices ( = ",&
            size(ifreq_measured_data),")"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(nfreq_synthetic_data.le.0) then
       call add(errmsg,2,"number frequencies of the synthetic data is less or equal to zero",myname)
       return
    end if
    if(nfreq_synthetic_data /= size(ifreq_synthetic_data)) then
       write(errstr,*) "number frequencies of the synthetic data ( = ",nfreq_synthetic_data,&
            ") does not match the size of the corresponding vector of frequency indices ( = ",&
            size(ifreq_synthetic_data),")"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    if(read_corrections) then
       sdata_file_base = 'corr_'
       allocate(this%scdata(this%ndata))
       s_or_sc_data => this%scdata
    else ! read_corrections
       sdata_file_base = 'synthetics_'
       allocate(this%sdata(this%ndata))
       s_or_sc_data => this%sdata
    end if ! read_corrections
!
    ! build the inverse maps to vectors ifreq_measured_data and ifreq_synthetic_data, e.g.:
    !   array map_ifreq_to_synthetic_index maps a frequency index jf to the corresponding index i in array ifreq_synthetic_data, 
    !   i.e. if ifreq_synthetic_data(i) = jf, we define map_ifreq_to_synthetic_index(jf) = i
    ! initiating the maps with value -1 at least causes the program to crash, if at any point of the routine the 
    ! map should be used for frequency indices which are NOT contained in vectors
    ! ifreq_measured_data or ifreq_synthetic_data, respectively.
    allocate(map_ifreq_to_synthetic_index(maxval(ifreq_synthetic_data))); map_ifreq_to_synthetic_index(:) = -1
    do j = 1,nfreq_synthetic_data
       map_ifreq_to_synthetic_index(ifreq_synthetic_data(j)) = j
    end do
    allocate(map_ifreq_to_measured_index(maxval(ifreq_measured_data))); map_ifreq_to_measured_index(:) = -1
    do j = 1,nfreq_measured_data
       map_ifreq_to_measured_index(ifreq_measured_data(j)) = j
    end do
!
    ! get the sdata normalization factors, if requested
    if(normalize_sdata) then
       call getNormalization(this%dmspace,(/ (j, j=1,this%ndata) /),errmsg,norm_sdata=normalization_factors)
       if(.level.errmsg == 2) goto 1
    end if ! normalize_sdata
!
    ! keep track on which rows were added to the kernel matrix (to avoid rows which are not initiated, might still be equal to zero on exit)
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
!
    nullify(event_filter,station_comp_filter,sdata_cmplx)
    evid_old = ''; staname_old = ''
    nullify(pifreq,pcomp,pimre,pwdata)
    nullify(idata_path,comp,idata_path_comp)
!
    ! loop over paths
    ipath = 0
    do while(nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp))
       ipath = ipath + 1
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       write(errstr,*) "reading synthetic data for ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
       call add(errmsg2,0,errstr,myname)
!
       ! read in event_filter for this path, if necessary
       if(do_event_filtering .and. evid_old /= evid) then
          if(associated(event_filter)) deallocate(event_filter)
          errmsg3 = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu,event_filter,ndata=nfreq_measured_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) call print(errmsg3)
          if(.level.errmsg3 == 2) then
             call add(errmsg2,2,"error reading event filter from formatted file '"//trim(path_event_filter)//&
                  "filter_"//trim(evid)//"'",myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") reading event filter"
             call add(errmsg,2,errstr,myname)
             goto 2
          endif
          call dealloc(errmsg3)
          ! remember from which event the filter values are
          evid_old = evid
       end if ! do_event_filtering
!
       do icomp = 1,size(comp)
!
          ! read in synthetic data for this path and component
          if(associated(sdata_cmplx)) deallocate(sdata_cmplx)
          errmsg3 = readAsciiData(trim(path_synthetic_data)//trim(sdata_file_base)//trim(evid)//"_"//trim(staname)//"_"//&
               trim(comp(icomp)),lu,sdata_cmplx,ndata=nfreq_synthetic_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) call print(errmsg3)
          if(.level.errmsg3 == 2) then
             call add(errmsg2,2,"error reading synthetic data from formatted file '"//trim(path_synthetic_data)//&
                  trim(sdata_file_base)//trim(evid)//"_"//trim(staname)//"_"//trim(comp(icomp))//"' for component '"//&
                  trim(comp(icomp))//"'",myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") reading synthetic data for component '"//trim(comp(icomp))//"'"
             call add(errmsg,2,errstr,myname)
             goto 2
          endif
          call dealloc(errmsg3)
!
          ! read in station_comp_filter for this path and component
          if(do_station_filtering) then
             if(associated(station_comp_filter)) deallocate(station_comp_filter)
             errmsg3 = readAsciiData(trim(path_station_filter)//"filter_"//trim(staname)//"_"//trim(comp(icomp)),&
                  lu,station_comp_filter,ndata=nfreq_measured_data)
             call addTrace(errmsg3,myname)
             if(.level.errmsg3 /= 0) call print(errmsg3)
             if(.level.errmsg3 == 2) then
                call add(errmsg2,2,"error reading station filter from formatted file '"//trim(path_station_filter)//&
                     "filter_"//trim(staname)//"_"//trim(comp(icomp))//"'",myname)
                call print(errmsg2)
                write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                     ") reading station filter for component '"//trim(comp(icomp))//"'"
                call add(errmsg,2,errstr,myname)
                goto 2
             endif
             call dealloc(errmsg3)
          end if ! do_station_filtering
!
          ! get info on data samples for this component and current path
          if(associated(pcomp)) deallocate(pcomp)
          if(associated(pifreq)) deallocate(pifreq); nullify(pifreq)
          if(associated(pimre)) deallocate(pimre); nullify(pimre)
          if(associated(pwdata)) deallocate(pwdata); nullify(pwdata)
          if(associated(idata_path_comp)) deallocate(idata_path_comp)
          allocate(pcomp(1)); pcomp(1) = comp(icomp)
          idata_path_comp => getIndxDataSamples(this%dmspace,ifreq=pifreq,comp=pcomp,imre=pimre,wdata=pwdata,indx_in=idata_path)
          ! idata_path_comp should be associated here (since comp(icomp) is in data space for this path), don't check. 
          ! otherwise function nextPathDataModelSpaceInfo is corrupt
          write(errstr,*) "there are ",size(idata_path_comp),&
               " data samples in data space for this path and component '"//trim(comp(icomp))//"'"
          call add(errmsg2,0,errstr,myname)
!
          ! check, if for any reasons the current values of synthetic data were already defined
          if(any(idata_added(idata_path_comp))) then
             write(errstr,*) "there are ",count(idata_added(idata_path_comp))," data sample values which are ",&
                  "already defined. data_model_space_info object may be corrupt"
             call add(errmsg2,2,errstr,myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
             call add(errmsg,2,errstr,myname)
             goto 2
          end if
!
          ! neutralize the data weights, if they are to be ignored
          if(ignore_wdata) pwdata = 1.0
          ! in case of sdata normalization, multiply normalization factors into the weights pwdata in order to account for normalization
          if(normalize_sdata) pwdata = pwdata*normalization_factors(idata_path_comp)
!
          ! filter the synthetic data values, overwriting the values of sdata (it's a local copy)
          ! compute the filtered values at ALL frequencies, even if some might not be needed for this path and comp.
          ! Alternatively, you would need the total different ifreq values for this path and comp. Please note that 
          ! pifreq contains every frequency index twice (in a sensible setup): one for the 'im' and one for the 're' 
          ! data sample
          if(do_event_filtering) sdata_cmplx = sdata_cmplx * &
               event_filter(map_ifreq_to_measured_index(ifreq_synthetic_data))
          if(do_station_filtering) sdata_cmplx = sdata_cmplx * &
               station_comp_filter(map_ifreq_to_measured_index(ifreq_synthetic_data))
!
          ! finally set values in this%sdata
          where(pimre == 'im')
             s_or_sc_data(idata_path_comp) = pwdata * imag(sdata_cmplx(map_ifreq_to_synthetic_index(pifreq)))
             idata_added(idata_path_comp) = .true.
          end where
          where(pimre == 're')
             s_or_sc_data(idata_path_comp) = pwdata * real(sdata_cmplx(map_ifreq_to_synthetic_index(pifreq)))
             idata_added(idata_path_comp) = .true.
          end where
!
       end do ! icomp
!
    end do ! ipath
!
    ! check if there are any rows which were not filled
    if(.not.all(idata_added)) then
       write(errstr,*) "maybe some data samples occur more than once in the data space (not allowed!), as ",&
            count(.not.idata_added)," synthetic data values were not set"
            !", namely these ones (by index) :",pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(pimre)) deallocate(pimre)
    if(associated(pwdata)) deallocate(pwdata)
    if(associated(idata_path_comp)) deallocate(idata_path_comp)
    if(associated(normalization_factors)) deallocate(normalization_factors)
!
    if(associated(event_filter)) deallocate(event_filter)
    if(associated(station_comp_filter)) deallocate(station_comp_filter)
    if(allocated(idata_added)) deallocate(idata_added)
!
    if(associated(sdata_cmplx)) deallocate(sdata_cmplx)
!
    if(allocated(map_ifreq_to_synthetic_index)) deallocate(map_ifreq_to_synthetic_index)
    if(allocated(map_ifreq_to_measured_index)) deallocate(map_ifreq_to_measured_index)
!
    ! if code comes here, return normally
    return
!
    ! if there is an error inside the loop on nextPath, need to reset the iterator (and deallocate all the pointers, etc.)
2   ldummy = nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp,reset=.true.)
    goto 1
  end subroutine readSyntheticDataSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief read in measured data files for system right-hand-side (not parallelized)
!! \details Path-wise construction. Normalize by measured data normalization factors, if requested.
!! \param this kernel system
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param path_measured_data path where measured data files are
!! \param lu file unit
!! \param errmsg error message
!! \param ignore_data_weights optional logical flag whether to ignore the data weights contained in data space
!! \param apply_mdata_normalization optional logical flag whether to normalize the measured data values by respective normalization factors contained in data space
!
  subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
       path_measured_data,lu,errmsg,ignore_data_weights,apply_mdata_normalization)
    ! incoming
    type (kernel_linear_system) :: this
    integer, dimension(:) ::  ifreq_measured_data
    integer :: nfreq_measured_data,lu
    character(len=*) :: path_measured_data
    logical, optional :: ignore_data_weights,apply_mdata_normalization
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=40) :: myname='readMeasuredDataSerialKernelLinearSystem'
    character(len=400)  :: errstr
    type (error_message) :: errmsg2,errmsg3
    integer :: j
    logical, dimension(:), allocatable :: idata_added
    ! data space
    integer :: ipath,icomp
    character(len=character_length_evid) :: evid,evid_old
    character(len=character_length_staname) :: staname,staname_old
    character(len=character_length_component), dimension(:), pointer :: comp,pcomp
    character(len=2), dimension(:), pointer :: pimre
    real, dimension(:), pointer :: pwdata
    integer, dimension(:), pointer :: pifreq,idata_path,idata_path_comp
    logical :: ignore_wdata,normalize_mdata,ldummy
    real, dimension(:), pointer :: normalization_factors
    ! measured data
    complex, dimension(:), pointer :: mdata_cmplx
    integer, dimension(:), allocatable :: map_ifreq_to_measured_index
!
    nullify(comp,pcomp,pimre,pwdata,pifreq,idata_path,idata_path_comp,normalization_factors,mdata_cmplx)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before reading in synthetic data",myname)
       return
    end if
    if(associated(this%mdata)) then
       call add(errmsg,2,"measured data vector seems to be read in already. please deallocate linear system first",myname)
       return
    end if
!
    call add(errmsg,0,"incoming path of measured data: '"//trim(path_measured_data)//"'",myname)
!
    if(present(ignore_data_weights)) then
       ignore_wdata = ignore_data_weights
    else
       ! by default do not ignore, but apply the data weights!
       ignore_wdata = .false.
    end if
    if(ignore_wdata) then
       call add(errmsg,0,"WILL IGNORE data weights defined in data space",myname)
    else
       call add(errmsg,0,"will apply data weights defined in data space normally (will not ignore them)",myname)
    end if
!
    if(present(apply_mdata_normalization)) then
       normalize_mdata = apply_mdata_normalization
    else
       ! by default do not normalize the scdata values
       normalize_mdata = .false.
    end if
    if(normalize_mdata) then
       call add(errmsg,0,"WILL NORMALIZE measured data values as defined in data space",myname)
    else
       call add(errmsg,0,"will NOT normalize measured data values",myname)
    end if
!
    if(nfreq_measured_data.le.0) then
       call add(errmsg,2,"number frequencies of the measured data is less or equal to zero",myname)
       return
    end if
    if(nfreq_measured_data /= size(ifreq_measured_data)) then
       write(errstr,*) "number frequencies of the measured data ( = ",nfreq_measured_data,&
            ") does not match the size of the corresponding vector of frequency indices ( = ",&
            size(ifreq_measured_data),")"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    allocate(this%mdata(this%ndata))
!
    ! build the inverse maps to vector ifreq_measured_data, e.g.:
    !   array map_ifreq_to_measured_index maps a frequency index jf to the corresponding index i in array ifreq_measured_data, 
    !   i.e. if ifreq_measured_data(i) = jf, we define map_ifreq_to_measured_index(jf) = i
    ! initiating the map with value -1 at least causes the program to crash, if at any point of the routine the 
    ! map should be used for frequency indices which are NOT contained in vector ifreq_measured_data.
    allocate(map_ifreq_to_measured_index(maxval(ifreq_measured_data))); map_ifreq_to_measured_index(:) = -1
    do j = 1,nfreq_measured_data
       map_ifreq_to_measured_index(ifreq_measured_data(j)) = j
    end do
!
    ! get the mdata normalization factors, if requested
    if(normalize_mdata) then
       call getNormalization(this%dmspace,(/ (j, j=1,this%ndata) /),errmsg,norm_mdata=normalization_factors)
       if(.level.errmsg == 2) goto 1
    end if ! normalize_mdata
!
    ! keep track on which rows were added to the kernel matrix (to avoid rows which are not initiated, might still be equal to zero on exit)
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
!
    evid_old = ''; staname_old = ''
    nullify(pifreq,pcomp,pimre,pwdata)
    nullify(idata_path_comp)
!
    ! loop over paths
    ipath = 0
    do while(nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp))
       ipath = ipath + 1
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       write(errstr,*) "reading measured data for ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
       call add(errmsg2,0,errstr,myname)
!
       do icomp = 1,size(comp)
!
          ! read in measured data for this path and component
          if(associated(mdata_cmplx)) deallocate(mdata_cmplx)
          errmsg3 = readAsciiData(trim(path_measured_data)//"data_"//trim(evid)//"_"//trim(staname)//"_"//trim(comp(icomp)),&
                  lu,mdata_cmplx,ndata=nfreq_measured_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) call print(errmsg3)
          if(.level.errmsg3 == 2) then
             call add(errmsg2,2,"error reading measured data from formatted file '"//trim(path_measured_data)//&
                  "data_"//trim(evid)//"_"//trim(staname)//trim(comp(icomp))//"' for component'"//&
                  trim(comp(icomp))//"'",myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//&
                  ") reading measured data for component '"//trim(comp(icomp))//"'"
             call add(errmsg,2,errstr,myname)
             goto 2
          end if
          call dealloc(errmsg3)
!
          ! get info on data samples for this component and current path
          if(associated(pcomp)) deallocate(pcomp)
          if(associated(pifreq)) deallocate(pifreq); nullify(pifreq)
          if(associated(pimre)) deallocate(pimre); nullify(pimre)
          if(associated(pwdata)) deallocate(pwdata); nullify(pwdata)
          if(associated(idata_path_comp)) deallocate(idata_path_comp)
          allocate(pcomp(1)); pcomp(1) = comp(icomp)
          idata_path_comp => getIndxDataSamples(this%dmspace,ifreq=pifreq,comp=pcomp,imre=pimre,wdata=pwdata,indx_in=idata_path)
          ! idata_path_comp should be associated here (since comp(icomp) is in data space for this path), don't check. 
          ! otherwise function nextPathDataModelSpaceInfo is corrupt
          write(errstr,*) "there are ",size(idata_path_comp),&
               " data samples in data space for this path and component '"//trim(comp(icomp))//"'"
          call add(errmsg2,0,errstr,myname)
!
          ! check, if for any reasons the current values of measured data were already defined
          if(any(idata_added(idata_path_comp))) then
             write(errstr,*) "there are ",count(idata_added(idata_path_comp))," measured data sample values which are ",&
                  "already defined. data_model_space_info object may be corrupt"
             call add(errmsg2,2,errstr,myname)
             call print(errmsg2)
             write(errstr,*) "error in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
             call add(errmsg,2,errstr,myname)
             goto 2
          end if
!
          ! neutralize the data weights, if they are to be ignored
          if(ignore_wdata) pwdata = 1.0
          ! in case of mdata normalization, multiply normalization factors into the weights pwdata in order to account for normalization
          if(normalize_mdata) pwdata = pwdata*normalization_factors(idata_path_comp)
!
          ! finally set values in this%mdata
          where(pimre == 'im')
             this%mdata(idata_path_comp) = pwdata * imag(mdata_cmplx(map_ifreq_to_measured_index(pifreq)))
             idata_added(idata_path_comp) = .true.
          end where
          where(pimre == 're')
             this%mdata(idata_path_comp) = pwdata * real(mdata_cmplx(map_ifreq_to_measured_index(pifreq)))
             idata_added(idata_path_comp) = .true.
          end where
!
       end do ! icomp
!
    end do ! ipath
!
    ! check if there are any rows which were not filled
    if(.not.all(idata_added)) then
       write(errstr,*) "maybe some data samples occur more than once in the data space (not allowed!), as ",&
            count(.not.idata_added)," measured data values were not set"
            !", namely these ones (by index) :",pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(pimre)) deallocate(pimre)
    if(associated(pwdata)) deallocate(pwdata)
    if(associated(idata_path_comp)) deallocate(idata_path_comp)
    if(associated(normalization_factors)) deallocate(normalization_factors)
!
    if(allocated(idata_added)) deallocate(idata_added)
    if(associated(mdata_cmplx)) deallocate(mdata_cmplx)
    if(allocated(map_ifreq_to_measured_index)) deallocate(map_ifreq_to_measured_index)
!
    ! if code comes here, return normally
    return
!
    ! if there is an error inside the loop on nextPath, need to reset the iterator (and deallocate all the pointers, etc.)
2   ldummy = nextPathDataModelSpaceInfo(this%dmspace,evid,staname,indx=idata_path,all_comp=comp,reset=.true.)
    goto 1
  end subroutine readMeasuredDataSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief define one right hand side vector of the linear system by data residual mdata-sdata
!! \param this kernel linear system
!! \param errmsg error message
  subroutine setRhsAsDataResidualKernelLinearSystem(this,errmsg)
    type (kernel_linear_system) :: this
    type (error_message) :: errmsg
    ! local
    character(len=38) :: myname = 'setRhsAsDataResidualKernelLinearSystem'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first "//&
            "and read in synthetic and measured data",myname)
       return
    end if
    if(.not.(associated(this%mdata).and.associated(this%sdata))) then
       call add(errmsg,2,"synthetic data or measured data not yet read in. required for computing residuals",myname)
       return
    end if
    if(associated(this%rhs)) then
       call add(errmsg,1,"right hand side(s) were already allocated, deallocating now before setting to new ones. "//&
            "also deallocate any existing solution vector(s) of the system",myname)
       deallocate(this%rhs)
       if(associated(this%sol)) deallocate(this%sol)
       this%transpose = ''
    end if
!
    allocate(this%rhs(this%nrow,1)); this%nrhs = 1
    this%rhs(:,1) = 0.
!
    this%rhs(1:this%ndata,1) = this%mdata - this%sdata
!
    this%transpose = 'N'
  end subroutine setRhsAsDataResidualKernelLinearSystem
!------------------------------------------------------------------------
!> \brief define one right hand side vector of the linear system by corrected data residual mdata-sdata-scdata
!! \param this kernel linear system
!! \param errmsg error message
  subroutine setRhsAsCorrectedDataResidualKernelLinearSystem(this,errmsg)
    type (kernel_linear_system) :: this
    type (error_message) :: errmsg
    ! local
    character(len=47) :: myname = 'setRhsAsCorrectedDataResidualKernelLinearSystem'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first "//&
            "and read in synthetic (correction) and measured data",myname)
       return
    end if
    if(.not.(associated(this%mdata).and.associated(this%sdata).and.associated(this%scdata))) then
       call add(errmsg,2,"measured data or synthetic data or synthetics corrections not yet read in. "//&
            "required for computing corrected residuals",myname)
       return
    end if
    if(associated(this%rhs)) then
       call add(errmsg,1,"right hand side(s) were already allocated, deallocating now before setting to new ones",myname)
       deallocate(this%rhs)
       if(associated(this%sol)) deallocate(this%sol)
       this%transpose = ''
    end if
!
    allocate(this%rhs(this%nrow,1)); this%nrhs = 1
    this%rhs(:,1) = 0.
!
    this%rhs(1:this%ndata,1) = this%mdata - this%sdata - this%scdata
!
    this%transpose = 'N'
  end subroutine setRhsAsCorrectedDataResidualKernelLinearSystem
!------------------------------------------------------------------------
!> \brief define right hand side vector(s) of the linear system manually
!! \param this kernel linear system
!! \param errmsg error message
  subroutine setRhsKernelLinearSystem(this,transpose,rhs,errmsg)
    type (kernel_linear_system) :: this
    character(len=*) :: transpose
    real, dimension(:,:) :: rhs
    type (error_message) :: errmsg
    ! local
    character(len=24) :: myname = 'setRhsKernelLinearSystem'
    character(len=400) :: errstr
    integer :: nrhs
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before setting right hand side(s)",myname)
       return
    end if
    if(associated(this%rhs)) then
       call add(errmsg,1,"right hand side(s) were already allocated, deallocating now before setting to new ones. "//&
            "also deallocating any existing solution(s) of this system",myname)
       deallocate(this%rhs)
       if(associated(this%sol)) deallocate(this%sol)
       this%transpose = ''
       this%nrhs = 0
    end if
!
    if(size(rhs,1) <= 0) then
       call add(errmsg,2,"length of incoming right-hand-side vectors is (less than) zero",myname)
       return
    end if
    nrhs = size(rhs,2)
    if(nrhs <= 0) then
       call add(errmsg,2,"number of incoming right-hand-side vectors is (less than) zero",myname)
       return
    end if
!
    select case(transpose)
    case ('N')
       if(size(rhs,1) /= this%nrow) then
          write(errstr,*) "incoming length of right-hand-side-vectors is ",size(rhs,1),&
               ", but must equal the number of rows ",this%nrow," of the system (not transposed)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       this%transpose = 'N'
    case ('T')
       if(size(rhs,1) /= this%ncol) then
          write(errstr,*) "incoming length of right-hand-side-vectors is ",size(rhs,1),&
               ", but must equal the number of columns ",this%ncol," of the system (transposed)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       this%transpose = 'T'
    case default
       call add(errmsg,2,"Incoming character '"//trim(transpose)//"' does not indicate transpose state. "//&
            "Must be either 'N' or 'T'.",myname)
       return
    end select
    allocate(this%rhs(size(rhs,1),nrhs))
    this%rhs = rhs
    this%nrhs = nrhs
  end subroutine setRhsKernelLinearSystem
!------------------------------------------------------------------------
!> \brief define solution vector(s) of the linear system manually (for forward multiplication)
!! \param this kernel linear system
!! \param errmsg error message
  subroutine setSolKernelLinearSystem(this,transpose,sol,errmsg)
    type (kernel_linear_system) :: this
    character(len=*) :: transpose
    real, dimension(:,:) :: sol
    type (error_message) :: errmsg
    ! local
    character(len=24) :: myname = 'setSolKernelLinearSystem'
    character(len=400) :: errstr
    integer :: nrhs
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before setting solution(s)",myname)
       return
    end if
    if(associated(this%sol)) then
       call add(errmsg,1,"solution(s) were already allocated, deallocating now before setting to new one(s). "//&
            "also deallocating any existing right-hand-side(s) of this system",myname)
       deallocate(this%sol)
       if(associated(this%rhs)) deallocate(this%rhs)
       this%transpose = ''
       this%nrhs = 0
    end if
!
    if(size(sol,1) <= 0) then
       call add(errmsg,2,"length of incoming solution vector(s) is (less than) zero",myname)
       return
    end if
    nrhs = size(sol,2)
    if(nrhs <= 0) then
       call add(errmsg,2,"number of incoming solution vectors is (less than) zero",myname)
       return
    end if
!
    select case(transpose)
    case ('N')
       if(size(sol,1) /= this%ncol) then
          write(errstr,*) "incoming length of solution vector(s) is ",size(sol,1),&
               ", but must equal the number of columns ",this%ncol," of the system (not transposed)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       this%transpose = 'N'
    case ('T')
       if(size(sol,1) /= this%nrow) then
          write(errstr,*) "incoming length of solution vector(s) is ",size(sol,1),&
               ", but must equal the number of rows ",this%nrow," of the system (transposed)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       this%transpose = 'T'
    case default
       call add(errmsg,2,"Incoming character '"//trim(transpose)//"' does not indicate transpose state. "//&
            "Must be either 'N' or 'T'.",myname)
       return
    end select
    allocate(this%sol(size(sol,1),nrhs))
    this%sol = sol
    this%nrhs = nrhs
  end subroutine setSolKernelLinearSystem
!------------------------------------------------------------------------
!> \brief add ONE regularization column (the n'th one) to the system
!! \param 
!
  subroutine setColumnRegKernelLinearSystem(this,v,n,errmsg)
    type (kernel_linear_system) :: this
    real, dimension(:) :: v
    integer :: n
    type (error_message) :: errmsg
    ! local
    character(len=30) :: myname = 'addColumnRegKernelLinearSystem'
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "before setting any column regularization",myname)
       return
    end if
    if(n.gt.this%ncolreg .or. n.lt.1) then
       write(errstr,*) "requested position of incoming regularizing column ",n,&
            " is invalid: must be between 1 and the number of allocated regularizing columns = ",this%ncolreg
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(v) /= this%nrow) then
       write(errstr,*) "length of incoming regularizing column ",size(v)," does not match number of rows ",&
            this%nrow," of the kernel system"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%K(:,this%nmval+n) = v
  end subroutine setColumnRegKernelLinearSystem
!------------------------------------------------------------------------
!> \brief solve kernel linear system (not parallelized)
!! \details Use the matrix K and the residual vector as the right hand side to form a 
!!  linear system. Compute the solution of that system, which is set as this%sol
!! \param this kernel system
!! \param errmsg error message
!! \return error message
!
  subroutine solveSerialKernelLinearSystem(this,errmsg,overwrite_A)
    ! incoming
    type (kernel_linear_system) :: this
    logical, optional :: overwrite_A
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=29) :: myname='solveSerialKernelLinearSystem'
    type (error_message) :: errmsg2
    ! linear system
    type (serial_linear_system) :: LSE
    real, dimension(:,:), pointer :: b,solution,residual
!
    nullify(b,solution,residual)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first "//&
            "and read in kernel matrix and set any right hand sides before system can be solved",myname)
       return
    end if
    if(.not.associated(this%rhs)) then
       call add(errmsg,2,"it seems, that right-hand-side vector(s) were not yet set",myname)
       return
    endif
    if(associated(this%sol)) then
       call add(errmsg,2,"solution vector(s) of the system are already defined. this should not happen when "//&
            "properly solving the kernel system",myname)
       return
    end if
!
    ! make a copy of rhs, not to overwrite it
    allocate(b(size(this%rhs,1),size(this%rhs,2))); b = this%rhs
    nullify(solution,residual)
!
    call createSerialLinearSystem(LSE,this%nrow,this%ncol,1,this%K,b,transpose=this%transpose)
    errmsg2 = solveLeastSquaresSerialLinearSystem(LSE,solution,residual,overwrite_A)
    call addTrace(errmsg2,myname)
    call print(errmsg2)
    if(.level.errmsg2/=0) then
       call add(errmsg,.level.errmsg2,"error or warning in routine 'solveLeastSquaresLinearSystem'",myname)
    endif
    if(.level.errmsg2==2) then
       deallocate(b)
       if(associated(solution)) deallocate(solution)
       if(associated(residual)) deallocate(residual)
       call dealloc(errmsg2)
       call dealloc(LSE)
       return
    endif
    call dealloc(errmsg2)
!
    this%sol => solution
    nullify(solution)
!
    deallocate(b)
    call dealloc(LSE)
  end subroutine solveSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief multiply kernel matrix with given solution vector(s) and yield rhs(s)
!! \details can be used for computing approximate data for some given model
!! \param this kernel system
!! \param errmsg error message
!! \return error message
!
  subroutine multiplyForwardKernelLinearSystem(this,errmsg)
    ! incoming
    type (kernel_linear_system) :: this
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=33) :: myname='multiplyForwardKernelLinearSystem'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"the linear system is not yet initiated. call initiateSerialKernelLinearSystem first, "//&
            "read in kernel matrix and define some solution vector(s), before doing a forward multiplication",myname)
       return
    end if
    if(.not.associated(this%sol)) then
       call add(errmsg,2,"solution vector(s) not yet defined",myname)
       return
    endif
    if(associated(this%rhs)) then
       call add(errmsg,2,"right-hand-side vector(s) of the system are already defined. this should not happen "//&
            "when properly doing a forward multiplication",myname)
       return
    end if
    ! when setting this%sol, then this%rhs is defined and, hence, here it can be trusted
    ! that when this%sol is defined, then this%sol has size (this%ncol,this%rhs) if 'N' or (this%nrow,this%rhs) if 'T'
!
    select case(this%transpose)
    case ('N'); allocate(this%rhs(this%nrow,this%nrhs)) 
    case ('T'); allocate(this%rhs(this%ncol,this%nrhs)) 
    end select
!
    ! BLAS level 3: general matrix-matrix multiplication, also with one solution, (i.e. vector)
    ! this%rhs = this%K * this%solution (or this%rhs = this%K^T * this%solution, in case this%transpose = 'T')
    ! IN THE LONG RUN: MOVE THE BLAS CODE TO MODULE serialLinearSystem ? WOULD  NEED TO CREATE A (linked) OBJECT
    ! OF serial_linear_system HERE FOR THE PURPOSE OF MULTIPLICATION
    call SGEMM(this%transpose , 'N' , size(this%rhs,1) , size(this%rhs,2) , size(this%sol,1) , 1.0 , this%K , &
         size(this%K,1) , this%sol , size(this%sol,1) , 0.0 , this%rhs , size(this%rhs,1))
!
  end subroutine multiplyForwardKernelLinearSystem
!------------------------------------------------------------------------
!> \brief deallocate kernel linear system
!! \param this kernel linear system to be deallocated
!
  subroutine deallocateKernelLinearSystem(this)
    type (kernel_linear_system) :: this
    if (associated(this%K)) deallocate(this%K)
    if (associated(this%rhs)) deallocate(this%rhs)
    if (associated(this%sol)) deallocate(this%sol)
    this%nrow = 0; this%ndata = 0; this%nrowreg = 0 
    this%ncol = 0; this%nmval = 0; this%ncolreg = 0
    this%nrhs = 0
    this%transpose = ''
    if (associated(this%mdata)) deallocate(this%mdata)
    if (associated(this%sdata)) deallocate(this%sdata)
    if (associated(this%scdata)) deallocate(this%scdata)
    call dealloc(this%dmspace)
    this%initiated = .false.
  end subroutine deallocateKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to kernel matrix
!! \param this kernel linear system
!! \param K pointer to kernel matrix of this kernel linear system
!
  function getMatrixKernelLinearSystem(this) result(K)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:,:), pointer :: K
    K => this%K
  end function getMatrixKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to measured data vector
!! \param this kernel linear system
!! \param mdata pointer to measured data vector of this kernel linear system
!
  function getMeasuredDataKernelLinearSystem(this) result(mdata)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:), pointer :: mdata
    mdata => this%mdata
  end function getMeasuredDataKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to synthetic data vector
!! \param this kernel linear system
!! \param sdata pointer to synthetic data vector of this kernel linear system
!
  function getSyntheticDataKernelLinearSystem(this) result(sdata)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:), pointer :: sdata
    sdata => this%sdata
  end function getSyntheticDataKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to synthetic corrections data vector
!! \param this kernel linear system
!! \param sdata pointer to synthetic data vector of this kernel linear system
!
  function getSyntheticCorrectionDataKernelLinearSystem(this) result(scdata)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:), pointer :: scdata
    scdata => this%scdata
  end function getSyntheticCorrectionDataKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to right-hand-sides
!! \param this kernel linear system
!! \param rhs pointer to right-hand-side vector(s) of this kernel linear system
!
  function getRhsKernelLinearSystem(this) result(rhs)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:,:), pointer :: rhs
    rhs => this%rhs
  end function getRhsKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return sum of squares of components of residual vector
!! \param this kernel linear system
!! \param iostat optional integer which tells if result can be trusted (ios=0) or if there was an error (ios/=0), e.g. no residual associated
!! \param misfit sum of squares of components of residual vector
!
  function getMisfitKernelLinearSystem(this,iostat) result(misfit)
    type (kernel_linear_system), intent(in) :: this
    real :: misfit
    integer, optional, intent(out) :: iostat
    if(.not.(associated(this%sdata).and.associated(this%mdata))) then
       misfit = 0.
       if(present(iostat)) iostat = 1 ! error, do not use result
       return
    endif
    misfit = sum((this%mdata-this%sdata)**2) ! return the sum of squares of the differences data - synthetics
    if(present(iostat)) iostat = 0 ! ok, you can trust result
  end function getMisfitKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return sum of squares of components of residual vector
!! \param this kernel linear system
!! \param iostat optional integer which tells if result can be trusted (ios=0) or if there was an error (ios/=0), e.g. no residual associated
!! \param misfit sum of squares of components of residual vector
!
  function getCorrectedMisfitKernelLinearSystem(this,iostat) result(misfit)
    type (kernel_linear_system), intent(in) :: this
    real :: misfit
    integer, optional, intent(out) :: iostat
    if(.not.(associated(this%mdata).and.associated(this%sdata).and.associated(this%scdata))) then
       misfit = 0.
       if(present(iostat)) iostat = 1 ! error, do not use result
       return
    endif
    misfit = sum((this%mdata-this%sdata-this%scdata)**2) ! return the sum of squares of the corrected residual, i.e. the differences data - synthetics - synthetic_correction
    if(present(iostat)) iostat = 0 ! ok, you can trust result
  end function getCorrectedMisfitKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return pointer to solution of kernel linear system (i.e. model uptdate)
!! \param this kernel linear system
!! \param sol pointer to solution of kernel linear system (i.e. model uptdate)
!
  function getSolutionKernelLinearSystem(this) result(sol)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:,:), pointer :: sol
    sol => this%sol
  end function getSolutionKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return number of rows nrow of this kernel linear system
  function getNrowKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%nrow
  end function getNrowKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return number of model values nmval of this kernel linear system
  function getNmvalKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%nmval
  end function getNmvalKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return number of data samples ndata (first rows) of this kernel linear system
  function getNdataKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%ndata
  end function getNdataKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return number of columns ncol of this kernel linear system
  function getNcolKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%ncol
  end function getNcolKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return number of row regularization conditions nrowreg which were added to this kernel linear system
  function getNrowregKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%nrowreg
  end function getNrowregKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return the sum of the column vectors of this%K 
  function getSumColumnsKernelLinearSystem(this) result(v)
    type (kernel_linear_system), intent(in) :: this
    real, dimension(:), pointer :: v
    integer :: j
    nullify(v)
    if(.not.associated(this%K)) return
    allocate(v(this%nrow))
    v(:) = 0.
    do j = 1,this%ncol
       v = v + this%K(1:this%nrow,j)
    end do ! j
  end function getSumColumnsKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return the data and model space info object
  function getDmspaceKernelLinearSystem(this) result(dmspace)
    type (kernel_linear_system), intent(in) :: this
    type (data_model_space_info) :: dmspace
    dmspace = this%dmspace
  end function getDmspaceKernelLinearSystem
!------------------------------------------------------------------------
!> \brief return logical flag, whether this object is initiated or not (i.e. contains sensible values of ndata,nmval,nrowreg,ncolreg,nrow,ncol,dmspace
  function isInitiatedKernelLinearSystem(this) result(initiated)
    type (kernel_linear_system), intent(in) :: this
    logical :: initiated
    initiated = this%initiated
  end function isInitiatedKernelLinearSystem
!!$!------------------------------------------------------------------------
!!$!> \brief multiply incoming matrix with complete (regularized) kernel system matrix
!!$!! \details if transpose is present and 'T', then this%K^T*A = solution will be computed, 
!!$!!  otherwise if transpose is 'N' or not present this%K*A = solution will be computed.
!!$!!  All matrix dimensions must be valid.
!!$!! \param this kernel linear system
!!$!! \param A matrix by which this%K should be multiplied
!!$!! \param solution solution array
!!$!! \param errmsg error message
!!$!! \param transpose optional character to indicate whether K^T*A should be computed ('T') or K*A ('N')
!!$!
!!$  subroutine multiplyMatrixKernelLinearSystem(this,A,solution,errmsg,transpose)
!!$    type (kernel_linear_system) :: this
!!$    real, dimension(:,:) :: A,solution
!!$    character(len=*), optional :: transpose
!!$    ! local
!!$    character(len=32) :: myname = 'multiplyMatrixKernelLinearSystem'
!!$    logical :: compute_transpose
!!$    integer :: LD
!!$!
!!$    call addTrace(errmsg,myname)
!!$!
!!$    compute_transpose = .false.
!!$    if(present(transpose)) then
!!$       select case(transpose)
!!$       case ('T')
!!$          compute_transpose = .true.
!!$       case('N')
!!$          ! OK do nothing
!!$       case default
!!$          call add(errmsg,2,"incoming transpose flag '"//trim(transpose)//"' invalid: must be either 'N' or 'T'",myname)
!!$       end select
!!$    end if
!!$!
!!$    if(compute_transpose) then
!!$       
!!$    else
!!$    end if
!!$  end subroutine multiplyMatrixKernelLinearSystem
!
end module kernelLinearSystem
