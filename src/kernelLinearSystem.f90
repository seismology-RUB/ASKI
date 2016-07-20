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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! TODO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! need ONE routine whicH is called from outside which
!!  - sets up the complete linear system, requireing as input:
!!    * system is parallel or not
!!    * data model space
!!    * spacial regularization parameters (if any)
!!    * type of right hand side (which misfit?, only data minus synt?, additional terms on r.h.s.?)
!!    * 
!!
!! internally, there should be several routines doing specific tasks, like
!!    * setting up kernel matrix
!!    * adding spacial regularization equations (?)
!!    * reading in data
!!    * reading in synthetics
!!    * 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!> \brief set up and handle kernel matrix system
!!
!! \details Given a data_sample_info object and general prequesits
!!  of the inversion and iteration step, the kernel matrix is set up
!!  and solved using module linearSystem. The result is interpreted
!!  and an overall weighted sensitivity may be computed. 
!!
!! \author Florian Schumacher
!! \date May 2013
!
module kernelLinearSystem
!
  use linearSystem
  use dataModelSpaceInfo
  use modelParametrization
  use seismicEvent
  use seismicStation
  use asciiDataIO
  use componentTransformation
  use spectralWaveformKernel
  use vectorPointer
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateKernelLinearSystem; end interface
  interface operator (.KM.); module procedure getMatrixKernelLinearSystem; end interface
  interface operator (.sol.); module procedure getSolutionKernelLinearSystem; end interface
  interface operator (.rhs.); module procedure getRhsKernelLinearSystem; end interface
  interface operator (.nrow.); module procedure getNrowKernelLinearSystem; end interface
  interface operator (.ndata.); module procedure getNdataKernelLinearSystem; end interface
  interface operator (.ncol.); module procedure getNcolKernelLinearSystem; end interface
  interface operator (.nsmooth.); module procedure getNsmoothKernelLinearSystem; end interface
  interface operator (.md.); module procedure getMeasuredDataKernelLinearSystem; end interface
  interface operator (.sd.); module procedure getSyntheticDataKernelLinearSystem; end interface
!
!> \brief kernel linear system object
  type kernel_linear_system
     private
     real, dimension(:,:), pointer :: K => null() !< kernel system matrix, including smoothing conditions (if any)
     real, dimension(:,:), pointer :: rhs => null() !< right-hand-side vector(s) of kernel system
     real, dimension(:,:), pointer :: sol => null() !< solution vector array (model updates) (has size(ncol,nrhs))
     integer :: nrow = 0 !< (local) number of rows of kernel linear system
     integer :: ndata = 0 !< (local) number of rows reserved for kernels or synthetic data (first rows of the system)
     integer :: nsmooth = 0 !< (local) number rows for smoothing conditions (rows) that were added to the system (last rows of the system)
     integer :: ncol = 0 !< (local) number of columns = unknowns of kernel linear system
     integer :: nrhs = 0 !< number of right-hand-side vectors (columns of array rhs)
     real, dimension(:), pointer :: mdata => null() !< vector of measured data according to data space
     real, dimension(:), pointer :: sdata => null() !< vector of synthetic data according to data space
  end type kernel_linear_system
!
contains
!------------------------------------------------------------------------
!> \brief allocate system matrix, according to number of data, parameters and smoothing conditions
!! \param this kernel linear system object
!! \param ndata number of data in data space, which will be used to fill kernel matrix
!! \param nparam number of parameters in model space, which will be used to fill kernel matrix
!! \param nsmooth number of smoothing conditions that will be added to the system
!! \param errmsg error message
!
  subroutine allocateMatrixSerialKernelLinearSystem(this,ndata,nparam,nsmooth,errmsg)
    type (kernel_linear_system) :: this
    integer :: ndata,nparam,nsmooth
    type (error_message) :: errmsg
    ! local
    integer :: nrow,ncol,status
    character(len=400) :: errstr
    character(len=38) :: myname = 'allocateMatrixSerialKernelLinearSystem'
!
    call addTrace(errmsg,myname)
    if(associated(this%K).or.associated(this%sol)) then
       if(associated(this%K)) deallocate(this%K)
       if(associated(this%sol)) deallocate(this%sol)
       this%ncol = 0
       this%nsmooth = 0
       if(.not.(associated(this%mdata) .or. associated(this%sdata))) this%ndata = 0
       call add(errmsg,1,"kernel matrix already allocated: deallocating it now",myname)
    end if
!
    if(ndata .le. 0 .or. nparam .le. 0) then
       write(errstr,*) "requested number of data (rows) = ",ndata,", or number of parameters (columns) = ",nparam,&
            " invalid: must be positive numbers"
       call add(errmsg,2,errstr,myname)
       return
    end if
    ! make it consistend with this%mdata or this%sdata, if already defined
    if(this%ndata > 0) then
       if(ndata /= this%ndata) then
          write(errstr,*) "requested number of data (rows) = ",ndata,&
               ", is inconsistend with existing number of data of this system = ",this%ndata
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    if(nsmooth < 0) then
       write(errstr,*) "requested number of smoothing conditions = ",nsmooth," must not be negative"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nrow = ndata + nsmooth
    ncol = nparam
!
    allocate(this%K(nrow,ncol),stat=status)
    if(status/=0) then
       write(errstr,*) "could not allocate kernel matrix for ",nrow," rows and ",ncol," columns, allocate status = ",status
       call add(errmsg,2,errstr,myname)
       return
    end if
    this%nrow = nrow
    this%ndata = ndata
    this%nsmooth = nsmooth
    this%ncol = ncol
  end subroutine allocateMatrixSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief setup kernel matrix (not parallelized)
!! \details loop over paths and read in only kernel frequencies which are needed for current path, 
!!  then place the kernel values for those data samples in the respective rows of the kernel matrix,
!!  which are in same order as the data samples
!! \param this kernel system
!! \param path_event_filter path were event filter files are
!! \param path_station_filter path were station filter files are
!! \param df_measured_data frequency step of measured data
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param path_sensitivity_kernels path where sensitivity kernel files are
!! \param comptrans component_transformation object for the stations in use
!! \param ntot_invgrid total number of inversion grid cells of the inversion grid in use (used to check kernel files)
!! \param lu1 first file unit
!! \param lu2 second file unit
!! \param dmspace data model space definition
!! \param errmsg error message
!
  subroutine readMatrixSerialKernellLinearSystem(this,path_event_filter,path_station_filter,df_measured_data,&
       nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,comptrans,ntot_invgrid,&
       lu1,lu2,dmspace,errmsg)
    ! incoming
    type (kernel_linear_system) :: this
    character(len=*) :: path_event_filter,path_station_filter,path_sensitivity_kernels
    real :: df_measured_data
    integer, dimension(:) ::  ifreq_measured_data
    integer :: nfreq_measured_data,ntot_invgrid,lu1,lu2
    type (component_transformation) :: comptrans
    type (data_model_space_info) :: dmspace
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=400)  :: errstr
    character(len=35) :: myname = 'readMatrixSerialKernellLinearSystem'
    type (error_message) :: errmsg2,errmsg3
    double precision, dimension(:,:), pointer :: trans_coef
    integer :: j,ndata_count,ndata_div_50
    logical, dimension(:), allocatable :: idata_added
    ! data space
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    integer :: ipath,npath,ncomp_all,idata,icomp,ndata_path_jf,ndata_path_jf_count
    character(len=character_length_evid) :: evid
    character(len=character_length_evid), dimension(:), pointer :: pevid
    character(len=character_length_staname) :: staname
    character(len=character_length_staname), dimension(:), pointer :: pstaname
    character(len=character_length_component), dimension(:), pointer :: comp_all,pcomp
    character(len=2), dimension(:), pointer :: pimre
    integer, dimension(:), pointer :: idata_path_jf,idata_path_jf_comp,pifreq
    ! model space
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name
    integer :: iparam_pmtrz,nparam_pmtrz
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: idx_param,idx_cell,pcell
    type (integer_vector_pointer), dimension(:,:), allocatable :: idx_param_cell
    ! filtering
    complex, dimension(:), pointer :: event_filter,station_comp_filter
    integer :: ifilter
    complex :: filter_value
    ! kernels
    character(len=500) :: file_kernel
    type (spectral_waveform_kernel) :: kernel
    integer :: jf
    complex, dimension(:,:,:), pointer :: kernel_values
!
    call addTrace(errmsg,myname)
    call add(errmsg,0,"incoming path of source filters : '"//trim(path_event_filter)//"'",myname)
    call add(errmsg,0,"incoming path of station filters: '"//trim(path_station_filter)//"'",myname)
    call add(errmsg,0,"incoming path of spectral sensitivity kernels: '"//trim(path_sensitivity_kernels)//"'",myname)
!
    if(.not.associated(this%K)) then
       call add(errmsg,2,"kernel matrix not allocated yet, allocate first before reading",myname)
       return
    end if
    if(this%ndata /= .ndata.dmspace) then
       write(errstr,*) "number of data samples in data space = ",.ndata.dmspace," does not match number of data = ",&
            this%ndata," for which this kernel matrix was allocated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(this%ncol /= .nparam.dmspace) then
       write(errstr,*) "number of model parameters in model space = ",.nparam.dmspace," does not match number of columns = ",&
            this%ndata," for which this kernel matrix was allocated"
       call add(errmsg,2,errstr,myname)
       return
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
    parametrization = .pmtrz.dmspace
    ! get all paths contained in data space
    paths => getPathsDataModelSpaceInfo(dmspace)
    if(.not.associated(paths)) then
       call add(errmsg,2,"no paths were returned by data_model_space_info object",myname)
       return
    end if
    npath = size(paths,2)
    ! get all different model parameters contained in model space, and corresponding invgrid cell indices
    nparam_pmtrz = numberOfParamModelParametrization(parametrization)
    allocate(idx_param_cell(2,nparam_pmtrz))
    iparam_pmtrz = 0
    do while(nextParamModelParametrization(parametrization,param_name))
       iparam_pmtrz = iparam_pmtrz + 1
       allocate(pparam(1)); pparam(1) = param_name
       idx_param => getIndxModelParam(dmspace,param=pparam,cell=pcell)
       call associateVectorPointer(idx_param_cell(1,iparam_pmtrz),idx_param)
       call associateVectorPointer(idx_param_cell(2,iparam_pmtrz),pcell)
       nullify(idx_param,pcell)
       deallocate(pparam)
    end do !  while nextParamModelParametrization
!
    nullify(event_filter,station_comp_filter)
    evid = ''; staname = ''
    nullify(pevid,pstaname,pifreq,pcomp,pimre)
    nullify(idata_path_jf,idata_path_jf_comp,trans_coef)
!
    comp_all => allComp(dmspace)
    if(.not.associated(comp_all)) then
       ! actually cannot be, as then .ndata.dmspace == 0
       call add(errmsg,2,"no components contained in data space, i.e. no data samples",myname)
       goto 1
    end if
    ncomp_all = size(comp_all)
!
    ! keep track on which rows were added to the kernel matrix (to avoid rows which are not initiated, might still be equal to zero on exit)
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
    ndata_count = 0
    ndata_div_50 = this%ndata / 50 + 1 ! give at most 50 messages of progress
!
    ! loop over paths
    do ipath = 1,npath
!
       ! read in event_filter for this path, if necessary
       if(evid /= paths(1,ipath)) then
          evid = paths(1,ipath)
          if(associated(event_filter)) deallocate(event_filter)
          errmsg2 = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu1,event_filter,ndata=nfreq_measured_data)
          call addTrace(errmsg2,myname)
          if(.level.errmsg2 /= 0) then
             write(errstr,*) "error or warning in path ("//trim(paths(1,ipath))//","//trim(paths(2,ipath))//&
                  ") reading source filter from ascii file '"//trim(path_event_filter)//&
                  "filter_"//trim(evid)
                call add(errmsg,2,trim(errstr),myname)
                call print(errmsg2)
                goto 1
             endif
             call dealloc(errmsg2)
       end if
       ! if necessary, get new transformation coefficients
       if(staname /= paths(2,ipath)) then
          staname = paths(2,ipath)
          ! get new transformation coefficients
          if(associated(trans_coef)) deallocate(trans_coef)
          trans_coef => transform(comptrans,(/'CX','CY','CZ'/),comp_all,staname) ! trans_coef should have size (ncomp_all,3)
          if(.not.associated(trans_coef)) then
             write(errstr,*) "'"//comp_all//"', "
             call add(errmsg,2,"there are no transformation coefficients for station '"//&
                  trim(staname)//"' and components "//trim(errstr),myname)
             goto 1
          end if
       end if
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       call add(errmsg2,0,"adding kernel values to system for path ("//trim(evid)//","//trim(staname)//")",myname)
!
       ! open kernel file
       file_kernel = trim(path_sensitivity_kernels)//"spectral_kernel_"//trim(parametrization)//&
            "_"//trim(evid)//"_"//trim(staname)
       call initialReadSpectralWaveformKernel(kernel,file_kernel,lu1,errmsg2)
       if(.level.errmsg2 /= 0) then
          call add(errmsg,2,"error or warning in path ("//trim(evid)//","//trim(staname)//&
                  ") initially reading kernel file '"//trim(file_kernel)//"'",myname)
          call print(errmsg2)
          goto 1
       end if
       call add(errmsg2,0,"successfully openend kernel file '"//trim(file_kernel)//"'",myname)
!
       ! check content of kernel file
       if(parametrization /= .pmtrz.kernel) then
          call add(errmsg2,2,"the current parametrization '"//trim(parametrization)//&
               "' differs from the parametrization '"//trim(.pmtrz.kernel)//"' contained in the kernel file",myname)
          call print(errmsg2)
          write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(ntot_invgrid /= .ntot.kernel) then
          write(errstr,*) "number of kernel values (",.ntot.kernel,") differs from current number of inversion grid cells (",&
               ntot_invgrid,"), which suggests that kernels were computed w.r.t. an inversion grid different from that in use!"
          call add(errmsg2,2,errstr,myname)
          call print(errmsg2)
          write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if( (df_measured_data-.df.kernel)/df_measured_data > 1.e-4) then
          write(errstr,*) "the frequency step of the frequencies in the kernel file (",.df.kernel,&
               ") differs significantly from the frequency step of the measured data (",df_measured_data,&
               "), which suggests that the kernels were computed w.r.t. a different frequency discretization"
          call add(errmsg2,2,errstr,myname)
          call print(errmsg2)
          write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       ! iterate over all frequencies contained in the kernel
       do while (nextFrequencySpectralWaveformKernel(kernel,jf))
!
          ! check if this frequency is contained in data space for the current path
          ! if so, use the respective indices of data samples to fill the respective rows of the kernel matrix
          ! if not, cycle to next frequency
          if(associated(pevid)) deallocate(pevid)
          if(associated(pstaname)) deallocate(pstaname)
          if(associated(pifreq)) deallocate(pifreq)
          if(associated(idata_path_jf)) deallocate(idata_path_jf)
          allocate(pevid(1),pstaname(1),pifreq(1))
          pevid(1) = evid; pstaname(1) = staname; pifreq(1) = jf
          idata_path_jf => getIndxDataSamples(dmspace,evid=pevid,staname=pstaname,ifreq=pifreq)
          if(.not.associated(idata_path_jf)) then
             write(errstr,*) "frequency index ",jf,", of this kernel file is not contained in data space"
             call add(errmsg2,0,errstr,myname)
             cycle
          else
             ndata_path_jf_count = 0 ! count the number of data samples which are added for this path and frequency
             ndata_path_jf = size(idata_path_jf)
             write(errstr,*) "there are ",ndata_path_jf," data samples in data space for this path and frequency index ",jf
             call add(errmsg2,0,errstr,myname)
          end if
!
          ! find correct index of filter value for this frequency
          ifilter = -1
          do j = 1,nfreq_measured_data
             if(ifreq_measured_data(j) == jf) then
                ifilter = j
                exit
             end if
          end do ! j
          if(ifilter .le. 0) then
             write(errstr,*) "frequency index ",jf,", of this kernel file is not contained in the ",&
                  "vector of frequency indices of measured data = ",ifreq_measured_data
             call add(errmsg2,2,errstr,myname)
             call print(errmsg2)
             write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          ! read kernel for this frequency
          call readSpectralWaveformKernel(kernel,jf,errmsg2)
          if(.level.errmsg2 /= 0) then
             write(errstr,*) "error or warning in path ("//trim(evid)//","//trim(staname)//&
                  ") reading frequency ",jf," of kernel"
             call add(errmsg,2,errstr,myname)
             call print(errmsg2)
             goto 1
          end if
          write(errstr,*) "successfully read frequency ",jf," of kernel"
          call add(errmsg2,0,errstr,myname)
          kernel_values => getValuesSpectralWaveformKernel(kernel)
!
          do icomp = 1,ncomp_all
             if(associated(pcomp)) deallocate(pcomp)
             if(associated(pimre)) deallocate(pimre); nullify(pimre)
             if(associated(idata_path_jf_comp)) deallocate(idata_path_jf_comp)
             allocate(pcomp(1)); pcomp(1) = comp_all(icomp)
             idata_path_jf_comp => getIndxDataSamples(dmspace,comp=pcomp,imre=pimre,indx_in=idata_path_jf)
!print *, "sample:: evid = '"//trim(evid)//"', staname = '"//trim(staname)//"', jf = ",jf,&
!", comp = '"//trim(comp_all(icomp))//"', idata_path_jf = ",idata_path_jf,", idata_path_jf_comp = ",idata_path_jf_comp
             if(.not.associated(idata_path_jf_comp)) then
                call add(errmsg2,0,"there are no data samples in data space for this path and frequency and component '"//&
                     trim(comp_all(icomp))//"'",myname)
                cycle
             else
                write(errstr,*) "there are ",size(idata_path_jf_comp),&
                     " data samples in data space for this path and frequency and component '"//trim(comp_all(icomp))//"'"
                call add(errmsg2,0,errstr,myname)
             end if
!
             ! read in station_comp_filter for this path and component
             if(associated(station_comp_filter)) deallocate(station_comp_filter)
             errmsg3 = readAsciiData(trim(path_station_filter)//"filter_"//trim(staname)//"_"//trim(comp_all(icomp)),&
                  lu2,station_comp_filter,ndata=nfreq_measured_data)
             call addTrace(errmsg3,myname)
             if(.level.errmsg3 /= 0) then
                call add(errmsg2,2,"error or warning reading station filter from ascii file '"//trim(path_station_filter)//&
                     "filter_"//trim(staname)//"_"//trim(comp_all(icomp))//"'",myname)
                call print(errmsg3)
                call print(errmsg2)
                write(errstr,*) "error or warning in path ("//trim(evid)//","//trim(staname)//&
                     ") reading station filter for component '"//trim(comp_all(icomp))//"'"
                call add(errmsg,2,errstr,myname)
                goto 1
             endif
             call dealloc(errmsg3)
             ! define filter_value for this path, frequency and component
             filter_value = event_filter(ifilter)*station_comp_filter(ifilter)
!
             do idata=1,size(idata_path_jf_comp) ! in a sensible setup, size(idata_path_jf_comp) == 2, one 'im' and one 're' sample
!print *, "sample:: evid = '"//trim(evid)//"', staname = '"//trim(staname)//"', jf = ",jf,&
!", comp = '"//trim(comp_all(icomp))//"', imre = '"//pimre(idata)//"'"
!                
                ! check, if for any reasons this row was already filled
                if(idata_added(idata_path_jf_comp(idata))) then
                   write(errstr,*) idata_path_jf_comp(idata),"'th row of kernel matrix was already filled with kernel values. ",&
                        "data_model_space_info object may be corrupt"
                   call add(errmsg2,2,errstr,myname)
                   call print(errmsg2)
                   write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
                   call add(errmsg,2,errstr,myname)
                   goto 1                   
                end if
!
                ! now, for this particular row, loop on all parameters of this parametrization and see
                ! if they occur in the model space. If yes, the corresponding model parameter indices and 
                ! inversion grid cell indices are contained in vector pointer array idx_param_cell, that was initiated above
                do iparam_pmtrz = 1,nparam_pmtrz
                   idx_param => getVectorPointer(idx_param_cell(1,iparam_pmtrz))
                   idx_cell => getVectorPointer(idx_param_cell(2,iparam_pmtrz))
                   if(associated(idx_param)) then
                      ! transform the CX,CY,CZ components of the kernel to current component (by matmul with transformation coefficients)
                      ! filter with event_filter and station_comp_filter to bring syntetics to the form of the data
                      !
                      select case (pimre(idata))
                      case('im')
                         this%K(idata_path_jf_comp(idata),idx_param) = imag(filter_value* ( &
                              matmul( kernel_values(idx_cell,:,iparam_pmtrz) , real(trans_coef(icomp,:)) ) &
                              ))
                      case('re')
                         this%K(idata_path_jf_comp(idata),idx_param) = real(filter_value* ( &
                              matmul( kernel_values(idx_cell,:,iparam_pmtrz) , real(trans_coef(icomp,:)) ) &
                              ))
                      end select
                   end if
                end do ! iparam_pmtrz
!
                ! remember this row was filled
                ndata_count = ndata_count + 1
                idata_added(idata_path_jf_comp(idata)) = .true.
                if(mod(ndata_count,ndata_div_50)==0) write(*,*) "read serial kernel matrix row ", &
                     ndata_count," out of",this%ndata,"; ",real(ndata_count)/real(this%ndata)*100.,"%"
!
             end do ! idata
!
          end do ! icomp
!
       end do ! while(nextFrequencySpectralWaveformKernel(kernel,jf))
!
       call finalReadSpectralWaveformKernel(kernel)
       call dealloc(errmsg2)
    end do ! ipath
!
    ! check if there are any rows which were not filled
    if(.not.all(idata_added)) then
       write(errstr,*) "maybe some data samples occur more than once in the data space (not allowed!), as ",&
            count(.not.idata_added),&
            " rows of the kernel matrix were not filled with values, namely these ones (by index) :",&
            pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(paths)) deallocate(paths)
    if(associated(pevid)) deallocate(pevid)
    if(associated(pstaname)) deallocate(pstaname)
    if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(comp_all)) deallocate(comp_all)
    if(associated(pimre)) deallocate(pimre)
    if(associated(idata_path_jf)) deallocate(idata_path_jf)
    if(associated(idata_path_jf_comp)) deallocate(idata_path_jf_comp)
!
    if(allocated(idx_param_cell)) then
       do iparam_pmtrz=1,nparam_pmtrz
          call dealloc(idx_param_cell(1,iparam_pmtrz))
          call dealloc(idx_param_cell(2,iparam_pmtrz))
       end do
       deallocate(idx_param_cell)
    end if
!
    if(associated(event_filter)) deallocate(event_filter)
    if(associated(station_comp_filter)) deallocate(station_comp_filter)
    if(associated(trans_coef)) deallocate(trans_coef)
    if(allocated(idata_added)) deallocate(idata_added)
!
  end subroutine readMatrixSerialKernellLinearSystem
!------------------------------------------------------------------------
!> \brief read in synthetic data files for system right-hand-side (not parallelized)
!! \details Row-wise construction, assume that the data samples
!!  in dmspace object are in an intelligent order, preventing unnecessary opening/closing
!!  of files. 
!! \param this kernel system
!! \param path_event_filter path were event filter files are
!! \param path_station_filter path were station filter files are
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param nfreq_synthetic_data number of frequencies of synthetic data (must equal length of vector ifreq_synthetic_data)
!! \param ifreq_synthetic_data vector containing all indices of frequencies of synthetic data (used to interpret filter files)
!! \param path_synthetic_data path where synthetic data files are
!! \param comptrans component_transformation object for the stations in use
!! \param lu file unit
!! \param dmspace data model space definition
!! \param errmsg error message
!
  subroutine readSyntheticDataSerialKernelLinearSystem(this,path_event_filter,path_station_filter,&
       nfreq_measured_data,ifreq_measured_data,nfreq_synthetic_data,ifreq_synthetic_data,&
       path_synthetic_data,comptrans,lu,dmspace,errmsg)
    ! incoming
    type (kernel_linear_system) :: this
    character(len=*) :: path_event_filter,path_station_filter,path_synthetic_data
    integer, dimension(:) ::  ifreq_measured_data,ifreq_synthetic_data
    integer :: nfreq_measured_data,nfreq_synthetic_data,lu
    type (component_transformation) :: comptrans
    type (data_model_space_info) :: dmspace
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=41) :: myname='readSyntheticDataSerialKernelLinearSystem'
    character(len=400)  :: errstr
    type (error_message) :: errmsg2,errmsg3
    double precision, dimension(:,:), pointer :: trans_coef
    integer :: j
    logical, dimension(:), allocatable :: idata_added
    ! data space
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    integer :: ndata,ipath,npath,icomp,ncomp_all,ndata_path_comp
    character(len=character_length_evid) :: evid
    character(len=character_length_evid), dimension(:), pointer :: pevid
    character(len=character_length_staname) :: staname
    character(len=character_length_staname), dimension(:), pointer :: pstaname
    character(len=character_length_component), dimension(:), pointer :: comp_all,pcomp
    character(len=2), dimension(:), pointer :: pimre
    integer, dimension(:), pointer :: pifreq,idata_path_comp
    ! filtering
    complex, dimension(:), pointer :: event_filter,station_comp_filter
    integer, dimension(:), allocatable :: map_ifreq_to_synthetic_index,map_ifreq_to_measured_index
    ! synthetic data
    complex, dimension(:,:), pointer :: sdata_all_comp
    complex, dimension(:), allocatable :: sdata_transform_filter
!
    call addTrace(errmsg,myname)
    call add(errmsg,0,"incoming path of source filters : '"//trim(path_event_filter)//"'",myname)
    call add(errmsg,0,"incoming path of station filters: '"//trim(path_station_filter)//"'",myname)
    call add(errmsg,0,"incoming path of synthetic data: '"//trim(path_synthetic_data)//"'",myname)
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
    if(associated(this%sdata)) then
       deallocate(this%sdata)
       if(.not.(associated(this%K) .or. associated(this%mdata))) this%ndata = 0
       call add(errmsg,1,"synthetic data vector already allocated: deallocating it now",myname)
    end if
!
    ndata = .ndata.dmspace
    if(ndata .le. 0) then
       write(errstr,*) "requested number of data  = ",ndata," (contained in data space) invalid: must be a positive number"
       call add(errmsg,2,errstr,myname)
       return
    end if
    ! make it consistend with this%mdata or this%K, if already defined
    if(this%ndata > 0) then
       if(ndata /= this%ndata) then
          write(errstr,*) "requested number of data = ",ndata,&
               ", (contained in data space) is inconsistend with existing number of data of this system = ",this%ndata
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    this%ndata = ndata
    allocate(this%sdata(ndata))
!
    ! get all paths contained in data space
    paths => getPathsDataModelSpaceInfo(dmspace)
    if(.not.associated(paths)) then
       call add(errmsg,2,"no paths were returned by data_model_space_info object",myname)
       return
    end if
    npath = size(paths,2)
!
    nullify(event_filter,station_comp_filter)
    evid = ''; staname = ''
    nullify(pevid,pstaname,pifreq,pcomp,pimre)
    nullify(idata_path_comp,trans_coef)
!
    comp_all => allComp(dmspace)
    if(.not.associated(comp_all)) then
       ! actually cannot be, as then .ndata.dmspace == 0
       call add(errmsg,2,"no components contained in data space, i.e. no data samples",myname)
       goto 1
    end if
    ncomp_all = size(comp_all)
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
    ! keep track on which rows were added to the kernel matrix (to avoid rows which are not initiated, might still be equal to zero on exit)
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
!
    ! loop over paths
    do ipath = 1,npath
!
       ! read in event_filter for this path, if necessary
       if(evid /= paths(1,ipath)) then
          evid = paths(1,ipath)
          if(associated(event_filter)) deallocate(event_filter)
          errmsg2 = readAsciiData(trim(path_event_filter)//"filter_"//trim(evid),lu,event_filter,ndata=nfreq_measured_data)
          call addTrace(errmsg2,myname)
          if(.level.errmsg2 /= 0) then
             write(errstr,*) "error or warning in path ("//trim(paths(1,ipath))//","//trim(paths(2,ipath))//&
                  ") reading source filter from ascii file '"//trim(path_event_filter)//&
                  "filter_"//trim(evid)
             call add(errmsg,2,trim(errstr),myname)
             call print(errmsg2)
             goto 1
          endif
          call dealloc(errmsg2)
       end if
       ! if necessary, get new transformation coefficients
       if(staname /= paths(2,ipath)) then
          staname = paths(2,ipath)
          ! get new transformation coefficients
          if(associated(trans_coef)) deallocate(trans_coef)
          trans_coef => transform(comptrans,(/'CX','CY','CZ'/),comp_all,staname) ! trans_coef should have size (ncomp_all,3)
          if(.not.associated(trans_coef)) then
             write(errstr,*) "'"//comp_all//"', "
             call add(errmsg,2,"there are no transformation coefficients for station '"//&
                  trim(staname)//"' and components "//trim(errstr),myname)
             goto 1
          end if
       end if
!
       ! read in synthetic data for this path
       if(associated(sdata_all_comp)) deallocate(sdata_all_comp)
       errmsg2 = readAsciiData(trim(path_synthetic_data)//"synthetics_"//trim(evid)//"_"//trim(staname),&
                  lu,sdata_all_comp,3,nrow=nfreq_synthetic_data)
       call addTrace(errmsg2,myname)
       if(.level.errmsg2 /= 0) then
          call print(errmsg2)
          call add(errmsg,2,"error or warning reading synthetic data from ascii file '"//trim(path_synthetic_data)//&
               "synthetics_"//trim(evid)//"_"//trim(staname)//"' in path ("//trim(evid)//","//trim(staname)//")",myname)
          goto 1
       end if
       call dealloc(errmsg2)
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       call add(errmsg2,0,"treating path ("//trim(evid)//","//trim(staname)//")",myname)
!
       do icomp = 1,ncomp_all
!
          ! check if data space contains samples for this component and current path
          if(associated(pevid)) deallocate(pevid)
          if(associated(pstaname)) deallocate(pstaname)
          if(associated(pcomp)) deallocate(pcomp)
          if(associated(pifreq)) deallocate(pifreq); nullify(pifreq)
          if(associated(pimre)) deallocate(pimre); nullify(pimre)
          if(associated(idata_path_comp)) deallocate(idata_path_comp)
          allocate(pevid(1)); allocate(pstaname(1)); allocate(pcomp(1))
          pevid(1) = evid; pstaname(1) = staname; pcomp(1) = comp_all(icomp)
          idata_path_comp => getIndxDataSamples(dmspace,evid=pevid,staname=pstaname,ifreq=pifreq,comp=pcomp,imre=pimre)
          if(.not.associated(idata_path_comp)) then
             write(errstr,*) "component '",trim(comp_all(icomp)),"' is not contained in data space for this path"
             call add(errmsg2,0,errstr,myname)
             cycle
          else
             ndata_path_comp = size(idata_path_comp)
             write(errstr,*) "there are ",ndata_path_comp,&
                  " data samples in data space for this path and component '"//trim(comp_all(icomp))//"'"
             call add(errmsg2,0,errstr,myname)
          end if

          ! read in station_comp_filter for this path and component
          if(associated(station_comp_filter)) deallocate(station_comp_filter)
          errmsg3 = readAsciiData(trim(path_station_filter)//"filter_"//trim(staname)//"_"//trim(comp_all(icomp)),&
               lu,station_comp_filter,ndata=nfreq_measured_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) then
             call add(errmsg2,2,"error or warning reading station filter from ascii file '"//trim(path_station_filter)//&
                  "filter_"//trim(staname)//"_"//trim(comp_all(icomp))//"'",myname)
             call print(errmsg3)
             call print(errmsg2)
             write(errstr,*) "error or warning in path ("//trim(evid)//","//trim(staname)//&
                  ") reading station filter for component '"//trim(comp_all(icomp))//"'"
             call add(errmsg,2,errstr,myname)
             goto 1
          endif
          call dealloc(errmsg3)
!
          ! allocate temporary data array sdata_transform_filter to perform transformation and filtering on required data samples
          if(allocated(sdata_transform_filter)) deallocate(sdata_transform_filter)
          allocate(sdata_transform_filter(ndata_path_comp))
!
          ! initiate sdata_transform_filter by transforming to current component
          sdata_transform_filter(:) = matmul(sdata_all_comp(map_ifreq_to_synthetic_index(pifreq),:) , trans_coef(icomp,:) )
!
          ! filter the transformed sdata
          sdata_transform_filter = sdata_transform_filter * &
               event_filter(map_ifreq_to_measured_index(pifreq)) * &
               station_comp_filter(map_ifreq_to_measured_index(pifreq))
!
          ! check, if for any reasons the current values of synthetic data were already defined
          if(any(idata_added(idata_path_comp))) then
             write(errstr,*) "there are ",count(idata_added(idata_path_comp))," synthetic data sample values which are ",&
                  "already defined. data_model_space_info object may be corrupt"
             call add(errmsg2,2,errstr,myname)
             call print(errmsg2)
             write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
             call add(errmsg,2,errstr,myname)
             goto 1                   
          end if
!
          ! finally set values in this%sdata
          where(pimre == 'im')
             this%sdata(idata_path_comp) = imag(sdata_transform_filter)
             idata_added(idata_path_comp) = .true.
          elsewhere
             this%sdata(idata_path_comp) = real(sdata_transform_filter)
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
            count(.not.idata_added),&
            " synthetic data values were not set, namely these ones (by index) :",&
            pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(paths)) deallocate(paths)
    if(associated(pevid)) deallocate(pevid)
    if(associated(pstaname)) deallocate(pstaname)
    if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(pimre)) deallocate(pimre)
    if(associated(comp_all)) deallocate(comp_all)
    if(associated(idata_path_comp)) deallocate(idata_path_comp)
!
    if(associated(event_filter)) deallocate(event_filter)
    if(associated(station_comp_filter)) deallocate(station_comp_filter)
    if(associated(trans_coef)) deallocate(trans_coef)
    if(allocated(idata_added)) deallocate(idata_added)
!
    if(associated(sdata_all_comp)) deallocate(sdata_all_comp)
    if(allocated(sdata_transform_filter)) deallocate(sdata_transform_filter)
!
    if(allocated(map_ifreq_to_synthetic_index)) deallocate(map_ifreq_to_synthetic_index)
    if(allocated(map_ifreq_to_measured_index)) deallocate(map_ifreq_to_measured_index)
!
  end subroutine readSyntheticDataSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief read in synthetic data files for system right-hand-side (not parallelized)
!! \details Row-wise construction, assume that the data samples
!!  in dmspace object are in an intelligent order, preventing unnecessary opening/closing
!!  of files. 
!! \param this kernel system
!! \param nfreq_measured_data number of frequencies of measured data (must equal length of vector ifreq_measured_data)
!! \param ifreq_measured_data vector containing all indices of frequencies of measured data (used to interpret filter files)
!! \param path_measured_data path where measured data files are
!! \param lu file unit
!! \param dmspace data model space definition
!! \param errmsg error message
!
  subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
       path_measured_data,lu,dmspace,errmsg)
    ! incoming
    type (kernel_linear_system) :: this
    integer, dimension(:) ::  ifreq_measured_data
    integer :: nfreq_measured_data,lu
    character(len=*) :: path_measured_data
    type (data_model_space_info) :: dmspace
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=40) :: myname='readMeasuredDataSerialKernelLinearSystem'
    character(len=400)  :: errstr
    type (error_message) :: errmsg2,errmsg3
    integer :: j
    logical, dimension(:), allocatable :: idata_added
    ! data space
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    integer :: ndata,ipath,npath,icomp,ncomp_all,ndata_path_comp
    character(len=character_length_evid) :: evid
    character(len=character_length_evid), dimension(:), pointer :: pevid
    character(len=character_length_staname) :: staname
    character(len=character_length_staname), dimension(:), pointer :: pstaname
    character(len=character_length_component), dimension(:), pointer :: comp_all,pcomp
    character(len=2), dimension(:), pointer :: pimre
    integer, dimension(:), pointer :: pifreq,idata_path_comp
    ! measured data
    complex, dimension(:), pointer :: mdata_cmplx
    integer, dimension(:), allocatable :: map_ifreq_to_measured_index
!
    call addTrace(errmsg,myname)
    call add(errmsg,0,"incoming path of measured data: '"//trim(path_measured_data)//"'",myname)
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
    if(associated(this%mdata)) then
       deallocate(this%mdata)
       if(.not.(associated(this%K) .or. associated(this%sdata))) this%ndata = 0
       call add(errmsg,1,"measured data vector already allocated: deallocating it now",myname)
    end if
!
    ndata = .ndata.dmspace
    if(ndata .le. 0) then
       write(errstr,*) "requested number of data  = ",ndata," (contained in data space) invalid: must be a positive number"
       call add(errmsg,2,errstr,myname)
       return
    end if
    ! make it consistend with this%mdata or this%K, if already defined
    if(this%ndata > 0) then
       if(ndata /= this%ndata) then
          write(errstr,*) "requested number of data = ",ndata,&
               ", (contained in data space) is inconsistend with existing number of data of this system = ",this%ndata
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    this%ndata = ndata
    allocate(this%mdata(ndata))
!
    ! get all paths contained in data space
    paths => getPathsDataModelSpaceInfo(dmspace)
    if(.not.associated(paths)) then
       call add(errmsg,2,"no paths were returned by data_model_space_info object",myname)
       return
    end if
    npath = size(paths,2)
!
    evid = ''; staname = ''
    nullify(pevid,pstaname,pifreq,pcomp,pimre)
    nullify(idata_path_comp)
!
    comp_all => allComp(dmspace)
    if(.not.associated(comp_all)) then
       ! actually cannot be, as then .ndata.dmspace == 0
       call add(errmsg,2,"no components contained in data space, i.e. no data samples",myname)
       goto 1
    end if
    ncomp_all = size(comp_all)
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
    ! keep track on which rows were added to the kernel matrix (to avoid rows which are not initiated, might still be equal to zero on exit)
    allocate(idata_added(this%ndata))
    idata_added(:) = .false.
!
    ! loop over paths
    do ipath = 1,npath
       evid = paths(1,ipath); staname = paths(2,ipath)
!
       ! for every path, use a new error message (otherwise, list of messages will get too long)
       call new(errmsg2,myname)
       call add(errmsg2,0,"treating path ("//trim(evid)//","//trim(staname)//")",myname)
!
       do icomp = 1,ncomp_all
!
          ! check if data space contains samples for this component and current path
          if(associated(pevid)) deallocate(pevid)
          if(associated(pstaname)) deallocate(pstaname)
          if(associated(pcomp)) deallocate(pcomp)
          if(associated(pifreq)) deallocate(pifreq); nullify(pifreq)
          if(associated(pimre)) deallocate(pimre); nullify(pimre)
          if(associated(idata_path_comp)) deallocate(idata_path_comp)
          allocate(pevid(1)); allocate(pstaname(1)); allocate(pcomp(1))
          pevid(1) = evid; pstaname(1) = staname; pcomp(1) = comp_all(icomp)
          idata_path_comp => getIndxDataSamples(dmspace,evid=pevid,staname=pstaname,ifreq=pifreq,comp=pcomp,imre=pimre)
          if(.not.associated(idata_path_comp)) then
             write(errstr,*) "component '",trim(comp_all(icomp)),"' is not contained in data space for this path"
             call add(errmsg2,0,errstr,myname)
             cycle
          else
             ndata_path_comp = size(idata_path_comp)
             write(errstr,*) "there are ",ndata_path_comp,&
                  " data samples in data space for this path and component '"//trim(comp_all(icomp))//"'"
             call add(errmsg2,0,errstr,myname)
          end if
!
          ! read in measured data for this path
          if(associated(mdata_cmplx)) deallocate(mdata_cmplx)
          errmsg3 = readAsciiData(trim(path_measured_data)//"data_"//trim(evid)//"_"//trim(staname)//"_"//trim(comp_all(icomp)),&
                  lu,mdata_cmplx,ndata=nfreq_measured_data)
          call addTrace(errmsg3,myname)
          if(.level.errmsg3 /= 0) then
             call print(errmsg3)
             call add(errmsg2,2,"error or warning reading measured data from ascii file '"//trim(path_measured_data)//&
                  "data_"//trim(evid)//"_"//trim(staname)//trim(comp_all(icomp))//"'",myname)
             call print(errmsg2)
             write(errstr,*) "error or warning in path ("//trim(evid)//","//trim(staname)//&
                  ") reading measured data for component '"//trim(comp_all(icomp))//"'"
             call add(errmsg,2,errstr,myname)
          goto 1
       end if
       call dealloc(errmsg3)
!
       ! check, if for any reasons the current values of synthetic data were already defined
       if(any(idata_added(idata_path_comp))) then
          write(errstr,*) "there are ",count(idata_added(idata_path_comp))," measured data sample values which are ",&
               "already defined. data_model_space_info object may be corrupt"
          call add(errmsg2,2,errstr,myname)
          call print(errmsg2)
          write(errstr,*) "error or warning in ",ipath,"'th path ("//trim(evid)//","//trim(staname)//")"
          call add(errmsg,2,errstr,myname)
          goto 1                   
       end if
!
       ! finally set values in this%mdata
       where(pimre == 'im')
          this%mdata(idata_path_comp) = imag(mdata_cmplx(map_ifreq_to_measured_index(pifreq)))
          idata_added(idata_path_comp) = .true.
       elsewhere
          this%mdata(idata_path_comp) = real(mdata_cmplx(map_ifreq_to_measured_index(pifreq)))
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
            count(.not.idata_added),&
            " synthetic data values were not set, namely these ones (by index) :",&
            pack( (/ (j,j=1,this%ndata) /) , .not.idata_added )
       call add(errmsg,2,errstr,myname)
    end if
!
    ! clean up
!
1   if(associated(paths)) deallocate(paths)
    if(associated(pevid)) deallocate(pevid)
    if(associated(pstaname)) deallocate(pstaname)
    if(associated(pifreq)) deallocate(pifreq)
    if(associated(pcomp)) deallocate(pcomp)
    if(associated(pimre)) deallocate(pimre)
    if(associated(comp_all)) deallocate(comp_all)
    if(associated(idata_path_comp)) deallocate(idata_path_comp)
!
    if(allocated(idata_added)) deallocate(idata_added)
    if(associated(mdata_cmplx)) deallocate(mdata_cmplx)
    if(allocated(map_ifreq_to_measured_index)) deallocate(map_ifreq_to_measured_index)
!
  end subroutine readMeasuredDataSerialKernelLinearSystem
!------------------------------------------------------------------------
!> \brief define one right hand side vector of the linear system by data residual
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
    if(.not.associated(this%K)) then
       call add(errmsg,2,"Kernel matrix is not yet allocated. Please allocate before setting residuals",myname)
       return
    end if
    if(.not.(associated(this%mdata).and.associated(this%sdata))) then
       call add(errmsg,2,"synthetic data or measured data not yet read in. required for computing residuals",myname)
       return
    end if
    if(associated(this%rhs)) then
       call add(errmsg,1,"right hand side(s) were already allocated, deallocating now before settin to new ones",myname)
       deallocate(this%rhs)
    end if
!
    allocate(this%rhs(this%nrow,1)); this%nrhs = 1
    this%rhs(:,1) = 0.
!
    this%rhs(1:this%ndata,1) = this%mdata - this%sdata
  end subroutine setRhsAsDataResidualKernelLinearSystem
!------------------------------------------------------------------------
!> \brief solve kernel linear system (not parallelized)
!! \details Use the matrix K and the residual vector as the right hand side to form a 
!!  linear system. Compute the solution of that system, which is set as dparam
!! \param this kernel system
!! \param errmsg error message
!! \return error message
!
  subroutine solveSerialKernellLinearSystem(this,errmsg)
    ! incoming
    type (kernel_linear_system) :: this
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=34) :: myname='solveSequentialKernellLinearSystem'
    type (error_message) :: errmsg2
    ! linear system
    type (linear_system) :: LSE
    real, dimension(:,:), pointer :: b,solution,residual
!
    call addTrace(errmsg,myname)
!
    if(.not.associated(this%K)) then
       call add(errmsg,2,"system matrix K not associated",myname)
       return
    endif
    if(.not.associated(this%rhs)) then
       call add(errmsg,2,"right-hand-side vector(s) not associated",myname)
       return
    endif
!
    ! make a copy of rhs, not to overwrite it
    allocate(b(this%nrow,this%nrhs)); b = this%rhs
    nullify(solution,residual)
!
    call createLinearSystem(LSE,this%nrow,this%ncol,1,this%K,b)
    errmsg2 = solveLeastSquaresLinearSystem(LSE,solution,residual,override_A=.false.)
    call addTrace(errmsg2,myname)
    if(.level.errmsg2/=0) then
       call print(errmsg2)
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
    if(associated(this%sol)) deallocate(this%sol)
    this%sol => solution
    nullify(solution)
!
    deallocate(b)
    call dealloc(LSE)
  end subroutine solveSerialKernellLinearSystem
!------------------------------------------------------------------------
!> \brief deallocate kernel linear system
!! \param this kernel linear system to be deallocated
!
  subroutine deallocateKernelLinearSystem(this)
    type (kernel_linear_system) :: this
    if (associated(this%K)) deallocate(this%K)
    if (associated(this%rhs)) deallocate(this%rhs)
    if (associated(this%sol)) deallocate(this%sol)
    this%nrow = 0; this%ndata = 0; this%nsmooth = 0; this%ncol = 0; this%nrhs = 0
    if (associated(this%mdata)) deallocate(this%mdata)
    if (associated(this%sdata)) deallocate(this%sdata)
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
!> \brief return pointer to solution of kernel linear system (i.e. model uptdate)
!! \param this kernel linear system
!! \param dparam pointer to solution of kernel linear system (i.e. model uptdate)
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
!> \brief return number of smoothing conditions nsmooth which were added to this kernel linear system
  function getNsmoothKernelLinearSystem(this) result(n)
    type (kernel_linear_system), intent(in) :: this
    integer :: n
    n = this%nsmooth
  end function getNsmoothKernelLinearSystem
!
end module kernelLinearSystem
