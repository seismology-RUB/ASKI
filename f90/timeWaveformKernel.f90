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
!> \brief module which computes integrated time waveform sensitivity kernels
!!
!! \details By inversely Fourier transforming spectral waveform sensitivity kernels
!!  (or by direct computation from time domain wavefields?! -> not supported yet), time domain
!!  waveform sensitivity kernels are computed or written/read to/from file by this module.
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module timeWaveformKernel
!
  use modelParametrization
  use spectralWaveformKernel
  use kernelDisplacement, only: length_ID_kernel_displacement
  use kernelGreenTensor, only: length_ID_kernel_green_tensor
  use discreteFourierTransform
  use flexibleType
  use streamAccess
  use realloc
  use errorMessage
!
  implicit none
!
  private :: privateReadMetaInfoTimeWaveformKernel
!
  interface dealloc; module procedure deallocateTimeWaveformKernel; end interface
  interface operator (.pmtrz.); module procedure getParametrizationTimeWaveformKernel; end interface
  interface operator (.ntot.); module procedure getNtotKernelTimeWaveformKernel; end interface
  interface operator (.ncomp.); module procedure getNcompTimeWaveformKernel; end interface
  interface operator (.comp.); module procedure getCompTimeWaveformKernel; end interface
  interface operator (.tzero.); module procedure getT0TimeWaveformKernel; end interface
  interface operator (.dt.); module procedure getDtTimeWaveformKernel; end interface
  interface operator (.nt.); module procedure getNtTimeWaveformKernel; end interface
  interface operator (.jt.); module procedure getJtTimeWaveformKernel; end interface
!
!< derived type containing specifications of a time waveform kernel object
  type time_waveform_kernel
     private
     ! model parametrization
     character(len=character_length_pmtrz) :: parametrization = '' !< model parametrization of kernel (also defines number of parameters nparam, below)
     integer :: nparam = 0 !< nparam = number of model parameters for which there are kernel values in this object (size of second dimension of kernel array)
     character(len=character_length_param), dimension(:), pointer :: param => null() !< vector of length nparam indicating the actual names of model parameters for which there are kernel values in this object
     integer :: nparam_in_file = 0 !< number of model parameters contained in kernel file, ONLY USED FOR READING KERNEL FILES!
     integer, dimension(:), pointer :: param_indx_in_file => null() !< mapping ONLY USED FOR READING KERNEL FILES! dim(nparam) . param_indx_in_file(iparam) = index of param used to generate data sets at the time writing the kernel file by stream access
     integer :: ntot_kernel = 0 !< total number of kernel values (either inversion grid cells onto which kernels are pre-integrated, or number of wavefield points, dependent on flag this%on_invgrid)
     ! temporal discretization
     real :: t0 = 0. !< global time shift of times defined by jt*dt, IS SUPPOSED TO BE USED (only), IF INDICES jt GET TOO LARGE FOR INTEGER(4)!!
     real :: dt = -1. !< time step of time discretization
     integer :: nt = 0 !< total number of time samples
     integer, dimension(:), pointer :: jt => null() !< time indices in order of subgroups of kernel file; times compute as t = t0+jt*dt [s]
     integer :: jtcur = -1 !< indicates current time index jt (i.e. time t=t0+jtcur*dt) for which kernel values are contained in array kernel
     ! source
     character(len=length_ID_kernel_displacement) :: kd_id = '' !< .id.kd for kernel_displacement object kd from which this kernel was created
     ! receiver
     character(len=length_ID_kernel_green_tensor) :: kgt_id = '' !< .id.kgt for kernel_green_tensor object kgt from which this kernel was created
     integer :: ncomp = 0 !< number of components for which kernel values are contained in array kernel (2nd dimension)
     character(len=character_length_component), dimension(:), pointer :: comp => null() !< the ncomp component names
     integer :: ncomp_in_file = 0 !< number of receiver components contained in kernel file, ONLY USED FOR READING KERNEL FILES!
     integer, dimension(:), pointer :: comp_indx_in_file => null() !< mapping ONLY USED FOR READING KERNEL FILES! dim(ncomp) . comp_indx_in_file(icomp) = index of comp used to generate data sets at the time writing the kernel file by stream access
     ! kernel values
     logical :: on_invgrid = .true. !< indicates whether this%kernel holds pre-integrated values on inversion grid (.true. , default) or on wavefield points (.false.)
     real, dimension(:,:,:), pointer :: kernel => null()  !< pre-integrated kernel values on invgrid (or plain kernel values on wp) for one time step dim(ntot_kernel,nparam,ncomp)     ! file handling
     ! file handling
     integer :: kernel_stat = -1 !< indicating status of kernel object: -1 = closed and not initiated, 0 = closed and initiated, 1 = open to write, 2 = open to read
     type (file_stream_access) :: fsa
     type (group_stream_access) :: root
  end type time_waveform_kernel
!
contains
!
!------------------------------------------------------------------------
!> \brief initiate time waveform kernel and allocate kernel array
!! \param this time waveform kernel
!! \param parametrization valid model parametrization, will be the parametrization of the kernel
!! \param ntot_invgrid total number of inversion grid cells (kernel values)
!! \param errmsg error message
!
  subroutine initiateTimeWaveformKernel(this,parametrization,param,ntot_kernel,comp,errmsg,kernel_on_wp)
    type (time_waveform_kernel) :: this
    character(len=*) :: parametrization
    integer :: ntot_kernel
    character(len=*), dimension(:) :: comp,param
    type (error_message) :: errmsg
    logical, optional :: kernel_on_wp
    ! local
    character(len=400) :: errstr
    character(len=26) :: myname = 'initiateTimeWaveformKernel'
    logical :: on_invgrid
    integer :: icomp,ncomp,iparam
!
    call addTrace(errmsg,myname)
!
    if(this%kernel_stat /= -1) then
       call add(errmsg,2,"kernel object is already initiated, please deallocate first",myname)
       return
    end if
!
    if(.not.validModelParametrization(parametrization)) then
       call add(errmsg,2,"incoming model parametrization '"//trim(parametrization)//&
            "' is not valid, valid parametrizations(parameters) are "//all_valid_pmtrz_param,myname)
       return
    end if
!
    if(size(param) <= 0) then
       call add(errmsg,2,"there are no incoming model parameters",myname)
       return
    end if
    do iparam = 1,size(param)
       if(.not.validParamModelParametrization(parametrization,param(iparam))) then
          write(errstr,*) iparam,"'th incoming model parameter '",trim(param(iparam)),&
               "' not a valid parameter of parametrization '",trim(parametrization),&
               "'. valid parametrizations (parameters) are "//all_valid_pmtrz_param
          call add(errmsg,2,errstr,myname)
          return
       end if
       ! from the second parameter on, check if the parameter already occured in the list
       if(iparam >=2) then
          if(any(param(1:iparam-1)==param(iparam))) then
             write(errstr,*) iparam,"'th incoming model parameter '",trim(param(iparam)),&
                  "' already occured in the vector of parameters. There must not be duplicate parameter requests."
             call add(errmsg,2,errstr,myname)
             return
          end if ! any duplicate
       end if ! iparam >=2
    end do ! iparam
!
    on_invgrid = .true. ! by default
    if(present(kernel_on_wp)) then
       on_invgrid = .not. kernel_on_wp
    end if
!
    if(ntot_kernel <= 0) then
       if(on_invgrid) then
          write(errstr,*) "incoming number of inversion grid cells ( = ",ntot_kernel,") must be positive!"
       else
          write(errstr,*) "incoming number of wavefield points ( = ",ntot_kernel,") must be positive!"
       end if
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ncomp = size(comp)
    if(ncomp <= 0) then
       write(errstr,*) "incoming number of station components ",ncomp," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(.not.allValidComponents(comp,i_invalid=icomp)) then
       write(errstr,*) icomp,"'th incoming component '"//trim(comp(icomp))//"' not valid. Valid components are '"//&
            all_valid_components//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    do icomp = 2,ncomp
       if(any(comp(1:icomp-1)==comp(icomp))) then
          write(errstr,*) icomp,"'th requested component '",trim(comp(icomp)),"' occurs more than once in the list of the ",&
               ncomp," requested components. there must not be duplicate components in the request"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end do ! icomp
!
    this%parametrization = parametrization
    this%nparam = size(param)
    allocate(this%param(this%nparam))
    this%param = param
    this%ntot_kernel = ntot_kernel
    this%on_invgrid = on_invgrid
    this%ncomp = ncomp
    allocate(this%comp(this%ncomp))
    this%comp = comp
!
    this%jtcur = -1
    allocate(this%kernel(this%ntot_kernel,this%nparam,this%ncomp))
!
    if(this%on_invgrid) then
       write(errstr,*) "initiated pre-integrated time waveform kernel on ",this%ntot_kernel," inversion grid cells, ",&
            this%nparam," parameters of parametrization '",trim(this%parametrization),"', and ",this%ncomp,&
            " station components"
    else
       write(errstr,*) "initiated plain time waveform kernel on ",this%ntot_kernel," wavefield points, ",&
            this%nparam," parameters of parametrization '",trim(this%parametrization),"', and ",this%ncomp,&
            " station components"
    end if
    call add(errmsg,0,errstr,myname)
!
    this%kernel_stat = 0
  end subroutine initiateTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this time waveform kernel
!
  subroutine deallocateTimeWaveformKernel(this)
    type (time_waveform_kernel) :: this
    this%parametrization = ''
    this%nparam = 0
    if(associated(this%param)) deallocate(this%param)
    this%nparam_in_file = 0
    if(associated(this%param_indx_in_file)) deallocate(this%param_indx_in_file)
    this%ntot_kernel = 0
    this%dt = -1.
    this%nt = 0
    if(associated(this%jt)) deallocate(this%jt)
    this%jtcur = -1
    this%kd_id = ''
    this%kgt_id = ''
    this%ncomp = 0
    if(associated(this%comp)) deallocate(this%comp)
    this%ncomp_in_file = 0
    if(associated(this%comp_indx_in_file)) deallocate(this%comp_indx_in_file)
    this%on_invgrid = .true.
    if(associated(this%kernel)) deallocate(this%kernel)
    if(this%kernel_stat>0) then
       call clearGroupTree(this%root)
       call dealloc(this%fsa)
    end if
    this%kernel_stat = -1
  end subroutine deallocateTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief compute time domain kernel from spectral kernel and write to file for all times
!! \details Using module discreteFourierTransform and the spectral waveform kernel
!!  
!! \param 
!
  subroutine transformSpectralToTimeWaveformKernel(skernel,dt,jt,file_tkernel,lu,errmsg,filter,t0)
    type (spectral_waveform_kernel) :: skernel
    real :: dt
    integer, dimension(:) :: jt
    character(len=*) :: file_tkernel
    integer :: lu
    type (error_message) :: errmsg
    ! optional
    complex, dimension(:,:), optional :: filter
    real, optional :: t0
    ! local
    character(len=37) :: myname = 'transformSpectralToTimeWaveformKernel'
    character(len=400) :: errstr
    type (error_message) :: errmsg2
    integer, dimension(:), pointer :: jf
    character(len=character_length_component), dimension(:), pointer :: comp
    character(len=character_length_param), dimension(:), pointer :: param
    type (time_waveform_kernel) :: tkernel
    character(len=length_ID_kernel_displacement) :: kd_id
    character(len=length_ID_kernel_green_tensor) :: kgt_id
    type (discrete_fourier_transform) :: DFT
    integer :: icomp,ncomp,it,nt,nf1,nf2,ifreq,frequency_count,nparam,n,n_shift
    real, dimension(:,:), allocatable :: tkernel_kernel_all_t
    complex, dimension(:,:), pointer :: skernel_vals
    complex, dimension(:), allocatable :: skernel_vals_filtered_reshaped
    logical :: apply_filter
!
    nullify(jf,comp,param,skernel_vals)
!
    call addTrace(errmsg,myname)
!
! first check incoming values and sizes of arrays
!
    if(dt < 0.) then
       write(errstr,*) "incoming time step dt = ",dt," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
    nt = size(jt)
    if(nt <= 0) then
       write(errstr,*) "incoming number of time indices ",nt," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! check content of incoming spectral kernel skernel
!
    kd_id = .kdID.skernel
    if(kd_id == '') then
       call add(errmsg,2,"no kernel_displacement ID contained in spectral kernel object, "//&
            "which should have been opened to read before",myname)
       return
    end if
    kgt_id = .kgtID.skernel
    if(kgt_id == '') then
       call add(errmsg,2,"no kernel_green_tensor ID contained in spectral kernel object, "//&
            "which should have been opened to read before",myname)
       return
    end if
!
    comp => .comp.skernel
    if(.not.associated(comp)) then
       call add(errmsg,2,"no station components contained in spectral kernel object, "//&
            "which should have been opened to read before",myname)
       return
    end if
    ncomp = size(comp)
    param => .param.skernel
    if(.not.associated(param)) then
       call add(errmsg,2,"no model parameters contained in spectral kernel object, "//&
            "which should have been opened to read before",myname)
       return
    end if
    nparam = size(param)
!
    jf => .jf.skernel
    if(.not.associated(jf)) then
       call add(errmsg,2,"no frequency indices contained in spectral kernel object, "//&
            "which should have been opened to read",myname)
       return
    end if
    apply_filter = present(filter)
    if(apply_filter) then
       if(size(filter,1) /= size(jf)) then
          write(errstr,*) "number of incoming filter values per component ",size(filter,1),&
               " does not match number of frequency indices ",size(jf)," contained in the spectral kernel file"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(filter,2) /= ncomp) then
          write(errstr,*) "number of incoming filter values per frequency ",size(filter,2),&
               " does not match number of station components ",ncomp," contained in the spectral kernel"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ! check if frequency discretization of spectral kernel is suitable for inverse fourier transform 
    ! (jf must be nf1:nf2 for some nf1,nf2)
    nf1 = minval(jf)
    nf2 = maxval(jf)
    do ifreq = nf1+1,nf2-1
       if(.not.any(jf==ifreq)) then
          write(errstr,*) "frequency indices of spectral kernel do not meet the requirements: "//&
               "all indices in range jf_min,jf_max = ",nf1,nf2," must be present, but index ",ifreq," is not"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end do
!
! now compute time waveform kernel
!
    ! initiate tkernel
    call initiateTimeWaveformKernel(tkernel,.pmtrz.skernel,param,.ntot.skernel,comp,errmsg,.onwp.skernel)
    if(.level.errmsg == 2) goto 1
    ! IN OTHER computeTimeWaveformKernel ROUTINES (if kernels are computed conventionally, i.e. analogous to computeSpectralWaveformKernel)
    ! PARAMETERS LIKE DT AND T0 ARE ALWAYS CHECKED BY THE RESPECTIVE VALUES CONTAINED IN timeKernelDisplacement, timeKernelGreenTensor ETC
    ! AND THE FIRST TIME SAMPLE FOR WHICH KERNELS ARE COMPUTED, DEFINES DF AND T0
    ! AS THIS ROUTINE IS NOT CONVENTIONAL IN THAT SENSE, THESE VALUES ARE DEFINED MANUALLY NOW!!!
    tkernel%dt = dt
    if(present(t0)) then
       tkernel%t0 = t0
    else
       tkernel%t0 = 0.
    end if
    tkernel%kd_id = kd_id
    tkernel%kgt_id = kgt_id
!
    ! open tkernel file already now to write 
    ! (simply in order to test, if file can be opened, open file is not needed before writing time kernel in loop over it below)
    call initialWriteTimeWaveformKernel(tkernel,file_tkernel,lu,errmsg,nt)
    if(.level.errmsg == 2) goto 1
!
    ! initiate discrete fourier transform
    call initiateInverseDFT(DFT,.df.skernel,nf1,nf2,tkernel%dt*jt+tkernel%t0,errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! in order to be able to do the Fourier Transformation:
    ! allocate tkernel values for ALL times (assuming that nt < nfreq)
    ! (ALTERNATIVELY you need to read and store spectral kernel values for ALL frequencies and keep those in memory, which can be more efficient if nt > nfreq)
    ! IN CASE YOU FIGURE THAT THIS NEEDS TOO MUCH MEMORY: YOU NEED TO READ ALL FREQUENCIES nt TIMES (for all times). THEN YOU NEED TO REWRITE THIS FOR TIME-WISE COMPUTATION AND WRITING)
    allocate(tkernel_kernel_all_t(nt,tkernel%ntot_kernel*ncomp*nparam))
    tkernel_kernel_all_t = 0.
!
    allocate(skernel_vals_filtered_reshaped(.ntot.skernel*nparam*ncomp))
    n_shift = .ntot.skernel*nparam
!
    ! loop on all frequencies of skernel (was checked: all indices from nf1,nf2)
    frequency_count = 0
    !write(*,*) "in transformSpectralToTimeWaveformKernel: DFT on the fly for ",.nfreq.skernel," frequencies now"
    do while(nextFrequencySpectralWaveformKernel(skernel,ifreq))
       frequency_count = frequency_count + 1
       !write(*,*) frequency_count,"'th frequency, index ",ifreq," = ",ifreq*.df.skernel," Hz"
!
       ! use a new error message for each frequency, as otherwise warning messages would get appended
       ! to errmsg (too many?!)
       call new(errmsg2,myname)
!
       ! read skernel values of that frequency from kernel file
       call readSpectralWaveformKernel(skernel,ifreq,errmsg2)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "error in readSpectralWaveformKernel reading frequency index ",ifreq
          call add(errmsg,2,errstr,myname)
          call print(errmsg2)
          call dealloc(errmsg2)
          goto 1
       end if
!
       ! reshape, and (if needed) filter the spectral kernel values at all components
       n = 0
       do icomp=1,ncomp
          skernel_vals => getValuesByCompSpectralWaveformKernel(skernel,comp(icomp)) ! must have size(.ntot.skernel,nparam) !
          if(apply_filter) then
             skernel_vals_filtered_reshaped(n+1:n+n_shift) = reshape( skernel_vals*filter(frequency_count,icomp) , (/n_shift/))
          else
             skernel_vals_filtered_reshaped(n+1:n+n_shift) = reshape( skernel_vals , (/n_shift/))
          end if
          n = n + n_shift
       end do ! icomp
!
       ! transform to time domain on the fly, immediately adding the ifreq-summand of the DFT to
       ! the values in tkernel_kernel_all_t
       call transformSpectraOnTheFlyInverseDFT(DFT,skernel_vals_filtered_reshaped,ifreq,&
            tkernel_kernel_all_t,errmsg2,add2traces=.true.)
       if(.level.errmsg2 /= 0) call print(errmsg2)
       if(.level.errmsg2 == 2) then
          !write(errstr,*) "error in transformSpectraOnTheFlyInverseDFT at frequency index ",ifreq
          call add(errmsg,2,errstr,myname)
          call print(errmsg2)
          call dealloc(errmsg2)
          goto 1
       end if
!
       call dealloc(errmsg2)
    end do ! while(nextFrequencySpectralWaveformKernel(skernel,ifreq))
!
    ! loop on all times and write time kernel file
    !write(*,*) "in transformSpectralToTimeWaveformKernel: writing time kernel for ",nt," time steps"
    do it = 1,nt
       !write(*,*) it,"'th time step, index ",jt(it)," = ",jt(it)*tkernel%dt," s"
       ! set values in tkernel%kernel and tkernel%jtcur manually
       tkernel%kernel = reshape(tkernel_kernel_all_t(it,:),(/tkernel%ntot_kernel,nparam,ncomp/))
       tkernel%jtcur = jt(it)
!
       ! write time kernel for this time step
       call new(errmsg2,myname)
       call writeTimeWaveformKernel(tkernel,errmsg2)
       if(.level.errmsg2 /= 0) call print(errmsg2)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "error in writeTimeWaveformKernel at time index ",it
          call add(errmsg,2,errstr,myname)
          call print(errmsg2)
          call dealloc(errmsg2)
          goto 1
       end if
       call dealloc(errmsg2)
    end do ! it
!
    ! close tkernel file
    call finalWriteTimeWaveformKernel(tkernel,errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! clean up
1   if(allocated(tkernel_kernel_all_t)) deallocate(tkernel_kernel_all_t)
    if(allocated(skernel_vals_filtered_reshaped)) deallocate(skernel_vals_filtered_reshaped)
    call dealloc(DFT)
    call deallocateTimeWaveformKernel(tkernel)
  end subroutine transformSpectralToTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given parameter
!
  function getValuesByCompTimeWaveformKernel(this,comp) result(p)
    type (time_waveform_kernel) :: this
    character(len=*) :: comp
    real, dimension(:,:), pointer :: p
    integer :: icomp
    nullify(p)
    if(this%kernel_stat == -1) return
    if(.not.associated(this%kernel)) return !! should be redundant with if(this%kernel_stat == -1) return
    do icomp = 1,this%ncomp
       if(this%comp(icomp) == comp) then
          p => this%kernel(:,:,icomp)
          return
       end if
    end do ! icomp
  end function getValuesByCompTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given component and parameter
!
  function getValuesByParamCompTimeWaveformKernel(this,param,comp) result(p)
    type (time_waveform_kernel) :: this
    character(len=*) :: param,comp
    real, dimension(:), pointer :: p
    integer :: jcomp,jcomp_requested,jparam,jparam_requested
!
    nullify(p)
!
    if(this%kernel_stat == -1) return
    if(.not.associated(this%kernel)) return !! should be redundant with if(this%kernel_stat == -1) return
!
    ! find index of requested parameter
    jparam_requested = -1
    do jparam=1,this%nparam
       if(this%param(jparam)==param) then
          jparam_requested = jparam
          exit
       end if
    end do ! jparam
    if(jparam_requested == -1) return
!
    ! find index of requested component
    jcomp_requested = -1
    do jcomp=1,this%ncomp
       if(this%comp(jcomp)==comp) then
          jcomp_requested = jcomp
          exit
       end if
    end do ! jcomp
    if(jcomp_requested == -1) return
!    
    p => this%kernel(:,jparam_requested,jcomp_requested)
  end function getValuesByParamCompTimeWaveformKernel
! !------------------------------------------------------------------------
! !> \brief return kernel values for given index of parameter
! !
!   function getValuesByIndexTimeWaveformKernel(this,i) result(p) !! ROUTINE IS DEPRECIATED (would need to refer by index to the current sorting and content of vector this%param, so it's safer to use above routine getValuesByParamSpectralWaveformKernel)
!     type (time_waveform_kernel) :: this
!     integer :: i
!     real, dimension(:,:), pointer :: p
!     nullify(p)
!     if(.not.associated(this%kernel)) return
!     if(i<1 .or. i>this%nparam) return
!     p => this%kernel(:,:,i)
!   end function getValuesByIndexTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to all kernel values
!
  function getAllValuesTimeWaveformKernel(this) result(p)
    type (time_waveform_kernel) :: this
    real, dimension(:,:,:), pointer :: p
    p => this%kernel
  end function getAllValuesTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return parametrization of kernel
!
  function getParametrizationTimeWaveformKernel(this) result(pmtrz)
    type (time_waveform_kernel), intent(in) :: this
    character(len=character_length_pmtrz) :: pmtrz
    pmtrz = this%parametrization
  end function getParametrizationTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of inversion grid cells of kernel
!
  function getNtotKernelTimeWaveformKernel(this) result(ntot)
    type (time_waveform_kernel), intent(in) :: this
    integer :: ntot
    ntot = this%ntot_kernel
  end function getNtotKernelTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of components of kernel
!
  function getNcompTimeWaveformKernel(this) result(ncomp)
    type (time_waveform_kernel), intent(in) :: this
    integer :: ncomp
    ncomp = this%ncomp
  end function getNcompTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of components of kernel
!
  function getCompTimeWaveformKernel(this) result(comp)
    type (time_waveform_kernel), intent(in) :: this
    character(len=character_length_component), dimension(:), pointer :: comp
    comp => this%comp
  end function getCompTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return t0 of kernel
!
  function getT0TimeWaveformKernel(this) result(t0)
    type (time_waveform_kernel), intent(in) :: this
    real :: t0
    t0 = this%t0
  end function getT0TimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return dt of kernel
!
  function getDtTimeWaveformKernel(this) result(dt)
    type (time_waveform_kernel), intent(in) :: this
    real :: dt
    dt = this%dt
  end function getDtTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of time samples of kernel
!
  function getNtTimeWaveformKernel(this) result(nt)
    type (time_waveform_kernel), intent(in) :: this
    integer :: nt
    nt = this%nt
  end function getNtTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to array of time indices of kernel
!
  function getJtTimeWaveformKernel(this) result(jt)
    type (time_waveform_kernel), intent(in) :: this
    integer, dimension(:), pointer :: jt
    jt => this%jt
  end function getJtTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief open stream access file to write kernel values, one time step at a time
!! \param 
!
  subroutine initialWriteTimeWaveformKernel(this,filename,lu,errmsg,nt)
    type (time_waveform_kernel) :: this
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    integer, optional :: nt
    ! local
    integer :: ios
    character(len=400) :: errstr
    character(len=30) :: myname = 'initialWriteTimeWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    if(this%kernel_stat == -1) then
       call add(errmsg,2,"time_waveform_kernel object not initiated yet",myname)
       return
    end if
!
    if(this%kernel_stat > 0) then
       write(errstr,*) "this waveform kernel object is already opened to read or to write. ",&
            "please call final read/write or deallocate first"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ios = createFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to write, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%kernel_stat = 1
!
    if(present(nt)) then
       ! it is more efficient to create root like this, knowing the number of subgroups to come
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=nt,maxdset=4)
    else
       ! otherwise assume enough subgroups, which is still more efficient than creating root without subrgoup information
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=150,maxdset=4)
    end if
  end subroutine initialWriteTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief write current values of waveform kernel to file
!! \details It is important that the time samples are written one after another, 
!!  i.e. first time index 1, then index 2, then ... etc., because the time samples
!!  are just appended to the stream access file and when reading the file again, 
!!  it is assumed that the time samples were written to file in time index order.
!! \param 
!
  subroutine writeTimeWaveformKernel(this,errmsg)
    type (time_waveform_kernel) :: this
    type (error_message) :: errmsg
    ! local
    integer :: iparam,icomp,maxdataset
    type (group_stream_access) :: time
    type (data_stream_access) :: dset
    character(len=23) :: myname = 'writeTimeWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    ! if file is not open to write, return
    if(this%kernel_stat /= 1) then
       call add(errmsg,2,"file not yet opened to write, call initialWriteTimeWaveformKernel first",myname)
       return
    end if
!
    ! if there is an indicator that there are no sensible kernel values stored, just return doing nothing
    if(this%jtcur == -1) then
       call add(errmsg,2,"no sensible kernel values were computed, cannot write any values",myname)
       return
    end if
!
    maxdataset = this%ncomp*this%nparam
    call createGroupStreamAccess(time,'Time',this%jtcur,maxsubgroup=0,maxdset=maxdataset)
    do icomp = 1,this%ncomp
       do iparam = 1,this%nparam
          call new(dset,1,(/ this%ntot_kernel /),T_REAL)
          call writeDatasetVectorStreamAccess(dset,this%fsa,this%kernel(:,iparam,icomp))
          call addDatasetStreamAccess(time,dset)
          call dealloc(dset)
       end do ! iparam
    end do ! icomp
    call addSubgroupStreamAccess(this%root,time)
    call dealloc(time)
!
    ! remember all time samples written to file (in their order) in array this%jt by appending the current time to array this%jt
    this%jt => reallocate(this%jt,this%nt+1)
    this%nt = this%nt + 1
    this%jt(this%nt) = this%jtcur
  end subroutine writeTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief close stream access file
!! \param this time waveform kernel
!! \param errmsg error message
!! \param lu return file unit of this waveform kernel file, it can be added to file unit handler
!
  subroutine finalWriteTimeWaveformKernel(this,errmsg,lu)
    type (time_waveform_kernel) :: this
    type (error_message) :: errmsg
    integer, optional :: lu
    ! local
    character(len=28) :: myname = 'finalWriteTimeWaveformKernel'
    type (data_stream_access) :: dset
    type (flexible), dimension(:), allocatable :: ft 
    integer :: iparam,icomp,on_invgrid_int
!
    call addTrace(errmsg,myname)
!
    ! if file is not open to write, return
    if(this%kernel_stat /= 1) then
       call add(errmsg,2,"file not yet opened to write, call initialWriteSpectralWaveformKernel first",myname)
       return
    end if
!
    ! write general information to first dataset in group root
    allocate(ft(10))
    call new(dset,1,(/ 10 /),T_FLEXIBLE)
    ft(1) = this%parametrization
    ft(2) = this%nparam
    ft(3) = this%ntot_kernel
    ft(4) = this%dt
    ft(5) = this%nt
    ft(6) = this%t0
    ft(7) = this%kd_id
    ft(8) = this%kgt_id
    ft(9) = this%ncomp
    if(this%on_invgrid) then
       on_invgrid_int = 1
    else
       on_invgrid_int = 0
    end if
    ft(10) = on_invgrid_int
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write parameter names to second dataset in group root
    allocate(ft(this%nparam))
    call new(dset,1,(/this%nparam/),T_FLEXIBLE)
    ! use a do loop in order to assure that character lengths are properly defined
    ! in any case, ft = this%comp is not defined (since "=" is an overlayn operator for type_flexible)
    do iparam = 1,this%nparam
       ft(iparam) = this%param(iparam)
    end do ! iparam
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write components to third dataset in group root
    allocate(ft(this%ncomp))
    call new(dset,1,(/ this%ncomp /),T_FLEXIBLE)
    ! use a do loop in order to assure that character lengths are properly defined
    ! in any case, ft = this%comp is not defined (since "=" is an overlayn operator for type_flexible)
    do icomp = 1,this%ncomp
       ft(icomp) = this%comp(icomp)
    end do ! icomp
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write time indices to fourth dataset in group root
    if(this%nt > 0) then
       call new(dset,1,(/ this%nt /),T_INTEGER)
       call writeDatasetVectorStreamAccess(dset,this%fsa,this%jt)
       call addDatasetStreamAccess(this%root,dset)
       call dealloc(dset)
    else
       call add(errmsg,1,"No kernel values have been written to file! "//&
            "File only contains meta information",myname)
    end if
!
    ! write file content information
    call writeGroupStreamAccess(this%root,this%fsa)
!
    ! clean up
    call clearGroupTree(this%root)
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call dealloc(this%fsa)
    this%kernel_stat = 0
  end subroutine finalWriteTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief get meta information from a kernel file
!
  subroutine readMetaInfoTimeWaveformKernel(filename,lu,errmsg,parametrization,nparam,ntot_kernel,&
         dt,t0,nt,kd_id,kgt_id,ncomp,on_invgrid,param,comp,jt)
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    character(len=character_length_pmtrz), optional :: parametrization
    integer, optional :: nparam,ntot_kernel,nt,ncomp
    real, optional :: dt,t0
    character(len=length_ID_kernel_displacement), optional :: kd_id
    character(len=length_ID_kernel_green_tensor), optional :: kgt_id
    logical, optional :: on_invgrid
    character(len=character_length_param), dimension(:), pointer, optional :: param
    character(len=character_length_component), dimension(:), pointer, optional :: comp
    integer, dimension(:), pointer, optional :: jt
    ! local
    character(len=400) :: errstr
    character(len=30) :: myname = 'readMetaInfoTimeWaveformKernel'
    type (file_stream_access) :: fsa
    type (group_stream_access) :: root
    integer :: ios
    character(len=character_length_pmtrz) :: parametrization_loc
    integer :: nparam_loc,ntot_kernel_loc,nt_loc,ncomp_loc
    real :: dt_loc,t0_loc
    character(len=length_ID_kernel_displacement) :: kd_id_loc
    character(len=length_ID_kernel_green_tensor) :: kgt_id_loc
    logical :: on_invgrid_loc
    character(len=character_length_param), dimension(:), pointer :: param_loc
    character(len=character_length_component), dimension(:), pointer :: comp_loc
    integer, dimension(:), pointer :: jt_loc
!
    nullify(param_loc,comp_loc,jt_loc)
!
    call addTrace(errmsg,myname)
    if(present(param)) nullify(param)
    if(present(comp)) nullify(comp)
    if(present(jt)) nullify(jt)
!
    ! open stream access file
    ios = openFileStreamAccess(fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to read, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    else
       write(errstr,*) "successfully opened file '"//trim(filename)//"' to read; iostat = ",ios
       call add(errmsg,0,errstr,myname)
   end if
!
    ! read content information of stream access file
    call readGroupStreamAccess(root,fsa)
!
    call privateReadMetaInfoTimeWaveformKernel(fsa,root,parametrization_loc,nparam_loc,ntot_kernel_loc,&
         dt_loc,t0_loc,nt_loc,kd_id_loc,kgt_id_loc,ncomp_loc,on_invgrid_loc,param_loc,comp_loc,jt_loc,errmsg,myname)
    if(.level.errmsg == 2) goto 1
!
    if(present(parametrization)) parametrization = parametrization_loc
    if(present(nparam)) nparam = nparam_loc
    if(present(ntot_kernel)) ntot_kernel = ntot_kernel_loc
    if(present(dt)) dt = dt_loc
    if(present(t0)) t0 = t0_loc
    if(present(nt)) nt = nt_loc
    if(present(kd_id)) kd_id = kd_id_loc
    if(present(kgt_id)) kgt_id = kgt_id_loc
    if(present(ncomp)) ncomp = ncomp_loc
    if(present(on_invgrid)) on_invgrid = on_invgrid_loc
    if(present(param)) then
       param => param_loc
       nullify(param_loc)
    end if
    if(present(comp)) then
       comp => comp_loc
       nullify(comp_loc)
    end if
    if(present(jt)) then
       jt => jt_loc
       nullify(jt_loc)
    end if
!
1   if(associated(param_loc)) deallocate(param_loc)
    if(associated(comp_loc)) deallocate(comp_loc)
    if(associated(jt_loc)) deallocate(jt_loc)
    call clearGroupTree(root)
    call dealloc(fsa)
  end subroutine readMetaInfoTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief private routine to read all meta information from opened kernel file
!
  subroutine privateReadMetaInfoTimeWaveformKernel(fsa,root,parametrization,nparam_in_file,ntot_kernel,&
         dt,t0,nt,kd_id,kgt_id,ncomp_in_file,on_invgrid,param_in_file,comp_in_file,jt_in_file,errmsg,myname)
    type (file_stream_access) :: fsa
    type (group_stream_access) :: root
    character(len=character_length_pmtrz) :: parametrization
    integer:: nparam_in_file,ntot_kernel,nt,ncomp_in_file
    real :: dt,t0
    character(len=length_ID_kernel_displacement) :: kd_id
    character(len=length_ID_kernel_green_tensor) :: kgt_id
    logical :: on_invgrid
    character(len=character_length_param), dimension(:), pointer :: param_in_file
    character(len=character_length_component), dimension(:), pointer :: comp_in_file
    integer, dimension(:), pointer :: jt_in_file
    type (error_message) :: errmsg
    character(len=*) :: myname
    ! local
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    type (flexible), dimension(:), pointer :: ft
    integer :: on_invgrid_int,iparam_in_file,icomp_in_file,ifreq
    character(len=400) :: errstr
!
    nullify(group,dset,ft)
!
    nullify(param_in_file,comp_in_file,jt_in_file)
!
    ! read first info vector
    call traversePathStreamAccess(root,0,(/ 1 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fsa,ft)
    if(.not.associated(ft)) then
       write(errstr,*) "there is no first information vector. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(ft)/=10) then
       write(errstr,*) "size of first information vector is ",size(ft),", must be 10. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    parametrization = ft(1)
    nparam_in_file = ft(2)
    ntot_kernel = ft(3)
    dt = ft(4)
    nt = ft(5)
    t0 = ft(6)
    kd_id = ft(7)
    kgt_id = ft(8)
    ncomp_in_file = ft(9)
    on_invgrid_int = ft(10)
    deallocate(ft)
!    
    if(nparam_in_file <= 0) then
       write(errstr,*) "number of elastic parameters as of information vector = ",nparam_in_file,&
            "; must be > 0, file is inconsistent"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    if(ncomp_in_file <= 0) then
       write(errstr,*) "number of receiver components as of information vector = ",ncomp_in_file,&
            "; must be > 0, file is inconsistent"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    if(nt < 0) then
       write(errstr,*) "number of time samples as of information vector = ",nt,"; must be >= 0 , file is incosistent"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    if(on_invgrid_int /= 0 .and. on_invgrid_int /= 1) then
       write(errstr,*) "inconsistency: invgrid/wp flag is ",on_invgrid_int,&
            "; must be either 1 or 0"
       call add(errmsg,2,errstr,myname)
       return
    end if
    on_invgrid = on_invgrid_int == 1
!
    ! read second info vector containing the elastic parameter names
    call traversePathStreamAccess(root,0,(/ 2 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fsa,ft)
    if(.not.associated(ft)) then
       write(errstr,*) "there is no vector of elastic parameter names. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(ft)/=nparam_in_file) then
       write(errstr,*) "information on parameter names inconsistent: number of parameters ",size(ft),&
            " differs from expected number ",nparam_in_file
       call add(errmsg,2,errstr,myname)
       return
    end if
    allocate(param_in_file(nparam_in_file))
    do iparam_in_file = 1,nparam_in_file
       param_in_file(iparam_in_file) = ft(iparam_in_file)
    end do ! iparam_in_file
    deallocate(ft)
!
    ! read third info vector containing the receiver component names
    call traversePathStreamAccess(root,0,(/ 3 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fsa,ft)
    if(.not.associated(ft)) then
       write(errstr,*) "there is no vector of receiver components. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(ft)/=ncomp_in_file) then
       write(errstr,*) "information on components inconsistent: number of components ",size(ft),&
            " differs from expected number ",ncomp_in_file
       call add(errmsg,2,errstr,myname)
       return
    end if
    allocate(comp_in_file(ncomp_in_file))
    do icomp_in_file = 1,ncomp_in_file
       comp_in_file(icomp_in_file) = ft(icomp_in_file)
    end do ! icomp_in_file
    deallocate(ft)
!
    ! read fourth info vector containing the time indices (if it was written at all)
    ! if nt == 0, it means that there is no data contained in the file, only meta-information contained in file
    ! this module actually permitts to produce such files (although they might by somewhat useless)
    if(nt>0) then
       call traversePathStreamAccess(root,0,(/ 4 /),group,dset)
       call readDatasetVectorStreamAccess(dset,fsa,jt_in_file)
       if(.not.associated(jt_in_file)) then
          write(errstr,*) "there is no vector of time indices. there is a problem with the file"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(jt_in_file)/=nt) then
          write(errstr,*) "information on time steps inconsistent: number of time indices ",size(jt_in_file),&
               " differs from expected number ",nt,"; dt= ",dt!,", time indices = ",jt_in_file !! this could cause errstr to overflow (too many characters)
          write(*,*) "ERROR in ",trim(myname),": time indices = ",jt_in_file
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(any(jt_in_file<0)) then
          write(errstr,*) "there are ",count(jt_in_file<0)," time indices in the kernel file which are < 0"
          call add(errmsg,2,errstr,myname)
          return
       end if
       do ifreq = 2,nt
          if(any(jt_in_file(1:ifreq-1)==jt_in_file(ifreq))) then
             write(errstr,*) ifreq,"'th time index contained in kernel file occurs more than once. ",&
                  "There must not be duplicated time indices"
             call add(errmsg,2,errstr,myname)
             return
          end if
       end do
    end if ! nt>0
!
  end subroutine privateReadMetaInfoTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief open stream access file to read kernel values, one time step at a time
!! \param 
!
  subroutine initialReadTimeWaveformKernel(this,filename,lu,errmsg)
    type (time_waveform_kernel) :: this
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    integer :: ios
    character(len=character_length_pmtrz) :: parametrization
    character(len=length_ID_kernel_displacement) :: kd_id
    character(len=length_ID_kernel_green_tensor) :: kgt_id
    integer :: nparam_in_file,ntot_kernel,ncomp_in_file,nt,iparam,iparam_in_file,icomp,icomp_in_file
    real :: dt,t0
    logical :: on_invgrid
    integer, dimension(:), pointer :: jt_in_file,comp_indx_in_file,param_indx_in_file
    character(len=character_length_param), dimension(:), pointer :: param_in_file
    character(len=character_length_component), dimension(:), pointer :: comp_in_file
    character(len=400) :: errstr
    character(len=29) :: myname = 'initialReadTimeWaveformKernel'
!
    nullify(jt_in_file,comp_indx_in_file,param_indx_in_file,param_in_file,comp_in_file)
!
    call addTrace(errmsg,myname)
!
    if(this%kernel_stat == -1) then
       call add(errmsg,2,"spectral_waveform_kernel object not initiated yet",myname)
       return
    end if
!
    if(this%kernel_stat > 0) then
       write(errstr,*) "this waveform kernel object is already opened to read or to write. ",&
            "please call final read/write or deallocate first"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! open stream access file
    ios = openFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to read, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    else
       write(errstr,*) "successfully opened file '"//trim(filename)//"' to read; iostat = ",ios
       call add(errmsg,0,errstr,myname)
    end if
!
    ! read content information of stream access file
    call readGroupStreamAccess(this%root,this%fsa)
!
    call privateReadMetaInfoTimeWaveformKernel(this%fsa,this%root,parametrization,nparam_in_file,ntot_kernel,&
         dt,t0,nt,kd_id,kgt_id,ncomp_in_file,on_invgrid,param_in_file,comp_in_file,jt_in_file,errmsg,myname)
    if(.level.errmsg == 2) goto 1
!
    if(parametrization /= this%parametrization) then
       write(errstr,*) "the kernel in this file is of model parameterization '",trim(parametrization),&
            "', which does not equal the parametrization '",trim(this%parametrization),&
            "' for which this object was initiated."
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    if(.not.(on_invgrid .eqv. this%on_invgrid)) then
       if(on_invgrid) then
          write(errstr,*) "inconsistency: values in kernel file are on inversion grid cells, but initially ",&
               "values on wavefield points were requested for this object"
       else
          write(errstr,*) "inconsistency: values in kernel file are on wavefield_points, but initially ",&
               "values on inversion grid cells were requested for this object"
       end if
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    if(ntot_kernel /= this%ntot_kernel) then
       if(on_invgrid) then
          write(errstr,*) "number of inversion grid cells contained in kernel file = ",ntot_kernel,&
               " differs from number of inversion grid cells for which this object was initiated = ",&
               this%ntot_kernel
       else
          write(errstr,*) "number of wavefield points contained in kernel file = ",ntot_kernel,&
               " differs from number of wavefield points for which this object was initiated = ",&
               this%ntot_kernel
       end if
    end if
!
    ! check if all parameter names for which this kernel object was initiated (i.e. which are requested to be
    ! read from the kernel file) are also actually contained in the kernel file.
    ! memorize which parameters should be read from file (and their order and position within the file)
    allocate(param_indx_in_file(this%nparam))
    param_indx_in_file = -1
    do iparam = 1,this%nparam
       do iparam_in_file = 1,nparam_in_file
          if(this%param(iparam)==param_in_file(iparam_in_file)) then
             ! store index in vector param (using index offset 1) which can be used below for traversePathStreamAccess
             param_indx_in_file(iparam) = iparam_in_file
             exit ! if an index is found, exit the loop on iparam_in_file
          end if
       end do ! iparam_in_file
    end do ! iparam
    if(any(param_indx_in_file == -1)) then
       do iparam = 1,this%nparam
          if(param_indx_in_file(iparam)==-1) then
             write(errstr,*) "model parameter '",trim(this%param(iparam)),"' was not found in the file"
             call add(errmsg,2,errstr,myname)
          end if
       end do
       write(errstr,*) count(param_indx_in_file == -1)," out of ",this%nparam,&
            " model parameters for which this object was initiated were not found in the file"
       write(*,*) "ERROR in initialReadTimeWaveformKernel: model parameters for which this object was initiated = ",&
            this%param,"; model parameters contained in kernel file = ",param_in_file
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    deallocate(param_in_file)
!
    ! check if all components for which this kernel object was initiated (i.e. which are requested to be
    ! read from the kernel file) are also actually contained in the kernel file.
    ! memorize which components should be read from file (and their order and position within the file)
    allocate(comp_indx_in_file(this%ncomp))
    comp_indx_in_file = -1
    do icomp = 1,this%ncomp
       do icomp_in_file = 1,ncomp_in_file
          if(this%comp(icomp)==comp_in_file(icomp_in_file)) then
             ! store index in vector comp (using index offset 1) which can be used below for traversePathStreamAccess
             comp_indx_in_file(icomp) = icomp_in_file
             exit ! if an index is found, exit the loop on iparam_in_file
          end if
       end do ! icomp_in_file
    end do ! icomp
    if(any(comp_indx_in_file == -1)) then
       do icomp = 1,this%ncomp
          if(comp_indx_in_file(icomp)==-1) then
             write(errstr,*) "receiver component '",trim(this%comp(icomp)),"' was not found in the file"
             call add(errmsg,2,errstr,myname)
          end if
       end do
       write(errstr,*) count(comp_indx_in_file == -1)," out of ",this%ncomp,&
            " receiver components for which this object was initiated were not found in the file"
       write(*,*) "ERROR in initialReadTimeWaveformKernel: receiver components for which this object was initiated = ",&
            this%comp,"; receiver components contained in kernel file = ",comp_in_file
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    deallocate(comp_in_file)
!
    ! if everything went alright, return here after initiating everything and setting flag this%kernel_stat = 2 (file is open to read)
    if(nt>0) then
       ! initiate object with the correct frequency content (this is not done in initiateSpectralWaveformKernel above)
       this%t0 = t0
       this%dt = dt
       this%nt = nt
       this%jt => jt_in_file
    end if
    ! initiate id's of source and receiver, to which the kernel values are associated
    ! (this is not done in initiateSpectralWaveformKernel above)
    this%kd_id = kd_id
    this%kgt_id = kgt_id
    this%ncomp_in_file = ncomp_in_file
    this%comp_indx_in_file => comp_indx_in_file
    this%nparam_in_file = nparam_in_file
    this%param_indx_in_file => param_indx_in_file
    this%kernel_stat = 2
    return
!
    ! if there went anything wrong in intial read, close file here and revert all initiation that's been done so far
1   if(associated(comp_indx_in_file)) deallocate(comp_indx_in_file)
    if(associated(param_indx_in_file)) deallocate(param_indx_in_file)
    if(associated(jt_in_file)) deallocate(jt_in_file)
    if(associated(comp_in_file)) deallocate(comp_in_file)
    if(associated(param_in_file)) deallocate(param_in_file)
    call finalReadTimeWaveformKernel(this)
  end subroutine initialReadTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief read values of time waveform kernel for some specific time (index jt)
!! \param 
!
  subroutine readTimeWaveformKernel(this,jt,errmsg)
    type (time_waveform_kernel) :: this
    integer :: jt
    type (error_message) :: errmsg
    ! local
    type (group_stream_access), pointer :: group,tgroup
    type (data_stream_access), pointer :: dset
    real, dimension(:), pointer :: d
    integer :: it,iparam,iparam_in_file,icomp,icomp_in_file,igroup,idset
    logical :: traverse_path_from_root
    character(len=400) :: errstr
    character(len=22) :: myname = 'readTimeWaveformKernel'
!
    nullify(group,tgroup,dset,d)
!
    call addTrace(errmsg,myname)
!
    if(this%kernel_stat /= 2) then
       call add(errmsg,2,"this waveform kernel object was not yet initiated to read. "//&
            "call initialReadTimeWaveformKernel first",myname)
       return
    end if
!
    igroup = 0
    do it = 1,this%nt
       if(this%jt(it) == jt) then
          igroup = it
          exit
       end if
    end do
    if(igroup == 0) then
       write(errstr,*) "incoming time index ",jt,", not contained in file. times contained in file are: dt,t0 = ",&
            this%dt,this%t0!,"; time indices = ",this%jt !! this could cause errstr to overflow (too many characters)
       write(*,*) "ERROR in readTimeWaveformKernel: time indices = ",this%jt
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    traverse_path_from_root = .true.
    do icomp = 1,this%ncomp
       icomp_in_file = this%comp_indx_in_file(icomp)
       do iparam = 1,this%nparam
          iparam_in_file = this%param_indx_in_file(iparam)
          idset = (icomp_in_file-1)*this%nparam_in_file + iparam_in_file
          if(traverse_path_from_root) then
             ! remember time index group to make calls traversePathStreamAccess faster below
             call traversePathStreamAccess(this%root,1,(/ igroup,idset /),tgroup,dset)
             ! only do this in the very first iteration of this loops, so set flag to false
             traverse_path_from_root = .false.
          else ! traverse_path_from_root
             call traversePathStreamAccess(tgroup,0,(/ idset /),group,dset)
          end if ! traverse_path_from_root
          call readDatasetVectorStreamAccess(dset,this%fsa,d)
          if(.not.associated(d)) then
             write(errstr,*) "there is no dataset in stream_access file for component ",icomp," of ",iparam,&
                  "'th parameter '",trim(getParamFromIndexModelParametrization(this%parametrization,iparam)),&
                  "'; kernel file may be corrupt"
             call add(errmsg,2,errstr,myname)
             return
          end if
          if(size(d) /= this%ntot_kernel) then
             write(errstr,*) "number of kernel values ",size(d)," for component ",icomp," = '"//trim(this%comp(icomp))//&
                  "' of ",iparam,&
                  "'th parameter '",trim(getParamFromIndexModelParametrization(this%parametrization,iparam)),&
                  "' differs from number of inversion grid cells ",this%ntot_kernel,", i.e. kernel file may be corrupt"
             call add(errmsg,2,errstr,myname)
             return
          end if
          this%kernel(:,iparam,icomp) = d
          deallocate(d)
       end do ! iparam
    end do ! icomp
!
    this%jtcur = jt
  end subroutine readTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief close stream access file
!! \param this time waveform kernel
!! \param lu return file unit of this waveform kernel file, it can be added to file unit handler
!
  subroutine finalReadTimeWaveformKernel(this,lu)
    type (time_waveform_kernel) :: this    
    integer, optional :: lu
    call clearGroupTree(this%root)
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call dealloc(this%fsa)
    this%kernel_stat = 0
  end subroutine finalReadTimeWaveformKernel
!
end module timeWaveformKernel
