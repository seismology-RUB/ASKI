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
!> \brief module which computes integrated time waveform sensitivity kernels
!!
!! \details By inversely Fourier transforming spectral waveform sensitivity kernels
!!  (or by direct computation from time domain wavefields?! -> not supported yet), time domain
!!  waveform sensitivity kernels are computed or written/read to/from file by this module.
!!
!! \author Florian Schumacher
!! \date August 2013
!
module timeWaveformKernel
!
  use modelParametrization
  use componentTransformation
  use spectralWaveformKernel
  use discreteFourierTransform
  use flexibleType
  use streamAccess
  use realloc
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateTimeWaveformKernel; end interface
  interface getValuesTimeWaveformKernel
     module procedure getValuesByParamTimeWaveformKernel
     module procedure getValuesByCompParamTimeWaveformKernel
     module procedure getValuesByIndexTimeWaveformKernel
     module procedure getAllValuesTimeWaveformKernel
  end interface getValuesTimeWaveformKernel
  interface operator (.pmtrz.); module procedure getParametrizationTimeWaveformKernel; end interface
  interface operator (.ntot.); module procedure getNtotInvgridTimeWaveformKernel; end interface
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
     integer :: nparam = 0 !< nparam = numberOfParamModelParametrization(parametrization), size of last dimension of kernel
     integer :: ntot_invgrid = 0 !< total number of inversion grid cells on which kernels are given
     ! station component
     integer :: ncomp = 0 !< number of components for which kernel values are contained in array kernel (2nd dimension)
     character(len=character_length_component), dimension(:), pointer :: comp => null() !< the ncomp component names
     ! temporal discretization
     real :: t0 = 0. !< global time shift of times defined by jt*dt, IS SUPPOSED TO BE USED (only), IF INDICES jt GET TOO LARGE FOR INTEGER(4)!!
     real :: dt = -1. !< time step of time discretization
     integer :: nt = 0 !< total number of time samples
     integer, dimension(:), pointer :: jt => null() !< time indices in order of subgroups of kernel file; times compute as t = t0+jt*dt [s]
     integer :: jtcur = -1 !< indicates current time index jt (i.e. time t=t0+jtcur*dt) for which kernel values are contained in array kernel
     ! kernel values
     real, dimension(:,:,:), pointer :: kernel => null() !< kernel values for one frequency dim(ntot_invgrid,ncomp,nparam)
     ! file handling
     integer :: filestat = 0 !< indicating closed (0), open to write (1), open to read (2)
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
  subroutine initiateTimeWaveformKernel(this,parametrization,ntot_invgrid,comp,errmsg)
    type (time_waveform_kernel) :: this
    character(len=*) :: parametrization
    integer :: ntot_invgrid
    character(len=*), dimension(:) :: comp
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=26) :: myname = 'initiateTimeWaveformKernel'
    integer :: icomp,ncomp
!
    call addTrace(errmsg,myname)
!
    if(.not.validModelParametrization(parametrization)) then
       call add(errmsg,2,"incoming model parametrization '"//trim(parametrization)//&
            "' is not valid, valid parametrizations(parameters) are "//all_valid_pmtrz_param,myname)
       return
    end if
    if(ntot_invgrid <= 0) then
       write(errstr,*) "incoming number of inversion grid cells ( = ",ntot_invgrid,") must be positive!"
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
    do icomp = 1,ncomp
       if(.not.validComponent(comp(icomp))) then
          write(errstr,*) icomp,"'th incoming station component '"//trim(comp(icomp))//"' (out of ",&
               ncomp,")is invalid"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end do ! icomp
!
    if(associated(this%kernel)) then
       call add(errmsg,1,"time waveform kernel is already initiated, deallocating it now before initiating new one",myname)
       call deallocateTimeWaveformKernel(this)
    end if
!
    this%parametrization = parametrization
    this%nparam = numberOfParamModelParametrization(this%parametrization)
    this%ntot_invgrid = ntot_invgrid
    this%ncomp = ncomp
    allocate(this%comp(this%ncomp))
    this%comp = comp
    this%jtcur = -1
!
    allocate(this%kernel(this%ntot_invgrid,3,this%nparam))
!
    write(errstr,*) "initiated time waveform kernel for ",this%ntot_invgrid," inversion grid cells, ",this%nparam,&
         " parameters of parametrization '",trim(this%parametrization),"' and ",ncomp," station components"
    call add(errmsg,0,errstr,myname)
  end subroutine initiateTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this time waveform kernel
!
  subroutine deallocateTimeWaveformKernel(this)
    type (time_waveform_kernel) :: this
    this%parametrization = ''
    this%nparam = 0
    this%ntot_invgrid = 0
    this%ncomp = 0
    if(associated(this%comp)) deallocate(this%comp)
    this%dt = -1.
    this%nt = 0
    if(associated(this%jt)) deallocate(this%jt)
    this%jtcur = -1
    if(associated(this%kernel)) deallocate(this%kernel)
    if(this%filestat/=0) then
       call clearGroupTree(this%root)
       call dealloc(this%fsa)
       this%filestat = 0
    end if
  end subroutine deallocateTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief compute time domain kernel from spectral kernel and write to file for all times
!! \details Using module discreteFourierTransform and the spectral waveform kernel
!!  
!! \param 
!
  subroutine transformSpectralToTimeWaveformKernel(skernel,comp,comptrans,staname,&
       dt,jt,file_tkernel,lu,errmsg,filter,t0)
    type (spectral_waveform_kernel) :: skernel
    character(len=*), dimension(:) :: comp
    type (component_transformation) :: comptrans
    character(len=*) :: staname
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
    type (time_waveform_kernel) :: tkernel
    type (discrete_fourier_transform) :: DFT
    integer, dimension(:), pointer :: jf
    integer :: icomp,ncomp,it,nt,nf1,nf2,ifreq,frequency_count,ipar,n,n_shift
    real, dimension(:,:), allocatable :: tkernel_kernel_all_t
    double precision, dimension(:,:), pointer :: trans_coef
    complex, dimension(:,:,:), pointer :: skernel_vals
    complex, dimension(:), allocatable :: skernel_vals_trans
    double complex, dimension(:,:), allocatable :: vals_tmp
    logical :: apply_filter
!
    call addTrace(errmsg,myname)
!
! first check incoming values and sizes of arrays
!
    ncomp = size(comp)
    if(ncomp <= 0) then
       write(errstr,*) "incoming number of station components ",ncomp," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
    do icomp = 1,ncomp
       if(.not.validComponent(comp(icomp))) then
          write(errstr,*) icomp,"'th incoming station component '"//trim(comp(icomp))//"' (out of ",&
               ncomp,")is invalid"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end do ! icomp
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
               " does not match number of incoming station components ",ncomp
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
    ! get all transformation coefficients for station staname
    ! interchange C_in and C_out here, as we below need the transposed transformation matrix
    trans_coef => transform(comptrans,comp,(/'CX','CY','CZ'/),staname)
    if(.not.associated(trans_coef)) then
       write(errstr,*) "'"//comp//"', "
       call add(errmsg,2,"there are no transformation coefficients for station '"//&
            trim(staname)//"' and components "//trim(errstr),myname)
       goto 1
    end if
!
    ! initiate tkernel
    call initiateTimeWaveformKernel(tkernel,.pmtrz.skernel,.ntot.skernel,comp,errmsg)
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
    ! allocate tkernel values for ALL times (otherwise you need to read ALL spectral values nt times, this is inefficient. 
    ! IN CASE YOU FIGURE THAT THIS NEEDS TOO MUCH MEMORY: YOU NEED TO REWRITE THIS FOR TIME-WISE COMPUTATION AND WRITING)
    allocate(tkernel_kernel_all_t(nt,tkernel%ntot_invgrid*ncomp*tkernel%nparam))
    tkernel_kernel_all_t = 0.
!
    allocate(skernel_vals_trans(.ntot.skernel*ncomp*numberOfParamModelParametrization(.pmtrz.skernel)),&
         vals_tmp(.ntot.skernel,ncomp))
    n_shift = .ntot.skernel*ncomp
!
    ! loop on all frequencies of skernel (was checked: all indices from nf1,nf2)
    frequency_count = 0
    write(*,*) "in transformSpectralToTimeWaveformKernel: DFT on the fly for ",.nfreq.skernel," frequencies now"
    do while(nextFrequencySpectralWaveformKernel(skernel,ifreq))
       frequency_count = frequency_count + 1
       write(*,*) frequency_count,"'th frequency, index ",ifreq," = ",ifreq*.df.skernel," Hz"
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
       skernel_vals => getValuesSpectralWaveformKernel(skernel)
!
       ! apply component transform, reshaping, and (if needed) filtering to spectral kernel values
       n = 0
       do ipar = 1,numberOfParamModelParametrization(.pmtrz.skernel)
          ! transform (ntot_invgrid,3)-matrix to (ntot_invgrid,ncomp)-matrix for this parameter
          vals_tmp = matmul(skernel_vals(:,:,ipar),trans_coef)
          ! apply filter values of current frequency to columns of vals_tmp (filters are component dependent!)
          if(apply_filter) then
             do icomp=1,ncomp
                vals_tmp(:,icomp) = vals_tmp(:,icomp) * filter(frequency_count,icomp)
             end do ! icomp
          end if ! apply_filter
          skernel_vals_trans(n+1:n+n_shift) = reshape( vals_tmp , (/n_shift/))
          n = n + n_shift
       end do ! ipar
!
       ! transform to time domain on the fly, immediately adding the ifreq-summand of the DFT to
       ! the values in tkernel_kernel_all_t
       call transformSpectraOnTheFlyInverseDFT(DFT,skernel_vals_trans,ifreq,&
            tkernel_kernel_all_t,errmsg2,add2traces=.true.)
       if(.level.errmsg2 /= 0) call print(errmsg2)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "error in transformSpectraOnTheFlyInverseDFT at frequency index ",ifreq
          call add(errmsg,2,errstr,myname)
          call dealloc(errmsg2)
          goto 1
       end if
!
       call dealloc(errmsg2)
    end do ! while(nextFrequencySpectralWaveformKernel(skernel,ifreq))
!
    ! loop on all times and write time kernel file
    write(*,*) "in transformSpectralToTimeWaveformKernel: writing time kernel for ",nt," time steps"
    do it = 1,nt
       write(*,*) it,"'th time step, index ",jt(it)," = ",jt(it)*tkernel%dt," s"
       ! set values in tkernel%kernel and tkernel%jtcur manually
       tkernel%kernel = reshape(tkernel_kernel_all_t(it,:),(/tkernel%ntot_invgrid,ncomp,tkernel%nparam/))
       tkernel%jtcur = jt(it)
!
       ! write time kernel for this time step
       call writeTimeWaveformKernel(tkernel)
    end do ! it
!
    ! close tkernel file
    call finalWriteTimeWaveformKernel(tkernel,errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! clean up
1   if(allocated(tkernel_kernel_all_t)) deallocate(tkernel_kernel_all_t)
    if(allocated(skernel_vals_trans)) deallocate(skernel_vals_trans)
    if(allocated(vals_tmp)) deallocate(vals_tmp)
    call dealloc(DFT)
    call deallocateTimeWaveformKernel(tkernel)
  end subroutine transformSpectralToTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given parameter
!
  function getValuesByParamTimeWaveformKernel(this,param) result(p)
    type (time_waveform_kernel) :: this
    character(len=*) :: param
    real, dimension(:,:), pointer :: p
    nullify(p)
    if(.not.associated(this%kernel)) return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    p => this%kernel(:,:,indexOfParamModelParametrization(this%parametrization,param))
  end function getValuesByParamTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given component and parameter
!
  function getValuesByCompParamTimeWaveformKernel(this,comp,param) result(p)
    type (time_waveform_kernel) :: this
    character(len=*) :: param,comp
    real, dimension(:), pointer :: p
    integer :: jcomp
    nullify(p)
    if(.not.associated(this%kernel)) return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    do jcomp=1,this%ncomp
       if(this%comp(jcomp)==comp) then
          p => this%kernel(:,jcomp,indexOfParamModelParametrization(this%parametrization,param))
          exit
       end if
    end do ! jcomp
  end function getValuesByCompParamTimeWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given index of parameter
!
  function getValuesByIndexTimeWaveformKernel(this,i) result(p)
    type (time_waveform_kernel) :: this
    integer :: i
    real, dimension(:,:), pointer :: p
    nullify(p)
    if(.not.associated(this%kernel)) return
    if(i<1 .or. i>this%nparam) return
    p => this%kernel(:,:,i)
  end function getValuesByIndexTimeWaveformKernel
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
  function getNtotInvgridTimeWaveformKernel(this) result(ntot)
    type (time_waveform_kernel), intent(in) :: this
    integer :: ntot
    ntot = this%ntot_invgrid
  end function getNtotInvgridTimeWaveformKernel
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
    ios = createFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to write, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%filestat = 1
!
    if(present(nt)) then
       ! it is more efficient to create root like this, knowing the number of subgroups to come
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=nt,maxdset=3)
    else
       ! otherwise assume enough subgroups, which is still more efficient than creating root without subrgoup information
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=150,maxdset=3)
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
  subroutine writeTimeWaveformKernel(this)
    type (time_waveform_kernel) :: this
    ! local
    integer :: iparam,icomp,maxdataset
    type (group_stream_access) :: time
    type (data_stream_access) :: dset
!
    ! if file is not open to write, return
    if(this%filestat /= 1) return
!
    ! if there is an indicator that there are no sensible kernel values stored, just return doing nothing
    if(this%jtcur == -1) return
    if(this%ncomp == 0) return
!
    maxdataset = this%ncomp*this%nparam
    call createGroupStreamAccess(time,'Time',this%jtcur,maxsubgroup=0,maxdset=maxdataset)
    do iparam = 1,this%nparam
       do icomp = 1,this%ncomp
          call new(dset,1,(/ this%ntot_invgrid /),T_REAL)
          call writeDatasetVectorStreamAccess(dset,this%fsa,this%kernel(:,icomp,iparam))
          call addDatasetStreamAccess(time,dset)
          call dealloc(dset)
       end do ! icomp
    end do ! iparam
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
    integer :: icomp
!
    ! write general information to first dataset in group root
    allocate(ft(7))
    call new(dset,1,(/ 7 /),T_FLEXIBLE)
    ft(1) = this%parametrization
    ft(2) = this%nparam
    ft(3) = this%ntot_invgrid
    ft(4) = this%ncomp
    ft(5) = this%dt
    ft(6) = this%nt
    ft(7) = this%t0
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write components to second dataset in group root
    if(this%ncomp > 0) then
       allocate(ft(this%ncomp))
       call new(dset,1,(/ this%ncomp /),T_FLEXIBLE)
       do icomp = 1,this%ncomp
          ft(icomp) = this%comp(icomp)
       end do ! icomp
       call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
       call addDatasetStreamAccess(this%root,dset)
       deallocate(ft); call dealloc(dset)
    else
       call add(errmsg,2,"No station components defined for this time kernel, "//&
            "i.e. kernel was not initiated yet",myname)
       goto 1
    end if
!
    ! write time indices to third dataset in group root
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
1   call clearGroupTree(this%root)
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call dealloc(this%fsa)
    this%filestat = 0
  end subroutine finalWriteTimeWaveformKernel
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
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    type (flexible), dimension(:), pointer :: ft
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_component), dimension(:), allocatable :: comp
    integer :: nparam,ntot_invgrid,ncomp,icomp,nt
    real :: dt,t0
    integer, dimension(:), pointer :: jt
    character(len=400) :: errstr
    character(len=29) :: myname = 'initialReadTimeWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    ! open stream access file
    ios = openFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to read, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! read content information of stream access file
    call readGroupStreamAccess(this%root,this%fsa)
!
    ! read first info vector
    call traversePathStreamAccess(this%root,0,(/ 1 /),group,dset)
    call readDatasetVectorStreamAccess(dset,this%fsa,ft)
    if(size(ft)/=7) then
       write(errstr,*) "size of first information vector is ",size(ft),", must be 7. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    parametrization = ft(1)
    nparam = ft(2)
    ntot_invgrid = ft(3)
    ncomp = ft(4)
    dt = ft(5)
    nt = ft(6)
    t0 = ft(7)
    deallocate(ft)
    if(numberOfParamModelParametrization(parametrization) /= nparam) then
       write(errstr,*) "number of parameters ",nparam,", contained in this file is inconsistent with its parametrization '"&
            //trim(parametrization)//"' which has ",numberOfParamModelParametrization(parametrization),&
            " parameters, i.e. module modelParametrization has been modified since creation of the file"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    ! read second info vector containing the components (if it was written at all)
    if(ncomp > 0) then
       call traversePathStreamAccess(this%root,0,(/ 2 /),group,dset)
       call readDatasetVectorStreamAccess(dset,this%fsa,ft)
       if(size(ft)/=ncomp) then
          write(errstr,*) "info on components inconsistent: number of components ",size(ft),&
               " differs from expected number ",ncomp
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       allocate(comp(ncomp))
       do icomp=1,ncomp
          comp(icomp) = ft(icomp)
       end do ! icomp
       deallocate(ft)
    else
       write(errstr,*) "number of components of time kernel = ",ncomp,", must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1       
    end if

    ! read third info vector containing the time sample values (if it was written at all)
    ! if nt == 0, it means that there is no data contained in the file
    if(nt > 0) then
       call traversePathStreamAccess(this%root,0,(/ 3 /),group,dset)
       call readDatasetVectorStreamAccess(dset,this%fsa,jt)
       if(size(jt)/=nt) then
          write(errstr,*) "information on time samples inconsistent: number of time indices ",size(jt),&
               " differs from expected number ",nt,"; dt,t0 = ",dt,t0
          write(*,*) "ERROR initialReadTimeWaveformKernel: time indices = ",this%jt
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       this%t0 = t0
       this%dt = dt
       this%nt = nt
       this%jt => jt
    end if
!
    ! finally initiate timeWaveformKernel object
    call initiateTimeWaveformKernel(this,parametrization,ntot_invgrid,comp,errmsg)
    if(.level.errmsg ==2) goto 1
    deallocate(comp)    
!
    ! if everything went alright, return here, after setting flag this%filestat = 2 (file is open to read)
    this%filestat = 2
    return
!
    ! if there went anything wrong in intial read, close file here and revert all initiation that's been done so far
1   call deallocateTimeWaveformKernel(this)
    if(allocated(comp)) deallocate(comp)
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
    integer :: it,iparam,icomp,igroup,idset
    character(len=400) :: errstr
    character(len=22) :: myname = 'readTimeWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    if(this%filestat /= 2) then
       call add(errmsg,2,"there was no file opened to read yet, call initialReadTimeWaveformKernel first",myname)
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
            this%dt,this%t0!,"; time indices = ",this%jt
       write(*,*) "ERROR in readTimeWaveformKernel: time indices = ",this%jt
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    idset = 0
    do iparam = 1,this%nparam
       do icomp = 1,this%ncomp
          idset = idset + 1
          if(idset == 1) then
             ! remember time index group to make calls traversePathStreamAccess faster below
             call traversePathStreamAccess(this%root,1,(/ igroup,idset /),tgroup,dset)
          else
             call traversePathStreamAccess(tgroup,0,(/ idset /),group,dset)
          end if
          call readDatasetVectorStreamAccess(dset,this%fsa,d)
          if(size(d) /= this%ntot_invgrid) then
             write(errstr,*) "number of kernel values ",size(d)," for component ",icomp," = '"//trim(this%comp(icomp))//&
                  "' of ",iparam,"'th parameter '",trim(getParamFromIndexModelParametrization(this%parametrization,iparam)),&
                  "' differs from number of inversion grid cells ",this%ntot_invgrid,", i.e. kernel file may be corrupt"
             call add(errmsg,2,errstr,myname)
             return
          end if
          this%kernel(:,icomp,iparam) = d
          deallocate(d)
       end do ! icomp
    end do ! iparam
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
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call deallocateTimeWaveformKernel(this)
  end subroutine finalReadTimeWaveformKernel
!
end module timeWaveformKernel
