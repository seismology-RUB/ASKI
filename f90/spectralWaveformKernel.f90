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
!> \brief module which computes integrated spectral waveform sensitivity kernels
!!
!! \details Accessing the generic modules kernelDisplacement, kernelGreenTensor and
!!  kernelReferenceModel, this module computes waveform sensitivity kernels on the 
!!  wavefield points at the given frequencies for parametrizations consistent with 
!!  module modelParametrization and integrates them onto the inversion grid using 
!!  an instance of module integrationWeights.
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module spectralWaveformKernel
!
  use modelParametrization
  use kernelReferenceModel
  use kernelDisplacement
  use kernelGreenTensor
  use componentTransformation, only: character_length_component
  use integrationWeights
  use flexibleType
  use streamAccess
  use realloc
  use mathConstants
  use errorMessage
!
  implicit none
!
  private :: privateReadMetaInfoSpectralWaveformKernel,privateComputeSpectralWaveformKernel,&
       computeOnWpIsoLameSpectralWaveformKernel,computeOnWpIsoVelocitySpectralWaveformKernel,&
       integrateSpectralWaveformKernel
!
  interface dealloc; module procedure deallocateSpectralWaveformKernel; end interface
  interface onInvgrid; module procedure getOnInvgridSpectralWaveformKernel; end interface
  interface onWp; module procedure getNotOnInvgridSpectralWaveformKernel; end interface
  interface computeSpectralWaveformKernel
     module procedure computeOnInvgridSpectralWaveformKernel
     module procedure computeOnWpSpectralWaveformKernel
  end interface computeSpectralWaveformKernel
  interface operator (.pmtrz.); module procedure getParametrizationSpectralWaveformKernel; end interface
  interface operator (.param.); module procedure getParamSpectralWaveformKernel; end interface
  interface operator (.ntot.); module procedure getNtotKernelSpectralWaveformKernel; end interface
  interface operator (.df.); module procedure getDfSpectralWaveformKernel; end interface
  interface operator (.nfreq.); module procedure getNfreqSpectralWaveformKernel; end interface
  interface operator (.jf.); module procedure getJfSpectralWaveformKernel; end interface
  interface operator (.ncomp.); module procedure getNcompSpectralWaveformKernel; end interface
  interface operator (.comp.); module procedure getCompSpectralWaveformKernel; end interface
  interface operator (.kdID.); module procedure getKdIdSpectralWaveformKernel; end interface
  interface operator (.kgtID.); module procedure getKgtIdSpectralWaveformKernel; end interface
  interface operator (.oninvgrid.); module procedure getOnInvgridSpectralWaveformKernel; end interface
  interface operator (.onwp.); module procedure getNotOnInvgridSpectralWaveformKernel; end interface

!
  type spectral_waveform_kernel
     private
     ! model parametrization
     character(len=character_length_pmtrz) :: parametrization = '' !< model parametrization of kernel (also defines number of parameters nparam, below)
     integer :: nparam = 0 !< nparam = number of model parameters for which there are kernel values in this object (size of second dimension of kernel array)
     character(len=character_length_param), dimension(:), pointer :: param => null() !< vector of length nparam indicating the actual names of model parameters for which there are kernel values in this object
     integer :: nparam_in_file = 0 !< number of model parameters contained in kernel file, ONLY USED FOR READING KERNEL FILES!
     integer, dimension(:), pointer :: param_indx_in_file => null() !< mapping ONLY USED FOR READING KERNEL FILES! dim(nparam) . param_indx_in_file(iparam) = index of param used to generate data sets at the time writing the kernel file by stream access
     integer :: ntot_kernel = 0 !< total number of kernel values (either inversion grid cells onto which kernels are pre-integrated, or number of wavefield points, dependent on flag this%on_invgrid)
     ! spectral discretization
     real :: df = -1. !< df
     integer :: nfreq = 0 !< number of frequencies, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null() !< frequency indices in order of subgroups of kernel file; frequencies compute as f = jf*df [Hz]
     integer :: jfcur = -1 !< indicates current frequency f=jfcur*df for which kernel values are contained in array kernel
     ! source
     character(len=length_ID_kernel_displacement) :: kd_id = '' !< .id.kd for kernel_displacement object kd from which this kernel was created
     ! receiver
     character(len=length_ID_kernel_green_tensor) :: kgt_id = '' !< .id.kgt for kernel_green_tensor object kgt from which this kernel was created
     integer :: ncomp !< number of receiver components for which there are kernel values
     character(len=character_length_component), dimension(:), pointer :: comp => null() !< vector of length ncomp giving the actual names of the receiver compoments (this vector defines the order of values in array kernel)
     integer :: ncomp_in_file = 0 !< number of receiver components contained in kernel file, ONLY USED FOR READING KERNEL FILES!
     integer, dimension(:), pointer :: comp_indx_in_file => null() !< mapping ONLY USED FOR READING KERNEL FILES! dim(ncomp) . comp_indx_in_file(icomp) = index of comp used to generate data sets at the time writing the kernel file by stream access
     ! kernel values
     logical :: on_invgrid = .true. !< indicates whether this%kernel holds pre-integrated values on inversion grid (.true. , default) or on wavefield points (.false.)
     complex, dimension(:,:,:), pointer :: kernel => null() !< pre-integrated kernel values on invgrid (or plain kernel values on wp) for one frequency dim(ntot_kernel,nparam,ncomp)     ! file handling
     integer :: kernel_stat = -1 !< indicating status of kernel object: -1 = closed and not initiated, 0 = closed and initiated, 1 = open to write, 2 = open to read
     type (file_stream_access) :: fsa
     type (group_stream_access) :: root
  end type spectral_waveform_kernel
!
contains
!
!------------------------------------------------------------------------
!> \brief initiate spectral waveform kernel and allocate kernel array
!! \param this spectral waveform kernel
!! \param parametrization valid model parametrization, will be the parametrization of the kernel
!! \param ntot_invgrid total number of inversion grid cells (kernel values)
!! \param errmsg error message
!
  subroutine initiateSpectralWaveformKernel(this,parametrization,param,ntot_kernel,comp,errmsg,kernel_on_wp)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: parametrization
    integer :: ntot_kernel
    character(len=*), dimension(:) ::  comp,param
    type (error_message) :: errmsg
    logical, optional :: kernel_on_wp
    ! local
    character(len=400) :: errstr
    character(len=30) :: myname = 'initiateSpectralWaveformKernel'
    logical :: on_invgrid
    integer :: icomp,iparam
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
    if(size(comp) <= 0) then
       call add(errmsg,2,"there are no incoming components",myname)
       return
    end if
    if(.not.allValidComponents(comp,i_invalid=icomp)) then
       write(errstr,*) icomp,"'th incoming component '"//trim(comp(icomp))//"' not valid. Valid components are '"//&
            all_valid_components//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    do icomp = 2,size(comp)
       if(any(comp(1:icomp-1)==comp(icomp))) then
          write(errstr,*) icomp,"'th requested component '",trim(comp(icomp)),"' occurs more than once in the list of the ",&
               size(comp)," requested components. there must not be duplicate components in the request"
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
    this%ncomp = size(comp)
    allocate(this%comp(this%ncomp))
    this%comp = comp
!
    this%jfcur = -1
    allocate(this%kernel(this%ntot_kernel,this%nparam,this%ncomp))
!
    if(this%on_invgrid) then
       write(errstr,*) "initiated pre-integrated spectral waveform kernel on ",this%ntot_kernel," inversion grid cells, ",&
            this%nparam," parameters of parametrization '",trim(this%parametrization),"', and ",this%ncomp,&
            " station components"
    else
       write(errstr,*) "initiated plain spectral waveform kernel on ",this%ntot_kernel," wavefield points, ",&
            this%nparam," parameters of parametrization '",trim(this%parametrization),"', and ",this%ncomp,&
            " station components"
    end if
    call add(errmsg,0,errstr,myname)
!
    this%kernel_stat = 0
  end subroutine initiateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this spectral waveform kernel
!
  subroutine deallocateSpectralWaveformKernel(this)
    type (spectral_waveform_kernel) :: this
    this%parametrization = ''
    this%nparam = 0
    if(associated(this%param)) deallocate(this%param)
    this%nparam_in_file = 0
    if(associated(this%param_indx_in_file)) deallocate(this%param_indx_in_file)
    this%ntot_kernel = 0
    this%df = -1.
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
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
  end subroutine deallocateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief invgrid fork of computeSpectralWaveformKernel: private forwarding to joint routine privateComputeSpectralWaveformKernel
!! \details this mechanism is implemented, in order to have an explicit interface 
!!  for calls to computeSpectralWaveformKernel from outside (without any optional dummy variables).
!!  Still, all forks of computeSpectralWaveformKernel need mainly to do the same thing, thus it is feasible to have 
!!  ONE JOINT internal routine privateComputeSpectralWaveformKernel (which takes optional arguments, etc.)
!!  to process all forks of computeSpectralWaveformKernel
  subroutine computeOnInvgridSpectralWaveformKernel(this,kd,kgt,krm,uf_mdata,intw,errmsg)
    type (spectral_waveform_kernel) :: this
    type (kernel_displacement) :: kd
    type (kernel_green_tensor) :: kgt
    type (kernel_reference_model) :: krm
    real :: uf_mdata
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    character (len=38) :: myname = 'computeOnInvgridSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
    if(this%kernel_stat == -1) then
       call add(errmsg,2,"spectral_waveform_kernel object not initiated yet",myname)
       return
    end if
    if(this%kernel_stat == 2) then
       call add(errmsg,2,&
            "spectral_waveform_kernel object is opened to read. computing kernel values is not permitted in this case",&
            myname)
       return
    end if
!
    if(.not.this%on_invgrid) then
       call add(errmsg,2,"spectral_waveform_kernel object was initiated for plain kernel values on wavefield points; "//&
            "cannot compute kernel on inversion grid in this case",myname)
       return
    end if
!
    call privateComputeSpectralWaveformKernel(this,kd,kgt,krm,uf_mdata,errmsg,myname,intw=intw)
  end subroutine computeOnInvgridSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief wp fork of computeSpectralWaveformKernel: private forwarding to joint routine privateComputeSpectralWaveformKernel
!! \details this mechanism is implemented, in order to have an explicit interface 
!!  for calls to computeSpectralWaveformKernel from outside (without any optional dummy variables).
!!  Still, all forks of computeSpectralWaveformKernel need mainly to do the same thing, thus it is feasible to have 
!!  ONE JOINT internal routine privateComputeSpectralWaveformKernel (which takes optional arguments, etc.)
!!  to process all forks of computeSpectralWaveformKernel
  subroutine computeOnWpSpectralWaveformKernel(this,kd,kgt,krm,uf_mdata,errmsg)
    type (spectral_waveform_kernel) :: this
    type (kernel_displacement) :: kd
    type (kernel_green_tensor) :: kgt
    type (kernel_reference_model) :: krm
    real :: uf_mdata
    type (error_message) :: errmsg
    ! local
    character (len=33) :: myname = 'computeOnWpSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
    if(this%kernel_stat == -1) then
       call add(errmsg,2,"spectral_waveform_kernel object not initiated yet",myname)
       return
    end if
    if(this%kernel_stat == 2) then
       call add(errmsg,2,&
            "spectral_waveform_kernel object is opened to read. computing kernel values is not permitted in this case",&
            myname)
       return
    end if
!
    if(this%on_invgrid) then
       call add(errmsg,2,"spectral_waveform_kernel object was initiated for pre-integrated kernel values on inversion "//&
            "grid; cannot compute kernel on wavefield points in this case",myname)
       return
    end if
!
    call privateComputeSpectralWaveformKernel(this,kd,kgt,krm,uf_mdata,errmsg,myname)
  end subroutine computeOnWpSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute spectral waveform kernel at current frequency and pre-integrate on inversion grid
!! \details This is a generic routine which, dependent on the parametrization of this kernel
!!  calls specific subroutines to compute the kernel values
!! \param this spectral waveform kernel
!! \param kd kernel displacement (must contain displacement values at same frequency as kgt)
!! \param kgt kernel green tensor (must contain green tensor values at same frequency as kd)
!! \param krm kernel reference model object
!! \param uf_mdata unit factor of the measured data
!! \param intw optional integration weights which connect the wavefield points on which kd,kgt are given with inversion grid on which kernel is computed (IN CASE this%on_invgrid)
!! \param errmsg error message
!
  subroutine privateComputeSpectralWaveformKernel(this,kd,kgt,krm,uf_mdata,errmsg,myname,intw)
    type (spectral_waveform_kernel) :: this
    type (kernel_displacement) :: kd
    type (kernel_green_tensor) :: kgt
    type (kernel_reference_model) :: krm
    real :: uf_mdata
    type (integration_weights), optional :: intw
    type (error_message) :: errmsg
    character(len=*) :: myname
    ! local
    real :: df_kd,uf_u,uf_ustr,df_kgt,uf_g,uf_gstr,uf_param,uf_intw
    double precision, dimension(:), allocatable :: uf_equation
    integer :: jfcur_kd,jfcur_kgt,nwp
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    double complex, dimension(:,:,:), allocatable :: kernel_on_wp
    integer :: iparam
    character(len=400) :: errstr
!
    ! check if kd and kgt objects have same (positive) df and same jfcur (i.e. contain values for the very same frequency)
    df_kd = .df.kd; df_kgt = .df.kgt
    if(df_kd .le. 0. .or. df_kgt .le. 0.) then
       write(errstr,*) "frequency step df of kernel displacement ( = ",df_kd,&
            ") or frequency step of kernel green tensor ( = ",df_kgt,") invalid: must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if( abs(df_kd-df_kgt)/df_kd > 1.e-4 ) then
       write(errstr,*) "frequency step df of kernel displacement ( = ",df_kd,&
            ") differs from frequency step df of kernel green tensor ( = ",df_kgt,") by more "//&
            "than 0.01 percent"
       call add(errmsg,2,errstr,myname)
       return
    end if
    jfcur_kd = .jfcur.kd; jfcur_kgt = .jfcur.kgt
    if(jfcur_kd < 0 .or. jfcur_kgt < 0) then
       write(errstr,*) "current frequency index of kernel displacement ( = ",jfcur_kd,&
            ") or current frequency index of kernel green tensor ( = ",jfcur_kgt,") invalid: must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(jfcur_kd /= jfcur_kgt) then
       write(errstr,*) "current frequency index of kernel displacement ( = ",jfcur_kd,&
            ") differs from current frequency index of kernel green tensor ( = ",jfcur_kgt,")"
       call add(errmsg,2,errstr,myname)
       return
    end if
    ! check if df of kd and kgt is consistent with my df
    if(this%df<0.) then
       ! if this is the first computed frequency, just use the df of kd and kgt
       this%df = df_kd
    else
       ! if this is not the first computed frequency, check if df is consistent
       if( (this%df-df_kd)/this%df > 1.e-4 ) then
          write(errstr,*) "frequency step df of kernel displacement and kernel green tensor ( = ",df_kd,&
               ") differs from frequency step df of this waveform kernel ( = ",this%df,")"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
    ! check if the id of kd is still the same that is memorized by this object
    if(this%kd_id == '') then
       ! if this is the first computed frequency, just use the id of kd
       this%kd_id = .id.kd
    else
       if(this%kd_id /= .id.kd) then
          write(errstr,*) "ID of incoming kernel displacement object = '",trim(.id.kd),&
               "' differs from the ID memorized by this waveform kernel object = '",trim(this%kd_id),"'"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
    ! check if the id of kgt is still the same that is memorized by this object
    if(this%kgt_id == '') then
       ! if this is the first computed frequency, just use the id of kgt
       this%kgt_id = .id.kgt
    else
       if(this%kgt_id /= .id.kgt) then
          write(errstr,*) "ID of incoming kernel green tensor object = '",trim(.id.kgt),&
               "' differs from the ID memorized by this waveform kernel object = '",trim(this%kgt_id),"'"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ! check if kgt is initiated for the same components (and order) as this object was initiated
    if(.not.hasExactComponentsKernelGreenTensor(kgt,this%comp)) then
       write(errstr,*) "kernel green tensor object is not initiated for components '",this%comp//",",&
            "' (in this order). kgt object and this kernel object must be initiated for the very same components"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! test incoming unit factor of measured data (must be strictly positive)
    if(uf_mdata <= 0) then
       write(errstr,*) "incoming unit factor of measured data = ",uf_mdata,"; must be strictly positive!"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! get wavefields and strains
    call getKernelDisplacement(kd,u,errmsg)
    if(.level.errmsg == 2) goto 1
    call getStrainsKernelDisplacement(kd,ustr,errmsg)
    if(.level.errmsg == 2) goto 1
    call getKernelGreenTensor(kgt,g,errmsg)
    if(.level.errmsg == 2) goto 1
    if(size(g,3)/=this%ncomp) then
       write(errstr,*) "number of force components returned by routine getKernelGreenTensor = ",&
            size(g,3),"; expected value = ",this%ncomp
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    call getStrainsKernelGreenTensor(kgt,gstr,errmsg)
    if(.level.errmsg == 2) goto 1
    if(size(gstr,3)/=this%ncomp) then
       write(errstr,*) "number of force components returned by routine getKernelGreenTensor = ",&
            size(gstr,3),"; expected value = ",this%ncomp
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    nwp = size(u,1)
    if(nwp /= size(g,1)) then
       write(errstr,*) "number of kernel displacement values ",nwp," differs from number of green tensor values ",size(g,1)
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    ! check if number of wavefield points is consistent with this kernel object in case of kernel on wavefield points
    if(.not.present(intw)) then
       if(nwp/=this%ntot_kernel) then
          write(errstr,*) "number of kernel displacement values ",nwp,&
               " differs from number of wavefield points for which this kernel was initiated: ",this%ntot_kernel
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
    end if
!
    ! get the unit factors for u, ustr, g, gstr from the kernel(Displacment/GreenTensor) modules
    call getUnitFactorKernelDisplacement(kd,uf_u,errmsg)
    if(.level.errmsg == 2) return
    call getUnitFactorStrainsKernelDisplacement(kd,uf_ustr,errmsg)
    if(.level.errmsg == 2) return
    call getUnitFactorKernelGreenTensor(kgt,uf_g,errmsg)
    if(.level.errmsg == 2) return
    call getUnitFactorStrainsKernelGreenTensor(kgt,uf_gstr,errmsg)
    if(.level.errmsg == 2) return
!
    ! communicate this frequency index to computing routines below via this object 
    this%jfcur = jfcur_kd
!
    ! Compute kernel factors which accounts for the rest of the kernel equation, outside the kernel expressions, 
    ! i.e. pre-integration, values of inverted model (i.e. need a vector of factors) and measured data
    ! These factors account for (1) unit factor of measured data (dividing by this factor), 
    ! (2) integration weights factor (if not OnWp), (3) model parameters factor
    ! Transfer these factors to the routines which compute the kernels. It should be accounted for inside the routines:
    !   first multiplying ALL factors (which might neturalize or at least compensate each other) 
    !   before multiplying to kernel values. This tries to avoid to lose precision or resulting 
    !   in zero or infty values which might happen when you multiply the kernel values by one factor after another
    allocate(uf_equation(this%nparam))
    if(present(intw)) then
       call getUnitFactorIntegrationWeights(intw,uf_intw,errmsg)
       if(.level.errmsg == 2) return
       uf_equation = dble(uf_intw) / dble(uf_mdata)
    else
       uf_equation = 1.d0 / dble(uf_mdata)
    end if
    do iparam = 1,this%nparam
       uf_param = getUnitFactorOfParamModelParametrization(this%parametrization,this%param(iparam)) ! do not check output by flag ios , since at this point this%param sould contain valid parameter names
       uf_equation(iparam) = uf_equation(iparam)*uf_param
    end do ! iparam
!
    allocate(kernel_on_wp(size(u,1),this%nparam,this%ncomp))
    select case(this%parametrization)
    case('isoLameSI')
       call computeOnWpIsoLameSpectralWaveformKernel(this,u,uf_u,ustr,uf_ustr,g,uf_g,gstr,uf_gstr,uf_equation,kernel_on_wp,errmsg)
    case('isoVelocitySI','isoVelocity1000')
       call computeOnWpIsoVelocitySpectralWaveformKernel(this,krm,u,uf_u,ustr,uf_ustr,g,uf_g,gstr,uf_gstr,uf_equation,&
            kernel_on_wp,errmsg)
       ! ADD YOUR NEW PARAMETRIZATION HERE, AND ADD A RESPECTIVE SUBROUTINE BELOW
       ! case('your_new_parametrization'); call computeYourNewParametrizationSpectralWaveformKernel(this,u,uf_u,ustr,uf_ustr,g,uf_g,gstr,uf_gstr,uf_equation,kernel_on_wp,errmsg)
    case default
       write(errstr,*) "there are no routines implemented to compute spectral waveform kernels of parametrization '",&
            trim(this%parametrization),"'"
       call add(errmsg,2,errstr,myname)
    end select
    if(.level.errmsg == 2) then
       deallocate(kernel_on_wp)
       goto 1
    end if
!
    if(present(intw)) then
       ! in this case it is assumed that this%on_invgrid is .true.
       ! integrate kernel values onto inversion grid
       call integrateSpectralWaveformKernel(this,kernel_on_wp,intw,errmsg,myname)
       if(.level.errmsg == 2) then
          deallocate(kernel_on_wp)
          goto 1
       end if
    else ! present(intw)
       ! in this case it is assumed that this%on_invgrid is .false.
       ! assign plain kernel values on wavefield points to this%kernel (finally converting to single precision here)
       this%kernel = kernel_on_wp
    end if ! present(intw)
    deallocate(kernel_on_wp)
!
1   if(associated(u)) deallocate(u)
    if(associated(ustr)) deallocate(ustr)
    if(associated(g)) deallocate(g)
    if(associated(gstr)) deallocate(gstr)
    if(allocated(uf_equation)) deallocate(uf_equation)
  end subroutine privateComputeSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute kernel values for parametrizations of type 'isoLame'
!! \param this spectral waveform kernel
!! \param u displacement values as returned by kernelDisplacement object
!! \param uf_u unit factor of kernel displacement
!! \param ustr displacement strains as returned by kernelDisplacement object
!! \param uf_ustr unit factor of kernel displacement strains
!! \param g green tensor values as returned by kernelGreenTensor object
!! \param uf_g unit factor of kernel green tensor
!! \param gstr green tensor strains as returned by kernelGreenTensor object
!! \param uf_gstr unit factor of kernel green tensor strains
!! \param uf_equation array of unit factors, one factor for each prameter in this%param, containing all remaining unit factors of the kernel equation (integration weights, model perturbation, data residual)
!! \param intw integration weights as given to computeSpectralWaveformKernel
!! \param errmsg error message
!
  subroutine computeOnWpIsoLameSpectralWaveformKernel(this,u,uf_u,ustr,uf_ustr,g,uf_g,gstr,uf_gstr,uf_equation,kernel_on_wp,errmsg)
    type (spectral_waveform_kernel) :: this
    real :: uf_u,uf_ustr,uf_g,uf_gstr
    double precision, dimension(:) :: uf_equation
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    double complex, dimension(:,:,:), intent(inout) :: kernel_on_wp
    type (error_message) :: errmsg
    ! local
    integer :: icomp,iparam
    real :: omega
    double precision :: uf_multiplied
    character(len=40) :: myname = 'computeOnWpIsoLameSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    omega = 2.*mc_pi*(this%jfcur*this%df)
!
    do iparam = 1,this%nparam

       select case(this%param(iparam))
       case ('rho')
          ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
          uf_multiplied = uf_u*uf_g*uf_equation(iparam)
          ! compute kernel on wavefield points for parameter rho
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = uf_multiplied*(omega*omega)*(dcmplx(u(:,1))*dcmplx(g(:,1,icomp)) + &
                  dcmplx(u(:,2))*dcmplx(g(:,2,icomp)) + dcmplx(u(:,3))*dcmplx(g(:,3,icomp)))
          end do ! icomp

       case ('lambda')
          ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
          uf_multiplied = uf_ustr*uf_gstr*uf_equation(iparam)
          ! compute kernel on wavefield points for parameter lambda
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = -uf_multiplied * (dcmplx(ustr(:,1))+dcmplx(ustr(:,2))+dcmplx(ustr(:,3))) * &
                  (dcmplx(gstr(:,1,icomp))+dcmplx(gstr(:,2,icomp))+dcmplx(gstr(:,3,icomp)))
          end do ! iparam

       case ('mu')
          ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
          uf_multiplied = uf_ustr*uf_gstr*uf_equation(iparam)
          ! compute kernel on wavefield points for third parameter mu
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = uf_multiplied * ( &
                  -2.d0*(dcmplx(ustr(:,1))*dcmplx(gstr(:,1,icomp)) + dcmplx(ustr(:,2))*dcmplx(gstr(:,2,icomp)) + &
                  dcmplx(ustr(:,3))*dcmplx(gstr(:,3,icomp))) &
                  -4.d0*(dcmplx(ustr(:,4))*dcmplx(gstr(:,4,icomp)) + dcmplx(ustr(:,5))*dcmplx(gstr(:,5,icomp)) + &
                  dcmplx(ustr(:,6))*dcmplx(gstr(:,6,icomp)))  )
          end do ! icomp
       end select ! this%param(iparam)

    end do ! iparam
!
  end subroutine computeOnWpIsoLameSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute kernel values for parametrizations of type 'isoVelocity'
!! \param this spectral waveform kernel
!! \param krm kernel reference model
!! \param u displacement values as returned by kernelDisplacement object
!! \param uf_u unit factor of kernel displacement
!! \param ustr displacement strains as returned by kernelDisplacement object
!! \param uf_ustr unit factor of kernel displacement strains
!! \param g green tensor values as returned by kernelGreenTensor object
!! \param uf_g unit factor of kernel green tensor
!! \param gstr green tensor strains as returned by kernelGreenTensor object
!! \param uf_gstr unit factor of kernel green tensor strains
!! \param uf_equation array of unit factors, one factor for each prameter in this%param, containing all remaining unit factors of the kernel equation (integration weights, model perturbation, data residual)
!! \param intw integration weights as given to computeSpectralWaveformKernel
!! \param errmsg error message
!
  subroutine computeOnWpIsoVelocitySpectralWaveformKernel(this,krm,u,uf_u,ustr,uf_ustr,g,uf_g,gstr,uf_gstr,uf_equation,&
       kernel_on_wp,errmsg)
    type (spectral_waveform_kernel) :: this
    type (kernel_reference_model) :: krm
    real :: uf_u,uf_ustr,uf_g,uf_gstr
    double precision, dimension(:) :: uf_equation
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    double complex, dimension(:,:,:), intent(inout) :: kernel_on_wp
    type (error_message) :: errmsg
    ! local
    integer :: icomp,iparam,nwp
    double precision :: omega
    double precision :: uf_multiplied,uf_rho,uf_vp,uf_vs
    real, dimension(:), pointer :: rho,vp,vs
    double complex, dimension(:,:), allocatable :: krho,klambda,kmu
    character(len=400) :: errstr
    character(len=44) :: myname = 'computeOnWpIsoVelocitySpectralWaveformKernel'
    logical :: any_rho,any_vp,any_vs
!
    call addTrace(errmsg,myname)
    nwp = size(u,1)
!
    any_rho = any(this%param == 'rho')
    any_vp = any(this%param == 'vp')
    any_vs = any(this%param == 'vs')
!
    omega = 2.d0*mc_pid*(this%jfcur*dble(this%df))
!
    ! FIRST GET ALL REQUIRED KERNEL REFERENCE MODEL VALUES
!
    ! rho values of kernel reference model are only required for vp and vs kernels
    if(any_vp .or. any_vs) then
       rho => getModelValuesKernelReferenceModel(krm,this%parametrization,'rho')
       if(.not.associated(rho)) then
          call add(errmsg,2,"no model values for 'isoVelocity'-'rho' in kernel reference model",myname)
          this%jfcur = -1
          goto 1
       end if
       if(size(rho) /= nwp) then
          write(errstr,*) "number of kernel reference model 'rho' values = ",size(rho),&
               " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          goto 1
       end if
       uf_rho = getUnitFactorOfParamModelParametrization(this%parametrization,'rho')
    end if ! any_vp .or. any_vs
!
    ! vp values of kernel reference model are only required for rho and vp kernels
    if(any_rho .or. any_vp) then
       vp => getModelValuesKernelReferenceModel(krm,this%parametrization,'vp')
       if(.not.associated(vp)) then
          call add(errmsg,2,"no model values for 'isoVelocity'-'vp' in kernel reference model",myname)
          this%jfcur = -1
          goto 1
       end if
       if(size(vp) /= nwp) then
          write(errstr,*) "number of kernel reference model 'vp' values = ",size(vp),&
               " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          goto 1
       end if
       uf_vp = getUnitFactorOfParamModelParametrization(this%parametrization,'vp')
    end if ! any_rho .or. any_vp
!
    ! vs values of kernel reference model are only required for rho and vs kernels
    if(any_rho .or. any_vs) then
       vs => getModelValuesKernelReferenceModel(krm,this%parametrization,'vs')
       if(.not.associated(vs)) then
          call add(errmsg,2,"no model values for 'isoVelocity'-'vs' in kernel reference model",myname)
          this%jfcur = -1
          goto 1
       end if
       if(size(vs) /= nwp) then
          write(errstr,*) "number of kernel reference model 'vs' values = ",size(vs),&
               " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          goto 1
       end if
       uf_vs = getUnitFactorOfParamModelParametrization(this%parametrization,'vs')
    end if ! any_rho .or. any_vs
!
    ! SECONDLY, COMPUTE ALL REQUIRED KERNELS FOR PARAMETRIZATION ISOLAME
!
    ! rho-isoLame kernel is only required for rho-isoVelocity kernel
    if(any_rho) then
       allocate(krho(nwp,this%ncomp))
       ! rho (of parametrization isoLame), account here for unit factors of krho kernel (uf_u*uf_g)
       ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
       uf_multiplied = uf_u*uf_g
       do icomp = 1,this%ncomp
          krho(:,icomp) = (uf_multiplied * omega*omega) * &
               (dcmplx(u(:,1))*dcmplx(g(:,1,icomp)) + dcmplx(u(:,2))*dcmplx(g(:,2,icomp)) + dcmplx(u(:,3))*dcmplx(g(:,3,icomp)))
       end do
    end if ! any_rho
!
    ! lambda-isoLame kernel is required for ALL rho-, vp- and vs-isoVelocity kernels 
    ! (so do not put an if clause, must be done in any case)
    allocate(klambda(nwp,this%ncomp))
    ! account here for unit factor of klambda kernel (uf_ustr*uf_gstr)
    ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
    uf_multiplied = uf_ustr*uf_gstr
    do icomp = 1,this%ncomp
       klambda(:,icomp) = -uf_multiplied * (dcmplx(ustr(:,1))+dcmplx(ustr(:,2))+dcmplx(ustr(:,3))) * &
            (dcmplx(gstr(:,1,icomp))+dcmplx(gstr(:,2,icomp))+dcmplx(gstr(:,3,icomp)))
    end do

    ! mu-isoLame kernel is required only for rho- and vs-isoVelocity kernels
    if(any_rho .or. any_vs) then
       allocate(kmu(nwp,this%ncomp))
       ! account here for unit factor of kmu kernel (uf_ustr*uf_gstr)
       ! multiply all involved unit factors in order to allow for them to compensate each other (may result in numerically more stable multiplication in kernel formulas below)
       uf_multiplied = uf_ustr*uf_gstr
       do icomp = 1,this%ncomp
          kmu(:,icomp) = uf_multiplied * ( &
               -2.d0*(dcmplx(ustr(:,1))*dcmplx(gstr(:,1,icomp)) + dcmplx(ustr(:,2))*dcmplx(gstr(:,2,icomp)) + &
               dcmplx(ustr(:,3))*dcmplx(gstr(:,3,icomp))) &
               -4.d0*(dcmplx(ustr(:,4))*dcmplx(gstr(:,4,icomp)) + dcmplx(ustr(:,5))*dcmplx(gstr(:,5,icomp)) + &
               dcmplx(ustr(:,6))*dcmplx(gstr(:,6,icomp)))  )
       end do
    end if ! any_rho .or. any_vs
!
    ! FINALLY, COMPUTE isoVelocity KERNELS VIA LINEARIZED RELATIONS BETWEEN PARAMETRIZATION RHO,LAMBDA,MU AND RHO,VP,VS
    ! (derived e.g. by total derivative)
!
    do iparam = 1,this%nparam
       select case(this%param(iparam))

       case ('rho')
          ! compute kernel on wavefield points for isoVelocity parameter rho
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = uf_equation(iparam) * ( krho(:,icomp) + &
                  ((uf_vp*dble(vp(:)))**2-2.d0*(uf_vs*dble(vs(:)))**2) * klambda(:,icomp) + &
                  ((uf_vs*dble(vs(:)))**2) * kmu(:,icomp) )
          end do 

       case ('vp')
          ! compute kernel on wavefield points for isoVelocity parameter vp
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = uf_equation(iparam) * &
                  ( 2.d0*uf_rho*dble(rho(:))*uf_vp*dble(vp(:)) * klambda(:,icomp) )
          end do

       case ('vs')
          ! compute kernel on wavefield points for isoVelocity parameter vs
          do icomp = 1,this%ncomp
             kernel_on_wp(:,iparam,icomp) = uf_equation(iparam) * &
                  (2.d0*uf_rho*dble(rho(:))*uf_vs*dble(vs(:))*(kmu(:,icomp)-2.d0*klambda(:,icomp)))
          end do

       end select
    end do ! iparam

1   if(allocated(krho)) deallocate(krho)
    if(allocated(klambda)) deallocate(klambda)
    if(allocated(kmu)) deallocate(kmu)
    if(associated(rho)) deallocate(rho)
    if(associated(vp)) deallocate(vp)
    if(associated(vs)) deallocate(vs)
!
  end subroutine computeOnWpIsoVelocitySpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief given kernel values on wavefield points, do integration onto inversion grid
!! \param this spectral waveform kernel
!! \param kernel_on_wp kernel values computed on wavefield points, as returned by private routines of this module
!! \param intw integration weights to integrate values onto inversion grid cells
!! \param errmsg error message
!! \param myname name of routine which called this routine integrateSpectralWaveformKernel (error will look like raised there)
!
  subroutine integrateSpectralWaveformKernel(this,kernel_on_wp,intw,errmsg,myname)
    type (spectral_waveform_kernel) :: this
    double complex, dimension(:,:,:) :: kernel_on_wp
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    character(len=*) :: myname
    ! local
    integer :: nwp_kernel,icell,iparam,icomp
    integer, dimension(:), pointer :: wp_idx
    real, dimension(:), pointer :: weight
    character(len=400) :: errstr
!
    ! IN THIS PRIVATE ROUTINE IT IS ASSUMED, THAT this%on_invgrid == .true. , i.e. this%ntot_kernel means 
    ! number of inversion grid cells (which must be consistent with object intw)
!
    nwp_kernel = size(kernel_on_wp,1)
!
    do icell = 1,this%ntot_kernel
       wp_idx => intw.wpidx.icell
       weight => intw.weight.icell
       if(.not.associated(wp_idx)) then
          write(errstr,*) "there are no integration weights for inversion grid cell ",icell,&
               ", inversion grid cell index seems to be out of range"
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          return
       end if
       if(any(wp_idx>nwp_kernel)) then
          write(errstr,*) "some wavefield point indices of integration weights for inversion grid cell ",icell,&
               ", exceed the maximum number of wavefield points for which this kernel was initiated. ",&
               "Integration weights object seems to be corrupt"
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          return
       end if
       do icomp = 1,this%ncomp
          do iparam = 1,this%nparam
             this%kernel(icell,iparam,icomp) = sum(dble(weight)*kernel_on_wp(wp_idx,iparam,icomp))
          end do ! iparam
       end do ! icomp
    end do ! icell
  end subroutine integrateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given component
!
  function getValuesByCompSpectralWaveformKernel(this,comp) result(p)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: comp
    complex, dimension(:,:), pointer :: p
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
  end function getValuesByCompSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given parameter
!
  function getValuesByParamSpectralWaveformKernel(this,param) result(p)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: param
    complex, dimension(:,:), pointer :: p
    integer :: iparam
    nullify(p)
    if(this%kernel_stat == -1) return
    if(.not.associated(this%kernel)) return !! should be redundant with if(this%kernel_stat == -1) return
    do iparam = 1,this%nparam
       if(this%param(iparam) == param) then
          p => this%kernel(:,iparam,:)
          return
       end if
    end do ! iparam
  end function getValuesByParamSpectralWaveformKernel
! !------------------------------------------------------------------------
! !> \brief return kernel values for given index of parameter
! !
!   function getValuesByIndexSpectralWaveformKernel(this,i) result(p) !! ROUTINE IS DEPRECIATED (would need to refer by index to the current sorting and content of vector this%param, so it's safer to use above routine getValuesByParamSpectralWaveformKernel)
!     type (spectral_waveform_kernel) :: this
!     integer :: i
!     complex, dimension(:,:), pointer :: p
!     nullify(p)
!     if(.not.associated(this%kernel)) return
!     if(i<1 .or. i>this%nparam) return
!     p => this%kernel(:,i,:)
!   end function getValuesByIndexSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to all kernel values
!
  function getAllValuesSpectralWaveformKernel(this) result(p)
    type (spectral_waveform_kernel) :: this
    complex, dimension(:,:,:), pointer :: p
    p => this%kernel
  end function getAllValuesSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return parametrization of kernel
!
  function getParametrizationSpectralWaveformKernel(this) result(pmtrz)
    type (spectral_waveform_kernel), intent(in) :: this
    character(len=character_length_pmtrz) :: pmtrz
    pmtrz = this%parametrization
  end function getParametrizationSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to array of parameter names
!
  function getParamSpectralWaveformKernel(this) result(param)
    type (spectral_waveform_kernel), intent(in) :: this
    character(len=character_length_param), dimension(:), pointer :: param
    param => this%param
  end function getParamSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of inversion grid cells of kernel
!
  function getNtotKernelSpectralWaveformKernel(this) result(ntot)
    type (spectral_waveform_kernel), intent(in) :: this
    integer :: ntot
    ntot = this%ntot_kernel
  end function getNtotKernelSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return df of kernel
!
  function getDfSpectralWaveformKernel(this) result(df)
    type (spectral_waveform_kernel), intent(in) :: this
    real :: df
    df = this%df
  end function getDfSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of frequencies of kernel
!
  function getNfreqSpectralWaveformKernel(this) result(nfreq)
    type (spectral_waveform_kernel), intent(in) :: this
    integer :: nfreq
    nfreq = this%nfreq
  end function getNfreqSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to array of frequency indices of kernel
!
  function getJfSpectralWaveformKernel(this) result(jf)
    type (spectral_waveform_kernel), intent(in) :: this
    integer, dimension(:), pointer :: jf
    jf => this%jf
  end function getJfSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return number of receiver components of kernel
!
  function getNcompSpectralWaveformKernel(this) result(ncomp)
    type (spectral_waveform_kernel), intent(in) :: this
    integer :: ncomp
    ncomp = this%ncomp
  end function getNcompSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return pointer to array of receiver components of kernel
!
  function getCompSpectralWaveformKernel(this) result(comp)
    type (spectral_waveform_kernel), intent(in) :: this
    character(len=character_length_component), dimension(:), pointer :: comp
    comp => this%comp
  end function getCompSpectralWaveformKernel
!-------------------------------------------------------------------------
!> \brief Get ID of kernel displacement object, which is associated with this kernel
!
  function getKdIdSpectralWaveformKernel(this) result(id)
    type (spectral_waveform_kernel), intent(in) :: this
    character(len=length_ID_kernel_displacement) :: id
    id = this%kd_id
  end function getKdIdSpectralWaveformKernel
!-------------------------------------------------------------------------
!> \brief Get ID of kernel green tensor object, which is associated with this kernel
!
  function getKgtIdSpectralWaveformKernel(this) result(id)
    type (spectral_waveform_kernel), intent(in) :: this
    character(len=length_ID_kernel_green_tensor) :: id
    id = this%kgt_id
  end function getKgtIdSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return .not.this%on_invgrid (indicating that the values of this are on wavefield points)
!
  function getNotOnInvgridSpectralWaveformKernel(this) result(l)
    type (spectral_waveform_kernel), intent(in) :: this
    logical :: l
    if(this%kernel_stat == -1) then
       ! cannot really indicate that object is not yet initiated, but return false in that case
       l = .false.
    else
       l = .not.(this%on_invgrid)
    end if
  end function getNotOnInvgridSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return this%on_invgrid (indicating that the values of this are on inversion grid, pre-integrated)
!
  function getOnInvgridSpectralWaveformKernel(this) result(l)
    type (spectral_waveform_kernel), intent(in) :: this
    logical :: l
    if(this%kernel_stat == -1) then
       ! cannot really indicate that object is not yet initiated, but return false in that case
       l = .false.
    else
       l = this%on_invgrid
    end if
  end function getOnInvgridSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief open stream access file to write kernel values, one frequency at a time
!! \param 
!
  subroutine initialWriteSpectralWaveformKernel(this,filename,lu,errmsg,nfreq)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    integer, optional :: nfreq
    ! local
    integer :: ios
    character(len=400) :: errstr
    character(len=34) :: myname = 'initialWriteSpectralWaveformKernel'
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
    ios = createFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to write, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%kernel_stat = 1
!
    if(present(nfreq)) then
       ! it is more efficient to create root like this, knowing the number of subgroups to come
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=nfreq,maxdset=4)
    else
       ! otherwise assume enough subgroups, which is still more efficient than creating root without subrgoup information
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=150,maxdset=4)
    end if
  end subroutine initialWriteSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief write current values of waveform kernel to file
!! \details It is important that the frequencies are written one after another, 
!!  i.e. first frequency index 1, then 2, then ... etc., because the frequencies
!!  are just appended to the stream access file and when reading the file again, 
!!  it is assumed that the frequencies were written to file in index order.
!! \param 
!
  subroutine writeSpectralWaveformKernel(this,errmsg)
    type (spectral_waveform_kernel) :: this
    type (error_message) :: errmsg
    ! local
    integer :: iparam,icomp,maxdataset
    type (group_stream_access) :: frequency
    type (data_stream_access) :: dset
    character(len=27) :: myname = 'writeSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    ! if file is not open to write, return
    if(this%kernel_stat /= 1) then
       call add(errmsg,2,"file not yet opened to write, call initialWriteSpectralWaveformKernel first",myname)
       return
    end if
!
    ! if there is an indicator that there are no sensible kernel values stored, just return doing nothing
    if(this%jfcur == -1) then
       call add(errmsg,2,"no sensible kernel values were computed, cannot write any values",myname)
       return
    end if
!
    maxdataset = this%ncomp*this%nparam
    call createGroupStreamAccess(frequency,'Frequency',this%jfcur,maxsubgroup=0,maxdset=maxdataset)
    do icomp = 1,this%ncomp
       do iparam = 1,this%nparam
          call new(dset,1,(/ this%ntot_kernel /),T_COMPLEX)
          call writeDatasetVectorStreamAccess(dset,this%fsa,this%kernel(:,iparam,icomp))
          call addDatasetStreamAccess(frequency,dset)
          call dealloc(dset)
       end do ! iparam
    end do ! icomp
    call addSubgroupStreamAccess(this%root,frequency)
    call dealloc(frequency)
!
    ! remember all frequencies written to file (in their order) in array this%jf by appending the current frequency to array this%jf
    this%jf => reallocate(this%jf,this%nfreq+1)
    this%nfreq = this%nfreq + 1
    this%jf(this%nfreq) = this%jfcur
  end subroutine writeSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief close stream access file
!! \param this spectral waveform kernel
!! \param errmsg error message
!! \param lu return file unit of this waveform kernel file, it can be added to file unit handler
!
  subroutine finalWriteSpectralWaveformKernel(this,errmsg,lu)
    type (spectral_waveform_kernel) :: this
    type (error_message) :: errmsg
    integer, optional :: lu
    ! local
    character(len=32) :: myname = 'finalWriteSpectralWaveformKernel'
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
    allocate(ft(9))
    call new(dset,1,(/ 9 /),T_FLEXIBLE)
    ft(1) = this%parametrization
    ft(2) = this%nparam
    ft(3) = this%ntot_kernel
    ft(4) = this%df
    ft(5) = this%nfreq
    ft(6) = this%kd_id
    ft(7) = this%kgt_id
    ft(8) = this%ncomp
    if(this%on_invgrid) then
       on_invgrid_int = 1
    else
       on_invgrid_int = 0
    end if
    ft(9) = on_invgrid_int
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
    call new(dset,1,(/this%ncomp/),T_FLEXIBLE)
    ! use a do loop in order to assure that character lengths are properly defined
    ! in any case, ft = this%comp is not defined (since "=" is an overlayn operator for type_flexible)
    do icomp = 1,this%ncomp
       ft(icomp) = this%comp(icomp)
    end do ! icomp
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write frequency indices to fourth dataset in group root (if any frequencies)
    if(this%nfreq > 0) then
       call new(dset,1,(/ this%nfreq /),T_INTEGER)
       call writeDatasetVectorStreamAccess(dset,this%fsa,this%jf)
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
  end subroutine finalWriteSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief get meta information from a kernel file
!
  subroutine readMetaInfoSpectralWaveformKernel(filename,lu,errmsg,parametrization,nparam,ntot_kernel,&
         df,nfreq,kd_id,kgt_id,ncomp,on_invgrid,param,comp,jf)
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    character(len=character_length_pmtrz), optional :: parametrization
    integer, optional :: nparam,ntot_kernel,nfreq,ncomp
    real, optional :: df
    character(len=length_ID_kernel_displacement), optional :: kd_id
    character(len=length_ID_kernel_green_tensor), optional :: kgt_id
    logical, optional :: on_invgrid
    character(len=character_length_param), dimension(:), pointer, optional :: param
    character(len=character_length_component), dimension(:), pointer, optional :: comp
    integer, dimension(:), pointer, optional :: jf
    ! local
    character(len=400) :: errstr
    character(len=34) :: myname = 'readMetaInfoSpectralWaveformKernel'
    type (file_stream_access) :: fsa
    type (group_stream_access) :: root
    integer :: ios
    character(len=character_length_pmtrz) :: parametrization_loc
    integer :: nparam_loc,ntot_kernel_loc,nfreq_loc,ncomp_loc
    real :: df_loc
    character(len=length_ID_kernel_displacement) :: kd_id_loc
    character(len=length_ID_kernel_green_tensor) :: kgt_id_loc
    logical :: on_invgrid_loc
    character(len=character_length_param), dimension(:), pointer :: param_loc
    character(len=character_length_component), dimension(:), pointer :: comp_loc
    integer, dimension(:), pointer :: jf_loc
!
    call addTrace(errmsg,myname)
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
    call privateReadMetaInfoSpectralWaveformKernel(fsa,root,parametrization_loc,nparam_loc,ntot_kernel_loc,&
         df_loc,nfreq_loc,kd_id_loc,kgt_id_loc,ncomp_loc,on_invgrid_loc,param_loc,comp_loc,jf_loc,errmsg,myname)
    if(.level.errmsg == 2) goto 1
!
    if(present(parametrization)) parametrization = parametrization_loc
    if(present(nparam)) nparam = nparam_loc
    if(present(ntot_kernel)) ntot_kernel = ntot_kernel_loc
    if(present(df)) df = df_loc
    if(present(nfreq)) nfreq = nfreq_loc
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
    if(present(jf)) then
       jf => jf_loc
       nullify(jf_loc)
    end if
!
1   if(associated(param_loc)) deallocate(param_loc)
    if(associated(comp_loc)) deallocate(comp_loc)
    if(associated(jf_loc)) deallocate(jf_loc)
    call clearGroupTree(root)
    call dealloc(fsa)
  end subroutine readMetaInfoSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief private routine to read all meta information from opened kernel file
!
  subroutine privateReadMetaInfoSpectralWaveformKernel(fsa,root,parametrization,nparam_in_file,ntot_kernel,&
         df,nfreq,kd_id,kgt_id,ncomp_in_file,on_invgrid,param_in_file,comp_in_file,jf_in_file,errmsg,myname)
    type (file_stream_access) :: fsa
    type (group_stream_access) :: root
    character(len=character_length_pmtrz) :: parametrization
    integer:: nparam_in_file,ntot_kernel,nfreq,ncomp_in_file
    real :: df
    character(len=length_ID_kernel_displacement) :: kd_id
    character(len=length_ID_kernel_green_tensor) :: kgt_id
    logical :: on_invgrid
    character(len=character_length_param), dimension(:), pointer :: param_in_file
    character(len=character_length_component), dimension(:), pointer :: comp_in_file
    integer, dimension(:), pointer :: jf_in_file
    type (error_message) :: errmsg
    character(len=*) :: myname
    ! local
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    type (flexible), dimension(:), pointer :: ft
    integer :: on_invgrid_int,iparam_in_file,icomp_in_file,ifreq
    character(len=400) :: errstr
!
    nullify(param_in_file,comp_in_file,jf_in_file)
!
    ! read first info vector
    call traversePathStreamAccess(root,0,(/ 1 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fsa,ft)
    if(.not.associated(ft)) then
       write(errstr,*) "there is no first information vector. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(ft)/=9) then
       write(errstr,*) "size of first information vector is ",size(ft),", must be 9. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    parametrization = ft(1)
    nparam_in_file = ft(2)
    ntot_kernel = ft(3)
    df = ft(4)
    nfreq = ft(5)
    kd_id = ft(6)
    kgt_id = ft(7)
    ncomp_in_file = ft(8)
    on_invgrid_int = ft(9)
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
    if(nfreq < 0) then
       write(errstr,*) "number of frequencies as of information vector = ",nfreq,"; must be >= 0 , file is incosistent"
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
    ! read fourth info vector containing the frequency values (if it was written at all)
    ! if nfreq == 0, it means that there is no data contained in the file, only meta-information contained in file
    ! this module actually permitts to produce such files (although they might by somewhat useless)
    if(nfreq>0) then
       call traversePathStreamAccess(root,0,(/ 4 /),group,dset)
       call readDatasetVectorStreamAccess(dset,fsa,jf_in_file)
       if(.not.associated(jf_in_file)) then
          write(errstr,*) "there is no vector of frequency indices. there is a problem with the file"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(jf_in_file)/=nfreq) then
          write(errstr,*) "information on frequencies inconsistent: number of frequency indices ",size(jf_in_file),&
               " differs from expected number ",nfreq,"; df= ",df!,", frequency indices = ",jf_in_file !! this could cause errstr to overflow (too many characters)
          write(*,*) "ERROR in ",trim(myname),": frequency indices = ",jf_in_file
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(any(jf_in_file<0)) then
          write(errstr,*) "there are ",count(jf_in_file<0)," frequency indices in the kernel file which are < 0"
          call add(errmsg,2,errstr,myname)
          return
       end if
       do ifreq = 2,nfreq
          if(any(jf_in_file(1:ifreq-1)==jf_in_file(ifreq))) then
             write(errstr,*) ifreq,"'th frequency index contained in kernel file occurs more than once. ",&
                  "There must not be duplicated frequency indices"
             call add(errmsg,2,errstr,myname)
             return
          end if
       end do
    end if ! nfreq>0
!
  end subroutine privateReadMetaInfoSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief open stream access file and read all meta information
!
  subroutine initialReadSpectralWaveformKernel(this,filename,lu,errmsg)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    integer :: ios
    character(len=character_length_pmtrz) :: parametrization
    character(len=length_ID_kernel_displacement) :: kd_id
    character(len=length_ID_kernel_green_tensor) :: kgt_id
    integer :: nparam_in_file,ntot_kernel,ncomp_in_file,nfreq,iparam,iparam_in_file,icomp,icomp_in_file
    real :: df
    logical :: on_invgrid
    integer, dimension(:), pointer :: jf_in_file,comp_indx_in_file,param_indx_in_file
    character(len=character_length_param), dimension(:), pointer :: param_in_file
    character(len=character_length_component), dimension(:), pointer :: comp_in_file
    character(len=400) :: errstr
    character(len=33) :: myname = 'initialReadSpectralWaveformKernel'
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
    call privateReadMetaInfoSpectralWaveformKernel(this%fsa,this%root,parametrization,nparam_in_file,ntot_kernel,&
         df,nfreq,kd_id,kgt_id,ncomp_in_file,on_invgrid,param_in_file,comp_in_file,jf_in_file,errmsg,myname)
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
       write(*,*) "ERROR in initialReadSpectralWaveformKernel: model parameters for which this object was initiated = ",&
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
       write(*,*) "ERROR in initialReadSpectralWaveformKernel: receiver components for which this object was initiated = ",&
            this%comp,"; receiver components contained in kernel file = ",comp_in_file
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    deallocate(comp_in_file)
!
    ! if everything went alright, return here after initiating everything and setting flag this%kernel_stat = 2 (file is open to read)
    if(nfreq>0) then
       ! initiate object with the correct frequency content (this is not done in initiateSpectralWaveformKernel above)
       this%df = df
       this%nfreq = nfreq
       this%jf => jf_in_file
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
    ! if there went anything wrong in intial read, close file here and deallocate everything
1   if(associated(comp_indx_in_file)) deallocate(comp_indx_in_file)
    if(associated(param_indx_in_file)) deallocate(param_indx_in_file)
    if(associated(jf_in_file)) deallocate(jf_in_file)
    if(associated(comp_in_file)) deallocate(comp_in_file)
    if(associated(param_in_file)) deallocate(param_in_file)
    call finalReadSpectralWaveformKernel(this)
  end subroutine initialReadSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief read values of waveform kernel for some specific frequency (index jf)
!! \param 
!
  subroutine readSpectralWaveformKernel(this,jf,errmsg)
    type (spectral_waveform_kernel) :: this
    integer :: jf
    type (error_message) :: errmsg
    ! local
    type (group_stream_access), pointer :: group,fgroup
    type (data_stream_access), pointer :: dset
    complex, dimension(:), pointer :: d
    integer :: ifreq,iparam,iparam_in_file,icomp,icomp_in_file,igroup,idset
    logical :: traverse_path_from_root
    character(len=400) :: errstr
    character(len=26) :: myname = 'readSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    if(this%kernel_stat /= 2) then
       call add(errmsg,2,"this waveform kernel object was not yet initiated to read. "//&
            "call initialReadSpectralWaveformKernel first",myname)
       return
    end if
!
    igroup = 0
    do ifreq = 1,this%nfreq
       if(this%jf(ifreq) == jf) then
          igroup = ifreq
          exit
       end if
    end do
    if(igroup == 0) then
       write(errstr,*) "incoming frequency index ",jf,", not contained in kernel file" !,", frequency indices = ",jf !! this could cause errstr to overflow (too many characters)
          write(*,*) "ERROR in readSpectralWaveformKernel: frequency step contained in file = ",this%df,&
               "; frequency indices contained in kernel file = ",this%jf
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
             ! remember frequency group to make calls traversePathStreamAccess faster below
             call traversePathStreamAccess(this%root,1,(/ igroup,idset /),fgroup,dset)
             ! only do this in the very first iteration of this loops, so set flag to false
             traverse_path_from_root = .false.
          else ! traverse_path_from_root
             call traversePathStreamAccess(fgroup,0,(/ idset /),group,dset)
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
       end do ! icomp
    end do ! iparam
!
    this%jfcur = jf
  end subroutine readSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief close stream access file
!! \param this spectral waveform kernel
!! \param lu return file unit of this waveform kernel file, it can be added to file unit handler
!
  subroutine finalReadSpectralWaveformKernel(this,lu)
    type (spectral_waveform_kernel) :: this    
    integer, optional :: lu
    call clearGroupTree(this%root)
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call dealloc(this%fsa)
    this%kernel_stat = 0
  end subroutine finalReadSpectralWaveformKernel
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencySpectralWaveformKernel(this,jf) result(next)
    type (spectral_waveform_kernel) :: this
    integer :: jf
    logical :: next
    integer :: call_count = 0
    save call_count
    !
    if (call_count == this%nfreq) then
       call_count = 0
       next = .false.; return
    endif
    call_count = call_count+1
    jf = this%jf(call_count)
    next = .true.
  end function nextFrequencySpectralWaveformKernel
!
end module spectralWaveformKernel
