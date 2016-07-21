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
!> \brief module which computes integrated spectral waveform sensitivity kernels
!!
!! \details Accessing the generic modules kernelDisplacement, kernelGreenTensor and
!!  kernelReferenceModel, this module computes waveform sensitivity kernels on the 
!!  wavefield points at the given frequencies for parametrizations consistent with 
!!  module modelParametrization and integrates them onto the inversion grid using 
!!  an instance of module integrationWeights.
!!
!! \author Florian Schumacher
!! \date April 2013
!
module spectralWaveformKernel
!
  use modelParametrization
  use kernelReferenceModel
  use kernelDisplacement
  use kernelGreenTensor
  use integrationWeights
  use flexibleType
  use streamAccess
  use realloc
  use mathConstants
  use errorMessage
!
  implicit none
!
  private :: computeIsoLameSpectralWaveformKernel,computeIsoVelocitySpectralWaveformKernel,&
       integrateSpectralWaveformKernel
!
  interface dealloc; module procedure deallocateSpectralWaveformKernel; end interface
  interface getValuesSpectralWaveformKernel
     module procedure getValuesByParamSpectralWaveformKernel
     module procedure getValuesByIndexSpectralWaveformKernel
     module procedure getAllValuesSpectralWaveformKernel
  end interface
  interface operator (.pmtrz.); module procedure getParametrizationSpectralWaveformKernel; end interface
  interface operator (.ntot.); module procedure getNtotInvgridSpectralWaveformKernel; end interface
  interface operator (.df.); module procedure getDfSpectralWaveformKernel; end interface
  interface operator (.nfreq.); module procedure getNfreqSpectralWaveformKernel; end interface
  interface operator (.jf.); module procedure getJfSpectralWaveformKernel; end interface
!
  type spectral_waveform_kernel
     private
     ! model parametrization
     character(len=character_length_pmtrz) :: parametrization = '' !< model parametrization of kernel (also defines number of parameters nparam, below)
     integer :: nparam = 0 !< nparam = numberOfParamModelParametrization(parametrization), size of last dimension of kernel
     integer :: ntot_invgrid = 0 !< total number of inversion grid cells onto which kernels are integrated
     ! spectral discretization
     real :: df = -1. !< df
     integer :: nfreq = 0 !< number of frequencies, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null() !< frequency indices in order of subgroups of kernel file; frequencies compute as f = jf*df [Hz]
     integer :: jfcur = -1 !< indicates current frequency f=jfcur*df for which kernel values are contained in array kernel
     ! kernel values
     complex, dimension(:,:,:), pointer :: kernel => null() !< kernel values for one frequency dim(ntot_invgrid,ncomp=3,nparam)
     ! file handling
     integer :: filestat = 0 !< indicating closed (0), open to write (1), open to read (2)
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
  subroutine initiateSpectralWaveformKernel(this,parametrization,ntot_invgrid,errmsg)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: parametrization
    integer :: ntot_invgrid
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=30) :: myname = 'initiateSpectralWaveformKernel'
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
    this%parametrization = parametrization
    this%nparam = numberOfParamModelParametrization(this%parametrization)
    this%ntot_invgrid = ntot_invgrid
!
    if(associated(this%kernel)) deallocate(this%kernel)
    this%jfcur = -1
    allocate(this%kernel(this%ntot_invgrid,3,this%nparam))
!
    write(errstr,*) "initiated spectral waveform kernel for ",this%ntot_invgrid," inversion grid cells, ",this%nparam,&
         " parameters of parametrization '",trim(this%parametrization),"'"
    call add(errmsg,0,errstr,myname)
  end subroutine initiateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this spectral waveform kernel
!
  subroutine deallocateSpectralWaveformKernel(this)
    type (spectral_waveform_kernel) :: this
    this%parametrization = ''
    this%nparam = 0
    this%ntot_invgrid = 0
    this%df = -1.
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
    if(associated(this%kernel)) deallocate(this%kernel)
    if(this%filestat/=0) then
       call clearGroupTree(this%root)
       call dealloc(this%fsa)
       this%filestat = 0
    end if
  end subroutine deallocateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute spectral waveform kernel at current frequency
!! \details This is a generic routine which, dependent on the parametrization of this kernel
!!  calls specific subroutines to compute the kernel values
!! \param this spectral waveform kernel
!! \param kd kernel displacement (must contain displacement values at same frequency as kgt)
!! \param kgt kernel green tensor (must contain green tensor values at same frequency as kd)
!! \param krm kernel reference model object
!! \param intw integration weights which connect the wavefield points on which kd,kgt are given with inversion grid on which kernel is computed
!! \param errmsg error message
!
  subroutine computeSpectralWaveformKernel(this,kd,kgt,krm,intw,errmsg)
    type (spectral_waveform_kernel) :: this
    type (kernel_displacement) :: kd
    type (kernel_green_tensor) :: kgt
    type (kernel_reference_model) :: krm
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    real :: df_kd,df_kgt
    integer :: jfcur_kd,jfcur_kgt,nwp
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    character(len=400) :: errstr
    character(len=29) :: myname = 'computeSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
    if(.not.associated(this%kernel)) then
       call add(errmsg,2,"spectral_waveform_kernel object not initiated yet",myname)
       return
    end if
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
!
    ! get wavefields and strains
    u => getKernelDisplacement(kd)
    ustr => getStrainsKernelDisplacement(kd)
    g => getKernelGreenTensor(kgt)
    gstr => getStrainsKernelGreenTensor(kgt)
    if(.not.associated(u) .or. .not.associated(ustr) .or. &
         .not.associated(g) .or. .not.associated(gstr)) then
       call add(errmsg,2,'displacement, green tensors or their strains not properly returned from '//&
            'modules kernelDisplacement and kernelGreenTensor',myname)
       return
    end if
    nwp = size(u,1)
    if(nwp /= size(g,1)) then
       write(errstr,*) "number of kernel displacement values ",nwp," differs from number of green tensor values ",size(g,1)
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! communicate this frequency index to computing routines below via my derived data type
    this%jfcur = jfcur_kd
!
    select case(this%parametrization)
       case('isoLame'); call computeIsoLameSpectralWaveformKernel(this,u,ustr,g,gstr,intw,errmsg)
       case('isoVelocity'); call computeIsoVelocitySpectralWaveformKernel(this,krm,u,ustr,g,gstr,intw,errmsg)
       ! ADD YOUR PARAMETRIZATION HERE, AND ADD THE RESPECTIVE SUBROUTINE BELOW
       ! case('your_new_parametrization'); call computeYourNewParametrizationSpectralWaveformKernel(this,u,ustr,g,gstr,intw,errmsg)
    end select
  end subroutine computeSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute kernel values for parametrization 'isoLame'
!! \param this spectral waveform kernel
!! \param u displacement values as returned by kernelDisplacement object
!! \param ustr displacement strains as returned by kernelDisplacement object
!! \param g green tensor values as returned by kernelGreenTensor object
!! \param gstr green tensor strains as returned by kernelGreenTensor object
!! \param intw integration weights as given to computeSpectralWaveformKernel
!! \param errmsg error message
!
  subroutine computeIsoLameSpectralWaveformKernel(this,u,ustr,g,gstr,intw,errmsg)
    type (spectral_waveform_kernel) :: this
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    integer :: icomp
    real :: omega
    complex, dimension(:,:,:), allocatable :: kernel_on_wp
    character(len=36) :: myname = 'computeIsoLameSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    omega = 2.*mc_pi*(this%jfcur*this%df)
!
    allocate(kernel_on_wp(size(u,1),3,this%nparam))
!
    ! compute kernel on wavefield points for first parameter rho
    do icomp = 1,3
       kernel_on_wp(:,icomp,1) = (omega*omega)*(u(:,1)*g(:,1,icomp) + u(:,2)*g(:,2,icomp) + u(:,3)*g(:,3,icomp))
    end do
!
    ! compute kernel on wavefield points for second parameter lambda
    do icomp = 1,3
       kernel_on_wp(:,icomp,2) = -(ustr(:,1)+ustr(:,2)+ustr(:,3))*(gstr(:,1,icomp)+gstr(:,2,icomp)+gstr(:,3,icomp))
    end do
!
    ! compute kernel on wavefield points for third parameter mu
    do icomp = 1,3
       kernel_on_wp(:,icomp,3) = -2.*(ustr(:,1)*gstr(:,1,icomp) + ustr(:,2)*gstr(:,2,icomp) + ustr(:,3)*gstr(:,3,icomp)) &
            -4.*(ustr(:,4)*gstr(:,4,icomp) + ustr(:,5)*gstr(:,5,icomp) + ustr(:,6)*gstr(:,6,icomp))
    end do
!
    ! integrate kernel values onto inversion grid
    call integrateSpectralWaveformKernel(this,kernel_on_wp,intw,errmsg,myname)
!
    deallocate(kernel_on_wp)
  end subroutine computeIsoLameSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief compute kernel values for parametrization 'isoVelocity'
!! \param this spectral waveform kernel
!! \param krm kernel reference model
!! \param u displacement values as returned by kernelDisplacement object
!! \param ustr displacement strains as returned by kernelDisplacement object
!! \param g green tensor values as returned by kernelGreenTensor object
!! \param gstr green tensor strains as returned by kernelGreenTensor object
!! \param intw integration weights as given to computeSpectralWaveformKernel
!! \param errmsg error message
!
  subroutine computeIsoVelocitySpectralWaveformKernel(this,krm,u,ustr,g,gstr,intw,errmsg)
    type (spectral_waveform_kernel) :: this
    type (kernel_reference_model) :: krm
    complex, dimension(:,:), pointer :: ustr,u
    complex, dimension(:,:,:), pointer :: gstr,g
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    integer :: icomp,nwp
    real :: omega
    real, dimension(:), pointer :: rho,vp,vs
    complex, dimension(:,:,:), allocatable :: kernel_on_wp
    complex, dimension(:,:), allocatable :: krho,klambda,kmu
    character(len=400) :: errstr
    character(len=40) :: myname = 'computeIsoVelocitySpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
    nwp = size(u,1)
!
    rho => getModelValuesKernelReferenceModel(krm,this%parametrization,'rho')
    if(.not.associated(rho)) then
       call add(errmsg,2,"no model values for 'isoVelocity'-'rho' in kernel reference model",myname)
       this%jfcur = -1
       return
    end if
    if(size(rho) /= nwp) then
       write(errstr,*) "number of kernel reference model 'rho' values = ",size(rho),&
            " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
       call add(errmsg,2,errstr,myname)
       this%jfcur = -1
       return
    end if
    vp => getModelValuesKernelReferenceModel(krm,this%parametrization,'vp')
    if(.not.associated(vp)) then
       call add(errmsg,2,"no model values for 'isoVelocity'-'vp' in kernel reference model",myname)
       this%jfcur = -1
       return
    end if
    if(size(vp) /= nwp) then
       write(errstr,*) "number of kernel reference model 'vp' values = ",size(vp),&
            " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
       call add(errmsg,2,errstr,myname)
       this%jfcur = -1
       return
    end if
    vs => getModelValuesKernelReferenceModel(krm,this%parametrization,'vs')
    if(.not.associated(vs)) then
       call add(errmsg,2,"no model values for 'isoVelocity'-'vs' in kernel reference model",myname)
       this%jfcur = -1
       return
    end if
    if(size(vs) /= nwp) then
       write(errstr,*) "number of kernel reference model 'vs' values = ",size(vs),&
            " differs from number of wavefield values (displacement, green tensor, strains) = ",nwp
       call add(errmsg,2,errstr,myname)
       this%jfcur = -1
       return
    end if
!
    omega = 2.*mc_pi*(this%jfcur*this%df)
!
    allocate(kernel_on_wp(nwp,3,this%nparam))
!
    ! first compute kernels for parametrization isoLame
!
    allocate(krho(nwp,3),klambda(nwp,3),kmu(nwp,3))
    ! rho (of parametrization isoLame)
    do icomp = 1,3
       krho(:,icomp) = (omega*omega)*(u(:,1)*g(:,1,icomp) + u(:,2)*g(:,2,icomp) + u(:,3)*g(:,3,icomp))
    end do
    ! lambda
    do icomp = 1,3
       klambda(:,icomp) = -(ustr(:,1)+ustr(:,2)+ustr(:,3))*(gstr(:,1,icomp)+gstr(:,2,icomp)+gstr(:,3,icomp))
    end do
    ! mu
    do icomp = 1,3
       kmu(:,icomp) = -2.*(ustr(:,1)*gstr(:,1,icomp) + ustr(:,2)*gstr(:,2,icomp) + ustr(:,3)*gstr(:,3,icomp)) &
            -4.*(ustr(:,4)*gstr(:,4,icomp) + ustr(:,5)*gstr(:,5,icomp) + ustr(:,6)*gstr(:,6,icomp))
    end do
!
    ! then compute isoVelocity kernels via linearized relations between parametrization rho,lambda,mu and rho,vp,vs
    ! (derived e.g. by total derivative)
!
    ! compute kernel on wavefield points for first isoVelocity parameter rho
    do icomp = 1,3
       kernel_on_wp(:,icomp,1) = krho(:,icomp)+(vp(:)**2-2.*vs(:)**2)*klambda(:,icomp)+vs(:)**2*kmu(:,icomp)
    end do
    !
    ! compute kernel on wavefield points for second isoVelocity parameter vp
    do icomp = 1,3
       kernel_on_wp(:,icomp,2) = 2.*rho(:)*vp(:)*klambda(:,icomp)
    end do
    !
    ! compute kernel on wavefield points for third isoVelocity parameter vs
    do icomp = 1,3
       kernel_on_wp(:,icomp,3) = 2.*rho(:)*vs(:)*(kmu(:,icomp)-2.*klambda(:,icomp))
    end do
!
    deallocate(krho,klambda,kmu)
!
    ! integrate kernel values onto inversion grid
    call integrateSpectralWaveformKernel(this,kernel_on_wp,intw,errmsg,myname)
!
    deallocate(kernel_on_wp)
  end subroutine computeIsoVelocitySpectralWaveformKernel
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
    complex, dimension(:,:,:) :: kernel_on_wp
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    character(len=*) :: myname
    ! local
    integer :: icell,iparam,icomp
    integer, dimension(:), pointer :: wp_idx
    real, dimension(:), pointer :: weight
    character(len=400) :: errstr
!
    do icell = 1,this%ntot_invgrid
       wp_idx => intw.wpidx.icell
       weight => intw.weight.icell
       if(.not.associated(wp_idx)) then
          write(errstr,*) "there are no integration weights for inversion grid cell ",icell,&
               ", inversion grid cell index seems to be out of range"
          call add(errmsg,2,errstr,myname)
          this%jfcur = -1
          return
       end if
       do iparam = 1,this%nparam
          do icomp = 1,3
             this%kernel(icell,icomp,iparam) = sum(weight*kernel_on_wp(wp_idx,icomp,iparam))
          end do ! icomp
       end do ! iparam
    end do ! icell
  end subroutine integrateSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given parameter
!
  function getValuesByParamSpectralWaveformKernel(this,param) result(p)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: param
    complex, dimension(:,:), pointer :: p
    nullify(p)
    if(.not.associated(this%kernel)) return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    p => this%kernel(:,:,indexOfParamModelParametrization(this%parametrization,param))
  end function getValuesByParamSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief return kernel values for given index of parameter
!
  function getValuesByIndexSpectralWaveformKernel(this,i) result(p)
    type (spectral_waveform_kernel) :: this
    integer :: i
    complex, dimension(:,:), pointer :: p
    nullify(p)
    if(.not.associated(this%kernel)) return
    if(i<1 .or. i>this%nparam) return
    p => this%kernel(:,:,i)
  end function getValuesByIndexSpectralWaveformKernel
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
!> \brief return number of inversion grid cells of kernel
!
  function getNtotInvgridSpectralWaveformKernel(this) result(ntot)
    type (spectral_waveform_kernel), intent(in) :: this
    integer :: ntot
    ntot = this%ntot_invgrid
  end function getNtotInvgridSpectralWaveformKernel
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
    ios = createFileStreamAccess(this%fsa,lu,filename)
    if(ios/=0) then
       write(errstr,*) "could not open file '"//trim(filename)//"' to write, raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%filestat = 1
!
    if(present(nfreq)) then
       ! it is more efficient to create root like this, knowing the number of subgroups to come
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=nfreq,maxdset=2)
    else
       ! otherwise assume enough subgroups, which is still more efficient than creating root without subrgoup information
       call createGroupStreamAccess(this%root,'Root',0,maxsubgroup=150,maxdset=2)
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
  subroutine writeSpectralWaveformKernel(this)
    type (spectral_waveform_kernel) :: this
    ! local
    integer :: iparam,icomp,maxdataset
    type (group_stream_access) :: frequency
    type (data_stream_access) :: dset
!
    ! if file is not open to write, return
    if(this%filestat /= 1) return
!
    ! if there is an indicator that there are no sensible kernel values stored, just return doing nothing
    if(this%jfcur == -1) return
!
    maxdataset = 3*this%nparam
    call createGroupStreamAccess(frequency,'Frequency',this%jfcur,maxsubgroup=0,maxdset=maxdataset)
    do iparam = 1,this%nparam
       do icomp = 1,3
          call new(dset,1,(/ this%ntot_invgrid /),T_COMPLEX)
          call writeDatasetVectorStreamAccess(dset,this%fsa,this%kernel(:,icomp,iparam))
          call addDatasetStreamAccess(frequency,dset)
          call dealloc(dset)
       end do ! icomp
    end do ! iparam
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
!
    ! write general information to first dataset in group root
    allocate(ft(5))
    call new(dset,1,(/ 5 /),T_FLEXIBLE)
    ft(1) = this%parametrization
    ft(2) = this%nparam
    ft(3) = this%ntot_invgrid
    ft(4) = this%df
    ft(5) = this%nfreq
    call writeDatasetVectorStreamAccess(dset,this%fsa,ft)
    call addDatasetStreamAccess(this%root,dset)
    deallocate(ft); call dealloc(dset)
!
    ! write frequency indices to second dataset in group root
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
    this%filestat = 0
  end subroutine finalWriteSpectralWaveformKernel
!------------------------------------------------------------------------
!> \brief open stream access file to read kernel values, one frequency at a time
!! \param 
!
  subroutine initialReadSpectralWaveformKernel(this,filename,lu,errmsg)
    type (spectral_waveform_kernel) :: this
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    integer :: ios
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    type (flexible), dimension(:), pointer :: ft
    character(len=character_length_pmtrz) :: parametrization
    integer :: nparam,ntot_invgrid,nfreq
    real :: df
    integer, dimension(:), pointer :: jf
    character(len=400) :: errstr
    character(len=33) :: myname = 'initialReadSpectralWaveformKernel'
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
    if(size(ft)/=5) then
       write(errstr,*) "size of first information vector is ",size(ft),", must be 5. there is a problem with the file"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    parametrization = ft(1)
    nparam = ft(2)
    ntot_invgrid = ft(3)
    df = ft(4)
    nfreq = ft(5)
    deallocate(ft)
    if(numberOfParamModelParametrization(parametrization) /= nparam) then
       write(errstr,*) "number of parameters ",nparam,", contained in this file is inconsistent with its parametrization '"&
            //trim(parametrization)//"' which has ",numberOfParamModelParametrization(parametrization),&
            " parameters. This could indicate that module modelParametrization has been modified since createion of the file"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    ! read second info vector containing the frequency values (if it was written at all)
    ! if nfreq == 0, it means that there is no data contained in the file
    if(nfreq>0) then
       call traversePathStreamAccess(this%root,0,(/ 2 /),group,dset)
       call readDatasetVectorStreamAccess(dset,this%fsa,jf)
       if(size(jf)/=nfreq) then
          write(errstr,*) "information on frequencies inconsistent: number of frequency indices ",size(jf),&
               " differs from expected number ",nfreq,"; df= ",df!,", frequency indices = ",jf !! this could cause errstr to overflow (too many characters)
          write(*,*) "ERROR in initialReadSpectralWaveformKernel: frequency indices = ",this%jf
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       this%df = df
       this%nfreq = nfreq
       this%jf => jf
    end if
!
    ! finally initiate spectralWaveformKernel object
    call initiateSpectralWaveformKernel(this,parametrization,ntot_invgrid,errmsg)
    if(.level.errmsg ==2) goto 1
!
    ! if everything went alright, return here, after setting flag this%filestat = 2 (file is open to read)
    this%filestat = 2
    return
!
    ! if there went anything wrong in intial read, close file here and revert all initiation that's been done so far
1   call deallocateSpectralWaveformKernel(this)
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
    integer :: ifreq,iparam,icomp,igroup,idset
    character(len=400) :: errstr
    character(len=26) :: myname = 'readSpectralWaveformKernel'
!
    call addTrace(errmsg,myname)
!
    if(this%filestat /= 2) then
       call add(errmsg,2,"there was no file opened to read yet, call initialReadSpectralWaveformKernel first",myname)
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
       write(errstr,*) "incoming frequency index ",jf,", not contained in file. frequencies contained in file are: df = ",&
            this%df!,", frequency indices = ",jf !! this could cause errstr to overflow (too many characters)
          write(*,*) "ERROR in readSpectralWaveformKernel: frequency indices = ",this%jf
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    idset = 0
    do iparam = 1,this%nparam
       do icomp = 1,3
          idset = idset + 1
          if(idset == 1) then
             ! remember frequency group to make calls traversePathStreamAccess faster below
             call traversePathStreamAccess(this%root,1,(/ igroup,idset /),fgroup,dset)
          else
             call traversePathStreamAccess(fgroup,0,(/ idset /),group,dset)
          end if
          call readDatasetVectorStreamAccess(dset,this%fsa,d)
          if(size(d) /= this%ntot_invgrid) then
             write(errstr,*) "number of kernel values ",size(d)," for component ",icomp," of ",iparam,&
                  "'th parameter '",trim(getParamFromIndexModelParametrization(this%parametrization,iparam)),&
                  "' differs from number of inversion grid cells ",this%ntot_invgrid,", i.e. kernel file may be corrupt"
             call add(errmsg,2,errstr,myname)
             return
          end if
          this%kernel(:,icomp,iparam) = d
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
    if(present(lu)) lu = getFileUnitStreamAccess(this%fsa)
    call deallocateSpectralWaveformKernel(this)
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
