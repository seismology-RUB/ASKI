! ==============================================================================
!  Earth model module for use with kernel computation
! ==============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   ASKI version 1.2 is distributed in the hope that it
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!------------------------------------------------------------------------------
!   Deals with earth model information at external radial nodes.
!   Needed for kernel computation.
!-----------------------------------------------------------------------------
module geminiKernelReferenceModel
    use modelParametrization
    use errorMessage
    use locatePoint
    use string
    implicit none
    interface dealloc; module procedure deallocGeminikernelReferenceModel; end interface
    type gemini_kernel_reference_model
       character (len = max_length_string) :: name              !< name of earth model
       character (len = max_length_string) :: attmode           !< attenuation mode string
       integer :: aniflag                                       !< isotropic (0), transversely isotropic (1)
       double precision :: fref                                 !< reference frequency
       double precision :: rearth                               !< earth radius
       integer :: nnod                                          !< number of nodes
       integer :: ng                                            !< number of lateral grid nodes
       real, dimension(:), allocatable :: rnod                  !< radial nodes in kilometers
       complex, dimension(:,:), allocatable :: zelcon           !< complex A,C,F,L,N,Kappa,Mue
       real, dimension(:), allocatable :: rho                   !< density
       real, dimension(:), allocatable :: qkinv                 !< inverse Qkappa, qk = 1/Qk
       real, dimension(:), allocatable :: qminv                 !< inverse Qmu, qm = 1/Qm
    end type gemini_kernel_reference_model
!
 contains
!----------------------------------------------------------------------------
!  Read gemini earth model from ASCII file (used by kernelReferenceModel)
!
    subroutine readGeminiKernelReferenceModel(this,lu,filename,errmsg)
    type (gemini_kernel_reference_model) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=20) :: myname = 'readGeminiKernelReferenceModel'
    integer :: ierr
!
    call addTrace(errmsg,myname)
    call readWithoutNgGeminiKernelReferenceModel(this,lu,filename,errmsg,.false.)  ! do not close file
!
!  read number of lateral grid nodes
!    
    read(lu,*,iostat = ierr) this%ng
    if (ierr /= 0) then
        call add(errmsg,2,'gemini earth model files does not have line with ng',myname)
        return
    endif
    close(lu)
    end subroutine readGeminiKernelReferenceModel
!----------------------------------------------------------------------------
!  Read gemini earth model from ASCII file without ng
!  Used by computeKernelWavefield and computeKernelGreenTensor
!
    subroutine readWithoutNgGeminiKernelReferenceModel(this,lu,filename,errmsg,closeflag)
    type (gemini_kernel_reference_model) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    logical :: closeflag
    character (len=20) :: myname = 'readGeminiKernelReferenceModel'
    integer :: i,j,ierr,nnod
    real, dimension(14) :: zelcon_tmp
!
    call addTrace(errmsg,myname)
    open(lu,file = filename, status = 'old', iostat = ierr)
    if (ierr /= 0) then
        call add(errmsg,2,'cannot open '+filename,myname)
        return
    endif
    read(lu,'(a)') this%name
    read(lu,'(a)') this%attmode
    read(lu,*) this%aniflag,this%fref,this%rearth
    read(lu,*) nnod
    this%nnod = nnod
    this%ng = 0                      ! ng not known yet
!
!  allocate space for parameters and nodes and read
!    
    allocate(this%rnod(nnod),this%rho(nnod),this%qkinv(nnod),this%qminv(nnod))
    allocate(this%zelcon(nnod,7))
    do i = 1,nnod
       read(lu,*) this%rnod(i),this%rho(i),this%qkinv(i),this%qminv(i)
       read(lu,*) zelcon_tmp(:)
       do j = 1,7
          this%zelcon(i,j) = cmplx(zelcon_tmp(2*j-1),zelcon_tmp(2*j))
       end do
    enddo
    if (closeflag) close(lu)
    end subroutine readWithoutNgGeminiKernelReferenceModel
!--------------------------------------------------------------------------------------
!  Write Gemini earth model adding final line with number of lateral grid nodes
!  Used by computeKernelWavefield and computeKernelGreenTensor
!
    subroutine writeGeminiKernelReferenceModel(this,lu,filename,ng,errmsg)
    type (gemini_kernel_reference_model) :: this
    integer :: lu
    character (len=*) :: filename
    integer :: ng
    type (error_message) :: errmsg
    integer :: i,j,ierr
    character (len=21) :: myname = 'writeGeminiKernelReferenceModel'
!
    call addTrace(errmsg,myname)
    open(lu,file = filename,iostat = ierr)
    if (ierr /= 0) then
        call add(errmsg,2,'cannot open '+filename,myname)
        return
    endif
    write(lu,'(a)') trim(this%name)
    write(lu,'(a)') trim(this%attmode)
    write(lu,'(i4,2e15.6)') this%aniflag,this%fref,this%rearth
    write(lu,'(i5)') this%nnod
    do i = 1,this%nnod
       write(lu,'(4e15.6)') this%rnod(i),this%rho(i),this%qkinv(i),this%qminv(i)
       write(lu,'(14e15.6)') (this%zelcon(i,j),j = 1,7)
    enddo
    write(lu,'(i7)') ng
    close(lu)
    end subroutine writeGeminiKernelReferenceModel
!--------------------------------------------------------------------------------------------
! Deallocate gemini_kernel_reference_model
!
    subroutine deallocGeminiKernelReferenceModel(this)
    type (gemini_kernel_reference_model) :: this
    if (allocated(this%rnod)) deallocate(this%rnod)
    if (allocated(this%rho)) deallocate(this%rho)
    if (allocated(this%qkinv)) deallocate(this%qkinv)
    if (allocated(this%qminv)) deallocate(this%qminv)
    if (allocated(this%zelcon)) deallocate(this%zelcon)
    end subroutine deallocGeminiKernelReferenceModel
!----------------------------------------------------------------------------
!> \brief Get index of node closest to given depth in km
!
    function getNodeIndexFromDepthGeminiKernelReferenceModel(this,zs) result(jr)
    type (gemini_kernel_reference_model) :: this
    real :: zs,rs
    integer :: jr
!
    rs = this%rearth-zs
    jr = locate(rs,this%nnod,this%rnod)              ! index of node below rs
    if (jr == 0) jr = 1                              ! use deepest node if rs is below it
    if (jr < this%nnod) then                         ! take node closest to rs
       if (abs(rs-this%rnod(jr)) > abs(rs-this%rnod(jr+1))) jr = jr+1
    endif
    end function getNodeIndexFromDepthGeminiKernelReferenceModel
!------------------------------------------------------------------------
!> \brief Get model values for some particular parameter of some particular parametrization
!! \param pmtrz some model parametrization, as defined in module modelParametrization (e.g. "isoVelocity1000")
!! \param param model parameter, which must be one of parametrization pmtrz
!! \param model_values pointer to array of model values on wavefield points. Nullified if there is a problem
!
    function getModelValuesWPGeminiKernelReferenceModel(this,pmtrz,param) result(model_values)
    type (gemini_kernel_reference_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
!
    nullify(model_values)
    if (.not.validModelParametrization(pmtrz)) return
    if (.not.validParamModelParametrization(pmtrz,param)) return
    select case(pmtrz)
       case('isoVelocity1000')
          select case(param)
             case('rho'); model_values => getDensityWPGeminiKernelReferenceModel(this) ! [g/cm^3]
             case('vp'); model_values => getPVelocityWPGeminiKernelReferenceModel(this) ! [km/s]
             case('vs'); model_values => getSVelocityWPGeminiKernelReferenceModel(this) ! [km/s]
          end select
       case('isoVelocitySI')
          select case(param)
             case('rho'); model_values => getDensityWPGeminiKernelReferenceModel(this) ! [g/cm^3]
             case('vp'); model_values => getPVelocityWPGeminiKernelReferenceModel(this) ! [km/s] 
             case('vs'); model_values => getSVelocityWPGeminiKernelReferenceModel(this) ! [km/s]
          end select
          model_values = model_values * 1.0e3 ! account for unit factor 1000.0 to get Kg/m^3 and m/s
       case('isoLameSI')
          select case(param)
             case('rho')
                model_values => getDensityWPGeminiKernelReferenceModel(this)     ! [g/cm^3]
                model_values = model_values * 1.0e3                    ! [kg/m^3]                 
             case('lambda')
                model_values => getLambdaWPGeminiKernelReferenceModel(this)    ! [GPa]
                model_values = model_values * 1.0e9                          ! [Pa]                 
             case('mu')
                model_values => getMuWPGeminiKernelReferenceModel(this)         ! [GPa]
                model_values = model_values * 1.0e9                   ! [Pa]                 
          end select
    end select
    end function getModelValuesWPGeminiKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get density of earth model at wavefield points
!
    function getDensityWPGeminiKernelReferenceModel(this) result(rho)
    type (gemini_kernel_reference_model) :: this
    real, dimension(:), pointer :: rho
    integer :: j
!
    allocate(rho(this%ng*this%nnod))
    do j = 1,this%nnod
        rho((j-1)*this%ng+1:j*this%ng) = this%rho(j)
    enddo
    end function getDensityWPGeminiKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get P velocity of earth model at wavefield points
!
    function getPVelocityWPGeminiKernelReferenceModel(this) result(vp)
    type (gemini_kernel_reference_model) :: this
    real, dimension(:), pointer :: vp
    integer :: j
!
    allocate(vp(this%ng*this%nnod))
    do j = 1,this%nnod
        vp((j-1)*this%ng+1:j*this%ng) = sqrt( real(this%zelcon(j,6)+4.d0*this%zelcon(j,7)/3.d0)/this%rho(j) )
    enddo
    end function getPVelocityWPGeminiKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get S velocity of earth model at wavefield points
!
    function getSVelocityWPGeminiKernelReferenceModel(this) result(vs)
    type (gemini_kernel_reference_model) :: this
    real, dimension(:), pointer :: vs
    integer :: j
!
    allocate(vs(this%ng*this%nnod))
    do j = 1,this%nnod
        vs((j-1)*this%ng+1:j*this%ng) = sqrt( real(this%zelcon(j,7))/this%rho(j) )
    enddo
    end function getSVelocityWPGeminiKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get lambda of earth model at wavefield points
!
    function getLambdaWPGeminiKernelReferenceModel(this) result(lambda)
    type (gemini_kernel_reference_model) :: this
    real, dimension(:), pointer :: lambda
    integer :: j
!
    allocate(lambda(this%ng*this%nnod))
    do j = 1,this%nnod
        lambda((j-1)*this%ng+1:j*this%ng) = real(this%zelcon(j,6)-2.d0*this%zelcon(j,7)/3.d0)
     enddo
     end function getLambdaWPGeminiKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get mu of earth model at wavefield points
!
    function getMuWPGeminiKernelReferenceModel(this) result(mu)
    type (gemini_kernel_reference_model) :: this
    real, dimension(:), pointer :: mu
    integer :: j
!
    allocate(mu(this%ng*this%nnod))
    do j = 1,this%nnod
        mu((j-1)*this%ng+1:j*this%ng) = real(this%zelcon(j,7))
     enddo
     end function getMuWPGeminiKernelReferenceModel
!-----------------------------------------------------------------------------
 end module
