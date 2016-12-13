! ===============================================================================
!  GEMINI specific module for generic module kernelGreenTensor
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!-----------------------------------------------------------------------------
!   Provides GEMINI specific interface to routines required by the
!   kernelGreenTensor module to read kernel Green tensors, store them, and later
!   access them when computing waveform kernels.
!-----------------------------------------------------------------------------
module geminiKernelGreenTensor
    use errorMessage
    use fileUnitHandler
    use string
    implicit none
    interface dealloc; module procedure deallocGeminiKernelGreenTensor; end interface
    interface operator (.id.); module procedure getIdGeminiKernelGreenTensor; end interface
    interface operator (.df.); module procedure getDfGeminiKernelGreenTensor; end interface
    interface operator (.nf.); module procedure getNfGeminiKernelGreenTensor; end interface
    interface operator (.jfcur.); module procedure getJfcurGeminiKernelGreenTensor; end interface
    interface operator (.disp.); module procedure getGeminiKernelGreenTensor; end interface
    interface operator (.strain.); module procedure getStrainsGeminiKernelGreenTensor; end interface
    integer, parameter :: length_staname_gkgt = 5 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
    type gemini_kernel_green_tensor
       character (len=length_staname_gkgt) :: id               ! id for kernel Green tensor
       integer :: luf                                          ! file unit for frequency kgt file
       integer :: nf1                                          ! first frequency index
       integer :: nf                                           ! number of frequencies
       real :: df                                              ! frequency stepping
       real :: sigma                                           ! negative imaginary part of angular frequency
       real :: kfcur                                           ! current frequency index
       integer :: ng                                           ! total number of wavefield points
       integer :: nnod                                         ! total number of depth nodes
       character (len=max_length_string) :: kgtbase            ! basename for single-frequency kd files
       complex, dimension(:,:,:), pointer :: g => null()       ! g(wp,comp,force) for current frequency at wp
       complex, dimension(:,:,:), pointer :: gstr => null()    ! gstr(wp,1-6,force) for current frequency at wp
    end type gemini_kernel_green_tensor
!
 contains
!-----------------------------------------------------------------------------------
!  Read basic information needed for all frequencies
!  fuh:          file unit handler
!  filename:     kernel Green tensor meta file
!  errmsg:       error message
!
    subroutine initialReadGeminiKernelGreenTensor(this,fuh,filename,errmsg)
    type (gemini_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    integer :: ierr,lu,ntot
    character (len=34) :: myname = 'initialReadGeminiKernelGreenTensor'
!
    call addTrace(errmsg,myname)
!
! construct kernel Green tensor file base name from name of meta file
! and obtain a file unit from handler for later reading of single-frequency
! kernel Green tensor files  
!
    this%kgtbase = filename
    this%luf = get(fuh)
!
!  open kernel wavefield meta file
!
    lu = get(fuh)
    open(lu,file = this%kgtbase+'.meta',status = 'old',iostat = ierr)
    if (ierr /= 0) then
        call add(errmsg,2,filename+' can not be opened',myname)
        call add(fuh,lu)
        return
    endif
!
!  read two header lines and close file
!
    read(lu,'(a)') this%id
    read(lu,*) this%df,this%nf1,this%nf,this%ng,this%nnod,this%sigma
    close(lu); call add(fuh,lu)
!
!  allocate space for kernels here
!
    ntot = this%ng*this%nnod
    allocate(this%g(ntot,3,3))
    allocate(this%gstr(ntot,6,3))
    end subroutine initialReadGeminiKernelGreenTensor
!----------------------------------------------------------------------------------------------------
!  Deallocate object
!
    subroutine deallocGeminiKernelGreenTensor(this,fuh)
    type (gemini_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    if (associated(this%g)) deallocate(this%g)
    if (associated(this%gstr)) deallocate(this%gstr)
    end subroutine deallocGeminiKernelGreenTensor
!--------------------------------------------------------------------------------------
!  Read Green tensor spectra file for a selected frequency
!  modified to new convention that f = kf*df
!
    subroutine readFrequencyGeminiKernelGreenTensor(this,kf,errmsg)
    type (gemini_kernel_green_tensor) :: this
    integer :: kf
    type (error_message) :: errmsg
    integer :: ic,n,jr,ios
    character (len=3) :: ckf
    character (len=36) :: myname = 'readFrequencyGeminiKernelGreenTensor'
!
    call addTrace(errmsg,myname)
!
    write(ckf,'(i3.3)') kf        ! ASKI uses f = kf*df
    open(this%luf,file = this%kgtbase+'.'+ckf, iostat = ios, form = 'unformatted')
    if (ios /= 0) then
       call add(errmsg,2,'Could not open output file: '+this%kgtbase+'.'+ckf,myname)
       return
    endif
!
    this%kfcur = kf
    do jr = 1,this%nnod
       do n = 1,3
          do ic = 1,3
             read(this%luf) this%g((jr-1)*this%ng+1:jr*this%ng,ic,n)
          enddo
          do ic = 1,6
             read(this%luf) this%gstr((jr-1)*this%ng+1:jr*this%ng,ic,n)
          enddo
       enddo
    enddo
    close(this%luf)
    end subroutine readFrequencyGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get ID
!
    function getIdGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    character(len=length_staname_gkgt) :: res
    res = this%id
    end function getIdGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency spacing df
!
    function getDfGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: res
    res = this%df
    end function getDfGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
    function getNfGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nf
    end function getNfGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency index kfcur
!
    function getJfcurGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%kfcur
    end function getJfcurGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get strains
!
    function getStrainsGeminiKernelGreenTensor(this) result(str)
    type (gemini_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: str
    str => this%gstr
    end function getStrainsGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsGeminiKernelGreenTensor(this,uf_gstr)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: uf_gstr
    ! strains of gemini kernel green tensor are in units (1.0e-3 nm)/km (per N etc...), 
    ! i.e. according to length, they have the dimension-less unit 1.0e-15
    ! This results in very large numbers, which must be scaled by a factor of 1.0e-15 in order to remove this unit
    uf_gstr = 1.0e-15
  end subroutine getUnitFactorStrainsGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get green_tensors
!
    function getGeminiKernelGreenTensor(this) result(u)
    type (gemini_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: u
    u => this%g
    end function getGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of green tensor
!
  subroutine getUnitFactorGeminiKernelGreenTensor(this,uf_g)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: uf_g
    ! gemini kernel green tensor are in units (1.0e-3 nm) (per N etc...), 
    ! i.e. according to length, they have the unit 1.0e-12 m
    ! In order to get SI unit m, green tensor values must be scaled by a factor of 1.0e-12
    uf_g = 1.0e-12
  end subroutine getUnitFactorGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!! modified to new convention that f = kf*df instead of (jf-1)*df before
!! therefore kf = jf-1
!
    function nextFrequencyGeminiKernelGreenTensor(this,kf) result(next)
    type (gemini_kernel_green_tensor) :: this
    integer :: jf,kf
    logical :: next
    integer :: call_count = 0
    save call_count
!
    if (call_count == this%nf) then
        next = .false.
        call_count = 0
        return
    endif
    call_count = call_count+1
    jf = this%nf1+call_count-1
    kf = jf-1
    next = .true.
    end function nextFrequencyGeminiKernelGreenTensor
!
 end module
