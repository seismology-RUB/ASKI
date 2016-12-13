! ===============================================================================
!  GEMINI specific module for generic module kernelDisplacement.f90
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
!   kernelDisplacement module to read kernel wavefields, store them, and later
!   access them when computing waveform kernels.
!-----------------------------------------------------------------------------
 module geminiKernelDisplacement
    use errorMessage
    use fileUnitHandler
    use string
    implicit none
    interface dealloc; module procedure deallocGeminiKernelDisplacement; end interface
    interface operator (.id.); module procedure getIdGeminiKernelDisplacement; end interface
    interface operator (.nf.); module procedure getNfGeminiKernelDisplacement; end interface
    interface operator (.df.); module procedure getDfGeminiKernelDisplacement; end interface
    interface operator (.jfcur.); module procedure getJfcurGeminiKernelDisplacement; end interface
    interface operator (.disp.); module procedure getGeminiKernelDisplacement; end interface
    interface operator (.strain.); module procedure getStrainsGeminiKernelDisplacement; end interface
    integer, parameter :: length_evid_gkd = 13 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
    type gemini_kernel_displacement
       private
       character (len=length_evid_gkd) :: id                   ! event id for kernel displacement
       integer :: luf                                          ! file unit for frequency kd file
       real :: df                                              ! frequency stepping
       real :: kfcur                                           ! current frequency index: f = kf*df
       integer :: nf1                                          ! first freqeuncy index
       integer :: nf                                           ! total number of frequencies
       real :: sigma                                           ! negative imaginary part of angular frequency
       integer :: ng                                           ! total number of surface wavefield points
       integer :: nnod                                         ! nuber of depth nodes
       character (len=max_length_string) :: kdbase             ! basename for single-frequency kd files
       complex, dimension(:,:), pointer :: u => null()         ! 2D array for displacement on grid
       complex, dimension(:,:), pointer :: ustr => null()      ! 2D array for strain components on grid
    end type gemini_kernel_displacement
!
 contains
!-----------------------------------------------------------------------------------
!  Read basic information needed for all frequencies
!  fuh:          file unit handler
!  filename:     kernel wavefield meta file
!  errmsg:       error message
!
    subroutine initialReadGeminiKernelDisplacement(this,fuh,filename,errmsg)
    type (gemini_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    integer :: ierr,lu,ntot
    character (len=35) :: myname = 'initialReadGeminiKernelDisplacement'
!
    call addTrace(errmsg,myname)
!
! construct kernel displacement file base name from name of meta file
! and obtain file unit form handler for later reading of single-frequency
! kernel displacement files  
!
    this%kdbase = filename
    this%luf = get(fuh)
!
!  open kernel wavefield meta file
!
    lu = get(fuh)
    open(lu,file = this%kdbase+'.meta',status = 'old',iostat = ierr)
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
!  allocate space for displacement and strain
!
    ntot = this%ng*this%nnod
    allocate(this%u(ntot,3))
    allocate(this%ustr(ntot,6))
    end subroutine initialReadGeminiKernelDisplacement
!----------------------------------------------------------------------------------------------------
!  Deallocate object
!
    subroutine deallocGeminiKernelDisplacement(this,fuh)
    type (gemini_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    if (associated(this%u)) deallocate(this%u)
    if (associated(this%ustr)) deallocate(this%ustr)
    end subroutine deallocGeminiKernelDisplacement
!-----------------------------------------------------------------------------------------------------
!  Read displacement spectra file for different frequencies
!  modified to new convention that f = kf*df
!
    subroutine readFrequencyGeminiKernelDisplacement(this,kf,errmsg)
    type (gemini_kernel_displacement) :: this
    integer :: kf
    type (error_message) :: errmsg
    integer :: ic,jr,ios
    character (len=3) :: ckf
    character (len=37) :: myname = 'readFrequencyGeminiKernelDisplacement'
!
    call addTrace(errmsg,myname)
!
    write(ckf,'(i3.3)') kf        ! ASKI uses f = kf*df
    open(this%luf,file = this%kdbase+'.'+ckf, iostat = ios, form = 'unformatted')
    if (ios /= 0) then
       call add(errmsg,2,'Could not open output file: '+this%kdbase+'.'+ckf,myname)
       return
    endif
!
    this%kfcur = kf
    do jr = 1,this%nnod
       do ic = 1,3
          read(this%luf) this%u((jr-1)*this%ng+1:jr*this%ng,ic)
       enddo
       do ic = 1,6
          read(this%luf) this%ustr((jr-1)*this%ng+1:jr*this%ng,ic)
       enddo
    enddo
    close(this%luf)
    end subroutine readFrequencyGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get ID
!
    function getIdGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    character(len=length_evid_gkd) :: res
    res = this%id
    end function getIdGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency spacing df
!
    function getDfGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: res
    res = this%df
    end function getDfGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
    function getNfGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nf
    end function getNfGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index kfcur
!
    function getJfcurGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%kfcur
    end function getJfcurGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get strains
!
    function getStrainsGeminiKernelDisplacement(this) result(str)
    type (gemini_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: str
    str => this%ustr
    end function getStrainsGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsGeminiKernelDisplacement(this,uf_ustr)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: uf_ustr
    ! strains of gemini kernel displacement are in units nm/km (per N or per Nm etc...), 
    ! i.e. according to length, they have the dimension-less unit 1.0e-12
    ! This results in very large numbers, which must be scaled by a factor of 1.0e-12 in order to remove this unit
    uf_ustr = 1.0e-12 
  end subroutine getUnitFactorStrainsGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get displacements
!
    function getGeminiKernelDisplacement(this) result(u)
    type (gemini_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: u
    u => this%u
    end function getGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of displacements
!
  subroutine getUnitFactorGeminiKernelDisplacement(this,uf_u)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: uf_u
    ! gemini kernel displacements are in units nm (per N or per Nm etc...), i.e. according to length, 
    ! they have the unit 1.0e-9 m
    ! In order to get SI unit m, displacement values must be scaled by a factor of 1.0e-9
    uf_u = 1.0e-9
  end subroutine getUnitFactorGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!!  modified for new convention that f = kf*df instead of f = (jf-1)*df
!!  hence kf = jf-1 or jf = kf+1
!
    function nextFrequencyGeminiKernelDisplacement(this,kf) result(next)
    type (gemini_kernel_displacement) :: this
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
    end function nextFrequencyGeminiKernelDisplacement
!
 end module
