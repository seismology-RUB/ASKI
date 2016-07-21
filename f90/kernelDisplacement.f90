!----------------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Wrapper module for kernel displacement
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  The interfaces to method GEMINI still need to be adapted to the new ASKI version,
!!  so far only method SPECFEM3D is supported. 
!!  Displacement and strain arrays ordered consistently with ordering of wavefield points.
!
module kernelDisplacement
!!$  use geminiKernelDisplacement
  use specfem3dKernelDisplacement
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateKernelDisplacement; end interface
  interface operator (.id.); module procedure getIdKernelDisplacement; end interface
  interface operator (.df.); module procedure getDfKernelDisplacement; end interface
  interface operator (.nf.); module procedure getNfKernelDisplacement; end interface
  interface operator (.jfcur.); module procedure getJfcurKernelDisplacement; end interface
  interface operator (.disp.); module procedure getKernelDisplacement; end interface
  interface operator (.strain.); module procedure getStrainsKernelDisplacement; end interface
!< fork here to specific method by associating some pointer (only ONE!)
  type kernel_displacement
     private
!!$     type (gemini_kernel_displacement), pointer :: gemini_kd => null()         ! gemini specific kernel displacement object
     type (specfem3d_kernel_displacement), pointer :: specfem_kd => null()       ! Specfem specific kernel displacement object
  end type kernel_displacement
!
  integer, parameter :: length_ID_kernel_displacement = 13 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
!
contains
!------------------------------------------------------------------------
!> \brief Initiate kernel displacement
!
  subroutine initiateKernelDisplacement(this,method,fuh,filename,errmsg)
    type (kernel_displacement) :: this
    character(len=*) :: method
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=26) :: myname = 'initiateKernelDisplacement'
!
    call addTrace(errmsg,myname)
    select case (method)
!!$    case('GEMINI')
!!$       allocate(this%gemini_kd)
!!$       errmsg = initialReadGeminiKernelDisplacement(this%gemini_kd,fuh,filename)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
    case('SPECFEM3D')
       allocate(this%specfem_kd)
       call initialReadSpecfem3dKernelDisplacement(this%specfem_kd,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case default
       call add(errmsg,2,"Invalid forward computation method '"//trim(method)//"'",myname)
       return
    end select
  end subroutine initiateKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Deallocate complete kernel_displacement object
!
  subroutine deallocateKernelDisplacement(this,fuh)
    type (kernel_displacement) :: this
    type (file_unit_handler) :: fuh
!
!!$    if (associated(this%gemini_kd)) then
!!$       call deallocGeminiKernelDisplacement(this%gemini_kd,fuh)
!!$    else if (associated(this%specfem_kd)) then
       call dealloc(this%specfem_kd,fuh)
!!$    endif
  end subroutine deallocateKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Read in kernel displacement for one frequency
!
  subroutine readFrequencyKernelDisplacement(this,jf,errmsg)
    type (kernel_displacement) :: this
    integer :: jf
    type (error_message) :: errmsg
    character (len=31) :: myname = 'readFrequencyKernelDisplacement'
!
    call addTrace(errmsg,myname)
!!$    if (associated(this%gemini_kd)) then
!!$       errmsg = readFrequencyGeminiKernelDisplacement(this%gemini_kd,jf)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
!!$    else if (associated(this%specfem_kd)) then
       call readFrequencySpecfem3dKernelDisplacement(this%specfem_kd,jf,errmsg)
       if (.level.errmsg == 2) return
!!$    endif
  end subroutine readFrequencyKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get ID
!
  function getIdKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    character(len=length_ID_kernel_displacement) :: res
!
    res = ''
!!$    if (associated(this%gemini_kd)) then
!!$       res = .id.(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       res = .id.(this%specfem_kd)
!!$    endif
  end function getIdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    real :: res
!
    res = 0.
!!$    if (associated(this%gemini_kd)) then
!!$       res = .df.(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       res = .df.(this%specfem_kd)
!!$    endif
  end function getDfKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
  function getNfKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    integer :: res
!
    res = 0
!!$    if (associated(this%gemini_kd)) then
!!$       res = .nf.(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       res = .nf.(this%specfem_kd)
!!$    endif
  end function getNfKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    integer :: res
!
    res = 0
!!$    if (associated(this%gemini_kd)) then
!!$       res = .jfcur.(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       res = .jfcur.(this%specfem_kd)
!!$    endif
  end function getJfcurKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get pointer to strains
!
  function getStrainsKernelDisplacement(this) result(ustr)
    type (kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: ustr
!!$    if (associated(this%gemini_kd)) then
!!$       ustr => getStrainsGeminiKernelDisplacement(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       ustr => .strain.(this%specfem_kd)
!!$    endif
  end function getStrainsKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get pointer to displacements
!
  function getKernelDisplacement(this) result(u)
    type (kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: u
!!$    if (associated(this%gemini_kd)) then
!!$       u => getGeminiKernelDisplacement(this%gemini_kd)
!!$    else if (associated(this%specfem_kd)) then
       u => .disp.(this%specfem_kd)
!!$    endif
  end function getKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyKernelDisplacement(this,jf) result(next)
    type (kernel_displacement) :: this
    integer :: jf
    logical :: next
!!$    if (associated(this%gemini_kd)) then
!!$       next = nextFrequencyGeminiKernelDisplacement(this%gemini_kd,jf)
!!$    else if (associated(this%specfem_kd)) then
       next = nextFrequencySpecfem3dKernelDisplacement(this%specfem_kd,jf)
!!$    endif
  end function nextFrequencyKernelDisplacement
!
end module kernelDisplacement
