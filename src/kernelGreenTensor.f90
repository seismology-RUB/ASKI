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
!> \brief Wrapper module for kernel Green tensor
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  The interfaces to method GEMINI still need to be adapted to the new ASKI version,
!!  so far only method SPECFEM3D is supported. 
!!  Green functions and strain arrays ordered consistently with ordering of wavefield points 
!
module kernelGreenTensor
!!$  use geminiKernelGreenTensor
  use specfem3dKernelGreenTensor
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateKernelGreenTensor; end interface
  interface operator (.id.); module procedure getIdKernelGreenTensor; end interface
  interface operator (.df.); module procedure getDfKernelGreenTensor; end interface
  interface operator (.nf.); module procedure getNfKernelGreenTensor; end interface
  interface operator (.jfcur.); module procedure getJfcurKernelGreenTensor; end interface
  interface operator (.disp.); module procedure getKernelGreenTensor; end interface
  interface operator (.strain.); module procedure getStrainsKernelGreenTensor; end interface
!< fork here to specific method by associating some pointer (only ONE!)
  type kernel_green_tensor
     private
!!$     type (gemini_kernel_green_tensor), pointer :: gemini_kgt => null()         !< gemini specific kernel displacement object
     type (specfem3d_kernel_green_tensor), pointer :: specfem_kgt => null()       !< specfem3d specific kernel displacement object
  end type kernel_green_tensor
!
  integer, parameter :: length_ID_kernel_green_tensor = 13 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
!
contains
!------------------------------------------------------------------------
!> \brief Initiate kernel Green tensor
!
  subroutine initiateKernelGreenTensor(this,method,fuh,filename,errmsg)
    type (kernel_green_tensor) :: this
    character(len=*) :: method
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=25) :: myname = 'initiateKernelGreenTensor'
!
    call addTrace(errmsg,myname)
    select case (method)
!!$    case('GEMINI')
!!$       allocate(this%gemini_kgt)
!!$       errmsg = initialReadGeminiKernelGreenTensor(this%gemini_kgt,fuh,filename)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
    case('SPECFEM3D')
       allocate(this%specfem_kgt)
       call initialReadSpecfem3dKernelGreenTensor(this%specfem_kgt,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case default
       call add(errmsg,2,"Invalid forward computation method '"//trim(method)//"'",myname)
       return
    end select
  end subroutine initiateKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Deallocate complete kernel_green_tensor object
!
  subroutine deallocateKernelGreenTensor(this,fuh)
    type (kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
!
!!$    if (associated(this%gemini_kgt)) then
!!$       call deallocGeminiKernelGreenTensor(this%gemini_kgt,fuh)
!!$    else if (associated(this%specfem_kgt)) then
       call dealloc(this%specfem_kgt,fuh)
!!$    endif
  end subroutine deallocateKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Read in kernel displacement for one frequency
!
  subroutine readFrequencyKernelGreenTensor(this,jf,errmsg)
    type (kernel_green_tensor) :: this
    integer :: jf
    type (error_message) :: errmsg
    character (len=30) :: myname = 'readFrequencyKernelGreenTensor'
!
    call addTrace(errmsg,myname)
!!$    if (associated(this%gemini_kgt)) then
!!$       errmsg = readFrequencyGeminiKernelGreenTensor(this%gemini_kgt,jf)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
!!$    else if (associated(this%specfem_kgt)) then
       call readFrequencySpecfem3dKernelGreenTensor(this%specfem_kgt,jf,errmsg)
       if (.level.errmsg == 2) return
!!$    endif
  end subroutine readFrequencyKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get ID
!
  function getIdKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    character(len=length_ID_kernel_green_tensor) :: res
!
    res = ''
!!$    if (associated(this%gemini_kgt)) then
!!$       res = .id.(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       res = .id.(this%specfem_kgt)
!!$    endif
  end function getIdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    real :: res
!
    res = 0.
!!$    if (associated(this%gemini_kgt)) then
!!$       res = .df.(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       res = .df.(this%specfem_kgt)
!!$    endif
  end function getDfKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
  function getNfKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    integer :: res
!
    res = 0
!!$    if (associated(this%gemini_kgt)) then
!!$       res = .nf.(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       res = .nf.(this%specfem_kgt)
!!$    endif
  end function getNfKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency fcur
!
  function getJfcurKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    integer :: res
!
    res = 0
!!$    if (associated(this%gemini_kgt)) then
!!$       res = .jfcur.(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       res = .jfcur.(this%specfem_kgt)
!!$    endif
  end function getJfcurKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get pointer to strains
!
  function getStrainsKernelGreenTensor(this) result(gstr)
    type (kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: gstr
!!$    if (associated(this%gemini_kgt)) then
!!$       gstr => getStrainsGeminiKernelGreenTensor(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       gstr => .strain.(this%specfem_kgt)
!!$    endif
  end function getStrainsKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get pointer to displacements
!
  function getKernelGreenTensor(this) result(g)
    type (kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: g
!!$    if (associated(this%gemini_kgt)) then
!!$       g => getGeminiKernelGreenTensor(this%gemini_kgt)
!!$    else if (associated(this%specfem_kgt)) then
       g => .disp.(this%specfem_kgt)
!!$    endif
  end function getKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyKernelGreenTensor(this,jf) result(next)
    type (kernel_green_tensor) :: this
    integer :: jf
    logical :: next
!!$    if (associated(this%gemini_kgt)) then
!!$       next = nextFrequencyGeminiKernelGreenTensor(this%gemini_kgt,jf)
!!$    else if (associated(this%specfem_kgt)) then
       next = nextFrequencySpecfem3dKernelGreenTensor(this%specfem_kgt,jf)
!!$    endif
  end function nextFrequencyKernelGreenTensor
!
end module kernelGreenTensor
