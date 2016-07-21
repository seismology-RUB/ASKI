!----------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   and Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Wrapper module for kernel displacement
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  Displacement and strain arrays are assumed to be ordered consistently with ordering of wavefield points.
!
module kernelDisplacement
  use geminiKernelDisplacement
  use specfem3dKernelDisplacement
  use nexdKernelDisplacement
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateKernelDisplacement; end interface
  interface operator (.id.); module procedure getIdKernelDisplacement; end interface
  interface operator (.df.); module procedure getDfKernelDisplacement; end interface
  interface operator (.nf.); module procedure getNfKernelDisplacement; end interface
  interface operator (.jfcur.); module procedure getJfcurKernelDisplacement; end interface
!< fork here to specific method by associating some pointer (only ONE!)
  type kernel_displacement
     private
     type (gemini_kernel_displacement), pointer :: gemini_kd => null()         ! gemini specific kernel displacement object
     type (specfem3d_kernel_displacement), pointer :: specfem3d_kd => null()       ! Specfem specific kernel displacement object
     type (nexd_kernel_displacement), pointer :: nexd_kd => null()       ! NEXD specific kernel displacement object
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
    case('GEMINI')
       allocate(this%gemini_kd)
       call initialReadGeminiKernelDisplacement(this%gemini_kd,fuh,filename,errmsg)
       if (.level.errmsg == 2) then
          call addTrace(errmsg,myname); return
       endif
    case('SPECFEM3D')
       allocate(this%specfem3d_kd)
       call initialReadSpecfem3dKernelDisplacement(this%specfem3d_kd,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case('NEXD')
       allocate(this%nexd_kd)
       call initialReadNexdKernelDisplacement(this%nexd_kd,fuh,filename,errmsg)
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
    if (associated(this%gemini_kd)) call deallocGeminiKernelDisplacement(this%gemini_kd,fuh)
    if (associated(this%specfem3d_kd)) call dealloc(this%specfem3d_kd,fuh)
    if (associated(this%nexd_kd)) call dealloc(this%nexd_kd,fuh)
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
    if (associated(this%gemini_kd)) then
       call readFrequencyGeminiKernelDisplacement(this%gemini_kd,jf,errmsg)
       if (.level.errmsg == 2) then
          call addTrace(errmsg,myname); return
       endif
    else if (associated(this%specfem3d_kd)) then
       call readFrequencySpecfem3dKernelDisplacement(this%specfem3d_kd,jf,errmsg)
       if (.level.errmsg == 2) return
    else if (associated(this%nexd_kd)) then
       call readFrequencyNexdKernelDisplacement(this%nexd_kd,jf,errmsg)
       if (.level.errmsg == 2) return
    endif
  end subroutine readFrequencyKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get ID
!
  function getIdKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    character(len=length_ID_kernel_displacement) :: res
!
    res = ''
    if (associated(this%gemini_kd)) then
       res = .id.(this%gemini_kd)
    else if (associated(this%specfem3d_kd)) then
       res = .id.(this%specfem3d_kd)
    else if (associated(this%nexd_kd)) then
       res = .id.(this%nexd_kd)
    endif
  end function getIdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    real :: res
!
    res = 0.
    if (associated(this%gemini_kd)) then
       res = .df.(this%gemini_kd)
    else if (associated(this%specfem3d_kd)) then
       res = .df.(this%specfem3d_kd)
    else if (associated(this%nexd_kd)) then
       res = .df.(this%nexd_kd)
    endif
  end function getDfKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
  function getNfKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    integer :: res
!
    res = 0
    if (associated(this%gemini_kd)) then
       res = .nf.(this%gemini_kd)
    else if (associated(this%specfem3d_kd)) then
       res = .nf.(this%specfem3d_kd)
    else if (associated(this%nexd_kd)) then
       res = .nf.(this%nexd_kd)
    endif
  end function getNfKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurKernelDisplacement(this) result(res)
    type (kernel_displacement), intent(in) :: this
    integer :: res
!
    res = 0
    if (associated(this%gemini_kd)) then
       res = .jfcur.(this%gemini_kd)
    else if (associated(this%specfem3d_kd)) then
       res = .jfcur.(this%specfem3d_kd)
    else if (associated(this%nexd_kd)) then
       res = .jfcur.(this%nexd_kd)
    endif
  end function getJfcurKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get pointer to strains
!! \details Allocate here for the return array ustr (allocation can also be done in the
!!  method-specific realization of kernelDisplacement). This is done in order to be most flexible,
!!  not enforcing a method to store a pointer of for ustr inside its kernelDisplacement object. 
!!  Routines calling getStrainsKernelDisplacement must care for deallocation of pointer ustr!
!
  subroutine getStrainsKernelDisplacement(this,ustr,errmsg)
    type (kernel_displacement) :: this
    complex, dimension(:,:), pointer :: ustr
    type (error_message) :: errmsg
    character (len=28) :: myname = 'getStrainsKernelDisplacement'
    ! local
    complex, dimension(:,:), pointer :: ustr_local
!
    if (associated(this%gemini_kd)) then
       ustr_local => getStrainsGeminiKernelDisplacement(this%gemini_kd)
       if(.not.associated(ustr_local)) then
          call add(errmsg,2,"no kernel displacement strains were returned by function getStrainsGeminiKernelDisplacement",myname)
          return
       end if
    else if (associated(this%specfem3d_kd)) then
       ustr_local => getStrainsSpecfem3dKernelDisplacement(this%specfem3d_kd)
       if(.not.associated(ustr_local)) then
          call add(errmsg,2,"no kernel displacement strains were returned by function getStrainsSpecfem3dKernelDisplacement",myname)
          return
       end if
    else if (associated(this%nexd_kd)) then
       ustr_local => getStrainsNexdKernelDisplacement(this%nexd_kd)
       if(.not.associated(ustr_local)) then
          call add(errmsg,2,"no kernel displacement strains were returned by function getStrainsNexdKernelDisplacement",myname)
          return
       end if
    endif
    allocate(ustr(size(ustr_local,1),size(ustr_local,2)))
    ustr = ustr_local
  end subroutine getStrainsKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsKernelDisplacement(this,uf_ustr,errmsg)
    type (kernel_displacement) :: this
    real :: uf_ustr
    type (error_message) :: errmsg
    if (associated(this%gemini_kd)) then
       call getUnitFactorStrainsGeminiKernelDisplacement(this%gemini_kd,uf_ustr)
    else if (associated(this%specfem3d_kd)) then
       call getUnitFactorStrainsSpecfem3dKernelDisplacement(this%specfem3d_kd,uf_ustr,errmsg)
    else if (associated(this%nexd_kd)) then
       call getUnitFactorStrainsNexdKernelDisplacement(this%nexd_kd,uf_ustr)
    end if
  end subroutine getUnitFactorStrainsKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get pointer to displacements
!! \details Allocate here for the return array u (allocation can also be done in the
!!  method-specific realization of kernelDisplacement). This is done in order to be most flexible,
!!  not enforcing a method to store a pointer of for u inside its kernelDisplacement object. 
!!  Routines calling getKernelDisplacement must care for deallocation of pointer u!
!
  subroutine getKernelDisplacement(this,u,errmsg)
    type (kernel_displacement) :: this
    complex, dimension(:,:), pointer :: u
    type (error_message) :: errmsg
    character (len=21) :: myname = 'getKernelDisplacement'
    ! local
    complex, dimension(:,:), pointer :: u_local
!
    if (associated(this%gemini_kd)) then
       u_local => getGeminiKernelDisplacement(this%gemini_kd)
       if(.not.associated(u_local)) then
          call add(errmsg,2,"no kernel displacement values were returned by function getGeminiKernelDisplacement",myname)
          return
       end if
    else if (associated(this%specfem3d_kd)) then
       u_local => getSpecfem3dKernelDisplacement(this%specfem3d_kd)
       if(.not.associated(u_local)) then
          call add(errmsg,2,"no kernel displacement values were returned by function getSpecfem3dKernelDisplacement",myname)
          return
       end if
    else if (associated(this%nexd_kd)) then
       u_local => getNexdKernelDisplacement(this%nexd_kd)
       if(.not.associated(u_local)) then
          call add(errmsg,2,"no kernel displacement values were returned by function getNexdKernelDisplacement",myname)
          return
       end if
    endif
    allocate(u(size(u_local,1),size(u_local,2)))
    u = u_local
  end subroutine getKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of displacements
!
  subroutine getUnitFactorKernelDisplacement(this,uf_u,errmsg)
    type (kernel_displacement) :: this
    real :: uf_u
    type (error_message) :: errmsg
    if (associated(this%gemini_kd)) then
       call getUnitFactorGeminiKernelDisplacement(this%gemini_kd,uf_u)
    else if (associated(this%specfem3d_kd)) then
       call getUnitFactorSpecfem3dKernelDisplacement(this%specfem3d_kd,uf_u,errmsg)
    else if (associated(this%nexd_kd)) then
       call getUnitFactorNexdKernelDisplacement(this%nexd_kd,uf_u)
    end if
  end subroutine getUnitFactorKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyKernelDisplacement(this,jf) result(next)
    type (kernel_displacement) :: this
    integer :: jf
    logical :: next
    if (associated(this%gemini_kd)) then
       next = nextFrequencyGeminiKernelDisplacement(this%gemini_kd,jf)
    else if (associated(this%specfem3d_kd)) then
       next = nextFrequencySpecfem3dKernelDisplacement(this%specfem3d_kd,jf)
    else if (associated(this%nexd_kd)) then
       next = nextFrequencyNexdKernelDisplacement(this%nexd_kd,jf)
    endif
  end function nextFrequencyKernelDisplacement
!
end module kernelDisplacement
