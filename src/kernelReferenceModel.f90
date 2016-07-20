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
!> \brief Wrapper module for kernelReferenceModel
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  The interfaces to method GEMINI still need to be adapted to the new ASKI version,
!!  so far only method SPECFEM3D is supported. 
!
module kernelReferenceModel
  use modelParametrization
!!$  use geminiEarthModel
  use specfem3dKernelReferenceModel
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocKernelReferenceModel; end interface
  type kernel_reference_model
!!$     type (gemini_earth_model), pointer :: gemini_em => null()
     type (specfem3d_kernel_reference_model), pointer :: specfem_em => null()
  end type kernel_reference_model
!
contains
!-------------------------------------------------------------------
!> \brief  Create kernel reference model
!
  subroutine createKernelReferenceModel(this,method,fuh,filename,errmsg)
    type (kernel_reference_model) :: this
    character(len=*) :: method
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=26) :: myname = 'createKernelReferenceModel'
    !
    call addTrace(errmsg,myname)
    select case (method)
!!$    case('GEMINI')
!!$       allocate(this%gemini_em)
!!$       errmsg = readSAGeminiEarthModel(this%gemini_em,fuh,filename)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
    case('SPECFEM3D')
       allocate(this%specfem_em)
       call readSpecfem3dKernelReferenceModel(this%specfem_em,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case default
       call add(errmsg,2,'Invalid forward computation method',myname)
       return
    end select
  end subroutine createKernelReferenceModel
!------------------------------------------------------------------------
!> \brief Deallocate kernel reference model
!
  subroutine deallocKernelReferenceModel(this)
    type (kernel_reference_model) :: this
!!$    if (associated(this%gemini_em)) call dealloc(this%gemini_em)
    if (associated(this%specfem_em)) call dealloc(this%specfem_em)
  end subroutine deallocKernelReferenceModel
!------------------------------------------------------------------------
!> \brief Get model values for some particular parameter of some particular parametrization
!! \param pmtrz some model parametrization, as defined in module modelParametrization (e.g. "isoVelocity")
!! \param param model parameter, which must be one of parametrization pmtrz
!! \param model_values pointer to array of model values on wavefield points. Nullified if there is a problem
!
  function getModelValuesKernelReferenceModel(this,pmtrz,param) result(model_values)
    type (kernel_reference_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
    nullify(model_values)
    if(.not.validModelParametrization(pmtrz)) return
    if(.not.validParamModelParametrization(pmtrz,param)) return
!!$    if (associated(this%gemini_em)) then
!!$       model_values => getModelValuesWPGeminiEarthModel(this%gemini_em,pmtrz,param)
!!$    else if (associated(this%specfem_em)) then
       model_values => getModelValuesSpecfem3dKernelReferenceModel(this%specfem_em,pmtrz,param)
!!$    endif
  end function getModelValuesKernelReferenceModel
!
end module kernelReferenceModel
