!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   and Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!----------------------------------------------------------------------------
!> \brief Wrapper module for kernelReferenceModel
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  The interfaces to method GEMINI still need to be adapted to the new ASKI version,
!!  so far only method SPECFEM3D is supported. 
!
module kernelReferenceModel
  use modelParametrization
  use geminiKernelReferenceModel
  use specfem3dKernelReferenceModel
  use nexdKernelReferenceModel
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateKernelReferenceModel; end interface
  type kernel_reference_model
     type (gemini_kernel_reference_model), pointer :: gemini_krm => null()
     type (specfem3d_kernel_reference_model), pointer :: specfem3d_krm => null()
     type (nexd_kernel_reference_model), pointer :: nexd_krm => null()
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
    case('GEMINI')
       allocate(this%gemini_krm)
       call readGeminiKernelReferenceModel(this%gemini_krm,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case('SPECFEM3D')
       allocate(this%specfem3d_krm)
       call readSpecfem3dKernelReferenceModel(this%specfem3d_krm,get(fuh),filename,errmsg)
       call undo(fuh)
       if (.level.errmsg == 2) return
    case('NEXD')
       allocate(this%nexd_krm)
       call readNexdKernelReferenceModel(this%nexd_krm,get(fuh),filename,errmsg)
       call undo(fuh)
       if (.level.errmsg == 2) return
    case default
       call add(errmsg,2,'Invalid forward computation method',myname)
       return
    end select
  end subroutine createKernelReferenceModel
!------------------------------------------------------------------------
!> \brief Deallocate kernel reference model
!
  subroutine deallocateKernelReferenceModel(this)
    type (kernel_reference_model) :: this
    if (associated(this%gemini_krm)) call dealloc(this%gemini_krm)
    if (associated(this%specfem3d_krm)) call dealloc(this%specfem3d_krm)
    if (associated(this%nexd_krm)) call dealloc(this%nexd_krm)
  end subroutine deallocateKernelReferenceModel
!------------------------------------------------------------------------
!> \brief Get model values for some particular parameter of some particular parametrization
!! \details When writing a method-specific version of this routine, allocate space for
!!  model_values. Some realizations of kernelReferenceModel may not store model values on wavefield points
!!  in the object but construct them on the fly here. Routines calling getModelValuesKernelReferenceModel 
!!  should care for deallocation of pointer model_values!
!! \param pmtrz some model parametrization, as defined in module modelParametrization (e.g. "isoVelocitySI")
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
    if (associated(this%gemini_krm)) then
       model_values => getModelValuesWPGeminiKernelReferenceModel(this%gemini_krm,pmtrz,param)
    else if (associated(this%specfem3d_krm)) then
       model_values => getModelValuesSpecfem3dKernelReferenceModel(this%specfem3d_krm,pmtrz,param)
    else if (associated(this%nexd_krm)) then
       model_values => getModelValuesNexdKernelReferenceModel(this%nexd_krm,pmtrz,param)
    endif
  end function getModelValuesKernelReferenceModel
!
end module kernelReferenceModel
