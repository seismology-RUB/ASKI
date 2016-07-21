!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief Module dealing with kernel reference model for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module specfem3dKernelReferenceModel
  use specfem3dForASKIFiles
  use modelParametrization
  use errorMessage
  implicit none
  interface dealloc; module procedure deallocSpecfem3dKernelReferenceModel; end interface
  type specfem3d_kernel_reference_model
     private
     integer :: specfem_version = -1 !< 1 = SPECFEM3D_GLOBE, 2 = SPEECFEM3D_Cartesian
     real, dimension(:), pointer :: rho => null()                    !< rho, only isotropic model supported at the moment!
     real, dimension(:), pointer :: vp => null()                    !< vp, only isotropic model supported at the moment!
     real, dimension(:), pointer :: vs => null()                    !< vs, only isotropic model supported at the moment!
  end type specfem3d_kernel_reference_model
!
contains
!----------------------------------------------------------------------------
!> \brief  Read in specfem earth model
!
  subroutine readSpecfem3dKernelReferenceModel(this,lu,filename,errmsg)
    type (specfem3d_kernel_reference_model) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    character (len=33) :: myname = 'readSpecfem3dKernelReferenceModel'
!
    call addTrace(errmsg,myname)
!
    call readSpecfem3dForASKIMainFile(filename,lu,errmsg,specfem_version=this%specfem_version,&
         rho=this%rho,vp=this%vp,vs=this%vs)
    if(.level.errmsg == 2) goto 2
!
    select case(this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files,version_cartesian_specfem3d_for_ASKI_files) ! OK, do nothing
    case default
       write(errstr,*) "readSpecfem3dForASKIMainFile returned specfem_version = ",this%specfem_version,&
            "; this module only supports specfem versions ",version_globe_specfem3d_for_ASKI_files," = SPECFEM3D_GLOBE and ",&
            version_cartesian_specfem3d_for_ASKI_files," = SPECFEM3D_Cartesian ."
       call add(errmsg,2,errstr,myname)
       goto 2
    end select
!
    if(.not.(associated(this%rho).and.associated(this%vp).and.associated(this%vs))) then
       call add(errmsg,2,"readSpecfem3dForASKIMainFile did not return all of rho,vp,vs vectors",myname)
       goto 2
    end if
!
    if(.not.(size(this%rho)==size(this%vp) .and. size(this%rho)==size(this%vs))) then
       call add(errmsg,2,"readSpecfem3dForASKIMainFile returned vectors rho,vp,vs which are not all of same size",myname)
       goto 2
    end if
!
    ! if code comes here, return normally.
    return
!
2   call deallocSpecfem3dKernelReferenceModel(this)
    return
  end subroutine readSpecfem3dKernelReferenceModel
!--------------------------------------------------------------------------------------------
!> \brief Deallocate specfem3d_kernel_reference_model
!
  subroutine deallocSpecfem3dKernelReferenceModel(this)
    type (specfem3d_kernel_reference_model) :: this
    this%specfem_version = -1
    if (associated(this%rho)) deallocate(this%rho)
    if (associated(this%vp)) deallocate(this%vp)
    if (associated(this%vs)) deallocate(this%vs)
  end subroutine deallocSpecfem3dKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get earth model values on wavefield points
!
  function getModelValuesSpecfem3dKernelReferenceModel(this,pmtrz,param) result(model_values)
    type (specfem3d_kernel_reference_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
!
    nullify(model_values)
    if(this%specfem_version <= 0) return ! this%specfem_version == -1 indicates that the model was not yet read from file
    if(.not.validModelParametrization(pmtrz)) return
    if(.not.validParamModelParametrization(pmtrz,param)) return
!
    select case(pmtrz)
       case('isoVelocitySI')
          allocate(model_values(size(this%rho)))
          select case(param)
             case('rho'); model_values = this%rho
             case('vp'); model_values = this%vp
             case('vs'); model_values = this%vs
          end select
          ! in case of SPECFEM3D_GLOBE, model_values are in g/cm^3 or km/s , so need to multiply here by 1000.0 to yield SI units
          if(this%specfem_version == version_globe_specfem3d_for_ASKI_files) model_values = model_values * 1.0e3
!
       case('isoVelocity1000')
          allocate(model_values(size(this%rho)))
          select case(param)
             case('rho'); model_values = this%rho
             case('vp'); model_values = this%vp
             case('vs'); model_values = this%vs
          end select
          ! in case of SPECFEM3D_Cartesian, model_values are in Kg/m^3 or m/s , so need to divide here by 1000.0 to yield isoVelocity1000 values
          if(this%specfem_version == version_cartesian_specfem3d_for_ASKI_files) model_values = model_values * 1.0e-3
       case('isoLameSI')
          ! could compute mu, lambda from rho,vp,vs and return values
    end select
  end function getModelValuesSpecfem3dKernelReferenceModel
!
end module specfem3dKernelReferenceModel
