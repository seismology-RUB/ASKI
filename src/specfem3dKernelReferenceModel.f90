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
!> \brief Module dealing with kernel reference model for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date June 2013
!
module specfem3dKernelReferenceModel
  use modelParametrization
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocSpecfem3dKernelReferenceModel; end interface
  type specfem3d_kernel_reference_model
     private
     character (len=80) :: name                                      !< name of earth model
     real, dimension(:), pointer :: rho => null()                    !< rho
     real, dimension(:), pointer :: vph => null()                    !< vph
     real, dimension(:), pointer :: vpv => null()                    !< vpv
     real, dimension(:), pointer :: vsh => null()                    !< vsh
     real, dimension(:), pointer :: vsv => null()                    !< vsv
  end type specfem3d_kernel_reference_model
!
  integer, parameter :: length_ID_specfem3d_kernel_reference_model = 13 !< change this number CONSISTENTLY with length_ASKI_output_ID in SPECFEM codes
!
contains
!--------------------------------------------------------------------------------------------
!> \brief Deallocate specfem3d_kernel_reference_model
!
  subroutine deallocSpecfem3dKernelReferenceModel(this)
    type (specfem3d_kernel_reference_model) :: this
    if (associated(this%rho)) deallocate(this%rho)
    if (associated(this%vph)) deallocate(this%vph)
    if (associated(this%vpv)) deallocate(this%vpv)
    if (associated(this%vsh)) deallocate(this%vsh)
    if (associated(this%vsv)) deallocate(this%vsv)
  end subroutine deallocSpecfem3dKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get density of earth model at wavefield points
!
  function getDensitySpecfem3dKernelReferenceModel(this) result(rho)
    type (specfem3d_kernel_reference_model) :: this
    real, dimension(:), pointer :: rho
    !
    rho => this%rho
  end function getDensitySpecfem3dKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get P velocity of earth model at wavefield points
!
  function getPVelocitySpecfem3dKernelReferenceModel(this) result(vp)
    type (specfem3d_kernel_reference_model) :: this
    real, dimension(:), pointer :: vp
    !
    vp => this%vph
  end function getPVelocitySpecfem3dKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get S velocity of earth model at wavefield points
!
  function getSVelocitySpecfem3dKernelReferenceModel(this) result(vs)
    type (specfem3d_kernel_reference_model) :: this
    real, dimension(:), pointer :: vs
    !
    vs => this%vsh
  end function getSVelocitySpecfem3dKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get S velocity of earth model at wavefield points
!
  function getModelValuesSpecfem3dKernelReferenceModel(this,pmtrz,param) result(model_values)
    type (specfem3d_kernel_reference_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
    nullify(model_values)
    if(.not.validModelParametrization(pmtrz)) return
    if(.not.validParamModelParametrization(pmtrz,param)) return
    select case(pmtrz)
       case('isoVelocity')
          select case(param)
             case('rho'); model_values => this%rho
             case('vp'); model_values => this%vph
             case('vs'); model_values => this%vsh
          end select
       case('isoLame')
          ! could compute mu, lambda from rho,vp,vs and return values
    end select
  end function getModelValuesSpecfem3dKernelReferenceModel
!----------------------------------------------------------------------------
!> \brief  Read in specfem earth model
!
  subroutine readSpecfem3dKernelReferenceModel(this,fuh,filename,errmsg)
    type (specfem3d_kernel_reference_model) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    integer :: lu,ier,ntot_wp,nf,specfem_version,type_invgrid,nproc,len_id
    character(len=length_ID_specfem3d_kernel_reference_model) :: id
    integer, dimension(:), allocatable :: jf,nb
    double precision :: df
    integer :: NGLLX,NGLLY,NGLLZ,ncell,nnb_total
    real, dimension(:), allocatable :: x,y,z,jacobian
    character (len=33) :: myname = 'readSpecfem3dKernelReferenceModel'
!
    call addTrace(errmsg,myname)
!
    lu = get(fuh)
    open(unit=lu,file=filename,status='old',action='read',access='stream',form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
       goto 1
    endif

    read(lu) specfem_version
    select case(specfem_version)
    case(1,2)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "specfem version = ",specfem_version," is not supported: 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select

    read(lu) len_id
    if(len_id/=length_ID_specfem3d_kernel_reference_model) then
       write(errstr,*) "length of ASKI ID is ",len_id,&
            ", only length ",length_ID_specfem3d_kernel_reference_model," supported here. Please update this code accordingly"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if

    read(lu) id,nproc,type_invgrid,ntot_wp,df,nf

    select case(type_invgrid)
    case(2,3,4)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "type of inversion grid = ",type_invgrid," is not supported by SPECFEM3D"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select

    ! check if specfem_version and type_invgrid are a valid match
    if(specfem_version == 1) then
       select case(type_invgrid)
       case(1,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_invgrid," is not supported by SPECFEM3D_GLOBE"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(specfem_version == 2) then
       select case(type_invgrid)
       case(2,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_invgrid," is not supported by SPECFEM3D_Cartesian"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(nf<1 .or. ntot_wp<1) then
       write(errstr,*) "ntot_wp,df,nf =",ntot_wp,df,nf,"; ntot_wp and nf should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    allocate(jf(nf))
    read(lu) jf
    allocate(x(ntot_wp),y(ntot_wp),z(ntot_wp))
    read(lu) x
    read(lu) y
    read(lu) z
    select case(type_invgrid)
    case(4)
       read(lu) NGLLX,NGLLY,NGLLZ
       allocate(jacobian(ntot_wp))
       read(lu) jacobian
       read(lu) ncell,nnb_total
       if(ncell /= ntot_wp/(NGLLX*NGLLY*NGLLZ) .or. ncell<1) then
          write(errstr,*) "ncell = ",ncell," does not compute as ntot_wp/(NGLLX*NGLLY*NGLLZ) = ",&
               ntot_wp/(NGLLX*NGLLY*NGLLZ)," or it is negative -> file is inconsistent!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(nnb_total<ncell) then
          write(errstr,*) "total number of integers defining the neighbours ",nnb_total,&
               " should be at least the number of cells ",ncell," -> file is inconsistent!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       allocate(nb(nnb_total))
       read(lu) nb
    end select
    allocate(this%rho(ntot_wp),this%vph(ntot_wp),this%vsh(ntot_wp))
    read(lu) this%rho
    read(lu) this%vph
    read(lu) this%vsh

1   close(lu)
    call add(fuh,lu)
    if(allocated(jf)) deallocate(jf)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)
    if(allocated(jacobian)) deallocate(jacobian)
    if(allocated(nb)) deallocate(nb)
  end subroutine readSpecfem3dKernelReferenceModel
!
end module specfem3dKernelReferenceModel
