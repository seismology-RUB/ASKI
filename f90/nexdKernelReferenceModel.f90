!----------------------------------------------------------------------------
!   Copyright 2016 Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module dealing with kernel reference model for method NEXD
!! \details The code of this module is based on module specfem3dKernelReferenceModel
!!
!! \author Christian Ullisch
!! \date Nov 2015
!
module nexdKernelReferenceModel
  use modelParametrization
  use errorMessage
  implicit none
  interface dealloc; module procedure deallocNexdKernelReferenceModel; end interface
  type nexd_kernel_reference_model
     private
     real, dimension(:), pointer :: rho => null()                   !< rho
     real, dimension(:), pointer :: vp => null()                    !< vp
     real, dimension(:), pointer :: vs => null()                    !< vs
  end type nexd_kernel_reference_model
!
contains
!--------------------------------------------------------------------------------------------
!> \brief Deallocate nexd_kernel_reference_model
!
  subroutine deallocNexdKernelReferenceModel(this)
    type (nexd_kernel_reference_model) :: this
    if (associated(this%rho)) deallocate(this%rho)
    if (associated(this%vp)) deallocate(this%vp)
    if (associated(this%vs)) deallocate(this%vs)
  end subroutine deallocNexdKernelReferenceModel
!------------------------------------------------------------------------------
!> \brief Get earth model values at wavefield points
!
  function getModelValuesNexdKernelReferenceModel(this,pmtrz,param) result(model_values)
    type (nexd_kernel_reference_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
    nullify(model_values)
    if(.not.associated(this%rho)) return ! model was not yet read in
    if(.not.validModelParametrization(pmtrz)) return
    if(.not.validParamModelParametrization(pmtrz,param)) return
    select case(pmtrz)
       case('isoVelocitySI')
          allocate(model_values(size(this%rho)))
          select case(param)
             case('rho'); model_values = this%rho
             case('vp'); model_values = this%vp
             case('vs'); model_values = this%vs
          end select
       case('isoVelocity1000')
          allocate(model_values(size(this%rho)))
          select case(param)
             case('rho'); model_values = this%rho
             case('vp'); model_values = this%vp
             case('vs'); model_values = this%vs
          end select
          model_values = model_values * 1.0e-3 ! scale from SI units to m/s or g/cm^3
       case('isoLameSI')
          ! could compute mu, lambda from rho,vp,vs and return values
    end select
  end function getModelValuesNexdKernelReferenceModel
!----------------------------------------------------------------------------
!> \brief  Read in specfem earth model
!
  subroutine readNexdKernelReferenceModel(this,lu,filename,errmsg)
    type (nexd_kernel_reference_model) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    integer :: ier,ntot_wp,nf,nproc
    integer, dimension(:), allocatable :: jf
    real :: df
    real, dimension(:), allocatable :: x,y,z
    character (len=28) :: myname = 'readNexdKernelReferenceModel'
!
    call addTrace(errmsg,myname)
!
    open(unit=lu,file=filename,status='old',action='read',access='stream',form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
       goto 1
    endif

    read(lu) nproc,ntot_wp,df,nf

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
    allocate(this%rho(ntot_wp),this%vp(ntot_wp),this%vs(ntot_wp))
    read(lu) this%rho
    read(lu) this%vp
    read(lu) this%vs

1   close(lu)
    if(allocated(jf)) deallocate(jf)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)
  end subroutine readNexdKernelReferenceModel
!
end module nexdKernelReferenceModel
