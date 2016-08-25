!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module dealing with wavefield points for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module specfem3dWavefieldPoints
  use specfem3dForASKIFiles
  use errorMessage
  implicit none
  interface dealloc; module procedure deallocateSpecfem3dWavefieldPoints; end interface
  type specfem3d_wavefield_points
     private
     integer :: specfem_version = -1 !< 1 = SPECFEM3D_GLOBE, 2 = SPEECFEM3D_Cartesian
     integer :: nwp = 0 !< number of wavefield points
     real, dimension(:), pointer :: x => null()        ! x coordinates of points
     real, dimension(:), pointer :: y => null()        ! y coordinates of points
     real, dimension(:), pointer :: z => null()        ! z coordinates of points
  end type specfem3d_wavefield_points
!
contains
!------------------------------------------------------------------
!> \brief Read in wavefield point information from file
!
  subroutine readSpecfem3dWavefieldPoints(this,lu,filename,errmsg)
    type (specfem3d_wavefield_points) :: this
    integer :: lu
    character(len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    character(len=28) :: myname = 'readSpecfem3dWavefieldPoints'
!
    call addTrace(errmsg,myname)
!
    call readSpecfem3dForASKIMainFile(filename,lu,errmsg,specfem_version=this%specfem_version,&
         nwp=this%nwp,x=this%x,y=this%y,z=this%z)
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
    if(this%nwp <= 0) then
       write(errstr,*) "readSpecfem3dForASKIMainFile returned number of wavefield points = ",this%nwp,". Must be positive"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    if(.not.(associated(this%x).and.associated(this%y).and.associated(this%z))) then
       call add(errmsg,2,"readSpecfem3dForASKIMainFile did not return all of x,y,z vectors",myname)
       goto 2
    end if
!
    if(size(this%x)/=this%nwp .or. size(this%y)/=this%nwp .or. size(this%z)/=this%nwp) then
       write(errstr,*) "readSpecfem3dForASKIMainFile returned vectors x,y,z which do not all have expected size ",this%nwp
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    ! if code comes here, return normally.
    return
!
    ! if there was an error above, destruct this object before returning
2   call deallocateSpecfem3dWavefieldPoints(this)
    return
  end subroutine readSpecfem3dWavefieldPoints
!--------------------------------------------------------------------------
!> \brief Deallocate specfem3d_wavefield_points
!
  subroutine deallocateSpecfem3dWavefieldPoints(this)
    type (specfem3d_wavefield_points ) :: this
    this%specfem_version = -1
    this%nwp = 0
    if (associated(this%x)) deallocate(this%x)
    if (associated(this%y)) deallocate(this%y)
    if (associated(this%z)) deallocate(this%z)
  end subroutine deallocateSpecfem3dWavefieldPoints
!-----------------------------------------------------------------------------
!> \brief Get total number of wavefield points
!
  function getNtotSpecfem3dWavefieldPoints(this) result(ntot)
    type (specfem3d_wavefield_points), intent(in) :: this
    integer :: ntot
    ntot = this%nwp
  end function getNtotSpecfem3dWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the specfem wavefield points
! get arrays of wavefield points as are, inversion grid must deal with them
  subroutine getSpecfem3dWavefieldPoints(this,x,y,z)
    type (specfem3d_wavefield_points), intent(in) :: this
    real, dimension(:), pointer :: x,y,z
!
    nullify(x,y,z)
    if(associated(this%x)) then
       allocate(x(this%nwp),y(this%nwp),z(this%nwp))
       x = this%x
       y = this%y
       z = this%z
    end if
  end subroutine getSpecfem3dWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the unit factor of specfem wavefield point coordinates
  subroutine getUnitFactorSpecfem3dWavefieldPoints(this,uf_wp,errmsg)
    type (specfem3d_wavefield_points), intent(in) :: this
    real :: uf_wp
    type (error_message) :: errmsg
    select case(this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files)
       ! SPECFEM3D_GLOBE: wavefield points in the *.main output file are in units of km , i.e. unit factor is 1000
       uf_wp = 1.0e3
    case(version_cartesian_specfem3d_for_ASKI_files)
       ! SPEECFEM3D_Cartesia: ASKI for SPECFEM3D_Cartesian assumes SI units
       uf_wp = 1.0
    case default
       ! wavefield points were not yet read from file, do not know the unit factor
       call add(errmsg,2,"wavefield points were not yet read from file, the unit factor is unknown",&
            "getUnitFactorSpecfem3dWavefieldPoints")
       uf_wp = 0
       return
    end select
  end subroutine getUnitFactorSpecfem3dWavefieldPoints
!
end module specfem3dWavefieldPoints
