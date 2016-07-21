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
!> \brief Module dealing with wavefield points for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date June 2013
!
module specfem3dWavefieldPoints
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateSpecfem3dWavefieldPoints; end interface
  interface operator (.csys.); module procedure getCoordinateSystemSpecfem3dWavefieldPoints; end interface
  interface operator (.R.); module procedure getRadiusSpecfem3dWavefieldPoints; end interface
  type specfem3d_wavefield_points
     private
     integer :: specfem_version = 0 !< 1 = SPECFEM3D_GLOBE, 2 = SPEECFEM3D_Cartesian
     integer :: type_inversion_grid = 0 !< type of inversion grid
     character (len=1) :: csys                         ! coordinate system (spherical, cartesian)
     real :: R                                         ! Radius of sphere in m (in case of csys == 'S')
     real :: clon,ctheta                 ! center of kernel chunk in degrees ('S') or m ('C')
		                                    ! ('S'): clon=longitude, ctheta=colatitude; ('C'): clon=y-coordinate, ctheta=x-coordinate of center
     real :: cdepth                                    ! center of depth in m ('S'), center of z in m ('C')
     real :: wlon,wlat                                 ! extension of kernel area in degrees ('S') or m ('C')
		                                                  ! ('C'): wlon=y-extension, wlat=x-extension
     real :: wdepth                                    ! extension of depth in m('S'), extension of z in m ('C')
     real :: gam                                       ! rotation on kernel chunk anticlockwise
     real, dimension(:), pointer :: x => null()        ! x coordinates of points
     real, dimension(:), pointer :: y => null()        ! y coordinates of points
     real, dimension(:), pointer :: z => null()        ! z coordinates of points
  end type specfem3d_wavefield_points
!
  integer, parameter :: length_ID_specfem3d_wavefield_points = 13 !< change this number CONSISTENTLY with length_ASKI_output_ID in SPECFEM codes
!
contains
!------------------------------------------------------------------
!> \brief Read in wavefield point information from file
!
  subroutine readSpecfem3dWavefieldPoints(this,fuh,filename,errmsg)
    type (specfem3d_wavefield_points) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    character (len = 28) :: myname = 'readSpecfem3dWavefieldPoints'
    integer :: ntot_wp,ier,lu,nf,nproc,len_id
    character(len=length_ID_specfem3d_wavefield_points) :: id
    integer, dimension(:), allocatable :: jf
    double precision :: df
    !
    call addTrace(errmsg,myname)
    !
    lu = get(fuh)
    open(unit=lu,file=filename,status='old',action='read',access='stream',form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
       goto 1
    endif

! THIS IS FOR COMPATABILITY WITH OLD INVERSION GRID, ONLY FOR EXAMPLE test/new_specfem3d
    this%R = 0.; this%csys = 'C'
    this%ctheta=0.;this%clon=0.;this%cdepth=155.0;this%wlat=128.0;this%wlon=128.0;this%wdepth=128.0

    read(lu) this%specfem_version
    select case(this%specfem_version)
    case(1,2)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "specfem version = ",this%specfem_version," is not supported: 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select

    read(lu) len_id
    if(len_id/=length_ID_specfem3d_wavefield_points) then
       write(errstr,*) "length of ASKI ID is ",len_id,&
            ", only length ",length_ID_specfem3d_wavefield_points," supported here. Please update this code accordingly"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if

    read(lu) id,nproc,this%type_inversion_grid,ntot_wp,df,nf

    select case(this%type_inversion_grid)
    case(2,3,4)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "type of inversion grid = ",this%type_inversion_grid," is not supported by SPECFEM3D"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select

    ! check if this%specfem_version and this%type_inversion_grid are a valid match
    if(this%specfem_version == 1) then
       select case(this%type_inversion_grid)
       case(1,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",this%type_inversion_grid," is not supported by SPECFEM3D_GLOBE"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(this%specfem_version == 2) then
       select case(this%type_inversion_grid)
       case(2,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",this%type_inversion_grid," is not supported by SPECFEM3D_Cartesian"
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
    allocate(this%x(ntot_wp),this%y(ntot_wp),this%z(ntot_wp))
    read(lu) this%x
    read(lu) this%y
    read(lu) this%z

1   close(lu)
    call add(fuh,lu)
    if(allocated(jf)) deallocate(jf)
  end subroutine readSpecfem3dWavefieldPoints
!--------------------------------------------------------------------------
!> \brief Deallocate specfem3d_wavefield_points
!
  subroutine deallocateSpecfem3dWavefieldPoints(this)
    type (specfem3d_wavefield_points ) :: this
    this%specfem_version = 0
    this%type_inversion_grid = 0
    if (associated(this%x)) deallocate(this%x)
    if (associated(this%y)) deallocate(this%y)
    if (associated(this%z)) deallocate(this%z)
  end subroutine deallocateSpecfem3dWavefieldPoints
!--------------------------------------------------------------------------------
!> \brief Get coordinate system
!
  function getCoordinateSystemSpecfem3dWavefieldPoints(this) result(csys)
    type (specfem3d_wavefield_points), intent(in) :: this
    character (len=1) :: csys
    csys = this%csys
  end function getCoordinateSystemSpecfem3dWavefieldPoints
!--------------------------------------------------------------------------------
!> \brief Get radius
!
  function getRadiusSpecfem3dWavefieldPoints(this) result(R)
    type (specfem3d_wavefield_points), intent(in) :: this
    real :: R
    R = this%R
  end function getRadiusSpecfem3dWavefieldPoints
!-----------------------------------------------------------------------------
!> \brief Get total number of wavefield points
!
  function getNtotSpecfem3dWavefieldPoints(this) result(ntot)
    type (specfem3d_wavefield_points), intent(in) :: this
    integer :: ntot
    if(associated(this%x)) then
       ntot = size(this%x)
    else
       ntot = 0
    endif
  end function getNtotSpecfem3dWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the specfem wavefield poinst
! get arrays of wavefield points as are, inversion grid must deal with them
  subroutine getSpecfem3dWavefieldPoints(this,x,y,z)
    type (specfem3d_wavefield_points), intent(in) :: this
    real, dimension(:), pointer :: x,y,z
!
    nullify(x,y,z)
    if(associated(this%x)) then
       allocate(x(size(this%x)),y(size(this%x)),z(size(this%x)))
       x = this%x
       y = this%y
       z = this%z
    end if
  end subroutine getSpecfem3dWavefieldPoints
!
end module specfem3dWavefieldPoints
