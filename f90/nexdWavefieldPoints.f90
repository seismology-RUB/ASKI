!----------------------------------------------------------------------------
!   Copyright 2016 Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module dealing with wavefield points for method NEXD
!! \details The code of this module is based on module specfem3dWavefieldPoints
!!
!! \author Christian Ullisch
!! \date November 2015
!
module nexdWavefieldPoints
  use errorMessage
  implicit none
  interface dealloc; module procedure deallocateNexdWavefieldPoints; end interface
  type nexd_wavefield_points
     private
     real, dimension(:), pointer :: x => null()        ! x coordinates of points
     real, dimension(:), pointer :: y => null()        ! y coordinates of points
     real, dimension(:), pointer :: z => null()        ! z coordinates of points
  end type nexd_wavefield_points
!
contains
!------------------------------------------------------------------
!> \brief Read in wavefield point information from file
!
  subroutine readNexdWavefieldPoints(this,lu,filename,errmsg)
    type (nexd_wavefield_points) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    character (len = 23) :: myname = 'readNexdWavefieldPoints'
    integer :: ntot_wp,ier,nf,nproc
    integer, dimension(:), allocatable :: jf
    real :: df

    call addTrace(errmsg,myname)

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
    allocate(this%x(ntot_wp),this%y(ntot_wp),this%z(ntot_wp))
    read(lu) this%x
    read(lu) this%y
    read(lu) this%z

1   close(lu)
    if(allocated(jf)) deallocate(jf)
  end subroutine readNexdWavefieldPoints
!--------------------------------------------------------------------------
!> \brief Deallocate nexd_wavefield_points
!
  subroutine deallocateNexdWavefieldPoints(this)
    type (nexd_wavefield_points ) :: this
    if (associated(this%x)) deallocate(this%x)
    if (associated(this%y)) deallocate(this%y)
    if (associated(this%z)) deallocate(this%z)
  end subroutine deallocateNexdWavefieldPoints
!-----------------------------------------------------------------------------
!> \brief Get total number of wavefield points
!
  function getNtotNexdWavefieldPoints(this) result(ntot)
    type (nexd_wavefield_points), intent(in) :: this
    integer :: ntot
    if(associated(this%x)) then
       ntot = size(this%x)
    else
       ntot = 0
    endif
  end function getNtotNexdWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the specfem wavefield poinst
! get arrays of wavefield points as are, inversion grid must deal with them
  subroutine getNexdWavefieldPoints(this,x,y,z)
    type (nexd_wavefield_points), intent(in) :: this
    real, dimension(:), pointer :: x,y,z
!
    nullify(x,y,z)
    if(associated(this%x)) then
       allocate(x(size(this%x)),y(size(this%x)),z(size(this%x)))
       x = this%x
       y = this%y
       z = this%z
    end if
  end subroutine getNexdWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the unit factor of nexd wavefield point coordinates
  subroutine getUnitFactorNexdWavefieldPoints(this,uf_wp,errmsg)
    type (nexd_wavefield_points), intent(in) :: this
    real :: uf_wp
    type (error_message) :: errmsg
!
    ! wavefield points were not yet read from file, do not know the unit factor
    if(associated(this%x)) then
       ! NEXD assumes SI units, so unit factor of coordinates of wavefield points is always 1.0
       uf_wp = 1.0
    else
       ! return erroneous unit factor and raise error
       uf_wp = 0.0
       call add(errmsg,2,"wavefield points were not yet read from file, the unit factor is unknown",&
            "getUnitFactorNexdWavefieldPoints")
    end if
  end subroutine getUnitFactorNexdWavefieldPoints
!
end module nexdWavefieldPoints
