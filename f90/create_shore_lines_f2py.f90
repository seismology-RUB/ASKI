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
!> \brief fortran module to generate create_shore_lines_f2py.so by f2py for import of create_shore_lines_f2py.shore_f2py in create_shore_lines.py
!!
!! \details This is a wrapper about the inversion grid module that can be compiled
!!  with f2py and is imported in create_gmt_shore_lines.py in order to interface with
!!  the ASKI inversion grid. 
!!
!! \author Florian Schumacher
!! \date May 2016
!! 
!! \copyright &copy; GNU Public License
!
module shore_f2py

  use create_shore_lines_f2pyVar

contains
!
!---------------------------------------------------------------------------------
!> \brief create inversion grid
!! \param type_name character string indicating inversion grid type
!! \param parfile character string containing filename of inversion grid parameter file
!! \param path character string containing path for local invgrid files
!! \param recreate logical indicating whether to recreate inversion grid or not
!! \param istatus integer containing status on return (0 = success, 1 = error)
!
  subroutine set_up_invgrid(type_name,parfile,path,recreate,istatus)
    ! incoming
    character(len=*), intent(in) :: type_name,parfile,path
    logical, intent(in) :: recreate
    ! returning
    integer, intent(out) :: istatus
    ! local
    integer :: lu
!
    call new(errmsg,'set_up_invgrid')
    lu = 11
    call createInversionGrid(invgrid,type_name,parfile,path,lu,errmsg,recreate)
    if(.level.errmsg /= 0) call print(errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! at this point, everything went alright, so indicate success
    istatus = 0
    invgrid_created = .true.
    return
!
1   call dealloc(errmsg)
    ! indicate error
    istatus = 1
  end subroutine set_up_invgrid
!
!---------------------------------------------------------------------------------
!> \brief deallocate inversion grid
!
  subroutine deallocate_invgrid()
    if(invgrid_created) call dealloc(invgrid)
  end subroutine deallocate_invgrid
!
!---------------------------------------------------------------------------------
!> \brief ask the inversion grid, whether a point on Earth's surface is located laterally inside of the inversion grid's bounds
!! \param lat latitude (deg) of point on surface to be ckecked
!! \param lon longitude (deg) of point on surface to be checked
!! \param istatus integer indicating whether the point is inside ( = 0 , meaning success) or outside ( = 1 , meaning error)
!! \param ichunk if istatus==0 AND invgrid type is chunksInversionGrid, then ichunk contains the chunk index of that point; otherwise ichunk is -1
!
  subroutine point_is_inside_inversion_grid(lat,lon,istatus,ichunk)
    ! incoming
    real, intent(in) :: lat,lon
    ! returning
    integer, intent(out) :: istatus,ichunk
    ! local
    logical :: is_inside
    real :: c1,c2,c3
!
    c1 = lat
    c2 = lon
    c3 = 0
!    
    is_inside = pointInsideInversionGrid(invgrid,c1,c2,c3,'station',ichunk=ichunk)
!
    if(is_inside) then
       ! istatus = 0 means success, i.e. point is inside 
       istatus = 0
    else
       istatus = 1
       ! if you return istatus=1, then also set ichunk = -1 (no sensible value to return for ichunk)
       ichunk = -1 
    end if
  end subroutine point_is_inside_inversion_grid
!
!---------------------------------------------------------------------------------
!> \brief transform coordinates to vtk projection
!! \param lat 
!! \param n bound of rank-1 arrays lat,lon,x_vtk,y_vtk,z_vtk . Without indicating intend, f2py makes the argument n optional. This way it works to return explicit shapes back to python (compare http://scicomp.stackexchange.com/a/6987)
!
  subroutine transform_to_vtk_projection(lat,lon,istatus,x_vtk,y_vtk,z_vtk,n)
    ! optional
    integer :: n
    ! incoming
    real, intent(in), dimension(n) :: lat,lon
    ! returning
    integer, intent(out) :: istatus
    real, intent(out), dimension(n) :: x_vtk,y_vtk,z_vtk
    ! local
    real, dimension(n) :: c1,c2,c3
!
    call new(errmsg,'transform_to_vtk_projection')
    c1 = lat
    c2 = lon
    c3 = 0
    call transformToVtkInversionGrid(invgrid,c1,c2,c3,'station',errmsg)
    if(.level.errmsg /= 0) call print(errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! at this point, everything went alright, so return correct values and indicate success
    x_vtk = c1
    y_vtk = c2
    z_vtk = c3
    istatus = 0
    return
!
1   call dealloc(errmsg)
    ! indicate error
    x_vtk = 0
    y_vtk = 0
    z_vtk = 0
    istatus = 1
  end subroutine transform_to_vtk_projection
end module shore_f2py
