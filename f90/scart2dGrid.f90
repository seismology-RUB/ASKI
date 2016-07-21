!--------------------------------------------------------------------------
!	Copyright 2015 Wolfgang Friederich
!
!	This file is part of Gemini II.
!
!	Gemini II is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 2 of the License, or
!	any later version.
!
!	Gemini II is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!--------------------------------------------------------------------------------
!  Module to implement a regular discretization on a plane
!--------------------------------------------------------------------------------
 module scart2dGrid
    use mathConstants
    use flexibleType
    use errorMessage
    implicit none
    interface dealloc; module procedure deallocScart2dGrid; end interface
    interface operator (.ng.); module procedure getNGridPointsScart2dGrid; end interface
    interface operator (.nx.); module procedure getNXScart2dGrid; end interface
    interface operator (.ny.); module procedure getNYScart2dGrid; end interface
    interface operator (.wx.); module procedure getWXScart2dGrid; end interface
    interface operator (.wy.); module procedure getWYScart2dGrid; end interface
    interface operator (.gam.); module procedure getGamScart2dGrid; end interface
    interface operator (.xs.); module procedure getXArrayScart2dGrid; end interface
    interface operator (.ys.); module procedure getYArrayScart2dGrid; end interface
    type scart2d_grid
        private
        integer :: ny                                 ! number of grid points y-direction
        integer :: nx                                 ! number of gridpoints x-direction
        real :: gam                                   ! counterclockwise rotation of grid in degrees
        real :: wx,wy                                 ! extension of grid (m)
        real, dimension(:), pointer :: xs => null()   ! x pointing south (m)
        real, dimension(:), pointer :: ys => null()   ! y pointing east (m)
    end type
!
 contains
!----------------------------------------------------------------------------------
!> \brief Create dummy chunk_cubed_sphere object for cartesian applications
!> \param gam Counterclockwise rotation of grid in degrees
!> \param width_x Width in x-direction in meters
!> \param width_y Width in y-direction in meters 
!> \param nx Number of grid points in x-direction
!> \param ny Number of grid points in y-direction
!
    subroutine createScart2dGrid(this,gam,wx,wy,nx,ny,errmsg)
    type (scart2d_grid) :: this
    type (error_message) :: errmsg
    real :: wx,wy,gam
    integer :: nx,ny,ng,i,j,cnt
    real :: dx,dy,x,y
    real, dimension(3,3) :: tm
    character (len=17) :: myname = 'createScart2dGrid'
!
    call addTrace(errmsg,myname)
    if (nx < 2 .or. ny < 2) then
        call add(errmsg,2,'You need at least two grid points per dimension',myname)
        return
    endif
!
    this%gam = gam
    this%wy = wy; this%wx = wx
    this%ny = ny; this%nx = nx
!
    dx = wx/(nx-1)                                          ! width / number of intervals
    dy = wy/(ny-1)                                          ! width / number of intervals
    ng = nx*ny                                              ! total number of grid points
!
    allocate(this%xs(ng),this%ys(ng))
    call rotationMatrixLocalToGlobalScart2dGrid(this,tm)
    cnt = 0
    do i = 1,ny
        y = -0.5*wy+(i-1)*dy
        do j = 1,nx
            cnt = cnt+1
            x = -0.5*wx+(j-1)*dx
            this%xs(cnt) = x*tm(1,1)+y*tm(2,1)
            this%ys(cnt) = x*tm(1,2)+y*tm(2,2)
        enddo
    enddo
    end subroutine createScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Deallocate scart2dGrid
!
    subroutine deallocScart2dGrid(this)
    type (scart2d_grid) :: this
    if (associated(this%xs)) deallocate(this%xs)
    if (associated(this%ys)) deallocate(this%ys)
    end subroutine deallocScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Rotation matrix R to transform coordinates when gam != 0
!!  h_n = R_nk e_k, where h_n local basis vectors and e_k global ones
!
    subroutine rotationMatrixLocalToGlobalScart2dGrid(this,tm)
    type (scart2d_grid), intent(in) :: this
    real, dimension(:,:) :: tm
    real :: cosgam,singam
!
    cosgam = cos(this%gam*mc_deg2rad)
    singam = sin(this%gam*mc_deg2rad)
    tm(1,1) = cosgam
    tm(1,2) = singam
    tm(2,1) = -singam
    tm(2,2) = cosgam
    end subroutine rotationMatrixLocalToGlobalScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get total number of grid points
!
    function getNGridPointsScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    integer :: res
    res = this%nx*this%ny
    end function getNGridPointsScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get nx
!
    function getNXScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    integer :: res
    res = this%nx
    end function getNXScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get ny
!
    function getNYScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    integer :: res
    res = this%ny
    end function getNYScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get wx
!
    function getWXScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    real :: res
    res = this%wx
    end function getWXScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get wy
!
    function getWYScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    real :: res
    res = this%wy
    end function getWYScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get xs-array
!
    function getXArrayScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%xs
    end function getXArrayScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get ys-array
!
    function getYArrayScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%ys
    end function getYArrayScart2dGrid
!-----------------------------------------------------------------------------------
!> \brief Get xs-array
!
    function getGamScart2dGrid(this) result(res)
    type (scart2d_grid), intent(in) :: this
    real :: res
    res = this%gam
    end function getGamScart2dGrid
!
 end module
