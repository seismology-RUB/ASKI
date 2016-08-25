!----------------------------------------------------------------------------
!	Copyright 2016 Wolfgang Friederich
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
!  Module to implement a discretization on the surface of the unit sphere
!  using the "cubed-sphere" algorithm (see explanations in cubedSphere.tex)
!--------------------------------------------------------------------------------
 module chunkCubedSphere
    use mathConstants
    use flexibleType
    use errorMessage
    implicit none
!    interface new
!        module procedure createCartesianChunkCubedSphere
!        module procedure createSphericalChunkCubedSphere
!    end interface
    interface dealloc; module procedure deallocChunkCubedSphere; end interface
    interface operator (.ng.); module procedure getNGridPointsChunkCubedSphere; end interface
    interface operator (.nlat.); module procedure getNlatChunkCubedSphere; end interface
    interface operator (.nlon.); module procedure getNlonChunkCubedSphere; end interface
    interface operator (.wlat.); module procedure getWlatChunkCubedSphere; end interface
    interface operator (.wlon.); module procedure getWlonChunkCubedSphere; end interface
    interface operator (.clat.); module procedure getClatChunkCubedSphere; end interface
    interface operator (.clon.); module procedure getClonChunkCubedSphere; end interface
    interface operator (.gam.); module procedure getGamChunkCubedSphere; end interface
    interface operator (.xg.); module procedure getXArrayChunkCubedSphere; end interface
    interface operator (.yg.); module procedure getYArrayChunkCubedSphere; end interface
    interface operator (.zg.); module procedure getZArrayChunkCubedSphere; end interface
    interface operator (.th.); module procedure getThArrayChunkCubedSphere; end interface
    interface operator (.ph.); module procedure getPhArrayChunkCubedSphere; end interface
    type chunk_cubed_sphere
        private
        integer :: nlon                               ! number of grid points in E/W-direction, y-direction
        integer :: nlat                               ! number of gridpoints in N/S-direction, x-direction
        real :: ctheta,clon                           ! center of grid in degrees, (0,0) for cartesian applications
        real :: gam                                   ! counterclockwise rotation of grid in degrees
        real :: wlat,wlon                             ! extension of grid in degrees  (NS,x and WE,y for gam = 0)
        real, dimension(:), pointer :: xg => null()   ! x on global unit sphere 
        real, dimension(:), pointer :: yg => null()   ! y on global unit sphere
        real, dimension(:), pointer :: zg => null()   ! z on global unit sphere
        real, dimension(:), pointer :: th => null()   ! theta in rad on global unit sphere 
        real, dimension(:), pointer :: ph => null()   ! phi in rad on global unit sphere
        real, dimension(:,:), pointer :: tm => null() ! transformation matrix from global to local basis vectors
    end type
!
 contains
!----------------------------------------------------------------------------------
!> \brief Create chunk_cubed_sphere object for spherical applications
!> \param clon Longitude of center of grid area in degrees
!> \param clat Latitude of center of grid area in degrees
!> \param gam Counterclockwise rotation of grid in degrees
!> \param ang_width_lon Angular width in longitudinal direction in degrees 
!> \param ang_width_lat Angular width in latitudinal direction in degrees 
!> \param nlon Number of grid points in longitudinal direction
!> \param nlat Number of grid points in latitudinal direction
!! use equiangular grid and not equidistant grid on tangential plane
!
    subroutine createChunkCubedSphere(this,clon,clat,gam,wlon,wlat,nlon,nlat,errmsg)
    type (chunk_cubed_sphere) :: this
    type (error_message) :: errmsg
    real :: clat,clon,gam,wlon,wlat
    integer :: nlon,nlat,ng,i,j,k,cnt
    character (len=31) :: myname = 'createSphericalChunkCubedSphere'
    real :: ang_width_lat,ang_width_lon,x,y,dlat,dlon,d,xg,yg,zg
    real, dimension(3) :: xs
    real, dimension(3,3) :: tm
!
    call addTrace(errmsg,myname)
    if (nlat < 2 .or. nlon < 2) then
        call add(errmsg,2,'You need at least two grid points per dimension',myname)
        return
    endif                        
!
    this%ctheta = 90.-clat; this%clon = clon; this%gam = gam
    this%wlon = wlon; this%wlat = wlat
    this%nlon = nlon; this%nlat = nlat
    call rotationMatrixLocalToGlobalChunkCubedSphere(this,tm)
    allocate(this%tm(3,3))
    this%tm = tm
!
    ang_width_lat = wlat*mc_deg2rad
    ang_width_lon = wlon*mc_deg2rad
!
    dlat = ang_width_lat/(nlat-1)                           ! width / number of intervals
    dlon = ang_width_lon/(nlon-1)                           ! width / number of intervals
    ng = nlon*nlat                                          ! total number of grid points
!
    allocate(this%xg(ng),this%yg(ng),this%zg(ng),this%th(ng),this%ph(ng))
    cnt = 0
    do i = 1,this%nlon                                      ! global numbering second in y/east-direction
        y = tan(-0.5*ang_width_lon+(i-1)*dlon)
        do j = 1,this%nlat                                  ! global numbering first in x/south-direction
            cnt = cnt+1
            x = tan(-0.5*ang_width_lat+(j-1)*dlat)
            d = sqrt(1.+y*y+x*x)                            ! distance of point on tangential plane to origin
            xs(1) = x/d; xs(2) = y/d; xs(3) = 1.0/d         ! projection onto unit sphere: x/d = xs/1
            ! the xs(j) coordinates refer to basis vectors h(i) with h(3) || radial direction at (clat,clon),
            ! h(1) || S-direction at (clat,clon) for gam=0 and h(2) || local east direction for gam=0.
            ! If gam /= 0, h(1) and h(2) are rotated around the local radial direction. In cubedSphere.tex
            ! equations are derived that give expressions for the h(i) in terms of global cartesion basis vectors e(n)
            ! where e(3) points towards the north pole, e(1) towards equatior at lon=0 and e(2) towards the equator at lon=90.
            ! The relation is h(i) = R(i,n) e(n). Thus, a grid point described by radius vector p can be expressed as
            ! p = xs(i) h(i) = xg(n) e(n). With h(i) = R(i,n)e(n) we find: p = xs(i)R(i,n)e(n) = xg(n)e(n) and therefore
            ! xg(n) = xs(i)R(i,n).
            xg = 0.; yg = 0.; zg = 0. 
            do k = 1,3
                xg = xg+xs(k)*tm(k,1)        ! transformation into global cartesian system
                yg = yg+xs(k)*tm(k,2)        ! with x-axis pointing through (lat=0,lon=0)
                zg = zg+xs(k)*tm(k,3)        ! z-axis through north pole and y-axis through (0,90)
            enddo
            this%xg(cnt) = xg
            this%yg(cnt) = yg
            this%zg(cnt) = zg
            this%th(cnt) = acos(zg)
            this%ph(cnt) = atan2(yg,xg)
        enddo
    enddo
    end subroutine createChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Deallocate ChunkCubedSphere
!
    subroutine deallocChunkCubedSphere(this)
    type (chunk_cubed_sphere) :: this
    if (associated(this%xg)) deallocate(this%xg)
    if (associated(this%yg)) deallocate(this%yg)
    if (associated(this%zg)) deallocate(this%zg)
    if (associated(this%th)) deallocate(this%th)
    if (associated(this%ph)) deallocate(this%ph)
    if (associated(this%tm)) deallocate(this%tm)
    end subroutine deallocChunkCubedSphere
!---------------------------------------------------------------------------------
!> \brief Get essential data of chunkCubedSphere in a flexible array
!
    subroutine getEssentialDataChunkCubedSphere(this,ft)
    type (chunk_cubed_sphere), intent(in) :: this
    type (flexible), dimension(:) :: ft
    ft(1) = 'S'; ft(2) = this%nlon; ft(3) = this%nlat
    ft(4) = this%clon; ft(5) = 90.-this%ctheta; ft(6) = this%gam
    ft(7) = this%wlat; ft(8) = this%wlon
    end subroutine getEssentialDataChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get total number of grid points
!
    function getNGridPointsChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    integer :: res
    res = this%nlat*this%nlon
    end function getNGridPointsChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get nlat
!
    function getNlatChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    integer :: res
    res = this%nlat
    end function getNlatChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get nlon
!
    function getNlonChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    integer :: res
    res = this%nlon
    end function getNlonChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get wlat
!
    function getWlatChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real :: res
    res = this%wlat
    end function getWlatChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get wlon
!
    function getWlonChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real :: res
    res = this%wlon
    end function getWlonChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get clat
!
    function getClatChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real :: res
    res = 90.-this%ctheta
    end function getClatChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get clon
!
    function getClonChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real :: res
    res = this%clon
    end function getClonChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get gamma
!
    function getGamChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real :: res
    res = this%gam
    end function getGamChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get pointer to array of xg-values
!
    function getXArrayChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%xg
    end function getXArrayChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get pointer to array of yg-values
!
    function getYArrayChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%yg
    end function getYArrayChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get pointer to array of zg-values
!
    function getZArrayChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%zg
    end function getZArrayChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get pointer to array of th-values
!
    function getThArrayChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%th
    end function getThArrayChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get pointer to array of zg-values
!
    function getPhArrayChunkCubedSphere(this) result(res)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:), pointer :: res
    res => this%ph
    end function getPhArrayChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Rotation matrix R as in cubedSphere.tex
!!  h_n = R_nk e_k, where h_n local basis vectors and e_k global ones
!
    subroutine rotationMatrixLocalToGlobalChunkCubedSphere(this,tm)
    type (chunk_cubed_sphere), intent(in) :: this
    real, dimension(:,:) :: tm
    real :: costheta,cosphi,cosgam,sintheta,sinphi,singam
!
    cosgam = cos(this%gam*mc_deg2rad)
    singam = sin(this%gam*mc_deg2rad)
    costheta = cos(this%ctheta*mc_deg2rad)
    cosphi = cos(this%clon*mc_deg2rad)
    sintheta = sin(this%ctheta*mc_deg2rad)
    sinphi = sin(this%clon*mc_deg2rad)
!
    tm(1,1) = cosgam*costheta*cosphi-singam*sinphi
    tm(1,2) = cosgam*costheta*sinphi+singam*cosphi
    tm(1,3) = -cosgam*sintheta
    tm(2,1) = -singam*costheta*cosphi-cosgam*sinphi
    tm(2,2) = -singam*costheta*sinphi+cosgam*cosphi
    tm(2,3) = singam*sintheta
    tm(3,1) = sintheta*cosphi
    tm(3,2) = sintheta*sinphi
    tm(3,3) = costheta
    end subroutine rotationMatrixLocalToGlobalChunkCubedSphere
!-----------------------------------------------------------------------------------
!> \brief Get cell index of some point on sphere
!!  x,y,z: cartesian coordinates on unit sphere
!!
!!  Note: Cells are areas defined by four grid points at the corners
!!  Example:
!!        ------------------------------->  y,ip,lon
!!   |
!!   |    +-----+-----+-----+-----+-----+
!!   |    |  1  |  3  |  5  |  7  |  9  |
!!   |    +-----+-----+-----+-----+-----+       
!!   V    |  2  |  4  |  6  |  8  | 10  |
!!  x,jp, +-----+-----+-----+-----+-----+
!!  lat
!!
!!  Numbers in the cells are cell indices.
!!  Global cell index: jgl = (ip-1)*(nx-1)+jp 
!!  where ip,jp are point indices at the upper left corner of the cell.
!!  For given global index jgl, the cell row in y-direction is:
!!  ry = (jgl-1)/(nx-1)+1 
!!  and the global point index of the upper left corner is:
!!  iul = jgl+ry-1  
!----------------------------------------------------------
    function getCellIndexChunkCubedSphere(this,x,y,z) result(jgl)
    type (chunk_cubed_sphere) :: this
    real :: x,y,z
    integer :: k,jgl,jp,ip
    real :: xig,yig,ang_width_lat,ang_width_lon,dlat,dlon
    real, dimension(3) :: xl
!
!  transform to local spherical coordinates with z-axis pointing 
!  through clon,clat
!
    do k = 1,3
       xl(k) = this%tm(k,1)*x+this%tm(k,2)*y+this%tm(k,3)*z
    enddo
!
!  point not on spherical chunk (negative local z)
!
    if (xl(3) < 0.0) then
       jgl = -1; return
    endif
!
!  project back to cube face on tangential plane
!
    xig = xl(1)/xl(3)
    yig = xl(2)/xl(3)
!
    ang_width_lat = this%wlat*mc_deg2rad
    ang_width_lon = this%wlon*mc_deg2rad
    dlat = ang_width_lat/(this%nlat-1)
    dlon = ang_width_lon/(this%nlon-1)
!
!  point inside grid ?
!
    if (xig < -0.5*ang_width_lat .or. xig >= 0.5*ang_width_lat) then
        jgl = -1; return
    endif 
    if (yig < -0.5*ang_width_lon .or. yig >= 0.5*ang_width_lon) then
        jgl = -1; return
    endif 
!
!  indices of upper left grid point (see figure) and cell index
!
    jp = int((atan(xig)+0.5*ang_width_lat)/dlat)+1
    ip = int((atan(yig)+0.5*ang_width_lon)/dlon)+1
    jgl = (ip-1)*(this%nlat-1)+jp
    end function getCellIndexChunkCubedSphere
!------------------------------------------------------------------------------------
!> \brief return indices of corners of next cell in grid
!!        array with four indices of corner points
!
    function nextCellChunkCubedSphere(this,idx) result(next)
    type (chunk_cubed_sphere) :: this
    integer, dimension(:) :: idx
    logical :: next
    integer :: ry,ill
    integer :: cnt = 0
    save :: cnt
!
    cnt = cnt+1
    ry = (cnt-1)/(this%nlat-1)+1                          ! number of row in y-direction
    if (ry > this%nlon-1) then
        next = .false.
        cnt = 0
        return
    endif
!
    ill = cnt+ry-1                                      ! index of lower-left point of cell
    idx = (/ ill,ill+1,ill+1+this%nlat,ill+this%nlat /)
    next = .true.
    end function nextCellChunkCubedSphere
!
 end module
