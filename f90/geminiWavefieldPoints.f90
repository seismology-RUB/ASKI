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
!-----------------------------------------------------------------
!> \brief Module dealing with Gemini wavefield points
!! \details For cartesian applications, associate x-coordinate with lat
!! and y-coordinate with lon whne specifying nx,ny and wx and wy
!----------------------------------------------------------------
 module geminiWavefieldPoints
    use chunkCubedSphere
    use scart2dGrid
    use fileUnitHandler
    use inputParameter
    use externalRadialNodes
    implicit none
    interface dealloc; module procedure deallocGeminiWavefieldPoints; end interface
    interface operator (.csys.); module procedure getCoordinateSystemGeminiWavefieldPoints; end interface
    interface operator (.ntot.); module procedure getNtotGeminiWavefieldPoints; end interface
    type gemini_wavefield_points
        private
         character (len=1) :: csys
         integer :: nwp                                    ! total number of wavefield points
         real, dimension(:), pointer :: xg => null()       ! x-values in sphere scaled with rearth or on halfspace
         real, dimension(:), pointer :: yg => null()       ! y-values in sphere scaled with rearth or on halfspace
         real, dimension(:), pointer :: zg => null()       ! z-values in sphere scaled with rearth or in halfspace
    end type
!
 contains
!-----------------------------------------------------------------------
!> \brief Read in parameters (not points) for wavefield points
!! Construct wavefield points
!
    subroutine readGeminiWavefieldPoints(this,fuh,filename,errmsg)
    type (gemini_wavefield_points) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    ! local
    character (len = 25) :: myname = 'readGeminiWavefieldPoints'
    integer :: nlat,nlon,lu,i,nz,ng
    real :: clon,clat,gam,wlon,wlat
    double precision, dimension(:), pointer :: rnod => null()
    double precision :: rearth
    real, dimension(:), pointer :: xg1,yg1,zg1
    type (chunk_cubed_sphere) :: ccs
    type (scart2d_grid) :: scg
    type (external_radial_nodes) :: exnodes
    type (input_parameter) :: inpar
    character (len=80), dimension(13) :: para_keywords
    data para_keywords/'WAVEFIELD_POINTS_CLON','WAVEFIELD_POINTS_CLAT',&
          & 'WAVEFIELD_POINTS_WLON','WAVEFIELD_POINTS_WLAT',&
          & 'WAVEFIELD_POINTS_NLON','WAVEFIELD_POINTS_NLAT',&
          & 'WAVEFIELD_POINTS_ROT',&
          & 'EXTERNAL_NODES_NBLOCKS','EXTERNAL_NODES_NNOD','EXTERNAL_NODES_DR',&
          & 'EXTERNAL_NODES_SHIFT','EXTERNAL_NODES_REARTH',&
          & 'COORDINATE_SYSTEM'/
!-----------------------------------------------------------------------------------
    lu = get(fuh)
    call addTrace(errmsg,myname)
    call createKeywordsInputParameter(inpar,para_keywords)
    call readSubroutineInputParameter(inpar,lu,filename,errmsg)
    if (.level.errmsg == 2) goto 1
!
!  read out parameters into variables
!
    this%csys = sval(inpar,'COORDINATE_SYSTEM')
    nlat = ival(inpar,'WAVEFIELD_POINTS_NLAT')
    nlon = ival(inpar,'WAVEFIELD_POINTS_NLON')
    clat = rval(inpar,'WAVEFIELD_POINTS_CLAT')
    clon = rval(inpar,'WAVEFIELD_POINTS_CLON')
    wlat = rval(inpar,'WAVEFIELD_POINTS_WLAT')
    wlon = rval(inpar,'WAVEFIELD_POINTS_WLON')
    gam  = rval(inpar,'WAVEFIELD_POINTS_ROT')
!
!  create chunk cubed sphere object
!
    if (this%csys == 'S') then
        call createChunkCubedSphere(ccs,clon,clat,gam,wlon,wlat,nlon,nlat,errmsg)
        ng = .ng.ccs
    else if (this%csys == 'C') then
        call createScart2dGrid(scg,gam,wlat,wlon,nlat,nlon,errmsg)
        ng = .ng.scg
    else 
        call add(errmsg,2,'invalid coordinate system',myname)
    endif
    if (.level.errmsg == 2) goto 1
!
!  create external radial nodes (unit is km)
!
    call createExternalRadialNodes(exnodes,lu,filename,errmsg)
    if (.level.errmsg == 2) goto 1
!
!  calculate depths, ordered bottom up
!
    nz = getNnodExternalRadialNodes(exnodes)
    rearth = getRearthExternalRadialNodes(exnodes)
    call getRadiiExternalRadialNodes(exnodes,rnod,errmsg)
    if (.level.errmsg == 2) goto 1
!
!  calculate wavefield points
!
    this%nwp = ng*nz
    allocate(this%xg(ng*nz),this%yg(ng*nz),this%zg(ng*nz))
    if (this%csys == 'S') then
        xg1 => getXArrayChunkCubedSphere(ccs)
        yg1 => getYArrayChunkCubedSphere(ccs)
        zg1 => getZArrayChunkCubedSphere(ccs)
        do i = 1,nz
            this%xg((i-1)*ng+1:i*ng) = xg1*rnod(i)
            this%yg((i-1)*ng+1:i*ng) = yg1*rnod(i)
            this%zg((i-1)*ng+1:i*ng) = zg1*rnod(i)
        enddo
    else if (this%csys == 'C') then
        xg1 => getXArrayScart2dGrid(scg)
        yg1 => getYArrayScart2dGrid(scg)
        do i = 1,nz
            this%xg((i-1)*ng+1:i*ng) = xg1
            this%yg((i-1)*ng+1:i*ng) = yg1
            this%zg((i-1)*ng+1:i*ng) = real(rearth-rnod(i))*1.e3  ! convert to meters
        enddo
    else
        call add(errmsg,2,'invalid coordinate system',myname)
    endif
!
1   call dealloc(ccs); call dealloc(scg)
    call dealloc(exnodes)
    call dealloc(inpar)
    if (associated(rnod)) deallocate(rnod)
    call add(fuh,lu)
    end subroutine readGeminiWavefieldPoints
!--------------------------------------------------------------------------
!> \brief Deallocate gemini_wavefield_points
!
    subroutine deallocGeminiWavefieldPoints(this)
    type (gemini_wavefield_points ) :: this
    if (associated(this%xg)) deallocate(this%xg)
    if (associated(this%xg)) deallocate(this%yg)
    if (associated(this%zg)) deallocate(this%zg)
    end subroutine deallocGeminiWavefieldPoints
!--------------------------------------------------------------------------------
!> \brief Get coordinate system
!
    function getCoordinateSystemGeminiWavefieldPoints(this) result(res)
    type (gemini_wavefield_points), intent(in) :: this
    character (len=1) :: res
    res = this%csys
    end function getCoordinateSystemGeminiWavefieldPoints
!-----------------------------------------------------------------------------
!> \brief get total number of wavefield points
!
    function getNtotGeminiWavefieldPoints(this) result(res)
    type (gemini_wavefield_points), intent(in) :: this
    integer :: res
    res = this%nwp
    end function getNtotGeminiWavefieldPoints
!---------------------------------------------------------------------
!> \brief get all wavefield points
!! Allocate space here for returned wavefield points
!! Other realizations of wavefieldPoints may not store the points
!! in the object but construct them on the fly here
!! Calling routines should care for deallocation
!
  subroutine getGeminiWavefieldPoints(this,c1,c2,c3)
    type (gemini_wavefield_points), intent(in) :: this
    real, dimension(:), pointer :: c1,c2,c3
    integer :: n
    if(associated(this%xg)) then
       n = size(this%xg)
    else
       n = 0
       nullify(c1,c2,c3)
       return
    end if
    allocate(c1(n),c2(n),c3(n))
    c1 = this%xg
    c2 = this%yg
    c3 = this%zg
  end subroutine getGeminiWavefieldPoints
!---------------------------------------------------------------------------
!> \brief Get the unit factor of gemini wavefield point coordinates
  subroutine getUnitFactorGeminiWavefieldPoints(this,uf_wp,errmsg)
    type (gemini_wavefield_points), intent(in) :: this
    real :: uf_wp
    type (error_message) :: errmsg
    if(associated(this%xg)) then
       if (this%csys == 'S') then
          ! gemini wavefield point coordinates have unit km , so return 1000 here
          uf_wp = 1.0e3
       else if (this%csys == 'C') then
          ! gemini wavefield point coordinates have unit m , so return 1.0 here
          uf_wp = 1.0
       end if
    else
       call add(errmsg,2,"wavefield points were not yet read from file, the unit factor is unknown",&
            "getUnitFactorGeminiWavefieldPoints")
       uf_wp = 0
    end if
  end subroutine getUnitFactorGeminiWavefieldPoints
!
 end module
