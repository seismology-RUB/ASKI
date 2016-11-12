! ===============================================================================
!  GEMINI specific module for ASKI generic module wavefieldPoints.f90
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!-----------------------------------------------------------------------------
!   Provides GEMINI specific interface to routines required by the
!   ASKI wavefieldPoints module to construct wavefield points, store them, and later
!   access them when computing waveform sensitivity kernels.
!-----------------------------------------------------------------------------
module geminiWavefieldPoints
    use chunkCubedSphere
    use fileUnitHandler
    use inputParameter
    use externalRadialNodes
    implicit none
    interface dealloc; module procedure deallocGeminiWavefieldPoints; end interface
    type gemini_wavefield_points
       private
       integer :: nwp                                    ! total number of wavefield points
       real, dimension(:), pointer :: xg => null()       ! x-values in sphere scaled with rearth 
       real, dimension(:), pointer :: yg => null()       ! y-values in sphere scaled with rearth
       real, dimension(:), pointer :: zg => null()       ! z-values in sphere scaled with rearth
    end type gemini_wavefield_points
!
 contains
!-----------------------------------------------------------------------
!  Read parameters (not points) for wavefield points and construct them
!
    subroutine readGeminiWavefieldPoints(this,fuh,filename,errmsg)
    type (gemini_wavefield_points) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len = 25) :: myname = 'readGeminiWavefieldPoints'
    integer :: nlat,nlon,lu,i,nnod,ng
    real :: clon,clat,gam,wlon,wlat
    double precision, dimension(:), pointer :: rnod => null()
    double precision :: rearth
    real, dimension(:), pointer :: xg1,yg1,zg1
    type (chunk_cubed_sphere) :: ccs
    type (external_radial_nodes) :: exnodes
    type (input_parameter) :: inpar
    character (len=80), dimension(12) :: para_keywords
    data para_keywords/'WAVEFIELD_POINTS_CLON','WAVEFIELD_POINTS_CLAT',&
          & 'WAVEFIELD_POINTS_WLON','WAVEFIELD_POINTS_WLAT',&
          & 'WAVEFIELD_POINTS_NLON','WAVEFIELD_POINTS_NLAT',&
          & 'WAVEFIELD_POINTS_ROT',&
          & 'EXTERNAL_NODES_NBLOCKS','EXTERNAL_NODES_NNOD','EXTERNAL_NODES_DR',&
          & 'EXTERNAL_NODES_SHIFT','EXTERNAL_NODES_REARTH'/
!-----------------------------------------------------------------------------------
    call addTrace(errmsg,myname)
    call createKeywordsInputParameter(inpar,para_keywords)
    lu = get(fuh)
    call readSubroutineInputParameter(inpar,lu,filename,errmsg)
    if (.level.errmsg == 2) then
       call dealloc(inpar); call add(fuh,lu); return
    endif
    call add(fuh,lu)
!
!  read out parameters into variables
!
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
    call createChunkCubedSphere(ccs,clon,clat,gam,wlon,wlat,nlon,nlat,errmsg)
    if (.level.errmsg == 2) then
       call dealloc(inpar); return
    endif
    ng = .ng.ccs
!
!  create external radial nodes (unit is km)
!
    call createExternalRadialNodes(exnodes,lu,filename,errmsg)
    if (.level.errmsg == 2) then
       call dealloc(inpar); call dealloc(ccs); return
    endif
!
!  calculate depths, ordered bottom up
!
    nnod = getNnodExternalRadialNodes(exnodes)
    rearth = getRearthExternalRadialNodes(exnodes)
    rnod => getPointerDoubleRadiiExternalRadialNodes(exnodes)
!
!  construct wavefield points
!
    this%nwp = ng*nnod
    allocate(this%xg(ng*nnod),this%yg(ng*nnod),this%zg(ng*nnod))
    xg1 => getXArrayChunkCubedSphere(ccs)
    yg1 => getYArrayChunkCubedSphere(ccs)
    zg1 => getZArrayChunkCubedSphere(ccs)
    do i = 1,nnod
       this%xg((i-1)*ng+1:i*ng) = xg1*rnod(i)
       this%yg((i-1)*ng+1:i*ng) = yg1*rnod(i)
       this%zg((i-1)*ng+1:i*ng) = zg1*rnod(i)
    enddo
!
    call dealloc(ccs)
    call dealloc(exnodes)
    call dealloc(inpar)
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
!  gemini wavefield point coordinates have unit km , so return 1000 here
!  
    subroutine getUnitFactorGeminiWavefieldPoints(this,uf_wp,errmsg)
    type (gemini_wavefield_points), intent(in) :: this
    real :: uf_wp
    type (error_message) :: errmsg
!
    if (associated(this%xg)) then
       uf_wp = 1.0e3
    else
       call add(errmsg,2,"wavefield points were not yet read from file, the unit factor is unknown",&
            "getUnitFactorGeminiWavefieldPoints")
       uf_wp = 0.0
    endif
    end subroutine getUnitFactorGeminiWavefieldPoints
!
 end module
