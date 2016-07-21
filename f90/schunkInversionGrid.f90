!----------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.0.
!
!   ASKI version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief simple spherical inversion grid using one chunk of a cubed sphere
!!
!! \details The inversion domain that can be covered by a simple chunk inversion grid
!!  is one chunk of a cubed sphere with specific maximum depth, center location and 
!!  rotation about the local vertical. From outside the module (e.g. from module wavefieldPoints)
!!  always global geographic cartesian coordinates in [km] are assumed, with 3-axis pointing towards 
!!  the north pole, 1-axis towards the equator at lon=0 and 2-axis towards the equator at lon=90 degrees.
!!  This coordinate system is referred to as "Global". Internally, the module handles two more sets of 
!!  coordinate system: 
!!  One is referred to as "Local", or "LOCAL_CURV" and is essentially the same as "Global", except
!!  that the chunk is not centered at clat,clon but at the north pole and the chunk is not rotated about 
!!  the local vertical by gam.
!!  The third coordinate system is referred to as "Scart" or "LOCAL_FLAT" which is essentially the same
!!  as "Local", except that the curvature is removed, such that the original chunk now is projected to 
!!  a flat cartesian block in [km] centered at the north pole. 
!!  For sake of vtk plotting ONLY, there are two more coordinate projections which are variants of "LOCAL_CURV"
!!  and "LOCAL_FLAT", namely "LOCAL_NORTH_CURV" and "LOCAL_NORTH_FLAT". The "north" refers to the fact
!!  that in those cases the chunk is kept rotated about the local vertical by gam and this means that the
!!  northing is along a coordinate axis (namely the negative x-axis, i.e. x really points to south in those cases,
!!  which can be nice for plotting instead of the gam-rotation removed)
!!
!! \author Wolfgang Friederich and Florian Schumacher
!! \date Okt 2015
!
module schunkInversionGrid
  !
  use inputParameter
  use scartInversionGrid
  use vectorPointer
  use mathConstants
  use errorMessage
  implicit none
  type schunk_inversion_grid
     private
     logical :: is_defined = .false. !< flag indicating the correct definition (initialization) of the object (i.e. all following values)
     real :: clat,clon !< center of grid on sphere
     real :: rmax      !< radius at top of grid
     real :: gam       !< angle in degrees of counter-clockwise rotation about local z-axis
     real, dimension(3,3) :: tm_global2local   !< h(i) = tm(i,n) e(n) where h(i) are cartesian basis vectors attached to the grid
     real, dimension(2,2) :: tm_gam !< simple rotation of angle gam anti-clockwise about the Z axis
     !< located at clat,clon and rotated by gam and e(n) are geographic cartesian
     !< basis vectors with e(3) pointing towards the north pole, e(1) towards the equator at lon=0
     !< and e(2) towards the equator at lon=90.
     character(len=16) :: vtk_projection = '' !< 'GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT', indicating the vtk coordinates for routines getGeometryVtkSchunkInversionGrid and transformToVtkSchunkInversionGrid
     logical :: apply_vtk_coords_scaling_factor
     real :: vtk_coords_scaling_factor
     type (scart_inversion_grid) :: scart
  end type schunk_inversion_grid
  !
contains
!------------------------------------------------------------------------
!> \brief logical return whether this schunk_inversion_grid is able to transform points (of given coords type) to vtk plot projection
!
  function canTransformToVtkPointsOutsideSchunkInversionGrid(this,coords_type) result(l)
    type(schunk_inversion_grid) :: this
    character(len=*) :: coords_type
    logical :: l
    ! the schunk_inversion_grid has capability to transform any points to vtk,
    ! provided the object is defined
    l = this%is_defined
  end function canTransformToVtkPointsOutsideSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief get unit factor of the volume element
!! \param this schunk inversion grid
!! \param uf_vol unit factor of volume element (return value of this subroutine)
!! \param errmsg error message
!
  subroutine getUnitFactorOfVolumeElementSchunkInversionGrid(this,uf_vol,errmsg)
    type (schunk_inversion_grid) :: this
    real :: uf_vol
    type (error_message) :: errmsg
    character (len=47) :: myname = 'getUnitFactorOfVolumeElementSchunkInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) call add(errmsg,1,"be aware that the inversion grid not yet defined; "//&
         "however, the unit factor of the volume element can be correctly computed at this point",myname)
!
    ! The schunk inversion grid assumes units for its spatial extension in km !
    ! This is independent of the unit of wavefield points
    uf_vol = 1.0e9
  end subroutine getUnitFactorOfVolumeElementSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief create simple chunk inversion grid
!! \details the parameter file given, must contain all necessary parameters to define
!!  an object of this type. 
!! \param this simple chunk inversion grid
!! \param parfile filename of parameter file containing definintion of this inversion grid
!! \param lu file unit to use for reading and writing files
!! \param errmsg error message
!
  subroutine createSchunkInversionGrid(this,parfile,lu,errmsg)
    type(schunk_inversion_grid) :: this
    character (len=*) :: parfile
    integer :: lu,ios,nblock,nlay,ilay,ibl,iz
    integer, dimension(:), pointer :: nlay_block,nlat,nlon
    real :: wlat,wlon,wx,wy
    real :: ctheta,costheta,cosphi,cosgam,sintheta,sinphi,singam
    real, dimension(:), pointer :: thickness
    real, dimension(:), pointer :: z
    type (error_message) :: errmsg
    character (len=400) :: errstr
    character (len=25) :: myname = 'createSchunkInversionGrid'
    type (input_parameter) :: inpar
    character (len=80), dimension(15) :: inpar_keys
    data inpar_keys/'SCHUNK_INVGRID_CLAT','SCHUNK_INVGRID_CLON','SCHUNK_INVGRID_RMAX','SCHUNK_INVGRID_WLAT',&
         'SCHUNK_INVGRID_WLON','SCHUNK_INVGRID_ROT','SCHUNK_INVGRID_NREF_BLOCKS','SCHUNK_INVGRID_NLAY',&
         'SCHUNK_INVGRID_THICKNESS','SCHUNK_INVGRID_NLAT','SCHUNK_INVGRID_NLON','VTK_PROJECTION',&
         'SCALE_VTK_COORDS','VTK_COORDS_SCALING_FACTOR','VTK_GEOMETRY_TYPE'/
    !
    call addTrace(errmsg,myname)
    if(this%is_defined) then
       call add(errmsg,1,"this object is already defined, deallocating it now before creating new one",myname)
       call deallocateSchunkInversionGrid(this)
    end if
    !
    call createKeywordsInputParameter(inpar,inpar_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) goto 1
    !
    !   define clat,clon,rmin,wlat,wlon,gam
    !
    this%clat = rval(inpar,'SCHUNK_INVGRID_CLAT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_CLAT' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_CLAT')//"'",myname)
       goto 1
    end if
    !
    this%clon = rval(inpar,'SCHUNK_INVGRID_CLON',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_CLON' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_CLON')//"'",myname)
       goto 1
    end if
    !
    this%rmax = rval(inpar,'SCHUNK_INVGRID_RMAX',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_RMAX' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_RMAX')//"'",myname)
       goto 1
    end if
    !
    wlat = rval(inpar,'SCHUNK_INVGRID_WLAT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_WLAT' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_WLAT')//"'",myname)
       goto 1
    end if
    if(wlat <= 0) then
       write(errstr,*) "SCHUNK_INVGRID_WLAT = ",wlat," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    !
    wlon = rval(inpar,'SCHUNK_INVGRID_WLON',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_WLON' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_WLON')//"'",myname)
       goto 1
    end if
    if(wlon <= 0) then
       write(errstr,*) "SCHUNK_INVGRID_WLON = ",wlon," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    !
    this%gam = rval(inpar,'SCHUNK_INVGRID_ROT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCHUNK_INVGRID_ROT' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_ROT')//"'",myname)
       goto 1
    end if
    !
    !  rotation matrix
    !
    ctheta = 90.-this%clat
    cosgam = cos(this%gam*mc_deg2rad)
    singam = sin(this%gam*mc_deg2rad)
    costheta = cos(ctheta*mc_deg2rad)
    cosphi = cos(this%clon*mc_deg2rad)
    sintheta = sin(ctheta*mc_deg2rad)
    sinphi = sin(this%clon*mc_deg2rad)
    !
    this%tm_global2local(1,1) = cosgam*costheta*cosphi-singam*sinphi
    this%tm_global2local(1,2) = cosgam*costheta*sinphi+singam*cosphi
    this%tm_global2local(1,3) = -cosgam*sintheta
    this%tm_global2local(2,1) = -singam*costheta*cosphi-cosgam*sinphi
    this%tm_global2local(2,2) = -singam*costheta*sinphi+cosgam*cosphi
    this%tm_global2local(2,3) = singam*sintheta
    this%tm_global2local(3,1) = sintheta*cosphi
    this%tm_global2local(3,2) = sintheta*sinphi
    this%tm_global2local(3,3) = costheta
    !
    this%tm_gam(1,1) = cosgam
    this%tm_gam(1,2) = -singam
    this%tm_gam(2,1) = singam
    this%tm_gam(2,2) = cosgam
    !
    ! define nblock,nlay_block,nlay,z
    !
    nblock = ival(inpar,'SCHUNK_INVGRID_NREF_BLOCKS',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read integer value for 'SCHUNK_INVGRID_NREF_BLOCKS' from '"//&
            trim(inpar.sval.'SCHUNK_INVGRID_NREF_BLOCKS')//"'",myname)
       goto 1
    end if
    if(nblock<1) then
       call add(errmsg,2,"'SCHUNK_INVGRID_NREF_BLOCKS' must be greater than zero",myname)
       goto 1
    end if
    !
    nlay_block => ivecp(inpar,'SCHUNK_INVGRID_NLAY',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCHUNK_INVGRID_NREF_BLOCKS = ",nblock,&
            " integers for 'SCHUNK_INVGRID_NLAY' from '"//trim(inpar.sval.'SCHUNK_INVGRID_NLAY')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(nlay_block < 1)) then
       call add(errmsg,2,"all values of SCHUNK_INVGRID_NLAY must be greater than zero",myname)
       goto 1
    end if
    nlay = sum(nlay_block)    
    !
    thickness => rvecp(inpar,'SCHUNK_INVGRID_THICKNESS',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCHUNK_INVGRID_NREF_BLOCKS = ",nblock,&
            " integers for 'SCHUNK_INVGRID_THICKNESS' from '"//trim(inpar.sval.'SCHUNK_INVGRID_THICKNESS')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(thickness <= 0.)) then
       call add(errmsg,2,"all values of SCHUNK_INVGRID_THICKNESS must be greater than zero",myname)
       goto 1
    end if
    !
    allocate(z(nlay+1))
    z(1) = this%rmax
    iz = 1
    do ibl = 1,nblock
       do ilay = 1,nlay_block(ibl)
          iz = iz + 1
          z(iz) = z(iz-1) - thickness(ibl)
       end do ! ilay
    end do ! ibl
    !
    ! define nlat,nlon
    !
    nlat => ivecp(inpar,'SCHUNK_INVGRID_NLAT',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCHUNK_INVGRID_NREF_BLOCKS = ",nblock,&
            " integers for 'SCHUNK_INVGRID_NLAT' from '"//trim(inpar.sval.'SCHUNK_INVGRID_NLAT')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(nlat < 1)) then
       call add(errmsg,2,"all values of SCHUNK_INVGRID_NLAT must be greater than zero",myname)
       goto 1
    end if
    !
    nlon => ivecp(inpar,'SCHUNK_INVGRID_NLON',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCHUNK_INVGRID_NREF_BLOCKS = ",nblock,&
            " integers for 'SCHUNK_INVGRID_NLON' from '"//trim(inpar.sval.'SCHUNK_INVGRID_NLON')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(nlon < 1)) then
       call add(errmsg,2,"all values of SCHUNK_INVGRID_NLON must be greater than zero",myname)
       goto 1
    end if
    !
    ! define vtk_projection,apply_vtk_coords_scaling_factor,vtk_coords_scaling_factor
    !
    this%vtk_projection = sval(inpar,'VTK_PROJECTION',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"keyword 'VTK_PROJECTION' is not present in parameter file",myname)
       goto 1
    end if
    select case(this%vtk_projection)
    case('GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT')
       ! OK, so do nothing
    case default
       call add(errmsg,2,"value '"//trim(this%vtk_projection)//"' of keyword 'VTK_PROJECTION' is invalid; "//&
            "must be one of 'GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT'",myname)
       goto 1
    end select
    this%apply_vtk_coords_scaling_factor = lval(inpar,'SCALE_VTK_COORDS',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read logical value for 'SCALE_VTK_COORDS' from '"//&
            trim(inpar.sval.'SCALE_VTK_COORDS')//"'",myname)
       goto 1
    end if
    if(this%apply_vtk_coords_scaling_factor) then
       this%vtk_coords_scaling_factor = rval(inpar,'VTK_COORDS_SCALING_FACTOR',iostat=ios)
       if(ios /= 0) then
          call add(errmsg,2,"could not read real value for 'VTK_COORDS_SCALING_FACTOR' from '"//&
               trim(inpar.sval.'VTK_COORDS_SCALING_FACTOR')//"'",myname)
          goto 1
       end if
    end if
    !
    !  create scartInversionGrid object
    !
    wx = 2.0*tan(0.5*wlat*mc_deg2rad)*this%rmax
    wy = 2.0*tan(0.5*wlon*mc_deg2rad)*this%rmax
    call createFromValuesScartInversionGrid(this%scart,nlat,nlon,wx,wy,0.,0.,0.,0.,nblock,nlay_block,z, &
         .false.,.false.,1.0,inpar.sval.'VTK_GEOMETRY_TYPE',errmsg)
    !   call createFromValuesScartInversionGrid(this%scart,nx,ny,wx,wy,cx,cy,zmax,rot,nblock,nlay_block,z,
    !       use_local_coords_for_vtk,apply_vtk_coords_scaling_factor,vtk_coords_scaling_factor)
    !
    if (.level.errmsg == 2) goto 1
    ! if everything was ok, indicate so and return
    this%is_defined = .true.
    if(associated(thickness)) deallocate(thickness)
    if(associated(nlay_block)) deallocate(nlay_block)
    if(associated(z)) deallocate(z)
    if(associated(nlat)) deallocate(nlat)
    if(associated(nlon)) deallocate(nlon)
    return
    !
    ! if there went anything wrong (i.e. there was a goto here), destroy whatever was created so far
1   call deallocateSchunkInversionGrid(this)
    if(associated(thickness)) deallocate(thickness)
    if(associated(nlay_block)) deallocate(nlay_block)
    if(associated(z)) deallocate(z)
    if(associated(nlat)) deallocate(nlat)
    if(associated(nlon)) deallocate(nlon)
    !
  end subroutine createSchunkInversionGrid
!--------------------------------------------------------------------
!> \brief Deallocate createSchunkInversionGrid
!
  subroutine deallocateSchunkInversionGrid(this)
    type (schunk_inversion_grid) :: this
    call deallocateScartInversionGrid(this%scart)
    this%is_defined = .false.
  end subroutine deallocateSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief return overall number of invgrid cells, if any
!
  function getNcellSchunkInversionGrid(this) result(ncell)
    type (schunk_inversion_grid), intent(in) :: this
    integer :: ncell
    if(this%is_defined) then
       ncell = getNcellScartInversionGrid(this%scart)
    else
       ncell = 0
    end if
  end function getNcellSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief transform wp,event or station coords to x,y,z coords for vtk application
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station' and that c1,c2,c3 are associated and have all same length.
!!  Also it can be assured that this schunk inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here).
!!  in case coords_type is 'wp', c1,c2,c3 are expected to be global cartesian x,y,z coordinates for spherical application
!!  in case coords_type is 'event', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and depth (in km)
!!  in case coords_type is 'station', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and altitude (in m)
!!  altitude of stations is ignored on return (in order to plot stations on surface of spherical chunk)
!! \param this schunk inversion grid
!! \param c1 vector or first coordinate (contains vtk x-values on exit)
!! \param c2 vector or second coordinate (contains vtk y-values on exit)
!! \param c3 vector or third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!! \param errmsg error message
!! \param uf_wp optional unit of wavefield points (recommended to be used along with coords_type == 'wp'). If not present, wavefield points are assumed to be in km
!
  subroutine transformToVtkSchunkInversionGrid(this,c1,c2,c3,coords_type,errmsg,uf_wp)
    type (schunk_inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
    real, optional :: uf_wp
    real, dimension(:), allocatable :: lat,lon,z
    real :: uf_correction
!
    call addTrace(errmsg,'transformToVtkSchunkInversionGrid')
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",'transformToVtkSchunkInversionGrid')
       return
    end if
!
    select case(coords_type)
    case('station','event')
       ! when using spherical inversion grids, the event and station coordinates
       ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
       ! so in case of coords_type = 'event' or 'station', transform the incoming c1,c2,c3 to respective
       ! global Cartesian coordinates first, then execute the transformations below
!
       allocate(lat(size(c1)),lon(size(c1)),z(size(c1)))
       lat = c1; lon = c2; z = c3
!
       ! first define unit vectors (c1,c2,c3) in global Cartesian coordinates pointing to the direction of the event/station location
       c1 = cos(lat*mc_deg2rad)*cos(lon*mc_deg2rad)
       c2 = cos(lat*mc_deg2rad)*sin(lon*mc_deg2rad)
       c3 = sin(lat*mc_deg2rad)
       ! then scale the unit vectors (c1,c2,c3) according to the incoming third coordinate z
       select case(coords_type)
       case('station')
          ! for station coordinates, the incoming 3rd coordinate is altitude in meters, which here is ignored: so scale to surface of the earth
          c1 = c1 * this%rmax
          c2 = c2 * this%rmax
          c3 = c3 * this%rmax
       case('event')
          ! for event coordinates, the incoming 3rd coordinate is depth in km
          c1 = c1 * (this%rmax-z)
          c2 = c2 * (this%rmax-z)
          c3 = c3 * (this%rmax-z)
       end select
!
       deallocate(lat,lon,z)
    case('wp')
       ! if the optional argument uf_wp is given, check if wavefield point coordinates are given in km,
       ! otherwise transform to km first, then execute the transformations below
       if(present(uf_wp)) then
          if(uf_wp /= 1.0e3) then
             ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
             uf_correction = uf_wp * 1.0e-3
             c1 = c1 * uf_correction
             c2 = c2 * uf_correction
             c3 = c3 * uf_correction
          end if ! uf_wp /= 1.0e3
       end if ! present(uf_wp)
    end select
!
    ! in case of this%vtk_projection == 'GLOBAL', no further transformations are required (except scaling)
    ! for the other cases, apply further transformations to local chunk or flat scart coordinates, respectively
    select case(this%vtk_projection)
    case('LOCAL_CURV')
       call transformVectorGlobalToLocalSchunkInversionGrid(this,c1,c2,c3,size(c1))
    case('LOCAL_FLAT')
       call transformVectorGlobalToLocalSchunkInversionGrid(this,c1,c2,c3,size(c1))
       call transformVectorLocalToScartSchunkInversionGrid(this,c1,c2,c3,size(c1))
    case('LOCAL_NORTH_CURV')
       call transformVectorGlobalToLocalSchunkInversionGrid(this,c1,c2,c3,size(c1))
       call transformVectorGammaRotationSchunkInversionGrid(this,c1,c2,size(c1))
    case('LOCAL_NORTH_FLAT')
       call transformVectorGlobalToLocalSchunkInversionGrid(this,c1,c2,c3,size(c1))
       call transformVectorLocalToScartSchunkInversionGrid(this,c1,c2,c3,size(c1))
       call transformVectorGammaRotationSchunkInversionGrid(this,c1,c2,size(c1))
    end select
!
    if(this%apply_vtk_coords_scaling_factor) then
       c1 = c1 * this%vtk_coords_scaling_factor
       c2 = c2 * this%vtk_coords_scaling_factor
       c3 = c3 * this%vtk_coords_scaling_factor
    end if
  end subroutine transformToVtkSchunkInversionGrid
!-------------------------------------------------------------------------------------------
!> \brief Assuming local coordinates frame, perform rotation about the Z-axis
!
  subroutine transformVectorGammaRotationSchunkInversionGrid(this,c1,c2,n,inverse)
    type (schunk_inversion_grid) :: this
    real, dimension(n) :: c1,c2
    integer :: n
    logical, optional :: inverse
    real, dimension(n) :: x1,x2
    logical :: inverse_rotation
    if(present(inverse)) then
       inverse_rotation = inverse
    else
       inverse_rotation = .false.
    end if
    x1 = c1; x2 = c2
    if(inverse_rotation) then
       ! apply rotation with transposed matrix tm_gam^T
       c1 = this%tm_gam(1,1)*x1+this%tm_gam(2,1)*x2
       c2 = this%tm_gam(1,2)*x1+this%tm_gam(2,2)*x2
    else
       c1 = this%tm_gam(1,1)*x1+this%tm_gam(1,2)*x2
       c2 = this%tm_gam(2,1)*x1+this%tm_gam(2,2)*x2
    end if
  end subroutine transformVectorGammaRotationSchunkInversionGrid
!-------------------------------------------------------------------------------------------
!> \brief Transform geographic cartesian coordinates to chunk-centered cartesian coordinates
!! \details h(i) = tm(i,n)e(n). Because of xc(i)h(i) = xg(n)e(n) we get xc(i)tm(i,n)e(n)=xg(n)e(n).
!! xg(n) = xc(i)tm(i,n), xg = xc*M = M^T*xc, xc = (M^T)^(-1)*xg = M*xg, xc(i) = tm(i,n)xg(n). 
!
  subroutine transformVectorGlobalToLocalSchunkInversionGrid(this,c1,c2,c3,n)
    type (schunk_inversion_grid) :: this
    real, dimension(n) :: c1,c2,c3
    integer :: n
    real, dimension(n) :: x1,x2,x3
    x1 = c1; x2 = c2; x3 = c3
    c1 = this%tm_global2local(1,1)*x1+this%tm_global2local(1,2)*x2+this%tm_global2local(1,3)*x3
    c2 = this%tm_global2local(2,1)*x1+this%tm_global2local(2,2)*x2+this%tm_global2local(2,3)*x3
    c3 = this%tm_global2local(3,1)*x1+this%tm_global2local(3,2)*x2+this%tm_global2local(3,3)*x3
  end subroutine transformVectorGlobalToLocalSchunkInversionGrid
!-------------------------------------------------------------------------------------------
!> \brief Transform chunk-centered cartesian coordinates to geographic cartesian coordinates
!! \details h(i) = tm(i,n)e(n). Because of xc(i)h(i) = xg(n)e(n) we get xc(i)tm(i,n)e(n)=xg(n)e(n).
!! xg(n) = xc(i)tm(i,n), xg = xc*M = M^T*xc. 
!
  subroutine transformVectorLocalToGlobalSchunkInversionGrid(this,c1,c2,c3,n)
    type (schunk_inversion_grid) :: this
    real, dimension(n) :: c1,c2,c3
    integer :: n
    real, dimension(n) :: x1,x2,x3
    x1 = c1; x2 = c2; x3 = c3
    c1 = this%tm_global2local(1,1)*x1+this%tm_global2local(2,1)*x2+this%tm_global2local(3,1)*x3
    c2 = this%tm_global2local(1,2)*x1+this%tm_global2local(2,2)*x2+this%tm_global2local(3,2)*x3
    c3 = this%tm_global2local(1,3)*x1+this%tm_global2local(2,3)*x2+this%tm_global2local(3,3)*x3
  end subroutine transformVectorLocalToGlobalSchunkInversionGrid
!-------------------------------------------------------------------------------------------
!> \brief Transform chunk-centered cartesian coordinates to geographic cartesian coordinates
!! \details h(i) = tm(i,n)e(n). Because of xc(i)h(i) = xg(n)e(n) we get xc(i)tm(i,n)e(n)=xg(n)e(n).
!! xg(n) = xc(i)tm(i,n), xg = xc*M = M^T*xc. 
!
  subroutine transformPointLocalToGlobalSchunkInversionGrid(this,c1,c2,c3)
    type (schunk_inversion_grid) :: this
    real :: c1,c2,c3
    real :: x1,x2,x3
    x1 = c1; x2 = c2; x3 = c3
    c1 = this%tm_global2local(1,1)*x1+this%tm_global2local(2,1)*x2+this%tm_global2local(3,1)*x3
    c2 = this%tm_global2local(1,2)*x1+this%tm_global2local(2,2)*x2+this%tm_global2local(3,2)*x3
    c3 = this%tm_global2local(1,3)*x1+this%tm_global2local(2,3)*x2+this%tm_global2local(3,3)*x3
  end subroutine transformPointLocalToGlobalSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief Transform local scartInversionGrid coordinates to chunk-centered cartesian schunkInversionGrid coordinates
!! \details Local schunkInversionGrid transforms to local scartInversionGrid by projecting schunk points to
!! tangential plane at r and then expand width of plane to maximum value at rmax 
!
  subroutine transformVectorScartToLocalSchunkInversionGrid(this,c1,c2,c3,n)
    type (schunk_inversion_grid) :: this
    integer :: n
    real, dimension(n) :: c1,c2,c3
    real, dimension(n) :: d
    ! adjust scart x,y-coordinates to local width of tangential plane (*r/rmax)
    c1 = c1*(c3/this%rmax)
    c2 = c2*(c3/this%rmax)
    ! project these points onto to spherical shell of radius c3
    d = c3*sqrt((c1/c3)**2+(c2/c3)**2+1.0)
    c1 = (c1/d)*c3
    c2 = (c2/d)*c3
    c3 = (c3/d)*c3
  end subroutine transformVectorScartToLocalSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief Transform scartInversionGrid coordinates to local chunk-centered cartesian coordinates
!! \details Local schunkInversionGrid transforms to local scartInversionGrid by projecting schunk points to
!! tangential plane at r and then expand width of plane to maximum value at rmax 
!
  subroutine transformPointScartToLocalSchunkInversionGrid(this,c1,c2,c3)
    type (schunk_inversion_grid) :: this
    real :: c1,c2,c3
    real :: d
    ! adjust scart x,y-coordinates to local width of tangential plane (*z/rmax)
    c1 = c1*(c3/this%rmax)
    c2 = c2*(c3/this%rmax)
    ! project these points onto to spherical shell of radius x3
    d = c3*sqrt((c1/c3)**2+(c2/c3)**2+1.0)
    c1 = (c1/d)*c3
    c2 = (c2/d)*c3
    c3 = (c3/d)*c3
  end subroutine transformPointScartToLocalSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief Transform local chunk-centered cartesian coordinates to scartInversionGrid coordinates
!! \details Local local scartInversionGrid transforms to local schunkInversionGrid by projecting scart points to
!! sphere at r
!
  subroutine transformVectorLocalToScartSchunkInversionGrid(this,c1,c2,c3,n)
    type (schunk_inversion_grid) :: this
    integer :: n
    real, dimension(n) :: c1,c2,c3
    real, dimension(n) :: r
    ! project these points onto to tangential plane (xt/r = xg1/xg3)
    r = c3*sqrt((c1/c3)**2+(c2/c3)**2+1.0)
    c1 = (c1/c3)*r
    c2 = (c2/c3)*r
    c3 = r
    ! adjust scart x,y-coordinates to local width of tangential plane (*rmax/z)
    c1 = c1*(this%rmax/c3)
    c2 = c2*(this%rmax/c3)
  end subroutine transformVectorLocalToScartSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief return geometry information on cells for vtk output
!
  subroutine getGeometryVtkSchunkInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,&
       cell_indx_req,indx_map_out)
    type (schunk_inversion_grid) :: this
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    integer :: geometry_type
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
!
   call addTrace(errmsg,'getGeometryVtkSchunkInversionGrid')
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",'getGeometryVtkSchunkInversionGrid')
       return
    endif
!
    call getGeometryVtkScartInversionGrid(this%scart,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,&
       cell_indx_req,indx_map_out)
    if(.level.errmsg == 2) return
!
    ! in case of this%vtk_projection == 'LOCAL_FLAT', no further transformations are required (except scaling)
    ! for the other cases, apply further transformations to local chunk or global coordinates, respectively
    select case(this%vtk_projection)
    case('LOCAL_NORTH_FLAT')
       call transformVectorGammaRotationSchunkInversionGrid(this,points(1,:),points(2,:),size(points,2))
    case('LOCAL_CURV')
       call transformVectorScartToLocalSchunkInversionGrid(this,points(1,:),points(2,:),points(3,:),size(points,2))
    case('LOCAL_NORTH_CURV')
       call transformVectorScartToLocalSchunkInversionGrid(this,points(1,:),points(2,:),points(3,:),size(points,2))
       call transformVectorGammaRotationSchunkInversionGrid(this,points(1,:),points(2,:),size(points,2))
    case('GLOBAL')
       call transformVectorScartToLocalSchunkInversionGrid(this,points(1,:),points(2,:),points(3,:),size(points,2))
       call transformVectorLocalToGlobalSchunkInversionGrid(this,points(1,:),points(2,:),points(3,:),size(points,2))
    end select
!
    if(this%apply_vtk_coords_scaling_factor) then
       points = points*this%vtk_coords_scaling_factor
    end if
  end subroutine getGeometryVtkSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief return indices of all face neighbours for all (optionally only subset of) cells
!!  \details If boundary_conditions is set, there should be a special form of nb_idx returned, which is used
!!  e.g. to define special smoothing conditions:
!!   'no_nb_inner_bnd': with some neighbours removed (which should not be smoothed with on internal invgrid 
!!                      boundaries). FOR schunkInversionGrid GRIDS, THERE IS NO POSSIBILITY YET TO DEFINE INTERNAL
!!                      INVGRID BOUNDARIES!! SO IF THIS IS VALUE IS SET, NO NEIGHBOURS ARE RETURNED! 
!!                      In the future might introduce this functionality via modul scartInversionGrid. It should
!!                      be easy to realize this by introducing certain layers below which there should be an internal boundary.
!!                      (might be even easier to define a set of refinement blocks via inidices, below which there should be 
!!                      an internal boundary)
!!   'extra_nbs_outer_bnd': additional fake neighbours (having cell index 0) in case of zero boundary conditions 
!!                          on outer invgrid boundaries.
!!   'extra_nbs_outer_bnd_except_free_surface': same as 'extra_nbs_outer_bnd', but not applied for cells on free surfaces.
!!   '','standard': no special boundary handling, i.e. standard functionality (exactly the geometrical face neighbours are given)
!! \param this schunk inversion grid
!! \param nb_idx pointer to array of length this%ncell which contains vector pointer to face neighbour indices
!! \param boundary_conditions optional string indicating type of boundary conditions for which neighbours should be returned
!
  subroutine getIndicesFaceNeighboursSchunkInversionGrid(this,nb_idx,boundary_conditions)
    type (schunk_inversion_grid) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
    character(len=*), optional :: boundary_conditions
!
    call getIndicesFaceNeighboursScartInversionGrid(this%scart,nb_idx,boundary_conditions)
  end subroutine getIndicesFaceNeighboursSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief for each cell return indices of wavefield points contained in that cell
!! \param this schunk inversion grid
!! \param x vector or first coordinate of wavefield points
!! \param y vector or second coordinate of wavefield points
!! \param z vector or third coordinate of wavefield points
!! \param uf_wp unit factor of wavefield points
!! \param wp_idx pointer to array of length .ncell.this; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideSchunkInversionGrid(this,x,y,z,uf_wp,wp_idx,errmsg)
    type (schunk_inversion_grid) :: this
    real, dimension(:), intent(in) :: x,y,z
    real :: uf_wp
    type (integer_vector_pointer), dimension(:), pointer :: wp_idx
    type (error_message) :: errmsg
    real, dimension(:), allocatable :: xl,yl,zl
    real :: uf_correction
!
    call addTrace(errmsg,'locateWpInsideSchunkInversionGrid')
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",'locateWpInsideSchunkInversionGrid')
       return
    end if
!
    allocate(xl(size(x)),yl(size(x)),zl(size(x)))
!
    if(uf_wp /= 1.0e3) then
       ! account for wavefield points which are not given in km (e.g. using SPECFEM3D_Cartesian (or NEXD) with SI
       ! units for spherical applications)

       ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
       uf_correction = uf_wp * 1.0e-3
       xl = x * uf_correction
       yl = y * uf_correction
       zl = z * uf_correction
    else
       xl = x
       yl = y
       zl = z
    end if
!
    call transformVectorGlobalToLocalSchunkInversionGrid(this,xl,yl,zl,size(x))
    call transformVectorLocalToScartSchunkInversionGrid(this,xl,yl,zl,size(x))
    call locateWpInsideScartInversionGrid(this%scart,xl,yl,zl,wp_idx,errmsg)
    deallocate(xl,yl,zl)
  end subroutine locateWpInsideSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in cell icell to standard cell and compute their jacobian
!! \param this schunk inversion grid
!! \param icell global inversion grid cell index
!! \param x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param uf_wp unit factor of wavefield points
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param type_standard_cell defines the shape of the standard cell, select specific routine dependent on type (4=Tetrahedron,6=Hexahedron)
!! \param errmsg error message
!
  subroutine transformToStandardCellSchunkInversionGrid(this,icell,x,y,z,uf_wp,jacobian,type_standard_cell,errmsg)
    type (schunk_inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    real :: uf_wp
    type (error_message) :: errmsg
    real :: xmin,xmax,ymin,ymax,zmin,zmax,xm,ym,zm,d,uf_correction
!
    call addTrace(errmsg,'transformToStandardCellSchunkInversionGrid')
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",'transformToStandardCellSchunkInversionGrid')
       return
    end if
    if(type_standard_cell == -1) then
       call add(errmsg,2,"Incoming value of type_standard_cell is -1, indicating a request for total "//&
            "integration weights (e.g. used by integration weights of type 6) instead of jacobian values. "//&
            "This functionality is not supported by inversion grids of type schunkInversionGrid",&
            'transformToStandardCellSchunkInversionGrid')
       return
    end if
!
    ! account for wavefield points which are not given in km (e.g. using SPECFEM3D_Cartesian (or NEXD) with SI
    ! units for spherical applications)
    if(uf_wp /= 1.0e3) then
       ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
       uf_correction = uf_wp * 1.0e-3
       x = x * uf_correction
       y = y * uf_correction
       z = z * uf_correction
    end if
!
    ! transform global wavefield point coordinates to local scart coordinates
    ! then call corresponding scart routine
    call transformVectorGlobalToLocalSchunkInversionGrid(this,x,y,z,size(x))
    call transformVectorLocalToScartSchunkInversionGrid(this,x,y,z,size(x))
    call transformToStandardCellScartInversionGrid(this%scart,icell,x,y,z,jacobian,type_standard_cell,errmsg)
    call getCellBoundariesScartInversiongrid(this%scart,icell,xmin,xmax,ymin,ymax,zmin,zmax)
    zm = 0.5*(zmin+zmax)
    xm = 0.5*(xmin+xmax)
    ym = 0.5*(ymin+ymax)
    d = zm*sqrt(1.+(xm/zm)**2+(ym/zm)**2)
    ! scale scart cell down to its true radius (z/rmax)
    jacobian = jacobian*(zm/this%rmax)**2
    ! correct for strectching due to projection to tangential plane
    ! mapping equations for points on tangential plane with r=z:
    ! x = r tan(theta)cos(phi), y = tan(theta)*sin(phi), z = r
    ! Det J = 1/cos(theta)**3 * r**2*sin(theta), cos(theta) = z/d
    jacobian = jacobian*(zm/d)**3
  end subroutine transformToStandardCellSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of inversion grid cell
!! \param this schunk inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param volume volume of cell icell
!
  subroutine getVolumeCellSchunkInversionGrid(this,icell,volume)
    type (schunk_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax,d
    real :: xm,ym,zm
!
    ! in routine getVolumeCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell is valid
!
    call getCellBoundariesScartInversionGrid(this%scart,icell,xmin,xmax,ymin,ymax,zmin,zmax)
    zm = 0.5*(zmin+zmax)
    xm = 0.5*(xmin+xmax)
    ym = 0.5*(ymin+ymax)
    d = zm*sqrt(1.+(xm/zm)**2+(ym/zm)**2)
!
    volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)*(zm/d)**3
  end subroutine getVolumeCellSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief get center of inversion grid cell
!! \param this schunk inversion grid
!! \param icell index of inversion grid for which center should be returned
!! \param c1 first coordinate of center of cell icell
!! \param c2 second coordinate of center of cell icell
!! \param c3 third coordinate of center of cell icell
!! \param coords_type 'wp','event','station'; optional request
!
  subroutine getCenterCellSchunkInversionGrid(this,icell,c1,c2,c3,coords_type)
    type (schunk_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: c1,c2,c3
    character(len=*), optional :: coords_type
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax
    double precision :: r,theta,phi
!
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell and coords_type are valid
!
    call getCellBoundariesScartInversionGrid(this%scart,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    c1 = 0.5*(xmax+xmin)
    c2 = 0.5*(ymax+ymin)
    c3 = 0.5*(zmax+zmin)
!
    call transformPointScartToLocalSchunkInversionGrid(this,c1,c2,c3)
    call transformPointLocalToGlobalSchunkInversionGrid(this,c1,c2,c3)
!
    ! when using spherical inversion grids, the event and station coordinates
    ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
    ! so in case of coords_type = 'event' or 'station', transform the outgoing c1,c2,c3 
    ! further (which at this point are in global Cartesian coordinates, i.e. 'wp'-form)
    select case(coords_type)
    case('station','event')
       ! r = sqrt(x^2+y^2+z^2)
       ! theta = arccos(z/r) \in [0,pi] (d.h. 90deg-theta rechnen!)
       ! phi = atan2(x,y)
       r = dsqrt( dble(c1)*dble(c1) + dble(c2)*dble(c2) + dble(c3)*dble(c3) )
       theta = acos(dble(c3)/r)
       phi = atan2(dble(c2),dble(c1))
       
       ! define the first two coordinates (c1 = latitude in degrees (-90<=c1<=90), c2 = lon in degrees (0<=lon<=360))
       c1 = real( (0.5d0*mc_pid - theta) / mc_deg2radd )
       c2 = real( phi / mc_deg2radd )

       ! define third coordinate (different definition for 'station' than for 'event'):
       ! 'station': put point onto surface of the earth (ignoring altitude [m], which c3 actually means in that case)
       ! 'event': c3 is depth [km]
       select case(coords_type)
       case('station')
          c3 = this%rmax
       case('event')
          c3 = real( dble(this%rmax)-r )
       end select ! coords_type
    end select ! coords_type
!
  end subroutine getCenterCellSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of inversion grid cell
!! \param this schunk inversion grid
!! \param icell index of inversion grid for which radius should be returned
!! \param radius radius of cell icell
!
  subroutine getRadiusCellSchunkInversionGrid(this,icell,radius)
    type (schunk_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax,d
!
    ! in routine getRadiusCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell is valid
!
    call getCellBoundariesScartInversiongrid(this%scart,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    d = zmax*sqrt(1.+(xmax/zmax)**2+(ymax/zmax)**2)
    radius = zmax/d*0.5*sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin))
  end subroutine getRadiusCellSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief function answers whether a given point inside the inversion grid domain
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station'
!!  Also it can be assured that this schunk inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here).
!!  in case coords_type is 'wp', c1,c2,c3 are expected to be global cartesian x,y,z coordinates for spherical application
!!  in case coords_type is 'event', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and depth (in km)
!!  in case coords_type is 'station', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and altitude (in m)
!!  altitude of stations is ignored on return (in order to plot stations on surface of spherical chunk)
!! \param this schunk inversion grid
!! \param c1 first coordinate
!! \param c2 second coordinate
!! \param c3 third coordinate
!! \param coords_type 'wp','event','station'
!! \param uf_wp unit factor of wavefield points; optional (recommended to be used along with 'wp')
!
  function pointInsideSchunkInversionGrid(this,c1,c2,c3,coords_type,uf_wp) result(l)
    type (schunk_inversion_grid) :: this
    real, intent(in) :: c1,c2,c3
    character(len=*) :: coords_type
    real, optional :: uf_wp
    logical :: l
    real, dimension(1) :: c1_copy,c2_copy,c3_copy
    real :: lat,lon,z,uf_correction
!
    l = .false.
!
    if(.not.this%is_defined) return
!
    c1_copy(1) = c1
    c2_copy(1) = c2
    c3_copy(1) = c3
!
    select case(coords_type)
    case('station','event')
       ! when using spherical inversion grids, the event and station coordinates
       ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
       ! so in case of coords_type = 'event' or 'station', transform the incoming c1,c2,c3 to respective
       ! global Cartesian coordinates first, then execute the transformations below
!
       lat = c1_copy(1); lon = c2_copy(1); z = c3_copy(1)
!
       ! first define unit vector (c1,c2,c3) in global Cartesian coordinates pointing to the direction of the event/station location
       c1_copy(1) = cos(lat*mc_deg2rad)*cos(lon*mc_deg2rad)
       c2_copy(1) = cos(lat*mc_deg2rad)*sin(lon*mc_deg2rad)
       c3_copy(1) = sin(lat*mc_deg2rad)
       ! then scale the unit vector (c1,c2,c3) according to the incoming third coordinate z
       select case(coords_type)
       case('station')
          ! for station coordinates, the incoming 3rd coordinate is altitude in meters, which here is ignored: so scale to surface of the earth
          c1_copy(1) = c1_copy(1) * this%rmax
          c2_copy(1) = c2_copy(1) * this%rmax
          c3_copy(1) = c3_copy(1) * this%rmax
       case('event')
          ! for event coordinates, the incoming 3rd coordinate is depth in km
          c1_copy(1) = c1_copy(1) * (this%rmax-z)
          c2_copy(1) = c2_copy(1) * (this%rmax-z)
          c3_copy(1) = c3_copy(1) * (this%rmax-z)
       end select
    case('wp')
       ! if the optional argument uf_wp is given, check if wavefield point coordinates are given in km,
       ! otherwise transform to km first, then execute the transformations below
       if(present(uf_wp)) then
          if(uf_wp /= 1.0e3) then
             ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
             uf_correction = uf_wp * 1.0e-3
             c1_copy(1) = c1_copy(1) * uf_correction
             c2_copy(1) = c2_copy(1) * uf_correction
             c3_copy(1) = c3_copy(1) * uf_correction
          end if ! uf_wp /= 1.0e3
       end if ! present(uf_wp)
    end select
!
    call transformVectorGlobalToLocalSchunkInversionGrid(this,c1_copy,c2_copy,c3_copy,1)
    call transformVectorLocalToScartSchunkInversionGrid(this,c1_copy,c2_copy,c3_copy,1)
    l = pointInsideScartInversionGrid(this%scart,c1_copy(1),c2_copy(1),c3_copy(1),coords_type)
  end function pointInsideSchunkInversionGrid
!
end module schunkInversionGrid
