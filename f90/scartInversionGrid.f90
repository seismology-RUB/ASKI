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
!> \brief simple layered Cartesian inversion grid
!!
!! \details The inversion domain that can be covered with a simple layered Cartesian inversion grid
!!   is a cuboid of certain width which may be rotated about the local central z-axis about some angle and then 
!!   translated by a certain vector. From outside the module, always absolute coordinate values
!!   are assumed (e.g. from module wavefieldPoints), except for the parameter file which creates
!!   a scart inversion grid. Internally, however, the cell boundaries are stored in terms of a centralized
!!   coordinate system (cuboid is centered ad x=y=zmay=0 and no rotation is applied)
!!
!! \author Florian Schumacher
!! \date Okt 2015
!
module scartInversionGrid
!
  use inputParameter
  use vectorPointer
  use mathConstants
  use errorMessage
!
  implicit none
!
  private :: defineFaceNeighboursScartInversionGrid,getBoundariesCoveringCellsScartInversionGrid,&
       transformVectorGlobalToLocalScartInversionGrid,transformPointLocalToGlobalScartInversionGrid,&
       locatePointInCellScartInversionGrid,locateCoordinateScartInversionGrid,&
       inverseCellIndexScartInversionGrid,intGeometryTypeScartInversionGrid !,getCellBoundariesScartInversiongrid
!
  type scart_inversion_grid
     private
     logical :: is_defined = .false. !< flag indicating the correct definition (initialization) of the object (i.e. all following values)
!
     ! GEOMETRY OF INVERSION GRID DOMAIN
     real :: cx,cy,zmax !< translation vector
     real :: rotation_angle !< angle in degrees of counter-clockwise rotation about local z-axis (before translating)
     real, dimension(2,2) :: Mrot !< matrix of rotation in XY plane from real coordinates back to reference frame (i.e. CLOCKWISE rotation)
!
     ! INVERSION GRID CELLS, coordinates x,y,z given in LOCAL reference frame (cuboid is centered at x=y=zmax=0, no rotation)
!
     ! REFINEMENT BLOCKS
     integer :: nblock !< number of refinement blocks with constant lateral resolution
     integer, dimension(:), pointer :: nlay_block => null() !< for each block, contains number of layers
     integer, dimension(:), pointer :: iblock => null() !< for each layer, contains index of block (inverse mapping to nlay_block)
     integer :: nlay !< total number of layers
!
     ! CELLS
     integer :: ncell !< overall number of inversion grid cells
     integer, dimension(:), pointer :: ncell_layer_cum => null() !< for each layer, cumulative number of cells in all above layers
     ! layer definition
     real, dimension(:), pointer :: z => null() !< z values of layer boundaries in decreasing order (size = nlay+1)
     ! refinement block specific lateral definitions of cells in the layers
     integer, dimension(:), pointer :: nx => null(), ny => null() !< for each block, contains number of invgrid cells in x (y) direction
     type (real_vector_pointer), dimension(:), pointer :: x => null(), y => null() !< for each block, contains bounding coordinates of cells in increasing order
!
     ! NEIGHBOURS
     integer, dimension(:,:), pointer :: face_neighbours_xy => null() !< dim(4,ncell), xmin,xmax,ymin,ymax neighb. indx for each cell (set -1 for no neighb)
     type (integer_vector_pointer), dimension(:,:), pointer :: face_neighbours_z => null() !< dim(2,ncell), zmin,zmax neighb. indices each cell
!
     ! SPECIFICATION FOR VTK OUTPUT
     logical :: use_local_coords_for_vtk
     logical :: apply_vtk_coords_scaling_factor
     real :: vtk_coords_scaling_factor
     integer :: vtk_geometry_type_int = -1 !< type of vtk geometry:  0 = volumetric cells , 1 = cell center points
     ! in the future: there could be flags in parameter file like: DONT_SMOOTH_LAYER_BOUNDARIES, 
     ! or SMOOTHING_BOUNDARY_CONDITIONS which could be taken into account here, and memorized for better handling 
     ! of smoothing conditions in calls to certain routines below
  end type scart_inversion_grid
!
contains
!------------------------------------------------------------------------
!> \brief logical return whether this scart_inversion_grid is able to transform points (of given coords type) to vtk plot projection
!
  function canTransformToVtkPointsOutsideScartInversionGrid(this,coords_type) result(l)
    type(scart_inversion_grid) :: this
    character(len=*) :: coords_type
    logical :: l
    ! the scart_inversion_grid has capability to transform any points to vtk,
    ! provided the object is defined
    l = this%is_defined
  end function canTransformToVtkPointsOutsideScartInversionGrid
!------------------------------------------------------------------------
!> \brief get unit factor of the volume element
!! \param this scart inversion grid
!! \param uf_wp unit factor of wavefield points
!! \param uf_vol unit factor of volume element (return value of this subroutine)
!! \param errmsg error message
!
  subroutine getUnitFactorOfVolumeElementScartInversionGrid(this,uf_wp,uf_vol,errmsg)
    type (scart_inversion_grid) :: this
    real :: uf_wp,uf_vol
    type (error_message) :: errmsg
    character (len=46) :: myname = 'getUnitFactorOfVolumeElementScartInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) call add(errmsg,1,"be aware that the inversion grid not yet defined; "//&
         "however, the unit factor of the volume element can be correctly computed at this point",myname)
!
    ! The scart inversion grid does not assume specific units for its spatial extension but simply 
    ! assumes that they are the same as for the wavefield point coordinates (which are located inside the
    ! inversion grid cells only according to their pure numerical values). Hence, for this 3D volumetric
    ! inversion grid, the volume element has a unit which is the cube of the unit of the wavefield points.
    uf_vol = uf_wp*uf_wp*uf_wp
  end subroutine getUnitFactorOfVolumeElementScartInversionGrid
!------------------------------------------------------------------------
!> \brief map vtk geometry type names to integers
!
  function intGeometryTypeScartInversionGrid(vtk_geometry_type_str) result(vtk_geometry_type_int)
    character(len=*), intent(in) :: vtk_geometry_type_str
    integer :: vtk_geometry_type_int
    select case(trim(vtk_geometry_type_str))
    case('CELLS')
       vtk_geometry_type_int = 0
    case('CELL_CENTERS')
       vtk_geometry_type_int = 1
    case default
       vtk_geometry_type_int = -1
    end select
  end function intGeometryTypeScartInversionGrid
!------------------------------------------------------------------------
!> \brief create simple Cartesian inversion grid
!! \details construct by passing all necessary values
!
  subroutine createFromValuesScartInversionGrid(this,nx,ny,wx,wy,cx,cy,zmax,rot,nblock,nlay_block,z,&
    & use_local_coords_for_vtk,apply_vtk_coords_scaling_factor,vtk_coords_scaling_factor,vtk_geometry_type,errmsg)
    type (scart_inversion_grid) :: this
    integer :: nblock
    integer, dimension(:) :: nlay_block,nx,ny
    real :: wx,wy,cx,cy,zmax,rot
    real, dimension(:) :: z
    logical :: use_local_coords_for_vtk
    logical :: apply_vtk_coords_scaling_factor
    real :: vtk_coords_scaling_factor
    character(len=*) :: vtk_geometry_type
    type (error_message) :: errmsg
    integer :: ibl,ilay,ix,iy
    real, dimension(:), pointer :: rp
    character (len=400) :: errstr
    character (len=34) :: myname = 'createFromValuesScartInversionGrid'
!
    nullify(rp)
    call addTrace(errmsg,myname)
    if(this%is_defined) then
       call add(errmsg,1,"this object is already defined, deallocating it now before creating new one",myname)
       call deallocateScartInversionGrid(this)
    end if
    this%cx = cx
    this%cy = cy
    this%zmax = zmax
    if(wx <= 0) then
       write(errstr,*) "SCART_INVGRID_WX = ",wx," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(wy <= 0) then
       write(errstr,*) "SCART_INVGRID_WY = ",wy," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    this%rotation_angle = rot
    ! define Mrot from rotation_angle, which rotates CLOCKWISEly from real XY to reference coordinates 
    ! (i.e. inverse rotation from rotated block back to non-rotated)
    this%Mrot(1,1) = cos(this%rotation_angle*mc_deg2rad)
    this%Mrot(1,2) = sin(this%rotation_angle*mc_deg2rad)
    this%Mrot(2,1) = -sin(this%rotation_angle*mc_deg2rad)
    this%Mrot(2,2) = cos(this%rotation_angle*mc_deg2rad)
!
    this%nblock = nblock
    if(this%nblock<1) then
       call add(errmsg,2,"'SCART_INVGRID_NREF_BLOCKS' must be greater than zero",myname)
       goto 1
    end if
!
    if (size(nlay_block) /= nblock) then
       call add(errmsg,2,"size of nlay_block /= nblock",myname)
       goto 1
    endif
    allocate(this%nlay_block(nblock))
    this%nlay_block = nlay_block
    if(any(this%nlay_block < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NLAY must be greater than zero",myname)
       goto 1
    end if
    this%nlay = sum(this%nlay_block)
!
    allocate(this%iblock(this%nlay))
    ilay = 0
    do ibl = 1,this%nblock
       this%iblock((ilay+1):(ilay+this%nlay_block(ibl))) = ibl
       ilay = ilay + this%nlay_block(ibl)
    end do ! ibl
!
    if (size(z) /= this%nlay+1) then
       call add(errmsg,2,"size of z inconsistent with nlay",myname)
       goto 1
    endif
    allocate(this%z(this%nlay+1))
    this%z = z
 !
    if (size(nx) /= nblock) then
       call add(errmsg,2,"size of nx /= nblock",myname)
       goto 1
    endif
    allocate(this%nx(nblock))
    this%nx = nx
    if(any(this%nx < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NX must be greater than zero",myname)
       goto 1
    end if
    allocate(this%x(this%nblock))
    do ibl = 1,this%nblock
       allocate(rp(this%nx(ibl)+1))
       rp(1) = -0.5*wx
       do ix = 2,this%nx(ibl)
          rp(ix) = rp(1) + ((ix-1)*wx)/this%nx(ibl)
       end do ! ix
       rp(this%nx(ibl)+1) = 0.5*wx
       call associateVectorPointer(this%x(ibl),rp)
       nullify(rp)
    end do ! ibl
!
    if (size(ny) /= nblock) then
       call add(errmsg,2,"size of ny /= nblock",myname)
       goto 1
    endif
    allocate(this%ny(nblock))
    this%ny = ny
    if(any(this%ny < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NY must be greater than zero",myname)
       goto 1
    end if
    allocate(this%y(this%nblock))
    do ibl = 1,this%nblock
       allocate(rp(this%ny(ibl)+1))
       rp(1) = -0.5*wy
       do iy = 2,this%ny(ibl)
          rp(iy) = rp(1) + ((iy-1)*wy)/this%ny(ibl)
       end do ! iy
       rp(this%ny(ibl)+1) = 0.5*wy
       call associateVectorPointer(this%y(ibl),rp)
       nullify(rp)
    end do ! ibl
!
    ! define number of inversion grid cells
    allocate(this%ncell_layer_cum(this%nlay))
    this%ncell_layer_cum(1) = 0
    do ilay = 2,this%nlay
       this%ncell_layer_cum(ilay) = this%ncell_layer_cum(ilay-1) + &
            this%nx(this%iblock(ilay-1))*this%ny(this%iblock(ilay-1))
    end do ! ilay
    this%ncell = this%ncell_layer_cum(this%nlay) + &
            this%nx(this%iblock(this%nlay))*this%ny(this%iblock(this%nlay))
!
    this%use_local_coords_for_vtk = use_local_coords_for_vtk
    this%apply_vtk_coords_scaling_factor = apply_vtk_coords_scaling_factor
    this%vtk_coords_scaling_factor = vtk_coords_scaling_factor
    this%vtk_geometry_type_int = intGeometryTypeScartInversionGrid(vtk_geometry_type)
    if(this%vtk_geometry_type_int < 0) then
       call add(errmsg,2,"incoming vtk geometry type '"//trim(vtk_geometry_type)//"' is invalid",myname)
       goto 1
    end if
!
    ! define face neigbours for all cells
    ! in the future: there could be flags in parameter file like: DONT_SMOOTH_LAYER_BOUNDARIES, 
    ! or SMOOTHING_BOUNDARY_CONDITIONS which could be taken into account here, and memorized for better handling of smoothing conditions
    call defineFaceNeighboursScartInversionGrid(this)
!
    ! if everything was ok, indicate so and return
    this%is_defined = .true.
    return
!
    ! if there went anything wrong (i.e. there was a goto here), destroy whatever was created so far
1   call deallocateScartInversionGrid(this)
  end subroutine createFromValuesScartInversionGrid
!------------------------------------------------------------------------
!> \brief create simple Cartesian inversion grid
!! \details the parameter file given, must contain all necessary parameters to define
!!  an object of this type. 
!! \param this simple Cartesian inversion grid
!! \param parfile filename of parameter file containing definintion of this inversion grid
!! \param lu file unit to use for reading and writing files
!! \param errmsg error message
!
  subroutine createScartInversionGrid(this,parfile,lu,errmsg)
    type (scart_inversion_grid) :: this
    character(len=*) :: parfile
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=24) :: myname = 'createScartInversionGrid'
    integer :: ios,ilay,iblock,iz,ix,iy
    real :: wx,wy
    real, dimension(:), pointer :: thickness,rp
    ! parfile
    type (input_parameter) :: inpar
    character (len=80), dimension(15) :: inpar_keys
    data inpar_keys/'SCART_INVGRID_CX','SCART_INVGRID_CY','SCART_INVGRID_ZMAX','SCART_INVGRID_WX',&
         'SCART_INVGRID_WY','SCART_INVGRID_ROT','SCART_INVGRID_NREF_BLOCKS','SCART_INVGRID_NLAY',&
         'SCART_INVGRID_THICKNESS','SCART_INVGRID_NX','SCART_INVGRID_NY','USE_LOCAL_INVGRID_COORDS_FOR_VTK',&
         'SCALE_VTK_COORDS','VTK_COORDS_SCALING_FACTOR','VTK_GEOMETRY_TYPE'/
!
    nullify(thickness,rp)
    call addTrace(errmsg,myname)
    if(this%is_defined) then
       call add(errmsg,1,"this object is already defined, deallocating it now before creating new one",myname)
       call deallocateScartInversionGrid(this)
    end if
!
    call createKeywordsInputParameter(inpar,inpar_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) return
!
    ! define cx,cy,zmax,wx,wy,rotation_angle
!
    this%cx = rval(inpar,'SCART_INVGRID_CX',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_CX' from '"//&
            trim(inpar.sval.'SCART_INVGRID_CX')//"'",myname)
       goto 1
    end if
!
    this%cy = rval(inpar,'SCART_INVGRID_CY',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_CY' from '"//&
            trim(inpar.sval.'SCART_INVGRID_CY')//"'",myname)
       goto 1
    end if
!
    this%zmax = rval(inpar,'SCART_INVGRID_ZMAX',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_ZMAX' from '"//&
            trim(inpar.sval.'SCART_INVGRID_ZMAX')//"'",myname)
       goto 1
    end if
!
    wx = rval(inpar,'SCART_INVGRID_WX',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_WX' from '"//&
            trim(inpar.sval.'SCART_INVGRID_WX')//"'",myname)
       goto 1
    end if
    if(wx <= 0) then
       write(errstr,*) "SCART_INVGRID_WX = ",wx," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    wy = rval(inpar,'SCART_INVGRID_WY',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_WY' from '"//&
            trim(inpar.sval.'SCART_INVGRID_WY')//"'",myname)
       goto 1
    end if
    if(wy <= 0) then
       write(errstr,*) "SCART_INVGRID_WY = ",wy," must be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    this%rotation_angle = rval(inpar,'SCART_INVGRID_ROT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'SCART_INVGRID_ROT' from '"//&
            trim(inpar.sval.'SCART_INVGRID_ROT')//"'",myname)
       goto 1
    end if
!
    ! define Mrot from rotation_angle, which rotates CLOCKWISEly from real XY to reference coordinates 
    ! (i.e. inverse rotation from rotated block back to non-rotated)
    this%Mrot(1,1) = cos(this%rotation_angle*mc_deg2rad)
    this%Mrot(1,2) = sin(this%rotation_angle*mc_deg2rad)
    this%Mrot(2,1) = -sin(this%rotation_angle*mc_deg2rad)
    this%Mrot(2,2) = cos(this%rotation_angle*mc_deg2rad)
!
    ! define nblock,nlay_block,iblock,nlay,z
!
    this%nblock = ival(inpar,'SCART_INVGRID_NREF_BLOCKS',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read integer value for 'SCART_INVGRID_NREF_BLOCKS' from '"//&
            trim(inpar.sval.'SCART_INVGRID_NREF_BLOCKS')//"'",myname)
       goto 1
    end if
    if(this%nblock<1) then
       call add(errmsg,2,"'SCART_INVGRID_NREF_BLOCKS' must be greater than zero",myname)
       goto 1
    end if
!
    this%nlay_block => ivecp(inpar,'SCART_INVGRID_NLAY',this%nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCART_INVGRID_NREF_BLOCKS = ",this%nblock,&
            " integers for 'SCART_INVGRID_NLAY' from '"//trim(inpar.sval.'SCART_INVGRID_NLAY')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(this%nlay_block < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NLAY must be greater than zero",myname)
       goto 1
    end if
    this%nlay = sum(this%nlay_block)
!
    allocate(this%iblock(this%nlay))
    ilay = 0
    do iblock = 1,this%nblock
       this%iblock((ilay+1):(ilay+this%nlay_block(iblock))) = iblock
       ilay = ilay + this%nlay_block(iblock)
    end do ! iblock
!
    thickness => rvecp(inpar,'SCART_INVGRID_THICKNESS',this%nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCART_INVGRID_NREF_BLOCKS = ",this%nblock,&
            " integers for 'SCART_INVGRID_THICKNESS' from '"//trim(inpar.sval.'SCART_INVGRID_THICKNESS')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(thickness <= 0.)) then
       call add(errmsg,2,"all values of SCART_INVGRID_THICKNESS must be greater than zero",myname)
       goto 1
    end if
!
    allocate(this%z(this%nlay+1))
    this%z(1) = 0.
    iz = 1
    do iblock = 1,this%nblock
       do ilay = 1,this%nlay_block(iblock)
          iz = iz + 1
          this%z(iz) = this%z(iz-1) - thickness(iblock)
       end do ! ilay
    end do ! iblock
    deallocate(thickness)
!
    ! define nx,ny,x,y
!
    this%nx => ivecp(inpar,'SCART_INVGRID_NX',this%nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCART_INVGRID_NREF_BLOCKS = ",this%nblock,&
            " integers for 'SCART_INVGRID_NX' from '"//trim(inpar.sval.'SCART_INVGRID_NX')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(this%nx < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NX must be greater than zero",myname)
       goto 1
    end if
!
    allocate(this%x(this%nblock))
    do iblock = 1,this%nblock
       allocate(rp(this%nx(iblock)+1))
       rp(1) = -0.5*wx
       do ix = 2,this%nx(iblock)
          rp(ix) = rp(1) + ((ix-1)*wx)/this%nx(iblock)
       end do ! ix
       rp(this%nx(iblock)+1) = 0.5*wx
       call associateVectorPointer(this%x(iblock),rp)
       nullify(rp)
    end do ! iblock
!
    this%ny => ivecp(inpar,'SCART_INVGRID_NY',this%nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read SCART_INVGRID_NREF_BLOCKS = ",this%nblock,&
            " integers for 'SCART_INVGRID_NY' from '"//trim(inpar.sval.'SCART_INVGRID_NY')//"'"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(this%ny < 1)) then
       call add(errmsg,2,"all values of SCART_INVGRID_NY must be greater than zero",myname)
       goto 1
    end if
!
    allocate(this%y(this%nblock))
    do iblock = 1,this%nblock
       allocate(rp(this%ny(iblock)+1))
       rp(1) = -0.5*wy
       do iy = 2,this%ny(iblock)
          rp(iy) = rp(1) + ((iy-1)*wy)/this%ny(iblock)
       end do ! iy
       rp(this%ny(iblock)+1) = 0.5*wy
       call associateVectorPointer(this%y(iblock),rp)
       nullify(rp)
    end do ! iblock
!
    ! define number of inversion grid cells
    allocate(this%ncell_layer_cum(this%nlay))
    this%ncell_layer_cum(1) = 0
    do ilay = 2,this%nlay
       this%ncell_layer_cum(ilay) = this%ncell_layer_cum(ilay-1) + &
            this%nx(this%iblock(ilay-1))*this%ny(this%iblock(ilay-1))
    end do ! iblock
    this%ncell = this%ncell_layer_cum(this%nlay) + &
            this%nx(this%iblock(this%nlay))*this%ny(this%iblock(this%nlay))
!
    ! define use_local_coords_for_vtk,apply_vtk_coords_scaling_factor,vtk_coords_scaling_factor
!
    this%use_local_coords_for_vtk = lval(inpar,'USE_LOCAL_INVGRID_COORDS_FOR_VTK',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read logical value for 'USE_LOCAL_INVGRID_COORDS_FOR_VTK' from '"//&
            trim(inpar.sval.'USE_LOCAL_INVGRID_COORDS_FOR_VTK')//"'",myname)
       goto 1
    end if
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
    this%vtk_geometry_type_int = intGeometryTypeScartInversionGrid(inpar.sval.'VTK_GEOMETRY_TYPE')
    if(this%vtk_geometry_type_int < 0) then
       call add(errmsg,2,"vtk geometry type '"//trim(inpar.sval.'VTK_GEOMETRY_TYPE')//"' is invalid",myname)
       goto 1
    end if
!
    ! define face neigbours for all cells
    ! in the future: there could be flags in parameter file like: DONT_SMOOTH_LAYER_BOUNDARIES, 
    ! or SMOOTHING_BOUNDARY_CONDITIONS which could be taken into account here, and memorized for better handling of smoothing conditions
    call defineFaceNeighboursScartInversionGrid(this)
!
    ! if everything was ok, indicate so and return
    this%is_defined = .true.
    return
!
    ! if there went anything wrong (i.e. there was a goto here), destroy whatever was created so far
1   call deallocateScartInversionGrid(this)
    if(associated(thickness)) deallocate(thickness)
!
  end subroutine createScartInversionGrid
!------------------------------------------------------------------------
!> \brief define neighbours this%face_neighbours of the inversion grid cells
!
  subroutine defineFaceNeighboursScartInversionGrid(this)
    type (scart_inversion_grid) :: this
    ! local
    integer :: iblock,iz,iy,ix,icell,nz
    integer :: nx_below,nx,nx_above
    integer :: ny_below,ny,ny_above
    integer :: j,jx,jx1,jx2,jy,jy1,jy2
    integer :: ncell_cum_below,ncell_cum_above
    integer, dimension(:), pointer :: ip
    real :: xmin,xmax,ymin,ymax
    real, dimension(:), pointer :: x_below,x,x_above,y_below,y,y_above
    logical :: layer_is_bottom,layer_is_top
!
    nullify(ip)
!
    allocate(this%face_neighbours_xy(4,this%ncell),this%face_neighbours_z(2,this%ncell))
!
    ! iterate over blocks, as for "inner" layers of a block (i.e. everything between top and bottom layer of a block), 
    ! the z-neighbours can be easily defined
!
    icell = 0
!
    do iblock = 1,this%nblock
!
       if(iblock<this%nblock) then
          x_below => getVectorPointer(this%x(iblock+1)); nx_below = this%nx(iblock+1)
          y_below => getVectorPointer(this%y(iblock+1)); ny_below = this%ny(iblock+1)
          ncell_cum_below = this%ncell_layer_cum(sum(this%nlay_block(1:iblock))+1)
       end if
       x => getVectorPointer(this%x(iblock)); nx = this%nx(iblock)
       y => getVectorPointer(this%y(iblock)); ny = this%ny(iblock)
       if(iblock>1) then
          x_above => getVectorPointer(this%x(iblock-1)); nx_above = this%nx(iblock-1)
          y_above => getVectorPointer(this%y(iblock-1)); ny_above = this%ny(iblock-1)
          ! store number of all cells above, except for the above layer, in order for the following line of code to work below
          !    ip(j) = ncell_cum_above + (jy-1)*nx_above + jx
          ncell_cum_above = this%ncell_layer_cum(sum(this%nlay_block(1:iblock-1)))
       end if
       nz = this%nlay_block(iblock)
!
       do iz = 1,nz
          layer_is_top = iz==1
          layer_is_bottom = iz==nz
!
          do iy = 1,ny
             do ix = 1,nx
                icell = icell + 1
!
                ! xmin-neighbour
                if(ix==1) then
                   this%face_neighbours_xy(1,icell) = -1
                else
                   this%face_neighbours_xy(1,icell) = icell-1
                end if
                ! xmax-neighbour
                if(ix==nx) then
                   this%face_neighbours_xy(2,icell) = -1
                else
                   this%face_neighbours_xy(2,icell) = icell+1
                end if
!
                ! ymin-neighbour
                if(iy==1) then
                   this%face_neighbours_xy(3,icell) = -1
                else
                   this%face_neighbours_xy(3,icell) = icell-nx
                end if
                ! ymax-neighbour
                if(iy==ny) then
                   this%face_neighbours_xy(4,icell) = -1
                else
                   this%face_neighbours_xy(4,icell) = icell+nx
                end if
!
                ! zmin-neighbour(s)
                if(layer_is_bottom) then
                   ! search in top layer of block below (if any) for zmin neighbours
                   if(iblock/=this%nblock) then
                      ! get xmin,xmax,ymin,ymax of current cell
                      xmin = x(ix); xmax = x(ix+1)
                      ymin = y(iy); ymax = y(iy+1)
!
                      ! find covering cells in x-direction in top layer of block below
                      call getBoundariesCoveringCellsScartInversionGrid(this,xmin,xmax,x_below,nx_below+1,jx1,jx2)
                      ! find covering cells in y-direction in top layer of block below
                      call getBoundariesCoveringCellsScartInversionGrid(this,ymin,ymax,y_below,ny_below+1,jy1,jy2)
!
                      allocate(ip((jx2-jx1)*(jy2-jy1)))
                      j = 0
                      do jy = jy1,jy2-1
                         do jx = jx1,jx2-1
                            j = j+1
                            ip(j) = ncell_cum_below + (jy-1)*nx_below + jx
                         end do ! jx
                      end do ! jy
                      call associateVectorPointer(this%face_neighbours_z(1,icell),ip)
                      nullify(ip)
                   ! else, if this is the bottom layer of the bottom block, there are no zmin neighbours, 
                   ! so do nothing (as vector pointer is nullified already)
                   end if
                else
                   ! if this is not the bottom layer of the block, there is a layer below with same x,y-resolution
                   ! so there is only one zmin neighbour, directly below
                   allocate(ip(1)); ip(1) = icell + nx*ny
                   call associateVectorPointer(this%face_neighbours_z(1,icell),ip)
                   nullify(ip)
                end if ! layer_is_bottom
!
                ! zmax-neighbour(s)
                if(layer_is_top) then
                   ! search in bottom layer of above block (if any) for zmax neighbours
                   if(iblock/=1) then
                      ! get xmin,xmax,ymin,ymax of current cell
                      xmin = x(ix); xmax = x(ix+1)
                      ymin = y(iy); ymax = y(iy+1)
!
                      ! find covering cells in x-direction in top layer of above block
                      call getBoundariesCoveringCellsScartInversionGrid(this,xmin,xmax,x_above,nx_above+1,jx1,jx2)
                      ! find covering cells in y-direction in top layer of above block
                      call getBoundariesCoveringCellsScartInversionGrid(this,ymin,ymax,y_above,ny_above+1,jy1,jy2)
!
                      allocate(ip((jx2-jx1)*(jy2-jy1)))
                      j = 0
                      do jy = jy1,jy2-1
                         do jx = jx1,jx2-1
                            j = j+1
                            ip(j) = ncell_cum_above + (jy-1)*nx_above + jx
                         end do ! jx
                      end do ! jy
                      call associateVectorPointer(this%face_neighbours_z(2,icell),ip)
                      nullify(ip)
                   end if
                   ! else, if this is the top layer of the first block, there are no zmax neighbours, 
                   ! so do nothing (as vector pointer is nullified already)
                else
                   ! if this is not the top layer of the block, there is a layer above with same x,y-resolution
                   ! so there is only one zmax neighbour, directly above
                   allocate(ip(1)); ip(1) = icell - nx*ny
                   call associateVectorPointer(this%face_neighbours_z(2,icell),ip)
                   nullify(ip)
                end if ! layer_is_top
!
             end do ! ix
          end do ! iy
       end do ! iz
    end do ! iblock
  end subroutine defineFaceNeighboursScartInversionGrid
!------------------------------------------------------------------------
!> \brief for a coordinate interval, select from given boundaries those which cover the interval
!! \details given a vector of ncoords coordinates coords (assuming increasing order), find
!!  indices ncov1<ncov2 such that [coords(ncov1),coords(ncov2)] is the smalles interval with boundaries
!!  in coords which completely overlaps with [min,max]
!! \param this scart inversion grid
!! \param min lower bound of interval [min,max] which should be covered
!! \param max upper bound of interval [min,max] which should be covered
!! \param coords array of boundaries to select from for covering interval
!! \param ncoords size of coords
!! \param[out] ncov1 largest lower bound of an interval covering [min,max]
!! \param[out] ncov2 smallest upper bound of an interval covering [min,max]
!
  subroutine getBoundariesCoveringCellsScartInversionGrid(this,min,max,coords,ncoords,ncov1,ncov2)
    type (scart_inversion_grid) :: this
    real :: min,max
    integer :: ncoords,ncov1,ncov2
    real, dimension(ncoords) :: coords
    ! local
    real :: w
!
    w = coords(ncoords) - coords(1)
!
    ! first guess of ncov1
    ncov1 = floor( (min - (-0.5*w)) / (w/(ncoords-1)) ) + 1
!
    ! if necessary, correct the first guess of ncov1: as above computation operations differ from 
    ! operations for constructing coords, the result may be wrong if min is very close to a coords entry
    if(coords(ncov1) > min) then
       ! there will be a segmentation fault, if ncov1 gets negative, but in that case, coords and min violate
       ! conventions (and this module is intrinsically erroneous), so leave that risk here in order to stop
       ! the program and make aware of that problem
       do while(coords(ncov1) > min)
          ncov1 = ncov1 - 1
       end do
    elseif(coords(ncov1+1) <= min) then
       do while(coords(ncov1+1) <= min)
          ncov1 = ncov1 + 1
       end do
    end if
!
    ! first guess of ncov2
    ncov2 = ceiling( (max - (-0.5*w)) / (w/(ncoords-1)) ) + 1
    ! if necessary, correct the first guess of ncov2: as above computation operations differ from 
    ! operations for constructing coords, the result may be wrong if min is very close to a coords entry
    if(coords(ncov2) < max) then
       ! there will be a segmentation fault, if ncov1 gets negative, but in that case, coords and min violate
       ! conventions (and this module is intrinsically erroneous), so leave that risk here in order to stop
       ! the program and make aware of that problem
       do while(coords(ncov2) < max)
          ncov2 = ncov2 + 1
       end do
    elseif(coords(ncov2-1) >= max) then
       do while(coords(ncov2-1) >= max)
          ncov2 = ncov2 - 1
       end do
    end if
  end subroutine getBoundariesCoveringCellsScartInversionGrid
!------------------------------------------------------------------------
!> \brief deallocate simple layered Cartesian inversion grid
!! \param this simple layered Cartesian inversion grid
!
  subroutine deallocateScartInversionGrid(this)
    type (scart_inversion_grid) :: this
    integer :: n
    if(associated(this%nlay_block)) deallocate(this%nlay_block)
    if(associated(this%iblock)) deallocate(this%iblock)
    if(associated(this%ncell_layer_cum)) deallocate(this%ncell_layer_cum)
    if(associated(this%z)) deallocate(this%z)
    if(associated(this%nx)) deallocate(this%nx)
    if(associated(this%ny)) deallocate(this%ny)
    if(associated(this%x)) then 
       do n = 1,size(this%x)
          call dealloc(this%x(n))
       end do ! n
       deallocate(this%x)
    end if
    if(associated(this%y)) then 
       do n = 1,size(this%y)
          call dealloc(this%y(n))
       end do ! n
       deallocate(this%y)
    end if
    if(associated(this%face_neighbours_xy)) deallocate(this%face_neighbours_xy)
    if(associated(this%face_neighbours_z)) then
       do n = 1,size(this%face_neighbours_z,2)
          call dealloc(this%face_neighbours_z(1,n))
          call dealloc(this%face_neighbours_z(2,n))
       end do ! n
       deallocate(this%face_neighbours_z)
    end if
    this%is_defined = .false.
  end subroutine deallocateScartInversionGrid
!------------------------------------------------------------------------
!> \brief return overall number of invgrid cells, if any
!
  function getNcellScartInversionGrid(this) result(ncell)
    type (scart_inversion_grid), intent(in) :: this
    integer :: ncell
    if(this%is_defined) then
       ncell = this%ncell
    else
       ncell = 0
    end if
  end function getNcellScartInversionGrid
!------------------------------------------------------------------------
!> \brief transform global to local coordinates of many points
!
  subroutine transformVectorGlobalToLocalScartInversionGrid(this,x,y,z,n)
    type (scart_inversion_grid) :: this
    integer :: n
    real, dimension(n) :: x,y,z
    real, dimension(n) :: xg,yg
    xg = x; yg = y
    ! shift poins by -cx,-cy,-zmax
    xg = xg - this%cx
    yg = yg - this%cy
    z = z - this%zmax
    ! rotate back to local xy by Mrot
    x = this%Mrot(1,1)*xg + this%Mrot(1,2)*yg
    y = this%Mrot(2,1)*xg + this%Mrot(2,2)*yg
  end subroutine transformVectorGlobalToLocalScartInversionGrid
!------------------------------------------------------------------------
!> \brief transform local to global coordinates of one point
!
  subroutine transformPointLocalToGlobalScartInversionGrid(this,x,y,z)
    type (scart_inversion_grid) :: this
    real :: x,y,z
    real :: xl,yl,zl
    xl = x; yl = y; zl = z
    ! rotate xl,yl back to real xy coords by transpose(Mrot)
    x = this%Mrot(1,1)*xl + this%Mrot(2,1)*yl
    y = this%Mrot(1,2)*xl + this%Mrot(2,2)*yl
    ! then shift point by cx,cy,zmax
    x = x + this%cx
    y = y + this%cy
    z = zl + this%zmax
  end subroutine transformPointLocalToGlobalScartInversionGrid
!------------------------------------------------------------------------
!> \brief transform wp,event or station coords to x,y,z coords for vtk application
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station' and that c1,c2,c3 are associated and have all same length
!!  also it can be assured that this scart inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here)
!! \param this scart inversion grid
!! \param c1 vector or first coordinate (contains vtk x-values on exit)
!! \param c2 vector or second coordinate (contains vtk y-values on exit)
!! \param c3 vector or third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!! \param errmsg error message
!
  subroutine transformToVtkScartInversionGrid(this,c1,c2,c3,coords_type,errmsg)
    type (scart_inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
!
    call addTrace(errmsg,'transformToVtkScartInversionGrid')
!
    ! no need of selecting coords_type: always do the same, since in case of using the 
    ! scart inversion grid the wavefield points as well as event and station coordinates
    ! are expected as x,y,z coordinates (in that order)
!
    if(this%use_local_coords_for_vtk) call transformVectorGlobalToLocalScartInversionGrid(this,c1,c2,c3,size(c1))
    if(this%apply_vtk_coords_scaling_factor) then
       c1 = c1 * this%vtk_coords_scaling_factor
       c2 = c2 * this%vtk_coords_scaling_factor
       c3 = c3 * this%vtk_coords_scaling_factor
    end if
  end subroutine transformToVtkScartInversionGrid
!------------------------------------------------------------------------
!> \brief return geometry information on cells for vtk output
!
  subroutine getGeometryVtkScartInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
       errmsg,cell_indx_req,indx_map_out)
    type (scart_inversion_grid) :: this
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    integer :: geometry_type
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
    ! local
    character(len=32) :: myname = 'getGeometryVtkScartInversionGrid'
    character(len=400) :: errstr
    integer, dimension(:), pointer :: indx_map
    integer, dimension(:), allocatable :: indx_map_tmp
    integer :: ilay,ix,nx,iy,ny,ncell,ncell_return,global_cell_count,icell_con,ipoints,i
    real, dimension(:), pointer :: x,y
    real :: xmin,xmax,ymin,ymax,zmin,zmax
    logical :: select_cell_indices
!
    nullify(indx_map,x,y)
!
    call addTrace(errmsg,myname)
    nullify(points,cell_connectivity,cell_type,cell_indx_out)
    if(present(indx_map_out)) nullify(indx_map_out)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    ! only if cell_indx_req is present and there are any indices in valid range, select those specific cells.
    ! otherwise return no cells (nullified pointers)
    if(present(cell_indx_req)) then
       if(any(cell_indx_req .ge. 1 .and. cell_indx_req .le. this%ncell)) then
          select_cell_indices = .true.
          ncell_return = 0
          ! define mapping indx_map_tmp which maps invgrid index to index in array cell_indx_req
          allocate(indx_map_tmp(this%ncell)); indx_map_tmp(:) = -1 ! indicate missing invgrid indices by -1
          do i = 1,size(cell_indx_req)
             if(cell_indx_req(i)>=1 .and. cell_indx_req(i)<=this%ncell) then
                ! only if the index did not occurr so far (i.e. value -1) , take it.
                ! this way, always the FIRST occuring index (in case of duplicate indices) is returned, 
                ! the others are lost. But in this case, you should not use the indx_map anyway
                if(indx_map_tmp(cell_indx_req(i)) == -1) then
                   indx_map_tmp(cell_indx_req(i)) = i
                   ncell_return = ncell_return + 1
                end if
             endif
          enddo ! i
       else
          write(errstr,*) "there are no valid indices among requested cell indices; indices must be between 1 and ",&
               this%ncell
          call add(errmsg,2,errstr,myname)
          return
       endif
    else
       select_cell_indices = .false.
       ! define index mapping anyway (although it is trivially the identiy in this case), 
       ! to provide optional output indx_map whenever requested (in case .not.present(cell_indx_req), indx_map is just the identity)
       allocate(indx_map_tmp(this%ncell))
       indx_map_tmp = (/ (i,i=1,this%ncell) /)
       ncell_return = this%ncell
    endif
!
    ! allocate for number of cells to be returned
    ! as all cells are hexahedra (for now), cell_connectivity((8+1)*ncell), 8 points + number of points
    ! for multiple cell types (e.g. hexa AND tets in invgrid), adjust here
    select case(this%vtk_geometry_type_int)
    case(0) ! CELLS
       allocate(points(3,8*ncell_return),cell_connectivity((8+1)*ncell_return),cell_type(ncell_return),&
            cell_indx_out(ncell_return),indx_map(ncell_return))
       geometry_type = 0
    case(1) ! CELL CENTERS
       allocate(points(3,ncell_return),cell_indx_out(ncell_return),indx_map(ncell_return))
       geometry_type = 1
    end select
!
    ncell = 0 ! counts the number of cells which will be returned (index of cell_type array)
    ipoints = 0 ! counts the number of points which will be returned (first index of points array)
    icell_con = 0 ! counts the number of entries in the cell_connectivity array
    global_cell_count = 0 ! globally counts the number of cells (from 1 to this%ncell)
!
    do ilay = 1,this%nlay
       nx = this%nx(this%iblock(ilay))
       x => getVectorPointer(this%x(this%iblock(ilay)))
       ny = this%ny(this%iblock(ilay))
       y => getVectorPointer(this%y(this%iblock(ilay)))
!
       zmin = this%z(ilay+1); zmax = this%z(ilay)
!
       do iy = 1,ny
          ymin = y(iy); ymax = y(iy+1)
          do ix = 1,nx
             xmin = x(ix); xmax = x(ix+1)
             global_cell_count = global_cell_count + 1
!
             if(select_cell_indices) then
                if(indx_map_tmp(global_cell_count) <= 0 ) cycle
             endif
!
             ncell=ncell+1
!
             select case(this%vtk_geometry_type_int)
             case(0) ! VOLUMETRIC CELLS
                cell_indx_out(ncell) = global_cell_count ! store the corresponding inversion grid cell index for this cell
                indx_map(ncell) = indx_map_tmp(global_cell_count) ! store the corresponding index in icoming array cell_indx_req (if present)
!
                cell_type(ncell) = 12 ! cell type is Hexahedron (compare vtk source code file vtkCellType.h)
                icell_con=icell_con+1; cell_connectivity(icell_con) = 8 ! first store number of point indices to come for this cell
!
                ! define the required points for the cell type of this cell (Hexahedron)
                ipoints=ipoints+1; points(:,ipoints) = (/ xmin,ymin,zmin /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1 ! indices in vtk files have offset 0, so always store ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmax,ymin,zmin /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmax,ymax,zmin /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmin,ymax,zmin /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmin,ymin,zmax /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmax,ymin,zmax /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmax,ymax,zmax /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
                ipoints=ipoints+1; points(:,ipoints) = (/ xmin,ymax,zmax /)
                icell_con=icell_con+1; cell_connectivity(icell_con) = ipoints -1
!
             case(1) ! CELL CENTER POINTS
                cell_indx_out(ncell) = global_cell_count ! store the corresponding inversion grid cell index for this cell
                indx_map(ncell) = indx_map_tmp(global_cell_count) ! store the corresponding index in icoming array cell_indx_req (if present)
!
                ! store the center point of this cell
                ipoints=ipoints+1; points(:,ipoints) = (/ 0.5*(xmax+xmin),0.5*(ymax+ymin),0.5*(zmax+zmin) /)
             end select ! case(this%vtk_geometry_type_int)
!
          end do ! ix
       end do ! iy
    end do ! ilay
!
    if(present(indx_map_out)) then
       indx_map_out => indx_map
    else
       deallocate(indx_map)
    endif
    deallocate(indx_map_tmp)
!
    if(.not.this%use_local_coords_for_vtk) then
       ! transform all points to global coordinates (laterally rotated and shifted)
       do ipoints = 1,size(points,2)
          call transformPointLocalToGlobalScartInversionGrid(this,points(1,ipoints),points(2,ipoints),points(3,ipoints))
       end do
    end if
    if(this%apply_vtk_coords_scaling_factor) then
       points = points*this%vtk_coords_scaling_factor
    end if
!
  end subroutine getGeometryVtkScartInversionGrid
!------------------------------------------------------------------------
!> \brief return indices of all face neighbours for all (optionally only subset of) cells
!!  \details If boundary_conditions is set, there should be a special form of nb_idx returned, which is used
!!  e.g. to define special smoothing conditions:
!!   'no_nb_inner_bnd': with some neighbours removed (which should not be smoothed with on internal invgrid 
!!                      boundaries). FOR scartInversionGrid GRIDS, THERE IS NO POSSIBILITY YET TO DEFINE INTERNAL
!!                      INVGRID BOUNDARIES!! SO IF THIS IS VALUE IS SET, NO NEIGHBOURS ARE RETURNED! 
!!                      In the future might introduce this functionality here. It should
!!                      be easy to realize this by introducing certain layers below which there should be an internal boundary.
!!                      (might be even easier to define a set of refinement blocks via inidices, below which there should be 
!!                      an internal boundary)
!!   'extra_nbs_outer_bnd': additional fake neighbours (having cell index 0) in case of zero boundary conditions 
!!                          on outer invgrid boundaries.
!!   'extra_nbs_outer_bnd_except_free_surface': same as 'extra_nbs_outer_bnd', but not applied for cells on free surfaces.
!!   '','standard': no special boundary handling, i.e. standard functionality (exactly the geometrical face neighbours are given)
!! \param this scart inversion grid
!! \param nb_idx pointer to array of length this%ncell which contains vector pointer to face neighbour indices
!! \param boundary_conditions optional string indicating type of boundary conditions for which neighbours should be returned
!
  subroutine getIndicesFaceNeighboursScartInversionGrid(this,nb_idx,boundary_conditions)
    type (scart_inversion_grid) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
    character(len=*), optional :: boundary_conditions
    ! local
    integer :: icell,nnb,nnbzmin,nnbzmax,nnbxy,inb
    integer, dimension(:), pointer :: nb,nb_zmin,nb_zmax
    character(len=39) :: bnd_cond
!
    nullify(nb,nb_zmin,nb_zmax)
!
    nullify(nb_idx)
    if(.not.this%is_defined) return
!
    if(present(boundary_conditions)) then
       select case(boundary_conditions)
       case('')
          bnd_cond = 'standard'
       case('extra_nbs_outer_bnd','extra_nbs_outer_bnd_except_free_surface','standard')
          bnd_cond = boundary_conditions
       case default
          return
       end select
    else
       ! By setting bnd_cond = 'standard' , any specific boundary handling is ignored below, i.e. 
       ! exactly the geometrical neighbours are given. This is the standard functionality
       bnd_cond = 'standard'
    end if
!
    allocate(nb_idx(this%ncell))
    nullify(nb)
    do icell = 1,this%ncell
       nb_zmin => getVectorPointer(this%face_neighbours_z(1,icell))
       nb_zmax => getVectorPointer(this%face_neighbours_z(2,icell))
!
       nnbzmin = 0; nnbzmax = 0
       if(associated(nb_zmin)) nnbzmin = size(nb_zmin)
       if(associated(nb_zmax)) nnbzmax = size(nb_zmax)
       nnbxy = count(this%face_neighbours_xy(:,icell) > 0)
!
       ! Dependent on the type of boundary condition, define the neighbour index vector nb
       select case(bnd_cond)
       case ('standard')
          ! give excatly the geometrical face neighbours of this cell
          nnb = nnbzmin + nnbzmax + nnbxy
          if(nnb > 0) then
             allocate(nb(nnb))
             inb = 0
             if(nnbxy > 0) then
                nb(inb+1:inb+nnbxy) = pack(this%face_neighbours_xy(:,icell),this%face_neighbours_xy(:,icell) > 0)
                inb = inb+nnbxy
             end if
             if(nnbzmin > 0) then
                nb(inb+1:inb+nnbzmin) = nb_zmin
                inb = inb + nnbzmin
             end if
             if(nnbzmax > 0) then
                nb(inb+1:inb+nnbzmax) = nb_zmax
                inb = inb + nnbzmax
             end if
          end if ! nnb > 0
!
       case ('extra_nbs_outer_bnd','extra_nbs_outer_bnd_except_free_surface')
          ! always has 4 lateral neighbours, even if on an outer boundary
          nnb = 4 
          ! if has neighbours below, add those, otherwise add one artificial
          if(nnbzmin > 0) then
             nnb = nnb + nnbzmin
          else
             nnb = nnb + 1
          end if
          ! If has neighbours above, add those, otherwise add one artificial.
          if(nnbzmax > 0) then
             nnb = nnb + nnbzmax
          else
             ! FOR NEIGHBOURS ABOVE THIS CELL, ONLY ADD ARTIFICIAL NEIGHBOUR FOR 'extra_nbs_outer_bnd'.
             ! IN CASE OF 'extra_nbs_outer_bnd_except_free_surface', ONLY GEOMETRICAL NEIGHBOURS HERE!
             if(bnd_cond=='extra_nbs_outer_bnd') nnb = nnb + 1
          end if
          allocate(nb(nnb))
          ! initialize vector of neighbour indices by -1 indicating artificial neighbours
          nb = -1
          ! Then overwrite (partially) by all actual geometrical neighbours.
          ! First add lateral face neighbours. As this%face_neighbours_xy is -1 for non existing neighbours, 
          ! it is consistent with what we want here, so simply write to first 4 entries in vector nb.
          nb(1:4) = this%face_neighbours_xy(:,icell)
          inb = 4
          ! Then add z neighbours. If there are no z neighbours (above, below), do nothing, the rest of the
          ! vector nb keeps the value -1 then.
          if(nnbzmin > 0) then
             nb(inb+1:inb+nnbzmin) = nb_zmin
             inb = inb + nnbzmin
          end if
          if(nnbzmax > 0) then
             nb(inb+1:inb+nnbzmax) = nb_zmax
             inb = inb + nnbzmax
          end if
       end select ! case bnd_cond
!
       ! Finally associate the neighbour index vector pointer for this cell.
       ! Even if nb was not allocated above, it works here (nothing happens).
       call associateVectorPointer(nb_idx(icell),nb)
       nullify(nb)
    end do ! icell
!
  end subroutine getIndicesFaceNeighboursScartInversionGrid
!------------------------------------------------------------------------
!> \brief for each cell return indices of wavefield points contained in that cell
!! \param this scart inversion grid
!! \param x vector or first coordinate of wavefield points
!! \param y vector or second coordinate of wavefield points
!! \param z vector or third coordinate of wavefield points
!! \param wp_idx pointer to array of length .ncell.this; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideScartInversionGrid(this,x,y,z,wp_idx,errmsg)
    type (scart_inversion_grid) :: this
    real, dimension(:), intent(in) :: x,y,z
    type (integer_vector_pointer), dimension(:), pointer :: wp_idx
    type (error_message) :: errmsg
    ! local
    character(len=32) :: myname = 'locateWpInsideScartInversionGrid'
    character(len=400) :: errstr
    real, dimension(:), allocatable :: xl,yl,zl
    integer, dimension(:), allocatable :: wp_nidx
    integer :: iwp,nwp,icell
    integer, dimension(:), pointer :: idx
!
    nullify(idx)
!
    call addTrace(errmsg,myname)
    nullify(wp_idx)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    allocate(wp_idx(this%ncell))
    allocate(wp_nidx(this%ncell)); wp_nidx(:) = 0
!
    ! in routine locateWpInsideInversionGrid of module inversionGrid it was already assured
    ! that x,y,z contain values and are all of same length!
    nwp = size(x)
    allocate(xl(nwp),yl(nwp),zl(nwp))
    xl = x; yl = y; zl = z
    call transformVectorGlobalToLocalScartInversionGrid(this,xl,yl,zl,nwp)
!
    do iwp = 1,nwp
       icell = locatePointInCellScartInversionGrid(this,xl(iwp),yl(iwp),zl(iwp))
       if(icell<0) cycle
!
       ! if current number of wavefield points of this cell wp_nidx(icell) is also 
       ! size of pointer in wp_idx(icell), reallocate
       if(mod(wp_nidx(icell),500) == 0) then
          idx => getVectorPointer(wp_idx(icell))
          idx => reallocate(idx,wp_nidx(icell)+500)
          call associateVectorPointer(wp_idx(icell),idx)
       end if
!
       ! add this wavefield point index iwp to list of indices of cell icell
       wp_nidx(icell) = wp_nidx(icell) + 1
       call fillVectorPointer(wp_idx(icell),(/iwp/),wp_nidx(icell))
    end do ! iwp
!
    ! reallocate wp_idx pointer for all cells to the actual number of wp indices found
    do icell = 1,this%ncell
       if(wp_nidx(icell) /= 0) then
          idx => getVectorPointer(wp_idx(icell))
          idx => reallocate(idx,wp_nidx(icell))
          call associateVectorPointer(wp_idx(icell),idx)
       end if
    end do ! icell
!
    if(nwp-sum(wp_nidx) > 0) then
       write(errstr,*) nwp-sum(wp_nidx)," wavefield points (out of ",nwp,&
            ") could not be located inside the inversion grid"
       call add(errmsg,1,errstr,myname)
    end if
!
    deallocate(wp_nidx,xl,yl,zl)
  end subroutine locateWpInsideScartInversionGrid
!------------------------------------------------------------------------
!> \brief for given point in LOCAL coordinates find index of inversion grid cell containing this point
!! \param this scart inverison grid
!! \param x x-coordinate of point
!! \param y y-coordinate of point
!! \param z z-coordinate of point
!! \param icell index of inversion grid cell containing point (x,y,z)
!! \return index of inversion grid cell containing the given point
!
  function locatePointInCellScartInversionGrid(this,x,y,z) result(icell)
    type (scart_inversion_grid) :: this
    real :: x,y,z
    integer :: icell
    ! local
    integer :: ix,iy,iz,iblock
!
    icell = -1
!
    iz = locateCoordinateScartInversionGrid(this%z,this%nlay+1,z,is_ascending=.false.)
    if(iz<0) return
!
    iblock = this%iblock(iz)
!
    ix = locateCoordinateScartInversionGrid(getVectorPointer(this%x(iblock)),this%nx(iblock)+1,x,is_ascending=.true.)
    if(ix<0) return
!
    iy = locateCoordinateScartInversionGrid(getVectorPointer(this%y(iblock)),this%ny(iblock)+1,y,is_ascending=.true.)
    if(iy<0) return
!
    icell = this%ncell_layer_cum(iz) + (iy-1)*this%nx(iblock) + ix
  end function locatePointInCellScartInversionGrid
!------------------------------------------------------------------------
!> \brief return index of interval in coords array which contains c
!> \details regardless of array coords being ascending or descending, two neighbouring entries
!!  coords(j),coords(j+1) of array coords define an interval (which defined a cell in this coordinate direction)
!!  which ALWAYS contains the left value coords(j) and does NOT contain the right value coords(j+1), 
!!  EXCEPT for the right-most interval defined by coords(n-1),coords(n) which always contains BOTH boundary 
!!  values coords(n-1) and coords(n)
!! \param coords array of coordinates (x,y,z), is assumed to be either in ascending or descending order!
!! \param n size of coords
!! \param c coordinate to be located inside intervals defined by values in coords
!! \param is_ascending indicates whether coords is in ascending order or not
!! \param i if i=-1 on exit, c is outside range of coords. Otherwise 1<=i<=(n-1) is index of interval containing c
! 
  function locateCoordinateScartInversionGrid(coords,n,c,is_ascending) result(i)
    integer :: n
    real, dimension(n) :: coords
    real :: c
    logical :: is_ascending
    ! returning
    integer :: i
    ! local
    integer :: jl,jr,jm
!
    if(is_ascending) then
!
       if(c < coords(1) .or. c > coords(n)) then
          i = -1
          return
       end if
       jl = 1; jr = n
       do while(jr-jl>1)
          jm = (jr+jl)/2
          if(c >= coords(jm)) then
             jl = jm
          else
             jr = jm
          end if
       end do ! while(jr-jl>1)
!
    else ! is_ascending
!
       if(c > coords(1) .or. c < coords(n)) then
          i = -1
          return
       end if
       jl = 1; jr = n
       do while(jr-jl>1)
          jm = (jr+jl)/2
          if(c <= coords(jm)) then
             jl = jm
          else
             jr = jm
          end if
       end do ! while(jr-jl>1)
!
    end if ! is_ascending
!
    i = jl
  end function locateCoordinateScartInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in cell icell to standard cell and compute their jacobian
!! \param this scart inversion grid
!! \param icell global inversion grid cell index
!! \param x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param type_standard_cell defines on return the shape of the standard cell (specific integration weights routine can be chosen): (4=Tetrahedron,6=Hexahedron). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param errmsg error message
!
  subroutine transformToStandardCellScartInversionGrid(this,icell,x,y,z,jacobian,type_standard_cell,errmsg)
    type (scart_inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    type (error_message) :: errmsg
    ! local
    character (len=41) :: myname = 'transformToStandardCellScartInversionGrid'
    character(len=400) :: errstr
    integer :: np
    real :: xmin,xmax,ymin,ymax,zmin,zmax
    double precision :: dx,dy,dz
!
    call addTrace(errmsg,myname)
!
    if(type_standard_cell == -1) then
       call add(errmsg,2,"Incoming value of type_standard_cell is -1, indicating a request for total "//&
            "integration weights (e.g. used by integration weights of type 6) instead of jacobian values. "//&
            "This functionality is not supported by inversion grids of type scartInversionGrid",myname)
       return
    end if
!
    ! in routine transformToStandardCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined, icell is valid and that x,y,z contain values and are all of same length!
!
    np = size(x)
    call transformVectorGlobalToLocalScartInversionGrid(this,x,y,z,np)
!
    call getCellBoundariesScartInversiongrid(this,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    dx = dble(xmax)-dble(xmin); dy = dble(ymax)-dble(ymin); dz = dble(zmax)-dble(zmin)
    ! transform x,y,z to hexahedral standard cell [-1,1]^3
    x = real( (2.d0*x -dble(xmax)-dble(xmin))/dx )
    y = real( (2.d0*y -dble(ymax)-dble(ymin))/dy )
    z = real( (2.d0*z -dble(zmax)-dble(zmin))/dz )
!
    if(any(x<-1.) .or. any(x>1.) .or. any(y<-1.) .or. any(y>1.) .or. any(z<-1.) .or. any(z>1.)) then
       write(*,*) "ERROR in transformToStandardCellScartInversionGrid :  xmin, xmax, x(:)  =  ",xmin,xmax,x
       write(*,*) "ERROR in transformToStandardCellScartInversionGrid :  ymin, ymax, y(:)  =  ",ymin,ymax,y
       write(*,*) "ERROR in transformToStandardCellScartInversionGrid :  zmin, zmax, z(:)  =  ",zmin,zmax,z
       write(errstr,*) "these wavefield points (supposed to be contained in cell ",icell,&
            ") transform to outside the hexahedral standard cell [-1,1]^3"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    jacobian(:) = 0.125*dx*dy*dz
!
    type_standard_cell = 6
!
  end subroutine transformToStandardCellScartInversionGrid
!------------------------------------------------------------------------
!! \brief for given cell index icell (assumed to be valid!!) return boundary values of cell
!
  subroutine getCellBoundariesScartInversiongrid(this,icell,xmin,xmax,ymin,ymax,zmin,zmax)
    type (scart_inversion_grid) :: this
    integer :: icell
    real :: xmin,xmax,ymin,ymax,zmin,zmax
    ! local
    integer :: ix,iy,iz,iblock
    real, dimension(:), pointer :: rp
!
    nullify(rp)
!
    call inverseCellIndexScartInversionGrid(this,icell,ix,iy,iz,iblock)
!
    rp => getVectorPointer(this%x(iblock))
    xmin = rp(ix); xmax = rp(ix+1)
!
    rp => getVectorPointer(this%y(iblock))
    ymin = rp(iy); ymax = rp(iy+1)
!
    zmin = this%z(iz+1); zmax = this%z(iz)
  end subroutine getCellBoundariesScartInversiongrid
!------------------------------------------------------------------------
!! \brief for given cell index icell (assumed to be valid!!) return x,y,z,block-indices
!
  subroutine inverseCellIndexScartInversionGrid(this,icell,ix,iy,iz,iblock)
    type (scart_inversion_grid) :: this
    integer :: icell,ix,iy,iz,iblock
    ! local
    integer :: jl,jr,jm,jcell_remain
!
    ! locate icell in array ncell_layer_cum -> yield iz (and iblock by iz)
    if(icell>this%ncell_layer_cum(this%nlay)) then
       iz = this%nlay
    else
       jl = 1; jr = this%nlay
       do while(jr-jl>1)
          jm = (jr+jl)/2
          if(icell > this%ncell_layer_cum(jm)) then
             jl = jm
          else
             jr = jm
          end if
       end do ! while(jr-jl>1)
       iz = jl
    end if
    iblock = this%iblock(iz)
!
    ! compute ix,iy from mod operations on layer
    jcell_remain = icell - this%ncell_layer_cum(iz)
!
    ix = mod(jcell_remain,this%nx(iblock))
    if(ix==0) then
       ix = this%nx(iblock)
       iy = jcell_remain/this%nx(iblock)
    else
       iy = jcell_remain/this%nx(iblock) + 1
    end if
!
  end subroutine inverseCellIndexScartInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param volume volume of cell icell
!
  subroutine getVolumeCellScartInversionGrid(this,icell,volume)
    type (scart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax
!
    ! in routine getVolumeCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell is valid
!
    call getCellBoundariesScartInversiongrid(this,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    volume = (xmax-xmin)*(ymax-ymin)*(zmax-zmin)
  end subroutine getVolumeCellScartInversionGrid
!------------------------------------------------------------------------
!> \brief get center of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which center should be returned
!! \param xc first coordinate of center of cell icell
!! \param yc second coordinate of center of cell icell
!! \param zc third coordinate of center of cell icell
!
  subroutine getCenterCellScartInversionGrid(this,icell,xc,yc,zc)
    type (scart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: xc,yc,zc
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax
!
    ! no need of selecting coords_type: 
    ! always do the same, since in case of using the 
    ! scart inversion grid the wavefield points as well as event and station coordinates
    ! are expected as x,y,z coordinates (in that order)
!
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell is valid
!
    call getCellBoundariesScartInversiongrid(this,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    xc = 0.5*(xmax+xmin)
    yc = 0.5*(ymax+ymin)
    zc = 0.5*(zmax+zmin)
!
    call transformPointLocalToGlobalScartInversionGrid(this,xc,yc,zc)
  end subroutine getCenterCellScartInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which radius should be returned
!! \param radius radius of cell icell
!
  subroutine getRadiusCellScartInversionGrid(this,icell,radius)
    type (scart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
    ! local
    real :: xmin,xmax,ymin,ymax,zmin,zmax
!
    ! in routine getRadiusCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined and icell is valid
!
    call getCellBoundariesScartInversiongrid(this,icell,xmin,xmax,ymin,ymax,zmin,zmax)
!
    radius = 0.5*sqrt((xmax-xmin)*(xmax-xmin) + (ymax-ymin)*(ymax-ymin) + (zmax-zmin)*(zmax-zmin))
  end subroutine getRadiusCellScartInversionGrid
!------------------------------------------------------------------------
!> \brief function answers whether a given point inside the inversion grid domain
!! \param c1 first coordinate
!! \param c2 second coordinate
!! \param c3 third coordinate
!! \param coords_type 'wp','event','station'
!
  function pointInsideScartInversionGrid(this,c1,c2,c3,coords_type) result(l)
    type (scart_inversion_grid) :: this
    real, intent(in) :: c1,c2,c3
    character(len=*) :: coords_type
    logical :: l
    real :: c1_copy,c2_copy,c3_copy
    real, dimension(:), pointer :: cx,cy
!
    nullify(cx,cy)
!
    l = .false.
!
    if(.not.this%is_defined) return
!
    c1_copy = c1
    c2_copy = c2
    c3_copy = c3
!
    ! no need of selecting coords_type: always do the same, since in case of using the 
    ! scart inversion grid the wavefield points as well as event and station coordinates
    ! are expected as x,y,z coordinates (in that order)
!
    call transformVectorGlobalToLocalScartInversionGrid(this,(/c1_copy/),(/c2_copy/),(/c3_copy/),1)
!
    ! check depth (z) -cordinate first
    if(c3_copy > this%z(1) .or. c3_copy < this%z(this%nlay+1)) return
!
    ! now check the first refinement block for lateral coordinates (actually this is sufficient, 
    ! assuming all refinement blocks have the same lateral coverage
    ! IF YOU ADD A CONSTRUCTOR TO THIS MODULE WHICH CREATES INVERSION GRID WITH DIFFERENT LATERAL EXTENT
    ! IN THE REFINEMENT BLOCKS, YOU NEED TO LOOP ON ALL REFINEMENT BLOCKS HERE!!
    cx => getVectorPointer(this%x(1))
    if(c1_copy < cx(1) .or. c1_copy > cx(this%nx(1))) return
!
    cy => getVectorPointer(this%y(1))
    if(c2_copy < cy(1) .or. c2_copy > cy(this%ny(1))) return

    ! if the code comes here, all the checks above were successfull, so return true
    l = .true.
  end function pointInsideScartInversionGrid
!
end module scartInversionGrid
