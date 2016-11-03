!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!
!#################################################################################################
!##  
!##  A T T E N T I O N  ! !
!##  
!##  THIS MODULE HAS KEPT GROWING AND BECAME SOME KIND OF MONSTEROUS CONSTRUCT, WHICH IS
!##  PROBABLY CHALLENGING TO EXTEND/MAINTAIN !
!##  
!##  F. Schumacher, July 2016
!##  
!#################################################################################################
!
!> \brief multi-chunk spherical inversion grid, semi-regular base cells, supports sophisticated base cell refinement
!!
!! \details The chunks inversion grid consists of 1, 2, 3 or 6 chunks of a cubed sphere. No "center cube"
!!  supported, i.e. covering the whole globe is not possible. The chunk width is always 90 degrees, except for 1-chunk grids
!!  where smaller chunks are allowed. The inversion grid cells are equi-angularly distributed,
!!  resulting in cells which are more or less evenly distributed in size (by contrast to the 
!!  schunk inversion grid, where the lateral projections of the cells onto the tangential plane
!!  have an equi-distant distribution on the plane). 
!!  2-chunk grids consist of two neighbouring 90-degree chunks, 3 chunk grids of three chunks 
!!  which are ALL neighbours of each other (i.e. essentially half of the Earth, see manuals for image
!!  of chunk positioning or below this doxygen comment block in code file chunksInversionGrid.f90).
!!  Inside the chunks, the inversion grid cells are constructed on the basis of regularly distributed
!!  base cells in the fashion of the schunk inversion grid (or scart inversion grid), having a certain 
!!  number of refinement blocks in depth inside which a fixed regular lateral resolution and depth resolution 
!!  can be chosen. Thereafter, by certain mechanisms, the base cells can be locally refined (so far only a toy
!!  method of random cell refinement is implemented for illustration , other methods can now be added easily, 
!!  e.g. based on ray coverage of the data).
!!  From outside the module (e.g. from module wavefieldPoints)
!!  always global Cartesian coordinates in [km] are assumed as "wp" coordinates ("wavefield points"), 
!!  with 3-axis pointing towards the north pole, 1-axis towards the equator at lon=0 and 2-axis 
!!  towards the equator at lon=90 degrees. "event" coordinates are lat,lon (degrees) and depth in km, 
!!  "station" coordinates are lat,lon (degrees) and altitude in m.
!!  This coordinate system is referred to as "global".
!!  <BR> 
!!  Internally, the module handles two other sets of  coordinate systems:
!!  "flat" = non-dimensionalized values on the tangential planes of each chunk, w.r.t. an 
!!  Earth radius of 1 (values x,y,z in type chunks_inversion_grid refer to these coordinates, with z being
!!  RADIUS w.r.t. Earth radius 1). These coordinates are referred to (by arrays cell_ixmin, cell_ixmax,...)
!!  as defining coordinates for inversion grid cells in all chunks, and in this coordinate system all chunks
!!  "overlay" each other. By chunk-individual rotations (inverse, i.e. transposed, matrices of Mrot_local2flat), 
!!  the individual chunks are rotated to their position in "local" coordinate system, where all chunks are curved 
!!  and attached to each other, as illustrated on images in manuals (and immediately after this doxygen comment 
!!  block in code file chunksInversionGrid.f90). In "local" reference frame, the first chunk (chunkd index 1) is 
!!  centered on the north pole of Earth and the counter-clockwise azimuthal rotation about local z-axis by az_rot is
!!  NOT applied. In order to go from "local" to "global": apply matrix Mrot_azimuth (only acting on x and y) and 
!!  thereafter the inverse (i.e. transpose) rotation to Mrot_global2local (constant rotation valid for all chunks).
!!  <BR>
!!  For sake of vtk plotting ONLY, the coordinate projections are referred to as follows: "GLOBAL" is "wp" 
!!  Cartesian coordinates; "LOCAL_CURV" is "GLOBAL" but with chunk 1 centered at the north pole (and not at 
!!  clon,clat) and no azimuthal rotation about local z-axis by az_rot applied; "LOCAL_FLAT" is "LOCAL_CURV" but 
!!  with curvature removed and the individual chunks simply shifted in the x-y plane and positioned as illustrated 
!!  on images in manuals (and immediately after this doxygen 
!!  comment block ended in code file chunksInversionGrid.f90); "LOCAL_NORTH_CURV" and "LOCAL_NORTH_FLAT" are
!!  variants of "LOCAL_CURV" and "LOCAL_FLAT", respectively, where "NORTH" refers to the fact that in those cases 
!!  the inversion grid is kept rotated about the local vertical by az_rot and this means that the northing (as seen
!!  at the center point of the first chunk) is along a coordinate axis (namely the negative x-axis, i.e. x really 
!!  points to south in those cases, which can be nice for plotting and cutting plots in paraview).
!!  <BR>
!!  
!!  <BR>
!!  For programmer's orientation: The inversion grid cells are categorized as "internal" and "external" cells. 
!!  Outside this module (in executables and as vtk file data name etc.), normally only the terminology "cell index" 
!!  is used, which is equivalent to external call index! External cells
!!  are all cells that are used from the ASKI package as actual inversion grid cells (i.e. on which kernels are
!!  pre-integrated and which are referred to by data model space info objects, etc.). ADDITIONALLY there can be 
!!  cells for internal use only, which are not visible to the rest of the ASKI software package ("internal" cells).
!!  Originally, this construct was chosen to realize base cell refinement in a multi hierachical way (e.g. refining
!!  cells recursively several times) and only the lowest level of refined cells that do not have subcells should
!!  be used as external inversion grid cells. All other parent cells should be kept as internal cells, in order 
!!  to conveniently search for neighbours (going 1 hierarchical refinement level up, etc.). When eventually 
!!  programming the first base cell refinement method in July 2016 (Florian Schumacher), it was noted that such
!!  multiple levels of cell refinement complicate the module structure significantly (because cells are not handled
!!  recursively but iteratively/serially in arrays). So Florian decided to have only ONE hierarchy level of cell 
!!  refinement, namely subcells of base cells. Thus, the concept of "external" and "internal" cells is restricted
!!  only to the case that base cells can lose their "external" status when they are refined. All contained subcells
!!  are automatically external (as they are not refined themselves). If a base cell is not refined, it stays 
!!  external. Please note, however, that ALL cells have a reference as an INTERNAL cell and their internal cell
!!  index is their position in arrays cell_ixmin, cell_ixmax, ... etc. For ALL operations within this module
!!  chunksInversionGrid, ALWAYS the internal cell index is used, and only for communication to the rest of the
!!  ASKI package the external cell index is used (realized by mapping icell_internal (mapping external cell index
!!  to internal cell index) as well as mapping icell_external (mapping internal cell index to external cell index)).
!!  <BR>
!!  The internal cell indexing is as follows: It is assumed that the internal cell index is sorted such that all 
!!  base cells come first, i.e. have indices in range 1,..,ncell_base. WHAT IS REFERRED TO AS BASE CELL INDEX (also
!!  in extecutables, e.g. chunksInvgrid2vtk) IS THE INTERNAL CELL INDEX OF THE BASE CELLS. 
!!  The base cells of the first chunk constitute 
!!  the first portion of internal cells and they are sorted analogously to cells in inversion grids of type 
!!  schunk_inversion_grid or type scart_inversion_grid (i.e. outer loop z from (decreasing, i.e. from top to 
!!  bottom), intermediate loop on y (increasing), inner loop on x (increasing)). If there is more than one chunk, 
!!  the base cells of the other chunks (chunk 2, chunk 3, chunk 4, chunk 5, chunk 6) are appended to that list of
!!  base cells. The indexing of the base cells in the other chunks is equal to the order of base cells in the first 
!!  chunk. Note that for multi-chunk grids, ALL chunks are identical in their flat reference frame w.r.t. their
!!  base cells: inversion grids of more than 1 chunk MUST have wlat = wlon = 90.0 (for all chunks). For 2-chunk 
!!  inversion grids, the two (neighbouring) chunks are identical in terms of base cells, regardless of values 
!!  CHUNKS_INVGRID_BASE_NLAT and CHUNKS_INVGRID_BASE_NLON (chunk 2 is just shifted to the left) and the depth
!!  refinement is valid for the whole inversion grid.
!!  For inversion grids with 3 or more chunks, additionally CHUNKS_INVGRID_BASE_NLAT must equal
!!  CHUNKS_INVGRID_BASE_NLON for each depth layer. This is required, because at the boundary between chunk 2 and 
!!  chunk 3, for instance, the lon direction of chunk 2 equals the lat direction of chunk 3 (and if 
!!  CHUNKS_INVGRID_BASE_NLAT was not equal to CHUNKS_INVGRID_BASE_NLON, we would have two different lateral 
!!  base cell resolutions on the boundary, which physically would not be sensible). Therefore, also for nchunk >= 3
!!  all chunks are identical in terms of their base cells in their flat reference frame (and just rotated to their
!!  location in local reference frame). This is the reason, why the definition of the base cells of chunks >= 2
!!  is just a duplicate of the base cells of chunk 1 (i.e. the first ncell_base entries of arrays cell_ixmin,
!!  cell_ixmax, .. etc. repeat themselves nchunk-times). In case there is no base cell refinement, the base cells
!!  make up the inversion grid (all base cells are external cells) and the external indexing equals the internal
!!  cell index. In case there is an additional base cell refinement, the internal indexing of the base cells is NOT 
!!  changed and newly created cells (subcells of the base cells) are just appended to the list of internal cells
!!  (i.e. arrays cell_ixmin, cell_ixmax, .. etc. are just extended by new subcells during their creation).
!!  However, since a base cell is no longer an external cell when it is refined, it must be removed from the 
!!  external cell index. NOTE, THAT EXTERNAL CELL INDEXING MUST BE CONNECTED (from 1, .., to ncell) AND MUST NOT
!!  HAVE GAPS, I.E. MISSING INDICES! One way to achieve that, is to 
!!  replace the base cell (in terms of its external index) by the first subcell and append all other subcells
!!  (in terms of their external index) to the existing set of external inversion grid cells. Alternatively, the 
!!  external indices of all base cells (with higher external index) must be shifted some way (decreased by 1 in 
!!  order to fill the gap or increased in order to make room for the new subcells). The advantage of making room
!!  for the new cells and "squeezing" them into the place where the original base cell was, is that the (external)
!!  cell index remains ascending from base cell to base cell (which is actually nice for orientation by the user,
!!  however, not relevant for the functionality of ASKI). Therefore, this approach is preferred by the exemplary
!!  method EXPERIMENTAL_RANDOM_SUBDIVISION, for instance (as an example in the code of subroutine 
!!  doExpRandCellRefinementChunksInversionGrid, the method of replacing the base cell only and appending the rest
!!  of the new subcells is also given, in a commented code section -> this is slightly more convenient to implement
!!  / might have slightly better performance since no shifting operations are required).
!!
!! \author Florian Schumacher
!! \date July 2016
!
!  DISTRIBUTION OF CHUNKS in chunksInversionGrid (LOCAL_FLAT projection here, i.e. x points down, y points right, 
!  z goes out of the screen):
!      +---+
!      | 3 |
!  +---+---+---+---+
!  | 2 | 1 | 5 | 4 |
!  +---+---+---+---+
!      | 6 |
!      +---+
!
module chunksInversionGrid
!
  use inputParameter
  use vectorPointer
  use realloc
  use mathConstants
  use errorMessage
!
  implicit none
!
  private :: constructNewChunksInversionGrid,readChunksInversionGrid,writeChunksInversionGrid,&
       readCheckParFileChunksInversionGrid,transformVectorFlatToGlobalChunksInversionGrid,&
       transformVectorGlobalToFlatChunksInversionGrid,&
       transformFlatToVtkChunksInversionGrid,computeTransformsChunksInversionGrid,&
       locateCoordinateChunksInversionGrid,sortCellBoundaryCoordinatesChunksInversionGrid,&
       findFaceNeighboursOfBaseCellsChunksInversionGrid,findFaceNeighboursOfRefinedCellsChunksInversionGrid,&
       computeCellRadiiChunksInversionGrid,privateGetGeometryVtkChunksInversionGrid,&
       doCellRefinementChunksInversionGrid,doExpRandCellRefinementChunksInversionGrid,&
       addCoordinatesChunksInversionGrid,getTouchingBoundaryOfNbChunkChunksInversionGrid,&
       locateWpLaterallyInsideChunksInversionGrid
!
  interface is_refined; module procedure containsRefinedCellsChunksInversionGrid ; end interface
  interface operator (.nbase.); module procedure getNcellBaseChunksInversionGrid; end interface
!
  !< derived type of complete inversion grid information (geometry and all cells)
  type chunks_inversion_grid
     private
     logical :: is_defined = .false. !< flag indicating the correct definition (initialization) of the object (i.e. all following values)
!
! GENERAL GEOMETRY INFORMATION OF CHUNKS
!
     integer :: nchunk = 0 !< number of chunks
!
     real :: clat = 0. !< Geographic latitude of center of the first (or only) chunk in degrees
     real :: clon = 0. !< Geographic longitude of center of the first (or only) chunk in degrees
     real :: wlat = 0. !< width of chunk(s) parallel to latitude in degrees (must be 90.0 for NCHUNK > 1)
     real :: wlon = 0. !< width of chunk(s) parallel to longitude in degrees (must be 90.0 for NCHUNK > 1)
!
     real :: rmax = 0. !< maximum radius of sphere [usually km]
     real :: az_rot = 0. !< angle in degrees of counter-clockwise azimuthal rotation about local z-axis
!
     double precision, dimension(2,2) :: Mrot_azimuth !< azimuthal rotation by angle rot in the x-y plane
     double precision, dimension(3,3) :: Mrot_global2local !< rotation from global cartesian to local curved chunks
     double precision, dimension(:,:,:), pointer :: Mrot_local2flat => null() !< rotation from local curved chunks to flat reference chunk; (3,3,nchunk)-size array containing for each chunk different transformation matrices
     double precision, dimension(:,:,:), pointer :: Mrot_global2flat => null() !< rotation from global Cartesian chunks to flat reference chunk (i.e. Mrot_local2flat times Mrot_global2local); (3,3,nchunk)-size array containing for each chunk different transformation matrices
!
! INVERSION GRID CELLS, OUTSIDE COMMUNICATION
!
     integer :: ncell = 0 !< number of inversion grid cells in use for the grid (lowest level of subcells), this number is used outside this module
     integer, dimension(:), pointer :: icell_internal => null() !< index map (size ncell): icell_internal(icell) = internal_index_of_cell_icell , where 1 <= icell <= ncell
     integer, dimension(:), pointer :: icell_external => null() !< inverse map to icell_internal (size ncell_internal): if an internal cell is not external, the index value is -1 (always ncell_internal >= ncell)
!
! INTERNAL INVERSION GRID CELLS, BASE CELLS AND SUB(SUB..)CELLS
!
     integer :: ncell_internal = 0 !< total number of internal cells, i.e. including full hierarchy of refined cells starting with base cells, including all subcells, subsubcells...
     integer :: ncell_base = 0 !< number of base cells; it is assumed that the internal cell index is sorted such that all base cells come first, i.e. have indices in range 1..ncell_base
!
     integer :: nx = 0 !< number of x-coordinates in table x
     double precision, dimension(:), pointer :: x => null() !< table of unsorted x-coordinates in uncurved local chunk (w.r.t. radius 1); cells refer to index in this array for their bounding coordinates
     integer, dimension(:), pointer :: sorted2ix => null() !< mapping from sorted coordinate index to index in x: If the array x was sorted (increasingly), then a coordinate x(i) (denoting i = sorted2ix(j)) would be at position j in this sorted array of x-coordinates
     integer, dimension(:), pointer :: ix2sorted => null() !< ix2sorted is the inverse mapping of sorted2ix
     integer :: ny = 0 !< number of y-coordinates in table y
     double precision, dimension(:), pointer :: y => null() !< table of unsorted y-coordinates in uncurved local chunk (w.r.t. radius 1); cells refer to index in this array for their bounding coordinates
     integer, dimension(:), pointer :: sorted2iy => null() !< mapping from sorted coordinate index to index in y: If the array y was sorted (increasingly), then a coordinate y(i) (denoting i = sorted2iy(j)) would be at position j in this sorted array of y-coordinates
     integer, dimension(:), pointer :: iy2sorted => null() !< iy2sorted is the inverse mapping of sorted2iy
     integer :: nz = 0 !< number of z-coordinates in table z
     double precision, dimension(:), pointer :: z => null() !< table of unsorted z-coordinates in uncurved local chunk (z = radius, same unit as rmax); cells refer to index in this array for their bounding coordinates
     integer, dimension(:), pointer :: sorted2iz => null() !< mapping from sorted coordinate index to index in z: If the array z was sorted (increasingly), then a coordinate z(i) (denoting i = sorted2iz(j)) would be at position j in this sorted array of z-coordinates
     integer, dimension(:), pointer :: iz2sorted => null() !< iz2sorted is the inverse mapping of sorted2iz

     ! detailed geometry information on each internal cell (lateral and depth coverage)
     integer, dimension(:), pointer :: cell_ichunk => null() !< for each internal cell, contains the chunk index
     integer, dimension(:), pointer :: chunk_ncell_internal => null() !< for each chunk, contains the number of internal inversion grid cells (size of respective vectors in chunk_icell_internal
     type (integer_vector_pointer), dimension(:), pointer :: chunk_icell_internal => null() !< for each chunk i, a vector of length chunk_ncell_internal(i) giving all internal inversion grid cell indices of cells contained in chunk i
     integer, dimension(:), pointer :: cell_ixmin => null() !< (ncell_internal)-array of indizes of minimum x-coordinates (refers to position in array x)
     integer, dimension(:), pointer :: cell_ixmax => null() !< (ncell_internal)-array of indizes of maximum x-coordinates (refers to position in array x)
     integer, dimension(:), pointer :: cell_iymin => null() !< (ncell_internal)-array of indizes of minimum y-coordinates (refers to position in array y)
     integer, dimension(:), pointer :: cell_iymax => null() !< (ncell_internal)-array of indizes of maximum y-coordinates (refers to position in array y)
     integer, dimension(:), pointer :: cell_izmin => null() !< (ncell_internal)-array of indizes of minimum z-coordinates (refers to position in array z)
     integer, dimension(:), pointer :: cell_izmax => null() !< (ncell_internal)-array of indizes of maximum z-coordinates (refers to position in array z)
!
     ! cell center coordinates and radius of cell
     double precision, dimension(:,:), pointer :: cell_center => null() !< (3,ncell_internal)-array of cell center coordinates, global cartesian (x,y,z) per point
     double precision, dimension(:), pointer :: cell_radius => null() !< for each internal cell, its radius
!
     ! neighbours of all internal cells
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_xmin => null() !< for each internal cell, a vector of internal cell indices is stored which indicates all face neighbours of that cell on its xmin face
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_xmax => null() !< for each internal cell,indicates all face neighbours on its xmax face
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_ymin => null() !< for each internal cell,indicates all face neighbours on its ymin face
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_ymax => null() !< for each internal cell,indicates all face neighbours on its ymax face
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_zmin => null() !< for each internal cell,indicates all face neighbours on its zmin face
     type (integer_vector_pointer), dimension(:), pointer :: cell_nb_zmax => null() !< for each internal cell,indicates all face neighbours on its zmax face
!
     ! inheritance of cells, cell hierarchy
     integer, dimension(:), pointer :: cell_parent => null() !< for each internal cell: internal index of parent cell (this cell is subcell of parent cell). For base cells (i.e. no parent), value is -1
     type (integer_vector_pointer), dimension(:), pointer :: cell_subcells => null() !< for each base cell: vector of internal cell indices of subcells. If no subcells: points to null()
!
! SPECIFICATIONS FOR VTK OUTPUT
!
     character(len=16) :: vtk_projection = "" !< 'GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT', defining the vtk transformation
     double precision, dimension(:,:), pointer :: vtk_flat_shift => null() !< (2,nchunk)-array containing for each chunk a lateral shift vector (used for the "*FLAT" vtk projections)
     logical :: apply_vtk_coords_scaling_factor !< flag indicating whether coordinates for vtk output are scaled by vtk_coords_scaling_factor 
     double precision :: vtk_coords_scaling_factor !< factor by which coordinates for vtk output are scaled if apply_vtk_coords_scaling_factor is .true.
     integer :: vtk_geometry_type_int = -1 !< type of vtk geometry:  0 = volumetric cells , 1 = cell center points
!
! OTHER
!
     logical :: contains_refined_cells = .false. !< flag defining whether there are only base cells (.false.) or there was some cell refinement done (.true.)
!
     ! in the future: there could be flags in parameter file like: DONT_SMOOTH_LAYER_BOUNDARIES, 
     ! or SMOOTHING_BOUNDARY_CONDITIONS which could be taken into account here, and memorized for better handling 
     ! of smoothing conditions in calls to certain routines below
  end type chunks_inversion_grid
!
contains
!------------------------------------------------------------------------
!> \brief logical return whether this type of inversion grid is able to transform points of a specific coords type to vtk plot projection
!
  function canTransformToVtkPointsOutsideChunksInversionGrid(this,coords_type) result(l)
    type(chunks_inversion_grid) :: this
    character(len=*) :: coords_type
    logical :: l
    l = .false.
    ! the capability of transforming points to vtk is ONLY dependent on this%vtk_projection, 
    ! but not dependent on coords type. So do not select coords_type here, simply ignore it
    ! projections of type 'LOCAL_FLAT' and' LOCAL_NORTH_FLAT' are NOT supported for transformation
    ! to vtk, since it must be known in which chunk the points are (not given, will raise error)
    ! all other projections are possible
    if(.not.this%is_defined) return
    select case(this%vtk_projection)
    case('GLOBAL','LOCAL_CURV','LOCAL_NORTH_CURV')
       l = .true.
    end select
  end function canTransformToVtkPointsOutsideChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get unit factor of the volume element
!! \param this chunks inversion grid
!! \param uf_vol unit factor of volume element (return value of this subroutine)
!! \param errmsg error message
!
  subroutine getUnitFactorOfVolumeElementChunksInversionGrid(this,uf_vol,errmsg)
    type (chunks_inversion_grid) :: this
    real :: uf_vol
    type (error_message) :: errmsg
    character (len=47) :: myname = 'getUnitFactorOfVolumeElementChunksInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) call add(errmsg,1,"be aware that the inversion grid not yet defined; "//&
         "however, the unit factor of the volume element can be correctly computed at this point",myname)
!
    ! The chunks inversion grid assumes units for its spatial extension in km !
    ! This is independent of the unit of wavefield points
    uf_vol = 1.0e9
  end subroutine getUnitFactorOfVolumeElementChunksInversionGrid
!------------------------------------------------------------------------
!> \brief map vtk geometry type names to integers
!
  function intGeometryTypeChunksInversionGrid(vtk_geometry_type_str) result(vtk_geometry_type_int)
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
  end function intGeometryTypeChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get logical value whether this chunks inversion grid was refined or only contains base cells
!
  function containsRefinedCellsChunksInversionGrid(this) result(l)
    type(chunks_inversion_grid), intent(in) :: this
    logical :: l
    l = this%contains_refined_cells
  end function containsRefinedCellsChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get logical value whether this chunks inversion grid was refined or only contains base cells
!
  function getNcellBaseChunksInversionGrid(this) result(n)
    type(chunks_inversion_grid), intent(in) :: this
    integer :: n
    n = this%ncell_base
  end function getNcellBaseChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get chunk indices from internal cell indices; return null() if there are any invalid indices or inversion grid is not yet defined
!
  function getChunkIndexOfBaseCellsChunksInversionGrid(this,icell_base_in) result(p)
    type(chunks_inversion_grid), intent(in) :: this
    integer, dimension(:), intent(in) :: icell_base_in
    integer, dimension(:), pointer :: p
!
    nullify(p)
    if(.not.this%is_defined) return
    if(size(icell_base_in)<=0) return
    if(any(icell_base_in<1 .or. icell_base_in>this%ncell_base)) return
!
    allocate(p(size(icell_base_in)))
    p = this%cell_ichunk(icell_base_in)
  end function getChunkIndexOfBaseCellsChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get chunk indices from external cell indices; return null() if there are any invalid indices or inversion grid is not yet defined
!
  function getChunkIndexOfExternalCellsChunksInversionGrid(this,icell_external_in) result(p)
    type(chunks_inversion_grid), intent(in) :: this
    integer, dimension(:), intent(in) :: icell_external_in
    integer, dimension(:), pointer :: p
!
    nullify(p)
    if(.not.this%is_defined) return
    if(size(icell_external_in)<=0) return
    if(any(icell_external_in<1 .or. icell_external_in>this%ncell)) return
!
    allocate(p(size(icell_external_in)))
    p = this%cell_ichunk(this%icell_internal(icell_external_in))
  end function getChunkIndexOfExternalCellsChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get external cell indices from base cell indices; return null() if there are any invalid indices or inversion grid is not yet defined
!
  function getExternalCellIndexOfBaseCellsChunksInversionGrid(this,icell_base_in) result(p)
    type(chunks_inversion_grid), intent(in) :: this
    integer, dimension(:), intent(in) :: icell_base_in
    integer, dimension(:), pointer :: p
!
    nullify(p)
    if(.not.this%is_defined) return
    if(size(icell_base_in)<=0) return
    if(any(icell_base_in<1 .or. icell_base_in>this%ncell_base)) return
!
    allocate(p(size(icell_base_in)))
    p = this%icell_external(icell_base_in)
  end function getExternalCellIndexOfBaseCellsChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get external subcell indices for given base cell index; return null() if base cell index is invalid or inversion grid is not yet defined
!
  function getSubcellsOfBaseCellChunksInversionGrid(this,icell_base_in) result(p)
    type(chunks_inversion_grid), intent(in) :: this
    integer :: icell_base_in
    integer, dimension(:), pointer :: p
    ! local
    integer, dimension(:), pointer :: sub
!
    nullify(p)
    if(.not.this%is_defined) return
    if(icell_base_in<1 .or. icell_base_in>this%ncell_base) return
    sub => getVectorPointer(this%cell_subcells(icell_base_in))
    if(.not.associated(sub)) return
!
    allocate(p(size(sub)))
    p = this%icell_external(sub)
  end function getSubcellsOfBaseCellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief create chunks inversion grid
!
  subroutine createChunksInversionGrid(this,parfile,path,lu,errmsg,recreate)
    type(chunks_inversion_grid) :: this
    character (len=*) :: parfile,path
    integer :: lu
    type (error_message) :: errmsg
    logical, optional :: recreate
    ! local
    character(len=25) :: myname = 'createChunksInversionGrid'
    logical :: recreate_invgrid,file_exists
    ! parfile
    type (input_parameter) :: inpar
    character (len=80), dimension(19) :: inpar_keys
    data inpar_keys/'CHUNKS_INVGRID_GEOM_NCHUNK',&
         'CHUNKS_INVGRID_GEOM_RMAX','CHUNKS_INVGRID_GEOM_CLAT','CHUNKS_INVGRID_GEOM_CLON',&
         'CHUNKS_INVGRID_GEOM_WLAT','CHUNKS_INVGRID_GEOM_WLON','CHUNKS_INVGRID_GEOM_ROT',&
         'CHUNKS_INVGRID_BASE_NREF_BLOCKS','CHUNKS_INVGRID_BASE_NLAY','CHUNKS_INVGRID_BASE_THICKNESS',&
         'CHUNKS_INVGRID_BASE_NLAT','CHUNKS_INVGRID_BASE_NLON','CHUNKS_INVGRID_CREF_METHOD',&
         'CHUNKS_INVGRID_CREF_PARAMETERS','CHUNKS_INVGRID_FILE','VTK_PROJECTION',&
         'VTK_GEOMETRY_TYPE','SCALE_VTK_COORDS','VTK_COORDS_SCALING_FACTOR'/
!
    call addTrace(errmsg,myname)
    if(this%is_defined) then
       call add(errmsg,1,"this object is already defined, deallocating it now before creating new one",myname)
       call deallocateChunksInversionGrid(this)
    end if
!
    call createKeywordsInputParameter(inpar,inpar_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) return
!
    if(present(recreate)) then
       recreate_invgrid = recreate
    else
       ! by default do not recreate but try to read in an existing inversion grid file
       recreate_invgrid = .false.
    end if
    if(.not.recreate_invgrid) then
       inquire(file = trim(path)//trim(inpar.sval.'CHUNKS_INVGRID_FILE'), exist = file_exists) 
       recreate_invgrid = .not.file_exists
    end if
!
    if(recreate_invgrid) then
       call constructNewChunksInversionGrid(this,path,inpar,lu,errmsg)
    else
       call readChunksInversionGrid(this,path,inpar,lu,errmsg)
    end if
    if(.level.errmsg == 2) return
!
  end subroutine createChunksInversionGrid
!--------------------------------------------------------------------
!> \brief read chunks inversion grid parameter file, check availability and basic consistency
  subroutine readCheckParFileChunksInversionGrid(inpar,nchunk,clat,clon,wlat,wlon,rmax,az_rot,nblock,nlay_block,&
       thickness,nlat_block,nlon_block,cref_param,vtk_geometry_type_int,apply_vtk_scaling,vtk_factor,errmsg,myname)
    ! incoming
    type (input_parameter) :: inpar
    character(len=*) :: myname
    ! returning
    integer :: nchunk,nblock,iblock,vtk_geometry_type_int
    real :: clat,clon,wlat,wlon,rmax,az_rot
    integer, dimension(:), pointer :: nlay_block,nlat_block,nlon_block
    real, dimension(:), pointer :: thickness,cref_param
    logical :: apply_vtk_scaling
    double precision :: vtk_factor
    type (error_message) :: errmsg
    ! local
    integer :: ios
    character(len=400) :: errstr
    character(len=500) :: char_string
!
    nullify(nlay_block,nlat_block,nlon_block,thickness,cref_param)
!
    nchunk = ival(inpar,'CHUNKS_INVGRID_GEOM_NCHUNK',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'CHUNKS_INVGRID_GEOM_NCHUNK' = '"//trim(inpar.sval.&
            'CHUNKS_INVGRID_GEOM_NCHUNK')//"' of chunks invgrid parfile is not a valid integer value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    select case(nchunk)
    case(1,2,3,6)
       ! OK, do nothing
    case default
       write(errstr,*) "parameter 'CHUNKS_INVGRID_GEOM_NCHUNK' = ",nchunk,&
            " of chunks invgrid parfile is not valid; must be one either 1, 2, 3, or 6"
       call add(errmsg,2,trim(errstr),myname)
      return
    end select
!
    clat = rval(inpar,'CHUNKS_INVGRID_GEOM_CLAT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_CLAT' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_CLAT')//"'",myname)
       return
    end if
    if(clat > 90.0 .or. clat < -90.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_CLAT = ",clat,"; must be between -90.0 and 90.0"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    clon = rval(inpar,'CHUNKS_INVGRID_GEOM_CLON',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_CLON' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_CLON')//"'",myname)
       return
    end if
    if(clon > 360.0 .or. clon < -360.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_CLON = ",clon,"; must be between -360.0 and 360.0"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    rmax = rval(inpar,'CHUNKS_INVGRID_GEOM_RMAX',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_RMAX' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_RMAX')//"'",myname)
       return
    end if
    if(rmax <= 0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_RMAX = ",rmax,"; must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    wlat = rval(inpar,'CHUNKS_INVGRID_GEOM_WLAT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_WLAT' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_WLAT')//"'",myname)
       return
    end if
    if(wlat <= 0 .or. wlat > 90.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_WLAT = ",wlat,"; must be > 0 and <= 90.0"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(nchunk>1 .and. wlat/=90.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_WLAT = ",wlat,"; must be 90.0 in case of nchunk > 1"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    wlon = rval(inpar,'CHUNKS_INVGRID_GEOM_WLON',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_WLON' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_WLON')//"'",myname)
       return
    end if
    if(wlon <= 0 .or. wlon > 90.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_WLON = ",wlon,"; must be > 0 and <= 90.0"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(nchunk>1 .and. wlon/=90.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_WLON = ",wlon,"; must be 90.0 in case of nchunk > 1"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    az_rot = rval(inpar,'CHUNKS_INVGRID_GEOM_ROT',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read real value for 'CHUNKS_INVGRID_GEOM_ROT' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_GEOM_ROT')//"'",myname)
       return
    end if
    if(az_rot < -360.0 .or. az_rot > 360.0) then
       write(errstr,*) "CHUNKS_INVGRID_GEOM_ROT = ",az_rot,", i.e. more than one complete rotation. ",&
            "check if this is intended"
       call add(errmsg,1,errstr,myname)
       return
    end if
!
    nblock = ival(inpar,'CHUNKS_INVGRID_BASE_NREF_BLOCKS',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read integer value for 'CHUNKS_INVGRID_BASE_NREF_BLOCKS' from '"//&
            trim(inpar.sval.'CHUNKS_INVGRID_BASE_NREF_BLOCKS')//"'",myname)
       return
    end if
    if(nblock<1) then
       write(errstr,*) "CHUNKS_INVGRID_BASE_NREF_BLOCKS = ",nblock,"; must be greater than zero"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nlay_block => ivecp(inpar,'CHUNKS_INVGRID_BASE_NLAY',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read CHUNKS_INVGRID_BASE_NREF_BLOCKS = ",nblock,&
            " integers for 'CHUNKS_INVGRID_BASE_NLAY' from '"//trim(inpar.sval.'CHUNKS_INVGRID_BASE_NLAY')//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(any(nlay_block < 1)) then
       call add(errmsg,2,"all values of CHUNKS_INVGRID_BASE_NLAY must be greater than zero",myname)
       return
    end if
!
    thickness => rvecp(inpar,'CHUNKS_INVGRID_BASE_THICKNESS',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read CHUNKS_INVGRID_BASE_NREF_BLOCKS = ",nblock,&
            " integers for 'CHUNKS_INVGRID_BASE_THICKNESS' from '"//trim(inpar.sval.'CHUNKS_INVGRID_BASE_THICKNESS')//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(any(thickness <= 0.)) then
       call add(errmsg,2,"all values of CHUNKS_INVGRID_BASE_THICKNESS must be greater than zero",myname)
       return
    end if
!
    nlat_block => ivecp(inpar,'CHUNKS_INVGRID_BASE_NLAT',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read CHUNKS_INVGRID_BASE_NREF_BLOCKS = ",nblock,&
            " integers for 'CHUNKS_INVGRID_BASE_NLAT' from '"//trim(inpar.sval.'CHUNKS_INVGRID_BASE_NLAT')//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(any(nlat_block < 1)) then
       call add(errmsg,2,"all values of CHUNKS_INVGRID_BASE_NLAT must be greater than zero",myname)
       return
    end if
!
    nlon_block => ivecp(inpar,'CHUNKS_INVGRID_BASE_NLON',nblock,iostat=ios)
    if(ios /= 0) then
       write(errstr,*) "could not read CHUNKS_INVGRID_BASE_NREF_BLOCKS = ",nblock,&
            " integers for 'CHUNKS_INVGRID_BASE_NLON' from '"//trim(inpar.sval.'CHUNKS_INVGRID_BASE_NLON')//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(any(nlon_block < 1)) then
       call add(errmsg,2,"all values of CHUNKS_INVGRID_BASE_NLON must be greater than zero",myname)
       return
    end if
!
    if(nchunk > 2) then
       ! check whether vector nlat_block equals vector nlon_block
       ! cannot permit differences in that case, since at some chunk boundaries base cells do not meet
       ! those in the neighbouring chunk (e.g. the boundary between chunk 2 and 3)
       do iblock = 1,nblock
          if(nlat_block(iblock) /= nlon_block(iblock)) then
             write(errstr,*) "in ",iblock,"'th refinement block: NLAT = ",nlat_block(iblock)," and NLON = ",&
                  nlon_block(iblock),"; in case of NCHUNK > 2, NLAT must equal NLON in all blocks"
             call add(errmsg,2,errstr,myname)
             return             
          end if
       end do ! iblock
    end if
!
    char_string = sval(inpar,'CHUNKS_INVGRID_CREF_METHOD',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"keyword 'CHUNKS_INVGRID_CREF_METHOD' is not present in parameter file",myname)
       return
    end if
    select case(trim(char_string))
    case('NONE')
       ! OK, no values from CHUNKS_INVGRID_CREF_PARAMETERS required, so do nothing
    case('EXPERIMENTAL_RANDOM_SUBDIVISION')
       ! check if there is at least one real value in vector CHUNKS_INVGRID_CREF_PARAMETERS and read this value
       cref_param => rvecp(inpar,'CHUNKS_INVGRID_CREF_PARAMETERS',2,iostat=ios)
       if(ios /= 0) then
          write(errstr,*) "could not read two real values of 'CHUNKS_INVGRID_CREF_PARAMETERS' from string '"//&
               trim(inpar.sval.'CHUNKS_INVGRID_CREF_PARAMETERS')//"'"
          call add(errmsg,2,errstr,myname)
          return
       end if
!##################################
! ADD YOUR REFINEMENT METHOD HERE
!##################################
    !case ('ADD_NEW_REFINEMENT_METHOD_HERE')
    !   ! If required, read in some numbers given by keyword 'CHUNKS_INVGRID_CREF_PARAMETERS'.
    !   ! Alternatively, you may need to introduce another keyword to read in some filenames
    !   ! or introduce some naming convention for files containing additional information.
    !   cref_param => rvecp(inpar,'CHUNKS_INVGRID_CREF_PARAMETERS',number_of_values_you_expect,iostat=ios)
    case default
       call add(errmsg,2,"value '"//trim(char_string)//"' of keyword 'CHUNKS_INVGRID_CREF_METHOD' is invalid; "//&
            "must be one of 'NONE', 'EXPERIMENTAL_RANDOM_SUBDIVISION'",myname)
       return
    end select
!
    char_string = sval(inpar,'VTK_PROJECTION',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"keyword 'VTK_PROJECTION' is not present in parameter file",myname)
       return
    end if
    select case(trim(char_string))
    case('GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT')
       ! OK, so do nothing
    case default
       call add(errmsg,2,"value '"//trim(char_string)//"' of keyword 'VTK_PROJECTION' is invalid; "//&
            "must be one of 'GLOBAL', 'LOCAL_CURV', 'LOCAL_FLAT', 'LOCAL_NORTH_CURV', 'LOCAL_NORTH_FLAT'",myname)
       return
    end select
!
    char_string = sval(inpar,'VTK_GEOMETRY_TYPE',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"keyword 'VTK_GEOMETRY_TYPE' is not present in parameter file",myname)
       return
    end if
    vtk_geometry_type_int = intGeometryTypeChunksInversionGrid(char_string)
    if(vtk_geometry_type_int < 0) then
       call add(errmsg,2,"vtk geometry type '"//trim(char_string)//"' is invalid",myname)
       return
    end if
!
    apply_vtk_scaling = lval(inpar,'SCALE_VTK_COORDS',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read logical value for 'SCALE_VTK_COORDS' from '"//&
            trim(inpar.sval.'SCALE_VTK_COORDS')//"'",myname)
       return
    end if
    if(apply_vtk_scaling) then
       vtk_factor = rval(inpar,'VTK_COORDS_SCALING_FACTOR',iostat=ios)
       if(ios /= 0) then
          call add(errmsg,2,"could not read real value for 'VTK_COORDS_SCALING_FACTOR' from '"//&
               trim(inpar.sval.'VTK_COORDS_SCALING_FACTOR')//"'",myname)
          return
       end if
    end if
  end subroutine readCheckParFileChunksInversionGrid
!--------------------------------------------------------------------
!> \brief construct a new chunks inversion grid
  subroutine constructNewChunksInversionGrid(this,path,inpar,lu,errmsg)
    type(chunks_inversion_grid) :: this
    character (len=*) :: path
    type (input_parameter) :: inpar
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=31) :: myname = 'constructNewChunksInversionGrid'
    character(len=400) :: errstr
    integer :: ios
    logical :: file_exists
    integer :: nblock,ilay,ibl,ix,iy,iz,icell,ichunk,ncell_base_1chunk,ishift
    integer, dimension(:), pointer :: nlay_block,nlat,nlon,ip
    double precision :: wx_half,wy_half,alpha_min,alpha_max,walpha
    real, dimension(:), pointer :: thickness,cref_param
    integer, dimension(:), allocatable :: ixmin_base_1chunk,ixmax_base_1chunk,iymin_base_1chunk,iymax_base_1chunk,&
         izmin_base_1chunk,izmax_base_1chunk,index_in_x,index_in_y
    double precision, dimension(:), allocatable :: c_add
    integer, dimension(:,:,:,:), allocatable :: icell_base_1chunk
    character (len=500) :: invgrid_filename,cref_method
!
    nullify(nlay_block,nlat,nlon,ip,thickness,cref_param)
    call addTrace(errmsg,myname)
    this%is_defined = .false.
!
    invgrid_filename = trim(path)//trim(inpar.sval.'CHUNKS_INVGRID_FILE')
    ! delete old binary inversion grid file, if existing
    inquire(file = invgrid_filename, exist = file_exists)
    if(file_exists) then
       open(unit=lu,file=invgrid_filename,form='UNFORMATTED',access='STREAM',status='OLD',action='READ',iostat=ios)
       close(unit=lu,status='DELETE')
    end if
!
    ! GET VALUES FROM PARAMETER FILE
!
    call readCheckParFileChunksInversionGrid(inpar,this%nchunk,this%clat,this%clon,this%wlat,this%wlon,this%rmax,&
         this%az_rot,nblock,nlay_block,thickness,nlat,nlon,cref_param,this%vtk_geometry_type_int,&
         this%apply_vtk_coords_scaling_factor,this%vtk_coords_scaling_factor,errmsg,myname)
    if(.level.errmsg == 2) goto 2
!
    cref_method = trim(inpar.sval.'CHUNKS_INVGRID_CREF_METHOD')
    this%vtk_projection = trim(inpar.sval.'VTK_PROJECTION')
    call computeTransformsChunksInversionGrid(this)
!
    ! NOW DEFINE z-COORDINATE DISTRIBUTION OF ALL INVERSION GRID BASE CELLS 
    ! since the inversion grid is defined in layers of base cells, the complete z coordinate 
    ! vector can be defined in advance now and referred to below
    ! Note that this%z will have DECREASING order (starting with Earth's radius = surface going down, decreasing radius).
    this%nz = sum(nlay_block)+1
    allocate(this%z(this%nz))
    this%z(1) = dble(this%rmax)
    iz = 1
    do ibl = 1,nblock
       do ilay = 1,nlay_block(ibl)
          iz = iz + 1
          this%z(iz) = this%z(iz-1) - dble(thickness(ibl))
          if(this%z(iz) < 0) then
             write(errstr,*) iz-1,"'th depth layer of base cells bottoms below the center of Earth"
             call add(errmsg,2,errstr,myname)
             goto 2
          end if
       end do ! ilay
    end do ! ibl
!
    ! NOW PREPARE x-COORDINATES AND y-COORDINATES IN A LOCAL FLAT CHUNK, NORMALIZED BY this%rmax, 
    ! i.e. the tangential plane touches a sphere of radius 1 at the north pole z=1
    if(this%nchunk == 1) then
       wx_half = tan(0.5d0*dble(this%wlat)*mc_deg2radd)
       wy_half = tan(0.5d0*dble(this%wlon)*mc_deg2radd)
    else
       ! for multiple chunks, only wlat=wlon=90.0 is permitted (checked above), so use the exact analytical values here
       wx_half = 1.d0
       wy_half = 1.d0
    end if
!
    ! start with just two coordinates that are used in every refienement block (outer boundaries of inversion grid chunk)
    ! inside the loop on refinement blocks (below), add the rest of the base cell boundary coordinates (if not yet 
    ! exising due to the lateral refinement if a previous block)
    this%nx = 2 
    this%ny = 2
    allocate(this%x(this%nx),this%y(this%ny))
    ! for all refinement blocks of base cells, use x(1),x(2),y(1),y(2) for the outer-most coordinates of a chunk 
    ! (all blocks definitley have those coordinates in common, this way it is also assured that the cells are 
    ! exactly aligned at the outer lateral boundaries of a chunk)
    this%x(1) = -wx_half; this%x(2) = wx_half
    this%y(1) = -wy_half; this%y(2) = wy_half
!
    ncell_base_1chunk = sum(nlon * nlat * nlay_block) ! number of base cells in a chunk (the same for all chunks)
    allocate(ixmin_base_1chunk(ncell_base_1chunk),ixmax_base_1chunk(ncell_base_1chunk),&
         iymin_base_1chunk(ncell_base_1chunk),iymax_base_1chunk(ncell_base_1chunk),&
         izmin_base_1chunk(ncell_base_1chunk),izmax_base_1chunk(ncell_base_1chunk))
    allocate(icell_base_1chunk(maxval(nlat),maxval(nlon),maxval(nlay_block),nblock))
    icell_base_1chunk = -1
!
    ! loop on all blocks, layers in a block, and lateral cells and define this%x,this%y as well as 
    ! temporary index arrays ixmin,ixmax,iymin,iymax,izmin,izmax, which define one chunk (can be 
    ! reused for base cells for all the chunks)
    icell = 0
    iz = 0
    do ibl = 1,nblock
       ! define x-coordinates NOT EQUIDISTANT on the tangential plane (as done in schunkInversionGrid), 
       ! BUT EQUIANGULARLY distributed
       walpha = dble(this%wlat)*mc_deg2radd
       alpha_min = -0.5d0*walpha
       ! create temporary array c_add containing all x-boundary values and "add" it to this%x
       ! by routine addCoordinatesChunksInversionGrid
       if(allocated(c_add)) deallocate(c_add)
       allocate(c_add(nlat(ibl)+1))
       c_add(1) = this%x(1) ! xmin of chunk -> coordinate should be located as existing at index 1
       do ix = 2,nlat(ibl)
          c_add(ix) = tan( alpha_min + ((ix-1)*walpha)/nlat(ibl) )
       end do ! ix
       c_add(nlat(ibl)+1) = this%x(2) ! xmax of chunk -> coordinate should be located as existing at index 2
       if(allocated(index_in_x)) deallocate(index_in_x)
       allocate(index_in_x(nlat(ibl)+1))
       call addCoordinatesChunksInversionGrid(this%x,this%nx,c_add,nlat(ibl)+1,index_in_x)
!
       ! define y-coordinates NOT EQUIDISTANT on the tangential plane (as done in schunkInversionGrid), 
       ! BUT EQUIANGULARLY distributed
       walpha = dble(this%wlon)*mc_deg2radd
       alpha_min = -0.5d0*walpha
       ! create temporary array c_add containing all x-boundary values and "add" it to this%x
       ! by routine addCoordinatesChunksInversionGrid
       if(allocated(c_add)) deallocate(c_add)
       allocate(c_add(nlon(ibl)+1))
       c_add(1) = this%y(1) ! ymin of chunk -> coordinate should be located as existing at index 1
       do iy = 2,nlon(ibl)
          c_add(iy) = tan( alpha_min + ((iy-1)*walpha)/nlon(ibl) )
       end do ! iy
       c_add(nlon(ibl)+1) = this%y(2) ! ymax of chunk -> coordinate should be located as existing at index 2
       if(allocated(index_in_y)) deallocate(index_in_y)
       allocate(index_in_y(nlon(ibl)+1))
       call addCoordinatesChunksInversionGrid(this%y,this%ny,c_add,nlon(ibl)+1,index_in_y)
!
       do ilay = 1,nlay_block(ibl)
          iz = iz + 1
          do iy = 1,nlon(ibl)
             do ix = 1,nlat(ibl)
                icell = icell+1
!
                ixmin_base_1chunk(icell) = index_in_x(ix)
                ixmax_base_1chunk(icell) = index_in_x(ix+1)
!
                iymin_base_1chunk(icell) = index_in_y(iy)
                iymax_base_1chunk(icell) = index_in_y(iy+1)
!
                izmin_base_1chunk(icell) = iz + 1 ! array this%z has decreasing order (radius: surface -> center)r
                izmax_base_1chunk(icell) = iz
!
                ! for convenience, remember this cell index in dependence of ix,iy,ilay,ibl
                ! (used below in findFaceNeighboursOfBaseCellsChunksInversionGrid)
                icell_base_1chunk(ix,iy,ilay,ibl) = icell
             end do ! ix
          end do ! iy
       end do ! ilay
    end do ! ibl
    if(allocated(c_add)) deallocate(c_add)
    if(allocated(index_in_x)) deallocate(index_in_x)
    if(allocated(index_in_y)) deallocate(index_in_y)
!
    ! NOW DEFINE THE BASE CELLS OF THE CHUNKS
! 
    this%ncell_base = this%nchunk*ncell_base_1chunk
    this%ncell_internal = this%ncell_base
    this%ncell = this%ncell_internal
!
    allocate(this%cell_ichunk(this%ncell_internal),this%icell_internal(this%ncell),&
         this%icell_external(this%ncell_internal),&
         this%cell_ixmin(this%ncell_internal),this%cell_ixmax(this%ncell_internal),&
         this%cell_iymin(this%ncell_internal),this%cell_iymax(this%ncell_internal),&
         this%cell_izmin(this%ncell_internal),this%cell_izmax(this%ncell_internal),&
         this%cell_parent(this%ncell_internal),this%cell_subcells(this%ncell_internal))
!
    ! parents (all base cells, i.e. cell_parents = -1
    this%cell_parent(1:this%ncell_internal) = -1
!
    ! at the moment there are no subcells, i.e. the particular pointers in array cell_subcells stay unassociated
    ! (as they are, initially)
!
    ! since there are only base cells (before doing the cell refinement), the mappings icell_internal and icell_external are identities
    this%icell_internal(1:this%ncell) = (/ (icell,icell=1,this%ncell) /)
    this%icell_external(1:this%ncell_internal) = (/ (icell,icell=1,this%ncell_internal) /)
!
    ! define the cell boundary index arrays cell_i(x/y/z)(min/max) and array cell_ichunk
    do ichunk = 1,this%nchunk
       ishift = (ichunk-1)*ncell_base_1chunk
       this%cell_ichunk(ishift+1:ishift+ncell_base_1chunk) = ichunk
       this%cell_ixmin(ishift+1:ishift+ncell_base_1chunk) = ixmin_base_1chunk(1:ncell_base_1chunk)
       this%cell_ixmax(ishift+1:ishift+ncell_base_1chunk) = ixmax_base_1chunk(1:ncell_base_1chunk)
       this%cell_iymin(ishift+1:ishift+ncell_base_1chunk) = iymin_base_1chunk(1:ncell_base_1chunk)
       this%cell_iymax(ishift+1:ishift+ncell_base_1chunk) = iymax_base_1chunk(1:ncell_base_1chunk)
       this%cell_izmin(ishift+1:ishift+ncell_base_1chunk) = izmin_base_1chunk(1:ncell_base_1chunk)
       this%cell_izmax(ishift+1:ishift+ncell_base_1chunk) = izmax_base_1chunk(1:ncell_base_1chunk)
    end do ! ichunk
!
    ! Define this%chunk_ncell_internal and the mapping this%chunk_icell_internal, defining for each chunk
    ! a vector of internal cell indices of all cells contained in that chunk.
    allocate(this%chunk_ncell_internal(this%nchunk),this%chunk_icell_internal(this%nchunk))
    ! since there are only base cells (before doing the cell refinement), every chunk contains the same number of cells, 
    ! namely ncell_base_1chunk and the index arrays this%chunk_icell_internal(:) all equal (1,..,ncell_base_1chunk)
    this%chunk_ncell_internal(:) = ncell_base_1chunk
    do ichunk = 1,this%nchunk
       allocate(ip(ncell_base_1chunk))
       ip = (/(icell,icell=1,ncell_base_1chunk)/)
       call associateVectorPointer(this%chunk_icell_internal(ichunk),ip)
       nullify(ip)
    end do ! ichunk
!
    ! NOW DO CELL REFINEMENT
    !   extend arrays icell_internal,x,y,z,cell_ichunk,chunk_icell_internal(:),cell_i(x,y,z)(min,max),cell_center,cell_radius
    !   extend and modify array icell_external (mark all base cells as NOT external when they are refined)
    !   modify array chunk_ncell_internal (increase the numbers of internal cells for each chunk by the refined subcells)
    call doCellRefinementChunksInversionGrid(this,cref_method,cref_param)
!
    ! NOW SORT THE x,y,z COORDINATES
    call sortCellBoundaryCoordinatesChunksInversionGrid(this,sum(nlay_block)+1)
!
    ! NOW DEFINE NEIGHBOUR ARRAY FOR ALL BASE CELLS
    call findFaceNeighboursOfBaseCellsChunksInversionGrid(this,nblock,nlay_block,nlon,nlat,ncell_base_1chunk,&
         icell_base_1chunk)
!
    ! NOW DEFINE NEIGHBOUR ARRAY FOR ALL REFINED CELLS
    call findFaceNeighboursOfRefinedCellsChunksInversionGrid(this)
!
    ! NOW DEFINE CELL CENTERS IN LOCAL FLAT COORDINATES
    allocate(this%cell_center(3,this%ncell_internal))
    ! cell centers are laterally defined w.r.t. angle!!, so apply atan and tan operations here
    do icell = 1,this%ncell_internal
       ! x coordinate of cell center in tangential plane
       alpha_min = atan(this%x(this%cell_ixmin(icell)))
       alpha_max = atan(this%x(this%cell_ixmax(icell)))
       this%cell_center(1,icell) = tan(0.5d0*(alpha_min+alpha_max))
       ! y coordinate of cell center in tangential plane
       alpha_min = atan(this%y(this%cell_iymin(icell)))
       alpha_max = atan(this%y(this%cell_iymax(icell)))
       this%cell_center(2,icell) = tan(0.5d0*(alpha_min+alpha_max))
       ! z coordinate of cell center (centered w.r.t. depth)
       this%cell_center(3,icell) = 0.5d0*(this%z(this%cell_izmin(icell)) + this%z(this%cell_izmax(icell)))
    end do ! icell
!
    ! NOW DEFINE CELL RADII
    call computeCellRadiiChunksInversionGrid(this)
!
    ! write binary inversion grid file
    ! since the old invgrid file was deleted above, set open_status = 'NEW' here
    call writeChunksInversionGrid(this,invgrid_filename,lu,'NEW',errmsg)
    if(.level.errmsg == 2) goto 2
!
    ! cleaning up and returning, if everything went alright
!
    this%is_defined = .true.
!
1   call dealloc(inpar)
    if(associated(thickness)) deallocate(thickness)
    if(associated(nlay_block)) deallocate(nlay_block)
    if(associated(nlat)) deallocate(nlat)
    if(associated(nlon)) deallocate(nlon)
    if(allocated(ixmin_base_1chunk)) deallocate(ixmin_base_1chunk)
    if(allocated(ixmax_base_1chunk)) deallocate(ixmax_base_1chunk)
    if(allocated(iymin_base_1chunk)) deallocate(iymin_base_1chunk)
    if(allocated(iymax_base_1chunk)) deallocate(iymax_base_1chunk)
    if(allocated(izmin_base_1chunk)) deallocate(izmin_base_1chunk)
    if(allocated(izmax_base_1chunk)) deallocate(izmax_base_1chunk)
    if(allocated(icell_base_1chunk)) deallocate(icell_base_1chunk)
    return
!
    ! destroy everything newly created before returning
2   call deallocateChunksInversionGrid(this)
    goto 1
  end subroutine constructNewChunksInversionGrid
!--------------------------------------------------------------------
  subroutine computeTransformsChunksInversionGrid(this)
    type (chunks_inversion_grid) :: this
    double precision :: ctheta,costheta,cosphi,cosazmth,sintheta,sinphi,sinazmth,x_shift,y_shift
    integer :: ichunk
!
    ! (allocate and) compute here:
    ! this%Mrot_azimuth
    ! this%Mrot_global2local
    ! this%Mrot_local2flat
    ! this%Mrot_global2flat
    ! this%vtk_flat_shift
!
    ! precomputation of some parameters
    ctheta = 90.d0-dble(this%clat)
    cosazmth = cos(dble(this%az_rot)*mc_deg2radd)
    sinazmth = sin(dble(this%az_rot)*mc_deg2radd)
    costheta = cos(ctheta*mc_deg2radd)
    cosphi = cos(dble(this%clon)*mc_deg2radd)
    sintheta = sin(ctheta*mc_deg2radd)
    sinphi = sin(dble(this%clon)*mc_deg2radd)
!
    this%Mrot_azimuth(1,1) = cosazmth
    this%Mrot_azimuth(1,2) = -sinazmth
    this%Mrot_azimuth(2,1) = sinazmth
    this%Mrot_azimuth(2,2) = cosazmth
!
    this%Mrot_global2local(1,1) = cosazmth*costheta*cosphi-sinazmth*sinphi
    this%Mrot_global2local(1,2) = cosazmth*costheta*sinphi+sinazmth*cosphi
    this%Mrot_global2local(1,3) = -cosazmth*sintheta
    this%Mrot_global2local(2,1) = -sinazmth*costheta*cosphi-cosazmth*sinphi
    this%Mrot_global2local(2,2) = -sinazmth*costheta*sinphi+cosazmth*cosphi
    this%Mrot_global2local(2,3) = sinazmth*sintheta
    this%Mrot_global2local(3,1) = sintheta*cosphi
    this%Mrot_global2local(3,2) = sintheta*sinphi
    this%Mrot_global2local(3,3) = costheta
!
    allocate(this%Mrot_local2flat(3,3,this%nchunk))
    this%Mrot_local2flat(:,:,:) = 0.d0
    do ichunk = 1,this%nchunk
       select case(ichunk)
       case(1) ! chunk 1 is already in the correct position after rotation global2local, so local2flat is identity
          this%Mrot_local2flat(1,1,ichunk) = 1.d0
          this%Mrot_local2flat(2,2,ichunk) = 1.d0
          this%Mrot_local2flat(3,3,ichunk) = 1.d0
       case(2) ! chunk 2 touches chunk 1 on its "left" side (looking from above onto unrotated chunk 1)
          this%Mrot_local2flat(1,1,ichunk) = 1.d0
          this%Mrot_local2flat(2,3,ichunk) = 1.d0
          this%Mrot_local2flat(3,2,ichunk) = -1.d0
       case(3) ! chunk 3 touches chunk 1 on its "top" side (looking from above onto unrotated chunk 1)
          this%Mrot_local2flat(1,3,ichunk) = 1.d0
          this%Mrot_local2flat(2,2,ichunk) = 1.d0
          this%Mrot_local2flat(3,1,ichunk) = -1.d0
       case(4) ! chunk 4 lies opposite of chunk 1
          this%Mrot_local2flat(1,1,ichunk) = 1.d0
          this%Mrot_local2flat(2,2,ichunk) = -1.d0
          this%Mrot_local2flat(3,3,ichunk) = -1.d0
       case(5) ! chunk 5 lies opposite of chunk 2
          this%Mrot_local2flat(1,1,ichunk) = 1.d0
          this%Mrot_local2flat(2,3,ichunk) = -1.d0
          this%Mrot_local2flat(3,2,ichunk) = 1.d0
       case(6) ! chunk 6 lies opposite of chunk 3
          this%Mrot_local2flat(1,3,ichunk) = -1.d0
          this%Mrot_local2flat(2,2,ichunk) = 1.d0
          this%Mrot_local2flat(3,1,ichunk) = 1.d0
       end select ! ichunk
    end do ! ichunk
!
    allocate(this%Mrot_global2flat(3,3,this%nchunk))
    do ichunk = 1,this%nchunk
       this%Mrot_global2flat(:,:,ichunk) = matmul(this%Mrot_local2flat(:,:,ichunk),this%Mrot_global2local(:,:))
    end do ! ichunk
!
    x_shift = dble(this%rmax)*dble(this%wlat)*mc_deg2radd
    y_shift = dble(this%rmax)*dble(this%wlon)*mc_deg2radd
    allocate(this%vtk_flat_shift(2,this%nchunk))
    do ichunk = 1,this%nchunk
       select case(ichunk)
       case(1)
          this%vtk_flat_shift(1,ichunk) = 0.d0
          this%vtk_flat_shift(2,ichunk) = 0.d0
       case(2)
          this%vtk_flat_shift(1,ichunk) = 0.d0
          this%vtk_flat_shift(2,ichunk) = -y_shift
       case(3)
          this%vtk_flat_shift(1,ichunk) = -x_shift
          this%vtk_flat_shift(2,ichunk) = 0.d0
       case(4)
          this%vtk_flat_shift(1,ichunk) = 0.d0
          this%vtk_flat_shift(2,ichunk) = y_shift + y_shift
       case(5)
          this%vtk_flat_shift(1,ichunk) = 0.d0
          this%vtk_flat_shift(2,ichunk) = y_shift
       case(6)
          this%vtk_flat_shift(1,ichunk) = x_shift
          this%vtk_flat_shift(2,ichunk) = 0.d0
       end select
    end do ! ichunk
  end subroutine computeTransformsChunksInversionGrid
!--------------------------------------------------------------------
!> \brief do (recusrive) cell refinement of base cells
  subroutine doCellRefinementChunksInversionGrid(this,cref_method,cref_param)
    type (chunks_inversion_grid) :: this
    character(len=*) :: cref_method
    real, dimension(:), pointer :: cref_param
    !type (error_message) :: errmsg
    ! local
    !character(len=35) :: myname = 'doCellRefinementChunksInversionGrid'
!
    ! in readCheckParFileChunksInversionGrid it was already checked that the value of cref_method 
    ! is valid and that there exist according values in cref_param. DO NOT CHECK THIS AGAIN HERE!
!
    select case(cref_method)
    case ('NONE')
       return
    case ('EXPERIMENTAL_RANDOM_SUBDIVISION')
       call srand(int(cref_param(1)))
       call doExpRandCellRefinementChunksInversionGrid(this,int(cref_param(2)))
!##################################
! ADD YOUR REFINEMENT METHOD HERE
!##################################
    !case ('ADD_NEW_REFINEMENT_METHOD_HERE')
    !   call requiredSubroutinesForNewRefinementMethod(this)
    end select
!
  end subroutine doCellRefinementChunksInversionGrid
!--------------------------------------------------------------------
!> \brief do cell refinement for Your New Refinement Method
  !subroutine requiredSubroutinesForNewRefinementMethod(this)
    !type (chunks_inversion_grid) :: this
!
!##################################
! ADD YOUR REFINEMENT METHOD HERE
!##################################
!
    ! THERE IS NO MULTIPLE REFINEMENT, i.e. MAKE SURE YOU ONLY HAVE ONE REFINEMENT HIERARCHY: 
    ! BASE CELLS CAN HAVE SUBCELLS, SUBCELLS  C A N N O T  HAVE SUBCELLS (otherwise, this module might become inconsistent)
    ! This means, that a refined base cell becomes internal only (no external cell anymore) and all its subcells are external.
    ! You should NOT implement a recursive subroutine here, for larger grids you may easily have stack-overflow. 
    ! If you want to have a recursive scheme, you should implmement it iteratively (compare 
    ! doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid contained in doExpRandCellRefinementChunksInversionGrid)
    ! 
!
    ! FOR EACH BASE CELL (index range 1,..,this%ncell_base) CHOOSE WHETHER TO REFINE IT OR NOT AND (if yes), DEFINE A SET OF SUBCELLS.
    ! FOR EACH BASE CELL YOU NEED TO:
    ! - allocate and set its pointer in this%cell_subcells (if base cell is refined, containing internal cell indices of the new subcells), leave unassociated otherwise
    ! - if refined, remove base cell index from this%icell_internal (you shoud overwrite it by internal index of first subcell)
    ! - set this%icell_external(icell_base) = -1 (base cell is not external anymore)
!
    ! EACH NEW SUBCELL IS A NEW EXTERNAL CELL OF THIS INVERSION GRID.
    ! FOR EACH SUBCELL YOU NEED TO:
    ! - increase counter this%ncell_internal -> internal cell indices of new subcells are just successors of this%ncell_internal, i.e. append them to the existing cell indices
    ! - increase counter this%ncell (number of external cells, be aware that this%ncell additionally decreases by 1 due to the base cell not being external anymore)
    ! - if necessary, add new cell boundary values to this%x,this%y,this%z (MAKE SURE YOU DO NOT ADD DUPLICATE VALUES!!, better use routine addCoordinatesChunksInversionGrid)
    ! - append a value to this%cell_ixmin,this%cell_ixmax,this%cell_iymin,... (subcell boundaries, indices in this%x,this%y,this%z)
    ! - append a value to this%cell_ichunk (index of chunk in which the cell is located: same index as that of base cell)
    ! - for its chunk index 'ichunk', increase this%chunk_ncell_internal(ichunk) and append the internal cell index to this%chunk_icell_internal(ichunk)
    ! - append a value to this%cell_parent (base cell index in which the subcell is contained)
    ! - add a value to this%icell_internal (first subcell of the base cell should overwrite the base cell's entry in this%icell_internal, since the base cell is not external anymore
    ! - add a value to this%icell_external (defining the external cell index of the subcell, make sure you account for the index shift due to the base cell not being external anymore, or simply define as inverse mapping of this%icell_internal, compare doExpRandCellRefinementChunksInversionGrid)
!
  !end subroutine requiredSubroutinesForNewRefinementMethod
!--------------------------------------------------------------------
!> \brief do cell refinement for experimental random method
  subroutine doExpRandCellRefinementChunksInversionGrid(this,nmax_subcells)
    type (chunks_inversion_grid) :: this
    integer :: nmax_subcells
    ! local
    integer :: isub,isub_start,isub_end,nsub_added,nsub,nsub_external
    integer :: ic,ix,nx,iy,ny,iz,nz,nxnynz
    integer :: icell_base,icell_base_external,ichunk,j
    logical,  dimension(:), pointer :: isub_is_external
    double precision, dimension(:), pointer :: dpp,xmin_sub,xmax_sub,ymin_sub,ymax_sub,zmin_sub,zmax_sub
    double precision, dimension(:), allocatable :: x,y,z
    double precision :: xmin,xmax,ymin,ymax,zmin,zmax,delta
    integer, dimension(:), pointer :: ip
    integer, dimension(:), allocatable :: int_array
!
    nullify(isub_is_external,xmin_sub,xmax_sub,ymin_sub,ymax_sub,zmin_sub,zmax_sub)
!
    ! loop on all base cells and do refinement
    do icell_base = 1,this%ncell_base
!
       !write(*,*) "BASE CELL ",icell_base," out of ",this%ncell_base
!
       if(associated(xmin_sub)) deallocate(xmin_sub)
       if(associated(xmax_sub)) deallocate(xmax_sub)
       if(associated(ymin_sub)) deallocate(ymin_sub)
       if(associated(ymax_sub)) deallocate(ymax_sub)
       if(associated(zmin_sub)) deallocate(zmin_sub)
       if(associated(zmax_sub)) deallocate(zmax_sub)
       nsub = 0
       call doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid()
       !! RECURSION DOES NOT WORK PROPERLY (see below)
       !! call doRecursiveExpRandCellRefinementOfBaseCellChunksInversionGrid(&
       !!      this%x(this%cell_ixmin(icell_base)),this%x(this%cell_ixmax(icell_base)),xmin_sub,xmax_sub,&
       !!      this%y(this%cell_ixmin(icell_base)),this%y(this%cell_ixmax(icell_base)),ymin_sub,ymax_sub,&
       !!      this%z(this%cell_ixmin(icell_base)),this%z(this%cell_ixmax(icell_base)),zmin_sub,zmax_sub,&
       !!      0,parent)
!
       ! all returning pointers x**_sub,y**_sub,z**_sub are assumed to have the same size on return
       ! (or be not associated).
       ! the size of these arrays is stored in variable nsub
!
       if(.not.associated(xmin_sub)) cycle
!
       ! reallocate arrays in 'this' to add all subcells of current base cell
       this%cell_ixmin => reallocate(this%cell_ixmin,this%ncell_internal+nsub)
       this%cell_ixmax => reallocate(this%cell_ixmax,this%ncell_internal+nsub)
       this%cell_iymin => reallocate(this%cell_iymin,this%ncell_internal+nsub)
       this%cell_iymax => reallocate(this%cell_iymax,this%ncell_internal+nsub)
       this%cell_izmin => reallocate(this%cell_izmin,this%ncell_internal+nsub)
       this%cell_izmax => reallocate(this%cell_izmax,this%ncell_internal+nsub)
       this%cell_ichunk => reallocate(this%cell_ichunk,this%ncell_internal+nsub)
       this%cell_parent => reallocate(this%cell_parent,this%ncell_internal+nsub)
       this%cell_subcells => extendArrayVectorPointer(this%cell_subcells,this%ncell_internal+nsub)
!
       allocate(int_array(nsub))
       call addCoordinatesChunksInversionGrid(this%x,this%nx,xmin_sub,nsub,int_array)
       this%cell_ixmin(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       call addCoordinatesChunksInversionGrid(this%x,this%nx,xmax_sub,nsub,int_array)
       this%cell_ixmax(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       call addCoordinatesChunksInversionGrid(this%y,this%ny,ymin_sub,nsub,int_array)
       this%cell_iymin(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       call addCoordinatesChunksInversionGrid(this%y,this%ny,ymax_sub,nsub,int_array)
       this%cell_iymax(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       call addCoordinatesChunksInversionGrid(this%z,this%nz,zmin_sub,nsub,int_array)
       this%cell_izmin(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       call addCoordinatesChunksInversionGrid(this%z,this%nz,zmax_sub,nsub,int_array)
       this%cell_izmax(this%ncell_internal+1:this%ncell_internal+nsub) = int_array
       deallocate(int_array)
!
       ! update this%cell_parent
       this%cell_parent(this%ncell_internal+1:this%ncell_internal+nsub) = this%icell_internal(icell_base)
!
       ! update this%cell_subcells
       allocate(ip(nsub))
       ip = (/ (j+this%ncell_internal,j=1,nsub) /)
       call associateVectorPointer(this%cell_subcells(icell_base),ip)
       nullify(ip)
!
       ichunk = this%cell_ichunk(icell_base)
       this%cell_ichunk(this%ncell_internal+1:this%ncell_internal+nsub) = ichunk
       ip => getVectorPointer(this%chunk_icell_internal(ichunk))
       ip => reallocate(ip,this%chunk_ncell_internal(ichunk)+nsub)
       ip(this%chunk_ncell_internal(ichunk)+1:this%chunk_ncell_internal(ichunk)+nsub) = &
            (/(j,j=this%ncell_internal+1,this%ncell_internal+nsub)/)
       call associateVectorPointer(this%chunk_icell_internal(ichunk),ip)
       nullify(ip)
       this%chunk_ncell_internal(ichunk) = this%chunk_ncell_internal(ichunk) + nsub
!
!
! THE FOLLOWING CODE REPLACES THE BASE CELL BY THE FIRST NEW SUBCELL (IN TERMS OF EXTERNAL CELL INDEX)
! AND APPENDS ALL OTHER SUBCELLS TO THE MAXIMUM EXTERNAL CELL INDEX -> YIELDS "SCATTERED" INDEX DISTRIBUTION
       ! ! there were ONLY external subcells returned by routine
       ! ! doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid,
       ! ! so add new cells as external cells and mark this base cell as INTERNAL only
       ! allocate(int_array(nsub))
       ! int_array = (/(j+this%ncell_internal,j=1,nsub)/)
       ! this%icell_internal => reallocate(this%icell_internal,this%ncell+nsub-1)
       ! this%icell_external => reallocate(this%icell_external,this%ncell_internal+nsub)
       ! this%icell_external(this%ncell_internal+1:this%ncell_internal+nsub) = -1
       ! ! in terms of external cells, replace the current base cell (formerly external) by the first subcell
       ! ! (that is external) instead of removing the external cell (this way, the rest of the array 
       ! ! this%icell_external does not need to be shifted)
       ! this%icell_external(int_array(1)) = this%icell_external(icell_base)
       ! this%icell_internal(this%icell_external(int_array(1))) = int_array(1)
       ! this%icell_external(icell_base) = -1
       ! ! define external cells for the rest of the new external cells
       ! this%icell_internal(this%ncell+1:this%ncell+nsub-1) = int_array(2:nsub)
       ! this%icell_external(int_array(2:nsub)) = (/(j,j=this%ncell+1,this%ncell+nsub-1)/)
       ! deallocate(int_array)
!
! NICER ALTERNATIVE TO ABOVE COMMENTED CODE SECTION (YIELDS NICER INVERSION GRID EXTERNAL CELL INDEX DISTRIBUTION):
! SHIFT THE EXTERNAL CELL INDEX AND "SQUEEZE" THE NEW CELLS INTO THE INDEX POSITION WHERE THE BASE CELL WAS
       ! There were ONLY external subcells (nsub many) returned by routine
       ! doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid.
       ! This base cell, however, is NOT internal anymore (therfore, -1 in next line)
       this%icell_internal => reallocate(this%icell_internal,this%ncell+nsub-1)
       this%icell_external => reallocate(this%icell_external,this%ncell_internal+nsub)
       this%icell_external(this%ncell_internal+1:this%ncell_internal+nsub) = -1 ! safety: initialize to -1
       icell_base_external = this%icell_external(icell_base)
       ! shift existing external indices (only those "behind" icell_base_external, if there are any)
       if(icell_base_external < this%ncell) then
          this%icell_internal(icell_base_external+nsub:this%ncell+nsub-1) = &
               this%icell_internal(icell_base_external+1:this%ncell)
          this%icell_external(this%icell_internal(icell_base_external+nsub:this%ncell+nsub-1)) = &
               (/ (j,j=icell_base_external+nsub,this%ncell+nsub-1) /)
       end if
       ! Insert the new external subcells into the range of existing external cell indices, 
       ! starting with the index that the base cell had so far (icell_base_external), i.e. replacing this
       ! entry by nsub ones.
       this%icell_internal(icell_base_external:icell_base_external+nsub-1) = &
            (/ (j,j=this%ncell_internal+1,this%ncell_internal+nsub) /)
       this%icell_external(this%ncell_internal+1:this%ncell_internal+nsub) = &
            (/ (j,j=icell_base_external,icell_base_external+nsub-1) /)
       ! Mark base cell as internal only
       this%icell_external(icell_base) = -1
!
!
       ! eventually, increase the counters of internal and external cells
       this%ncell_internal = this%ncell_internal + nsub
       this%ncell = this%ncell + nsub - 1 ! -1 here, because base cell is no external cell anymore
!
       this%contains_refined_cells = .true.
    end do  ! icell_base
!
  contains
!
    subroutine doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid()
      ! do random coin flip
!       if(rand() < 0.5) then
!          ! do not subdivide further
! write(*,*) "'NO' COIN FLIP => NO REFINING OF BASE CELL"
!          return
!       end if

!       ! yes, subdivide cell
! write(*,*) "'YES' COIN FLIP => REFINING BASE CELL"
!
      ! nsub should have been initiated to 0 before calling this subroutine
      !nsub = 0
!
      ! SUBDIVIDE BASE CELL: ARTIFICIALLY TREAT BASE CELL AS ONE NEW SUBCELL THAT WAS ADDED AND HAS TO 
      ! BE REFINED IN THE FIRST EXECUTION OF THE FOLLOWING LOOP (will only loop on one cell, namely the base cell)
      ! (by this mechanism, the execution of the following loop (designed for subcells) works also for the base cell)
      nsub_added = 1
!
100   isub_start = nsub-nsub_added+1
      isub_end = nsub
      nsub_added = 0
      do isub = isub_start,isub_end
!
         ! first check, if this cell is do be subdivided AT ALL (if not, the the following operations do not need
         ! to be done)
         nx = int(rand()*3)+1
         ny = int(rand()*3)+1
         nz = int(rand()*3)+1
         nxnynz = nx*ny*nz
         if(nxnynz <= 1 .or. nsub+nxnynz>=nmax_subcells) then
            !write(*,*) "NO REFINING OF CELL",isub
            ! memorize that this is an external cell (reduce arrays xmin_sub,xmax_sub,... below in order to 
            ! only contain external cells)
            if(isub>0) then 
               ! do not return base cell as one of the external subcells, i.e. in case of isub==0, there is
               ! no need to do anything
               isub_is_external(isub) = .true.
            end if
            cycle
         end if
         !write(*,*) "YES REFINING OF CELL",isub," BY CELLS ",nsub+1," ... ",nsub+nxnynz
!
         ! IF THE CODE COMES HERE, THIS (SUB)CELL IS TO BE SUBDIVIDED
!
         ! SUBDIVIDE THE isub'th SUBCELL: ALL SUBCELLS CREATED NOW WILL BE APPENDED TO TOTAL LIST OF SUBCELLS
!
         if(isub==0) then
            ! this is the base cell
            xmin = this%x(this%cell_ixmin(icell_base))
            xmax = this%x(this%cell_ixmax(icell_base))
            ymin = this%y(this%cell_iymin(icell_base))
            ymax = this%y(this%cell_iymax(icell_base))
            zmin = this%z(this%cell_izmin(icell_base))
            zmax = this%z(this%cell_izmax(icell_base))
         else
            ! this is a subcell of the base cell that was created before in this subroutine
            xmin = xmin_sub(isub)
            xmax = xmax_sub(isub)
            ymin = ymin_sub(isub)
            ymax = ymax_sub(isub)
            zmin = zmin_sub(isub)
            zmax = zmax_sub(isub)
          end if
!
         allocate(x(nx+1))
         x(1) = xmin
         x(nx+1) = xmax
         delta = (xmax-xmin)/dble(nx)
         do ix = 2,nx
            x(ix) = x(ix-1)+delta
         end do ! ix
!
         allocate(y(ny+1))
         y(1) = ymin
         y(ny+1) = ymax
         delta = (ymax-ymin)/dble(ny)
         do iy = 2,ny
            y(iy) = y(iy-1)+delta
         end do ! iy
!
         allocate(z(nz+1))
         z(1) = zmin
         z(nz+1) = zmax
         delta = (zmax-zmin)/dble(nz)
         do iz = 2,nz
            z(iz) = z(iz-1)+delta
         end do ! iz
!
         ! subdivide this cell into nx*ny*nz subcells, extend arrays that memorize subcell boundaries
         xmin_sub => reallocate(xmin_sub,nsub+nxnynz)
         xmax_sub => reallocate(xmax_sub,nsub+nxnynz)
         ymin_sub => reallocate(ymin_sub,nsub+nxnynz)
         ymax_sub => reallocate(ymax_sub,nsub+nxnynz)
         zmin_sub => reallocate(zmin_sub,nsub+nxnynz)
         zmax_sub => reallocate(zmax_sub,nsub+nxnynz)
         isub_is_external => reallocate(isub_is_external,nsub+nxnynz)
         isub_is_external(nsub+1:nsub+nxnynz) = .false. ! initiate to .false. (is set to .true. above if cell is external)
         ic = 0
         do iz = 1,nz
            do iy = 1,ny
               do ix = 1,nx
                  ic = ic + 1
                  xmin_sub(nsub+ic) = x(ix)
                  xmax_sub(nsub+ic) = x(ix+1)
                  ymin_sub(nsub+ic) = y(iy)
                  ymax_sub(nsub+ic) = y(iy+1)
                  zmin_sub(nsub+ic) = z(iz)
                  zmax_sub(nsub+ic) = z(iz+1)
               end do ! ix
            end do ! iy
         end do ! iz
!
         deallocate(x,y,z)
!
         ! increase counters of total number of subscells
         nsub = nsub+nxnynz
         ! increase counters of subscells added in this loop
         nsub_added = nsub_added+nxnynz
      end do ! isub
!
      ! if there were any subdivisions done inside the loop, invoke "recursion" by 
      ! going back and iterating over those subdivided cells, checking if themselves need to be subdivided
      if(nsub_added >0) then
         !write(*,*) "RE-EXECUTE THE LOOP"
         goto 100
      end if
!
      if(nsub>0) then
         nsub_external = count(isub_is_external)
         ! reduce array xmin_sub,xmax_sub,... in order to only contain external cells
         if(nsub_external < nsub) then
            ! some of the nsub subcells are NOT external (i.e. these cells were themselves subdivided)
            ! remove those from the return arrays
!
            ! xmin
            allocate(dpp(nsub_external))
            dpp = pack(xmin_sub,isub_is_external)
            deallocate(xmin_sub)
            xmin_sub => dpp
            nullify(dpp)
            ! xmax
            allocate(dpp(nsub_external))
            dpp = pack(xmax_sub,isub_is_external)
            deallocate(xmax_sub)
            xmax_sub => dpp
            nullify(dpp)
            ! ymin
            allocate(dpp(nsub_external))
            dpp = pack(ymin_sub,isub_is_external)
            deallocate(ymin_sub)
            ymin_sub => dpp
            nullify(dpp)
            ! ymax
            allocate(dpp(nsub_external))
            dpp = pack(ymax_sub,isub_is_external)
            deallocate(ymax_sub)
            ymax_sub => dpp
            nullify(dpp)
            ! zmin
            allocate(dpp(nsub_external))
            dpp = pack(zmin_sub,isub_is_external)
            deallocate(zmin_sub)
            zmin_sub => dpp
            nullify(dpp)
            ! zmax
            allocate(dpp(nsub_external))
            dpp = pack(zmax_sub,isub_is_external)
            deallocate(zmax_sub)
            zmax_sub => dpp
            nullify(dpp)
!
            nsub = nsub_external
         else ! nsub_external < nsub
            ! no need to do anything, as all nsub subcells are external.
            ! even if nsub is still 0, nsub_external must be 0, too, and arrays xmin_sub,... are not associated. 
            ! So everything alright.
         end if ! nsub_external < nsub
      end if ! nsub>0
      if(associated(isub_is_external)) deallocate(isub_is_external)
    end subroutine doIterativeExpRandCellRefinementOfBaseCellChunksInversionGrid

!#################################################################################################
! THE FOLLOWING RECURSIVE SUBROUTINE CAUSES SEGMENTATION FAULTS THAT I CANNOT TRACE BACK (expecially
! for large inversion grids/a lot of refinement recursions, probably stack overflow?!)
! ALTERNATIVELY, I IMPLEMENTED THE ABOVE ITERATIVE SUBROUTINE 
!#################################################################################################
!     recursive subroutine doRecursiveExpRandCellRefinementOfBaseCellChunksInversionGrid(&
!          xmin_base,xmax_base,xmin_sub,xmax_sub,&
!          ymin_base,ymax_base,ymin_sub,ymax_sub,&
!          zmin_base,zmax_base,zmin_sub,zmax_sub,base_index,parent_sub)
!       double precision, intent(in) :: xmin_base,xmax_base,ymin_base,ymax_base,zmin_base,zmax_base
!       double precision, dimension(:), pointer, intent(inout) :: xmin_sub,xmax_sub,ymin_sub,ymax_sub,zmin_sub,zmax_sub
!       integer, intent(in) :: base_index
!       integer, dimension(:), pointer, intent(inout) :: parent_sub
!       ! local
!       integer :: nsub,ic,ix,nx,iy,ny,iz,nz,nxnynz
!       double precision :: delta
!       double precision, dimension(:), allocatable :: x,y,z
! !
!       ! do random coin flip
!       if(rand() < 0.5) then
! !write(*,*) "YES COIN FLIP => REFINING CELL",base_index
!          ! yes, subdivide cell
!          if(associated(parent_sub)) then
!             !write(*,*) "associated(parent_sub)"
!             nsub = size(parent_sub)
!          else
!             !write(*,*) "NOT associated(parent_sub)"
!             nsub = 0
!          end if
! !write(*,*) "base_index, size(parent_sub) = ",base_index, nsub
! !
!          ! first check, if this base cell already contains enough subcells (termination criterion of recursion)
!          if(nsub >= 81) then
! !write(*,*) "nsub >= 81, so returning"
! write(*,*) "NO REFINING OF CELL",base_index
!             return
!          else
!          end if
! !
!          nx = int(rand()*3)+1
!          ny = int(rand()*3)+1
!          nz = int(rand()*3)+1
!          nxnynz = nx*ny*nz
!          if(nxnynz <= 1) then
! write(*,*) "NO REFINING OF CELL",base_index
!             return
!          end if
! write(*,*) "YES REFINING OF CELL",base_index," BY CELLS ",nsub+1," ... ",nsub+nxnynz
! !
!          allocate(x(nx+1))
!          x(1) = xmin_base
!          x(nx+1) = xmax_base
!          delta = (xmax_base-xmin_base)/dble(nx)
!          do ix = 2,nx
!             x(ix) = x(ix-1)+delta
!          end do ! ix
! !
!          allocate(y(ny+1))
!          y(1) = ymin_base
!          y(ny+1) = ymax_base
!          delta = (ymax_base-ymin_base)/dble(ny)
!          do iy = 2,ny
!             y(iy) = y(iy-1)+delta
!          end do ! iy
! !
!          allocate(z(nz+1))
!          z(1) = zmin_base
!          z(nz+1) = zmax_base
!          delta = (zmax_base-zmin_base)/dble(nz)
!          do iz = 2,nz
!             z(iz) = z(iz-1)+delta
!          end do ! iz
! !
!          ! subdivide this cell into nx*ny*nz subcells, extend outgoing arrays
!          xmin_sub => reallocate(xmin_sub,nsub+nxnynz)
!          xmax_sub => reallocate(xmax_sub,nsub+nxnynz)
!          ymin_sub => reallocate(ymin_sub,nsub+nxnynz)
!          ymax_sub => reallocate(ymax_sub,nsub+nxnynz)
!          zmin_sub => reallocate(zmin_sub,nsub+nxnynz)
!          zmax_sub => reallocate(zmax_sub,nsub+nxnynz)
!          ic = 0
!          do iz = 1,nz
!             do iy = 1,ny
!                do ix = 1,nx
!                   ic = ic + 1
!                   xmin_sub(nsub+ic) = x(ix)
!                   xmax_sub(nsub+ic) = x(ix+1)
!                   ymin_sub(nsub+ic) = y(iy)
!                   ymax_sub(nsub+ic) = y(iy+1)
!                   zmin_sub(nsub+ic) = z(iz)
!                   zmax_sub(nsub+ic) = z(iz+1)
!                end do ! ix
!             end do ! iy
!          end do ! iz
! !
!          deallocate(x,y,z)
! !
!          ! extend the array parents
!          parent_sub => reallocate(parent_sub,nsub+nxnynz)
!          parent_sub(nsub+1:nsub+nxnynz) = base_index
! !
!          ! invoke recursion
!          ic = 0
!          do iz = 1,nz
!             do iy = 1,ny
!                do ix = 1,nx
!                   ic = ic + 1
!                   call doRecursiveExpRandCellRefinementOfBaseCellChunksInversionGrid(&
!                        xmin_sub(nsub+ic),xmax_sub(nsub+ic),xmin_sub,xmax_sub,&
!                        ymin_sub(nsub+ic),ymax_sub(nsub+ic),ymin_sub,ymax_sub,&
!                        zmin_sub(nsub+ic),zmax_sub(nsub+ic),zmin_sub,zmax_sub,nsub+ic,parent_sub)
!                end do ! ix
!             end do ! iy
!          end do ! iz
! !         
!       else
!          ! do not subdivide further
! write(*,*) "NO REFINING OF CELL",base_index
!          return
!       end if
! !
!     end subroutine doRecursiveExpRandCellRefinementOfBaseCellChunksInversionGrid
!
  end subroutine doExpRandCellRefinementChunksInversionGrid
!--------------------------------------------------------------------
!> \brief add new coordinates in new_c to vector c, but detect existing ones (create NO duplicates!)
!! \param c pointer to array of existing coordinates (e.g. this%x,this%y,this%z): WILL BE EXTENDED IF NECESSARY
!! \param nc size of c (e.g. this%nx,this%ny,this%nz): WILL BE INCREASED IF NECESSARY
!! \param new_c array of new coordinates that should be added to c: COORDINATES ALREADY EXISTING IN c WILL NOT BE ADDED TO c!
!! \param nnew_c size of array new_c
!! \param index_new_c has size nnew_c; on return, index_new_c(i) equals the position of new_c(i) in array c, i.e. c(index_new_c(i)) = new_c(i) (for previously existing coordinates, index_new_c(i) <= nc, nc as on enter)
  subroutine addCoordinatesChunksInversionGrid(c,nc,new_c,nnew_c,&
            index_new_c)
    double precision, dimension(:), pointer, intent(inout) :: c
    integer, intent(inout) :: nc
    integer, intent(in) :: nnew_c
    double precision, dimension(nnew_c), intent(in) :: new_c
    integer, dimension(nnew_c), intent(out) :: index_new_c
    ! local
    integer :: nadd_c,inew_c,ic
    double precision, dimension(nnew_c) :: add_c
    logical :: coordinate_existed
!
    nadd_c = 0
!
    do inew_c = 1,nnew_c
!
       ! check existing coordinates, if new_c(inew_c) is among them
       coordinate_existed = .false.
       do ic = 1,nc
          if(c(ic)==new_c(inew_c)) then
             index_new_c(inew_c) = ic
             coordinate_existed = .true.
             exit
          end if
       end do ! ic
       if(coordinate_existed) cycle
!
       ! additionally, loop over all coordinates newly added so far
       do ic=1,nadd_c ! if nadd_c == 0, this loop does not do anything
          if(add_c(ic)==new_c(inew_c)) then
             index_new_c(inew_c) = nc+ic ! the new coordinates will be appended, so its index ic is shifted by nc
             coordinate_existed = .true.
             exit
          end if
       end do ! ic
       if(coordinate_existed) cycle
!
       ! if new_c(inew_c) is not among the existing coordinates nor among newly added ones, remember to add it
       nadd_c = nadd_c + 1
       add_c(nadd_c) = new_c(inew_c)
       index_new_c(inew_c) = nc+nadd_c
    end do ! inew_c
!
    ! if there were new coordinates found, that were not already contained
    ! in vector c, then extend c and append those new  coordinates; also update nc
    if(nadd_c > 0) then
       c => reallocate(c,nc+nadd_c)
       c(nc+1:nc+nadd_c) = add_c(1:nadd_c)
       nc = nc+nadd_c
    end if
  end subroutine addCoordinatesChunksInversionGrid
!--------------------------------------------------------------------
!> \brief define the sorting mappings of x,y,z coordinates of cell boundaries
  subroutine sortCellBoundaryCoordinatesChunksInversionGrid(this,nz_base)
    type (chunks_inversion_grid) :: this
    integer :: nz_base
    ! local
    integer :: j_sorted,ix,iy,iz
    double precision :: x_current,y_current,z_current
!
    allocate(this%sorted2ix(this%nx))
!
    ! the order of the boundary coordinates of the base cells are known to some extend: 
    !  - first entry in this%x is most left coordinate of the chunk (xmin chunk boundary)
    !  - second entry in this%x is most "right" x coordinate of the chunk (xmax chunk boundary)
    !  - third to (nx_base)'th entries are sorted in blocks (according to refinement blocks), i.e. an insertion sort is not the worst thing here
    this%sorted2ix(1) = 1
    this%sorted2ix(this%nx) = 2
    ! run insertion sort algorithm on the rest of the points
    do ix = 3,this%nx
       x_current = this%x(ix)
!
       ! we assume at this point that sorted2ix(1:ix-2) refers to sorted coordinates in this%x
       ! (additionally, the last index is known, which is 2, which is not of importance here)
       ! at first, put the current coordinate to the next position, i.e. to index ix-1
       ! thereafter, insert it into the sorted array (size ix-2) on its left
!
       j_sorted = ix-1
!
       do while(x_current < this%x(this%sorted2ix(j_sorted-1)))
          this%sorted2ix(j_sorted) = this%sorted2ix(j_sorted-1) ! move the left point by 1 to the right
          j_sorted = j_sorted - 1 ! go one position to the left and check again
       end do
       ! when the do loop exits, j_sorted is the sorting index for the current coordinate
       this%sorted2ix(j_sorted) = ix
    end do ! ix
!
    ! after completing the sorting for sorted2ix, define the inverse mapping ix2sorted
    allocate(this%ix2sorted(this%nx))
    do ix = 1,this%nx
       this%ix2sorted(this%sorted2ix(ix)) = ix
    end do ! ix
!
    ! NOW DO THE VERY SAME THING FOR y AS WAS DONE FOR x ABOVE
!
    allocate(this%sorted2iy(this%ny))
    this%sorted2iy(1) = 1 ! the first value in this%y is most "left" y coordinate (ymin chunk boundary)
    this%sorted2iy(this%ny) = 2 ! the second value in this%y is most "right" y coordinate (ymax chunk boundary)
    do iy = 3,this%ny
       y_current = this%y(iy)
!
       j_sorted = iy-1
!
       do while(y_current < this%y(this%sorted2iy(j_sorted-1)))
          this%sorted2iy(j_sorted) = this%sorted2iy(j_sorted-1) ! move the left point by 1 to the right
          j_sorted = j_sorted - 1 ! go one position to the left and check again
       end do
       this%sorted2iy(j_sorted) = iy
    end do ! iy
!
    ! after completing the sorting for sorted2iy, define the inverse mapping iy2sorted
    allocate(this%iy2sorted(this%ny))
    do iy = 1,this%ny
       this%iy2sorted(this%sorted2iy(iy)) = iy
    end do ! iy
!
    ! now define the sorting for the z coordinate, this%z should be sorted in DECREASING order, since
    ! it was build (at least for the base cells) from surface to bottom of inversion grid (representing
    ! radius, i.e. decreasing values)
!
    allocate(this%sorted2iz(this%nz))
    ! the first nz_base values in this%z are sorted in decreasing order (by construction of the base cells)
    ! place the last nz_base-1 values in the beginning of the sorted array (in reversed order, i.e. increasing coordinates)
    ! and place the first value (maximum radius value) at the very end.
    ! afterwards (if there are any) sort the remaining refined values from right into the 
    ! leading sorted block (in the same manner as x and y above, however here sorting in DECREASING order)
    this%sorted2iz(1:nz_base-1) = (/ (j_sorted,j_sorted=nz_base,2,-1) /)
    this%sorted2iz(this%nz) = 1
    do iz = nz_base+1,this%nz ! this loop does nothing, if nz_base >= this%nz
       z_current = this%z(iz)
!
       j_sorted = iz-1
!
       do while(z_current < this%z(this%sorted2iz(j_sorted-1)))
          this%sorted2iz(j_sorted) = this%sorted2iz(j_sorted-1)
          j_sorted = j_sorted -1
       end do
       this%sorted2iz(j_sorted) = iz
    end do ! iz
!
    ! after completing the sorting for sorted2iz, define the inverse mapping iz2sorted
    allocate(this%iz2sorted(this%nz))
    do iz = 1,this%nz
       this%iz2sorted(this%sorted2iz(iz)) = iz
    end do ! iz
  end subroutine sortCellBoundaryCoordinatesChunksInversionGrid
!--------------------------------------------------------------------
!> \brief find face neighbours of base cells
  subroutine findFaceNeighboursOfBaseCellsChunksInversionGrid(this,nblock,nlay_block,nlon,nlat,ncell_base_1chunk,&
         icell_base_1chunk)
    type (chunks_inversion_grid) :: this
    integer :: nblock,ncell_base_1chunk
    integer, dimension(nblock) :: nlay_block,nlon,nlat
    integer, dimension(:,:,:,:) :: icell_base_1chunk
    ! local
    integer :: ibl,ilay,icell,ichunk,ichunk_nb,icell_shift,icell_shift_nb,reverse_orientation,n_nb,i_nb
    integer :: ix,ix_nb,ix_nb_1,ix_nb_2,nx,nx_above,nx_below,iy,iy_nb,iy_nb_1,iy_nb_2,ny,ny_above,ny_below
    logical :: layer_is_bottom,layer_is_top
    character(len=4) :: bound_nb
    integer, dimension(:), pointer :: idx_nb

    nullify(idx_nb)

    ! NOW DEFINE NEIGHBOUR ARRAY FOR ALL BASE CELLS, LOOP ON 1 CHUNK AND DEFINE FOR ALL CHUNS (simple index shift)

    allocate(this%cell_nb_xmin(this%ncell_internal),this%cell_nb_xmax(this%ncell_internal),&
         this%cell_nb_ymin(this%ncell_internal),this%cell_nb_ymax(this%ncell_internal),&
         this%cell_nb_zmin(this%ncell_internal),this%cell_nb_zmax(this%ncell_internal))

    icell = 0

    do ibl = 1,nblock
       do ilay = 1,nlay_block(ibl)
          layer_is_top = ilay==1
          layer_is_bottom = ilay==nlay_block(ibl)
          ny = nlon(ibl)
          if(ibl/=nblock) ny_below = nlon(ibl+1)
          if(ibl/=1) ny_above = nlon(ibl-1)
          do iy = 1,ny
             nx = nlat(ibl)
             if(ibl/=nblock) nx_below = nlat(ibl+1)
             if(ibl/=1) nx_above = nlat(ibl-1)
             do ix = 1,nx
                ! index of the current cell in first (or only) chunk (indices for n'th chunk shifted by (n-1)*ncell_base_1chunk)
                icell = icell + 1

         ! xmin neighbours of this cell (define for all chunks by simple index shift)

                if(ix==1) then 
                   ! for xmin neighbours on the current xmin boundary, check if there is a neighbouring chunk

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      call getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,'xmin',ichunk_nb,bound_nb,&
                           reverse_orientation) ! reverse_orientation is 0 (no reverse) or 1 (reverse); use this in summation operations below to correctly compute indices on reverted boundaries
                      if(ichunk_nb > 0) then
                         select case(bound_nb)
                         case('xmin')
                            ix_nb = 1
                            iy_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                         case('xmax')
                            ix_nb = nx
                            iy_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                         case('ymin')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                            iy_nb = 1
                         case('ymax')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                            iy_nb = ny
                         end select
                         icell_shift_nb = (ichunk_nb-1)*ncell_base_1chunk
                         allocate(idx_nb(1))
                         idx_nb(1) = icell_base_1chunk(ix_nb,iy_nb,ilay,ibl) + icell_shift_nb
                         call associateVectorPointer(this%cell_nb_xmin(icell+icell_shift),idx_nb)
                         nullify(idx_nb)                         
                      !else ! ichunk_nb > 0
                         ! if there is no touching boundary, there is no xmin neighbour cell;
                         ! in this case, nothing needs to be done, since the pointer
                         ! this%cell_nb_xmin(icell+icell_shift) must stay unassociated (as it is initially)
                      end if ! ichunk_nb > 0
                   end do ! ichunk

                else ! ix==1
                   ! this base cell has exactly 1 xmin neighbour, which is contained in the same chunk (so use same icell_shift for base cell and neighbour cell)

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix-1,iy,ilay,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_xmin(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end if ! ix==1

         ! xmax neighbours of this cell (define for all chunks by simple index shift)

                if(ix == nx) then
                   ! for xmax neighbours on the current xmax boundary, check if there is a neighbouring chunk

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      call getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,'xmax',ichunk_nb,bound_nb,&
                           reverse_orientation) ! reverse_orientation is 0 (no reverse) or 1 (reverse); use this in summation operations below to correctly compute indices on reverted boundaries
                      if(ichunk_nb > 0) then
                         select case(bound_nb)
                         case('xmin')
                            ix_nb = 1
                            iy_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                         case('xmax')
                            ix_nb = nx
                            iy_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                         case('ymin')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                            iy_nb = 1
                         case('ymax')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = (1-reverse_orientation)*iy + reverse_orientation*(ny-iy+1)
                            iy_nb = ny
                         end select
                         icell_shift_nb = (ichunk_nb-1)*ncell_base_1chunk
                         allocate(idx_nb(1))
                         idx_nb(1) = icell_base_1chunk(ix_nb,iy_nb,ilay,ibl) + icell_shift_nb
                         call associateVectorPointer(this%cell_nb_xmax(icell+icell_shift),idx_nb)
                         nullify(idx_nb)                         
                      !else ! ichunk_nb > 0
                         ! if there is no touching boundary, there is no xmax neighbour cell;
                         ! in this case, nothing needs to be done, since the pointer
                         ! this%cell_nb_xmax(icell+icell_shift) must stay unassociated (as it is initially)
                      end if ! ichunk_nb > 0
                   end do ! ichunk                   

                else ! ix == nx
                   ! this base cell has exactly 1 xmax neighbour, which is contained in the same chunk (so use same icell_shift for base cell and neighbour cell)

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix+1,iy,ilay,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_xmax(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end if ! ix == nx

         ! ymin neighbours of this cell (define for all chunks by simple index shift)

                if(iy==1) then 
                   ! for ymin neighbours on the current ymin boundary, check if there is a neighbouring chunk

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      call getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,'ymin',ichunk_nb,bound_nb,&
                           reverse_orientation) ! reverse_orientation is 0 (no reverse) or 1 (reverse); use this in summation operations below to correctly compute indices on reverted boundaries
                      if(ichunk_nb > 0) then
                         select case(bound_nb)
                         case('xmin')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = 1
                            iy_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                         case('xmax')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = nx
                            iy_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                         case('ymin')
                            ix_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                            iy_nb = 1
                         case('ymax')
                            ix_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                            iy_nb = ny
                         end select
                         icell_shift_nb = (ichunk_nb-1)*ncell_base_1chunk
                         allocate(idx_nb(1))
                         idx_nb(1) = icell_base_1chunk(ix_nb,iy_nb,ilay,ibl) + icell_shift_nb
                         call associateVectorPointer(this%cell_nb_ymin(icell+icell_shift),idx_nb)
                         nullify(idx_nb)                         
                      !else ! ichunk_nb > 0
                         ! if there is no touching boundary, there is no ymin neighbour cell;
                         ! in this case, nothing needs to be done, since the pointer
                         ! this%cell_nb_ymin(icell+icell_shift) must stay unassociated (as it is initially)
                      end if ! ichunk_nb > 0
                   end do ! ichunk

                else ! iy==1
                   ! this base cell has exactly 1 ymin neighbour, which is contained in the same chunk (so use same icell_shift for base cell and neighbour cell)

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix,iy-1,ilay,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_ymin(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end if ! iy==1

         ! ymax neighbours of this cell (define for all chunks by simple index shift)

                if(iy == ny) then
                   ! for ymax neighbours on the current ymax boundary, check if there is a neighbouring chunk

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      call getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,'ymax',ichunk_nb,bound_nb,&
                           reverse_orientation) ! reverse_orientation is 0 (no reverse) or 1 (reverse); use this in summation operations below to correctly compute indices on reverted boundaries
                      if(ichunk_nb > 0) then
                         select case(bound_nb)
                         case('xmin')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = 1
                            iy_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                         case('xmax')
                            ! in case of nchunk > 1 (which is the case here, since there exists a neighbouring chunk), 
                            ! it is a priorily assured that always nx==ny, so no problem to match x-indices with y-indices
                            ix_nb = nx
                            iy_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                         case('ymin')
                            ix_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                            iy_nb = 1
                         case('ymax')
                            ix_nb = (1-reverse_orientation)*ix + reverse_orientation*(nx-ix+1)
                            iy_nb = ny
                         end select
                         icell_shift_nb = (ichunk_nb-1)*ncell_base_1chunk
                         allocate(idx_nb(1))
                         idx_nb(1) = icell_base_1chunk(ix_nb,iy_nb,ilay,ibl) + icell_shift_nb
                         call associateVectorPointer(this%cell_nb_ymax(icell+icell_shift),idx_nb)
                         nullify(idx_nb)                         
                      !else ! ichunk_nb > 0
                         ! if there is no touching boundary, there is no ymax neighbour cell;
                         ! in this case, nothing needs to be done, since the pointer
                         ! this%cell_nb_ymax(icell+icell_shift) must stay unassociated (as it is initially)
                      end if ! ichunk_nb > 0
                   end do ! ichunk                   

                else ! iy == ny
                   ! this base cell has exactly 1 ymax neighbour, which is contained in the same chunk (so use same icell_shift for base cell and neighbour cell)

                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix,iy+1,ilay,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_ymax(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end if ! iy == ny

         ! zmin neighbours of this cell (define for all chunks by simple index shift)
         ! by construction of the chunks, all z-neighbours of base cells are contained in the same chunk 
         ! (so use same icell_shift for base cell and neighbour cell)

         ! z-coordinate is RADIUS, but the loop is going from Earth's surface downwards;
         ! hence, the zmin neighbour is BELOW the current cell
                if(layer_is_bottom) then

                   ! search in top layer of block below (if any) for zmin neighbours
                   ! find there a covering range of cells in x- and y- direction, dependent on nx,ny,ix,iy, ...
                   if(ibl/=nblock) then
                      ix_nb_1 = ((ix-1)*nx_below)/nx + 1
                      ix_nb_2 = nx_below - ((nx-ix)*nx_below)/nx
                      iy_nb_1 = ((iy-1)*ny_below)/ny + 1
                      iy_nb_2 = ny_below - ((ny-iy)*ny_below)/ny
                      n_nb = (ix_nb_2-ix_nb_1+1)*(iy_nb_2-iy_nb_1+1)
                      if(ix_nb_2 < ix_nb_1 .or. iy_nb_2 < iy_nb_1) then
                         ! THIS IS ACTUALLY BAD PRACTICE, BETTER RAISE A REAL ERROR HERE (well, hope this will not be necessary...)
                         write(*,*) "ERROR 1 IN findFaceNeighboursOfBaseCellsChunksInversionGrid: ",&
                              "THIS ERROR SHOULD NOT OCCUR: number of zmin neigbhour base cells = ",n_nb
                         stop
                      end if
                      do ichunk = 1,this%nchunk
                         icell_shift = (ichunk-1)*ncell_base_1chunk
                         allocate(idx_nb(n_nb))
                         i_nb = 0
                         do iy_nb = iy_nb_1,iy_nb_2
                            do ix_nb = ix_nb_1,ix_nb_2
                               i_nb = i_nb + 1
                               idx_nb(i_nb) = icell_base_1chunk(ix_nb,iy_nb,1,ibl+1) + icell_shift
                            end do ! iy_nb
                         end do ! ix_nb
                         call associateVectorPointer(this%cell_nb_zmin(icell+icell_shift),idx_nb)
                         nullify(idx_nb)
                      end do ! ichunk

                   ! else, if this is the bottom layer of the bottom block, there are no zmin neighbours, 
                   ! so do nothing (as vector pointer is nullified already)
                   end if ! ibl/=nblock

                else ! layer_is_bottom

                   ! this base cell has exactly 1 zmin neighbour
                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix,iy,ilay+1,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_zmin(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end  if ! layer_is_bottom

         ! zmax neighbours of this cell (define for all chunks by simple index shift)
         ! by construction of the chunks, all z-neighbours of base cells are contained in the same chunk 
         ! (so use same icell_shift for base cell and neighbour cell)

         ! z-coordinate is RADIUS, but the loop is going from Earth's surface downwards;
         ! hence, the zmax neighbour is ABOVE the current cell

                if(layer_is_top) then

                   ! search in bottom layer of block above (if any) for zmax neighbours
                   ! find there a covering range of cells in x- and y- direction, dependent on nx,ny,ix,iy, ...
                   if(ibl/=1) then
                      ix_nb_1 = ((ix-1)*nx_above)/nx + 1
                      ix_nb_2 = nx_above - ((nx-ix)*nx_above)/nx
                      iy_nb_1 = ((iy-1)*ny_above)/ny + 1
                      iy_nb_2 = ny_above - ((ny-iy)*ny_above)/ny
                      n_nb = (ix_nb_2-ix_nb_1+1)*(iy_nb_2-iy_nb_1+1)
                      if(ix_nb_2 < ix_nb_1 .or. iy_nb_2 < iy_nb_1) then
                         ! THIS IS ACTUALLY BAD PRACTICE, BETTER RAISE A REAL ERROR HERE (well, hope this will not be necessary...)
                         write(*,*) "ERROR 2 IN findFaceNeighboursOfBaseCellsChunksInversionGrid: ",&
                              "THIS ERROR SHOULD NOT OCCUR: number of zmax neigbhour base cells = ",n_nb
                         stop
                      end if
                      do ichunk = 1,this%nchunk
                         icell_shift = (ichunk-1)*ncell_base_1chunk
                         allocate(idx_nb(n_nb))
                         i_nb = 0
                         do iy_nb = iy_nb_1,iy_nb_2
                            do ix_nb = ix_nb_1,ix_nb_2
                               i_nb = i_nb + 1
                               idx_nb(i_nb) = icell_base_1chunk(ix_nb,iy_nb,nlay_block(ibl-1),ibl-1) + icell_shift
                            end do ! iy_nb
                         end do ! ix_nb
                         call associateVectorPointer(this%cell_nb_zmax(icell+icell_shift),idx_nb)
                         nullify(idx_nb)
                      end do ! ichunk

                   ! else, if this is the top layer of the top block, there are no zmax neighbours, 
                   ! so do nothing (as vector pointer is nullified already)
                   end if ! ibl/=1

                else ! layer_is_top

                   ! this base cell has exactly 1 zmax neighbour
                   do ichunk = 1,this%nchunk
                      icell_shift = (ichunk-1)*ncell_base_1chunk
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_base_1chunk(ix,iy,ilay-1,ibl) + icell_shift
                      call associateVectorPointer(this%cell_nb_zmax(icell+icell_shift),idx_nb)
                      nullify(idx_nb)
                   end do ! ichunk

                end  if ! layer_is_top

             end do ! ix
          end do ! iy
       end do ! ilay
    end do ! ibl
!
  end subroutine findFaceNeighboursOfBaseCellsChunksInversionGrid
!--------------------------------------------------------------------
!> \brief get index and boundary of LATERALLY touching neighbouring chunk (if any)
  subroutine getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,bound,ichunk_nb,bound_nb,reverse_orientation)
    type (chunks_inversion_grid) :: this
    integer :: ichunk,ichunk_nb,reverse_orientation
    character(len=*), intent(in) :: bound
    character(len=4), intent(out) :: bound_nb
!
    ichunk_nb = -1
    bound_nb = ''
    reverse_orientation = 0
!
    select case(this%nchunk)

    case(1) ! this%nchunk

       ! there are no neighbouring chunks (regardless of input ichunk,bound)
       return

    case(2) ! this%nchunk

       select case(ichunk)

       case(1) ! ichunk
          ! in this case, there is only one boundary with a neighbour, namely ymin
          select case(bound)
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'ymax'
          end select ! bound

       case(2) ! ichunk
          ! in this case, there is only one boundary with a neighbour, namely ymax
          select case(bound)
          case('ymax')
             ichunk_nb = 1
             bound_nb = 'ymin'
          end select ! bound

       end select ! ichunk

    case(3) ! this%nchunk
       
       select case(ichunk)

       case(1) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'xmax'
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'ymax'
          end select ! bound

       case(2) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'ymin'
          case('ymax')
             ichunk_nb = 1
             bound_nb = 'ymin'
          end select ! bound

       case(3) ! ichunk
          select case(bound)
          case('xmax')
             ichunk_nb = 1
             bound_nb = 'xmin'
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'xmin'
          end select ! bound

       end select ! ichunk

    case(6) ! this%nchunk

       select case(ichunk)

       case(1) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'xmax'
          case('xmax')
             ichunk_nb = 6
             bound_nb = 'xmin'
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'ymax'
          case('ymax')
             ichunk_nb = 5
             bound_nb = 'ymin'
          end select ! bound

       case(2) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'ymin'
          case('xmax')
             ichunk_nb = 6
             bound_nb = 'ymin'
             reverse_orientation = 1
          case('ymin')
             ichunk_nb = 4
             bound_nb = 'ymax'
          case('ymax')
             ichunk_nb = 1
             bound_nb = 'ymin'
          end select ! bound

       case(3) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 4
             bound_nb = 'xmin'
             reverse_orientation = 1
          case('xmax')
             ichunk_nb = 1
             bound_nb = 'xmin'
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'xmin'
          case('ymax')
             ichunk_nb = 5
             bound_nb = 'xmin'
             reverse_orientation = 1
          end select ! bound

       case(4) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'xmin'
             reverse_orientation = 1
          case('xmax')
             ichunk_nb = 6
             bound_nb = 'xmax'
             reverse_orientation = 1
          case('ymin')
             ichunk_nb = 5
             bound_nb = 'ymax'
          case('ymax')
             ichunk_nb = 2
             bound_nb = 'ymin'
          end select ! bound

       case(5) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 3
             bound_nb = 'ymax'
             reverse_orientation = 1
          case('xmax')
             ichunk_nb = 6
             bound_nb = 'ymax'
          case('ymin')
             ichunk_nb = 1
             bound_nb = 'ymax'
          case('ymax')
             ichunk_nb = 4
             bound_nb = 'ymin'
          end select ! bound

       case(6) ! ichunk
          select case(bound)
          case('xmin')
             ichunk_nb = 1
             bound_nb = 'xmax'
          case('xmax')
             ichunk_nb = 4
             bound_nb = 'xmax'
             reverse_orientation = 1
          case('ymin')
             ichunk_nb = 2
             bound_nb = 'xmax'
             reverse_orientation = 1
          case('ymax')
             ichunk_nb = 5
             bound_nb = 'xmax'
          end select ! bound

       end select ! ichunk

    end select ! this%nchunk
  end subroutine getTouchingBoundaryOfNbChunkChunksInversionGrid
!--------------------------------------------------------------------
!> \brief find face neighbours of refined cells
  subroutine findFaceNeighboursOfRefinedCellsChunksInversionGrid(this)
    type (chunks_inversion_grid) :: this
    ! local
    integer :: icell_base,ibound,isub,nsub,inb,nnb
    integer, dimension(:), pointer :: idx_subcells,idx_nb,ixmin_sub,ixmax_sub,iymin_sub,iymax_sub,izmin_sub,izmax_sub
    character(len=4) :: bound,bound_nb
    integer, dimension(:), pointer :: nb_base,idx_sub_base_at_bound,idx_sub_nb_base_at_bound
    integer :: icell_nb_base,ichunk,ichunk_nb,reverse_orientation,c2_min,c2_max,c3_min,c3_max
    logical :: base_cell_is_refined,nb_base_cell_is_refined
    logical, dimension(:), allocatable :: condition_c1,condition_c1_c2_c3
!
    nullify(idx_subcells,idx_nb,ixmin_sub,ixmax_sub,iymin_sub,iymax_sub,izmin_sub,izmax_sub)
    nullify(nb_base,idx_sub_base_at_bound,idx_sub_nb_base_at_bound)
!
    do icell_base = 1,this%ncell_base
       idx_subcells => getVectorPointer(this%cell_subcells(icell_base))
       if(.not.associated(idx_subcells)) cycle
!
       nsub = size(idx_subcells)
       if(associated(ixmin_sub)) deallocate(ixmin_sub)
       if(associated(ixmax_sub)) deallocate(ixmax_sub)
       if(associated(iymin_sub)) deallocate(iymin_sub)
       if(associated(iymax_sub)) deallocate(iymax_sub)
       if(associated(izmin_sub)) deallocate(izmin_sub)
       if(associated(izmax_sub)) deallocate(izmax_sub)
       allocate(ixmin_sub(nsub),ixmax_sub(nsub),iymin_sub(nsub),iymax_sub(nsub),izmin_sub(nsub),izmax_sub(nsub))
       ixmin_sub = this%cell_ixmin(idx_subcells)
       ixmax_sub = this%cell_ixmax(idx_subcells)
       iymin_sub = this%cell_iymin(idx_subcells)
       iymax_sub = this%cell_iymax(idx_subcells)
       izmin_sub = this%cell_izmin(idx_subcells)
       izmax_sub = this%cell_izmax(idx_subcells)
!
       if(allocated(condition_c1)) deallocate(condition_c1)
       if(allocated(condition_c1_c2_c3)) deallocate(condition_c1_c2_c3)
       allocate(condition_c1(nsub),condition_c1_c2_c3(nsub))
       do isub = 1,nsub
          !find xmin neighbours of cell isub
          condition_c1(:) = ixmax_sub(:) == ixmin_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = y, c3 = z
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%iy2sorted(iymin_sub(isub)),this%iy2sorted(iymax_sub(isub)),& ! c2 test cell
                  this%iz2sorted(izmin_sub(isub)),this%iz2sorted(izmax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%iy2sorted(iymin_sub),this%iy2sorted(iymax_sub),& ! c2 array
                  this%iz2sorted(izmin_sub),this%iz2sorted(izmax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_xmin(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for xmin neighbours
!
          !find xmax neighbours of cell isub
          condition_c1(:) = ixmin_sub(:) == ixmax_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = y, c3 = z
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%iy2sorted(iymin_sub(isub)),this%iy2sorted(iymax_sub(isub)),& ! c2 test cell
                  this%iz2sorted(izmin_sub(isub)),this%iz2sorted(izmax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%iy2sorted(iymin_sub),this%iy2sorted(iymax_sub),& ! c2 array
                  this%iz2sorted(izmin_sub),this%iz2sorted(izmax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_xmax(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for xmax neighbours
!
          !find ymin neighbours of cell isub
          condition_c1(:) = iymax_sub(:) == iymin_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = x, c3 = z
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%ix2sorted(ixmin_sub(isub)),this%ix2sorted(ixmax_sub(isub)),& ! c2 test cell
                  this%iz2sorted(izmin_sub(isub)),this%iz2sorted(izmax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%ix2sorted(ixmin_sub),this%ix2sorted(ixmax_sub),& ! c2 array
                  this%iz2sorted(izmin_sub),this%iz2sorted(izmax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_ymin(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for ymin neighbours
!
          !find ymax neighbours of cell isub
          condition_c1(:) = iymin_sub(:) == iymax_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = x, c3 = z
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%ix2sorted(ixmin_sub(isub)),this%ix2sorted(ixmax_sub(isub)),& ! c2 test cell
                  this%iz2sorted(izmin_sub(isub)),this%iz2sorted(izmax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%ix2sorted(ixmin_sub),this%ix2sorted(ixmax_sub),& ! c2 array
                  this%iz2sorted(izmin_sub),this%iz2sorted(izmax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_ymax(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for ymax neighbours
!
          !find zmin neighbours of cell isub
          condition_c1(:) = izmax_sub(:) == izmin_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = x, c3 = y
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%ix2sorted(ixmin_sub(isub)),this%ix2sorted(ixmax_sub(isub)),& ! c2 test cell
                  this%iy2sorted(iymin_sub(isub)),this%iy2sorted(iymax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%ix2sorted(ixmin_sub),this%ix2sorted(ixmax_sub),& ! c2 array
                  this%iy2sorted(iymin_sub),this%iy2sorted(iymax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_zmin(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for zmin neighbours
!
          !find zmax neighbours of cell isub
          condition_c1(:) = izmin_sub(:) == izmax_sub(isub) ! select those indices which satisfy condition on first coordinate
          if(any(condition_c1)) then
             ! HERE: c2 = x, c3 = y
             call findOverlapping2DFacesChunksInversionGrid(&
                  this%ix2sorted(ixmin_sub(isub)),this%ix2sorted(ixmax_sub(isub)),& ! c2 test cell
                  this%iy2sorted(iymin_sub(isub)),this%iy2sorted(iymax_sub(isub)),& ! c3 test cell
                  nsub,&
                  this%ix2sorted(ixmin_sub),this%ix2sorted(ixmax_sub),& ! c2 array
                  this%iy2sorted(iymin_sub),this%iy2sorted(iymax_sub),& ! c3 array
                  condition_c1_c2_c3,condition_c1=condition_c1)
             nnb = count(condition_c1_c2_c3)
             if(nnb>0) then
                allocate(idx_nb(nnb))
                idx_nb = pack( idx_subcells , condition_c1_c2_c3 )
                call associateVectorPointer(this%cell_nb_zmax(idx_subcells(isub)),idx_nb)
                nullify(idx_nb)
             end if
          end if ! any(condition_c1) for zmax neighbours
       end do ! isub
    end do ! icell_base
!
    if(associated(ixmin_sub)) deallocate(ixmin_sub)
    if(associated(ixmax_sub)) deallocate(ixmax_sub)
    if(associated(iymin_sub)) deallocate(iymin_sub)
    if(associated(iymax_sub)) deallocate(iymax_sub)
    if(associated(izmin_sub)) deallocate(izmin_sub)
    if(associated(izmax_sub)) deallocate(izmax_sub)
    if(allocated(condition_c1)) deallocate(condition_c1)
!
    ! re-loop over all base cells, updating neighbour arrays of all subcells (and possibly the base cell 
    ! itself) w.r.t. subcells of neighbouring base cells; pay special attention if the neighbouring base 
    ! cell is in a different chunk (possibly inverse orientation of lateral coordinates)
    do icell_base = 1,this%ncell_base
       base_cell_is_refined = associated(getVectorPointer(this%cell_subcells(icell_base)))
       ichunk = this%cell_ichunk(icell_base)
!
       ! Loop on all 4 lateral boundaries of this base cell and update the neighbours of this base cell 
       ! (if not refined) or the neighbours of subcells of this base cell (if refined)
       ! For lateral neighbours, there is special treatment required if neighbouring base cells are in 
       ! a different chunk
       do ibound = 1,4
          select case (ibound)
          case(1)
             bound = 'xmin'
             bound_nb = 'xmax'
             nb_base => getVectorPointer(this%cell_nb_xmin(icell_base))
          case(2)
             bound = 'ymin'
             bound_nb = 'ymax'
             nb_base => getVectorPointer(this%cell_nb_ymin(icell_base))
          case(3)
             bound = 'xmax'
             bound_nb = 'xmin'
             nb_base => getVectorPointer(this%cell_nb_xmax(icell_base))
          case(4)
             bound = 'ymax'
             bound_nb = 'ymin'
             nb_base => getVectorPointer(this%cell_nb_ymax(icell_base))
          end select
!
          if(associated(nb_base)) then
             ! Yes, there is a neighbouring base cell
!
             if(size(nb_base)/=1) then
                write(*,*) "ERROR 3 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: base cell ",icell_base,&
                     " has ",size(nb_base)," neighbouring base cells at boundary '",bound,"'. This violates the ",&
                     "internal assumption of module chunksInversionGrid that in case of NCHUNK > 2, NLAT must ",&
                     "equal NLON in all blocks. Hence, module chunksInversionGrid is inconsistent!"
                stop
             end if
             icell_nb_base = nb_base(1)
!
             ! test all subcells of neighbouring base cell if they have a touching face with (the subcells of) this base cell
!
             ! if neighbour is contained in a different chunk, we need to select boundary (and orientation) of neighbour 
             ! base cell (otherwise the boundary of the neighbour cell was defined above by bound_nb already)
             ichunk_nb = this%cell_ichunk(icell_nb_base)
             if(ichunk == ichunk_nb) then
                ! both base cells are in the same chunk: need to consider the opposite boundary of the neighbour,
                ! the value of bound_nb was set above
                if(associated(idx_sub_nb_base_at_bound)) deallocate(idx_sub_nb_base_at_bound)
                call getSubcellsAtBoundaryChunksInversionGrid(icell_nb_base,bound_nb,idx_sub_nb_base_at_bound)
                reverse_orientation = 0
             else
                ! overwrite previously set value of bound_nb by correct boundary of neighbouring chunk
                ! and get reversre_orientation value (0 or 1)
                call getTouchingBoundaryOfNbChunkChunksInversionGrid(this,ichunk,bound,ichunk_nb,bound_nb,reverse_orientation)
                if(associated(idx_sub_nb_base_at_bound)) deallocate(idx_sub_nb_base_at_bound)
                call getSubcellsAtBoundaryChunksInversionGrid(icell_nb_base,bound_nb,idx_sub_nb_base_at_bound)
             end if ! ichunk == ichunk_nb
             nb_base_cell_is_refined = associated(idx_sub_nb_base_at_bound)
!
             if(base_cell_is_refined) then
                ! get subcells of icell_base which are located on boundary bound
                if(associated(idx_sub_base_at_bound)) deallocate(idx_sub_base_at_bound)
                call getSubcellsAtBoundaryChunksInversionGrid(icell_base,bound,idx_sub_base_at_bound)
                if(.not.associated(idx_sub_base_at_bound)) then
                   write(*,*) "ERROR 4 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: base cell ",icell_base,&
                        " has subcells, but none of them is located at boundary '",bound,&
                        "' -> module chunksInversionGrid is inconsistent!"
                   stop
                end if
                ! loop on all subcells located on boundary bound and update their 'bound'-neighbours array
                do isub = 1,size(idx_sub_base_at_bound)
                   if(nb_base_cell_is_refined) then
                      ! call subroutine findOverlapping2DFacesChunksInversionGrid accounting for reverse_orientation.
                      ! The boundaries c1_min,... (on the current base cell side) can be defined according to the value of 'bound'.
                      ! The boundaries of the subcells on the neighbouring side (for which overlap should be checked, 
!
                      ! always choose c3 to be z, as we iterate over x,y boundaries only
                      c3_min = this%iz2sorted(this%cell_izmin(idx_sub_base_at_bound(isub)))
                      c3_max = this%iz2sorted(this%cell_izmax(idx_sub_base_at_bound(isub)))
                      select case(bound)
                      case('xmin','xmax')
                         ! on the 'x'-sides of the base cell, we need to test 'y' coordinates for overlap, so c2 is y
                         c2_min = this%iy2sorted(this%cell_iymin(idx_sub_base_at_bound(isub)))
                         c2_max = this%iy2sorted(this%cell_iymax(idx_sub_base_at_bound(isub)))
                      case('ymin','ymax')
                         ! on the 'y'-sides of the base cell, we need to test 'x' coordinates for overlap, so c2 is x
                         c2_min = this%ix2sorted(this%cell_ixmin(idx_sub_base_at_bound(isub)))
                         c2_max = this%ix2sorted(this%cell_ixmax(idx_sub_base_at_bound(isub)))
                      end select ! bound
!
                      if(allocated(condition_c1_c2_c3)) deallocate(condition_c1_c2_c3)
                      allocate(condition_c1_c2_c3(size(idx_sub_nb_base_at_bound)))
!
                      select case(bound_nb)
                      case('xmin','xmax')
                         ! on the 'x'-sides of the neighbouring base cell, we need to test 'y' coordinates for overlap, so c2 is y
                         ! c3 is z (everywhere, so also here)
                         call findOverlapping2DFacesChunksInversionGrid(&
                              c2_min,c2_max,c3_min,c3_max,size(idx_sub_nb_base_at_bound),&
                              this%iy2sorted(this%cell_iymin(idx_sub_nb_base_at_bound)),& ! c2 min array
                              this%iy2sorted(this%cell_iymax(idx_sub_nb_base_at_bound)),& ! c2 max array
                              this%iz2sorted(this%cell_izmin(idx_sub_nb_base_at_bound)),& ! c3 min array -> always z
                              this%iz2sorted(this%cell_izmax(idx_sub_nb_base_at_bound)),& ! c3 max array -> always z
                              condition_c1_c2_c3,c2_reverse_order=reverse_orientation)
                      case('ymin','ymax')
                         ! on the 'y'-sides of the neighbouring base cell, we need to test 'x' coordinates for overlap, so c2 is x
                         ! c3 is z (everywhere, so also here)
                         call findOverlapping2DFacesChunksInversionGrid(&
                              c2_min,c2_max,c3_min,c3_max,size(idx_sub_nb_base_at_bound),&
                              this%ix2sorted(this%cell_ixmin(idx_sub_nb_base_at_bound)),& ! c2 min array
                              this%ix2sorted(this%cell_ixmax(idx_sub_nb_base_at_bound)),& ! c2 max array
                              this%iz2sorted(this%cell_izmin(idx_sub_nb_base_at_bound)),& ! c3 min array -> always z
                              this%iz2sorted(this%cell_izmax(idx_sub_nb_base_at_bound)),& ! c3 max array -> always z
                              condition_c1_c2_c3,c2_reverse_order=reverse_orientation)
                      end select ! bound_nb
!
                      nnb = count(condition_c1_c2_c3)
                      if(nnb>0) then
                         allocate(idx_nb(nnb))
                         idx_nb = pack( idx_sub_nb_base_at_bound , condition_c1_c2_c3 )
                      else
                         write(*,*) "ERROR 5 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external_index, contained in base cell ",&
                              icell_base,") does not find any '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell (which has one refined base-cell neighbour) ",&
                              "-> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                   else ! nb_base_cell_is_refined
                      ! the only neighbour to be added is the neighbouring base cell
                      allocate(idx_nb(1))
                      idx_nb(1) = icell_nb_base
                   end if ! nb_base_cell_is_refined
!
                   ! update the 'bound'-neighbours array by idx_nb
                   select case(ibound)
                   case(1)
                      if(associated(getVectorPointer(this%cell_nb_xmin(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 6 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined",&
                              " -> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_xmin(idx_sub_base_at_bound(isub)),idx_nb)
                   case(2)
                      if(associated(getVectorPointer(this%cell_nb_ymin(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 7 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined",&
                              " -> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_ymin(idx_sub_base_at_bound(isub)),idx_nb)
                   case(3)
                      if(associated(getVectorPointer(this%cell_nb_xmax(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 8 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined",&
                              " -> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_xmax(idx_sub_base_at_bound(isub)),idx_nb)
                   case(4)
                      if(associated(getVectorPointer(this%cell_nb_ymax(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 9 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined",&
                              " -> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_ymax(idx_sub_base_at_bound(isub)),idx_nb)
                   end select ! ibound
                   nullify(idx_nb)
                end do ! isub
                ! remove neighbour of current base cell, since the base cell is not external anymore (was refined)
                ! -> only keep external neighbours for simplicity
                select case(ibound)
                case(1); call dealloc(this%cell_nb_xmin(icell_base))
                case(2); call dealloc(this%cell_nb_ymin(icell_base))
                case(3); call dealloc(this%cell_nb_xmax(icell_base))
                case(4); call dealloc(this%cell_nb_ymax(icell_base))
                end select ! ibound
                nullify(nb_base)
!
             else ! base_cell_is_refined
!
                if(nb_base_cell_is_refined) then
                   ! modify 'bound'-neighbours of this base cell:
                   ! remove neighbouring base cell index (since neighbouring base cell is not external, since it was refined)
                   ! and add all subcell indices at the touching boundary
                   if(size(idx_sub_nb_base_at_bound)==1) then
                      ! directly overwrite value in respective array this%cell_nb_'bound'(icell_base)
                      nb_base(1) = idx_sub_nb_base_at_bound(1)
                   else
                      nb_base => reallocate(nb_base,size(idx_sub_nb_base_at_bound))
                      nb_base = idx_sub_nb_base_at_bound
                      ! reallocation has changed the memory location of the neighbour array and deallocated the
                      ! old memory, so need to re-associated this%cell_nb_'bound'(icell_base)
                      select case(ibound)
                      case(1); call associateVectorPointer(this%cell_nb_xmin(icell_base),nb_base)
                      case(2); call associateVectorPointer(this%cell_nb_ymin(icell_base),nb_base)
                      case(3); call associateVectorPointer(this%cell_nb_xmax(icell_base),nb_base)
                      case(4); call associateVectorPointer(this%cell_nb_ymax(icell_base),nb_base)
                      end select ! ibound
                   end if
                else ! nb_base_cell_is_refined
                   ! the only neighbour to be added would be the neighbouring base cell, which already is
                   ! set as neighbour, so nothing to do here
                end if ! nb_base_cell_is_refined
             end if ! base_cell_is_refined

          else ! associated(nb_base)
             ! since there is no 'bound'-neighbour of this base cell, there is nothing to do at this point
          end if ! associated(nb_base)
       end do ! ibound, xmin,ymin,xmax,ymax
!
       ! account for the vertical boundaries zmin,zmax
       do ibound = 5,6
          select case (ibound)
          case(5)
             bound = 'zmin'
             bound_nb = 'zmax'
             nb_base => getVectorPointer(this%cell_nb_zmin(icell_base))
          case(6)
             bound = 'zmax'
             bound_nb = 'zmin'
             nb_base => getVectorPointer(this%cell_nb_zmax(icell_base))
          end select
!
          if(associated(nb_base)) then
             ! Yes, there is at least one neighbouring base cell

             ! test all subcells of neighbouring base cell(s) if they have a touching face with (the subcells of) this base cell
!
             ! in case of vertical neighbours, there can be more than one neighbouring base cell, so collect
             ! potential subcells
             if(associated(idx_sub_nb_base_at_bound)) deallocate(idx_sub_nb_base_at_bound)
             nsub = 0
             nb_base_cell_is_refined = .false.
             do inb = 1,size(nb_base)
                if(associated(idx_nb)) deallocate(idx_nb)
                call getSubcellsAtBoundaryChunksInversionGrid(nb_base(inb),bound_nb,idx_nb)
                if(associated(idx_nb)) then
                   ! base cell nb_base(inb) was refined, so append idx to idx_sub_nb_base_at_bound
                   nnb = size(idx_nb)
                   idx_sub_nb_base_at_bound => reallocate(idx_sub_nb_base_at_bound,nsub+nnb)
                   idx_sub_nb_base_at_bound(nsub+1:nsub+nnb) = idx_nb
                   nsub = nsub+nnb
                   ! remember that there is at least one neighbouring base cell which contains subcells
                   ! (otherwise we do not need to test for overlapping faces below, but simply can set
                   ! the base cells as neighbours)
                   nb_base_cell_is_refined = .true. 
                else
                   ! base cell nb_base(inb) is NOT refined, so append base cell index nb_base(inb) to idx_sub_nb_base_at_bound
                   idx_sub_nb_base_at_bound => reallocate(idx_sub_nb_base_at_bound,nsub+1)
                   idx_sub_nb_base_at_bound(nsub+1) = nb_base(inb)
                   nsub = nsub+1
                end if
             end do ! inb
             if(associated(idx_nb)) deallocate(idx_nb)
!
             if(base_cell_is_refined) then
                ! get subcells of icell_base which are located on boundary bound
                call getSubcellsAtBoundaryChunksInversionGrid(icell_base,bound,idx_sub_base_at_bound)
                if(.not.associated(idx_sub_base_at_bound)) then
                   write(*,*) "ERROR 10 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: base cell ",icell_base,&
                        " has subcells, but none of them is located at boundary '",bound,&
                        "' -> module chunksInversionGrid is inconsistent!"
                   stop
                end if
                ! loop on all subcells located on boundary bound and update their 'bound'-neighbours array
                do isub = 1,size(idx_sub_base_at_bound)
                   if( (size(idx_sub_base_at_bound)==1.and.(.not.nb_base_cell_is_refined .or. size(nb_base)==1)) .or. &
                        size(idx_sub_nb_base_at_bound)==1 ) then
                      ! in these cases, we do not need to search for overlapping 2D faces:
                      ! 1) if there is just one subcell and only base cells as neighbours (or one base cell which is allowed 
                      !    to be refined), assign vector idx_sub_nb_base_at_bound as neighbours
                      ! 2) if there is just one neighbour on the face (either base cell or subcell), it 
                      !    automatically overlaps with all subcells here
                      !
                      ! in these cases, do allocate new memory for idx_nb and DO NOT simply 
                      ! idx_nb => idx_sub_nb_base_at_bound, since this must be done for ALL subcells isub 
                      ! (otherwise, this operation only works for the first subcell isub==1).
                      allocate(idx_nb(size(idx_sub_nb_base_at_bound)))
                      idx_nb = idx_sub_nb_base_at_bound
                   else ! (size(idx_sub_base_at_bound)==1.and.(.not.nb_base_cell_is_refined ....
                      ! we need to search for overlapping 2D faces, because even if all neighbouring base cells 
                      ! (there is more than one!) are not refined, it is not clear whether all of them overlap 
                      ! with all subcells of this base cell
!
                      if(allocated(condition_c1_c2_c3)) deallocate(condition_c1_c2_c3)
                      allocate(condition_c1_c2_c3(size(idx_sub_nb_base_at_bound)))
!
                      ! call subroutine findOverlapping2DFacesChunksInversionGrid, with c3 = x, c3 = y 
                      ! (c1 is always z, as we are on a z-boundary)
                      ! There is no reverse orientation here, since in z direction we are always in the same chunks
                      ! with same orientations of x- and y- coordinates
                      call findOverlapping2DFacesChunksInversionGrid(&
                           this%ix2sorted(this%cell_ixmin(idx_sub_base_at_bound(isub))),& ! c2 min
                           this%ix2sorted(this%cell_ixmax(idx_sub_base_at_bound(isub))),& ! c2 max
                           this%iy2sorted(this%cell_iymin(idx_sub_base_at_bound(isub))),& ! c3 min
                           this%iy2sorted(this%cell_iymax(idx_sub_base_at_bound(isub))),& ! c3 max
                           size(idx_sub_nb_base_at_bound),&
                           this%ix2sorted(this%cell_ixmin(idx_sub_nb_base_at_bound)),& ! c2 min array
                           this%ix2sorted(this%cell_ixmax(idx_sub_nb_base_at_bound)),& ! c2 max array
                           this%iy2sorted(this%cell_iymin(idx_sub_nb_base_at_bound)),& ! c3 min array
                           this%iy2sorted(this%cell_iymax(idx_sub_nb_base_at_bound)),& ! c3 max array
                           condition_c1_c2_c3)
                      nnb = count(condition_c1_c2_c3)
                      if(nnb>0) then
                         allocate(idx_nb(nnb))
                         idx_nb = pack( idx_sub_nb_base_at_bound , condition_c1_c2_c3 )
                      else
                         write(*,*) "ERROR 11 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") does not find any '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell (which has base-cell neighbour(s)) ",&
                              "-> module chunksInversionGrid is inconsistent!"
                         stop
                      end if
                   end if ! (size(idx_sub_base_at_bound)==1.and.(.not.nb_base_cell_is_refined ....
!
                   ! update the 'bound'-neighbours array by idx_nb
                   select case(ibound)
                   case(5)
                      if(associated(getVectorPointer(this%cell_nb_zmin(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 12 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined -> module ",&
                              "chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_zmin(idx_sub_base_at_bound(isub)),idx_nb)
                   case(6)
                      if(associated(getVectorPointer(this%cell_nb_zmax(idx_sub_base_at_bound(isub))))) then
                         write(*,*) "ERROR 13 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: subcell ",&
                              this%icell_external(idx_sub_base_at_bound(isub))," (external index, contained in base cell ",&
                              icell_base,") already has '",bound,"' neighbours, although it is located on the '",bound,&
                              "' face of the base cell. This should not yet have been defined -> module ",&
                              "chunksInversionGrid is inconsistent!"
                         stop
                      end if
                      call associateVectorPointer(this%cell_nb_zmax(idx_sub_base_at_bound(isub)),idx_nb)
                   end select ! ibound
                   nullify(idx_nb)
                end do ! isub
                ! remove neighbour of current base cell, since the base cell is not external anymore (was refined)
                ! -> only keep external neighbours for simplicity
                select case(ibound)
                case(5); call dealloc(this%cell_nb_zmin(icell_base))
                case(6); call dealloc(this%cell_nb_zmax(icell_base))
                end select ! ibound
                nullify(nb_base)
!
             else ! base_cell_is_refined
!
                if(nb_base_cell_is_refined) then
                   if(size(nb_base)==1) then
                      ! if there is only one refined neighbouring base cell, all subcells are neighbours
                      ! of this (unrefined) base cell
                      ! -> Replace 'bound'-neighbours of this base cell by complete array idx_sub_nb_base_at_bound
                      ! remove neighbouring base cell indices (if there are external neighbouring base cells, 
                      ! their indices are contained in , so we get the correct result)
                      if(size(idx_sub_nb_base_at_bound)==1) then
                         ! directly overwrite value in respective array this%cell_nb_'bound'(icell_base)
                         nb_base(1) = idx_sub_nb_base_at_bound(1)
                      else
                         nb_base => reallocate(nb_base,size(idx_sub_nb_base_at_bound))
                         nb_base = idx_sub_nb_base_at_bound
                         ! reallocation has changed the memory location of the neighbour array and deallocated the
                         ! old memory, so need to re-associated this%cell_nb_'bound'(icell_base)
                         select case(ibound)
                         case(5); call associateVectorPointer(this%cell_nb_zmin(icell_base),nb_base)
                         case(6); call associateVectorPointer(this%cell_nb_zmax(icell_base),nb_base)
                         end select ! ibound
                      end if
                   else ! size(nb_base)==1
                      ! we need to search for overlapping 2D faces, because there can be subcells which lie laterally 
                      ! outside the face of this base cell
!
                      if(allocated(condition_c1_c2_c3)) deallocate(condition_c1_c2_c3)
                      allocate(condition_c1_c2_c3(size(idx_sub_nb_base_at_bound)))
!
                      ! call subroutine findOverlapping2DFacesChunksInversionGrid, with c3 = x, c3 = y 
                      ! (c1 is always z, as we are on a z-boundary)
                      ! There is no reverse orientation here, since in z direction we are always in the same chunks
                      ! with same orientations of x- and y- coordinates
                      call findOverlapping2DFacesChunksInversionGrid(&
                           this%ix2sorted(this%cell_ixmin(icell_base)),& ! c2 min  this is the base cell
                           this%ix2sorted(this%cell_ixmax(icell_base)),& ! c2 max
                           this%iy2sorted(this%cell_iymin(icell_base)),& ! c3 min
                           this%iy2sorted(this%cell_iymax(icell_base)),& ! c3 max
                           size(idx_sub_nb_base_at_bound),&
                           this%ix2sorted(this%cell_ixmin(idx_sub_nb_base_at_bound)),& ! c2 min array
                           this%ix2sorted(this%cell_ixmax(idx_sub_nb_base_at_bound)),& ! c2 max array
                           this%iy2sorted(this%cell_iymin(idx_sub_nb_base_at_bound)),& ! c3 min array
                           this%iy2sorted(this%cell_iymax(idx_sub_nb_base_at_bound)),& ! c3 max array
                           condition_c1_c2_c3)
                      nnb = count(condition_c1_c2_c3)
                      if(nnb>0) then
                         allocate(idx_nb(nnb))
                         idx_nb = pack( idx_sub_nb_base_at_bound , condition_c1_c2_c3 )
                      else
                         write(*,*) "ERROR 14 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: base cell ",&
                              icell_base," has no '",bound,"' neighbours in terms of neighbouring subcells, ",&
                              "even though there is at least one neighbouring base cell which are refined. ",&
                              "-> module chunksInversionGrid is inconsistent"
                         stop
                      end if
                      ! replace the 'bound'-neighbours of this base cell by idx_nb, deallocating the previous 
                      ! neighbour(s) (that were base cells only) -> only keep external neighbours for simplicity
                      select case(ibound)
                      case(5)
                         call dealloc(this%cell_nb_zmin(icell_base))
                         call associateVectorPointer(this%cell_nb_zmin(icell_base),idx_nb)
                      case(6)
                         call dealloc(this%cell_nb_zmax(icell_base))
                         call associateVectorPointer(this%cell_nb_zmax(icell_base),idx_nb)
                      end select ! ibound
                      nullify(idx_nb)
                      nullify(nb_base)
                   end if ! size(nb_base)==1
!
                else ! nb_base_cell_is_refined
                   ! the only neighbour(s) to be added would be the neighbouring base cell(s), which already is/are
                   ! set as neighbour(s), so nothing to do here
                end if ! nb_base_cell_is_refined
             end if ! base_cell_is_refined
!
          else ! associated(nb_base)
             ! since there is no 'bound'-neighbour of this base cell, there is nothing to do at this point
          end if ! associated(nb_base)
       end do ! ibound zmin,zmax
!
       ! CHECK FOR CONSISTENCY OF THIS SUBROUTINE: REFINED BASE CELLS SHOULD NOT HAVE ASSIGNED ANY NEIGHBOURS ANYMORE
       if(base_cell_is_refined) then
          if(associated(getVectorPointer(this%cell_nb_xmin(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has xmin neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
          if(associated(getVectorPointer(this%cell_nb_xmax(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has xmax neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
          if(associated(getVectorPointer(this%cell_nb_ymin(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has ymin neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
          if(associated(getVectorPointer(this%cell_nb_ymax(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has ymax neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
          if(associated(getVectorPointer(this%cell_nb_zmin(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has zmin neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
          if(associated(getVectorPointer(this%cell_nb_zmax(icell_base)))) then
             write(*,*) "ERROR 17 in findFaceNeighboursOfRefinedCellsChunksInversionGrid: after subcell neighbour ",&
                  "definition, base cell ",icell_base," still has zmax neighbours, even though it has subcells",&
                  "-> module chunksInversionGrid is inconsistent"
             stop
          end if
       end if ! base_cell_is_refined
    end do ! icell_base
!
    if(associated(idx_sub_base_at_bound)) deallocate(idx_sub_base_at_bound)
    if(associated(idx_sub_nb_base_at_bound)) deallocate(idx_sub_nb_base_at_bound)
    if(allocated(condition_c1_c2_c3)) deallocate(condition_c1_c2_c3)
!
!------------------------------------------------------------------------
  contains
!
    subroutine findOverlapping2DFacesChunksInversionGrid(&
         c2_min_testcell_sorted,c2_max_testcell_sorted,c3_min_testcell_sorted,c3_max_testcell_sorted,&
         n,c2_min_array_sorted,c2_max_array_sorted,c3_min_array_sorted,c3_max_array_sorted,&
         condition_c1_c2_c3,condition_c1,c2_reverse_order)
      integer, intent(in) :: n,c2_min_testcell_sorted,c2_max_testcell_sorted,c3_min_testcell_sorted,c3_max_testcell_sorted
      integer, dimension(n), intent(in) :: c2_min_array_sorted,c2_max_array_sorted,&
           c3_min_array_sorted,c3_max_array_sorted
      logical, dimension(n), intent(out) :: condition_c1_c2_c3
      logical, dimension(n), optional, intent(in) :: condition_c1
      integer, optional, intent(in) :: c2_reverse_order
      ! local
      logical :: cond1_present,c2_is_reverse_order
      integer :: i
!
      condition_c1_c2_c3(:) = .false.
!
      cond1_present = present(condition_c1)
      if(present(c2_reverse_order)) then
         c2_is_reverse_order = c2_reverse_order /= 0
      else
         c2_is_reverse_order = .false.
      end if
!
      if(c2_is_reverse_order) then
         do i=1,n
            if(cond1_present) then
               if(.not.condition_c1(i)) cycle
            end if
            if(c2_min_array_sorted(i) > c2_max_testcell_sorted) then
            if(c2_max_array_sorted(i) < c2_min_testcell_sorted) then
            if(c3_min_array_sorted(i) < c3_max_testcell_sorted) then
            if(c3_max_array_sorted(i) > c3_min_testcell_sorted) then
               ! this if-construct defines
               ! condition_c1_c2_c3(i) = (condition_c1(i), if present) .and. &
               !    c2_min_array_sorted(i) > c2_max_testcell_sorted (due to reverse order) .and. &
               !    c2_max_array_sorted(i) < c2_min_testcell_sorted (due to reverse order) .and. &
               !    c3_min_array_sorted(i) < c3_max_testcell_sorted .and. &
               !    c3_max_array_sorted(i) > c3_min_testcell_sorted
               ! for reasons of performance, test one condition after another instead of every time 
               ! evaluating all of them. 
               ! if all above conditions are fulfilled, set condition_c1_c2_c3(i) to .true.
               condition_c1_c2_c3(i) = .true.
            end if
            end if
            end if
            end if
         end do ! i
!
      else ! c2_is_reverse_order
!
         do i=1,n
            if(cond1_present) then
               if(.not.condition_c1(i)) cycle
            end if
            if(c2_min_array_sorted(i) < c2_max_testcell_sorted) then
            if(c2_max_array_sorted(i) > c2_min_testcell_sorted) then
            if(c3_min_array_sorted(i) < c3_max_testcell_sorted) then
            if(c3_max_array_sorted(i) > c3_min_testcell_sorted) then
               ! this if-construct defines
               ! condition_c1_c2_c3(i) = (condition_c1(i), if present) .and. &
               !    c2_min_array_sorted(i) < c2_max_testcell_sorted .and. &
               !    c2_max_array_sorted(i) > c2_min_testcell_sorted .and. &
               !    c3_min_array_sorted(i) < c3_max_testcell_sorted .and. &
               !    c3_max_array_sorted(i) > c3_min_testcell_sorted
               ! for reasons of performance, test one condition after another instead of every time 
               ! evaluating all of them. 
               ! if all above conditions are fulfilled, set condition_c1_c2_c3(i) to .true.
               condition_c1_c2_c3(i) = .true.
            end if
            end if
            end if
            end if
         end do ! i
      end if ! c2_is_reverse_order
    end subroutine findOverlapping2DFacesChunksInversionGrid
!
    subroutine getSubcellsAtBoundaryChunksInversionGrid(icell,bnd,idx_return)
      !type (chunks_inversion_grid), intent(in) :: this
      integer, intent(in) :: icell
      character(len=*) :: bnd
      integer, dimension(:), pointer, intent(out) :: idx_return
      ! local
      integer, dimension(:), pointer :: subcells,idx_bnd
      integer :: n_subcells,n_return
      logical, dimension(:), allocatable :: map
!
      nullify(idx_return)
!
      subcells => getVectorPointer(this%cell_subcells(icell))
      ! return null pointer if this cell is not refined
      if(.not.associated(subcells)) return
      n_subcells = size(subcells)
!
      select case(bnd)
      case('xmin'); idx_bnd => this%cell_ixmin
      case('xmax'); idx_bnd => this%cell_ixmax
      case('ymin'); idx_bnd => this%cell_iymin
      case('ymax'); idx_bnd => this%cell_iymax
      case('zmin'); idx_bnd => this%cell_izmin
      case('zmax'); idx_bnd => this%cell_izmax
      case default
         write(*,*) "ERROR 15 in getSubcellsAtBoundaryChunksInversionGrid: incoming boundary code = '",bnd,&
              "', must be one of 'xmin','xmax','ymin','ymax','zmin','zmax'"
         stop
      end select
!
      allocate(map(n_subcells))
      map = (idx_bnd(subcells) == idx_bnd(icell))
      n_return = count(map)
!
      if(n_return == 0) then
         deallocate(map)
         return
      else
         allocate(idx_return(n_return))
         idx_return = pack( subcells , map )
         deallocate(map)
      end if
    end subroutine getSubcellsAtBoundaryChunksInversionGrid
  end subroutine findFaceNeighboursOfRefinedCellsChunksInversionGrid
!--------------------------------------------------------------------
!> \brief compute cell radii of all cells
  subroutine computeCellRadiiChunksInversionGrid(this)
    type (chunks_inversion_grid) :: this
    ! local
    double precision :: xmin,xmax,ymin,ymax,zmin,zmax,rmin,rmax,dx,dy,dz
    integer :: icell
!
    allocate(this%cell_radius(this%ncell_internal))
!
    do icell = 1,this%ncell_internal
       ! get the coordinates of two corners of the cell: (xmin,ymin,zmin) , (xmax,ymax,zmax)
       xmin = this%x(this%cell_ixmin(icell))
       xmax = this%x(this%cell_ixmax(icell))
       ymin = this%y(this%cell_iymin(icell))
       ymax = this%y(this%cell_iymax(icell))
       zmin = this%z(this%cell_izmin(icell))
       zmax = this%z(this%cell_izmax(icell))
!
       ! transform to CURVED local chunk; compare "FIRST STEP" in routine transformVectorFlatToGlobalChunksInversionGrid
       ! at the same time, buid the difference vector (dx,dy,dz) pointing from one corner to the other (in "curved" chunk frame)
!
       rmin = dsqrt(xmin**2+ymin**2+1.d0)
       rmax = dsqrt(xmax**2+ymax**2+1.d0)
!
       dx = xmax*zmax/rmax - xmin*zmin/rmin
       dy = ymax*zmax/rmax - ymin*zmin/rmin
       dz = zmax/rmax - zmin/rmin
!
       ! finally, the cell radius is the half of the length of the diagonal vector (pointing from the min to the max corner in "curved" chunk frame)
       this%cell_radius(icell) = 0.5d0 * dsqrt(dx**2+dy**2+dz**2)
    end do ! icell
  end subroutine computeCellRadiiChunksInversionGrid
!--------------------------------------------------------------------
!> \brief read chunks inversion grid file which was created before
  subroutine readChunksInversionGrid(this,path,inpar,lu,errmsg)
    type (chunks_inversion_grid) :: this
    character(len=*) :: path
    type (input_parameter) :: inpar
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=600) :: filename
    character(len=23) :: myname = 'readChunksInversionGrid'
    integer :: ios,ichunk,icell,size_ip,n
    integer, dimension(:), pointer :: ip
    ! par file
    integer :: nchunk,nblock
    real :: clat,clon,wlat,wlon,rmax,az_rot
    integer, dimension(:), pointer :: nlay_block,nlat_block,nlon_block
    real, dimension(:), pointer :: thickness,cref_param
!
    nullify(ip,nlay_block,nlat_block,nlon_block,thickness,cref_param)
    call addTrace(errmsg,myname)
!
    this%is_defined = .false.
!
    ! GET VALUES FROM PARAMETER FILE
!
    call readCheckParFileChunksInversionGrid(inpar,nchunk,clat,clon,wlat,wlon,rmax,az_rot,nblock,nlay_block,thickness,&
         nlat_block,nlon_block,cref_param,this%vtk_geometry_type_int,this%apply_vtk_coords_scaling_factor,&
         this%vtk_coords_scaling_factor,errmsg,myname)
    if(.level.errmsg == 2) goto 2
!
    this%vtk_projection = trim(inpar.sval.'VTK_PROJECTION')
!
    filename = trim(path)//trim(inpar.sval.'CHUNKS_INVGRID_FILE')
!
    open(unit=lu,file=filename,form='UNFORMATTED',access='STREAM',status='OLD',action='READ',iostat=ios)
    if(ios/=0) then
       call add(errmsg,2,"could not open binary file '"//trim(filename)//"'to read",myname)
    else
       call add(errmsg,0,"opened binary file '"//trim(filename)//"' to read",myname)
    end if
!
    ! read the content of the inversion grid file and check consistency with the current state of the invgrid parfile
    ! issue error messages, telling the values in the binary file
    ! this tries to assure the user to know what he is doing
!
    read(lu,iostat=ios) this%nchunk,this%clat,this%clon,this%wlat,this%wlon,this%rmax,this%az_rot
    if(ios/=0) then
       write(errstr,*) "could not read nchunk,clat,clon,wlat,wlon,rmax,az_rot: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    if(this%nchunk /= nchunk) then
       write(errstr,*) "nchunk = ",this%nchunk," in inversion grid file is inconsistent with its value in parfile = ",&
            nchunk,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%clat /= clat) then
       write(errstr,*) "clat = ",this%clat," in inversion grid file is inconsistent with its value in parfile = ",&
            clat,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%clon /= clon) then
       write(errstr,*) "clon = ",this%clon," in inversion grid file is inconsistent with its value in parfile = ",&
            clon,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%wlat /= wlat) then
       write(errstr,*) "wlat = ",this%wlat," in inversion grid file is inconsistent with its value in parfile = ",&
            wlat,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%wlon /= wlon) then
       write(errstr,*) "wlon = ",this%wlon," in inversion grid file is inconsistent with its value in parfile = ",&
            wlon,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%rmax /= rmax) then
       write(errstr,*) "rmax = ",this%rmax," in inversion grid file is inconsistent with its value in parfile = ",&
            rmax,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%az_rot /= az_rot) then
       write(errstr,*) "az_rot = ",this%az_rot," in inversion grid file is inconsistent with its value in parfile = ",&
            az_rot,"; you may want to recreate the inversion grid file ",&
            "(e.g. by running initBasics -recr or (re)moving the binary inversion grid file"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
   call computeTransformsChunksInversionGrid(this)
!
    read(lu,iostat=ios) this%nx
    if(ios/=0) then
       write(errstr,*) "could not read nx: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%nx < 1) then
       write(errstr,*) "value of nx in file = ",this%nx,"; must be strictly positive, file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%x(this%nx))
    read(lu,iostat=ios) this%x
    if(ios/=0) then
       write(errstr,*) "could not read x-coordinates array (nx = ",this%nx,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%sorted2ix(this%nx))
    read(lu,iostat=ios) this%sorted2ix
    if(ios/=0) then
       write(errstr,*) "could not read sorting mapping of x-coordinates (nx = ",this%nx,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%ix2sorted(this%nx))
    read(lu,iostat=ios) this%ix2sorted
    if(ios/=0) then
       write(errstr,*) "could not read inverse sorting mapping of x-coordinates (nx = ",this%nx,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!    
    read(lu,iostat=ios) this%ny
    if(ios/=0) then
       write(errstr,*) "could not read ny: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%ny < 1) then
       write(errstr,*) "value of ny in file = ",this%ny,"; must be strictly positive, file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%y(this%ny))
    read(lu,iostat=ios) this%y
    if(ios/=0) then
       write(errstr,*) "could not read y-coordinates array (ny = ",this%ny,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%sorted2iy(this%ny))
    read(lu,iostat=ios) this%sorted2iy
    if(ios/=0) then
       write(errstr,*) "could not read sorting mapping of y-coordinates (ny = ",this%ny,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%iy2sorted(this%ny))
    read(lu,iostat=ios) this%iy2sorted
    if(ios/=0) then
       write(errstr,*) "could not read inverse sorting mapping of y-coordinates (ny = ",this%ny,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    read(lu,iostat=ios) this%nz
    if(ios/=0) then
       write(errstr,*) "could not read nz: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%nz < 1) then
       write(errstr,*) "value of nz in file = ",this%nz,"; must be strictly positive, file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%z(this%nz))
    read(lu,iostat=ios) this%z
    if(ios/=0) then
       write(errstr,*) "could not read z-coordinates array (nz = ",this%nz,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%sorted2iz(this%nz))
    read(lu,iostat=ios) this%sorted2iz
    if(ios/=0) then
       write(errstr,*) "could not read sorting mapping of z-coordinates (nz = ",this%nz,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%iz2sorted(this%nz))
    read(lu,iostat=ios) this%iz2sorted
    if(ios/=0) then
       write(errstr,*) "could not read inverse sorting mapping of z-coordinates (nz = ",this%nz,"): raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    read(lu,iostat=ios) this%ncell,this%ncell_internal,this%ncell_base
    if(ios/=0) then
       write(errstr,*) "could not read ncell,ncell_internal,ncell_base: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%ncell < 1) then
       write(errstr,*) "value of ncell in file = ",this%ncell,"; must be strictly positive, file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(this%ncell_internal < 1) then
       write(errstr,*) "value of ncell_internal in file = ",this%ncell_internal,&
            "; must be strictly positive, file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    allocate(this%icell_internal(this%ncell),this%icell_external(this%ncell_internal))
    read(lu,iostat=ios) this%icell_internal,this%icell_external
    if(ios/=0) then
       write(errstr,*) "could not read arrays icell_internal,icell_external: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(any(this%icell_internal < 1 .or. this%icell_internal > this%ncell_internal)) then
       write(errstr,*) count(this%icell_internal < 1 .or. this%icell_internal > this%ncell_internal),&
            " values of array icell_internal are either < 1 , or > ncell_internal = ",this%ncell_internal,&
            "; file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(count(this%icell_external /= -1) /= this%ncell) then
       write(errstr,*) count(this%icell_internal /= -1)," values of array icell_external are different from -1. ",&
            "however, this number must equal ncell = ",this%ncell,"; file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(.not.all(this%icell_external == -1 .or. (this%icell_external >= 1 .and. this%icell_external <= this%ncell))) then
       write(errstr,*) "there are values of array icell_external, which are neither -1 nor in permitted range ",&
            "between 1 and ncell = ",this%ncell,"; file may be corrupt"
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    ! could further check if the arrays icell_internal and icell_external are consistent, e.g. by checking whether
    ! the mappings are injective (except value -1 in array icell_external occurring (ncell_internal-ncell+1)-times)
    ! and icell_external(icell_internal(i)) == i for all 1<=i<=ncell
!
    allocate(this%cell_ichunk(this%ncell_internal),&
         this%cell_ixmin(this%ncell_internal),this%cell_ixmax(this%ncell_internal),&
         this%cell_iymin(this%ncell_internal),this%cell_iymax(this%ncell_internal),&
         this%cell_izmin(this%ncell_internal),this%cell_izmax(this%ncell_internal))
    read(lu,iostat=ios) this%cell_ichunk,this%cell_ixmin,this%cell_ixmax,this%cell_iymin,this%cell_iymax,this%cell_izmin,&
         this%cell_izmax
    if(ios/=0) then
       write(errstr,*) "could not read arrays cell_ichunk,cell_ixmin,cell_ixmax,cell_iymin,cell_iymax,cell_izmin: ",&
            "raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    allocate(this%chunk_ncell_internal(this%nchunk),this%chunk_icell_internal(this%nchunk))
    read(lu,iostat=ios) this%chunk_ncell_internal
    if(ios/=0) then
       write(errstr,*) "could not read array chunk_ncell_internal: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    if(any(this%chunk_ncell_internal < 0)) then
       write(*,*) "ERROR 16 in readChunksInversionGrid: there are negative values in vector chunk_ncell_internal = ",&
            this%chunk_ncell_internal
       call add(errmsg,2,"there are negative values in vector chunk_ncell_internal",myname)
       goto 2
    end if
    do ichunk = 1,this%nchunk
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of chunk_icell_internal of chunk ",ichunk," (out of ",this%nchunk,&
               "): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip /= this%chunk_ncell_internal(ichunk)) then
          write(errstr,*) "size of chunk_icell_internal of chunk ",ichunk," (out of ",this%nchunk,&
               ") equals ",size_ip,". This is inconsistent with the expected (or erroneous) value chunk_ncell_internal = ",&
               this%chunk_ncell_internal(ichunk)
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not vector chunk_icell_internal of chunk ",ichunk," (out of ",this%nchunk,&
               "), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%chunk_icell_internal(ichunk),ip)
       nullify(ip)
    end do ! ichunk
!
    allocate(this%cell_center(3,this%ncell_internal),this%cell_radius(this%ncell_internal))
    read(lu,iostat=ios) this%cell_center,this%cell_radius
    if(ios/=0) then
       write(errstr,*) "could not read arrays cell_center,cell_radius: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
!
    allocate(this%cell_parent(this%ncell_internal),this%cell_subcells(this%ncell_internal))
    read(lu,iostat=ios) this%cell_parent
    if(ios/=0) then
       write(errstr,*) "could not read array cell_parent: raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of subcells-vector of cell ",icell," (out of ",this%ncell_internal,&
               "): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 2) then
          write(errstr,*) "size of subcells-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 2; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read subcells-vector of cell ",icell," (out of ",this%ncell_internal,&
               "), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_subcells(icell),ip)
       nullify(ip)
    end do ! icell
!
    ! xmin neighbours
    allocate(this%cell_nb_xmin(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of xmin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of xmin-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read xmin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_xmin(icell),ip)
       nullify(ip)
    end do ! icell
    ! xmax neighbours
    allocate(this%cell_nb_xmax(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of xmax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of xmax-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read xmax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_xmax(icell),ip)
       nullify(ip)
    end do ! icell
    ! ymin neighbours
    allocate(this%cell_nb_ymin(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of ymin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of ymin-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read ymin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_ymin(icell),ip)
       nullify(ip)
    end do ! icell
    ! ymax neighbours
    allocate(this%cell_nb_ymax(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of ymax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of ymax-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read ymax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_ymax(icell),ip)
       nullify(ip)
    end do ! icell
    ! zmin neighbours
    allocate(this%cell_nb_zmin(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of zmin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of zmin-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read zmin-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_zmin(icell),ip)
       nullify(ip)
    end do ! icell
    ! zmax neighbours
    allocate(this%cell_nb_zmax(this%ncell_internal))
    do icell = 1,this%ncell_internal
       read(lu,iostat=ios) size_ip
       if(ios/=0) then
          write(errstr,*) "could not read size of zmax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"): raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(size_ip == 0) cycle
       if(size_ip < 1) then
          write(errstr,*) "size of zmax-neighbours-vector of cell ",icell," (out of ",this%ncell_internal,&
               ") has value ",size_ip,": must be at least 1; file may be corrupt"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       allocate(ip(size_ip))
       read(lu,iostat=ios) ip
       if(ios/=0) then
          write(errstr,*) "could not read zmax-neighbours-vector of cell ",icell," (out of ",&
               this%ncell_internal,"), which has size ",size_ip,": raised iostat = ",ios
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       call associateVectorPointer(this%cell_nb_zmax(icell),ip)
       nullify(ip)
    end do ! icell
    ! contains_refeined_cells
    read(lu,iostat=ios) n
    if(ios/=0) then
       write(errstr,*) "could not read integer indicating the logical 'contains_refined_cells': raised iostat = ",ios
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    select case(n)
    case (0)
       this%contains_refined_cells = .false.
    case (1)
       this%contains_refined_cells = .true.
    case default
       write(errstr,*) "integer indicating the logical 'contains_refined_cells' has value ",n,": expect either 0 or 1"
       call add(errmsg,2,errstr,myname)
       goto 2
    end select
!
    ! if code comes here, everything was read in and is consistent; so, inversion grid is correctly defined
    this%is_defined = .true.
!
1   close(lu)
    if(associated(nlay_block)) deallocate(nlay_block)
    if(associated(nlat_block)) deallocate(nlat_block)
    if(associated(nlon_block)) deallocate(nlon_block)
    if(associated(thickness)) deallocate(thickness)
    return
!
    ! in case of an error, deallocate everything so far defined before returning
2   call deallocateChunksInversionGrid(this)
    goto 1
  end subroutine readChunksInversionGrid
!--------------------------------------------------------------------
!> \brief write chunks inversion grid file
  subroutine writeChunksInversionGrid(this,filename,lu,open_status,errmsg)
    type (chunks_inversion_grid) :: this
    character(len=*) :: filename,open_status
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=24) :: myname = 'writeChunksInversionGrid'
    integer :: ios,ichunk,icell,zero,one
    integer, dimension(:), pointer :: ip
!
    nullify(ip)
    call addTrace(errmsg,myname)
!
    zero = 0
    one = 1
!
    open(unit=lu,file=filename,form='UNFORMATTED',access='STREAM',status=trim(open_status),action='WRITE',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "could not open binary file '"//trim(filename)//"'to write, raised iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    write(lu) this%nchunk,this%clat,this%clon,this%wlat,this%wlon,this%rmax,this%az_rot
    write(lu) this%nx,this%x,this%sorted2ix,this%ix2sorted
    write(lu) this%ny,this%y,this%sorted2iy,this%iy2sorted
    write(lu) this%nz,this%z,this%sorted2iz,this%iz2sorted
    write(lu) this%ncell,this%ncell_internal,this%ncell_base
    write(lu) this%icell_internal,this%icell_external
    write(lu) this%cell_ichunk,this%cell_ixmin,this%cell_ixmax,this%cell_iymin,this%cell_iymax,this%cell_izmin,&
         this%cell_izmax
    write(lu) this%chunk_ncell_internal
    do ichunk=1,this%nchunk
       ip => getVectorPointer(this%chunk_icell_internal(ichunk))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! ichunk
    write(lu) this%cell_center,this%cell_radius
    write(lu) this%cell_parent
    !subcells
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_subcells(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! xmin neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_xmin(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! xmax neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_xmax(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! ymin neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_ymin(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! ymax neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_ymax(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! zmin neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_zmin(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    ! zmax neighbours
    do icell = 1,this%ncell_internal
       ip => getVectorPointer(this%cell_nb_zmax(icell))
       if(associated(ip)) then
          write(lu) size(ip),ip
       else
          write(lu) zero
       end if
    end do ! icell
    if(this%contains_refined_cells) then
       write(lu) one
    else
       write(lu) zero
    end if
    close(lu)
  end subroutine writeChunksInversionGrid
!--------------------------------------------------------------------
!> \brief deallocate chunks inversion grid
!
  subroutine deallocateChunksInversionGrid(this)
    type (chunks_inversion_grid) :: this
    integer :: ichunk,icell
    this%nchunk = 0
    this%clat = 0.; this%clon = 0.
    this%wlat = 0.; this%wlon = 0.
    this%rmax = 0.
    this%az_rot = 0.
    if(associated(this%Mrot_local2flat)) deallocate(this%Mrot_local2flat)
    if(associated(this%Mrot_global2flat)) deallocate(this%Mrot_global2flat)
    if(associated(this%vtk_flat_shift)) deallocate(this%vtk_flat_shift)
    this%ncell = 0; this%ncell_internal = 0; this%ncell_base = 0
    if(associated(this%icell_internal)) deallocate(this%icell_internal)
    if(associated(this%icell_external)) deallocate(this%icell_external)
    if(associated(this%x)) deallocate(this%x)
    if(associated(this%sorted2ix)) deallocate(this%sorted2ix)
    if(associated(this%ix2sorted)) deallocate(this%ix2sorted)
    this%nx = 0
    if(associated(this%y)) deallocate(this%y)
    if(associated(this%sorted2iy)) deallocate(this%sorted2iy)
    if(associated(this%iy2sorted)) deallocate(this%iy2sorted)
    this%ny = 0
    if(associated(this%z)) deallocate(this%z)
    if(associated(this%sorted2iz)) deallocate(this%sorted2iz)
    if(associated(this%iz2sorted)) deallocate(this%iz2sorted)
    this%nz = 0
    if(associated(this%cell_ichunk)) deallocate(this%cell_ichunk)
    if(associated(this%chunk_ncell_internal)) deallocate(this%chunk_ncell_internal)
    if(associated(this%chunk_icell_internal)) then
       do ichunk = 1,size(this%chunk_icell_internal)
          call dealloc(this%chunk_icell_internal(ichunk))
       end do ! ichunk
       deallocate(this%chunk_icell_internal)
    end if
    if(associated(this%cell_ixmin)) deallocate(this%cell_ixmin)
    if(associated(this%cell_ixmax)) deallocate(this%cell_ixmax)
    if(associated(this%cell_iymin)) deallocate(this%cell_iymin)
    if(associated(this%cell_iymax)) deallocate(this%cell_iymax)
    if(associated(this%cell_izmin)) deallocate(this%cell_izmin)
    if(associated(this%cell_izmax)) deallocate(this%cell_izmax)
    if(associated(this%cell_center)) deallocate(this%cell_center)
    if(associated(this%cell_radius)) deallocate(this%cell_radius)
    if(associated(this%cell_nb_xmin)) then
       do icell = 1,size(this%cell_nb_xmin)
          call dealloc(this%cell_nb_xmin(icell))
       end do ! icell
       deallocate(this%cell_nb_xmin)
    end if
    if(associated(this%cell_nb_xmax)) then
       do icell = 1,size(this%cell_nb_xmax)
          call dealloc(this%cell_nb_xmax(icell))
       end do ! icell
       deallocate(this%cell_nb_xmax)
    end if
    if(associated(this%cell_nb_ymin)) then
       do icell = 1,size(this%cell_nb_ymin)
          call dealloc(this%cell_nb_ymin(icell))
       end do ! icell
       deallocate(this%cell_nb_ymin)
    end if
    if(associated(this%cell_nb_ymax)) then
       do icell = 1,size(this%cell_nb_ymax)
          call dealloc(this%cell_nb_ymax(icell))
       end do ! icell
       deallocate(this%cell_nb_ymax)
    end if
    if(associated(this%cell_nb_zmin)) then
       do icell = 1,size(this%cell_nb_zmin)
          call dealloc(this%cell_nb_zmin(icell))
       end do ! icell
       deallocate(this%cell_nb_zmin)
    end if
    if(associated(this%cell_nb_zmax)) then
       do icell = 1,size(this%cell_nb_zmax)
          call dealloc(this%cell_nb_zmax(icell))
       end do ! icell
       deallocate(this%cell_nb_zmax)
    end if
    if(associated(this%cell_parent)) deallocate(this%cell_parent)
    if(associated(this%cell_subcells)) then
       do icell = 1,size(this%cell_subcells)
          call dealloc(this%cell_subcells(icell))
       end do ! icell
       deallocate(this%cell_subcells)
    end if
    this%vtk_geometry_type_int = - 1
    this%is_defined = .false.
  end subroutine deallocateChunksInversionGrid
!------------------------------------------------------------------------
!> \brief return overall number of invgrid cells, if any
!
  function getNcellChunksInversionGrid(this) result(ncell)
    type (chunks_inversion_grid), intent(in) :: this
    integer :: ncell
    if(this%is_defined) then
       ncell = this%ncell
    else
       ncell = 0
    end if
  end function getNcellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief for a given chunk, transform flat tangential plane (and depth) coords to global cartesian coordinates
!! \details First project planar coordinates down to spherical chunk (centered on North pole, no azimuthal rotation).
!!  Then rotate the whole chunk to the original center and apply azimutal rotation (by inverse of matrices Mrot_global2flat)
!! \param this chunks inversion grid
!! \param x1 vector of x coordinate
!! \param x2 vector of y coordinate
!! \param x3 vector of z coordinate
!! \param idx_chunk vector of chunk indices, defining for each point in which chunk it lies (dependent on ichunk, the correct rotation matrix Mrot_global2flat is chosen)
!! \param n length of vectors x1,x2,x3
  subroutine transformVectorFlatToGlobalChunksInversionGrid(this,x1,x2,x3,idx_chunk,n)
    type (chunks_inversion_grid) :: this
    integer :: n
    double precision, dimension(n) :: x1,x2,x3
    integer, dimension(n) :: idx_chunk
    ! local
    double precision, dimension(n) :: r,x1_tmp,x2_tmp,x3_tmp
!
    ! ###  FIRST STEP: GENERATE CURVATURE OF SPHERICAL CHUNK  ###
!
    ! INTUITIVE MATH
    ! Note that coordinates in this%x,this%y correspond to tangential plane coords w.r.t. a sphere of radius 1, 
    ! and x3 contains radius r (in unit of distance, e.g. m or km) defining the actual spherical radius of that 
    ! point inside Earth
    !
    ! 1) scale x1 and x2 by this%rmax. Then x1,x2,x3 constitute global Cartesian coordinates in a cuboidal
    !    chunk (not spherical chunk, i.e. no curvature at all, but in units of distance)
    !      x1 = x1*this%rmax
    !      x2 = x2*this%rmax
    !
    ! 2) transform cuboidal chunk to a "cone" with a flat surface (same surface as cuboid's surface), but 
    !    where x and y coordinates correspond to coordinates in a tangential plane w.r.t. a sphere of radius x3
    !      x1 = x1*(x3/this%rmax)
    !      x2 = x2*(x3/this%rmax)
    !
    ! 3) project the points in this "cone" down to the surface of a sphere of radius x3 in order to form
    !    a curved spherical chunk. therefore, compute the norm r of the vector (x1,x2,x3), normalize the vector
    !    by its norm r and then scale it by x3
    !      r = sqrt(x1**2 + x2**2 + x3**2)
    !      x1 = x1*(x3/r)
    !      x2 = x2*(x3/r)
    !      x3 = x3*(x3/r)
!
    ! NOW AN  E Q U I V A L E N T  COMPUTATION WHICH IS MUCH MORE PERFORMANT
    ! (accounuting for everything cancelling out in the above steps, hence avoiding redundant computations)
    ! The euklidean distance sqrt(x1**2+x2**2+1.d0) w.r.t. sphere of radius 1 (below) might numerically be 
    ! much more stable than r = sqrt(x1**2 + x2**2 + x3**2) (above) for coords in distance unit (e.g. for coords 
    ! in m the values get very very large, so the numbers might become hard to represent even in double precision?!)
!
    r = dsqrt(x1**2+x2**2+1.d0) 
    x1 = x1*x3/r
    x2 = x2*x3/r
    x3 = x3/r
!
!
    ! ###  SECOND STEP: APPLY ROTATION TO ORIGINAL CHUNK CENTER, INCLUDING AZIMUTHAL ROTATION  ###
!
    ! local copy of coords in unrotated spherical chunk
    x1_tmp = x1; x2_tmp = x2; x3_tmp = x3
!
    ! apply a rotation which is inverse to Mrot_global2flat by using the transposed matrix
    x1 = this%Mrot_global2flat(1,1,idx_chunk)*x1_tmp+this%Mrot_global2flat(2,1,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(3,1,idx_chunk)*x3_tmp
    x2 = this%Mrot_global2flat(1,2,idx_chunk)*x1_tmp+this%Mrot_global2flat(2,2,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(3,2,idx_chunk)*x3_tmp
    x3 = this%Mrot_global2flat(1,3,idx_chunk)*x1_tmp+this%Mrot_global2flat(2,3,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(3,3,idx_chunk)*x3_tmp
  end subroutine transformVectorFlatToGlobalChunksInversionGrid
!------------------------------------------------------------------------
!> \brief for a given chunk, transform global cartesian coordinates to flat tangential plane (and depth) coords
!! \details First project planar coordinates down to spherical chunk (centered on North pole, no azimuthal rotation).
!!  Then rotate the whole chunk to the original center and apply azimutal rotation (by inverse of matrices Mrot_global2flat)
!! \param this chunks inversion grid
!! \param x1 vector of x coordinate
!! \param x2 vector of y coordinate
!! \param x3 vector of z coordinate
!! \param idx_chunk vector of chunk indices, defining for each point in which chunk it lies (dependent on ichunk, the correct rotation matrix Mrot_global2flat is chosen)
!! \param n length of vectors x1,x2,x3
  subroutine transformVectorGlobalToFlatChunksInversionGrid(this,x1,x2,x3,idx_chunk,n)
    type (chunks_inversion_grid) :: this
    integer :: n
    double precision, dimension(n) :: x1,x2,x3
    integer, dimension(n) :: idx_chunk
    ! local
    double precision, dimension(n) :: x1_tmp,x2_tmp,x3_tmp
!
    ! local copy of global Cartesian coords
    x1_tmp = x1; x2_tmp = x2; x3_tmp = x3
!
    ! apply rotation from global Cartesian to flat reference chunk (individual rotation per chunk)
    x1 = this%Mrot_global2flat(1,1,idx_chunk)*x1_tmp+this%Mrot_global2flat(1,2,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(1,3,idx_chunk)*x3_tmp
    x2 = this%Mrot_global2flat(2,1,idx_chunk)*x1_tmp+this%Mrot_global2flat(2,2,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(2,3,idx_chunk)*x3_tmp
    x3 = this%Mrot_global2flat(3,1,idx_chunk)*x1_tmp+this%Mrot_global2flat(3,2,idx_chunk)*x2_tmp+&
         this%Mrot_global2flat(3,3,idx_chunk)*x3_tmp
!
    ! project the lateral x1,x2 coords to a tangential plane corresponding to radius 1.0 (not possible if z=0, 
    ! but this would correspond either to the center of the earth or a chunk of infinite width, which both should
    ! not occurr)
    where(x3/= 0.d0)
       x1 = x1/x3
       x2 = x2/x3
    elsewhere
       x1 = 0.d0
       x2 = 0.d0
    end where
!
    ! the new z coordinate is the radius of the original point.
    ! since x1 -> x1/x3  and x2 -> x2/x3, you can rewrite the radius by pulling out x3 from the sqrt:
    x3 = x3*dsqrt(x1**2+x2**2+1.d0)
  end subroutine transformVectorGlobalToFlatChunksInversionGrid
!------------------------------------------------------------------------
!> \brief internal transform from flat inversion grid coords to vtk projection
!! \details for all "FLAT_*" projection types: do not return the projections to
!!  the tangential plane, but apply atan/(0.25 pi) first in order to produce equiDISTANT
!!  FLAT plots from equiANGULAR cell distributions. this is only done for vtk plotting
!!  and does not affect the general equiangular distribution of the cells
  subroutine transformFlatToVtkChunksInversionGrid(this,x1,x2,x3,idx_chunk,n)
    type (chunks_inversion_grid) :: this
    integer :: n
    double precision, dimension(n) :: x1,x2,x3
    integer, dimension(n) :: idx_chunk
    ! local
    double precision, dimension(n) :: r,x1_tmp,x2_tmp,x3_tmp
!
    select case(this%vtk_projection)
    case('GLOBAL')
       call transformVectorFlatToGlobalChunksInversionGrid(this,x1,x2,x3,idx_chunk,n)
    case('LOCAL_CURV','LOCAL_NORTH_CURV')
       ! curv down the flat reference chunks
       r = dsqrt(x1**2+x2**2+1.d0) 
       x1 = x1*x3/r
       x2 = x2*x3/r
       x3 = x3/r
!
       ! rotate the flat reference chunks to local position (by inverse (transposed) rotation local2flat^T)
       x1_tmp = x1; x2_tmp = x2; x3_tmp = x3
       x1 = this%Mrot_local2flat(1,1,idx_chunk)*x1_tmp+this%Mrot_local2flat(2,1,idx_chunk)*x2_tmp+&
            this%Mrot_local2flat(3,1,idx_chunk)*x3_tmp
       x2 = this%Mrot_local2flat(1,2,idx_chunk)*x1_tmp+this%Mrot_local2flat(2,2,idx_chunk)*x2_tmp+&
            this%Mrot_local2flat(3,2,idx_chunk)*x3_tmp
       x3 = this%Mrot_local2flat(1,3,idx_chunk)*x1_tmp+this%Mrot_local2flat(2,3,idx_chunk)*x2_tmp+&
            this%Mrot_local2flat(3,3,idx_chunk)*x3_tmp
!
       ! for LOCAL_NORTH_CURV, additionally apply the azimuthal rotation to the x- and y- coordinates
       if(this%vtk_projection == 'LOCAL_NORTH_CURV') then
          x1_tmp = x1; x2_tmp = x2
          x1 = this%Mrot_azimuth(1,1)*x1_tmp+this%Mrot_azimuth(1,2)*x2_tmp
          x2 = this%Mrot_azimuth(2,1)*x1_tmp+this%Mrot_azimuth(2,2)*x2_tmp
       end if
    case('LOCAL_FLAT','LOCAL_NORTH_FLAT')
       ! first transform lateral coordinates of flat reference chunk (tangential plane coordinates)
       ! to arclength w.r.t to the radius at the surface. This is done to display an equidistant 
       ! cell distribution instead of the stretched tangential plane projections
       x1 = atan(x1)*dble(this%rmax)
       x2 = atan(x2)*dble(this%rmax)
!
       ! then shift the individual flat reference chunks to their correct lateral position
       x1 = x1 + this%vtk_flat_shift(1,idx_chunk)
       x2 = x2 + this%vtk_flat_shift(2,idx_chunk)
!
       ! for LOCAL_NORTH_FLAT, additionally apply the azimuthal rotation to the x- and y- coordinates
       if(this%vtk_projection == 'LOCAL_NORTH_FLAT') then
          x1_tmp = x1; x2_tmp = x2
          x1 = this%Mrot_azimuth(1,1)*x1_tmp+this%Mrot_azimuth(1,2)*x2_tmp
          x2 = this%Mrot_azimuth(2,1)*x1_tmp+this%Mrot_azimuth(2,2)*x2_tmp
       end if
    end select
!
    if(this%apply_vtk_coords_scaling_factor) then
       x1 = x1 * this%vtk_coords_scaling_factor
       x2 = x2 * this%vtk_coords_scaling_factor
       x3 = x3 * this%vtk_coords_scaling_factor
    end if
  end subroutine transformFlatToVtkChunksInversionGrid
!------------------------------------------------------------------------
!> \brief transform wp,event or station coords to x,y,z coords for vtk application
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station' and that c1,c2,c3 are associated and have all same length.
!!  Also it can be assured that this chunks inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here).
!!  in case coords_type is 'wp', c1,c2,c3 are expected to be global cartesian x,y,z coordinates for spherical application
!!  in case coords_type is 'event', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and depth (in km)
!!  in case coords_type is 'station', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and altitude (in m)
!!  altitude of stations is ignored on return (in order to plot stations on surface of spherical chunk)
!! \param this chunks inversion grid
!! \param c1 vector or first coordinate (contains vtk x-values on exit)
!! \param c2 vector or second coordinate (contains vtk y-values on exit)
!! \param c3 vector or third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!! \param errmsg error message
!! \param uf_wp optional unit of wavefield points (recommended to be used along with coords_type == 'wp'). If not present, wavefield points are assumed to be in km
!
  subroutine transformToVtkChunksInversionGrid(this,c1,c2,c3,coords_type,errmsg,uf_wp)
    type (chunks_inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
    real, optional :: uf_wp
    ! local
    character(len=33) :: myname = 'transformToVtkChunksInversionGrid'
    character(len=400) :: errstr
    integer :: np
    double precision, dimension(:), allocatable :: c1d,c2d,c3d
    double precision, dimension(:), allocatable :: lat,lon,z
    double precision, dimension(:), allocatable :: x1_tmp,x2_tmp,x3_tmp
    integer, dimension(:), allocatable :: idx_chunk
    double precision :: uf_correction
!
    call addTrace(errmsg,myname)
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    np = size(c1)
!
    ! for internal use, transform to double precision first
    allocate(c1d(np),c2d(np),c3d(np))
    c1d = dble(c1); c2d = dble(c2); c3d = dble(c3)
!
    select case(coords_type)
    case('station','event')
       ! when using spherical inversion grids, the event and station coordinates
       ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
       ! so in case of coords_type = 'event' or 'station', transform the incoming c1,c2,c3 to respective
       ! global Cartesian coordinates first, then execute the transformations below
!
       allocate(lat(np),lon(np),z(np))
       lat = c1d; lon = c2d; z = c3d
!
       ! first define unit vectors (c1d,c2d,c3d) in global Cartesian coordinates pointing to the direction of the event/station location
       c1d = cos(lat*mc_deg2radd)*cos(lon*mc_deg2radd)
       c2d = cos(lat*mc_deg2radd)*sin(lon*mc_deg2radd)
       c3d = sin(lat*mc_deg2radd)
       ! then scale the unit vectors (c1d,c2d,c3d) according to the incoming third coordinate z
       select case(coords_type)
       case('station')
          ! for station coordinates, the incoming 3rd coordinate is altitude in meters, which here is ignored: so scale to surface of the earth
          c1d = c1d * dble(this%rmax)
          c2d = c2d * dble(this%rmax)
          c3d = c3d * dble(this%rmax)
       case('event')
          ! for event coordinates, the incoming 3rd coordinate is depth in km
          c1d = c1d * (dble(this%rmax)-z)
          c2d = c2d * (dble(this%rmax)-z)
          c3d = c3d * (dble(this%rmax)-z)
       end select
!
       deallocate(lat,lon,z)
!
    case('wp')
       ! if the optional argument uf_wp is present, check if wavefield point coordinates are given in km
       ! and if necessary transform to km first before executing the transformations below
       ! if uf_wp is not present, IT IS ASSUMED THAT WAVEFIELD POINT COORDINATES ARE IN km !
       if(present(uf_wp)) then
          if(uf_wp /= 1.0e3) then
             ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
             uf_correction = dble(uf_wp) * 1.0d-3
             c1d = c1d * uf_correction
             c2d = c2d * uf_correction
             c3d = c3d * uf_correction
          end if ! uf_wp /= 1.0e3
       end if ! present(uf_wp)
    end select
!
    ! at this point, variables c1d,c2d,c3d are global Cartesian x,y,z coordinates
    ! in case of this%vtk_projection == 'GLOBAL', no further transformations are required except scaling
    ! for the other cases, apply further transformations to local chunk or flat scart coordinates, respectively
!
    select case(this%vtk_projection)
    case('GLOBAL')
       if(this%apply_vtk_coords_scaling_factor) then
          c1d = c1d * this%vtk_coords_scaling_factor
          c2d = c2d * this%vtk_coords_scaling_factor
          c3d = c3d * this%vtk_coords_scaling_factor
       end if
    case('LOCAL_CURV','LOCAL_NORTH_CURV')
       allocate(x1_tmp(np),x2_tmp(np),x3_tmp(np))
       ! rotate the global coordinates to local curved chunks (by rotation global2local)
       x1_tmp = c1d; x2_tmp = c2d; x3_tmp = c3d
       c1d = this%Mrot_global2local(1,1)*x1_tmp+this%Mrot_global2local(1,2)*x2_tmp+this%Mrot_global2local(1,3)*x3_tmp
       c2d = this%Mrot_global2local(2,1)*x1_tmp+this%Mrot_global2local(2,2)*x2_tmp+this%Mrot_global2local(2,3)*x3_tmp
       c3d = this%Mrot_global2local(3,1)*x1_tmp+this%Mrot_global2local(3,2)*x2_tmp+this%Mrot_global2local(3,3)*x3_tmp
!
       ! for LOCAL_NORTH_CURV, additionally apply the azimuthal rotation to the x- and y- coordinates
       if(this%vtk_projection == 'LOCAL_NORTH_CURV') then
          x1_tmp = c1d; x2_tmp = c2d
          c1d = this%Mrot_azimuth(1,1)*x1_tmp+this%Mrot_azimuth(1,2)*x2_tmp
          c2d = this%Mrot_azimuth(2,1)*x1_tmp+this%Mrot_azimuth(2,2)*x2_tmp
       end if
!
       if(this%apply_vtk_coords_scaling_factor) then
          c1d = c1d * this%vtk_coords_scaling_factor
          c2d = c2d * this%vtk_coords_scaling_factor
          c3d = c3d * this%vtk_coords_scaling_factor
       end if
!
       deallocate(x1_tmp,x2_tmp,x3_tmp)
    case('LOCAL_FLAT','LOCAL_NORTH_FLAT')
       allocate(idx_chunk(np))
       call locateWpLaterallyInsideChunksInversionGrid(this,c1d,c2d,c3d,idx_chunk,np)
       if(any(idx_chunk == -1)) then
          write(errstr,*) count(idx_chunk == -1)," out of ",np," incoming points cannot be located laterally ",&
               "inside any of the ",this%nchunk," chunks; inversion grids of type 'chunksInversionGrid' cannot ",&
               "transform such points"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       call transformVectorGlobalToFlatChunksInversionGrid(this,c1d,c2d,c3d,idx_chunk,np)
       call transformFlatToVtkChunksInversionGrid(this,c1d,c2d,c3d,idx_chunk,np)
       deallocate(idx_chunk)
    end select
!
    ! transform back to single precision output
    c1 = real(c1d)
    c2 = real(c2d)
    c3 = real(c3d)
!
    ! deallocate internal double precision variables
1   deallocate(c1d,c2d,c3d)
  end subroutine transformToVtkChunksInversionGrid
!------------------------------------------------------------------------
!> \brief return geometry information on cells for vtk output
!! \details Dependent on the geometry type defined by the parameter file of this chunks inversion grid (i.e. CELLS or CELL_CENTERS),
!!  define an array of point coordinates and (if geometry type is CELLS) an array defining cells of a
!!  certain cell_type as in vtk file format for unstructured grid.
!!  We need the additional index mapping indx_map here (confer invgrid_vtk_file%req_indx) for the case of 
!!  present(cell_indx_req), as this routine (or rather privateGetGeometryVtkChunksInversionGrid) may throw away invalid and duplicate 
!!  invgrid cell indices. in order to be able to still use the cell geometry information with data vectors 
!!  of same length (and order) as cell_indx_req, indx_map maps the vtk cell index of the returned vtk cells 
!!  to the original position in array cell_indx_req. also, if routine privateGetGeometryVtkChunksInversionGrid does not 
!!  preserve the order of vtk cells as requested by cell_indx_req, map indx_map will reconstruct this order 
!!  (and be the identity otherwise)
!! \param this inversion grid
!! \param geometry_type integer returning the type of geometry: 0 = volumetric cells, 1 = cell center points
!! \param points array of nodes which cell_connectivity array refers to
!! \param cell_connectivity array contianing indices of points (indices having zero offset) as required by vtk format
!! \param cell_types vtk cell types of cells as of vtk convention
!! \param cell_indx_out array of same size as the number of vtk cells (same order), which contains the invgrid cell index of a corresponding vtk cell
!! \param cell_indx_req optional incoming array defining invgrid cell indices for which vtk cells should be returned only
!! \param indx_map if cell_indx_req present, then indx_map(cell_indx_req(i)) = i for all valid and non duplicate cell_indx_req(i) , otherwise identity
!
  subroutine getGeometryVtkChunksInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,& 
       cell_indx_req,indx_map_out)
    type (chunks_inversion_grid) :: this
    integer :: geometry_type
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
    ! local
    character(len=33) :: myname = 'getGeometryVtkChunksInversionGrid'
    character(len=400) :: errstr
    integer, dimension(:), allocatable :: icell_internal_req
    integer :: j
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    if(present(cell_indx_req)) then
       ! incoming array cell_indx_req contains external cell indices, transforming them now to internal ones
       allocate(icell_internal_req(size(cell_indx_req)))
       ! make sure to use only valid internal cell indices as input for privateGetGeometryVtkChunksInversionGrid
       where(cell_indx_req>=1 .and. cell_indx_req<=this%ncell)
          icell_internal_req = this%icell_internal(cell_indx_req)
       elsewhere
          icell_internal_req = -1
       end where
       if(all(icell_internal_req == -1)) then
          write(errstr,*) "there are no valid cell indices among the requested ones; cell indices must "//&
               "be between 1 and ",this%ncell
          call add(errmsg,2,errstr,myname)
          deallocate(icell_internal_req)
          return
       end if
    else
       ! in this case, all external cells are requested, by default
       allocate(icell_internal_req(this%ncell))
       icell_internal_req = (/ (this%icell_internal(j), j=1,this%ncell) /)
    end if
!
    call privateGetGeometryVtkChunksInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,& 
       icell_internal_req,indx_map_out)
!
    ! array cell_indx_out contains internal cell indices, transforming them now to external ones
    cell_indx_out = this%icell_external(cell_indx_out)
!
    deallocate(icell_internal_req)
  end subroutine getGeometryVtkChunksInversionGrid
!------------------------------------------------------------------------
!> \brief do the same as routine getGeometryVtkChunksInversionGrid but for base cells only
!! \details Dependent on the geometry type defined by the parameter file of this chunks inversion grid (i.e. CELLS or CELL_CENTERS),
!!  define an array of point coordinates and (if geometry type is CELLS) an array defining cells of a
!!  certain cell_type as in vtk file format for unstructured grid.
!!  We need the additional index mapping indx_map here (confer invgrid_vtk_file%req_indx) for the case of 
!!  present(cell_indx_req), as this routine (or rather privateGetGeometryVtkChunksInversionGrid) may throw away invalid and duplicate 
!!  invgrid cell indices. in order to be able to still use the cell geometry information with data vectors 
!!  of same length (and order) as cell_indx_req, indx_map maps the vtk cell index of the returned vtk cells 
!!  to the original position in array cell_indx_req. also, if routine privateGetGeometryVtkChunksInversionGrid does not 
!!  preserve the order of vtk cells as requested by cell_indx_req, map indx_map will reconstruct this order 
!!  (and be the identity otherwise)
!! \param this inversion grid
!! \param geometry_type integer returning the type of geometry: 0 = volumetric cells, 1 = cell center points
!! \param points array of nodes which cell_connectivity array refers to
!! \param cell_connectivity array contianing indices of points (indices having zero offset) as required by vtk format
!! \param cell_types vtk cell types of cells as of vtk convention
!! \param cell_indx_out array of same size as the number of vtk cells (same order), which contains the invgrid cell index of a corresponding vtk cell
!! \param cell_indx_req optional incoming array defining invgrid cell indices for which vtk cells should be returned only
!! \param indx_map if cell_indx_req present, then indx_map(cell_indx_req(i)) = i for all valid and non duplicate cell_indx_req(i) , otherwise identity
!
  subroutine getBaseCellGeometryVtkChunksInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,& 
       cell_indx_req,indx_map_out)
    type (chunks_inversion_grid) :: this
    integer :: geometry_type
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
    ! local
    character(len=41) :: myname = 'getBaseCellGeometryVtkChunksInversionGrid'
    character(len=400) :: errstr
    integer, dimension(:), allocatable :: icell_internal_req
    integer :: j
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    if(present(cell_indx_req)) then
       ! incoming array cell_indx_req contains external cell indices, transforming them now to internal ones
       allocate(icell_internal_req(size(cell_indx_req)))
       ! make sure to use only valid internal cell indices as input for privateGetGeometryVtkChunksInversionGrid
       ! note that the base cells cover the first this%ncell_base internal cell indices, i.e. only mark invalid ones
       where(cell_indx_req>=1 .and. cell_indx_req<=this%ncell_base)
          icell_internal_req = cell_indx_req
       elsewhere
          icell_internal_req = -1
       end where
       if(all(icell_internal_req == -1)) then
          write(errstr,*) "there are no valid cell indices among the requested ones; cell indices must "//&
               "be between 1 and the number of base cells = ",this%ncell_base
          call add(errmsg,2,errstr,myname)
          deallocate(icell_internal_req)
          return
       end if
    else
       ! in this case, all external cells are requested, by default
       allocate(icell_internal_req(this%ncell_base))
       icell_internal_req = (/ (j, j=1,this%ncell_base) /)
    end if
!
    call privateGetGeometryVtkChunksInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,& 
       icell_internal_req,indx_map_out)
!
    ! array cell_indx_out contains internal cell indices, which is OK here, since we are returning base cells:
    ! if the inversion grid is refined, the internal index of base cells is ordered as 1,..,this%ncell_base, so OK for external use
    ! if the inversion grid is not refined, then internal and external indices are equal anyway
!
    deallocate(icell_internal_req)
  end subroutine getBaseCellGeometryVtkChunksInversionGrid
!------------------------------------------------------------------------
!> \brief private routine: return geometry information on cells for vtk output
!! \param this inversion grid
!! \param geometry_type integer returning the type of geometry: 0 = volumetric cells, 1 = cell center points
!! \param points array of nodes which cell_connectivity array refers to
!! \param cell_connectivity array contianing indices of points (indices having zero offset) as required by vtk format
!! \param cell_types vtk cell types of cells as of vtk convention
!! \param cell_indx_out array of INTERNAL CELL INDEX, same size as the number of vtk cells (same order), which contains the invgrid cell index of a corresponding vtk cell
!! \param cell_indx_req incoming array of INTERNAL CELL INDICES defining those invgrid cells for which vtk output should be returned only
!! \param indx_map indx_map(cell_indx_req(i)) = i for all valid and non duplicate cell_indx_req(i)
!
  subroutine privateGetGeometryVtkChunksInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,& 
       cell_indx_req,indx_map_out)
    type (chunks_inversion_grid) :: this
    integer :: geometry_type
    integer, dimension(:) :: cell_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
    ! local
    character(len=40) :: myname = 'privateGetGeometryVtkChunksInversionGrid'
    character(len=400) :: errstr
    integer,  parameter :: nrealloc = 400000 ! nrealloc = 400,000 corresponds to arrays points_dble(3,nrealloc) and ichunk_corner(nrealloc) consuming approx. 10 MB memory (this reallocation step should be OK for normal hardware)
    double precision, dimension(:,:), pointer :: points_dble
    double precision :: ix_dble,iy_dble,iz_dble
    integer, dimension(:), pointer :: indx_map
    logical, dimension(:), allocatable :: cell_occurred
    logical :: select_cell_indices
    integer :: ncell_return,ipoint,npoint_return,n,i,ishift
    integer :: icell,icell_internal,ix,iy,iz,ichunk,point_index_of_corner
    integer, dimension(:,:,:,:), allocatable :: idx_point
    integer, dimension(:), pointer :: ichunk_corner
    logical :: use_idx_point
!
    nullify(indx_map,points_dble,ichunk_corner)
!
    call addTrace(errmsg,myname)
    nullify(points,cell_connectivity,cell_type,cell_indx_out)
    if(present(indx_map_out)) nullify(indx_map_out)
!
! DEFINE THOSE INVGRID CELLS FOR WHICH THE GEOMETRY IS RETURNED
!
    ! only if there are any indices in valid range, select those cells and remove duplicate ones.
    ! otherwise return no cells (nullified pointers)
    if(any(cell_indx_req .ge. 1 .and. cell_indx_req .le. this%ncell_internal)) then
       select_cell_indices = .true.
       ncell_return = 0
       ! define mapping indx_map_tmp which maps invgrid index to index in array cell_indx_req
       allocate(indx_map(this%ncell_internal),cell_indx_out(this%ncell_internal)) ! allocate for maximum number of cells returned (reallocate later, if necessary)
       allocate(cell_occurred(this%ncell_internal)); cell_occurred(:) = .false.
       do i = 1,size(cell_indx_req)
          if(cell_indx_req(i)>=1 .and. cell_indx_req(i)<=this%ncell_internal) then
             ! only if this valid index did not occurr so far, take it.
             ! this way, always the FIRST occuring index (in case of duplicate indices) is returned, 
             ! the others are lost. But in this case, you should not use the indx_map anyway
             if(.not.cell_occurred(cell_indx_req(i))) then
                ncell_return = ncell_return + 1
                indx_map(ncell_return) = i
                cell_indx_out(ncell_return) = cell_indx_req(i)
                cell_occurred(cell_indx_req(i)) = .true.
             end if
          endif
       enddo ! i
       deallocate(cell_occurred)
       if(ncell_return < this%ncell_internal) then
          indx_map => reallocate(indx_map,ncell_return)
          cell_indx_out => reallocate(cell_indx_out,ncell_return)
       end if
    else ! any(cell_indx_req .ge. 1 .and. cell_indx_req .le. this%ncell_internal)
       write(errstr,*) "there are no valid internal indices among the requested ones; internal cell indices must "//&
            "be between 1 and ",this%ncell_internal
       call add(errmsg,2,errstr,myname)
       return
    endif ! any(cell_indx_req .ge. 1 .and. cell_indx_req .le. this%ncell_internal)
!
! LOOP ON ALL CELLS TO BE RETURNED AND DEFINE POINTS AND CELL CONNECTIVITY
!
    select case(this%vtk_geometry_type_int)
    case(0) ! CELLS
       ! TRY TO ALLOCATE ARRAY idx_point (can be very very large for large refined inversion grids)
       allocate(idx_point(this%nx,this%ny,this%nz,this%nchunk),stat=i)
       use_idx_point = (i == 0)
       if(use_idx_point) then
          idx_point = -1
       else
          write(*,*) "WARNING in getGeometryVtkChunksInversionGrid: too many x,y,z coordinates, cannot ",&
               "allocate array for memorizing already collected point indices (would be ",&
               4*this%nx*this%ny*this%nz*this%nchunk," bytes)-> APPLYING A DIFFERENT METHOD FOR ",&
               "THAT PURPOSE WHICH IS COMPUTATIONALLY EXPENSIVE (so the program might slow down now significantly!)"
       end if
!
       if(this%contains_refined_cells) then
          ! the number this%nx*this%ny*this%nz*this%nchunk is much too large and not sensible to be used
          ! for pre-allocation of points_dble. In case of cell refinement, it makes sense to reallocate
          ! the arrays
          allocate(points_dble(3,nrealloc),ichunk_corner(nrealloc))
       else
          ! allocate for maximum number of vtk points (the following is the exact number in case of NO cell 
          ! refinement); in the end, reallocate array points for the actual number of found corners
          allocate(points_dble(3,this%nx*this%ny*this%nz*this%nchunk),ichunk_corner(this%nx*this%ny*this%nz*this%nchunk))
       end if
       ! need to track the chunk index for each corner, since only (a subset of) external cells are processed here
       ichunk_corner(:) = -1
!
       allocate(cell_connectivity((8+1)*ncell_return),cell_type(ncell_return))
!       ! so far, all cells are hexahedra, hence cell_type equals 12 everywhere (by vtk standard)
       cell_type(:) = 12
!
       ishift = 0 ! remember the current position in array cell_connectivity
       npoint_return = 0
       do icell = 1,ncell_return
          ! first entry in cell_connectivity for this cell: number of indices to come (always = 8)
          cell_connectivity(1+ishift) = 8
!
          icell_internal = cell_indx_out(icell)
          ichunk = this%cell_ichunk(icell_internal)
!
          ! 1st corner xmin,ymin,zmin
          ix = this%cell_ixmin(icell_internal)
          iy = this%cell_iymin(icell_internal)
          iz = this%cell_izmin(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(2+ishift) = point_index_of_corner
!
          ! 2nd corner xmax,ymin,zmin
          ix = this%cell_ixmax(icell_internal)
          iy = this%cell_iymin(icell_internal)
          iz = this%cell_izmin(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(3+ishift) = point_index_of_corner
!
          ! 3rd corner xmax,ymax,zmin
          ix = this%cell_ixmax(icell_internal)
          iy = this%cell_iymax(icell_internal)
          iz = this%cell_izmin(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(4+ishift) = point_index_of_corner
!
          ! 4th corner xmin,ymax,zmin
          ix = this%cell_ixmin(icell_internal)
          iy = this%cell_iymax(icell_internal)
          iz = this%cell_izmin(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(5+ishift) = point_index_of_corner
!
          ! 5th corner xmin,ymin,zmax
          ix = this%cell_ixmin(icell_internal)
          iy = this%cell_iymin(icell_internal)
          iz = this%cell_izmax(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(6+ishift) = point_index_of_corner
!
          ! 6th corner xmax,ymin,zmax
          ix = this%cell_ixmax(icell_internal)
          iy = this%cell_iymin(icell_internal)
          iz = this%cell_izmax(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(7+ishift) = point_index_of_corner
!
          ! 7th corner xmax,ymax,zmax
          ix = this%cell_ixmax(icell_internal)
          iy = this%cell_iymax(icell_internal)
          iz = this%cell_izmax(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(8+ishift) = point_index_of_corner
!
          ! 8th corner xmin,ymax,zmax
          ix = this%cell_ixmin(icell_internal)
          iy = this%cell_iymax(icell_internal)
          iz = this%cell_izmax(icell_internal)
          call addCornerIfNotAlreadyAdded() ! in this subroutine, the value of index point_index_of_corner is set
          cell_connectivity(9+ishift) = point_index_of_corner
!
          ishift = ishift + 9
       end do ! icell
!
       if(allocated(idx_point)) deallocate(idx_point)
!
       if(.not.use_idx_point) then
          ! convert integer indices still contained in array points_dble to actual point coordinates
          do ipoint = 1,npoint_return
             ix = int(real(points_dble(1,ipoint)))
             points_dble(1,ipoint) = this%x(ix)
             iy = int(real(points_dble(2,ipoint)))
             points_dble(2,ipoint) = this%y(iy)
             iz = int(real(points_dble(3,ipoint)))
             points_dble(3,ipoint) = this%z(iz)
          end do ! ipoint
       end if
!
       geometry_type = 0
!
    case(1) ! CELL CENTERS
       ! do not define arrays cell_connectivity and cell_type in this case. Only return points as cell centers
       allocate(points_dble(3,ncell_return),ichunk_corner(ncell_return))
       do icell = 1,ncell_return
          icell_internal = this%icell_internal(cell_indx_out(icell))
          points_dble(1,icell) = this%cell_center(1,icell_internal)
          points_dble(2,icell) = this%cell_center(2,icell_internal)
          points_dble(3,icell) = this%cell_center(3,icell_internal)
          ichunk_corner(icell) = this%cell_ichunk(icell_internal)
       end do ! icell
       npoint_return = ncell_return
       geometry_type = 1
    end select ! this%vtk_geometry_type
!
! TRANSFORM POINTS ARRAY AND CONVERT TO SINGLE PRECISION
!
    ! transform the flat point coordinates to VTK projection
    call transformFlatToVtkChunksInversionGrid(this,points_dble(1,1:npoint_return),points_dble(2,1:npoint_return),&
         points_dble(3,1:npoint_return),ichunk_corner(1:npoint_return),npoint_return)
!
    allocate(points(3,npoint_return))
    points = real(points_dble(:,1:npoint_return))
!
    if(associated(ichunk_corner)) deallocate(ichunk_corner)
    if(associated(points_dble)) deallocate(points_dble)
!
    if(present(indx_map_out)) then
       indx_map_out => indx_map
       nullify(indx_map)
    else
       deallocate(indx_map)
    endif
!
  contains
!
    subroutine addCornerIfNotAlreadyAdded()
      ! NO LOCAL DECLARATION OF VARIABLES, USE THE (global) ONES DECLARED IN subroutine getGeometryVtkChunksInversionGrid
!
      ! check if this corner is added to points_dble (dependent on whether idx_point is used or not); 
      ! if yes, define the existing index point_index_of_corner and return
!
      if(use_idx_point) then
         if(idx_point(ix,iy,iz,ichunk) > 0) then
            point_index_of_corner = idx_point(ix,iy,iz,ichunk) -1
            return
         end if
      else ! use_idx_point
         ix_dble = dble(ix)
         iy_dble = dble(iy)
         iz_dble = dble(iz)
         do ipoint = 1,npoint_return
            if(ichunk_corner(ipoint) == ichunk) then
               if(points_dble(1,ipoint)==ix_dble) then
                  if(points_dble(2,ipoint)==iy_dble) then
                     if(points_dble(3,ipoint)==iz_dble) then
                        point_index_of_corner = ipoint-1 ! indices in vtk files have offset 0, so always store ipoint-1
                        return
                     end if
                  end if
               end if
            end if
         end do ! ipoint
      end if ! use_idx_point
!
      ! if subroutine comes to this line of code, the current cell corner is NOT yet contained in points_dble, so add it 
!
      ! for either case (use_idx_point or not) increase global counter of points in points_dble and reallocate if necessary
      npoint_return = npoint_return + 1
      n = size(ichunk_corner)
      if(npoint_return > n) then
         ! reallocate arrays points_dble, ichunk_corner
         points_dble => reallocate(points_dble,3,n+nrealloc)
         ichunk_corner => reallocate(ichunk_corner,n+nrealloc)
      end if
!
      if(use_idx_point) then
         ! store the actual coordinates in array points_dble
         points_dble(1,npoint_return) = this%x(ix)
         points_dble(2,npoint_return) = this%y(iy)
         points_dble(3,npoint_return) = this%z(iz)
         ! memorize the index of this corner in array points_dble
         idx_point(ix,iy,iz,ichunk) = npoint_return
      else ! use_idx_point
         ! if the loop above did not find the corner as an existing point in array points_dble, 
         ! add the integer indices (in arrays this%x,this%y,this%z) of this point as dbles 
         ! to array points_dble
         points_dble(1,npoint_return) = ix_dble
         points_dble(2,npoint_return) = iy_dble
         points_dble(3,npoint_return) = iz_dble
      end if ! use_idx_point
!
      ! for either case (use_idx_point or not), store the chunk index into ichunk_corner and define point_index_of_corner as the new index just added
      ichunk_corner(npoint_return) = ichunk
      point_index_of_corner = npoint_return - 1  ! indices in vtk files have offset 0, so always store ipoint-1
    end subroutine addCornerIfNotAlreadyAdded
  end subroutine privateGetGeometryVtkChunksInversionGrid
!------------------------------------------------------------------------
!> \brief return indices of all face neighbours for all (optionally only subset of) cells
!!  \details If boundary_conditions is set, there should be a special form of nb_idx returned, which is used
!!  e.g. to define special smoothing condigions:
!!   'no_nb_inner_bnd': with some neighbours removed (which should not be smoothed with on internal invgrid 
!!                      boundaries). FOR chunksInversionGrid GRIDS, THERE IS NO POSSIBILITY YET TO DEFINE INTERNAL
!!                      INVGRID BOUNDARIES!! SO IF THIS IS VALUE IS SET, NO NEIGHBOURS ARE RETURNED! 
!!                      In the future might introduce actual functionality here. A simple way to introduce this is
!!                      introducing certain layers of base cells below which there should be an internal boundary.
!!                      (of even easier, to define a set of refinement blocks of base cells, below which there should be 
!!                      an internal boundary). 
!!   'extra_nbs_outer_bnd': additional fake neighbours (having cell index 0) in case of zero boundary conditions 
!!                          on outer invgrid boundaries.
!!   'extra_nbs_outer_bnd_except_free_surface': same as 'extra_nbs_outer_bnd', but not applied for cells on free surfaces.
!!   '','standard': no special boundary handling, i.e. standard functionality (exactly the geometrical face neighbours are given)
!! \param this chunks inversion grid
!! \param nb_idx pointer to array of length this%ncell which contains vector pointer to face neighbour indices
!! \param boundary_conditions optional string indicating type of boundary conditions for which neighbours should be returned
!
  subroutine getIndicesFaceNeighboursChunksInversionGrid(this,nb_idx,boundary_conditions)
    type (chunks_inversion_grid) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
    character(len=*), optional :: boundary_conditions
    ! local
    integer :: icell,nnb_xmin,nnb_xmax,nnb_ymin,nnb_ymax,nnb_zmin,nnb_zmax,nnb,inb
    integer, dimension(:), pointer :: nb_xmin,nb_xmax,nb_ymin,nb_ymax,nb_zmin,nb_zmax,nb
    character(len=39) :: bnd_cond
!
    nullify(nb_xmin,nb_xmax,nb_ymin,nb_ymax,nb_zmin,nb_zmax,nb)
!
    nullify(nb_idx)
!
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
!
    do icell = 1,this%ncell
       ! we need to be careful here to account only for external cells as neighbours. Note that vectors 
       ! this%cell_nb_xmin contain neighbour indices of internal cells!

       ! find the number of EXTERNAL CELLS which are neighbours in a particular direction of space
       nnb_xmin = 0; nnb_xmax = 0; nnb_ymin = 0; nnb_ymax = 0; nnb_zmin = 0; nnb_zmax = 0

       nb_xmin => getVectorPointer(this%cell_nb_xmin(this%icell_internal(icell)))
       if(associated(nb_xmin)) nnb_xmin = count(this%icell_external(nb_xmin(:))/=-1)

       nb_xmax => getVectorPointer(this%cell_nb_xmax(this%icell_internal(icell)))
       if(associated(nb_xmax)) nnb_xmax = count(this%icell_external(nb_xmax(:))/=-1)

       nb_ymin => getVectorPointer(this%cell_nb_ymin(this%icell_internal(icell)))
       if(associated(nb_ymin)) nnb_ymin = count(this%icell_external(nb_ymin(:))/=-1)

       nb_ymax => getVectorPointer(this%cell_nb_ymax(this%icell_internal(icell)))
       if(associated(nb_ymax)) nnb_ymax = count(this%icell_external(nb_ymax(:))/=-1)

       nb_zmin => getVectorPointer(this%cell_nb_zmin(this%icell_internal(icell)))
       if(associated(nb_zmin)) nnb_zmin = count(this%icell_external(nb_zmin(:))/=-1)

       nb_zmax => getVectorPointer(this%cell_nb_zmax(this%icell_internal(icell)))
       if(associated(nb_zmax)) nnb_zmax = count(this%icell_external(nb_zmax(:))/=-1)
!
       ! Dependent on the type of boundary condition, define the neighbour index vector nb
       select case(bnd_cond)
!
       case ('extra_nbs_outer_bnd','extra_nbs_outer_bnd_except_free_surface','standard')
!
          ! First define the total number of neighbours (including artificial) for current cell
          select case(bnd_cond)
          case ('standard')
             ! give excatly the geometrical face neighbours of this cell
             nnb = nnb_xmin + nnb_xmax + nnb_ymin + nnb_ymax + nnb_zmin + nnb_zmax
          case ('extra_nbs_outer_bnd','extra_nbs_outer_bnd_except_free_surface')
             nnb = 0
             ! For each direction, test if there are neighbours. If yes, add those, otherwise add one artificial.
             if(nnb_xmin > 0) then
                nnb = nnb + nnb_xmin
             else
                nnb = nnb + 1
             end if
             if(nnb_xmax > 0) then
                nnb = nnb + nnb_xmax
             else
                nnb = nnb + 1
             end if
             if(nnb_ymin > 0) then
                nnb = nnb + nnb_ymin
             else
                nnb = nnb + 1
             end if
             if(nnb_ymax > 0) then
                nnb = nnb + nnb_ymax
             else
                nnb = nnb + 1
             end if
             if(nnb_zmin > 0) then
                nnb = nnb + nnb_zmin
             else
                nnb = nnb + 1
             end if
             if(nnb_zmax > 0) then
                nnb = nnb + nnb_zmax
             else
                ! FOR NEIGHBOURS ABOVE THIS CELL, ONLY ADD ARTIFICIAL NEIGHBOUR FOR 'extra_nbs_outer_bnd'.
                ! IN CASE OF 'extra_nbs_outer_bnd_except_free_surface', ONLY GEOMETRICAL NEIGHBOURS HERE!
                if(bnd_cond=='extra_nbs_outer_bnd') nnb = nnb + 1
             end if
          end select ! bnd_cond
!
          ! Secondly, build the vector of neighbour indices for this cell, length nnb
          if(nnb > 0) then
             allocate(nb(nnb))
             ! initialize vector of neighbour indices by -1 , indicating artificial neighbours
             nb(:) = -1
             ! Then overwrite by all actual geometrical neighbours. If some tailing indices stay -1, this 
             ! is intended to indicate artificial neighbours (controlled by increasing the number nnb 
             ! in select statement above).
             inb = 0 ! counter on all entries in vector nb
             if(nnb_xmin > 0) then
                nb(inb+1:inb+nnb_xmin) = pack( this%icell_external(nb_xmin(:)) , this%icell_external(nb_xmin(:))/=-1 )
                inb = inb + nnb_xmin
             end if
             if(nnb_xmax > 0) then
                nb(inb+1:inb+nnb_xmax) = pack( this%icell_external(nb_xmax(:)) , this%icell_external(nb_xmax(:))/=-1 )
                inb = inb + nnb_xmax
             end if
             if(nnb_ymin > 0) then
                nb(inb+1:inb+nnb_ymin) = pack( this%icell_external(nb_ymin(:)) , this%icell_external(nb_ymin(:))/=-1 )
                inb = inb + nnb_ymin
             end if
             if(nnb_ymax > 0) then
                nb(inb+1:inb+nnb_ymax) = pack( this%icell_external(nb_ymax(:)) , this%icell_external(nb_ymax(:))/=-1 )
                inb = inb + nnb_ymax
             end if
             if(nnb_zmin > 0) then
                nb(inb+1:inb+nnb_zmin) = pack( this%icell_external(nb_zmin(:)) , this%icell_external(nb_zmin(:))/=-1 )
                inb = inb + nnb_zmin
             end if
             if(nnb_zmax > 0) then
                nb(inb+1:inb+nnb_zmax) = pack( this%icell_external(nb_zmax(:)) , this%icell_external(nb_zmax(:))/=-1 )
                inb = inb + nnb_zmax
             end if
          end if ! nnb > 0
!
       ! case ('no_nb_inner_bnd')
       !    either write completely new case statement here, or integrate 'no_nb_inner_bnd' 
       !    in the select statements above
       end select ! case bnd_cond
!
       ! Finally associate the neighbour index vector pointer for this cell.
       ! Even if nb was not allocated above (and at this point is still not associated), it works here (nothing happens).
       call associateVectorPointer(nb_idx(icell),nb); nullify(nb)
    end do ! icell
  end subroutine getIndicesFaceNeighboursChunksInversionGrid
!------------------------------------------------------------------------
!> \brief for each cell return indices of wavefield points contained in that cell
!! \param this chunks inversion grid
!! \param x vector or first coordinate of wavefield points
!! \param y vector or second coordinate of wavefield points
!! \param z vector or third coordinate of wavefield points
!! \param uf_wp unit factor of wavefield points
!! \param wp_idx pointer to array of length this%ncell; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideChunksInversionGrid(this,x,y,z,uf_wp,wp_idx,errmsg)
    type (chunks_inversion_grid) :: this
    real, dimension(:), intent(in) :: x,y,z
    real :: uf_wp
    type (integer_vector_pointer), dimension(:), pointer :: wp_idx
    type (error_message) :: errmsg
    ! local
    character(len=33) :: myname = 'locateWpInsideChunksInversionGrid'
    character(len=400) :: errstr
    double precision, dimension(:), allocatable :: r,xd,yd,zd,xd_flat,yd_flat,zd_flat
    logical, dimension(:), allocatable :: mask_r,mask_ichunk
    integer, dimension(:), allocatable :: ichunk_tmp,wp_ichunk
    integer, dimension(:), pointer :: idx
    integer, dimension(:), allocatable :: idx_r,idx_flat,nwp_chunk,iwp_chunk,wp_nidx
    integer :: iwp,iwp_flat,nwp,nwp_flat,n,nr
    type(integer_vector_pointer), dimension(:,:), allocatable :: test_vectors_wp
    integer :: ichunk,nwp_to_be_located,icell,icell_internal,icell_external
    integer, dimension(:), pointer :: jzmax_wp,jxmax_wp,jymax_wp,iglob_wp,chunk_icell_internal
    integer :: jxmin,jxmax,jymin,jymax,jzmin,jzmax
    double precision :: uf_correction
!
    nullify(idx,jzmax_wp,jxmax_wp,jymax_wp,iglob_wp,chunk_icell_internal)
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
    ! that x,y,z are not empty and are all of same length!
    nwp = size(x)
!
    ! convert to double precision for internal use
    allocate(xd(nwp),yd(nwp),zd(nwp))
    if(uf_wp /= 1.0e3) then
       ! account for wavefield points which are not given in km (e.g. using SPECFEM3D_Cartesian (or NEXD) with SI
       ! units for spherical applications)

       ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
       uf_correction = dble(uf_wp) * 1.0d-3
       xd = dble(x) * uf_correction
       yd = dble(y) * uf_correction
       zd = dble(z) * uf_correction
    else
       xd = dble(x)
       yd = dble(y)
       zd = dble(z)
    end if
!
    ! First check which points are in radial range of this inversion grid
    allocate(r(nwp),mask_r(nwp))
    r = dsqrt( xd**2 + yd**2 + zd**2)
    mask_r = r <= dble(this%rmax) .and. r >= this%z(this%sorted2iz(1))
    deallocate(r)
    nr = count(mask_r)
    if(nr == 0) goto 1
    allocate(idx_r(nr)); idx_r = pack( (/ (iwp,iwp=1,nwp) /) , mask_r )
    deallocate(mask_r)
!
    ! Locate those points laterally inside the chunks, which are in depth range of this inversion grid
    ! allocate ichunk_tmp with size nr, since for intent(out) dummy variables, array-section actual arguments 
    ! with vector subscripts are not allowed (could NOT use wp_ichunkg(idx_r) here!!).
    ! In order to select correct indices, allocate wp_ichunk having the total size nwp (just like xd,yd,zd).
    allocate(ichunk_tmp(nr),wp_ichunk(nwp))
    call locateWpLaterallyInsideChunksInversionGrid(this,xd(idx_r),yd(idx_r),zd(idx_r),ichunk_tmp,nr)
    wp_ichunk = -1; wp_ichunk(idx_r) = ichunk_tmp
    deallocate(ichunk_tmp)
!
    ! Transform those points which could be laterally located in some chunk, to that flat reference chunk.
    ! This allows to locate those points inside the inversion grid cells.
    ! Those points are defined as the ones for which wp_ichunk /= -1 .
    allocate(mask_ichunk(nr))
    mask_ichunk = wp_ichunk(idx_r) /= -1
    nwp_flat = count(mask_ichunk)
    if(nwp_flat == 0) goto 1
    allocate(idx_flat(nwp_flat)); idx_flat = pack( idx_r , mask_ichunk )
    deallocate(idx_r,mask_ichunk)
    allocate(xd_flat(nwp_flat),yd_flat(nwp_flat),zd_flat(nwp_flat))
    xd_flat = xd(idx_flat)
    yd_flat = yd(idx_flat)
    zd_flat = zd(idx_flat)
    deallocate(xd,yd,zd)
    call transformVectorGlobalToFlatChunksInversionGrid(this,xd_flat,yd_flat,zd_flat,wp_ichunk(idx_flat),nwp_flat)
!
    ! now sort the wavefield points by chunk and pre-compute indices of sorted cell boundaries beneath which the wavefield points are located
!
    allocate(nwp_chunk(this%nchunk),iwp_chunk(this%nchunk),test_vectors_wp(4,this%nchunk))
    do ichunk = 1,this%nchunk
       n = count(wp_ichunk(idx_flat) == ichunk) ! wp_ichunk can have values 1..this%nchunk (locateWpLaterally checks only existing chunks!)
       nwp_chunk(ichunk) = n
       call allocateVectorPointer(test_vectors_wp(1,ichunk),n)
       call allocateVectorPointer(test_vectors_wp(2,ichunk),n)
       call allocateVectorPointer(test_vectors_wp(3,ichunk),n)
       call allocateVectorPointer(test_vectors_wp(4,ichunk),n)
    end do ! ichunk
!
    iwp_chunk(:) = 0
    do iwp_flat = 1,nwp_flat
       ichunk = wp_ichunk(idx_flat(iwp_flat))
!
       ! locate point between all possible coordinates of cell boundaries
       jzmax = locateCoordinateChunksInversionGrid(zd_flat(iwp_flat),this%z(this%sorted2iz),this%nz,is_ascending=.true.)
       jxmax = locateCoordinateChunksInversionGrid(xd_flat(iwp_flat),this%x(this%sorted2ix),this%nx,is_ascending=.true.)
       jymax = locateCoordinateChunksInversionGrid(yd_flat(iwp_flat),this%y(this%sorted2iy),this%ny,is_ascending=.true.)
!
       ! increase counter of wavefield points contained in chunk ichunk
       iwp_chunk(ichunk) = iwp_chunk(ichunk) + 1
       n = iwp_chunk(ichunk)
       ! memorize all values of j[x,y,z]max and global wp index for this wavefield point
       call fillVectorPointer(test_vectors_wp(1,ichunk),(/jzmax/),n) ! jzmax
       call fillVectorPointer(test_vectors_wp(2,ichunk),(/jxmax/),n) ! jxmax
       call fillVectorPointer(test_vectors_wp(3,ichunk),(/jymax/),n) ! jymax
       call fillVectorPointer(test_vectors_wp(4,ichunk),(/idx_flat(iwp_flat)/),n) ! global index of wavefield point
    end do ! iwp_flat
    deallocate(wp_ichunk,idx_flat,iwp_chunk)
!
    ! now loop over all chunks
    do ichunk = 1,this%nchunk

       jzmax_wp => getVectorPointer(test_vectors_wp(1,ichunk))
       jxmax_wp => getVectorPointer(test_vectors_wp(2,ichunk))
       jymax_wp => getVectorPointer(test_vectors_wp(3,ichunk))
       iglob_wp => getVectorPointer(test_vectors_wp(4,ichunk))

       ! loop over all external cells of this chunk, for each cell find all wavefield points contained in it 
       ! (test only wavefield points contained in this chunk)

       ! Use a counter of the number of wavefield points which were not yet located.
       ! If this counter becomes zero, exit loop on cells of this chunk.
       nwp_to_be_located = nwp_chunk(ichunk)
!
       chunk_icell_internal => getVectorPointer(this%chunk_icell_internal(ichunk))
       do icell = 1,this%chunk_ncell_internal(ichunk)
!
          ! if in this chunk, there are no (more) wavefield points to be located, exit this loop on the cells
          if(nwp_to_be_located <= 0) exit
!
          icell_internal = chunk_icell_internal(icell)
          icell_external = this%icell_external(icell_internal)
          ! only treat external cells here, do not locate wavefield points inside internal invgrid cells
          if(icell_external < 0) cycle
!
          ! get sorted indices of boundary coordinates of this cell
          jxmin = this%ix2sorted(this%cell_ixmin(this%icell_internal(icell)))
          jxmax = this%ix2sorted(this%cell_ixmax(this%icell_internal(icell)))
          jymin = this%iy2sorted(this%cell_iymin(this%icell_internal(icell)))
          jymax = this%iy2sorted(this%cell_iymax(this%icell_internal(icell)))
          jzmin = this%iz2sorted(this%cell_izmin(this%icell_internal(icell)))
          jzmax = this%iz2sorted(this%cell_izmax(this%icell_internal(icell)))
!
          do iwp = 1,nwp_chunk(ichunk)

             ! if this point is already located, cycle 
             if(iglob_wp(iwp) < 0) cycle ! (memorize that this point is located by setting iglob_wp = -1 after locating it) 

             ! A wavefield point is defined to be contained in the current cell if for each dimension x,y,z
             ! the point's coordinate is <= the upper cell boundary and > the lower cell boundary.
             ! Points which lie on the lower coordinate boundary of the inversion grid are located in the
             ! boundary cell (for which its lower boundary coincides with the lower boundary of the inversion grid)

             ! The following if-statements are, however, formulated in a negated way, i.e. they test whether
             ! a point does NOT fulfill the above conditions and, hence, is rejected. 
             ! For example testing coordinate z:
             !
             !   if((jzmax_wp(iwp) <= jzmin .and. jzmin > 1) .or. jzmax_wp(iwp) > jzmax) cycle
             !
             ! Note that jzmax_wp(iwp) is the smallest cell z-coordinate which is larger or equal to the
             ! z-coordinate of the wavefield point. Therefore a point is rejected if jzmax_wp(iwp) > jzmax
             ! or if jzmax_wp(iwp) <= jzmin in case that we are not testing a cell at the z-min boundary of the
             ! inversion grid (because in that case a point should lie in the cell even if jzmax_wp(iwp) == jzmin == 1).
             ! 
             ! If a point "survives" all of the following x,y,z if-tests (and is not rejected), it is contained in the current cell.

             ! The following if-statements of the form 
             !   if((jzmax_wp(iwp) <= jzmin .and. jzmin > 1) .or. jzmax_wp(iwp) > jzmax) cycle
             ! are subdivided into their constituents. I.e. there are equivalently two successive if-statements of the form
             !    if(jzmax_wp(iwp) <= jzmin) then
             !       if(jzmin > 1) cycle
             !    end if
             !    if(jzmax_wp(iwp) > jzmax) cycle
             ! Florian Schumacher does not know whether this is more computationally efficient compared with testing
             ! the composite statement (probably depends on how the compiler structures such statements)
             ! The idea is, e.g., only to evaluate the second if-statement if the first one was true (instead of
             ! always evaluating both and looking at the logical comparison).

             if(jzmax_wp(iwp) <= jzmin) then
                if(jzmin > 1) cycle
             end if
             if(jzmax_wp(iwp) > jzmax) cycle

             if(jxmax_wp(iwp) <= jxmin) then
                if(jxmin > 1) cycle
             end if
             if(jxmax_wp(iwp) > jxmax) cycle

             if(jymax_wp(iwp) <= jymin) then
                if(jymin > 1) cycle
             end if
             if(jymax_wp(iwp) > jymax) cycle

             ! If code comes here, the current wavefield point is contained in the current cell

             ! if number of wavefield points already located into this cell, wp_nidx(icell_external), equals
             ! the allocated  size of pointer in wp_idx(icell_external), reallocate 
             ! (does also work in the beginning, when wp_nidx(icell_external)==0)
             if(mod(wp_nidx(icell_external),500) == 0) then
                idx => getVectorPointer(wp_idx(icell_external))
                idx => reallocate(idx,wp_nidx(icell_external)+500)
                call associateVectorPointer(wp_idx(icell_external),idx)
             end if
!
             ! add this wavefield point index iwp to list of indices of cell icell_external
             wp_nidx(icell_external) = wp_nidx(icell_external) + 1
             call fillVectorPointer(wp_idx(icell_external),(/iglob_wp(iwp)/),wp_nidx(icell_external))
!
             ! decrease counter of all points to be located inside this chunk
             nwp_to_be_located = nwp_to_be_located - 1
!
             ! memorize that this point was located in some cell by setting iglob_wp = -1 here
             iglob_wp(iwp) = -1
          end do ! iwp
!
          ! reallocate wp_index array for this cell to correct size, if there were any points located in this cell 
          ! and if there are less valid entries than its actual size
          if(wp_nidx(icell_external) /= 0) then
             idx => getVectorPointer(wp_idx(icell_external))
             if(wp_nidx(icell_external) /= size(idx)) then
                idx => reallocate(idx,wp_nidx(icell_external))
                call associateVectorPointer(wp_idx(icell_external),idx)
             end if
          end if

       end do ! icell

    end do ! ichunk
!
1   if(nwp-sum(wp_nidx) > 0) then
       write(errstr,*) nwp-sum(wp_nidx)," wavefield points (out of ",nwp,&
            ") could not be located inside the inversion grid"
       call add(errmsg,1,errstr,myname)
    end if
!
    if(allocated(r)) deallocate(r)
    if(allocated(xd)) deallocate(xd)
    if(allocated(yd)) deallocate(yd)
    if(allocated(zd)) deallocate(zd)
    if(allocated(xd_flat)) deallocate(xd_flat)
    if(allocated(yd_flat)) deallocate(yd_flat)
    if(allocated(zd_flat)) deallocate(zd_flat)
    if(allocated(mask_r)) deallocate(mask_r)
    if(allocated(mask_ichunk)) deallocate(mask_ichunk)
    if(allocated(ichunk_tmp)) deallocate(ichunk_tmp)
    if(allocated(wp_ichunk)) deallocate(wp_ichunk)
    if(allocated(idx_r)) deallocate(idx_r)
    if(allocated(idx_flat)) deallocate(idx_flat)
    if(allocated(nwp_chunk)) deallocate(nwp_chunk)
    if(allocated(iwp_chunk)) deallocate(iwp_chunk)
    if(allocated(wp_nidx)) deallocate(wp_nidx)
    if(allocated(test_vectors_wp)) then
       do n = 1,size(test_vectors_wp,1)
          do nr = 1,size(test_vectors_wp,2)
             call dealloc(test_vectors_wp(n,nr))
          end do ! nr
       end do ! n
       deallocate(test_vectors_wp)
    end if
  end subroutine locateWpInsideChunksInversionGrid
!------------------------------------------------------------------------
  subroutine locateWpLaterallyInsideChunksInversionGrid(this,x1,x2,x3,idx_chunk,n)
    type (chunks_inversion_grid) :: this
    integer :: n
    double precision, dimension(n) :: x1,x2,x3
    integer, dimension(n), intent(out) :: idx_chunk
    ! local
    integer :: i,ichunk
    double precision :: check_lat,check_lon
    double precision :: x1_tmp,x2_tmp,x3_tmp
!
    ! after rotating to local curved coordinates, a point (x1,x2,x3) lies laterally in the chunk if
    !   1) x3 >= 0  AND
    !   2) (abs(x1)/x3) <= tan(wlat/2)  AND
    !   3) (abs(x2)/x3) <= tan(wlon/2)
    ! the implementation below tries to avoid any unnecessary computations and tries to check things
    ! as soon as possible.
!
    if(this%wlat == 90.0) then
       check_lat = 1.d0 ! = tan(0.5 * 90.0 * mc_deg2radd)
    else
       check_lat = dtan(0.5d0*dble(this%wlat)*mc_deg2radd)
    end if
    if(this%wlon == 90.0) then
       check_lon = 1.d0 ! = tan(0.5 * 90.0 * mc_deg2radd)
    else
       check_lon = dtan(0.5d0*dble(this%wlon)*mc_deg2radd)
    end if
!
    do i = 1,n
       idx_chunk(i) = -1
       do ichunk = 1,this%nchunk
          x3_tmp = this%Mrot_global2flat(3,1,ichunk)*x1(i)+this%Mrot_global2flat(3,2,ichunk)*x2(i)+&
               this%Mrot_global2flat(3,3,ichunk)*x3(i)
          if(x3_tmp >= 0.d0) then
             x1_tmp = this%Mrot_global2flat(1,1,ichunk)*x1(i)+this%Mrot_global2flat(1,2,ichunk)*x2(i)+&
                  this%Mrot_global2flat(1,3,ichunk)*x3(i)
             if(abs(x1_tmp) > x3_tmp*check_lat) cycle
             x2_tmp = this%Mrot_global2flat(2,1,ichunk)*x1(i)+this%Mrot_global2flat(2,2,ichunk)*x2(i)+&
                  this%Mrot_global2flat(2,3,ichunk)*x3(i)
             if(abs(x2_tmp) <= x3_tmp*check_lon) then
                idx_chunk(i) = ichunk
                exit ! point is located laterally in this chunk, so exit loop over ichunk and continue with next point, i.e. next index in loop over i
             end if
          end if
       end do ! ichunk
    end do ! i
  end subroutine locateWpLaterallyInsideChunksInversionGrid
!------------------------------------------------------------------------
!> \brief return index of interval + 1 in coords array which contains c
!> \details regardless of array coords being ascending or descending, two neighbouring entries
!!  coords(j),coords(j+1) of array coords define an interval (which defined a cell in this coordinate direction)
!!  which ALWAYS contains the left value coords(j) and does NOT contain the right value coords(j+1), 
!!  EXCEPT for the right-most interval defined by coords(n-1),coords(n) which always contains BOTH boundary 
!!  values coords(n-1) and coords(n)
!!  The code does NOT test whether the point is at all inside the range of the incoming intervals, but assumes
!!  it is contained in any interval. E.g. if the point is left of the second interval boundary, it is assumed
!!  that the point is in the first interval (without testing the left boundary)
!! \param c coordinate to be located inside intervals defined by values in coords
!! \param coords array of coordinates (x,y,z), is assumed to be either in ascending or descending order!
!! \param n size of coords
!! \param is_ascending indicates whether coords is in ascending order or not
!! \param i always in range 2,...,n indicating the minimal index in the array, for which the array value 
!!        is strictly larger than the value to locate (or i=n if c is larger or equal to coords(n-1)
! 
  function locateCoordinateChunksInversionGrid(c,coords,n,is_ascending) result(i)
    integer :: n
    double precision, dimension(n) :: coords
    double precision :: c
    logical :: is_ascending
    ! returning
    integer :: i
    ! local
    integer :: jl,jr,jm
!
    if(is_ascending) then
!
       if(c < coords(2)) then
          i = 2
          return
       end if
       jl = 2; jr = n
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
       if(c > coords(2)) then
          i = 2
          return
       end if
       jl = 2; jr = n
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
    i = jl+1
  end function locateCoordinateChunksInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of wavefield points contained in cell icell to standard cell and compute their jacobian
!! \param this chunks inversion grid
!! \param icell global inversion grid cell index
!! \param wp-x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param wp-y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param wp-z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param uf_wp unit factor of wavefield points
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param type_standard_cell defines on return the shape of the standard cell (specific integration weights routine can be chosen): (4=Tetrahedron,6=Hexahedron). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param errmsg error message
!
  subroutine transformToStandardCellChunksInversionGrid(this,icell,x,y,z,uf_wp,jacobian,type_standard_cell,errmsg)
    type (chunks_inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    real :: uf_wp
    type (error_message) :: errmsg
    ! local
    character(len=42) :: myname = 'transformToStandardCellChunksInversionGrid'
    character (len=400) :: errstr
    integer :: n,icell_internal
    double precision, dimension(:), allocatable :: xd,yd,zd
    integer, dimension(:), allocatable :: ichunk
    double precision :: dx_half,dy_half,dz_half,cx,cy,cz,uf_correction
!
    call addTrace(errmsg,myname)
    if(.not.this%is_defined) then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    if(type_standard_cell == -1) then
       call add(errmsg,2,"Incoming value of type_standard_cell is -1, indicating a request for total "//&
            "integration weights (e.g. used by integration weights of type 6) instead of jacobian values. "//&
            "This functionality is not supported by inversion grids of type chunksInversionGrid",myname)
       return
    end if
!
    ! assume that incoming icell is in range 1,..,this%ncell (was checked by module inversionGrid)
    icell_internal = this%icell_internal(icell)
!
    cx = 0.5d0 * ( this%x(this%cell_ixmax(icell_internal)) + this%x(this%cell_ixmin(icell_internal)) )
    cy = 0.5d0 * ( this%y(this%cell_iymax(icell_internal)) + this%y(this%cell_iymin(icell_internal)) )
    cz = 0.5d0 * ( this%z(this%cell_izmax(icell_internal)) + this%z(this%cell_izmin(icell_internal)) )
!
    dx_half = 0.5d0 * ( this%x(this%cell_ixmax(icell_internal)) - this%x(this%cell_ixmin(icell_internal)) )
    dy_half = 0.5d0 * ( this%y(this%cell_iymax(icell_internal)) - this%y(this%cell_iymin(icell_internal)) )
    dz_half = 0.5d0 * ( this%z(this%cell_izmax(icell_internal)) - this%z(this%cell_izmin(icell_internal)) )
!
    ! assume that all incoming vectors x,y,z are of same size ( > 0); this is checked in module inversionGrid
    n = size(x)
    allocate(xd(n),yd(n),zd(n),ichunk(n))
    if(uf_wp /= 1.0e3) then
       ! account for wavefield points which are not given in km (e.g. using SPECFEM3D_Cartesian (or NEXD) with SI
       ! units for spherical applications)

       ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
       uf_correction = dble(uf_wp) * 1.0d-3
       xd = dble(x) * uf_correction
       yd = dble(y) * uf_correction
       zd = dble(z) * uf_correction
    else
       xd = dble(x)
       yd = dble(y)
       zd = dble(z)
    end if
    ichunk(:) = this%cell_ichunk(icell_internal)
    call transformVectorGlobalToFlatChunksInversionGrid(this,xd,yd,zd,ichunk,n)
!
    ! transform to hexahedral standard cell [-1,1]^3
    x = real( (xd-cx)/dx_half )
    y = real( (yd-cy)/dy_half )
    z = real( (zd-cz)/dz_half )
!
    ! check if the incoming points were actually located in cell icell and transform correctly into the inside of standard cell icell
    if(any(x<-1.) .or. any(x>1.) .or. any(y<-1.) .or. any(y>1.) .or. any(z<-1.) .or. any(z>1.)) then
       write(errstr,*) "these wavefield points (supposed to be contained in cell ",icell,&
            ") transform to outside the hexahedral standard cell [-1,1]^3"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
    jacobian = real(  dx_half*dy_half*dz_half * (zd*zd)/(xd*xd+yd*yd+1.d0)**1.5  )
!
    type_standard_cell = 6
!
1   deallocate(xd,yd,zd,ichunk)
  end subroutine transformToStandardCellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of inversion grid cell
!! \param this chunks inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param volume volume of cell icell
!
  subroutine getVolumeCellChunksInversionGrid(this,icell,volume)
    type (chunks_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    ! local
    integer :: icell_internal
    double precision :: dx,dy,dz,cx,cy,cz
!
    ! assume that incoming icell is in range 1,..,this%ncell and that this%is_defined (was checked by module inversionGrid)
    icell_internal = this%icell_internal(icell)
!
    cx = 0.5d0 * ( this%x(this%cell_ixmax(icell_internal)) + this%x(this%cell_ixmin(icell_internal)) )
    cy = 0.5d0 * ( this%y(this%cell_iymax(icell_internal)) + this%y(this%cell_iymin(icell_internal)) )
    cz = 0.5d0 * ( this%z(this%cell_izmax(icell_internal)) + this%z(this%cell_izmin(icell_internal)) )
!
    dx = this%x(this%cell_ixmax(icell_internal)) - this%x(this%cell_ixmin(icell_internal))
    dy = this%y(this%cell_iymax(icell_internal)) - this%y(this%cell_iymin(icell_internal))
    dz = this%z(this%cell_izmax(icell_internal)) - this%z(this%cell_izmin(icell_internal))
!
    ! volume = vol_standard_cell(=8)  * jacobian_at_center
    ! for jacobian see above routine transformToStandardCellChunksInversionGrid
    volume = real(  dx*dy*dz * (cz*cz)/(cx*cx+cy*cy+1.d0)**1.5  )
  end subroutine getVolumeCellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get center of inversion grid cell
!! \param this chunks inversion grid
!! \param icell index of inversion grid for which center should be returned
!! \param c1 first coordinate of center of cell icell
!! \param c2 second coordinate of center of cell icell
!! \param c3 third coordinate of center of cell icell
!! \param coords_type 'wp','event','station'; optional request
!
  subroutine getCenterCellChunksInversionGrid(this,icell,c1,c2,c3,coords_type)
    type (chunks_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: c1,c2,c3
    character(len=*), optional :: coords_type
    ! local
    integer :: icell_internal
    double precision :: r,theta,phi
    double precision, dimension(1) :: c1d,c2d,c3d
!
    ! assume that incoming icell is in range 1,..,this%ncell and that this%is_defined (was checked by module inversionGrid)
    icell_internal = this%icell_internal(icell)
!
    ! in array this%cell_center, global Cartesian coordinates of the cell centers are stored (i.e. 'wp'-type coordinates)
    c1d(1) = this%cell_center(1,icell_internal)
    c2d(1) = this%cell_center(2,icell_internal)
    c3d(1) = this%cell_center(3,icell_internal)
!
    ! transform the "flat" coordinates of the cell centers to global Cartesian
    call transformVectorFlatToGlobalChunksInversionGrid(this,c1d,c2d,c3d,(/this%cell_ichunk(icell_internal)/),1)
!
    ! when using spherical inversion grids, the event and station coordinates
    ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
    ! so in case of coords_type = 'event' or 'station', transform the outgoing c1,c2,c3 
    ! further (which at this point are in global Cartesian coordinates, i.e. 'wp'-form)
    select case(coords_type)
    case('wp')
       ! convert to single precision before return
       c1 = real(c1d(1))
       c2 = real(c2d(1))
       c3 = real(c3d(1))
    case('station','event')
       ! r = sqrt(x^2+y^2+z^2)
       ! theta = arccos(z/r) \in [0,pi] (d.h. 90deg-theta rechnen!)
       ! phi = atan2(x,y)
       r = dsqrt( c1d(1)*c1d(1) + c2d(1)*c2d(1) + c3d(1)*c3d(1) )
       theta = acos(c3d(1)/r)
       phi = atan2(c2d(1),c1d(1))
       
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
  end subroutine getCenterCellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of inversion grid cell
!! \param this chunks inversion grid
!! \param icell index of inversion grid for which radius should be returned
!! \param radius radius of cell icell
!
  subroutine getRadiusCellChunksInversionGrid(this,icell,radius)
    type (chunks_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
!
    ! assume that incoming icell is in range 1,..,this%ncell and that this%is_defined (was checked by module inversionGrid)
    radius = this%cell_radius(this%icell_internal(icell))
  end subroutine getRadiusCellChunksInversionGrid
!------------------------------------------------------------------------
!> \brief function answers whether a given point inside the inversion grid domain
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station'
!!  Also it can be assured that this chunks inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here).
!!  in case coords_type is 'wp', c1,c2,c3 are expected to be global cartesian x,y,z coordinates for spherical application
!!  in case coords_type is 'event', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and depth (in km)
!!  in case coords_type is 'station', c1,c2,c3 are expected to be geographical lat (c1), lon (c2) and altitude (in m)
!!  altitude of stations is ignored on return (in order to plot stations on surface of spherical chunk)
!! \param this chunks inversion grid
!! \param c1 first coordinate
!! \param c2 second coordinate
!! \param c3 third coordinate
!! \param coords_type 'wp','event','station'
!! \param uf_wp unit factor of wavefield points; optional (recommended to be used along with 'wp')
!! \param ichunk on return, ichunk contains the index of the chunk in which the point is located; optional
!! \param l logical indicating whether the given point inside the inversion grid domain
!! \return logical indicating whether the given point inside the inversion grid domain
!
  function pointInsideChunksInversionGrid(this,c1,c2,c3,coords_type,uf_wp,ichunk) result(l)
    type (chunks_inversion_grid) :: this
    real, intent(in) :: c1,c2,c3
    character(len=*) :: coords_type
    real, optional :: uf_wp
    integer, optional :: ichunk
    logical :: l
    ! local
    double precision, dimension(1) :: c1_copy,c2_copy,c3_copy
    double precision :: lat,lon,z,uf_correction
    integer, dimension(1) :: wp_ichunk
!
    l = .false.
    if(present(ichunk)) ichunk = -1
!
    if(.not.this%is_defined) return
!
    ! create a double precision copy in an array of dimension 1 (in order to use the point in array routines below)
    c1_copy(1) = dble(c1)
    c2_copy(1) = dble(c2)
    c3_copy(1) = dble(c3)
!
    select case(coords_type)
    case('station','event')
       ! when using spherical inversion grids, the event and station coordinates
       ! are assumed to be given in spherical coordinates in degrees (and depth / altitude, respectively)
       ! so in case of coords_type = 'event' or 'station', transform the incoming c1,c2,c3 to respective
       ! global Cartesian coordinates first, then execute the search operations below (which work on wavefield 
       ! point coordinates)
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
          ! check here, if this depth exceeds the maximum depth of the inversion grid (if so, return immediately)
          if( this%rmax-z < this%z(this%sorted2iz(1)) ) return
!
          c1_copy(1) = c1_copy(1) * (this%rmax-z)
          c2_copy(1) = c2_copy(1) * (this%rmax-z)
          c3_copy(1) = c3_copy(1) * (this%rmax-z)
       end select
!
       ! Now check, if the point is laterally contained in any of the chunks. Therefore, the above scaling
       ! to the correct depth is actually NOT NECESSARY. Stick to it anyway, in order not to assume too much
       ! about the routine locateWpLaterallyInsideChunksInversionGrid.
       call locateWpLaterallyInsideChunksInversionGrid(this,c1_copy,c2_copy,c3_copy,wp_ichunk,1)
       l = wp_ichunk(1) > 0
       if(present(ichunk) .and. l) ichunk = wp_ichunk(1)
!
    case('wp')
       ! if the optional argument uf_wp is given, check if wavefield point coordinates are given in km,
       ! otherwise transform to km first, then execute the tests below
       if(present(uf_wp)) then
          if(uf_wp /= 1.0e3) then
             ! multiplying by uf_wp brings the coordinates to SI units of m and then dividing by 1000 brings them to km
             uf_correction = dble(uf_wp) * 1.0d-3
             c1_copy(1) = c1_copy(1) * uf_correction
             c2_copy(1) = c2_copy(1) * uf_correction
             c3_copy(1) = c3_copy(1) * uf_correction
          end if ! uf_wp /= 1.0e3
       end if ! present(uf_wp)
!
       ! First check if point is in depth range of this inversion Grid
       z = dsqrt( c1_copy(1)**2 + c2_copy(1)**2 + c3_copy(1)**2)
       if(real(z) > this%rmax .or. z < this%z(this%sorted2iz(1)) ) return
!
       ! If the point is in correct depth range, check if it is laterally located inside any of the chunks
       call locateWpLaterallyInsideChunksInversionGrid(this,c1_copy,c2_copy,c3_copy,wp_ichunk,1)
       l = wp_ichunk(1) > 0
       if(present(ichunk) .and. l) ichunk = wp_ichunk(1)
    end select
  end function pointInsideChunksInversionGrid
!
end module chunksInversionGrid
