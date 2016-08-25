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
!> \brief generic inversion grid module which forks to specific inversion grid modules
!!
!! \details ASKI supports several types of inversion grids:
!! 
!!  type_name_inversion_grid = 'scartInversionGrid', module scartInversionGrid
!!     ASKI internal, method independent simple Cartesian inversion grid containing of layeres of 
!!     regular hexahedra, layer-dependent refinement possible
!!
!!  type_name_inversion_grid = 'ecartInversionGrid', module ecartInversionGrid
!!     ASKI internal, external Cartesian inversion grid provided e.g. by CUBIT which may contain tetrahedra,
!!      hexahedra are not fully supported yet
!!
!!  type_name_inversion_grid = 'specfem3dInversionGrid', module specfem3dInversionGrid
!!     method dependent inversion grid for METHOD = SPECFEM3D 
!!     use SPECFEM elements as inversion grid:
!!     use spectral elements as the inversion grid cells and read in the geometry including information 
!!     on neighbour cells, using ALL GLL points in an element as wavefield points, reading in their jacobian
!!
!!  type_name_inversion_grid = 'schunkInversionGrid', module schunkInversionGrid
!!     ASKI internal, method independent spherical inverison grid using one chunk, based on simple Cartesian 
!!     inversion grid, adding curvature projections using the concept of a cubed sphere.
!!     Only suitable for chunks of small width, since some simplifications and approximations are assumed.
!!     For proper spherical application, especially for larger chunks, use chunksInversionGrid, which , 
!!     however, is more expensive to set up than schunkInversionGrid. 
!!
!!  type_name_inversion_grid = 'chunksInversionGrid', module chunksInversionGrid
!!     ASKI internal, method independent spherical inverison grid using (one or several) chunks in 
!!     the concept of a cubed sphere and variable cell refinement. 
!!
!! \author Florian Schumacher
!! \date Okt 2015
!
module inversionGrid
!
  use chunksInversionGrid
  use schunkInversionGrid
  use scartInversionGrid
  use ecartInversionGrid
  use specfem3dInversionGrid
  use errorMessage
!
  implicit none
!
  interface getInvgrid
     module procedure getScartInversionGrid
     module procedure getEcartInversionGrid
     module procedure getSpecfem3dInversionGrid
     module procedure getSchunkInversionGrid
     module procedure getChunksInversionGrid
  end interface getInvgrid
  interface dealloc; module procedure deallocateInversionGrid; end interface
  interface operator (.ncell.); module procedure getNcellInversionGrid; end interface
!
  integer, parameter :: character_length_type_inversion_grid = 22
  character(len=114), parameter :: all_valid_types_inversion_grid = &
       "'scartInversionGrid', 'ecartInversionGrid', 'specfem3dInversionGrid', 'schunkInversionGrid', 'chunksInversionGrid'"
!
  type inversion_grid
     private
     character(len=character_length_type_inversion_grid) :: type_name_inversion_grid = '' !< indicates which pointer below is allocated (can have values as in string all_valid_types_inversion_grid)
     type (scart_inversion_grid), pointer :: scart => null()
     type (ecart_inversion_grid), pointer :: ecart => null()
     type (specfem3d_inversion_grid), pointer :: specfem3d => null()
     type (schunk_inversion_grid), pointer :: schunk => null()
     type (chunks_inversion_grid), pointer :: chunks => null()
  end type inversion_grid
!
contains
!------------------------------------------------------------------------
!> \brief check if some name is a valid inversion grid type, supported by this module
!! \details As there are some forward method dependent types of inversion grids and 
!!  inversion grid dependent types of integration weights, this routine can, optionally
!!  check compatibility of all those types. This routine at least requires a name of
!!  an inversion grid type, the validity of which will be evaluated
!
  function validTypeInversionGrid(name,method,intw_type,err) result(l)
    character(len=*), intent(in) :: name
    logical :: l
    character(len=*), optional :: method
    integer, optional :: intw_type
    type (error_message), optional :: err
    ! local
    character(len=22) :: myname = 'validTypeInversionGrid'
!
    if(present(err)) call addTrace(err,myname)
    l = .true.
!
    select case(name)
    case('scartInversionGrid','ecartInversionGrid','specfem3dInversionGrid','schunkInversionGrid','chunksInversionGrid')
       ! ok, do nothing
    case default
       l = .false.
       if(present(err)) call add(err,2,"incoming type of inversion grid '"//trim(name)//&
            "' is invalid: valid types are "//all_valid_types_inversion_grid,myname)
    end select
!
    ! dependent on type of inversion grid, check if the method is incompatible with this inversion grid
    if(present(method)) then
       select case(name)
       case('specfem3dInversionGrid')
          if(method /= "SPECFEM3D") then
             l = .false.
             if(present(err)) call add(err,2,"type 'specfem3dInversionGrid' must be used with forward "//&
                  "method 'SPECFEM3D' only",myname)
          end if
       end select
    end if
!
    ! Dependent on type of inversion grid, check if the type of integration weights is incompatible with this 
    ! inversion grid.
    ! The only integration weights type to check here, is type 6 (external weights), which is only compatible
    ! with specfem3dInversionGrid. For all other inversion types, return .false.
    if(present(intw_type)) then
       select case(name)
       case('scartInversionGrid','ecartInversionGrid','schunkInversionGrid','chunksInversionGrid')
          if(intw_type == 6) then
             l = .false.
             if(present(err)) call add(err,2,"type 6 integration weights (external) are not supported by "//&
                  "inversion grid types 'scartInversionGrid','ecartInversionGrid','chunksInversionGrid'",myname)
          end if
       end select
    end if
  end function validTypeInversionGrid
!------------------------------------------------------------------------
  function canTransformToVtkPointsOutsideInversionGrid(this,coords_type) result(l)
    type (inversion_grid) :: this
    character(len=*) :: coords_type
    logical :: l
    l = .false.
    if(this%type_name_inversion_grid == '') return
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       l = canTransformToVtkPointsOutsideScartInversionGrid(this%scart,coords_type)
    case('ecartInversionGrid')
       l = canTransformToVtkPointsOutsideEcartInversionGrid(this%ecart,coords_type)
    case('specfem3dInversionGrid')
       l = canTransformToVtkPointsOutsideSpecfem3dInversionGrid(this%specfem3d,coords_type)
    case('schunkInversionGrid')
       l = canTransformToVtkPointsOutsideSchunkInversionGrid(this%schunk,coords_type)
    case('chunksInversionGrid')
       l = canTransformToVtkPointsOutsideChunksInversionGrid(this%chunks,coords_type)
    end select
  end function canTransformToVtkPointsOutsideInversionGrid
!------------------------------------------------------------------------
!> \brief get unit factor of the volume element dependent on the dimensions of the specific inversion grid
!! \param this inversion grid
!! \param uf_wp unit factor of wavefield points
!! \param uf_vol unit factor of volume element (return value of this subroutine)
!! \param errmsg error message
  subroutine getUnitFactorOfVolumeElementInversionGrid(this,uf_wp,uf_vol,errmsg)
    type (inversion_grid) :: this
    real :: uf_wp,uf_vol
    type (error_message) :: errmsg
    character(len=41) :: myname = "getUnitFactorOfVolumeElementInversionGrid"
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
!
    if(this%type_name_inversion_grid == '') then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    if(uf_wp <= 0) then
       write(errstr,*) "incoming unit factor of wavefield points = ",uf_wp,&
            " is not strictly positive. This is not supported by ASKI, there seems to be some problem."
       call add(errmsg,2,errstr,myname)
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call getUnitFactorOfVolumeElementScartInversionGrid(this%scart,uf_wp,uf_vol,errmsg)
    case('ecartInversionGrid')
       call getUnitFactorOfVolumeElementEcartInversionGrid(this%ecart,uf_wp,uf_vol,errmsg)
    case('specfem3dInversionGrid')
       call getUnitFactorOfVolumeElementSpecfem3dInversionGrid(this%specfem3d,uf_wp,uf_vol,errmsg)
    case('schunkInversionGrid')
       ! For schunk inversion grids, uf_wp is not needed to compute the unit factor of the volume element
       call getUnitFactorOfVolumeElementSchunkInversionGrid(this%schunk,uf_vol,errmsg)
    case('chunksInversionGrid')
       ! For chunks inversion grids, uf_wp is not needed to compute the unit factor of the volume element
       call getUnitFactorOfVolumeElementChunksInversionGrid(this%chunks,uf_vol,errmsg)
    end select
  end subroutine getUnitFactorOfVolumeElementInversionGrid
!------------------------------------------------------------------------
!> \brief get pointer to scart_inversion_grid of this inversion_grid
  subroutine getScartInversionGrid(this,p)
    type (inversion_grid), intent(in) :: this
    type (scart_inversion_grid), pointer, intent(out) :: p
    p => this%scart
  end subroutine getScartInversionGrid
!------------------------------------------------------------------------
!> \brief get pointer to ecart_inversion_grid of this inversion_grid
  subroutine getEcartInversionGrid(this,p)
    type (inversion_grid), intent(in) :: this
    type (ecart_inversion_grid), pointer, intent(out) :: p
    p => this%ecart
  end subroutine getEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get pointer to specfem3d_inversion_grid of this inversion_grid
  subroutine getSpecfem3dInversionGrid(this,p)
    type (inversion_grid), intent(in) :: this
    type (specfem3d_inversion_grid), pointer, intent(out) :: p
    p => this%specfem3d
  end subroutine getSpecfem3dInversionGrid
!------------------------------------------------------------------------
!> \brief get pointer to schunk_inversion_grid of this inversion_grid
  subroutine getSchunkInversionGrid(this,p)
    type (inversion_grid), intent(in) :: this
    type (schunk_inversion_grid), pointer, intent(out) :: p
    p => this%schunk
  end subroutine getSchunkInversionGrid
!------------------------------------------------------------------------
!> \brief get pointer to chunks_inversion_grid of this inversion_grid
  subroutine getChunksInversionGrid(this,p)
    type (inversion_grid), intent(in) :: this
    type (chunks_inversion_grid), pointer, intent(out) :: p
    p => this%chunks
  end subroutine getChunksInversionGrid
!------------------------------------------------------------------------
!> \brief create inversion grid of specific type
!! \details Each specific type inversion grid module may use own specifications, 
!!  as defined in parameter file parfile. If recreate is present, it controls if the 
!!  inversion grid should be newly created from specifications in parfile, or may be read 
!!  in from file(s) that were previously created (if existend). Use this flag, if parfile
!!  has changed but this is not recognized by the files previously created.
!! \param this inversion grid
!! \param type_name name of inversion grid type to be created
!! \param parfile filename of parameter file/ any file(s) for creation of specific grid
!! \param path path where to find/write own files (usually path of current iteration step)
!! \param lu file unit to be used for reading and writing files
!! \param errmsg error message
!! \param recreate if set, it controls if inversion grid should be newly created, or can be read in
!
  subroutine createInversionGrid(this,type_name,parfile,path,lu,errmsg,recreate)
    type (inversion_grid) :: this
    character(len=*) :: type_name,parfile,path
    integer :: lu
    type (error_message) :: errmsg
    logical, optional :: recreate
    ! local
    character(len=19) :: myname = 'createInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(.not.validTypeInversionGrid(type_name)) then
       call add(errmsg,2,"inversion grid type '"//trim(type_name)//"' not supported",myname)
       return
    end if
!
    if(this%type_name_inversion_grid /= '') then
       call add(errmsg,1,"inversion grid of type '"//trim(this%type_name_inversion_grid)//&
            "' is already created. deallocating it now, before creating new grid",myname)
       call deallocateInversionGrid(this)
    end if
!
    select case(type_name)
    case('scartInversionGrid')
       allocate(this%scart)
       call createScartInversionGrid(this%scart,parfile,lu,errmsg)
       ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
       if(.level.errmsg/=2) this%type_name_inversion_grid = 'scartInversionGrid'
    case('ecartInversionGrid')
       allocate(this%ecart)
       call createEcartInversionGrid(this%ecart,parfile,path,lu,errmsg,recreate)
       ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
       if(.level.errmsg/=2) this%type_name_inversion_grid = 'ecartInversionGrid'
    case('specfem3dInversionGrid')
       allocate(this%specfem3d)
       call createSpecfem3dInversionGrid(this%specfem3d,parfile,path,lu,errmsg)
       ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
       if(.level.errmsg/=2) this%type_name_inversion_grid = 'specfem3dInversionGrid'
    case('schunkInversionGrid')
       allocate(this%schunk)
       call createSchunkInversionGrid(this%schunk,parfile,lu,errmsg)
       ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
       if(.level.errmsg/=2) this%type_name_inversion_grid = 'schunkInversionGrid'
    case('chunksInversionGrid')
       allocate(this%chunks)
       call createChunksInversionGrid(this%chunks,parfile,path,lu,errmsg,recreate)
       ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
       if(.level.errmsg/=2) this%type_name_inversion_grid = 'chunksInversionGrid'
    end select
  end subroutine createInversionGrid
!------------------------------------------------------------------------
!> \brief deallocate inversion grid object
!
  subroutine deallocateInversionGrid(this)
    type (inversion_grid) :: this
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call deallocateScartInversionGrid(this%scart)
       deallocate(this%scart)
    case('ecartInversionGrid')
       call deallocateEcartInversionGrid(this%ecart)
       deallocate(this%ecart)
    case('specfem3dInversionGrid')
       call deallocateSpecfem3dInversionGrid(this%specfem3d)
       deallocate(this%specfem3d)
    case('schunkInversionGrid')
       call deallocateSchunkInversionGrid(this%schunk)
       deallocate(this%schunk)
    case('chunksInversionGrid')
       call deallocateChunksInversionGrid(this%chunks)
       deallocate(this%chunks)
    end select
    this%type_name_inversion_grid = ''
  end subroutine deallocateInversionGrid
!------------------------------------------------------------------------
!> \brief get total number of inversion grid cells
!
  function getNcellInversionGrid(this) result(ncell)
    type (inversion_grid), intent(in) :: this
    integer :: ncell
!
    ncell = 0
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       ncell =  getNcellScartInversionGrid(this%scart)
    case('ecartInversionGrid')
       ncell =  getNcellEcartInversionGrid(this%ecart)
    case('specfem3dInversionGrid')
       ncell =  getNcellSpecfem3dInversionGrid(this%specfem3d)
    case('schunkInversionGrid')
       ncell =  getNcellSchunkInversionGrid(this%schunk)
    case('chunksInversionGrid')
       ncell =  getNcellChunksInversionGrid(this%chunks)
    end select
  end function getNcellInversionGrid
!------------------------------------------------------------------------
!> \brief transform coordinates (of certain kind) into vtk inversion grid coordinates
!! \param this inversion grid
!! \param c1 vector of first coordinate (contains vtk x-values on exit)
!! \param c2 vector of second coordinate (contains vtk y-values on exit)
!! \param c3 vector of third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!! \param uf_wp optional unit of wavefield points (recommended to be used along with coords_type == 'wp')
!
  subroutine transformToVtkInversionGrid(this,c1,c2,c3,coords_type,errmsg,uf_wp)
    type (inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
    real, optional :: uf_wp
    ! local
    character (len=400) :: errstr
    character (len=27) :: myname = 'transformToVtkInversionGrid'
    integer :: size_c1
!
    call addTrace(errmsg,myname)
!
    select case(coords_type)
    case('wp','event','station')
       ! ok, do nothing
    case default
       call add(errmsg,2,"incoming coordinate type '"//trim(coords_type)//&
            "' not supported: must be one of 'wp','event','station'",myname)
       return
    end select
!
    if(size(c1)==0 .or. size(c2)==0 .or. size(c3)==0) then
       call add(errmsg,2,"there are no incoming coordinates",myname)
       return
    end if
!
    size_c1 = size(c1)
    if(size(c2) /= size_c1 .or. size(c3) /= size_c1) then
       call add(errmsg,2,"coordinate vectors must all be of same size",myname)
       return
    end if
!
    if(present(uf_wp)) then
       if(uf_wp <= 0) then
          write(errstr,*) "incoming unit factor of wavefield points = ",uf_wp,&
               " is not strictly positive. This is not supported by ASKI, there seems to be some problem."
          call add(errmsg,2,errstr,myname)
       end if
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       ! For scart inversion grids, uf_wp is not required to transform wavefield points to vtk coordinates
       call transformToVtkScartInversionGrid(this%scart,c1,c2,c3,coords_type,errmsg)
    case('ecartInversionGrid')
       ! For ecart inversion grids, uf_wp is not required to transform wavefield points to vtk coordinates
       call transformToVtkEcartInversionGrid(this%ecart,c1,c2,c3,coords_type,errmsg)
    case('specfem3dInversionGrid')
       ! For specfem3d inversion grids, uf_wp is not required to transform wavefield points to vtk coordinates
       call transformToVtkSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,coords_type,errmsg)
    case('schunkInversionGrid')
       call transformToVtkSchunkInversionGrid(this%schunk,c1,c2,c3,coords_type,errmsg,uf_wp)
    case('chunksInversionGrid')
       call transformToVtkChunksInversionGrid(this%chunks,c1,c2,c3,coords_type,errmsg,uf_wp)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine transformToVtkInversionGrid
!------------------------------------------------------------------------
!> \brief get inversion grid geometry information for vtk files
!! \details Dependent on the type of inversion grid cells (i.e. hexahedra, tetrahedra, ...)
!!  define an array of point coordinates and an array defining cells of a
!!  certain cell_type as in vtk file format for unstructured grid (for geometry type 1 = cell center points, only define the point array).
!!  We need the additional index mapping indx_map here (confer invgrid_vtk_file%req_indx) for the case of 
!!  present(cell_indx_req), as the specific subroutines below may throw away invalid and duplicate 
!!  invgrid cell indices. in order to be able to still use the cell geometry information with data vectors 
!!  of same length (and order) as cell_indx_req, indx_map maps the vtk cell index of the returned vtk cells 
!!  to the original position in array cell_indx_req. also, if a specific subroutine below (for some type of invgrid) does not 
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
!! \param cell_type_requested optional string for requesting a specific cell type only. At the moment only for chunksInversionGrid the values 'EXTERNAL_CELLS' and 'BASE_CELLS' are supported
!
  subroutine getGeometryVtkInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,errmsg,&
       cell_indx_req,indx_map,cell_type_requested)
    type (inversion_grid) :: this
    integer, dimension(:), optional :: cell_indx_req
    character(len=*), optional :: cell_type_requested
    ! outgoing
    integer :: geometry_type
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map
    type (error_message) :: errmsg
    ! local
    character(len=27) :: myname = 'getGeometryVtkInversionGrid'
!
    call addTrace(errmsg,myname)
    nullify(points,cell_connectivity,cell_type,cell_indx_out)
    if(present(indx_map)) nullify(indx_map)
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       if(present(cell_type_requested)) call add(errmsg,1,"inversion grids of type 'scartInversionGrid' do not support the "//&
            "presence of optional argument cell_type_requested: it is ignored here",myname)
       call getGeometryVtkScartInversionGrid(this%scart,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
            errmsg,cell_indx_req,indx_map)

    case('ecartInversionGrid')
       if(present(cell_type_requested)) call add(errmsg,1,"inversion grids of type 'ecartInversionGrid' do not support the "//&
            "presence of optional argument cell_type_requested: it is ignored here",myname)
       call getGeometryVtkEcartInversionGrid(this%ecart,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
            errmsg,cell_indx_req,indx_map)

    case('specfem3dInversionGrid')
       if(present(cell_type_requested)) call add(errmsg,1,"inversion grids of type 'specfem3dInversionGrid' do not support the "//&
            "presence of optional argument cell_type_requested: it is ignored here",myname)
       call getGeometryVtkSpecfem3dInversionGrid(this%specfem3d,geometry_type,points,cell_connectivity,cell_type,&
            cell_indx_out,errmsg,cell_indx_req,indx_map)

    case('schunkInversionGrid')
       if(present(cell_type_requested)) call add(errmsg,1,"inversion grids of type 'schunkInversionGrid' do not support the "//&
            "presence of optional argument cell_type_requested: it is ignored here",myname)
       call getGeometryVtkSchunkInversionGrid(this%schunk,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
            errmsg,cell_indx_req,indx_map)

    case('chunksInversionGrid')
       if(present(cell_type_requested)) then
          select case(cell_type_requested)
          case('EXTERNAL_CELLS')
             call getGeometryVtkChunksInversionGrid(this%chunks,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
                  errmsg,cell_indx_req,indx_map)
          case('BASE_CELLS')
             call getBaseCellGeometryVtkChunksInversionGrid(this%chunks,geometry_type,points,cell_connectivity,cell_type,&
                  cell_indx_out,errmsg,cell_indx_req,indx_map)
          case default
             call add(errmsg,2,"inversion grids of type chunksInversionGrid do not support a cell type request '"//&
                  trim(cell_type_requested)//"', only 'EXTERNAL' and 'BASE_CELLS' are supported",myname)
             return
          end select ! cell_type_requested
       else ! present(cell_type_requested)
          call getGeometryVtkChunksInversionGrid(this%chunks,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
               errmsg,cell_indx_req,indx_map)
       end if ! present(cell_type_requested)

    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select ! this%type_name_inversion_grid
  end subroutine getGeometryVtkInversionGrid
!------------------------------------------------------------------------
!> \brief get for all cells the indices of their face neighbour cells
!!  \details If boundary_conditions is set, there should be a special form of nb_idx returned, which is used
!!  e.g. to define special smoothing condigions:
!!   'no_nb_inner_bnd': with some neighbours removed (which should not be smoothed with on internal invgrid 
!!                      boundaries). FOR ALL INVERSIONGRIDS, THERE IS NO POSSIBILITY YET TO DEFINE INTERNAL
!!                      INVGRID BOUNDARIES!! SO IF THIS IS VALUE IS SET, NO NEIGHBOURS ARE RETURNED! IN THE FUTURE
!!                      MIGHT INTRODUCE ACTUAL FUNCTIONALITY IN THE SUBMODULES
!!   'extra_nbs_outer_bnd': additional fake neighbours (having cell index 0) in case of zero boundary conditions 
!!                          on outer invgrid boundaries.
!!   'extra_nbs_outer_bnd_except_free_surface': same as 'extra_nbs_outer_bnd', but not applied for cells on free surfaces.
!!   '','standard': no special boundary handling, i.e. standard functionality (exactly the geometrical face neighbours are given)
!! \param this inversion grid
!! \param nb_idx pointer to array of length .ncell.this defining neighbours; if invgrid not defined yet, nullified on exit
!! \param boundary_conditions optional string indicating type of boundary conditions for which neighbours should be returned
!
  subroutine getIndicesFaceNeighboursInversionGrid(this,nb_idx,boundary_conditions)
    type (inversion_grid), intent(in) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
    character(len=*), optional :: boundary_conditions
!
    nullify(nb_idx)
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call getIndicesFaceNeighboursScartInversionGrid(this%scart,nb_idx,boundary_conditions)
    case('ecartInversionGrid')
       ! For ecart inversion grids, boundary condition functionality is not supported. Indicate this
       ! by returning nullified nb_idx. Allow for values boundary_conditions == '','standard' which is defined 
       ! to be equivalent to .not.present(boundary_conditions)
       if(present(boundary_conditions)) then
          select case(boundary_conditions)
          case('','standard') ! OK do nothing, this is supported standard functionality
          case default; return
          end select
       end if
       call getIndicesFaceNeighboursEcartInversionGrid(this%ecart,nb_idx)
    case('specfem3dInversionGrid')
       ! For specfem3d inversion grids, boundary condition functionality is not supported. Indicate this
       ! by returning nullified nb_idx. Allow for values boundary_conditions == '','standard' which is defined 
       ! to be equivalent to .not.present(boundary_conditions)
       if(present(boundary_conditions)) then
          select case(boundary_conditions)
          case('','standard') ! OK do nothing, this is supported standard functionality
          case default; return
          end select
       end if
       call getIndicesFaceNeighboursSpecfem3dInversionGrid(this%specfem3d,nb_idx)
    case('schunkInversionGrid')
       call getIndicesFaceNeighboursSchunkInversionGrid(this%schunk,nb_idx,boundary_conditions)
    case('chunksInversionGrid')
       call getIndicesFaceNeighboursChunksInversionGrid(this%chunks,nb_idx,boundary_conditions)
    end select
  end subroutine getIndicesFaceNeighboursInversionGrid
!------------------------------------------------------------------------
!> \brief get for all cells the wavefield point indices of points contained in that cell
!! \param this inversion grid
!! \param c1 vector of first coordinate of wavefield points
!! \param c2 vector of second coordinate of wavefield points
!! \param c3 vector of third coordinate of wavefield points
!! \param uf_wp unit factor of wavefield points
!! \param wp_idx pointer to array of length .ncell.this; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideInversionGrid(this,c1,c2,c3,uf_wp,wp_idx,errmsg)
    type (inversion_grid) :: this
    real, dimension(:), intent(in) :: c1,c2,c3
    real :: uf_wp
    type (integer_vector_pointer), dimension(:), pointer :: wp_idx
    type (error_message) :: errmsg
    ! local
    character(len=27) :: myname = 'locateWpInsideInversionGrid'
    integer :: size_c1
!
    nullify(wp_idx)
!
    if(size(c1)==0 .or. size(c2)==0 .or. size(c3)==0) then
       call add(errmsg,2,"there are no incoming wavefield point coordinates",myname)
       return
    end if
!
    size_c1 = size(c1)
    if(size(c2) /= size_c1 .or. size(c3) /= size_c1) then
       call add(errmsg,2,"coordinate vectors must all be of same size",myname)
       return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       ! For scart inversion grids, uf_wp is not required to locate wavefield points inside inversion grid cells
       call locateWpInsideScartInversionGrid(this%scart,c1,c2,c3,wp_idx,errmsg)
    case('ecartInversionGrid')
       ! For ecart inversion grids, uf_wp is not required to locate wavefield points inside inversion grid cells
       call locateWpInsideEcartInversionGrid(this%ecart,c1,c2,c3,wp_idx,errmsg)
    case('specfem3dInversionGrid')
       ! For specfem3d inversion grids, uf_wp is not required to locate wavefield points inside inversion grid cells
       call locateWpInsideSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,wp_idx,errmsg)
    case('schunkInversionGrid')
       call locateWpInsideSchunkInversionGrid(this%schunk,c1,c2,c3,uf_wp,wp_idx,errmsg)
    case('chunksInversionGrid')
       call locateWpInsideChunksInversionGrid(this%chunks,c1,c2,c3,uf_wp,wp_idx,errmsg)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine locateWpInsideInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in cell icell to standard cell and compute their jacobian
!! \param this inversion grid
!! \param icell index of inversion grid cell which contains wavefield points c1,c2,c3
!! \param c1 vector of first coordinate (contains x-values in standard cell on exit)
!! \param c2 vector of second coordinate (contains y-values in standard cell on exit)
!! \param c3 vector of third coordinate (contains z-values in standard cell on exit)
!! \param uf_wp unit factor of wavefield points
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param type_standard_cell defines on return the shape of the standard cell (specific integration weights routine can be chosen): (4=Tetrahedron,6=Hexahedron). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param errmsg error message
!
  subroutine transformToStandardCellInversionGrid(this,icell,c1,c2,c3,uf_wp,jacobian,type_standard_cell,errmsg)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: c1,c2,c3,jacobian
    real :: uf_wp
    type (error_message) :: errmsg
    ! local
    character (len=36) :: myname = 'transformToStandardCellInversionGrid'
    character(len=400) :: errstr
    integer :: size_c1,ncell
!
    call addTrace(errmsg,myname)
!
    if(this%type_name_inversion_grid == '') then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    if(uf_wp <= 0) then
       write(errstr,*) "incoming unit factor of wavefield points = ",uf_wp,&
            " is not strictly positive. This is not supported by ASKI, there seems to be some problem."
       call add(errmsg,2,errstr,myname)
    end if
!
    ncell = getNcellInversionGrid(this)
    if(icell<1 .or. icell>ncell) then
       write(errstr,*) "incoming inversion grid cell index ",icell," is out of range: must be between 1 and ",ncell
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    if(size(c1)==0 .or. size(c2)==0 .or. size(c3)==0) then
       call add(errmsg,2,"there are no incoming wavefield point coordinates",myname)
       return
    end if
!
    size_c1 = size(c1)
    if(size(c2) /= size_c1 .or. size(c3) /= size_c1) then
       call add(errmsg,2,"coordinate vectors must all be of same size",myname)
       return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       ! for scart inversion grids, the unit factor of the wavefield points is not relevant for transforming
       ! wavefield points to standard cells, since the inversion grid does not assume specific units for its 
       ! spatial extension but simply assumes that they are the same as for the wavefield point coordinates, 
       ! which are located inside the inversion grid cells only according to their pure numerical values. 
       call transformToStandardCellScartInversionGrid(this%scart,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('ecartInversionGrid')
       ! for ecart inversion grids, the unit factor of the wavefield points is not relevant for transforming
       ! wavefield points to standard cells, since the inversion grid does not assume specific units for its 
       ! spatial extension but simply assumes that they are the same as for the wavefield point coordinates, 
       ! which are located inside the inversion grid cells only according to their pure numerical values. 
       call transformToStandardCellEcartInversionGrid(this%ecart,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('specfem3dInversionGrid')
       ! For specfem3d inversion grids, the unit factor of the wavefield points is not relevant for transforming
       ! wavefield points to standard cells, since the points defining a cell are (a subset of) the wavefield points
       call transformToStandardCellSpecfem3dInversionGrid(this%specfem3d,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('schunkInversionGrid')
       call transformToStandardCellSchunkInversionGrid(this%schunk,icell,c1,c2,c3,uf_wp,jacobian,type_standard_cell,errmsg)
    case('chunksInversionGrid')
       call transformToStandardCellChunksInversionGrid(this%chunks,icell,c1,c2,c3,uf_wp,jacobian,type_standard_cell,errmsg)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine transformToStandardCellInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param volume volume of cell icell
!! \param errmsg error message
!
  subroutine getVolumeCellInversionGrid(this,icell,volume,errmsg)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    type (error_message) :: errmsg
    ! local
    character (len=26) :: myname = 'getVolumeCellInversionGrid'
    character(len=400) :: errstr
    integer :: ncell
!
    call addTrace(errmsg,myname)
!
    if(this%type_name_inversion_grid == '') then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    ncell = getNcellInversionGrid(this)
    if(icell<1 .or. icell>ncell) then
       write(errstr,*) "incoming inversion grid cell index ",icell," is out of range: must be between 1 and ",ncell
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call getVolumeCellScartInversionGrid(this%scart,icell,volume)
    case('ecartInversionGrid')
       call getVolumeCellEcartInversionGrid(this%ecart,icell,volume,errmsg)
    case('specfem3dInversionGrid')
       call getVolumeCellSpecfem3dInversionGrid(this%specfem3d,icell,volume)
    case('schunkInversionGrid')
       call getVolumeCellSchunkInversionGrid(this%schunk,icell,volume)
    case('chunksInversionGrid')
       call getVolumeCellChunksInversionGrid(this%chunks,icell,volume)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine getVolumeCellInversionGrid
!------------------------------------------------------------------------
!> \brief get center point of inversion grid cell
!! \details this point should be given in coordinates c1,c2,c3 as understood as
!!  wavefield point or event or station coordinates. It should represent some sort of point center
!!  of inversion grid cell icell, w.r.t. which quantities on cells (e.g. model
!!  values, kernel values etc.) will be assigned to in space, if a point 
!!  description of the respective quantity is required.
!! \param this inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param c1 first coordinate of center of cell icell
!! \param c2 second coordinate of center of cell icell
!! \param c3 third coordinate of center of cell icell
!! \param errmsg error message
!! \param coords_type 'wp','event','station'; optional request
!
  subroutine getCenterCellInversionGrid(this,icell,c1,c2,c3,errmsg,coords_type)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    real :: c1,c2,c3
    type (error_message) :: errmsg
    character(len=*), optional :: coords_type
    ! local
    character (len=26) :: myname = 'getCenterCellInversionGrid'
    character(len=400) :: errstr
    character(len=7) :: coords_type_tested
    integer :: ncell
!
    call addTrace(errmsg,myname)
!
    if(present(coords_type)) then
       select case(coords_type)
       case('wp','event','station')
          ! ok, do nothing
       case default
          call add(errmsg,2,"incoming coordinate type '"//trim(coords_type)//&
               "' not supported: must be one of 'wp','event','station'",myname)
          return
       end select
       coords_type_tested = coords_type
    else
       coords_type_tested = 'wp'
    end if
!
    if(this%type_name_inversion_grid == '') then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    ncell = getNcellInversionGrid(this)
    if(icell<1 .or. icell>ncell) then
       write(errstr,*) "incoming inversion grid cell index ",icell," is out of range: must be between 1 and ",ncell
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call getCenterCellScartInversionGrid(this%scart,icell,c1,c2,c3)
    case('ecartInversionGrid')
       call getCenterCellEcartInversionGrid(this%ecart,icell,c1,c2,c3)
    case('specfem3dInversionGrid')
       call getCenterCellSpecfem3dInversionGrid(this%specfem3d,icell,c1,c2,c3,coords_type_tested)
    case('schunkInversionGrid')
       call getCenterCellSchunkInversionGrid(this%schunk,icell,c1,c2,c3,coords_type_tested)
    case('chunksInversionGrid')
       call getCenterCellChunksInversionGrid(this%chunks,icell,c1,c2,c3,coords_type_tested)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine getCenterCellInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of inversion grid cell
!! \details this value should represent some sort of radial expansion (as seen from 
!!  the cell center). it is used to interpolate quantities on cells (e.g. model
!!  values, kernel values etc.) to points in space and should, henceforth, describe
!!  the spacial range of influence of cell icell
!! \param icell index of inversion grid for which radius should be returned
!! \param radius radius of cell icell
!! \param errmsg error message
!
  subroutine getRadiusCellInversionGrid(this,icell,radius,errmsg)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
    type (error_message) :: errmsg
    ! local
    character (len=26) :: myname = 'getRadiusCellInversionGrid'
    character(len=400) :: errstr
    integer :: ncell
!
    call addTrace(errmsg,myname)
!
    if(this%type_name_inversion_grid == '') then
       call add(errmsg,2,"inversion grid not yet defined",myname)
       return
    end if
!
    ncell = getNcellInversionGrid(this)
    if(icell<1 .or. icell>ncell) then
       write(errstr,*) "incoming inversion grid cell index ",icell," is out of range: must be between 1 and ",ncell
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       call getRadiusCellScartInversionGrid(this%scart,icell,radius)
    case('ecartInversionGrid')
       call getRadiusCellEcartInversionGrid(this%ecart,icell,radius)
    case('specfem3dInversionGrid')
       call getRadiusCellSpecfem3dInversionGrid(this%specfem3d,icell,radius)
    case('schunkInversionGrid')
       call getRadiusCellSchunkInversionGrid(this%schunk,icell,radius)
    case('chunksInversionGrid')
       call getRadiusCellChunksInversionGrid(this%chunks,icell,radius)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine getRadiusCellInversionGrid
!------------------------------------------------------------------------
!> \brief function answers whether a given point inside the inversion grid domain
!! \param c1 first coordinate
!! \param c2 second coordinate
!! \param c3 third coordinate
!! \param coords_type 'wp','event','station'
!! \param uf_wp unit factor of wavefield points; optional (recommended to be used along with 'wp')
!! \param ichunk on return, ichunk contains the index of the chunk in which the point is located; optional (only sensible inversion grids that have "chunks")
!
  function pointInsideInversionGrid(this,c1,c2,c3,coords_type,uf_wp,ichunk) result(l)
    type (inversion_grid) :: this
    real, intent(in) :: c1,c2,c3
    character(len=*) :: coords_type
    real, optional :: uf_wp
    integer, optional :: ichunk
    logical :: l
!
    l = .false.
    if(present(ichunk)) ichunk = -1
!
    select case(coords_type)
    case('wp','event','station')
       ! ok, do nothing
    case default
       return
    end select
!
    if(present(uf_wp)) then
       if(uf_wp <= 0) return
    end if
!
    select case(this%type_name_inversion_grid)
    case('scartInversionGrid')
       ! For scart inversion grids, uf_wp is not required to check if a point is inside the inversion grid.
       ! The scart inversion grid does not assume specific units for its spatial extension but simply 
       ! assumes that they are the same as for the wavefield point coordinates (which are located inside the
       ! inversion grid cells only according to their pure numerical values). 
       l = pointInsideScartInversionGrid(this%scart,c1,c2,c3,coords_type)
    case('ecartInversionGrid')
       !l = pointInsideEcartInversionGrid(this%ecart,c1,c2,c3,coords_type,uf_wp) ! this routine does not yet exist
    case('specfem3dInversionGrid')
       !l = pointInsideSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,coords_type,uf_wp) ! this routine does not yet exist
    case('schunkInversionGrid')
       l = pointInsideSchunkInversionGrid(this%schunk,c1,c2,c3,coords_type,uf_wp)
       if(present(ichunk)) ichunk = 1
    case('chunksInversionGrid')
       l = pointInsideChunksInversionGrid(this%chunks,c1,c2,c3,coords_type,uf_wp,ichunk)
    case default
       return
    end select
  end function pointInsideInversionGrid
!!$!------------------------------------------------------------------------
!!$!> \brief transform coordinates (of certain kind) into vtk inversion grid coordinates
!!$!! \param this inversion grid
!!$!! \param c1 vector of first coordinate (contains vtk x-values on exit)
!!$!! \param c2 vector of second coordinate (contains vtk y-values on exit)
!!$!! \param c3 vector of third coordinate (contains vtk z-values on exit)
!!$!! \param coords_type 'wp','event','station'
!!$!
!!$  subroutine transformCoordinatesInversionGrid(this,c1,c2,c3,coords_type_in,coords_type_out,errmsg)
!!$    type (inversion_grid) :: this
!!$    real, dimension(:), intent(inout) :: c1,c2,c3
!!$    character(len=*) :: coords_type_in,coords_type_out
!!$    type (error_message) :: errmsg
!!$    ! local
!!$    character (len=33) :: myname = 'transformCoordinatesInversionGrid'
!!$    integer :: size_c1
!!$!
!!$    call addTrace(errmsg,myname)
!!$!
!!$    select case(coords_type_in)
!!$    case('wp','event','station')
!!$       ! ok, do nothing
!!$    case default
!!$       call add(errmsg,2,"incoming coordinate type '"//trim(coords_type)//&
!!$            "' not supported: must be one of 'wp','event','station'",myname)
!!$       return
!!$    end select
!!$!
!!$    select case(coords_type_out)
!!$    case('wp','event','station')
!!$       ! ok, do nothing
!!$    case default
!!$       call add(errmsg,2,"requested outgoing coordinate type '"//trim(coords_type)//&
!!$            "' not supported: must be one of 'wp','event','station'",myname)
!!$       return
!!$    end select
!!$!
!!$    if(size(c1)==0 .or. size(c2)==0 .or. size(c3)==0) then
!!$       call add(errmsg,2,"there are no incoming coordinates",myname)
!!$       return
!!$    end if
!!$!
!!$    size_c1 = size(c1)
!!$    if(size(c2) /= size_c1 .or. size(c3) /= size_c1) then
!!$       call add(errmsg,2,"coordinate vectors must all be of same size",myname)
!!$       return
!!$    end if
!!$!
!!$    select case(this%type_name_inversion_grid)
!!$    !case('chunksInversionGrid')
!!$    !   call transformCoordinatesChunksInversionGrid(this%chunks,c1,c2,c3,coords_type,errmsg)
!!$    case('scartInversionGrid')
!!$       call transformCoordinatesScartInversionGrid(this%scart,c1,c2,c3,coords_type,errmsg)
!!$    case('ecartInversionGrid')
!!$       call transformCoordinatesEcartInversionGrid(this%ecart,c1,c2,c3,coords_type,errmsg)
!!$    case('specfem3dInversionGrid')
!!$       call transformCoordinatesSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,coords_type,errmsg)
!!$    case('schunkInversionGrid')
!!$       call transformCoordinatesSchunkInversionGrid(this%schunk,c1,c2,c3,coords_type,errmsg)
!!$    case default
!!$       call add(errmsg,2,"inversion grid not yet defined",myname)
!!$    end select
!!$  end subroutine transformCoordinatesInversionGrid
!
end module inversionGrid
