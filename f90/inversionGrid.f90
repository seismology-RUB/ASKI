!----------------------------------------------------------------------------
!   Copyright 2013 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 0.3.
!
!   ASKI version 0.3 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 0.3 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 0.3.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief generic inversion grid module which forks to specific inversion grid modules
!!
!! \details ASKI supports several types of inversion grids for FORWARD_METHOD = SPECFEM3D:
!! 
!!  type_name_inversion_grid = 'ccsInversionGrid', module ccsInversionGrid -> NOT SUPPORTED YET!!
!!     ASKI internal, method independent spherical inverison grid using (one or several) chunks in 
!!     the concept of a cubed sphere
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
!!  For FORWARD_METHOD = GEMINI (to be incorporated properly in the futuer), the following inversion grid
!!  types are supported:
!!  ...
!!
!! \author Florian Schumacher
!! \date Jul 2013
!
module inversionGrid
!
  !use ccsInversionGrid ! TODO: module for spherical inversion grid based on the chunk-cubed-sphere concept
  use scartInversionGrid
  use ecartInversionGrid
  use specfem3dInversionGrid
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateInversionGrid; end interface
  interface operator (.ncell.); module procedure getNcellInversionGrid; end interface
!
  integer, parameter :: character_length_type_inversion_grid = 22
  character(len=68), parameter :: all_valid_types_inversion_grid = &
       "'scartInversionGrid', 'ecartInversionGrid', 'specfem3dInversionGrid'"
  !character(len=88), parameter :: all_valid_types_inversion_grid = "'ccsInversionGrid', 'scartInversionGrid', 'ecartInversionGrid', 'specfem3dInversionGrid'"
!
  type inversion_grid
     private
     character(len=character_length_type_inversion_grid) :: type_name_inversion_grid = '' !< indicates which pointer below is allocated, compare getTypeNameInversionGrid
     !type (ccs_inversion_grid), pointer :: ccs => null()
     type (scart_inversion_grid), pointer :: scart => null()
     type (ecart_inversion_grid), pointer :: ecart => null()
     type (specfem3d_inversion_grid), pointer :: specfem3d => null()
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
    case('scartInversionGrid','ecartInversionGrid','specfem3dInversionGrid')!,'ccsInversionGrid')
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
    ! dependent on type of inversion grid, check if the type of integration weights is incompatible with this inversion grid
    if(present(intw_type)) then
       select case(name)
       case('scartInversionGrid','ecartInversionGrid','ccsInversionGrid')
          if(intw_type == 6) then
             l = .false.
             if(present(err)) call add(err,2,"type 6 integration weights (external) are not supported by "//&
                  "inversion grid types 'scartInversionGrid','ecartInversionGrid','ccsInversionGrid'",myname)
          end if
       end select
    end if
  end function validTypeInversionGrid
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
    !case('ccsInversionGrid')
    !   allocate(this%ccs)
    !   call createCcsInversionGrid(this%ccs,parfile,lu,errmsg,recreate)
    !   ! only if creation was not erroneous, indicate valid creation by setting this%type_name_inversion_grid
    !   if(.level.errmsg/=2) this%type_name_inversion_grid = 'ccsInversionGrid'
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
    end select
  end subroutine createInversionGrid
!------------------------------------------------------------------------
!> \brief deallocate inversion grid object
!
  subroutine deallocateInversionGrid(this)
    type (inversion_grid) :: this
    select case(this%type_name_inversion_grid)
    !case('ccsInversionGrid')
    !   call deallocateCcsInversionGrid(this%ccs)
    !   deallocate(this%ccs)
    case('scartInversionGrid')
       call deallocateScartInversionGrid(this%scart)
       deallocate(this%scart)
    case('ecartInversionGrid')
       call deallocateEcartInversionGrid(this%ecart)
       deallocate(this%ecart)
    case('specfem3dInversionGrid')
       call deallocateSpecfem3dInversionGrid(this%specfem3d)
       deallocate(this%specfem3d)
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
    !case('ccsInversionGrid')
    !   ncell =  getNcellCcsInversionGrid(this%ccs)
    case('scartInversionGrid')
       ncell =  getNcellScartInversionGrid(this%scart)
    case('ecartInversionGrid')
       ncell =  getNcellEcartInversionGrid(this%ecart)
    case('specfem3dInversionGrid')
       ncell =  getNcellSpecfem3dInversionGrid(this%specfem3d)
    end select
  end function getNcellInversionGrid
!------------------------------------------------------------------------
!> \brief transform coordinates (of certain kind) into vtk inversion grid coordinates
!! \param this inversion grid
!! \param c1 vector of first coordinate (contains vtk x-values on exit)
!! \param c2 vector of second coordinate (contains vtk y-values on exit)
!! \param c3 vector of third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!
  subroutine transformToVtkInversionGrid(this,c1,c2,c3,coords_type,errmsg)
    type (inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
    ! local
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
    select case(this%type_name_inversion_grid)
    !case('ccsInversionGrid')
    !   call transformToVtkCcsInversionGrid(this%ccs,c1,c2,c3,coords_type,errmsg)
    case('scartInversionGrid')
       call transformToVtkScartInversionGrid(this%scart,c1,c2,c3,coords_type,errmsg)
    case('ecartInversionGrid')
       call transformToVtkEcartInversionGrid(this%ecart,c1,c2,c3,coords_type,errmsg)
    case('specfem3dInversionGrid')
       call transformToVtkSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,coords_type,errmsg)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine transformToVtkInversionGrid
!------------------------------------------------------------------------
!> \brief get inversion grid geometry information for vtk files
!! \details Dependent on the type of inversion grid cells (i.e. hexahedra, tetrahedra, ...)
!!  define an array of point coordinates and an array defining cells of a
!!  certain cell_type as in vtk file format for unstructured grid.
!!  We need the additional index mapping indx_map here (confer invgrid_vtk_file%req_indx) for the case of 
!!  present(cell_indx_req), as the specific subroutines belowroutine may throw away invalid and duplicate 
!!  invgrid cell indices. in order to be able to still use the cell geometry information with data vectors 
!!  of same length (and order) as cell_indx_req, indx_map maps the vtk cell index of the returned vtk cells 
!!  to the original position in array cell_indx_req. also, if routine getGeometryVtkInversionGrid does not 
!!  preserve the order of vtk cells as requested by cell_indx_req, map indx_map will reconstruct this order 
!!  (and be the identity otherwise)
!! \param this inversion grid
!! \param points array of nodes which cell_connectivity array refers to
!! \param cell_connectivity array contianing indices of points (indices having zero offset) as required by vtk format
!! \param cell_types vtk cell types of cells as of vtk convention
!! \param cell_indx_out array of same size as the number of vtk cells (same order), which contains the invgrid cell index of a corresponding vtk cell
!! \param cell_indx_req optional incoming array defining invgrid cell indices for which vtk cells should be returned only
!! \param indx_map if cell_indx_req present, then indx_map(cell_indx_req(i)) = i for all valid and non duplicate cell_indx_req(i) , otherwise identity
!
  subroutine getGeometryVtkInversionGrid(this,points,cell_connectivity,cell_type,cell_indx_out,errmsg,cell_indx_req,indx_map)
    type (inversion_grid) :: this
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map
    type (error_message) :: errmsg
    ! local
    character(len=27) :: myname = 'getGeometryVtkInversionGrid'
!
    call addTrace(errmsg,myname)
!
    select case(this%type_name_inversion_grid)
    !case('ccsInversionGrid')
    !   call getGeometryVtkCcsInversionGrid(this%ccs,points,cell_connectivity,cell_type,cell_indx_out,errmsg,cell_indx_req,indx_map)
    case('scartInversionGrid')
       call getGeometryVtkScartInversionGrid(this%scart,points,cell_connectivity,cell_type,cell_indx_out,&
            errmsg,cell_indx_req,indx_map)
    case('ecartInversionGrid')
       call getGeometryVtkEcartInversionGrid(this%ecart,points,cell_connectivity,cell_type,cell_indx_out,&
            errmsg,cell_indx_req,indx_map)
    case('specfem3dInversionGrid')
       call getGeometryVtkSpecfem3dInversionGrid(this%specfem3d,points,cell_connectivity,cell_type,&
            cell_indx_out,errmsg,cell_indx_req,indx_map)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select    
  end subroutine getGeometryVtkInversionGrid
!------------------------------------------------------------------------
!> \brief get for all cells the indices of their face neighbour cells
!! \param nb_idx pointer to array of length .ncell.this; if invgrid not defined yet, nullified on exit
!
  subroutine getIndicesFaceNeighboursInversionGrid(this,nb_idx)
    type (inversion_grid), intent(in) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
!
    nullify(nb_idx)
    select case(this%type_name_inversion_grid)
    !case('ccsInversionGrid')
    !   call getIndicesFaceNeighboursCcsInversionGrid(this%ccs,nb_idx)
    case('scartInversionGrid')
       call getIndicesFaceNeighboursScartInversionGrid(this%scart,nb_idx)
    case('ecartInversionGrid')
       call getIndicesFaceNeighboursEcartInversionGrid(this%ecart,nb_idx)
    case('specfem3dInversionGrid')
       call getIndicesFaceNeighboursSpecfem3dInversionGrid(this%specfem3d,nb_idx)
    end select
  end subroutine getIndicesFaceNeighboursInversionGrid
!------------------------------------------------------------------------
!> \brief get for all cells the wavefield point indices of points contained in that cell
!! \param this inversion grid
!! \param c1 vector of first coordinate of wavefield points
!! \param c2 vector of second coordinate of wavefield points
!! \param c3 vector of third coordinate of wavefield points
!! \param wp_idx pointer to array of length .ncell.this; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideInversionGrid(this,c1,c2,c3,wp_idx,errmsg)
    type (inversion_grid) :: this
    real, dimension(:), intent(in) :: c1,c2,c3
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
    !case('ccsInversionGrid')
    !   call locateWpInsideCcsInversionGrid(this%ccs,c1,c2,c3,wp_idx,errmsg)
    case('scartInversionGrid')
       call locateWpInsideScartInversionGrid(this%scart,c1,c2,c3,wp_idx,errmsg)
    case('ecartInversionGrid')
       call locateWpInsideEcartInversionGrid(this%ecart,c1,c2,c3,wp_idx,errmsg)
    case('specfem3dInversionGrid')
       call locateWpInsideSpecfem3dInversionGrid(this%specfem3d,c1,c2,c3,wp_idx,errmsg)
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
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param type_standard_cell defines the shape of the standard cell, select specific routine dependent on type (4=Tetrahedron,6=Hexahedron)
!! \param errmsg error message
!
  subroutine transformToStandardCellInversionGrid(this,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: c1,c2,c3,jacobian
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
    !case('ccsInversionGrid')
    !   call transformToStandardCellCcsInversionGrid(this%ccs,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('scartInversionGrid')
       call transformToStandardCellScartInversionGrid(this%scart,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('ecartInversionGrid')
       call transformToStandardCellEcartInversionGrid(this%ecart,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
    case('specfem3dInversionGrid')
       call transformToStandardCellSpecfem3dInversionGrid(this%specfem3d,icell,c1,c2,c3,jacobian,type_standard_cell,errmsg)
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
    !case('ccsInversionGrid')
    !   call getVolumeCellCcsInversionGrid(this%ccs,icell,volume,errmsg)
    case('scartInversionGrid')
       call getVolumeCellScartInversionGrid(this%scart,icell,volume,errmsg)
    case('ecartInversionGrid')
       call getVolumeCellEcartInversionGrid(this%ecart,icell,volume,errmsg)
    case('specfem3dInversionGrid')
       call getVolumeCellSpecfem3dInversionGrid(this%specfem3d,icell,volume,errmsg)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine getVolumeCellInversionGrid
!------------------------------------------------------------------------
!> \brief get center point of inversion grid cell
!! \details this point should be given in coordinates c1,c2,c3 as understood as
!!  wavefield point coordinates (compare routines locateWpInsideInversionGrid,
!!  transformToVtkInversionGrid etc.) It should represent some sort of point center
!!  of inversion grid cell icell, w.r.t. which quantities on cells (e.g. model
!!  values, kernel values etc.) will be assigned to in space, if a point 
!!  description of the respective quantity is required.
!! \param this inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param c1 first coordinate of center of cell icell
!! \param c2 second coordinate of center of cell icell
!! \param c3 third coordinate of center of cell icell
!! \param errmsg error message
!
  subroutine getCenterCellInversionGrid(this,icell,c1,c2,c3,errmsg)
    type (inversion_grid) :: this
    integer, intent(in) :: icell
    real :: c1,c2,c3
    type (error_message) :: errmsg
    ! local
    character (len=26) :: myname = 'getCenterCellInversionGrid'
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
    !case('ccsInversionGrid')
    !   call getCenterCellCcsInversionGrid(this%ccs,icell,c1,c2,c3,errmsg)
    case('scartInversionGrid')
       call getCenterCellScartInversionGrid(this%scart,icell,c1,c2,c3,errmsg)
    case('ecartInversionGrid')
       call getCenterCellEcartInversionGrid(this%ecart,icell,c1,c2,c3,errmsg)
    case('specfem3dInversionGrid')
       call getCenterCellSpecfem3dInversionGrid(this%specfem3d,icell,c1,c2,c3,errmsg)
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
    !case('ccsInversionGrid')
    !   call getRadiusCellCcsInversionGrid(this%ccs,icell,radius,errmsg)
    case('scartInversionGrid')
       call getRadiusCellScartInversionGrid(this%scart,icell,radius,errmsg)
    case('ecartInversionGrid')
       call getRadiusCellEcartInversionGrid(this%ecart,icell,radius,errmsg)
    case('specfem3dInversionGrid')
       call getRadiusCellSpecfem3dInversionGrid(this%specfem3d,icell,radius,errmsg)
    case default
       call add(errmsg,2,"inversion grid not yet defined",myname)
    end select
  end subroutine getRadiusCellInversionGrid
!
end module inversionGrid
