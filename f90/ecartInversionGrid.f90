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
!
!#########################################################################################
!#########################################################################################
!##    W  A  R  N  I  N  G 
!## 
!## SUBROUTINE createFaceNeighboursEcartInversionGrid DOES NOT YET WORK CORRECTLY, EVEN 
!## FOR TET4 CELLS! THERE WERE SOME TEST SETUPS, WHERE THE NEIGHBOUR SEARCH WORKED OUT 
!## FINE, FOR OTHER SETUPS IT DID NOT WORK OUT (in general, too few neighbours were found,
!## so maybe a matter of threashold epsilon stuff in the decision processes?!)
!#########################################################################################
!#########################################################################################
!
!> \brief external Cartesian inversion grid
!!
!! \details An inverison grid of type ecartInversionGrid consists of a collection of an arbitrary number of
!!  tetrahedral (tet4) and hexahedral (hex8) cells of arbitrarily distorted shapes. 
!!  The cells are defined via textfiles containing gridpoints in Cartesian x,y,z coordinates, and files defining
!!  the cells by refering to a collection of point indices. These files may be produced by any external meshing software
!!  (such as Cubit, confer module cubit2ASKIecartInversionGrid.py) as long as they meet the required format.
!!   the nodes coordinates text file(s) must be of following format:
!!      + the first line contains a single integer value, indicating the number of lines to 
!!        come (i.e. the number of points)
!!      + each following line contains 3 floating point numbers (separated by white space) 
!!        defining Cartesian X Y Z coordinates of a point
!!   the cell connectivity text files must be of following format:
!!     + the first line contains a single integer value, indicating the number of lines to 
!!       come (i.e. the number of cells)
!!     + each following line contains n integer numbers (separated by white space, n = 4 in case of 
!!       tet4 cells, n = 8 in case of hex8 cells), which define the corners (or higher order mid-edge 
!!       points etc.) of the cell and correspond to the point indices in the respective nodes coordinates 
!!       file, where the lowest point index is 1, corresponding to the second line (first point) in 
!!       the nodes coordinates file.
!!   THE ORDER OF THE POINTS IS ASSUMED TO CORRESPOND TO THE VTK CELL CONVENTIONS!
!!   in case one of the cell connectivity files does not exist, or their first line containing value 0, no cells
!!   of the respective type will be created
!!
!! \author Florian Schumacher
!! \date Jul 2013
!
module ecartInversionGrid
!
  use inputParameter
  use vectorPointer
  use realloc
  use errorMessage
!
  implicit none
!
  private :: createPointCellEcartInversionGrid,createFaceNeighboursEcartInversionGrid,&
       writeFaceNeighbourEcartInversionGrid,readFaceNeighbourEcartInversionGrid,&
       locateWpInsideTet4CellEcartInversionGrid,vectorProductEcartInversionGrid,&
       transformToStandardTetEcartInversionGrid,getVolumeTet4EcartInversionGrid,&
       locateWpInsideHex8CellEcartInversionGrid,transformToStandardHexEcartInversionGrid,&
       getVolumeHex8EcartInversionGrid
!
  type ecart_inversion_grid
     private
     logical :: is_defined = .false. !< flag indicating the correct definition (initialization) of the object (i.e. all following values)
!
     ! POINTS
     integer :: npoint !< overall number of base-points of the cells
     real, dimension(:,:), pointer :: point => null() !< (3,npoint)-array containing x,y,z coordinates of points
!
     ! CELLS
     integer :: ncell !< overall number of inversion grid cells
     integer, dimension(:), pointer :: ctype => null() !< for each cell, indicating its type (4 = tet4 , 8 = hex8,...)
     integer, dimension(:), pointer :: icell_type => null() !< for each cell, the local index in type-specific arrays cell_tet4,cell_hex8, ...
     type (integer_vector_pointer), dimension(:), pointer :: face_neighbour => null() !< for each cell, contains global cell indices of its neighbours
     ! tet4
     integer :: ncell_tet4 !< number of tet4 cells
     integer, dimension(:,:), pointer :: cell_tet4 => null() !< (4,ncell_tet4)-array (if any tet4 cells), containing cell connectivity / point indices
     integer, dimension(:), pointer :: iglob_tet4 => null() !< for each tet4 cell, contains its corresponding global cell index (index of array icell_type)
     ! hex8
     integer :: ncell_hex8 !< number of hex8 cells
     integer, dimension(:,:), pointer :: cell_hex8 => null() !< (8,ncell_hex8)-array (if any hex8 cells), containing cell connectivity / point indices
     integer, dimension(:), pointer :: iglob_hex8 => null() !< for each hex8 cell, contains its corresponding global cell index (index of array icell_type)
     ! further cell types possible here, just add
     !integer :: ncell_new_type
     !integer, dimension(:,:), pointer :: cell_new_type => null() !< (??,ncell_new_type)-array
!
     ! COORDINATES SPECIFICATION FOR VTK OUTPUT
     logical :: apply_vtk_coords_scaling_factor
     real :: vtk_coords_scaling_factor
     integer :: vtk_geometry_type_int = -1 !< type of vtk geometry:  0 = volumetric cells , 1 = cell center points
     ! in the future: there could be flags in parameter file like: DONT_SMOOTH_LAYER_BOUNDARIES, 
     ! or SMOOTHING_BOUNDARY_CONDITIONS which could be taken into account here, and memorized for better handling 
     ! of smoothing conditions in calls to certain routines below
  end type ecart_inversion_grid
!
contains
!------------------------------------------------------------------------
!> \brief logical return whether this ecart_inversion_grid is able to transform points (of given coords type) to vtk plot projection
!
  function canTransformToVtkPointsOutsideEcartInversionGrid(this,coords_type) result(l)
    type(ecart_inversion_grid) :: this
    character(len=*) :: coords_type
    logical :: l
    ! the ecart_inversion_grid has capability to transform any points to vtk,
    ! provided the object is defined
    l = this%is_defined
  end function canTransformToVtkPointsOutsideEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get unit factor of the volume element
!! \param this ecart inversion grid
!! \param uf_wp unit factor of wavefield points
!! \param uf_vol unit factor of volume element (return value of this subroutine)
!! \param errmsg error message
!
  subroutine getUnitFactorOfVolumeElementEcartInversionGrid(this,uf_wp,uf_vol,errmsg)
    type (ecart_inversion_grid) :: this
    real :: uf_wp,uf_vol
    type (error_message) :: errmsg
    character (len=46) :: myname = 'getUnitFactorOfVolumeElementEcartInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) call add(errmsg,1,"be aware that the inversion grid not yet defined; "//&
         "however, the unit factor of the volume element can be correctly computed at this point",myname)
!
    ! The ecart inversion grid does not assume specific units for its spatial extension but simply 
    ! assumes that they are the same as for the wavefield point coordinates (which are located inside the
    ! inversion grid cells only according to their pure numerical values). Hence, for this 3D volumetric
    ! inversion grid, the volume element has a unit which is the cube of the unit of the wavefield points.
    uf_vol = uf_wp*uf_wp*uf_wp
  end subroutine getUnitFactorOfVolumeElementEcartInversionGrid
!------------------------------------------------------------------------
!> \brief map vtk geometry type names to integers
!
  function intGeometryTypeEcartInversionGrid(vtk_geometry_type_str) result(vtk_geometry_type_int)
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
  end function intGeometryTypeEcartInversionGrid
!------------------------------------------------------------------------
!> \brief create external Cartesian inversion grid
!! \details the parameter file given, must contain all necessary parameters to define
!!  an object of this type. 
!! \param this external Cartesian inversion grid
!! \param parfile filename of parameter file containing definintion of this inversion grid
!! \param path path where to find/write own files (usually path of current iteration step)
!! \param lu file unit to use for reading and writing files
!! \param errmsg error message
!
  subroutine createEcartInversionGrid(this,parfile,path,lu,errmsg,recreate)
    type (ecart_inversion_grid) :: this
    character(len=*) :: parfile,path
    integer :: lu
    type (error_message) :: errmsg
    logical, optional :: recreate
    ! local
    character(len=24) :: myname = 'createEcartInversionGrid'
    logical :: recreate_neighbours,file_exists,is_binary
    integer :: ios
    ! parfile
    type (input_parameter) :: inpar
    character (len=80), dimension(11) :: inpar_keys
    data inpar_keys/'ECART_INVGRID_FILE_NODES_HEX8', 'SCALE_VTK_COORDS', 'ECART_INVGRID_FILE_NODES_TET4', &
         'ECART_INVGRID_USE_NODES_COMMON', 'ECART_INVGRID_FILE_CELLS_TET4', &
         'ECART_INVGRID_FILE_NEIGHBOURS_IS_BINARY', 'ECART_INVGRID_FILE_NEIGHBOURS', 'ECART_INVGRID_FILE_NODES_COMMON', &
         'VTK_COORDS_SCALING_FACTOR', 'VTK_GEOMETRY_TYPE', 'ECART_INVGRID_FILE_CELLS_HEX8'/
!
    call addTrace(errmsg,myname)
    if(this%is_defined) then
       call add(errmsg,1,"this object is already defined, deallocating it now before creating new one",myname)
       call deallocateEcartInversionGrid(this)
    end if
!
    call createKeywordsInputParameter(inpar,inpar_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) return
!
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
    this%vtk_geometry_type_int = intGeometryTypeEcartInversionGrid(inpar.sval.'VTK_GEOMETRY_TYPE')
    if(this%vtk_geometry_type_int < 0) then
       call add(errmsg,2,"vtk geometry type '"//trim(inpar.sval.'VTK_GEOMETRY_TYPE')//"' is invalid",myname)
       goto 1
    end if
!
    is_binary = lval(inpar,'ECART_INVGRID_FILE_NEIGHBOURS_IS_BINARY',iostat=ios)
    if(ios /= 0) then
       call add(errmsg,2,"could not read logical value for 'ECART_INVGRID_FILE_NEIGHBOURS_IS_BINARY' from '"//&
            trim(inpar.sval.'ECART_INVGRID_FILE_NEIGHBOURS_IS_BINARY')//"'",myname)
       goto 1
    end if
!
    ! create points and cell connectivity arrays
    call createPointCellEcartInversionGrid(this,inpar,path,lu,errmsg)
    if(.level.errmsg == 2) goto 1
!
    ! if recreate is present and true or file ECART_INVGRID_FILE_NEIGHBOURS does not exist, create neighbours
    ! anew and write to file ECART_INVGRID_FILE_NEIGHBOURS, otherwise just read in neighbours from file 
    ! ECART_INVGRID_FILE_NEIGHBOURS
!
    recreate_neighbours = .false.
    if(present(recreate)) then
       if(recreate) recreate_neighbours = .true.
    end if
!
    inquire(file = trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NEIGHBOURS'), exist = file_exists)
!
    if(recreate_neighbours .or. .not.file_exists) then
       call createFaceNeighboursEcartInversionGrid(this,errmsg)
       if(.level.errmsg == 2) goto 1
       call writeFaceNeighbourEcartInversionGrid(this,&
            trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NEIGHBOURS'),is_binary,file_exists,lu,errmsg)
    else
       call readFaceNeighbourEcartInversionGrid(this,trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NEIGHBOURS'),&
            is_binary,lu,errmsg)
    end if
!
    ! if everything was ok, indicate so and return, otherwise destroy whatever was created so far
1   this%is_defined = .level.errmsg /= 2
    if(.not.this%is_defined) call deallocateEcartInversionGrid(this)
!
  end subroutine createEcartInversionGrid
!------------------------------------------------------------------------
!> \brief create node and cell connectivity for external Cartesian inversion grid from nodes and cell connectivity files
!! \param this external Cartesian inversion grid
!! \param inpar input parameter object containing ecart parameter file
!! \param path path where to find/write own files (usually path of current iteration step)
!! \param lu file unit to read ponts and cell connectivity files
!! \param errmsg error message
!
  subroutine createPointCellEcartInversionGrid(this,inpar,path,lu,errmsg)
    type (ecart_inversion_grid) :: this
    type (input_parameter) :: inpar
    character(len=*) :: path
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=33) :: myname = 'createPointCellEcartInversionGrid'
    integer :: ios,ipoint,npoint_tet4,npoint_hex8,icell,point_shift,ncell
    logical :: read_tet4,read_hex8,use_nodes_common,file_exists
    character(len=1) :: dummy_char
!
    call addTrace(errmsg,myname)
!
    use_nodes_common = lval(inpar,'ECART_INVGRID_USE_NODES_COMMON',iostat=ios)
    if(ios/=0) then
       call add(errmsg,2,"could not read logical value for 'ECART_INVGRID_USE_NODES_COMMON' from '"//&
            trim(inpar.sval.'ECART_INVGRID_USE_NODES_COMMON')//"'",myname)
       return
    end if
!
! FIRST CHECK WHICH CELL CONNECTIVITY FILES EXIST AND CONTAIN A POSITIVE NUMBER OF CELLS
!
    ! read tet4 files?
    read_tet4 = .false.; this%ncell_tet4 = 0
    inquire(file = trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4'), exist = file_exists)
    if(file_exists) then
       ! open tet4 cell file
       open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4'),form='formatted',&
            status='old',action='read',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open tet4 cell file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"' to read",myname)
          return
       end if
       ! read number of cells from first line
       read(lu,*,iostat=ios) this%ncell_tet4
       close(lu)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not read number of tet4 cells from tet4 cell file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"'",myname)
          return
       end if
       if(this%ncell_tet4 <= 0) then
          this%ncell_tet4 = 0
       else
          read_tet4 = .true.
       end if
    else ! file_exists
       call add(errmsg,1,"tet4 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"' "//&
            "does not exist, hence assume 0 tet4 cells for this inversion grid",myname)
    end if ! file_exists
!
    ! read hex8 files?
    read_hex8 = .false.; this%ncell_hex8 = 0
    inquire(file = trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8'), exist = file_exists)
    if(file_exists) then
       ! open hex8 cell file
       open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8'),form='formatted',&
            status='old',action='read',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open hex8 cell file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"' to read",myname)
          return
       end if
       ! read number of cells from first line
       read(lu,*,iostat=ios) this%ncell_hex8
       close(lu)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not read number of hex8 cells from hex8 cell file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"'",myname)
          return
       end if
       if(this%ncell_hex8 <= 0) then
          this%ncell_hex8 = 0
       else
          read_hex8 = .true.
       end if
    else ! file_exists
       call add(errmsg,1,"hex8 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"' "//&
            "does not exist, hence assume 0 hex8 cells for this inversion grid",myname)
    end if ! file_exists
!
    if(.not.read_tet4 .and. .not.read_hex8) then
       call add(errmsg,2,"there are no cells to create",myname)
       return
    else
       this%ncell = this%ncell_tet4 + this%ncell_hex8
    end if
!
! READ IN NODES COORDINATES, EITHER FROM COMMON FILE, OR FROM INDIVIDUAL FILES (ONLY FOR CELL TYPES FOR WHICH read_TYPE == .true.)
!
    if(use_nodes_common) then
       ! open common nodes file
       open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_COMMON'),form='formatted',&
            status='old',action='read',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open common nodes file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_NODES_COMMON')//"' to read",myname)
          return
       end if
       ! read number of points from first line
       read(lu,*,iostat=ios) this%npoint
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not read number of nodes from common nodes file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_NODES_COMMON')//"'",myname)
          return
       end if
       if(this%npoint <= 0) then
          close(lu)
          write(errstr,*) "integer ",this%npoint," on first line in common nodes file '"//trim(path)//&
               trim(inpar.sval.'ECART_INVGRID_FILE_NODES_COMMON')//"' must be positive (inidcates number of lines to come)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       ! allocate array this%point
       allocate(this%point(3,this%npoint))
       ! iterate over the next npoint lines and read in coordinates
       do ipoint = 1,this%npoint
          read(lu,*,iostat=ios) this%point(:,ipoint)
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read Cartesian X Y Z coordinates for ",ipoint,"'th point from line ",ipoint+1,&
                  " of common nodes file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_COMMON')//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
       end do ! ipoint
       close(lu)
!
    else ! use_nodes_common
       this%npoint = 0
!
       ! read in tet4 nodes file, if read_tet4
!
       if(read_tet4) then
          ! open nodes file and get number of points
          open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_TET4'),form='formatted',&
               status='old',action='read',iostat=ios)
          if(ios/=0) then
             close(lu)
             call add(errmsg,2,"could not open tet4 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_TET4')//"' to read",myname)
             return
          end if
          ! read number of points from first line
          read(lu,*,iostat=ios) npoint_tet4
          if(ios/=0) then
             close(lu)
             call add(errmsg,2,"could not read number of nodes from tet4 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_TET4')//"'",myname)
             return
          end if
          if(npoint_tet4 <= 0) then
             close(lu)
             write(errstr,*) "integer ",npoint_tet4," on first line in tet4 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_TET4')//"' must be positive (inidcates number of lines to come)"
             call add(errmsg,2,errstr,myname)
             return
          end if
          ! reallocate this%point (will be allocated, if not associated)
          this%point => reallocate(this%point,3,this%npoint+npoint_tet4)
          ! iterate over the next npoint_tet4 lines and read in coordinates
          do ipoint = 1,npoint_tet4
             read(lu,*,iostat=ios) this%point(:,this%npoint+ipoint)
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read Cartesian X Y Z coordinates for ",ipoint,"'th point from line ",ipoint+1,&
                     " of tet4 nodes file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_TET4')//"'"
                call add(errmsg,2,errstr,myname)
                return
             end if
          end do ! ipoint
          close(lu)
          ! shift counter this%npoint
          this%npoint = this%npoint + npoint_tet4
       end if ! read_tet4
!
       ! read in hex8 nodes file, if read_hex8
!
       if(read_hex8) then
          ! open nodes file and get number of points
          open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_HEX8'),form='formatted',&
               status='old',action='read',iostat=ios)
          if(ios/=0) then
             close(lu)
             call add(errmsg,2,"could not open hex8 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_HEX8')//"' to read",myname)
             return
          end if
          ! read number of points from first line
          read(lu,*,iostat=ios) npoint_hex8
          if(ios/=0) then
             close(lu)
             call add(errmsg,2,"could not read number of nodes from hex8 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_HEX8')//"'",myname)
             return
          end if
          if(npoint_hex8 <= 0) then
             close(lu)
             write(errstr,*) "integer ",npoint_hex8," on first line in hex8 nodes file '"//trim(path)//&
                  trim(inpar.sval.'ECART_INVGRID_FILE_NODES_HEX8')//"' must be positive (inidcates number of lines to come)"
             call add(errmsg,2,errstr,myname)
             return
          end if
          ! reallocate this%point (will be allocated, if not associated)
          this%point => reallocate(this%point,3,this%npoint+npoint_hex8)
          ! iterate over the next npoint_hex8 lines and read in coordinates
          do ipoint = 1,npoint_hex8
             read(lu,*,iostat=ios) this%point(:,this%npoint+ipoint)
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read Cartesian X Y Z coordinates for ",ipoint,"'th point from line ",ipoint+1,&
                     " of hex8 nodes file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_NODES_HEX8')//"'"
                call add(errmsg,2,errstr,myname)
                return
             end if
          end do ! ipoint
          close(lu)
          this%npoint = this%npoint + npoint_hex8
       end if ! read_hex8
!
    end if ! use_nodes_common
!
! READ IN CELL CONNECTIVITY FILES FOR CELL TYPES FOR WHICH NCELL_TYPE > 0
!
    allocate(this%ctype(this%ncell),this%icell_type(this%ncell))
    ! cell counter for arrays ctype,icell_type
    ncell = 0
    ! in case of .not.use_nodes_common, keep track of range in array this%point in order to properly shift point indices contained in cell files
    point_shift = 0
!
    ! tet4 cell file
    if(read_tet4) then
       allocate(this%cell_tet4(4,this%ncell_tet4),this%iglob_tet4(this%ncell_tet4))
       this%ctype(ncell+1:ncell+this%ncell_tet4) = 4
       ! open tet4 cell file
       open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4'),form='formatted',&
            status='old',action='read')
       ! skip first line (already read above)
       read(lu,"(a1)") dummy_char
       ! iterate over the next this%ncell_tet4 lines and read in array cell_tet4
       do icell = 1,this%ncell_tet4
          read(lu,*,iostat=ios) this%cell_tet4(:,icell)
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read 4 node indices for ",icell,"'th tet4 cell from line ",icell+1,&
                  " of tet4 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
          if(use_nodes_common) then
             if(any(this%cell_tet4(:,icell)<1 .or. this%cell_tet4(:,icell)>this%npoint)) then
                close(lu)
                write(errstr,*) "invalid node indices for ",icell,"'th tet4 cell from line ",icell+1,&
                  " of tet4 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"'"
                call add(errmsg,2,errstr,myname)
             end if
          else ! use_nodes_common
             ! if use individual node files, the cell files contain node indices w.r.t. the individual node files
             if(any(this%cell_tet4(:,icell)<1 .or. this%cell_tet4(:,icell)>npoint_tet4)) then
                close(lu)
                write(errstr,*) "invalid node indices for ",icell,"'th tet4 cell from line ",icell+1,&
                  " of tet4 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_TET4')//"'"
                call add(errmsg,2,errstr,myname)
             end if
             ! shift the node indices by the starting index in global node array this%point
             ! (which at this point in the code is the value of point_shift)
             this%cell_tet4(:,icell) = this%cell_tet4(:,icell) + point_shift
          end if ! use_nodes_common
          this%icell_type(ncell+icell) = icell
          this%iglob_tet4(icell) = ncell+icell
       end do ! icell
       close(lu)
       ! increase shift in global node array this%point
       point_shift = point_shift + npoint_tet4
       ! increase global cell counter for arrays ctype,icell_type
       ncell = ncell + this%ncell_tet4
    end if ! read_tet4
!
    ! hex8 cell file
    if(read_hex8) then
       allocate(this%cell_hex8(8,this%ncell_hex8),this%iglob_hex8(this%ncell_hex8))
       this%ctype(ncell+1:ncell+this%ncell_hex8) = 8
       ! open hex8 cell file
       open(unit=lu,file=trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8'),form='formatted',&
            status='old',action='read')
       ! skip first line (already read above)
       read(lu,"(a1)") dummy_char
       ! iterate over the next this%ncell_hex8 lines and read in array cell_hex8
       do icell = 1,this%ncell_hex8
          read(lu,*,iostat=ios) this%cell_hex8(:,icell)
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read 8 node indices for ",icell,"'th hex8 cell from line ",icell+1,&
                  " of hex8 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
          if(use_nodes_common) then
             if(any(this%cell_hex8(:,icell)<1 .or. this%cell_hex8(:,icell)>this%npoint)) then
                close(lu)
                write(errstr,*) "invalid node indices for ",icell,"'th hex8 cell on line ",icell+1,&
                  " of hex8 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"'"
                call add(errmsg,2,errstr,myname)
             end if
          else ! use_nodes_common
             ! if use individual node files, the cell files contain node indices w.r.t. the individual node files
             if(any(this%cell_hex8(:,icell)<1 .or. this%cell_hex8(:,icell)>npoint_hex8)) then
                close(lu)
                write(errstr,*) "invalid node indices for ",icell,"'th hex8 cell on line ",icell+1,&
                  " of hex8 cell file '"//trim(path)//trim(inpar.sval.'ECART_INVGRID_FILE_CELLS_HEX8')//"'"
                call add(errmsg,2,errstr,myname)
             end if
             ! shift the node indices by the starting index in global node array this%point
             ! (which at this point in the code is the value of point_shift)
             this%cell_hex8(:,icell) = this%cell_hex8(:,icell) + point_shift
          end if
          this%icell_type(ncell+icell) = icell
          this%iglob_hex8(icell) = ncell+icell
       end do ! icell
       close(lu)
       ! increase shift in global node array this%point
       point_shift = point_shift + npoint_hex8
       ! increase global cell counter for arrays ctype,icell_type
       ncell = ncell + this%ncell_tet4
    end if ! read_hex8
!
  end subroutine createPointCellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief create neighbours this%face_neighbour of the inversion grid cells
!
  subroutine createFaceNeighboursEcartInversionGrid(this,errmsg)
    type (ecart_inversion_grid) :: this
    type (error_message) :: errmsg
    ! local
    character(len=38) :: myname = 'createFaceNeighboursEcartInversionGrid'
    character(len=400) :: errstr
    integer :: icell,icell_tet4,jcell_tet4,iface,inb,nnb_potential,nnb_add,nnb,j
    real :: eps,origin_distance_tet4_test
    real, dimension(:,:), allocatable :: outward_normal_tet4
    real, dimension(:), allocatable :: dot,origin_distance_tet4
    logical, dimension(:), allocatable :: is_potential_neighbour_face
    real, dimension(3) :: p1,p2,p3,p4
    integer, dimension(:), pointer :: idx
    integer, dimension(:), allocatable :: idx_add
!!$    ! LAPACK SGESVD
!!$    real, dimension(5,6) :: A_tet4
!!$    real, dimension(5,7) :: Ab_tet4
!!$    integer :: LWORK_SGESVD_A_tet4,LWORK_SGESVD_Ab_tet4,INFO
!!$    integer :: rank_A_tet4,rank_Ab_tet4
!!$    real, dimension(:),allocatable :: WORK_SGESVD_A_tet4,WORK_SGESVD_Ab_tet4
!!$    real, dimension(5) :: S_A_tet4,S_Ab_tet4
!!$    real, dimension(1,1) :: U,VT
    ! Triangle Triangle Test TTT
    real, dimension(3) :: p1_TTT,q1_TTT,r1_TTT,p2_TTT,q2_TTT,r2_TTT,n_TTT
!
    nullify(idx)
!
! BEWARE: SO FAR, THIS ROUTINE ONLY FINDS TET4 NEIGHBOURS OF TET4 CELLS! HEX8 CELLS ARE COMPLETELY IGNORED FOR NOW (complicated
    if(this%ncell_hex8>0) &
         call add(errmsg,1,"there are 8-node hexahedral cells in this inversion grid: "//&
         "for this cell type, the neighbour search is not yet supported and those "//&
         "cells will simply have no neighbours and will not be a neighbour of any other cell",myname)
!
!
    ! as this routine is private, it should have been assured that this inversion grid was successfully created
    ! before calling this routine, so need no checks here
    allocate(this%face_neighbour(this%ncell))
!
    eps = epsilon(1.0)
!!$print *, "eps = ",eps
!
    ! for all tet4 cells check if the nodes are correctly oriented and the tet is not degenerate
    ! and compute all outward normal vectors of all 4 faces
    allocate(outward_normal_tet4(4*this%ncell_tet4,3),origin_distance_tet4(4*this%ncell_tet4))
    iface = 0
    do icell_tet4 = 1,this%ncell_tet4
       ! get all nodes of this thetrahedron
       p1 = this%point(:,this%cell_tet4(1,icell_tet4))
       p2 = this%point(:,this%cell_tet4(2,icell_tet4))
       p3 = this%point(:,this%cell_tet4(3,icell_tet4))
       p4 = this%point(:,this%cell_tet4(4,icell_tet4))
!
       ! check correct orientation of nodes
       if(checkOrientationTet4EcartInversionGrid(p1,p2,p3,p4) < 3*eps) then
          if(checkOrientationTet4EcartInversionGrid(p1,p2,p3,p4) <= -3*eps) then
             write(errstr,*) this%iglob_tet4(icell_tet4),&
                  "'th inversion grid cell (4-node tetrahedral) has wrong node indexing: ",&
                  checkOrientationTet4EcartInversionGrid(p1,p2,p3,p4)
          else
             write(errstr,*) this%iglob_tet4(icell_tet4),&
                  "'th inversion grid cell (4-node tetrahedral) is degenerate (points are coplanar): ",&
                  checkOrientationTet4EcartInversionGrid(p1,p2,p3,p4)
          end if
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       ! first face (opposite node p1, outward normal is normalized vector product of vectors p3-p2 and p4-p2)
       iface = iface + 1
       outward_normal_tet4(iface,:) = vectorProductEcartInversionGrid(p3-p2,p4-p2,normal=.true.)
       origin_distance_tet4(iface) = sum(outward_normal_tet4(iface,:)*p2)
!
       ! second face (opposite node p2, outward normal is normalized vector product of vectors p4-p1 and p3-p1)
       iface = iface + 1
       outward_normal_tet4(iface,:) = vectorProductEcartInversionGrid(p4-p1,p3-p1,normal=.true.)
       origin_distance_tet4(iface) = sum(outward_normal_tet4(iface,:)*p1)
!
       ! third face (opposite node p3, outward normal is normalized vector product of vectors p2-p1 and p4-p1)
       iface = iface + 1
       outward_normal_tet4(iface,:) = vectorProductEcartInversionGrid(p2-p1,p4-p1,normal=.true.)
       origin_distance_tet4(iface) = sum(outward_normal_tet4(iface,:)*p1)
!
       ! fourth face (opposite node p4, outward normal is normalized vector product of vectors p3-p1 and p2-p1)
       iface = iface + 1
       outward_normal_tet4(iface,:) = vectorProductEcartInversionGrid(p3-p1,p2-p1,normal=.true.)
       origin_distance_tet4(iface) = sum(outward_normal_tet4(iface,:)*p1)
    end do ! icell_tet4
!
! prepare values for search below
!
!##################################################################################
!##################################################################################
! THE FOLLOWING DOES  N O T  WORK!!, (see below
!!$    ! conduct a workspace query for SGESVD operations on A_tet4
!!$    A_tet4 = 0.
!!$    allocate(WORK_SGESVD_A_tet4(1)); LWORK_SGESVD_A_tet4 = -1
!!$    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
!!$    call SGESVD( 'N', 'N', 5 , 6, A_tet4, 5, S_A_tet4, U, 1, VT, 1, WORK_SGESVD_A_tet4, LWORK_SGESVD_A_tet4, INFO )
!!$    if(INFO/=0) then
!!$       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to A returned INFO = ',INFO
!!$       call add(errmsg,2,errstr,myname)
!!$       goto 1
!!$    endif ! INFO/=0
!!$    LWORK_SGESVD_A_tet4 = WORK_SGESVD_A_tet4(1)
!!$    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A_tet4 is ',LWORK_SGESVD_A_tet4
!!$    call add(errmsg,0,errstr,myname)
!!$    deallocate(WORK_SGESVD_A_tet4); allocate(WORK_SGESVD_A_tet4(LWORK_SGESVD_A_tet4))
!!$!
!!$    ! conduct a workspace query for SGESVD operations on Ab_tet4
!!$    Ab_tet4 = 0.
!!$    allocate(WORK_SGESVD_Ab_tet4(1)); LWORK_SGESVD_Ab_tet4 = -1
!!$    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
!!$    call SGESVD( 'N', 'N', 5 , 7, Ab_tet4, 5, S_A_tet4, U, 1, VT, 1, WORK_SGESVD_Ab_tet4, LWORK_SGESVD_Ab_tet4, INFO )
!!$    if(INFO/=0) then
!!$       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to Ab_tet4 returned INFO = ',INFO
!!$       call add(errmsg,2,errstr,myname)
!!$       goto 1
!!$    endif ! INFO/=0
!!$    LWORK_SGESVD_Ab_tet4 = WORK_SGESVD_Ab_tet4(1)
!!$    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to Ab_tet4 is ',LWORK_SGESVD_A_tet4
!!$    call add(errmsg,0,errstr,myname)
!!$    deallocate(WORK_SGESVD_Ab_tet4); allocate(WORK_SGESVD_Ab_tet4(LWORK_SGESVD_Ab_tet4))
!##################################################################################
!##################################################################################
!
    allocate(dot(4*this%ncell_tet4),is_potential_neighbour_face(4*this%ncell_tet4))
!
    ! now loop on all cells and chech each face
    do icell = 1,this%ncell
       select case(this%ctype(icell))
       case(4)
          ! initiate index vector of neighbours idx for this cell
          nnb = 0; nullify(idx)
          !
          icell_tet4 = this%icell_type(icell)
          ! get all nodes of this thetrahedron
          p1 = this%point(:,this%cell_tet4(1,icell_tet4))
          p2 = this%point(:,this%cell_tet4(2,icell_tet4))
          p3 = this%point(:,this%cell_tet4(3,icell_tet4))
          p4 = this%point(:,this%cell_tet4(4,icell_tet4))
          do iface = 1,4
             ! compute the origin distance of the current face plane (is tested below against origin distances
             ! of potential neighbour face planes)
             if(iface==1) then
                origin_distance_tet4_test = sum(outward_normal_tet4((icell_tet4-1)*4+iface,:)*p2)
             else
                origin_distance_tet4_test = sum(outward_normal_tet4((icell_tet4-1)*4+iface,:)*p1)
             end if
!
             ! compute the scalar product of the outward normal of this face with the outward normals of ALL faces, 
             ! only those faces, where this scalar product is -1, and where the origin distance of the face plane
             ! is the same up to sign (in fact must have different sign, as the normals differ by sign) are potential neighbour faces
!
             ! SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
             call SGEMV('N',4*this%ncell_tet4,3,1.,outward_normal_tet4,4*this%ncell_tet4,&
                  outward_normal_tet4((icell_tet4-1)*4+iface,:),1,0.,dot,1)
!
!!$if(icell==3) then
!!$print *, "cell 1, face = ",iface,", dot product with cell 6 faces:"
!!$print *, "   cell 6 face 1: ",dot(21)
!!$print *, "   cell 6 face 2: ",dot(22)
!!$print *, "   cell 6 face 3: ",dot(23)
!!$print *, "   cell 6 face 4: ",dot(24)
!!$print *, "cell 3, face = ",iface,", dot product with cell 5 faces:"
!!$print *, "   cell 5 face 1: ",dot(17)
!!$print *, "   cell 5 face 2: ",dot(18)
!!$print *, "   cell 5 face 3: ",dot(19)
!!$print *, "   cell 5 face 4: ",dot(20)
!!$print *, "cell 3, face = ",iface,", dot product with cell 12 faces:"
!!$print *, "   cell 12 face 1: ",dot(45)
!!$print *, "   cell 12 face 2: ",dot(46)
!!$print *, "   cell 12 face 3: ",dot(47)
!!$print *, "   cell 12 face 4: ",dot(48)
!!$end if
!
             ! is_potential_neighbour_face = dot == - 1. .and. origin_distance_tet4_test == -origin_distance_tet4((icell_tet4-1)*4+iface)
             ! use 3* machine-epsilon as treshold here, as dot product included 3 multiplications
             is_potential_neighbour_face = abs(dot + 1.) < 3*eps .and. &
                  abs(origin_distance_tet4 + origin_distance_tet4_test) < 3*eps
             nnb_potential = count(is_potential_neighbour_face)
!print *, "cell ",icell_tet4,", face ",iface,", potential neighbours ",nnb_potential
             if(nnb_potential == 0) cycle
!
!
!##################################################################################
!##################################################################################
! THE FOLLOWING DOES  N O T  WORK!!,
! as the set of points under consideration is actually NOT the 2D-convex hulls
! since it cannot be assured that the solution vectors contain strictly positive
! values. For triangles in the same plane, there can be a neighbour detection
! even if the triangles do not overlap (in this case the linear system Ax=b might
! have solutions, but negative ones!)
!##################################################################################
!##################################################################################
!
!!$             ! now select from potential neighbours the real neighbours by removing those, which are no neighbours.
!!$             ! do this by investigating, if the linear system Ax=b has a solution, where
!!$             ! 
!!$             !     ( p1 | p2 | p3 | -q1 | -q2 | -q3 )
!!$             ! A = ( 1. | 1. | 1. |  0. |  0. |  0. ) , b = ( 0. , 0., 0., 1., 1.) , 
!!$             !     ( 0. | 0. | 0. |  1. |  1. |  1. )
!!$             !
!!$             ! p1,p2,p3 are column vectors of corners of this current face, q1,q2,q3 are column vectors of corners of potential neighbour face
!!$             !
!!$             ! IF SYSTEM Ax=b HAS A SOLUTION, IT MEANS THAT THE 2D-CONVEX HULLS CREATED BY THE CORNERS OF BOTH FACES (i.e. the face surfaces)
!!$             ! HAVE COMMON POINTS, i.e. OVERLAP, HENCE THE TWO FACES ARE NEIGHBOURING FACES
!!$             !
!!$             ! use the following criterion: Ax=b has a solution if and oly if rank(A)==rank(A,b), where A,b is the matrix A with additional column b
!!$!
!!$             allocate(idx_add(nnb_potential))
!!$             idx_add = pack( (/ (j,j=1,4*this%ncell_tet4) /) , is_potential_neighbour_face )
!!$             nnb_add = nnb_potential
!!$             do inb = 1,nnb_potential
!!$                ! initiate values for A_tet4, Ab_tet4 which are independent of any points
!!$                A_tet4(4,:) = (/ 1., 1., 1., 0., 0., 0. /)
!!$                A_tet4(5,:) = (/ 0., 0., 0., 1., 1., 1. /)
!!$                Ab_tet4(4,1:6) = (/ 1., 1., 1., 0., 0., 0. /) ! 4th,5th row of A
!!$                Ab_tet4(5,1:6) = (/ 0., 0., 0., 1., 1., 1. /)
!!$                Ab_tet4(:,7) = (/ 0., 0., 0., 1., 1. /) ! rhs-vector b in column 7 of Ab_tet4
!!$                ! put corner points of the current face into matrices A and Ab, rows 1,2,3 / columns 1,2,3
!!$                select case(iface)
!!$                case(1)
!!$                   A_tet4(1:3,1) = p2; Ab_tet4(1:3,1) = p2
!!$                   A_tet4(1:3,2) = p4; Ab_tet4(1:3,2) = p4
!!$                   A_tet4(1:3,3) = p3; Ab_tet4(1:3,3) = p3
!!$                case(2)
!!$                   A_tet4(1:3,1) = p1; Ab_tet4(1:3,1) = p1
!!$                   A_tet4(1:3,2) = p3; Ab_tet4(1:3,2) = p3
!!$                   A_tet4(1:3,3) = p4; Ab_tet4(1:3,3) = p4
!!$                case(3)
!!$                   A_tet4(1:3,1) = p1; Ab_tet4(1:3,1) = p1
!!$                   A_tet4(1:3,2) = p4; Ab_tet4(1:3,2) = p4
!!$                   A_tet4(1:3,3) = p2; Ab_tet4(1:3,3) = p2
!!$                case(4)
!!$                   A_tet4(1:3,1) = p1; Ab_tet4(1:3,1) = p1
!!$                   A_tet4(1:3,2) = p2; Ab_tet4(1:3,2) = p2
!!$                   A_tet4(1:3,3) = p3; Ab_tet4(1:3,3) = p3
!!$                end select ! case iface
!!$                ! put corner points of this potential neighbour face into matrices A and Ab, rows 1,2,3 / columns 4,5,6
!!$                select case(mod(idx_add(inb),4))
!!$                case(1) ! means face 1
!!$                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
!!$                   A_tet4(1:3,4) = - this%point(:,this%cell_tet4(2,jcell_tet4))
!!$                   A_tet4(1:3,5) = - this%point(:,this%cell_tet4(4,jcell_tet4))
!!$                   A_tet4(1:3,6) = - this%point(:,this%cell_tet4(3,jcell_tet4))
!!$                   Ab_tet4(1:3,4) = A_tet4(1:3,4)
!!$                   Ab_tet4(1:3,5) = A_tet4(1:3,5)
!!$                   Ab_tet4(1:3,6) = A_tet4(1:3,6)
!!$                case(2) ! means face 2
!!$                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
!!$                   A_tet4(1:3,4) = - this%point(:,this%cell_tet4(1,jcell_tet4))
!!$                   A_tet4(1:3,5) = - this%point(:,this%cell_tet4(3,jcell_tet4))
!!$                   A_tet4(1:3,6) = - this%point(:,this%cell_tet4(4,jcell_tet4))
!!$                   Ab_tet4(1:3,4) = A_tet4(1:3,4)
!!$                   Ab_tet4(1:3,5) = A_tet4(1:3,5)
!!$                   Ab_tet4(1:3,6) = A_tet4(1:3,6)
!!$                case(3) ! means face 3
!!$                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
!!$                   A_tet4(1:3,4) = - this%point(:,this%cell_tet4(1,jcell_tet4))
!!$                   A_tet4(1:3,5) = - this%point(:,this%cell_tet4(4,jcell_tet4))
!!$                   A_tet4(1:3,6) = - this%point(:,this%cell_tet4(2,jcell_tet4))
!!$                   Ab_tet4(1:3,4) = A_tet4(1:3,4)
!!$                   Ab_tet4(1:3,5) = A_tet4(1:3,5)
!!$                   Ab_tet4(1:3,6) = A_tet4(1:3,6)
!!$                case(0) ! means face 4, as in that case the index on all faces is a multiple of 4
!!$                   jcell_tet4 = idx_add(inb)/4 ! tet4 cell index of cell, to which this face belongs
!!$                   A_tet4(1:3,4) = - this%point(:,this%cell_tet4(1,jcell_tet4))
!!$                   A_tet4(1:3,5) = - this%point(:,this%cell_tet4(2,jcell_tet4))
!!$                   A_tet4(1:3,6) = - this%point(:,this%cell_tet4(3,jcell_tet4))
!!$                   Ab_tet4(1:3,4) = A_tet4(1:3,4)
!!$                   Ab_tet4(1:3,5) = A_tet4(1:3,5)
!!$                   Ab_tet4(1:3,6) = A_tet4(1:3,6)
!!$                end select
!!$!
!!$                ! compute rank of A_tet4
!!$                call SGESVD( 'N', 'N', 5 , 6, A_tet4, 5, S_A_tet4, U, 1, VT, 1, WORK_SGESVD_A_tet4, LWORK_SGESVD_A_tet4, INFO )
!!$                if(INFO/=0) then
!!$                   j = mod(idx_add(inb),4); if(j==0) j=4
!!$                   write(errstr,*) "computation of singular values of A_tet4 failed comparing ",iface,"'th face of ",&
!!$                        icell,"'th cell and ",j,"'th face of ",this%iglob_tet4(jcell_tet4),&
!!$                        "'th cell: SGESVD returned INFO = ",INFO
!!$                   call add(errmsg,2,errstr,myname)
!!$                   goto 1
!!$                endif ! INFO/=0
!!$                ! apply the same criterion as matlab-2013-a does
!!$                !rank_A_tet4 = count(S_A_tet4 >= S_A_tet4(1)*7.*eps)
!!$                rank_A_tet4 = count(S_A_tet4 >= S_A_tet4(1)*7.e-5)
!!$!
!!$                ! compute rank of Ab_tet4
!!$                call SGESVD( 'N', 'N', 5 , 7, Ab_tet4, 5, S_Ab_tet4, U, 1, VT, 1, WORK_SGESVD_Ab_tet4, LWORK_SGESVD_Ab_tet4, INFO )
!!$                if(INFO/=0) then
!!$                   j = mod(idx_add(inb),4); if(j==0) j=4
!!$                   write(errstr,*) "computation of singular values of Ab_tet4 failed comparing ",iface,"'th face of ",&
!!$                        icell,"'th cell and ",j,"'th face of ",this%iglob_tet4(jcell_tet4),&
!!$                        "'th cell: SGESVD returned INFO = ",INFO
!!$                   call add(errmsg,2,errstr,myname)
!!$                   goto 1
!!$                endif ! INFO/=0
!!$                ! apply the same criterion as matlab-2013-a does:
!!$                !rank_Ab_tet4 = count(S_Ab_tet4 >= S_Ab_tet4(1)*7.*eps)
!!$                rank_Ab_tet4 = count(S_Ab_tet4 >= S_Ab_tet4(1)*7.e-5)
!!$!
!!$if(icell==3) then
!!$if(idx_add(inb)==17) then
!!$print *, "cell 3, face ",iface,", SVD with cell 5 face 1:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==18) then
!!$print *, "cell 3, face ",iface,", SVD with cell 5 face 2:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==19) then
!!$print *, "cell 3, face ",iface,", SVD with cell 5 face 3:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==20) then
!!$print *, "cell 3, face ",iface,", SVD with cell 5 face 4:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==45) then
!!$print *, "cell 3, face ",iface,", SVD with cell 12 face 1:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==46) then
!!$print *, "cell 3, face ",iface,", SVD with cell 12 face 2:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==47) then
!!$print *, "cell 3, face ",iface,", SVD with cell 12 face 3:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$if(idx_add(inb)==48) then
!!$print *, "cell 3, face ",iface,", SVD with cell 12 face 4:"
!!$print *, rank_A_tet4," S_A",S_A_tet4
!!$print *, rank_Ab_tet4," S_Ab",S_Ab_tet4
!!$end if
!!$endif ! icell== 3
!!$!print *, "S_A",S_A_tet4
!!$!print *, "S_Ab",S_Ab_tet4
!!$                if(rank_A_tet4 == rank_Ab_tet4 .and. rank_A_tet4 > 0 .and. rank_A_tet4 <=4) then
!!$                   idx_add(inb) = this%iglob_tet4(jcell_tet4)
!!$                   ! check if this neighbour cell was already found (might happen, if edge neighbours are detected)
!!$                   if(associated(idx)) then
!!$                      if(any(idx == idx_add(inb))) then
!!$                         idx_add(inb) = -1
!!$                         nnb_add = nnb_add - 1
!!$                      end if
!!$                   end if
!!$                else
!!$                   idx_add(inb) = -1
!!$                   nnb_add = nnb_add - 1
!!$                end if
!!$             end do ! inb
!##################################################################################
!##################################################################################
! THE ABOVE DOES  N O T  WORK!!, see above
!##################################################################################
!##################################################################################
!
!
!##################################################################################
! INSTEAD
! follow the approach of testing signed areas of triangles which are spanned by 
! one point of the first tri face and an edge of the other tri face
! as described in section 4 of paper
!   "Faster Triangle-Triangle Intersection Tests", Olivier Devillers, Philippe Guige, 
!   INRIA Sophia Antipolis research report RR 4488, June 2002, http://hal.inria.fr/inria-00072100,
!   PDF: https://hal.inria.fr/inria-00072100/file/RR-4488.pdf
! in short here referred to as Tri Tri Test = TTT, the nomenclature p,q,r etc.
! is used as in the paper
! however, we modify this algorithm a little, in order to exclude edge- and point- 
! neighbours, in order to soley detect face-neighbours
!##################################################################################
             ! define variables related to current face (p1,q1,r1), as needed for procedure in paper
             select case(iface)
             case(1)
                n_TTT = p1
                p1_TTT = p2
                q1_TTT = p4
                r1_TTT = p3
             case(2)
                n_TTT = p2
                p1_TTT = p1
                q1_TTT = p3
                r1_TTT = p4
             case(3)
                n_TTT = p3
                p1_TTT = p1
                q1_TTT = p4
                r1_TTT = p2
             case(4)
                n_TTT = p4
                p1_TTT = p1
                q1_TTT = p2
                r1_TTT = p3
             end select
!
             allocate(idx_add(nnb_potential))
             idx_add = pack( (/ (j,j=1,4*this%ncell_tet4) /) , is_potential_neighbour_face )
             nnb_add = nnb_potential
!!$if(icell==56 .and. iface==3) then
!!$print *, "##########################################################################"
!!$print *, "##########################################################################"
!!$print *,"icell==56, iface ",iface,"; n_TT = ",n_TTT,", idx_add = ",idx_add
!!$print *, "p1 = ",p1_TTT
!!$print *, "q1 = ",q1_TTT
!!$print *, "r1 = ",r1_TTT
!!$end if
!
             ! loop on all potential neighbour tri faces (for which already was assured that they are 
             ! in the same plane and have the correct "outside" direction") and check if the potential triangle (tri2)
             ! overlaps with triangle p1_TTT,q1_TTT,r1_TTT (tri1)
             do inb = 1,nnb_potential
!!$             do inb = 1,1!nnb_potential,nnb_potential !! FS FS
!!$             do inb = nnb_potential,nnb_potential !! FS FS
!
                ! select the respective corners of tri2, BUT with changed orientation:
                ! for the algorithm in the paper, we need an anticlockwise orientation in the SAME
                ! plane and from the same viewpoint as tri1 (p1,q1,r1), which is oriented anticlockwisely from the viewpoint of 
                ! the tet4-node opposite tri1 (stored in n_TTT) (i.e. looking from the inside of the tetrahedron)
                ! so here for tri2 we need an anticlockwise orientation from viewpoint from n_TTT as well, which is on the OUTSIDE 
                ! side of tri2 (p2,q2,r2), i.e. change the orientation of tri2, compared with normal orientation
                select case(mod(idx_add(inb),4))
                case(1) ! means face 1
                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
                   p2_TTT = this%point(:,this%cell_tet4(2,jcell_tet4))
                   q2_TTT = this%point(:,this%cell_tet4(3,jcell_tet4)) !, note the different orientation compared with tri1, interchange two indices
                   r2_TTT = this%point(:,this%cell_tet4(4,jcell_tet4))
!!$if(icell==56 .and. iface==1) then
!!$print *, "   inb = ",inb," is face 1"
!!$print *, "   p2 = ",p2_TTT
!!$print *, "   q2 = ",q2_TTT
!!$print *, "   r2 = ",r2_TTT
!!$endif
                case(2) ! means face 2
                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
                   p2_TTT = this%point(:,this%cell_tet4(1,jcell_tet4))
                   q2_TTT = this%point(:,this%cell_tet4(4,jcell_tet4))
                   r2_TTT = this%point(:,this%cell_tet4(3,jcell_tet4))
!!$if(icell==56 .and. iface==1 .and.idx_add(inb)==170) then
!!$print *, "   inb = ",inb," is face 2"
!!$print *, "   p2 = ",p2_TTT
!!$print *, "   q2 = ",q2_TTT
!!$print *, "   r2 = ",r2_TTT
!!$endif
                case(3) ! means face 3
                   jcell_tet4 = idx_add(inb)/4 + 1 ! tet4 cell index of cell, to which this face belongs
                   p2_TTT = this%point(:,this%cell_tet4(1,jcell_tet4))
                   q2_TTT = this%point(:,this%cell_tet4(2,jcell_tet4))
                   r2_TTT = this%point(:,this%cell_tet4(4,jcell_tet4))
                case(0) ! means face 4, as in that case the index on all faces is a multiple of 4
                   jcell_tet4 = idx_add(inb)/4 ! tet4 cell index of cell, to which this face belongs
                   p2_TTT = this%point(:,this%cell_tet4(1,jcell_tet4))
                   q2_TTT = this%point(:,this%cell_tet4(3,jcell_tet4))
                   r2_TTT = this%point(:,this%cell_tet4(2,jcell_tet4))
!!$if(icell==56 .and. iface==3 .and.idx_add(inb)==172) then
!!$print *, "   inb = ",inb," is face 4"
!!$print *, "   p2 = ",p2_TTT
!!$print *, "   q2 = ",q2_TTT
!!$print *, "   r2 = ",r2_TTT
!!$endif
                end select
!
             select case(locateP1TTTEcartInversionGrid(p1_TTT,p2_TTT,q2_TTT,r2_TTT,n_TTT))
                ! after call to locateP1TTTEcartInversionGrid, p2_TTT,q2_TTT,r2_TTT are permuted circularly in order to 
                ! meet the conditions for decision trees below
             case(-1)
                ! if p1_TTT could not be located properly, the algorithm does not work
                ! in theory this should not happen!
                j = mod(idx_add(inb),4); if(j==0) j=4
                write(errstr,*) "Triangle Triangle Test failed comparing ",iface,"'th face of ",&
                     icell,"'th cell and ",j,"'th face of ",this%iglob_tet4(jcell_tet4),&
                     "'th cell: could not find a start of the algorithm (should not happen in theory)"
                call add(errmsg,2,errstr,myname)
                goto 1
             case(0)
!!$if(icell==56 .and. iface==3 .and.idx_add(inb)==172) then
!!$print *, "locate P1 TTT is 0 and"
!!$print *, "new p1 = ",p1_TTT
!!$print *, "new p2 = ",p2_TTT
!!$print *, "new q2 = ",q2_TTT
!!$print *, "new r2 = ",r2_TTT
!!$endif
                ! locate = 0 means that p1_TTT is strictly inside tri2, hence tri1 and tri2 overlap, so are truely neighbours
                ! so mark it as a neighbour
                idx_add(inb) = this%iglob_tet4(jcell_tet4)
                if(associated(idx)) then
                   if(any(idx == idx_add(inb))) then
                      idx_add(inb) = -1
                      nnb_add = nnb_add - 1
                   end if
                end if
             case(1)
!!$if(icell==56 .and. iface==3 .and.idx_add(inb)==172) then
!!$print *, "locate P1 TTT is 1 and"
!!$print *, "new p1 = ",p1_TTT
!!$print *, "new p2 = ",p2_TTT
!!$print *, "new q2 = ",q2_TTT
!!$print *, "new r2 = ",r2_TTT
!!$endif
                ! use decision tree "p1_TTT in R1" (confer paper, Figure 9) to evaluate if tris overlap
                if(overlapR1TTTEcartInversionGrid(p1_TTT,q1_TTT,r1_TTT,p2_TTT,r2_TTT,n_TTT)) then
!print *, "and overlap R1 is true"
                   idx_add(inb) = this%iglob_tet4(jcell_tet4)
                   if(associated(idx)) then
                      if(any(idx == idx_add(inb))) then
                         idx_add(inb) = -1
                         nnb_add = nnb_add - 1
                      end if
                   end if
                else
!print *, "and overlap R1 is false"
                   ! if not, remove this neighbour tri face from list of potential neighbours
                   idx_add(inb) = -1    
                   nnb_add = nnb_add - 1
                end if
             case(2)
!!$if(icell==56 .and. iface==3 .and.idx_add(inb)==172) then
!!$print *, "locate P1 TTT is 2 and"
!!$print *, "new p1 = ",p1_TTT
!!$print *, "new p2 = ",p2_TTT
!!$print *, "new q2 = ",q2_TTT
!!$print *, "new r2 = ",r2_TTT
!!$endif
                ! use decision tree "p1_TTT in in R2" (confer paper, Figure 10) to evaluate if tris overlap
                if(overlapR2TTTEcartInversionGrid(p1_TTT,q1_TTT,r1_TTT,p2_TTT,q2_TTT,r2_TTT,n_TTT)) then
!print *, "and overlap R2 is true"
                   idx_add(inb) = this%iglob_tet4(jcell_tet4)
                   if(associated(idx)) then
                      if(any(idx == idx_add(inb))) then
                         idx_add(inb) = -1
                         nnb_add = nnb_add - 1
                      end if
                   end if
                else
!print *, "and overlap R2 is false"
                   ! if not, remove this neighbour tri face from list of potential neighbours
                   idx_add(inb) = -1    
                   nnb_add = nnb_add - 1
                end if
             end select
!
!if(icell==56 .and. iface==3.and.inb==2) stop
          end do ! inb
!
!print *, "nnb_add = ",nnb_add
          if(nnb_add>0) then
                idx => reallocate(idx,nnb+nnb_add)
                if(nnb_add<nnb_potential) then
                   idx(nnb+1:nnb+nnb_add) = pack(idx_add,idx_add>0)
                else
                   idx(nnb+1:nnb+nnb_add) = idx_add
                end if
                nnb = nnb + nnb_add
             end if
             deallocate(idx_add)
!if(icell==56 .and. iface==3) stop
          end do ! iface
          call associateVectorPointer(this%face_neighbour(icell),idx); nullify(idx)
       case(8); cycle
       end select
    end do ! icell
!
1   if(allocated(outward_normal_tet4)) deallocate(outward_normal_tet4)
    if(allocated(origin_distance_tet4)) deallocate(origin_distance_tet4)
    if(allocated(dot)) deallocate(dot)
    if(allocated(is_potential_neighbour_face)) deallocate(is_potential_neighbour_face)
!!$    if(allocated(WORK_SGESVD_A_tet4)) deallocate(WORK_SGESVD_A_tet4)
!!$    if(allocated(WORK_SGESVD_Ab_tet4)) deallocate(WORK_SGESVD_Ab_tet4)
    if(allocated(idx_add)) deallocate(idx_add)
  end subroutine createFaceNeighboursEcartInversionGrid
!------------------------------------------------------------------------
!> \brief find a permutation of p2,q2,r2 to start the search algorithm as in paper by Devillers/Guigue
!! \details Confer paper by Devillers/Guigue, Figure 6: identify the location of p1 as either 
!!  "inside" tri2, or in region R1 (including boundary to region "+++"), or in region R1
!!  including boundary to R1 and boundary to region "+-+", hence including corner r2
!!  if necessary, p2,q2,r2 are permuted circularly on exit (just renaming the points, actually)
!!  and a situation as in Figure 6 of the paper is constructed
  function locateP1TTTEcartInversionGrid(p1,p2,q2,r2,n) result(iflag)
    real, dimension(3) :: p1,p2,q2,r2,n
    integer :: iflag
    ! local
    integer :: j_times
    real, dimension(3) :: tmp
!
    iflag = -1
!
    ! conduct at most 3 searches, permuting the points after each iteration if not successfull
    ! the fourth permutation would regain the start configuration of this loop
    do j_times = 1,3
       if(checkOrientationTet4EcartInversionGrid(p1,p2,q2,n) > 0) then
          ! R1, R2, inside, "+-+"
          if(checkOrientationTet4EcartInversionGrid(p1,q2,r2,n) > 0) then
             ! R1, inside
             if(checkOrientationTet4EcartInversionGrid(p1,r2,p2,n) > 0) then
                ! inside
                iflag = 0; return
             else
                ! R1 including boundary to inside ("+++")
                iflag = 1; return
             end if
          else
             ! R2, "+-+"
             if(checkOrientationTet4EcartInversionGrid(p1,r2,p2,n) <= 0) then
                ! R2 including boundary to R1 and boundary to "+-+", hence including corner r2
                iflag = 2; return
             end if
          end if
       end if
       ! if loop comes here, p1 was not yet located, so permute points p2,q2,r2 circularly
       ! (if not for the 3rd time) and do this again 
       if(j_times < 3) then
          tmp=r2; r2=q2; q2=p2; p2=tmp
       end if
    end do ! j_times
  end function locateP1TTTEcartInversionGrid
!------------------------------------------------------------------------
!> \brief go through decision tree Figure 9, in paper by Devillers/Guigue
!! \details note, however, that this routine is also called, if p1 lies on the boundary of
!!  "+++" and R1 (excluding! poins r1,p2)
!!  and it is not completely analogous to decision tree, as here all "=" cases in inequalities
!!  were switched (i.e. ">=" / "<" decisions became ">" / "<=" decisions) in order to prevent the 
!!  algorithm to find edge and point neighbours
  function overlapR1TTTEcartInversionGrid(p1,q1,r1,p2,r2,n) result(overlap)
    real, dimension(3) :: p1,q1,r1,p2,r2,n
    logical :: overlap
    real, parameter :: eps = epsilon(1.0)
!
    ! decision I (confer paper)
    if(checkOrientationTet4EcartInversionGrid(r2,p2,q1,n) > 0) then
!
       ! IIa
       if(checkOrientationTet4EcartInversionGrid(r2,p1,q1,n) <= 3*eps) goto 1
!
       ! IIIa
       if(checkOrientationTet4EcartInversionGrid(p1,p2,q1,n) > 3*eps) goto 2
!
       ! IVa
       if(checkOrientationTet4EcartInversionGrid(p1,p2,r1,n) <= 3*eps) goto 1
!
       ! V
       if(checkOrientationTet4EcartInversionGrid(q1,r1,p2,n) <= 3*eps) goto 1
!
       ! if routine comes here, we have a positive outcome in decision V, so return overlap=true
       goto 2
!
    else ! I
!
       ! IIb
       if(checkOrientationTet4EcartInversionGrid(r2,p2,r1,n) <= 3*eps) goto 1!then;print *,"R1 TTT IIb";goto 1;endif !! FS FS
!
       ! IIIb
       if(checkOrientationTet4EcartInversionGrid(q1,r1,r2,n) <= 3*eps) goto 1
!
       ! IVb -> ATTENTION, WRONG CONDITION IN THE PAPER (inequality signs must be the other way round, as in IIIb, compare text int he paper, page 13, end of 1st paragraph)
       if(checkOrientationTet4EcartInversionGrid(p1,p2,r1,n) <= 3*eps) goto 1
!
       ! if routine comes here, we have a positive outcome in decision IVb, so return overlap=true
       goto 2
!
    end if ! I
!
1   overlap = .false.; return
2   overlap = .true.; return
  end function overlapR1TTTEcartInversionGrid
!------------------------------------------------------------------------
!> \brief go through decision tree Figure 10, in paper by Devillers/Guigue
!! \details note, however, that this routine is also called, if p1=r1
!!  and it is not completely analogous to decision tree, as here all "=" cases in inequalities
!!  were switched (i.e. ">=" / "<" decisions became ">" / "<=" decisions) in order to prevent the 
!!  algorithm to find edge and point neighbours
  function overlapR2TTTEcartInversionGrid(p1,q1,r1,p2,q2,r2,n) result(overlap)
    real, dimension(3) :: p1,q1,r1,p2,q2,r2,n
    logical :: overlap
    real, parameter :: eps = epsilon(1.0)
!
    ! decision I
    if(checkOrientationTet4EcartInversionGrid(r2,p2,q1,n) > 0) then
!
       ! IIa
       if(checkOrientationTet4EcartInversionGrid(q2,r2,q1,n) > 0) then
!
          ! IIIa
          if(checkOrientationTet4EcartInversionGrid(p1,p2,q1,n) > 0) then
             ! IVa
             if(checkOrientationTet4EcartInversionGrid(p1,q2,q1,n) >= -3*eps) goto 1
             goto 2 ! positive decision in IVa, return true
!
          else ! IIIa
             ! IVb
             if(checkOrientationTet4EcartInversionGrid(p1,p2,r1,n) <= 3*eps) goto 1
!
             ! Va
             if(checkOrientationTet4EcartInversionGrid(r2,p2,r1,n) <= 3*eps) goto 1
             goto 2 ! positive decision in Va, return true
          end if ! IIIa
!
       else ! IIa
!
          ! IIIb
          if(checkOrientationTet4EcartInversionGrid(p1,q2,q1,n) >= -3*eps) goto 1
!
          ! IVc
          if(checkOrientationTet4EcartInversionGrid(q2,r2,r1,n) <= 3*eps) goto 1
!
          ! Vb
          if(checkOrientationTet4EcartInversionGrid(q1,r1,q2,n) <= 3*eps) goto 1
          goto 2 ! positive decision in Vb, return true
!
       end if ! IIa
!
    else ! I
!
       ! IIb
       if(checkOrientationTet4EcartInversionGrid(r2,p2,r1,n) <= 3*eps) goto 1
!
       ! IIIc
       if(checkOrientationTet4EcartInversionGrid(q1,r1,r2,n) > 0) then
          ! IVd
          if(checkOrientationTet4EcartInversionGrid(r1,p1,p2,n) <= 3*eps) goto 1
          goto 2 ! positive decision in IVd, return true
!
       else ! IIIc
          ! IVe
          if(checkOrientationTet4EcartInversionGrid(q1,r1,q2,n) <= 3*eps) goto 1
!
          ! Vc
          if(checkOrientationTet4EcartInversionGrid(q2,r2,r1,n) <= 3*eps) goto 1
          goto 2 ! positive decision in Vc, return true
       end if ! IIIc
!
    end if ! I
!
1   overlap = .false.; return
2   overlap = .true.; return
  end function overlapR2TTTEcartInversionGrid
!------------------------------------------------------------------------
!! \param checks orientation of tet4 cell
!! \details computes the determinant
!!           ( p2(1) p3(1) p4(1) p1(1) )
!!           ( p2(2) p3(2) p4(2) p1(2) )       ( p2(1)-p1(1)  p3(1)-p1(1)  p4(1)-p1(1) )
!! ori = det ( p2(3) p3(3) p4(3) p1(3) ) = det ( p2(2)-p1(2)  p3(2)-p1(2)  p4(2)-p1(2) )
!!           (  1     1     1     1    )       ( p2(3)-p1(3)  p3(3)-p1(3)  p4(3)-p1(3) )
!! ori > 0 : p1,p2,p3 have counter-clockwise orientation in their plane from viewpoint of p4
!! ori < 0 : p1,p2,p3 have clockwise orientation in their plane from viewpoint of p4
!! ori = 0 : p1,p2,p3,p4 are coplanar
!
  function checkOrientationTet4EcartInversionGrid(p1,p2,p3,p4) result(ori)
    real, dimension(3) :: p1,p2,p3,p4
    real :: ori
    ori = ( (p2(1)-p1(1)) * (p3(2)-p1(2)) * (p4(3)-p1(3)) + &
            (p3(1)-p1(1)) * (p4(2)-p1(2)) * (p2(3)-p1(3)) + &
            (p4(1)-p1(1)) * (p2(2)-p1(2)) * (p3(3)-p1(3)) ) &
        - ( (p4(1)-p1(1)) * (p3(2)-p1(2)) * (p2(3)-p1(3)) + &
            (p2(1)-p1(1)) * (p4(2)-p1(2)) * (p3(3)-p1(3)) + &
            (p3(1)-p1(1)) * (p2(2)-p1(2)) * (p4(3)-p1(3)) )
  end function checkOrientationTet4EcartInversionGrid
!------------------------------------------------------------------------
!> \brief write face neighbours of external Cartesian inversion grid to file
!! \param this external Cartesian inversion grid
!! \param filename filename of neighbours file
!! \param is_binary logical indicating whether file is binary or not
!! \param lu file unit to write file
!! \param errmsg error message
!
  subroutine writeFaceNeighbourEcartInversionGrid(this,filename,is_binary,file_exists,lu,errmsg)
    type (ecart_inversion_grid) :: this
    character(len=*) :: filename
    logical :: is_binary,file_exists
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=36) :: myname = 'writeFaceNeighbourEcartInversionGrid'
    character (len=7) :: open_status
    integer :: ios,icell
    integer, dimension(:), pointer :: nb
    integer, parameter :: zero = 0
!
    nullify(nb)
    call addTrace(errmsg,myname)
!
    if(file_exists) then
       open_status = 'REPLACE'
    else
       open_status = 'NEW'
    end if
!
    if(is_binary) then
       open(unit=lu,file=filename,form='unformatted',access='stream',status=open_status,action='write',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open binary neighbours file '"//trim(filename)//"' to write",myname)
          return
       end if
       write(lu) this%ncell
       do icell = 1,this%ncell
          nb => getVectorPointer(this%face_neighbour(icell))
          if(associated(nb)) then
             write(lu) size(nb),nb
          else
             write(lu) zero
          end if
       end do ! icell
!
    else ! is_binary
!
       open(unit=lu,file=filename,form='formatted',status=open_status,action='write',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open ascii neighbours file '"//trim(filename)//"' to write",myname)
          return
       end if
       write(lu,*) this%ncell
       do icell = 1,this%ncell
          nb => getVectorPointer(this%face_neighbour(icell))
          if(associated(nb)) then
             write(lu,*) size(nb),nb
          else
             write(lu,*) zero
          end if
       end do ! icell
!
    end if ! is_binary
!
    close(lu)
  end subroutine writeFaceNeighbourEcartInversionGrid
!------------------------------------------------------------------------
!> \brief read external Cartesian inversion grid to file
!! \param this external Cartesian inversion grid
!! \param filename filename of inversion grid file
!! \param is_binary logical indicating whether file is binary or not
!! \param lu file unit to read file
!! \param errmsg error message
!
  subroutine readFaceNeighbourEcartInversionGrid(this,filename,is_binary,lu,errmsg)
    type (ecart_inversion_grid) :: this
    character(len=*) :: filename
    logical :: is_binary
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=35) :: myname = 'readFaceNeighbourEcartInversionGrid'
    integer :: ios,icell,ncell,nnb
    integer, dimension(:), pointer :: nb
    character(len=400) :: line
!
    nullify(nb)
    call addTrace(errmsg,myname)
    allocate(this%face_neighbour(this%ncell))
!
    if(is_binary) then
       open(unit=lu,file=filename,form='unformatted',access='stream',status='old',action='read',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open binary neighbours file '"//trim(filename)//"' to read",myname)
          return
       end if
       read(lu,iostat=ios) ncell
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not read number of cells as first entry in binary neighbours file '"&
               //trim(filename)//"'",myname)
          return
       end if
       if(ncell/=this%ncell) then
          close(lu)
          write(errstr,*) "number of cells ",ncell," contained in binary neighbours file '"//trim(filename)//&
               "' does not match number of cells ",this%ncell," constructed from cell connectivity files"
          call add(errmsg,2,errstr,myname)
          return
       end if
       do icell = 1,ncell
          read(lu,iostat=ios) nnb
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read number of neighbours of ",icell,&
                  "'th cell from binary neighbours file '"//trim(filename)//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
          if(nnb>0) then
             allocate(nb(nnb))
             read(lu,iostat=ios) nb
             if(ios/=0) then
                close(lu); deallocate(nb)
                write(errstr,*) "could not read ",nnb," neighbour indices of ",icell,&
                     "'th cell from binary neighbours file '"//trim(filename)//"'"
                call add(errmsg,2,errstr,myname)
                return
             end if
             if(any(nb<1 .or. nb>this%ncell)) then
                close(lu)
                write(errstr,*) "there are neighbour indices out of range for ",icell,&
                     "'th cell; must be >= 1 and <= ncell =",this%ncell
                call add(errmsg,2,errstr,myname)
                return
             end if
             call associateVectorPointer(this%face_neighbour(icell),nb)
             nullify(nb)
          end if
       end do ! icell
!
    else ! is_binary
!
       open(unit=lu,file=filename,form='formatted',status='old',action='read',iostat=ios)
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not open ascii neighbours file '"//trim(filename)//"' to read",myname)
          return
       end if
       read(lu,*,iostat=ios) ncell
       if(ios/=0) then
          close(lu)
          call add(errmsg,2,"could not read number of cells as first entry in ascii neighbours file '"&
               //trim(filename)//"'",myname)
          return
       end if
       if(ncell/=this%ncell) then
          close(lu)
          write(errstr,*) "number of cells ",ncell," contained in ascii neighbours file '"//trim(filename)//&
               "' does not match number of cells ",this%ncell," constructed from cell connectivity files"
          call add(errmsg,2,errstr,myname)
          return
       end if
       do icell = 1,ncell
          read(lu,"(a400)",iostat=ios) line
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read ",icell,"'th line from ascii neighbours file '"//trim(filename)//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
          read(line,*,iostat=ios) nnb
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read number of neighbours of ",icell,&
                  "'th cell from ascii neighbours file '"//trim(filename)//"'"
             call add(errmsg,2,errstr,myname)
             return
          end if
          if(nnb>0) then
             allocate(nb(nnb))
             read(line,*,iostat=ios) nnb,nb
             if(ios/=0) then
                close(lu); deallocate(nb)
                write(errstr,*) "could not read ",nnb," neighbour indices of ",icell,&
                     "'th cell from ascii neighbours file '"//trim(filename)//"'"
                call add(errmsg,2,errstr,myname)
                return
             end if
             if(any(nb<1 .or. nb>this%ncell)) then
                close(lu)
                write(errstr,*) "there are neighbour indices out of range for ",icell,&
                     "'th cell; must be >= 1 and <= ncell =",this%ncell
                call add(errmsg,2,errstr,myname)
                return
             end if
             call associateVectorPointer(this%face_neighbour(icell),nb)
             nullify(nb)
          end if
       end do ! icell
    end if ! is_binary
  end subroutine readFaceNeighbourEcartInversionGrid
!------------------------------------------------------------------------
!> \brief deallocate external layered Cartesian inversion grid
!! \param this external layered Cartesian inversion grid
!
  subroutine deallocateEcartInversionGrid(this)
    type (ecart_inversion_grid) :: this
    integer :: n
    this%npoint = 0
    if(associated(this%point)) deallocate(this%point)
    this%ncell = 0
    if(associated(this%ctype)) deallocate(this%ctype)
    if(associated(this%icell_type)) deallocate(this%icell_type)
    if(associated(this%face_neighbour)) then
       do n = 1,size(this%face_neighbour)
          call dealloc(this%face_neighbour(n))
       end do
       deallocate(this%face_neighbour)
    end if
    this%ncell_tet4 = 0
    if(associated(this%cell_tet4)) deallocate(this%cell_tet4)
    if(associated(this%iglob_tet4)) deallocate(this%iglob_tet4)
    this%ncell_hex8 = 0
    if(associated(this%cell_hex8)) deallocate(this%cell_hex8)
    if(associated(this%iglob_hex8)) deallocate(this%iglob_hex8)
    this%vtk_geometry_type_int = -1
    this%is_defined = .false.
  end subroutine deallocateEcartInversionGrid
!------------------------------------------------------------------------
!> \brief return overall number of invgrid cells, if any
!
  function getNcellEcartInversionGrid(this) result(ncell)
    type (ecart_inversion_grid), intent(in) :: this
    integer :: ncell
    ncell = this%ncell
  end function getNcellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief transform wp,event or station coords to x,y,z coords for vtk application
!! \details in module inversioGrid it was already checked, if coords_type is one of
!!  'wp','event','station' and that c1,c2,c3 are associated and have all same length
!!  also it can be assured that this ecart inversion grid is properly defined (otherwise
!!  inversionGrid module would not fork here)
!! \param this ecart inversion grid
!! \param c1 vector or first coordinate (contains vtk x-values on exit)
!! \param c2 vector or second coordinate (contains vtk y-values on exit)
!! \param c3 vector or third coordinate (contains vtk z-values on exit)
!! \param coords_type 'wp','event','station'
!! \param errmsg error message
!
  subroutine transformToVtkEcartInversionGrid(this,c1,c2,c3,coords_type,errmsg)
    type (ecart_inversion_grid) :: this
    real, dimension(:), intent(inout) :: c1,c2,c3
    character(len=*) :: coords_type
    type (error_message) :: errmsg
!
    call addTrace(errmsg,'transformToVtkEcartInversionGrid')
!
    ! no need of selecting coords_type: always do the same, since in case of using the 
    ! ecart inversion grid the wavefield points as well as event and station coordinates
    ! are expected as x,y,z coordinates (in that order)
!
    if(this%apply_vtk_coords_scaling_factor) then
       c1 = c1 * this%vtk_coords_scaling_factor
       c2 = c2 * this%vtk_coords_scaling_factor
       c3 = c3 * this%vtk_coords_scaling_factor
    end if
  end subroutine transformToVtkEcartInversionGrid
!------------------------------------------------------------------------
!> \brief return geometry information on cells for vtk output
!
  subroutine getGeometryVtkEcartInversionGrid(this,geometry_type,points,cell_connectivity,cell_type,cell_indx_out,&
       errmsg,cell_indx_req,indx_map_out)
    type (ecart_inversion_grid) :: this
    integer :: geometry_type
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: cell_connectivity,cell_type,cell_indx_out
    integer, dimension(:), pointer, optional :: indx_map_out
    type (error_message) :: errmsg
    ! local
    character(len=32) :: myname = 'getGeometryVtkEcartInversionGrid'
    character(len=400) :: errstr
    integer :: size_cell_indx_req,i,ncell_return,icell,jcell,ncell_con
    logical, dimension(:), allocatable :: valid_non_duplicate_cell_indx_req
    real :: xc,yc,zc
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
    if(present(cell_indx_req)) then
       ! if cell_indx_req is present and there are any indices in valid range, select those specific cells
       ! and remove duplicate indices from cell_indx_req. otherwise return no cells (nullified pointers)
       ! HOWEVER, the general order of the returned cell indices is not changed, only gaps of invalid/duplicate
       ! indices in cell_indx_req are closed! (this approach is different from e.g. scartInversionGrid, where
       ! valid/non-duplicate indices in cell_indx_req are always returned sorted w.r.t. the internal cell order)
       size_cell_indx_req = size(cell_indx_req)
       allocate(valid_non_duplicate_cell_indx_req(size_cell_indx_req))
!
       ! first select valid indices
       valid_non_duplicate_cell_indx_req = cell_indx_req .ge. 1 .and. cell_indx_req .le. this%ncell
!
       ! then make all duplicate indices also invalid
       do i = 1,size_cell_indx_req-1
          if(valid_non_duplicate_cell_indx_req(i)) then
             where(cell_indx_req(i+1:size_cell_indx_req) == cell_indx_req(i))
                valid_non_duplicate_cell_indx_req(i+1:size_cell_indx_req) = .false.
             end where
          end if
       end do ! i
!
       ncell_return = count(valid_non_duplicate_cell_indx_req)
       if(ncell_return == 0) then
          deallocate(valid_non_duplicate_cell_indx_req)
          write(errstr,*) "there are no valid indices among requested cell indices; indices must be between 1 and ",&
               this%ncell
          call add(errmsg,2,errstr,myname)
          return
       endif
!
       ! define arrays cell_indx_out, indx_map_out
       allocate(cell_indx_out(ncell_return))
       cell_indx_out = pack(cell_indx_req , valid_non_duplicate_cell_indx_req)
       if(present(indx_map_out)) then
          allocate(indx_map_out(ncell_return))
          indx_map_out = pack( (/ (i,i=1,size_cell_indx_req) /) , valid_non_duplicate_cell_indx_req)
       end if
!
       deallocate(valid_non_duplicate_cell_indx_req)
!
    else ! present(cell_indx_req)
       ! return all cells in the order of inversion grid cells, hence arrays cell_indx_out
       ! and indx_map_out are trivial
       allocate(cell_indx_out(this%ncell))
       cell_indx_out = (/ (i,i=1,this%ncell) /)
       if(present(indx_map_out)) then
          allocate(indx_map_out(this%ncell))
          indx_map_out = cell_indx_out
       end if
       ncell_return = this%ncell
!
    end if ! present(cell_indx_req)
!
    select case(this%vtk_geometry_type_int)
    case(0) ! CELLS
       ! return all points, regardless of cell_idx_req
       allocate(points(3,this%npoint)); points = this%point
!
       ! define cell_type array
       allocate(cell_type(ncell_return))
       where(this%ctype(cell_indx_out) == 4)
          cell_type = 10
       end where
       where(this%ctype(cell_indx_out) == 8)
          cell_type = 12
       end where
!
       ! define cell_connectivity by looping over cell_indx_out and adding respective point indices
       ncell_con = ncell_return + 4*count(cell_type == 10) + 8*count(cell_type == 12)
       allocate(cell_connectivity(ncell_con))
!
       ! loop on cell_indx_out
       ncell_con = 0 ! reset counter on entries in connectivity array
       do icell = 1,ncell_return; jcell = cell_indx_out(icell)
          select case(this%ctype(jcell))
          case(4)
             ! first store number of point indices to come for this cell
             ncell_con=ncell_con+1; cell_connectivity(ncell_con) = 4
             ! define the required points for the cell type of this cell (4-node tetrahedron)
             ! indices in vtk files have offset 0, so always store ipoint -1
             cell_connectivity(ncell_con+1:ncell_con+4) = this%cell_tet4(:,this%icell_type(jcell)) - 1
             ncell_con = ncell_con + 4
          case(8)
             ! first store number of point indices to come for this cell
             ncell_con=ncell_con+1; cell_connectivity(ncell_con) = 8
             ! define the required points for the cell type of this cell (8-node hexahedron)
             ! indices in vtk files have offset 0, so always store ipoint -1
             cell_connectivity(ncell_con+1:ncell_con+8) = this%cell_hex8(:,this%icell_type(jcell)) - 1
             ncell_con = ncell_con + 8
          end select
       end do ! icell
       geometry_type = 0
!
    case(1) ! CELL CENTERS
       ! do not define arrays cell_connectivity and cell_type in this case. Only return points as cell centers
       allocate(points(3,ncell_return))
       do icell = 1,ncell_return
          call getCenterCellEcartInversionGrid(this,cell_indx_out(icell),xc,yc,zc)
          points(1,icell) = xc
          points(2,icell) = yc
          points(3,icell) = zc
       end do ! icell
       geometry_type = 1
    end select ! this%vtk_geometry_type
!
    if(this%apply_vtk_coords_scaling_factor) then
       points = points*this%vtk_coords_scaling_factor
    end if
!
  end subroutine getGeometryVtkEcartInversionGrid
!------------------------------------------------------------------------
!> \brief return indices of all face neighbours for all (optionally only subset of) cells
!! \param this ecart inversion grid
!! \param nb_idx pointer to array of length this%ncell which contains
!
  subroutine getIndicesFaceNeighboursEcartInversionGrid(this,nb_idx)
    type (ecart_inversion_grid) :: this
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx
    ! local
    integer :: icell
    integer, dimension(:), pointer :: nb,nb2
!
    nullify(nb,nb2)
    nullify(nb_idx)
    if(.not.this%is_defined) return
!
    allocate(nb_idx(this%ncell))
    do icell = 1,this%ncell
       nb => getVectorPointer(this%face_neighbour(icell))
       if(associated(nb)) then
          allocate(nb2(size(nb))); nb2 = nb
          call associateVectorPointer(nb_idx(icell),nb2)
          nullify(nb,nb2)
       end if
    end do ! icell
!
  end subroutine getIndicesFaceNeighboursEcartInversionGrid
!------------------------------------------------------------------------
!> \brief for each cell return indices of wavefield points contained in that cell
!! \param this ecart inversion grid
!! \param x vector or first coordinate of wavefield points
!! \param y vector or second coordinate of wavefield points
!! \param z vector or third coordinate of wavefield points
!! \param wp_idx pointer to array of length this%ncell; if invgrid not defined yet, nullified on exit
!! \param errmsg error message
!
  subroutine locateWpInsideEcartInversionGrid(this,x,y,z,wp_idx,errmsg)
    type (ecart_inversion_grid) :: this
    real, dimension(:), intent(in) :: x,y,z
    type (integer_vector_pointer), dimension(:), pointer :: wp_idx
    type (error_message) :: errmsg
    ! local
    character(len=32) :: myname = 'locateWpInsideEcartInversionGrid'
    character(len=400) :: errstr
    integer :: nwp,nwp_to_locate,i,icell
    integer, dimension(:), pointer :: idx
    integer, dimension(:), pointer :: idx_to_locate,ptmp
    real, dimension(:,:), allocatable :: wp_coords
!
    nullify(idx,idx_to_locate,ptmp)
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
!
    ! in routine locateWpInsideInversionGrid of module inversionGrid it was already assured
    ! that x,y,z contain values and are all of same length!
    nwp = size(x)
!
    ! start with all points not having been located inside cells (obviously)
    nwp_to_locate = nwp
    allocate(idx_to_locate(nwp_to_locate))
    idx_to_locate = (/ (i,i=1,nwp_to_locate) /)
    allocate(wp_coords(nwp_to_locate,3))
    wp_coords(:,1) = x; wp_coords(:,2) = y; wp_coords(:,3) = z
!
    ! for each cell, find its containing wavefield points and remove their indices from idx_to_locate
    do icell = 1,this%ncell
       ! if all wavefield points already have been located in other cells, the remaining cells are empty, 
       ! hence there is nothing left to do
       if(nwp_to_locate == 0) exit
!
       ! find all indices (among idx_to_locate) for which the wavefield points are inside this cell
       select case(this%ctype(icell))
       case(4); idx=>locateWpInsideTet4CellEcartInversionGrid(this,this%icell_type(icell),wp_coords(idx_to_locate,:))
       case(8); idx=>locateWpInsideHex8CellEcartInversionGrid(this,this%icell_type(icell),wp_coords(idx_to_locate,:))
          write(errstr,*) "for 8-node hexahedral invgrid cells (as cell ",icell,&
               "), localization of wavefield points is not supported yet"
          call add(errmsg,2,errstr,myname); return
       end select
!
       if(associated(idx)) then
          ! as returned by locateWpInside***CellEcartInversionGrid, idx contains indices w.r.t. the rows of wp_coords(idx_to_locate,:),
          ! hence need to assign to idx the correct wavefield point indices, which are idx_to_locate(idx)
          call allocateVectorPointer(wp_idx(icell),size(idx))
          call fillVectorPointer(wp_idx(icell),idx_to_locate(idx),1)
!
          ! remove the located wavefield point indices from array idx_to_locate
          nwp_to_locate = nwp_to_locate - size(idx)
          if(nwp_to_locate > 0) then
             idx_to_locate(idx) = -1
             allocate(ptmp(nwp_to_locate)); ptmp = pack(idx_to_locate,idx_to_locate>0)
             deallocate(idx_to_locate); idx_to_locate => ptmp; nullify(ptmp)
          else
             deallocate(idx_to_locate)
          end if
!
          deallocate(idx)
       end if ! associated(idx)
    end do ! icell
!
    if(nwp_to_locate > 0) then
       write(errstr,*) nwp_to_locate," wavefield points (out of ",nwp,&
            ") could not be located inside the inversion grid"
       call add(errmsg,1,errstr,myname)
    end if
!
  end subroutine locateWpInsideEcartInversionGrid
!------------------------------------------------------------------------
!> \brief for given tet4 inversion grid cell, select indices of points containing the cell
!! \param this ecart inverison grid
!! \param icell_tet4 index in this%cell_tet4 corresponding to the cell under consideration
!! \param wp_coords (nwp,3)-array of x,y,z coordinates of wavefield points to be located
!! \return idx array of indices (index of first dimension of wp_coords) of points contained in this cell
!
  function locateWpInsideTet4CellEcartInversionGrid(this,icell_tet4,wp_coords) result(idx)
    type (ecart_inversion_grid) :: this
    integer :: icell_tet4
    real, dimension(:,:) :: wp_coords
    integer, dimension(:), pointer :: idx
    ! local
    integer, dimension(:), pointer :: idx_wp_inside,ptmp
    real, dimension(:), allocatable :: dot
    real :: dot_correction
    integer :: nwp_inside,i,iface
    real, dimension(3) :: outward_direction,p1,p2,p3,p4
!
! BEWARE: THIS ROUTINE MIGHT WELL BE OF NOT VERY GOOD PERFORMANCE! (in terms of memory access in huge array wp_coords, etc)
! PLEASE DON'T HESITATE TO IMPROVE!
!
    nullify(idx_wp_inside,ptmp)
    nullify(idx)
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_tet4(1,icell_tet4))
    p2 = this%point(:,this%cell_tet4(2,icell_tet4))
    p3 = this%point(:,this%cell_tet4(3,icell_tet4))
    p4 = this%point(:,this%cell_tet4(4,icell_tet4))
!
    nwp_inside = size(wp_coords,1)
    allocate(idx_wp_inside(nwp_inside),dot(nwp_inside))
    idx_wp_inside = (/ (i,i=1,nwp_inside) /)
!
!!$print *, "tet4 cell ",icell_tet4," p1,p2,p3,p4 = "
!!$print *, p1
!!$print *, p2
!!$print *, p3
!!$print *, p4
!
    ! FOR ALL FACES DO:
    !
    !   * SELECT ALL THOSE INDICES, FOR WHICH WAVEFIELD POINTS LIE ON THE "INSIDE"-SIDE OF THE FACE:
    !     compare two vectors: 
    !     1) vector n, orthogonal on face pointing outward
    !     2) vector x' = x - p , where vector x is a wavefield point and p lies on the face
    !     if dot(n,x') > 0 , x is strictly on the "outside", i.e. that side of the face, where n points to
    !     FOR SIMPLICITY, do not explicitely compute x', but simplify dot(n,x') = dot(n,x) - dot(n,p)
    !
    !   * THROW AWAY ALL OTHER INDICES
    !
    ! THEN RETURN REMAINING INDICES, IF ANY
!
    do iface = 1,4
!
       ! compute outward pointing vector, which is orthogonal on current face (not necessary to normalize)
       select case(iface)
       case(1)
          ! first face (opposite node p1, outward normal is vector product of vectors p3-p2 and p4-p2)
          outward_direction = vectorProductEcartInversionGrid(p3-p2,p4-p2)
          ! dot_correction = dot(n,p), p2 is on the face (see comment above do loop)
          dot_correction = sum(outward_direction*p2)
!!$print *, "face 1: outward_direction = ",outward_direction," dot_correction = ",dot_correction
       case(2)
          ! second face (opposite node p2, outward normal is vector product of vectors p4-p1 and p3-p1)
          outward_direction = vectorProductEcartInversionGrid(p4-p1,p3-p1)
          ! dot_correction = dot(n,p), p1 is on the face (see comment above do loop)
          dot_correction = sum(outward_direction*p1)
!!$print *, "face 2: outward_direction = ",outward_direction," dot_correction = ",dot_correction
       case(3)
          ! third face (opposite node p3, outward normal is vector product of vectors p2-p1 and p4-p1)
          outward_direction = vectorProductEcartInversionGrid(p2-p1,p4-p1)
          ! dot_correction = dot(n,p), p1 is on the face (see comment above do loop)
          dot_correction = sum(outward_direction*p1)
!!$print *, "face 3: outward_direction = ",outward_direction," dot_correction = ",dot_correction
       case(4)
          ! fourth face (opposite node p4, outward normal is vector product of vectors p3-p1 and p2-p1)
          outward_direction = vectorProductEcartInversionGrid(p3-p1,p2-p1)
          ! dot_correction = dot(n,p), p1 is on the face (see comment above do loop)
          dot_correction = sum(outward_direction*p1)
!!$print *, "face 4: outward_direction = ",outward_direction," dot_correction = ",dot_correction
       end select
!
       ! compute dot products of outward_direction with all (remaining) wavefield point vectors by a general matrix-vector product
!
!!$do i = 1,size(idx_wp_inside)
!!$   print *, "wp_coords(idx_wp_inside(i),:) = ",wp_coords(idx_wp_inside(i),:)
!!$end do
       ! SUBROUTINE SGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
       call SGEMV('N',nwp_inside,3,1.,wp_coords(idx_wp_inside,:),nwp_inside,outward_direction,1,0.,dot,1)
       !dot = matmul(wp_coords(idx_wp_inside,:),outward_direction)
!
       ! dot(n,x') = dot(n,x) - dot(n,p), see comment above do loop
       dot = dot - dot_correction
!!$print *, "dot = ",dot
!
       ! throw away all indices for which dot>0 (the respective points are located strictly outside the tet)
       nwp_inside = count(dot<=0)
       if(nwp_inside>0) then
          allocate(ptmp(nwp_inside)); ptmp = pack(idx_wp_inside,dot<=0)
          deallocate(idx_wp_inside); idx_wp_inside => ptmp
          deallocate(dot); allocate(dot(nwp_inside))
       else
          goto 1
       end if
!
    end do ! iface
!
    ! if program comes here, nwp_inside>0 (as otherwise there was a call to goto 1 above)
    ! and idx_wp_inside contains all indices to be returned 
    idx => idx_wp_inside; nullify(idx_wp_inside)
!!$print *, "idx = ",idx
!
1   if(associated(idx_wp_inside)) deallocate(idx_wp_inside)
    if(allocated(dot)) deallocate(dot)
!!$print *, ""
  end function locateWpInsideTet4CellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief return the vector product of two vectors in R^3
  function vectorProductEcartInversionGrid(a,b,normal) result(c)
    real, dimension(3) :: a,b,c
    logical, optional :: normal
    real :: norm
    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)
    if(present(normal)) then
       if(normal) then
          norm = sqrt(c(1)*c(1) + c(2)*c(2) + c(3)*c(3))
          c(1) = c(1)/norm
          c(2) = c(2)/norm
          c(3) = c(3)/norm
       end if
    end if
  end function vectorProductEcartInversionGrid
!------------------------------------------------------------------------
!> \brief for given hex8 inversion grid cell, select indices of points containing the cell
!! \param this ecart inverison grid
!! \param icell_hex8 index in this%cell_hex8 corresponding to the cell under consideration
!! \param wp_coords (nwp,3)-array of x,y,z coordinates of wavefield points to be located
!! \return idx array of indices (index of first dimension of wp_coords) of points contained in this cell
!
  function locateWpInsideHex8CellEcartInversionGrid(this,icell_hex8,wp_coords) result(idx)
    type (ecart_inversion_grid) :: this
    integer :: icell_hex8
    real, dimension(:,:) :: wp_coords
    integer, dimension(:), pointer :: idx
    nullify(idx)
  end function locateWpInsideHex8CellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in cell icell to standard cell and compute their jacobian
!! \param this ecart inversion grid
!! \param icell global inversion grid cell index
!! \param x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param type_standard_cell defines on return the shape of the standard cell (specific integration weights routine can be chosen): (4=Tetrahedron,6=Hexahedron). If ON INPUT type_standard_cell=-1, then instead of jacobian values actual integration weights are requested
!! \param errmsg error message
!
  subroutine transformToStandardCellEcartInversionGrid(this,icell,x,y,z,jacobian,type_standard_cell,errmsg)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    integer :: type_standard_cell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    type (error_message) :: errmsg
    ! local
    character (len=41) :: myname = 'transformToStandardCellEcartInversionGrid'
!
    call addTrace(errmsg,myname)
!
    if(type_standard_cell == -1) then
       call add(errmsg,2,"Incoming value of type_standard_cell is -1, indicating a request for total "//&
            "integration weights (e.g. used by integration weights of type 6) instead of jacobian values. "//&
            "This functionality is not supported by inversion grids of type ecartInversionGrid",myname)
       return
    end if
!
    ! in routine transformToStandardCellInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined, icell is valid and that x,y,z contain values and are all of same length!
!
    select case (this%ctype(icell))
    case(4)
       type_standard_cell = 4
       call transformToStandardTetEcartInversionGrid(this,icell,x,y,z,jacobian,errmsg)
    case(8)
       type_standard_cell = 6
       call transformToStandardHexEcartInversionGrid(this,icell,x,y,z,jacobian,errmsg)
    end select
!
  end subroutine transformToStandardCellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in tet4 cell icell_tet4 to standard tet and compute their jacobian
!! \param this ecart inversion grid
!! \param icell global cell index of tet4 under consideration
!! \param x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!
  subroutine transformToStandardTetEcartInversionGrid(this,icell,x,y,z,jacobian,errmsg)
    type (ecart_inversion_grid) :: this
    integer :: icell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    type (error_message) :: errmsg
    ! local
    character (len=40) :: myname = 'transformToStandardTetEcartInversionGrid'
    character(len=400) :: errstr
    integer :: np
    real :: det_123,det_124,det_134,det_234
    real, dimension(:), allocatable :: det_X12,det_X13,det_X14,det_X23,det_X24,det_X34
    real :: D0
    real, dimension(:), allocatable :: D2,D3,D4
    real, dimension(3) :: p1,p2,p3,p4
!
    call addTrace(errmsg,myname)
!
    ! in routine transformVectorGlobalToLocalInversionGrid of module inversionGrid it was already assured
    ! that this%is_defined, icell is valid and that x,y,z contain values and are all of same length!
    np = size(x)
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_tet4(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_tet4(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_tet4(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_tet4(4,this%icell_type(icell)))
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! for every wavefield point X=x,y,z compute its representation in coordinates X'=xi,eta,zeta in
! the standard tetrahedron (defined by nodes q1=(0,0,0), q2=(1,0,0), q3=(0,1,0), q4=(0,0,1) )
! where the real nodes pj are mapped to qj (i.e. p1 <-> q1 , ... , p4 <-> q4) and the transformation
! from X' to X shall be
! X(xi,eta,zeta) = p1 + (p2-p1)*xi + (p3-p1)*eta + (p4-p1)*zeta
!                = (1-xi-eta-zeta)*p1 + xi*p2 + eta*p3 + zeta*p4
!
! hence, xi,eta,zeta are related to the BARYCENTRIC COORDINATES bi of this tetrahedron, where an
! arbitrary point X inside the tetrahedron is expressed as X = b1*p1 + b2*p2 + b3*p3 + b4*p4 with b1 + b2 + b3 + b4 = 1, 
! in this case b1=(1-xi-eta-zeta), b2=xi, b3=eta, b4=zeta
! 
! compute the barycentric coordinates bi by solving the linear system
! ( p1(1) p2(1) p3(1) p4(1) )   (b1)   (y)
! ( p1(2) p2(2) p3(2) p4(2) )   (b2)   (x)
! ( p1(3) p2(3) p3(3) p4(3) ) * (b3) = (z)
! (  1     1     1     1    )   (b4)   (1)
! using Cramers rule, which says that the solution may be computed as bi = Di/D0, where:
!
!          ( p1(1) p2(1) p3(1) p4(1) )
!          ( p1(2) p2(2) p3(2) p4(2) )
! D0 = det ( p1(3) p2(3) p3(3) p4(3) )
!          (  1     1     1     1    )
!
!          (  x   p2(1) p3(1) p4(1) )
!          (  y   p2(2) p3(2) p4(2) )
! D1 = det (  z   p2(3) p3(3) p4(3) ) , actuall we don't need D1, as we don't need to calculate b1
!          (  1    1     1     1    )
!
!          ( p1(1)  x   p3(1) p4(1) )
!          ( p1(2)  y   p3(2) p4(2) )
! D2 = det ( p1(3)  z   p3(3) p4(3) )
!          (  1     1    1     1    )
!
!          ( p1(1) p2(1)  x   p4(1) )
!          ( p1(2) p2(2)  y   p4(2) )
! D3 = det ( p1(3) p2(3)  z   p4(3) )
!          (  1     1     1    1    )
!
!          ( p1(1) p2(1) p3(1)  x   )
!          ( p1(2) p2(2) p3(2)  y   )
! D4 = det ( p1(3) p2(3) p3(3)  z   )
!          (  1     1     1     1   )
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! D0  
!
    ! in order to compute D0, we only need 4 sub-determinants (each 3x3): det_ijk = det(M_4l), whereby the 3x3 matrix M_4l
    ! is the so-called "minor" submatrix of above matrix, yield by cuting out row 4 and column l, with l /= i,j,k
    det_123 = p1(1)*p2(2)*p3(3) + p1(2)*p2(3)*p3(1) + p1(3)*p2(1)*p3(2) - &
         ( p1(3)*p2(2)*p3(1) + p1(1)*p2(3)*p3(2) + p1(2)*p2(1)*p3(3))
    det_124 = p1(1)*p2(2)*p4(3) + p1(2)*p2(3)*p4(1) + p1(3)*p2(1)*p4(2) - &
         ( p1(3)*p2(2)*p4(1) + p1(1)*p2(3)*p4(2) + p1(2)*p2(1)*p4(3))
    det_134 = p1(1)*p3(2)*p4(3) + p1(2)*p3(3)*p4(1) + p1(3)*p3(1)*p4(2) - &
         ( p1(3)*p3(2)*p4(1) + p1(1)*p3(3)*p4(2) + p1(2)*p3(1)*p4(3))
    det_234 = p2(1)*p3(2)*p4(3) + p2(2)*p3(3)*p4(1) + p2(3)*p3(1)*p4(2) - &
         ( p2(3)*p3(2)*p4(1) + p2(1)*p3(3)*p4(2) + p2(2)*p3(1)*p4(3))
!
    ! use Laplace's formula on the 4'th row (omit multiplying by 1 here, of course)
    D0 = -det_234 + det_134 - det_124 + det_123
!
    ! check if tet is degenerate or has negative jacobian ( jacobian = -D0, see comment box on bottom of routine)
    if(D0 >= 0.) then
       write(errstr,*) "jacobian = ",-D0,", <= 0, means that invgrid cell ",icell," (4-node tetrahedral) is degenerate "//&
            "(jacobian = 0), or its nodes are indexed in a wrong order (negative jacobian)"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
!
! D2,D3,D4
!
    allocate(det_X12(np),det_X13(np),det_X14(np),det_X23(np),det_X24(np),det_X34(np))
    allocate(D2(np),D3(np),D4(np))
    ! in order to compute D2,..,D4, we only need 6 sub-determinants (each 3x3): det_Xij with x,y,z in first column (which 'X' stands for),
    ! pi in second and pj in third column. Then: 
    ! use Laplace's formula on the 4'th row (omit multiplying by 1 here, of course)
    ! AND USE THE RULE: CHANING TWO COLUMNS OF THE MATRIX SWAPS THE SIGN OF THE DETERMINANT!
    det_X12 = ( x*(p1(2)*p2(3)) + y*(p1(3)*p2(1)) + z*(p1(1)*p2(2)) ) - &
         ( z*(p1(2)*p2(1)) + x*(p1(3)*p2(2)) + y*(p1(1)*p2(3)) )
    det_X13 = ( x*(p1(2)*p3(3)) + y*(p1(3)*p3(1)) + z*(p1(1)*p3(2)) ) - &
         ( z*(p1(2)*p3(1)) + x*(p1(3)*p3(2)) + y*(p1(1)*p3(3)) )
    det_X14 = ( x*(p1(2)*p4(3)) + y*(p1(3)*p4(1)) + z*(p1(1)*p4(2)) ) - &
         ( z*(p1(2)*p4(1)) + x*(p1(3)*p4(2)) + y*(p1(1)*p4(3)) )
    det_X23 = ( x*(p2(2)*p3(3)) + y*(p2(3)*p3(1)) + z*(p2(1)*p3(2)) ) - &
         ( z*(p2(2)*p3(1)) + x*(p2(3)*p3(2)) + y*(p2(1)*p3(3)) )
    det_X24 = ( x*(p2(2)*p4(3)) + y*(p2(3)*p4(1)) + z*(p2(1)*p4(2)) ) - &
         ( z*(p2(2)*p4(1)) + x*(p2(3)*p4(2)) + y*(p2(1)*p4(3)) )
    det_X34 = ( x*(p3(2)*p4(3)) + y*(p3(3)*p4(1)) + z*(p3(1)*p4(2)) ) - &
         ( z*(p3(2)*p4(1)) + x*(p3(3)*p4(2)) + y*(p3(1)*p4(3)) )
!
    !D1 = - det_234 + det_X34 - det_X24 + det_X23
!
    !D2= - det_X34 + det_134 - det_1X4 + det_1X3  rule: chaning two columns of the matrix swaps the sign of the determinant!
    D2 = - det_X34 + det_134 + det_X14 - det_X13 
!
    !D3= - det_2X4 + det_1X4 - det_124 + det_12X  need do swap columns twice in det_12X (to yield det_X12), so sign stays ('+')
    D3 =   det_X24 - det_X14 - det_124 + det_X12
!
    !D4= - det_23X + det_13X - det_12X + det_123
    D4 = - det_X23 + det_X13 - det_X12 + det_123
!
    ! xi = b2 = D2/D0 , see comment box above
    x = D2/D0
    ! eta = b3 = D3/D0 , see comment box above
    y = D3/D0
    ! zeta = b4 = D4/D0 , see comment box above
    z = D4/D0
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! JACOBIAN
!
! having now computed xi,eta,zeta (contained in x,y,z), we still need the jacobian of the transformation
! from standard cell to real coordinates, which is an affine linear transform, hence its jacobian is constant
!
! IN FACT, THE JACOBIAN WAS ALREADY COMPUTED ABOVE: IT EQUALS -D0 , AS EXPLAINED IN THE FOLLOWING:
!
! the transformation is:
! X(xi,eta,zeta) = p1 + (p2-p1)*xi + (p3-p1)*eta + (p4-p1)*zeta
!
! its jacobian matrix is:
! J = (  p2-p1 | p3-p1 | p4-p1  ) , where pj-p1 are column vectors (of length 3)
!
! transforming the matrix related to D0, det(J) can be computed as follows
! (note that the determinant does not change, when adding a multiple
! of some column to some other column):
!
! D0 =  det ( p1 | p2 | p3 | p4 ) = det ( p1 | p2-p1 | p3-p1 | p4-p1 ) =  - det (p2-p1 | p3-p1 | p4-p1) = - det(J)
!           (  1 |  1 |  1 |  1 )       (  1 |   0   |   0   |   0   )
!
! THUS, the jacobian equals -D0
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
    jacobian = -D0
! 
    deallocate(det_X12,det_X13,det_X14,det_X23,det_X24,det_X34,D2,D3,D4)
  end subroutine transformToStandardTetEcartInversionGrid
!------------------------------------------------------------------------
!> \brief transform given coordinates of points contained in hex8 cell icell_hex8 to standard tet and compute their jacobian
!! \param this ecart inversion grid
!! \param icell global cell index of hex8 under consideration
!! \param x vector of global x coordinate (contains x-values in standard cell on exit)
!! \param y vector of global y coordinate (contains y-values in standard cell on exit)
!! \param z vector of global z coordinate (contains z-values in standard cell on exit)
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!
  subroutine transformToStandardHexEcartInversionGrid(this,icell,x,y,z,jacobian,errmsg)
    type (ecart_inversion_grid) :: this
    integer :: icell
    real, dimension(:), intent(inout) :: x,y,z,jacobian
    type (error_message) :: errmsg
    ! local
    character (len=40) :: myname = 'transformToStandardHexEcartInversionGrid'
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
    write(errstr,*) "for 8-node hexahedral invgrid cells (as cell ",icell,"), this operation is not supported yet"
    call add(errmsg,2,errstr,myname)
  end subroutine transformToStandardHexEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which volume should be returned
!! \param volume volume of cell icell
!! \param errmsg error message
!
  subroutine getVolumeCellEcartInversionGrid(this,icell,volume,errmsg)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    type (error_message) :: errmsg
    ! local
    character (len=31) :: myname = 'getVolumeCellEcartInversionGrid'
!
    call addTrace(errmsg,myname)
!
    select case (this%ctype(icell))
    case(4)
       call getVolumeTet4EcartInversionGrid(this,icell,volume,errmsg)
    case(8)
       call getVolumeHex8EcartInversionGrid(this,icell,volume,errmsg)
    end select

  end subroutine getVolumeCellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of tet4 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which volume should be returned
!! \param volume volume of tet4 cell icell
!! \param errmsg error message
!
  subroutine getVolumeTet4EcartInversionGrid(this,icell,volume,errmsg)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    type (error_message) :: errmsg
    ! local
    character (len=31) :: myname = 'getVolumeTet4EcartInversionGrid'
    character(len=400) :: errstr
    real, dimension(3) :: p1,p2,p3,p4
!
    call addTrace(errmsg,myname)
    ! in routine getVolumeCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_tet4(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_tet4(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_tet4(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_tet4(4,this%icell_type(icell)))
!
    ! volume is a six'th of volume of a corresponding parallelepiped which computes as the triple 
    ! product of the vectors (p2-p1),(p3-p1),(p4-p1)
    volume = ( (p2(1)-p1(1))*(p3(2)-p1(2))*(p4(3)-p1(3)) + &
               (p2(2)-p1(2))*(p3(3)-p1(3))*(p4(1)-p1(1)) + &
               (p2(3)-p1(3))*(p3(1)-p1(1))*(p4(2)-p1(2)) ) &
           - ( (p2(3)-p1(3))*(p3(2)-p1(2))*(p4(1)-p1(1)) + &
               (p2(1)-p1(1))*(p3(3)-p1(3))*(p4(2)-p1(2)) + &
               (p2(2)-p1(2))*(p3(1)-p1(1))*(p4(3)-p1(3)) )
    volume = volume*(1./6.)
!
    ! if the edges (p2-p1),(p3-p1),(p4-p1) do not form a right-handed system (i.e. the triple product is negative),
    ! then the nodes of the tetrahedron are indexed in a wrong way
    if(volume <= 0) then
       write(errstr,*) "cell volume = ",volume,", <= 0, means that invgrid cell ",icell," (4-node tetrahedral) is degenerate "//&
            "(volume = 0), or its nodes are indexed in a wrong order (negative jacobian)"
       call add(errmsg,2,errstr,myname)
       return
    end if
  end subroutine getVolumeTet4EcartInversionGrid
!------------------------------------------------------------------------
!> \brief get volume of hex8 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which volume should be returned
!! \param volume volume of hex8 cell icell
!! \param errmsg error message
!
  subroutine getVolumeHex8EcartInversionGrid(this,icell,volume,errmsg)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: volume
    type (error_message) :: errmsg
    ! local
    character (len=31) :: myname = 'getVolumeHex8EcartInversionGrid'
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
    write(errstr,*) "for 8-node hexahedral invgrid cells (as cell ",icell,"), volume computation is not supported yet"
    call add(errmsg,2,errstr,myname)
  end subroutine getVolumeHex8EcartInversionGrid
!------------------------------------------------------------------------
!> \brief get center of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which center should be returned
!! \param xc first coordinate of center of cell icell
!! \param yc second coordinate of center of cell icell
!! \param zc third coordinate of center of cell icell
!
  subroutine getCenterCellEcartInversionGrid(this,icell,xc,yc,zc)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: xc,yc,zc
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! no need of selecting coords_type: 
    ! always do the same, since in case of using the 
    ! ecart inversion grid the wavefield points as well as event and station coordinates
    ! are expected as x,y,z coordinates (in that order)
!
    select case (this%ctype(icell))
    case(4)
       call getCenterTet4EcartInversionGrid(this,icell,xc,yc,zc)
    case(8)
       call getCenterHex8EcartInversionGrid(this,icell,xc,yc,zc)
    end select
  end subroutine getCenterCellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get center of tet4 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which center should be returned
!! \param xc first coordinate of center of cell icell
!! \param yc second coordinate of center of cell icell
!! \param zc third coordinate of center of cell icell
!
  subroutine getCenterTet4EcartInversionGrid(this,icell,xc,yc,zc)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: xc,yc,zc
    ! local
    real, dimension(3) :: p1,p2,p3,p4
!
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_tet4(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_tet4(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_tet4(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_tet4(4,this%icell_type(icell)))
!
    ! return as the center the barycenter
    xc = 0.25*(p1(1)+p2(1)+p3(1)+p4(1))
    yc = 0.25*(p1(2)+p2(2)+p3(2)+p4(2))
    zc = 0.25*(p1(3)+p2(3)+p3(3)+p4(3))
  end subroutine getCenterTet4EcartInversionGrid
!------------------------------------------------------------------------
!> \brief get center of hex8 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which center should be returned
!! \param xc first coordinate of center of cell icell
!! \param yc second coordinate of center of cell icell
!! \param zc third coordinate of center of cell icell
!
  subroutine getCenterHex8EcartInversionGrid(this,icell,xc,yc,zc)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: xc,yc,zc
    ! local
    real, dimension(3) :: p1,p2,p3,p4,p5,p6,p7,p8
!
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_hex8(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_hex8(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_hex8(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_hex8(4,this%icell_type(icell)))
    p5 = this%point(:,this%cell_hex8(5,this%icell_type(icell)))
    p6 = this%point(:,this%cell_hex8(6,this%icell_type(icell)))
    p7 = this%point(:,this%cell_hex8(7,this%icell_type(icell)))
    p8 = this%point(:,this%cell_hex8(8,this%icell_type(icell)))
!
    ! return as the center the "barycenter"
    xc = 0.125*(p1(1)+p2(1)+p3(1)+p4(1)+p5(1)+p6(1)+p7(1)+p8(1))
    yc = 0.125*(p1(2)+p2(2)+p3(2)+p4(2)+p5(2)+p6(2)+p7(2)+p8(2))
    zc = 0.125*(p1(3)+p2(3)+p3(3)+p4(3)+p5(3)+p6(3)+p7(3)+p8(3))
  end subroutine getCenterHex8EcartInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of inversion grid cell
!! \param this inversion grid
!! \param icell index of inversion grid for which radius should be returned
!! \param radius radius of cell icell
!
  subroutine getRadiusCellEcartInversionGrid(this,icell,radius)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
!
    select case (this%ctype(icell))
    case(4)
       call getRadiusTet4EcartInversionGrid(this,icell,radius)
    case(8)
       call getRadiusHex8EcartInversionGrid(this,icell,radius)
    end select
  end subroutine getRadiusCellEcartInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of tet4 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which radius should be returned
!! \param radius radius of tet4 cell icell
!
  subroutine getRadiusTet4EcartInversionGrid(this,icell,radius)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
    ! local
    real, dimension(3) :: p1,p2,p3,p4
    real :: xc,yc,zc
!
    ! in routine getRadiusCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_tet4(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_tet4(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_tet4(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_tet4(4,this%icell_type(icell)))
!
    ! compute the barycenter
    xc = 0.25*(p1(1)+p2(1)+p3(1)+p4(1))
    yc = 0.25*(p1(2)+p2(2)+p3(2)+p4(2))
    zc = 0.25*(p1(3)+p2(3)+p3(3)+p4(3))
!
    ! return as radius the maximum distance of a corner point from the barycenter
    ! IN THE FUTURE: maybe radius of circumball is better suited here
    radius = sqrt( (p1(1)-xc)*(p1(1)-xc) + (p1(2)-yc)*(p1(2)-yc) + (p1(3)-zc)*(p1(3)-zc) )
    radius = max(radius, sqrt( (p2(1)-xc)*(p2(1)-xc) + (p2(2)-yc)*(p2(2)-yc) + (p2(3)-zc)*(p2(3)-zc) ) )
    radius = max(radius, sqrt( (p3(1)-xc)*(p3(1)-xc) + (p3(2)-yc)*(p3(2)-yc) + (p3(3)-zc)*(p3(3)-zc) ) )
    radius = max(radius, sqrt( (p4(1)-xc)*(p4(1)-xc) + (p4(2)-yc)*(p4(2)-yc) + (p4(3)-zc)*(p4(3)-zc) ) )
  end subroutine getRadiusTet4EcartInversionGrid
!------------------------------------------------------------------------
!> \brief get radius of hex8 inversion grid cell
!! \param this inversion grid
!! \param icell global index of inversion grid cell for which radius should be returned
!! \param radius radius of hex8 cell icell
!
  subroutine getRadiusHex8EcartInversionGrid(this,icell,radius)
    type (ecart_inversion_grid) :: this
    integer, intent(in) :: icell
    real :: radius
    ! local
    real, dimension(3) :: p1,p2,p3,p4,p5,p6,p7,p8
    real :: xc,yc,zc
!
    ! in routine getCenterCellInversionGrid of module inversionGrid it was already assured
    ! icell is valid
!
    ! get all nodes of this thetrahedron
    p1 = this%point(:,this%cell_hex8(1,this%icell_type(icell)))
    p2 = this%point(:,this%cell_hex8(2,this%icell_type(icell)))
    p3 = this%point(:,this%cell_hex8(3,this%icell_type(icell)))
    p4 = this%point(:,this%cell_hex8(4,this%icell_type(icell)))
    p5 = this%point(:,this%cell_hex8(5,this%icell_type(icell)))
    p6 = this%point(:,this%cell_hex8(6,this%icell_type(icell)))
    p7 = this%point(:,this%cell_hex8(7,this%icell_type(icell)))
    p8 = this%point(:,this%cell_hex8(8,this%icell_type(icell)))
!
    ! compute the "barycenter"
    xc = 0.125*(p1(1)+p2(1)+p3(1)+p4(1)+p5(1)+p6(1)+p7(1)+p8(1))
    yc = 0.125*(p1(2)+p2(2)+p3(2)+p4(2)+p5(2)+p6(2)+p7(2)+p8(2))
    zc = 0.125*(p1(3)+p2(3)+p3(3)+p4(3)+p5(3)+p6(3)+p7(3)+p8(3))
!
    ! return as radius the maximum distance of a corner point from the "barycenter"
    radius = sqrt( (p1(1)-xc)*(p1(1)-xc) + (p1(2)-yc)*(p1(2)-yc) + (p1(3)-zc)*(p1(3)-zc) )
    radius = max(radius, sqrt( (p2(1)-xc)*(p2(1)-xc) + (p2(2)-yc)*(p2(2)-yc) + (p2(3)-zc)*(p2(3)-zc) ) )
    radius = max(radius, sqrt( (p3(1)-xc)*(p3(1)-xc) + (p3(2)-yc)*(p3(2)-yc) + (p3(3)-zc)*(p3(3)-zc) ) )
    radius = max(radius, sqrt( (p4(1)-xc)*(p4(1)-xc) + (p4(2)-yc)*(p4(2)-yc) + (p4(3)-zc)*(p4(3)-zc) ) )
    radius = max(radius, sqrt( (p5(1)-xc)*(p5(1)-xc) + (p5(2)-yc)*(p5(2)-yc) + (p5(3)-zc)*(p5(3)-zc) ) )
    radius = max(radius, sqrt( (p6(1)-xc)*(p6(1)-xc) + (p6(2)-yc)*(p6(2)-yc) + (p6(3)-zc)*(p6(3)-zc) ) )
    radius = max(radius, sqrt( (p7(1)-xc)*(p7(1)-xc) + (p7(2)-yc)*(p7(2)-yc) + (p7(3)-zc)*(p7(3)-zc) ) )
    radius = max(radius, sqrt( (p8(1)-xc)*(p8(1)-xc) + (p8(2)-yc)*(p8(2)-yc) + (p8(3)-zc)*(p8(3)-zc) ) )
!
  end subroutine getRadiusHex8EcartInversionGrid
!
end module ecartInversionGrid
