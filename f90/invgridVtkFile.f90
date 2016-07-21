!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief module to write data living on invgrids to vkt output
!!
!! \details The inversion grid geometry (as type USTRUCTURED_GRID) 
!!  alongside with scalar data living on the inversion grid is written
!!  to binary or ascii vtk files. Complex data is handled as 2 component 
!!  scalar float data. As an option, multiple files containing data
!!  w.r.t. some index (frequency, time) may be written having the same
!!  file base name followed by an index, in order to be considered by 
!!  Paraview as a sequence of data.
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module invgridVtkFile
!
  use inversionGrid
  use errorMessage
!
  implicit none
!
  ! by allowing calls to only certain routines, it is tried to assure that 
  ! a file is open when writing to it, etc. (this way, certain security checks
  ! may be omitted, as unintended calls to auxilliary routines are forbidden)
  private
  public :: invgrid_vtk_file,init,writeData,writeInvgrid,dealloc
!
  interface init; module procedure initiateInvgridVtkFile; end interface
  interface writeData
     module procedure writeRealDataInvgridVtkFile
     !module procedure writeIntegerDataInvgridVtkFile
     module procedure writeComplexDataInvgridVtkFile
  end interface writeData
  interface writeInvgrid
     module procedure writeInvgridVtkFile
  end interface writeInvgrid
  interface dealloc; module procedure deallocateInvgridVtkFile; end interface
!
!> \brief general file, geometry and cell information of vtk file
  type invgrid_vtk_file
     private
     ! general file information
     character (len=300) :: filename = '' !< (base) file name of vtk file (WITHOUT '.vtk')
     logical is_ascii !< indicating ascii (true) or binary (false) format of vtk file
     character (len=200) :: title = '' !< second line of vtk file
     !
     ! other general information
     integer :: geometry_type = -1 !< type of vtk geometry:  0 = volumetric cells , 1 = cell center points
     integer :: ntot_invgrid = 0 !< total number of cells in inversion grid
     integer :: ndata = 0 !< number of data handled in this vtk file (i.e. number of vtk cells or cell center points)
     integer, dimension(:), pointer :: invgrid_cell_indx => null() !< mapping of size(ndata): for each datum i, invgrid_cell_indx(i) is the corresponding invgrid cell index
     !
     ! vtk points
     integer :: npoints = 0 !< number of points in vtk file
     real, dimension(:,:), pointer :: points => null() !< POINTS geometry for UNSTRUCTURED_GRID
     !
     ! vtk cells
     integer :: ncells = 0 !< number of inversion grid cells handled in this vtk file
     integer, dimension(:), pointer :: cell_connectivity => null() !< array of point indices defining the vtk cells
     integer, dimension(:), pointer :: cell_type => null() !< vtk cell type for each cell (usually all equal)
     !
     ! we need the additional index mapping "req_indx" here for the case of present(cell_indx_req), as routine 
     ! getGeometryVtkInversionGrid may throw away invalid and duplicate invgrid cell indices. in order to be able to 
     ! still use the cell geometry information with data vectors of same length (and order) as cell_indx_req, req_indx
     ! maps the vtk cell index of the returned vtk cells to the original position in array cell_indx_req.
     ! also, if routine getGeometryVtkInversionGrid does not preserve the order of vtk cells as requested by cell_indx_req, 
     ! map req_indx will reconstruct this order (and be the identity otherwise)
     integer, dimension(:), pointer :: req_indx => null() !< if present(cell_indx_req), maps invgrid_cell_indx to cell_indx_req, otherwise identity
  end type invgrid_vtk_file
!
contains
!------------------------------------------------------------------------
!> \brief initiate naming and geometry structure of vtk file
!! \details Define base filename, format (ascii /binary) and header title of vtk file
!!  and get the cell geometry from module inversionGrid.
!! \param this invgrid vtk file
!! \param invgrid inversion grid
!! \param filename vtk base name (without '.vtk'). Will optionally be concatenated with a filename extension
!! \param vtk_format 'ASCII' or 'BINARY' indicating the vtk file format
!! \param vtk_title optional second line of vtk file (by default 'data on inversion grid')
!! \param errmsg error message
!! \return error message
!
  subroutine initiateInvgridVtkFile(this,invgrid,filename,vtk_format,errmsg,vtk_title,cell_indx_req)
    ! incoming
    type (invgrid_vtk_file) :: this
    type (inversion_grid) :: invgrid
    character(len=*) :: filename,vtk_format
    character(len=*), optional :: vtk_title
    integer, dimension(:), optional :: cell_indx_req
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=18) :: myname = 'initInvgridVtkFile'
    character (len=400) :: errstr
    integer :: num_cell_indx_req
    !
    call addTrace(errmsg,myname)
    !
    if(this%filename/= '') then
       call add(errmsg,1,"this object was already initiated. deallocating it now, before initiating anew",myname)
       call deallocateInvgridVtkFile(this)
    endif
    !
    if(trim(vtk_format) /= 'ASCII' .and. trim(vtk_format) /= 'BINARY') then
       errstr = "incoming vtk_format string '"//trim(vtk_format)//"' not valid, must be either 'ASCII' or 'BINARY'"
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    this%is_ascii = (trim(vtk_format) == 'ASCII')
    if(trim(filename) == '') then
       call add(errmsg,2,'incoming filename is empty string',myname)
       return
    else
       this%filename = trim(filename)
    endif
    if(present(cell_indx_req)) then 
       num_cell_indx_req = size(cell_indx_req)
       if(num_cell_indx_req .le. 0) then
          write(errstr,*) "number of incoming requested invgrid cells (",num_cell_indx_req,&
               ") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
    end if
    nullify(this%points,this%cell_connectivity,this%cell_type,this%invgrid_cell_indx,this%req_indx)
    this%geometry_type = -1
    call getGeometryVtkInversionGrid(invgrid,this%geometry_type,this%points,this%cell_connectivity,&
         this%cell_type,this%invgrid_cell_indx,errmsg,cell_indx_req,this%req_indx)
    if(.level.errmsg == 2) return
    
    select case (this%geometry_type)
    case(0) ! VOLUMETRIC CELLS
       if(.not.(associated(this%points).and.associated(this%cell_connectivity).and.associated(this%cell_type).and.&
            associated(this%invgrid_cell_indx))) then
          call add(errmsg,2,"although a volumetric cell geometry should be created, not all required components "//&
               "are returned by routine getGeometryVtkInversionGrid",myname)
       end if
       this%npoints = size(this%points,2)
       this%ncells = size(this%cell_type)
       this%ndata = this%ncells
       if(present(vtk_title)) then
          this%title = trim(vtk_title)
       else
          this%title = 'data on inversion grid cells'
       endif
    case(1) ! CELL CENTER POINTS
       if(.not.(associated(this%points).and.associated(this%invgrid_cell_indx))) then
          call add(errmsg,2,"although a point geometry (cell centers) should be created, not all required components "//&
               "are returned by routine getGeometryVtkInversionGrid",myname)
       end if
       this%npoints = size(this%points,2)
       this%ncells = 0
       this%ndata = this%npoints
       if(present(vtk_title)) then
          this%title = trim(vtk_title)
       else
          this%title = 'data on inversion grid cell center points'
       endif
    case default
       write(errstr,*) "geometry type ",this%geometry_type," returned by routine getGeometryVtkInversionGrid is not ",&
            "supported. supported types are: 0 (cells), 1 (cell center points)"
       call add(errmsg,2,trim(errstr),myname)
       return
    end select

    if(present(cell_indx_req)) then
       if(this%ndata/=num_cell_indx_req) then
          write(errstr,*) "num of vtk cells / cell centers returned by getGeometryVtkInversionGrid (",this%ndata,&
               ") not equal to num of invgrid cells originally requested to be used (",num_cell_indx_req,&
               ") => there were duplicate or invalid indices in requested index array cell_indx_req, "//&
               "which may indicate inconsistencies of files"
          call add(errmsg,2,trim(errstr),myname)
          return
       endif
    endif
    this%ntot_invgrid = .ncell.invgrid
  end subroutine initiateInvgridVtkFile
!------------------------------------------------------------------------
!> \brief open vtk file to write
!! \param this invgrid vtk file
!! \param lu file unit
!! \param file_index optional index of file (will be appended to filename base)
!! \param overwrite_in optional logical to indicate behaviour in case file exists, by default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine openInvgridVtkFile(this,lu,filename_extension,errmsg,overwrite_in)
    ! incoming
    type (invgrid_vtk_file) :: this
    integer :: lu
    character (len=*) :: filename_extension
    logical, optional :: overwrite_in
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=18) :: myname = 'openInvgridVtkFile'
    character (len=400) :: vtk_file
    integer :: ios
    logical :: overwrite,file_exists
    character (len=7) :: open_status
    !
    call addTrace(errmsg,myname)
    !
    if(present(overwrite_in)) then
       overwrite = overwrite_in
    else
       ! by default, do not overwrite
       overwrite = .False.
    endif
    !
    ! create filename from basename and (possibly empty) filename extension plus '.vtk' extension
    vtk_file = trim(this%filename)//trim(filename_extension)//'.vtk'
    !
    ! check if file exists, set open_status as required
    inquire(file=trim(vtk_file),exist=file_exists)
    if(file_exists) then
       if(overwrite) then
          ! give warning and replace existing file
          call add(errmsg,1,"file '"//trim(vtk_file)//"' exists and will be overwritten",myname)
          open_status = 'REPLACE'
       else ! overwrite
          ! raise error
          call add(errmsg,2,"file '"//trim(vtk_file)//"' exists. Please rename or indicate 'overwrite'",myname)
          return
       endif ! overwrite
    else ! file_exists
       open_status = 'NEW'
    endif ! file_exists
    !
    ! open file. According to value of open_status, an existing file is overwritten
    if(this%is_ascii) then
       open(unit=lu,file=trim(vtk_file),form='FORMATTED',status=trim(open_status),action='WRITE',iostat=ios)
       if(ios/=0) call add(errmsg,2,"could not open ascii file '"//trim(vtk_file)//"'to write",myname)
    else ! this%is_ascii
       open(unit=lu,file=trim(vtk_file),form='UNFORMATTED',access='STREAM',status=trim(open_status),action='WRITE',&
            convert='BIG_ENDIAN',iostat=ios)
       if(ios/=0) call add(errmsg,2,"could not open binary file '"//trim(vtk_file)//"'to write",myname)
    endif ! this%is_ascii
  end subroutine openInvgridVtkFile
!------------------------------------------------------------------------
!> \brief write vtk header, points and cell structure to open vtk file
!! \param this invgrid vtk file
!! \param lu file unit of file (MUST BE OPENED ALREADY!)
!! \param errmsg error message
!! \return error message
!
  subroutine writeHeaderGeometryInvgridVtkFile(this,lu,errmsg)
    type (invgrid_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    ! local
    integer :: ios
    character(len=33) :: myname = 'writeHeaderGeometryInvgridVtkFile'
    character(len=25) :: vtk_dataset_type
    character (len=1), parameter :: eol_char = char(10)
    logical :: err
!
    call addTrace(errmsg,myname)
!
    select case(this%geometry_type)
    case(0) ! VOLUMETRIC CELLS
       vtk_dataset_type = 'DATASET UNSTRUCTURED_GRID'
    case(1) ! CELL CENTER POINTS
       vtk_dataset_type = 'DATASET POLYDATA'
    end select
!
    ! remember with err if there was an error somewhere
    err = .false.
    if(this%is_ascii) then
       ! vkt Header
       write(unit=lu,fmt='(a)',iostat=ios) '# vtk DataFile Version 3.1'  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) trim(this%title)              ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'ASCII'                       ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) trim(vtk_dataset_type)        ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing vtk Header',myname)
          return
       endif
    else ! this%is_ascii
       ! vtk Header
       write(unit=lu,iostat=ios) '# vtk DataFile Version 3.1'//eol_char  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) trim(this%title)//eol_char              ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'BINARY'//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) trim(vtk_dataset_type)//eol_char        ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing vtk Header',myname)
          return
       endif
    endif ! this%is_ascii

    select case(this%geometry_type)
    case(0) ! VOLUMETRIC CELLS
       call writePointsGeometryInvgridVtkFile(this,lu,errmsg,myname)
       call writeCellsGeometryInvgridVtkFile(this,lu,errmsg,myname)
    case(1) ! CELL CENTER POINTS
       call writePointsGeometryInvgridVtkFile(this,lu,errmsg,myname)
       call writeVerticesGeometryInvgridVtkFile(this,lu,errmsg,myname)
    end select
  end subroutine writeHeaderGeometryInvgridVtkFile
!------------------------------------------------------------------------
  subroutine writePointsGeometryInvgridVtkFile(this,lu,errmsg,myname)
    type (invgrid_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    character(len=*) :: myname
    ! logical
    integer :: ios
    character (len=500) :: string
    logical :: err
    character (len=1), parameter :: eol_char = char(10)
!
    err = .false.
    if(this%is_ascii) then
       ! POINTS
       write(unit=lu,fmt='(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(3e14.6e2)',iostat=ios) this%points                     ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
          return
       endif
    else ! this%is_ascii
       ! POINTS
       write(string,'(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'
       write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%points                             ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
          return
       endif
    end if ! this%is_ascii
  end subroutine writePointsGeometryInvgridVtkFile
!------------------------------------------------------------------------
  subroutine writeVerticesGeometryInvgridVtkFile(this,lu,errmsg,myname)
    type (invgrid_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    character(len=*) :: myname
    ! logical
    integer :: ios,i
    character (len=500) :: string
    logical :: err
    character (len=1), parameter :: eol_char = char(10)
!
    err = .false.
    if(this%is_ascii) then
       ! VERTICES
       write(unit=lu,fmt='(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints     ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(2i12)',iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)            ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    else ! this%is_ascii
       ! VERTICES
       write(string,'(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)          ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    end if ! this%is_ascii
  end subroutine writeVerticesGeometryInvgridVtkFile
!------------------------------------------------------------------------
  subroutine writeCellsGeometryInvgridVtkFile(this,lu,errmsg,myname)
    type (invgrid_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    character(len=*) :: myname
    ! logical
    integer :: ios
    character (len=500) :: string
    logical :: err
    character (len=1), parameter :: eol_char = char(10)
!
    err = .false.
    if(this%is_ascii) then
       ! CELL CONNECTIVITY
       write(unit=lu,fmt='(a,2i12)',iostat=ios)'CELLS ',this%ncells,size(this%cell_connectivity)  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(i12)',iostat=ios) this%cell_connectivity                               ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELLS',myname)
          return
       endif
       ! CELL TYPES
       write(unit=lu,fmt='(a,i12)',iostat=ios) 'CELL_TYPES ',this%ncells  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(i12)',iostat=ios) this%cell_type               ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELL_TYPES',myname)
          return
       endif
    else ! this%is_ascii
       ! CELL CONNECTIVITY
       write(string,'(a,2i12)',iostat=ios)'CELLS ',this%ncells,size(this%cell_connectivity)
       write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%cell_connectivity                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELLS',myname)
          return
       endif
       ! CELL TYPES
       write(string,'(a,i12)',iostat=ios)'CELL_TYPES ',this%ncells
       write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%cell_type                          ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELL_TYPES',myname)
          return
       endif
    end if ! this%is_ascii
  end subroutine writeCellsGeometryInvgridVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and one component float scalar valued cell data to vtk file
!! \details The number of incoming real data values must match the number this%ndata
!!  of this invgrid_vtk_file object, because here scalar CELL_DATA (i.e. one scalar value per cell in case 
!!  this%geometry_type == 0) or scalar POINT_DATA (one value per inversion grid cell center point in case 
!!  this%geometry_type == 1) is added to the vtk file. 
!!  First a file is opened calling openInvgridVtkFile and header and point and cell 
!!  geometry information is written to that file calling writeHeaderGeometryInvgridVtkFile.
!!  Then the incoming data values are added to the vtk file as scalar valued float cell data.
!!  If data_indx_is_invgrid is present and true, incoming values in data are expected for every invgrid cell (in 
!!  invgrid cell order), otherwise incoming values in data are expected for every cell index as in cell_indx_req
!!  when initiating this. Id cell_indx_req was not present when initiating this, the requested indices are
!!  simply all invgrid cell (in invgrid cell order)
!! \param this invgrid vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param data_indx_is_invgrid optional logical indicating if there are incoming data for ALL invgrid cells or not
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param fname_extension optional character string as file name extension IN FRONT OF file index: this%filename//fname_extension//file_indx
!! \param errmsg error message
!! \return error message
!
  subroutine writeRealDataInvgridVtkFile(this,lu,data,errmsg,data_name,file_index,data_indx_is_invgrid,overwrite,&
       fname_extension)
    ! incoming
    type (invgrid_vtk_file) :: this
    integer :: lu
    real, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: data_indx_is_invgrid,overwrite
    character (len=*), optional :: fname_extension
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=27) :: myname = 'writeRealDataInvgridVtkFile'
    character(len=10) :: vtk_data_type
    character (len=500) :: string
    character (len=500) :: filename_extension
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata_expect
    logical :: err,map_indx_all
!
    call addTrace(errmsg,myname)
!
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this invgrid_vtk_file object is not initiated yet',myname)
       return
    endif
!
    select case(this%geometry_type)
    case(0) ! VOLUMETRIC CELLS
       ! set variables for the vtk file line defining the type of data
       vtk_data_type = 'CELL_DATA'
    case(1) ! CELL CENTER POINTS
       ! set variables for the vtk file line defining the type of data
       vtk_data_type = 'POINT_DATA'
    end select
!
    ! by default, expect as many data as cells are handled in this vtk file
    ndata_expect = this%ndata
    map_indx_all = .false.
!
    ! check if optionally data for ALL cells are to be maped to those defined in this vtk file (can be substet)
    if(present(data_indx_is_invgrid)) then
       if(data_indx_is_invgrid) then
          ndata_expect = this%ntot_invgrid
          map_indx_all = .true.
       end if
    end if
!
    ! check size of incoming data
    if(size(data) /= ndata_expect) then
       write(errstr,*) 'number of incoming data ( =',size(data),&
            ') does not match number of expected values ( =',ndata_expect,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
!
    filename_extension = ''
    if(present(fname_extension)) then
       filename_extension = trim(filename_extension)//trim(fname_extension)
    end if
    if(present(file_index)) then
       write(string,"('_',i6.6)") file_index
       filename_extension = trim(filename_extension)//trim(string)
    endif
    ! open vtk file to write
    call openInvgridVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
!
    ! write header and geometry information to file
    call writeHeaderGeometryInvgridVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
!
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata    ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,fmt='(e14.6e2)', iostat=ios) data(this%invgrid_cell_indx) ; err = err.or.(ios/=0)
       else
          write(unit=lu,fmt='(e14.6e2)', iostat=ios) data(this%req_indx)          ; err = err.or.(ios/=0)
       end if
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing '//trim(vtk_data_type),myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,iostat=ios) data(this%invgrid_cell_indx)             ; err = err.or.(ios/=0)
       else
          write(unit=lu,iostat=ios) data(this%req_indx)                      ; err = err.or.(ios/=0)
       end if
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing '//trim(vtk_data_type),myname)
          close(lu)
          return
       endif
    endif
    !
    ! close file
    close(lu)
  end subroutine writeRealDataInvgridVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and two component float scalar valued cell data to vtk file
!! \details The number of incoming complex data values must match the number this%ndata
!!  of this invgrid_vtk_file object, because here two component scalar CELL_DATA (i.e. one complex scalar value per 
!!  cell in case this%geometry_type == 0) or two component scalar POINT_DATA (one complex value per inversion grid 
!!  cell center point in case this%geometry_type == 1) is added to the vtk file. 
!!  First a file is opened calling openInvgridVtkFile and header and point and cell 
!!  geometry information is written to that file calling writeHeaderGeometryInvgridVtkFile.
!!  Then the incoming data values are added to the vtk file as two component scalar valued float cell data.
!!  If data_indx_is_invgrid is present and true, incoming values in data are expected for every invgrid cell (in 
!!  invgrid cell order), otherwise incoming values in data are expected for every cell index as in cell_indx_req
!!  when initiating this. Id cell_indx_req was not present when initiating this, the requested indices are
!!  simply all invgrid cell (in invgrid cell order)
!! \param this invgrid vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param data_indx_is_invgrid optional logical indicating if there are incoming data for ALL invgrid cells or not
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param fname_extension optional character string as file name extension IN FRONT OF file index: this%filename//fname_extension//file_indx
!! \param errmsg error message
!! \return error message
!
  subroutine writeComplexDataInvgridVtkFile(this,lu,data,errmsg,data_name,file_index,data_indx_is_invgrid,overwrite,&
       fname_extension)
    ! incoming
    type (invgrid_vtk_file) :: this
    integer :: lu
    complex, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: data_indx_is_invgrid,overwrite
    character (len=*), optional :: fname_extension
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=30) :: myname = 'writeComplexDataInvgridVtkFile'
    character(len=10) :: vtk_data_type
    character (len=500) :: string
    character (len=500) :: filename_extension
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata_expect
    logical :: err,map_indx_all
!
    call addTrace(errmsg,myname)
!
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this invgrid_vtk_file object is not initiated yet',myname)
       return
    endif
!
    select case(this%geometry_type)
    case(0) ! VOLUMETRIC CELLS
       ! set variables for the vtk file line defining the type of data
       vtk_data_type = 'CELL_DATA'
    case(1) ! CELL CENTER POINTS
       ! set variables for the vtk file line defining the type of data
       vtk_data_type = 'POINT_DATA'
    end select

    ! by default, expect as many data as cells are handled in this vtk file
    ndata_expect = this%ndata
    map_indx_all = .false.
!
    ! check if optionally data for ALL cells are to be maped to those defined in this vtk file (can be substet)
    if(present(data_indx_is_invgrid)) then
       if(data_indx_is_invgrid) then
          ndata_expect = this%ntot_invgrid
          map_indx_all = .true.
       end if
    end if
!
    ! check size of incoming data
    if(size(data) /= ndata_expect) then
       write(errstr,*) 'number of incoming data ( =',size(data),&
            ') does not match number of expected values ( =',ndata_expect,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
!
    filename_extension = ''
    if(present(fname_extension)) then
       filename_extension = trim(filename_extension)//trim(fname_extension)
    end if
    if(present(file_index)) then
       write(string,"('_',i6.6)") file_index
       filename_extension = trim(filename_extension)//trim(string)
    endif
    ! open vtk file to write
    call openInvgridVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and geometry information to file
    call writeHeaderGeometryInvgridVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
!
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata         ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 2'
       else
          string = 'SCALARS data float 2'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                            ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'                  ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,fmt='(2e14.6e2)', iostat=ios) data(this%invgrid_cell_indx) ; err = err.or.(ios/=0)
       else
          write(unit=lu,fmt='(2e14.6e2)', iostat=ios) data(this%req_indx)          ; err = err.or.(ios/=0)
       end if
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELL_DATA',myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) trim(vtk_data_type)//' ',this%ndata
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 2'
       else
          string = 'SCALARS data float 2'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,iostat=ios) data(this%invgrid_cell_indx)             ; err = err.or.(ios/=0)
       else
          write(unit=lu,iostat=ios) data(this%req_indx)                      ; err = err.or.(ios/=0)
       end if
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing CELL_DATA',myname)
          close(lu)
          return
       endif
    endif
    !
    ! close file
    close(lu)
  end subroutine writeComplexDataInvgridVtkFile
!------------------------------------------------------------------------
!> \brief open file, write only header and geometry without data to vtk file
!! \param this invgrid vtk file
!! \param lu file unit
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeInvgridVtkFile(this,lu,errmsg,overwrite)
    ! incoming
    type (invgrid_vtk_file) :: this
    integer :: lu
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=19) :: myname = 'writeInvgridVtkFile'
    character (len=1), parameter :: eol_char = char(10)
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this invgrid_vtk_file object is not initiated yet',myname)
       return
    endif
    !
    ! open vtk file to write
    call openInvgridVtkFile(this,lu,'_invgrid',errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and geometry information to file
    call writeHeaderGeometryInvgridVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! close file
    close(lu)
  end subroutine writeInvgridVtkFile
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this invgrid_vtk_file object
!
  subroutine deallocateInvgridVtkFile(this)
    type (invgrid_vtk_file) :: this
    this%filename = ''
    this%title = ''
    this%geometry_type = -1
    this%ntot_invgrid = 0
    this%ndata = 0
    if(associated(this%points)) deallocate(this%points)
    this%npoints = 0
    if(associated(this%cell_connectivity)) deallocate(this%cell_connectivity)
    if(associated(this%cell_type)) deallocate(this%cell_type)
    this%ncells = 0
    if(associated(this%invgrid_cell_indx)) deallocate(this%invgrid_cell_indx)
    if(associated(this%req_indx)) deallocate(this%req_indx)
  end subroutine deallocateInvgridVtkFile
!
end module invgridVtkFile
