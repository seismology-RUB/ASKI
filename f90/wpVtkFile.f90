!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief module to write data on wavefield points to vkt output
!!
!! \details A vtk POLYDATA dataset consisting of wavefield points as VERTICES
!!  along with any scalar POINT DATA living on wavefield points is written
!!  to binary or ascii vtk files. Complex data is handled as 2 component 
!!  scalar float data. As an option, multiple files containing data
!!  w.r.t. some index (frequency, time) may be written having the same
!!  file base name followed by an index, in order to be considered by 
!!  Paraview as a sequence of data.
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module wpVtkFile
!
  use wavefieldPoints
  use inversionGrid
  use errorMessage
!
  implicit none
!
  ! by allowing calls to only certain routines, it is tried to assure that e.g. 
  ! a file is open when writing to it, etc. (this way, certain security checks
  ! may be omitted, as unintended calls to auxilliary routines are forbidden)
  private
  public :: wp_vtk_file,init,writeData,writeWp,dealloc
!
  interface init; module procedure initiateWpVtkFile; end interface
  interface writeData
     module procedure writeRealDataWpVtkFile
     !module procedure writeIntegerDataWpVtkFile
     module procedure writeComplexDataWpVtkFile
  end interface writeData
  interface writeWp; module procedure writeWpVtkFile; end interface
  interface dealloc; module procedure deallocateWpVtkFile; end interface
!
!> \brief general file, geometry and point information of vtk file
  type wp_vtk_file
     private
     character (len=300) :: filename = '' !< (base) file name of vtk file (WITHOUT '.vtk')
     logical :: is_ascii !< indicating ascii (true) or binary (false) format of vtk file
     character (len=100) :: title = '' !< second line of vtk file
     integer :: npoints !< number of vtk points (can be a only a part of all wavefield points)
     real, dimension(:,:), pointer :: points => null() !< POINTS geometry for POLYDATA VERTICES
     integer, dimension(:), pointer :: wp_indx => null() !< for each vtk point, the corresponding wavefield point index
     integer :: ntot_wp !< total number of wavefield points
     !
     ! we need the additional index mapping "req_indx" here for the case of present(wp_indx_req), as routine 
     ! getVtkWavefieldPoints may throw away invalid and duplicate wavefield point indices. in order to be able to 
     ! still use the points geometry information with data vectors of same length (and order) as wp_indx_req, req_indx
     ! maps the vtk point index of the selected vtk points to the original position in array wp_indx_req
     ! also, if routine getVtkWavefieldPoints would not presever the order of vtk points as requested by wp_indx_req, 
     ! map req_indx will reconstruct this order
     integer, dimension(:), pointer :: req_indx => null() !< if present(wp_indx_req), maps wp_indx_req to wp_indx, otherwise identity
  end type wp_vtk_file
!
contains
!------------------------------------------------------------------------
!> \brief initiate naming and geometry structure of vtk file
!! \details Define base filename, format (ascii /binary) and header title of vtk file
!!  and get the wavefield points coordinates from module wavefieldPoints.
!! \param this wavefield point vtk file
!! \param wp wavefield points
!! \param invgrid inversion grid
!! \param filename vtk base name (without '.vtk'). Will optionally be concatenated with a filename extension
!! \param vtk_format 'ASCII' or 'BINARY' indicating the vtk file format
!! \param vtk_title optional second line of vtk file (by default 'data on wavefield points')
!! \param errmsg error message
!! \return error message
!
  subroutine initiateWpVtkFile(this,wp,invgrid,filename,vtk_format,errmsg,vtk_title,wp_indx_req)
    ! incoming
    type (wp_vtk_file) :: this
    type (wavefield_points) :: wp
    type (inversion_grid) :: invgrid
    character(len=*) :: filename,vtk_format
    character(len=*), optional :: vtk_title
    integer, dimension(:), optional :: wp_indx_req
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=17) :: myname = 'initiateWpVtkFile'
    character (len=400) :: errstr
    integer :: num_wp_indx_req
    !
    call addTrace(errmsg,myname)
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
    if(present(vtk_title)) then
       this%title = trim(vtk_title)
    else
       this%title = 'data on wavefield points'
    endif
    if(present(wp_indx_req)) then
       num_wp_indx_req = size(wp_indx_req)
       if(num_wp_indx_req .le. 0) then
          write(errstr,*) "number of incoming requested wavefield points (",num_wp_indx_req,&
               ") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return       
       end if
    end if
    nullify(this%points,this%wp_indx,this%req_indx)
    call getVtkWavefieldPoints(wp,invgrid,this%points,this%wp_indx,errmsg,wp_indx_req,this%req_indx)
    if(.level.errmsg == 2) return
    if(associated(this%points)) then
       this%npoints = size(this%points,2)
    else
       this%npoints = 0
    endif
    if(present(wp_indx_req)) then
       if(this%npoints/=num_wp_indx_req) then
          write(errstr,*) "num of vtk points returned by getVtkWavefieldPoints (",this%npoints,&
               ") not equal to num of wavefield points originally requested to be used (",num_wp_indx_req,&
               ") => there were duplicate or invalid indices in requested index array wp_indx_req, which "//&
               "may indicate inconsistencies of files"
          call add(errmsg,2,trim(errstr),myname)
          return
       endif
    endif
    this%ntot_wp = .ntot.wp
  end subroutine initiateWpVtkFile
!------------------------------------------------------------------------
!> \brief open vtk file to write
!! \param this wavefield points vtk file
!! \param lu file unit
!! \param file_index optional index of file (will be appended to filename base)
!! \param overwrite_in optional logical to indicate behaviour in case file exists, by default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine openWpVtkFile(this,lu,filename_extension,errmsg,overwrite_in)
    ! incoming
    type (wp_vtk_file) :: this
    integer :: lu
    character (len=*) :: filename_extension
    logical, optional :: overwrite_in
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=13) :: myname = 'openWpVtkFile'
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
  end subroutine openWpVtkFile
!------------------------------------------------------------------------
!> \brief write vtk header, points and vertices structure to open vtk file
!! \param this wavefield points vtk file
!! \param lu file unit of file (MUST BE OPENED ALREADY!)
!! \param errmsg error message
!! \return error message
!
  subroutine writeHeaderGeometryWpVtkFile(this,lu,errmsg)
    type (wp_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    ! local
    integer :: ios,i
    character(len=28) :: myname = 'writeHeaderGeometryWpVtkFile'
    character (len=500) :: string
    character (len=1), parameter :: eol_char = char(10)
    logical :: err
    !
    call addTrace(errmsg,myname)
    !
    ! remember with err if there was an error somewhere
    err = .false.
    if(this%is_ascii) then
       ! vkt Header
       write(unit=lu,fmt='(a)',iostat=ios) '# vtk DataFile Version 3.1'  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) trim(this%title)              ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'ASCII'                       ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'DATASET POLYDATA'            ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing vtk Header',myname)
          return
       endif
       ! POINTS
       write(unit=lu,fmt='(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(3e14.6e2)',iostat=ios) this%points                     ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
          return
       endif
       ! VERTICES
       write(unit=lu,fmt='(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints     ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(2i12)',iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)            ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    else ! this%is_ascii
       ! vtk Header
       write(unit=lu,iostat=ios) '# vtk DataFile Version 3.1'//eol_char  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) trim(this%title)//eol_char              ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'BINARY'//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'DATASET POLYDATA'//eol_char            ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing vtk Header',myname)
          return
       endif
       ! POINTS
       write(string,'(a,i12,a)',iostat=ios) 'POINTS ',this%npoints,' float'
       write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%points                             ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
          return
       endif
       ! VERTICES
       write(string,'(a,2i12)',iostat=ios)'VERTICES ',this%npoints,2*this%npoints
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) (/ ((/1,i-1/),i=1,this%npoints) /)          ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    endif ! this%is_ascii
  end subroutine writeHeaderGeometryWpVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and one component float scalar valued point data to vtk file
!! \details The number of incoming real data values must match the number of points this%npoints
!!  of this wp_vtk_file object, because here scalar POINT_DATA (i.e. one scalar value per point)
!!  is added to the vtk file. 
!!  First a file is opened calling openWpVtkFile and header and point and vertex 
!!  geometry information is written to that file calling writeHeaderGeometryWpVtkFile.
!!  Then the incoming data values are added to the vtk file as scalar valued float point data.
!!  If data_indx_is_wp is present and true, incoming values in data are expected for every wavefield point (in 
!!  wavefield point order), otherwise incoming values in data are expected for every wp index as in wp_indx_req
!!  when initiating this. If wp_indx_req was not present when initiating this, the requested indices are
!!  simply all wp indices (in wp order)
!! \param this wavefield points vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param data_indx_is_wp optional logical indicating if there are incoming data for ALL wavefield points or not
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeRealDataWpVtkFile(this,lu,data,errmsg,data_name,file_index,data_indx_is_wp,overwrite)
    ! incoming
    type (wp_vtk_file) :: this
    integer :: lu
    real, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: data_indx_is_wp,overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=22) :: myname = 'writeRealDataWpVtkFile'
    character (len=500) :: string
    character (len=100) :: filename_extension
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata_expect
    logical :: err,map_indx_all
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this wp_vtk_file object is not initiated yet',myname)
       return
    endif
    ! check if number of data is correct (writing point data here)
    ndata_expect = this%npoints; map_indx_all = .false.
    if(present(data_indx_is_wp)) then
       if(data_indx_is_wp) then
          ndata_expect = this%ntot_wp
          map_indx_all = .true.
       end if
    end if
    if(size(data) /= ndata_expect) then
       write(errstr,*) 'number of incoming data ( =',size(data),&
            ') does not match number of expected values ( =',ndata_expect,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    !
    if(present(file_index)) then
       write(filename_extension,"('_',i6.6)") file_index
    else
       filename_extension = ''
    endif
    ! open vtk file to write
    call openWpVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and geometry information to file
    call writeHeaderGeometryWpVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) 'POINT_DATA ',this%npoints        ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,fmt='(e14.6e2)', iostat=ios) data(this%wp_indx)         ; err = err.or.(ios/=0)
       else
          write(unit=lu,fmt='(e14.6e2)', iostat=ios) data(this%req_indx)        ; err = err.or.(ios/=0)
       end if
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) 'POINT_DATA ',this%npoints
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,iostat=ios) data(this%wp_indx)                       ; err = err.or.(ios/=0)
       else
          write(unit=lu,iostat=ios) data(this%req_indx)                      ; err = err.or.(ios/=0)
       end if
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    endif
    !
    ! close file
    close(lu)
  end subroutine writeRealDataWpVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and two component float scalar valued cell data to vtk file
!! \details The number of incoming complex data values must match the number of points this%npoints
!!  of this wp_vtk_file object, because here two component scalar POINT_DATA (i.e. two scalar float values per point)
!!  are added to the vtk file. 
!!  First a file is opened calling openWpVtkFile and header and point and vertex 
!!  geometry information is written to that file calling writeHeaderGeometryWpVtkFile.
!!  Then the incoming data values are added to the vtk file as two component scalar valued float cell data.
!!  If data_indx_is_wp is present and true, incoming values in data are expected for every wavefield point (in 
!!  wavefield point order), otherwise incoming values in data are expected for every wp index as in wp_indx_req
!!  when initiating this. If wp_indx_req was not present when initiating this, the requested indices are
!!  simply all wp indices (in wp order)
!! \param this wavefield points vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param data_indx_is_wp optional logical indicating if there are incoming data for ALL wavefield points or not
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeComplexDataWpVtkFile(this,lu,data,errmsg,data_name,file_index,data_indx_is_wp,overwrite)
    ! incoming
    type (wp_vtk_file) :: this
    integer :: lu
    complex, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: data_indx_is_wp,overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=25) :: myname = 'writeComplexDataWpVtkFile'
    character (len=500) :: string
    character (len=100) :: filename_extension
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata_expect
    logical :: err,map_indx_all
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this wp_vtk_file object is not initiated yet',myname)
       return
    endif
    ! check if number of data is correct (writing point data here)
    ndata_expect = this%npoints; map_indx_all = .false.
    if(present(data_indx_is_wp)) then
       if(data_indx_is_wp) then
          ndata_expect = this%ntot_wp
          map_indx_all = .true.
       end if
    end if
    if(size(data) /= ndata_expect) then
       write(errstr,*) 'number of incoming data ( =',size(data),&
            ') does not match number of expected values ( =',ndata_expect,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    !
    if(present(file_index)) then
       write(filename_extension,"('_',i6.6)") file_index
    else
       filename_extension = ''
    endif
    ! open vtk file to write
    call openWpVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and geometry information to file
    call writeHeaderGeometryWpVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) 'POINT_DATA ',this%npoints        ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 2'
       else
          string = 'SCALARS data float 2'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,fmt='(2e14.6e2)', iostat=ios) data(this%wp_indx)        ; err = err.or.(ios/=0)
       else
          write(unit=lu,fmt='(2e14.6e2)', iostat=ios) data(this%req_indx)       ; err = err.or.(ios/=0)
       end if
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) 'POINT_DATA ',this%npoints
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 2'
       else
          string = 'SCALARS data float 2'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       if(map_indx_all) then
          write(unit=lu,iostat=ios) data(this%wp_indx)                       ; err = err.or.(ios/=0)
       else
          write(unit=lu,iostat=ios) data(this%req_indx)                      ; err = err.or.(ios/=0)
       end if
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    endif
    !
    ! close file
    close(lu)
  end subroutine writeComplexDataWpVtkFile
!------------------------------------------------------------------------
!> \brief open file, write only header and geometry without data to vtk file
!! \param this wavefield points vtk file
!! \param lu file unit
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeWpVtkFile(this,lu,errmsg,overwrite)
    ! incoming
    type (wp_vtk_file) :: this
    integer :: lu
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=14) :: myname = 'writeWpVtkFile'
    character (len=1), parameter :: eol_char = char(10)
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this wp_vtk_file object is not initiated yet',myname)
       return
    endif
    !
    ! open vtk file to write
    call openWpVtkFile(this,lu,'_wp',errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and geometry information to file
    call writeHeaderGeometryWpVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! close file
    close(lu)
  end subroutine writeWpVtkFile
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this wp_vtk_file object
!
  subroutine deallocateWpVtkFile(this)
    type (wp_vtk_file) :: this
    this%filename = ''
    if(associated(this%points)) deallocate(this%points)
    if(associated(this%wp_indx)) deallocate(this%wp_indx)
    if(associated(this%req_indx)) deallocate(this%req_indx)
  end subroutine deallocateWpVtkFile
!
end module wpVtkFile
