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
!> \brief module to write data on event,station coordinates or pahts to vkt output
!!
!! \details A vtk POLYDATA dataset consisting of event or station coordinate points 
!!  as VERTICES along with any scalar POINT DATA living on the event or station coordinate 
!!  points is written to binary or ascii vtk files. Alternatively, any data living on paths, 
!!  defined as vtk LINES, can be written to vtk files.
!!  Complex data is handled as 2 component scalar float data. As an option, multiple 
!!  files containing data w.r.t. some index (frequency, time) may be written having the same
!!  file base name followed by an index, in order to be considered by Paraview as a sequence of data.
!!
!! \author Florian Schumacher
!! \date March 2013
!
module eventStationVtkFile
!
  use seismicEventList
  use seismicEvent
  use seismicNetwork
  use seismicStation
  use inversionGrid
  use fileUnitHandler
  use errorMessage
!
  implicit none
!
  ! by allowing calls to only certain routines, it is tried to assure that e.g. 
  ! a file is open when writing to it, etc. (this way, certain security checks
  ! may be omitted, as unintended calls to auxilliary routines are forbidden)
  private
  public :: event_station_vtk_file,init,writeData,writeEvents,writeStations,writePaths,dealloc
!
  interface init
     module procedure initiateEventsVtkFile
     module procedure initiateStationsVtkFile
     module procedure initiatePathsVtkFile
  end interface init
  interface writeData
     module procedure writeRealDataEventStationVtkFile
     !module procedure writeRealVectorDataEventStationVtkFile ! for multidimensional data (suitable on paths e.g. for focussing coefficients per path, for all components), can be (1?),2,3 dimensional
     module procedure writeComplexDataEventStationVtkFile
  end interface writeData
  interface writeEvents; module procedure writeEventsVtkFile; end interface
  interface writeStations; module procedure writeStationsVtkFile; end interface
  interface writePaths; module procedure writePathsVtkFile; end interface
  interface dealloc; module procedure deallocateEventStationVtkFile; end interface
!
!> \brief general file, geometry and cell information of vtk file
  type event_station_vtk_file
     private
     integer :: type_of_file = 0 !< 1 = events vtk file, 2 = stations vtk file, 3 = paths vtk file
     character (len=300) :: filename = '' !< (base) file name of vtk file (WITHOUT '.vtk')
     logical :: is_ascii !< indicating ascii (true) or binary (false) format of vtk file
     character (len=100) :: title = '' !< second line of vtk file
     integer :: nstat = 0 !< number of stations
     integer :: nev = 0 !< number of events
     integer :: npath = 0 !< number of paths
     integer :: npoint = 0 !< number of points (size of second dimension of array points)
     real, dimension(:,:), pointer :: points => null() !< POINTS geometry for POLYDATA VERTICES, coordinates of stations or events according to type_of_file
     integer, dimension(:,:), pointer :: lines => null() !< LINES geometry for POLYDATA LINES, array of point indices defining the paths
  end type event_station_vtk_file
!
contains
!------------------------------------------------------------------------
!> \brief initiate naming and geometry structure of events vtk file
!! \details Define base filename, format (ascii /binary) and header title of vtk file
!!  and get set the respective point coordinates
!! \param this event station vtk file
!! \param events seismic event list
!! \param invgrid inversion grid
!! \param filename vtk base name (without '.vtk'). Will optionally be concatenated with a filename extension
!! \param vtk_format 'ASCII' or 'BINARY' indicating the vtk file format
!! \param vtk_title optional second line of vtk file (by default 'data on event coordinates')
!! \param errmsg error message
!! \return error message
!
  subroutine initiateEventsVtkFile(this,events,invgrid,filename,vtk_format,errmsg,vtk_title)
    ! incoming
    type (event_station_vtk_file) :: this
    type (seismic_event_list) :: events
    type (inversion_grid) :: invgrid
    character(len=*) :: filename,vtk_format
    character(len=*), optional :: vtk_title
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=21) :: myname = 'initiateEventsVtkFile'
    character (len=400) :: errstr
    type (seismic_event) :: event
    integer :: ipoint
    real, dimension(:), allocatable :: c1,c2,c3
!
    call addTrace(errmsg,myname)
    if(trim(this%filename) /= '') call deallocateEventStationVtkFile(this)
!
    this%type_of_file = 1
    this%nev = .nev.events
    if(this%nev .le. 0) then
       write(errstr,*) "there are ",this%nev," events in incoming event list. there must be a positive number of events"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
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
       this%title = 'data on event coordinates'
    endif
!
    allocate(c1(this%nev),c2(this%nev),c3(this%nev))
    ipoint = 0
    do while (nextEventSeismicEventList(events,event))
       ipoint = ipoint + 1
       c1(ipoint) = .slat.event
       c2(ipoint) = .slon.event
       c3(ipoint) = .sdepth.event
    end do ! while (nextEvent)
!
    call transformToVtkInversionGrid(invgrid,c1,c2,c3,'event',errmsg)
    if(.level.errmsg==2) goto 1
!
    this%npoint = this%nev
    allocate(this%points(3,this%npoint))
    this%points(1,:) = c1
    this%points(2,:) = c2
    this%points(3,:) = c3
!
1   if(allocated(c1)) deallocate(c1)
    if(allocated(c2)) deallocate(c2)
    if(allocated(c3)) deallocate(c3)
  end subroutine initiateEventsVtkFile
!------------------------------------------------------------------------
!> \brief initiate naming and geometry structure of stations vtk file
!! \details Define base filename, format (ascii /binary) and header title of vtk file
!!  and get set the respective point coordinates
!! \param this event station vtk file
!! \param stations seismic network
!! \param invgrid inversion grid
!! \param filename vtk base name (without '.vtk'). Will optionally be concatenated with a filename extension
!! \param vtk_format 'ASCII' or 'BINARY' indicating the vtk file format
!! \param vtk_title optional second line of vtk file (by default 'data on event coordinates')
!! \param errmsg error message
!! \return error message
!
  subroutine initiateStationsVtkFile(this,stations,invgrid,filename,vtk_format,errmsg,vtk_title)
    ! incoming
    type (event_station_vtk_file) :: this
    type (seismic_network) :: stations
    type (inversion_grid) :: invgrid
    character(len=*) :: filename,vtk_format
    character(len=*), optional :: vtk_title
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=21) :: myname = 'initiateEventsVtkFile'
    character (len=400) :: errstr
    type (seismic_station) :: station
    integer :: ipoint
    real, dimension(:), allocatable :: c1,c2,c3
!
    call addTrace(errmsg,myname)
    if(trim(this%filename) /= '') call deallocateEventStationVtkFile(this)
!
    this%type_of_file = 2
    this%nstat = .nstat.stations
    if(this%nstat .le. 0) then
       write(errstr,*) "there are ",this%nstat," stations in incoming network. there must be a positive number of stations"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
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
       this%title = 'data on station coordinates'
    endif
!
    allocate(c1(this%nstat),c2(this%nstat),c3(this%nstat))
    ipoint = 0
    do while (nextStationSeismicNetwork(stations,station))
       ipoint = ipoint + 1
       c1(ipoint) = .lat.station
       c2(ipoint) = .lon.station
       c3(ipoint) = .alt.station
    end do ! while (nextStation)
!
    call transformToVtkInversionGrid(invgrid,c1,c2,c3,'station',errmsg)
    if(.level.errmsg==2) goto 1
!
    this%npoint = this%nstat
    allocate(this%points(3,this%npoint))
    this%points(1,:) = c1
    this%points(2,:) = c2
    this%points(3,:) = c3
!
1   if(allocated(c1)) deallocate(c1)
    if(allocated(c2)) deallocate(c2)
    if(allocated(c3)) deallocate(c3)
  end subroutine initiateStationsVtkFile
!------------------------------------------------------------------------
!> \brief initiate naming and geometry structure of paths vtk file
!! \details Define base filename, format (ascii /binary) and header title of vtk file
!!  and get set the respective point coordinates
!! \param this event station vtk file
!! \param events seismic event list
!! \param stations seismic network
!! \param paths array of dimension(2,npath) which contains event IDs (first entry) and station names (second entry)
!! \param invgrid inversion grid
!! \param filename vtk base name (without '.vtk'). Will optionally be concatenated with a filename extension
!! \param vtk_format 'ASCII' or 'BINARY' indicating the vtk file format
!! \param vtk_title optional second line of vtk file (by default 'data on event coordinates')
!! \param errmsg error message
!! \return error message
!
  subroutine initiatePathsVtkFile(this,events,stations,paths,invgrid,filename,vtk_format,errmsg,vtk_title)
    ! incoming
    type (event_station_vtk_file) :: this
    type (seismic_event_list) :: events
    type (seismic_network) :: stations
    character(len=*), dimension(:,:) :: paths
    type (inversion_grid) :: invgrid
    character(len=*) :: filename,vtk_format
    character(len=*), optional :: vtk_title
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=20) :: myname = 'initiatePathsVtkFile'
    character (len=400) :: errstr
    type (error_message) :: errmsg2
    type (seismic_event) :: event
    type (seismic_station) :: station
    integer :: ipoint,ipath,nev,iev,nstat,istat
    real, dimension(:), allocatable :: ev_c1,ev_c2,ev_c3,stat_c1,stat_c2,stat_c3
    integer, dimension(:), allocatable :: ipoint_of_station,ipoint_of_event
    real, dimension(:,:), pointer :: points_tmp
!
    nullify(points_tmp)
    call addTrace(errmsg,myname)
    if(trim(this%filename) /= '') call deallocateEventStationVtkFile(this)
!
    this%type_of_file = 3
    this%npath = size(paths,2)
    if(this%npath .le. 0 .or. size(paths,1)/=2) then
       write(errstr,*) "there are ",this%npath," paths in incoming  path array (must be positive), and ",size(paths,1),&
            " names per path (must be 2: 'eventID','station name' pair)"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
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
       this%title = 'data on paths'
    endif
!
    nev = .nev.events
    if(nev .le. 0) then
       write(errstr,*) "there are ",nev," events in incoming event list. there must be a positive number of events"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    nstat = .nstat.stations
    if(nstat .le. 0) then
       write(errstr,*) "there are ",nstat," stations in incoming network. there must be a positive number of stations"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    ! FIRST TRANSFORM ALL STATION AND EVENT COORDINATE TO VTK FRAME BY INVERSION GRID MODULE
!
    allocate(stat_c1(nstat),stat_c2(nstat),stat_c3(nstat))
    ipoint = 0
    do while (nextStationSeismicNetwork(stations,station))
       ipoint = ipoint + 1
       stat_c1(ipoint) = .lat.station
       stat_c2(ipoint) = .lon.station
       stat_c3(ipoint) = .alt.station
    end do ! while (nextStation)
!
    call transformToVtkInversionGrid(invgrid,stat_c1,stat_c2,stat_c3,'station',errmsg)
    if(.level.errmsg==2) goto 1
!
    allocate(ev_c1(nev),ev_c2(nev),ev_c3(nev))
    ipoint = 0
    do while (nextEventSeismicEventList(events,event))
       ipoint = ipoint + 1
       ev_c1(ipoint) = .slat.event
       ev_c2(ipoint) = .slon.event
       ev_c3(ipoint) = .sdepth.event
    end do ! while (nextEvent)
!
    call transformToVtkInversionGrid(invgrid,ev_c1,ev_c2,ev_c3,'event',errmsg)
    if(.level.errmsg==2) goto 1
!
    ! SECONDLY, SELECT THOSE STATIONS AND EVENTS THAT ARE CONTAINED IN THE PATHS AT ALL
    ! AND STORE THEIR COORDINATES IN ARRAY THIS%POINTS. REMEMBER THEIR POINT INDICES 
    ! IN ORDER TO CORRECTLY CONNECT THE LINES GEOMETRY BELOW
!
    allocate(points_tmp(3,nstat+nev))
    ipoint = 0
!
    ! first iterate through the events
    allocate(ipoint_of_event(nev)); ipoint_of_event = -1
    iev = 0
    do while (nextEventSeismicEventList(events,event))
       iev = iev + 1
       if(any(paths(1,:)==.evid.event)) then
          ipoint = ipoint + 1
          ! add the event coordinates to the (temporary) point array
          points_tmp(:,ipoint) = (/ ev_c1(iev),ev_c2(iev),ev_c3(iev) /)
          ! connect this point to the event index in the event list
          ipoint_of_event(iev) = ipoint
       end if
    end do

    ! then iterate through the stations, keep the value of ipoint!! (simply continue)
    allocate(ipoint_of_station(nstat)); ipoint_of_station = -1
    istat = 0
    do while (nextStationSeismicNetwork(stations,station))
       istat = istat + 1
       if(any(paths(2,:)==.staname.station)) then
          ipoint = ipoint + 1
          ! add the station coordinates to the (temporary) point array
          points_tmp(:,ipoint) = (/ stat_c1(istat),stat_c2(istat),stat_c3(istat) /)
          ! connect this point to the event index in the event list
          ipoint_of_station(istat) = ipoint
       end if
    end do

    ! finally define the actual points array this%points
    if(ipoint == nstat+nev) then
       this%npoint = nstat+nev
       this%points => points_tmp
       nullify(points_tmp)
    else
       this%npoint = ipoint
       allocate(this%points(3,this%npoint))
       this%points = points_tmp(:,1:ipoint)
       deallocate(points_tmp)
    end if
!
    ! THIRDLY, DEFINE THE LINES GEOMETRY: FOR EACH GIVEN PATH CONNECT THE TWO POINTS CORRESPONDING TO 
    ! THE RESPECTIVE STATION AND EVENT
!
    allocate(this%lines(3,this%npath))
    this%lines(1,:) = 2 ! always have 2 point indices following defining a line
    do ipath = 1,this%npath
       ! first search for the event index of this event in the incoming event list
       errmsg2 = searchEventidSeismicEventList(events,paths(1,ipath),iev=iev)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "event ID '"//trim(paths(1,ipath))//"' of ",ipath,"'th path was not found in event list "//&
               "(the corresponding station name of this path is '"//trim(paths(2,ipath))//"'"
          call add(errmsg,2,trim(errstr),myname)
          call deallocateEventStationVtkFile(this)
          goto 1
       end if
       call dealloc(errmsg2)
!
       ! define geometry connectivity array LINE for this path, event coordinates are first point
       this%lines(2,ipath) = ipoint_of_event(iev)-1
!
       ! then search for the station index of this station in the incoming station list
       errmsg2 =  searchStationNameSeismicNetwork(stations,paths(2,ipath),istat=istat)
       if(.level.errmsg2 == 2) then
          write(errstr,*) "station name '"//trim(paths(2,ipath))//"' of ",ipath,"'th path was not found in station list"//&
               "(the corresponding event ID of this path is '"//trim(paths(1,ipath))//"'"
          call add(errmsg,2,trim(errstr),myname)
          call deallocateEventStationVtkFile(this)
          goto 1
       end if
       call dealloc(errmsg2)
!
       ! define geometry connectivity array LINE for this path, station coordinates are second point
       this%lines(3,ipath) = ipoint_of_station(istat)-1
    end do ! ipath
!
1   if(allocated(ev_c1)) deallocate(ev_c1)
    if(allocated(ev_c2)) deallocate(ev_c2)
    if(allocated(ev_c3)) deallocate(ev_c3)
    if(allocated(stat_c1)) deallocate(stat_c1)
    if(allocated(stat_c2)) deallocate(stat_c2)
    if(allocated(stat_c3)) deallocate(stat_c3)
    if(allocated(ipoint_of_station)) deallocate(ipoint_of_station)
    if(allocated(ipoint_of_event)) deallocate(ipoint_of_event)
  end subroutine initiatePathsVtkFile
!------------------------------------------------------------------------
!> \brief open vtk file to write
!! \param this event station vtk file
!! \param lu file unit
!! \param file_index optional index of file (will be appended to filename base)
!! \param overwrite_in optional logical to indicate behaviour in case file exists, by default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine openEventStationVtkFile(this,lu,filename_extension,errmsg,overwrite_in)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    character (len=*) :: filename_extension
    logical, optional :: overwrite_in
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=23) :: myname = 'openEventStationVtkFile'
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
  end subroutine openEventStationVtkFile
!------------------------------------------------------------------------
!> \brief write vtk header and points structure to open vtk file
!! \param this event station vtk file
!! \param lu file unit of file (MUST BE OPENED ALREADY!)
!! \param errmsg error message
!! \return error message
!
  subroutine writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    type (event_station_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    ! local
    integer :: ios
    character(len=38) :: myname = 'writeHeaderGeometryEventStationVtkFile'
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
       write(unit=lu,fmt='(a,i12,a)',iostat=ios) 'POINTS ',this%npoint,' float'  ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(3e14.6e2)',iostat=ios) this%points                     ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
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
       write(string,'(a,i12,a)',iostat=ios) 'POINTS ',this%npoint,' float'
       write(unit=lu,iostat=ios) trim(string)//eol_char                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%points                             ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINTS',myname)
          return
       endif
    endif ! this%is_ascii
  end subroutine writeHeaderPointsEventStationVtkFile
!------------------------------------------------------------------------
!> \brief write POLYDATA VERTICES structure (one vertex for each defined point) to open vtk file
!! \param this event station vtk file
!! \param lu file unit of file (MUST BE OPENED ALREADY!)
!! \param errmsg error message
!! \return error message
!
  subroutine writeVerticesEventStationVtkFile(this,lu,errmsg)
    type (event_station_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    ! local
    integer :: ios,i
    character(len=32) :: myname = 'writeVerticesEventStationVtkFile'
    character (len=500) :: string
    character (len=1), parameter :: eol_char = char(10)
    logical :: err
    !
    call addTrace(errmsg,myname)
    !
    ! remember with err if there was an error somewhere
    err = .false.
    if(this%is_ascii) then
       ! VERTICES
       write(unit=lu,fmt='(a,2i12)',iostat=ios)'VERTICES ',this%npoint,2*this%npoint     ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(i12)',iostat=ios) (/ ((/1,i-1/),i=1,this%npoint) /)            ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    else ! this%is_ascii
       ! VERTICES
       write(string,'(a,2i12)',iostat=ios)'VERTICES ',this%npoint,2*this%npoint
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) (/ ((/1,i-1/),i=1,this%npoint) /)          ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing VERTICES',myname)
          return
       endif
    endif ! this%is_ascii
  end subroutine writeVerticesEventStationVtkFile
!------------------------------------------------------------------------
!> \brief write POLYDATA LINES structure (basically content of this%path array) to open vtk file
!! \param this event station vtk file
!! \param lu file unit of file (MUST BE OPENED ALREADY!)
!! \param errmsg error message
!! \return error message
!
  subroutine writeLinesEventStationVtkFile(this,lu,errmsg)
    type (event_station_vtk_file) :: this
    integer :: lu
    type(error_message) :: errmsg
    ! local
    integer :: ios
    character(len=29) :: myname = 'writeLinesEventStationVtkFile'
    character (len=500) :: string
    character (len=1), parameter :: eol_char = char(10)
    logical :: err
    !
    call addTrace(errmsg,myname)
    !
    ! remember with err if there was an error somewhere
    err = .false.
    if(this%is_ascii) then
       ! LINES
       write(unit=lu,fmt='(a,2i12)',iostat=ios)'LINES ',this%npath,size(this%lines)   ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(3i12)',iostat=ios) this%lines                              ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing LINES',myname)
          return
       endif
    else ! this%is_ascii
       ! LINES
       write(string,'(a,2i12)',iostat=ios)'LINES ',this%npath,size(this%lines)
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) this%lines                                  ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing LINES',myname)
          return
       endif
    endif ! this%is_ascii
  end subroutine writeLinesEventStationVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and one component float scalar valued point (cell) data to vtk file
!! \details The number of incoming real data values must match the number of points (i.e. this%nev,this%nstat,
!!  respectively) or number of cells (lines, i.e. this%npath), dependent on this%type_of_file,
!!  because here scalar POINT_DATA (CELL_DATA) (i.e. one scalar value per point (cell))
!!  is added to the vtk file. 
!!  First a file is opened calling openEventStationVtkFile and header and point and vertex (lines)
!!  geometry information is written to that file calling writeHeaderPointsEventStationVtkFile, 
!!  writeVerticesEventStationVtkFile,writeLinesEventStationVtkFile.
!!  Then the incoming data values are added to the vtk file as scalar valued float point data.
!! \param this event station vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeRealDataEventStationVtkFile(this,lu,data,errmsg,data_name,file_index,overwrite)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    real, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=32) :: myname = 'writeRealDataEventStationVtkFile'
    character (len=500) :: string
    character (len=100) :: filename_extension
    character (len=10) :: dataset_attribute ! either "POINT_DATA", or "CELL_DATA"
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata
    logical :: err
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this event_station_vtk_file object is not initiated yet',myname)
       return
    endif
    !
    ! check if number of incoming data is correct, depending on type of file
    select case(this%type_of_file)
    case(1)
       ndata = this%nev
       write(errstr,*) 'number of data ( =',size(data),') does not match number of events ( =',this%nev,')'
    case(2)
       ndata = this%nstat
       write(errstr,*) 'number of data ( =',size(data),') does not match number of stations ( =',this%nstat,')'
    case(3)
       ndata = this%npath
       write(errstr,*) 'number of data ( =',size(data),') does not match number of paths ( =',this%npath,')'
    end select
    if(size(data) /= ndata) then
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
    call openEventStationVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and points information to file
    call writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write vertices or lines to file, depending on type of file
    select case(this%type_of_file)
    case(1,2)
       dataset_attribute = 'POINT_DATA'
       call writeVerticesEventStationVtkFile(this,lu,errmsg)
       if(.level.errmsg == 2) then
          close(lu)
          return
       endif
    case(3)
       dataset_attribute = 'CELL_DATA'
       call writeLinesEventStationVtkFile(this,lu,errmsg)
       if(.level.errmsg == 2) then
          close(lu)
          return
       endif
    end select
    !
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) trim(dataset_attribute)//' ',ndata        ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(e14.6e2)', iostat=ios) data                          ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) trim(dataset_attribute)//' ',ndata
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) data                                        ; err = err.or.(ios/=0)
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
  end subroutine writeRealDataEventStationVtkFile
!------------------------------------------------------------------------
!> \brief open file, write header and geometry and two component scalar valued float data to vtk file
!! \details The number of incoming complex data values must match the number of points (i.e. this%nev,this%nstat,
!!  respectively) or number of cells (lines, i.e. this%npath), dependent on this%type_of_file,
!!  because here scalar POINT_DATA (CELL_DATA) (i.e. one scalar value per point (cell))
!!  is added to the vtk file. 
!!  First a file is opened calling openEventStationVtkFile and header and point and vertex (lines)
!!  geometry information is written to that file calling writeHeaderPointsEventStationVtkFile, 
!!  writeVerticesEventStationVtkFile,writeLinesEventStationVtkFile.
!!  Then the incoming data values are added to the vtk file as two component scalar valued float data.
!! \param this event station vtk file
!! \param lu file unit
!! \param data data values to be added to vtk file
!! \param data_name optional name of the data (by default 'data')
!! \param file_index optional index of file (will be appended to filename base)
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeComplexDataEventStationVtkFile(this,lu,data,errmsg,data_name,file_index,overwrite)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    complex, dimension(:) :: data
    character (len=*), optional :: data_name
    integer, optional :: file_index
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character (len=400) :: errstr
    character(len=35) :: myname = 'writeComplexDataEventStationVtkFile'
    character (len=500) :: string
    character (len=100) :: filename_extension
    character (len=10) :: dataset_attribute ! either "POINT_DATA", or "CELL_DATA"
    character (len=1), parameter :: eol_char = char(10)
    integer :: ios,ndata
    logical :: err
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this event_station_vtk_file object is not initiated yet',myname)
       return
    endif
    !
    ! check if number of incoming data is correct, depending on type of file
    select case(this%type_of_file)
    case(1)
       ndata = this%nev
       write(errstr,*) 'number of data ( =',size(data),') does not match number of events ( =',this%nev,')'
    case(2)
       ndata = this%nstat
       write(errstr,*) 'number of data ( =',size(data),') does not match number of stations ( =',this%nstat,')'
    case(3)
       ndata = this%npath
       write(errstr,*) 'number of data ( =',size(data),') does not match number of paths ( =',this%npath,')'
    end select
    if(size(data) /= ndata) then
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
    call openEventStationVtkFile(this,lu,trim(filename_extension),errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and points information to file
    call writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write vertices or lines to file, depending on type of file
    select case(this%type_of_file)
    case(1,2)
       dataset_attribute = 'POINT_DATA'
       call writeVerticesEventStationVtkFile(this,lu,errmsg)
       if(.level.errmsg == 2) then
          close(lu)
          return
       endif
    case(3)
       dataset_attribute = 'CELL_DATA'
       call writeLinesEventStationVtkFile(this,lu,errmsg)
       if(.level.errmsg == 2) then
          close(lu)
          return
       endif
    end select
    !
    ! write data to file
    err = .false.
    if(this%is_ascii) then
       write(unit=lu,fmt='(a,i12)',iostat=ios) trim(dataset_attribute)//' ',ndata        ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 2'
       else
          string = 'SCALARS data float 2'
       endif
       write(unit=lu,fmt='(a)',iostat=ios) trim(string)                         ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(a)',iostat=ios) 'LOOKUP_TABLE default'               ; err = err.or.(ios/=0)
       write(unit=lu,fmt='(2e14.6e2)', iostat=ios) data                         ; err = err.or.(ios/=0)
       if(err) then ! if any of the above ios were /= 0
          call add(errmsg,2,'there was an error writing POINT_DATA',myname)
          close(lu)
          return
       endif
    else
       write(string,fmt='(a,i12)',iostat=ios) trim(dataset_attribute)//' ',ndata
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       if(present(data_name)) then 
          string = 'SCALARS '//trim(data_name)//' float 1'
       else
          string = 'SCALARS data float 1'
       endif
       write(unit=lu,iostat=ios) trim(string)//eol_char                      ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) 'LOOKUP_TABLE default'//eol_char            ; err = err.or.(ios/=0)
       write(unit=lu,iostat=ios) data                                        ; err = err.or.(ios/=0)
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
  end subroutine writeComplexDataEventStationVtkFile
!------------------------------------------------------------------------
!> \brief open file, write only header and geometry without data to vtk file
!! \param this event station vtk file
!! \param lu file unit
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeEventsVtkFile(this,lu,errmsg,overwrite)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=17) :: myname = 'writeEventVtkFile'
    character (len=1), parameter :: eol_char = char(10)
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this event_station_vtk_file object is not initiated yet',myname)
       return
    endif
    if(this%type_of_file /= 1) then
       call add(errmsg,2,'it appears that this event_station_vtk_file object was not initiated with event information',myname)
       return
    end if
    !
    ! open vtk file to write
    call openEventStationVtkFile(this,lu,'_events',errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and points information to file
    call writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write vertices to file
    call writeVerticesEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! close file
    close(lu)
  end subroutine writeEventsVtkFile
!------------------------------------------------------------------------
!> \brief open file, write only header and geometry without data to vtk file
!! \param this event station vtk file
!! \param lu file unit
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writeStationsVtkFile(this,lu,errmsg,overwrite)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=20) :: myname = 'writeStationsVtkFile'
    character (len=1), parameter :: eol_char = char(10)
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this event_station_vtk_file object is not initiated yet',myname)
       return
    endif
    if(this%type_of_file /= 2) then
       call add(errmsg,2,'it appears that this event_station_vtk_file object was not initiated with station information',myname)
       return
    end if
    !
    ! open vtk file to write
    call openEventStationVtkFile(this,lu,'_stations',errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and points information to file
    call writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write vertices to file
    call writeVerticesEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! close file
    close(lu)
  end subroutine writeStationsVtkFile
!------------------------------------------------------------------------
!> \brief open file, write only header and geometry without data to vtk file
!! \param this event station vtk file
!! \param lu file unit
!! \param overwrite_in optional logical to indicate behaviour in case file exists. By default, no overwrite
!! \param errmsg error message
!! \return error message
!
  subroutine writePathsVtkFile(this,lu,errmsg,overwrite)
    ! incoming
    type (event_station_vtk_file) :: this
    integer :: lu
    logical, optional :: overwrite
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=17) :: myname = 'writePathsVtkFile'
    character (len=1), parameter :: eol_char = char(10)
    !
    call addTrace(errmsg,myname)
    !
    ! check if object is initiated already
    if(this%filename== '') then
       call add(errmsg,2,'it appears that this event_station_vtk_file object is not initiated yet',myname)
       return
    endif
    if(this%type_of_file /= 3) then
       call add(errmsg,2,'it appears that this event_station_vtk_file object was not initiated with path information',myname)
       return
    end if
    !
    ! open vtk file to write
    call openEventStationVtkFile(this,lu,'_paths',errmsg,overwrite)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write header and points information to file
    call writeHeaderPointsEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! write LINES to file
    call writeLinesEventStationVtkFile(this,lu,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       return
    endif
    !
    ! close file
    close(lu)
  end subroutine writePathsVtkFile
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this event_station_vtk_file object
!
  subroutine deallocateEventStationVtkFile(this)
    type (event_station_vtk_file) :: this
    this%type_of_file = 0
    this%filename = ''
    this%title = ''
    this%nev = 0
    this%nstat = 0
    this%npath = 0
    this%npoint = 0
    if(associated(this%points)) deallocate(this%points)
    if(associated(this%lines)) deallocate(this%lines)
  end subroutine deallocateEventStationVtkFile
!
end module eventStationVtkFile
