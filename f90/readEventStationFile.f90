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
!> \brief read content of ASKI event file and ASKI station file
!!
!! \details define how to read content of ASKI event file and ASKI station file
!!  and return arrays of seismic event, seismic station objects accordingly.
!!  different formats of station and event files could be implemented here
!!
!! \author Florian Schumacher
!! \date March 2013
!
module readEventStationFile
!
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use dateTime
  use errorMessage
!
  implicit none
!
  interface createEventsFromEventFile
     module procedure createEventsFromStandardEventFile
     !module procedure createEventsFromSpecificTypeEventFile ! see template below.
  end interface createEventsFromEventFile
  interface createEventListFromEventFile
     module procedure createEventListFromStandardEventFile
     !module procedure createEventListFromSpecificTypeEventFile ! see template below.
  end interface createEventListFromEventFile
  interface createStationsFromStationFile
     module procedure createStationsFromStandardStationFile
     !module procedure createStationsFromSpecificTypeStationFile ! see template below.
  end interface createStationsFromStationFile
  interface createStationListFromStationFile
     module procedure createStationListFromStandardStationFile
     !module procedure createNetworkFromSpecificTypeStationFile ! see template below.
  end interface createStationListFromStationFile
!
contains
!
!------------------------------------------------------------------------
!> \brief read in standard event list text file and return seismic events objects
!! \param event_file file name of file containing standard event list
!! \param lu file unit to open event file
!! \param events list of seismic event objects created from events present in event file
!! \param errmsg error message
!! \return error message
!
  subroutine createEventsFromStandardEventFile(event_file,lu,events,errmsg)
    ! incoming
    character(len=*) :: event_file
    integer :: lu
    ! returning
    type (seismic_event), dimension(:), pointer :: events
    type (error_message) :: errmsg
    ! local
    character(len=21) :: myname = 'readStandardEventFile'
    character(len=400) :: errstr
    character(len=1) :: csys
    integer :: ios,iev,nev
    character(len=256) :: line
    character(len=13) :: evid
    character(len=25) :: otime_string
    type(date_time) :: otime
    real :: slat,slon,sdepth,mag
    integer :: istyp
    real, dimension(3) :: f
    real, dimension(6) :: M
!
    call addTrace(errmsg,myname)
!
    ! open event file
    open(unit=lu, file=trim(event_file) , status = 'old', action  = 'read', iostat = ios)
    if (ios /= 0) then
       call add(errmsg,2,"could not open event file '"//trim(event_file)//"'",myname)
       return
    else
       call add(errmsg,0,"successfully opened event file '"//trim(event_file)//"'",myname)
    endif
!
    ! check first entry of first line of event file, must be 'C' or 'S', defining the coordinate system of events
    read(lu,"(a256)",iostat = ios) line
    if (ios /= 0) then
       call add(errmsg,2,"could not read the single character defining the coordinate system "//&
            "from first line of event file",myname)
       return
    endif
    read(line,*,iostat=ios) csys
    if(ios /= 0) then
       call add(errmsg,2,"could not read single character defining the coordinate system from first line '"&
            //trim(line)//"' of event file",myname)
       return
    end if
    select case(csys)
    case ('C','S') ! ok, do nothing
    case default
       call add(errmsg,2,"character '"//csys//"' on first line of event file defining the coordinate "//&
            "system is not valid, must bei either 'C' or 'S'",myname)
       return
    end select
!
    ! count number of non-empty lines below first line, interpreted as event definitions (one event per line)
    nev = 0
    do while(ios == 0)
       read(lu,"(a256)",iostat = ios) line
       if(ios /= 0) exit

       if( len_trim(line) > 0 ) nev = nev + 1
    enddo
    close(lu)
    write(errstr,"(a,i4)") "number of non-empty lines found in event file (below first line): ",nev
    call add(errmsg,0,trim(errstr),myname)
!
    allocate(events(nev))
!    
    ! reads in lines of the form:  eventid  origintime  slat  slon  sdepth  istyp  mag  mom/frce
    open(unit=lu, file=trim(event_file), status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) then
       call add(errmsg,2,"could not re-open event file",myname)
       return
    end if
    read(lu,"(a256)",iostat = ios) line
    if (ios /= 0) then
       call add(errmsg,2,"could not read first line of event file",myname)
       return
    endif
    iev = 0
    do while(ios == 0)
       read(lu,"(a256)",iostat = ios) line
       if( ios /= 0 ) exit
       if( len_trim(line) > 0 ) then
          iev = iev + 1
          line = trim(line)
          ! create seismic_event object from line
          read(line,*) evid,otime_string,slat,slon,sdepth,mag,istyp
          call createFromFullTimestringDateTime(otime,otime_string)
          if(istyp==0) then
             read(line,*) evid,otime_string,slat,slon,sdepth,mag,istyp,f
             call createForceMomentSeismicEvent(events(iev),evid,csys,slat,slon,sdepth,mag,otime,istyp,f)
          elseif(istyp==1) then
             read(line,*) evid,otime_string,slat,slon,sdepth,mag,istyp,M
             call createForceMomentSeismicEvent(events(iev),evid,csys,slat,slon,sdepth,mag,otime,istyp,M)
          else
             call createSeismicEvent(events(iev),evid,csys,slat,slon,sdepth,mag,otime)
          end if ! istyp==0
       endif ! len_trim(line) > 0
    enddo ! while(ios == 0)
    close(lu)
  end subroutine createEventsFromStandardEventFile
!------------------------------------------------------------------------
!> \brief read in standard event list text file and return seismic event list object
!! \param event_file file name of file containing standard event list
!! \param lu file unit to open event file
!! \param events list of seismic event objects created from events present in event file
!! \param errmsg error message
!! \return error message
!
  subroutine createEventListFromStandardEventFile(event_file,lu,name,event_list,errmsg)
    ! incoming
    character(len=*) :: event_file,name
    integer :: lu
    ! returning
    type (seismic_event_list) :: event_list
    type (error_message) :: errmsg
    ! local
    type (seismic_event), dimension(:), pointer :: events
!
    call addTrace(errmsg,'createEventListFromStandardEventFile')
    call createEventsFromStandardEventFile(event_file,lu,events,errmsg)
    if(.level.errmsg == 2) return
    if(.not.associated(events)) then
       call add(errmsg,2,'there are no events returned by createEventsFromStandardEventFile',&
            'createEventListFromStandardEventFile')
       return
    end if
    call createSeismicEventList(event_list,name,.csys.events(1),events)
    nullify(events)
  end subroutine createEventListFromStandardEventFile
!------------------------------------------------------------------------
!> \brief read in standard station list text file and return seismic stations objects
!! \param station_file file name of file containing standard station list
!! \param lu file unit to open station file
!! \param stations list of seismic station objects created from stations present in station file
!! \param errmsg error message
!! \return error message
!
  subroutine createStationsFromStandardStationFile(station_file,lu,stations,errmsg)
    ! incoming
    character(len=*) :: station_file
    integer :: lu
    ! returning
    type (seismic_station), dimension(:), pointer :: stations
    type (error_message) :: errmsg
    ! local
    character(len=23) :: myname = 'readStandardStationFile'
    character(len=400) :: errstr
    character(len=1) :: csys
    integer :: ios,istat,nstat
    character(len=256) :: line
    character (len=5) :: staname
    character (len=6) :: netw
    real :: lat, lon, alt
!
    call addTrace(errmsg,myname)
!
    ! open station file
    open(unit=lu, file=trim(station_file) , status = 'old', action  = 'read', iostat = ios)
    if (ios /= 0) then
       call add(errmsg,2,'could not open stations file "'//trim(station_file)//'"',myname)
       return
    else
       call add(errmsg,0,"successfully opened station file '"//trim(station_file)//"'",myname)
    endif
!
    ! check first entry of first line of station file, must be 'C' or 'S', defining the coordinate system of events
    read(lu,"(a256)",iostat = ios) line
    if (ios /= 0) then
       call add(errmsg,2,"could not read the single character defining the coordinate system "//&
            "from first line of station file",myname)
       return
    endif
    read(line,*,iostat=ios) csys
    if(ios /= 0) then
       call add(errmsg,2,"could not read single character defining the coordinate system from first line '"&
            //trim(line)//"' of station file",myname)
       return
    end if
    select case(csys)
    case ('C','S') ! ok, do nothing
    case default
       call add(errmsg,2,"character '"//csys//"' on first line of event file defining the coordinate "//&
            "system is not valid, must bei either 'C' or 'S'",myname)
       return
    end select
!
    ! count number of non-empty lines below first line, interpreted as station definitions (one station per line)
    nstat = 0
    do while(ios == 0)
       read(lu,"(a256)",iostat = ios) line
       if(ios /= 0) exit

       if( len_trim(line) > 0 ) nstat = nstat + 1
    enddo
    close(lu)
    write(errstr,"(a,i4)") "number of non-empty lines found in stations file (below first line): ",nstat
    call add(errmsg,0,trim(errstr),myname)
!
    allocate(stations(nstat))
!
    ! reads in lines of the form:  station_name  network_name  latitude  longitude  elevation
    open(unit=lu, file=trim(station_file), status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) then
       call add(errmsg,2,"could not re-open station file",myname)
       return
    end if
    read(lu,"(a256)",iostat = ios) line
    if (ios /= 0) then
       call add(errmsg,2,"could not read first line of station file",myname)
       return
    endif
    istat = 0
    do while(ios == 0)
       read(lu,"(a256)",iostat = ios) line
       if( ios /= 0 ) exit
       if( len_trim(line) > 0 ) then
          istat = istat + 1
          line = trim(line)
          ! create seismic_station object from line
          read(line,*) staname, netw, lat, lon, alt
          call createSeismicStation(stations(istat),csys,trim(staname),lat,lon,alt,netcode=trim(netw))
       endif
    enddo
    close(lu)
  end subroutine createStationsFromStandardStationFile
!------------------------------------------------------------------------
!> \brief read in standard station list text file and return seismic stations objects
!! \param station_file file name of file containing standard station list
!! \param lu file unit to open station file
!! \param stations list of seismic station objects created from stations present in station file
!! \param errmsg error message
!! \return error message
!
  subroutine createStationListFromStandardStationFile(station_file,lu,name,network,errmsg)
    ! incoming
    character(len=*) :: station_file,name
    integer :: lu
    ! returning
    type (seismic_network) :: network
    type (error_message) :: errmsg
    ! local
    type (seismic_station), dimension(:), pointer :: stations
!
    call addTrace(errmsg,'createStationListFromStandardStationFile')
    call createStationsFromStandardStationFile(station_file,lu,stations,errmsg)
    if(.level.errmsg == 2) return
    if(.not.associated(stations)) then
       call add(errmsg,2,'there are no stations returned by createStationsFromStandardStationFile',&
            'createStationListFromStandardStationFile')
       return
    end if
    call createSeismicNetwork(network,name,.csys.stations(1),stations)
    nullify(stations)
  end subroutine createStationListFromStandardStationFile
!------------------------------------------------------------------------
!!$!> \brief check specific type of event list and call respective read function
!!$!! \param event_file file name of file containing standard event list
!!$!! \param lu file unit to open event file
!!$!! \param file_type specific type of event file
!!$!! \param events list of seismic event objects created from events present in event file
!!$!! \param errmsg error message
!!$!! \return error message
!!$!
!!$  subroutine createEventsFromSpecificTypeEventFile(event_file,lu,file_type,events,errmsg)
!!$    character(len=*) :: event_file
!!$    integer :: lu
!!$    character(len=*) :: file_type
!!$    type (seismic_event), dimension(:), pointer :: events
!!$    type (error_message) :: errmsg
!!$!
!!$    call addTrace(errmsg,'readSpecificTypeEventFile')
!!$    select case(file_type)
!!$    case('type1')
!!$       errmsg = createEventsFromType1EventFile(event_file,lu,events)
!!$    case('type2')
!!$       ....
!!$
!!$    case default
!!$       call add(errmsg,2,"event file type '"//trim(file_type)//"' not supported",'readSpecificTypeEventFile')
!!$       nullify(events)
!!$       return
!!$    end select
!!$  end subroutine createEventsFromSpecificTypeEventFile
!
end module readEventStationFile
