!----------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!-------------------------------------------------------------
!  Module describing a list of seismic events
!-------------------------------------------------------------
 module seismicEventList
	use seismicEvent
	use errorMessage
	implicit none
	interface createSeismicEventList
		module procedure createCopySeismicEventList
		module procedure createLinkSeismicEventList
	end interface createSeismicEventList
	interface dealloc; module procedure deallocSeismicEventList; end interface
	interface operator (.nev.); module procedure getNevSeismicEventList; end interface
	interface operator (.csys.); module procedure getCoordinateSystemSeismicEventList; end interface
	interface operator (.event.); module procedure getSelectedEventSeismicEventList; end interface
	interface operator (.events.); module procedure getAllEventsSeismicEventList; end interface
	type seismic_event_list
		private
		character (len=1) :: csys
		integer :: nev
		character (len=132) :: name
		type (seismic_event), dimension(:), pointer :: event => null()
	end type seismic_event_list
!
 contains
!---------------------------------------------------------------
!> \brief Create a seismic network object by copying incoming events
!
	subroutine createCopySeismicEventList(this,name,csys,nev,event)
	type (seismic_event_list) :: this
	character (len=*) :: name,csys
	integer :: nev
	type (seismic_event), dimension(:) :: event
!
	this%name = name
	this%csys = csys
	this%nev = nev
	allocate(this%event(nev))
	this%event(1:nev) = event(1:nev)
	end subroutine createCopySeismicEventList
!---------------------------------------------------------------
!> \brief Create a seismic network object by pointing to incoming events
!
	subroutine createLinkSeismicEventList(this,name,csys,event)
	type (seismic_event_list) :: this
	character (len=*) :: name,csys
	type (seismic_event), dimension(:), pointer :: event
!
	this%name = name
	this%csys = csys
	if(associated(event)) then
		this%nev = size(event)
		this%event => event
	else
		this%nev = 0
		nullify(this%event)
	endif
	end subroutine createLinkSeismicEventList
!----------------------------------------------------------------
!> \brief Deallocate a seismic network object
!
	subroutine deallocSeismicEventList(this)
	type (seismic_event_list) :: this
	integer :: j
	if (associated(this%event)) then
		do j = 1,this%nev
			call dealloc(this%event(j))
		enddo
		deallocate(this%event)
	endif
	end subroutine deallocSeismicEventList
!----------------------------------------------------------------
!> \brief Get coordinate system of events
!
	function getCoordinateSystemSeismicEventList(this) result(res)
	type (seismic_event_list), intent(in) :: this
	character (len=1) :: res
	res = this%csys
	end function getCoordinateSystemSeismicEventList
!----------------------------------------------------------------
!> \brief Get number of events
!
	function getNevSeismicEventList(this) result(res)
	type (seismic_event_list), intent(in) :: this
	integer :: res
	res = this%nev
	end function getNevSeismicEventList
!-----------------------------------------------------------------
!> \brief Get selected seismic event
!
	function getSelectedEventSeismicEventList(this,j) result(res)
	type (seismic_event_list), intent(in) :: this
	integer, intent(in) :: j
	type (seismic_event) :: res
	res = this%event(j)
	end function getSelectedEventSeismicEventList
!-----------------------------------------------------------------
!> \brief Get all seismic events as pointer to array
!
	function getAllEventsSeismicEventList(this) result(res)
	type (seismic_event_list), intent(in) :: this
	type (seismic_event), dimension(:), pointer :: res
	res => this%event
	end function getAllEventsSeismicEventList
!-------------------------------------------------------------------------
!> \brief Return next event
!> \param f Start with event index f
!> \param s go through events with step s
!> \param f End with event index f
!
	function nextEventSeismicEventList(this,event,f,s,l,reset) result(next)
	type (seismic_event_list), intent(in) :: this
	type (seismic_event) :: event
	integer, optional :: f,s,l
	logical, optional :: reset
	logical :: next
	integer :: first,step,last,current,call_count = 0
	save :: call_count
!
	if(present(reset)) then
		if(reset) then
			call_count = 0
			next = .false.
			return
		end if
	endif
!
	if (present(f)) then; first = f; else; first = 1; endif
	if (present(s)) then; step = s; else; step = 1; endif
	if (present(l)) then; last = l; else; last = this%nev; endif
	call_count = call_count+1
	current = first+(call_count-1)*step
	if (current > last) then
		next = .false.
		call_count = 0
	else if (current < 1) then
		next = .false.
		call_count = 0
	else
		event = this%event(current)
		next = .true.
	endif
	end function nextEventSeismicEventList
!--------------------------------------------------------------------------------------
!> \brief Print out seismic network
!
	subroutine printSeismicEventList(this)
	type (seismic_event_list) :: this
	integer :: j
	write(6,'(a15,a6,4a12,a25,a15)') 'Event Id','CS','Latitude','Longitude','Depth','Magnitude','Origin time','Force/Moment'
	do j = 1,this%nev
		call printSeismicEvent(this%event(j))
	enddo
	end subroutine printSeismicEventList
!--------------------------------------------------------------------------------
!> \brief Search event list for eventid, return event as output argument
!
	function searchEventidSeismicEventList(this,eventid,event,iev) result(errmsg)
	type (seismic_event_list) :: this
	character (len=*) :: eventid
	type (seismic_event), optional :: event
	integer, optional :: iev
	type (error_message) :: errmsg
	character (len=29) :: myname = 'searchEventidSeismicEventList'
	integer :: idx,j
!
	idx = 0
	call new(errmsg,myname)
	do j = 1, this%nev
		!if (index(trim(.evid.(this%event(j))),trim(eventid)) > 0) then ! old version: this is not robust!
		if (adjustl(trim(.evid.(this%event(j)))) == adjustl(trim(eventid))) then
			idx = j; exit
		endif
	enddo
	if (idx == 0) then
		call new(errmsg,2,'Event '//eventid//' not found',myname)
		return
	endif
	if(present(event)) event = this%event(idx)
	if(present(iev)) iev = idx
	end function searchEventidSeismicEventList
!
 end module
