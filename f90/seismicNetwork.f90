!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!-------------------------------------------------------------
!  Module describing a seismic network
!-------------------------------------------------------------
 module seismicNetwork
	use seismicStation
	use errorMessage
	implicit none
	interface createSeismicNetwork
		module procedure createCopySeismicNetwork
		module procedure createLinkSeismicNetwork
	end interface createSeismicNetwork
	interface dealloc; module procedure deallocSeismicNetwork; end interface
	interface operator (.csys.); module procedure getCoordinateSystemSeismicNetwork; end interface
	interface operator (.nstat.); module procedure getNstatSeismicNetwork; end interface
	interface operator (.station.); module procedure getSelectedStationSeismicNetwork; end interface
	interface operator (.stations.); module procedure getAllStationsSeismicNetwork; end interface
	type seismic_network
		private
		integer :: nstat
		character (len=1) :: csys
		character (len=6) :: netcode
		type (seismic_station), dimension(:), pointer :: station => null()
	end type
!
 contains
!---------------------------------------------------------------
!> \brief Create a seismic network object by copying incoming stations
!
	subroutine createCopySeismicNetwork(this,netcode,csys,nstat,station)
	type (seismic_network) :: this
	character (len=*) :: netcode,csys
	integer :: nstat
	type (seismic_station), dimension(:) :: station
!
	this%nstat = nstat
	this%netcode = netcode
	this%csys = csys
	allocate(this%station(nstat))
	this%station(1:nstat) = station(1:nstat)
	end subroutine createCopySeismicNetwork
!---------------------------------------------------------------
!> \brief Create a seismic network object by pointing to incoming stations
!
	subroutine createLinkSeismicNetwork(this,netcode,csys,station)
	type (seismic_network) :: this
	character (len=*) :: netcode,csys
	type (seismic_station), dimension(:), pointer :: station
!
	this%netcode = netcode
	this%csys = csys
	if(associated(station)) then
		this%nstat = size(station)
		this%station => station
	else
		this%nstat = 0
		nullify(this%station)
	endif
	end subroutine createLinkSeismicNetwork
!----------------------------------------------------------------
!> \brief Deallocate a seismic network object
!
	subroutine deallocSeismicNetwork(this)
	type (seismic_network) :: this
	if (associated(this%station)) deallocate(this%station)
	end subroutine deallocSeismicNetwork
!----------------------------------------------------------------
!> \brief Get coordinate system
!
	function getCoordinateSystemSeismicNetwork(this) result(res)
	type (seismic_network), intent(in) :: this
	character (len=1) :: res
	res = this%csys
	end function getCoordinateSystemSeismicNetwork
!----------------------------------------------------------------
!> \brief Get number of stations
!
	function getNstatSeismicNetwork(this) result(res)
	type (seismic_network), intent(in) :: this
	integer :: res
	res = this%nstat
	end function getNstatSeismicNetwork
!-----------------------------------------------------------------
!> \brief Get selected seismic station
!
	function getSelectedStationSeismicNetwork(this,j) result(res)
	type (seismic_network), intent(in) :: this
	integer, intent(in) :: j
	type (seismic_station) :: res
	res = this%station(j)
	end function getSelectedStationSeismicNetwork
!-----------------------------------------------------------------
!> \brief Get all seismic stations as pointer to array
!
	function getAllStationsSeismicNetwork(this) result(res)
	type (seismic_network), intent(in) :: this
	type (seismic_station), dimension(:), pointer :: res
	res => this%station
	end function getAllStationsSeismicNetwork
!---------------------------------------------------------------------
!> \brief Get station index by name and operating time
!! returns 0 in case of no match
!
	function getStationIndexByNameOptimeSeismicNetwork(this,staname,top) result(res)
	type (seismic_network) :: this
	character (len=*) :: staname
	type (date_time) :: top
	integer :: res
	type (seismic_station) :: si
	integer :: j
!
	res = 0
	do j = 1,this%nstat
		si = this%station(j)
		if(trim(adjustl(staname)) == trim(adjustl(.staname.si))) then
			if ((.opstart.si) > top .or. (.opend.si) < top) then      ! is station operating ?
				res = 0
				print *, 'Out of operating range ',staname
				cycle
			endif
			res = j
			exit
		endif		
	enddo
	end function getStationIndexByNameOptimeSeismicNetwork
!-------------------------------------------------------------------------
!> \brief Return next station
!> \param f Start with event index f
!> \param s go through events with step s
!> \param f End with event index f
!
	function nextStationSeismicNetwork(this,si,f,s,l,reset) result(next)
	type (seismic_network), intent(in) :: this
	type (seismic_station) :: si
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
	if (present(l)) then; last = l; else; last = this%nstat; endif
	call_count = call_count+1
	current = first+(call_count-1)*step
	if (current > last) then
		next = .false.
		call_count = 0
	else if (current < 1) then
		next = .false.
		call_count = 0
	else
		si = this%station(current)
		next = .true.
	endif
	end function nextStationSeismicNetwork
!--------------------------------------------------------------------------------
!> \brief Search station list for stationid, return station as output argument
!
	function searchStationNameSeismicNetwork(this,staname,station,istat) result(errmsg)
	type (seismic_network) :: this
	character (len=*) :: staname
	type (seismic_station), optional :: station
	integer, optional :: istat
	type (error_message) :: errmsg
	character (len=31) :: myname = 'searchStationNameSeismicNetwork'
	integer :: idx,j
!
	idx = 0
	call new(errmsg,myname)
	do j = 1, this%nstat
		if (adjustl(trim(.staname.(this%station(j)))) == adjustl(trim(staname))) then
			idx = j; exit
		endif
	enddo
	if (idx == 0) then
		call new(errmsg,2,'Station '//trim(staname)//' not found',myname)
		return
	endif
	if(present(station)) station = this%station(idx)
	if(present(istat)) istat = idx
	end function searchStationNameSeismicNetwork
!--------------------------------------------------------------------------------------
!> \brief Print out seismic network
!
	subroutine printSeismicNetwork(this)
	type (seismic_network) :: this
	type (seismic_station) :: si
	integer :: j
	write(6,'(3a10,3a12,2a18)') 'Station','Net','IName','Latitude','Longitude','Altitude','Opstart','Opend'
	do j = 1,this%nstat
		si= this%station(j)
		write(6,'(3a10,2f12.4,f12.2,2a18)') &
			.staname.si,.netcode.si,.inst.si,.lat.si,.lon.si,.alt.si, &
			convertToTimestringDateTime(.opstart.si),convertToTimestringDateTime(.opend.si)
	enddo
	end subroutine printSeismicNetwork
!
 end module
