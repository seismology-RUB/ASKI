!--------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 1.0.
!
!   Gemini II and ASKI version 1.0 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 1.0 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!------------------------------------------------------------------------
!> \brief  Module for a seismic event
!------------------------------------------------------------------------
 module seismicEvent
	use mathConstants
	use dateTime
	use errorMessage
	use flexibleType
	implicit none
	interface new
		module procedure createSeismicEvent
		module procedure createForceMomentSeismicEvent
	end interface
	interface dealloc; module procedure deallocSeismicEvent; end interface
	interface operator (.otime.); module procedure getOriginTimeSeismicEvent; end interface
	interface operator (.slat.); module procedure getLatitudeSeismicEvent; end interface
	interface operator (.slon.); module procedure getLongitudeSeismicEvent; end interface
	interface operator (.slatrad.); module procedure getRadLatitudeSeismicEvent; end interface
	interface operator (.slonrad.); module procedure getRadLongitudeSeismicEvent; end interface
	interface operator (.sdepth.); module procedure getDepthSeismicEvent; end interface
	interface operator (.styp.); module procedure getTypeSeismicEvent; end interface
	interface operator (.mag.); module procedure getMagnitudeSeismicEvent; end interface
	interface operator (.momten.); module procedure getMomentTensorSeismicEvent; end interface
	interface operator (.force.); module procedure getForceSeismicEvent; end interface
	interface operator (.evid.); module procedure getEventIdSeismicEvent; end interface
	interface operator (.csys.); module procedure getCoordinateSystemSeismicEvent; end interface
	interface operator (.year.); module procedure getYearSeismicEvent; end interface
	interface operator (.doy.); module procedure getDoySeismicEvent; end interface
	integer, parameter :: character_length_evid = 13
	type seismic_event
		private
		integer :: istyp                                   ! Source type: 0 = force, 1 = moment tensor, -1: not specified
		character (len=character_length_evid) :: evid      ! Event ID (e.g. 2006.10.2977 or 061113_141238)
		character (len=1) :: csys                          ! Character code for coordinate system (C,S)
		real :: slat,slon,sdepth                           ! Lat and lon and depth of event (in degrees/km) or x,y,z in meters
		real :: mag                                        ! Magnitude of event
		type (date_time) :: otime                          ! origin time
		real, dimension(:), pointer :: momten => null()    ! moment tensor (depending on istyp)
		real, dimension(:), pointer :: force => null()     ! force (depending on istyp)
	end type
!
 contains
!--------------------------------------------------------------------------------------
!> \brief Create a seismic event object
!
	subroutine createSeismicEvent(this,evid,csys,slat,slon,sdepth,mag,otime)
	type (seismic_event) :: this
	character (len=*) :: evid,csys
	real :: slat,slon,sdepth,mag
	type (date_time) :: otime
!
	this%evid = evid; this%csys = csys; this%slat = slat; this%slon = slon
	this%sdepth = sdepth; this%mag = mag; this%otime = otime
	this%istyp = -1
	end subroutine createSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Create a seismic event object with force/moment information
!
	subroutine createForceMomentSeismicEvent(this,evid,csys,slat,slon,sdepth,mag,otime,istyp,momfor)
	type (seismic_event) :: this
	character (len=*) :: evid,csys
	real :: slat,slon,sdepth,mag
	type (date_time) :: otime
	integer :: istyp
	real, dimension(:) :: momfor
!
	this%evid = evid; this%csys = csys; this%slat = slat; this%slon = slon
	this%sdepth = sdepth; this%mag = mag; this%otime = otime
	if (istyp == 0) then
		allocate(this%force(3)); this%force = momfor
	else if(istyp == 1) then
		allocate(this%momten(6)); this%momten = momfor
	endif
	this%istyp = istyp
	end subroutine createForceMomentSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Create seismic event from flexible vector
!! \details given a vector of type flexible, this routine creates a seismic event object
!!  recovering its members from the flexible type vector in the same manner as the members
!!  are returned by function getMembersSeismicEvent
!
	subroutine createFromFlexibleSeismicEvent(this,members)
	type (seismic_event) :: this
	type (flexible), dimension(:) :: members
	character(len=25) :: string
	integer :: i
	call deallocSeismicEvent(this)
	this%evid = members(1)
	this%csys = members(2)
	this%slat = members(3)
	this%slon = members(4)
	this%sdepth = members(5)
	this%mag = members(6)
	string = members(7)
	call createFromFullTimestringDateTime(this%otime,string)
	this%istyp = members(8)
	if (this%istyp == 0) then
		allocate(this%force(3))
		do i=1,3; this%force(i) = members(8+i); enddo
	elseif (this%istyp == 1) then
		allocate(this%momten(6))
		do i=1,6; this%momten(i) = members(8+i); enddo
	endif
	end subroutine createFromFlexibleSeismicEvent
!---------------------------------------------------------------------------
!> \brief True copy of a seismic event object
!
	subroutine copySeismicEvent(this,copy)
	type (seismic_event) :: this,copy
	copy = this
	if (associated(this%force)) then
		allocate(copy%force(3)); copy%force = this%force
	endif
	if (associated(this%momten)) then
		allocate(copy%momten(6)); copy%momten = this%momten
	endif
	end subroutine copySeismicEvent
!---------------------------------------------------------------------------
!> \brief Deallocate seismic event object
!
	subroutine deallocSeismicEvent(this)
	type (seismic_event) :: this
	if (associated(this%momten)) deallocate(this%momten)
	if (associated(this%force)) deallocate(this%force)
	end subroutine deallocSeismicEvent
!---------------------------------------------------------------------------
!> \brief Extend a pointer array of seismic event objects
!! If object contains pointers to array, the copied ones should be unlinked
!! and the no more used ones deallocated (example miniSeed.f90)
!
	function extendArraySeismicEvent(array,n) result(newarray)
	type (seismic_event), dimension(:), pointer :: array
	type (seismic_event), dimension(:), pointer :: newarray
	integer :: n,nold,i
!
	allocate(newarray(n))
	if (.not. associated(array)) return
	nold = min(size(array),n)
	newarray(1:nold) = array(1:nold)
	do i = 1,nold
		call unlinkSeismicEvent(array(i))
	enddo
	do i = nold+1,size(array)
		call deallocSeismicEvent(array(i))
	enddo
	deallocate(array)
	end function extendArraySeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get cartesian coordinates of epicenter on unit sphere
!
	subroutine getCartesianEpicenterSeismicEvent(this,x,y,z)
	type (seismic_event) :: this
	real :: x,y,z
	x = cos(mc_deg2rad*(90.-this%slat))*cos(mc_deg2rad*this%slon)
	y = cos(mc_deg2rad*(90.-this%slat))*sin(mc_deg2rad*this%slon)
	z = cos(mc_deg2rad*(90.-this%slat))
	end subroutine getCartesianEpicenterSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get event id
!
	function getCoordinateSystemSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	character (len=1) :: res
	res = this%csys
	end function getCoordinateSystemSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get depth of event
!
	function getDepthSeismicEvent(this) result(d)
	type (seismic_event), intent(in) :: this
	real :: d
	d = this%sdepth
	end function getDepthSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get doy of event
!
	function getDoySeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	integer :: res
	res = .doy.(this%otime)
	end function getDoySeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get event id
!
	function getEventIdSeismicEvent(this) result(evid)
	type (seismic_event), intent(in) :: this
	character (len=character_length_evid) :: evid
	evid = this%evid
	end function getEventIdSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get force of event
!
	function getForceSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	real, dimension(:), pointer :: res
	res => this%force
	end function getForceSeismicEvent
!-----------------------------------------------------------------------------------
!> \brief Create a hypocenter line from event information
!
	function getHypolineSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	character (len=146) :: res
	character (len=25) :: tm
	integer :: yy,mo,dd,hh,mm,is,ns,i
	real :: ss
!
	tm = .timestring.(this%otime)
	read(tm,'(i4,i2,i2,1x,i2,i2,i2,1x,i9)') yy,mo,dd,hh,mm,is,ns
	ss = is+1.e-9*ns
	if (this%istyp < 0 .or. this%istyp > 1) then
		write(res,'(a13,1x,i4.4,4(1x,i2.2),f6.2,f9.4,f10.4,f8.3,f5.1)') &
		    & this%evid,yy,mo,dd,hh,mm,ss,this%slat,this%slon, &
		    & this%sdepth,this%mag
	else if (this%istyp == 0) then                    ! force
		write(res,'(a13,1x,i5,1x,i4.4,4(1x,i2.2),f6.2,f9.4,f10.4,f8.3,f5.1,3e12.3)') &
		    & this%evid,this%istyp,yy,mo,dd,hh,mm,ss,this%slat,this%slon, &
		    & this%sdepth,this%mag,(this%force(i),i=1,3)
	else if (this%istyp == 1) then                    ! moment
		write(res,'(a13,1x,i5,1x,i4.4,4(1x,i2.2),f6.2,f9.4,f10.4,f8.3,f5.1,6e12.3)') &
		    & this%evid,this%istyp,yy,mo,dd,hh,mm,ss,this%slat,this%slon, &
		    & this%sdepth,this%mag,(this%momten(i),i=1,6)
	endif
	end function
!--------------------------------------------------------------------------------------
!> \brief Get latitude of event
!
	function getLatitudeSeismicEvent(this) result(lat)
	type (seismic_event), intent(in) :: this
	real :: lat
	lat = this%slat
	end function getLatitudeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get latitude of event in rad
!
	function getRadLatitudeSeismicEvent(this) result(lat)
	type (seismic_event), intent(in) :: this
	real :: lat
	lat = this%slat*mc_deg2rad
	end function getRadLatitudeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get longitude of event
!
	function getLongitudeSeismicEvent(this) result(lon)
	type (seismic_event), intent(in) :: this
	real :: lon
	lon = this%slon
	end function getLongitudeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get longitude of event in rad
!
	function getRadLongitudeSeismicEvent(this) result(lon)
	type (seismic_event), intent(in) :: this
	real :: lon
	lon = this%slon*mc_deg2rad
	end function getRadLongitudeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get magnitude of event
!
	function getMagnitudeSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	real :: res
	res = this%mag
	end function getMagnitudeSeismicEvent
!-------------------------------------------------------------------------------------
!> \brief Get members of event type
!
	function getMembersSeismicEvent(this) result(res)
	type (seismic_event) :: this
	type (flexible), dimension(:), pointer :: res
	integer :: i
!
	if (this%istyp < 0 .or. this%istyp > 1) then
		allocate(res(8))
	else if (this%istyp == 0) then
		allocate(res(11))
	else if (this%istyp == 1) then
		allocate(res(14))
	endif
	res(1) = trim(this%evid)
	res(2) = this%csys
	res(3) = this%slat
	res(4) = this%slon
	res(5) = this%sdepth
	res(6) = this%mag
	res(7) = convertToFullTimestringDateTime(this%otime)
	res(8) = this%istyp
	if (this%istyp == 0) then
		do i = 1,3; res(8+i) = this%force(i); enddo
	else if (this%istyp == 1) then
		do i = 1,6; res(8+i) = this%momten(i); enddo
	endif
	end function getMembersSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get moment tensor of event
!
	function getMomentTensorSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	real, dimension(:), pointer :: res
	res => this%momten
	end function getMomentTensorSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get origin time
!
	function getOriginTimeSeismicEvent(this) result(otime)
	type (seismic_event), intent(in) :: this
	type (date_time) :: otime
	otime = this%otime
	end function getOriginTimeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get type of event
!
	function getTypeSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	integer :: res
	res = this%istyp
	end function getTypeSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Get year of event
!
	function getYearSeismicEvent(this) result(res)
	type (seismic_event), intent(in) :: this
	integer :: res
	res = .year.(this%otime)
	end function getYearSeismicEvent
!---------------------------------------------------------------------------
!> \brief Is moment tensor available ?
!
	logical function isMomentSeismicEvent(this)
	type (seismic_event) :: this
	isMomentSeismicEvent = (this%istyp == 1)
	end function isMomentSeismicEvent
!----------------------------------------------------------------------------
!> \brief Print seismic event
!
	subroutine printSeismicEvent(this)
	type (seismic_event) :: this
!
	if (this%istyp  < 0 .or. this%istyp > 1) then
		write(6,'(a13,2x,a4,2x,4f12.3,2x,a25)') trim(this%evid),this%csys,this%slat,this%slon,this%sdepth,this%mag, &
	                                               convertToFullTimestringDateTime(this%otime)
	else if (this%istyp == 0) then
		write(6,'(a13,2x,a4,2x,4f12.3,2x,a25,3e15.3)') &
			trim(this%evid),this%csys,this%slat,this%slon,this%sdepth,this%mag, &
                     convertToFullTimestringDateTime(this%otime),this%force
	else if (this%istyp == 1) then
		write(6,'(a13,2x,a4,2x,4f12.3,2x,a25,6e15.5)') &
			trim(this%evid),this%csys,this%slat,this%slon,this%sdepth,this%mag, &
                     convertToFullTimestringDateTime(this%otime),this%momten
	endif
	end subroutine printSeismicEvent
!-------------------------------------------------------------
!> \brief set new epicenter coordinates
!
	subroutine setEpicenterSeismicEvent(this,slat,slon)
	type (seismic_event) :: this
	real :: slat,slon
	this%slat = slat; this%slon = slon
	end subroutine setEpicenterSeismicEvent
!--------------------------------------------------------------------------------------
!> \brief Nullify pointers
!
	subroutine unlinkSeismicEvent(this)
	type (seismic_event) :: this
	if (associated(this%momten)) nullify(this%momten)
	if (associated(this%force)) nullify(this%force)
	end subroutine unlinkSeismicEvent
!------------------------------------------------------------------------------
 end module
