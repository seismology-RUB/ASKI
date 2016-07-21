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
!----------------------------------------------------------------
!> \brief Module with time utilities
!
 module dateTime
	use timeUtils
	use realloc
	implicit none
	interface new
		module procedure createDateTime
		module procedure createYMDDateTime
		module procedure createFromFullTimestringDateTime
	end interface
	interface addDateTime
		module procedure addFsNsDateTime
		module procedure addDPDateTime
	end interface
	interface subtractDateTime
		module procedure subtractFsNsDateTime
		module procedure subtractDPDateTime
	end interface
	interface operator (>); module procedure isLaterDateTime; end interface 
	interface operator (<); module procedure isEarlierDateTime; end interface 
	interface operator (==); module procedure isEqualDateTime; end interface
	interface operator (-); module procedure differenceDateTime; end interface
	interface operator (.year.); module procedure getYearDateTime; end interface
	interface operator (.month.); module procedure getMonthDateTime; end interface
	interface operator (.day.); module procedure getDayDateTime; end interface
	interface operator (.doy.); module procedure getDoyDateTime; end interface
	interface operator (.tfs.); module procedure getTfsDateTime; end interface
	interface operator (.tns.); module procedure getTnsDateTime; end interface
	interface operator (.tanfdp.); module procedure getTanfDPDateTime; end interface
	interface operator (.timestring.); module procedure convertToFullTimestringDateTime; end interface
	interface operator (.plus.); module procedure addDPDateTime; end interface
	interface operator (.minus.); module procedure subtractDPDateTime; end interface
	interface operator (.sdate.); module procedure convertToSdateDateTime; end interface
	interface operator (.eod.); module procedure endOfDayDateTime; end interface
	type date_time
		private
		integer :: year     !  year (yyyy)
		integer :: doy      !  day of year
		integer :: hh       !  hour
		integer :: mm       !  minute
		integer :: ss       !  seconds
		integer :: ns       !  nanoseconds
	end type
!
 contains
!---------------------------------------------------------------------------------
!> \brief Direct constructor
!
	subroutine createDateTime(this,year,doy,hh,mm,ss,ns)
	type (date_time) :: this
	integer :: year,doy,hh,mm,ss,ns
	this = date_time(year,doy,hh,mm,ss,ns)
	end subroutine createDateTime
!---------------------------------------------------------------------------------
!> \brief Constructor from year, month, day etc
!
	subroutine createYMDDateTime(this,year,month,day,hh,mm,ss,ns)
	type (date_time) :: this
	integer :: year,month,day,hh,mm,ss,ns,doy
	doy = getDayOfYearTimeUtils(year,month,day)
	this = date_time(year,doy,hh,mm,ss,ns)
	end subroutine createYMDDateTime
!---------------------------------------------------------------------------------
!> \brief Constructor from yyyymmdd_hhmmss_nnnnnnnnn string
!
	subroutine createFromFullTimestringDateTime(this,timestring)
	type (date_time) :: this
	character (len=*) :: timestring
	integer :: year,month,day,hh,mm,ss,doy,ns
	if (len_trim(timestring) < 25) then
		ns = 0
		read(timestring,'(i4,i2,i2,1x,i2,i2,i2)') year,month,day,hh,mm,ss
	else
		read(timestring,'(i4,i2,i2,1x,i2,i2,i2,1x,i9)') year,month,day,hh,mm,ss,ns
	endif
	doy = getDayOfYearTimeUtils(year,month,day)
	this = date_time(year,doy,hh,mm,ss,ns)
	end subroutine createFromFullTimestringDateTime
!----------------------------------------------------------------------
!> \brief Add some time interval given as double precision
!> \param this Date_time object
!> \param dtdb double precision change in time 
!
	function addDPDateTime(this,dbdt) result(newtm)
	type (date_time), intent(in) :: this
	double precision, intent(in) :: dbdt
	type (date_time) :: newtm
	integer :: dt,dtns
	dt = int(dbdt)
	dtns = int(1000000000*(dbdt-dt))
	newtm = changeDateTime(this,dt,dtns,'add')
	end function addDPDateTime
!----------------------------------------------------------------------
!> \brief Add some time interval (wrapper for changeDateTime)
!> \param this Date_time object
!> \param dt change in time in full seconds
!> \param dtns change in time in full nanoseconds
!
	function addFsNsDateTime(this,dt,dtns) result(newtm)
	type (date_time) :: this,newtm
	integer :: dt,dtns
	newtm = changeDateTime(this,dt,dtns,'add')
	end function addFsNsDateTime
!--------------------------------------------------------------------
!> \brief  build a yyyy/doy string
!
	function buildYearDoyStringDateTime(this) result(string)
	type (date_time) :: this
	character (len=8) :: string
	write(string,'(i4.4,a1,i3.3)') this%year,'/',this%doy
	end function buildYearDoyStringDateTime
!---------------------------------------------------------------------------------
!> \brief Modify time by some amount of seconds and nanoseconds (either positive or negative)
!> \param this Date_time object
!> \param dt change in time in full seconds
!> \param dtns change in time in full nanoseconds
!> \param action either 'add' or 'sub'
!
	function changeDateTime(this,dt,dtns,action) result(newtm)
	type (date_time) :: this,newtm
	integer :: dt,dtns
	integer :: year,doy,hh,mm,ss,dth,day,hour,minute,sec,ndy,vorz=0,dthns,ns
	character (len=3) :: action
!
	dth = iabs(dt)
	dthns = iabs(dtns)
	if (dth == 0 .and. dthns == 0) then
		newtm = copyDateTime(this)
		return
	endif
!
	select case (action)
	case ('sub'); vorz = -1
	case ('add'); vorz = +1
	case default
		print *,'changeDateTime: unknown action'
		stop
	end select

	day = dth/86400; dth = dth-day*86400
	hour = dth/3600; dth = dth-hour*3600
	minute = dth/60
	sec = dth-minute*60
	year = this%year
	doy = this%doy+vorz*day
	hh = this%hh+vorz*hour
	mm = this%mm+vorz*minute
	ss = this%ss+vorz*sec
	ns = this%ns+vorz*dthns
	if (action == 'sub') then
		if (ns < 0) then; ns = ns+1000000000; ss = ss-1; endif
		if (ss < 0) then; ss = ss+60; mm = mm-1; endif	
		if (mm < 0) then; mm = mm+60; hh = hh-1; endif
		if (hh < 0) then; hh = hh+24; doy = doy-1; endif
		do while(doy <= 0)
			year = year-1
			ndy = getNdaysYearTimeUtils(year)
			doy = doy+ndy
		enddo
		call createDateTime(newtm,year,doy,hh,mm,ss,ns)
	else if (action == 'add') then
		if (ns > 999999999) then; ns = ns-1000000000; ss = ss+1; endif
		if (ss > 59) then; ss = ss-60; mm = mm+1; endif
		if (mm > 59) then; mm = mm-60; hh = hh+1; endif
		if (hh > 23) then; hh = hh-24; doy = doy+1; endif
		ndy = getNdaysYearTimeUtils(year)
		do while(doy > ndy)
			doy = doy-ndy
			year = year+1
			ndy = getNdaysYearTimeUtils(year)
		enddo
		call createDateTime(newtm,year,doy,hh,mm,ss,ns)
	endif
	end function changeDateTime
!--------------------------------------------------------------------
!> \brief  return a yyyymmdd_hhmmss time string
!
	function convertToTimestringDateTime(this) result(timestring)
	type (date_time), intent(in) :: this
	character (len=15) :: timestring
	integer :: month,day
!
	call monthDayFromDayOfYearTimeUtils(this%year,this%doy,month,day)
	write(timestring,'(i4.4,2i2.2,a1,3i2.2)') this%year,month,day,'_',this%hh,this%mm,this%ss
	end function convertToTimestringDateTime
!--------------------------------------------------------------------
!> \brief return a yyyymmdd_hhmmss_nnnnnnnnn time string
!
	function convertToFullTimestringDateTime(this) result(timestring)
	type (date_time), intent(in) :: this
	character (len=25) :: timestring
	integer :: month,day
!
	call monthDayFromDayOfYearTimeUtils(this%year,this%doy,month,day)
	write(timestring,'(i4.4,2i2.2,a1,3i2.2,a1,i9.9)') this%year,month,day,'_',this%hh,this%mm,this%ss,'_',this%ns
	end function convertToFullTimestringDateTime
!--------------------------------------------------------------------
!> \brief  return a yyyymmdd (sdate) string
!
	function convertToSdateDateTime(this) result(sdate)
	type (date_time), intent(in) :: this
	character (len=8) :: sdate
	integer :: month,day
!
	call monthDayFromDayOfYearTimeUtils(this%year,this%doy,month,day)
	write(sdate,'(i4.4,2i2.2)') this%year,month,day
	end function convertToSdateDateTime
!--------------------------------------------------------------------
!> \brief  convert to a SEED time string (yyyy,ddd,hh,mm,ss.ffffff)
!
	function convertToSeedTimestringDateTime(this) result(timestring)
	type (date_time), intent(in) :: this
	character (len=24) :: timestring
	write(timestring,'(i4.4,a1,i3.3,a1,i2.2,a1,i2.2,a1,i2.2,a1,i6.6)') &
	& this%year,',',this%doy,',',this%hh,',',this%mm,',',this%ss,'.',nint(this%ns/1000.)
	end function convertToSeedTimestringDateTime
!--------------------------------------------------------------------
!> \brief copy a date_time object
!
	function copyDateTime(this) result(newtm)
	type (date_time) :: this,newtm
	newtm = this
	end function copyDateTime
!---------------------------------------------------------------------
!> \brief Absolute time difference expressed by a date_time object
!> \par 
!! Warning: Do not use year here, always keep it at zero, largest unit is days
!<
	function differenceDateTime(this,tm) result(tdiff)
	type (date_time), intent(in) :: this,tm
	type (date_time) :: tdiff
	integer :: vorz=0,ns,ss,mm,hh,doy,year,ndy,dyear
!
	if (this == tm) then; call createDateTime(tdiff,0,0,0,0,0,0); return; endif
!
	if (this > tm) then; vorz = +1; else; vorz = -1; endif
	ns = vorz*(this%ns-tm%ns)
	ss = vorz*(this%ss-tm%ss)
	mm = vorz*(this%mm-tm%mm)
	hh = vorz*(this%hh-tm%hh)
	doy = vorz*(this%doy-tm%doy)
	year = max(this%year,tm%year)
	dyear = vorz*(this%year-tm%year)
	if (ns < 0) then; ns = ns+1000000000; ss = ss-1; endif
	if (ss < 0) then; ss = ss+60; mm = mm-1; endif	
	if (mm < 0) then; mm = mm+60; hh = hh-1; endif
	if (hh < 0) then; hh = hh+24; doy = doy-1; endif
	do while(dyear > 0)
		year = year-1
		dyear = dyear-1
		ndy = getNdaysYearTimeUtils(year)
		doy = doy+ndy
	enddo
	call createDateTime(tdiff,dyear,doy,hh,mm,ss,ns)
	end function differenceDateTime
!-------------------------------------------------
!> \brief Get end of day as date_time object
!
	function endOfDayDateTime(this) result(eod)
	type (date_time), intent(in) :: this
	type (date_time) :: eod
	eod = date_time(this%year,this%doy,23,59,59,999999999)
	end function endOfDayDateTime
!---------------------------------------------------------------------------
!> \brief Extend a pointer array of date_time objects
!! If object contains pointers to array, the copied ones should be unlinked
!! and the no more used ones deallocated (example miniSeed.f90)
!
	function extendArrayDateTime(array,n) result(newarray)
	type (date_time), dimension(:), pointer :: array
	type (date_time), dimension(:), pointer :: newarray
	integer :: n,nold
!
	allocate(newarray(n))
	if (.not. associated(array)) return
	nold = min(size(array),n)
	newarray(1:nold) = array(1:nold)
	deallocate(array)
	end function extendArrayDateTime
!-------------------------------------------------------------------
!> \brief Get year
!
	integer function getYearDateTime(this)
	type (date_time), intent(in) :: this
	getYearDateTime = this%year
	end function getYearDateTime
!--------------------------------------------------------------------
!> \brief  get month
!
	integer function getMonthDateTime(this)
	type (date_time), intent(in) :: this
	integer :: month,day
!
	call monthDayFromDayOfYearTimeUtils(this%year,this%doy,month,day)
	getMonthDateTime = month
	end function getMonthDateTime
!--------------------------------------------------------------------
!> \brief  get mday
!
	integer function getDayDateTime(this)
	type (date_time), intent(in) :: this
	integer :: month,day
!
	call monthDayFromDayOfYearTimeUtils(this%year,this%doy,month,day)
	getDayDateTime = day
	end function getDayDateTime
!----------------------------------------------------------
!> \brief Get day of year
!
	integer function getDoyDateTime(this)
	type (date_time), intent(in) :: this
	getDoyDateTime = this%doy
	end function getDoyDateTime
!----------------------------------------------------------
!> \brief Get members of date time object in integer array
!
	function getMembersDateTime(this) result(res)
	type (date_time) :: this
	integer, dimension(:), pointer :: res
	allocate(res(6))
	res = (/ this%year,this%doy,this%hh,this%mm,this%ss,this%ns /)
	end function getMembersDateTime
!---------------------------------------------------
!> \brief Get seconds after midnight in double precision
!
	double precision function getTanfDPDateTime(this)
	type (date_time), intent(in) :: this
	getTanfDPDateTime = this%hh*3600.d0+this%mm*60.d0+dble(this%ss)+1.d-9*dble(this%ns)
	end function getTanfDPDateTime
!---------------------------------------------------
!> \brief Get time in full seconds after midnight
!
	integer function getTfsDateTime(this)
	type (date_time), intent(in) :: this
	getTfsDateTime = tsecFromHmsTimeUtils(this%hh,this%mm,this%ss)
	end function getTfsDateTime
!---------------------------------------------------
!> \brief Get nanoseconds part of time
!
	integer function getTnsDateTime(this)
	type (date_time), intent(in) :: this
	getTnsDateTime = this%ns
	end function getTnsDateTime
!---------------------------------------------------------------------------------
!> \brief Calculate number of intervals of length dt contained between this and tm
!> \param this,tm date_time objects
!> \param dt Interval length
!> \param err Absolute fitting error in seconds
!> \return Number of intervals of length \a dt
!
	integer function intervalsBetweenDateTime(this,tm,dt,err)
	type (date_time) :: this,tm,tdiff
	double precision :: dt,ssd,nsd,err
	integer :: ks,kns,kr
	double precision :: rs,rns,rest
!
	tdiff = differenceDateTime(this,tm)
	ssd = tdiff%ss+60.d0*(tdiff%mm+60.d0*(tdiff%hh+24.d0*tdiff%doy))
	nsd = dble(tdiff%ns)
	if (ssd/dt > huge(ks)) then
		print *,'Number of intervals too large to be represented by normal integer'
		stop
	endif
!
!  number of intervals
!
	ks = int(ssd/dt)                       ! sampling intervals fitting into ss-part
	rs = ssd-ks*dt                         ! remainder
	kns = int(nsd*1.d-9/dt)                ! sampling intervals fitting into ns-part
	rns = nsd*1.d-9-kns*dt                 ! remainder in seconds
	rest = rs+rns                          ! total remainder from ss and ns parts
	kr = nint(rest/dt)                     ! sampling intervals fitting into total remainder
	err = abs(rest-dble(kr)*dt)            ! remainder/(sampling interval) to estimate error
	intervalsBetweenDateTime = ks+kns+kr
	end function intervalsBetweenDateTime
!---------------------------------------------------------------------
!> \brief Is this earlier than comp ?
!
	function isEarlierDateTime(this,comp) result(earlier)
	type (date_time), intent(in) :: this,comp
	logical :: earlier
!
	if (isEqualDateTime(this,comp)) then; earlier = .false.; return; endif
!
	if (this%year < comp%year) then; earlier = .true.
	else if (this%year > comp%year) then; earlier = .false.
	else
		if (this%doy < comp%doy) then; earlier = .true.
		else if (this%doy > comp%doy) then; earlier = .false.
		else
			if (this%hh < comp%hh) then; earlier = .true.
			else if (this%hh > comp%hh) then; earlier = .false.
			else
				if (this%mm < comp%mm) then; earlier = .true.
				else if (this%mm > comp%mm) then; earlier = .false.
				else
					if (this%ss < comp%ss) then; earlier = .true.
					else if (this%ss > comp%ss) then; earlier = .false.
					else
						if (this%ns < comp%ns) then; earlier = .true.
						else if (this%ns >= comp%ns) then; earlier = .false.
						endif
					endif
				endif
			endif
		endif
	endif
	end function isEarlierDateTime
!---------------------------------------------------------------------
!> \brief Is this later than comp ?
!
	function isLaterDateTime(this,comp) result(later)
	type (date_time), intent(in) :: this,comp
	logical :: later
!
	later = .true.
	if (isEqualDateTime(this,comp)) then; later = .false.; endif
	if (isEarlierDateTime(this,comp)) then; later = .false.; endif
	end function isLaterDateTime
!---------------------------------------------------------------------
!> \brief Is comp equal to this ?
!
	function isEqualDateTime(this,comp) result(equal)
	type (date_time), intent(in) :: this,comp
	logical :: equal
!
	equal = (comp%year == this%year .and. comp%doy == this%doy &
	       & .and. comp%hh == this%hh .and. comp%mm == this%mm &
	       & .and. comp%ss == this%ss .and. comp%ns == this%ns)
	end function isEqualDateTime
!--------------------------------------------------------------------
!> \brief  print a data_time object
!
	subroutine printDateTime(this)
	type (date_time) :: this
	write(*,'(i4.4,i3.3,a1,3i2.2,a1,i9.9)') this%year,this%doy,'_',this%hh,this%mm,this%ss,'_',this%ns
	end subroutine printDateTime
!----------------------------------------------------------------------
!> \brief Subtract some time interval (wrapper for changeDateTime)
!> \param this Date_time object
!> \param dt change in time in full seconds
!> \param dtns change in time in full nanoseconds
!
	function subtractFsNsDateTime(this,dt,dtns) result(newtm)
	type (date_time) :: this,newtm
	integer :: dt,dtns
	newtm = changeDateTime(this,dt,dtns,'sub')
	end function subtractFsNsDateTime
!----------------------------------------------------------------------
!> \brief Subtract some time interval given as double precision
!> \param this Date_time object
!> \param dtdb double precision change in time 
!
	function subtractDPDateTime(this,dbdt) result(newtm)
	type (date_time), intent(in) :: this
	double precision, intent(in) :: dbdt
	type (date_time) :: newtm
	integer :: dt,dtns
	dt = int(dbdt)
	dtns = int(1000000000*(dbdt-dt))
	newtm = changeDateTime(this,dt,dtns,'sub')
	end function subtractDPDateTime
!
 end module dateTime
