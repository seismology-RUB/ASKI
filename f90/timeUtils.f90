!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 1.2.
!
!   Gemini II and ASKI version 1.2 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 1.2 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------
!> \brief Module with time utilities
!
 module timeUtils
	implicit none
!
 contains
!-----------------------------------------------------------
!> \brief Convert intelligently from two-digits year to full year
!! Works from 1951 to 2050
!
	function convertYear2ToYear4TimeUtils(yy) result(fy)
	integer :: yy,fy
	if (yy < 100) then                    ! yy is really year2
		if (yy > 50) then
			fy = yy+1900
		else
			fy = yy+2000
		endif
	else
		fy = yy
	endif
	end function convertYear2ToYear4TimeUtils
!-----------------------------------------------------------
!> \brief Convert intelligently from full year to two-digits year
!
	function convertYear4ToYear2TimeUtils(fy) result(yy)
	integer :: yy,fy
	if (fy >= 2000) then
		yy = fy-2000
	else
		yy = fy-1900
	endif
	end function convertYear4ToYear2TimeUtils
!----------------------------------------------------------
!> \brief  calculate day of year
!> \param year Year
!> \param month Month
!> \param day Day
!
	integer function getDayOfYearTimeUtils(year,month,day) result(doy)
	integer :: year,month,day,i
	integer, dimension(12) :: dm
!
	dm = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
	if (isLeapYearTimeUtils(year)) dm(2) = 29
	doy = 0
	do i=1,month-1
		doy = doy+dm(i)
	enddo
	doy = doy+day
	end function getDayOfYearTimeUtils
!-------------------------------------------------------------------
!> \brief Get number of days in year
!> \param year Year in yyyy form
!
	integer function getNdaysYearTimeUtils(year)
	integer :: year
	if (isLeapYearTimeUtils(year)) then
		getNdaysYearTimeUtils = 366
	else
		getNdaysYearTimeUtils = 365
	endif
	end function getNdaysYearTimeUtils
!---------------------------------------------------------------------------------
!> \brief Get a hhmmss string from time in seconds after midnight
!> \param t Time in seconds
!> \result hhmmss output string giving hour, minute and seconds
! 
	function hmsStringFromTsecTimeUtils(t) result(hhmmss)
	integer :: t,th
	character (len=6) :: hhmmss
	integer :: hour,minute,sec
!
	hour = t/3600; th = t-hour*3600
	minute = th/60; th = th-minute*60.
	sec = th
	write(hhmmss,'(i2.2,i2.2,i2.2)') hour,minute,sec
	end function hmsStringFromTsecTimeUtils
!------------------------------------------------------------------
!> \brief Get hh, mm, ss from time in full seconds after midnight
!> \param tfs Time in full seconds after midnight
!> \param hh Hours
!> \param mm Minutes
!> \param ss Seconds
!
	subroutine hmsFromTsecTimeUtils(tfs,hh,mm,ss)
	integer :: tfs,th,hh,mm,ss
	th = tfs
	hh = th/3600; th = th-hh*3600
	mm = th/60;
	ss = th-mm*60
	end subroutine hmsFromTsecTimeUtils
!-------------------------------------------------------------------
! > \brief Is year a leap year ?
! > \param year Year in yyyy form
!
	logical function isLeapYearTimeUtils(year)
	integer :: year
	logical :: schalt
!
	schalt = (mod(year,4) == 0)
	if (schalt) then
		if (mod(year,100) == 0 .and. mod(year,400) /= 0) schalt = .false.
	endif
	isLeapYearTimeUtils = schalt
	end function isLeapYearTimeUtils
!----------------------------------------------------------
!> \brief  calculate month and day from year and day of year
!> \param year Year
!> \param doy Day of year
!> \param month Month (output)
!> \param mday Day (output)
!
	subroutine monthDayFromDayOfYearTimeUtils(year,doy,month,mday)
	integer :: year,doy,month,mday
	integer, dimension(12), target :: daysum, daysumschalt
	integer, dimension(:), pointer :: daysump
	integer :: jl,ju,jm
	daysum = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334 /)
	daysumschalt = (/ 0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335 /)
	if (isLeapYearTimeUtils(year)) then; daysump => daysumschalt; else; daysump => daysum; endif
	jl=0
	ju=13
	do while (ju-jl > 1)
		jm=(ju+jl)/2
		if(doy > daysump(jm)) then; jl=jm; else; ju=jm; endif
	enddo
	month = jl
	mday = doy-daysump(jl)
	end subroutine monthDayFromDayOfYearTimeUtils
!------------------------------------------------------------------
!> \brief get time in seconds from hh, mm, ss
!> \param hh Hours
!> \param mm Minutes
!> \param ss Seconds
!
	integer function tsecFromHmsTimeUtils(hh,mm,ss)
	integer :: hh,mm,ss
	tsecFromHmsTimeUtils = ss+60*(mm+60*hh)
	end function tsecFromHmsTimeUtils
!
 end module
