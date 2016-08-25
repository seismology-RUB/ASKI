!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------
!> \brief Module providing functions for string handling
!--------------------------------------------------------
module string
    use errorMessage
    implicit none
    interface operator(+)
       module procedure concatenateString
    end interface
    interface operator(.equal.)
       module procedure equalString
    end interface
    integer, parameter :: max_length_string = 400
    integer, parameter :: max_length_sep_string = 10
!
contains
!-------------------------------------------------------
!> \brief Count words contained in a string
!! Spliting is done where a "separator string" appears.
!! Leading and trailing blanks are always ignored.
!! Separator string is taken "as is".
!
    function countWordsString(s,sps,errmsg) result(res)
    character(len=*) :: s
    character(len=*) :: sps
    integer :: res
    type (error_message) :: errmsg
    character(len=17) :: myname = 'countWordsString' 
    integer :: lsp,nlb,nsc,icrel,ic
    character(len=max_length_string) :: errstr
!
    res = 0
    call addTrace(errmsg,myname)
    if (len(sps) > max_length_sep_string) then
       write(errstr,*) 'Separator string is longer than ',max_length_sep_string,' characters'
       call add(errmsg,2,errstr,myname)
       call print(errmsg); return
    endif
    if (len_trim(adjustl(s)) == 0) then
       call add(errmsg,2,'Input string either all blank or empty',myname)
       call print(errmsg); return
    endif
!
!  count words LEFT of separator string
!  add single word right of last separator string (if present) later
!
    lsp = len(sps)
    nlb = len_trim(s)-len_trim(adjustl(s))    ! # of leading blanks
    nsc = 0
    icrel = index(trim(s(nlb+1:)),sps)
    ic = icrel+nlb
    do while (icrel /= 0 .and. ic < len_trim(s))  
       if (icrel > 1) nsc = nsc+1                ! do not count if sep string is right at beginning
       icrel = index(trim(s(ic+lsp:)),sps)       ! search behind word+sps
       ic = ic+lsp+icrel-1
    enddo
    res = nsc+1
    end function countWordsString
!-------------------------------------------------------
!> \brief  Split a character string into multiple words
!! Splitting is done where a "separator string" appears.
!! Leading and trailing blanks are always ignored.
!! Returned words are stripped of leading blanks.
!! Separator string is taken "as is".
!> \param s input string
!> \param c split string
!> \param errmsg error message object
!> \return allocated pointer array with words
!
    function getWordsString(s,sps,errmsg) result(res)
    character(len=*) :: s
    character(len=*) :: sps
    character(len=max_length_string), dimension(:), pointer :: res
    type (error_message) :: errmsg
    character(len=15) :: myname = 'getWordsString' 
    integer :: ic,icrel,iwb,nsc,nword,lsp,nlb
!
    res => null()
    call addTrace(errmsg,myname)
    nword = countWordsString(s,sps,errmsg)
    if (.level.errmsg == 2) return
    allocate(res(nword))
!
!  now set words
!
    lsp = len(sps)
    nlb = len_trim(s)-len_trim(adjustl(s))    ! # of leading blanks
    iwb = nlb+1                               ! index of word begin
    icrel = index(trim(s(nlb+1:)),sps)
    ic = icrel+nlb
    nsc = 0
    do while (icrel /= 0)
       if (icrel > max_length_string+1) then
          call add(errmsg,2,'there are words longer than max_length_string',myname)
          call print(errmsg); return
       endif
       if (icrel > 1) then
          nsc = nsc+1
          res(nsc) = adjustl(s(iwb:ic-1))
       endif
       icrel = index(trim(s(ic+lsp:)),sps)
       iwb = ic+lsp                           ! next word starts here
       ic = ic+lsp+icrel-1
    enddo
    if (nsc < nword) res(nsc+1) = adjustl(s(iwb:len_trim(s))) 
    end function getWordsString
!-----------------------------------------------------------
!> \brief Append string
!! Input string is modified
!! Trailing blanks are deleted
!! Leading blanks are kept
!
    subroutine appendString(s,a,errmsg)
    character(len=*) :: s
    character(len=*) :: a
    type (error_message) :: errmsg
    character(len=max_length_string) :: errstr
    character(len=12) :: myname = 'appendString'
!
    if (len_trim(adjustl(s))+len_trim(a) > max_length_string) then
       write(errstr,*) "Concatenated string's length exceeds ",max_length_string," characters."
       call add(errmsg,2,errstr,myname)
       call print(errmsg); return
    endif
    if (len_trim(adjustl(s)) == 0) then
       call add(errmsg,1,'Input string either all blank or empty',myname)
       call print(errmsg); return
    endif
    if (len_trim(adjustl(a)) == 0) then
       call add(errmsg,1,'Append string either all blank or empty',myname)
       call print(errmsg); return
    endif
!
    s = trim(s)//trim(a)
    end subroutine appendString
!-----------------------------------------------------------
!> \brief Concatenate two strings
!! Trailing blanks of string 1 are deleted
!! Leading blanks of both strings are kept.
!! Return of error meassge is omitted here because requirement for intent(in)
!> \return concatenated string
!
    function concatenateString(s,a) result(res)
    character(len=*), intent(in) :: s
    character(len=*), intent(in) :: a
    character(len=max_length_string) :: res
!
    if (len_trim(adjustl(s))+len_trim(a) > max_length_string) then
       print *, "Concatenated string's length exceeds ",max_length_string," characters."
       return
    endif
!
    res = trim(s)//trim(a)
    end function concatenateString
!--------------------------------------------------------------
!> \brief Are two strings equal
!
    function equalString(s1,s2) result(res)
    character(len=*), intent(in) :: s1,s2
    logical :: res
    res = (trim(adjustl(s1)) == trim(adjustl(s2)))
    end function equalString
!--------------------------------------------------------------
!> \brief Get string ahead of last separator
!! Eliminate leading blanks from result
!! Separator string is taken "as is"
!! Optionally include separator string at end of result
!
    function getAheadLastSeparatorString(s,sps,errmsg,addsepatend) result(res)
    character(len=*) :: s,sps
    type (error_message) :: errmsg
    logical, optional :: addsepatend
    character(len=max_length_string) :: res
    character(len=27) :: myname = 'getAheadLastSeparatorString'
    character(len=400) :: errstr
    integer :: nlb,icrel,iwe
!
    call addTrace(errmsg,myname)
    if (len(sps) > max_length_sep_string) then
       write(errstr,*) 'Separator string is longer than ',max_length_sep_string,' characters'
       call add(errmsg,2,errstr,myname)
       call print(errmsg); return
    endif
    if (len_trim(adjustl(s)) == 0) then
       call add(errmsg,2,'Input string either all blank or empty',myname)
       call print(errmsg); return
    endif
!
    nlb = len_trim(s)-len_trim(adjustl(s))       ! # of leading blanks
    icrel = index(trim(s(nlb+1:)),sps,.true.)    ! search last occurence of sps
    if (icrel < 1) then
       call add(errmsg,1,'Result string is empty, return empty string',myname)
       res = ''; return                          ! return empty string
    endif
    if (present(addsepatend)) then
       iwe = icrel+nlb                              ! end of word
    else
       iwe = icrel+nlb-1
    endif
    if (iwe-nlb > max_length_string) then
       call add(errmsg,1,'Result string is longer than max_length_string, return truncated string',myname)
    endif
    res = s(nlb+1:min(iwe,max_length_string+nlb))
    end function getAheadLastSeparatorString
!--------------------------------------------------------------
!> \brief Get string behind last separator
!! Eliminate leading blanks from result
!! Separator string is taken "as is"
!
    function getBehindLastSeparatorString(s,sps,errmsg) result(res)
    character(len=*) :: s,sps
    type (error_message) :: errmsg
    character(len=max_length_string) :: res
    character(len=28) :: myname = 'getBehindLastSeparatorString'
    character(len=400) :: errstr
    integer :: nlb,icrel,iwb,lsp
!
    call addTrace(errmsg,myname)
    if (len(sps) > max_length_sep_string) then
       write(errstr,*) 'Separator string is longer than ',max_length_sep_string,' characters'
       call add(errmsg,2,errstr,myname)
       call print(errmsg); return
    endif
    if (len_trim(adjustl(s)) == 0) then
       call add(errmsg,2,'Input string either all blank or empty',myname)
       call print(errmsg); return
    endif
!
    lsp = len(sps)
    nlb = len_trim(s)-len_trim(adjustl(s))       ! # of leading blanks
    icrel = index(trim(s(nlb+1:)),sps,.true.)    ! search last occurence of sps
    iwb = icrel+nlb+lsp                          ! begin of word
    if (iwb > len(trim(s))) then
       call add(errmsg,1,'Result string is empty, return empty string',myname)
       res = ''
       return
    endif
    if (len_trim(s)-iwb+1 > max_length_string) then
       call add(errmsg,1,'Result string is longer than max_length_string, return truncated string',myname)
    endif
    res = s(iwb:min(len_trim(s),max_length_string+iwb-1))
    end function getBehindLastSeparatorString
!-----------------------------------------------------------
end module string
