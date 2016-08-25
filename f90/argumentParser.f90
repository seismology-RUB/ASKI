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
!-----------------------------------------------------------------
!> \brief Parse command line arguments
!
module argumentParser
    use string
    use realloc
    implicit none
    interface init; module procedure initArgumentParser; end interface
    interface addPosarg; module procedure addPositionalArgumentParser; end interface
    interface addOption; module procedure addOptionArgumentParser; end interface
    interface parse; module procedure parseArgumentParser; end interface
    interface usage; module procedure usageArgumentParser; end interface
    interface document; module procedure documentArgumentParser; end interface
    interface dealloc; module procedure deallocArgumentParser; end interface
    interface operator (.ival.); module procedure getIvalArgumentParser; end interface
    interface operator (.rval.); module procedure getRvalArgumentParser; end interface
    interface operator (.dval.); module procedure getDvalArgumentParser; end interface
    interface operator (.sval.); module procedure getSvalArgumentParser; end interface
    interface operator (.ivec.); module procedure getIvecArgumentParser; end interface
    interface operator (.rvec.); module procedure getRvecArgumentParser; end interface
    interface operator (.dvec.); module procedure getDvecArgumentParser; end interface
    interface operator (.svec.); module procedure getSvecArgumentParser; end interface
    interface operator (.optset.); module procedure isOptionSetArgumentParser; end interface
    interface operator (.errmsg.); module procedure getErrmsgArgumentParser; end interface
    integer, parameter :: max_length_option = 10
    integer, parameter :: max_length_posnam = 20
    integer, parameter :: max_length_argtyp = 5
    integer, parameter :: num_argument_parser = 5                                    ! expected number of posargs and opts
    type argument_parser
       private
       character(len=max_length_string) :: progname                                  ! program name
       character(len=max_length_string) :: description                               ! program description
       integer :: npos                                                               ! number of valid pos args
       integer :: nopt                                                               ! number of valid options
!
       character(len=max_length_posnam), dimension(:), pointer :: posnam => null()   ! name of positional arguments
       character(len=max_length_string), dimension(:), pointer :: posmsg => null()   ! posarg help message
       character(len=max_length_argtyp), dimension(:), pointer :: postyp => null()   ! positional argument types
!
       character(len=max_length_option), dimension(:), pointer :: option => null()   ! option strings
       character(len=max_length_string), dimension(:), pointer :: optmsg => null()   ! option help message
       character(len=max_length_argtyp), dimension(:), pointer :: opttyp => null()   ! option argument types
       character(len=max_length_string), dimension(:), pointer :: optdef => null()   ! option default value (char!!)
       logical, dimension(:), pointer :: hasarg => null()                            ! option has argument or not
!
       character(len=max_length_string), dimension(:), pointer :: posarg => null()   ! values of positional arguments
       character(len=max_length_string), dimension(:), pointer :: optarg => null()   ! values of optional arguments
       logical, dimension(:), pointer :: optset => null()                            ! option is set or not
!
       character(len=max_length_argtyp), dimension(8) :: valid_argtyp                ! valid argument types
       type (error_message) :: errmsg                                                ! collects errors
    end type
!
    type value_argument_parser
       integer :: ival
       real :: rval
       double precision :: dval
       character(len=max_length_string) :: sval
       integer, dimension(:), pointer :: ivec => null()
       real, dimension(:), pointer :: rvec => null()
       double precision, dimension(:), pointer :: dvec => null()
       character(len=max_length_string), dimension(:), pointer :: svec => null() 
    end type
!
contains
!----------------------------------------------------------------
!> \brief Initialization of argumentParser
!
    subroutine initArgumentParser(this,progname,description)
    type (argument_parser) :: this
    character(len=*) :: progname,description
    this%progname = progname
    this%description = description
    this%valid_argtyp = (/'ival','sval','rval','dval','ivec','svec','rvec','dvec'/)
    allocate(this%posnam(num_argument_parser))
    allocate(this%posmsg(num_argument_parser))
    allocate(this%postyp(num_argument_parser))
    allocate(this%option(num_argument_parser))
    allocate(this%optmsg(num_argument_parser))
    allocate(this%optdef(num_argument_parser))
    allocate(this%opttyp(num_argument_parser))
    allocate(this%hasarg(num_argument_parser))
    allocate(this%posarg(num_argument_parser))
    allocate(this%optarg(num_argument_parser))
    allocate(this%optset(num_argument_parser))
    this%nopt = 0
    this%npos = 0
    call new(this%errmsg,'initArgumentParser')
    end subroutine initArgumentParser
!---------------------------------------------------------------
!> \brief Deallocate argumentParser object
!
    subroutine deallocArgumentParser(this)
    type (argument_parser) :: this
    if (associated(this%posnam)) deallocate(this%posnam)
    if (associated(this%posmsg)) deallocate(this%posmsg)
    if (associated(this%postyp)) deallocate(this%postyp)
    if (associated(this%option)) deallocate(this%option)
    if (associated(this%optmsg)) deallocate(this%optmsg)
    if (associated(this%optdef)) deallocate(this%optdef)
    if (associated(this%opttyp)) deallocate(this%opttyp)
    if (associated(this%hasarg)) deallocate(this%hasarg)
    if (associated(this%posarg)) deallocate(this%posarg)
    if (associated(this%optarg)) deallocate(this%optarg)
    if (associated(this%optset)) deallocate(this%optset)
    call dealloc(this%errmsg)
    end subroutine deallocArgumentParser
!---------------------------------------------------------------
!> \brief Add positional argument
!
    subroutine addPositionalArgumentParser(this,posnam,postyp,posmsg)
    type (argument_parser) :: this
    character(len=*) :: posnam,posmsg,postyp
    character(len=27) :: myname = 'addPositionalArgumentParser'
!
!  check length of posnam, posmsg and validity of postyp
!
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    if (len_trim(posnam) > max_length_posnam) then
       call add(this%errmsg,2,'name of positional argument too long',myname)
       call print(this%errmsg); return
    endif
    if (len_trim(posmsg) > max_length_string) then
       call add(this%errmsg,2,'message for positional argument too long',myname)
       call print(this%errmsg); return
    endif
    if (.not. typeValidArgumentParser(this,postyp)) then
       call add(this%errmsg,2,'type of positional argument "'//trim(posnam)+'" is invalid',myname)
       call print(this%errmsg); return
    endif
    if (mod(this%npos,num_argument_parser) == 0) then
       this%posnam => reallocate(this%posnam,this%npos+num_argument_parser)
       this%posmsg => reallocate(this%posmsg,this%npos+num_argument_parser)
       this%postyp => reallocate(this%postyp,this%npos+num_argument_parser)
    endif
    this%npos = this%npos+1
    this%posnam(this%npos) = posnam
    this%posmsg(this%npos) = posmsg
    this%postyp(this%npos) = postyp
    end subroutine addPositionalArgumentParser
!---------------------------------------------------------------
!> \brief Add optional argument
!
    subroutine addOptionArgumentParser(this,option,hasarg,optmsg,opttyp,optdef)
    type (argument_parser) :: this
    character(len=*) :: option,optmsg
    character(len=*), optional :: opttyp,optdef
    logical :: hasarg
    character(len=25) :: myname = 'addOptionalArgumentParser'
!
    if (hasarg) then
       if (.not.present(opttyp) .or. .not.present(optdef)) then
          call add(this%errmsg,2,'Option "'//trim(option)//'" has an argument. You must specify type and default value',myname)
          call print(this%errmsg); return
       endif
    endif
!
!  check length of option, optmsg and validity of opttyp
!
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    if (len_trim(option) > max_length_option) then
       call add(this%errmsg,2,'name of option too long',myname)
       call print(this%errmsg); return
    endif
    if (len_trim(optmsg) > max_length_string) then
       call add(this%errmsg,2,'message for option too long',myname)
       call print(this%errmsg); return
    endif
    if (len_trim(optdef) > max_length_string) then
       call add(this%errmsg,2,'default value for option too long',myname)
       call print(this%errmsg); return
    endif
    if (hasarg .and. .not. typeValidArgumentParser(this,opttyp)) then
       call add(this%errmsg,2,'type of argument for option "'//trim(option)+'" is invalid',myname)
       call print(this%errmsg); return
    endif
!
    if (mod(this%nopt,num_argument_parser) == 0) then
       this%option => reallocate(this%option,this%nopt+num_argument_parser)
       this%optmsg => reallocate(this%optmsg,this%nopt+num_argument_parser)
       this%optdef => reallocate(this%optdef,this%nopt+num_argument_parser)
       this%opttyp => reallocate(this%opttyp,this%nopt+num_argument_parser)
       this%hasarg => reallocate(this%hasarg,this%nopt+num_argument_parser)
    endif
    this%nopt = this%nopt+1
    this%option(this%nopt) = option
    this%optmsg(this%nopt) = optmsg
    this%hasarg(this%nopt) = hasarg
    if (present(optdef)) then
       this%optdef(this%nopt) = optdef
    else
       this%optdef(this%nopt) = 'none'
    endif
    if (present(opttyp)) then
       this%opttyp(this%nopt) = opttyp
    else
       this%opttyp(this%nopt) = 'none'
    endif
    end subroutine addOptionArgumentParser
!---------------------------------------------------------------
!> \brief Parse arguments
!
    subroutine parseArgumentParser(this)
    type (argument_parser) :: this
    character(len=19) :: myname = 'parseArgumentParser'
    character(len=max_length_string) :: arg,errstr
    integer :: nargs,jarg,jpos,jopt
    character(len=max_length_option) :: option
!
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    allocate(this%optset(this%nopt),this%optarg(this%nopt),this%posarg(this%npos))
    this%optset = .false.
    this%optarg = 'none'
    this%posarg = 'none'
!
    nargs = iargc()
    if (nargs == 0 .and. this%npos > 0) goto 3
    jarg = 1
    jpos = 0
    do while (jarg <= nargs)
       call getarg(jarg,arg)
       if (arg(1:1) == '-') then          ! argument is an option string
          if (len_trim(arg) > max_length_option) then
             call add(this%errmsg,2,'option string is too long, truncated',myname)
             return
          endif
          option = arg(1:len_trim(arg))
          jopt = 1
          do while(jopt <= this%nopt)             ! search if among defined options
             if (option.equal.(this%option(jopt))) then
                this%optset(jopt) = .true.
                if (this%hasarg(jopt)) then
                   jarg = jarg+1
                   if (jarg > nargs) then
                      call add(this%errmsg,2,'no argument for option: '+option,myname)
                      return
                   endif
                   call getarg(jarg,this%optarg(jopt))  ! get optional argument if hasarg
                endif
                exit
             endif
             jopt = jopt+1
          enddo
          if (jopt > this%nopt) then
             call add(this%errmsg,2,'invalid option "'//trim(option)//'" found on command line',myname)
             return
          endif
       else                                ! positional argument
          jpos = jpos+1
          if(jpos>this%npos) then
             write(errstr,*) 'expect only ',this%npos,' positional argument(s), but at least one more is given'
             call add(this%errmsg,2,errstr,myname)
             return
          end if
          this%posarg(jpos) = arg
       endif
       jarg = jarg+1
    enddo
!
!  checks
!
    if (jpos < this%npos) goto 2
    return
!
2   call add(this%errmsg,2,'Not enough positional arguments',myname)
    return
3   call add(this%errmsg,2,'No arguments at all',myname)
    end subroutine parseArgumentParser
!---------------------------------------------------------------------------
!> \brief type of argument valid?
!
    function typeValidArgumentParser(this,argtyp) result(res)
    type (argument_parser) :: this
    character(len=*) :: argtyp
    logical :: res
    integer :: j
!
    res = .false.
    if (len_trim(argtyp) > max_length_argtyp) return
    j = 1
    do while(j <= size(this%valid_argtyp))
       if (argtyp .equal. (this%valid_argtyp(j))) then
          res = .true.
          exit
       endif
       j = j+1
    enddo
    end function typeValidArgumentParser
!---------------------------------------------------------------------------
!> \brief print usage message
!
    subroutine usageArgumentParser(this)
    type (argument_parser) :: this
    integer :: j
!
    print '(80(1h-))'
    print '(30x,a9)','| USAGE |'
    print '(30x,9(1h-))'
    print '(a9,4x,a)','Program: ',trim(this%progname)
    print '(a13,a)','Description: ',trim(this%description)
    print *,''
    print '(a)','Positional arguments:'
    do j = 1,this%npos
       print '(a24,a9,a5,a16,a)',trim(this%posnam(j)),': Type = ',trim(this%postyp(j)),', Description = ',trim(this%posmsg(j))
    enddo
    print *,''
    print '(a)','Optional arguments:'
    do j = 1,this%nopt
       print '(a14,a9,a5,a16,a,a13,a,a1)',trim(this%option(j)),': Type = ',trim(this%opttyp(j)),&
            ', Description = ',trim(this%optmsg(j)),', (Default = ',trim(this%optdef(j)),')'
    enddo
    print '(80(1h-))'
    end subroutine usageArgumentParser
!-----------------------------------------------------------
!> \brief Get error messages
!
    function getErrmsgArgumentParser(this) result(res)
    type (argument_parser), intent(in) :: this
    type (error_message) :: res
    res = this%errmsg
    end function getErrmsgArgumentParser
!------------------------------------------------------------
!> \brief get value of positional or optional argument
!
    function getValueArgumentParser(this,argnam,typ) result(res)
    type (argument_parser) :: this
    character(len=*) :: argnam
    character(len=*) :: typ
    type (value_argument_parser) :: res
    character(len=22) :: myname = 'getValueArgumentParser'
    character(len=max_length_string) :: errstr
    integer :: jopt,jpos,j,nword
    logical :: optset,ispos,isopt
!
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    if (argnam(1:1) == '-') goto 2       ! optional argument
!
!  check if argnam is a positional argument and if type fits
!
    ispos = isPositionalArgumentParser(this,argnam,jpos)
    if (.not. ispos) then
       write(errstr,*) 'Positional argument "',trim(argnam),'" does not exist' 
       call add(this%errmsg,2,errstr,myname)
       call print(this%errmsg)
    endif
    if (.not. (this%postyp(jpos).equal.typ)) then
       write(errstr,*) 'Positional argument "',trim(argnam),'" is not of type "',trim(typ),'"' 
       call add(this%errmsg,2,errstr,myname)
       call print(this%errmsg)
       return
    endif
    errstr = this%posarg(jpos)
    goto 3
!
!  check if argnam is an optional argument and if type fits
!
 2  isopt = isOptionalArgumentParser(this,argnam,jopt,optset)
    if (.not. isopt) then
       write(errstr,*) 'Optional argument "',trim(argnam),'" does not exist' 
       call add(this%errmsg,2,errstr,myname)
       call print(this%errmsg)
       return
    endif
    if (.not.(this%opttyp(jopt).equal.typ)) then
       write(errstr,*) 'Optional argument "',trim(argnam),'" is not of type "',trim(typ),'"' 
       call add(this%errmsg,2,errstr,myname)
       call print(this%errmsg)
       return
    endif
    if (.not.optset) then
       errstr = this%optdef(jopt)
    else
       errstr = this%optarg(jopt)
    endif
!
!  read out values
!
 3  if (index(typ,'val') > 0) then
       select case (trim(typ))
       case ('ival'); read(errstr,*) res%ival
       case ('rval'); read(errstr,*) res%rval
       case ('sval'); read(errstr,'(a)') res%sval
       case ('dval'); read(errstr,*) res%dval
       end select
    else if(index(typ,'vec') > 0) then
       nword = countWordsString(errstr,' ',this%errmsg)
       if (.level.this%errmsg == 2) return
       select case (trim(typ))
       case ('ivec')
          allocate(res%ivec(nword))
          read(errstr,*) (res%ivec(j),j = 1,nword)
       case ('rvec')
          allocate(res%rvec(nword))
          read(errstr,*) (res%rvec(j),j = 1,nword)
       case ('dvec')
          allocate(res%dvec(nword))
          read(errstr,*) (res%dvec(j),j = 1,nword)
       case ('svec')
          res%svec => getWordsString(errstr,' ',this%errmsg)
          if (.level.this%errmsg == 2) then; nullify(res%svec); return; endif
       end select
    endif
    end function getValueArgumentParser
!----------------------------------------------------------------------------------
!> \brief Check if provided argument name is positional argument and, if so, return index
!
    function isPositionalArgumentParser(this,argnam,jpos) result(res)
    type (argument_parser) :: this
    character(len=*) :: argnam
    integer :: jpos
    logical :: res
    integer :: j
!
    res = .false.
    jpos = 0
    do j = 1,this%npos
       if (argnam.equal.this%posnam(j)) then
          res = .true.
          jpos = j
          exit
       endif
    enddo
    end function isPositionalArgumentParser
!----------------------------------------------------------------------------------
!> \brief Check if provided argument name is optional argument 
!!  and, if so, return index and flag whether option is set
!
    function isOptionalArgumentParser(this,argnam,jopt,optset) result(res)
    type (argument_parser) :: this
    character(len=*) :: argnam
    integer :: jopt
    logical optset
    logical :: res
    integer :: j
!
    res = .false.
    optset = .false.
    jopt = 0
    do j = 1,this%nopt
       if (argnam.equal.this%option(j)) then
          res = .true.
          jopt = j
          if (this%optset(j)) optset = .true.
          exit
       endif
    enddo
    end function isOptionalArgumentParser
!---------------------------------------------------------------------------
!> \brief print actual commandline contents
!
    subroutine documentArgumentParser(this)
    type (argument_parser) :: this
    integer :: j
!
    print '(80(1h-))'
    print '(25x,41a)','| DOCUMENTATION OF COMMANDLINE CONTENTS |'
    print '(25x,41(1h-))'
    print '(a9,4x,a)','Program: ',trim(this%progname)
    print '(a13,a)','Description: ',trim(this%description)
    print *,''
    print '(a)','Positional arguments:'
    do j = 1,this%npos
       print '(a24,a9,a5,a10,a)',trim(this%posnam(j)),': Type = ',trim(this%postyp(j)),', Value = ',trim(this%posarg(j))
    enddo
    print *,''
    print '(a)','Optional arguments:'
    do j = 1,this%nopt
       if (.not.this%optset(j)) then
          print '(a14,a9,a5,a12,a)',trim(this%option(j)),': Type = ',trim(this%opttyp(j)),', Default = ',trim(this%optdef(j))
       else
          print '(a14,a9,a5,a12,a)',trim(this%option(j)),': Type = ',trim(this%opttyp(j)),', Value   = ',trim(this%optarg(j))
       endif
    enddo
    print '(80(1h-))'
    end subroutine documentArgumentParser
!-----------------------------------------------------------
!> \brief Get integer argument value
!
    function getIvalArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    integer :: res
    character(len=21) :: myname = 'getIvalArgumentParser'
    type (value_argument_parser) :: val
!
    res = -999
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'ival')
    res = val%ival
    end function getIvalArgumentParser
!-----------------------------------------------------------
!> \brief Get real argument value
!
    function getRvalArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    real :: res
    character(len=21) :: myname = 'getRvalArgumentParser'
    type (value_argument_parser) :: val
!
    res = -999.9
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'rval')
    res = val%rval
    end function getRvalArgumentParser
!-----------------------------------------------------------
!> \brief Get double argument value
!
    function getDvalArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    double precision :: res
    character(len=21) :: myname = 'getDvalArgumentParser'
    type (value_argument_parser) :: val
!
    res = -999.9d0
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'dval')
    res = val%dval
    end function getDvalArgumentParser
!-----------------------------------------------------------
!> \brief Get character argument value
!
    function getSvalArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    character(len=max_length_string) :: res
    character(len=21) :: myname = 'getSvalArgumentParser'
    type (value_argument_parser) :: val
!
    res = 'none'
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'sval')
    res = val%sval
    end function getSvalArgumentParser
!-----------------------------------------------------------
!> \brief Get integer vector
!
    function getIvecArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    integer, dimension(:), pointer :: res
    character(len=21) :: myname = 'getIvecArgumentParser'
    type (value_argument_parser) :: val
!
    nullify(res)
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'ivec')
    res => val%ivec
    end function getIvecArgumentParser
!-----------------------------------------------------------
!> \brief Get real vector
!
    function getRvecArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    real, dimension(:), pointer :: res
    character(len=21) :: myname = 'getRvecArgumentParser'
    type (value_argument_parser) :: val
!
    nullify(res)
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'rvec')
    res => val%rvec
    end function getRvecArgumentParser
!-----------------------------------------------------------
!> \brief Get double vector
!
    function getDvecArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    double precision, dimension(:), pointer :: res
    character(len=21) :: myname = 'getDvecArgumentParser'
    type (value_argument_parser) :: val
!
    nullify(res)
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'dvec')
    res => val%dvec
    end function getDvecArgumentParser
!-----------------------------------------------------------
!> \brief Get character vector
!
    function getSvecArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    character(len=max_length_string), dimension(:), pointer :: res
    character(len=21) :: myname = 'getSvecArgumentParser'
    type (value_argument_parser) :: val
!
    nullify(res)
    if (.level.this%errmsg == 2) return
    call addTrace(this%errmsg,myname)
    val = getValueArgumentParser(this,argnam,'svec')
    res => val%svec
    end function getSvecArgumentParser
!--------------------------------------------------------
!> \brief Is option set
!
    function isOptionSetArgumentParser(this,argnam) result(res)
    type (argument_parser), intent(in) :: this
    character(len=*), intent(in) :: argnam
    logical :: res,isopt,optset
    integer :: jopt
    character(len=max_length_string) :: errstr
    character(len=25) :: myname = 'isOptionSetArgumentParser'
!
    res = .false.
    call addTrace(this%errmsg,myname)
    isopt = isOptionalArgumentParser(this,argnam,jopt,optset)
    if (.not. isopt) then
       write(errstr,*) 'Optional argument "',trim(argnam),'" does not exist' 
       call add(this%errmsg,2,errstr,myname)
       call print(this%errmsg)
       return
    endif
    res = optset
    end function isOptionSetArgumentParser
!
end module argumentParser
