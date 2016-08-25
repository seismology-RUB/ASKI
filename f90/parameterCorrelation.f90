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
!> \brief defines correlations between certain model parameters of modelParametrization
!!
!! \details So far, there is only a quick hack done for constructing objects (for reasons of testing):
!!  In the parameter correlation file as given by the main parameter file, assume only ONE line, 
!!  containing for parameter vs the parameters of  rho , vp (in that order, space separated).
!!  Please finish implementing this module by introduce a correlation matrix ("reciprocally-symmetric" 
!   with diagonal elements equal to 1) which is constructed from a file.
!!  The correlation coefficients should be accessed by functions using names of parameters to comunicate
!!  
!!
!! \author Florian Schumacher
!! \date Nov 2015
!! 
!! \copyright &copy; GNU Public License
!
module parameterCorrelation
!
  use modelParametrization
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateParameterCorrelation; end interface
  interface correlateAnyParameters
     module procedure atAllAnyParameterCorrelation
     module procedure byIndexAnyParameterCorrelation
     module procedure byNameAnyParameterCorrelation     
  end interface correlateAnyParameters
  interface correlateParameters
     module procedure byIndexParameterCorrelation
     module procedure byNameParameterCorrelation     
  end interface correlateParameters
!interface ; module procedure ; end interface
!interface operator (..); module procedure ; end interface
!
  type parameter_correlation
     private
     character(len=50) :: mode = '' !< either 'CORRELATE_KERNELS', (or for future implementation 'PASSIVE_UPDATE', or whatever)
!
     character(len=character_length_pmtrz) :: model_parametrization
     logical, dimension(:,:), pointer :: correlate => null() !< indicates for each model parameter p_i, whether it correlates to other parameters q_j, i.e. whether p_i = c1*q_1 + ... + cn*q_n
     real, dimension(:,:), pointer :: corr_coef => null() !< correlation coefficients
     !! FS FS ADD PARAMETER CORRELATION HERE
     ! introduce a correlation matrix ("reciprocally-symmetric with diagonal elements equal to 1")
     ! which is constructed from a file (so far only a quick hack below)
     ! the correlation coefficients should be accessed by functions using names of parameters to comunicate
  end type parameter_correlation
!
contains
!
!------------------------------------------------------------------------
!> \brief constructor
!! \param 
!
  subroutine createParameterCorrelation(this,mode,filename,lu,parametrization,errmsg)
    type (parameter_correlation) :: this
    character(len=*) :: mode,filename,parametrization
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=26) :: myname = 'createParameterCorrelation'
    character(len=400) :: errstr,line
    integer :: nparam,iparam
    integer :: ios
!
    call addTrace(errmsg,myname)
    ! if object was already allocated, deallocate now before creating a new one
    call deallocateParameterCorrelation(this)
!
    select case(mode)
    case('CORRELATE_KERNELS','NONE')! 'PASSIVE_UPDATE', ... 
       this%mode = mode
    case default
       call add(errmsg,2,"incoming correlation mode '"//trim(mode)//"' not supported",myname)
       return
    end select
!
    ! in case of no parameter correlation, do not read the file and do not store coefficients
    if(this%mode == 'NONE') return
!
    open(unit=lu,file=trim(filename),form='FORMATTED',status='OLD',action='READ',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "could not open file, iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
!
    if(.not.validModelParametrization(parametrization)) then
       call add(errmsg,2,"incoming model parametrization '"//trim(parametrization)//"' is not valid",myname)
       goto 2
    end if
!
    ! allocate for this model parametrization
    nparam = numberOfParamModelParametrization(parametrization)
    allocate(this%correlate(nparam,nparam),this%corr_coef(nparam,nparam))
    this%correlate = .false.
    this%corr_coef = 0.
    this%model_parametrization = parametrization
!
    ! ! read parameter which is to be correlated
    ! read(lu,"(a400)",iostat=ios) line
    ! if(ios/=0) goto 2
    ! if(validParamModelParametrization(parametrization,trim(adjustl(trim(line))))) then
    !    iparam = indexOfParamModelParametrization(parametrization,param)
    !    ! read line containing the parameter names of those parameters to which parameter iparam correlates
    !    read(lu,"(a400)",iostat=ios) line
    !    if(ios/=0) goto 2
       
    !    ! FS FS GET NUMBER OF PARAMETER AND PARAMETERS FROM STRING (by string module?)
       
    !    ! read line containing the correlation coefficients of the  parameters (same number as line above)
    !    read(lu,"(a400)",iostat=ios) line
    !    if(ios/=0) goto 2
       
    ! else
    !    call add(errmsg,"model parameter '"//trim(line)//"' is not valid",myname)
    !    goto 2
    ! end if

    ! do while(ios/=0)
    !    ! now read the rest of the blocks ...
    ! end do ! while ios/=0
    
    !! FS FS QUICK HACK: ASSUME ONLY ONE LINE, CONTAINING FOR vs THE PARAMETERS OF   rho,  vp  (in that order, space separated)
     read(lu,"(a400)",iostat=ios) line
     if(ios/=0) goto 2
     this%correlate(3,1:2) = .true. ! leave the diagonal if this%correlate .false. is actually sensible!!
     read(line,*,iostat=ios) this%corr_coef(3,1),this%corr_coef(3,2)
     ! this%corr_coef(3,3) = 0.0 ! leaving the own coefficient 0.0 is actually sensible, since the coefficients
     ! are meant in the context of ADDING to the own unweighted kernel. 
!
1   return
2   call add(errmsg,2,"object could not be created",myname)
    call deallocateParameterCorrelation(this)
  end subroutine createParameterCorrelation
!------------------------------------------------------------------------
  subroutine deallocateParameterCorrelation(this)
    type (parameter_correlation) :: this
    if (associated(this%correlate)) deallocate(this%correlate)
    if (associated(this%corr_coef)) deallocate(this%corr_coef)
    this%mode = ''
    this%model_parametrization = ''
  end subroutine deallocateParameterCorrelation
!------------------------------------------------------------------------
!> \brief check whether there are any correlations at all
!! \return logical that indicates whether there are any correlations at all
!
  function atAllAnyParameterCorrelation(this) result(l)
    type (parameter_correlation) :: this
    logical :: l
    l = .false.
    select case(this%mode)
    case('NONE','')
       return
    end select
    if(.not.associated(this%correlate)) return
    l = any(this%correlate)
  end function atAllAnyParameterCorrelation
!------------------------------------------------------------------------
!> \brief check whether a specific parameter should be correlated to others
!! \details define parameter by their index (consistent with module modelParametrization)
!! \return logical that indicates whether iparam'th parameter should be correlated to other parameters
!
  function byIndexAnyParameterCorrelation(this,iparam) result(l)
    type (parameter_correlation) :: this
    integer :: iparam
    logical :: l
    l = .false.
    select case(this%mode)
    case('NONE','')
       return
    end select
    if(.not.associated(this%correlate)) return
    if(iparam < 1 .or. iparam > size(this%correlate)) return
    l = any(this%correlate(iparam,:))
  end function byIndexAnyParameterCorrelation
!------------------------------------------------------------------------
!> \brief check whether a specific parameter should be correlated to others
!! \details define parameter by their name (consistent with module modelParametrization)
!! \return logical that indicates whether iparam'th parameter should be correlated to other parameters
!
  function byNameAnyParameterCorrelation(this,name_param) result(l)
    type (parameter_correlation) :: this
    character(len=*) :: name_param
    logical :: l
    l = byIndexAnyParameterCorrelation(this,indexOfParamModelParametrization(this%model_parametrization,name_param))
  end function byNameAnyParameterCorrelation
!------------------------------------------------------------------------
!> \brief check whether a specific parameter should be correlated to others
!! \details define parameter by their index (consistent with module modelParametrization)
!! \return logical that indicates whether iparam'th parameter should be correlated to other parameters
!
  function byIndexParameterCorrelation(this,iparam_main,iparam_corr,c_corr) result(l)
    type (parameter_correlation) :: this
    integer :: iparam_main,iparam_corr
    real :: c_corr
    logical :: l
    l = .false.; c_corr = 0.
    select case(this%mode)
    case('NONE','')
       return
    end select
    if(.not.associated(this%correlate)) return
    if(iparam_main < 1 .or. iparam_main > size(this%correlate)) return
    if(iparam_corr < 1 .or. iparam_corr > size(this%correlate)) return
    l = this%correlate(iparam_main,iparam_corr)
    if(l) c_corr = this%corr_coef(iparam_main,iparam_corr)
  end function byIndexParameterCorrelation
!------------------------------------------------------------------------
!> \brief check whether a specific parameter should be correlated to others
!! \details define parameter by their name (consistent with module modelParametrization)
!! \return logical that indicates whether iparam'th parameter should be correlated to other parameters
!
  function byNameParameterCorrelation(this,name_param_main,name_param_corr,c_corr) result(l)
    type (parameter_correlation) :: this
    character(len=*) :: name_param_main,name_param_corr
    real :: c_corr
    logical :: l
    integer :: iparam_main,iparam_corr
    iparam_main = indexOfParamModelParametrization(this%model_parametrization,name_param_main)
    iparam_corr = indexOfParamModelParametrization(this%model_parametrization,name_param_corr)
    l = byIndexParameterCorrelation(this,iparam_main,iparam_corr,c_corr)
  end function byNameParameterCorrelation
!------------------------------------------------------------------------
!> \brief for a given model parameter, iterate over all parameters which are correlated to it
!! \param  model parametrization
!! \param param this parameter correlation
!! \param 
!! \return logical value which is false if there is no next model parameter
!
  logical function nextCorrelationParameter(this,name_param_main,name_param_corr,c_corr)
    type (parameter_correlation) :: this
    character(len=*) :: name_param_main,name_param_corr
    real, optional :: c_corr
    integer :: iparam_corr,iparam
    integer :: count = 0
    integer :: iparam_main = 0
    save :: count,iparam_main
!
    ! COUNT IS NOT A CALL COUNT, but rather is the index iparam_coror of the latest found parameter
    ! which is to be correlated to param_main (or 0 on first call), hence this index can also jump from
    ! one call to another:  
    !   e.g. if for a given main parameter the parameters 1 , 2 and 5 correlate,
    !   then count will have values 1, 2 and 5 on exits of this routine
!
    if(iparam_main==0) then
       ! if this is the first call, memorize the main parameter (index)
       iparam_main = indexOfParamModelParametrization(this%model_parametrization,name_param_main)
    else
       ! otherwise, if the main parameter differs from the one of previous calls, return false and set everything back to start
       if(iparam_main /= indexOfParamModelParametrization(this%model_parametrization,name_param_main)) goto 100
    end if
!
    ! now find the next parameter which is to be correlated to the param_main:
    ! on entry, count is the index iparam_corr found at the last call (or 0 at first call)
    ! so start loop on parameter indices at count+1
    iparam_corr = -1
    do iparam = count+1,numberOfParamModelParametrization(this%model_parametrization)
       if(this%correlate(iparam_main,iparam)) then
          iparam_corr = iparam
          exit
       end if
    end do
    ! if there was no correlation index found, set everything back to start
    if(iparam_corr == -1) goto 100
!
    ! otherwise store the relevant values
    name_param_corr = getParamFromIndexModelParametrization(this%model_parametrization,iparam_corr)
    if(present(c_corr)) c_corr = this%corr_coef(iparam_main,iparam_corr)
    ! memorize the index of the parameter to be correlated. 
    count = iparam_corr
!
    ! if routine comes here, everything went alright, so return memorizing
    nextCorrelationParameter = .true.
    return
!
    ! if there is a goto 100, the outome is negative, so set everything back to start
100 count = 0
    iparam_main = 0
    name_param_corr = ''
    if(present(c_corr)) c_corr = 0.0
    nextCorrelationParameter = .false.
  end function nextCorrelationParameter
!
end module parameterCorrelation
