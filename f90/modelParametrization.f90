!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief define (an)elastic material properties
!!
!! \details organizes (an)elastic parametrizations of eart models by defining
!!  the number of parameters n used in a parametrization and the names of the parameters,
!!  whereby each parameter is assigned a specific index between 1 and n
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module modelParametrization
!
  implicit none
!
  integer, parameter :: character_length_pmtrz = 15 !< can be used outside module to declare character strings for parametrization
  integer, parameter :: character_length_param = 6 !< can be used outside module to declare character strings for parameters
  ! ADD YOUR PARAMETRIZATION HERE
  character(len=87) :: all_valid_pmtrz_param = &
       "'isoLameSI' ('rho','lambda','mu') ; 'isoVelocitySI','isoVelocity1000' ('rho','vp','vs')" !< for user/help output only
!
contains
!
!##############################################################################################
!#  O N L Y  IN THE FOLLOWING TWO FUNCTIONS, 
!#   - numberOfParamModelParametrization
!#   - getParamFromIndexModelParametrization
!#   - getUnitFactorOfParamModelParametrization
!#  YOU NEED TO DEFINE NEW PARAMETRIZATIONS (or modify existing ones). Additionally 
!#  you may want to modify info string all_valid_pmtrz_param
!#  above (just used for user/help output). 
!#  THERE IS NO NEED TO MODIFY THIS MODULE ANYWHERE ESLE
!#
!#  Additionally, of course, you will have to define new routines for ASKI kernel
!#  calculation in module spectralWaveformKernel and assure compatibility with module
!#  kernelReferenceModel and method specific reference model modules.
!##############################################################################################
!
!> \brief returns number of model parameters for given parametrization
!! \param parametrization model parametrization
!! \param n number of model parameters
!! \return number of model parameters for given parametrization. returns n=-1
!!  in order to indicate that incoming model parametrization is not supported by this module
!
  function numberOfParamModelParametrization(parametrization) result(n)
    character(len=*), intent(in) :: parametrization
    integer :: n
    ! ADD YOUR PARAMETRIZATION HERE
    ! 
    ! N O T I C E : the name of a parametrization must NOT be an empty character string '' and should not contain any blanks!
    !               the number of parameters n must be strictly larger than zero: n > 0 !
    !
    select case(parametrization)
    case('isoLameSI')
       n = 3
    case('isoVelocitySI')
       n = 3
    case('isoVelocity1000')
       n = 3
       !case('your_new_parametrization')
       !	n = number_of_your_new_parameters
    case default
       ! RETURN A NEGATIVE NUMBER TO INDICATE AN INVALID PARAMETRIZATION
       n = -1 
    end select
  end function numberOfParamModelParametrization
!------------------------------------------------------------------------
!> \brief for given parameter index of a given parametrization get the name of the parameter
!! \param parametrization model parametrization
!! \param idx index of param in given parametrization
!! \param param name of model parameter of the given index in given parametrization
!! \return name of parameter having given index in given parametrization. returns empty string
!!  if incoming parametrization is invalid, or incoming index is invalid for given parametrization
!
  function getParamFromIndexModelParametrization(parametrization,idx) result(param)
    character(len=*) :: parametrization
    integer :: idx
    character(len=character_length_param):: param
    ! ADD YOUR PARAMETRIZATION HERE
    ! 
    ! N O T I C E : names of parameters must NOT be empty character strings '' and must not contain any blanks!
    !               the indexing of parameters MUST begin with idx=1 and end with idx=numberOfParamModelParametrization(parametrization)
    !
    select case(parametrization)
    case('isoLameSI')
       select case(idx)
       case (1); param = 'rho'
       case (2); param = 'lambda'
       case (3); param = 'mu'
       case default; param = ''
       end select
    case('isoVelocitySI','isoVelocity1000')
       select case(idx)
       case (1); param = 'rho'
       case (2); param = 'vp'
       case (3); param = 'vs'
       case default; param = ''
       end select
    !case('your_new_parametrization')
       !select case(idx)
       !case (1); param = 'name_of_first_new_parameter'
       ! ...
       !case ( numberOfParamModelParametrization('your_new_parametrization') ); param = 'name_of_last_new_parameter'
       !case default; param = ''
       !end select
    case default; param = ''
    end select
  end function getParamFromIndexModelParametrization
!------------------------------------------------------------------------
!> \brief for given parameter get its unit factor which brings its values to SI units
!! \param parametrization model parametrization
!! \param param model parameter
!! \param ios optional status integer indicating whether the incoming parameter (of parametrization) is valid (then ios = 0 on return), i.e. whether the result can be trusted, or not  (then ios < 0 on return, has the same value as factor then.)
!! \param factor unit factor of given parameter; has negative value if param is not a valid parameter of parametrization
!! \return unit factor of given parameter. In case of an error return a negative value:
!!    -1 if incoming parametrization is invalid
!!    -2 if incoming parametrization is valid, but incoming parameter is not valid in the given parametrization
!!    -3 if incoming parameter is an empty string (or contains blanks only)
!
  function getUnitFactorOfParamModelParametrization(parametrization,param,ios) result(factor)
    character(len=*) :: parametrization,param
    integer, optional :: ios
    real :: factor
    ! ADD YOUR PARAMETRIZATION HERE
    ! 
    ! N O T I C E : names of parameters must NOT be empty character strings '' and must not contain any blanks!
    !               the indexing of parameters MUST begin with idx=1 and end with idx=numberOfParamModelParametrization(parametrization)
    !
    select case(parametrization)
    case('isoLameSI')
       ! parametrization 'isoLameSI' has SI units, i.e. all unit factors are 1.0
       select case(param)
       case ('rho'); factor = 1.0 ! unit of model values of 'rho' is kg/m^3
       case ('lambda'); factor = 1.0 ! unit of model values of 'lambda' is Pa = N/m^2
       case ('mu'); factor = 1.0 ! unit of model values of 'mu' is Pa = N/m^2
       case (''); goto 3
       case default; goto 2
       end select
    case('isoVelocitySI')
       ! parametrization 'isoVelocitySI' has SI units, i.e. all unit factors are 1.0
       select case(param)
       case ('rho'); factor = 1.0 ! unit of model values of 'rho' is kg/m^3
       case ('vp'); factor = 1.0 ! unit of model values of 'vp' is m/s
       case ('vs'); factor = 1.0 ! unit of model values of 'vs' is m/s
       case (''); goto 3
       case default; goto 2
       end select
    case('isoVelocity1000')
       ! parametrization 'isoVelocitySI' has units with a unit factor of 1000, 
       ! i.e. the actual model values must be multiplied by a factor of 1000 to get values w.r.t. SI units
       select case(param)
       case ('rho'); factor = 1.0e3 ! unit of model values of 'rho' is g/cm^3
       case ('vp'); factor = 1.0e3 ! unit of model values of 'vp' is km/s
       case ('vs'); factor = 1.0e3  ! unit of model values of 'vs' is km/s
       case (''); goto 3
       case default; goto 2
       end select
    !case('your_new_parametrization')
       !select case(param)
       !case ('name_of_first_new_parameter'); factor = ...
       ! ...
       !case ('name_of_last_new_parameter'); factor = ...
       !case (''); goto 3
       !case default; goto 2
       !end select
    case default; goto 1
    end select
!
    ! if code comes here, eveything went OK, so return
    if(present(ios)) ios = 0
    return
!
    ! if incoming parametrization is invalid, return factor = -1
1   factor = -1.0; goto 4
!
    ! if incoming parametrization is valid, but the given parameter is invalid in this parametrization AND IS NO EMPTY STRING , return factor = -2
2   factor = -2.0; goto 4
!
    ! if incoming parametrization is valid, but the given parameter is an empty string, return factor = -3
    !   we can check this here, because the emtpy string '' is a possible return value of function 
    !   getParamFromIndexModelParametrization (for invalid indices) that a user could have called before giving it to this function
3   factor = -3.0; goto 4
!
    ! in case of any error, finally set ios to factor (if ios is present) and return
4   if(present(ios)) then
       ios = int(factor)
    end if
    return
  end function getUnitFactorOfParamModelParametrization
!
!###############################################################################################
!#  END OF FUNCTIONS, which need to be modified in case of adding new parameters to this module
!###############################################################################################
!
!> \brief for given parameter get its index in the given parametrization
!! \param parametrization model parametrization
!! \param param model parameter
!! \param idx index of param in given parametrization (negative if invalid)
!! \return index of given parameter in given parametrization. In case of an error
!!  return a negative value:
!!    -1 if incoming parametrization is invalid
!!    -2 if incoming parametrization is valid, but incoming parameter is not valid in the given parametrization
!!    -3 if incoming parameter is an empty string (or contains blanks only)
!
  function indexOfParamModelParametrization(parametrization,param) result(idx)
    character(len=*) :: parametrization,param
    integer :: idx
    integer :: i,n
    ! return index=-3, if the incoming parameter is empty (a priori not valid)
    !   need to check this here, because the emtpy string '' is a possible return value
    !   of function getParamFromIndexModelParametrization (for invalid indices)
    if(trim(param) == '') then
       idx = -3 
       return
    endif
    ! return index=-1, if the parametrization does not exist
    n = numberOfParamModelParametrization(parametrization)
    if(n<0) then
       idx = -1
       return
    endif
    ! return index=-2, if the parametrization exists, but param=/'' is not a valid parameter
    idx = -2
    do i = 1,n
       if(getParamFromIndexModelParametrization(parametrization,i)==param) then
          idx = i
          return
       endif
    enddo
  end function indexOfParamModelParametrization
!------------------------------------------------------------------------
!> \brief check if a given model parametrization is supported by this module
!! \param parametrization model parametrization
!! \return logical value if parametrization is supported by this module
!
  function validModelParametrization(parametrization) result(l)
    character(len=*) :: parametrization
    logical :: l
    l = numberOfParamModelParametrization(parametrization) > 0
  end function validModelParametrization
!------------------------------------------------------------------------
!> \brief check if given parameter name is defined for given model parametrization
!! \param parametrization model parametrization
!! \param param model parameter
!! \return logical value if param is actually a parameter of given model parametrization
!
  function validParamModelParametrization(parametrization,param) result(l)
    character(len=*) :: parametrization,param
    logical :: l
    l = .not.(indexOfParamModelParametrization(parametrization,param) < 0)
  end function validParamModelParametrization
!------------------------------------------------------------------------
!> \brief iterate over parameters of a given model parametrization
!! \details the order of iteration MUST be consistent with function indexOfParamModelParametrization
!!  which actually is the case by construction, since in fucntion indexOfParamModelParametrization
!!  there is a loop on all parameters from 1 to numberOfParam
!! \param parametrization model parametrization
!! \param param model parameter
!! \return logical value which is false if there is no next model parameter
!
  logical function nextParamModelParametrization(parametrization,param,iparam,reset)
    character(len=*), intent(in) :: parametrization
    character(len=character_length_param), optional :: param
    integer, optional :: iparam
    logical, optional :: reset
    integer :: call_count = 0
    save :: call_count
!
    ! if this iterator is to be reset, do so
    if(present(reset)) then
       if(reset) goto 1
    end if
!
    ! increase counter
    call_count = call_count+1
!
    ! if the counter rises above the upper bound, reset this iterator
    if(call_count > numberOfParamModelParametrization(parametrization)) goto 1
!
    ! otherwise set the optional variables, if present, and indicate success
    if(present(param)) param = getParamFromIndexModelParametrization(parametrization,call_count)
    if(present(iparam)) iparam = call_count
    nextParamModelParametrization = .true.
! 
    ! IF FUNCTION COMES HERE, RETURN NORMALLY
    return
!
    ! RESET THE ITERATOR
1   call_count = 0
    if(present(param)) param = ''
    if(present(iparam)) iparam = -1
    nextParamModelParametrization = .false.
    return
  end function nextParamModelParametrization
!
end module modelParametrization
