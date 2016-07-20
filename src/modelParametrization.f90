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
!> \brief define (an)elastic material properties
!!
!! \details organizes (an)elastic parametrizations of eart models by defining
!!  the number of parameters n used in a parametrization and the names of the parameters,
!!  whereby each parameter is assigned a specific index between 1 and n
!!
!! \author Florian Schumacher
!! \date Mar 2013
!
module modelParametrization
!
  implicit none
!
  integer, parameter :: character_length_pmtrz = 11 !< can be used outside module to declare character strings for parametrization
  integer, parameter :: character_length_param = 6 !< can be used outside module to declare character strings for parameters
  ! ADD YOUR PARAMETRIZATION HERE
  character(len=65) :: all_valid_pmtrz_param = "'isoLame' ('rho','lambda','mu') ; 'isoVelocity' ('rho','vp','vs')" !< for user/help output only
!
contains
!
!##############################################################################################
!#  O N L Y  IN THE FOLLOWING TWO FUNCTIONS, 
!#   - numberOfParamModelParametrization
!#   - getParamFromIndexModelParametrization
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
    case('isoLame')
       n = 3
    case('isoVelocity')
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
    ! N O T I C E : names of parameters must NOT be empty character strings '' and should not contain any blanks!
    !               the indexing of parameters MUST begin with idx=1 and end with idx=numberOfParamModelParametrization(parametrization)
    !
    select case(parametrization)
    case('isoLame')
       select case(idx)
       case (1); param = 'rho'
       case (2); param = 'lambda'
       case (3); param = 'mu'
       case default; param = ''
       end select
    case('isoVelocity')
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
!! \param parametrization model parametrization
!! \param param model parameter
!! \return logical value which is false if there is no next model parameter
!
  logical function nextParamModelParametrization(parametrization,param,iparam)
    character(len=*), intent(in) :: parametrization
    character(len=character_length_param), optional :: param
    integer, optional :: iparam
    integer :: call_count = 0
    save :: call_count
    call_count = call_count+1
    if(call_count > numberOfParamModelParametrization(parametrization)) then
       call_count = 0
       param = ''
       nextParamModelParametrization = .false.
       return
    endif
    if(present(param)) param = getParamFromIndexModelParametrization(parametrization,call_count)
    if(present(iparam)) iparam = call_count
    nextParamModelParametrization = .true.
  end function nextParamModelParametrization
!
end module modelParametrization
