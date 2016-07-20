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
!> \brief holds basic requirements for inversion process
!!
!! \details Basic requirements for an inversion procedure which are
!!  independent of a specific iteration step, such as the main parameter file,
!!  events and stations, station componenent transformation coefficients, etc.
!!  are supervised by this module. 
!!
!! \author Florian Schumacher
!! \date March 2013
!
module inversionBasics
!
  use inputParameter
  use readEventStationFile
  use seismicEventList
  use seismicNetwork
  use componentTransformation
  use errorMessage
  use modelParametrization
!
  implicit none
!
  interface dealloc; module procedure deallocateInversionBasics; end interface
  interface init; module procedure initiateInversionBasics; end interface
  interface mapIfreq2ArrayIndex
     module procedure mapOneIfreq2ArrayIndexInversionBasics
     module procedure mapManyIfreq2ArrayIndexInversionBasics
  end interface mapIfreq2ArrayIndex
  interface operator (.iterpath.); module procedure getCurrentIterPathInversionBasics; end interface
  interface operator (.inpar.); module procedure getInputParameterInversionBasics; end interface
  interface operator (.evlist.); module procedure getEventListInversionBasics; end interface
  interface operator (.statlist.); module procedure getStationListInversionBasics; end interface
  interface operator (.comptrans.); module procedure getComponentTransformationInversionBasics; end interface
  interface operator (.ifreq.); module procedure getMeasuredDataFrequencyIndicesInversionBasics; end interface
!
!> \brief type contains all basic requirements for an inversion procedure
  type inversion_basics
     private
     character(len=350) :: iter_path !< absolute path to directory of iteration step
     type (input_parameter) :: inpar !< content of main parameter file
     type (seismic_event_list) :: event_list !< list of all events involved in this inversion
     type (seismic_network) :: station_list !< list of all stations involved in this inversion
     type (component_transformation) :: cmptr !< component transformation
     integer, dimension(:), pointer :: ifreq => null() !< measured data frequency indices
  end type inversion_basics
!
contains
!
!------------------------------------------------------------------------
!> \brief initiate basic inversion requirements
!! \details Main parfile is read. Everything needed for 
!!  the inversion, independet of each iteration step is calculated and 
!!  saved (if not done yet), or read in (if existent).
!! \param this inversion_basics object
!! \param parfile main parfile of the inversion
!! \param lu file unit used for reading / writing of files
!! \return errmsg error message
!
  subroutine initiateInversionBasics(this,parfile,lu,errmsg)
    ! incoming
    type (inversion_basics) :: this
    character(len=*) :: parfile
    integer :: lu
    ! returning
    type (error_message) :: errmsg
    ! local
    integer :: j,nf,i_partest,ios
    real :: r_partest
    character (len=80), dimension(15) :: main_inpar_keys
    character(len=400) :: errstr
    character(len=23) :: myname = 'initiateInversionBasics'
    ! keywords for input parameter
    data main_inpar_keys/'CURRENT_ITERATION_STEP', 'DEFAULT_VTK_FILE_FORMAT', &
         'FILE_EVENT_LIST', 'FILE_STATION_LIST', 'FORWARD_METHOD', 'ITERATION_STEP_PATH', 'MAIN_PATH_INVERSION', &
         'MEASURED_DATA_FREQUENCY_STEP', 'MEASURED_DATA_NUMBER_OF_FREQ', 'MEASURED_DATA_INDEX_OF_FREQ', &
         'MODEL_PARAMETRIZATION', 'PARFILE_ITERATION_STEP','PATH_MEASURED_DATA', 'PATH_EVENT_FILTER', &
         'PATH_STATION_FILTER'/
!------------------------------------------------------------------------
!  subroutine starts here
!
    call addTrace(errmsg,myname)
!------------------------------------------------------------------------
!  read input parameters
!
    call createKeywordsInputParameter(this%inpar,main_inpar_keys)
    call readSubroutineInputParameter(this%inpar,lu,trim(parfile),errmsg)
    if (.level.errmsg == 2) return
!------------------------------------------------------------------------
!  check consistency of entries in main parfile!!!
!THIS IS ONLY HELPFUL, IF INDEED ALL VALUES ARE CORRECTLY SET IN THE PARFILE. "PRELIMINARY" PARFILES CANNOT BE USED THEN WITH iversionBasics
!ADVANTAGE: then other routines, which rely on inversionBasics (like iteration step basics) can be sure, that all values are valid!
!
    ! check consistency of entries in main parfile
    ! during creation of object this%inpar, there was already checked if the required keywords are present
    ! now check correct type of values here and some special values
!
    i_partest = ival(this%inpar,'CURRENT_ITERATION_STEP',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'CURRENT_ITERATION_STEP' = '"//trim(this%inpar.sval.'CURRENT_ITERATION_STEP')//&
            "' of main parfile is not a valid integer value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    select case(this%inpar.sval.'DEFAULT_VTK_FILE_FORMAT')
    case('ASCII','BINARY')
       ! valid, do nothing
    case default
       write(errstr,*) "parameter 'DEFAULT_VTK_FILE_FORMAT' = '"//trim(this%inpar.sval.'DEFAULT_VTK_FILE_FORMAT')//&
            "' of main parfile is not valid; must be either 'ASCII' or 'BINARY'"
       call add(errmsg,2,trim(errstr),myname)
       return       
    end select
!
    r_partest = rval(this%inpar,'MEASURED_DATA_FREQUENCY_STEP',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'MEASURED_DATA_FREQUENCY_STEP' = '"//trim(this%inpar.sval.&
            'MEASURED_DATA_FREQUENCY_STEP')//"' of main parfile is not a valid real value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(r_partest .le. 0.) then
       write(errstr,*) "parameter 'MEASURED_DATA_FREQUENCY_STEP' = ",r_partest,&
            " of main parfile is not valid; must be positive"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    nf = ival(this%inpar,'MEASURED_DATA_NUMBER_OF_FREQ',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'MEASURED_DATA_NUMBER_OF_FREQ' = '"//trim(this%inpar.sval.&
            'MEASURED_DATA_NUMBER_OF_FREQ')//"' of main parfile is not a valid integer value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(nf .le. 0) then
       write(errstr,*) "parameter 'MEASURED_DATA_NUMBER_OF_FREQ' = ",nf,&
            " of main parfile is not valid; must be positive"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    this%ifreq => ivecp(this%inpar,'MEASURED_DATA_INDEX_OF_FREQ',nf,iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'MEASURED_DATA_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
            'MEASURED_DATA_INDEX_OF_FREQ')//"' is not a vector of ",nf," integers"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(any(this%ifreq < 0)) then
       write(errstr,*) "entries of vector 'MEASURED_DATA_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
            'MEASURED_DATA_INDEX_OF_FREQ')//"' must not be negative"
       call add(errmsg,2,trim(errstr),myname)
       deallocate(this%ifreq)
       return
    end if
    do j = 1,nf
       if(any(this%ifreq(j+1:nf) == this%ifreq(j))) then
          write(errstr,*) "vector 'MEASURED_DATA_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
               'MEASURED_DATA_INDEX_OF_FREQ')//"' must not contain multiple entries"
          call add(errmsg,2,trim(errstr),myname)
          deallocate(this%ifreq)
          return
       end if
    end do ! j
!
    if(.not.validModelParametrization(this%inpar.sval.'MODEL_PARAMETRIZATION')) then
       write(errstr,*) "parameter 'MODEL_PARAMETRIZATION' = '"//trim(this%inpar.sval.'MODEL_PARAMETRIZATION')//&
            "' of main parfile is not valid; valid parametrizations and parameters are: "//trim(all_valid_pmtrz_param)
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    write(this%iter_path,"(2a,i3.3,a)") trim(this%inpar.sval.'MAIN_PATH_INVERSION'),trim(this%inpar.sval.'ITERATION_STEP_PATH'), &
         this%inpar.ival.'CURRENT_ITERATION_STEP','/'
!------------------------------------------------------------------------
!  read event file, create event list and write events.vtk
!
    call createEventListFromEventFile(this%inpar.sval.'FILE_EVENT_LIST',lu,&
         'ASKI_events',this%event_list,errmsg)
    if (.level.errmsg == 2) return
!------------------------------------------------------------------------
!  read station file, create station list and write stations.vtk
!
    call createStationListFromStationFile(this%inpar.sval.'FILE_STATION_LIST',lu,&
         'ASKI_stations',this%station_list,errmsg)
    if (.level.errmsg == 2) return
!
    if(.csys.this%event_list /= .csys.this%station_list) then
       call add(errmsg,2,"the coordinate system '"//.csys.this%event_list//"' used in event list differs from '"&
            //.csys.this%station_list//"' used in station list",myname)
       return
    end if
!------------------------------------------------------------------------
!  create component transformation object and check value of COORDINATE_SYSTEM
!
    if(.csys.this%station_list=='S') then
       call createComponentTransformation(this%cmptr,.stations.this%station_list)
    elseif(.csys.this%station_list=='C') then
       call createComponentTransformation(this%cmptr)
    else
       call add(errmsg,2,"coordinate system '"//.csys.this%station_list//&
            "' of station list is neither 'S' nor 'C': cannot create component transformation",myname)
       return
    endif
  end subroutine initiateInversionBasics
!------------------------------------------------------------------------
!> \brief for one given frequency index ifreq_in return its index position in array this%ifreq
!! \param this iteration step basics
!! \param ifreq_in integer frequency index for which the array position in this%ifreq should be returned
!! \param idx integer such that this%ifreq(idx)==ifreq_in, or -1 if no entry of this%ifreq has value ifreq_in
!! \return position index of ifreq_in in array this%ifreq, or -1 if no entry of this%ifreq has value ifreq_in
!
  function mapOneIfreq2ArrayIndexInversionBasics(this,ifreq_in) result(idx)
    type (inversion_basics), intent(in) :: this
    integer :: ifreq_in,idx
    ! local
    integer :: j
    idx = -1
    if(.not.associated(this%ifreq)) return
    do j=1,size(this%ifreq)
       if(this%ifreq(j) == ifreq_in) then
          idx = j
          exit
       end if
    end do ! j
  end function mapOneIfreq2ArrayIndexInversionBasics
!------------------------------------------------------------------------
!> \brief for a vector of given frequency indices ifreq_in return their index positions in array this%ifreq
!! \param this iteration step basics
!! \param ifreq_in vector of integer frequency indices for which the array positions in this%ifreq should be returned
!! \param idx of same length as ifreq_in and this%ifreq(idx(i))==ifreq_in(i), or not associated if any ifreq_in(i) is not in this%ifreq
!! \return position indices of values in ifreq_in in array this%ifreq, or not associated if any ifreq_in(i) is not in this%ifreq
!
  function mapManyIfreq2ArrayIndexInversionBasics(this,ifreq_in) result(idx)
    type (inversion_basics), intent(in) :: this
    integer, dimension(:) :: ifreq_in
    integer, dimension(:), pointer :: idx
    ! local
    integer :: size_ifreq_in,i,j
!
    idx => null()
    if(.not.associated(this%ifreq)) return
    size_ifreq_in = size(ifreq_in)
    if(size_ifreq_in <= 0) return
!
    allocate(idx(size_ifreq_in))
    idx = -1
    do i=1,size_ifreq_in
       do j=1,size(this%ifreq)
          if(this%ifreq(j) == ifreq_in(i)) then
             idx(i) = j
             exit
          end if
       end do ! j
       if(idx(i)==-1) then
          deallocate(idx); nullify(idx)
          return
       end if
    end do ! i
  end function mapManyIfreq2ArrayIndexInversionBasics
!------------------------------------------------------------------------
!> \brief deallocate inversion basics object
!! \param this inversion_basics object
!
  subroutine deallocateInversionBasics(this)
    type (inversion_basics) :: this
    this%iter_path = ''
    call dealloc(this%inpar)
    call dealloc(this%event_list)
    call dealloc(this%station_list)
    call dealloc(this%cmptr)
    if(associated(this%ifreq)) deallocate(this%ifreq)
  end subroutine deallocateInversionBasics
!------------------------------------------------------------------------
!> \brief get absolute path to directory of current iteration step
!! \param this inversion_basics object
!! \param iter_path iteration step path
!! \return absolute path to directory of current iteration step
!
  function getCurrentIterPathInversionBasics(this) result(iter_path)
    type (inversion_basics), intent(in) :: this
    character(len=350) :: iter_path
    iter_path = this%iter_path
  end function getCurrentIterPathInversionBasics
!------------------------------------------------------------------------
!> \brief get input_parameter object contained in inversion_basics object
!! \param this inversion_basics object
!! \return inpath input_parameter object contained in this 
!
  function getInputParameterInversionBasics(this) result(inpar)
    type (inversion_basics), intent(in) :: this
    type (input_parameter) :: inpar
    inpar = this%inpar
  end function getInputParameterInversionBasics
!------------------------------------------------------------------------
!> \brief get event_list object contained in inversion_basics object
!! \param this inversion_basics object
!! \return event_list seismic_event_list object contained in this 
!
  function getEventListInversionBasics(this) result(event_list)
    type (inversion_basics), intent(in) :: this
    type (seismic_event_list) :: event_list
    event_list = this%event_list
  end function getEventListInversionBasics
!------------------------------------------------------------------------
!> \brief get station_list object contained in inversion_basics object
!! \param this inversion_basics object
!! \return station_list seismic_network object contained in this 
!
  function getStationListInversionBasics(this) result(station_list)
    type (inversion_basics), intent(in) :: this
    type (seismic_network) :: station_list
    station_list = this%station_list
  end function getStationListInversionBasics
!------------------------------------------------------------------------
!> \brief get component transformation contained in inversion_basics object
!! \param this inversion_basics object
!! \param cmptr component transformation
!! \return iteration step path contained in this
!
  function getComponentTransformationInversionBasics(this) result(cmptr)
    type (inversion_basics), intent(in) :: this
    type (component_transformation) :: cmptr
    cmptr = this%cmptr
  end function getComponentTransformationInversionBasics
!------------------------------------------------------------------------
!> \brief get measured data frequency indices contained in inversion_basics object
!! \param this inversion_basics object
!! \param ifreq frequency indices
!! \return measured data frequency indices contained in this
!
  function getMeasuredDataFrequencyIndicesInversionBasics(this) result(ifreq)
    type (inversion_basics), intent(in) :: this
    integer, dimension(:), pointer :: ifreq
    ifreq => this%ifreq
  end function getMeasuredDataFrequencyIndicesInversionBasics
!
end module inversionBasics
