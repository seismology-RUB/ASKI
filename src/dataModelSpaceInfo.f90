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
!> \brief Define and hold information about a set of data and a set of model parameters
!!
!! \details This module defines a type containing all information about the data samples 
!!  and model parameters, i.e. data space and model space, which should be used for a 
!!  certain operation like setting up the kernel matrix, doing kernel focussing, etc.
!!  A data sample is defined by a station, a station component, an event, a frequency index
!!  and distinction between imaginary or real part (of the complex valued data).  
!!  A model parameter is defined by a parametrization (e.g. 'isoVelocity', 
!!  allowing for parameters like 'rho', 'vp', 'vs'), the respective parameter (e.g. 'vs') 
!!  and an inversion grid cell).
!!  Only indices refering to properties defined in a different object (iteration step info, 
!!  inversion grid) are hold in this module, like station/event index, inversion grid cell 
!!  index etc. 
!!  For parallelized application, these data samples and model parameters may be defined and 
!!  used locally and do not need to represent the whole data and model space.
!!
!! \author Florian Schumacher
!! \date April 2013
!------------------------------------------------------------------------
module dataModelSpaceInfo
!
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use componentTransformation
  use modelParametrization
  use integrationWeights
  use errorMessage
!
  implicit none
!
  ! by disallowing calls to certain routines, it is tried to reduce the efford of 
  ! checking validity of parameters and assume correct values in those routines (must
  ! be checked outside the routines somewhere else in the module)
  !public
  private :: addDataSamplesDataModelSpaceInfo,addModelParametersDataModelSpaceInfo
!
  interface dealloc; module procedure deallocateDataModelSpaceInfo; end interface
  interface getArrayEvid; module procedure getArrayEvidDataSamplesDataModelSpaceInfo; end interface
  interface getArrayStaname; module procedure getArrayStanameDataSamplesDataModelSpaceInfo; end interface
  interface getIndxDataSamples; module procedure getIndicesDataSamplesDataModelSpaceInfo; end interface
  interface getIndxModelParam; module procedure getIndicesModelParametersDataModelSpaceInfo; end interface
  interface allComp; module procedure getAllDifferentCompDataSamplesDataModelSpaceInfo; end interface
! FS FS
! instead of having ONE routine which returnes a vector of indices of data samples (model parameters)
! which have certain properties, we can have ONE routine (well two actually: one for data samples, one for model parameters)
! which gets ALL optional parameters staname,evid,comp,ifreq,imre (param,cell respectively)
! and returnes a vector of indices of data samples (model parameters) which have all the requested properties
! (or a pointer to null, if there is no data sample (model parameter) having the requested properties)
! (could be realized like: get arrays of all properties, then use identity array (indices) and "where" statement
! and finally pack statement (using count))
! FS FS
  interface operator (.ndata.); module procedure getNdataDataModelSpaceInfo; end interface
  interface operator (.nparam.); module procedure getNparamDataModelSpaceInfo; end interface
  interface operator (.evid.); module procedure getEvidDataSampleDataModelSpaceInfo; end interface
  interface operator (.staname.); module procedure getStanameDataSampleDataModelSpaceInfo; end interface
  interface operator (.comp.); module procedure getCompDataSampleDataModelSpaceInfo; end interface
  interface operator (.ifreq.); module procedure getIfreqDataSampleDataModelSpaceInfo; end interface
  interface operator (.imre.); module procedure getImreDataSampleDataModelSpaceInfo; end interface
  interface operator (.pmtrz.); module procedure getParametrizationDataModelSpaceInfo; end interface
  interface operator (.param.); module procedure getParamModelParameterDataModelSpaceInfo; end interface
  interface operator (.cell.); module procedure getCellModelParameterDataModelSpaceInfo; end interface
!
!> \brief meta info defining rows (data space) and columns (model space) of inversion matrix
  type data_model_space_info
     private
     ! data space
     integer :: ndata = 0 !< total number of data in data space (size of data_space_info)

     character(len=character_length_evid), dimension(:), pointer :: evid => null() !< event ID
     character(len=character_length_staname), dimension(:), pointer :: staname => null() !< name of station
     character(len=character_length_component), dimension(:), pointer :: comp => null() !< data component (valid components defined by module componentTransformation)
     integer, dimension(:), pointer :: ifreq => null() !< frequency index of all frequencies (either nf1 <= ifreq <= nf2 for equally spaced, or 1 <= ifreq <= nfreq for arbitrary freq.)
     character(len=2), dimension(:), pointer :: imre => null() !< imaginary part (imre='im') or real part (imre='re')

     ! model space
     integer :: nparam = 0 !< total number of parameters in model space (size of model_space_info)
     character(len=character_length_pmtrz) :: parametrization = '' !< defines the parametrization of the model space (no multiple parametrizations allowed)

     character(len=character_length_param), dimension(:), pointer :: param => null() !< e.g. 'lambda','vs',... (dependent on parametrization defined in data_model_space_info object)
     integer, dimension(:), pointer :: cell => null() !< index of inversion grid cell which corresponds to this parameter
  end type data_model_space_info
!
  integer, parameter :: max_character_length_evid_staname = max(character_length_evid,character_length_staname)
!
contains
!------------------------------------------------------------------------
!> \brief add data samples to existing ones of this data model space info object
!! \details This routine is private, hence cannot be called from outside the module. 
!!  Hence, this routine ASSUMES all incoming values to be valid, there is no check done here!
!!  Furthermore it is assumed that all incoming arrays have the same size, except in case of 
!!  all_combinations=.true. when all incoming values are combined to nev*nstat*ncomp*nfreq*nimre
!!  data samples. So all routines of this module which call this routine have to check validity of
!!  the values and need to assure the correct sizes of the arrays in the respective cases.
!! \param this data_model_space_info object to which data samples are added
!! \param evids array of event IDs
!! \param stanames array of station names
!! \param comp array of components
!! \param ifreq array of frequency indices
!! \param imre array containing values 'im' or 're' indicating imaginary or real part of the data sample
!! \param all_combinations optional logical to indicate whether data samples for all nev*nstat*ncomp*nfreq*nimre
!!  combinations of the incoming arrays are added or not. If not present (or .false.) all incoming arrays are assumed to have same size!
!
  subroutine addDataSamplesDataModelSpaceInfo(this,evids,stanames,comp,ifreq,imre,all_combinations)
    ! incoming
    type (data_model_space_info) :: this
    character(len=*), dimension(:) :: evids,stanames
    character(len=*), dimension(:) :: comp
    integer, dimension(:) :: ifreq
    character(len=*), dimension(:) :: imre
    logical, optional :: all_combinations
    ! local
    integer :: nstat,nev,ncomp,nfreq,nimre
    integer :: jstat,jev,jcomp,jfreq,jimre
    integer :: ndata_new,jdata
    logical :: add_all_combinations
    character(len=character_length_evid), dimension(:), pointer :: evid_tmp
    character(len=character_length_staname), dimension(:), pointer :: staname_tmp
    character(len=character_length_component), dimension(:), pointer :: comp_tmp
    integer, dimension(:), pointer :: ifreq_tmp
    character(len=2), dimension(:), pointer :: imre_tmp
    !
    if(present(all_combinations)) then
       add_all_combinations = all_combinations
    else
       add_all_combinations = .false.
    end if
!
    nev = size(evids); nstat = size(stanames); ncomp = size(comp); nfreq = size(ifreq); nimre = size(imre)
!
    if(add_all_combinations) then
       ndata_new = nev * nstat * ncomp * nfreq * nimre
    else
       ndata_new = nev ! it is assumed that all incoming arrays have the same size!!
    end if
!
    if(ndata_new .le. 0) return
!
    ! reallocate
    evid_tmp => this%evid
    staname_tmp => this%staname
    comp_tmp => this%comp
    ifreq_tmp => this%ifreq
    imre_tmp => this%imre
    allocate(this%evid(this%ndata+ndata_new))
    allocate(this%staname(this%ndata+ndata_new))
    allocate(this%comp(this%ndata+ndata_new))
    allocate(this%ifreq(this%ndata+ndata_new))
    allocate(this%imre(this%ndata+ndata_new))
    if(this%ndata > 0) then
       this%evid(1:this%ndata) = evid_tmp
       this%staname(1:this%ndata) = staname_tmp
       this%comp(1:this%ndata) = comp_tmp
       this%ifreq(1:this%ndata) = ifreq_tmp
       this%imre(1:this%ndata) = imre_tmp
    end if
    if(associated(evid_tmp)) deallocate(evid_tmp)
    if(associated(staname_tmp)) deallocate(staname_tmp)
    if(associated(comp_tmp)) deallocate(comp_tmp)
    if(associated(ifreq_tmp)) deallocate(ifreq_tmp)
    if(associated(imre_tmp)) deallocate(imre_tmp)
!
    ! add samples
    if(add_all_combinations) then
!
       jdata = 0
       do jev = 1,nev
          do jstat = 1,nstat
             do jcomp = 1,ncomp
                do jfreq = 1,nfreq
                   do jimre = 1,nimre
                      jdata = jdata + 1
                      this%evid(this%ndata+jdata) = evids(jev)
                      this%staname(this%ndata+jdata) = stanames(jstat)
                      this%comp(this%ndata+jdata) = comp(jcomp)
                      this%ifreq(this%ndata+jdata) = ifreq(jfreq)
                      this%imre(this%ndata+jdata) = imre(jimre)
                   end do ! jimre
                end do ! jfreq
             end do ! jcomp
          end do ! jstat
       end do ! jev
!
    else ! add_all_combinations
!
       do jdata = 1,ndata_new
          this%evid(this%ndata+jdata) = evids(jdata)
          this%staname(this%ndata+jdata) = stanames(jdata)
          this%comp(this%ndata+jdata) = comp(jdata)
          this%ifreq(this%ndata+jdata) = ifreq(jdata)
          this%imre(this%ndata+jdata) = imre(jdata)
       end do ! jdata
!
    endif ! add_all_combinations
!
    this%ndata = this%ndata + ndata_new
  end subroutine addDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief add model parameters to data model space info by combination of all input items
!! \details This routine is private, hence cannot be called from outside the module. 
!!  Hence, this routine ASSUMES all incoming values to be valid, there is no check done here!
!!  Furthermore it is assumed that all incoming arrays have the same size, except in case of 
!!  all_combinations=.true. when all incoming values are combined to nparam*ncell
!!  model parameters. So all routines of this module which call this routine have to check validity of
!!  the values and need to assure the correct sizes of the arrays in the respective cases.
!! \param this data_model_space_info object to which model parameters are added
!! \param param array containing names of parameter (defines this%param)
!! \param cell array of indices of inversion grid cells
!! \param all_combinations optional logical to indicate whether model parameters for all nparam*ncell
!!  combinations of the incoming arrays are added or not. If not present (or .false.) all incoming arrays are assumed to have same size!
!
  subroutine addModelParametersDataModelSpaceInfo(this,param,cell,all_combinations)
    ! incoming
    type (data_model_space_info) :: this
    character(len=*), dimension(:) :: param
    integer, dimension(:) :: cell
    logical, optional :: all_combinations
    ! local
    integer :: nparam,ncell
    integer :: jparam,jcell
    integer :: nparam_new,iparam
    logical :: add_all_combinations
    character(len=character_length_param), dimension(:), pointer :: param_tmp
    integer, dimension(:), pointer :: cell_tmp
    !
    if(present(all_combinations)) then
       add_all_combinations = all_combinations
    else
       add_all_combinations = .false.
    end if
    !
    nparam = size(param); ncell = size(cell)
!
    if(add_all_combinations) then
       nparam_new = nparam * ncell
    else
       nparam_new = nparam ! it is assumed that all incoming arrays have the same size!!       
    end if
!
    if(nparam_new .le. 0) return
!
    ! reallocate
    param_tmp => this%param
    cell_tmp => this%cell
    allocate(this%param(this%nparam+nparam_new))
    allocate(this%cell(this%nparam+nparam_new))
    if(this%nparam > 0) then
       this%param(1:this%nparam) = param_tmp
       this%cell(1:this%nparam) = cell_tmp
    end if
    if(associated(param_tmp)) deallocate(param_tmp)
    if(associated(cell_tmp)) deallocate(cell_tmp)
!
    ! add parameters
    if(add_all_combinations) then
!
       iparam = 0
       do jparam = 1,nparam
          do jcell = 1,ncell
             iparam = iparam + 1
             this%param(this%nparam+iparam) = param(jparam)
             this%cell(this%nparam+iparam) = cell(jcell)
          end do ! jcell
       end do ! jparam
!
    else ! add_all_combinations
!
       do iparam = 1,nparam_new
          this%param(this%nparam+iparam) = param(iparam)
          this%cell(this%nparam+iparam) = cell(iparam)
       end do ! iparam
!
    endif ! add_all_combinations
!
    this%nparam = this%nparam + nparam_new
  end subroutine addModelParametersDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief from data model space input file, a data_model_space object is filled
!! \details here, simply the two subroutines createDataSamplesFromFileDataModelSpaceInfo
!!  and createModelParametersFromFileDataModelSpaceInfo are called
!! \param this data_model_space_info object to be created
!! \param eventlist list of all seismic events, to check validity of event IDs
!! \param stationlist list of all seismic stations, to check validity of station names
!! \param ifreq_valid array of valid frequency indices
!! \param parametrization all parameters are assumed to be of this parametrization 
!! \param ntot_invgrid total number of inversion grid cells to check valid range of cell indices in file
!! \param intw integration weights object which containes wavefield point localization inside invgrid
!! \param filename file name to use
!! \param lu file unit to use
!! \param errmsg error message
!
  subroutine createFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,&
       parametrization,ntot_invgrid,intw,filename,lu,errmsg)
    ! incoming
    type (data_model_space_info) :: this
    type (seismic_event_list) :: eventlist
    type (seismic_network) :: stationlist
    integer, dimension(:) :: ifreq_valid
    character(len=*) :: parametrization
    integer :: ntot_invgrid
    type (integration_weights) :: intw
    character(len=*) :: filename
    integer :: lu
    ! outgoing
    type (error_message) :: errmsg
    ! local
    character(len=32) :: myname = 'createFromFileDataModelSpaceInfo'
!
    call addTrace(errmsg,myname)
!
    call createDataSamplesFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,filename,lu,errmsg)
    if(.level.errmsg == 2) return
!
    call createModelParametersFromFileDataModelSpaceInfo(this,parametrization,ntot_invgrid,intw,filename,lu,errmsg)
    if(.level.errmsg == 2) return
  end subroutine createFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief create data samples from the DATA SAMPLES block of the data model space input file
!! \details This routine searches for the line "DATA SAMPLES" in the given file and reads in the
!!  block following this line. 
!!  Format of 'DATA SAMPLES' block:
!!     line of form: 'PATHS value', where value is either 'ALL' (all paths for a given set of events and stations
!!        are used), or 'SPECIFIC' (a specific definition of paths as a series of event/station pairs follows below)
!!     If 'PATHS ALL', the next two lines are of form 'nev evid_1 ... evid_n' and 'nstat staname_1 ... staname_n', defining the set
!!        of event IDs and station names which are combined to all possible paths.
!!     line of form: 'COMPONENTS value', where value is either 'ALL' (for all paths, the same components are used) or 'SPECIFIC'
!!        (only allowed if 'PATHS SPECIFIC', for each path a specific set of components may be defined)
!!     If 'COMPONENTS ALL', the next line is of form 'ncomp comp_1 ... comp_n' defining the components for all paths.
!!     line of form: 'FREQUENCIES value', where value is either 'ALL' (for all paths, the same frequency indices are used) or 'SPECIFIC'
!!        (only allowed if 'PATHS SPECIFIC', for each path a specific set of frequency indices may be defined)
!!     If 'FREQUENCIES ALL', the next line is of form 'nfreq ifreq_1 ... ifreq_n' defining the frequency indices for all paths.
!!     line of form: 'IMRE value', where value is either 'ALL' (for all paths, the same set of imaginary/real parts are used) or 'SPECIFIC'
!!        (only allowed if 'PATHS SPECIFIC', for each path a specific set of imaginary/real parts may be defined)
!!     If 'IMRE ALL', the next line is of form 'nimre imre_1 ... imre_n' defining imaginary (i.e. imre_i = 'im') or real parts 
!!        (imre_i = 're') for all paths.
!!    If 'PATHS SPECIFIC', the following line must contain the number of paths npaths which should be taken, followed by 
!!        npahts blocks of lines, each defining the path and the data samples for that path. These blocks constist of at least one
!!        line containing the eventid/stationname pair 'evid staname'. For each keyword 'COMPONENTS', 'FREQUENCIES' and 'IMRE' -if 'SPECIFIC'- 
!!        one line is added to such a block of lines, in the same form as the line following value 'ALL' (see above, e.g. 'nfreq ifreq_1 ... ifreq_n'),
!!        defining the specific components, frequencies or set of imaginary/real parts for that path. 
!! \param this data model space info
!! \param eventlist list of all seismic events, to check validity of event IDs
!! \param stationlist list of all seismic stations, to check validity of station names
!! \param ifreq_valid array of valid frequency indices
!! \param filename file name to use
!! \param lu file unit, with which the input file is already open
!! \param errmsg error message
!
  subroutine createDataSamplesFromFileDataModelSpaceInfo(this,eventlist,stationlist,ifreq_valid,filename,lu,errmsg)
    type (data_model_space_info) :: this
    type (seismic_event_list) :: eventlist
    type (seismic_network) :: stationlist
    integer, dimension(:) :: ifreq_valid
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    type (error_message) :: errmsg2
    character(len=43) :: myname = 'createDataSamplesFromFileDataModelSpaceInfo'
    character(len=500) :: line
    integer :: ios,iline,nvalid,ndata_before
    character(len=40) :: key,val
    integer :: nev,jev,nstat,jstat,npath,jpath,ncomp,jcomp,nfreq,jfreq,nimre,jimre
    character(len=character_length_evid) :: evid
    character(len=character_length_evid), dimension(:), allocatable :: evids,evids_tmp
    character(len=character_length_staname) :: staname
    character(len=character_length_staname), dimension(:), allocatable :: stanames,stanames_tmp
    character(len=character_length_component), dimension(:), allocatable :: comp,comp_tmp
    integer, dimension(:), allocatable :: ifreq,ifreq_tmp
    character(len=2), dimension(:), allocatable :: imre,imre_tmp
    logical, dimension(:), allocatable :: valid
    logical :: specific_paths,specific_comp,specific_freq,specific_imre,line_data_samples_found
    !
    call addTrace(errmsg,myname)
    iline = 0
    ndata_before = this%ndata
    !
    open(unit=lu,file=trim(filename),status='old',form='formatted',action='read',iostat=ios)
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not open file '"//trim(filename)//"', iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    !
    ! TODO in the future:
    ! if introducing some header to dataModelSpaceInfo files, read in header here, and react accordingly
    !
    ! parse through whole file until line "DATA SAMPLES" and start reading in block from there
    line_data_samples_found = .false.
    do while(ios==0)
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) exit
       if(line == "DATA SAMPLES") then
          line_data_samples_found = .true.
          exit
       endif
    enddo
    if(.not.line_data_samples_found) then
       write(errstr,*) "could not find line 'DATA SAMPLES' in file '"//trim(filename)//"', searching until line ",iline-1
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    write(errstr,*) "start reading 'DATA SAMPLES' block from line ",iline," of file '"//trim(filename)//"'"
    call add(errmsg,0,trim(errstr),myname)
    !
    ! PATHS
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword PATHS followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'PATHS') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in DATA SAMPLES block not supported. 'PATHS' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
!
    select case(trim(val))
    case('ALL')
       specific_paths = .false.
       !
       ! event IDs
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of events nev as first value on that line followed by nev eventIDs"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) nev
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of events as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(nev .le. 0) then
          close(lu)
          write(errstr,*) "number of events ",nev," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(evids_tmp(nev))
       read(line,*,iostat=ios) nev,evids_tmp
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read ",nev," event IDs starting from second value of line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid event IDs
       nvalid = nev
       allocate(valid(nvalid))
       do jev = 1,nev
          errmsg2 = searchEventidSeismicEventList(eventlist,evids_tmp(jev))
          valid(jev) = .level.errmsg2 == 0
          if(.not. valid(jev)) then
             nvalid = nvalid - 1
             write(errstr,*) jev,"'th event id '",trim(evids_tmp(jev)),"' on line ",iline," is not contained in event list, ",&
                  "hence it is excluded!"
             call add(errmsg,1,trim(errstr),myname)
          endif
          call dealloc(errmsg2)
       enddo ! jev
       ! if there were invalid event IDs, exclude them. if there are no valid event IDs, return
       if(nvalid<nev) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid event ID on line ",iline,", hence no data samples are created!"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(evids(nvalid))
          evids = pack(evids_tmp,valid)
       else
          allocate(evids(nev))
          evids = evids_tmp
       end if
       deallocate(evids_tmp,valid)
       !
       ! station names
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of stations nstat as first value on that line, followed by nstat station names"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) nstat
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of stations as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(nstat .le. 0) then
          close(lu)
          write(errstr,*) "number of stations ",nstat," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(stanames_tmp(nstat))
       read(line,*,iostat=ios) nstat,stanames_tmp
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read ",nstat," station names starting from second value of line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid station names
       nvalid = nstat
       allocate(valid(nvalid))
       do jstat = 1,nstat
          errmsg2 = searchStationNameSeismicNetwork(stationlist,stanames_tmp(jstat))
          valid(jstat) = .level.errmsg2 == 0
          if(.not. valid(jstat)) then
             nvalid = nvalid - 1
             write(errstr,*) jstat,"'th station name '",trim(stanames_tmp(jstat)),"' on line ",iline,&
                  " is not contained in station list, hence it is excluded!"
             call add(errmsg,1,trim(errstr),myname)
          endif
          call dealloc(errmsg2)
       enddo ! jstat
       ! if there were invalid station names, exclude them. if there are no valid station names, return
       if(nvalid<nstat) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid station name on line ",iline,", hence no data samples are created!"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(stanames(nvalid))
          stanames = pack(stanames_tmp,valid)
       else
          allocate(stanames(nstat))
          stanames = stanames_tmp
       end if
       deallocate(stanames_tmp,valid)
    case('SPECIFIC') ! val
       specific_paths = .true.
    case default
       write(errstr,*) "'PATHS' specification '"//trim(val)//"' on line ",iline," in DATA SAMPLES block not supported."//&
            " 'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val	
    !
    ! COMPONENTS
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword COMPONENTS followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'COMPONENTS') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in DATA SAMPLES block not supported. 'COMPONENTS' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('ALL')
       specific_comp = .false.
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of components ncomp as first value on that line, followed by ncomp components"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) ncomp
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of components as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(ncomp .le. 0) then
          close(lu)
          write(errstr,*) "number of components ",ncomp," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(comp_tmp(ncomp))
       read(line,*,iostat=ios) ncomp,comp_tmp
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read ",ncomp," components starting from second value of line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid components
       nvalid = ncomp
       allocate(valid(nvalid))
       do jcomp = 1,ncomp
          valid(jcomp) = validComponent(comp_tmp(jcomp))
          if(.not.valid(jcomp)) then
             nvalid = nvalid - 1
             write(errstr,*) jcomp,"'th component '",trim(comp_tmp(jcomp)),"' on line ",iline, &
                  " is not valid, hence it is excluded! Valid components are: "//trim(all_valid_components)
             call add(errmsg,1,trim(errstr),myname)
          endif
       enddo ! jcomp
       ! if there were invalid components, exclude them. if there are no valid components, return
       if(nvalid<ncomp) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid components on line ",iline,", hence no data samples are created!"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(comp(nvalid))
          comp = pack(comp_tmp,valid)
       else
          allocate(comp(ncomp))
          comp = comp_tmp
       end if
       deallocate(comp_tmp,valid)
    case('SPECIFIC') ! val
       if(.not.specific_paths) then
          call add(errmsg,2,"'COMPONENTS SPECIFIC' is only allowed in case of 'PATHS SPECIFIC'",myname)
          close(lu)
          return
       endif
       specific_comp = .true.
    case default
       write(errstr,*) "'COMPONENTS' specification '"//trim(val)//"' on line ",iline,"in  DATA SAMPLES block not supported. "//&
            "'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val
    !
    ! FREQUENCIES
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword FREQUENCIES followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'FREQUENCIES') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in DATA SAMPLES block not supported. 'FREQUENCIES' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('ALL')
       specific_freq = .false.
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of frequencies nfreq as first value on that line, followed by nfreq frequency indices"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) nfreq
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of frequencies as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(nfreq .le. 0) then
          close(lu)
          write(errstr,*) "number of frequencies ",nfreq," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(ifreq_tmp(nfreq))
       read(line,*,iostat=ios) nfreq,ifreq_tmp
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid frequencies
       nvalid = nfreq
       allocate(valid(nvalid))
       do jfreq = 1,nfreq
          valid(jfreq) = any(ifreq_valid == ifreq_tmp(jfreq))
          if(.not.valid(jfreq)) then
             nvalid = nvalid - 1
             write(errstr,*) jfreq,"'th frequency index '",ifreq_tmp(jfreq),"' on line ",iline, &
                  " is invalid (valid indices: ",ifreq_valid,"), hence it is excluded!"
             call add(errmsg,1,trim(errstr),myname)
          endif
       enddo ! jfreq
       ! if there were invalid frequency indices, exclude them. if there are no valid frequencies, return
       if(nvalid<nfreq) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid frequency index on line ",iline,", hence no data samples are created!"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(ifreq(nvalid))
          ifreq = pack(ifreq_tmp,valid)
       else
          allocate(ifreq(nfreq))
          ifreq = ifreq_tmp
       end if
       deallocate(ifreq_tmp,valid)
    case('SPECIFIC') ! val
       if(.not.specific_paths) then
          call add(errmsg,2,"'FREQUENCIES SPECIFIC' is only allowed in case of 'PATHS SPECIFIC'",myname)
          close(lu)
          return			
       endif
       specific_freq = .true.
    case default
       write(errstr,*) "'FREQUENCIES' specification '"//trim(val)//"' on line ",iline,"in  DATA SAMPLES block not supported. "//&
            "'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val
    !
    ! IMRE
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword IMRE followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'IMRE') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in DATA SAMPLES block not supported. 'IMRE' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('ALL')
       specific_imre = .false.
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of imaginary/real part values nimre as first value on that line, followed by nimre of such values"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) nimre
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of imaginary/real part values as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(nimre .le. 0) then
          close(lu)
          write(errstr,*) "number of imaginary/real part values ",nimre," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(imre_tmp(nimre))
       read(line,*,iostat=ios) nimre,imre_tmp
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid imre values
       nvalid = nimre
       allocate(valid(nvalid))
       do jimre = 1,nimre
          select case(imre_tmp(jimre))
          case('im','re')
             valid(jimre) = .true.
          case default
             nvalid = nvalid - 1
             valid(jimre) = .false.
             write(errstr,*) jimre,"'th imaginary/real part value '",trim(imre_tmp(jimre)),&
                  "' on line ",iline," is not supported, hence it is excluded! Must be either 'im' or 're'."
             call add(errmsg,1,trim(errstr),myname)
          end select
       enddo ! jimre
       ! if there were invalid frequency indices, exclude them. if there are no valid frequencies, return
       if(nvalid<nimre) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid imaginary/real part value on line ",iline,", hence no data samples are created!"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(imre(nvalid))
          imre = pack(imre_tmp,valid)
       else
          allocate(imre(nimre))
          imre = imre_tmp
       end if
       deallocate(imre_tmp,valid)
    case('SPECIFIC') ! val
       if(.not.specific_paths) then
          call add(errmsg,2,"'IMRE SPECIFIC' is only allowed in case of 'PATHS SPECIFIC'",myname)
          close(lu)
          return			
       endif
       specific_imre = .true.
    case default
       write(errstr,*) "'IMRE' specification '"//trim(val)//"' on line ",iline,"in  DATA SAMPLES block not supported. "//&
            "'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val
!
! now iterate over lines defining specific paths (if there are any)
!
    if(specific_paths) then
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of paths as first (and only) value on that line"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) npath
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of paths as first (and only) value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(npath .le. 0) then
          close(lu)
          write(errstr,*) "number of paths ",npath," (first (and only) value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       do jpath = 1,npath
          read(lu,"(a)",iostat=ios) line; iline = iline+1
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                  "; expected ",jpath,"'th eventID/stationname pair (out of ",npath,") on that line"
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          read(line,*,iostat=ios) evid,staname
          if(ios/=0) then
             write(errstr,*) "could not read 'eventID stationName' pair on line",iline,&
                  ", hence excluding this path!, iostat = ",ios
             call add(errmsg,1,trim(errstr),myname)
             ! skip all remaining lines which belong to this path and cycle
             if(specific_comp) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_freq) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_imre) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             cycle
          end if ! ios/=0
          ! check if the event ID is valid
          errmsg2 = searchEventidSeismicEventList(eventlist,evid)
          if(.level.errmsg2 /= 0) then
             write(errstr,*) "event id '",trim(evid),"' (first entry) on line ",iline," is not contained in event list",&
                  ", hence excluding this path!"
             call add(errmsg,1,trim(errstr),myname)
             call dealloc(errmsg2)
             ! skip all remaining lines which belong to this path and cycle
             if(specific_comp) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_freq) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_imre) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             cycle
          end if
          call dealloc(errmsg2)
          ! check if the station name is valid
          errmsg2 = searchStationNameSeismicNetwork(stationlist,staname)
          if(.level.errmsg2 /= 0) then
             write(errstr,*) "station name '",trim(staname),"' (second entry) on line ",iline," is not contained in ",&
                  "station list, hence excluding this path!"
             call add(errmsg,1,trim(errstr),myname)
             call dealloc(errmsg2)
             ! skip all remaining lines which belong to this path and cycle
             if(specific_comp) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_freq) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             if(specific_imre) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             cycle
          end if
          call dealloc(errmsg2)
          !
          ! if there are specific components for each path, the next line contains those
          if(specific_comp) then
             read(lu,"(a)",iostat=ios) line; iline = iline+1
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                     "; expected number of components ncomp as first value on that line, followed by ncomp components"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             read(line,*,iostat=ios) ncomp
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read number of components as first value on line",iline,", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(ncomp .le. 0) then
                close(lu)
                write(errstr,*) "number of components ",ncomp," (first value on line",iline,") must be positive"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(allocated(comp)) deallocate(comp)
             allocate(comp_tmp(ncomp))
             read(line,*,iostat=ios) ncomp,comp_tmp
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read ",ncomp," components starting from second value of line",iline,", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             ! check if there are invalid components
             nvalid = ncomp
             allocate(valid(nvalid))
             do jcomp = 1,ncomp
                valid(jcomp) = validComponent(comp_tmp(jcomp))
                if(.not.valid(jcomp)) then
                   nvalid = nvalid - 1
                   write(errstr,*) jcomp,"'th component '",trim(comp_tmp(jcomp)),"' on line ",iline, &
                        " is not valid, hence it is excluded! Valid components are: "//trim(all_valid_components)
                   call add(errmsg,1,trim(errstr),myname)
                endif
             enddo ! jcomp
             ! if there were invalid components, exclude them. if there are no valid components, cycle to next path
             if(nvalid<ncomp) then
                if(nvalid .le. 0) then
                   write(errstr,*) "no valid components on line ",iline,", hence excluding this path!"
                   call add(errmsg,1,trim(errstr),myname)
                   ! skip all remaining lines which belong to this path and cycle
                   if(specific_freq) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
                   if(specific_imre) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
                   deallocate(comp_tmp,valid)
                   cycle
                end if
                allocate(comp(nvalid))
                comp = pack(comp_tmp,valid)
             else
                allocate(comp(ncomp))
                comp = comp_tmp
             end if
             deallocate(comp_tmp,valid)
          endif ! specific_comp
          !
          ! if there are specific frequencies for each path, the next line contains those
          if(specific_freq) then
             read(lu,"(a)",iostat=ios) line; iline = iline+1
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                     "; expected number of frequencies nfreq as first value on that line, followed by nfreq frequency indices"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             read(line,*,iostat=ios) nfreq
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read number of frequency indices as first value on line",iline,", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(nfreq .le. 0) then
                close(lu)
                write(errstr,*) "number of frequency indices ",nfreq," (first value on line",iline,") must be positive"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(allocated(ifreq)) deallocate(ifreq)
             allocate(ifreq_tmp(nfreq))
             read(line,*,iostat=ios) nfreq,ifreq_tmp
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read ",nfreq," frequency indices starting from second value of line",iline,&
                     ", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             ! check if there are invalid frequencies
             nvalid = nfreq
             allocate(valid(nvalid))
             do jfreq = 1,nfreq
                valid(jfreq) = any(ifreq_valid == ifreq_tmp(jfreq))
                if(.not.valid(jfreq)) then
                   nvalid = nvalid - 1
                   write(errstr,*) jfreq,"'th frequency index '",ifreq_tmp(jfreq),"' on line ",iline, &
                        " is invalid (valid indices: ",ifreq_valid,"), hence it is excluded!"
                   call add(errmsg,1,trim(errstr),myname)
                endif
             enddo ! jfreq
             ! if there were invalid frequencies, exclude them. if there are no valid frequencies, cycle to next path
             if(nvalid<nfreq) then
                if(nvalid .le. 0) then
                   write(errstr,*) "no valid frequency index on line ",iline,", hence excluding this path!"
                   call add(errmsg,1,trim(errstr),myname)
                   ! skip all remaining lines which belong to this path and cycle
                   if(specific_imre) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
                   deallocate(ifreq_tmp,valid)
                   cycle
                end if
                allocate(ifreq(nvalid))
                ifreq = pack(ifreq_tmp,valid)
             else
                allocate(ifreq(nfreq))
                ifreq = ifreq_tmp
             end if
             deallocate(ifreq_tmp,valid)
          endif ! specific_freq
          !
          ! if there are specific imre's for each path, the next line contains those
          if(specific_imre) then
             read(lu,"(a)",iostat=ios) line; iline = iline+1
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                     "; expected number of imaginary/real part values nimre as first value on that line, "//&
                     "followed by nimre such values"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             read(line,*,iostat=ios) nimre
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read number of imaginary/real part values as first value on line",iline,&
                     ", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(nimre .le. 0) then
                close(lu)
                write(errstr,*) "number of imaginary/real part values ",nimre," (first value on line",iline,") must be positive"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(allocated(imre)) deallocate(imre)
             allocate(imre_tmp(nimre))
             read(line,*,iostat=ios) nimre,imre_tmp
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read ",nimre," frequency indices starting from second value of line",iline,&
                     ", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             ! check if there are invalid imre values
             nvalid = nimre
             allocate(valid(nvalid))
             do jimre = 1,nimre
                select case(imre_tmp(jimre))
                case('im','re') 
                   valid(jimre) = .true.
                case default
                   nvalid = nvalid - 1
                   valid(jimre) = .false.
                   write(errstr,*) jimre,"'th imaginary/real part value '",trim(imre_tmp(jimre)),&
                        "' on line ",iline," is not supported, hence it is excluded! Must be either 'im' or 're'."
                   call add(errmsg,1,trim(errstr),myname)
                end select
             enddo ! jimre
             ! if there were invalid imre values, exclude them. if there are no valid imre values, cycle to next path
             if(nvalid<nimre) then
                if(nvalid .le. 0) then
                   write(errstr,*) "no valid imaginary/real part value on line ",iline,", hence excluding this path!"
                   call add(errmsg,1,trim(errstr),myname)
                   ! as there are no remaining lines which belong to this path, just cycle
                   deallocate(imre_tmp,valid)
                   cycle
                end if
                allocate(imre(nvalid))
                imre = pack(imre_tmp,valid)
             else
                allocate(imre(nimre))
                imre = imre_tmp
             end if
             deallocate(imre_tmp,valid)
          endif ! specific_imre
          !
          ! now, for this path all information on components, frequencies, imre's are gathered, so finally add data samples
          call addDataSamplesDataModelSpaceInfo(this,(/evid/),(/staname/),comp,ifreq,imre,all_combinations=.true.)
          if(specific_comp) deallocate(comp)
          if(specific_freq) deallocate(ifreq)
          if(specific_imre) deallocate(imre)
       enddo ! jpath
       !
       ! if arrays comp,ifreq,imre were allocated above (before jpath-loop), deallocate here
       if(.not.specific_comp) deallocate(comp)
       if(.not.specific_freq) deallocate(ifreq)
       if(.not.specific_imre) deallocate(imre)
    else ! specific_paths
       call addDataSamplesDataModelSpaceInfo(this,evids,stanames,comp,ifreq,imre,all_combinations=.true.)
       deallocate(evids,stanames,comp,ifreq,imre)
    endif ! specific_paths
!
    close(lu)
    write(errstr,*) this%ndata-ndata_before," data samples were added"
    call add(errmsg,0,trim(errstr),myname)
  end subroutine createDataSamplesFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief create model parameters from the MODEL PARAMETERS block of the data model space input file
!! \details This routine searches for the line "MODEL PARAMETERS" in the given file and reads in the
!!  block following this line. 
!!  Format of 'MODEL PARAMETERS' block:
!!     line of form: 'INVERSION_GRID_CELLS value', where value is either 'ALL' (all inversion 
!!                   grid cells are taken) or 'SPECIFIC' (specific definition of set of invgrid cells following below)
!!     line of form: 'PARAMETERS value', where value is either 'ALL' (all inversion grid cells 
!!                   are taken) or 'SPECIFIC' ('SPECIFIC' only allowed if 'INVERSION_GRID_CELLS SPECIFIC')
!!     If 'PARAMETERS SPECIFIC' the next line must be of form 'nparam param_1 ... param_n', 
!!        defining the parameters of all inversion grid cells (assumed to be parameters of given parametrization).
!!     If 'INVERSION_GRID_CELLS SPECIFIC', the following line must contain the number of cells ncell which should be taken, 
!!        followed by ncell blocks of lines, each defining an inversion grid cell. These blocks consist of at least 
!!        one line containing the inversion grid cell index (if 'PARAMETERS ALL'). If 'PARAMETERS SPECIFIC', an additional 
!!        line of the form 'nparam param_1 ... param_n' defines the parameters to be used for this 
!!        specific inversion grid cell
!! \param this data model space info
!! \param parametrization all parameters are assumed to be of this parametrization 
!! \param ntot_invgrid total number of inversion grid cells to check valid range of cell indices in file
!! \param intw integration weights object which is used to check if inversion grid cells are empty (hence invalid)
!! \param filename file name to use
!! \param lu file unit, with which the input file is already open
!! \param errmsg error message
!
  subroutine createModelParametersFromFileDataModelSpaceInfo(this,parametrization,ntot_invgrid,intw,filename,lu,errmsg)
    type (data_model_space_info) :: this
    character(len=*) :: parametrization
    integer :: ntot_invgrid
    type (integration_weights) :: intw
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=47) :: myname = 'createModelParametersFromFileDataModelSpaceInfo'
    character(len=500) :: line
    integer :: ios,iline,nvalid,nparam_before
    character(len=40) :: key,val
    integer :: ncell,icell,jcell,nparam,jparam
    character(len=character_length_param), dimension(:), allocatable :: param,param_tmp
    logical, dimension(:), allocatable :: valid
    logical :: specific_cells,specific_param,line_model_parameters_found,raised_warning_invalid_cells
    integer, dimension(:), pointer :: indx
!
    call addTrace(errmsg,myname)
    iline = 0
    nparam_before = this%nparam
!
    ! first of all make sure that the current data_mode_space_info object does not already
    ! contain any model parameters of a parametrization other than requested here
    if(trim(this%parametrization) /= '' .and. trim(this%parametrization) /= trim(parametrization)) then
       call add(errmsg,2,"the model space already contains model parameters of parametrization '"//&
            trim(this%parametrization)//"', which is different from the parametrization '"//trim(parametrization)//&
            "' of the model parameters to be added. There are no multiple parametrizations allowed.",myname)
       return
    end if
    this%parametrization = parametrization
!
    open(unit=lu,file=trim(filename),status='old',form='formatted',action='read',iostat=ios)
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not open file '"//trim(filename)//"', iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    !
    ! TODO in the future:
    ! if introducing some header to dataModelSpaceInfo files, read in header here, and react accordingly
    !
    ! parse through whole file until line "MODEL PARAMETERS" and start reading in block from there
    line_model_parameters_found = .false.
    do while(ios==0)
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) exit
       if(line == "MODEL PARAMETERS") then
          line_model_parameters_found = .true.
          exit
       endif
    enddo
    if(.not.line_model_parameters_found) then
       write(errstr,*) "could not find line 'MODEL PARAMETERS' in file '"//trim(filename)//"', searching until line ",iline-1
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    write(errstr,*) "start reading 'MODEL PARAMETERS' block from line ",iline," of file '"//trim(filename)//"'"
    call add(errmsg,0,trim(errstr),myname)
!
    ! INVERSION_GRID_CELLS
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword INVERSION_GRID_CELLS followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'INVERSION_GRID_CELLS') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline,&
            " in MODEL PARAMETERS block not supported. 'INVERSION_GRID_CELLS' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('ALL')
       specific_cells = .false.
    case('SPECIFIC')
       specific_cells = .true.
    case default
       write(errstr,*) "'INVERSION_GRID_CELLS' specification '"//trim(val)//"' on line ",iline,&
            " in MODEL PARAMETERS block not supported. 'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val	
!
    ! PARAMETERS
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword PARAMETERS followed by ALL or SPECIFIC"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'PARAMETERS') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in MODEL PARAMETERS block not supported. 'PARAMETERS' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('ALL')
       specific_param = .false.
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of parameters nparam as first value on that line, followed by nparam parameter names"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) nparam
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of parameters as first value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(nparam .le. 0) then
          close(lu)
          write(errstr,*) "number of parameters ",nparam," (first value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       allocate(param_tmp(nparam))
       read(line,*,iostat=ios) nparam,param_tmp
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read ",nparam," parameter names starting from second value of line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       !write(*,*) iline,"#"//trim(line)//"#"
       ! check if there are invalid parameters
       nvalid = nparam
       allocate(valid(nvalid))
       do jparam = 1,nparam
          valid(jparam) = validParamModelParametrization(parametrization,param_tmp(jparam))
          if(.not.valid(jparam)) then
             nvalid = nvalid - 1
             write(errstr,*) jparam,"'th parameter '",trim(param_tmp(jparam)),"' on line ",iline, &
                  " is not a valid parameter of parametrization '",trim(parametrization),&
                  "', hence it is excluded! Valid parametrizations (parameters) are: ",all_valid_pmtrz_param
             call add(errmsg,1,trim(errstr),myname)
          endif
       enddo ! jparam
       ! if there were invalid parameters, exclude them. if there were no valid parameters, return
       if(nvalid<nparam) then
          if(nvalid .le. 0) then
             write(errstr,*) "no valid parameter of parametrization '",trim(parametrization),"' on line ",iline,&
                  ", hence no model parameters are created! Valid parametrizations (parameters) are: ",all_valid_pmtrz_param
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          end if
          allocate(param(nvalid))
          param = pack(param_tmp,valid)
       else
          allocate(param(nparam))
          param = param_tmp
       end if
       deallocate(param_tmp,valid)
!
    case('SPECIFIC') ! val
!
       if(.not.specific_cells) then
          call add(errmsg,2,"'PARAMETERS SPECIFIC' is only allowed in case of 'INVERSION_GRID_CELLS SPECIFIC'",myname)
          close(lu)
          return
       endif
       specific_param = .true.
!
    case default ! val
!
       write(errstr,*) "'PARAMETERS' specification '"//trim(val)//"' on line ",iline," in MODEL PARAMETERS block not supported."//&
            " 'ALL' or 'SPECIFIC' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select ! val
!
! now iterate over lines defining specific invgrid cells (if there are any)
!
    if(specific_cells) then
       raised_warning_invalid_cells = .false.
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read line",iline,", iostat = ",ios,&
               "; expected number of inversion grid cells as first (and only) value on that line"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       read(line,*,iostat=ios) ncell
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read number of inversion grid cells as first (and only) value on line",iline,", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       if(ncell .le. 0) then
          close(lu)
          write(errstr,*) "number of invgrid cells ",ncell," (first (and only) value on line",iline,") must be positive"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       do jcell = 1,ncell
          read(lu,"(a)",iostat=ios) line; iline = iline+1
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                  "; expected ",jcell,"'th invgrid grid cell index (out of ",ncell,") on that line"
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          read(line,*,iostat=ios) icell
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read ",jcell,"'th invgrid grid cell index (out of ",ncell,") on line",&
                  iline,", iostat = ",ios
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          if((icell < 1 .or. icell > ntot_invgrid) .or. emptyCell(intw,icell)) then
             if(.not.raised_warning_invalid_cells) then
                raised_warning_invalid_cells = .true.
                write(errstr,*) "there are specific inversion grid cell indices which are either out of range ",&
                     "(total number of inversion grid cells is ",ntot_invgrid,&
                     ") or for which the cell is empty, hence they are excluded!"
                call add(errmsg,1,trim(errstr),myname)
             end if
             ! before going to the next invgrid cell, we need to skip the next line in case of specific_param
             if(specific_param) then; read(lu,"(a)",iostat=ios) line; iline=iline+1; endif
             cycle
          end if
          if(specific_param) then
             read(lu,"(a)",iostat=ios) line; iline = iline+1
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                     "; expected number of parameters nparam as first value on that line, followed by nparam parameter names"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             read(line,*,iostat=ios) nparam
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read number of parameters as first value on line",iline,", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             if(nparam .le. 0) then
                close(lu)
                write(errstr,*) "number of parameters ",nparam," (first value on line",iline,") must be positive"
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             allocate(param_tmp(nparam))
             read(line,*,iostat=ios) nparam,param_tmp
             if(ios/=0) then
                close(lu)
                write(errstr,*) "could not read ",nparam," parameter names starting from second value of line",&
                     iline,", iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                return
             end if
             ! check if there are invalid parameters
             nvalid = nparam
             allocate(valid(nvalid))
             do jparam = 1,nparam
                valid(jparam) = validParamModelParametrization(parametrization,param_tmp(jparam))
                if(.not.valid(jparam)) then
                   nvalid = nvalid - 1
                   write(errstr,*) jparam,"'th parameter '",trim(param_tmp(jparam)),"' on line ",iline, &
                        " is not a valid parameter of parametrization '",trim(parametrization),&
                        "', hence it is excluded! Valid parametrizations (parameters) are: ",all_valid_pmtrz_param
                   call add(errmsg,1,trim(errstr),myname)
                endif
             enddo ! jparam
             ! if there were invalid parameters, exclude them. if there were no valid parameters, cycle
             if(nvalid<nparam) then
                if(nvalid .le. 0) then
                   write(errstr,*) "no valid parameter of parametrization '",trim(parametrization),"' on line ",iline,&
                        ", hence excluding this path!"
                   call add(errmsg,1,trim(errstr),myname)
                   ! go to next invgrid cell
                   deallocate(param_tmp,valid)
                   cycle
                end if
                allocate(param(nvalid))
                param = pack(param_tmp,valid)
             else
                allocate(param(nparam))
                param = param_tmp
             end if
             deallocate(param_tmp,valid)
          endif ! specific_param
          !
          ! now, for this cell all information on parameters are gathered, so finally add model parameters
          call addModelParametersDataModelSpaceInfo(this,param,(/icell/),all_combinations=.true.)
          if(specific_param) deallocate(param)
       enddo ! jcell
       !
       ! if array param was allocated above (before jcell-loop), deallocate here
       if(.not.specific_param) deallocate(param)
    else ! specific_cells
       indx => getFilledCells(intw)
       if(.not.associated(indx)) then
          call add(errmsg,2,'there are only empty cells in inversion grid, says the integration weights object',myname)
          deallocate(indx)
          close(lu)
          return
       end if
       call addModelParametersDataModelSpaceInfo(this,param,indx,all_combinations=.true.)
       deallocate(indx)
       deallocate(param)
    endif ! specific_cells
!
    close(lu)
    write(errstr,*) this%nparam-nparam_before," model parameters were added"
    call add(errmsg,0,trim(errstr),myname)
  end subroutine createModelParametersFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return paths contained in this data_model_space_info object
!! \param 
!
  function getPathsDataModelSpaceInfo(this) result(paths)
    type (data_model_space_info) :: this
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    ! local
    character(len=character_length_evid), dimension(:), allocatable :: evid
    character(len=character_length_staname), dimension(:), allocatable :: staname
    integer :: idata,npaths
    logical, dimension(:), allocatable :: mask_paths
!
    nullify(paths)
    if(this%ndata == 0) return
!
    ! make local copies of evid and staname to modify
    allocate(evid(this%ndata),staname(this%ndata))
    evid = this%evid
    staname = this%staname
!
    do idata = 1,this%ndata
       if(evid(idata)/='' .and. staname(idata)/='') then
          where(evid(idata+1:this%ndata)==evid(idata) .and. staname(idata+1:this%ndata)==staname(idata))
             evid(idata+1:this%ndata) = ''
             staname(idata+1:this%ndata) = ''
          end where
       end if
    end do ! idata
!
    allocate(mask_paths(this%ndata))
    mask_paths = evid/= '' .and. staname/=''
    npaths = count(mask_paths)
    if(npaths==0) goto 1
!
    allocate(paths(2,npaths))
    paths(1,:) = pack(evid,mask_paths)
    paths(2,:) = pack(staname,mask_paths)
!
1   deallocate(mask_paths,evid,staname)
  end function getPathsDataModelSpaceInfo
!
!! FS FS
! FUNCTION:
! write data model space as text file (with MODEL PARAMETERS SPECIFIC, DATA SAMPLES SPECIFIC)
!
! THE TEXT FILE FORMAT IS THE ONLY INTERFACE TO THIS MODULE!
! later on (if needed), people can write routines like "addDataSampleFromValues" (having optional "allCombinations" flag
! which is handled in the same way as in modules addDataSamplesDataModelSpaceInfo,addModelParametersDataModelSpaceInfo)
! where incoming values are checked (in the same way as in routines "add*FromFile") and inside which routines 
! addDataSamplesDataModelSpaceInfo,addModelParametersDataModelSpaceInfo are called. 
!! FS FS
!
!------------------------------------------------------------------------
!> \brief deallocate data_model_space_info object
!! \param this data_model_space_info object
!
  subroutine deallocateDataModelSpaceInfo(this)
    type (data_model_space_info) :: this
    if(associated(this%evid)) deallocate(this%evid)
    if(associated(this%staname)) deallocate(this%staname)
    if(associated(this%comp)) deallocate(this%comp)
    if(associated(this%ifreq)) deallocate(this%ifreq)
    if(associated(this%imre)) deallocate(this%imre)
    this%ndata = 0
    if(associated(this%param)) deallocate(this%param)
    if(associated(this%cell)) deallocate(this%cell)
    this%nparam = 0
    this%parametrization = ''
  end subroutine deallocateDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get total number of data samples
!! \param this data model space info
!! \param ndata this%ndata
!! \return number of data samples 
!
  function getNdataDataModelSpaceInfo(this) result(ndata)
    type (data_model_space_info), intent(in) :: this
    integer :: ndata
    ndata = this%ndata
  end function getNdataDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get total number of model parameters
!! \param this data model space info
!! \param nparam this%nparam
!! \return number of model parameters
!
  function getNparamDataModelSpaceInfo(this) result(nparam)
    type (data_model_space_info), intent(in) :: this
    integer :: nparam
    nparam = this%nparam
  end function getNparamDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get event ID of i-th data sample
!! \param this data model space info
!! \param i index of data sample
!! \param evid event ID this%data_space_info(i)%evid
!! \return event ID of i-th data sample
!
  function getEvidDataSampleDataModelSpaceInfo(this,i) result(evid)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=character_length_evid) :: evid
    if(i<1 .or. i>this%ndata) then
       evid = ''
    else
       evid = this%evid(i)
    endif
  end function getEvidDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get all event IDs as pointer to array
!! \param this data model space info
!! \param evid pointer to event ID array
!! \return event IDs as pointer to array
!
  function getArrayEvidDataSamplesDataModelSpaceInfo(this) result(evid)
    type (data_model_space_info), intent(in) :: this
    character(len=character_length_evid), dimension(:), pointer :: evid
    evid => this%evid
  end function getArrayEvidDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get station name of i-th data sample
!! \param this data model space info
!! \param i index of data sample
!! \param staname station name this%data_space_info(i)%staname
!! \return station name of i-th data sample
!
  function getStanameDataSampleDataModelSpaceInfo(this,i) result(staname)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=character_length_staname) :: staname
    if(i<1 .or. i>this%ndata) then
       staname = ''
    else
       staname = this%staname(i)
    endif
  end function getStanameDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get all station names as pointer to array
!! \param this data model space info
!! \param staname pointer to station name array
!! \return station names as pointer to array
!
  function getArrayStanameDataSamplesDataModelSpaceInfo(this) result(staname)
    type (data_model_space_info), intent(in) :: this
    character(len=character_length_staname), dimension(:), pointer :: staname
    staname => this%staname
  end function getArrayStanameDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get station component of i-th data sample
!! \param this data model space info
!! \param i index of data sample
!! \param comp component this%data_space_info(i)%comp
!! \return station component of i-th data sample
!
  function getCompDataSampleDataModelSpaceInfo(this,i) result(comp)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=character_length_component) :: comp
    if(i<1 .or. i>this%ndata) then
       comp = ''
    else
       comp = this%comp(i)
    endif
  end function getCompDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get all different occurring station components
!! \param this data model space info
!! \param comp component this%data_space_info(i)%comp
!! \return station component of i-th data sample
!
  function getAllDifferentCompDataSamplesDataModelSpaceInfo(this) result(comp)
    type (data_model_space_info), intent(in) :: this
    character(len=character_length_component), dimension(:), pointer :: comp
    character(len=character_length_component), dimension(:), allocatable :: comp_tmp
    integer :: idata,ncomp
    logical, dimension(:), allocatable :: mask
!
    nullify(comp)
    if(this%ndata == 0) return
!
    ! make a local copy of this%comp to modify
    allocate(comp_tmp(this%ndata)); comp_tmp = this%comp
!
    do idata = 1,this%ndata
       if(comp_tmp(idata)/='') then
          where(comp_tmp(idata+1:this%ndata)==comp_tmp(idata))
             comp_tmp(idata+1:this%ndata) = ''
          end where
       end if
    end do ! idata
!
    allocate(mask(this%ndata))
    mask = comp_tmp/=''
    ncomp = count(mask)
    if (ncomp==0) goto 1
!
    allocate(comp(ncomp))
    comp = pack(comp_tmp,mask)
!
1   deallocate(mask,comp_tmp)
  end function getAllDifferentCompDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get frequency index of i-th data sample
!! \param this data model space info
!! \param i index of data sample
!! \param ifreq frequency index this%data_space_info(i)%ifreq
!! \return frequency index of i-th data sample
!
  function getIfreqDataSampleDataModelSpaceInfo(this,i) result(ifreq)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    integer :: ifreq
    if(i<1 .or. i>this%ndata) then
       ifreq = -1
    else
       ifreq = this%ifreq(i)
    endif
  end function getIfreqDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get imre of i-th data sample
!! \param this data model space info
!! \param i index of data sample
!! \param imre 'im' or 're' as in this%data_space_info(i)%imre
!! \return imre value of i-th data sample
!
  function getImreDataSampleDataModelSpaceInfo(this,i) result(imre)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=2) :: imre
    if(i<1 .or. i>this%ndata) then
       imre = ''
    else
       imre = this%imre(i)
    endif
  end function getImreDataSampleDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get parametrization of this model space
!! \param this data model space info
!! \param pmtrz this%parametrization
!! \return parametrization of this model space
!
  function getParametrizationDataModelSpaceInfo(this) result(pmtrz)
    type (data_model_space_info), intent(in) :: this
    character(len=character_length_pmtrz) :: pmtrz
    pmtrz = this%parametrization
  end function getParametrizationDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get parameter of i-th model parameter
!! \param this data model space info
!! \param i index of model parameter
!! \param param this%model_space_info(i)%param
!! \return parameter of i-th model parameter
!
  function getParamModelParameterDataModelSpaceInfo(this,i) result(param)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=character_length_param) :: param
    if(i<1 .or. i>this%nparam) then
       param = ''
    else
       param = this%param(i)
    endif
  end function getParamModelParameterDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief inversion grid cell index of i-th model parameter
!! \param this data model space info
!! \param i index of model parameter
!! \param cell invgrid cell index this%model_space_info(i)%cell
!! \return inversion grid cell index of i-th model parameter
!
  function getCellModelParameterDataModelSpaceInfo(this,i) result(cell)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    integer :: cell
    if(i<1 .or. i>this%nparam) then
       cell = -1
    else
       cell = this%cell(i)
    endif
  end function getCellModelParameterDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return indices of data samples which have certain properties
!! \details On Enter: optional arrays evid,staname,comp,ifreq,imre define properties of data 
!!  samples for which the indices should be returned, whereby any properties in the same array (e.g. two
!!  evenIDs S001,S002 contained in array evid) are connected by an "or" condition (i.e. the index of a
!!  data sample is returned, if it has eventID S001 OR S002) and properties in different arrays (e.g.
!!  station name R001 in array staname and component CY in array comp) are connected by an "and" condition (i.e.
!!  the index of a data sample is returned, if it has station name R001 AND component CY).
!!  If such an optional array is not present, or points to null() (e.g. ifreq not present), this property is not
!!  used to condition the data samples for which indices are returned (i.e. indices for data samples with arbitrary
!!  frequency indices are returned).
!!  Optionally, array indx_in will a priorily restrict the indices of data samples among which are searched for the 
!!  given properties (e.g. if you want to look for specific data samples having certain properties among the first 
!!  100 data samples, set indx_in = (/ (i,i=1,100) /) . The indices returned in array indx, of course will then be 
!!  a subset of array indx_in)
!!  On Exit: index array indx contains all indices of data samples contained in this data_model_space_info object, which
!!  satisfy the conditions given by arrays evid,staname,comp,ifreq,imre,indx_in. If there are no data samples which 
!!  satisfy the conditions, it points to nul(). For any array of evid,staname,comp,ifreq,imre, 
!!  which is present on enter, the array is reallocated inside this subroutine and contains on exit the respective values
!!  corresponding to indices in array indx. It is not touched, however, and points to the same array as it did on enter, if
!!  indx => null() on exit (i.e. no valid indices found)
!!  For proper reallocation to work, the lengths of the incoming characters MUST be:
!!    character_length_evid for array evid
!!    character_length_staname for array staname
!!    character_length_component for array comp
!!    2 for array imre
!!  Example: 
!!    On enter: array evid = ("S001","S002")
!!    On exit: indx = (1,2,5,6,7) and evid = ("S001","S001","S002","S002","S002")
!!    This means that data samples 1 and 2 have eventID S001, and data samples 5,6,7 have eventID S002
!! \param this data model space
!! \param evid optional pointer to array of event IDs
!! \param staname optional pointer to array of station names
!! \param comp optional pointer to array of components
!! \param ifreq optional pointer to array of frequency indices
!! \param imre optional pointer to array containing values 'im' or 're' indicating imaginary or real part of the data sample
!! \param indx_in optional pointer to index array of data samples among which will be searched only
!! \param indx array of indices of data samples which satisfy all above constraints (if present)
!! \return array of indices of data samples which satisfy given properties
!
  function getIndicesDataSamplesDataModelSpaceInfo(this,evid,staname,comp,ifreq,imre,indx_in) result(indx)
    type (data_model_space_info) :: this
    character(len=*), dimension(:), pointer, optional :: evid,staname,comp,imre
    integer, dimension(:), pointer, optional :: ifreq,indx_in
    integer, dimension(:), pointer :: indx
    ! local
    integer :: nindx_in,nindx_search,nindx_found
    integer :: nevid,nstaname,ncomp,nfreq,nimre
    integer :: j
    integer, dimension(:), allocatable :: indx_in_valid
    logical, dimension(:), allocatable :: valid,valid_tmp
    logical :: use_indx_in,search_evid,search_staname,search_comp,search_ifreq,search_imre
    character(len=character_length_evid), dimension(:), pointer :: evid_tmp
    character(len=character_length_staname), dimension(:), pointer :: staname_tmp
    character(len=character_length_component), dimension(:), pointer :: comp_tmp
    integer, dimension(:), pointer :: ifreq_tmp
    character(len=2), dimension(:), pointer :: imre_tmp
!
    nullify(indx)
!
    ! define which properties condition the output
    search_evid = .false.
    if(present(evid)) then
       if(associated(evid)) then
          search_evid = .true.
          nevid = size(evid)
       end if
    end if
    search_staname = .false.
    if(present(staname)) then
       if(associated(staname)) then
          search_staname = .true.
          nstaname = size(staname)
       end if
    end if
    search_comp = .false.
    if(present(comp)) then
       if(associated(comp)) then
          search_comp = .true.
          ncomp = size(comp)
       end if
    end if
    search_ifreq = .false.
    if(present(ifreq)) then
       if(associated(ifreq)) then
          search_ifreq = .true.
          nfreq = size(ifreq)
       end if
    end if
    search_imre = .false.
    if(present(imre)) then
       if(associated(imre)) then
          search_imre = .true.
          nimre = size(imre)
       end if
    end if
!
    ! define indices among which will be searched for data samples having the given properties
    use_indx_in = .false.
    nindx_search = this%ndata
    if(present(indx_in)) then
       if(associated(indx_in)) then
          use_indx_in = .true.
          ! first check if there are indices out of range in incoming array indx_in and choose as indx_in_valid all valid indices
          nindx_in = size(indx_in)
          if(nindx_in == 0) return
          allocate(valid(nindx_in))
          valid = (indx_in>0).and.(indx_in.le.this%ndata)
          nindx_search = count(valid)
          if(nindx_search .le. 0) return
          allocate(indx_in_valid(nindx_search))
          if(nindx_search == nindx_in) then
             indx_in_valid = indx_in
          else
             indx_in_valid = pack(indx_in,valid)
          end if
          deallocate(valid)
       end if ! associated(indx_in)
    end if ! present(indx_in)
!
    ! find valid indices to return, i.e. find all data samples which satisfy the given conditions (if any)
    allocate(valid(nindx_search))
    valid = .true.
!
    ! find all data samples with required event IDs
    if(search_evid) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,nevid
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%evid(indx_in_valid) == evid(j) )
          else
             valid_tmp = valid_tmp .or. ( this%evid == evid(j) )
          end if
       end do ! j
       ! update all valid data samples
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no data samples satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all data samples with required station name
    if(search_staname) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,nstaname
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%staname(indx_in_valid) == staname(j) )
          else
             valid_tmp = valid_tmp .or. ( this%staname == staname(j) )
          end if
       end do ! j
       ! update all valid data samples
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no data samples satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all data samples with required components
    if(search_comp) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,ncomp
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%comp(indx_in_valid) == comp(j) )
          else
             valid_tmp = valid_tmp .or. ( this%comp == comp(j) )
          end if
       end do ! j
       ! update all valid data samples
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no data samples satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all data samples with required frequency indices
    if(search_ifreq) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,nfreq
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%ifreq(indx_in_valid) == ifreq(j) )
          else
             valid_tmp = valid_tmp .or. ( this%ifreq == ifreq(j) )
          end if
       end do ! j
       ! update all valid data samples
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no data samples satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all data samples with required im/re values
    if(search_imre) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,nimre
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%imre(indx_in_valid) == imre(j) )
          else
             valid_tmp = valid_tmp .or. ( this%imre == imre(j) )
          end if
       end do ! j
       ! update all valid data samples
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if there are no data samples satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! prepare all valid indices for return
    nindx_found = count(valid)
    if(nindx_found .le. 0) return ! actually this cannot be, because of lines "if(.not.any(valid)) goto 1" above
    allocate(indx(nindx_found))
    if(use_indx_in) then
       indx = pack(indx_in_valid,valid)
    else
       indx = pack( (/ (j,j=1,this%ndata) /),valid)
    end if
!
    ! reallocate property arrays for return, if were present on enter
    ! ALSO REALLOCATE IN CASE OF INCOMING NULLIFIED POINTERS!, this way you can get the values of some property even 
    ! if this property was not used to constrain the index search
    if(present(evid)) then
       if(associated(evid)) deallocate(evid)
       allocate(evid_tmp(nindx_found))
       evid_tmp = this%evid(indx)
       evid => evid_tmp
    end if
    if(present(staname)) then
       if(associated(staname)) deallocate(staname)
       allocate(staname_tmp(nindx_found))
       staname_tmp = this%staname(indx)
       staname => staname_tmp
    end if
    if(present(comp)) then
       if(associated(comp)) deallocate(comp)
       allocate(comp_tmp(nindx_found))
       comp_tmp = this%comp(indx)
       comp => comp_tmp
    end if
    if(present(ifreq)) then
       if(associated(ifreq)) deallocate(ifreq)
       allocate(ifreq_tmp(nindx_found))
       ifreq_tmp = this%ifreq(indx)
       ifreq => ifreq_tmp
    end if
    if(present(imre)) then
       if(associated(imre)) deallocate(imre)
       allocate(imre_tmp(nindx_found))
       imre_tmp = this%imre(indx)
       imre => imre_tmp
    end if
!
1   if(allocated(indx_in_valid)) deallocate(indx_in_valid)
    if(allocated(valid)) deallocate(valid)
  end function getIndicesDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return indices of data samples which have certain properties
!! \details handling is completely analogous to routine getIndicesDataSamplesDataModelSpaceInfo
!! \param this data model space
!! \param param optional pointer to array of names of model parameters
!! \param cell optional pointer to array of inversion grid cell indices
!! \param indx_in optional pointer to index array of model parameters among which will be searched only
!! \return array of indices of model parameters which satisfy given properties
!
  function getIndicesModelParametersDataModelSpaceInfo(this,param,cell,indx_in) result(indx)
    type (data_model_space_info) :: this
    character(len=*), dimension(:), pointer, optional :: param
    integer, dimension(:), pointer, optional :: cell,indx_in
    integer, dimension(:), pointer :: indx
    ! local
    integer :: nindx_in,nindx_search,nindx_found
    integer :: nparam,ncell
    integer :: j
    integer, dimension(:), allocatable :: indx_in_valid
    logical, dimension(:), allocatable :: valid,valid_tmp
    logical :: use_indx_in,search_param,search_cell
    character(len=character_length_param), dimension(:), pointer :: param_tmp
    integer, dimension(:), pointer :: cell_tmp
!
    nullify(indx)
!
    ! define which properties condition the output
    search_param = .false.
    if(present(param)) then
       if(associated(param)) then
          search_param = .true.
          nparam = size(param)
       end if
    end if
    search_cell = .false.
    if(present(cell)) then
       if(associated(cell)) then
          search_cell = .true.
          ncell = size(cell)
       end if
    end if
!
    ! define indices among which will be searched for model parameters having the given properties
    use_indx_in = .false.
    nindx_search = this%nparam
    if(present(indx_in)) then
       if(associated(indx_in)) then
          use_indx_in = .true.
          ! first check if there are indices out of range in incoming array indx_in and choose as indx_in_valid all valid indices
          nindx_in = size(indx_in)
          if(nindx_in == 0) return
          allocate(valid(nindx_in))
          valid = (indx_in>0).and.(indx_in.le.this%nparam)
          nindx_search = count(valid)
          if(nindx_search .le. 0) return
          allocate(indx_in_valid(nindx_search))
          if(nindx_search == nindx_in) then
             indx_in_valid = indx_in
          else
             indx_in_valid = pack(indx_in,valid)
          end if
          deallocate(valid)
       end if ! associated(indx_in)
    end if ! present(indx_in)
!
    ! find valid indices to return, i.e. find all model parameters which satisfy the given conditions (if any)
    allocate(valid(nindx_search))
    valid = .true.
!
    ! find all model parameters with required parameter name
    if(search_param) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,nparam
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%param(indx_in_valid) == param(j) )
          else
             valid_tmp = valid_tmp .or. ( this%param == param(j) )
          end if
       end do ! j
       ! update all valid model parameters
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no model parameters satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all model parameters with required inversion grid cell index
    if(search_cell) then
       allocate(valid_tmp(nindx_search))
       valid_tmp = .false.
       do j = 1,ncell
          if(use_indx_in) then
             valid_tmp = valid_tmp .or. ( this%cell(indx_in_valid) == cell(j) )
          else
             valid_tmp = valid_tmp .or. ( this%cell == cell(j) )
          end if
       end do ! j
       ! update all valid model parameters
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if there are no model parameters satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! prepare all valid indices for return
    nindx_found = count(valid)
    if(nindx_found .le. 0) return ! actually this cannot be, because of lines "if(.not.any(valid)) goto 1" above
    allocate(indx(nindx_found))
    if(use_indx_in) then
       indx = pack(indx_in_valid,valid)
    else
       indx = pack( (/ (j,j=1,this%nparam) /),valid)
    end if
!
    ! reallocate property arrays for return, if were present on enter
    ! ALSO REALLOCATE IN CASE OF INCOMING NULLIFIED POINTERS!, this way you can get the values of some property even 
    ! if this property was not used to constrain the index search
    if(present(param)) then
       if(associated(param)) deallocate(param)
       allocate(param_tmp(nindx_found))
       param_tmp = this%param(indx)
       param => param_tmp
    end if
    if(present(cell)) then
       if(associated(cell)) deallocate(cell)
       allocate(cell_tmp(nindx_found))
       cell_tmp = this%cell(indx)
       cell => cell_tmp
    end if
!
1   if(allocated(indx_in_valid)) deallocate(indx_in_valid)
    if(allocated(valid)) deallocate(valid)
  end function getIndicesModelParametersDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief print all data samples and model parameters
!! \param this data model space info
!
  subroutine printDataModelSpaceInfo(this)
    type (data_model_space_info) :: this
    integer :: i
    write(*,*) "#############################################################"
    write(*,*) "DATA-MODEL-SPACE-INFO"
    write(*,*) this%ndata," data samples:"
    write(*,*) "evid    staname   comp    ifreq    imre"
    write(*,*) "-------------------------------------------------------------"
    do i = 1,this%ndata
       write(*,*) "'"//trim(this%evid(i))//"'  ","'"//trim(this%staname(i))//"'  ",&
            "'"//trim(this%comp(i))//"'",this%ifreq(i),"'"//trim(this%imre(i))//"'"
    enddo
    write(*,*) ""
    write(*,*) this%nparam," model parameters of parametrization '"//trim(this%parametrization)//"' :"
    write(*,*) "param   cell "
    write(*,*) "-------------------------------------------------------------"
    do i = 1,this%nparam
       write(*,*) "'"//trim(this%param(i))//"'",this%cell(i)
    end do
    write(*,*) "#############################################################"
  end subroutine printDataModelSpaceInfo
!------------------------------------------------------------------------
end module dataModelSpaceInfo
