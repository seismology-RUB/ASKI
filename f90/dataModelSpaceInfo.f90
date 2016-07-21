!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

!############################################################################
!## IN THE LONG RUN: improve this module !
!## 
!## do not generate the data model space over and over again from scratch
!## 
!## instead: in the beginning of every iteration step, generate one (binary) 
!## file by a separate program which contains ALL information. Assume that
!## all values are checked (do not check for validity). 
!## Additionally: remembber in this file information like paths, npaths,
!## which receiver components, frequencies are contained in which paths, 
!## all different comp, all different ifreq, ...
!## 
!## The best thing would be to memorize all (or a lot of) information required by 
!## the programs like: all data indices per path ....
!## 
!## Alternatively a more convenient access to this information should be
!## implemented in this module.
!## 
!## This could allow for quicker access in the application of the programs
!############################################################################

!> \brief Define and hold information about a set of data and a set of model values
!!
!! \details This module defines a type containing all information about the data samples 
!!  and model values, i.e. data space and model space, which should be used for a 
!!  certain operation like setting up the kernel matrix, doing kernel focussing, etc.
!!  A data sample is defined by a station, a station component, an event, a frequency index
!!  and distinction between imaginary or real part (of the complex valued data).  
!!  A model value is defined by a parametrization (e.g. 'isoVelocity', 
!!  allowing for parameters like 'rho', 'vp', 'vs'), the respective parameter (e.g. 'vs') 
!!  and an inversion grid cell.
!!  Only indices refering to properties defined in a different object (iteration step info, 
!!  inversion grid) are hold in this module, like station/event index, inversion grid cell 
!!  index etc. 
!!  For parallelized application, these data samples and model values may be defined and 
!!  used locally and do not need to represent the whole data and model space.
!!
!! \author Florian Schumacher
!! \date July 2015
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
  use realloc
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateDataModelSpaceInfo; end interface
!!$  interface anyParam; module procedure anyModelParameterDataModelSpaceInfo; end interface
  interface getArrayEvid; module procedure getArrayEvidDataSamplesDataModelSpaceInfo; end interface
  interface getArrayStaname; module procedure getArrayStanameDataSamplesDataModelSpaceInfo; end interface
  interface getIndxDataSamples
     module procedure getIndicesDataSamplesDataModelSpaceInfo
     module procedure getIndicesDataSamplesInOtherDataModelSpaceInfo
  end interface getIndxDataSamples
  interface getIndxModelValues
     module procedure getIndicesModelValuesDataModelSpaceInfo
     module procedure getIndicesModelValuesInOtherDataModelSpaceInfo
  end interface getIndxModelValues
  interface subtractModelValues
     module procedure subtractModelValuesDataModelSpaceInfo
     module procedure subtractModelValuesCopyDataModelSpaceInfo
  end interface subtractModelValues
  interface sameDataSamples; module procedure sameDataSamplesDataModelSpaceInfo; end interface
  interface sameModelValues; module procedure sameModelValuesDataModelSpaceInfo; end interface
  interface getNormalization; module procedure getNormalizationDataSamplesDataModelSpaceInfo; end interface
  interface allComp; module procedure getAllDifferentCompDataSamplesDataModelSpaceInfo; end interface
  interface allIfreq; module procedure getAllDifferentIfreqDataSamplesDataModelSpaceInfo; end interface
  interface allParam; module procedure getAllDifferentParamModelValuesDataModelSpaceInfo; end interface
!!$  interface getArrayComp; module procedure getArrayCompDataSamplesDataModelSpaceInfo; end interface
!!$  interface getArrayIfreq; module procedure getArrayIfreqDataSamplesDataModelSpaceInfo; end interface
!!$  interface getArrayImre; module procedure getArrayImreDataSamplesDataModelSpaceInfo; end interface
!!$  interface getArrayCell; module procedure getArrayCellModelValuesDataModelSpaceInfo; end interface

! FS FS
! instead of having ONE routine which returnes a vector of indices of data samples (model parameters)
! which have certain properties, we can have ONE routine (well two actually: one for data samples, one for model parameters)
! which gets ALL optional parameters staname,evid,comp,ifreq,imre (param,cell respectively)
! and returnes a vector of indices of data samples (model parameters) which have all the requested properties
! (or a pointer to null, if there is no data sample (model parameter) having the requested properties)
! (could be realized like: get arrays of all properties, then use identity array (indices) and "where" statement
! and finally pack statement (using count))
!!$  interface getIndicesParam; module procedure getIndicesParameterDataModelSpaceInfo; end interface
! FS FS
  interface operator (.ndata.); module procedure getNdataDataModelSpaceInfo; end interface
  interface operator (.nmval.); module procedure getNmvalDataModelSpaceInfo; end interface
  interface operator (.npath.); module procedure getNpathDataModelSpaceInfo; end interface
  interface operator (.evid.); module procedure getEvidDataSampleDataModelSpaceInfo; end interface
  interface operator (.staname.); module procedure getStanameDataSampleDataModelSpaceInfo; end interface
  interface operator (.comp.); module procedure getCompDataSampleDataModelSpaceInfo; end interface
  interface operator (.ifreq.); module procedure getIfreqDataSampleDataModelSpaceInfo; end interface
  interface operator (.imre.); module procedure getImreDataSampleDataModelSpaceInfo; end interface
  interface operator (.pmtrz.); module procedure getParametrizationDataModelSpaceInfo; end interface
  interface operator (.param.); module procedure getParamModelValueDataModelSpaceInfo; end interface
  interface operator (.cell.); module procedure getCellModelValueDataModelSpaceInfo; end interface
!
!> \brief meta info defining rows (data space) and columns (model space) of inversion matrix
  type data_model_space_info
     private
     ! data space
     integer :: ndata = 0 !< total number of data in data space (size of evid,staname,comp,ifreq,imre,wdata,normalization*)

     character(len=character_length_evid), dimension(:), pointer :: evid => null() !< event ID
     character(len=character_length_staname), dimension(:), pointer :: staname => null() !< name of station
     character(len=character_length_component), dimension(:), pointer :: comp => null() !< data component (valid components defined by module componentTransformation)
     integer, dimension(:), pointer :: ifreq => null() !< frequency index of all frequencies (either nf1 <= ifreq <= nf2 for equally spaced, or 1 <= ifreq <= nfreq for arbitrary freq.)
     character(len=2), dimension(:), pointer :: imre => null() !< imaginary part (imre='im') or real part (imre='re')
     real, dimension(:), pointer :: wdata => null() !< weight coefficient for each data sample ( 0.0 < wdata <= 1.0 ), will be used to scale equations in kernel linear system
     real, dimension(:), pointer :: normalization_mdata => null() !< normalization factors for measured data (one factor for each datum in this data space)
     real, dimension(:), pointer :: normalization_sdata => null() !< normalization factors for synthetic data (one factor for each datum in this data space)
     real, dimension(:), pointer :: normalization_kernel => null() !< normalization factors for kernels (complete row of kernel matrix; one factor for each datum in this data space; also used for path-specific synthetic corrections, since those are computed as scalar product of kernel matrix row and model vector)

     ! model space
     integer :: nmval = 0 !< total number of model values in model space (size of param,cell)
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
!! \details This routine ASSUMES all incoming values to be valid, there is no check done here!
!!  It is usually called from routines within this module only and not directly called from outside the
!!  module. If you do so, make sure yourself that the incoming values of evid, staname,comp,ifreq,... are valid.
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
  subroutine addDataSamplesDataModelSpaceInfo(this,evids,stanames,comp,ifreq,imre,all_combinations,wdata)
    ! incoming
    type (data_model_space_info) :: this
    character(len=*), dimension(:) :: evids,stanames
    character(len=*), dimension(:) :: comp
    integer, dimension(:) :: ifreq
    character(len=*), dimension(:) :: imre
    real, dimension(:), optional :: wdata !< always assumed to provide one weight per incoming DATA SAMPLE! (also in case all_combinations=.true., in the correct order)
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
    real, dimension(:), pointer :: wdata_tmp
!
    nullify(evid_tmp,staname_tmp,comp_tmp,ifreq_tmp,imre_tmp,wdata_tmp)
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
    ! IT IS ASSUMED THAT WDATA HAS ALSO SIZE NDATA_NEW AND THAT THE ORDER OF THE VALUES IS CONSISTENT WITH THE LOOPS ABOVE
!
    wdata_tmp => this%wdata
    allocate(this%wdata(this%ndata+ndata_new))
    if(this%ndata > 0) this%wdata(1:this%ndata) = wdata_tmp
    if(associated(wdata_tmp)) deallocate(wdata_tmp)
    if(present(wdata)) then
       this%wdata(this%ndata+1:this%ndata+ndata_new) = wdata(:)
    else
       this%wdata(this%ndata+1:this%ndata+ndata_new) = 1.0
    end if
!
    this%ndata = this%ndata + ndata_new
  end subroutine addDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief add model values to data model space info by combination of all input items
!! \details This routine ASSUMES all incoming values to be valid, there is no check done here!
!!  It is usually called from routines within this module only and not directly called from outside the
!!  module. If you do so, make sure yourself that the incoming values of param and cell are valid.
!!  Furthermore it is assumed that all incoming arrays have the same size, except in case of 
!!  all_combinations=.true. when all incoming values are combined to nparam*ncell
!!  model values. So all routines of this module which call this routine have to check validity of
!!  the values and need to assure the correct sizes of the arrays in the respective cases.
!! \param this data_model_space_info object to which model values are added
!! \param param array containing names of parameter (defines this%param)
!! \param cell array of indices of inversion grid cells
!! \param all_combinations optional logical to indicate whether model values for all nparam*ncell
!!  combinations of the incoming arrays are added or not. If not present (or .false.) all incoming arrays are assumed to have same size!
!
  subroutine addModelValuesDataModelSpaceInfo(this,param,cell,all_combinations)
    ! incoming
    type (data_model_space_info) :: this
    character(len=*), dimension(:) :: param
    integer, dimension(:) :: cell
    logical, optional :: all_combinations
    ! local
    integer :: nparam,ncell
    integer :: jparam,jcell
    integer :: nmval_new,iparam
    logical :: add_all_combinations
    character(len=character_length_param), dimension(:), pointer :: param_tmp
    integer, dimension(:), pointer :: cell_tmp
!
    nullify(param_tmp,cell_tmp)
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
       nmval_new = nparam * ncell
    else
       nmval_new = nparam ! it is assumed that all incoming arrays have the same size!!       
    end if
!
    if(nmval_new .le. 0) return
!
    ! reallocate
    param_tmp => this%param
    cell_tmp => this%cell
    allocate(this%param(this%nmval+nmval_new))
    allocate(this%cell(this%nmval+nmval_new))
    if(this%nmval > 0) then
       this%param(1:this%nmval) = param_tmp
       this%cell(1:this%nmval) = cell_tmp
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
             this%param(this%nmval+iparam) = param(jparam)
             this%cell(this%nmval+iparam) = cell(jcell)
          end do ! jcell
       end do ! jparam
!
    else ! add_all_combinations
!
       do iparam = 1,nmval_new
          this%param(this%nmval+iparam) = param(iparam)
          this%cell(this%nmval+iparam) = cell(iparam)
       end do ! iparam
!
    endif ! add_all_combinations
!
    this%nmval = this%nmval + nmval_new
  end subroutine addModelValuesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief subtract model values contained in that from model values contained in this: "this minus that"
!! \param this data and model space from which certain model values will be removed
!! \param that contains model values which will be removed from this (data space ignored)

!
  subroutine subtractModelValuesDataModelSpaceInfo(this,that)
    type (data_model_space_info) :: this,that
    ! local
    integer, dimension(:), allocatable :: indx,indx2
    integer :: nmval_reduced,imval
    character(len=character_length_param), dimension(:), pointer :: param
    integer, dimension(:), pointer :: cell
!
    nullify(param,cell)
!
    if(this%nmval == 0) return
    if(that%nmval == 0) return
    if(this%parametrization /= that%parametrization) return
!
    ! start with all model parameters in this
    allocate(indx(this%nmval))
    indx = (/ (imval,imval=1,this%nmval) /)
!
    ! now detect those model values (indices) which should be removed and mark them by -1 in array indx
    do imval = 1,that%nmval
       where(this%param == that%param(imval) .and. this%cell == that%cell(imval))
          indx = -1
       end where
    end do ! imval
    nmval_reduced = count(indx/=-1)
!
    ! if no model parameters are to be removed, just return
    if(nmval_reduced == this%nmval) goto 1
!
    ! if all model parameters are to be removed, remove model space from object this
    if(nmval_reduced == 0) then
       this%nmval = 0
       deallocate(this%param,this%cell)
       goto 1
    end if
!
    ! if only some model parameters are to be removed: allocate, define and assign new arrays for param and cell
    allocate(indx2(nmval_reduced),param(nmval_reduced),cell(nmval_reduced))
    indx2 = pack(indx,indx/=-1)
    param = this%param(indx2)
    cell = this%cell(indx2)
    deallocate(this%param,this%cell)
    this%param => param
    this%cell => cell
    this%nmval = nmval_reduced
    nullify(param,cell)
!
1   if(allocated(indx)) deallocate(indx)
    if(allocated(indx2)) deallocate(indx2)
  end subroutine subtractModelValuesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief "result_copy = this without that" in terms of model values
!! \param this data and model space from which certain model values will be removed
!! \param that contains model values which will be removed from this (data space ignored)
!! \param result_copy, object will be newly created and contains the result on exit, 'this' is unchanged
!
  subroutine subtractModelValuesCopyDataModelSpaceInfo(this,that,result_copy)
    type (data_model_space_info) :: this,that,result_copy
!
    call copyDataModelSpaceInfo(result_copy,this)
!
    call subtractModelValuesDataModelSpaceInfo(result_copy,that)
!
  end subroutine subtractModelValuesCopyDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief from data model space input file, a data_model_space object is filled
!! \details here, simply the two subroutines createDataSamplesFromFileDataModelSpaceInfo
!!  and createModelValuesFromFileDataModelSpaceInfo are called
!! \param this data_model_space_info object to be created
!! \param eventlist list of all seismic events, to check validity of event IDs
!! \param stationlist list of all seismic stations, to check validity of station names
!! \param ifreq_valid array of valid frequency indices
!! \param parametrization all model values are assumed to be of this parametrization 
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
    call createModelValuesFromFileDataModelSpaceInfo(this,parametrization,ntot_invgrid,intw,filename,lu,errmsg)
    if(.level.errmsg == 2) return
  end subroutine createFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief create data samples from the DATA SAMPLES block of the data model space input file
!! \details This routine searches for the line "DATA SAMPLES" in the given file and reads in the
!!  block following this line. 
!!  Format of 'DATA SAMPLES' block:
!!     line of form: 'WEIGHTING value', where value is one of 'NONE', 'BY_PATH', 'BY_FREQUENCY', 'BY_PATH_AND_FREQUENCY'
!!        in case of 'NONE', all data weights are internally set to 1.0 (i.e. no actual weighting is performed)
!!        values 'BY_PATH' and 'BY_PATH_AND_FREQUENCY' are only allowed in case of 'PATH SPECIFIC', in which case
!!        the pairs 'evid staname' in a path block is expected to be followed by one number > 0.0 and <= 1.0  ;
!!        the frequency dependent weighting values (in cases 'BY_FREQUENCY', 'BY_PATH_AND_FREQUENCY') are defined
!!        in a separate line following the lines of form 'nfreq ifreq_1 ... ifreq_n' (for either case of 'FREQUENCIES ALL' 
!!        or 'FREQUENCIES SPECIFIC'. this separate lines have themselves the form 'nfreq w_1 ... w_n' defining nfreq weights
!!        in range > 0.0 and <= 1.0. 
!!        in case 'BY_PATH_AND_FREQUENCY', both weights (for path and frequency) are MULTIPLIED for each data sample
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
    integer :: n,nev,jev,nstat,jstat,npath,jpath,ncomp,jcomp,nfreq,jfreq,nimre,jimre
    character(len=character_length_evid) :: evid
    character(len=character_length_evid), dimension(:), allocatable :: evids,evids_tmp
    character(len=character_length_staname) :: staname
    character(len=character_length_staname), dimension(:), allocatable :: stanames,stanames_tmp
    character(len=character_length_component), dimension(:), allocatable :: comp,comp_tmp
    integer, dimension(:), allocatable :: ifreq,ifreq_tmp
    real, dimension(:), allocatable :: weight_freq,weight_freq_tmp,weight
    real :: weight_path
    character(len=2), dimension(:), allocatable :: imre,imre_tmp
    logical, dimension(:), allocatable :: valid
    logical :: specific_paths,specific_comp,specific_freq,specific_imre,line_data_samples_found,&
         weight_by_path,weight_by_freq
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
    ! WEIGHTING
    ! 
    read(lu,"(a)",iostat=ios) line; iline = iline+1
    if(ios/=0) then
       close(lu)
       write(errstr,*) "could not read line",iline,", iostat = ",ios,&
            "; expected line containing keyword 'WEIGHTING' followed by one of 'NONE','BY_FREQUENCY','BY_PATH',"//&
            "'BY_PATH_AND_FREQUENCY'"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    read(line,*,iostat=ios) key,val
    if(trim(key) /= 'WEIGHTING') then
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in DATA SAMPLES block not supported. 'WEIGHTING' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    select case(trim(val))
    case('NONE')
       weight_by_path = .false.
       weight_by_freq = .false.
    case('BY_FREQUENCY')
       weight_by_path = .false.
       weight_by_freq = .true.
    case('BY_PATH')
       weight_by_path = .true.
       weight_by_freq = .false.
    case('BY_PATH_AND_FREQUENCY')
       weight_by_path = .true.
       weight_by_freq = .true.
    case default
       write(errstr,*) "'WEIGHTING' specification '"//trim(val)//"' on line ",iline," in DATA SAMPLES block not supported."//&
            " One of 'NONE','BY_FREQUENCY','BY_PATH','BY_PATH_AND_FREQUENCY' expected."
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    end select
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
       if(weight_by_path) then
          close(lu)
          write(errstr,*) "WEIGHTING BY_PATHS or WEIGHTING BY_PATHS_AND_FREQUENCY is not allowed for PATHS ALL, "//&
               "only for PATHS SPECIFIC"
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
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
       if(ios/=0) then
          close(lu)
          write(errstr,*) "could not read ",nfreq," frequency indices starting from second value of line",iline,&
               ", iostat = ",ios
          call add(errmsg,2,trim(errstr),myname)
          return
       end if
       !write(*,*) iline,"#"//trim(line)//"#"
       if(weight_by_freq) then
          read(lu,"(a)",iostat=ios) line; iline = iline+1
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                  "; expected number of frequencies nfreq as first value on that line, "//&
                  "followed by nfreq frequency dependent data weights (> 0 and <= 1.0)"
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          read(line,*,iostat=ios) n
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read number of frequency weights as first value on line",iline,", iostat = ",ios
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          if(n/= nfreq) then
             close(lu)
             write(errstr,*) "number of frequencies ",n," (first value on line",iline,") must be the same as in "//&
                  "line before (",nfreq,")"
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
          allocate(weight_freq_tmp(n))
          read(line,*,iostat=ios) n,weight_freq_tmp
          if(ios/=0) then
             close(lu)
             write(errstr,*) "could not read ",n," frequency dependent data weights starting from second value of line",&
                  iline,", iostat = ",ios
             call add(errmsg,2,trim(errstr),myname)
             return
          end if
       else ! weight_by_freq
          ! for simplicity of the operations below (no if clauses necessary), setup a fake weights array of constant 1.0
          allocate(weight_freq_tmp(nfreq))
          weight_freq_tmp = 1.0
       end if ! weight_by_freq
       ! check if there are invalid frequencies (or invalid weights)
       nvalid = nfreq
       allocate(valid(nvalid))
       do jfreq = 1,nfreq
          valid(jfreq) = any(ifreq_valid == ifreq_tmp(jfreq)) .and. &
               weight_freq_tmp(jfreq) > 0.0 .and. weight_freq_tmp(jfreq) <= 1.0
          if(.not.valid(jfreq)) then
             nvalid = nvalid - 1
             write(errstr,*) jfreq,"'th frequency index '",ifreq_tmp(jfreq),"' on line ",iline, &
                  " or its weight (",weight_freq_tmp(jfreq),") is invalid, hence it is excluded!"
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
          allocate(weight_freq(nvalid))
          weight_freq = pack(weight_freq_tmp,valid)
          nfreq = nvalid
       else
          allocate(ifreq(nfreq))
          ifreq = ifreq_tmp
          allocate(weight_freq(nfreq))
          weight_freq = weight_freq_tmp
       end if
       deallocate(ifreq_tmp,weight_freq_tmp,valid)
       !
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
          if(weight_by_path) then
             read(line,*,iostat=ios) evid,staname,weight_path
             if(ios/=0) then
                write(errstr,*) "could not read 'eventID stationName weight' on line",iline,&
                     ", hence excluding this path!, iostat = ",ios
             else ! only if you could read weight, check if it is valid
                if(weight_path <= 0.0 .or. weight_path > 1.0) then
                   write(errstr,*) "path dependent data weight ",weight_path," (third entry) on line ",iline, &
                        " is invalid (must be > 0.0 and <= 1.0), hence this path is excluded!"
                   ios = -1
                end if
             end if
          else
             read(line,*,iostat=ios) evid,staname
             if(ios/=0) write(errstr,*) "could not read 'eventID stationName' pair on line",iline,&
                  ", hence excluding this path!, iostat = ",ios
             weight_path = 1.0 ! for simplicity of operations below, define a fake weight (no if clauses necessary below)
          end if
          if(ios/=0) then
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
             if(weight_by_freq) then
                read(lu,"(a)",iostat=ios) line; iline = iline+1
                if(ios/=0) then
                   close(lu)
                   write(errstr,*) "could not read line",iline,", iostat = ",ios,&
                        "; expected number of frequencies nfreq as first value on that line, "//&
                        "followed by nfreq frequency dependent data weights (> 0 and <= 1.0)"
                   call add(errmsg,2,trim(errstr),myname)
                   return
                end if
                read(line,*,iostat=ios) n
                if(ios/=0) then
                   close(lu)
                   write(errstr,*) "could not read number of frequency weights as first value on line",iline,", iostat = ",ios
                   call add(errmsg,2,trim(errstr),myname)
                   return
                end if
                if(n/= nfreq) then
                   close(lu)
                   write(errstr,*) "number of frequencies ",n," (first value on line",iline,") must be the same as in "//&
                        "line before (",nfreq,")"
                   call add(errmsg,2,trim(errstr),myname)
                   return
                end if
                allocate(weight_freq_tmp(n))
                read(line,*,iostat=ios) n,weight_freq_tmp
                if(ios/=0) then
                   close(lu)
                   write(errstr,*) "could not read ",n," frequency dependent data weights starting from second value of line",&
                        iline,", iostat = ",ios
                   call add(errmsg,2,trim(errstr),myname)
                   return
                end if
             else ! weight_by_freq
                ! for simplicity of the operations below (no if clauses necessary), setup a fake weights array of constant 1.0
                allocate(weight_freq_tmp(nfreq))
                weight_freq_tmp = 1.0
             end if ! weight_by_freq
             ! check if there are invalid frequencies
             nvalid = nfreq
             allocate(valid(nvalid))
             do jfreq = 1,nfreq
                valid(jfreq) = any(ifreq_valid == ifreq_tmp(jfreq)) .and. &
                     weight_freq_tmp(jfreq) > 0.0 .and. weight_freq_tmp(jfreq) <= 1.0
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
                allocate(weight_freq(nvalid))
                weight_freq = pack(weight_freq_tmp,valid)
                nfreq = nvalid
             else
                allocate(ifreq(nfreq))
                ifreq = ifreq_tmp
                allocate(weight_freq(nfreq))
                weight_freq = weight_freq_tmp
             end if
             deallocate(ifreq_tmp,weight_freq_tmp,valid)
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
          ! now, for this path all information on components, frequencies, imre's and weights are gathered, so finally add data samples
          ! define weight array
          allocate(weight(ncomp*nfreq*nimre))
          n = 0
          do jcomp = 1,ncomp
             do jfreq = 1,nfreq
                do jimre = 1,nimre
                   n = n + 1
                   weight(n) = weight_path*weight_freq(jfreq)
                end do ! jimre
             end do ! jfreq
          end do ! jcomp
          call addDataSamplesDataModelSpaceInfo(this,(/evid/),(/staname/),comp,ifreq,imre,wdata=weight,all_combinations=.true.)
          deallocate(weight)
          if(specific_comp) deallocate(comp)
          if(specific_freq) deallocate(ifreq,weight_freq)
          if(specific_imre) deallocate(imre)
       enddo ! jpath
       !
       ! if arrays comp,ifreq,imre were allocated above (before jpath-loop), deallocate here
       if(.not.specific_comp) deallocate(comp)
       if(.not.specific_freq) deallocate(ifreq,weight_freq)
       if(.not.specific_imre) deallocate(imre)
    else ! specific_paths
       allocate(weight(nev*nstat*ncomp*nfreq*nimre))
       n = 0
       do jev = 1,nev
          do jstat = 1,nstat
             do jcomp = 1,ncomp
                do jfreq = 1,nfreq
                   do jimre = 1,nimre
                      n = n + 1
                      weight(n) = weight_freq(jfreq)
                   end do ! jimre
                end do ! jfreq
             end do ! jcomp
          end do ! jstat
       end do ! jev
       call addDataSamplesDataModelSpaceInfo(this,evids,stanames,comp,ifreq,imre,wdata=weight,all_combinations=.true.)
       deallocate(weight,weight_freq)
       deallocate(evids,stanames,comp,ifreq,imre)
    endif ! specific_paths
!
    close(lu)
    write(errstr,*) this%ndata-ndata_before," data samples were added"
    call add(errmsg,0,trim(errstr),myname)
  end subroutine createDataSamplesFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief create model values from the MODEL VALUES block of the data model space input file
!! \details This routine searches for the line "MODEL VALUES" in the given file and reads in the
!!  block following this line. 
!!  Format of 'MODEL VALUES' block:
!!     line of form: 'INVERSION_GRID_CELLS value', where value is either 'ALL' (all inversion 
!!                   grid cells are taken) or 'SPECIFIC' (specific definition of set of invgrid cells following below)
!!     line of form: 'PARAMETERS value', where value is either 'ALL' (same parameters as used for all inversion grid
!!                   cells) or 'SPECIFIC' ('SPECIFIC' only allowed if 'INVERSION_GRID_CELLS SPECIFIC', specific
!!                   definition of elastic parameters for each inversion grid cell)
!!     If 'PARAMETERS ALL' the next line must be of form 'nparam param_1 ... param_n', 
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
  subroutine createModelValuesFromFileDataModelSpaceInfo(this,parametrization,ntot_invgrid,intw,filename,lu,errmsg)
    type (data_model_space_info) :: this
    character(len=*) :: parametrization
    integer :: ntot_invgrid
    type (integration_weights) :: intw
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=43) :: myname = 'createModelValuesFromFileDataModelSpaceInfo'
    character(len=500) :: line
    integer :: ios,iline,nvalid,nmval_before
    character(len=40) :: key,val
    integer :: ncell,icell,jcell,nparam,jparam
    character(len=character_length_param), dimension(:), allocatable :: param,param_tmp
    logical, dimension(:), allocatable :: valid
    logical :: specific_cells,specific_param,line_model_values_found,raised_warning_invalid_cells
    integer, dimension(:), pointer :: indx
!
    nullify(indx)
!
    call addTrace(errmsg,myname)
    iline = 0
    nmval_before = this%nmval
!
    ! first of all make sure that the current data_mode_space_info object does not already
    ! contain any model values of a parametrization other than requested here
    if(trim(this%parametrization) /= '' .and. trim(this%parametrization) /= trim(parametrization)) then
       call add(errmsg,2,"the model space already contains model values of parametrization '"//&
            trim(this%parametrization)//"', which is different from the parametrization '"//trim(parametrization)//&
            "' of the model values to be added. There are no multiple parametrizations allowed.",myname)
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
    ! parse through whole file until line "MODEL VALUES" and start reading in block from there
    line_model_values_found = .false.
    do while(ios==0)
       read(lu,"(a)",iostat=ios) line; iline = iline+1
       if(ios/=0) exit
       if(line == "MODEL VALUES") then
          line_model_values_found = .true.
          exit
       endif
    enddo
    if(.not.line_model_values_found) then
       write(errstr,*) "could not find line 'MODEL VALUES' in file '"//trim(filename)//"', searching until line ",iline-1
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    write(errstr,*) "start reading 'MODEL VALUES' block from line ",iline," of file '"//trim(filename)//"'"
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
            " in MODEL VALUES block not supported. 'INVERSION_GRID_CELLS' expected."
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
            " in MODEL VALUES block not supported. 'ALL' or 'SPECIFIC' expected."
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
       write(errstr,*) "keyword '"//trim(key)//"' on line ",iline," in MODEL VALUES block not supported. 'PARAMETERS' expected."
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
       allocate(valid(nparam))
       nvalid = nparam
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
                  ", hence no model values are created! Valid parametrizations (parameters) are: ",all_valid_pmtrz_param
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
       write(errstr,*) "'PARAMETERS' specification '"//trim(val)//"' on line ",iline," in MODEL VALUES block not supported."//&
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
             allocate(valid(nparam))
             nvalid = nparam
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
                        ", hence excluding this inversion grid cell!"
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
          call addModelValuesDataModelSpaceInfo(this,param,(/icell/),all_combinations=.true.)
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
       call addModelValuesDataModelSpaceInfo(this,param,indx,all_combinations=.true.)
       deallocate(indx)
       deallocate(param)
    endif ! specific_cells
!
    close(lu)
    write(errstr,*) this%nmval-nmval_before," model values were added"
    call add(errmsg,0,trim(errstr),myname)
  end subroutine createModelValuesFromFileDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return paths contained in this data_model_space_info object
!! \param idata optional index array of data sample indices defining a data substet for which the contained paths should be returned only
!
  function getPathsDataModelSpaceInfo(this,idata_in) result(paths)
    type (data_model_space_info) :: this
    integer, dimension(:), optional :: idata_in
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    ! local
    character(len=character_length_evid), dimension(:), allocatable :: evid
    character(len=character_length_staname), dimension(:), allocatable :: staname
    integer :: idata,npaths,ndata
    logical, dimension(:), allocatable :: mask_paths
!
    nullify(paths)
!
    if(this%ndata == 0) return
!
    ! make local copies of evid and staname to modify
    if(present(idata_in)) then
       if(size(idata_in) <1) return
       allocate(mask_paths(size(idata_in)))
       mask_paths = idata_in >= 1 .and. idata_in <= this%ndata
       ndata = count(mask_paths)
       if(ndata < 1) then
          deallocate(mask_paths)
          return
       end if
       allocate(evid(ndata),staname(ndata))
       evid = this%evid(pack(idata_in,mask_paths))
       staname = this%staname(pack(idata_in,mask_paths))
       deallocate(mask_paths)
    else
       ndata = this%ndata
       allocate(evid(ndata),staname(ndata))
       evid = this%evid
       staname = this%staname
    end if
!
    do idata = 1,ndata
       if(evid(idata)/='' .and. staname(idata)/='') then
          where(evid(idata+1:ndata)==evid(idata) .and. staname(idata+1:ndata)==staname(idata))
             evid(idata+1:ndata) = ''
             staname(idata+1:ndata) = ''
          end where
       end if
    end do ! idata
!
    allocate(mask_paths(ndata))
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
!------------------------------------------------------------------------
!> \brief create normalization factors for this data space
!! \param this data_model_space_info object
!
  subroutine createDataNormalizationDataModelSpaceInfo(this,normalization_type,errmsg,mdata,sdata)
    type (data_model_space_info) :: this
    character(len=*) normalization_type
    type (error_message) :: errmsg
    real, dimension(:), optional :: mdata,sdata
    ! local
    character(len=400) :: errstr
    character(len=47) :: myname = 'createDataNormalizationDataModelSpaceInfo'
    character(len=93) :: valid_normalization_types = &
         "'maxamp_mdata_by_paths', 'maxamp_mdata_by_paths_and_frequency', 'scale_maxamp_mdata_by_paths'"
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    character(len=character_length_component), dimension(:), pointer :: all_comp_path
    character(len=2), dimension(:), pointer :: pimre
    character(len=character_length_staname), dimension(:), pointer :: pstaname
    character(len=character_length_evid), dimension(:), pointer :: pevid
    integer, dimension(:), pointer :: pifreq,idx_data_path,idx_data_comp_jf,all_ifreq_path
    integer, dimension(:,:), allocatable :: idx_data
    real, dimension(:), allocatable :: amp_mdata,amp_sdata
    real :: global_maxamp_mdata,maxamp_mdata,maxamp_sdata
    integer :: i,ipath,jf,nf,npath
    integer :: ipath_global_maxamp_mdata,ifreq_global_maxamp_mdata
    integer, dimension(1) :: maxloc_amp
!
    nullify(paths,all_comp_path,pimre,pstaname,pevid,pifreq,idx_data_path,idx_data_comp_jf,all_ifreq_path)
!
    call addTrace(errmsg,myname)
!
    if(associated(this%normalization_mdata) .or. associated(this%normalization_sdata) .or. &
         associated(this%normalization_kernel)) then
       call add(errmsg,1,"there are already normalization factors defined; deallocating them now before "//&
            "defining new ones",myname)
       if(associated(this%normalization_mdata)) deallocate(this%normalization_mdata)
       if(associated(this%normalization_sdata)) deallocate(this%normalization_sdata)
       if(associated(this%normalization_kernel)) deallocate(this%normalization_kernel)
    end if
!
    ! for each normalization type, different prequisites and conditions must be met. So check here in a select case statement
    select case(normalization_type)
    case ('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency')
       ! check if data space is already defined (cannot define normalization otherwise for these types of normalization)
       if(this%ndata == 0) then
          call add(errmsg,2,"data space is not yet defined; However, this is a prequisite to define normalizations of "//&
               "types 'maxamp_mdata_by_paths' or 'maxamp_mdata_by_paths_and_frequency'",myname)
          return          
       end if
       ! check if BOTH mdata and sdata are present; both are required for these types of normalization
       if(.not.present(mdata) .or. .not.present(sdata)) then
          call add(errmsg,2,"for normalization types 'maxamp_mdata_by_paths' or 'maxamp_mdata_by_paths_and_frequency', "//&
               "BOTH mdata AND sdata vectors must be provided",myname)
          return
       end if
       ! check if the incoming size of mdata and sdata equal the number of data samples
       ! THE VECTORS mdata AND sdata ARE ASSUMED TO BE ORDERED AS DEFINED IN THIS data_model_space_info OBJECT!
       if(size(mdata) /= this%ndata) then
          write(errstr,*) "incoming number of measured data values = ",size(mdata),&
               " does not match the number of data samples = ",this%ndata," contained in this data space object"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(sdata) /= this%ndata) then
          write(errstr,*) "incoming number of synthetic data values = ",size(sdata),&
               " does not match the number of data samples = ",this%ndata," contained in this data space object"
          call add(errmsg,2,errstr,myname)
          return
       end if
    case ('scale_maxamp_mdata_by_paths')
       ! check if data space is already defined (cannot define normalization otherwise for these types of normalization)
       if(this%ndata == 0) then
          call add(errmsg,2,"data space is not yet defined; However, this is a prequisite to define normalizations of "//&
               "type 'scale_maxamp_mdata_by_paths'",myname)
          return          
       end if
       ! check if mdata is present; only mdata required for this types of normalization
       if(.not.present(mdata)) then
          call add(errmsg,2,"for normalization type 'scale_maxamp_mdata_by_paths', mdata vector must be provided",myname)
          return
       end if
       ! check if the incoming size of mdata equals the number of data samples
       ! THE VECTOR mdata IS ASSUMED TO BE ORDERED AS DEFINED IN THIS data_model_space_info OBJECT!
       if(size(mdata) /= this%ndata) then
          write(errstr,*) "incoming number of measured data values = ",size(mdata),&
               " does not match the number of data samples = ",this%ndata," contained in this data space object"
          call add(errmsg,2,errstr,myname)
          return
       end if
    !add your normalization type here and check everything you require
    !case ('your_new_type')
    case default
       call add(errmsg,2,"incoming normalization type '"//trim(normalization_type)//&
            "' not supported; supported types are "//trim(valid_normalization_types),myname)
       return
    end select ! normalization_type
!
    allocate(this%normalization_mdata(this%ndata),this%normalization_sdata(this%ndata),&
         this%normalization_kernel(this%ndata))
!
    ! PREPARE INDEX MAPPINGS ETC., IF NECESSARY
!
!
    ! COMPUTE NORMALIZATIONS
    select case(normalization_type)
    case('maxamp_mdata_by_paths','maxamp_mdata_by_paths_and_frequency','scale_maxamp_mdata_by_paths')
       paths => getPathsDataModelSpaceInfo(this)
       npath = size(paths,2)
       ! loop on all paths
       do ipath = 1,npath
          if(associated(idx_data_path)) deallocate(idx_data_path)
          if(associated(pstaname)) deallocate(pstaname)
          if(associated(pevid)) deallocate(pevid)
          allocate(pevid(1)); pevid(1) = paths(1,ipath)
          allocate(pstaname(1)); pstaname(1) = paths(2,ipath)
          idx_data_path => getIndxDataSamples(this,staname=pstaname,evid=pevid)
!
          all_comp_path => getAllDifferentCompDataSamplesDataModelSpaceInfo(this,idx_data_path)
          if(size(all_comp_path)>1) then
             write(errstr,*) ipath,"'th path (evid ='",trim(paths(1,ipath)),"', staname = '",trim(paths(2,ipath)),&
                  "') has ",size(all_comp_path)," components. So far, only 1 component suppported!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
          ! IDEA FOR THE FUTURE: first loop on frequencies, and then loop on components (append for all components,
          ! i.e. for all frequency test all components)
!
          all_ifreq_path => getAllDifferentIfreqDataSamplesDataModelSpaceInfo(this,idx_data_path)
          nf = size(all_ifreq_path)
          allocate(idx_data(nf,2),amp_mdata(nf))
          do jf = 1,nf
             if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
             if(associated(pifreq)) deallocate(pifreq)
             if(associated(pimre)) deallocate(pimre)
             allocate(pifreq(1)); pifreq(1) = all_ifreq_path(jf)
             idx_data_comp_jf => getIndxDataSamples(this,imre=pimre,ifreq=pifreq,indx_in=idx_data_path)
             if(size(pimre)/= 2 .or. (.not.any(pimre=='im')) .or. (.not.any(pimre=='re'))) then
                write(errstr,*) "inconsistent data set: for component '",trim(all_comp_path(1)),"', freq. indx ",&
                     all_ifreq_path(jf),", path (evid,staname) = ('",trim(paths(1,ipath)),"','",trim(paths(2,ipath)),&
                     "') there are ",size(pimre)," data samples with im/re = ","'"//pimre//&
                     "'.  Must be exactly one 'im', one 're'"
                call add(errmsg,2,errstr,myname)
                goto 1
             end if
             do i = 1,2
                if(pimre(i) == 're') then
                   idx_data(jf,1) = idx_data_comp_jf(i)
                else
                   idx_data(jf,2) = idx_data_comp_jf(i)
                end if
             end do ! i
          end do ! jf
          amp_mdata(:) = sqrt(mdata(idx_data(:,1))*mdata(idx_data(:,1)) + mdata(idx_data(:,2))*mdata(idx_data(:,2)))
!
          maxamp_mdata = maxval(amp_mdata)
          if(ipath==1) then
             global_maxamp_mdata = maxamp_mdata
             ipath_global_maxamp_mdata = 1
             maxloc_amp = maxloc(amp_mdata)
             ifreq_global_maxamp_mdata = all_ifreq_path(maxloc_amp(1))
          else
             if(maxamp_mdata > global_maxamp_mdata) then
                global_maxamp_mdata = maxamp_mdata
                ipath_global_maxamp_mdata = ipath
                maxloc_amp = maxloc(amp_mdata)
                ifreq_global_maxamp_mdata = all_ifreq_path(maxloc_amp(1))
             end if
          end if
!
          select case(normalization_type)
          case('maxamp_mdata_by_paths')
             ! First, here divide through by maxamp of path, multiply by global_maxamp_mdata after loop on paths!
             maxamp_mdata = maxval(amp_mdata)
             this%normalization_mdata(idx_data_path) = 1./maxamp_mdata
             maxloc_amp = maxloc(amp_mdata);i=maxloc_amp(1)
!write(*,*) "path ",ipath,": '",trim(paths(1,ipath)),"','",trim(paths(2,ipath)),"' maxamp_mdata = ",maxamp_mdata,&
!" at ifreq = ",i,"; w = ",1./maxamp_mdata
!
             allocate(amp_sdata(nf))
             amp_sdata(:) = sqrt(sdata(idx_data(:,1))*sdata(idx_data(:,1)) + sdata(idx_data(:,2))*sdata(idx_data(:,2)))
!
             !maxamp_sdata = maxval(amp_sdata) ! THIS YIELDS WRONG RESULTS, DO NOT NORMALIZE AT THE MAXIMUM AMPLITUDE, BUT AT THAT AMPLITUDE WHERE THE DATA IS MAXIMAL (array index i)
             maxamp_sdata = amp_sdata(i)
             this%normalization_sdata(idx_data_path) = 1./maxamp_sdata
             this%normalization_kernel(idx_data_path) = 1./maxamp_sdata ! kernels are also handled as sdata, in this case
             !maxloc_amp = maxloc(amp_sdata);j=maxloc_amp(1)
!!$if(i/=j) then
!!$write(*,*) "path ",ipath,": '",trim(paths(1,ipath)),"','",trim(paths(2,ipath)),"' maxamp_mdata = ",maxamp_mdata,&
!!$" at ifreq = ",i,"; maxamp_sdata = ",maxamp_sdata," at DIFFERENT ifreq = ",j
!!$                stop
!!$end if
!write(*,*) "path ",ipath,": '",trim(paths(1,ipath)),"','",trim(paths(2,ipath)),"' maxamp_sdata = ",maxamp_sdata,&
!"maxval,maxloc ",maxval(amp_sdata),j!" at ifreq = ",j,"; amp_sdata(i_maxamp_mdata)=",amp_sdata(i)!,"; w = ",1./maxamp_sdata
          case('maxamp_mdata_by_paths_and_frequency')
             allocate(amp_sdata(nf))
             amp_sdata(:) = sqrt(sdata(idx_data(:,1))*sdata(idx_data(:,1)) + sdata(idx_data(:,2))*sdata(idx_data(:,2)))
!
             ! First, here divide through by maxamp of path, multiply by global_maxamp_mdata after loop on paths! 
             do jf = 1,nf
                this%normalization_mdata(idx_data(jf,:)) = 1./amp_mdata(jf)
                this%normalization_sdata(idx_data(jf,:)) = 1./amp_sdata(jf)
                this%normalization_kernel(idx_data(jf,:)) = 1./amp_sdata(jf) ! kernels are also handled as sdata, in this case
             end do ! jf
          case('scale_maxamp_mdata_by_paths')
             ! for this type of normalization, all factors (measured, synthetic, kernel) are the SAME for one path: 
             ! simply this path block of the linear system is scaled by the path specific 1/maxamp_mdata 
             this%normalization_mdata(idx_data_path) = 1./maxamp_mdata
             this%normalization_sdata(idx_data_path) = 1./maxamp_mdata
             this%normalization_kernel(idx_data_path) = 1./maxamp_mdata
          end select ! normalization_type
!
          if(allocated(idx_data)) deallocate(idx_data)
          if(allocated(amp_mdata)) deallocate(amp_mdata)
          if(allocated(amp_sdata)) deallocate(amp_sdata)
       end do ! ipath
!
       ! after loop on path, hence finding global_maxamp_mdata, multiply the normalization factors by this number:
       !    TECHNICALLY, DO NOT DO THIS, BUT INSTEAD SCALE THE MAXIMUM DATA TO 1.0 ! (i.e. throw away factor global_maxamp_mdata)
       !    this is done, since for (very) large systems, the absolute numbers might get so large, that single precision makes 
       !    problems in program solveCglsKernelSystem. We get an equivalent kernel system not actually multiplying by 
       !    global_maxamp_mdata, also consistent with  regularization, since regularization scaling types will use maxval of 
       !    kernel matrix for scaling.
       !this%normalization_mdata = this%normalization_mdata * global_maxamp_mdata
       !this%normalization_sdata = this%normalization_sdata * global_maxamp_mdata
       !this%normalization_kernel = this%normalization_kernel * global_maxamp_mdata
!
!!$write(*,*) "global mdata amplitude = ",global_maxamp_mdata," at path '",trim(paths(1,ipath_global_maxamp_mdata)),&
!!$"','",trim(paths(2,ipath_global_maxamp_mdata)),"' and frequency index ",ifreq_global_maxamp_mdata
!
    end select ! normalization_type
!
!do i = 1,size(this%normalization_mdata)
!write(*,*) this%normalization_mdata(i),this%normalization_sdata(i),this%normalization_kernel(i)
!end do
!
    ! clean up
2   if(associated(paths)) deallocate(paths)
    if(associated(all_comp_path)) deallocate(all_comp_path)
    if(associated(all_ifreq_path)) deallocate(all_ifreq_path)
    if(associated(idx_data_path)) deallocate(idx_data_path)
    if(associated(pstaname)) deallocate(pstaname)
    if(associated(pevid)) deallocate(pevid)
    if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
    if(associated(pifreq)) deallocate(pifreq)
    if(associated(pimre)) deallocate(pimre)
!
    ! if routine comes here, everything went alright, so return
    return
!
    ! if an error has occurred after allocating memory for normalization factors, 
    ! destroy everything partially created before returning
1   if(associated(this%normalization_mdata)) deallocate(this%normalization_mdata)
    if(associated(this%normalization_sdata)) deallocate(this%normalization_sdata)
    if(associated(this%normalization_kernel)) deallocate(this%normalization_kernel)
    goto 2
  end subroutine createDataNormalizationDataModelSpaceInfo


!! FS FS
! FUNCTION:
! write data model space as text file (with MODEL VALUES SPECIFIC, DATA SAMPLES SPECIFIC)

! THE TEXT FILE FORMAT IS THE ONLY INTERFACE TO THIS MODULE!
! later on (if needed), people can write routines like "addDataSampleFromValues" (having optional "allCombinations" flag
! which is handled in the same way as in modules addDataSamplesDataModelSpaceInfo,addModelValuesDataModelSpaceInfo)
! where incoming values are checked (in the same way as in routines "add*FromFile") and inside which routines 
! addDataSamplesDataModelSpaceInfo,addModelValuesDataModelSpaceInfo are called. 
!! FS FS


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
    if(associated(this%wdata)) deallocate(this%wdata)
    if(associated(this%normalization_mdata)) deallocate(this%normalization_mdata)
    if(associated(this%normalization_sdata)) deallocate(this%normalization_sdata)
    if(associated(this%normalization_kernel)) deallocate(this%normalization_kernel)
    this%ndata = 0
    if(associated(this%param)) deallocate(this%param)
    if(associated(this%cell)) deallocate(this%cell)
    this%nmval = 0
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
!! \param nmval this%nmval
!! \return number of model values
!
  function getNmvalDataModelSpaceInfo(this) result(nmval)
    type (data_model_space_info), intent(in) :: this
    integer :: nmval
    nmval = this%nmval
  end function getNmvalDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get total number of data paths
!! \param this data model space info
!! \param path number of data paths
!! \return number of data paths
!
  function getNpathDataModelSpaceInfo(this) result(npath)
    type (data_model_space_info), intent(in) :: this
    integer :: npath
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    nullify(paths)
    npath = 0
    paths => getPathsDataModelSpaceInfo(this)
    if(.not.associated(paths)) return
    npath = size(paths,2)
    deallocate(paths)
  end function getNpathDataModelSpaceInfo
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
!! \param idata_in optional index array defining subset of data samples among which to search only
!! \param comp vector of all diffferent components
!! \return all different occurring station components
!
  function getAllDifferentCompDataSamplesDataModelSpaceInfo(this,idata_in) result(comp)
    type (data_model_space_info), intent(in) :: this
    integer, dimension(:), optional :: idata_in
    character(len=character_length_component), dimension(:), pointer :: comp
    ! local
    character(len=character_length_component), dimension(:), allocatable :: comp_tmp
    integer :: idata,ncomp,ndata
    logical, dimension(:), allocatable :: mask
!
    nullify(comp)
    if(this%ndata == 0) return
!
    ! make a local copy of this%comp to modify
    if(present(idata_in)) then
       if(size(idata_in) <1) return
       allocate(mask(size(idata_in)))
       mask = idata_in >= 1 .and. idata_in <= this%ndata
       ndata = count(mask)
       if(ndata < 1) then
          deallocate(mask)
          return
       end if
       allocate(comp_tmp(ndata))
       comp_tmp = this%comp(pack(idata_in,mask))
       deallocate(mask)
    else
       ndata = this%ndata
       allocate(comp_tmp(ndata))
       comp_tmp = this%comp
    end if
!
    do idata = 1,ndata-1
       if(comp_tmp(idata)/='') then
          where(comp_tmp(idata+1:ndata)==comp_tmp(idata))
             comp_tmp(idata+1:ndata) = ''
          end where
       end if
    end do ! idata
!
    allocate(mask(ndata))
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
!> \brief get all different occurring frequencies
!! \param this data model space info
!! \param idata_in optional index array defining subset of data samples among which to search only
!! \param ifreq vector of all different frequency indices
!! \return all different occurring frequencies
!
  function getAllDifferentIfreqDataSamplesDataModelSpaceInfo(this,idata_in) result(ifreq)
    type (data_model_space_info), intent(in) :: this
    integer, dimension(:), optional :: idata_in
    integer, dimension(:), pointer :: ifreq
    ! local
    integer, dimension(:), allocatable :: ifreq_tmp
    integer :: idata,nfreq,ndata
    logical, dimension(:), allocatable :: mask
!
    nullify(ifreq)
    if(this%ndata == 0) return
!
    ! make a local copy of this%ifreq to modify
    if(present(idata_in)) then
       if(size(idata_in) <1) return
       allocate(mask(size(idata_in)))
       mask = idata_in >= 1 .and. idata_in <= this%ndata
       ndata = count(mask)
       if(ndata < 1) then
          deallocate(mask)
          return
       end if
       allocate(ifreq_tmp(ndata))
       ifreq_tmp = this%ifreq(pack(idata_in,mask))
       deallocate(mask)
    else
       ndata = this%ndata
       allocate(ifreq_tmp(ndata))
       ifreq_tmp = this%ifreq
    end if
!
    do idata = 1,ndata-1
       if(ifreq_tmp(idata)/=-1) then
          where(ifreq_tmp(idata+1:ndata)==ifreq_tmp(idata))
             ifreq_tmp(idata+1:ndata) = -1
          end where
       end if
    end do ! idata
!
    allocate(mask(ndata))
    mask = ifreq_tmp/=-1
    nfreq = count(mask)
    if (nfreq==0) goto 1
!
    allocate(ifreq(nfreq))
    ifreq = pack(ifreq_tmp,mask)
!
1   deallocate(mask,ifreq_tmp)
  end function getAllDifferentIfreqDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get all different parameter names
!! \param this data model space info
!! \param imval_in optional index array defining subset of model values among which to search only
!! \param param component this%data_space_info(i)%comp
!! \return all different occurring frequencies
!
  function getAllDifferentParamModelValuesDataModelSpaceInfo(this,imval_in) result(param)
    type (data_model_space_info), intent(in) :: this
    integer, dimension(:), optional :: imval_in
    character(len=character_length_param), dimension(:), pointer :: param
    ! local
    character(len=character_length_param), dimension(:), allocatable :: param_tmp
    integer :: iparam,nparam,nparam_out
    logical, dimension(:), allocatable :: mask
!
    nullify(param)
    if(this%nmval == 0) return
!
    ! make a local copy of this%param to modify
    if(present(imval_in)) then
       if(size(imval_in) <1) return
       allocate(mask(size(imval_in)))
       mask = imval_in >= 1 .and. imval_in <= this%nmval
       nparam = count(mask)
       if(nparam < 1) then
          deallocate(mask)
          return
       end if
       allocate(param_tmp(nparam))
       param_tmp = this%param(pack(imval_in,mask))
       deallocate(mask)
    else
       nparam = this%nmval
       allocate(param_tmp(nparam))
       param_tmp = this%param
    end if
!
    do iparam = 1,nparam-1
       if(param_tmp(iparam)/='') then
          where(param_tmp(iparam+1:nparam)==param_tmp(iparam))
             param_tmp(iparam+1:nparam) = ''
          end where
       end if
    end do ! iparam
!
    allocate(mask(nparam))
    mask = param_tmp/=''
    nparam_out = count(mask)
    if (nparam_out==0) goto 1
!
    allocate(param(nparam_out))
    param = pack(param_tmp,mask)
!
1   deallocate(mask,param_tmp)
  end function getAllDifferentParamModelValuesDataModelSpaceInfo
!!$!------------------------------------------------------------------------
!!$!> \brief get all components as pointer to array
!!$!! \param this data model space info
!!$!! \param comp pointer to components array
!!$!! \return components as pointer to array
!!$!
!!$	function getArrayCompDataSamplesDataModelSpaceInfo(this) result(comp)
!!$	type (data_model_space_info), intent(in) :: this
!!$	character(len=character_length_component), dimension(:), pointer :: comp
!!$	integer :: i
!!$	if(this%ndata==0) then
!!$		comp => null()
!!$	else
!!$		allocate(comp(this%ndata))
!!$		do i = 1,this%ndata
!!$			comp(i) = this%data_space_info(i)%comp
!!$		enddo
!!$	endif
!!$	end function getArrayCompDataSamplesDataModelSpaceInfo
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
!!$!------------------------------------------------------------------------
!!$!> \brief get all frequency indices as pointer to array
!!$!! \param this data model space info
!!$!! \param ifreq pointer to frequency indices
!!$!! \return frequency indices as pointer to array
!!$!
!!$	function getArrayIfreqDataSamplesDataModelSpaceInfo(this) result(ifreq)
!!$	type (data_model_space_info), intent(in) :: this
!!$	integer, dimension(:), pointer :: ifreq
!!$	integer :: i
!!$	if(this%ndata==0) then
!!$		ifreq => null()
!!$	else
!!$		allocate(ifreq(this%ndata))
!!$		do i = 1,this%ndata
!!$			ifreq(i) = this%data_space_info(i)%ifreq
!!$		enddo
!!$	endif
!!$	end function getArrayIfreqDataSamplesDataModelSpaceInfo
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
!!$!------------------------------------------------------------------------
!!$!> \brief get all imre values as pointer to array
!!$!! \param this data model space info
!!$!! \param imre pointer to array of imre values
!!$!! \return imre values as pointer to array
!!$!
!!$	function getArrayImreDataSamplesDataModelSpaceInfo(this) result(imre)
!!$	type (data_model_space_info), intent(in) :: this
!!$	character(len=2), dimension(:), pointer :: imre
!!$	integer :: i
!!$	if(this%ndata==0) then
!!$		imre => null()
!!$	else
!!$		allocate(imre(this%ndata))
!!$		do i = 1,this%ndata
!!$			imre(i) = this%data_space_info(i)%imre
!!$		enddo
!!$	endif
!!$	end function getArrayImreDataSamplesDataModelSpaceInfo
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
!> \brief get parameter of i-th model value
!! \param this data model space info
!! \param i index of model value
!! \param param this%param(i)
!! \return parameter of i-th model parameter
!
  function getParamModelValueDataModelSpaceInfo(this,i) result(param)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    character(len=character_length_param) :: param
    if(i<1 .or. i>this%nmval) then
       param = ''
    else
       param = this%param(i)
    endif
  end function getParamModelValueDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief inversion grid cell index of i-th model parameter
!! \param this data model space info
!! \param i index of model parameter
!! \param cell invgrid cell index this%cell(i)
!! \return inversion grid cell index of i-th model parameter
!
  function getCellModelValueDataModelSpaceInfo(this,i) result(cell)
    type (data_model_space_info), intent(in) :: this
    integer, intent(in) :: i
    integer :: cell
    if(i<1 .or. i>this%nmval) then
       cell = -1
    else
       cell = this%cell(i)
    endif
  end function getCellModelValueDataModelSpaceInfo
!!$!------------------------------------------------------------------------
!!$!> \brief get all invgrid cell indices of model parameters as pointer to array of integer
!!$!! \param this data model space info
!!$!! \param cell pointer to array of cell index values of all model parameters defined in this
!!$!! \return invgrid cell index values of all model parameters as pointer to array
!!$!
!!$	function getArrayCellModelValuesDataModelSpaceInfo(this) result(cell)
!!$	type (data_model_space_info), intent(in) :: this
!!$	integer, dimension(:), pointer :: cell
!!$	integer :: i
!!$	if(this%nmval==0) then
!!$		cell => null()
!!$	else
!!$		allocate(cell(this%nmval))
!!$		do i = 1,this%nmval
!!$			cell(i) = this%model_space_info(i)%cell
!!$		enddo
!!$	endif
!!$	end function getArrayCellModelValuesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief get normalization factors for data samples
!! \param this data model space info
!! \param indx_in array of indices of data samples for which the normalization factors should be returned
!! \param errmsg error message object
!! \param norm_mdata optional pointer to array; if successfull, will be allocated on exit and contains the requested measured data normalization factors
!! \param norm_sdata optional pointer to array; if successfull, will be allocated on exit and contains the requested synthetic data normalization factors
!! \param norm_kernel optional pointer to array; if successfull, will be allocated on exit and contains the requested kernel normalization factors
  subroutine getNormalizationDataSamplesDataModelSpaceInfo(this,indx_in,errmsg,norm_mdata,norm_sdata,norm_kernel)
    type (data_model_space_info) :: this
    integer, dimension(:) :: indx_in
    type (error_message) :: errmsg
    real, dimension(:), pointer, optional :: norm_mdata,norm_sdata,norm_kernel
    ! local
    character(len=400) :: errstr
    character(len=45) :: myname = 'getNormalizationDataSamplesDataModelSpaceInfo'
    integer :: nindx_in
!
    call addTrace(errmsg,myname)
!
    if(.not.(present(norm_mdata).or.present(norm_sdata).or.present(norm_kernel))) then
       call add(errmsg,1,"no output arrays present on input",myname)
       return
    end if
    if(present(norm_mdata)) nullify(norm_mdata)
    if(present(norm_sdata)) nullify(norm_sdata)
    if(present(norm_kernel)) nullify(norm_kernel)
    if(this%ndata == 0) then
       call add(errmsg,1,"no data samples in this data space yet",myname)
       return
    end if
    if(size(indx_in)==0) then
       call add(errmsg,1,"no indizes of data samples given on input",myname)
       return
    else
       nindx_in = size(indx_in)
    end if
!
    if(present(norm_mdata)) then
       if(.not.associated(this%normalization_mdata)) then
          call add(errmsg,2,"normalization factors for mdata are requested, but not yet defined in this object",myname)
          return
       elseif(size(this%normalization_mdata)/=this%ndata) then
          write(errstr,*) "normalization factors for mdata inconsistent: there are ",size(this%normalization_mdata),&
               " factors but ",this%ndata," data samples (must be the same)"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
    if(present(norm_sdata)) then
       if(.not.associated(this%normalization_sdata)) then
          call add(errmsg,2,"normalization factors for sdata are requested, but not yet defined in this object",myname)
          return
       elseif(size(this%normalization_sdata)/=this%ndata) then
          write(errstr,*) "normalization factors for sdata inconsistent: there are ",size(this%normalization_sdata),&
               " factors but ",this%ndata," data samples (must be the same)"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
    if(present(norm_kernel)) then
       if(.not.associated(this%normalization_kernel)) then
          call add(errmsg,2,"normalization factors for kernels are requested, but not yet defined in this object",myname)
          return
       elseif(size(this%normalization_kernel)/=this%ndata) then
          write(errstr,*) "normalization factors for kernel inconsistent: there are ",size(this%normalization_kernel),&
               " factors but ",this%ndata," data samples (must be the same)"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    if(any(indx_in<1 .or. indx_in>this%ndata)) then
       write(errstr,*) "there are incoming requested indices of data samples which are <1 or >ndata = ",this%ndata
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! if code comes here, finally everyhting is checked (I hope), so return the requested values
!
    if(present(norm_mdata)) then
       allocate(norm_mdata(nindx_in))
       norm_mdata = this%normalization_mdata(indx_in)
    end if
    if(present(norm_sdata)) then
       allocate(norm_sdata(nindx_in))
       norm_sdata = this%normalization_sdata(indx_in)
    end if
    if(present(norm_kernel)) then
       allocate(norm_kernel(nindx_in))
       norm_kernel = this%normalization_kernel(indx_in)
    end if
  end subroutine getNormalizationDataSamplesDataModelSpaceInfo
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
!!  satisfy the conditions, it points to null(). For any array of evid,staname,comp,ifreq,imre, 
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
!! \param success optional logical value indicating on exit if there are matching data samples or not. NO OUTPUT IS PRODUCED, INCOMING POINTERS ARE NOT MODIFIED. 
!! \param indx array of indices of data samples which satisfy all above constraints (if present)
!! \return array of indices of data samples which satisfy given properties
!
  function getIndicesDataSamplesDataModelSpaceInfo(this,evid,staname,comp,ifreq,imre,wdata,indx_in,success) result(indx)
    type (data_model_space_info) :: this
    character(len=*), dimension(:), pointer, optional :: evid,staname,comp,imre
    integer, dimension(:), pointer, optional :: ifreq,indx_in
    integer, dimension(:), pointer :: indx
    real, dimension(:), pointer, optional :: wdata
    logical, optional :: success
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
    nullify(evid_tmp,staname_tmp,comp_tmp,ifreq_tmp,imre_tmp)
!
    nullify(indx)
    if(present(wdata)) nullify(wdata)
    if(present(success)) success = .false.
    if(this%ndata == 0) return
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
          if(nindx_search .le. 0) goto 1
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
    if(nindx_found .le. 0) goto 1 ! actually this cannot be, because of lines "if(.not.any(valid)) goto 1" above
    ! if only success status is requested, return now without allocating return values
    if(present(success)) then
       success = .true.
       goto 1
    end if
    allocate(indx(nindx_found))
    if(use_indx_in) then
       indx = pack(indx_in_valid,valid)
    else
       indx = pack( (/ (j,j=1,this%ndata) /),valid)
    end if
!
    ! reallocate property arrays for return, if were present on enter
    ! ALSO REALLOCATE IN CASE OF INCOMING NULLIFIED POINTERS!, this way you can get the values of some property even 
    ! if this property was not used to constrain the index search (this mechanism is always used for wdata)
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
    if(present(wdata)) then
       allocate(wdata(nindx_found))
       wdata = this%wdata(indx)
    end if
!
1   if(allocated(indx_in_valid)) deallocate(indx_in_valid)
    if(allocated(valid)) deallocate(valid)
  end function getIndicesDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return indices of model values which have certain properties
!! \details handling is completely analogous to routine getIndicesDataSamplesDataModelSpaceInfo
!! \param this data model space
!! \param param optional pointer to array of names of model parameters
!! \param cell optional pointer to array of inversion grid cell indices
!! \param indx_in optional pointer to index array of model values among which will be searched only
!! \return array of indices of model values which satisfy the given properties
!
  function getIndicesModelValuesDataModelSpaceInfo(this,param,cell,indx_in) result(indx)
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
    nullify(param_tmp,cell_tmp)
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
    ! define indices among which will be searched for model values having the given properties
    use_indx_in = .false.
    nindx_search = this%nmval
    if(present(indx_in)) then
       if(associated(indx_in)) then
          use_indx_in = .true.
          ! first check if there are indices out of range in incoming array indx_in and choose as indx_in_valid all valid indices
          nindx_in = size(indx_in)
          if(nindx_in == 0) return
          allocate(valid(nindx_in))
          valid = (indx_in>0).and.(indx_in.le.this%nmval)
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
    ! find valid indices to return, i.e. find all model values which satisfy the given conditions (if any)
    allocate(valid(nindx_search))
    valid = .true.
!
    ! find all model values with required parameter name
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
       ! update all valid model values
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if already now there are no model values satisfying the conditions, return
       if(.not.any(valid)) goto 1
    end if
!
    ! find all model values with required inversion grid cell index
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
       ! update all valid model values
       valid = valid .and. valid_tmp
       deallocate(valid_tmp)
       ! if there are no model values satisfying the conditions, return
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
       indx = pack( (/ (j,j=1,this%nmval) /),valid)
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
  end function getIndicesModelValuesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief returns indices of this' model values as contained in that
!! \details For all model values with index j contained in this (1<=j<=this%nmval), 
!!  it is determined if the modal value is also contained in that. If so, the respective 
!!  index position k at which the model value is held in object that is returned 
!!  by array indx as: indx(j) = k. If j is not contained in that, indx(j) = -1.
!!  If model spaces of this or that are empty, or their
!!  parametrizations do not match, or the model spaces are disjoint, then indx is not
!!  associated on exit.
!! \param this contains model space to be located inside that
!! \param that contains model space in which this is located
!! \return indx contains indices of this' model values as contained in that (if so)
!
  function getIndicesModelValuesInOtherDataModelSpaceInfo(this,that) result(indx)
    type (data_model_space_info) :: this,that
    integer, dimension(:), pointer :: indx
    ! local
    integer :: imval
!
    nullify(indx)
!
    if(this%nmval == 0) return
    if(that%nmval == 0) return
    if(this%parametrization /= that%parametrization) return
!
    allocate(indx(this%nmval))
    indx(:) = -1
!
    do imval = 1,that%nmval
       where(this%param == that%param(imval) .and. this%cell == that%cell(imval))
          indx = imval
       end where
    end do ! imval
!
    if(all(indx==-1)) then
       deallocate(indx)
       nullify(indx)
    end if
  end function getIndicesModelValuesInOtherDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief returns indices of this' data samples as contained in that
!! \details For all data samples with index j contained in this (1<=j<=this%ndata), 
!!  it is determined if the data sample is also contained in that. If so, the respective 
!!  index position k at which the data sample is held in object that is returned 
!!  by array indx as: indx(j) = k. If j is not contained in that, indx(j) = -1.
!!  If data spaces of this or that are empty, or the data spaces are disjoint, 
!!  then indx is not associated on exit.
!! \param this contains data space to be located inside that
!! \param that contains data space in which this is located
!! \return indx contains indices of this' data samples as contained in that (if so)
!
  function getIndicesDataSamplesInOtherDataModelSpaceInfo(this,that) result(indx)
    type (data_model_space_info) :: this,that
    integer, dimension(:), pointer :: indx
    ! local
    integer :: idata
!
    nullify(indx)
!
    if(this%ndata == 0) return
    if(that%ndata == 0) return
!
    allocate(indx(this%ndata))
    indx(:) = -1
!
    do idata = 1,that%ndata
       where(this%evid == that%evid(idata) .and. &
            this%staname == that%staname(idata) .and. &
            this%comp == that%comp(idata) .and. &
            this%ifreq == that%ifreq(idata) .and. &
            this%imre == that%imre(idata))
          indx = idata
       end where
    end do ! idata
!
    if(all(indx==-1)) then
       deallocate(indx)
       nullify(indx)
    end if
  end function getIndicesDataSamplesInOtherDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief set model parametrization of model space
!! \details If the model parametrization is already defined and different from the incoming requested one, 
!!  i.e. this%parametrization == pmtrz, or incoming pmtrz is not valid, then an error is returned.
!!  Otherwise this function simply defines this%parametrization = pmtrz.
!! \param this data model space info for which parametrization of model space should be defined
!! \param pmtrz requested parametrization
!! \param errmsg error message
  subroutine setParametrizationDataModelSpaceInfo(this,pmtrz,errmsg)
    type (data_model_space_info) :: this
    character(len=*) :: pmtrz
    type (error_message) :: errmsg
!
    if(this%parametrization /= '') then
       if(this%parametrization /= pmtrz) then
          call add(errmsg,2,"incoming model space is already defined for parametrization '"//trim(this%parametrization)//&
               "' which is different from the requested one '"//trim(pmtrz)//"'",'setParametrizationDataModelSpaceInfo')
          return
       else
          return
       end if
    end if
!
    if(.not.validModelParametrization(pmtrz)) then
       call add(errmsg,2,"incoming requested model parametrization '"//trim(pmtrz)//&
            "' is invalid",'setParametrizationDataModelSpaceInfo')
       return
    end if
!
    this%parametrization = pmtrz
  end subroutine setParametrizationDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief create this as a one-to-one copy of object that
!! \param this data model space to be created as a copy of that
!! \param that data model space to be copied to this
!
  subroutine copyDataModelSpaceInfo(this,that,idata1,idata2,ipath1,ipath2,imval1,imval2)
    type (data_model_space_info) :: this,that
    integer, optional :: idata1,idata2,ipath1,ipath2,imval1,imval2
!
    ! first decide which type of data sample request: by sample index or path index
!
    ! if both types are present do not copy any data samples (not supported in this case), but go to check whether to copy model values 
    if( (present(idata1).or.present(idata2)) .and. (present(ipath1).or.present(ipath2)) ) goto 1
!
    if( (present(idata1).or.present(idata2)) ) then
       call copyDataSamplesByIndicesDataModelSpaceInfo(this,that,idata1,idata2)
    elseif(  (present(ipath1).or.present(ipath2)) ) then
       call copyDataSamplesByPathsDataModelSpaceInfo(this,that,ipath1,ipath2)
    else
       ! copy all data samples
       call copyDataSamplesByIndicesDataModelSpaceInfo(this,that,1,that%ndata)
    end if
!
    ! then copy model values (not yet different types of requests, so no need for decision here)
!
1   call copyModelValuesByIndicesDataModelSpaceInfo(this,that,imval1,imval2)

  end subroutine copyDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine copyDataSamplesByIndicesDataModelSpaceInfo(this,that,idata1,idata2)
    type (data_model_space_info) :: this,that
    integer, optional :: idata1,idata2
    ! local
    integer :: i,idata_start,idata_end,ndata
    integer, dimension(:), allocatable :: indx_data

    if(present(idata1)) then
       if(idata1 < 1 .or. idata1 > that%ndata) return
       idata_start = idata1
    else
       idata_start = 1
    end if
    if(present(idata2)) then
       if(idata2 < idata_start .or. idata2 > that%ndata) return
       idata_end = idata2
    else
       idata_end = that%ndata
    end if
!
    ndata = idata_end -idata_start + 1
    allocate(indx_data(ndata))
    indx_data = (/ (i,i=idata_start,idata_end) /)
!
    call copyDataSamplesByIndexDataModelSpaceInfo(this,that,ndata,indx_data)
    deallocate(indx_data)
  end subroutine copyDataSamplesByIndicesDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine copyDataSamplesByPathsDataModelSpaceInfo(this,that,ipath1,ipath2)
    type (data_model_space_info) :: this,that
    integer, optional :: ipath1,ipath2
    ! local
    integer :: n,ndata,npath,ipath,ipath_start,ipath_end
    character(len=character_length_evid), dimension(:), pointer :: evid
    character(len=character_length_staname), dimension(:), pointer :: staname
    integer, dimension(:), pointer :: indx_data,indx
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
!
    nullify(evid,staname,indx_data,indx,paths)
!
    paths => getPathsDataModelSpaceInfo(that)
    if(.not.associated(paths)) return
    npath = size(paths,2)
!
    if(present(ipath1)) then
       if(ipath1 < 1 .or. ipath1 > npath) return
       ipath_start = ipath1
    else
       ipath_start = 1
    end if
    if(present(ipath2)) then
       if(ipath2 < ipath_start .or. ipath2 > npath) return
       ipath_end = ipath2
    else
       ipath_end = npath
    end if
!
    if(ipath_start == 1 .and. ipath_end == npath) then
       call copyDataSamplesByIndicesDataModelSpaceInfo(this,that,1,that%ndata)
       return
    end if
!
    nullify(indx_data,evid,staname,indx)
    ndata = 0
    ! loop on all paths and collect data sample indices for each path in array indx_data (extending it throughout the loop)
    do ipath = ipath_start,ipath_end
       if(associated(evid)) deallocate(evid)
       if(associated(staname)) deallocate(staname)
       if(associated(indx)) deallocate(indx)
       allocate(evid(1),staname(1))
       evid(1) = paths(1,ipath); staname(1) = paths(2,ipath)
       indx => getIndicesDataSamplesDataModelSpaceInfo(that,evid,staname)
       if(.not.associated(indx)) then
          ! THIS SHOULD ACTUALLY NOT HAPPEN!! OTHERWISE MODULE IS INCONSISTENT
          ! TRY TO INDICATE BY DEALLOCATING RETURNING OBJECT this
          call deallocateDataModelSpaceInfo(this)
          goto 1
       end if
       n = size(indx)
       if(n==0) then
          ! THIS SHOULD ACTUALLY NOT HAPPEN!! OTHERWISE MODULE IS INCONSISTENT
          ! TRY TO INDICATE BY DEALLOCATING RETURNING OBJECT this
          call deallocateDataModelSpaceInfo(this)
          goto 1
       end if
       indx_data => reallocate(indx_data,ndata+n)
       indx_data(ndata+1:ndata+n) = indx
       ndata = ndata + n
    end do ! ipath
!
    call copyDataSamplesByIndexDataModelSpaceInfo(this,that,ndata,indx_data)    
!
    ! clean up
1   if(associated(indx_data)) deallocate(indx_data)
    if(associated(evid)) deallocate(evid)
    if(associated(staname)) deallocate(staname)
    if(associated(indx)) deallocate(indx)
  end subroutine copyDataSamplesByPathsDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine copyDataSamplesByIndexDataModelSpaceInfo(this,that,n,indx_data)
    type (data_model_space_info) :: this,that
    integer :: n
    integer, dimension(n) :: indx_data
!
    ! check whether values in indx_data are in valid range 1,...,that%ndata
    if(n <= 0) return
    if(any(indx_data < 1 .or. indx_data > that%ndata)) return
!
    ! ADD THE INCOMING n DATA SAMPLES TO THOSE PRESENT IN this
!
    this%evid => reallocate(this%evid,this%ndata+n)
    this%evid(this%ndata+1:this%ndata+n) = that%evid(indx_data)

    this%staname => reallocate(this%staname,this%ndata+n)
    this%staname(this%ndata+1:this%ndata+n) = that%staname(indx_data)

    this%comp => reallocate(this%comp,this%ndata+n)
    this%comp(this%ndata+1:this%ndata+n) = that%comp(indx_data)

    this%ifreq => reallocate(this%ifreq,this%ndata+n)
    this%ifreq(this%ndata+1:this%ndata+n) = that%ifreq(indx_data)

    this%imre => reallocate(this%imre,this%ndata+n)
    this%imre(this%ndata+1:this%ndata+n) = that%imre(indx_data)

    this%wdata => reallocate(this%wdata,this%ndata+n)
    this%wdata(this%ndata+1:this%ndata+n) = that%wdata(indx_data)

    ! THE FOLLOWING CAN LEAD TO INCONSISTENCIES, IF YOU ADD ONE DATA SPACE WITH NORMALIZATION, AND ANOTHER WITHOUT!!
    if(associated(that%normalization_mdata)) then
       this%normalization_mdata => reallocate(this%normalization_mdata,this%ndata+n)
       this%normalization_mdata(this%ndata+1:this%ndata+n) = that%normalization_mdata(indx_data)
    end if
    if(associated(that%normalization_sdata)) then
       this%normalization_sdata => reallocate(this%normalization_sdata,this%ndata+n)
       this%normalization_sdata(this%ndata+1:this%ndata+n) = that%normalization_sdata(indx_data)
    end if
    if(associated(that%normalization_kernel)) then
       this%normalization_kernel => reallocate(this%normalization_kernel,this%ndata+n)
       this%normalization_kernel(this%ndata+1:this%ndata+n) = that%normalization_kernel(indx_data)
    end if

    this%ndata = this%ndata + n
  end subroutine copyDataSamplesByIndexDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine copyModelValuesByIndicesDataModelSpaceInfo(this,that,imval1,imval2)
    type (data_model_space_info) :: this,that
    integer, optional :: imval1,imval2
    ! local
    integer :: i,imval_start,imval_end,nmval
    integer, dimension(:), allocatable :: indx_mval

    if(present(imval1)) then
       if(imval1 < 1 .or. imval1 > that%nmval) return
       imval_start = imval1
    else
       imval_start = 1
    end if
    if(present(imval2)) then
       if(imval2 < imval_start .or. imval2 > that%nmval) return
       imval_end = imval2
    else
       imval_end = that%nmval
    end if
!
    nmval = imval_end -imval_start + 1
    allocate(indx_mval(nmval))
    indx_mval = (/ (i,i=imval_start,imval_end) /)
!
    call copyModelValuesByIndexDataModelSpaceInfo(this,that,nmval,indx_mval)
    deallocate(indx_mval)
  end subroutine copyModelValuesByIndicesDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine copyModelValuesByIndexDataModelSpaceInfo(this,that,n,indx_mval)
    type (data_model_space_info) :: this,that
    integer :: n
    integer, dimension(n) :: indx_mval
!
    ! check whether values in indx_mval are in valid range 1,...,that%nmval
    if(n <= 0) return
    if(any(indx_mval < 1 .or. indx_mval > that%nmval)) return
!
    if(this%parametrization /= '') then
       ! in case there are already model parameters defined in this, make sure to append model values
       ! of the same model parametrization, otherwise deallocate returning object and return
       if(that%parametrization /= this%parametrization) then
          call deallocateDataModelSpaceInfo(this)
          return
       end if
    else
       this%parametrization = that%parametrization
    end if
!
    ! ADD THE INCOMING n MODEL VALUES TO THOSE PRESENT IN this
!
    this%param => reallocate(this%param,this%nmval+n)
    this%param(this%nmval+1:this%nmval+n) = that%param(indx_mval)

    this%cell => reallocate(this%cell,this%nmval+n)
    this%cell(this%nmval+1:this%nmval+n) = that%cell(indx_mval)

    this%nmval = this%nmval + n    
  end subroutine copyModelValuesByIndexDataModelSpaceInfo
!------------------------------------------------------------------------
  subroutine mapCellIndicesDataModelSpaceInfo(this,param_name,cell_in,indx_out,dmspace_indx_out)
    ! incoming
    type (data_model_space_info) :: this
    character(len=*) :: param_name
    integer, dimension(:) :: cell_in
    ! returning
    integer, dimension(:), pointer :: indx_out,dmspace_indx_out
    ! local
    integer :: ncell_in,nindx_out,nindx_out_tmp,icell
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: idx_dmspace,pcell
    integer, dimension(:), allocatable :: indx_out_tmp,dmspace_indx_out_tmp
    logical, dimension(:), allocatable :: map
!
    nullify(pparam,idx_dmspace,pcell)
!
    nullify(indx_out,dmspace_indx_out)
!
    if(.not.validParamModelParametrization(this%parametrization,param_name)) return
    ncell_in = size(cell_in)
    if(ncell_in < 1) return
!
    allocate(pparam(1)); pparam(1) = param_name
    allocate(pcell(ncell_in)); pcell = cell_in
    idx_dmspace => getIndicesModelValuesDataModelSpaceInfo(this,param=pparam,cell=pcell)
    if(.not.associated(idx_dmspace)) goto 1
    nindx_out_tmp = size(idx_dmspace)
!
    allocate(indx_out_tmp(nindx_out_tmp),dmspace_indx_out_tmp(nindx_out_tmp),map(nindx_out_tmp))
    indx_out_tmp = -1; dmspace_indx_out_tmp = -1
!
    do icell = 1,ncell_in
       where(pcell == cell_in(icell))
          indx_out_tmp = icell
          dmspace_indx_out_tmp = idx_dmspace
       end where
    end do
    map = indx_out_tmp > 0
    nindx_out = count(map)
!
    if(nindx_out == 0) goto 1
!
    allocate(indx_out(nindx_out),dmspace_indx_out(nindx_out))
    indx_out = pack(indx_out_tmp,map)
    dmspace_indx_out = pack(dmspace_indx_out_tmp,map)
!
1   if(associated(pparam)) deallocate(pparam)
    if(associated(pcell)) deallocate(pcell)
    if(allocated(indx_out_tmp)) deallocate(indx_out_tmp)
    if(allocated(dmspace_indx_out_tmp)) deallocate(dmspace_indx_out_tmp)
    if(allocated(map)) deallocate(map)
  end subroutine mapCellIndicesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return logical indicating whether that contains the VERY SAME data samples (same order) as this
!! \param this first data model space info object
!! \param that second data model space info object
!! \param check_normalization optional logical indicating whether to check equality of all %normalization_* vectors (default: NO check)
!! \param check_weights optional logical indicating whether to check equality of weight of data samples (default: NO check)
!! \param l logical indicating on return whether that contains the VERY SAME data samples (same order) as this
!! \return logical indicating whether that contains the VERY SAME data samples (same order) as this
!
  function sameDataSamplesDataModelSpaceInfo(this,that,check_normalization,check_weights) result(l)
    type (data_model_space_info), intent(in) :: this
    type (data_model_space_info), intent(in) :: that
    logical :: l
    logical, optional :: check_normalization,check_weights
    ! local
    logical :: this_associated,that_associated
!
    ! initiate return variable to false, only indicate true before return when conditions are met
    l = .false.
!
    if(this%ndata /= that%ndata) return
!
    ! technically, this and that contain the very same data space if they both do not contain any data values
    if(this%ndata == 0) goto 1
!
    if(any(this%evid /= that%evid)) return
    if(any(this%staname /= that%staname)) return
    if(any(this%comp /= that%comp)) return
    if(any(this%ifreq /= that%ifreq)) return
    if(any(this%imre /= that%imre)) return
!
    if(present(check_weights)) then
       if(check_weights) then
          if(any(this%wdata /= that%wdata)) return
       end if
    end if
!
    if(present(check_normalization)) then
       if(check_normalization) then
          this_associated = associated(this%normalization_mdata)
          that_associated = associated(that%normalization_mdata)
          if(.not. (this_associated.eqv.that_associated)) return
          if(this_associated) then
             if(any(this%normalization_mdata /= this%normalization_mdata)) return
             if(any(this%normalization_sdata /= this%normalization_sdata)) return
             if(any(this%normalization_kernel /= this%normalization_kernel)) return
          end if
       end if
    end if
!
1   l = .true.
  end function sameDataSamplesDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief return logical indicating whether that contains the VERY SAME model values (same order) as this
!! \param this first data model space info object
!! \param that second data model space info object
!! \param l logical indicating on return whether that contains the VERY SAME model values (same order) as this
!! \return logical indicating whether that contains the VERY SAME model values (same order) as this
!
  function sameModelValuesDataModelSpaceInfo(this,that) result(l)
    type (data_model_space_info), intent(in) :: this
    type (data_model_space_info), intent(in) :: that
    logical :: l
!
    ! initiate return variable to false, only indicate true before return when conditions are met
    l = .false.
!
    if(this%nmval /= that%nmval) return
!
    ! technically, this and that contain the very same model space if they both do not contain any model values
    if(this%nmval == 0) goto 1
!
    if(this%parametrization /= that%parametrization) return

    if(any(this%param /= that%param)) return
    if(any(this%cell /= that%cell)) return
!
1   l = .true.
  end function sameModelValuesDataModelSpaceInfo
!!$!------------------------------------------------------------------------
!!$	subroutine mapInvgridSmoothingIndicesToDataModelSpaceInfo(this,invgrid,nb_idx_in,nb_w_in,center_w_in,&
!!$	  nb_idx_return,nb_w_return,center_w_return)
!!$	! incoming
!!$	type (data_model_space_info) :: this
!!$	type (inversion_grid) :: invgrid
!!$	! return
!!$	real, dimension(:), pointer :: center_w_in,center_w_return
!!$	type (integer_vector_pointer), dimension(:), pointer :: nb_idx_in,nb_idx_return
!!$	type (real_vector_pointer), dimension(:), pointer :: nb_w_in,nb_w_return
!!$	!type (error_message) :: errmsg
!!$	! local
!!$	!real, dimension(:), pointer :: center_w
!!$	!type (integer_vector_pointer), dimension(:), pointer :: nb_idx
!!$	!type (real_vector_pointer), dimension(:), pointer :: nb_w
!!$	integer, dimension(:), pointer :: nb_idx_p
!!$	integer, dimension(:), allocatable :: idx
!!$	!character(len=52) :: myname = 'transformInvgridSmoothingIndicesToDataModelSpaceInfo'
!!$	integer, dimension(:), allocatable :: isoVelocity_rho,isoVelocity_vp,isoVelocity_vs,&
!!$	   isoElasticity_rho,isoElasticity_lambda,isoElasticity_mu
!!$	logical :: any_isoElasticity,any_isoVelocity
!!$	integer :: ntot_invgrid,iparam,icell
!!$!
!!$	!call new(errmsg,myname)
!!$	nullify(center_w_return,nb_idx_return,nb_w_return)
!!$!
!!$	! assure incoming Smoothing constraints exist
!!$	if(.not.(associated(center_w_in).and.associated(nb_idx_in).and.associated(nb_w_in))) return
!!$	! assure they have the correct size
!!$	ntot_invgrid = .ntot.invgrid
!!$	if(.not.(size(center_w_in)==ntot_invgrid .and. size(nb_idx_in)==ntot_invgrid .and. size(nb_w_in)==ntot_invgrid)) return
!!$!
!!$	any_isoElasticity = anyModelPmtrizationDataModelSpaceInfo(this,'isoElasticity')
!!$	any_isoVelocity = anyModelPmtrizationDataModelSpaceInfo(this,'isoVelocity')
!!$!
!!$	if(.not.any_isoElasticity .and. .not.any_isoVelocity) then
!!$		!call add(errmsg,2,'no parameters are of parametrization isoElasticity or isoVelocity',myname)
!!$		return
!!$	endif
!!$	if(any_isoElasticity) then
!!$		allocate(isoElasticity_rho(ntot_invgrid),isoElasticity_lambda(ntot_invgrid),isoElasticity_mu(ntot_invgrid))
!!$		isoElasticity_rho(:) = -1; isoElasticity_lambda(:) = -1; isoElasticity_mu(:) = -1
!!$	endif
!!$	if(any_isoVelocity) then
!!$		allocate(isoVelocity_rho(ntot_invgrid),isoVelocity_vp(ntot_invgrid),isoVelocity_vs(ntot_invgrid))
!!$		isoVelocity_rho(:) = -1; isoVelocity_vp(:) = -1; isoVelocity_vs(:) = -1
!!$	endif
!!$!
!!$	! first check for all invgrid cells, which parameters are in the model space on this cell
!!$	do iparam = 1,this%nmval
!!$		select case(this%model_space_info(iparam)%pmtrization)
!!$		case ('isoElasticity')
!!$			select case(this%model_space_info(iparam)%param)
!!$			case ('rho'); isoElasticity_rho(this%model_space_info(iparam)%cell) = iparam
!!$			case ('lambda'); isoElasticity_lambda(this%model_space_info(iparam)%cell) = iparam
!!$			case ('mu'); isoElasticity_mu(this%model_space_info(iparam)%cell) = iparam
!!$			end select
!!$		case ('isoVelocity')
!!$			select case(this%model_space_info(iparam)%param)
!!$			case ('rho'); isoVelocity_rho(this%model_space_info(iparam)%cell) = iparam
!!$			case ('vp'); isoVelocity_vp(this%model_space_info(iparam)%cell) = iparam
!!$			case ('vs'); isoVelocity_vs(this%model_space_info(iparam)%cell) = iparam
!!$			end select
!!$		end select
!!$	enddo ! iparam
!!$!
!!$	!call getLateralLaplaceSmoothingInversionGrid(invgrid,nb_idx,nb_w,center_w)
!!$	allocate(center_w_return(this%nmval),nb_idx_return(this%nmval),nb_w_return(this%nmval))
!!$	center_w_return(:) = 0.
!!$!
!!$	! for every model parameter, change the invgrid index of its neighbours to their respective model parameters index
!!$	do iparam = 1,this%nmval
!!$		icell = this%model_space_info(iparam)%cell
!!$		nb_idx_p => getVectorPointer(nb_idx_in(icell))
!!$		if(.not.associated(nb_idx_p)) cycle
!!$!
!!$		! for this specific model parameters, select the correct indices of all neighbouring model parameters of the same kind
!!$		! if there are none, leave this smoothing condition empty
!!$		select case(this%model_space_info(iparam)%pmtrization)
!!$		case ('isoElasticity')
!!$			select case(this%model_space_info(iparam)%param)
!!$			case ('rho')
!!$				if(all(isoElasticity_rho(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoElasticity_rho(nb_idx_p)
!!$				endif
!!$			case ('lambda')
!!$				if(all(isoElasticity_lambda(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoElasticity_lambda(nb_idx_p)
!!$				endif
!!$			case ('mu')
!!$				if(all(isoElasticity_mu(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoElasticity_mu(nb_idx_p)
!!$				endif
!!$			end select
!!$		case ('isoVelocity')
!!$			select case(this%model_space_info(iparam)%param)
!!$			case ('rho')
!!$				if(all(isoVelocity_rho(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoVelocity_rho(nb_idx_p)
!!$				endif				
!!$			case ('vp')
!!$				if(all(isoVelocity_vp(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoVelocity_vp(nb_idx_p)
!!$				endif				
!!$			case ('vs')
!!$				if(all(isoVelocity_vs(nb_idx_p) > 0)) then
!!$					allocate(idx(size(nb_idx_p)))
!!$					idx = isoVelocity_vs(nb_idx_p)
!!$				endif				
!!$			end select
!!$		end select
!!$!
!!$		! if for this model parameter, for all inversion grid neighbours a parameter of the same kind (e.g. vp) is present,
!!$		! then define a smoothing condition
!!$		if(allocated(idx)) then
!!$			! set center weight
!!$			center_w_return(iparam) = center_w_in(icell)
!!$			! set indices (of model parameters)
!!$			call allocateVectorPointer(nb_idx_return(iparam),size(idx))
!!$			call fillVectorPointer(nb_idx_return(iparam),idx,1)
!!$			! set smoothing weights
!!$			call allocateVectorPointer(nb_w_return(iparam),size(nb_idx_p))
!!$			call fillVectorPointer(nb_w_return(iparam),getVectorPointer(nb_w_in(icell)),1)
!!$			! deallocate idx
!!$			deallocate(idx)
!!$		endif
!!$	enddo ! iparam
!!$	end subroutine mapInvgridSmoothingIndicesToDataModelSpaceInfo
!-------------------------------------------------------------------------
!> \brief iterate over paths
!! \details return evid and staname of next path. optionally, return array of all data sample indices of that path,
!!  all different components of that path (and all different frequencies of that path?! does it make sense?)
!! \param this data_model_space_info object which is searched
!! \param evid event ID of next path
!! \param staname station name of next path
!! \param indx optional pointer to array of indices of data samples of the next path. will be reallocated here, will be deallocated if next = .false.
!! \param all_comp optional pointer to array of all different component names of the next path. will be reallocated here, will be deallocated if next = .false.
!! \param all_ifreq optional pointer to array of all different frequency indices of the next path. will be reallocated here, will be deallocated if next = .false.
!! \param reset optional logical indicating whether this iterator should be reset (and, hence, everything deallocated etc.) It is sensible to call this function with reset = .true. when leaving a loop before it is finished
!! \param ipath_start optional integer indicating path index at which the loop should start (only recognized in the FIRST call, if invalid or start/end inconsistent then loop terminates)
!! \param ipath_end optional integer indicating path index at which the loop should end (only recognized in the FIRST call, if invalid or start/end inconsistent then loop terminates)
!! \param next logical indicating whether there is a next path or not (if not, the return values are dummies and not meaningful)
!! \return logical indicating whether there is a next path or not
!
  function nextPathDataModelSpaceInfo(this,evid,staname,indx,all_comp,all_ifreq,reset,ipath_start,ipath_end) result(next)
    ! incoming
    type (data_model_space_info) :: this
    ! returning
    character(len=character_length_evid) :: evid
    character(len=character_length_staname) :: staname
    integer, dimension(:), pointer, optional :: indx
    character(len=character_length_component), dimension(:), pointer, optional :: all_comp
    integer, dimension(:), pointer, optional :: all_ifreq
    logical, optional :: reset
    integer, optional :: ipath_start,ipath_end
    logical :: next
    ! save
    integer :: ipath = 0
    integer :: ipath_start_internal,ipath_end_internal
    character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
    save :: ipath,paths,ipath_start_internal,ipath_end_internal
    ! local
    integer, dimension(:), pointer :: indx_tmp
    character(len=character_length_evid), dimension(:), pointer :: evid_tmp
    character(len=character_length_staname), dimension(:), pointer :: staname_tmp
!
    nullify(indx_tmp,evid_tmp,staname_tmp)
!
    ! if this iterator is to be reset, do so
    if(present(reset)) then
       if(reset) goto 1
    end if
!
    ! in the first iteration, get paths and setup iteration range
    ! ONLY HERE THE OPTIONAL INDICES ipath_start,ipath_end ARE ACCOUNTED FOR
    if(ipath == 0) then
       paths => getPathsDataModelSpaceInfo(this)
       if(.not.associated(paths)) goto 1
       ! define start index of path iteration (if not present, start at path 1)
       if(present(ipath_start)) then
          if(ipath_start < 1) goto 1
          ipath_start_internal = ipath_start
       else
          ipath_start_internal = 1
       end if
       ! define end index of path iteration (if not present, end at last path)
       if(present(ipath_end)) then
          if(ipath_end > size(paths,2)) goto 1
          ipath_end_internal = ipath_end
       else
          ipath_end_internal = size(paths,2)
       end if
       ! in case of incoming start/end indices, they could be inconsistent, so check
       if(ipath_start_internal > ipath_end_internal) goto 1
       ! initiate the counter according to ipath_start_internal (will be incremented below)
       ipath = ipath_start_internal - 1
    end if
!
    ! increment counter
    ipath = ipath+1
!
    ! if last path is exeeded, reset the iterator
    if(ipath > ipath_end_internal) goto 1
!
    evid = paths(1,ipath)
    staname = paths(2,ipath)
!
    ! if any of the optional pointers is present, the indices of data samples of this path are required
    if(present(indx) .or. present(all_comp) .or. present(all_ifreq)) then
       allocate(evid_tmp(1)); evid_tmp(1) = evid
       allocate(staname_tmp(1)); staname_tmp(1) = staname       
       indx_tmp => getIndicesDataSamplesDataModelSpaceInfo(this,evid_tmp,staname_tmp)
       ! indx_tmp should be associated here, otherwise this module is corrupt
       if(associated(evid_tmp)) deallocate(evid_tmp)
       if(associated(staname_tmp)) deallocate(staname_tmp)
    end if
!
    ! first treat all_comp and all_ifreq, using indices in indx_tmp (before deallocating/nullifying it)
    if(present(all_comp)) then
       if(associated(all_comp)) deallocate(all_comp)
       all_comp => getAllDifferentCompDataSamplesDataModelSpaceInfo(this,indx_tmp)
    end if
    if(present(all_ifreq)) then
       if(associated(all_ifreq)) deallocate(all_ifreq)
       all_ifreq => getAllDifferentIfreqDataSamplesDataModelSpaceInfo(this,indx_tmp)
    end if
!
    ! if the indices are requested as output, point to indx_tmp and nullify indx_tmp, otherwise deallocate indx_tmp
    if(present(indx)) then
       if(associated(indx)) deallocate(indx)
       indx => indx_tmp
       nullify(indx_tmp)
    else
       if(associated(indx_tmp)) deallocate(indx_tmp); nullify(indx_tmp)
    end if
!
    ! IF FUNCTION COMES HERE, RETURN NORMALLY
    next = .true.
    return
!
    ! IF goto 1 RESET THE ITERATOR BEFORE RETURNING
1   ipath = 0
    if(associated(paths)) deallocate(paths); nullify(paths)
    evid = ''; staname = ''
    if(present(indx)) then
       if(associated(indx)) deallocate(indx); nullify(indx)
    end if
    if(present(all_comp)) then
       if(associated(all_comp)) deallocate(all_comp); nullify(all_comp)
    end if
    if(present(all_ifreq)) then
       if(associated(all_ifreq)) deallocate(all_ifreq); nullify(all_ifreq)
    end if
    next = .false.
    return
  end function nextPathDataModelSpaceInfo
!------------------------------------------------------------------------
!> \brief print all data samples and model values
!! \param this data model space info
!
  subroutine printDataModelSpaceInfo(this)
    type (data_model_space_info) :: this
    integer :: i
    write(*,*) "#############################################################"
    write(*,*) "DATA-MODEL-SPACE-INFO"
    write(*,*) this%ndata," data samples:"
    write(*,*) "evid    staname   comp    ifreq    imre    weight"
    write(*,*) "-------------------------------------------------------------"
    do i = 1,this%ndata
       write(*,*) "'"//trim(this%evid(i))//"'  ","'"//trim(this%staname(i))//"'  ",&
            "'"//trim(this%comp(i))//"'",this%ifreq(i),"'"//trim(this%imre(i))//"'",this%wdata(i)
    enddo
    write(*,*) ""
    write(*,*) this%nmval," model values of parametrization '"//trim(this%parametrization)//"' :"
    write(*,*) "param   cell "
    write(*,*) "-------------------------------------------------------------"
    do i = 1,this%nmval
       write(*,*) "'"//trim(this%param(i))//"'",this%cell(i)
    end do
    write(*,*) "#############################################################"
  end subroutine printDataModelSpaceInfo
!------------------------------------------------------------------------
end module dataModelSpaceInfo
