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
!> \brief program for transforming time-domain measured data to ASKI-conform frequency-domain measured data files
!!
!! \details At the moment, two forms of input data are supported: plain text trace files and *.su files.
!!  For both cases, a parameter file has to be provided (filename of parameter file as argument of options
!!  -txt or -su, respectively) which contains some specifications about the location of the files, their naming
!!  and their content. See templates of the respective parameter file for further documentation on its content.
!!
!! \author Florian Schumacher
!! \date January 2015

program createMeasuredData

  use argumentParser
  use string
  use inversionBasics
  use dataModelSpaceInfo
  use inputParameter
  use asciiDataIO
  use dataSu
  use seismicEventList
  use seismicEvent
  use seismicNetwork
  use seismicStation
  use discreteFourierTransform
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: argparse
  character(len=max_length_string) :: main_parfile,filename
  logical :: input_data_is_txt,input_data_is_su

  type (file_unit_handler) :: fuh

  type (inversion_basics) :: invbasics

  type (data_model_space_info) :: dmspace
  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname
  character(len=character_length_component), dimension(:), pointer :: comp_path

  type (seismic_event) :: event
  type (seismic_station) :: station
  character(len=character_length_component) :: comp

  character(len=max_length_string) :: path_su_files

  type (input_parameter) :: inpar
  character (len=80), dimension(5) :: inpar_keys_txt
  data inpar_keys_txt/'DT', 'NSTEP', 'COLUMN_OF_TRACE', 'FILE_DATA_MODEL_SPACE_INFO',&
       'PATH_TXT_TRACES'/
  real :: dt
  integer :: nstep,column_of_trace
  character(len=max_length_string) :: parfile_txt,path_txt_traces,dmspace_file

  type (data_su) :: su_file
  integer :: su_sampint

  logical :: apply_hanning_taper
  real :: portion_hanning_taper

  type (discrete_fourier_transform) :: DFT
  real :: df
  integer :: nfreq
  integer, dimension(:), pointer :: jf
  real, dimension(:,:), pointer :: traces
  real, dimension(:), pointer :: trace
  complex, dimension(:), allocatable :: spectrum

  type (error_message) :: errmsg
  character(len=18) :: prog_name = 'createMeasuredData'

  logical :: terminate_program,next,file_exists,initiated_DFT
  integer :: ios,icomp,ncomp,istat,nstat,nev

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  terminate_program = .false.

!------------------------------------------------------------------------
!  processing arguments and options

  call init(argparse,prog_name,'Fourier transform of time-domain data to ASKI-conform frequency-domain measured data; '//&
       'frequency discretization of measured data as defined in main parfile')

  ! define optional arguments
  call addOption(argparse,'-txt',.true.,&
       'input data files are plain text: parameter file with details on location, naming and content of txt '//&
       'data traces (exactly one of -txt, -su must be set)','sval','parfile_data_txt')
  call addOption(argparse,'-su',.true.,&
       "input data files are seismic unix: argument is path where the input files named eventID_COMP.su can be found "//&
       "(exactly one of -txt, -su must be set)",&
       'sval','./data_su/')
  call addOption(argparse,'-htaper',.true.,&
       "apply cosign-hanning taper to time-domain traces before Fourier Transform. Argument gives portion of "//&
       "the end of the time-series (between 0.0 and 1.0) to which taper is applied",'rval','0.05')
  !call addPosarg(argparse,'dmspace_file','sval',&
  !     'data-model-space info file which defines the data paths and data components')
  call addPosarg(argparse,'main_parfile','sval',&
       'Main parameter file of ASKI inversion; defines PATH_MEASURED_DATA and frequency discretization')
  call parse(argparse)
  if (.level.(.errmsg.argparse) == 2) then
     call print(.errmsg.argparse)
     call usage(argparse)
     goto 1
  end if

  ! get values of positional arguments
  !dmspace_file = argparse.sval.'dmspace_file'
  main_parfile = argparse.sval.'main_parfile'
  if (.level.(.errmsg.argparse) == 2) then
     call usage(argparse)
     goto 1
  end if

  ! check options
  input_data_is_txt = argparse.optset.'-txt'
  input_data_is_su = argparse.optset.'-su'
  apply_hanning_taper = argparse.optset.'-htaper'
  if (.level.(.errmsg.argparse) == 2) then
     call usage(argparse)
     goto 1
  end if

  ! EXACTLY ONE of the options -txt, -su MUST BE SET. IF THAT IS NOT THE CASE, RAISE AN ERROR
  if(input_data_is_txt .eqv. input_data_is_su) then
     write(*,*) "ERROR: EXACTLY ONE of the options -txt, -su must be set"
     call usage(argparse)
     goto 1
  end if

  ! get values of optional arguments
  if(input_data_is_txt) parfile_txt = argparse.sval.'-txt'
  if(input_data_is_su) path_su_files = argparse.sval.'-su'
  if(apply_hanning_taper) portion_hanning_taper = argparse.rval.'-htaper'
  if (.level.(.errmsg.argparse) == 2) goto 1

  ! print the values used
  call document(argparse)
  write(*,*) ""

!------------------------------------------------------------------------
!  other preliminary stuff

  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)

!------------------------------------------------------------------------
!  setup basics

  write(*,*) "initiating inversion basics"
  call new(errmsg,prog_name)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,prog_name)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "done"
  write(*,*) ""

!------------------------------------------------------------------------
!  do the transformation: fork to subroutines, dependent on the kind of input data
!  assume that EXACTLY ONE of the possibilities is logically .true.

  if(input_data_is_txt) then
     call transformMeasuredDataTxt()
  elseif(input_data_is_su) then
     call transformMeasuredDataSu()
  end if
  ! this is actually redundant for now; keept it here in case of adding more code below this line in the future
  if(terminate_program) goto 1 

!------------------------------------------------------------------------
!  clean up

1 call dealloc(invbasics)
  call dealloc(errmsg)
  call dealloc(argparse)
  call dealloc(fuh)

  write(*,*) "good bye"
  write(*,*) ""

  ! program is finished here, so stop (stop command probably not necessary)
  stop

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains

  subroutine transformMeasuredDataTxt()
    write(*,*) "transforming txt traces to ASKI spectra now"
    write(*,*) ""

    ! read in parameter file
    call createKeywordsInputParameter(inpar,inpar_keys_txt)
    write(*,*) "reading parameter file '",trim(parfile_txt),"'"
    call new(errmsg,prog_name)
    call readSubroutineInputParameter(inpar,get(fuh),parfile_txt,errmsg)
    call undo(fuh)
    call addTrace(errmsg,prog_name)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)

    path_txt_traces = inpar.sval.'PATH_TXT_TRACES'
    write(*,*) "   PATH_TXT_TRACES = ",trim(path_txt_traces)
    dmspace_file = inpar.sval.'FILE_DATA_MODEL_SPACE_INFO'
    write(*,*) "   FILE_DATA_MODEL_SPACE_INFO = ",trim(dmspace_file)

    dt = rval(inpar,'DT',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: could not read real value for key 'DT' from its value string '"//&
            trim(inpar.sval.'DT')//"' in parameter file '"//trim(parfile_txt)//"'"
       call dealloc(inpar)
       goto 2
    end if
    write(*,*) "   DT = ",dt

    nstep = ival(inpar,'NSTEP',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: could not read integer value for key 'NSTEP' from its value string '"//&
            trim(inpar.sval.'NSTEP')//"' in  parameter file '"//trim(parfile_txt)//"'"
       call dealloc(inpar)
       goto 2
    end if
    write(*,*) "   NSTEP = ",nstep

    column_of_trace = ival(inpar,'COLUMN_OF_TRACE',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: could not read integer value for key 'COLUMN_OF_TRACE' from its value string '"//&
            trim(inpar.sval.'COLUMN_OF_TRACE')//"' in  parameter file '"//trim(parfile_txt)//"'"
       call dealloc(inpar)
       goto 2
    end if
    write(*,*) "   COLUMN_OF_TRACE = ",column_of_trace

    write(*,*) "done"
    write(*,*) ""

    call dealloc(inpar)

    !  get the data paths and components for which measured 
    write(*,*) "creating data space from file '",trim(dmspace_file),"'"
    call new(errmsg,prog_name)
    call createDataSamplesFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
         .ifreq.invbasics,trim(dmspace_file),get(fuh),errmsg)
    call undo(fuh)
    call addTrace(errmsg,prog_name)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)
    write(*,*) "   there are ",.ndata.dmspace," data samples from ",.npath.dmspace," data paths"
    write(*,*) "done"
    write(*,*) ""

    ! initiate DFT object
    df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
    jf => .ifreq.invbasics
    nfreq = ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')
    write(*,*) "initiating discrete Fourier transform operations"
    write(*,*) "   DF = ",df
    write(*,*) "   number of frequencies = ",nfreq
    write(*,*) "   frequencies [Hz] = ",df*jf
    call new(errmsg,prog_name)
    if(apply_hanning_taper) then
       write(*,*) "   applying cosign-hanning taper to the last ",portion_hanning_taper*100.0,&
            " percent of the time-series"
       call initiateForwardDFT(DFT,dt,0,nstep-1,jf*df,errmsg,hanning_taper=portion_hanning_taper)
    else
       write(*,*) "   do not apply any taper to the time-series"
       call initiateForwardDFT(DFT,dt,0,nstep-1,jf*df,errmsg)
    end if
    call addTrace(errmsg,prog_name)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
    write(*,*) "done"
    write(*,*) ""

    allocate(spectrum(nfreq))

    ! iterate over all paths
    write(*,*) "iterating over all paths now, writing output files '",&
         trim((.inpar.invbasics).sval.'PATH_MEASURED_DATA'),"data_EVID_STANAME_COMP'"
    write(*,*) ""
    do while(nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path))
       ! if there is a next path, then comp must be associated here, otherwise module dataModelSpaceInfo is corrupt, so trust the output here and do not check again if comp is associated
       ncomp = size(comp_path)
       write(*,*) trim("   path '" + evid + "', '" + staname + "';")," present components: ","'"//comp_path//"', "
       do icomp = 1,ncomp
          ! read in this trace
          filename = path_txt_traces + evid + "_" + staname + "_" + comp_path(icomp) + ".txt"
          if(associated(traces)) deallocate(traces)
          errmsg = readAsciiData(filename,get(fuh),traces,column_of_trace,nrow=nstep)
          call undo(fuh)
          call addTrace(errmsg,prog_name)
          if (.level.errmsg /= 0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)

          ! Fourier transform of this trace
          call new(errmsg,prog_name)
          call transformForwardDFT(DFT,traces(:,column_of_trace),spectrum,errmsg)
          call addTrace(errmsg,prog_name)
          if (.level.errmsg /= 0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)

          ! write this spectrum file
          filename = ((.inpar.invbasics).sval.'PATH_MEASURED_DATA') + &
               "data_" + evid + "_" + staname + "_" + comp_path(icomp)
          errmsg = writeAsciiData(filename,get(fuh),spectrum)
          call undo(fuh)
          call addTrace(errmsg,prog_name)
          if (.level.errmsg /= 0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
       end do ! icomp
    end do ! while nextPath

    write(*,*) ""
    write(*,*) "loop over all paths done, transformation of txt traces to ASKI spectra completed"
    write(*,*) ""

    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   if(associated(traces)) deallocate(traces)
    if(allocated(spectrum)) deallocate(spectrum)
    call dealloc(DFT)
    call dealloc(dmspace)
    return

    ! due to an error, terminate program after return
2   terminate_program = .true.
    ! reset the iterator in order to account for those cases where the above loop was exited before it 
    ! finished normally
    next = nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path,reset=.true.)
    goto 1
  end subroutine transformMeasuredDataTxt

!------------------------------------------------------------------------

  subroutine transformMeasuredDataSu()
    nstat = .nstat.(.statlist.invbasics)
    nev = .nev.(.evlist.invbasics)
    write(*,*) "transforming Seismic Unix filesto ASKI spectra now"
    write(*,*) "   assuming one file per receiver component in path '",trim(path_su_files),"' named 'EVID_COMP.su'"
    write(*,*) "   for all ",nev," shots, each containing ",nstat," traces for all receivers "
    write(*,*) "   ordered as ASKI receiver list, all having the same time discretization"
    write(*,*) ""

    initiated_DFT = .false.

    write(*,*) "iterating over all events and possible components now, writing output files '",&
         trim((.inpar.invbasics).sval.'PATH_MEASURED_DATA'),"data_EVID_STANAME_COMP'"
    write(*,*) ""

    do while(nextEventSeismicEventList(.evlist.invbasics,event))
       evid = .evid.event
       do while(nextComponent(comp))

          filename = path_su_files + evid + "_" + comp + ".su"
          inquire(file=filename,exist=file_exists)

          if(file_exists) then
             write(*,*) "   reading and processing su file '",trim(filename),"' now"
             call readDataSu(su_file,get(fuh),filename)
             call undo(fuh)
             if(.n.su_file /= nstat) then
                write(*,*) "ERROR: su file contains ",.n.su_file," data traces, but expecting ",nstat
                goto 2
             end if

             if(.not.initiated_DFT) then
                su_sampint = su_file.sampint.1
                nstep = su_file.numsamp.1
                ! initiate DFT object
                df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
                jf => .ifreq.invbasics
                nfreq = ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')
                write(*,*) "   initiating discrete Fourier transform operations for"
                write(*,*) "      DT = ",su_sampint * 1.e-6
                write(*,*) "      NSTEP = ",nstep
                write(*,*) "      DF = ",df
                write(*,*) "      number of frequencies = ",nfreq
                write(*,*) "      frequencies [Hz] = ",df*jf
                call new(errmsg,prog_name)
                if(apply_hanning_taper) then
                   write(*,*) "      applying cosign-hanning taper to the last ",portion_hanning_taper*100.0,&
                        " percent of the time-series"
                   call initiateForwardDFT(DFT,su_sampint * 1.e-6,0,nstep-1,jf*df,errmsg,&
                        hanning_taper=portion_hanning_taper)
                else
                   write(*,*) "      do not apply any taper to the time-series"
                   call initiateForwardDFT(DFT,su_sampint * 1.e-6,0,nstep-1,jf*df,errmsg)
                end if
                call addTrace(errmsg,prog_name)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) goto 2
                call dealloc(errmsg)
                allocate(spectrum(nfreq))
                initiated_DFT = .true.
             end if ! .not.initiated_DFT

             ! now loop on all receivers
             istat = 0
             do while(nextStationSeismicNetwork(.statlist.invbasics,station))
                staname = .staname.station
                istat = istat+1

                ! check if current su trace has the same time discretization as the previous traces
                if(su_sampint /= (su_file.sampint.istat)) then
                   write(*,*) "ERROR: ",istat,"'th trace ('",trim(staname),"') in su file has sampint = ",&
                        su_file.sampint.istat,", but expecting the same number for all traces, namely ",su_sampint
                end if
                if(nstep /= (su_file.numsamp.istat)) then
                   write(*,*) "ERROR: ",istat,"'th trace ('",trim(staname),"') in su file has numsamp = ",&
                        su_file.numsamp.istat,", but expecting the same number for all traces, namely ",nstep
                end if

                trace => su_file.traceieee.istat

                ! Fourier transform of this trace
                call new(errmsg,prog_name)
                call transformForwardDFT(DFT,trace,spectrum,errmsg)
                call addTrace(errmsg,prog_name)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) goto 2
                call dealloc(errmsg)

                ! write this spectrum file
                filename = ((.inpar.invbasics).sval.'PATH_MEASURED_DATA') + &
                     "data_" + evid + "_" + staname + "_" + comp
                errmsg = writeAsciiData(filename,get(fuh),spectrum)
                call undo(fuh)
                call addTrace(errmsg,prog_name)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) goto 2
                call dealloc(errmsg)
             end do ! nextStation

          else ! file_exists

             write(*,*) "   su file '",trim(filename),"' for event '",trim(evid),"' and component '",&
                  trim(comp),"' does not exist, this file is skipped"
             cycle
          end if ! file_exists

       end do ! nextComponent
    end do ! nextEvent

    write(*,*) ""
    write(*,*) "loop over all events, components and receivers done, transformation of all present su files ",&
         "to ASKI spectra completed"
    write(*,*) ""

    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   if(allocated(spectrum)) deallocate(spectrum)
    call dealloc(DFT)
    return

    ! due to an error, terminate program after return
2   terminate_program = .true.
    ! reset the iterators in order to account for those cases where the above loops were exited before they finished normally
    next = nextStationSeismicNetwork(.statlist.invbasics,station,reset=.true.)
    next = nextComponent(comp,reset=.true.)
    next = nextEventSeismicEventList(.evlist.invbasics,event,reset=.true.)
    goto 1
  end subroutine transformMeasuredDataSu

end program createMeasuredData
