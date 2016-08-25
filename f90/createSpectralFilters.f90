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
program createSpectralFilters
  use inversionBasics
  use seismicEventList
  use seismicEvent
  use discreteFourierTransform
  use complexKernelFrequency
  use inputParameter
  use asciiDataIO
  use argumentParser
  use string
  use mathConstants
  use fileUnitHandler
  use errorMessage

  implicit none

  ! command line
  type (argument_parser) :: ap
  character(len=max_length_string) :: str,main_parfile,eventf_parfile!,stationf_parfile

  ! basics
  type (error_message) :: errmsg
  character(len=21) :: prog_name = 'createSpectralFilters'
  type (file_unit_handler) :: fuh
  type (inversion_basics) :: invbasics
  
  ! event filter
  type (input_parameter) :: eventf_inpar
  character (len=80), dimension(23) :: eventf_inpar_keys
  data eventf_inpar_keys/'CREATE_FILTERS_BY_SOURCE_WAVELET', 'CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH', &
       'EVENT_IDS', 'NUMBER_OF_EVENTS', 'CONVOLVE_WITH_OTHER_WAVELET', 'MULTIPLY_BY_OTHER_SPECTRAL_FILTER', 'STF_FILE', &
       'STF_COLUMN_OF_TRACE', 'STF_DT', 'STF_NSTEP', 'BTW_LOW_PASS_APPLY', 'BTW_LOW_PASS_ORDER', &
       'BTW_LOW_PASS_FC', 'BTW_HIGH_PASS_APPLY', 'BTW_HIGH_PASS_ORDER', 'BTW_HIGH_PASS_FC', 'CONV_WAVEL_FILE', &
       'CONV_WAVEL_COLUMN_OF_TRACE', 'CONV_WAVEL_DT', 'CONV_WAVEL_NSTEP', 'DECONVOLVE_WAVEL_INSTEAD', &
       'CONV_SPEC_FILE', 'DECONVOLVE_SPEC_INSTEAD'/
  logical :: create_filters_by_source_wavelet, create_filters_by_spectral_butterworth,convolve_with_other_wavelet,&
       multiply_by_other_spectral_filter,btw_low_pass_apply,btw_high_pass_apply,deconvolve_wavelet_instead,&
       divide_by_spectral_filter_instead,use_all_event_ids
  character(len=500), dimension(:), pointer :: event_ids
  character(len=max_length_string) :: eventf_stf_file,eventf_convw_file,eventf_convs_file,eventf_filebase
  integer :: eventf_nev,eventf_stf_col,eventf_stf_nstep,eventf_convw_col,eventf_convw_nstep,eventf_btw_lp_o,eventf_btw_hp_o
  real :: eventf_stf_dt,eventf_convw_dt,eventf_btw_lp_fc,eventf_btw_hp_fc
  type (seismic_event) :: event
  real, dimension(:,:), pointer :: eventf_stf,eventf_convw
  complex, dimension(:), allocatable :: eventf_filter_final,eventf_convw_spec
  complex, dimension(:), pointer :: eventf_convs

  ! ! station filter
  ! type (input_parameter) :: stationf_inpar
  ! character (len=80), dimension() :: stationf_inpar_keys

  ! discrete Fourier transform
  type (discrete_fourier_transform) :: DFT
  integer, dimension(:), pointer :: mdata_jf
  real :: mdata_df
  integer :: mdata_nf
  complex, dimension(:), allocatable :: f

  ! other stuff
  logical :: print_usage_and_stop,terminate_program,create_event_filters,create_station_filters
  integer :: ios,j,n,iev
  double complex :: filter_val,omega_over_omega_c
  double complex :: pi_i = mc_cid*mc_pid


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(event_ids,eventf_stf,eventf_convw,eventf_convs,mdata_jf)

  terminate_program = .false.

!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call init(ap,prog_name,"Generate spectral filter files as used in ASKI from source wavelets or as Butterworth "//&
       "high-/low-/bandpass filters (de)convolved by other wavelets or spectral filters")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-eventf',.true.,"create event filters (describing the SOURCE contribution to the "//&
       "filter which brings the synthetic displacement wavefield to the measured data); argument of -eventf is the "//&
       "parameter file for event filters",'sval','')
  call addOption(ap,'-stationf',.true.,"NOT SUPPORTED YET!! create station filters (describing the STATION COMPONENT "//&
       "contributions to the filter which brings the synthetic displacement wavefield to the measured data); "//&
       "argument of -stationf is the parameter file for station filters",'sval','')
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) goto 3
!
  str = ap.sval."main_parfile"
  main_parfile = str
  if (.level.(.errmsg.ap) == 2) goto 3
!
  print_usage_and_stop = .false.
!
  ! -eventf
  create_event_filters = ap.optset."-eventf"
  if(create_event_filters) then
     str = ap.sval."-eventf"
     eventf_parfile = str
     print_usage_and_stop = print_usage_and_stop .or. .level.(.errmsg.ap) == 2
  end if
!
  ! -stationf
  create_station_filters = ap.optset."-stationf"
  if(create_station_filters) then
     write(*,*) "ERROR: this program does not yet support the option -stationf , sorry"
     print_usage_and_stop = .true.
     ! str = ap.sval."-stationf"
     ! stationf_parfile = str
     ! print_usage_and_stop = print_usage_and_stop .or. .level.(.errmsg.ap) == 2
  end if
!
  if(.not.(create_event_filters .or. create_station_filters)) then
     write(*,*) "ERROR: at least one of options -eventf , -stationf must be set!"
     print_usage_and_stop = .true.
  end if
!
  print_usage_and_stop = print_usage_and_stop .or. .level.(.errmsg.ap) == 2
!
  if(print_usage_and_stop) goto 3
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,10)
!
!------------------------------------------------------------------------
!  get basics
!
   ! setup inversion basics
   call new(errmsg,prog_name)
   call init(invbasics,main_parfile,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   mdata_df = rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP')
   mdata_jf => .ifreq.invbasics
   mdata_nf = ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')
   str = (.inpar.invbasics).sval.'FORWARD_METHOD'
   allocate(f(mdata_nf))
   do j = 1,mdata_nf
      f(j) = getComplexKernelFrequency(str,mdata_df,mdata_jf(j))
   end do
!
   write(*,*) "spectral filters will be created at ",mdata_nf," frequencies according to"
   write(*,*) "  df = ",mdata_df
   write(*,*) "  frequency indices = ",mdata_jf
   if(methodHasComplexKernelFrequency(str)) then
      write(*,*) "for forward method '",trim(str),"', this corresponds to the following complex frequencies:"
      write(*,*) "  ",f
   else
      write(*,*) "corresponding to the real-valued frequencies [Hz]:"
      write(*,*) "  ",real(f)
   end if
   write(*,*) ""
!
!------------------------------------------------------------------------
!  read parfiles
!
  if(create_event_filters) then
     allocate(eventf_filter_final(mdata_nf))
     call create_event_filters_subroutine()
  end if
  if(terminate_program) goto 1
!
  ! if(create_station_filters) call stationf_subroutine()
  ! if(terminate_program) goto 1
!
!
!------------------------------------------------------------------------
!  clean up
!
   write(*,*) "good bye"
!
1  call dealloc(invbasics)
   call dealloc(fuh)
   call dealloc(ap)
   call dealloc(DFT)
   call dealloc(errmsg)
   if(allocated(f)) deallocate(f)
   if(allocated(eventf_filter_final)) deallocate(eventf_filter_final)
!
   stop
!
3  if(.level.(.errmsg.ap)>=1) call print(.errmsg.ap)
   call usage(ap)
   goto 1
!
!------------------------------------------------------------------------
!
contains
!
!------------------------------------------------------------------------
!
  subroutine create_event_filters_subroutine()
    nullify(eventf_stf,eventf_convw,eventf_convs)
!
    write(*,*) "creating event filters now"
    write(*,*) ""
!
    ! read in parameter file and retrieve all relevant information from it
    call eventf_read_parfile()
    if(terminate_program) goto 1
!
    if(create_filters_by_source_wavelet) then
       write(*,*) "creating event filters by source wavelet"
       write(*,*) "  ... given by column ",eventf_stf_col," in file '",trim(eventf_stf_file),"'"
       write(*,*) "  ... having DT, NSTEP = ",eventf_stf_dt,eventf_stf_nstep
!
       ! read in source wavelet as defined in parfile
       errmsg = readAsciiData(eventf_stf_file,get(fuh),eventf_stf,ncol=eventf_stf_col,nrow=eventf_stf_nstep)
       call undo(fuh)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       if(.not.associated(eventf_stf)) then
          write(*,*) "ERROR: even though there was no error reading the source wavelet from file, there are no ",&
               "values returned. This should not happen!"
          goto 2
       end if
       call dealloc(errmsg)
!
       ! transform source wavelet to spectral domain
       ! create an individual DFT object, since the other wavelet(s) involved are allowed to have different
       ! time sampling
       call new(errmsg,prog_name)
       call initiateForwardDFT(DFT,eventf_stf_dt,0,eventf_stf_nstep-1,f,errmsg)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
!
       call transformForwardDFT(DFT,eventf_stf(:,eventf_stf_col),eventf_filter_final,errmsg)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       call dealloc(errmsg)
       call dealloc(DFT)
       deallocate(eventf_stf)
    end if ! create_filters_by_source_wavelet
!
    if(create_filters_by_spectral_butterworth) then
       write(*,*) "creating event filters as spectral butterworth filter"
       if(btw_low_pass_apply) &
            write(*,*) "  ... using low pass filter (order,corner-freq.): ",eventf_btw_lp_o,eventf_btw_lp_fc
       if(btw_high_pass_apply) &
            write(*,*) "  ... using high pass filter (order,corner-freq.): ",eventf_btw_hp_o,eventf_btw_hp_fc
!
       ! create butterworth filter(s)
       do j = 1,mdata_nf
          filter_val = (1.d0,0.d0)
!
          if(btw_low_pass_apply) then
             omega_over_omega_c = dcmplx(f(j))/dcmplx(eventf_btw_lp_fc) ! can use f/fc = (2pi f)/(2pi fc)
             do n = 1,eventf_btw_lp_o
                filter_val = -filter_val*mc_cid/(omega_over_omega_c-exp(pi_i*(dble(n)-0.5d0)/dble(eventf_btw_lp_o)))
             end do ! n
          end if ! btw_low_pass_apply
!
          if(btw_high_pass_apply) then
             omega_over_omega_c = dcmplx(f(j))/dcmplx(eventf_btw_hp_fc) ! can use f/fc here since 1/2pi cancels out
             do n = 1,eventf_btw_hp_o
                filter_val = &
                     filter_val*omega_over_omega_c/(omega_over_omega_c-exp(pi_i*(dble(n)-0.5d0)/dble(eventf_btw_hp_o)))
             end do ! n
          end if ! btw_high_pass_apply
!
          eventf_filter_final(j) = cmplx(filter_val)
       end do ! j
!
    end if ! create_filters_by_spectral_butterworth
!
    if(convolve_with_other_wavelet) then
       if(deconvolve_wavelet_instead) then
          write(*,*) "will deconvolve the event filters by another wavelet"
       else
          write(*,*) "will convolve the event filters by another wavelet"
       end if
       write(*,*) "  ... given by column ",eventf_convw_col," in file '",trim(eventf_convw_file),"'"
       write(*,*) "  ... having DT, NSTEP = ",eventf_convw_dt,eventf_convw_nstep
!
       ! read in (de)convolve wavelet as defined in parfile
       errmsg = readAsciiData(eventf_convw_file,get(fuh),eventf_convw,ncol=eventf_convw_col,nrow=eventf_convw_nstep)
       call undo(fuh)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       if(.not.associated(eventf_convw)) then
          write(*,*) "ERROR: even though there was no error reading the (de)convolve-wavelet from file, there are no ",&
               "values returned. This should not happen!"
          goto 2
       end if
       call dealloc(errmsg)
!
       ! transform (de)convolve wavelet to spectral domain
       ! create an individual DFT object, since the other wavelet(s) involved are allowed to have different
       ! time sampling
       call new(errmsg,prog_name)
       call initiateForwardDFT(DFT,eventf_convw_dt,0,eventf_convw_nstep-1,f,errmsg)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
!
       allocate(eventf_convw_spec(mdata_nf))
       call transformForwardDFT(DFT,eventf_convw(:,eventf_convw_col),eventf_convw_spec,errmsg)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       call dealloc(errmsg)
       call dealloc(DFT)
       deallocate(eventf_convw)
!
       ! do (de)convolution by (dividing)multiplying eventf_filter_final by eventf_convw_spec
       if(deconvolve_wavelet_instead) then
          eventf_filter_final(:) = eventf_filter_final(:) / eventf_convw_spec(:)
       else
          eventf_filter_final(:) = eventf_filter_final(:) * eventf_convw_spec(:)
       end if
!
       deallocate(eventf_convw_spec)
    end if ! convolve_with_other_wavelet
!
    if(multiply_by_other_spectral_filter) then
       if(divide_by_spectral_filter_instead) then
          write(*,*) "will divide the event filters through by another spectral filter"
       else
          write(*,*) "will multiply the event filters by another spectral filter"
       end if
       write(*,*) "  ... given by file '",trim(eventf_convs_file),"'"
!
       ! read in spectrum to be multiplied as defined in parfile
       errmsg = readAsciiData(eventf_convs_file,get(fuh),eventf_convs,ndata=mdata_nf)
       call undo(fuh)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       if(.not.associated(eventf_convs)) then
          write(*,*) "ERROR: even though there was no error reading the spectrum from file, there are no ",&
               "values returned. This should not happen!"
          goto 2
       end if
       call dealloc(errmsg)
!
       ! do multiplication (division) of other spectral filter that was just read in
       if(divide_by_spectral_filter_instead) then
          eventf_filter_final(:) = eventf_filter_final(:) / eventf_convs(:)
       else
          eventf_filter_final(:) = eventf_filter_final(:) * eventf_convs(:)
       end if
       deallocate(eventf_convs)
    end if ! multiply_by_other_spectral_filter
    write(*,*) ""
!
    write(*,*) "the final spectral filter that will be written to file equals (lines correspond to frequencies):"
    do j = 1,mdata_nf
       write(*,*) eventf_filter_final(j)
    end do ! j
    write(*,*) ""
!
    ! write filter files
!
    eventf_filebase = trim((.inpar.invbasics).sval.'PATH_EVENT_FILTER')//'filter_'
    write(*,*) "writing ",eventf_nev," event filter files (all the same content as printed above)"
    write(*,*) "  ... with file basename: '",trim(eventf_filebase),"'"
    write(*,*) "  ... extended by the following event ID's defined in parfile:"
    if(use_all_event_ids) then
       write(*,*) "      ALL events:"
    else
       write(*,*) "      selected events:"
    end if
    do iev = 1,eventf_nev
       write(*,*) "        ",trim(event_ids(iev))
    end do ! iev
!
    ! loop on all events, set filename of filter file, write filter eventf_filter_final to file
    do iev = 1,eventf_nev
       errmsg = writeAsciiData(trim(eventf_filebase)//trim(event_ids(iev)),get(fuh),eventf_filter_final)
       call undo(fuh)
       if(.level.errmsg/=0) call print(errmsg)
       if(.level.errmsg==2) goto 2
       call dealloc(errmsg)
    end do ! iev
!
    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   call dealloc(DFT)
    if(associated(eventf_stf)) deallocate(eventf_stf)
    if(associated(eventf_convw)) deallocate(eventf_convw)
    if(allocated(eventf_convw_spec)) deallocate(eventf_convw_spec)
    if(associated(eventf_convs)) deallocate(eventf_convs)
    write(*,*) ""
    return
!
    ! due to an error, terminate program after return
2   terminate_program = .true.
    call dealloc(errmsg)
    goto 1
  end subroutine create_event_filters_subroutine
!
!------------------------------------------------------------------------
!
  subroutine eventf_read_parfile()
    call createKeywordsInputParameter(eventf_inpar,eventf_inpar_keys)
    write(*,*) "reading parameter file '",trim(eventf_parfile),"'"
    call new(errmsg,prog_name)
    call readSubroutineInputParameter(eventf_inpar,get(fuh),eventf_parfile,errmsg)
    call undo(fuh)
    call addTrace(errmsg,prog_name)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    create_filters_by_source_wavelet = lval(eventf_inpar,'CREATE_FILTERS_BY_SOURCE_WAVELET',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CREATE_FILTERS_BY_SOURCE_WAVELET'),"' of keyword '",&
            "CREATE_FILTERS_BY_SOURCE_WAVELET' in parameter file is no valid logical"
       goto 2
    end if
! 
   create_filters_by_spectral_butterworth = lval(eventf_inpar,'CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH'),"' of keyword '",&
            "CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH' in parameter file is no valid logical"
       goto 2
    end if
!
    if(create_filters_by_source_wavelet .eqv. create_filters_by_spectral_butterworth) then
       write(*,*) "ERROR: CREATE_FILTERS_BY_SOURCE_WAVELET and CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH must not have ",&
            "the same logical value. EXACTLY one of them must be .true."
       goto 2
    end if
!
    if(create_filters_by_source_wavelet) then
       eventf_stf_file = eventf_inpar.sval.'STF_FILE'
       eventf_stf_col = ival(eventf_inpar,'STF_COLUMN_OF_TRACE',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'STF_COLUMN_OF_TRACE'),"' of keyword '",&
               "STF_COLUMN_OF_TRACE' in parameter file is no valid integer"
          goto 2
       end if
       eventf_stf_dt = rval(eventf_inpar,'STF_DT',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'STF_DT'),"' of keyword '",&
               "STF_DT' in parameter file is no valid real value"
          goto 2
       end if
       eventf_stf_nstep = ival(eventf_inpar,'STF_NSTEP',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'STF_NSTEP'),"' of keyword '",&
               "STF_NSTEP' in parameter file is no valid integer"
          goto 2
       end if
    end if ! create_filters_by_source_wavelet
!
    if(create_filters_by_spectral_butterworth) then
       btw_low_pass_apply = lval(eventf_inpar,'BTW_LOW_PASS_APPLY',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_LOW_PASS_APPLY'),"' of keyword '",&
               "BTW_LOW_PASS_APPLY' in parameter file is no valid logical"
          goto 2
       end if
       btw_high_pass_apply = lval(eventf_inpar,'BTW_HIGH_PASS_APPLY',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_HIGH_PASS_APPLY'),"' of keyword '",&
               "BTW_HIGH_PASS_APPLY' in parameter file is no valid logical"
          goto 2
       end if
       if((.not.btw_high_pass_apply).and.(.not.btw_low_pass_apply)) then
          write(*,*) "ERROR: since 'CREATE_FILTERS_BY_SPECTRAL_BUTTERWORTH = .true.' in parfile,",&
               "at least one of BTW_LOW_PASS_APPLY or BTW_HIGH_PASS_APPLY must be set to .true. (both are .false.)"
          goto 2
       end if
!
       if(btw_low_pass_apply) then
          eventf_btw_lp_o = ival(eventf_inpar,'BTW_LOW_PASS_ORDER',iostat=ios)
          if(ios/=0) then
             write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_LOW_PASS_ORDER'),"' of keyword '",&
                  "BTW_LOW_PASS_ORDER' in parameter file is no valid integer"
             goto 2
          end if
          if(eventf_btw_lp_o <= 0) then
             write(*,*) "ERROR: the integer value ",eventf_btw_lp_o," of keyword '",&
                  "BTW_LOW_PASS_ORDER' must be strictly greater than zero"
             goto 2
          end if
          eventf_btw_lp_fc = rval(eventf_inpar,'BTW_LOW_PASS_FC',iostat=ios)
          if(ios/=0) then
             write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_LOW_PASS_FC'),"' of keyword '",&
                  "BTW_LOW_PASS_FC' in parameter file is no valid real value"
             goto 2
          end if
          if(eventf_btw_lp_fc <= 0) then
             write(*,*) "ERROR: the real value ",eventf_btw_lp_fc," of keyword '",&
                  "BTW_LOW_PASS_FC' must be strictly greater than zero"
             goto 2
          end if
       end if ! btw_low_pass_apply
!
       if(btw_high_pass_apply) then
          eventf_btw_hp_o = ival(eventf_inpar,'BTW_HIGH_PASS_ORDER',iostat=ios)
          if(ios/=0) then
             write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_HIGH_PASS_ORDER'),"' of keyword '",&
                  "BTW_HIGH_PASS_ORDER' in parameter file is no valid integer"
             goto 2
          end if
          if(eventf_btw_hp_o <= 0) then
             write(*,*) "ERROR: the integer value ",eventf_btw_hp_o," of keyword '",&
                  "BTW_HIGH_PASS_ORDER' must be strictly greater than zero"
             goto 2
          end if
          eventf_btw_hp_fc = rval(eventf_inpar,'BTW_HIGH_PASS_FC',iostat=ios)
          if(ios/=0) then
             write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'BTW_HIGH_PASS_FC'),"' of keyword '",&
                  "BTW_HIGH_PASS_FC' in parameter file is no valid real value"
             goto 2
          end if
          if(eventf_btw_hp_fc <= 0) then
             write(*,*) "ERROR: the real value ",eventf_btw_hp_fc," of keyword '",&
                  "BTW_HIGH_PASS_FC' must be strictly greater than zero"
             goto 2
          end if
       end if ! btw_high_pass_apply
!
       if(btw_low_pass_apply .and. btw_high_pass_apply) then
          if(eventf_btw_hp_fc > eventf_btw_lp_fc) then
             write(*,*) "ERROR: the corner frequency of the butterworth high pass filter ",&
                  eventf_btw_hp_fc," is greater than the corner frequency of the low pass filter ",&
                  eventf_btw_lp_fc
             goto 2
          end if
       end if
    end if ! create_filters_by_spectral_butterworth
!
    use_all_event_ids = (eventf_inpar.sval.'EVENT_IDS') == 'ALL'
    if(use_all_event_ids) then
       eventf_nev = .nev.(.evlist.invbasics)
       allocate(event_ids(eventf_nev))
       iev = 0
       do while(nextEventSeismicEventList(.evlist.invbasics,event))
          iev = iev+1
          event_ids(iev) = .evid.event
       end do ! while nextEvent
    else ! use_all_event_ids
       eventf_nev = ival(eventf_inpar,'NUMBER_OF_EVENTS',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'NUMBER_OF_EVENTS'),"' of keyword '",&
               "NUMBER_OF_EVENTS' in parameter file is no valid integer value"
          goto 2
       end if
       if(eventf_nev <= 0) then
          write(*,*) "ERROR: the number of events ",eventf_nev," in the parfile must be strictly greater than 0"
          goto 2
       end if
       event_ids => svecp(eventf_inpar,'EVENT_IDS',eventf_nev,iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'EVENT_IDS'),"' of keyword '",&
               "EVENT_IDS' in parameter file is no valid white-space separated vector of ",eventf_nev,&
               " character strings"
          goto 2
       end if
       if(.not.associated(event_ids)) then
          write(*,*) "ERROR: even though there was no error when reading EVENT_IDS from parfile, there is no ",&
               "list of event ids returned... this should not happen"
          goto 2
       end if
       ! check the events
       do iev = 1,eventf_nev
          errmsg = searchEventidSeismicEventList(.evlist.invbasics,trim(event_ids(iev)))
          if(.level.errmsg/=0) then
             call print(errmsg)
             write(*,*) "ERROR: ",iev,"'th event ID '",trim(event_ids(j)),"' of list of event IDs in parfile is ",&
                  "no valid event ID, i.e. it is not present in the ASKI event list"
             goto 2
          end if
          call dealloc(errmsg)
       end do ! iev
    end if ! use_all_event_ids
!
    convolve_with_other_wavelet = lval(eventf_inpar,'CONVOLVE_WITH_OTHER_WAVELET',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CONVOLVE_WITH_OTHER_WAVELET'),"' of keyword '",&
            "CONVOLVE_WITH_OTHER_WAVELET' in parameter file is no valid logical"
       goto 2
    end if
    if(convolve_with_other_wavelet) then
       eventf_convw_file = eventf_inpar.sval.'CONV_WAVEL_FILE'
       eventf_convw_col = ival(eventf_inpar,'CONV_WAVEL_COLUMN_OF_TRACE',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CONV_WAVEL_COLUMN_OF_TRACE'),"' of keyword '",&
               "CONV_WAVEL_COLUMN_OF_TRACE' in parameter file is no valid integer"
          goto 2
       end if
       eventf_convw_dt = rval(eventf_inpar,'CONV_WAVEL_DT',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CONV_WAVEL_DT'),"' of keyword '",&
               "CONV_WAVEL_DT' in parameter file is no valid real value"
          goto 2
       end if
       eventf_convw_nstep = ival(eventf_inpar,'CONV_WAVEL_NSTEP',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'CONV_WAVEL_NSTEP'),"' of keyword '",&
               "CONV_WAVEL_NSTEP' in parameter file is no valid integer"
          goto 2
       end if
       deconvolve_wavelet_instead = lval(eventf_inpar,'DECONVOLVE_WAVEL_INSTEAD',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'DECONVOLVE_WAVEL_INSTEAD'),"' of keyword '",&
               "DECONVOLVE_WAVEL_INSTEAD' in parameter file is no valid logical"
          goto 2
       end if
    end if ! convolve_with_other_wavelet
!
    multiply_by_other_spectral_filter = lval(eventf_inpar,'MULTIPLY_BY_OTHER_SPECTRAL_FILTER',iostat=ios)
    if(ios/=0) then
       write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'MULTIPLY_BY_OTHER_SPECTRAL_FILTER'),"' of keyword '",&
            "MULTIPLY_BY_OTHER_SPECTRAL_FILTER' in parameter file is no valid logical"
       goto 2
    end if
    if(multiply_by_other_spectral_filter) then
       eventf_convs_file = eventf_inpar.sval.'CONV_SPEC_FILE'
       divide_by_spectral_filter_instead = lval(eventf_inpar,'DECONVOLVE_SPEC_INSTEAD',iostat=ios)
       if(ios/=0) then
          write(*,*) "ERROR: value '",trim(eventf_inpar.sval.'DECONVOLVE_SPEC_INSTEAD'),"' of keyword '",&
               "DECONVOLVE_SPEC_INSTEAD' in parameter file is no valid logical"
          goto 2
       end if
    end if ! multiply_by_other_spectral_filter
!
    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   call dealloc(eventf_inpar)
    return
!
    ! due to an error, terminate program after return
2   terminate_program = .true.
    call dealloc(errmsg)
    goto 1
  end subroutine eventf_read_parfile
!
!------------------------------------------------------------------------
!
  ! subroutine stationf_subroutine()    
  ! end subroutine stationf_subroutine
!
end program createSpectralFilters
