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
program spec2timeKernels
  use inversionBasics
  use iterationStepBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use dataModelSpaceInfo
  use componentTransformation
  use spectralWaveformKernel
  use timeWaveformKernel
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,skernel_filebase,tkernel_filebase,str
  character(len=max_length_string), dimension(:), pointer :: str_vec

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg,errmsg2
  character(len=16) :: myname = 'spec2timeKernels'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics
  
  type (data_model_space_info) :: dmspace
  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  ! spectral kernel
  type (spectral_waveform_kernel) :: skernel
  integer, dimension(:), pointer :: jf
  real :: df_skernel
  integer :: ntot_kernel
  ! time kernel
  real :: t0,dt
  integer :: nwin,iwin,njt,ncomp,jcomp
  integer, dimension(:), pointer :: nt1,nt2
  integer, dimension(:), allocatable :: jt
  character(len=character_length_component), dimension(:), pointer :: comp_path
  character(len=16) :: comp_path_write
  integer :: iparam
  character(len=character_length_param), dimension(:), pointer :: param

  ! filtering
  integer :: nfreq_measured
  real :: df_measured
  integer, dimension(:), pointer :: map_jf_2_filter_index
  complex, dimension(:), pointer :: event_filter,station_comp_filter
  complex, dimension(:,:), allocatable :: filter

  logical :: use_dmspace,compute_one_path,kernels_on_wp,&
       apply_event_filter,apply_station_filter,apply_filter,&
       stop_after_command_line,next,terminate_program,&
       check_nt1_nt2,check_ipath
  integer :: npath,ipath1,ipath2,ipath,lu,j

  nullify(str_vec,jf,nt1,nt2,comp_path,param,map_jf_2_filter_index,event_filter,station_comp_filter)
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Compute a set of time kernel files FROM EXISTING spectral kernel files (by inverse "//&
       "Fourier transform). Supports both types of spectral kernel files: pre-integrated or on wavefield points "//&
       "(flag -wp). The set of kernels can be characterized in two different ways: 'way 1' "//&
       "(one path only) and 'way 2' (from dmspace file). Time "//&
       "kernels are computed at times t = t0 + jt*dt , nt1 <= jt <= nt2 .")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-evid',.true.,"(way 1) defines the event id of the one path (must belong to an event in"//&
       " main event list), mandatory for 'way 1'",'sval','')
  call addOption(ap,'-staname',.true.,"(way 1) defines the station name of the one path (must belong to a "//&
        "station in main station list, mandatory for 'way 1'",'sval','')
  call addOption(ap,'-comp',.true.,"(way 1) vector of receiver components for which the kernel should be "//&
       "computed, mandatory for 'way 1'; valid components are: "//all_valid_components,'svec','')
  call addOption(ap,'-param',.true.,"(way 1) vector of parameter names for which the kernel should be "//&
       "computed, mandatory for 'way 1'",'svec','')
  call addOption(ap,'-dmspace',.true.,"(way 2) data model space input file to define a set of paths (and "//&
       "respective receiver components) as well as model parameters, mandatory for 'way 2'",'sval','')
  call addOption(ap,'-ipath1',.true.,"(way 2) first index of the path loop, optional "//&
       "for 'way 2'",'ival','1')
   call addOption(ap,'-ipath2',.true.,"(way 2) last index of the path loop, optional for 'way 2'",'ival',&
        'max_number_of_paths')
  call addOption(ap,"-dt",.true.,"(mandatory) global time step of time discretization; must be > 0","rval","")
  call addOption(ap,"-nt1",.true.,"(mandatory) vector of length nwin, giving starting indices of nwin time "//&
       "windows (must have same length as vector given by -nt2)","ivec","")
  call addOption(ap,"-nt2",.true.,"(mandatory) vector of length nwin, giving end indices of nwin time "//&
       "windows (must have same length as vector given by -nt1)","ivec","")
  call addOption(ap,"-t0",.true.,"(optional) global time shift which is added to all times defined by dt,nt1,nt2",&
       "rval","0.0")
   call addOption(ap,'-wp',.false.,"(optional) If set, then 'ON-WP' spectral kernel files are transformed. If not set, "//&
        "normal kernel files (pre-integrated) are transformed.")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   stop_after_command_line = .false.
!
   main_parfile = ap.sval.'main_parfile'
!
   kernels_on_wp = ap.optset.'-wp'
!
   use_dmspace = ap.optset.'-dmspace'
   compute_one_path = (ap.optset.'-evid') .and. (ap.optset.'-staname') .and. (ap.optset.'-comp') .and. &
        (ap.optset.'-param')
!
  if(.not. (use_dmspace.or.compute_one_path)) then
     write(*,*) "ERROR: please use either one path ('way 1') or data model space file ('way 2')"
     stop_after_command_line = .true.
  end if
!
  if(use_dmspace .and. compute_one_path) then
     write(*,*) "ERROR: please use either one path or data model space file"
     stop_after_command_line = .true.
  end if
!
  if(compute_one_path .and. ((ap.optset.'-ipath1').or.(ap.optset.'-ipath2'))) then
     write(*,*) "ERROR: -ipath1, -ipath2 only to be used together with -dmspace"
     stop_after_command_line = .true.
  end if
!
  if(use_dmspace .and. ((ap.optset.'-evid').or. (ap.optset.'-staname').or.(ap.optset.'-comp').or.(ap.optset.'-param'))) then
     write(*,*) "ERROR: -evid, -staname, -comp, -param only to be used when NOT using -dmspce"
     stop_after_command_line = .true.
  end if
!
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   if(stop_after_command_line) then
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) then; call print(errmsg); endif
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   apply_event_filter = lval(.inpar.invbasics,'APPLY_EVENT_FILTER')
   apply_station_filter = lval(.inpar.invbasics,'APPLY_STATION_FILTER')
   apply_filter = apply_event_filter .or. apply_station_filter
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg)
   if (.level.errmsg /= 0) then; call print(errmsg); endif
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!------------------------------------------------------------------------
!  other stuff to do beforehand
!
   if(use_dmspace) then
      stop_after_command_line = .false.
!
      dmspace_file = ap.sval.'-dmspace'
      call new(errmsg,myname)
      call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
           .ifreq.iterbasics,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',.ncell.(.invgrid.iterbasics),&
           .intw.iterbasics,trim(dmspace_file),get(fuh),errmsg)
      call undo(fuh)
      call addTrace(errmsg,myname)
      if (.level.errmsg /= 0) call print(errmsg)
      !call print(errmsg)
      if (.level.errmsg == 2) goto 1
      if(.ndata.dmspace == 0 .or. .nmval.dmspace == 0) then
         write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace,&
              "model values contained in the data_model_space_info object. using -dmspace there must BOTH, ",&
              "data samples and model parameters, be defined. Printing now the error message returned from ",&
              "createFromFileDataModelSpaceInfo"
         call print(errmsg)
         goto 1
      end if
      call dealloc(errmsg)
      npath = .npath.dmspace ! cannot be 0, since .ndata.dmspace /= 0
!
      ! define model parameters
      ! param should be returned associated, since .nmval.dmspace /= 0 (assured above)
      param => getAllDifferentParamModelValuesDataModelSpaceInfo(dmspace)
!
      check_ipath = .true.
      ! define ipath1
      if(.not.(ap.optset.'-ipath1')) then
         ipath1 = 1
      else
         ipath1 = ap.ival.'-ipath1'
         if (.level.(.errmsg.ap) == 2) then
            call print(.errmsg.ap)
            call usage(ap)
            goto 1
         end if
         if(ipath1<1) then
            write(*,*) "ERROR: first path index = ",ipath1,", given by -ipath1 must be at least 1"
            stop_after_command_line = .true.
            check_ipath = .false.
         end if
      end if
      ! define ipath2
      if(.not.(ap.optset.'-ipath1')) then
         ipath2 = npath
      else
         ipath2 = ap.ival.'-ipath2'
         if (.level.(.errmsg.ap) == 2) then
            call print(.errmsg.ap)
            call usage(ap)
            goto 1
         end if
         if(ipath2>npath) then
            write(*,*) "ERROR: last path index = ",ipath2,", given by -ipath2 is greater than maximum number of "//&
                 "paths (",npath,") contained in data-model-space info defined by file '"//trim(dmspace_file)//"'"
            stop_after_command_line = .true.
            check_ipath = .false.
         end if
      end if
      ! check if ipath1<=ipath2
      if(check_ipath) then
         if(ipath1>ipath2) then
            write(*,*) "ERROR: first path index ( = ",ipath1,") is greater than last path index ( = ",ipath2,")"
            stop_after_command_line = .true.
         end if
      end if
!
   else ! use_dmspace
!
      stop_after_command_line = .false.
!
      ! previously checked:
      ! if program comes here: compute_one_path = ap.optset.'-evid' .and. ap.optset.'-staname' .and. ap.optset.'-comp' .and. ap.optset.'-param' == .true.
!
      npath = 1; ipath1 = 1; ipath2 = 1
!
      ! -evid
      ! check if event ID is valid
      str = ap.sval.'-evid'
      evid = str
      errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid)
      if(.level. errmsg/=0) then
         write(*,*) "event ID '"//trim(evid)//"' (input string of option -evid) is not contained in event list"
         stop_after_command_line = .true.
      end if
      call dealloc(errmsg)
      ! -stname
      ! check if station name is valid
      str = ap.sval.'-staname'
      staname = str
      errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname)
      if(.level. errmsg/=0) then
         write(*,*) "station name '"//trim(staname)//"' (input string of option -stname) is not contained in station list"
         stop_after_command_line = .true.
      end if
      call dealloc(errmsg)
!
      ! comp
      str_vec => ap.svec.'-comp'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
      if(.not.associated(str_vec)) then
         write(*,*) "ERROR: for some reason, there is no list of station components returned by argument parser, "//&
              "even though there was no error parsing argument -comp. This is strange..."
         write(*,*) ""
         call usage(ap)
         goto 1
      end if
      allocate(comp_path(size(str_vec)))
      do jcomp = 1,size(str_vec)
         comp_path(jcomp) = str_vec(jcomp)
      end do
      deallocate(str_vec)
      if(.not.allValidComponents(comp_path,i_invalid=jcomp)) then
         write(*,*) "ERROR: ",jcomp,"'th component '"//trim(comp_path(jcomp))//"' not valid. Valid components are '"//&
              all_valid_components//"'"
         stop_after_command_line = .true.
      end if
!
      ! param
      str_vec => ap.svec.'-param'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
      if(.not.associated(str_vec)) then
         write(*,*) "ERROR: for some reason, there is no list of model parameters returned by argument parser, "//&
              "even though there was no error parsing argument -param. This is strange..."
         write(*,*) ""
         call usage(ap)
         goto 1
      end if
      allocate(param(size(str_vec)))
      do iparam = 1,size(str_vec)
         param(iparam) = str_vec(iparam)
         if(.not.validParamModelParametrization((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param(iparam))) then
            write(*,*) "ERROR: ",iparam,"'th given model parameter '",trim(param(iparam)),"' not valid in current "//&
                 "parametrization '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//"'. Valid "//&
                 "parameterizations(parameters are ",all_valid_pmtrz_param
            stop_after_command_line = .true.
         end if
      end do ! iparam
      deallocate(str_vec)
!
   end if ! use_dmspace
!
   ! now treat options regarding time kernel computation
!
   ! t0
   if(.not.(ap.optset.'-t0')) then
      t0 = 0.0
   else
      t0 = ap.rval.'-t0'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         stop_after_command_line = .true.
      end if
   end if
!
   ! dt
   if(.not.(ap.optset.'-dt')) then
      write(*,*) "ERROR: please indicate the time step '-dt'"
      stop_after_command_line = .true.
   else
      dt = ap.rval.'-dt'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         stop_after_command_line = .true.
      end if
      if(dt .le. 0.) then
         write(*,*) "ERROR: global time step indicated by '-dt' =",dt,". This value must be strictly positive!"
         stop_after_command_line = .true.
      end if
   end if
!
   check_nt1_nt2 = .true.
   ! nt1
   if(.not.(ap.optset.'-nt1')) then
      write(*,*) "ERROR: please indicate starting indices of nwin>0 time windows by option '-nt1'"
      stop_after_command_line = .true.
      check_nt1_nt2 = .false.
   else
      nt1 => ap.ivec.'-nt1'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         check_nt1_nt2 = .false.
         stop_after_command_line = .true.
      end if
   end if
   ! nt2
   if(.not.(ap.optset.'-nt2')) then
      write(*,*) "ERROR: please indicate end indices of nwin>0 time windows by option '-nt2'"
      stop_after_command_line = .true.
      check_nt1_nt2 = .false.
   else
      nt2 => ap.ivec.'-nt2'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         check_nt1_nt2 = .false.
         stop_after_command_line = .true.
      end if
   end if
   ! check if nt2 >= nt1 for all time windows
   if(check_nt1_nt2) then
      nwin = size(nt1)
      if(size(nt2)/=nwin .or. nwin <= 0) then
         write(*,*) "ERROR: size of vector -nt1 = ",nwin," and size of vector -nt2 =",size(nt2)," must be equal and > 0"
         stop_after_command_line = .true.
      else
         do iwin=1,nwin
            if(nt1(iwin)>nt2(iwin)) then
               write(*,*) "ERROR: for ",iwin,"'th time window (out of ",nwin,") nt1,nt2 do not fulfill nt1 <= "//&
                    "nt2 :  nt1,nt2 =",nt1(iwin),nt2(iwin)
               stop_after_command_line = .true.
            end if
         end do ! iwin
      end if ! size(nt2)/=nwin .or. nwin <= 0
   end if ! check_nt1_nt2
!
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   if(stop_after_command_line) then
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
!
   call document(ap)
   write(*,*) ""
!
   ! now define vector of all time indices (all time windows)
   ! check whether time windows to overlap or not, remove duplicate indices
   njt = sum( (nt2-nt1) + 1)
   allocate(jt(njt))
!
   ! put all indices of first time window into vector jt
   njt = nt2(1)-nt1(1) + 1
   jt(1:njt) = (/ (j,j=nt1(1),nt2(1)) /)
   ! afterwards loop on the rest of the time windows, only adding time indices which are not yet present in jt
   do iwin = 2,nwin
      do j = nt1(iwin),nt2(iwin)
         if(.not.any(jt(1:njt)==j)) then
            njt = njt + 1
            jt(njt) = j
         end if
      end do ! j
   end do ! iwin
   ! now, value njt indicates sensible values in array jt!
   if(njt < size(jt)) then
      write(*,*) "WARNING: time windows overlap: there are ",size(jt)-njt,&
           " duplicate time indices (out of a total of ",size(jt),"), which are ignored "//&
           "(no duplicate kernel values written to time kernel file)"
   end if
!
!------------------------------------------------------------------------
!  transform spectral to time kernels
!
   ! set base of all time and spectral kernel filenames
   if(kernels_on_wp) then
      tkernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'time_kernel_ON-WP_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
      skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_ON-WP_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
      ntot_kernel = .ntot.(.wp.iterbasics)
   else
      tkernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'time_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
      skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'spectral_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
      ntot_kernel = .ncell.(.invgrid.iterbasics)
   end if
!
   nfreq_measured = (.inpar.invbasics).ival.'MEASURED_DATA_NUMBER_OF_FREQ'
   df_measured = (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'
   nullify(event_filter,station_comp_filter,map_jf_2_filter_index)
!
   terminate_program = .false.
!
  if(use_dmspace) call spec2TimeKernelsDmspace()
  if(terminate_program) goto 1 
!
  if(compute_one_path) call spec2TimeKernelsOnePath()
  if(terminate_program) goto 1 
!------------------------------------------------------------------------
!  clean up
!
1  if(associated(str_vec)) deallocate(str_vec)
   if(associated(nt1)) deallocate(nt1)
   if(associated(nt2)) deallocate(nt2)
   if(allocated(jt)) deallocate(jt)
   if(associated(comp_path)) deallocate(comp_path)
   if(associated(param)) deallocate(param)
   if(associated(map_jf_2_filter_index)) deallocate(map_jf_2_filter_index)
   if(associated(event_filter)) deallocate(event_filter)
   if(associated(station_comp_filter)) deallocate(station_comp_filter)
   if(allocated(filter)) deallocate(filter)
   call dealloc(errmsg); call dealloc(errmsg2)
   call dealloc(dmspace)
   call dealloc(skernel)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
!
  stop
!
contains
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine spec2TimeKernelsDmspace()
!
   write(*,*) "will now compute time '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//"' kernels of parameters ",&
        param//','," for ",ipath2-ipath1+1," path(s), transforming from spectral kernels (assuming them to exist)"
   if(kernels_on_wp) then
      write(*,*) "will transform plain kernel values on wavefield points"
   else
      write(*,*) "will transform pre-integrated kernel values on inversion grid cells"
   end if
   write(*,*) "for ",nwin," time windows and a total of ",njt," time samples: "
   write(*,*) "   time indices jt = ",jt(1:njt)
   write(*,*) "   i.e. times t = t0 + jt*dt = ",t0+jt(1:njt)*dt
   write(*,*) ""
   write(*,*) " path index | event ID        | station name    | station components         |"
   write(*,*) "------------+-----------------+-----------------+----------------------------|"
!
   ! incremet index ipath manually, for output on screen
   ipath = ipath1-1
   ! loop on all required paths
   do while(nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path,ipath_start=ipath1,ipath_end=ipath2))
      ipath = ipath + 1
      ncomp = size(comp_path)
!
      ! use a new error message for each path
      call new(errmsg,myname)
!
      ! if there is a next path, then comp_path must be associated here, otherwise module dataModelSpaceInfo 
      ! is corrupt, so trust the output here and do not check again if comp_path is associated

      write(comp_path_write,*) comp_path//","

      write(*,"(i12,a,a15,a,a15,2a)") ipath," | ","'"//trim(evid)//"'"," | ","'"//trim(staname)//"'"," |",trim(comp_path_write)
      !write(*,*) "# COMPUTE TIME KERNEL for ",ipath,"'th path: event '",trim(evid),"', station '",trim(staname),"'"
!
      call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param,ntot_kernel,&
           comp_path,errmsg,kernels_on_wp)
      if(.level.errmsg == 2) then
         call print(errmsg)
         goto 2
      end if
      call initialReadSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg)
      if(.level.errmsg == 2) then
         call print(errmsg)
         goto 2
      end if
!
      ! check if frequency discretization of this spectral kernel file is consistent with
      ! frequency discretization of the filters (i.e. the main frequency discretization):
      jf => .jf.skernel
      if(.not.associated(jf)) then
         write(*,*) "ERROR: no frequencies contained in kernel file"
         goto 2
      end if
      df_skernel = .df.skernel
      if( abs(df_measured-df_skernel)/df_measured > 1.e-4 ) then
         write(*,*) "ERROR: frequency step of spectral kernel ( = ",df_skernel,&
              ") differs from frequency step of measured data / filters ( = ",df_measured,&
              ") by more than 0.01 percent"
         goto 2
      end if
!
      ! if required, read in filter values for this event and all components of this station
      ! and select required filter values for this kernel
      if(apply_filter) then
         if(associated(map_jf_2_filter_index)) deallocate(map_jf_2_filter_index)
         map_jf_2_filter_index => mapIfreq2ArrayIndex(invbasics,jf)
         if(.not.associated(map_jf_2_filter_index)) then
            write(*,*) "ERROR: some frequency indices contained in kernel file are not contained "//&
                 "in frequency indices of measured data / filters"
            write(*,*) "jf_kernel = ",jf
            write(*,*) "jf_measured = ",.ifreq.invbasics
            goto 2
         end if
!
         if(allocated(filter)) deallocate(filter)
         allocate(filter(.nfreq.skernel,ncomp)) ! assuming .nfreq.skernel == size(jf) (which I can assume, this is the rules of the module)
         filter = (1.,0.)
!
         if(apply_event_filter) then
            if(associated(event_filter)) deallocate(event_filter)
            errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_EVENT_FILTER')//"filter_"//trim(evid),get(fuh),&
                 event_filter,ndata=nfreq_measured)
            call undo(fuh)
            call addTrace(errmsg2,myname)
            if(.level.errmsg2 /= 0) call print(errmsg2)
            if(.level.errmsg2 == 2) goto 2
            call dealloc(errmsg2)
!
            do jcomp = 1,ncomp
               filter(:,jcomp) = filter(:,jcomp)*event_filter(map_jf_2_filter_index)
            end do
         end if ! apply_event_filter
!
         if(apply_station_filter) then
            do jcomp = 1,ncomp
               if(associated(station_comp_filter)) deallocate(station_comp_filter)
               errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_STATION_FILTER')//"filter_"//trim(staname)//&
                    "_"//trim(comp_path(jcomp)),get(fuh),station_comp_filter,ndata=nfreq_measured)
               call undo(fuh)
               call addTrace(errmsg2,myname)
               if(.level.errmsg2 /= 0) call print(errmsg2)
               if(.level.errmsg2 == 2) goto 2
               call dealloc(errmsg2)
!
               filter(:,jcomp) = event_filter(map_jf_2_filter_index)*station_comp_filter(map_jf_2_filter_index)
            end do ! jcomp
         end if ! apply_station_filter
      end if ! apply_filter
!
      if(apply_filter) then
         call transformSpectralToTimeWaveformKernel(skernel,dt,jt(1:njt),&
              trim(tkernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg,filter=filter,t0=t0)
      else
         call transformSpectralToTimeWaveformKernel(skernel,dt,jt(1:njt),&
              trim(tkernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg,t0=t0)
      end if
      call undo(fuh)
      if(.level.errmsg == 2) then
         call print(errmsg)
         goto 2
      end if
!
      call finalReadSpectralWaveformKernel(skernel,lu)
      call add(fuh,lu)
      call dealloc(skernel)
!
      if(.level.errmsg /= 0) call print(errmsg)
      call dealloc(errmsg)
   end do ! while(nextPath)
!
   ! clean up before returning
1  if(allocated(filter)) deallocate(filter)
   if(associated(event_filter)) deallocate(event_filter)
   if(associated(station_comp_filter)) deallocate(station_comp_filter)
!
   ! if code comes here, return normally
   return
!
   ! due to an error, terminate program after return
2   terminate_program = .true.
   ! if there was an error inside the loop on all paths, need to reset the iterator (and deallocate all
   ! arrays involved) before returning
   next = nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path,reset=.true.)
   call dealloc(skernel)
   call dealloc(errmsg); call dealloc(errmsg2)
   goto 1
 end subroutine spec2TimeKernelsDmspace
!
!-----------------------------------------------------------------------------------------------------------------
!
 subroutine spec2TimeKernelsOnePath()
!
   write(*,*) "will now compute time '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//"' kernels of parameters ",&
        param//','," for one path, transforming from spectral kernels (assuming them to exist):"
   write(*,*)"event ID = ",trim(evid)
   write(*,*)"station name = ",trim(staname)," , receiver components = ",comp_path//","
   if(kernels_on_wp) then
      write(*,*) "will transform plain kernel values on wavefield points"
   else
      write(*,*) "will transform pre-integrated kernel values on inversion grid cells"
   end if
   write(*,*) "for ",nwin," time windows and a total of ",njt," time samples: "
   write(*,*) "   time indices jt = ",jt(1:njt)
   write(*,*) "   i.e. times t = t0 + jt*dt = ",t0+jt(1:njt)*dt
   write(*,*) ""
!
   ! use one error message for everything
   call new(errmsg,myname)
!
   call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param,ntot_kernel,&
        comp_path,errmsg,kernels_on_wp)
   if(.level.errmsg == 2) then
      call print(errmsg)
      goto 2
   end if
   call initialReadSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg)
   if(.level.errmsg == 2) then
      call print(errmsg)
      goto 2
   end if
!
   ! check if frequency discretization of this spectral kernel file is consistent with
   ! frequency discretization of the filters (i.e. the main frequency discretization):
   jf => .jf.skernel
   if(.not.associated(jf)) then
      write(*,*) "ERROR: no frequencies contained in kernel file"
      goto 2
   end if
   df_skernel = .df.skernel
   if( abs(df_measured-df_skernel)/df_measured > 1.e-4 ) then
      write(*,*) "ERROR: frequency step of spectral kernel ( = ",df_skernel,&
           ") differs from frequency step of measured data / filters ( = ",df_measured,&
           ") by more than 0.01 percent"
      goto 2
   end if
!
   ! if required, read in filter values for this event and all components of this station
   ! and select required filter values for this kernel
   if(apply_filter) then
      if(associated(map_jf_2_filter_index)) deallocate(map_jf_2_filter_index)
      map_jf_2_filter_index => mapIfreq2ArrayIndex(invbasics,jf)
      if(.not.associated(map_jf_2_filter_index)) then
         write(*,*) "ERROR: some frequency indices contained in kernel file are not contained "//&
              "in frequency indices of measured data / filters"
         write(*,*) "jf_kernel = ",jf
         write(*,*) "jf_measured = ",.ifreq.invbasics
         goto 2
      end if
!
      if(allocated(filter)) deallocate(filter)
      allocate(filter(.nfreq.skernel,ncomp)) ! assuming .nfreq.skernel == size(jf) (which I can assume, this is the rules of the module)
      filter = (1.,0.)
!
      if(apply_event_filter) then
         if(associated(event_filter)) deallocate(event_filter)
         errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_EVENT_FILTER')//"filter_"//trim(evid),get(fuh),&
              event_filter,ndata=nfreq_measured)
         call undo(fuh)
         call addTrace(errmsg2,myname)
         if(.level.errmsg2 /= 0) call print(errmsg2)
         if(.level.errmsg2 == 2) goto 2
         call dealloc(errmsg2)
!
         do jcomp = 1,ncomp
            filter(:,jcomp) = filter(:,jcomp)*event_filter(map_jf_2_filter_index)
         end do
      end if ! apply_event_filter
!
      if(apply_station_filter) then
         do jcomp = 1,ncomp
            if(associated(station_comp_filter)) deallocate(station_comp_filter)
            errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_STATION_FILTER')//"filter_"//trim(staname)//&
                 "_"//trim(comp_path(jcomp)),get(fuh),station_comp_filter,ndata=nfreq_measured)
            call undo(fuh)
            call addTrace(errmsg2,myname)
            if(.level.errmsg2 /= 0) call print(errmsg2)
            if(.level.errmsg2 == 2) goto 2
            call dealloc(errmsg2)
!
            filter(:,jcomp) = event_filter(map_jf_2_filter_index)*station_comp_filter(map_jf_2_filter_index)
         end do ! jcomp
      end if ! apply_station_filter
   end if ! apply_filter
!
   if(apply_filter) then
      call transformSpectralToTimeWaveformKernel(skernel,dt,jt(1:njt),&
           trim(tkernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg,filter=filter,t0=t0)
   else
      call transformSpectralToTimeWaveformKernel(skernel,dt,jt(1:njt),&
           trim(tkernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg,t0=t0)
   end if
   call undo(fuh)
   if(.level.errmsg == 2) then
      call print(errmsg)
      goto 2
   end if
!
   call finalReadSpectralWaveformKernel(skernel,lu)
   call add(fuh,lu)
   call dealloc(skernel)
!
   if(.level.errmsg /= 0) call print(errmsg)
   call dealloc(errmsg)
!
   ! clean up before returning
1  if(allocated(filter)) deallocate(filter)
   if(associated(event_filter)) deallocate(event_filter)
   if(associated(station_comp_filter)) deallocate(station_comp_filter)
!
   ! if code comes here, return normally
   return
!
   ! due to an error, terminate program after return
2  terminate_program = .true.
   call dealloc(skernel)
   call dealloc(errmsg); call dealloc(errmsg2)
   goto 1
 end subroutine spec2TimeKernelsOnePath
!
end program spec2timeKernels
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'    spec2timeKernels [-h] [-evid event_id -stname station_name -comp "ncomp comp" -param "nparam param"]'
!   print *,'       [-dmspce dmspace_file [-ipath1 path1] [-ipath2 path2]]'
!   print *,'       -dt time_step -nwin ntime_windows -nt1 nt1_string -nt2 nt2_string [-t0 tzero] [-wp] main_parfile'
!   print *,''
!   print *,'Description:'
!   print *,'Compute time_kernel files FROM EXISTING spectral_kernel files (by inverse Fourier transform), assuming'
!   print *,'that the spectral kernel files were already produced by program computeKernels.'
!   print *,'Can either handle pre-integrated kernel values on inversion grid cells, or original kernel values on'
!   print *,'wavefield points (flag -wp): MAKE SURE to produce same type of spectral kernel files before!'
!   print *,'The time kernels are computed for the time windows as specified by options t0,dt,nwin,nt1,nt2:'
!   print *,'The total set of time samples consists of nwin time windows. nt1,nt2 are strings containing nwin integers'
!   print *,'defining the start (nt1) and end (nt2) of each time window by an index. Times compute as t = t0 + jt*dt'
!   print *,'for nt1 <= jt <= nt2.'
!   print *,''
!   print *,'THERE ARE TWO POSSIBLE WAYS TO DEFINE A SET OF KERNELS THAT ARE TRANSFORMED:'
!   print *,'  (way 1) compute kernel for only one path, defined by eventID and station name and component(s) using'
!   print *,'          options -evid  , -stname , -comp and -param'
!   print *,''
!   print *,'  (way 2) use flag -dmspce in connection with optional range definition of the path index (flags -ipath1 -ipath2)'
!   print *,'          in order to define a subset of the paths contained in the given data-model-space description.'
!   print *,'          If the upper (lower) limit of the path range is not defined (i.e. -ipath1 (-ipath2) not set), the maximum'
!   print *,'          (minimum) possible value is used.'
!   print *,'          Then, all kernels for the defined range of paths in the given data-model-space description are computed'
!   print *,'          at the respective receiver components and the model parameters defined in the data-model-space description'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print this help message and do not do anything else'
!   print *,''
!   print *,'-wp    : if set, then plain kernel values on WAVEFIELD POINTS are handled. Otherwise (if not set),'
!   print *,'         the standard pre-integrated kernels on inversion grid cells are produced.'
!   print *,''
!   print *,'SET OF KERNELS TO BE TRANSFORMED'
!   print *,'  (way 1)'
!   print *,'    -evid event_id : defines the event id of the one path (must belong to an event in main event list)'
!   print *,'    -stname station_name   : defines the station name of the one path (must belong to a station in '//&
!        'main station list)'
!   print *,'    -comp "n comp1..compn" : defines a list of n receiver components for which the kernels should be computed '
!   print *,'                         (the string must start with the number n, followed by n component names). '
!   print *,'                         Valid component names are: '//all_valid_components//' (way 1)'
!   print *,'    -param "n param1..paramn" : defines a list of n parameter names (e.g. "2 vp vs") for which the kernels should be '
!   print *,'                         computed (the string must start with the number n, followed by n valid parameter names). '
!   print *,'                         Valid parameter names are: '//all_valid_pmtrz_param//' (way 1)'
!   print *,'  (way 2)'
!   print *,"    -dmspce  : data model space input file to define a set of paths"
!   print *,'    -ipath1 path1 (optional) : path1 is the first index of the path loop. By default, path1 = 1'
!   print *,'    -ipath2 path2 (optional) : path2 is the last index of the path loop. By default, path2 = max_number_of_paths '//&
!        'as to data-model-space description'
!   print *,''
!   print *,'TIME WAVEFORM KERNELS'
!   print *,'-t0 tzero  : optional global time shift which is added to all times defined by dt,nt1,nt2 (default t0=0)'
!   print *,'-dt time_step  : global time step of time discretization'
!   print *,'-nwin ntime_windows  : integer number of time windows'
!   print *,'-nt1 nt1_string  : string containing ntime_windows space separated integers defining nt1 for each time window'
!   print *,'-nt2 nt2_string  : string containing ntime_windows space separated integers defining nt2 for each time window'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
