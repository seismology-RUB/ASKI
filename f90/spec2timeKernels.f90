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
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=300) :: parfile,string,dmspace_file,skernel_filebase,tkernel_filebase

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg,errmsg2
  character(len=16) :: myname = 'spec2timeKernels'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics
  
  type (data_model_space_info) :: dmspace
  character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  ! spectral kernel
  type (spectral_waveform_kernel) :: skernel
  integer, dimension(:), pointer :: jf
  real :: df_skernel
  ! time kernel
  real :: t0,dt
  integer :: nwin,iwin,njt,ncomp,jcomp
  integer, dimension(:), allocatable :: nt1,nt2,jt
  character(len=character_length_component), dimension(:), allocatable :: comp

  ! filtering
  integer :: nfreq_measured
  real :: df_measured
  integer, dimension(:), pointer :: map_jf_2_filter_index
  complex, dimension(:), pointer :: event_filter,station_comp_filter
  complex, dimension(:,:), allocatable :: filter

  logical :: use_dmspace,compute_one_path,&
       stop_after_command_line,check_nt1_nt2,check_ipath
  integer :: npath,ipath1,ipath2,ipath,lu,ios,j

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
   ! process command line
   call new(cl,12,1,'h evid stname dmspce ipath1 ipath2 t0 dt nwin nt1 nt2 comp',&
        '0 1 1 1 1 1 1 1 1 1 1 1',printhelp)
   parfile = clManarg(cl,1)
!
   use_dmspace = clOptset(cl,4)
   compute_one_path = clOptset(cl,2) .and. clOptset(cl,3)
!
   stop_after_command_line = .false.
!
   if(.not. (use_dmspace.or.compute_one_path)) then
      print *, "use either one path or data model space file"
      stop_after_command_line = .true.
   end if
   if(use_dmspace .and. compute_one_path) then
      print *, "use either one path or data model space file"
      stop_after_command_line = .true.
   end if
   if(compute_one_path .and. (clOptset(cl,5).or.clOptset(cl,6))) then
      print *, "-ipath1, -ipath2 only to be used together with -dmspce"
      stop_after_command_line = .true.
   end if
!
   if(stop_after_command_line) then
      print *, ""
      stop_after_command_line = .true.
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
   call init(invbasics,trim(parfile),get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) then; call print(errmsg); endif
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
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
      dmspace_file = trim(clOptarg(cl,4))
      call new(errmsg,myname)
      call createDataSamplesFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
           .ifreq.iterbasics,trim(dmspace_file),get(fuh),errmsg)
      call undo(fuh)
      call addTrace(errmsg,myname)
      !if (.level.errmsg /= 0) call print(errmsg)
      call print(errmsg)
      if (.level.errmsg == 2) goto 1
      call dealloc(errmsg)
      paths => getPathsDataModelSpaceInfo(dmspace)
      if(.not.associated(paths))then
         print *, "no paths returned by data model space info object"
         goto 1
      end if
      npath = size(paths,2)
      call dealloc(dmspace)
!
      check_ipath = .true.
      ! define ipath1
      if(.not.clOptset(cl,5)) then
         ipath1 = 1
      else
         string = clOptarg(cl,5)
         read(string,*,iostat=ios) ipath1
         if(ios/=0) then
            print *, "there was an error reading integer index path1 from '-ipath1' input string '"//trim(string)//"'"
            stop_after_command_line = .true.
            check_ipath = .false.
         else
            if(ipath1<1) then
               write(*,*) "path1 ( = ",ipath1,") read from '-ipath1' input string '", &
                    trim(string),"' must be at least 1"
               stop_after_command_line = .true.
               check_ipath = .false.
            end if
         end if ! ios/=0
      end if
      ! define ipath2
      if(.not.clOptset(cl,6)) then
         ipath2 = npath
      else
         string = clOptarg(cl,6)
         read(string,*,iostat=ios) ipath2
         if(ios/=0) then
            print *, "there was an error reading integer index path2 from '-ipath2' input string '"//trim(string)//"'"
            stop_after_command_line = .true.
            check_ipath = .false.
         else
            if(ipath2>npath) then
               write(*,*) "path2 ( = ",ipath2,") read from '-ipath2' input string '", &
                    trim(string),"' is greater than maximum number of paths (",npath,") contained in data-model-space info "//&
                    "defined by file '"//trim(dmspace_file)//"'"
               stop_after_command_line = .true.
               check_ipath = .false.
            end if
         end if ! ios/=0
      end if
      ! check if ipath1<=ipath2
      if(check_ipath) then
         if(ipath1>ipath2) then
            write(*,*) "path1 ( = ",ipath1,") is greater than path2 ( = ",ipath2,")"
            stop_after_command_line = .true.
         end if
      end if
!
      if(stop_after_command_line) then
         print *, ""
         call printhelp
         goto 1
      end if
!
   else ! use_dmspace
!
      stop_after_command_line = .false.
!
      ! previously checked:
      ! if program comes here:  compute_one_path == (clOptset(cl,2) .and. clOptset(cl,3)) == .true.
!
      npath = 1; ipath1 = 1; ipath2 = 1
      allocate(paths(2,npath))
!
      ! check if event ID is valid
      paths(1,1) = clOptarg(cl,2)
      errmsg = searchEventidSeismicEventList(.evlist.invbasics,paths(1,1))
      if(.level. errmsg/=0) then
         write(*,*) "event ID '"//trim(paths(1,1))//"' (input string of option -evid) is not contained in event list"
         stop_after_command_line = .true.
      end if
      call dealloc(errmsg)
      ! check if station name is valid
      paths(2,1) = clOptarg(cl,3)
      errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,paths(2,1))
      if(.level. errmsg/=0) then
         write(*,*) "station name '"//trim(paths(2,1))//"' (input string of option -stname) is not contained in station list"
         stop_after_command_line = .true.
      end if
      call dealloc(errmsg)
!
      if(stop_after_command_line) then
         print *, ""
         call printhelp
         goto 1
      end if
   end if ! use_dmspace
!
   ! now treat options regarding time kernel computation
   stop_after_command_line = .false.
!
   ! t0
   if(.not.clOptset(cl,7)) then
      t0 = 0.0         
   else
      string = clOptarg(cl,7)
      read(string,*,iostat=ios) t0
      if(ios/=0) then
         print *, "there was an error reading real valued global time shift from '-t0' "//&
              "input string '"//trim(string)//"'"
         stop_after_command_line = .true.
      end if
   end if
!
   ! dt
   if(.not.clOptset(cl,8)) then
      print *, "please indicate the time step '-dt'"
      stop_after_command_line = .true.
   else
      string = clOptarg(cl,8)
      read(string,*,iostat=ios) dt
      if(ios/=0) then
         print *, "there was an error reading real valued global time step from '-dt' input string '"&
              //trim(string)//"'"
         stop_after_command_line = .true.
      end if
   end if
!
   ! nwin
   if(.not.clOptset(cl,9)) then
      print *, "please indicate the number of time windows '-nwin'"
      stop_after_command_line = .true.
   else
      string = clOptarg(cl,9)
      read(string,*,iostat=ios) nwin
      if(ios/=0) then
         print *, "there was an error reading integer number of time windows from '-nwin' input string '"&
              //trim(string)//"'"
         stop_after_command_line = .true.
      else
         if(nwin<0) then
            print *, "integer number of time windows ",nwin," read from '-nwin' input string '"&
                 //trim(string)//"' must not be negative"
            stop_after_command_line = .true.
         elseif(nwin>0) then
            check_nt1_nt2 = .true.
            ! nt1
            if(.not.clOptset(cl,10)) then
               print *, "please indicate the numbers nt1 by option '-nt1'"
               stop_after_command_line = .true.
               check_nt1_nt2 = .false.
            else
               string = clOptarg(cl,10)
               allocate(nt1(nwin))
               read(string,*,iostat=ios) nt1
               if(ios/=0) then
                  print *, "there was an error reading ",nwin," integer values of nt1 from '-nt1' input string '"&
                       //trim(string)//"'"
                  stop_after_command_line = .true.
                  check_nt1_nt2 = .false.
               end if
            end if
            ! nt2
            if(.not.clOptset(cl,11)) then
               print *, "please indicate the numbers nt2 by option '-nt2'"
               stop_after_command_line = .true.
               check_nt1_nt2 = .false.
            else
               string = clOptarg(cl,11)
               allocate(nt2(nwin))
               read(string,*,iostat=ios) nt2
               if(ios/=0) then
                  print *, "there was an error reading ",nwin," integer values of nt2 from '-nt2' input string '"&
                       //trim(string)//"'"
                  stop_after_command_line = .true.
                  check_nt1_nt2 = .false.
               end if
            end if
            ! check if nt2 >= nt1 for all iwin
            if(check_nt1_nt2) then
               do iwin=1,nwin
                  if(nt1(iwin)>nt2(iwin)) then
                     print *, "for ",iwin,"'th time window (out of ",nwin,") nt1,nt2 do not fulfill nt1 <= nt2 :"//&
                          "  nt1,nt2 =",nt1(iwin),nt2(iwin)
                     stop_after_command_line = .true.
                  end if
               end do ! iwin
            end if ! check_nt1_nt2
!
         else ! if(nwin<0) elseif(nwin>0)
            ! in this case nwin == 0
            print *, "number of time windows read from '-nwin' input string '"&
                 //trim(string)//"' is zero, in which case nothing is to be computed"
            stop_after_command_line = .true.
         end if ! if(nwin<0) elseif(nwin>0)
      end if ! ios/=0
   end if ! .not.clOptset(cl,9) ! nwin
!
   ! comp
   if(clOptset(cl,12)) then
      string = clOptarg(cl,12)
      read(string,*,iostat=ios) ncomp
      if(ios/=0) then
         write(*,*) "could not read integer number of components as first word of '-comp' input string '"&
              //trim(string)//"'"
         stop_after_command_line = .true.
         goto 10
      end if
      if(ncomp.le.0) then
         write(*,*) "number of components ",ncomp," (first word of '-comp' input string '"&
              //trim(string)//"' must be positive"
         stop_after_command_line = .true.
         goto 10
      end if
      allocate(comp(ncomp))
      read(string,*,iostat=ios) ncomp,comp
      if(ios/=0) then
         write(*,*) "could not read ",ncomp," components from '-comp' input string '"//trim(string)//"'"
         stop_after_command_line = .true.
         goto 10
      end if
      do jcomp=1,ncomp
         if(.not.validComponent(comp(jcomp))) then
            print *, jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
                 all_valid_components//"'"
            stop_after_command_line = .true.
            goto 10
         end if ! .not.validComponent(comp(jcomp))
      end do ! jcomp
   else ! clOptset(cl,12)
      print *, "please indicate the components '-comp'"
      stop_after_command_line = .true.
   end if ! clOptset(cl,12)
!
10 if(stop_after_command_line) then
      print *, ""
      call printhelp
      goto 1
   end if
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
      print *, " WARNING, time windows overlap: there are ",size(jt)-njt,&
           " duplicate time indices (out of a total of ",size(jt),"), which are ignored "//&
           "(no duplicate kernel values written to time kernel file)"
   end if
!
!------------------------------------------------------------------------
!  compute kernels
!
   print *, "time '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//"' kernels for ",&
        ipath2-ipath1+1," paths will be computed now: "
   print *, " path index | event ID        | station name    |"
   print *, "------------+-----------------+-----------------+"
   do ipath = ipath1,ipath2
      write(*,"(i12,a,a15,a,a15,a)") ipath," | ","'"//trim(paths(1,ipath))//"'"," | ","'"//trim(paths(2,ipath))//"'"," |"
   end do ! ipath
   print *, ""
!
   ! base of all time kernel filenames
   tkernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
        'time_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
   ! base of all spectral kernel filenames
   skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
        'spectral_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
!
   nfreq_measured = (.inpar.invbasics).ival.'MEASURED_DATA_NUMBER_OF_FREQ'
   df_measured = (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'
   nullify(event_filter,station_comp_filter,map_jf_2_filter_index)
!
   ! loop on all required paths
   do ipath = ipath1,ipath2
      ! use a new error message for each path
      call new(errmsg,myname)
!
      evid = paths(1,ipath)
      staname = paths(2,ipath)
!
      write(*,*) "# COMPUTE TIME KERNEL for ",ipath,"'th path: event '",trim(evid),"', station '",trim(staname),"'"
!
      call initialReadSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg)
      if(.level.errmsg == 2) then
         call print(errmsg)
         goto 1
      end if
!
      ! check if frequency discretization of this spectral kernel file is consistent with
      ! frequency discretization of the filters:
      jf => .jf.skernel 
      if(.not.associated(jf)) then
         write(*,*) "ERROR: no frequencies contained in kernel file"
         goto 1
      end if
      df_skernel = .df.skernel
      if( abs(df_measured-df_skernel)/df_measured > 1.e-4 ) then
         write(*,*) "ERROR: frequency step of spectral kernel ( = ",df_skernel,&
              ") differs from frequency step of measured data / filters ( = ",df_measured,&
              ") by more than 0.01 percent"
         goto 1
      end if
      if(associated(map_jf_2_filter_index)) deallocate(map_jf_2_filter_index)
      map_jf_2_filter_index => mapIfreq2ArrayIndex(invbasics,jf)
      if(.not.associated(map_jf_2_filter_index)) then
         write(*,*) "ERROR: some frequency indices contained in kernel file are not contained "//&
              "in frequency indices of measured data / filters"
         write(*,*) "jf_kernel = ",jf
         write(*,*) "jf_measured = ",.ifreq.invbasics
         goto 1
      end if
!
      ! now read in filter values for this event and all requested components of this station
      ! and select required filter values for this kernel
!
      if(associated(event_filter)) deallocate(event_filter)
      errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_EVENT_FILTER')//"filter_"//trim(evid),get(fuh),&
           event_filter,ndata=nfreq_measured)
      call undo(fuh)
      call addTrace(errmsg2,myname)
      if(.level.errmsg2 /= 0) call print(errmsg2)
      if(.level.errmsg2 == 2) goto 1
      call dealloc(errmsg2)
!
      allocate(filter(.nfreq.skernel,ncomp)) ! assuming .nfreq.skernel == size(jf) (which I can assume, this is the rules of the module)
      do jcomp = 1,ncomp         
         if(associated(station_comp_filter)) deallocate(station_comp_filter)
         errmsg2 = readAsciiData(trim((.inpar.invbasics).sval.'PATH_STATION_FILTER')//"filter_"//trim(staname)//&
              "_"//trim(comp(jcomp)),get(fuh),station_comp_filter,ndata=nfreq_measured)
         call undo(fuh)
         call addTrace(errmsg2,myname)
         if(.level.errmsg2 /= 0) call print(errmsg2)
         if(.level.errmsg2 == 2) goto 1
         call dealloc(errmsg2)
!
         filter(:,jcomp) = event_filter(map_jf_2_filter_index)*station_comp_filter(map_jf_2_filter_index)
      end do ! jcomp
!
      call transformSpectralToTimeWaveformKernel(skernel,comp,.comptrans.invbasics,staname,dt,jt(1:njt),&
           trim(tkernel_filebase)//trim(evid)//'_'//trim(staname),get(fuh),errmsg,filter,t0)
      call undo(fuh)
      if(.level.errmsg == 2) then
         call print(errmsg)
         goto 1
      end if
!
      call finalReadSpectralWaveformKernel(skernel,lu)
      call add(fuh,lu)
!
      if(.level.errmsg /= 0) call print(errmsg)
      if(.level.errmsg == 2) goto 1
      call dealloc(errmsg)
   end do ! ipath
!------------------------------------------------------------------------
!  clean up
!
1  if(associated(paths)) deallocate(paths)
   if(allocated(nt1)) deallocate(nt1)
   if(allocated(nt2)) deallocate(nt2)
   if(allocated(comp)) deallocate(comp)
   if(allocated(jt)) deallocate(jt)
   if(allocated(filter)) deallocate(filter)
   if(associated(map_jf_2_filter_index)) deallocate(map_jf_2_filter_index)
   if(associated(event_filter)) deallocate(event_filter)
   if(associated(station_comp_filter)) deallocate(station_comp_filter)
   call dealloc(errmsg); call dealloc(errmsg2)
   call dealloc(skernel)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
end program spec2timeKernels
!
!-----------------------------------------------------------------------------------------------------------------
!
 subroutine printhelp
   use componentTransformation, only: all_valid_components
   print '(50(1h-))'
   print *,'    spec2timeKernels [-h] [-evid event_id -stname station_name] [-dmspce dmspace_file [-ipath1 path1]'
   print *,'       [-ipath2 path2]] [-t0 tzero] -dt time_step -nwin ntime_windows -nt1 nt1_string -nt2 nt2_string'
   print *,'       -comp "ncomp comp" main_parfile'
   print *,''
   print *,'Description:'
   print *,'Compute time_kernel files from existing spectral_kernel files (by inverse Fourier transform) for the'
   print *,'time windows as specified by options t0,dt,nwin,nt1,nt2:'
   print *,'The total set of time samples consists of nwin time windows. nt1,nt2 are strings containing nwin integers'
   print *,'defining the start (nt1) and end (nt2) of each time window by an index. Times compute as t = t0 + jt*dt'
   print *,'whereby nt1 <= jt <= nt2.'
   print *,''
   print *,'THERE ARE TWO POSSIBLE WAYS TO DEFINE A SET OF KERNELS THAT ARE TRANSFORMED:'
   print *,'  (way 1) compute kernel for only one path, defined by eventID and station name using options -evid and -stname'
   print *,''
   print *,'  (way 2) use flag -dmspace in connection with optional range definition of the path index (flags -ipath1 -ipath2)'
   print *,'          in order to define a subset of the paths contained in the given data-model-space description.'
   print *,'          If the upper (lower) limit of the path range is not defined (i.e. -ipath1 (-ipath2) not set), the maximum'
   print *,'          (minimum) possible value is used.'
   print *,'          Then, all kernels for the defined range of paths in the given data-model-space description are computed.'
   print *,''
   print *,'Arguments:'
   print *,''
   print *,"    parfile: main parameter file of inversion"
   print *,''
   print *,'Options:'
   print *,''
   print *,'-h     : print help'
   print *,''
   print *,'SET OF KERNELS TO BE TRANSFORMED'
   print *,'  (way 1)'
   print *,'    -evid event_id : defines the event id of the one path (must belong to an event in main event list)'
   print *,'    -stname station_name   : defines the station name of the one path (must belong to a station in '//&
        'main station list)'
   print *,'  (way 2)'
   print *,"    -dmspce  : data model space input file to define a set of paths"
   print *,'    -ipath1 path1  : path1 is the first index of the path loop. By default, path1 = 1'
   print *,'    -ipath2 path2  : path2 is the last index of the path loop. By default, path2 = max_number_of_paths as to '//&
        'data-model-space description'
   print *,''
   print *,'TIME WAVEFORM KERNELS'
   print *,'-t0 tzero  : optional global time shift which is added to all times defined by dt,nt1,nt2 (default t0=0)'
   print *,'-dt time_step  : global time step of time discretization'
   print *,'-nwin ntime_windows  : integer number of time windows'
   print *,'-nt1 nt1_string  : string containing ntime_windows space separated integers defining nt1 for each time window'
   print *,'-nt2 nt2_string  : string containing ntime_windows space separated integers defining nt2 for each time window'
  print *,"-comp : ncomp, number of components follwing; comp, components (one of '"//trim(all_valid_components)//"')"
   print '(50(1h-))'
   return
 end subroutine printhelp
