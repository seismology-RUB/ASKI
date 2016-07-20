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
program computeKernels
  use inversionBasics
  use iterationStepBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use dataModelSpaceInfo
  use kernelDisplacement
  use kernelGreenTensor
  use spectralWaveformKernel
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=300) :: parfile,string,dmspace_file,&
       skernel_filebase,kd_filebase,kgt_filebase

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=14) :: myname = 'computeKernels'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics
  
  type (data_model_space_info) :: dmspace
  character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  type (kernel_displacement) :: kd
  type (kernel_green_tensor) :: kgt

  type (spectral_waveform_kernel) :: skernel

  logical :: use_dmspace,compute_one_path
  integer :: npath,ipath1,ipath2,ipath,jf,lu,ios

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
   ! process command line
   ! here, more options to indicate a set of paths could be added
   call new(cl,6,1,'h evid stname dmspce ipath1 ipath2','0 1 1 1 1 1',printhelp)
   parfile = clManarg(cl,1)
!
   use_dmspace = clOptset(cl,4)
   compute_one_path = clOptset(cl,2) .and. clOptset(cl,3)
!
   if(.not. (use_dmspace.or.compute_one_path)) then
      print *, "use either one path or data model space file"
      call printhelp
      stop
   end if
!
   if(use_dmspace .and. compute_one_path) then
      print *, "use either one path or data model space file"
      call printhelp
      stop
   end if
!
   if(compute_one_path .and. (clOptset(cl,5).or.clOptset(cl,6))) then
      print *, "-ipath1, -ipath2 only to be used together with -dmspce"
      call printhelp
      stop
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
      ! define ipath1
      if(.not.clOptset(cl,5)) then
         ipath1 = 1
      else
         string = clOptarg(cl,5)
         read(string,*,iostat=ios) ipath1
         if(ios/=0) then
            print *, "there was an error reading integer index path1 from '-ipath1' input string '"//trim(string)//"'"
            print *, ""
            call printhelp
            goto 1
         end if
         if(ipath1<1) then
            write(*,*) "path1 ( = ",ipath1,") read from '-ipath1' input string '", &
                 trim(string),"' must be at least 1"
            print *, ""
            call printhelp
            goto 1         
         end if
      end if
      ! define ipath2
      if(.not.clOptset(cl,6)) then
         ipath2 = npath
      else
         string = clOptarg(cl,6)
         read(string,*,iostat=ios) ipath2
         if(ios/=0) then
            print *, "there was an error reading integer index path2 from '-ipath2' input string '"//trim(string)//"'"
            print *, ""
            call printhelp
            goto 1
         end if
         if(ipath2>npath) then
            write(*,*) "path2 ( = ",ipath2,") read from '-ipath2' input string '", &
                 trim(string),"' is greater than maximum number of paths (",npath,") contained in data-model-space info "//&
                 "defined by file '"//trim(dmspace_file)//"'"
            print *, ""
            call printhelp
            goto 1         
         end if
      end if
      ! check if ipath1<=ipath2
      if(ipath1>ipath2) then
         write(*,*) "path1 ( = ",ipath1,") is greater than path2 ( = ",ipath2,")"
         print *, ""
         call printhelp
         goto 1
      end if
!
   else ! use_dmspace
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
         print *, ""
         call printhelp
         goto 1
      end if
      call dealloc(errmsg)
      ! check if station name is valid
      paths(2,1) = clOptarg(cl,3)
      errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,paths(2,1))
      if(.level. errmsg/=0) then
         write(*,*) "station name '"//trim(paths(2,1))//"' (input string of option -stname) is not contained in station list"
         print *, ""
         call printhelp
         goto 1
      end if
      call dealloc(errmsg)
!
   end if ! use_dmspace
!------------------------------------------------------------------------
!  compute kernels
!
   print *, "spectral '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//"' kernels for ",&
        ipath2-ipath1+1," paths will be computed now: "
   print *, " path index | event ID        | station name    |"
   print *, "------------+-----------------+-----------------+"
   do ipath = ipath1,ipath2
      write(*,"(i12,a,a15,a,a15,a)") ipath," | ","'"//trim(paths(1,ipath))//"'"," | ","'"//trim(paths(2,ipath))//"'"," |"
   end do ! ipath
   print *, ""
!
   ! base of all kernel displacement filenames
   kd_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//'kernel_displ_'
   ! base of all kernel green tensor filenames
   kgt_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_GREEN_TENSORS')//'kernel_gt_'
   ! base of all spectral kernel filenames
   skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
        'spectral_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
!
   ! use variables evid / staname to memorize the event ID / station name for which current kd / kgt objects are initiated
   ! indicate by setting evid / staname to '' that no kd,kgt objects have been initiated yet
   evid = ''; staname = ''
   ! loop on all required paths
   do ipath = ipath1,ipath2
!
      write(*,*) "# COMPUTE SPECTRAL KERNEL for ",ipath,"'th path: event '",trim(paths(1,ipath)),"', station '",trim(paths(2,ipath)),"'"
!
      ! check if for this path the previous kd or kgt files can be used (i.e. left open). Otherwise, close old ones and open new ones
      ! kd object still valid?
      if(evid /= paths(1,ipath)) then
         ! close old file if open
         if(evid/='') call dealloc(kd,fuh)
         evid = paths(1,ipath)
         ! initiate kernel displacement
         call new(errmsg,myname)
         call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
              trim(kd_filebase)//trim(evid),errmsg)
         if (.level.errmsg /= 0) then; call print(errmsg); endif
         if (.level.errmsg == 2) goto 2
         call dealloc(errmsg)
      end if ! evid /= paths(1,ipath)
      ! kgt object still valid?
      if(staname /= paths(2,ipath)) then
         ! close old file if open
         if(staname/='') call dealloc(kgt,fuh)
         staname = paths(2,ipath)
         ! initiate kernel Green tensor
         call new(errmsg,myname)
         call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
              trim(kgt_filebase)//trim(staname),errmsg)
         if (.level.errmsg /= 0) then; call print(errmsg); endif
         if (.level.errmsg == 2) goto 2
         call dealloc(errmsg)
      end if ! staname /= paths(2,ipath)
!
      ! initiate and inital write spectral waveform kernel
      call new(errmsg,myname)
      call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION', &
           .ncell.(.invgrid.iterbasics),errmsg)
      if (.level.errmsg /= 0) then; call print(errmsg); endif
      if (.level.errmsg == 2) goto 2
      call dealloc(errmsg)

      call new(errmsg,myname)
      call initialWriteSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),&
           get(fuh),errmsg,nfreq=.nf.kd)
      if (.level.errmsg /= 0) then; call print(errmsg); endif
      if (.level.errmsg == 2) goto 2
      call dealloc(errmsg)
!
      ! now loop on all frequencies
      do while (nextFrequencyKernelDisplacement(kd,jf))
         write(*,*) "#     frequency index ",jf," , frequency ",(.df.kd)*jf," [Hz]"
         ! read current frequency of kernel displacement
         call new(errmsg,myname)
         call readFrequencyKernelDisplacement(kd,jf,errmsg)
         if (.level.errmsg /= 0) then; call print(errmsg); endif
         if (.level.errmsg == 2) goto 2
         call dealloc(errmsg)
!
         ! read current frequency of kernel Green tensor
         call new(errmsg,myname)
         call readFrequencyKernelGreenTensor(kgt,jf,errmsg)
         if (.level.errmsg /= 0) then; call print(errmsg); endif
         if (.level.errmsg == 2) goto 2
         call dealloc(errmsg)
!                  
         ! compute spectral waveform kernel
         call new(errmsg,myname)
         call computeSpectralWaveformKernel(skernel,kd,kgt,.krm.iterbasics,.intw.iterbasics,errmsg)
         if (.level.errmsg /= 0) then; call print(errmsg); endif
         if (.level.errmsg == 2) goto 2
         call dealloc(errmsg)
!
         ! write spectral waveform kernel
         call writeSpectralWaveformKernel(skernel)
      end do ! while (nextFrequencyKernelDisplacement(kd,jf))
!
      ! final write and deallocate velocity kernel
      call new(errmsg,myname)
      call finalWriteSpectralWaveformKernel(skernel,errmsg,lu)
      if (.level.errmsg /= 0) then; call print(errmsg); endif
      if (.level.errmsg == 2) goto 2
      call dealloc(errmsg)
      call add(fuh,lu)
      call dealloc(skernel)
!
   end do ! ipath
!
2  if(evid/='') call dealloc(kgt,fuh)
   if(staname/='') call dealloc(kd,fuh)
!------------------------------------------------------------------------
!  clean up
!
1  if(associated(paths)) deallocate(paths)
   call dealloc(skernel)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
 end program computeKernels
!
!-----------------------------------------------------------------------------------------------------------------
!
 subroutine printhelp
   print '(50(1h-))'
   print *,'    computeKernels [-h] [-evid event_id] [-stname station_name] [-dmspce dmspace_file] [-ipath1 path1] '//&
        '[-ipath2 path2] main_parfile'
   print *,''
   print *,'THERE ARE TWO POSSIBLE WAYS TO DEFINE A SET OF SPECTRAL KERNELS THAT ARE COMPUTED:'
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
   print *,'-evid event_id : defines the event id of the one path (must belong to an event in main event list) (way 1)'
   print *,''
   print *,'-stname station_name   : defines the station name of the one path (must belong to a station in main station '//&
        'list) (way 1)'
   print *,''
   print *,"-dmspce  : data model space input file to define a set of paths (way 2)"
   print *,''
   print *,'-ipath1 path1  : path1 is the first index of the path loop. By default, path1 = 1 (way 2)'
   print *,''
   print *,'-ipath2 path2  : path2 is the last index of the path loop. By default, path2 = max_number_of_paths as to '//&
        'data-model-space description (way 2)'
   print '(50(1h-))'
   return
 end subroutine printhelp
