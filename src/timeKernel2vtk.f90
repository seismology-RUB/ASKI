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
program timeKernel2vtk
  use inversionBasics
  use iterationStepBasics
  use modelParametrization
  use timeWaveformKernel
  use invgridVtkFile
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=132) :: parfile,string

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=14) :: myname = 'timeKernel2vtk'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (time_waveform_kernel) :: kernel

  type (invgrid_vtk_file), dimension(:,:), allocatable :: vtkFile

  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  integer :: ios,lu
  integer :: nparam,jparam,ncomp,jcomp,nwin,iwin,njt,ijt
  character(len=character_length_param), dimension(:), allocatable :: param
  character(len=character_length_param) :: one_param
  character(len=character_length_pmtrz) :: parametrization
  integer, dimension(:), allocatable :: nt1,nt2,jt
  character(len=character_length_component), dimension(:), allocatable :: comp
  character(len=character_length_component), dimension(:), pointer :: comp_kernel

  character(len=400) :: kernel_file,vtkFile_base,vtkFile_title,vtkFile_data_name
  real, dimension(:), pointer :: k
  integer, dimension(:), pointer :: cells_filled

  real :: dt,t0

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
   ! process command line
   call new(cl,8,1,'h evid stname param nwin nt1 nt2 comp','0 1 1 1 1 1 1 1',printhelp)
   parfile = clManarg(cl,1)
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,150)
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,parfile,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) stop
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) stop
   call dealloc(errmsg)
!
   parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
!------------------------------------------------------------------------
!  processing of command line
!
   ! handle -evid
   if(clOptset(cl,2)) then
      evid = clOptarg(cl,2)
      errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid)
      if(.level. errmsg/=0) then
         write(*,*) "event ID '"//trim(evid)//"' (input string of option -evid) is not contained in event list"
         print *, ""
         call printhelp
         goto 1
      end if
      call dealloc(errmsg)
   else
      print *, "please indicate '-evid'"
      print *, ""
      call printhelp
      goto 1
   end if
!
   ! handle -stname
   if(clOptset(cl,3)) then
      staname = clOptarg(cl,3)
      errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname)
      if(.level. errmsg/=0) then
         write(*,*) "station name '"//trim(staname)//"' (input string of option -stname) is not contained in station list"
         print *, ""
         call printhelp
         goto 1
      end if
      call dealloc(errmsg)
   else
      print *, "please indicate '-stname'"
      print *, ""
      call printhelp
      goto 1
   end if
!
   ! handle -param
   if(clOptset(cl,4)) then
      string = clOptarg(cl,4)
      if(trim(string) == 'all') then
         nparam = numberOfParamModelParametrization(parametrization)
         allocate(param(nparam))
         do while(nextParamModelParametrization(parametrization,one_param))
            param(jparam) = one_param
         end do
      else ! trim(string) == 'all'
         read(string,*,iostat=ios) nparam
         if(ios/=0) then
            write(*,*) "could not read integer number of parameters as first word of '-param' input string '"&
                 //trim(string)//"'"
            print *, ""
            call printhelp
            goto 1
         end if
         if(nparam.le.0) then
            write(*,*) "number of parameters ",nparam," (first word of '-param' input string '"&
                 //trim(string)//"' must be positive"
            print *, ""
            call printhelp
            goto 1
         end if
         allocate(param(nparam))
         read(string,*,iostat=ios) nparam,param
         if(ios/=0) then
            write(*,*) "could not read ",nparam," parameters from '-param' input string '"//trim(string)//"'"
            print *, ""
            call printhelp
            goto 1
         end if
         do jparam = 1,nparam
            if(.not.validParamModelParametrization(parametrization,param(jparam))) then
               write(*,*) jparam,"'th parameter '"//trim(param(jparam))//"' is not valid in parametrization '"//&
                    trim(parametrization)//"'. Valid Parametrizations (Parameters) are: "//&
                    all_valid_pmtrz_param
               print *, ""
               call printhelp
               goto 1
            end if
         end do
      end if ! trim(string) == 'all'
   else ! clOptset(cl,4)
      print *, "please indicate -param"
      print *, ""
      call printhelp
      goto 1
   end if ! clOptset(cl,4)
!
   ! nwin
   if(.not.clOptset(cl,5)) then
      print *, "please indicate the number of time windows '-nwin'"
      print *, ""
      call printhelp
      goto 1
   else
      string = clOptarg(cl,5)
      read(string,*,iostat=ios) nwin
      if(ios/=0) then
         print *, "there was an error reading integer number of time windows from '-nwin' input string '"&
              //trim(string)//"'"
         print *, ""
         call printhelp
         goto 1
      end if
      if(nwin<0) then
         print *, "integer number of time windows ",nwin," read from '-nwin' input string '"&
              //trim(string)//"' must not be negative"
         print *, ""
         call printhelp
         goto 1
      elseif(nwin>0) then
         ! nt1
         if(.not.clOptset(cl,6)) then
            print *, "please indicate the numbers nt1 by option '-nt1'"
            print *, ""
            call printhelp
            goto 1
         else
            string = clOptarg(cl,6)
            allocate(nt1(nwin))
            read(string,*,iostat=ios) nt1
            if(ios/=0) then
               print *, "there was an error reading ",nwin," integer values of nt1 from '-nt1' input string '"&
                    //trim(string)//"'"
               print *, ""
               call printhelp
               goto 1
            end if
         end if
         ! nt2
         if(.not.clOptset(cl,7)) then
            print *, "please indicate the numbers nt2 by option '-nt2'"
            print *, ""
            call printhelp
            goto 1
         else
            string = clOptarg(cl,7)
            allocate(nt2(nwin))
            read(string,*,iostat=ios) nt2
            if(ios/=0) then
               print *, "there was an error reading ",nwin," integer values of nt2 from '-nt2' input string '"&
                    //trim(string)//"'"
               print *, ""
               call printhelp
               goto 1
            end if
         end if
         ! check if nt2 >= nt1 for all iwin
         do iwin=1,nwin
            if(nt1(iwin)>nt2(iwin)) then
               print *, "for ",iwin,"'th time window (out of ",nwin,") nt1,nt2 do not fulfill nt1 <= nt2 :"//&
                    "  nt1,nt2 =",nt1(iwin),nt2(iwin)
               print *, ""
               call printhelp
               goto 1
            end if
         end do ! iwin
!
      else ! if(nwin<0) elseif(nwin>0)
         ! in this case nwin == 0
         print *, "number of time windows read from '-nwin' input string '"&
              //trim(string)//"' is zero, in which case nothing is to be computed"
         print *, ""
         call printhelp
         goto 1
      end if ! if(nwin<0) elseif(nwin>0)
   end if ! .not.clOptset(cl,5) ! nwin
!
   ! handle -comp
   if(clOptset(cl,8)) then
      string = clOptarg(cl,8)
      read(string,*,iostat=ios) ncomp
      if(ios/=0) then
         write(*,*) "could not read integer number of components as first word of '-comp' input string '"&
              //trim(string)//"'"
         print *, ""
         call printhelp
         goto 1
      end if
      if(ncomp.le.0) then
         write(*,*) "number of components ",ncomp," (first word of '-comp' input string '"&
              //trim(string)//"' must be positive"
         print *, ""
         call printhelp
         goto 1
      end if
      allocate(comp(ncomp))
      read(string,*,iostat=ios) ncomp,comp
      if(ios/=0) then
         write(*,*) "could not read ",ncomp," components from '-comp' input string '"//trim(string)//"'"
         print *, ""
         call printhelp
         goto 1
      end if
      do jcomp=1,ncomp
         if(.not.validComponent(comp(jcomp))) then
            print *, jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
                 all_valid_components//"'"
            print *, ""
            call printhelp
            goto 1
         end if ! .not.validComponent(comp(jcomp))
      end do ! comp
   else ! clOptset(cl,8)
      print *, "please indicate '-comp'"
      print *, ""
      call printhelp
      goto 1
   end if ! clOptset(cl,8)
!
   ! now define vector of all time indices (all time windows)
   ! check whether time windows to overlap or not, remove duplicate indices
   njt = sum( (nt2-nt1) + 1)
   allocate(jt(njt))
!
   ! put all indices of first time window into vector jt
   njt = nt2(1)-nt1(1) + 1
   jt(1:njt) = (/ (ijt,ijt=nt1(1),nt2(1)) /)
   ! afterwards loop on the rest of the time windows, only adding time indices which are not yet present in jt
   do iwin = 2,nwin
      do ijt = nt1(iwin),nt2(iwin)
         if(.not.any(jt(1:njt)==ijt)) then
            njt = njt + 1
            jt(njt) = ijt
         end if
      end do ! ijt
   end do ! iwin
   ! now, value njt indicates sensible values in array jt!
   if(njt < size(jt)) then
      print *, " WARNING, time windows overlap: there are ",size(jt)-njt,&
           " duplicate time indices (out of a total of ",size(jt),"), which are ignored "//&
           "(no duplicate kernel vtk files will be produced)"
   end if
!
!------------------------------------------------------------------------
!  program starting here
!
   kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
        'time_kernel_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

   print *,"timeKernel2vtk: kernel file '"//trim(kernel_file)//"'"

   call new(errmsg,myname)
   call initialReadTimeWaveformKernel(kernel,kernel_file,get(fuh),errmsg)
   if (.level.errmsg /= 0) then; call print(errmsg); endif
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)

   if(.pmtrz.kernel /= parametrization) then
      write(*,*) "ERROR in timeKernel2vtk: parametrization '"//trim(.pmtrz.kernel)//&
           "' of kernel file differs from parametrization '"//&
           trim(parametrization)//"' currently set in main parameter file"
      goto 1
   end if

   dt = .dt.kernel
   t0 = .tzero.kernel
   comp_kernel => .comp.kernel
   if(.not.associated(comp_kernel)) then
      write(*,*) "ERROR in timeKernel2vtk: no components contained in kernel file"
      goto 1
   end if
   do jcomp = 1,ncomp
      if(.not.any(comp_kernel == comp(jcomp))) then
         write(*,*) "ERROR in timeKernel2vtk: requested component '"//trim(comp(jcomp))//&
              "' is not contained in kernel file; components contained in kernel file: '"//(comp_kernel)//"'"
         goto 1         
      end if
   end do ! jcomp

   print *,"contains"
   print *,"   parametrization '"//trim(parametrization)//"'"
   print *,"   nt, dt, t0 = ",.nt.kernel,dt,t0
   print *,"   components = '"//.comp.kernel//"'"
   print *,"will be converted to vtk files for"
   print *,"   parameters ","'"//param//"',"
   print *,"   components ","'"//comp//"',"
   print *,"   time indices ",jt(1:njt)
   print *,""

   allocate(vtkFile(ncomp,nparam))

   cells_filled => getFilledCells(.intw.iterbasics)

   ! loop on time samples
   do ijt = 1,njt
      call new(errmsg,myname)
      call readTimeWaveformKernel(kernel,jt(ijt),errmsg)
      if (.level.errmsg /= 0) then; call print(errmsg); endif
      if (.level.errmsg == 2) goto 1
      call dealloc(errmsg)

      ! handle all parameters
      do jparam=1,nparam

         ! loop on all components
         do jcomp=1,ncomp
            k => getValuesTimeWaveformKernel(kernel,comp(jcomp),param(jparam))
            if(.not.associated(k)) then
               write(*,*) "no '"//trim(param(jparam))//"' sensitivity values for component '"//trim(comp(jcomp))//&
                    "'contained in kernel"
               goto 1
            end if

            if(ijt==1) then
               ! initiate vtk file
               write(vtkFile_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
               write(vtkFile_title,*) trim(comp(jcomp)),"-component of time ",trim(param(jparam)),&
                    '-'//trim(parametrization)//' Kernel at time ',t0+jt(ijt)*dt,' s on inversion grid'
               call new(errmsg,myname)
               call init(vtkfile(jcomp,jparam),.invgrid.iterbasics,trim(vtkFile_base),&
                    trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtkFile_title),&
                    cell_indx_req=cells_filled)
               if (.level.errmsg /= 0) then; call print(errmsg); endif
               if (.level.errmsg == 2) goto 1
               call dealloc(errmsg)
               print *,"timeKernel2vtk: creating vtk files with basename '"//trim(vtkFile_base)//"'"
            end if ! ijt==1
            ! write kernel values to vtk file
            write(vtkFile_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'
            call new(errmsg,myname)
            call writeData(vtkfile(jcomp,jparam),get(fuh),k(cells_filled),&
                 errmsg,data_name=trim(vtkFile_data_name),file_index=jt(ijt))!,overwrite=.true.)
            call undo(fuh)
            if (.level.errmsg /= 0) then; call print(errmsg); endif
            if (.level.errmsg == 2) goto 1
            call dealloc(errmsg)

            print *,"timeKernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                 "-"//trim(parametrization)//" Kernel, time index ",jt(ijt),", time ",t0+jt(ijt)*dt," s"
         end do ! jcomp

      end do ! jparam
   end do ! ijt

   call finalReadTimeWaveformKernel(kernel,lu)
   call add(fuh,lu)
   call dealloc(kernel)
! 
!------------------------------------------------------------------------
!  clean up
!
1  if(allocated(vtkfile)) then
      do jparam=1,nparam
         do jcomp=1,ncomp
            call dealloc(vtkFile(jcomp,jparam))
         end do ! jcomp
      end do ! jparam
      deallocate(vtkFile)
   end if
   if(allocated(param)) deallocate(param)
   if(allocated(comp)) deallocate(comp)
   if(allocated(nt1)) deallocate(nt1)
   if(allocated(nt2)) deallocate(nt2)
   if(allocated(jt)) deallocate(jt)
   if(associated(cells_filled)) deallocate(cells_filled)
   call dealloc(kernel)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
   print *, "timeKernel2vtk: good bye"
end program timeKernel2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  use componentTransformation
  print '(50(1h-))'
  print *,'      timeKernel2vtk [-h] -evid event_id -stname station_name -param "nparam param" -nwin ntime_windows '
  print *,'          -nt1 nt1_string -nt2 nt2_string -comp "ncomp comp"  main_parfile'
  print *,''
  print *,"    main_parfile: main parameter file of inversion"
  print *,''
  print *,'Options:'
  print *,''
  print *,'-h     : print help'
  print *,''
   print *,'-evid event_id : defines the event id of the one path (must belong to an event in main event list)'
   print *,''
   print *,'-stname station_name   : defines the station name of the one path (must belong to a station in main station '//&
        'list)'
  print *,''
  print *,"-param : nparam, number of parameters following; param, parameters "//&
       "(e.g. '2 vp vs', or '3 rho lambda mu')"
  print *,''
   print *,'-nwin ntime_windows  : integer number of time windows'
  print *,''
   print *,'-nt1 nt1_string  : string containing ntime_windows space separated integers defining nt1 (start time index) '//&
        'for each time window'
  print *,''
   print *,'-nt2 nt2_string  : string containing ntime_windows space separated integers defining nt2 (end time index) '//&
        'for each time window'
  print *,''
  print *,"-comp : ncomp, number of components follwing; comp, components (one of '"//&
       trim(all_valid_components)//"')"
  print *,''
  print '(50(1h-))'
  return
end subroutine printhelp
