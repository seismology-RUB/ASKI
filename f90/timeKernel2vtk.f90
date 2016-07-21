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
program timeKernel2vtk
  use inversionBasics
  use iterationStepBasics
  use modelParametrization
  use timeWaveformKernel
  use invgridVtkFile
  use wpVtkFile
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,str
   character(len=max_length_string), dimension(:), pointer :: str_vec

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=14) :: myname = 'timeKernel2vtk'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (time_waveform_kernel) :: kernel

   type (invgrid_vtk_file), dimension(:,:), allocatable :: ig_vtk
   type (wp_vtk_file), dimension(:,:), allocatable :: wp_vtk

  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  integer :: lu
  integer :: nparam,jparam,ncomp,jcomp,nwin,iwin,njt,ijt,ntot_kernel
  character(len=character_length_param), dimension(:), allocatable :: param
  character(len=character_length_pmtrz) :: parametrization
  integer, dimension(:), pointer :: nt1,nt2
  integer, dimension(:), allocatable :: jt
  character(len=character_length_component), dimension(:), allocatable :: comp

  character(len=400) :: kernel_file,vtkFile_base,vtkFile_title,vtkFile_data_name
  real, dimension(:), pointer :: k
  integer, dimension(:), pointer :: cells_filled,wp_inside

  real :: dt,t0

  logical :: kernel_on_wp

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Produces vtk files from one EXISTING binary time sensitivity kernel file for a selection "//&
       "of time steps defined by vectors of starting and end indices -nt1 , -nt2 ")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-evid',.true.,"(mandatory) defines the event id of the one path (must belong to an event in"//&
       " main event list)",'sval','')
  call addOption(ap,'-staname',.true.,"(mandatory) defines the station name of the one path (must belong to a "//&
        "station in main station list",'sval','')
  call addOption(ap,'-comp',.true.,"(mandatory) vector of receiver components for which the kernel should be "//&
       "computed; valid components are: "//all_valid_components,'svec','')
  call addOption(ap,'-param',.true.,"(mandatory) vector of parameter names for which the kernel should be "//&
       "computed",'svec','')
  call addOption(ap,"-nt1",.true.,"(mandatory) vector of length nwin, giving starting indices of nwin time "//&
       "windows (must have same length as vector given by -nt2)","ivec","")
  call addOption(ap,"-nt2",.true.,"(mandatory) vector of length nwin, giving end indices of nwin time "//&
       "windows (must have same length as vector given by -nt1)","ivec","")
  call addOption(ap,"-wp",.false.,"if set, time waveform kernel files '...ON-WP...' (kernels on wavefield "//&
       "points) will be read (assuming they exist). If not set, normal time kernel files (pre-integrated) are "//&
       "read.")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  ! handle -comp
  if(ap.optset.'-comp') then
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
     ncomp = size(str_vec)
     allocate(comp(ncomp))
     do jcomp = 1,ncomp
        comp(jcomp) = str_vec(jcomp)
     end do
     deallocate(str_vec)
     if(.not.allValidComponents(comp,i_invalid=jcomp)) then
        write(*,*) "ERROR: ",jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
             all_valid_components//"'"
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
  end if
!
   ! handle -wp
   kernel_on_wp = ap.optset.'-wp'
!
   main_parfile = ap.sval.'main_parfile'
!
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,150)
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,main_parfile,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
!------------------------------------------------------------------------
!  processing of the rest of the arguments given by argument parser
!
   ! check if event ID is valid
   str = ap.sval.'-evid'
   evid = str
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid)
   if(.level. errmsg/=0) then
      write(*,*) "ERROR: event ID '"//trim(evid)//"' (input string of option -evid) is not contained in event list"
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
   call dealloc(errmsg)
!
   ! check if station name is valid
   str = ap.sval.'-staname'
   staname = str
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname)
   if(.level. errmsg/=0) then
      write(*,*) "ERROR: station name '"//trim(staname)//"' (input string of option -stname) is not contained in "//&
           "station list"
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
   call dealloc(errmsg)
!
   ! handle -param
   if(ap.optset.'-param') then
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
      nparam = size(str_vec)
      allocate(param(nparam))
      do jparam = 1,nparam
         param(jparam) = str_vec(jparam)
         if(.not.validParamModelParametrization(parametrization,param(jparam))) then
            write(*,*) jparam,"'th given model parameter '",trim(param(jparam)),"' not valid. Valid "//&
                 "parameterizations(parameters are ",all_valid_pmtrz_param
            write(*,*) ""
            call usage(ap)
            goto 1
         end if
      end do ! jparam
      deallocate(str_vec)
   end if
!
   ! nt1
   if(.not.(ap.optset.'-nt1')) then
      write(*,*) "ERROR: please indicate starting indices of nwin>0 time windows by option '-nt1'"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      nt1 => ap.ivec.'-nt1'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
   end if
   ! nt2
   if(.not.(ap.optset.'-nt2')) then
      write(*,*) "ERROR: please indicate end indices of nwin>0 time windows by option '-nt2'"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      nt2 => ap.ivec.'-nt2'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
   end if
   ! check if nt2 >= nt1 for all time windows
   nwin = size(nt1)
   if(size(nt2)/=nwin .or. nwin <= 0) then
      write(*,*) "ERROR: size of vector -nt1 = ",nwin," and size of vector -nt2 =",size(nt2)," must be equal and > 0"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      do iwin=1,nwin
         if(nt1(iwin)>nt2(iwin)) then
            write(*,*) "ERROR: for ",iwin,"'th time window (out of ",nwin,") nt1,nt2 do not fulfill nt1 <= "//&
                 "nt2 :  nt1,nt2 =",nt1(iwin),nt2(iwin)
            write(*,*) ""
            call usage(ap)
            goto 1
         end if
      end do ! iwin
   end if ! size(nt2)/=nwin .or. nwin <= 0
!
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
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
      write(*,*) " WARNING, time windows overlap: there are ",size(jt)-njt,&
           " duplicate time indices (out of a total of ",size(jt),"), which are ignored "//&
           "(no duplicate kernel vtk files will be produced)"
   end if
!
!------------------------------------------------------------------------
!  program starting here
!
   if(kernel_on_wp) then
      write(*,*)"timeKernel2vtk: converting plain '"//trim(parametrization)//"'-kernel values on wavefield points to vtk for"
      kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'time_kernel_ON-WP_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)
      ntot_kernel = .ntot.(.wp.iterbasics)
   else ! kernel_on_wp
      write(*,*)"timeKernel2vtk: converting pre-integrated '"//trim(parametrization)//"'-kernel values on inversion grid ",&
           "cells to vtk for"
      kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
           'time_kernel_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)
      ntot_kernel = .ncell.(.invgrid.iterbasics)
   end if ! kernel_on_wp
   write(*,*)"  event '",trim(evid),"'"
   write(*,*)"  station '",trim(staname),"' ; components ","'"//comp//"',"
   write(*,*)"  parameters ","'"//param//"',"
   write(*,*) ""
   write(*,*)"kernel file '"//trim(kernel_file)//"'"

   call new(errmsg,myname)
   call initiateTimeWaveformKernel(kernel,parametrization,param,ntot_kernel,comp,errmsg,kernel_on_wp=kernel_on_wp)
   if (.level.errmsg == 2) then
      call print(errmsg)
      goto 1
   end if
   call initialReadTimeWaveformKernel(kernel,kernel_file,get(fuh),errmsg)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)

   dt = .dt.kernel
   t0 = .tzero.kernel

   write(*,*)"contains:"
   write(*,*)"   parametrization '"//trim(parametrization)//"'"
   write(*,*)"   nt, dt, t0 = ",.nt.kernel,dt,t0
   write(*,*)"will be converted to vtk files for ",nwin," time windows and a total of ",njt," time samples: "
   write(*,*)"   time indices jt = ",jt(1:njt)
   write(*,*)"   i.e. times t = t0 + jt*dt = ",t0+jt(1:njt)*dt
   write(*,*)""

   if(kernel_on_wp) then
      allocate(wp_vtk(nparam,ncomp))
      wp_inside => getWpInside(.intw.iterbasics)
      if(.not.associated(wp_inside)) then
         write(*,*) "ERROR: there are wavefield points which are located inside the inversion grid"
         goto 1
      end if
   else
      allocate(ig_vtk(nparam,ncomp))
      cells_filled => getFilledCells(.intw.iterbasics)
      if(.not.associated(cells_filled)) then
         write(*,*) "ERROR: there are only empty inversion grid cells (or no cells at all)"
         goto 1
      end if
   end if

   ! loop on time samples
   do ijt = 1,njt
      ! use one error message for each time samples (better debugging in case of an error)
      call new(errmsg,myname)

      call readTimeWaveformKernel(kernel,jt(ijt),errmsg)
      if (.level.errmsg /= 0) call print(errmsg)
      if (.level.errmsg == 2) goto 1

      ! handle all parameters
      do jcomp=1,ncomp

         ! loop on all components
         do jparam=1,nparam
            k => getValuesByParamCompTimeWaveformKernel(kernel,param(jparam),comp(jcomp))
            if(.not.associated(k)) then
               write(*,*) "no '"//trim(param(jparam))//"' sensitivity values for component '"//trim(comp(jcomp))//&
                    "'contained in kernel. THIS ERROR SHOULD NOT OCCUR!!"
               goto 1
            end if

            if(ijt==1) then
               ! initiate vtk file
               write(vtkFile_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
               write(vtkFile_title,*) trim(comp(jcomp)),"-component of time ",trim(param(jparam)),&
                    '-'//trim(parametrization)//' Kernel at time ',t0+jt(ijt)*dt,' s on inversion grid'
               call new(errmsg,myname)
               if(kernel_on_wp) then
                  call init(wp_vtk(jparam,jcomp),.wp.iterbasics,.invgrid.iterbasics,trim(vtkFile_base),&
                       trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtkFile_title),&
                       wp_indx_req=wp_inside)
               else
                  call init(ig_vtk(jparam,jcomp),.invgrid.iterbasics,trim(vtkFile_base),&
                       trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtkFile_title),&
                       cell_indx_req=cells_filled)
               end if
               if (.level.errmsg /= 0) call print(errmsg)
               if (.level.errmsg == 2) goto 1
               write(*,*)"timeKernel2vtk: creating vtk files with basename '"//trim(vtkFile_base)//"'"
            end if ! ijt==1
            ! write kernel values to vtk file
            write(vtkFile_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'
            call new(errmsg,myname)
            if(kernel_on_wp) then
               call writeData(wp_vtk(jparam,jcomp),get(fuh),k(wp_inside),&
                    errmsg,data_name=trim(vtkFile_data_name),file_index=jt(ijt))!,overwrite=.true.)
            else
               call writeData(ig_vtk(jparam,jcomp),get(fuh),k(cells_filled),&
                    errmsg,data_name=trim(vtkFile_data_name),file_index=jt(ijt))!,overwrite=.true.)
            end if
            call undo(fuh)
            if (.level.errmsg /= 0) call print(errmsg)
            if (.level.errmsg == 2) goto 1

            write(*,*)"timeKernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                 "-"//trim(parametrization)//" Kernel, time index ",jt(ijt),", time ",t0+jt(ijt)*dt," s"
         end do ! jparam

      end do ! jcomp
      call dealloc(errmsg)
   end do ! ijt

   call finalReadTimeWaveformKernel(kernel,lu)
   call add(fuh,lu)
   call dealloc(kernel)
! 
!------------------------------------------------------------------------
!  clean up
!
   write(*,*) "timeKernel2vtk: good bye"
!
1  if(allocated(ig_vtk)) then
      do jcomp=1,ncomp
         do jparam=1,nparam
            call dealloc(ig_vtk(jparam,jcomp))
         end do ! jparam
      end do ! jcomp
      deallocate(ig_vtk)
   end if
   if(allocated(wp_vtk)) then
      do jcomp=1,ncomp
         do jparam=1,nparam
            call dealloc(wp_vtk(jparam,jcomp))
         end do ! jparam
      end do ! jcomp
      deallocate(wp_vtk)
   end if
   if(allocated(param)) deallocate(param)
   if(allocated(comp)) deallocate(comp)
   if(associated(nt1)) deallocate(nt1)
   if(associated(nt2)) deallocate(nt2)
   if(allocated(jt)) deallocate(jt)
   if(associated(cells_filled)) deallocate(cells_filled)
   if(associated(wp_inside)) deallocate(wp_inside)
   call dealloc(errmsg)
   call dealloc(kernel)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
end program timeKernel2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   use componentTransformation
!   print '(50(1h-))'
!   print *,'Program timeKernel2vtk writes time sensitivity kernels as vtk files for specific paths, parameters,'
!   print *,'components and time steps. Can either handle pre-integrated kernel values on inversion grid cells, or'
!   print *,'original kernel values on wavefield points (flag -wp). '
!   print *,'IN ANY CASE binary time kernel files must be produced (for the requested time samples) by program'
!   print *,'spec2timeKernel BEFORE using timeKernel2vtk.'
!   print *,'Usage:'
!   print *,''
!   print *,'      timeKernel2vtk [-h] -evid event_id -stname station_name -comp "ncomp comp" -param "nparam param"'
!   print *,'          -nwin ntime_windows -nt1 nt1_string -nt2 nt2_string [-wp] main_parfile'
!   print *,''
!   print *,"    main_parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!    print *,'-evid event_id : defines the event id of the one path (must belong to an event in main event list)'
!    print *,''
!    print *,'-stname station_name   : defines the station name of the one path (must belong to a station in main station '//&
!         'list)'
!   print *,''
!   print *,"-comp : ncomp, number of components follwing; comp, components (one of '"//&
!        trim(all_valid_components)//"')"
!   print *,''
!   print *,"-param : nparam, number of parameters following; param, parameters "//&
!        "(e.g. '2 vp vs', or '3 rho lambda mu')"
!   print *,''
!    print *,'-nwin ntime_windows  : integer number of time windows'
!   print *,''
!    print *,'-nt1 nt1_string  : string containing ntime_windows space separated integers defining nt1 (start time index) '//&
!         'for each time window'
!   print *,''
!    print *,'-nt2 nt2_string  : string containing ntime_windows space separated integers defining nt2 (end time index) '//&
!         'for each time window'
!   print *,''
!   print *,"-wp     : if set, then the original kernel values on the WAVEFIELD POINTS are handled, which will be"
!   print *,"          read from time waveform kernel files '...ON-WP...' (assuming that they exist, i.e. they were"
!   print *,"          produced before by program spec2timeKernel)."
!   print *,"          Otherwise, pre-integrated time kernel values are transformed to vtk which are read from"
!   print *,"          the regular time kernel files (not containing 'ON-WP' in their filenames)."
!   print '(50(1h-))'
!   return
! end subroutine printhelp
