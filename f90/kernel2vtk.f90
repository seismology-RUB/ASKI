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
program kernel2vtk
  use inversionBasics
  use iterationStepBasics
  use modelParametrization
  use componentTransformation
  use spectralWaveformKernel
  use invgridVtkFile
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=132) :: parfile,string

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=10) :: myname = 'kernel2vtk'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (spectral_waveform_kernel) :: kernel

  type (invgrid_vtk_file), dimension(:,:), allocatable :: vtkFile

  character(len=character_length_evid) :: evid
  character(len=character_length_staname) :: staname

  integer :: ios,lu,istat
  integer :: nparam,jparam,nfreq,jfreq,ncomp,jcomp
  character(len=character_length_param), dimension(:), allocatable :: param
  character(len=character_length_param) :: one_param
  character(len=character_length_pmtrz) :: parametrization
  integer, dimension(:), pointer :: all_ifreq
  integer, dimension(:), allocatable :: ifreq
  character(len=character_length_component), dimension(:), allocatable :: comp

  character(len=400) :: kernel_file,vtkFile_base,vtkFile_title,vtkFile_data_name
  complex, dimension(:,:), pointer :: k
  double precision, dimension(:,:), pointer :: trans_coef
  integer, dimension(:), pointer :: cells_filled

   real :: df

   external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
   ! process command line
   call new(cl,6,1,'h evid stname param ifreq comp','0 1 1 1 1 1',printhelp)
   parfile = clManarg(cl,1)
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,150)
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,trim(parfile),get(fuh),errmsg)
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
         stop
      end if
      call dealloc(errmsg)
   else
      print *, "please indicate '-evid'"
      print *, ""
      call printhelp
      stop
   end if
   ! handle -stname
   if(clOptset(cl,3)) then
      staname = clOptarg(cl,3)
      errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname,istat=istat)
      if(.level. errmsg/=0) then
         write(*,*) "station name '"//trim(staname)//"' (input string of option -stname) is not contained in station list"
         print *, ""
         call printhelp
         stop
      end if
      call dealloc(errmsg)
   else
      print *, "please indicate '-stname'"
      print *, ""
      call printhelp
      stop
   end if
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
            stop
         end if
         if(nparam.le.0) then
            write(*,*) "number of parameters ",nparam," (first word of '-param' input string '"&
                 //trim(string)//"' must be positive"
            print *, ""
            call printhelp
            stop
         end if
         allocate(param(nparam))
         read(string,*,iostat=ios) nparam,param
         if(ios/=0) then
            write(*,*) "could not read ",nparam," parameters from '-param' input string '"//trim(string)//"'"
            print *, ""
            call printhelp
            stop
         end if
         do jparam = 1,nparam
            if(.not.validParamModelParametrization(parametrization,param(jparam))) then
               write(*,*) jparam,"'th parameter '"//trim(param(jparam))//"' is not valid in parametrization '"//&
                    trim(parametrization)//"'. Valid Parametrizations (Parameters) are: "//&
                    all_valid_pmtrz_param
               print *, ""
               call printhelp
               stop
            end if
         end do
      end if ! trim(string) == 'all'
   else ! clOptset(cl,4)
      print *, "please indicate -param"
      print *, ""
      call printhelp
      stop
   end if ! clOptset(cl,4)
!
   ! handle -ifreq
   if(clOptset(cl,5)) then
      string = clOptarg(cl,5)
      if(trim(string) == 'all') then
         all_ifreq => .ifreq.iterbasics
         if(.not.associated(all_ifreq)) then
            write(*,*) "no iteration step specific frequencies defined"
            print *, ""
            call printhelp
            stop
         end if
         nfreq = size(all_ifreq)
         allocate(ifreq(nfreq))
         ifreq = all_ifreq
      else ! trim(string) == 'all'
         read(string,*,iostat=ios) nfreq
         if(ios/=0) then
            write(*,*) "could not read integer number of frequencies as first word of '-ifreq' input string '"&
                 //trim(string)//"'"
            print *, ""
            call printhelp
            stop
         end if
         if(nfreq.le.0) then
            write(*,*) "number of frequencies ",nfreq," (first word of '-ifreq' input string '"&
                 //trim(string)//"' must be positive"
            print *, ""
            call printhelp
            stop
         end if
         allocate(ifreq(nfreq))
         read(string,*,iostat=ios) nfreq,ifreq
         if(ios/=0) then
            write(*,*) "could not read ",nfreq," frequency indices from '-ifreq' input string '"//trim(string)//"'"
            print *, ""
            call printhelp
            stop
         end if
         all_ifreq => .ifreq.iterbasics
         do jfreq=1,nfreq
            if(.not.any(all_ifreq==ifreq(jfreq))) then
               write(*,*) jfreq,"'th frequency index ",ifreq(jfreq),&
                    " is not an interation step specific frequency as defined by iteration step specific parameter file"
               print *, ""
               call printhelp
               stop
            end if
         end do
      end if ! trim(string) == 'all'
   else ! clOptset(cl,5)
      print *, "please indicate '-ifreq'"
      print *, ""
      call printhelp
      stop
   end if ! clOptset(cl,5)
!
   ! handle -comp
   if(clOptset(cl,6)) then
      string = clOptarg(cl,6)
      read(string,*,iostat=ios) ncomp
      if(ios/=0) then
         write(*,*) "could not read integer number of components as first word of '-comp' input string '"&
              //trim(string)//"'"
         print *, ""
         call printhelp
         stop
      end if
      if(ncomp.le.0) then
         write(*,*) "number of components ",ncomp," (first word of '-comp' input string '"&
              //trim(string)//"' must be positive"
         print *, ""
         call printhelp
         stop
      end if
      allocate(comp(ncomp))
      read(string,*,iostat=ios) ncomp,comp
      if(ios/=0) then
         write(*,*) "could not read ",ncomp," components from '-comp' input string '"//trim(string)//"'"
         print *, ""
         call printhelp
         stop
      end if
      do jcomp=1,ncomp
         if(.not.validComponent(comp(jcomp))) then
            print *, jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
                 all_valid_components//"'"
            print *, ""
            call printhelp
            stop
         end if ! .not.validComponent(comp(jcomp))
      end do ! comp
      trans_coef => transform(.comptrans.invbasics,(/'CX','CY','CZ'/),comp,istat)
      if(.not.associated(trans_coef)) then
         print *, "the computation of the transformation coefficients to transform CX,CY,CZ to ",comp," was not successful"
         print *, ""
         stop
      end if ! .not.associated(trans_coef)
   else ! clOptset(cl,6)
      print *, "please indicate '-comp'"
      print *, ""
      call printhelp
      stop
   end if ! clOptset(cl,6)
!------------------------------------------------------------------------
!  program starting here
!
   kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
        'spectral_kernel_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

   call new(errmsg,myname)
   call initialReadSpectralWaveformKernel(kernel,kernel_file,get(fuh),errmsg)
   if (.level.errmsg /= 0) then; call print(errmsg); endif
   if (.level.errmsg == 2) stop
   call dealloc(errmsg)

   if(.pmtrz.kernel /= parametrization) then
      write(*,*) "parametrization '"//trim(.pmtrz.kernel)//"' of kernel file differs from parametrization '"//&
           trim(parametrization)//"' currently set in main parameter file"
      stop
   end if

   df = .df.kernel
   if( df<0. .or. abs(((.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP') -df)/df > 1.e-4 ) then
      write(*,*) "frequency step of kernel ",df," differs from frequency step of measured data ",&
           (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'," by more than 0.01 percent"
      stop
   end if

   print *,"kernel2vtk: converting '"//trim(parametrization)//"'-kernel to vtk for"
   print *,"  parameters ","'"//param//"',"
   print *,"  frequency indices (corresponding to df = ",df,")  ",ifreq
   print *,"  components ","'"//comp//"',"
   print *,""
   
   allocate(vtkFile(ncomp,nparam))

   cells_filled => getFilledCells(.intw.iterbasics)

   ! loop on frequencies
   do jfreq = 1,nfreq
      call new(errmsg,myname)
      call readSpectralWaveformKernel(kernel,ifreq(jfreq),errmsg)
      if (.level.errmsg /= 0) then; call print(errmsg); endif
      if (.level.errmsg == 2) stop
      call dealloc(errmsg)

      ! handle all parameters
      do jparam=1,nparam
         k => getValuesSpectralWaveformKernel(kernel,param(jparam))
         if(.not.associated(k)) then
            write(*,*) "no '"//trim(param(jparam))//"' sensitivity values contained in kernel"
            stop
         end if

         ! loop on all components
         do jcomp=1,ncomp
            if(jfreq==1) then
               ! initiate vtk file
               write(vtkFile_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
               write(vtkFile_title,*) trim(comp(jcomp)),"-component of spectral ",trim(param(jparam)),&
                    '-'//trim(parametrization)//' Kernel at frequency ',ifreq(jfreq)*df,' Hz on inversion grid'
               call new(errmsg,myname)
               call init(vtkfile(jcomp,jparam),.invgrid.iterbasics,trim(vtkFile_base),&
                    trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtkFile_title),&
                    cell_indx_req=cells_filled)
               if (.level.errmsg /= 0) then; call print(errmsg); endif
               if (.level.errmsg == 2) stop
               call dealloc(errmsg)
               print *,"kernel2vtk: creating vtk files with basename '"//trim(vtkFile_base)//"'"
            end if ! jfreq==1
            ! write kernel values to vtk file
            write(vtkFile_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'
            call new(errmsg,myname)
            call writeData(vtkfile(jcomp,jparam),get(fuh),matmul(k(cells_filled,:),real(trans_coef(jcomp,:))),&
                 errmsg,data_name=trim(vtkFile_data_name),file_index=ifreq(jfreq))
            call undo(fuh)
            if (.level.errmsg /= 0) then; call print(errmsg); endif
            if (.level.errmsg == 2) stop
            call dealloc(errmsg)

            print *,"kernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                 "-"//trim(parametrization)//" Kernel, frequency index ",ifreq(jfreq),", frequency ",ifreq(jfreq)*df," Hz"
         end do ! jcomp

      end do ! jparam
   end do ! jfreq

   call finalReadSpectralWaveformKernel(kernel,lu)
   call add(fuh,lu)
   call dealloc(kernel)
! 
!------------------------------------------------------------------------
!  clean up
!
   deallocate(param,ifreq,comp)
   do jparam=1,nparam
      do jcomp=1,ncomp
         call dealloc(vtkfile(jcomp,jparam))
      end do ! jcomp
   end do ! jparam
   if(associated(cells_filled)) deallocate(cells_filled)
   if(associated(trans_coef)) deallocate(trans_coef)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
   print *, "kernel2vtk: good bye"
end program kernel2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  use componentTransformation
  print '(50(1h-))'
  print *,'                   kernel2vtk [-h] -evid event_id -stname station_name -param "nparam param" '//&
       '-ifreq ["nfreq ifreq","all"] -comp "ncomp comp"  main_parfile'
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
  print *,"-ifreq     : nfreq, number of frequency indices following; ifreq, frequency indices "//&
       "(e.g. '4  3 4 6 10','2  1 8')"
  print *,''
  print *,"-comp : ncomp, number of components follwing; comp, components (one of '"//&
       trim(all_valid_components)//"')"
  print '(50(1h-))'
  return
end subroutine printhelp
