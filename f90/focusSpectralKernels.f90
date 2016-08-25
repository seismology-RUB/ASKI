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
program focusSpectralKernels
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use kernelFocus
  use invgridVtkFile
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,basename
  character(len=max_length_string), dimension(:), pointer :: str_vec
  integer :: nparam_foc
  character(len=character_length_param), dimension(:), allocatable :: param_foc
  ! real :: fintense ! option for simple Backus Gilbert focussing (not used below)

  type (file_unit_handler) :: fuh
  integer :: lu1,lu2

  type (error_message) :: errmsg
  character(len=19) :: myname = 'focusSpectraKernels'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace_full,mspace_focus

  ! kernel focus stuff
  type (kernel_focus) :: kfocus
  real, dimension(:), pointer :: km
  real, dimension(:,:), pointer :: kf
  real, dimension(:,:), pointer :: a
  integer :: nfoc,ifoc

  ! vtk file stuff
  type (invgrid_vtk_file) :: invgrid_vtk
  character(len=300) :: vtk_filename
  character(len=100) :: vtk_title,vtk_data_name
  character(len=character_length_param) :: param_name
  character(len=character_length_param), dimension(:), pointer :: param
  integer, dimension(:), pointer :: indx,cell

  integer :: icell,ncell_focus,iparam
  real :: c1,c2,c3
  integer, dimension(:), allocatable :: icell_focus
  character(len=character_length_pmtrz) :: parametrization

  logical :: stop_after_command_line

  real :: c1_foc_min, c1_foc_max, c2_foc_min, c2_foc_max, c3_foc_min, c3_foc_max

  nullify(str_vec,km,kf,a,param,indx,cell)

!------------------------------------------------------------------------
!  preliminary processing
!
  stop_after_command_line = .false.
!
  call init(ap,myname,"compute Backus-Gilbert focussing of sensitivity kernels on a defined focussing region in "//&
       "the model space")
  call addPosarg(ap,"dmspace_file","sval","data model space input file which defines the rows and columns of "//&
       "the kernel matrix")
  call addPosarg(ap,"outfile_base","sval","base name of output files RELATIVE to iteration step output "//&
       "directory (used for all files, with suitable extensions)")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-c1min",.true.,"(mandatory) minimum first coordinate of focussing subvolume (Cartesian "//&
       "X, or spherical lat in deg (-90<=lat<=90))","rval","")
  call addOption(ap,"-c1max",.true.,"(mandatory) maximum first coordinate of focussing subvolume (Cartesian "//&
       "X, or spherical lat in deg (-90<=lat<=90))","rval","")
  call addOption(ap,"-c2min",.true.,"(mandatory) minimum second coordinate of focussing subvolume "//&
       "(Cartesian Y or spherical lon in deg (0<=lon<=360))","rval","")
  call addOption(ap,"-c2max",.true.,"(mandatory) maximum second coordinate of focussing subvolume "//&
       "(Cartesian Y or spherical lon in deg (0<=lon<=360))","rval","")
  call addOption(ap,"-c3min",.true.,"(mandatory) minimum third coordinate of focussing subvolume "//&
       "(Cartesian Z or depth in km)","rval","")
  call addOption(ap,"-c3max",.true.,"(mandatory) maximum third coordinate of focussing subvolume "//&
       "(Cartesian Z or depth in km)","rval","")
  call addOption(ap,"-param",.true.,"(mandatory) list of model parameters which will be focussed on","svec","")
  ! call addOption(ap,"-ffactor",.true.,"(optional) a factor controlling the intensity of the focussing","rval","1.0") ! option for simple Backus Gilbert focussing (not used below)
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  if(.not. ((ap.optset.'-c1min').and.(ap.optset.'-c1max').and.&
       (ap.optset.'-c2min').and.(ap.optset.'-c2max').and.&
       (ap.optset.'-c3min').and.(ap.optset.'-c3max'))) then
     write(*,*) "ERROR: please indicate all of -c1min, -c1max, -c2min, -c2max, -c3min, -c3max"
     stop_after_command_line = .true.
  end if
  if(.not.(ap.optset.'-param')) then
     write(*,*) "ERROR: please indicate -param"
     stop_after_command_line = .true.
  end if
!
  if(stop_after_command_line) then
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  dmspace_file = ap.sval.'dmspace_file'
  basename = ap.sval.'outfile_base'
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  ! handle -c(123)_(min,max)
  c1_foc_min = ap.rval.'-c1min'
  c1_foc_max = ap.rval.'-c1max'
  c2_foc_min = ap.rval.'-c2min'
  c2_foc_max = ap.rval.'-c2max'
  c3_foc_min = ap.rval.'-c3min'
  c3_foc_max = ap.rval.'-c3max'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  ! handle -param
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
  nparam_foc = size(str_vec)
  allocate(param_foc(nparam_foc))
  do iparam = 1,nparam_foc
     param_foc(iparam) = str_vec(iparam)
  end do ! iparam
  deallocate(str_vec)
!
  ! ! handle -ffactor ! option for simple Backus Gilbert focussing (not used below)
  ! fintense = ap.rval.'-ffactor'
  ! if (.level.(.errmsg.ap) == 2) then
  !    call print(.errmsg.ap)
  !    call usage(ap)
  !    goto 1
  ! end if
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,150)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(main_parfile),get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if(.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if(.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  basename = trim(.iterpath.iterbasics)//trim((.inpar.iterbasics).sval.'PATH_OUTPUT_FILES')//trim(basename)
!
  parametrization = sval(.inpar.invbasics,'MODEL_PARAMETRIZATION')
  ! check if all requested focussing model parameters are valid
  do iparam = 1,nparam_foc
     if(.not.validParamModelParametrization(parametrization,param_foc(iparam))) then
        write(*,*) "focusSpectraKernels ### ERROR : ",iparam,"'th incoming focussing model parameter '",&
             trim(param_foc(iparam)),"' is not valid in parametrization '",trim(parametrization),"'"
        goto 1
     end if
  end do ! iparam

!
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace_full,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,parametrization,&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if(.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  actual program starts
!
  write(*,*) "focusSpectraKernels ### WELCOME TO KERNEL FOCUSSING"
  write(*,*) ""
  write(*,*) "focusSpectraKernels ### complete data and model space defined by file '"//trim(dmspace_file)//"'"
  write(*,*) "focusSpectraKernels ### containing ",.ndata.dmspace_full," data samples and ",.nmval.dmspace_full," model values"
  write(*,*) ""
  write(*,*) "focusSpectraKernels ### try to focus sensitivity for the model parameters ","'"//param_foc//"'"," in the subvolume"
  write(*,*) c1_foc_min," < C1 < ",c1_foc_max
  write(*,*) c2_foc_min," < C2 < ",c2_foc_max
  write(*,*) c3_foc_min," < C3 < ",c3_foc_max
  write(*,*) ""
  write(*,*) "focusSpectraKernels ### writing output files with base filename '"//trim(basename)//"'"
  write(*,*) ""
!
!------------------------------------------------------------------------
!  define focussing area
!
  ncell_focus = 0
  allocate(icell_focus(.ncell.(.invgrid.iterbasics)))
!
  do icell = 1,.ncell.(.invgrid.iterbasics)
     ! get cell center of current cell
     call new(errmsg,myname)
     call getCenterCellInversionGrid(.invgrid.iterbasics,icell,c1,c2,c3,errmsg,coords_type='event')
     if(.level.errmsg /= 0) call print(errmsg)
     if(.level.errmsg == 2) goto 1
     call dealloc(errmsg)
!
     if(c1 .lt. c1_foc_max .and. c1 .gt. c1_foc_min .and. &
          c2 .lt. c2_foc_max .and. c2 .gt. c2_foc_min .and. &
          c3 .lt. c3_foc_max .and. c3 .gt. c3_foc_min) then
        ncell_focus = ncell_focus + 1
        icell_focus(ncell_focus) = icell
     end if
  end do ! icell
!
  if(ncell_focus == 0) then
     write(*,*) "focusSpectraKernels ### ERROR: no inversion grid cells were found in focussing subvolume"
     goto 1
  end if
!
  ! create focussing model subspace
  call new(errmsg,myname)
  call setParametrizationDataModelSpaceInfo(mspace_focus,parametrization,errmsg)
  call addModelValuesDataModelSpaceInfo(mspace_focus,param_foc,icell_focus(1:ncell_focus),all_combinations=.true.)
  write(*,*) "focusSpectraKernels ### sucessfully found ",.nmval.mspace_focus," model values to focus on:"
  write(*,*) ncell_focus," inversion grid cells for each of the ",nparam_foc," parameters : ","'"//param_foc//"'; "
  ! write focussing subvolume as sub-inversion-grid to vtk file
  ! initiate vtk file
  vtk_filename = trim(basename)//'_fvol'
  vtk_title = 'focussing subvolume of inversion grid'
  call new(errmsg,myname)
  call init(invgrid_vtk,.invgrid.iterbasics,vtk_filename,sval(.inpar.invbasics,'DEFAULT_VTK_FILE_FORMAT'),&
       errmsg,vtk_title,cell_indx_req=icell_focus(1:ncell_focus))
  if(.level.errmsg == 2) then; call print(errmsg); goto 1; endif
  ! write values to file
  write(*,*) "focusSpectraKernels ### writing focussing subvolume as vtk file"
  write(*,*) ""
  call writeInvgrid(invgrid_vtk,get(fuh),errmsg,overwrite=.true.)
  call undo(fuh)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(invgrid_vtk)
  call dealloc(errmsg)
!
  ! initiate kernel focus object
  write(*,*) "focusSpectraKernels ### initiating kernel focus object"
  call new(errmsg,myname)
  lu1 = get(fuh)
  lu2 = get(fuh)
  call initiateKernelFocus(kfocus,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
       ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
       .ncell.(.invgrid.iterbasics),.pcorr.invbasics,lu1,lu2,dmspace_full,errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'))
  call add(fuh,lu1); call add(fuh,lu2)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  ! ! compute simple Backus Gilbert focussing
  ! write(*,*) "focusSpectraKernels ### compute simple Backus Gilbert focussing"
  ! call new(errmsg,myname)
  ! call computeSimpleBackusGilbertKernelFocus(kfocus,dmspace_full,mspace_focus,errmsg,fintense=fintense)
  ! if(.level.errmsg /= 0) call print(errmsg)
  ! if(.level.errmsg == 2) goto 1
  ! call dealloc(errmsg)
  ! write(*,*) ""
!
  ! compute original Backus Gilbert focussing
  write(*,*) "focusSpectraKernels ### compute original Backus Gilbert focussing"
  call new(errmsg,myname)
  call computeOriginalBackusGilbertKernelFocus(kfocus,dmspace_full,mspace_focus,.invgrid.iterbasics,errmsg)
  !call print(errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  ! get focussing coefficients
  a => .a.kfocus
!
  ! get mean sensitivity
  call new(errmsg,myname)
  call getMeanKernelFocus(kfocus,km,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  ! get focussed sensitivity
  nfoc = .nfoc.kfocus
  call new(errmsg,myname)
  call getKernelFocus(kfocus,kf,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  if(nfoc /= size(kf,2)) then
     write(*,*) "focusSpectraKernels ### ERROR: number of returned focussings ",size(kf,2),&
          " differs from the expected number ",nfoc
     goto 1
  end if
!
!
!------------------------------------------------------------------------
!  write results to files with basename
!
  ! write focussing coefficients
  write(*,*) "focusSpectraKernels ### writing focussing coefficients to file"
  errmsg = writeAsciiData(trim(basename)//'_coef.txt',get(fuh),a)
  call undo(fuh)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  nullify(param,indx,cell)
!
  ! for each parameter in data model space, write mean and focussed sensitivity as vtk files on inversion grid
  do while(nextParamModelParametrization(parametrization,param_name))
     if(associated(param)) deallocate(param)
     if(associated(cell)) deallocate(cell)
     if(associated(indx)) deallocate(indx)
     allocate(param(1)); param(1) = param_name
!
     indx => getIndxModelValues(dmspace_full,param,cell)
     if(.not.associated(indx)) then
        write(*,*) "focusSpectraKernels ### no kernel values for model parameter '"//trim(param_name)//"' in full model space"
        cycle
     end if
!
     ! WRITE MEAN KERNEL VALUES TO FILE
!
     call new(errmsg,myname)
!
     ! initiate vtk file
     vtk_filename = trim(basename)//'_kmean_'//trim(parametrization)//'-'//trim(param_name)
     vtk_title = trim(parametrization)//'-'//trim(param_name)//' mean kernel values on inversion grid'
     call init(invgrid_vtk,.invgrid.iterbasics,vtk_filename,sval(.inpar.invbasics,'DEFAULT_VTK_FILE_FORMAT'),&
          errmsg,vtk_title,cell_indx_req=cell)
     if(.level.errmsg == 2) then; call print(errmsg); goto 1; endif
!
     ! write values to file
     vtk_data_name = trim(parametrization)//'-'//trim(param_name)//'_kmean'
     write(*,*) "focusSpectraKernels ### writing mean kernel values of parameter '"//trim(param_name)//&
          "' to file '"//trim(vtk_filename)//"'.vtk"
     call writeData(invgrid_vtk,get(fuh),km(indx),errmsg,data_name=trim(vtk_data_name),overwrite=.true.)
     call undo(fuh)
     if(.level.errmsg /= 0) call print(errmsg)
     if(.level.errmsg == 2) goto 1
     call dealloc(invgrid_vtk)
     call dealloc(errmsg)
!
     ! WRITE FOCUSSED KERNEL VALUES TO FILE
!
     call new(errmsg,myname)
!
     ! initiate vtk file
     vtk_filename = trim(basename)//'_kfoc_'//trim(parametrization)//'-'//trim(param_name)
     vtk_title = trim(parametrization)//'-'//trim(param_name)//' focussed kernel values on inversion grid'
     call init(invgrid_vtk,.invgrid.iterbasics,vtk_filename,sval(.inpar.invbasics,'DEFAULT_VTK_FILE_FORMAT'),&
          errmsg,vtk_title,cell_indx_req=cell)
     if(.level.errmsg == 2) then; call print(errmsg); goto 1; endif
!
     write(*,*) "focusSpectraKernels ### writing focussed kernel values of parameter '"//trim(param_name)//"'"
     ! write values to file
     do ifoc = 1,nfoc
        write(*,*) "focusSpectraKernels ###    focal volume ",ifoc
        vtk_data_name = trim(parametrization)//'-'//trim(param_name)//'_kfoc'
        call writeData(invgrid_vtk,get(fuh),kf(indx,ifoc),errmsg,data_name=trim(vtk_data_name),file_index=ifoc,overwrite=.true.)
        call undo(fuh)
        if(.level.errmsg == 2) then; call print(errmsg); goto 1; endif
     end do ! ifoc
     call dealloc(invgrid_vtk)
     if(.level.errmsg /= 0) call print(errmsg)
     call dealloc(errmsg)
  end do ! while (nextParam)
  write(*,*) ""
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "focusSpectraKernels ### good bye"
!
1  call dealloc(kfocus)
   call dealloc(dmspace_full); call dealloc(mspace_focus)
!
   if(associated(km)) deallocate(km)
   if(associated(kf)) deallocate(kf)
!
   if(associated(param)) deallocate(param)
   if(associated(indx)) deallocate(indx)
   if(associated(cell)) deallocate(cell)
!
   if(associated(str_vec)) deallocate(str_vec)
   if(allocated(param_foc)) deallocate(param_foc)
   if(allocated(icell_focus)) deallocate(icell_focus)
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
   call dealloc(errmsg)
end program focusSpectralKernels
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'                   focusSpectralKernels [-h] [-ffac focus_intensity] -param "nparam param" '//&
!        'dmspace_file outfile_base parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    dmspace_file: data model space input file which defines data and model space"
!   print *,"    outfile_base: base name of output files relative to output directory "//&
!        "(will be used for all files, with suitable extensions)"
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print *,"-param : nparam, number of parameters following; param, parameters "//&
!        "(e.g. '2 vp vs', or '3 rho lambda mu')"
!   print *,''
!   print *,"-ffac : focus_intensity is a factor controlling the intensity of the focussing (default = 1.0)"
!   print *,''
!   print '(50(1h-))'
!   return
! end subroutine printhelp
