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
program solveKernelSystem
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use linearModelSmoothing
  use kernelLinearSystem
  use modelParametrization
  use kernelInvertedModel
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=300) :: parfile,dmspace_file,outfile,string

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=17) :: myname = 'solveKernelSystem'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  type (kernel_linear_system) :: KLSE
  type (linear_model_smoothing) :: lmsmooth
  logical :: add_smoothing,outfile_is_relative
  character(len=character_length_scaling_type) :: smoothing_scaling_type
  integer :: nscale_val
  real, dimension(:), allocatable :: smoothing_scaling_values

  real :: misfit
  character(len=100) :: sm_boundary_condition
  real, dimension(:,:), pointer :: dparam

  type (kernel_inverted_model) :: kim_up,kim_ref
  character(len=character_length_param) :: param_name
  character(len=character_length_param), dimension(:), pointer :: pparam
  integer, dimension(:), pointer :: idx_dmspace,pcell

  integer :: ios,lu_out,lu1,lu2
  logical :: outfile_exists

  external printhelp
!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call new(cl,6,3,'h smooth scltyp sclval smbnd odir','0 0 1 1 1 0',printhelp)
  parfile = clManarg(cl,3)
  dmspace_file = clManarg(cl,1)
  outfile = clManarg(cl,2)
!
  ! smooth
  add_smoothing = clOptset(cl,2)
!
  if((clOptset(cl,3).or.clOptset(cl,4).or.clOptset(cl,5)) .and. .not.add_smoothing) then
     write(*,*) "scaling_type, scaling_values, or boundary_conditions for smoothing must only be given when -smooth is set"
     call printhelp
     stop
  end if
!
  ! -scltyp
  smoothing_scaling_type = ''
  if(clOptset(cl,3)) smoothing_scaling_type = clOptarg(cl,3)
!
  ! -sclval
  nscale_val = -1
  if(clOptset(cl,4)) then
     string = clOptarg(cl,4)
     read(string,*,iostat=ios) nscale_val
     if(ios/=0) then
        write(*,*) "could not read integer number of scaling values for smoothing as first entry of "//&
             "-sclval input '"//trim(string)//"'"
        call printhelp
        stop
     end if
     if(nscale_val<1) then
        write(*,*) "number of scaling values (first entry of -sclval input '"//trim(string)//"' must be positive"
        call printhelp
        stop
     end if
     allocate(smoothing_scaling_values(nscale_val))
     read(string,*,iostat=ios) nscale_val,smoothing_scaling_values
     if(ios/=0) then
        write(*,*) "could not read ",nscale_val," real numbers (scaling values) for smoothing from "//&
             "-sclval input '"//trim(string)//"'"
        call printhelp
        stop
     end if
  end if
!  
  ! smbnd
! IGNORE, FOR NOW
!  if(clOptset(cl,5)) then
!     sm_boundary_condition = clOptarg(cl,5)
!  else
!     sm_boundary_condition = 'continuous'
!  end if
!
  outfile_is_relative = clOptset(cl,6)
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
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
  if(outfile_is_relative) outfile = trim(.iterpath.iterbasics)//&
       trim((.inpar.iterbasics).sval.'PATH_OUTPUT_FILES')//trim(outfile)
!------------------------------------------------------------------------
!  check if output files already exist
!
  write(*,*) "base name of output files will be '"//trim(outfile)//"'"
!
  ! check if output text file exists
  inquire(file=trim(outfile)//'_out.txt',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "output text file '"//trim(outfile)//"_out.txt' already exists. Please (re)move it."
     stop
  end if
  ! check if kernel inverted model update file exists
  inquire(file=trim(outfile)//'_up.kim',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "model uptdate output file '"//trim(outfile)//"_up.kim' already exists. Please (re)move it."
     stop
  end if
  ! check if kernel inverted model new file exists
  inquire(file=trim(outfile)//'_new.kim',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "new inverted model output file '"//trim(outfile)//"_new.kim' already exists. Please (re)move it."
     stop
  end if
  ! also check if vtk output files already exist ?! (-> very extensive, as we would have to assume here naming of vtk files as done by kernelInvertedModel module)
!------------------------------------------------------------------------
!  open output textfile to write
!
  lu_out = get(fuh)
  open(unit=lu_out,file=trim(outfile)//"_out.txt",status='unknown',form='formatted',action='write',iostat=ios)
  if(ios/=0) then
     write(*,*) "could not open output text file '"//trim(outfile)//"_out.txt' to write. Raised iostat = ",ios
     close(lu_out)
     stop
  end if
  write(lu_out,*) ""; write(lu_out,*) ""
  write(lu_out,*) "welcome to ASKI (put more general information about everything here)"; write(lu_out,*) ""
  write(lu_out,*) "inverting data now"; write(lu_out,*) ""
  write(lu_out,*) "base name of output files will be '"//trim(outfile)//"'"; write(lu_out,*) ""
!------------------------------------------------------------------------
!  setup data model space info object
!
  write(*,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
  write(lu_out,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
!
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,sval(.inpar.invbasics,'MODEL_PARAMETRIZATION'),&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!
  write(*,*) "there are ",.ndata.dmspace," data samples and ",.nparam.dmspace," model parameters"
  write(lu_out,*) "there are ",.ndata.dmspace," data samples and ",.nparam.dmspace," model parameters"; write(lu_out,*) ""
!------------------------------------------------------------------------
!  create smoothing conditions
!
  if(add_smoothing) then
!
     call new(errmsg,myname)
     if(smoothing_scaling_type=='' .and. nscale_val==-1) then
        write(*,*) "creating linear smoothing constraints, no scaling"
        write(lu_out,*) "creating linear smoothing constraints, no scaling"
        call createNeighbourAverageLinearModelSmoothing(lmsmooth,.invgrid.iterbasics,dmspace,errmsg)
     elseif(smoothing_scaling_type/='' .and. nscale_val==-1) then
        write(*,*) "creating linear smoothing constraints, with scaling '"//trim(smoothing_scaling_type)//"'"
        write(lu_out,*) "creating linear smoothing constraints, with scaling '"//trim(smoothing_scaling_type)//"'"
        call createNeighbourAverageLinearModelSmoothing(lmsmooth,.invgrid.iterbasics,dmspace,errmsg,&
             scaling_type=smoothing_scaling_type)
     elseif(smoothing_scaling_type/='' .and. nscale_val>0) then
        write(*,*) "creating linear smoothing constraints, with scaling '"//trim(smoothing_scaling_type)//&
             "' and scaling values ",smoothing_scaling_values
        write(lu_out,*) "creating linear smoothing constraints, with scaling '"//trim(smoothing_scaling_type)//&
             "' and scaling values ",smoothing_scaling_values
        call createNeighbourAverageLinearModelSmoothing(lmsmooth,.invgrid.iterbasics,dmspace,errmsg,&
             scaling_type=smoothing_scaling_type,scaling_values=smoothing_scaling_values)
     end if
     !if (.level.errmsg /= 0) call print(errmsg)
     call print(errmsg)
     if (.level.errmsg == 2) then; close(lu_out); stop; endif
     call dealloc(errmsg)
!
     write(*,*) "there were ",.neq.lmsmooth," linear smoothing constraints created"
     write(lu_out,*) "there were ",.neq.lmsmooth," linear smoothing constraints created"; write(lu_out,*) ""
  end if
!
!------------------------------------------------------------------------
!  read in kernel matrix
!
  write(*,*) "allocating kernel linear system now"
  call new(errmsg,myname)
  call allocateMatrixSerialKernelLinearSystem(KLSE,.ndata.dmspace,.nparam.dmspace,.neq.lmsmooth,errmsg)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!
  write(*,*) "reading in kernel matrix now"
  call new(errmsg,myname)
  lu1 = get(fuh)
  lu2 = get(fuh)
  call readMatrixSerialKernellLinearSystem(KLSE,sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       sval(.inpar.invbasics,'PATH_STATION_FILTER'),rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
       ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
       .comptrans.invbasics,.ncell.(.invgrid.iterbasics),lu1,lu2,dmspace,errmsg)
  call add(fuh,lu1); call add(fuh,lu2)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!------------------------------------------------------------------------
!  read in measured and synthetic data, compute difference residual and misfit
!
  write(*,*) "reading in measured data now"
  call new(errmsg,myname)
  call readMeasuredDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),dmspace,errmsg)
  call undo(fuh)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!
  write(*,*) "reading in synthetic data now"
  call new(errmsg,myname)
  call readSyntheticDataSerialKernelLinearSystem(KLSE,sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       sval(.inpar.invbasics,'PATH_STATION_FILTER'),ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
       ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),.comptrans.invbasics,get(fuh),dmspace,errmsg)
  call undo(fuh)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!
  write(*,*) "setting right-hand-side to data residuals now"
  call new(errmsg,myname)
  call setRhsAsDataResidualKernelLinearSystem(KLSE,errmsg)
  !if (.level.errmsg /= 0) call print(errmsg)
  call print(errmsg)
  if (.level.errmsg == 2) then; close(lu_out); stop; endif
  call dealloc(errmsg)
!
   misfit = getMisfitKernelLinearSystem(KLSE,iostat=ios)
   if(ios/=0) then
      write(*,*) "there was an error computing the misfit, raised iostat = ",ios
      write(lu_out,*) "there was an error computing the misfit, raised iostat = ",ios
      close(lu_out); stop
   end if
   write(*,*) "the misfit (i.e. sum of squares of components of residual vector) is ",misfit
   write(lu_out,*) "the misfit (i.e. sum of squares of components of residual vector) is ",misfit; write(lu_out,*) ""
!------------------------------------------------------------------------
!  add smoothing constraints to kernel linear system
!
   if(add_smoothing) then
      write(*,*) "adding linear smoothing constrains to kernel linear system now"
      call new(errmsg,myname)
      call addToKernelLinearSystemLinearModelSmoothing(lmsmooth,KLSE,dmspace,errmsg)
      !if (.level.errmsg /= 0) call print(errmsg)
      call print(errmsg)
      if (.level.errmsg == 2) then; close(lu_out); stop; endif
         call dealloc(errmsg)
   end if
!------------------------------------------------------------------------
!  solve kernel linear system now
!
   write(*,*) "solve kernel linear system now, having ",.nrow.KLSE," rows and ",.ncol.KLSE," columns"
   write(lu_out,*) "solve kernel linear system now, having ",.nrow.KLSE," rows and ",.ncol.KLSE," columns"
   call new(errmsg,myname)
   call solveSerialKernellLinearSystem(KLSE,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
   dparam => .sol.KLSE
   if(associated(dparam)) then
      write(*,*) "there are a total of ",size(dparam,1)," model update values in the solution vector"
      write(lu_out,*) "there are a total of ",size(dparam,1)," model update values in the solution vector"
   else
      write(*,*) "after solving linear system: solution vector not associated! There must have gone sth. wrong"
      write(lu_out,*) "after solving linear system: solution vector not associated! There must have gone sth. wrong"
      close(lu_out)
      stop 
   end if
!
   write(lu_out,*) "successfully solved linear system"; write(lu_out,*) ""
!------------------------------------------------------------------------
!  retrieve model update from solution vector
!
   write(*,*) "creating model objects from solution vector and reference model and computing new model"
!
! create kernel_inverted_model object model update
!
   call new(errmsg,myname)
   call init(kim_up,.pmtrz.dmspace,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
   do while(nextParamModelParametrization(.pmtrz.dmspace,param_name))
      if(associated(pparam)) deallocate(pparam)
      if(associated(pcell)) deallocate(pcell); nullify(pcell)
      if(associated(idx_dmspace)) deallocate(idx_dmspace)
      allocate(pparam(1)); pparam(1) = param_name
      idx_dmspace => getIndxModelParam(dmspace,param=pparam,cell=pcell)
      if(associated(idx_dmspace)) then
         write(*,*) "there are ",size(idx_dmspace)," '"//trim(param_name)//"' model update values in the solution vector"
         write(lu_out,*) "there are ",size(idx_dmspace)," '"//trim(param_name)//"' model update values in the solution vector"
         call new(errmsg,myname)
         call addValuesKernelInvertedModel(kim_up,param_name,pcell,dparam(idx_dmspace,1),errmsg)
         if (.level.errmsg /= 0) call print(errmsg)
         !call print(errmsg)
         if (.level.errmsg == 2) then; close(lu_out); stop; endif
         call dealloc(errmsg)
         deallocate(idx_dmspace)
      else
         write(*,*) "there are no '"//trim(param_name)//"' model update values in the solution vector"
         write(lu_out,*) "there are no '"//trim(param_name)//"' model update values in the solution vector"
      end if
   end do ! while next param
   write(lu_out,*) ""
!
! write update to file(s) now, as model object kim_up will be modified later on when adding kim_ref to kim_up (to get new model)
!
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim_up,trim(outfile)//'_up.kim',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim_up,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile)//'_up',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
! create kernel_inverted_model object reference model 
!
   call new(errmsg,myname)
   call interpolateKernelReferenceToKernelInvertedModel(kim_ref,.krm.iterbasics,.pmtrz.dmspace,&
        .invgrid.iterbasics,.intw.iterbasics,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
! create kernel_inverted_model object new model = reference + update by modifying object kim_up
! print new max/min model values
! write new model to files
!
   write(*,*) "updating model now, computing new model values"
   write(lu_out,*) "updating model now, computing new model values"
   call new(errmsg,myname)
   call summateInstancesKernelInvertedModel(kim_up,kim_ref,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
   do while(nextParamModelParametrization(.pmtrz.kim_up,param_name))
         write(*,*) "max('"//trim(param_name)//"') = ",maxValue(kim_up,param_name),&
              "min('"//trim(param_name)//"') = ",minValue(kim_up,param_name)
         write(lu_out,*) "max('"//trim(param_name)//"') = ",maxValue(kim_up,param_name),&
              "min('"//trim(param_name)//"') = ",minValue(kim_up,param_name)
   end do
   write(lu_out,*) ""
!
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim_up,trim(outfile)//'_new.kim',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim_up,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),trim(outfile)//'_new',get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) then; close(lu_out); stop; endif
   call dealloc(errmsg)
!
   write(lu_out,*) "successfully written all output files"; write(lu_out,*) ""
!
!------------------------------------------------------------------------
!  clean up
!
   write(lu_out,*) "good bye"; write(lu_out,*) ""
   close(lu_out); call add(fuh,lu_out)
!
   if(allocated(smoothing_scaling_values)) deallocate(smoothing_scaling_values)
   if(associated(pparam)) deallocate(pparam)
   if(associated(idx_dmspace)) deallocate(idx_dmspace)
   if(associated(pcell)) deallocate(pcell)
!
   call dealloc(KLSE)
   call dealloc(dmspace)
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
!   
   call dealloc(kim_up); call dealloc(kim_ref)
!
end program solveKernelSystem
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  print '(50(1h-))'
  print *,'    solveKernelSystem [-h] [-smooth] [-scltyp smoothing_scaling_type] '
  print *,'        [-sclval smoothing_scaling_values] [-smbnd] [-odir] dmspace_file outfile_base parfile'
  print *,''
  print *,'Arguments:'
  print *,''
  print *,"    dmspace_file: data model space input file which defines data and model space"
  print *,"    outfile_base: base name of output files (will be used for all files, with suitable extensions)"
  print *,"    parfile: main parameter file of inversion"
  print *,''
  print *,'Options:'
  print *,''
  print *,'-h     : print help'
  print *,''
  print *,'-smooth  : indicates if linear smoothing constraints (average neighbour) are added to the kernel linear system'
  print *,''
  print *,"-scltyp  : type of scaling of smoothing constraints, at the moment only 'absmax_per_param,overall_factor' "//&
       "(requires one value of -sclval) allowed."
  print *,''
  print *,'-sclval  : dependent on -scltyp, smoothing_scaling_values is of form "n val1 .. valn" giving the number of '//&
       'values n followed by n real numbers'
  print *,''
  print *,'-smbnd  : is ignored at the moment, will in the future define the way (non-existing) neighbours are treated '//&
       'at the boundaries of the inversion grid'
  print *,''
  print *,'-odir  : if set, outfile_base will be assumed relatively to iteration step output files directory'
  print '(50(1h-))'
  return
end subroutine printhelp
