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
program exportKim

  use inversionBasics
  use iterationStepBasics
  use kernelInvertedModel
  use modelParametrization
  use inversionGrid
  use vectorPointer
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=300) :: parfile,outfile,kim_file

  type (file_unit_handler) :: fuh

  type (error_message) :: errmsg
  character(len=9) :: myname = 'exportKim'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  ! kernel inverted model
  type (kernel_inverted_model) :: kim
  character(len=character_length_pmtrz) :: pmtrz
  character(len=character_length_param) :: param
  integer :: nparam_pmtrz,nval
  real, dimension(:), pointer :: model_values
  integer, dimension(:), pointer :: indx

  ! output file
  character(len=400) :: second_line
  integer :: lu,ios
  integer :: icell,nnb
  real :: c1,c2,c3,r
  type (integer_vector_pointer), dimension(:), pointer :: nb_idx

  logical :: outfile_exists,stop_after_command_line

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call new(cl,3,1,'h kim o','0 1 1',printhelp)
!
  stop_after_command_line = .false.
!
  parfile = clManarg(cl,1)
!
  ! kim_file
  if(clOptset(cl,2)) then
     kim_file = clOptarg(cl,2)
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -kim"
  end if
!
  ! outfile
  if(clOptset(cl,3)) then
     outfile = clOptarg(cl,3)
     ! check if output model file exists
     inquire(file=outfile,exist=outfile_exists)
     if(outfile_exists) then
        write(*,*) "output file '"//trim(outfile)//"' already exists. Please (re)move it"
        stop_after_command_line = .true.
     end if
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -o"
  end if
!
  if(stop_after_command_line) then
     write(*,*) ""
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
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully initated inversion basics"
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully initated iteration step basics"
!
!------------------------------------------------------------------------
!  read kernel inverted model file
!
  call new(errmsg,myname)
  call readFileKernelInvertedModel(kim,kim_file,get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully read kernel inverted model from file '"//trim(kim_file)//"'"
!
!------------------------------------------------------------------------
!  export inversion grid und kernel inverted model information to text file
!
  pmtrz = .pmtrz.kim
  if(.not.validModelParametrization(pmtrz)) then
     write(*,*) "kernel inverted model has invalid parametrization '"//trim(pmtrz)//"'"
     goto 1
  end if
  nparam_pmtrz = numberOfParamModelParametrization(pmtrz)
  write(second_line,*) nparam_pmtrz
  do while(nextParamModelParametrization(pmtrz,param))
     second_line = trim(second_line)//'   '//trim(param)
  end do
!
  call getIndicesFaceNeighboursInversionGrid(.invgrid.iterbasics,nb_idx)
  if(.not.associated(nb_idx)) then
     write(*,*) "no neighbour information returned from inversion grid, this means "//&
          "inversion grid is not yet defined"
     goto 1
  end if
!
  lu = get(fuh)
  open(unit=lu,file=outfile,form='formatted',status='new',action='write',iostat=ios)
  if(ios/=0) then
     write(*,*) "cannot open output file '"//trim(outfile)//"' to write, raised iostat = ",ios
     close(lu)
     goto 1
  end if
!
  write(lu,*) trim(pmtrz)
!
  write(lu,*) trim(second_line)
!
  write(lu,*) .ncell.(.invgrid.iterbasics)
!
  do icell = 1,.ncell.(.invgrid.iterbasics)
     call new(errmsg,myname)
     call getCenterCellInversionGrid(.invgrid.iterbasics,icell,c1,c2,c3,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu); goto 1; endif
     call dealloc(errmsg)
!     
     call new(errmsg,myname)
     call getRadiusCellInversionGrid(.invgrid.iterbasics,icell,r,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu); goto 1; endif
     call dealloc(errmsg)
!
     indx => getVectorPointer(nb_idx(icell))
     if(associated(indx)) then
        nnb = size(indx)
        write(lu,*) c1,c2,c3,r,nnb,indx
     else
        nnb = 0
        write(lu,*) c1,c2,c3,r,nnb
     end if
  end do ! icell
!
  do while(nextParamModelParametrization(pmtrz,param))
     write(lu,*) param
     indx => getIndx(kim,param)
     if(associated(indx)) then
        nval = size(indx)
        write(lu,*) nval
        write(lu,*) indx
        model_values => getVal(kim,param)
        write(lu,*) model_values
     else
        nval = 0
        write(lu,*) nval
     end if
  end do
!
  close(lu)
  call add(fuh,lu)
  write(*,*) "successfully wrote output file '"//trim(outfile)//"'"
!
!------------------------------------------------------------------------
!  clean up
!
1 write(*,*) "good bye"
!
  if(associated(nb_idx)) then
     do nnb=1,size(nb_idx)
        call dealloc(nb_idx(nnb))
     end do
     deallocate(nb_idx)
  end if
!
  call dealloc(invbasics); call dealloc(iterbasics)
  call dealloc(fuh)
  call dealloc(cl)
!   
  call dealloc(kim)
!
end program exportKim
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  print '(50(1h-))'
  print *,'Program exportKim produces a text file containing information of cell centers and'
  print *,'radii (i.e. rough expansion) of all inversion grid cells and all cell neighbours,'
  print *,'as well as model values and respective invgrid cell indices for each model parameter,'
  print *,'as contained in a given kernelInvertedModel file. This text file may be used by any'
  print *,'forward method to define the simulation model for the next iteration step, or by any'
  print *,'tool handling final models.'
  print *,'The format of the produced text file is as follows:'
  print *,'The first 3 lines contain:'
  print *,'   model_parametrization'
  print *,'   number_of_parameters  name_param_1 ... name_param_n'
  print *,'   number_of_invgrid_cells'
  print *,'The next number_of_invgrid_cells lines contain for each inversion grid cell:'
  print *,'  c1 c2 c3 r nnb inb1 ... inbn'
  print *,'whereby c1,c2,c3 are the 3 coordinates of the cell center (in wavefield point coords),'
  print *,'r is the cell radius (i.e. rough expansion of cell), nnb is the number of cell neighbours'
  print *,'and inb1,...,inbn are their nnb cell indices (if nnb>0, otherwise line ends on nnb=0).'
  print *,'Then number_of_parameters blocks are following in the file, each having the following format:'
  print *,'   name_param'
  print *,'   nval'
  print *,'   cell_indx'
  print *,'   model_values'
  print *,'whereby name_param is the name of the parameter to which the following model values belong,'
  print *,'nval is the number of model values following, cell_indx is a vector of nval cell indices to'
  print *,'which the model values belong (space separated on one line) and model_values is an vector of'
  print *,'the actual nval model values (space separated on one line)'
  print *,''
  print *,''
  print *,'Usage:        exportKim [-h] -kim kernel_inverted_mode_file -o outfile parfile'
  print *,''
  print *,'Arguments:'
  print *,''
  print *,"    parfile: main parameter file of inversion"
  print *,''
  print *,'Mandatory options:'
  print *,''
  print *,'-kim   : kernel_inverted_mode_file is the binary file containing the inverted model, which is to be exported'
  print *,''
  print *,'-o     : outfile is the text output file to which inversion grid and model information will be written'
  print *,''
  print *,'Optional options:'
  print *,''
  print *,'-h     : print help'
  print *,''
  print '(50(1h-))'
  return
end subroutine printhelp
