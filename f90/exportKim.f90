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
program exportKim

  use inputParameter
  use kernelInvertedModel
  use modelParametrization
  use inversionGrid
  use vectorPointer
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,outfile_txt,outfile_vtk,kim_file

  type (file_unit_handler) :: fuh

  type (error_message) :: errmsg
  character(len=9) :: myname = 'exportKim'

  ! parameter file(s)
  type (input_parameter) :: main_inpar
  character (len=80), dimension(5) :: main_inpar_keys
  data main_inpar_keys/'CURRENT_ITERATION_STEP', 'DEFAULT_VTK_FILE_FORMAT','ITERATION_STEP_PATH', &
       'MAIN_PATH_INVERSION', 'PARFILE_ITERATION_STEP'/
  integer :: i_iter
  character(len=350) :: iter_path
  type (input_parameter) :: iter_inpar
  character (len=80), dimension(2) :: iter_inpar_keys
  data iter_inpar_keys/'TYPE_INVERSION_GRID', 'PARFILE_INVERSION_GRID'/

  ! inversion grid
  type (inversion_grid) :: invgrid

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

  logical :: outfile_exists,stop_after_command_line,export_txt,export_vtk,terminate_program

  nullify(model_values,indx,nb_idx)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"exports the given .kim file to vtk (-ovtk) and/or to a text file (-otxt, containing cell "//&
       "centers, radii, neighbours as well as all model values, may be used for communicating inverted model to "//&
       "forward method)")
  call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
  call addOption(ap,"-kim",.true.,"(mandatory) binary '.kim' file containing the inverted model, which is to be "//&
       "exported","sval","")
  call addOption(ap,"-otxt",.true.,"(optional) file name of output text file. If not set, no text file will be "//&
       "produced","sval","")
  call addOption(ap,"-ovtk",.true.,"(optional) file base of output vtk files. If not set, no vtk output will be "//&
       "produced","sval","")
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
  ! kim_file
  if(ap.optset.'-kim') then
     kim_file = ap.sval.'-kim'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -kim"
  end if
!
  ! outfile txt
  export_txt = ap.optset.'-otxt'
  if(export_txt) then
     outfile_txt = ap.sval.'-otxt'
     ! check if output model file exists
     inquire(file=outfile_txt,exist=outfile_exists)
     if(outfile_exists) then
        write(*,*) "output txt file '"//trim(outfile_txt)//"' already exists. Please (re)move it"
        stop_after_command_line = .true.
     end if
  end if
!
  ! outfile vtk
  export_vtk = ap.optset.'-ovtk'
  if(export_vtk) outfile_vtk = ap.sval.'-ovtk'
!
  if(.not.(export_txt .or. export_vtk)) then
     stop_after_command_line = .true.
     write(*,*) "none of -otxt , -ovtk is indicated, so there is nothing to do for this program"
  end if
!
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  if(stop_after_command_line) then
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics (manually, since this program might be needed when not all 
!  prequisites of inversionBasics/iterationStepBasics are ready yet)
!
  ! read in main parameter file (only those keys which are needed here)
  call new(errmsg,myname)
  call createKeywordsInputParameter(main_inpar,main_inpar_keys)
  call readSubroutineInputParameter(main_inpar,get(fuh),trim(main_parfile),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully read main parameter file"
!
  i_iter = ival(main_inpar,'CURRENT_ITERATION_STEP',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: parameter 'CURRENT_ITERATION_STEP' = '"//trim(main_inpar.sval.'CURRENT_ITERATION_STEP')//&
          "' of main parfile is not a valid integer value"
     goto 1
  end if
  write(iter_path,"(2a,i3.3,a)") trim(main_inpar.sval.'MAIN_PATH_INVERSION'),&
       trim(main_inpar.sval.'ITERATION_STEP_PATH'),i_iter,'/'
!
  select case(main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT')
  case('BINARY','ASCII')
     ! OK do nothing
  case default
     write(*,*) "ERROR: parameter 'DEFAULT_VTK_FILE_FORMAT' = '"//trim(main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT')//&
          "' of main parfile is not valid; must be either 'ASCII' or 'BINARY'"
     goto 1       
  end select
!
  ! read in iteration step parameter file (only those keys which are needed here)
  call new(errmsg,myname)
  call createKeywordsInputParameter(iter_inpar,iter_inpar_keys)
  call readSubroutineInputParameter(iter_inpar,get(fuh),trim(iter_path)//&
       trim(main_inpar.sval.'PARFILE_ITERATION_STEP'),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully read iteration step parameter file"
!
  call new(errmsg,myname)
  call createInversionGrid(invgrid,iter_inpar.sval.'TYPE_INVERSION_GRID',&
       trim(iter_path)//trim(iter_inpar.sval.'PARFILE_INVERSION_GRID'),iter_path,&
       get(fuh),errmsg,recreate=.false.)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully set up inversion grid"
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
!------------------------------------------------------------------------
!  export kernel inverted model information to vtk file(s)
!
  if(export_vtk) then
     call new(errmsg,myname)
     call writeVtkKernelInvertedModel(kim,invgrid,trim(main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT'),&
          outfile_vtk,get(fuh),errmsg,overwrite=.false.)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
     write(*,*) "successfully written kim to vtk files, basename '"//trim(outfile_vtk)//"'"
  end if
!------------------------------------------------------------------------
!  export inversion grid und kernel inverted model information to text file
!
  terminate_program = .false.
  if(export_txt) call exportKimToTextfile()
  if(terminate_program) goto 1
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "good bye"
!
1 if(associated(nb_idx)) then
     do nnb=1,size(nb_idx)
        call dealloc(nb_idx(nnb))
     end do
     deallocate(nb_idx)
  end if
!
  call dealloc(errmsg)
  call dealloc(main_inpar); call dealloc(iter_inpar)
  call dealloc(invgrid)
  call dealloc(fuh)
  call dealloc(ap)
!   
  call dealloc(kim)
!

contains

subroutine exportKimToTextfile()

  pmtrz = .pmtrz.kim
  if(.not.validModelParametrization(pmtrz)) then
     write(*,*) "ERROR: kernel inverted model has invalid parametrization '"//trim(pmtrz)//"'"
     goto 1
  end if
  nparam_pmtrz = numberOfParamModelParametrization(pmtrz)
  write(second_line,*) nparam_pmtrz
  do while(nextParamModelParametrization(pmtrz,param))
     second_line = trim(second_line)//'   '//trim(param)
  end do
!
  call getIndicesFaceNeighboursInversionGrid(invgrid,nb_idx)
  if(.not.associated(nb_idx)) then
     write(*,*) "ERROR: no neighbour information returned from inversion grid, this means "//&
          "inversion grid is not yet defined"
     goto 1
  end if
!
  lu = get(fuh)
  open(unit=lu,file=outfile_txt,form='formatted',status='new',action='write',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: cannot open output file '"//trim(outfile_txt)//"' to write, raised iostat = ",ios
     close(lu)
     goto 1
  end if
!
  write(*,*) "export model now:"
!
  write(*,*) "   ",trim(pmtrz)
  write(lu,*) trim(pmtrz)
!
  write(*,*) "   ",trim(second_line)
  write(lu,*) trim(second_line)
!
  write(*,*) "   ",.ncell.invgrid," inversion grid cells"
  write(lu,*) .ncell.invgrid
!
  do icell = 1,.ncell.invgrid
     call new(errmsg,myname)
     call getCenterCellInversionGrid(invgrid,icell,c1,c2,c3,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) then; close(lu); goto 1; endif
     call dealloc(errmsg)
!     
     call new(errmsg,myname)
     call getRadiusCellInversionGrid(invgrid,icell,r,errmsg)
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
        write(*,*) "   parameter ",trim(param),": there are model values on ",nval," cells"
        if(nval<.ncell.invgrid) then
           write(*,*) "     WARNING: number of model values of this parameter is smaller than number ",&
                "of inversion grid cells; make sure you know what you are doing"
        end if
        write(lu,*) nval
        write(lu,*) indx
        model_values => getVal(kim,param)
        write(lu,*) model_values
     else
        nval = 0
        write(*,*) "   parameter ",trim(param),": THERE ARE NO MODEL VALUES"
        write(lu,*) nval
     end if
  end do
!
  close(lu)
  call add(fuh,lu)
  write(*,*) "successfully wrote output file '"//trim(outfile_txt)//"'"

  ! if subroutine comes here, everything went OK, so return to main program
  return
!
  ! if an error has occurred make the main program terminate after return to main program
1 terminate_program = .true.
  return
end subroutine exportKimToTextfile

end program exportKim
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'Program exportKim produces a text file (option -otxt) containing information of cell centers and'
!   print *,'radii (i.e. rough expansion) of all inversion grid cells and all cell neighbours,'
!   print *,'as well as model values and respective invgrid cell indices for each model parameter,'
!   print *,'as contained in a given kernelInvertedModel file. This text file may be used by any'
!   print *,'forward method to define the simulation model for the next iteration step, or by any'
!   print *,'tool handling final models.'
!   print *,'The format of the produced text file is as follows:'
!   print *,'The first 3 lines contain:'
!   print *,'   model_parametrization'
!   print *,'   number_of_parameters  name_param_1 ... name_param_n'
!   print *,'   number_of_invgrid_cells'
!   print *,'The next number_of_invgrid_cells lines contain for each inversion grid cell:'
!   print *,'  c1 c2 c3 r nnb inb1 ... inbn'
!   print *,'whereby c1,c2,c3 are the 3 coordinates of the cell center (in wavefield point coords),'
!   print *,'r is the cell radius (i.e. rough expansion of cell), nnb is the number of cell neighbours'
!   print *,'and inb1,...,inbn are their nnb cell indices (if nnb>0, otherwise line ends on nnb=0).'
!   print *,'Then number_of_parameters blocks are following in the file, each having the following format:'
!   print *,'   name_param'
!   print *,'   nval'
!   print *,'   cell_indx'
!   print *,'   model_values'
!   print *,'whereby name_param is the name of the parameter to which the following model values belong,'
!   print *,'nval is the number of model values following, cell_indx is a vector of nval cell indices to'
!   print *,'which the model values belong (space separated on one line) and model_values is an vector of'
!   print *,'the actual nval model values (space separated on one line)'
!   print *,''
!   print *,'Additionally (or alternatively) the program converts the .kim file to vtk files (option -ovtk)'
!   print *,''
!   print *,'The two options -otxt , -ovtk can be used independently of each other'
!   print *,''
!   print *,'Usage:        exportKim [-h] -kim kernel_inverted_mode_file -otxt outfile_txt -ovtk outfile_vtk parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Mandatory options:'
!   print *,''
!   print *,'-kim   : kernel_inverted_mode_file is the binary file containing the inverted model, which is to be exported'
!   print *,''
!   print *,'-otxt     : outfile_txt is the text output file to which inversion grid and model information will be written'
!   print *,''
!   print *,'-ovtk     : outfile_vtk is the output file base to which standard vtk files of the .kim file will be written'
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print '(50(1h-))'
!   return
! end subroutine printhelp
