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
program krm2kim

  use inversionBasics
  use iterationStepBasics
  use kernelInvertedModel
  use kernelReferenceModel
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,outfile,krm_file
  logical :: use_alternative_krm

  type (file_unit_handler) :: fuh

  type (error_message) :: errmsg
  character(len=7) :: myname = 'krm2kim'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (kernel_inverted_model) :: kim
  type (kernel_reference_model) :: krm

  logical :: outfile_exists

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"interpolate kernel reference model onto inversion grid and produce a .kim file of it")
  call addPosarg(ap,"outfile_base","sval","outfile_base")
  call addPosarg(ap,"main_parfile","sval","main_parfile")
  call addOption(ap,"-krm",.true.,"Filename of alternative kernel reference model. If set, this krm "//&
       "is used INSTEAD of the one of the current iteration step.","sval","")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
  main_parfile = ap.sval.'main_parfile'
  outfile = ap.sval.'outfile_base'
!
  use_alternative_krm = ap.optset.'-krm'
  if(use_alternative_krm) krm_file = ap.sval.'-krm'
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
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  check if output files already exist
!
  ! check if output model file exists
  inquire(file=trim(outfile)//'.kim',exist=outfile_exists)
  if(outfile_exists) then
     write(*,*) "inverted model file '"//trim(outfile)//".kim' already exists. Please (re)move it."
     goto 1
  end if
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(main_parfile),get(fuh),errmsg)
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
!------------------------------------------------------------------------
!  create kernel reference model on invgrid
!
  if(use_alternative_krm) then

     call new(errmsg,myname)
     call createKernelReferenceModel(krm,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh,krm_file,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)

     call new(errmsg,myname)
     call interpolateKernelReferenceToKernelInvertedModel(kim,krm,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',&
          .invgrid.iterbasics,.intw.iterbasics,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)

     call dealloc(krm)

  else

     call new(errmsg,myname)
     call interpolateKernelReferenceToKernelInvertedModel(kim,.krm.iterbasics,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',&
          .invgrid.iterbasics,.intw.iterbasics,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)

  end if
  write(*,*) "successfully interpolated reference model to inversion grid"
!
!------------------------------------------------------------------------
!  write kernel inverted model to output file(s)
!
  call new(errmsg,myname)
  call writeFileKernelInvertedModel(kim,trim(outfile)//'.kim',get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call new(errmsg,myname)
  call writeVtkKernelInvertedModel(kim,.invgrid.iterbasics,&
       (.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT',outfile,get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully written all output files"
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "good bye"
!
1 call dealloc(invbasics); call dealloc(iterbasics)
  call dealloc(fuh)
  call dealloc(ap)
   
  call dealloc(kim)
!
end program krm2kim
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'                   krm2kim [-h] [-krm kernel_reference_mode_file] outfile parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    outfile: basename of output model files - will additionally be written as vtk"
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print *,'-krm   : if set, then instead of using the kernel reference model file defined by the iteration step parfile,'
!   print *,'         the given file kernel_reference_mode_file is used to read in the kernel reference model.'
!   print *,''
!   print '(50(1h-))'
!   return
! end subroutine printhelp
