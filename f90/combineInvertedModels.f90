!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
program combineInvertedModels
   use inversionBasics
   use iterationStepBasics
   use kernelInvertedModel
   use argumentParser
   use string
   use fileUnitHandler
   use errorMessage

   implicit none

   type (argument_parser) :: ap
   character(len=max_length_string) :: main_parfile,kim1_file,kim2_file,kim_outfile!,string

   type (file_unit_handler) :: fuh
   type (error_message) :: errmsg
   character(len=21) :: myname = 'combineInvertedModels'

   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics

   type (kernel_inverted_model) :: kim1,kim2

   logical :: compute_relative
   real :: c1,c2


!------------------------------------------------------------------------
!  preliminary processing
!
   call init(ap,myname,'Compute linear combination   coef1 * kim1 + coef2 * kim2   of two models kim1, kim2')
   call addPosarg(ap,'kim1_file','sval','file of first kernel inverted model')
   call addPosarg(ap,'kim2_file','sval','file of second kernel inverted model')
   call addPosarg(ap,'kim_outfile','sval','output file basename of linear combination of kim1 and kim2')
   call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
   call addOption(ap,'-c1',.true.,"real number which is the coefficient of kim1 in the linear combination",'rval','1.0')
   call addOption(ap,'-c2',.true.,"real number which is the coefficient of kim2 in the linear combination",'rval','1.0')
   call addOption(ap,'-rel',.false.,"if set, the program computes relative to kim1, i.e. (coef1*kim1 + coef2*kim2) / kim1")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   kim1_file = ap.sval.'kim1_file'
   kim2_file = ap.sval.'kim2_file'
   main_parfile = ap.sval.'main_parfile'
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   c1 = ap.rval.'-c1'
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   c2 = ap.rval.'-c2'
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   compute_relative = ap.optset.'-rel'
!
   call document(ap)
   write(*,*) ""
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
   write(*,*) "setting up basics"
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
   write(*,*) ""
!
!------------------------------------------------------------------------
!  read kernel inverted model files kim1_file,kim2_file 
!
   write(*,*) "reading in kernel inverted model files for kim1,kim2"
   call new(errmsg,myname)
   call readFileKernelInvertedModel(kim1,kim1_file,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   call new(errmsg,myname)
   call readFileKernelInvertedModel(kim2,kim2_file,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
   write(*,*) ""

!
!------------------------------------------------------------------------
!  linearly combine models and create new model 
!
   write(*,*) "computing linear combination"
   call new(errmsg,myname)
   call summateInstancesKernelInvertedModel(kim1,kim2,errmsg,c1=c1,c2=c2,relative=compute_relative)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
   write(*,*) ""
!
!------------------------------------------------------------------------
!  write combined model (contained in object kim1, as it was overwritten/modified in routine
!  summateInstancesKernelInvertedModel) to file
!
   write(*,*) "writing output files"
   call new(errmsg,myname)
   call writeFileKernelInvertedModel(kim1,trim(kim_outfile)//'.kim',get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   call new(errmsg,myname)
   call writeVtkKernelInvertedModel(kim1,.invgrid.iterbasics,&
        trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),kim_outfile,get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
   write(*,*) ""
   write(*,*) "good bye"
!
!------------------------------------------------------------------------
!  clean up
!
1  call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(kim1); call dealloc(kim2)
   call dealloc(fuh)
   call dealloc(ap)
!
end program combineInvertedModels
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'Compute linear combination   coef1 * kim1 + coef2 * kim2   of two models kim1, kim2'
!   print *,''
!   print *,'Usage:'
!   print *,''
!   print *,'    combineInvertedModels [-h] [-c1 coef1] [-c2 coef2] [-rel] kim1_file kim2_file kim_ outfile parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    kim1_file: file of first kernel inverted model"
!   print *,''
!   print *,"    kim1_file: file of first kernel inverted model"
!   print *,''
!   print *,"    kim_outfile: output file basename of linear combination of kim1 and kim2"
!   print *,''
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-c1     : coef1 gives a real number which is the coefficient of kim1 in the linear combination'
!   print *,''
!   print *,'-c2     : coef2 gives a real number which is the coefficient of kim2 in the linear combination'
!   print *,''
!   print *,'-rel     : if set, the program computes relative to kim1:   (coef1*kim1 + coef2*kim2) / kim1'
!   print *,''
!   print *,'-h     : print help'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
