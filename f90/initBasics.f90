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
program initBasics
  use inversionBasics
  use iterationStepBasics
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

   type (argument_parser) :: ap
   character(len=max_length_string) :: main_parfile

   type (file_unit_handler) :: fuh
   type (error_message) :: errmsg
   character(len=10) :: myname = 'initBasics'

   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics

   logical :: recreate_files

!------------------------------------------------------------------------
!  preliminary processing
!
   call init(ap,myname,"Initiating and testing all basic requirements for ASKI programs (parameter files, event "//&
        "and station list, inversion grid, wavefield points, integration weights, reference model)")
   call addPosarg(ap,"main_parfile","sval","Main parameter file of inversion")
   call addOption(ap,"-recr",.false.,"If set, existing files (inversion grid, integration weights, vtk files) "//&
        "will be recreated. If not set, existing files will only be read in!")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   main_parfile = ap.sval.'main_parfile'
   recreate_files = ap.optset.'-recr'
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
!  setup basics
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),get(fuh),errmsg)
   call undo(fuh)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg,recreate_existing_files=recreate_files)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!------------------------------------------------------------------------
!  clean up
!
1  call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
end program initBasics
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'                   initBasics [-h] [-recr] parfile'
!   print *,''
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-recr  : if set, existing files will be recreated with current parfile specifications (overwrites existing files'
!   print *,'         will be read in and not newly created).'
!   print *,'         affected files are:'
!   print *,'           - inversion_grid and integration_weigts and all related .vtk files (including wavefield_points .vtk files)'
!   print *,'           - stations.vtk, events.vtk'
!   print *,'           - kernel reference model files (on wavefield points, and as interpolation on inversion grid)'
!   print *,'         IF NOT SET (default), EXISTING FILES (especially inversion_grid, integration_weights) WILL BE READ IN ONLY,'
!   print *,'         REGARDLESS OF ANY CHANGES IN PARFILES! So, if you change the specification of inversion grid, or integration'
!   print *,'         weights (or stations, events which is only relevant for station.vtk events.vtk), after having already'
!   print *,'         created the respective files, you should set -recr'
!   print *,''
!   print *,'-h     : print help'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
