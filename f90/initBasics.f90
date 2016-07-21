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
program initBasics
  use inversionBasics
  use iterationStepBasics
  use commandLine
  use fileUnitHandler
  use errorMessage

  implicit none

   type (cmdLine) :: cl
   character(len=132) :: parfile

   type (file_unit_handler) :: fuh
   type (error_message) :: errmsg
   character(len=10) :: myname = 'initBasics'

   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics

   logical :: recreate_files

   external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
   ! process command line
   call new(cl,2,1,'h recr','0 0',printhelp)
   recreate_files = clOptset(cl,2)
   parfile = clManarg(cl,1)
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
   if (.level.errmsg == 2) stop
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg,recreate_existing_files=recreate_files)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) stop
   call dealloc(errmsg)
!------------------------------------------------------------------------
!  clean up
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(cl)
end program initBasics
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  print '(50(1h-))'
  print *,'                   initBasics [-h] [-recr] parfile'
  print *,''
  print *,"    parfile: main parameter file of inversion"
  print *,''
  print *,'Options:'
  print *,''
  print *,'-recr  : if set, existing files will be recreated with current parfile specifications (overwrites existing files'
  print *,'         will be read in and not newly created).'
  print *,'         affected files are:'
  print *,'           - inversion_grid and integration_weigts and all related .vtk files (including wavefield_points .vtk files)'
  print *,'           - stations.vtk, events.vtk'
  print *,'           - kernel reference model files (on wavefield points, and as interpolation on inversion grid)'
  print *,'         IF NOT SET (default), EXISTING FILES (especially inversion_grid, integration_weights) WILL BE READ IN ONLY,'
  print *,'         REGARDLESS OF ANY CHANGES IN PARFILES! So, if you change the specification of inversion grid, or integration'
  print *,'         weights (or stations, events which is only relevant for station.vtk events.vtk), after having already'
  print *,'         created the respective files, you should set -recr'
  print *,''
  print *,'-h     : print help'
  print '(50(1h-))'
  return
end subroutine printhelp
