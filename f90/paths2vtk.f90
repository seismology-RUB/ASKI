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
program paths2vtk

  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use eventStationVtkFile
  use argumentParser
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=400) :: main_parfile,dmspace_file,outfile

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=9) :: myname = 'paths2vtk'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths

  type (event_station_vtk_file) :: evstat_vtk
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,'plot all paths contained in the data space definition as vtk lines')
  ! define optional arguments
  ! call addOption(ap,'-smooth',.false.,&
  !      'indicates if linear smoothing constraints (average neighbour) are added to the kernel linear system')
  ! call addOption(ap,'-scltyp',.true.,&
  !      "type of scaling of smoothing constraints, at the moment only 'absmax_per_param,overall_factor' "//&
  !      "(requires one value of -sclval) allowed.",'sval','absmax_per_param,overall_factor')
  ! define positional arguments
  call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
  call addPosarg(ap,'outfile','sval','output file (basename) of vtk file(s)')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if

  ! get values of positional arguments
  main_parfile = ap.sval.'main_parfile'
  dmspace_file = ap.sval.'dmsi_file'
  outfile = ap.sval.'outfile'

  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  setup data space info object
!
  write(*,*) "creating data space info from file '"//trim(dmspace_file)//"'"
!
  call new(errmsg,myname)
  call createDataSamplesFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  write(*,*) "there are ",.ndata.dmspace," data samples from ",.npath.dmspace," paths"
  write(*,*) ""
!
  ! get all the paths contained in the data space
  paths => getPathsDataModelSpaceInfo(dmspace)
  if(.not.associated(paths)) then
     write(*,*) "ERROR: no paths contained in data space, data space seems to be empty"
     goto 1
  end if
!
  ! first of all output one vtk file (without data), containing all paths of the data model space
  call new(errmsg,myname)
  call init(evstat_vtk,.evlist.invbasics,.statlist.invbasics,paths,.invgrid.iterbasics,trim(outfile),&
       trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  write(*,*) "writing paths vtk file with base-filename '"//trim(outfile)//"'"
  write(*,*) ""
  call writePaths(evstat_vtk,get(fuh),errmsg,overwrite =.false.)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
1 call dealloc(dmspace)
  call dealloc(iterbasics)
  call dealloc(invbasics)
  call dealloc(fuh)
  call dealloc(errmsg)
  write(*,*) "good bye"; write(*,*) ""
end program paths2vtk
