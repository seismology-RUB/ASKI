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
program invgrid2vtk
  use inversionGrid
  use invgridVtkFile
  use commandLine
  use errorMessage

  implicit none

  type (cmdLine) :: cl
  character(len=300) :: string,type_inversion_grid,inpar_inversion_grid,outbase,filebase,vtk_title,invgrid_path

  type (error_message) :: errmsg
  character(len=11) :: myname = 'invgrid2vtk'

  type (inversion_grid) :: invgrid
  type (invgrid_vtk_file) :: invgrid_vtk
  character(len=6) :: vtk_format

  logical :: write_all_neighbours_vtk,write_selected_neighbours_vtk
  integer :: ncell_nb
  integer, dimension(:), allocatable :: icell_nb
  type (integer_vector_pointer), dimension(:), pointer :: nb_idx
  integer, dimension(:), pointer :: nb

  integer :: ios,j
  logical :: stop_after_command_line,overwrite_output

  external printhelp

!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call new(cl,8,0,'h igtype igpar o overwr nb bin igpath','0 1 1 1 0 1 0 1',printhelp)
!
  stop_after_command_line = .false.
!
  ! type_inversion_grid
  if(clOptset(cl,2)) then
     type_inversion_grid = clOptarg(cl,2)
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igtype"
  end if
!
  ! inpar_inversion_grid
  if(clOptset(cl,3)) then
     inpar_inversion_grid = clOptarg(cl,3)
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpar"
  end if
!
  ! invgrid_path
  if(clOptset(cl,8)) then
     invgrid_path = clOptarg(cl,8)
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpath"
  end if
!
  ! outbase
  if(clOptset(cl,4)) then
     outbase = clOptarg(cl,4)
  else
     outbase = 'inversion_grid'
  end if
!
  ! overwrite_output
  overwrite_output = clOptset(cl,5)
!
  write_all_neighbours_vtk = .false.; write_selected_neighbours_vtk = .false.
  if(clOptset(cl,6)) then
     string = clOptarg(cl,6)
     if(string=='all') then
        write_all_neighbours_vtk = .true.
     else
        read(string,*,iostat=ios) ncell_nb
        if(ios/=0) then
           stop_after_command_line = .true.
           write(*,*) "ERROR: could not read integer number of cells as first word of -nb input string '"&
                //trim(string)//"'"
        else
           allocate(icell_nb(ncell_nb))
           read(string,*,iostat=ios) ncell_nb,icell_nb
           if(ios/=0) then
              stop_after_command_line = .true.
              write(*,*) "ERROR: could not read ",ncell_nb," integer cell indices from -nb input string '"&
                   //trim(string)//"' starting with second word"
           end if ! ios/=0 read ncell_nb,icell_nb
        end if ! ios/=0 read ncell_nb
     end if ! string=='all'
     write_selected_neighbours_vtk = .true.
  end if ! clOptset(cl,6)
!
  if(clOptset(cl,7)) then
     vtk_format = 'BINARY'
  else
     vtk_format = 'ASCII'
  end if
!
  if(stop_after_command_line) then
     write(*,*) ""
     call printhelp
     stop
  end if
!
!------------------------------------------------------------------------
!  create inversion grid
!
  call new(errmsg,myname)
  call createInversionGrid(invgrid,type_inversion_grid,inpar_inversion_grid,invgrid_path,11,errmsg)
  if(.level.errmsg/=0) call print(errmsg)
  !call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
  if(write_all_neighbours_vtk) then
     ncell_nb = .ncell.invgrid
     allocate(icell_nb(ncell_nb))
     icell_nb = (/ (j,j=1,ncell_nb) /)
  end if
!
!------------------------------------------------------------------------
!  write basic vtk file
!
  write(*,*) "writing inversion grid with invgrid index as data to vtk file with basename '"//trim(outbase)//"'"
!
  call new(errmsg,myname)
  call init(invgrid_vtk,invgrid,outbase,vtk_format,errmsg,'inversion grid cell indices')
  if(.level.errmsg/=0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  call new(errmsg,myname)
  call writeData(invgrid_vtk,11,(/(real(j),j=1,.ncell.invgrid)/),errmsg,data_name='cell_index',&
       overwrite=overwrite_output)
  if(.level.errmsg/=0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)  
  call dealloc(invgrid_vtk)
!
!------------------------------------------------------------------------
!  write neighbour files, if any
!
  if(write_all_neighbours_vtk .or. write_selected_neighbours_vtk) then
     write(*,*) "writing neigbour vtk files for ",ncell_nb,&
          " invgrid cells with basename '"//trim(outbase)//"_nb######'"
!
     call getIndicesFaceNeighboursInversionGrid(invgrid,nb_idx)
     if(.not.associated(nb_idx)) then
        write(*,*) "ERROR: no neighbour indices returned by getIndicesFaceNeighboursInversionGrid"
        stop
     end if
!
     do j = 1,ncell_nb
        if(icell_nb(j)<1 .or. icell_nb(j)>.ncell.invgrid) then
           write(*,*) j,"'th invgrid index ",icell_nb(j),&
                " out of range, must be between 1 and .ncell.invgrid = ",.ncell.invgrid
           cycle
        end if
!
        nb => getVectorPointer(nb_idx(icell_nb(j)))
        if(.not.associated(nb)) then
           write(*,*) "there are no neighbours of invgrid cell ",icell_nb(j)
           cycle
        end if
!
        write(filebase,"(a,i6.6)") trim(outbase)//"_nb",icell_nb(j)
!
        call new(errmsg,myname)
        write(vtk_title,"('inversion grid neighbours of cell',i6)") icell_nb(j)
        call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title=vtk_title,cell_indx_req=nb)
        if(.level.errmsg/=0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 1
        call dealloc(errmsg)
!
        call new(errmsg,myname)
        call writeData(invgrid_vtk,11,real(nb),errmsg,data_name='cell_index',&
             overwrite=overwrite_output)
        if(.level.errmsg/=0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 1
        call dealloc(errmsg)  
        call dealloc(invgrid_vtk)
        
     end do ! j
  end if
!
!------------------------------------------------------------------------
!  clean up
!
1 if(allocated(icell_nb)) deallocate(icell_nb)
  if(associated(nb_idx)) then
     do j = 1,size(nb_idx)
        call dealloc(nb_idx(j))
     end do ! j
     deallocate(nb_idx)
  end if
  call dealloc(invgrid_vtk)
  call dealloc(errmsg)
  call dealloc(invgrid)
  call dealloc(cl)
end program invgrid2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
subroutine printhelp
  print '(50(1h-))'
  print *,'Usage:'
  print *,''
  print *,'    invgrid2vtk -igtype invgrid_type -igpar invgrid_parfile -igpath invgrid_path [-o outbase] [-overwr] '//&
       '[-nb neighbours] [-bin] [-h]'
  print *,''
  print *,'Mandatory options:'
  print *,''
  print *,'-igtype  : invgrid_type is TYPE_INVERSION_GRID as in ASKI iteration step parameter file'
  print *,''
  print *,'-igpar : invgrid_parfile is PARFILE_INVERSION_GRID as in ASKI iteration step parameter file'
  print *,''
  print *,'-igpath: invgrid_path is treated as current iteration step path, used for inversion grid to write/read own files'
  print *,''
  print *,'Optional options:'
  print *,''
  print *,'-o     : outbase is output base name (default is inversion_grid)'
  print *,''
  print *,'-overwr : if set, existing output files will be overwritten'
  print *,''
  print *,'-nb : neighbours is either "all" or "N idx1 idx2 .. idxN", i.e. the number of cells N followed by '
  print *,'      N cell indices, indicating a set of cells the neighbours of which will be written as vtk files'
  print *,'      IN THE FUTURE: may be helpful to also accept ranges like "20:40"'
  print *,''
  print *,'-bin  : if set, the output vtk files will be binary, otherwise they will be ascii'
  print *,''
  print *,'-h     : print this help'
  print '(50(1h-))'
  return
end subroutine printhelp
