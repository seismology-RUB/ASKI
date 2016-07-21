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
program invgrid2vtk
  use inversionGrid
  use invgridVtkFile
  use argumentParser
  use string
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: type_inversion_grid,inpar_inversion_grid,invgrid_path,outbase,filebase,str
  character(len=200) :: vtk_title

  type (error_message) :: errmsg
  character(len=11) :: myname = 'invgrid2vtk'

  type (inversion_grid) :: invgrid
  type (invgrid_vtk_file) :: invgrid_vtk
  character(len=6) :: vtk_format

  logical :: write_all_neighbours_vtk,write_selected_neighbours_vtk
  integer :: ncell_nb,nnb_actual
  integer, dimension(:), pointer :: icell_nb
  type (integer_vector_pointer), dimension(:), pointer :: nb_idx
  integer, dimension(:), pointer :: nb,nb_actual
  real, dimension(:), allocatable :: nnb_total
  character(len=39) :: bnd_cond

  integer :: j,icell
  logical :: stop_after_command_line,overwrite_output,recreate_invgrid,allocate_nb_actual

  nullify(icell_nb,nb_idx,nb,nb_actual)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Create vtk file(s) of the given inversion grid (useful to look at, to see if the "//&
       "specifications are correct)")
  call addOption(ap,"-igtype",.true.,"(mandatory) TYPE_INVERSION_GRID as in ASKI iteration step parameter file",&
       "sval","")
  call addOption(ap,"-igpar",.true.,"(mandatory) PARFILE_INVERSION_GRID as in ASKI iteration step parameter file",&
       "sval","")
  call addOption(ap,"-igpath",.true.,"(mandatory) is treated as current iteration step path, used for some "//&
       "inversion grids to write/read own files","sval","")
  call addOption(ap,"-o",.true.,"(optional) output flile base name","sval","inversion_grid")
  call addOption(ap,"-overwr",.false.,"(optional) if set, existing output files will be overwritten")
  call addOption(ap,'-nb',.true.,"(optional) explicit vector of cell indices for which the neighbours will be "//&
       "written as vtk files. Options -nb and -all_nb must not be set at the same time",'ivec','')
  call addOption(ap,'-all_nb',.false.,"(optional) if set, for all invgrid cells the neighbours will be written "//&
       "as vtk files. Options -nb and -all_nb must not be set at the same time")
  call addOption(ap,"-bin",.false.,"(optional) if set, the output vtk files will be binary, otherwise they will "//&
       "be ascii")
  call addOption(ap,"-recr",.false.,"(optional) if set, possible existing local inversion grid files will be "//&
       "recreated")
  call addOption(ap,"-bndcond",.true.,"(optional) type of boundary conditions ('extra_nbs_outer_bnd', "//&
       "'extra_nbs_outer_bnd_except_free_surface')","sval","standard")
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
  ! type_inversion_grid
  if(ap.optset.'-igtype') then
     type_inversion_grid = ap.sval.'-igtype'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igtype"
  end if
!
  ! inpar_inversion_grid
  if(ap.optset.'-igpar') then
     inpar_inversion_grid = ap.sval.'-igpar'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpar"
  end if
!
  ! invgrid_path
  if(ap.optset.'-igpath') then
     invgrid_path = ap.sval.'-igpath'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpath"
  end if
!
  write_all_neighbours_vtk = ap.optset.'-all_nb'
  write_selected_neighbours_vtk = ap.optset.'-nb'
  ! either one of -nb , -nb_all can be set, not both at the same time
  if(write_all_neighbours_vtk .and. write_selected_neighbours_vtk) then
     stop_after_command_line = .true.
     write(*,*) "ERROR: Options -nb and -all_nb must not be set at the same time!"
  end if
  ! -bndcond can only be set if there is any neighbour output requested
  if(ap.optset.'-bndcond') then
     if(.not.(write_all_neighbours_vtk .or. write_selected_neighbours_vtk)) then
        stop_after_command_line = .true.
        write(*,*) "ERROR: -bndcond can only be set if some neighbours are requested (by -nb or -all_nb)"
     end if
  end if
!
  if(stop_after_command_line) then
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  ! outbase
  outbase = ap.sval.'-o'
!
  ! overwrite_output
  overwrite_output = ap.optset.'-overwr'
!
  if(write_selected_neighbours_vtk) then
     icell_nb => ap.ivec.'-nb'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(icell_nb)) then
        write(*,*) "ERROR: for some reason, there is no list of cell indeices returned by argument parser, "//&
             "even though there was no error parsing argument -nb. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
     ncell_nb = size(icell_nb)
  end if ! write_selected_neighbours_vtk
!
  ! -bin
  if(ap.optset.'-bin') then
     vtk_format = 'BINARY'
  else
     vtk_format = 'ASCII'
  end if
!
  ! recreate_invgrid
  recreate_invgrid = ap.optset.'-recr'
!
  ! boundary conditions
  str = ap.sval.'-bndcond'
  bnd_cond = str
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
!------------------------------------------------------------------------
!  create inversion grid
!
  call new(errmsg,myname)
  call createInversionGrid(invgrid,type_inversion_grid,inpar_inversion_grid,invgrid_path,11,errmsg,&
       recreate=recreate_invgrid)
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
  if(write_selected_neighbours_vtk) then
     if(any(icell_nb<1 .or. icell_nb>.ncell.invgrid)) then
        write(*,*) "ERROR: among the invgrid cell indices requested by option -nb, there are",&
             count(icell_nb<1 .or. icell_nb>.ncell.invgrid)," invalid indices which are either < 1 or greater ",&
             "than the maximum number of cells of this inversion grid = ",.ncell.invgrid
        goto 1
     end if
  end if ! write_selected_neighbours_vtk
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
     write(*,*) "writing neigbour vtk files for ",size(icell_nb),&
          " invgrid cells with basename '"//trim(outbase)//"_nb######'"
     write(*,*) "boundary condition type = '",trim(bnd_cond),"'"
!
     call getIndicesFaceNeighboursInversionGrid(invgrid,nb_idx,boundary_conditions=bnd_cond)
     if(.not.associated(nb_idx)) then
        write(*,*) "ERROR: no neighbour indices returned by getIndicesFaceNeighboursInversionGrid. It might be that "//&
             "for this type of inversion grid, the requested boundary conditions are not supported ('standard' is "//&
             "supported by all inversion grids)"
        goto 1
     end if
!
     allocate(nnb_total(ncell_nb))
     nullify(nb_actual)
     do j = 1,ncell_nb
        icell = icell_nb(j)
!
        nb => getVectorPointer(nb_idx(icell))
        if(.not.associated(nb)) then
           write(*,*) "there are no neighbours of invgrid cell ",icell
           nnb_total(j) = 0.
           cycle
        end if
!
        ! remember the total number of neighbours of this cell (including artificial neighbours)
        nnb_total(j) = real(size(nb))
!
        if(all(nb==-1)) then
           write(*,*) "there are only artificial neighbours of invgrid cell ",icell
           cycle
        end if
!
        allocate_nb_actual = any(nb==-1)
        if(allocate_nb_actual) then
           nnb_actual = count(nb /= -1)
           allocate(nb_actual(nnb_actual))
           nb_actual = pack(nb,nb /= -1)
        else
           nb_actual => nb
        end if
!
        ! NOW WRITE vtk FILE CONTAINING ALL ACTUAL NEIGHBOURS OF CURRENT CELL (as data, write cell indices on cells)
!
        write(filebase,"(a,i6.6)") trim(outbase)//"_nb",icell
!
        call new(errmsg,myname)
        write(vtk_title,"('inversion grid neighbours of cell',i6)") icell
        call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title=vtk_title,cell_indx_req=nb_actual)
        if(.level.errmsg/=0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 1
!
        call writeData(invgrid_vtk,11,real(nb_actual),errmsg,data_name='cell_index',&
             overwrite=overwrite_output)
        if(.level.errmsg/=0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 1
        call dealloc(errmsg)  
        call dealloc(invgrid_vtk)
!
        if(allocate_nb_actual) then
           deallocate(nb_actual)
           allocate_nb_actual = .false.
        end if
        nullify(nb_actual)
!        
     end do ! j

     ! FINALLY WRITE vtk FILE CONTAINING ALL REQUESTED CELLS icell_nb , AS DATA WRITE TOTAL NUMBER OF NEIGHBOURS, INCLUDING ARTIFICIAL
     filebase = trim(outbase)//"_total_nnb_incl_artificial"
!
     call new(errmsg,myname)
     vtk_title = "total number of neighbours including artificial ones on boundaries"
     call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title=vtk_title,cell_indx_req=icell_nb)
     if(.level.errmsg/=0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
!
     call writeData(invgrid_vtk,11,nnb_total,errmsg,data_name='nnb_total_incl_artif',overwrite=overwrite_output)
     if(.level.errmsg/=0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
     call dealloc(invgrid_vtk)

  end if ! write_all_neighbours_vtk .or. write_selected_neighbours_vtk
!
!------------------------------------------------------------------------
!  clean up
!
1 if(allocate_nb_actual) then
     if(associated(nb_actual)) deallocate(nb_actual)
  end if
  if(allocated(nnb_total)) deallocate(nnb_total)
  if(associated(icell_nb)) deallocate(icell_nb)
  if(associated(nb_idx)) then
     do j = 1,size(nb_idx)
        call dealloc(nb_idx(j))
     end do ! j
     deallocate(nb_idx)
  end if
  call dealloc(invgrid_vtk)
  call dealloc(errmsg)
  call dealloc(invgrid)
  call dealloc(ap)
end program invgrid2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'Usage:'
!   print *,''
!   print *,'    invgrid2vtk -igtype invgrid_type -igpar invgrid_parfile -igpath invgrid_path [-o outbase] [-overwr] '//&
!        '[-nb neighbours] [-bin] [-recr] [-h]'
!   print *,''
!   print *,'Mandatory options:'
!   print *,''
!   print *,'-igtype  : invgrid_type is TYPE_INVERSION_GRID as in ASKI iteration step parameter file'
!   print *,''
!   print *,'-igpar : invgrid_parfile is PARFILE_INVERSION_GRID as in ASKI iteration step parameter file'
!   print *,''
!   print *,'-igpath: invgrid_path is treated as current iteration step path, used for inversion grid to write/read own files'
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'-o     : outbase is output base name (default is inversion_grid)'
!   print *,''
!   print *,'-overwr : if set, existing output files will be overwritten'
!   print *,''
!   print *,'-nb : neighbours is either "all" or "N idx1 idx2 .. idxN", i.e. the number of cells N followed by '
!   print *,'      N cell indices, indicating a set of cells the neighbours of which will be written as vtk files'
!   print *,'      IN THE FUTURE: may be helpful to also accept ranges like "20:40"'
!   print *,''
!   print *,'-bin  : if set, the output vtk files will be binary, otherwise they will be ascii'
!   print *,''
!   print *,'-recr  : if set, the existing inversion grid file(s) (if any existing, and if type creates any) will be '
!   print *,'         recreated with current invgrid parfile specifications'
!   print *,'         if not set (default), existing inversion grid will be read in only, REGARDLESS OF ANY CHANGES IN THE PARFILE!'
!   print *,''
!   print *,'-h     : print this help'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
