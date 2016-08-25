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
program chunksInvgrid2vtk
  use inversionGrid
  use chunksInversionGrid
  use invgridVtkFile
  use argumentParser
  use string
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: inpar_inversion_grid,invgrid_path,outbase,filebase
  character(len=200) :: vtk_title

  type (error_message) :: errmsg
  character(len=17) :: myname = 'chunksInvgrid2vtk'

  type (inversion_grid) :: invgrid
  type (chunks_inversion_grid), pointer :: chunks_invgrid
  type (invgrid_vtk_file) :: invgrid_vtk
  character(len=6) :: vtk_format

  logical :: write_base_cells_icell_base,write_base_cells_ichunk,write_cells_ichunk,&
       write_selected_subcells,write_all_subcells,do_anything_at_all

  integer :: ncell_base_for_subs
  integer, dimension(:), pointer :: icell_base_for_subs,idx

  integer :: j,icell_base
  logical :: stop_after_command_line,overwrite_output,recreate_invgrid,stop_program

  nullify(icell_base_for_subs,idx,chunks_invgrid)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Create special vtk file(s) of a chunksInversiongrid-type inversion grid; "//&
       "use this executable complementary to invgrid2vtk, especially in case of base cell refinement.")
  call addOption(ap,"-igpar",.true.,"(mandatory) PARFILE_INVERSION_GRID as in ASKI iteration step parameter file",&
       "sval","")
  call addOption(ap,"-igpath",.true.,"(mandatory) is treated as current iteration step path, used for the "//&
       "inversion grid to write/read own file(s)","sval","")
  call addOption(ap,"-o",.true.,"(optional) output file base name","sval","inversion_grid")
  call addOption(ap,'-base',.false.,"(optional) if set, produce a vtk file containing the inversion grid's base "//&
       "cells, with assigned base cell indices")
  call addOption(ap,'-base_ichk',.false.,"(optional) if set, produces a vtk file containing the inversion grid's base "//&
       "cells, with assigned chunk indices")
  call addOption(ap,'-ichunk',.false.,"(optional) if set, produces a vtk file containing the inversion grid cells, "//&
       "with assigned cell indices")
  call addOption(ap,'-subs',.true.,"(optional) explicit vector of base cell indices for which subcell "//&
       "vtk files (with cell index) will be produced (one file per base cell, contains base cell if not refined)",'ivec','')
  call addOption(ap,'-all_subs',.false.,"(optional) if set, same as for -subs will be done for ALL base cells")
  call addOption(ap,"-overwr",.false.,"(optional) if set, existing output files will be overwritten")
  call addOption(ap,"-bin",.false.,"(optional) if set, the output vtk files will be binary, otherwise they will "//&
       "be ascii")
  call addOption(ap,"-recr",.false.,"(optional) if set, possible existing local inversion grid files will be "//&
       "recreated")
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
  do_anything_at_all = .false.
!
  write_base_cells_icell_base = ap.optset.'-base'
  do_anything_at_all = do_anything_at_all .or. write_base_cells_icell_base
!
  write_base_cells_ichunk = ap.optset.'-base_ichk'
  do_anything_at_all = do_anything_at_all .or. write_base_cells_ichunk
!
  write_cells_ichunk = ap.optset.'-ichunk'
  do_anything_at_all = do_anything_at_all .or. write_cells_ichunk
!
  write_selected_subcells = ap.optset.'-subs'
  write_all_subcells = ap.optset.'-all_subs'
  ! either one of -subs , -all_subs can be set, not both at the same time
  if(write_selected_subcells .and. write_all_subcells) then
     stop_after_command_line = .true.
     write(*,*) "ERROR: Options -subs and -all_subs must not be set at the same time!"
  end if
  do_anything_at_all = do_anything_at_all .or. (write_selected_subcells.or.write_all_subcells)
!
  if(.not.do_anything_at_all) then
     stop_after_command_line = .true.
     write(*,*) "ERROR: there is no option set that requests any vtk files to be produced, so there is nothing to do"
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
  if(write_selected_subcells) then
     icell_base_for_subs => ap.ivec.'-subs'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(icell_base_for_subs)) then
        write(*,*) "ERROR: for some reason, there is no list of base cell indices returned by argument parser, "//&
             "even though there was no error parsing argument -subs. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
     ncell_base_for_subs = size(icell_base_for_subs)
  end if ! write_selected_subcells
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
  ! this executable assumes inversion grid type chunksInversionGrid
  call createInversionGrid(invgrid,'chunksInversionGrid',inpar_inversion_grid,invgrid_path,11,errmsg,&
       recreate=recreate_invgrid)
  if(.level.errmsg/=0) call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
  call getInvgrid(invgrid,chunks_invgrid)
!
  if(write_selected_subcells) then
     if(any(icell_base_for_subs<1 .or. icell_base_for_subs>.nbase.chunks_invgrid)) then
        write(*,*) "ERROR: among the base cell indices requested by option -subs, there are",&
             count(icell_base_for_subs<1 .or. icell_base_for_subs>.nbase.chunks_invgrid),&
             " invalid indices which are either < 1 or greater ",&
             "than the maximum number of base cells of this chunks inversion grid = ",.nbase.chunks_invgrid
        goto 1
     end if
  end if ! write_selected_subcells
!
  if(write_all_subcells) then
     ncell_base_for_subs = .nbase.chunks_invgrid
     allocate(icell_base_for_subs(ncell_base_for_subs))
     icell_base_for_subs = (/ (j,j=1,ncell_base_for_subs) /)
  end if ! write_all_subcells
!
  if(is_refined(chunks_invgrid)) then
     write(*,*) "successfully created chunks inversion grid with ",.nbase.chunks_invgrid," base cells and ",&
          .ncell.invgrid," inversion grid cells after base cell refinement"
  else
     write(*,*) "successfully created chunks inversion grid with ",.ncell.invgrid," inversion grid cells and NO "//&
          "cell refinement (i.e. only base cells)"
  end if
  write(*,*) ""
  write(*,*) "all produced output vtk files will have basename '"//trim(outbase)//"'"
  write(*,*) ""
!
!------------------------------------------------------------------------
!  produce vtk output files
!
  stop_program = .false.
!
  if(write_base_cells_icell_base) call write_vtk_base_cells_with_ibase()
  if(stop_program) goto 1
!
  if(write_base_cells_ichunk) call write_vtk_base_cells_with_ichunk()
  if(stop_program) goto 1
!
  if(write_cells_ichunk) call write_vtk_external_cells_with_ichunk()
  if(stop_program) goto 1
!
  if(write_selected_subcells .or. write_all_subcells) call write_vtk_subcell_files_for_base_cells()
  if(stop_program) goto 1
!
!------------------------------------------------------------------------
!  clean up
!
  write(*,*) "good bye"
!
1 if(associated(icell_base_for_subs)) deallocate(icell_base_for_subs)
  nullify(chunks_invgrid)
  call dealloc(invgrid)
  call dealloc(invgrid_vtk)
  call dealloc(errmsg)
  call dealloc(ap)
!
  stop
!
!------------------------------------------------------------------------
!
contains
!
!------------------------------------------------------------------------
!

  subroutine write_vtk_base_cells_with_ibase()
    write(*,*) "producing vtk file: base cells, assigned with base cell index"
!
    if(.not.is_refined(chunks_invgrid)) &
         write(*,*) "WARNING: as the chunks inversion grid in use has no base cell refinement, it consists of ",&
         "base cells only. Hence, this operation produces a standard inversion grid vtk file (with cell indices)"
!
    filebase = trim(outbase)//"_base_cells"
!
    call new(errmsg,myname)
    call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title='base cell index on chunks_inversion_grid base cells',&
         cell_type_requested='BASE_CELLS')
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    call new(errmsg,myname)
    call writeData(invgrid_vtk,11,(/(real(j),j=1,.nbase.chunks_invgrid)/),errmsg,&
         data_name='base_cell_index',overwrite=overwrite_output)
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)  
    call dealloc(invgrid_vtk)
!
    ! clean up
1   if(associated(idx)) deallocate(idx)
    call dealloc(invgrid_vtk)
    call dealloc(errmsg)
    write(*,*) ""
    return
!
2   stop_program = .true.
    goto 1
  end subroutine write_vtk_base_cells_with_ibase
!
!------------------------------------------------------------------------
!
  subroutine write_vtk_base_cells_with_ichunk()
    write(*,*) "producing vtk file: base cells, assigned with chunk index"
!
    if(.not.is_refined(chunks_invgrid)) &
         write(*,*) "WARNING: as the chunks inversion grid in use has no base cell refinement, it consists of ",&
         "base cells only. Hence, this operation produces the same file as option -ichunk"
!
    filebase = trim(outbase)//"_ichunk_on_base_cells"
!
    call new(errmsg,myname)
    call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title='ichunk on chunks_inversion_grid base cells',&
         cell_type_requested='BASE_CELLS')
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
!
    nullify(idx)
!write(*,*) "###1 .nbase.chunks_invgrid = ",.nbase.chunks_invgrid
    idx => getChunkIndexOfBaseCellsChunksInversionGrid(chunks_invgrid,(/(j,j=1,.nbase.chunks_invgrid)/))
    if(.not.associated(idx)) then
       write(*,*) "ERROR: no chunk indices returned for ",.nbase.chunks_invgrid," base cell indices: ",&
            "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
       goto 2
    end if
    if(size(idx) /= .nbase.chunks_invgrid) then
       write(*,*) "ERROR: there were ",size(idx)," chunk indices returned for ",.nbase.chunks_invgrid,&
            " base cell indices: ",&
            "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
       goto 2
    end if
    call new(errmsg,myname)
    call writeData(invgrid_vtk,11,real(idx),errmsg,data_name='chunk_index',overwrite=overwrite_output)
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)  
    call dealloc(invgrid_vtk)
!
    ! clean up
1   if(associated(idx)) deallocate(idx)
    call dealloc(invgrid_vtk)
    call dealloc(errmsg)
    write(*,*) ""
    return
!
2   stop_program = .true.
    goto 1
  end subroutine write_vtk_base_cells_with_ichunk
!
!------------------------------------------------------------------------
!
  subroutine write_vtk_external_cells_with_ichunk()
    write(*,*) "producing vtk file: inversion grid cells, assigned with chunk index"
!
    filebase = trim(outbase)//"_ichunk"
!
    call new(errmsg,myname)
    call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title='chunk index on chunks_inversion_grid cells')
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    nullify(idx)
    idx => getChunkIndexOfExternalCellsChunksInversionGrid(chunks_invgrid,(/(j,j=1,.ncell.invgrid)/))
    if(.not.associated(idx)) then
       write(*,*) "ERROR: no chunk indices returned for ",.ncell.invgrid," inversion grid cell indices: ",&
            "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
       goto 2
    end if
    if(size(idx) /= .ncell.invgrid) then
       write(*,*) "ERROR: there were ",size(idx)," chunk indices returned for ",.ncell.invgrid,&
            " inversion grid cell indices: ",&
            "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
       goto 2
    end if
    call new(errmsg,myname)
    call writeData(invgrid_vtk,11,real(idx),errmsg,data_name='chunk_index',overwrite=overwrite_output)
    if(.level.errmsg/=0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)  
    call dealloc(invgrid_vtk)
!
    ! clean up
1   if(associated(idx)) deallocate(idx)
    call dealloc(invgrid_vtk)
    call dealloc(errmsg)
    write(*,*) ""
    return
!
2   stop_program = .true.
    goto 1
  end subroutine write_vtk_external_cells_with_ichunk
!
!------------------------------------------------------------------------
!
  subroutine write_vtk_subcell_files_for_base_cells()
    write(*,*) "producing one vtk file for ",ncell_base_for_subs," (selected) base cells: subcells assigned with ",&
         "cell index (if base cell is not refined, vtk file contains base cell only)"
!
    do j = 1,ncell_base_for_subs
       icell_base = icell_base_for_subs(j)
       write(filebase,"(a,i6.6)") trim(outbase)//"_subcells_",icell_base
!
       idx => getSubcellsOfBaseCellChunksInversionGrid(chunks_invgrid,icell_base)
       if(associated(idx)) then
          ! this base cell is refined
          write(*,*) "   base cell index ",icell_base," , ",size(idx)," subcells"
!
          call new(errmsg,myname)
          write(vtk_title,"('subcells of base cell',i6)") icell_base
          call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title=vtk_title,cell_indx_req=idx)
          if(.level.errmsg/=0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
!
          call writeData(invgrid_vtk,11,real(idx),errmsg,data_name='cell_index',&
               overwrite=overwrite_output)
          if(.level.errmsg/=0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)  
          call dealloc(invgrid_vtk)
!          
          deallocate(idx)
!
       else ! associated(idx)
          ! this base cell is not refined, produce vtk file containing base cell only
          write(*,*) "   base cell index ",icell_base," , no subcells"
!
          call new(errmsg,myname)
          write(vtk_title,"('unrefined base cell',i6)") icell_base
          call init(invgrid_vtk,invgrid,filebase,vtk_format,errmsg,vtk_title=vtk_title,&
               cell_indx_req=(/icell_base/),cell_type_requested='BASE_CELLS')
          if(.level.errmsg/=0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
!
          ! get external cell index of this base cell
          idx => getExternalCellIndexOfBaseCellsChunksInversionGrid(chunks_invgrid,(/icell_base/))
          if(.not.associated(idx)) then
             write(*,*) "ERROR: no external cell index returned for base cell ",icell_base," (out of max. ",&
                  .nbase.chunks_invgrid,"): ",&
                  "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
             goto 2
          end if
          if(size(idx) /= 1) then
             write(*,*) "ERROR: there were ",size(idx)," external cell indices returned for 1 base cell index: ",&
                  "THIS ERROR SHOULD NOT OCCUR, MODULE chunksInversionGrid SEEMS INCONSISTENT"
             goto 2
          end if
!
          call writeData(invgrid_vtk,11,(/real(idx(1))/),errmsg,data_name='cell_index',&
               overwrite=overwrite_output)
          if(.level.errmsg/=0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)  
          call dealloc(invgrid_vtk)
!          
          deallocate(idx)
!          
       end if ! associated(idx)
    end do ! j
!
    ! clean up
1   if(associated(idx)) deallocate(idx)
    call dealloc(invgrid_vtk)
    call dealloc(errmsg)
    write(*,*) ""
    return
!
2   stop_program = .true.
    goto 1
  end subroutine write_vtk_subcell_files_for_base_cells
!
end program chunksInvgrid2vtk
