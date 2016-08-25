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
program createShorelines

  use inversionGrid
  use inputParameter
  use vectorPointer
  use realloc
  use argumentParser
  use fileUnitHandler
  use errorMessage

  implicit none

  ! argument parser
  type (argument_parser) :: ap
  character(len=400) :: main_parfile,outfile_base,filename_GSHHS_bin

  ! parameter files
  type (input_parameter) :: main_inpar
  character (len=80), dimension(5) :: main_inpar_keys
  data main_inpar_keys/'CURRENT_ITERATION_STEP', 'DEFAULT_VTK_FILE_FORMAT','ITERATION_STEP_PATH', &
       'MAIN_PATH_INVERSION', 'PARFILE_ITERATION_STEP'/
  integer :: i_iter
  character(len=350) :: iter_path
  type (input_parameter) :: iter_inpar
  character (len=80), dimension(2) :: iter_inpar_keys
  data iter_inpar_keys/'TYPE_INVERSION_GRID', 'PARFILE_INVERSION_GRID'/
  type (input_parameter) :: invgrid_inpar
  character (len=80), dimension(1) :: invgrid_inpar_keys
  data invgrid_inpar_keys/'VTK_PROJECTION'/

  ! inversion grid
  type (inversion_grid) :: invgrid
  character(len=350) :: invgrid_type,vtk_projection

  ! gshhs file
  integer :: gshhs_fu,npol,np_pol,np_seg,np_seg_vtk,ipoint,ichunk_now,ichunk_previous
  integer, parameter :: points_alloc_chunk = 1000
  integer, parameter :: lines_alloc_chunk = 1200
  integer, dimension(11) :: pol_header
  integer, dimension(:,:), allocatable :: lonlat_pol
  real, dimension(:), allocatable :: plon_seg,plat_seg,plon_vtk,plat_vtk
  real, dimension(:), allocatable :: c1,c2,c3

  ! vtk file
  integer :: np_vtk,nlines_vtk,size_lines_vtk,vtk_fu,ios
  logical :: file_exists,vtk_file_is_ascii,err
  character(len=403) :: vtk_filename
  character (len=7) :: open_status
  character (len=500) :: char_string
  real, dimension(:,:), pointer :: points
  integer, dimension(:), pointer :: lines
  character (len=1), parameter :: eol_char = char(10)

  ! other
  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=16) :: myname = 'createShorelines'
  integer :: j
  logical :: stop_program

  nullify(points,lines)
!
!------------------------------------------------------------------------
!  preliminary processing
!
  stop_program = .false.
!
  call init(ap,myname,'create shorelines vtk file with data from native binary GSHHS file, projected on '//&
       'current inversion grid surface')
  ! define positional arguments
  call addPosarg(ap,'filename_GSHHS_bin','sval','GSHHS shore line file (native binary format)')
  call addPosarg(ap,'outfile_base','sval','base name of output file')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
  ! get values of positional arguments
  filename_GSHHS_bin = ap.sval.'filename_GSHHS_bin'
  outfile_base = ap.sval.'outfile_base'
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) goto 1
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
  write(*,*) ""
  write(*,*) "WELCOME TO GSHHS SHORE LINE CREATION FOR ASKI!"
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
  write(*,*) "   successfully read main parameter file '",trim(main_parfile),"'"
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
  case('BINARY')
     vtk_file_is_ascii = .false.
  case('ASCII')
     vtk_file_is_ascii = .true.
  case default
     write(*,*) "ERROR: parameter 'DEFAULT_VTK_FILE_FORMAT' = '"//trim(main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT')//&
          "' of main parfile is not valid; must be either 'ASCII' or 'BINARY'"
     stop_program = .true.
  end select
  write(*,*) "      current iteration step = ",i_iter
  write(*,*) "      default vtk file format (will be used here) = ",trim(main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT')
!
  ! read in iteraton step parameter file (only those keys which are needed here)
  call new(errmsg,myname)
  call createKeywordsInputParameter(iter_inpar,iter_inpar_keys)
  call readSubroutineInputParameter(iter_inpar,get(fuh),trim(iter_path)//&
       trim(main_inpar.sval.'PARFILE_ITERATION_STEP'),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "   successfully read iteration step parameter file '",trim(iter_path)//&
       trim(main_inpar.sval.'PARFILE_ITERATION_STEP'),"'"
!
  invgrid_type = iter_inpar.sval.'TYPE_INVERSION_GRID'
  write(*,*) "      inversion grid type = ",trim(invgrid_type)
!
  ! read in inversion grid parameter file (only those keys which are needed here)
  call new(errmsg,myname)
  call createKeywordsInputParameter(invgrid_inpar,invgrid_inpar_keys)
  call readSubroutineInputParameter(invgrid_inpar,get(fuh),trim(iter_path)//&
       trim(iter_inpar.sval.'PARFILE_INVERSION_GRID'),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "   successfully read inversion grid parameter file '",trim(iter_path)//&
       trim(iter_inpar.sval.'PARFILE_INVERSION_GRID'),"'"
!
  vtk_projection = invgrid_inpar.sval.'VTK_PROJECTION'
!
  ! create inversion grid
  call new(errmsg,myname)
  call createInversionGrid(invgrid,invgrid_type,&
       trim(iter_path)//trim(iter_inpar.sval.'PARFILE_INVERSION_GRID'),iter_path,&
       get(fuh),errmsg,recreate=.false.)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "      vtk projection type of inversion grid = ",trim(vtk_projection)
!
!------------------------------------------------------------------------
!  read shorelines, select (parts of) segments, transform to vtk inversion grid coords, write vtk file
!
  ! first estimate a priori boundaries of the target region w,s,e,n from the centers (and width(s)?!) 
  ! of the given inversion grid

  if(stop_program) goto 1
  write(*,*) ""
!
  write(*,*) "READING GSHHS SHORELINE FILE '"//trim(filename_GSHHS_bin)//&
       "' NOW, SELECTING LINE SEGMENTS IN RANGE OF CURRENT INVERSION GRID"
  call getShoreLinesNativeBinaryGSHHS()
  if(stop_program) goto 1
  write(*,*) ""

  write(*,*) "WRITING SHORELINE VTK FILE (basename '"//trim(outfile_base)//"') NOW"
  call writeVtkFileShoreLines()
  if(stop_program) goto 1
!
!------------------------------------------------------------------------
!  clean up
!
1 call dealloc(errmsg)
  call dealloc(ap)
  call dealloc(fuh)
  if(allocated(lonlat_pol)) deallocate(lonlat_pol)
  if(allocated(plon_seg)) deallocate(plon_seg)
  if(allocated(plat_seg)) deallocate(plat_seg)
  if(allocated(plon_vtk)) deallocate(plon_vtk)
  if(allocated(plat_vtk)) deallocate(plat_vtk)
  if(associated(points)) deallocate(points)
  if(associated(lines)) deallocate(lines)
  if(allocated(c1)) deallocate(c1)
  if(allocated(c2)) deallocate(c2)
  if(allocated(c3)) deallocate(c3)

  write(*,*) ""
  write(*,*) "Good Bye"
  write(*,*) ""

  stop

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine getShoreLinesNativeBinaryGSHHS()
  ! open shore line file (native binary format) to read
  gshhs_fu = get(fuh)
  open(unit=gshhs_fu,file=trim(filename_GSHHS_bin),form='UNFORMATTED',access='STREAM',action='READ',&
       convert='BIG_ENDIAN',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR in getShoreLinesNativeBinaryGSHHS: could not open shore line native binary file '"//&
          trim(filename_GSHHS_bin)//"'"
     stop_program = .true.
     goto 1
  endif

  ! allocate arrays
  allocate(points(3,points_alloc_chunk),lines(lines_alloc_chunk))
  np_vtk = 0
  nlines_vtk = 0
  size_lines_vtk = 0

  ! read in file, looping on all contained polygons

  npol = 0
  do while(.true.)
     ! read header of current polygon (if there is any)
     read(gshhs_fu,iostat=ios) pol_header
     if(ios /= 0) exit

     if(pol_header(1)/=npol) then
        ! first entry of header does not coincide with polygon ID = index starting from 0
        write(*,*) "ERROR in getShoreLinesNativeBinaryGSHHS: the header of the ",npol+1,&
             "'th polygon says: polygon ID = ",pol_header(1),"; expectig value ",npol," => unexpected file format!"
        stop_program = .true.
        goto 1
     end if

     ! second entry of header is number of points in polygon
     np_pol = pol_header(2)
     if(np_pol <= 0) then
        write(*,*) "ERROR in getShoreLinesNativeBinaryGSHHS: the header of the ",npol,&
             "'th polygon says: number of points in polygon = ",np_pol,&
             "; expecting number > 0 => unexpected file format!"
        stop_program = .true.
        goto 1
     end if

     ! read np_pol lat/lon pairs from file
     if(allocated(lonlat_pol)) deallocate(lonlat_pol)
     allocate(lonlat_pol(2,np_pol))
     read(gshhs_fu,iostat=ios) lonlat_pol
     if(ios /= 0) then
        write(*,*) "ERROR in getShoreLinesNativeBinaryGSHHS: reading the ",npol,&
             "'th polygon points: could not read ",np_pol," lon/lat pairs from file (as expected by polygon header) ",&
             "file ended in the middle of this polygon => unexpected file format!"
        stop_program = .true.
        goto 1
     end if
     ! memorize that the npol'th polygon was read from file
     npol = npol + 1

     call selectSplitPolygons()

  end do ! loop on content of file

  write(*,*) "   found ",npol," polygons in native binary GSHHS file"
  if(npol>0) then
     write(*,*) "   from these, there was a total of ",np_vtk," points selected to be inside the inversion grid and a total of ",&
          nlines_vtk," line segments that will be written to vtk (possible some polygons were split up/truncated/removed ",&
          "in order to fit inside the inversion grid domain and not to intersect chunk boundaries in an unwanted manner)"
  end if

1 close(gshhs_fu)
  call add(fuh,gshhs_fu)
end subroutine getShoreLinesNativeBinaryGSHHS
!
!------------------------------------------------------------------------
!
subroutine selectSplitPolygons()
  ! loop on all points of this polygon, check whether points are inside the inversion grid, or polygon needs to 
  ! be split up into more than one vtk line segment
  if(allocated(plat_seg)) deallocate(plat_seg)
  if(allocated(plon_seg)) deallocate(plon_seg)
  if(allocated(plat_vtk)) deallocate(plat_vtk)
  if(allocated(plon_vtk)) deallocate(plon_vtk)
  allocate(plat_seg(np_pol),plon_seg(np_pol),plat_vtk(np_pol),plon_vtk(np_pol))
  np_seg = 0
  np_seg_vtk = 0
  ichunk_now = -1
  do ipoint = 1,np_pol
     ichunk_previous = ichunk_now

     ! the np_pol lon/lat pairs read from file represent micro degrees, so multiply by 1.e-6 here to get degrees
     if(pointInsideInversionGrid(invgrid,1.e-6*lonlat_pol(2,ipoint),1.e-6*lonlat_pol(1,ipoint),0.0,&
          'station',ichunk=ichunk_now)) then
        ! YES, point is inside the inverison grid volume:

        ! if this is not the first point of the vtk segment inside the inversion grid, check if it
        ! crosses a chunk border which is discontinuous in the current inversion grid's vtk projection
        ! so that the vtk segment needs to be split here
        if(np_seg >= 1) then
           if(splitVtkSegment()) then
              ! if vtk segment needs to be splitted here, check whether there was a sensible vtk line found 
              ! so far (at least 2 points)
              if(np_seg >= 2) then
                 ! if there was a sensible vtk line found so far, memorize it
                 plat_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plat_seg(1:np_seg)
                 plon_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plon_seg(1:np_seg)
                 ! store line connectivity in array lines, reallocate it if too small to add this segment
                 if(size_lines_vtk+np_seg+1 > size(lines)) lines => reallocate(lines,size(lines)+lines_alloc_chunk)
                 lines(size_lines_vtk+1) = np_seg
                 ! in global counting of vtk points, the next index is np_vtk+np_seg_vtk+1
                 lines(size_lines_vtk+2:size_lines_vtk+1+np_seg) = &
                      (/ (j-1 , j=np_vtk+np_seg_vtk+1 , np_vtk+np_seg_vtk+np_seg) /) 
                 ! increase counters
                 size_lines_vtk = size_lines_vtk + 1 + np_seg
                 nlines_vtk = nlines_vtk + 1
                 np_seg_vtk = np_seg_vtk+np_seg
              end if ! np_seg >= 2
              ! create a new empty segment in which the current point will be stored (as first point)
              np_seg = 0 ! simply point to the first entry of arrays plon_seg,plat_seg, re-use the same arrays
           end if ! split_vtk_segment
        end if ! np_seg >= 1

        ! store point in the current vtk segment (existing one or the one just newly created)
        np_seg = np_seg + 1
        plat_seg(np_seg) = 1.e-6*lonlat_pol(2,ipoint)
        plon_seg(np_seg) = 1.e-6*lonlat_pol(1,ipoint)

     else ! pointInsideInversionGrid
        ! NO, point is not inside the inversion grid volume:

        ! check whether there was a sensible vtk line found so far (at least 2 points)
        if(np_seg >= 2) then
           ! if there was a sensible vtk line found so far, memorize it
           plat_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plat_seg(1:np_seg)
           plon_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plon_seg(1:np_seg)
           ! store line connectivity in array lines, reallocate it if too small to add this segment
           if(size_lines_vtk+np_seg+1 > size(lines)) lines => reallocate(lines,size(lines)+lines_alloc_chunk)
           lines(size_lines_vtk+1) = np_seg
           ! in global counting of vtk points, the next index is np_vtk+np_seg_vtk+1
           lines(size_lines_vtk+2:size_lines_vtk+1+np_seg) = &
                (/ (j-1 , j=np_vtk+np_seg_vtk+1 , np_vtk+np_seg_vtk+np_seg) /) 
           ! increase counters
           size_lines_vtk = size_lines_vtk + 1 + np_seg
           nlines_vtk = nlines_vtk + 1
           np_seg_vtk = np_seg_vtk+np_seg
        end if ! np_seg >= 2

        ! create a new empty vtk segment for rest of the GSHHS segment (possibly GSHHS segment enters inversion grid again)
        np_seg = 0

        ! in any case, make the current chunk index invalid (this point is not inside the inversion grid)
        ichunk_now = -1

     end if ! pointInsideInversionGrid

  end do ! ipoint

  ! after looping on all points of the segment, memorize (the tail part of) it for vtk output
  ! ONLY IN CASE there are at least 2 points in it.
  ! this could be only the tail part of the segment, if some segment points were outside the inversion grid
  if(np_seg >= 2) then
     ! if there was a sensible vtk line found so far, memorize it
     plat_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plat_seg(1:np_seg)
     plon_vtk(np_seg_vtk+1:np_seg_vtk+np_seg) = plon_seg(1:np_seg)
     ! store line connectivity in array lines, reallocate it if too small to add this segment
     if(size_lines_vtk+np_seg+1 > size(lines)) lines => reallocate(lines,size(lines)+lines_alloc_chunk)
     lines(size_lines_vtk+1) = np_seg
     ! in global counting of vtk points, the next index is np_vtk+np_seg_vtk+1
     lines(size_lines_vtk+2:size_lines_vtk+1+np_seg) = &
          (/ (j-1 , j=np_vtk+np_seg_vtk+1 , np_vtk+np_seg_vtk+np_seg) /) 
     ! increase counters
     size_lines_vtk = size_lines_vtk + 1 + np_seg
     nlines_vtk = nlines_vtk + 1
     np_seg_vtk = np_seg_vtk+np_seg
  end if ! np_seg >= 2

  ! finally, if any points were found above, transform them from plat_vtk,plon_vtk to inversion grid vtk 
  ! projection and add them to array points
  if(np_seg_vtk > 0) then
     if(allocated(c1)) deallocate(c1)
     if(allocated(c2)) deallocate(c2)
     if(allocated(c3)) deallocate(c3)
     allocate(c1(np_seg_vtk),c2(np_seg_vtk),c3(np_seg_vtk))
     c1 = plat_vtk(1:np_seg_vtk)
     c2 = plon_vtk(1:np_seg_vtk)
     c3 = 0.0
     call transformToVtkInversionGrid(invgrid,c1,c2,c3,'station',errmsg)

     if(np_vtk+np_seg_vtk > size(points,2)) points => reallocate(points,3,size(points,2)+points_alloc_chunk)
     points(1,np_vtk+1:np_vtk+np_seg_vtk) = c1
     points(2,np_vtk+1:np_vtk+np_seg_vtk) = c2
     points(3,np_vtk+1:np_vtk+np_seg_vtk) = c3
     np_vtk = np_vtk + np_seg_vtk
  end if ! np_seg_vtk > 0
end subroutine selectSplitPolygons
!
!------------------------------------------------------------------------
!
function splitVtkSegment() result(split)
  logical :: split
!
  ! initiate default on return: do not split
  split = .false.
!
  if(ichunk_previous <= 0 .or. ichunk_now <= 0) return
  if(ichunk_previous == ichunk_now) return
!
  ! in the following, assume that ichunk_previous /= ichunk_now and that both are > 0
  select case(invgrid_type)
  case('chunksInversionGrid')
     select case(vtk_projection)
     case('LOCAL_FLAT','LOCAL_NORTH_FLAT')
        !
        ! DISTRIBUTION OF CHUNKS IN chunksInversionGrid (LOCAL_FLAT projection here):
        !
        !     +---+
        !     | 3 |
        ! +---+---+---+---+
        ! | 2 | 1 | 5 | 4 |
        ! +---+---+---+---+
        !     | 6 |
        !     +---+
        !
        select case(ichunk_previous)
        case(1)
           ! only split, if ichunk_now == 4
           split = ichunk_now == 4
        case(2,3,6)
           ! split unless ichunk_now == 1
           split = ichunk_now /= 1
        case(4)
           ! split unless ichunk_now == 5
           split = ichunk_now /= 5
        case(5)
           ! split unless ichunk_now == 1 or ichunk_now == 4
           split = (ichunk_now /= 1).or.(ichunk_now /= 4)
        end select ! ichunk_previous
     case default
        return
     end select ! case(vtk_projection)
  case default
     return
  end select ! case(invgrid_type)
end function splitVtkSegment
!
!------------------------------------------------------------------------
!
subroutine writeVtkFileShoreLines()
  if(np_vtk == 0 .or. size_lines_vtk == 0) then
     write(*,*) "ERROR in writeVtkFileShoreLines: no vtk points or lines defined"
     stop_program = .true.
     return
  end if

  ! OPEN VTK FILE
  vtk_filename = trim(outfile_base)//'.vtk'
  write(*,*) "   output vtk file will be '"//trim(vtk_filename)//"'"

  ! check if file exists, set open_status as required
  inquire(file=vtk_filename,exist=file_exists)
  if(file_exists) then
     open_status = 'REPLACE'
     write(*,*) "   output vtk file exists and will be overwritten"
  else ! file_exists
     open_status = 'NEW'
     write(*,*) "   output vtk file does not exist yet and will newly created"
  endif ! file_exists

  ! BINARY OR ASCII FILE?
  vtk_file_is_ascii = ( main_inpar.sval.'DEFAULT_VTK_FILE_FORMAT') == 'ASCII'
  if(vtk_file_is_ascii) then
     write(*,*) "   output vtk file will be ASCII"
  else
     write(*,*) "   output vtk file will be BINARY"
  end if

  vtk_fu = get(fuh)
  if(vtk_file_is_ascii) then
     open(unit=vtk_fu,file=trim(vtk_filename),form='FORMATTED',status=trim(open_status),action='WRITE',iostat=ios)
     if(ios/=0) then
        write(*,*) "ERROR in writeVtkFileShoreLines: could not open vtk output file to write, iostat = ",ios
        stop_program = .true.
        return
     end if
  else
     open(unit=vtk_fu,file=trim(vtk_filename),form='UNFORMATTED',access='STREAM',status=trim(open_status),action='WRITE',&
          convert='BIG_ENDIAN',iostat=ios)
     if(ios/=0) then
        write(*,*) "ERROR in writeVtkFileShoreLines: could not open vtk output file to write, iostat = ",ios
        stop_program = .true.
        return
     end if
  end if

  ! WRITE vtk FILE HEADER

  ! remember with err if there was an error somewhere
  err = .false.
  if(vtk_file_is_ascii) then
     write(unit=vtk_fu,fmt='(a)',iostat=ios) '# vtk DataFile Version 3.1'  ; err = err.or.(ios/=0)
     write(unit=vtk_fu,fmt='(a)',iostat=ios) 'GSHHS Shorelines, projected / cut to inversion grid' ; err = err.or.(ios/=0)
     write(unit=vtk_fu,fmt='(a)',iostat=ios) 'ASCII'                       ; err = err.or.(ios/=0)
     write(unit=vtk_fu,fmt='(a)',iostat=ios) 'DATASET POLYDATA'            ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing vtk Header"
        stop_program = .true.
        return
     endif
  else ! vtk_file_is_ascii
     write(unit=vtk_fu,iostat=ios) '# vtk DataFile Version 3.1'//eol_char  ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) 'GSHHS Shorelines, projected / cut to inversion grid'//eol_char ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) 'BINARY'//eol_char                      ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) 'DATASET POLYDATA'//eol_char            ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing vtk Header"
        stop_program = .true.
        return
     endif
  endif ! vtk_file_is_ascii

  ! WRITE POINTS TO VTK FILE

  ! remember with err if there was an error somewhere
  err = .false.
  if(vtk_file_is_ascii) then
     write(unit=vtk_fu,fmt='(a,i12,a)',iostat=ios) 'POINTS ',np_vtk,' float'  ; err = err.or.(ios/=0)
     write(unit=vtk_fu,fmt='(3e14.6e2)',iostat=ios) points(:,1:np_vtk)        ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing POINTS"
        stop_program = .true.
        return
     endif
  else ! vtk_file_is_ascii
     write(char_string,'(a,i12,a)',iostat=ios) 'POINTS ',np_vtk,' float'
     write(unit=vtk_fu,iostat=ios) trim(char_string)//eol_char                  ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) points(:,1:np_vtk)                      ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) eol_char                                ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing POINTS"
        stop_program = .true.
        return
     endif
  endif ! vtk_file_is_ascii

  ! WRITE LINES TO VTK FILE

  ! remember with err if there was an error somewhere
  err = .false.
  if(vtk_file_is_ascii) then
     ! LINES
     write(unit=vtk_fu,fmt='(a,2i12)',iostat=ios)'LINES ',nlines_vtk,size_lines_vtk   ; err = err.or.(ios/=0)
     write(unit=vtk_fu,fmt='(i12)',iostat=ios) lines(1:size_lines_vtk)               ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing LINES"
        stop_program = .true.
        return
     endif
  else ! vtk_file_is_ascii
     ! LINES
     write(char_string,'(a,2i12)',iostat=ios)'LINES ',nlines_vtk,size_lines_vtk     ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) trim(char_string)//eol_char                      ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) lines(1:size_lines_vtk)                     ; err = err.or.(ios/=0)
     write(unit=vtk_fu,iostat=ios) eol_char                                    ; err = err.or.(ios/=0)
     if(err) then ! if any of the above ios were /= 0
        write(*,*) "ERROR in writeVtkFileShoreLines: there was an error writing LINES"
        stop_program = .true.
        return
     endif
  endif ! vtk_file_is_ascii

  close(vtk_fu)
  call add(fuh,vtk_fu)
end subroutine writeVtkFileShoreLines

end program createShorelines
