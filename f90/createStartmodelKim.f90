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
program createStartmodelKim
  use inversionGrid
  use modelParametrization
  use kernelInvertedModel
  use argumentParser
  use string
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: outbase,model_type,model_file,&
       inversion_grid_type,inversion_grid_inpar,inversion_grid_path,str
       
  character(len=character_length_pmtrz) :: model_pmtrz

  type (error_message) :: errmsg
  character(len=19) :: myname = 'createStartmodelKim'

  type (inversion_grid) :: invgrid
  character(len=6) :: vtk_format

  type (kernel_inverted_model) :: kim

  integer :: lu
  logical :: stop_after_command_line

!------------------------------------------------------------------------
!  preliminary processing
!
   call init(ap,myname,"Create a file of type kernel_inverted_model (.kim) containing pre-defined values on the "//&
       "inversion grid. Can be used to create a start model for the full waveform inversion process.")
   call addOption(ap,"-igtype",.true.,"(mandatory) TYPE_INVERSION_GRID as in ASKI iteration step parameter file",&
        "sval","")
   call addOption(ap,"-igpar",.true.,"(mandatory) PARFILE_INVERSION_GRID as in ASKI iteration step parameter file",&
        "sval","")
   call addOption(ap,"-igpath",.true.,"(mandatory) is treated as current iteration step path, used for some "//&
        "inversion grids to write/read own files","sval","")
   call addOption(ap,"-mpmtrz",.true.,"(mandatory) model parametrization of the model which is to be created "//&
        "(must be consistent with content of model input file, see -mfile)","sval","")
   call addOption(ap,"-mtype",.true.,"(mandatory) Defines type of model file and the interpolation type. "//&
        "Supported types: 1D_linear, 3D_structured","sval","")
   call addOption(ap,"-mfile",.true.,"(mandatory) model input file of type defined by -mtype","sval","")
   call addOption(ap,"-o",.true.,"(optional) output base name","sval","start_model")
   call addOption(ap,"-bin",.false.,"if set, the output vtk files will be binary, otherwise they will be ascii")
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
  ! inversion_grid_type
  if(ap.optset.'-igtype') then
     inversion_grid_type = ap.sval.'-igtype'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igtype"
  end if
!
  ! inversion_grid_inpar
  if(ap.optset.'-igpar') then
     inversion_grid_inpar = ap.sval.'-igpar'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpar"
  end if
!
  ! inversion_grid_path
  if(ap.optset.'-igpath') then
     inversion_grid_path = ap.sval.'-igpath'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -igpath"
  end if
!
  ! model_pmtrz
  if(ap.optset.'-mpmtrz') then
     str = ap.sval.'-mpmtrz'
     model_pmtrz = str
     if(.not.validModelParametrization(model_pmtrz)) then
        stop_after_command_line = .true.
        write(*,*) "ERROR: model parametrization '"//trim(model_pmtrz)//"' read from -mpmtrz input string "//&
             "is not valid, valid parametrizations (parameters) are "//all_valid_pmtrz_param
     end if
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -mpmtrz"
  end if
!
  ! model_type
  if(ap.optset.'-mtype') then
     model_type = ap.sval.'-mtype'
     select case(model_type)
        case('1D_linear','3D_structured')
           ! ok, do nothing
        case default
           stop_after_command_line = .true.
           write(*,*) "ERROR: model type '"//trim(model_type)//"' read from -mtype input string "//&
                "is not supported, supported types are '1D_linear'and '3D_structured'"
     end select
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -mtype"
  end if
!
  ! model_file
  if(ap.optset.'-mfile') then
     model_file = ap.sval.'-mfile'
  else
     stop_after_command_line = .true.
     write(*,*) "ERROR: please indicate -mfile"
  end if
!
  ! outbase
  outbase = ap.sval.'-o'
!
  if(ap.optset.'-bin') then
     vtk_format = 'BINARY'
  else
     vtk_format = 'ASCII'
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
  lu = 11
!
!------------------------------------------------------------------------
!  create inversion grid
!
  call new(errmsg,myname)
  call createInversionGrid(invgrid,inversion_grid_type,inversion_grid_inpar,inversion_grid_path,lu,errmsg)
  !if(.level.errmsg/=0) call print(errmsg)
  call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  create kim of current model type from model file
!
  call new(errmsg,myname)
  select case(model_type)
  case('1D_linear')
     call createStartmodel1DLinearKim(model_file,lu,invgrid,model_pmtrz,kim,errmsg)
  case('3D_structured')
     call createStartmodel3DstructuredKim(model_file,lu,invgrid,model_pmtrz,kim,errmsg)
  end select
  !if(.level.errmsg/=0) call print(errmsg)
  call print(errmsg)
  if(.level.errmsg==2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  write kim and vtk files
!
  call new(errmsg,myname)
  call writeFileKernelInvertedModel(kim,trim(outbase)//'.kim',lu,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call new(errmsg,myname)
  call writeVtkKernelInvertedModel(kim,invgrid,vtk_format,outbase,lu,errmsg,overwrite=.true.)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "successfully written all output files with base name '"//trim(outbase)//"'"
!
!------------------------------------------------------------------------
!  clean up
!
1 call dealloc(errmsg)
  call dealloc(invgrid)
  call dealloc(kim)
  call dealloc(ap)
!
!------------------------------------------------------------------------
!  subprograms
!
contains
!
  subroutine createStartmodel1DLinearKim(model_file,lu,invgrid,model_pmtrz,kim,errmsg)
    use inversionGrid
    use kernelInvertedModel
    use modelParametrization,only : character_length_param
    use errorMessage
    implicit none
!
    character(len=*) :: model_file,model_pmtrz
    integer :: lu
    type (inversion_grid) :: invgrid
    type (kernel_inverted_model) :: kim
    type (error_message) :: errmsg
    ! local
    character(len=27) :: myname = 'createStartmodel1DLinearKim'
    type (error_message) :: errmsg2
    integer :: icell,iparam,i
    real, dimension(3) :: cc
    ! kim
    integer, dimension(:), allocatable :: indx_kim
    real, dimension(:,:), allocatable :: val_kim
    integer :: nval_kim
    real :: w1,w2
    ! 1D model file
    integer :: nparam,icoord,nval
    character(len=character_length_param), dimension(:), pointer :: param
    real, dimension(:), pointer :: coords
    real, dimension(:,:), pointer :: val
!
    nullify(param,coords,val)
    call addTrace(errmsg,myname)
!
    call read1DModelFile(model_file,lu,nparam,param,icoord,nval,coords,val,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call init(kim,model_pmtrz,errmsg)
    if(.level.errmsg == 2) goto 1
!
    allocate(indx_kim(.ncell.invgrid),val_kim(.ncell.invgrid,nparam))
!
    ! interpolate model values for all cells
    nval_kim = 0
    do icell = 1,.ncell.invgrid
       call new(errmsg2,myname)
       call getCenterCellInversionGrid(invgrid,icell,cc(1),cc(2),cc(3),errmsg2)
       if(.level.errmsg2 == 2) goto 1
       call dealloc(errmsg2)
!
       ! locate coordinate cc(icoord) of current inversion grid cell inside all interpolation coordinats
       i = locateCoordinate(coords,nval,cc(icoord))
       if(i<0) cycle
       ! compute interpolation weights
       w2 = (cc(icoord)-coords(i))/(coords(i+1)-coords(i))
       w1 = 1. - w2
!
       ! interpolate model values
       nval_kim = nval_kim + 1
       indx_kim(nval_kim) = icell
       val_kim(nval_kim,:) = w1*val(i,:) + w2*val(i+1,:)
    end do ! icell
!
    ! add model values to kernel inverted model
    if(nval_kim > 0) then
       do iparam = 1,nparam
          call addValuesKernelInvertedModel(kim,param(iparam),indx_kim(1:nval_kim),val_kim(1:nval_kim,iparam),errmsg)
          if(.level.errmsg == 2) goto 1
       end do ! iparam
    end if

1   call dealloc(errmsg2)
    if(allocated(indx_kim)) deallocate(indx_kim)
    if(allocated(val_kim)) deallocate(val_kim)
  end subroutine createStartmodel1DLinearKim
!
!------------------------------------------------------------------------
!
  subroutine createStartmodel3DstructuredKim(model_file,lu,invgrid,model_pmtrz,kim,errmsg)
    use inversionGrid
    use kernelInvertedModel
    use modelParametrization,only : character_length_param !character_length_param = 6
    use errorMessage
    implicit none
!
    character(len=*) :: model_file,model_pmtrz
    integer :: lu
    type (inversion_grid) :: invgrid
    type (kernel_inverted_model) :: kim
    type (error_message) :: errmsg
    ! local
    character(len=31) :: myname = 'createStartmodel3DstructuredKim'
    type (error_message) :: errmsg2
    integer :: icell,iparam,i,j,k
    real, dimension(3) :: cc
    ! kim
    integer, dimension(:), allocatable :: indx_kim
    real, dimension(:,:), allocatable :: val_kim
    integer :: nval_kim
    real :: w1,w2,w3,w4,w5,w6,w7,w8,x,y,z
    ! 3D_structured model file
    integer :: nparam,icoord,nval
    character(len=character_length_param), dimension(:), pointer :: param
    real, dimension(:), pointer :: Xcoord,Ycoord,Zcoord
    real, dimension(:,:,:,:), pointer :: val
    integer :: nx,ny,nz
!
    nullify(param,Xcoord,Ycoord,Zcoord,val)
    call addTrace(errmsg,myname)
!
    call read3DstructuredModelFile(model_file,lu,nparam,param,icoord,nval,Xcoord,Ycoord,Zcoord,nx,ny,nz,val,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call init(kim,model_pmtrz,errmsg)
    if(.level.errmsg == 2) goto 1
!
    allocate(indx_kim(.ncell.invgrid),val_kim(.ncell.invgrid,nparam))
!
    ! interpolate model values for all cells
    nval_kim = 0
    do icell = 1,.ncell.invgrid
       call new(errmsg2,myname)
       call getCenterCellInversionGrid(invgrid,icell,cc(1),cc(2),cc(3),errmsg2)
       if(.level.errmsg2 == 2) goto 1
       call dealloc(errmsg2)
!
       ! locate coordinate cc(icoord) of current inversion grid cell inside all interpolation coordinats
       i = locateCoordinate(Zcoord,nz,cc(3))
       j = locateCoordinate(Ycoord,ny,cc(2))
       k = locateCoordinate(Xcoord,nx,cc(1))
       if(i<0) cycle
       if(j<0) cycle
       if(k<0) cycle
       ! compute interpolation weights
        x = ((cc(1)-Xcoord(k))/(Xcoord(k+1)-Xcoord(k)))
        y = ((cc(2)-Ycoord(j))/(Ycoord(j+1)-Ycoord(j)))
        z = ((cc(3)-Zcoord(i))/(Zcoord(i+1)-Zcoord(i)))
        w1 = (1-x)*(1-y)*(1-z)
        w2 = x*(1-y)*(1-z)
        w3 = (1-x)*y*(1-z)
        w4 = (1-x)*(1-y)*z
        w5 = x*(1-y)*z
        w6 = (1-x)*y*z
        w7 = x*y*(1-z)
        w8 = x*y*z
!
       ! interpolate model values
       nval_kim = nval_kim + 1
       indx_kim(nval_kim) = icell
       val_kim(nval_kim,:) = w1*val(k,j,i,:) + w2*val(k+1,j,i,:)+w3*val(k,j+1,i,:)+w4*val(k,j,i+1,:)+ &
                             w5*val(k+1,j,i+1,:)+w6*val(k,j+1,i+1,:)+w7*val(k+1,j+1,i,:)+w8*val(k+1,j+1,i+1,:)
    end do ! icell
!
    ! add model values to kernel inverted model
    if(nval_kim > 0) then
       do iparam = 1,nparam
          call addValuesKernelInvertedModel(kim,param(iparam),indx_kim(1:nval_kim),val_kim(1:nval_kim,iparam),errmsg)
          if(.level.errmsg == 2) goto 1
       end do ! iparam
    end if

1   call dealloc(errmsg2)
    if(allocated(indx_kim)) deallocate(indx_kim)
    if(allocated(val_kim)) deallocate(val_kim)
  end subroutine createStartmodel3DstructuredKim
!
!------------------------------------------------------------------------
!
  subroutine read1DModelFile(filename,lu,nparam,param,icoord,nval,coords,val,errmsg)
    use modelParametrization,only : character_length_param
    character(len=*) :: filename
    integer :: lu,nparam,icoord,nval
    character(len=character_length_param), dimension(:), pointer :: param
    real, dimension(:), pointer :: coords
    real, dimension(:,:), pointer :: val
    type (error_message) :: errmsg
    ! local
    character(len=15) :: myname = 'read1DModelFile'
    character(len=400) :: errstr
    integer :: ncol,ios,ival
!
    call addTrace(errmsg,myname)
!
    open(unit=lu,file=filename,status='unknown',form='formatted',action='read',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "ERROR: could not open 1D model file '"//trim(filename)//"' to read"
       call add(errmsg,2,errstr,myname); return
    end if
!
    read(lu,*,iostat=ios) nval,ncol
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read number of model value lines and number of columns "//&
            "from first line of 1D model file '"//trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
    if(nval<2 .or. ncol<2) then
       close(lu)
       write(errstr,*) "ERROR: number of model value lines ",nval," and number of columns ",ncol,&
            "on first line of 1D model file '"//trim(filename)//"' must at least be 2"
       call add(errmsg,2,errstr,myname); return
    end if
!
    nparam = ncol-1
    allocate(param(nparam))
    read(lu,*,iostat=ios) icoord,param
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read index of interpolation coordinate and ",nparam,&
            " parameter names from second line of 1D model file '"//trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
    if(icoord<1 .or. icoord>3) then
       close(lu)
       write(errstr,*) "ERROR: index of interpolation coordinate ",icoord,&
            "on second line of 1D model file '"//trim(filename)//"' must be 1,2 or 3"
       call add(errmsg,2,errstr,myname); return
    end if
!
    allocate(coords(nval),val(nval,nparam))
    do ival = 1,nval
       read(lu,*,iostat=ios) coords(ival),val(ival,:)
       if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read ",nparam+1," real values from ",ival+2,&
            "'th line of 1D model file '"//trim(filename)//"' (which contains the ",ival,&
            "'th set of interpolation coordinate and model values)"
       call add(errmsg,2,errstr,myname); return
       end if
    end do ! ival
!
    close(lu)
  end subroutine read1DModelFile
!
!------------------------------------------------------------------------
!
  subroutine read3DstructuredModelFile(filename,lu,nparam,param,icoord,nval,Xcoord,Ycoord,Zcoord,nx,ny,nz,val,errmsg)
    use modelParametrization,only : character_length_param
    use string, only : countWordsString
    character(len=*) :: filename
    integer :: lu,nparam,icoord,nval
    character(len=character_length_param), dimension(:), pointer :: param
    real, dimension(:), pointer :: Xcoord, Ycoord, Zcoord
    real, dimension(:,:,:,:), pointer :: val
    integer :: nx,ny,nz
    type (error_message) :: errmsg
    ! local
    character(len=25) :: myname = 'read3DstructuredModelFile'
    character(len=400) :: errstr
    real :: minX,minY,minZ,maxX,maxY,maxZ 
    integer :: ios,i,j,k
    character(len=10) :: s
    character(len=1) :: sps
!
    call addTrace(errmsg,myname)
!
    open(unit=lu,file=filename,status='unknown',form='formatted',action='read',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "ERROR: could not open 3D_structured model file '"//trim(filename)//"' to read"
       call add(errmsg,2,errstr,myname); return
    end if
!
!   Reads the first line of the model file, number of model pionts (nx, ny, nz)
    read(lu,*,iostat=ios) nx,ny,nz
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read number of model pionts in X-, Y- & Z-direction "//&
            "from first line of 3D_structured model file '"//trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
!
!   Reads the second line of the model file, minimal coordinates (minX, minY, minZ)
    read(lu,*,iostat=ios) minX,minY,minZ
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read minmal coordinate values in X-, Y- & Z-direction "//&
            "from second line of 3D_structured model file '"//trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
!
!   Reads the third line of the model file, maximanl coordinates (maxX, maxY, maxZ)
    read(lu,*,iostat=ios) maxX,maxY,maxZ
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read maximal coordinate values in X-, Y- & Z-direction "//&
            "from third line of 3D_structured model file '"//trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
!
    nval = (nx)*(ny)*(nz)
!
!   Reads parameters (vp, vs, rho)
    read(lu,fmt="(a)",iostat=ios) s
    sps=' '
    nparam = countWordsString(s,sps,errmsg) ! reads parameters in 4th line of the 3D_structured file and counts given parameters
    if(.level.errmsg == 2 .or. nparam < 1) then
       close(lu)
       write(errstr,*) "ERROR: could not read any ' '-separated words from 4th line of 3D_structured file '"//&
            trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
    allocate(param(nparam))
    read(s,*,iostat=ios) param
    if(ios/=0) then
       close(lu)
       write(errstr,*) "ERROR: could not read ",nparam, " parameter names from 4th line of 3D_structured file '"&
           //trim(filename)//"'"
       call add(errmsg,2,errstr,myname); return
    end if
!
!   Reads the model values from the model file and creates coordinates (from minX, maxX, nX,...)
    allocate(val(nx,ny,nz,nparam))
    allocate(Xcoord(nx))
    allocate(Ycoord(ny))
    allocate(Zcoord(nz))
    do i= 1,nx
       Xcoord(i)=minX+(i*(maxX-minX)/(nx-1))-(maxX-minX)/(nx-1) ! All possible coordinates in X-direction
       do j=1,ny
          Ycoord(j)=minY+(j*(maxY-minY)/(ny-1))-(maxY-minY)/(ny-1) ! All possible coordinates in Y-direction
          do k=1,nz
             Zcoord(k)=minZ+(k*(maxZ-minZ)/(nz-1))-(maxZ-minZ)/(nz-1) ! All possible coordinates in Z-direction
             read(lu,*,iostat=ios) val(i,j,k,:)
             if(ios/=0) then
             close(lu)
             write(errstr,*) "ERROR: could not read real values from 3D_structured model file '"//trim(filename)//"', "//&
                  "please check that the total number of values in the 3D_structued model file is equal to ",nval,""
              call add(errmsg,2,errstr,myname); return
             end if
           end do !k
        end do !j
    end do ! i
!
    close(lu)
  end subroutine read3DstructuredModelFile
!
!------------------------------------------------------------------------
!
  function locateCoordinate(coords,n,c) result(i)
    integer :: n
    real, dimension(n) :: coords
    real :: c
    ! returning
    integer :: i
    ! local
    integer :: jl,jr,jm
    logical :: is_ascending
!
    is_ascending = coords(n) > coords(1)
!
    if(is_ascending) then
!
       if(c < coords(1) .or. c > coords(n)) then
          i = -1
          return
       end if
       jl = 1; jr = n
       do while(jr-jl>1)
          jm = (jr+jl)/2
          if(c >= coords(jm)) then
             jl = jm
          else
             jr = jm
          end if
       end do ! while(jr-jl>1)
!
    else ! is_ascending
!
       if(c > coords(1) .or. c < coords(n)) then
          i = -1
          return
       end if
       jl = 1; jr = n
       do while(jr-jl>1)
          jm = (jr+jl)/2
          if(c <= coords(jm)) then
             jl = jm
          else
             jr = jm
          end if
       end do ! while(jr-jl>1)
!
    end if ! is_ascending
!
    i = jl
  end function locateCoordinate

end program createStartmodelKim
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'Usage:'
!   print *,''
!   print *,'    createStartmodelKim -igtype invgrid_type -igpar invgrid_parfile -igpath invgrid_path '
!   print *,'                        -mpmtrz model_pmtrz -mtype model_type -mfile model_file [-o outbase] [-bin]'
!   print *,''
!   print *,'Mandatory options:'
!   print *,''
!   print *,'-igtype  : invgrid_type is TYPE_INVERSION_GRID as in ASKI iteration step parameter file'
!   print *,''
!   print *,'-igpar : invgrid_parfile is PARFILE_INVERSION_GRID as in ASKI iteration step parameter file'
!   print *,''
!   print *,'-igpath: invgrid_path is treated as current iteration step path, used for inversion grid to write/read own files'
!   print *,''
!   print *,'-mpmtrz: model_pmtrz is the model parametrization of the model which is to be created (must be consistent'
!   print *,'         with content of model_file, see -mfile)'
!   print *,''
!   print *,"-mtype: model_type is a string defining the type of model_file and the interpolation type:"
!   print *,"        '1D_linear' : linear interpolation of values between coordinates given in 1D model file"
!   print *,"        '3D_structured': trilinear interpolation of values between coordinates given in 3D_structured model file"
!   print *,''
!   print *,'        1D model file must have the following format:'
!   print *,'           line 1:   nval  ncol'
!   print *,'              nval = number of model values /interpolation coordinates to come'
!   print *,'              ncol = number of columns to read from file, ncol=1+N, where N is number of parameters'
!   print *,'           line 2:   icoord  param1 ... paramN'
!   print *,'              icoord = index of inversion grid / wavefield point coordinate for which the interpolation'
!   print *,'                       should be applied (either 1,2 or 3)'
!   print *,'              param1 ... paramN   N parameter names associated with the values in the columns below'
!   print *,'           lines 3...nval+2 :   coord val1 ... valN'
!   print *,'              these following nval lines contain the interpolation coordinate and N model values for the '
!   print *,'              respective parameters, as defined by line 2'
!   print *,"              The values 'coord' are assumed to be be strictly monotonical(!) and can be either increasing"
!   print *,'              or decreasing. You may choose this monotonicity at your will dependent on the coordinate for'
!   print *,'              which the interpolation should be done (which, dependent on your inversion grid type may be'
!   print *,'              depth or positive z-value...)'
!   print *,''
!   print *,"        3D_structured model file must have the following format:"
!   print *,'           line 1:   nx  ny  nz'
!   print *,'                     nx,ny and nz are the number of model points in X-, Y- and Z-direction'
!   print *,'           line 2:   minX  minY  minZ'
!   print *,'                     minX, minY and minZ are the smallest coordinates in the model'
!   print *,'           line 3:   maxX  maxY  maxZ'
!   print *,'                     maxX, maxY and maxZ are the highest coordinates in the model'
!   print *,"               nval, the total number of vaules in the model is given by nval = nx*ny*nz"
!   print *,'           line 4:   param1 ... paramN'
!   print *,'               param1 ... paramN   N parameter names associated with the values in the columns below'
!   print *,'           line 5+:  nval model values'
!   print *,'                     The model values have to be sorted like:'
!   print *,'                     DO i=minX,maxX+1'
!   print *,'                       DO j=minY,maxY+1'
!   print *,'                         DO k=minZ,maxZ+1'
!   print *,'                           READ model value'
!   print *,'                         END DO'
!   print *,'                       END DO'
!   print *,'                     END DO'
!   print *,''
!   print *,'-mfile: model_file is the model input file of type defined by model_type'
!   print *,''
!   print *,'Optional options:'
!   print *,''
!   print *,'-o     : outbase is output base name (default is start_model)'
!   print *,''
!   print *,'-bin : if set, the output vtk files will be binary, otherwise they will be ascii'
!   print *,''
!   print *,'-h     : print this help'
!   print '(50(1h-))'
!   return
! end subroutine printhelp
