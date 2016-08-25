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
program kgt2vtk
  use inversionBasics
  use iterationStepBasics
  use kernelGreenTensor
  use seismicNetwork
  use wpVtkFile
  use argumentParser
  use componentTransformation, only: all_valid_components
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  type (error_message) :: errmsg
  type (file_unit_handler) :: fuh
  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  character(len=max_length_string) :: main_parfile

  character(len=max_length_string) :: staname
  integer, dimension(:), pointer :: ifreq,ifreq_iterbasics
  character(len=max_length_string), dimension(:), pointer :: scomp,ucomp
  logical :: use_all_ucomp,use_selected_ucomp,use_all_ifreq,use_selected_ifreq
  integer :: nscomp,nucomp,iscomp,iucomp,nfreq,jfreq
  real :: df

  type (kernel_green_tensor) :: kgt
  complex, dimension(:,:,:), pointer :: kgt_ustr,kgt_u
  character(len=max_length_string) :: kgt_file
  integer :: nwp,un,en

  type (wp_vtk_file), dimension(:,:), allocatable :: wp_vtk
  complex, dimension(:), allocatable :: data
  character(len=max_length_string) :: vtk_file_base,vtk_file_title,vtk_file_data_name

  character (len=7) :: myname = 'kgt2vtk'

  nullify(ifreq,ifreq_iterbasics,scomp,ucomp,kgt_ustr,kgt_u)

!------------------------------------------------------------------------
!  definition and basic processing of command line
!
  call init(ap,myname,'Extract kernel green tensor spectra to vtk files for certain wavefield and strain components and '//&
  'frequencies')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-staname',.true.,"defines the station name of the kernel green tensor object. This option must be set.",&
       'sval','')
  call addOption(ap,'-scomp',.true.,"explicit vector of station components which can be '"//trim(all_valid_components)//&
       "'. This option must be set.",'svec','')
  call addOption(ap,'-ifreq',.true.,"explicit vector of frequency indices at which the wavefield output should be extracted. "//&
       "Exactly one of options -ifreq , -all_ifreq must be set",'ivec','')
  call addOption(ap,'-all_ifreq',.false.,"if set, all frequency indices are used. Exactly one of options -ifreq , "//&
       "-all_ifreq must be set")
  call addOption(ap,'-ucomp',.true.,"explicit vector of wavefield components which can be 'ux', 'uy', 'uz' (denoting underived "//&
       "x,y,z components of wavefield) and 'exx', 'eyy', 'ezz', 'eyz', 'exz', 'exy' (denoting the strain components). "//&
       "Exactly one of options -ucomp , -all_ucomp must be set'",'svec','')
  call addOption(ap,'-all_ucomp',.false.,"if set, all wavefield components are (3 underived and 6 strain components). "//&
       "Exactly one of options -ucomp , -all_ucomp must be set")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  use_all_ucomp = ap.optset.'-all_ucomp'
  use_selected_ucomp = ap.optset.'-ucomp'
  if(use_all_ucomp .eqv. use_selected_ucomp) then
     write(*,*) "ERROR: exactly ONE of the options -ucomp and -all_ucomp must be set!"
     call usage(ap)
     goto 1
  end if
!
  use_all_ifreq = ap.optset.'-all_ifreq'
  use_selected_ifreq = ap.optset.'-ifreq'
  if(use_all_ifreq .eqv. use_selected_ifreq) then
     write(*,*) "ERROR: exactly ONE of the options -ifreq and -all_ifreq must be set!"
     call usage(ap)
     goto 1
  end if
!
  if(.not.(ap.optset.'-staname')) then
     write(*,*) "ERROR: option -staname must be set!"
     call usage(ap)
     goto 1
  end if
  staname = ap.sval.'-staname'
!
  if(.not.(ap.optset.'-scomp')) then
     write(*,*) "ERROR: option -scomp must be set!"
     call usage(ap)
     goto 1
  end if
  scomp => ap.svec.'-scomp'
  if(.not.(associated(scomp))) then
     write(*,*) "ERROR: no components could be read from the argument of -scomp"
     call usage(ap)
     goto 1
  end if
  nscomp = size(scomp)
  do iscomp = 1,nscomp
     if(.not.validComponent(scomp(iscomp))) then
        write(*,*) "ERROR: ",iscomp,"'th station component '",trim(scomp(iscomp)),&
             "' of the -scomp argument string is not one of '",trim(all_valid_components),"'"
        call usage(ap)
        goto 1
     end if
  end do ! iscomp
!------------------------------------------------------------------------
!  setup basics
!
  call new(fuh,20)
!
  call new(errmsg,myname)
  call init(invbasics,trim(main_parfile),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!  setup iteration step basics
!
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  ifreq_iterbasics => .ifreq.iterbasics
  df = .df.invbasics
!------------------------------------------------------------------------
!  detailed processing of command line arguments
!
  if(use_selected_ucomp) then
     ucomp => ap.svec.'-ucomp'
     if(.not.(associated(ucomp))) then
        write(*,*) "ERROR: no components could be read from the argument of -ucomp"
        call usage(ap)
        goto 1
     end if
     nucomp = size(ucomp)
     do iucomp = 1,nucomp
        select case (ucomp(iucomp))
        case('ux','uy','uz','exx','eyy','ezz','eyz','exz','exy')
           ! OK, do nothing
        case default
           write(*,*) "ERROR: ",iucomp,"'th wavefield component '",trim(ucomp(iucomp)),&
                "' of the -ucomp argument string is not one of 'ux','uy','uz','exx','eyy','ezz','eyz','exz','exy'"
           call usage(ap)
           goto 1
        end select
     end do ! iucomp
  else ! use_selected_ucomp
     nucomp = 9
     allocate(ucomp(nucomp))
     ucomp = (/'ux ','uy ','uz ','exx','eyy','ezz','eyz','exz','exy'/)
  end if ! use_selected_ucomp
!
  if(use_selected_ifreq) then
     ifreq => ap.ivec.'-ifreq'
     if(.not.(associated(ifreq))) then
        write(*,*) "ERROR: no frequency indices could be read from the argument of -ifreq"
        call usage(ap)
        goto 1
     end if
     nfreq = size(ifreq)
     do jfreq = 1,nfreq
        if(.not.any(ifreq(jfreq)==ifreq_iterbasics)) then
           write(*,*) "ERROR: ",jfreq,"'th frequency index ",ifreq(jfreq)," of the -ifreq argument string is not "//&
                "contained in the frequency indices of the current iteration step: ",ifreq_iterbasics
           call usage(ap)
           goto 1
        end if
     end do ! iucomp
  else ! use_selected_ifreq
     nfreq = size(ifreq_iterbasics)
     allocate(ifreq(nfreq))
     ifreq = ifreq_iterbasics
  end if ! use_selected_ifreq
!
  errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname)
  if(.level. errmsg/=0) then
     write(*,*) "ERROR: station name '"//trim(staname)//"' (input string of option -staname) is not contained in station list"
     goto 1
  end if
  call dealloc(errmsg)
!
  call document(ap)
  write(*,*) ""
!------------------------------------------------------------------------
!  prepare for the loop below
!
  kgt_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_GREEN_TENSORS')//&
       'kernel_gt_'//trim(staname)
  write(*,*) "OPEN KERNEL GREEN TENSOR FILE '",trim(kgt_file),"' TO READ"
  call new(errmsg,myname)
  call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh,kgt_file,&
       scomp,.comptrans.invbasics,staname,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  allocate(wp_vtk(nucomp,nscomp))
  nwp = .ntot.(.wp.iterbasics)
  allocate(data(nwp))
  
!------------------------------------------------------------------------
!  now loop on all frequencies and components, read kernel displ and write
!
  do jfreq = 1,nfreq
     write(*,*) "PROCESSING FREQUENCY INDEX ",ifreq(jfreq)

     call new(errmsg,myname) ! use one error message for the next few calls

     call readFrequencyKernelGreenTensor(kgt,ifreq(jfreq),errmsg)
     if (.level.errmsg == 2) then; call print(errmsg); goto 1; endif

     call getKernelGreenTensor(kgt,kgt_u,errmsg)
     if (.level.errmsg == 2) then; call print(errmsg); goto 1; endif
     if(.not.associated(kgt_u)) then
        write(*,*) "ERROR: no kernel green tensor values were returned, this error should not have occurred!"
        if (.level.errmsg /= 0) call print(errmsg)
        goto 1
     end if
     ! kgt_u: dimension 1 = wavefield points, dimension 2 = wavefield component, dimension 3 = force component of green function
     if(size(kgt_u,1) /= nwp) then
        write(*,*) "ERROR: kernel green tensor object has ",size(kgt_u,1)," values, but there are ",nwp," wavefield points"
        if (.level.errmsg /= 0) call print(errmsg)
        goto 1
     end if

     call getStrainsKernelGreenTensor(kgt,kgt_ustr,errmsg)
     if (.level.errmsg == 2) then; call print(errmsg); goto 1; endif
     ! kgt_ustr: dimension 1 = wavefield points, dimension 2 = strain component, dimension 3 = force component of green function
     if(.not.associated(kgt_ustr)) then
        write(*,*) "ERROR: no kernel green tensor strain values were returned, this error should not have occurred!"
        if (.level.errmsg /= 0) call print(errmsg)
        goto 1
     end if
     if(size(kgt_ustr,1) /= nwp) then
        write(*,*) "ERROR: kernel green tensor strains have ",size(kgt_u,1)," values, but there are ",nwp," wavefield points"
        if (.level.errmsg /= 0) call print(errmsg)
        goto 1
     end if

     if (.level.errmsg /= 0) call print(errmsg)
     call dealloc(errmsg)

     ! loop on all station components
     do iscomp=1,nscomp

        ! loop on all wavefield components
        do iucomp=1,nucomp
           un = -1; en = -1
           select case(ucomp(iucomp))
           case ('ux'); un = 1
           case ('uy'); un = 2
           case ('uz'); un = 3
           case ('exx'); en = 1
           case ('eyy'); en = 2
           case ('ezz'); en = 3
           case ('eyz'); en = 4
           case ('exz'); en = 5
           case ('exy'); en = 6
           end select

           ! in case that this is an underived wavefield component, get kernel displacement at correct component
           if(un > 0) data = kgt_u(:,un,iscomp)
           ! in case that this is a strain component, get strain at correct component
           if(en > 0) data = kgt_ustr(:,en,iscomp)

           ! finally write vtk file

           if(jfreq==1) then
              write(vtk_file_base,"(a,'_',a,'_',a)") trim(kgt_file),trim(scomp(iscomp)),trim(ucomp(iucomp))
              ! initiate vtk file
              write(vtk_file_title,*) trim(ucomp(iucomp)),"-component of spectral kernel green tensor, receiver component ",&
                   trim(scomp(iscomp))," at frequency ",ifreq(jfreq)*df,' Hz on wavefield points'
              call new(errmsg,myname)
              call init(wp_vtk(iucomp,iscomp),.wp.iterbasics,.invgrid.iterbasics,trim(vtk_file_base),&
                   trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtk_file_title))
              if (.level.errmsg /= 0) call print(errmsg)
              if (.level.errmsg == 2) goto 1
              call dealloc(errmsg)
              print *,"    creating vtk files with basename '"//trim(vtk_file_base)//"' (plus extension for each frequency index)"
           end if ! jfreq==1

           ! write kgt values to vtk file
           write(vtk_file_data_name,*) trim(ucomp(iucomp)),'-kgt(',trim(scomp(iscomp)),')'
           call new(errmsg,myname)
           call writeData(wp_vtk(iucomp,iscomp),get(fuh),data,errmsg,data_name=trim(vtk_file_data_name),file_index=ifreq(jfreq))
           call undo(fuh)
           if (.level.errmsg /= 0) call print(errmsg)
           if (.level.errmsg == 2) goto 1
           call dealloc(errmsg)

        end do ! iucomp
     end do ! iscomp

     if(associated(kgt_ustr)) deallocate(kgt_ustr)
     if(associated(kgt_u)) deallocate(kgt_u)

  end do ! jfreq

!------------------------------------------------------------------------
!  clean up before terminating the program
!
1   call dealloc(iterbasics); call dealloc(invbasics)
    call dealloc(fuh)
    call dealloc(errmsg)
    call dealloc(ap)
    if(associated(ifreq)) deallocate(ifreq)
    if(associated(ucomp)) deallocate(ucomp)
    if(associated(kgt_ustr)) deallocate(kgt_ustr)
    if(associated(kgt_u)) deallocate(kgt_u)
    if(allocated(wp_vtk)) then
       do iscomp = 1,nscomp
          do iucomp = 1,nucomp
             call dealloc(wp_vtk(iucomp,iscomp))
          end do ! iucomp
       end do ! iscomp
       deallocate(wp_vtk)
    end if
    if(allocated(data)) deallocate(data)

end program kgt2vtk
