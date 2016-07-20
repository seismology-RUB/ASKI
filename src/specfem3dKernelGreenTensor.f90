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
!> \brief Module dealing with kernel green tensor for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date June 2013
!
module specfem3dKernelGreenTensor
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateSpecfem3dKernelGreenTensor; end interface
  interface operator (.id.); module procedure getStanameSpecfem3dKernelGreenTensor; end interface
  interface operator (.df.); module procedure getDfSpecfem3dKernelGreenTensor; end interface
  interface operator (.nf.); module procedure getNfSpecfem3dKernelGreenTensor; end interface
  interface operator (.jfcur.); module procedure getJfcurSpecfem3dKernelGreenTensor; end interface
  interface operator (.disp.); module procedure getSpecfem3dKernelGreenTensor; end interface
  interface operator (.strain.); module procedure getStrainsSpecfem3dKernelGreenTensor; end interface
!
  integer, parameter :: length_ID_specfem3d_kernel_green_tensor = 13 !< change this number CONSISTENTLY with length_ASKI_output_ID in SPECFEM codes
!
!> \brief type containing all information on the kernel green tensor files
  type specfem3d_kernel_green_tensor
     private
     character (len=498) :: basename = '' !< file base for all files of this kernel green tensor object
     character (len=length_ID_specfem3d_kernel_green_tensor), dimension(3) :: staname = ''  !< station name
     integer :: nproc = -1 !< number of parallel processors by which this file was created (safety feature, to check with the frequency files)
     integer :: specfem_version = -1 !< for safety reasons, remember this (check with other files)
     double precision :: df = -1.0d0   !< frequency step
     integer :: nfreq = 0  !< number of frequency indices, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null()     !< frequency indices
     integer :: jfcur = -1  !< indicates current frequency index for which wavefield values are contained in arrays g,gstr
     integer :: ntot = 0   !< total number of wavefield points (size(g,1), size(gstr,1))
     complex, dimension(:,:,:), pointer :: g => null()   !< dimension 1: ntot, dimension 2: wavefield component 1,..,3, dimension 3: force component 1,..,3
     complex, dimension(:,:,:), pointer :: gstr => null() !< dimension 1: ntot, dimension 2: strain component 1,..,6, dimension 3: force component 1,..,3
     integer :: lu = 0 !< file unit used for this object (hold this for the whole time between initialRead and deallocate)
  end type specfem3d_kernel_green_tensor
!
contains
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed from "main" file of this kernel green tensor object
!
  subroutine initialReadSpecfem3dKernelGreenTensor(this,fuh,basename,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    integer :: ier,type_invgrid_X,len_id
    integer :: j,j_file,specfem_version,nproc,type_invgrid,ntot,nf
    character(len=2) :: extension
    double precision :: df
    integer, dimension(:), allocatable :: jf
    character (len=37) :: myname = 'initialReadSpecfem3dKernelGreenTensor'
!
    call addTrace(errmsg,myname)
!
    this%lu = get(fuh)
    this%basename = basename
!
    ! FIRST FILE, X-direction green function
!
    open(unit=this%lu,file=trim(this%basename)//'_X.main',status='old',action='read',access='stream',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(this%basename)//"_X.main' cannot be opened to read",myname)
       goto 1
    else
       call add(errmsg,0,"Successfully opened file '"//trim(this%basename)//"_X.main' to read",myname)
    endif

    read(this%lu) this%specfem_version
    select case(this%specfem_version)
    case(1,2)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "specfem version = ",this%specfem_version,&
            " is not supported: 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select

    read(this%lu) len_id
    if(len_id/=length_ID_specfem3d_kernel_green_tensor) then
       write(errstr,*) "length of ASKI ID is ",len_id,&
            ", only length ",length_ID_specfem3d_kernel_green_tensor," supported here. Please update this code accordingly"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if

    read(this%lu) this%staname(1),this%nproc,type_invgrid_X,this%ntot,this%df,this%nfreq

    ! check if specfem_version and type_invgrid_X are a valid match
    if(this%specfem_version == 1) then
       select case(type_invgrid_X)
       case(1,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_invgrid_X," is not supported by SPECFEM3D_GLOBE"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(this%specfem_version == 2) then
       select case(type_invgrid_X)
       case(2,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_invgrid_X," is not supported by SPECFEM3D_Cartesian"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if

    if(this%nfreq<1 .or. this%ntot<1) then
       write(errstr,*) "ntot,df,nf =",this%ntot,this%df,this%nfreq,"; ntot and nf should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    allocate(this%jf(this%nfreq))
    read(this%lu) this%jf

    ! everything went alright reading in information from first main file (X-direction green function)
    close(this%lu)
!
    ! NOW CHECK, IF THE OTHER MAIN FILES (Y-,Z-direction green functions) CONTAIN THE SAME INFORMATION!
    do j_file = 1,2
       select case(j_file)
       case(1); extension='_Y'
       case(2); extension='_Z'
       end select

       open(unit=this%lu,file=trim(this%basename)//extension//'.main',status='old',action='read',access='stream',&
            form='unformatted',iostat=ier)
       if (ier /= 0) then
          call add(errmsg,2,"File '"//trim(this%basename)//extension//".main' cannot be opened to read",myname)
          goto 1
       else
          call add(errmsg,0,"Successfully opened file '"//trim(this%basename)//extension//".main' to read",myname)
       endif

       read(this%lu) specfem_version
       if(specfem_version/=this%specfem_version) then
          write(errstr,*) "specfem_version ",specfem_version," contained in file '"//extension//&
               ".main' does not match specfem_version ",this%specfem_version," of file '_X.main'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       read(this%lu) len_id
       if(len_id/=length_ID_specfem3d_kernel_green_tensor) then
          write(errstr,*) "length of ASKI ID is ",len_id,&
               ", only length ",length_ID_specfem3d_kernel_green_tensor,&
               " supported here. Please update this code accordingly"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       read(this%lu) this%staname(j_file+1),nproc,type_invgrid,ntot,df,nf

       if(nproc/=this%nproc) then
          write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc," of _X.main file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(type_invgrid/=type_invgrid_X) then
          write(errstr,*) "type of inversion grid ",type_invgrid," contained in this file does not match type ",&
               type_invgrid_X," in _X.main file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(ntot/=this%ntot) then
          write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot," of _X.main file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(df/=this%df) then
          write(errstr,*) "df ",df," contained in this file does not match df ",this%df," of _X.main file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(nf/=this%nfreq) then
          write(errstr,*) "nf ",nf," contained in this file does not match nf ",this%nfreq," of _X.main file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       allocate(jf(nf))
       read(this%lu) jf
       do j=1,nf
          if(.not.any(jf==this%jf(j))) then
             write(errstr,*) "frequency index ",this%jf(j),", contained in file _X.main, "//&
                  "is NOT contained in this file!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if 
       end do ! j
       deallocate(jf)

       ! everything went alright reading in information from second or third main file
       close(this%lu)

    end do ! j_file

    ! finally allocate for spectra
    allocate(this%g(this%ntot,3,3),this%gstr(this%ntot,6,3))

    ! if program comes here, there was no error reading the files and jf is deallocated, so return
    return

    ! if there is an error reading the files, deallocate everything before return
1   close(this%lu)
    if(allocated(jf)) deallocate(jf)
    call deallocateSpecfem3dKernelGreenTensor(this,fuh)
  end subroutine initialReadSpecfem3dKernelGreenTensor
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocateSpecfem3dKernelGreenTensor(this,fuh)
    type (specfem3d_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    this%basename = ''
    this%staname = ''
    this%nproc = -1
    this%specfem_version = -1
    this%df = -1.0d0
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
    this%ntot = 0
    if(associated(this%g)) deallocate(this%g)
    if(associated(this%gstr)) deallocate(this%gstr)
    call add(fuh,this%lu); this%lu = 0
  end subroutine deallocateSpecfem3dKernelGreenTensor
!--------------------------------------------------------------------------------------
!> \brief read specific files of kernel green tensor spectrum for the requested frequency (if existend)
!! \param jf frequency index for which kernel green tensor should be read in (defines filename)
!
  subroutine readFrequencySpecfem3dKernelGreenTensor(this,jf,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    integer :: jf
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=39) :: myname = 'readFrequencySpecfem3dKernelGreenTensor'
    character(len=509) :: filename
    integer :: j_file,ier,len_id,nproc,specfem_version,ntot,jf_file
    double precision :: df
    character(len=length_ID_specfem3d_kernel_green_tensor) :: id
    character(len=2) :: extension
    !
    call addTrace(errmsg,myname)
!
    ! check validity of frequency index
!
    if(.not.any(this%jf == jf)) then
       write(errstr,*) "incoming frequency index ",jf,", not contained in this kernelGreenTensor object. "
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! open specific files for this frequency
!
    do j_file = 1,3
       select case(j_file)
       case(1); extension='_X'
       case(2); extension='_Y'
       case(3); extension='_Z'
       end select
       
       ! open specific file for this frequency
       write(filename,"(a,i6.6)") trim(this%basename)//extension//".jf",jf
       open(unit=this%lu,file=filename,status='old',action='read',access='stream',&
            form='unformatted',iostat=ier)
       if (ier /= 0) then
          call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
          close(this%lu)
          return
       else
          call add(errmsg,0,"Successfully opened file '"//trim(filename)//"' to read",myname)
       endif

       read(this%lu) specfem_version
       if(specfem_version/=this%specfem_version) then
          write(errstr,*) "specfem_version ",specfem_version," contained in this file does not match specfem_version ",&
               this%specfem_version," of main files"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if

       read(this%lu) len_id
       if(len_id/=length_ID_specfem3d_kernel_green_tensor) then
          write(errstr,*) "length of ASKI ID is ",len_id,&
               ", only length ",length_ID_specfem3d_kernel_green_tensor," supported here. Please update this code accordingly"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if

       read(this%lu) id,nproc,ntot,df,jf_file
       if(id/=this%staname(j_file)) then
          call add(errmsg,2,"ID contained in this file '"//trim(id)//"' does not match ID '"//&
               trim(this%staname(j_file))//"' of main files",myname)
          close(this%lu)
          return
       end if
       if(nproc/=this%nproc) then
          write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc," of main files"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if
       if(ntot/=this%ntot) then
          write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot," of main files"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if
       if(df/=this%df) then
          write(errstr,*) "df ",df," contained in this file does not match df ",this%df," of main files"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if
       if(jf_file/=jf) then
          write(errstr,*) "requested frequency index ",jf," (also contained in filename) does not match the frequency index ",&
               jf_file,"actually contained in the file"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if

       ! make this%jf invalid before changing values in arrays this%g,this%gstr
       ! in case of returning on error, it can be no longer assured what is contained in the object
       ! when returning on success (end of subroutine), set this%jfcur properly
       this%jfcur = -1
       ! read spectrum from file
       read(this%lu) this%g(:,:,j_file),this%gstr(:,:,j_file)

       close(this%lu)

    end do ! j_file

    this%jfcur = jf
  end subroutine readFrequencySpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfSpecfem3dKernelGreenTensor(this) result(res)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    real :: res
    res = real(this%df)
  end function getDfSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequency indices nfreq
!
  function getNfSpecfem3dKernelGreenTensor(this) result(res)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nfreq
  end function getNfSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurSpecfem3dKernelGreenTensor(this) result(res)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%jfcur
  end function getJfcurSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get strains
!
  function getStrainsSpecfem3dKernelGreenTensor(this) result(gstr)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: gstr
    gstr => this%gstr
  end function getStrainsSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get Green tensor
!
  function getSpecfem3dKernelGreenTensor(this) result(g)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: g
    g => this%g
  end function getSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get station name
!
  function getStanameSpecfem3dKernelGreenTensor(this) result(res)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    character(len=length_ID_specfem3d_kernel_green_tensor) :: res
    res = this%staname(1)
  end function getStanameSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencySpecfem3dKernelGreenTensor(this,jf) result(next)
    type (specfem3d_kernel_green_tensor) :: this
    integer :: jf
    logical :: next
    integer :: call_count = 0
    save call_count
    !
    if (call_count == this%nfreq) then
       call_count = 0
       next = .false.; return
    endif
    call_count = call_count+1
    jf = this%jf(call_count)
    next = .true.
  end function nextFrequencySpecfem3dKernelGreenTensor
!
end module specfem3dKernelGreenTensor
