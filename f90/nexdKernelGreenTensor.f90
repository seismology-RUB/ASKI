!----------------------------------------------------------------------------
!   Copyright 2016 Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module dealing with kernel green tensor for method NEXD
!! \details The code of this module is based on module specfem3dKernelGreenTensor
!!
!! \author Christian Ullisch
!! \date Nov 2015
!
module nexdKernelGreenTensor
  use componentTransformation
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateNexdKernelGreenTensor; end interface
  interface operator (.id.); module procedure getStanameNexdKernelGreenTensor; end interface
  interface operator (.df.); module procedure getDfNexdKernelGreenTensor; end interface
  interface operator (.nf.); module procedure getNfNexdKernelGreenTensor; end interface
  interface operator (.jfcur.); module procedure getJfcurNexdKernelGreenTensor; end interface
!
  integer, parameter :: length_ID_nexd_kernel_green_tensor = 13 !< change this number CONSISTENTLY with module kernelGreenTensor
!
!> \brief type containing all information on the kernel green tensor files
  type nexd_kernel_green_tensor
     private
     ! basic information
     character (len=498) :: basename = '' !< file base for all files of this kernel green tensor object
     character (len=length_ID_nexd_kernel_green_tensor) :: staname = '' !< station name (string before '_' in the ID contained in the file
     integer :: nproc = -1 !< number of parallel processors by which this file was created (safety feature, to check with the frequency files)
     integer :: lu = -1 !< file unit used for this object (hold this for the whole time between initialRead and deallocate)
!
     ! frequency content
     real :: df = -1.   !< frequency step
     integer :: nfreq = 0  !< number of frequency indices, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null()     !< frequency indices
     integer :: jfcur = -1  !< indicates current frequency index for which wavefield values are contained in arrays g,gstr
!
     ! ACTUAL COMPONENTS WHICH ARE HANDLED IN THIS OBJECT (can be any subset of all available components)
     integer :: ncomp = 0 !< number of actual components which are managed by this object
     character (len=character_length_component), dimension(:), pointer :: comp => null() !< list of actually contained components of this object (subset of available components), this represents the content of arrays g,gstr
!
     ! green tensor values and strains for requested force components (3rd dimensions according to actually contained components comp)
     integer :: ntot = 0   !< total number of wavefield points (size(g,1), size(gstr,1))
     complex, dimension(:,:,:), pointer :: g => null()   !< dimension 1: ntot, dimension 2: wavefield component 1,..,3, dimension 3: requested force component 1,..,ncomp
     complex, dimension(:,:,:), pointer :: gstr => null() !< dimension 1: ntot, dimension 2: strain component 1,..,6, dimension 3: requested force component 1,..,ncomp
  end type nexd_kernel_green_tensor
!
contains
!----------------------------------------------------------------------------------------------------
!> \brief check which components are available at all
!
  subroutine getAvailableComponentsNexdKernelGreenTensor(basename,lu,comp_available,errmsg)
    character(len=*) :: basename
    integer :: lu
    character (len=character_length_component), dimension(:), pointer :: comp_available
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=43) :: myname = 'getAvailableComponentsNexdKernelGreenTensor'
    integer :: ncomp_available,icomp,ier
!
    call addTrace(errmsg,myname)
    nullify(comp_available)
!
    open(unit=lu,file=trim(basename)//'.comp',status='old',action='read',&
         form='formatted',iostat=ier)
    if(ier /= 0) then
       call add(errmsg,2,"File '"//trim(basename)//".comp' cannot be opened to read",myname)
       goto 1
    else
       call add(errmsg,0,"Successfully opened file '"//trim(basename)//".comp' to read",myname)
    endif
!
    read(lu,*,iostat=ier) ncomp_available
    if(ier /= 0) then
       call add(errmsg,2,"could not read number of components from first line of *.comp file",myname)
       goto 1
    end if
    if(ncomp_available <= 0) then
       write(errstr,*) "number of components on first line of *.comp file = ",ncomp_available,&
            ". must be strictly positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    else
       write(errstr,*) "number on first line of *.comp file equals ",ncomp_available,&
            "; expect that many lines of components following"
       call add(errmsg,0,errstr,myname)
    end if
!
    allocate(comp_available(ncomp_available))
    do icomp = 1,ncomp_available
       read(lu,*,iostat=ier) comp_available(icomp)
       if(ier /= 0) then
          write(errstr,*) "could not read ",icomp,"'th component from ",icomp+1,"'th line of *.comp file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if ! ier/=0
       if(.not.validComponent(comp_available(icomp))) then
          write(errstr,*) icomp,"'th component on ",icomp+1,"'th line of *.comp file ('",trim(comp_available(icomp)),&
               "') is invalid! valid components are: ",all_valid_components
          call add(errmsg,2,errstr,myname)
          goto 1
       end if ! validComponent
       if(icomp > 1) then
          if(any(comp_available(1:icomp-1) == comp_available(icomp))) then
             write(errstr,*) icomp,"'th component on ",icomp+1,"'th line of *.comp file ('",trim(comp_available(icomp)),&
                  "') did already occur on the previous lines. There seems to be an inconsistency of the *.comp file"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if ! any(..)
       end if ! icomp>0
    end do ! icomp
!
    ! if code comes here, everything went alright, so simply return
    return
!
1   close(lu)
    if(associated(comp_available)) deallocate(comp_available)
  end subroutine getAvailableComponentsNexdKernelGreenTensor
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed from "main" file of this kernel green tensor object
!
  subroutine initialReadNexdKernelGreenTensor(this,fuh,basename,staname,comp,errmsg)
    type (nexd_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename,staname
    character(len=*), dimension(:) :: comp
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=508) :: filename
    integer :: ier,icomp
    integer :: j,nproc,ntot,nf
    real :: df
    integer, dimension(:), allocatable :: jf
    character (len=32) :: myname = 'initialReadNexdKernelGreenTensor'
!
    ! IN MODULE kernelGreenTensor IT WAS CHECKED THAT comp IS NOT EMPTY AND ALL ENTRIES ARE 
    ! VALID COMPONENTS AND THAT ALL ENTRIES IN VECTOR comp ARE AVAILABLE COMPONENTS FOR THIS OBJECT
!
    call addTrace(errmsg,myname)
    write(errstr,*) "initiating Nexd Green tensor object for station '",trim(staname),"' for components '"!,comp//",","'"
    do j = 1,size(comp)
       errstr = trim(errstr)//trim(comp(j))//","
    end do ! j
    errstr = trim(errstr)//"'"
    call add(errmsg,0,errstr,myname)
!
    this%lu = get(fuh)
    this%basename = basename
!
    ! in module kernelGreenTensor, there was already a call to 
    ! getAvailableComponentsNexdKernelGreenTensor which reads in file this%basename//'.comp' which 
    ! contains the components of Green functions which were written from NEXD and are 
    ! availabe to be read in. In module kernelGreenTensor it was already decided for which components
    ! this object should be initiated (entries in incoming vector comp)
    ! Here, there is no need anymore to decide if those can be transformed to the requested components 
    ! (which are not known to module nexdKernelGreenTensor)
    this%ncomp = size(comp)
    allocate(this%comp(this%ncomp))
    this%comp = comp
!
    ! FIRST FILE, green function component this%comp(1)
!
    filename = trim(this%basename)//'_'//trim(this%comp(1))//'.main'
    open(unit=this%lu,file=filename,status='old',action='read',access='stream',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
       goto 1
    else
       call add(errmsg,0,"Successfully opened file '"//trim(filename)//"' to read",myname)
    endif

    read(this%lu) this%nproc,this%ntot,this%df,this%nfreq

    this%staname = staname

    if(this%nfreq<1 .or. this%ntot<1 .or. this%df<0) then
       write(errstr,*) "ntot,df,nf =",this%ntot,this%df,this%nfreq,"; ntot, df and nf should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if

    allocate(this%jf(this%nfreq))
    read(this%lu) this%jf

    if(any(this%jf<1)) then
       call add(errmsg,2,"there are frequency indices which are <= 0",myname)
       goto 1
    end if
    do j = 2,this%nfreq
       if(any(this%jf(1:j-1)==this%jf(j))) then
          write(errstr,*) j,"'th frequency index read from file = ",this%jf(j),&
               " occurrs twice in vector of frequency indices. Duplicate frequency ",&
               " indices inicate an inconsistency of the files."
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
    end do ! j

    ! everything went alright reading in information from first component main file
    close(this%lu)
!
    ! NOW CHECK, IF THE MAIN FILES OF THE OTHER COMPONENTS CONTAIN THE SAME INFORMATION!
    do icomp = 2,this%ncomp ! if this%ncomp == 1, this loop does nothing (so works for all cases)

       filename = trim(this%basename)//'_'//trim(this%comp(icomp))//'.main'
       open(unit=this%lu,file=filename,status='old',action='read',access='stream',&
            form='unformatted',iostat=ier)
       if (ier /= 0) then
          call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
          goto 1
       else
          call add(errmsg,0,"Successfully opened file '"//trim(filename)//"' of component '"//trim(this%comp(icomp))//&
               "'to read",myname)
       endif

       read(this%lu) nproc,ntot,df,nf

       if(nproc/=this%nproc) then
          write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc,&
               " of the .main file of component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(ntot/=this%ntot) then
          write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot,&
               " of .main file of component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(df/=this%df) then
          write(errstr,*) "df ",df," contained in this file does not match df ",this%df,&
               " of .main file of component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(nf/=this%nfreq) then
          write(errstr,*) "nf ",nf," contained in this file does not match nf ",this%nfreq,&
               " of .main file of component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       allocate(jf(nf))
       read(this%lu) jf
       do j=1,nf
          if(.not.any(jf==this%jf(j))) then
             write(errstr,*) "frequency index ",this%jf(j),", contained in file of component '",trim(this%comp(1)),&
                  "' is NOT contained in this file!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if 
       end do ! j
       deallocate(jf)

       ! everything went alright reading in information from the other main files
       close(this%lu)

    end do ! icomp

    ! finally allocate for spectra
    allocate(this%g(this%ntot,3,this%ncomp),this%gstr(this%ntot,6,this%ncomp))

    ! if program comes here, there was no error reading the files and jf is deallocated, so return
    return

    ! if there is an error reading the files, deallocate everything before return
1   close(this%lu)
    if(allocated(jf)) deallocate(jf)
    call deallocateNexdKernelGreenTensor(this,fuh)
  end subroutine initialReadNexdKernelGreenTensor
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocateNexdKernelGreenTensor(this,fuh)
    type (nexd_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    this%basename = ''
    this%staname = ''
    this%nproc = -1
    this%df = -1.0d0
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
    this%ncomp = 0
    if(associated(this%comp)) deallocate(this%comp)
    this%ntot = 0
    if(associated(this%g)) deallocate(this%g)
    if(associated(this%gstr)) deallocate(this%gstr)
    call add(fuh,this%lu); this%lu = -1
  end subroutine deallocateNexdKernelGreenTensor
!--------------------------------------------------------------------------------------
!> \brief read specific files of kernel green tensor spectrum for the requested frequency (if existing)
!! \param jf frequency index for which kernel green tensor should be read in (defines filename)
!
  subroutine readFrequencyNexdKernelGreenTensor(this,jf,errmsg)
    type (nexd_kernel_green_tensor) :: this
    integer :: jf
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=34) :: myname = 'readFrequencyNexdKernelGreenTensor'
    character(len=512) :: filename
    integer :: icomp,ier,nproc,ntot,jf_file
    real :: df
!
    call addTrace(errmsg,myname)
!
    ! check if this was initiated
    if(this%lu == -1) then
       call add(errmsg,2,"Green tensor object not yet initiated, call initialReadNexdKernelGreenTensor first"&
            ,myname)
       return
    end if
!
    ! check validity of frequency index
!
    if(.not.any(this%jf == jf)) then
       write(errstr,*) "incoming frequency index ",jf,", not contained in this kernelGreenTensor object."
       call add(errmsg,2,errstr,myname)
       return
    else
       write(errstr,*) "will read kernel green tensor for frequency index ",jf
       call add(errmsg,0,errstr,myname)
    end if
!
    ! make this%jfcur invalid before changing values in arrays this%g,this%gstr.
    ! in case of returning on error, it can be no longer assured what is contained in the object.
    ! when returning on success (end of subroutine), set this%jfcur properly
    this%jfcur = -1
!
    ! open specific files for this frequency
!
    do icomp = 1,this%ncomp
       
       ! open specific file for this frequency
       write(filename,"(a,i6.6)") trim(this%basename)//"_"//trim(this%comp(icomp))//".jf",jf
       open(unit=this%lu,file=filename,status='old',action='read',access='stream',&
            form='unformatted',iostat=ier)
       if (ier /= 0) then
          call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
          close(this%lu)
          return
       else
          call add(errmsg,0,"Successfully opened file '"//trim(filename)//"' to read",myname)
       endif

       !read(this%lu) id,nproc,ntot,df,jf_file
       read(this%lu) nproc,ntot,df,jf_file

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
               jf_file," which is actually contained in the file"
          call add(errmsg,2,errstr,myname)
          close(this%lu)
          return
       end if

       ! read spectrum from file
       read(this%lu) this%g(:,:,icomp),this%gstr(:,:,icomp)

       close(this%lu)

    end do ! icomp

    this%jfcur = jf
  end subroutine readFrequencyNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfNexdKernelGreenTensor(this) result(res)
    type (nexd_kernel_green_tensor), intent(in) :: this
    real :: res
    res = real(this%df)
  end function getDfNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequency indices nfreq
!
  function getNfNexdKernelGreenTensor(this) result(res)
    type (nexd_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nfreq
  end function getNfNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurNexdKernelGreenTensor(this) result(res)
    type (nexd_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%jfcur
  end function getJfcurNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get strains at specific components. Raise error if incoming components do not match those for which this tensor object was initiated
!
  subroutine getStrainsNexdKernelGreenTensor(this,gstr,errmsg)
    type (nexd_kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: gstr
    type (error_message) :: errmsg
    character (len=31) :: myname = 'getStrainsNexdKernelGreenTensor'
!
    call addTrace(errmsg,myname)
    nullify(gstr)
!
    if(this%jfcur == -1) then
       call add(errmsg,2,"this object does not contain any sensible values yet. "//&
            "call readFrequencyNexdKernelGreenTensor first",myname)
       return
    end if
!
    gstr => this%gstr
  end subroutine getStrainsNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsNexdKernelGreenTensor(this,uf_gstr,errmsg)
    type (nexd_kernel_green_tensor) :: this
    real :: uf_gstr
    type (error_message) :: errmsg
    if(this%staname /= '') then
       ! NEXD assumes SI units, so unit factor of strains is always 1.0
       uf_gstr = 1.0
    else
       ! object was not yet initiated, raise error
       uf_gstr = 0.
       call add(errmsg,2,"nexd_kernel_green_tensor object not yet initiated",&
            "getUnitFactorStrainsNexdKernelGreenTensor")
    end if
  end subroutine getUnitFactorStrainsNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get Green tensor at specific components. Raise error if incoming components do not match those for which this tensor object was initiated
!
  subroutine getNexdKernelGreenTensor(this,g,errmsg)
    type (nexd_kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: g
    type (error_message) :: errmsg
    character (len=24) :: myname = 'getNexdKernelGreenTensor'
!
    call addTrace(errmsg,myname)
    nullify(g)
!
    if(this%jfcur == -1) then
       call add(errmsg,2,"this object does not contain any sensible values yet. "//&
            "call readFrequencyNexdKernelGreenTensor first",myname)
       return
    end if
!
    g => this%g
  end subroutine getNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorNexdKernelGreenTensor(this,uf_g,errmsg)
    type (nexd_kernel_green_tensor) :: this
    real :: uf_g
    type (error_message) :: errmsg
    if(this%staname /= '') then
       ! NEXD assumes SI units, so unit factor of green tensor is always 1.0
       uf_g = 1.0
    else
       ! object was not yet initiated, raise error
       uf_g = 0.
       call add(errmsg,2,"nexd_kernel_green_tensor object not yet initiated",&
            "getUnitFactorNexdKernelGreenTensor")
    end if
  end subroutine getUnitFactorNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get station name
!
  function getStanameNexdKernelGreenTensor(this) result(res)
    type (nexd_kernel_green_tensor), intent(in) :: this
    character(len=length_ID_nexd_kernel_green_tensor) :: res
    res = this%staname
  end function getStanameNexdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyNexdKernelGreenTensor(this,jf) result(next)
    type (nexd_kernel_green_tensor) :: this
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
  end function nextFrequencyNexdKernelGreenTensor
!
end module nexdKernelGreenTensor
