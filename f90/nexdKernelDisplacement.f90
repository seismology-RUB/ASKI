!----------------------------------------------------------------------------
!   Copyright 2015 Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module dealing with kernel displacements for method NEXD
!! \details The code of this module is based on module specfem3dKernelDisplacement
!!
!! \author Christian Ullisch
!! \date Nov 2015
!
module nexdKernelDisplacement
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateNexdKernelDisplacement; end interface
  interface operator (.id.); module procedure getEventIDNexdKernelDisplacement; end interface
  interface operator (.df.); module procedure getDfNexdKernelDisplacement; end interface
  interface operator (.nf.); module procedure getNfNexdKernelDisplacement; end interface
  interface operator (.jfcur.); module procedure getJfcurNexdKernelDisplacement; end interface
!
  integer, parameter :: length_ID_nexd_kernel_displacement = 13 !< change this number CONSISTENTLY with length_ASKI_output_ID in SPECFEM codes
!
!> \brief type containing all information on the kernel displacement files
  type nexd_kernel_displacement
     private
     character (len=500) :: basename = '' !< file base for all files of this kernel displacement object
     integer :: nproc = -1 !< number of parallel processors by which this file was created (safety feature, to check with the frequency files)
     real :: df = -1.   !< frequency step
     integer :: nfreq = 0  !< number of frequency indices, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null()     !< frequency indices
     integer :: jfcur = -1  !< indicates current frequency index for which wavefield values are contained in arrays u,ustr
     integer :: ntot = 0   !< total number of wavefield points (size(u,1), size(ustr,1))
     complex, dimension(:,:), pointer :: u => null()   !< dimension 1: number of wavefield points, dimension 2: wavefield component 1,..,3
     complex, dimension(:,:), pointer :: ustr => null() !< dimension 1: number of wavefield points, dimension 2: strain component 1,..,6
     integer :: lu = 0 !< file unit used for this object (hold this for the whole time between initialRead and deallocate)
  end type nexd_kernel_displacement
!
contains
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed from "main" file of this kernel displacement object
!
  subroutine initialReadNexdKernelDisplacement(this,fuh,basename,errmsg)
    type (nexd_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    integer :: ier
    character (len=33) :: myname = 'initialReadNexdKernelDisplacement'
!
    call addTrace(errmsg,myname)
!
    this%lu = get(fuh)
    this%basename = basename
!
    open(unit=this%lu,file=trim(this%basename)//'.main',status='old',action='read',access='stream',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(this%basename)//".main' cannot be opened to read",myname)
       goto 1
    endif

    read(this%lu) this%nproc,this%ntot,this%df,this%nfreq

    if(this%nfreq<1 .or. this%ntot<1) then
       write(errstr,*) "ntot,df,nf =",this%ntot,this%df,this%nfreq,"; ntot and nf should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    allocate(this%jf(this%nfreq))
    read(this%lu) this%jf

    close(this%lu)

    ! finally allocate for spectra
    allocate(this%u(this%ntot,3),this%ustr(this%ntot,6))

    ! if program comes here, there was no error reading the file, so return
    return

    ! if there is an error reading the file, deallocate everything before return
1   close(this%lu)
    call deallocateNexdKernelDisplacement(this,fuh)
  end subroutine initialReadNexdKernelDisplacement
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocateNexdKernelDisplacement(this,fuh)
    type (nexd_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    this%basename = ''
    this%nproc = -1
    this%df = -1.0d0
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
    this%ntot = 0
    if(associated(this%u)) deallocate(this%u)
    if(associated(this%ustr)) deallocate(this%ustr)
    call add(fuh,this%lu); this%lu = 0
  end subroutine deallocateNexdKernelDisplacement
!-----------------------------------------------------------------------------------------------------
!> \brief read specific file of displacement spectrum for the requested frequency (if existend)
!! \param jf frequency index for which kernel displacement should be read in (defines filename)
!
  subroutine readFrequencyNexdKernelDisplacement(this,jf,errmsg)
    type (nexd_kernel_displacement) :: this
    integer :: jf
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=35) :: myname = 'readFrequencyNexdKernelDisplacement'
    character(len=509) :: filename
    integer :: ier,nproc,ntot,jf_file
    real :: df
!
    call addTrace(errmsg,myname)
!
    ! check validity of frequency index
!
    if(.not.any(this%jf == jf)) then
       write(errstr,*) "incoming frequency index ",jf,", not contained in this kernelDisplacement object. "
       call add(errmsg,2,errstr,myname)
       return
    else
       write(errstr,*) "will read kernel displacement for frequency index ",jf
       call add(errmsg,0,errstr,myname)
    end if
!
    ! open specific file for this frequency
    write(filename,"(a,i6.6)") trim(this%basename)//".jf",jf
    open(unit=this%lu,file=filename,status='old',action='read',access='stream',&
         form='unformatted',iostat=ier)
    if (ier /= 0) then
       call add(errmsg,2,"File '"//trim(filename)//"' cannot be opened to read",myname)
       close(this%lu)
       return
    else
       call add(errmsg,0,"Successfully opened file '"//trim(filename)//"' to read",myname)
    endif
!
    read(this%lu) nproc,ntot,df,jf_file
    if(nproc/=this%nproc) then
       write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc," of main file"
       call add(errmsg,2,errstr,myname)
       close(this%lu)
       return
    end if
    if(ntot/=this%ntot) then
       write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot," of main file"
       call add(errmsg,2,errstr,myname)
       close(this%lu)
       return
    end if
    if(df/=this%df) then
       write(errstr,*) "df ",df," contained in this file does not match df ",this%df," of main file"
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

    ! make this%jf invalid before changing values in arrays this%u,this%ustr
    ! in case of returning on error (well, crashing the program in this case..), it can be no longer assured what is contained in the object
    ! when returning on success (end of subroutine), set this%jfcur properly
    this%jfcur = -1
    ! read spectrum from file
    read(this%lu) this%u,this%ustr

    close(this%lu)
!
    this%jfcur = jf
  end subroutine readFrequencyNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfNexdKernelDisplacement(this) result(res)
    type (nexd_kernel_displacement), intent(in) :: this
    real :: res
    res = this%df
  end function getDfNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequency indices nfreq
!
  function getNfNexdKernelDisplacement(this) result(res)
    type (nexd_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nfreq
  end function getNfNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurNexdKernelDisplacement(this) result(res)
    type (nexd_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%jfcur
  end function getJfcurNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get strains
!
  function getStrainsNexdKernelDisplacement(this) result(ustr)
    type (nexd_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: ustr
    ustr => this%ustr
  end function getStrainsNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsNexdKernelDisplacement(this,uf_ustr)
    type (nexd_kernel_displacement), intent(in) :: this
    real :: uf_ustr
    ! NEXD assumes SI units, so unit factor of strains is always 1.0
    uf_ustr = 1.0
  end subroutine getUnitFactorStrainsNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get displacements
!
  function getNexdKernelDisplacement(this) result(u)
    type (nexd_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: u
    u => this%u
  end function getNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of displacements
!
  subroutine getUnitFactorNexdKernelDisplacement(this,uf_u)
    type (nexd_kernel_displacement), intent(in) :: this
    real :: uf_u
    ! NEXD assumes SI units, so unit factor of displacement is always 1.0
    uf_u = 1.0
  end subroutine getUnitFactorNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get event ID
!
  function getEventIDNexdKernelDisplacement(this) result(res)
    type (nexd_kernel_displacement), intent(in) :: this
    character(len=length_ID_nexd_kernel_displacement) :: res
    ! AT THE MOMENT, RETURN A FAKE ID, THIS FEATURE IS NOT SUPPORTED AT THE MOMENT, BUT SHOULD WORK LIKE THIS IN ASKI
    res = 'fakeEventNEXD'
  end function getEventIDNexdKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyNexdKernelDisplacement(this,jf) result(next)
    type (nexd_kernel_displacement) :: this
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
  end function nextFrequencyNexdKernelDisplacement
!
end module nexdKernelDisplacement
