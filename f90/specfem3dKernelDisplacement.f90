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
!> \brief Module dealing with kernel displacements for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module specfem3dKernelDisplacement
  use specfem3dForASKIFiles
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateSpecfem3dKernelDisplacement; end interface
  interface operator (.id.); module procedure getEventIDSpecfem3dKernelDisplacement; end interface
  interface operator (.df.); module procedure getDfSpecfem3dKernelDisplacement; end interface
  interface operator (.nf.); module procedure getNfSpecfem3dKernelDisplacement; end interface
  interface operator (.jfcur.); module procedure getJfcurSpecfem3dKernelDisplacement; end interface
!
!> \brief type containing all information on the kernel displacement files
  type specfem3d_kernel_displacement
     private
     character (len=500) :: basename = '' !< file base for all files of this kernel displacement object
     character (len=length_ID_specfem3d_for_ASKI_files) :: eventid = ''  !< event id
     integer :: nproc = -1 !< number of parallel processors by which this file was created (safety feature, to check with the frequency files)
     integer :: specfem_version = -1 !< for safety reasons, remember this (check with other files)
     double precision :: df = -1.0d0   !< frequency step
     integer :: nfreq = 0  !< number of frequency indices, i.e. size(jf)
     integer, dimension(:), pointer :: jf => null()     !< frequency indices
     integer :: jfcur = -1  !< indicates current frequency index for which wavefield values are contained in arrays u,ustr
     integer :: ntot = 0   !< total number of wavefield points (size(u,1), size(ustr,1))
     complex, dimension(:,:), pointer :: u => null()   !< dimension 1: number of wavefield points, dimension 2: wavefield component 1,..,3
     complex, dimension(:,:), pointer :: ustr => null() !< dimension 1: number of wavefield points, dimension 2: strain component 1,..,6
     integer :: lu = 0 !< file unit used for this object (hold this for the whole time between initialRead and deallocate)
  end type specfem3d_kernel_displacement
!
contains
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed from "main" file of this kernel displacement object
!
  subroutine initialReadSpecfem3dKernelDisplacement(this,fuh,basename,errmsg)
    type (specfem3d_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=38) :: myname = 'initialReadSpecfem3dKernelDisplacement'
!
    call addTrace(errmsg,myname)
!
    this%lu = get(fuh)
    this%basename = basename
!
    call readSpecfem3dForASKIMainFile(trim(this%basename)//'.main',this%lu,errmsg,&
         specfem_version=this%specfem_version,file_ID=this%eventid,nproc=this%nproc,&
         nwp=this%ntot,df=this%df,nf=this%nfreq,jf=this%jf)
    if(.level.errmsg == 2) goto 2
!
    select case(this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files,version_cartesian_specfem3d_for_ASKI_files)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "readSpecfem3dForASKIMainFile returned specfem_version = ",this%specfem_version,&
            "; this module only supports specfem versions ",version_globe_specfem3d_for_ASKI_files," = SPECFEM3D_GLOBE and ",&
            version_cartesian_specfem3d_for_ASKI_files," = SPECFEM3D_Cartesian ."
       call add(errmsg,2,errstr,myname)
       goto 2
    end select
!
    ! if program comes here, there was no error reading the file, so return
    return
!
    ! if there is an error reading the file, deallocate everything before return
2   call deallocateSpecfem3dKernelDisplacement(this,fuh)
  end subroutine initialReadSpecfem3dKernelDisplacement
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
  subroutine deallocateSpecfem3dKernelDisplacement(this,fuh)
    type (specfem3d_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    this%basename = ''
    this%eventid = ''
    this%nproc = -1
    this%specfem_version = -1
    this%df = -1.0d0
    this%nfreq = 0
    if(associated(this%jf)) deallocate(this%jf)
    this%jfcur = -1
    this%ntot = 0
    if(associated(this%u)) deallocate(this%u)
    if(associated(this%ustr)) deallocate(this%ustr)
    call add(fuh,this%lu); this%lu = 0
  end subroutine deallocateSpecfem3dKernelDisplacement
!-----------------------------------------------------------------------------------------------------
!> \brief read specific file of displacement spectrum for the requested frequency (if existend)
!! \param jf frequency index for which kernel displacement should be read in (defines filename)
!
  subroutine readFrequencySpecfem3dKernelDisplacement(this,jf,errmsg)
    type (specfem3d_kernel_displacement) :: this
    integer :: jf
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=40) :: myname = 'readFrequencySpecfem3dKernelDisplacement'
    character(len=509) :: filename
    integer :: nproc,specfem_version,ntot,jf_file
    double precision :: df
    character(len=length_ID_specfem3d_for_ASKI_files) :: id
!
    call addTrace(errmsg,myname)
!
    if(this%lu == -1) then
       call add(errmsg,2,"object was not yet initiated, call initialReadSpecfem3dKernelDisplacement first",myname)
       return
    end if
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
    ! make this%jfcur invalid before changing values in arrays this%u,this%ustr.
    ! in case of returning on error, it can be no longer assured what is contained in the object.
    ! when returning on success (end of subroutine), set this%jfcur properly
    this%jfcur = -1
!
    ! read specific file for this frequency
    write(filename,"(a,i6.6)") trim(this%basename)//".jf",jf
    if(associated(this%u)) deallocate(this%u)
    if(associated(this%ustr)) deallocate(this%ustr)
    call readSpecfem3dForASKISpectralWavefieldFile(filename,this%lu,errmsg,specfem_version=specfem_version,&
         file_ID=id,nproc=nproc,nwp=ntot,df=df,jfcur=jf_file,u=this%u,ustr=this%ustr)
!
    ! check content of frequency file
    if(specfem_version/=this%specfem_version) then
       write(errstr,*) "specfem_version ",specfem_version," contained in this file does not match specfem_version ",&
            this%specfem_version," of main file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(id/=this%eventid) then
       call add(errmsg,2,"ID contained in this file '"//trim(id)//"' does not match ID '"//&
            trim(this%eventid)//"' of main file",myname)
       return
    end if
    if(nproc/=this%nproc) then
       write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc," of main file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(ntot/=this%ntot) then
       write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot," of main file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(df/=this%df) then
       write(errstr,*) "df ",df," contained in this file does not match df ",this%df," of main file"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(jf_file/=jf) then
       write(errstr,*) "requested frequency index ",jf," (also contained in filename) does not match the frequency index ",&
            jf_file,"actually contained in the file"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! check spectra
    if(.not.(associated(this%u).and.associated(this%ustr))) then
       call add(errmsg,2,"readSpecfem3dForASKISpectralWavefieldFile did not return both of u, ustr",myname)
       return
    end if
    if(size(this%u,1)/=this%ntot .or. size(this%u,2)/= 3) then
       write(errstr,*) "readSpecfem3dForASKISpectralWavefieldFile returned u of size ",&
            size(this%u,1),size(this%u,2),"; expected size ",this%ntot,"   3"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(this%ustr,1)/=this%ntot .or. size(this%ustr,2)/= 6) then
       write(errstr,*) "readSpecfem3dForASKISpectralWavefieldFile returned ustr of size ",&
            size(this%ustr,1),size(this%ustr,2),"; expected size ",this%ntot,"   6"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    this%jfcur = jf
  end subroutine readFrequencySpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfSpecfem3dKernelDisplacement(this) result(res)
    type (specfem3d_kernel_displacement), intent(in) :: this
    real :: res
    res = real(this%df)
  end function getDfSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequency indices nfreq
!
  function getNfSpecfem3dKernelDisplacement(this) result(res)
    type (specfem3d_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nfreq
  end function getNfSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index jfcur
!
  function getJfcurSpecfem3dKernelDisplacement(this) result(res)
    type (specfem3d_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%jfcur
  end function getJfcurSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get strains
!
  function getStrainsSpecfem3dKernelDisplacement(this) result(ustr)
    type (specfem3d_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: ustr
    ustr => this%ustr
  end function getStrainsSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsSpecfem3dKernelDisplacement(this,uf_ustr,errmsg)
    type (specfem3d_kernel_displacement) :: this
    real :: uf_ustr
    type(error_message) :: errmsg
    select case (this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files)
       uf_ustr = 1.0
    case(version_cartesian_specfem3d_for_ASKI_files)
       uf_ustr = 1.0
    case default ! object not yet initiated
       uf_ustr = 0.
       call add(errmsg,2,"specfem3d_kernel_displacement object not yet initiated. don't know the specfem version",&
            "getUnitFactorStrainsSpecfem3dKernelDisplacement")
    end select
  end subroutine getUnitFactorStrainsSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get displacements
!
  function getSpecfem3dKernelDisplacement(this) result(u)
    type (specfem3d_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: u
    u => this%u
  end function getSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of displacements
!
  subroutine getUnitFactorSpecfem3dKernelDisplacement(this,uf_u,errmsg)
    type (specfem3d_kernel_displacement) :: this
    real :: uf_u
    type(error_message) :: errmsg
    select case (this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files)
       uf_u = 1.0
    case(version_cartesian_specfem3d_for_ASKI_files)
       uf_u = 1.0
    case default ! object not yet initiated
       uf_u = 0.
       call add(errmsg,2,"specfem3d_kernel_displacement object not yet initiated. don't know the specfem version",&
            "getUnitFactorSpecfem3dKernelDisplacement")
    end select
  end subroutine getUnitFactorSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get event ID
!
  function getEventIDSpecfem3dKernelDisplacement(this) result(res)
    type (specfem3d_kernel_displacement), intent(in) :: this
    character(len=length_ID_specfem3d_for_ASKI_files) :: res
    res = this%eventid
  end function getEventIDSpecfem3dKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencySpecfem3dKernelDisplacement(this,jf) result(next)
    type (specfem3d_kernel_displacement) :: this
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
  end function nextFrequencySpecfem3dKernelDisplacement
!
end module specfem3dKernelDisplacement
