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
!> \brief Module dealing with kernel green tensor for methods SPECFEM3D (Cartesian and GLOBE)
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module specfem3dKernelGreenTensor
  use specfem3dForASKIFiles
  use componentTransformation
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateSpecfem3dKernelGreenTensor; end interface
  interface operator (.id.); module procedure getStanameSpecfem3dKernelGreenTensor; end interface
  interface operator (.df.); module procedure getDfSpecfem3dKernelGreenTensor; end interface
  interface operator (.nf.); module procedure getNfSpecfem3dKernelGreenTensor; end interface
  interface operator (.jfcur.); module procedure getJfcurSpecfem3dKernelGreenTensor; end interface
!
!> \brief type containing all information on the kernel green tensor files
  type specfem3d_kernel_green_tensor
     private
     ! basic information
     character (len=498) :: basename = '' !< file base for all files of this kernel green tensor object
     character (len=length_ID_specfem3d_for_ASKI_files) :: staname = '' !< station name (string before '_' in the ID contained in the file
     integer :: nproc = -1 !< number of parallel processors by which this file was created (safety feature, to check with the frequency files)
     integer :: specfem_version = -1 !< for safety reasons, remember this (check with other files)
     integer :: lu = -1 !< file unit used for this object (hold this for the whole time between initialRead and deallocate)
!
     ! frequency content
     double precision :: df = -1.0d0   !< frequency step
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
  end type specfem3d_kernel_green_tensor
!
contains
!----------------------------------------------------------------------------------------------------
!> \brief check which components are available at all
!
  subroutine getAvailableComponentsSpecfem3dKernelGreenTensor(basename,lu,comp_available,errmsg)
    character(len=*) :: basename
    integer :: lu
    character (len=character_length_component), dimension(:), pointer :: comp_available
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=48) :: myname = 'getAvailableComponentsSpecfem3dKernelGreenTensor'
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
  end subroutine getAvailableComponentsSpecfem3dKernelGreenTensor
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed from "main" file of this kernel green tensor object
!
  subroutine initialReadSpecfem3dKernelGreenTensor(this,fuh,basename,staname,comp,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename,staname
    character(len=*), dimension(:) :: comp
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character (len=508) :: filename
    integer :: icomp,type_invgrid_1
    integer :: j,specfem_version,nproc,type_invgrid,ntot,nf
    double precision :: df
    integer, dimension(:), pointer :: jf
    character (len=length_ID_specfem3d_for_ASKI_files) :: staname_comp
    character (len=37) :: myname = 'initialReadSpecfem3dKernelGreenTensor'
!
    nullify(jf)
!
    ! IN MODULE kernelGreenTensor IT WAS CHECKED THAT comp IS NOT EMPTY AND ALL ENTRIES ARE 
    ! VALID COMPONENTS AND THAT ALL ENTRIES IN VECTOR comp ARE AVAILABLE COMPONENTS FOR THIS OBJECT
!
    call addTrace(errmsg,myname)
    write(errstr,*) "initiating Specfem3d Green tensor object for station '",trim(staname),"' for components '"!,comp//",","'"
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
    ! getAvailableComponentsSpecfem3dKernelGreenTensor which reads in file this%basename//'.comp' which 
    ! contains the components of Green functions which were written from SPECFEM3D and are 
    ! availabe to be read in. In module kernelGreenTensor it was already decided for which components
    ! this object should be initiated (entries in incoming vector comp)
    ! Here, there is no need anymore to decide if those can be transformed to the requested components 
    ! (which are not known to module specfem3dKernelGreenTensor)
    this%ncomp = size(comp)
    allocate(this%comp(this%ncomp))
    this%comp = comp
!
    ! FIRST FILE, green function component this%comp(1)
!
    filename = trim(this%basename)//'_'//trim(this%comp(1))//'.main'
    call readSpecfem3dForASKIMainFile(filename,this%lu,errmsg,&
         specfem_version=this%specfem_version,file_ID=staname_comp,nproc=this%nproc,&
         type_inversion_grid=type_invgrid_1,nwp=this%ntot,df=this%df,nf=this%nfreq,jf=this%jf)
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
    if(staname_comp /= trim(staname)//'_'//trim(this%comp(1))) then
       write(errstr,*) "ASKI output ID read from file = '",trim(staname_comp),&
            "' has not the expected value 'staname_comp' = '",trim(staname)//'_'//trim(this%comp(1)),&
            "' referring to station name '",trim(staname),"' and component '",trim(this%comp(1)),&
            "' of this Green function. By convention, ASKI output ID must be of this form."
       call add(errmsg,2,errstr,myname)
       goto 2
    end if
    this%staname = staname
!
    ! NOW CHECK, IF THE MAIN FILES OF THE OTHER COMPONENTS CONTAIN THE SAME INFORMATION!
    do icomp = 2,this%ncomp ! if this%ncomp == 1, this loop does nothing (so works for all cases)
!
       filename = trim(this%basename)//'_'//trim(this%comp(icomp))//'.main'
       call readSpecfem3dForASKIMainFile(filename,this%lu,errmsg,&
            specfem_version=specfem_version,file_ID=staname_comp,nproc=nproc,&
            type_inversion_grid=type_invgrid,nwp=ntot,df=df,nf=nf,jf=jf)
       if(.level.errmsg == 2) goto 2
!
       if(specfem_version/=this%specfem_version) then
          write(errstr,*) "specfem_version ",specfem_version," contained in this file does not match specfem_version ",&
               this%specfem_version," contained in .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
!
       if(staname_comp /= trim(staname)//'_'//trim(this%comp(icomp))) then
          write(errstr,*) "ASKI output ID read from file = '",trim(staname_comp),&
               "' has not the expected value 'staname_comp' = '",trim(staname)//'_'//trim(this%comp(icomp)),&
               "' referring to station name '",trim(staname),"' and component '",trim(this%comp(icomp)),&
               "' of this Green function. By convention, ASKI output ID must be of this form."
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
!
       if(nproc/=this%nproc) then
          write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc,&
               " of the .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(type_invgrid/=type_invgrid_1) then
          write(errstr,*) "type of inversion grid ",type_invgrid," contained in this file does not match type ",&
               type_invgrid_1," in .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(ntot/=this%ntot) then
          write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot,&
               " of .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(df/=this%df) then
          write(errstr,*) "df ",df," contained in this file does not match df ",this%df,&
               " of .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
       if(nf/=this%nfreq) then
          write(errstr,*) "nf ",nf," contained in this file does not match nf ",this%nfreq,&
               " of .main file of first component '",trim(this%comp(1)),"'"
          call add(errmsg,2,errstr,myname)
          goto 2
       end if
!
       do j=1,nf
          if(.not.any(jf==this%jf(j))) then
             write(errstr,*) "frequency index ",this%jf(j),", contained in file of first component '",&
                  trim(this%comp(1)),"' is NOT contained in this file!"
             call add(errmsg,2,errstr,myname)
             goto 2
          end if 
       end do ! j
       deallocate(jf)
!
    end do ! icomp
!
    if(associated(jf)) deallocate(jf)
!
    ! Finally allocate for spectra. 
    ! By contrast to module specfem3dKernelDisplacement, the arrays this%g, this%gstr are filled by deep 
    ! copies in readFrequencySpecfem3dKernelGreenTensor, not by linking
    allocate(this%g(this%ntot,3,this%ncomp),this%gstr(this%ntot,6,this%ncomp))
!
    ! if program comes here, there was no error reading the files and jf is deallocated, so return
    return
!
    ! if there is an error reading the files, deallocate everything before return
2   if(associated(jf)) deallocate(jf)
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
    this%ncomp = 0
    if(associated(this%comp)) deallocate(this%comp)
    this%ntot = 0
    if(associated(this%g)) deallocate(this%g)
    if(associated(this%gstr)) deallocate(this%gstr)
    call add(fuh,this%lu); this%lu = -1
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
    character(len=512) :: filename
    integer :: icomp,nproc,specfem_version,ntot,jf_file
    double precision :: df
    character(len=length_ID_specfem3d_for_ASKI_files) :: staname_comp
    complex, dimension(:,:), pointer :: g_comp,gstr_comp
!
    nullify(g_comp,gstr_comp)
!
    call addTrace(errmsg,myname)
!
    ! check if this was initiated
    if(this%lu == -1) then
       call add(errmsg,2,"Green tensor object not yet initiated, call initialReadSpecfem3dKernelGreenTensor first",myname)
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
    ! read specific files for this frequency
!
    do icomp = 1,this%ncomp
       
       ! read file for this component
       write(filename,"(a,i6.6)") trim(this%basename)//"_"//trim(this%comp(icomp))//".jf",jf
       if(associated(g_comp)) deallocate(g_comp)
       if(associated(gstr_comp)) deallocate(gstr_comp)
       call readSpecfem3dForASKISpectralWavefieldFile(filename,this%lu,errmsg,specfem_version=specfem_version,&
            file_ID=staname_comp,nproc=nproc,nwp=ntot,df=df,jfcur=jf_file,u=g_comp,ustr=gstr_comp)
!
       ! check content of file
       if(specfem_version/=this%specfem_version) then
          write(errstr,*) "specfem_version ",specfem_version," contained in this file does not match specfem_version ",&
               this%specfem_version," of main files"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(staname_comp/=trim(this%staname)//'_'//trim(this%comp(icomp))) then
          write(errstr,*) "ASKI output ID read from file = '",trim(staname_comp),&
               "' has not the expected value 'staname_comp' = '",trim(this%staname)//'_'//trim(this%comp(icomp)),&
               "' referring to station name '",trim(this%staname),"' and component '",trim(this%comp(icomp)),&
               "' of this Green function. By convention, ASKI output ID must be of this form."
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(nproc/=this%nproc) then
          write(errstr,*) "nproc ",nproc," contained in this file does not match nproc ",this%nproc," of main files"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(ntot/=this%ntot) then
          write(errstr,*) "ntot ",ntot," contained in this file does not match ntot ",this%ntot," of main files"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(df/=this%df) then
          write(errstr,*) "df ",df," contained in this file does not match df ",this%df," of main files"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(jf_file/=jf) then
          write(errstr,*) "requested frequency index ",jf," (also contained in filename) does not match the frequency index ",&
               jf_file," which is actually contained in the file"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       ! check spectra
       if(.not.(associated(g_comp).and.associated(gstr_comp))) then
          call add(errmsg,2,"readSpecfem3dForASKISpectralWavefieldFile did not return both of g, gstr",myname)
          goto 1
       end if
       if(size(g_comp,1)/=this%ntot .or. size(g_comp,2)/= 3) then
          write(errstr,*) "readSpecfem3dForASKISpectralWavefieldFile returned g of size ",&
               size(g_comp,1),size(g_comp,2),"; expected size ",this%ntot,"   3"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(size(gstr_comp,1)/=this%ntot .or. size(gstr_comp,2)/= 6) then
          write(errstr,*) "readSpecfem3dForASKISpectralWavefieldFile returned gstr of size ",&
               size(gstr_comp,1),size(gstr_comp,2),"; expected size ",this%ntot,"   6"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       this%g(:,:,icomp) = g_comp
       this%gstr(:,:,icomp) = gstr_comp
       deallocate(g_comp,gstr_comp)
!
    end do ! icomp
!
    ! if code comes here, everything went alright, so indicate the current frequency contained in this object
    this%jfcur = jf
!
1   if(associated(g_comp)) deallocate(g_comp)
    if(associated(gstr_comp)) deallocate(gstr_comp)
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
!> \brief Get strains at specific components. Raise error if incoming components do not match those for which this tensor object was initiated
!
  subroutine getStrainsSpecfem3dKernelGreenTensor(this,gstr,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: gstr
    type (error_message) :: errmsg
    character (len=36) :: myname = 'getStrainsSpecfem3dKernelGreenTensor'
!
    nullify(gstr)
!
    call addTrace(errmsg,myname)
!
    if(this%jfcur == -1) then
       call add(errmsg,2,"this object does not contain any sensible values yet. "//&
            "call readFrequencySpecfem3dKernelGreenTensor first",myname)
       return
    end if
!
    gstr => this%gstr
  end subroutine getStrainsSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsSpecfem3dKernelGreenTensor(this,uf_gstr,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    real :: uf_gstr
    type(error_message) :: errmsg
    select case (this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files)
       ! presumably there is a bug (?) in form of a nearly planar wave arriving from the deep earth
       ! (observed in a 1chunk simulation). The amplitude and waveform of that wave is INDEPENDENT of 
       ! the exciting seismic source. When applying a single force with 1 N, the resulting aplitudes
       ! of the seismic waves are VERY small, whereby the strange artefact "wave" becomes highly 
       ! significant and contaminates the actual waveforms tremendously.
       ! (Florian Schumacher is still unsure whether this is a "bug" which is known to the developer of SPECFEM3D)
       ! In order to account for this phenomenom, we boost the amplitude from 1 Newton to 1e12 Newton.
       ! We do NOT divide by 1e12 before writing ASKI output to file, but instead apply a unit factor of 1e-12 
       ! in the ASKI main package for kernel green tensor outpu. 
       ! Since the amplitude of the strange artefact signal 
       ! is independent of the applied seismic source, this removes the artefact from the waveforms.
       uf_gstr = 1.e-12
    case(version_cartesian_specfem3d_for_ASKI_files)
       uf_gstr = 1.0
    case default ! object not yet initiated
       uf_gstr = 0.
       call add(errmsg,2,"specfem3d_kernel_green_tensor object not yet initiated. don't know the specfem version",&
            "getUnitFactorStrainsSpecfem3dKernelGreenTensor")
    end select
  end subroutine getUnitFactorStrainsSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get Green tensor at specific components. Raise error if incoming components do not match those for which this tensor object was initiated
!
  subroutine getSpecfem3dKernelGreenTensor(this,g,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: g
    type (error_message) :: errmsg
    character (len=29) :: myname = 'getSpecfem3dKernelGreenTensor'
!
    nullify(g)
!
    call addTrace(errmsg,myname)
!
    if(this%jfcur == -1) then
       call add(errmsg,2,"this object does not contain any sensible values yet. "//&
            "call readFrequencySpecfem3dKernelGreenTensor first",myname)
       return
    end if
!
    g => this%g
  end subroutine getSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of green functions
!
  subroutine getUnitFactorSpecfem3dKernelGreenTensor(this,uf_g,errmsg)
    type (specfem3d_kernel_green_tensor) :: this
    real :: uf_g
    type(error_message) :: errmsg
    select case (this%specfem_version)
    case(version_globe_specfem3d_for_ASKI_files)
       ! presumably there is a bug (?) in form of a nearly planar wave arriving from the deep earth
       ! (observed in a 1chunk simulation). The amplitude and waveform of that wave is INDEPENDENT of 
       ! the exciting seismic source. When applying a single force with 1 N, the resulting aplitudes
       ! of the seismic waves are VERY small, whereby the strange artefact "wave" becomes highly 
       ! significant and contaminates the actual waveforms tremendously.
       ! (Florian Schumacher is still unsure whether this is a "bug" which is known to the developer of SPECFEM3D)
       ! In order to account for this phenomenom, we boost the amplitude from 1 Newton to 1e12 Newton.
       ! We do NOT divide by 1e12 before writing ASKI output to file, but instead apply a unit factor of 1e-12 
       ! in the ASKI main package for kernel green tensor output. 
       ! Since the amplitude of the strange artefact signal 
       ! is independent of the applied seismic source, this removes the artefact from the waveforms.
       uf_g = 1.e-12
    case(version_cartesian_specfem3d_for_ASKI_files)
       uf_g = 1.0
    case default ! object not yet initiated
       uf_g = 0.
       call add(errmsg,2,"specfem3d_kernel_green_tensor object not yet initiated. don't know the specfem version",&
            "getUnitFactorSpecfem3dKernelGreenTensor")
    end select
  end subroutine getUnitFactorSpecfem3dKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get station name
!
  function getStanameSpecfem3dKernelGreenTensor(this) result(res)
    type (specfem3d_kernel_green_tensor), intent(in) :: this
    character(len=length_ID_specfem3d_for_ASKI_files) :: res
    res = this%staname
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
