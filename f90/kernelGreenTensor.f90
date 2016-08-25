!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   and Christian Ullisch (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Wrapper module for kernel Green tensor
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  Green functions and strain arrays ordered consistently with ordering of wavefield points 
!
module kernelGreenTensor
  use geminiKernelGreenTensor
  use specfem3dKernelGreenTensor
  use nexdKernelGreenTensor
  use complexKernelFrequency
  use componentTransformation
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateKernelGreenTensor; end interface
  interface operator (.id.); module procedure getIdKernelGreenTensor; end interface
  interface operator (.df.); module procedure getDfKernelGreenTensor; end interface
  interface operator (.nf.); module procedure getNfKernelGreenTensor; end interface
  interface operator (.jfcur.); module procedure getJfcurKernelGreenTensor; end interface
  interface operator (.fc.); module procedure getComplexFrequencyKernelGreenTensor; end interface
!
!< fork here to specific method by associating some pointer (exactly ONE!)
  type kernel_green_tensor
     private
     integer :: ncomp = 0 !< number of components for which this green tensor object is initiated (could be only a few green functions, not the full tensor)
     character (len=character_length_component), dimension(:), pointer :: comp => null() !< vector of the ncomp components for which this object is initiated
     logical :: apply_component_transform = .false. !< logical, indicating whether to apply trans_coef or simply copy in get routines
     double precision, dimension(:,:), pointer :: trans_coef => null()            !< dependent on flag apply_component_transform, this pointer is associated (containing the transformation matrix) or not
!
     type (gemini_kernel_green_tensor), pointer :: gemini_kgt => null()         !< gemini specific kernel green tensor object
     type (specfem3d_kernel_green_tensor), pointer :: specfem3d_kgt => null()     !< specfem3d specific kernel green tensor object
     type (nexd_kernel_green_tensor), pointer :: nexd_kgt => null()     !< nexd specific kernel green tensor object
  end type kernel_green_tensor
!
  integer, parameter :: length_ID_kernel_green_tensor = 13 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
!
contains
!------------------------------------------------------------------------
!> \brief Initiate kernel Green tensor
!
  subroutine initiateKernelGreenTensor(this,method,fuh,filename,comp,comptrans,staname,errmsg)
    type (kernel_green_tensor) :: this
    character(len=*) :: method
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    character(len=*), dimension(:) :: comp
    type (component_transformation) :: comptrans
    character(len=*) :: staname
    type (error_message) :: errmsg
    ! local
    character (len=25) :: myname = 'initiateKernelGreenTensor'
    character(len=400) :: errstr
    integer :: ncomp,ncomp_available,icomp,i
    character (len=character_length_component), dimension(:), pointer :: comp_available
!
    nullify(comp_available)
    call addTrace(errmsg,myname)
    write(errstr,*) "initiating Green tensor object for station '",trim(staname),"' for components '"!,trim(comp)//",","'"
    do i = 1,size(comp)
       errstr = trim(errstr)//trim(comp(i))//","
    end do ! i
    errstr = trim(errstr)//"'"
    call add(errmsg,0,errstr,myname)
!
    ! check if this object is already initiated. If so, return an error. It could be problematic to simply deallocate 
    ! it here and just give a warning, since fuh might be inconsistent. The programmer/user should fix this problem 
    ! in the routines calling here.
    if(this%ncomp /= 0) then
       call add(errmsg,2,"This kernel_green_tensor_object is already initiated. Please deallocate this object first.",myname)
       return
    end if
!
    ! check if vector comp is valid
    ncomp = size(comp)
    if(ncomp <= 0) then
       call add(errmsg,2,"there are no incoming components",myname)
       return
    end if ! ncomp <= 0
    if(.not.allValidComponents(comp,i_invalid=icomp)) then
       write(errstr,*) icomp,"'th incoming component '"//trim(comp(icomp))//"' not valid. Valid components are '"//&
            all_valid_components//"'"
       call add(errmsg,2,errstr,myname)
       return
    end if ! .not.allValid
    do icomp = 2,ncomp
       if(any(comp(1:icomp-1)==comp(icomp))) then
          write(errstr,*) icomp,"'th requested component '",trim(comp(icomp)),"' occurs more than once in the list of the ",&
               ncomp," requested components. there must not be duplicate components in the request"
       end if
    end do ! icomp
!
    select case (method)
    case('GEMINI')
       allocate(this%gemini_kgt)
       call initialReadGeminiKernelGreenTensor(this%gemini_kgt,fuh,filename,errmsg)
       if (.level.errmsg == 2) then
          call addTrace(errmsg,myname); return
       endif
       ! if this object is initated for CX CY CZ (in exactly this order), we don't need
       ! to  initiate trans_coef! in this case, we can simply copy values below, as
       ! this%gemini_kgt always provides complete tensor with components CX,CY,CZ (in that order)
       ! so check this here first
       if(ncomp /= 3) then
          this%apply_component_transform = .true.
       else
          this%apply_component_transform = comp(1)/='CX' .or. comp(2)/='CY' .or. comp(3)/='CZ'
       end if
       if(this%apply_component_transform) this%trans_coef => transform(comptrans,(/'CX','CY','CZ'/),comp,staname)
!
    case('SPECFEM3D')
       allocate(this%specfem3d_kgt)
       call getAvailableComponentsSpecfem3dKernelGreenTensor(filename,get(fuh),comp_available,errmsg)
       call undo(fuh)
       if(.not.associated(comp_available)) return
       write(errstr,*) "Green tensor components available by SPECFEM3D are : '",comp_available//",","'"
       call add(errmsg,0,errstr,myname)
       if(.not.allValidComponents(comp_available,i_invalid=icomp)) then
          write(errstr,*) icomp,"'th component '"//trim(comp(icomp))//&
               "' available by SPECFEM3D not valid. Valid components are '"//all_valid_components//"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if ! .not.allValid

       ! if there are more than 3 available components, probably the user does not know what he's doing, so raise an error
       ncomp_available = size(comp_available)
       if(ncomp_available > 3) then
          write(errstr,*) "number of components available by SPECFEM3D = ",ncomp_available,&
               ". In practice, any number > 3 does not make sense. Please make sure of what you're doing"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       ! check if requested components are all contained in vector comp_available
       if(allComponentsContained(comp,comp_available)) then
          ! if so, simply initiate specfem3d_kernel_green_tensor object with requested components comp
          ! otherwise a component transformation must be applied
          call initialReadSpecfem3dKernelGreenTensor(this%specfem3d_kgt,fuh,filename,staname,comp,errmsg)
          if (.level.errmsg == 2) goto 1
          this%apply_component_transform = .false.
       else ! allComponentsContained
          ! in case we need to apply a component transformation, raise a warning and ask the user to 
          ! doublecheck if the requested components can really be produced from the available ones 
          ! (this check would be very intense to do here, so leave it to the user, he should know what he's doing)
          if(ncomp_available < 3) then
             write(errstr,*) "there are only ",ncomp_available," component(s) available by SPECFEM3D, namely: '",&
                  comp_available//",",&
                  "'. PLEASE BE SURE THAT THOSE REALLY CAN BE TRANSFORMED WITHOUT LOSS OF INFORMATION INTO ",&
                  "THE REQUESTED COMPONENTS, WHICH ARE: '",comp//",","'. IF NOT, YOU SHOULD NOT TRUST THE RESULTS!"
             call add(errmsg,1,errstr,myname)
          else
             write(errstr,*) "there are ",ncomp_available," component(s) available by SPECFEM3D, namely: '",comp_available//",",&
                  "'. please be sure that those really can be transformed without loss of information into ",&
                  "the requested components, which are: '",comp//",","'. if not, you should not trust the results!"
             call add(errmsg,1,errstr,myname)
          end if
          call initialReadSpecfem3dKernelGreenTensor(this%specfem3d_kgt,fuh,filename,staname,comp_available,errmsg)
          if (.level.errmsg == 2) goto 1
          this%apply_component_transform = .true.
          this%trans_coef => transform(comptrans,comp_available,comp,staname)
       end if ! allComponentsContained
!
    case('NEXD')
       allocate(this%nexd_kgt)
       call getAvailableComponentsNexdKernelGreenTensor(filename,get(fuh),comp_available,errmsg)
       call undo(fuh)
       if(.not.associated(comp_available)) return
       write(errstr,*) "Green tensor components available by NEXD are : '",comp_available//",","'"
       call add(errmsg,0,errstr,myname)
       if(.not.allValidComponents(comp_available,i_invalid=icomp)) then
          write(errstr,*) icomp,"'th component '"//trim(comp(icomp))//&
               "' available by NEXD not valid. Valid components are '"//all_valid_components//"'"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if ! .not.allValid

       ! if there are more than 3 available components, probably the user does not know what he's doing, so raise an error
       ncomp_available = size(comp_available)
       if(ncomp_available > 3) then
          write(errstr,*) "number of components available by NEXD = ",ncomp_available,&
               ". In practice, any number > 3 does not make sense. Please make sure of what you're doing"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if

       ! check if requested components are all contained in vector comp_available
       if(allComponentsContained(comp,comp_available)) then
          ! if so, simply initiate nexd_kernel_green_tensor object with requested components comp
          ! otherwise a component transformation must be applied
          call initialReadNexdKernelGreenTensor(this%nexd_kgt,fuh,filename,staname,comp,errmsg)
          if (.level.errmsg == 2) goto 1
          this%apply_component_transform = .false.
       else ! allComponentsContained
          ! in case we need to apply a component transformation, raise a warning and ask the user to 
          ! doublecheck if the requested components can really be produced from the available ones 
          ! (this check would be very intense to do here, so leave it to the user, he should know what he's doing)
          if(ncomp_available < 3) then
             write(errstr,*) "there are only ",ncomp_available," component(s) available by NEXD, namely: '",&
                  comp_available//",",&
                  "'. PLEASE BE SURE THAT THOSE REALLY CAN BE TRANSFORMED WITHOUT LOSS OF INFORMATION INTO ",&
                  "THE REQUESTED COMPONENTS, WHICH ARE: '",comp//",","'. IF NOT, YOU SHOULD NOT TRUST THE RESULTS!"
             call add(errmsg,1,errstr,myname)
          else
             write(errstr,*) "there are ",ncomp_available," component(s) available by NEXD, namely: '",comp_available//",",&
                  "'. please be sure that those really can be transformed without loss of information into ",&
                  "the requested components, which are: '",comp//",","'. if not, you should not trust the results!"
             call add(errmsg,1,errstr,myname)
          end if
          call initialReadNexdKernelGreenTensor(this%nexd_kgt,fuh,filename,staname,comp_available,errmsg)
          if (.level.errmsg == 2) goto 1
          this%apply_component_transform = .true.
          this%trans_coef => transform(comptrans,comp_available,comp,staname)
       end if ! allComponentsContained
!
    case default
       call add(errmsg,2,"Invalid forward computation method '"//trim(method)//"'",myname)
       return
    end select
!
    this%ncomp = ncomp
    allocate(this%comp(this%ncomp))
    this%comp = comp
!
    return
!
1   if(associated(comp_available)) deallocate(comp_available)
  end subroutine initiateKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Deallocate complete kernel_green_tensor object
!
  subroutine deallocateKernelGreenTensor(this,fuh)
    type (kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
!
    if (associated(this%gemini_kgt)) call dealloc(this%gemini_kgt,fuh)
    if (associated(this%specfem3d_kgt)) call dealloc(this%specfem3d_kgt,fuh)
    if (associated(this%nexd_kgt)) call dealloc(this%nexd_kgt,fuh)
    if(associated(this%comp)) deallocate(this%comp)
    this%ncomp = 0
    if(associated(this%trans_coef)) deallocate(this%trans_coef)
  end subroutine deallocateKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief check if this kernel_green_tensor object has the same components (and order) as the given ones
!
  function hasExactComponentsKernelGreenTensor(this,comp) result(l)
    type (kernel_green_tensor) :: this
    character(len=*), dimension(:) :: comp
    logical :: l
    ! local
    integer :: icomp,ncomp
!
    ! initiate to false. then simply return if any check is not successfull
    l = .false.
    ncomp = size(comp)
!
    if(ncomp /= this%ncomp) return
!
    do icomp = 1,ncomp ! this even works if ncomp == 0 (then this loop does nothing)
       if(comp(icomp) /= this%comp(icomp)) return
    end do ! icomp
!
    l = .true.
  end function hasExactComponentsKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Read in kernel green tensor for one frequency
!
  subroutine readFrequencyKernelGreenTensor(this,jf,errmsg)
    type (kernel_green_tensor) :: this
    integer :: jf
    type (error_message) :: errmsg
    character (len=30) :: myname = 'readFrequencyKernelGreenTensor'
!
    call addTrace(errmsg,myname)
    if (associated(this%gemini_kgt)) then
       call readFrequencyGeminiKernelGreenTensor(this%gemini_kgt,jf,errmsg)
       if (.level.errmsg == 2) return
    else if (associated(this%specfem3d_kgt)) then
       call readFrequencySpecfem3dKernelGreenTensor(this%specfem3d_kgt,jf,errmsg)
       if (.level.errmsg == 2) return
    else if (associated(this%nexd_kgt)) then
       call readFrequencyNexdKernelGreenTensor(this%nexd_kgt,jf,errmsg)
       if (.level.errmsg == 2) return
    endif
  end subroutine readFrequencyKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get ID
!
  function getIdKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    character(len=length_ID_kernel_green_tensor) :: res
!
    res = ''
    if (associated(this%gemini_kgt)) then
       res = .id.(this%gemini_kgt)
    else if (associated(this%specfem3d_kgt)) then
       res = .id.(this%specfem3d_kgt)
    else if (associated(this%nexd_kgt)) then
       res = .id.(this%nexd_kgt)
    endif
  end function getIdKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency step df
!
  function getDfKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    real :: res
!
    res = 0.
    if (associated(this%gemini_kgt)) then
       res = .df.(this%gemini_kgt)
    else if (associated(this%specfem3d_kgt)) then
       res = .df.(this%specfem3d_kgt)
    else if (associated(this%nexd_kgt)) then
       res = .df.(this%nexd_kgt)
    endif
  end function getDfKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
  function getNfKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    integer :: res
!
    res = 0
    if (associated(this%gemini_kgt)) then
       res = .nf.(this%gemini_kgt)
    else if (associated(this%specfem3d_kgt)) then
       res = .nf.(this%specfem3d_kgt)
    else if (associated(this%nexd_kgt)) then
       res = .nf.(this%nexd_kgt)
    endif
  end function getNfKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency fcur
!
  function getJfcurKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    integer :: res
!
    res = 0
    if (associated(this%gemini_kgt)) then
       res = .jfcur.(this%gemini_kgt)
    else if (associated(this%specfem3d_kgt)) then
       res = .jfcur.(this%specfem3d_kgt)
    else if (associated(this%nexd_kgt)) then
       res = .jfcur.(this%nexd_kgt)
    endif
  end function getJfcurKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get complex frequency, return negative real frequency if not yet initiated
!
  function getComplexFrequencyKernelGreenTensor(this) result(res)
    type (kernel_green_tensor), intent(in) :: this
    complex :: res
    ! local
    real :: df
    integer :: jf
!
    res = cmplx(-1.0)
    if (associated(this%gemini_kgt)) then
       df = .df.(this%gemini_kgt)
       jf = .jfcur.(this%gemini_kgt)
       res = getComplexKernelFrequency('GEMINI',df,jf)
    else if (associated(this%specfem3d_kgt)) then
       df = .df.(this%specfem3d_kgt)
       jf = .jfcur.(this%specfem3d_kgt)
       res = getComplexKernelFrequency('SPECFEM3D',df,jf)
    else if (associated(this%nexd_kgt)) then
       df = .df.(this%nexd_kgt)
       jf = .jfcur.(this%nexd_kgt)
       res = getComplexKernelFrequency('NEXD',df,jf)
    endif
  end function getComplexFrequencyKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get pointer to green tensor strains.
!! \details Allocate here for the return array gstr (allocation can also be done in the
!!  method-specific realization of kernelGreenTensor). This is done in order to be most flexible,
!!  not enforcing a method to store a pointer of for gstr inside its kernelGreenTensor object. 
!!  Routines calling getStrainsKernelGreenTensor must care for deallocation of pointer gstr!
!
  subroutine getStrainsKernelGreenTensor(this,gstr,errmsg)
    type (kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: gstr
    type (error_message) :: errmsg
    character (len=27) :: myname = 'getStrainsKernelGreenTensor'
    ! local
    complex, dimension(:,:,:), pointer :: gstr_local
    integer :: icomp_out,icomp_local
!
    nullify(gstr_local)
!
    call addTrace(errmsg,myname)
    nullify(gstr)
!
    if (associated(this%gemini_kgt)) then
       ! gstr_local must not be deallocated here!
       gstr_local => getStrainsGeminiKernelGreenTensor(this%gemini_kgt)
       if(associated(gstr_local)) then
          ! gstr: dimension 1 = wavefield points, dimension 2 = strain component, dimension 3 = force component of green function
          allocate(gstr(size(gstr_local,1),size(gstr_local,2),this%ncomp))
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation do the transformation by hand; 
             ! matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                ! since we defined this%trans_coef => transform(comptrans,(/'CX','CY','CZ'/),comp,.id.this%gemini_kgt), 
                ! we know here that size(this%trans_coef,2) = 3 corresponds to CX,CY,CZ in 3rd dimension of gstr_local
                ! Therefore, do the matrix-vector multiplication of the transformation by hand as:
                gstr(:,:,icomp_out) = this%trans_coef(icomp_out,1)*gstr_local(:,:,1) + &
                     this%trans_coef(icomp_out,2)*gstr_local(:,:,2) + this%trans_coef(icomp_out,3)*gstr_local(:,:,3)
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the object was initiated for components 'CX','CY','CZ' which are returned be GEMINI in that order.
             ! so here, we simply need to copy (component transformation would be identity)
             gstr = gstr_local
          end if ! this%apply_component_transform
       else ! associated(gstr_local)
          call add(errmsg,2,"no Green tensor strains were returned by function getStrainsGeminiKernelGreenTensor",myname)
          return
       end if ! associated(gstr_local)

    else if (associated(this%specfem3d_kgt)) then
       ! gstr_local must not be deallocated here!
       call getStrainsSpecfem3dKernelGreenTensor(this%specfem3d_kgt,gstr_local,errmsg)
       if(associated(gstr_local)) then
          allocate(gstr(size(gstr_local,1),size(gstr_local,2),this%ncomp))
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation, do the transformation by hand; 
             ! matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                gstr(:,:,icomp_out) = (0.,0.)
                do icomp_local = 1,size(gstr_local,3)
                   gstr(:,:,icomp_out) = gstr(:,:,icomp_out) + &
                        this%trans_coef(icomp_out,icomp_local)*gstr_local(:,:,icomp_local)
                end do ! icomp_local
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the SPECFEM3D object was initiated for the requested components this%comp (in that order)
             ! so here, we simply need to copy (as component transformation would be the identity)
             gstr = gstr_local
          end if ! this%apply_component_transform
       else ! associated(gstr_local)
          call add(errmsg,2,"no Green tensor strains were returned by function getStrainsSpecfem3dKernelGreenTensor",myname)
          return
       end if ! associated(gstr_local)

    else if (associated(this%nexd_kgt)) then
       ! gstr_local must not be deallocated here!
       call getStrainsNexdKernelGreenTensor(this%nexd_kgt,gstr_local,errmsg)
       if(associated(gstr_local)) then
          allocate(gstr(size(gstr_local,1),size(gstr_local,2),this%ncomp))
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation, do the transformation by hand; 
             ! matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                gstr(:,:,icomp_out) = (0.,0.)
                do icomp_local = 1,size(gstr_local,3)
                   gstr(:,:,icomp_out) = gstr(:,:,icomp_out) + &
                        this%trans_coef(icomp_out,icomp_local)*gstr_local(:,:,icomp_local)
                end do ! icomp_local
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the NEXD object was initiated for the requested components this%comp (in that order)
             ! so here, we simply need to copy (as component transformation would be the identity)
             gstr = gstr_local
          end if ! this%apply_component_transform
       else ! associated(gstr_local)
          call add(errmsg,2,"no Green tensor strains were returned by function getStrainsNexdKernelGreenTensor",myname)
          return
       end if ! associated(gstr_local)

    else ! associated this%***_kgt
       call add(errmsg,2,"kernel_green_tensor object not yet initiated",myname)
       return
    endif ! associated this%***_kgt
  end subroutine getStrainsKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsKernelGreenTensor(this,uf_gstr,errmsg)
    type (kernel_green_tensor) :: this
    real :: uf_gstr
    type (error_message) :: errmsg
    if (associated(this%gemini_kgt)) then
       call getUnitFactorStrainsGeminiKernelGreenTensor(this%gemini_kgt,uf_gstr)
    else if (associated(this%specfem3d_kgt)) then
       call getUnitFactorStrainsSpecfem3dKernelGreenTensor(this%specfem3d_kgt,uf_gstr,errmsg)
    else if (associated(this%nexd_kgt)) then
       call getUnitFactorStrainsNexdKernelGreenTensor(this%nexd_kgt,uf_gstr,errmsg)
    else
       call add(errmsg,2,"kernel_green_tensor object not yet initiated","getUnitFactorStrainsKernelGreenTensor")
       return
    end if
  end subroutine getUnitFactorStrainsKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get values of green tensor
!! \details Allocate here for the return array g (allocation can also be done in the
!!  method-specific realization of kernelGreenTensor). This is done in order to be most flexible,
!!  not enforcing a method to store a pointer of for g inside its kernelGreenTensor object. 
!!  Routines calling getKernelGreenTensor must care for deallocation of pointer g!
!
  subroutine getKernelGreenTensor(this,g,errmsg)
    type (kernel_green_tensor) :: this
    complex, dimension(:,:,:), pointer :: g
    type (error_message) :: errmsg
    character (len=20) :: myname = 'getKernelGreenTensor'
    ! local
    complex, dimension(:,:,:), pointer :: g_local
    integer :: icomp_out,icomp_local
!
    nullify(g_local)
!
    call addTrace(errmsg,myname)
    nullify(g)
!
    if (associated(this%gemini_kgt)) then
       ! g_local must not be deallocated here!
       g_local => getGeminiKernelGreenTensor(this%gemini_kgt)
       if(associated(g_local)) then
          ! g_local : dimension 1 = wavefield points, dimension 2 = wavefield component, dimension 3 = force component of green function
          allocate(g(size(g_local,1),size(g_local,2),this%ncomp)) 
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation for GEMINI in this module,
             ! do the transformation by hand, matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                ! since we defined this%trans_coef => transform(comptrans,(/'CX','CY','CZ'/),comp,.id.this%gemini_kgt), 
                ! we know here that size(this%trans_coef,2) = 3 corresponds to CX,CY,CZ in 3rd dimension of g_local
                ! Therefore, do the matrix-vector multiplication of the transformation by hand as:
                g(:,:,icomp_out) = this%trans_coef(icomp_out,1)*g_local(:,:,1) + &
                     this%trans_coef(icomp_out,2)*g_local(:,:,2) + this%trans_coef(icomp_out,3)*g_local(:,:,3)
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the object was initiated for components 'CX','CY','CZ' which are returned be GEMINI in that order.
             ! so here, we simply need to copy (component transformation would be identity)
             g = g_local
          end if ! this%apply_component_transform
       else ! associated(this%gemini_kgt)
          call add(errmsg,2,"no Green tensor values were returned by function getGeminiKernelGreenTensor",myname)
          return
       end if ! associated(this%gemini_kgt)

    else if (associated(this%specfem3d_kgt)) then
       ! g_local must not be deallocated here!
       call getSpecfem3dKernelGreenTensor(this%specfem3d_kgt,g_local,errmsg)
       if(associated(g_local)) then
          allocate(g(size(g_local,1),size(g_local,2),this%ncomp)) 
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation, do the transformation by hand; 
             ! matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                g(:,:,icomp_out) = (0.,0.)
                do icomp_local = 1,size(g_local,3)
                   g(:,:,icomp_out) = g(:,:,icomp_out) + &
                        this%trans_coef(icomp_out,icomp_local)*g_local(:,:,icomp_local)
                end do ! icomp_local
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the SPECFEM3D object was initiated for the requested components this%comp (in that order)
             ! so here, we simply need to copy (as component transformation would be the identity)
             g = g_local
          end if ! this%apply_component_transform
       else ! associated(this%specfem3d_kgt)
          call add(errmsg,2,"no Green tensor strains were returned by function getSpecfem3dKernelGreenTensor",myname)
          return
       end if ! associated(this%specfem3d_kgt)

    else if (associated(this%nexd_kgt)) then
       ! g_local must not be deallocated here!
       call getNexdKernelGreenTensor(this%nexd_kgt,g_local,errmsg)
       if(associated(g_local)) then
          allocate(g(size(g_local,1),size(g_local,2),this%ncomp)) 
          if(this%apply_component_transform) then
             ! in case we need to apply the component transformation, do the transformation by hand; 
             ! matmul is not possible because of more than two dimensions
             do icomp_out = 1,this%ncomp
                g(:,:,icomp_out) = (0.,0.)
                do icomp_local = 1,size(g_local,3)
                   g(:,:,icomp_out) = g(:,:,icomp_out) + &
                        this%trans_coef(icomp_out,icomp_local)*g_local(:,:,icomp_local)
                end do ! icomp_local
             end do ! icomp_out
          else ! this%apply_component_transform
             ! in this case, the NEXD object was initiated for the requested components this%comp (in that order)
             ! so here, we simply need to copy (as component transformation would be the identity)
             g = g_local
          end if ! this%apply_component_transform
       else ! associated(this%nexd_kgt)
          call add(errmsg,2,"no Green tensor strains were returned by function getNexdKernelGreenTensor",myname)
          return
       end if ! associated(this%nexd_kgt)

    else ! associated this%***_kgt
       call add(errmsg,2,"kernel_green_tensor object not yet initiated",myname)
       return
    endif ! associated(this%***_kgt)
  end subroutine getKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of green functions
!
  subroutine getUnitFactorKernelGreenTensor(this,uf_g,errmsg)
    type (kernel_green_tensor) :: this
    real :: uf_g
    type (error_message) :: errmsg
    if (associated(this%gemini_kgt)) then
       call getUnitFactorGeminiKernelGreenTensor(this%gemini_kgt,uf_g)
    else if (associated(this%specfem3d_kgt)) then
       call getUnitFactorSpecfem3dKernelGreenTensor(this%specfem3d_kgt,uf_g,errmsg)
    else if (associated(this%nexd_kgt)) then
       call getUnitFactorNexdKernelGreenTensor(this%nexd_kgt,uf_g,errmsg)
    else
       call add(errmsg,2,"kernel_green_tensor object not yet initiated","getUnitFactorKernelGreenTensor")
       return
    end if
  end subroutine getUnitFactorKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!
  function nextFrequencyKernelGreenTensor(this,jf) result(next)
    type (kernel_green_tensor) :: this
    integer :: jf
    logical :: next
    if (associated(this%gemini_kgt)) then
       next = nextFrequencyGeminiKernelGreenTensor(this%gemini_kgt,jf)
    else if (associated(this%specfem3d_kgt)) then
       next = nextFrequencySpecfem3dKernelGreenTensor(this%specfem3d_kgt,jf)
    else if (associated(this%nexd_kgt)) then
       next = nextFrequencyNexdKernelGreenTensor(this%nexd_kgt,jf)
    endif
  end function nextFrequencyKernelGreenTensor
!
end module kernelGreenTensor
