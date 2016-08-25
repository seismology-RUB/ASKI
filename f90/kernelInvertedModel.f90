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
!> \brief this module manages the inverted Earth model on the inversion grid
!!
!! \details At the end of each iteration step of the ASKI kernel inversion, the
!!  resulting updated (inverted) model can computed (by adding perturbations to
!!  the starting model of that iteration step) and written to file by this module.
!!  Parmetrizations other than isotropic velocities are not supported yet. It may
!!  be in general difficult to extend ASKI to other parametrizations. 
!!
!! \author Florian Schumacher
!! \date Nov 2015
!
module kernelInvertedModel
!
  use modelParametrization
  use vectorPointer
  use inversionGrid
  use wavefieldPoints
  use integrationWeights
  use kernelReferenceModel
  use invgridVtkFile
  use dataModelSpaceInfo
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateKernelInvertedModel; end interface
  interface init; module procedure initiateKernelInvertedModel; end interface
!  interface add; module procedure addValuesKernelInvertedModel; end interface ! "add" is in conflict with some other interface
  interface getVal; module procedure getValuesKernelInvertedModel; end interface
  interface getIndx; module procedure getIndicesKernelInvertedModel; end interface
  interface maxValue; module procedure getMaxValueKernelInvertedModel; end interface
  interface minValue; module procedure getMinValueKernelInvertedModel; end interface
  interface sum
     module procedure summateInstancesKernelInvertedModel
     module procedure summateKernelReferenceAndKernelInvertedModel
  end interface sum
!
  interface operator (.pmtrz.); module procedure getParametrizationKernelInvertedModel; end interface
!
!> \brief type which holds model values (on a selection of inversion grid cells) for some specific model parametrization
  type kernel_inverted_model
     private
     character(len=character_length_pmtrz) :: parametrization = '' !< model parametrization in use; empty string = object not initialized
     !
     type (real_vector_pointer), dimension(:), pointer :: model_values => null() !< 1D model value arrays, one array per parameter
     !
     type (integer_vector_pointer), dimension(:), pointer :: indx => null() !< for each parameter, array of inversion grid cell indices to which model_values belong
  end type kernel_inverted_model
!
contains
!------------------------------------------------------------------------
!> \brief initiate kernel inverted model object
!! \details just allocate and initiate the inversion grid cell index array
!!  and define parametrization. Specific values for parameters are added by 
!!  routine addValuesKernelInvertedModel
!! \param this kernel inverted model object
!! \param parametrization parametrization of model (at the moment 'isoVelocity', 'isoElasticity' supported)
!! \param errmsg error message
!! \return error message
!
  subroutine initiateKernelInvertedModel(this,parametrization,errmsg)
    type (kernel_inverted_model) :: this
    character(len=*) :: parametrization
    type (error_message) :: errmsg
    character(len=27) :: myname = 'initiateKernelInvertedModel'
    !
    call addTrace(errmsg,myname)
    !
    if(trim(this%parametrization) /= '') then
       call add(errmsg,1,'this object was already allocated. deallocating it now before initiating anew',myname)
       call deallocateKernelInvertedModel(this)
    endif
    !
    if(.not.validModelParametrization(parametrization)) then
       call add(errmsg,2,"parametrization '"//trim(parametrization)//"' not supported, "//&
            "valid parametrizations (parameters) are: "//all_valid_pmtrz_param,myname)
       return
    end if
    this%parametrization = parametrization
    allocate(this%model_values(numberOfParamModelParametrization(parametrization)),&
         this%indx(numberOfParamModelParametrization(parametrization)))
  end subroutine initiateKernelInvertedModel
!------------------------------------------------------------------------
!> \brief deallocate kernel inverted model object
!! \param this kernel inverted model object
!
  subroutine deallocateKernelInvertedModel(this)
    type (kernel_inverted_model) :: this
    integer :: i
    this%parametrization = ''
    if(associated(this%model_values)) then
       do i = 1,size(this%model_values)
          call dealloc(this%model_values(i))
       enddo
       deallocate(this%model_values)
    endif
    if(associated(this%indx)) then
       do i = 1,size(this%indx)
          call dealloc(this%indx(i))
       enddo
       deallocate(this%indx)
    endif
  end subroutine deallocateKernelInvertedModel
!------------------------------------------------------------------------
!> \brief allocate and add model values for one specific parameter
!! \param this kernel inverted model object
!! \param param parameter to be added
!! \param indx array of inversion grid cell indices for which model values will be added
!! \param modelval_param array of model values for that specific parameter
!! \param errmsg error message
!! \return error message
!
  subroutine addValuesKernelInvertedModel(this,param,indx,modelval_param,errmsg)
    type (kernel_inverted_model) :: this
    character(len=*) :: param
    integer, dimension(:) :: indx
    real, dimension(:) :: modelval_param
    type (error_message) :: errmsg
    character(len=300) :: errstr
    character(len=28) :: myname = 'addValuesKernelInvertedModel'
    integer :: iparam,size_indx,size_modelval_param
    !
    call addTrace(errmsg,myname)
    !
    if(trim(this%parametrization) == '') then
       call add(errmsg,2,"object not initiated yet. please call initiateKernelInvertedModel first",myname)
    endif
    !
    if(.not.validParamModelParametrization(this%parametrization,param)) then
       call add(errmsg,2,"'"//trim(param)//"' is an invalid parameter in parametrization '"//&
            trim(this%parametrization)//"'. Valid parametrizations (parameters) : "//all_valid_pmtrz_param,myname)
       return
    endif
    iparam = indexOfParamModelParametrization(this%parametrization,param)
    !
    size_indx = size(indx); size_modelval_param = size(modelval_param)
    if(.not.size_indx>0) then
       call add(errmsg,2,'no incoming inversion grid cell indices',myname)
       return
    endif
    if(any(indx<1)) then
       call add(errmsg,2,'there are invalid, i.e. non-positive, incoming inversion grid cell indices',myname)
       return
    endif
    if(.not.size_modelval_param>0) then
       call add(errmsg,2,"no incoming model parameter values",myname)
       return
    endif
    if(size_indx/=size_modelval_param) then
       write(errstr,*) "number of incoming inversion grid indices (",size_indx, &
            ") does not match number of incoming model parameter values",size_modelval_param,")"
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    !
    if(associated(getVectorPointer(this%model_values(iparam)))) then
       call dealloc(this%model_values(iparam))
       call add(errmsg,1,"there were already model values for '"//trim(this%parametrization)//"' parameter '"//&
            trim(param)//"'. overwriting the old values now",myname)
    endif
    if(associated(getVectorPointer(this%indx(iparam)))) call dealloc(this%indx(iparam))
    !
    call allocateVectorPointer(this%indx(iparam),size_indx)
    call fillVectorPointer(this%indx(iparam),indx,1)
    call allocateVectorPointer(this%model_values(iparam),size_modelval_param)
    call fillVectorPointer(this%model_values(iparam),modelval_param,1)
    !
  end subroutine addValuesKernelInvertedModel
!------------------------------------------------------------------------
  function getParametrizationKernelInvertedModel(this) result(pmtrz)
    type (kernel_inverted_model), intent(in) :: this
    character(len=character_length_pmtrz) :: pmtrz
    pmtrz = this%parametrization
  end function getParametrizationKernelInvertedModel
!------------------------------------------------------------------------
  function getValuesKernelInvertedModel(this,param) result(p)
    type (kernel_inverted_model) :: this
    character(len=*) :: param
    real, dimension(:), pointer :: p
    integer :: iparam
    nullify(p)
    if(trim(this%parametrization)=='') return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    iparam = indexOfParamModelParametrization(this%parametrization,param)
    p => getVectorPointer(this%model_values(iparam))
  end function getValuesKernelInvertedModel
!------------------------------------------------------------------------
  function getIndicesKernelInvertedModel(this,param) result(p)
    type (kernel_inverted_model) :: this
    character(len=*) :: param
    integer, dimension(:), pointer :: p
    integer :: iparam
    nullify(p)
    if(trim(this%parametrization)=='') return
    iparam = indexOfParamModelParametrization(this%parametrization,param)
    if(.not.iparam>0) return
    p => getVectorPointer(this%indx(iparam))
  end function getIndicesKernelInvertedModel
!------------------------------------------------------------------------
  function getModelValuesVectorPointerKernelInvertedModel(this) result(vp)
    type (kernel_inverted_model) :: this
    type (real_vector_pointer), dimension(:), pointer :: vp
    nullify(vp)
    vp => this%model_values
  end function getModelValuesVectorPointerKernelInvertedModel
!------------------------------------------------------------------------
  function getIndexVectorPointerKernelInvertedModel(this) result(vp)
    type (kernel_inverted_model) :: this
    type (integer_vector_pointer), dimension(:), pointer :: vp
    nullify(vp)
    vp => this%indx
  end function getIndexVectorPointerKernelInvertedModel
!------------------------------------------------------------------------
  function getMaxValueKernelInvertedModel(this,param) result(max)
    type (kernel_inverted_model) :: this
    character(len=*) :: param
    real :: max
    integer :: iparam
    real, dimension(:), pointer :: p
    nullify(p)
    max = -huge(1.0) ! initiate to invalid value
    if(trim(this%parametrization)=='') return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    iparam = indexOfParamModelParametrization(this%parametrization,param)
    p => getVectorPointer(this%model_values(iparam))
    if(.not.associated(p)) return
    max = maxval(p)
  end function getMaxValueKernelInvertedModel
!------------------------------------------------------------------------
  function getMinValueKernelInvertedModel(this,param) result(min)
    type (kernel_inverted_model) :: this
    character(len=*) :: param
    real :: min
    integer :: iparam
    real, dimension(:), pointer :: p
    nullify(p)
    min = huge(1.0) ! initiate to invalid value
    if(trim(this%parametrization)=='') return
    if(.not.validParamModelParametrization(this%parametrization,param)) return
    iparam = indexOfParamModelParametrization(this%parametrization,param)
    p => getVectorPointer(this%model_values(iparam))
    if(.not.associated(p)) return
    min = minval(p)
  end function getMinValueKernelInvertedModel
!------------------------------------------------------------------------
!> \brief update this by adding model values of that, if there are common parameters
!! \param this model object which will be modified/updated (changed on output)
!! \param that model object which will be added to this (not changed on output)
!! \param c1 optional factor with which this is multiplied
!! \param c2 optional factor with which that is multiplied
!! \param relative indicates whether after (weighted) summation, the resulting values should be divided through by values of this
!! \param errmsg error message
!! \return error message
!
  subroutine summateInstancesKernelInvertedModel(this,that,errmsg,c1,c2,relative)
    type (kernel_inverted_model) :: this,that
    real, optional :: c1,c2
    logical, optional :: relative
    type (error_message) :: errmsg
    ! local
    character(len=35) :: myname = 'summateInstancesKernelInvertedModel'
    integer :: iparam,i,j,max_icell_in_this_indx,ival_append,sze,n
    integer, dimension(:), allocatable :: this_modelval_map,this_indx_append,indx_tmp
    real, dimension(:), allocatable :: this_model_values_append
    integer, dimension(:), pointer :: this_indx,that_indx
    real, dimension(:), pointer :: this_model_values,that_model_values
    real :: val
    logical :: compute_relative
!
    nullify(this_indx,that_indx,this_model_values,that_model_values)
!
    call addTrace(errmsg,myname)
!
    if(this%parametrization=='') then
       call add(errmsg,2,'model object, which is to be updated, is not initiated yet',myname)
       return
    endif
    if(that%parametrization=='') then
       call add(errmsg,2,'model object, which is to be added, is not initiated yet',myname)
       return
    endif
    if(this%parametrization/=that%parametrization) then
       call add(errmsg,2,"incoming model objects have different parametrizations. cannot add '"//&
            trim(that%parametrization)//"' parameters to '"//trim(this%parametrization)//"' parameters",myname)
       return
    endif
!
    if(present(relative)) then
       compute_relative = relative
    else
       compute_relative = .false.
    endif
!
    ! loop over all parameters of this (that has same parameters, as this and that have same parametrization)
    do iparam = 1,numberOfParamModelParametrization(this%parametrization)
       ! if there are no model values in that for current parameter, there is nothing to be added, 
       ! so check if c1 and  cycle
       if(.not.associated(getVectorPointer(that%model_values(iparam)))) then
          this_model_values => getVectorPointer(this%model_values(iparam))
          if(compute_relative) then
             if(present(c1)) then
                this_model_values = c1
             else
                this_model_values = 1.
             endif
          else
             if(present(c1)) this_model_values = c1*this_model_values
          endif
          cycle
       endif
!
       ! if there are no model values yet in this for current parameter, interpret them as zero 
       ! and associate model values of that as the new model values of this (if not compute_relative!)
       if(.not.associated(getVectorPointer(this%model_values(iparam)))) then
          ! if compute_relative, DO NOTHING (cannot "relate" somenthing to nothing, would be "infinity")
          if(compute_relative) cycle
!
          call addValuesKernelInvertedModel(this,getParamFromIndexModelParametrization(this%parametrization,iparam),&
               getVectorPointer(that%indx(iparam)),getVectorPointer(that%model_values(iparam)),errmsg)
          if(.level.errmsg == 2) return
!
          if(present(c2)) then
             ! then that_model_values should have been scaled by c2 before adding to this
             ! hence, multiply now (after adding) this_model_values by c2
             this_model_values => getVectorPointer(this%model_values(iparam))
             this_model_values = c2*this_model_values
          endif
!
          cycle
       endif
!
       ! if loop comes here, there are model values in both objects (this,that) for current parameter.
       ! have to assure here to summate values only at the same inversion grid cells.
       ! append values for inversion grid cells which are new in this. just use value of that (assuming value 0. if it does not exist)
!
       ! create index mapping for model values in this to be able to compare with values in that
       this_indx => getVectorPointer(this%indx(iparam))
       max_icell_in_this_indx = maxval(this_indx)
       allocate(this_modelval_map(max_icell_in_this_indx)); this_modelval_map = -1
       do i = 1,size(this_indx)
          ! for each invgrid cell indx j present in this%indx, set modelval_indx(j) = i
          ! if modelval_indx(j) = -1 for some invgrid cell indx j, it indicates that j is not present in this%indx
          ! as this%indx>0 (checked in addValues), there is no problem with accessing array this_modelval_map
!
          this_modelval_map(this_indx(i)) = i
       enddo
       !
       ! now loop over all model values in that and decide if (and to which index) it is to be added to this model values, 
       ! or if it is to be appended to this model values because the respective inversion grid cell is not present in this
       that_indx => getVectorPointer(that%indx(iparam))
       that_model_values => getVectorPointer(that%model_values(iparam))
       this_model_values => getVectorPointer(this%model_values(iparam))
       allocate(this_model_values_append(size(that_model_values)),this_indx_append(size(that_indx)))
       ival_append = 0
       do j = 1,size(that_indx)
          if(that_indx(j)<=max_icell_in_this_indx) then
             i = this_modelval_map(that_indx(j))
             if(i>0) then
                ! inversion grid cell index that_indx(j) is also present in this (at modelval index i), 
                ! so simply sum values
                val = that_model_values(j)
                if(present(c2)) val = c2*val
                if(present(c1)) then
                   val = val + c1*this_model_values(i)
                else
                   val = val + this_model_values(i)
                endif
                if(compute_relative) val = val/this_model_values(i)
                ! overwrite this_model_values(i)
                this_model_values(i) = val
                !this_model_values(i) = this_model_values(i) + that_model_values(j)
!
                ! remember i'th entry of array this_model_vales was treated here. 
                ! if, in the end, there were any invgrid cell indices which were only
                ! present in this, but not in that, (see below) then you have to treat those according
                ! to present(c1), compute_relative, so remeber cell that_indx(j) was treated here
                this_modelval_map(that_indx(j)) = -1
                cycle
             endif
          endif
          ! if loop comes here, either that_indx(i)>max_icell_in_this_indx, or modelval_indx(that_indx(i)) negative
!
          ! if compute_relative, here is also the problem that we cannot relate something to nothing
          ! (no model value in this, but model value in that)
          if(compute_relative) cycle ! go to next index j

          ! otherwise, append that model value
          ival_append = ival_append + 1
          this_indx_append(ival_append) = that_indx(j)
          if(present(c2)) then
             this_model_values_append(ival_append) = c2*that_model_values(j)
          else
             this_model_values_append(ival_append) = that_model_values(j)
          endif
       enddo ! j
!
       ! check if there were any invgrid cell indices which were only
       ! present in this, but not present in that (i.e. which were not treated in j-loop above)
       ! for all those indices i, this_modelval_map(i) is still positive.
       ! this is important in order to respect present(c1) and compute_relative correctly.
       if(any(this_modelval_map>0) .and. (compute_relative.or.present(c1))) then
          n = count(this_modelval_map>0); allocate(indx_tmp(n))
          indx_tmp = pack(this_modelval_map,this_modelval_map>0)
!
          if(compute_relative) then
             if(present(c1)) then
                this_model_values(indx_tmp) = c1
             else
                this_model_values(indx_tmp) = 1.
             endif
          else
             if(present(c1)) this_model_values(indx_tmp) = c1*this_model_values(indx_tmp)
          endif
!
          deallocate(indx_tmp)
       endif
!
       ! reallocate this_indx,this_model_values if there were any values appended
       if(ival_append > 0) then
          sze = size(this_indx)
          this_indx => reallocate(this_indx,sze+ival_append)
          call associateVectorPointer(this%indx(iparam),this_indx)
          call fillVectorPointer(this%indx(iparam),this_indx_append(1:ival_append),sze+1)
!
          sze = size(this_model_values)
          this_model_values => reallocate(this_model_values,sze+ival_append)
          call associateVectorPointer(this%model_values(iparam),this_model_values)
          call fillVectorPointer(this%model_values(iparam),this_model_values_append(1:ival_append),sze+1)
       endif
!
       ! clean up
       deallocate(this_modelval_map,this_model_values_append,this_indx_append)
    enddo ! iparam
  end subroutine summateInstancesKernelInvertedModel
!------------------------------------------------------------------------
!> \brief add kernel reference model (interpolated on invgrid) to this
!! \param this inverted model object to which reference model will be added
!! \param krm kernel reference model which will be added to this inverted model object
!! \param invgrid inversion grid on which the inverted model values live (indices stored in this are assumed to relate to invgrid)
!! \param intw integration weights object which contain localizations of wavefield points (relating to krm) in invgrid
!! \param errmsg error message
!! \return error message
!
  subroutine summateKernelReferenceAndKernelInvertedModel(this,krm,invgrid,intw,errmsg)
    type (kernel_inverted_model) :: this
    type (kernel_reference_model) :: krm
    type (inversion_grid) :: invgrid
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    type (kernel_inverted_model) :: kim_krm
    character(len=44) :: myname = 'summateKernelReferenceAndKernelInvertedModel'
    !
    call addTrace(errmsg,myname)
    !
    if(this%parametrization=='') then
       call add(errmsg,2,'model object, which is to be updated, is not initiated yet',myname)
       return
    endif
    !
    call interpolateKernelReferenceToKernelInvertedModel(kim_krm,krm,this%parametrization,invgrid,intw,errmsg)
    if(.level.errmsg == 2) then
       call deallocateKernelInvertedModel(kim_krm)
       return
    endif
    !
    call summateInstancesKernelInvertedModel(this,kim_krm,errmsg)
    if(.level.errmsg == 2) then
       call deallocateKernelInvertedModel(kim_krm)
       return
    endif
    !
    call deallocateKernelInvertedModel(kim_krm)
  end subroutine summateKernelReferenceAndKernelInvertedModel
!------------------------------------------------------------------------
!> \brief interpolate kernel reference model (on wavefield points) to inversion grid
!! \details using average "integration weights" (TYPE_INTEGRATION_WEIGHTS = 0), the kernel reference model values
!!  for requested parametrization (given by value parametrization), which live on the wavefield points, are interpolated
!!  (averaged) in each inversion grid cell. The result ist added as model values to a newly created instance
!!  of type kernel_inverted_model. As there is no interface to module kernelReferenceModel which works with
!!  character strings of parametrization and parameters (as suggested to have in TODO_waveformInversion.txt), here
!!  we have to manually select a supported type of parametrization. 
!!  IN THE FUTURE: just check if incoming value for parametrization is valid in terms of this module (otherwise we cannot 
!!  create a new instance of type kernel_inverted_model here) and request the respective values by character strings from
!!  module kernelReferenceModel.
!! \param krm kernel reference model which will be added to this inverted model object
!! \param invgrid inversion grid on which the inverted model values live (indices stored in this are assumed to relate to invgrid)
!! \param intw integration weights object which contain localizations of wavefield points (relating to krm) in invgrid
!! \param errmsg error message
!! \return error message
!
  subroutine interpolateKernelReferenceToKernelInvertedModel(this,krm,parametrization,invgrid,intw,errmsg)
    type (kernel_inverted_model) :: this
    type (kernel_reference_model) :: krm
    character(len=*) :: parametrization
    type (inversion_grid) :: invgrid
    type (integration_weights) :: intw
    type (error_message) :: errmsg
    ! local
    character(len=47) :: myname = 'interpolateKernelReferenceToKernelInvertedModel'
    character(len=300) :: errstr
    character(len=character_length_param) :: param_name
    real, dimension(:), pointer :: model_values_wp
    real, dimension(:), allocatable :: model_values_invgrid
    integer, dimension(:), allocatable :: indx_invgrid
    integer, dimension(:), pointer :: wp_idx
    real, dimension(:), pointer :: weights
    integer :: icell,ncell,ntot_invgrid
!
    nullify(model_values_wp,wp_idx,weights)
!
    call addTrace(errmsg,myname)
!
    call initiateKernelInvertedModel(this,parametrization,errmsg)
    if(.level.errmsg == 2) then
       call deallocateKernelInvertedModel(this)
       goto 1
    endif
    !
    ntot_invgrid = .ncell.invgrid
    !
    ! loop over all parameters 
    do while(nextParamModelParametrization(this%parametrization,param_name))
       !
       model_values_wp => getModelValuesKernelReferenceModel(krm,this%parametrization,param_name)
       if(.not.associated(model_values_wp)) then
          call add(errmsg,1,"there are no kernel reference model values for '"//&
               trim(this%parametrization)//"' parameter '"//trim(param_name)//"'",myname)
          cycle
       end if
       !
       allocate(model_values_invgrid(ntot_invgrid),indx_invgrid(ntot_invgrid))
       ncell = 0 ! count valid entries in arrays model_values_invgrid,indx_invgrid (skip invalid invgrid cells)
       do icell = 1,ntot_invgrid
          ! check if this cell is empty. if yes, cycle
          if(emptyCell(intw,icell)) cycle
          !
          ! get wp indices for this cell
          wp_idx => intw.wpidx.icell
          if(.not.associated(wp_idx)) then
             write(errstr,*) "pointer to array of wavefield point indices for invgrid cell ",icell,&
                  " is not associated, which usually means that invgrid index is out of range. This could "//&
                  "mean your integration weights are corrupted or not compatible with this inversion grid"
             call add(errmsg,2,trim(errstr),myname)
             goto 1
          endif
          !
          ! get weights for this cell
          weights => intw.weight.icell
          if(.not.associated(weights)) then
             write(errstr,*) "pointer to array of weights for invgrid cell ",icell,&
                  " is not associated, which usually means that invgrid index is out of range. This could "//&
                  "mean your integration weights are corrupted or not compatible with this inversion grid"
             call add(errmsg,2,trim(errstr),myname)
             goto 1
          endif
          !
          if(size(wp_idx) /= size(weights)) then
             write(errstr,*) "there are ",size(wp_idx)," wavefield point indices but ",size(weights),&
                  " weights for invgrid cell ",icell,"; this means your integration weights are corrupted"
             call add(errmsg,2,trim(errstr),myname)
             goto 1
          end if
          !
          ! append averaged model value to array model_values_invgrid and memorize this invgrid cell index
          ncell = ncell + 1
          !model_values_invgrid(ncell) = sum(model_values_wp(wp_idx))/size(wp_idx) ! OLD IMPLEMENTATION: JUST AVERAGE (depreciated, better use integration weights)
          model_values_invgrid(ncell) = sum(weights*model_values_wp(wp_idx))/sum(weights)
          indx_invgrid(ncell) = icell
       enddo ! icell
       !
       ! now add all valid model values found (and respective invgrid indices) as model values of current parameter to this object
       call addValuesKernelInvertedModel(this,param_name,indx_invgrid(1:ncell),model_values_invgrid(1:ncell),errmsg)
       deallocate(model_values_invgrid,indx_invgrid)
       if(.level.errmsg == 2) goto 1
       !
       if(associated(model_values_wp)) deallocate(model_values_wp)
    end do ! while (nextParam)
    !
1   if(associated(model_values_wp)) deallocate(model_values_wp)
    if(allocated(model_values_invgrid)) deallocate(model_values_invgrid)
    if(allocated(indx_invgrid)) deallocate(indx_invgrid)
  end subroutine interpolateKernelReferenceToKernelInvertedModel
!------------------------------------------------------------------------
!> \brief read kernel inverted model from file
!! \param this kernel inverted model object
!! \param lu file unit
!! \param filename filename of model file
!! \param errmsg error message
!! \return error message
!
  subroutine readFileKernelInvertedModel(this,filename,lu,errmsg)
    type (kernel_inverted_model) :: this
    integer :: lu
    character(len=*) :: filename
    type (error_message) :: errmsg
    character(len=300) :: errstr
    character(len=27) :: myname = 'readFileKernelInvertedModel'
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name
    integer :: ios,nval
    integer, dimension(:), allocatable :: indx
    real, dimension(:), allocatable :: val
    !
    call addTrace(errmsg,myname)
    !
    ! open file to read
    open(unit=lu,file=filename,form='unformatted',status='old',action='read',access='stream',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "opening unformatted file '"//trim(filename)//"' to read raised iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    call add(errmsg,0,"successfully opened file '"//trim(filename)//"' to read",myname)
    !
    read(lu,iostat=ios) parametrization
    if(ios/=0) then
       write(errstr,*) 'trying to read parametrization from file raised iostat = ',ios
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    call initiateKernelInvertedModel(this,parametrization,errmsg)
    if(.level.errmsg == 2) then
       close(lu)
       call deallocateKernelInvertedModel(this)
       return
    endif
    call add(errmsg,0,"model parametrization contained in file is '"//trim(parametrization)//"'",myname)
    !
    do while (ios==0)
       read(lu,iostat=ios) param_name,nval
       if(ios==0) then
          if(nval>0) then
             allocate(indx(nval),val(nval))
             read(lu,iostat=ios) indx,val
             if(ios/=0) then
                write(errstr,*) 'trying to read ',nval,"'",trim(this%parametrization),"'-'",trim(param_name),&
                     "' model values and their indices raised iostat = ",ios
                call add(errmsg,2,trim(errstr),myname)
                close(lu)
                return
             endif
             call addValuesKernelInvertedModel(this,param_name,indx,val,errmsg)
             if(.level.errmsg == 2) then
                close(lu)
                return
             endif
             deallocate(indx,val)
          else
             write(errstr,*) "negative number of model values (",nval,") in file"
             call add(errmsg,2,trim(errstr),myname)
             close(lu)
             return
          endif
       endif
    enddo ! while (ios==0)
    !
    close(lu)
  end subroutine readFileKernelInvertedModel
!------------------------------------------------------------------------
!> \brief write kernel inverted model to file
!! \param this kernel inverted model object
!! \param lu file unit
!! \param filename filename of model file
!! \param errmsg error message
!! \return error message
!
  subroutine writeFileKernelInvertedModel(this,filename,lu,errmsg)
    type (kernel_inverted_model) :: this
    integer :: lu
    character(len=*) :: filename
    type (error_message) :: errmsg
    character(len=300) :: errstr
    character(len=28) :: myname = 'writeFileKernelInvertedModel'
    character(len=character_length_param) :: param_name
    integer :: ios,iparam
    logical :: file_exists
    character(len=7) :: open_status
    !
    call addTrace(errmsg,myname)
    !
    if(trim(this%parametrization) == '') then
       call add(errmsg,2,"object not initiated yet. please call initiateKernelInvertedModel first",myname)
       return
    endif
    !
    ! open file to write
    inquire(file=trim(filename),exist=file_exists)
    if(file_exists) then
       open_status = 'replace'
       write(errstr,*) "file '"//trim(filename)//"' already exists, replacing it now!"
       call add(errmsg,1,trim(errstr),myname)
    else
       open_status = 'new'
    endif
    open(unit=lu,file=trim(filename),form='unformatted',status=trim(open_status),action='write',access='stream',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "opening unformatted file '"//trim(filename)//"' to write raised iostat = ",ios
       call add(errmsg,2,trim(errstr),myname)
       close(lu)
       return
    endif
    call add(errmsg,0,"successfully opened file '"//trim(filename)//"' to write",myname)
    !
    ! write object to file
    write(lu) this%parametrization
    do while(nextParamModelParametrization(this%parametrization,param_name))
       iparam = indexOfParamModelParametrization(this%parametrization,param_name)
       if(associated(getVectorPointer(this%model_values(iparam))) .and. &
            associated(getVectorPointer(this%indx(iparam)))) then
          write(lu) param_name,size(getVectorPointer(this%indx(iparam)))
          write(lu) getVectorPointer(this%indx(iparam)),getVectorPointer(this%model_values(iparam))
       endif
    enddo ! while (nextParam)
    close(lu)
  end subroutine writeFileKernelInvertedModel
!------------------------------------------------------------------------
!> \brief write model object as vtk files
!! \details for each model parameter, write an own file using values and invgrid indices for that parameter
!! \param this kernel inverted model object
!! \param invgrid inversion grid on which the model values live (indices stored in this are assumed to relate to invgrid)
!! \param lu file unit
!! \param filebase basename for vtk output files (will be appended with parametrization and parameters)
!! \param errmsg error message
!! \return error message
!
  subroutine writeVtkKernelInvertedModel(this,invgrid,vtk_format,basename,lu,errmsg,overwrite)
    type (kernel_inverted_model) :: this
    type (inversion_grid) :: invgrid
    character(len=*) :: vtk_format,basename
    integer :: lu
    logical, optional :: overwrite
    type (error_message) :: errmsg
    type (invgrid_vtk_file) :: invgrid_vtk
    character(len=300) :: vtk_filename
    character(len=100) :: vtk_title,vtk_data_name
    character(len=27) :: myname = 'writeVtkKernelInvertedModel'
    integer :: iparam
    logical :: overwrite_file
    character(len=character_length_param) :: param_name
    integer, dimension(:), pointer :: indx
    real, dimension(:), pointer :: model_values
!
    nullify(indx,model_values)
!
    call addTrace(errmsg,myname)
    if(present(overwrite)) then
       overwrite_file = overwrite
    else
       overwrite_file = .false.
    end if
!
    if(trim(this%parametrization) == '') then
       call add(errmsg,2,"object not initiated yet. please call initiateKernelInvertedModel first",myname)
       return
    endif
!
    do while(nextParamModelParametrization(this%parametrization,param_name))
       iparam = indexOfParamModelParametrization(this%parametrization,param_name)
       model_values => getVectorPointer(this%model_values(iparam))
       indx => getVectorPointer(this%indx(iparam))
       ! create vtk file only, if there are any model values for this parameter
       if(associated(model_values) .and. associated(indx)) then
!
          ! initiate vtk file
          vtk_filename = trim(basename)//'_'//trim(this%parametrization)//'-'//trim(param_name)
          vtk_title = trim(this%parametrization)//'-'//trim(param_name)//' model values on inversion grid'
          call init(invgrid_vtk,invgrid,vtk_filename,vtk_format,errmsg,trim(vtk_title),cell_indx_req=indx)
          if(.level.errmsg == 2) return
!
          call add(errmsg,0,"creating vtk file with basename '"//trim(vtk_filename)//"'",myname)
!
          ! write model values to file
          vtk_data_name = trim(this%parametrization)//'-'//trim(param_name)
          call writeData(invgrid_vtk,lu,model_values,errmsg,data_name=trim(vtk_data_name),overwrite=overwrite_file)
          if(.level.errmsg == 2) return
!
          call dealloc(invgrid_vtk)
       endif
    enddo ! while (nextParam)
  end subroutine writeVtkKernelInvertedModel
!------------------------------------------------------------------------
!> \brief create kernel inverted model object as a copy of another object
  subroutine copyKernelInvertedModel(this,that)!,param,icell) ! in the future: introduce optional variables defining a subset of model parameters or a subset of inversion grid cells (or anything sensible)
    type (kernel_inverted_model) :: this,that
    ! local
    integer :: size_that,i
    real, dimension(:), pointer :: rpthis,rpthat
    integer, dimension(:), pointer :: ipthis,ipthat
!
    nullify(rpthis,rpthat,ipthis,ipthat)
!
    if(this%parametrization /= '') call deallocateKernelInvertedModel(this)
    if(that%parametrization == '') return
!
    ! parametrization
    this%parametrization = that%parametrization
!
    ! model_values
    if(associated(that%model_values)) then
       size_that = size(that%model_values)
       allocate(this%model_values(size_that))
       do i = 1,size_that
          rpthat => getVectorPointer(that%model_values(i))
          if(associated(rpthat)) then
             allocate(rpthis(size(rpthat)))
             rpthis = rpthat
             call associateVectorPointer(this%model_values(i),rpthis)
             nullify(rpthis)
          end if ! associated(rpthat)
       end do ! i
    end if ! associated(that%model_values)
!
    ! indx
    if(associated(that%indx)) then
       size_that = size(that%indx)
       allocate(this%indx(size_that))
       do i = 1,size_that
          ipthat => getVectorPointer(that%indx(i))
          if(associated(ipthat)) then
             allocate(ipthis(size(ipthat)))
             ipthis = ipthat
             call associateVectorPointer(this%indx(i),ipthis)
             nullify(ipthis)
          end if ! associated(ipthat)
       end do ! i
    end if ! associated(that%indx)
  end subroutine copyKernelInvertedModel
!------------------------------------------------------------------------
!> \brief pack a given vector of model values to a kernel_inverted_model object according to a given model space info object
  subroutine packVectorToKernelInvertedModel(this,model_vector,mspace,errmsg)
    type (kernel_inverted_model) :: this
    real, dimension(:) :: model_vector
    type (data_model_space_info) :: mspace
    type (error_message) :: errmsg
    ! local
    character(len=31) :: myname = 'packVectorToKernelInvertedModel'
    character(len=400) :: errstr
    character(len=character_length_param) :: param_name
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: idx_mspace,pcell
!
    nullify(pparam,idx_mspace,pcell)
!
    call addTrace(errmsg,myname)
!
    call init(this,.pmtrz.mspace,errmsg)
    if (.level.errmsg == 2) return
!
    if(size(model_vector) /= .nmval.mspace) then
       write(errstr,*) "size of incoming vector of model values = ",size(model_vector),&
            " differs from number of model parameters contained in incoming model space info object = ",&
            .nmval.mspace,", hence cannot interpret the model vector"
       call add(errmsg,2,errstr,myname)
       return       
    end if
!
    if(size(model_vector) == 0) then
       call add(errmsg,1,"incoming vector of model values is empty, hence created empty model object",myname)
       return
    end if
!
    nullify(pparam,pcell,idx_mspace)
    do while(nextParamModelParametrization(.pmtrz.mspace,param_name))
       if(associated(pparam)) deallocate(pparam); nullify(pparam)
       if(associated(pcell)) deallocate(pcell); nullify(pcell)
       if(associated(idx_mspace)) deallocate(idx_mspace)
       allocate(pparam(1)); pparam(1) = param_name
       idx_mspace => getIndxModelValues(mspace,param=pparam,cell=pcell)
       if(associated(idx_mspace)) then
          write(errstr,*) "there are ",size(idx_mspace)," '"//trim(param_name)//"' values in the model vector"
          call add(errmsg,0,errstr,myname)
          call addValuesKernelInvertedModel(this,param_name,pcell,model_vector(idx_mspace),errmsg)
          if (.level.errmsg == 2) goto 1
          deallocate(idx_mspace)
       else
          write(errstr,*) "there are no '"//trim(param_name)//"' values in the model vector"
          call add(errmsg,0,errstr,myname)
       end if
    end do ! while next param
!
    ! clean up
1   if(associated(pparam)) deallocate(pparam)
    if(associated(pcell)) deallocate(pcell)
    if(associated(idx_mspace)) deallocate(idx_mspace)
  end subroutine packVectorToKernelInvertedModel
!------------------------------------------------------------------------
!> \brief unpack a given kernel_inverted_model object to a vector of model values according to a given model space info object
  subroutine unpackToVectorKernelInvertedModel(model_vector,this,mspace,errmsg)
    type (kernel_inverted_model) :: this
    real, dimension(:), pointer :: model_vector
    type (data_model_space_info) :: mspace
    type (error_message) :: errmsg
    ! local
    character(len=33) :: myname = 'unpackToVectorKernelInvertedModel'
    character(len=400) :: errstr
    logical, dimension(:), allocatable :: initiated
    character(len=character_length_param) :: param_name
    integer, dimension(:), pointer :: pcell,idx_mspace,idx_kim
    real, dimension(:), pointer :: pval
!
    nullify(pcell,idx_mspace,idx_kim,pval)
!
    call addTrace(errmsg,myname)
    nullify(model_vector)
!
    if(.nmval.mspace == 0) then
       call add(errmsg,1,"no model parameters defined in the incoming model space, hence returning empty model vector",myname)
       return
    end if
    if(.pmtrz.mspace == '') then
       call add(errmsg,1,"incoming model space undefined, hence returning empty model vector",myname)
       return
    end if
    if(trim(this%parametrization) == '') then
       call add(errmsg,2,"this kernel_inverted_model object is not initiated yet and does not contain any values",myname)
       return
    end if
    if(trim(this%parametrization) /= trim(.pmtrz.mspace)) then
       write(errstr,*) "parametrization '"//trim(.pmtrz.mspace)//"' of incoming model space differs from parametrization '"&
            //trim(this%parametrization)//"' of this kernel_inverted_model object"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! allocate return vector and vector for checks for the requested number of values
    allocate(model_vector(.nmval.mspace),initiated(.nmval.mspace))
    initiated(:) = .false.
!
     do while (nextParamModelParametrization(this%parametrization,param_name))
        pcell => getIndicesKernelInvertedModel(this,param_name)
        pval => getValuesKernelInvertedModel(this,param_name)
!
        if(.not.(associated(pcell) .and. associated(pval))) cycle
!
        if(associated(idx_kim)) deallocate(idx_kim)
        if(associated(idx_mspace)) deallocate(idx_mspace)
        call mapCellIndicesDataModelSpaceInfo(mspace,param_name,pcell,idx_kim,idx_mspace)
!
        if(.not.(associated(idx_kim).and.associated(idx_mspace))) cycle
!
        if(any(initiated(idx_mspace))) then
           write(errstr,*) count(initiated(idx_mspace))," values already initiated when trying to set '"//trim(param_name)//&
                "' values: THIS ERROR SHOULD NOT HAPPEN!"
           call add(errmsg,2,errstr,myname)
           goto 2
        end if
!
        ! idx_mspace should contain index values from 1 to .nmval.mspace, i.e. the next assignment should work 
        ! (otherwise module dataModelSpaceInfo, or its routine mapCellIndicesDataModelSpaceInfo is inconsistent)
        model_vector(idx_mspace) = pval(idx_kim)
!
        ! remember that these indices have been set
        initiated(idx_mspace) = .true.
     end do ! while nextParam
     if(associated(idx_kim)) deallocate(idx_kim)
     if(associated(idx_mspace)) deallocate(idx_mspace)
!
     if(any(.not.initiated)) then
        write(errstr,*) count(.not.initiated)," of the requested values in the model vector could not be defined, this "//&
             "means that there are model values in the model space which are not contained in this kernel_inverted_model_object"
        call add(errmsg,2,errstr,myname)
        goto 2
     end if
!
    ! clean up and return
1   if(allocated(initiated)) deallocate(initiated)
    return
!
    ! If the code comes here, there was an error after allocating the return variable, hence deallocate it before cleaning up and returning
2   deallocate(model_vector); nullify(model_vector)
    if(associated(idx_kim)) deallocate(idx_kim)
    if(associated(idx_mspace)) deallocate(idx_mspace)
    goto 1
  end subroutine unpackToVectorKernelInvertedModel
!	
end module kernelInvertedModel
