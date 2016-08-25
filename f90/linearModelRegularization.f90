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
!> \brief define and handle linear equations, which as smoothing or damping or similar conditions 
!!  are added to a kernel_linear_system object
!!
!! \details as the regularization equations will generally be very sparse, avoid storing zeros 
!!  by using vector_pointer structures and only store relevant coefficients of an equation
!!  (along with die indices of the corresponding variables of that equation)
!!
!! \author Florian Schumacher
!! \date Okt 2015
!
module linearModelRegularization
!
  use modelParametrization
  use dataModelSpaceInfo
  use inversionGrid
  use vectorPointer
  use kernelLinearSystem
  use errorMessage
!
  implicit none
!
  interface init; module procedure initiateLinearModelRegularization; end interface
  interface dealloc; module procedure deallocateLinearModelRegularization; end interface
  interface addSmoothing; module procedure addNeighbourAverageSmoothingLinearModelRegularization; end interface
  interface addDamping; module procedure addDampingLinearModelRegularization; end interface
  interface getEquations
     module procedure getEquationsMatricesLinearModelRegularization
     module procedure getEquationsVectorPointersLinearModelRegularization
  end interface getEquations
  interface operator (.neq.); module procedure getNeqtotalLinearModelRegularization ; end interface
  interface operator (.requires.); module procedure requiresQuantityLinearModelRegularization ; end interface
!
  integer, parameter :: character_length_regscal_type = 31
!
!> \brief type contains indices and values which define equations in a linear system 
  type linear_model_regularization
     private
     logical :: initiated = .false. !< flag to indicate whether initiateLinearModelRegularization was successfully called
     type (data_model_space_info) :: mspace !< model space, the model values of which are used to define equation arrays
     integer :: nparam !< number of different parameter names in mspace: nparam = size(param)
     character(len=character_length_param), dimension(:), pointer :: param => null() !< all different parameter names contained in mspace
     ! scaling
     character(len=character_length_regscal_type) :: scaling_type = '' !< type of scaling supported by this model
     ! equations
     integer :: neq_total = 0 !< total number of equations, neq_total == sum(neq)
     integer, dimension(:), pointer :: neq => null() !< number of equations per model parameter
     type (integer_vector_pointer), dimension(:,:), pointer :: idx => null() !< indices of variables with non-zero coefficient (per model parameter)
     type (real_vector_pointer), dimension(:,:), pointer :: coef => null() !< coefficients corresponding to indices in idx (per model parameter)
     real, dimension(:,:), pointer :: rhs => null() !< rhs values of equations (per model parameter)
  end type linear_model_regularization
!
contains
!
!------------------------------------------------------------------------
!> \brief add simple Laplace smoothing (average over neighbours) to possibly existing regularization equations
!! \param this linear model smoothing conditions
!! \param dmspace data model space info object
!! \param errmsg error message
!! \param scaling_type optional string defining a method to scale the regularization equations dependent on values of the kernel matrix
  subroutine initiateLinearModelRegularization(this,dmspace,errmsg,scaling_type)
    type (linear_model_regularization) :: this
    type (data_model_space_info) :: dmspace
    type (error_message) :: errmsg
    character(len=*), optional :: scaling_type
    ! local
    character(len=33) :: myname = 'initiateLinearModelRegularization'
    character(len=character_length_param), dimension(:), pointer :: param_dmspace
    character(len=character_length_param) :: param_name
    integer :: iparam_this
!
    nullify(param_dmspace)
!
    call addTrace(errmsg,myname)
!
    if(this%initiated) then
       call add(errmsg,2,"Linear model regularization equations are already defined (initiated). "//&
            "Please deallocate object before creating a new one",myname)
       return
    end if
!
    if(.nmval.dmspace == 0) then
       call add(errmsg,2,"there are no parameters in model space",myname)
       return
    end if
!
    if(present(scaling_type)) then
       select case(scaling_type)
       case('none','') ! no scaling
          this%scaling_type = ''
       case('absmax_per_param,overall_factor')
          this%scaling_type = 'absmax_per_param,overall_factor'
       case('absmax_per_param,param_factors')
          this%scaling_type = 'absmax_per_param,param_factors'
       !case('your_scaling_type')
       !   process your type here and define any expected values this%scaling_values below in routines
       !   addSmoothing, addDamping
       case default
          call add(errmsg,2,"scaling_type '"//trim(scaling_type)//"' not supported. Supported at the moment are: "//&
               "'none', 'absmax_per_param,overall_factor', 'absmax_per_param,param_factors'",myname)
          return
       end select ! case scaling_type
       call add(errmsg,0,"incoming scaling type: '"//trim(scaling_type)//"'",myname)
    else
       call add(errmsg,0,"NO incoming scaling type",myname)
    end if
!
    param_dmspace => null()
    param_dmspace => allParam(dmspace)
    if(.not.associated(param_dmspace)) then
       call add(errmsg,2,"THIS ERROR SHOULD NOT OCCUR! no parameter names returned by model space object",myname)
       goto 2
    end if
!
    ! now define vector this%param containing exactly the same values as param_dmspace, but with the values
    ! in conventional order as defined by module modelParametrization (with possible gaps caused by missing 
    ! parameters). dataModelSpaceInfo module returns values param_dmspace possibly in mixed order.
    ! Find the correct permutation of array param_dmspace now, using the parameter name iterator
    this%nparam = size(param_dmspace)
    allocate(this%param(this%nparam))
    iparam_this = 0
    do while(nextParamModelParametrization(.pmtrz.dmspace,param_name))
       if(any(param_dmspace == param_name)) then
          iparam_this = iparam_this + 1
          this%param(iparam_this) = param_name
       end if
    end do ! while nextParam
    if(iparam_this < this%nparam) then
       call add(errmsg,2,"There were duplicate parameter names returned as 'allParam' by data model space info "//&
            "module. The module dataModelSpaceInfo is corrupt.",myname)
       goto 2
    end if
!
    ! only memorize the model space, i.e. do not copy the data space (by setting first and last data sample 
    ! index to 0 and allow for all model value indices)
    call copyDataModelSpaceInfo(this%mspace,dmspace,idata1=0,idata2=0)
!
    allocate(this%neq(this%nparam))
!
    this%initiated = .true.
!
1   if(associated(param_dmspace)) deallocate(param_dmspace)
    return
!
2   call deallocateLinearModelRegularization(this)
    goto 1
  end subroutine initiateLinearModelRegularization
!------------------------------------------------------------------------
!> \brief add simple Laplace smoothing (average over neighbours) to possibly existing regularization equations
!! \details boundary_conditions can have values: 
!!  '' : empty string is equivalent to .not.present(boundary_conditions)
!!  'zero_all_outer_bnd' : on all outer boundaries, apply zero boundary smoothing conditions
!!  'zero_burried_outer_bnd,cont_free_surface' : on burried outer boundaries, apply zero boundary smoothing 
!!        conditions, on free surfaces apply continuous smoothing boundary conditions
!! \param this linear model smoothing conditions
!! \param invgrid inversion grid used to get neighbours of cells
!! \param errmsg error message
!! \param scaling_values optional vector of numbers used for scaling, corresponding to the type of scaling (type specific)
!! \param boundary_conditions optional character string indicating whether to define different smoothing conditions for free surfaces and/or burried invgrid boundaries
!
  subroutine addNeighbourAverageSmoothingLinearModelRegularization(this,invgrid,errmsg,scaling_values,&
       boundary_conditions,neq_added)
    type (linear_model_regularization) :: this
    type (inversion_grid) :: invgrid
    type (error_message) :: errmsg
    real, dimension(:), optional :: scaling_values
    character(len=*), optional :: boundary_conditions
    integer, optional :: neq_added
    ! locally allocated smoothing equations
    integer :: neq_total_smooth
    integer, dimension(:), allocatable :: neq_smooth
    type (integer_vector_pointer), dimension(:,:), allocatable :: idx_smooth
    type (real_vector_pointer), dimension(:,:), allocatable :: coef_smooth
    real, dimension(:,:), allocatable :: rhs_smooth
    type (integer_vector_pointer), dimension(:,:), pointer :: tmp_idx
    type (real_vector_pointer), dimension(:,:), pointer :: tmp_coef
    real, dimension(:,:), pointer :: tmp_rhs
    ! local other
    character(len=53) :: myname = 'addNeighbourAverageSmoothingLinearModelRegularization'
    character(len=400) :: errstr
    real, dimension(:), allocatable :: scaling_values_local
    type (integer_vector_pointer), dimension(:), pointer :: idx_invgrid_nb
    integer :: nidx_invgrid_nb,imval,iparam,iparam_pmtrz,nparam_pmtrz,nmval_mspace,&
         ieq,ieq_shift,nidx_param_cell,n,n_artif,max_icell_param,max_neq
    real :: scaling_value
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: idx_param_cell,pcell,nb_idx_p,eq_idx
    integer, dimension(:), allocatable :: pcell_inverse,nb_idx_valid
    real, dimension(:), pointer :: eq_coef
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name
    character(len=100) :: smoothing_boundary_conditions,invgrid_boundary_conditions
!
    nullify(tmp_idx,tmp_coef,tmp_rhs,idx_invgrid_nb,pparam,idx_param_cell,pcell,nb_idx_p,eq_idx,eq_coef)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"Regularization object not yet initiated. Initiate before adding any regularization "//&
            "conditions",myname)
       return
    end if
!
    if(this%neq_total /= 0) then
       call add(errmsg,1,"Some linear model regularization condititions were already defined, possibly damping "//&
            "conditions. If you did not request damping, there could be some problem. Adding smoothing "//&
            "conditions to the existing equations now.",myname)
    end if
!
    parametrization = .pmtrz.(this%mspace)
    nparam_pmtrz = numberOfParamModelParametrization(parametrization)
!
    select case(this%scaling_type)
    case('')
       ! in this case, do not expect any scaling values. If there are any scaling values given nonetheless, 
       ! raise a warning that those are ignored and that there is no scaling now
       if(present(scaling_values)) then
          call add(errmsg,1,"Regularization object was initiated WITHOUT any scaling. However, there are incoming "//&
               "scaling values. Those will be ignored, there will be no scaling!",myname)
       end if
!
    case('absmax_per_param,overall_factor')
       if(.not.present(scaling_values)) then
          call add(errmsg,2,"Regularization object was initiated for scaling of type 'absmax_per_param,overall_factor', "//&
               "hence expect one scaling value here. However, there are no incoming scaling values.",myname)
          return
       end if
       if(size(scaling_values)>1) then
          call add(errmsg,1,"Regularization object was initiated for scaling of type 'absmax_per_param,overall_factor', "//&
               "hence expect one scaling value here. Since there are more than one value given, will use only the first "//&
               "one here ignoring the others!",myname)
       end if
       if(size(scaling_values)>=1) then
          if(scaling_values(1)<= 0.) then
             call add(errmsg,2,"object initiated for scaling type 'absmax_per_param,overall_factor'. "//&
                  "The first incoming scaling value, however, is <= 0.0 . This is not supported.",myname)
             return
          end if
          ! The scaling values do not necessarily need to be remembered by this object.
          ! They are incorporated in the smoothing equations and are (for now) not requried to be known afterwards.
          ! In the future (if necessary) the incoming scaling values could also be added to the derived type 
          ! linear_model_regularization. Note that there will be different such scaling values for each type of 
          ! regularization (smoothing, damping etc., possibly also depth dependent values)
          allocate(scaling_values_local(1)) 
          scaling_values_local(1:1) = scaling_values(1:1)
       else
          call add(errmsg,2,"in case of scaling_type = 'absmax_per_param,overall_factor', "//&
               "scaling_values must at least contain one value (the first is used)",myname)
          return                   
       end if
!
    case('absmax_per_param,param_factors')
       if(.not.present(scaling_values)) then
          call add(errmsg,2,"Regularization object was initiated for scaling of type 'absmax_per_param,param_factors', "//&
               "hence expect one scaling value pere model parameter here. However, there are no incoming scaling values.",&
               myname)
          return
       end if
       if(size(scaling_values)>nparam_pmtrz) then
          write(errstr,*) "Regularization object was initiated for scaling of type 'absmax_per_param,param_factors', "//&
               "hence expect ",nparam_pmtrz," values, one for each model parameter of current parametrizaton '",&
               trim(parametrization)//"' (in conventional order). Since there are more values given here, will use "//&
               "the first ones ignoring the others!"
          call add(errmsg,1,errstr,myname)
       end if
       if(size(scaling_values)>=nparam_pmtrz) then
          ! The scaling values do not necessarily need to be remembered by this object.
          ! They are incorporated in the smoothing equations and are (for now) not requried to be known afterwards.
          ! In the future (if necessary) the incoming scaling values could also be added to the derived type 
          ! linear_model_regularization. Note that there will be different such scaling values for each type of 
          ! regularization (smoothing, damping etc., possibly also depth dependent values)
          allocate(scaling_values_local(this%nparam))
          scaling_values_local = 0.0
          do while(nextParamModelParametrization(parametrization,param_name,iparam_pmtrz))
             if(any(this%param==param_name)) then
                where(this%param==param_name)
                   scaling_values_local = scaling_values(iparam_pmtrz)
                end where
             end if
          end do ! while(nextParam)
          if(any(scaling_values_local <= 0)) then
             call add(errmsg,2,"object initiated for scaling type 'absmax_per_param,param_factors'. However, "//&
                  "after processing there remain scaling values being <= 0.0 . This is not supported. Make sure "//&
                  "you know which values need to be provided.",myname)
             return
          end if
       else
          write(errstr,*) "in case of scaling_type = 'absmax_per_param,param_factors', scaling_values must contain at least "//&
               "as much values as there are model parameters (in conventional order) for the current model parametrization '",&
               trim(parametrization),"' (",nparam_pmtrz,"). values given for unused parameters will be ignored (just fill the "//&
               "vector by zeros or somehting)."
          call add(errmsg,2,errstr,myname)
          return
       end if
    end select ! this%scaling_type
!
    ! Variable invgrid_boundary_conditions is the only mechanism here to control the smoothing boundary condition below.
    ! It is assumed, that artificial neighbours beyond bounadries are indicated by neighbour index -1 and that
    ! in no other situation neighbour indices are negative.
    ! For zero smoothing boundary conditions, simply all neighbour indices == -1 are accounted for as artificial 
    ! neighbours and their number is added to the number of actual valid neighbours to define the smoothing weights. 
    ! Continuity boundary conditions are automatically accounted for on outer boundaries (and in the future possibly
    ! on innder boundaries), since in this case simply some neighbours are left out and are missing here. 
    ! Variable smoothing_boundary_conditions is unused below (for now...)
    smoothing_boundary_conditions  = 'continuous'
    invgrid_boundary_conditions  = 'standard'
    if(present(boundary_conditions)) then
       select case(boundary_conditions)
       case('','standard','continuous')
          ! This is defined to be equivalent to .not.present(boundary_conditions)
          ! So do nothing, keep values 'continuous' and 'standard'
       case('zero_all_outer_bnd')
          smoothing_boundary_conditions = boundary_conditions
          invgrid_boundary_conditions = 'extra_nbs_outer_bnd'
       case('zero_burried_outer_bnd,cont_free_surface')
          smoothing_boundary_conditions = boundary_conditions
          invgrid_boundary_conditions = 'extra_nbs_outer_bnd_except_free_surface'
       case default
          call add(errmsg,2,"incoming smoothing boundary condition flag '"//trim(boundary_conditions)//&
               "' not supported. Supported flags are: 'zero_all_outer_bnd', 'zero_burried_outer_bnd,cont_free_surface'",myname)
          return
       end select
       call add(errmsg,0,"incoming smoothing boundary condition flag: '"//trim(boundary_conditions)//&
            "'; using neighbour boundary condition '"//trim(invgrid_boundary_conditions)//"' for inversion grid",myname)
    else
       call add(errmsg,0,"NO incoming smoothing boundary condition flag, will apply standard boundary conditions",myname)
    endif
!
    call getIndicesFaceNeighboursInversionGrid(invgrid,idx_invgrid_nb,boundary_conditions=invgrid_boundary_conditions)
    if(.not.associated(idx_invgrid_nb)) then
       ! Only raise warning, do not create any smoothing. However, it is not very likely that
       ! there are actually no neighbours, so it might be better to raise an error here (?!)
       call add(errmsg,1,"Inversion grid did not return any neighbours. This may suggest that the inversion grid was "//&
            "not properly created, or that the neighbour boundary conditions '"//trim(invgrid_boundary_conditions)//&
            "' are not supported for this inversion grid",myname)
       return
    end if
    nidx_invgrid_nb = size(idx_invgrid_nb)
!
    ! allocate for maximal number of smoothing equations and known number of parameters
    nmval_mspace = .nmval.(this%mspace)
    allocate(neq_smooth(this%nparam),idx_smooth(nmval_mspace,this%nparam),coef_smooth(nmval_mspace,this%nparam),&
         rhs_smooth(nmval_mspace,this%nparam))
!
    neq_smooth(:) = 0
    do iparam = 1,this%nparam
!
       ! define here the correct scaling value to be used in loop below, 
       select case(this%scaling_type)
       case(''); scaling_value = 1.0 ! no scaling, hence set to 1.0
       case('absmax_per_param,overall_factor'); scaling_value = scaling_values_local(1)
       case('absmax_per_param,param_factors'); scaling_value = scaling_values_local(iparam)
       end select
!
       ! count the equations added for this parameter
       ieq = 0
!
       ! for this parameter name (e.g. vp) get all model parameter indices and invgrid cell indices
       if(associated(pparam)) deallocate(pparam)
       if(associated(pcell)) deallocate(pcell); nullify(pcell)
       if(associated(idx_param_cell)) deallocate(idx_param_cell)
       allocate(pparam(1)); pparam(1) = this%param(iparam)
       idx_param_cell => getIndxModelValues(this%mspace,param=pparam,cell=pcell)
       if(.not.associated(idx_param_cell)) then
          write(errstr,*) "There are no model values in the model for parameter '",trim(this%param(iparam)),&
               "' although there should be by definition. This error should not occur!!"
          call add(errmsg,2,errstr,myname)
          return
       end if
       nidx_param_cell = size(idx_param_cell)
!
       ! define the inverse mapping pcell_inverse of array pcell which maps an inversion grid cell index to its position in array pcell
       ! i.e.  pcell_inverse(i) = j <=> pcell(j) = i
       if(allocated(pcell_inverse)) deallocate(pcell_inverse)
       max_icell_param = maxval(pcell)
       allocate(pcell_inverse(max_icell_param)); pcell_inverse = -1
       do n = 1,nidx_param_cell
          pcell_inverse(pcell(n)) = n
       end do ! i

       ! now loop on all model values of current parameter name this%param(iparam)
       do imval = 1,nidx_param_cell

          ! get invgrid cell neighbours of current model value
          nb_idx_p => getVectorPointer(idx_invgrid_nb(pcell(imval)))
          if(.not.associated(nb_idx_p)) cycle
!
          ! Memorize how many artificial neighbours this cell has.
          ! Variable invgrid_boundary_conditions is the only mechanism here to control the smoothing boundary condition.
          ! It is assumed, that artificial neighbours beyond bounadries are indicated by neighbour index -1 and that
          ! in no other situation neighbour indices are negative.
          ! For zero smoothing boundary conditions, simply all neighbour indices == -1 returned by the inversion grid 
          ! module (passing boundary flag invgrid_boundary_conditions) are accounted for as artificial neighbours and their 
          ! number is added to the number of actual valid neighbours to define the smoothing weights. 
          ! Continuity boundary conditions are automatically accounted for on outer boundaries (and in the future possibly
          ! on innder boundaries), since in this case simply some neighbours are left out and are missing here 
          ! so that n_artif should be 0 in this case (having no effect when defining the coefficients below).
          n_artif = count(nb_idx_p == -1)
!
          ! only use inversion grid cell neighbours which are also in model space for current parameter name

          ! first of all, throw away artificial neighbours and possible inversion grid indices which are larger than 
          ! max_icell_param. If there are no such cells (e.g. only artificial), we cannot set up a smoothing equation
          n = count(nb_idx_p <= max_icell_param .and. nb_idx_p > 0)
          if(n==0) cycle
          if(allocated(nb_idx_valid)) deallocate(nb_idx_valid)
          allocate(nb_idx_valid(n))
          nb_idx_valid = pack(nb_idx_p , nb_idx_p <= max_icell_param  .and. nb_idx_p > 0)
!
          n = count(pcell_inverse( nb_idx_valid ) > 0)
          ! If there are no valid neighbours for this parameter, we cannot set up a smoothing equation, even if there are artificial cells.
          if(n==0) cycle
!
          ! finally add new equation
          ieq = ieq+1
!
          ! set indices of variables in equation
          allocate(eq_idx(n+1))
          eq_idx(1:n) = idx_param_cell(pack(pcell_inverse(nb_idx_valid),pcell_inverse(nb_idx_valid) > 0)) ! mspace index of valid neighbours
          eq_idx(n+1) = idx_param_cell(imval) ! mspace index of center cell
          call associateVectorPointer(idx_smooth(ieq,iparam),eq_idx); nullify(eq_idx)
!
          ! set equation coefficients
          allocate(eq_coef(n+1))
          eq_coef(1:n) = (1./real(n+n_artif))*scaling_value ! coefficient of neighbour cells, account for zero smoothing boundary conditions by '+n_artif'
          eq_coef(n+1) = -1.*scaling_value ! smoothing coefficient of center cell
          call associateVectorPointer(coef_smooth(ieq,iparam),eq_coef); nullify(eq_coef)
!
          ! set rhs of equation
          rhs_smooth(ieq,iparam) = 0.!*scaling_value
!
       end do ! imval
!
       neq_smooth(iparam) = ieq
!
    end do ! iparam
!
    neq_total_smooth = sum(neq_smooth)
!
    if(neq_total_smooth == 0) then
       call add(errmsg,1,"there were no smoothing equations added to the regularization object",myname)
       goto 1
    end if

    ! (re)allocate for the additional regularization equations (smoothing equations) found above
    if(this%neq_total == 0) then
       allocate(this%neq(this%nparam))
       this%neq = 0
       max_neq = maxval(neq_smooth)
       allocate(this%neq(this%nparam),this%idx(max_neq,this%nparam),this%coef(max_neq,this%nparam),&
            this%rhs(max_neq,this%nparam))
    else
       ! Need to re-allocate first dimension of equation arrays to the maximum size per parameter AFTER
       ! adding the smoothing equations.
       ! Since it might happen that smoothing equations were only found for parameters which do not have
       ! (a lot of) equations yet, we need a further check whether we need to reallocate at all
       max_neq = maxval(this%neq+neq_smooth)
       if(max_neq > size(this%idx,1)) then
          ! In this case we need to reallocate! (assume that the first dimensions of this%idx,this%coef,this%rhs
          ! all have the same size!)
          ! First, allocate temporary arrays of new size
          allocate(tmp_idx(max_neq,this%nparam),tmp_coef(max_neq,this%nparam),tmp_rhs(max_neq,this%nparam))
          ! Second, copy all contents from this to temporary arrays (or link pointers, respectively)
          do iparam = 1,this%nparam
             do ieq = 1,this%neq(iparam)
                eq_idx => getVectorPointer(this%idx(ieq,iparam))
                call associateVectorPointer(tmp_idx(ieq,iparam),eq_idx)
                eq_coef => getVectorPointer(this%coef(ieq,iparam))
                call associateVectorPointer(tmp_coef(ieq,iparam),eq_coef)
                tmp_rhs(ieq,iparam) = this%rhs(ieq,iparam)
             end do ! ieq
          end do ! iparam
          ! Third, deallocate arrays from this (attention: do not call dealloc on the vector pointers, since they
          ! were linked to the temporary arrays, just deallocate this%idx, this%coef)
          deallocate(this%idx,this%coef,this%rhs)
          ! Fourth, let the temporary arrays point to this
          this%idx => tmp_idx
          this%coef => tmp_coef
          this%rhs => tmp_rhs
          ! Fifth, let temporary arrays point to null() since they must not be deallocated below before return
          nullify(tmp_idx,tmp_coef,tmp_rhs)
       end if
    end if

    ! Finally add the above found smoothig equations to the existing regularization equations
    do iparam = 1,this%nparam
       ieq_shift = this%neq(iparam)
       do ieq = 1,neq_smooth(iparam)
          eq_idx => getVectorPointer(idx_smooth(ieq,iparam))
          call associateVectorPointer(this%idx(ieq_shift+ieq,iparam),eq_idx)
          eq_coef => getVectorPointer(coef_smooth(ieq,iparam))
          call associateVectorPointer(this%coef(ieq_shift+ieq,iparam),eq_coef)
          this%rhs(ieq_shift+ieq,iparam) = rhs_smooth(ieq,iparam)
       end do ! ieq
    end do ! iparam
    ! deallocate local smoothing equation arrays (attention: do not call dealloc on the vector pointers, since they
    ! were linked to the arrays in this, just deallocate idx_smooth,coef_smooth,rhs_smooth)
    deallocate(idx_smooth,coef_smooth,rhs_smooth)

    ! Update the number of valid equations per parameter and the total number of equations
    this%neq = this%neq + neq_smooth
    this%neq_total = this%neq_total + neq_total_smooth
!
    if(present(neq_added)) neq_added = neq_total_smooth
!
    ! clean up
1   if(allocated(neq_smooth)) deallocate(neq_smooth)
    if(allocated(idx_smooth)) then
       do iparam = 1,size(idx_smooth,2)
          do ieq = 1,size(idx_smooth,1)
             call dealloc(idx_smooth(ieq,iparam))
          end do ! ieq
       end do ! iparam
       deallocate(idx_smooth)
    end if
    if(allocated(coef_smooth)) then
       do iparam = 1,size(coef_smooth,2)
          do ieq = 1,size(coef_smooth,1)
             call dealloc(coef_smooth(ieq,iparam))
          end do ! ieq
       end do ! iparam
       deallocate(coef_smooth)
    end if
    if(allocated(rhs_smooth)) deallocate(rhs_smooth)
    if(associated(idx_invgrid_nb)) then
       do n = 1,size(idx_invgrid_nb)
          call dealloc(idx_invgrid_nb(n))
       end do
       deallocate(idx_invgrid_nb)
    end if
    if(associated(pparam)) deallocate(pparam)
    if(associated(pcell)) deallocate(pcell)
    if(associated(idx_param_cell)) deallocate(idx_param_cell)
    if(allocated(pcell_inverse)) deallocate(pcell_inverse)
    if(allocated(nb_idx_valid)) deallocate(nb_idx_valid)
    if(allocated(scaling_values_local)) deallocate(scaling_values_local)
  end subroutine addNeighbourAverageSmoothingLinearModelRegularization
!------------------------------------------------------------------------
!> \brief add damping equations to possibly existing regularization equations
!! \param this linear model smoothing conditions
!! \param errmsg error message
!! \param scaling_values optional vector of numbers used for scaling, corresponding to the type of scaling (type specific)
!
  subroutine addDampingLinearModelRegularization(this,errmsg,scaling_values,neq_added)
    type (linear_model_regularization) :: this
    type (error_message) :: errmsg
    real, dimension(:), optional :: scaling_values
    integer, optional :: neq_added
    ! locally allocated damping equations
    integer :: neq_total_damp
    integer, dimension(:), allocatable :: neq_damp
    type (integer_vector_pointer), dimension(:,:), allocatable :: idx_damp
    type (real_vector_pointer), dimension(:,:), allocatable :: coef_damp
    real, dimension(:,:), allocatable :: rhs_damp
    type (integer_vector_pointer), dimension(:,:), pointer :: tmp_idx
    type (real_vector_pointer), dimension(:,:), pointer :: tmp_coef
    real, dimension(:,:), pointer :: tmp_rhs
    ! local other
    character(len=35) :: myname = 'addDampingLinearModelRegularization'
    character(len=400) :: errstr
    real, dimension(:), allocatable :: scaling_values_local
    real :: scaling_value
    integer :: nparam_pmtrz,iparam_pmtrz,iparam,imval,ieq_shift,ieq,max_neq,nmval_mspace
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: midx_param,eq_idx
    real, dimension(:), pointer :: eq_coef
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name
!
    nullify(tmp_idx,tmp_coef,tmp_rhs,pparam,midx_param,eq_idx,eq_coef)
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"Regularization object not yet initiated. Initiate before adding any regularization "//&
            "conditions",myname)
       return
    end if
!
    if(this%neq_total /= 0) then
       call add(errmsg,1,"Some linear model regularization condititions were already defined, possibly smoothing "//&
            "conditions. If you did not request smoothing, there could be some problem. Adding damping "//&
            "conditions to the existing equations now.",myname)
    end if
!
    parametrization = .pmtrz.(this%mspace)
    nparam_pmtrz = numberOfParamModelParametrization(parametrization)
!
    select case(this%scaling_type)
    case('') 
       ! in this case, do not expect any scaling values. If there are any scaling values given nonetheless, 
       ! raise a warning that those are ignored and that there is no scaling now
       if(present(scaling_values)) then
          call add(errmsg,1,"Regularization object was initiated WITH NO scaling. However, there are incoming "//&
               "scaling values. Those will be ignored, there will be no scaling!",myname)
       end if
!
    case('absmax_per_param,overall_factor')
       if(.not.present(scaling_values)) then
          call add(errmsg,2,"Regularization object was initiated for scaling of type 'absmax_per_param,overall_factor', "//&
               "hence expect one scaling value here. However, there are no incoming scaling values.",myname)
          return
       end if
       if(size(scaling_values)>1) then
          call add(errmsg,1,"Regularization object was initiated for scaling of type 'absmax_per_param,overall_factor', "//&
               "hence expect one scaling value here. Since there is more than one value given, will use only the first "//&
               "one here ignoring the others!",myname)
       end if
       if(size(scaling_values)>=1) then
          if(scaling_values(1)<= 0.) then
             call add(errmsg,2,"object initiated for scaling type 'absmax_per_param,overall_factor'. "//&
                  "The first incoming scaling value, however, is <= 0.0 . This is not supported.",myname)
             return
          end if
          ! The scaling values do not necessarily need to be remembered by this object.
          ! They are incorporated in the damping equations and are (for now) not requried to be known afterwards.
          ! In the future (if necessary) the incoming scaling values could also be added to the derived type 
          ! linear_model_regularization. Note that there will be different such scaling values for each type of 
          ! regularization (smoothing, damping etc., possibly also depth dependent values)
          allocate(scaling_values_local(1)) 
          scaling_values_local(1:1) = scaling_values(1:1)
       else
          call add(errmsg,2,"in case of scaling_type = 'absmax_per_param,overall_factor', "//&
               "scaling_values must at least contain one value (the first is used)",myname)
          return                   
       end if
!
    case('absmax_per_param,param_factors')
       if(.not.present(scaling_values)) then
          call add(errmsg,2,"Regularization object was initiated for scaling of type 'absmax_per_param,param_factors', "//&
               "hence expect one scaling value pere model parameter here. However, there are no incoming scaling values.",&
               myname)
          return
       end if
       if(size(scaling_values)>nparam_pmtrz) then
          write(errstr,*) "Regularization object was initiated for scaling of type 'absmax_per_param,param_factors', "//&
               "hence expect ",nparam_pmtrz," values, one for each model parameter of current parametrizaton '",&
               trim(parametrization)//"' (in conventional order). Since there are more values given here, will use "//&
               "the first ones ignoring the others!"
          call add(errmsg,1,errstr,myname)
       end if
       if(size(scaling_values)>=nparam_pmtrz) then
          ! The scaling values do not necessarily need to be remembered by this object.
          ! They are incorporated in the damping equations and are (for now) not requried to be known afterwards.
          ! In the future (if necessary) the incoming scaling values could also be added to the derived type 
          ! linear_model_regularization. Note that there will be different such scaling values for each type of 
          ! regularization (smoothing, damping etc., possibly also depth dependent values)
          allocate(scaling_values_local(this%nparam))
          scaling_values_local = 0.0
          do while(nextParamModelParametrization(parametrization,param_name,iparam_pmtrz))
             if(any(this%param==param_name)) then
                where(this%param==param_name)
                   scaling_values_local = scaling_values(iparam_pmtrz)
                end where
             end if
          end do ! while(nextParam)
          if(any(scaling_values_local <= 0)) then
             call add(errmsg,2,"object initiated for scaling type 'absmax_per_param,param_factors'. However, "//&
                  "after processing there remain scaling values being <= 0.0 . This is not supported. Make sure "//&
                  "you know which values need to be provided.",myname)
             goto 1
          end if
       else
          write(errstr,*) "in case of scaling_type = 'absmax_per_param,param_factors', scaling_values must contain at least "//&
               "as much values as there are model parameters (in conventional order) for the current model parametrization '",&
               trim(parametrization),"' (",nparam_pmtrz,"). values given for unused parameters will be ignored (just fill the "//&
               "vector by zeros or somehting)."
          call add(errmsg,2,errstr,myname)
          return
       end if
    end select ! this%scaling_type
!
    ! allocate for maximal number of damping equations and known number of parameters
    nmval_mspace = .nmval.(this%mspace)
    allocate(neq_damp(this%nparam),idx_damp(nmval_mspace,this%nparam),coef_damp(nmval_mspace,this%nparam),&
         rhs_damp(nmval_mspace,this%nparam))
!
    neq_damp(:) = 0
    do iparam = 1,this%nparam
       ! define here the correct scaling value to be used in loop below, 
       select case(this%scaling_type)
       case(''); scaling_value = 1.0 ! no scaling, hence set to 1.0
       case('absmax_per_param,overall_factor'); scaling_value = scaling_values_local(1)
       case('absmax_per_param,param_factors'); scaling_value = scaling_values_local(iparam)
       end select
!
       ! for this parameter name (e.g. vp) get all model parameter indices
       if(associated(pparam)) deallocate(pparam)
       if(associated(midx_param)) deallocate(midx_param)
       allocate(pparam(1)); pparam(1) = this%param(iparam)
       midx_param => getIndxModelValues(this%mspace,param=pparam)
       if(.not.associated(midx_param)) then
          write(errstr,*) "There are no model values in the model for parameter '",trim(this%param(iparam)),&
               "' although there should be by definition. This error should not occur!!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
!
       ! now loop on all model values of current parameter name this%param(iparam)
       do imval = 1,size(midx_param)

          ! Add a damping equation for each model value in the model space associated with this parameter name.
          ! Regardless of the relation of this model value to the other model values, there is a damping equation
          ! (not like with smoothing equations which are only sensible if there are neighbouring cells of the 
          ! same parameter contained in the model space as well)

          ! set indices of variables in equation
          allocate(eq_idx(1)) ! equation has just one entry different from zero
          eq_idx(1) = midx_param(imval) ! model space index of this model value
          call associateVectorPointer(idx_damp(imval,iparam),eq_idx); nullify(eq_idx)

          ! set equation coefficients, add the damping coefficient on the diagonal of the damping matrix
          allocate(eq_coef(1))
          eq_coef(1) = 1.0*scaling_value
          call associateVectorPointer(coef_damp(imval,iparam),eq_coef); nullify(eq_coef)

          ! set rhs of equation
          rhs_damp(imval,iparam) = 0.
!
       end do ! imval
!
       neq_damp(iparam) = size(midx_param)
!
    end do ! iparam
!
    neq_total_damp = sum(neq_damp)
!
    if(neq_total_damp == 0) then
       call add(errmsg,1,"there were no damping equations added to the regularization object",myname)
       goto 1
    end if
!
    ! (re)allocate for the additional regularization equations (damping) found above
    if(this%neq_total == 0) then
       allocate(this%neq(this%nparam))
       this%neq = 0
       max_neq = maxval(neq_damp)
       allocate(this%neq(this%nparam),this%idx(max_neq,this%nparam),this%coef(max_neq,this%nparam),&
            this%rhs(max_neq,this%nparam))
    else
       ! Need to re-allocate first dimension of equation arrays to the maximum size per parameter AFTER
       ! adding the damping equations.
       ! Since it might happen that damping equations were only found for parameters which do not have
       ! (a lot of) equations yet, we need a further check whether we need to reallocate at all
       max_neq = maxval(this%neq+neq_damp)
       if(max_neq > size(this%idx,1)) then
          ! In this case we need to reallocate! (assume that the first dimensions of this%idx,this%coef,this%rhs
          ! all have the same size!)
          ! First, allocate temporary arrays of new size
          allocate(tmp_idx(max_neq,this%nparam),tmp_coef(max_neq,this%nparam),tmp_rhs(max_neq,this%nparam))
          ! Second, copy all contents from this to temporary arrays (or link pointers, respectively)
          do iparam = 1,this%nparam
             do ieq = 1,this%neq(iparam)
                eq_idx => getVectorPointer(this%idx(ieq,iparam))
                call associateVectorPointer(tmp_idx(ieq,iparam),eq_idx)
                eq_coef => getVectorPointer(this%coef(ieq,iparam))
                call associateVectorPointer(tmp_coef(ieq,iparam),eq_coef)
                tmp_rhs(ieq,iparam) = this%rhs(ieq,iparam)
             end do ! ieq
          end do ! iparam
          ! Third, deallocate arrays from this (attention: do not call dealloc on the vector pointers, since they
          ! were linked to the temporary arrays, just deallocate this%idx, this%coef)
          deallocate(this%idx,this%coef,this%rhs)
          ! Fourth, let the temporary arrays point to this
          this%idx => tmp_idx
          this%coef => tmp_coef
          this%rhs => tmp_rhs
          ! Fifth, let temporary arrays point to null() since they must not be deallocated below before return
          nullify(tmp_idx,tmp_coef,tmp_rhs)
       end if
    end if
!
    ! Finally add the above found damping equations to the existing regularization equations
    do iparam = 1,this%nparam
       ieq_shift = this%neq(iparam)
       do ieq = 1,neq_damp(iparam)
          eq_idx => getVectorPointer(idx_damp(ieq,iparam))
          call associateVectorPointer(this%idx(ieq_shift+ieq,iparam),eq_idx)
          eq_coef => getVectorPointer(coef_damp(ieq,iparam))
          call associateVectorPointer(this%coef(ieq_shift+ieq,iparam),eq_coef)
          this%rhs(ieq_shift+ieq,iparam) = rhs_damp(ieq,iparam)
       end do ! ieq
    end do ! iparam
    ! deallocate local damping equation arrays (attention: do not call dealloc on the vector pointers, since they
    ! were linked to the arrays in this, just deallocate idx_damp,coef_damp,rhs_damp)
    deallocate(idx_damp,coef_damp,rhs_damp)

    ! Update the number of valid equations per parameter and the total number of equations
    this%neq = this%neq + neq_damp
    this%neq_total = this%neq_total + neq_total_damp
!
    if(present(neq_added)) neq_added = neq_total_damp
!
    ! clean up
1   if(allocated(neq_damp)) deallocate(neq_damp)
    if(allocated(idx_damp)) then
       do iparam = 1,size(idx_damp,2)
          do imval = 1,size(idx_damp,1)
             call dealloc(idx_damp(imval,iparam))
          end do ! imval
       end do ! iparam
       deallocate(idx_damp)
    end if
    if(allocated(coef_damp)) then
       do iparam = 1,size(coef_damp,2)
          do imval = 1,size(coef_damp,1)
             call dealloc(coef_damp(imval,iparam))
          end do ! imval
       end do ! iparam
       deallocate(coef_damp)
    end if
    if(allocated(rhs_damp)) deallocate(rhs_damp)
    if(associated(pparam)) deallocate(pparam)
    if(associated(midx_param)) deallocate(midx_param)
    if(allocated(scaling_values_local)) deallocate(scaling_values_local)
  end subroutine addDampingLinearModelRegularization
!------------------------------------------------------------------------
!> \brief add smoothing conditions as equations to incoming linear system
!! \param this linear model smoothing conditions
!! \param KLSE kernel linear system of equations
!! \param errmsg error message
!
  subroutine addToKernelLinearSystemLinearModelRegularization(this,KLSE,errmsg)
    type (linear_model_regularization) :: this
    type (kernel_linear_system) :: KLSE
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=48) :: myname = 'addToKernelLinearSystemLinearModelRegularization'
    real, dimension(:), allocatable :: eq_scale
    integer :: iparam,ieq,ieq_tot,ndata
    real :: rscale
    integer, dimension(:), pointer :: midx_param,idx
    real, dimension(:), pointer :: coef
    character(len=character_length_param), dimension(:), pointer :: pparam
    real, dimension(:,:), pointer :: K,rhs
!
    nullify(midx_param,idx,coef,pparam,K,rhs)
!
    call addTrace(errmsg,myname)
!
    if(.not.(this%initiated)) then
       call add(errmsg,2,"regularization object not yet initiated",myname)
       return
    end if
    if(this%neq_total == 0) then
       call add(errmsg,2,"no smoothing conditions defined yet",myname)
       return
    end if
    if(.not.isInitiated(KLSE)) then
       write(errstr,*) "incoming kernel linear system is not yet initiated (i.e. it has no valid specifications)"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(.nrowreg.KLSE /= this%neq_total) then
       write(errstr,*) "kernel linear system is allocated for ",.nrowreg.KLSE," model regularization conditions, but here ",&
            "are ",this%neq_total," conditions defined"
       call add(errmsg,2,errstr,myname)
       return
    end if
    rhs => .rhs.KLSE
    if(.not.associated(rhs)) then
       call add(errmsg,2,"right-hand-side vector(s) of kernel linear system not defined yet. "//&
            "Please set befor adding smoothing",myname)
       return
    end if
    if(.not.sameModelValues(this%mspace,.dmspace.KLSE)) then
       call add(errmsg,2,"the kernel linear system was not set up for the VERY same model space "//&
            "(requires very same order of model values here)",myname)
       return
    end if
!
    K => .KM.KLSE
    ndata = .ndata.KLSE
!
    allocate(eq_scale(this%neq_total))
    ! dependent on this%scaling_type, define FOR EACH smoothing equation a scaling factor in array eq_scale
    ! scaling values in eq_scale are in following order : do while(nextParam(iparam)); do ieq=i,neq(iparam); enddo; enddo
    select case(this%scaling_type)
    case('')
       eq_scale(:) = 1.
    case('absmax_per_param,overall_factor','absmax_per_param,param_factors')
       ieq_tot = 0
       do iparam = 1,this%nparam
!
          if(associated(midx_param)) deallocate(midx_param)
          if(associated(pparam)) deallocate(pparam)
          allocate(pparam(1)); pparam(1) = this%param(iparam)
          midx_param => getIndxModelValues(this%mspace,param=pparam)
          if(.not.associated(midx_param)) then
             write(errstr,*) "There are no model values in the model for parameter '",trim(this%param(iparam)),&
                  "' although there should be by definition. This error should not occur!!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          ! the following line computes max(abs(K(1:ndata,midx_param))) , but (just a guess!, not actually tested!) performes better 
          ! than max(abs(K(1:ndata,midx_param))), especially for large K (Florian Schumacher, June 2013)
          rscale = max(  sign( maxval(K(1:ndata,midx_param)) , 1. )  ,  sign( minval(K(1:ndata,midx_param)) , 1. )  )
!
          eq_scale(ieq_tot+1:ieq_tot+this%neq(iparam)) = rscale
          ieq_tot = ieq_tot + this%neq(iparam)
       end do ! iparam
       !case('your_scaling_type')
          ! process your type here and define respecive values this%scaling_type
    end select ! this%type_scaling
!
    ! initiate smoothing conditions with zero values
    K(ndata+1:ndata+this%neq_total,:) = 0.
    rhs(ndata+1:ndata+this%neq_total,:) = 0.
    ! now actually add (non-zero coefficients of) scaling equations to Kernel system
    ieq_tot = 0
    do iparam = 1,this%nparam
!
       do ieq = 1,this%neq(iparam)
          idx => getVectorPointer(this%idx(ieq,iparam))
          coef => getVectorPointer(this%coef(ieq,iparam))
!
          ieq_tot = ieq_tot + 1
          K(ndata+ieq_tot,idx) = eq_scale(ieq_tot)*coef
          rhs(ndata+ieq_tot,:) = eq_scale(ieq_tot)*this%rhs(ieq,iparam)
       end do ! ieq
!
    end do ! while (nextParam)
!
1   if(associated(pparam)) deallocate(pparam)
    if(associated(midx_param)) deallocate(midx_param)
    if(allocated(eq_scale)) deallocate(eq_scale)
  end subroutine addToKernelLinearSystemLinearModelRegularization
!------------------------------------------------------------------------
  function requiresQuantityLinearModelRegularization(this,quantity) result(l)
    type (linear_model_regularization), intent(in) :: this
    character(len=*), intent(in) :: quantity
    logical :: l
    l = .false.
    if(this%scaling_type == '') return
    select case(quantity)
    case ('minmaxval_K_per_param')
       l = (this%scaling_type) == 'absmax_per_param,overall_factor' .or. &
            (this%scaling_type) == 'absmax_per_param,param_factors'
    end select
  end function requiresQuantityLinearModelRegularization
!------------------------------------------------------------------------
  subroutine getEquationsVectorPointersLinearModelRegularization(this,eq_indx,eq_coef,eq_rhs,errmsg,&
       ieq_start,ieq_end,minmaxval_K_per_param)
    ! incoming
    type (linear_model_regularization) :: this
    integer, optional :: ieq_start,ieq_end
    real, dimension(:,:), optional :: minmaxval_K_per_param
    ! returning
    type (integer_vector_pointer), dimension(:), pointer :: eq_indx
    type (real_vector_pointer), dimension(:), pointer :: eq_coef
    real, dimension(:), pointer :: eq_rhs
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=51) :: myname = 'getEquationsVectorPointersLinearModelRegularization'
    integer :: ieq1,ieq2,ieq,j,iparam_pmtrz,iparam
    integer, dimension(:), pointer :: idx
    real, dimension(:), pointer :: coef
    real, dimension(:), allocatable :: eq_scale
!
    nullify(idx,coef)
!
    call addTrace(errmsg,myname)
    nullify(eq_indx,eq_coef,eq_rhs)
!
    if(.not.(this%initiated)) then
       call add(errmsg,2,"regularization object not yet initiated",myname)
       return
    end if
!
    if(this%neq_total == 0) then
       call add(errmsg,1,"there are no smoothing equations defined in this smoothing object",myname)
       return
    end if
!
    ieq1 = 1
    if(present(ieq_start)) then
       if(ieq_start >0 .and. ieq_start <= this%neq_total) then
          ieq1 = ieq_start
       else
          write(errstr,*) "incoming requested first index of smoothing equation ",ieq_start,&
               " must be >= 1 and cannot be larger than the total number of constraints in this object (",this%neq_total,")"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ieq2 = this%neq_total
    if(present(ieq_end)) then
       if(ieq_end >=ieq1 .and. ieq_end <= this%neq_total) then
          ieq2 = ieq_end
       else
          write(errstr,*) "incoming requested last index of smoothing equation ",ieq_end,&
               " must be >= first index (",ieq1,") and cannot be larger than the total number of constraints ",&
               "in this object (",this%neq_total,")"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ! is minmaxval_K_per_param present if required?
    select case(this%scaling_type)
    case('absmax_per_param,overall_factor','absmax_per_param,param_factors')
       if(.not.present(minmaxval_K_per_param)) then
          call add(errmsg,2,"optional input minmaxval_K_per_param must be present for scaling types "//&
               "'absmax_per_param,overall_factor','absmax_per_param,param_factors'",myname)
          return
       end if
       if(size(minmaxval_K_per_param,1) /= 2) then
          write(errstr,*) "first dimension of incoming optional array minmaxval_K_per_param is of size ",&
               size(minmaxval_K_per_param,1),"; must be of size 2"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(minmaxval_K_per_param,2) /= numberOfParamModelParametrization(.pmtrz.(this%mspace))) then
          write(errstr,*) "second dimension of incoming optional array minmaxval_K_per_param is of size ",&
               size(minmaxval_K_per_param,2),"; must be as large as the number of model parameters of "//&
               "parametrization '",trim(.pmtrz.(this%mspace)),"' which is ",&
               numberOfParamModelParametrization(.pmtrz.(this%mspace))
          call add(errmsg,2,errstr,myname)
          return
       end if
    end select
!
    ! Define scaling factor for each and every regularization equation, even if some equations are not returned.
    allocate(eq_scale(this%neq_total))
    select case(this%scaling_type)
    case('')
       eq_scale(:) = 1.
    case('absmax_per_param,overall_factor','absmax_per_param,param_factors')
       ! Even if in these cases the additional scaling factor is constant for all model parameters, define 
       ! an array eq_scale (one value for each equation), becaus for other scaling methods, this might be 
       ! different. Sustain maximum flexibility here.
       ieq = 0
       do iparam = 1,this%nparam
          iparam_pmtrz = indexOfParamModelParametrization(.pmtrz.(this%mspace),this%param(iparam))
          if(iparam_pmtrz<=0) then
             write(errstr,*) iparam,"'th parameter '",trim(this%param(iparam)),"' seems not to be a valid "//&
                  "parameter of parametrization '",trim(.pmtrz.(this%mspace)),"', index returned = ",&
                  iparam_pmtrz,". THIS ERROR SHOULD NOT OCCUR!!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          ! do the same operation as in routine addToKernelLinearSystemLinearModelSmoothing to define equation 
          ! scaling for scaling types 'absmax_per_param,overall_factor','absmax_per_param,param_factors'
          eq_scale(ieq+1:ieq+this%neq(iparam)) = &
               max(  sign( minmaxval_K_per_param(2,iparam_pmtrz) , 1. )  ,  &
                     sign( minmaxval_K_per_param(1,iparam_pmtrz) , 1. )  )
          ieq = ieq + this%neq(iparam)
       end do ! iparam
    end select ! this%scaling_type
!
    allocate(eq_indx(ieq2-ieq1+1),eq_coef(ieq2-ieq1+1),eq_rhs(ieq2-ieq1+1))
!
    ieq = 0
    do iparam = 1,this%nparam
!
       ! if the requested starting equation index ieq1 is not in this parameter, 
       ! shift global equation index ieq to end of indices of this parameter and cycle to next parameter
       if(ieq+this%neq(iparam) < ieq1) then
          ieq = ieq+this%neq(iparam)
          cycle
       end if
!
       do j = 1,this%neq(iparam)
          ! if loop comes here, starting index ieq1 must be in this parameter
          ! so cycle until ieq == ieq1 is the requested starting parameter
          if(ieq +1 < ieq1) then
             ieq = ieq + 1
             cycle
          end if
!
          ! if ieq == ieq2 this routine has finished its job
          if(ieq +1 > ieq2) goto 1
!
          ! if loops come here, ieq2 <= ieq <= ieq2, so add smoothing equations to SA,Sb
          idx => getVectorPointer(this%idx(j,iparam))
          coef => getVectorPointer(this%coef(j,iparam))
!
          ieq = ieq + 1
          call allocateVectorPointer(eq_indx(ieq- ieq1+1),size(idx))
          call fillVectorPointer(eq_indx(ieq- ieq1+1),idx,1)
          call allocateVectorPointer(eq_coef(ieq- ieq1+1),size(coef))
          call fillVectorPointer(eq_coef(ieq- ieq1+1),eq_scale(ieq)*coef,1)
          eq_rhs(ieq- ieq1+1) = eq_scale(ieq)*this%rhs(j,iparam)
       end do ! j
!
    end do ! iparam

1   if(allocated(eq_scale)) deallocate(eq_scale)
  end subroutine getEquationsVectorPointersLinearModelRegularization
!------------------------------------------------------------------------
  subroutine getEquationsMatricesLinearModelRegularization(this,SA,Sb,errmsg,ieq_start,ieq_end,&
       minmaxval_K_per_param)
    ! incoming
    type (linear_model_regularization) :: this
    integer, optional :: ieq_start,ieq_end
    real, dimension(:,:), optional :: minmaxval_K_per_param
    ! returning
    real, dimension(:,:), pointer :: SA
    real, dimension(:), pointer :: Sb
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=45) :: myname = 'getEquationsMatricesLinearModelRegularization'
    integer :: ieq1,ieq2,nmval,ieq,j,iparam_pmtrz,iparam
    integer, dimension(:), pointer :: idx
    real, dimension(:), pointer :: coef
    real, dimension(:), allocatable :: eq_scale
!
    nullify(idx,coef)
!
    call addTrace(errmsg,myname)
    nullify(SA,Sb)
!
    if(.not.(this%initiated)) then
       call add(errmsg,2,"regularization object not yet initiated",myname)
       return
    end if
!
    if(this%neq_total == 0) then
       call add(errmsg,1,"there are no smoothing equations defined in this smoothing object",myname)
       return
    end if
    nmval = .nmval.(this%mspace)
    if(nmval <= 0) then
       call add(errmsg,1,"there are no model parameters in the model space object for which "//&
            "this smoothing object was constructed. THIS ERROR SHOULD NOT OCCUR",myname)
       return
    end if
!
    ieq1 = 1
    if(present(ieq_start)) then
       if(ieq_start >0 .and. ieq_start <= this%neq_total) then
          ieq1 = ieq_start
       else
          write(errstr,*) "incoming requested first index of smoothing equation ",ieq_start,&
               " must be >= 1 and cannot be larger than the total number of constraints in this object (",this%neq_total,")"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ieq2 = this%neq_total
    if(present(ieq_end)) then
       if(ieq_end >=ieq1 .and. ieq_end <= this%neq_total) then
          ieq2 = ieq_end
       else
          write(errstr,*) "incoming requested last index of smoothing equation ",ieq_end,&
               " must be >= first index (",ieq1,") and cannot be larger than the total number of constraints ",&
               "in this object (",this%neq_total,")"
          call add(errmsg,2,errstr,myname)
          return
       end if
    end if
!
    ! is minmaxval_K_per_param present if required?
    select case(this%scaling_type)
    case('absmax_per_param,overall_factor','absmax_per_param,param_factors')
       if(.not.present(minmaxval_K_per_param)) then
          call add(errmsg,2,"optional input minmaxval_K_per_param must be present for scaling types "//&
               "'absmax_per_param,overall_factor','absmax_per_param,param_factors'",myname)
          return
       end if
       if(size(minmaxval_K_per_param,1) /= 2) then
          write(errstr,*) "first dimension of incoming optional array minmaxval_K_per_param is of size ",&
               size(minmaxval_K_per_param,1),"; must be of size 2"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(minmaxval_K_per_param,2) /= numberOfParamModelParametrization(.pmtrz.(this%mspace))) then
          write(errstr,*) "second dimension of incoming optional array minmaxval_K_per_param is of size ",&
               size(minmaxval_K_per_param,2),"; must be as large as the number of model parameters of "//&
               "parametrization '",trim(.pmtrz.(this%mspace)),"' which is ",&
               numberOfParamModelParametrization(.pmtrz.(this%mspace))
          call add(errmsg,2,errstr,myname)
          return
       end if
    end select
!
    ! Define scaling factor for each and every equation, even if some equations are not returned.
    allocate(eq_scale(this%neq_total))
    select case(this%scaling_type)
    case('')
       eq_scale(:) = 1.
    case('absmax_per_param,overall_factor','absmax_per_param,param_factors')
       ! Even if in these cases the additional scaling factor is constant for all model parameters, define 
       ! an array eq_scale (one value for each equation), becaus for other scaling methods, this might be 
       ! different. Sustain maximum flexibility here.
       ieq = 0
       do iparam = 1,this%nparam
          iparam_pmtrz = indexOfParamModelParametrization(.pmtrz.(this%mspace),this%param(iparam))
          if(iparam_pmtrz<=0) then
             write(errstr,*) iparam,"'th parameter '",trim(this%param(iparam)),"' seems not to be a valid "//&
                  "parameter of parametrization '",trim(.pmtrz.(this%mspace)),"', index returned = ",&
                  iparam_pmtrz,". THIS ERROR SHOULD NOT OCCUR!!"
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
!
          ! do the same operationo as in routine addToKernelLinearSystemLinearModelSmoothing to define equation 
          ! scaling for scaling type 'absmax_per_param,overall_factor','absmax_per_param,param_factors'
          eq_scale(ieq+1:ieq+this%neq(iparam)) = &
               max(  sign( minmaxval_K_per_param(2,iparam_pmtrz) , 1. )  ,  &
                     sign( minmaxval_K_per_param(1,iparam_pmtrz) , 1. )  )
          ieq = ieq + this%neq(iparam)
       end do ! iparam
    end select ! this%scaling_type
!
    allocate(SA(ieq2-ieq1+1,nmval),Sb(ieq2-ieq1+1))
    SA(:,:) = 0.; Sb(:) = 0.
!
    ieq = 0
    do iparam = 1,this%nparam
!
       ! if the requested starting equation index ieq1 is not in this parameter, 
       ! shift global equation index ieq to end of indices of this parameter and cycle to next parameter
       if(ieq+this%neq(iparam) < ieq1) then
          ieq = ieq+this%neq(iparam)
          cycle
       end if
!
       do j = 1,this%neq(iparam)
          ! if loop comes here, starting index ieq1 must be in this parameter
          ! so cycle until ieq == ieq1 is the requested starting parameter
          if(ieq +1 < ieq1) then
             ieq = ieq + 1
             cycle
          end if
!
          ! if ieq == ieq2 this routine has finished its job
          if(ieq +1 > ieq2) goto 1

          ! if loops come here, ieq2 <= ieq <= ieq2, so add smoothing equations to SA,Sb
          idx => getVectorPointer(this%idx(j,iparam))
          coef => getVectorPointer(this%coef(j,iparam))
!
          ieq = ieq + 1
          SA(ieq- ieq1+1,idx) = eq_scale(ieq)*coef
          Sb(ieq- ieq1+1) = eq_scale(ieq)*this%rhs(j,iparam)
       end do ! j
!
    end do ! iparam
!
1   if(allocated(eq_scale)) deallocate(eq_scale)
  end subroutine getEquationsMatricesLinearModelRegularization
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this linear model smoothing conditions
  function getNeqtotalLinearModelRegularization(this) result(neq)
    type (linear_model_regularization), intent(in) :: this
    integer :: neq
    neq = this%neq_total
  end function getNeqtotalLinearModelRegularization
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this linear model smoothing conditions
  function getScalingTypeLinearModelRegularization(this) result(scltyp)
    type (linear_model_regularization), intent(in) :: this
    character(len=character_length_regscal_type) :: scltyp
    scltyp = this%scaling_type
  end function getScalingTypeLinearModelRegularization
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this linear model smoothing conditions
!
  subroutine deallocateLinearModelRegularization(this)
    type (linear_model_regularization) :: this
    this%initiated = .false.
    call dealloc(this%mspace)
    if(associated(this%param)) deallocate(this%param)
    this%nparam = 0
    this%scaling_type = ''
    this%neq_total = 0
    if(associated(this%neq)) deallocate(this%neq)
    if(associated(this%idx)) deallocate(this%idx)
    if(associated(this%coef)) deallocate(this%coef)
    if(associated(this%rhs)) deallocate(this%rhs)
  end subroutine deallocateLinearModelRegularization
end module linearModelRegularization

