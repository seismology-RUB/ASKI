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
!> \brief define and handle linear equations, which as smoothing conditions are added to a  kernel_linear_system object
!!
!! \details as the smoothing equations will generally be very sparse, avoid storing zeros 
!!  by using vector_pointer structures and only store relevant coefficients of an equation
!!  (along with die indices of the corresponding variables of that equation)
!!
!! \author Florian Schumacher
!! \date June 2013
!
module linearModelSmoothing
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
  interface dealloc; module procedure deallocateLinearModelSmoothing; end interface
  interface operator (.neq.); module procedure getNeqtotalLinearModelSmoothing ; end interface
!
  integer, parameter :: character_length_scaling_type = 50
!
!> \brief type contains indices and values which define equations in a linear system 
  type linear_model_smoothing
     private
     ! equations
     integer :: neq_total = 0 !< total number of equations (sum(ncond))
     character(len=character_length_pmtrz) :: parametrization = '' !< model parametrization, the parameters of which are used to define equation arrays
     integer, dimension(:), pointer :: neq => null() !< number of equations per model parameter name
     type (integer_vector_pointer), dimension(:,:), pointer :: idx => null() !< indices of variables with non-zero coefficient (per model parameter name)
     type (real_vector_pointer), dimension(:,:), pointer :: coef => null() !< coefficients corresponding to indices in idx (per model parameter name)
     real, dimension(:,:), pointer :: rhs => null() !< rhs values of equations (per model parameter name)
     ! scaling
     character(len=character_length_scaling_type) :: scaling_type = '' !< type of scaling this model understands
     real, dimension(:), pointer :: scaling_values => null() !< some numbers used for scaling, corresponding to the type of scaling (type specific)
  end type linear_model_smoothing
!
contains
!
!------------------------------------------------------------------------
!> \brief create simple Laplace smoothing (average over neighbours)
!! \param this linear model smoothing conditions
!
  subroutine createNeighbourAverageLinearModelSmoothing(this,invgrid,dmspace,errmsg,scaling_type,scaling_values,boundary_conditions)
    type (linear_model_smoothing) :: this
    type (inversion_grid) :: invgrid
    type (data_model_space_info) :: dmspace
    type (error_message) :: errmsg
    character(len=*), optional :: scaling_type
    real, dimension(:), optional :: scaling_values
    character(len=*), optional :: boundary_conditions !< IGNORE FOR NOW
    ! local
    character(len=42) :: myname = 'createNeighbourAverageLinearModelSmoothing'
    type (integer_vector_pointer), dimension(:), pointer :: idx_invgrid_nb
    integer :: nidx_invgrid_nb,iparam,iparam_pmtrz,nparam_pmtrz,nparam_dmspace,ieq,nidx_param_cell,n
    character(len=character_length_param) :: param_name
    character(len=character_length_param), dimension(:), pointer :: pparam
    integer, dimension(:), pointer :: idx_param_cell,pcell,nb_idx_p,eq_idx
    integer, dimension(:), allocatable :: pcell_inverse
    real, dimension(:), pointer :: eq_coef
!
    call addTrace(errmsg,myname)
!
    if(this%neq_total /= 0) then
       call add(errmsg,1,"linear model smoothing condititions were already defined. "//&
            "deallocating now before creating new ones",myname)
       call deallocateLinearModelSmoothing(this)
    end if
!
    nparam_dmspace = .nparam.dmspace
    if(nparam_dmspace == 0) then
       call add(errmsg,2,"there are no parameters in model space",myname)
       return
    end if
!
    if(present(scaling_type)) then
       select case(scaling_type)
       case('none') ! do nothing, pass
       case('absmax_per_param,overall_factor')
          if(present(scaling_values)) then
             if(size(scaling_values)>0) then
                allocate(this%scaling_values(1))
                this%scaling_values(1) = scaling_values(1)
             else
                call add(errmsg,2,"in case of scaling_type = 'absmax_per_param,overall_factor', "//&
                     "scaling_values must at least contain one value (the first one is used)",myname)
                return                   
             end if
          else
             call add(errmsg,2,"in case of scaling_type = 'absmax_per_param,overall_factor', "//&
                  "scaling_values must at least contain one value (the first one is used)",myname)
             return
          end if
          this%scaling_type = 'absmax_per_param,overall_factor'
       !case('your_scaling_type')
       !   process your type here and define respecive values this%scaling_type and this%scaling_values
       case default
          call add(errmsg,2,"scaling_type '"//trim(scaling_type)//"' not supported",myname)
       end select ! case scaling_type
    end if
!
    ! FOR NOW, IGNORE boundary_condition!
    ! Later it can control, if on inversion grid boundaries, there should be continuity assumed, or a certain value (e.g. 0 )
    !if(present(boundary_condition)) then
    !   ...
    !endif
!
    call getIndicesFaceNeighboursInversionGrid(invgrid,idx_invgrid_nb)
    if(.not.associated(idx_invgrid_nb)) then
       call add(errmsg,1,"Inversion grid did not return any neighbours. This may suggests that the inversion grid was "//&
            "not properly created",myname)
       return
    end if
    nidx_invgrid_nb = size(idx_invgrid_nb)
!
    ! allocate for maximal number of parametrization parameters
    this%parametrization = .pmtrz.dmspace
    nparam_pmtrz = numberOfParamModelParametrization(this%parametrization)
    allocate(this%neq(nparam_pmtrz),this%idx(nparam_dmspace,nparam_pmtrz),this%coef(nparam_dmspace,nparam_pmtrz),&
         this%rhs(nparam_dmspace,nparam_pmtrz))
!
    this%neq(:) = 0
    do while (nextParamModelParametrization(this%parametrization,param_name,iparam_pmtrz))
!
       ! count the equations added for this parameter
       ieq = 0
!
       ! for this parameter name (e.g. vp) get all invgrid cell indices contained in the model space
       if(associated(pparam)) deallocate(pparam)
       if(associated(pcell)) deallocate(pcell); nullify(pcell)
       if(associated(idx_param_cell)) deallocate(idx_param_cell)
       allocate(pparam(1)); pparam(1) = param_name
       idx_param_cell => getIndxModelParam(dmspace,param=pparam,cell=pcell)
       if(.not.associated(idx_param_cell)) cycle
       nidx_param_cell = size(idx_param_cell)
!
       ! define the inverse mapping pcell_inverse of array pcell which maps an inversion grid cell index to its position in array pcell
       ! i.e.  pcell_inverse(i) = j <=> pcell(j) = i
       if(allocated(pcell_inverse)) deallocate(pcell_inverse)
       allocate(pcell_inverse(maxval(pcell))); pcell_inverse = -1
       do n = 1,nidx_param_cell
          pcell_inverse(pcell(n)) = n
       end do ! i
       
       ! now loop on all model parameters of current param_name
       do iparam = 1,nidx_param_cell
          ! get invgrid cell neighbours of current model parameter
          nb_idx_p => getVectorPointer(idx_invgrid_nb(pcell(iparam)))
          if(.not.associated(nb_idx_p)) cycle
!
          ! only use inversion grid cell neighbours which are also in model space for current parameter name
          n = count(pcell_inverse(nb_idx_p) > 0)
          if(n==0) cycle
!
          ! finally add new equation
          ieq = ieq+1
          ! set indices of variables in equation
          allocate(eq_idx(n+1))
          eq_idx(1:n) = idx_param_cell(pack(pcell_inverse(nb_idx_p),pcell_inverse(nb_idx_p) > 0)) ! dmspace index of valid neighbours
          eq_idx(n+1) = idx_param_cell(iparam) ! dmspace index of center cell
          call associateVectorPointer(this%idx(ieq,iparam_pmtrz),eq_idx); nullify(eq_idx)
          ! set equation coefficients
          allocate(eq_coef(n+1))
          eq_coef(1:n) = 1./real(n) ! coefficient of neighbour cells
          eq_coef(n+1) = -1. ! coefficient of center cell
          call associateVectorPointer(this%coef(ieq,iparam_pmtrz),eq_coef); nullify(eq_coef)
          ! set rhs of equation
          this%rhs(ieq,iparam_pmtrz) = 0.
!
       end do ! iparam
!
       this%neq(iparam_pmtrz) = ieq
!
    end do ! iparam_pmtrz
!
    this%neq_total = sum(this%neq)
!
    ! clean up
    if(this%neq_total == 0) call deallocateLinearModelSmoothing(this)
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
  end subroutine createNeighbourAverageLinearModelSmoothing
!------------------------------------------------------------------------
!> \brief add smoothing conditions as equations to incoming linear system
!! \param this linear model smoothing conditions
!! \param KLSE kernel linear system of equations
!! \param dmspace data model space info for which kernel linear system an smoothing was defined
!! \param errmsg error message
!
  subroutine addToKernelLinearSystemLinearModelSmoothing(this,KLSE,dmspace,errmsg)
    type (linear_model_smoothing) :: this
    type (kernel_linear_system) :: KLSE
    type (data_model_space_info) :: dmspace
    type (error_message) :: errmsg
    ! local
    character(len=400) :: errstr
    character(len=43) :: myname = 'addToKernelLinearSystemLinearModelSmoothing'
    real, dimension(:), allocatable :: eq_scale
    integer :: iparam_pmtrz,ieq,ieq_tot,ndata
    real :: rscale
    integer, dimension(:), pointer :: idx_dmspace,idx
    real, dimension(:), pointer :: coef
    character(len=character_length_param) :: param_name
    character(len=character_length_param), dimension(:), pointer :: pparam
    real, dimension(:,:), pointer :: K,rhs
!
    call addTrace(errmsg,myname)
!
    if(this%neq_total == 0) then
       call add(errmsg,2,"no smoothing conditions defined yet",myname)
       return
    end if
    if(.nsmooth.KLSE /= this%neq_total) then
       write(errstr,*) "kernel linear system is allocated for ",.nsmooth.KLSE," smoothing conditions, but here ",&
            "are ",this%neq_total," smoothing conditions defined"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(.not.associated(.rhs.KLSE)) then
       call add(errmsg,2,"right-hand-side vector(s) of kernel linear system not defined yet. "//&
            "Please set befor adding smoothing",myname)
       return
    end if
    if(this%parametrization /= .pmtrz.dmspace) then
       call add(errmsg,2,"model parametrization '"//trim(this%parametrization)//&
            "' of smoothing differs from model parametrization '"//trim(.pmtrz.dmspace)//"' of model space",myname)
       return
    end if
!
    K => .KM.KLSE
    ndata = .ndata.KLSE
    rhs => .rhs.KLSE
!
    allocate(eq_scale(this%neq_total))
    ! dependent on this%scaling_type, define FOR EACH smoothing equation a scaling factor in array eq_scale
    ! scaling values in eq_scale are in following order : do while(nextParam(iparam)); do ieq=i,neq(iparam); enddo; enddo
    select case(this%scaling_type)
    case('')
       eq_scale(:) = 1.
    case('absmax_per_param,overall_factor')
       ieq_tot = 0
       do while (nextParamModelParametrization(this%parametrization,param_name,iparam_pmtrz))
          if(this%neq(iparam_pmtrz) == 0) cycle
!
          if(associated(idx_dmspace)) deallocate(idx_dmspace)
          if(associated(pparam)) deallocate(pparam)
          allocate(pparam(1)); pparam(1) = param_name
          idx_dmspace => getIndxModelParam(dmspace,param=pparam)
          if(.not.associated(idx_dmspace)) then
             ! actually this cannot be, by construction of equations
             call add(errmsg,2,"for model parameters of parameter '"//trim(param_name)//&
                  "' there are smoothing constraints, although '"//trim(param_name)//&
                  "' is not contained in model space. This suggests that smoothing was defined w.r.t. a different model space",&
                  myname)
             goto 1
          end if
!
          ! the following line computes max(abs(K(1:ndata,idx_dmspace))) , but (just a guess!, not actually tested!) performes better 
          ! than max(abs(K(1:ndata,idx_dmspace))), especially for large K (Florian Schumacher, June 2013)
          rscale = max(  sign( maxval(K(1:ndata,idx_dmspace)) , 1. )  ,  sign( minval(K(1:ndata,idx_dmspace)) , 1. )  )
!
          eq_scale(ieq_tot+1:ieq_tot+this%neq(iparam_pmtrz)) = rscale*this%scaling_values(1)
          ieq_tot = ieq_tot + this%neq(iparam_pmtrz)
       end do ! iparam
       !case('your_scaling_type')
          ! process your type here and define respecive values this%scaling_type and this%scaling_values
    end select ! this%type_scaling
!
    ! initiate smoothing conditions with zero values
    K(ndata+1:ndata+this%neq_total,:) = 0.
    rhs(ndata+1:ndata+this%neq_total,:) = 0.
    ! now actually add (non-zero coefficients of) scaling equations to Kernel system
    ieq_tot = 0
    do while (nextParamModelParametrization(this%parametrization,param_name,iparam_pmtrz))
!
       do ieq = 1,this%neq(iparam_pmtrz)
          idx => getVectorPointer(this%idx(ieq,iparam_pmtrz))
          coef => getVectorPointer(this%coef(ieq,iparam_pmtrz))
!
          ieq_tot = ieq_tot + 1
          K(ndata+ieq_tot,idx) = eq_scale(ieq_tot)*coef
          rhs(ndata+ieq_tot,:) = eq_scale(ieq_tot)*this%rhs(ieq,iparam_pmtrz)
       end do ! ieq
!
    end do ! while (nextParam)
!
1   if(associated(pparam)) deallocate(pparam)
    if(associated(idx_dmspace)) deallocate(idx_dmspace)
    if(allocated(eq_scale)) deallocate(eq_scale)
  end subroutine addToKernelLinearSystemLinearModelSmoothing
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this linear model smoothing conditions
  function getNeqtotalLinearModelSmoothing(this) result(neq)
    type (linear_model_smoothing), intent(in) :: this
    integer :: neq
    neq = this%neq_total
  end function getNeqtotalLinearModelSmoothing
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this linear model smoothing conditions
!
  subroutine deallocateLinearModelSmoothing(this)
    type (linear_model_smoothing) :: this
    if(associated(this%neq)) deallocate(this%neq)
    if(associated(this%idx)) deallocate(this%idx)
    if(associated(this%coef)) deallocate(this%coef)
    if(associated(this%rhs)) deallocate(this%rhs)
    this%neq_total = 0
    this%parametrization = ''
    this%scaling_type = ''
    if(associated(this%scaling_values)) deallocate(this%scaling_values)
  end subroutine deallocateLinearModelSmoothing
end module linearModelSmoothing
