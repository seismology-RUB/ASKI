!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief compute data weights which focus overall sensitivity
!!
!! \details By using a method inspired by G.Backus, F.Gilbert "The Resolving
!!  Power of Gross Earth Data", 1968, weighting coefficients for data samples are computed
!!  such that the overall sensitivity w.r.t. these data weights is large within
!!  a certain subset of the inversion grid, and is close to zero outside that subset. 
!!
!! \author Florian Schumacher
!! \date Nov 2015
!! 
!! \copyright &copy; GNU Public License
!
module kernelFocus
!
  use dataModelSpaceInfo
  use kernelLinearSystem
  use parameterCorrelation
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateKernelFocus; end interface
  interface operator (.a.); module procedure getCoefficientsKernelFocus; end interface
  interface operator (.nfoc.); module procedure getNfocKernelFocus; end interface
!
!> \brief general definition of focus
     type kernel_focus
        private
        logical :: initiated = .false. !< indicates whether object is initiated (i.e. KLSE is set up) or not
        type (kernel_linear_system) :: KLSE !< contains sensitivity kernel matrix for complete data space and complete (unfocussed) model space
!!$        type (data_model_space_info), pointer :: dmspace => null() !< complete data space and complete (unfocussed) model space
!!$        integer :: ndata !< number of data samples in complete data space (.ndata.dmspace)
!!$        integer :: nmval !< number of model values in complete model space (.nmval.dmspace)
!!$        real, dimension(:,:), pointer :: K => null() !< pure kernel matrix without any regularization, for complete data space and complete model space. size(ndata,nmval)
        integer :: nfoc = 0 !< number of focal subvolumes = number of focussing coefficient vectors a
        real, dimension(:,:), pointer :: a => null() !< focussing coefficient vector(s), size(ndata,nfoc)
     end type kernel_focus
!
contains
!------------------------------------------------------------------------
!> \brief 
!! \details 
!!  
!! \param 
!
  subroutine initiateKernelFocus(this,df_measured_data,nfreq_measured_data,ifreq_measured_data,&
       path_sensitivity_kernels,ntot_invgrid,pcorr,lu1,lu2,dmspace,errmsg,&
       apply_event_filter,path_event_filter,apply_station_filter,path_station_filter)
    ! incoming
    type (kernel_focus) :: this
    character(len=*) :: path_sensitivity_kernels
    real :: df_measured_data
    integer, dimension(:) ::  ifreq_measured_data
    integer :: nfreq_measured_data,ntot_invgrid,lu1,lu2
    type (parameter_correlation) :: pcorr
    type (data_model_space_info) :: dmspace
    logical, optional :: apply_event_filter,apply_station_filter
    character(len=*), optional :: path_event_filter,path_station_filter
    ! returning
    type (error_message) :: errmsg
    ! local
    character(len=19) :: myname = 'initiateKernelFocus'
!
    call addTrace(errmsg,myname)
    if(.level.errmsg == 2) return
!
    call initiateSerialKernelLinearSystem(this%KLSE,dmspace,0,0,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call readMatrixSerialKernelLinearSystem(this%KLSE,df_measured_data,&
       nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,&
       pcorr,lu1,lu2,errmsg,apply_event_filter=apply_event_filter,path_event_filter=path_event_filter,&
       apply_station_filter=apply_station_filter,path_station_filter=path_station_filter)
    if(.level.errmsg == 2) goto 1
!
    ! if routine comes here, everything went OK and the object is initiated. So return
    this%initiated = .true.
    return
!
1   call deallocateKernelFocus(this)
  end subroutine initiateKernelFocus
!------------------------------------------------------------------------
!> \brief focus sensitivity on a subset of model values, as suggested by chapter 7 of dissertation by Florian Schumacher 2014 (motivated by Backus and Gilbert, but simplified)
!! \details DOES NOT WORK!!
!!  
!! \param this kernel focus object to be created
!! \param dmspace complete data model space as communicated to initateKernelFocus
!! \param fmspace focussing model subspace, defining focal area in model space which the data (i.e. its sensitivity) should be focussed on
!! \param errmsg error message
!
  subroutine computeSimpleBackusGilbertKernelFocus(this,dmspace,fmspace,errmsg,fintense)
    type (kernel_focus) :: this
    type (data_model_space_info) :: dmspace,fmspace
    type (error_message) :: errmsg
    real, optional :: fintense
    ! local
    character(len=37) :: myname = 'computeSimpleBackusGilbertKernelFocus'
    character(len=400) :: errstr
    type (data_model_space_info) :: dmspace2
    type (kernel_linear_system) :: KLSE_focus
    real, dimension(:), pointer :: k
    real, dimension(:,:), allocatable :: rhs
    real, dimension(:,:), pointer :: sol
!!$real, dimension(:,:), pointer :: KM !! FS FS REMOVE, FOR SETTING KERNEL MATRIX TO ABSOLUTE VALUES
!
    call addTrace(errmsg,myname)
    if(.not.this%initiated) then
       call add(errmsg,2,"kernel focussing object not yet initiated",myname)
       return
    end if
!
    ! dmspace - fmspace = dmspace2
    call subtractModelValues(dmspace,fmspace,dmspace2)
    if(.nmval.dmspace2 == .nmval.dmspace) then
       call add(errmsg,2,"incoming focussing model subspace is not contained in original model space, "//&
            "or both model spaces are empty",myname)
       goto 1
    end if
!
    call initiateSerialKernelLinearSystem(KLSE_focus,dmspace2,0,1,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call copyMatrixKernelLinearSystem(KLSE_focus,this%KLSE,errmsg)
    if(.level.errmsg == 2) goto 1
!!$KM => .KM.KLSE_focus !! FS FS REMOVE, FOR SETTING KERNEL MATRIX TO ABSOLUTE VALUES
!!$KM = abs(KM) !! FS FS REMOVE, FOR SETTING KERNEL MATRIX TO ABSOLUTE VALUES
!
    k => getSumColumns(this%KLSE)
    if(.not.associated(k)) then
       call add(errmsg,2,"error computing sum of columns of full kernel matrix",myname)
       goto 1
    end if
    if(present(fintense)) then
       k = k*fintense
    end if
    call setColumnRegularization(KLSE_focus,k,1,errmsg)
    if(.level.errmsg == 2) goto 1
!
    allocate(rhs(.ncol.KLSE_focus,1))
    rhs(:,:) = 0.
    if(present(fintense)) then
       rhs(.ncol.KLSE_focus,1) = fintense
    else
       rhs(.ncol.KLSE_focus,1) = 1.
    end if
    call setRhsKernelLinearSystem(KLSE_focus,'T',rhs,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call solveSerialKernelLinearSystem(KLSE_focus,errmsg)
    if(.level.errmsg == 2) goto 1
!
    sol => .sol.KLSE_focus
    if(.not.associated(sol)) then
       call add(errmsg,2,"there is no solution in kernel system object after solving it",myname)
       goto 1
    end if
    if(size(sol,2) /= 1) then
       write(errstr,*) "there should be exactly one solution vector, instead there are ",size(sol,2)
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(size(sol,1) /= .nrow.this%KLSE) then
       write(errstr,*) "length ",size(sol,1)," of solution vector does not match number of rows ",&
            .nrow.this%KLSE," of original kernel system"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    this%nfoc = 1
    allocate(this%a(size(sol,1),1))
    this%a = sol
!
1   call dealloc(dmspace2)
    call dealloc(KLSE_focus)
    if(associated(k)) deallocate(k)
  end subroutine computeSimpleBackusGilbertKernelFocus
!------------------------------------------------------------------------
!> \brief focus sensitivity on a subset of model values, as suggested originally by the Backus and Gilbert paper
!! \details 
!!  
!! \param this kernel focus object to be created
!! \param dmspace complete data model space as communicated to initateKernelFocus
!! \param fmspace focussing model subspace, defining focal area in model space which the data (i.e. its sensitivity) should be focussed on
!! \param errmsg error message
!
  subroutine computeOriginalBackusGilbertKernelFocus(this,dmspace,fmspace,invgrid,errmsg)
    use serialLinearSystem
    type (kernel_focus) :: this
    type (data_model_space_info) :: dmspace,fmspace
    type (inversion_grid) :: invgrid
    type (error_message) :: errmsg,errmsg2
    ! local
    character(len=37) :: myname = 'computeSimpleBackusGilbertKernelFocus'
    character(len=400) :: errstr
    type (data_model_space_info) :: dmspace_outside
    real, dimension(:,:), pointer :: K
    integer :: ndata,K_nmval,nmval_outside,npacked
    integer :: icell,ipacked,icol,irow
    integer, dimension(:), pointer :: indx_param_outside,indx,indxcell
    real :: mean_volume_cell_outside,vol
    real, dimension(:), allocatable :: K_icol
    real, dimension(:,:), allocatable :: C
    type (kernel_linear_system) :: KLSE_focus
    ! serial packed linear system
    type (serial_linear_system) :: LSE
    real, dimension(:), pointer :: A_packed
    real, dimension(:,:), pointer :: rhs
    real, dimension(:,:), pointer :: sol
!
    call addTrace(errmsg,myname)
    if(.not.this%initiated) then
       call add(errmsg,2,"kernel focussing object not yet initiated",myname)
       return
    end if
!
    ! dmspace - fmspace -> dmspace_outside
    call subtractModelValues(dmspace,fmspace,dmspace_outside)
    if(.nmval.dmspace_outside == .nmval.dmspace) then
       call add(errmsg,2,"incoming focussing model subspace is not contained in original model space, "//&
            "or both model spaces are empty",myname)
       goto 1
    end if
!
    ! GET ALL MODEL VALUES INDICES (i.e. indices of columns of kernel matrix) WHICH CORRESPOND TO  O U T S I D E  FOCUSSING VOLUME
    indx_param_outside => getIndxModelValues(dmspace_outside,dmspace)
    if(.not.associated(indx_param_outside)) then
       call add(errmsg,2,"no indices returned for model space outside focussing area",myname)
       goto 1
    end if
    nmval_outside = size(indx_param_outside)
    write(errstr,*) "number of model parameters outside focussing volume: ",nmval_outside
write(*,*) "ORIG B&G FOCUSSING: number of model parameters outside focussing volume: ",nmval_outside
    call add(errmsg,0,errstr,myname)
!
    ! GET TOTAL VOLUME OF CELLS OUTSIDE FOCUSSING VOLUME, USED TO CORRECT PREINTEGRATED KERNEL SUMMATION FOR APPROXIMATE BACKUS&GILBERT INTEGRAL FORMULAS
    mean_volume_cell_outside = 0.
    nullify(indxcell)
    indx => getIndxModelValues(dmspace_outside,cell=indxcell)
    if(.not.associated(indx)) then
       call add(errmsg,2,"ERROR NUMBER 1",myname)
       goto 1
    end if
    if(.not.associated(indxcell)) then
       call add(errmsg,2,"ERROR NUMBER 1b",myname)
       goto 1
    end if
    if(size(indxcell) /= nmval_outside) then
       call add(errmsg,2,"ERROR NUMBER 2",myname)
       goto 1
    end if
    deallocate(indx)
write(*,*) "ORIG B&G FOCUSSING: calculating mean invgrid cell volume outside focussing volume"
    call new(errmsg2,myname)
    do icell = 1,nmval_outside
       !! subroutine getVolumeCellInversionGrid(this,icell,volume,errmsg)
       call getVolumeCellInversionGrid(invgrid,indxcell(icell),vol,errmsg2)
       if(.level.errmsg2 /= 0) then
          call print(errmsg2)
          write(errstr,*) "error getting colume of invgrid cell ",indxcell(icell)
          call add(errmsg,2,errstr,myname)
          return
       end if
       mean_volume_cell_outside = mean_volume_cell_outside + vol 
    end do ! icell
    mean_volume_cell_outside = mean_volume_cell_outside / real(nmval_outside)
    call dealloc(errmsg2)
    if(mean_volume_cell_outside <= 0.) then
       call add(errmsg,2,"ERROR NUMBER 3",myname)
       goto 1
    end if
!
    write(errstr,*) "mean invgrid cell volume outside focussing volume: ",mean_volume_cell_outside
write(*,*) "mean invgrid cell volume outside focussing volume: ",mean_volume_cell_outside
    call add(errmsg,0,errstr,myname)
!!$    ! INVERT THE VALUE (need to divide by that later on)
!!$    mean_volume_cell_outside = 1./mean_volume_cell_outside  !! DO NOT DIVIDE KERNEL INTEGRALS BY VOLUME, INSTEAD MULTIPLY THE WHOLE MATRIX EQUATION BY mean_volume_cell_outside (see below)
!
!
!
!
!!########### PERFORM THE MULTIPLICATION K*K^T MANUALLY,  T H I S  P E R F O R M S  V E R Y  B A D L Y
!!$    ndata = .ndata.(this%KLSE)
!!$    K_nmval = .nmval.(this%KLSE)
!!$    K => .KM.(this%KLSE)
!!$!
!!$    ! ALLOCATE AND DEFINE UPPER-DIAGONAL PACKED BACKUS&GILBERT FOCUSSING MATRIX, 
!!$    ! FOCUSSING MATRIX HAS SIZE (ndata+1)-times-(ndata+1) (additional Lagrange parameter)
!!$    ! HENCE, PACKED VECTOR HAS SIZE (ndata+2)*(ndata+1)/2 = sum_{i=1}^(ndata+1) i
!!$    npacked = ((ndata+2)*(ndata+1))/2
!!$    allocate(A_packed(npacked))
!!$write(*,*) "ORIG B&G FOCUSSING: SETTING UP PACKED FOCUSSING MATRIX, total number of elements: ",npacked
!!$    ! packed vector contains upper diagonal B&G focussing matrix COLUMEWISELY!
!!$    ipacked = 0
!!$    ! first do the columns containing the integrals over kernel products
!!$    allocate(K_icol(nmval_outside))
!!$    do icol = 1,ndata
!!$       K_icol = K(icol,indx_param_outside)
!!$       do irow = 1,icol
!!$          ipacked = ipacked + 1
!!$          A_packed(ipacked) = sum(K(irow,indx_param_outside)*K_icol)
!!$!
!!$          ! WRITE TO SCREEN WHICH ENTRY
!!$          !if(mod(ipacked,npacked/100)==0) then
!!$          if(mod(ipacked,10000)==0) then
!!$             write(*,*) "ORIG B&G FOCUSSING: SET UP PACKED FOCUSSING MATRIX, ENTRY ",ipacked," OUT OF ",npacked
!!$          end if
!!$       end do ! irow
!!$    end do ! icol
!!$    ! then do the last column (of size ndata) EXCEPT LAST ENTRY!
!!$    do irow = 1,ndata
!!$       ipacked = ipacked + 1
!!$       A_packed(ipacked) = sum(K(irow,1:K_nmval))
!!$!
!!$       ! WRITE TO SCREEN WHICH ENTRY
!!$       !if(mod(ipacked,npacked/100)==0) then
!!$       if(mod(ipacked,10000)==0) then
!!$          write(*,*) "ORIG B&G FOCUSSING: SET UP PACKED FOCUSSING MATRIX, ENTRY ",ipacked," OUT OF ",npacked
!!$       end if
!!$    end do ! irow
!!$    ! then do last entry (which is zero)
!!$    ipacked = ipacked + 1
!!$    A_packed(ipacked) = 0.
!
!
!
!!########### PERFORM THE MULTIPLICATION K*K^T BY LEVEL3-BLAS ROUTINES
write(*,*) 'ORIG B&G FOCUSSING: copying kernel matrix, K" = K( : , j\notin J_F) which misses the focussing columns'
    call initiateSerialKernelLinearSystem(KLSE_focus,dmspace_outside,0,0,errmsg)
    if(.level.errmsg == 2) goto 1
!
    call copyMatrixKernelLinearSystem(KLSE_focus,this%KLSE,errmsg)
    if(.level.errmsg == 2) goto 1
!
    K => .KM.KLSE_focus
    ndata = .ndata.(KLSE_focus)
    K_nmval = .nmval.(KLSE_focus)
!
!
write(*,*) "ORIG B&G FOCUSSING: SETTING UP PACKED FOCUSSING MATRIX, total number of elements: ",((ndata+2)*(ndata+1))/2
write(*,*) 'ORIG B&G FOCUSSING: perform K"*(K"^T) by the level 3 BLAS routine SSYRK now'
!
    allocate(C(ndata,ndata))
    !! SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
    !!!call SSYRK ( 'U', 'N', ndata, K_nmval, mean_volume_cell_outside, K, ndata, 0.0, C, ndata )
    call SSYRK ( 'U', 'N', ndata, K_nmval, 1.0, K, ndata, 0.0, C, ndata ) !! DO NOT DIVIDE KERNEL INTEGRALS BY VOLUME, INSTEAD MULTIPLY THE WHOLE MATRIX EQUATION BY mean_volume_cell_outside (see below)

    ! ALLOCATE AND DEFINE UPPER-DIAGONAL PACKED BACKUS&GILBERT FOCUSSING MATRIX, 
    ! FOCUSSING MATRIX HAS SIZE (ndata+1)-times-(ndata+1) (additional Lagrange parameter)
    ! HENCE, PACKED VECTOR HAS SIZE (ndata+2)*(ndata+1)/2 = sum_{i=1}^(ndata+1) i
    npacked = ((ndata+2)*(ndata+1))/2
    allocate(A_packed(npacked))
!
write(*,*) "ORIG B&G FOCUSSING: done, copying to packed vector now"
    ! COPY 2-rank OUTPUT MATRIX REFERENCED BY SSYRK TO 1-rank PACKED MATRIX VECTOR A_packed
    ipacked = 0
    ! packed vector contains upper diagonal B&G focussing matrix COLUMEWISELY!
    ! first do the columns containing the integrals over kernel products, i.e. upper triangular matrix of C
    do icol = 1,ndata
       do irow = 1,icol
          ipacked = ipacked + 1
          A_packed(ipacked) = C(irow,icol)
       end do ! irow
    end do ! icol
write(*,*) "ORIG B&G FOCUSSING: done, defining last column of packed vector now"
    ! then do the last column (of size ndata+1) EXCEPT LAST ENTRY!
    K => .KM.this%KLSE
    ! FOR PROPER ORIGINAL B&G FOCUSSING WITH MORE THAN 1 MODEL PARAMETER (e.g. vs AND vp), DO NOT USE  1:K_nmval  BELOW, BUT
    ! AN INDEX VECTOR DEFINING ALL COLUMNS WHICH CORRESPOND TO THE FOCUSSED PARAMETER (e.g. all invgrid cells of the focussing parameter vs)
    ! FOR NOW, ASSUME ONE PARAMETER (e.g. only vs) IN COMPLETE DATA SPACE
    K_nmval = .nmval.(this%KLSE)
    do irow = 1,ndata
       ipacked = ipacked + 1
       !!A_packed(ipacked) = sum(K(irow,1:K_nmval))
       A_packed(ipacked) = mean_volume_cell_outside*sum(K(irow,1:K_nmval))  !! DO NOT DIVIDE KERNEL INTEGRALS BY VOLUME, INSTEAD MULTIPLY THE WHOLE MATRIX EQUATION BY mean_volume_cell_outside
    end do ! irow
    ! then do last entry (which is zero)
    ipacked = ipacked + 1
    A_packed(ipacked) = 0.
write(*,*) "ORIG B&G FOCUSSING: done"

    ! DEFINE rhs OF PACKED SYSTEM
    allocate(rhs(ndata+1,1))
    rhs(:,:) = 0.
    !!rhs(ndata+1,1) = 1.
    rhs(ndata+1,1) = mean_volume_cell_outside  !! DO NOT DIVIDE KERNEL INTEGRALS BY VOLUME, INSTEAD MULTIPLY THE WHOLE MATRIX EQUATION BY mean_volume_cell_outside
!    
write(*,*) "ORIG B&G FOCUSSING: creating packed serial linear system"
    !! subroutine createPackedSerialLinearSystem(this,nrow,ncol,nrhs,A_packed,UPLO,b)
    call createPackedSerialLinearSystem(LSE,ndata+1,ndata+1,1,A_packed,'U',rhs)
!
write(*,*) "ORIG B&G FOCUSSING: solving packed serial linear system now"
    !! function solvePackedSymmetricSerialLinearSystem(this,x,overwrite_A) result(errmsg)
    errmsg2 = solvePackedSymmetricSerialLinearSystem(LSE,sol,overwrite_A=.true.)
    if(.level.errmsg2 /= 0) call print(errmsg2)
    if(.level.errmsg2 == 2) then
       call add(errmsg,2,"ERROR NUMBER 3",myname)
       goto 1
    end if
    this%nfoc = 1
    allocate(this%a(ndata,1))
    this%a(:,1) = sol(1:ndata,1)
!
write(*,*) "ORIG B&G FOCUSSING: done"
1   call dealloc(dmspace_outside)
    call dealloc(LSE)
    if(associated(indx_param_outside)) deallocate(indx_param_outside)
    if(associated(indx)) deallocate(indx)
    if(associated(indxcell)) deallocate(indxcell)
    if(associated(sol)) deallocate(sol)
    if(associated(rhs)) deallocate(rhs)
    if(allocated(K_icol)) deallocate(K_icol)
    if(allocated(C)) deallocate(C)
  end subroutine computeOriginalBackusGilbertKernelFocus
!------------------------------------------------------------------------
!> \brief get focussed sensitivity kernel values
!! \details for given coefficients (or coefficient matrices) a=this%a and full sensitivity matrix K=.KM.KLSE, return
!!  the focussed sensitivity vector(s) K^T*a
!! \param this kernel focus object
!
  subroutine getKernelFocus(this,kf,errmsg)
    type (kernel_focus) :: this
    real, dimension(:,:), pointer :: kf
    type (error_message) :: errmsg
    ! local
    character(len=14) :: myname = 'getKernelFocus'
    real, dimension(:,:), pointer :: KM
    integer :: m,k,n
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"kernel focussing object not yet initiated",myname)
       return
    end if
    if(.not.associated(this%a)) then
       call add(errmsg,2,"there were no focussing coefficients computed yet",myname)
       return
    end if
!
    ! IN THE FUTURE: CALL A ROUTINE IN kernelLinearSystem TO DO THE MATRIX PRODUCT
    !SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
    KM => .KM.(this%KLSE)
    if(.not.associated(KM)) then
       call add(errmsg,2,"Although kernel focussing object initiated, there is no kernel matrix of "//&
            "full kernel system. This should not happen, modules are corrupt",myname)
       return
    end if
    ! as we use KM^T in the multiplication, m,k must represent those numbers such that KM^T is an m-by-k matrix
    m = size(KM,2)
    k = size(KM,1)
    ! as we use this%a (not transposed), n must represent a number such that this%a is a k-by-n matrix
    n = this%nfoc
    ! allocate the solution array
    nullify(kf)
    allocate(kf(m,n))
    call SGEMM('T','N',m,n,k,1.,KM,k,this%a,m,0.,kf,m)
  end subroutine getKernelFocus
!------------------------------------------------------------------------
  subroutine getMeanKernelFocus(this,k,errmsg)
    type (kernel_focus) :: this
    real, dimension(:), pointer :: k
    type (error_message) :: errmsg
    ! local
    character(len=18) :: myname = 'getMeanKernelFocus'
    real, dimension(:,:), pointer :: KM
    integer :: nrow,i
!
    call addTrace(errmsg,myname)
!
    if(.not.this%initiated) then
       call add(errmsg,2,"kernel focussing object not yet initiated",myname)
       return
    end if
!
    ! IN THE FUTURE: CALL A ROUTINE IN kernelLinearSystem TO GET THE SUM OF ROWS
    KM => .KM.(this%KLSE)
    if(.not.associated(KM)) then
       call add(errmsg,2,"Although kernel focussing object initiated, there is no kernel matrix of "//&
            "full kernel system. This should not happen, modules are corrupt",myname)
       return
    end if
!
    nrow = size(KM,1)
!
    nullify(k)
    allocate(k(size(KM,2)))
    k = 0.
    do i = 1,nrow
       k = k + KM(i,:)
    end do ! i
    k = k*(1./nrow)
  end subroutine getMeanKernelFocus
!------------------------------------------------------------------------
!> \brief get focussing coefficients this%a
!! \param this kernel focus
!! \param a focussing coefficients this%a
!
  function getCoefficientsKernelFocus(this) result(a)
    type (kernel_focus), intent(in) :: this
    real, dimension(:,:), pointer :: a
    a => this%a
  end function getCoefficientsKernelFocus
!------------------------------------------------------------------------
!> \brief get number of focal subvolumes
!! \param this kernel focus
!! \param nfoc number of focal subvolumes
!
  function getNfocKernelFocus(this) result(nfoc)
    type (kernel_focus), intent(in) :: this
    integer :: nfoc
    nfoc = this%nfoc
  end function getNfocKernelFocus
!------------------------------------------------------------------------
!> \brief deallocate object
!! \param this kernel focus object to be deallocated
!
  subroutine deallocateKernelFocus(this)
    type (kernel_focus) :: this
    call dealloc(this%KLSE)
    if(associated(this%a)) deallocate(this%a)
    this%nfoc = 0
    this%initiated = .false.
  end subroutine deallocateKernelFocus
!
end module kernelFocus
