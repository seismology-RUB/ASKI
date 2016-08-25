!----------------------------------------------------------------------------
!	Copyright 2016 Wolfgang Friederich
!
!	This file is part of Gemini II.
!
!	Gemini II is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 2 of the License, or
!	any later version.
!
!	Gemini II is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!> \brief Gemini earth model
!---------------------------------------------
 module geminiEarthModel
    use modelParametrization
    use streamAccess
    use errorMessage
    use fileUnitHandler
    use locatePoint
    use flexibleType
    implicit none
    interface dealloc; module procedure deallocGeminiEarthModel; end interface
    interface operator (.rnod.); module procedure getSelectedRadialNodeGeminiEarthModel; end interface
    interface operator (.znod.); module procedure getSelectedDepthNodeGeminiEarthModel; end interface
    interface operator (.nnod.); module procedure getNnodGeminiEarthModel; end interface
    interface operator (.rearth.); module procedure getEarthRadiusGeminiEarthModel; end interface
    interface operator (.flexpack.); module procedure packFlexibleGeminiEarthModel; end interface
    interface operator (.rpack.); module procedure packDoubleGeminiEarthModel; end interface
    interface operator (.zpack.); module procedure packDoubleComplexGeminiEarthModel; end interface
    type gemini_earth_model
        character (len=80) :: name                                     !< name of earth model
        integer :: nnod                                                !< number of nodes
        integer :: ng                                                  !< number of lateral grid points (for 3D)
        integer :: dampflag                                            !< damping mode (see flnm.h)
        integer :: aniflag                                             !< isotropic (0), transversely isotropic (1)
        double precision :: fref                                       !< reference frequency
        double precision :: rearth                                     !< earth radius
        double precision, dimension(:), pointer :: znod => null()      !< depth nodes in meters
        double precision, dimension(:), pointer :: rnod => null()      !< depth nodes in kilometers
        complex, dimension(:), pointer :: zpa => null()                !< A
        complex, dimension(:), pointer :: zpc => null()                !< C
        complex, dimension(:), pointer :: zpf => null()                !< F
        complex, dimension(:), pointer :: zpl => null()                !< L
        complex, dimension(:), pointer :: zpn => null()                !< N
        complex, dimension(:), pointer :: zkap => null()               !< kappa
        complex, dimension(:), pointer :: zmue => null()               !< mue
        real, dimension(:), pointer :: rho => null()                   !< rho
        real, dimension(:), pointer :: qk => null()                    !< qk = 1/Qk
        real, dimension(:), pointer :: qm => null()                    !< qm = 1/Qm
    end type
!
 contains
!----------------------------------------------------------------------------
!> \brief  Read in gemini earth model from cartesianGreenDisplacement output
!
    subroutine readSAGeminiEarthModel(this,fuh,filename,errmsg)
    type (gemini_earth_model) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    integer :: ierr
    type (file_stream_access) :: fda
    type (group_stream_access) :: root
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    double precision, dimension(:,:), pointer :: rp
    double complex, dimension(:,:), pointer :: zp
    type (flexible), dimension(:), pointer :: ft
    character (len=22) :: myname = 'readSAGeminiEarthModel'
!
    nullify(group,dset,rp,zp,ft)
!
    call addTrace(errmsg,myname)
    ierr = openFileStreamAccess(fda,get(fuh),filename)
    if (ierr /= 0) then
        call add(errmsg,2,'cannot open '//trim(filename),myname)
        call undo(fuh)
        return
    endif
    call readGroupStreamAccess(root,fda)
    call traversePathStreamAccess(root,0,(/ 1, 0 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fda,ft)
    call traversePathStreamAccess(root,0,(/ 2, 0 /),group,dset)
    call readDatasetDouble2DArrayStreamAccess(dset,fda,rp)
    call traversePathStreamAccess(root,0,(/ 3, 0 /),group,dset)
    call readDatasetDoubleComplex2DArrayStreamAccess(dset,fda,zp)
    call createFromPackedGeminiEarthModel(this,ft,rp,zp)
!
!    print *,size(rp,1),size(rp,2),size(zp,1),size(zp,2)
!
!  clean up
!
    deallocate(ft,rp,zp)
    call clearGroupTree(root); call dealloc(fda); call undo(fuh)
!
    end subroutine readSAGeminiEarthModel
!--------------------------------------------------------------------------------------------
!> \brief Create gemini earth model from packed arrays
!
    subroutine createFromPackedGeminiEarthModel(this,ft,rp,zp)
    type (gemini_earth_model) :: this
    type (flexible), dimension(:) :: ft
    double precision, dimension(:,:) :: rp
    double complex, dimension(:,:) :: zp
    integer :: nnod
!
    this%name = ft(1); this%dampflag = ft(2); this%aniflag = ft(3) 
    this%fref = ft(4); this%rearth = ft(5); this%nnod = ft(6)
    if (size(ft) == 7) then; this%ng = ft(7); else; this%ng = 1; endif
    nnod = this%nnod
!
    allocate(this%znod(nnod),this%rnod(nnod),this%rho(nnod),this%qk(nnod),this%qm(nnod))
    allocate(this%zpa(nnod),this%zpc(nnod),this%zpf(nnod),this%zpl(nnod))
    allocate(this%zpn(nnod),this%zkap(nnod),this%zmue(nnod))
!
    this%rnod = rp(:,1)
    this%znod = (this%rearth-rp(:,1))*1.d3
    this%rho = rp(:,2)
    this%qk = rp(:,3)
    this%qm = rp(:,4)
    this%zpa = zp(:,1)
    this%zpc = zp(:,2)
    this%zpf = zp(:,3)
    this%zpl = zp(:,4)
    this%zpn = zp(:,5)
    this%zkap = zp(:,6)
    this%zmue = zp(:,7)
    end subroutine createFromPackedGeminiEarthModel
!--------------------------------------------------------------------------------------------
!> \brief Deallocate gemini_earth_model
!
    subroutine deallocGeminiEarthModel(this)
    type (gemini_earth_model) :: this
    if (associated(this%znod)) deallocate(this%znod)
    if (associated(this%rho)) deallocate(this%rho)
    if (associated(this%qk)) deallocate(this%qk)
    if (associated(this%qm)) deallocate(this%qm)
    if (associated(this%zpa)) deallocate(this%zpa)
    if (associated(this%zpc)) deallocate(this%zpc)
    if (associated(this%zpf)) deallocate(this%zpf)
    if (associated(this%zpl)) deallocate(this%zpl)
    if (associated(this%zpn)) deallocate(this%zpn)
    if (associated(this%zkap)) deallocate(this%zkap)
    if (associated(this%zmue)) deallocate(this%zmue)
    end subroutine deallocGeminiEarthModel
!------------------------------------------------------------------------
!> \brief Get model values for some particular parameter of some particular parametrization
!! \param pmtrz some model parametrization, as defined in module modelParametrization (e.g. "isoVelocity1000")
!! \param param model parameter, which must be one of parametrization pmtrz
!! \param model_values pointer to array of model values on wavefield points. Nullified if there is a problem
!
    function getModelValuesWPGeminiEarthModel(this,pmtrz,param) result(model_values)
    type (gemini_earth_model) :: this
    character(len=*) :: pmtrz,param
    real, dimension(:), pointer :: model_values
    nullify(model_values)
    if(.not.validModelParametrization(pmtrz)) return
    if(.not.validParamModelParametrization(pmtrz,param)) return
    select case(pmtrz)
       case('isoVelocity1000')
          select case(param)
             case('rho'); model_values => getDensityWPGeminiEarthModel(this) ! [g/cm^3]
             case('vp'); model_values => getPVelocityWPGeminiEarthModel(this) ! [km/s]
             case('vs'); model_values => getSVelocityWPGeminiEarthModel(this) ! [km/s]
          end select
       case('isoVelocitySI')
          select case(param)
             case('rho'); model_values => getDensityWPGeminiEarthModel(this) ! [g/cm^3]
             case('vp'); model_values => getPVelocityWPGeminiEarthModel(this) ! [km/s] 
             case('vs'); model_values => getSVelocityWPGeminiEarthModel(this) ! [km/s]
          end select
          model_values = model_values * 1.0e3 ! account for unit factor 1000.0 of the returned values to get Kg/m^3 and m/s
       case('isoLameSI')
          ! could compute mu, lambda from rho,vp,vs and return values
    end select
    end function getModelValuesWPGeminiEarthModel
!------------------------------------------------------------------------------
!> \brief Get density of earth model at wavefield points
!
    function getDensityWPGeminiEarthModel(this) result(rho)
    type (gemini_earth_model) :: this
    real, dimension(:), pointer :: rho
    integer :: j
!
    allocate(rho(this%ng*this%nnod))
    do j = 1,this%nnod
        rho((j-1)*this%ng+1:j*this%ng) = this%rho(j)
    enddo
    end function getDensityWPGeminiEarthModel
!------------------------------------------------------------------------------
!> \brief Get P velocity of earth model at wavefield points
!
    function getPVelocityWPGeminiEarthModel(this) result(vp)
    type (gemini_earth_model) :: this
    real, dimension(:), pointer :: vp
    integer :: j
!
    allocate(vp(this%ng*this%nnod))
    do j = 1,this%nnod
        vp((j-1)*this%ng+1:j*this%ng) = sqrt( real(this%zkap(j)+4.d0*this%zmue(j)/3.d0)/this%rho(j) )
    enddo
    end function getPVelocityWPGeminiEarthModel
!------------------------------------------------------------------------------
!> \brief Get S velocity of earth model at wavefield points
!
    function getSVelocityWPGeminiEarthModel(this) result(vs)
    type (gemini_earth_model) :: this
    real, dimension(:), pointer :: vs
    integer :: j
!
    allocate(vs(this%ng*this%nnod))
    do j = 1,this%nnod
        vs((j-1)*this%ng+1:j*this%ng) = sqrt( real(this%zmue(j))/this%rho(j) )
    enddo
    end function getSVelocityWPGeminiEarthModel
!---------------------------------------------------------------------------------------
!> \brief Get number of radial nodes
!
    function getNnodGeminiEarthModel(this) result(nnod)
    type (gemini_earth_model), intent(in) :: this
    integer :: nnod
    nnod = this%nnod
    end function getNnodGeminiEarthModel
!---------------------------------------------------------------------------------------
!> \brief Get earth radius in km
!
    function getEarthRadiusGeminiEarthModel(this) result(re)
    type (gemini_earth_model), intent(in) :: this
    double precision :: re
    re = this%rearth
    end function getEarthRadiusGeminiEarthModel
!---------------------------------------------------------------------------------------
!> \brief Return pointer to radial nodes in km
!
    function getRnodGeminiEarthModel(this) result(rnod)
    type (gemini_earth_model) :: this
    double precision, dimension(:), pointer :: rnod
    rnod => this%rnod
    end function getRnodGeminiEarthModel 
!---------------------------------------------------------------------------------------
!> \brief Get selected node radius in km
!
    function getSelectedRadialNodeGeminiEarthModel(this,jr) result(rnod)
    type (gemini_earth_model), intent(in) :: this
    integer, intent(in) :: jr
    double precision :: rnod
    rnod = this%rearth-this%znod(jr)*1.d-3
    end function getSelectedRadialNodeGeminiEarthModel
!---------------------------------------------------------------------------------------
!> \brief Get selected node depth in m
!
    function getSelectedDepthNodeGeminiEarthModel(this,jr) result(znod)
    type (gemini_earth_model), intent(in) :: this
    integer, intent(in) :: jr
    double precision :: znod
    znod = this%znod(jr)
    end function getSelectedDepthNodeGeminiEarthModel
!-------------------------------------------------------------------------------
!> \brief Pack flexible part
!
    function packFlexibleGeminiEarthModel(this) result(ft)
    type (gemini_earth_model), intent(in) :: this
    type (flexible), dimension(:), pointer :: ft
    allocate(ft(7))
    ft(1) = this%name; ft(2) = this%dampflag; ft(3) = this%aniflag
    ft(4) = this%fref; ft(5) = this%rearth; ft(6) = this%nnod
    ft(7) = this%ng
    end function packFlexibleGeminiEarthModel
!--------------------------------------------------------------------------------
!> \brief pack double precision part of model into 2D array
!
    function packDoubleGeminiEarthModel(this) result(rp)
    type (gemini_earth_model), intent(in) :: this
    double precision, dimension(:,:), pointer :: rp
    allocate(rp(this%nnod,4))
    rp(:,1) = this%rearth-this%znod*1.d-3
    rp(:,2) = this%rho
    rp(:,3) = this%qk
    rp(:,4) = this%qm
    end function packDoubleGeminiEarthModel
!--------------------------------------------------------------------------------
!> \brief pack double complex part of model into 2D array
!
    function packDoubleComplexGeminiEarthModel(this) result(zp)
    type (gemini_earth_model), intent(in) :: this
    double complex, dimension(:,:), pointer :: zp
    allocate(zp(this%nnod,7))
    zp(:,1) = this%zpa
    zp(:,2) = this%zpc
    zp(:,3) = this%zpf
    zp(:,4) = this%zpl
    zp(:,5) = this%zpn
    zp(:,6) = this%zkap
    zp(:,7) = this%zmue
    end function packDoubleComplexGeminiEarthModel
!----------------------------------------------------------------------------
!> \brief Get index of node nclosest to given depth in meters
!
    function radialNodeIndexGeminiEarthModel(this,zs) result(jr)
    type (gemini_earth_model) :: this
    double precision :: zs
    integer :: jr
    jr = locate(zs,this%nnod,this%znod)
    if (jr > 1) then
        if ((this%znod(jr-1)-zs) .lt. (zs-this%znod(jr))) jr = jr-1
    endif
    end function radialNodeIndexGeminiEarthModel
!----------------------------------------------------------------------
!> \brief Set this%ng
!
    subroutine setNgGeminiEarthModel(this,ng)
    type (gemini_earth_model) :: this
    integer :: ng
    this%ng = ng
    end subroutine setNgGeminiEarthModel
!------------------------------------------------------------------------------
!> \brief Write model to stream-access file
!
    subroutine writeGeminiEarthModel(this,lu,filename,errmsg)
    type (gemini_earth_model) :: this
    integer :: lu
    character(len=*) :: filename
    type (error_message) :: errmsg
    character(len=21) :: myname = 'writeGeminiEarthModel'
    type (file_stream_access) :: fda 
    type (data_stream_access) :: dset
    type (group_stream_access) :: root
    integer :: ios
!
    call addTrace(errmsg,myname)
    ios = createFileStreamAccess(fda,lu,filename)
    if (ios /= 0) then
       call add(errmsg,2,'Open file error: '//trim(filename),myname)
       return
    endif
!
    call createGroupStreamAccess(root,'Root',0,maxsubgroup = 0,maxdset = 3)
    call createDatasetStreamAccess(dset,1,(/ 7 /),T_FLEXIBLE)
    call writeDatasetVectorStreamAccess(dset,fda,.flexpack.this)
    call addDatasetStreamAccess(root,dset)
    call dealloc(dset)
    !        
    call createDatasetStreamAccess(dset,2,(/ this%nnod,4 /),T_DOUBLE)
    call writeDatasetDouble2DArrayStreamAccess(dset,fda,.rpack.this)
    call addDatasetStreamAccess(root,dset)
    call dealloc(dset)
    !        
    call createDatasetStreamAccess(dset,2,(/ this%nnod,7 /),T_DOUBLE_COMPLEX)
    call writeDatasetDoubleComplex2DArrayStreamAccess(dset,fda,.zpack.this)
    call addDatasetStreamAccess(root,dset)
    call dealloc(dset)
    !
    call writeGroupStreamAccess(root,fda)
    call clearGroupTree(root); call dealloc(fda)
    end subroutine writeGeminiEarthModel
!
 end module
