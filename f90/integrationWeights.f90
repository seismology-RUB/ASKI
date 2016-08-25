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
!> \brief Module computing, writing, reading, dealing with integration weights
!!
!! \details The integration weights constitute linear integration rules, which integrate functions 
!!  living on the wavefield points (such as sensitivity kernels) onto the inversion grid cells.
!!  Different types of integration weights are supported by this module:
!!    0 = all weights the same, weight = 1/number_of_points_in_box, i.e. no integration, just build average
!!    1 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 1
!!    2 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 2, i.e. approximation order 3 (?)
!!    3 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 3, i.e. approximation order 4 (?)
!!    4 = for each cell compute the highest possible order of SDI integration (trying types 3,2,1 (in that order) until computation was successful)
!!    5 = average of function values, multiplied with volume of box (i.e. linear integration)
!!    6 = external integration weights, to be used along with a suitable inversion grid (e.g. specfem3dInversionGrid)
!
module integrationWeights
!
  use wavefieldPoints
  use inversionGrid
  use vectorPointer
  use realloc
  use mathConstants
!
  implicit none
!
  private :: computeWeightAverage,computeWeightLinearIntegration,computeWeightLinearSDI,&
       computeWeightLinearSDIHex,computeWeightLinearSDITet,computeWeightQuadraticSDI,&
       computeWeightQuadraticSDIHex,computeWeightQuadraticSDITet,computeWeightCubicSDI,&
       computeWeightCubicSDIHex,computeWeightCubicSDITet
!
  interface dealloc; module procedure deallocateIntegrationWeights; end interface
  interface getFilledCells; module procedure getIndicesFilledCellsIntegrationWeights; end interface
  interface getEmptyCells; module procedure getIndicesEmptyCellsIntegrationWeights; end interface
  interface emptyCell; module procedure cellIsEmptyIntegrationWeights; end interface
  interface anyCellEmpty; module procedure anyCellIsEmptyIntegrationWeights; end interface
  interface getWpInside; module procedure getIndicesWpInsideInvgridIntegrationWeights; end interface
  interface getWpOutside; module procedure getIndicesWpOutsideInvgridIntegrationWeights; end interface
  interface anyWpOutside; module procedure anyWpIsOutsideInvgridIntegrationWeights; end interface
  interface operator (.wpidx.); module procedure getSelectedWpIdxIntegrationWeights; end interface
  interface operator (.weight.); module procedure getSelectedIntegrationWeights; end interface
  interface operator (.ncell.); module procedure getNtotInvgridIntegrationWeights; end interface
  interface operator (.nwp.); module procedure getNtotWpIntegrationWeights; end interface
!
  type integration_weights
     private
     integer :: ntot_wp = 0
     type (integer_vector_pointer), dimension(:), pointer :: wp_idx => null() !< Global indices of wavefield points located in given inversion cell
     logical, dimension(:), pointer :: empty_cell => null() !< indicates for each invgrid cell if it is empty
     type (real_vector_pointer), dimension(:), pointer :: weight => null() !< Integration weights for wavefield points in given inversion cell
     integer, dimension(:), pointer :: err_level_weights => null() !< stores for each cell error level of weight computation
     integer, dimension(:), pointer :: type_weights => null() !< stores for each cell the type of integration weights which are used
     real :: unit_factor = 0 !< unit factor of integration weights
  end type integration_weights
!
contains
!-------------------------------------------------------------------------------
!> \brief check if intw_type is valid, i.e. supported by this module
  function validTypeIntegrationWeights(itype,err) result(l)
    integer :: itype
    logical :: l
    type (error_message), optional :: err
    character(len=400) :: errstr
    if(present(err)) call addTrace(err,'validTypeIntegrationWeights')
    select case(itype)
    case(0,1,2,3,4,5,6)
    ! Type of integration weights: 
    !   0 = all weights the same, weight = 1/number_of_points_in_box, i.e. no integration, just build average
    !   1 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 1
    !   2 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 2, i.e. approximation order 3 (?)
    !   3 = Scattered Data Integration, as in D. Levin [1999], polynomial degree 3, i.e. approximation order 4 (?)
    !   4 = for each cell compute the highest possible order of SDI integration (trying degrees 3,2,1 , in that order). if not even SDI degree 1 is successful, choose weights of type 5
    !   5 = average of function values, multiplied with volume of box (i.e. some sort of linear integration)
    !   6 = external integration weights, to be used along with a suitable inversion grid (e.g. specfem3dInversionGrid)
       l = .true.
    case default
       l = .false.
       if(present(err)) then
          write(errstr,*) "incoming type ",itype," of integration weights is not supported: "//&
               "valid types are 0, 1, 2, 3, 4, 5, 6"
          call add(err,2,errstr,'validTypeIntegrationWeights')
       end if
    end select
  end function validTypeIntegrationWeights
!-------------------------------------------------------------------------------
!> \brief compute integration weights
!! \details First go through all wavefield points (given in inversion grid coordinates
!!  through function getTransformedWavefieldPoints) and locate them inside a specific inversion
!!  grid cell. Give warnings if there are wavefield points which lie outside the
!!  inversion grid. Then go through the inversion grid cells and compute the requested
!!  type of integration weights (defined by intw_type). Give warnings if there are inversion
!!  grid cells which do not contain any wavefield points. For those inversion grid cells,
!!  define an artificial wavefield point (wavefield point index 1 ) with zero integration weight
!!  such that the integral of arbitrary functions (e.g. kernels) over that invgrid cell are
!!  always zero, allowing for holes, caves (e.g. tunnels) inside a regular inversion grid.
!! \param this integration weights
!! \param intw_type integration weights type code
!! \param wp wavefield points
!! \param invgrid inversion grid
!! \param errmsg error message
!! \return error message
!
  subroutine createIntegrationWeights(this,intw_type,wp,invgrid,errmsg)
    ! incoming
    type (integration_weights) :: this
    integer :: intw_type
    type (wavefield_points) :: wp
    type (inversion_grid) :: invgrid
    ! outgoing
    type (error_message) :: errmsg
    ! local
    type (error_message) :: errmsg2
    character (len=400) :: errstr
    character (len=24) :: myname = 'createIntegrationWeights'
    real, dimension(:), pointer :: c1,c2,c3
    integer, dimension(:), pointer :: idx
    real, dimension(:), pointer :: weight,jacobian
    real, dimension(:), allocatable :: x,y,z
    real :: volume,uf_intw,uf_wp
    double precision :: sum_intw
    integer :: ncell_invgrid,icell,cnt_empty_cells,size_idx,type_standard_cell
!
    nullify(c1,c2,c3,idx,weight,jacobian)
    call addTrace(errmsg,myname)
    if(this%ntot_wp /= 0) then
       call add(errmsg,1,"there are already integration weights defined, "//&
            "deallocating now before creating new ones",myname)
       call deallocateIntegrationWeights(this)
    end if
!
!  compute unit factor of integration weights (only dependent on unit factor of wavefield points and on type integration weights)
!
    call getUnitFactorWavefieldPoints(wp,uf_wp,errmsg)
    if(.level.errmsg == 2) goto 1
    if(uf_wp <= 0) then
       write(errstr,*) "The unit factor of wavefield points = ",uf_wp,&
            " is not strictly positive. This is not supported by ASKI, there seems to be some problem."
       call add(errmsg,2,trim(errstr),myname)
       goto 1
    end if
    select case(intw_type)
    case ( 0 )
       ! When building the average, no actual integration is done, hence there is no volume element and no 
       ! unit to account for. In the concept of ASKI, this can be achieved by setting this%unit_factor to 
       ! 1.0 , i.e. the neutral value of multiplication
       this%unit_factor = 1.0
    case ( 1, 2, 3, 4, 5, 6 )
       ! For types 1, 2, 3, 4, a transformation to a standard cell and jacobian values
       ! must be provided by the inversion grid module. Even if those are all weights for 3D integration, 
       ! the unit of the volume element does not necessarily need to be the cube of the unit of the 
       ! wavefield point coordinates:
       ! E.g. the inversion grid could allow for different units of wavefield
       ! points, while mandating an inversion grid specific unit of the volume element (spherical grids will
       ! e.g. use km^3, but could allow for wavefield points given in m etc).
       ! Weights of type 5 will use the cell volume which is provided by the inversion grid module (and 
       ! should naturally have the unit of the volume element).
       ! Weights of type 6 will be provided directly by the inversion grid module, hence the unit factor
       ! of the volume element should in any case be provided by the inversion grid.
       ! Hence, we ask the inversion grid object, to which value this%unit_factor should be set here.
       call getUnitFactorOfVolumeElementInversionGrid(invgrid,uf_wp,uf_intw,errmsg)
       if(.level.errmsg == 2) goto 1
       if(uf_intw <= 0) then
          write(errstr,*) "The unit factor of integration weights returned by inversion grid = ",uf_intw,&
               " is not strictly positive. This is not supported by ASKI, there seems to be some problem."
          call add(errmsg,2,trim(errstr),myname)
          goto 1
       end if
       this%unit_factor = uf_intw
    end select ! case intw_type
!
!  locate wavefield points inside inversion grid cells
!
    call getWavefieldPoints(wp,c1,c2,c3)
    if(.not.associated(c1)) then
       call add(errmsg,2,"there were no wavefield points returned by 'getWavefieldPoints'",myname)
       goto 1
    end if
    this%ntot_wp = .ntot.wp
!
    call locateWpInsideInversionGrid(invgrid,c1,c2,c3,uf_wp,this%wp_idx,errmsg)
    if(.level.errmsg == 2) goto 1
    if(.not.associated(this%wp_idx)) then
       call add(errmsg,2,"wavefield points could not be located inside inversion grid",myname)
       goto 1
    end if
    ncell_invgrid = .ncell.invgrid
    if(size(this%wp_idx) /= ncell_invgrid) then
       write(errstr,*) "number of returned cells with wavefield points localizations ",size(this%wp_idx),&
            " does not match number of inversion grid cells ",ncell_invgrid
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
! create this%empty_cell, this%err_level_weights, this%type_weights and compute this%weights
!
    allocate(this%empty_cell(ncell_invgrid),this%err_level_weights(ncell_invgrid),&
         this%type_weights(ncell_invgrid),this%weight(ncell_invgrid))
    this%err_level_weights = 2 ! initiate all cells to erroneous, then empty cells will always be erroneous
    nullify(jacobian)
!
    cnt_empty_cells = 0
    sum_intw = 0.d0
    do icell = 1,ncell_invgrid
!
       idx => getVectorPointer(this%wp_idx(icell))
       this%empty_cell(icell) = .not.associated(idx)
!
       if(this%empty_cell(icell)) then
!
          cnt_empty_cells = cnt_empty_cells + 1
          ! in case of this cell being empty, define one artificial wavefield point contained in that cell 
          ! (wp index 1, which always exists) and a corresponding integration weight = 0.
          ! this is done to be able to do integration operations, e.g. to compute kernels, for ALL inversion 
          ! grid cells regardless of a cell being empty or not. Array this%empty_cell then indicates how these 
          ! integration operations should be interpreted (note that integration of any function over empty 
          ! cells always yields value 0.)
          allocate(idx(1)); idx(1) = 1
          call associateVectorPointer(this%wp_idx(icell),idx)
          allocate(weight(1)); weight(1) = 0.
          call associateVectorPointer(this%weight(icell),weight)
!
       else ! this%empty_cell(icell)
          ! if this cell is not empty, compute integration weights
!
          size_idx = size(idx)
!
          select case(intw_type)
          ! for weights of type 0, we do not need to transform the wavefield points to a standard cell
          case ( 0 ) ! intw_type
             weight => computeWeightAverage(size_idx)
             this%err_level_weights(icell) = 0
             this%type_weights(icell) = 0
!
          ! for weights of type 1,2,3,4 first transform the wavefield points to a standard cell before computing the weights
          case ( 1 , 2 , 3 , 4 ) ! intw_type
             ! first transform wavefield points to standard inversion grid cell
             if(allocated(x)) deallocate(x)
             if(allocated(y)) deallocate(y)
             if(allocated(z)) deallocate(z)
             if(associated(jacobian)) deallocate(jacobian)
             allocate(x(size_idx),y(size_idx),z(size_idx),jacobian(size_idx))
             x = c1(idx); y = c2(idx); z = c3(idx)
             type_standard_cell = 0
             call new(errmsg2,myname)
             call transformToStandardCellInversionGrid(invgrid,icell,x,y,z,uf_wp,jacobian,type_standard_cell,errmsg2)
             if(.level.errmsg2/=0) call print(errmsg2)
             if(.level.errmsg2==2) then
                write(errstr,*) "error transforming wavefield points contained in cell ",icell," to standard cell"
                call add(errmsg,2,trim(errstr),myname)
                goto 1
             endif
             ! assume all arrays to be associated and have same length on exit, if no error
!
             ! then compute weights of specific type
             select case(intw_type)
             case ( 1 )
                call new(errmsg2,myname)
                weight => computeWeightLinearSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                this%err_level_weights(icell) = .level.errmsg2
                this%type_weights(icell) = 1
                call addTrace(errmsg2,myname)
                if(icell==1) then
                    ! the first time, print error message
                   call print(errmsg2)
                else
                   ! all other times print errmsg2 only if there is a warning or error
                   if(.level.errmsg2 /= 0) call print(errmsg2)
                end if
                if(.level.errmsg2 == 2) then
                   write(errstr,*) "error computing type 1 weights of inversion grid cell ",icell
                   call add(errmsg,2,errstr,myname)
                   goto 1
                endif
                call dealloc(errmsg2)
             case ( 2 )
                call new(errmsg2,myname)
                weight => computeWeightQuadraticSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                this%err_level_weights(icell) = .level.errmsg2
                this%type_weights(icell) = 2
                call addTrace(errmsg2,myname)
                if(icell==1) then ! the first time, print error message
                   call print(errmsg2)
                else ! all other times print errmsg2 only if there is a warning or error
                   if(.level.errmsg2 /= 0) call print(errmsg2)
                endif
                if(.level.errmsg2 == 2) then
                   write(errstr,*) "error computing type 2 weights of inversion grid cell ",icell
                   call add(errmsg,2,errstr,myname)
                   goto 1
                endif
                call dealloc(errmsg2)
             case ( 3 )
                call new(errmsg2,myname)
                weight => computeWeightCubicSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                this%err_level_weights(icell) = .level.errmsg2
                this%type_weights(icell) = 3
                call addTrace(errmsg2,myname)
                if(icell==1) then ! the first time, print error message
                   call print(errmsg2)
                else ! all other times print errmsg2 only if there is a warning or error
                   if(.level.errmsg2 /= 0) call print(errmsg2)
                endif
                if(.level.errmsg2 == 2) then
                   write(errstr,*) "error computing type 3 weights of inversion grid cell ",icell
                   call add(errmsg,2,errstr,myname)
                   goto 1
                endif
                call dealloc(errmsg2)
             case ( 4 )
                ! first try cubic SDI (type 3)
                call new(errmsg2,myname)
                weight => computeWeightCubicSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                if(.level.errmsg2 == 0) then
                   this%type_weights(icell) = 3
                else
                   ! if not worked out, try quadratic SDI (type 2)
                   if(associated(weight)) deallocate(weight)
                   call dealloc(errmsg2)
                   call new(errmsg2,myname)
                   weight => computeWeightQuadraticSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                   if(.level.errmsg2 == 0) then
                      this%type_weights(icell) = 2
                   else
                      ! if still not worked out, try linear SDI (type 1)
                      if(associated(weight)) deallocate(weight)
                      call dealloc(errmsg2)
                      call new(errmsg2,myname)
                      weight => computeWeightLinearSDI(type_standard_cell,size_idx,x,y,z,jacobian,errmsg2)
                      if(.level.errmsg2 == 0) then
                         this%type_weights(icell) = 1
                      else
                         ! if not worked out, set to linear weights (constant weights = volum/number_of_points )
                         if(associated(weight)) deallocate(weight)
                         call dealloc(errmsg2)
                         call new(errmsg2,myname)
                         call getVolumeCellInversionGrid(invgrid,icell,volume,errmsg2)
                         if(.level.errmsg2==2) then
                            write(errstr,*) "error getting volume of inversion grid cell ",icell
                            call add(errmsg,2,trim(errstr),myname)
                            goto 1
                         endif
                         weight => computeWeightLinearIntegration(volume,size_idx)
                         this%type_weights(icell) = 5
                      endif
                   endif
                endif
                this%err_level_weights(icell) = .level.errmsg2
                if(icell==1) then ! the first time, print error message
                   call print(errmsg2)
                else
                   if(.level.errmsg2 /= 0) call print(errmsg2)
                end if
                if(.level.errmsg2 == 2) then
                   write(errstr,*) "error computing type 4 weights of inversion grid cell ",icell
                   call add(errmsg,2,trim(errstr),myname)
                   goto 1
                endif
                call dealloc(errmsg2)
             end select
!
          ! for type 5 integration weights, we dont need to transform wavefield points to standard 
          ! cell, we just need the volume of the cell
          case ( 5 ) ! intw_type
             call new(errmsg2,myname)
             call getVolumeCellInversionGrid(invgrid,icell,volume,errmsg2)
             if(.level.errmsg2/=0) call print(errmsg2)
             if(.level.errmsg2==2) then
                write(errstr,*) "error getting volume of inversion grid cell ",icell
                call add(errmsg,2,trim(errstr),myname)
                goto 1
             endif
             call dealloc(errmsg2)
             weight => computeWeightLinearIntegration(volume,size_idx)
             this%err_level_weights(icell) = 0
             this%type_weights(icell) = 5
!
          ! for type 6 integration weights, call transformToStandardCell routine, but use return values
          ! of jacobian as weights directly
          ! this may be used for external integration weights, handled by some specific inversion grid module
          ! (e.g. specfem3dInversionGrid)
          case ( 6 ) ! intw_type
             ! first transform wavefield points to standard inversion grid cell, but indicating 
             ! by type_standard_cell = -1 that on exit jacobian should contain the total weights already
             if(allocated(x)) deallocate(x)
             if(allocated(y)) deallocate(y)
             if(allocated(z)) deallocate(z)
             if(associated(jacobian)) deallocate(jacobian)
             allocate(x(size_idx),y(size_idx),z(size_idx),jacobian(size_idx))
             x = c1(idx); y = c2(idx); z = c3(idx)
             ! by indicating type_standard_cell = -1 on enter, jacobian is supposed to contain the total weights on return and not only the jacobian values
             type_standard_cell = -1
             call new(errmsg2,myname)
             call transformToStandardCellInversionGrid(invgrid,icell,x,y,z,uf_wp,jacobian,type_standard_cell,errmsg2)
             if(.level.errmsg2/=0) call print(errmsg2)
             if(.level.errmsg2==2) then
                write(errstr,*) "error transforming wavefield points contained in cell ",icell," to standard cell"
                call add(errmsg,2,trim(errstr),myname)
                goto 1
             endif
             weight => jacobian
             nullify(jacobian)
             this%err_level_weights(icell) = 0
             this%type_weights(icell) = 6
!
          case default
             write(errstr,*) 'TYPE_INTEGRATION_WEIGHTS = ',intw_type,' not supported'
             call add(errmsg,2,errstr,myname)
             goto 1
          end select ! intw_type
!
          sum_intw = sum_intw + sum(dble(weight))
          ! finally store conputed weights
          call associateVectorPointer(this%weight(icell),weight)
          nullify(weight)
!
       end if ! this%empty_cell(icell)
!
    end do ! icell
!
    write(errstr,*) cnt_empty_cells," inversion grid cells ( out of ",ncell_invgrid,&
         ") do not contain any wavefield points"
    call add(errmsg,0,errstr,myname)
    if(count(this%err_level_weights/=0 .and. .not. this%empty_cell) > 0 ) then
       write(errstr,*) 'for ',count(this%err_level_weights/=0 .and. .not. this%empty_cell),' out of ',&
            ncell_invgrid-cnt_empty_cells,&
            ' inversion grid cells, the computation of the integration weights was not successful (error or warning)'
       call add(errmsg,1,errstr,myname)
    else
       write(errstr,*) 'no errors or warnings occurred computing the actual integration weights for ',&
            ncell_invgrid-cnt_empty_cells," inversion grid cells"
       call add(errmsg,0,errstr,myname)
    endif
    write(errstr,*) 'Total sum of all weights in ',ncell_invgrid-cnt_empty_cells,' invgrid cells = ',sum_intw
    select case (intw_type)
    case (0)
       errstr = trim(errstr)//'; should equal the total number of non-empty inversion grid cells'
    case (1,2,3,4,5,6)
       errstr = trim(errstr)//'; should equal the total volume of all non-empty inversion grid cells'
    end select
    call add(errmsg,0,errstr,myname)
!
    ! clean up
2   if(associated(c1)) deallocate(c1)
    if(associated(c2)) deallocate(c2)
    if(associated(c3)) deallocate(c3)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)
    if(associated(jacobian)) deallocate(jacobian)
!
    return
!
    ! in case there was an error constructing this object, deallocate everything before returning
1   call deallocateIntegrationWeights(this)
    goto 2
  end subroutine createIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief weights for average computation (no integration)
!! \details All weights are the same, namely 1/number_of_points_in_box. 
!! I.e. no integration is done, just build average
!! \param np number of wavefield point in this invgrid cell
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightAverage(np) result(weight)
    integer :: np
    real, dimension(:), pointer :: weight
    allocate(weight(np))
    weight(:) = 1./np ! just build the avarage value (no integration!)
  end function computeWeightAverage
!------------------------------------------------------------------------------------
!> \brief weights for linear integration
!! \details All weights are the same, building the average of the function values 
!!  multiplied by the volume of the inversion grid cell (i.e. linear integration)
!! \param np number of wavefield point in this invgrid cell
!! \param volume volume of inversion grid cell
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightLinearIntegration(volume,np) result(weight)
    ! incoming
    integer :: np
    real :: volume
    ! outgoing
    real, dimension(:), pointer :: weight
!
    allocate(weight(np))
    weight(:) = volume/np ! avarage function value multiplied by cell volume
  end function computeWeightLinearIntegration
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 1
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 1 taken here
!!  Confer document integrationWeights.pdf
!! \param type_standard_cell defines the shape of the standard cell, select specific routine dependent on type (4=Tetrahedron,6=Hexahedron)
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightLinearSDI(type_standard_cell,np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: type_standard_cell,np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local
    character (len=22) :: myname = 'computeWeightLinearSDI'
    character (len=400) :: errstr
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    select case (type_standard_cell) ! (4=Tetrahedron,6=Hexahedron)
    case( 4 )
       weight => computeWeightLinearSDITet(np,xi,eta,zeta,jacobian,errmsg)
    case( 6 )
       weight => computeWeightLinearSDIHex(np,xi,eta,zeta,jacobian,errmsg)
    case default
       write(errstr,*) "type of standard inversion grid cell ",type_standard_cell," is not supported"
       call add(errmsg,2,errstr,myname)
    end select
  end function computeWeightLinearSDI
!------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 1, hexahedral cell [-1,1]^3
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 1 taken here
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightLinearSDIHex(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local, scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: h,cxc,cyc,czc
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(4) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcube
    real, dimension(np,4) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(4,4) :: A
    integer, dimension(4) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=25) :: myname = 'computeWeightLinearSDIHex'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real :: cx0,cx1,cy0,cy1,cz0,cz1
    real :: sumw_std_cell
    real, parameter :: third = 0.333333333
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell [-1,1]^3) will be subdivided into nh*nh*nh subcubes
!
    ! OLD METHOD, I THINK MISLEADING, NUMERICALLY NOT STABLE
    ! np = 0,..,15 => nh = 1
    ! np = 16,..,215 => nh = 2
    ! np = 216,..,511 => nh = 3
    ! np = 512,..,999 => nh = 4
    !
    ! BETTER: 
    ! nh = max(floor(0.5*np**third),1)
    ! this tries to assure, that there are at least 4 (or maximum available) points within one subcube, 
    ! assuming a uniform distribution of points within cell
    ! (otherwise the damping by eta_SDI might be so strong, that the 4x4 matrix A is close to singular)
    nh = max(floor((np/4.)**third),1)
    h = 2./nh ! 2., as we use [-1,1]^3 as the 'standard' invgrid cell in which we do the integration
!
    write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh,", since np = ",np
    call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 4 points, as then
    ! matrix E_SDI cannot have full rank J=4, which is assumed in the theory
    ! a violation to that could most probably result in matrices A being close to singular, and
    ! the integration weights being erroneous
    if(np < 4) then
       write(errstr,*) 'number of points in this cell is smaller than 4: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
!
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,4)),U(1,1),VT(1,1),E_SDI_COPY(np,4))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 4, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 4, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,4))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need these variables inside loop to check if every A is regular
    allocate(S(4),U(1,1),VT(1,1),A_COPY(4,4))
!
    ! loop on all subcubes of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
             ! calculate boundary coordinates of subcube cx0,cx1,cy0,cy1,cz0,cz1
             cx0 = -1. + (ixsubcell -1)*h
             cx1 = cx0 + h
             cy0 = -1. + (iysubcell -1)*h
             cy1 = cy0 + h
             cz0 = -1. + (izsubcell -1)*h
             cz1 = cz0 + h
             ! calculate subcube center cxc,cyc,czc
             cxc = 0.5*(cx0 + cx1)
             cyc = 0.5*(cy0 + cy1)
             czc = 0.5*(cz0 + cz1)
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcube
             ! the basis of the space of all polynomials of degree lower than 1 which is used here is:
             ! {1 , x , y , z}
             c_SDI(1) = (cx1-cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "1"
             c_SDI(2) = 0.5*(cx1*cx1-cx0*cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "x"
             c_SDI(3) = 0.5*(cx1-cx0)*(cy1*cy1-cy0*cy0)*(cz1-cz0) ! integral of "y"
             c_SDI(4) = 0.5*(cx1-cx0)*(cy1-cy0)*(cz1*cz1-cz0*cz0) ! integral of "z"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,4
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',4, 4, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,4)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                ! do workspace query for SSYSV
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',4, 1,    A,  4, IPIV, c_SDI, 4, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 4 , 4, A_COPY, 4, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 4 , 4, A_COPY, 4, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(4) / (S(1)*4)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: call SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',4, 1,    A,  4, IPIV, c_SDI, 4, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightLinearSDIHex: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,4
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcube being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell [-1,1]^3 should equal its volume ( = 8. )
    if(abs(sumw_std_cell-8.)>0.08) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 8.)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-8.)>0.08
       if(abs(sumw_std_cell-8.)>0.004) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 8.)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-8.)>0.004
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 8.)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-8.)>0.004
    endif ! abs(sumw_std_cell-8.)>0.08
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightLinearSDIHex
!------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 1, tetrahedral cell {0 <= x,y,z ; x+y+z<=1}
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 1 taken here
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightLinearSDITet(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local, scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: cxc,cyc,czc,h
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(4) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcell
    real, dimension(np,4) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(4,4) :: A
    integer, dimension(4) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=25) :: myname = 'computeWeightLinearSDITet'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real :: sumw_std_cell
    real, parameter :: sixth = 0.16666666666666
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell {0 <= x,y,z ; x+y+z<=1}) will be subdivided into nh*nh*nh subcells
!
    ! for now, use nh = 1 , as a subdivision of the standard cell into subcells is much more complicated than in the hex case!!
    nh = 1
    h = 1.
    ! as the cell center will be c = (0.25,0.25,0.25), use the following value for h: twice the mean value of the 4 distances of c to all cell corners
    ! MAKES NO DIFFERENCE... (same results, as using 1.)
    !h = 2 * 0.25 * ( sqrt(0.25**2 + 0.25**2 + 0.25**2) + 3*sqrt(0.75**2 + 0.25**2 + 0.25**2) ) != 1.46 
!
    !write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh
    !call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 4 points, as then
    ! matrix E_SDI cannot have full rank J=4, which is assumed in the theory
    ! a violation to that could most probably result in matrices A being close to singular, and
    ! the integration weights being erroneous
    if(np < 4) then
       write(errstr,*) 'number of points in this cell is smaller than 4: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
!
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,4)),U(1,1),VT(1,1),E_SDI_COPY(np,4))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 4, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 4, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,4))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need these variables inside loop to check if every A is regular
    allocate(S(4),U(1,1),VT(1,1),A_COPY(4,4))
!
    ! loop on all subcells of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             ! calculate subcell center cxc,cyc,czc
             cxc = 0.25
             cyc = 0.25
             czc = 0.25
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcell
             ! the basis of the space of all polynomials of degree lower than 1 which is used here is:
             ! {1 , x , y , z}
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             c_SDI(1) = sixth ! integral of "1"
             c_SDI(2) = 1./24. ! integral of "x"
             c_SDI(3) = c_SDI(2) ! integral of "y"
             c_SDI(4) = c_SDI(2) ! integral of "z"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,4
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',4, 4, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,4)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                ! do workspace query for SSYSV
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',4, 1,    A,  4, IPIV, c_SDI, 4, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 4 , 4, A_COPY, 4, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 4 , 4, A_COPY, 4, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(4) / (S(1)*4)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: call SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',4, 1,    A,  4, IPIV, c_SDI, 4, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightLinearSDITet: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,4
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcell being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell {0 <= x,y,z ; x+y+z<=1} should equal its volume ( = 1/6 )
    if(abs(sumw_std_cell-sixth)>(sixth/100.)) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 1/6)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-sixth)>(sixth/100.)
       if(abs(sumw_std_cell-sixth)>(sixth/2000.)) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 1/6)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-sixth)>(sixth/2000.)
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 1/6)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-sixth)>(sixth/2000.)
    endif ! abs(sumw_std_cell-sixth)>(sixth/100.)
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightLinearSDITet
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 2
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 2 taken here (consistent with approximation order 3)
!!  Confer document integrationWeights.pdf
!! \param type_standard_cell defines the shape of the standard cell, select specific routine dependent on type (4=Tetrahedron,6=Hexahedron)
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightQuadraticSDI(type_standard_cell,np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: type_standard_cell,np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local
    character (len=25) :: myname = 'computeWeightQuadraticSDI'
    character (len=400) :: errstr
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    select case (type_standard_cell) ! (4=Tetrahedron,6=Hexahedron)
    case( 4 )
       weight => computeWeightQuadraticSDITet(np,xi,eta,zeta,jacobian,errmsg)
    case( 6 )
       weight => computeWeightQuadraticSDIHex(np,xi,eta,zeta,jacobian,errmsg)
    case default
       write(errstr,*) "type of standard inversion grid cell ",type_standard_cell," is not supported"
       call add(errmsg,2,errstr,myname)
    end select
  end function computeWeightQuadraticSDI
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 2, hexahedral cell [-1,1]^3
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 2 taken here (consistent with approximation order 3)
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightQuadraticSDIHex(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local, scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: h,cxc,cyc,czc
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(10) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcube
    real, dimension(np,10) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(10,10) :: A
    integer, dimension(10) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=28) :: myname = 'computeWeightQuadraticSDIHex'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real :: cx0,cx0cx0,cx1,cx1cx1,cy0,cy0cy0,cy1,cy1cy1,cz0,cz0cz0,cz1,cz1cz1
    real, parameter :: third = 0.333333333
    real :: sumw_std_cell
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell [-1,1]^3) will be subdivided into nh*nh*nh subcubes
!
    ! OLD METHOD, I THINK MISLEADING, NUMERICALLY NOT STABLE
    ! np = 0,..,15 => nh = 1
    ! np = 16,..,215 => nh = 2
    ! np = 216,..,511 => nh = 3
    ! np = 512,..,999 => nh = 4
    !
    ! BETTER: 
    ! nh = max(floor(0.5*np**third),1)
    ! this tries to assure, that there are at least 10 (or maximum available) points within one subcube, 
    ! assuming a uniform distribution of points within cell
    ! (otherwise the damping by eta_SDI might be so strong, that the 10x10 matrix A is close to singular)
    nh = max(floor((np/10.)**third),1)
    h = 2./nh ! 2., as we use [-1,1]^3 as the 'standard' invgrid cell in which we do the integration
!
    write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh
    call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 10 points, as then
    ! matrix E_SDI cannot have full rank J=10, which is assumed in the theory
    ! a violation to that could most probably result in matrices A being close to singular, and
    ! the integration weights being erroneous
    if(np < 10) then
       write(errstr,*) 'number of points in this cell is smaller than 10: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
    E_SDI(:,5) = xi*xi
    E_SDI(:,6) = xi*eta
    E_SDI(:,7) = xi*zeta
    E_SDI(:,8) = eta*eta
    E_SDI(:,9) = eta*zeta
    E_SDI(:,10) = zeta*zeta
!
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,10)),U(1,1),VT(1,1),E_SDI_COPY(np,10))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 10, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 10, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,10))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need these variables inside loop to check if every A is regular
    allocate(S(10),U(1,1),VT(1,1),A_COPY(10,10))
!
    ! loop on all subcubes of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
             ! calculate boundary coordinates of subcube cx0,cx1,cy0,cy1,cz0,cz1
             cx0 = -1. + (ixsubcell -1)*h
             cx1 = cx0 + h
             cy0 = -1. + (iysubcell -1)*h
             cy1 = cy0 + h
             cz0 = -1. + (izsubcell -1)*h
             cz1 = cz0 + h
             ! calculate subcube center cxc,cyc,czc
             cxc = 0.5*(cx0 + cx1)
             cyc = 0.5*(cy0 + cy1)
             czc = 0.5*(cz0 + cz1)
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcube
             ! the basis of the space of all polynomials of degree lower than 2 which is used here is:
             ! {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2}
             ! calculate squared values cx0^2,... in advance to reduce numbers of flops
             cx0cx0 = cx0*cx0
             cx1cx1 = cx1*cx1
             cy0cy0 = cy0*cy0
             cy1cy1 = cy1*cy1
             cz0cz0 = cz0*cz0
             cz1cz1 = cz1*cz1
             c_SDI(1) = (cx1-cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "1"
             c_SDI(2) = 0.5*(cx1cx1-cx0cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "x"
             c_SDI(3) = 0.5*(cx1-cx0)*(cy1cy1-cy0cy0)*(cz1-cz0) ! integral of "y"
             c_SDI(4) = 0.5*(cx1-cx0)*(cy1-cy0)*(cz1cz1-cz0cz0) ! integral of "z"
             c_SDI(5) = third*(cx1cx1*cx1-cx0cx0*cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "x^2"
             c_SDI(6) = 0.25*(cx1cx1-cx0cx0)*(cy1cy1-cy0cy0)*(cz1-cz0) ! integral of "xy"
             c_SDI(7) = 0.25*(cx1cx1-cx0cx0)*(cy1-cy0)*(cz1cz1-cz0cz0) ! integral of "xz"
             c_SDI(8) = third*(cx1-cx0)*(cy1cy1*cy1-cy0cy0*cy0)*(cz1-cz0) ! integral of "y^2"
             c_SDI(9) = 0.25*(cx1-cx0)*(cy1cy1-cy0cy0)*(cz1cz1-cz0cz0) ! integral of "yz"
             c_SDI(10) = third*(cx1-cx0)*(cy1-cy0)*(cz1cz1*cz1-cz0cz0*cz0) ! integral of "z^2"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,10
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',10, 10, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,10)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                ! do workspace query for SSYSV
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',10, 1,    A,  10, IPIV, c_SDI, 10, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 10 , 10, A_COPY, 10, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 10 , 10, A_COPY, 10, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(10) / (S(1)*10)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: call SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',10, 1,    A,  10, IPIV, c_SDI, 10, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightQuadraticSDI: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,10
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcube being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell [-1,1]^3 should equal its volume ( = 8. )
    if(abs(sumw_std_cell-8.)>0.08) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 8.)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-8.)>0.08
       if(abs(sumw_std_cell-8.)>0.004) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 8.)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-8.)>0.004
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 8.)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-8.)>0.004
    endif ! abs(sumw_std_cell-8.)>0.08
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightQuadraticSDIHex
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 2, tetrahedral cell {0 <= x,y,z ; x+y+z<=1}
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 2 taken here (consistent with approximation order 3)
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightQuadraticSDITet(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local, scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: h,cxc,cyc,czc
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(10) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcell
    real, dimension(np,10) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(10,10) :: A
    integer, dimension(10) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=28) :: myname = 'computeWeightQuadraticSDITet'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real :: sumw_std_cell
    real, parameter :: sixth = 0.16666666666666
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell {0 <= x,y,z ; x+y+z<=1}) will be subdivided into nh*nh*nh subcells
!
    ! for now, use nh = 1 , as a subdivision of the standard cell into subcells is much more complicated than in the hex case!!
    nh = 1
    h = 1.
    ! as the cell center will be c = (0.25,0.25,0.25), use the following value for h: twice the mean value of the 4 distances of c to all cell corners
    ! MAKES NO DIFFERENCE... (same results, as using 1.)
    !h = 2 * 0.25 * ( sqrt(0.25**2 + 0.25**2 + 0.25**2) + 3*sqrt(0.75**2 + 0.25**2 + 0.25**2) ) != 1.46 
!
    !write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh
    !call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 10 points, as then
    ! matrix E_SDI cannot have full rank J=10, which is assumed in the theory
    ! a violation to that could most probably result in matrices A being close to singular, and
    ! the integration weights being erroneous
    if(np < 10) then
       write(errstr,*) 'number of points in this cell is smaller than 10: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
    E_SDI(:,5) = xi*xi
    E_SDI(:,6) = xi*eta
    E_SDI(:,7) = xi*zeta
    E_SDI(:,8) = eta*eta
    E_SDI(:,9) = eta*zeta
    E_SDI(:,10) = zeta*zeta
!
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,10)),U(1,1),VT(1,1),E_SDI_COPY(np,10))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 10, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 10, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,10))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need these variables inside loop to check if every A is regular
    allocate(S(10),U(1,1),VT(1,1),A_COPY(10,10))
!
    ! loop on all subcells of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             ! calculate subcell center cxc,cyc,czc
             cxc = 0.25
             cyc = 0.25
             czc = 0.25
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcell
             ! the basis of the space of all polynomials of degree lower than 2 which is used here is:
             ! {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2}
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             c_SDI(1) = sixth     ! integral of "1"
             c_SDI(2) = 1./24.    ! integral of "x"
             c_SDI(3) = c_SDI(2)  ! integral of "y"
             c_SDI(4) = c_SDI(2)  ! integral of "z"
             c_SDI(5) = 1./60.    ! integral of "x^2"
             c_SDI(6) = 1./120.   ! integral of "xy"
             c_SDI(7) = c_SDI(6)  ! integral of "xz"
             c_SDI(8) = c_SDI(5)  ! integral of "y^2"
             c_SDI(9) = c_SDI(6)  ! integral of "yz"
             c_SDI(10) = c_SDI(5) ! integral of "z^2"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,10
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',10, 10, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,10)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                ! do workspace query for SSYSV
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',10, 1,    A,  10, IPIV, c_SDI, 10, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 10 , 10, A_COPY, 10, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 10 , 10, A_COPY, 10, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(10) / (S(1)*10)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: call SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',10, 1,    A,  10, IPIV, c_SDI, 10, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightQuadraticSDI: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,10
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcell being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell {0 <= x,y,z ; x+y+z<=1} should equal its volume ( = 1/6 )
    if(abs(sumw_std_cell-sixth)>(sixth/100.)) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 1/6)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-sixth)>(sixth/100.)
       if(abs(sumw_std_cell-sixth)>(sixth/2000.)) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 1/6)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-sixth)>(sixth/2000.)
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 1/6)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-sixth)>(sixth/2000.)
    endif ! abs(sumw_std_cell-sixth)>(sixth/100.)
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightQuadraticSDITet
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 3
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 1 taken here (consistent with approximation order 4)
!!  Confer document integrationWeights.pdf
!! \param type_standard_cell defines the shape of the standard cell, select specific routine dependent on type (4=Tetrahedron,6=Hexahedron)
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightCubicSDI(type_standard_cell,np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: type_standard_cell,np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local
    character (len=21) :: myname = 'computeWeightCubicSDI'
    character (len=400) :: errstr
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    select case (type_standard_cell) ! (4=Tetrahedron,6=Hexahedron)
    case( 4 )
       weight => computeWeightCubicSDITet(np,xi,eta,zeta,jacobian,errmsg)
    case( 6 )
       weight => computeWeightCubicSDIHex(np,xi,eta,zeta,jacobian,errmsg)
    case default
       write(errstr,*) "type of standard inversion grid cell ",type_standard_cell," is not supported"
       call add(errmsg,2,errstr,myname)
    end select
  end function computeWeightCubicSDI
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 3, hexahedral cell [-1,1]^3
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 3 taken here (consistent with approximation order 4)
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightCubicSDIHex(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: h,cxc,cyc,czc
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(20) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcube
    real, dimension(np,20) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(20,20) :: A
    integer, dimension(20) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=24) :: myname = 'computeWeightCubicSDIHex'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real :: cx0,cx0cx0,cx0cx0cx0,cx1,cx1cx1,cx1cx1cx1, &
         cy0,cy0cy0,cy0cy0cy0,cy1,cy1cy1,cy1cy1cy1, &
         cz0,cz0cz0,cz0cz0cz0,cz1,cz1cz1,cz1cz1cz1
    real, parameter :: third = 0.3333333333333
    real, parameter :: sixth = 0.1666666666666
    real :: sumw_std_cell
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell [-1,1]^3) will be subdivided into nh*nh*nh subcubes
!
    ! OLD METHOD, I THINK MISLEADING, NUMERICALLY NOT STABLE
    ! np = 0,..,15 => nh = 1
    ! np = 16,..,215 => nh = 2
    ! np = 216,..,511 => nh = 3
    ! np = 512,..,999 => nh = 4
    !
    ! BETTER: 
    ! nh = max(floor((np/20.)**third),1)
    ! this tries to assure, that there are at least 20 (or maximum available) points within one subcube 
    ! (otherwise the damping by eta_SDI might be so strong, that the 20x20 matrix A is close to singular)
    nh = max(floor((np/20.)**third),1)
    h = 2./nh ! 2., as we use [-1,1]^3 as the 'standard' invgrid cell in which we do the integration
!
    write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh
    call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 20 points, as then
    ! matrix E_SDI cannot have full rank J=20, which is assumed in the theory
    ! a violation to that could most probably result in matrix A being close to singular, and
    ! the integration weights being erroneous
    if(np < 20) then
       write(errstr,*) 'number of points in this cell is smaller than 20: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2 ..
    ! ..x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
    E_SDI(:,5) = xi*xi
    E_SDI(:,6) = xi*eta
    E_SDI(:,7) = xi*zeta
    E_SDI(:,8) = eta*eta
    E_SDI(:,9) = eta*zeta
    E_SDI(:,10) = zeta*zeta
    E_SDI(:,11) = E_SDI(:,5)*xi
    E_SDI(:,12) = E_SDI(:,5)*eta
    E_SDI(:,13) = E_SDI(:,5)*zeta
    E_SDI(:,14) = E_SDI(:,6)*eta
    E_SDI(:,15) = E_SDI(:,6)*zeta
    E_SDI(:,16) = E_SDI(:,7)*zeta
    E_SDI(:,17) = E_SDI(:,8)*eta
    E_SDI(:,18) = E_SDI(:,8)*zeta
    E_SDI(:,19) = E_SDI(:,9)*zeta
    E_SDI(:,20) = E_SDI(:,10)*zeta
!	
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,20)),U(1,1),VT(1,1),E_SDI_COPY(np,20))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 20, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 20, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,20))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need those variables inside loop to check if matrices A are regular
    allocate(S(20),U(1,1),VT(1,1),A_COPY(20,20))
!
    ! loop on all subcubes of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
             ! calculate boundary coordinates of subcube cx0,cx1,cy0,cy1,cz0,cz1
             cx0 = -1. + (ixsubcell -1)*h
             cx1 = cx0 + h
             cy0 = -1. + (iysubcell -1)*h
             cy1 = cy0 + h
             cz0 = -1. + (izsubcell -1)*h
             cz1 = cz0 + h
             ! calculate subcube center cxc,cyc,czc
             cxc = 0.5*(cx0 + cx1)
             cyc = 0.5*(cy0 + cy1)
             czc = 0.5*(cz0 + cz1)
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcube
             ! the basis of the space of all polynomials of degree lower than 3 which is used here is:
             ! {  1,    x,   y,     z, x^2,   xy,  xz,  y^2,   yz, z^2, ..
             !..x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3}
             ! calculate powers cx0^2,cx0^3... in advance to reduce numbers of flops
             cx0cx0 = cx0*cx0; cx0cx0cx0 = cx0cx0*cx0
             cx1cx1 = cx1*cx1; cx1cx1cx1 = cx1cx1*cx1
             cy0cy0 = cy0*cy0; cy0cy0cy0 = cy0cy0*cy0
             cy1cy1 = cy1*cy1; cy1cy1cy1 = cy1cy1*cy1
             cz0cz0 = cz0*cz0; cz0cz0cz0 = cz0cz0*cz0
             cz1cz1 = cz1*cz1; cz1cz1cz1 = cz1cz1*cz1
             c_SDI(1) = (cx1-cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "1"
             c_SDI(2) = 0.5*(cx1cx1-cx0cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "x"
             c_SDI(3) = 0.5*(cx1-cx0)*(cy1cy1-cy0cy0)*(cz1-cz0) ! integral of "y"
             c_SDI(4) = 0.5*(cx1-cx0)*(cy1-cy0)*(cz1cz1-cz0cz0) ! integral of "z"
             c_SDI(5) = third*(cx1cx1cx1-cx0cx0cx0)*(cy1-cy0)*(cz1-cz0) ! integral of "x^2"
             c_SDI(6) = 0.25*(cx1cx1-cx0cx0)*(cy1cy1-cy0cy0)*(cz1-cz0) ! integral of "xy"
             c_SDI(7) = 0.25*(cx1cx1-cx0cx0)*(cy1-cy0)*(cz1cz1-cz0cz0) ! integral of "xz"
             c_SDI(8) = third*(cx1-cx0)*(cy1cy1cy1-cy0cy0cy0)*(cz1-cz0) ! integral of "y^2"
             c_SDI(9) = 0.25*(cx1-cx0)*(cy1cy1-cy0cy0)*(cz1cz1-cz0cz0) ! integral of "yz"
             c_SDI(10) = third*(cx1-cx0)*(cy1-cy0)*(cz1cz1cz1-cz0cz0cz0) ! integral of "z^2"
             c_SDI(11) = 0.25*(cx1cx1cx1*cx1-cx0cx0cx0*cx0)*(cy1-cy0)*(cz1-cz0)! integral of "x^3"
             c_SDI(12) = sixth*(cx1cx1cx1-cx0cx0cx0)*(cy1cy1-cy0cy0)*(cz1-cz0) ! integral of "x^2y"
             c_SDI(13) = sixth*(cx1cx1cx1-cx0cx0cx0)*(cy1-cy0)*(cz1cz1-cz0cz0) ! integral of "x^2z"
             c_SDI(14) = sixth*(cx1cx1-cx0cx0)*(cy1cy1cy1-cy0cy0cy0)*(cz1-cz0) ! integral of "xy^2"
             c_SDI(15) = 0.125*(cx1cx1-cx0cx0)*(cy1cy1-cy0cy0)*(cz1cz1-cz0cz0) ! integral of "xyz"
             c_SDI(16) = sixth*(cx1cx1-cx0cx0)*(cy1-cy0)*(cz1cz1cz1-cz0cz0cz0) ! integral of "xz^2"
             c_SDI(17) = 0.25*(cx1-cx0)*(cy1cy1cy1*cy1-cy0cy0cy0*cy0)*(cz1-cz0) ! integral of "y^3"
             c_SDI(18) = sixth*(cx1-cx0)*(cy1cy1cy1-cy0cy0cy0)*(cz1cz1-cz0cz0) ! integral of "y^2z"
             c_SDI(19) = sixth*(cx1-cx0)*(cy1cy1-cy0cy0)*(cz1cz1cz1-cz0cz0cz0) ! integral of "yz^2"
             c_SDI(20) = 0.25*(cx1-cx0)*(cy1-cy0)*(cz1cz1cz1*cz1-cz0cz0cz0*cz0) ! integral of "z^3"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,20
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',10, 10, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,10)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',20, 1,    A,  20, IPIV, c_SDI, 20, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 20 , 20, A_COPY, 20, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 20 , 20, A_COPY, 20, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(20) / (S(1)*20)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',20, 1,    A,  20, IPIV, c_SDI, 20, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightCubicSDI: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,10
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcube ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcube being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell [-1,1]^3 should equal its volume ( = 8. )
!
    !########################################################################################
    ! SHOULD HERE (IN COMPARISON TO computeWeightQuadraticSDI) THRESHOLDS LOWER THAN 1%,0.05% BE USED ???? (-> higher order = higher precision)
    !########################################################################################
    if(abs(sumw_std_cell-8.)>0.08) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 8.)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-8.)>0.08
       if(abs(sumw_std_cell-8.)>0.004) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 8.)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-8.)>0.004
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 8.)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-8.)>0.004
    endif ! abs(sumw_std_cell-8.)>0.08
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightCubicSDIHex
!------------------------------------------------------------------------------------
!> \brief weights for Scattered Data Integration, degree 3, tetrahedral cell {0 <= x,y,z ; x+y+z<=1}
!! \details Scattered Data Integration, as in paper by 
!!  David Levin, Journal of Computational and Applied Mathematics 112 (1999) 181-187
!!  polynomial degree 3 taken here (consistent with approximation order 4)
!!  Confer document integrationWeights.pdf
!! \param np number of incoming points, i.e. size of xi,eta,zeta,jacobian
!! \param xi x coordinates of wavefield points inside standard cell
!! \param eta y coordinates of wavefield points inside standard cell
!! \param zeta z coordinates of wavefield points inside standard cell
!! \param jacobian jacobian of transformation from standard cell to real coordinate cell (to be multiplied to standard weights)
!! \param errmsg error message
!! \param weight vector of integration weights
!! \return integration weights
!
  function computeWeightCubicSDITet(np,xi,eta,zeta,jacobian,errmsg) result(weight)
    ! incoming
    integer :: np
    real, dimension(np) :: xi,eta,zeta,jacobian
    ! outgoing
    real, dimension(:), pointer :: weight
    type (error_message) :: errmsg
    ! local scattered data integration method
    ! nomenclature "_SDI" corresponds to SDI paper by D. Levin
    real :: h,cxc,cyc,czc
    integer :: nh
    real, dimension(np) :: eta_SDI ! actually this is 1/eta, so diagonal entries of D^-1 as in paper by Levin
    real, dimension(20) :: c_SDI ! vector that contains the integrals of the polynomial basis over the current subcell
    real, dimension(np,20) :: E_SDI,DinvE_SDI ! matrices corresponding to "E" and "D^-1*E" (see paper and/or documentation)
    real, dimension(20,20) :: A
    integer, dimension(20) :: IPIV ! output of LAPACK routine
    real, dimension(:),allocatable :: WORK,WORK_SGESVD
    real, dimension(:), allocatable :: S
    real, dimension(:,:), allocatable :: U,VT,E_SDI_COPY,A_COPY
    integer :: INFO,LWORK,LWORK_SGESVD
    real :: sv_ratio_E_SDI,sv_ratio_A,eps
    ! local other
    character (len=400) :: errstr
    character (len=24) :: myname = 'computeWeightCubicSDITet'
    integer :: ixsubcell,iysubcell,izsubcell,ibase,itmp
    real, parameter :: sixth = 0.1666666666666
    real :: sumw_std_cell
!
    call addTrace(errmsg,myname)
    nullify(weight)
!
    ! define nh (and, hence, h). This invgrid cell (or rather the 'standard' cell {0 <= x,y,z ; x+y+z<=1}) will be subdivided into nh*nh*nh subcells
!
    ! for now, use nh = 1 , as a subdivision of the standard cell into subcells is much more complicated than in the hex case!!
    nh = 1
    h = 1.
    ! as the cell center will be c = (0.25,0.25,0.25), use the following value for h: twice the mean value of the 4 distances of c to all cell corners
    ! MAKES NO DIFFERENCE... (same results, as using 1.)
    !h = 2 * 0.25 * ( sqrt(0.25**2 + 0.25**2 + 0.25**2) + 3*sqrt(0.75**2 + 0.25**2 + 0.25**2) ) != 1.46 
!
    !write(errstr,*) 'cell will be subdivided into nh^3 subcells, nh = ',nh
    !call add(errmsg,0,trim(errstr),myname)
!
    ! give warning, if this invgrid cell contains less than 20 points, as then
    ! matrix E_SDI cannot have full rank J=20, which is assumed in the theory
    ! a violation to that could most probably result in matrix A being close to singular, and
    ! the integration weights being erroneous
    if(np < 20) then
       write(errstr,*) 'number of points in this cell is smaller than 20: weights can likely be erroneous!'
       call add(errmsg,1,trim(errstr),myname)
    endif
!
    allocate(weight(np))
    weight(:) = 0.
!
    ! set up matrix E with E_ij = p_j(v_i), where v=(xi,eta,zeta) is the vector of x,y,z coordinates in standard cell
    ! here, p_j is in the polinomial basis {1 , x , y , z , x^2 , xy , xz , y^2 , yz , z^2 ..
    ! ..x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3}
    E_SDI(:,1) = 1.
    E_SDI(:,2) = xi
    E_SDI(:,3) = eta
    E_SDI(:,4) = zeta
    E_SDI(:,5) = xi*xi
    E_SDI(:,6) = xi*eta
    E_SDI(:,7) = xi*zeta
    E_SDI(:,8) = eta*eta
    E_SDI(:,9) = eta*zeta
    E_SDI(:,10) = zeta*zeta
    E_SDI(:,11) = E_SDI(:,5)*xi
    E_SDI(:,12) = E_SDI(:,5)*eta
    E_SDI(:,13) = E_SDI(:,5)*zeta
    E_SDI(:,14) = E_SDI(:,6)*eta
    E_SDI(:,15) = E_SDI(:,6)*zeta
    E_SDI(:,16) = E_SDI(:,7)*zeta
    E_SDI(:,17) = E_SDI(:,8)*eta
    E_SDI(:,18) = E_SDI(:,8)*zeta
    E_SDI(:,19) = E_SDI(:,9)*zeta
    E_SDI(:,20) = E_SDI(:,10)*zeta
!	
    ! estimate if E_SDI has full rank by looking at the ratio of largest and smallest singular value
    allocate(S(min(np,20)),U(1,1),VT(1,1),E_SDI_COPY(np,20))
    E_SDI_COPY = E_SDI
    ! do workspace query first
    allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
    !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
    call SGESVD( 'N', 'N', np , 20, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'workspace query failed: LAPACK routine SGESVD applied to E_SDI returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    LWORK_SGESVD = WORK_SGESVD(1)
    write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to E_SDI is ',LWORK_SGESVD
    call add(errmsg,0,trim(errstr),myname)
    deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
    ! now compute singular values
    call SGESVD( 'N', 'N', np , 20, E_SDI_COPY, np, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
    if(INFO/=0) then
       write(errstr,*) 'computation of singular values of E_SDI failed: SGESVD returned INFO = ',INFO
       call add(errmsg,2,trim(errstr),myname)
       return
    endif ! INFO/=0
    ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
    sv_ratio_E_SDI = S(size(S)) / (S(1)*max(np,20))
    eps = epsilon(1.0)
    deallocate(U,VT,WORK_SGESVD,E_SDI_COPY,S)
    if(sv_ratio_E_SDI < eps) then
       write(errstr,*) 'matrix E_SDI is (close to) singular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' < machine_epsilon = ',eps
       call add(errmsg,1,trim(errstr),myname)
    else
       write(errstr,*) 'matrix E_SDI seems to be regular, as sv_ratio_E_SDI = ',sv_ratio_E_SDI,&
            ' >= machine_epsilon = ',eps
       call add(errmsg,0,trim(errstr),myname)
    endif
!
    ! need those variables inside loop to check if matrices A are regular
    allocate(S(20),U(1,1),VT(1,1),A_COPY(20,20))
!
    ! loop on all subcells of this invgrid cell
    do ixsubcell = 1,nh
       do iysubcell = 1,nh
          do izsubcell = 1,nh		
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             ! calculate subcell center cxc,cyc,czc
             cxc = 0.25
             cyc = 0.25
             czc = 0.25
!
             ! calculate vector containing the integrals of the polynomial basis over the current subcell
             ! the basis of the space of all polynomials of degree lower than 3 which is used here is:
             ! {  1,    x,   y,     z, x^2,   xy,  xz,  y^2,   yz, z^2, ..
             !..x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3}
! ATTENTION: THIS "SUBCELL" IS THE  W H O L E   TETRAHEDRON!!! IF YOU CHOOSE ANY REAL SUBCELLS, ADJUST ACCORDINGLY!!
             c_SDI(1) = sixth     ! integral of "1"
             c_SDI(2) = 1./24.    ! integral of "x"
             c_SDI(3) = c_SDI(2)  ! integral of "y"
             c_SDI(4) = c_SDI(2)  ! integral of "z"
             c_SDI(5) = 1./60.    ! integral of "x^2"
             c_SDI(6) = 1./120.   ! integral of "xy"
             c_SDI(7) = c_SDI(6)  ! integral of "xz"
             c_SDI(8) = c_SDI(5)  ! integral of "y^2"
             c_SDI(9) = c_SDI(6)  ! integral of "yz"
             c_SDI(10) = c_SDI(5) ! integral of "z^2"
             c_SDI(11) = c_SDI(6) ! integral of "x^3"
             c_SDI(12) = 1./360.  ! integral of "x^2y"
             c_SDI(13) = c_SDI(12)! integral of "x^2z"
             c_SDI(14) = c_SDI(12)! integral of "xy^2"
             c_SDI(15) = 1./720.  ! integral of "xyz"
             c_SDI(16) = c_SDI(12)! integral of "xz^2"
             c_SDI(17) = c_SDI(6) ! integral of "y^3"
             c_SDI(18) = c_SDI(12)! integral of "y^2z"
             c_SDI(19) = c_SDI(12)! integral of "yz^2"
             c_SDI(20) = c_SDI(6) ! integral of "z^3"
!
             ! calculate entries of diagonal matrix D^-1
             eta_SDI = 0.5*exp(-((cxc-xi)*(cxc-xi) + (cyc-eta)*(cyc-eta) + (czc-zeta)*(czc-zeta)) / (h*h))
! 
             ! set up matrix D^-1*E with (D^-1*E)_ij = eta(norm(cxc-xi_i))*p_j(xi_i)
             do ibase = 1,20
                DinvE_SDI(:,ibase) = E_SDI(:,ibase)*eta_SDI
             enddo ! ibase
! 
             ! compute matrix A = E^T*D^-1*E of the linear system that is to be solved
             A = matmul(transpose(E_SDI),DinvE_SDI)
             !ALTERNATIVELY: call SGEMM(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) ! on exit, C contains the result
             !call SGEMM('T','N',10, 10, np,1.,E_SDI,np,DinvE_SDI,np, 0., A,10)
!
             A_COPY = A
!
             ! in the very first iteration, carry out workspace query in order to allocate WORK array of optimal size
             if(ixsubcell==1 .and. iysubcell==1 .and. izsubcell==1) then
                allocate(WORK(1)); LWORK = -1
                call SSYSV(  'U',20, 1,    A,  20, IPIV, c_SDI, 20, WORK, LWORK, INFO)
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SSYSV returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK = WORK(1)
                deallocate(WORK); allocate(WORK(LWORK))
!
                ! do workspace query for SGESVD
                allocate(WORK_SGESVD(1)); LWORK_SGESVD = -1
                !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
                call SGESVD( 'N', 'N', 20 , 20, A_COPY, 20, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
                if(INFO/=0) then
                   write(errstr,*) 'workspace query failed; SGESVD applied to A returned INFO = ',INFO
                   call add(errmsg,2,trim(errstr),myname)
                   return
                endif ! INFO/=0
                LWORK_SGESVD = WORK_SGESVD(1)
                deallocate(WORK_SGESVD); allocate(WORK_SGESVD(LWORK_SGESVD))
!
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SSYSV is ',LWORK
                call add(errmsg,0,trim(errstr),myname)
                write(errstr,*) 'optimal size of WORK array for LAPACK routine SGESVD applied to A is ',LWORK_SGESVD
                call add(errmsg,0,trim(errstr),myname)
             endif ! ixsubcell==iysubcell==izsubcell==1
!
             ! estimate if A has full rank by looking at the ratio of largest and smallest singular value
             !SYNTAX: call SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
             call SGESVD( 'N', 'N', 20 , 20, A_COPY, 20, S, U, 1, VT, 1, WORK_SGESVD, LWORK_SGESVD, INFO )
             if(INFO/=0) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SGESVD failed, it returned INFO = ',INFO
                call add(errmsg,2,trim(errstr),myname)
                return
             endif ! INFO/=0
             ! sv_ratio = smallest singular value / ( largest singular value * number of singular values )  (as MATLAB does)
             sv_ratio_A = S(20) / (S(1)*20)
             ! give immediate warning, if for this subcell sv_ratio_A is small (i.e. this A seems to be (close to) singular)
             if(sv_ratio_A < eps) then
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : matrix A is (close to) singular, as sv_ratio_A = ',sv_ratio_A,&
                     ' < machine_ epsilon = ',eps
                call add(errmsg,1,trim(errstr),myname)
             endif ! sv_ratio_A < eps
!
             ! solve linear system 
             !SYNTAX: SSYSV(UPLO, N, NRHS, A, LDA, IPIV, B,    LDB, WORK, LWORK, INFO)
             call SSYSV(  'U',20, 1,    A,  20, IPIV, c_SDI, 20, WORK, LWORK, INFO)
             if(INFO /= 0) then
                ! UNCOMMENT FOR DETAILED DEBUGGING:
                !print *, 'in computeWeightCubicSDI: ixsubcell,iysubcell,izsubcell = ', &
                !   ixsubcell,iysubcell,izsubcell,'; INFO = ',INFO
                !print *, 'np,xtmin,xtmax,ytmin,ytmax,zb,dz,R=',np,xtmin,xtmax,ytmin,ytmax,zb,dz,R
                !print *, 'nh,h=',nh,h
                !print *, '####################  A= ###################################'
                !do itmp = 1,10
                !	print *, A(itmp,:)
                !end do
                !print *, '#####################  DinvE_SDI= ##################################'
                !do itmp = 1,np
                !	print *, DinvE_SDI(itmp,:)
                !end do
                !print *, 'xt,yt,zig,xi,eta,zeta,eta_SDI : '
                !do itmp = 1,np
                !	print *,xt(itmp),yt(itmp),zig(itmp),xi(itmp),eta(itmp),zeta(itmp),eta_SDI(itmp)
                !enddo
                write(errstr,*) 'subcell ix,iy,iz = ',ixsubcell,iysubcell,izsubcell,&
                     ' : SSYSV failed, it returned INFO = ',INFO,&
                     '. If there was no warning about the A of this subcell being (close to) singular, it seems to be regular'
                call add(errmsg,2,trim(errstr),myname)
                return
             endif
!
             ! update weights
             weight = weight + matmul(DinvE_SDI,c_SDI)
             !ALTERNATIVELY: call SGEMV(TRANS, M, N, ALPHA, A LDA,X,INCX,BETA,Y,INCY) ! on exit, Y contains the result
          enddo ! izsubcell
       enddo ! iysubcell
    enddo ! ixsubcell
    deallocate(U,VT,A_COPY,S)
    deallocate(WORK,WORK_SGESVD)
!
    ! CHECK IF COMPUTED WEIGHTS ARE OK
    sumw_std_cell = sum(weight)
    ! the sum of the weights for the standard invgrid cell {0 <= x,y,z ; x+y+z<=1} should equal its volume ( = 1/6 )
!
    !########################################################################################
    ! SHOULD HERE (IN COMPARISON TO computeWeightQuadraticSDI) THRESHOLDS LOWER THAN 1%,0.05% BE USED ???? (-> higher order = higher precision)
    !########################################################################################
    if(abs(sumw_std_cell-sixth)>(sixth/100.)) then
       ! UNCOMMENT FOR DETAILED DEBUGGING:
       !write(*,*) "ERRONEOUS CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
       !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
       !do itmp = 1,np
       !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
       !enddo ! itmp
       write(errstr,*) 'weights may be erroneous: sum(weights_standard_cell) = ',sumw_std_cell,&
            ' differs by more than 1 percent from expected value ( = 1/6)'
       call add(errmsg,1,trim(errstr),myname)
    else ! abs(sumw_std_cell-sixth)>(sixth/100.)
       if(abs(sumw_std_cell-sixth)>(sixth/2000.)) then
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "INEXACT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          write(errstr,*) 'weights may be inexact: sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by between 0.05 and 1 percent from expected value ( = 1/6)'
          call add(errmsg,1,trim(errstr),myname)
       else ! abs(sumw_std_cell-sixth)>(sixth/2000.)
          ! UNCOMMENT FOR DETAILED DEBUGGING:
          !write(*,*) "CORRECT CELL: icell,sum(weight),sum(weight*jacobian),np,sv_max,sv_min = ", &
          !  icell,sum(weight),sum(weight*jacobian),np,S(1),S(size(S))
          !do itmp = 1,np
          !	write(*,*) xi(itmp),eta(itmp),zeta(itmp),weight(itmp),jacobian(itmp)
          !enddo ! itmp
          !call add(errmsg,2,'STOP, found what I looked for',myname)
          !return
          write(errstr,*) 'integration weights look ok, as sum(weights_standard_cell) = ',sumw_std_cell,&
               ' differs by less than 0.05 percent from expected value ( = 1/6)'
          call add(errmsg,0,trim(errstr),myname)
       endif ! abs(sumw_std_cell-sixth)>(sixth/2000.)
    endif ! abs(sumw_std_cell-sixth)>(sixth/100.)
!
    ! account for change of volume caused by transformation from real world to standard inversion grid cell
    weight = weight*jacobian
  end function computeWeightCubicSDITet
!------------------------------------------------------------------------------------
!> \brief Deallocate integration weights
!
  subroutine deallocateIntegrationWeights(this)
    type (integration_weights) :: this
    integer :: j
    if (associated(this%wp_idx)) then
       do j = 1,size(this%wp_idx)
          call deallocIntegerVectorPointer(this%wp_idx(j))
       enddo
    endif
    if (associated(this%weight)) then
       do j = 1,size(this%weight)
          call deallocRealVectorPointer(this%weight(j))
       enddo
    endif
    if (associated(this%empty_cell)) deallocate(this%empty_cell)
    if (associated(this%err_level_weights)) deallocate(this%err_level_weights)
    if (associated(this%type_weights)) deallocate(this%type_weights)
    this%ntot_wp = 0
    this%unit_factor = 0
  end subroutine deallocateIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get total number of inversion grid cells, used to create this integration weights object
!
  function getNtotInvgridIntegrationWeights(this) result(n)
    type (integration_weights), intent(in) :: this
    integer :: n
    if(associated(this%wp_idx)) then
       n = size(this%wp_idx)
    else
       n = 0
    end if
  end function getNtotInvgridIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get total number of inversion grid cells, used to create this integration weights object
!
  function getNtotWpIntegrationWeights(this) result(n)
    type (integration_weights), intent(in) :: this
    integer :: n
    n = this%ntot_wp
  end function getNtotWpIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get indices of wavefield points in given inversion box
!
  function getSelectedWpIdxIntegrationWeights(this,icell) result(idx)
    type (integration_weights), intent(in) :: this
    integer, intent(in) :: icell
    integer, dimension(:), pointer :: idx
!
    nullify(idx)
    if(.not.associated(this%wp_idx)) return
!
    if(icell >= 1 .and. icell <= size(this%wp_idx)) &
         idx => getVectorPointer(this%wp_idx(icell))
  end function getSelectedWpIdxIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get integration weights for wavefield points in given inversion box
!
  function getSelectedIntegrationWeights(this,icell) result(weight)
    type (integration_weights), intent(in) :: this
    integer, intent(in) :: icell
    real, dimension(:), pointer :: weight
!
    nullify(weight)
    if(.not.associated(this%weight)) return
!
    if(icell >= 1 .and. icell <= size(this%weight)) &
         weight => getVectorPointer(this%weight(icell))
  end function getSelectedIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get integration weights for wavefield points in given inversion box
!
  function getAtWavefieldPointsIntegrationWeights(this) result(w)
    type (integration_weights) :: this
    real, dimension(:), pointer :: w
    integer :: i
    real, dimension(:), pointer :: weight
    integer, dimension(:), pointer :: idx
!
    nullify(weight,idx)
    nullify(w)
    if(this%ntot_wp==0 .or. .not.associated(this%wp_idx)) return
!
    allocate(w(this%ntot_wp))
    w(:) = 0.
    do i = 1,size(this%wp_idx)
       if(.not.this%empty_cell(i)) then
          idx => getVectorPointer(this%wp_idx(i))
          weight => getVectorPointer(this%weight(i))
          w(idx) = weight(:)
       endif
    enddo ! i
  end function getAtWavefieldPointsIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get cell index for each wavefield points (if not located, set -1)
!
  function getCellIdexAtWavefieldPointsIntegrationWeights(this) result(icell)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: icell
    integer :: i
    integer, dimension(:), pointer :: idx
!
    nullify(idx)
    nullify(icell)
    if(this%ntot_wp==0 .or. .not.associated(this%wp_idx)) return
!
    allocate(icell(this%ntot_wp))
    icell = -1
    do i = 1,size(this%wp_idx) ! i is inversion grid cell index
       if(.not.this%empty_cell(i)) then
          idx => getVectorPointer(this%wp_idx(i))
          icell(idx) = i
       end if
    end do ! i
  end function getCellIdexAtWavefieldPointsIntegrationWeights
!------------------------------------------------------------------------------------
  function getNwpPerBoxIntegrationWeights(this) result(n)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: n
    integer :: i
    integer, dimension(:), pointer :: idx
!
    nullify(idx)
    nullify(n)
    if(.not.associated(this%wp_idx)) return
!
    allocate(n(size(this%wp_idx)))
    do i = 1,size(this%wp_idx)
       if(this%empty_cell(i)) then
          n(i) = 0
       else
          idx => getVectorPointer(this%wp_idx(i))
          n(i) = size(idx)
       endif
    enddo ! i
  end function getNwpPerBoxIntegrationWeights
!------------------------------------------------------------------------------------
  function getErrLevelIntegrationWeights(this) result(err_level)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: err_level
    err_level => this%err_level_weights
  end function getErrLevelIntegrationWeights
!------------------------------------------------------------------------------------
  function getTypeIntegrationWeights(this) result(typ)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: typ
    typ => this%type_weights
  end function getTypeIntegrationWeights
!------------------------------------------------------------------------------------
  function getIndicesEmptyCellsIntegrationWeights(this) result(indx)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: indx
    integer :: n,i
!
    indx => null()
    if(.not.associated(this%empty_cell)) return
!
    n = count(this%empty_cell)
    if(n>0) then
       allocate(indx(n))
       indx = pack( (/ (i,i=1,size(this%empty_cell)) /), this%empty_cell)
    endif
  end function getIndicesEmptyCellsIntegrationWeights
!------------------------------------------------------------------------------------
  function getIndicesFilledCellsIntegrationWeights(this) result(indx)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: indx
    integer :: n,i
!
    indx => null()
    if(.not.associated(this%empty_cell)) return
!
    n = count(.not.this%empty_cell)
    if(n>0) then
       allocate(indx(n))
       indx = pack( (/ (i,i=1,size(this%empty_cell)) /), .not.this%empty_cell)
    endif
  end function getIndicesFilledCellsIntegrationWeights
!------------------------------------------------------------------------------------
  function cellIsEmptyIntegrationWeights(this,indx) result(b)
    type (integration_weights) :: this
    integer :: indx
    logical :: b
    ! when returning only one logical value, you cannot distinguish cases where indx is out of range. 
    ! for this case, the answer shall be always true (which is the "bad" case, meaning the cell is empty, i.e. should not be used)
    b = .true.
    if(.not.associated(this%wp_idx)) return
    if(indx>=1 .and. indx<=size(this%wp_idx)) then
       b = this%empty_cell(indx)
    endif
  end function cellIsEmptyIntegrationWeights
!------------------------------------------------------------------------------------
  function anyCellIsEmptyIntegrationWeights(this) result(b)
    type (integration_weights) :: this
    logical :: b
    if(.not.associated(this%empty_cell)) then
       b = .true.
    else
       b = any(this%empty_cell)
    end if
  end function anyCellIsEmptyIntegrationWeights
!------------------------------------------------------------------------------------
  function getIndicesWpInsideInvgridIntegrationWeights(this) result(indx_return)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: indx_return
    integer :: icell,n
    integer, dimension(:), pointer :: idx
    integer, dimension(:), allocatable :: indx_inside
!
    nullify(idx)
    nullify(indx_return)
    if(.not.associated(this%wp_idx) .or. this%ntot_wp==0) return
!
    allocate(indx_inside(this%ntot_wp))
    ! first, set no wp to be inside
    indx_inside = 0
    do icell = 1,size(this%wp_idx)
       if(.not.this%empty_cell(icell)) then
          ! if wp inside, set indx_inside to the respective index
          idx => getVectorPointer(this%wp_idx(icell))
          indx_inside(idx) = idx
       endif
    enddo ! icell
!
    n = count(indx_inside .gt. 0)
    if(n>0) then
       allocate(indx_return(n))
       indx_return = pack( indx_inside, indx_inside .gt. 0)
    endif
    deallocate(indx_inside)
  end function getIndicesWpInsideInvgridIntegrationWeights
!------------------------------------------------------------------------------------
  function getIndicesWpOutsideInvgridIntegrationWeights(this) result(indx_return)
    type (integration_weights) :: this
    integer, dimension(:), pointer :: indx_return
    integer :: icell,n
    integer, dimension(:), pointer :: idx
    integer, dimension(:), allocatable :: indx_outside
!
    nullify(idx)
    nullify(indx_return)
    if(.not.associated(this%wp_idx) .or. this%ntot_wp==0) return
!
    allocate(indx_outside(this%ntot_wp))
    ! first, set all wp to be outside
    indx_outside = (/ (n,n=1,this%ntot_wp) /)
    do icell = 1,size(this%wp_idx)
       if(.not.this%empty_cell(icell)) then
          ! if wp inside, set indx_outside to zero
          idx => getVectorPointer(this%wp_idx(icell))
          indx_outside(idx) = 0
       endif
    enddo ! icell
!
    n = count(indx_outside .gt. 0)
    if(n>0) then
       allocate(indx_return(n))
       indx_return = pack( indx_outside, indx_outside .gt. 0)
    endif
    deallocate(indx_outside)
  end function getIndicesWpOutsideInvgridIntegrationWeights
!------------------------------------------------------------------------------------
  function anyWpIsOutsideInvgridIntegrationWeights(this) result(b)
    type (integration_weights) :: this
    logical :: b
    integer :: icell
    integer, dimension(:), pointer :: idx
    logical, dimension(:), allocatable :: indx_outside
!
    nullify(idx)
    b = .true.
    if(this%ntot_wp == 0 .or. .not.associated(this%wp_idx)) return
!
    allocate(indx_outside(this%ntot_wp))
    ! first, set all wp to be outside
    indx_outside = .true.
    do icell = 1,size(this%wp_idx)
       if(.not.this%empty_cell(icell)) then
          ! if wp inside, set indx_outside false
          idx => getVectorPointer(this%wp_idx(icell))
          indx_outside(idx) = .false.
       endif
    enddo ! icell
    b = any(indx_outside)
    deallocate(indx_outside)
  end function anyWpIsOutsideInvgridIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Get unit factor of integration weights
!
  subroutine getUnitFactorIntegrationWeights(this,uf_intw,errmsg)
    type (integration_weights) :: this
    real :: uf_intw
    type (error_message) :: errmsg
    if(this%unit_factor == 0) then
       call add(errmsg,2,"the integration weights object is not yet created; unit factor is still unknown",&
            "getUnitFactorIntegrationWeights")
       uf_intw = 0
    else
       uf_intw = this%unit_factor
    end if
  end subroutine getUnitFactorIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Write integration weights to unformatted file
!
  subroutine writeIntegrationWeights(this,lu,filename,errmsg)
    type (integration_weights) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    character(len=400) :: errstr
    integer :: ierr,j
    integer, dimension(:), pointer :: pidx
    real, dimension(:), pointer :: pw
    character (len=23) :: myname = 'writeIntegrationWeights'
!
    nullify(pidx,pw)
    call addTrace(errmsg,myname)
!
    if(.not.associated(this%wp_idx)) then
       call add(errmsg,2,"integration weights not defined, so nothing to write",myname)
       return
    end if
!
    open(lu,file = trim(filename),form = 'unformatted',iostat = ierr)
    if (ierr > 0) then
       call add(errmsg,2,trim(filename)//' could not be opened',myname)
       return
    endif
!
!  Loop over inversion grid cells
!
    write(lu) size(this%wp_idx)
    do j = 1,size(this%wp_idx)
       pidx => getVectorPointer(this%wp_idx(j))
       if(.not.associated(pidx)) then
          write(errstr,*) "for inversion grid cell ",j,", there are no wavefield point indices present, "//&
               "meaning that there was a problem in locaĺizing the wavefield points inside the inversion grid "//&
               "(even if cell is empty, there should be one artificial point)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       pw => getVectorPointer(this%weight(j))
       if(.not.associated(pw)) then
          write(errstr,*) "for inversion grid cell ",j,", there are no integration weights present, "//&
               "meaning that there was a problem computing the weights (even if cell is empty, "//&
               "there should be one artificial weight corresponding to an artificial point)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       if(size(pidx)/=size(pw)) then
          write(errstr,*) "for inversion grid cell ",j,", there are ",size(pidx)," wavefield point indices, but ",&
               size(pw)," weights (should be equal!)"
          call add(errmsg,2,errstr,myname)
          return
       end if
       ! two write statements to allow for appropriate allocation in readIntegrationWeights
       write(lu) j,size(pidx)
       write(lu) pidx,pw
    enddo
    write(lu) this%unit_factor,this%ntot_wp
    write(lu) size(this%empty_cell)
    write(lu) this%empty_cell
    write(lu) size(this%err_level_weights)
    write(lu) this%err_level_weights
    write(lu) size(this%type_weights)
    write(lu) this%type_weights
    close(lu)
  end subroutine writeIntegrationWeights
!------------------------------------------------------------------------------------
!> \brief Read integration weights from file
!
  subroutine readIntegrationWeights(this,lu,filename,errmsg)
    type (integration_weights) :: this
    integer :: lu
    character (len=*) :: filename
    type (error_message) :: errmsg
    integer :: ierr,j,j_read,size_wp_idx,size_pidx,size_array
    integer, dimension(:), pointer :: pidx
    real, dimension(:), pointer :: pw
    real :: unit_factor
    integer :: ntot_wp
    character (len=22) :: myname = 'readIntegrationWeights'
    character (len=400) :: errstr
!
    nullify(pidx,pw)
    call addTrace(errmsg,myname)
!
    open(lu,file = trim(filename),form = 'unformatted',action='read',status='old',iostat = ierr)
    if (ierr > 0) then
       call add(errmsg,2,trim(filename)//' could not be opened',myname)
       return
    endif
!
    read(lu) size_wp_idx
    if (size_wp_idx <= 0) then
       write(errstr,*) 'number of inversion grid cells = ',size_wp_idx,' is invalid'
       call add(errmsg,2,errstr,myname)
       return
    endif
    !print *,'in readIntegrationWeights: size_wp_idx=',size_wp_idx
    allocate(this%wp_idx(size_wp_idx),this%weight(size_wp_idx))
!
    do j = 1,size_wp_idx
       read(lu) j_read,size_pidx
       if (j_read/=j) then
          write(errstr,*) "current inversion grid cell index read from file = ",j_read,", but expected ",j
          call add(errmsg,2,errstr,myname)
          return
       endif
       if (size_pidx <= 0) then
          write(errstr,*) "size of wp_idx vector of inversion grid cell ",j," is invalid: ",size_pidx
          call add(errmsg,2,errstr,myname)
          return
       endif
       allocate(pidx(size_pidx),pw(size_pidx))
       read(lu) pidx,pw
       call associateIntegerVectorPointer(this%wp_idx(j),pidx)
       call associateRealVectorPointer(this%weight(j),pw)
    enddo
    read(lu) unit_factor,ntot_wp
    if(unit_factor <= 0) then
       write(errstr,*) "unit factor of integration weights read from file = ",unit_factor,&
            "; must be > 0 , hence integration weights file seems to be inconsistent"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(ntot_wp <= 0) then
       write(errstr,*) "number of wavefield points read from file = ",ntot_wp,&
            "; is expected to be > 0 , hence integration weights file seems to be inconsistent"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    this%unit_factor = unit_factor
    this%ntot_wp = ntot_wp
    ! this%empty_cell
    read(lu) size_array
    if (size_array /= size_wp_idx) then
       write(errstr,*) 'size of empty_cells array ( = ',size_array,')does not equal number of inversion grid cells ',&
            '( = size wp_idx = ',size_wp_idx,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    allocate(this%empty_cell(size_array))
    read(lu) this%empty_cell
    ! this%err_level_weights
    read(lu) size_array
    if (size_array /= size_wp_idx) then
       write(errstr,*) 'size of err_level_weights array ( = ',size_array,')does not equal number of inversion grid cells ',&
            '( = size wp_idx = ',size_wp_idx,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    allocate(this%err_level_weights(size_array))
    read(lu) this%err_level_weights
    ! this%empty_type_weights
    read(lu) size_array
    if (size_array /= size_wp_idx) then
       write(errstr,*) 'size of type_weights array ( = ',size_array,')does not equal number of inversion grid cells ',&
            '( = size wp_idx = ',size_wp_idx,')'
       call add(errmsg,2,trim(errstr),myname)
       return
    endif
    allocate(this%type_weights(size_array))
    read(lu) this%type_weights
    close(lu)
  end subroutine readIntegrationWeights
!
end module integrationWeights
