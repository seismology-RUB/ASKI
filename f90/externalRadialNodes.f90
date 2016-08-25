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
!-----------------------------------------------------------------
!> \brief Module dealing with external radial nodes for Gemini II 
!
!  Definition of external radial nodes (in terms of depth):
!  For an arbitrary number of EXTERNAL_NODES_NBLOCKS blocks of nodes,
!  the vectors EXTERNAL_NODES_NNOD and EXTERNAL_NODES_DR (both of length EXTERNAL_NODES_NBLOCKS)
!  define the subdivision of each block by its number of nodes EXTERNAL_NODES_NNOD(i) contained 
!  in such a block, having constant spacing EXTERNAL_NODES_DR(i). The first node is always
!  at zero depth (surface) and does not count in EXTERNAL_NODES_NNOD. It can be shifted to greater depths.
!  Units: Length unit is always km.
!
!  Example: 
!      EXTERNAL_NODES_NBLOCKS =  3
!      EXTERNAL_NODES_NNOD =  8   5   2
!      EXTERNAL_NODES_DR =  5. 10. 20.
!  this means:
!  3 blocks of NODES (one below the other starting with the first BELOW the surface) with different spacing:
!     8 nodes with spacing of 5 km between surface and 40 km depth (5,10,15,20,25,30,35,40) plus one at surface
!     5 nodes with spacing of of 10 km between 40 km and 90 km depth (50,60,70,80,90)
!     2 layers with spacing of 20 km between 90 km and 130 km depth (110,130)
!----------------------------------------------------------------
module externalRadialNodes
    use errorMessage
    use inputParameter
    implicit none
    interface dealloc; module procedure deallocExternalRadialNodes; end interface
    interface operator (.nnod.); module procedure getNnodExternalRadialNodes; end interface
    interface getRadiiExternalRadialNodes
       module procedure getDoubleRadiiExternalRadialNodes
       module procedure getRealRadiiExternalRadialNodes
    end interface getRadiiExternalRadialNodes
    type external_radial_nodes
       private
       integer :: nblocks                                              ! total number of nodes
       integer, dimension(:), pointer :: nnodblocks => null()          ! array with nodes per block
       double precision, dimension(:), pointer :: drblocks => null()   ! array with node distance per block
       double precision :: shift                                       ! shift of "surface radius"
       double precision :: rearth                                      ! earth radius taken as surface
       double precision, dimension(:), pointer :: rnod => null()       ! nodes radii
       integer :: nnod                                                 ! number of radii
     end type external_radial_nodes
!
contains
!-------------------------------------------------------------------------------
!  \brief Read out parameters from file, nodes are created on the fly upon request
!
    subroutine createExternalRadialNodes(this,lu,parfile,errmsg)
    type (external_radial_nodes) :: this
    integer :: lu
    character (len=*) :: parfile
    type (error_message) :: errmsg
    ! local
    type (input_parameter) :: inpar
    integer :: ios,n,ib,i
    character (len=25) :: myname = 'createExternalRadialNodes'
    character (len=80), dimension(5) :: par_keys
    data par_keys/'EXTERNAL_NODES_NBLOCKS','EXTERNAL_NODES_NNOD',&
         'EXTERNAL_NODES_DR','EXTERNAL_NODES_SHIFT','EXTERNAL_NODES_REARTH'/
!
    call addTrace(errmsg,myname)
    call createKeywordsInputParameter(inpar,par_keys)
    call readSubroutineInputParameter(inpar,lu,parfile,errmsg)
    if (.level.errmsg == 2) then; call dealloc(inpar); return
    endif
    this%nblocks = ival(inpar,'EXTERNAL_NODES_NBLOCKS')
    this%nnodblocks => ivecp(inpar,'EXTERNAL_NODES_NNOD',this%nblocks,ios)
    if (ios /= 0) then
       call add(errmsg,2,'problems reading out EXTERNAL_NODES_NNOD',myname)
       goto 1
    endif
    this%drblocks => dvecp(inpar,'EXTERNAL_NODES_DR',this%nblocks,ios)
    if (ios /= 0) then
       call add(errmsg,2,'problems reading out EXTERNAL_NODES_DR',myname)
       goto 1
    endif
    this%shift = dval(inpar,'EXTERNAL_NODES_SHIFT')
    this%rearth = dval(inpar,'EXTERNAL_NODES_REARTH')
    call dealloc(inpar)
!
!  calculate node radii
!
    this%nnod = sum(this%nnodblocks)+1                ! +1 is surface node
    allocate(this%rnod(this%nnod))
    this%rnod(this%nnod) = this%rearth-this%shift
    n = this%nnod
    do ib = 1,this%nblocks
       do i = 1,this%nnodblocks(ib)
          this%rnod(n-i) = this%rnod(n)-i*this%drblocks(ib)
       enddo
       n = n-this%nnodblocks(ib)
    enddo
!
    if (any(this%rnod < 0.d0)) then
       call add(errmsg,2,'Negative radii, check your input',myname)
       goto 1
    endif
    return
!
1   if (associated(this%drblocks)) deallocate(this%drblocks)
    if (associated(this%nnodblocks)) deallocate(this%nnodblocks)
    if (associated(this%rnod)) deallocate(this%rnod)
    call dealloc(inpar)
    end subroutine createExternalRadialNodes
!------------------------------------------------------------------------
!  Deallocate object
!
    subroutine deallocExternalRadialNodes(this)
    type (external_radial_nodes) :: this
    if (associated(this%nnodblocks)) deallocate(this%nnodblocks)
    if (associated(this%drblocks)) deallocate(this%drblocks)
    end subroutine deallocExternalRadialNodes
!------------------------------------------------------------------------
!  Get total number of nodes
!
    integer function getNnodExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    res = this%nnod
    end function getNnodExternalRadialNodes
!------------------------------------------------------------------------
!  Get rearth
!
    double precision function getRearthExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    res = this%rearth
    end function getRearthExternalRadialNodes
!--------------------------------------------------------------------------
!  Return allocated double precision pointer to radii
!
    subroutine getDoubleRadiiExternalRadialNodes(this,rnod,errmsg)
    type (external_radial_nodes), intent(in) :: this
    double precision, dimension(:), pointer :: rnod
    type (error_message) :: errmsg
    character (len=33) :: myname = 'getDoubleRadiiExternalRadialNodes'
!
    call addTrace(errmsg,myname)
    allocate(rnod(this%nnod))
    rnod = this%rnod
    end subroutine getDoubleRadiiExternalRadialNodes
!--------------------------------------------------------------------------
!  Return double precision pointer to radii
!
    function getPointerDoubleRadiiExternalRadialNodes(this) result(res)
    type (external_radial_nodes), intent(in) :: this
    double precision, dimension(:), pointer :: res
    res => this%rnod
    end function getPointerDoubleRadiiExternalRadialNodes
!--------------------------------------------------------------------------
!  Return allocated real pointer to radii
!
    subroutine getRealRadiiExternalRadialNodes(this,rnod,errmsg)
    type (external_radial_nodes), intent(in) :: this
    real, dimension(:), pointer :: rnod
    type (error_message) :: errmsg
    character (len=31) :: myname = 'getRealRadiiExternalRadialNodes'
!
    call addTrace(errmsg,myname)
    allocate(rnod(this%nnod))
    rnod = sngl(this%rnod)
    end subroutine getRealRadiiExternalRadialNodes
!-------------------------------------------------------------------
!> \brief Get radius of uppermost node
!
    subroutine getDoubleRadiusTopNodeExternalRadialNodes(this,rtop)
    type (external_radial_nodes), intent(in) :: this
    double precision :: rtop
    rtop = this%rearth-this%shift
    end subroutine getDoubleRadiusTopNodeExternalRadialNodes
!
end module externalRadialNodes
