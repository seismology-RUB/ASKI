!----------------------------------------------------------------------------
!   Copyright 2013 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   and Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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
!> \brief Module to represent point grid where wavefield is sampled
!!
!! \details Generic module which forks to forward-method specific subroutines.
!!  The interfaces to method GEMINI still need to be adapted to the new ASKI version,
!!  so far only method SPECFEM3D is supported. 
!
module wavefieldPoints
  use inversionGrid
!!$  use geminiWavefieldPoints
  use specfem3dWavefieldPoints
  use realloc
  use errorMessage
  use fileUnitHandler
  implicit none
  interface dealloc; module procedure deallocateWavefieldPoints; end interface
  interface operator (.ntot.); module procedure getNtotWavefieldPoints; end interface
  type wavefield_points
     private
!!$     type (gemini_wavefield_points), pointer :: gemini_wp => null()
     type (specfem3d_wavefield_points), pointer :: specfem_wp => null()
  end type wavefield_points
!
contains
!---------------------------------------------------------------------------------------
!> \brief Create a wavefield_point object from gemini_wavefield_points
!
  subroutine createWavefieldPoints(this,method,fuh,filename,errmsg)
    type (wavefield_points) :: this
    character(len=*) :: method
    type (file_unit_handler) :: fuh
    character (len=*) :: filename
    type (error_message) :: errmsg
    character (len=21) :: myname = 'createWavefieldPoints'
!
    call addTrace(errmsg,myname)
    select case (method)
!!$    case('GEMINI')
!!$       allocate(this%gemini_wp)
!!$       errmsg = readGeminiWavefieldPoints(this%gemini_wp,fuh,filename)
!!$       if (.level.errmsg == 2) then
!!$          call addTrace(errmsg,myname); return
!!$       endif
    case('SPECFEM3D')
       allocate(this%specfem_wp)
       call readSpecfem3dWavefieldPoints(this%specfem_wp,fuh,filename,errmsg)
       if (.level.errmsg == 2) return
    case default
       call add(errmsg,2,"Invalid forward computation method '"//trim(method)//"'",myname)
       return
    end select
  end subroutine createWavefieldPoints
!--------------------------------------------------------------------------------
!> \brief Deallocate wavefield points
!
  subroutine deallocateWavefieldPoints(this)
    type (wavefield_points), intent(in) :: this
!!$    if (associated(this%gemini_wp)) then
!!$       call deallocGeminiWavefieldPoints(this%gemini_wp)
!!$    else if(associated(this%specfem_wp)) then
    if(associated(this%specfem_wp)) & 
         call dealloc(this%specfem_wp)
!!$    endif
  end subroutine deallocateWavefieldPoints
!-------------------------------------------------------------------
!> \brief Get total number of wavefield points
!
  function getNtotWavefieldPoints(this) result(res)
    type (wavefield_points), intent(in) :: this
    integer :: res
!!$    if (associated(this%gemini_wp)) then
!!$       res = getNtotGeminiWavefieldPoints(this%gemini_wp)
!!$    else if(associated(this%specfem_wp)) then
       res = getNtotSpecfem3dWavefieldPoints(this%specfem_wp)
!!$    endif
  end function getNtotWavefieldPoints
!---------------------------------------------------------------------
!> \brief Get all wavefield points
!
  subroutine getWavefieldPoints(this,c1,c2,c3)
    type (wavefield_points), intent(in) :: this
    real, dimension(:), pointer :: c1,c2,c3
    !
!!$    if (associated(this%gemini_wp)) then
!!$       call getGeminiWavefieldPoints(this%gemini_wp,c1,c2,c3)
!!$    else if(associated(this%specfem_wp)) then
       call getSpecfem3dWavefieldPoints(this%specfem_wp,c1,c2,c3)
!!$    endif
  end subroutine getWavefieldPoints
!------------------------------------------------------------------------
!> \brief Get point-geometry information of wavefield points for vtk file output
!
  subroutine getVtkWavefieldPoints(this,invgrid,points,wp_indx_out,errmsg,wp_indx_req,indx_map_out)
    ! incoming
    type (wavefield_points) :: this
    type (inversion_grid) :: invgrid
    integer, dimension(:), optional :: wp_indx_req
    ! outgoing
    real, dimension(:,:), pointer :: points
    integer, dimension(:), pointer :: wp_indx_out
    type (error_message) :: errmsg
    integer, dimension(:), pointer, optional :: indx_map_out
    ! local
    character(len=21) :: myname = 'getVtkWavefieldPoints'
    real, dimension(:), pointer :: c1,c2,c3
    integer :: nwp,npoints,i
    logical, dimension(:), allocatable :: wp_indx_already_found
!
    call addTrace(errmsg,myname)
    nullify(points,wp_indx_out)
    if(present(indx_map_out)) nullify(indx_map_out)
!
    ! get all wavefield points
    call getWavefieldPoints(this,c1,c2,c3)
    ! transform wavefield points into inversion grid vtk coordinates for plotting in the same way, as the inversion grid is plotted
    call transformToVtkInversionGrid(invgrid,c1,c2,c3,'wp',errmsg)
    if(.level.errmsg==2) goto 1
    nwp = size(c1)
!
    ! now check requested indices wp_indx_req and create the respective index array wp_indx_out
    if(present(wp_indx_req)) then
       ! select valid indices for which points will be returned. remove duplicate indices,
       ! where only the first occurring wavefield point index is returned, the later occurring ones others are lost
       !
       ! allocate for maximally possible number of points
       allocate(wp_indx_out(size(wp_indx_req)))
       if(present(indx_map_out)) allocate(indx_map_out(size(wp_indx_req)))
       allocate(wp_indx_already_found(nwp)); wp_indx_already_found(:) = .false.
       npoints = 0
       ! loop over icoming requested wavefield point indices and select only valid ones and remove duplicate ones
       do i = 1,size(wp_indx_req)
          if(wp_indx_req(i) .ge. 1 .and. wp_indx_req(i) .le. nwp) then
             if(.not.wp_indx_already_found(wp_indx_req(i))) then
                wp_indx_already_found(wp_indx_req(i)) = .true.
                npoints = npoints + 1
                wp_indx_out(npoints) = wp_indx_req(i)
                if(present(indx_map_out)) indx_map_out(npoints) = i
             endif
          endif
       enddo
       deallocate(wp_indx_already_found)
       ! if there were no valid indices in wp_indx_req, return nullified pointers
       if(npoints == 0) then
          ! return no points (nullified pointer)
          deallocate(wp_indx_out); nullify(wp_indx_out)
          return
       endif
       ! reallocate if necessary
       if(npoints<size(wp_indx_req)) then
          wp_indx_out => reallocate(wp_indx_out,npoints)
          if(present(indx_map_out)) indx_map_out => reallocate(indx_map_out,npoints)
       end if
    else ! present(wp_indx_req)
       ! return all points
       allocate(wp_indx_out(nwp))
       wp_indx_out = (/ (i,i=1,nwp) /)
       if(present(indx_map_out)) then
          allocate(indx_map_out(nwp))
          indx_map_out = (/ (i,i=1,nwp) /)
       end if
    endif ! present(wp_indx_req)
    !
    ! return the requested points
    allocate(points(3,size(wp_indx_out)))
    points(1,:) = c1(wp_indx_out)
    points(2,:) = c2(wp_indx_out)
    points(3,:) = c3(wp_indx_out)
!
1   if(associated(c1)) deallocate(c1)
    if(associated(c2)) deallocate(c2)
    if(associated(c3)) deallocate(c3)
    !
  end subroutine getVtkWavefieldPoints
!
 end module wavefieldPoints
