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

!##################################################
! THIS IS A HARD-CODED EXECUTABLE, REQUIRING THIS
! SOURCE FILE TO BE ADAPTED (and re-compiled) 
! BEFORE USING IT!
!##################################################

program addSpikeCheckerToKim
  use kernelInvertedModel
  use errorMessage
  use inversionGrid

  implicit none

  character (len=400) :: filename_kim,filebase_kim_checker,invgrid_parfile
  type (kernel_inverted_model) :: kim,kim_checker_abs,kim_checker_rel
  real, dimension(:), pointer :: model_values_abs,model_values_rel

  type (inversion_grid) :: invgrid

  type (error_message) :: errmsg

  integer :: ilat,ncell_lat,ilon,ncell_lon,idepth,ncell_depth,ncell_total,icell
  integer :: nchecker_lat,ichecker_lat,nchecker_lon,ichecker_lon,nchecker_depth,ichecker_depth
  integer :: j,sum_mod_ichecker
  integer, dimension(:,:), allocatable :: ibchecker_lat,ibchecker_lon,ibchecker_depth
  real :: percent_anomaly,sign_checker_anomaly
  logical :: first_anomaly_is_positive

  ! HARDCODED QUICK-AND-DIRTY IMPLEMENTATION FOR schunkInversionGrid-type INVERSION GRIDS FOR WHICH 
  ! kim-MODELS ARE STORED STRAIGHTLY (all cell-index arrays sorted from 1 to ncell_invgrid)
  ! SHOULD ALSO WORK FOR chunkInversionGird-type INVERSION GRIDS OF 1 CHUNK WITHOUT CELL REFINEMENT

  ! AT THE MOMENT, THE CHECKERBOARD MODEL IS ONLY CREATED FOR vs AND THE GIVEN .kim FILE WHICH SHOULD
  ! BE SUPERIMPOSED BY A CHECKERBOARD IS ASSUMED TO HAVE PARAMETRIZATION isoVelocity1000 OR isoVelocitySI

  nullify(model_values_abs,model_values_rel)



!#######################
! BEGIN MANUAL INTITATION OF VARIABLES
!#######################

  first_anomaly_is_positive = .false.

  invgrid_parfile = '../schunk_invgrid_parfile'
  !invgrid_parfile = '../chunks_invgrid_parfile'

  filename_kim = 'krm_on_invgrid.kim'
  filebase_kim_checker = 'krm_checker_'


  ! manual specification of regular distribution of inversion grid cells (must coincide with actual inversion grid)
  ncell_lat = 51
  ncell_lon = 66
  ncell_depth = 20
  ncell_total = ncell_lat*ncell_lon*ncell_depth

  percent_anomaly = 5.0


  ! LATERAL RESOLUTION

  ! numbers of checkers in lat and lon direction
  nchecker_lat = 4
  nchecker_lon = 5
  allocate(ibchecker_lat(nchecker_lat,2),ibchecker_lon(nchecker_lon,2))

  ! for each checker in lat-direction, the minimum- and maximum indices of the inversion grid cells IN LAT DIRECTION (so, not global invgrid cell index, but only index running from 1 to ncell_lat)
  ibchecker_lat(1,1) = 2 ! minimum index
  ibchecker_lat(1,2) = 8 ! maximum index
  ibchecker_lat(2,1) = 16
  ibchecker_lat(2,2) = 22
  ibchecker_lat(3,1) = 30
  ibchecker_lat(3,2) = 36
  ibchecker_lat(4,1) = 44
  ibchecker_lat(4,2) = 50

  ibchecker_lon(1,1) = 2 ! minimum index
  ibchecker_lon(1,2) = 8 ! maximum index
  ibchecker_lon(2,1) = 16
  ibchecker_lon(2,2) = 22
  ibchecker_lon(3,1) = 30
  ibchecker_lon(3,2) = 36
  ibchecker_lon(4,1) = 44
  ibchecker_lon(4,2) = 50
  ibchecker_lon(5,1) = 58
  ibchecker_lon(5,2) = 64


  ! DEPTH RESOLUTION

  nchecker_depth = 5
  allocate(ibchecker_depth(nchecker_depth,2))

  ibchecker_depth(1,1) = 3 ! minimum index
  ibchecker_depth(1,2) = 4 ! maximum index
  ibchecker_depth(2,1) = 7
  ibchecker_depth(2,2) = 8
  ibchecker_depth(3,1) = 13
  ibchecker_depth(3,2) = 14
  ibchecker_depth(4,1) = 17
  ibchecker_depth(4,2) = 17
  ibchecker_depth(5,1) = 19
  ibchecker_depth(5,2) = 19


!#######################
! END OF MANUAL INTITATION OF VARIABLES
!#######################



  call new(errmsg,'addSpikeCheckerToKim')
  call createInversionGrid(invgrid,'schunkInversionGrid',invgrid_parfile,'',11,errmsg,recreate=.false.)
  !call createInversionGrid(invgrid,'chunksInversionGrid',invgrid_parfile,'../',11,errmsg,recreate=.false.)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call new(errmsg,'addSpikeCheckerToKim')
  call readFileKernelInvertedModel(kim,filename_kim,11,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call copyKernelInvertedModel(kim_checker_abs,kim)

  ! initiate relative checker update by setting all values to 0.0 (later on, only add percentages within the checkers)
  call copyKernelInvertedModel(kim_checker_rel,kim)
  model_values_rel => getVal(kim_checker_rel,'rho')
  if(.not.associated(model_values_rel)) then
     write(*,*) "ERROR: model_values_rel rho not associated"
     goto 1
  end if
  model_values_rel = 0.
  model_values_rel => getVal(kim_checker_rel,'vp')
  if(.not.associated(model_values_rel)) then
     write(*,*) "ERROR: model_values_rel vp not associated"
     goto 1
  end if
  model_values_rel = 0.
  model_values_rel => getVal(kim_checker_rel,'vs')
  if(.not.associated(model_values_rel)) then
     write(*,*) "ERROR: model_values_rel vs not associated"
     goto 1
  end if
  model_values_rel = 0.


  ! NOW LOOP ON ALL INVERSION GRID CELLS, ASSUMING THE MODEL VALUES ARE SORTED FOR INDICES 1,..,ncell
  ! FIND OUT IF THE CURRENT CELL IS WITHIN A CHECKER, IF YES THEN MANIPULATE THE MODE VALUES
  model_values_abs => getVal(kim_checker_abs,'vs')
  model_values_rel => getVal(kim_checker_rel,'vs')
  if(.not.associated(model_values_abs)) then
     write(*,*) "ERROR: model_values_abs vs not associated"
     goto 1
  end if
  if(size(model_values_abs) /= ncell_total) then
     write(*,*) "ERROR: model_values_abs vs have size ",size(model_values_abs),", but expected ",ncell_total
     goto 1
  end if
  icell = 0
  do idepth = 1,ncell_depth
     do ilon = 1,ncell_lon
        do ilat = 1,ncell_lat
           icell = icell+1
           if(cellInChecker()) then
              model_values_abs(icell) = model_values_abs(icell) * (1.0 + sign_checker_anomaly*percent_anomaly*0.01)
              model_values_rel(icell) = sign_checker_anomaly*percent_anomaly
           end if
        end do ! ilat
     end do ! ilon
  end do ! idepth


  ! WRITE ABSOLUTE AND RELATIVE CHECKER MODEL TO .kim FILES
  call new(errmsg,'addSpikeCheckerToKim')
  call writeFileKernelInvertedModel(kim_checker_abs,trim(filebase_kim_checker)//'abs.kim',11,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call new(errmsg,'addSpikeCheckerToKim')
  call writeFileKernelInvertedModel(kim_checker_rel,trim(filebase_kim_checker)//'rel.kim',11,errmsg)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  ! WRITE ABSOLUTE AND RELATIVE CHECKER MODEL TO Vtk FILES
  call new(errmsg,'addSpikeCheckerToKim')
  call writeVtkKernelInvertedModel(kim_checker_abs,invgrid,'BINARY',trim(filebase_kim_checker)//'abs',11,&
       errmsg,overwrite=.true.)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)

  call new(errmsg,'addSpikeCheckerToKim')
  call writeVtkKernelInvertedModel(kim_checker_rel,invgrid,'BINARY',trim(filebase_kim_checker)//'rel',11,&
       errmsg,overwrite=.true.)
  if(.level.errmsg /= 0) call print(errmsg)
  if(.level.errmsg == 2) goto 1
  call dealloc(errmsg)


  ! CLEAN UP

1 if(allocated(ibchecker_lat)) deallocate(ibchecker_lat)
  if(allocated(ibchecker_lon)) deallocate(ibchecker_lon)
  if(allocated(ibchecker_depth)) deallocate(ibchecker_depth)
  call dealloc(kim)
  call dealloc(errmsg)

  stop

contains

  function cellInChecker() result(l)
    logical :: l
    l = .false.

    ichecker_depth = 0
    do j = 1,nchecker_depth
       if(idepth < ibchecker_depth(j,1)) return
       if(idepth <= ibchecker_depth(j,2)) then
          ichecker_depth = j
          exit
       end if
    end do
    if(ichecker_depth == 0) return

    ichecker_lon = 0
    do j = 1,nchecker_lon
       if(ilon < ibchecker_lon(j,1)) return
       if(ilon <= ibchecker_lon(j,2)) then
          ichecker_lon = j
          exit
       end if
    end do
    if(ichecker_lon == 0) return

    ichecker_lat = 0
    do j = 1,nchecker_lat
       if(ilat < ibchecker_lat(j,1)) return
       if(ilat <= ibchecker_lat(j,2)) then
          ichecker_lat = j
          exit
       end if
    end do
    if(ichecker_lat == 0) return

    ! first find out if this checker cell has positive or negative model anomaly
    ! the following sum is an indicator of the sign of anomaly:
    sum_mod_ichecker = mod(ichecker_depth,2) + mod(ichecker_lon,2) + mod(ichecker_lat,2)
    ! this value being odd or even can decide about the sign. this way, you get a 3D checker pattern
    if( mod(sum_mod_ichecker,2) == 1 .eqv. first_anomaly_is_positive) then
       sign_checker_anomaly = 1.0
    else
       sign_checker_anomaly = -1.0
    end if

    l = .true.
  end function cellInChecker

end program addSpikeCheckerToKim
