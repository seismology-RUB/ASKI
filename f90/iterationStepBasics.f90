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
!> \brief hold basic requirements for a specific iteration step
!!
!! \details Basic requirements for a specific iteration step, such as the iteration 
!!  specific parfile, wavefield points, inversion grid, kernel reference model, 
!!  integration weights, etc. are supervised by this module
!!
!! \author Florian Schumacher
!! \date March 2013
!
module iterationStepBasics
!
  use inversionBasics
  use inputParameter
  use wavefieldPoints
  use wpVtkFile
  use inversionGrid
  use invgridVtkFile
  use eventStationVtkFile
  use integrationWeights
  use kernelReferenceModel
  use kernelInvertedModel
  use fileUnitHandler
  use errorMessage
!
  implicit none
!
  interface dealloc; module procedure deallocateIterationStepBasics; end interface
  interface init; module procedure initiateIterationStepBasics; end interface
  interface mapIfreq2ArrayIndex
     module procedure mapOneIfreq2ArrayIndexIterationStepBasics
     module procedure mapManyIfreq2ArrayIndexIterationStepBasics
  end interface mapIfreq2ArrayIndex
  interface operator (.iterpath.); module procedure getIterationStepPathIterationStepBasics; end interface
  interface operator (.inpar.); module procedure getInputParameterIterationStepBasics; end interface
  interface operator (.wp.); module procedure getWavefieldPointsIterationStepBasics; end interface
  interface operator (.invgrid.); module procedure getInversionGridIterationStepBasics; end interface
  interface operator (.intw.); module procedure getIntegrationWeightsIterationStepBasics; end interface
  interface operator (.krm.); module procedure getKernelReferenceModelIterationStepBasics; end interface
  interface operator (.nf.); module procedure getNumberOfFrequenciesIterationStepBasics; end interface
  interface operator (.ifreq.); module procedure getFrequencyIndicesIterationStepBasics; end interface
!
!> \brief type contains all basic requirements for a specific iteration step
  type iteration_step_basics
     private
     character(len=350) :: iter_path !< absolute iteration step path (starting with inverison main path)
     type (input_parameter) :: inpar !< parfile content
     type (wavefield_points) :: wp !< wavefield points used for inversion
     type (inversion_grid) :: invgrid !< inversion grid used for inversion
     type (integration_weights) :: intw !< integration weights used for inversion
     type (kernel_reference_model) :: krm !< reference model used in this iteration step
     integer, dimension(:), pointer :: ifreq => null() !< iteration step frequency indices
  end type iteration_step_basics
!
contains
!
!------------------------------------------------------------------------
!> \brief initiate basic iteration step requirements
!! \details Pathfile and parfile of this iteration step are read. Every
!!  basic information needed for an iteration step is read or created.
!! \param this interation step basics
!! \param invbasics inversion basics
!! \return errmsg error message
!
  subroutine initiateIterationStepBasics(this,invbasics,fuh,errmsg,recreate_existing_files)
    ! incoming
    type (iteration_step_basics) :: this
    type (inversion_basics) :: invbasics
    type (file_unit_handler) :: fuh
    logical, optional :: recreate_existing_files
    ! returning
    type (error_message) :: errmsg
    ! local
       ! integration weights
    logical :: intw_existed
       ! vtk files
    integer, dimension(:), pointer :: idata,cells_empty,cells_filled,wp_outside,wp_inside
    real, dimension(:), pointer :: rdata
    type (invgrid_vtk_file) :: invgrid_vtk
    type (wp_vtk_file) :: wp_vtk
    type (event_station_vtk_file) :: evstat_vtk
       ! reference model
    character(len=character_length_pmtrz) :: parametrization
    logical :: kim_krm_existed
    type (kernel_inverted_model) :: kim_krm
    character(len=character_length_param) :: param_name
       ! other
    logical :: recreate_files
    integer :: i_partest
    integer :: j,n,ios,nf
    real :: volume
    character (len=80), dimension(14) :: iter_inpar_keys
    type (error_message) :: errmsg2
    character(len=400) :: errstr
    character(len=27) :: myname = 'initiateIterationStepBasics'
    !
    ! keywords for input parameter
    data iter_inpar_keys/'FILE_WAVEFIELD_POINTS', 'FILE_INTEGRATION_WEIGHTS', 'TYPE_INVERSION_GRID', &
         'PARFILE_INVERSION_GRID', 'FILE_KERNEL_REFERENCE_MODEL', 'ITERATION_STEP_NUMBER_OF_FREQ', &
         'ITERATION_STEP_INDEX_OF_FREQ', 'PATH_KERNEL_DISPLACEMENTS', 'PATH_KERNEL_GREEN_TENSORS', &
         'PATH_OUTPUT_FILES', 'PATH_SENSITIVITY_KERNELS', 'PATH_SYNTHETIC_DATA', 'TYPE_INTEGRATION_WEIGHTS', &
         'FILEBASE_BASIC_STATS'/
!
!  function starts here
!
    call addTrace(errmsg,myname)
    if(present(recreate_existing_files)) then
       recreate_files = recreate_existing_files
    else
       recreate_files = .false.
    end if
    this%iter_path = .iterpath.invbasics
!
!  read input parameters
!
    call createKeywordsInputParameter(this%inpar,iter_inpar_keys)
    call readSubroutineInputParameter(this%inpar,get(fuh),trim(this%iter_path)//&
         trim((.inpar.invbasics).sval.'PARFILE_ITERATION_STEP'),errmsg)
    call undo(fuh)
    if (.level.errmsg == 2) return
!
    ! check consistency of entries in iteration step specific parfile
    ! during creation of object this%inpar, there was already checked if the required keywords are present
    ! now check correct type of values here (except character types) and consistencies
!
    nf = ival(this%inpar,'ITERATION_STEP_NUMBER_OF_FREQ',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'ITERATION_STEP_NUMBER_OF_FREQ' = '"//trim(this%inpar.sval.&
            'ITERATION_STEP_NUMBER_OF_FREQ')//"' of iteration step parfile is not a valid integer value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(nf .le. 0) then
       write(errstr,*) "parameter 'ITERATION_STEP_NUMBER_OF_FREQ' = ",nf,&
            " of iteration step parfile is not valid; must be positive"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(nf > ((.inpar.invbasics).ival.'MEASURED_DATA_NUMBER_OF_FREQ')) then
       write(errstr,*) "parameter 'ITERATION_STEP_NUMBER_OF_FREQ' = ",nf,&
            " of iteration step parfile must not be larger than 'MEASURED_DATA_NUMBER_OF_FREQ' = ",&
            (.inpar.invbasics).ival.'MEASURED_DATA_NUMBER_OF_FREQ'," of main parfile"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
!
    this%ifreq => ivecp(this%inpar,'ITERATION_STEP_INDEX_OF_FREQ',nf,iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'ITERATION_STEP_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
            'ITERATION_STEP_INDEX_OF_FREQ')//"' is not a vector of ",nf," integers"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(any(this%ifreq < 0)) then
       write(errstr,*) "entries of vector 'ITERATION_STEP_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
            'ITERATION_STEP_INDEX_OF_FREQ')//"' must not be negative"
       call add(errmsg,2,trim(errstr),myname)
       deallocate(this%ifreq)
       return
    end if
    do j = 1,nf
       if(any(this%ifreq(j+1:nf) == this%ifreq(j))) then
          write(errstr,*) "vector 'ITERATION_STEP_INDEX_OF_FREQ' = '"//trim(this%inpar.sval.&
               'ITERATION_STEP_INDEX_OF_FREQ')//"' must not contain multiple entries"
          call add(errmsg,2,trim(errstr),myname)
          deallocate(this%ifreq)
          return
       end if
       if(.not.any((.ifreq.invbasics)==this%ifreq(j))) then
          write(errstr,*) "frequency index ",this%ifreq(j)," contained in vector 'ITERATION_STEP_INDEX_OF_FREQ' = '"&
               //trim(this%inpar.sval.'ITERATION_STEP_INDEX_OF_FREQ')//&
               "' is not contained in global frequency index vector 'MEASURED_DATA_INDEX_OF_FREQ' of main parfile"
          call add(errmsg,2,trim(errstr),myname)
          deallocate(this%ifreq)
          return
       end if
    end do ! j
!
    i_partest = ival(this%inpar,'TYPE_INTEGRATION_WEIGHTS',iostat=ios)
    if(ios/=0) then
       write(errstr,*) "parameter 'TYPE_INTEGRATION_WEIGHTS' = '"//trim(this%inpar.sval.'TYPE_INTEGRATION_WEIGHTS')//&
            "' of iteration step parfile is not a valid integer value"
       call add(errmsg,2,trim(errstr),myname)
       return
    end if
    if(.not.validTypeIntegrationWeights(i_partest,err=errmsg)) return
!
    ! check dependent validity of inversion grid, wavefield points (i.e. FORWARD_METHOD) and integration weights
    if(.not.validTypeInversionGrid(this%inpar.sval.'TYPE_INVERSION_GRID',&
         method=(.inpar.invbasics).sval.'FORWARD_METHOD',&
         intw_type=this%inpar.ival.'TYPE_INTEGRATION_WEIGHTS',err=errmsg)) return
!
! create wavefield points
!
    call createWavefieldPoints(this%wp,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
         trim(this%iter_path)//trim(this%inpar.sval.'FILE_WAVEFIELD_POINTS'),errmsg)
    if (.level.errmsg == 2) return
!
! create inversion grid
!
    call createInversionGrid(this%invgrid,this%inpar.sval.'TYPE_INVERSION_GRID',&
         trim(this%iter_path)//trim(this%inpar.sval.'PARFILE_INVERSION_GRID'),this%iter_path,&
         get(fuh),errmsg,recreate=recreate_files)
    call undo(fuh)
    if (.level.errmsg == 2) return
!
!  create integration weights (if not existent yet, or if recreate_files=True)
!
    inquire(file=trim(this%iter_path)//trim(this%inpar.sval.'FILE_INTEGRATION_WEIGHTS'), exist=intw_existed)
    if(intw_existed .and. .not.recreate_files) then
       call readIntegrationWeights(this%intw,get(fuh),trim(this%iter_path)//&
            trim(this%inpar.sval.'FILE_INTEGRATION_WEIGHTS'),errmsg)
       call undo(fuh)
       if (.level.errmsg == 2) return
    else ! intw_existed .and. .not.recreate_files
       ! create integration weights
       call createIntegrationWeights(this%intw,this%inpar.ival.'TYPE_INTEGRATION_WEIGHTS',this%wp,this%invgrid,errmsg)
       if (.level.errmsg == 2) return
       !
       ! write integration weights
       call writeIntegrationWeights(this%intw,get(fuh),trim(this%iter_path)//&
            trim(this%inpar.sval.'FILE_INTEGRATION_WEIGHTS'),errmsg)
       call undo(fuh)
       if (.level.errmsg == 2) return
    endif ! intw_existed .and. invgrid_existed
!
!  check consistency of integration weights and inversion grid, wavefield points
!
    if(.ncell.this%intw /= .ncell.this%invgrid) then
       write(errstr,*) "number of invgrid cells contained in integration weights = ",.ncell.this%intw,&
            " differs from actual number of cells of inversion grid  = ",.ncell.this%invgrid,&
            " -> files are inconsistent, so recreate files by running 'initBasics -recr ...'"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(.nwp.this%intw /= .ntot.this%wp) then
       write(errstr,*) "number of wavefield points contained in integration weights = ",.nwp.this%intw,&
            " differs from actual number wavefield points = ",.ntot.this%wp,&
            " -> files are inconsistent, so recreate files by running 'initBasics -recr ...'"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
! write vtk files, if not existendyou have new inversion grid or new integration weights or recreate_files=True
!
    if(.not.intw_existed .or. recreate_files) then
       !
       ! in general, write all vtk files only for invgrid cells containing wavefield points, and wavefield points inside the inversion grid
       cells_filled => getFilledCells(this%intw)
       if(.not.associated(cells_filled)) then
          call add(errmsg,2,"there are only empty inversion grid cells",myname)
          return
       end if
       wp_inside => getWpInside(this%intw)
       if(.not.associated(wp_inside)) then
          call add(errmsg,2,"there no wavefield points inside the inversion grid",myname)
          return
       end if
       !
       ! wavefield points
       call init(wp_vtk,this%wp,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS'),&
            trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,'Wavefield points',&
            wp_indx_req=wp_inside)
       if (.level.errmsg == 2) return
       ! write wavefield points tp vtk file
       call writeWp(wp_vtk,get(fuh),errmsg,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(wp_vtk)
       !
       ! if there are any wavefield points outside the inversion grid, write additionally all the wavefield points
       ! (including those outside the inversion grid) and (separately) only the points outside
       if(anyWpOutside(this%intw)) then
          ! all wavefield points
          call init(wp_vtk,this%wp,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_all',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
               vtk_title='Wavefield points including points outside invgrid')
          if (.level.errmsg == 2) return
          ! write wavefield points tp vtk file
          call writeWp(wp_vtk,get(fuh),errmsg,overwrite=.true.)
          call undo(fuh)
          if (.level.errmsg == 2) return
          call dealloc(wp_vtk)
          !
          ! wavefield points outside invgrid
          wp_outside => getWpOutside(this%intw)
          call init(wp_vtk,this%wp,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_outside',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
               vtk_title='Wavefield points including points outside invgrid',wp_indx_req=wp_outside)
          if (.level.errmsg == 2) return
          ! write wavefield points tp vtk file
          call writeWp(wp_vtk,get(fuh),errmsg,overwrite=.true.)
          call undo(fuh)
          if (.level.errmsg == 2) return
          call dealloc(wp_vtk)
          if(associated(wp_outside)) deallocate(wp_outside)
          !
       endif ! anyWpOutside(this%intw)
       !
       ! inversion grid
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS'),&
            trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title='inversion grid',&
            cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeInvgrid(invgrid_vtk,get(fuh),errmsg,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       !  if there are any empty inversion grid cells, write additionally the whole inversion grid (including those empty cells)
       !  and (seperately) the empty cells
       if(anyCellEmpty(this%intw)) then
          ! total inversion grid (including empty cells)
          call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_total',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),&
               errmsg,vtk_title='inversion grid including empty cells')
          if (.level.errmsg == 2) return
          call writeInvgrid(invgrid_vtk,get(fuh),errmsg,overwrite=.true.)
          call undo(fuh)
          if (.level.errmsg == 2) return
          call dealloc(invgrid_vtk)
          !
          ! empty cells
          cells_empty => getEmptyCells(this%intw)
          call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_empty_cells',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
               vtk_title='empty inversion grid cells',cell_indx_req=cells_empty)
          if (.level.errmsg == 2) return
          call writeInvgrid(invgrid_vtk,get(fuh),errmsg,overwrite=.true.)
          call undo(fuh)
          if (.level.errmsg == 2) return
          call dealloc(invgrid_vtk)
          if(associated(cells_empty)) deallocate(cells_empty)
          !
       endif ! anyCellEmpty(this%intw)
       !
       ! write inversion grid cell indices to invgridVtkFile
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_invgrid_cell_indices',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='cell indices of inversion grid',cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeData(invgrid_vtk,get(fuh),real(cells_filled),errmsg,data_name='cell_index',overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       ! write events.vtk (depdent on (.inpar.invbasics).lval.'USE_LOCAL_INVGRID_COORDS_FOR_VTK', 
       ! this is invgrid dependent, so use invgrid file name as base name)
       call init(evstat_vtk,.evlist.invbasics,this%invgrid,trim(this%iter_path)//&
            trim(this%inpar.sval.'FILEBASE_BASIC_STATS'),&
            trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title='events')
       if (.level.errmsg == 2) return
       call writeEvents(evstat_vtk,get(fuh),errmsg,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(evstat_vtk)
       !
       ! write stations.vtk (depdent on (.inpar.invbasics).lval.'USE_LOCAL_INVGRID_COORDS_FOR_VTK', 
       ! this is invgrid dependent, so use invgrid file name as base name)
       call init(evstat_vtk,.statlist.invbasics,this%invgrid,trim(this%iter_path)//&
            trim(this%inpar.sval.'FILEBASE_BASIC_STATS'),&
            trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title='stations')
       if (.level.errmsg == 2) return
       call writeStations(evstat_vtk,get(fuh),errmsg,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(evstat_vtk)
       !
       ! write number of wavefield points per invgrid box to invgridVtkFile
       idata => getNwpPerBoxIntegrationWeights(this%intw)
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_nwp_per_cell',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='Number of wavefield Points per inversion grid cell',cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeData(invgrid_vtk,get(fuh),real(idata),errmsg,data_name='number_of_wavefield_points',&
            data_indx_is_invgrid=.true.,overwrite=.true.)
       call undo(fuh)
       if(associated(idata)) deallocate(idata)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       ! write sum of integration weights per invgrid cell (should correspond to volume) to invgridVtkFile
       ! also compute the sum of weights for empty cells, as it is easiest to input all data (of size .ncell.invgrid)
       ! and control to leave out empty cells by setting cell_indx_req=cells_filled (and then data_indx_is_invgrid=.true.)
       allocate(rdata(.ncell.(this%invgrid)))
       do j = 1,.ncell.(this%invgrid)
          rdata(j) = sum((this%intw).weight.j)
       end do
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_sumw_per_cell',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='Sum of integration weights per (~ volume of) inversion grid cell',cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeData(invgrid_vtk,get(fuh),rdata,errmsg,data_name='sum_of_integration_weights',&
            data_indx_is_invgrid=.true.,overwrite=.true.)
       call undo(fuh)
       deallocate(rdata)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       ! write cell volume per invgrid cell to invgridVtkFile
       ! get the volume for all available cells, and set cell_indx_req respectively. Then set data_indx_is_invgrid=.false.
       allocate(rdata(.ncell.(this%invgrid)),idata(.ncell.(this%invgrid)))
       n = 0
       do j = 1,.ncell.(this%invgrid)
          call new(errmsg2,myname)
          call getVolumeCellInversionGrid(this%invgrid,j,volume,errmsg2)
          if(.level.errmsg2 == 0) then
             n = n+1
             idata(n) = j
             rdata(n) = volume
          end if
          call dealloc(errmsg2)
       end do ! j
       if(n > 0) then
          call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_cell_vol',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
               vtk_title='Volume of inversion grid cell',cell_indx_req=idata(1:n))
          if (.level.errmsg == 2) return
          call writeData(invgrid_vtk,get(fuh),rdata(1:n),errmsg,data_name='cell_volume',&
               data_indx_is_invgrid=.false.,overwrite=.true.)
          call undo(fuh)
          deallocate(rdata,idata)
          if (.level.errmsg == 2) return
          call dealloc(invgrid_vtk)
       else
          ! raise warning, no volumes can be computed
          call add(errmsg,1,"cell volume could not be computed for any inversion grid cell",myname)
          deallocate(rdata,idata)
       end if
       if(n < .ncell.(this%invgrid)) then
          ! raise warning, for some cells the volume could not be computed
          write(errstr,*) "for ",.ncell.(this%invgrid)-n," out of ",.ncell.(this%invgrid),&
               " inversion grid cells, their volume could not be computed"
          call add(errmsg,1,errstr,myname)
       end if
       !
       ! write integration weights on wavefield points
       rdata => getAtWavefieldPointsIntegrationWeights(this%intw)
       call init(wp_vtk,this%wp,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_intw_on_wp',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='Integration weights on wavefield points',wp_indx_req=wp_inside)
       if (.level.errmsg == 2) return
       call writeData(wp_vtk,get(fuh),rdata,errmsg,'integration_weights',data_indx_is_wp=.true.,overwrite=.true.)
       call undo(fuh)
       deallocate(rdata)
       if (.level.errmsg == 2) return
       call dealloc(wp_vtk)
       !
       ! write type of weights (per invgrid cell)
       idata => getTypeIntegrationWeights(this%intw)
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_intw_type',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='type of integration weights',cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeData(invgrid_vtk,get(fuh),real(idata),errmsg,data_name='type_of_weights',&
            data_indx_is_invgrid=.true.,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       ! write error level returned from computation of weights (per invgrid cell)
       idata => getErrLevelIntegrationWeights(this%intw)
       call init(invgrid_vtk,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_intw_err_level',trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,&
            vtk_title='Error level returned from computation of weights',cell_indx_req=cells_filled)
       if (.level.errmsg == 2) return
       call writeData(invgrid_vtk,get(fuh),real(idata),errmsg,data_name='weights_error_level',&
            data_indx_is_invgrid=.true.,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(invgrid_vtk)
       !
       if(associated(cells_filled)) deallocate(cells_filled)
       if(associated(wp_inside)) deallocate(wp_inside)
    endif ! .not.intw_existed .or. .not.invgrid_existed
    !
    !  create kernel reference model
    !
    call createKernelReferenceModel(this%krm,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
         trim(this%iter_path)//trim(this%inpar.sval.'FILE_KERNEL_REFERENCE_MODEL'),errmsg)
    if (.level.errmsg == 2) return
    !
    ! write kernel reference model interpolated on inversion grid as kim file and vtk files
    !(filename base is file kernel reference model + 'on invgrid' !)
    !
    inquire(file=trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')//'_krm_on_invgrid.kim',&
         exist=kim_krm_existed)
    if(.not.kim_krm_existed .or. recreate_files) then
       parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
       !
       ! kreate kernel inverted model object by interpolating kernel reference model to inversion grid
       call interpolateKernelReferenceToKernelInvertedModel(kim_krm,this%krm,&
            parametrization,this%invgrid,this%intw,errmsg)
       if (.level.errmsg == 2) return
       ! write kernel inverted model file
       call writeFileKernelInvertedModel(kim_krm,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
            //'_krm_on_invgrid.kim',get(fuh),errmsg)
       call undo(fuh)
       if (.level.errmsg == 2) return
       ! write kernel inverted model vtk files
       call writeVtkKernelInvertedModel(kim_krm,this%invgrid,(.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT',&
            trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')//'_krm_on_invgrid.kim',&
            get(fuh),errmsg,overwrite=.true.)
       call undo(fuh)
       if (.level.errmsg == 2) return
       call dealloc(kim_krm)
       !
       ! write kernel reference model on wavefield points as wpVtkFile (filename base is file kernel reference model + 'on wp' !)
       !
       do while(nextParamModelParametrization(parametrization,param_name))
          !
          ! get model values for this parameter
          rdata => getModelValuesKernelReferenceModel(this%krm,parametrization,param_name)
          if(.not.associated(rdata)) then
             call add(errmsg,1,"there are no kernel reference model values for '"//&
                  trim(parametrization)//"' parameter '"//trim(param_name)//"'",myname)
             cycle
          end if
          call init(wp_vtk,this%wp,this%invgrid,trim(this%iter_path)//trim(this%inpar.sval.'FILEBASE_BASIC_STATS')&
               //'_krm_on_wp_'//trim(parametrization)//"-"//trim(param_name),&
               trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),&
               errmsg,vtk_title=trim(parametrization)//'-'//trim(param_name)//' model values on wavefield points')
          if (.level.errmsg == 2) return
          call writeData(wp_vtk,get(fuh),rdata,errmsg,trim(parametrization)//'-'//trim(param_name),&
               data_indx_is_wp=.true.,overwrite=.true.)
          call undo(fuh)
          if (.level.errmsg == 2) return
          call dealloc(wp_vtk)
       end do ! while(nextParam)
    end if ! .not.kim_krm_existed .or. recreate_files
!
  end subroutine initiateIterationStepBasics
!------------------------------------------------------------------------
!> \brief deallocate iteration step basics object
!! \param this iteration step basics
!
  subroutine deallocateIterationStepBasics(this)
    type (iteration_step_basics) :: this
    call dealloc(this%inpar)
    call dealloc(this%wp)
    call dealloc(this%invgrid)
    call dealloc(this%intw)
    call dealloc(this%krm)
    if(associated(this%ifreq)) deallocate(this%ifreq)
  end subroutine deallocateIterationStepBasics
!------------------------------------------------------------------------
!> \brief for one given frequency index ifreq_in return its index position in array this%ifreq
!! \param this iteration step basics
!! \param ifreq_in integer frequency index for which the array position in this%ifreq should be returned
!! \param idx integer such that this%ifreq(idx)==ifreq_in, or -1 if no entry of this%ifreq has value ifreq_in
!! \return position index of ifreq_in in array this%ifreq, or -1 if no entry of this%ifreq has value ifreq_in
!
  function mapOneIfreq2ArrayIndexIterationStepBasics(this,ifreq_in) result(idx)
    type (iteration_step_basics), intent(in) :: this
    integer :: ifreq_in,idx
    ! local
    integer :: j
    idx = -1
    if(.not.associated(this%ifreq)) return
    do j=1,size(this%ifreq)
       if(this%ifreq(j) == ifreq_in) then
          idx = j
          exit
       end if
    end do ! j
  end function mapOneIfreq2ArrayIndexIterationStepBasics
!------------------------------------------------------------------------
!> \brief for a vector of given frequency indices ifreq_in return their index positions in array this%ifreq
!! \param this iteration step basics
!! \param ifreq_in vector of integer frequency indices for which the array positions in this%ifreq should be returned
!! \param idx of same length as ifreq_in and this%ifreq(idx(i))==ifreq_in(i), or not associated if any ifreq_in(i) is not in this%ifreq
!! \return position indices of values in ifreq_in in array this%ifreq, or not associated if any ifreq_in(i) is not in this%ifreq
!
  function mapManyIfreq2ArrayIndexIterationStepBasics(this,ifreq_in) result(idx)
    type (iteration_step_basics), intent(in) :: this
    integer, dimension(:) :: ifreq_in
    integer, dimension(:), pointer :: idx
    ! local
    integer :: size_ifreq_in,i,j
!
    idx => null()
    if(.not.associated(this%ifreq)) return
    size_ifreq_in = size(ifreq_in)
    if(size_ifreq_in <= 0) return
!
    allocate(idx(size_ifreq_in))
    idx = -1
    do i=1,size_ifreq_in
       do j=1,size(this%ifreq)
          if(this%ifreq(j) == ifreq_in(i)) then
             idx(i) = j
             exit
          end if
       end do ! j
       if(idx(i)==-1) then
          deallocate(idx); nullify(idx)
          return
       end if
    end do ! i
  end function mapManyIfreq2ArrayIndexIterationStepBasics
!------------------------------------------------------------------------
!> \brief get iteration step path
!! \param this iteration_step_basics object
!! \param iterpath iteration step path
!! \return iteration step path contained in this
!
  function getIterationStepPathIterationStepBasics(this) result(iterpath)
    type (iteration_step_basics), intent(in) :: this
    character(len=350) :: iterpath
    iterpath = this%iter_path
  end function getIterationStepPathIterationStepBasics
!------------------------------------------------------------------------
!> \brief get input_parameter object contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \param input_parameter object
!! \return input_parameter object contained in this
!
  function getInputParameterIterationStepBasics(this) result(inpar)
    type (iteration_step_basics), intent(in) :: this
    type (input_parameter) :: inpar
    inpar = this%inpar
  end function getInputParameterIterationStepBasics
!------------------------------------------------------------------------
!> \brief get wavefield points contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \return wp wavefield points contained in this
!
  function getWavefieldPointsIterationStepBasics(this) result(wp)
    type (iteration_step_basics), intent(in) :: this
    type (wavefield_points) :: wp
    wp = this%wp
  end function getWavefieldPointsIterationStepBasics
!------------------------------------------------------------------------
!> \brief get inversion grid contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \return invgrid inversion grid contained in this
!
  function getInversionGridIterationStepBasics(this) result(invgrid)
    type (iteration_step_basics), intent(in) :: this
    type (inversion_grid) :: invgrid
    invgrid = this%invgrid
  end function getInversionGridIterationStepBasics
!------------------------------------------------------------------------
!> \brief get integration weights contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \return intw integration weights contained in this
!
  function getIntegrationWeightsIterationStepBasics(this) result(intw)
    type (iteration_step_basics), intent(in) :: this
    type (integration_weights) :: intw
    intw = this%intw
  end function getIntegrationWeightsIterationStepBasics
!------------------------------------------------------------------------
!> \brief get kernel reference model contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \param krm kernel_reference_mode object
!! \return kernel reference model contained in this
!
  function getKernelReferenceModelIterationStepBasics(this) result(krm)
    type (iteration_step_basics), intent(in) :: this
    type (kernel_reference_model) :: krm
    krm = this%krm
  end function getKernelReferenceModelIterationStepBasics
!------------------------------------------------------------------------
!> \brief get frequency indices contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \param ifreq frequency indices
!! \return frequency indices contained in this
!
  function getFrequencyIndicesIterationStepBasics(this) result(ifreq)
    type (iteration_step_basics), intent(in) :: this
    integer, dimension(:), pointer :: ifreq
    ifreq => this%ifreq
  end function getFrequencyIndicesIterationStepBasics
!------------------------------------------------------------------------
!> \brief get number of frequency indices contained in iteration_step_basics object
!! \param this iteration_step_basics object
!! \param nf number of frequency indices (size(this%ifreq))
!! \return number of frequency indices contained in this
!
  function getNumberOfFrequenciesIterationStepBasics(this) result(nf)
    type (iteration_step_basics), intent(in) :: this
    integer :: nf
    if(associated(this%ifreq)) then
       nf = size(this%ifreq)
    else
       nf = 0
    end if
  end function getNumberOfFrequenciesIterationStepBasics
!
end module iterationStepBasics
