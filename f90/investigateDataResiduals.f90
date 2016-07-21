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
program investigateDataResiduals

  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use kernelLinearSystem
  use eventStationVtkFile
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,outdir_stats_files,outfile_vtk,vtk_filebase,outfile

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=24) :: myname = 'investigateDataResiduals'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace,dmspace_path_comp

  type (kernel_linear_system) :: KLSE
  real, dimension(:), pointer :: mdata,sdata,csdata
  real, dimension(:), allocatable :: res,res_norm,cres,cres_norm

  character(len=character_length_component), dimension(:), pointer :: pcomp,all_comp
  character(len=2), dimension(:), pointer :: pimre
  character(len=character_length_staname) :: staname
  character(len=character_length_staname), dimension(:), pointer :: pstaname
  character(len=character_length_evid) :: evid
  character(len=character_length_evid), dimension(:), pointer :: pevid
  character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths,paths_comp_jf
  integer, dimension(:), pointer :: pifreq,ifreq,idx_data_path,idx_data_comp_jf
  
  real, dimension(:,:), allocatable :: rarray2d

  ! VTK OUTPUT VARIABELS
  type (event_station_vtk_file) :: evstat_vtk,vtk_res_imre,vtk_res_ap
  integer, dimension(:,:), allocatable :: idx_data
  real, dimension(:), allocatable :: rdata_vtk_tmp1,rdata_vtk_tmp2
  complex, dimension(:), allocatable :: cdata_vtk
  character(len=200) :: vtk_title

  logical :: produce_vtk_files,stop_program,next
  real :: df
  integer :: i,icomp,jf,ipath,nf,npath_comp_jf,lu
!
!------------------------------------------------------------------------
!  preliminary processing
!
  stop_program = .false.
!
  call init(ap,myname,'Get some statistics about the data residuals of a given data (sub)set')
  ! define optional arguments
  ! call addOption(ap,'-smooth',.false.,&
  !      'indicates if linear smoothing constraints (average neighbour) are added to the kernel linear system')
  ! call addOption(ap,'-scltyp',.true.,&
  !      "type of scaling of smoothing constraints, at the moment only 'absmax_per_param,overall_factor' "//&
  !      "(requires one value of -sclval) allowed.",'sval','absmax_per_param,overall_factor')
  ! define positional arguments
  call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
  call addPosarg(ap,'outdir_stats_files','sval','output directory where residual files per path and component will written')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-ovtk',.true.,"base name of output vtk files. If not set, no vtk files will be produced",'sval','dataset')
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if

  ! get values of positional arguments
  main_parfile = ap.sval.'main_parfile'
  dmspace_file = ap.sval.'dmsi_file'
  outdir_stats_files = ap.sval.'outdir_stats_files'

  ! check if vtk files should be produced
  produce_vtk_files = ap.optset.'-ovtk'
  if(produce_vtk_files) outfile_vtk = ap.sval.'-ovtk'

  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
!
  df = (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  nf = .nf.iterbasics
  ifreq => .ifreq.iterbasics
!
!------------------------------------------------------------------------
!  setup data space info object
!
  write(*,*) "creating data space info from file '"//trim(dmspace_file)//"'"
!
  call new(errmsg,myname)
  call createDataSamplesFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  write(*,*) "there are ",.ndata.dmspace," data samples from ",.npath.dmspace," paths"
  write(*,*) ""
!
!------------------------------------------------------------------------
!  for all paths and components contained in the data space, produce statistics/residual files
!
  write(*,*) "PRODUCING FILES CONTAINING DATA RESIDUAL STATISTICS NOW IN DIRECTORY"
  write(*,*) "'",trim(outdir_stats_files),"'"
  write(*,*) ""
  call produceDataResidualsFiles()
  if(stop_program) goto 1
!
!------------------------------------------------------------------------
!  produce the vtk files, if requested
!
  if(produce_vtk_files) then
     write(*,*) "PRODUCING VTK FILES NOW"
     write(*,*) ""
     call produceDataResidualsVtkFiles()
     if(stop_program) goto 1
  end if
!
!------------------------------------------------------------------------
!  clean up and stop program
!
1 call dealloc(errmsg)
  call dealloc(ap)
  call dealloc(fuh)
  call dealloc(iterbasics)
  call dealloc(invbasics)
  call dealloc(dmspace)
!
  stop
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
contains
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine produceDataResidualsFiles()
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
     allocate(rarray2d(nf,10),cres(2*nf))
  else
     allocate(rarray2d(nf,5))
  end if
  allocate(idx_data(nf,2),res(2*nf))

  ! write textual file explaining the content of all files which will be written in the following
  lu = get(fuh)
  open(unit=lu,file = trim(outdir_stats_files)//"/content_info.txt", form = 'formatted', status = 'unknown', action ='write')
  write(lu,*) "the files 'residuals_evid_staname_comp' in this directory '"//trim(outdir_stats_files)//"'"
  write(lu,*) "contain residual quantities for all paths taken from data model space info file"
  write(lu,*) "'"//trim(dmspace_file)//"'"
  write(lu,*) "and all ",nf," frequencies used in current iteration ",(.inpar.invbasics).ival.'CURRENT_ITERATION_STEP'," :"
  write(lu,*) "df = ",df,"; frequency indices jf (such that frequency = jf*df) = "
  write(lu,*) ifreq
  write(lu,*) ""
  write(lu,*) "each row corresponds to a different frequency (in order of these frequency indices)"
  write(lu,*) "each column corresponds to a different residual quantity:"
  write(lu,*) "column 1: real part of residual 'data - synthetics'"
  write(lu,*) "column 2: imaginary part of residual 'data - synthetics'"
  write(lu,*) "column 3: absolute value of complex residual 'data - synthetics' (i.e. amplitude)"
  write(lu,*) "column 4: absolute value of complex residual 'data - synthetics' divided by absolute "//&
       "value of complex data (i.e. normalized amplitude)"
  write(lu,*) "column 5: atan2 value of complex residual 'data - synthetics' (i.e. phase)"
  write(lu,*) "columns 6-10: the same quantities as in colums 1-5 for the corrected residuals 'data - synthetics - "//&
       "synth_corrections' (in case of path specific models)"
  write(lu,*) ""
  close(lu)
  call undo(fuh)

  ! loop on all paths
  do while(nextPathDataModelSpaceInfo(dmspace,evid,staname,indx=idx_data_path,all_comp=all_comp))
     ! at this point, there is a next path for which idx_data_path and all_comp are associated 
     ! (otherwise module dataModelSpaceInfo is corrupt)

     ! subloop on all components of this path
     do icomp = 1,size(all_comp)
        
        ! create temporary data space object of this evid,staname,component and ALL frequencies of this iteration and BOTH im/re
        call addDataSamplesDataModelSpaceInfo(dmspace_path_comp,(/evid/),(/staname/),&
             (/all_comp(icomp)/),ifreq,(/'im','re'/),all_combinations=.true.)
        if(.ndata.dmspace_path_comp /= 2*nf) then
           write(*,*) "ERROR: temporary data space for one path and one component contains ",.ndata.dmspace_path_comp,&
                "data samples, expected to be ",2*nf,". Module dataModelSpaceInfo is corrupt"
           goto 3
        end if

        ! read in measured and synthetic (and corrections) in temporary KLSE object
        call new(errmsg,myname)
        call initiateSerialKernelLinearSystem(KLSE,dmspace_path_comp,0,0,errmsg)
        if (.level.errmsg == 2) then
           call print(errmsg)
           goto 3
        end if
        ! re-use same error message for reading in measured and synthetic data (more comprehensive if an error ocurrs)
        call readMeasuredDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg)
        call undo(fuh)
        if (.level.errmsg == 2) then
           call print(errmsg)
           goto 3
        end if
        mdata => .md.KLSE
        call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
             ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
             trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
             apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
             path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
             apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
             path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'))
        call undo(fuh)
        if (.level.errmsg == 2) then
           call print(errmsg)
           goto 3
        end if
        sdata => .sd.KLSE
        res = mdata-sdata
        ! IN CASE OF PATH SPECIFIC MODELS, ACCOUNT FOR SYNTHETIC CORRECTIONS
        if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
           !subroutine readSyntheticDataSerialKernelLinearSystem(this,path_event_filter,path_station_filter,&
           !  nfreq_measured_data,ifreq_measured_data,nfreq_synthetic_data,ifreq_synthetic_data,&
           !  path_synthetic_data,comptrans,lu,dmspace,errmsg)
           call new(errmsg,myname)
           call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
                ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
                ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
                ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
                trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
                apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
                path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
                apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
                path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),read_synthetic_corrections=.true.)
           call undo(fuh)
           if (.level.errmsg /= 0) call print(errmsg)
           !call print(errmsg)
           if (.level.errmsg == 2) goto 3
           csdata => .scd.KLSE
           cres = mdata-sdata-csdata
        else
           if (.level.errmsg /= 0) call print(errmsg)
        end if ! USE_PATH_SPECIFIC_MODELS
        call dealloc(errmsg)

        ! get index mappings for correct order of frequencies (as defined by iter basics) 
        do jf=1,nf
           if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
           if(associated(pifreq)) deallocate(pifreq)
           if(associated(pimre)) deallocate(pimre)
           allocate(pifreq(1)); pifreq(1) = ifreq(jf)
           idx_data_comp_jf => getIndxDataSamples(dmspace_path_comp,imre=pimre,ifreq=pifreq)
           if(.not.associated(idx_data_comp_jf)) then
              write(*,*) "ERROR: no data samples returned for component '",trim(all_comp(icomp)),&
                   "', freq. indx ",ifreq(jf)," and path (evid,staname) = ('",trim(evid),&
                   "','",trim(staname),"'), although there should be some; ",&
                   "this ERROR should actually not occurr!! please contact the developer"
              goto 3
           end if

           ! there should be exactly two different data samples now, one with imre = "im", the other with "re"
           ! if not, raise error (for the future: could support only one of them being present, in this case
           ! only imre residuals can be written)
           if(size(pimre)/= 2 .or. (.not.any(pimre=='im')) .or. (.not.any(pimre=='re'))) then
              write(*,*) "ERROR: inconsistent data set: for component '",trim(all_comp(icomp)),"', freq. indx ",&
                   ifreq(jf)," and path (evid,staname) = ('",trim(evid),"','",trim(staname),&
                   "') there are ",size(pimre)," data samples (must be 2) with im/re = ","'"//pimre//"'"
              goto 3
           end if
           do i = 1,2
              if(pimre(i) == 're') then
                 idx_data(jf,1) = idx_data_comp_jf(i)
              else
                 idx_data(jf,2) = idx_data_comp_jf(i)
              end if
           end do
        end do ! jf

        ! compute the required quantities, write in a real matrix (columnwise, each row corresponds to the correct frequency) 
       ! DATA RESIDUALS, REAL PART
        rarray2d(:,1) = res(idx_data(:,1))
        ! DATA RESIDUALS, IMAGINARY PART
        rarray2d(:,2) = res(idx_data(:,2))
        ! DATA RESIDUALS, AMPLITUDE
        rarray2d(:,3) = sqrt(res(idx_data(:,1))*res(idx_data(:,1)) + res(idx_data(:,2))*res(idx_data(:,2)))
        ! DATA RESIDUALS, NORMALIZED AMPLITUDE
        rarray2d(:,4) = rarray2d(:,3)/&
             sqrt(mdata(idx_data(:,1))*mdata(idx_data(:,1)) + mdata(idx_data(:,2))*mdata(idx_data(:,2)))
        ! DATA RESIDUALS, PHASE
        rarray2d(:,5) = atan2(res(idx_data(:,2)),res(idx_data(:,1)))

        ! CORRECTED QUANTITIES, (if path specific)
        if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
           ! DATA RESIDUALS, REAL PART
           rarray2d(:,6) = cres(idx_data(:,1))
           ! DATA RESIDUALS, IMAGINARY PART
           rarray2d(:,7) = cres(idx_data(:,2))
           ! DATA RESIDUALS, AMPLITUDE
           rarray2d(:,8) = sqrt(cres(idx_data(:,1))*cres(idx_data(:,1)) + cres(idx_data(:,2))*cres(idx_data(:,2)))
           ! DATA RESIDUALS, NORMALIZED AMPLITUDE
           rarray2d(:,9) = rarray2d(:,8)/&
                sqrt(mdata(idx_data(:,1))*mdata(idx_data(:,1)) + mdata(idx_data(:,2))*mdata(idx_data(:,2)))
           ! DATA RESIDUALS, PHASE
           rarray2d(:,10) = atan2(cres(idx_data(:,2)),cres(idx_data(:,1)))
        end if

        ! write by module asciiDataIO to output dir, give naming by residuals_evid_staname_comp
        outfile = trim(outdir_stats_files)//"/residuals_"//trim(evid)//"_"//trim(staname)//"_"//&
             trim(trim(all_comp(icomp)))
        errmsg = writeAsciiData(outfile,get(fuh),rarray2d)
        call undo(fuh)
        if(.level.errmsg /= 0) call print(errmsg)
        if(.level.errmsg == 2) goto 3
        call dealloc(errmsg)

        call dealloc(KLSE)
        call dealloc(dmspace_path_comp)
     end do ! icomp

  end do ! while(nextPath)

1 deallocate(rarray2d,idx_data,res)
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) deallocate(cres)
  call dealloc(errmsg)
  call dealloc(dmspace_path_comp)
  call dealloc(KLSE)
  if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
  if(associated(pifreq)) deallocate(pifreq)
  if(associated(pimre)) deallocate(pimre)
!
  ! if code comes here, return normally
  return
!
  ! in case there was some error, return and indicate to stop the main program
2 stop_program = .true.
  goto 1
!
  ! in case there was an error inside the loop on all paths, reset the iterator first (deallocating all arrays)
3 next = nextPathDataModelSpaceInfo(dmspace,evid,staname,indx=idx_data_path,all_comp=all_comp,reset=.true.)
  goto 2
end subroutine produceDataResidualsFiles
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine produceDataResidualsVtkFiles()
!------------------------------------------------------------------------
!  read in measured and synthetic data (and synthetic corrections data)
!
  write(*,*) "reading in measured data now"
  !subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !  path_measured_data,lu,dmspace,errmsg)
  call new(errmsg,myname)
  call initiateSerialKernelLinearSystem(KLSE,dmspace,0,0,errmsg)
  if (.level.errmsg == 2) then
     call print(errmsg)
     goto 2
  end if
  ! re-use same error message for reading in measured and synthetic data (more comprehensive if an error ocurrs)
  call readMeasuredDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg == 2) then
     call print(errmsg)
     goto 2
  end if
  mdata => .md.KLSE
  write(*,*) ""
!
  write(*,*) "reading in synthetic data now"
  !subroutine readSyntheticsSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !  nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
  !  apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !  ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
  call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
       ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
       trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
       apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
       path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
       apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
       path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'))
  call undo(fuh)
  if (.level.errmsg == 2) then
     call print(errmsg)
     goto 2
  end if
  sdata => .sd.KLSE
  write(*,*) ""
!
  write(*,*) "computing (normalized) residuals now"
  allocate(res(.ndata.dmspace),res_norm(.ndata.dmspace))
  res = mdata-sdata
  res_norm = abs(res)/abs(mdata)
  write(*,*) ""
!
  ! IN CASE OF PATH SPECIFIC MODELS, ACCOUNT FOR SYNTHETIC CORRECTIONS
  if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
     write(*,*) "reading in synthetics correction data now"
     !subroutine readSyntheticsSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
     !  nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
     !  apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
     !  ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
     call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
          ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
          trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
          apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
          path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
          apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
          path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),read_synthetic_corrections=.true.)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 2
     write(*,*) ""
!
     write(*,*) "computing corrected (normalized) residuals now"
     csdata => .scd.KLSE
     allocate(cres(.ndata.dmspace),cres_norm(.ndata.dmspace))
     cres = mdata-sdata-csdata
     cres_norm = abs(cres)/abs(mdata)
     write(*,*) ""
  else
     if (.level.errmsg /= 0) call print(errmsg)
  end if ! USE_PATH_SPECIFIC_MODELS
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  and output the results
!
  ! OUTPUT PATH VTK FILES

  ! first of all output one vtk file (without data), containing all paths of the data model space

  ! get all the paths contained in the data space
  paths => getPathsDataModelSpaceInfo(dmspace)
  if(.not.associated(paths)) then
     write(*,*) "ERROR: no paths contained in data space, data space seems to be empty"
     goto 2
  end if

  call new(errmsg,myname)
  call init(evstat_vtk,.evlist.invbasics,.statlist.invbasics,paths,.invgrid.iterbasics,trim(outfile_vtk),&
       trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 2
  call writePaths(evstat_vtk,get(fuh),errmsg,overwrite =.false.)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  if (.level.errmsg == 2) goto 2
  call dealloc(errmsg)
  if(associated(paths)) deallocate(paths)
!
  ! then output vtk files for all components and all frequencies, each containing all data samples (i.e. paths and im/re) for each combination of component and frequency
!
  all_comp => allComp(dmspace)
  if(.not.associated(all_comp)) then
     write(*,*) "ERROR: no components (i.e. no data samples) contained in data model space info"
     goto 2
  end if
!
  do icomp = 1,size(all_comp)
     do jf = 1,nf

        ! FOR EVERY COMPONENT AND FREQUENCY, THERE WILL BE ONE VTK FILE (one file for each of 
        ! (amplitude/phase) residuals, (amplitude/phase) normalized residuals, ...)

        if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
        if(associated(pcomp)) deallocate(pcomp)
        if(associated(pifreq)) deallocate(pifreq)
        allocate(pcomp(1),pifreq(1))
        pcomp(1) = all_comp(icomp)
        pifreq(1) = ifreq(jf)
!
        idx_data_comp_jf => getIndxDataSamples(dmspace,comp=pcomp,ifreq=pifreq)
!
        ! now get all the paths which hold data samples of this component and this frequency
        ! (in general NOT ALL paths, hence, we need to construct the vtk files separately inside this loop)
        ! loop on all paths and get indices of im and re
!
        ! get all the paths for which the current component
        if(associated(paths_comp_jf)) deallocate(paths_comp_jf)
        paths_comp_jf => getPathsDataModelSpaceInfo(dmspace,idata_in=idx_data_comp_jf)
        if(.not.associated(paths_comp_jf)) then
           write(*,*) "ERROR: no paths returned for component '",trim(all_comp(icomp)),&
                "' and frequency index ",ifreq(jf)," although there should be",&
                "; this ERROR should actually not occurr!! please contact the developer"
           if(.not.associated(idx_data_comp_jf)) then
              write(*,*) "idx_data_comp_jf was already not associated"
           else
              write(*,*) "there were ",size(idx_data_comp_jf)," data samples in which there were paths to be searched"
           end if
           goto 2
        end if
        npath_comp_jf = size(paths_comp_jf,2)
!
        ! initiate the vtk files now
        vtk_filebase = trim(outfile_vtk)//"_residuals_im-re_"//trim(all_comp(icomp))
        write(*,*) "writing data residuals (imaginary/real parts) on paths for component '",&
             trim(all_comp(icomp)),"' and frequency index ",ifreq(jf)

        write(vtk_title,*) "data residuals (imaginary/real parts) at frequency ",ifreq(jf)*df," Hz on paths"
        call new(errmsg,myname)
        call init(vtk_res_imre,.evlist.invbasics,.statlist.invbasics,paths_comp_jf,.invgrid.iterbasics,vtk_filebase,&
             trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title)
        if (.level.errmsg /= 0) call print(errmsg)
        if (.level.errmsg == 2) goto 2
        call dealloc(errmsg)
!
        vtk_filebase = trim(outfile_vtk)//"_residuals_amp-phase_"//trim(all_comp(icomp))
        write(*,*) "writing data residuals (amplitudes and pases) on paths for component '",&
             trim(all_comp(icomp)),"' and frequency index ",ifreq(jf)

        write(vtk_title,*) "data residuals (amplitudes and pases) at frequency ",ifreq(jf)*df," Hz on paths"
        call new(errmsg,myname)
        call init(vtk_res_ap,.evlist.invbasics,.statlist.invbasics,paths_comp_jf,.invgrid.iterbasics,vtk_filebase,&
             trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title)
        if (.level.errmsg /= 0) call print(errmsg)
        if (.level.errmsg == 2) goto 2
        call dealloc(errmsg)
!
        ! .. initiate more vtk files for other quantities (a bit lengthy, but cannot be realized otherwise)
!
        ! collect index mappings for the the data of the vtk files
        allocate(idx_data(npath_comp_jf,2))
        do ipath = 1,npath_comp_jf
           if(associated(idx_data_path)) deallocate(idx_data_path)
           if(associated(pimre)) deallocate(pimre)
           if(associated(pstaname)) deallocate(pstaname)
           if(associated(pevid)) deallocate(pevid)
           allocate(pevid(1)); pevid(1) = paths_comp_jf(1,ipath)
           allocate(pstaname(1)); pstaname(1) = paths_comp_jf(2,ipath)
!
           idx_data_path => getIndxDataSamples(dmspace,staname=pstaname,evid=pevid,imre=pimre,indx_in=idx_data_comp_jf)
           if(.not.associated(idx_data_path)) then
              write(*,*) "ERROR: no data samples returned for component '",trim(all_comp(icomp)),&
                   "' and frequency index ",ifreq(jf)," and path (evid,staname) = ('",trim(paths_comp_jf(1,ipath)),&
                   "','",trim(paths_comp_jf(2,ipath)),"'), although there should be some; ",&
                   "this ERROR should actually not occurr!! please contact the developer"
              goto 2
           end if
!
           ! there should be exactly two different data samples now, one with imre = "im", the other with "re"
           ! if not, raise error (for the future: could support only one of them being present, in this case
           ! only imre residuals can be written)
           if(size(pimre)/= 2 .or. (.not.any(pimre=='im')) .or. (.not.any(pimre=='re'))) then
              write(*,*) "ERROR: inconsistent data set: for component '",trim(all_comp(icomp)),"', freq. indx ",&
                   ifreq(jf),", path (evid,staname) = ('",trim(paths_comp_jf(1,ipath)),"','",trim(paths_comp_jf(2,ipath)),&
                   "') there are ",size(pimre)," data samples (must be 2) with im/re = ","'"//pimre//"'"
              goto 2
           end if
           do i = 1,2
              if(pimre(i) == 're') then
                 idx_data(ipath,1) = idx_data_path(i)
              else
                 idx_data(ipath,2) = idx_data_path(i)
              end if
           end do
        end do ! ipath
!
        allocate(cdata_vtk(npath_comp_jf),rdata_vtk_tmp1(npath_comp_jf),rdata_vtk_tmp2(npath_comp_jf))
!
        ! write residuals imaginary/real part, need to convert to complex data array
        ! for proper writing to 2-component-vector float vtk files
        rdata_vtk_tmp1 = res(idx_data(:,1))
        rdata_vtk_tmp2 = res(idx_data(:,2))
        cdata_vtk(:) = cmplx( rdata_vtk_tmp1 , rdata_vtk_tmp2)
        call new(errmsg,myname)
        call writeData(vtk_res_imre,get(fuh),cdata_vtk,errmsg,"residuals_imre",file_index=ifreq(jf))
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        if (.level.errmsg == 2) goto 2
        call dealloc(errmsg)

        ! compute and write residuals in amplitudes / phases, need to convert to complex data array
        ! for proper writing to 2-component-vector float vtk files
        rdata_vtk_tmp1 = sqrt(res(idx_data(:,1))*res(idx_data(:,1)) + res(idx_data(:,2))*res(idx_data(:,2))) ! amplitude
        rdata_vtk_tmp2 = atan2(res(idx_data(:,2)),res(idx_data(:,1))) ! phase
        cdata_vtk(:) = cmplx( rdata_vtk_tmp1 , rdata_vtk_tmp2)
        call new(errmsg,myname)
        call writeData(vtk_res_ap,get(fuh),cdata_vtk,errmsg,"residuals_ampl_phase",file_index=ifreq(jf))
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        if (.level.errmsg == 2) goto 2
        call dealloc(errmsg)

        deallocate(idx_data,cdata_vtk,rdata_vtk_tmp1,rdata_vtk_tmp2)
     end do ! jf
  end do ! icomp

1 call dealloc(errmsg)
  call dealloc(KLSE)
  if(associated(paths)) deallocate(paths)
  if(associated(all_comp)) deallocate(all_comp)
  if(allocated(res)) deallocate(res)
  if(allocated(res_norm)) deallocate(res_norm)
  if(allocated(cres)) deallocate(cres)
  if(allocated(cres_norm)) deallocate(cres_norm)
  if(associated(idx_data_comp_jf)) deallocate(idx_data_comp_jf)
  if(associated(pcomp)) deallocate(pcomp)
  if(associated(pifreq)) deallocate(pifreq)
  if(associated(paths_comp_jf)) deallocate(paths_comp_jf)
  if(associated(idx_data_path)) deallocate(idx_data_path)
  if(associated(pimre)) deallocate(pimre)
  if(associated(pstaname)) deallocate(pstaname)
  if(associated(pevid)) deallocate(pevid)
  if(allocated(idx_data)) deallocate(idx_data)
  if(allocated(cdata_vtk)) deallocate(cdata_vtk)
  if(allocated(rdata_vtk_tmp1)) deallocate(rdata_vtk_tmp1)
  if(allocated(rdata_vtk_tmp2)) deallocate(rdata_vtk_tmp2)

  ! if code comes here, return normally
  return

  ! in case there was some error, return and indicate to stop the main program
2 stop_program = .true.
  goto 1

end subroutine produceDataResidualsVtkFiles

end program investigateDataResiduals
