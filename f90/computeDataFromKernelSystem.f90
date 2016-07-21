!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------

program computeDataFromKernelSystem

  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use kernelInvertedModel
  use kernelLinearSystem
  use asciiDataIO
  use argumentParser
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=400) :: main_parfile,dmspace_file,outdir,kim_file,outfile

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace,dmspace_path
  character(len=max_character_length_evid_staname), dimension(:,:), pointer :: paths
  character(len=character_length_component), dimension(:), pointer :: all_comp,pcomp
  character(len=2), dimension(:), pointer :: pimre
  integer, dimension(:), pointer :: pifreq,idx_data_path_comp,idx_data_path_comp_jf
  integer :: ipath,npath,icomp,ncomp,i

  type (kernel_inverted_model) :: kim_in,kim_krm
  real, dimension(:,:), allocatable :: model_vector
  real, dimension(:), pointer :: tmp

  type (kernel_linear_system) :: KLSE
  real, dimension(:), pointer :: sdata,scdata
  real, dimension(:,:), pointer :: rhs
  real, dimension(:), allocatable :: mdata
  integer :: lu1,lu2

  complex, dimension(:), allocatable :: cvec_mdata
  integer :: jf_mdata,nfreq_mdata
  integer, dimension(:), pointer :: ifreq_mdata

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=27) :: myname = 'computeDataFromKernelSystem'

  nullify(paths,all_comp,pcomp,pimre,pifreq,idx_data_path_comp,idx_data_path_comp_jf,tmp,sdata,scdata,rhs,ifreq_mdata)
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Compute the 'measured' data d in the equation  d-s(-sc) = K (m1 - mref)  for given model m1")
  ! define positional arguments
  call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
  call addPosarg(ap,'kim_file','sval','kernel-inverted-model file which contains model m1')
  call addPosarg(ap,'outdir_data','sval',"directory where the new 'measured' data files are written to")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  ! get values of positional arguments
  dmspace_file = ap.sval.'dmsi_file'
  kim_file = ap.sval.'kim_file'
  outdir = ap.sval.'outdir_data'
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) goto 1
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  write(*,*) "SETTING UP INVERSION BASICS, main parfile = '",trim(main_parfile),"'"
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,main_parfile,get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  write(*,*) ""
!
  write(*,*) "SETTING UP ITERATION STEP BASICS"
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  nfreq_mdata = (.inpar.invbasics).ival.'MEASURED_DATA_NUMBER_OF_FREQ'
  ifreq_mdata => .ifreq.invbasics
  allocate(cvec_mdata(nfreq_mdata))
  write(*,*) "there are ",nfreq_mdata,"  frequencies of measured data with indices ",ifreq_mdata
  write(*,*) ""
!
!------------------------------------------------------------------------
!  setup data model space info object
!
  write(*,*) "SETTING UP DATA MODEL SPACE, dmspace_file = '",trim(dmspace_file),"'"
  call new(errmsg,myname)
  call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
       .ifreq.iterbasics,sval(.inpar.invbasics,'MODEL_PARAMETRIZATION'),&
       .ncell.(.invgrid.iterbasics),.intw.iterbasics,&
       trim(dmspace_file),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) "created data model space containing ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"
  paths => getPathsDataModelSpaceInfo(dmspace)
  if(.not.associated(paths)) then
     write(*,*) "ERROR: no paths contained in data space, data space seems to be empty"
     goto 1
  end if
  npath = size(paths,2)
  write(*,*) "data space contains ",npath," data paths"
  write(*,*) ""
!
!------------------------------------------------------------------------
!  read in model values for which the data should be computed (kernel-inverted-model object)
!
  write(*,*) "READ IN MODEL VALUES FOR WHICH THE DATA SHOULD BE COMPUTED, kim-object from file '",trim(kim_file),"'"
  call new(errmsg,myname)
  call readFileKernelInvertedModel(kim_in,kim_file,get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  write(*,*) "INTERPOLATING (GLOBAL) KERNEL REFERENCE MODEL TO INVERSION GRID AND SUBSTRACTING FROM THE REQUESTED KIM-MODEL"
!  subroutine interpolateKernelReferenceToKernelInvertedModel(this,krm,parametrization,invgrid,intw,errmsg)
  call new(errmsg,myname)
  call interpolateKernelReferenceToKernelInvertedModel(kim_krm,.krm.iterbasics,.pmtrz.dmspace,.invgrid.iterbasics,&
       .intw.iterbasics,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!  subroutine summateInstancesKernelInvertedModel(this,that,errmsg,c1,c2,relative)
  call new(errmsg,myname)
  call summateInstancesKernelInvertedModel(kim_in,kim_krm,errmsg,c1=1.0,c2=-1.0)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
  write(*,*) ""
!
  call new(errmsg,myname)
  call unpackToVectorKernelInvertedModel(tmp,kim_in,dmspace,errmsg)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  if(.not.associated(tmp)) then
     write(*,*) "ERROR: model vector returned from routine 'unpackToVectorKernelInvertedModel' is empty"
     goto 1
  end if
  if(size(tmp) /= .nmval.dmspace) then
     write(*,*) "ERROR: model vector returned from routine 'unpackToVectorKernelInvertedModel' has size ",size(tmp),&
          ", but requested values (number of model values in model space) is ",.nmval.dmspace
     goto 1
  end if
  allocate(model_vector(size(tmp),1))
  model_vector(:,1) = tmp
  deallocate(tmp)
!
  write(*,*) "LOOPING ON ALL PAHTS NOW"
  write(*,*) ""
! LOOP ON PATHS
  do ipath = 1,npath
     write(*,*) "path ",ipath,"; event '",trim(paths(1,ipath)),"' ,  station '",trim(paths(2,ipath)),"'"
!
     call copyDataModelSpaceInfo(dmspace_path,dmspace,ipath1=ipath,ipath2=ipath)
!
!  subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
     call new(errmsg,myname)
     call initiateSerialKernelLinearSystem(KLSE,dmspace_path,0,0,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
!
     write(*,*) "   read kernel matrix for this path"
!  subroutine readMatrixSerialKernelLinearSystem(this,df_measured_data,&
!       nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,ntot_invgrid,pcorr,&
!       lu1,lu2,errmsg,apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
!       ignore_data_weights,apply_kernel_normalization)
     lu1 = get(fuh); lu2 = get(fuh)
     call new(errmsg,myname)
     call readMatrixSerialKernelLinearSystem(KLSE,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
          ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
          .ncell.(.invgrid.iterbasics),.pcorr.invbasics,lu1,lu2,errmsg,&
          apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
          path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
          apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
          path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
          ignore_data_weights=.true.)
     call add(fuh,lu1); call add(fuh,lu2)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
!
!  subroutine setSolKernelLinearSystem(this,transpose,sol,errmsg)
     call new(errmsg,myname)
     call setSolKernelLinearSystem(KLSE,'N',model_vector,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
!
     write(*,*) "   multiply kernel matrix of this path with model (update) vector"
!  subroutine multiplyForwardKernelLinearSystem(this,errmsg)
     call new(errmsg,myname)
     call multiplyForwardKernelLinearSystem(KLSE,errmsg)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
     rhs => .rhs.KLSE
!
     call new(errmsg,myname)
     call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
          ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
          trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
          apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
          path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
          apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
          path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
          ignore_data_weights=.true.)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     call dealloc(errmsg)
     sdata => .sd.KLSE
!
     allocate(mdata(.ndata.dmspace_path))
     mdata = rhs(:,1) + sdata
!
     if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
        call new(errmsg,myname)
        call readSyntheticDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
             ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
             ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ'),&
             ivec(.inpar.iterbasics,'ITERATION_STEP_INDEX_OF_FREQ',ival(.inpar.iterbasics,'ITERATION_STEP_NUMBER_OF_FREQ')),&
             trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SYNTHETIC_DATA'),get(fuh),errmsg,&
             apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
             path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
             apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
             path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'),&
             ignore_data_weights=.true.,read_synthetic_corrections=.true.)
        call undo(fuh)
        if (.level.errmsg /= 0) call print(errmsg)
        !call print(errmsg)
        if (.level.errmsg == 2) goto 1
        call dealloc(errmsg)
        scdata => .scd.KLSE
        mdata = mdata + scdata
     end if ! USE_PATH_SPECIFIC_MODELS
!
!  computed mdata as rhs + sdata (+ scdata , if path specific)
!  write mdata to files (select components, frequencies; sort frequencies (set non existing ones to zero)
!
     if(lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')) then
        write(*,*) "   writing now to files:  'data' =  K * (kim_in - krm) + synthetics + synthetic_corrections"
     else
        write(*,*) "   writing now to files:  'data' =  K * (kim_in - krm) + synthetics"
     end if
!
     if(associated(all_comp)) deallocate(all_comp)
     all_comp => allComp(dmspace_path)
     if(.not.associated(all_comp)) then
        write(*,*) "ERROR: no components returned for path (evid,staname) = ('",trim(paths(1,ipath)),&
             "','",trim(paths(2,ipath)),"'), although there should be some; ",&
             "this ERROR should actually not occurr!! please contact the developer"
        goto 1
     end if
     ncomp = size(all_comp)
!
     write(*,*) "   looping on all ",ncomp," component(s) of this path:"
     ! subloop on all components of this path
     do icomp = 1,ncomp
        write(*,*) "      component '",trim(all_comp(icomp)),"'"
        if(associated(idx_data_path_comp)) deallocate(idx_data_path_comp)
        if(associated(pcomp)) deallocate(pcomp)
        allocate(pcomp(1)); pcomp(1) = all_comp(icomp)
        idx_data_path_comp => getIndxDataSamples(dmspace_path,comp=pcomp)
        if(.not.associated(idx_data_path_comp)) then
           write(*,*) "ERROR: no data samples returned for component '",trim(all_comp(icomp)),&
                "', although there should be some; this ERROR should actually not occurr!! please contact the developer"
           goto 1
        end if
!
        do jf_mdata = 1,nfreq_mdata
           cvec_mdata(jf_mdata) = (0.0,0.0)
!
           if(associated(idx_data_path_comp_jf)) deallocate(idx_data_path_comp_jf)
           if(associated(pifreq)) deallocate(pifreq)
           if(associated(pimre)) deallocate(pimre)
           allocate(pifreq(1)); pifreq(1) = ifreq_mdata(jf_mdata)
           idx_data_path_comp_jf => getIndxDataSamples(dmspace_path,imre=pimre,ifreq=pifreq,indx_in=idx_data_path_comp)
           if(.not.associated(idx_data_path_comp_jf)) then
              write(*,*) "WARNING: no data samples returned for measured data frequency index ",ifreq_mdata(jf_mdata),&
                   ", hence setting the data values in the output file in line ",jf_mdata," to (0.0,0.0)"
              cycle
           end if
           ! there should be exactly two different data samples now, one with imre = "im", the other with "re"
           ! if not, raise warning that the missing data sample is set to (real) 0.0
           select case(size(pimre))
           case (1)
              write(*,*) "WARNING: there is only one data sample for frequency index ",ifreq_mdata(jf_mdata)," with im/re = ",&
                   "'"//pimre//"'; setting the data value of the other one to 0.0"
           case (2)
              if(all(pimre=='im') .or. all(pimre=='re')) then
                 write(*,*) "ERROR: inconsistent data set: there are 2 data samples with im/re = ",&
                      "'"//pimre//"'; must be of different kind"
              end if
           case default
              write(*,*) "ERROR: inconsistent data set: there are ",size(pimre)," data samples for frequency index ",&
                   ifreq_mdata(jf_mdata)," with im/re = ","'"//pimre//"'; must be maximally 2 and all must be of different kind"
              goto 1
           end select
           do i = 1,size(pimre)
              if(pimre(i) == 're') then
                 cvec_mdata(jf_mdata) = cvec_mdata(jf_mdata) + cmplx(mdata(idx_data_path_comp_jf(i)),0.0)
              else
                 cvec_mdata(jf_mdata) = cvec_mdata(jf_mdata) + cmplx(0.0,mdata(idx_data_path_comp_jf(i)))
              end if
           end do ! i
        end do ! jf_mdata
!
        ! write complex vector of data values by module asciiDataIO to output dir, give naming by convention for measured data
        outfile = trim(outdir)//"/data_"//trim(paths(1,ipath))//"_"//trim(paths(2,ipath))//"_"//&
             trim(trim(all_comp(icomp)))
        write(*,*) "      writing output file '",trim(outfile),"'"
        errmsg = writeAsciiData(outfile,get(fuh),cvec_mdata)
        call undo(fuh)
        if(.level.errmsg /= 0) call print(errmsg)
        if(.level.errmsg == 2) goto 1
        call dealloc(errmsg)
!
     end do ! icomp
!
     deallocate(mdata)
     call dealloc(dmspace_path)
     call dealloc(KLSE)
  end do ! ipath
!  END OF LOOP ON PAHTS
!
  write(*,*) ""; write(*,*) "good bye"; write(*,*) ""
1 call dealloc(iterbasics); call dealloc(invbasics)
  call dealloc(dmspace); call dealloc(dmspace_path)
  call dealloc(KLSE)
  call dealloc(kim_in); call dealloc(kim_krm)
  call dealloc(ap)
  if(allocated(model_vector)) deallocate(model_vector)
  if(allocated(mdata)) deallocate(mdata)
  if(associated(paths)) deallocate(paths)
  if(associated(all_comp)) deallocate(all_comp)
  if(associated(idx_data_path_comp)) deallocate(idx_data_path_comp)
  if(associated(idx_data_path_comp_jf)) deallocate(idx_data_path_comp_jf)
  if(associated(pcomp)) deallocate(pcomp)
  if(associated(pifreq)) deallocate(pifreq)
  if(associated(pimre)) deallocate(pimre)
  if(allocated(cvec_mdata)) deallocate(cvec_mdata)
!
end program computeDataFromKernelSystem
