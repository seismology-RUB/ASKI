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
program computeMisfit
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use kernelLinearSystem
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=13) :: myname = 'computeMisfit'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  type (kernel_linear_system) :: KLSE

  real, dimension(:), pointer :: mdata,sdata

  integer, dimension(:), pointer :: idx,ifreq
  integer, dimension(:), pointer :: ifreq_all
  integer :: nf,j

  real :: misfit

  integer :: ios

  logical :: restrict_to_frequencies

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"reads in measured and synthetic data characterized by data model space info file and "//&
       "computes the data misfit")
  call addPosarg(ap,'dmspace_file','sval',"data model space input file which defines dataspace for which the "//&
       "data misfit is calculated")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-jf',.true.,"vector of frequency indices. If set, additionally to the misfit of the whole "//&
       "dataset defined by dmspace_file, the misfit is computed for data subsets restricted to these "//&
       "individual frequencies",'ivec','')
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
  dmspace_file = ap.sval.'dmspace_file'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  ! is jf set?
  restrict_to_frequencies = ap.optset.'-jf'
  if(restrict_to_frequencies) then
     ifreq_all => ap.ivec.'-jf'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(ifreq_all)) then
        write(*,*) "ERROR: for some reason, there is no list of frequency indices returned by argument parser, "//&
             "even though there was no error parsing argument -jf. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
  end if
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(main_parfile),get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
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
  write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"
!
!------------------------------------------------------------------------
!  read in measured and synthetic data, compute difference residual and misfit
!
  write(*,*) "reading in measured data now"
  call new(errmsg,myname)
  call initiateSerialKernelLinearSystem(KLSE,dmspace,0,0,errmsg)
  if (.level.errmsg == 2) then
     call print(errmsg)
     goto 1
  end if
  !subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !  path_measured_data,lu,dmspace,errmsg)
  call readMeasuredDataSerialKernelLinearSystem(KLSE,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg == 2) then
     call print(errmsg)
     goto 1
  end if
!
  write(*,*) "reading in synthetic data now"
  !subroutine readSyntheticsSerialKernelLinearSystem(this,path_event_filter,path_station_filter,&
  !  nfreq_measured_data,ifreq_measured_data,nfreq_synthetic_data,ifreq_synthetic_data,&
  !  path_synthetic_data,comptrans,lu,dmspace,errmsg)
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
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  compute misfit(s)
!
   write(*,*) ""
   write(*,*) "the misfits are computed as sum of squares of components of the residual data vector "
   write(*,*) "(measured minus synthetic data), dependent on the data(sub)set"
   write(*,*) ""
!
   ! get vector of synthetic and measured data samples
   sdata => .sd.KLSE
   mdata => .md.KLSE
!
   ! total misfit, i.e. whole dataset
   misfit = getMisfitKernelLinearSystem(KLSE,iostat=ios)
   if(ios/=0) then
      write(*,*) "MISFIT OF COMPLETE DATASET:    "//"error computing misfit, raised iostat = ",ios
      goto 1
   end if
   write(*,*) "MISFIT OF COMPLETE DATASET:    ",misfit,"  NORMALIZED: ",misfit/sum(mdata**2)
!
   ! in case of data subsets restricted to individual frequencies, compute individual misfits
   if(restrict_to_frequencies) then
      do j = 1,nf
         if(associated(ifreq)) deallocate(ifreq)
         allocate(ifreq(1)); ifreq(1) = ifreq_all(j)
         if(associated(idx)) deallocate(idx)
!
         idx => getIndxDataSamples(dmspace,ifreq=ifreq)
!
         if(.not.associated(idx)) then
            write(*,*) "MISFIT OF DATA SUBSET FOR jf =",ifreq_all(j)," : "//"requested frequency index not contained in dataset"
         else
            misfit = sum( (mdata(idx)-sdata(idx))**2 )
            write(*,*) "MISFIT OF DATA SUBSET FOR jf =",ifreq_all(j)," : ",misfit,"  NORMALIZED: ",misfit/sum(mdata(idx)**2)
         end if
!
      end do ! j
      if(associated(ifreq)) deallocate(ifreq)
      if(associated(idx)) deallocate(idx)
   end if
!
!------------------------------------------------------------------------
!  clean up
!
   write(*,*) "good bye"
!
1  call dealloc(KLSE)
   call dealloc(dmspace)
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)
!
   if(associated(ifreq_all)) deallocate(ifreq_all)
!
end program computeMisfit
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'    computeMisfit [-h] [-jf "nf jf1..jfn"] dmspace_file parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    dmspace_file: data model space input file which defines data and model space"
!   print *,"    parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print *,'-jf    : "nf jf1..jfn" defines a set of nf frequency indices. Additionally to the misfit of the whole dataset '
!   print *,'         defined by dmspace_file, the misfit is computed for data subsets restricted to these individual frequencies'
!   print *,''
!   print '(50(1h-))'
!   return
! end subroutine printhelp
