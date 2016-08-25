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
!----------------------------------------------------------------------------
!## THIS IS A HARD-CODED PROGRAM WHICH IS PROBABLY SUITABLE FOR DEVELOPERS USE ONLY
!----------------------------------------------------------------------------
program computeIterationMisfits
  use inversionBasics
  use dataModelSpaceInfo
  use kernelLinearSystem
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: parfile,dmspace_file,sdata_subpath,sdata_path,&
       outfile_misfit,outfile_misfit_normalized

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=24) :: myname = 'computeIterqationMisfits'

  type (inversion_basics) :: invbasics
  !type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  type (kernel_linear_system) :: KLSE_mdata
  type (kernel_linear_system) :: KLSE_sdata

  real, dimension(:), pointer :: mdata,sdata

  integer, dimension(:), pointer :: idx,ifreq,ifreq_all
  integer, dimension(:), allocatable :: jiter
  integer :: nf,jf,niter,j

  real, dimension(:,:), allocatable :: iteration_misfits,iteration_misfits_normalized
  real :: misfit,sum_mdata_2

  logical :: restrict_to_frequencies

  nullify(mdata,sdata,idx,ifreq,ifreq_all)

!------------------------------------------------------------------------
!  HARDCODED PARAMETERS, PLEASE ADJUST HERE AS NECESSARY:
!
  ! iteration steps over which should be looped below
  niter = 13
  allocate(jiter(niter))
  jiter = (/1,2,3,4,5,6,7,8,9,10,11,12,13/)
!
  ! subpath of synthetic data in each iteration_step directory, in measured_data discretization (and measured_data file formats)
  sdata_subpath = 'synthetic_data_stf_rickr_50_T_0.25/'
!
  ! outfiles
  outfile_misfit = "/data/Kernel/misfit.dat"
  outfile_misfit_normalized = "/data/Kernel/misfit_normalized.dat"
!
!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"DEVELOPERS TOOL, CONTAINS HARD-CODED PARAMETERS, PLEASE ADJUST IN CODE FILE AND RE-COMPILE! "//&
       "computes the non-linear misfit for a set of iteration steps (provided there were additional forward "//&
       "calculations done with the new models)")
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
     stop
  end if
!
  parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     stop
  end if
  dmspace_file = ap.sval.'dmspace_file'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     stop
  end if
!
  ! -jf
  restrict_to_frequencies = ap.optset.'-jf'
  if(restrict_to_frequencies) then
     ifreq_all => ap.ivec.'-jf'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        stop
     end if
     if(.not.associated(ifreq_all)) then
        write(*,*) "ERROR: for some reason, there is no list of frequency indices returned by argument parser, "//&
             "even though there was no error parsing argument -jf. This is strange..."
        write(*,*) ""
        call usage(ap)
        stop
     end if
  else
     nf = 0
  end if
!
  call document(ap)
  write(*,*) ""
!
  ! creat file unit handler  
  call createFileUnitHandler(fuh,100)
!
  allocate(iteration_misfits(niter,nf+1))
  iteration_misfits = -1.0
  allocate(iteration_misfits_normalized(niter,nf+1))
  iteration_misfits_normalized = -1.0
!
!------------------------------------------------------------------------
!  setup basics
!
  ! setup inversion basics
  call new(errmsg,myname)
  call init(invbasics,trim(parfile),get(fuh),errmsg)
  call undo(fuh)
  call addTrace(errmsg,myname)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
!------------------------------------------------------------------------
!  setup data space info object
!
  write(*,*) "creating data space info from file '"//trim(dmspace_file)//"'"
!
  call new(errmsg,myname)
  call createDataSamplesFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
!!$       .ifreq.iterbasics,trim(dmspace_file),get(fuh),errmsg)
       .ifreq.invbasics,trim(dmspace_file),get(fuh),errmsg) !! FS FS TMP
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
  write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"
!
!------------------------------------------------------------------------
!  read in measured data
!
  write(*,*) "reading in measured data now"
  call new(errmsg,myname)
  call initiateSerialKernelLinearSystem(KLSE_mdata,dmspace,0,0,errmsg)
  if (.level.errmsg == 2) then
     call print(errmsg)
     stop
  end if
  !subroutine readMeasuredDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !  path_measured_data,lu,dmspace,errmsg)
  call readMeasuredDataSerialKernelLinearSystem(KLSE_mdata,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
       ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
       sval(.inpar.invbasics,'PATH_MEASURED_DATA'),get(fuh),errmsg)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
  mdata => .md.KLSE_mdata
  sum_mdata_2 = real( sum(dble(mdata)**2) )
!
   write(*,*) ""
   write(*,*) "the misfits are computed as sum of squares of components of the residual data vector "
   write(*,*) "(measured minus synthetic data), dependent on the data(sub)set(s)"
   write(*,*) ""
!
!------------------------------------------------------------------------
!  loop on all iteration steps indicated at the top, 
!  and synthetic data, compute difference residual and misfit
!
  do j = 1,niter
!
     write(*,*) "ITERATION ",jiter(j)
     write(*,*) ""
!
     write(sdata_path,"(2a,i3.3,2a)") trim((.inpar.invbasics).sval.'MAIN_PATH_INVERSION'),&
          trim((.inpar.invbasics).sval.'ITERATION_STEP_PATH'),jiter(j),'/',trim(sdata_subpath)
!
     ! read in synthetic data
!
     write(*,*) "reading in synthetic data of now, in the form of ASKI measured data from path ITER_PATH/"//&
          trim(sdata_subpath)//" :"
     write(*,*) "  '"//trim(sdata_path)//"'"
     call new(errmsg,myname)
     call initiateSerialKernelLinearSystem(KLSE_sdata,dmspace,0,0,errmsg)
     if (.level.errmsg == 2) then
        call print(errmsg)
        stop
     end if
     call readMeasuredDataSerialKernelLinearSystem(KLSE_sdata,ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
          ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
          trim(sdata_path),get(fuh),errmsg)
     call undo(fuh)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) stop
     call dealloc(errmsg)
     sdata => .md.KLSE_sdata
!
     !  compute misfit(s)
!
     ! total misfit, i.e. whole dataset
     misfit = real( sum( (dble(mdata)-dble(sdata))**2 ) )
     iteration_misfits(j,1) = misfit
     iteration_misfits_normalized(j,1) = real( dble(misfit)/dble(sum_mdata_2) )
     write(*,*) "MISFIT OF COMPLETE DATASET:    ",misfit,"  NORMALIZED: ",real( dble(misfit)/dble(sum_mdata_2) )
!
     ! in case of data subsets restricted to individual frequencies, compute individual misfits
     if(restrict_to_frequencies) then
        do jf = 1,nf
           if(associated(ifreq)) deallocate(ifreq)
           allocate(ifreq(1)); ifreq(1) = ifreq_all(jf)
           if(associated(idx)) deallocate(idx)
!
           idx => getIndxDataSamples(dmspace,ifreq=ifreq)
!
           if(.not.associated(idx)) then
              write(*,*) "MISFIT OF DATA SUBSET FOR jf =",ifreq_all(jf),&
                   " : requested frequency index not contained in dataset"
           else
              misfit = real( sum( (dble(mdata(idx))-dble(sdata(idx)))**2 ) )
              iteration_misfits(j,jf+1) = misfit
              iteration_misfits_normalized(j,jf+1) = real( dble(misfit)/sum(dble(mdata(idx))**2) )
              write(*,*) "MISFIT OF DATA SUBSET FOR jf =",ifreq_all(jf)," : ",misfit,"  NORMALIZED: ",&
                   real( dble(misfit)/sum(dble(mdata(idx))**2) )
           end if
!
        end do ! jf
     end if
     write(*,*) ""
     call dealloc(KLSE_sdata)
!
  end do ! jiter
!
!------------------------------------------------------------------------
!  write out the two arrays iteration_misfits,iteration_misfits_normalized
!
  write(*,*) "writing misfit values now to file '"//trim(outfile_misfit)//"'"
  errmsg = writeAsciiData(trim(outfile_misfit),get(fuh),iteration_misfits)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
  write(*,*) "writing normalized misfit values now to file '"//trim(outfile_misfit_normalized)//"'"
  errmsg = writeAsciiData(trim(outfile_misfit_normalized),get(fuh),iteration_misfits_normalized)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) stop
  call dealloc(errmsg)
!
  write(*,*) ""
!
!------------------------------------------------------------------------
!  clean up
!
   write(*,*) "good bye"
!
   if(associated(ifreq_all)) deallocate(ifreq_all)
   if(allocated(jiter)) deallocate(jiter)
   if(allocated(iteration_misfits)) deallocate(iteration_misfits)
   if(allocated(iteration_misfits_normalized)) deallocate(iteration_misfits_normalized)
   if(associated(idx)) deallocate(idx)
   if(associated(ifreq)) deallocate(ifreq)
!
   call dealloc(KLSE_mdata); call dealloc(KLSE_sdata)
!
   call dealloc(dmspace)
!
   call dealloc(invbasics)
   call dealloc(fuh)
   call dealloc(ap)
!
end program computeIterationMisfits
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'    computeIterationMisfits [-h] [-jf "nf jf1..jfn"] dmspace_file parfile'
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
