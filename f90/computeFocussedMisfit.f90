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
program computeFocussedMisfit
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use kernelLinearSystem
  use asciiDataIO
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,foc_coef_file

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=21) :: myname = 'computeFocussedMisfit'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace

  integer :: nfocvol,ifocvol
  real, dimension(:,:), pointer :: foc_coef

  type (kernel_linear_system) :: KLSE

  real, dimension(:), pointer :: mdata,sdata

  real, dimension(:), allocatable :: fmisfit
  real :: misfit
  integer :: ios

  nullify(foc_coef,mdata,sdata)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,'compute focussed misfit of a dataset applying the focussing coefficients produced by program '//&
       'focusSpectralKernels')
  call addPosarg(ap,'dmspace_file','sval','data model space input file which defines data and model space')
  call addPosarg(ap,'foc_coef_file','sval','text file containing the focussing coefficients, as written by program '//&
       'focusSpectralKernels')
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  dmspace_file = ap.sval.'dmspace_file'
  foc_coef_file = ap.sval.'foc_coef_file'
  main_parfile = ap.sval.'main_parfile'
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
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
!  read in focussing coefficients from file foc_coef_file
!
  !! FOR NOW, ONE COLUMN IS ASSUMED IN COEFFICIENTS FILE (i.e. ONE focussing volume)
  !! IN THE FUTURE: READ NUMBER OF COLUMNS (i.e. focussing volums) FROM COMMAND LINE OPTION
  nfocvol = 1
  errmsg = readAsciiData(foc_coef_file,get(fuh),foc_coef,nfocvol)
  call undo(fuh)
  if (.level.errmsg /= 0) call print(errmsg)
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  if(.not.associated(foc_coef)) then
     write(*,*) "ERROR: array of focussing coefficients read from file '",&
          trim(foc_coef_file),"' is not associated"
     if(.level.errmsg == 0) call print(errmsg)
     goto 1
  end if
  call dealloc(errmsg)
!
  if(size(foc_coef,1)/=.ndata.dmspace) then
     write(*,*) "ERROR: number of focussing coefficients fread from file ",size(foc_coef,1),&
          " does not equal the number of data samples ",.ndata.dmspace," in data space"
     goto 1
  end if
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
  !  path_measured_data,lu,errmsg,ignore_data_weights,apply_mdata_normalization)
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
  ! subroutine readSyntheticDataSerialKernelLinearSystem(this,nfreq_measured_data,ifreq_measured_data,&
  !      nfreq_synthetic_data,ifreq_synthetic_data,path_synthetic_data,lu,errmsg,&
  !      apply_event_filter,path_event_filter,apply_station_filter,path_station_filter,&
  !      ignore_data_weights,apply_sdata_normalization,read_synthetic_corrections)
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
!  compute focussed misfit(s)
!
   write(*,*) ""
   write(*,*) "the focussed misfits are computed in different ways (for testing, in order to see what is sensible to look at)"
   write(*,*) ""
!
   allocate(fmisfit(nfocvol))
   sdata => .sd.KLSE
   mdata => .md.KLSE
   if(.not.associated(sdata)) then
      write(*,*) "ERROR: synthetic data was not read successfully"
      goto 1
   end if
   if(.not.associated(mdata)) then
      write(*,*) "ERROR: measured data was not read successfully"
      goto 1
   end if
!
   misfit = getMisfitKernelLinearSystem(KLSE,iostat=ios)
   if(ios/=0) then
      write(*,*) "ERROR: computing regular misfit raised iostat = ",ios
      goto 1
   end if
   write(*,*) "REGULAR (non-focussed) MISFIT:    "
   write(*,*) misfit
   write(*,*) "normalized misfit, i.e.  misfit/(mdata**2):"
   write(*,*) misfit/sum(mdata**2)
   write(*,*) ""
!
   ! BUILD SIMPLE WEIGHTED SUM OF DATA RESIDUALS
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * (mdata-sdata) )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF DATA RESIDUALS: "
   write(*,*) fmisfit
   write(*,*) "including normalization of residuals by absolute measured data values:"
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * ((mdata-sdata)/abs(mdata)) )
   end do ! ifocvol
   write(*,*) fmisfit
!
   ! BUILD WEIGHTED SUM OF ABSOLUTE DATA RESIDUALS
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * abs(mdata-sdata) )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF ABSOLUTE DATA RESIDUALS: "
   write(*,*) fmisfit
   write(*,*) "including normalization of residuals by absolute measured data values:"
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * abs(mdata-sdata)/abs(mdata) )
   end do ! ifocvol
   write(*,*) fmisfit
!
   ! BUILD WEIGHTED SUM OF SQUARED DATA RESIDUALS
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * (mdata-sdata)**2 )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF SQUARED DATA RESIDUALS: "
   write(*,*) fmisfit
   write(*,*) "including normalization of residuals by absolute measured data values:"
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * ((mdata-sdata)/abs(mdata))**2 )
   end do ! ifocvol
   write(*,*) fmisfit
!
   ! BUILD WEIGHTED SUM OF ABSOLUTE DATA RESIDUALS WITH ABSOLUTE WEIGHTS
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( abs(foc_coef(:,ifocvol)) * abs(mdata-sdata) )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF ABSOLUTE DATA RESIDUALS WITH ABSOLUTE WEIGHTS: "
   write(*,*) fmisfit
   write(*,*) "including normalization of residuals by absolute measured data values:"
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( abs(foc_coef(:,ifocvol)) * abs(mdata-sdata)/abs(mdata) )
   end do ! ifocvol
   write(*,*) fmisfit
!
   ! BUILD WEIGHTED SUM OF SQUARED DATA RESIDUALS WITH SQUARED WEIGHTS
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( (foc_coef(:,ifocvol))**2 * (mdata-sdata)**2 )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF SQUARED DATA RESIDUALS WITH SQUARED WEIGHTS: "
   write(*,*) fmisfit
   write(*,*) "including normalization of residuals by absolute measured data values:"
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( (foc_coef(:,ifocvol))**2 * ((mdata-sdata)/abs(mdata))**2 )
   end do ! ifocvol
   write(*,*) fmisfit
!
   ! BUILD WEIGHTED SUM OF SQUARED DATA RESIDUALS, DIVIDED BY REGULAR (non-weighted) MISFIT
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * (mdata-sdata)**2 )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF SQUARED DATA RESIDUALS, DIVIDED BY REGULAR MISFIT: "
   write(*,*) fmisfit * (1./misfit)
   write(*,*) ""
!
   ! BUILD WEIGHTED SUM OF SQUARED DATA RESIDUALS, DIVIDED BY WEIGHTED SUM OF SQUARED MDATA VALUES
   do ifocvol=1,nfocvol
      fmisfit(ifocvol) = sum( foc_coef(:,ifocvol) * (mdata-sdata)**2 ) / sum( foc_coef(:,ifocvol) * (mdata)**2 )
   end do ! ifocvol
   write(*,*) ""
   write(*,*) "WEIGHTED SUM OF SQUARED DATA RESIDUALS, DIVIDED BY WEIGHTED SUM OF SQUARED MDATA VALUES: "
   write(*,*) fmisfit
   write(*,*) ""
!
!------------------------------------------------------------------------
!  clean up
!
   write(*,*) ""
   write(*,*) "good bye"
!
1  if(associated(foc_coef)) deallocate(foc_coef)
   if(allocated(fmisfit)) deallocate(fmisfit)
!
   call dealloc(KLSE)
   call dealloc(dmspace)
!
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(errmsg)
   call dealloc(fuh)
   call dealloc(ap)
!
end program computeFocussedMisfit
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   print '(50(1h-))'
!   print *,'    computeFocussedMisfit [-h] dmspace_file foc_coef_file main_parfile'
!   print *,''
!   print *,'Arguments:'
!   print *,''
!   print *,"    dmspace_file: data model space input file which defines data and model space"
!   print *,"    foc_coef_file: text file containing the focussing coefficients, as written by program focusSpectralKernels"
!   print *,"    main_parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print '(50(1h-))'
!   return
! end subroutine printhelp

