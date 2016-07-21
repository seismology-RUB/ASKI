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
program kernel2vtk
   use inversionBasics
   use iterationStepBasics
   use modelParametrization
   use componentTransformation
   use spectralWaveformKernel
   use kernelDisplacement
   use kernelGreenTensor
   use kernelReferenceModel
   use invgridVtkFile
   use wpVtkFile
   use argumentParser
   use string
   use fileUnitHandler
   use errorMessage

   implicit none

   type (argument_parser) :: ap
   character(len=max_length_string) :: main_parfile,str
   character(len=max_length_string), dimension(:), pointer :: str_vec


   type (file_unit_handler) :: fuh
   type (error_message) :: errmsg
   character(len=10) :: myname = 'kernel2vtk'

   type (inversion_basics) :: invbasics
   type (iteration_step_basics) :: iterbasics

   type (kernel_displacement) :: kd
   type (kernel_green_tensor) :: kgt
   type (kernel_reference_model) :: pskrm

   type (spectral_waveform_kernel) :: kernel

   type (invgrid_vtk_file), dimension(:,:), allocatable :: ig_vtk
   type (wp_vtk_file), dimension(:,:), allocatable :: wp_vtk

   character(len=character_length_evid) :: evid
   character(len=character_length_staname) :: staname

   integer :: lu
   integer :: nparam,jparam,nfreq,jfreq,ncomp,jcomp
   character(len=character_length_param), dimension(:), allocatable :: param
   character(len=character_length_pmtrz) :: parametrization
   integer, dimension(:), pointer :: all_ifreq,ifreq
   character(len=character_length_component), dimension(:), allocatable :: comp

   character(len=400) :: kernel_file,kd_file,kgt_file,vtk_file_base,vtk_file_title,vtk_file_data_name,errstr
   complex, dimension(:,:), pointer :: k
   integer, dimension(:), pointer :: cells_filled,wp_inside

   logical :: use_all_ifreq,use_selected_ifreq
   logical :: kernels_on_wp,path_specific,terminate_program,kernel_file_exists,kernel_file_is_on_invgrid

   real :: df

   nullify(str_vec,all_ifreq,ifreq,k,cells_filled,wp_inside)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"Writes spectral sensitivity kernels as vtk files for specific paths, parameters, "//&
       "components and frequencies. Can either handle pre-integrated kernel values on inversion grid cells, or "//&
       "original kernel values on wavefield points (flag -wp).")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-evid',.true.,"(mandatory) defines the event id of the one path (must belong to an event in"//&
       " main event list)",'sval','')
  call addOption(ap,'-staname',.true.,"(mandatory) defines the station name of the one path (must belong to a "//&
        "station in main station list",'sval','')
  call addOption(ap,'-comp',.true.,"(mandatory) vector of receiver components for which the kernel should be "//&
       "computed; valid components are: "//all_valid_components,'svec','')
  call addOption(ap,'-param',.true.,"(mandatory) vector of parameter names for which the kernel should be "//&
       "computed",'svec','')
  call addOption(ap,'-ifreq',.true.,"explicit vector of frequency indices at which the wavefield output should "//&
       "be extracted. Exactly one of options -ifreq , -all_ifreq must be set",'ivec','')
  call addOption(ap,'-all_ifreq',.false.,"if set, all frequency indices are used. Exactly one of options "//&
       "-ifreq , -all_ifreq must be set")
  call addOption(ap,"-wp",.false.,"if set, the original kernels on the WAVEFIELD POINTS are produced. If not "//&
       "set, pre-integrated kernel values on inversion grid cells are produced")
!
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if
!
  use_all_ifreq = ap.optset.'-all_ifreq'
  use_selected_ifreq = ap.optset.'-ifreq'
  if(use_all_ifreq .eqv. use_selected_ifreq) then
     write(*,*) "ERROR: exactly ONE of the options -ifreq and -all_ifreq must be set!"
     call usage(ap)
     goto 1
  end if
  if(.not.( (ap.optset.'-evid').and. (ap.optset.'-staname').and.(ap.optset.'-comp').and.(ap.optset.'-param'))) then
     write(*,*) "ERROR: all of the options -evid, -staname, -comp, -param must be set!"
  end if
!
  ! handle -comp
  if(ap.optset.'-comp') then
     str_vec => ap.svec.'-comp'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     if(.not.associated(str_vec)) then
        write(*,*) "ERROR: for some reason, there is no list of station components returned by argument parser, "//&
             "even though there was no error parsing argument -comp. This is strange..."
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
     ncomp = size(str_vec)
     allocate(comp(ncomp))
     do jcomp = 1,ncomp
        comp(jcomp) = str_vec(jcomp)
     end do
     deallocate(str_vec)
     if(.not.allValidComponents(comp,i_invalid=jcomp)) then
        write(*,*) "ERROR: ",jcomp,"'th component '"//trim(comp(jcomp))//"' not valid. Valid components are '"//&
             all_valid_components//"'"
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
  end if
!
   ! handle -wp
   kernels_on_wp = ap.optset.'-wp'
!
   main_parfile = ap.sval.'main_parfile'
!
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   ! creat file unit handler  
   call createFileUnitHandler(fuh,150)
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,trim(main_parfile),get(fuh),errmsg)
   call undo(fuh)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   ! setup iteration step basics
   call new(errmsg,myname)
   call init(iterbasics,invbasics,fuh,errmsg)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
!------------------------------------------------------------------------
!  processing of the rest of the arguments given by argument parser
!
   ! check if event ID is valid
   str = ap.sval.'-evid'
   evid = str
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   errmsg = searchEventidSeismicEventList(.evlist.invbasics,evid)
   if(.level. errmsg/=0) then
      write(*,*) "ERROR: event ID '"//trim(evid)//"' (input string of option -evid) is not contained in event list"
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
   call dealloc(errmsg)
!
   ! check if station name is valid
   str = ap.sval.'-staname'
   staname = str
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
   errmsg = searchStationNameSeismicNetwork(.statlist.invbasics,staname)
   if(.level. errmsg/=0) then
      write(*,*) "ERROR: station name '"//trim(staname)//"' (input string of option -stname) is not contained in "//&
           "station list"
      write(*,*) ""
      call usage(ap)
      goto 1
   end if
   call dealloc(errmsg)
!
   ! handle -param
   if(ap.optset.'-param') then
      str_vec => ap.svec.'-param'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
      if(.not.associated(str_vec)) then
         write(*,*) "ERROR: for some reason, there is no list of model parameters returned by argument parser, "//&
              "even though there was no error parsing argument -param. This is strange..."
         write(*,*) ""
         call usage(ap)
         goto 1
      end if
      nparam = size(str_vec)
      allocate(param(nparam))
      do jparam = 1,nparam
         param(jparam) = str_vec(jparam)
         if(.not.validParamModelParametrization(parametrization,param(jparam))) then
            write(*,*) jparam,"'th given model parameter '",trim(param(jparam)),"' not valid. Valid "//&
                 "parameterizations(parameters are ",all_valid_pmtrz_param
            write(*,*) ""
            call usage(ap)
            goto 1
         end if
      end do ! jparam
      deallocate(str_vec)
   end if
!
   all_ifreq => .ifreq.iterbasics
   if(use_all_ifreq) then
      nfreq = size(all_ifreq)
      allocate(ifreq(nfreq))
      ifreq = all_ifreq
   end if ! use_all_ifreq
   if(use_selected_ifreq) then
      ifreq => ap.ivec."-ifreq"
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
      if(.not.associated(ifreq)) then
         write(*,*) "ERROR: for some reason, there is no vector of frequency indices returned by argument "//&
              "parser, even though there was no error parsing argument -ifreq. This is strange..."
         write(*,*) ""
         call usage(ap)
         goto 1
      end if
      nfreq = size(ifreq)
      do jfreq=1,nfreq
         if(.not.any(all_ifreq==ifreq(jfreq))) then
            write(*,*) "ERROR: ",jfreq,"'th frequency index ",ifreq(jfreq),&
                " is not an interation step specific frequency as defined by iteration step specific parameter file"
            print *, ""
            call usage(ap)
            goto 1
         end if
      end do
   end if ! use_selected_ifreq
!
   call document(ap)
   write(*,*) ""
!
!------------------------------------------------------------------------
!  program starting here
!
   terminate_program = .false.

   if(kernels_on_wp) then
      allocate(wp_vtk(nparam,ncomp))
      call kernel_on_wp_2_vtk()
      if(terminate_program) goto 1
   else ! kernels_on_wp
      allocate(ig_vtk(nparam,ncomp))
      call kernel_on_invgrid_2_vtk()
      if(terminate_program) goto 1
   end if ! kernels_on_wp
! 
!------------------------------------------------------------------------
!  clean up
!
   print *, "kernel2vtk: good bye"

1  if(allocated(param)) deallocate(param)
   if(allocated(comp)) deallocate(comp)
   if(associated(ifreq)) deallocate(ifreq)
   if(associated(str_vec)) deallocate(str_vec)
   if(allocated(wp_vtk)) then
      do jcomp=1,ncomp
         do jparam=1,nparam
            call dealloc(wp_vtk(jparam,jcomp))
         end do ! jcomp
      end do ! jparam
      deallocate(wp_vtk)
   end if ! allocated(wp_vtk)
   if(allocated(ig_vtk)) then
      do jcomp=1,ncomp
         do jparam=1,nparam
            call dealloc(ig_vtk(jparam,jcomp))
         end do ! jcomp
      end do ! jparam
      deallocate(ig_vtk)
   end if ! allocated(ig_vtk)
   call dealloc(invbasics); call dealloc(iterbasics)
   call dealloc(fuh)
   call dealloc(ap)

   stop
!
!###########################################################################################
!###########################################################################################
!
contains

  subroutine kernel_on_wp_2_vtk()

    df = (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'
    print *,"kernel2vtk: converting plain '"//trim(parametrization)//"'-kernel values on wavefield points to vtk for"
    print *,"  parameters ","'"//param//"',"
    print *,"  frequency indices (corresponding to df = ",df,")  ",ifreq
    print *,"  components ","'"//comp//"',"
    print *,""
    
    ! First check, if the spectral waveform kernel file exists. 
    ! If so, open and see if it contains kernel on wavefield points
    ! Otherwise, recompute kernel on wavefield points

    kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
         'spectral_kernel_ON-WP_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

    inquire(file=kernel_file,exist=kernel_file_exists)

    if(kernel_file_exists) then
       ! check if existing kernel file contains kernel values on wavefield points
       ! if so read in the kernel file, if not recompute
       call new(errmsg,myname)
       call readMetaInfoSpectralWaveformKernel(kernel_file,get(fuh),errmsg,on_invgrid=kernel_file_is_on_invgrid)
       call undo(fuh)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; return; endif
       call dealloc(errmsg)
       if(kernel_file_is_on_invgrid) then
          print *, "kernel file '",trim(kernel_file),"' exists, but contains pre-integrated kernel values on inversion ",&
               "grid cells. Recomputing plain kernel values on wavefield points before writing vtk files."
          print *, ""
          call recompute_kernel_on_wp_2_vtk()
       else ! kernel_file_is_on_invgrid
          print *, "kernel file '",trim(kernel_file),"' exists, and contains plain kernel values on wavefield points. ",&
               "Reading kernel values from file before writing vtk files."
          print *, ""
          call read_kernel_on_wp_2_vtk()
       end if ! kernel_file_is_on_invgrid
    else ! kernel_file_exists
       print *, "kernel file '",trim(kernel_file),"' does not exists. ",&
            "Recomputing plain kernel values on wavefield points before writing vtk files."
       print *, ""
       call recompute_kernel_on_wp_2_vtk()
    end if ! kernel_file_exists

    if(terminate_program) return ! this command is somewhat redudant for now, but keep it in case that code will be added later on
  end subroutine kernel_on_wp_2_vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine read_kernel_on_wp_2_vtk()

    kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
         'spectral_kernel_ON-WP_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

    call new(errmsg,myname)
    call initiateSpectralWaveformKernel(kernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION', &
         param,.ntot.(.wp.iterbasics),comp,errmsg,kernel_on_wp=.true.)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    call new(errmsg,myname)
    call initialReadSpectralWaveformKernel(kernel,kernel_file,get(fuh),errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    df = .df.kernel
    if( df<0. .or. abs(((.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP') -df)/df > 1.e-4 ) then
       write(*,*) "frequency step of kernel ",df," differs from frequency step of measured data ",&
            (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'," by more than 0.01 percent"
       terminate_program = .true.
       goto 1
    end if

    wp_inside => getWpInside(.intw.iterbasics)
    if(.not.associated(wp_inside)) then
       write(*,*) "there are wavefield points which are located inside the inversion grid"
       terminate_program = .true.
       goto 1
    end if

    ! loop on frequencies
    do jfreq = 1,nfreq
       call new(errmsg,myname)
       call readSpectralWaveformKernel(kernel,ifreq(jfreq),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
       call dealloc(errmsg)

       ! loop on all components
       do jcomp=1,ncomp
          k => getValuesByCompSpectralWaveformKernel(kernel,comp(jcomp))
          if(.not.associated(k)) then
             write(*,*) "no '"//trim(param(jparam))//"' sensitivity values contained in kernel"
             terminate_program = .true.
             goto 1
          end if
          ! at this point, k should have size(.ntot.(.wp.iterbasics),nparam)
          ! otherwise, module spectralWaveformKernel is corrupt

          ! handle all parameters for which this kernel was initiated (i.e. all for which vtk files are requested)
          do jparam=1,nparam
             if(jfreq==1) then
                write(vtk_file_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
                ! initiate vtk file
                write(vtk_file_title,*) trim(comp(jcomp)),"-component of spectral ",trim(param(jparam)),&
                     '-'//trim(parametrization)//' Kernel at frequency ',ifreq(jfreq)*df,' Hz on wavefield points'
                call new(errmsg,myname)
                call init(wp_vtk(jparam,jcomp),.wp.iterbasics,.invgrid.iterbasics,trim(vtk_file_base),&
                     trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtk_file_title),&
                     wp_indx_req=wp_inside)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
                call dealloc(errmsg)
                print *,"kernel2vtk: creating vtk files with basename '"//trim(vtk_file_base)//"'"
             end if ! jfreq==1

             ! write kernel values to vtk file
             write(vtk_file_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'

             call new(errmsg,myname)
             call writeData(wp_vtk(jparam,jcomp),get(fuh),k(wp_inside,jparam),&
                  errmsg,data_name=trim(vtk_file_data_name),file_index=ifreq(jfreq))
             call undo(fuh)
             if (.level.errmsg /= 0) call print(errmsg)
             if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
             call dealloc(errmsg)

             print *,"kernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                  "-"//trim(parametrization)//" Kernel, frequency index ",ifreq(jfreq),", frequency ",ifreq(jfreq)*df," Hz"
          end do ! jparam

       end do ! jcomp
    end do ! jfreq
    print *, ""

    ! CLEAN UP LOCALLY ALLOCATED STUFF
1   call finalReadSpectralWaveformKernel(kernel,lu)
    call add(fuh,lu)
    call dealloc(kernel)
    if(associated(wp_inside)) deallocate(wp_inside)
  end subroutine read_kernel_on_wp_2_vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine recompute_kernel_on_wp_2_vtk()
    ! in case we need to recompute the kernel values, we need to know if we are in path specific mode or not
    ! (otherwise this information is not needed)
    path_specific = lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')

    if(path_specific) then
       kd_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//&
            'kernel_displ_'//trim(evid)//'_'//trim(staname)
       kgt_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_GREEN_TENSORS')//&
            'kernel_gt_'//trim(staname)//'_'//trim(staname)
       
       call new(errmsg,myname)
       call createKernelReferenceModel(pskrm,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(.iterpath.invbasics)//trim(sval(.inpar.iterbasics,'PATH_KERNEL_REFERENCE_MODELS'))//&
            'krm_'//trim(evid)//'_'//trim(staname),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
       call dealloc(errmsg)

    else ! path_specific
         
       kd_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//&
            'kernel_displ_'//trim(evid)
       kgt_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_GREEN_TENSORS')//&
            'kernel_gt_'//trim(staname)
    end if ! path_specific

    call new(errmsg,myname)
    call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh,kd_file,errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    call new(errmsg,myname)
    call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh,kgt_file,comp,&
         .comptrans.invbasics,staname,errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    df = .df.kd
    if( df<0. .or. abs(((.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP') -df)/df > 1.e-4 ) then
       write(*,*) "frequency step of kernel displacement ",df," differs from frequency step of measured data ",&
            (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'," by more than 0.01 percent"
       terminate_program = .true.
       goto 1
    end if
    df = .df.kgt
    if( df<0. .or. abs(((.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP') -df)/df > 1.e-4 ) then
       write(*,*) "frequency step of kernel green tensor ",df," differs from frequency step of measured data ",&
            (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'," by more than 0.01 percent"
       terminate_program = .true.
       goto 1
    end if

    wp_inside => getWpInside(.intw.iterbasics)
    if(.not.associated(wp_inside)) then
       write(*,*) "there are wavefield points which are located inside the inversion grid"
       terminate_program = .true.
       goto 1
    end if

    call new(errmsg,myname)
    call initiateSpectralWaveformKernel(kernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION', &
         param,.ntot.(.wp.iterbasics),comp,errmsg,kernel_on_wp=.true.)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    ! need base kernel file for vtk output
    kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
         'spectral_kernel_ON-WP_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

    ! loop on frequencies
    do jfreq = 1,nfreq
       call new(errmsg,myname) ! use one error message for all calls at this frequency, especially when omitting to write the frequencies to screen!
       write(errstr,*) "computing spectral waveform kernel for freqeuency index ",ifreq(jfreq)
       call add(errmsg,0,errstr,myname)

       call readFrequencyKernelDisplacement(kd,ifreq(jfreq),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif

       call readFrequencyKernelGreenTensor(kgt,ifreq(jfreq),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif

       if (path_specific) then
          call computeSpectralWaveformKernel(kernel,kd,kgt,pskrm,.ufmdata.invbasics,errmsg)
       else
          call computeSpectralWaveformKernel(kernel,kd,kgt,.krm.iterbasics,.ufmdata.invbasics,errmsg)
       endif
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif

       ! loop on all components
       do jcomp=1,ncomp
          k => getValuesByCompSpectralWaveformKernel(kernel,comp(jcomp))
          if(.not.associated(k)) then
             write(*,*) "no sensitivity values for component '"//trim(comp(jcomp))//"' contained in kernel"
             terminate_program = .true.
             goto 1
          end if
          ! at this point, k should have size(.ntot.(.wp.iterbasics),nparam)
          ! otherwise, module spectralWaveformKernel is corrupt

          ! handle all parameters for which this kernel was initiated (i.e. all for which vtk files are requested)
          do jparam=1,nparam
             if(jfreq==1) then
                write(vtk_file_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
                ! initiate vtk file
                write(vtk_file_title,*) trim(comp(jcomp)),"-component of spectral ",trim(param(jparam)),&
                     '-'//trim(parametrization)//' Kernel at frequency ',ifreq(jfreq)*df,' Hz on wavefield points'
                call new(errmsg,myname)
                call init(wp_vtk(jparam,jcomp),.wp.iterbasics,.invgrid.iterbasics,trim(vtk_file_base),&
                     trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtk_file_title),&
                     wp_indx_req=wp_inside)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
                call dealloc(errmsg)
                print *,"kernel2vtk: creating vtk files with basename '"//trim(vtk_file_base)//"'"
             end if ! jfreq==1

             ! write kernel values to vtk file
             write(vtk_file_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'

             call new(errmsg,myname)
             call writeData(wp_vtk(jparam,jcomp),get(fuh),k(wp_inside,jparam),&
                  errmsg,data_name=trim(vtk_file_data_name),file_index=ifreq(jfreq))
             call undo(fuh)
             if (.level.errmsg /= 0) call print(errmsg)
             if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
             call dealloc(errmsg)

             print *,"kernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                  "-"//trim(parametrization)//" Kernel, frequency index ",ifreq(jfreq),", frequency ",ifreq(jfreq)*df," Hz"
          end do ! jparam

       end do ! jcomp
    end do ! jfreq
    print *, ""

    ! CLEAN UP LOCALLY ALLOCATED STUFF
1   if(path_specific) call dealloc(pskrm)
    call dealloc(kernel)
    call dealloc(kd,fuh)
    call dealloc(kgt,fuh)
    if(associated(wp_inside)) deallocate(wp_inside)
  end subroutine recompute_kernel_on_wp_2_vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine kernel_on_invgrid_2_vtk()

    kernel_file = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
         'spectral_kernel_'//trim(parametrization)//'_'//trim(evid)//'_'//trim(staname)

    call new(errmsg,myname)
    call initiateSpectralWaveformKernel(kernel,parametrization,param,.ncell.(.invgrid.iterbasics),comp,errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    call new(errmsg,myname)
    call initialReadSpectralWaveformKernel(kernel,kernel_file,get(fuh),errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
    call dealloc(errmsg)

    df = .df.kernel
    if( df<0. .or. abs(((.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP') -df)/df > 1.e-4 ) then
       write(*,*) "frequency step of kernel ",df," differs from frequency step of measured data ",&
            (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'," by more than 0.01 percent"
       terminate_program = .true.
       goto 1
    end if

    cells_filled => getFilledCells(.intw.iterbasics)
    if(.not.associated(cells_filled)) then
       write(*,*) "there are only empty inversion grid cells (or no cells at all)"
       terminate_program = .true.
       goto 1
    end if

    print *,"kernel2vtk: converting pre-integrated '"//trim(parametrization)//"'-kernel on invgrid to vtk for"
    print *,"  parameters ","'"//param//"',"
    print *,"  frequency indices (corresponding to df = ",df,")  ",ifreq
    print *,"  components ","'"//comp//"',"
    print *,""

    ! loop on frequencies
    do jfreq = 1,nfreq
       call new(errmsg,myname)
       call readSpectralWaveformKernel(kernel,ifreq(jfreq),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
       call dealloc(errmsg)

       ! loop on all components
       do jcomp=1,ncomp
          k => getValuesByCompSpectralWaveformKernel(kernel,comp(jcomp))
          if(.not.associated(k)) then
             write(*,*) "no '"//trim(param(jparam))//"' sensitivity values contained in kernel"
             terminate_program = .true.
             goto 1
          end if
          ! at this point, k should have size(.ncell.(.invgrid.iterbasics),nparam)
          ! otherwise, module spectralWaveformKernel is corrupt

          ! handle all parameters for which this kernel was initiated (i.e. all for which vtk files are requested)
          do jparam=1,nparam
             if(jfreq==1) then
                write(vtk_file_base,"(a,'_',a,'_',a)") trim(kernel_file),trim(param(jparam)),trim(comp(jcomp))
                ! initiate vtk file
                write(vtk_file_title,*) trim(comp(jcomp)),"-component of spectral ",trim(param(jparam)),&
                     '-'//trim(parametrization)//' Kernel at frequency ',ifreq(jfreq)*df,' Hz on inversion grid'
                call new(errmsg,myname)
                call init(ig_vtk(jparam,jcomp),.invgrid.iterbasics,trim(vtk_file_base),&
                     trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtk_file_title),&
                     cell_indx_req=cells_filled)
                if (.level.errmsg /= 0) call print(errmsg)
                if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
                call dealloc(errmsg)
                print *,"kernel2vtk: creating vtk files with basename '"//trim(vtk_file_base)//"'"
             end if ! jfreq==1

             ! write kernel values to vtk file
             write(vtk_file_data_name,*) trim(comp(jcomp)),'_',trim(param(jparam)),'-kernel'

             call new(errmsg,myname)
             call writeData(ig_vtk(jparam,jcomp),get(fuh),k(cells_filled,jparam),&
                  errmsg,data_name=trim(vtk_file_data_name),file_index=ifreq(jfreq))
             call undo(fuh)
             if (.level.errmsg /= 0) call print(errmsg)
             if (.level.errmsg == 2) then; terminate_program = .true.; goto 1; endif
             call dealloc(errmsg)

             print *,"kernel2vtk: "//trim(comp(jcomp))//"-component of "//trim(param(jparam))// &
                  "-"//trim(parametrization)//" Kernel, frequency index ",ifreq(jfreq),", frequency ",ifreq(jfreq)*df," Hz"
          end do ! jparam

       end do ! jcomp
    end do ! jfreq
    print *, ""

    ! CLEAN UP LOCALLY ALLOCATED STUFF
1   call finalReadSpectralWaveformKernel(kernel,lu)
    call add(fuh,lu)
    call dealloc(kernel)
    if(associated(cells_filled)) deallocate(cells_filled)

  end subroutine kernel_on_invgrid_2_vtk

end program kernel2vtk
!
!-----------------------------------------------------------------------------------------------------------------
!
! subroutine printhelp
!   use componentTransformation
!   print '(50(1h-))'
!   print *,'Program kernel2vtk writes spectral sensitivity kernels as vtk files for specific paths, parameters,'
!   print *,'components and frequencies. Can either handle pre-integrated kernel values on inversion grid cells, or'
!   print *,'original kernel values on wavefield points (flag -wp).'
!   print *,'Usage:'
!   print *,'                   kernel2vtk [-h] -evid event_id -stname station_name -comp "ncomp comp" -param "nparam param" '
!   print *,'                                   -ifreq ["nfreq ifreq","all"] [-wp] main_parfile'
!   print *,''
!   print *,"    main_parfile: main parameter file of inversion"
!   print *,''
!   print *,'Options:'
!   print *,''
!   print *,'-h     : print help'
!   print *,''
!   print *,'-evid event_id : defines the event id of the one path (must belong to an event in main event list)'
!   print *,''
!   print *,'-stname station_name   : defines the station name of the one path (must belong to a station in main station '//&
!        'list)'
!   print *,''
!   print *,"-comp : ncomp, number of components follwing; comp, components (one of '"//&
!        trim(all_valid_components)//"')"
!   print *,''
!   print *,"-param : nparam, number of parameters following; param, parameters "//&
!        "(e.g. '2 vp vs', or '3 rho lambda mu')"
!   print *,''
!   print *,"-ifreq     : nfreq, number of frequency indices following; ifreq, frequency indices "//&
!        "(e.g. '4  3 4 6 10','2  1 8')"
!   print *,''
!   print *,"-wp     : if set, then the original kernels on the WAVEFIELD POINTS are produced."
!   print *,"          If the kernel files contain pre-integrated values on inversion grid, kernel values are "//&
!        "RECALCULATED on the wavefield points"
!   print *,"          Otherwise, if the kernel files already contain values on wavefield points, they are simply "//&
!        "written to vtk"
!   print '(50(1h-))'
!   return
! end subroutine printhelp
