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
program computeKernels
  use inversionBasics
  use iterationStepBasics
  use seismicEvent
  use seismicEventList
  use seismicStation
  use seismicNetwork
  use dataModelSpaceInfo
  use kernelDisplacement
  use kernelGreenTensor
  use kernelReferenceModel
  use spectralWaveformKernel
  use componentTransformation
  use modelParametrization
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: str,main_parfile,dmspace_file,&
       skernel_filebase,kd_filebase,kgt_filebase,psmodelpath
  character(len=max_length_string), dimension(:), pointer :: str_vec

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=400) :: errstr
  character(len=14) :: myname = 'computeKernels'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics
  
  type (data_model_space_info) :: dmspace
  integer :: npath,ipath,ipath1,ipath2
  character(len=character_length_evid) :: evid,evid_opened
  character(len=character_length_staname) :: staname,staname_opened
  integer :: icomp
  character(len=character_length_component), dimension(:), pointer :: comp_path
  character(len=16) :: comp_path_write
  integer :: iparam
  character(len=character_length_param), dimension(:), pointer :: param

  type (kernel_displacement) :: kd
  type (kernel_green_tensor) :: kgt
  type (kernel_reference_model) :: pskrm

  type (spectral_waveform_kernel) :: skernel

  logical :: kernels_on_wp,use_dmspace,compute_one_path,path_specific,next,terminate_program
  integer :: jf,lu
  character(len=character_length_pmtrz) :: parametrization

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  PROGRAM STARTS HERE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  nullify(str_vec,comp_path,param)

  terminate_program = .false.

!------------------------------------------------------------------------
!  preliminary processing
!
  ! process command line
  call init(ap,myname,"Compute a set of kernels. The set can be characterized in two different ways: 'way 1' "//&
       "(one path only, using options -evid , -staname , -comp and -param) and 'way 2' (from data model space "//&
       "file take all model parameters and paths (with respective components), using options -dmspace , "//&
       "-ipath1 , -ipath2)")
  call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
  call addOption(ap,'-evid',.true.,"(way 1) defines the event id of the one path (must belong to an event in"//&
       " main event list), mandatory for 'way 1'",'sval','')
  call addOption(ap,'-staname',.true.,"(way 1) defines the station name of the one path (must belong to a "//&
        "station in main station list, mandatory for 'way 1'",'sval','')
  call addOption(ap,'-comp',.true.,"(way 1) vector of receiver components for which the kernel should be "//&
       "computed, mandatory for 'way 1'; valid components are: "//all_valid_components,'svec','')
  call addOption(ap,'-param',.true.,"(way 1) vector of parameter names for which the kernel should be "//&
       "computed, mandatory for 'way 1'",'svec','')
  call addOption(ap,'-dmspace',.true.,"(way 2) data model space input file to define a set of paths (and "//&
       "respective receiver components) as well as model parameters, mandatory for 'way 2'",'sval','')
  call addOption(ap,'-ipath1',.true.,"(way 2) first index of the path loop, optional "//&
       "for 'way 2'",'ival','1')
   call addOption(ap,'-ipath2',.true.,"(way 2) last index of the path loop, optional for 'way 2'",'ival',&
        'max_number_of_paths')
   call addOption(ap,'-wp',.false.,"if set, then plain kernel values on WAVEFIELD POINTS are produced. Otherwise"//&
        " (if not set), pre-integrated kernels on inversion grid cells are computed")
!
   call parse(ap)
   if (.level.(.errmsg.ap) == 2) then
      call print(.errmsg.ap)
      call usage(ap)
      goto 1
   end if
!
   main_parfile = ap.sval.'main_parfile'
!
   use_dmspace = ap.optset.'-dmspace'
   compute_one_path = (ap.optset.'-evid') .and. (ap.optset.'-staname') .and. (ap.optset.'-comp') .and. &
        (ap.optset.'-param')
!
  if(.not. (use_dmspace.or.compute_one_path)) then
     write(*,*) "ERROR: please use either one path ('way 1') or data model space file ('way 2')"
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  if(use_dmspace .and. compute_one_path) then
     write(*,*) "ERROR: please use either one path or data model space file"
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  if(compute_one_path .and. ((ap.optset.'-ipath1').or.(ap.optset.'-ipath2'))) then
     write(*,*) "ERROR: -ipath1, -ipath2 only to be used together with -dmspace"
     write(*,*) ""
     call usage(ap)
     goto 1
  end if
!
  if(use_dmspace .and. ((ap.optset.'-evid').or. (ap.optset.'-staname').or.(ap.optset.'-comp').or.(ap.optset.'-param'))) then
     write(*,*) "ERROR: -evid, -staname, -comp, -param only to be used when NOT using -dmspce"
     write(*,*) ""
     call usage(ap)
     goto 1
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
     allocate(comp_path(size(str_vec)))
     do icomp = 1,size(str_vec)
        comp_path(icomp) = str_vec(icomp)
     end do
     deallocate(str_vec)
     if(.not.allValidComponents(comp_path,i_invalid=icomp)) then
        write(*,*) "ERROR: ",icomp,"'th component '"//trim(comp_path(icomp))//"' not valid. Valid components are '"//&
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
  if (.level.errmsg /= 0) then; call print(errmsg); endif
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!
  ! setup iteration step basics
  call new(errmsg,myname)
  call init(iterbasics,invbasics,fuh,errmsg)
  if (.level.errmsg /= 0) then; call print(errmsg); endif
  !call print(errmsg)
  if (.level.errmsg == 2) goto 1
  call dealloc(errmsg)
!------------------------------------------------------------------------
!  other stuff to do beforehand
!
  path_specific = lval(.inpar.iterbasics,'USE_PATH_SPECIFIC_MODELS')
  psmodelpath = ''
  if(path_specific) psmodelpath = sval(.inpar.iterbasics,'PATH_KERNEL_REFERENCE_MODELS')
!
  parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
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
     allocate(param(size(str_vec)))
     do iparam = 1,size(str_vec)
        param(iparam) = str_vec(iparam)
        if(.not.validParamModelParametrization(parametrization,param(iparam))) then
           write(*,*) "ERROR: ",iparam,"'th given model parameter '",trim(param(iparam)),"' not valid in current "//&
                "parameterization '"//trim(parametrization)//"'. Valid "//&
                "parameterizations(parameters are ",all_valid_pmtrz_param
           write(*,*) ""
           call usage(ap)
           goto 1
        end if
     end do ! iparam
     deallocate(str_vec)
  end if
!
  if(use_dmspace) then
     dmspace_file = ap.sval.'-dmspace'
     if (.level.(.errmsg.ap) == 2) then
        call print(.errmsg.ap)
        call usage(ap)
        goto 1
     end if
     call new(errmsg,myname)
     call createFromFileDataModelSpaceInfo(dmspace,.evlist.invbasics,.statlist.invbasics,&
          .ifreq.iterbasics,parametrization,.ncell.(.invgrid.iterbasics),.intw.iterbasics,&
          trim(dmspace_file),get(fuh),errmsg)
     call undo(fuh)
     call addTrace(errmsg,myname)
     if (.level.errmsg /= 0) call print(errmsg)
     !call print(errmsg)
     if (.level.errmsg == 2) goto 1
     if(.ndata.dmspace == 0 .or. .nmval.dmspace == 0) then
        write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace,&
             "model values contained in the data_model_space_info object. using -dmspace there must BOTH, ",&
             "data samples and model parameters, be defined. Printing now the error message returned from ",&
             "createFromFileDataModelSpaceInfo"
        call print(errmsg)
        goto 1
     end if
     call dealloc(errmsg)
     npath = .npath.dmspace ! cannot be 0, since .ndata.dmspace /= 0
!
     ! define model parameters
     ! param should be returned associated, since .nmval.dmspace /= 0 (assured above)
     param => getAllDifferentParamModelValuesDataModelSpaceInfo(dmspace)
!
     ! define ipath1
     if(.not.(ap.optset.'-ipath1')) then
        ipath1 = 1
     else
        ipath1 = ap.ival.'-ipath1'
        if (.level.(.errmsg.ap) == 2) then
           call print(.errmsg.ap)
           call usage(ap)
           goto 1
        end if
        if(ipath1<1) then
           write(*,*) "ERROR: first path index = ",ipath1,", given by -ipath1 must be at least 1"
           write(*,*) ""
           call usage(ap)
           goto 1         
        end if
     end if
     ! define ipath2
     if(.not.(ap.optset.'-ipath2')) then
        ipath2 = npath
     else
        ipath2 = ap.ival.'-ipath2'
        if (.level.(.errmsg.ap) == 2) then
           call print(.errmsg.ap)
           call usage(ap)
           goto 1
        end if
        if(ipath2>npath) then
           write(*,*) "ERROR: last path index = ",ipath2,", given by -ipath2 is greater than maximum number of "//&
                "paths (",npath,") contained in data-model-space info defined by file '"//trim(dmspace_file)//"'"
           write(*,*) ""
           call usage(ap)
           goto 1         
        end if
     end if
     ! check if ipath1<=ipath2
     if(ipath1>ipath2) then
        write(*,*) "ERROR: first path index ( = ",ipath1,") is greater than last path index ( = ",ipath2,")"
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
!
  else ! use_dmspace
!
     ! previously checked:
     ! if program comes here: compute_one_path = ap.optset.'-evid' .and. ap.optset.'-staname' .and. ap.optset.'-comp' .and. ap.optset.'-param' == .true.
!
     npath = 1
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
        write(*,*) "ERROR: station name '"//trim(staname)//"' (input string of option -stname) is not contained in station list"
        write(*,*) ""
        call usage(ap)
        goto 1
     end if
     call dealloc(errmsg)
!
  end if ! use_dmspace
!
  call document(ap)
  write(*,*) ""
!
!------------------------------------------------------------------------
!  compute kernels
!
  ! base of all kernel displacement filenames
  kd_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_DISPLACEMENTS')//'kernel_displ_'
  ! base of all kernel green tensor filenames
  kgt_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_KERNEL_GREEN_TENSORS')//'kernel_gt_'
  ! base of all spectral kernel filenames
  if(kernels_on_wp) then
     skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
          'spectral_kernel_ON-WP_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
  else
     skernel_filebase = trim(.iterpath.invbasics)//trim((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')//&
          'spectral_kernel_'//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//'_'
  end if
!
  if(use_dmspace) call computeKernelsDmspace()
  if(terminate_program) goto 1 
!
  if(compute_one_path) call computeKernelsOnePath()
  if(terminate_program) goto 1 
!
!------------------------------------------------------------------------
!  clean up
!
1 if(associated(str_vec)) deallocate(str_vec)
  if(associated(comp_path)) deallocate(comp_path)
  if(associated(param)) deallocate(param)
  call dealloc(skernel)
  if (path_specific) call dealloc(pskrm)
  call dealloc(dmspace)
  call dealloc(invbasics); call dealloc(iterbasics)
  call dealloc(fuh)
  call dealloc(ap)
!
  stop
!
contains
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine computeKernelsDmspace()
    if(path_specific) then
       write(*,*) "will now compute spectral '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//&
            "' kernels of parameters ",param//','," in path-specific mode for ",ipath2-ipath1+1," path(s) "
       write(*,*)"Kernel displacements will be read from files: '",trim(kd_filebase)//"EVENTID_STANAME'"
       write(*,*)"Kernel Green tensors will be read from files: '",trim(kgt_filebase)//"EVENTID_STANAME'"
       write(*,*)"Kernel reference models will be taken seperately for each path from files: '",&
            trim(.iterpath.invbasics)//trim(psmodelpath)//"krm_EVENTID_STANAME'"
    else
       write(*,*) "will now compute spectral '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//&
            "' kernels of parameters ",param//','," for ",ipath2-ipath1+1," path(s)"
       write(*,*)"Kernel displacements will be read from files: '",trim(kd_filebase)//"EVENTID'"
       write(*,*)"Kernel Green tensors will be read from files: '",trim(kgt_filebase)//"STANAME'"
    end if
    if(kernels_on_wp) then
       write(*,*) "will compute plain kernel values on wavefield points"
    else
       write(*,*) "will compute pre-integrated kernel values on inversion grid cells"
    end if
    write(*,*) ""
    write(*,*) " path index | event ID        | station name    | station components         |"
    write(*,*) "------------+-----------------+-----------------+----------------------------|"
!
    ! use variables evid_opened / staname_opened to memorize the event ID / station name for which current kd / kgt objects are initiated
    ! indicate by setting evid_opened / staname_opened to '' that no kd,kgt objects have been initiated yet
    evid_opened = ''; staname_opened = ''
    ! incremet index ipath manually, for output on screen
    ipath = ipath1-1
    ! loop on all paths
    do while(nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path,ipath_start=ipath1,ipath_end=ipath2))
       ipath = ipath + 1
!
       ! if there is a next path, then comp_path must be associated here, otherwise module dataModelSpaceInfo 
       ! is corrupt, so trust the output here and do not check again if comp_path is associated
!
       write(comp_path_write,*) comp_path//","
       write(*,"(i12,a,a15,a,a15,2a)") ipath," | ","'"//trim(evid)//"'"," | ","'"//trim(staname)//"'"," | ",trim(comp_path_write)
       !write(*,*) "# COMPUTE SPECTRAL KERNEL for ",ipath,"'th path: event '",trim(paths(1,ipath)),"', station '",trim(paths(2,ipath)),"'"
!
       ! check if for this path the previous kd or kgt files can be used (i.e. left open). Otherwise, close old ones and open new ones
       ! kd object still valid?
       if(path_specific .or. evid_opened /= evid) then    ! always do initiate if path_specific !
          ! close old file if open
          if(evid_opened/='') call dealloc(kd,fuh)
          ! initiate kernel displacement
          call new(errmsg,myname)
          if (path_specific) then
             call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
                  trim(kd_filebase)//trim(evid)//'_'//trim(staname),errmsg)
             !write(*,*)'Kernel displacement from file: ',trim(kd_filebase)//trim(evid)//'_'//trim(staname)
          else
             call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
                  trim(kd_filebase)//trim(evid),errmsg)
          endif
          if (.level.errmsg /= 0) then; call print(errmsg); endif
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
          evid_opened = evid
       end if ! path_specific .or. evid_opened /= evid
       ! kgt object still valid?
       if(path_specific .or. staname_opened /= staname) then    ! always do initiate if path_specific !
          ! close old file if open
          if(staname_opened /= '') call dealloc(kgt,fuh)
          ! initiate kernel Green tensor
          call new(errmsg,myname)
          if (path_specific) then
             call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
                  trim(kgt_filebase)//trim(evid)//'_'//trim(staname),comp_path,.comptrans.invbasics,staname,errmsg)
             !write(*,*)'Kernel Green tensor from file: ',trim(kgt_filebase)//trim(evid)//'_'//trim(staname)
          else
             call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
                  trim(kgt_filebase)//trim(staname),comp_path,.comptrans.invbasics,staname,errmsg)
          endif
          if (.level.errmsg /= 0) then; call print(errmsg); endif
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
          staname_opened = staname
       end if ! path_specific .or. staname_opened /= staname
!
       ! initiate and inital write spectral waveform kernel
       call new(errmsg,myname)
       if(kernels_on_wp) then
          call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param, &
               .ntot.(.wp.iterbasics),comp_path,errmsg,kernel_on_wp=kernels_on_wp)
       else
          call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param, &
               .ncell.(.invgrid.iterbasics),comp_path,errmsg)
       end if
       if (.level.errmsg /= 0) then; call print(errmsg); endif
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)

       call new(errmsg,myname)
       call initialWriteSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),&
            get(fuh),errmsg,nfreq=.nf.kd)
       if (.level.errmsg /= 0) then; call print(errmsg); endif
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
!
!  if path_specific, read in path_specific kernel reference model
!
       if (path_specific) then
          call new(errmsg,myname)
          call createKernelReferenceModel(pskrm,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
               trim(.iterpath.invbasics)//trim(psmodelpath)//'krm_'//trim(evid)//'_'//trim(staname),errmsg)
          if (.level.errmsg /= 0) then; call print(errmsg); endif
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
          !write(*,*)'Kernel reference model from file: ', &
          !     trim(.iterpath.invbasics)//trim(psmodelpath)//'krm_'//trim(evid)//'_'//trim(staname)
       endif
!
       ! now loop on all frequencies
       do while (nextFrequencyKernelDisplacement(kd,jf))
          call new(errmsg,myname) ! use one error message for all calls at this frequency, especially when omitting to write the frequencies to screen!
          write(errstr,*) "computing spectral waveform kernel for freqeuency index ",jf
          call add(errmsg,0,errstr,myname)
          !write(*,*) "#     frequency index ",jf," , frequency ",(.df.kd)*jf," [Hz]"

          ! read current frequency of kernel displacement
          call readFrequencyKernelDisplacement(kd,jf,errmsg)
          if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
          ! read current frequency of kernel Green tensor
          call readFrequencyKernelGreenTensor(kgt,jf,errmsg)
          if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!                  
          ! compute spectral waveform kernel
          if (path_specific) then
             if(kernels_on_wp) then
                call computeOnWpSpectralWaveformKernel(skernel,kd,kgt,pskrm,.ufmdata.invbasics,errmsg)
             else
                call computeOnInvgridSpectralWaveformKernel(skernel,kd,kgt,pskrm,.ufmdata.invbasics,.intw.iterbasics,errmsg)
             end if
          else ! path_specific
             if(kernels_on_wp) then
                call computeOnWpSpectralWaveformKernel(skernel,kd,kgt,.krm.iterbasics,.ufmdata.invbasics,errmsg)
             else
                call computeOnInvgridSpectralWaveformKernel(skernel,kd,kgt,.krm.iterbasics,.ufmdata.invbasics,&
                     .intw.iterbasics,errmsg)
             end if
          endif ! path_specific
          if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
         ! write spectral waveform kernel
          call writeSpectralWaveformKernel(skernel,errmsg)
          if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
          if (.level.errmsg /= 0) then; call print(errmsg); endif ! print warning, if any
          call dealloc(errmsg)
       end do ! while (nextFrequencyKernelDisplacement(kd,jf))
!
       ! final write and deallocate velocity kernel
       call new(errmsg,myname)
       call finalWriteSpectralWaveformKernel(skernel,errmsg,lu)
       if (.level.errmsg /= 0) then; call print(errmsg); endif
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
       call add(fuh,lu)
       call dealloc(skernel)
!
       if (path_specific) call dealloc(pskrm)
!
       !write(*,*) ""
    end do ! while nextPath
!
    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   if(evid_opened/='') call dealloc(kgt,fuh)
    if(staname_opened/='') call dealloc(kd,fuh)
    return
!   
    ! due to an error, terminate program after return
2   terminate_program = .true.
    ! reset the iterator in order to account for those cases where the above loop was exited before it 
    ! finished normally (pass all_comp to the function in order to deallocate it)
    next = nextPathDataModelSpaceInfo(dmspace,evid,staname,all_comp=comp_path,reset=.true.)
    call dealloc(skernel)
    if (path_specific) call dealloc(pskrm)
    goto 1
  end subroutine computeKernelsDmspace
!
!-----------------------------------------------------------------------------------------------------------------
!
  subroutine computeKernelsOnePath()
    if(path_specific) then
       write(*,*) "will now compute spectral '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//&
            "' kernel of parameters ",param//','," in path-specific mode for one path:"
       write(*,*)"event ID = ",trim(evid)
       write(*,*)"station name = ",trim(staname)," , receiver components = ",comp_path//","
       write(*,*)"Kernel displacements will be read from file: '",trim(kd_filebase)//trim(evid)//"_"//trim(staname)//"'"
       write(*,*)"Kernel Green tensors will be read from file: '",trim(kgt_filebase)//trim(evid)//"_"//trim(staname)//"'"
       write(*,*)"Kernel reference model will be taken from files '",&
            trim(.iterpath.invbasics)//trim(psmodelpath)//"krm_"//trim(evid)//"_"//trim(staname)//"'"
    else
       write(*,*) "will now compute spectral '"//trim((.inpar.invbasics).sval.'MODEL_PARAMETRIZATION')//&
            "' kernels of parameters ",param//','," for one path:"
       write(*,*)"event ID = ",trim(evid)
       write(*,*)"station name = ",trim(staname)," , receiver components = ",comp_path//","
       write(*,*)"Kernel displacements will be read from file: '",trim(kd_filebase)//trim(evid)//"'"
       write(*,*)"Kernel Green tensors will be read from file: '",trim(kgt_filebase)//trim(staname)//"'"
    end if
    if(kernels_on_wp) then
       write(*,*) "will compute plain kernel values on wavefield points"
    else
       write(*,*) "will compute pre-integrated kernel values on inversion grid cells"
    end if
!
    ! initiate kernel displacement
    call new(errmsg,myname)
    if (path_specific) then
       call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(kd_filebase)//trim(evid)//'_'//trim(staname),errmsg)
       !write(*,*)'Kernel displacement from file: ',trim(kd_filebase)//trim(evid)//'_'//trim(staname)
    else
       call initiateKernelDisplacement(kd,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(kd_filebase)//trim(evid),errmsg)
    endif
    if (.level.errmsg /= 0) then; call print(errmsg); endif
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    ! initiate kernel Green tensor
    call new(errmsg,myname)
    if (path_specific) then
       call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(kgt_filebase)//trim(evid)//'_'//trim(staname),comp_path,.comptrans.invbasics,staname,errmsg)
       !write(*,*)'Kernel Green tensor from file: ',trim(kgt_filebase)//trim(evid)//'_'//trim(staname)
    else
       call initiateKernelGreenTensor(kgt,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(kgt_filebase)//trim(staname),comp_path,.comptrans.invbasics,staname,errmsg)
    endif
    if (.level.errmsg /= 0) then; call print(errmsg); endif
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    ! initiate and inital write spectral waveform kernel
    call new(errmsg,myname)
    if(kernels_on_wp) then
       call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param, &
            .ntot.(.wp.iterbasics),comp_path,errmsg,kernel_on_wp=.true.)
    else
       call initiateSpectralWaveformKernel(skernel,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',param, &
            .ncell.(.invgrid.iterbasics),comp_path,errmsg)
    end if
    if (.level.errmsg /= 0) then; call print(errmsg); endif
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    call new(errmsg,myname)
    call initialWriteSpectralWaveformKernel(skernel,trim(skernel_filebase)//trim(evid)//'_'//trim(staname),&
         get(fuh),errmsg,nfreq=.nf.kd)
    if (.level.errmsg /= 0) then; call print(errmsg); endif
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
!
    !  if path_specific, read in path_specific kernel reference model
    if (path_specific) then
       call new(errmsg,myname)
       call createKernelReferenceModel(pskrm,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh, &
            trim(.iterpath.invbasics)//trim(psmodelpath)//'krm_'//trim(evid)//'_'//trim(staname),errmsg)
       if (.level.errmsg /= 0) then; call print(errmsg); endif
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
       !write(*,*)'Kernel reference model from file: ', &
          !     trim(.iterpath.invbasics)//trim(psmodelpath)//'krm_'//trim(evid)//'_'//trim(staname)
    endif
!
    ! now loop on all frequencies
    do while (nextFrequencyKernelDisplacement(kd,jf))
       call new(errmsg,myname) ! use one error message for all calls at this frequency, especially when omitting to write the frequencies to screen!
       write(errstr,*) "computing spectral waveform kernel for freqeuency index ",jf
       call add(errmsg,0,errstr,myname)
       !write(*,*) "#     frequency index ",jf," , frequency ",(.df.kd)*jf," [Hz]"

       ! read current frequency of kernel displacement
       call readFrequencyKernelDisplacement(kd,jf,errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
       ! read current frequency of kernel Green tensor
       call readFrequencyKernelGreenTensor(kgt,jf,errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!                  
       ! compute spectral waveform kernel
       if (path_specific) then
          if(kernels_on_wp) then
             call computeOnWpSpectralWaveformKernel(skernel,kd,kgt,pskrm,.ufmdata.invbasics,errmsg)
          else
             call computeOnInvgridSpectralWaveformKernel(skernel,kd,kgt,pskrm,.ufmdata.invbasics,.intw.iterbasics,errmsg)
          end if
       else ! path_specific
          if(kernels_on_wp) then
             call computeOnWpSpectralWaveformKernel(skernel,kd,kgt,.krm.iterbasics,.ufmdata.invbasics,errmsg)
       else
             call computeOnInvgridSpectralWaveformKernel(skernel,kd,kgt,.krm.iterbasics,.ufmdata.invbasics,&
                  .intw.iterbasics,errmsg)
          end if
       end if ! path_specific
       if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
       ! write spectral waveform kernel
       call writeSpectralWaveformKernel(skernel,errmsg)
       if (.level.errmsg == 2) then; call print(errmsg); goto 2; endif
!
       if (.level.errmsg /= 0) then; call print(errmsg); endif ! print warning, if any
       call dealloc(errmsg)
    end do ! while (nextFrequencyKernelDisplacement(kd,jf))
!
    ! final write and deallocate velocity kernel
    call new(errmsg,myname)
    call finalWriteSpectralWaveformKernel(skernel,errmsg,lu)
    if (.level.errmsg /= 0) then; call print(errmsg); endif
    if (.level.errmsg == 2) goto 2
    call dealloc(errmsg)
    call add(fuh,lu)
!
    ! if subroutine comes here, deallocate stuff that was allocated here and return normally
1   call dealloc(kgt,fuh)
    call dealloc(kd,fuh)
    if (path_specific) call dealloc(pskrm)
    call dealloc(skernel)
    return
!   
    ! due to an error, terminate program after return
2   terminate_program = .true.
    goto 1
  end subroutine computeKernelsOnePath

 end program computeKernels
!
!-----------------------------------------------------------------------------------------------------------------
!
  ! subroutine printhelp
  !   use componentTransformation,only: all_valid_components
  !   print '(50(1h-))'
  !   print *,'    computeKernels [-h] [-wp] [-evid event_id] [-stname station_name] [-comp "n comp1 .. compn"] '//&
  !        '[-param "n param1 .. paramn"][-dmspce dmspace_file] [-ipath1 path1] [-ipath2 path2] main_parfile'
  !   print *,''
  !   print *,'THERE ARE TWO POSSIBLE WAYS TO DEFINE A SET OF SPECTRAL KERNELS THAT ARE COMPUTED:'
  !   print *,'  (way 1) compute kernel for only one path, defined by eventID and station name and component(s) using'
  !   print *,'          options -evid , -stname , -comp and -param'
  !   print *,''
  !   print *,'  (way 2) use flag -dmspce in connection with optional range definition of the path index (flags -ipath1 -ipath2)'
  !   print *,'          in order to define a subset of the paths contained in the given data-model-space description.'
  !   print *,'          If the upper (lower) limit of the path range is not defined (i.e. -ipath1 (-ipath2) not set), the maximum'
  !   print *,'          (minimum) possible value is used.'
  !   print *,'          Then, all kernels for the defined range of paths in the given data-model-space description are computed'
  !   print *,'          at the respective receiver components and the model parameters defined in the data-model-space description'
  !   print *,''
  !   print *,'Arguments:'
  !   print *,''
  !   print *,"    parfile: main parameter file of inversion"
  !   print *,''
  !   print *,'Options:'
  !   print *,''
  !   print *,'-h     : print help'
  !   print *,''
  !   print *,'-wp    : if set, then plain kernel values on WAVEFIELD POINTS are produced. Otherwise (if not set),'
  !   print *,'          pre-integrated kernels on inversion grid cells are computed'
  !   print *,''
  !   print *,'-evid event_id : defines the event id of the one path (must belong to an event in main event list) (way 1)'
  !   print *,''
  !   print *,'-stname station_name   : defines the station name of the one path (must belong to a station in main station '//&
  !        'list) (way 1)'
  !   print *,''
  !   print *,'-comp "n comp1..compn" : defines a list of n receiver components for which the kernels should be computed '
  !   print *,'                         (the string must start with the number n, followed by n component names). '
  !   print *,'                         Valid component names are: '//all_valid_components//' (way 1)'
  !   print *,''
  !   print *,'-param "n param1..paramn" : defines a list of n parameter names (e.g. "2 vp vs") for which the kernels should be '
  !   print *,'                         computed (the string must start with the number n, followed by n valid parameter names). '
  !   print *,'                         Valid parameter names are: '//all_valid_pmtrz_param//' (way 1)'
  !   print *,''
  !   print *,"-dmspce  : data model space input file to define a set of paths (way 2)"
  !   print *,''
  !   print *,'-ipath1 path1  : path1 is the first index of the path loop. By default, path1 = 1 (way 2)'
  !   print *,''
  !   print *,'-ipath2 path2  : path2 is the last index of the path loop. By default, path2 = max_number_of_paths as to '//&
  !        'data-model-space description (way 2)'
  !   print '(50(1h-))'
  !   return
  ! end subroutine printhelp
