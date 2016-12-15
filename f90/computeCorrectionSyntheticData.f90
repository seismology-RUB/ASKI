!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   and Wolfgang Friederich (Ruhr-Universitaet Bochum Germany)
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
!-----------------------------------------------------------------------
!  computeCorrectionSyntheticData
!
!  Calculate corrections to synthetics owing to difference 
!  of path specific m_p and global reference model M,
!  m is 3D model to be determined. For each path p we have:
!  d_p-s_p = K_p*(m-m_p) = K_p*(m-M+M-m_p) => d-s+K_p(m_p-M) = K_p(m-M)
!  Hence, instead of synthetics s_p, we use s_p-K_p(m_p-M) as new synthetics.
!  The correction term computed here is therefore: c_p = -K_p(m_p-M)
!
!  What about frequency dependence of earth model (dispersion)? 
!  In principle, the krm model should be produced from the nm-model on the fly according to frequency
!  when computing the kernel. Then the kernel would be correct. The same applies here for the correction
!  of the synthetics. Read in the node model, convert it to krm and then to kim.
!  But we may argure that the dispersion correction for model differences is second order, since it
!  would be e.g. dmu = dmu_r*(1+eps) = dmur+dmur*eps where eps is the dispersion correction.
!  For this reason, we ignore dispersion here.
!---------------------------------------------------------------
program correctionSyntheticData
    use errorMessage
    use argumentParser
    use fileUnitHandler
    use inversionBasics
    use iterationStepBasics
    use dataModelSpaceInfo
    use kernelInvertedModel
    use kernelReferenceModel
    use spectralWaveformKernel
    use modelParametrization
    use parameterCorrelation
    use realloc
    use inputParameter
    use asciiDataIO
    use string
    implicit none
    type (argument_parser) :: ap
    type (error_message) :: errmsg
    type (file_unit_handler) :: fuh
    type (inversion_basics) :: invbasics
    type (iteration_step_basics) :: iterbasics
    type (data_model_space_info) :: dms
    type (kernel_inverted_model) :: kimglob,kimps
    type (kernel_reference_model) :: krmps
    type (spectral_waveform_kernel) :: kernel
    integer :: lu,j,jf,jfreq,nfreq,ic,nc,iparam,nparam_dms,nparam_kernel
    integer, dimension(:), pointer :: ifreq
    real :: df_mdata
    real, dimension(:), pointer :: model_values
    complex, dimension(:,:), pointer :: kernel_values
    complex, dimension(:,:), allocatable :: corr
    character(len=character_length_component), dimension(:), pointer :: comp_path
    character(len=character_length_pmtrz) :: parametrization
    character(len=character_length_param) :: param_name,param_name2
    character(len=character_length_param), dimension(:), pointer :: param_dms,param_kernel,param_name_p
    integer, dimension(:), pointer :: icell_p,mvindex_p
    real :: c_correlation
    character(len=character_length_evid) :: evid
    character(len=character_length_staname) :: staname
    logical :: next
    character (len=max_length_string) :: main_parfile,dmsifile
    character (len=400) :: psrmpath,pskrmfile,file_kernel,corrfile
    character (len=30) :: myname = 'computeCorrectionSyntheticData'
!----------------------------------------------------------------
    nullify(ifreq,model_values,kernel_values,comp_path,param_dms,param_kernel,param_name_p,icell_p,mvindex_p)
!----------------------------------------------------------------
    call init(ap,myname,'Compute corrections to synthetic data due to change from path to global reference model')
    call addPosarg(ap,'dmsi_file','sval','Data-model-space-info file')
    call addPosarg(ap,'main_parfile','sval','Main parameter file of inversion')
    call parse(ap)
    main_parfile = ap.sval.'main_parfile'
    dmsifile = ap.sval.'dmsi_file'
    if (.level.(.errmsg.ap) == 2) then
       call print(.errmsg.ap)
       call usage(ap)
       goto 1
    end if
    call document(ap); call dealloc(ap)
!
    call new(fuh,20)
!------------------------------------------------------------------------
!  setup basics
!
    call new(errmsg,myname)
    call init(invbasics,trim(main_parfile),get(fuh),errmsg)
    call undo(fuh)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)
!
!  setup iteration step basics
!
    call new(errmsg,myname)
    call init(iterbasics,invbasics,fuh,errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)
!
    ifreq => .ifreq.iterbasics
    nfreq = .nf.iterbasics
!---------------------------------------------------------------
! create data model space object from dmsi-file
!
    call new(errmsg,myname)
    ! call createDataSamplesFromFileDataModelSpaceInfo(dms,.evlist.invbasics,.statlist.invbasics,&
    !      .ifreq.iterbasics,trim(dmsifile),get(fuh),errmsg)
    call createFromFileDataModelSpaceInfo(dms,.evlist.invbasics,.statlist.invbasics,&
         .ifreq.iterbasics,(.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',.ncell.(.invgrid.iterbasics),&
         .intw.iterbasics,trim(dmsifile),get(fuh),errmsg)
    call undo(fuh)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)
!
    param_dms => getAllDifferentParamModelValuesDataModelSpaceInfo(dms)
    if(.not.associated(param_dms)) then
       write(*,*) "ERROR: model space does not contain any model values"
       goto 1
    end if
    nparam_dms = size(param_dms)
!
    ! check if there is any parameter correlation to be done for any parameter.
    ! if so, find out the ADDITIONAL parameters which are not contained in the model space
    ! but for which kernel values will be needed for correlation. create a vector
    ! of parameter names containing param_dms AND the additional ones. This
    ! vector will be used to read kernel values.
    nparam_kernel = nparam_dms
    allocate(param_kernel(nparam_kernel))
    param_kernel = param_dms
    do iparam = 1,nparam_dms
       do while (nextCorrelationParameter(.pcorr.invbasics,param_dms(iparam),param_name2))
          if(.not.any(param_kernel == param_name2)) then
             param_kernel => reallocate(param_kernel,nparam_kernel+1)
             param_kernel(nparam_kernel+1) = param_name2
             nparam_kernel = nparam_kernel + 1
          end if
       end do ! while(nextCorrelationParameter)
    end do ! iparam
!----------------------------------------------------------------------
!  convert global krm to kernelInvertedModel
!
    call new(errmsg,myname)
    call interpolateKernelReferenceToKernelInvertedModel(kimglob,.krm.iterbasics,&
         (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION',&
         .invgrid.iterbasics,.intw.iterbasics,errmsg)
    if (.level.errmsg /= 0) call print(errmsg)
    if (.level.errmsg == 2) goto 1
    call dealloc(errmsg)
!---------------------------------------------------------------------
!  Loop over paths
!
    parametrization = (.inpar.invbasics).sval.'MODEL_PARAMETRIZATION'
    psrmpath = (.iterpath.iterbasics)+((.inpar.iterbasics).sval.'PATH_KERNEL_REFERENCE_MODELS')
    write(*,*) " path index | event ID        | station name    |"
    write(*,*) "------------+-----------------+-----------------+"
!
    do while(nextPathDataModelSpaceInfo(dms,evid,staname,all_comp=comp_path))
!
       write(*,"(i12,a,a15,a,a15,a)") j," | ","'"+evid+"'"," | ","'"+staname+"'"," |"
       pskrmfile = psrmpath+'krm_'+evid+'_'+staname
!
       nc = size(comp_path)
       allocate(corr(nfreq,nc))
!
!  read path kernel reference model
!
       call new(errmsg,myname)
       call createKernelReferenceModel(krmps,(.inpar.invbasics).sval.'FORWARD_METHOD',fuh,pskrmfile,errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
!
!  convert path krm to kernelInvertedModel
!
       call new(errmsg,myname)
       call interpolateKernelReferenceToKernelInvertedModel(kimps,krmps,parametrization,&
            .invgrid.iterbasics,.intw.iterbasics,errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
       call dealloc(krmps)
!
!  take difference of path specific and global kim and overwrite kimps:
!  m_p =  (m_p-m_glob), correction is c_p = -K_p*(m_p-M)
!
       call new(errmsg,myname)
       call summateInstancesKernelInvertedModel(kimps,kimglob,errmsg,+1.0,-1.0)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
!
!  read path specific kernel
!
       file_kernel = (.iterpath.iterbasics)+((.inpar.iterbasics).sval.'PATH_SENSITIVITY_KERNELS')+&
            "spectral_kernel_"+parametrization+"_"+evid+"_"+staname
       call new(errmsg,myname)
       call initiateSpectralWaveformKernel(kernel,parametrization,param_kernel,.ncell.(.invgrid.iterbasics),comp_path,&
            errmsg,kernel_on_wp=.false.)
       if (.level.errmsg == 2) then
          call print(errmsg)
          goto 2
       end if
       ! re-use same error message for initial read (might be more comprehensive in case there is an error)
       call initialReadSpectralWaveformKernel(kernel,file_kernel,get(fuh),errmsg)
       if (.level.errmsg /= 0) call print(errmsg)
       if (.level.errmsg == 2) goto 2
       call dealloc(errmsg)
!
!  check content of kernel file
!
       df_mdata = (.inpar.invbasics).rval.'MEASURED_DATA_FREQUENCY_STEP'
       if ( abs(df_mdata-.df.kernel) > 1.e-4*df_mdata) then
          write(*,*) "ERROR: the frequency step of the frequencies in the kernel file (",.df.kernel,&
               ") differs by more than 0.01 percent from the frequency step of the measured data (",df_mdata,&
               "), which suggests that the kernels were computed w.r.t. a different frequency discretization"
          goto 2
       endif
!
!  Loop over frequency indices listed in ITERATION_STEP_INDEX_OF_FREQ
!
       do jfreq = 1,nfreq
          jf = ifreq(jfreq)
!
          call new(errmsg,myname)
          call readSpectralWaveformKernel(kernel,jf,errmsg)
          if(.level.errmsg /= 0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
          call dealloc(errmsg)
!
!  loop over all parameters which are contained in the model space (and for which the kernel values were read in)
!  kernel_values(icell,comp), model_values(icell)
!
          corr(jfreq,:) = (0.d0,0.d0)
          nullify(param_name_p,icell_p,mvindex_p)
          do iparam = 1,nparam_dms
             ! for this model parameter, get all invgrid cell indices icell_p
             ! it is assumed that for the very same model space the Kernel matrix is solved later on
             ! so, this assures that the corrections are computed for the final kernel matrix
             if(associated(param_name_p)) deallocate(param_name_p)
             if(associated(icell_p)) deallocate(icell_p)
             if(associated(mvindex_p)) deallocate(mvindex_p)
             nullify(param_name_p,icell_p,mvindex_p)
             allocate(param_name_p(1)); param_name_p(1) = param_dms(iparam)
             mvindex_p => getIndxModelValues(dms,param=param_name_p,cell=icell_p)
             ! At this point, mvindex_p should be associated (as param(iparam) is contained in model space).
             ! Otherwise, module dataModelSpaceInfo is corrupt. Don't check here again,
!
             ! get pointer to model and kernel values at ALL inversion grid cells
             model_values => getValuesKernelInvertedModel(kimps,param_dms(iparam))
             kernel_values => getValuesByParamSpectralWaveformKernel(kernel,param_dms(iparam))
             ! check for any inconsistencies
             if(.not.(associated(model_values) .and. associated(kernel_values))) then
                write(*,*) "ERROR at parameter ",param_name," and frequency index ",jf,&
                     ": model values or kernel values not defined"
             end if
             if(size(model_values) /= size(kernel_values,1)) then
                write(*,*) "ERROR at parameter ",param_name," and frequency index ",jf,&
                     ": model and kernel values have different size:"
                write(*,*) "  size(model_values) (ncell) = ",size(model_values)
                write(*,*) "  size(kernel_values,1) (ncell) = ",size(kernel_values,1)
                write(*,*) "  size(kernel_values,2) (ncomp of this path) = ",size(kernel_values,2)
                goto 2
             end if

             ! now compute correction term for this model parameter param_name at the correct set of inversion grid cells icell_p
             do ic = 1,nc
                corr(jfreq,ic) = corr(jfreq,ic)-sum(kernel_values(icell_p,ic)*model_values(icell_p))
             enddo ! ic

             ! if other model parameters shoud be correlated to the current model parameter param_name, loop over those
             if(correlateAnyParameters(.pcorr.invbasics,param_name)) then

                ! technically we could also skip the previous "if(correlateAnyParameters(.pcorr.invbasics,param_name))"
                do while (nextCorrelationParameter(.pcorr.invbasics,param_name,param_name2,c_correlation))
                   ! get the kernel values of parameter param_name2
                   kernel_values => getValuesByParamSpectralWaveformKernel(kernel,param_name2)

                   ! then ADDITIONALLY substract those kernel values (weighted by correlation coefficient c_correlation) from the 
                   ! corrections that are already computed. This way, the correct synthetic corrections term is derived, 
                   ! which arises in the finally used kernel linear system (which is assumed to be set up using the very
                   ! same model space and the very same parameter correlation (if any))
                   do ic = 1,nc
                      corr(jfreq,ic) = corr(jfreq,ic)-c_correlation*sum(kernel_values(icell_p,ic)*model_values(icell_p))
                   enddo ! ic

                end do ! while next correlation parameter
             end if ! correlateAnyParameters

          enddo ! iparam
       enddo ! jfreq
       call finalReadSpectralWaveformKernel(kernel,lu)
       call add(fuh,lu)
       call dealloc(kernel)
!
!  write correction to output files, one file per component
!
       do ic = 1,nc
          corrfile = (.iterpath.iterbasics)+((.inpar.iterbasics).sval.'PATH_SYNTHETIC_DATA')+&
               'corr_'+evid+'_'+staname+'_'+comp_path(ic)
          errmsg = writeComplexVectorAsciiDataIO(corrfile,get(fuh),corr(:,ic))
          call undo(fuh)
          if(.level.errmsg /= 0) call print(errmsg)
          if (.level.errmsg == 2) goto 2
       end do ! ic
!
!  clean up before dealing with next path
!
       call dealloc(kimps)
       deallocate(corr)
    end do ! while(nextPath)
!
!  clean up before terminating the program
!
1   call dealloc(errmsg)
    call dealloc(fuh)    
    call dealloc(invbasics); call dealloc(iterbasics)
    call dealloc(dms)
    call dealloc(kimglob); call dealloc(kimps); call dealloc(krmps)
    call dealloc(kernel)
    if(associated(param_dms)) deallocate(param_dms)
    if(associated(param_kernel)) deallocate(param_kernel)
    if(associated(icell_p)) deallocate(icell_p)
    if(associated(param_name_p)) deallocate(param_name_p)
    if(associated(mvindex_p)) deallocate(mvindex_p)
    if(allocated(corr)) deallocate(corr)
!
    ! terminate program if code comes here
    stop
!
    ! reset the iterator in order to account for those cases where the above loop was exited before it 
    ! finished normally (pass all_comp to the function in order to deallocate it)
2   next = nextPathDataModelSpaceInfo(dms,evid,staname,all_comp=comp_path,reset=.true.)
    goto 1
end program correctionSyntheticData
