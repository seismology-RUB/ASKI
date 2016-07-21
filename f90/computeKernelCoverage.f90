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
program computeKernelCoverage
  use inversionBasics
  use iterationStepBasics
  use dataModelSpaceInfo
  use vectorPointer
  use kernelLinearSystem
  use invgridVtkFile
  use modelParametrization
  use argumentParser
  use string
  use fileUnitHandler
  use errorMessage

  implicit none

  type (argument_parser) :: ap
  character(len=max_length_string) :: main_parfile,dmspace_file,vtk_outfile_base

  integer, dimension(:), pointer :: jf1_in,jf2_in,jf_window,jf_all
  integer :: nwin,iwin,nf_window
  logical, dimension(:), allocatable :: map_jf_window

  type (file_unit_handler) :: fuh
  type (error_message) :: errmsg
  character(len=21) :: myname = 'computeKernelCoverage'

  type (inversion_basics) :: invbasics
  type (iteration_step_basics) :: iterbasics

  type (data_model_space_info) :: dmspace
  character(len=character_length_pmtrz) :: parametrization
  character(len=character_length_param), dimension(:), pointer :: param,pparam
  integer :: nparam,iparam
  integer, dimension(:), pointer :: idx_mval,pcell,idx_dsamp
  type (integer_vector_pointer), dimension(:,:), allocatable :: idx_mval_cell

  type (kernel_linear_system) :: KLSE
  integer :: lu1,lu2,ndata_K,nmval_K,icol
  real, dimension(:,:), pointer :: K
  double precision, dimension(:), allocatable :: colsum_K,max_colsum_K_param

  type (invgrid_vtk_file), dimension(:), allocatable :: ig_vtk
  character (len=200) :: vtk_file_title,vtk_file_data_name,vtk_filename_extension
  logical :: overwrite_output

  nullify(jf1_in,jf2_in,jf_window,jf_all,param,pparam,idx_mval,pcell,idx_dsamp,K)

!------------------------------------------------------------------------
!  preliminary processing
!
  call init(ap,myname,"summate the absolute values of the column vectors of a given kernel matrix for nwin frequency windows")
  call addPosarg(ap,'dmspace_file','sval',"Data-model-space-info file to set up the kernel linear system")
  call addPosarg(ap,'vtk_outfile_base','sval',"base name for output files (vtk output on inversion grid only)")
  call addPosarg(ap,'main_parfile','sval',"'Main parameter file of inversion")
  call addOption(ap,'-jf1',.true.,"(mandatory) vector of starting indices of nwin frequency windows "//&
       "(must have the same length, nwin, as -jf2, jf1 <= jf2 for all windows)",'ivec','')
  call addOption(ap,'-jf2',.true.,"(mandatory) vector of end indices of nwin frequency windows "//&
       "(must have the same length, nwin, as -jf1, jf1 <= jf2 for all windows)",'ivec','')
  call addOption(ap,'-overwr',.false.,"(optional) indicate to overwrite the file output (NO overwrite by default )")
  call parse(ap)
  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if

  ! get values of positional arguments
  dmspace_file = ap.sval.'dmspace_file'
  vtk_outfile_base = ap.sval.'vtk_outfile_base'
  main_parfile = ap.sval.'main_parfile'

  overwrite_output = ap.optset.'-overwr'

  if (.level.(.errmsg.ap) == 2) then
     call print(.errmsg.ap)
     call usage(ap)
     goto 1
  end if

   ! jf1
   if(.not.(ap.optset.'-jf1')) then
      write(*,*) "ERROR: please indicate starting indices of nwin>0 frequency windows by option '-jf1'"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      jf1_in => ap.ivec.'-jf1'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
   end if
   ! jf2
   if(.not.(ap.optset.'-jf2')) then
      write(*,*) "ERROR: please indicate end indices of nwin>0 frequency windows by option '-jf2'"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      jf2_in => ap.ivec.'-jf2'
      if (.level.(.errmsg.ap) == 2) then
         call print(.errmsg.ap)
         call usage(ap)
         goto 1
      end if
   end if
   ! check if jf2 >= jf1 for all time windows
   nwin = size(jf1_in)
   if(size(jf2_in)/=nwin .or. nwin <= 0) then
      write(*,*) "ERROR: size of vector -jf1 = ",nwin," and size of vector -jf2 =",size(jf2_in)," must be equal and > 0"
      write(*,*) ""
      call usage(ap)
      goto 1
   else
      do iwin=1,nwin
         if(jf1_in(iwin)>jf2_in(iwin)) then
            write(*,*) "ERROR: for ",iwin,"'th frequency window (out of ",nwin,") jf1,jf2 do not fulfill jf1 <= "//&
                 "jf2 :  jf1,jf2 =",jf1_in(iwin),jf2_in(iwin)
            write(*,*) ""
            call usage(ap)
            goto 1
         end if
      end do ! iwin
   end if ! size(jf2)/=nwin .or. nwin <= 0
!
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
   call createFileUnitHandler(fuh,150)
!
   ! setup inversion basics
   call new(errmsg,myname)
   call init(invbasics,main_parfile,get(fuh),errmsg)
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
   jf_all => .ifreq.iterbasics
   allocate(map_jf_window(size(jf_all)))
!
!------------------------------------------------------------------------
!  setup data model space info and read in Kernel matrix
!
   write(*,*) "creating data model space info from file '"//trim(dmspace_file)//"'"
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
   write(*,*) "there are ",.ndata.dmspace," data samples and ",.nmval.dmspace," model values"
   if(.ndata.dmspace <= 0) then
      write(*,*) "ERROR: there are no data samples in the data space"
      goto 1
   end if
   if(.nmval.dmspace <= 0) then
      write(*,*) "ERROR: there are no model values in the model space"
      goto 1
   end if
   write(*,*) ""
!
   parametrization = .pmtrz.dmspace
!
   param => allParam(dmspace)
   nparam = size(param)
!
   write(*,*) "initiating kernel linear system now"
   !   subroutine initiateSerialKernelLinearSystem(this,dmspace,nrowreg,ncolreg,errmsg)
   call new(errmsg,myname)
   call initiateSerialKernelLinearSystem(KLSE,dmspace,0,0,errmsg)
   !if (.level.errmsg /= 0) call print(errmsg)
   call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   write(*,*) "reading in kernel matrix now"
   !subroutine readMatrixSerialKernelLinearSystem(this,path_event_filter,path_station_filter,df_measured_data,&
   !  nfreq_measured_data,ifreq_measured_data,path_sensitivity_kernels,comptrans,ntot_invgrid,&
   !  lu1,lu2,dmspace,errmsg)
   call new(errmsg,myname)
   lu1 = get(fuh)
   lu2 = get(fuh)
   call readMatrixSerialKernelLinearSystem(KLSE,rval(.inpar.invbasics,'MEASURED_DATA_FREQUENCY_STEP'),&
        ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ'),&
        ivec(.inpar.invbasics,'MEASURED_DATA_INDEX_OF_FREQ',ival(.inpar.invbasics,'MEASURED_DATA_NUMBER_OF_FREQ')),&
        trim(.iterpath.invbasics)//sval(.inpar.iterbasics,'PATH_SENSITIVITY_KERNELS'),&
        .ncell.(.invgrid.iterbasics),.pcorr.invbasics,lu1,lu2,errmsg,&
        apply_event_filter=lval(.inpar.invbasics,'APPLY_EVENT_FILTER'),&
        path_event_filter=sval(.inpar.invbasics,'PATH_EVENT_FILTER'),&
        apply_station_filter=lval(.inpar.invbasics,'APPLY_STATION_FILTER'),&
        path_station_filter=sval(.inpar.invbasics,'PATH_STATION_FILTER'))
   call add(fuh,lu1); call add(fuh,lu2)
   if (.level.errmsg /= 0) call print(errmsg)
   !call print(errmsg)
   if (.level.errmsg == 2) goto 1
   call dealloc(errmsg)
!
   K => .KM.KLSE
   ndata_K = .ndata.KLSE
   nmval_K = .nmval.KLSE
   allocate(colsum_K(nmval_K))
!
!------------------------------------------------------------------------
!  For reasons of normalization, compute the absolute sum of the complete column vectors and take 
!  maximum for each model parameter.
!  Simultaneously, remeber for each model parameter its model value and invgrid cell indices for use below
!
   do icol = 1,nmval_K
      colsum_K(icol) = sum( abs( dble( K(1:ndata_K,icol) ) ) )
   end do ! icol
!
   allocate(max_colsum_K_param(nparam),idx_mval_cell(2,nparam))
   do iparam = 1,nparam
      if(associated(pparam)) deallocate(pparam); nullify(pparam)
      if(associated(pcell)) deallocate(pcell); nullify(pcell)
      if(associated(idx_mval)) deallocate(idx_mval)
      allocate(pparam(1)); pparam(1) = param(iparam)
      idx_mval => getIndxModelValues(dmspace,param=pparam,cell=pcell)
      if(associated(idx_mval)) then
         max_colsum_K_param(iparam) = maxval(colsum_K(idx_mval))
         call associateVectorPointer(idx_mval_cell(1,iparam),idx_mval)
         call associateVectorPointer(idx_mval_cell(2,iparam),pcell)
         nullify(idx_mval,pcell)
         deallocate(pparam)
      else
         write(*,*) "ERROR: there are no model values for ",iparam,"'th parameter '",trim(param(iparam)),&
              "' (in total ",nparam," parameters contained in model space). This error should not occur!!"
         goto 1
      end if
   end do ! iparam
!
!------------------------------------------------------------------------
!  Initiate one vtk file object for each model parameter, since different model parameter might be associated
!  with a different (sub)set of inversion grid cell indices.
!  The vtk files are used for data output at all frequency windows, using a distinct filename extension
!
   allocate(ig_vtk(nparam))
   do iparam = 1,nparam
      pcell => getVectorPointer(idx_mval_cell(2,iparam))
      vtk_file_title = "normalized "//trim(parametrization)//"-"//trim(param(iparam))//&
           "-kernel coverage for specific frequency window"
      call new(errmsg,myname)
      call init(ig_vtk(iparam),.invgrid.iterbasics,trim(vtk_outfile_base),&
           trim((.inpar.invbasics).sval.'DEFAULT_VTK_FILE_FORMAT'),errmsg,vtk_title=trim(vtk_file_title),&
           cell_indx_req=pcell)
      if (.level.errmsg /= 0) call print(errmsg)
      !call print(errmsg)
      if (.level.errmsg == 2) goto 1
      call dealloc(errmsg)
   end do ! iparam
!
!------------------------------------------------------------------------
!  iterate over all frequency windows and compute the absolute sum of the (sub)column vectors of the kernel matrix
!
   do iwin = 1,nwin

      write(*,*) "processing frequency window ",iwin," out of ",nwin,":   frequency indices range from ",&
           jf1_in(iwin)," to ",jf2_in(iwin)

      ! create a vector containing frequency indices of iwin'th window which are valid, i.e. which are defined
      ! in the iter parfile, i.e. contained in jf_all
      map_jf_window = jf_all >= jf1_in(iwin) .and. jf_all <= jf2_in(iwin)
      nf_window = count(map_jf_window)
      if(nf_window > 0) then
         if(associated(jf_window)) deallocate(jf_window)
         allocate(jf_window(nf_window))
         jf_window = pack(jf_all,map_jf_window)
         write(*,*) "     there are ",nf_window," valid frequency indices in this window : ",jf_window
      else
         write(*,*) "     there are no valid frequency indices in this window (i.e. which are used in this iteration step)"
         cycle
      end if

      ! find all data samples which belong to this frequency window
      if(associated(idx_dsamp)) deallocate(idx_dsamp)
      idx_dsamp => getIndxDataSamples(dmspace,ifreq=jf_window)
      if(associated(idx_dsamp)) then
         write(*,*) "     there are ",size(idx_dsamp)," data samples in the data space associated with these frequencies"
         do icol = 1,nmval_K
            colsum_K(icol) = sum( abs( dble( K(idx_dsamp,icol) ) ) )
         end do ! icol

      else
         write(*,*) "     there are no data samples in the data space associated with these frequencies"
         cycle
      end if

      ! loop on all parameters and produce vtk output file for each parameter at this frequency window
      do iparam = 1,nparam
         idx_mval => getVectorPointer(idx_mval_cell(1,iparam))

         ! normalize the column sum values for this parameter
         colsum_K(idx_mval) = colsum_K(idx_mval) / max_colsum_K_param(iparam)

         ! write column sum values for this parameter to vtk file
         vtk_file_data_name = "normalized_coverage"
         write(vtk_filename_extension,"(a,i6.6,'-',i6.6)") "_"//trim(parametrization)//"-"//trim(param(iparam))//"_",&
              int(minval(jf_window)),int(maxval(jf_window))
         call new(errmsg,myname)
         call writeData(ig_vtk(iparam),get(fuh),real(colsum_K(idx_mval)),errmsg,data_name=trim(vtk_file_data_name),&
              overwrite=overwrite_output,fname_extension=trim(vtk_filename_extension))
         call undo(fuh)
         if (.level.errmsg /= 0) call print(errmsg)
         !call print(errmsg)
         if (.level.errmsg == 2) goto 1
         call dealloc(errmsg)
      end do ! iparam

   end do ! iwin
!
!------------------------------------------------------------------------
!  clean up
!
1  if(associated(jf1_in)) deallocate(jf1_in)
   if(associated(jf2_in)) deallocate(jf2_in)
   if(associated(jf_window)) deallocate(jf_window)
   if(associated(param)) deallocate(param)
   if(associated(idx_dsamp)) deallocate(idx_dsamp)
   if(allocated(map_jf_window)) deallocate(map_jf_window)
   if(allocated(idx_mval_cell)) then
      do iparam = 1,size(idx_mval_cell,2)
         call dealloc(idx_mval_cell(1,iparam))
         call dealloc(idx_mval_cell(2,iparam))
      end do ! iparam
      deallocate(idx_mval_cell)
   end if
   if(allocated(colsum_K)) deallocate(colsum_K)
   if(allocated(max_colsum_K_param)) deallocate(max_colsum_K_param)
   if(allocated(ig_vtk)) then
      do iparam = 1,size(ig_vtk)
         call dealloc(ig_vtk(iparam))
      end do
      deallocate(ig_vtk)
   end if
   call dealloc(iterbasics)
   call dealloc(invbasics)
   call dealloc(dmspace)
   call dealloc(KLSE)
   call dealloc(ap)
   call dealloc(errmsg)
   call dealloc(fuh)

end program computeKernelCoverage
