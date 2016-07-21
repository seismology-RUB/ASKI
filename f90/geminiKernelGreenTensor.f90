!--------------------------------------------------------------------------
!	Copyright 2015 Wolfgang Friederich
!
!	This file is part of Gemini II.
!
!	Gemini II is free software: you can redistribute it and/or modify
!	it under the terms of the GNU General Public License as published by
!	the Free Software Foundation, either version 2 of the License, or
!	any later version.
!
!	Gemini II is distributed in the hope that it will be useful,
!	but WITHOUT ANY WARRANTY; without even the implied warranty of
!	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!	GNU General Public License for more details.
!
!	You should have received a copy of the GNU General Public License
!	along with Gemini II.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!---------------------------------------------------------------
!  Module dealing with Gemini kernel Green tensor
!---------------------------------------------------------------
 module geminiKernelGreenTensor
    use seismicStation
    use streamAccess
    use errorMessage
    use vectorPointer
    use fileUnitHandler
    implicit none
    interface dealloc; module procedure deallocGeminiKernelGreenTensor; end interface
    interface operator (.id.); module procedure getIdGeminiKernelGreenTensor; end interface
    interface operator (.df.); module procedure getDfGeminiKernelGreenTensor; end interface
    interface operator (.nf.); module procedure getNfGeminiKernelGreenTensor; end interface
    interface operator (.jfcur.); module procedure getJfcurGeminiKernelGreenTensor; end interface
    interface operator (.disp.); module procedure getGeminiKernelGreenTensor; end interface
    interface operator (.strain.); module procedure getStrainsGeminiKernelGreenTensor; end interface
    type gemini_kernel_green_tensor
        character (len=80) :: id                                        ! id for kernel Green tensor
        integer :: nf1,nf2                                              ! frequency index range
        real :: df                                                      ! frequency stepping
        real :: kfcur                                                   ! current frequency index
        integer :: numtasks                                             ! number of parallel processes for gt-calculation
        integer :: ng                                                   ! total number of wavefield points
        integer :: nnod                                                 ! total number of depth nodes
        complex, dimension(:,:,:), pointer :: g => null()               ! g(wp,comp,force) for current frequency at wp
        complex, dimension(:,:,:), pointer :: gstr => null()            ! gstr(wp,1-6,force) for current frequency at wp
        type (seismic_station) :: si                                    ! station information
        type (file_stream_access), dimension(:), pointer :: fda => null()    ! file pointers to all files
        type (group_stream_access), dimension(:), pointer :: root => null()  ! root group
    end type
!
 contains
!-----------------------------------------------------------------------------------
!> \brief Read in Green tensor and Green tensor gradient for one frequency
!
    subroutine initialReadGeminiKernelGreenTensor(this,fuh,basename,errmsg)
    type (gemini_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename
    type (error_message) :: errmsg
    type (file_stream_access) :: fdaini
    type (group_stream_access) :: rootini
    type (data_stream_access), pointer :: dset
    type (group_stream_access), pointer :: group
    type (flexible), dimension(:), pointer :: ft
    integer :: ierr,j,lu,ntot
    character (len=400) :: filename
    character (len=3) :: crank
    character (len=34) :: myname = 'initialReadGeminiKernelGreenTensor'
!
    nullify(dset,group,ft)
!
    call addTrace(errmsg,myname)
    lu = get(fuh)
    filename = trim(basename)//'_'//'000'
    ierr = openFileStreamAccess(fdaini,lu,trim(filename))
    if (ierr /= 0) then
        call add(errmsg,2,trim(filename)//' can not be opened',myname)
        call add(fuh,lu); return
    endif
!
!  read in group tree
!
    call readGroupStreamAccess(rootini,fdaini)
!    print *,'group tree read in from file: ',trim(filename)
!
!  read out header info: nf1,nf2,df,ntot,numtasks
!
    call traversePathStreamAccess(rootini,0,(/ 1,0 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fdaini,ft)
    this%nf1 = ft(1); this%nf2 = ft(2); this%df = ft(3); this%id = ft(4) 
    this%ng = ft(5); this%numtasks = ft(6); this%nnod = ft(7)
    deallocate(ft)
    this%kfcur = 0
!    print *,'numtasks =',this%numtasks
!
!  open all other files needed and read in group tree
!
    allocate(this%root(this%numtasks),this%fda(this%numtasks))
    this%root(1) = rootini; this%fda(1) = fdaini
    call dealloc(rootini)
    do j = 2,this%numtasks
        write(crank,'(i3.3)') j-1
        filename = trim(basename)//'_'//crank
        lu = get(fuh)
        ierr = openFileStreamAccess(this%fda(j),lu,trim(filename))
        if (ierr /= 0) then
            call add(errmsg,2,trim(filename)//' can not be opened',myname)
            call add(fuh,lu); return
        endif
        call readGroupStreamAccess(this%root(j),this%fda(j))
    enddo
!
!  allocate space for kernels here
!  makes clearFrequency obsolete
!
    ntot = this%ng*this%nnod
    allocate(this%g(ntot,3,3))
    allocate(this%gstr(ntot,6,3))
    end subroutine initialReadGeminiKernelGreenTensor
!--------------------------------------------------------------------------------------
!> \brief Repeated read of Green tensor spectra for different frequency
!! modified to new convention that f = kf*df instead of (jf-1)*df before
!! therefore kf = jf-1
!
    subroutine readFrequencyGeminiKernelGreenTensor(this,kf,errmsg)
    type (gemini_kernel_green_tensor) :: this
    integer :: jf,kf
    type (error_message) :: errmsg
    integer :: j,ic,n,jr,nnod
    complex, dimension(:), pointer :: d
    type (data_stream_access), pointer :: dset
    type (group_stream_access), pointer :: group
    character (len=36) :: myname = 'readFrequencyGeminiKernelGreenTensor'
!
    nullify(d,dset,group)
!
    call addTrace(errmsg,myname)
!
    jf = kf+1
    if (jf < this%nf1 .or. jf > this%nf2) then
        call add(errmsg,2,'Frequency index out of range',myname)
        return
    endif
!
    this%kfcur = kf
    nnod = this%nnod
!
    do j = 1,this%numtasks
        do jr = j,nnod,this%numtasks
            do n = 1,3
                do ic = 1,9
                    call traversePathStreamAccess(this%root(j),1, &
                        & (/ jf-this%nf1+1,ic+(n-1)*9+(jr-j)/this%numtasks*27 /),group,dset)
                    if (ic <= 3 ) then
                        call readDatasetVectorStreamAccess(dset,this%fda(j),d)
                        this%g((jr-1)*this%ng+1:jr*this%ng,ic,n) = d
!                        call fillComplexVectorPointer(this%g(ic,n),d,(jr-1)*this%ng+1)
                        deallocate(d)
                    else
                        call readDatasetVectorStreamAccess(dset,this%fda(j),d)
                        this%gstr((jr-1)*this%ng+1:jr*this%ng,ic-3,n) = d
!                        call fillComplexVectorPointer(this%gstr(ic-3,n),d,(jr-1)*this%ng+1)
                        deallocate(d)
                    endif
                enddo
            enddo
        enddo
    enddo
    end subroutine readFrequencyGeminiKernelGreenTensor
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
    subroutine deallocGeminiKernelGreenTensor(this,fuh)
    type (gemini_kernel_green_tensor) :: this
    type (file_unit_handler) :: fuh
    integer :: j
    if (associated(this%root)) then
        do j = 1, this%numtasks
            call clearGroupTree(this%root(j))
        enddo
        deallocate(this%root)
    endif
    if (associated(this%fda)) then
        do j = 1, this%numtasks
            call add(fuh,getFileUnitStreamAccess(this%fda(j)))
            call dealloc(this%fda(j))
        enddo
        deallocate(this%fda)
    endif
    if (associated(this%g)) deallocate(this%g)
    if (associated(this%gstr)) deallocate(this%gstr)
    end subroutine deallocGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get ID
!
    function getIdGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    character(len=80) :: res
    res = this%id
    end function getIdGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get frequency spacing df
!
    function getDfGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: res
    res = this%df
    end function getDfGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
    function getNfGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nf2-this%nf1+1
    end function getNfGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get current frequency index kfcur
!
    function getJfcurGeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%kfcur
    end function getJfcurGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get strains
!
    function getStrainsGeminiKernelGreenTensor(this) result(str)
    type (gemini_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: str
    str => this%gstr
    end function getStrainsGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsGeminiKernelGreenTensor(this,uf_gstr)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: uf_gstr
    ! strains of gemini kernel green tensor are in units (1.0e-3 nm)/km (per N etc...), 
    ! i.e. according to length, they have the dimension-less unit 1.0e-15
    ! This results in very large numbers, which must be scaled by a factor of 1.0e-15 in order to remove this unit
    uf_gstr = 1.0e-15
  end subroutine getUnitFactorStrainsGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get green_tensors
!
    function getGeminiKernelGreenTensor(this) result(u)
    type (gemini_kernel_green_tensor), intent(in) :: this
    complex, dimension(:,:,:), pointer :: u
    u => this%g
    end function getGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get unit factor of green tensor
!
  subroutine getUnitFactorGeminiKernelGreenTensor(this,uf_g)
    type (gemini_kernel_green_tensor), intent(in) :: this
    real :: uf_g
    ! gemini kernel green tensor are in units (1.0e-3 nm) (per N etc...), 
    ! i.e. according to length, they have the unit 1.0e-12 m
    ! In order to get SI unit m, green tensor values must be scaled by a factor of 1.0e-12
    uf_g = 1.0e-12
  end subroutine getUnitFactorGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!! modified to new convention that f = kf*df instead of (jf-1)*df before
!! therefore kf = jf-1
!
    function nextFrequencyGeminiKernelGreenTensor(this,kf) result(next)
    type (gemini_kernel_green_tensor) :: this
    integer :: jf,kf
    logical :: next
    integer :: call_count = 0
    save call_count
!
    if (call_count == (this%nf2-this%nf1+1)) then
        next = .false.
        call_count = 0
        return
    endif
    call_count = call_count+1
    jf = this%nf1+call_count-1
    kf = jf-1
    next = .true.
    end function nextFrequencyGeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get first frequency nf1
!
    function getNf1GeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nf1
    end function getNf1GeminiKernelGreenTensor
!-------------------------------------------------------------------------
!> \brief Get last frequency nf2
!
    function getNf2GeminiKernelGreenTensor(this) result(res)
    type (gemini_kernel_green_tensor), intent(in) :: this
    integer :: res
    res = this%nf2
    end function getNf2GeminiKernelGreenTensor
!
 end module
