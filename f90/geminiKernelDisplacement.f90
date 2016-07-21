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
!------------------------------------------------------------------------------
!> \brief Module dealing with Gemini kernel displacements
!------------------------------------------------------------------------------
 module geminiKernelDisplacement
    use streamAccess
    use errorMessage
    use vectorPointer
    use fileUnitHandler
    implicit none
    interface dealloc; module procedure deallocGeminiKernelDisplacement; end interface
    interface operator (.id.); module procedure getIdGeminiKernelDisplacement; end interface
    interface operator (.nf.); module procedure getNfGeminiKernelDisplacement; end interface
    interface operator (.df.); module procedure getDfGeminiKernelDisplacement; end interface
    interface operator (.jfcur.); module procedure getJfcurGeminiKernelDisplacement; end interface
    interface operator (.disp.); module procedure getGeminiKernelDisplacement; end interface
    interface operator (.strain.); module procedure getStrainsGeminiKernelDisplacement; end interface
    interface operator (.ng.); module procedure getNgGeminiKernelDisplacement; end interface
    type gemini_kernel_displacement
        private
        character (len=80) :: id                                ! id for kernel displacement
        integer :: nf1,nf2                                      ! frequency index range
        real :: df                                              ! frequency stepping
        real :: kfcur                                           ! current frequency index: f = kf*df
        integer :: ng                                           ! total number of surface wavefield points
        integer :: nnod                                         ! nuber of depth nodes
        complex, dimension(:,:), pointer :: u                   ! array of 3 vector pointers for displacement
        complex, dimension(:,:), pointer :: ustr                ! array of 6 vector pointers for strain components on grid
        integer :: numtasks                                     ! number of parallel tasks that wrote kernel displacements
        type (file_stream_access), dimension(:), pointer :: fda         ! file pointers to all files
        type (group_stream_access), dimension(:), pointer :: root       ! root groups of all files
    end type
    integer, parameter :: length_ID_gemini_kernel_displacement = 13 !< change this number CONSISTENTLY with respective numbers in submodules (specfem3d,gemini,etc) 
!
 contains
!-----------------------------------------------------------------------------------
!> \brief Read in basic information needed for all frequencies and open parallel files
!
    subroutine initialReadGeminiKernelDisplacement(this,fuh,basename,errmsg)
    type (gemini_kernel_displacement) :: this
    type (file_unit_handler) :: fuh
    character (len=*) :: basename
    type (error_message) :: errmsg
    type (file_stream_access) :: fdaini
    type (group_stream_access) :: rootini
    type (data_stream_access), pointer :: dset
    type (group_stream_access), pointer :: group
    type (flexible), dimension(:), pointer :: ft
    integer :: j,ierr,lu,ntot
    character (len=400) :: filename
    character (len=3) :: crank
    character (len=35) :: myname = 'initialReadGeminiKernelDisplacement'
!
    call addTrace(errmsg,myname)
    lu = get(fuh)
    filename = trim(basename)//'_'//'000'
    ierr = openFileStreamAccess(fdaini,lu,trim(filename))
    if (ierr /= 0) then
        call add(errmsg,2,trim(filename)//' can not be opened',myname)
        call add(fuh,lu)
        return
    endif
!
!  read in group tree
!
    call readGroupStreamAccess(rootini,fdaini)
!    print *,'group tree read in from file: ',trim(filename)
!
!  read out header info: nf1,nf2,df,istyp,ntot,numtasks,nnod
!
    call traversePathStreamAccess(rootini,0,(/ 1,0 /),group,dset)
    call readDatasetVectorStreamAccess(dset,fdaini,ft)
    this%nf1 = ft(1); this%nf2 = ft(2); this%df = ft(3); this%id = ft(4)
    this%ng = ft(5); this%numtasks = ft(6); this%nnod = ft(7)
    deallocate(ft)
    this%kfcur = 0
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
            call add(fuh,lu)
            return
        endif
        call readGroupStreamAccess(this%root(j),this%fda(j))
    enddo
!
!  allocate space for displacements and strains
!  makes clearFrequency obsolete
!
    ntot = this%ng*this%nnod
    allocate(this%u(ntot,3))
    allocate(this%ustr(ntot,6))
    end subroutine initialReadGeminiKernelDisplacement
!----------------------------------------------------------------------------------------------------
!> \brief Deallocate object
!
    subroutine deallocGeminiKernelDisplacement(this,fuh)
    type (gemini_kernel_displacement) :: this
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
    if (associated(this%u)) deallocate(this%u)
    if (associated(this%ustr)) deallocate(this%ustr)
    end subroutine deallocGeminiKernelDisplacement
!-----------------------------------------------------------------------------------------------------
!> \brief Repeated read of displacement spectra file for different frequencies
!! modified to new convention that f = kf*df instead of (jf-1)*df before
!! therefore kf = jf-1
!
    subroutine readFrequencyGeminiKernelDisplacement(this,kf,errmsg)
    type (gemini_kernel_displacement) :: this
    integer :: jf,kf
    type (error_message) :: errmsg
    integer :: j,ic,jr,nnod
    complex, dimension(:), pointer :: d
    type (data_stream_access), pointer :: dset
    type (group_stream_access), pointer :: group
    character (len=37) :: myname = 'readFrequencyGeminiKernelDisplacement'
!
    call addTrace(errmsg,myname)
!
! check validity of frequency index
!
    jf = kf+1
    if (jf < this%nf1 .or. jf > this%nf2) then
        call add(errmsg,2,'Frequency index out of range ',myname)
        return
    endif
    this%kfcur = kf
    nnod = this%nnod
!
    do j = 1,this%numtasks
        do jr = j,nnod,this%numtasks
            do ic = 1,9
                call traversePathStreamAccess(this%root(j),1,(/ jf-this%nf1+1,ic+(jr-j)/this%numtasks*9 /),group,dset)
                if (ic <= 3 ) then
                    call readDatasetVectorStreamAccess(dset,this%fda(j),d)
                    this%u((jr-1)*this%ng+1:jr*this%ng,ic) = d
                    deallocate(d)
                else
                    call readDatasetVectorStreamAccess(dset,this%fda(j),d)
                    this%ustr((jr-1)*this%ng+1:jr*this%ng,ic-3) = d
                    deallocate(d)
                endif
            enddo
        enddo
    enddo
    end subroutine readFrequencyGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get ID
!
    function getIdGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    character(len=80) :: res
    res = this%id
    end function getIdGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get frequency spacing df
!
    function getDfGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: res
    res = this%df
    end function getDfGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of frequencies nf
!
    function getNfGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nf2-this%nf1+1
    end function getNfGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get current frequency index kfcur
!
    function getJfcurGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%kfcur
    end function getJfcurGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get strains
!
    function getStrainsGeminiKernelDisplacement(this) result(str)
    type (gemini_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: str
    str => this%ustr
    end function getStrainsGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of strains
!
  subroutine getUnitFactorStrainsGeminiKernelDisplacement(this,uf_ustr)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: uf_ustr
    ! strains of gemini kernel displacement are in units nm/km (per N or per Nm etc...), 
    ! i.e. according to length, they have the dimension-less unit 1.0e-12
    ! This results in very large numbers, which must be scaled by a factor of 1.0e-12 in order to remove this unit
    uf_ustr = 1.0e-12 
  end subroutine getUnitFactorStrainsGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get displacements
!
    function getGeminiKernelDisplacement(this) result(u)
    type (gemini_kernel_displacement), intent(in) :: this
    complex, dimension(:,:), pointer :: u
    u => this%u
    end function getGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get unit factor of displacements
!
  subroutine getUnitFactorGeminiKernelDisplacement(this,uf_u)
    type (gemini_kernel_displacement), intent(in) :: this
    real :: uf_u
    ! gemini kernel displacements are in units nm (per N or per Nm etc...), i.e. according to length, 
    ! they have the unit 1.0e-9 m
    ! In order to get SI unit m, displacement values must be scaled by a factor of 1.0e-9
    uf_u = 1.0e-9
  end subroutine getUnitFactorGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Iterator over frequencies
!!  modified for new convention that f = kf*df instead of f = (jf-1)*df
!!  hence kf = jf-1 or jf = kf+1
!
    function nextFrequencyGeminiKernelDisplacement(this,kf) result(next)
    type (gemini_kernel_displacement) :: this
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
    end function nextFrequencyGeminiKernelDisplacement
!-------------------------------------------------------------------------
!  Old stuff, needed anymore ??
!-------------------------------------------------------------------------
!> \brief Get displacements on a horizontal slice at depth index jr and component ic
!
    function getHsliceDispGeminiKernelDisplacement(this,jr,ic) result(res)
    type (gemini_kernel_displacement) :: this
    integer :: jr,ic
    complex, dimension(:), pointer :: res
    if (jr < 1 .or. jr > this%nnod) then
        print *,'<getHsliceDispGeminiKernelDisplacement>: Invalid depth index: ',jr 
        res => null(); return
    endif
    if (ic < 1 .or. ic > 3) then
        print *,'<getHsliceDispGeminiKernelDisplacement>: Invalid component: ',ic 
        res => null(); return
    endif
    res => this%u((jr-1)*this%ng+1:jr*this%ng,ic)
    end function getHsliceDispGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get strains on a horizontal slice at depth index jr and component ic
!
    function getHsliceStrainsGeminiKernelDisplacement(this,jr,ic) result(res)
    type (gemini_kernel_displacement) :: this
    integer :: jr,ic
    complex, dimension(:), pointer :: res
    if (jr < 1 .or. jr > this%nnod) then
        print *,'<getHsliceStrainsGeminiKernelDisplacement>: Invalid depth index: ',jr 
        res => null(); return
    endif
    if (ic < 1 .or. ic > 6) then
        print *,'<getHsliceStrainsGeminiKernelDisplacement>: Invalid component: ',ic 
        res => null(); return
    endif
    res => this%ustr((jr-1)*this%ng+1:jr*this%ng,ic)
    end function getHsliceStrainsGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get number of grid points per layer (this%ng)
!
    function getNgGeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%ng
    end function getNgGeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get first frequency nf1
!
    function getNf1GeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nf1
    end function getNf1GeminiKernelDisplacement
!-------------------------------------------------------------------------
!> \brief Get last frequency nf2
!
    function getNf2GeminiKernelDisplacement(this) result(res)
    type (gemini_kernel_displacement), intent(in) :: this
    integer :: res
    res = this%nf2
    end function getNf2GeminiKernelDisplacement
!
 end module
