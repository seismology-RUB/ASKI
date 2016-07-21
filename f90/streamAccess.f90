!--------------------------------------------------------------------------
!   Copyright 2015 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of Gemini II.
!   This file is part of ASKI version 1.0.
!
!   Gemini II and ASKI version 1.0 are free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   Gemini II and ASKI version 1.0 are distributed in the hope that they
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with Gemini II and ASKI version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!----------------------------------------------------------------
!> \brief Handle structured stream input and output

!> \author Wolfgang Friederich

!> \par Description
!>  Support for stream access file reading and writing.
!!  File tree starts with a root group that may contain further groups
!!  and datasets. Each subgroup can contain subgroups and datasets
!!  itself. Idea is to create groups and datasets and to build up a
!!  tree under the root group which contains links to all subgroups
!!  and datasets. Datasets are written to file, file positions are stored
!!  in the group and dataset objects. Finally, the complete tree
!!  with information about the data is also written to file in a way
!!  that by reading the file the complete information and data can
!!  be retrieved.
!<---------------------------------------------------------------
 module streamAccess
    use flexibleType
    use kindDefinitions
    implicit none
    interface createGroupStreamAccess
        module procedure createPostallocatedGroupStreamAccess
        module procedure createPreallocatedGroupStreamAccess
    end interface
    interface new
        module procedure createDatasetStreamAccess
    end interface
    interface dealloc
        module procedure deallocFileStreamAccess
        module procedure nullifyGroupStreamAccess
        module procedure fullyDeallocDatasetStreamAccess
    end interface
    interface clearGroupTree; module procedure recursiveDeallocGroupStreamAccess; end interface
    interface writeDatasetVectorStreamAccess
        module procedure writeDatasetIntegerVectorStreamAccess
        module procedure writeDatasetRealVectorStreamAccess
        module procedure writeDatasetDoubleVectorStreamAccess
        module procedure writeDatasetComplexVectorStreamAccess
        module procedure writeDatasetDoubleComplexVectorStreamAccess
        module procedure writeDatasetCharVectorStreamAccess
        module procedure writeDatasetFlexibleVectorStreamAccess
    end interface
    interface writeDataset2DArrayStreamAccess
        module procedure writeDatasetReal2DArrayStreamAccess
        module procedure writeDatasetDouble2DArrayStreamAccess
        module procedure writeDatasetComplex2DArrayStreamAccess
        module procedure writeDatasetDoubleComplex2DArrayStreamAccess
    end interface
    interface readDatasetVectorStreamAccess
        module procedure readDatasetIntegerVectorStreamAccess
        module procedure readDatasetRealVectorStreamAccess
        module procedure readDatasetDoubleVectorStreamAccess
        module procedure readDatasetComplexVectorStreamAccess
        module procedure readDatasetDoubleComplexVectorStreamAccess
        module procedure readDatasetCharVectorStreamAccess
        module procedure readDatasetFlexibleVectorStreamAccess
    end interface
    interface readDataset2DArrayStreamAccess
        module procedure readDatasetReal2DArrayStreamAccess
        module procedure readDatasetDouble2DArrayStreamAccess
        module procedure readDatasetComplex2DArrayStreamAccess
        module procedure readDatasetDoubleComplex2DArrayStreamAccess
    end interface
    interface operator (.ndataset.); module procedure getNdsGroupStreamAccess; end interface
    interface operator (.nsubgroup.); module procedure getNsubGroupStreamAccess; end interface
    interface operator (.groupname.); module procedure getNameGroupStreamAccess; end interface
    interface operator (.grouptag.); module procedure getTagGroupStreamAccess; end interface
!
    type file_stream_access
        private
        integer :: lu                            !< Fortran identifier of file
        integer (longint) :: current_file_pos    !< current position of file pointer
        character (len=132) :: filename          !< name of file
    end type
!
    type data_stream_access
        private
        integer (longint) :: filepos_start                 !< position where data start
        integer :: rank                                    !< dimensionality of data array
        integer, dimension(:), pointer :: dims =>  null()  !< number of data
        integer :: datatype                                !< integer coded primitive type of data
    end type
!
    type group_stream_access
        private
        integer (longint) :: filepos      !< position where to find subgroup and data information in file
        integer :: tag                    !< positive integer number indexing the group in some way (root has zero!)
        character (len=132) :: name       !< name specifying group properties
        type (group_stream_access), pointer :: parent => null()  !< pointer to parent
        type (group_stream_access), dimension(:), pointer :: subgroup => null()   !< pointer to subgroups
        integer :: maxsubgroup            !< maximum number of subgroups that shall be added to this group
                                          !! optional argument in createGroupStreamAccess
                                          !! if .not.present(maxsubgroup) in createGroupStreamAccess, it is set to -1 and is not used 
                                          !<
        integer :: lastsubgroup           !< index of subgroup that was added last (internal use)
        type (data_stream_access), dimension(:), pointer :: dataset => null()     !< pointer to datasets
        integer :: maxdset                !< maximum number of datasets that shall be added to this group
                                          !! optional argument in createGroupStreamAccess
                                          !! if .not.present(maxdset) in createGroupStreamAccess, it is set to -1 and is not used
                                          !<
        integer :: lastdset               !< index of dataset that was added last (internal use)
    end type
!
    logical :: verboseStreamAccess = .false.
!
 contains
!---------------------------------------------------------------------------
!> \brief Create a new stream access file object. Open file for writing
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
    integer function createFileStreamAccess(this,lu,filename) result(ios)
    type (file_stream_access) :: this
    integer :: lu
    character (len=*) :: filename
    logical :: exflag
!
    inquire(file=filename, exist = exflag)   ! delete file if it exists
    if (exflag) then
        open(lu,file=filename,access = 'stream', form='unformatted')
        close(lu,status = 'delete')
    endif
    open(lu,file=filename,form='unformatted',access='stream',status='new',iostat=ios)
    if (ios /= 0) return
    this%lu = lu
    this%current_file_pos = 1+bit_size(this%current_file_pos)/8      !  leave space for a long integer
    this%filename = filename
    end function createFileStreamAccess
!-----------------------------------------------------------------
!> \brief Open an existing stream access file for reading
!> \param this file_stream_access object
!> \param lu Fortran file identifier
!> \param filename Name of stream access file
!
    integer function openFileStreamAccess(this,lu,filename) result(ios)
    type (file_stream_access) :: this
    integer :: lu
    character (len=*) :: filename
    open(lu,file=filename,form='unformatted',access='stream',status='old',iostat=ios)
    if (ios /= 0) return
    this%lu = lu
    read(lu,pos = 1) this%current_file_pos
    end function openFileStreamAccess
!-----------------------------------------------------------------
!> \brief Deallocate a direct access file object and close file
!> \param this file_stream_access object
!
    subroutine deallocFileStreamAccess(this)
    type (file_stream_access) :: this
    close(this%lu)
    end subroutine deallocFileStreamAccess
!----------------------------------------------------------------
!> \brief Get the logical unit of the file
!> \param this file_stream_access object
!
    integer function getFileUnitStreamAccess(this)
    type (file_stream_access) :: this
    getFileUnitStreamAccess = this%lu
    end function getFileUnitStreamAccess
!----------------------------------------------------------------
!> \brief Get the current file position
!> \param this file_stream_access object
!
    function getCurrentFilePositionStreamAccess(this) result(filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    filepos = this%current_file_pos
    end function getCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Get the file name
!> \param this file_stream_access object
!
    function getFileNameStreamAccess(this) result(fname)
    type (file_stream_access) :: this
    character(len=132) :: fname
    fname = this%filename
    end function getFileNameStreamAccess
!----------------------------------------------------------------
!> \brief Increment current file position
!> \param this file_stream_access object
!> \param n number of bytes to move
!
    function IncrementCurrentFilePositionStreamAccess(this,n) result(filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    integer :: n
    this%current_file_pos = this%current_file_pos + n
    filepos = this%current_file_pos
    end function IncrementCurrentFilePositionStreamAccess
!----------------------------------------------------------------
!> \brief Set current file position
!> \param this file_stream_access object
!> \param filepos desired position of file pointer
!
    subroutine setCurrentFilePositionStreamAccess(this,filepos)
    type (file_stream_access) :: this
    integer (longint) :: filepos
    this%current_file_pos = filepos
    end subroutine setCurrentFilePositionStreamAccess
!---------------------------------------------------------------------
!> \brief Create a group. Contents to be filled in later.
!> \param this group_stream_access object
!> \param name name of the group
!> \param tag an identifier of the group
!> \par
!> Maximal number of subgroups and datasets need not to be known a priori, 
!! arrays of subgroups and datasets reallocated frequently
!>
!
    subroutine createPostallocatedGroupStreamAccess(this,name,tag)
    type (group_stream_access) :: this
    character (len=*) :: name
    integer :: tag
!
    this%name = name
    this%tag = tag
    this%maxsubgroup = -1
    this%lastsubgroup = -1
    this%maxdset = -1
    this%lastdset = -1
    if (verboseStreamAccess) then
        print *,'create group named ',trim(this%name),' Tag = ',this%tag
    endif
!
    end subroutine createPostallocatedGroupStreamAccess
!---------------------------------------------------------------------
!> \brief Create a group. Contents to be filled in later. 
!> \param this group_stream_access object
!> \param name name of the group
!> \param tag an identifier of the group
!> \param maxsubgroup number of subgroups that are accounted for in array preallocation
!> \param maxdset number of datasets that are accounted for in array preallocation
!> \par
!> Maximal number of subgroups and datasets are know a priori, 
!! arrays of subgroups and datasets are allocated in advance 
!! to avoid frequent reallocation.
!>
!
    subroutine createPreallocatedGroupStreamAccess(this,name,tag,maxsubgroup,maxdset)
    type (group_stream_access) :: this
    character (len=*) :: name
    integer :: tag
    integer :: maxsubgroup
    integer :: maxdset
    integer :: ierr
!
    this%name = name
    this%tag = tag
    if(maxsubgroup > 0) then
        this%maxsubgroup = maxsubgroup
        allocate(this%subgroup(1:this%maxsubgroup), stat=ierr)
        if(ierr /= 0) stop "allocate error"
    else
        this%maxsubgroup = 0
    endif
    this%lastsubgroup = 0
    if(maxdset > 0) then
        this%maxdset = maxdset
        allocate(this%dataset(1:this%maxdset), stat=ierr)
        if(ierr /= 0) stop "allocate error"
    else
        this%maxdset = 0
    endif
    this%lastdset = 0
    if (verboseStreamAccess) then
        print *,'create group named ',trim(this%name),' Tag = ',this%tag, & 
                        ' preallocated ',this%maxsubgroup,' subgroups and ',this%maxdset,' datasets'
    endif
!
    end subroutine createPreallocatedGroupStreamAccess
!---------------------------------------------------------------------
!> \brief Recursively deallocate group including subgroups and datasets
!> \param this group_stream_access object
!
    recursive subroutine recursiveDeallocGroupStreamAccess(this)
    type (group_stream_access) :: this
    integer :: j
    if (associated(this%subgroup)) then
        do j=1,size(this%subgroup)
            call recursiveDeallocGroupStreamAccess(this%subgroup(j))
        enddo
        deallocate(this%subgroup)
    endif
    if (associated(this%dataset)) then
        do j=1,size(this%dataset)
            call fullyDeallocDatasetStreamAccess(this%dataset(j))
        enddo
        deallocate(this%dataset)
    endif
    end subroutine recursiveDeallocGroupStreamAccess
!-----------------------------------------------------------------------
!> \brief Only nullify pointers in group object. Do not deallocate.
!> \param this group_stream_access object
!
    subroutine nullifyGroupStreamAccess(this)
    type (group_stream_access) :: this
    if (associated(this%parent)) nullify(this%parent)
    if (associated(this%subgroup)) nullify(this%subgroup)
    if (associated(this%dataset)) nullify(this%dataset)
    end subroutine nullifyGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a subgroup to this group.
!> \param this group_stream_access object
!> \param group subgroup object to be incorporated into the tree
!> \par
!> The subgroup must be completely filled with contents.
!! Only a SHALLOW copy of the group object into the tree is performed.
!! Avoid copying the same group into different trees, because
!! clearing one tree leads to a deallocation of memory still associated
!! with pointers in other trees. Clearing these trees will then fail. 
!! After calling this routine \a group may be deallocated using dealloc. 
!<
    subroutine addSubgroupStreamAccess(this,group)
    type (group_stream_access), target :: this
    type (group_stream_access) :: group
    integer :: n,j
!
    if (associated(this%subgroup)) then
        n = size(this%subgroup)
    else
        n = 0
    endif
    if (this%maxsubgroup == -1) then  ! case of frequent reallocation
        this%subgroup => reallocateGroupStreamAccess(this%subgroup,n+1)
        this%subgroup(n+1) = group
        if (verboseStreamAccess) then
            print *,'add group named ',trim(this%subgroup(n+1)%name),' as ',n+1,' th subgroup of group ', &
                  & trim(this%name),' with tag = ',this%tag
        endif
        this%subgroup(n+1)%parent => this
    !
    !  redirect the parent member of group's subgroups from group to this%subgroup(n+1)
    !
        if (associated(group%subgroup)) then
            do j=1,size(group%subgroup)
                group%subgroup(j)%parent => this%subgroup(n+1)
            enddo
        endif
    else  ! case of preallocation, only reallocate when maxsubgroup is reached
        if (this%lastsubgroup == n) &
            this%subgroup => reallocateGroupStreamAccess(this%subgroup,n + this%maxsubgroup/10 + 1) !"+1" for case maxsubgroup < 10
        this%subgroup(this%lastsubgroup+1) = group
        if (verboseStreamAccess) then
            print *,'add group named ',trim(this%subgroup(this%lastsubgroup+1)%name),' as ',this%lastsubgroup+1, &
                  & ' th subgroup of group ',trim(this%name),' with tag = ',this%tag
        endif
        this%subgroup(this%lastsubgroup+1)%parent => this
    !
    !  redirect the parent member of group's subgroups from group to this%subgroup(lastsubgroup+1)
    !
        if (group%lastsubgroup > 0) then
            do j=1,group%lastsubgroup
                group%subgroup(j)%parent => this%subgroup(this%lastsubgroup+1)
            enddo
        endif
        this%lastsubgroup = this%lastsubgroup + 1
    endif
    end subroutine addSubgroupStreamAccess
!----------------------------------------------------------------------
!> \brief Add a dataset to some existing group
!> \param this group_stream_access object
!> \param dset dataset object
!> \par 
!> The dataset object must contain all information. This implies
!! that the data have already been written to file
!<
    subroutine addDatasetStreamAccess(this,dset)
    type (group_stream_access) :: this
    type (data_stream_access) :: dset
    integer :: n
!
    if (associated(this%dataset)) then
        n = size(this%dataset)
    else
        n = 0
    endif
    if (this%maxdset == -1) then
        this%dataset => reallocateDatasetStreamAccess(this%dataset,n+1)
        call createDatasetStreamAccess(this%dataset(n+1),dset%rank,dset%dims,dset%datatype)
        this%dataset(n+1)%filepos_start = dset%filepos_start
        if (verboseStreamAccess) then
            print *,'add ',n+1,' th dataset with size ',this%dataset(n+1)%dims,' to group ',trim(this%name),' with tag ',this%tag
        endif
    else
        if (this%lastdset == n) &
            this%dataset => reallocateDatasetStreamAccess(this%dataset,n + this%maxdset/10 + 1) !"+1" for case maxdset < 10
        call createDatasetStreamAccess(this%dataset(this%lastdset+1),dset%rank,dset%dims,dset%datatype)
        this%dataset(this%lastdset+1)%filepos_start = dset%filepos_start
        if (verboseStreamAccess) then
            print *,'add ',this%lastdset+1,' th dataset with size ',this%dataset(this%lastdset+1)%dims,' to group ', &
                & trim(this%name),' with tag ',this%tag
        endif
        this%lastdset = this%lastdset + 1
    endif
    end subroutine addDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of datasets contained in group
!> \param this group_stream_access object
!
    integer function getNdsGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    if (associated(this%dataset)) then
        if (this%maxdset == -1) then
            getNdsGroupStreamAccess = size(this%dataset)
        else
            getNdsGroupStreamAccess = this%lastdset
        endif
    else
        getNdsGroupStreamAccess = 0
    endif
    end function getNdsGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get number of subgroups contains in group
!> \param this group_stream_access object
!
    integer function getNsubGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    if (associated(this%subgroup)) then
        if (this%maxsubgroup == -1) then
            getNsubGroupStreamAccess = size(this%subgroup)
        else
            getNsubGroupStreamAccess = this%lastsubgroup
        endif
    else
        getNsubGroupStreamAccess = 0
    endif
    end function getNsubGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get tag of group
!> \param this group_stream_access object
!
    integer function getTagGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    getTagGroupStreamAccess = this%tag
    end function getTagGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get name of group
!> \param this group_stream_access object
!
    character (len=132) function getNameGroupStreamAccess(this)
    type (group_stream_access), intent(in) :: this
    getNameGroupStreamAccess = this%name
    end function getNameGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Get a pointer to the parent group
!> \param this group_stream_access object
!> \return a pointer to the parent group
!
    function getParentStreamAccess(this) result(parent)
    type (group_stream_access) :: this
    type (group_stream_access), pointer :: parent
    parent => this%parent
    end function getParentStreamAccess
!----------------------------------------------------------------------
!> \brief Get pointer to selected dataset in group
!> \param this group_stream_access object
!> \param k index of selected dataset in group
!> \return a pointer to a data_stream_access object
!
    function getDatasetSelectedStreamAccess(this,k) result(dset)
    type (group_stream_access) :: this
    type (data_stream_access), pointer:: dset
    integer :: k
!
    dset => this%dataset(k)
    end function getDatasetSelectedStreamAccess
!----------------------------------------------------------------------
!> \brief Recursively write the group data to file
!> \param this file_stream_access object
!> \param fda file_stream_access object
!
    recursive subroutine writeGroupStreamAccess(this,fda)
    type (group_stream_access) :: this
    type (file_stream_access) :: fda
    integer :: nsub,nds,j,lu,move
    integer (longint) :: filepos
!
    if (.nsubgroup.this > 0) then    ! first work through my subgroups
        do j=1,.nsubgroup.this
            call writeGroupStreamAccess(this%subgroup(j),fda)
        enddo
    endif
!
!  now deal with myself
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
!
    this%filepos = filepos           ! position where infos about my group start
                                     ! this info will later be written to file by my parent group
!
!  if I am the root group I write now the first position of the file
!  specifiying the position where info about the root group's contents has been written
!
    if (.not. associated(this%parent)) write(lu,pos = 1) filepos
!
    nsub = getNsubGroupStreamAccess(this)
    nds = getNdsGroupStreamAccess(this)
!
    write(lu,pos=filepos) this%name,this%tag,nsub,nds   ! write name, tag etc
    if (verboseStreamAccess) then
        print *,'write group: ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
    endif
    filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
!
    if (nsub > 0) then                                  ! write subgroup record info
        write(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
        move = nsub*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    endif
    if (nds > 0) then                                   ! write dataset record info
        write(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
        move = nds*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    endif
    end subroutine writeGroupStreamAccess
!------------------------------------------------------------------------------------
!> \brief Recursively read group data from file into memory starting with root group
!> \param this group_stream_access object
!> \param fda file from which group data read
!
    recursive subroutine readGroupStreamAccess(this,fda,readInfoDset)
    type (group_stream_access), target :: this
    type (file_stream_access), target :: fda
    logical, optional :: readInfoDset
    logical :: readInfoDsetFlag
    integer :: lu,nsub,nds,j,move
    integer (longint) :: filepos
!
    if(present(readInfoDset)) then; readInfoDsetFlag=readInfoDset; else; readInfoDsetFlag = .true.;endif
!
    lu = getFileUnitStreamAccess(fda)
    filepos = getCurrentFilePositionStreamAccess(fda)
    this%filepos = filepos
    read(lu,pos=filepos) this%name,this%tag,nsub,nds      ! read first group record
    if (verboseStreamAccess) then
        print *,'read group ',trim(this%name),' with tag ',this%tag,' nsub = ',nsub,' nds = ',nds
    endif
    filepos = IncrementCurrentFilePositionStreamAccess(fda,len(this%name)+3*kind(1))
    this%maxsubgroup = nsub; this%lastsubgroup = nsub
    this%maxdset = nds ; this%lastdset = nds
!
    if (nsub > 0) then                                   ! read subgroup record info
        allocate(this%subgroup(nsub))
        do j=1,nsub
            this%subgroup(j)%parent => this
        enddo
        read(lu,pos=filepos) (this%subgroup(j)%filepos,j=1,nsub)
        move = nsub*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    else
        this%subgroup => null()
    endif
    if (nds > 0) then                                    ! read dataset record info
        allocate(this%dataset(nds))
        read(lu,pos=filepos) (this%dataset(j)%filepos_start,j=1,nds)
        move = nds*bit_size(filepos)/8
        filepos = IncrementCurrentFilePositionStreamAccess(fda,move)
    else
        this%dataset => null()
    endif
    if (nsub > 0) then                                   !  read in subgroups
        do j=1,nsub
            call setCurrentFilePositionStreamAccess(fda,this%subgroup(j)%filepos)
            call readGroupStreamAccess(this%subgroup(j),fda,readInfoDsetFlag)   ! read in subgroups
        enddo
    endif
    if (nds > 0 .and. readInfoDsetFlag) then
        do j = 1,nds
!            call setCurrentFilePositionStreamAccess(fda,this%dataset(j)%filepos_start)
            call readInfoDataStreamAccess(this%dataset(j),fda)
        enddo
    endif
    end subroutine readGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a group array
!> \param p pointer to an array of groups
!> \param n number of group objects to be allocated
!> \return newp pointer to the new array of groups
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
    function reallocateGroupStreamAccess(p, n) result(newp)
    type (group_stream_access), pointer, dimension(:) :: p, newp
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(newp(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    newp(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocateGroupStreamAccess
!----------------------------------------------------------------------
!> \brief Traverse group tree along some given path down to some dataset.
!> \param this Path refers to this group
!> \param depth gives level of the group containing the dataset as seen from \a this
!> \param path an index array of indices specifiying the path through the group tree, 
!!  the last index identifies the data set in the group.
!<
!> \param group output pointer to group the index array points to
!> \param dset a pointer to the dataset specified by the path
!> \par
!!  If \a group does not exist, return a null pointer.
!!  If \a dset does not exist, return a null pointer.
!!  The group pointer is returned to allow fast access to other datasets in the group.
!<
    recursive subroutine traversePathStreamAccess(this,depth,path,group,dset)
    type (group_stream_access), target :: this
    integer :: depth
    integer, dimension(:) :: path
    integer, dimension(:), allocatable :: path_copy
    type (group_stream_access), pointer :: group
    type (data_stream_access), pointer :: dset
    integer :: j,depth2
!
    allocate(path_copy(size(path)))
    path_copy = path
    if (depth > 0 ) then
        if (associated(this%subgroup)) then
            j = path_copy(1)                         ! take index from path
            path_copy = cshift(path_copy,1)          ! put first index at end
            depth2 = depth-1                         ! decrease depth
            call traversePathStreamAccess(this%subgroup(j),depth2,path_copy,group,dset)
        else
            group => null()       ! depth can't be reached
            dset => null()
        endif
    else
        j = path_copy(1)
        group => this
        if (associated(this%dataset) .and. j > 0) then
            if (size(this%dataset) >= j) then
                dset => this%dataset(j)
            else
                dset => null()
            endif
        else
            dset => null()
        endif
    endif
    deallocate(path_copy)
    end subroutine traversePathStreamAccess
!----------------------------------------------------------------------
!> \brief Create a dataset
!> \param this a data_stream_access object
!> \param rank rank of data array
!> \param dims array of integers specifiying size in each rank
!> \param datatype type of data in array using conventions in \link primitiveTypeEncoding.f90 primitiveTypeEncoding \endlink
!> \par
!> Specify rank, dims and datatype but nothing else
!! which is done later on writing to file
!<
    subroutine createDatasetStreamAccess(this,rank,dims,datatype)
    type (data_stream_access) :: this
    integer :: rank
    integer, dimension(:) :: dims
    integer :: datatype
!
    this%rank = rank
    this%datatype = datatype
    allocate(this%dims(rank))
    this%dims = dims
    if (verboseStreamAccess) then
        print *,'create dataset of size ',this%dims
    endif
!
    end subroutine createDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Fully deallocate a dataset
!> \param this data_stream_access object
!
    subroutine fullyDeallocDatasetStreamAccess(this)
    type (data_stream_access) :: this
    if (associated(this%dims)) deallocate(this%dims)
    end subroutine fullyDeallocDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Nullify the pointers contained in a dataset object
!> \param this data_stream_access object
!
    subroutine nullifyDatasetStreamAccess(this)
    type (data_stream_access) :: this
    if (associated(this%dims)) nullify(this%dims)
    end subroutine nullifyDatasetStreamAccess
!----------------------------------------------------------------------
!> \brief Reallocate a dataset pointer array
!> \param p pointer to an array of datasets
!> \param n number of dataset objects to be allocated
!> \return newp pointer to the new array of datasets
!> \par
!> The contents of the old array \a p is copied into the new one.
!! Afterwards, the old array \a p is deallocated.
!<
!--------------------------------------------------------------------
    function reallocateDatasetStreamAccess(p, n) result(newp)
    type (data_stream_access), pointer, dimension(:) :: p, newp
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(newp(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    newp(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocateDatasetStreamAccess
!---------------------------------------------------------------------
!> \brief Read info part of dataset (without data)
!
    subroutine readInfoDataStreamAccess(this,fda)
    type (data_stream_access) :: this
    type (file_stream_access) :: fda
    integer :: lu
    integer (longint) :: filepos
!
    lu = getFileUnitStreamAccess(fda)
!    filepos = getCurrentFilePositionStreamAccess(fda)
!    this%filepos_start = filepos
    filepos = this%filepos_start
    read(lu,pos = filepos) this%rank                  ! read dataset info
    allocate(this%dims(this%rank))
    read(lu,pos = filepos+kind(1)) this%dims,this%datatype
    call setCurrentFilePositionStreamAccess(fda,filepos+kind(1)*(1+this%rank+1))
!    filepos = IncrementCurrentFilePositionStreamAccess(fda,kind(1)*(1+this%rank+1))
    end subroutine readInfoDataStreamAccess
!----------------------------------------------------------------------
!> \brief Return size of data vector in data set
!
    function getVectorSizeDataStreamAccess(this) result(n)
    type (data_stream_access) :: this
    integer :: n
    n = this%dims(1)
    end function getVectorSizeDataStreamAccess
!----------------------------------------------------------------------
!> \brief  Write an integer vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d integer data array
!
    subroutine writeDatasetIntegerVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    integer, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1+size(d))*kind(1))
    end subroutine writeDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a real vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d real data array
!
    subroutine writeDatasetRealVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.0))
    end subroutine writeDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data array
!
    subroutine writeDatasetDoubleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.d0))
    end subroutine writeDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d complex data array
!
    subroutine writeDatasetComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.0))
    end subroutine writeDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double complex vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double complex data array
!
    subroutine writeDatasetDoubleComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.d0))
    end subroutine writeDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a character vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d character data array
!
    subroutine writeDatasetCharVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    character(len=*), dimension(:) :: d
    integer :: clen
    integer :: lu
    integer (longint) :: filepos
!
    clen = len(d(1))
    if (clen > 80) then
        print *,'writeDatasetCharVectorStreamAccess: string length greater 80 !'
        stop
    endif
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,clen,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1)+size(d)*clen)
    end subroutine writeDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a flexible type vector to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d flexible data array
!
    subroutine writeDatasetFlexibleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    type (flexible), dimension(:) :: d
    integer :: lu,nbytes,j
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,T_FLEXIBLE
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    do j=1,size(d)
        call writeSAFlexibleType(d(j),lu,filepos,nbytes)
        filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
    enddo
    end subroutine writeDatasetFlexibleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Write a real 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d real data 2D-array
!
    subroutine writeDatasetReal2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.0))
    end subroutine writeDatasetReal2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data 2D-array
!
    subroutine writeDatasetDouble2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*kind(1.d0))
    end subroutine writeDatasetDouble2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a complex 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d complex data 2D-array
!
    subroutine writeDatasetComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.0))
    end subroutine writeDatasetComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Write a double complex 2D array to file
!> \param this data_stream_access object
!> \param fda file to be written
!> \param d double precision data 2D-array
!
    subroutine writeDatasetDoubleComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:,:) :: d
    integer :: lu
    integer (longint) :: filepos
!
    filepos = getCurrentFilePositionStreamAccess(fda)
    lu = getFileUnitStreamAccess(fda)
    this%filepos_start = filepos
    write(lu,pos=filepos) this%rank,this%dims,this%datatype,d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1)+size(d)*2*kind(1.d0))
    end subroutine writeDatasetDoubleComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read an integer vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d integer data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetIntegerVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    integer, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetIntegerVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_INTEGER) &
        print *, "WARNING in readDatasetIntegerVectorStreamAccess: this dataset is NOT ",&
        "of type integer, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1))
    end subroutine readDatasetIntegerVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a real vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d real data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetRealVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetRealVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_REAL) &
        print *, "WARNING in readDatasetRealVectorStreamAccess: this dataset is NOT ",&
        "of type real, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.0))
    end subroutine readDatasetRealVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE) &
        print *, "WARNING in readDatasetDoubleVectorStreamAccess: this dataset is NOT ",&
        "of type double, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.d0))
    end subroutine readDatasetDoubleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetComplexVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_COMPLEX) &
        print *, "WARNING in readDatasetComplexVectorStreamAccess: this dataset is NOT ",&
        "of type complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.0))
    end subroutine readDatasetComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double complex vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleComplexVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleComplexVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE_COMPLEX) &
        print *, "WARNING in readDatasetDoubleComplexVectorStreamAccess: this dataset is NOT ",&
        "of type double complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.d0))
    end subroutine readDatasetDoubleComplexVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a character vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d character data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetCharVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    character (len=80), dimension(:), pointer :: d
    integer :: clen,i,lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype,clen
    else
        read(lu,pos=filepos) rank,dims,datatype,clen
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetCharVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_CHAR) &
        print *, "WARNING in readDatasetCharVectorStreamAccess: this dataset is NOT ",&
        "of type character, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+2)*kind(1))
    allocate(d(this%dims(1)))
    read(lu,pos=filepos) (d(i)(1:clen),i=1,size(d))
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*clen)
    end subroutine readDatasetCharVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a flexible type vector from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d flexible data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetFlexibleVectorStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    type (flexible), dimension(:), pointer :: d
    integer :: lu,nbytes,j
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(1)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(1))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetFlexibleVectorStreamAccess: values rank,dims(1),datatype = ",&
            rank,dims(1),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%datatype = ",this%rank,this%dims(1),this%datatype," !! This means ",&
            "that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_FLEXIBLE) &
        print *, "WARNING in readDatasetFlexibleVectorStreamAccess: this dataset is NOT ",&
        "of flexible type, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1)))
    do j=1,size(d)
        call readSAFlexibleType(d(j),lu,filepos,nbytes)
        filepos = IncrementCurrentFilePositionStreamAccess(fda,nbytes)
    enddo
    end subroutine readDatasetFlexibleVectorStreamAccess
!----------------------------------------------------------------------
!> \brief Read a real 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetReal2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    real, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetReal2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_REAL) &
        print *, "WARNING in readDatasetReal2DArrayStreamAccess: this dataset is NOT ",&
        "of type real, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.0))
    end subroutine readDatasetReal2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDouble2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double precision, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDouble2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE) &
        print *, "WARNING in readDatasetDouble2DArrayStreamAccess: this dataset is NOT ",&
        "of type double, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*kind(1.d0))
    end subroutine readDatasetDouble2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a complex 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d complex data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    complex, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetComplex2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_COMPLEX) &
        print *, "WARNING in readDatasetComplex2DArrayStreamAccess: this dataset is NOT ",&
        "of type complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.0))
    end subroutine readDatasetComplex2DArrayStreamAccess
!----------------------------------------------------------------------
!> \brief Read a double complex 2D array from file
!> \param this data_stream_access object
!> \param fda file from which data are read
!> \param d double precision data array pointer (output, allocated in the subroutine)
!
    subroutine readDatasetDoubleComplex2DArrayStreamAccess(this,fda,d)
    type (data_stream_access) :: this
    type (file_stream_access), target :: fda
    double complex, dimension(:,:), pointer :: d
    integer :: lu
    integer (longint) :: filepos
    integer :: rank,datatype
    integer :: dims(2)
!
    lu = getFileUnitStreamAccess(fda)
    filepos = this%filepos_start
    call setCurrentFilePositionStreamAccess(fda,filepos)
    if(.not.associated(this%dims)) then
        allocate(this%dims(2))
        read(lu,pos=filepos) this%rank,this%dims,this%datatype
    else
        read(lu,pos=filepos) rank,dims,datatype
        if(rank/=this%rank.or.dims(1)/=this%dims(1).or.dims(2)/=this%dims(2).or.datatype/=this%datatype) &
            print *, "WARNING in readDatasetDoubleComplex2DArrayStreamAccess: values rank,dims(1),dims(2),datatype = ",&
            rank,dims(1),dims(2),datatype," read from file do not match allocated data_stream_access object ",&
            "dset%rank,dset%dims(1),dset%dims(2),dset%datatype = ",this%rank,this%dims(1),this%dims(2),&
            this%datatype," !! This means that dataset information held by root group differs from file content."
    endif
    if(this%datatype /= T_DOUBLE_COMPLEX) &
        print *, "WARNING in readDatasetDoubleComplex2DArrayStreamAccess: this dataset is NOT ",&
        "of type double complex, but of type (using primitive type encoding): ",this%datatype,". ",&
        "Do not trust the data returned by this subroutine!"
    filepos = IncrementCurrentFilePositionStreamAccess(fda,(1+this%rank+1)*kind(1))
    allocate(d(this%dims(1),this%dims(2)))
    read(lu,pos=filepos) d
    filepos = IncrementCurrentFilePositionStreamAccess(fda,size(d)*2*kind(1.d0))
    end subroutine readDatasetDoubleComplex2DArrayStreamAccess
!--------------------------------------------------------------------------------
!> \brief Print Info about a group
!
    subroutine printGroupInfoStreamAccess(this)
    type (group_stream_access) :: this
    write(*,'(80(1h-))')
    write(*,'(a,i4,a,i5,a,i5,a,a,a)') '|  Group: ',this%tag,'|  N-Subgroups: ',&
         .nsubgroup.this,'|  N-Datasets: ',.ndataset.this,'|  Name: ',trim(this%name),'|'
    write(*,'(80(1h-))')
    end subroutine printGroupInfoStreamAccess
!--------------------------------------------------------------------------------
!> \brief Print Info about a dataset
!
    subroutine printDsetInfoStreamAccess(this,idx)
    type (data_stream_access) :: this
    integer :: idx
    character(len=29) :: fmt
    write(fmt,'(a25,i1,a3)') '(a2,i7,a2,i5,a2,i5,a2,1x,',this%rank,'i5)'
    write(*,fmt) '| ',idx,' |',this%rank,' |',this%datatype,' |',this%dims
    end subroutine printDsetInfoStreamAccess
!----------------------------------------------------------------------
!> \brief Recursively write the file hierarchy to screen
!> \param this file_stream_access object
!> \param fda file_stream_access object
!
    recursive subroutine printGroupStreamAccess(this)
    type (group_stream_access) :: this
    integer :: nds,j
!
!  first deal with myself
!
    call printGroupInfoStreamAccess(this)
    nds = getNdsGroupStreamAccess(this)
    if (nds > 0) then
       write(*,'(a11,2a7,a12)') '| Dataset |',' Rank |',' Type |',' Dimensions |'
       do j = 1,nds
          call printDsetInfoStreamAccess(this%dataset(j),j)
       enddo
    endif
!
    if (.nsubgroup.this > 0) then    ! work through my subgroups
        do j=1,.nsubgroup.this
            call printGroupStreamAccess(this%subgroup(j))
        enddo
    endif
    end subroutine printGroupStreamAccess
!
 end module
