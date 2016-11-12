!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.2.
!
!   ASKI version 1.2 is free software: you can
!   redistribute it and/or modify it under the terms of the GNU
!   General Public License as published by the Free Software
!   Foundation, either version 2 of the License, or (at your option) 
!   any later version.
!
!   ASKI version 1.2 is distributed in the hope that it
!   will be useful, but WITHOUT ANY WARRANTY; without even the implied
!   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.2.
!--------------------------------------------------------------
!  enlarge or downsize an array in memory
!--------------------------------------------------------------
 module realloc
    interface reallocate
        module procedure reallocate_real
        module procedure reallocate_double_precision
        module procedure reallocate_complex
        module procedure reallocate_char
        module procedure reallocate_int
        module procedure reallocate_logical
        module procedure reallocate_double_precision_2d
        module procedure reallocate_integer_2d
        module procedure reallocate_real_2d
        module procedure reallocate_logical_2d
        module procedure reallocate_complex_2d
    end interface reallocate
 contains
!--------------------------------------------------------------------
    function reallocate_real(p, n)               ! reallocate REAL
    real, pointer, dimension(:) :: p, reallocate_real
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_real(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_real(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_real
!--------------------------------------------------------------------
    function reallocate_double_precision(p, n)    ! reallocate DOUBLE PRECISION
    double precision, pointer, dimension(:) :: p, reallocate_double_precision
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_double_precision(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_double_precision(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_double_precision
!--------------------------------------------------------------------
    function reallocate_complex(p, n)               ! reallocate COMPLEX
    complex, pointer, dimension(:) :: p, reallocate_complex
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_complex(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_complex(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_complex
!--------------------------------------------------------------------
    function reallocate_char(p, n)               ! reallocate char
    character (len=*), pointer, dimension(:) :: p
    character (len=len(p)), pointer, dimension(:) :: reallocate_char
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_char(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_char(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_char
!--------------------------------------------------------------------
    function reallocate_int(p, n)               ! reallocate int
    integer, pointer, dimension(:) :: p, reallocate_int
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_int(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_int(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_int
!---------------------------------------------------------------------------
    function reallocate_logical(p, n)               ! reallocate logical
    logical, pointer, dimension(:) :: p, reallocate_logical
    integer, intent(in) :: n
    integer :: nold, ierr
    allocate(reallocate_logical(1:n), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p), n)
    reallocate_logical(1:nold) = p(1:nold)
    deallocate(p)
    end function reallocate_logical
!--------------------------------------------------------------------
!     2D
!--------------------------------------------------------------------
    function reallocate_double_precision_2d(p, n, k)    ! reallocate DOUBLE PRECISION 2D
    double precision, pointer, dimension(:,:) :: p, reallocate_double_precision_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_double_precision_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_double_precision_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_double_precision_2d
!-----------------------------------------------------------------------
    function reallocate_real_2d(p, n, k)    ! reallocate real 2D
    real, pointer, dimension(:,:) :: p, reallocate_real_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_real_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_real_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_real_2d
!-----------------------------------------------------------------------
    function reallocate_integer_2d(p, n, k)    ! reallocate integer 2D
    integer, pointer, dimension(:,:) :: p, reallocate_integer_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_integer_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_integer_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_integer_2d
!-----------------------------------------------------------------------
    function reallocate_logical_2d(p, n, k)    ! reallocate logical 2D
    logical, pointer, dimension(:,:) :: p, reallocate_logical_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_logical_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_logical_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_logical_2d
!-----------------------------------------------------------------------
    function reallocate_complex_2d(p, n, k)    ! reallocate complex 2D
    complex, pointer, dimension(:,:) :: p, reallocate_complex_2d
    integer, intent(in) :: n,k
    integer :: nold, kold, ierr
    allocate(reallocate_complex_2d(1:n,1:k), stat=ierr)
    if(ierr /= 0) stop "allocate error"
    if(.not. associated(p)) return
    nold = min(size(p,1), n)
    kold = min(size(p,2), k)
    reallocate_complex_2d(1:nold,1:kold) = p(1:nold,1:kold)
    deallocate(p)
    end function reallocate_complex_2d
!
 end module realloc
