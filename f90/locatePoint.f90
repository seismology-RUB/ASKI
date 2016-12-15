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
!----------------------------------------------------------------------------
!--------------------------------------------------------------
!  module with locate functions
!--------------------------------------------------------------------
!  For given value of xx and an ordered sequence x(i), i=1...n
!  find the index j such that xx lies between x(j) and x(j+1).
!  If xx sits exactly on one of the x(i), j is chosen such that
!  xx sits on the greater one of x(j) and x(j+1). This means:
!  x(i) ascending:  x(j) < xx <= x(j+1)
!  x(i) descending  x(j) >= xx > x(j+1)
!-------------------------------------------------------
 module locatePoint
    implicit none
    interface locate
        module procedure locate_real
        module procedure locate_double_precision
        module procedure locate_int
    end interface locate
 contains
!--------------------------------------------------------------
!  real array
!
    integer function locate_real(xx,n,x)
    integer :: n
    real :: xx,x(n)
    integer :: lidx,ridx,midx
!
!  start with lower and upper index bounds 0 and n+1
!
    lidx = 0
    ridx = n+1
    do while (ridx-lidx .gt. 1)
        midx = (lidx+ridx)/2
!
!  x(i) is an ascending sequence
!
        if (x(n) .gt. x(1)) then
            if (xx .gt. x(midx)) then
                lidx = midx
            else
                ridx = midx
            endif
!
!  x(i) is a descending sequence
!
        else
            if (xx .gt. x(midx)) then
                ridx = midx
            else
                lidx = midx
            endif
        endif
    enddo
    locate_real = lidx
    end function locate_real
!----------------------------------------------------------------
!  double precision array
!
    integer function locate_double_precision(xx,n,x)
    integer :: n
    double precision :: xx,x(n)
    integer :: lidx,ridx,midx
!
!  start with lower and upper index bounds 0 and n+1
!
    lidx = 0
    ridx = n+1
    do while (ridx-lidx .gt. 1)
        midx = (lidx+ridx)/2
!
!  x(i) is an ascending sequence
!
        if (x(n) .gt. x(1)) then
            if (xx .gt. x(midx)) then
                lidx = midx
            else
                ridx = midx
            endif
!
!  x(i) is a descending sequence
!
        else
            if (xx .gt. x(midx)) then
                ridx = midx
            else
                lidx = midx
            endif
        endif
    enddo
    locate_double_precision = lidx
    end function locate_double_precision
!--------------------------------------------------------------
!  integer array
!
    integer function locate_int(xx,n,x)
    integer :: n
    integer :: xx,x(n)
    integer :: lidx,ridx,midx
!
!  start with lower and upper index bounds 0 and n+1
!
    lidx = 0
    ridx = n+1
    do while (ridx-lidx .gt. 1)
        midx = (lidx+ridx)/2
!
!  x(i) is an ascending sequence
!
        if (x(n) .gt. x(1)) then
            if (xx .gt. x(midx)) then
                lidx = midx
            else
                ridx = midx
            endif
!
!  x(i) is a descending sequence
!
        else
            if (xx .gt. x(midx)) then
                ridx = midx
            else
                lidx = midx
            endif
        endif
    enddo
    locate_int = lidx
    end function locate_int
!
 end module locatePoint
