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
!--------------------------------------------------
!  Module with mathematical constants
!--------------------------------------------------
module mathConstants
    real, parameter :: pi = 3.141592653589793                          ! pi
    complex, parameter :: mc_ci = (0., 1.)                             ! complex i
    real, parameter :: mc_pi = 3.141592653589793                       ! pi single
    real, parameter :: mc_two_pi =  2.0 * mc_pi                        ! 2*pi single
    real, parameter :: mc_deg2rad = 3.141592653589793/180.             ! degree to radian
!
!  doubles
!
    double precision, parameter :: mc_pid = 3.141592653589793          ! pi double precision
    double precision, parameter :: mc_two_pid =  2.d0 * mc_pid         ! 2*pi double precision
    double precision, parameter :: mc_ed  = 2.718281828459045          ! e double precision
    double precision, parameter :: mc_deg2radd = 3.141592653589793d0/180.d0   ! degree to radian
    double complex, parameter :: mc_cid = (0.d0, 1.d0)                 ! double complex i
!
end module mathConstants 
