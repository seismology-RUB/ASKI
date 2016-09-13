!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
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

! compile by:
! gfortran -O3 -o xgoertzel goertzel.f90

program goertzel

  implicit none

  ! time series
  character(len=400) :: seisfile
  real, dimension(:), allocatable :: s
  integer :: NSTEP,jt
  double precision :: DT
  real :: t_dummy

  ! frequency discretization
  integer, dimension(:), allocatable :: ASKI_jf
  integer :: ASKI_nf,jf
  double precision :: ASKI_df

  ! conventional DFT
  complex(kind=kind(1.d0)), dimension(:,:), allocatable :: efactors
  complex(kind=kind(1.d0)), dimension(:), allocatable :: spectrum_dble

  ! Goertzel DFT
  complex(kind=kind(1.d0)) :: spectral_value
  double precision, dimension(:), allocatable :: Wr
  double precision, dimension(:), pointer :: U0,U1,U2,U2_tmp
  double precision :: x

  ! Hanning taper
  logical :: apply_hanning_taper
  real :: hanning_taper_tail
  double precision :: wtaper
  integer :: ntaper_start
  double precision, dimension(:), allocatable :: taper_values

  integer :: ios,lu

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI

  integer, parameter :: number_of_repetitions = 50000
  integer :: ibench
  real :: cpu_start,cpu_end,duration_summation,duration_goertzel

  nullify(U0,U1,U2,U2_tmp)


  ! exemplary time series 
  ! (measured data from "cross borehole" example inversion, event SB04, station RC02, component CY
  ! as in GJI 2016 ASKI paper Figure 7, upper-right, and in Florian's dissertation Figure 5.4 (b) )

  seisfile = './data_SB04_OUTPUT_FILES--NETW.RC02.FXY.semd'
  DT = 2.114d-4
  NSTEP = 1184
  lu = 11
  open(unit=lu,file=seisfile,form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     write(*,*) "ERROR: could not open file '",trim(seisfile),"' to read on unit ",lu,"; raised iostat = ",ios
     stop
  end if
  allocate(s(1:NSTEP))
  do jt = 1,NSTEP
     read(lu,*,iostat=ios) t_dummy,s(jt)
     if(ios/=0) then
        write(*,*) "ERROR: could not read in ",jt,"'th line of seismogram file '",trim(seisfile),"'"
        stop
     end if
  end do ! jt
  write(*,*) "TIME SERIES:"
  write(*,*) "seismogram file: '",trim(seisfile),"'"
  write(*,*) "DT = ",DT
  write(*,*) "NSTEP = ",NSTEP
  write(*,*) ""


  ! Choose whether to apply a Hanning taper (and to which tail portion (between 0.0 and 1.0) it should be applied).
  ! Note that for the exemplary seismogram used in the following, a Hanning taper of 5 percent was appplied
  ! in the "cross borehole" inversion example (hard-coded in executable transformSpecfem3dCartesianMeasuredData
  ! of ASKI extension package SPECFEM3D_Cartesian_for_ASKI, used to produce file measured_data--data_SB04_RC02_CY)

  apply_hanning_taper = .true.

  if(apply_hanning_taper) then
     hanning_taper_tail = 0.05
     wtaper = DT*dble(NSTEP-1)*dble(hanning_taper_tail)
     ntaper_start = NSTEP - int(wtaper/DT) + 1 ! 2 <= ntaper_start <= NSTEP+1
     if(ntaper_start<=NSTEP) then
        allocate(taper_values(ntaper_start:NSTEP))
        do jt = ntaper_start,NSTEP
           taper_values(jt) = 0.5d0*(1.d0-dcos(PI*DT*dble(NSTEP-jt)/wtaper))
        end do ! jt
     else
        apply_hanning_taper = .false.
     end if
  else ! apply_hanning_taper
     ntaper_start = NSTEP+1
  end if ! apply_hanning_taper
  if(apply_hanning_taper) then
     write(*,*) "APPLYING HANNING TAPER ON TAILING ",NSTEP-ntaper_start+1," SAMPLES OF TIME SERIES (i.e. ",&
          hanning_taper_tail*100.," percent) BEFORE FOURIER TRANSFORM"
     write(*,*) ""
  else
     write(*,*) "APPLYING NOE HANNING TAPER BEFORE FOURIER TRANSFORM"
     write(*,*) ""
  end if


  ! frequency discretization
  ! (as in "cross borehole" example inversion)

  ASKI_nf = 21
  allocate(ASKI_jf(ASKI_nf))
  ASKI_jf = (/18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38/)
  ASKI_df = 2.d0
  !ASKI_df = 1/(NSTEP*DT)
  write(*,*) "FREQUENCY DISCRETIZATION:"
  write(*,*) "df = ",ASKI_df
  write(*,*) "nf = ",ASKI_nf
  write(*,*) "jf = ",ASKI_jf
  write(*,*) ""



  ! ++++ CONVENTIONAL DFT BY EXPLICIT SUM ++++

  write(*,*) "CONVENTIONAL DFT BY EXPLICIT SUMMATION:"
  

  ! defining efactors (including line element DT of integration) before time loop
  allocate(efactors(ASKI_nf,NSTEP))
  do jt = 1,NSTEP
     do jf = 1,ASKI_nf
        efactors(jf,jt) = exp(-(0.d0,1.d0)*TWO_PI*(ASKI_jf(jf)*ASKI_df)*(dble(jt-1)*DT)) * DT
     end do ! jf
  end do ! jt
  if(apply_hanning_taper) then
     do jt = ntaper_start,NSTEP
        efactors(:,jt) = efactors(:,jt) * taper_values(jt)
     end do ! jt
  end if

  ! time loop
  allocate(spectrum_dble(ASKI_nf))

  write(*,*) "will repeat the time loop ",number_of_repetitions," times for benchmarking (adapt in the code, if necessary)"
  call cpu_time(cpu_start)
  do ibench = 1,number_of_repetitions

     spectrum_dble = (0.d0,0.d0)
     do jt = 1,NSTEP
        do jf = 1,ASKI_nf
           spectrum_dble(jf) = spectrum_dble(jf) + s(jt)*efactors(jf,jt)
        end do ! jf
     end do ! jt

  end do ! ibench
  call cpu_time(cpu_end)
  duration_summation = cpu_end-cpu_start
  write(*,*) "total time [s] = ",duration_summation," ; i.e. per repition = ",(duration_summation)/real(number_of_repetitions)

  write(*,*) "resulting spectrum:"
  do jf = 1,ASKI_nf
     write(*,*) cmplx(spectrum_dble(jf))
  end do ! jf
  write(*,*) ""



  ! ! ++++ Goertzel's ALGORITHM ++++
  ! ! ++++ original recursion from N down to 1 (not suitabel for on-the-fly Fourier transform) ++++

  ! write(*,*) "DFT BY ORIGINAL Goertzel's ALGORITHM:"

  ! ! defining constants before time loop:
  ! allocate(Wr(ASKI_nf),U0(ASKI_nf),U1(ASKI_nf),U2(ASKI_nf))
  ! do jf = 1,ASKI_nf
  !    x = -TWO_PI*ASKI_jf(jf)*ASKI_df*DT
  !    Wr(jf) = 2.d0 * cos(x)
  ! end do ! jf
  ! U1(:) = 0.d0
  ! U2(:) = 0.d0

  ! ! time loop
  ! do jt = NSTEP,2,-1
     
  !    if(apply_hanning_taper .and. jt >= ntaper_start) then
  !       do jf = 1,ASKI_nf
  !          U0(jf) = s(jt)*taper_values(jt) + Wr(jf)*U1(jf) - U2(jf)
  !       end do ! jf
  !    else
  !       do jf = 1,ASKI_nf
  !          U0(jf) = s(jt) + Wr(jf)*U1(jf) - U2(jf)
  !       end do ! jf
  !    end if
  !    ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
  !    U2_tmp => U2
  !    U2 => U1
  !    U1 => U0
  !    U0 => U2_tmp ! use the memory of U2 (pointed to by U2_tmp) for setting the values U0 in next time step

  ! end do ! jt

  ! write(*,*) "resulting spectrum:"
  ! do jf = 1,ASKI_nf
  !    x = -TWO_PI*ASKI_jf(jf)*ASKI_df*DT
  !    if(apply_hanning_taper .and. ntaper_start <= 1) then
  !       write(*,*) cmplx( DT*( s(1)*taper_values(1)+U1(jf)*0.5d0*Wr(jf)-U2(jf) ) , DT*sin(x)*U1(jf) )
  !    else
  !       write(*,*) cmplx( DT*( s(1)+U1(jf)*0.5d0*Wr(jf)-U2(jf) ) , DT*sin(x)*U1(jf) )
  !    end if
  ! end do ! jf
  ! write(*,*) ""



  ! ++++ Goertzel's ALGORITHM 
  ! ++++ reverse time by accounting for reversal property of Fourier transform and correcting phase ++++

  write(*,*) "DFT BY Goertzel ALGORITHM correcting for reverse time series (by Fourier transform rules):"

  if(.not.allocated(Wr)) allocate(Wr(ASKI_nf))
  if(.not.associated(U0)) allocate(U0(ASKI_nf))
  if(.not.associated(U1)) allocate(U1(ASKI_nf))
  if(.not.associated(U2)) allocate(U2(ASKI_nf))
  ! compute at NEGATIVE frequencies (accounting for time-reversal property, i.e. chosing -x here
  ! compared to (commented) original Goertzel algorithm above):
  do jf = 1,ASKI_nf
     x = TWO_PI*ASKI_jf(jf)*ASKI_df*DT
     Wr(jf) = 2.d0 * cos(x)
  end do ! jf

  ! time loop
  write(*,*) "will repeat the time loop ",number_of_repetitions," times for benchmarking (adapt in the code, if necessary)"
  call cpu_time(cpu_start)
  do ibench = 1,number_of_repetitions

     U1(:) = 0.d0
     U2(:) = 0.d0
     do jt = 1,NSTEP-1
     
        if(apply_hanning_taper .and. jt >= ntaper_start) then
           do jf = 1,ASKI_nf
              U0(jf) = s(jt)*taper_values(jt) + Wr(jf)*U1(jf) - U2(jf)
           end do ! jf
        else
           do jf = 1,ASKI_nf
              U0(jf) = s(jt) + Wr(jf)*U1(jf) - U2(jf)
           end do ! jf
        end if
        !   THE FOLLOWING ALTERNATIVE CODE LINE IS  E X T R E M E L Y  SLOW COMPARED WITH THE ABOVE FREQUENCY LOOPS
        !   (especially if compiler flags -O1 , -O2 or -O3 are used)
        ! U0 = (Wr*U1-U2) + s(jt)  or
        ! U0 = (Wr*U1-U2) + s(jt)*taper_values(jt)

        ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
        U2_tmp => U2
        U2 => U1
        U1 => U0
        U0 => U2_tmp ! use the memory of U2 (pointed to by U2_tmp) for setting the values U0 in next time step

     end do ! jt

  end do ! ibench
  call cpu_time(cpu_end)
  duration_goertzel = cpu_end-cpu_start
  write(*,*) "total time [s] = ",duration_goertzel," ; i.e. per repition = ",(duration_goertzel)/real(number_of_repetitions)
  write(*,*) "Goertzel's algorithm was ",duration_summation/duration_goertzel," times faster than the explicit summation"

  write(*,*) "resulting spectrum:"
  do jf = 1,ASKI_nf
     x = TWO_PI*ASKI_jf(jf)*ASKI_df*DT
     if(apply_hanning_taper) then
        spectral_value = dcmplx( DT*( s(NSTEP)*taper_values(NSTEP)+0.5d0*Wr(jf)*U1(jf)-U2(jf) ) , DT*sin(x)*U1(jf) )
     else
        spectral_value = dcmplx( DT*( s(NSTEP)+0.5d0*Wr(jf)*U1(jf)-U2(jf) ) , DT*sin(x)*U1(jf) )
     end if
     spectral_value = spectral_value*exp(-(0.d0,1.d0)*TWO_PI*(ASKI_jf(jf)*ASKI_df)*(NSTEP-1)*DT) ! shift the phase back to the right by record length T = (NSTEP-1)*DT
     write(*,*) cmplx(spectral_value)
  end do ! jf
  write(*,*) ""


  ! clean up
  if(allocated(s)) deallocate(s)
  if(allocated(taper_values)) deallocate(taper_values)
  if(allocated(ASKI_jf)) deallocate(ASKI_jf)
  if(allocated(efactors)) deallocate(efactors)
  if(allocated(spectrum_dble)) deallocate(spectrum_dble)
  if(allocated(Wr)) deallocate(Wr)
  if(associated(U0)) deallocate(U0)
  if(associated(U1)) deallocate(U1)
  if(associated(U2)) deallocate(U2)
end program goertzel
