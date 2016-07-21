!----------------------------------------------------------------------------
!   Copyright 2015 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.0.
!
!   ASKI version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief module for discrete Fourier transform
!!
!! \details In particular developed for transform from time to frequency domain
!!  in case of a sparse frequency spectrum (few frequency samples, in general
!!  not equidistant). A "direct" discrete Fourier transform is applied, explicitely
!!  computing coefficients efactors(omega=2*pi*f,t) = exp ( +- i * omega * t ) with
!!  which a Rieman-sum type of integration is conducted, (hence, very basic after all)
!!  FOR THE FUTURE it may also be helpful, if some fast Fourier transform algorithms
!!  are incorporated into the module, for faster transformation in case of generating
!!  time_waveform_kernel objects.
!!
!! \author Florian Schumacher
!! \date August 2013
!
module discreteFourierTransform
!
  use errorMessage
  use mathConstants
!
  implicit none
!
  private :: computeEfactorsDFT
!
  interface initiateForwardDFT
     module procedure initiateRealFreqForwardDFT
     module procedure initiateComplexFreqForwardDFT
  end interface initiateForwardDFT
  interface transformForwardDFT
     module procedure transformTracesForwardDFT
     module procedure transformOneTraceForwardDFT
  end interface transformForwardDFT
  interface transformInverseDFT
     module procedure transformSpectraInverseDFT
     module procedure transformSpectraOnTheFlyInverseDFT
     module procedure transformOneSpectrumInverseDFT
     !module procedure transformOneSpectrumOnTheFlyInverseDFT
  end interface transformInverseDFT
  interface dealloc
     module procedure deallocateDFT
  end interface dealloc
!
!< collection of complex exponential coefficients (and meta infor about them) with which time series (spectra) are discretely "convolved"
  type discrete_fourier_transform
     private
     logical :: is_defined = .false. !< logical value whether the object has been initiated yet or not
     !
     ! spectral domain specifications
     ! df,nf1,nf2 are only used for inverse transform!
     double precision :: df !< frequency step of spectra
     integer :: nf1 !< first frequency index defining starting frequency f1 = nf1*df of spectra (positive, and nf2 > nf1)
     integer :: nf2 !< second frequency index defining end frequency f2 = nf2*df of spectra (positive, and nf2 > nf1)
     integer :: nf !< number of frequencies (lengt of f, number of rows of efactors), nf = nf2-nf1+1 for inverse transform
     double complex, dimension(:), pointer :: f => null() !< frequencies for which (rows of) efactors are computed (order related to order of rows of efactors)
     !
     ! time domain specifications
     ! dt,nt1,nt2 are only used for forward transform!
     double precision :: dt !< time step of time series
     integer :: nt1 !< first time index defining start time of the time series t1 = nt1*dt (can be negative, but nt2 > nt1)
     integer :: nt2 !< second time index defining end time of the time series t2 = nt2*dt (can be negative, but nt2 > nt1)
     integer :: nt !< total number of time samples (length of t, columns of efactors)
     double precision, dimension(:), pointer :: t => null() !< times for which (columns of) efactors are computed (order related to order of columns of efactors)
     !
     ! complex exponential coefficients used to compute fourie transforms
     logical :: forward !< logical indicating whether efactors contains coefficients for forward or inverse Fourier transform
     double complex, dimension(:,:), pointer :: efactors => null() !< complex exponential coefficients for forward or inverse Fourier transform
  end type discrete_fourier_transform
!
contains
!
!------------------------------------------------------------------------
!> \brief initiate complex exponential coefficients for direct forward Fourier transform, real frequencies
!! \details simply calls subroutine initiateComplexFreqForwardDFT with argument cmplx(f)
!! \param this DFT object
!! \param dt time step of the time series which will be transformed
!! \param nt1 first time index defining starting time t1 = nt1*dt of the time series which will be transformed
!! \param nt2 second time index defining end time t2 = nt2*dt of the time series which will be transformed
!! \param f vector of real frequency values, at which the spectrum of the Fourier transform should be evaluated
!! \param errmsg error message
!! \param hanning_taper value between 0. and 1., defining a tail portion of the timeseries which are tapered by a cos hanning taper
!
  subroutine initiateRealFreqForwardDFT(this,dt,nt1,nt2,f,errmsg,hanning_taper)
    type (discrete_fourier_transform) :: this
    real :: dt
    integer :: nt1,nt2
    real, dimension(:) :: f
    type (error_message) :: errmsg
    real, optional :: hanning_taper
    ! local
    character(len=26) :: myname = 'initiateRealFreqForwardDFT'
!
    call addTrace(errmsg,myname)
!
    call initiateComplexFreqForwardDFT(this,dt,nt1,nt2,cmplx(f),errmsg,hanning_taper)
  end subroutine initiateRealFreqForwardDFT
!------------------------------------------------------------------------
!> \brief initiate complex exponential coefficients for direct forward Fourier transform
!! \details Explicitely compute coefficients of the form exp(-i*omega*t)*dt for numeric
!!  integration of  time_series(t)*exp(-i*omega*t) over all t. These coefficients are computed
!!  for some given (few) frequencies and are used in routines in interface transformForwardDFT
!! \param this DFT object
!! \param dt time step of the time series which will be transformed
!! \param nt1 first time index defining starting time t1 = nt1*dt of the time series which will be transformed
!! \param nt2 second time index defining end time t2 = nt2*dt of the time series which will be transformed
!! \param f vector of real frequency values, at which the spectrum of the Fourier transform should be evaluated
!! \param errmsg error message
!! \param hanning_taper value between 0. and 1., defining a tail portion of the timeseries which are tapered by a cos hanning taper
!
  subroutine initiateComplexFreqForwardDFT(this,dt,nt1,nt2,f,errmsg,hanning_taper)
    type (discrete_fourier_transform) :: this
    real :: dt
    integer :: nt1,nt2
    complex, dimension(:) :: f
    type (error_message) :: errmsg
    real, optional :: hanning_taper
    ! local
    character(len=29) :: myname = 'initiateComplexFreqForwardDFT'
    character(len=400) :: errstr
    logical :: return_after_check,apply_taper
    double precision :: wtaper
    integer :: ntaper,it
!
    call addTrace(errmsg,myname)
!
    ! first check incoming values
    return_after_check = .false.
!
    if(dt <= 0) then
       write(errstr,*) "incoming timestep dt = ",dt," must be positive"
       call add(errmsg,2,errstr,myname)
       return_after_check = .true.
    end if
    if(nt2 <= nt1) then
       write(errstr,*) "incoming time indices nt1,nt2 = ",nt1,nt2," invalid, must be nt1 < nt2"
       call add(errmsg,2,errstr,myname)
       return_after_check = .true.
    end if
!
    if(size(f) <= 0) then
       call add(errmsg,2,"there are no incoming frequencies",myname)
       return_after_check = .true.
    end if
!
    if(present(hanning_taper)) then
       if(hanning_taper < 0. .or. hanning_taper > 1.) then
          write(errstr,*) "incoming percentage ",hanning_taper," of tapering the tail of time "//&
               "series must be between 0. and 1."
          call add(errmsg,2,errstr,myname)
          return_after_check = .true.          
       end if
       if(hanning_taper == 0.) then
          apply_taper = .false.
       else
          apply_taper = .true.
       end if
    else
       apply_taper = .false.
    end if
!
    if(return_after_check) return
!
    if(this%is_defined) then
       call add(errmsg,1,"DFT object already initiated, deallocating it now before defining new one",myname)
       call deallocateDFT(this)
    end if
!
    ! now start defining object
!
    this%nt1 = nt1; this%nt2 = nt2; this%nt = nt2-nt1+1
    this%dt = dble(dt)
    allocate(this%t(this%nt))
    this%t = (/ (dble(it)*this%dt, it=nt1,nt2) /)
!
    this%nf = size(f)
    allocate(this%f(this%nf))
    this%f = dcmplx(f)
!
    this%forward = .true.
    allocate(this%efactors(this%nf,this%nt))
    call computeEfactorsDFT(this)
!
    ! multiply by volume element dt for integration over time
    this%efactors = this%efactors*this%dt
!
    ! apply taper
    if(apply_taper) then
       ! checked above:
       ! - hanning_taper is present
       ! - hanning_taper > 0.
       ! - hanning_taper <= 1.
       wtaper = this%dt*dble(this%nt-1)*dble(hanning_taper) ! hence, 0. < wtaper <= t2-t1 , whereby ti = dt*nti (start/end time of time series)
       ntaper = wtaper/this%dt ! 0 <= ntaper <= nt-1

       if(ntaper>0) then
          do it = this%nt-ntaper+1,this%nt
             this%efactors(:,it) = this%efactors(:,it)* &
                  (0.5d0*(1.d0-dcos(mc_pid*this%dt*dble(this%nt-it)/wtaper)))
          end do ! it
       end if ! ntaper>0
    end if ! apply_taper
!
    this%is_defined = .true.
!    
  end subroutine initiateComplexFreqForwardDFT
!------------------------------------------------------------------------
!> \brief initiate complex exponential coefficients for direct inverse Fourier transform
!! \details Explicitely compute coefficients of the form exp(i*omega*t)*dt for numeric
!!  integration of  time_series(t)*exp(i*omega*t) over all t. These coefficients are computed
!!  for equidistant frequencies in a sufficiently large frequency window, which must have been properly
!!  chosen in the forward simulations! The coefficients are computed for an equidistant time window. They
!!  are used in routines in interface transformInverseDFT
!! \param this DFT object
!! \param df frequency step of the frequency window in which the spectra are given (used to compute this%f)
!! \param nf1 first frequency index defining starting frequency f1 = nf1*df of the spectrum which will be transformed
!! \param nf2 second frequency index defining end frequency f2 = nf2*df of the spectrum which will be transformed
!! \param t vector of time values, at which the time series of the inverse Fourier transform should be evaluated
!! \param errmsg error message
!
  subroutine initiateInverseDFT(this,df,nf1,nf2,t,errmsg)
    type (discrete_fourier_transform) :: this
    real :: df
    integer :: nf1,nf2
    real, dimension(:) :: t
    type (error_message) :: errmsg
    ! local
    character(len=18) :: myname = 'initiateInverseDFT'
    character(len=400) :: errstr
    logical :: return_after_check
    integer :: ifreq
!
    call addTrace(errmsg,myname)
!
    ! first check incoming values
    return_after_check = .false.
!
    if(df <= 0) then
       write(errstr,*) "incoming frequency step df = ",df," must be positive"
       call add(errmsg,2,errstr,myname)
       return_after_check = .true.
    end if
    if(nf1 < 0 .or. nf2 < 0 .or. nf2 <= nf1) then
       write(errstr,*) "incoming frequency indices nf1,nf2 = ",nf1,nf2," invalid, must both not be negative and nf1 < nf2"
       call add(errmsg,2,errstr,myname)
       return_after_check = .true.
    end if
!
    if(size(t) <= 0) then
       call add(errmsg,2,"there are no incoming time samples",myname)
       return_after_check = .true.
    end if
!
    if(return_after_check) return
!
    if(this%is_defined) then
       call add(errmsg,1,"DFT object already initiated, deallocating it now before defining new one",myname)
       call deallocateDFT(this)
    end if
!
    ! now start defining object
!
    this%nf1 = nf1; this%nf2 = nf2; this%nf = nf2-nf1+1
    this%df = dble(df)
    allocate(this%f(this%nf))
    this%f = (/ (dble(ifreq)*this%df, ifreq=nf1,nf2) /)
!
    this%nt = size(t)
    allocate(this%t(this%nt))
    this%t = dble(t)
!
    this%forward = .false.
    allocate(this%efactors(this%nf,this%nt))
    call computeEfactorsDFT(this)
!
    ! multiply by volume element df for integration over frequency
    this%efactors = this%efactors*dble(df)
!
    this%is_defined = .true.
!    
  end subroutine initiateInverseDFT
!------------------------------------------------------------------------
!> \brief compute coefficients of the form exp(+-i*omega*t) to be used in forward and inverse DFT
!! \param this DFT object
!
  subroutine computeEfactorsDFT(this)
    type (discrete_fourier_transform) :: this
    ! local
    double complex :: plus_minus_i_2pi
    integer :: it,ifreq
    if(this%forward) then
       plus_minus_i_2pi = -mc_cid*mc_two_pid
    else
       plus_minus_i_2pi = mc_cid*mc_two_pid
    end if
!
! definintion of pure efactors: no tapering, no multiplication by dt or df
! to be used in both, forward and inverse transformation
    do it = 1,this%nt
       do ifreq = 1,this%nf
          this%efactors(ifreq,it) = &
               cdexp(plus_minus_i_2pi*this%f(ifreq)*this%t(it))
       end do ! ifreq
    end do ! it
!
  end subroutine computeEfactorsDFT
!------------------------------------------------------------------------
!> \brief transform more than one time series to frequency domain
!! \details Using the complex exponential coefficients defined in routine initateForwardDFT
!!  each incoming trace is transformed by a matrix vector multiplication efactors*trace
!!  which essentially computes the sum of the time samples weighted by respective exponential
!!  coefficients, which also include (simple constant) integeration weights, such that the 
!!  actual Fourier integral is computed.
!! \param this DFT object
!! \param traces incoming (nt,ntrace)-array of traces
!! \param spectra (nf,ntrace)-array containing allocated space for the spectra of the transformed traces
!! \param errmsg error message
!
  subroutine transformTracesForwardDFT(this,traces,spectra,errmsg)
    type (discrete_fourier_transform) :: this
    real, dimension(:,:) :: traces
    complex, dimension(:,:) :: spectra
    type (error_message) :: errmsg
    ! local
    character(len=25) :: myname = 'transformTracesForwardDFT'
    character(len=400) :: errstr
    integer :: ntrace,nt_in,nf_in
!!$integer :: itrace,it,ifreq
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"object not initiated yet, call initiateForwardDFT first",myname)
       return
    end if
    if(.not.this%forward) then
       call add(errmsg,2,"object was not initiated for forward but inverse Fourier transform",myname)
       return
    end if
!
    nt_in = size(traces,1)
    if(nt_in /= this%nt) then
       write(errstr,*) "incoming number of time samples of traces ",nt_in,&
            " does not equal number of time samples ",this%nt," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    ntrace = size(traces,2)
    if(ntrace <= 0) then
       write(errstr,*) "incoming number of traces ",ntrace," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nf_in = size(spectra,1)
    if(nf_in /= this%nf) then
       write(errstr,*) "incoming number of frequency samples of spectra ",nf_in,&
            " does not equal number of frequencies ",this%nf," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(spectra,2) /= ntrace) then
       write(errstr,*) "incoming number of traces ",ntrace," does not equal incoming number of spectra ",&
            size(spectra,2)
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! the type and kind of the result of matmul follow the usual type and kind promotion rules, as for the * operator,
    ! says http://gcc.gnu.org/onlinedocs/gfortran/MATMUL.html
    spectra = matmul(this%efactors,traces)
!!$    do itrace = 1,ntrace
!!$       do ifreq = 1,this%nf
!!$          spectra(ifreq,itrace) = (0.,0.)
!!$          do it = 1,this%nt
!!$             spectra(ifreq,itrace) = spectra(ifreq,itrace) + cmplx(this%efactors(ifreq,it))*traces(it,itrace)
!!$          end do ! it
!!$       end do ! ifreq
!!$    end do ! itrace
  end subroutine transformTracesForwardDFT
!------------------------------------------------------------------------
!> \brief transform exactly one time series to frequency domain
!! \details Using the complex exponential coefficients defined in routine initateForwardDFT
!!  the incoming trace is transformed by a matrix vector multiplication efactors*trace
!!  which essentially computes the sum of the time samples weighted by respective exponential
!!  coefficients, which also include (simple constant) integeration weights, such that the 
!!  actual Fourier integral is computed.
!! \param this DFT object
!! \param trace incoming (nt)-array containing time series
!! \param spectrum (nf)-array containing allocated specfe for the spectrum of the transformed trace
!! \param errmsg error message
!
  subroutine transformOneTraceForwardDFT(this,trace,spectrum,errmsg)
    type (discrete_fourier_transform) :: this
    real, dimension(:) :: trace
    complex, dimension(:) :: spectrum
    type (error_message) :: errmsg
    ! local
    character(len=27) :: myname = 'transformOneTraceForwardDFT'
    character(len=400) :: errstr
    integer :: nt_in,nf_in
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"object not initiated yet, call initiateForwardDFT first",myname)
       return
    end if
    if(.not.this%forward) then
       call add(errmsg,2,"object was not initiated for forward but inverse Fourier transform",myname)
       return
    end if
!
    nt_in = size(trace)
    if(nt_in /= this%nt) then
       write(errstr,*) "incoming number of time samples of trace ",nt_in,&
            " does not equal number of time samples ",this%nt," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nf_in = size(spectrum)
    if(nf_in /= this%nf) then
       write(errstr,*) "incoming number of frequency samples of spectrum ",nf_in,&
            " does not equal number of frequencies ",this%nf," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! the type and kind of the result of matmul follow the usual type and kind promotion rules, as for the * operator,
    ! says http://gcc.gnu.org/onlinedocs/gfortran/MATMUL.html
    spectrum = matmul(this%efactors,trace)
  end subroutine transformOneTraceForwardDFT
!------------------------------------------------------------------------
!> \brief transform more than one frequency spectrum to time domain
!! \details Using the complex exponential coefficients defined in routine initiateInverseDFT
!!  each incoming spectrum is transformed by a matrix vector multiplication efactors^T*spectrum
!!  which essentially computes the sum of the spectral values weighted by respective exponential
!!  coefficients, which also include (simple constant) integeration weights, such that the 
!!  actual inverse Fourier integral is computed.
!!  It is assumed here, that the original time domain signal of the spectra was real! This assumption
!!  enables to treat the negative frequencies by the complex conjugate values of the corresponding positive
!!  frequency!
!! \param this DFT object
!! \param spectra incoming (nf,ntrace)-array of spectra
!! \param traces incoming (nt,ntrace)-array containing allocated space for the transformed traces
!! \param errmsg error message
!
  subroutine transformSpectraInverseDFT(this,spectra,traces,errmsg)
    type (discrete_fourier_transform) :: this
    real, dimension(:,:) :: traces
    complex, dimension(:,:) :: spectra
    type (error_message) :: errmsg
    ! local
    character(len=25) :: myname = 'transformTracesForwardDFT'
    character(len=400) :: errstr
    integer :: nspectrum,nt_in,nf_in
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"object not initiated yet, call initiateInverseDFT first",myname)
       return
    end if
    if(this%forward) then
       call add(errmsg,2,"object was not initiated for inverse but forward Fourier transform",myname)
       return
    end if
!
    nf_in = size(spectra,1)
    if(nf_in /= this%nf) then
       write(errstr,*) "incoming number of frequency samples of spectra ",nf_in,&
            " does not equal number of frequencies ",this%nf," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    nspectrum = size(spectra,2)
    if(nspectrum <= 0) then
       write(errstr,*) "incoming number of spectra ",nspectrum," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nt_in = size(traces,1)
    if(nt_in /= this%nt) then
       write(errstr,*) "incoming number of time samples of traces ",nt_in,&
            " does not equal number of time samples ",this%nt," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(traces,2) /= nspectrum) then
       write(errstr,*) "incoming number of spectra ",nspectrum," does not equal incoming number of traces ",&
            size(traces,2)
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! the type and kind of the result of matmul follow the usual type and kind promotion rules, as for the * operator,
    ! says http://gcc.gnu.org/onlinedocs/gfortran/MATMUL.html
    traces = real( matmul(transpose(this%efactors),spectra) + matmul(transpose(conjg(this%efactors)),conjg(spectra)) )
  end subroutine transformSpectraInverseDFT
!------------------------------------------------------------------------
!> \brief compute one summand of on-the-fly transformation of more than one frequency spectrum to time domain
!! \details Using the complex exponential coefficients defined in routine initiateInverseDFT, for
!!  every incoming spectral value (at frequency (index) jf) the jf'th summand of the discrete fourier 
!!  transformation is computed by the multiplication of column-vector efactors(jf,:)^T with row-vector 
!!  spectral_values, yielding a matrix. These matrices sum up (over all frequencies) to matrix traces as
!!  returned by routine transformSpectraInverseDFT.
!!  It is assumed here, that the original time domain signal of the spectrum was real! This assumption
!!  enables to treat the negative frequencies by the complex conjugate values of the corresponding positive
!!  frequency!
!! \param this DFT object
!! \param spectra_one_freq incoming (nspectra)-array of spectral values of all spectra at frequency (index) jf
!! \param jf frequency index at which the incoming spectral values spectra_one_freq are given (nf1<=jf<=nf2)
!! \param traces (nt,nspectra)-array containing allocated space for the (frequency summand of the) transformed traces
!! \param errmsg error message
!! \param add optional logical indicating whether the jf'th summand of the DFT shoud be added to incoming values contained in traces or not
!
  subroutine transformSpectraOnTheFlyInverseDFT(this,spectra_one_freq,jf,traces,errmsg,add2traces)
    type (discrete_fourier_transform) :: this
    complex, dimension(:) :: spectra_one_freq
    integer :: jf
    real, dimension(:,:) :: traces
    type (error_message) :: errmsg
    ! optional
    logical, optional :: add2traces
    ! local
    character(len=34) :: myname = 'transformSpectraOnTheFlyInverseDFT'
    character(len=400) :: errstr
    integer :: nspectrum,nt_in,ispectrum,it,jf_shifted
    logical :: add_to_traces
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"object not initiated yet, call initiateInverseDFT first",myname)
       return
    end if
    if(this%forward) then
       call add(errmsg,2,"object was not initiated for inverse but forward Fourier transform",myname)
       return
    end if
!
    if(jf<this%nf1 .or. jf>this%nf2) then
       write(errstr,*) "incoming frequency index ",jf," is out of range of indices nf1,nf2 = ",&
            this%nf1,this%nf2," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nspectrum = size(spectra_one_freq)
    if(nspectrum <= 0) then
       write(errstr,*) "incoming number of spectral frequency samples ",nspectrum," at frequency index ",&
            jf," must be positive"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nt_in = size(traces,1)
    if(nt_in /= this%nt) then
       write(errstr,*) "incoming number of time samples of traces ",nt_in,&
            " does not equal number of time samples ",this%nt," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
    if(size(traces,2) /= nspectrum) then
       write(errstr,*) "incoming number of spectra ",nspectrum," does not equal incoming number of traces ",&
            size(traces,2)
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    add_to_traces = .false.
    if(present(add2traces)) add_to_traces = add2traces
!
    ! shift back index jf to local index in first dimension of array this%efactors (which is indexed from 1 to nf = nf2-nf1+1)
    jf_shifted = jf-this%nf1+1
!
    if(add_to_traces) then
       do ispectrum = 1,nspectrum
          do it = 1,this%nt
             traces(it,ispectrum) = traces(it,ispectrum) + real( this%efactors(jf_shifted,it)*spectra_one_freq(ispectrum) + &
                  conjg(this%efactors(jf_shifted,it))*conjg(spectra_one_freq(ispectrum)) )
          end do ! it
       end do ! ispectrum

    else ! add_to_traces

       do ispectrum = 1,nspectrum
          do it = 1,this%nt
             traces(it,ispectrum) = real( this%efactors(jf_shifted,it)*spectra_one_freq(ispectrum) + &
                  conjg(this%efactors(jf_shifted,it))*conjg(spectra_one_freq(ispectrum)) )
          end do ! it
       end do ! ispectrum
    end if ! add_to_traces
  end subroutine transformSpectraOnTheFlyInverseDFT
!------------------------------------------------------------------------
!> \brief transform exactly one frequency spectrum to time domain
!! \details Using the complex exponential coefficients defined in routine initiateInverseDFT
!!  the incoming spectrum is transformed by a matrix vector multiplication efactors^T*spectrum
!!  which essentially computes the sum of the spectral values weighted by respective exponential
!!  coefficients, which also include (simple constant) integeration weights, such that the 
!!  actual inverse Fourier integral is computed.
!!  It is assumed here, that the original time domain signal of the spectrum was real! This assumption
!!  enables to treat the negative frequencies by the complex conjugate values of the corresponding positive
!!  frequency!
!! \param this DFT object
!! \param spectrum incoming (nf)-array of spectrum
!! \param trace incoming (nt)-array containing allocated space for the transformed trace
!! \param errmsg error message
!
  subroutine transformOneSpectrumInverseDFT(this,spectrum,trace,errmsg)
    type (discrete_fourier_transform) :: this
    real, dimension(:) :: trace
    complex, dimension(:) :: spectrum
    type (error_message) :: errmsg
    ! local
    character(len=30) :: myname = 'transformOneSpectrumInverseDFT'
    character(len=400) :: errstr
    integer :: nt_in,nf_in
!
    call addTrace(errmsg,myname)
!
    if(.not.this%is_defined) then
       call add(errmsg,2,"object not initiated yet, call initiateInverseDFT first",myname)
       return
    end if
    if(this%forward) then
       call add(errmsg,2,"object was not initiated for inverse but forward Fourier transform",myname)
       return
    end if
!
    nf_in = size(spectrum)
    if(nf_in /= this%nf) then
       write(errstr,*) "incoming number of frequency samples of spectrum ",nf_in,&
            " does not equal number of frequencies ",this%nf," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    nt_in = size(trace)
    if(nt_in /= this%nt) then
       write(errstr,*) "incoming number of time samples of trace ",nt_in,&
            " does not equal number of time samples ",this%nt," for which DFT has been initiated"
       call add(errmsg,2,errstr,myname)
       return
    end if
!
    ! the type and kind of the result of matmul follow the usual type and kind promotion rules, as for the * operator,
    ! says http://gcc.gnu.org/onlinedocs/gfortran/MATMUL.html
    trace = real( matmul(transpose(this%efactors),spectrum) + matmul(transpose(conjg(this%efactors)),conjg(spectrum)) )
  end subroutine transformOneSpectrumInverseDFT
!------------------------------------------------------------------------
!> \brief deallocate object of type discrete_fourier_transform
!! \param this DFT object
!
  subroutine deallocateDFT(this)
    type (discrete_fourier_transform) :: this
    if(associated(this%f)) deallocate(this%f)
    if(associated(this%efactors)) deallocate(this%efactors)
    this%is_defined = .false.
  end subroutine deallocateDFT
!
end module discreteFourierTransform
