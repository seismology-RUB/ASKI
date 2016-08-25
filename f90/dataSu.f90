!----------------------------------------------------------------------------
!   Copyright 2016 Nils Müller (Ruhr-Universitaet Bochum, Germany)
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
!-------------------------------------------------------------
!   ____    __   ____   __    ___  __  __ 
!  (  _ \  /__\ (_  _) /__\  / __)(  )(  )
!   )(_) )/(__)\  )(  /(__)\ \__ \ )(__)( 
!  (____/(__)(__)(__)(__)(__)(___/(______)
!
!  - reads SeismicUnix files
!  - writes SeismicUnix files
!
!  subroutines:
!
!  allocDataSu(this,n)				-	allocates memory for dataSu object [this] to store [n] traces.
!  addTraceDataSu(this,ni,trace,sampint,...	-	adds trace to slot [ni] in dataSu object [this]. Array containing
!							trace values [trace] are passed by reference. Sampling Interval [sampint]
!							has to be given. All other items are optional.
!  delTraceDataSu(this,ni)			-	deletes trace in slot [ni] in dataSu object [this].
!  writeDataSu(this,lu,filename,endian)		-	writes data stored in dataSu object [this] to file [filename]. The io-unit
!							is given by [lu]. Endianness of the file is given by the optional argument
!							[endian] (not present or 0 = system standard, 1 = convert to little_endian,
!							2 = convert to big_endian).
!  readDataSu(this,lu,filename,endian)		-	reads data stored in SeismicUnix file [filename] The io-unit is given by [lu].
!							Endianness of the file is given by the optional argument [endian] (not present
!							or 0 = system standard, 1 = convert to little_endian, 2 = convert to big_endian).
!  deallocDataSu(this)				-	deallocates all previously allocated memory for dataSu object [this].
!
!  functions (operators are defined for each function):
!
!  getNumTrsDataSu(this)			-	Yields number of traces in dataSu object [this].
!  get[..]DataSu(this,ni)			-	Yields item [..] (see definition of trace_su object below) of trace [ni] in
!							dataSu object [this]
!  
!-------------------------------------------------------------
module dataSu
!
  implicit none
!
!-------------------------------------------------------------
!
  type trace_su
    private
    logical          :: allocated_by_readDataSu  = .false.  ! this is for correctly deallocating/disassociating the trace
    integer (kind=4) :: ni           ! SU number of trace
    integer (kind=2) :: trid         ! SU trace identification code
    integer (kind=4) :: recele       ! SU receiver group elevation
    integer (kind=4) :: surfele      ! SU surface elevation at source
    integer (kind=4) :: srcdepth     ! SU source depth below surface
    integer (kind=2) :: elepow       ! SU
    integer (kind=2) :: coordpow     ! SU
    integer (kind=2) :: timepow      ! SU
    integer (kind=4) :: xsrccoord    ! SU source coordinate x
    integer (kind=4) :: ysrccoord    ! SU source coordinate y
    integer (kind=4) :: xgrpcoord    ! SU group coordinate x
    integer (kind=4) :: ygrpcoord    ! SU group coordinate y
    integer (kind=2) :: coordunits   ! SU coordinate units
    integer (kind=2) :: srccorr      ! SU source static correction in milliseconds
    integer (kind=2) :: grpcorr      ! SU group static correction in milliseconds
    integer (kind=2) :: recdelay     ! SU delay recording time
    integer (kind=4) :: numsamp      ! SU number of samples in this trace
    integer (kind=4) :: sampint      ! SU sample interval in microseconds
    integer (kind=2) :: recyear      ! SU year data recorded
    integer (kind=2) :: doy          ! SU day of year
    integer (kind=2) :: hod          ! SU hour of minute
    integer (kind=2) :: moh          ! SU minute of hour
    integer (kind=2) :: som          ! SU second of minute
    integer (kind=2) :: timecode     ! SU time basis code
    integer (kind=2) :: trweight     ! SU trace weighting factor
    real (kind=4) :: d1              ! SU (commonly?) used for sampling time (dt) in seconds
    integer (kind=2) :: valunits     ! SU trace value measurement unit
    integer (kind=2) :: srctype      ! SU source type/orientation
    real (kind=4), dimension(:), pointer :: traceieee => null() ! SU 
  end type trace_su

  type data_su
    private
    integer (kind=4) :: n            ! number of traces
    type (trace_su), allocatable, dimension(:) :: trace_su
  end type data_su
!
  interface operator (.n.); module procedure getNumTrsDataSu; end interface
  interface operator (.ni.); module procedure getNumTrDataSu; end interface
  interface operator (.trid.); module procedure getTrIdDataSu; end interface
  interface operator (.recele.); module procedure getRecEleDataSu; end interface
  interface operator (.surfele.); module procedure getSurfEleDataSu; end interface
  interface operator (.srcdepth.); module procedure getSrcDepthDataSu; end interface
  interface operator (.elepow.); module procedure getElePowDataSu; end interface
  interface operator (.coordpow.); module procedure getCoordPowDataSu; end interface
  interface operator (.timepow.); module procedure getTimePowDataSu; end interface
  interface operator (.xsrccoord.); module procedure getXSrcCoordDataSu; end interface
  interface operator (.ysrccoord.); module procedure getYSrcCoordDataSu; end interface
  interface operator (.xgrpcoord.); module procedure getXGrpCoordDataSu; end interface
  interface operator (.ygrpcoord.); module procedure getYGrpCoordDataSu; end interface
  interface operator (.coordunits.); module procedure getCoordUnitsDataSu; end interface
  interface operator (.srccorr.); module procedure getSrcCorrDataSu; end interface
  interface operator (.grpcorr.); module procedure getGrpCorrDataSu; end interface
  interface operator (.recdelay.); module procedure getRecDelayDataSu; end interface
  interface operator (.numsamp.); module procedure getNumsampDataSu; end interface
  interface operator (.recyear.); module procedure getRecYearDataSu; end interface
  interface operator (.sampint.); module procedure getSampIntDataSu; end interface
  interface operator (.doy.); module procedure getDOYDataSu; end interface
  interface operator (.hod.); module procedure getHODDataSu; end interface
  interface operator (.moh.); module procedure getMOHDataSu; end interface
  interface operator (.som.); module procedure getSOMDataSu; end interface
  interface operator (.timecode.); module procedure getTimeCodeDataSu; end interface
  interface operator (.trweight.); module procedure getTrWeightDataSu; end interface
  interface operator (.done.); module procedure getD1DataSu; end interface     
  interface operator (.valunits.); module procedure getValUnitsDataSu; end interface
  interface operator (.srctype.); module procedure getSrcTypeDataSu; end interface
  interface operator (.traceieee.); module procedure getTraceIEEEDataSu; end interface
!
  contains
!
!-------------------------------------------------------------
!
  subroutine allocDataSu(this,n)
!
    type (data_su), intent(inout) :: this
    integer (kind=4), intent(in) :: n
    integer (kind=1) :: status
!
    if (allocated(this%trace_su)) call deallocDataSu(this)
    allocate(this%trace_su(n),stat=status)
    if (status /= 0) stop "createDataSu_mod_dataSu: unable to allocate traces"
    this%n = n
!
  end subroutine allocDataSu
!
!-------------------------------------------------------------
!
  subroutine addTraceDataSu(this,ni,trace,sampint,trid_opt,recele_opt,surfele_opt,srcdepth_opt,elepow_opt, &
                          & coordpow_opt,timepow_opt,xsrccoord_opt,ysrccoord_opt,xgrpcoord_opt,ygrpcoord_opt, &
                          & coordunits_opt,srccorr_opt,grpcorr_opt,recdelay_opt,recyear_opt,doy_opt,hod_opt, &
                          & moh_opt,som_opt,timecode_opt,trweight_opt,d1_opt,valunits_opt,srctype_opt)
!
    type (data_su), intent(inout) :: this
    integer (kind=4), intent(in) :: ni
    real (kind=4), intent(in), dimension(:), target :: trace
    integer (kind=2), intent(in) :: sampint
    integer (kind=2), intent(in), optional :: trid_opt
    integer (kind=4), intent(in), optional :: recele_opt
    integer (kind=4), intent(in), optional :: surfele_opt
    integer (kind=4), intent(in), optional :: srcdepth_opt
    integer (kind=2), intent(in), optional :: elepow_opt
    integer (kind=2), intent(in), optional :: coordpow_opt
    integer (kind=2), intent(in), optional :: timepow_opt
    integer (kind=4), intent(in), optional :: xsrccoord_opt
    integer (kind=4), intent(in), optional :: ysrccoord_opt
    integer (kind=4), intent(in), optional :: xgrpcoord_opt
    integer (kind=4), intent(in), optional :: ygrpcoord_opt
    integer (kind=2), intent(in), optional :: coordunits_opt
    integer (kind=2), intent(in), optional :: srccorr_opt
    integer (kind=2), intent(in), optional :: grpcorr_opt
    integer (kind=2), intent(in), optional :: recdelay_opt
    integer (kind=2), intent(in), optional :: recyear_opt
    integer (kind=2), intent(in), optional :: doy_opt
    integer (kind=2), intent(in), optional :: hod_opt
    integer (kind=2), intent(in), optional :: moh_opt
    integer (kind=2), intent(in), optional :: som_opt
    integer (kind=2), intent(in), optional :: timecode_opt
    integer (kind=2), intent(in), optional :: trweight_opt
    real (kind=4), intent(in), optional :: d1_opt
    integer (kind=2), intent(in), optional :: valunits_opt
    integer (kind=2), intent(in), optional :: srctype_opt
!
    if (ni .gt. this%n) stop "addTraceDataSu_mod_dataSu: trace no. greater than maximum no. of traces in file"
!
    this%trace_su(ni)%sampint = sampint
!
    if (present(trid_opt)) then
      this%trace_su(ni)%trid = trid_opt
    else
      this%trace_su(ni)%trid = 0
    endif
!
    if (present(recele_opt)) then
      this%trace_su(ni)%recele = recele_opt
    else
      this%trace_su(ni)%recele = 0
    endif
!
    if (present(surfele_opt)) then
      this%trace_su(ni)%surfele = surfele_opt
    else
      this%trace_su(ni)%surfele = 0
    endif
!
    if (present(srcdepth_opt)) then
      this%trace_su(ni)%srcdepth = srcdepth_opt
    else
      this%trace_su(ni)%srcdepth = 0
    endif
!
    if (present(elepow_opt)) then
      this%trace_su(ni)%elepow = elepow_opt
    else
      this%trace_su(ni)%elepow = 1
    endif
!
    if (present(coordpow_opt)) then
      this%trace_su(ni)%coordpow = coordpow_opt
    else
      this%trace_su(ni)%coordpow = 1
    endif
!
    if (present(timepow_opt)) then
      this%trace_su(ni)%timepow = timepow_opt
    else
      this%trace_su(ni)%timepow = 1
    endif
!
    if (present(xsrccoord_opt)) then
      this%trace_su(ni)%xsrccoord = xsrccoord_opt
    else
      this%trace_su(ni)%xsrccoord = 0
    endif
!
    if (present(ysrccoord_opt)) then
      this%trace_su(ni)%ysrccoord = ysrccoord_opt
    else
      this%trace_su(ni)%ysrccoord = 0
    endif
!
    if (present(xgrpcoord_opt)) then
      this%trace_su(ni)%xgrpcoord = xgrpcoord_opt
    else
      this%trace_su(ni)%xgrpcoord = 0
    endif
!
    if (present(ygrpcoord_opt)) then
      this%trace_su(ni)%ygrpcoord = ygrpcoord_opt
    else
      this%trace_su(ni)%ygrpcoord = 0
    endif
!
    if (present(coordunits_opt)) then
      this%trace_su(ni)%coordunits = coordunits_opt
    else
      this%trace_su(ni)%coordunits = 1
    endif
!
    if (present(srccorr_opt)) then
      this%trace_su(ni)%srccorr = srccorr_opt
    else
      this%trace_su(ni)%srccorr = 0
    endif
!
    if (present(grpcorr_opt)) then
      this%trace_su(ni)%grpcorr = grpcorr_opt
    else
      this%trace_su(ni)%grpcorr = 0
    endif
!
    if (present(recdelay_opt)) then
      this%trace_su(ni)%recdelay = recdelay_opt
    else
      this%trace_su(ni)%recdelay = 0
    endif
!
    if (present(recyear_opt)) then
      this%trace_su(ni)%recyear = recyear_opt
    else
      this%trace_su(ni)%recyear = 0
    endif
!
    if (present(doy_opt)) then
      this%trace_su(ni)%doy = doy_opt
    else
      this%trace_su(ni)%doy = 0
    endif
!
    if (present(hod_opt)) then
      this%trace_su(ni)%hod = hod_opt
    else
      this%trace_su(ni)%hod = 0
    endif
!
    if (present(moh_opt)) then
      this%trace_su(ni)%moh = moh_opt
    else
      this%trace_su(ni)%moh = 0
    endif
!
    if (present(som_opt)) then
      this%trace_su(ni)%som = som_opt
    else
      this%trace_su(ni)%som = 0
    endif
!
    if (present(timecode_opt)) then
      this%trace_su(ni)%timecode = timecode_opt
    else
      this%trace_su(ni)%timecode = 1
    endif
!
    if (present(trweight_opt)) then
      this%trace_su(ni)%trweight = trweight_opt
    else
      this%trace_su(ni)%trweight = 0
    endif
!
    if (present(d1_opt)) then
      this%trace_su(ni)%d1 = d1_opt
    else
      this%trace_su(ni)%d1 = 0.
    endif
!
    if (present(valunits_opt)) then
      this%trace_su(ni)%valunits = valunits_opt
    else
      this%trace_su(ni)%valunits = 0
    endif
!
    if (present(srctype_opt)) then
      this%trace_su(ni)%srctype = srctype_opt
    else
      this%trace_su(ni)%srctype = 0
    endif
!
    this%trace_su(ni)%ni = ni
!
    this%trace_su(ni)%numsamp = size(trace)
!
    this%trace_su(ni)%traceieee => trace
!
    this%trace_su(ni)%allocated_by_readDataSu = .false.
!
  end subroutine addTraceDataSu
!
!-------------------------------------------------------------
!
  subroutine delTraceDataSu(this,ni)
!
    type (data_su), intent(inout) :: this
    integer (kind=4), intent(in) :: ni
!
    if(this%trace_su(ni)%allocated_by_readDataSu) deallocate(this%trace_su(ni)%traceieee)
    this%trace_su(ni)%traceieee => null()
!
    this%trace_su(ni)%allocated_by_readDataSu = .false.
!
  end subroutine delTraceDataSu
!
!-------------------------------------------------------------
!
  subroutine writeDataSu(this,lu,filename,endian)
!
    type (data_su), intent(inout) :: this
    integer (kind=4), intent(in) :: lu
    integer (kind=4), intent(in), optional :: endian
    character (len=*), intent(in) :: filename
    integer (kind=4) :: ni,ti,shift
    integer (kind=1) :: status
!
    if (present(endian) .and. endian .eq. 1) then
      open(unit=lu, file=filename, status="replace", access='stream', action="write", convert="little_endian", iostat=status)
      if (status /= 0) stop "writeDataSu_mod_dataSu: Unable to open output-file"
    elseif (present(endian) .and. endian .eq. 2) then
      open(unit=lu, file=filename, status="replace", access='stream', action="write", convert="big_endian", iostat=status)
      if (status /= 0) stop "writeDataSu_mod_dataSu: Unable to open output-file"
    elseif (.not. present(endian) .or. present(endian) .and. endian .eq. 0) then
      open(unit=lu, file=filename, status="replace", access='stream', action="write", iostat=status)
      if (status /= 0) stop "writeDataSu_mod_dataSu: Unable to open output-file"
    else
      stop "writeDataSu_mod_dataSu: endian parameter has to be equal to 0,1 or 2"
    end if
!
    shift = 0
    do ni=1,this%n
      if (associated(this%trace_su(ni)%traceieee)) then
        write(lu,pos=13+shift) this%trace_su(ni)%ni
        write(lu,pos=29+shift) this%trace_su(ni)%trid
        write(lu,pos=41+shift) this%trace_su(ni)%recele
        write(lu,pos=45+shift) this%trace_su(ni)%surfele
        write(lu,pos=49+shift) this%trace_su(ni)%srcdepth
        write(lu,pos=69+shift) this%trace_su(ni)%elepow
        write(lu,pos=71+shift) this%trace_su(ni)%coordpow
        write(lu,pos=73+shift) this%trace_su(ni)%xsrccoord
        write(lu,pos=77+shift) this%trace_su(ni)%ysrccoord
        write(lu,pos=81+shift) this%trace_su(ni)%xgrpcoord
        write(lu,pos=85+shift) this%trace_su(ni)%ygrpcoord
        write(lu,pos=89+shift) this%trace_su(ni)%coordunits
        write(lu,pos=99+shift) this%trace_su(ni)%srccorr
        write(lu,pos=101+shift) this%trace_su(ni)%grpcorr
        write(lu,pos=109+shift) this%trace_su(ni)%recdelay
        if(this%trace_su(ni)%numsamp > 32767) then ! su format requires "unsigned short" for numsamp
           write(lu,pos=115+shift) int(this%trace_su(ni)%numsamp-65536,2)
        else
           write(lu,pos=115+shift) int(this%trace_su(ni)%numsamp,2)
        end if
        if(this%trace_su(ni)%sampint > 32767) then ! su format requires "unsigned short" for sampint
           write(lu,pos=117+shift) int(this%trace_su(ni)%sampint-65536,2)
        else
           write(lu,pos=117+shift) int(this%trace_su(ni)%sampint,2)
        end if
        write(lu,pos=157+shift) this%trace_su(ni)%recyear
        write(lu,pos=159+shift) this%trace_su(ni)%doy
        write(lu,pos=161+shift) this%trace_su(ni)%hod
        write(lu,pos=163+shift) this%trace_su(ni)%moh
        write(lu,pos=165+shift) this%trace_su(ni)%som
        write(lu,pos=167+shift) this%trace_su(ni)%timecode
        write(lu,pos=169+shift) this%trace_su(ni)%trweight
        write(lu,pos=181+shift) this%trace_su(ni)%d1
        write(lu,pos=203+shift) this%trace_su(ni)%valunits
        write(lu,pos=215+shift) this%trace_su(ni)%timepow
        write(lu,pos=217+shift) this%trace_su(ni)%srctype
        do ti=1,this%trace_su(ni)%numsamp
          write(lu,pos=241+shift+(ti-1)*4) this%trace_su(ni)%traceieee(ti)
        enddo
        shift = shift + 240 + this%trace_su(ni)%numsamp*4
      else
        stop "writeDataSu_mod_dataSu: data_su structure contains no data for at least one trace"
      endif
    enddo
!
    close (lu, status="keep", iostat=status)
    if (status /= 0) stop "mod dataSu: Unable to close output-file"
!
  end subroutine writeDataSu
!
!-------------------------------------------------------------
!
  subroutine readDataSu(this,lu,filename,endian)
!
    type(data_su), intent(inout) :: this
    integer (kind=4) :: lu
    integer (kind=4), intent(in), optional :: endian
    character (len=*), intent(in) :: filename
    integer (kind=4) :: ni,ti,shift,numsamp
    integer (kind=2) :: numsamp_read,sampint_read
    integer (kind=1) :: status
!
    if (present(endian) .and. endian .eq. 1) then
      open(unit=lu, file=filename, status="old", access='stream', action="read", convert="little_endian", iostat=status)
      if (status /= 0) stop "readDataSu_mod_dataSu: Unable to open output-file"
    elseif (present(endian) .and. endian .eq. 2) then
      open(unit=lu, file=filename, status="old", access='stream', action="read", convert="big_endian", iostat=status)
      if (status /= 0) stop "readDataSu_mod_dataSu: Unable to open output-file"
    elseif (.not. present(endian) .or. present(endian) .and. endian .eq. 0) then
      open(unit=lu, file=filename, status="old", access='stream', action="read", iostat=status)
      if (status /= 0) stop "readDataSu_mod_dataSu: Unable to open output-file"
    else
      stop "readDataSu_mod_dataSu: endian parameter has to be equal to 0,1 or 2"
    end if
!
    this%n = 0
    shift = 0
    do
      shift = shift + 115
      read(lu, pos=shift,iostat=status) numsamp_read
      if(numsamp_read < 0) then 
         numsamp = numsamp_read + 65536
      else
         numsamp = numsamp_read
      end if
      if (status .ne. 0) exit
      this%n = this%n + 1
      shift = shift + 125 + 4*numsamp
    enddo
!
    if (allocated(this%trace_su)) call deallocDataSu(this)
    allocate(this%trace_su(this%n),stat=status)
    if (status /= 0) stop "readDataSu_mod_dataSu: unable to allocate traces"
!
    shift = 0
    do ni=1,this%n
      read(lu,pos=13+shift) this%trace_su(ni)%ni
      read(lu,pos=29+shift) this%trace_su(ni)%trid
      read(lu,pos=41+shift) this%trace_su(ni)%recele
      read(lu,pos=45+shift) this%trace_su(ni)%surfele
      read(lu,pos=49+shift) this%trace_su(ni)%srcdepth
      read(lu,pos=69+shift) this%trace_su(ni)%elepow
      read(lu,pos=71+shift) this%trace_su(ni)%coordpow
      read(lu,pos=73+shift) this%trace_su(ni)%xsrccoord
      read(lu,pos=77+shift) this%trace_su(ni)%ysrccoord
      read(lu,pos=81+shift) this%trace_su(ni)%xgrpcoord
      read(lu,pos=85+shift) this%trace_su(ni)%ygrpcoord
      read(lu,pos=89+shift) this%trace_su(ni)%coordunits
      read(lu,pos=99+shift) this%trace_su(ni)%srccorr
      read(lu,pos=101+shift) this%trace_su(ni)%grpcorr
      read(lu,pos=109+shift) this%trace_su(ni)%recdelay
      read(lu,pos=115+shift) numsamp_read
      if(numsamp_read < 0) then 
         this%trace_su(ni)%numsamp = numsamp_read + 65536
      else
         this%trace_su(ni)%numsamp = numsamp_read
      end if
      read(lu,pos=117+shift) sampint_read
      if(sampint_read < 0) then 
         this%trace_su(ni)%sampint = sampint_read + 65536
      else
         this%trace_su(ni)%sampint = sampint_read
      end if
      read(lu,pos=157+shift) this%trace_su(ni)%recyear
      read(lu,pos=159+shift) this%trace_su(ni)%doy
      read(lu,pos=161+shift) this%trace_su(ni)%hod
      read(lu,pos=163+shift) this%trace_su(ni)%moh
      read(lu,pos=165+shift) this%trace_su(ni)%som
      read(lu,pos=167+shift) this%trace_su(ni)%timecode
      read(lu,pos=169+shift) this%trace_su(ni)%trweight
      read(lu,pos=181+shift) this%trace_su(ni)%d1
      read(lu,pos=203+shift) this%trace_su(ni)%valunits
      read(lu,pos=215+shift) this%trace_su(ni)%timepow
      read(lu,pos=217+shift) this%trace_su(ni)%srctype
      allocate(this%trace_su(ni)%traceieee(this%trace_su(ni)%numsamp),stat=status)
      if (status /= 0) stop "readDataSu_mod_dataSu: unable to allocate traces"
      this%trace_su(ni)%allocated_by_readDataSu = .true.
      do ti=1,this%trace_su(ni)%numsamp
        read(lu,pos=241+shift+(ti-1)*4) this%trace_su(ni)%traceieee(ti)
      enddo
      shift = shift + 240 + this%trace_su(ni)%numsamp*4
    enddo
!
    close (lu, status="keep", iostat=status)
    if (status /= 0) stop "mod dataSu: Unable to close input-file"
!
  end subroutine readDataSu
!
!-------------------------------------------------------------
!
  subroutine deallocDataSu(this)
!
    type (data_su), intent(inout) :: this
    integer (kind=4) :: ni
    integer (kind=1) :: status
!
    if (allocated(this%trace_su)) then
      do ni=1,this%n
        call delTraceDataSu(this,ni)
      enddo
      deallocate(this%trace_su,stat=status)
      if (status /= 0) stop "deallocDataSu_mod_dataSu: unable to deallocate traces"
    endif
  end subroutine deallocDataSu
!
!-------------------------------------------------------------
!
  integer (kind = 4) function getNumTrsDataSu(this)
    type (data_su), intent(in) :: this
    getNumTrsDataSu = this%n
  end function getNumTrsDataSu
!
  integer (kind = 4) function getNumTrDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getNumTrDataSu = this%trace_su(ni)%ni
  end function getNumTrDataSu
!
  integer (kind = 2) function getTrIdDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getTrIdDataSu = this%trace_su(ni)%trid
  end function getTrIdDataSu
!
  integer (kind = 4) function getRecEleDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getRecEleDataSu = this%trace_su(ni)%recele
  end function getRecEleDataSu
!
  integer (kind = 4) function getSurfEleDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getSurfEleDataSu = this%trace_su(ni)%surfele
  end function getSurfEleDataSu
!
  integer (kind = 4) function getSrcDepthDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getSrcDepthDataSu = this%trace_su(ni)%srcdepth
  end function getSrcDepthDataSu
!
  integer (kind = 2) function getElePowDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getElePowDataSu = this%trace_su(ni)%elepow
  end function getElePowDataSu
!
  integer (kind = 2) function getCoordPowDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getCoordPowDataSu = this%trace_su(ni)%coordpow
  end function getCoordPowDataSu
!
  integer (kind = 2) function getTimePowDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind=4), intent(in) :: ni
    getTimePowDataSu = this%trace_su(ni)%timepow
  end function getTimePowDataSu
!
  integer (kind = 4) function getXSrcCoordDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getXSrcCoordDataSu = this%trace_su(ni)%xsrccoord
  end function getXSrcCoordDataSu
!
  integer (kind = 4) function getYSrcCoordDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getYSrcCoordDataSu = this%trace_su(ni)%ysrccoord
  end function getYSrcCoordDataSu
!
  integer (kind = 4) function getXGrpCoordDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getXGrpCoordDataSu = this%trace_su(ni)%xgrpcoord
  end function getXGrpCoordDataSu
!
  integer (kind = 4) function getYGrpCoordDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getYGrpCoordDataSu = this%trace_su(ni)%ygrpcoord
  end function getYGrpCoordDataSu
!
  integer (kind = 2) function getCoordUnitsDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getCoordUnitsDataSu = this%trace_su(ni)%coordunits
  end function getCoordUnitsDataSu
!
  integer (kind = 2) function getSrcCorrDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getSrcCorrDataSu = this%trace_su(ni)%srccorr
  end function getSrcCorrDataSu
!
  integer (kind = 2) function getGrpCorrDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getGrpCorrDataSu = this%trace_su(ni)%grpcorr
  end function getGrpCorrDataSu
!
  integer (kind = 2) function getRecDelayDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getRecDelayDataSu = this%trace_su(ni)%recdelay
  end function getRecDelayDataSu
!
  integer (kind = 4) function getNumsampDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getNumsampDataSu = this%trace_su(ni)%numsamp
  end function getNumsampDataSu
!
  integer (kind = 4) function getSampIntDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getSampIntDataSu = this%trace_su(ni)%sampint
  end function getSampIntDataSu
!
  integer (kind = 2) function getRecYearDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getRecYearDataSu = this%trace_su(ni)%recyear
  end function getRecYearDataSu
!
  integer (kind = 2) function getDOYDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getDOYDataSu = this%trace_su(ni)%doy
  end function getDOYDataSu
!
  integer (kind = 2) function getHODDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getHODDataSu = this%trace_su(ni)%hod
  end function getHODDataSu
!
  integer (kind = 2) function getMOHDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getMOHDataSu = this%trace_su(ni)%moh
  end function getMOHDataSu
!
  integer (kind = 2) function getSOMDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getSOMDataSu = this%trace_su(ni)%som
  end function getSOMDataSu
!
  integer (kind = 2) function getTimeCodeDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getTimeCodeDataSu = this%trace_su(ni)%timecode
  end function getTimeCodeDataSu
!
  integer (kind = 2) function getTrWeightDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getTrWeightDataSu = this%trace_su(ni)%trweight
  end function getTrWeightDataSu
!
  real (kind=4) function getD1DataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getD1DataSu = this%trace_su(ni)%d1
  end function getD1DataSu
!
  integer (kind = 2) function getValUnitsDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getValUnitsDataSu = this%trace_su(ni)%valunits
  end function getValUnitsDataSu
!
  integer (kind = 2) function getSrcTypeDataSu(this,ni)
    type (data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    getSrcTypeDataSu = this%trace_su(ni)%srctype
  end function getSrcTypeDataSu
!
  function getTraceIEEEDataSu(this,ni) result(p)
    type(data_su), intent(in) :: this
    integer (kind = 4), intent(in) :: ni
    real (kind = 4), dimension(:), pointer :: p
    p => this%trace_su(ni)%traceieee
  end function getTraceIEEEDataSu
!
!-------------------------------------------------------------
!
end module dataSu
