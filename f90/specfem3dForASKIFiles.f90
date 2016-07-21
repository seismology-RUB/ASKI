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

module specfem3dForASKIFiles

  use vectorPointer
  use errorMessage

  implicit none

  integer, parameter :: length_ID_specfem3d_for_ASKI_files = 13 !< change this number CONSISTENTLY with length_ASKI_output_ID in SPECFEM3D (GLOBE/Cartesian) codes
  integer, parameter :: version_globe_specfem3d_for_ASKI_files = 1
  integer, parameter :: version_cartesian_specfem3d_for_ASKI_files = 2

contains

!------------------------------------------------------------------------
!> \brief read in specfem3d for ASKI main file
!! \param filename file name
!! \param lu file unit
!! \param errmsg error message
!! \param specfem_version specfem version, to be returned optionally
!! \param len_file_ID length of file ID character, to be returned optionally
!! \param file_ID file ID, to be returned optionally
!! \param nproc number of procs, to be returned optionally
!! \param type_inversion_grid type of inversion grid, to be returned optionally
!! \param nwp number of wavefield points, to be returned optionally
!! \param df frequency step, to be returned optionally
!! \param nf number of frequencies, to be returned optionally
!! \param jf frequency indices, to be returned optionally
!! \param x x-coordinates of wavefield points, to be returned optionally
!! \param y y-coordinates of wavefield points, to be returned optionally
!! \param z z-coordinates of wavefield points, to be returned optionally
!! \param rho density reference model, to be returned optionally
!! \param vp vp reference model, to be returned optionally
!! \param vs vs reference model, to be returned optionally
!! \param ngllx ngllx (only if specfem3dInversionGrid), to be returned optionally
!! \param nglly nglly (only if specfem3dInversionGrid), to be returned optionally
!! \param ngllz ngllz (only if specfem3dInversionGrid), to be returned optionally
!! \param jacobian jacobians on wavefield points (only if specfem3dInversionGrid), to be returned optionally
!! \param ncell number of invgrid cells (only if specfem3dInversionGrid), to be returned optionally
!! \param nb_idx invgrid cell neighbours (only if specfem3dInversionGrid), to be returned optionally
  subroutine readSpecfem3dForASKIMainFile(filename,lu,errmsg,&
       specfem_version,len_file_ID,file_ID,nproc,type_inversion_grid,nwp,df,nf,jf,&
       x,y,z,rho,vp,vs,ngllx,nglly,ngllz,jacobian,ncell,nb_idx)
    ! incoming
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! optional return
    integer, optional :: specfem_version,len_file_ID,nproc,type_inversion_grid,nwp,nf,ngllx,nglly,ngllz,ncell
    character(len=length_ID_specfem3d_for_ASKI_files), optional :: file_ID
    double precision, optional :: df
    integer, dimension(:), pointer, optional :: jf
    real, dimension(:), pointer, optional :: x,y,z,rho,vp,vs,jacobian
    type (integer_vector_pointer), dimension(:), pointer, optional :: nb_idx
    ! local
    integer :: idummy
    double precision :: ddummy
    character(len=length_ID_specfem3d_for_ASKI_files) :: ID_dummy
    integer, dimension(:), pointer :: ip_dummy
    real, dimension(:), pointer :: rp_dummy
    type (integer_vector_pointer), dimension(:), pointer :: nb_idx_dummy
    integer :: npresent,ier,nnb_total,nnb_count,icell,&
         specfem_version_local,type_inversion_grid_local,nf_local,nwp_local,&
         ngllx_local,nglly_local,ngllz_local,ncell_local
    character(len=28) :: myname = 'readSpecfem3dForASKIMainFile'
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
!
    ! Find out the total number of optional parameters that are requested from this file (optional arguments after errmsg)
    ! Decrease this number every time a requested parameter was read from file, return ones this number reaches 0
    npresent = 0
    if(present(specfem_version)) npresent = npresent + 1 
    if(present(len_file_ID)) npresent = npresent + 1 
    if(present(file_ID)) npresent = npresent + 1 
    if(present(nproc)) npresent = npresent + 1 
    if(present(type_inversion_grid)) npresent = npresent + 1 
    if(present(nwp)) npresent = npresent + 1 
    if(present(df)) npresent = npresent + 1 
    if(present(nf)) npresent = npresent + 1 
    if(present(jf)) npresent = npresent + 1 
    if(present(x)) npresent = npresent + 1 
    if(present(y)) npresent = npresent + 1 
    if(present(z)) npresent = npresent + 1 
    if(present(rho)) npresent = npresent + 1 
    if(present(vp)) npresent = npresent + 1 
    if(present(vs)) npresent = npresent + 1 
    if(present(ngllx)) npresent = npresent + 1 
    if(present(nglly)) npresent = npresent + 1 
    if(present(ngllz)) npresent = npresent + 1 
    if(present(jacobian)) npresent = npresent + 1 
    if(present(ncell)) npresent = npresent + 1 
    if(present(nb_idx)) npresent = npresent + 1 
!
    if(npresent == 0) then
       call add(errmsg,1,"there are no incoming parameters that should be read from file, the call to this "//&
            "routine is redundant",myname)
       return
    end if
!
    ! open file
    open(unit=lu,file=filename,status='old',action='read',access='stream',form='unformatted',iostat=ier)
    if (ier /= 0) then
       write(errstr,*) "File '"//trim(filename)//"' cannot be opened to read, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    else
       call add(errmsg,0,"successfully opened file '"//trim(filename)//"' to read",myname)
    endif
!
    ! specfem_version
    read(lu,iostat=ier) specfem_version_local
    if(ier/=0) then
       write(errstr,*) "could not read specfem version from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    select case(specfem_version_local)
    case(1,2)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "specfem version = ",specfem_version_local,&
            " is not supported: 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select
    if(present(specfem_version)) then
       specfem_version = specfem_version_local
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! len_file_ID
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read length of file ID from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(idummy/=length_ID_specfem3d_for_ASKI_files) then
       write(errstr,*) "length of ASKI file ID is ",idummy,", only length ",&
            length_ID_specfem3d_for_ASKI_files," supported here. Please update this code accordingly."
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(len_file_ID)) then
       len_file_ID = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! file_ID
    read(lu,iostat=ier) ID_dummy
    if(ier/=0) then
       write(errstr,*) "could not read file ID from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(file_ID)) then
       file_ID = ID_dummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! nproc
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read numper of procs from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(nproc)) then
       nproc = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! type_inversion_grid
    read(lu,iostat=ier) type_inversion_grid_local
    if(ier/=0) then
       write(errstr,*) "could not read inversion grid type from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    ! check if specfem_version and type_inversion_grid are a valid match
    if(specfem_version_local == 1) then
       select case(type_inversion_grid_local)
       case(1,3,4,5)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_inversion_grid_local," is not supported by SPECFEM3D_GLOBE"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(specfem_version_local == 2) then
       select case(type_inversion_grid_local)
       case(2,3,4)
          ! OK, so pass doing nothing
       case default
          write(errstr,*) "type of inversion grid = ",type_inversion_grid_local," is not supported by SPECFEM3D_Cartesian"
          call add(errmsg,2,errstr,myname)
          goto 1
       end select
    end if
    if(present(type_inversion_grid)) then
       type_inversion_grid = type_inversion_grid_local
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! nwp
    read(lu,iostat=ier) nwp_local
    if(ier/=0) then
       write(errstr,*) "could not read total number of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(nwp_local < 1) then
       write(errstr,*) "number of wavefield points =",nwp_local,", should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(nwp)) then
       nwp = nwp_local
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! df
    read(lu,iostat=ier) ddummy
    if(ier/=0) then
       write(errstr,*) "could not read frequency step from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(df)) then
       df = ddummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! nf
    read(lu,iostat=ier) nf_local
    if(ier/=0) then
       write(errstr,*) "could not read number of frequencies from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(nf)) then
       nf = nf_local
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! jf
    allocate(ip_dummy(nf_local))
    read(lu,iostat=ier) ip_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nf_local," frequency indices from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(any(ip_dummy<0)) then
       call add(errmsg,2,"there are frequency indices in the file, which are < 0",myname)
       goto 1
    end if
    do idummy = 2,nf_local
       if(any(ip_dummy(1:idummy-1)==ip_dummy(idummy))) then
          write(errstr,*) idummy,"'th frequency index read from file = ",ip_dummy(idummy),&
               " occurrs twice in vector of frequency indices. Duplicate frequency ",&
               " indices inicate an inconsistency of the file."
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
    end do ! idummy
    if(present(jf)) then
       jf => ip_dummy
       nullify(ip_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
    else
       deallocate(ip_dummy)
    end if
!
    ! x
    allocate(rp_dummy(nwp_local))
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," x-coordinates of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(x)) then
       x => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
       ! since at this point there was no return, y must be read, too. so already allocate here for reading y
       allocate(rp_dummy(nwp_local))
    !else
       ! else: do nothing! 
       ! Do not deallocate here, as before for jf, but reuse rp_dummy for reading y.
       ! This way unnecessary (de)allocations are saved
    end if
!
    ! y
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," y-coordinates of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(y)) then
       y => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
       ! since at this point there was no return, z must be read, too. so already allocate here for reading z
       allocate(rp_dummy(nwp_local))
    !else
       ! else: do nothing! 
       ! Do not deallocate here, as before for jf, but reuse rp_dummy for reading z.
       ! This way unnecessary (de)allocations are saved
    end if
!
    ! z
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," z-coordinates of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(z)) then
       z => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
       ! since at this point there was no return, rho must be read, too. so already allocate here for reading rho
       allocate(rp_dummy(nwp_local))
    !else
       ! else: do nothing! 
       ! Do not deallocate here, as before for jf, but reuse rp_dummy for reading rho.
       ! This way unnecessary (de)allocations are saved
    end if
!
    ! rho
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," rho values from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(rho)) then
       rho => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
       ! since at this point there was no return, vp must be read, too. so already allocate here for reading vp
       allocate(rp_dummy(nwp_local))
    !else
       ! else: do nothing! 
       ! Do not deallocate here, as before for jf, but reuse rp_dummy for reading vp.
       ! This way unnecessary (de)allocations are saved
    end if
!
    ! vp
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," vp values from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(vp)) then
       vp => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
       ! since at this point there was no return, vs must be read, too. so already allocate here for reading vs
       allocate(rp_dummy(nwp_local))
    !else
       ! else: do nothing! 
       ! Do not deallocate here, as before for jf, but reuse rp_dummy for reading vs.
       ! This way unnecessary (de)allocations are saved
    end if
!
    ! vs
    read(lu,iostat=ier) rp_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," vs values from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(vs)) then
       vs => rp_dummy
       nullify(rp_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
    else
       deallocate(rp_dummy)
    end if
!
    select case(type_inversion_grid_local)
    case(4) ! specfem3dInversionGrid
       ! ngllx
       read(lu,iostat=ier) ngllx_local
       if(ier/=0) then
          write(errstr,*) "could not read ngllx from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(ngllx)) then
          ngllx = ngllx_local
          npresent = npresent -1
          if(npresent == 0) goto 1
       end if
!
       ! nglly
       read(lu,iostat=ier) nglly_local
       if(ier/=0) then
          write(errstr,*) "could not read nglly from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(nglly)) then
          nglly = nglly_local
          npresent = npresent -1
          if(npresent == 0) goto 1
       end if
!
       ! ngllz
       read(lu,iostat=ier) ngllz_local
       if(ier/=0) then
          write(errstr,*) "could not read ngllz from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(ngllz)) then
          ngllz = ngllz_local
          npresent = npresent -1
          if(npresent == 0) goto 1
       end if
!
       ! jacobian
       allocate(rp_dummy(nwp_local))
       read(lu,iostat=ier) rp_dummy
       if(ier/=0) then
          write(errstr,*) "could not read ",nwp_local," jacobian values from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(jacobian)) then
          jacobian => rp_dummy
          nullify(rp_dummy)
          npresent = npresent -1
          if(npresent == 0) goto 1
       else
          deallocate(rp_dummy)
       end if
!
       ! ncell
       read(lu,iostat=ier) ncell_local
       if(ier/=0) then
          write(errstr,*) "could not read ncell from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(ncell_local<1 .or. ncell_local /= nwp_local/(ngllx_local*nglly_local*ngllz_local)) then
          write(errstr,*) "ncell = ",ncell_local," does not compute as ntot_wp/(NGLLX*NGLLY*NGLLZ) = ",&
               nwp_local/(ngllx_local*nglly_local*ngllz_local)," or it is negative -> file is inconsistent!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(ncell)) then
          ncell = ncell_local
          npresent = npresent -1
          if(npresent == 0) goto 1
       end if
!
       ! nb_idx
       read(lu,iostat=ier) nnb_total
       if(ier/=0) then
          write(errstr,*) "could not read total number of integers defining invgrid neighbours from file, raised iostat = ",ier
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(nnb_total<ncell_local) then
          write(errstr,*) "total number of integers defining the neighbours ",nnb_total,&
               " should be at least the number of cells ",ncell_local," -> file is inconsistent!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       allocate(nb_idx_dummy(ncell_local))
       nnb_count = 0
       do icell = 1,ncell_local
          read(lu,iostat=ier) idummy
          if(ier/=0) then
             write(errstr,*) "could not read number of neighbours for ",icell,"'th cell from file, raised iostat = ",ier
             call add(errmsg,2,errstr,myname)
             goto 1
          end if
          nnb_count = nnb_count + 1
          if(idummy>0) then
             allocate(ip_dummy(idummy))
             read(lu,iostat=ier) ip_dummy
             if(ier/=0) then
                write(errstr,*) "could not read neighbours for ",icell,"'th cell from file, raised iostat = ",ier
                call add(errmsg,2,errstr,myname)
                goto 1
             end if
             call associateVectorPointer(nb_idx_dummy(icell),ip_dummy); nullify(ip_dummy)
             nnb_count = nnb_count + idummy
          end if
       end do ! icell
       if(nnb_count /= nnb_total) then
          write(errstr,*) "the total number of integers in main file containing neighbour information ",nnb_count,&
               " does not match the expected number of integers ",nnb_total," -> file is inconsistent!"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
       if(present(nb_idx)) then
          nb_idx => nb_idx_dummy
          nullify(nb_idx_dummy)
          npresent = npresent -1
          if(npresent == 0) goto 1
       else
          do icell = 1,size(nb_idx_dummy)
             call dealloc(nb_idx_dummy(icell))
          end do
          deallocate(nb_idx_dummy)
       end if
    case default
       if(present(ngllx).or.present(nglly).or.present(ngllz).or.present(jacobian).or.present(ncell).or.present(nb_idx)) then
          write(errstr,*) "ngllx,nglly,ngllz,jacobian,ncell or nb_idx were requested to be read from this main file, although ",&
               "it does not contain any specfem3d inversion grid information"
          call add(errmsg,2,errstr,myname)
          goto 1
       end if
    end select ! case type_inversion_grid_local
!
    ! the code should actually not come here, since above npresent should have reached 0 at some point and
    ! goto 1 statement should have been called
    if(npresent /= 0) then
       write(errstr,*) "after reading the complete file, there are still ",npresent,&
            " requested parameters which were not yet returned. This number must be 0 at this point and this ",&
            "error should not occur!!"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
1   close(lu)
    if(associated(ip_dummy)) deallocate(ip_dummy)
    if(associated(rp_dummy)) deallocate(rp_dummy)
    if(associated(nb_idx_dummy)) then
       do icell = 1,size(nb_idx_dummy)
          call dealloc(nb_idx_dummy(icell))
       end do
    end if
  end subroutine readSpecfem3dForASKIMainFile
!
!------------------------------------------------------------------------
!> \brief read in specfem3d for ASKI main file
!! \param filename file name
!! \param lu file unit
!! \param errmsg error message
!! \param specfem_version specfem version, to be returned optionally
!! \param len_file_ID length of file ID character, to be returned optionally
!! \param file_ID file ID, to be returned optionally
!! \param nproc number of procs, to be returned optionally
!! \param nwp number of wavefield points, to be returned optionally
!! \param df frequency step, to be returned optionally
!! \param jf frequency index for which spectra are contained in the file, to be returned optionally
!! \param u displacement spectrum, to be returned optionally
!! \param ustr strain spectrum, to be returned optionally
  subroutine readSpecfem3dForASKISpectralWavefieldFile(filename,lu,errmsg,specfem_version,len_file_ID,file_ID,&
       nproc,nwp,df,jfcur,u,ustr)
    character(len=*) :: filename
    integer :: lu
    type (error_message) :: errmsg
    ! optional return
    integer, optional :: specfem_version,len_file_ID,nproc,nwp,jfcur
    character(len=length_ID_specfem3d_for_ASKI_files), optional :: file_ID
    double precision, optional :: df
    complex, dimension(:,:), pointer, optional :: u,ustr
    ! local
    integer :: idummy
    double precision :: ddummy
    character(len=length_ID_specfem3d_for_ASKI_files) :: ID_dummy
    complex, dimension(:,:), pointer :: u_dummy
    integer :: npresent,ier,nwp_local
    character(len=41) :: myname = 'readSpecfem3dForASKISpectralWavefieldFile'
    character(len=400) :: errstr
!
    call addTrace(errmsg,myname)
!
    ! Find out the total number of optional parameters that are requested from this file (optional arguments after errmsg)
    ! Decrease this number every time a requested parameter was read from file, return ones this number reaches 0
    npresent = 0
    if(present(specfem_version)) npresent = npresent + 1 
    if(present(len_file_ID)) npresent = npresent + 1 
    if(present(file_ID)) npresent = npresent + 1 
    if(present(nproc)) npresent = npresent + 1 
    if(present(nwp)) npresent = npresent + 1 
    if(present(df)) npresent = npresent + 1 
    if(present(jfcur)) npresent = npresent + 1 
    if(present(u)) npresent = npresent + 1 
    if(present(ustr)) npresent = npresent + 1 
!
    if(npresent == 0) then
       call add(errmsg,1,"there are no incoming parameters that should be read from file, the call to this "//&
            "routine is redundant",myname)
       return
    end if
!
    ! open file
    open(unit=lu,file=filename,status='old',action='read',access='stream',form='unformatted',iostat=ier)
    if (ier /= 0) then
       write(errstr,*) "File '"//trim(filename)//"' cannot be opened to read, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    else
       call add(errmsg,0,"successfully opened file '"//trim(filename)//"' to read",myname)
    endif
!
    ! specfem_version
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read specfem version from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    select case(idummy)
    case(1,2)
       ! OK, so pass doing nothing
    case default
       write(errstr,*) "specfem version = ",idummy,&
            " is not supported: 1 = SPECFEM3D_GLOBE, 2 = SPECFEM3D_Cartesian"
       call add(errmsg,2,errstr,myname)
       goto 1
    end select
    if(present(specfem_version)) then
       specfem_version = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! len_file_ID
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read length of file ID from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(idummy/=length_ID_specfem3d_for_ASKI_files) then
       write(errstr,*) "length of ASKI file ID is ",idummy,", only length ",&
            length_ID_specfem3d_for_ASKI_files," supported here. Please update this code accordingly."
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(len_file_ID)) then
       len_file_ID = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! file_ID
    read(lu,iostat=ier) ID_dummy
    if(ier/=0) then
       write(errstr,*) "could not read file ID from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(file_ID)) then
       file_ID = ID_dummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! nproc
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read numper of procs from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(nproc)) then
       nproc = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! nwp
    read(lu,iostat=ier) nwp_local
    if(ier/=0) then
       write(errstr,*) "could not read total number of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(nwp_local < 1) then
       write(errstr,*) "number of wavefield points =",nwp_local,", should be positive"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(nwp)) then
       nwp = nwp_local
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! df
    read(lu,iostat=ier) ddummy
    if(ier/=0) then
       write(errstr,*) "could not read frequency step from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(df)) then
       df = ddummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! jfcur
    read(lu,iostat=ier) idummy
    if(ier/=0) then
       write(errstr,*) "could not read total number of wavefield points from file, raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(jfcur)) then
       jfcur = idummy
       npresent = npresent -1
       if(npresent == 0) goto 1
    end if
!
    ! u
    allocate(u_dummy(nwp_local,3))
    read(lu,iostat=ier) u_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," times 3 values of spectral displacement from file, ",&
            "raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(u)) then
       u => u_dummy
       nullify(u_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
    else
       deallocate(u_dummy)
    end if
!
    ! ustr
    allocate(u_dummy(nwp_local,6))
    read(lu,iostat=ier) u_dummy
    if(ier/=0) then
       write(errstr,*) "could not read ",nwp_local," times 6 values of spectral strains from file, ",&
            "raised iostat = ",ier
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
    if(present(ustr)) then
       ustr => u_dummy
       nullify(u_dummy)
       npresent = npresent -1
       if(npresent == 0) goto 1
    else
       deallocate(u_dummy)
    end if
!
    ! the code should actually not come here, since above npresent should have reached 0 at some point and
    ! goto 1 statement should have been called
    if(npresent /= 0) then
       write(errstr,*) "after reading the complete file, there are still ",npresent,&
            " requested parameters which were not yet returned. This number must be 0 at this point and this ",&
            "error should not occur!!"
       call add(errmsg,2,errstr,myname)
       goto 1
    end if
!
1   close(lu)
    if(associated(u_dummy)) deallocate(u_dummy)
  end subroutine readSpecfem3dForASKISpectralWavefieldFile
  
end module specfem3dForASKIFiles
