!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of ASKI version 1.1.
!
!   ASKI version 1.1 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   ASKI version 1.1 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with ASKI version 1.1.  If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!> \brief Handles definition of and transformation between different components of orthogonal bases (e.g. directions of data components)
!!
!! \details In this module, components of different local (on the sphere) or global (orthogonal) basis systems 
!!  (e.g. components, i.e. local trace directions of seismic data) may be defined (named) along with their 
!!  relation to other components, i.e. transformations between them. The transformation coefficients, i.e. 
!!  (parts of) the transformation matrix, from one set of components C_in to another C_out, is computed by
!!  transforming C_in to global Cartesian coordinates CX,CY,CZ and then CX,CY,CZ to C_out
!!
!! \author Florian Schumacher
!! \date August 2012
!
module componentTransformation
!
    use seismicStation
    use mathConstants
!
    implicit none
!
    private :: indxComponent,componentByIndex,number_of_components
!
    interface createComponentTransformation
       module procedure setCoefficientsCartesianComponentTransformation
       module procedure setCoefficientsSphericalComponentTransformation
    end interface createComponentTransformation
    interface transform
       module procedure getCoefficientsCartesianComponentTransformation
       module procedure getCoefficientsByIndexComponentTransformation
       module procedure getCoefficientsByNameComponentTransformation
    end interface transform
    interface dealloc; module procedure deallocateCoefficientsComponentTransformation; end interface
!
    !> \brief stores transformation coefficients between data components
    type component_transformation
       private
       character(len=1) :: csys = '' !< character indicating coordinate system in which the components are used
       integer :: nstat = 0 !< number of stations (i.e. size of 3rd dimension of arrays C2XYZ,XYZ2C) (also indicates if object was created already)
       double precision, dimension(:,:,:), pointer :: C2XYZ => null() !< array containing the coefficients needed to transform components to X,Y,Z
       double precision, dimension(:,:,:), pointer :: XYZ2C => null() !< array containing the coefficients needed to transform X,Y,Z to components
       character(len=character_length_staname), dimension(:), pointer :: staname => null() !< in case csys='S', contains respective station names (nstat many)
    end type component_transformation
!
    integer, parameter :: character_length_component = 4
    character(len=24), parameter :: all_valid_components = 'CX,CY,CZ,N,S,E,W,UP,DOWN' !< string containing all valid components (for output)
    integer, parameter :: number_of_components = 9
!
contains
!------------------------------------------------------------------------
!> \brief map component character to internal index
!! \param C component character
!! \param i index
!! \return index (will be -1 if character not supported)
!
  function indxComponent(C) result(i)
    character(len=*), intent(in) :: C
    integer :: i
    select case (trim(C))
    case ('CX') ! Cartesian X (X axis through equator and zero meridian)
       i = 1
    case ('CY') ! Cartesian Y (Y axis through equator and 90 deg East meridian)
       i = 2
    case ('CZ') ! Cartesian Z (Z axis through north pole)
       i = 3
    case ('N') ! local north
       i = 4
    case ('S') ! local south
       i = 5
    case ('E') ! local east
       i = 6
    case ('W') ! local west
       i = 7
    case ('UP') ! local up
       i = 8
    case ('DOWN') ! local down
       i = 9
    ! EXAMPLE
    ! case ('HX') ! "half X" - defines the transformation CX <-> CX/2
    !    i = 10
    case default
       i = -1
    end select
  end function indxComponent
!------------------------------------------------------------------------
!> \brief map internal index to component character MUST BE CONSISTENT WITH function indxComponent
!! \param i index
!! \param C component character
!! \return character (will be '' if index out of range)
!
  function componentByIndex(i) result(C)
    integer, intent(in) :: i
    character(len=character_length_component) :: C
    select case (i)
    case (1) ! Cartesian X (X axis through equator and zero meridian)
       C = 'CX'
    case (2) ! Cartesian Y (Y axis through equator and 90 deg East meridian)
       C ='CY' 
    case (3) ! Cartesian Z (Z axis through north pole)
       C ='CZ' 
    case (4) ! local north
       C ='N' 
    case (5) ! local south
       C = 'S'
    case (6) ! local east
       C = 'E'
    case (7) ! local west
       C = 'W'
    case (8) ! local up
       C = 'UP'
    case (9) ! local down
       C = 'DOWN'
    ! EXAMPLE
    ! case (10) ! "half X" - defines the transformation CX <-> CX/2
    !    C = 'HX'
    case default
       C = ''
    end select
  end function componentByIndex
!------------------------------------------------------------------------
!> \brief iterator over all supported components
  logical function nextComponent(C,reset)
    character(len=character_length_component) :: C
    logical, optional :: reset
    ! local
    integer :: call_count = 0
    save :: call_count
!
    ! if this iterator is to be reset, do so
    if(present(reset)) then
       if(reset) goto 1
    end if
!
    ! increase counter
    call_count = call_count+1
!
    ! if the counter rises above the upper bound, reset this iterator
    if(call_count > number_of_components) goto 1
!
    ! otherwise set the return variable and indicate success
    C = componentByIndex(call_count)
    nextComponent = .true.
!
    ! IF FUNCTION COMES HERE, RETURN NORMALLY
    return
!
    ! RESET THE ITERATOR
1   call_count = 0
    C = ''
    nextComponent = .false.
    return
  end function nextComponent
!------------------------------------------------------------------------
!> \brief create transformation coefficients for Cartesian case
!! \param this component transformation object
!
  subroutine setCoefficientsCartesianComponentTransformation(this)
    type (component_transformation) :: this
    if(this%nstat /= 0) call deallocateCoefficientsComponentTransformation(this)
    this%csys = 'C'
    allocate(this%C2XYZ(3,number_of_components,1),this%XYZ2C(3,number_of_components,1))
    this%C2XYZ = 0.d0 ; this%XYZ2C = 0.d0
    !
    ! Set coefficients (i.e. entries of transformation matrices) to transform a specific (set of) component(s) to CX,CY,CZ.
    ! Any specific collection of columns of array this%C2XYZ forms a transformation matrix from those corresponding 
    ! components to CX,CY,CZ
    this%C2XYZ(:,indxComponent('CX'),1) = (/ 1.d0,0.d0,0.d0 /)
    this%C2XYZ(:,indxComponent('CY'),1) = (/ 0.d0,1.d0,0.d0 /)
    this%C2XYZ(:,indxComponent('CZ'),1) = (/ 0.d0,0.d0,1.d0 /)
    this%C2XYZ(:,indxComponent('N'),1) = (/ -1.d0,0.d0,0.d0 /)
    this%C2XYZ(:,indxComponent('S'),1) = (/ 1.d0,0.d0,0.d0 /)
    this%C2XYZ(:,indxComponent('E'),1) = (/ 0.d0,1.d0,0.d0 /)
    this%C2XYZ(:,indxComponent('W'),1) = (/ 0.d0,-1.d0,0.d0 /) ! e.g. the W-coordinate contributes with factor -1. to CY and with factor 0. to CX and CZ
    this%C2XYZ(:,indxComponent('UP'),1) = (/ 0.d0,0.d0,1.d0 /)
    this%C2XYZ(:,indxComponent('DOWN'),1) = (/ 0.d0,0.d0,-1.d0 /)
    ! EXAMPLE
    !    this%C2XYZ(:,indxComponent('HX'),1) = (/ 2.d0,0.d0,0.d0 /) ! the HX=CX/2-coordinate contributes with factor 2. to CX and  with factor 0. to CY and CZ
    !
    ! Set coefficients (i.e. entries of transformation matrices) to transform CX,CY,CZ to a specific (set of) component(s).
    ! Any specific collection of columns of array this%XYZ2C forms THE TRANSPOSE of a transformation matrix from CX,CY,CZ to 
    ! those corresponding components.
    ! Store in transposed from (i.e. CX,CY,CZ in first dimension), as for transformation CX,CY,CZ, -> C you always need all values
    ! of CX,CY,CZ for a specific component C (more efficient when accessing the array this%XYZ2C, compare transpose() statement 
    ! inside matmul in getCoefficients.. routines)
    ! IN CASE OF ORTHOGONAL TRANSFORMATIONS, the values of a specific column of arrays this%C2XYZ and this%XYZ2C will
    ! coincide, as in this case the inverse transformation matrix is its transposed!
    ! However, the coefficients this%XYZ2C (which correspond to the inverse of the transformation by this%C2XYZ) 
    ! are kept seperately from this%C2XYZ in order to also support non-orthogonal basis transformations in this module 
    ! (compare example component 'HX'=CX/2 ). 
    this%XYZ2C(:,indxComponent('CX'),1) = (/ 1.d0,0.d0,0.d0 /)
    this%XYZ2C(:,indxComponent('CY'),1) = (/ 0.d0,1.d0,0.d0 /)
    this%XYZ2C(:,indxComponent('CZ'),1) = (/ 0.d0,0.d0,1.d0 /)
    this%XYZ2C(:,indxComponent('N'),1) = (/ -1.d0,0.d0,0.d0 /)
    this%XYZ2C(:,indxComponent('S'),1) = (/ 1.d0,0.d0,0.d0 /)
    this%XYZ2C(:,indxComponent('E'),1) = (/ 0.d0,1.d0,0.d0 /)
    this%XYZ2C(:,indxComponent('W'),1) = (/ 0.d0,-1.d0,0.d0 /) ! e.g. W = 0.*X -1.*Y + 0.*Z
    this%XYZ2C(:,indxComponent('UP'),1) = (/ 0.d0,0.d0,1.d0 /)
    this%XYZ2C(:,indxComponent('DOWN'),1) = (/ 0.d0,0.d0,-1.d0 /)
    ! EXAMPLE
    !    this%XYZ2C(:,indxComponent('HX'),1) = (/ 0.5d0,0.d0,0.d0 /) ! X/2 = 0.5*X + 0.*Y + 0.*Z
!
    this%nstat = 1
  end subroutine setCoefficientsCartesianComponentTransformation
!------------------------------------------------------------------------
!> \brief create transformation coefficients for spherical case
!! \details here, also the cartesian case is supported: dependent on the coordinate system
!!  of the first incoming seismic station (.csys.stations(1)), subroutine
!!  setCoefficientsCartesianComponentTransformation is called instead. This tries to assure a
!!  proper behaviour of this module, if the program does not know (or does not cares about) the 
!!  coordinate system of the setup.
!! \param this component transformation object
!! \param stations array of seismic station objects for which the coefficients will be computed
!
  subroutine setCoefficientsSphericalComponentTransformation(this,stations)
    ! incoming
    type (component_transformation) :: this
    type (seismic_station), dimension(:) :: stations
    ! local
    double precision :: lat,lon,sinlat,sinlon,coslat,coslon
    integer :: istat,nstat
    !
    nstat = size(stations)
    if(nstat == 0) return
    !
    if(.csys.stations(1) /= 'S') then
       call setCoefficientsCartesianComponentTransformation(this)
       return
    endif
    !
    if(this%nstat /= 0) call deallocateCoefficientsComponentTransformation(this)
    this%csys = 'S'
    allocate(this%C2XYZ(3,number_of_components,nstat),this%XYZ2C(3,number_of_components,nstat),this%staname(nstat))
    this%C2XYZ = 0.d0 ; this%XYZ2C = 0.d0
    !
    do istat = 1,nstat
       this%staname(istat) = .staname.stations(istat)
       lat = dble(.lat.stations(istat))*mc_deg2radd
       lon = dble(.lon.stations(istat))*mc_deg2radd
       !
       sinlat = dsin(lat); sinlon = dsin(lon)
       coslat = dcos(lat); coslon = dcos(lon)
       !
       ! Set coefficients (i.e. entries of transformation matrices) to transform a specific (set of) component(s) to CX,CY,CZ.
       ! Any specific collection of columns of array this%C2XYZ forms a transformation matrix from those corresponding 
       ! components to CX,CY,CZ
       this%C2XYZ(:,indxComponent('CX'),istat) = (/ 1.d0,0.d0,0.d0 /)
       this%C2XYZ(:,indxComponent('CY'),istat) = (/ 0.d0,1.d0,0.d0 /)
       this%C2XYZ(:,indxComponent('CZ'),istat) = (/ 0.d0,0.d0,1.d0 /)
       this%C2XYZ(:,indxComponent('N'),istat) = (/ -sinlat*coslon,-sinlat*sinlon,coslat /)
       this%C2XYZ(:,indxComponent('S'),istat) = (/ sinlat*coslon,sinlat*sinlon,-coslat /)
       this%C2XYZ(:,indxComponent('E'),istat) = (/ -sinlon,coslon,0.d0 /)
       this%C2XYZ(:,indxComponent('W'),istat) = (/ sinlon,-coslon,0.d0 /)
       this%C2XYZ(:,indxComponent('UP'),istat) = (/ coslat*coslon,coslat*sinlon,sinlat /)
       this%C2XYZ(:,indxComponent('DOWN'),istat) = (/- coslat*coslon,-coslat*sinlon,-sinlat /)
       ! EXAMPLE
       !    this%C2XYZ(:,indxComponent('HX'),istat) = (/ 2.d0,0.d0,0.d0 /)
       !
       ! Set coefficients (i.e. entries of transformation matrices) to transform CX,CY,CZ to a specific (set of) component(s).
       ! Any specific collection of columns of array this%XYZ2C forms THE TRANSPOSE of a transformation matrix from CX,CY,CZ to 
       ! those corresponding components.
       ! Store in transposed from (i.e. CX,CY,cZ in first dimension), as for transformation CX,CY,CZ, -> C you always need all values
       ! of CX,CY,CZ for a specific component C (more efficient when accessing the array this%XYZ2C, compare transpose() statement 
       ! inside matmul in getCoefficients.. routines)
       ! IN CASE OF ORTHOGONAL TRANSFORMATIONS, the values of a specific column of arrays this%C2XYZ and this%XYZ2C will
       ! coincide, as in this case the inverse of the transformation matrix is its transposed!
       ! However, the coefficients this%XYZ2C (which correspond to the inverse of the transformation by this%C2XYZ) 
       ! are kept seperately from this%C2XYZ in order to also support non-orthogonal basis transformations in this module
       ! (compare example component 'HX' = CX/2 ).
       this%XYZ2C(:,indxComponent('CX'),istat) = (/ 1.d0,0.d0,0.d0 /)
       this%XYZ2C(:,indxComponent('CY'),istat) = (/ 0.d0,1.d0,0.d0 /)
       this%XYZ2C(:,indxComponent('CZ'),istat) = (/ 0.d0,0.d0,1.d0 /)
       this%XYZ2C(:,indxComponent('N'),istat) = this%C2XYZ(:,indxComponent('N'),istat) !(/ -sinlat*coslon,-sinlat*sinlon,coslat /)
       this%XYZ2C(:,indxComponent('S'),istat) = this%C2XYZ(:,indxComponent('S'),istat) !(/ sinlat*coslon,sinlat*sinlon,-coslat /)
       this%XYZ2C(:,indxComponent('E'),istat) = this%C2XYZ(:,indxComponent('E'),istat) !(/ -sinlon,coslon,0.d0 /)
       this%XYZ2C(:,indxComponent('W'),istat) = this%C2XYZ(:,indxComponent('W'),istat) !(/ sinlon,-coslon,0.d0 /)
       this%XYZ2C(:,indxComponent('UP'),istat) = this%C2XYZ(:,indxComponent('UP'),istat) !(/ coslat*coslon,coslat*sinlon,sinlat /)
       this%XYZ2C(:,indxComponent('DOWN'),istat) = this%C2XYZ(:,indxComponent('DOWN'),istat) !(/ -coslat*coslon,-coslat*sinlon,-sinlat /)
       ! EXAMPLE
       !    this%XYZ2C(:,indxComponent('HX'),istat) = (/ 0.5d0,0.d0,0.d0 /)
    enddo ! istat
    this%nstat = nstat
  end subroutine setCoefficientsSphericalComponentTransformation
!------------------------------------------------------------------------
!> \brief get coefficients to express one component in 3 others, Cartesian case only
!! \details the user has to take care of choosing sensible (numbers of) components
!!  to transform!
!! \param this component transformation object
!! \param C_in array n_in of components in which C_out components shall be expressed
!! \param C_out array of n_out component which shall be expressed in C_out components
!! \param coef n_out-by-n_in coefficient matrix such that C_out = coef * C_in
!! \return coefficient matrix to transform components C_in to components C_out
!
  function getCoefficientsCartesianComponentTransformation(this,C_in,C_out) result(coef)
    ! incoming
    type (component_transformation) :: this
    character(len=*), dimension(:) :: C_in,C_out
    ! outgoing
    double precision, dimension(:,:), pointer :: coef
    ! local
    integer :: n_in,n_out,i
    integer, dimension(:), allocatable :: indx_in,indx_out
    !
    coef => null()
    if(this%nstat == 0) return
    !
    ! if coefficients were actually defined for spherical case, return
    if(this%csys == 'S') return
    !
    ! if there are no incoming or outgoing coefficients, return
    n_in = size(C_in); n_out = size(C_out)
    if(n_in == 0 .or. n_out == 0) return
    !
    allocate(indx_in(n_in),indx_out(n_out))
    do i = 1,n_in
       indx_in(i) = indxComponent(C_in(i))
    enddo ! i
    do i = 1,n_out
       indx_out(i) = indxComponent(C_out(i))
    enddo ! i
    ! if there are any invalid incoming or outgoing coefficients, return
    if(any(indx_in<1) .or. any(indx_out<1)) then
       deallocate(indx_in,indx_out)
       return
    endif
    !
    allocate(coef(n_out,n_in))
    coef = matmul(transpose(this%XYZ2C(:,indx_out,1)),this%C2XYZ(:,indx_in,1))
    !
    deallocate(indx_in,indx_out)
  end function getCoefficientsCartesianComponentTransformation
!------------------------------------------------------------------------
!> \brief get coefficients of istat'th station to express one set of component in some others
!! \details the user has to take care of choosing sensible (numbers of) components
!!  to transform! Here, Cartesian case is also supported (then istat is ignored)
!! \param this component transformation object
!! \param C_in array n_in of components in which C_out components shall be expressed
!! \param C_out array of n_out component which shall be expressed in C_out components
!! \param istat index of requested station (as transformation coefficients might be different for each station)
!! \param coef n_out-by-n_in coefficient matrix such that C_out = coef * C_in
!! \return coefficient matrix to transform components C_in to components C_out
!
  function getCoefficientsByIndexComponentTransformation(this,C_in,C_out,istat) result(coef)
    ! incoming
    type (component_transformation) :: this
    character(len=*), dimension(:) :: C_in,C_out
    integer :: istat
    ! outgoing
    double precision, dimension(:,:), pointer :: coef
    ! local
    integer :: n_in,n_out,i
    integer, dimension(:), allocatable :: indx_in,indx_out
    !
    coef => null()
    if(this%nstat == 0) return
    !
    ! if coefficients were defined for Cartesian case, call respective function
    if(this%csys == 'C') then
       coef => getCoefficientsCartesianComponentTransformation(this,C_in,C_out)
       return
    elseif(this%csys == 'S') then
       ! if station index is out of range, return
       if(istat>this%nstat .or. istat<1) return
    else
       return
    endif
    !
    ! if there are no incoming or outgoing coefficients, return
    n_in = size(C_in); n_out = size(C_out)
    if(n_in == 0 .or. n_out == 0) return
    !
    allocate(indx_in(n_in),indx_out(n_out))
    do i = 1,n_in
       indx_in(i) = indxComponent(C_in(i))
    enddo ! i
    do i = 1,n_out
       indx_out(i) = indxComponent(C_out(i))
    enddo ! i
    ! it there are any invalid incoming or outgoing coefficients, return
    if(any(indx_in<1) .or. any(indx_out<1)) then
       deallocate(indx_in,indx_out)
       return
    endif
    !
    allocate(coef(n_out,n_in))
    coef = matmul(transpose(this%XYZ2C(:,indx_out,istat)),this%C2XYZ(:,indx_in,istat))
    !
    deallocate(indx_in,indx_out)
  end function getCoefficientsByIndexComponentTransformation
!------------------------------------------------------------------------
!> \brief get coefficients of station "staname" to express one set of component in some others
!! \details the user has to take care of choosing sensible (numbers of) components
!!  to transform! Here, Cartesian case is also supported (then istat is ignored)
!! \param this component transformation object
!! \param C_in array n_in of components in which C_out components shall be expressed
!! \param C_out array of n_out component which shall be expressed in C_out components
!! \param staname station name of requested station (as transformation coefficients might be different for each station)
!! \param coef n_out-by-n_in coefficient matrix such that C_out = coef * C_in
!! \return coefficient matrix to transform components C_in to components C_out
!
  function getCoefficientsByNameComponentTransformation(this,C_in,C_out,staname) result(coef)
    ! incoming
    type (component_transformation) :: this
    character(len=*), dimension(:) :: C_in,C_out
    character(len=*) :: staname
    ! outgoing
    double precision, dimension(:,:), pointer :: coef
    ! local
    integer :: istat
!
    coef => null()
    if(this%nstat == 0) return
!
    ! if coefficients were defined for Cartesian case, call respective function
    if(this%csys == 'C') then
       coef => getCoefficientsCartesianComponentTransformation(this,C_in,C_out)
       return
    elseif(this%csys == 'S') then
       ! for incoming station name, find respective index and call getCoefficientsByIndexComponentTransformation
       do istat = 1,this%nstat
          if(this%staname(istat)==staname) then
             coef => getCoefficientsByIndexComponentTransformation(this,C_in,C_out,istat)
             return
          end if
       end do ! istat
       ! if staname was not found in array this%staname, coef just points to null() on exit
    end if
!
  end function getCoefficientsByNameComponentTransformation
!------------------------------------------------------------------------
!> \brief check if component character is supported by this module
!! \param C component character
!! \param l logical value if C is supported by this module
!! \return logical value if C is supported by this module
!
  function validComponent(C) result(l)
    character(len=*) :: C
    logical :: l
    l = indxComponent(C)>0
  end function validComponent
!------------------------------------------------------------------------
!> \brief check if all of the component names are supported by this module
!! \param C vector of component names
!! \param i_invalid optional integer, returns the smmallest index for which a component in C is invalid
!! \param l logical value if all entries in C are supported by this module
!! \return logical value if all entries in C are supported by this module
!
  function allValidComponents(C,i_invalid) result(l)
    character(len=*), dimension(:) :: C
    integer, optional :: i_invalid
    logical :: l
    ! local
    integer :: icomp
!
    ! initiate negative result
    l = .false.
    if(present(i_invalid)) i_invalid = -1
!
    ! technically, if no component is present, then all present components are valid
    ! however, this might not be intended by the program calling this routine, so return .false.
    if(size(C)<=0) return
!
    do icomp=1,size(C)
        if(indxComponent(C(icomp))==-1) then
           if(present(i_invalid)) i_invalid = icomp
           return
        end if ! indxComponent(C(icomp))==-1 , i.e. C(icomp) is not valid
     end do ! icomp
!
     ! if function comes here, all components are valid, so return .true.
     l = .true.
  end function allValidComponents
!------------------------------------------------------------------------
!> \brief check if all of the incoming component names are contained in given components list
!! \param C_in vector of component names to be checked
!! \param C_list vector of components
!! \param indx_map optional mapping: indx_map(indx_in_C_in) = indx_in_C_list
!! \param l logical value if all entries in C_in are contained in given list of components C_list
!! \return logical value if all entries in C_in are contained in given list of components C_list
!
  function allComponentsContained(C_in,C_list,indx_map) result(l)
    character(len=*), dimension(:) :: C_in,C_list
    integer, dimension(:), pointer, optional :: indx_map
    logical :: l
    ! local
    integer :: icomp,jcomp
!
    ! initiate negative result
    l = .false.
    if(present(indx_map)) nullify(indx_map)
!
    ! technically, if no component is present, then all present components are contained in the test list.
    ! however, this might not be intended by the program calling this routine, so return .false. .
    ! also return false, if the test list holds no components
    if(size(C_in)<=0 .or. size(C_list)<=0) return
!
    do icomp=1,size(C_in)
        if(.not.any(C_list==C_in(icomp))) return
     end do ! icomp
!
     ! if function comes here, all components in C_in are contained in given list of components C_list, so return .true.
     l = .true.
!
     if(present(indx_map)) then
        allocate(indx_map(size(C_in)))
        do icomp = 1,size(C_in)
           do jcomp = 1,size(C_list)
              if(C_in(icomp)==C_list(jcomp)) then
                 indx_map(icomp) = jcomp
                 exit
              end if
           end do ! jcomp
        end do ! icomp
     end if
  end function allComponentsContained
!------------------------------------------------------------------------
! TODO 
! routines to add here:
!   - check if components_in or components_out are independent (check that transformation matrices have full rank) (unnecessary?? user should know what he is doing!)
!------------------------------------------------------------------------
!> \brief deallocate transformation arrays
!! \param this component transformation object
!
  subroutine deallocateCoefficientsComponentTransformation(this)
    type (component_transformation) :: this
    this%csys = ''
    if(associated(this%C2XYZ)) deallocate(this%C2XYZ)
    if(associated(this%XYZ2C)) deallocate(this%XYZ2C)
    if(associated(this%staname)) deallocate(this%staname)
    this%nstat = 0
  end subroutine deallocateCoefficientsComponentTransformation
!
end module componentTransformation
