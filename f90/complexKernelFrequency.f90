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
!> \brief Manage possible complex frequencies for wavefields, kernels, data and filters
!!
!! \details For specific (frequency domain) forward methods, this module manages the knowledge
!!  about possible imaginary parts of frequencies used in ASKI for wavefields and kernels, as
!!  well as (synthetic) data and filters. By default, frequencies in ASKI are real-valued 
!!  and defined by f = df*jf , with the frequency step df and frequency index jf being defined
!!  in the ASKI parameter files.
!!
!! \author Florian Schumacher
!! \date Feb 2016
!! 
!! \copyright &copy; GNU Public License
!
module complexKernelFrequency
!
  use mathConstants, only:mc_pi
!
implicit none
!
contains
!
!------------------------------------------------------------------------
!> \brief ask if a method uses complex kernel frequencies at all
!! \param method forward method used by ASKI
!! \param l logical indicating whether the forward method uses complex frequencies or not
!
  function methodHasComplexKernelFrequency(method) result(l)
    ! incoming
    character(len=*) :: method
    ! returning
    logical :: l
!
    select case (method)
    case('GEMINI')
       l = .true.
    case('SPECFEM3D')
       l = .false.
    case('NEXD')
       l = .false.
    case default
       ! If the method is invalid/not defined here, do not use complex frequencies.
       ! In this case, the default interpretation of real-valued frequencies should
       ! be used everywhere.
       l = .false.
    end select ! method
  end function methodHasComplexKernelFrequency
!
!------------------------------------------------------------------------
!> \brief get the actual complex frequency dependent on real frequency step and index
!! \param method forward method used by ASKI
!! \param df real-valued frequency step as defined by ASKI main parameter file
!! \param jf frequency index of the requested frequency
!! \param f resulting complex-valued frequency
!
  function getComplexKernelFrequency(method,df,jf) result(f)
    ! incoming
    character(len=*) :: method
    real :: df
    integer :: jf
    ! returning
    complex :: f
!
    select case (method)
    case('GEMINI')
       ! in case of Gemini, there is a constant imaginary part of  -5*df/2pi
       f = cmplx( df*jf , -2.5*df/mc_pi )
    case('SPECFEM3D','NEXD')
       ! SPECFEM3D and NEXD use standard real-valued frequencies df*jf
       f = cmplx( df*jf )
    case default
       ! If the method is invalid/not defined here, do not use complex frequencies.
       ! In this case, the default interpretation of real-valued frequencies should
       ! be used everywhere.
       f = cmplx( df*jf )
    end select
  end function getComplexKernelFrequency
!
end module complexKernelFrequency
