!------------------------------------------------------------------------------
!> Settings and parameters for the LINTRAN v2.0 radiative transfer software
!------------------------------------------------------------------------------
module lintran_header
   implicit none

! Normal settings
real, parameter :: taua_split = 0.4
real, parameter :: tautot_max = 30.d0 !15.0
real, parameter :: fourier_tolerance = 1.d-4
real, parameter :: gs_tolerance = 1.d-4
integer, parameter ::gs_maxiter = 25            ! max numbers of gauss seidel iterations
 
!*** Settings for cloudy scenes
real, parameter :: taus_split_cld = 0.1
logical, parameter ::split_double_cld = .true.   ! do double scattering analytical
integer, parameter ::interpolation_cld = 1          ! Interpolation method for intensity in layer
integer, parameter ::solver_cld = 2                 ! 1=Gauss-Seidel, 2=LU decomposition

!*** Settings for clear-sky scenes
real, parameter :: taus_split_clr = 0.02
logical, parameter ::split_double_clr = .false.
integer, parameter ::interpolation_clr = 0          ! Interpolation method for intensity in layer
integer, parameter ::solver_clr = 1                 ! 1=Gauss-Seidel, 2=LU decomposition
                                                ! 0=averaging, 1=linear
!------------------------------------------------------------------------------

! Paranoid settings
! real, parameter :: taua_split = 0.4
! real, parameter :: tautot_max = 15.0
! real, parameter :: fourier_tolerance = 0.0
! real, parameter :: gs_tolerance = 1.d-4
! integer, parameter ::gs_maxiter = 25            ! max numbers of gauss seidel iterations
!  
! !*** Settings for cloudy scenes
! real, parameter :: taus_split_cld = 0.01
! logical, parameter ::split_double_cld = .true.   ! do double scattering analytical
! integer, parameter ::interpolation_cld = 1          ! Interpolation method for intensity in layer
! integer, parameter ::solver_cld = 2                 ! 1=Gauss-Seidel, 2=LU decomposition
! 
! !*** Settings for clear-sky scenes
! real, parameter :: taus_split_clr = 0.01
! logical, parameter ::split_double_clr = .true.
! integer, parameter ::interpolation_clr = 1          ! Interpolation method for intensity in layer
! integer, parameter ::solver_clr = 2                 ! 1=Gauss-Seidel, 2=LU decomposition
!                                                 ! 0=averaging, 1=linear
! !------------------------------------------------------------------------------

end module lintran_header
