!> \file lintran_constants_module.f90
!! Module file for mathematical, numeric and organizational constants.

!> Module for mathematical, numeric and organizational constants.
!!
!! This module only contants member values with the 'parameter' keyword.
!! Therefore, they are no variables. The values are either mathematical
!! constants, values that assert numeric stability or array indices that are
!! preferrably referred by their name instead of their numbers.
module lintran_constants_module

  implicit none

  ! Error handling.
  ! The status number is a combination of zero or more warning flags. Errors are also considered
  ! as warning, in the sense: 'Warning: your program has died. Buy a new computer'. Warning flags
  ! co-exist in a bitwise manner. For now, we reserve the lowest hexadecimal digit for deadly
  ! warning, also referred to as errors. Three others are for warnings. This is just a small number
  ! that fits in an integer no matter whether it is signed or not, so no compiler flags have
  ! to be set. This restricts Lintran to 4 deadly warning and 12 normal warnings.
  integer, parameter :: lintran_errormask = Z"000F" !< Mask to ignore warnings and keep the errors, using bitwise and operator. This constant is not used in Lintran, but can be imported by the user.
  integer, parameter :: errorflag_allocation = Z"0001" !< Error flag for a memory error on allocating a heap array.
  integer, parameter :: errorflag_deallocation = Z"0002" !< Error flag for a memory error on deallocating a heap array.
  integer, parameter :: errorflag_stokes = Z"0004" !< Error flag for an non-supported number of Stokes parameters.
  integer, parameter :: warningflag_oddstreams = Z"0010" !< Warning code for an odd number of streams. Continuing with even number, rounding down.
  integer, parameter :: warningflag_dgnonconvergence = Z"0020" !< Warning flag for a non-convergence in calculating the double Gaussian stream grid.
  integer, parameter :: warningflag_instablegrid = Z"0040" !< Warning flag for a double Gaussian grid that will lead to numeric instability.
  integer, parameter :: warningflag_vulnerablegrid = Z"0080" !< Warning flag for a double Gaussian grid in which some solar or viewing geometry cannot be protected against numeric instability.
  integer, parameter :: warningflag_pseudosphericalerror = Z"0100" !< Warning code for setting pseudo-spherical geometry without altitude information. Continuing with plane-parallel geometry.
  integer, parameter :: warningflag_displaecedsurface = Z"0200" !< Altitude grid given with surface not at zero meters altitude.
  integer, parameter :: warningflag_unpreparedderivatives = Z"0400" !< Warning code for calculating derivatives when initialized for forward calculation only. Continuing without derivatives.
  integer, parameter :: warningflag_instableangle = Z"0800" !< Warning code for adapting an effective angle for numeric stability.
  integer, parameter :: warningflag_instablematrix = Z"1000" !< Warning code for numeric problems with the matrix.
  integer, parameter :: warningflag_gsnonconvergence = Z"2000" !< Warning code for non-convergence in Gauss-Seidel.

  ! Mathematics.
  real, parameter :: pi = 2.0*asin(1.0) !< Perimeter of a circle divided by the diameter, a mathematical constant.

  ! Assertion of stability.
  real, parameter :: minimum_mudiff = 1.0e-6 !< Minimum difference between different angle cosines for numeric stability.
  real, parameter :: lowlimit_tau = 1.0e-20 !< Minimum optical depth in one layer for numeric stability.
  real, parameter :: gs_protect_zero = 1.0e-20 !< Protection against division by a small number in convergence check for Gauss-Seidel.
  real, parameter :: protect_division = 1.0e-80 !< Protection against division by an almost-zero number.
  real, parameter :: protect_exponent = -log(protect_division) !< Highest exponent for which an exponential function can be properly evaluated.

  ! Convergence of double-Gaussian grid.
  integer, parameter :: dg_max_iterations = 1000 !< Maximum number of iterations for creating double-Gaussian stream grid.
  real, parameter :: dg_tolerance = 1.0e-14 !< Convergence tolerance for creating double-Gaussian grid.

  ! Stokes parameters.
  integer, parameter :: nst_max = 4 !< Maximum number of Stokes parameters.
  integer, parameter :: ist_i = 1 !< Index of Stokes parameter I.
  integer, parameter :: ist_q = 2 !< Index of Stokes parameter Q.
  integer, parameter :: ist_u = 3 !< Index of Stokes parameter U.
  integer, parameter :: ist_v = 4 !< Index of Stokes parameter V.
  integer, parameter :: ist_sun = ist_i !< Stokes parameter in which the sun shines.
  integer, parameter :: ist_therm = ist_i !< Stokes parameter in which thermal emission radiates.
  integer :: nst_fake !< An unknown number of Stokes parameters, to make the interface working under gfortran. This cannot be a parameter, though its value is completely irrelevant. It is only to convince gfortran that the interface does not yet know the number of Stokes parameters.

  ! What happens with the number of Stokes parameters.
  logical, dimension(nst_max), parameter :: legal_nst = (/.true.,.false.,.true.,.true./) !< Flags for each number of Stokes parameters from 1 to 4 whether or not it is possible to run Lintran with that.
  ! For the rest, the second value does not matter, because that is illegal. But, we want
  ! to index them with nst, not some chain of indices.
  logical, dimension(nst_max), parameter :: pol_nst = (/.false.,.false.,.true.,.true./) !< Flags for each number of Stokes parameters, including the illegal 2, whether or not we call this polarized (or vector) radiative transport.
  logical, dimension(nst_max), parameter :: fullpol_nst = (/.false.,.false.,.false.,.true./) !< Flags for each number of Stokes parameters, including the illegal 2, whether or not full polarization is included, where some symmetry relationships are broken.
  integer, dimension(nst_max), parameter :: nscat_nst = (/1,0,4,6/) !< Number of independent matrix elements for the phase matrix belonging to the number of Stokes parameter. A bullshit number is inserted for the illegal two.

  ! Independent matrix elements.
  integer, parameter :: iscat_a1 = 1 !< Diagonal matrix element at Stokes parameter I.
  integer, parameter :: iscat_b1 = 2 !< Off-diagonal matrix elements at Stokes parameters I and Q. Both have the same sign.
  integer, parameter :: iscat_a2 = 3 !< Diagonal matrix element at Stokes parameter Q.
  integer, parameter :: iscat_a3 = 4 !< Diagonal matrix element at Stokes parameter U.
  integer, parameter :: iscat_b2 = 5 !< Off-diagonal matrix elements at Stokes parameters U and V. Sign convention as in de Haan et al. (1987), positive from V to U and negative from U to V.
  integer, parameter :: iscat_a4 = 6 !< Diagonal matrix element at Stokes parameter V.

  real, dimension(nst_max), parameter :: flip_v = (/1.0,1.0,1.0,-1.0/) !< Array used to flip the sign of V-polarization, used for zipping a normal and adjoint field. Only relevant for V-polarization, but it is always created for four Stokes parameters, because it is a constant.

  ! Indices for down and up.
  integer, parameter :: idn = 0 !< Up/down index for downward radiation.
  integer, parameter :: iup = 1 !< Up/down index for upward radiation.

  ! Matrix solvers.
  integer, parameter :: solver_gauss_seidel = 1 !< Solver identifier for Gauss-Seidel.
  integer, parameter :: solver_lu_decomposition = 2 !< Solver identifier for LU-decomposition.

end module lintran_constants_module
