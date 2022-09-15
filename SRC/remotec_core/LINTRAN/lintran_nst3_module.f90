!> \file lintran_nst3_module.f90
!! Module file the reduced vector radiative transport.

!> Module the reduced vector radiative transport.
!!
!! This module defines all the radiative transfer functions and their perturbations to the
!! reduced vector version. That is done just be defining the number of Stokes parameters to three
!! inside the module and include the implementations of the functions where there will be
!! references to the number of Stokes parameters. The implementations are programmed as generic
!! functions, so they have generic names and just assume that someone defines the number of
!! Stokes parameters. The implementation files also contain the routine names, otherwise we
!! could not include a whole file. This implies that when importing a routine, it should be
!! renamed to specify for which number of Stokes parameters the function has been implemented.
module lintran_nst3_module

  implicit none

  integer, parameter :: nst = 3 !< Number of Stokes parameters that is applied if any routine is imported from this module. For the rest, all routines do the same.

contains

  ! The implementations are separated in different parts, the unperturbed radiative transport
  ! and a perturbation for each type of parameter that is perturbed. The modular structure
  ! is a bit damaged, becuase all different parts of the radiative transfer end up in the
  ! same module, but the files are still separated in what should have been a module. If we
  ! keps these modules separate, there would be an enormous number of files and in each
  ! files all these comment lines would have been repeated.
  include "lintran_transfer_implementation.f90"
  include "lintran_perturbation_scat_implementation.f90"
  include "lintran_perturbation_tau_implementation.f90"
  include "lintran_perturbation_bdrf_implementation.f90"
  include "lintran_perturbation_emi_implementation.f90"
  include "lintran_perturbation_therm_implementation.f90"

end module lintran_nst3_module
