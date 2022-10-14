!> \file lintran_interface_module.f90
!! Module file for blueprints of all radiative transfer and perturbation routines.

!> Module for blueprints of all radiative transfer and perturbation routines.
!!
!! This module only contains the inputs and the outputs of all transfer and perturbation
!! routines. These are used to create procedure pointers to routines that look like these.
!! The selection only determines the number of Stokes parameters. Doing this on compile
!! time increases the speed.
module lintran_interface_module

  implicit none

  !> This is an abstract interface.
  !!
  !! This just contains function blueprints without implementation. Compare this to
  !! an h-file in C++. It is needed that the types and shapes of the arguments are
  !! known before the implementation is known. This allows functions to be defined
  !! before it is decided what it exactly does. In the interfaces, dimension sizes
  !! zero are given for any dimension. This is because Fortran does not recognize
  !! dimension sizes are a type specification, but it does recognize the rank as
  !! a type specification. The dimension size is not important, so a zero is okay.
  !! A colon as assumed dimension does not work in this interface. In a test program,
  !! it resulted in uninitialized values.
  interface

    ! Unperturbed radiative transfer routines. {{{

    !> Blueprint for direct_1_elem.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine direct_1_elem_interface(nlay,exist_therm,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity. As this is direct single-scattering, this can only be if the viewing Stokes parameter is the one in which thermal emission takes place.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    end subroutine direct_1_elem_interface ! }}}

    !> Blueprint for direct_2.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine direct_2_interface(nstrhalf,nlay,exist_therm,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    end subroutine direct_2_interface ! }}}

    !> Blueprint for transmit.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine transmit_interface(nstrhalf,nlay,flag_reverse,intens_in,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      logical, intent(in) :: flag_reverse !< True to reverse transmission directions, useful when the input is a response function.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: intens_in !< Input intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(out) :: res !< The resulting field.

    end subroutine transmit_interface ! }}}

    !> Blueprint for transmit_sp.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine transmit_sp_interface(nstrhalf,nlay,nsplit_total,flag_reverse,intens_in,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: flag_reverse !< True to reverse transmission directions, useful when the input is a response function.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_in !< Input intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    end subroutine transmit_sp_interface ! }}}

    !> Blueprint for source_response_1.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine source_response_1_interface(nstrhalf,nlay,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_emi !< Flag for if fluorescent emission exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(out) :: res !< The resulting field.

    end subroutine source_response_1_interface ! }}}

    !> Blueprint for source_response_1_sp.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine source_response_1_sp_interface(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_emi !< Flag for if fluorescent emission exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    end subroutine source_response_1_sp_interface ! }}}

    !> Blueprint for source_response_1_int.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine source_response_1_int_interface(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_emi !< Flag for if fluorescent emission exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    end subroutine source_response_1_int_interface ! }}}

    !> Blueprint for source_response_2_sp_addition.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine source_response_2_sp_addition_interface(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_therm,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(inout) :: res !< The resulting field added to the original values.

    end subroutine source_response_2_sp_addition_interface ! }}}

    !> Blueprint for scatter_1_sp.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine scatter_1_sp_interface(nstrhalf,nlay,nsplit_total,exist_bdrf,intens_in,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: intens_in !< Input intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out), target :: res !< The resulting field.

    end subroutine scatter_1_sp_interface ! }}}

    !> Blueprint for scatter_1_int.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine scatter_1_int_interface(nstrhalf,nlay,nsplit_total,exist_bdrf,intens_in,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_in !< Input intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    end subroutine scatter_1_int_interface ! }}}

    !> Blueprint for innerproduct.
    !!
    !! The implementation code is in lintran_transfer_implementation.f90.
    subroutine innerproduct_interface(nstrhalf,nsplit_total,fld1,fld2,res) ! {{{

      use lintran_constants_module, only: nst_fake

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld1 !< One field for the inner product.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld2 !< The other field for the inner product.
      real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    end subroutine innerproduct_interface ! }}}

    ! }}}

    ! Perturbation routines for the scattering. {{{
   
    !> Blueprint for perturbation_scat_direct_1_elem.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_direct_1_elem_interface(nlay,nlay_deriv,ilay_deriv,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay), intent(inout) :: res !< Resulting derivative with respect to the single-scattering C-parameter for the relevant Stokes parameters. The result of this routine will be added to the input.

    end subroutine perturbation_scat_direct_1_elem_interface ! }}}

    !> Blueprint for perturbation_scat_direct_2.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_direct_2_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res_fwd_0,res_bck_0,res_fwd_v,res_bck_v) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_0 !< Resulting derivative with respect to C-parameters for forward scattering from the solar direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_0 !< Resulting derivative with respect to C-parameters for backward scattering from the solar direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_v !< Resulting derivative with respect to C-parameters for forward scattering to the viewing direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_v !< Resulting derivative with respect to C-parameters for backward scattering to the viewing direction. The result of this routine will be added to the input.

    end subroutine perturbation_scat_direct_2_interface ! }}}

    !> Blueprint for perturbation_scat_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_source_response_1_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    end subroutine perturbation_scat_source_response_1_interface ! }}}

    !> Blueprint for perturbation_scat_source_response_1_sp.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_source_response_1_sp_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    end subroutine perturbation_scat_source_response_1_sp_interface ! }}}

    !> Blueprint for perturbation_scat_source_response_1_int.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_source_response_1_int_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    end subroutine perturbation_scat_source_response_1_int_interface ! }}}

    !> Blueprint for perturbation_scat_source_response_2_int.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_source_response_2_int_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_therm,fld,tlc,res_fwd,res_bck,res_fwd_x,res_bck_x) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
      real, dimension(nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    end subroutine perturbation_scat_source_response_2_int_interface ! }}}

    !> Blueprint for perturbation_scat_scatter_1_int.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_scatter_1_int_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res_fwd,res_bck) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.

    end subroutine perturbation_scat_scatter_1_int_interface ! }}}

    !> Blueprint for perturbation_scat_zip.
    !!
    !! The implementation code is in lintran_perturbation_scat_implementation.f90.
    subroutine perturbation_scat_zip_interface(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,nlay_deriv_base,ilay_deriv_base,intens,intens_adj,tlc,chain_c,res_fwd,res_bck,chain_base,res_base) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv_c !< Number of layers for which derivatives with respect to C are calculated.
      integer, dimension(nlay_deriv_c), intent(in) :: ilay_deriv_c !< Indices of layers for which derivatives with respect to C are calculated.
      integer, intent(in) :: nlay_deriv_base !< Number of layers for which derivatives with respect to basic properties (including optical depth) are calculated.
      integer, dimension(nlay_deriv_base), intent(in) :: ilay_deriv_base !< Indices of layers for which derivatives with respect to basic properties (including optical depth) are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay), intent(in) :: chain_c !< Chain rule factors for derivatives with respect to C.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.
      real, dimension(nlay), intent(in) :: chain_base !< Chain rule factors for derivatives with respect to the optical depth, the basic property that is relevant.
      real, dimension(nlay_deriv_base), intent(inout) :: res_base !< Resulting derivative with respect to the optical depth, the basic property that is relevant. The result of this routine will be added to the input.

    end subroutine perturbation_scat_zip_interface ! }}}

    ! }}}

    ! Perturbation routines for the optical depth. {{{

    !> Blueprint for perturbation_tau_direct_1_elem.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_direct_1_elem_interface(nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth, for the relevant viewing Stokes parameter. The result of this routine will be added to the input.

    end subroutine perturbation_tau_direct_1_elem_interface ! }}}

    !> Blueprint for perturbation_tau_direct_2_elem.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_direct_2_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_direct_2_interface ! }}}

    !> Blueprint for perturbation_tau_mutualtransmit.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_mutualtransmit_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: source !< Intensity that is at the interface of the perturbed layer.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: destination !< Field with which the inner product is taken. This field includes an inverse transmission, so that the perturbation can be done on only one layer.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_mutualtransmit_interface ! }}}

    !> Blueprint for perturbation_tau_mutualtransmit_sp.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_mutualtransmit_sp_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Intensity that is at the interface of a perturbed sublayer.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken. This field includes an inverse transmission, so that the perturbation can be done on only one sublayer.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_mutualtransmit_sp_interface ! }}}

    !> Blueprint for perturbation_tau_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_source_response_1_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_source_response_1_interface ! }}}

    !> Blueprint for perturbation_tau_source_response_1_sp.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_source_response_1_sp_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_source_response_1_sp_interface ! }}}

    !> Blueprint for perturbation_tau_source_response_1_int.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_source_response_1_int_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_source_response_1_int_interface ! }}}

    !> Blueprint for perturbation_tau_source_response_2_int.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_source_response_2_int_interface(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_therm,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_source_response_2_int_interface ! }}}

    !> Blueprint for perturbation_tau_scatter_1_int.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_tau_scatter_1_int_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_tau_scatter_1_int_interface ! }}}

    !> Blueprint for perturbation_abs_zip.
    !!
    !! The implementation code is in lintran_perturbation_tau_implementation.f90.
    subroutine perturbation_abs_zip_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,intens,intens_adj,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_abs_zip_interface ! }}}

    ! }}}

    ! Perturbation routines for the BDRF. {{{
   
    !> Blueprint for perturbation_bdrf_direct_1_elem.
    !!
    !! The implementation code is in lintran_perturbation_bdrf_implementation.f90.
    subroutine perturbation_bdrf_direct_1_elem_interface(nlay,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, intent(inout) :: res !< Resulting derivative with respect to the surface reflection at single-scattering geometry, for the relevant Stokes parameters. The result of this routine will be added to the input.

    end subroutine perturbation_bdrf_direct_1_elem_interface ! }}}

    !> Blueprint for perturbation_bdrf_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_bdrf_implementation.f90.
    subroutine perturbation_bdrf_source_response_1_interface(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
      logical, intent(in) :: flag_source !< True for source. False for response.
      logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF in which the solar or viewing direction is involved. The result of this routine will be added to the input.

    end subroutine perturbation_bdrf_source_response_1_interface ! }}}

    !> Blueprint for perturbation_bdrf_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_bdrf_implementation.f90.
    subroutine perturbation_bdrf_scatter_1_interface(nstrhalf,nsplit_total,source,destination,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF among internal streams. The result of this routine will be added to the input.

    end subroutine perturbation_bdrf_scatter_1_interface ! }}}

    !> Blueprint for perturbation_bdrf_zip.
    !!
    !! The implementation code is in lintran_perturbation_bdrf_implementation.f90.
    subroutine perturbation_bdrf_zip_interface(nstrhalf,nsplit_total,intens,intens_adj,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nst_fake,nstrhalf,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF among internal streams. The result of this routine will be added to the input.

    end subroutine perturbation_bdrf_zip_interface ! }}}

    ! }}}

    ! Perturbation routines for the fluorescent emissivity. {{{

    !> Blueprint for perturbation_emi_direct_1_elem.
    !!
    !! The implementation code is in lintran_perturbation_emi_implementation.f90.
    subroutine perturbation_emi_direct_1_elem_interface(nlay,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, intent(inout) :: res !< Resulting derivative with respect to the surface emissivity directed towards the instrument, for the relevant viewing Stokes parameter. The result of this routine will be added to the input.

    end subroutine perturbation_emi_direct_1_elem_interface ! }}}

    !> Blueprint for perturbation_emi_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_emi_implementation.f90.
    subroutine perturbation_emi_source_response_1_interface(nstrhalf,nsplit_total,flag_source,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
      logical, intent(in) :: flag_source !< True for source. False for response.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nst_fake,nstrhalf), intent(out) :: res !< Resulting derivative with respect to surface emissivity. The result of this routine will be added to the input.

    end subroutine perturbation_emi_source_response_1_interface ! }}}

    ! }}}

    ! Perturbation routines for thermal emission. {{{

    !> Blueprint for perturbation_therm_direct_1_elem.
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_direct_1_elem_interface(nlay_deriv,ilay_deriv,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the thermal emissivity. The result of this routine will be added to the input.

    end subroutine perturbation_therm_direct_1_elem_interface ! }}}

    !> Blueprint for perturbation_therm_direct_2
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_direct_2_interface(nlay_deriv,ilay_deriv,tlc,res) ! {{{

      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_therm_direct_2_interface ! }}}

    !> Blueprint for perturbation_therm_source_response_1.
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_source_response_1_interface(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nlay !< Number of atmospheric layers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      real, dimension(nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_therm_source_response_1_interface ! }}}

    !> Blueprint for perturbation_therm_source_response_1_sp.
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_source_response_1_sp_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_therm_source_response_1_sp_interface ! }}}

    !> Blueprint for perturbation_therm_source_response_1_int.
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_source_response_1_int_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_therm_source_response_1_int_interface ! }}}

    !> Blueprint for perturbation_therm_source_response_2_int.
    !!
    !! The implementation code is in lintran_perturbation_therm_implementation.f90.
    subroutine perturbation_therm_source_response_2_int_interface(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_tlc_class

      implicit none

      ! Input and output.
      integer, intent(in) :: nstrhalf !< Half the number of streams.
      integer, intent(in) :: nsplit_total !< Total number of sublayers.
      integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
      integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
      logical, intent(in) :: flag_source !< True for source. False for response.
      real, dimension(nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
      type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
      real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    end subroutine perturbation_therm_source_response_2_int_interface ! }}}

    ! }}}

  end interface

end module lintran_interface_module
