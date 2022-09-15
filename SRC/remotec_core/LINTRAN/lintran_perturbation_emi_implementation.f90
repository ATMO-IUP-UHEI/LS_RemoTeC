!> \file lintran_perturbation_emi_implementation.f90
!! Implementation file for derivatives with respect to the fluorescent emissivity.
!! This file has a part of the contents from lintran_transfer_implementation, but here
!! derivatives are taken with respect to surface emissivity. Inner products are directly
!! calculated.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Calculates the derivative of direct single-scatterig with respect to surface emissivity.
  !!
  !! This routine calculates the derivative of direct_1_elem from lintran_transfer_module
  !! with respect to surface emissivity.
  subroutine perturbation_emi_direct_1_elem(nlay,tlc,res) ! {{{

    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, intent(inout) :: res !< Resulting derivative with respect to the surface emissivity directed towards the instrument, for the relevant viewing Stokes parameter. The result of this routine will be added to the input.

    ! Emission is only relevant in the Stokes parameter equal to the viewing Stokes parameter,
    ! because this is single-scattering. And the fluorescent emission event is counted as a
    ! scatter event. Therefore, fluorescently emitted light can only change Stokes parameter
    ! from the second scatter event.

    ! TLC-parameters for fluorescent emission:
    ! T = T_tot_v at bottom
    ! L = 1
    ! C = Ce_v_elem

    ! Differentiating will turn C into her derivative.
    res = res + tlc%prefactor_dir * tlc%t_tot_v(nlay) * tlc%diff_c_emi_v_elem_emi

  end subroutine perturbation_emi_direct_1_elem ! }}}

  !> Calculates the inner products with the derivative of a first order source or response.
  !!
  !! This routine calculates the derivative of source_response_1, source_response_1_sp or
  !! source_response_1_int from lintran_transfer_module with respect to surface emissivity
  !! to internal streams. Directly, the inner product with an existing field is taken, so only
  !! the field at the bottom is relevant.
  subroutine perturbation_emi_source_response_1(nstrhalf,nsplit_total,flag_source,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
    logical, intent(in) :: flag_source !< True for source. False for response.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf), intent(out) :: res !< Resulting derivative with respect to surface emissivity. The result of this routine will be added to the input.

    ! Pointers to TLC-parameters.
    real, dimension(:), pointer :: prefactor
    real, dimension(:,:), pointer :: diff_c_emi_emi

    ! Direction index.
    integer :: idir_bck

    ! Iterator.
    integer :: istr ! Over streams.

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_bck = iup
      ! This is the normal case, otherwise the call is nonsense.
      diff_c_emi_emi => tlc%diff_c_emi_nrm_emi
    else
      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_bck = idn
      ! This is the adjoint case, otherwise the call is nonsense
      diff_c_emi_emi => tlc%diff_c_emi_adj_emi
    endif

    ! This is only relevant for adjoint response or normal source. Otherwise,
    ! this routine should not be called.

    ! TLC-parameters for surface-emission:
    ! T = 1
    ! L = 1
    ! C = Ce

    ! Differentiating turns C into her derivative.

    do istr = 1,nstrhalf
      res(:,istr) = res(:,istr) + prefactor(istr) * diff_c_emi_emi(:,istr) * fld(:,istr,idir_bck,nsplit_total)
    enddo

  end subroutine perturbation_emi_source_response_1 ! }}}
