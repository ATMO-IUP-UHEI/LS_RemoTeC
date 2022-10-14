!> \file lintran_perturbation_bdrf_implementation.f90
!! Implementation file for derivatives with respect to the bi-directional reflectance function.
!! This file has a part of the contents from lintran_transfer_implementation, but here
!! derivatives are taken with respect to surface reflection (BDRF). Inner products are directly
!! calculated.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Calculates the derivative of direct single-scatterig with respect to surface reflection.
  !!
  !! This routine calculates the derivative of direct_1_elem from lintran_transfer_module
  !! with respect to surface reflection.
  subroutine perturbation_bdrf_direct_1_elem(nlay,tlc,res) ! {{{

    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, intent(inout) :: res !< Resulting derivative with respect to the surface reflection at single-scattering geometry, for the relevant Stokes parameters. The result of this routine will be added to the input.

    ! Like the perturbation of the single-scattering phase function, this is only relevant for
    ! the source Stokes parameter equal to that of the sun and the destination Stokes parameter
    ! equal to the viewing Stokes parameter. That is handled in the main lintran module.

    ! TLC-parameters for surface reflection:
    ! T = T_tot_0v at bottom
    ! L = 1
    ! C = Cb_srf_0v_elem

    ! Differentiating will turn C into her derivative.
    res = res + tlc%prefactor_dir * tlc%t_tot_0v(nlay) * tlc%diff_c_bck_srf_0v_elem_bdrf

  end subroutine perturbation_bdrf_direct_1_elem ! }}}

  !> Calculates the inner products with the derivative of a first order source or response.
  !!
  !! This routine calculates the derivative of source_response_1, source_response_1_sp or
  !! source_response_1_int from lintran_transfer_module with respect to surface reflection
  !! from the solar direction or to the viewing direction. Directly, the inner product
  !! with an existing field is taken, so only radiation at the bottom is relevant.
  subroutine perturbation_bdrf_source_response_1(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_v, fullpol_nst
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
    logical, intent(in) :: flag_source !< True for source. False for response.
    logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF in which the solar or viewing direction is involved. The result of this routine will be added to the input.

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: diff_c_bck_srf_x_bdrf

    ! Temporary result to enable Delta-flipping for the adjoint.
    real, dimension(nst,nstrhalf) :: temp_res

    ! Iterator.
    integer :: istr ! Over streams.

    ! Direction index.
    integer :: idir_bck

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_bck = iup
    else
      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_bck = idn
    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      diff_c_bck_srf_x_bdrf => tlc%diff_c_bck_srf_v_bdrf
    else
      ! Special direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      diff_c_bck_srf_x_bdrf => tlc%diff_c_bck_srf_0_bdrf
    endif

    ! TLC-parameters for surface reflection:
    ! T = T_tot_x at bottom
    ! L = 1
    ! C = Cb_srf_x

    ! Differentiating turns C into her derivative.

    do istr = 1,nstrhalf
      temp_res(:,istr) = prefactor(istr) * t_tot_x(nlay) * diff_c_bck_srf_x_bdrf(istr) * fld(:,istr,idir_bck,nsplit_total)
    enddo

    ! Apply Delta-flipping for the adjoint, only relevant for V-polarization.
    if (flag_adjoint .and. fullpol_nst(nst)) then

      temp_res(ist_v:ist_v,:) = -temp_res(ist_v:ist_v,:)

    endif

    ! Apply result.
    res = res + temp_res

  end subroutine perturbation_bdrf_source_response_1 ! }}}

  !> Calculates the inner products with the derivative of an internal scatter event.
  !!
  !! This routine calculates the derivative of scatter_1_sp or scatter_1_int (or scatter_1
  !! if it had existed) from lintran_transfer_module with respect to surface reflection.
  !! Directly, the inner product with an existing field is taken, so that derivatives for all
  !! layers can be handled at once. As scatter_1_sp and scatter_1_int also need an input field,
  !! this routine needs both an input field for the scatter event and a field with which the
  !! inner product is taken.
  subroutine perturbation_bdrf_scatter_1(nstrhalf,nsplit_total,source,destination,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF among internal streams. The result of this routine will be added to the input.

    ! This routine differs from the zip routine for the BDRF in the sense
    ! that both sides have the same character, normal or adjoint. Therefore,
    ! no zip prefactor is needed and the directions are up and down, where
    ! down is for the source side and up is for the destination side.

    ! Iterators.
    integer :: istr ! Over destination streams.
    integer :: istr_src ! Over source streams.
    integer :: ist ! Over destination Stokes parameters.

    ! TLC-parameters for surface reflection:
    ! T = 1
    ! L = 1
    ! C = Cb_srf

    ! Differentiating turns C into her derivative.

    ! The prefactor is prefactor_scat.

    ! Scalar      : None
    ! Array in s  : source
    ! Array in d  : destination
    ! Matrix(s,d) : prefactor_scat, diff_c

    do istr = 1,nstrhalf ! Destination stream.
      do istr_src = 1,nstrhalf ! Source stream.
        do ist = 1,nst ! Destination Stokes parameter.
          ! Array operation over source Stokes parameters.
          res(:,ist,istr_src,istr) = res(:,ist,istr_src,istr) + destination(ist,istr,iup,nsplit_total) * source(:,istr_src,idn,nsplit_total) * tlc%prefactor_scat(istr_src,istr) * tlc%diff_c_bck_srf_bdrf(istr_src,istr)
        enddo
      enddo
    enddo

  end subroutine perturbation_bdrf_scatter_1 ! }}}

  !> Takes the inner product of a normal and adjoint intensity field with a pertrubation of the radiative transfer operator.
  !!
  !! Perturbs the surface reflection in radiative transfer operator and performs the inner
  !! product of the adjoint field with the perturbation of the radiative transfer operator acting
  !! on the normal field. Corrections for the pseudo-forward method for the adjoint field are
  !! applied so that the inner product is correct.
  subroutine perturbation_bdrf_zip(nstrhalf,nsplit_total,intens,intens_adj,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, flip_v
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf), intent(inout) :: res !< Resulting derivative with respect to BDRF among internal streams. The result of this routine will be added to the input.

    ! Iterators.
    integer :: istr ! Over destination streams.
    integer :: istr_src ! Over source streams.
    integer :: ist ! Over destination Stokes parameters.

    ! TLC-parameters for surface reflection:
    ! T = 1
    ! L = 1
    ! C = C_srf

    ! Differentiating turns C into her derivative.
    ! We also have a prefactor and two intensity fields.

    ! Array in s  : intens
    ! Array in d  : prefactor_zip, intens_adj
    ! Matrix(s,d) : prefactor_scat, diff_c_srf

    ! Note that surface reflection is from down to up, and because up is
    ! down in the adjoint, both intens and intens_adj are down.
    do istr = 1,nstrhalf ! Destination stream.
      do istr_src = 1,nstrhalf ! Source stream.
        do ist = 1,nst ! Destination Stokes parameter.
          ! Array operation over source Stokes parameters.
          ! The V-flipper is indexed with ist, so it does no harm
          ! if there is no V-polarization.
          res(:,ist,istr_src,istr) = res(:,ist,istr_src,istr) + tlc%prefactor_zip(istr) * flip_v(ist) * intens_adj(ist,istr,idn,nsplit_total) * intens(:,istr_src,idn,nsplit_total) * tlc%prefactor_scat(istr_src,istr) * tlc%diff_c_bck_srf_bdrf(istr_src,istr)
        enddo
      enddo
    enddo

  end subroutine perturbation_bdrf_zip ! }}}
