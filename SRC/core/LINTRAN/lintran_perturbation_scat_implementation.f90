!> \file lintran_perturbation_scat_implementation.f90
!! Implementation file for derivatives with respect to scattering.
!! This file has a part of the contents from lintran_transfer_implementation, but here
!! derivatives are taken with respect to scattering C-parameters. These derivatives are
!! progressed to derivatives with respect to phase functions or single-scattering albedo.
!! Inner products are directly calculated.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Calculates the derivative of direct single-scatterig with respect to scattering parameters.
  !!
  !! This routine calculates the derivative of direct_1_elem from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers.
  subroutine perturbation_scat_direct_1_elem(nlay,nlay_deriv,ilay_deriv,tlc,res) ! {{{

    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay), intent(inout) :: res !< Resulting derivative with respect to the single-scattering C-parameter for the relevant Stokes parameters. The result of this routine will be added to the input.

    ! We differentiate by the phase matrix with
    ! rotation, so it is Z (not F). Their chain rules to the independent
    ! F-elements will be applied afterwards.

    ! For single-scattering, the viewing Stokes parameter is always equal to the
    ! destination Stokes parameter of the relevant matrix element. And the source
    ! Stokes parameter of the matrix element is always the Stokes parameter of the sun.
    ! So, this reduces the number of Stokes dimensions of the result to zero.

    ! TLC-parameters for backward scattering:
    ! T = T_tot_0v at interface ilay_pert-1
    ! L = L_0v
    ! C = Cb_0v_elem

    ! Differentiating turns C into her derivative. But for now, we
    ! differentiate with respect to C, so her derivative is one.
    res(ilay_deriv) = res(ilay_deriv) + tlc%prefactor_dir * tlc%t_tot_0v(ilay_deriv-1) * tlc%l_0v(ilay_deriv)

  end subroutine perturbation_scat_direct_1_elem ! }}}

  !> Calculates the derivative of direct double-scatterig with respect to scattering parameters.
  !!
  !! This routine calculates the derivative of direct_2 from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers.
  subroutine perturbation_scat_direct_2(nstrhalf,nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res_fwd_0,res_bck_0,res_fwd_v,res_bck_v) ! {{{

    use lintran_constants_module, only: nst_fake, ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_0 !< Resulting derivative with respect to C-parameters for forward scattering from the solar direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_0 !< Resulting derivative with respect to C-parameters for backward scattering from the solar direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_v !< Resulting derivative with respect to C-parameters for forward scattering to the viewing direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_v !< Resulting derivative with respect to C-parameters for backward scattering to the viewing direction. The result of this routine will be added to the input.

    ! Aggregate of prefactors and T- and L-parameters (the easy ones).
    real :: scalar_aggregated ! All dimensions are iterated.

    ! Iterators.
    integer :: ideriv ! Over differentiations.
    integer :: istr ! Over streams.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! TLC-parameters for (back - forward):
    ! T = T_tot_0v at interface ilay-1
    ! L = effmu_mv * (L_p0 - L_0v)
    ! C = Cb_0 * Cf_v

    ! TLC-parameters for (forward - back):
    ! T = T_tot_0v at interface ilay-1
    ! L = effmu_m0 * (L_pv - L_0v)
    ! C = Cf_0 * Cb_v

    ! Differentiating turns both C-parameters into their derivatives with the product rule.
    ! And the derivative of one individual C is one, because we differentiate with respect to C.

    do ideriv = 1,nlay_deriv
      ilay = ilay_deriv(ideriv) ! Layer that is differentiated.
      do istr = 1,nstrhalf

        ! (back - forward).
        scalar_aggregated = tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_0v(ilay-1) * tlc%effmu_mv(istr,ilay)*(tlc%l_p0(istr,ilay) - tlc%l_0v(ilay))
        ! The result is the aggregated scalar times one C-parameter if we differentiate with
        ! respect to the other C-parameter.

        ! C = Cb_0 * Cf_v.
        res_bck_0(:,istr,ilay) = res_bck_0(:,istr,ilay) + scalar_aggregated * tlc%c_fwd_v_nrm(:,istr,ilay)
        res_fwd_v(:,istr,ilay) = res_fwd_v(:,istr,ilay) + scalar_aggregated * tlc%c_bck_0_nrm(:,istr,ilay)

        ! (forward - back).
        ! Same joke, but switch sun and instrument in L-parameters and switch backward
        ! and forward in C.
        scalar_aggregated = tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_0v(ilay-1) * tlc%effmu_m0(istr,ilay)*(tlc%l_pv(istr,ilay) - tlc%l_0v(ilay))
        res_fwd_0(:,istr,ilay) = res_fwd_0(:,istr,ilay) + scalar_aggregated * tlc%c_bck_v_nrm(:,istr,ilay)
        res_bck_v(:,istr,ilay) = res_bck_v(:,istr,ilay) + scalar_aggregated * tlc%c_fwd_0_nrm(:,istr,ilay)

        if (exist_therm) then

          ! TLC-parameters for (thermal - back):
          ! T = T_v
          ! L = effmu_p * (L_v - L_pv)
          ! C = C_T * Cb_v

          ! TLC-parameters for (thermal - forward):
          ! T = T_v
          ! L = effmu_mv * (L_p - L_v)
          ! C = C_T * Cf_v

          ! Differentiation turns the atmospheric C into her derivative.
          scalar_aggregated = tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_v(ilay-1) * tlc%c_therm(ilay)
          res_bck_v(ist_therm,istr,ilay) = res_bck_v(ist_therm,istr,ilay) + scalar_aggregated * tlc%effmu_p(istr) * (tlc%l_v(ilay) - tlc%l_pv(istr,ilay))
          res_fwd_v(ist_therm,istr,ilay) = res_fwd_v(ist_therm,istr,ilay) + scalar_aggregated * tlc%effmu_mv(istr,ilay) * (tlc%l_p(istr,ilay) - tlc%l_v(ilay))

        endif
      enddo

    enddo

  end subroutine perturbation_scat_direct_2 ! }}}

  !> Calculates the inner products with the derivative of a first order source or response without layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1 from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_scat_source_response_1(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_v, fullpol_nst
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: flag_source !< True for source. False for response.
    logical, intent(in) :: flag_adjoint !< True for adjoint. False for normal run.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:,:), pointer :: l_px
    real, dimension(:,:), pointer :: l_mx

    ! Results that are being constructed from zero.
    ! This has to do with Delta-flipping for the adjoint.
    real, dimension(nst,nstrhalf,nlay) :: temp_res_fwd
    real, dimension(nst,nstrhalf,nlay) :: temp_res_bck

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Iterators.
    integer :: ideriv ! Over differentiations.
    integer :: istr ! Over streams.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_fwd = idn
      idir_bck = iup
    else
      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_fwd = iup
      idir_bck = idn
    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      l_px => tlc%l_pv
      l_mx => tlc%l_mv
    else
      ! Special direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      l_px => tlc%l_p0
      l_mx => tlc%l_m0
    endif

    ! TLC-parameters for backward scattering:
    ! T = T_tot_x at interface ilay-1
    ! L = L_px
    ! C = Cb_x

    ! TLC-parameters for forward scattering:
    ! T = T_tot_x at interface ilay-1 * T
    ! L = L_mx
    ! C = Cf_x

    ! Differentiating turns C into her derivative.

    do ideriv = 1,nlay_deriv
      ilay = ilay_deriv(ideriv) ! Layer that is differentiated.
      do istr = 1,nstrhalf
        ! The result has dimensions:
        ! 1: Stokes parameter (array operation).
        ! 2: Stream (istr)
        ! 3: Layer (ilay)
        temp_res_bck(:,istr,ilay) = prefactor(istr) * t_tot_x(ilay-1) * l_px(istr,ilay) * fld(:,istr,idir_bck,ilay-1)
        temp_res_fwd(:,istr,ilay) = prefactor(istr) * t_tot_x(ilay-1) * tlc%t(istr,ilay) * l_mx(istr,ilay) * fld(:,istr,idir_fwd,ilay)
      enddo
    enddo

    ! Interpret the result. We have to Delta-flip the result for the adjoint. And
    ! for the choice of sun and instrument, we have pointers.
    if (flag_adjoint .and. fullpol_nst(nst)) then

      ! Do Delta-flipping.
      temp_res_fwd(ist_v:ist_v,:,:) = -temp_res_fwd(ist_v:ist_v,:,:)
      temp_res_bck(ist_v:ist_v,:,:) = -temp_res_bck(ist_v:ist_v,:,:)

    endif

    ! Add results after the eventual necessary Delta-flipping.
    res_fwd_x = res_fwd_x + temp_res_fwd
    res_bck_x = res_bck_x + temp_res_bck

  end subroutine perturbation_scat_source_response_1 ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1_sp from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_scat_source_response_1_sp(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_v, fullpol_nst
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:,:), pointer :: l_sp_px
    real, dimension(:,:), pointer :: l_sp_mx

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Iterators.
    integer :: isplit ! Over split layers.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over streams.

    ! Derived index.
    integer :: ilay ! Atmospheric layers.

    ! Results that are being constructed from zero.
    ! This has to do with Delta-flipping for the adjoint.
    real, dimension(nst,nstrhalf,nlay) :: temp_res_fwd
    real, dimension(nst,nstrhalf,nlay) :: temp_res_bck

    ! Results of one sublayer.
    ! 0 is backward to low-index layer.
    ! 1 is forward to high-index layer.
    real, dimension(nstrhalf,0:1) :: sublayer_aggregated

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_fwd = idn
      idir_bck = iup
    else
      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_fwd = iup
      idir_bck = idn
    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_sp_px => tlc%l_sp_pv
      l_sp_mx => tlc%l_sp_mv
    else
      ! Special direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_sp_px => tlc%l_sp_p0
      l_sp_mx => tlc%l_sp_m0
    endif

    ! TLC-parameters for backward scattering:
    ! T = T_tot_x at interface above current sublayer
    ! L = L_px
    ! C = Cb_x

    ! Differentiating turns C into her derivative.

    ! TLC-parameters for forward scattering:
    ! T = T_tot_x at interface above current sublayer * T
    ! L = L_mx
    ! C = Cf_x

    ! Differentiating turns C into her derivative.

    ! Initialize temporary result to zero, because it is acquired by addition of sublayers
    temp_res_fwd = 0.0
    temp_res_bck = 0.0

    ! Start at topmost layer.
    isplit_start = 0
    ilay = 1
    do ideriv = 1,nlay_deriv

      ! Progress to the perturbed layer.
      isplit_start = isplit_start + sum(tlc%nsplit(ilay:ilay_deriv(ideriv)-1))
      ! Set layer index to perturbed layer.
      ilay = ilay_deriv(ideriv)
      ! Acquire last sublayer in layer. Note that ilay is now updated.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Acquire result of the first sublayer.
      sublayer_aggregated(:,0) = prefactor * t_tot_x(ilay-1) * l_sp_px(:,ilay)
      sublayer_aggregated(:,1) = prefactor * t_tot_x(ilay-1) * tlc%t_sp(:,ilay) * l_sp_mx(:,ilay)

      ! Apply sublayer result to all sublayers of this external layer.
      do isplit = isplit_start,isplit_end-1

        do istr = 1,nstrhalf
          temp_res_fwd(:,istr,ilay) = temp_res_fwd(:,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sublayer_aggregated(istr,1) * fld(:,istr,idir_fwd,isplit+1)
          temp_res_bck(:,istr,ilay) = temp_res_bck(:,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sublayer_aggregated(istr,0) * fld(:,istr,idir_bck,isplit)
        enddo

      enddo

    enddo

    ! Interpret the result. We have to Delta-flip the result for the adjoint. And
    ! for the choice of sun and instrument, we have pointers.
    if (flag_adjoint .and. fullpol_nst(nst)) then

      ! Do Delta-flipping.
      temp_res_fwd(ist_v:ist_v,:,:) = -temp_res_fwd(ist_v:ist_v,:,:)
      temp_res_bck(ist_v:ist_v,:,:) = -temp_res_bck(ist_v:ist_v,:,:)

    endif

    ! Add results after the eventual necessary Delta-flipping.
    res_fwd_x = res_fwd_x + temp_res_fwd
    res_bck_x = res_bck_x + temp_res_bck

  end subroutine perturbation_scat_source_response_1_sp ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting and coupling to interpolated field.
  !!
  !! This routine calculates the derivative of source_response_1_int from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_scat_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,fld,tlc,res_fwd_x,res_bck_x) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_v, fullpol_nst
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    ! Iterators.
    integer :: isplit ! Over sub-layers in one layer.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over streams.

    ! Derived index.
    integer :: ilay ! Atmospheric layers.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:), pointer :: l_int_sm_x
    real, dimension(:), pointer :: l_int_ac_x

    ! Results that are being constructed from zero.
    ! This has to do with Delta-flipping for the adjoint.
    real, dimension(nst,nstrhalf,nlay) :: temp_res_fwd
    real, dimension(nst,nstrhalf,nlay) :: temp_res_bck

    ! split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Recycled TLC-parameters for different split layers.
    real, dimension(nstrhalf,0:1) :: sublayer_aggregated

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_fwd = idn
      idir_bck = iup
    else
      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_fwd = iup
      idir_bck = idn
    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_int_sm_x => tlc%l_int_sm_v
      l_int_ac_x => tlc%l_int_ac_v
    else
      ! Special direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_int_sm_x => tlc%l_int_sm_0
      l_int_ac_x => tlc%l_int_ac_0
    endif

    ! TLC-parameters for backward scattering:
    ! T = T_tot_x at interface above current sublayer
    ! L = L_int_x
    ! C = Cb_x

    ! Differentiating turns C into her derivative.

    ! TLC-parameters for forward scattering:
    ! T = T_tot_x at interface above current sublayer
    ! L = L_int_x
    ! C = Cf_x

    ! Differentiating turns C into her derivative.

    ! Initialize temporary result to zero, because it is acquired by addition of sublayers.
    temp_res_fwd = 0.0
    temp_res_bck = 0.0

    ! Start at topmost layer.
    isplit_start = 0
    ilay = 1
    do ideriv = 1,nlay_deriv

      ! Progress to the perturbed layer.
      isplit_start = isplit_start + sum(tlc%nsplit(ilay:ilay_deriv(ideriv)-1))
      ! Set layer index to perturbed layer.
      ilay = ilay_deriv(ideriv)
      ! Acquire last sublayer in layer. Note that ilay is now updated.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Construct TLC-parameters for the sublayer, for both contributions of the interpolation.
      sublayer_aggregated(:,0) = prefactor * t_tot_x(ilay-1) * l_int_sm_x(ilay)
      sublayer_aggregated(:,1) = prefactor * t_tot_x(ilay-1) * l_int_ac_x(ilay)

      ! Apply sublayer result to all sublayers of this external layer.
      do isplit = isplit_start,isplit_end-1

        do istr = 1,nstrhalf
          ! Not the most impressive matmul, it sums two whole elements.
          temp_res_fwd(:,istr,ilay) = temp_res_fwd(:,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * matmul(fld(:,istr,idir_fwd,isplit:isplit+1) , sublayer_aggregated(istr,:))
          temp_res_bck(:,istr,ilay) = temp_res_bck(:,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * matmul(fld(:,istr,idir_bck,isplit:isplit+1) , sublayer_aggregated(istr,:))
        enddo

      enddo

    enddo

    ! Interpret the result. We have to Delta-flip the result for the adjoint. And
    ! for the choice of sun and instrument, we have pointers.
    if (flag_adjoint .and. fullpol_nst(nst)) then

      ! Do Delta-flipping.
      temp_res_fwd(ist_v:ist_v,:,:) = -temp_res_fwd(ist_v:ist_v,:,:)
      temp_res_bck(ist_v:ist_v,:,:) = -temp_res_bck(ist_v:ist_v,:,:)

    endif

    ! Add results after the eventual necessary Delta-flipping.
    res_fwd_x = res_fwd_x + temp_res_fwd
    res_bck_x = res_bck_x + temp_res_bck

  end subroutine perturbation_scat_source_response_1_int ! }}}

  !> Calculates the inner products with the derivative of a second order source or response with layer spltting and coupling to interpolated field.
  !!
  !! The corresponding non-derivative source_response_2_int is not needed and therefore does
  !! not exist in lintran_transfer_module. It is a source or response function where two
  !! scatter events in one layer are involved. One scatter event couples to an interpolated
  !! field and the other scatter events involves the solar or viewing direction. Directly, the
  !! inner product with an existing field is taken, so that derivatives for all layers can be
  !! handled at once.
  subroutine perturbation_scat_source_response_2_int(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_therm,fld,tlc,res_fwd,res_bck,res_fwd_x,res_bck_x) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_i, ist_u, ist_v, ist_therm, fullpol_nst
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_fwd_x !< Resulting derivative with respect to C-parameters for forward scattering involving the solar or viewing direction. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nstrhalf,nlay), intent(inout) :: res_bck_x !< Resulting derivative with respect to C-parameters for backward scattering involving the solar or viewing direction. The result of this routine will be added to the input.

    ! Iterators.
    integer :: isplit ! Over sub-layers in one layer.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over streams, usually destination streams.
    integer :: istr_src ! Over source streams.
    integer :: ist ! Over Stokes parameters.

    ! Derived index.
    integer :: ilay ! Atmospheric layers.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:,:), pointer :: effmu_px
    real, dimension(:,:), pointer :: effmu_mx
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:), pointer :: l_int_sm_x
    real, dimension(:), pointer :: l_int_ac_x
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x

    ! Temporary TLC-parameters, because C-parameters are asymmetric due to V-polarization.
    ! These are only needed if V-polarization is involved.
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_fwd_id
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_bck_id

    ! Results that are being constructed from zero.
    ! This has to do with Delta-flipping for the adjoint.
    real, dimension(nst,nstrhalf,nlay) :: temp_res_fwd
    real, dimension(nst,nstrhalf,nlay) :: temp_res_bck

    ! split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Results that may need to be transposed. For the response function, the
    ! phase function from i to d actually goes from d to i, so it must be
    ! transposed. In the routine, it is written in the perspective of the source.
    ! The result is added to the original result, so we should not transpose
    ! the entire result, but only the part that comes from this routine.
    real, dimension(nst,nst,nstrhalf,nstrhalf,nlay) :: res_fwd_untransposed
    real, dimension(nst,nst,nstrhalf,nstrhalf,nlay) :: res_bck_untransposed

    ! Aggregated sublayer in dimensions as if it is a C-parameter from i to d.
    ! First dimension: Stokes parameter i (intermediate).
    ! Second dimension: Stokes parameter d (destination).
    ! Third dimension: Stream i (intermediate).
    ! Fourth dimension: Stream d (destination).
    ! Fifth dimension: sm/ac.
    real, dimension(nst,nst,nstrhalf,nstrhalf,0:1) :: sublayer_aggregated_id
    ! We always need the destination direction, because we have to multiply with
    ! the field in the destination direction. We also always need the intermediate
    ! direction, because we resolve the phase function from i to d, or from x to i,
    ! both including i.

    ! Stokesless TLC-parameters (is all T and L and the prefactors).
    real, dimension(nstrhalf) :: stokesless_aggregated_i ! Direction i (intermediate).
    real, dimension(nstrhalf) :: stokesless_aggregated_d ! Direction d (destination).

    ! Set direction indices and pointers to TLC-parameters depending on
    ! source/response and normal/adjoint settings.
    if (flag_source) then

      prefactor => tlc%prefactor_src
      ! Forward-scattering direction is downward.
      idir_fwd = idn
      idir_bck = iup

      ! Asymmetric TLC-parameters.
      if (flag_adjoint) then
        c_fwd_x => tlc%c_fwd_v_adj
        c_bck_x => tlc%c_bck_v_adj
      else
        c_fwd_x => tlc%c_fwd_0_nrm
        c_bck_x => tlc%c_bck_0_nrm
      endif

    else

      prefactor => tlc%prefactor_resp
      ! Forward-scattering direction is upward.
      idir_fwd = iup
      idir_bck = idn

      ! Asymmetric TLC-parameters.
      if (flag_adjoint) then
        c_fwd_x => tlc%c_fwd_0_adj
        c_bck_x => tlc%c_bck_0_adj
      else
        c_fwd_x => tlc%c_fwd_v_nrm
        c_bck_x => tlc%c_bck_v_nrm
      endif

    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      effmu_mx => tlc%effmu_mv
      effmu_px => tlc%effmu_pv
      l_int_sm_x => tlc%l_int_sm_v
      l_int_ac_x => tlc%l_int_ac_v
    else
      ! Special direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      effmu_mx => tlc%effmu_m0
      effmu_px => tlc%effmu_p0
      l_int_sm_x => tlc%l_int_sm_0
      l_int_ac_x => tlc%l_int_ac_0
    endif

    ! Initialize these intermediate results.
    res_fwd_untransposed = 0.0
    res_bck_untransposed = 0.0

    ! Initialize temporary result to zero, because it is acquired by addition of sublayers.
    temp_res_fwd = 0.0
    temp_res_bck = 0.0

    ! Start at topmost layer.
    isplit_start = 0
    ilay = 1
    do ideriv = 1,nlay_deriv

      ! Progress to the perturbed layer.
      isplit_start = isplit_start + sum(tlc%nsplit(ilay:ilay_deriv(ideriv)-1))
      ! Set layer index to perturbed layer.
      ilay = ilay_deriv(ideriv)
      ! Acquire last sublayer in layer. Note that ilay is now updated.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Choose scattering from i to d (source) or from d to i (response).
      ! This transposition is necessary because of the asymmetry related to V-
      ! polarized radiation.

      ! Copy the data of the new layer.
      c_fwd_id = tlc%c_fwd(:,:,:,:,ilay)
      c_bck_id = tlc%c_bck(:,:,:,:,ilay)

      if (.not. flag_source .and. fullpol_nst(nst)) then

        ! Delta-flip with V.
        c_fwd_id(ist_v:ist_v,ist_i:ist_u,:,:) = -c_fwd_id(ist_v:ist_v,ist_i:ist_u,:,:)
        c_fwd_id(ist_i:ist_u,ist_v:ist_v,:,:) = -c_fwd_id(ist_i:ist_u,ist_v:ist_v,:,:)
        c_bck_id(ist_v:ist_v,ist_i:ist_u,:,:) = -c_bck_id(ist_v:ist_v,ist_i:ist_u,:,:)
        c_bck_id(ist_i:ist_u,ist_v:ist_v,:,:) = -c_bck_id(ist_i:ist_u,ist_v:ist_v,:,:)

      endif

      ! First, we will differentiate the second scatter event
      ! It can have any direction (fwd/bck).

      ! TLC-parameters for (forward - ?):
      ! The first scatter event is the analytic one and
      ! the second one couples to an interpolated field.

      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_mx(i) * (L_int_p(i) - L_int_x)
      ! C = Cf_x(i) * C?(i,d)

      ! Differentiating makes C? into nothing, but the differentiated C's Stokes
      ! parameters should be communicated to the other C.

      ! Scalar      : T_tot_x, L_x
      ! Array in i  : prefactor_intermediate, L, effmu_mx, Cf_x
      ! Array in d  : prefactor

      ! Aggregate everything and put the scalars with d, because there it is
      ! less crowded. First do the `sm' side.
      stokesless_aggregated_d = prefactor * t_tot_x(ilay-1)
      stokesless_aggregated_i = tlc%prefactor_intermediate * effmu_mx(:,ilay) * (tlc%l_int_sm_p(:,ilay) - l_int_sm_x(ilay))

      ! The missing TLC-parameters are Cf_x and C?.
      do istr = 1,nstrhalf ! Stream d.
        do istr_src = 1,nstrhalf ! Stream i.
          ! The Stokes parameters are not yet implemented.
          sublayer_aggregated_id(:,:,istr_src,istr,0) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
        enddo
      enddo

      ! Repeat joke for ac.
      ! Update the aggregate in which the `sm' parameter was and replace it with `ac'.
      stokesless_aggregated_i = tlc%prefactor_intermediate * effmu_mx(:,ilay) * (tlc%l_int_ac_p(:,ilay) - l_int_ac_x(ilay))

      ! The missing TLC-parameters are Cf_x and C?.
      do istr = 1,nstrhalf ! Stream d.
        do istr_src = 1,nstrhalf ! Stream i.
          ! The Stokes parameters are not yet implemented.
          sublayer_aggregated_id(:,:,istr_src,istr,1) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
        enddo
      enddo

      ! Apply result to field.
      do isplit = isplit_start,isplit_end-1

        do istr = 1,nstrhalf ! Stream i.
          do ist = 1,nst ! Stokes parameter i.
            ! Perturbation of i-d direction.
            ! The unperturbed scatter event is forward to x, so the net direction is
            ! equal to the perturbed scatter direction.
            res_fwd_untransposed(ist,:,istr,:,ilay) = res_fwd_untransposed(ist,:,istr,:,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_fwd,isplit:isplit+1),3) * c_fwd_x(ist,istr,ilay)
            res_bck_untransposed(ist,:,istr,:,ilay) = res_bck_untransposed(ist,:,istr,:,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_bck,isplit:isplit+1),3) * c_fwd_x(ist,istr,ilay)

            ! Perturbation of x-i direction.
            ! The unperturbed scatter event is from i to d. The net direction is equal
            ! to the this event. The perturbed scatter event is forward for the special
            ! direction.
            temp_res_fwd(ist,istr,ilay) = temp_res_fwd(ist,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(c_fwd_id(ist,:,istr,:) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_fwd,isplit:isplit+1),3))
            temp_res_fwd(ist,istr,ilay) = temp_res_fwd(ist,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(c_bck_id(ist,:,istr,:) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_bck,isplit:isplit+1),3))
          enddo
        enddo

      enddo

      ! TLC-parameters for (back - ?)
      ! The first scatter event is the analytic one and
      ! the second one couples to an interpolated field.

      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_px(i) * (L_int_x - T_px(i) * L_int_m(i))
      ! C = Cb_x(i) * C?(i,d)

      ! Differentiating does approximately the same as in the previous case.

      ! Scalar      : T_tot_x, L_int_x
      ! Array in i  : prefactor_intermediate, L_int_m, T_px, effmu_px, Cb_x
      ! Array in d  : prefactor

      ! We still have the same aggregate over d, so we keep that one.

      ! First `sm'.
      stokesless_aggregated_i = tlc%prefactor_intermediate * effmu_px(:,ilay) * (l_int_sm_x(ilay) - t_sp_x(ilay)*tlc%t_sp(:,ilay) * tlc%l_int_sm_m(:,ilay))

      ! The missing TLC-parameter are Cb_x and C?.
      do istr = 1,nstrhalf ! Stream d.
        do istr_src = 1,nstrhalf ! Stream i.
          ! The Stokes parameters are not yet implemented.
          sublayer_aggregated_id(:,:,istr_src,istr,0) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
        enddo
      enddo

      ! Repeat joke for ac.
      ! Update the aggregate in which the `sm' parameter was and replace it with `ac'.
      stokesless_aggregated_i = tlc%prefactor_intermediate * effmu_px(:,ilay) * (l_int_ac_x(ilay) - t_sp_x(ilay)*tlc%t_sp(:,ilay) * tlc%l_int_ac_m(:,ilay))

      ! The missing TLC-parameter are still Cb_x and C?.
      do istr = 1,nstrhalf ! Stream d.
        do istr_src = 1,nstrhalf ! Stream i.
          ! The Stokes parameters are not yet implemented.
          sublayer_aggregated_id(:,:,istr_src,istr,1) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
        enddo
      enddo

      ! Apply result to field.
      do isplit = isplit_start,isplit_end-1

        do istr = 1,nstrhalf ! Stream i.
          do ist = 1,nst ! Stokes parameter i.
            ! The unperturbed scatter event is backward, so the net direction is opposite to
            ! the perturbed scatter direction.
            res_fwd_untransposed(ist,:,istr,:,ilay) = res_fwd_untransposed(ist,:,istr,:,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_bck,isplit:isplit+1),3) * c_bck_x(ist,istr,ilay)
            res_bck_untransposed(ist,:,istr,:,ilay) = res_bck_untransposed(ist,:,istr,:,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_fwd,isplit:isplit+1),3) * c_bck_x(ist,istr,ilay)

            ! The unperturbed scatter event is from i to d and is in the opposite direction
            ! as the net (field) direction. The perturbed scatter event is backwards.
            temp_res_bck(ist,istr,ilay) = temp_res_bck(ist,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(c_bck_id(ist,:,istr,:) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_fwd,isplit:isplit+1),3))
            temp_res_bck(ist,istr,ilay) = temp_res_bck(ist,istr,ilay) + t_sp_x(ilay)**(isplit-isplit_start) * sum(c_fwd_id(ist,:,istr,:) * sum(sublayer_aggregated_id(ist,:,istr,:,:) * fld(:,:,idir_bck,isplit:isplit+1),3))
          enddo
        enddo

      enddo

      ! Thermal emission.
      if (exist_therm) then

        ! This only gives a contribution to the derivative of the C-parameter from i
        ! to d, because the first scatter event is the thermal emission itself.

        ! TLC-parameters for thermal emission:
        ! T = 1
        ! L = effmu_p * (L_int - L_int_p), all at intermediate direction
        ! C = C_T * C(i,d)
        ! These TLC-parameters are defined using the integral such that the
        ! `ac' side is defined as the layer where the just thermally emitted light
        ! goes to.

        ! However, thermal emission happens on both sides. This means that the
        ! destination light is going form `sm' to `ac' if the scatter event is
        ! forward and from `ac' to `sm' if the scatter event was backward.
        stokesless_aggregated_d = prefactor * tlc%c_therm(ilay)
        stokesless_aggregated_i = tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_sm_p(:,ilay))
        ! The missing TLC-parameter is only C?.
        do istr = 1,nstrhalf ! Stream d.
          do istr_src = 1,nstrhalf ! Stream i.
            ! The Stokes parameters are not yet implemented.
            ! For the intermediate direction, only Stokes parameter ist_therm is
            ! relevant. The rest will not be used.
            sublayer_aggregated_id(ist_therm,:,istr_src,istr,0) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
          enddo
        enddo
        ! Repeat joke for the `ac' side.
        stokesless_aggregated_i = tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_ac_p(:,ilay))
        ! The missing TLC-parameter is still only C?.
        do istr = 1,nstrhalf ! Stream d.
          do istr_src = 1,nstrhalf ! Stream i.
            ! The Stokes parameters are not yet implemented.
            ! For the intermediate direction, only Stokes parameter ist_therm is
            ! relevant. The rest will not be used.
            sublayer_aggregated_id(ist_therm,:,istr_src,istr,1) = stokesless_aggregated_d(istr) * stokesless_aggregated_i(istr_src)
          enddo
        enddo

        do istr = 1,nstrhalf ! Stream i.
          ! Stokes parameter i is ist_therm.
          ! If the sm side is the upper (low-indexed) one, the thermal emission was
          ! downward. Therefore, this TLC-aggregate can be used for the upper side
          ! in the combinations forward scattering and downward destination.
          ! Note that idir_fwd is downward radiation according to the intuition. For
          ! the adjoint, everything is reversed, but that will be okay. If we
          ! change one thing, sm / ac, or perturebed fwd / bck or destination
          ! bck / fwd, the indices switch between high and low.
          ! Combination 1 of 8: L = sm, perturbed = fwd, dest = idir_fwd.
          ! Then, we have the low indices.
          res_fwd_untransposed(ist_therm,:,istr,:,ilay) = res_fwd_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,0) * sum(fld(:,:,idir_fwd,isplit_start:isplit_end-1),3)
          ! Combination 2 of 8: L = sm, perturbed = fwd, dest = idir_bck.
          ! Then, we have the high indices.
          res_fwd_untransposed(ist_therm,:,istr,:,ilay) = res_fwd_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,0) * sum(fld(:,:,idir_bck,isplit_start+1:isplit_end),3)
          ! Combination 3 of 8: L = sm, perturbed = bck, dest = idir_fwd.
          ! Then, we have the high indices.
          res_bck_untransposed(ist_therm,:,istr,:,ilay) = res_bck_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,0) * sum(fld(:,:,idir_fwd,isplit_start+1:isplit_end),3)
          ! Combination 4 of 8: L = sm, perturbed = bck, dest = idir_bck.
          ! Then, we have the low indices.
          res_bck_untransposed(ist_therm,:,istr,:,ilay) = res_bck_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,0) * sum(fld(:,:,idir_bck,isplit_start:isplit_end-1),3)
          ! Combination 5 of 8: L = ac, perturbed = fwd, dest = idir_fwd.
          ! Then, we have the high indices.
          res_fwd_untransposed(ist_therm,:,istr,:,ilay) = res_fwd_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,1) * sum(fld(:,:,idir_fwd,isplit_start+1:isplit_end),3)

          ! Combination 6 of 8: L = ac, perturbed = fwd, dest = idir_bck.
          ! Then, we have the low indices.
          res_fwd_untransposed(ist_therm,:,istr,:,ilay) = res_fwd_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,1) * sum(fld(:,:,idir_bck,isplit_start:isplit_end-1),3)

          ! Combination 7 of 8: L = ac, perturbed = bck, dest = idir_fwd.
          ! Then, we have the low indices.
          res_bck_untransposed(ist_therm,:,istr,:,ilay) = res_bck_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,1) * sum(fld(:,:,idir_fwd,isplit_start:isplit_end-1),3)

          ! Combination 8 of 8: L = ac, perturbed = bck, dest = idir_bck.
          ! Then, we have the high indices.
          res_bck_untransposed(ist_therm,:,istr,:,ilay) = res_bck_untransposed(ist_therm,:,istr,:,ilay) + sublayer_aggregated_id(ist_therm,:,istr,:,1) * sum(fld(:,:,idir_bck,isplit_start+1:isplit_end),3)
        enddo

      endif

    enddo

    ! Transpose the result on i-d for response function if that is relevant.
    if (.not. flag_source .and. fullpol_nst(nst)) then
      res_fwd_untransposed(ist_i:ist_u,ist_v:ist_v,:,:,:) = -res_fwd_untransposed(ist_i:ist_u,ist_v:ist_v,:,:,:)
      res_fwd_untransposed(ist_v:ist_v,ist_i:ist_u,:,:,:) = -res_fwd_untransposed(ist_v:ist_v,ist_i:ist_u,:,:,:)
      res_bck_untransposed(ist_i:ist_u,ist_v:ist_v,:,:,:) = -res_bck_untransposed(ist_i:ist_u,ist_v:ist_v,:,:,:)
      res_bck_untransposed(ist_v:ist_v,ist_i:ist_u,:,:,:) = -res_bck_untransposed(ist_v:ist_v,ist_i:ist_u,:,:,:)
    endif
    ! Now the untransposed arrays are transposed.
    res_fwd = res_fwd + res_fwd_untransposed
    res_bck = res_bck + res_bck_untransposed

    ! Interpret the result. We have to Delta-flip the result for the adjoint. And
    ! for the choice of sun and instrument, we have pointers.
    if (flag_adjoint .and. fullpol_nst(nst)) then

      ! Do Delta-flipping
      temp_res_fwd(ist_v:ist_v,:,:) = -temp_res_fwd(ist_v:ist_v,:,:)
      temp_res_bck(ist_v:ist_v,:,:) = -temp_res_bck(ist_v:ist_v,:,:)

    endif

    ! Add results after the eventual necessary Delta-flipping.
    res_fwd_x = res_fwd_x + temp_res_fwd
    res_bck_x = res_bck_x + temp_res_bck

  end subroutine perturbation_scat_source_response_2_int ! }}}

  !> Calculates the inner products with the derivative of an interpolated internal scatter event.
  !!
  !! This routine calculates the derivative of scatter_1_int from lintran_transfer_module
  !! with respect to atmospheric scattering in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  !! As scatter_1_int also needs an input field, this routine needs both an input field for
  !! the scatter event and a field with which the inner product is taken.
  subroutine perturbation_scat_scatter_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res_fwd,res_bck) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Iterators.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over destination streams.
    integer :: istr_src ! Over source streams.
    integer :: ist ! Over destination Stokes parameters.

    ! Derived index.
    integer :: ilay ! Atmospheric layers.

    ! TLC-parameters for scattering:
    ! T = 1
    ! L = L_int_p(d)
    ! C = Cb or Cf

    ! Differentiating turns C into her derivative.

    ! We need to include the input intensity as well.

    ! Scalar      : None
    ! Array in s  : source
    ! Array in d  : L_int_p, destination
    ! Matrix(s,d) : prefactor_scat

    ! result = sum(array_d * matmul(array_s,matrix_sd))

    ! We will perform a loop over sublayers.

    ! Start at topmost layer.
    isplit_start = 0
    ilay = 1
    do ideriv = 1,nlay_deriv

      ! Progress to the perturbed layer.
      isplit_start = isplit_start + sum(tlc%nsplit(ilay:ilay_deriv(ideriv)-1))
      ! Set layer index to perturbed layer.
      ilay = ilay_deriv(ideriv)
      ! Acquire last sublayer in layer. Note that ilay is now updated.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      do istr = 1,nstrhalf
        do istr_src = 1,nstrhalf
          do ist = 1,nst

            ! Forward scattering, up-up, sm.
            ! Destination has low index, because it is up.
            ! Source has low index, because of `sm'.
            res_fwd(:,ist,istr_src,istr,ilay) = res_fwd(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * matmul(source(:,istr_src,iup,isplit_start:isplit_end-1) , destination(ist,istr,iup,isplit_start:isplit_end-1))
            ! Forward scattering, up-up, ac.
            ! Destination has low index, because it is up.
            ! Source has high index, because of `ac'.
            res_fwd(:,ist,istr_src,istr,ilay) = res_fwd(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * matmul(source(:,istr_src,iup,isplit_start+1:isplit_end) , destination(ist,istr,iup,isplit_start:isplit_end-1))
            ! Forward scattering, dn-dn, sm.
            ! Destination has high index, because it is dn.
            ! Source has high index, because of `sm'.
            res_fwd(:,ist,istr_src,istr,ilay) = res_fwd(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * matmul(source(:,istr_src,idn,isplit_start+1:isplit_end) , destination(ist,istr,idn,isplit_start+1:isplit_end))
            ! Forward scattering, dn-dn, ac.
            ! Destination has high index, because it is dn.
            ! Source has low index, because of `ac'.
            res_fwd(:,ist,istr_src,istr,ilay) = res_fwd(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * matmul(source(:,istr_src,idn,isplit_start:isplit_end-1) , destination(ist,istr,idn,isplit_start+1:isplit_end))
            ! Repeat joke for backward scattering, reversing the directions of the source
            ! Then, the directions and the indices can remain the same, since the
            ! link between the direction (dn/up) and a layer index is only there
            ! for the destination.
            res_bck(:,ist,istr_src,istr,ilay) = res_bck(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * matmul(source(:,istr_src,idn,isplit_start:isplit_end-1) , destination(ist,istr,iup,isplit_start:isplit_end-1))
            res_bck(:,ist,istr_src,istr,ilay) = res_bck(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * matmul(source(:,istr_src,idn,isplit_start+1:isplit_end) , destination(ist,istr,iup,isplit_start:isplit_end-1))
            res_bck(:,ist,istr_src,istr,ilay) = res_bck(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_sm_p(istr,ilay) * matmul(source(:,istr_src,iup,isplit_start+1:isplit_end) , destination(ist,istr,idn,isplit_start+1:isplit_end))
            res_bck(:,ist,istr_src,istr,ilay) = res_bck(:,ist,istr_src,istr,ilay) + tlc%prefactor_scat(istr_src,istr) * tlc%l_int_ac_p(istr,ilay) * matmul(source(:,istr_src,iup,isplit_start:isplit_end-1) , destination(ist,istr,idn,isplit_start+1:isplit_end))

          enddo
        enddo
      enddo

    enddo

  end subroutine perturbation_scat_scatter_1_int ! }}}

  !> Takes the inner product of a normal and adjoint intensity field with a pertrubation of the radiative transfer operator.
  !!
  !! Perturbs the scattering coefficient in radiative transfer operator and performs the inner
  !! product of the adjoint field with the perturbation of the radiative transfer operator acting
  !! on the normal field. Corrections for the pseudo-forward method for the adjoint field are
  !! applied so that the inner product is correct. Both fields are interpolated. Note that the
  !! radiative transfer operator is perturbed, not the matrix. And because scattering is
  !! perturbed, there is not such a thing as transmission. Therefore, the result of the pretrubed
  !! operator acting on the field does not do anything with the vertical grid, so the field
  !! remains interpolated. Chain rule coefficients convert derivatives with respect to the
  !! scattering coefficient into derivatives with respect to C-parameters and to the optical
  !! depths. This is the only derivative with respect to the optical depth that is in this module,
  !! just because it is the same routine as a derivative with respect to C. The rest of the
  !! derivatives with respect to the optical depth are in lintran_perturbation_tau_module.
  subroutine perturbation_scat_zip(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,nlay_deriv_base,ilay_deriv_base,intens,intens_adj,tlc,chain_c,res_fwd,res_bck,chain_base,res_base) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, flip_v
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay), intent(in) :: chain_c !< Chain rule factors for derivatives with respect to C.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_fwd !< Resulting derivative with respect to C-parameters for forward scattering among internal streams. The result of this routine will be added to the input.
    real, dimension(nst+0*nst_fake,nst+0*nst_fake,nstrhalf,nstrhalf,nlay), intent(inout) :: res_bck !< Resulting derivative with respect to C-parameters for backward scattering among internal streams. The result of this routine will be added to the input.
    real, dimension(nlay), intent(in) :: chain_base !< Chain rule factors for derivatives with respect to the optical depth, the basic property that is relevant.
    real, dimension(nlay_deriv_base), intent(inout) :: res_base !< Resulting derivative with respect to the optical depth, the basic property that is relevant. The result of this routine will be added to the input.

    ! Iterators.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over streams.
    integer :: ideriv ! Over differentiations.
    integer :: ist ! Over Stokes parameters.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! Administration of writing results for tau (base).
    logical :: deriv_base ! Flag for writing results for tau.
    integer :: ideriv_base ! Differentiation index for differentiation to tau.

    ! Temporary result.
    real, dimension(nst,nst) :: intensproduct

    ! Result without chain rule.
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: aggregated

    ! The derivative of the radiative transfer operator is the
    ! derivative of C times a new L that takes responsibility of the
    ! mutual interpolations of the normal and the adjoint fields.
    ! For the prefactors, there is a scatter event from the direction
    ! of the forward to the direction of the adjoint, plus a zip factor
    ! at the direction of the adjoint (prefactor_scat(s,d)*prefactor_zip(d))
    ! This is the same as prefactor_scat(d,s)*prefactor_zip(s), then the
    ! last scatter event had been done by the adjoint.

    ! We have to care for the fact that the adjoint field is still in
    ! the reverse direction.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Start at topmost layer.
    isplit_start = 0
    ilay = 1
    ideriv_base = 0
    do ideriv = 1,nlay_deriv_c

      ! Progress to the perturbed layer.
      isplit_start = isplit_start + sum(tlc%nsplit(ilay:ilay_deriv_c(ideriv)-1))
      ! Set layer index to perturbed layer.
      ilay = ilay_deriv_c(ideriv)
      ! Acquire last sublayer in layer. Note that ilay is now updated.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Check if the layer is included in the base derivatives.
      deriv_base = .false.
      if (ideriv_base .lt. nlay_deriv_base) then ! Prevents out-of-bounds exceptions.
        ! Check if the next base-derivative layer is indeed ilay.
        if (ilay .eq. ilay_deriv_base(ideriv_base+1)) then
          ideriv_base = ideriv_base + 1
          deriv_base = .true.
        endif
      endif

      do istr = 1,nstrhalf ! Destination stream.

        do istr_src = 1,nstrhalf ! Source stream.

          ! TLC-parameters for backward scattering:
          ! T = 1
          ! L = l_zip
          ! C = Cb

          ! Differentiating turns C into her derivative.

          ! Scalar      : l_zip
          ! Array in s  : intens
          ! Array in d  : intens_adj, prefactor_zip
          ! Matrix(s,d) : prefactor_scat

          ! For each sublayer, there will be a contribution with l_zip_sm for
          ! intens(isplit_low)*intesn_adj(isplit_low) + intens(isplit_high)*intens_adj(isplit_high)

          ! That is once the border interfaces and twice the internal interfaces.

          ! We will perform the matrix multiplication and thereby zip the
          ! dimension over the split layers within external layer ilay
          ! We will compute the intensity product per combination of Stokes
          ! parameters for each combination of directions. Then the common
          ! TLC-parameters are applied, including the derivative of C, that
          ! is there to collapse the Stokes parameters again.

          ! Note that backward scattering zips equal directions, because the
          ! adjoint is defined in a reverse direction. Forward scattering zips
          ! different directions. Layer split indices have no relationship with forward
          ! or backward scattering, only with sm/ac interpolation contribution, as zipping
          ! happens with interpolated fields on both sides.

          ! The product between the normal and adjoint intensities is split in some
          ! components, because we have to add up, two directions (up and down), and we
          ! have to handle the interior interfaces differently from the outer
          ! interfaces, because the interior count twice, because they
          ! contribute to two sublayers.

          ! Contribution for same-interface intensities.
          do ist = 1,nst ! Destination Stokes parameter.
            intensproduct(:,ist) = flip_v(ist) * (intens(:,istr_src,idn,isplit_start) * intens_adj(ist,istr,idn,isplit_start) + intens(:,istr_src,iup,isplit_start) * intens_adj(ist,istr,iup,isplit_start) + intens(:,istr_src,idn,isplit_end) * intens_adj(ist,istr,idn,isplit_end) + intens(:,istr_src,iup,isplit_end) * intens_adj(ist,istr,iup,isplit_end) + 2.0 * matmul(intens(:,istr_src,idn,isplit_start+1:isplit_end-1) , intens_adj(ist,istr,idn,isplit_start+1:isplit_end-1)) + 2.0 * matmul(intens(:,istr_src,iup,isplit_start+1:isplit_end-1) , intens_adj(ist,istr,iup,isplit_start+1:isplit_end-1)))
          enddo
          aggregated(:,:,istr_src,istr) = tlc%l_zip_sm(ilay) * tlc%prefactor_zip(istr) * tlc%prefactor_scat(istr_src,istr) * intensproduct

          ! Contribution for other-interface intensities.
          do ist = 1,nst ! Destination stream.
            intensproduct(:,ist) = flip_v(ist) * (matmul(intens(:,istr_src,idn,isplit_start:isplit_end-1) , intens_adj(ist,istr,idn,isplit_start+1:isplit_end)) + matmul(intens(:,istr_src,iup,isplit_start:isplit_end-1) , intens_adj(ist,istr,iup,isplit_start+1:isplit_end)) + matmul(intens(:,istr_src,idn,isplit_start+1:isplit_end) , intens_adj(ist,istr,idn,isplit_start:isplit_end-1)) + matmul(intens(:,istr_src,iup,isplit_start+1:isplit_end) , intens_adj(ist,istr,iup,isplit_start:isplit_end-1)))
          enddo
          aggregated(:,:,istr_src,istr) = aggregated(:,:,istr_src,istr) + tlc%l_zip_ac(ilay) * tlc%prefactor_zip(istr) * tlc%prefactor_scat(istr_src,istr) * intensproduct

        enddo

      enddo

      ! Apply chain rules.
      ! Because the zipping L-parameter integrates such that the integral becomes a
      ! derivative with respect to C*tau instead of just C, we have to apply chain
      ! rules for C and for tau.
      ! Furthermore, for the base (tau) integral, we have to get rid of the phase function
      ! (C-differentiation) and transform it first to a derivative to taus_s instead of
      ! C*tau, because that is the original perspective. Therefore, the tau-derivative
      ! needs the derivative of C with respect to the single-scattering albedo.
      res_bck(:,:,:,:,ilay) = res_bck(:,:,:,:,ilay) + chain_c(ilay) * aggregated
      if (deriv_base) then
        res_base(ideriv_base) = res_base(ideriv_base) + chain_base(ilay) * sum(aggregated * tlc%diff_c_bck_ssa(:,:,:,:,ilay))
      endif

      do istr = 1,nstrhalf ! Destination stream.

        do istr_src = 1,nstrhalf ! Source stream.

          ! TLC-parameters for forward scattering:
          ! The same, except for Cb that becomes Cf, so we will now contribute to
          ! the forward result. And the indices idn and iup now should be different
          ! for intens and intens_adj.

          ! Same layers.
          do ist = 1,nst ! Destination Stokes parameter.
            intensproduct(:,ist) = flip_v(ist) * (intens(:,istr_src,idn,isplit_start) * intens_adj(ist,istr,iup,isplit_start) + intens(:,istr_src,iup,isplit_start) * intens_adj(ist,istr,idn,isplit_start) + intens(:,istr_src,idn,isplit_end) * intens_adj(ist,istr,iup,isplit_end) + intens(:,istr_src,iup,isplit_end) * intens_adj(ist,istr,idn,isplit_end) + 2.0 * matmul(intens(:,istr_src,idn,isplit_start+1:isplit_end-1) , intens_adj(ist,istr,iup,isplit_start+1:isplit_end-1)) + 2.0 * matmul(intens(:,istr_src,iup,isplit_start+1:isplit_end-1) , intens_adj(ist,istr,idn,isplit_start+1:isplit_end-1)))
          enddo
          aggregated(:,:,istr_src,istr) = tlc%l_zip_sm(ilay) * tlc%prefactor_zip(istr) * tlc%prefactor_scat(istr_src,istr) * intensproduct

          ! Across layers.
          do ist = 1,nst ! Destination stream.
            intensproduct(:,ist) = flip_v(ist) * (matmul(intens(:,istr_src,idn,isplit_start:isplit_end-1) , intens_adj(ist,istr,iup,isplit_start+1:isplit_end)) + matmul(intens(:,istr_src,iup,isplit_start:isplit_end-1) , intens_adj(ist,istr,idn,isplit_start+1:isplit_end)) + matmul(intens(:,istr_src,idn,isplit_start+1:isplit_end) , intens_adj(ist,istr,iup,isplit_start:isplit_end-1)) + matmul(intens(:,istr_src,iup,isplit_start+1:isplit_end) , intens_adj(ist,istr,idn,isplit_start:isplit_end-1)))
          enddo
          aggregated(:,:,istr_src,istr) = aggregated(:,:,istr_src,istr) + tlc%l_zip_ac(ilay) * tlc%prefactor_zip(istr) * tlc%prefactor_scat(istr_src,istr) * intensproduct

        enddo

      enddo

      ! Apply chain rules, with the same remarks as for backward scattering.
      res_fwd(:,:,:,:,ilay) = res_fwd(:,:,:,:,ilay) + chain_c(ilay) * aggregated
      if (deriv_base) then
        res_base(ideriv_base) = res_base(ideriv_base) + chain_base(ilay) * sum(aggregated * tlc%diff_c_fwd_ssa(:,:,:,:,ilay))
      endif

    enddo

  end subroutine perturbation_scat_zip ! }}}
