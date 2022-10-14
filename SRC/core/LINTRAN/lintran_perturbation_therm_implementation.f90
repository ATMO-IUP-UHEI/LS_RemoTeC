!> \file lintran_perturbation_therm_implementation.f90
!! Implementation file for derivatives with respect to thermal emissivity.
!! This file has a part of the contents from lintran_transfer_implementation, but here
!! derivatives are taken with respect to thermal emissivity. Inner products are directly
!! calculated.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Calculates the derivative of direct single-scatterig with respect to thermal emissivity.
  !!
  !! This routine calculates the derivative of direct_1_elem from lintran_transfer_module
  !! with respect to thermal emissivity.
  subroutine perturbation_therm_direct_1_elem(nlay_deriv,ilay_deriv,tlc,res) ! {{{

    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the thermal emissivity. The result of this routine will be added to the input.

    ! Thermal emission is only relevant if the viewing Stokes parameter is equal to the
    ! thermally emitted Stokes parameter (I), because this is single-scattering. And the
    ! thermal emission event is counted as a scatter event. Therefore, thermally emitted
    ! light can only change Stokes parameter from the second scatter event.

    ! TLC-parameters for thermal emission:
    ! T = T_tot_v at interface ilay-1
    ! L = L_v
    ! C = C_T

    ! Differentiating will turn C into her derivative.
    ! But this derivative is handled later on, because chains to
    ! the Planck curve and ssa.
    res = res + tlc%prefactor_dir * tlc%t_tot_v(ilay_deriv-1) * tlc%l_v(ilay_deriv)

  end subroutine perturbation_therm_direct_1_elem ! }}}

  !> Calculates the derivative of direct double-scatterig with respect to thermal emissivity.
  !!
  !! This routine calculates the derivative of direct_2 from lintran_transfer_module
  !! with respect to the thermal emissivity in the selected layers.
  subroutine perturbation_therm_direct_2(nlay_deriv,ilay_deriv,tlc,res) ! {{{

    use lintran_constants_module, only: ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterator.
    integer :: ideriv ! Over differentiations.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! TLC-parameters for (thermal - back):
    ! T = T_v
    ! L = effmu_p * (L_v - L_pv)
    ! C = C_T * Cb_v

    ! Differentiation turns C_T into her derivative.

    ! TLC-parameters for (thermal - forward):
    ! T = T_v
    ! L = effmu_mv * (L_p - L_v)
    ! C = C_T * Cf_v

    ! Differentiation turns C_T into her derivative.

    ! The first term is (thermal - back), the second term is (thermal - forward).
    do ideriv = 1,nlay_deriv
      ilay = ilay_deriv(ideriv) ! Layer that is differentiated.
      res(ideriv) = res(ideriv) + tlc%prefactor_dir * sum(tlc%prefactor_intermediate * tlc%t_tot_v(ilay-1) * (tlc%effmu_p * (tlc%l_v(ilay) - tlc%l_pv(:,ilay)) * tlc%c_bck_v_nrm(ist_therm,:,ilay) + tlc%effmu_mv(:,ilay) * (tlc%l_p(:,ilay) - tlc%l_v(ilay)) * tlc%c_fwd_v_nrm(ist_therm,:,ilay)))
    enddo

  end subroutine perturbation_therm_direct_2 ! }}}

  !> Calculates the inner products with the derivative of a first order source or response without layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1 from lintran_transfer_module
  !! with respect to thermal emission in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_therm_source_response_1(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: flag_source !< True for source. False for response.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointer to TLC-parameter with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor

    ! Iterator.
    integer :: ideriv ! Over differentiations.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! Set direction indices.
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

    ! This is only relevant for adjoint response or normal source. Otherwise,
    ! this routine should not be called.

    ! TLC-parameters for thermal emission:
    ! T = 1
    ! L = L_p
    ! C = C_T

    ! Differentiation turns C into her derivative.
    do ideriv = 1,nlay_deriv
      ilay = ilay_deriv(ideriv) ! Layer that is differentiated.
      res(ideriv) = res(ideriv) + sum(prefactor * tlc%l_p(:,ilay) * (fld(ist_therm,:,idir_bck,ilay-1) + fld(ist_therm,:,idir_fwd,ilay)))
    enddo

  end subroutine perturbation_therm_source_response_1 ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1_sp from lintran_transfer_module
  !! with respect to thermal emission in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_therm_source_response_1_sp(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: flag_source !< True for source. False for response.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: ideriv ! Over perturbed layers.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointer to TLC-parameter with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Set direction indices.
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

      ! TLC-parameters for thermal emission:
      ! T = 1
      ! L = L_p
      ! C = C_T

      ! Differentiation turns C intro her derivative.

      res(ideriv) = res(ideriv) + sum(matmul(prefactor * tlc%l_sp_p(:,ilay) , fld(ist_therm,:,idir_bck,isplit_start:isplit_end-1)+fld(ist_therm,:,idir_fwd,isplit_start+1:isplit_end)))

    enddo

  end subroutine perturbation_therm_source_response_1_sp ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting and coupling to interpolated field.
  !!
  !! This routine calculates the derivative of source_response_1_int from lintran_transfer_module
  !! with respect to thermal emission in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_therm_source_response_1_int(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: flag_source !< True for source. False for response.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: ideriv ! Over perturbed layers.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointer to TLC-parameter with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Set direction indices.
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

      ! TLC-parameters for thermal emission:
      ! T = 1
      ! L = L_int
      ! C = C_T

      ! Differentiation turns C intro her derivative.

      ! The matmul is over the streams, the left sum is over iup and idn and the sum on
      ! the right is over split layers.
      res(ideriv) = res(ideriv) + tlc%l_int(ilay) * sum(matmul(prefactor , fld(ist_therm,:,:,isplit_start) + 2.0 * sum(fld(ist_therm,:,:,isplit_start+1:isplit_end-1),3) + fld(ist_therm,:,:,isplit_end)))

    enddo

  end subroutine perturbation_therm_source_response_1_int ! }}}

  !> Calculates the inner products with the derivative of a second order source or response with layer spltting and coupling to interpolated field.
  !!
  !! The corresponding non-derivative source_response_2_int is not needed and therefore does
  !! not exist in lintran_transfer_module. It is a source or response function where two
  !! scatter events in one layer are involved. One scatter event couples to an interpolated
  !! field and the other scatter event is the thermal emission itself. Directly, the
  !! inner product with an existing field is taken, so that derivatives for all layers can be
  !! handled at once.
  subroutine perturbation_therm_source_response_2_int(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,flag_source,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_i, ist_u, ist_v, ist_therm, fullpol_nst
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: flag_source !< True for source. False for response.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Pointer to TLC-parameter with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor

    ! Aggregate intermediate result.
    real, dimension(nst,nstrhalf), target :: sublayer_aggregated

    ! Pointers for folding Stokes and stream dimensions.
    real, dimension(:), pointer :: sublayer_aggregated_ptr
    real, dimension(:,:,:), pointer :: fld_ptr

    ! Temporary TLC-parameters, because C-parameters are asymmetric due to V-polarization.
    ! These are only needed if V-polarization is involved.
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_fwd_id
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_bck_id

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Iterators.
    integer :: ideriv ! Over differentiations.
    integer :: istr ! Over streams.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! The pointers are fixed. As we will be matmulling over all the split layers inside
    ! one layer (because we have no transmission term), we will set one pointer for the
    ! entire field, just folding the first two dimensions.
    fld_ptr(1:nst*nstrhalf,0:1,0:nsplit_total) => fld
    sublayer_aggregated_ptr(1:nst*nstrhalf) => sublayer_aggregated

    ! Set direction indices and prefactor, which depend on if it is the source
    ! or the response function.
    if (flag_source) then
      prefactor => tlc%prefactor_src
    else
      prefactor => tlc%prefactor_resp
    endif

    ! This is only relevant for adjoint response or normal source. Otherwise,
    ! this routine should not be called.

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

      ! TLC-parameters for thermal emission:
      ! T = 1
      ! L = effmu_p * (L_int - L_int_p), all at intermediate direction
      ! C = C_T * C(i,d)

      ! The emission is towards the `ac' side, with these TLC-parameters.

      ! Differentiation turns C_T into her derivative.

      ! Case 1: Backward scattering `sm' side.
      do istr = 1,nstrhalf
        ! The intermediate Stokes parameters is only ist_therm.
        ! The intermediate streams are matmulled.
        ! THe destination Stokes parameter is the array operation.
        ! The destination stream is looped.
        sublayer_aggregated(:,istr) = prefactor(istr) * matmul(c_bck_id(ist_therm,:,:,istr) , tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_sm_p(:,ilay)))
      enddo
      ! After backward scattering, the target radiation is directed towards the `sm' side.
      ! `sm' contribution, upward radiation with low indices, downward radiation with
      ! high indices.
      res(ideriv) = res(ideriv) + sum(matmul(sublayer_aggregated_ptr , fld_ptr(:,iup,isplit_start:isplit_end-1) + fld_ptr(:,idn,isplit_start+1:isplit_end))) ! Matmul over combined Stokes and stream dimensions. Sum over split layers.

      ! For the prose about dimensions, iterations and matmulling, read case 1. The new
      ! comments are only about things that are different from case 1.

      ! Case 2: Backward scattering `ac' side.
      do istr = 1,nstrhalf
        sublayer_aggregated(:,istr) = prefactor(istr) * matmul(c_bck_id(ist_therm,:,:,istr) , tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_ac_p(:,ilay)))
      enddo
      ! After backward scattering, the target radiation is directed towards the `sm' side.
      ! `ac' contribution, upward radiation with high indices, downward radiation with
      ! low indices.
      res(ideriv) = res(ideriv) + sum(matmul(sublayer_aggregated_ptr , fld_ptr(:,iup,isplit_start+1:isplit_end) + fld_ptr(:,idn,isplit_start:isplit_end-1)))

      ! Case 3: Forward scattering `sm' side.
      do istr = 1,nstrhalf
        sublayer_aggregated(:,istr) = prefactor(istr) * matmul(c_fwd_id(ist_therm,:,:,istr) , tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_sm_p(:,ilay)))
      enddo
      ! After forward scattering, the target radiation is directed towards the `ac' side.
      ! `sm' contribution, upward radiation with high indices, downward radiation with
      ! low indices.
      res(ideriv) = res(ideriv) + sum(matmul(sublayer_aggregated_ptr , fld_ptr(:,iup,isplit_start+1:isplit_end) + fld_ptr(:,idn,isplit_start:isplit_end-1)))

      ! Case 4: Forward scattering `ac' side.
      do istr = 1,nstrhalf
        sublayer_aggregated(:,istr) = prefactor(istr) * matmul(c_fwd_id(ist_therm,:,:,istr) , tlc%prefactor_intermediate * tlc%effmu_p * (tlc%l_int(ilay) - tlc%l_int_ac_p(:,ilay)))
      enddo
      ! After forward scattering, the target radiation is directed towards the `ac' side.
      ! `ac' contribution, upward radiation with low indices, downward radiation with
      ! high indices.
      res(ideriv) = res(ideriv) + sum(matmul(sublayer_aggregated_ptr , fld_ptr(:,iup,isplit_start:isplit_end-1) + fld_ptr(:,idn,isplit_start+1:isplit_end)))

    enddo

  end subroutine perturbation_therm_source_response_2_int ! }}}
