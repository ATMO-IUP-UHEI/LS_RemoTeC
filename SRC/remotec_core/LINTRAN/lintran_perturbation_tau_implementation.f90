!> \file lintran_perturbation_tau_implementation.f90
!! Implementation file for derivatives with respect to the optical depth.
!! This file has a part of the contents from lintran_transfer_implementation, but here
!! derivatives are taken with respect to optical depth. Inner products are directly
!! calculated.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Calculates the derivative of direct single-scatterig with respect to the optical depth.
  !!
  !! This routine calculates the derivative of direct_1_elem from lintran_transfer_module
  !! with respect to the optical depths in the selected layers.
  subroutine perturbation_tau_direct_1_elem(nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in), target :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth, for the relevant viewing Stokes parameter. The result of this routine will be added to the input.

    ! Aggregate of common TLC-parameters.
    real :: aggregated ! Aggregate of TLC-parameters for atmospheric scattering.
    real :: aggregated_therm ! Aggregate of TLC-parameters for thermal emission.

    ! Iterator.
    integer :: ilay ! Over scattering layers.

    ! Adminitration for diagonal terms, of which we do not know where they are.
    logical :: diag ! Diagonal term is included in result.
    integer :: ideriv_diag ! Index of diagonal term in result.

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end
    integer :: ideriv_end_therm ! Possibly different sensitivity for just the instrument.

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr

    ! Pointers to sensitive parts for thermal emission.
    real, dimension(:), pointer :: res_ptr_therm
    integer, dimension(:), pointer :: ilay_deriv_ptr_therm

    ! TLC-parameters for backward scattering:
    ! T = T_tot_0v at interface ilay-1
    ! L = L_0v at scattering layer
    ! C = Cb_0v at scattering layer

    ! Differentiation.
    ! There is a relative derivative of T and an absolute derivative of L.

    ! TLC-parameters for thermal emission:
    ! T = T_v at interface ilay-1
    ! L = L_v
    ! C = C_T

    ! Differentiation.
    ! adds relative derivative of T and absolute derivative of L.

    ideriv_diag = 0 ! Only relevant for thermal emission.
    do ilay = 1,nlay

      ! Set sensitivity pointers.
      res_ptr => res(1:tlc%ideriv_sensitive_0v(ilay))
      ilay_deriv_ptr => ilay_deriv(1:tlc%ideriv_sensitive_0v(ilay))

      ! TLC-parameters for single atmospheric scattering:
      ! T = T_tot_0v at interface ilay-1
      ! L = L_0v
      ! C = C_bck_0v_elem

      ! Differentation involves L and T.
      ! Derivative of L is diagonal under plane-parallel conditions.

      ! Derivative of T.
      aggregated = tlc%prefactor_dir * tlc%t_tot_0v(ilay-1) * tlc%c_bck_0v_elem(ilay)
      res_ptr = res_ptr + tlc%reldiff_t_tot_0v_tau(ilay_deriv_ptr,ilay-1) * tlc%l_0v(ilay) * aggregated

      ! Derivatives of T for thermal emission.

      if (exist_therm) then

        ! Reset pointers, because now only the instrument is relevant. The sun does not matter.
        res_ptr_therm => res(1:tlc%ideriv_sensitive_v(ilay))
        ilay_deriv_ptr_therm => ilay_deriv(1:tlc%ideriv_sensitive_v(ilay))
        ! TLC-parameters for thermal emission:
        ! T = T_v at interface ilay-1
        ! L = L_v
        ! C = C_T
        aggregated_therm = tlc%prefactor_dir * tlc%t_tot_v(ilay-1) * tlc%c_therm(ilay)
        res_ptr_therm = res_ptr_therm + tlc%reldiff_t_tot_v_tau(ilay_deriv_ptr_therm,ilay-1) * tlc%l_v(ilay) * aggregated_therm

      endif

      ! Derivative of L are diagonal for plane-parallel geometry.

      ! Check for diagonal terms.
      if (tlc%plane_parallel) then
        diag = .false.
        if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
          if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
            diag = .true.
            ideriv_diag = ideriv_diag + 1
          endif
        endif
        if (.not. diag) cycle
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
        if (exist_therm) ideriv_end_therm = ideriv_diag
      else
        ideriv_start = 1
        ideriv_end = tlc%ideriv_sensitive_0v(ilay)
        if (exist_therm) ideriv_end_therm = tlc%ideriv_sensitive_v(ilay)
      endif

      ! Reset pointers.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)
      if (exist_therm) then
        res_ptr_therm => res(ideriv_start:ideriv_end_therm)
        ilay_deriv_ptr_therm => ilay_deriv(ideriv_start:ideriv_end_therm)
      endif

      ! Apply derivatives of L.
      res_ptr = res_ptr + tlc%diff_l_0v_tau(ilay_deriv_ptr,ilay) * aggregated
      if (exist_therm) res_ptr_therm = res_ptr_therm + tlc%diff_l_v_tau(ilay_deriv_ptr_therm,ilay) * aggregated_therm

    enddo

    ! TLC-parameters for surface reflection:
    ! T = T_tot_0v at bottom
    ! L = 1
    ! C = Cb_srf_0v

    ! TLC-parameters for fluorescent emission:
    ! T = T_tot_v at bottom
    ! L = 1
    ! C = Ce_v

    ! Differentiation.
    ! adds relative derivative of T for both cases.
    res = res + tlc%prefactor_dir * (tlc%t_tot_0v(nlay) * tlc%reldiff_t_tot_0v_tau(ilay_deriv,nlay) * tlc%c_bck_srf_0v_elem + tlc%t_tot_v(nlay) * tlc%reldiff_t_tot_v_tau(ilay_deriv,nlay) * tlc%c_emi_v_elem)

  end subroutine perturbation_tau_direct_1_elem ! }}}

  !> Calculates the derivative of direct double-scatterig with respect to the optical depth.
  !!
  !! This routine calculates the derivative of direct_2 from lintran_transfer_module
  !! with respect to the optical depths in the selected layers.
  subroutine perturbation_tau_direct_2(nstrhalf,nlay,nlay_deriv,ilay_deriv,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: ist_therm
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

    ! Aggregate of common TLC-parameters.
    real, dimension(nstrhalf) :: arr_aggregated_i_fwd_bck ! Dimension: intermediate streams.
    real, dimension(nstrhalf) :: arr_aggregated_i_bck_fwd ! Dimension: intermediate streams.
    real, dimension(nstrhalf) :: arr_aggregated_i_therm_bck ! Dimension: intermediate streams.
    real, dimension(nstrhalf) :: arr_aggregated_i_therm_fwd ! Dimension: intermediate streams.

    ! Iterators
    integer :: ilay ! Over scattering layers
    integer :: istr ! Over streams

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr

    ! Pointers to sensitive parts for thermal emission.
    real, dimension(:), pointer :: res_ptr_therm
    integer, dimension(:), pointer :: ilay_deriv_ptr_therm

    ! Adminitration for diagonal terms, of which we do not know where they are.
    integer :: ideriv_diag ! Index of diagonal term in result.
    logical :: diag ! Diagonal term is included in result.

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end
    integer :: ideriv_end_therm ! Possibly different sensitivity for just the instrument.

    ideriv_diag = 0 ! Only relevant for thermal emission.
    do ilay = 1,nlay

      ! TLC-parameters for (back - forward):
      ! T = T_tot_0v at interface ilay-1
      ! L = effmu_mv * (L_p0 - L_0v)
      ! C = Cb_0 * Cf_v

      ! Differentiation.
      ! There will be a relative derivative of T and an absolute derivative of L.

      ! Common TLC-parameters: two prefactors, T and two Cs.
      ! Derivative parts:
      ! - Differentiate effmu_mv (with L_p0 - L_0v)
      ! - Differentiate L_p0 (with effmu_mv)
      ! - Differentiate L_0v (with -effmu_mv)
      ! And the relative derivative of T.

      ! When looping over scattering layers, the following dimensions remain:
      ! Scalar : prefactor_dir, T, L_0v
      ! Array in streams : prefactor_intermediate, Cb_0, Cf_v, effmu_mv, L_p0
      ! Array in perturbed layers: DT_v, DL_0v
      ! Matrix (stream,perturbed layer) Deffmu_mv, DL_p0

      ! The C-parameters have to be summed over the intermediate Stokes parameters,
      ! and the intermediate stream index remains as dimension, because that one
      ! interacts with T and L-parameters.

      ! For (forward - back), all TLC-parameters are the same, but sun and
      ! instrument are reversed.

      ! Set sensitivity pointers.
      res_ptr => res(1:tlc%ideriv_sensitive_0v(ilay))
      ilay_deriv_ptr => ilay_deriv(1:tlc%ideriv_sensitive_0v(ilay))

      ! Loop needed because of the way the Cs are multiplied.
      do istr = 1,nstrhalf
        ! (back - forward).
        arr_aggregated_i_bck_fwd(istr) = tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_0v(ilay-1) * dot_product(tlc%c_bck_0_nrm(:,istr,ilay) , tlc%c_fwd_v_nrm(:,istr,ilay))
        ! (forward - back).
        arr_aggregated_i_fwd_bck(istr) = tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_0v(ilay-1) * dot_product(tlc%c_bck_v_nrm(:,istr,ilay) , tlc%c_fwd_0_nrm(:,istr,ilay))
      enddo

      ! Derivative of T, also nondiagonal for plane-parallel geometry. The first term is
      ! back-forward. The second term is forward-back.
      res_ptr = res_ptr + tlc%reldiff_t_tot_0v_tau(ilay_deriv_ptr,ilay-1) * (sum(arr_aggregated_i_bck_fwd * tlc%effmu_mv(:,ilay) * (tlc%l_p0(:,ilay) - tlc%l_0v(ilay))) + sum(arr_aggregated_i_fwd_bck * (tlc%effmu_m0(:,ilay) * (tlc%l_pv(:,ilay) - tlc%l_0v(ilay)))))

      if (exist_therm) then

        ! Reset sensitivity pointers to just the instrument.
        res_ptr_therm => res(1:tlc%ideriv_sensitive_v(ilay))
        ilay_deriv_ptr_therm => ilay_deriv(1:tlc%ideriv_sensitive_v(ilay))

        ! TLC-parameters for (thermal - back):
        ! T = T_v
        ! L = effmu_p * (L_v - L_pv)
        ! C = C_T * Cb_v

        ! Differentiation.
        ! There will be a relative derivative of T and an absolute derivative of L.

        ! Common TLC-parameters: two prefactors, T, two Cs and effmu_p.
        ! Derivative parts:
        ! - Differentiate L_v
        ! - Differentiate L_pv (with minus sign)
        ! And the relative derivative of T.

        ! No loop required, becuase the thermal C-parameter does not depend on the
        ! intermediate stream. And we can apply the thermal Stokes parameter to the
        ! other C-parameter.
        arr_aggregated_i_therm_bck = tlc%prefactor_dir * tlc%t_tot_v(ilay-1) * tlc%c_therm(ilay) * tlc%prefactor_intermediate * tlc%effmu_p * tlc%c_bck_v_nrm(ist_therm,:,ilay)

        ! TLC-parameters for (thermal - forward):
        ! T = T_v
        ! L = effmu_mv * (L_p - L_v)
        ! C = C_T * Cf_v

        ! Differentiation.
        ! There will be a relative derivative of T and an absolute derivative of L.

        ! Common TLC-parameters: two prefactors, T, and two Cs.
        ! Derivative parts:
        ! - Differentiate effmu_mv (with L_p - L_v)
        ! - Differentiate L_p (with effmu_mv), only diagonal
        ! - Differentiate L_v (with -efmmu_mv)
        ! And the relative derivative of T.

        ! Apply same non-loop joke as for the backward counterpart, using the constant thermal Stokes parameter.
        arr_aggregated_i_therm_fwd = tlc%prefactor_dir * tlc%t_tot_v(ilay-1) * tlc%c_therm(ilay) * tlc%prefactor_intermediate * tlc%c_fwd_v_nrm(ist_therm,:,ilay)

        ! Now just do the derivatives of T.
        res_ptr_therm = res_ptr_therm + tlc%reldiff_t_tot_v_tau(ilay_deriv_ptr_therm,ilay-1) * (sum(arr_aggregated_i_therm_bck * (tlc%l_v(ilay) - tlc%l_pv(:,ilay))) + sum(arr_aggregated_i_therm_fwd * tlc%effmu_mv(:,ilay) * (tlc%l_p(:,ilay) - tlc%l_v(ilay))))

      endif

      ! Check for diagonal terms.
      if (tlc%plane_parallel .or. exist_therm) then
        diag = .false.
        if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
          if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
            diag = .true.
            ideriv_diag = ideriv_diag + 1
          endif
        endif
      endif

      ! Only diagonal terms are relevant for plane_parallel.
      if (tlc%plane_parallel .and. .not. diag) cycle

      ! Set rules for loop over ideriv.
      if (tlc%plane_parallel) then
        ! For plane-parallel, there must be diagonal index, otherwise, we would not be here.
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
        ideriv_end_therm = ideriv_diag
      else
        ideriv_start = 1
        ideriv_end = tlc%ideriv_sensitive_0v(ilay)
        ideriv_end_therm = tlc%ideriv_sensitive_v(ilay)
      endif

      ! Redefine sensitivity pointers.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)
      if (exist_therm) then
        res_ptr_therm => res(ideriv_start:ideriv_end_therm)
        ilay_deriv_ptr_therm => ilay_deriv(ideriv_start:ideriv_end_therm)
      endif

      ! Apply derivative of L. The first terms are back-forward. The end is forward-back.
      res_ptr = res_ptr + matmul(arr_aggregated_i_bck_fwd * (tlc%l_p0(:,ilay) - tlc%l_0v(ilay)) , tlc%diff_effmu_mv_tau(:,ilay_deriv_ptr,ilay)) + matmul(arr_aggregated_i_bck_fwd * tlc%effmu_mv(:,ilay) , tlc%diff_l_p0_tau(:,ilay_deriv_ptr,ilay)) - tlc%diff_l_0v_tau(ilay_deriv_ptr,ilay) * sum(arr_aggregated_i_bck_fwd * tlc%effmu_mv(:,ilay)) + matmul(arr_aggregated_i_fwd_bck * (tlc%l_pv(:,ilay) - tlc%l_0v(ilay)) , tlc%diff_effmu_m0_tau(:,ilay_deriv_ptr,ilay)) + matmul(arr_aggregated_i_fwd_bck * tlc%effmu_m0(:,ilay) , tlc%diff_l_pv_tau(:,ilay_deriv_ptr,ilay)) - tlc%diff_l_0v_tau(ilay_deriv_ptr,ilay) * sum(arr_aggregated_i_fwd_bck * tlc%effmu_m0(:,ilay))

      if (exist_therm) then

        ! Apply derivatives of L.
        res_ptr_therm = res_ptr_therm + tlc%diff_l_v_tau(ilay_deriv_ptr_therm,ilay) * sum(arr_aggregated_i_therm_bck) - matmul(arr_aggregated_i_therm_bck , tlc%diff_l_pv_tau(:,ilay_deriv_ptr_therm,ilay)) + matmul(arr_aggregated_i_therm_fwd * (tlc%l_p(:,ilay) - tlc%l_v(ilay)) , tlc%diff_effmu_mv_tau(:,ilay_deriv_ptr_therm,ilay)) - tlc%diff_l_v_tau(ilay_deriv_ptr_therm,ilay) * sum(arr_aggregated_i_therm_fwd * tlc%effmu_mv(:,ilay))

        ! The diagonal term.
        if (diag) res(ideriv_diag) = res(ideriv_diag) + sum(arr_aggregated_i_therm_fwd * tlc%diff_l_p_tau(:,ilay) * tlc%effmu_mv(:,ilay))

      endif

    enddo

  end subroutine perturbation_tau_direct_2 ! }}}

  !> Calculates the derivative of transmission with respect to the optical depth wihtout layer splitting.
  !!
  !! To simplify the differentiation, it is assumed the the radiation is at the interface just
  !! before the perturbed layer. In a transmission event through zero or more layers like
  !! in transmit from lintran_transfer_module, the radiation crosses multiple layers at once.
  !! To apply the perturbation, the first zero or more layers should be crossed unperturbedly
  !! with transmit. This routine handles the perturbed layer. The rest of the layers will
  !! be transmitted by transmit again. But that is done by applying a reverse transmisstion
  !! to the field with which the inner product is taken.
  subroutine perturbation_tau_mutualtransmit(nstrhalf,nlay,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: source !< Intensity that is at the interface of the perturbed layer.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: destination !< Field with which the inner product is taken. This field includes an inverse transmission, so that the perturbation can be done on only one layer.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over streams.

    ! Non-iterated index.
    integer :: ilay ! Atmospheric layers.

    ! This inner product is performed between two fields that both had been
    ! transmitted towards each other. This inner product perturbs the last layer
    ! before the fields meet each other.

    do ideriv = 1,nlay_deriv

      ilay = ilay_deriv(ideriv) ! Perturbed layer.

      ! We let the source and the destination meet at a perturbed layer.
      ! The destination is some adjoint field that was transmitted in the reverse.
      ! Therefore, it is no pseudo-forward adjoint, so the term destination is
      ! better than adjoint. It is the destination after a perturbed transmission
      ! through one layer.

      ! TLC-parameters for transmitting though one layer:
      ! T = T
      ! L = 1
      ! C = 1
      ! Differentiation adds the relative derivative of T.
      ! There is no prefactor.
      do istr = 1,nstrhalf
        res(ideriv) = res(ideriv) + tlc%t(istr,ilay)*tlc%reldiff_t_tau(istr) * (sum(source(:,istr,idn,ilay-1)*destination(:,istr,idn,ilay)) + sum(source(:,istr,iup,ilay)*destination(:,istr,iup,ilay-1)))
      enddo

    enddo

  end subroutine perturbation_tau_mutualtransmit ! }}}

  !> Calculates the derivative of transmission with respect to the optical depth wiht layer splitting.
  !!
  !! To simplify the differentiation, it is assumed the the radiation is at the interface just
  !! before a perturbed sublayer. In a transmission event through zero or more sublayers like
  !! in transmit_sp from lintran_transfer_module, the radiation crosses multiple sublayers at once.
  !! To apply the perturbation, the first zero or more sublayers should be crossed unperturbedly
  !! with transmit_sp. This routine handles the perturbed sublayer. The rest of the sublayers will
  !! be transmitted by transmit_sp again. But that is done by applying a reverse transmisstion
  !! to the field with which the inner product is taken.
  subroutine perturbation_tau_mutualtransmit_sp(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Intensity that is at the interface of a perturbed sublayer.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken. This field includes an inverse transmission, so that the perturbation can be done on only one sublayer.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ideriv ! Over differentiated layers.
    integer :: istr ! Over streams.

    ! Non-iterated index.
    integer :: ilay ! Atmospheric layers.

    ! Borders of a layers in split environment.
    integer :: isplit_start
    integer :: isplit_end

    ! This inner product is performed between two fields that both had been
    ! transmitted towards each other. This inner product perturbs the last sublayer
    ! before the fields meet each other.

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

      ! We let the source and the destination meet at a perturbed sublayer.
      ! The destination is some adjoint field that was transmitted in the reverse.
      ! Therefore, it is no pseudo-forward adjoint, so the term destination is
      ! better than adjoint. It is the destination after a perturbed transmission
      ! through one sublayer.

      ! TLC-parameters for transmitting though one sublayer:
      ! T = T
      ! L = 1
      ! C = 1
      ! Differentiation adds the relative derivative of T.
      ! There is no prefactor.
      ! Meanwhile, the sublayers are summed up.
      do istr = 1,nstrhalf
        res(ideriv) = res(ideriv) + tlc%t_sp(istr,ilay)*tlc%reldiff_t_sp_tau(istr,ilay) * (sum(source(:,istr,idn,isplit_start:isplit_end-1)*destination(:,istr,idn,isplit_start+1:isplit_end)) + sum(source(:,istr,iup,isplit_start+1:isplit_end)*destination(:,istr,iup,isplit_start:isplit_end-1)))
      enddo

    enddo

  end subroutine perturbation_tau_mutualtransmit_sp ! }}}

  !> Calculates the inner products with the derivative of a first order source or response without layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1 from lintran_transfer_module
  !! with respect to the opcial depth in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_tau_source_response_1(nstrhalf,nlay,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Pre-aggregated common TLC-parameters, in destination stream dimension.
    real, dimension(nstrhalf) :: arr_aggregated_bck
    real, dimension(nstrhalf) :: arr_aggregated_fwd

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: istr ! Over streams.

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr

    ! Adminitration for diagonal terms, of which we do not know where they are.
    integer :: ideriv_diag ! Index of diagonal term in result.
    logical :: diag ! Diagonal term is included in result.

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    integer, dimension(:), pointer :: ideriv_sensitive_x
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:,:), pointer :: l_px
    real, dimension(:,:), pointer :: l_mx
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: reldiff_t_tot_x_tau
    real, dimension(:,:,:), pointer :: diff_l_px_tau
    real, dimension(:,:,:), pointer :: diff_l_mx_tau

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
        c_bck_srf_x => tlc%c_bck_srf_v_adj
      else
        c_fwd_x => tlc%c_fwd_0_nrm
        c_bck_x => tlc%c_bck_0_nrm
        c_bck_srf_x => tlc%c_bck_srf_0_nrm
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
        c_bck_srf_x => tlc%c_bck_srf_0_adj
      else
        c_fwd_x => tlc%c_fwd_v_nrm
        c_bck_x => tlc%c_bck_v_nrm
        c_bck_srf_x => tlc%c_bck_srf_v_nrm
      endif

    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      ideriv_sensitive_x => tlc%ideriv_sensitive_v
      t_tot_x => tlc%t_tot_v
      l_px => tlc%l_pv
      l_mx => tlc%l_mv
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_v_tau
      diff_l_px_tau => tlc%diff_l_pv_tau
      diff_l_mx_tau => tlc%diff_l_mv_tau
    else
      ! Direction is the sun (0).
      ideriv_sensitive_x => tlc%ideriv_sensitive_0
      t_tot_x => tlc%t_tot_0
      l_px => tlc%l_p0
      l_mx => tlc%l_m0
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_0_tau
      diff_l_px_tau => tlc%diff_l_p0_tau
      diff_l_mx_tau => tlc%diff_l_m0_tau
    endif

    ideriv_diag = 0
    do ilay = 1,nlay

      ! Set sensitivity pointers.
      res_ptr => res(1:ideriv_sensitive_x(ilay))
      ilay_deriv_ptr => ilay_deriv(1:ideriv_sensitive_x(ilay))

      ! Start with derivatives of transmission terms using the original fields. This is
      ! even nondiagonal for plane-parallel geometry.

      ! TLC-parameters for backward scattering:
      ! T = T_tot_x at interface ilay-1
      ! L = L_px
      ! C = Cb_x

      ! TLC-parameters for forward scattering:
      ! T = T_tot_x at interface ilay-1 * T
      ! L = L_mx
      ! C = Cf_x

      ! Aggregate common TLC-parameters.
      do istr = 1,nstrhalf
        arr_aggregated_bck(istr) = prefactor(istr) * t_tot_x(ilay-1) * dot_product(c_bck_x(:,istr,ilay) , fld(:,istr,idir_bck,ilay-1))
        arr_aggregated_fwd(istr) = prefactor(istr) * t_tot_x(ilay-1) * tlc%t(istr,ilay) * dot_product(c_fwd_x(:,istr,ilay) , fld(:,istr,idir_fwd,ilay))
      enddo

      ! Differentiation now adds relative derivative of T.
      ! This is even off-diagonal for plane-parallel geometry.
      res_ptr = res_ptr + reldiff_t_tot_x_tau(ilay_deriv_ptr,ilay-1) * (dot_product(arr_aggregated_bck , l_px(:,ilay)) + dot_product(arr_aggregated_fwd , l_mx(:,ilay)))

      ! Check for diagonal terms.
      diag = .false.
      if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
        if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
          diag = .true.
          ideriv_diag = ideriv_diag + 1
        endif
      endif

      ! Only diagonal terms are relevant for plane_parallel.
      if (tlc%plane_parallel .and. .not. diag) cycle

      ! Set rules for loop over ideriv.
      if (tlc%plane_parallel) then
        ! For plane-parallel, there must be diagonal index, otherwise, we would not be here.
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
      else
        ideriv_start = 1
        ideriv_end = ideriv_sensitive_x(ilay)
      endif

      ! We repeat the TLC-parameters, mentioning the different
      ! differentiations.

      ! TLC-parameters for backward scattering:
      ! T = T_tot_x at interface ilay-1
      ! L = L_px
      ! C = Cb_x

      ! Differentiations.
      ! - only L

      ! TLC-parameters for forward scattering:
      ! T = T_tot_x at interface ilay-1 * T
      ! L = L_mx
      ! C = Cf_x

      ! Differenations.
      ! - L_mx (with nothing)
      ! - T (with L_mx)

      ! Redefine sensitivity pointers.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)

      ! Apply derivative of L.
      res_ptr = res_ptr + matmul(arr_aggregated_bck , diff_l_px_tau(:,ilay_deriv_ptr,ilay)) + matmul(arr_aggregated_fwd , diff_l_mx_tau(:,ilay_deriv_ptr,ilay))

      ! The derivative of T for in the forward part is only diagonal.
      ! Also, the thermal emission is only diagonal, because with thermal
      ! emission, there is no sun, so also no spherical Earth.

      ! So, only diagonal terms remain, even for pseudo-spherical geometry.
      if (.not. diag) cycle

      ! The diagonal derivative of the single-layer T in the internal stream for
      ! the forward scattering term.
      res(ideriv_diag) = res(ideriv_diag) + sum(arr_aggregated_fwd * tlc%reldiff_t_tau * l_mx(:,ilay))

      ! The thermal term is also diagonal.
      if (exist_therm) then

        ! TLC-parameters for thermal emission:
        ! T = 1
        ! L = L_p
        ! C = C_T

        ! Differentiation is on L_p.
        ! There is no pre-aggregated array.
        ! We will take the forward and backward emitted light at the same time.
        res(ideriv_diag) = res(ideriv_diag) + tlc%c_therm(ilay) * sum(prefactor * tlc%diff_l_p_tau(:,ilay) * (fld(ist_therm,:,idir_bck,ilay-1) + fld(ist_therm,:,idir_fwd,ilay)))

      endif

    enddo

    ! TLC-parameters for surface reflection:
    ! T = t_tot_x at bottom
    ! L = 1
    ! C = Cb_srf

    ! Add the relative derivative of the transmission term.
    if (exist_bdrf) then
      do istr = 1,nstrhalf
        res = res + prefactor(istr) * reldiff_t_tot_x_tau(ilay_deriv,nlay) * t_tot_x(nlay) * sum((c_bck_srf_x(:,istr) * fld(:,istr,idir_bck,nlay)))
      enddo
    endif

  end subroutine perturbation_tau_source_response_1 ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting.
  !!
  !! This routine calculates the derivative of source_response_1_sp from lintran_transfer_module
  !! with respect to the opcial depth in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_tau_source_response_1_sp(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: isplit ! Over sublayers.
    integer :: ideriv ! Over perturbed layers.

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr
    real, dimension(:,:), pointer :: sublayer_aggregated_for_l_bck_ptr ! Also reshapes.
    real, dimension(:,:), pointer :: sublayer_aggregated_for_l_fwd_ptr ! Also reshapes.
    real, dimension(:,:,:), pointer :: fld_ptr ! Only reshapes.

    ! Adminitration for diagonal terms, of which we do not know where they are.
    logical :: diag ! Diagonal term is included in result.
    integer :: ideriv_diag ! Index of diagonal term in result.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    integer, dimension(:), pointer :: ideriv_sensitive_x
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:,:), pointer :: l_sp_px
    real, dimension(:,:), pointer :: l_sp_mx
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: reldiff_t_tot_x_tau
    real, dimension(:,:), pointer :: reldiff_t_sp_x_tau
    real, dimension(:,:,:), pointer :: diff_l_sp_px_tau
    real, dimension(:,:,:), pointer :: diff_l_sp_mx_tau

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end

    ! Iterator.
    integer :: istr ! Over streams.

    ! Aggregate of TLC-parameters relevant for one sublayer. One for differentiating
    ! with respect to T, and one for differentiating with respect to L.
    ! The one for L is more complicated, because it combines different kinds of contributions
    ! in different terms of the derivatives in L. The one of T is more dynamic, because it gets
    ! an extra dynamic factor of a changing relative derivative of T for different sublayers.
    real, dimension(nst,nstrhalf,nlay_deriv), target :: sublayer_aggregated_for_l_bck
    real, dimension(nst,nstrhalf,nlay_deriv), target :: sublayer_aggregated_for_l_fwd
    ! We will use the original field for sublayer_aggregated_for_t, because different
    ! scattering layers are not contaminated.
    real, dimension(nst,nstrhalf) :: sublayer_aggregated_for_t_bck
    real, dimension(nst,nstrhalf) :: sublayer_aggregated_for_t_fwd

    ! An aggregate of TLC-parameters that can be re-used in a loop over perturbed layers.
    real, dimension(nst,nstrhalf) :: arr_aggregated_bck
    real, dimension(nst,nstrhalf) :: arr_aggregated_fwd

    ! Set fixed reshape pointer.
    fld_ptr(1:nst*nstrhalf,0:1,0:nsplit_total) => fld

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
        c_bck_srf_x => tlc%c_bck_srf_v_adj
      else
        c_fwd_x => tlc%c_fwd_0_nrm
        c_bck_x => tlc%c_bck_0_nrm
        c_bck_srf_x => tlc%c_bck_srf_0_nrm
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
        c_bck_srf_x => tlc%c_bck_srf_0_adj
      else
        c_fwd_x => tlc%c_fwd_v_nrm
        c_bck_x => tlc%c_bck_v_nrm
        c_bck_srf_x => tlc%c_bck_srf_v_nrm
      endif

    endif

    ! Direction indices and TLC-parameters. The direction that either
    ! corresponds to the sun or to the instrument is abbreviated with 'x'.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      ideriv_sensitive_x => tlc%ideriv_sensitive_v
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_sp_px => tlc%l_sp_pv
      l_sp_mx => tlc%l_sp_mv
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_v_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_v_tau
      diff_l_sp_px_tau => tlc%diff_l_sp_pv_tau
      diff_l_sp_mx_tau => tlc%diff_l_sp_mv_tau
    else
      ! Special direction is the sun (0).
      ideriv_sensitive_x => tlc%ideriv_sensitive_0
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_sp_px => tlc%l_sp_p0
      l_sp_mx => tlc%l_sp_m0
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_0_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_0_tau
      diff_l_sp_px_tau => tlc%diff_l_sp_p0_tau
      diff_l_sp_mx_tau => tlc%diff_l_sp_m0_tau
    endif

    isplit_start = 0
    ideriv_diag = 0
    do ilay = 1,nlay

      ! Acquire end sublayer index.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Set sensitivity pointers.
      res_ptr => res(1:ideriv_sensitive_x(ilay))
      ilay_deriv_ptr => ilay_deriv(1:ideriv_sensitive_x(ilay))

      ! Contribution of derivative of transmission term using the original fields,
      ! which is even nondiagonal for plane-parallel gemoetry.

      ! TLC-parameters for backward scattering:
      ! T = T_tot_x at interface above current sublayer
      ! L = L_px
      ! C = Cb_x

      ! TLC-parameters for forward scattering:
      ! T = T_tot_x at interface above current sublayer * T
      ! L = L_mx
      ! C = Cf_x

      do istr = 1,nstrhalf
        arr_aggregated_bck(:,istr) = prefactor(istr) * t_tot_x(ilay-1) * c_bck_x(:,istr,ilay)
        arr_aggregated_fwd(:,istr) = prefactor(istr) * t_tot_x(ilay-1) * tlc%t_sp(istr,ilay) * c_fwd_x(:,istr,ilay)
      enddo

      ! Add the non-derivative L-parameter to become th total sublayer aggregate for T,
      ! where everything except the derivatives of T are in.
      do istr = 1,nstrhalf
        sublayer_aggregated_for_t_bck(:,istr) = arr_aggregated_bck(:,istr) * l_sp_px(istr,ilay)
        sublayer_aggregated_for_t_fwd(:,istr) = arr_aggregated_fwd(:,istr) * l_sp_mx(istr,ilay)
      enddo

      ! Differentiation over T. This is even non-diagonal for plane-parallel geometry.
      do isplit = isplit_start,isplit_end-1

        res_ptr = res_ptr + (reldiff_t_tot_x_tau(ilay_deriv_ptr,ilay-1) + (isplit-isplit_start) * reldiff_t_sp_x_tau(ilay_deriv_ptr,ilay)) * t_sp_x(ilay)**(isplit-isplit_start) * (sum(sublayer_aggregated_for_t_bck * fld(:,:,idir_bck,isplit)) + sum(sublayer_aggregated_for_t_fwd * fld(:,:,idir_fwd,isplit+1)))

      enddo

      ! Check for diagonal terms.
      diag = .false.
      if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
        if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
          diag = .true.
          ideriv_diag = ideriv_diag + 1
        endif
      endif

      ! Further calcuation only necessary for pseudo-spherical calculation, or just the
      ! diagonal terms for plane-parallel.
      if (tlc%plane_parallel .and. .not. diag) then
        ! Progress to next layer and go out of the loop.
        isplit_start = isplit_end
        cycle
      endif

      ! Set rules for loop over ideriv.
      if (tlc%plane_parallel) then
        ! For plane-parallel, there must be diagonal index, otherwise, we would not be here.
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
      else
        ideriv_start = 1
        ideriv_end = ideriv_sensitive_x(ilay)
      endif

      ! Repeat TLC-parameters, considering the differentiation of L.

      ! TLC-parameters for backward scattering:
      ! T = T_tot_x at interface above current sublayer
      ! L = L_px
      ! C = Cb_x

      ! TLC-parameters for forward scattering:
      ! T = T_tot_x at interface above current sublayer * T
      ! L = L_mx
      ! C = Cf_x

      ! Redefine sensitivity pointers, including the TLC-aggregate for differentiation of L.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)
      sublayer_aggregated_for_l_bck_ptr(1:nst*nstrhalf,ideriv_start:ideriv_end) => sublayer_aggregated_for_l_bck(:,:,ideriv_start:ideriv_end) ! Also reshapes by merging Stokes parameters and streams.
      sublayer_aggregated_for_l_fwd_ptr(1:nst*nstrhalf,ideriv_start:ideriv_end) => sublayer_aggregated_for_l_fwd(:,:,ideriv_start:ideriv_end) ! Also reshapes by merging Stokes parameters and streams.

      ! Loop over perturbed layers for the sublayer results for L.
      do ideriv = ideriv_start,ideriv_end

        ! Calculate sublayer terms for L.
        ! They are all nondiagonal.
        do istr = 1,nstrhalf
          sublayer_aggregated_for_l_bck(:,istr,ideriv) = arr_aggregated_bck(:,istr) * diff_l_sp_px_tau(istr,ilay_deriv(ideriv),ilay)
          sublayer_aggregated_for_l_fwd(:,istr,ideriv) = arr_aggregated_fwd(:,istr) * diff_l_sp_mx_tau(istr,ilay_deriv(ideriv),ilay)
        enddo

      enddo

      ! Add diagonal term from the relative derivative of T.
      if (diag) then
        do istr = 1,nstrhalf
          sublayer_aggregated_for_l_fwd(:,istr,ideriv_diag) = sublayer_aggregated_for_l_fwd(:,istr,ideriv_diag) + tlc%reldiff_t_sp_tau(istr,ilay) * arr_aggregated_fwd(:,istr) * l_sp_mx(istr,ilay)
        enddo
      endif

      ! Apply results to all sublayers.
      do isplit = isplit_start,isplit_end-1

        res_ptr = res_ptr + t_sp_x(ilay)**(isplit-isplit_start) * (matmul(fld_ptr(:,idir_bck,isplit) , sublayer_aggregated_for_l_bck_ptr) + matmul(fld_ptr(:,idir_fwd,isplit+1) , sublayer_aggregated_for_l_fwd_ptr))

      enddo

      ! Diagonal contribution from thermal emission. It is separated because it does
      ! not have the same transmission term.
      if (exist_therm .and. diag) then

        ! TLC-parameters for thermal emission:
        ! T = 1
        ! L = L_p
        ! C = C_T

        ! Differentiation turns L into her derivative.

        ! This is independent of direction.
        ! We only have to fill in the thermal Stokes parameter.
        ! We do not use pre-aggregated results, because this happens only once.
        res(ideriv_diag) = res(ideriv_diag) + tlc%c_therm(ilay) * sum(matmul(prefactor * tlc%diff_l_sp_p_tau(:,ilay) , fld(ist_therm,:,idir_bck,isplit_start:isplit_end-1) + fld(ist_therm,:,idir_fwd,isplit_start+1:isplit_end)))

      endif

      ! Progress to next layer.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = t_tot_x at bottom
    ! L = 1
    ! C = Cb_srf

    ! Add the relative derivative of the transmission term.
    if (exist_bdrf) then
      do istr = 1,nstrhalf
        res = res + reldiff_t_tot_x_tau(ilay_deriv,nlay) * t_tot_x(nlay) * prefactor(istr) * sum(c_bck_srf_x(:,istr) * fld(:,istr,idir_bck,nsplit_total))
      enddo
    endif

  end subroutine perturbation_tau_source_response_1_sp ! }}}

  !> Calculates the inner products with the derivative of a first order source or response with layer spltting and coupling to interpolated field.
  !!
  !! This routine calculates the derivative of source_response_1_int from lintran_transfer_module
  !! with respect to the optical depth in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  subroutine perturbation_tau_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_bdrf,exist_therm,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Adminitration for diagonal terms, of which we do not know where they are.
    logical :: diag ! Diagonal term is included in result.
    integer :: ideriv_diag ! Index of diagonal term in result.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    integer, dimension(:), pointer :: ideriv_sensitive_x
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:), pointer :: l_int_sm_x
    real, dimension(:), pointer :: l_int_ac_x
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: reldiff_t_tot_x_tau
    real, dimension(:,:), pointer :: reldiff_t_sp_x_tau
    real, dimension(:,:), pointer :: diff_l_int_sm_x_tau
    real, dimension(:,:), pointer :: diff_l_int_ac_x_tau

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: istr ! Over streams.
    integer :: isplit ! Over sublayers.
    integer :: ideriv ! Over perturbed layers.

    ! Aggregate of TLC-parameters relevant for one sublayer. One for differentiating
    ! with respect to T, and one for differentiating with respect to L.
    ! The one for L is more complicated, because it combines different kinds of contributions
    ! in different terms of the derivatives in L. The one of T is more dynamic, because it gets
    ! an extra dynamic factor of a changing relative derivative of T for different sublayers.
    ! Added dimension: 0 for bck and 1 for fwd.
    ! Added dimension: 0 for sm and 1 for ac.
    real, dimension(nst,nstrhalf,0:1,0:1) :: sublayer_aggregated_for_t
    real, dimension(nst,nstrhalf,0:1,0:1,nlay_deriv), target :: sublayer_aggregated_for_l

    ! Aggregate of TLC-parameters that can be re-used in loops over perturbed layers.
    ! Added dimension: 0 for bck and 1 for fwd.
    real, dimension(nst,nstrhalf,0:1) :: arr_aggregated_d

    ! Reshape pointers.
    real, dimension(:), pointer :: fld_sublayer_1d
    real, dimension(:,:), pointer :: sublayer_aggregated_for_l_2d

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr
    ! For sublayer_aggregated_for_l, we use the reshape pointer for this purpose as well.

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
        c_bck_srf_x => tlc%c_bck_srf_v_adj
      else
        c_fwd_x => tlc%c_fwd_0_nrm
        c_bck_x => tlc%c_bck_0_nrm
        c_bck_srf_x => tlc%c_bck_srf_0_nrm
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
        c_bck_srf_x => tlc%c_bck_srf_0_adj
      else
        c_fwd_x => tlc%c_fwd_v_nrm
        c_bck_x => tlc%c_bck_v_nrm
        c_bck_srf_x => tlc%c_bck_srf_v_nrm
      endif

    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Special direction is the instrument (v).
      ideriv_sensitive_x => tlc%ideriv_sensitive_v
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_int_sm_x => tlc%l_int_sm_v
      l_int_ac_x => tlc%l_int_ac_v
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_v_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_v_tau
      diff_l_int_sm_x_tau => tlc%diff_l_int_sm_v_tau
      diff_l_int_ac_x_tau => tlc%diff_l_int_ac_v_tau
    else
      ! Special direction is the sun (0).
      ideriv_sensitive_x => tlc%ideriv_sensitive_0
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_int_sm_x => tlc%l_int_sm_0
      l_int_ac_x => tlc%l_int_ac_0
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_0_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_0_tau
      diff_l_int_sm_x_tau => tlc%diff_l_int_sm_0_tau
      diff_l_int_ac_x_tau => tlc%diff_l_int_ac_0_tau
    endif

    isplit_start = 0
    ideriv_diag = 0 ! Only relevant for plane-parallel mode and for thermal emission.
    do ilay = 1,nlay

      ! Acquire end sublayer index.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! Set sensitivity pointers.
      res_ptr => res(1:ideriv_sensitive_x(ilay))
      ilay_deriv_ptr => ilay_deriv(1:ideriv_sensitive_x(ilay))

      ! TLC-parameters for scattering in a layer:
      ! T = T_tot_x at interface above current sublayer
      ! L = L_int_x
      ! C = Cb_x or Cf_x

      do istr = 1,nstrhalf

        arr_aggregated_d(:,istr,0) = prefactor(istr) * t_tot_x(ilay-1) * c_bck_x(:,istr,ilay) ! For backward scattering
        arr_aggregated_d(:,istr,1) = prefactor(istr) * t_tot_x(ilay-1) * c_fwd_x(:,istr,ilay) ! For forward scattering

      enddo

      ! Construct aggregate for T.
      ! The last dimension index is for the interpolation contribution, 0=sm, 1=ac.
      ! This also translated the 0=bck and 1=fwd to idir_bck and idir_fwd, which may
      ! be the other way around. For sm/ac, it always good that sm=0 and ac=1.
      sublayer_aggregated_for_t(:,:,idir_bck,0) = arr_aggregated_d(:,:,0) * l_int_sm_x(ilay)
      ! Other possibilities sm/ac :,and fwd/back for C_x.
      sublayer_aggregated_for_t(:,:,idir_bck,1) = arr_aggregated_d(:,:,0) * l_int_ac_x(ilay)
      sublayer_aggregated_for_t(:,:,idir_fwd,0) = arr_aggregated_d(:,:,1) * l_int_sm_x(ilay)
      sublayer_aggregated_for_t(:,:,idir_fwd,1) = arr_aggregated_d(:,:,1) * l_int_ac_x(ilay)

      ! Apply perturbation of transmission term, also nondiagonal for plane-parallel case.
      do isplit = isplit_start,isplit_end-1

        res_ptr = res_ptr + t_sp_x(ilay)**(isplit-isplit_start) * (reldiff_t_tot_x_tau(ilay_deriv_ptr,ilay-1) + (isplit-isplit_start) * reldiff_t_sp_x_tau(ilay_deriv_ptr,ilay)) * sum(sublayer_aggregated_for_t * fld(:,:,:,isplit:isplit+1))

      enddo

      ! Check for diagonal term only necessary for plane parallel mode, and for
      ! thermal emission. because there is no part of the derivative of L that is
      ! always only diagonal, except for the one for thermal emission.
      if (tlc%plane_parallel .or. exist_therm) then
        diag = .false.
        if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
          if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
            diag = .true.
            ideriv_diag = ideriv_diag + 1
          endif
        endif
      endif

      ! For plane-parallel mode, only the diagonal term is relevant.
      if (tlc%plane_parallel .and. .not. diag) then
        ! Progress to next layer and get out of the loop.
        isplit_start = isplit_end
        cycle
      endif

      if (tlc%plane_parallel) then
        ! Directly set loop limits for ideriv.
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
      else
        ! Set loop limits to all possibly sensitive layers.
        ideriv_start = 1
        ideriv_end = ideriv_sensitive_x(ilay)
      endif

      ! Reset sensitivity pointers, including the reshape pointer.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)
      sublayer_aggregated_for_l_2d(1:4*nst*nstrhalf,ideriv_start:ideriv_end) => sublayer_aggregated_for_l(:,:,:,:,ideriv_start:ideriv_end) ! Also a reshape pointer.

      ! Loop over perturbed layers for the sublayer results for L.
      do ideriv = ideriv_start,ideriv_end

        ! Calculate sublayer terms for L.
        ! They are all nondiagonal.
        ! Also, the bck/fwd dimension is translated and the sm/ac dimension is left
        ! as it is, similarly to the aggregated sublayer for T.
        sublayer_aggregated_for_l(:,:,idir_bck,0,ideriv) = arr_aggregated_d(:,:,0) * diff_l_int_sm_x_tau(ilay_deriv(ideriv),ilay)
        ! Other possibilities sm/ac and fwd/bck for C_x
        sublayer_aggregated_for_l(:,:,idir_bck,1,ideriv) = arr_aggregated_d(:,:,0) * diff_l_int_ac_x_tau(ilay_deriv(ideriv),ilay)
        sublayer_aggregated_for_l(:,:,idir_fwd,0,ideriv) = arr_aggregated_d(:,:,1) * diff_l_int_sm_x_tau(ilay_deriv(ideriv),ilay)
        sublayer_aggregated_for_l(:,:,idir_fwd,1,ideriv) = arr_aggregated_d(:,:,1) * diff_l_int_ac_x_tau(ilay_deriv(ideriv),ilay)

      enddo

      do isplit = isplit_start,isplit_end-1

        ! Set reshape pointer to the field relevant for this sublayer. Fot T, no reshape
        ! pointer is needed, because it is no matmul. Summations work on multidimensional
        ! arrays, matmuls not. The 4 in the dimension folds up the direction and layer indices.
        fld_sublayer_1d(1:4*nst*nstrhalf) => fld(:,:,:,isplit:isplit+1)

        ! Get results.
        res_ptr = res_ptr + t_sp_x(ilay)**(isplit-isplit_start) * matmul(fld_sublayer_1d , sublayer_aggregated_for_l_2d)

      enddo

      ! Diagonal contribution from thermal emission. It is separated because it does
      ! not have the same transmission term.
      if (exist_therm .and. diag) then

        ! TLC-parameters for thermal emission:
        ! T = 1
        ! L = L_int
        ! C = C_T

        ! Differentiation of L_int is only diagonal.
        ! The matmul zips the streams.
        ! And Stokes parameters only exist in the fields and that
        ! must obey the Stokes parameter for thermal emission.
        ! The left sum sums upward and downward.
        ! The right sum sums split layers.
        res(ideriv_diag) = res(ideriv_diag) + tlc%diff_l_int_tau(ilay) * tlc%c_therm(ilay) * sum(matmul(prefactor , fld(ist_therm,:,:,isplit_start) + 2.0 * sum(fld(ist_therm,:,:,isplit_start+1:isplit_end-1),3) + fld(ist_therm,:,:,isplit_end)))

      endif

      ! Progress to next layer.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = t_tot_x at bottom
    ! L = 1
    ! C = Cb_srf

    ! Add the relative derivative of the transmission term.
    if (exist_bdrf) then
      do istr = 1,nstrhalf
        res = res + reldiff_t_tot_x_tau(ilay_deriv,nlay) * prefactor(istr) * t_tot_x(nlay) * sum(c_bck_srf_x(:,istr) * fld(:,istr,idir_bck,nsplit_total))
      enddo
    endif

  end subroutine perturbation_tau_source_response_1_int ! }}}

  !> Calculates the inner products with the derivative of a second order source or response with layer spltting and coupling to interpolated field.
  !!
  !! The corresponding non-derivative source_response_2_int is not needed and therefore does
  !! not exist in lintran_transfer_module. It is a source or response function where two
  !! scatter events in one layer are involved. One scatter event couples to an interpolated
  !! field and the other scatter events involves the solar or viewing direction. Directly, the
  !! inner product with an existing field is taken, so that derivatives for all layers can be
  !! handled at once.
  subroutine perturbation_tau_source_response_2_int(nstrhalf,nlay,nsplit_total,nlay_deriv,ilay_deriv,flag_source,flag_adjoint,exist_therm,fld,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_i, ist_u, ist_v, ist_therm, fullpol_nst
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: fld !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in), target :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout), target :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    integer, dimension(:), pointer :: ideriv_sensitive_x
    real, dimension(:), pointer :: prefactor
    real, dimension(:,:), pointer :: effmu_px
    real, dimension(:,:), pointer :: effmu_mx
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:), pointer :: l_int_sm_x
    real, dimension(:), pointer :: l_int_ac_x
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:,:), pointer :: diff_effmu_px_tau
    real, dimension(:,:,:), pointer :: diff_effmu_mx_tau
    real, dimension(:,:), pointer :: reldiff_t_tot_x_tau
    real, dimension(:,:), pointer :: reldiff_t_sp_x_tau
    real, dimension(:,:), pointer :: diff_l_int_sm_x_tau
    real, dimension(:,:), pointer :: diff_l_int_ac_x_tau

    ! Temporary TLC-parameters, because C-parameters are asymmetric due to V-polarization.
    ! These are only needed if V-polarization is involved.
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_fwd_id
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_bck_id

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Loop for ideriv. Must be from one to the sensitivity limit for pseudo-spherical
    ! case. Just ideriv_diag should be taken for plane-parallel mode.
    integer :: ideriv_start
    integer :: ideriv_end

    ! Iterators.
    integer :: ilay ! Atmospheric layers.
    integer :: ideriv ! Over perturbed layers.
    integer :: isplit ! Over sub-layers in one layer.
    integer :: istr ! Over streams, mostly destination
    integer :: istr_src ! Over source streams

    ! Pointers to sensitive parts.
    real, dimension(:), pointer :: res_ptr
    integer, dimension(:), pointer :: ilay_deriv_ptr
    ! For sublayer_aggregated_for_l, we use the reshape pointer for this purpose as well.

    ! Adminitration for diagonal terms, of which we do not know where they are.
    integer :: ideriv_diag ! Index of diagonal term in result.
    logical :: diag ! Diagonal term is included in result.

    ! Aggregate of TLC-parameters relevant for one sublayer. One for differentiating
    ! with respect to T, and one for differentiating with respect to L.
    ! The one for L is more complicated, because it combines different kinds of contributions
    ! in different terms of the derivatives in L. The one of T is more dynamic, because it gets
    ! an extra dynamic factor of a changing relative derivative of T for different sublayers.
    real, dimension(nst,nstrhalf,0:1,0:1) :: sublayer_aggregated_for_t
    real, dimension(nst,nstrhalf,0:1,0:1,nlay_deriv), target :: sublayer_aggregated_for_l

    ! Aggregates of TLC-parameters for re-use.
    real, dimension(nst,nstrhalf,0:1) :: arr_aggregated_i
    ! Further `advanced' constructions, because it is a complicated function.
    real, dimension(nst,nstrhalf,0:1) :: arr_aggregated_i_advanced_c1
    real, dimension(nst,nstrhalf,0:1) :: arr_aggregated_i_advanced_c2
    real, dimension(nst,nstrhalf,0:1) :: arr_aggregated_d_advanced
    real, dimension(nst,nstrhalf,0:1,0:1) :: arr_aggregated_d_advanced_4d
    real, dimension(nst,nstrhalf) :: arr_aggregated_d_therm ! 2D for thermal emission.

    ! Reshape pointers.
    real, dimension(:), pointer :: fld_sublayer_1d
    real, dimension(:,:), pointer :: sublayer_aggregated_for_l_2d

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
      ideriv_sensitive_x => tlc%ideriv_sensitive_v
      effmu_px => tlc%effmu_pv
      effmu_mx => tlc%effmu_mv
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_int_sm_x => tlc%l_int_sm_v
      l_int_ac_x => tlc%l_int_ac_v
      diff_effmu_px_tau => tlc%diff_effmu_pv_tau
      diff_effmu_mx_tau => tlc%diff_effmu_mv_tau
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_v_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_v_tau
      diff_l_int_sm_x_tau => tlc%diff_l_int_sm_v_tau
      diff_l_int_ac_x_tau => tlc%diff_l_int_ac_v_tau
    else
      ! Special direction is the sun (0).
      ideriv_sensitive_x => tlc%ideriv_sensitive_0
      effmu_px => tlc%effmu_p0
      effmu_mx => tlc%effmu_m0
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_int_sm_x => tlc%l_int_sm_0
      l_int_ac_x => tlc%l_int_ac_0
      diff_effmu_px_tau => tlc%diff_effmu_p0_tau
      diff_effmu_mx_tau => tlc%diff_effmu_m0_tau
      reldiff_t_tot_x_tau => tlc%reldiff_t_tot_0_tau
      reldiff_t_sp_x_tau => tlc%reldiff_t_sp_0_tau
      diff_l_int_sm_x_tau => tlc%diff_l_int_sm_0_tau
      diff_l_int_ac_x_tau => tlc%diff_l_int_ac_0_tau
    endif

    ideriv_diag = 0
    isplit_start = 0
    do ilay = 1,nlay

      ! Acquire end sublayer index.
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

      ! Set sensitivity pointers.
      res_ptr => res(1:ideriv_sensitive_x(ilay))
      ilay_deriv_ptr => ilay_deriv(1:ideriv_sensitive_x(ilay))

      ! The first scatter event is the analytic one and
      ! the second one couples to an interpolated field.

      ! TLC-parameters for (forward - ?):
      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_mx(i) * (L_int_p(i) - L_int_x), with interpolations
      ! C = Cf_x(i) * C?(i,d)

      ! Scalar      : T_tot_x, L_int_x
      ! Array in i  : prefactor_intermediate, L_int_p, effmu_mx, Cf_x
      ! Array in d  : prefactor
      ! Matrix(i,d) : Cb

      ! TLC-parameters for (back - ?):
      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_px(i) * (L_int_x - T_px(i) * L_int_m(i))
      ! C = Cb_x(i) * C?(i,d)

      ! Scalar      : T_tot_x, L_int_x
      ! Array in i  : prefactor_intermediate, L_int_m, T_px, effmu_px, Cb_x
      ! Array in d  : prefactor
      ! Matrix(i,d) : Cf

      ! Differentiations for (forward - ?):
      ! - Transmission, with the entire undifferentiated field
      ! - effmu_mx (with L_int_p - L_int_x)
      ! - L_int_p (with effmu_mx), diagonal
      ! - L_int_x (with -effmu_mx)

      ! Differentiations for (back - ?):
      ! - effmu_px (with L_int_x - T_px * L_int_m)
      ! - L_int_x (with effmu_px)
      ! - T_px (with -effmu_px * T_px * L_int_m), diagonal for p, normal for x
      ! - L_int_m (with -effmu_px * T_px), diagonal

      ! Construct aggregate of TLC-parameters for re-use.
      do istr = 1,nstrhalf
        arr_aggregated_i(:,istr,0) = t_tot_x(ilay-1) * tlc%prefactor_intermediate(istr) * c_bck_x(:,istr,ilay) ! With backward scattering.
        arr_aggregated_i(:,istr,1) = t_tot_x(ilay-1) * tlc%prefactor_intermediate(istr) * c_fwd_x(:,istr,ilay) ! With forward scattering.
      enddo

      ! Construct aggregate for T.
      ! The last dimension index is for the interpolation contribution, 0=sm, 1=ac.
      ! res = scalar * array_d * matmul(array_i,matrix_id)
      sublayer_aggregated_for_t = 0.0 ! Because we are summing in a loop over istr_src.
      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! `Source', this is direction i.
          sublayer_aggregated_for_t(:,istr,idir_bck,0) = sublayer_aggregated_for_t(:,istr,idir_bck,0) + prefactor(istr) * (effmu_mx(istr_src,ilay) * (tlc%l_int_sm_p(istr_src,ilay) - l_int_sm_x(ilay)) * matmul(arr_aggregated_i(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) + effmu_px(istr_src,ilay) * (l_int_sm_x(ilay) - tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay)) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)))
          ! Other possibilities sm/ac and fwd/back for C_id.
          sublayer_aggregated_for_t(:,istr,idir_bck,1) = sublayer_aggregated_for_t(:,istr,idir_bck,1) + prefactor(istr) * (effmu_mx(istr_src,ilay) * (tlc%l_int_ac_p(istr_src,ilay) - l_int_ac_x(ilay)) * matmul(arr_aggregated_i(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) + effmu_px(istr_src,ilay) * (l_int_ac_x(ilay) - tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay)) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)))
          sublayer_aggregated_for_t(:,istr,idir_fwd,0) = sublayer_aggregated_for_t(:,istr,idir_fwd,0) + prefactor(istr) * (effmu_mx(istr_src,ilay) * (tlc%l_int_sm_p(istr_src,ilay) - l_int_sm_x(ilay)) * matmul(arr_aggregated_i(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) + effmu_px(istr_src,ilay) * (l_int_sm_x(ilay) - tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay)) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)))
          sublayer_aggregated_for_t(:,istr,idir_fwd,1) = sublayer_aggregated_for_t(:,istr,idir_fwd,1) + prefactor(istr) * (effmu_mx(istr_src,ilay) * (tlc%l_int_ac_p(istr_src,ilay) - l_int_ac_x(ilay)) * matmul(arr_aggregated_i(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) + effmu_px(istr_src,ilay) * (l_int_ac_x(ilay) - tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay)) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)))
        enddo
      enddo

      ! Apply perturbation of transmission term.
      do isplit = isplit_start,isplit_end-1

        res_ptr = res_ptr + t_sp_x(ilay)**(isplit-isplit_start) * (reldiff_t_tot_x_tau(ilay_deriv_ptr,ilay-1) + (isplit-isplit_start) * reldiff_t_sp_x_tau(ilay_deriv_ptr,ilay)) * sum(sublayer_aggregated_for_t * fld(:,:,:,isplit:isplit+1))

      enddo

      ! Check for diagonal terms.
      diag = .false.
      if (ideriv_diag .ne. nlay_deriv) then ! To prevent out-of-bounds exception.
        if (ilay_deriv(ideriv_diag+1) .eq. ilay) then
          diag = .true.
          ideriv_diag = ideriv_diag + 1
        endif
      endif

      ! Derivatives of L, only relevant for pseudo-spherical or diagonal for plane-parallel.
      if (tlc%plane_parallel .and. .not. diag) then
        ! Progress to next layer and go out of the loop.
        isplit_start = isplit_end
        cycle
      endif

      do istr = 1,nstrhalf
        ! Aggregate case-specfic aggregates.
        ! The last dimension is associated with interpolation contributions sm/ac.
        ! Case 1: Differentiation of effmu_mx.
        arr_aggregated_i_advanced_c1(:,istr,0) = arr_aggregated_i(:,istr,1) * (tlc%l_int_sm_p(istr,ilay) - l_int_sm_x(ilay))
        arr_aggregated_i_advanced_c1(:,istr,1) = arr_aggregated_i(:,istr,1) * (tlc%l_int_ac_p(istr,ilay) - l_int_ac_x(ilay))
        ! Case 2: Differentiation of effmu_px.
        arr_aggregated_i_advanced_c2(:,istr,0) = arr_aggregated_i(:,istr,0) * (l_int_sm_x(ilay) - tlc%t_sp(istr,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr,ilay))
        arr_aggregated_i_advanced_c2(:,istr,1) = arr_aggregated_i(:,istr,0) * (l_int_ac_x(ilay) - tlc%t_sp(istr,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr,ilay))

        ! Aggregates in destination dimension. Now the last dimension is associated with involved
        ! scatter direction, 0=bck, 1=fwd.
        arr_aggregated_d_advanced(:,istr,0) = 0.0 ! Because there is summation in a loop over istr_src.
        arr_aggregated_d_advanced(:,istr,1) = 0.0 ! Because there is summation in a loop over istr_src.
        do istr_src = 1,nstrhalf
          arr_aggregated_d_advanced(:,istr,0) = arr_aggregated_d_advanced(:,istr,0) + effmu_px(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)) - effmu_mx(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) ! Net backward scattering
          arr_aggregated_d_advanced(:,istr,1) = arr_aggregated_d_advanced(:,istr,1) + effmu_px(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)) - effmu_mx(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) ! Net backward scattering
        enddo

        ! Four-dimensional aggregates. Here the third dimension is a net scatter direction,
        ! 0=bck, 1=fwd and the last dimension is an interpolation contribution, 0=sm, 1=ac.
        arr_aggregated_d_advanced_4d(:,istr,0,0) = 0.0 ! Again, there is a summation in a loop.
        arr_aggregated_d_advanced_4d(:,istr,1,0) = 0.0 ! Again, there is a summation in a loop.
        arr_aggregated_d_advanced_4d(:,istr,0,1) = 0.0 ! Again, there is a summation in a loop.
        arr_aggregated_d_advanced_4d(:,istr,1,1) = 0.0 ! Again, there is a summation in a loop.

        ! Add contributions for different source (i) streams
        do istr_src = 1,nstrhalf
          arr_aggregated_d_advanced_4d(:,istr,0,0) = arr_aggregated_d_advanced_4d(:,istr,0,0) + effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr))
          arr_aggregated_d_advanced_4d(:,istr,1,0) = arr_aggregated_d_advanced_4d(:,istr,1,0) + effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr))
          arr_aggregated_d_advanced_4d(:,istr,0,1) = arr_aggregated_d_advanced_4d(:,istr,0,1) + effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr))
          arr_aggregated_d_advanced_4d(:,istr,1,1) = arr_aggregated_d_advanced_4d(:,istr,1,1) + effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr))
        enddo

      enddo

      ! Set rules for loop over ideriv.
      if (tlc%plane_parallel) then
        ! For plane-parallel, there must be diagonal index, otherwise, we would not be here.
        ideriv_start = ideriv_diag
        ideriv_end = ideriv_diag
      else
        ideriv_start = 1
        ideriv_end = ideriv_sensitive_x(ilay)
      endif

      ! Set sensitivity pointers.
      res_ptr => res(ideriv_start:ideriv_end)
      ilay_deriv_ptr => ilay_deriv(ideriv_start:ideriv_end)
      sublayer_aggregated_for_l_2d(1:4*nst*nstrhalf,ideriv_start:ideriv_end) => sublayer_aggregated_for_l(:,:,:,:,ideriv_start:ideriv_end) ! Also a reshape pointer.

      ! Initialize the relevant information to zero, so it can be acquired through addition.
      sublayer_aggregated_for_l_2d = 0.0

      ! Loop over perturbed layers for the sublayer results for L.
      do ideriv = ideriv_start,ideriv_end

        ! Differentiations for (forward - ?)
        ! - Transmission, with the entire undifferentiated field
        ! - effmu_mx (with L_int_p - L_int_x)
        ! - L_int_p (with effmu_mx), diagonal
        ! - L_int_x (with -effmu_mx)

        ! Differentiations for (back - ?)
        ! - effmu_px (with L_int_x - T_px * L_int_m)
        ! - L_int_x (with effmu_px)
        ! - T_px (with -effmu_px * T_px * L_int_m), diagonal for p, normal for x
        ! - L_int_m (with -effmu_px * T_px), diagonal

        ! Calculate nondiagonal sublayer terms for L.
        do istr = 1,nstrhalf ! Stream d (destination).
          ! Contributions depending on source stream (only pseudo-spherical where derivatives
          ! of effmu_.x are nonzero).
          if (.not. tlc%plane_parallel) then
            do istr_src = 1,nstrhalf ! Stream i (intermediate).
              sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv) = sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv) + prefactor(istr) * (diff_effmu_mx_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c1(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)) + diff_effmu_px_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c2(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)))
              ! Other possibilities sm/ac and fwd/back for C_id.
              sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv) = sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv) + prefactor(istr) * (diff_effmu_mx_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c1(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) + diff_effmu_px_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c2(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)))
              sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv) = sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv) + prefactor(istr) * (diff_effmu_mx_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c1(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)) + diff_effmu_px_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c2(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)))
              sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv) = sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv) + prefactor(istr) * (diff_effmu_mx_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c1(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) + diff_effmu_px_tau(istr_src,ilay_deriv(ideriv),ilay) * matmul(arr_aggregated_i_advanced_c2(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)))
            enddo
          endif
          ! Contribution not coupled to the source stream.
          sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv) = sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv) + prefactor(istr) * (diff_l_int_sm_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced(:,istr,0) - reldiff_t_sp_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced_4d(:,istr,0,0))
          ! Other possibilities sm/ac and fwd/back for C_id.
          sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv) = sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv) + prefactor(istr) * (diff_l_int_ac_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced(:,istr,0) - reldiff_t_sp_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced_4d(:,istr,0,1))
          sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv) = sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv) + prefactor(istr) * (diff_l_int_sm_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced(:,istr,1) - reldiff_t_sp_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced_4d(:,istr,1,0))
          sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv) = sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv) + prefactor(istr) * (diff_l_int_ac_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced(:,istr,1) - reldiff_t_sp_x_tau(ilay_deriv(ideriv),ilay) * arr_aggregated_d_advanced_4d(:,istr,1,1))
        enddo

      enddo

      if (diag) then

        ! Diagonal terms.
        ! result = scalar * array_d * matmul(array_i,matrix_id)
        do istr = 1,nstrhalf
          do istr_src = 1,nstrhalf
            sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv_diag) = sublayer_aggregated_for_l(:,istr,idir_bck,0,ideriv_diag) + prefactor(istr) * (effmu_mx(istr_src,ilay) * tlc%diff_l_int_sm_p_tau(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) - tlc%reldiff_t_sp_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)) - tlc%diff_l_int_sm_m_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)))
            ! Other possibilities sm/ac and fwd/back for C_id
            sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv_diag) = sublayer_aggregated_for_l(:,istr,idir_bck,1,ideriv_diag) + prefactor(istr) * (effmu_mx(istr_src,ilay) * tlc%diff_l_int_ac_p_tau(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_bck_id(:,:,istr_src,istr)) - tlc%reldiff_t_sp_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)) - tlc%diff_l_int_ac_m_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_fwd_id(:,:,istr_src,istr)))
            sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv_diag) = sublayer_aggregated_for_l(:,istr,idir_fwd,0,ideriv_diag) + prefactor(istr) * (effmu_mx(istr_src,ilay) * tlc%diff_l_int_sm_p_tau(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) - tlc%reldiff_t_sp_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_sm_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)) - tlc%diff_l_int_sm_m_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)))
            sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv_diag) = sublayer_aggregated_for_l(:,istr,idir_fwd,1,ideriv_diag) + prefactor(istr) * (effmu_mx(istr_src,ilay) * tlc%diff_l_int_ac_p_tau(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,1) , c_fwd_id(:,:,istr_src,istr)) - tlc%reldiff_t_sp_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * tlc%l_int_ac_m(istr_src,ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)) - tlc%diff_l_int_ac_m_tau(istr_src,ilay) * effmu_px(istr_src,ilay) * tlc%t_sp(istr_src,ilay) * t_sp_x(ilay) * matmul(arr_aggregated_i(:,istr_src,0) , c_bck_id(:,:,istr_src,istr)))
          enddo
        enddo

      endif

      do isplit = isplit_start,isplit_end-1

        ! Set reshape pointer to the field relevant for this sublayer. Fot T, no reshape
        ! pointer is needed, because it is no matmul. Summations work on multidimensional
        ! arrays, matmuls not.
        fld_sublayer_1d(1:4*nst*nstrhalf) => fld(:,:,:,isplit:isplit+1)

        ! Get results.
        res_ptr = res_ptr + t_sp_x(ilay)**(isplit-isplit_start) * matmul(fld_sublayer_1d , sublayer_aggregated_for_l_2d)

      enddo

      ! Diagonal contribution from thermal emission. It is separated because it does
      ! not have the same transmission term.
      if (exist_therm .and. diag) then

        ! TLC-parameters for thermal emission:
        ! Emission towards the `ac'-side.
        ! T = 1
        ! L = effmu_p * (L_int - L_int_p), all at intermediate direction
        ! C = C_T * C(i,d)

        ! Scalar      : DL_int, C_T
        ! Array in i  : prefactor_intermediate, effmu_p, L_int_p
        ! Array in d  : prefactor
        ! Matrix(i,d) : C

        ! Differentiation turns the L-terms into their derivatives.

        ! result = scalar * array_layer * matmul(array_streams,matrix)

        ! The aggregate arrays are not to recycle results more often. It is just that
        ! one matmul in a line is more than enough.

        arr_aggregated_d_therm = 0.0 ! Because we will acquire this by addition.
        do istr = 1,nstrhalf
          do istr_src = 1,nstrhalf
            ! Case 1: Backward scattering `sm' side.
            ! After backward scattering, the target radiation is directed towards the `sm' side.
            ! `sm' contribution, upward radiation with low indices, downward radiation with
            ! high indices.
            arr_aggregated_d_therm(:,istr) = arr_aggregated_d_therm(:,istr) + tlc%c_therm(ilay) * prefactor(istr) * tlc%prefactor_intermediate(istr_src) * tlc%effmu_p(istr_src) * (tlc%diff_l_int_tau(ilay) - tlc%diff_l_int_sm_p_tau(istr_src,ilay)) * c_bck_id(ist_therm,:,istr_src,istr)
            ! Case 2: Forward scattering `ac' side.
            ! After forward scattering, the target radiation is directed towards the `ac' side.
            ! `ac' contribution, upward radiation with low indices, downward radiation with
            ! high indices.
            arr_aggregated_d_therm(:,istr) = arr_aggregated_d_therm(:,istr) + tlc%c_therm(ilay) * prefactor(istr) * tlc%prefactor_intermediate(istr_src) * tlc%effmu_p(istr_src) * (tlc%diff_l_int_tau(ilay) - tlc%diff_l_int_ac_p_tau(istr_src,ilay)) * c_fwd_id(ist_therm,:,istr_src,istr)
          enddo
        enddo

        ! As they are all for the same indices of the field, aggregate the results.
        res(ideriv_diag) = res(ideriv_diag) + sum(arr_aggregated_d_therm * sum(fld(:,:,iup,isplit_start:isplit_end-1) + fld(:,:,idn,isplit_start+1:isplit_end),3))

        ! Reset this array for new use.
        arr_aggregated_d_therm = 0.0 ! Because we will acquire this by addition.

        do istr = 1,nstrhalf
          do istr_src = 1,nstrhalf
            ! Case 3: Backward scattering `ac' side.
            ! After backward scattering, the target radiation is directed towards the `sm' side.
            ! `ac' contribution, upward radiation with high indices, downward radiation with
            ! low indices.
            arr_aggregated_d_therm(:,istr) = arr_aggregated_d_therm(:,istr) + tlc%c_therm(ilay) * prefactor(istr) * tlc%prefactor_intermediate(istr_src) * tlc%effmu_p(istr_src) * (tlc%diff_l_int_tau(ilay) - tlc%diff_l_int_ac_p_tau(istr_src,ilay)) * c_bck_id(ist_therm,:,istr_src,istr)
            ! Case 4: Forward scattering `sm' side.
            ! After forward scattering, the target radiation is directed towards the `ac' side.
            ! `sm' contribution, upward radiation with high indices, downward radiation with
            ! low indices.
            arr_aggregated_d_therm(:,istr) = arr_aggregated_d_therm(:,istr) + tlc%c_therm(ilay) * prefactor(istr) * tlc%prefactor_intermediate(istr_src) * tlc%effmu_p(istr_src) * (tlc%diff_l_int_tau(ilay) - tlc%diff_l_int_sm_p_tau(istr_src,ilay)) * c_fwd_id(ist_therm,:,istr_src,istr)
          enddo
        enddo

        ! As they are all for the same indices of the field, aggregate the results.
        res(ideriv_diag) = res(ideriv_diag) + sum(arr_aggregated_d_therm * sum(fld(:,:,iup,isplit_start+1:isplit_end) + fld(:,:,idn,isplit_start:isplit_end-1),3))

      endif

      ! Progress to the next sublayer.
      isplit_start = isplit_end

    enddo

    ! There is no surface reflection, because this is second order.

  end subroutine perturbation_tau_source_response_2_int ! }}}

  !> Calculates the inner products with the derivative of an interpolated internal scatter event.
  !!
  !! This routine calculates the derivative of scatter_1_int from lintran_transfer_module
  !! with respect to the optical depth in the selected layers. Directly, the inner product
  !! with an existing field is taken, so that derivatives for all layers can be handled at once.
  !! As scatter_1_int also needs an input field, this routine needs both an input field for
  !! the scatter event and a field with which the inner product is taken.
  subroutine perturbation_tau_scatter_1_int(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,source,destination,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: source !< Source for the scatter event.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: destination !< Field with which the inner product is taken.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: isplit ! Over sublayers.
    integer :: ideriv ! Over differentiations.
    integer :: istr ! Over streams (destination).
    integer :: istr_src ! Over source streams.

    ! Non-iterated derived index.
    integer :: ilay ! Atmospheric layers.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! TLC-parameters for backward scattering:
    ! T = 1
    ! L = L_int_p
    ! C = Cb

    ! The prefactor is prefactor_scat.

    ! Differentiating turns L into her derivative for ilay = ilay_pert.

    ! For forward scattering, only C changes to Cf.

    ! Array in s  : source
    ! Array in d  : L_int_p, destination
    ! Matrix(s,d) : prefactor_scat, C

    ! result = array_d*matmul(array_s,matrix_sd)

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

      ! Loop over sublayers.
      do isplit = isplit_start+1,isplit_end

        do istr = 1,nstrhalf ! Destination stream.
          do istr_src = 1,nstrhalf ! Source stream.
            res(ideriv) = res(ideriv) + tlc%prefactor_scat(istr_src,istr) * (tlc%diff_l_int_sm_p_tau(istr,ilay) * (sum(destination(:,istr,iup,isplit-1) * (matmul(source(:,istr_src,idn,isplit-1) , tlc%c_bck(:,:,istr_src,istr,ilay)) + matmul(source(:,istr_src,iup,isplit-1) , tlc%c_fwd(:,:,istr_src,istr,ilay)))) + sum(destination(:,istr,idn,isplit) * (matmul(source(:,istr_src,iup,isplit) , tlc%c_bck(:,:,istr_src,istr,ilay)) + matmul(source(:,istr_src,idn,isplit) , tlc%c_fwd(:,:,istr_src,istr,ilay))))) + tlc%diff_l_int_ac_p_tau(istr,ilay) * (sum(destination(:,istr,iup,isplit-1) * (matmul(source(:,istr_src,idn,isplit) , tlc%c_bck(:,:,istr_src,istr,ilay)) + matmul(source(:,istr_src,iup,isplit) , tlc%c_fwd(:,:,istr_src,istr,ilay)))) + sum(destination(:,istr,idn,isplit) * (matmul(source(:,istr_src,iup,isplit-1) , tlc%c_bck(:,:,istr_src,istr,ilay)) + matmul(source(:,istr_src,idn,isplit-1) , tlc%c_fwd(:,:,istr_src,istr,ilay))))))
          enddo
        enddo

      enddo

    enddo

  end subroutine perturbation_tau_scatter_1_int ! }}}

  !> Takes the inner product of a normal and adjoint intensity field with a pertrubation of the radiative transfer operator.
  !!
  !! Perturbs the extinction coefficient in radiative transfer operator and performs the inner
  !! product of the adjoint field with the perturbation of the radiative transfer operator acting
  !! on the normal field. Corrections for the pseudo-forward method for the adjoint field are
  !! applied so that the inner product is correct. Both fields are interpolated. Note that the
  !! radiative transfer operator is perturbed, not the matrix. Though transmission is about
  !! moving vertically, the perturbation acts on an infinitesimal part of the atmosphere, so
  !! after application of the perturbation of the radiative transfer equation, there is no
  !! movement, so the interpoalted field remains interpolated. A chain rule coefficient converts
  !! derivatives with respect to the extinction coefficient into derivatives with respect to
  !! the optical depth.
  subroutine perturbation_abs_zip(nstrhalf,nsplit_total,nlay_deriv,ilay_deriv,intens,intens_adj,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, flip_v
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    integer, intent(in) :: nlay_deriv !< Number of layers for which derivatives are calculated.
    integer, dimension(nlay_deriv), intent(in) :: ilay_deriv !< Indices of layers for which derivatives are calculated.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens !< Involved normal intensity field.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_adj !< Involved adjoint intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nlay_deriv), intent(inout) :: res !< Resulting derivative with respect to the optical depth. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ideriv ! Over differentiations.
    integer :: istr ! Over streams.

    ! Index that is not iterated.
    integer :: ilay ! Atmospheric layers.

    ! The derivative of the radiative transfer operator is the derivative
    ! of the transmission operator, transmitting only an infinitesimal
    ! distance. There is a prefactor for zipping.

    ! We have to care for the fact that the adjoint field is still in
    ! the reverse direction.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! TLC-parameters for zipping:
    ! T = 1
    ! L = l_zip
    ! C = 1

    ! Note that the zipping L-parameters already corrects for that we are not
    ! integrating over a scatter event.

    ! Differentiation adds reldiff_t_tau, because the 1 in T is actually
    ! an infinitesimal transmission term, so it has a derivative.

    ! Actually, the reldiff_t_tau compensates for a mu in the zipping prefactor
    ! and also creates the minus sign that should be here.

    ! Remember that upward and downward are switched in the adjoint.

    ! Scalar : l_zip
    ! Array  : prefactor_zip, reldiff_t_tau, intens, intens_adj

    ! For each sublayer, there will be a contribution with l_zip_sm for
    ! intens(isplit-1)*intesn_adj(isplit-1) + intens(isplit)*intens_adj(isplit).

    ! That is once the border interfaces and twice the internal interfaces.

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

      ! Contribution for same-interface intensities.
      ! The V-flipper is always four Stokes parameters long. If there are fewer Stokes
      ! parameters, it is not needed, but if we do flip_v(1:nst), then it will do nothing
      ! for less than four Stokes parameters and do the right thing for four Stokes parameters.
      do istr = 1,nstrhalf
        res(ideriv) = res(ideriv) + tlc%prefactor_zip(istr) * tlc%l_zip_sm(ilay) * tlc%reldiff_t_tau(istr) * (sum(intens(:,istr,idn,isplit_start) * flip_v(1:nst) * intens_adj(:,istr,iup,isplit_start)) + sum(intens(:,istr,iup,isplit_start) * flip_v(1:nst) * intens_adj(:,istr,idn,isplit_start)) + sum(intens(:,istr,idn,isplit_end) * flip_v(1:nst) * intens_adj(:,istr,iup,isplit_end)) + sum(intens(:,istr,iup,isplit_end) * flip_v(1:nst) * intens_adj(:,istr,idn,isplit_end)) + 2.0 * sum(matmul(flip_v(1:nst) , intens(:,istr,idn,isplit_start+1:isplit_end-1) * intens_adj(:,istr,iup,isplit_start+1:isplit_end-1))) + 2.0 * sum(matmul(flip_v(1:nst) , intens(:,istr,iup,isplit_start+1:isplit_end-1) * intens_adj(:,istr,idn,isplit_start+1:isplit_end-1))))

        ! Contribution for other-interface intensities
        res(ideriv) = res(ideriv) + tlc%prefactor_zip(istr) * tlc%l_zip_ac(ilay) * tlc%reldiff_t_tau(istr) * (sum(matmul(flip_v(1:nst) , intens(:,istr,idn,isplit_start:isplit_end-1) * intens_adj(:,istr,iup,isplit_start+1:isplit_end))) + sum(matmul(flip_v(1:nst) , intens(:,istr,iup,isplit_start:isplit_end-1) * intens_adj(:,istr,idn,isplit_start+1:isplit_end))) + sum(matmul(flip_v(1:nst) , intens(:,istr,idn,isplit_start+1:isplit_end) * intens_adj(:,istr,iup,isplit_start:isplit_end-1))) + sum(matmul(flip_v(1:nst) , intens(:,istr,iup,isplit_start+1:isplit_end) * intens_adj(:,istr,idn,isplit_start:isplit_end-1))))
      enddo

    enddo

  end subroutine perturbation_abs_zip ! }}}
