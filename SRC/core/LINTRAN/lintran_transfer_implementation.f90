!> \file lintran_transfer_implementation.f90
!! Implementation file for unperturbed radiative transfer functions.
!! This file contains all events that can happen to a photon from emission, transmission,
!! scattering and being detectred by a detector.

! The routines and their comments are indented one level, because you think you are inside
! a module.

  !> Direct single-scattering.
  !!
  !! This is a trajectory from the sun, via one scatter event directly to the instrument.
  subroutine direct_1_elem(nlay,exist_therm,tlc,res) ! {{{

    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity. As this is direct single-scattering, this can only be if the viewing Stokes parameter is the one in which thermal emission takes place.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    ! TLC-parameters for backward scattering:
    ! T = T_tot_0v at interface ilay-1
    ! L = L_0v
    ! C = Cb_0v_elem

    ! TLC-parameters for surface reflection:
    ! T = T_tot_0v at bottom
    ! L = 1
    ! C = Cb_srf_0v_elem

    ! TLC-parameters for fluorescent emission:
    ! T = T_tot_v at bottom
    ! L = 1
    ! C = Ce_v_elem

    ! TLC-parameters for thermal emission:
    ! T = T_v at interface ilay-1
    ! L = L_v
    ! C = C_T

    ! The first term is the atmosphere, the second is surface reflection
    ! and the final term is surface emission.
    res = res + tlc%prefactor_dir * (sum(tlc%t_tot_0v(0:nlay-1) * tlc%l_0v * tlc%c_bck_0v_elem) + tlc%t_tot_0v(nlay) * tlc%c_bck_srf_0v_elem + tlc%t_tot_v(nlay) * tlc%c_emi_v_elem)

    ! Add direct thermal emission if desired.
    if (exist_therm) res = res + tlc%prefactor_dir * sum(tlc%t_tot_v(0:nlay-1) * tlc%l_v * tlc%c_therm)

  end subroutine direct_1_elem ! }}}

  !> Direct double-scattering.
  !!
  !! This is a trajectory from the sun, scattering twice within the same atmospheric layer
  !! ending up in the instrument. Because double scattering within one layer is integrated
  !! at once, there is no intermediate field. Therefore, this is called direct scattering.
  subroutine direct_2(nstrhalf,nlay,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: ist_therm
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    logical, intent(in) :: exist_therm !< Flag for existence of thermal emissivity.
    real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.

    ! TLC-parameters for (back - forward):
    ! T = T_tot_0v at interface ilay-1
    ! L = effmu_mv * (L_p0 - L_0v)
    ! C = Cb_0 * Cf_v

    ! TLC-parameters for (forward - back):
    ! T = T_tot_0v at interface ilay-1
    ! L = effmu_m0 * (L_pv - L_0v)
    ! C = Cf_0 * Cb_v

    ! TLC-parameters for (thermal - back):
    ! T = T_v
    ! L = effmu_p * (L_v - L_pv)
    ! C = C_T * Cb_v

    ! TLC-parameters for (thermal - forward):
    ! T = T_v
    ! L = effmu_mv * (L_p - L_v)
    ! C = C_T * Cf_v

    do ilay = 1,nlay

      do istr = 1,nstrhalf ! Intermediate stream.

        ! The first term is (back - forward), the second term is (forward - back).
        ! Shared TLC-parameters are put outside the parentheses.
        ! Array operation: scalar * matmul(array_i,matrix_id), where i and d
        ! refer to intermediate and destination Stokes parameters
        res = res + tlc%prefactor_dir * tlc%prefactor_intermediate(istr) * tlc%t_tot_0v(ilay-1) * (tlc%effmu_mv(istr,ilay) * (tlc%l_p0(istr,ilay) - tlc%l_0v(ilay)) * dot_product(tlc%c_bck_0_nrm(:,istr,ilay) , tlc%c_fwd_v_nrm(:,istr,ilay)) + tlc%effmu_m0(istr,ilay) * (tlc%l_pv(istr,ilay) - tlc%l_0v(ilay)) * dot_product(tlc%c_fwd_0_nrm(:,istr,ilay) , tlc%c_bck_v_nrm(:,istr,ilay)))

      enddo

      ! Thermal emissivity has no array operation over intermediate Stokes parameters,
      ! Array operation is over intermediate streams. The intermediate Stokes parameter
      ! is I, because it is thermal emission.
      ! The first term is (thermal - back), the second term is (thermal - forward).
      if (exist_therm) res = res + tlc%prefactor_dir * tlc%t_tot_v(ilay-1) * tlc%c_therm(ilay) * sum(tlc%prefactor_intermediate * (tlc%effmu_p * (tlc%l_v(ilay) - tlc%l_pv(:,ilay)) * tlc%c_bck_v_nrm(ist_therm,:,ilay) + tlc%effmu_mv(:,ilay) * (tlc%l_p(:,ilay) - tlc%l_v(ilay)) * tlc%c_fwd_v_nrm(ist_therm,:,ilay)))

    enddo

  end subroutine direct_2 ! }}}

  !> Transmission through zero or more layers in original vertical grid.
  !!
  !! This routine lets radiation be transmitted from one interface to another, without
  !! changing directions. The next scatter operator will assume that the radiation will
  !! scatter in the next layer. All zero or more layers completely cross are handled by
  !! this routine. No layer splitting is applied, so it can only be applied for analytic
  !! (limited-scattering) calculations.
  subroutine transmit(nstrhalf,nlay,flag_reverse,intens_in,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    logical, intent(in) :: flag_reverse !< True to reverse transmission directions, useful when the input is a response function.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(in) :: intens_in !< Input intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(out) :: res !< The resulting field.

    ! Direction indices.
    integer :: idir_dn ! What the routine thinks is downwards.
    integer :: idir_up ! What the routine thinks is upwards.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.

    ! Temporary variable for the intensity that has been transmitted from
    ! previous layer up to the current level.
    real, dimension(nst,nstrhalf) :: transmitted

    ! Reverse indices when necessary.
    if (flag_reverse) then
      idir_dn = iup
      idir_up = idn
    else
      idir_dn = idn
      idir_up = iup
    endif

    ! Downward radiation, transmit from top of the atmosphere down to the
    ! surface. High in the atmosphere, layer indices are low, so the layer
    ! iterator goes from low to high.

    ! During the loop, we will first progress the current intensity through
    ! the next layer and then add the intensity from the input at the
    ! bottom of that layer.

    ! The source at the first layer has to be initialized.
    transmitted = intens_in(:,:,idir_dn,0)
    res(:,:,idir_dn,0) = transmitted

    do ilay = 1,nlay

      ! Apply extinction by layer ilay (between levels ilay-1 and ilay).
      do istr = 1,nstrhalf
        transmitted(:,istr) = transmitted(:,istr) * tlc%t(istr,ilay) ! Array operation over Stokes parameters.
      enddo

      ! Add source of level ilay (below the layer where we went through).
      transmitted = transmitted + intens_in(:,:,idir_dn,ilay)

      ! That is the result at level ilay.
      res(:,:,idir_dn,ilay) = transmitted

    enddo

    ! Now repeat the process in reverse direction for upward direction.
    transmitted = intens_in(:,:,idir_up,nlay)
    res(:,:,idir_up,nlay) = transmitted

    do ilay = nlay,1,-1

      ! Apply extinction by layer ilay (between levels ilay and ilay-1).
      do istr = 1,nstrhalf
        transmitted(:,istr) = transmitted(:,istr) * tlc%t(istr,ilay) ! Array operation over Stokes parameters.
      enddo

      ! Add source of level ilay-1 (above the layer where we went through).
      transmitted = transmitted + intens_in(:,:,idir_up,ilay-1)

      ! That is the output at level ilay-1.
      res(:,:,idir_up,ilay-1) = transmitted

    enddo

  end subroutine transmit ! }}}

  !> Transmission through zero or more layers in split vertical grid.
  !!
  !! This routine lets radiation be transmitted from one interface to another, without
  !! changing directions. The next scatter operator will assume that the radiation will
  !! scatter in the next sublayer. All zero or more sublayers completely cross are handled by
  !! this routine. Layer splitting is applied, so it can be applied for multi-scattering.
  subroutine transmit_sp(nstrhalf,nlay,nsplit_total,flag_reverse,intens_in,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    logical, intent(in) :: flag_reverse !< True to reverse transmission directions, useful when the input is a response function.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_in !< Input intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    ! Direction indices.
    integer :: idir_dn ! What the routine thinks is downwards.
    integer :: idir_up ! What the routine thinks is upwards.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: isplit ! Over sublayers in layer.

    ! Total split layer index.
    integer :: isplit_total

    ! Temporary variable for the intensity that has been transmitted from
    ! previous layer up to the current level.
    real, dimension(nst,nstrhalf) :: transmitted

    ! Reverse indices when necessary.
    if (flag_reverse) then
      idir_dn = iup
      idir_up = idn
    else
      idir_dn = idn
      idir_up = iup
    endif

    ! Downward radiation, transmit from top of the atmosphere down to the
    ! surface. High in the atmosphere, layer indices are low, so the layer
    ! iterator goes from low to high.

    ! During the loop, we will first progress the current intensity through
    ! the next layer and then add the intensity from the source at the
    ! bottom of that layer.

    ! The source at the first layer has to be initialized.
    transmitted = intens_in(:,:,idir_dn,0)
    res(:,:,idir_dn,0) = transmitted

    ! Initialize total sublayer counter.
    isplit_total = 0

    do ilay = 1,nlay

      do isplit = 1,tlc%nsplit(ilay)

        ! Count total split layer.
        isplit_total = isplit_total + 1

        ! Apply extinction by layer isplit_total (between levels isplit_total-1 and isplit_total).
        do istr = 1,nstrhalf
          transmitted(:,istr) = transmitted(:,istr) * tlc%t_sp(istr,ilay) ! Array operation over Stokes parameters.
        enddo

        ! Add source of level isplit_total (below the layer where we went through).
        transmitted = transmitted + intens_in(:,:,idir_dn,isplit_total)

        ! That is the output at level isplit_total.
        res(:,:,idir_dn,isplit_total) = transmitted

      enddo

    enddo

    ! Now repeat the process in reverse direction for upward direction.
    transmitted = intens_in(:,:,idir_up,nsplit_total)
    res(:,:,idir_up,nsplit_total) = transmitted

    isplit_total = nsplit_total

    do ilay = nlay,1,-1

      do isplit = 1,tlc%nsplit(ilay)

        ! Count layer (in reverse direction).
        isplit_total = isplit_total - 1

        ! Apply extinction by layer isplit_total+1 (between levels isplit_total+1 and isplit_total).
        do istr = 1,nstrhalf
          transmitted(:,istr) = transmitted(:,istr) * tlc%t_sp(istr,ilay) ! Array operation over Stokes parameters.
        enddo

        ! Add source of level isplit_total (above the layer where we went through).
        transmitted = transmitted + intens_in(:,:,idir_up,isplit_total)

        ! That is the output at level isplit_total.
        res(:,:,idir_up,isplit_total) = transmitted

      enddo

    enddo

  end subroutine transmit_sp ! }}}

  !> First-order source or response in original vertical grid.
  !!
  !! This routine handles the transmission from the sun, one scatter event and the intra-layer
  !! transmission to the first interface after the scatter event. The resulting field is saved.
  !! Also, this routine can handle exactly the opposite, the intra-layer transmission from the
  !! last interface to the last scatter event, that last scatter event and the transmission to
  !! the instrument. The latter is referred to as the response. In both cases, no layer splitting
  !! is applied, so it can only be applied for analytic (limited-scattering) calculations.
  subroutine source_response_1(nstrhalf,nlay,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nlay), intent(out) :: res !< The resulting field.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.

    ! Direction indices.
    integer :: idir_fwd
    integer :: idir_bck

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:,:), pointer :: l_px
    real, dimension(:,:), pointer :: l_mx
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: c_emi

    ! Intermediate result for thermal emission.
    ! There is no dimension on Stokes parameters, becuase
    ! thermal emission can only be one.
    real, dimension(nstrhalf) :: res_therm

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

    ! Emission, asymmetric for normal and adjoint. This differs by the sign of
    ! V-polarized light, because the adjoint flips it and the first C-parameter
    ! is taken responsible for this sign switch.
    if (flag_adjoint) then
      c_emi => tlc%c_emi_adj
    else
      c_emi => tlc%c_emi_nrm
    endif

    ! Instrument or sun oriented TLC-parameters. Those depend on whether it
    ! is the sun or the instrument where we are looking at.
    if (flag_source .eqv. flag_adjoint) then ! Adjoint source or normal response.
      ! Direction is the instrument (v).
      t_tot_x => tlc%t_tot_v
      l_px => tlc%l_pv
      l_mx => tlc%l_mv
    else
      ! Direction is the sun (0).
      t_tot_x => tlc%t_tot_0
      l_px => tlc%l_p0
      l_mx => tlc%l_m0
    endif

    ! Initialize at zero, the layer that will not be filled.
    res(:,:,idir_fwd,0) = 0.0

    do ilay = 1,nlay

      do istr = 1,nstrhalf ! Destination stream.

        ! TLC-parameters for backward scattering:
        ! T = T_tot_x at interface ilay-1
        ! L = L_px
        ! C = Cb_x

        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,ilay-1) = prefactor(istr) * t_tot_x(ilay-1) * l_px(istr,ilay) * c_bck_x(:,istr,ilay)

        ! TLC-parameters for forward scattering:
        ! T = T_tot_x at interface ilay-1 * T
        ! L = L_mx
        ! C = Cf_x

        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_fwd,ilay) = prefactor(istr) * t_tot_x(ilay-1) * tlc%t(istr,ilay) * l_mx(istr,ilay) * c_fwd_x(:,istr,ilay)

      enddo

      ! TLC-parameters for thermal emissivity: Only relevant for normal source or adjoint
      ! response. Both directions get the same contribution because of symnmetry.
      ! T = 1
      ! L = L
      ! C + C_T
      if (exist_therm) then
        res_therm = prefactor * tlc%l_p(:,ilay) * tlc%c_therm(ilay)
        res(ist_therm,:,idir_bck,ilay-1) = res(ist_therm,:,idir_bck,ilay-1) + res_therm
        res(ist_therm,:,idir_fwd,ilay) = res(ist_therm,:,idir_fwd,ilay) + res_therm
      endif

    enddo

    ! TLC-parameters for surface reflection:
    ! T = T_tot_x at bottom
    ! L = 1
    ! C = Cb_srf_x
    if (exist_bdrf) then
      do istr = 1,nstrhalf
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nlay) = prefactor(istr) * t_tot_x(nlay) * c_bck_srf_x(:,istr)
      enddo
    else
      res(:,:,idir_bck,nlay) = 0.0
    endif

    ! Emission, only relevant for normal source or adjoint response.
    ! The direction is backward.
    if (exist_emi) then

      ! TLC-parameters for surface-emission:
      ! T = 1
      ! L = 1
      ! C = Ce
      do istr = 1,nstrhalf
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nlay) = res(:,istr,idir_bck,nlay) + prefactor(istr) * c_emi(:,istr)
      enddo
    endif

  end subroutine source_response_1 ! }}}

  !> First-order source or response in split vertical grid.
  !!
  !! This routine handles the transmission from the sun, one scatter event and the intra-layer
  !! transmission to the first interface after the scatter event. The resulting field is saved.
  !! Also, this routine can handle exactly the opposite, the intra-layer transmission from the
  !! last interface to the last scatter event, that last scatter event and the transmission to
  !! the instrument. The latter is referred to as the response. In both cases, layer splitting
  !! is applied, so it should be used when performing interpolated (multi-scattering) calculations.
  subroutine source_response_1_sp(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: isplit ! Over sub-layers in one layer.

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
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: c_emi

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Intermediate result.
    ! There is no dimension on Stokes parameters, becuase
    ! thermal emission can only be one.
    real, dimension(nstrhalf) :: res_therm

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

    ! Emission, asymmetric for normal and adjoint. This differs by the sign of
    ! V-polarized light, because the adjoint flips it and the first C-parameter
    ! is taken responsible for this sign switch.
    if (flag_adjoint) then
      c_emi => tlc%c_emi_adj
    else
      c_emi => tlc%c_emi_nrm
    endif

    ! Initialize at zero, the layer that will not be filled.
    res(:,:,idir_fwd,0) = 0.0

    ! Initialize at the top of the first external layer.
    isplit_start = 0

    do ilay = 1,nlay

      ! Get other boundary of external layer.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      do istr = 1,nstrhalf ! Destination stream.

        ! TLC-parameters for backward scattering:
        ! T = T_tot_x at interface above current sublayer
        ! L = L_px
        ! C = Cb_x

        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,isplit_start) = prefactor(istr) * t_tot_x(ilay-1) * l_sp_px(istr,ilay) * c_bck_x(:,istr,ilay)

        ! TLC-parameters for forward scattering:
        ! T = T_tot_x at interface above current sublayer * T
        ! L = L_mx
        ! C = Cf_x

        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_fwd,isplit_start+1) = prefactor(istr) * t_tot_x(ilay-1) * tlc%t_sp(istr,ilay) * l_sp_mx(istr,ilay) * c_fwd_x(:,istr,ilay)

      enddo

      ! The subsequent sublayers have the same source, except that the
      ! parameter T gets another factor t_sp_x(ilay) for each next
      ! sublayer.
      do isplit = isplit_start+1,isplit_end-1 ! Excludes the first sublayer, because it has already been done.
        ! isplit is the index for backward scattering. The index for
        ! forward scattering will be isplit+1.
        ! Array operation over Stokes parameters and streams.
        res(:,:,idir_bck,isplit) = res(:,:,idir_bck,isplit-1) * t_sp_x(ilay)
        res(:,:,idir_fwd,isplit+1) = res(:,:,idir_fwd,isplit) * t_sp_x(ilay)
      enddo
      ! Note that in this analytic (non-interpolated) source or response function,
      ! contribution from different sublayers do not interfere with each other.

      ! TLC-parameters for thermal emissivity: Only relevant for normal source or adjoint
      ! response. Both directions get the same contribution because of symnmetry.
      ! This contribution is not subject to the t_sp_x loop just above.
      ! T = 1
      ! L = L
      ! C + C_T
      if (exist_therm) then
        res_therm = prefactor * tlc%l_sp_p(:,ilay) * tlc%c_therm(ilay)
        do isplit = isplit_start+1,isplit_end
          res(ist_therm,:,idir_bck,isplit-1) = res(ist_therm,:,idir_bck,isplit-1) + res_therm
          res(ist_therm,:,idir_fwd,isplit) = res(ist_therm,:,idir_fwd,isplit) + res_therm
        enddo
      endif

      ! Progress to the next external layer, dodging the sum-expression.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = T_tot_x at bottom
    ! L = 1
    ! C = Cb_srf_x
    if (exist_bdrf) then
      do istr = 1,nstrhalf ! Destination stream.
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nsplit_total) = prefactor(istr) * t_tot_x(nlay) * c_bck_srf_x(:,istr)
      enddo
    else
      res(:,:,idir_bck,nsplit_total) = 0.0
    endif

    ! Emission, only relevant for normal source or adjoint response.
    ! The direction is backward.
    if (exist_emi) then

      ! TLC-parameters for surface-emission:
      ! T = 1
      ! L = 1
      ! C = Ce
      do istr = 1,nstrhalf
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nsplit_total) = res(:,istr,idir_bck,nsplit_total) + prefactor(istr) * c_emi(:,istr)
      enddo

    endif

  end subroutine source_response_1_sp ! }}}

  !> First-order source or response in split vertical grid, coupling to interpolated field.
  !!
  !! This routine handles the transmission from the sun and one scatter event and end somewhere
  !! in the middle of the layer. The resulting field is saved. Also, and more commonly used, this
  !! routine can handle exactly the opposite, the last scatter event and the transmission from
  !! somewhere in the middle of a layer to the instrument. The latter is referred to as the
  !! response. The coupling to the `middle' of the layer is possible by interpolation. When
  !! performing a normal dot product (sum of the product of the overlapping elements), the
  !! vertical integral is automatically correctly performed.
  subroutine source_response_1_int(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_bdrf,exist_emi,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_therm
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: isplit ! Over sublayers in layer.

    ! Direction indices.
    integer :: idir_fwd ! For forward scattering.
    integer :: idir_bck ! For backward scattering.

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:), pointer :: l_int_sm_x
    real, dimension(:), pointer :: l_int_ac_x
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x
    real, dimension(:,:), pointer :: c_bck_srf_x
    real, dimension(:,:), pointer :: c_emi

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Each sublayer adds to the result an array. The next sublayer adds the
    ! same array with an extra transmission term. The misery is that these
    ! arrays overlap half in the result, so the array of the previous
    ! sublayer cannot be acquired from the result array. Therefore, it
    ! should be saved in a temporary array.
    real, dimension(nst,nstrhalf,0:1,0:1) :: res_sublayer

    ! Intermediate result for thermal emission.
    ! There is no dimension on Stokes parameters, becuase
    ! thermal emission can only be one.
    real, dimension(nstrhalf) :: res_therm

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

    ! Emission, asymmetric for normal and adjoint. This differs by the sign of
    ! V-polarized light, because the adjoint flips it and the first C-parameter
    ! is taken responsible for this sign switch.
    if (flag_adjoint) then
      c_emi => tlc%c_emi_adj
    else
      c_emi => tlc%c_emi_nrm
    endif

    ! Initialize at zero, so that different parts can be added.
    res = 0.0

    ! Initialize at the top of the first external layer.
    isplit_start = 0

    do ilay = 1,nlay

      ! Get other boundary of external layer.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      do istr = 1,nstrhalf ! Destination stream.

        ! TLC-parameters for backward scattering:
        ! T = T_tot_x at interface above current sublayer
        ! L = L_int_x
        ! C = Cb_x

        ! Same-layer parameters always for high layer (low index), because radiation
        ! communicates to sun or instrument via the top of the layer.

        ! Array operation over destination Stokes parameters.
        res_sublayer(:,istr,idir_bck,0) = prefactor(istr) * t_tot_x(ilay-1) * l_int_sm_x(ilay) * c_bck_x(:,istr,ilay)
        res_sublayer(:,istr,idir_bck,1) = prefactor(istr) * t_tot_x(ilay-1) * l_int_ac_x(ilay) * c_bck_x(:,istr,ilay)

        ! TLC-parameters for forward scattering:
        ! The same, except that Cb_x is now Cf_x. Also, the remark on same/across
        ! L-terms is the same as for backward scattering.

        ! Array operation over destination Stokes parameters.
        res_sublayer(:,istr,idir_fwd,0) = prefactor(istr) * t_tot_x(ilay-1) * l_int_sm_x(ilay) * c_fwd_x(:,istr,ilay)
        res_sublayer(:,istr,idir_fwd,1) = prefactor(istr) * t_tot_x(ilay-1) * l_int_ac_x(ilay) * c_fwd_x(:,istr,ilay)

      enddo

      ! Write first sublayer, after which there is a loop adding transmission terms
      ! for the subsequent sublayers.
      res(:,:,:,isplit_start:isplit_start+1) = res(:,:,:,isplit_start:isplit_start+1) + res_sublayer

      ! The subsequent sublayers have the same source, except that the
      ! parameter T gets another factor t_sp_x(ilay) for each next
      ! sublayer.
      do isplit = isplit_start+1,isplit_end-1
        res_sublayer = res_sublayer * t_sp_x(ilay)
        ! isplit is the index for backward scattering. The index for forward
        ! scattering will be isplit+1.
        res(:,:,:,isplit:isplit+1) = res(:,:,:,isplit:isplit+1) + res_sublayer
      enddo

      ! TLC-parameters for thermal emission:
      ! T = 1
      ! L = L_int
      ! C = C_T
      if (exist_therm) then
        res_therm = prefactor * tlc%l_int(ilay) * tlc%c_therm(ilay)
        do isplit = isplit_start+1,isplit_end
          ! For both sides and for both directions, the contribution is the same.
          res(ist_therm,:,idir_bck,isplit-1) = res(ist_therm,:,idir_bck,isplit-1) + res_therm
          res(ist_therm,:,idir_fwd,isplit-1) = res(ist_therm,:,idir_fwd,isplit-1) + res_therm
          res(ist_therm,:,idir_bck,isplit) = res(ist_therm,:,idir_bck,isplit) + res_therm
          res(ist_therm,:,idir_fwd,isplit) = res(ist_therm,:,idir_fwd,isplit) + res_therm
        enddo
      endif

      ! Progress to the next external layer, dodging the sum-expression.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = T_tot_x at bottom
    ! L = 1
    ! C = Cb_srf_x
    if (exist_bdrf) then
      do istr = 1,nstrhalf
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nsplit_total) = res(:,istr,idir_bck,nsplit_total) + prefactor(istr) * t_tot_x(nlay) * c_bck_srf_x(:,istr)
      enddo
    endif
    ! It is an addition, so nothing has to be done if exist_bdrf is false.

    ! Emission, only relevant for normal source or adjoint response.
    ! The direction is backward.
    if (exist_emi) then

      ! TLC-parameters for surface-emission:
      ! T = 1
      ! L = 1
      ! C = Ce
      do istr = 1,nstrhalf
        ! Array operation over destination Stokes parameters.
        res(:,istr,idir_bck,nsplit_total) = res(:,istr,idir_bck,nsplit_total) + prefactor(istr) * c_emi(:,istr)
      enddo

    endif

  end subroutine source_response_1_int ! }}}

  !> Second-order source or response in split vertical grid.
  !!
  !! This routine handles the transmission from the sun, two scatter events in one layer and
  !! the intra-layer transmission to the first interface after the scatter events. The resulting
  !! field is saved. Also, this routine can handle exactly the opposite, the intra-layer
  !! transmission from the last interface to the last scatter events, with the last two scatter events
  !! in the same layer, and the transmission to the instrument. The latter is referred to as the
  !! response. In both cases, layer splitting is applied, so it should be used when performing
  !! interpolated (multi-scattering) calculations. Because of the application of this field,
  !! the result is added to the pre-existing field.
  subroutine source_response_2_sp_addition(nstrhalf,nlay,nsplit_total,flag_source,flag_adjoint,exist_therm,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup, ist_i, ist_u, ist_v, ist_therm, fullpol_nst
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
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(inout) :: res !< The resulting field added to the original values.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over source streams.
    integer :: isplit ! Over sublayers in one layer.

    ! Pre-aggregated TLC-parameters that appear repeatedly. For the matrix, the
    ! dimension orders are switched, so that they can be interpreted with merged
    ! dimensions of Stokes parameters and streams with reshape pointers.
    real, dimension(nst,nstrhalf), target :: arr_aggregated_i ! Intermediate direction.
    real, dimension(nst,nstrhalf), target :: arr_aggregated_d ! Destination direction.
    real, dimension(nst,nstrhalf,nst,nstrhalf), target :: mat_aggregated ! Both directions.

    ! Result of the topmost sublayer of a layer.
    real, dimension(nst,nstrhalf,0:1), target :: res_sublayer

    ! Reshape pointers, where Stokes parameters and streams dimensions will be combined.
    real, dimension(:), pointer :: arr_i_ptr
    real, dimension(:), pointer :: arr_d_ptr
    real, dimension(:,:), pointer :: mat_ptr
    real, dimension(:), pointer :: res_ptr

    ! Reshape pointers for thermal contribution. Different pointers, because their dimensions
    ! are a bit different, because one Stokes parameter dimension vanishes into ist_therm.
    real, dimension(:), pointer :: arr_i_therm_ptr
    real, dimension(:), pointer :: arr_d_therm_ptr
    real, dimension(:,:), pointer :: mat_therm_ptr
    real, dimension(:), pointer :: res_therm_ptr

    ! Direction indices.
    integer :: idir_fwd ! For forward scattering.
    integer :: idir_bck ! For backward scattering.

    ! Pointers to TLC-parameters with the sun or the instrument.
    ! Here, x means either the sun or the instrument.
    real, dimension(:), pointer :: prefactor
    real, dimension(:,:), pointer :: effmu_mx
    real, dimension(:), pointer :: t_tot_x
    real, dimension(:), pointer :: t_sp_x
    real, dimension(:,:), pointer :: l_sp_mx
    real, dimension(:,:), pointer :: l_sp_px
    real, dimension(:,:), pointer :: l_1_sp_px
    real, dimension(:,:,:), pointer :: c_fwd_x
    real, dimension(:,:,:), pointer :: c_bck_x

    ! Temporary TLC-parameters, because C-parameters are asymmetric due to V-polarization.
    ! These are only needed if V-polarization is involved.
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_fwd_id
    real, dimension(nst,nst,nstrhalf,nstrhalf) :: c_bck_id

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Intermediate result for thermal emission.
    ! Now, we have a dimension on Stokes parameters, because there
    ! is a scatter event from the thermal emission to the destination
    ! direction.
    real, dimension(nst,nstrhalf), target :: res_therm

    ! Temporary arrays and matrices for the thermal part, where one Stokes
    ! parameter becomes fixed.
    real, dimension(nstrhalf), target :: arr_aggregated_i_therm ! No Stokes parameters here.
    real, dimension(nstrhalf,nst,nstrhalf), target :: mat_aggregated_therm ! Intermediate stream has no Stokes parameter.

    ! Note that this routine is the only routine that is designed to add the
    ! result to a pre-existing result. Therefor there is certainly no zero-
    ! initialization.

    ! Set fixed reshape pointers.
    arr_i_ptr(1:nst*nstrhalf) => arr_aggregated_i
    arr_d_ptr(1:nst*nstrhalf) => arr_aggregated_d
    mat_ptr(1:nst*nstrhalf,1:nst*nstrhalf) => mat_aggregated
    if (exist_therm) then
      arr_i_therm_ptr => arr_aggregated_i_therm ! It stays 1D, no reshape, but that is not trivial.
      arr_d_therm_ptr(1:nst*nstrhalf) => arr_aggregated_d ! It stays the same as fot the non-thermal array, but that is not trivial.
      mat_therm_ptr(1:nstrhalf,1:nst*nstrhalf) => mat_aggregated_therm ! Reshape from 3D to 2D.
      res_therm_ptr(1:nst*nstrhalf) => res_therm
    endif

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
      ! Special direction is the instrument (v).
      effmu_mx => tlc%effmu_mv
      t_tot_x => tlc%t_tot_v
      t_sp_x => tlc%t_sp_v
      l_sp_px => tlc%l_sp_pv
      l_sp_mx => tlc%l_sp_mv
      l_1_sp_px => tlc%l_1_sp_pv
    else
      ! Special direction is the sun (0).
      effmu_mx => tlc%effmu_m0
      t_tot_x => tlc%t_tot_0
      t_sp_x => tlc%t_sp_0
      l_sp_px => tlc%l_sp_p0
      l_sp_mx => tlc%l_sp_m0
      l_1_sp_px => tlc%l_1_sp_p0
    endif

    ! The prefactor is prefactor (d) times prefactor_intermediate (i)
    ! where 'prefactor' is either prefactor_src or prefactor_resp,
    ! depending on if we are the source or the response function.

    ! Initialize at the top of the first external layer.
    isplit_start = 0

    do ilay = 1,nlay

      ! Get other boundary of external layer.
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

      ! Net backward scattering.
      ! Set reshape pointer of the sublayer result to the 0 index, because
      ! backward scattering is always related to the top of the sublayer.
      res_ptr(1:nst*nstrhalf) => res_sublayer(:,:,0)

      ! TLC-parameters for (forward - back): {{{

      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_mx(i) * (L_pp(i,d) - L_px(d))
      ! C = Cf_x(i) * Cb(i,d)

      ! Scalar      : T_tot_x
      ! Array in i  : prefactor_intermediate, effmu_mx, Cf_x
      ! Array in d  : prefactor, -L_px
      ! Matrix(i,d) : L_pp, Cb

      ! The array over destination (d) is Stokes-parameter-independent, but will be
      ! represented with Stokes parameters (nst copies), so it has the same
      ! structure as the array over intermediate (i).

      ! Start with -L_px in destination direction.
      do istr = 1,nstrhalf ! General stream.
        arr_aggregated_i(:,istr) = tlc%prefactor_intermediate(istr) * c_fwd_x(:,istr,ilay) * effmu_mx(istr,ilay)
        arr_aggregated_d(:,istr) = -t_tot_x(ilay-1) * prefactor(istr) * l_sp_px(istr,ilay)
        ! Use the single stream index to re-arrange the dimensions of Cb
        ! istr is here the intermediate (i) stream.
        mat_aggregated(:,istr,:,:) = c_bck_id(:,:,istr,:)
      enddo

      ! result = array_d * matmul(array_i,matrix_id)
      res_ptr = arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! Reconstruct the array over d without her -L_px
      ! and add the L_pp to the aggregated matrix.

      do istr = 1,nstrhalf ! Destination stream.
        arr_aggregated_d(:,istr) = t_tot_x(ilay-1) * prefactor(istr)
        ! Adding this L-parameter needs another iterator.
        do istr_src = 1,nstrhalf ! Intermediate stream.
          mat_aggregated(:,istr_src,:,istr) = mat_aggregated(:,istr_src,:,istr) * tlc%l_sp_pp(istr_src,istr,ilay)
        enddo
      enddo

      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! }}}

      ! TLC-parameters for (back - forward) for i unequal to d: {{{

      ! T = T_tot_x at interface above current sublayer
      ! L = effmu_mp(i,d) * (L_px(i) - L_px(d))
      ! C = Cb_x(i), Cf(i,d)

      ! Scalar      : T_tot_x
      ! Array in i  : prefactor_intermediate, L_px, Cb_x
      ! Array in d  : prefactor, -L_px
      ! Matrix(i,d) : effmu_mp, Cf

      ! First do the -L_px part in the destination stream.
      do istr = 1,nstrhalf ! General stream.
        arr_aggregated_i(:,istr) = tlc%prefactor_intermediate(istr) * c_bck_x(:,istr,ilay)
        arr_aggregated_d(:,istr) = -t_tot_x(ilay-1) * prefactor(istr) * l_sp_px(istr,ilay)
        ! Use another iterator for the aggregated matrix. The original iterator
        ! is the d (destination) iterator and the new one because the i (intermediate)
        ! iterator, which is more like the source.
        do istr_src = 1,nstrhalf ! Intermediate stream
          ! Throw away the diagonal terms, because of zero-division challenges.
          if (istr_src .eq. istr) then
            mat_aggregated(:,istr_src,:,istr) = 0.0 ! We restrict to unequal streams, so equal streams are set to zero.
          else
            mat_aggregated(:,istr_src,:,istr) = tlc%effmu_mp(istr_src,istr) * c_fwd_id(:,:,istr_src,istr)
          endif

        enddo
        
      enddo

      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! Reconstruct the array over d without her -L_px
      ! and add the L_px to the array over i.

      do istr = 1,nstrhalf ! General stream.
        arr_aggregated_d(:,istr) = t_tot_x(ilay-1) * prefactor(istr)
        arr_aggregated_i(:,istr) = arr_aggregated_i(:,istr) * l_sp_px(istr,ilay)
      enddo

      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! }}}

      ! TLC-parameters for (back - forward) for i equal to d: {{{

      ! T = T_tot_x at interface above current sublayer
      ! L = L_1_px(d)
      ! C = Cb_x(i), Cf(i,d)

      ! There is only one L-term, so no pre-aggregation of arrays is
      ! required.

      ! In this case d = i, two streams become one, but two Stokes parameters
      ! remain two. Therefore, we refrain from using the reshape pointers here.
      do istr = 1,nstrhalf ! Both stream i and stream d.
        ! The matmul is over Stokes parameters.
        ! array_d = scalar * matmul(array_i , matrix_id), with addition
        res_sublayer(:,istr,0) = res_sublayer(:,istr,0) + prefactor(istr) * tlc%prefactor_intermediate(istr) * t_tot_x(ilay-1) * l_1_sp_px(istr,ilay) * matmul(c_bck_x(:,istr,ilay) , c_fwd_id(:,:,istr,istr))
      enddo

      ! }}}

      ! Net forward scattering.
      ! Set reshape pointer of the sublayer result to the 1 index.
      res_ptr(1:nst*nstrhalf) => res_sublayer(:,:,1)

      ! TLC-parameters for (back - back): {{{

      ! T = T_tot_x at interface above current sublayer * T(d)
      ! L = effmu_pp(i,d) * (L_mx(d) - L_px(i))
      ! C = Cb_x(i), Cb(i,d)

      ! Scalar      : T_tot_x
      ! Array in i  : prefactor_intermediate, -L_px, Cb_x
      ! Array in d  : prefactor, T, L_mx
      ! Matrix(i,d) : effmu_pp, Cb

      ! First do the L_mx in the destination direction.
      do istr = 1,nstrhalf ! General stream.
        arr_aggregated_i(:,istr) = tlc%prefactor_intermediate(istr) * c_bck_x(:,istr,ilay)
        arr_aggregated_d(:,istr) = t_tot_x(ilay-1) * prefactor(istr) * tlc%t_sp(istr,ilay) * l_sp_mx(istr,ilay)
        ! A new iterator for the matrix. The `source' iterator istr_src will be
        ! direction `i', the direction that is at the source side of the TLC-parameters
        ! in the matrix.
        do istr_src = 1,nstrhalf ! Intermediate stream.
          mat_aggregated(:,istr_src,:,istr) = tlc%effmu_pp(istr_src,istr) * c_bck_id(:,:,istr_src,istr)
        enddo
      enddo

      ! result = array_d * matmul(array_i,matrix_id).
      res_ptr = arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! Reconstruct the array over d without her L_mx
      ! and add the -L_px to the aggregated array over i.

      do istr = 1,nstrhalf ! General stream.
        arr_aggregated_d(:,istr) = t_tot_x(ilay-1) * prefactor(istr) * tlc%t_sp(istr,ilay)
        arr_aggregated_i(:,istr) = -arr_aggregated_i(:,istr) * l_sp_px(istr,ilay)
      enddo
        
      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! }}}

      ! TLC-parameters for (forward - forward): {{{

      ! T = T_tot_x at interface above current sublayer * T(d)
      ! L = effmu_mx(i) (L_mp(d,i) - L_mx(d))
      ! C = Cf_x(i), Cf(i,d)

      ! Scalar      : T_tot_x
      ! Array in i  : prefactor_intermediate, effmu_mx, Cf_x
      ! Array in d  : prefactor, T, -L_mx
      ! Matrix(i,d) : L_mp, Cf

      ! L_mp has to be transposed, because it must be on (d,i), while our matrix
      ! is at (i,d).

      ! First do the -L_mx in the destination array.
      do istr = 1,nstrhalf
        arr_aggregated_i(:,istr) = tlc%prefactor_intermediate(istr) * effmu_mx(istr,ilay) * c_fwd_x(:,istr,ilay)
        arr_aggregated_d(:,istr) = -t_tot_x(ilay-1) * prefactor(istr) * tlc%t_sp(istr,ilay) * l_sp_mx(istr,ilay)
        ! Use the single stream index to re-arrange the dimensions of Cf.
        mat_aggregated(:,istr,:,:) = c_fwd_id(:,:,istr,:)
      enddo

      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! Reconstruct the array over d without her -L_mx
      ! and add the L_mp to the aggregated matrix.

      do istr = 1,nstrhalf
        arr_aggregated_d(:,istr) = t_tot_x(ilay-1) * prefactor(istr) * tlc%t_sp(istr,ilay)
        ! Use a second iterator for direction `i'. That is the first dimension in the
        ! matrix, but in this case, it is the second dimension of the L-parameter.
        do istr_src = 1,nstrhalf
          mat_aggregated(:,istr_src,:,istr) = mat_aggregated(:,istr_src,:,istr) * tlc%l_sp_mp(istr,istr_src,ilay)
        enddo
      enddo

      ! result = array_d * matmul(array_i,matrix_id), but now with addition.
      res_ptr = res_ptr + arr_d_ptr * matmul(arr_i_ptr , mat_ptr)

      ! }}}

      ! All sublayers have the same source, except that the
      ! parameter T gets another factor t_sp_x(ilay) for each next
      ! sublayer, starting at exactly the res_sublayer for the first layer.
      do isplit = isplit_start,isplit_end-1
        ! isplit is the index for backward scattering. The index for
        ! forward scattering will be isplit+1.
        res(:,:,idir_bck,isplit) = res(:,:,idir_bck,isplit) + res_sublayer(:,:,0)
        res(:,:,idir_fwd,isplit+1) = res(:,:,idir_fwd,isplit+1) + res_sublayer(:,:,1)
        res_sublayer = res_sublayer * t_sp_x(ilay)
      enddo

      ! Thermal emission. This is not subject ot the loop above, where t_sp_x is
      ! added for each sublayer. As thermal emission is equal downward and upward,
      ! all combinations with scatter events apply equally for upward and downward
      ! destinations.
      if (exist_therm) then

        ! In this section, the intermediate Stokes parameter is always ist_therm.
        ! This means that a dimension changes in the aggregated array over i and
        ! over the matrix. The array over d does not change, but we will use the
        ! thermal pointer to avoid confusion.

        ! TLC-parameters for (thermal - back): ! {{{
        ! T = 1
        ! L = effmu_p(i) * (L_p(d) - L_pp(i,d))
        ! C = C_T * Cb(i,d)

        ! Scalar      : C_T
        ! Array in i  : prefactor_intermediate, effmu_p
        ! Array in d  : prefactor, L_p
        ! Matrix(i,d) : L_pp, Cb

        ! First, we take the L_p contribution.

        do istr = 1,nstrhalf ! General stream.
          arr_aggregated_d(:,istr) = tlc%c_therm(ilay) * prefactor(istr) * tlc%l_sp_p(istr,ilay)
          arr_aggregated_i_therm(istr) = tlc%prefactor_intermediate(istr) * tlc%effmu_p(istr)
          mat_aggregated_therm(istr,:,:) = c_bck_id(ist_therm,:,istr,:) ! Array operation over both destination Stokes parameter and destination stream. The istr is the source (intermediate) stream here.
        enddo

        ! result = array_d * matmul(array_i,matrix_id).
        res_therm_ptr = arr_d_therm_ptr * matmul(arr_i_therm_ptr , mat_therm_ptr)

        ! Reconstruct the array over d without her L_p
        ! and add the -L_pp to the aggregated matrix.

        do istr = 1,nstrhalf ! General stream.
          arr_aggregated_d(:,istr) = tlc%c_therm(ilay) * prefactor(istr)
          ! Adding this L-parameter needs another iterator.
          do istr_src = 1,nstrhalf ! Intermediate stream.
            mat_aggregated_therm(istr_src,:,istr) = -mat_aggregated_therm(istr_src,:,istr) * tlc%l_sp_pp(istr_src,istr,ilay)
          enddo
        enddo

        ! result = array_d * matmul(array_i,matrix_id), but now with addition.
        res_therm_ptr = res_therm_ptr + arr_d_therm_ptr * matmul(arr_i_therm_ptr , mat_therm_ptr)

        ! }}}

        ! TLC-parameters for (thermal - forward) for i unequal to d: {{{
        ! T = 1
        ! L = effmu_mp(i,d) * (L_p(i) - L_p(d))
        ! C = C_T * Cf(i,d)

        ! Scalar      : C_T
        ! Array in i  : prefactor_intermediate, L_p
        ! Array in d  : prefactor, L_p
        ! Matrix(i,d) : effmu_mp, Cf

        ! Let us first take the contribution with L_p in the intermediate array.
        ! The, we can keep our old array in d: C_T * prefactor.

        ! Construct the aggregated i array and the matrix.
        do istr = 1,nstrhalf ! General stream.
          arr_aggregated_i_therm(istr) = tlc%prefactor_intermediate(istr) * tlc%l_sp_p(istr,ilay)
          ! We need another iterator, because the matrix also contains a Stokesless contribution.
          do istr_src = 1,nstrhalf ! Intermediate stream.
            if (istr_src .eq. istr) then
              mat_aggregated_therm(istr_src,:,istr) = 0.0 ! We restrict to unequal streams, so equal streams are set to zero.
            else
              mat_aggregated_therm(istr_src,:,istr) = tlc%effmu_mp(istr_src,istr) * c_fwd_id(ist_therm,:,istr_src,istr)
            endif
          enddo
        enddo

        ! result = array_d * matmul(array_i,matrix_id), but now with addition.
        res_therm_ptr = res_therm_ptr + arr_d_therm_ptr * matmul(arr_i_therm_ptr , mat_therm_ptr)

        ! Reconstruct the array over i without her L_p
        ! and add the -L_p to the aggregated array over d.

        ! We can just use tlc%prefactor_intermediate as array over i.

        do istr = 1,nstrhalf
          arr_aggregated_d(:,istr) = -arr_aggregated_d(:,istr) * tlc%l_sp_p(istr,ilay)
        enddo

        ! result = array_d * matmul(array_i,matrix_id), but now with addition.
        res_therm_ptr = res_therm_ptr + arr_d_therm_ptr * matmul(tlc%prefactor_intermediate , mat_therm_ptr)

        ! }}}

        ! TLC-parameters for (thermal - forward) for i equal to d: {{{
        ! T = 1
        ! L = L_1(i)
        ! C = C_T * Cf(i,d)

        ! There is only one L-term, so no pre-aggregation of arrays is
        ! required.

        ! In this case d = i, two streams become one, but two Stokes parameters
        ! remain two. Therefore, we refrain from using the reshape pointers here.
        do istr = 1,nstrhalf ! Both stream i and stream d.
          res_therm(:,istr) = res_therm(:,istr) + prefactor(istr) * tlc%prefactor_intermediate(istr) * tlc%l_1_sp_p(istr,ilay) * tlc%c_therm(ilay) * c_fwd_id(ist_therm,:,istr,istr)
        enddo

        ! }}}

        ! Apply thermal emission to each sublayer.
        do isplit = isplit_start+1,isplit_end
          res(:,:,idir_bck,isplit-1) = res(:,:,idir_bck,isplit-1) + res_therm
          res(:,:,idir_fwd,isplit) = res(:,:,idir_fwd,isplit) + res_therm
        enddo

      endif

      ! Progress to the next external layer, dodging the sum-expression.
      isplit_start = isplit_end

    enddo

  end subroutine source_response_2_sp_addition ! }}}

  !> First order scattering in one sublayer in split vertical grid.
  !!
  !! The input radiation is transmitted into the sublayer it is facing, scattered in that
  !! sublayer and then transmitted to the next interface in the split grid. Because the
  !! split vertical grid is used, this routine can be used for multi-scattering.
  subroutine scatter_1_sp(nstrhalf,nlay,nsplit_total,exist_bdrf,intens_in,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in), target :: intens_in !< Input intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out), target :: res !< The resulting field.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over destination streams.
    integer :: istr_src ! Over source streams.
    integer :: isplit ! Over sublayer in layer.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Aggregated matrix that can be re-used for all split layers.
    real, dimension(nst,nstrhalf,nst,nstrhalf), target :: mat_aggregated

    ! Reshape pointers, with dimensions nst*nstrhalf.
    real, dimension(:,:), pointer :: mat_ptr
    real, dimension(:), pointer :: source_ptr
    real, dimension(:), pointer :: destination_ptr

    ! Set fixed reshape pointer.
    mat_ptr(1:nst*nstrhalf,1:nst*nstrhalf) => mat_aggregated

    ! Direction s is the source and direction d is the destination.

    ! The prefactor is prefactor_scat(s,d).

    ! Initialize at zero, the layer that will not be filled.
    res(:,:,idn,0) = 0.0

    ! Initialize at the top of the first external layer.
    isplit_start = 0

    do ilay = 1,nlay

      ! Get other boundary of external layer.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! TLC-parameters for backward scattering: {{{
      ! T = 1
      ! L = L_pp(s,d)
      ! C = Cb(s,d)

      ! Array in s  : intens_in
      ! Matrix(s,d) : prefactor_scat, L_pp, Cb

      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! Aggregate matrix that is the same for all sublayers
          ! The array operations are over all the Stokes parameters, with the
          ! source always left of the destination.
          mat_aggregated(:,istr_src,:,istr) = tlc%prefactor_scat(istr_src,istr) * tlc%l_sp_pp(istr_src,istr,ilay) * tlc%c_bck(:,:,istr_src,istr,ilay)

        enddo
      enddo

      ! Loop over sublayers so that the intensity field can be added to
      ! the TLC-aggregate.
      do isplit = isplit_start+1,isplit_end

        ! result = matmul(array_s,matrix_sd)
        ! This will be applied using reshape pointers.

        ! Set source and destination reshape pointers for this split layer (dn-up).
        source_ptr(1:nst*nstrhalf) => intens_in(:,:,idn,isplit-1)
        destination_ptr(1:nst*nstrhalf) => res(:,:,iup,isplit-1)
        destination_ptr = matmul(source_ptr , mat_ptr)

        ! Set source and destination reshape pointers for this split layer (up-dn).
        source_ptr(1:nst*nstrhalf) => intens_in(:,:,iup,isplit)
        destination_ptr(1:nst*nstrhalf) => res(:,:,idn,isplit)
        destination_ptr = matmul(source_ptr , mat_ptr)

      enddo

      ! }}}

      ! TLC-parameters for forward scattering: {{{
      ! T = T(s)
      ! L = L_mp(s,d)
      ! C = Cf(s,d)

      ! It is chosen to take T(s) * L_mp(s,d), which is the same as
      ! T(d) * L_mp(d,s), but the prefactor is already in dimension
      ! order (s,d), so we choose the same dimension array in L, so that
      ! no transpose is needed.

      ! Array in s  : intens_in, T
      ! Matrix(s,d) : prefactor_scat, L_mp, Cf

      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! Aggregate matrix that is the same for all sublayers.
          mat_aggregated(:,istr_src,:,istr) = tlc%prefactor_scat(istr_src,istr) * tlc%t_sp(istr_src,ilay) * tlc%l_sp_mp(istr_src,istr,ilay) * tlc%c_fwd(:,:,istr_src,istr,ilay)
          
        enddo
      enddo

      ! Loop over sublayers so that the intensity field can be added to
      ! the TLC-aggregate.
      do isplit = isplit_start+1,isplit_end

        ! result = matmul(array_s,matrix_sd).
        ! This will be applied using reshape pointers.

        ! Set source and destination reshape pointers for this split layer (up-up).
        source_ptr(1:nst*nstrhalf) => intens_in(:,:,iup,isplit)
        destination_ptr(1:nst*nstrhalf) => res(:,:,iup,isplit-1)
        destination_ptr = destination_ptr + matmul(source_ptr , mat_ptr)

        ! Set source and destination reshape pointers for this split layer (dn-dn).
        source_ptr(1:nst*nstrhalf) => intens_in(:,:,idn,isplit-1)
        destination_ptr(1:nst*nstrhalf) => res(:,:,idn,isplit)
        destination_ptr = destination_ptr + matmul(source_ptr , mat_ptr)

      enddo

      ! }}}

      ! Progress to the next external layer, dodging the sum-expression.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = 1
    ! L = 1
    ! C = Cb_srf

    ! Array in s  : intens_in
    ! Matrix(s,d) : prefactor_scat , Cb_srf

    ! Use the same pointers to merge Stokes parameters and streams.
    source_ptr(1:nst*nstrhalf) => intens_in(:,:,idn,nsplit_total)
    destination_ptr(1:nst*nstrhalf) => res(:,:,iup,nsplit_total)

    if (exist_bdrf) then
      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! Set aggregated matrix equal to the product of the TLC-parameters for
          ! surface reflection.
          mat_aggregated(:,istr_src,:,istr) = tlc%prefactor_scat(istr_src,istr) * tlc%c_bck_srf(:,:,istr_src,istr)

        enddo
      enddo

      ! result = matmul(array_s,matrix_sd), but now without sublayer loop
      destination_ptr = matmul(source_ptr , mat_ptr)
    else
      destination_ptr = 0.0
    endif

  end subroutine scatter_1_sp ! }}}

  !> First order scattering in one sublayer in split vertical grid, coupling to an interpolated field.
  !!
  !! The input radiation is scattered at its place, which is in the middle of a sublayer, then
  !! it is transmitted to the first interface. The input field is sampled such that the
  !! field in the middle is equal to the interpolated field.
  subroutine scatter_1_int(nstrhalf,nlay,nsplit_total,exist_bdrf,intens_in,tlc,res) ! {{{

    use lintran_constants_module, only: nst_fake, idn, iup
    use lintran_types_module, only: lintran_tlc_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: nsplit_total !< Total number of sublayers.
    logical, intent(in) :: exist_bdrf !< Flag for if surface reflection exist on this Fourier number.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: intens_in !< Input intensity field.
    type(lintran_tlc_class), intent(in) :: tlc !< Internal radiative-transfer building-block structure.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(out) :: res !< The resulting field.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over source streams.

    ! Split interface indices surrounding the perturbed layer.
    integer :: isplit_start
    integer :: isplit_end

    ! Direction s is the source and direction d is the destination.

    ! The prefactor is prefactor_scat(s,d)

    ! Initialize result to zero so different contributions can be acquired by addition.
    res = 0.0

    ! Initialize at the top of the first external layer.
    isplit_start = 0

    do ilay = 1,nlay

      ! Get other boundary of external layer.
      isplit_end = isplit_start + tlc%nsplit(ilay)

      ! TLC-parameters for atmospheric scattering:
      ! T = 1
      ! L = L_int_p(d)
      ! C = Cb(s,d) or Cf(s,d)

      ! Array in s  : intens_in
      ! Array in d  : L_int_p
      ! Matrix(s,d) : prefactor_scat, C

      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.

          ! Everything that ends up as upward radiation.
          ! This is a 2D array assignment on destination Stokes parameters and sublayers.
          ! The matmuls have to sum over both arguments' first dimensions, namely the
          ! source Stokes parameters, so the first one is transposed.
          res(:,istr,iup,isplit_start:isplit_end-1) = res(:,istr,iup,isplit_start:isplit_end-1) + tlc%prefactor_scat(istr_src,istr) * (tlc%l_int_sm_p(istr,ilay) * matmul(transpose(tlc%c_bck(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,idn,isplit_start:isplit_end-1)) + tlc%l_int_ac_p(istr,ilay) * matmul(transpose(tlc%c_bck(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,idn,isplit_start+1:isplit_end)) + tlc%l_int_sm_p(istr,ilay) * matmul(transpose(tlc%c_fwd(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,iup,isplit_start:isplit_end-1)) + tlc%l_int_ac_p(istr,ilay) * matmul(transpose(tlc%c_fwd(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,iup,isplit_start+1:isplit_end)))

          ! Everything that ends up as downard radiation.
          ! Same comment about matrix multiplication and array assignment.
          res(:,istr,idn,isplit_start+1:isplit_end) = res(:,istr,idn,isplit_start+1:isplit_end) + tlc%prefactor_scat(istr_src,istr) * (tlc%l_int_sm_p(istr,ilay) * matmul(transpose(tlc%c_bck(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,iup,isplit_start+1:isplit_end)) + tlc%l_int_ac_p(istr,ilay) * matmul(transpose(tlc%c_bck(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,iup,isplit_start:isplit_end-1)) + tlc%l_int_sm_p(istr,ilay) * matmul(transpose(tlc%c_fwd(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,idn,isplit_start+1:isplit_end)) + tlc%l_int_ac_p(istr,ilay) * matmul(transpose(tlc%c_fwd(:,:,istr_src,istr,ilay)) , intens_in(:,istr_src,idn,isplit_start:isplit_end-1)))

        enddo
      enddo

      ! Progress to the next external layer, dodging the sum-expression.
      isplit_start = isplit_end

    enddo

    ! TLC-parameters for surface reflection:
    ! T = 1
    ! L = 1
    ! C = Cb_srf

    ! Array in s  : intens_in
    ! Matrix(s,d) : prefactor_scat , Cb_srf

    if (exist_bdrf) then
      ! result = matmul(array_s,matrix_sd).
      do istr = 1,nstrhalf ! Destination stream.
        do istr_src = 1,nstrhalf ! Source stream.
          res(:,istr,iup,nsplit_total) = res(:,istr,iup,nsplit_total) + tlc%prefactor_scat(istr_src,istr) * matmul(intens_in(:,istr_src,idn,nsplit_total) , tlc%c_bck_srf(:,:,istr_src,istr))
        enddo
      enddo
    endif
    ! It is an addition. The result was already initialized to zero, because the
    ! calculation inside the atmosphere is a bit more complicated. Still, nothing needs
    ! to be done when exist_bdrf is false.

  end subroutine scatter_1_int ! }}}

  !> Calculates the inner product.
  !!
  !! This routine calculates the observable by multiplying an intensity field with a response
  !! funciton. The inner product is simply the sum of the products of the elements with the
  !! same indices. This routine can be used for both the original vertical grid and for the
  !! split vertical grid, as long as the two input fields are in the same grid.
  subroutine innerproduct(nstrhalf,nsplit_total,fld1,fld2,res) ! {{{

    use lintran_constants_module, only: nst_fake

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nsplit_total !< Total number of sublayers. Set equal to nlay for non-split environment.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld1 !< One field for the inner product.
    real, dimension(nst+0*nst_fake,nstrhalf,0:1,0:nsplit_total), intent(in) :: fld2 !< The other field for the inner product.
    real, intent(inout) :: res !< Result. The result of this routine will be added to the input.

    ! Regular inner product.
    res = res + sum(fld1*fld2)

  end subroutine innerproduct ! }}}
