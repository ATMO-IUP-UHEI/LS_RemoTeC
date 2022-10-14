!> file lintran_remotec_module
!! Module file for optimization for use in RemoTeC.

!> Module for optimization for use in RemoTeC.
!!
!! This module short-circuits a few options offered by Lintran and thereby gains speed in single
!! scattering. For RemoTeC, this is relevant for single-scattering, since the K-binning applied
!! in RemoTeC causes that many more single-scattering calculations are performed than multi-scattering
!! calculations.
module lintran_remotec_module

  interface

    subroutine remotec_single_scat_without_derivatives_interface(atm,self,rint) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_atmosphere, lintran_class

      implicit none

      ! Input and output.
      type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
      type(lintran_class), intent(inout) :: self !< Main Lintran instance.
      real, dimension(nst_fake,1), intent(out) :: rint !< Main output, the reflectance.

    end subroutine remotec_single_scat_without_derivatives_interface ! }}}

    subroutine remotec_single_scat_with_derivatives_interface(atm,set,self,rint,drv) ! {{{

      use lintran_constants_module, only: nst_fake
      use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

      implicit none

      ! Input and output.
      type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
      type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
      type(lintran_class), intent(inout) :: self !< Main Lintran instance.
      real, dimension(nst_fake,1), intent(out) :: rint !< Main output, the reflectance.
      type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.

    end subroutine remotec_single_scat_with_derivatives_interface ! }}}

  end interface

contains

  !> Sets the pointers to the functions as we want them.
  !!
  !! This is setting the pointers to both the correct routine wihtout derivatives and the
  !! routine with derivatives.
  subroutine remotec_single_scat_set_nst(nst,remotec_single_scat_without_derivatives,remotec_single_scat_with_derivatives) ! {{{

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    procedure(remotec_single_scat_without_derivatives_interface), intent(out), pointer :: remotec_single_scat_without_derivatives
    procedure(remotec_single_scat_with_derivatives_interface), intent(out), pointer :: remotec_single_scat_with_derivatives

    if (nst .eq. 1) then
      remotec_single_scat_without_derivatives => remotec_single_scat_without_derivatives_nst1
      remotec_single_scat_with_derivatives => remotec_single_scat_with_derivatives_nst1
    endif
    if (nst .eq. 3) then
      remotec_single_scat_without_derivatives => remotec_single_scat_without_derivatives_nst3
      remotec_single_scat_with_derivatives => remotec_single_scat_with_derivatives_nst3
    endif
    if (nst .eq. 4) then
      remotec_single_scat_without_derivatives => remotec_single_scat_without_derivatives_nst4
      remotec_single_scat_with_derivatives => remotec_single_scat_with_derivatives_nst4
    endif

  end subroutine remotec_single_scat_set_nst ! }}}

  !> Performs scalar single-scattering calculation without derivatives with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. This is optimized code, very quick, not
  !! too easy to read and not using the TLC-module, because it is written for the case
  !! that no multi-scattering is calculated with the same atmosphere.
  subroutine remotec_single_scat_without_derivatives_nst1(atm,self,rint) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, iscat_a1
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst = 1 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance. There is one Stokes parameter and one viewing geometry.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth.
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.

    ! Aggregate result, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.

    ! Iterator.
    integer :: ilay ! Over Atmospheric layers.

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      ! TLC-parameters for atmospheric scattering:
      ! T = T_tot_0v
      ! L = L_0v
      ! C = C_bck_ssg (=ssa * phase_ssg * rotmat)
      ! The prefactor is already merged with T.

      ! Calculate the only TLC-parameter that is used twice: Once in L and the other
      ! in the aggregate of T.
      t_0v = exp(-tau/effmu_0v)

      ! We now forget effmu_0v in L, but we will apply that factor later, just to keep
      ! aggregate_pt correct for surface reflection.
      rint(ist_i,1) = rint(ist_i,1) + aggregate_pt * (1.0-t_0v) * ssa * atm%phase_ssg(iscat_a1,ilay,1)

      ! Add new transmission term to aggregate_pt.
      aggregate_pt = aggregate_pt * t_0v

    enddo

    ! Apply forgotten factor (effmu_0v) and add up surface reflection.
    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    ! The prefactor and the correct T are now correctly aggregated.
    rint(ist_i,1) = rint(ist_i,1) * effmu_0v + aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1) * atm%bdrf_ssg(ist_i,ist_i,1)

  end subroutine remotec_single_scat_without_derivatives_nst1 ! }}}

  !> Performs single-scattering calculation without derivatives for three Stokes parameters with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. Also, all three Stokes parameters are calculated.
  !! Note that polarization has no effect on I for single-scattering, so when considering
  !! polarization but only interested in I, just use scalar single-scatteirng. This is optimized
  !! code, very quick, not too easy to read and not using the TLC-module, because it is written
  !! for the case that no multi-scattering is calculated with the same atmosphere.
  subroutine remotec_single_scat_without_derivatives_nst3(atm,self,rint) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, ist_q, ist_u, iscat_a1, iscat_b1
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst = 3 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance for three Stokes parameters. There is only one viewing geometry.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.
    real :: rotmat_c ! Rotation matrix element for the cosine with direction v, used for scattering from I to Q.
    real :: rotmat_s ! Rotation matrix element for the sine with direction v, used for scattering from I to U.

    ! Aggregate result, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.
    real :: aggregate_before_phase ! Everything, including the single-scattering albedo, but not yet the phase function.

    ! Iterator.
    integer :: ilay ! Over Atmospheric layers.

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)
    rotmat_c = self%geo%rotmat_c_v_allgeo(1)
    rotmat_s = self%geo%rotmat_s_v_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      ! TLC-parameters for atmospheric scattering:
      ! T = T_tot_0v
      ! L = L_0v
      ! C = C_bck_ssg (=ssa * phase_ssg * rotmat)
      ! The prefactor is already merged with T.

      ! Calculate the only TLC-parameter that is used twice: Once in L and the other
      ! in the aggregate of T.
      t_0v = exp(-tau/effmu_0v)

      ! We now forget effmu_0v in L, but we will apply that factor later, just to keep
      ! aggregate_pt correct for surface reflection.

      ! Everything except the phase (Z) matrix elements are stored in an intermediate result
      ! that can be recycled for different Stokes parameters.
      aggregate_before_phase = aggregate_pt * (1.0-t_0v) * ssa
      rint(ist_i,1) = rint(ist_i,1) + aggregate_before_phase * atm%phase_ssg(iscat_a1,ilay,1)
      rint(ist_q,1) = rint(ist_q,1) + aggregate_before_phase * rotmat_c * atm%phase_ssg(iscat_b1,ilay,1)
      rint(ist_u,1) = rint(ist_u,1) + aggregate_before_phase * rotmat_s * atm%phase_ssg(iscat_b1,ilay,1)

      ! Add new transmission term to aggregate_pt.
      aggregate_pt = aggregate_pt * t_0v

    enddo

    ! Apply forgotten factor (effmu_0v) and add up surface reflection.
    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    ! The prefactor and the correct T are now correctly aggregated.
    rint(:,1) = rint(:,1) * effmu_0v + aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1) * atm%bdrf_ssg(1:nst,ist_i,1)

  end subroutine remotec_single_scat_without_derivatives_nst3 ! }}}

  !> Performs single-scattering calculation without derivatives for four Stokes parameters with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. Also, all three Stokes parameters are calculated.
  !! Note that polarization has no effect on I for single-scattering, so when considering
  !! polarization but only interested in I, just use scalar single-scatteirng. Note that
  !! atmospheric scattering never yields V-polarized light, so only the BDRF can provide
  !! V-polarized light. This is optimized code, very quick, not too easy to read and not
  !! using the TLC-module, because it is written for the case that no multi-scattering is
  !! calculated with the same atmosphere.
  subroutine remotec_single_scat_without_derivatives_nst4(atm,self,rint) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, ist_q, ist_u, iscat_a1, iscat_b1, pol_nst
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst = 4 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance for four Stokes parameters. There is only one viewing geometry.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth.
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.
    real :: rotmat_c ! Rotation matrix element for the cosine with direction v, used for scattering from I to Q.
    real :: rotmat_s ! Rotation matrix element for the sine with direction v, used for scattering from I to U.

    ! Aggregate result, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.
    real :: aggregate_before_phase ! Everything, including the single-scattering albedo, but not yet the phase function.

    ! Iterator.
    integer :: ilay ! Over Atmospheric layers.

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)
    rotmat_c = self%geo%rotmat_c_v_allgeo(1)
    rotmat_s = self%geo%rotmat_s_v_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      ! TLC-parameters for atmospheric scattering:
      ! T = T_tot_0v
      ! L = L_0v
      ! C = C_bck_ssg (=ssa * phase_ssg * rotmat)
      ! The prefactor is already merged with T.

      ! Calculate the only TLC-parameter that is used twice: Once in L and the other
      ! in the aggregate of T.
      t_0v = exp(-tau/effmu_0v)

      ! We now forget effmu_0v in L, but we will apply that factor later, just to keep
      ! aggregate_pt correct for surface reflection.

      ! Everything except the phase (Z) matrix elements are stored in an intermediate result
      ! that can be recycled for different Stokes parameters.
      aggregate_before_phase = aggregate_pt * (1.0-t_0v) * ssa
      rint(ist_i,1) = rint(ist_i,1) + aggregate_before_phase * atm%phase_ssg(iscat_a1,ilay,1)
      rint(ist_q,1) = rint(ist_q,1) + aggregate_before_phase * rotmat_c * atm%phase_ssg(iscat_b1,ilay,1)
      rint(ist_u,1) = rint(ist_u,1) + aggregate_before_phase * rotmat_s * atm%phase_ssg(iscat_b1,ilay,1)

      ! Add new transmission term to aggregate_pt.
      aggregate_pt = aggregate_pt * t_0v

    enddo

    ! Apply forgotten factor (effmu_0v) and add up surface reflection.
    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    ! The prefactor and the correct T are now correctly aggregated.
    rint(:,1) = rint(:,1) * effmu_0v + aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1) * atm%bdrf_ssg(1:nst,ist_i,1)

  end subroutine remotec_single_scat_without_derivatives_nst4 ! }}}

  !> Performs scalar single-scattering calculation with derivatives with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. This is optimized code, very quick, not
  !! too easy to read and not using the TLC-module, because it is written for the case
  !! that no multi-scattering is calculated with the same atmosphere.
  subroutine remotec_single_scat_with_derivatives_nst1(atm,set,self,rint,drv) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, iscat_a1
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst = 1 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance. There is one Stokes parameter and one viewing geometry.
    type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Derivative with respect to tau and ssa.
    real, dimension(self%nlay) :: drint_tau
    real, dimension(self%nlay) :: drint_ssa

    ! Square of tau, for chain rule.
    real, dimension(self%nlay) :: tau_squared

    ! Dimension size.
    integer :: nlay_deriv_ph

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth.
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.
    real :: c_bck_ssg ! Single-scattering albedo times phase function for single-scattering geometry.

    ! Aggregate results, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.
    real :: aggregate_ptl ! Prefactor, transmission and layer-integral.
    real :: rint_lay ! Reflectance by scattering in this layer.

    ! Iterators.
    integer :: ilay ! Over Atmospheric layers.
    integer :: ilay_ph ! For phase derivatives.

    ! Administration of layer-index for derivatives with respect to phase function.
    nlay_deriv_ph = size(set%ilay_deriv_ph)
    ilay_ph = 1

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      tau_squared(ilay) = tau**2.0

      ! Calculate TLC-parameters for this layer.
      t_0v = exp(-tau/effmu_0v)
      ! Aggregate is prefactor * T * L.
      aggregate_ptl = aggregate_pt * effmu_0v * (1.0-t_0v)

      drint_ssa(ilay) = aggregate_ptl * atm%phase_ssg(iscat_a1,ilay,1)
      ! Get C-parameter, Z times ssa.
      c_bck_ssg = atm%phase_ssg(iscat_a1,ilay,1) * ssa

      rint_lay = aggregate_ptl * c_bck_ssg
      rint(ist_i,1) = rint(ist_i,1) + rint_lay
      ! Derivative with respect to single-scattering albedo.
      ! Derivative with respect to tau in layers above.
      drint_tau(1:ilay-1) = drint_tau(1:ilay-1) -1.0/effmu_0v * rint_lay
      if (ilay_ph .le. nlay_deriv_ph) then
        if (set%ilay_deriv_ph(ilay_ph) .eq. ilay) then
          drv%drint_phase_ssg(iscat_a1,ilay_ph,ist_i,1) = aggregate_ptl * ssa
          ilay_ph = ilay_ph + 1
        endif
      endif

      ! Define derivative of tau.
      ! Add new transmission term to aggregate_pt.
      ! For the currect layer, aggregate_pt will be equal to prefactor * T * dL/dtau.
      aggregate_pt = aggregate_pt * t_0v

      drint_tau(ilay) = aggregate_pt * c_bck_ssg ! Note that P*T is temporarily P*T*dL/dtau.

    enddo

    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    drv%drint_bdrf_ssg(ist_i,ist_i,ist_i,1) = aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1)
    rint_lay = drv%drint_bdrf_ssg(ist_i,ist_i,ist_i,1) * atm%bdrf_ssg(ist_i,ist_i,1)
    rint(ist_i,1) = rint(ist_i,1) + rint_lay

    drint_tau = drint_tau - rint_lay / effmu_0v

    ! Convert coordinates.
    drv%drint_taus(:,ist_i,1) = drint_ssa * atm%taua / tau_squared + drint_tau
    drv%drint_taua(:,ist_i,1) = -drint_ssa * atm%taus / tau_squared + drint_tau

  end subroutine remotec_single_scat_with_derivatives_nst1 ! }}}

  !> Performs single-scattering calculation with derivatives for three Stokes parameters with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. Also, all three Stokes parameters are calculated.
  !! Note that polarization has no effect on I for single-scattering, so when considering
  !! polarization but only interested in I, just use scalar single-scatteirng. This is optimized
  !! code, very quick, not too easy to read and not using the TLC-module, because it is written
  !! for the case that no multi-scattering is calculated with the same atmosphere.
  subroutine remotec_single_scat_with_derivatives_nst3(atm,set,self,rint,drv) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, ist_q, ist_u, iscat_a1, iscat_b1
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst = 3 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance for three Stokes parameters. There is only one viewing geometry.
    type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Derivative with respect to tau and ssa.
    real, dimension(self%nlay,nst) :: drint_tau
    real, dimension(self%nlay,nst) :: drint_ssa

    ! Square of tau, for chain rule.
    real, dimension(self%nlay) :: tau_squared

    ! Dimension size.
    integer :: nlay_deriv_ph

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth.
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.
    real :: rotmat_c ! Rotation matrix element for the cosine with direction v, used for scattering from I to Q.
    real :: rotmat_s ! Rotation matrix element for the sine with direction v, used for scattering from I to U.
    real, dimension(nst) :: c_bck_ssg ! Single-scattering albedo times phase function for single-scattering geometry.

    ! Aggregate results, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.
    real :: aggregate_ptl ! Prefactor, transmission and layer-integral.
    real, dimension(nst) :: rint_lay ! Reflectance by scattering in this layer.

    ! Iterators.
    integer :: ilay ! Over Atmospheric layers.
    integer :: ilay_ph ! For phase derivatives.
    integer :: ist ! Over Stokes parameters.

    ! Administration of layer-index for derivatives with respect to phase function.
    nlay_deriv_ph = size(set%ilay_deriv_ph)
    ilay_ph = 1

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)
    rotmat_c = self%geo%rotmat_c_v_allgeo(1)
    rotmat_s = self%geo%rotmat_s_v_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      tau_squared(ilay) = tau**2.0

      ! Calculate TLC-parameters for this layer.
      t_0v = exp(-tau/effmu_0v)
      ! Aggregate is prefactor * T * L.
      aggregate_ptl = aggregate_pt * effmu_0v * (1.0-t_0v)

      c_bck_ssg(ist_i) = atm%phase_ssg(iscat_a1,ilay,1)
      c_bck_ssg(ist_q) = rotmat_c * atm%phase_ssg(iscat_b1,ilay,1)
      c_bck_ssg(ist_u) = rotmat_s * atm%phase_ssg(iscat_b1,ilay,1)
      ! Now, C is still Z. We need a multiplication with the single-scattering albedo.
      drint_ssa(ilay,:) = aggregate_ptl * c_bck_ssg ! C is still Z.
      ! Update Z to C.
      c_bck_ssg = c_bck_ssg * ssa

      rint_lay = aggregate_ptl * c_bck_ssg
      rint(:,1) = rint(:,1) + rint_lay
      ! Derivative with respect to single-scattering albedo.
      ! Derivative with respect to tau in layers above.
      do ist = 1,nst
        drint_tau(1:ilay-1,ist) = drint_tau(1:ilay-1,ist) -1.0/effmu_0v * rint_lay(ist)
      enddo
      if (ilay_ph .le. nlay_deriv_ph) then
        if (set%ilay_deriv_ph(ilay_ph) .eq. ilay) then
          drv%drint_phase_ssg(iscat_a1,ilay_ph,ist_i,1) = aggregate_ptl * ssa
          drv%drint_phase_ssg(iscat_b1,ilay_ph,ist_q,1) = aggregate_ptl * ssa * rotmat_c
          drv%drint_phase_ssg(iscat_b1,ilay_ph,ist_u,1) = aggregate_ptl * ssa * rotmat_s
          ilay_ph = ilay_ph + 1
        endif
      endif

      ! Define derivative of tau.
      ! Add new transmission term to aggregate_pt.
      ! For the currect layer, aggregate_pt will be equal to prefactor * T * dL/dtau.
      aggregate_pt = aggregate_pt * t_0v

      drint_tau(ilay,:) = aggregate_pt * c_bck_ssg ! Note that P*T is temporarily P*T*dL/dtau.

    enddo

    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    do ist = 1,nst
      drv%drint_bdrf_ssg(ist,ist_i,ist,1) = aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1)
    enddo
    rint_lay = drv%drint_bdrf_ssg(ist_i,ist_i,ist_i,1) * atm%bdrf_ssg(1:nst,ist_i,1)
    rint(:,1) = rint(:,1) + rint_lay
    ! Misuse rint_lay, because a loop is needed for this 1D-2D vector operation.
    rint_lay = -rint_lay/effmu_0v
    do ist = 1,nst
      drint_tau(:,ist) = drint_tau(:,ist) + rint_lay(ist)
    enddo

    ! Convert coordinates.
    do ist = 1,nst
      drv%drint_taus(:,ist,1) = drint_ssa(:,ist) * atm%taua / tau_squared + drint_tau(:,ist)
      drv%drint_taua(:,ist,1) = -drint_ssa(:,ist) * atm%taus / tau_squared + drint_tau(:,ist)
    enddo

  end subroutine remotec_single_scat_with_derivatives_nst3 ! }}}

  !> Performs single-scattering calculation with derivatives for four Stokes parameters with settings from RemoTeC.
  !!
  !! These settings are a plane-parallel atmosphere, no fluorescent emission and,
  !! no thermal emission and normalized sun. Also, all three Stokes parameters are calculated.
  !! Note that polarization has no effect on I for single-scattering, so when considering
  !! polarization but only interested in I, just use scalar single-scatteirng. Note that
  !! atmospheric scattering never yields V-polarized light, so only the BDRF can provide
  !! V-polarized light. This is optimized code, very quick, not too easy to read and not
  !! using the TLC-module, because it is written for the case that no multi-scattering is
  !! calculated with the same atmosphere.
  subroutine remotec_single_scat_with_derivatives_nst4(atm,set,self,rint,drv) ! {{{

    use lintran_constants_module, only: nst_fake, ist_i, ist_q, ist_u, ist_v, iscat_a1, iscat_b1
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

    implicit none

    ! Needed for the dimensions.
    integer, parameter :: nst_atm = 3 ! Number of Stokes parameters reached by atmospheric scattering.
    integer, parameter :: nst = 4 ! Number of Stokes parameters.

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst+0*nst_fake,1), intent(out) :: rint !< Main output, the reflectance for four Stokes parameters. There is only one viewing geometry.
    type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.

    ! This single-scattering routine is limited to plane-parallel geometry and no delta-M
    ! and only one viewing geometry, and no emission.

    ! Derivative with respect to tau and ssa.
    real, dimension(self%nlay,nst_atm) :: drint_tau
    real, dimension(self%nlay,nst_atm) :: drint_ssa
    real, dimension(self%nlay) :: drint_v_tau

    ! Square of tau, for chain rule.
    real, dimension(self%nlay) :: tau_squared

    ! Dimension size.
    integer :: nlay_deriv_ph

    ! Atmospheric parmaeters.
    real :: tau ! Optical depth.
    real :: ssa ! Single-scattering albedo.

    ! TLC-parameters, used once per layer.
    real :: effmu_0v ! Effective direction cosine from geometry module (copied for local access).
    real :: t_0v ! Transmission through a layer in both solar and viewing direction.
    real :: rotmat_c ! Rotation matrix element for the cosine with direction v, used for scattering from I to Q.
    real :: rotmat_s ! Rotation matrix element for the sine with direction v, used for scattering from I to U.
    real, dimension(nst_atm) :: c_bck_ssg ! Single-scattering albedo times phase function for single-scattering geometry.

    ! Aggregate results, meaning changes, see comments whenever it is assigned.
    real :: aggregate_pt ! Prefactor and transmission.
    real :: aggregate_ptl ! Prefactor, transmission and layer-integral.
    real, dimension(nst) :: rint_lay ! Reflectance by scattering in this layer.

    ! Iterators.
    integer :: ilay ! Over Atmospheric layers.
    integer :: ilay_ph ! For phase derivatives.
    integer :: ist ! Over Stokes parameters.

    ! Administration of layer-index for derivatives with respect to phase function.
    nlay_deriv_ph = size(set%ilay_deriv_ph)
    ilay_ph = 1

    ! Copy often-used number to the local stack.
    effmu_0v = self%geo%effmu_0v_ppg_allgeo(1)
    rotmat_c = self%geo%rotmat_c_v_allgeo(1)
    rotmat_s = self%geo%rotmat_s_v_allgeo(1)

    ! Initialize results to zero.
    rint = 0.0

    ! TLC-parameters for single-scattering:
    ! T = T_tot_0v
    ! L = L_0v
    ! C = ssa * Z_ss

    ! The prefactor is prefactor_dir and exists in grd as prefactor_dir.
    ! We will also add the degeneracy factor, which exists in geometry as
    ! degenerate_ssg_allgeo(1).

    ! Transmission to top of the current layer. That starts with just one, so the
    ! aggregate of the prefactor with the transmission term starts with just the prefactor.
    aggregate_pt = self%grd%prefactor_dir * self%geo%degenerate_ssg_allgeo(1)
    do ilay = 1,self%nlay

      ! Get optical properties.
      tau = atm%taus(ilay) + atm%taua(ilay)
      ssa = atm%taus(ilay) / tau

      tau_squared(ilay) = tau**2.0

      ! Calculate TLC-parameters for this layer.
      t_0v = exp(-tau/effmu_0v)
      ! Aggregate is prefactor * T * L.
      aggregate_ptl = aggregate_pt * effmu_0v * (1.0-t_0v)

      c_bck_ssg(ist_i) = atm%phase_ssg(iscat_a1,ilay,1)
      c_bck_ssg(ist_q) = rotmat_c * atm%phase_ssg(iscat_b1,ilay,1)
      c_bck_ssg(ist_u) = rotmat_s * atm%phase_ssg(iscat_b1,ilay,1)
      ! Now, C is still Z. We need a multiplication with the single-scattering albedo.
      drint_ssa(ilay,:) = aggregate_ptl * c_bck_ssg(:) ! C is still Z.
      ! Update Z to C.
      c_bck_ssg = c_bck_ssg * ssa

      rint_lay(1:nst_atm) = aggregate_ptl * c_bck_ssg
      rint(1:nst_atm,1) = rint(1:nst_atm,1) + rint_lay(1:nst_atm)
      ! Derivative with respect to single-scattering albedo.
      ! Derivative with respect to tau in layers above.
      do ist = 1,nst_atm
        drint_tau(1:ilay-1,ist) = drint_tau(1:ilay-1,ist) -1.0/effmu_0v * rint_lay(ist)
      enddo
      if (ilay_ph .le. nlay_deriv_ph) then
        if (set%ilay_deriv_ph(ilay_ph) .eq. ilay) then
          drv%drint_phase_ssg(iscat_a1,ilay_ph,ist_i,1) = aggregate_ptl * ssa
          drv%drint_phase_ssg(iscat_b1,ilay_ph,ist_q,1) = aggregate_ptl * ssa * rotmat_c
          drv%drint_phase_ssg(iscat_b1,ilay_ph,ist_u,1) = aggregate_ptl * ssa * rotmat_s
          ilay_ph = ilay_ph + 1
        endif
      endif

      ! Define derivative of tau.
      ! Add new transmission term to aggregate_pt.
      ! For the currect layer, aggregate_pt will be equal to prefactor * T * dL/dtau.
      aggregate_pt = aggregate_pt * t_0v

      drint_tau(ilay,:) = aggregate_pt * c_bck_ssg ! Note that P*T is temporarily P*T*dL/dtau.

    enddo

    ! TLC-parmaeters for surface reflection:
    ! T = T_tot_0v
    ! L = 1
    ! C = BDRF * diff_c_bdrf
    do ist = 1,nst
      drv%drint_bdrf_ssg(ist,ist_i,ist,1) = aggregate_pt * self%geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(1)
    enddo
    rint_lay = drv%drint_bdrf_ssg(ist_i,ist_i,ist_i,1) * atm%bdrf_ssg(:,ist_i,1)
    rint(:,1) = rint(:,1) + rint_lay
    ! Misuse rint_lay, because a loop is needed for this 1D-2D vector operation.
    rint_lay = -rint_lay/effmu_0v
    do ist = 1,nst_atm
      drint_tau(:,ist) = drint_tau(:,ist) + rint_lay(ist)
    enddo
    drint_v_tau = rint_lay(ist_v)

    ! Convert coordinates.
    do ist = 1,nst_atm
      drv%drint_taus(:,ist,1) = drint_ssa(:,ist) * atm%taua / tau_squared + drint_tau(:,ist)
      drv%drint_taua(:,ist,1) = -drint_ssa(:,ist) * atm%taus / tau_squared + drint_tau(:,ist)
    enddo
    drv%drint_taus(:,ist_v,1) = drint_v_tau
    drv%drint_taua(:,ist_v,1) = drint_v_tau

  end subroutine remotec_single_scat_with_derivatives_nst4 ! }}}

end module lintran_remotec_module
