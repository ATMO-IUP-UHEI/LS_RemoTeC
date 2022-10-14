!> \file lintran_geometry_module.f90
!! Module file for the solar and viewing geometry.

!> Module for the solar and viewing geometry.
!!
!! This modules handles everything that does depend on the solar or viewing
!! geometry, but not on the optical properties of the atmosphere. The idea is
!! that radiative transfer is usually calculated for several wavelengths. For
!! different wavelengths, the same solar and viewing geometry applies, but the
!! optical properties of the atmosphere and the surface may be different. Therefore,
!! the calculations performed in the geometry module can be recycled for different
!! wavelengths in the same geometry.
module lintran_geometry_module

  implicit none

contains

  !> Constructs space for all properties that depend just on the geometry.
  !!
  !! No actual geometry is calculated. The space required to perform calculations
  !! for a geometry does not change, so the allocations are performed before any
  !! geometry is provided. It is possible to include different viewing geometries
  !! in one calculations. The number of viewing geometries should be fixed.
  subroutine geometry_init(nst,nstrhalf,nleg,nlay,ngeo,geo,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, pol_nst
    use lintran_types_module, only: lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    type(lintran_geometry_class), intent(out) :: geo !< Internal geometry structure.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Viewing zenith angles.
    allocate(geo%mu_v_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%azi_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Visibility ranges, dependent on viewing geometry. One of the solar
    ! geometry is scalar and need not be allocated. The only possibility for a
    ! non-trivial content is a pseudo-spherical geometry with solar or viewing
    ! zenith angle above 90 degress.
    allocate(geo%visible_v_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%visible_0v_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Sensitivity limits, for removing zero-calculations. The only possibility for a
    ! non-trivial content is a pseudo-spherical geometry with solar or viewing
    ! zenith angle above 90 degress.
    allocate(geo%ilay_sensitive_0(0:nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%ilay_sensitive_v_allgeo(0:nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Degeneracy array for all Stokes parameters. These degeneracy arrays also contain
    ! a division by the cosine of the viewing angle. Thereby, some other terms become
    ! geometry-independent. For single-scattering, a degeneracy with just the division
    ! by the cosine of the viewing angle is needed. That is Stokes-independent.
    allocate(geo%degenerate_allgeo(0:nleg,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%degenerate_ssg_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Effective direction cosines for plane-parallel geometry.
    allocate(geo%effmu_0v_ppg_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%effmu_p0_ppg(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%effmu_m0_ppg(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%effmu_pv_ppg_allgeo(nstrhalf,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%effmu_mv_ppg_allgeo(nstrhalf,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Relative derivatives of transmission terms with respect to the optical
    ! depth. These arrays are used to construct the transmission terms themselves
    ! as well, not only for the derivatives. The contents of the arrays are
    ! trivial for plane-parallel atmosphere, namely 0 and -1/mu_x, where x is sun,
    ! instrument or a combination of both. For pseudo-spherical calculations, these
    ! values become minus path lengths.
    ! First dimension: Perturbed layer.
    ! Second dimension: Interface until which is transmitted.
    allocate(geo%reldiff_t_tot_0_tau(nlay,0:nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%reldiff_t_tot_v_tau_allgeo(nlay,0:nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%reldiff_t_tot_0v_tau_allgeo(nlay,0:nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Single-scattering rotation matrix. Note that the solar (0) side also
    ! depends on the viewing geometry because of the angle of single scattering.
    if (pol_nst(nst)) then ! Only for vector.
      allocate(geo%rotmat_c_0_allgeo(ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%rotmat_c_v_allgeo(ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%rotmat_s_0_allgeo(ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%rotmat_s_v_allgeo(ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 4 ! Fake allocations.
    endif

    ! Genelarized spherical functions on sun or instrument geometry.
    allocate(geo%gsf_0_0(0:nleg,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%gsf_0_v_allgeo(0:nleg,0:nleg,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (pol_nst(nst)) then ! Only for vector.
      allocate(geo%gsf_p_0(0:nleg,0:nleg),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%gsf_p_v_allgeo(0:nleg,0:nleg,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%gsf_m_0(0:nleg,0:nleg),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(geo%gsf_m_v_allgeo(0:nleg,0:nleg,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,geo,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 4 ! Fake allocations.
    endif

    ! Derivatives with respect to the Bi-directional reflection function (BDRF).
    ! The Stokes parameter dimension vanishes in thes derivatives. This derivative
    ! contains just some geometrical terms.
    allocate(geo%diff_c_bck_srf_0_bdrf(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%diff_c_bck_srf_v_bdrf_allgeo(nstrhalf,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! For finite-element C-parameters that are scalar besides their dependence
    ! on the viewing geometry.
    allocate(geo%diff_c_bck_srf_0v_elem_bdrf_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(geo%diff_c_emi_v_elem_emi_allgeo(ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,geo,stat)
      return
    endif
    progression = progression + 1 ! }}}

  end subroutine geometry_init ! }}}

  !> Set the geometry and calculate related variables, except for transmission-related parameters.
  !!
  !! The solar and viewing geometries are set. All gemoetry-dependent except transmission-related
  !! terms are calculated as well. The transmission terms all depend on whether plane-parallel
  !! or pseudo-spherical geometry is applied. For numeric reasons, the viewing zenith angle cannot
  !! be too close to zero, so that is asserted and the viewing zenith angles may possibly be changed
  !! slightly.
  subroutine geometry_set_general(nst,nleg,ngeo,mu_0,mu_v,azi,grd,geo,stat) ! {{{

    use lintran_constants_module, only: ist_i, ist_q, ist_u, ist_v, ist_sun, pi, minimum_mudiff, warningflag_instableangle, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    real, intent(in) :: mu_0 !< Cosine of solar zenith angle.
    real, dimension(ngeo), intent(in) :: mu_v !< Cosines of viewing zenith angles.
    real, dimension(ngeo), intent(in) :: azi !< Azimuthal differences in radians, defined as in De Haan et al. (1987).
    type(lintran_grid_class), intent(in) :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: m ! Fourier number.
    integer :: igeo ! Over viewing geometries.

    ! Cosine of single-scattering angle
    real :: mu_ssg

    ! Save the geometry in the structure, so it is accessible for other routines in
    ! the geometry module.
    geo%mu_0 = mu_0
    geo%mu_v_allgeo = mu_v
    ! Only handy for the Nakajima wrapper. Inside Lintran, the azimuth need not be saved.
    geo%azi_allgeo = azi
    ! Protect mu_v for numbers too close to zero, because it appears in the degeneracy array.
    ! That does not depend on the pseudo-spherical or plane parallel geometry.
    do igeo = 1,ngeo
      if (abs(geo%mu_v_allgeo(igeo)) .lt. minimum_mudiff) then
        if (geo%mu_v_allgeo(igeo) .lt. 0.0) then
          geo%mu_v_allgeo(igeo) = -minimum_mudiff
        else
          geo%mu_v_allgeo(igeo) = minimum_mudiff
        endif
      endif
      if (mu_v(igeo) .ne. geo%mu_v_allgeo(igeo)) then
        stat = or(stat,warningflag_instableangle)
      endif
    enddo

    ! Pre-calculate degeneracy.The normal degeneracy is one for m=0.
    ! For m>0, both m and -m are calculated together and their azimuthal
    ! factors are exp(i*m*azi) and exp(i*(-m)*azi). Apart from that,
    ! m and -m have the same properties. These two azimuthal factors
    ! added up are 2*cos(m*azi).
    ! For U and V, there is an additional term i, which means that we end up
    ! with sines instead of cosines. And we have a zero for m=0. Moreover,
    ! for U and V, we have to undo a sign switch that we applied for upward
    ! radiation. This sign switch is done internally to make the phase matrices
    ! more symmetric, but for the observable, we have to return to the sign
    ! conventions.
    ! Though the sun is hardcoded to radiate light, not polarization, this is
    ! always the case. Still, we evaluate the if-clause for this, because Lintran
    ! is still mathematically healthy if the sun radiated polarization. But the
    ! degeneracy factors become different if it radiated U or V. This means that
    ! Lintran can be hacked to a polarization-radiating sun, for instance for a
    ! manual adjoint if somehow needed.
    ! Furthermore, there is a division by mu_v.
    ! In single-scattering geometry, it is only the division by mu_v.
    if (ist_sun .lt. ist_u) then ! I or Q. It is I, so this is true.
      do igeo = 1,ngeo
        geo%degenerate_allgeo(:,ist_i,igeo) = (/1.0,(2.0*cos(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        if (pol_nst(nst)) then ! Q and U.
          geo%degenerate_allgeo(:,ist_q,igeo) = (/1.0,(2.0*cos(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
          geo%degenerate_allgeo(:,ist_u,igeo) = (/0.0,(-2.0*sin(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        endif
        if (fullpol_nst(nst)) then ! V.
          geo%degenerate_allgeo(:,ist_v,igeo) = (/0.0,(-2.0*sin(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        endif
      enddo
    else
      do igeo = 1,ngeo
        geo%degenerate_allgeo(:,ist_i,igeo) = (/0.0,(-2.0*sin(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        if (pol_nst(nst)) then ! Q and U.
          geo%degenerate_allgeo(:,ist_q,igeo) = (/0.0,(-2.0*sin(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
          geo%degenerate_allgeo(:,ist_u,igeo) = (/-1.0,(-2.0*cos(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        endif
        if (fullpol_nst(nst)) then ! V.
          geo%degenerate_allgeo(:,ist_v,igeo) = (/-1.0,(-2.0*cos(m*azi(igeo)),m=1,nleg)/) / abs(geo%mu_v_allgeo(igeo))
        endif
      enddo
    endif
    ! For single-scattering, there is only a division by mu_v.
    geo%degenerate_ssg_allgeo = 1.0 / abs(geo%mu_v_allgeo) ! Array operation over geometries.

    ! Calculate single-scattering rotation matrix elements, using cos(2a) = 2*cos^2(a) - 1
    if (pol_nst(nst)) then
      do igeo = 1,ngeo
        mu_ssg = -geo%mu_0*geo%mu_v_allgeo(igeo) + sqrt(1.0-geo%mu_0**2.0)*sqrt(1.0-geo%mu_v_allgeo(igeo)**2.0)*cos(azi(igeo))
        ! If the single-scattering angle cosine is minus one, we have a problem, but in fact,
        ! the rotation matrices should not matter, because the phase function should be zero.
        if (mu_ssg .le. -1.0) then
          geo%rotmat_c_0_allgeo(igeo) = 0.0
          geo%rotmat_c_v_allgeo(igeo) = 0.0
        else
          ! Now, there is an additional problem for nadir (zero zenith angle for sun or
          ! instrument). That also gives a zero divided by zero quotient. The limit of
          ! the quotient is the cosine of the azimuth angle squared. And it is relevant,
          ! because the phase function is not zero. The scattering angle is just some angle.
          if (geo%mu_0 .ge. 1.0) then
            geo%rotmat_c_0_allgeo(igeo) = 2.0 * cos(azi(igeo))**2.0 - 1.0
          else
            geo%rotmat_c_0_allgeo(igeo) = 2.0 * (geo%mu_v_allgeo(igeo) + geo%mu_0*mu_ssg)**2.0 / ((1.0-mu_ssg**2.0) * (1.0-geo%mu_0**2.0)) - 1.0
          endif
          if (geo%mu_v_allgeo(igeo) .ge. 1.0) then
            geo%rotmat_c_v_allgeo(igeo) = 2.0 * cos(azi(igeo))**2.0 - 1.0
          else
            geo%rotmat_c_v_allgeo(igeo) = 2.0 * (geo%mu_0 + geo%mu_v_allgeo(igeo)*mu_ssg)**2.0 / ((1.0-mu_ssg**2.0) * (1.0-geo%mu_v_allgeo(igeo)**2.0)) - 1.0
          endif
        endif
        ! Protection against ones or minus ones that numerically become a little bit larger.
        ! This is not a physical cutoff.
        if (abs(geo%rotmat_c_0_allgeo(igeo)) .ge. 1.0) then
          geo%rotmat_s_0_allgeo(igeo) = 0.0
        else
          geo%rotmat_s_0_allgeo(igeo) = sqrt(1.0-geo%rotmat_c_0_allgeo(igeo)**2.0)
        endif
        if (abs(geo%rotmat_c_v_allgeo(igeo)) .ge. 1.0) then
          geo%rotmat_s_v_allgeo(igeo) = 0.0
        else
          geo%rotmat_s_v_allgeo(igeo) = sqrt(1.0-geo%rotmat_c_v_allgeo(igeo)**2.0)
        endif
        ! Make sure that the sine term gets the right sign. Note that the signs are a bit messed
        ! up because mu_0 and mu_v are defined postive and the formulas from Hovenier assume a
        ! different sign for mu_0 and mu_v, because one is upward and the other is downward.
        ! This issue created a sign switch for rotmat_s_0, and that is corrected to mix up
        ! the .lt. and the .gt. in the next lines.
        if (geo%mu_0 .ge. 1.0) then
          geo%rotmat_s_0_allgeo(igeo) = -2.0 * sin(azi(igeo)) * cos(azi(igeo))
        else
          if ((modulo(azi(igeo),2.0*pi) .gt. pi) .eqv. (geo%mu_v_allgeo(igeo) + geo%mu_0*mu_ssg .lt. 0.0)) geo%rotmat_s_0_allgeo(igeo) = -geo%rotmat_s_0_allgeo(igeo)
        endif
        if (geo%mu_v_allgeo(igeo) .ge. 1.0) then
          geo%rotmat_s_v_allgeo(igeo) = 2.0 * sin(azi(igeo)) * cos(azi(igeo))
        else
          if ((modulo(azi(igeo),2.0*pi) .gt. pi) .eqv. (geo%mu_0 + geo%mu_v_allgeo(igeo)*mu_ssg .gt. 0.0)) geo%rotmat_s_v_allgeo(igeo) = -geo%rotmat_s_v_allgeo(igeo)
        endif
      enddo
    endif

    ! Calculate phase derivatives for solar and viewing geometries.
    do m = 0,nleg
      call calculate_geometry_generalized_spherical_functions(nst,nleg,ngeo,m,geo)
    enddo

    ! Derivatives with respect to Bi-directional reflection function (BDRF)
    ! and emission. Note that absolute values of direction cosines are not
    ! necessary. If the direction cosine is negative, the surface is invisible.

    ! Surface reflection.
    ! C = BDRF * 4*mu_src*mu_dest

    ! Uniform emission field.
    ! C = emission * 4*pi*mu

    geo%diff_c_bck_srf_0_bdrf = 4.0 * grd%mu * geo%mu_0
    ! Do-loop required because 1D-2D chaos.
    do igeo = 1,ngeo
      geo%diff_c_bck_srf_v_bdrf_allgeo(:,igeo) = 4.0 * grd%mu * geo%mu_v_allgeo(igeo)
    enddo

    ! Finite-element derivatives of C with respect to albedo or emissivity.
    ! Array operation over geometries.
    geo%diff_c_bck_srf_0v_elem_bdrf_allgeo = 4.0 * geo%mu_0 * geo%mu_v_allgeo
    geo%diff_c_emi_v_elem_emi_allgeo = 4.0*pi * geo%mu_v_allgeo

  end subroutine geometry_set_general ! }}}

  !> Applies plane-parallel geometry.
  !!
  !! This routine calculates the path length factors. Because this is plane-parallel geometry,
  !! these path length factors are quite simple. With a minus signs, those are the relative
  !! derivatives of the cumulative transmission terms. From the relative derivatives of the
  !! cumulative transmission terms, all other transmission terms and related intregrals,
  !! in which the solar or viewing geometry is involved, will be calculated when the optical
  !! properties of the atmosphere are known. The flag for plane-parallel geometry is saved
  !! for optimization.
  subroutine geometry_set_plane_parallel(nstrhalf,nlay,ngeo,grd,geo,stat) ! {{{

    use lintran_constants_module, only: minimum_mudiff, warningflag_instableangle
    use lintran_types_module, only: lintran_grid_class, lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    type(lintran_grid_class), intent(in) :: grd !< Internal grid structure.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: ilay ! Over atmospheric layers, actually interfaces.
    integer :: igeo ! Over geometries.
    integer :: istr ! Over streams.

    ! For stability check warnings.
    real :: mu_x_original

    ! Check for solar or viewing zenith angles that are too close to the stream grid.
    ! Also, those angles should not be too close to zero (90 degrees, zero cosine).

    ! Solar zenith angle.
    mu_x_original = geo%mu_0
    ! Check for plus or minus internal streams.
    do istr = 1,nstrhalf
      if (abs(geo%mu_0 - grd%mu(istr)) .lt. minimum_mudiff) then
        if (geo%mu_0 .lt. grd%mu(istr)) then
          geo%mu_0 = grd%mu(istr) - minimum_mudiff
        else
          geo%mu_0 = grd%mu(istr) + minimum_mudiff
        endif
        exit
      endif
      if (abs(geo%mu_0 + grd%mu(istr)) .lt. minimum_mudiff) then
        if (geo%mu_0 .lt. -grd%mu(istr)) then
          geo%mu_0 = -grd%mu(istr) - minimum_mudiff
        else
          geo%mu_0 = -grd%mu(istr) + minimum_mudiff
        endif
        exit
      endif
    enddo
    ! Check for zero.
    if (abs(geo%mu_0) .lt. minimum_mudiff) then
      if (geo%mu_0 .lt. 0.0) then
        geo%mu_0 = -minimum_mudiff
      else
        geo%mu_0 = minimum_mudiff
      endif
    endif
    ! Raise warning if something is changed.
    if (mu_x_original .ne. geo%mu_0) then
      stat = or(stat,warningflag_instableangle)
    endif
    ! Effective viewing zenith angles.
    do igeo = 1,ngeo
      mu_x_original = geo%mu_v_allgeo(igeo)
      ! Initial check for the negative solar direction. Hereafter, this check
      ! is repeated each time the effective angle is moved. If this movement here
      ! enters the instability domain of another stream, it will be fixed there.
      ! For plane-parallel geometry, it is a bit strange, because a viewing direction
      ! in the negative solar direction gives zero signal, but there is a numeric
      ! instability in one of the terms.
      if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) then
        if (geo%mu_v_allgeo(igeo) .lt. -geo%mu_0) then
          geo%mu_v_allgeo(igeo) = -geo%mu_0 - minimum_mudiff
        else
          geo%mu_v_allgeo(igeo) = -geo%mu_0 + minimum_mudiff
        endif
      endif
      ! Check for plus or minus internal streams.
      do istr = 1,nstrhalf
        if (abs(geo%mu_v_allgeo(igeo) - grd%mu(istr)) .lt. minimum_mudiff) then
          if (geo%mu_v_allgeo(igeo) .lt. grd%mu(istr)) then
            geo%mu_v_allgeo(igeo) = grd%mu(istr) - minimum_mudiff
            ! Very unlikely: if mu_v moves close to -mu_0, there will be additional
            ! instability. While loop can be done at most twice. Therefore, two if-statements
            ! replace the while loop, because there is a razzia for while loops.
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
          else
            geo%mu_v_allgeo(igeo) = grd%mu(istr) + minimum_mudiff
            ! Very unlikely: if mu_v moves close to -mu_0, there will be additional
            ! instability. While loop can be done at most twice. Therefore, two if-statements
            ! replace the while loop, because there is a razzia for while loops.
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
          endif
          exit
        endif
        if (abs(geo%mu_v_allgeo(igeo) + grd%mu(istr)) .lt. minimum_mudiff) then
          if (geo%mu_v_allgeo(igeo) .lt. -grd%mu(istr)) then
            geo%mu_v_allgeo(igeo) = -grd%mu(istr) - minimum_mudiff
            ! Very unlikely: if mu_v moves close to -mu_0, there will be additional
            ! instability. While loop can be done at most twice. Therefore, two if-statements
            ! replace the while loop, because there is a razzia for while loops.
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
          else
            geo%mu_v_allgeo(igeo) = -grd%mu(istr) + minimum_mudiff
            ! Very unlikely: if mu_v moves close to -mu_0, there will be additional
            ! instability. While loop can be done at most twice. Therefore, two if-statements
            ! replace the while loop, because there is a razzia for while loops.
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
            if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
          endif
          exit
        endif
      enddo
      ! Check for zero. This was already done in the general initialization for the
      ! degeneracy terms, but it is error prone to remove this code, because if the
      ! assertion at the degeneracy is changed, the check here will be incomplete.
      if (abs(geo%mu_v_allgeo(igeo)) .lt. minimum_mudiff) then
        if (geo%mu_v_allgeo(igeo) .lt. 0.0) then
          geo%mu_v_allgeo(igeo) = -minimum_mudiff
          ! Very unlikely: if mu_v is close to -mu_0, there will be additional
          ! instability. While loop can be done at most twice. Therefore, two if-statements
          ! replace the while loop, because there is a razzia for while loops.
          if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
          if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) - minimum_mudiff
        else
          geo%mu_v_allgeo(igeo) = minimum_mudiff
          ! Very unlikely: if mu_v is close to -mu_0, there will be additional
          ! instability. While loop can be done at most twice. Therefore, two if-statements
          ! replace the while loop, because there is a razzia for while loops.
          if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
          if (abs(geo%mu_v_allgeo(igeo)+geo%mu_0) .lt. minimum_mudiff) geo%mu_v_allgeo(igeo) = geo%mu_v_allgeo(igeo) + minimum_mudiff
        endif
      endif
      ! Raise warning if something is changed.
      if (mu_x_original .ne. geo%mu_v_allgeo(igeo)) then
        stat = or(stat,warningflag_instableangle)
      endif
    enddo

    ! Relative derivatives of T with respect to the optical depth.
    do ilay = 0,nlay
      ! Perturbed layers above interface ilay, relative derivative is minus the path length.
      geo%reldiff_t_tot_0_tau(1:ilay,ilay) = -1.0 / geo%mu_0
      ! Perturbed layers below interface ilay. Transmission term insenstitve to those layers.
      geo%reldiff_t_tot_0_tau(ilay+1:nlay,ilay) = 0.0
      ! Transmission to interface ilay is sensitive to layers up (geometrically down) to ilay.
      geo%ilay_sensitive_0(ilay) = ilay

      ! Repeat joke for instrument, with a loop over geometries.
      do igeo = 1,ngeo
        geo%reldiff_t_tot_v_tau_allgeo(1:ilay,ilay,igeo) = -1.0 / geo%mu_v_allgeo(igeo)
        geo%reldiff_t_tot_v_tau_allgeo(ilay+1:nlay,ilay,igeo) = 0.0
        geo%ilay_sensitive_v_allgeo(ilay,igeo) = ilay
      enddo
    enddo

    ! All interfaces are visible if the sun is up. Otherwise, nothing is visible.
    if (geo%mu_0 .gt. 0.0) then
      geo%visible_0 = nlay
    else
      geo%visible_0 = 0
      geo%ilay_sensitive_0 = 0
    endif
    ! Repeat joke for instruments.
    do igeo = 1,ngeo
      if (geo%mu_v_allgeo(igeo) .gt. 0.0) then
        geo%visible_v_allgeo(igeo) = nlay
      else
        geo%visible_v_allgeo(igeo) = 0
        geo%ilay_sensitive_v_allgeo(:,igeo) = 0
      endif
    enddo

    ! Aggregate sun and instruments.
    do igeo = 1,ngeo
      geo%reldiff_t_tot_0v_tau_allgeo(:,:,igeo) = geo%reldiff_t_tot_0_tau + geo%reldiff_t_tot_v_tau_allgeo(:,:,igeo)
      geo%visible_0v_allgeo(igeo) = min(geo%visible_0,geo%visible_v_allgeo(igeo))
    enddo

    ! Just for optimization, to indicate that the world is a lot simpler
    ! than it could have been.
    geo%plane_parallel = .true.

    ! Set effective direction cosines with the geometry involved. They do not
    ! depend on the optical depth in plane-parallel geometry.
    geo%effmu_p0_ppg = (grd%mu * geo%mu_0) / (grd%mu + geo%mu_0)
    geo%effmu_m0_ppg = (grd%mu * geo%mu_0) / (grd%mu - geo%mu_0)
    do igeo = 1,ngeo
      geo%effmu_0v_ppg_allgeo(igeo) = (geo%mu_0 * geo%mu_v_allgeo(igeo)) / (geo%mu_0 + geo%mu_v_allgeo(igeo))
      geo%effmu_pv_ppg_allgeo(:,igeo) = (grd%mu * geo%mu_v_allgeo(igeo)) / (grd%mu + geo%mu_v_allgeo(igeo))
      geo%effmu_mv_ppg_allgeo(:,igeo) = (grd%mu * geo%mu_v_allgeo(igeo)) / (grd%mu - geo%mu_v_allgeo(igeo))
    enddo

  end subroutine geometry_set_plane_parallel ! }}}

  !> Applies pseudo-spherical geometry.
  !!
  !! This routine calculates the path length factors. Because this is pseudo-spherical geometry,
  !! these path length factors are quite complicated. With a minus signs, those are the relative
  !! derivatives of the cumulative transmission terms. From the relative derivatives of the
  !! cumulative transmission terms, all other transmission terms and related intregrals,
  !! in which the solar or viewing geometry is involved, will be calculated when the optical
  !! properties of the atmosphere are known.
  subroutine geometry_set_pseudo_spherical(nlay,ngeo,radius,rlev,geo,stat) ! {{{

    use lintran_constants_module, only: warningflag_displaecedsurface
    use lintran_types_module, only: lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    real, intent(in) :: radius !< Radius of the planet.
    real, dimension(0:nlay), intent(in) :: rlev !< Distance of interface levels from the center of the planeet.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure.
    integer, intent(inout) :: stat !< Error code.

    ! When the beam is directed to an interface, how close does it reach
    ! the planet center, including the part beyond the interface.
    real :: tangent ! Used for sun and instrument.

    ! Iterators.
    integer :: ilay ! Over layers, actually interfaces.
    integer :: igeo ! Over viewing geometries.

    ! Finding the layer of the tangent.
    integer :: ilay_tangent ! Index of the layer in which the tangent point is located.
    integer, dimension(1) :: ind ! Search index.

    ! Log if the atmosphere does not start exactly at the planet surface.
    if (rlev(nlay) .ne. radius) then
      stat = or(stat,warningflag_displaecedsurface)
    endif

    ! Throughout this calculation, a tangent point is used frequently. That point is
    ! acquired by taking the infinite line of the beam and take the point that is
    ! closest to the planet center. We will get right-angled triangles with a right
    ! angle at this point, using the center of the planet as well, because the distance
    ! from any interface point to the center of the planet is known.

    ! Set zero path to top of atmosphere. It will probably be overwritten. Only
    ! in the case of a fully invisible atmosphere, this value may prevent unidentified
    ! exceptions.
    geo%reldiff_t_tot_0_tau(:,0) = 0.0
    geo%reldiff_t_tot_v_tau_allgeo(:,0,:) = 0.0

    ! Start with the sun. The instrument will follow. {{{

    geo%visible_0 = nlay ! All layers are visible if not overwritten.

    ! Case 1: Sun higher than your nose (mu_0 higher than 0). {{{

    if (geo%mu_0 .gt. 0.0) then
      do ilay = 0,nlay
        tangent = rlev(ilay)*sqrt(1.0 - geo%mu_0**2.0)
        ! Draw two triangles, with two common points: the center of the planet and the
        ! point in the line that is closest to the planet, the tangent point
        ! The third point is where the beam crosses the interfaces around the perturbed layer.
        ! The path length is the difference between the lengths of that third point and the
        ! tangent point. The angle at the tangent point is ninety degrees, so Pythagoras can
        ! be used. Meanwhile, a minus sign is applied to convert path length to relative
        ! derivative. Therefore, the longer distance (from tangent point to interface with
        ! low index) has the minus sign. Also, there is a division by the vertical path
        ! length, so that it is a path length ampfilier, because tau is perturbed, not
        ! the extinction coefficient.
        geo%reldiff_t_tot_0_tau(1:ilay,ilay) = (sqrt(rlev(1:ilay)**2.0-tangent**2.0) - sqrt(rlev(0:ilay-1)**2.0-tangent**2.0)) / (rlev(0:ilay-1) - rlev(1:ilay))
        geo%reldiff_t_tot_0_tau(ilay+1:nlay,ilay) = 0.0 ! Layers are not crossed.
        ! Set sensitivity limit.
        geo%ilay_sensitive_0(ilay) = ilay
      enddo
    endif

    ! }}}

    ! Case 2: Sun is at nose altitude (mu_0 is 0). {{{

    if (geo%mu_0 .eq. 0.0) then
      ! Same algorithm, but the root sqrt(rlev(ilay)**2.0-tangent**2.0) is a sqrt(0.0)
      ! which must be protected from becoming a little negative root.
      geo%reldiff_t_tot_0_tau(:,0) = 0.0 ! All layers are not crossed. Ranges, like 1:ilay-1 will
      ! be ignored, but a single element, like ilay, will cause an out-of-bounds exception.
      do ilay = 1,nlay ! Interface 0 is skipped, see one line above.
        tangent = rlev(ilay)*sqrt(1.0 - geo%mu_0**2.0)
        ! Everything except perturbed layer ilay is the same.
        geo%reldiff_t_tot_0_tau(1:ilay-1,ilay) = (sqrt(rlev(1:ilay-1)**2.0-tangent**2.0) - sqrt(rlev(0:ilay-2)**2.0-tangent**2.0)) / (rlev(0:ilay-2) - rlev(1:ilay-1))
        ! Perturbed layer ilay has her sqrt(0.0) removed.
        geo%reldiff_t_tot_0_tau(ilay,ilay) = -sqrt(rlev(ilay-1)**2.0-tangent**2.0) / (rlev(ilay-1) - rlev(ilay))
        geo%reldiff_t_tot_0_tau(ilay+1:nlay,ilay) = 0.0 ! Layers are not crossed.
        ! Set sensitivity limit.
        geo%ilay_sensitive_0(ilay) = ilay
      enddo
    endif

    ! }}}

    ! Case 3: Sun is below your nose (mu_0 is lower than 0). {{{

    if (geo%mu_0 .lt. 0.0) then
      ! Different algorithm, because some layers are crossed more than once.
      ! The layer in which the tangent point is located will be entered at the
      ! top interface and also be left at the top interface, after which all the
      ! layers up until the interface are crossed again.
      do ilay = 0,nlay ! Now, interface 0 can be reached from below.
        tangent = rlev(ilay)*sqrt(1.0 - geo%mu_0**2.0)
        if (tangent .le. radius) then
          ! Interface is invisible. All further layers will be invisible as well,
          ! because they are even closer to the planet. Last fully visible layer is
          ! one lower than this index. If, for some strange reason, this already happens
          ! at interface 0, the entire model atmosphere is invisible. If only a part
          ! of a layer is visible, it is counted as invisible. This may raise an error
          ! but that is due to the transition between spherical and plane-parallel
          ! atmospheres.
          geo%visible_0 = max(0,ilay-1) ! We do not want a -1 here.
          ! Set sensitivity limit for invisible layers.
          geo%ilay_sensitive_0(ilay:nlay) = 0
          exit ! The rest is also invisible.
        else
          ! Find layer in which tangent point is located. The interface with this
          ! index is the lower interface of that layer, so that is the highest interface
          ! under the tangent point.
          ind = maxloc(rlev,mask=rlev .lt. tangent)
          ilay_tangent = ind(1) - 1 ! The minus one is because of the indexing.
          ! Protect and warn if tangent point is not in the model atmosphere. Do as if
          ! it is in the lowest layer.
          if (ilay_tangent .eq. 0) ilay_tangent = nlay

          ! Progress from the top of the atmosphere to the top of the tangent layer the
          ! same way as in case 1.
          geo%reldiff_t_tot_0_tau(1:ilay_tangent-1,ilay) = (sqrt(rlev(1:ilay_tangent-1)**2.0-tangent**2.0) - sqrt(rlev(0:ilay_tangent-2)**2.0-tangent**2.0)) / (rlev(0:ilay_tangent-2) - rlev(1:ilay_tangent-1))
          ! Tangent layer. Here, we draw a triangle with two equal sides by adding up
          ! two mirrored right-angled triangles to each other. The mirror line is from
          ! the center of the planet to the tangent point.
          geo%reldiff_t_tot_0_tau(ilay_tangent,ilay) = -2.0 * sqrt(rlev(ilay_tangent-1)**2.0-tangent**2.0) / (rlev(ilay_tangent-1) - rlev(ilay_tangent))
          ! The layers from the tangent layer to the interface are passed again, with
          ! the same path length as the path before the tangent layer.
          geo%reldiff_t_tot_0_tau(ilay+1:ilay_tangent-1,ilay) = 2.0 * geo%reldiff_t_tot_0_tau(ilay+1:ilay_tangent-1,ilay)
          ! Layers below the tangent layer are not crossed.
          geo%reldiff_t_tot_0_tau(ilay_tangent+1:nlay,ilay) = 0.0
          ! Set sensitivity limit.
          geo%ilay_sensitive_0(ilay) = ilay_tangent
        endif

      enddo

    endif

    ! }}}

    ! }}}

    ! Repeat joke for the instruments. {{{

    do igeo = 1,ngeo

      geo%visible_v_allgeo(igeo) = nlay ! All layers are visible if not overwritten.

      ! Case 1: Instrument higher than your nose (mu_v higher than 0). {{{

      if (geo%mu_v_allgeo(igeo) .gt. 0.0) then
        do ilay = 0,nlay
          tangent = rlev(ilay)*sqrt(1.0 - geo%mu_v_allgeo(igeo)**2.0)
          ! Draw two triangles, with two common points: the center of the planet and the
          ! point in the line that is closest to the planet, the tangent point
          ! The third point is where the beam crosses the interfaces around the perturbed layer.
          ! The path length is the difference between the lengths of that third point and the
          ! tangent point. The angle at the tangent point is ninety degrees, so Pythagoras can
          ! be used. Meanwhile, a minus sign is applied to convert path length to relative
          ! derivative. Therefore, the longer distance (from tangent point to interface with
          ! low index) has the minus sign. Also, there is a division by the vertical path
          ! length, so that it is a path length ampfilier, because tau is perturbed, not
          ! the extinction coefficient.
          geo%reldiff_t_tot_v_tau_allgeo(1:ilay,ilay,igeo) = (sqrt(rlev(1:ilay)**2.0-tangent**2.0) - sqrt(rlev(0:ilay-1)**2.0-tangent**2.0)) / (rlev(0:ilay-1) - rlev(1:ilay))
          geo%reldiff_t_tot_v_tau_allgeo(ilay+1:nlay,ilay,igeo) = 0.0 ! Layers are not crossed.
          ! Set sensitivity limit.
          geo%ilay_sensitive_v_allgeo(ilay,igeo) = ilay
        enddo
      endif

      ! }}}

      ! Case 2: Instrument is at nose altitude (mu_v is 0). {{{

      ! This case cannot apply, because mu_v is asserted not to be too close to zero. However,
      ! if that code is changed, this part of the code should be applied.
      if (geo%mu_v_allgeo(igeo) .eq. 0.0) then
        ! Same algorithm, but the root sqrt(rlev(ilay)**2.0-tangent**2.0) is a sqrt(0.0)
        ! which must be protected from becoming a little negative root.
        geo%reldiff_t_tot_v_tau_allgeo(:,0,igeo) = 0.0 ! All layers are not crossed. Ranges, like 1:ilay-1 will
        ! be ignored, but a single element, like ilay, will cause an out-of-bounds exception.
        do ilay = 1,nlay ! Interface 0 is skipped, see one line above.
          tangent = rlev(ilay)*sqrt(1.0 - geo%mu_v_allgeo(igeo)**2.0)
          ! Everything except perturbed layer ilay is the same.
          geo%reldiff_t_tot_v_tau_allgeo(1:ilay-1,ilay,igeo) = (sqrt(rlev(1:ilay-1)**2.0-tangent**2.0) - sqrt(rlev(0:ilay-2)**2.0-tangent**2.0)) / (rlev(0:ilay-2) - rlev(1:ilay-1))
          ! Perturbed layer ilay has her sqrt(0.0) removed.
          geo%reldiff_t_tot_v_tau_allgeo(ilay,ilay,igeo) = -sqrt(rlev(ilay-1)**2.0-tangent**2.0) / (rlev(ilay-1) - rlev(ilay))
          geo%reldiff_t_tot_v_tau_allgeo(ilay+1:nlay,ilay,igeo) = 0.0 ! Layers are not crossed.
          ! Set sensitivity limit.
          geo%ilay_sensitive_v_allgeo(ilay,igeo) = ilay
        enddo
      endif

      ! }}}

      ! Case 3: Instrument is below your nose (mu_v is lower than 0). {{{

      if (geo%mu_v_allgeo(igeo) .lt. 0.0) then
        ! Different algorithm, because some layers are crossed more than once.
        ! The layer in which the tangent point is located will be entered at the
        ! top interface and also be left at the top interface, after which all the
        ! layers up until the interface are crossed again.
        do ilay = 0,nlay ! Now, interface 0 can be reached from below.
          tangent = rlev(ilay)*sqrt(1.0 - geo%mu_v_allgeo(igeo)**2.0)
          if (tangent .lt. radius) then
            ! Interface is invisible. All further layers will be invisible as well,
            ! because they are even closer to the planet. Last fully visible layer is
            ! one lower than this index. If, for some strange reason, this already happens
            ! at interface 0, the entire model atmosphere is invisible. If only a part
            ! of a layer is visible, it is counted as invisible. This may raise an error
            ! but that is due to the transition between spherical and plane-parallel
            ! atmospheres.
            geo%visible_v_allgeo(igeo) = max(0,ilay-1) ! We do not want a -1 here.
            ! Set sensitivity limit for invisible layers.
            geo%ilay_sensitive_v_allgeo(ilay:nlay,igeo) = 0
            exit ! The rest is also invisible.
          else
            ! Find layer in which tangent point is located. The interface with this
            ! index is the lower interface of that layer, so that is the highest interface
            ! under the tangent point.
            ind = maxloc(rlev,mask=rlev .le. tangent)
            ilay_tangent = ind(1) - 1 ! The minus one is because of the indexing.
            ! Protect and warn if tangent point is not in the model atmosphere. Do as if
            ! it is in the lowest layer.
            if (ilay_tangent .eq. 0) ilay_tangent = nlay

            ! Progress from the top of the atmosphere to the top of the tangent layer the
            ! same way as in case 1.
            geo%reldiff_t_tot_v_tau_allgeo(1:ilay_tangent-1,ilay,igeo) = (sqrt(rlev(1:ilay_tangent-1)**2.0-tangent**2.0) - sqrt(rlev(0:ilay_tangent-2)**2.0-tangent**2.0)) / (rlev(0:ilay_tangent-2) - rlev(1:ilay_tangent-1))
            ! Tangent layer. Here, we draw a triangle with two equal sides by adding up
            ! two mirrored right-angled triangles to each other. The mirror line is from
            ! the center of the planet to the tangent point.
            geo%reldiff_t_tot_v_tau_allgeo(ilay_tangent,ilay,igeo) = -2.0 * sqrt(rlev(ilay_tangent-1)**2.0-tangent**2.0) / (rlev(ilay_tangent-1) - rlev(ilay_tangent))
            ! The layers from the tangent layer to the interface are passed again, with
            ! the same path length as the path before the tangent layer.
            geo%reldiff_t_tot_v_tau_allgeo(ilay+1:ilay_tangent-1,ilay,igeo) = 2.0 * geo%reldiff_t_tot_v_tau_allgeo(ilay+1:ilay_tangent-1,ilay,igeo)
            ! Layers below the tangent layer are not crossed.
            geo%reldiff_t_tot_v_tau_allgeo(ilay_tangent+1:nlay,ilay,igeo) = 0.0
            ! Set sensitivity limit.
            geo%ilay_sensitive_v_allgeo(ilay,igeo) = ilay_tangent
          endif

        enddo

      endif

      ! }}}

    enddo

    ! }}}

    ! Aggregate sun and instruments.
    do igeo = 1,ngeo
      geo%reldiff_t_tot_0v_tau_allgeo(:,:,igeo) = geo%reldiff_t_tot_0_tau + geo%reldiff_t_tot_v_tau_allgeo(:,:,igeo)
      geo%visible_0v_allgeo(igeo) = min(geo%visible_0,geo%visible_v_allgeo(igeo))
    enddo

    ! Just for optimization, indicating that the world is indeed
    ! as complicated as it could be.
    geo%plane_parallel = .false.

  end subroutine geometry_set_pseudo_spherical ! }}}

  !> Calculates the generalized spherical harmonics for the solar and viewing directions.
  !!
  !! In this function, generalized spherical harmonics for solar and viewing angles are calculated.
  !! These are required to calculate any phase matrix element in which one of these directions is
  !! involved, as well as their derivatives. Generalized spherical harmonics for internal streams
  !! are already handled in the grid module.
  subroutine calculate_geometry_generalized_spherical_functions(nst,nleg,ngeo,m,geo) ! {{{

    use lintran_constants_module, only: pol_nst
    use lintran_types_module, only: lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    integer, intent(in) :: m !< Fourier number.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure.

    integer :: ileg ! Legendre number.
    integer :: ibin ! Inplied do-loop iterator for binomial coefficient.
    real :: prefactor ! Common prefactor for initial Legendre polynomial.

    ! Building blocks of the p and m polynomials, who have their own recursive function.
    real :: temppolyp2_0
    real, dimension(ngeo) :: temppolyp2_v
    real :: temppolym2_0
    real, dimension(ngeo) :: temppolym2_v

    ! Recursive definition of the p and m building blocks, save history of the last
    ! two ones. The 'm' after the underscore stands for minus in the sense of 'l'-values,
    ! so that is the previous polynomial. The 'm' (or 'p') before the underscore stands
    ! for the type of function.
    real :: temppolyp2_m1_0
    real, dimension(ngeo) :: temppolyp2_m1_v
    real :: temppolyp2_m2_0
    real, dimension(ngeo) :: temppolyp2_m2_v
    real :: temppolym2_m1_0
    real, dimension(ngeo) :: temppolym2_m1_v
    real :: temppolym2_m2_0
    real, dimension(ngeo) :: temppolym2_m2_v

    ! Calculate the common prefactor.
    ! The prefactor is 2^(-m) * sqrt(2m over m), where the former is a
    ! power and the latter is the root of a binomial coefficient.
    prefactor = 2.0**(-m)*sqrt(product((/(float(m+ibin)/float(ibin),ibin=1,m)/)))

    ! Set everything at zero for l < m.
    geo%gsf_0_0(0:m-1,m) = 0.0
    geo%gsf_0_v_allgeo(0:m-1,m,:) = 0.0
    if (pol_nst(nst)) then ! Only for vector.
      ! Write a zero for Legendre coefficients 0 and 1, because that is always the
      ! case for gsf_p and gsf_m. In the module, Legendre numbers 0 and 1 are not
      ! calculated, so they are left uninitialized.
      geo%gsf_p_0(0:1,m) = 0.0
      geo%gsf_p_v_allgeo(0:1,m,:) = 0.0
      geo%gsf_m_0(0:1,m) = 0.0
      geo%gsf_m_v_allgeo(0:1,m,:) = 0.0
      ! If not already set to zero, write a zero for l < m.
      geo%gsf_p_0(2:m-1,m) = 0.0
      geo%gsf_p_v_allgeo(2:m-1,m,:) = 0.0
      geo%gsf_m_0(2:m-1,m) = 0.0
      geo%gsf_m_v_allgeo(2:m-1,m,:) = 0.0
    endif

    ! Set up the polynomials so that recursion can do the rest.

    ! Polynomial 0, set up ileg=m and ileg=m+1.
    geo%gsf_0_0(m,m) = prefactor * sqrt(1.0-geo%mu_0**2.0)**m
    geo%gsf_0_v_allgeo(m,m,:) = prefactor * sqrt(1.0-geo%mu_v_allgeo**2.0)**m
    if (m .ne. nleg) then
      geo%gsf_0_0(m+1,m) = sqrt(float(2*m+1))*geo%mu_0*geo%gsf_0_0(m,m)
      geo%gsf_0_v_allgeo(m+1,m,:) = sqrt(float(2*m+1))*geo%mu_v_allgeo*geo%gsf_0_v_allgeo(m,m,:)
    endif

    ! Full recursion for higher Legendre numbers.
    do ileg = m+2,nleg
      geo%gsf_0_0(ileg,m) = (float(2*ileg-1)*geo%mu_0*geo%gsf_0_0(ileg-1,m) - sqrt(float((ileg-1)**2-m**2))*geo%gsf_0_0(ileg-2,m))/sqrt(float(ileg**2-m**2))
      geo%gsf_0_v_allgeo(ileg,m,:) = (float(2*ileg-1)*geo%mu_v_allgeo*geo%gsf_0_v_allgeo(ileg-1,m,:) - sqrt(float((ileg-1)**2-m**2))*geo%gsf_0_v_allgeo(ileg-2,m,:))/sqrt(float(ileg**2-m**2))
    enddo

    if (pol_nst(nst)) then ! Only for vector.

      ! Polynomial p2 and m2, the lowest possible ileg is 2, or higher if m is higher
      ! than 2.

      ! This lowest ileg will be put in the _m2 polynomial and one higher in _m1.

      select case (m)
        case (0)
          temppolyp2_m2_0 = -0.25 * sqrt(6.0) * (1.0-geo%mu_0**2.0)
          temppolyp2_m2_v = -0.25 * sqrt(6.0) * (1.0-geo%mu_v_allgeo**2.0)
          temppolym2_m2_0 = -0.25 * sqrt(6.0) * (1.0-geo%mu_0**2.0)
          temppolym2_m2_v = -0.25 * sqrt(6.0) * (1.0-geo%mu_v_allgeo**2.0)
        case (1)
          temppolyp2_m2_0 = 0.5 * sqrt(1.0-geo%mu_0**2.0) * (1.0+geo%mu_0)
          temppolyp2_m2_v = 0.5 * sqrt(1.0-geo%mu_v_allgeo**2.0) * (1.0+geo%mu_v_allgeo)
          temppolym2_m2_0 = -0.5 * sqrt(1.0-geo%mu_0**2.0) * (1.0-geo%mu_0)
          temppolym2_m2_v = -0.5 * sqrt(1.0-geo%mu_v_allgeo**2.0) * (1.0-geo%mu_v_allgeo)
        case (2)
          temppolyp2_m2_0 = -0.25 * (1.0+geo%mu_0)**2.0
          temppolyp2_m2_v = -0.25 * (1.0+geo%mu_v_allgeo)**2.0
          temppolym2_m2_0 = -0.25 * (1.0-geo%mu_0)**2.0
          temppolym2_m2_v = -0.25 * (1.0-geo%mu_v_allgeo)**2.0
        case default
          temppolyp2_m2_0 = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-geo%mu_0**2.0)**(m-2) * (1.0+geo%mu_0)**2.0
          temppolyp2_m2_v = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-geo%mu_v_allgeo**2.0)**(m-2) * (1.0+geo%mu_v_allgeo)**2.0
          temppolym2_m2_0 = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-geo%mu_0**2.0)**(m-2) * (1.0-geo%mu_0)**2.0
          temppolym2_m2_v = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-geo%mu_v_allgeo**2.0)**(m-2) * (1.0-geo%mu_v_allgeo)**2.0
      endselect

      ! The ileg is set to two, except for the `case default', where ileg is set to m.
      ileg = max(m,2)

      ! Write results.
      geo%gsf_p_0(ileg,m) = 0.5 * (temppolym2_m2_0 + temppolyp2_m2_0)
      geo%gsf_p_v_allgeo(ileg,m,:) = 0.5 * (temppolym2_m2_v + temppolyp2_m2_v)
      geo%gsf_m_0(ileg,m) = 0.5 * (temppolym2_m2_0 - temppolyp2_m2_0)
      geo%gsf_m_v_allgeo(ileg,m,:) = 0.5 * (temppolym2_m2_v - temppolyp2_m2_v)

      if (m .ne. nleg) then

        ! Acquire next polynomials, for ileg+1. Polynomial _m2 is now the previous one
        ! and the next-to-previous one does not exist yet.
        temppolyp2_m1_0 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_0-2.0*m) * temppolyp2_m2_0) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolyp2_m1_v = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_v_allgeo-2.0*m) * temppolyp2_m2_v) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2_m1_0 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_0+2.0*m) * temppolym2_m2_0) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2_m1_v = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_v_allgeo+2.0*m) * temppolym2_m2_v) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))

        ! Write results for the next Legendre number.
        geo%gsf_p_0(ileg+1,m) = 0.5 * (temppolym2_m1_0 + temppolyp2_m1_0)
        geo%gsf_p_v_allgeo(ileg+1,m,:) = 0.5 * (temppolym2_m1_v + temppolyp2_m1_v)
        geo%gsf_m_0(ileg+1,m) = 0.5 * (temppolym2_m1_0 - temppolyp2_m1_0)
        geo%gsf_m_v_allgeo(ileg+1,m,:) = 0.5 * (temppolym2_m1_v - temppolyp2_m1_v)

      endif

      ! And process all ilegs that are higher. Note that the iterator ileg is
      ! one lower than the actual Legendre number, that is how the recursion is
      ! defined and the formulas are complicated.
      do ileg = max(m,2)+1,nleg-1

        temppolyp2_0 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_0-2.0*m) * temppolyp2_m1_0 - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolyp2_m2_0) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolyp2_v = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_v_allgeo-2.0*m) * temppolyp2_m1_v - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolyp2_m2_v) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2_0 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_0+2.0*m) * temppolym2_m1_0 - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolym2_m2_0) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2_v = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*geo%mu_v_allgeo+2.0*m) * temppolym2_m1_v - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolym2_m2_v) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))

        ! Write results. The ileg is still one lower than what it should be.
        geo%gsf_p_0(ileg+1,m) = 0.5 * (temppolym2_0 + temppolyp2_0)
        geo%gsf_p_v_allgeo(ileg+1,m,:) = 0.5 * (temppolym2_v + temppolyp2_v)
        geo%gsf_m_0(ileg+1,m) = 0.5 * (temppolym2_0 - temppolyp2_0)
        geo%gsf_m_v_allgeo(ileg+1,m,:) = 0.5 * (temppolym2_v - temppolyp2_v)

        ! Progress previous results for next recursive iteration.
        temppolyp2_m2_0 = temppolyp2_m1_0
        temppolyp2_m2_v = temppolyp2_m1_v
        temppolym2_m2_0 = temppolym2_m1_0
        temppolym2_m2_v = temppolym2_m1_v
        temppolyp2_m1_0 = temppolyp2_0
        temppolyp2_m1_v = temppolyp2_v
        temppolym2_m1_0 = temppolym2_0
        temppolym2_m1_v = temppolym2_v

      enddo

    endif

  end subroutine calculate_geometry_generalized_spherical_functions ! }}}

  !> Cleans up the geometry.
  !!
  !! All arrays that are allocated in geometry_init are cleaned up.
  subroutine geometry_close(nst,geo,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, pol_nst
    use lintran_types_module, only: lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    deallocate(geo%diff_c_emi_v_elem_emi_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%diff_c_bck_srf_0v_elem_bdrf_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%diff_c_bck_srf_v_bdrf_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%diff_c_bck_srf_0_bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then
      deallocate(geo%gsf_m_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%gsf_m_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%gsf_p_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%gsf_p_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(geo%gsf_0_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%gsf_0_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then ! Only for vector.
      deallocate(geo%rotmat_s_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%rotmat_s_0_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%rotmat_c_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(geo%rotmat_c_0_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(geo%reldiff_t_tot_0v_tau_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%reldiff_t_tot_v_tau_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%reldiff_t_tot_0_tau,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%effmu_mv_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%effmu_pv_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%effmu_m0_ppg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%effmu_p0_ppg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%effmu_0v_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%degenerate_ssg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%degenerate_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%ilay_sensitive_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%ilay_sensitive_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%visible_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%visible_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%azi_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(geo%mu_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine geometry_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(nst,progression,geo,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, pol_nst
    use lintran_types_module, only: lintran_geometry_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    type(lintran_geometry_class), intent(inout) :: geo !< Internal geometry structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (progression .ge. 30) deallocate(geo%diff_c_emi_v_elem_emi_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 29) deallocate(geo%diff_c_bck_srf_0v_elem_bdrf_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 28) deallocate(geo%diff_c_bck_srf_v_bdrf_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 27) deallocate(geo%diff_c_bck_srf_0_bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then
      if (progression .ge. 26) deallocate(geo%gsf_m_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 25) deallocate(geo%gsf_m_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 24) deallocate(geo%gsf_p_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 23) deallocate(geo%gsf_p_0,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 22) deallocate(geo%gsf_0_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 21) deallocate(geo%gsf_0_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then ! Only for vector.
      if (progression .ge. 20) deallocate(geo%rotmat_s_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 19) deallocate(geo%rotmat_s_0_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 18) deallocate(geo%rotmat_c_v_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 17) deallocate(geo%rotmat_c_0_allgeo,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 16) deallocate(geo%reldiff_t_tot_0v_tau_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 15) deallocate(geo%reldiff_t_tot_v_tau_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 14) deallocate(geo%reldiff_t_tot_0_tau,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 13) deallocate(geo%effmu_mv_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(geo%effmu_pv_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 11) deallocate(geo%effmu_m0_ppg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 10) deallocate(geo%effmu_p0_ppg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 9) deallocate(geo%effmu_0v_ppg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(geo%degenerate_ssg_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(geo%degenerate_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(geo%ilay_sensitive_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(geo%ilay_sensitive_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(geo%visible_0v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(geo%visible_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(geo%azi_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(geo%mu_v_allgeo,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail ! }}}

end module lintran_geometry_module
