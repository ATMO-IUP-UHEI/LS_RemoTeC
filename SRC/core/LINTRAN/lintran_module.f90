!> \file lintran_module.f90
!! Module file for the main interface of Lintran.

!> Module for the main interface of Lintran.
!!
!! This module initializes, provides geometries, performs calculations and closes
!! off Lintran. All routines that should be called by the external program are defined
!! in this module, some in polymorph form. Not all routines defined here should be called
!! by the user. Some are called directly by the interface routines, but are so close to
!! the interface that is seems logical to have them in this module as well. Those are
!! lintran_execute_without_derivatives and lintran_execute_with_derivatives.
module lintran_module

  implicit none

  !> Provides geometry to Lintran.
  !!
  !! Altitude information is needed when choosing for pseudo-spherical geometry.
  interface lintran_provide
    procedure lintran_provide_without_altitudes
    procedure lintran_provide_with_altitudes
  end interface lintran_provide

  !> Main calculation of Lintran.
  !!
  !! Whenever the structure for derivatives is given, they will be calculated. But
  !! then, Lintran must be initialized for supporting derivatives.
  interface lintran_calculate
    procedure lintran_calculate_without_derivatives
    procedure lintran_calculate_with_derivatives
  end interface lintran_calculate

  !> Nakajima wrapper.
  !!
  !! Calculates the multi-scattering contribution, using the original Nakajima method.
  !! This does not require a single-scattering geometry phase function in the atmosphere,
  !! but it requires the BDRF and fluorescent emissivity at single-scattering geometry.
  !! Nakajima multiscattering is equal to the total-scattering field, with truncated
  !! Legendre series for the single-scttering phase function, using &delta;-M minus the
  !! single-scattering field without &delta;-M with the same runcated phase function.
  !! Adding external precise single-scattering without &delta;-M will approximate the
  !! total radiance. However, using Lintran to calculate single-scattering with
  !! &delta;-M as well will give more precise results, but the administration and chain
  !! rules of all the derivatives for single-scattering will take relatively a lot of time.
  !! Nakajima multiscattering is useful when using an optimization scheme for multi-scattering.
  interface lintran_nakajima_multiscattering
    procedure lintran_nakajima_multiscattering_without_derivatives
    procedure lintran_nakajima_multiscattering_with_derivatives
  end interface lintran_nakajima_multiscattering

contains

  !> Initializes Lintran.
  !!
  !! This routine allocates all memory needed by Lintran. It needs all relevant dimension
  !! sizes, the number of streams, atmospheric layers and viewing geometries. During the
  !! initialization, all fixed information is already calculated. The user can choose whether
  !! or not derivatives can be calculated. If derivatives are turned on, the user is not forced
  !! to include derivatives for all runs, but if it is turned off, no derivatives can be
  !! calculated. When turning on derivatives, additional memory is calculated.
  subroutine lintran_init(nst,streams,nlay,ngeo,flag_derivatives,self,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, errorflag_stokes, warningflag_oddstreams, nst_max, legal_nst, nscat_nst
    use lintran_types_module, only: lintran_class
    use lintran_grid_module, only: grid_init
    use lintran_geometry_module, only: geometry_init
    use lintran_pixel_module, only: pixel_init
    use lintran_tlc_module, only: tlc_init
    use lintran_fields_module, only: fields_init

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: streams !< Number of streams. This must be an even number, so there is an equal number of upward and downward streams.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_class), intent(out) :: self !< Constructed main Lintran instance.
    integer, intent(out) :: stat !< Error code.

    ! Derived dimension sizes on local scope.
    integer :: nscat ! Number of independent matrix elements in a scattering matrix.
    integer :: nstrhalf ! Half number of streams.
    integer :: nleg ! Highest Legendre coefficient.

    ! Check for legal number of Stokes parameters without possibility of out-of-bounds exception.
    logical :: legal

    ! Error handling.
    integer :: progression

    progression = 0

    ! Initialize error flag.
    stat = 0

    ! Check number of Stokes parameters.
    if (nst .ge. 1 .and. nst .le. nst_max) then ! This prevents an out-of-bounds exception.
      legal = legal_nst(nst)
    else
      legal = .false.
    endif

    if (.not. legal) then
      ! Raise warning, actually an error, but Lintran knows only warnings.
      stat = or(stat,errorflag_stokes)
      call fail(progression,flag_derivatives,self,stat)
      return
    endif

    ! Look up the number of independent matrix elements in the scattering matrix.
    nscat = nscat_nst(nst)

    ! Check number of streams.
    if (modulo(streams,2) .eq. 1) then
      stat = or(stat,warningflag_oddstreams)
      ! The even number (rounded down) will automatically be applied by integer division.
    endif

    ! Calculate derived dimension sizes.
    nstrhalf = streams / 2 ! Half number of streams (auto-correcting if input streams is odd).
    nleg = 2*nstrhalf-1 ! Highest Legendre number (using corrected number of streams).

    ! Save dimension sizes.
    self%nst = nst
    self%nscat = nscat
    self%nstrhalf = nstrhalf
    self%nleg = nleg
    self%nlay = nlay
    self%ngeo = ngeo

    ! Initialize all modules.

    ! The grid module actually does everything it will ever need to do. It
    ! calculates the properties that do not depend on any parameter except
    ! for the dimensions selected above. Furthermore, the grid module pre-
    ! calculates all generalized spherical functions for the internal streams,
    ! so that phase matrix elements between internal streams can be calculated from
    ! the coefficients. Depending on the number of Stokes parameters, complicated
    ! generalized spherical functions are either or not calculated.
    call grid_init(nst,nstrhalf,nleg,flag_derivatives,self%grd,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail(progression,flag_derivatives,self,stat)
      return
    endif ! }}}
    progression = progression + 1

    ! The geometry module is for everything that depends on the solar or
    ! viewing angle, including the azimuthal difference. These properties
    ! vary from measurement to measurement, but usually do not vary between
    ! different wavelengths. During the initialization, the arrays are
    ! only shaped, not calculated, because no geometry is available yet.
    call geometry_init(nst,nstrhalf,nleg,nlay,ngeo,self%geo,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail(progression,flag_derivatives,self,stat)
      return
    endif ! }}}
    progression = progression + 1

    ! The pixel module is for everything that changes each calculation.
    ! Like for the geometry module, the arrays are only shaped, because
    ! nothing can be calculated yet.
    call pixel_init(nst,nscat,nstrhalf,nleg,nlay,ngeo,flag_derivatives,self%pix,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail(progression,flag_derivatives,self,stat)
      return
    endif ! }}}
    progression = progression + 1

    ! Shape arrays for terms appearing in the radiative transfer equation.
    ! Terms that do depend on the pixel module are allocated and shaped.
    ! Terms that only depend on members of the grid or the geometry are
    ! pointers. Then, they need not be looked after and will always be
    ! consistent with the current geometry. The pixel, geometry and
    ! grid modules are not accessible for the radiative transfer modules.
    ! Everything is done via the TLC module.
    call tlc_init(nst,nstrhalf,nlay,ngeo,flag_derivatives,self%grd,self%geo,self%pix,self%tlc,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail(progression,flag_derivatives,self,stat)
      return
    endif ! }}}
    progression = progression + 1

    ! Create intensity fields for analytic double-scattering. Those are the ones
    ! of which the size is known at this moment.
    call fields_init(nst,nstrhalf,nlay,flag_derivatives,self%fld,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail(progression,flag_derivatives,self,stat)
      return
    endif ! }}}
    progression = progression + 1

    ! Save derivatives flag.
    self%flag_derivatives = flag_derivatives

  end subroutine lintran_init ! }}}

  !> Provides a geometry to Lintran without giving altitude infromation.
  !!
  !! One solar direction and a number of viewing directions are provided to construct the
  !! geometry. The number of viewing geometries must be equal to the number of viewing
  !! geometries for which Lintran is initialized. In this form of lintran_provide, the
  !! optional arguments for the vertical dimensions are not given, so pseudo-spherical
  !! geometry cannot be applied. The user will be warned if he or she attempts that.
  subroutine lintran_provide_without_altitudes(mu_0,mu_v,azi,flag_pseudo_spherical,self,stat) ! {{{

    use lintran_constants_module, only: warningflag_pseudosphericalerror
    use lintran_types_module, only: lintran_class
    use lintran_geometry_module, only: geometry_set_general, geometry_set_plane_parallel

    implicit none

    ! Input and output.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance, defined as first argument, because it owns an important dimension size.
    real, intent(in) :: mu_0 !< Cosine of solar zenith angle.
    real, dimension(self%ngeo), intent(in) :: mu_v !< Cosines of viewing zenith angles.
    real, dimension(self%ngeo), intent(in) :: azi !< Azimuthal differences in radians, defined as in De Haan et al. (1987).
    logical, intent(in) :: flag_pseudo_spherical !< Flag for attempting pseudo-spherical geometry, though it is not possible in this form of lintran_provide.
    integer, intent(out) :: stat !< Error code.

    ! Initialize error flag.
    stat = 0

    if (flag_pseudo_spherical) then
      stat = or(stat,warningflag_pseudosphericalerror)
      ! This setting will be ignored in this implementation without altitudes.
    endif

    ! Set the geometry to the one of this pixel and calculate all
    ! stuff that can be calculated with just the geometry.
    call geometry_set_general(self%nst,self%nleg,self%ngeo,mu_0,mu_v,azi,self%grd,self%geo,stat)
    call geometry_set_plane_parallel(self%nstrhalf,self%nlay,self%ngeo,self%grd,self%geo,stat)

  end subroutine lintran_provide_without_altitudes ! }}}

  !> Provides a geometry to Lintran with giving altitude infromation.
  !!
  !! One solar direction and a number of viewing directions are provided to construct the
  !! geometry. The number of viewing geometries must be equal to the number of viewing
  !! geometries for which Lintran is initialized. In this form of lintran_provide, the
  !! optional arguments for the vertical dimensions are given, so pseudo-spherical
  !! geometry is possible. However, the user can still choose for plane-parallel geometry.
  !! In such a case, the optional altitude information is ignored.
  subroutine lintran_provide_with_altitudes(mu_0,mu_v,azi,flag_pseudo_spherical,radius,zlev,self,stat) ! {{{

    use lintran_types_module, only: lintran_class
    use lintran_geometry_module, only: geometry_set_general, geometry_set_plane_parallel, geometry_set_pseudo_spherical

    implicit none

    ! Input and output.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance, defined as first argument, because it owns some important dimension sizes.
    real, intent(in) :: mu_0 !< Cosine of solar zenith angle.
    real, dimension(self%ngeo), intent(in) :: mu_v !< Cosines of viewing zenith angles.
    real, dimension(self%ngeo), intent(in) :: azi !< Azimuthal differences in radians, defined as in De Haan et al. (1987).
    logical, intent(in) :: flag_pseudo_spherical !< Flag for pseudo-spherical geometry.
    real, intent(in) :: radius !< Radius of the planet in any unit.
    real, dimension(0:self%nlay), intent(in) :: zlev !< Altitudes of interfaces above surface in the same unit as radius. Index 0 is the top boundary of the top layer of the atmosphere. Index nlay is the bottom of the lowest layer, which should be at zero altitude.
    integer, intent(out) :: stat !< Error code.

    ! Initialize error flag.
    stat = 0

    ! Set the geometry to the one of this pixel and calculate all
    ! stuff that can be calculated with just the geometry
    ! In the case of pseudo-spherical geometry, the altitude grid is included,
    ! and is here converted to altitudes with respect to the center of the Earth,
    ! hereafter called rlev instead of zlev. In the case that no pseudo-spherical
    ! geometry is chosen, the altitude grid is dropped.
    call geometry_set_general(self%nst,self%nleg,self%ngeo,mu_0,mu_v,azi,self%grd,self%geo,stat)
    if (flag_pseudo_spherical) then
      call geometry_set_pseudo_spherical(self%nlay,self%ngeo,radius,zlev+radius,self%geo,stat)
    else
      call geometry_set_plane_parallel(self%nstrhalf,self%nlay,self%ngeo,self%grd,self%geo,stat)
    endif
    ! The polymorphism is not longer there in geometry_set, because it is not inside
    ! the interface outside Lintran.

  end subroutine lintran_provide_with_altitudes ! }}}

  !> Performs calculation without derivatives.
  !!
  !! This routine runs Lintran without derivatives. It is a form of the polymorph calculation
  !! and will be called if the derivatives structure is not passed.
  subroutine lintran_calculate_without_derivatives(execution,atm,set,self,rint,stat) ! {{{

    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_class

    implicit none

    ! Input and output.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(self%nst,self%ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries for the selected scattering orders.
    integer, intent(out) :: stat !< Error code.

    ! Initialize error flag.
    stat = 0

    ! Meanwhile, execution(2) is turned off if double-scattering is not split off, because
    ! that would mean nothing. Turning off execution(2) makes optimization easier.
    call lintran_execute_without_derivatives(self%nst,self%nscat,self%nstrhalf,self%nleg,self%nlay,self%ngeo,execution .and. (/.true.,set%split_double,.true./),atm,set,self,rint,stat)

  end subroutine lintran_calculate_without_derivatives ! }}}

  !> Performs calculation with derivatives.
  !!
  !! This routine runs Lintran with derivatives. It is a form of the polymorph calculation
  !! and will be called if the derivatives structure is passed.
  subroutine lintran_calculate_with_derivatives(execution,atm,set,self,rint,drv,stat) ! {{{

    use lintran_constants_module, only: warningflag_unpreparedderivatives
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

    implicit none

    ! Input and output.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(self%nst,self%ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries for the selected scattering orders.
    type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.
    integer, intent(out) :: stat !< Error code.

    ! Dimension sizes to pass to internal routine.
    integer :: nleg_pert ! A dimension size that will be derived.
    integer :: nlay_deriv_base ! A dimension size that will be derived.
    integer :: nlay_deriv_ph ! A dimension size that will be derived.

    ! Iterator.
    integer :: ilay ! Over atmospheric layers.

    ! Flags for layer to be differentiated in any way, for differentiating C.
    logical, dimension(self%nlay) :: differentiated

    ! Initialize error flag.
    stat = 0

    if (.not. self%flag_derivatives) then
      call lintran_calculate_without_derivatives(execution,atm,set,self,rint,stat)
      ! This routine is an interface call where the status flag is intent out. Therefore,
      ! this specific warning will be added at the end. Otherwise, it will be lost.
      stat = or(stat,warningflag_unpreparedderivatives)
      return
    endif

    ! Dimension sizes needed for declaring an intermediate result in the execution.
    nleg_pert = min(self%nleg,set%nleg_deriv)
    ! These pack-statements allow extra zeros in the array and they will be ignored.
    ! But the zeros can only be put in the end. The values that make sense must
    ! be put in front of the array.
    nlay_deriv_base = size(pack(set%ilay_deriv_base,set%ilay_deriv_base .ne. 0))
    nlay_deriv_ph = size(pack(set%ilay_deriv_ph,set%ilay_deriv_ph .ne. 0))

    if (nlay_deriv_base .eq. 0 .and. nlay_deriv_ph .eq. 0) then
      differentiated = .false.
    endif
    if (nlay_deriv_base .eq. 0 .and. nlay_deriv_ph .ne. 0) then
      differentiated = (/(any(set%ilay_deriv_ph(1:nlay_deriv_ph) .eq. ilay),ilay = 1,self%nlay)/)
    endif
    if (nlay_deriv_base .ne. 0 .and. nlay_deriv_ph .eq. 0) then
      differentiated = (/(any(set%ilay_deriv_base(1:nlay_deriv_base) .eq. ilay),ilay = 1,self%nlay)/)
    endif
    if (nlay_deriv_base .ne. 0 .and. nlay_deriv_ph .ne. 0) then
      differentiated = (/(any(set%ilay_deriv_base(1:nlay_deriv_base) .eq. ilay) .or. any(set%ilay_deriv_ph(1:nlay_deriv_ph) .eq. ilay),ilay = 1,self%nlay)/)
    endif

    ! Meanwhile, execution(2) is turned off if double-scattering is not split off, because
    ! that would mean nothing. Turning off execution(2) makes optimization easier.
    call lintran_execute_with_derivatives(self%nst,self%nscat,self%nstrhalf,self%nleg,self%nlay,self%ngeo,nlay_deriv_base,nlay_deriv_ph,count(differentiated),pack((/(ilay,ilay=1,self%nlay)/),differentiated),nleg_pert,execution .and. (/.true.,set%split_double,.true./),atm,set,self,rint,drv,stat)

  end subroutine lintran_calculate_with_derivatives ! }}}

  !> Performs original Nakajima multiscattering without derivatives.
  !!
  !! This is a wrapper for the normal calculations. But now, multi-scattering is calculated
  !! using &delta;-M, but the contribution of single scattering is limited to real single
  !! scattering. In &delta;-M, scattering in the forward peak is not counted as a scattering
  !! order, giving lower scattering orders than what would be real. Using Lintran's &delta;-M
  !! is, in principle, better then this Nakajima wrapper, but then single-scattering must also be
  !! calculated with Lintran. Performing a single-scattering calculation on your own results in
  !! in wrong results, unless the user watches very carefully and then, it comes down to a
  !! single-scattering calculation with Lintran. This is the form without derivatives, which
  !! is called when the derivatives structure is not passed.
  subroutine lintran_nakajima_multiscattering_without_derivatives(atm_inp,set_inp,self,rint,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, errorflag_deallocation, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_class

    implicit none

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm_inp !< Atmosphere provided by the user. Suffix _inp is added, because the Nakajima wrapper performs the desired calculation by making correction in the structure. The input structure remains unharmed.
    type(lintran_settings), intent(in) :: set_inp !< Calculation settings provided by the user. Suffix _inp is added, because the Nakajima wrapper performs the desired calculation by making correction in the structure. The input structure remains unharmed.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(self%nst,self%ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries.
    integer, intent(out) :: stat !< Error code.

    ! Radiance for a new calculation.
    real, dimension(self%nst,self%ngeo) :: rint_new

    ! Atmosphere and settings with possible hacks for Nakajima.
    type(lintran_atmosphere) :: atm
    type(lintran_settings) :: set

    ! Copy of dimension sizes.
    integer :: nst ! Number of Stokes parameters.
    integer :: nscat ! Number of relevant independent matrix elements in scattering matrix.
    integer :: nstrhalf ! Half the number of streams.
    integer :: nleg ! Highest Legendre coefficient.
    integer :: nlay ! Number of atmospheric layers.
    integer :: ngeo ! Number of viewing geometries.
    integer :: nleg_lay_cur ! Current number of Legendre coefficients in a layer, limited by the total number of Legendre coefficients.

    ! Iterators.
    integer :: ileg ! Over Legendre polynomials.
    integer :: ilay ! Over atmospheric layers.
    integer :: igeo ! Over viewing geometries.

    ! Generalized spherical functions at single-scattering geometry.
    ! These are the two numbers as lower index of the generalized spherical function, the m denotes a minus sign for the number on the right sie.
    real, dimension(0:self%nleg,self%ngeo) :: gsf_00
    real, dimension(0:self%nleg,self%ngeo) :: gsf_02
    real, dimension(0:self%nleg,self%ngeo) :: gsf_22
    real, dimension(0:self%nleg,self%ngeo) :: gsf_2m2 

    real, dimension(0:self%nleg) :: definitionconversion ! Array of 2*ileg+1 numbers.

    ! Single-scattering geometry.
    real :: uscat_ssg

    ! Forward peak for artifical delta-M.
    real :: forward_peak

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Initialize error flag.
    stat = 0

    ! Get dimension sizes on local scope.
    nst = self%nst
    nscat = self%nscat
    nstrhalf = self%nstrhalf
    nleg = self%nleg
    nlay = self%nlay
    ngeo = self%ngeo

    ! Copy atmosphere.
    allocate(atm%taus(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%taua(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%coefs(nscat,0:nleg+1,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%phase_ssg(nscat,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf(nst,nst,nstrhalf,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_0(nst,nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_v(nst,nst,nstrhalf,0:nleg,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_ssg(nst,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%emi(nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%emi_ssg(nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (atm_inp%thermal_emission) then
      allocate(atm%planck_curve(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 1 ! Fake allocation.
    endif
    allocate(atm%nleg_lay(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Copy atmosphere except for the phase function at single-scattering geometry. Meanwhile,
    ! the atmosphere may be encapsulated to the part that is used. Lintran used to cut off the
    ! used part of the atmosphere only in the pixel module, but that is because the arrays in
    ! the pixel module are allocated by Lintran, not by the user. This atmosphere is also
    ! allocated by Lintran, so the sizes must match. The number of atmospheric layers is not
    ! cut off, because truncating that will result in nonsense.
    atm%sun = atm_inp%sun
    atm%taus = atm_inp%taus
    atm%taua = atm_inp%taua
    atm%coefs = atm_inp%coefs(1:nscat,0:nleg+1,:) ! Possibly cut off non-used matrix elements and high Legendre coefficients.
    ! Single-scattering phase function not yet defined.
    atm%bdrf = atm_inp%bdrf(1:nst,1:nst,:,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%bdrf_0 = atm_inp%bdrf_0(1:nst,1:nst,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%bdrf_v = atm_inp%bdrf_v(1:nst,1:nst,:,0:nleg,1:ngeo) ! Possibly cut off non-used Stokes parameters, high Fourier coefficients and non-used viewing geometries.
    atm%bdrf_ssg = atm_inp%bdrf_ssg(1:nst,1:nst,1:ngeo) ! Possibly cut off non-used Stokes parameters and non-used viewing geometries.
    atm%emi = atm_inp%emi(1:nst,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%emi_ssg = atm_inp%emi_ssg(1:nst,1:ngeo) ! Possibly cut off non-used Stokes parameters and non-used viewing geometries.
    if (atm_inp%thermal_emission) atm%planck_curve = atm_inp%planck_curve
    atm%nleg_lay = atm_inp%nleg_lay
    atm%bdrf_only_0 = atm_inp%bdrf_only_0
    atm%emi_only_0 = atm_inp%emi_only_0
    atm%thermal_emission = atm_inp%thermal_emission

    ! Copy settings.
    allocate(set%execute_stokes(nst),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    set%execute_stokes = set_inp%execute_stokes(1:nst)
    set%taus_split = set_inp%taus_split
    set%taua_split = set_inp%taua_split
    set%tautot_max = set_inp%tautot_max
    set%fourier_tolerance = set_inp%fourier_tolerance
    set%split_double = set_inp%split_double
    set%interpolation = set_inp%interpolation
    ! Delta-M and some stuff about derivatives are not copied.
    set%solver = set_inp%solver
    set%gs_maxiter = set_inp%gs_maxiter
    set%gs_tolerance = set_inp%gs_tolerance

    ! We desire to do array operations over Legendre coefficients, so the term
    ! (2*ileg + 1) must be put into an array. And we will directly transform it
    ! into a floating point.
    definitionconversion = (/(float(2*ileg+1),ileg=0,nleg)/)

    do igeo = 1,ngeo

      ! Calculate the single-scattering angle. This is the only reason to save the azimuth
      ! on the geometry structure.
      uscat_ssg = -self%geo%mu_0*self%geo%mu_v_allgeo(igeo) + sqrt((1.0-self%geo%mu_0**2.0)*(1.0-self%geo%mu_v_allgeo(igeo)**2.0)) * cos(self%geo%azi_allgeo(igeo))

      ! An artifical truncated single-scattering phase function is constructed, using
      ! the phase coefficients after delta-M correction. This involves a conversion from
      ! Legendre coefficients to a normal phase function, so this involves the generalized
      ! spherical functions for single-scattering geometry.

      gsf_00(0,igeo) = 1.0
      ! Use recursive definition of generalized spherical functions.
      do ileg = 1,nleg

        if (ileg .eq. 1) then
          ! Recursive construction, knowing that for l-2, the polynomials
          ! does not yet exist.
          gsf_00(ileg,igeo) = float(2*ileg-1)*uscat_ssg*gsf_00(ileg-1,igeo)/sqrt(float(ileg**2))
        else
          ! Full recursion.
          gsf_00(ileg,igeo) = (float(2*ileg-1)*uscat_ssg*gsf_00(ileg-1,igeo) - sqrt(float((ileg-1)**2))*gsf_00(ileg-2,igeo))/sqrt(float(ileg**2))
        endif

      enddo

      if (pol_nst(nst)) then ! Only for vector.

        gsf_02(0:1,igeo) = 0.0
        gsf_02(2,igeo) = -0.25 * sqrt(6.0) * (1.0 - uscat_ssg**2.0)
        gsf_22(0:1,igeo) = 0.0
        gsf_22(2,igeo) = 0.25 * (1.0 + uscat_ssg)**2.0
        gsf_2m2(0:1,igeo) = 0.0
        gsf_2m2(2,igeo) = 0.25 * (1.0 - uscat_ssg)**2.0

        do ileg = 2,nleg-1

          ! Here, m=0.
          gsf_02(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg) * gsf_02(ileg,igeo) - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0)) * gsf_02(ileg-1,igeo))
          ! Here m=2.
          gsf_22(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-4.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg-4.0) * gsf_22(ileg,igeo) - (ileg+1)*sqrt((ileg**2.0-4.0)*(ileg**2.0-4.0)) * gsf_22(ileg-1,igeo))
          gsf_2m2(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-4.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg+4.0) * gsf_2m2(ileg,igeo) - (ileg+1)*sqrt((ileg**2.0-4.0)*(ileg**2.0-4.0)) * gsf_2m2(ileg-1,igeo))

        enddo

      endif

      ! Nakajima multiscattering is:
      ! + total-scattering with delta-M and truncated phase function
      ! - single-scattering without delta-M, with truncated phase function

      ! Total scattering with delta-M and truncated phase function.
      ! Construct truncated phase function.

      ! All existing elements are written, though in the physical world, only a1 and b1
      ! matter, because the sun radiates I.
      do ilay = 1,nlay
        if (atm%nleg_lay(ilay) .gt. nleg) then
          forward_peak = atm%coefs(iscat_a1,nleg+1,ilay) / (2.0*nleg+3.0)
        else
          forward_peak = 0.0
        endif
        ! Limit number of Legendre coefficients by the layer or by the model.
        nleg_lay_cur = min(nleg,atm%nleg_lay(ilay))
        atm%phase_ssg(iscat_a1,ilay,igeo) = dot_product(gsf_00(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a1,0:nleg_lay_cur,ilay) - definitionconversion(0:nleg_lay_cur) * forward_peak)
        if (pol_nst(nst)) then
          atm%phase_ssg(iscat_b1,ilay,igeo) = dot_product(gsf_02(0:nleg_lay_cur,igeo) , atm%coefs(iscat_b1,0:nleg_lay_cur,ilay)) ! Off-diagonal, so no subtraction of the forward peak.
          atm%phase_ssg(iscat_a2,ilay,igeo) = 0.5*(dot_product(gsf_22(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) + atm%coefs(iscat_a3,0:nleg_lay_cur,ilay) - 2.0 * definitionconversion(0:nleg_lay_cur) * forward_peak) + dot_product(gsf_2m2(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) - atm%coefs(iscat_a3,0:nleg_lay_cur,ilay)))
          atm%phase_ssg(iscat_a3,ilay,igeo) = 0.5*(dot_product(gsf_22(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) + atm%coefs(iscat_a3,0:nleg_lay_cur,ilay) - 2.0 * definitionconversion(0:nleg_lay_cur) * forward_peak) - dot_product(gsf_2m2(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) - atm%coefs(iscat_a3,0:nleg_lay_cur,ilay)))
        endif
        if (fullpol_nst(nst)) then
          atm%phase_ssg(iscat_b2,ilay,igeo) = dot_product(gsf_02(0:nleg_lay_cur,igeo) , atm%coefs(iscat_b2,0:nleg_lay_cur,ilay)) ! Off-diagonal, so no subtraction of the forward peak.
          atm%phase_ssg(iscat_a4,ilay,igeo) = dot_product(gsf_00(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a4,0:nleg_lay_cur,ilay) - definitionconversion(0:nleg_lay_cur) * forward_peak)
        endif
      enddo

    enddo

    set%deltam = .true.
    call lintran_execute_without_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,(/.true.,set%split_double,.true./),atm,set,self,rint,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif ! }}}

    ! Keep truncated phase function, but turn off delta-M.
    set%deltam = .false.
    call lintran_execute_without_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,(/.true.,.false.,.false./),atm,set,self,rint_new,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif ! }}}
    rint = rint - rint_new

    ! Clean up temporary heap information.
    deallocate(set%execute_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    if (atm_inp%thermal_emission) then
      deallocate(atm%planck_curve,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake deallocation.
    endif
    deallocate(atm%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%taua,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%taus,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1

  end subroutine lintran_nakajima_multiscattering_without_derivatives ! }}}

  !> Performs original Nakajima multiscattering with derivatives.
  !!
  !! This is a wrapper for the normal calculations. But now, multi-scattering is calculated
  !! using &delta;-M, but the contribution of single scattering is limited to real single
  !! scattering. In &delta;-M, scattering in the forward peak is not counted as a scattering
  !! order, giving lower scattering orders than what would be real. Using Lintran's &delta;-M
  !! is, in principle, better then this Nakajima wrapper, but then single-scattering must also be
  !! calculated with Lintran. Performing a single-scattering calculation on your own results in
  !! in wrong results, unless the user watches very carefully and then, it comes down to a
  !! single-scattering calculation with Lintran. This is the form with derivatives, which
  !! is called when the derivatives structure are passed.
  subroutine lintran_nakajima_multiscattering_with_derivatives(atm_inp,set_inp,self,rint,drv,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, errorflag_deallocation, warningflag_unpreparedderivatives, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class

    implicit none

    ! Input and output.
    type(lintran_atmosphere), intent(in) :: atm_inp !< Atmosphere provided by the user. Suffix _inp is added, because the Nakajima wrapper performs the desired calculation by making correction in the structure. The input structure remains unharmed.
    type(lintran_settings), intent(in) :: set_inp !< Calculation settings provided by the user. Suffix _inp is added, because the Nakajima wrapper performs the desired calculation by making correction in the structure. The input structure remains unharmed.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(self%nst,self%ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries.
    type(lintran_derivatives), intent(inout) :: drv !< Additional output, the structure of the derivatives.
    integer, intent(out) :: stat !< Error code.

    ! Dimension sizes to pass to internal routine.
    integer :: nleg_pert ! A dimension size that will be derived.
    integer :: nlay_deriv_base ! A dimension size that will be derived.
    integer :: nlay_deriv_ph ! A dimension size that will be derived.

    ! Flags for layer to be differentiated in any way, for differentiating C.
    logical, dimension(self%nlay) :: differentiated

    ! Radiance and derivatives for a new calculation.
    real, dimension(self%nst,self%ngeo) :: rint_new
    type(lintran_derivatives) :: drv_new

    ! Atmosphere and settings with possible hacks for Nakajima.
    type(lintran_atmosphere) :: atm
    type(lintran_settings) :: set

    ! Copy of dimension sizes.
    integer :: nst ! Number of Stokes parameters.
    integer :: nscat ! Number of relevant independent matrix elements in scattering matrix.
    integer :: nstrhalf ! Half the number of streams.
    integer :: nleg ! Highest Legendre coefficient.
    integer :: nlay ! Number of atmospheric layers.
    integer :: nleg_deriv ! Number of Legendre coefficients that are differentiated.
    integer :: ngeo ! Number of viewing geometries.
    ! Interpreted dimension size.
    integer :: nleg_pert_lay ! Number of Legendre polynomials perturbed in this layer.
    integer :: nleg_lay_cur ! Current number of Legendre coefficients in a layer, limited by the total number of Legendre coefficients.

    ! Whether or not differentiations take place at well.
    logical :: deriv_base
    logical :: deriv_ph

    ! Iterators.
    integer :: ileg ! Over Legendre polynomials.
    integer :: ilay ! Over atmospheric layers.
    integer :: ideriv ! Over differentiations.
    integer :: ist ! Over viewing Stokes parameters.
    integer :: igeo ! Over viewing geometries.

    ! Generalized spherical functions at single-scattering geometry.
    ! These are the two numbers as lower index of the generalized spherical function, the m denotes a minus sign for the number on the right sie.
    real, dimension(0:self%nleg,self%ngeo) :: gsf_00
    real, dimension(0:self%nleg,self%ngeo) :: gsf_02
    real, dimension(0:self%nleg,self%ngeo) :: gsf_22
    real, dimension(0:self%nleg,self%ngeo) :: gsf_2m2 

    real, dimension(0:self%nleg) :: definitionconversion ! Array of 2*ileg+1 numbers.

    ! Single-scattering geometry.
    real :: uscat_ssg

    ! Forward peak for artifical delta-M.
    real :: forward_peak

    ! The derivative with respect to single-scattering phase function may be absent, because
    ! for the user, there is no single-scattering. Internally, there is single-scattering,
    ! so its derivative should be calculated. But the user needn't do the administration.
    logical :: artificial_drint_phase_ssg ! Flag for that a field is allocated artificially.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Initialize error flag.
    stat = 0

    if (.not. self%flag_derivatives) then
      call lintran_nakajima_multiscattering_without_derivatives(atm_inp,set_inp,self,rint,stat)
      ! This routine is an interface call where the status flag is intent out. Therefore,
      ! this specific warning will be added at the end. Otherwise, it will be lost.
      stat = or(stat,warningflag_unpreparedderivatives)
      return
    endif

    ! Get dimension sizes on local scope.
    nst = self%nst
    nscat = self%nscat
    nstrhalf = self%nstrhalf
    nleg = self%nleg
    nlay = self%nlay
    ngeo = self%ngeo

    ! Dimension sizes needed for declaring an intermediate result in the execution.
    nleg_deriv = set_inp%nleg_deriv
    nleg_pert = min(self%nleg,nleg_deriv)
    ! These pack-statements allow extra zeros in the array and they will be ignored.
    ! But the zeros can only be put in the end. The values that make sense must
    ! be put in front of the array.
    nlay_deriv_base = size(pack(set_inp%ilay_deriv_base,set_inp%ilay_deriv_base .ne. 0))
    nlay_deriv_ph = size(pack(set_inp%ilay_deriv_ph,set_inp%ilay_deriv_ph .ne. 0))

    ! Whether differentiations take place.
    deriv_base = nlay_deriv_base .gt. 0
    deriv_ph = nlay_deriv_ph .gt. 0

    if (nlay_deriv_base .eq. 0 .and. nlay_deriv_ph .eq. 0) then
      differentiated = .false.
    endif
    if (nlay_deriv_base .eq. 0 .and. nlay_deriv_ph .ne. 0) then
      differentiated = (/(any(set_inp%ilay_deriv_ph(1:nlay_deriv_ph) .eq. ilay),ilay = 1,nlay)/)
    endif
    if (nlay_deriv_base .ne. 0 .and. nlay_deriv_ph .eq. 0) then
      differentiated = (/(any(set_inp%ilay_deriv_base(1:nlay_deriv_base) .eq. ilay),ilay = 1,nlay)/)
    endif
    if (nlay_deriv_base .ne. 0 .and. nlay_deriv_ph .ne. 0) then
      differentiated = (/(any(set_inp%ilay_deriv_base(1:nlay_deriv_base) .eq. ilay) .or. any(set_inp%ilay_deriv_ph(1:nlay_deriv_ph) .eq. ilay),ilay = 1,nlay)/)
    endif

    ! Copy atmosphere.
    allocate(atm%taus(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%taua(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%coefs(nscat,0:nleg+1,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%phase_ssg(nscat,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf(nst,nst,nstrhalf,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_0(nst,nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_v(nst,nst,nstrhalf,0:nleg,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%bdrf_ssg(nst,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%emi(nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(atm%emi_ssg(nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (atm_inp%thermal_emission) then
      allocate(atm%planck_curve(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 1 ! Fake allocation.
    endif
    allocate(atm%nleg_lay(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Copy atmosphere except for the phase function at single-scattering geometry. Meanwhile,
    ! the atmosphere may be encapsulated to the part that is used. Lintran used to cut off the
    ! used part of the atmosphere only in the pixel module, but that is because the arrays in
    ! the pixel module are allocated by Lintran, not by the user. This atmosphere is also
    ! allocated by Lintran, so the sizes must match. The number of atmospheric layers is not
    ! cut off, because truncating that will result in nonsense.
    atm%sun = atm_inp%sun
    atm%taus = atm_inp%taus
    atm%taua = atm_inp%taua
    atm%coefs = atm_inp%coefs(1:nscat,0:nleg+1,:) ! Possibly cut off non-used matrix elements and high Legendre coefficients.
    ! Single-scattering phase function not yet defined.
    atm%bdrf = atm_inp%bdrf(1:nst,1:nst,:,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%bdrf_0 = atm_inp%bdrf_0(1:nst,1:nst,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%bdrf_v = atm_inp%bdrf_v(1:nst,1:nst,:,0:nleg,1:ngeo) ! Possibly cut off non-used Stokes parameters, high Fourier coefficients and non-used viewing geometries.
    atm%bdrf_ssg = atm_inp%bdrf_ssg(1:nst,1:nst,1:ngeo) ! Possibly cut off non-used Stokes parameters and non-used viewing geometries.
    atm%emi = atm_inp%emi(1:nst,:,0:nleg) ! Possibly cut off non-used Stokes parameters and high Fourier coefficients.
    atm%emi_ssg = atm_inp%emi_ssg(1:nst,1:ngeo) ! Possibly cut off non-used Stokes parameters and non-used viewing geometries.
    if (atm_inp%thermal_emission) atm%planck_curve = atm_inp%planck_curve
    atm%nleg_lay = atm_inp%nleg_lay
    atm%bdrf_only_0 = atm_inp%bdrf_only_0
    atm%emi_only_0 = atm_inp%emi_only_0
    atm%thermal_emission = atm_inp%thermal_emission

    ! Copy settings.
    allocate(set%execute_stokes(nst),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    set%execute_stokes = set_inp%execute_stokes(1:nst)
    allocate(set%differentiate_stokes(nst),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_without_derivatives(progression,atm_inp%thermal_emission,atm,set,stat)
      return
    endif
    progression = progression + 1 ! }}}
    set%differentiate_stokes = set_inp%differentiate_stokes(1:nst)
    set%taus_split = set_inp%taus_split
    set%taua_split = set_inp%taua_split
    set%tautot_max = set_inp%tautot_max
    set%fourier_tolerance = set_inp%fourier_tolerance
    set%split_double = set_inp%split_double
    set%interpolation = set_inp%interpolation
    ! Delta-M and some stuff about derivatives are not copied.
    set%solver = set_inp%solver
    set%gs_maxiter = set_inp%gs_maxiter
    set%gs_tolerance = set_inp%gs_tolerance
    set%nleg_deriv = set_inp%nleg_deriv

    ! Settings of the derivatives.
    if (deriv_base) then
      allocate(set%ilay_deriv_base(nlay_deriv_base),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
      set%ilay_deriv_base = pack(set_inp%ilay_deriv_base,set_inp%ilay_deriv_base .ne. 0)
    else
      progression = progression + 1 ! Fake allocation.
    endif
    if (deriv_ph) then
      allocate(set%ilay_deriv_ph(nlay_deriv_ph),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
      set%ilay_deriv_ph = pack(set_inp%ilay_deriv_ph,set_inp%ilay_deriv_ph .ne. 0)
    else
      progression = progression + 1 ! Fake allocation.
    endif

    ! Construct temporary derivative field.
    if (deriv_base) then
      allocate(drv_new%drint_taus(nlay_deriv_base,nst,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(drv_new%drint_taua(nlay_deriv_base,nst,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
      if (atm_inp%thermal_emission) then
        allocate(drv_new%drint_planck_curve(nlay_deriv_base,nst,ngeo),stat=ierr) ! {{{
        if (ierr .ne. 0) then
          stat = or(stat,errorflag_allocation)
          call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
          return
        endif
        progression = progression + 1 ! }}}
      else
        progression = progression + 1 ! Fake allocation.
      endif
    else
      progression = progression + 3 ! Fake allocations.
    endif
    if (deriv_ph) then
      allocate(drv_new%drint_phase_ssg(nscat,nlay_deriv_ph,nst,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 1 ! Fake allocation.
    endif
    allocate(drv_new%drint_bdrf_ssg(nst,nst,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(drv_new%drint_emi(nst,nstrhalf,0:nleg,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(drv_new%drint_emi_ssg(nst,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! No need for non-single-scattering BDRF and Legendre coefficients, because this
    ! derivative is only used for single-scattering without delta-M.

    ! If not already there, allocate space for derivative with respect to single-scattering
    ! phase function in original field. If it did not exist beforehand, it will be destroyed
    ! afterwards.
    artificial_drint_phase_ssg = (.not. allocated(drv%drint_phase_ssg) .and. deriv_ph)

    ! Create field for derivatives with respect to single-scattering phase function,
    ! because Lintran thinks that single-scattering is calculated.
    if (artificial_drint_phase_ssg) then
      allocate(drv%drint_phase_ssg(nscat,nlay_deriv_ph,nst,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 1 ! Fake allocation.
    endif

    ! We desire to do array operations over Legendre coefficients, so the term
    ! (2*ileg + 1) must be put into an array. And we will directly transform it
    ! into a floating point.
    definitionconversion = (/(float(2*ileg+1),ileg=0,nleg)/)

    do igeo = 1,ngeo

      ! Calculate the single-scattering angle. This is the only reason to save the azimuth
      ! on the geometry structure.
      uscat_ssg = -self%geo%mu_0*self%geo%mu_v_allgeo(igeo) + sqrt((1.0-self%geo%mu_0**2.0)*(1.0-self%geo%mu_v_allgeo(igeo)**2.0)) * cos(self%geo%azi_allgeo(igeo))

      ! An artifical truncated single-scattering phase function is constructed, using
      ! the phase coefficients after delta-M correction. This involves a conversion from
      ! Legendre coefficients to a normal phase function, so this involves the generalized
      ! spherical functions for single-scattering geometry.

      gsf_00(0,igeo) = 1.0
      ! Use recursive definition of generalized spherical functions.
      do ileg = 1,nleg

        if (ileg .eq. 1) then
          ! Recursive construction, knowing that for l-2, the polynomials
          ! does not yet exist.
          gsf_00(ileg,igeo) = float(2*ileg-1)*uscat_ssg*gsf_00(ileg-1,igeo)/sqrt(float(ileg**2))
        else
          ! Full recursion.
          gsf_00(ileg,igeo) = (float(2*ileg-1)*uscat_ssg*gsf_00(ileg-1,igeo) - sqrt(float((ileg-1)**2))*gsf_00(ileg-2,igeo))/sqrt(float(ileg**2))
        endif

      enddo

      if (pol_nst(nst)) then ! Only for vector.

        gsf_02(0:1,igeo) = 0.0
        gsf_02(2,igeo) = -0.25 * sqrt(6.0) * (1.0 - uscat_ssg**2.0)
        gsf_22(0:1,igeo) = 0.0
        gsf_22(2,igeo) = 0.25 * (1.0 + uscat_ssg)**2.0
        gsf_2m2(0:1,igeo) = 0.0
        gsf_2m2(2,igeo) = 0.25 * (1.0 - uscat_ssg)**2.0

        do ileg = 2,nleg-1

          ! Here, m=0.
          gsf_02(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg) * gsf_02(ileg,igeo) - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0)) * gsf_02(ileg-1,igeo))
          ! Here m=2.
          gsf_22(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-4.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg-4.0) * gsf_22(ileg,igeo) - (ileg+1)*sqrt((ileg**2.0-4.0)*(ileg**2.0-4.0)) * gsf_22(ileg-1,igeo))
          gsf_2m2(ileg+1,igeo) = (1.0/(ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-4.0)))) * ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*uscat_ssg+4.0) * gsf_2m2(ileg,igeo) - (ileg+1)*sqrt((ileg**2.0-4.0)*(ileg**2.0-4.0)) * gsf_2m2(ileg-1,igeo))

        enddo

      endif

      ! Nakajima multiscattering is:
      ! + total-scattering with delta-M and truncated phase function
      ! - single-scattering without delta-M, with truncated phase function

      ! Total scattering with delta-M and truncated phase function.
      ! Construct truncated phase function.

      ! All existing elements are written, though in the physical world, only a1 and b1
      ! matter, because the sun radiates I.
      do ilay = 1,nlay
        ! Limit number of Legendre coefficients by the layer or by the model.
        nleg_lay_cur = min(nleg,atm%nleg_lay(ilay))
        if (atm%nleg_lay(ilay) .gt. nleg) then
          forward_peak = atm%coefs(iscat_a1,nleg+1,ilay) / (2.0*nleg+3.0)
        else
          forward_peak = 0.0
        endif
        atm%phase_ssg(iscat_a1,ilay,igeo) = dot_product(gsf_00(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a1,0:nleg_lay_cur,ilay) - definitionconversion(0:nleg_lay_cur) * forward_peak)
        if (pol_nst(nst)) then
          atm%phase_ssg(iscat_b1,ilay,igeo) = dot_product(gsf_02(0:nleg_lay_cur,igeo) , atm%coefs(iscat_b1,0:nleg_lay_cur,ilay)) ! Off-diagonal, so no subtraction of the forward peak.
          atm%phase_ssg(iscat_a2,ilay,igeo) = 0.5*(dot_product(gsf_22(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) + atm%coefs(iscat_a3,0:nleg_lay_cur,ilay) - 2.0 * definitionconversion(0:nleg_lay_cur) * forward_peak) + dot_product(gsf_2m2(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) - atm%coefs(iscat_a3,0:nleg_lay_cur,ilay)))
          atm%phase_ssg(iscat_a3,ilay,igeo) = 0.5*(dot_product(gsf_22(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) + atm%coefs(iscat_a3,0:nleg_lay_cur,ilay) - 2.0 * definitionconversion(0:nleg_lay_cur) * forward_peak) - dot_product(gsf_2m2(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a2,0:nleg_lay_cur,ilay) - atm%coefs(iscat_a3,0:nleg_lay_cur,ilay)))
        endif
        if (fullpol_nst(nst)) then
          atm%phase_ssg(iscat_b2,ilay,igeo) = dot_product(gsf_02(0:nleg_lay_cur,igeo) , atm%coefs(iscat_b2,0:nleg_lay_cur,ilay)) ! Off-diagonal, so no subtraction of the forward peak.
          atm%phase_ssg(iscat_a4,ilay,igeo) = dot_product(gsf_00(0:nleg_lay_cur,igeo) , atm%coefs(iscat_a4,0:nleg_lay_cur,ilay) - definitionconversion(0:nleg_lay_cur) * forward_peak)
        endif
      enddo

    enddo

    set%deltam = .true.
    call lintran_execute_with_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,nlay_deriv_base,nlay_deriv_ph,count(differentiated),pack((/(ilay,ilay=1,nlay)/),differentiated),nleg_pert,(/.true.,set%split_double,.true./),atm,set,self,rint,drv,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif ! }}}

    ! Keep truncated phase function, but turn off delta-M.
    set%deltam = .false.
    call lintran_execute_with_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,nlay_deriv_base,nlay_deriv_ph,count(differentiated),pack((/(ilay,ilay=1,nlay)/),differentiated),nleg_pert,(/.true.,.false.,.false./),atm,set,self,rint_new,drv_new,stat)
    ! Progress error. {{{
    if (and(stat,errorflag_allocation) .eq. errorflag_allocation) then
      call fail_nakajima_multiscattering_with_derivatives(progression,atm_inp%thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat)
      return
    endif ! }}}
    rint = rint - rint_new

    ! Aggregate derivatives. Note that the original derivatives structure may be too
    ! large, so indexing is needed to make sure the arrays have equal size.
    if (deriv_base) then
      drv%drint_taus(1:nlay_deriv_base,1:nst,1:ngeo) = drv%drint_taus(1:nlay_deriv_base,1:nst,1:ngeo) - drv_new%drint_taus
      drv%drint_taua(1:nlay_deriv_base,1:nst,1:ngeo) = drv%drint_taua(1:nlay_deriv_base,1:nst,1:ngeo) - drv_new%drint_taua
      if (atm_inp%thermal_emission) drv%drint_planck_curve(1:nlay_deriv_base,1:nst,1:ngeo) = drv%drint_planck_curve(1:nlay_deriv_base,1:nst,1:ngeo) - drv_new%drint_planck_curve
    endif

    ! In the case that dring_phase_ssg is not artifical, we have to subscript the array.
    ! Otherwise, it does not matter.
    if (deriv_ph) drv%drint_phase_ssg(1:nscat,1:nlay_deriv_ph,1:nst,1:ngeo) = drv%drint_phase_ssg(1:nscat,1:nlay_deriv_ph,1:nst,1:ngeo) - drv_new%drint_phase_ssg
    drv%drint_bdrf_ssg(1:nst,1:nst,1:nst,1:ngeo) = drv%drint_bdrf_ssg(1:nst,1:nst,1:nst,1:ngeo) - drv_new%drint_bdrf_ssg
    drv%drint_emi_ssg(1:nst,1:nst,1:ngeo) = drv%drint_emi_ssg(1:nst,1:nst,1:ngeo) - drv_new%drint_emi_ssg

    ! Translate derivative with respect to single-scattering phase function to
    ! derivatives with respect to phase coefficients using the generalized spherical
    ! functions.
    ! There are too many arrays with all different dimensions and differnet dimensions
    ! omitting. Even a magician would have difficulties matmulling and reshapes everything.
    if (deriv_ph) then
      do igeo = 1,ngeo
        do ist = 1,nst
          do ideriv = 1,nlay_deriv_ph
            ! This nleg_pert_lay already protects against higher numbers than nleg.
            nleg_pert_lay = min(nleg_pert,atm%nleg_lay(set%ilay_deriv_ph(ideriv)))
            drv%drint_coefs(iscat_a1,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_a1,0:nleg_pert_lay,ideriv,ist,igeo) + gsf_00(0:nleg_pert_lay,igeo) * drv%drint_phase_ssg(iscat_a1,ideriv,ist,igeo)
            if (pol_nst(nst)) then
              drv%drint_coefs(iscat_b1,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_b1,0:nleg_pert_lay,ideriv,ist,igeo) + gsf_02(0:nleg_pert_lay,igeo) * drv%drint_phase_ssg(iscat_b1,ideriv,ist,igeo)
              drv%drint_coefs(iscat_a2,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_a2,0:nleg_pert_lay,ideriv,ist,igeo) + 0.5*(gsf_22(0:nleg_pert_lay,igeo) * (drv%drint_phase_ssg(iscat_a2,ideriv,ist,igeo) + drv%drint_phase_ssg(iscat_a3,ideriv,ist,igeo)) + gsf_2m2(0:nleg_pert_lay,igeo) * (drv%drint_phase_ssg(iscat_a2,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a3,ideriv,ist,igeo)))
              drv%drint_coefs(iscat_a3,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_a3,0:nleg_pert_lay,ideriv,ist,igeo) + 0.5*(gsf_22(0:nleg_pert_lay,igeo) * (drv%drint_phase_ssg(iscat_a2,ideriv,ist,igeo) + drv%drint_phase_ssg(iscat_a3,ideriv,ist,igeo)) - gsf_2m2(0:nleg_pert_lay,igeo) * (drv%drint_phase_ssg(iscat_a2,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a3,ideriv,ist,igeo)))
            endif
            if (fullpol_nst(ist)) then
              drv%drint_coefs(iscat_b2,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_b2,0:nleg_pert_lay,ideriv,ist,igeo) + gsf_02(0:nleg_pert_lay,igeo) * drv%drint_phase_ssg(iscat_b2,ideriv,ist,igeo)
              drv%drint_coefs(iscat_a4,0:nleg_pert_lay,ideriv,ist,igeo) = drv%drint_coefs(iscat_a4,0:nleg_pert_lay,ideriv,ist,igeo) + gsf_00(0:nleg_pert_lay,igeo) * drv%drint_phase_ssg(iscat_a4,ideriv,ist,igeo)
            endif
          enddo
        enddo
      enddo

      ! Derivative with respect to delta-M parameter, add chain rule from the artificially created
      ! single-scattering phase function.
      if (nleg_deriv .gt. nleg) then
        do igeo = 1,ngeo
          do ist = 1,nst
            do ideriv = 1,nlay_deriv_ph
              if (atm%nleg_lay(set%ilay_deriv_ph(ideriv)) .le. nleg) cycle
              drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) = drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a1,ideriv,ist,igeo) * dot_product(gsf_00(:,igeo) , definitionconversion) / (2.0*nleg + 3.0)
              if (pol_nst(nst)) then
                drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) = drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a2,ideriv,ist,igeo) * dot_product(gsf_22(:,igeo) , definitionconversion) / (2.0*nleg + 3.0)
                drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) = drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a3,ideriv,ist,igeo) * dot_product(gsf_22(:,igeo) , definitionconversion) / (2.0*nleg + 3.0)
              endif
              if (fullpol_nst(nst)) then
                drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) = drv%drint_coefs(iscat_a1,nleg+1,ideriv,ist,igeo) - drv%drint_phase_ssg(iscat_a4,ideriv,ist,igeo) * dot_product(gsf_00(:,igeo) , definitionconversion) / (2.0*nleg + 3.0)
              endif
            enddo
          enddo
        enddo
      endif

      ! Clean derivative with respect to phase function because it has been transferred.
      drv%drint_phase_ssg = 0.0

    endif

    ! Clean up temporary heap information.
    if (artificial_drint_phase_ssg) then
      deallocate(drv%drint_phase_ssg,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake allocation.
    endif
    deallocate(drv_new%drint_emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(drv_new%drint_emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(drv_new%drint_bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    if (deriv_ph) then
      deallocate(drv_new%drint_phase_ssg,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake deallocation.
    endif
    if (deriv_base) then
      if (atm_inp%thermal_emission) then
        deallocate(drv_new%drint_planck_curve,stat=ierr)
        if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
        progression = progression - 1
      else
        progression = progression - 1 ! Fake allocation.
      endif
      deallocate(drv_new%drint_taua,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
      deallocate(drv_new%drint_taus,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 3 ! Fake allocations.
    endif
    if (deriv_ph) then
      deallocate(set%ilay_deriv_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake deallocation.
    endif
    if (deriv_base) then
      deallocate(set%ilay_deriv_base,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake deallocation.
    endif
    deallocate(set%differentiate_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(set%execute_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    if (atm_inp%thermal_emission) then
      deallocate(atm%planck_curve,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      progression = progression - 1
    else
      progression = progression - 1 ! Fake deallocation.
    endif
    deallocate(atm%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%taua,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1
    deallocate(atm%taus,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    progression = progression - 1

  end subroutine lintran_nakajima_multiscattering_with_derivatives ! }}}

  !> Calculation without derivatives.
  !!
  !! This routine is called by lintran_calculate_without_derivatives, but some organizational
  !! information is added as arguments. The actual calculations are performed in this routine.
  subroutine lintran_execute_without_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,execution,atm,set,self,rint,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, ist_i, ist_u, ist_sun, ist_therm, solver_gauss_seidel, solver_lu_decomposition, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_class, lintran_atmosphere, lintran_settings
    use lintran_pixel_module, only: pixel_set
    use lintran_tlc_module, only: calculate_fourier_independent_tlc_parameters, calculate_relevant_single_scattering_c_parameters, calculate_fourier_dependent_tlc_parameters, choose_viewing_geometry_for_fourier_independent_tlc_parameters, calculate_fourier_dependent_tlc_parameters_for_viewing_geometry
    use lintran_fields_module, only: fields_dynamic_init, fields_dynamic_close
    use lintran_matrix_module, only: matrix_int, gauss_seidel, lu_decompose, lu_solve
    use lintran_interface_module, only: direct_1_elem_interface, direct_2_interface, transmit_interface, transmit_sp_interface, source_response_1_interface, source_response_1_sp_interface, source_response_1_int_interface, source_response_2_sp_addition_interface, scatter_1_sp_interface, innerproduct_interface
    use lintran_nst1_module, only: direct_1_elem_nst1 => direct_1_elem, direct_2_nst1 => direct_2, transmit_nst1 => transmit, transmit_sp_nst1 => transmit_sp, source_response_1_nst1 => source_response_1, source_response_1_sp_nst1 => source_response_1_sp, source_response_1_int_nst1 => source_response_1_int, source_response_2_sp_addition_nst1 => source_response_2_sp_addition, scatter_1_sp_nst1 => scatter_1_sp, innerproduct_nst1 => innerproduct
    use lintran_nst3_module, only: direct_1_elem_nst3 => direct_1_elem, direct_2_nst3 => direct_2, transmit_nst3 => transmit, transmit_sp_nst3 => transmit_sp, source_response_1_nst3 => source_response_1, source_response_1_sp_nst3 => source_response_1_sp, source_response_1_int_nst3 => source_response_1_int, source_response_2_sp_addition_nst3 => source_response_2_sp_addition, scatter_1_sp_nst3 => scatter_1_sp, innerproduct_nst3 => innerproduct
    use lintran_nst4_module, only: direct_1_elem_nst4 => direct_1_elem, direct_2_nst4 => direct_2, transmit_nst4 => transmit, transmit_sp_nst4 => transmit_sp, source_response_1_nst4 => source_response_1, source_response_1_sp_nst4 => source_response_1_sp, source_response_1_int_nst4 => source_response_1_int, source_response_2_sp_addition_nst4 => source_response_2_sp_addition, scatter_1_sp_nst4 => scatter_1_sp, innerproduct_nst4 => innerproduct

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nscat !< Number of independent matrix element in scattering matrix.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst,ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: m ! Fourier index.
    integer :: ist ! Over Stokes parameters.
    integer :: igeo ! Over viewing geometries.

    ! Other dimension sizes.
    integer :: nsplit_total ! Number of layers after layer splitting.
    integer :: nmat ! Size of the matrix.
    integer :: nband ! Size of the matrix bands.

    ! Fourier convergence flag.
    logical, dimension(ngeo) :: fourier_converged
    real, dimension(nst,ngeo) :: rint_fourier ! Fourier-only part of the radiance.

    ! Flags for BDRF, fluorescence and thermal emissivity.
    ! They will possibly be turned off for Fourier numbers beyond 0
    ! or entirely for the thermal emissivity.
    logical :: exist_bdrf
    logical :: exist_emi
    logical :: exist_therm

    ! Results of current Fourier number (or single-scattering) and current viewing geometry
    ! and given viewing Stokes parameter.
    real :: rint_cur

    ! Function pointers. The interface of the functions aer the same for all polarization
    ! versions, so we just take the scalar ones, but that could have been anyone.
    procedure(direct_1_elem_interface), pointer :: direct_1_elem
    procedure(direct_2_interface), pointer :: direct_2
    procedure(transmit_interface), pointer :: transmit
    procedure(transmit_sp_interface), pointer :: transmit_sp
    procedure(source_response_1_interface), pointer :: source_response_1
    procedure(source_response_1_sp_interface), pointer :: source_response_1_sp
    procedure(source_response_1_int_interface), pointer :: source_response_1_int
    procedure(source_response_2_sp_addition_interface), pointer :: source_response_2_sp_addition
    procedure(scatter_1_sp_interface), pointer :: scatter_1_sp
    procedure(innerproduct_interface), pointer :: innerproduct

    ! Import atmospheric information with derivatives, therefore the .true.
    call pixel_set(nst,nscat,nstrhalf,nleg,nlay,ngeo,execution,.false.,atm,set,self%pix)

    ! Now the layer-splitting is done, the final dimension sizes
    ! can be acquired.
    if (execution(3)) then
      nsplit_total = sum(self%tlc%nsplit)
      nband = 4*nst*nstrhalf
      nmat = 2*nst*nstrhalf*(nsplit_total+1)

      ! Create interpolated fields for this total number of split layer.
      call fields_dynamic_init(nst,nstrhalf,nsplit_total,.false.,self%fld,stat)
      if (and(stat,errorflag_allocation) .eq. errorflag_allocation) return ! Progress error.
      ! No fail routine is required, because there is only one constructrion and that is also
      ! the only place where things can go wrong.
    endif

    ! Calculate TLC-parameters that do not depend on the Fourier number. Those
    ! that depend on the Fourier number will be updated each iteration of the Fourier
    ! loop.
    call calculate_fourier_independent_tlc_parameters(nstrhalf,nlay,ngeo,.false.,execution,set%split_double,set%interpolation,self%grd,self%geo,self%pix,self%tlc,stat)

    ! Set the pointers to the correct functions.
    if (fullpol_nst(nst)) then ! Vector 4.
      direct_1_elem => direct_1_elem_nst4
      direct_2 => direct_2_nst4
      transmit => transmit_nst4
      transmit_sp => transmit_sp_nst4
      source_response_1 => source_response_1_nst4
      source_response_1_sp => source_response_1_sp_nst4
      source_response_1_int => source_response_1_int_nst4
      source_response_2_sp_addition => source_response_2_sp_addition_nst4
      scatter_1_sp => scatter_1_sp_nst4
      innerproduct => innerproduct_nst4
    else if (pol_nst(nst)) then ! Vector 3.
      direct_1_elem => direct_1_elem_nst3
      direct_2 => direct_2_nst3
      transmit => transmit_nst3
      transmit_sp => transmit_sp_nst3
      source_response_1 => source_response_1_nst3
      source_response_1_sp => source_response_1_sp_nst3
      source_response_1_int => source_response_1_int_nst3
      source_response_2_sp_addition => source_response_2_sp_addition_nst3
      scatter_1_sp => scatter_1_sp_nst3
      innerproduct => innerproduct_nst3
    else ! Scalar.
      direct_1_elem => direct_1_elem_nst1
      direct_2 => direct_2_nst1
      transmit => transmit_nst1
      transmit_sp => transmit_sp_nst1
      source_response_1 => source_response_1_nst1
      source_response_1_sp => source_response_1_sp_nst1
      source_response_1_int => source_response_1_int_nst1
      source_response_2_sp_addition => source_response_2_sp_addition_nst1
      scatter_1_sp => scatter_1_sp_nst1
      innerproduct => innerproduct_nst1
    endif

    ! Save flag for whether thermal emission is taken into account. It will be set to false
    ! when a Fourier number above 0 is considered. Only use this flag for radiative transfer
    ! routines, not for chain rules or working out derivatives (not that we are calculating
    ! derivatives here, but still), because it will be set to false inside the Fourier loop
    ! if execution 2 or 3 is turned on.
    exist_therm = self%pix%thermal_emission

    ! Initialize aggregate results to zero that are aggregated inside
    ! the Fourier loop.
    rint = 0.0

    ! Finite-element single-scattering calculation.
    if (execution(1)) then

      do igeo = 1,ngeo

        ! Choose geometry for Fourier-independent TLC-parameters. This routine is called
        ! more often than maybe optimally needed, but it only sets pointers. It performs
        ! no calculations. The boolean argument is the flag for derivatives.
        call choose_viewing_geometry_for_fourier_independent_tlc_parameters(nlay,.false.,igeo,self%geo,self%tlc)

        do ist = 1,nst

          if (.not. set%execute_stokes(ist)) cycle

          ! Initialize results to zero.
          rint_cur = 0.0

          ! Calculate the relevant C-parameter for this viewing geometry and viewing Stokes parameter.
          call calculate_relevant_single_scattering_c_parameters(igeo,ist_sun,ist,self%geo,self%pix,self%tlc)

          ! Perform the calculation.
          call direct_1_elem(nlay,exist_therm .and. ist .eq. ist_therm,self%tlc,rint_cur)

          ! Apply single-scattering degeneracy factor.
          rint(ist,igeo) = rint(ist,igeo) + rint_cur * self%geo%degenerate_ssg_allgeo(igeo)

        enddo

      enddo

    endif

    ! Fourier loop only executed for double- and multi-scattering, because those
    ! are the only Fourier-dependent cases.
    if (execution(2) .or. execution(3)) then

      ! Fourier loop. This loop will be cut off if the relative contribution to
      ! the radiance is less than fourier_tolerance. This contribution is checked
      ! before applying the degeneracy factor (dependent on m and the azimuth angle)
      ! because that contribution can be zero for some angles, concealing the fact
      ! that the Fourier loop is actually not yet converged.
      fourier_converged = .false.
      rint_fourier = 0.0

      ! Existence of BDRF and fluorescent emissivity. Could be turned off after
      ! first Fourier number (0), depending on the atmosphere.
      exist_bdrf = .true.
      exist_emi = .true.
      do m = 0,nleg

        ! Turn off BDRF and/or emissivity at m=1 (after m=0) if only Fourier number 0 exists.
        ! Always turn off thermal emission for higher Fourier numbers.
        if (m .eq. 1) then
          if (self%pix%bdrf_only_0) exist_bdrf = .false.
          if (self%pix%emi_only_0) exist_emi = .false.
          exist_therm = .false.
        endif

        ! Update Fourier-dependent TLC-parameters to current Fourier number.
        ! Instrument-related TLC-parameters are calculated when the relevant
        ! viewing geometry is considered.
        call calculate_fourier_dependent_tlc_parameters(nst,nstrhalf,nleg,nlay,execution,m,ist_sun,exist_bdrf,exist_emi,self%grd,self%geo,self%pix,self%tlc)

        ! Prepare matrix if necessary.
        if (execution(3)) then
          call matrix_int(nst,nstrhalf,nlay,nsplit_total,m .eq. 0,set%solver .eq. solver_gauss_seidel,exist_bdrf,self%tlc,self%fld%mat,self%fld%shifts,stat) ! L (matrix).
          ! Perform LU-decomposition if necessary.
          if (set%solver .eq. solver_lu_decomposition) call lu_decompose(nmat,nband,self%fld%mat,self%fld%shifts,stat)
        endif

        if (set%split_double) then

          ! Analytic double-scattering.
          if (execution(2)) then

            ! <r | t s t s t | q> = D_2 + <R_1 | T | Q_1>
            call source_response_1(nstrhalf,nlay,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%a1) ! | Q_1> in A1.
            call transmit(nstrhalf,nlay,.false.,self%fld%a1,self%tlc,self%fld%a2) ! T | Q_1> in A2.

            ! Wait a moment, because the loop over geometries has to be opened. And
            ! for multi-scattering, also some geometry-independent stuff has to be done
            ! The geometry-dependent TLC-parameters are really calculated, so we want
            ! to minimize the number of calls to that routine.

          endif

          ! Multi-scattering.
          if (execution(3)) then

            ! L | I_2+> = | Q_2>
            ! E_3+ = <R_1 | I_2+>

            call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%i1) ! | Q_1> in I1.
            call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_1> in I2.
            call scatter_1_sp(nstrhalf,nlay,nsplit_total,exist_bdrf,self%fld%i2,self%tlc,self%fld%i1) ! S T | Q_1> in I1.
            call source_response_2_sp_addition(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_therm,self%tlc,self%fld%i1) ! | Q_2> (only same layer) added to in I1.

            ! Solve the matrix.
            if (set%solver .eq. solver_gauss_seidel) then
              call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_2 (total)> in I2 as first-guess for | I_2+>.
              call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2,stat) ! | I_2+> in I2.
            endif
            if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2) ! | I_2+> in I2.

            ! Okay, this was all the geometry-independent stuff for multi-scattering. Now
            ! the loop over geometries can be opened.

          endif

        else

          ! Multi-scattering.
          if (execution(3)) then

            ! L | I_1+> = | Q_1>
            ! E_2+ = <R_1 | I_1+>
          
            call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%i1) ! | Q_1> in I1.

            ! Solve the matrix.
            if (set%solver .eq. solver_gauss_seidel) then
              call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_1> in I2 as first-guess for | I_1+>.
              call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2,stat) ! | I_1+> in I2.
            endif
            if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2) ! | I_1+> in I2.

          endif

        endif

        ! Now, the geometry-dependent stuff is going to happen. Both the double-scattering
        ! and the multi-scattering will be resumed after the geometry-dependent
        ! TLC-parameters have been calculated for this geometry. Note that these
        ! are all Fourier-dependent TLC-parameters, so they are used only one iteration
        ! of this loop.
        do igeo = 1,ngeo

          ! Skip if this geometry is already converged.
          if (fourier_converged(igeo)) cycle
          ! Fourier convergence can only be achieved after doing the second Fourier component.
          ! This is because for Q and U, the second Fourier component is usually much more
          ! significant than the first, so convergence should not be concluded after evaluating
          ! the first Fourier component.
          if (m .ge. 2) fourier_converged(igeo) = .true. ! Will be set to false at the end of loop when not yet converged.

          ! Also set the pointers to Fourier-independent instrument-dependent
          ! TLC-parameters, because a new geometry is active now.
          ! The boolean is the flag for derivatives.
          call choose_viewing_geometry_for_fourier_independent_tlc_parameters(nlay,.false.,igeo,self%geo,self%tlc)

          do ist = 1,nst

            ! Skip Stokes parameter if not selected
            if (.not. set%execute_stokes(ist)) cycle

            ! Skip Fourier number zero for Stokes parameters on the other half as the sun.
            ! They have zero degeneracy, but may result is zero-division in Gauss-Seidel
            ! and they cost time, for just a zero.
            ! Greater than Q means U or V. This expression works for any number of
            ! Stokes parameters, because ist_u as parameter always exists.
            if (m .eq. 0 .and. ((ist_sun .ge. ist_u) .neqv. (ist .ge. ist_u))) cycle

            ! Initialize results to zero.
            rint_cur = 0.0

            ! Calculate the Fourier-dependent TLC-parameters in which the viewing
            ! geometry is involved.
            call calculate_fourier_dependent_tlc_parameters_for_viewing_geometry(nst,nstrhalf,nleg,nlay,igeo,m,ist,exist_bdrf,self%grd,self%geo,self%pix,self%tlc)

            if (set%split_double) then

              ! Resume double-scattering.
              if (execution(2)) then

                ! Overview of existing (instrument-independent) fields:
                ! | Q_1> in A1 (though not necessary).
                ! T | Q_1> in A2.

                call direct_2(nstrhalf,nlay,exist_therm,self%tlc,rint_cur)
                call source_response_1(nstrhalf,nlay,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%a1) ! <R_1 | in A1.
                call innerproduct(nstrhalf,nlay,self%fld%a1,self%fld%a2,rint_cur)

              endif

              ! Resume multi-scattering.
              if (execution(3)) then

                ! Overview of existing (instrument-independent) fields:
                ! | Q_2> in I1 (though not necessary).
                ! | I_2+> in I2.

                call source_response_1_int(nstrhalf,nlay,nsplit_total,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i1) ! <R_1 | in I1.
                call innerproduct(nstrhalf,nsplit_total,self%fld%i1,self%fld%i2,rint_cur)

              endif

            else

              ! Resume multi-scattering.
              if (execution(3)) then

                call source_response_1_int(nstrhalf,nlay,nsplit_total,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i1) ! <R_1 | in I1.
                call innerproduct(nstrhalf,nsplit_total,self%fld%i1,self%fld%i2,rint_cur)

              endif

            endif

            ! Count result of current Fourier number.
            rint_fourier(ist,igeo) = rint_fourier(ist,igeo) + rint_cur * self%geo%degenerate_allgeo(m,ist,igeo)

            ! Fourier convergence check. Uses the conventions of Lintran 1, except for the
            ! degeneracy factor and the Nakajima single-scattering difference.
            ! This implies that convergence is only checked for Stokes parameter I. the reference
            ! radiance is just the Fourier-part of the up-to-now calculated radiance. In contrast
            ! to Lintran 1, the Fourier part does not include some single-scattering difference
            ! from the Nakajima wrapper. There is always a calculation up to and including Fourier
            ! number 2, because that is often a very important one for polarization.
            if (ist .eq. ist_i .and. fourier_converged(igeo) .and. abs(rint_fourier(ist,igeo)) .gt. 0.0) then ! Against zero-division.
              if (abs(rint_cur / rint_fourier(ist,igeo)) .ge. set%fourier_tolerance) fourier_converged(igeo) = .false.
            endif

          enddo

        enddo

        ! Quit Fourier loop if all geometries are converged.
        if (all(fourier_converged)) exit

      enddo

      ! Add up the Fourier part to the eventual single-scattering part.
      rint = rint + rint_fourier

    endif

    if (execution(3)) call fields_dynamic_close(.false.,self%fld,stat)

  end subroutine lintran_execute_without_derivatives ! }}}

  !> Calculation with derivatives.
  !!
  !! This routine is called by lintran_calculate_with_derivatives, but some organizational
  !! information is added as arguments. The actual calculations are performed in this routine.
  subroutine lintran_execute_with_derivatives(nst,nscat,nstrhalf,nleg,nlay,ngeo,nlay_deriv_base,nlay_deriv_ph,nlay_deriv_c,ilay_deriv_c,nleg_pert,execution,atm,set,self,rint,drv,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation, ist_i, ist_q, ist_u, ist_v, ist_sun, ist_therm, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, solver_gauss_seidel, solver_lu_decomposition, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_class, lintran_atmosphere, lintran_settings, lintran_derivatives
    use lintran_pixel_module, only: pixel_set
    use lintran_tlc_module, only: calculate_fourier_independent_tlc_parameters, calculate_relevant_single_scattering_c_parameters, calculate_fourier_dependent_tlc_parameters, choose_viewing_geometry_for_fourier_independent_tlc_parameters, calculate_fourier_dependent_tlc_parameters_for_viewing_geometry
    use lintran_fields_module, only: fields_dynamic_init, fields_dynamic_close
    use lintran_matrix_module, only: matrix_int, gauss_seidel, lu_decompose, lu_solve
    use lintran_interface_module, only: direct_1_elem_interface, direct_2_interface, transmit_interface, transmit_sp_interface, source_response_1_interface, source_response_1_sp_interface, source_response_1_int_interface, source_response_2_sp_addition_interface, scatter_1_sp_interface, scatter_1_int_interface, innerproduct_interface, perturbation_scat_direct_1_elem_interface, perturbation_scat_direct_2_interface, perturbation_scat_source_response_1_interface, perturbation_scat_source_response_1_sp_interface, perturbation_scat_source_response_1_int_interface, perturbation_scat_source_response_2_int_interface, perturbation_scat_scatter_1_int_interface, perturbation_scat_zip_interface, perturbation_tau_direct_1_elem_interface, perturbation_tau_direct_2_interface, perturbation_tau_mutualtransmit_interface, perturbation_tau_mutualtransmit_sp_interface, perturbation_tau_source_response_1_interface, perturbation_tau_source_response_1_sp_interface, perturbation_tau_source_response_1_int_interface, perturbation_tau_source_response_2_int_interface, perturbation_tau_scatter_1_int_interface, perturbation_abs_zip_interface, perturbation_bdrf_direct_1_elem_interface, perturbation_bdrf_source_response_1_interface, perturbation_bdrf_scatter_1_interface, perturbation_bdrf_zip_interface, perturbation_emi_direct_1_elem_interface, perturbation_emi_source_response_1_interface, perturbation_therm_direct_1_elem_interface, perturbation_therm_direct_2_interface, perturbation_therm_source_response_1_interface, perturbation_therm_source_response_1_sp_interface, perturbation_therm_source_response_1_int_interface, perturbation_therm_source_response_2_int_interface
    use lintran_nst1_module, only: direct_1_elem_nst1 => direct_1_elem, direct_2_nst1 => direct_2, transmit_nst1 => transmit, transmit_sp_nst1 => transmit_sp, source_response_1_nst1 => source_response_1, source_response_1_sp_nst1 => source_response_1_sp, source_response_1_int_nst1 => source_response_1_int, source_response_2_sp_addition_nst1 => source_response_2_sp_addition, scatter_1_sp_nst1 => scatter_1_sp, scatter_1_int_nst1 => scatter_1_int, innerproduct_nst1 => innerproduct, perturbation_scat_direct_1_elem_nst1 => perturbation_scat_direct_1_elem, perturbation_scat_direct_2_nst1 => perturbation_scat_direct_2, perturbation_scat_source_response_1_nst1 => perturbation_scat_source_response_1, perturbation_scat_source_response_1_sp_nst1 => perturbation_scat_source_response_1_sp, perturbation_scat_source_response_1_int_nst1 => perturbation_scat_source_response_1_int, perturbation_scat_source_response_2_int_nst1 => perturbation_scat_source_response_2_int, perturbation_scat_scatter_1_int_nst1 => perturbation_scat_scatter_1_int, perturbation_scat_zip_nst1 => perturbation_scat_zip, perturbation_tau_direct_1_elem_nst1 => perturbation_tau_direct_1_elem, perturbation_tau_direct_2_nst1 => perturbation_tau_direct_2, perturbation_tau_mutualtransmit_nst1 => perturbation_tau_mutualtransmit, perturbation_tau_mutualtransmit_sp_nst1 => perturbation_tau_mutualtransmit_sp, perturbation_tau_source_response_1_nst1 => perturbation_tau_source_response_1, perturbation_tau_source_response_1_sp_nst1 => perturbation_tau_source_response_1_sp, perturbation_tau_source_response_1_int_nst1 => perturbation_tau_source_response_1_int, perturbation_tau_source_response_2_int_nst1 => perturbation_tau_source_response_2_int, perturbation_tau_scatter_1_int_nst1 => perturbation_tau_scatter_1_int, perturbation_abs_zip_nst1 => perturbation_abs_zip, perturbation_bdrf_direct_1_elem_nst1 => perturbation_bdrf_direct_1_elem, perturbation_bdrf_source_response_1_nst1 => perturbation_bdrf_source_response_1, perturbation_bdrf_scatter_1_nst1 => perturbation_bdrf_scatter_1, perturbation_bdrf_zip_nst1 => perturbation_bdrf_zip, perturbation_emi_direct_1_elem_nst1 => perturbation_emi_direct_1_elem, perturbation_emi_source_response_1_nst1 => perturbation_emi_source_response_1, perturbation_therm_direct_1_elem_nst1 => perturbation_therm_direct_1_elem, perturbation_therm_direct_2_nst1 => perturbation_therm_direct_2, perturbation_therm_source_response_1_nst1 => perturbation_therm_source_response_1, perturbation_therm_source_response_1_sp_nst1 => perturbation_therm_source_response_1_sp, perturbation_therm_source_response_1_int_nst1 => perturbation_therm_source_response_1_int, perturbation_therm_source_response_2_int_nst1 => perturbation_therm_source_response_2_int
    use lintran_nst3_module, only: direct_1_elem_nst3 => direct_1_elem, direct_2_nst3 => direct_2, transmit_nst3 => transmit, transmit_sp_nst3 => transmit_sp, source_response_1_nst3 => source_response_1, source_response_1_sp_nst3 => source_response_1_sp, source_response_1_int_nst3 => source_response_1_int, source_response_2_sp_addition_nst3 => source_response_2_sp_addition, scatter_1_sp_nst3 => scatter_1_sp, scatter_1_int_nst3 => scatter_1_int, innerproduct_nst3 => innerproduct, perturbation_scat_direct_1_elem_nst3 => perturbation_scat_direct_1_elem, perturbation_scat_direct_2_nst3 => perturbation_scat_direct_2, perturbation_scat_source_response_1_nst3 => perturbation_scat_source_response_1, perturbation_scat_source_response_1_sp_nst3 => perturbation_scat_source_response_1_sp, perturbation_scat_source_response_1_int_nst3 => perturbation_scat_source_response_1_int, perturbation_scat_source_response_2_int_nst3 => perturbation_scat_source_response_2_int, perturbation_scat_scatter_1_int_nst3 => perturbation_scat_scatter_1_int, perturbation_scat_zip_nst3 => perturbation_scat_zip, perturbation_tau_direct_1_elem_nst3 => perturbation_tau_direct_1_elem, perturbation_tau_direct_2_nst3 => perturbation_tau_direct_2, perturbation_tau_mutualtransmit_nst3 => perturbation_tau_mutualtransmit, perturbation_tau_mutualtransmit_sp_nst3 => perturbation_tau_mutualtransmit_sp, perturbation_tau_source_response_1_nst3 => perturbation_tau_source_response_1, perturbation_tau_source_response_1_sp_nst3 => perturbation_tau_source_response_1_sp, perturbation_tau_source_response_1_int_nst3 => perturbation_tau_source_response_1_int, perturbation_tau_source_response_2_int_nst3 => perturbation_tau_source_response_2_int, perturbation_tau_scatter_1_int_nst3 => perturbation_tau_scatter_1_int, perturbation_abs_zip_nst3 => perturbation_abs_zip, perturbation_bdrf_direct_1_elem_nst3 => perturbation_bdrf_direct_1_elem, perturbation_bdrf_source_response_1_nst3 => perturbation_bdrf_source_response_1, perturbation_bdrf_scatter_1_nst3 => perturbation_bdrf_scatter_1, perturbation_bdrf_zip_nst3 => perturbation_bdrf_zip, perturbation_emi_direct_1_elem_nst3 => perturbation_emi_direct_1_elem, perturbation_emi_source_response_1_nst3 => perturbation_emi_source_response_1, perturbation_therm_direct_1_elem_nst3 => perturbation_therm_direct_1_elem, perturbation_therm_direct_2_nst3 => perturbation_therm_direct_2, perturbation_therm_source_response_1_nst3 => perturbation_therm_source_response_1, perturbation_therm_source_response_1_sp_nst3 => perturbation_therm_source_response_1_sp, perturbation_therm_source_response_1_int_nst3 => perturbation_therm_source_response_1_int, perturbation_therm_source_response_2_int_nst3 => perturbation_therm_source_response_2_int
    use lintran_nst4_module, only: direct_1_elem_nst4 => direct_1_elem, direct_2_nst4 => direct_2, transmit_nst4 => transmit, transmit_sp_nst4 => transmit_sp, source_response_1_nst4 => source_response_1, source_response_1_sp_nst4 => source_response_1_sp, source_response_1_int_nst4 => source_response_1_int, source_response_2_sp_addition_nst4 => source_response_2_sp_addition, scatter_1_sp_nst4 => scatter_1_sp, scatter_1_int_nst4 => scatter_1_int, innerproduct_nst4 => innerproduct, perturbation_scat_direct_1_elem_nst4 => perturbation_scat_direct_1_elem, perturbation_scat_direct_2_nst4 => perturbation_scat_direct_2, perturbation_scat_source_response_1_nst4 => perturbation_scat_source_response_1, perturbation_scat_source_response_1_sp_nst4 => perturbation_scat_source_response_1_sp, perturbation_scat_source_response_1_int_nst4 => perturbation_scat_source_response_1_int, perturbation_scat_source_response_2_int_nst4 => perturbation_scat_source_response_2_int, perturbation_scat_scatter_1_int_nst4 => perturbation_scat_scatter_1_int, perturbation_scat_zip_nst4 => perturbation_scat_zip, perturbation_tau_direct_1_elem_nst4 => perturbation_tau_direct_1_elem, perturbation_tau_direct_2_nst4 => perturbation_tau_direct_2, perturbation_tau_mutualtransmit_nst4 => perturbation_tau_mutualtransmit, perturbation_tau_mutualtransmit_sp_nst4 => perturbation_tau_mutualtransmit_sp, perturbation_tau_source_response_1_nst4 => perturbation_tau_source_response_1, perturbation_tau_source_response_1_sp_nst4 => perturbation_tau_source_response_1_sp, perturbation_tau_source_response_1_int_nst4 => perturbation_tau_source_response_1_int, perturbation_tau_source_response_2_int_nst4 => perturbation_tau_source_response_2_int, perturbation_tau_scatter_1_int_nst4 => perturbation_tau_scatter_1_int, perturbation_abs_zip_nst4 => perturbation_abs_zip, perturbation_bdrf_direct_1_elem_nst4 => perturbation_bdrf_direct_1_elem, perturbation_bdrf_source_response_1_nst4 => perturbation_bdrf_source_response_1, perturbation_bdrf_scatter_1_nst4 => perturbation_bdrf_scatter_1, perturbation_bdrf_zip_nst4 => perturbation_bdrf_zip, perturbation_emi_direct_1_elem_nst4 => perturbation_emi_direct_1_elem, perturbation_emi_source_response_1_nst4 => perturbation_emi_source_response_1, perturbation_therm_direct_1_elem_nst4 => perturbation_therm_direct_1_elem, perturbation_therm_direct_2_nst4 => perturbation_therm_direct_2, perturbation_therm_source_response_1_nst4 => perturbation_therm_source_response_1, perturbation_therm_source_response_1_sp_nst4 => perturbation_therm_source_response_1_sp, perturbation_therm_source_response_1_int_nst4 => perturbation_therm_source_response_1_int, perturbation_therm_source_response_2_int_nst4 => perturbation_therm_source_response_2_int

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nscat !< Number of independent matrix element in scattering matrix.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    integer, intent(in) :: nlay_deriv_base !< Number of layers for which derivatives with respect to scattering and absorption are calculated. Their indices are hidden in the settings structure.
    integer, intent(in) :: nlay_deriv_ph !< Number of layers for which derivatives with respect to phase elements are calculated. Their indices are hidden in the settings structure.
    integer, intent(in) :: nlay_deriv_c !< Number of layers for which any derivatives are calculaed, needed for recycling intermediate results.
    integer, dimension(nlay_deriv_c), intent(in) :: ilay_deriv_c !< Indices of layers for which any derivatives are calculaed, needed for recycling intermediate results.
    integer, intent(in) :: nleg_pert !< Highest Legendre moment to be differentiated.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance.
    real, dimension(nst,ngeo), intent(out) :: rint !< Main output, the radiance for all viewing geometries.
    type(lintran_derivatives), intent(inout), target :: drv !< Additional output, the structure of the derivatives.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: m ! Fourier index.
    integer :: ileg ! Over Legendre coefficients.
    integer :: ist ! Over Stokes parameters.
    integer :: ist_src ! Over Source Stokes parameters.
    integer :: istr ! Over streams.
    integer :: istr_src ! Over source streams.
    integer :: ideriv ! Over differentiations.
    integer :: ilay ! Over atmospheric layers.
    integer :: igeo ! Over viewing geometries.

    ! Possibility to turn off derivatives of base and/or phase if the no layers are
    ! selected and segmentation faults are not desired.
    logical :: deriv_base
    logical :: deriv_ph
    logical :: deriv_c

    ! Other dimension sizes.
    integer :: nsplit_total ! Number of layers after layer splitting.
    integer :: nmat ! Size of the matrix.
    integer :: nband ! Size of the matrix bands.
    integer :: nleg_pert_lay ! Number of perturbed Legendre coefficients in current layer.

    ! Fourier convergence flag.
    logical, dimension(ngeo) :: fourier_converged
    real, dimension(nst,ngeo) :: rint_fourier ! Fourier-only part of the radiance.

    ! Flags for BDRF, fluorescence and thermal emissivity.
    ! They will possibly be turned off for Fourier numbers beyond 0
    ! or entirely for the thermal emissivity.
    logical :: exist_bdrf
    logical :: exist_emi
    logical :: exist_therm

    ! Results of current Fourier number (or single-scattering) and current viewing geometry.
    real :: rint_cur
    real, dimension(nlay_deriv_base) :: drint_ssa_cur
    real, dimension(nlay_deriv_base) :: drint_tau_cur
    real, dimension(nlay_deriv_base) :: drint_c_therm_cur
    real, dimension(nscat,0:nleg_pert,nlay_deriv_ph) :: drint_coefs_cur
    real, dimension(nlay) :: drint_c_ssg_cur
    real, dimension(:,:,:,:), pointer :: drint_bdrf_cur
    real, dimension(:,:), pointer :: drint_bdrf_0_cur
    real, dimension(:,:), pointer :: drint_bdrf_v_cur
    real :: drint_bdrf_ssg_cur
    real, dimension(:,:), pointer :: drint_emi_cur
    real :: drint_emi_ssg_cur

    ! Derivatives with respect to C, food for chain rules for phase functions
    ! and the single-scattering albedo.

    ! Temporary derivative with respect to the single-scattering phase function.
    real, dimension(nlay) :: drint_c_ssg
    real, dimension(nst,nstrhalf,nlay) :: drint_c_fwd_0
    real, dimension(nst,nstrhalf,nlay) :: drint_c_bck_0
    real, dimension(nst,nstrhalf,nlay) :: drint_c_fwd_v
    real, dimension(nst,nstrhalf,nlay) :: drint_c_bck_v
    real, dimension(nst,nst,nstrhalf,nstrhalf,nlay) :: drint_c_fwd
    real, dimension(nst,nst,nstrhalf,nstrhalf,nlay) :: drint_c_bck

    ! Derivative with respect to thermal C-parameter, to be progressed
    ! to the Planck curve and the optical depths.
    real, dimension(nlay_deriv_base,nst,ngeo) :: drint_c_therm

    ! Derivatives not yet applying the delta-M chain rules.
    real, dimension(nlay_deriv_base,nst,ngeo) :: drint_ssa_predm
    real, dimension(nlay_deriv_base,nst,ngeo) :: drint_tau_predm
    real, dimension(nscat,0:nleg_pert,nlay_deriv_ph,nst,ngeo) :: drint_coefs_predm
    real, dimension(nscat,nlay_deriv_ph,nst,ngeo) :: drint_phase_ssg_predm

    ! After delta-M correction, but still in tau/ssa coordinates.
    real, dimension(nlay_deriv_base,nst,ngeo) :: drint_ssa
    real, dimension(nlay_deriv_base,nst,ngeo) :: drint_tau

    ! Symmetry relationships for converting the nst by nst derivatives to phase matrices
    ! to derivatives with respect of coefficients.
    real :: sym_even
    real :: sym_odd

    ! Function pointers. The interface of the functions aer the same for all polarization
    ! versions, so we just take the scalar ones, but that could have been anyone.
    procedure(direct_1_elem_interface), pointer :: direct_1_elem
    procedure(direct_2_interface), pointer :: direct_2
    procedure(transmit_interface), pointer :: transmit
    procedure(transmit_sp_interface), pointer :: transmit_sp
    procedure(source_response_1_interface), pointer :: source_response_1
    procedure(source_response_1_sp_interface), pointer :: source_response_1_sp
    procedure(source_response_1_int_interface), pointer :: source_response_1_int
    procedure(source_response_2_sp_addition_interface), pointer :: source_response_2_sp_addition
    procedure(scatter_1_sp_interface), pointer :: scatter_1_sp
    procedure(scatter_1_int_interface), pointer :: scatter_1_int
    procedure(innerproduct_interface), pointer :: innerproduct
    procedure(perturbation_scat_direct_1_elem_interface), pointer :: perturbation_scat_direct_1_elem
    procedure(perturbation_scat_direct_2_interface), pointer :: perturbation_scat_direct_2
    procedure(perturbation_scat_source_response_1_interface), pointer :: perturbation_scat_source_response_1
    procedure(perturbation_scat_source_response_1_sp_interface), pointer :: perturbation_scat_source_response_1_sp
    procedure(perturbation_scat_source_response_1_int_interface), pointer :: perturbation_scat_source_response_1_int
    procedure(perturbation_scat_source_response_2_int_interface), pointer :: perturbation_scat_source_response_2_int
    procedure(perturbation_scat_scatter_1_int_interface), pointer :: perturbation_scat_scatter_1_int
    procedure(perturbation_scat_zip_interface), pointer :: perturbation_scat_zip
    procedure(perturbation_tau_direct_1_elem_interface), pointer :: perturbation_tau_direct_1_elem
    procedure(perturbation_tau_direct_2_interface), pointer :: perturbation_tau_direct_2
    procedure(perturbation_tau_mutualtransmit_interface), pointer :: perturbation_tau_mutualtransmit
    procedure(perturbation_tau_mutualtransmit_sp_interface), pointer :: perturbation_tau_mutualtransmit_sp
    procedure(perturbation_tau_source_response_1_interface), pointer :: perturbation_tau_source_response_1
    procedure(perturbation_tau_source_response_1_sp_interface), pointer :: perturbation_tau_source_response_1_sp
    procedure(perturbation_tau_source_response_1_int_interface), pointer :: perturbation_tau_source_response_1_int
    procedure(perturbation_tau_source_response_2_int_interface), pointer :: perturbation_tau_source_response_2_int
    procedure(perturbation_tau_scatter_1_int_interface), pointer :: perturbation_tau_scatter_1_int
    procedure(perturbation_abs_zip_interface), pointer :: perturbation_abs_zip
    procedure(perturbation_bdrf_direct_1_elem_interface), pointer :: perturbation_bdrf_direct_1_elem
    procedure(perturbation_bdrf_source_response_1_interface), pointer :: perturbation_bdrf_source_response_1
    procedure(perturbation_bdrf_scatter_1_interface), pointer :: perturbation_bdrf_scatter_1
    procedure(perturbation_bdrf_zip_interface), pointer :: perturbation_bdrf_zip
    procedure(perturbation_emi_direct_1_elem_interface), pointer :: perturbation_emi_direct_1_elem
    procedure(perturbation_emi_source_response_1_interface), pointer :: perturbation_emi_source_response_1
    procedure(perturbation_therm_direct_1_elem_interface), pointer :: perturbation_therm_direct_1_elem
    procedure(perturbation_therm_direct_2_interface), pointer :: perturbation_therm_direct_2
    procedure(perturbation_therm_source_response_1_interface), pointer :: perturbation_therm_source_response_1
    procedure(perturbation_therm_source_response_1_sp_interface), pointer :: perturbation_therm_source_response_1_sp
    procedure(perturbation_therm_source_response_1_int_interface), pointer :: perturbation_therm_source_response_1_int
    procedure(perturbation_therm_source_response_2_int_interface), pointer :: perturbation_therm_source_response_2_int

    ! Check if differentations are necessary.
    deriv_base = nlay_deriv_base .gt. 0
    deriv_ph = nlay_deriv_ph .gt. 0
    deriv_c = nlay_deriv_c .gt. 0

    ! Import atmospheric information with derivatives, therefore the .true.
    call pixel_set(nst,nscat,nstrhalf,nleg,nlay,ngeo,execution,.true.,atm,set,self%pix)

    ! Now the layer-splitting is done, the final dimension sizes
    ! can be acquired.
    if (execution(3)) then
      nsplit_total = sum(self%tlc%nsplit)
      nband = 4*nst*nstrhalf
      nmat = 2*nst*nstrhalf*(nsplit_total+1)

      ! Create interpolated fields for this total number of split layer.
      call fields_dynamic_init(nst,nstrhalf,nsplit_total,.true.,self%fld,stat)
      if (and(stat,errorflag_allocation) .eq. errorflag_allocation) return ! Progress error.
      ! No fail routine is required, because there is only one constructrion and that is also
      ! the only place where things can go wrong.
    endif

    ! Calculate TLC-parameters that do not depend on the Fourier number. Those
    ! that depend on the Fourier number will be updated each iteration of the Fourier
    ! loop.
    call calculate_fourier_independent_tlc_parameters(nstrhalf,nlay,ngeo,.true.,execution,set%split_double,set%interpolation,self%grd,self%geo,self%pix,self%tlc,stat,set%ilay_deriv_base(1:nlay_deriv_base))

    ! Set the pointers to the correct functions.
    if (fullpol_nst(nst)) then ! Vector 4.
      direct_1_elem => direct_1_elem_nst4
      direct_2 => direct_2_nst4
      transmit => transmit_nst4
      transmit_sp => transmit_sp_nst4
      source_response_1 => source_response_1_nst4
      source_response_1_sp => source_response_1_sp_nst4
      source_response_1_int => source_response_1_int_nst4
      source_response_2_sp_addition => source_response_2_sp_addition_nst4
      scatter_1_sp => scatter_1_sp_nst4
      scatter_1_int => scatter_1_int_nst4
      innerproduct => innerproduct_nst4
      perturbation_scat_direct_1_elem => perturbation_scat_direct_1_elem_nst4
      perturbation_scat_direct_2 => perturbation_scat_direct_2_nst4
      perturbation_scat_source_response_1 => perturbation_scat_source_response_1_nst4
      perturbation_scat_source_response_1_sp => perturbation_scat_source_response_1_sp_nst4
      perturbation_scat_source_response_1_int => perturbation_scat_source_response_1_int_nst4
      perturbation_scat_source_response_2_int => perturbation_scat_source_response_2_int_nst4
      perturbation_scat_scatter_1_int => perturbation_scat_scatter_1_int_nst4
      perturbation_scat_zip => perturbation_scat_zip_nst4
      perturbation_tau_direct_1_elem => perturbation_tau_direct_1_elem_nst4
      perturbation_tau_direct_2 => perturbation_tau_direct_2_nst4
      perturbation_tau_mutualtransmit => perturbation_tau_mutualtransmit_nst4
      perturbation_tau_mutualtransmit_sp => perturbation_tau_mutualtransmit_sp_nst4
      perturbation_tau_source_response_1 => perturbation_tau_source_response_1_nst4
      perturbation_tau_source_response_1_sp => perturbation_tau_source_response_1_sp_nst4
      perturbation_tau_source_response_1_int => perturbation_tau_source_response_1_int_nst4
      perturbation_tau_source_response_2_int => perturbation_tau_source_response_2_int_nst4
      perturbation_tau_scatter_1_int => perturbation_tau_scatter_1_int_nst4
      perturbation_abs_zip => perturbation_abs_zip_nst4
      perturbation_bdrf_direct_1_elem => perturbation_bdrf_direct_1_elem_nst4
      perturbation_bdrf_source_response_1 => perturbation_bdrf_source_response_1_nst4
      perturbation_bdrf_scatter_1 => perturbation_bdrf_scatter_1_nst4
      perturbation_bdrf_zip => perturbation_bdrf_zip_nst4
      perturbation_emi_direct_1_elem => perturbation_emi_direct_1_elem_nst4
      perturbation_emi_source_response_1 => perturbation_emi_source_response_1_nst4
      perturbation_therm_direct_1_elem => perturbation_therm_direct_1_elem_nst4
      perturbation_therm_direct_2 => perturbation_therm_direct_2_nst4
      perturbation_therm_source_response_1 => perturbation_therm_source_response_1_nst4
      perturbation_therm_source_response_1_sp => perturbation_therm_source_response_1_sp_nst4
      perturbation_therm_source_response_1_int => perturbation_therm_source_response_1_int_nst4
      perturbation_therm_source_response_2_int => perturbation_therm_source_response_2_int_nst4
    else if (pol_nst(nst)) then ! Vector 3.
      direct_1_elem => direct_1_elem_nst3
      direct_2 => direct_2_nst3
      transmit => transmit_nst3
      transmit_sp => transmit_sp_nst3
      source_response_1 => source_response_1_nst3
      source_response_1_sp => source_response_1_sp_nst3
      source_response_1_int => source_response_1_int_nst3
      source_response_2_sp_addition => source_response_2_sp_addition_nst3
      scatter_1_sp => scatter_1_sp_nst3
      scatter_1_int => scatter_1_int_nst3
      innerproduct => innerproduct_nst3
      perturbation_scat_direct_1_elem => perturbation_scat_direct_1_elem_nst3
      perturbation_scat_direct_2 => perturbation_scat_direct_2_nst3
      perturbation_scat_source_response_1 => perturbation_scat_source_response_1_nst3
      perturbation_scat_source_response_1_sp => perturbation_scat_source_response_1_sp_nst3
      perturbation_scat_source_response_1_int => perturbation_scat_source_response_1_int_nst3
      perturbation_scat_source_response_2_int => perturbation_scat_source_response_2_int_nst3
      perturbation_scat_scatter_1_int => perturbation_scat_scatter_1_int_nst3
      perturbation_scat_zip => perturbation_scat_zip_nst3
      perturbation_tau_direct_1_elem => perturbation_tau_direct_1_elem_nst3
      perturbation_tau_direct_2 => perturbation_tau_direct_2_nst3
      perturbation_tau_mutualtransmit => perturbation_tau_mutualtransmit_nst3
      perturbation_tau_mutualtransmit_sp => perturbation_tau_mutualtransmit_sp_nst3
      perturbation_tau_source_response_1 => perturbation_tau_source_response_1_nst3
      perturbation_tau_source_response_1_sp => perturbation_tau_source_response_1_sp_nst3
      perturbation_tau_source_response_1_int => perturbation_tau_source_response_1_int_nst3
      perturbation_tau_source_response_2_int => perturbation_tau_source_response_2_int_nst3
      perturbation_tau_scatter_1_int => perturbation_tau_scatter_1_int_nst3
      perturbation_abs_zip => perturbation_abs_zip_nst3
      perturbation_bdrf_direct_1_elem => perturbation_bdrf_direct_1_elem_nst3
      perturbation_bdrf_source_response_1 => perturbation_bdrf_source_response_1_nst3
      perturbation_bdrf_scatter_1 => perturbation_bdrf_scatter_1_nst3
      perturbation_bdrf_zip => perturbation_bdrf_zip_nst3
      perturbation_emi_direct_1_elem => perturbation_emi_direct_1_elem_nst3
      perturbation_emi_source_response_1 => perturbation_emi_source_response_1_nst3
      perturbation_therm_direct_1_elem => perturbation_therm_direct_1_elem_nst3
      perturbation_therm_direct_2 => perturbation_therm_direct_2_nst3
      perturbation_therm_source_response_1 => perturbation_therm_source_response_1_nst3
      perturbation_therm_source_response_1_sp => perturbation_therm_source_response_1_sp_nst3
      perturbation_therm_source_response_1_int => perturbation_therm_source_response_1_int_nst3
      perturbation_therm_source_response_2_int => perturbation_therm_source_response_2_int_nst3
    else ! Scalar.
      direct_1_elem => direct_1_elem_nst1
      direct_2 => direct_2_nst1
      transmit => transmit_nst1
      transmit_sp => transmit_sp_nst1
      source_response_1 => source_response_1_nst1
      source_response_1_sp => source_response_1_sp_nst1
      source_response_1_int => source_response_1_int_nst1
      source_response_2_sp_addition => source_response_2_sp_addition_nst1
      scatter_1_sp => scatter_1_sp_nst1
      scatter_1_int => scatter_1_int_nst1
      innerproduct => innerproduct_nst1
      perturbation_scat_direct_1_elem => perturbation_scat_direct_1_elem_nst1
      perturbation_scat_direct_2 => perturbation_scat_direct_2_nst1
      perturbation_scat_source_response_1 => perturbation_scat_source_response_1_nst1
      perturbation_scat_source_response_1_sp => perturbation_scat_source_response_1_sp_nst1
      perturbation_scat_source_response_1_int => perturbation_scat_source_response_1_int_nst1
      perturbation_scat_source_response_2_int => perturbation_scat_source_response_2_int_nst1
      perturbation_scat_scatter_1_int => perturbation_scat_scatter_1_int_nst1
      perturbation_scat_zip => perturbation_scat_zip_nst1
      perturbation_tau_direct_1_elem => perturbation_tau_direct_1_elem_nst1
      perturbation_tau_direct_2 => perturbation_tau_direct_2_nst1
      perturbation_tau_mutualtransmit => perturbation_tau_mutualtransmit_nst1
      perturbation_tau_mutualtransmit_sp => perturbation_tau_mutualtransmit_sp_nst1
      perturbation_tau_source_response_1 => perturbation_tau_source_response_1_nst1
      perturbation_tau_source_response_1_sp => perturbation_tau_source_response_1_sp_nst1
      perturbation_tau_source_response_1_int => perturbation_tau_source_response_1_int_nst1
      perturbation_tau_source_response_2_int => perturbation_tau_source_response_2_int_nst1
      perturbation_tau_scatter_1_int => perturbation_tau_scatter_1_int_nst1
      perturbation_abs_zip => perturbation_abs_zip_nst1
      perturbation_bdrf_direct_1_elem => perturbation_bdrf_direct_1_elem_nst1
      perturbation_bdrf_source_response_1 => perturbation_bdrf_source_response_1_nst1
      perturbation_bdrf_scatter_1 => perturbation_bdrf_scatter_1_nst1
      perturbation_bdrf_zip => perturbation_bdrf_zip_nst1
      perturbation_emi_direct_1_elem => perturbation_emi_direct_1_elem_nst1
      perturbation_emi_source_response_1 => perturbation_emi_source_response_1_nst1
      perturbation_therm_direct_1_elem => perturbation_therm_direct_1_elem_nst1
      perturbation_therm_direct_2 => perturbation_therm_direct_2_nst1
      perturbation_therm_source_response_1 => perturbation_therm_source_response_1_nst1
      perturbation_therm_source_response_1_sp => perturbation_therm_source_response_1_sp_nst1
      perturbation_therm_source_response_1_int => perturbation_therm_source_response_1_int_nst1
      perturbation_therm_source_response_2_int => perturbation_therm_source_response_2_int_nst1
    endif

    ! Initialize result and derivatives to zero. Some of those may not be touched, for instance
    ! the derivative with respect to single-scattering phase function if single-scattering
    ! is turned off. In such a case, it is possible to have those fields not allocated
    ! in the derivatives structure. If they are allocated, they will be set to zero.
    ! If fields in the derivatives are not allocated, but they are needed, a segmentation
    ! fault will occur, just like if any necessary field is not allocated.

    ! Therefore, we use the `allocated'-statement, because the user outside Lintran decides
    ! whether or not the fields are allocated. Though some of these should certainly be
    ! allocated, we decide to save the administrative misery and do the same for every field.

    ! Save thermal emission flag, to be used for single-scattering and Fourier number
    ! zero for multi-scattering. When higher Fourier numbers are calculated, the flag
    ! is automatically set to false. Do not use this flag for administrative actions
    ! such as chain rules for derivatives, because it ends up as false when execution 2
    ! or 3 is turned on.
    exist_therm = self%pix%thermal_emission

    rint = 0.0
    if (allocated(drv%drint_taus)) drv%drint_taus = 0.0
    if (allocated(drv%drint_taua)) drv%drint_taua = 0.0
    if (allocated(drv%drint_coefs)) drv%drint_coefs = 0.0
    if (allocated(drv%drint_phase_ssg)) drv%drint_phase_ssg = 0.0
    if (allocated(drv%drint_bdrf)) drv%drint_bdrf = 0.0
    if (allocated(drv%drint_bdrf_0)) drv%drint_bdrf_0 = 0.0
    if (allocated(drv%drint_bdrf_v)) drv%drint_bdrf_v = 0.0
    if (allocated(drv%drint_bdrf_ssg)) drv%drint_bdrf_ssg = 0.0
    if (allocated(drv%drint_emi)) drv%drint_emi = 0.0
    if (allocated(drv%drint_emi_ssg)) drv%drint_emi_ssg = 0.0
    if (allocated(drv%drint_planck_curve)) drv%drint_planck_curve = 0.0

    ! Initialize aggregate intermediate results to zero.
    drint_ssa_predm = 0.0
    drint_tau_predm = 0.0
    drint_coefs_predm = 0.0
    drint_phase_ssg_predm = 0.0
    drint_c_therm = 0.0

    ! Finite-element single-scattering calculation.
    if (execution(1)) then

      do igeo = 1,ngeo

        ! Choose geometry for Fourier-independent TLC-parameters. This routine is called
        ! more often than maybe optimally needed, but it only sets pointers. It performs
        ! no calculations. The boolean argument is the flag for derivatives.
        call choose_viewing_geometry_for_fourier_independent_tlc_parameters(nlay,.true.,igeo,self%geo,self%tlc)

        do ist = 1,nst

          if (.not. set%execute_stokes(ist)) cycle

          ! Initialize intermediate result for single-scattering to zero.
          drint_c_ssg = 0.0

          ! Inirialize results before applying degeneracy factor.
          rint_cur = 0.0
          drint_c_ssg_cur = 0.0
          drint_tau_cur = 0.0
          drint_bdrf_ssg_cur = 0.0
          drint_emi_ssg_cur = 0.0
          drint_c_therm_cur = 0.0

          ! Calculate the relevant C-parameter for this viewing geometry and viewing Stokes parameter.
          call calculate_relevant_single_scattering_c_parameters(igeo,ist_sun,ist,self%geo,self%pix,self%tlc)

          ! Perform the calculation.
          call direct_1_elem(nlay,exist_therm .and. ist .eq. ist_therm,self%tlc,rint_cur)

          ! Apply degeneracy factors.
          rint(ist,igeo) = rint(ist,igeo) + rint_cur * self%geo%degenerate_ssg_allgeo(igeo)

          ! Derivatives.
          if (set%differentiate_stokes(ist)) then

            if (deriv_c) call perturbation_scat_direct_1_elem(nlay,nlay_deriv_c,ilay_deriv_c,self%tlc,drint_c_ssg_cur)
            if (deriv_base) call perturbation_tau_direct_1_elem(nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),exist_therm .and. ist .eq. ist_therm,self%tlc,drint_tau_cur)
            call perturbation_bdrf_direct_1_elem(nlay,self%tlc,drint_bdrf_ssg_cur)
            call perturbation_emi_direct_1_elem(nlay,self%tlc,drint_emi_ssg_cur)
            if (exist_therm .and. deriv_base .and. ist .eq. ist_therm) call perturbation_therm_direct_1_elem(nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%tlc,drint_c_therm_cur)

            ! Apply degeneracy factors.
            drint_c_ssg = drint_c_ssg + drint_c_ssg_cur * self%geo%degenerate_ssg_allgeo(igeo)
            drint_tau_predm(:,ist,igeo) = drint_tau_predm(:,ist,igeo) + drint_tau_cur * self%geo%degenerate_ssg_allgeo(igeo)
            drv%drint_bdrf_ssg(ist_sun,ist,ist,igeo) = drv%drint_bdrf_ssg(ist_sun,ist,ist,igeo) + drint_bdrf_ssg_cur * self%geo%degenerate_ssg_allgeo(igeo)
            drv%drint_emi_ssg(ist,ist,igeo) = drv%drint_emi_ssg(ist,ist,igeo) + drint_emi_ssg_cur * self%geo%degenerate_ssg_allgeo(igeo)
            if (exist_therm) drint_c_therm(:,ist,igeo) = drint_c_therm(:,ist,igeo) + drint_c_therm_cur * self%geo%degenerate_ssg_allgeo(igeo)

            ! Convert derivatives with respect to C to derivatives with respect to
            ! single-scattering albedo and the phase matrix elements for single scattering.

            ! Single-scattering albedo.
            if (deriv_base) drint_ssa_predm(:,ist,igeo) = drint_ssa_predm(:,ist,igeo) + drint_c_ssg(set%ilay_deriv_base(1:nlay_deriv_base)) * self%tlc%diff_c_bck_0v_elem_ssa(set%ilay_deriv_base(1:nlay_deriv_base))

            ! Phase matrix.
            if (deriv_ph) then

              if (ist_sun .eq. ist_i) then

                if (ist .eq. ist_i) drint_phase_ssg_predm(iscat_a1,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a1,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))
                if (pol_nst(nst)) then
                  if (ist .eq. ist_q) drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_v_allgeo(igeo)
                  if (ist .eq. ist_u) drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_v_allgeo(igeo)
                endif
                ! Nothing to be done for V

              endif

              if (ist_sun .eq. ist_q) then

                ! Here we assume at least vector code. Otherwise, it is a problem that the
                ! sun is radiating Q instead of I. With vector, that is of course perfectly
                ! fine.
                if (ist .eq. ist_i) drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)
                if (ist .eq. ist_q) then
                  drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)*self%geo%rotmat_c_v_allgeo(igeo)
                  drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)*self%geo%rotmat_s_v_allgeo(igeo)
                endif
                if (ist .eq. ist_u) then
                  drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)*self%geo%rotmat_s_v_allgeo(igeo)
                  drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)*self%geo%rotmat_c_v_allgeo(igeo)
                endif
                if (fullpol_nst(nst)) then
                  if (ist .eq. ist_v) drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)
                endif

              endif

              if (ist_sun .eq. ist_u) then

                ! The same applies as for a sun radiating Q. This must be at least vector
                ! radiative transport.
                if (ist .eq. ist_i) drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b1,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)
                if (ist .eq. ist_q) then
                  drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)*self%geo%rotmat_c_v_allgeo(igeo)
                  drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)*self%geo%rotmat_s_v_allgeo(igeo)
                endif
                if (ist .eq. ist_u) then
                  drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a2,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_0_allgeo(igeo)*self%geo%rotmat_s_v_allgeo(igeo)
                  drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a3,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)*self%geo%rotmat_c_v_allgeo(igeo)
                endif
                if (fullpol_nst(nst)) then
                  if (ist .eq. ist_v) drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_0_allgeo(igeo)
                endif

              endif

              if (ist_sun .eq. ist_v) then

                ! Here, we can assume that full polarization is chosen.

                ! Nothing to be done for I
                if (ist .eq. ist_q) drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) - self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_s_v_allgeo(igeo)
                if (ist .eq. ist_u) drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_b2,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))*self%geo%rotmat_c_v_allgeo(igeo)
                if (ist .eq. ist_v) drint_phase_ssg_predm(iscat_a4,1:nlay_deriv_ph,ist,igeo) = drint_phase_ssg_predm(iscat_a4,1:nlay_deriv_ph,ist,igeo) + self%tlc%diff_c_ph(set%ilay_deriv_ph(1:nlay_deriv_ph)) * drint_c_ssg(set%ilay_deriv_ph(1:nlay_deriv_ph))

              endif

            endif

          endif

        enddo

      enddo

    endif

    ! Fourier loop is only done for double-scattering or multi-scattering. Most
    ! inside the loop has an if-clause for an execution, the entire administration,
    ! e.g. calculating Fourier-dependent TLC-parameters.
    if (execution(2) .or. execution(3)) then

      ! Fourier loop. This loop will be cut off if the relative contribution to
      ! the radiance is less than fourier_tolerance. This contribution is checked
      ! before applying the degeneracy factor (dependent on m and the azimuth angle)
      ! because that contribution can be zero for some angles, concealing the fact
      ! that the Fourier loop is actually not yet converged.
      fourier_converged = .false.
      rint_fourier = 0.0

      ! Existence of BDRF and fluorescent emissivity. Could be turned off after
      ! first Fourier number (0), depending on the atmosphere.
      exist_bdrf = .true.
      exist_emi = .true.
      do m = 0,nleg

        ! Turn off BDRF and/or emissivity at m=1 (after m=0) if only Fourier number 0 exists.
        ! Always turn off thermal emission for higher Fourier numbers.
        if (m .eq. 1) then
          if (self%pix%bdrf_only_0) exist_bdrf = .false.
          if (self%pix%emi_only_0) exist_emi = .false.
          exist_therm = .false.
        endif

        ! Update Fourier-dependent TLC-parameters to current Fourier number.
        ! Instrument-related TLC-parameters are calculated when the relevant
        ! viewing geometry is considered.
        call calculate_fourier_dependent_tlc_parameters(nst,nstrhalf,nleg,nlay,execution,m,ist_sun,exist_bdrf,exist_emi,self%grd,self%geo,self%pix,self%tlc)

        ! Prepare matrix if necessary.
        if (execution(3)) then
          call matrix_int(nst,nstrhalf,nlay,nsplit_total,m .eq. 0,set%solver .eq. solver_gauss_seidel,exist_bdrf,self%tlc,self%fld%mat,self%fld%shifts,stat) ! L (matrix).
          ! Perform LU-decomposition if necessary.
          if (set%solver .eq. solver_lu_decomposition) call lu_decompose(nmat,nband,self%fld%mat,self%fld%shifts,stat)
        endif

        if (set%split_double) then

          ! Analytic double-scattering.
          if (execution(2)) then

            ! All calculation that do not depend on the instrument are calculated in
            ! the beginning. The instrument may have a lot of duplicates, because of
            ! different Stokes parameters and different viewing geometries.

            ! <r | t s t s t | q> = D_2 + <R_1 | T | Q_1>
            call source_response_1(nstrhalf,nlay,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%a1) ! | Q_1> in A1.
            call transmit(nstrhalf,nlay,.false.,self%fld%a1,self%tlc,self%fld%a2) ! T | Q_1> in A2.

            ! Wait a moment, because the loop over geometries has to be opened. And
            ! for multi-scattering, also some geometry-independent stuff has to be done
            ! The geometry-dependent TLC-parameters are really calculated, so we want
            ! to minimize the number of calls to that routine.

          endif

          ! Multi-scattering.
          if (execution(3)) then

            ! L | I_2+> = | Q_2>
            ! E_3+ = <R_1 | I_2+>
            call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%i1) ! | Q_1> in I1.
            call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_1> in I2.
            call scatter_1_sp(nstrhalf,nlay,nsplit_total,exist_bdrf,self%fld%i2,self%tlc,self%fld%i1) ! S T | Q_1> in I1.
            call source_response_2_sp_addition(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_therm,self%tlc,self%fld%i1) ! | Q_2> (only same layer) added to in I1.

            ! Solve the matrix.
            if (set%solver .eq. solver_gauss_seidel) then
              call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_2 (total)> in I2 as first-guess for | I_2+>.
              call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2,stat) ! | I_2+> in I2.
            endif
            if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2) ! | I_2+> in I2.

            ! Acquire non-derivative T | R+_1> via | R+_1>.
            ! The intermediate result | R+_1> is not needed afterwards.
            ! These do not depend on the instrument and will be used in the end for the
            ! derivatives.
            call source_response_1_sp(nstrhalf,nlay,nsplit_total,.false.,.true.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%i3) ! R+_1> in I3.
            call transmit_sp(nstrhalf,nlay,nsplit_total,.true.,self%fld%i3,self%tlc,self%fld%i4) ! T | R+_1> in I4.

            ! Okay, this was all the geometry-independent stuff for multi-scattering. Now
            ! the loop over geometries can be opened.

          endif

        else

          ! Multi-scattering.
          if (execution(3)) then

            ! L | I_1+> = | Q_1>
            ! E_2+ = <R_1 | I_1+>
            call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.false.,exist_bdrf,exist_emi,exist_therm,self%tlc,self%fld%i1) ! | Q_1> in I1.

            ! Solve the matrix.
            if (set%solver .eq. solver_gauss_seidel) then
              call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i1,self%tlc,self%fld%i2) ! T | Q_1> in I2 as first-guess for | I_1+>.
              call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2,stat) ! | I_1+> in I2.
            endif
            if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i1,self%fld%i2) ! | I_1+> in I2.

          endif

        endif

        ! Now, the geometry-dependent stuff is going to happen. Both the double-scattering
        ! and the multi-scattering will be resumed after the geometry-dependent
        ! TLC-parameters have been calculated for this geometry. Note that these
        ! are all Fourier-dependent TLC-parameters, so they are used only one iteration
        ! of this loop.
        do igeo = 1,ngeo

          ! Skip if this geometry is already converged.
          if (fourier_converged(igeo)) cycle
          ! Fourier convergence can only be achieved after doing the second Fourier component.
          ! This is because for Q and U, the second Fourier component is usually much more
          ! significant than the first, so convergence should not be concluded after evaluating
          ! the first Fourier component.
          if (m .ge. 2) fourier_converged(igeo) = .true. ! Will be set to false at the end of loop when not yet converged.

          ! Also set the pointers to Fourier-independent instrument-dependent
          ! TLC-parameters, because a new geometry is active now.
          ! The boolean is the flag for derivatives.
          call choose_viewing_geometry_for_fourier_independent_tlc_parameters(nlay,.true.,igeo,self%geo,self%tlc)

          do ist = 1,nst

            ! Skip Stokes parameter if not selected
            if (.not. set%execute_stokes(ist)) cycle

            ! Skip Fourier number zero for Stokes parameters on the other half as the sun.
            ! They have zero degeneracy, but may result is zero-division in Gauss-Seidel
            ! and they cost time, for just a zero.
            if (m .eq. 0 .and. ((ist_sun .ge. ist_u) .neqv. (ist .ge. ist_u))) cycle

            ! Initialize results to zero.
            rint_cur = 0.0

            if (set%differentiate_stokes(ist)) then

              drint_tau_cur = 0.0
              drint_coefs_cur = 0.0
              drint_bdrf_cur => drv%drint_bdrf(:,:,:,:,m,ist,igeo)
              drint_bdrf_0_cur => drv%drint_bdrf_0(ist_sun,:,:,m,ist,igeo)
              drint_bdrf_v_cur => drv%drint_bdrf_v(:,ist,:,m,ist,igeo)
              drint_emi_cur => drv%drint_emi(:,:,m,ist,igeo)

              ! Initialize intermediate C-derivatives to zero.
              drint_c_fwd = 0.0
              drint_c_bck = 0.0
              drint_c_fwd_0 = 0.0
              drint_c_bck_0 = 0.0
              drint_c_fwd_v = 0.0
              drint_c_bck_v = 0.0

              ! Initialize intermediate C-parameter to zero, noting that there is also
              ! a single-scattering soluation.
              drint_c_therm_cur = 0.0

            endif

            ! Calculate the Fourier-dependent TLC-parameters in which the viewing
            ! geometry is involved.
            call calculate_fourier_dependent_tlc_parameters_for_viewing_geometry(nst,nstrhalf,nleg,nlay,igeo,m,ist,exist_bdrf,self%grd,self%geo,self%pix,self%tlc)

            if (set%split_double) then

              ! Resume double-scattering.
              if (execution(2)) then

                ! Overview of existing (instrument-independent) fields:
                ! | Q_1> in A1.
                ! T | Q_1> in A2.

                ! Complete forward run.
                call direct_2(nstrhalf,nlay,exist_therm,self%tlc,rint_cur)
                call source_response_1(nstrhalf,nlay,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%a1) ! <R_1 | in A1.
                call innerproduct(nstrhalf,nlay,self%fld%a1,self%fld%a2,rint_cur)

                ! Derivatives.

                if (set%differentiate_stokes(ist)) then

                  ! Case 1: Direct.
                  if (deriv_c) call perturbation_scat_direct_2(nstrhalf,nlay,nlay_deriv_c,ilay_deriv_c,exist_therm,self%tlc,drint_c_fwd_0,drint_c_bck_0,drint_c_fwd_v,drint_c_bck_v)
                  if (deriv_base) call perturbation_tau_direct_2(nstrhalf,nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),exist_therm,self%tlc,drint_tau_cur)
                  if (exist_therm .and. deriv_base) call perturbation_therm_direct_2(nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%tlc,drint_c_therm_cur)

                  ! Case 2: <R'_1 | T | Q_1> (using A2).
                  if (deriv_c) call perturbation_scat_source_response_1(nstrhalf,nlay,nlay_deriv_c,ilay_deriv_c,.false.,.false.,self%fld%a2,self%tlc,drint_c_fwd_v,drint_c_bck_v)
                  if (deriv_base) call perturbation_tau_source_response_1(nstrhalf,nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.false.,exist_bdrf,.false.,self%fld%a2,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nlay,.false.,.false.,self%fld%a2,self%tlc,drint_bdrf_v_cur)

                  ! Calculate non-derivative.
                  call transmit(nstrhalf,nlay,.true.,self%fld%a1,self%tlc,self%fld%a3) ! <R_1 | T in A3.

                  ! Case 3: <R_1 | T' | Q_1> (only for tau), using non-derivatives A2 and A3.
                  if (deriv_base) call perturbation_tau_mutualtransmit(nstrhalf,nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%a2,self%fld%a3,self%tlc,drint_tau_cur)

                  ! Case 4: <R_1 | T | Q'_1> (using A3).
                  if (deriv_c) call perturbation_scat_source_response_1(nstrhalf,nlay,nlay_deriv_c,ilay_deriv_c,.true.,.false.,self%fld%a3,self%tlc,drint_c_fwd_0,drint_c_bck_0)
                  if (deriv_base) call perturbation_tau_source_response_1(nstrhalf,nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.true.,.false.,exist_bdrf,exist_therm,self%fld%a3,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nlay,.true.,.false.,self%fld%a3,self%tlc,drint_bdrf_0_cur)
                  if (exist_emi) call perturbation_emi_source_response_1(nstrhalf,nlay,.true.,self%fld%a3,self%tlc,drint_emi_cur)
                  if (exist_therm .and. deriv_base) call perturbation_therm_source_response_1(nstrhalf,nlay,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.true.,self%fld%a3,self%tlc,drint_c_therm_cur)

                endif

              endif

              ! Resume multi-scattering.
              if (execution(3)) then

                ! Overview of existing (viewing-Stokes-parameter-independent) fields:
                ! | Q_2> in I1 (though not necessary).
                ! | I_2+> in I2.
                ! | R+_1> in I3 (though not necessary).
                ! T | R+_1> in I4.

                ! Complete forward run, response function, though not needed later on
                call source_response_1_int(nstrhalf,nlay,nsplit_total,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i3) ! <R_1 | in I3
                call innerproduct(nstrhalf,nsplit_total,self%fld%i3,self%fld%i2,rint_cur)

                ! Derivatives.
                if (set%differentiate_stokes(ist)) then

                  ! Case 1: <R_1' | I_2+> (using I2).
                  if (deriv_c) call perturbation_scat_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,.false.,.false.,self%fld%i2,self%tlc,drint_c_fwd_v,drint_c_bck_v)
                  if (deriv_base) call perturbation_tau_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.false.,exist_bdrf,.false.,self%fld%i2,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nsplit_total,.false.,.false.,self%fld%i2,self%tlc,drint_bdrf_v_cur)

                  ! Adjoint.
                  call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.true.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i3) ! <Q+_1 | in I3.

                  ! Solve the matrix.
                  if (set%solver .eq. solver_gauss_seidel) then
                    call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i3,self%tlc,self%fld%i1) ! <Q+_1 | T in I1 as first-guess for <I+_1+ |.
                    call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i3,self%fld%i1,stat) ! <I+_1+ | in I1.
                  endif
                  if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i3,self%fld%i1) ! <I+_1+ | in I1.

                  ! Case 2: <I+_1+ | (-DL) | I_2> (using I1 and I2).
                  if (deriv_c) call perturbation_scat_zip(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i2,self%fld%i1,self%tlc,self%pix%chain_zipscat_c,drint_c_fwd,drint_c_bck,self%pix%chain_zipscat_tau,drint_tau_cur)
                  if (deriv_base) call perturbation_abs_zip(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i2,self%fld%i1,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_zip(nstrhalf,nsplit_total,self%fld%i2,self%fld%i1,self%tlc,drint_bdrf_cur)

                  ! Case 3: <I+_1+ | R+'_2*> (using I1).
                  if (deriv_c) call perturbation_scat_source_response_2_int(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,.false.,.true.,exist_therm,self%fld%i1,self%tlc,drint_c_fwd,drint_c_bck,drint_c_fwd_0,drint_c_bck_0)
                  if (deriv_base) call perturbation_tau_source_response_2_int(nstrhalf,nlay,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.true.,exist_therm,self%fld%i1,self%tlc,drint_tau_cur)
                  if (exist_therm .and. deriv_base) call perturbation_therm_source_response_2_int(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,self%fld%i1,self%tlc,drint_c_therm_cur)

                  ! Case 4: <I+_1+ | S' T | R+_1>.
                  ! Using I1 and I4, the latter is a non-derivative: T | R+_1>.
                  ! Hereafter, <I_1+ | (I1) can be recycled.
                  if (deriv_c) call perturbation_scat_scatter_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,self%fld%i1,self%fld%i4,self%tlc,drint_c_fwd,drint_c_bck)
                  if (deriv_base) call perturbation_tau_scatter_1_int(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i1,self%fld%i4,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_scatter_1(nstrhalf,nsplit_total,self%fld%i1,self%fld%i4,self%tlc,drint_bdrf_cur)

                  ! Transform <I+_1+ |, which is no longer needed, into <I+_1+ | S_1 T
                  ! The unused field I3 will be used for the intermediate result <I+_1+ | S_1
                  ! after which I1 will be overwritten with the final result.
                  call scatter_1_int(nstrhalf,nlay,nsplit_total,exist_bdrf,self%fld%i1,self%tlc,self%fld%i3) ! <I+_1+ | S_1 in I3.
                  call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i3,self%tlc,self%fld%i1) ! <I+_1+ | S_1 T back in I1.

                  ! Case 5: <I+_1+ | S_1 T | R'+_1> (using I1).
                  if (deriv_c) call perturbation_scat_source_response_1_sp(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,.false.,.true.,self%fld%i1,self%tlc,drint_c_fwd_0,drint_c_bck_0)
                  if (deriv_base) call perturbation_tau_source_response_1_sp(nstrhalf,nlay,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.true.,exist_bdrf,exist_therm,self%fld%i1,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nsplit_total,.false.,.true.,self%fld%i1,self%tlc,drint_bdrf_0_cur)
                  if (exist_emi) call perturbation_emi_source_response_1(nstrhalf,nsplit_total,.false.,self%fld%i1,self%tlc,drint_emi_cur)
                  if (exist_therm .and. deriv_base) call perturbation_therm_source_response_1_sp(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,self%fld%i1,self%tlc,drint_c_therm_cur)

                  ! Case 6: <I+_1+ | S T' | R+_1> (using I4 and I1)
                  if (deriv_base) call perturbation_tau_mutualtransmit_sp(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i1,self%fld%i4,self%tlc,drint_tau_cur)

                endif

              endif

            else

              if (execution(3)) then

                ! Finish forward calculation.
                call source_response_1_int(nstrhalf,nlay,nsplit_total,.false.,.false.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i3) ! <R_1 | in I3.
                call innerproduct(nstrhalf,nsplit_total,self%fld%i3,self%fld%i2,rint_cur)

                ! Derivatives.
                if (set%differentiate_stokes(ist)) then

                  ! Case 1: <R_1' | I_1+> (using I2).
                  if (deriv_c) call perturbation_scat_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,.false.,.false.,self%fld%i2,self%tlc,drint_c_fwd_v,drint_c_bck_v)
                  if (deriv_base) call perturbation_tau_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.false.,exist_bdrf,.false.,self%fld%i2,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nsplit_total,.false.,.false.,self%fld%i2,self%tlc,drint_bdrf_v_cur)

                  ! Adjoint.
                  call source_response_1_sp(nstrhalf,nlay,nsplit_total,.true.,.true.,exist_bdrf,.false.,.false.,self%tlc,self%fld%i3) ! <Q+_1 | in I3.

                  ! Solve the matrix.
                  if (set%solver .eq. solver_gauss_seidel) then
                    call transmit_sp(nstrhalf,nlay,nsplit_total,.false.,self%fld%i3,self%tlc,self%fld%i1) ! <Q+_1 | T in I1 as first-guess for <I+_1+ |.
                    call gauss_seidel(nst,nstrhalf,nsplit_total,nmat,nband,set%gs_maxiter,set%gs_tolerance,self%fld%mat,self%fld%shifts,self%fld%i3,self%fld%i1,stat) ! <I+_1+ | in I1.
                  endif
                  if (set%solver .eq. solver_lu_decomposition) call lu_solve(nmat,nband,self%fld%mat,self%fld%shifts,self%fld%i3,self%fld%i1) ! <I+_1+ | in I1.

                  ! Case 2: <I+_1+ | (-DL) | I_1> (using I2 and I1).
                  if (deriv_c) call perturbation_scat_zip(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i2,self%fld%i1,self%tlc,self%pix%chain_zipscat_c,drint_c_fwd,drint_c_bck,self%pix%chain_zipscat_tau,drint_tau_cur)
                  if (deriv_base) call perturbation_abs_zip(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),self%fld%i2,self%fld%i1,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_zip(nstrhalf,nsplit_total,self%fld%i2,self%fld%i1,self%tlc,drint_bdrf_cur)

                  ! Case 3: <I+_1+ | R+'_1> (using I1).
                  if (deriv_c) call perturbation_scat_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_c,ilay_deriv_c,.false.,.true.,self%fld%i1,self%tlc,drint_c_fwd_0,drint_c_bck_0)
                  if (deriv_base) call perturbation_tau_source_response_1_int(nstrhalf,nlay,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,.true.,exist_bdrf,exist_therm,self%fld%i1,self%tlc,drint_tau_cur)
                  if (exist_bdrf) call perturbation_bdrf_source_response_1(nstrhalf,nlay,nsplit_total,.false.,.true.,self%fld%i1,self%tlc,drint_bdrf_0_cur)
                  if (exist_emi) call perturbation_emi_source_response_1(nstrhalf,nsplit_total,.false.,self%fld%i1,self%tlc,drint_emi_cur)
                  if (exist_therm .and. deriv_base) call perturbation_therm_source_response_1_int(nstrhalf,nsplit_total,nlay_deriv_base,set%ilay_deriv_base(1:nlay_deriv_base),.false.,self%fld%i1,self%tlc,drint_c_therm_cur)

                endif

              endif

            endif

            ! We are still in the loop over observed Stokes parameters. Derivatives only
            ! need to be prgoressed when calculating derivatives.

            ! First apply the degeneracy factor on the radiance itself. Then, we can open
            ! an if-clause for calculating the derivative.
            rint_fourier(ist,igeo) = rint_fourier(ist,igeo) + rint_cur * self%geo%degenerate_allgeo(m,ist,igeo)

            if (set%differentiate_stokes(ist)) then

              ! Convert nst by nst derivatives to derivative with respect to coefficients. {{{
              ! Single-scattering albedo.
              if (deriv_base) then
                do ideriv = 1,nlay_deriv_base
                  ilay = set%ilay_deriv_base(ideriv)
                  drint_ssa_cur(ideriv) = sum(drint_c_fwd_0(:,:,ilay) * self%tlc%diff_c_fwd_0_nrm_ssa(:,:,ilay)) + sum(drint_c_bck_0(:,:,ilay) * self%tlc%diff_c_bck_0_nrm_ssa(:,:,ilay)) + sum(drint_c_fwd_v(:,:,ilay) * self%tlc%diff_c_fwd_v_nrm_ssa(:,:,ilay)) + sum(drint_c_bck_v(:,:,ilay) * self%tlc%diff_c_bck_v_nrm_ssa(:,:,ilay)) + sum(drint_c_fwd(:,:,:,:,ilay) * self%tlc%diff_c_fwd_ssa(:,:,:,:,ilay)) + sum(drint_c_bck(:,:,:,:,ilay) * self%tlc%diff_c_bck_ssa(:,:,:,:,ilay))
                enddo
              endif
              ! Phase coefficients.
              if (deriv_ph) then
                do ideriv = 1,nlay_deriv_ph
                  ilay = set%ilay_deriv_ph(ideriv)
                  nleg_pert_lay = min(nleg_pert,self%pix%nleg_lay(ilay))

                  ! Initialize symmetry conditions. They change with ileg.
                  sym_even = 1.0
                  sym_odd = -1.0

                  ! Loop over Legendre moments and streams.
                  ! Array operation over target Stokes parameters and perturbed layers.
                  ! Individual expressions for matrix elements (nscat and nst by nst).

                  do ileg = m,nleg_pert_lay ! Legendre moment
                    do istr = 1,nstrhalf ! Stream

                      ! For derivatives with respect to C-parameters for the sun or the instrument,
                      ! it depends on which Stokes parameters are involved. That applies both
                      ! for the sun and for the instrument.
                      ! Sun
                      if (ist_sun .eq. ist_i) then
                        drint_coefs_cur(iscat_a1,ileg,ideriv) = drint_coefs_cur(iscat_a1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_i,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_0(ist_i,istr,ilay) * sym_even * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        if (pol_nst(nst)) then
                          drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_even * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_even * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_m(istr,ileg,m))
                        endif
                      endif
                      if (ist_sun .eq. ist_q) then
                        ! With a sun shining Q, we assume at least polarization.
                        drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_i,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_0(ist_i,istr,ilay) * sym_even * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        drint_coefs_cur(iscat_a2,ileg,ideriv) = drint_coefs_cur(iscat_a2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_even * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_even * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_m(istr,ileg,m))
                        drint_coefs_cur(iscat_a3,ileg,ideriv) = drint_coefs_cur(iscat_a3,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_odd * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_odd * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_p(istr,ileg,m))
                        if (fullpol_nst(nst)) then
                          drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (-drint_c_fwd_0(ist_v,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) - drint_c_bck_0(ist_v,istr,ilay) * sym_odd * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        endif
                      endif
                      if (ist_sun .eq. ist_u) then
                        ! With a sun shining U, we assume at least polarization.
                        drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_i,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_0(ist_i,istr,ilay) * sym_even * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        drint_coefs_cur(iscat_a2,ileg,ideriv) = drint_coefs_cur(iscat_a2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_even * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_even * self%geo%gsf_m_0(ileg,m) * self%grd%gsf_m(istr,ileg,m))
                        drint_coefs_cur(iscat_a3,ileg,ideriv) = drint_coefs_cur(iscat_a3,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_odd * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_odd * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_p(istr,ileg,m))
                        if (fullpol_nst(nst)) then
                          drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (-drint_c_fwd_0(ist_v,istr,ilay) * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) - drint_c_bck_0(ist_v,istr,ilay) * sym_odd * self%geo%gsf_p_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        endif
                      endif
                      if (ist_sun .eq. ist_v) then
                        ! Now, we can assume full polarization.
                        drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_q,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd_0(ist_u,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_0(ist_q,istr,ilay) * sym_odd * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_0(ist_u,istr,ilay) * sym_odd * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_p(istr,ileg,m))
                        drint_coefs_cur(iscat_a4,ileg,ideriv) = drint_coefs_cur(iscat_a4,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_0(ist_v,istr,ilay) * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_0(ist_v,istr,ilay) * sym_odd * self%geo%gsf_0_0(ileg,m) * self%grd%gsf_0(istr,ileg,m))
                      endif

                      ! Instrument.
                      if (ist .eq. ist_i) then
                        drint_coefs_cur(iscat_a1,ileg,ideriv) = drint_coefs_cur(iscat_a1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_i,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_i,istr,ilay) * sym_even * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                        if (pol_nst(nst)) then
                          drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_v(ist_q,istr,ilay) * sym_even * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_v(ist_u,istr,ilay) * sym_even * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m))
                        endif
                      endif
                      if (ist .eq. ist_q) then
                        ! With an instrument measuring Q, we assume at least polarization.
                        drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_i,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_i,istr,ilay) * sym_even * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                        drint_coefs_cur(iscat_a2,ileg,ideriv) = drint_coefs_cur(iscat_a2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_v(ist_q,istr,ilay) * sym_even * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_v(ist_u,istr,ilay) * sym_even * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m))
                        drint_coefs_cur(iscat_a3,ileg,ideriv) = drint_coefs_cur(iscat_a3,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_v(ist_q,istr,ilay) * sym_odd * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_v(ist_u,istr,ilay) * sym_odd * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m))
                        if (fullpol_nst(nst)) then
                          drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_v,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_v,istr,ilay) * sym_odd * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                        endif
                      endif
                      if (ist .eq. ist_u) then
                        ! With an instrument measuring U, we assume at least polarization.
                        drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_i,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_i,istr,ilay) * sym_even * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                        drint_coefs_cur(iscat_a2,ileg,ideriv) = drint_coefs_cur(iscat_a2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_v(ist_q,istr,ilay) * sym_even * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_v(ist_u,istr,ilay) * sym_even * self%geo%gsf_m_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m))
                        drint_coefs_cur(iscat_a3,ileg,ideriv) = drint_coefs_cur(iscat_a3,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck_v(ist_q,istr,ilay) * sym_odd * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck_v(ist_u,istr,ilay) * sym_odd * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m))
                        if (fullpol_nst(nst)) then
                          drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_v,istr,ilay) * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_v,istr,ilay) * sym_odd * self%geo%gsf_p_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                        endif
                      endif
                      if (ist .eq. ist_v) then
                        ! Now, we can assume full polarization.
                        drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (-drint_c_fwd_v(ist_q,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) - drint_c_fwd_v(ist_u,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m) - drint_c_bck_v(ist_q,istr,ilay) * sym_odd * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_m(istr,ileg,m) - drint_c_bck_v(ist_u,istr,ilay) * sym_odd * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_p(istr,ileg,m))
                        drint_coefs_cur(iscat_a4,ileg,ideriv) = drint_coefs_cur(iscat_a4,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd_v(ist_v,istr,ilay) * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck_v(ist_v,istr,ilay) * sym_odd * self%geo%gsf_0_v_allgeo(ileg,m,igeo) * self%grd%gsf_0(istr,ileg,m))
                      endif

                      ! Add contributions for phase functions from one internal stream to
                      ! another internal stream. These are relevant for all sun and instrument
                      ! related Stokes parameters.
                      do istr_src = 1,nstrhalf
                        drint_coefs_cur(iscat_a1,ileg,ideriv) = drint_coefs_cur(iscat_a1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd(ist_i,ist_i,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck(ist_i,ist_i,istr_src,istr,ilay) * sym_even * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        if (pol_nst(nst)) then
                          drint_coefs_cur(iscat_b1,ileg,ideriv) = drint_coefs_cur(iscat_b1,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd(ist_i,ist_q,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd(ist_i,ist_u,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd(ist_q,ist_i,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_fwd(ist_u,ist_i,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck(ist_i,ist_q,istr_src,istr,ilay) * sym_even * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck(ist_i,ist_u,istr_src,istr,ilay) * sym_even * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_q,ist_i,istr_src,istr,ilay) * sym_even * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck(ist_u,ist_i,istr_src,istr,ilay) * sym_even * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m))
                          drint_coefs_cur(iscat_a2,ileg,ideriv) = drint_coefs_cur(iscat_a2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd(ist_q,ist_q,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd(ist_q,ist_u,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd(ist_u,ist_q,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd(ist_u,ist_u,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_q,ist_q,istr_src,istr,ilay) * sym_even * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck(ist_q,ist_u,istr_src,istr,ilay) * sym_even * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_u,ist_q,istr_src,istr,ilay) * sym_even * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck(ist_u,ist_u,istr_src,istr,ilay) * sym_even * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m))
                          drint_coefs_cur(iscat_a3,ileg,ideriv) = drint_coefs_cur(iscat_a3,ileg,ideriv) + self%tlc%diff_c_ph(ilay) *(drint_c_fwd(ist_q,ist_q,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd(ist_q,ist_u,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_fwd(ist_u,ist_q,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd(ist_u,ist_u,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck(ist_q,ist_q,istr_src,istr,ilay) * sym_odd * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_q,ist_u,istr_src,istr,ilay) * sym_odd * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) + drint_c_bck(ist_u,ist_q,istr_src,istr,ilay) * sym_odd * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_u,ist_u,istr_src,istr,ilay) * sym_odd * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m))
                        endif
                        if (fullpol_nst(nst)) then
                          drint_coefs_cur(iscat_b2,ileg,ideriv) = drint_coefs_cur(iscat_b2,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (-drint_c_fwd(ist_q,ist_v,istr_src,istr,ilay) * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) - drint_c_fwd(ist_u,ist_v,istr_src,istr,ilay) * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_fwd(ist_v,ist_q,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_fwd(ist_v,ist_u,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m) - drint_c_bck(ist_q,ist_v,istr_src,istr,ilay) * sym_odd * self%grd%gsf_m(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) - drint_c_bck(ist_u,ist_v,istr_src,istr,ilay) * sym_odd * self%grd%gsf_p(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck(ist_v,ist_q,istr_src,istr,ilay) * sym_odd * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_m(istr,ileg,m) + drint_c_bck(ist_v,ist_u,istr_src,istr,ilay) * sym_odd * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_p(istr,ileg,m))
                          drint_coefs_cur(iscat_a4,ileg,ideriv) = drint_coefs_cur(iscat_a4,ileg,ideriv) + self%tlc%diff_c_ph(ilay) * (drint_c_fwd(ist_v,ist_v,istr_src,istr,ilay) * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m) + drint_c_bck(ist_v,ist_v,istr_src,istr,ilay) * sym_odd * self%grd%gsf_0(istr_src,ileg,m) * self%grd%gsf_0(istr,ileg,m))
                        endif
                      enddo
                    enddo

                    ! Switch even and odd.
                    sym_even = -sym_even
                    sym_odd = -sym_odd

                  enddo
                enddo
              endif

              ! }}}

              ! Count derivatives of current Fourier number.
              ! Derivatives with respect to ssa, tau and coefficients are subject to
              ! the chain rule for delta-M. Therefore, the `pre-dm' suffix is added.
              ! The derivative with respect to the phase function at single-scattering
              ! geometry is not there, because it is not involved in the Fourier loop.

              if (deriv_base) then
                ! Array operation over altitudes.
                drint_ssa_predm(:,ist,igeo) = drint_ssa_predm(:,ist,igeo) + drint_ssa_cur * self%geo%degenerate_allgeo(m,ist,igeo)
                drint_tau_predm(:,ist,igeo) = drint_tau_predm(:,ist,igeo) + drint_tau_cur * self%geo%degenerate_allgeo(m,ist,igeo)
              endif

              ! Repeat joke for coefficients.
              if (deriv_ph) then
                do ideriv = 1,nlay_deriv_ph
                  nleg_pert_lay = min(nleg_pert,self%pix%nleg_lay(set%ilay_deriv_ph(ideriv)))
                  drint_coefs_predm(:,m:nleg_pert_lay,ideriv,ist,igeo) = drint_coefs_predm(:,m:nleg_pert_lay,ideriv,ist,igeo) + drint_coefs_cur(:,m:nleg_pert_lay,ideriv) * self%geo%degenerate_allgeo(m,ist,igeo)
                enddo
              endif

              ! The BDRF is defined per Fourier mode, so it is not summed up. The
              ! pointers, however, may point to something that is not allocated in the
              ! case that there is no contribution. Therefore, the execution if-clauses
              ! are applied.
              if (execution(3) .and. exist_bdrf) drint_bdrf_cur = drint_bdrf_cur * self%geo%degenerate_allgeo(m,ist,igeo) ! It is a pointer.

              ! We are already in an if-execution-2-or-3 clause.
              ! Surface reflection.
              if (exist_bdrf) then
                drint_bdrf_0_cur = drint_bdrf_0_cur * self%geo%degenerate_allgeo(m,ist,igeo) ! It is a pointer.
                drint_bdrf_v_cur = drint_bdrf_v_cur * self%geo%degenerate_allgeo(m,ist,igeo) ! It is a pointer.
              endif

              ! Emissivity.
              if (exist_emi) drint_emi_cur = drint_emi_cur * self%geo%degenerate_allgeo(m,ist,igeo) ! It is a pointer.

              ! Thermal emissivity.
              if (exist_therm) drint_c_therm(:,ist,igeo) = drint_c_therm(:,ist,igeo) + drint_c_therm_cur * self%geo%degenerate_allgeo(m,ist,igeo)

            endif

            ! Fourier convergence check. Uses the conventions of Lintran 1, except for the
            ! degeneracy factor and the Nakajima single-scattering difference.
            ! This implies that convergence is only checked for Stokes parameter I. the reference
            ! radiance is just the Fourier-part of the up-to-now calculated radiance. In contrast
            ! to Lintran 1, the Fourier part does not include some single-scattering difference
            ! from the Nakajima wrapper. There is always a calculation up to and including Fourier
            ! number 2, because that is often a very important one for polarization.
            if (ist .eq. ist_i .and. fourier_converged(igeo) .and. abs(rint_fourier(ist,igeo)) .gt. 0.0) then ! Against zero-division.
              if (abs(rint_cur / rint_fourier(ist,igeo)) .ge. set%fourier_tolerance) fourier_converged(igeo) = .false.
            endif

          enddo

        enddo

        ! Quit Fourier loop if all geometries are converged.
        if (all(fourier_converged)) exit

      enddo

      ! Add up the Fourier part to the eventual single-scattering part.
      ! This is not needed for the derivatives, which can be stored directly into their
      ! end results. This is only to have the Fourier convergence criterion to be on
      ! the Fourier-only part of the radiance. The derivatives play no role in the Fourier
      ! convergence.
      rint = rint + rint_fourier

    endif

    ! Convert derivative with respect to thermal emission to derivatives with respect to
    ! single-scattering albedo and the Planck curve. This is done before handling delta-M.
    ! The thermal C-parameter was calculated after delta-M, so the calculation is done in the
    ! delta-M-corrected world. Note that we have to use the old thermal-emission flag again,
    ! because exist_therm is destroyed in the Fourier loop.
    if (self%pix%thermal_emission .and. deriv_base) then
      do igeo = 1,ngeo
        do ist = 1,nst
          if (set%execute_stokes(ist) .and. set%differentiate_stokes(ist)) then
            drint_ssa_predm(:,ist,igeo) = drint_ssa_predm(:,ist,igeo) + drint_c_therm(:,ist,igeo) * self%tlc%diff_c_therm_ssa(set%ilay_deriv_base(1:nlay_deriv_base))
            drv%drint_planck_curve(1:nlay_deriv_base,ist,igeo) = drint_c_therm(:,ist,igeo) * self%tlc%diff_c_therm_planck(set%ilay_deriv_base(1:nlay_deriv_base))
          endif
        enddo
      enddo
    endif

    ! Apply delta-M chain rules.
    if (set%deltam) then
      if (deriv_base) then
        do igeo = 1,ngeo
          do ist = 1,nst
            if (set%execute_stokes(ist) .and. set%differentiate_stokes(ist)) then
              drint_ssa(:,ist,igeo) = self%pix%chain_ssa_ssa(set%ilay_deriv_base(1:nlay_deriv_base)) * drint_ssa_predm(:,ist,igeo) + self%pix%chain_ssa_tau(set%ilay_deriv_base(1:nlay_deriv_base)) * drint_tau_predm(:,ist,igeo)
              drint_tau(:,ist,igeo) = self%pix%chain_tau_tau(set%ilay_deriv_base(1:nlay_deriv_base)) * drint_tau_predm(:,ist,igeo)
            else
              ! Set non-executed derivative to zero, because they will be progressed to the end
              ! results in an array operation regardless of whether the Stokes parameter was
              ! executed or not.
              drint_ssa(:,ist,igeo) = 0.0
              drint_tau(:,ist,igeo) = 0.0
            endif
          enddo
        enddo
      endif
      if (deriv_ph) then
        if (execution(2) .or. execution(3)) then
          do igeo = 1,ngeo
            do ist = 1,nst
              if (set%execute_stokes(ist) .and. set%differentiate_stokes(ist)) then
                do ideriv = 1,nlay_deriv_ph
                  nleg_pert_lay = min(nleg_pert,self%pix%nleg_lay(set%ilay_deriv_ph(ideriv)))
                  drv%drint_coefs(1:nscat,0:nleg_pert_lay,ideriv,ist,igeo) = self%pix%chain_ph_ph(set%ilay_deriv_ph(ideriv)) * drint_coefs_predm(:,0:nleg_pert_lay,ideriv,ist,igeo)
                enddo
              endif
            enddo
          enddo
        endif
        if (execution(1)) then
          do igeo = 1,ngeo
            do ist = 1,nst
              if (set%execute_stokes(ist) .and. set%differentiate_stokes(ist)) then
                do ideriv = 1,nlay_deriv_ph
                  drv%drint_phase_ssg(1:nscat,ideriv,ist,igeo) = self%pix%chain_ph_ph(set%ilay_deriv_ph(ideriv)) * drint_phase_ssg_predm(:,ideriv,ist,igeo)
                enddo
              endif
            enddo
          enddo
        endif
      endif

      ! The delta-M phase coefficient. Only the layers are written where all derivatives
      ! are calculated. Otherwise, this value remains zero.
      if (deriv_base .and. deriv_ph .and. set%nleg_deriv .gt. nleg) then
        do ideriv = 1,nlay_deriv_base
          ilay = set%ilay_deriv_base(ideriv)
          if (atm%nleg_lay(ilay) .gt. nleg) then
            do igeo = 1,ngeo
              do ist = 1,nst
                if (set%execute_stokes(ist) .and. set%differentiate_stokes(ist)) then
                  where (set%ilay_deriv_ph(1:nlay_deriv_ph) .eq. ilay) ! This is at most one element.
                    drv%drint_coefs(iscat_a1,nleg+1,1:nlay_deriv_ph,ist,igeo) = self%pix%chain_f_ssa(ilay) * drint_ssa_predm(ideriv,ist,igeo) + self%pix%chain_f_tau(ilay) * drint_tau_predm(ideriv,ist,igeo) + matmul(reshape(self%pix%chain_f_ph(:,:,ilay),(/(nleg+1)*nscat/)) , reshape(drint_coefs_predm(:,:,:,ist,igeo),(/(nleg+1)*nscat,nlay_deriv_ph/))) + matmul(self%pix%chain_f_ph_elem(:,ilay,igeo) , drint_phase_ssg_predm(:,:,ist,igeo))
                  endwhere
                endif
              enddo
            enddo
          endif
        enddo
      endif
    else
      ! No delta-M. Just rename the derivatives.
      drint_ssa = drint_ssa_predm
      drint_tau = drint_tau_predm
      if (execution(2) .or. execution(3)) then
        do ideriv = 1,nlay_deriv_ph
          nleg_pert_lay = min(nleg_pert,self%pix%nleg_lay(set%ilay_deriv_ph(ideriv)))
          drv%drint_coefs(1:nscat,0:nleg_pert_lay,ideriv,1:nst,1:ngeo) = drint_coefs_predm(:,0:nleg_pert_lay,ideriv,:,:)
        enddo
      endif
      if (execution(1)) drv%drint_phase_ssg(1:nscat,1:nlay_deriv_ph,1:nst,1:ngeo) = drint_phase_ssg_predm
    endif

    ! From here onwards, no execute Stokes or differentiate Stokes check is performed, because
    ! the observed Stokes parameter is taken as array operation.

    ! The interpretation is still as if the BDRF is Delta-flipped for
    ! U and V as destination directions, so that has to be undone. Not
    ! for single-scattering, since it does not use Fourier.
    if (pol_nst(nst)) then
      if (execution(2) .or. execution(3)) then
        drv%drint_bdrf_v(:,ist_u:nst,:,:,:,:) = -drv%drint_bdrf_v(:,ist_u:nst,:,:,:,:)
        drv%drint_bdrf_0(:,ist_u:nst,:,:,:,:) = -drv%drint_bdrf_0(:,ist_u:nst,:,:,:,:)

        ! Same joke for fluorescent emission. It is Delta-flipped for the
        ! Fourier cases, because it is upward radiation.
        drv%drint_emi(ist_u:nst,:,:,:,:) = -drv%drint_emi(ist_u:nst,:,:,:,:)
      endif

      ! This BDRF is only for execution 3.
      if (execution(3)) drv%drint_bdrf(:,ist_u:nst,:,:,:,:,:) = -drv%drint_bdrf(:,ist_u:nst,:,:,:,:,:)
    endif

    ! Transpose the derivative with respect to the BDRFs, so that the
    ! destination direction is the first dimension again. This is necessary
    ! because the BDRF is not perfectly symmetric.
    if (execution(1)) drv%drint_bdrf_ssg(1:nst,1:nst,1:nst,1:ngeo) = reshape(drv%drint_bdrf_ssg(1:nst,1:nst,1:nst,1:ngeo) , (/nst,nst,nst,ngeo/), order = (/2,1,3,4/))
    if (execution(2) .or. execution(3)) then
      drv%drint_bdrf_v(1:nst,1:nst,1:nstrhalf,0:nleg,1:nst,1:ngeo) = reshape(drv%drint_bdrf_v(1:nst,1:nst,1:nstrhalf,0:nleg,1:nst,1:ngeo) , (/nst,nst,nstrhalf,nleg+1,nst,ngeo/), order = (/2,1,3,4,5,6/))
      drv%drint_bdrf_0(1:nst,1:nst,1:nstrhalf,0:nleg,1:nst,1:ngeo) = reshape(drv%drint_bdrf_0(1:nst,1:nst,1:nstrhalf,0:nleg,1:nst,1:ngeo) , (/nst,nst,nstrhalf,nleg+1,nst,ngeo/), order = (/2,1,3,4,5,6/))
    endif
    if (execution(3)) drv%drint_bdrf(1:nst,1:nst,1:nstrhalf,1:nstrhalf,0:nleg,1:nst,1:ngeo) = reshape(drv%drint_bdrf(1:nst,1:nst,1:nstrhalf,1:nstrhalf,0:nleg,1:nst,1:ngeo) , (/nst,nst,nstrhalf,nstrhalf,nleg+1,nst,ngeo/), order = (/2,1,4,3,5,6,7/))

    ! There are bound elements in the BDRF, each with an individual
    ! derivative, but the derivatives are calculated assuming the
    ! correct symmetry relationship. So, those individual components
    ! make no sense, just the derivatives of the independent matrix
    ! elements. So, we will write a 0 if the first stream index is
    ! higher than the second. And for the same streams, we do the
    ! same with Stokes parameters. For the other elements, we write
    ! the correct aggregate.

    ! Both streams take the role of source and destination once.
    ! Moreover, we just transposed the matrix, so the meaning of
    ! source and destination is confusing. Now, it does not matter.
    if (execution(3)) then
      do istr = 1,nstrhalf
        do istr_src = 1,istr
          do ist = 1,nst
            do ist_src = 1,nst
              ! Either, the istr_src index is lower than istr,
              ! or they are equal and then, the ist_src is lower than ist.
              if (istr_src .eq. istr .and. ist .le. ist_src) cycle
              ! The one with low index left recieves the contribution of
              ! the one with low index right.
              if ((ist .eq. ist_u) .eqv. (ist_src .eq. ist_u)) then
                drv%drint_bdrf(ist_src,ist,istr_src,istr,:,:,:) = drv%drint_bdrf(ist_src,ist,istr_src,istr,:,:,:) + drv%drint_bdrf(ist,ist_src,istr,istr_src,:,:,:)
              else
                drv%drint_bdrf(ist_src,ist,istr_src,istr,:,:,:) = drv%drint_bdrf(ist_src,ist,istr_src,istr,:,:,:) - drv%drint_bdrf(ist,ist_src,istr,istr_src,:,:,:)
              endif
              drv%drint_bdrf(ist,ist_src,istr,istr_src,:,:,:) = 0.0
            enddo
          enddo
        enddo
      enddo
    endif

    ! Convert derivatives with respect to tau and ssa to derivatives with
    ! respect to taus and taua
    if (deriv_base) then
      do ideriv = 1,nlay_deriv_base
        drv%drint_taus(ideriv,1:nst,1:ngeo) = drint_ssa(ideriv,:,:) * self%pix%chain_taus_ssa(set%ilay_deriv_base(ideriv)) + drint_tau(ideriv,:,:)
        drv%drint_taua(ideriv,1:nst,1:ngeo) = drint_ssa(ideriv,:,:) * self%pix%chain_taua_ssa(set%ilay_deriv_base(ideriv)) + drint_tau(ideriv,:,:)
      enddo
    endif

    if (execution(3)) call fields_dynamic_close(.true.,self%fld,stat)

  end subroutine lintran_execute_with_derivatives ! }}}

  !> Cleans up Lintran.
  !!
  !! All arrays that are allocated in lintran_init are cleaned up.
  subroutine lintran_close(self,stat) ! {{{

    use lintran_types_module, only: lintran_class
    use lintran_grid_module, only: grid_close
    use lintran_geometry_module, only: geometry_close
    use lintran_pixel_module, only: pixel_close
    use lintran_tlc_module, only: tlc_close
    use lintran_fields_module, only: fields_close

    implicit none

    ! Input and output.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance to be deconstructed.
    integer, intent(out) :: stat !< Error code.

    ! Initialize error flag.
    stat = 0

    call fields_close(self%flag_derivatives,self%fld,stat)
    call tlc_close(self%nst,self%flag_derivatives,self%tlc,stat)
    call pixel_close(self%flag_derivatives,self%pix,stat)
    call geometry_close(self%nst,self%geo,stat)
    call grid_close(self%nst,self%flag_derivatives,self%grd,stat)

  end subroutine lintran_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(progression,flag_derivatives,self,stat) ! {{{

    use lintran_types_module, only: lintran_class
    use lintran_grid_module, only: grid_close
    use lintran_geometry_module, only: geometry_close
    use lintran_pixel_module, only: pixel_close
    use lintran_tlc_module, only: tlc_close
    use lintran_fields_module, only: fields_close

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was being initialized for derivatives.
    type(lintran_class), intent(inout) :: self !< Main Lintran instance to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    if (progression .ge. 5) call fields_close(flag_derivatives,self%fld,stat)
    if (progression .ge. 4) call tlc_close(self%nst,flag_derivatives,self%tlc,stat)
    if (progression .ge. 3) call pixel_close(flag_derivatives,self%pix,stat)
    if (progression .ge. 2) call geometry_close(self%nst,self%geo,stat)
    if (progression .ge. 1) call grid_close(self%nst,flag_derivatives,self%grd,stat)

  end subroutine fail ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail_nakajima_multiscattering_without_derivatives(progression,thermal_emission,atm,set,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_atmosphere, lintran_settings

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: thermal_emission !< Flag for having thermal emission in the input atmosphere.
    type(lintran_atmosphere), intent(inout) :: atm !< Temporary atmosphere to be cleaned up.
    type(lintran_settings), intent(inout) :: set !< Temporary settings to be cleaned up.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (progression .ge. 13) deallocate(set%execute_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(atm%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (thermal_emission) then
      if (progression .ge. 11) deallocate(atm%planck_curve,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 10) deallocate(atm%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 9) deallocate(atm%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(atm%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(atm%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(atm%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(atm%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(atm%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(atm%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(atm%taua,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(atm%taus,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail_nakajima_multiscattering_without_derivatives ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail_nakajima_multiscattering_with_derivatives(progression,thermal_emission,deriv_base,deriv_ph,artificial_drint_phase_ssg,atm,set,drv_new,drv,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: thermal_emission !< Flag for having thermal emission in the input atmosphere.
    logical, intent(in) :: deriv_base !< Flag for having any derivatives with respect to taus and taua.
    logical, intent(in) :: deriv_ph !< Flag for having any derivatives with respect to phase issues.
    logical, intent(in) :: artificial_drint_phase_ssg !< Flag for that a field is allocated artificially.
    type(lintran_atmosphere), intent(inout) :: atm !< Temporary atmosphere to be cleaned up.
    type(lintran_settings), intent(inout) :: set !< Temporary settings to be cleaned up.
    type(lintran_derivatives), intent(inout) :: drv_new !< Temporary derivatives to be cleaned up.
    type(lintran_derivatives), intent(inout) :: drv !< User derivatives, where possibly a temporary additional field was created.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (artificial_drint_phase_ssg) then
      if (progression .ge. 24) deallocate(drv%drint_phase_ssg,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 23) deallocate(drv_new%drint_emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 22) deallocate(drv_new%drint_emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 21) deallocate(drv_new%drint_bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (deriv_ph) then
      if (progression .ge. 20) deallocate(drv_new%drint_phase_ssg,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (deriv_base) then
      if (thermal_emission) then
        if (progression .ge. 19) deallocate(drv_new%drint_planck_curve,stat=ierr)
        if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      endif
      if (progression .ge. 18) deallocate(drv_new%drint_taua,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 17) deallocate(drv_new%drint_taus,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (deriv_ph) then
      if (progression .ge. 16) deallocate(set%ilay_deriv_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (deriv_base) then
      if (progression .ge. 15) deallocate(set%ilay_deriv_base,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 14) deallocate(set%differentiate_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 13) deallocate(set%execute_stokes,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(atm%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (thermal_emission) then
      if (progression .ge. 11) deallocate(atm%planck_curve,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 10) deallocate(atm%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 9) deallocate(atm%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(atm%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(atm%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(atm%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(atm%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(atm%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(atm%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(atm%taua,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(atm%taus,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail_nakajima_multiscattering_with_derivatives ! }}}

end module lintran_module
