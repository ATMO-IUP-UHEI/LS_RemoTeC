!> \file lintran_types_module.f90
!! Module file for all derived types.

!> Module for all derived types used in Lintran.
!!
!! The derived types include the interface and internal types. Internal
!! types have a suffix `_class'. Only the main internal class `lintran_class'
!! must be defined by the user, but the user need not do anything with it except
!! for passing it to all routines. The other non-internal types form the interface
!! and should be provided by the user, or it is output and it is passed back to
!! the user.
module lintran_types_module

  implicit none

  !> Structure describing the optical properties of an atmosphere.
  !!
  !! These contents should be provided by the user. It contains the optical properties
  !! of the atmosphere at one wavelength.
  !! <p>The optical properties are given as optical depths and phase coefficients,
  !! so not as wavelength, refractive indices and particle sizes. That means that
  !! the conversion from physical particles to optical properties (e.g. a Mie module)
  !! is not included in Lintran.
  !! <p>The atmosphere is defined from top to bottom. That means that the first layer
  !! defined is the highest in the atmosphere and the last layer ends at the surface.
  !! The surface itself is included in the model atmosphere, with the Bi-directional
  !! reflection function, and fluorescent emissivity.
  type lintran_atmosphere ! {{{

    real :: sun !< Solar irradiance. Make sure to use the same units as the surface emissivity and the planck curve.
    real, dimension(:), allocatable :: taus !< Scattering optical depth. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: taua !< Absorption optical depth. Dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: coefs !< Phase matrix Legendre coefficients in high definition, which means that Legendre coefficient 1 is in the range \f$ \{-3\dots 3\} \f$, not \f$ \{-1\dots 1\} \f$. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Legendre number (starting at zero). Third dimension: Atmospheric layer. Elements involving non-used Stokes parameters need not be included. At least the required Legendre coefficients must be given. That is \f$ \{0\dots N_s-1\} \f$. Here, \f$ N_s \f$ is the number of streams for which Lintran is initialized. It is possible to provide more Legendre coefficients. They will be ignored. This field need not be allocated for single-scattering.
    real, dimension(:,:,:), allocatable :: phase_ssg !< Scattering matrix at single-scattering geometry. Can be acquired with Legendre coefficients, according to de Haan et al.(1987), equations 68--73. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Atmospheric layer. Third dimension: Viewing geometry. Elements involving non-used Stokes parameters need not be included. This field need to be allocated for multi-scattering.
    real, dimension(:,:,:,:,:), allocatable :: bdrf !< Bi-directional reflection function for internal streams. First dimension: Destination Stokes parameter. Second dimension: Source Stokes parameter. Third dimension: Destination stream. Fourth dimension: Source stream. Fifth dimension: Fourier number (starting at zero). The streams are only defined in one hemisphere, so there are \f$ \frac{N_s}{2} \f$ streams, where \f$ N_s \f$ is the number of streams with which Lintran is initialized. The stream direction cosines can be found in &lt;lintran_instance&gt;%grd%mu when Lintran is initialized. At least Fourier numbers \f$ \{0\dots N_s-1 \} \f$ must be defined (\f$ N_s \f$ is the number of streams with which Lintran is initialized), or when bdrf_only_0 is set to true, only Fourier number 0 must be defined although the Fourier dimension should exist. Defining more Fourier numbers is possible, but they will be ignored. Note that the BDRF must have the correct symmetry relationship. When flipping source and destination streams and Stokes parameters, it should be the same, except when exactly one Stokes parameter is U (3), a minus sign is applied. An example: For Lambertian reflection, all elements for both Stokes parameter I (1) and Fourier index 0 are equal to the albedo, and all other elements are zero. This field need not be allocated for single-scattering.
    real, dimension(:,:,:,:), allocatable :: bdrf_0 !< Bi-directional reflection function with the solar angle. First dimension: Destination Stokes parameter. Second dimension: Source Stokes parameter. Third dimension: Destination stream. Fourth dimension: Fourier number. All properties for the internal BDRF also applies for this BDRF, but there is one internal stream dimension less. Therefore, there is also no longer a symmetry relationship. This field need not be allocated for single-scattering.
    real, dimension(:,:,:,:,:), allocatable :: bdrf_v !< Bi-directional reflection function with the viewing angles. First dimension: Destination Stokes parameter. Second dimension: Source Stokes parameter. Third dimension: Source stream. Fourth dimension: Fourier number. Fifth dimension: Viewing geometry. All properties for the internal BDRF also applies for this BDRF, but there is one internal stream dimension less. Therefore, there is also no longer a symmetry relationship. However, this field has an addition dimension for the viewing geometries. This field need not be allocated for single-scattering.
    real, dimension(:,:,:), allocatable :: bdrf_ssg !< Bi-directional reflection function at single-scattering geometry. This is the BDRF for reflection from the solar direction to the viewing directions, comparable with phase_ssg. First dimension: Destination Stokes parameter. Second dimension: Source Stokes parameter. Third dimension: Viewing geometry. This field need not be allocated for multi-scattering.
    real, dimension(:,:,:), allocatable :: emi !< Surface emissivity. First dimension: Stokes parameter of emitted light. Second dimension: Internal stream. Third dimension: Fourier number. For the Fourier dimension, the same rules apply as for the BDRF. This field need not be allocated for single-scattering.
    real, dimension(:,:), allocatable :: emi_ssg !< Surface emissivity in viewing geometries. First dimension: Stokes parameter of emitted light. Second dimension: Viewing geometry. This field need not be allocated for multi-scattering.
    real, dimension(:), allocatable :: planck_curve !< Position on the Planck curve at the simulated wavelength and the temperature of the atmospheric layer. Make sure to use the same units as the solar irradiance and the surface emissivity. Dimension: Atmospheric layer. This field need not be allocated if thermal emission is turned off.
    integer, dimension(:), allocatable :: nleg_lay !< Highest nonzero Legendre coefficient per layer. Just to save CPU time.
    logical :: bdrf_only_0 !< Flag for that the surface BDRF only has nonzero values for Fourier number 0. Just to save CPU time.
    logical :: emi_only_0 !< Flag for that the fluorescent emissivity only has nonzero values for Fourier number 0. Just to save CPU time.
    logical :: thermal_emission !< Flag for including thermal emission. When turned off, you save CPU time and need not allocate the Planck curve.

  end type lintran_atmosphere ! }}}

  !> Structure for calculation settings.
  !!
  !! This structure contains switches. Most of them trade between speed and precision.
  !! It also contains a selection of which Stokes parameters to calculate and for derivatives,
  !! it contains selection which derivatives should be calculated.
  !! <p>Pseudo-spherical geometry and the selection of single and/or multi
  !! scattering are not included in this structure, but are given as loose parameters when
  !! needed.
  type lintran_settings ! {{{

    ! This structure contains all the switches and settings for running Lintran
    ! It does not include the choice between single-scattering and multi-scattering,
    ! since that is done in two interfaces.
    logical, dimension(:), allocatable :: execute_stokes !< Flags for Stokes parameters to be calculate, for both normal run and derivatives. The order is I, Q, U, V. Non-implemented elements need not be included.
    logical, dimension(:), allocatable :: differentiate_stokes !< Flags for Stokes parameters for which derivatives should be calculated. The order is I, Q, U, V. Non-implemented elements need not be included.
    real :: taus_split !< Maximum scattering optical depth in one sublayer. Layers with higher scattering optical depth are automatically split into multiple sufficiently thin sublayers. Layers are split after eventual application of &delta;-M. Therefore, &delta;-M may reduce the number of sublayers.
    real :: taua_split !< Maximum absorption optical depth in one sublayer. Layers with higher absorption optical depth are automatically split into multiple sufficiently thin sublayers.
    real :: tautot_max !< Maximum total optical depth after which layer-splitting will be disabled. Layers below an optically thick cloud are considered less important, so layer-splitting can be turned off for those layers, increasing the speed. Like for taus_split, this limit is applied after eventual &delta;-M transformation.
    real :: fourier_tolerance !< Ignore higher Fourier numbers when the relative change of the result changes less than this number.
    logical :: split_double !< Switch for solving double-scattering analytically.
    integer :: interpolation !< Interpolation method for intensity inside layer. 0: Averaging. 1: Linear.
    logical :: deltam !< Apply &delta;-M to cut off the forward peak from the phase function.
    integer :: solver !< Choice of which matrix solver to use. 1: Gauss-Seidel. 2: LU-decomposition.
    integer :: gs_maxiter !< Maximum number of iterations for Gauss-Seidel solver. Not needed for LU-decomposition.
    real :: gs_tolerance !< Convergence tolerance for Gauss-Seidel. That is maximum relative change of the downward intensities at the bottom of the atmosphere and the upward intensities at the top of the atmosphere. Not needed for LU-decomposition.
    integer, dimension(:), allocatable :: ilay_deriv_base !< Layer indices for which the derivatives with respect to taus and taua are calculated. Note that layer 1 is at the top of the atmosphere. Need not be allocated when not taking derivatives. Also, leave unallocated if you do not want to differentiate with respect to taus and taua in any layer. Trailing zeros will be ignored, but do not put zeros before nonzero numbers.
    integer, dimension(:), allocatable :: ilay_deriv_ph !< Layer indices for which the derivatives with respect to the phase function is calculated. Note that layer 1 is at the top of the atmosphere. Need not be allocated when not taking derivatives. Also, leave unallocated if you do not want to differentiate with respect to phase functions in any layer. Trailing zeros will be ignored, but do not put zeros before nonzero numbers.
    integer :: nleg_deriv !< Highest Legendre number for which derivatives with respect to phase coefficients are calculated.

  end type lintran_settings ! }}}

  !> Structure for the derivatives.
  !!
  !! The dimensions of all fields are the same as in lintran_atmosphere,
  !! except that the dimensions for atmospheric layers are reduced to the layers
  !! that are selected for calculating the derivatives, and the Legendre coefficicents
  !! in drint_coefs are defined at 0:nleg_deriv, see lintran_settings. Also, each property
  !! gets an additional dimension over viewing Stokes parameters. Furthermore, all fields
  !! that did not yet have a dimension over viewing geometries get one. The viewing geometry
  !! will be the last dimension and the viewing Stokes parameters will be the dimension just
  !! before the one over viewing geometries.
  !! <p>If no layers are selected for differentiation with respect to certain properties
  !! (optical depths or phase functions), then those derivatives also need not be
  !! allocated. Also irrelevant fields, such as the single-scattering phase function when
  !! only calculating multi-scattering need not be allocated (like in lintran_atmosphere).
  !! Fields that are allocated, but did not need to be allocated, will be set to zero.
  type lintran_derivatives ! {{{

    real, dimension(:,:,:), allocatable :: drint_taus !< Derivatives with respect to scattering optical depth. First dimension: Differentiated layer. Second dimension: Viewing Stokes parameter. Third dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: drint_taua !< Derivatives with respect to absorption optical depth. First dimension: Differentiated layer. Second dimension: Viewing Stokes parameter. Third dimension: Viewing geometry.
    real, dimension(:,:,:,:,:), allocatable :: drint_coefs !< Derivatives with respect to Legendre coefficients of the phase matrix. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Legendre coefficient. Third dimension: Differentiated layer. Fourth dimension: Viewing Stokes parameter. Fifth dimension: Viewing geometry.
    real, dimension(:,:,:,:), allocatable :: drint_phase_ssg !< Derivatives with respect to the scattering matrix at single-scattering geometries. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Atmospheric layer. Third dimension: Viewing Stokes parameter. Fourth dimension: Viewing geometry.
    real, dimension(:,:,:,:,:,:,:), allocatable :: drint_bdrf !< Derivatives with respect to the bi-directional reflection function for internal streams. First dimension: Destination Stokes parameter of perturbed BDRF element. Second dimension: Source Stokes parameter of perturbed BDRF element. Third dimension: Destination stream of perturbed BDRF element. Fourth dimension: Source stream of perturbed BDRF element. Fifth dimension: Fourier number of perturbed BDRF element. Sixth dimension: Viewing Stokes parameter. Seventh dimension: Viewing geometry.
    real, dimension(:,:,:,:,:,:), allocatable :: drint_bdrf_0 !< Derivatives with respect to the bi-directional reflection function with the solar angle. First dimension: Destination Stokes parameter of perturbed BDRF element. Second dimension: Source Stokes parameter of perturbed BDRF element. Third dimension: Destination stream of perturbed BDRF element. Fourth dimension: Fourier number of perturbed BDRF element. Fifth dimension: Viewing Stokes parameter. Sixth dimension: Viewing geometry.
    real, dimension(:,:,:,:,:,:), allocatable :: drint_bdrf_v !< Derivatives with respect to the bi-directional reflection function with the viewing angle. First dimension: Destination Stokes parameter of perturbed BDRF element. Second dimension: Source Stokes parameter of perturbed BDRF element. Third dimension: Source stream of perturbed BDRF element. Fourth dimension: Fourier number of perturbed BDRF element. Fifth dimension: Viewing Stokes parameter. Sixth dimension: Viewing geometry.
    real, dimension(:,:,:,:), allocatable :: drint_bdrf_ssg !< Derivative with respect to the bi-directional reflection function at single-scattering geometry. First dimension: Destination Stokes parameter of perturbed BDRF element. Second dimension: Source Stokes parameter of perturbed BDRF element. Third dimension: Viewing Stokes parameter. Fourth dimension: Viewing geometry.
    real, dimension(:,:,:,:,:), allocatable :: drint_emi !< Derivative with respect to surface emissivity. First dimension: Stokes parameter of perturbed emissivity. Second dimension: Stream of perturbed emissivity. Third dimension: Fourier number of perturbed emissivity. Fourth dimension: Viewing Stokes parameter. Fifth dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: drint_emi_ssg !< Derivative with respect to surface emissivity towards the instrument. First dimension: Stokes parameter of perturbed emissivity. Second dimension: Viewing Stokes parameter. Third dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: drint_planck_curve !< Derivative with respect to the normalized position on the Planck curve. First dimension: Atmospheric layer. Second dimension: Viewing Stokes parameter. Third dimension: Viewing geometry.

  end type lintran_derivatives ! }}}

  !> Internal stream grid and related parameters.
  !!
  !! This structure contains all information that only depends on the
  !! internal stream grid. That internal stream grid is determined only by
  !! the number of streams. All contents of this structure do not depend
  !! on the atmosphere and not even on the solar and viewing geometries.
  type lintran_grid_class ! {{{

    ! Internal angle grid.
    real, dimension(:), allocatable :: mu !< Direction cosines of internal streams.
    real, dimension(:), allocatable :: wt !< Weights for internal streams when integrating over the zenith angle.

    ! Prefactors.
    real :: prefactor_dir !< Factor to be applied when considering scattering directly from the sun to the instrument, equal to \f$ \frac{1}{4\pi} \f$.
    real, dimension(:), allocatable :: prefactor_src !< Factor to be applied when considering a source function, equal to \f$ \frac{1}{4\pi\mu} \f$. Dimension: Internal stream of the source.
    real, dimension(:), allocatable :: prefactor_resp !< Factor to be applied for response functions, equal to \f$ \frac{a}{2} \f$. Dimension: Internal stream of the response.
    real, dimension(:,:), allocatable :: prefactor_scat !< Factor to be applied when light is scattered, equal to \f$ \frac{a_s}{2\mu_d} \f$, where \f$ a_s \f$ denotes the integration weight of the source direction and \f$ \mu_d \f$ denotes the direction cosine of the destination direction. First dimension: Source direction. Second dimension: Destination direction.
    real, dimension(:), allocatable :: prefactor_intermediate !< When multiple scatter events are handled at once, the prefactor of all scatter events together are obtained by multiplying the prefactor_scat from the source to the final destination by the prefactor_intermediate of all intermediate streams, equal to \f$ \frac{a}{2\mu} \f$.
    real, dimension(:), allocatable :: prefactor_zip !< Factor to be applied when performing an innerproduct with a normal field on one side and an adjoint field on the other side, equal to \f$ 2\pi a\mu \f$. First dimension: Internal stream. Second dimension: Viewing geometry.

    ! Effective direction cosines constructed from just internal stream angles.
    real, dimension(:,:), allocatable :: effmu_pp !< Effective direction cosine obtained by combining two angles: \f$ \frac{\mu_1\mu_2}{\mu_1+\mu_2} \f$.
    real, dimension(:,:), allocatable :: effmu_mp !< Effective direction cosine obtained by combining two angles of which the first negative: \f$ \frac{\mu_1\mu_2}{\mu_1-\mu_2} \f$.

    ! Relative derivatives of T with respect to tau.
    real, dimension(:), allocatable :: reldiff_t_tau !< Relative derivative of transmission term with respect to the optical depth.

    ! Generalized shperical functions in internal directions. These are used for creating the
    ! phase matrix itself.
    real, dimension(:,:,:), allocatable :: gsf_0 !< Generalized spherical harmonic of type 0, which are also known as associated Legendre polynomials, appearing on the diagonal spots of I and V. These are also used for scalar radiative transport. First dimension: stream. Second dimension: Legendre coefficient. Third dimension: Fourier number.
    real, dimension(:,:,:), allocatable :: gsf_p !< Generalized spherical harmonic of plus type appearing in the matrix on the diagonal spots of Q and U. First dimension: stream. Second dimension: Legendre coefficient. Third dimension: Fourier number.
    real, dimension(:,:,:), allocatable :: gsf_m !< Generalized spherical harmonic of minus type appearing in the matrix on the off-diagonal spots between Q and U. First dimension: stream. Second dimension: Legendre coefficient. Third dimension: Fourier number.

    ! Derivatives with respect to Bi-directional reflection function.
    real, dimension(:,:), allocatable :: diff_c_bck_srf_bdrf !< A factor that has to be applied to the BDRF to get a parameter that is comparable to \f$ \omega Z \f$, equal to \f$ 4\mu_1\mu_2 \f$. Stokes parameters are not involved here.

    ! Derivatives with respect to fluorescent emission for all Stokes parameters.
    ! The reason to include the Stokes parameter is that for the adjoint, there is a minus sign for V.
    real, dimension(:,:), allocatable :: diff_c_emi_nrm_emi !< A factor that has to be applied to the emissivity to get a parameter that is comparable to \f$ \omega Z \f$ for normal (non-adjoint) run. This is equal to \f$ 4\pi\mu \f$ and independent of Stokes parameter. First dimension: Stokes parameter. Second dimension: Stream.
    real, dimension(:,:), allocatable :: diff_c_emi_adj_emi_target !< Eventual different factor that has to be applied to the emissivity to get a parameter that is comparable to \f$ \omega Z \f$ for adjoint run. This is equal to \f$ 4\pi\mu \f$, but a minus sign is applied for V-polarization. First dimension: Stokes parameter. Second dimension: Stream.
    real, dimension(:,:), pointer :: diff_c_emi_adj_emi !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.

  end type lintran_grid_class ! }}}

  !> Solar and viewing geometry and related properties.
  !!
  !! This structure contains all information that depends on the solar
  !! or viewing geometry, but not on the atmosphere. The idea is that
  !! the geometry is constant for all wavelengths, while the atmospheric
  !! properties change with wavelength. These geometry settings can be
  !! initialized once for each spatial pixel, for all wavelengths.
  !! <p>Some variables have a suffix `allgeo'. Then, the parameter is
  !! defined for different viewing geometries. In many times, pointers
  !! exist that point to the correct geometry. The viewing geometry
  !! is always the last dimension.
  type lintran_geometry_class ! {{{

    ! All member variables with suffix _allgeo are defined for all
    ! viewing geometries. Those have an extra dimension over
    ! geometries added at the end. There are no pointers to the current
    ! geometry, because those only exist in the code from the TLC-module
    ! onwards. There, the code thinks there is only one geometry, but
    ! the geometry module is `before' the TLC module, so here, the geometries
    ! all co-exist.

    ! Sun and instrument geometry.
    real :: mu_0 !< Cosine of solar zenith angle.
    real, dimension(:), allocatable :: mu_v_allgeo !< Cosines of viewing zenith angles. Dimension: Viewing geometry.
    real, dimension(:), allocatable :: azi_allgeo !< Azimuthal difference, defined as in De Haan et al. (1987). Saved for the Nakajima wrapper. We could derive the azimuth from the degeneracy array, but that is clumsy. Dimension: Viewing geometry.

    ! Degeneracy factors for different Stokes parameters. Those are
    ! the same as for I, except that the cosine becomes a sine for U and V. And
    ! that for U and V, a minus sign is added to undo the sign switch that is
    ! applied for upward U and V radiation.
    real, dimension(:,:,:), allocatable :: degenerate_allgeo !< Factor to be applied for each Fourier number. For Stokes parameters I and Q, equal to 1 for Fourier number 0, and \f$ 2\cos(m\varphi) \f$ for Fourier numbers above 0 (m denotes the Fourier number). For U and V, it is 0 for Fourier number 0 and \f$ 2\sin(m\varphi) \f$ for Fourier numbers above 0. First dimension: Fourier number. Second dimension: Viewing Stokes parameter. Third dimension: Viewing geometry.
    real, dimension(:), allocatable :: degenerate_ssg_allgeo !< Factor to be applied on single-scattering results. There are no Fourier number there, but the division by \f$ \mu_v \f$ needs to applied, because that one is removed from the prefactors. The sign switch for U and V is not needed for single-scattering. Dimension: Viewing geometry.

    ! Administration on using pseudo spherical atmosphere or plane-parallel.
    ! If the latter is the case, operations can be shortened.
    logical :: plane_parallel !< Flag that saves whether plane-parallel (T) geometry or pseudo-spherical (F) geometry is chosen. This is used to speed up the plane-parallel calculations.

    ! Number of fully visible layers for the sun or the instrument.
    integer :: visible_0 !< Lowest layer (highest index) from which the sun is visible, only non-trivial for pseudo-spherical geometry.
    integer, dimension(:), allocatable :: visible_v_allgeo !< Lowest layer (highest index) visible for the instrument, only non-trivial for pseudo-spherical geometry. Dimension: Viewing geometry.
    integer, dimension(:), allocatable :: visible_0v_allgeo !< Lowest layer (highest index) both visible for the sun and the instrument, equal to the minimum of visible_0 and visisble_v_allgeo. Dimension: Viewing geometry.

    ! Layers to which transmission terms are sensitive, equal to index for
    ! plane-parallel geometry, possibly lower in the atmosphere (higher index)
    ! for pseudo-spherical geometry, or 0 for invisible layers.
    integer, dimension(:), allocatable :: ilay_sensitive_0 !< Lowest layer (highest index) that is passed through when light travels from sun to an interface. Only non-trivial for pseudo-spherical geometry. Dimension: Relevant interface.
    integer, dimension(:,:), allocatable :: ilay_sensitive_v_allgeo !< Lowest layer (highest index) that is passed through when light travels from an interface to the instrument. Only non-trivial for pseudo-spherical geometry. First dimension: Relevant interface. Second dimension: Viewing geometry.
    integer, dimension(:,:), allocatable :: ilay_sensitive_0v_allgeo !< Lowest layer (highest index) that is passed through when light travels from the sun via on interface to the instrument. Only non-trivial for pseudo-spherical geometry. Equal to the maximum of ilay_sensitive_0 and ilay_sensitive_v_allgeo. First dimension: Relevant interface. Second dimension: Viewing geometry.

    real, dimension(:), allocatable :: effmu_p0_ppg !< Effective direction cosine by aggregating the solar direction with a positive internal stream for plane-parallel geometry, therefore only dependent on geometry. Dimension: Internal stream.
    real, dimension(:), allocatable :: effmu_m0_ppg !< Effective direction cosine by aggregating the solar direction with a negative internal stream for plane-parallel geometry, therefore only dependent on geometry. Dimension: Internal stream.
    real, dimension(:,:), allocatable :: effmu_pv_ppg_allgeo !< Effective direction cosine by aggregating a viewing direction with a positive internal stream for plane-parallel geometry, therefore only dependent on geometry. First dimension: Internal stream. Second dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: effmu_mv_ppg_allgeo !< Effective direction cosine by aggregating a viewing direction with a negative internal stream for plane-parallel geometry, therefore only dependent on geometry. First dimension: Internal stream. Second dimension: Viewing geometry.
    real, dimension(:), allocatable :: effmu_0v_ppg_allgeo !< Effective direction cosine by aggregating the solar direction with a viewing direction for plane-parallel geometry, therefore only dependent on geometry. Dimension: Viewing geometry.

    ! Relative derivatives of T with respect to tau for solar or
    ! viewing geometries.
    real, dimension(:,:), allocatable :: reldiff_t_tot_0_tau !< Relative derivative of the cumulative transmission term along the solar direction from the top of the atmosphere to another interface, with respect to the optical depth of a layer. First dimension: Layer in which the optical depth is perturbed. Second dimension: Interface to which the transmission is done. For plane-parallel geometry, it contains only 0 and \f$ -\frac{1}{\mu_0} \f$. For pseudo-spherical geometry, it is much more complicated.
    real, dimension(:,:,:), allocatable :: reldiff_t_tot_v_tau_allgeo !< Relative derivative of the cumulative transmission term along the viewing direction from the top of the atmosphere to another interface, with respect to the optical depth of a layer. First dimension: Layer in which the optical depth is perturbed. Second dimension: Interface to which the transmission is done. Third dimension: Viewing geometry. For plane-parallel geometry, it contains only 0 and \f$ -\frac{1}{\mu_v} \f$. For pseudo-spherical geometry, it is much more complicated.
    real, dimension(:,:,:), allocatable :: reldiff_t_tot_0v_tau_allgeo !< Relative derivative of the combined transmission terms, from the sun to the interface and from the interface back to the instrument, equal to reldiff_t_tot_0_tau plus reldiff_t_tot_v_tau_allgeo. First dimension: Layer in which the optical depth is perturbed. Second dimension: Interface to which the transmission is done. Third dimension: Viewing geometry.

    ! Rotation matrix elements for single-scattering.
    ! Also, the sun side depends on the viewing geometry, since the single-scattering
    ! angle is involved.
    real, dimension(:), allocatable :: rotmat_c_0_allgeo !< Cosine term of the single-scattering rotation matrix on the sun side. Dimension: Viewing geometry.
    real, dimension(:), allocatable :: rotmat_c_v_allgeo !< Cosine term of the single-scattering rotation matrix on the instrument side. Dimension: Viewing geometry.
    real, dimension(:), allocatable :: rotmat_s_0_allgeo !< Sine term of the single-scattering rotation matrix on the sun side. Dimension: Viewing geometry.
    real, dimension(:), allocatable :: rotmat_s_v_allgeo !< Sine term of the single-scattering rotation matrix on the instrument side. Dimension: Viewing geometry.

    ! Generalized shperical functions in solar or viewing direction.
    real, dimension(:,:), allocatable :: gsf_0_0 !< Generalized spherical harmonic of type 0 sampled on solar angle. These are also known as associated Legendre polynomials, appearing on the diagonal spots of I and V. These are also used for scalar radiative transport. First dimension: Legendre coefficient. Second dimension: Fourier number.
    real, dimension(:,:,:), allocatable :: gsf_0_v_allgeo !< Generalized spherical harmonic of type 0 sampled on a viewing angle. These are also known as associated Legendre polynomials, appearing on the diagonal spots of I and V. These are also used for scalar radiative transport. First dimension: Legendre coefficient. Second dimension: Fourier number. Third dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: gsf_p_0 !< Generalized spherical harmonic of type plus sampled on the solar angle. These appear in the matrix on the diagonal spots of Q and U. First dimension: Legendre coefficient. Second dimension: Fourier number.
    real, dimension(:,:,:), allocatable :: gsf_p_v_allgeo !< Generalized spherical harmonic of type plus sampled on a viewing angle. These appear in the matrix on the diagonal spots of Q and U. First dimension: Legendre coefficient. Second dimension: Fourier number. Third dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: gsf_m_0 !< Generalized spherical harmonic of type minus sampled on the solar angle. These appear in the matrix on the off-diagonal spots between Q and U. First dimension: Legendre coefficient. Second dimension: Fourier number.
    real, dimension(:,:,:), allocatable :: gsf_m_v_allgeo !< Generalized spherical harmonic of type minus sampled on a viewing angle. These appear in the matrix on the off-diagonal spots between Q and U. First dimension: Legendre coefficient. Second dimension: Fourier number. Third dimension: Viewing geometry.

    ! Derivatives of C with respect to the Bi-directional reflection function (BDRF).
    real, dimension(:), allocatable :: diff_c_bck_srf_0_bdrf !< A factor that has to be applied to the BDRF to get a parameter that is comparable to \f$ \omega Z \f$, equal to \f$ 4\mu_0\mu \f$, because this is for the BDRF in which the solar direction is involved. It does not depend on the Stokes parameter. Dimension: stream.
    real, dimension(:,:), allocatable :: diff_c_bck_srf_v_bdrf_allgeo !< A factor that has to be applied to the BDRF to get a parameter that is comparable to \f$ \omega Z \f$, equal to \f$ 4\mu_v\mu \f$, because this is for the BDRF in which the viewing direction is involved. It does not depend on the Stokes parameter. First dimension: stream. Second dimension: Viewing geometry.

    ! Derivatives of finite-element C-parameters.
    real, dimension(:), allocatable :: diff_c_bck_srf_0v_elem_bdrf_allgeo !< A factor that has to be applied to the BDRF to get a parameter that is comparable to \f$ \omega Z \f$, equal to \f$ 4\mu_v\mu_0 \f$, because this is for the BDRF in single-scattering geometry. The `elem' stands for `finite element' method, which means that Fourier modes are not used. It does not depend on the Stokes parameter. Dimension: Viewing geometry.
    real, dimension(:), allocatable :: diff_c_emi_v_elem_emi_allgeo !< A factor that has to be applied to the emissivity to get a parameter that is comparable to \f$ \omega Z \f$, for emission directed straight to the instrument. Equal to \f$ 4\pi\mu_v \f$. It does not depend on the Stokes parameter. Dimension: Viewing geometry.

  end type lintran_geometry_class ! }}}

  !> Structure for atmospheric information of one wavelength pixel.
  !!
  !! This structure contains the atmospheric parameters that can change
  !! for each wavelength. This structure needs to be prepared for each
  !! wavelength. This contains the atmospheric input from lintran_atmosphere
  !! with eventual pre-processing. Also, parameters that need no pre-processing
  !! are copied to this structure, so that the entire pre-processed atmosphere
  !! is present in this structure. A few parameters are added to this structure.
  !! Those are needed only in the main Lintran module. For the core radiative
  !! transfer calculation, the lintran_tlc_class is created.
  type lintran_pixel_class ! {{{

    ! All input parameters that also may depend on the wavelength.
    real :: sun !< Solar irradiance.
    real, dimension(:), allocatable :: tau !< Extinction optical depth after eventual &delta;-M transformation. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: ssa !< Single-scattering albedo after eventual &delta;-M transformation. Dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: coefs !< Phase Legendre coefficients. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Legendre number. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: phase_ssg !< Phase function at single-scttering geometry. First dimension: Independent matrix element, The order is I&#8594;I, I&#8594;Q=Q&#8594;I, Q&#8594;Q, U&#8594;U, V&#8594;U=-U&#8594;V, V&#8594;V. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    real, dimension(:,:,:,:,:), allocatable :: bdrf !< Bi-directional reflection function for internal streams. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream. Fifth dimension: Fourier number.
    real, dimension(:,:,:,:), allocatable :: bdrf_0 !< Bi-directional reflection function in which the solar direction is involved. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Destination stream. Fourth dimension: Fourier number.
    real, dimension(:,:,:,:,:), allocatable :: bdrf_v !< Bi-directional reflection function in which a viewing direction is involved. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Fourier number. Fifth dimension: Viewing geometry
    real, dimension(:,:,:), allocatable :: bdrf_ssg !< Bi-directional reflection function at single-scattering geometry. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: emi !< Surface emissivity for internal streams, normalized with the solar irradiance. FIrst dimension: Stokes parameter of emitted light. Second dimension: Stream. Thrid dimension: Fourier number.
    real, dimension(:,:), allocatable :: emi_ssg !< Surface emissivity directed towards the instrument, normalized with the solar irradiance. First dimension: Stokes parameter. Second dimension: Viewing geometry.
    real, dimension(:), allocatable :: planck_curve !< Position on the Planck at the simulated wavelength and the temperature of the atmospheric layer. Dimension: Atmospheric layer. This field need not be allocated if thermal emission is turned off.

    ! Information for either-or-not calculating some less convenient stuff.
    integer, dimension(:), allocatable :: nleg_lay !< Highest nonzero Legendre coefficient per layer that is taken into account, so it will never be higher than the number of Legendre coefficients that exist in the calculation.
    logical :: bdrf_only_0 !< Flag for that the surface BDRF only has nonzero values for Fourier number 0.
    logical :: emi_only_0 !< Flag for that the fluorescent emissivity only has nonzero values for Fourier number 0.
    logical :: thermal_emission !< Flag for including thermal emission.

    ! Administrative array.
    integer, dimension(:), allocatable :: nsplit !< Number of sublayers in which a layer is split. Dimension: Atmospheric layer.

    ! Stuff only for derivatives.

    ! Chain rules for derivatives with delta-M.
    real, dimension(:), allocatable :: chain_ssa_ssa !< Chain rule from derivative with respect to preprocessed single-scattering albedo to derivative with respect to raw single-scattering albedo, only non-trivial with &delta;-M. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: chain_ssa_tau !< Chain rule from derivative with respect to preprocessed optical depth to derivative with respect to raw single-scattering albedo, only non-trivial with &delta;-M. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: chain_tau_tau !< Chain rule from derivative with respect to preprocessed optical depth to derivative with respect to raw optical depth, only non-trivial with &delta;-M. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: chain_ph_ph !< Chain rule from derivative with respect to preprocessed phase element to derivative with respect to raw phase element, only non-trivial with &delta;-M. Dimension: Atmospheric layer.

    ! Chain rules for differentiating to the forward peak.
    real, dimension(:), allocatable :: chain_f_ssa !< Chain rule from derivative with respect to preprocessed single-scattering albedo to derivative with respect to the forward-peak fraction, only non-trivial with &delta;-M. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: chain_f_tau !< Chain rule from derivative with respect to preprocessed optical depth to derivative with respect to the forward-peak fraction, only non-trivial with &delta;-M. Dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: chain_f_ph !< Chain rule from derivative with respect to preprocessed phase element to derivative with respect to the forward-peak fraction, only non-trivial with &delta;-M. First dimension: Independent matrix element. Second dimension: Legendre number. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: chain_f_ph_elem !< Chain rule from derivative with respect to preprocessed phase element in single-scattering geometry to derivative with respect to the forward-peak fraction, only non-trivial with &delta;-M. First dimension: Independent matrix element. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.

    ! Chain rule from scattering zip integrals to scattering parameters or tau.
    real, dimension(:), pointer :: chain_zipscat_c !< Chain rule from scattering zip-integral to a derivative with respect to C, a parameter that looks like \f$ \omega Z \f$. Dimension: Atmospheric layer.
    real, dimension(:), pointer :: chain_zipscat_tau !< Chain rule from scattering zip-integral to a derivative with respect to the optical depth. Dimension: Atmospheric layer.

    ! Chain rule for converting tau/ssa to taus/taua
    ! Chain rule from tau to taus or taua is a trivial one and will not be stored.
    real, dimension(:), allocatable :: chain_taus_ssa !< Chain rule for conversion with respect to the single-scattering albedo to a derivative with respect to scattering optical depth. Dimension: Atmospheric layer.
    real, dimension(:), allocatable :: chain_taua_ssa !< Chain rule for conversion with respect to the single-scattering albedo to a derivative with respect to absorption optical depth. Dimension: Atmospheric layer.

  end type lintran_pixel_class ! }}}

  !> Structure for building blocks of the radiative transfer equation.
  !!
  !! Lintran is split is two parts. In the first part, optical properties of the atmosphere
  !! are transformed into properties that describe what the atmospheric layers do to the
  !! radiation. These properties are named TLC-parameters (hence tlc_class). The T is
  !! for transmission. This is a factor that is applied each time a photon is transmitted
  !! from a layer interface to another one, under a specified zenith angle. The general shape
  !! of these parameters is \f$ e^{-\frac{\tau}{\mu}} \f$. The L is a layer-integral of
  !! the optical depth to integrate for scattering in a layer. The general shape of these
  !! parameters is \f$ \mu(1-e^{-\frac{\tau}{\mu}}) \f$. The C is for the part that is
  !! is put outside the integral over \f$ \tau \f$, therefore relatively constant. It has
  !! a general shape of \f$ \omega Z \f$. But other scatter events such as surface reflection
  !! and fluorescent emission (that is a scatter event in Lintran) is transformed such that
  !! there is a variable that takes the role of \f$ \omega Z \f$. Some direction cosines and
  !! integration weights are put in prefactors, counted as TLC-parameters though it is no
  !! T, L or C. These parameters should strictly be included in C, but to prexerve the
  !! desired symmetry relationships in C, these prefectors are set aside.
  !! <p> The second part of Lintran performs the actual radiative transfer equations.
  !! There, the administration of intensities, sources and response function play a
  !! central role. But becuse of the TLC-parameters, the formulas become a lot easier.
  !! The idea of Lintran is that the transfer, matrix and perturbation modules, those that work
  !! with intensity fields, only use the TLC class is used, not the grid, geometry and pixel
  !! classes. Furthermore, all arithmetic inside the transfer, matrix and perturbation modules
  !! is as simple as a few multiplications, additions and subtractions, with at the most a power
  !! for transmission through multiple identical sublayers.
  !! <p>Besides the TLC-parameters themselves, their derivatives are also included, because
  !! they are necessary for calculations with derivatives. When not running with derivatives
  !! many derivatives of TLC-parameters will not be calculated. However, some of them are
  !! also used as intermediate results for the actual TLC-parameters, so they will also be
  !! calculated for a run without derivatives.
  !! <p>Many members are pointers, possibly to grid or gemoetry members. In such a case,
  !! the values need not be calculated for each wavelength.
  type lintran_tlc_class ! {{{

    ! This module contains all building blocks from the radiative
    ! transport equation. This includes transmission terms (T), and
    ! an integral over the optical depth for the scattering term. The
    ! latter is split in an integral over the atmospheric layer (L) and
    ! a tau-independent constant (C).

    ! This is the only data needed by the radiative transfer code. This
    ! module translates the physical properties, such as 'phase' and 'tau',
    ! into what it means for radiative transfer, such as 'transmission'.

    ! Some layer integrals are quite complicated and may consist of a few
    ! easier building blocks. In that way, there is still some arithmetic
    ! in the transfer module. However, the derivation on which construction
    ! of TLC-parameters is taken, is just looked up from a dictrionary. The
    ! comments say: 'TLC-parameters for <event>', followd by the TLC-parameters
    ! as if they are learned by heart. Derivations of the integrals will
    ! take another document full of mathematics.

    ! It is possible to have different viewing geometries. Pointers will just
    ! point to the current viewing geometry. Allocatables will be calculated
    ! for all geometries and a pointer to the current geometry is added as well
    ! The allocatables will have one extra dimension, at the end, for the geometries.

    ! Pointer to number of split layers from pixel module, for administration
    ! of split layers in the transfer module.
    integer, dimension(:), pointer :: nsplit !< Pointer to number of sublayers for external layer. Dimension: Atmospheric layer.

    ! Administration on using pseudo spherical atmosphere or plane-parallel.
    ! If the latter is the case, operations can be shortened.
    logical, pointer :: plane_parallel !< Pointer to flag for plane-parallel geometry, for speed optimization. Comes from the geometry module.

    ! Administration on layers to which the transmission terms are sensitive to tau.
    integer, dimension(:), allocatable :: ideriv_sensitive_0 !< Differentiation index of lowest layer (highest index) that is passed through when light travels from sun to an interface. Only non-trivial for pseudo-spherical geometry. Dimension: Relevant interface.
    integer, dimension(:), pointer :: ideriv_sensitive_v !< Differentiation index of lowest layer (highest index) that is passed through when light travels from an interface to the instrument. Only non-trivial for pseudo-spherical geometry. Dimension: Relevant interface.
    integer, dimension(:,:), allocatable :: ideriv_sensitive_v_allgeo !< Same as ideriv_sensitive_v, but now for all viewing geometries. First dimension: Relevant interface. Second dimension: Viewing geometry.
    integer, dimension(:), pointer :: ideriv_sensitive_0v !< Differentiation index of lowest layer (highest index) that is passed through when light travels from the sun via on interface to the instrument. Only non-trivial for pseudo-spherical geometry. Equal to the maximum of ilay_sensitive_0 and ilay_sensitive_v_allgeo. Dimension: Relevant interface.
    integer, dimension(:,:), allocatable :: ideriv_sensitive_0v_allgeo !< Same as ideriv_sensitive_0v, but now for all viewing geometries. First dimension: Relevant interface. Second dimension: Viewing geometry.

    ! Integration prefactors that should not depend on the atmosphere.
    real, pointer :: prefactor_dir !< Pointer to prefactor for direct scattering, from geometry module.
    real, dimension(:), pointer :: prefactor_src !< Pointer to prefactor for source functions from grid module. Dimension: Stream.
    real, dimension(:), pointer :: prefactor_resp !< Pointer to prefactor for response function, from geometry module. Dimension: Stream.
    real, dimension(:,:), pointer :: prefactor_scat !< Pointer to prefactor for internal scattering, from grid module. First dimension: Source stream. Second dimension: Destination direction.
    real, dimension(:), pointer :: prefactor_intermediate !< Pointer to prefactor for inserting an intermediate direction for multiple scattering at once, from the grid module. Dimension: Stream.
    real, dimension(:), pointer :: prefactor_zip !< Pointer to prefactor for zipping normal and adjoint fields, from the geometry module. Dimension: Stream.

    ! Effective direction cosines, possible additional building blocks of (L).
    ! Those where the sun or the instrument is involved may depend on layer, because
    ! of pseudo-spherical geometry.
    real, dimension(:,:), allocatable :: effmu_p0 !< Effective direction cosine by aggregating the solar direction with a positive internal stream. First dimension: Internal stream. Second dimension: Atmospheric layer. This parameter depends on the atmospheric layer for pseudo-spherical geometry.
    real, dimension(:,:), allocatable :: effmu_m0 !< Effective direction cosine by aggregating the solar direction with a negative internal stream. First dimension: Internal stream. Second dimension: Atmospheric layer. This parameter depends on the atmospheric layer for pseudo-spherical geometry.
    real, dimension(:,:), pointer :: effmu_pv !< Effective direction cosine by aggregating the viewing direction with a positive internal stream. First dimension: Internal stream. Second dimension: Atmospheric layer. This parameter depends on the atmospheric layer for pseudo-spherical geometry.
    real, dimension(:,:,:), allocatable :: effmu_pv_allgeo !< Same as effmu_pv, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: effmu_mv !< Effective direction cosine by aggregating the viewing direction with a negative internal stream. First dimension: Internal stream. Second dimension: Atmospheric layer. This parameter depends on the atmospheric layer for pseudo-spherical geometry.
    real, dimension(:,:,:), allocatable :: effmu_mv_allgeo !< Same as effmu_mv, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: effmu_pp !< Effective direction cosines by aggregativng two positive stream directions. This is a pointer to a member of the grid class. First dimension: One stream. Second dimension: The other stream.
    real, dimension(:,:), pointer :: effmu_mp !< Effective direction cosines by aggregativng one negative and one positive stream direction. This is a pointer to a member of the grid class. First dimension: The negative stream. Second dimension: The positive stream.
    real, dimension(:), pointer :: effmu_p !< Effective direction cosine for just one positive stream. This is equal to the real direction cosine. Dimension: Stream.

    ! Transmission terms (T) for all layer above current layer.
    real, dimension(:), allocatable :: t_tot_0 !< Cumulative transmission from top of the atmosphere to an interface in the direction of the sun. Dimension: interface.
    real, dimension(:), pointer :: t_tot_v !< Cumulative transmission from top of the atmosphere to an interface (or the other way around) in the direction of the instrument. Dimension: interface.
    real, dimension(:,:), allocatable :: t_tot_v_allgeo !< Same as t_tot_v, but now for all viewing geometries. First dimension: Interface, Second dimension: Viewing geometry.
    real, dimension(:), pointer :: t_tot_0v !< Product of t_tot_0 and t_tot_v, for transmission first from the sun to an interface and then back to the instrument. Dimension: Interface.
    real, dimension(:,:), allocatable :: t_tot_0v_allgeo !< Same as t_tot_0v, but now for all viewing geometries: First dimension: Interface. Second dimension: Viewing geometry.

    ! Single-layer transmission term.
    real, dimension(:,:), allocatable :: t !< Transmission through one layer in an internal stream direction. First dimension: Stream. Second dimension: Atmospheric layer.

    ! transmission for split layers.
    real, dimension(:,:), allocatable :: t_sp !< Transmission in an internal stream through one sublayer of the specified external layer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:), allocatable :: t_sp_0 !< Transmission in solar direction through one sublayer of the specified external layer. Dimension: Atmospheric layer.
    real, dimension(:), pointer :: t_sp_v !< Transmission in viewing direction through one sublayer of the specified external layer. Dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: t_sp_v_allgeo !< Same as t_sp_v, but now for all viewing geometries. First dimension: Atmospheric layer. Second dimension: Viewing geometry.

    ! Layer-integrals (L).
    real, dimension(:), pointer :: l_0v !< Layer integral using an aggregate direction of the solar direction and the viewing direction. Dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_0v_allgeo !< Same as l_0v, but now for all viewing geometries. First dimension: Atmospheric layer. Second dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: l_p0 !< Layer integral using an aggregate direction of the solar direction and a positive internal stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_m0 !< Layer integral using an aggregate direction of the solar direction and a negative internal stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: l_pv !< Layer integral using an aggregate direction of the viewing direction and a positive internal stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_pv_allgeo !< Same as l_pv, but now for all viewing geometries. First dimension: Stream. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    real, dimension(:,:), pointer :: l_mv !< Layer integral using an aggregate direction of the viewing direction and a negative internal stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_mv_allgeo !< Same as l_mv, but now for all viewing geometries. First dimension: Stream. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    ! Single-stream L-parameters (for thermal emission).
    real, dimension(:,:), allocatable :: l_p !< Layer integral using just one positive internal stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:), pointer :: l_v !< Layer integral using just one instrument direction. Dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_v_allgeo !< Same as l_v but now for all viewing geometries. First dimension: Atmospheric layer. Second dimension: Viewing geometry.

    ! Layer-integrals for split layers.
    ! Most layer-integrals are needed for split environment. Actually
    ! there are more L-parameters needed for the split environment than
    ! for the normal situaion.
    real, dimension(:,:), allocatable :: l_sp_p0 !< Layer-integral using a positive stream and the solar direction, integrated over only a sublayer of the external layer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_sp_m0 !< Layer-integral using a negative stream and the solar direction, integrated over only a sublayer of the external layer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: l_sp_pv !< Layer-integral using a positive stream and the viewing direction, integrated over only a sublayer of the external layer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_sp_pv_allgeo !< Same as l_sp_pv, but now for all viewing geometries. First dimension: Stream. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    real, dimension(:,:), pointer :: l_sp_mv !< Layer-integral using a negative stream and the viewing direction, integrated over only a sublayer of the external layer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_sp_mv_allgeo !< Same as l_sp_mv, but now for all viewing geometries. First dimension: Stream. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: l_sp_pp !< Layer-integral using two positive streams, integrated over only a sublayer of the external layer. First dimension: One stream. Second dimension: The other stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_sp_mp !< Layer-integral using one negative and one positive stream, integrated over only a sublayer of the external layer. First dimension: The negative stream. Second dimension: The positive stream. Third dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_1_sp_p0 !< Layer integral, using the solar direction and a positive internal stream, but a linear term is added to the formula. Still, the integral is over a sublayer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: l_1_sp_pv !< Layer integral, using the viewing direction and a positive internal stream, but a linear term is added to the formula. Still, the integral is over a sublayer. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: l_1_sp_pv_allgeo !< Same as l_1_sp_pv, but now for all viewing geometries. First dimension: Stream. Second dimension: Atmospheric layer. Third dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: l_sp_p !< Layer-integral using just one positive stream direction in split enviroment. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_1_sp_p !< Layer integral, using just a positive internal stream, but a linear term is added to the formula. Still, the integral is over a sublayer. First dimension: Stream. Second dimension. Atmospheric layer.

    ! L-parameters for interpolated intensity fields.
    real, dimension(:,:), allocatable :: l_int_sm_p !< Same-layer contribution for split layer integral with an interpolated field, using a positive stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_int_sm_m !< Same-layer contribution for split layer integral with an interpolated field, using a negative stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:), allocatable :: l_int_sm_0 !< Same-layer contribution for split layer integral with an interpolated field, using the solar direction. Dimension: Atmospheric layer.
    real, dimension(:), pointer :: l_int_sm_v !< Same-layer contribution for split layer integral with an interpolated field, using the viewing direction. Dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_int_sm_v_allgeo !< Same as l_int_sm_v, but now for all viewing geometries. First dimension: Atmospheric layer. Second dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: l_int_ac_p !< Across-layer contribution for split layer integral with an interpolated field, using a positive stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_int_ac_m !< Across-layer contribution for split layer integral with an interpolated field, using a negative stream. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:), allocatable :: l_int_ac_0 !< Across-layer contribution for split layer integral with an interpolated field, using the solar direction. Dimension: Atmospheric layer.
    real, dimension(:), pointer :: l_int_ac_v !< Across-layer contribution for split layer integral with an interpolated field, using the viewing direction. Dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: l_int_ac_v_allgeo !< Across as l_int_ac_v, but now for all viewing geometries. First dimension: Atmospheric layer. Second dimension: Viewing geometry.
    ! L-terms that only involve interpolation. For interpolation. But this is always one half.
    real, dimension(:), allocatable :: l_int !< Integral of the any-side weight over the layer. For any healthy interpolation function, this is half the optical depth. Dimension: Atmospheric layer.

    ! Zipping L-parameters, needed for the zip term of the derivatives.
    real, dimension(:), allocatable :: l_zip_sm !< Same-layer contribution of the layer-integral performed when taking the inner product of a normal and adjoint field. Dimension: Stream.
    real, dimension(:), allocatable :: l_zip_ac !< Across-layer contribution of the layer-integral performed when taking the inner product of a normal and adjoint field. Dimension: Stream.

    ! Tau-constant scatter terms (C), including those only needed for the adjoint
    ! Fourier-dependent parameters are not defined for all geometries at once,
    ! because they will be overwritten for a new geometry. After that geometry
    ! they are no longer needed, because in the next Fourier number, they will
    ! be different anyway. The sun and the instrument has her own Stokes parameter,
    ! the sun has one fixed and the instrument has a loop. That Stokes parameter
    ! is not in the dimensions.
    real, dimension(:,:,:,:,:), allocatable :: c_fwd !< Scattering parameter for forward scattering from an internal stream to another one. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream. Fifth dimension: Atmospheric layer.
    real, dimension(:,:,:,:,:), allocatable :: c_bck !< Scattering parameter for backward scattering from an internal stream to another one. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream. Fifth dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_fwd_0_nrm !< Scattering parameter for forward scattering from the solar direction to an internal stream. The relevant solar Stokes parameter is already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_bck_0_nrm !< Scattering parameter for backward scattering from the solar direction to an internal stream. The relevant solar Stokes parameter is already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_fwd_v_nrm !< Scattering parameter for forward scattering from an internal stream to the instrument. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_bck_v_nrm !< Scattering parameter for backward scattering from an internal stream to the instrument. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_fwd_0_adj_target !< Eventual different scattering parameter for forward scattering from an internal stream to the solar direction, to be used for the adjoint. The relevant solar Stokes parameter is already applied. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_bck_0_adj_target !< Eventual different scattering parameter for backward scattering from an internal stream to the solar direction, to be used for the adjoint. The relevant solar Stokes parameter is already applied. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_fwd_v_adj_target !< Eventual different scattering parameter for forward scattering from the instrument to an internal stream, to be used for the adjoint. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: c_bck_v_adj_target !< Eventual different scattering parameter for backward scattering from the instrument to an internal stream, to be used for the adjoint. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), pointer :: c_fwd_0_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: c_bck_0_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: c_fwd_v_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: c_bck_v_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.

    ! Surface reflection.
    real, dimension(:,:,:,:), allocatable :: c_bck_srf !< Scattering parameter for surface reflection from an internal stream to another one. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream.
    real, dimension(:,:), allocatable :: c_bck_srf_0_nrm !< Scattering parameter for surface reflection from the solar direction to an internal stream. The relevant solar Stokes parameter is already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream.
    real, dimension(:,:), allocatable :: c_bck_srf_v_nrm !< Scattering parameter for surface reflection from an internal stream to the instrument. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Source Stokes parameter. Second dimension: Source stream.
    real, dimension(:,:), allocatable :: c_bck_srf_0_adj_target !< Eventual different scattering parameter for surface reflection from an internal stream to the solar direction, to be used for the adjoint. The relevant solar Stokes parameter is already applied. First dimension: Source Stokes parameter. Second dimension: Source stream.
    real, dimension(:,:), allocatable :: c_bck_srf_v_adj_target !< Eventual different scattering parameter for surface reflection from the instrument to an internal stream, to be used for the adjoint. The relevant viewing Stokes parameter and viewing geometry are already applied. First dimension: Destination Stokes parameter. Second dimension: Destination stream.
    real, dimension(:,:), pointer :: c_bck_srf_0_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.
    real, dimension(:,:), pointer :: c_bck_srf_v_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.

    ! Fluorescent emission.
    real, dimension(:,:), allocatable :: c_emi_nrm !< Parameter similar to C for fluorescent emission to internal stream to be used for normal (non-adjoint) calculation. First dimension: Stokes parameter. Second dimension: Stream.
    real, dimension(:,:), allocatable :: c_emi_adj_target !< Eventual different parameter similar to C for fluorescent emission to internal stream to be used for adjoint calculation. First dimension: Stokes parameter. Second dimension: Stream.
    real, dimension(:,:), pointer :: c_emi_adj !< Pointer for the adjoint, which can be a different target array or just the same as the normal one, depending on the symmetry.

    ! Thermal emission.
    real, dimension(:), allocatable :: c_therm !< Parameter similar to C for thermal emission. Dimension: Atmospheric layer.

    ! C-parameters in single-scattering geometry.
    real, dimension(:), allocatable :: c_bck_0v_elem !< Scattering parameter for scattering in single-scattering geometry for relevant Stokes parameters. Dimension: Atmospheric layer.
    real :: c_bck_srf_0v_elem !< Scattering parameter for surface reflection in single-scattering geometry for the relevant Stokes parameters. Scalar.
    real :: c_emi_v_elem !< Parameter similar to C for fluorescent emission directed towards the instrument for the relevant Stokes parameter. Scalar.

    ! Derivatives of TLC-parameters.

    ! Derivatives with respect to the optical depth.

    ! Due to the possible pseudo-spherical gemoetry, some TLC-parameters are defined
    ! on one layer and are dependent on the optical depth on another layer, making
    ! two layer-dimensions. In such a case, the perturbed layer is always just one
    ! dimension before the scattering layer.

    ! Derivatives of effective direction cosines.
    real, dimension(:,:,:), allocatable :: diff_effmu_p0_tau !< Derivative of effmu_p0 with respect to optical depth. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies.
    real, dimension(:,:,:), allocatable :: diff_effmu_m0_tau !< Derivative of effmu_m0 with respect to optical depth. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies.
    real, dimension(:,:,:), pointer :: diff_effmu_pv_tau !< Derivative of effmu_pv with respect to optical depth. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies.
    real, dimension(:,:,:,:), allocatable :: diff_effmu_pv_tau_allgeo !< Same as diff_effmu_pv_tau, but now for all viewing geometries. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies. Fourth dimension: Viewing geometry.
    real, dimension(:,:,:), pointer :: diff_effmu_mv_tau !< Derivative of effmu_mv with respect to optical depth. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies.
    real, dimension(:,:,:,:), allocatable :: diff_effmu_mv_tau_allgeo !< Same as diff_effmu_mv_tau, but now for all viewing geometries. First dimension: Stream. Second dimension: Layer where optical depth is perturbed. Third dimension: Layer for which the effective direction cosine applies. Fourth dimension: Viewing geometry.

    ! Relative derivatives of T.
    real, dimension(:), pointer :: reldiff_t_tau !< Relative derivative of t with respect to the optical depth. Dimension: Stream.
    real, dimension(:,:), pointer :: reldiff_t_tot_0_tau !< Relative derivative of t_tot_0 with respect to optical depth. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Only non-trivial for pseudo-spherical geometry.
    real, dimension(:,:), pointer :: reldiff_t_tot_v_tau !< Relative derivative of t_tot_v with respect to optical depth. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Only non-trivial for pseudo-spherical geometry.
    real, dimension(:,:), pointer :: reldiff_t_tot_0v_tau !< Relative derivative of t_tot_0v with respect to optical depth. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Only non-trivial for pseudo-spherical geometry.

    ! Relative derivatives of T in split environment.
    real, dimension(:,:), allocatable :: reldiff_t_sp_tau !< Relative derivative of t_sp with respect to optical depth. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: reldiff_t_sp_0_tau !< Relative derivative of t_sp_0 with respect to optical depth. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Only non-trivial for pseudo-spherical geometry.
    real, dimension(:,:), pointer :: reldiff_t_sp_v_tau !< Relative derivative of t_sp_v with respect to optical depth. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Only non-trivial for pseudo-spherical geometry.
    real, dimension(:,:,:), allocatable :: reldiff_t_sp_v_tau_allgeo !< Same as reldiff_t_sp_v, but now for all viewing geometries. First dimension: Layer where optical depth is perturbed. Second dimension: Layer through which is transmitted. Third dimension: Viewing geometry.

    ! Absolute derivatives of L.
    real, dimension(:,:,:), allocatable :: diff_l_p0_tau !< Derivative of l_p0 with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:), pointer :: diff_l_pv_tau !< Derivative of l_pv with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:,:), allocatable :: diff_l_pv_tau_allgeo !< Same as diff_l_pv_tau, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated. Fourth dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: diff_l_m0_tau !< Derivative of l_m0 with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:), pointer :: diff_l_mv_tau !< Derivative of l_mv with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:,:), allocatable :: diff_l_mv_tau_allgeo !< Same as diff_l_mv_tau, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated. Fourth dimension: Viewing geometry.
    real, dimension(:,:), pointer :: diff_l_0v_tau !< Derivative of l_0v with respect to otical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:,:), allocatable :: diff_l_0v_tau_allgeo !< Same as diff_l_0v_tau, but now for all viewing geometries. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated. Third dimension: Viewing geometry.
    real, dimension(:,:), pointer :: diff_l_p_tau !< Derivative of l_p with respect to the optical depth. This is equal to the transmission term. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), pointer :: diff_l_v_tau !< Derivative of l_v with respect to the optical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:,:), allocatable :: diff_l_v_tau_allgeo !< Same as diff_l_v_tau but now for all viewing geometries. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated. Third dimension: Viewing geometry.

    ! Absolute derivatives of L in split environemnt.
    real, dimension(:,:,:), allocatable :: diff_l_sp_p0_tau !< Derivative of l_sp_p0 with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:), pointer :: diff_l_sp_pv_tau !< Derivative of l_sp_pv with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:,:), allocatable :: diff_l_sp_pv_tau_allgeo !< Same as diff_l_sp_pv_tau, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated. Fourth dimension: Viewing geometry.
    real, dimension(:,:,:), allocatable :: diff_l_sp_m0_tau !< Derivative of l_sp_m0 with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:), pointer :: diff_l_sp_mv_tau !< Derivative of l_sp_mv with respect to otical depth. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated.
    real, dimension(:,:,:,:), allocatable :: diff_l_sp_mv_tau_allgeo !< Same as diff_l_sp_mv_tau, but now for all viewing geometries. First dimension: Internal stream. Second dimension: Layer on which optical depth is perturbed. Third dimension: Layer that is integrated. Fourth dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: diff_l_sp_p_tau !< Derivative of l_sp_p with respect to the optical depth. First dimension: Stream. Second dimension: Atmospheric layer.

    ! Absolute derivatives of interpolated L.
    real, dimension(:,:), allocatable :: diff_l_int_sm_p_tau !< Derivative of l_int_sm_p with respect to otical depth. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: diff_l_int_sm_m_tau !< Derivative of l_int_sm_m with respect to otical depth. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: diff_l_int_ac_p_tau !< Derivative of l_int_ac_p with respect to otical depth. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: diff_l_int_ac_m_tau !< Derivative of l_int_ac_m with respect to otical depth. First dimension: Stream. Second dimension: Atmospheric layer.
    real, dimension(:,:), allocatable :: diff_l_int_sm_0_tau !< Derivative of l_int_sm_0 with respect to otical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:), pointer :: diff_l_int_sm_v_tau !< Derivative of l_int_sm_v with respect to otical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:,:), allocatable :: diff_l_int_sm_v_tau_allgeo !< Same as diff_l_int_sm_v_tau , but now for all viewing geometries. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated. Third dimension: Viewing geometry.
    real, dimension(:,:), allocatable :: diff_l_int_ac_0_tau !< Derivative of l_int_ac_0 with respect to otical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:), pointer :: diff_l_int_ac_v_tau !< Derivative of l_int_ac_v with respect to otical depth. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated.
    real, dimension(:,:,:), allocatable :: diff_l_int_ac_v_tau_allgeo !< Same as diff_l_int_ac_v_tau , but now for all viewing geometries. First dimension: Layer on which optical depth is perturbed. Second dimension: Layer that is integrated. Third dimension: Viewing geometry.
    real, dimension(:), allocatable :: diff_l_int_tau !< Derivative of l_int with respect to the optical depth. Dimension: Atmospheric layer.

    ! Derivatives with respect to the single-scattering albedo.

    ! Derivatives of C.
    real, dimension(:,:,:,:,:), allocatable :: diff_c_fwd_ssa !< Derivative of c_fwd with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream. Fifth dimension: Atmospheric layer.
    real, dimension(:,:,:,:,:), allocatable :: diff_c_bck_ssa !< Derivative of c_bck with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Destination Stokes parameter. Third dimension: Source stream. Fourth dimension: Destination stream. Fifth dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_fwd_0_nrm_ssa !< Derivative of c_fwd_0_nrm with respect to the single-scattering albedo. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_bck_0_nrm_ssa !< Derivative of c_bck_0_nrm with respect to the single-scattering albedo. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_fwd_v_nrm_ssa !< Derivative of c_fwd_v_nrm with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_bck_v_nrm_ssa !< Derivative of c_bck_v_nrm with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_fwd_0_adj_ssa_target !< Eventual different derivative of c_fwd_0_adj with respect to the single-scattering albedo. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_bck_0_adj_ssa_target !< Eventual different derivative of c_bck_0_adj with respect to the single-scattering albedo. First dimension: Destination Stokes parameter. Second dimension: Destination stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_fwd_v_adj_ssa_target !< Eventual different derivative of c_fwd_v_adj with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), allocatable :: diff_c_bck_v_adj_ssa_target !< Eventual different derivative of c_bck_v_adj with respect to the single-scattering albedo. First dimension: Source Stokes parameter. Second dimension: Source stream. Third dimension: Atmospheric layer.
    real, dimension(:,:,:), pointer :: diff_c_fwd_0_adj_ssa !< Choice for the adjoint, just taking over the normal array or point to a new one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: diff_c_bck_0_adj_ssa !< Choice for the adjoint, just taking over the normal array or point to a new one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: diff_c_fwd_v_adj_ssa !< Choice for the adjoint, just taking over the normal array or point to a new one, depending on the symmetry.
    real, dimension(:,:,:), pointer :: diff_c_bck_v_adj_ssa !< Choice for the adjoint, just taking over the normal array or point to a new one, depending on the symmetry.

    ! Other derivatives.

    ! There is only a chain rule, the derivatives of C with respect to the
    ! phase function. That is the single-scattering albedo, so this derivative
    ! will point to the single-scattering albedo.
    real, dimension(:), pointer :: diff_c_ph !< Derivative of a C-parameter with respect to the phase function. Dimension: Atmospheric layer.

    ! Pointers to correct Fourier index.
    real, dimension(:,:), pointer :: diff_c_bck_srf_bdrf !< Derivative of c_bck_srf with respect to the bi-directional reflection function. Does not depend on the Stokes parameters. First dimension: Source stream. Second dimension: Destination stream.
    real, dimension(:), pointer :: diff_c_bck_srf_0_bdrf !< Derivative of c_bck_srf_0 with respect to the bi-directional reflection function. Does not depend on the Stokes parameters. Dimension: Destination stream.
    real, dimension(:), pointer :: diff_c_bck_srf_v_bdrf !< Derivative of c_bck_srf_v with respect to the bi-directional reflection function. Does not depend on the Stokes parameters. Dimension: Source stream.
    real, dimension(:,:), pointer :: diff_c_emi_nrm_emi !< Derivative of c_emi_nrm with respect to the emissivity. First dimension: Stokes parameter. Second dimension: Stream.
    real, dimension(:,:), pointer :: diff_c_emi_adj_emi !< Derivative of c_emi_adj with respect to the emissivity. First dimension: Stokes parameter. Second dimension: Stream.

    ! Derivatives of thermal emissivity.
    real, dimension(:), allocatable :: diff_c_therm_ssa !< Derivative of c_therm with respect to the single-scattering albedo.
    real, dimension(:), allocatable :: diff_c_therm_planck !< Derivative of c_therm with respect to the Planck point.

    ! Derivatives of finite-element C-parameters.
    real, dimension(:), allocatable :: diff_c_bck_0v_elem_ssa !< Derivative of c_bck_0v_elem with respect to single-scattering albedo. Dimension: Atmospheric layer.
    real, pointer :: diff_c_bck_srf_0v_elem_bdrf !< Derivative of c_bck_srf_0v_elem with respect to bi-direction reflection function. Does not depend on the Stokes parameters. Scalar.
    real, pointer :: diff_c_emi_v_elem_emi !< Derivative of c_emi_v_elem with respect to emissivity. Does not depend on the Stokes parameters. Scalar.

  end type lintran_tlc_class ! }}}

  !> Structure for intensity fields.
  !!
  !! This structure contains some intensity fields to work with. These
  !! can be the source, the response function and an intensity field.
  !! Because stuff is never as easy as it looks, it contains several
  !! fields where permanent or temporary stuff can be dumped. This
  !! stuff can vary to such an extent that it is not worthwhile to give
  !! meaningful names to the fields. In lintran_module, where the members
  !! of this structure are explicitly referred to, comments are used to
  !! indicate what is what.
  !! <p>Some arrays are analytic and contain a fixed number of atmospheric
  !! layers. Some are interpolated and benefit from layer splitting and
  !! therefore, the fields can grow larger or smaller. These fields are
  !! allocated and deallocated during the calculation.
  !! <p> All fields are four-dimensional. First dimension: Stokes parameter.
  !! Second dimension: Stream. Third dimension: down or up (0 or 1). Fourth
  !! dimension: interface (0 to number of layers).
  !! However, for the matrix, the Stokes and stream dimensions are combined,
  !! because some compilers do not allow eight dimensions. Therefore, the matrix
  !! shifts are also set up in that way to stay consistent with the matrix.
  type lintran_fields_class ! {{{


    ! Analytic fields.
    real, dimension(:,:,:,:), allocatable :: a1 !< Analytic field 1.
    real, dimension(:,:,:,:), allocatable :: a2 !< Analytic field 2.
    real, dimension(:,:,:,:), allocatable :: a3 !< Analytic field 3. Only used with derivatives.

    ! Interpolated fields.
    real, dimension(:,:,:,:), pointer :: i1 !< Interpolated field 1.
    real, dimension(:,:,:,:), pointer :: i2 !< Interpolated field 2.
    real, dimension(:,:,:,:), pointer :: i3 !< Interpolated field 3. Only used with derivatives.
    real, dimension(:,:,:,:), pointer :: i4 !< Interpolated field 4. Only used with derivatives.

    ! The matrix, in eight-dimensional band structure, where dimensions 1 and 2 and
    ! dimensions 5 and 6 are merged. But the physical meaning is eight dimensions.
    ! First four dimensions: Source, where the fourth dimension is only
    ! 0:1, for the two relevant layers
    ! Last four dimensions: Destination in full (0:nsplit_total) representation.
    real, dimension(:,:,:,:,:,:), pointer :: mat !< The matrix. First three dimensions: The source, where the Stokes and steam dimensions are meged. Last three dimensions: The destination, where also the Stokes and stream dimensions are merged. For the source, the layers are reduced to just 0 and 1.

    ! Matrix shifts (for use of banded matrix).
    integer, dimension(:,:,:), pointer :: shifts !< Because in the matrix, the source layer dimension is reduced to a 0--1 domain, simply reshaping it to two-dimensional does not yet give a two-dimensional matrix. For each destiation index, only a limited number of source indices are defined. The square matrix is applied by writing a number of zeros, then the defined values and continue with zeros again. The number of zeros to write at the beginning is stored in this array. The dimensions stand for the destination index, but it is defined in the original four-dimensional destination (Stokes, stream, down/up, interface). Similarly to the matrix, the Stokes and stream dimensions are combined. This is not necessary, but it is just to stay consistent with the matrix, the closest relative of this shifts array.

  end type lintran_fields_class ! }}}

  !> Main structure of Lintran.
  !!
  !! This structure contains all internal structures of Lintran and all dimension
  !! sizes that should be saved. When Lintran is initialized, or a geometry is provided
  !! to Lintran, the contents of this structure change. This structure also stores
  !! important intermediate results while performing a calculation.
  type lintran_class ! {{{

    ! Dimension sizes that Lintran has to remember from the initialization.
    integer :: nst !< Number of Stokes parameters. This can be 1 (scalar), 3 (polarization without V) or 4 (polarization with V).
    integer :: nscat !< Number of independent matrix element in scattering matrix, directly derived from the number of Stokes parameters.
    integer :: nstrhalf !< Half number of streams (\f$ \frac{N_s}{2} \f$). That is the number of streams for each upward and downward radiation. Lintran is initialized with a number of streams that must be even.
    integer :: nleg !< Highest Legendre number. That is equal to \f$ N_s-1 \f$. For higher Legendre numbers, the phase function cannot remain normalized with the given number of streams.
    integer :: nlay !< Number of atmospheric layers, without layer-splitting.
    integer :: ngeo !< Number of viewing geometries.

    ! Flag for whether or not Lintran has been initialized for calculating derivatives.
    logical :: flag_derivatives !< Flag for whether or not Lintran is initialized for calculating derivatives. If so, it is possible to calculate derivatives, though a calculation without derivatives remains possible. If not, only calculation without derivatives can be performed. Turning off this flag consumes less memory and may be a bit quicker for that reason.

    ! Objects from Lintran.
    type(lintran_grid_class) :: grd !< Internal grid structure. See lintran_grid_class.
    type(lintran_geometry_class) :: geo !< Internal geometry structure. See lintran_geometry_class.
    type(lintran_pixel_class) :: pix !< Internal pixel structure. See lintran_pixel_class.
    type(lintran_tlc_class) :: tlc !< Internal TLC structure. See lintran_tlc_class.
    type(lintran_fields_class) :: fld !< Internal field structure. See lintran_fields_class.

  end type lintran_class ! }}}

end module lintran_types_module
