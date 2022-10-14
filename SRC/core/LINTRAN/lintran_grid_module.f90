!> \file lintran_grid_module.f90
!! Module file for the internal stream grid.

!> Module for the internal stream grid.
!!
!! The internal grid of streams is constructed in this module. Furthermore,
!! everything that only depends on the internal stream grid is calculated as well.
!! The number of streams with which Lintran is initialized fully determines what
!! happens in this module.
module lintran_grid_module

  implicit none

contains

  !> Constructs the stream grid and everything that depends on just that.
  !!
  !! Everything that is calculated in this module depends only on the number
  !! of streams. The solar and viewing geometry and the atmosphere do not play
  !! a role here. This grid can be initialized once at the beginning, and used
  !! for all the calculations, for any geometry and any atmosphere. Only when
  !! changing the number of streams, the initialization must be redone.
  subroutine grid_init(nst,nstrhalf,nleg,flag_derivatives,grd,stat) ! {{{

    use lintran_constants_module, only: minimum_mudiff, pi, flip_v, warningflag_instablegrid, warningflag_vulnerablegrid, errorflag_allocation, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_grid_class), intent(out), target :: grd !< Internal grid structure.
    integer, intent(inout) :: stat !< Error code.

    ! Iterators.
    integer :: m ! Fourier number.
    integer :: istr ! Over streams.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    allocate(grd%mu(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%wt(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! grd%prefactor_dir is scalar.
    allocate(grd%prefactor_src(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%prefactor_resp(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%prefactor_scat(nstrhalf,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%prefactor_intermediate(nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%effmu_pp(nstrhalf,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%effmu_mp(nstrhalf,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%gsf_0(nstrhalf,0:nleg,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (pol_nst(nst)) then ! Only for vector.
      allocate(grd%gsf_p(nstrhalf,0:nleg,0:nleg),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,grd,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(grd%gsf_m(nstrhalf,0:nleg,0:nleg),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,grd,stat)
        return
      endif
      progression = progression + 1 ! }}}
    else
      progression = progression + 2 ! Fake allocations.
    endif
    allocate(grd%diff_c_bck_srf_bdrf(nstrhalf,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(grd%diff_c_emi_nrm_emi(nst,nstrhalf),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(nst,progression,flag_derivatives,grd,stat)
      return
    endif
    progression = progression + 1 ! }}}
    if (fullpol_nst(nst)) then ! With V, forward and adjoint lose some symmetry.
      allocate(grd%diff_c_emi_adj_emi_target(nst,nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,grd,stat)
        return
      endif
      progression = progression + 1 ! }}}
      grd%diff_c_emi_adj_emi => grd%diff_c_emi_adj_emi_target
    else
      progression = progression + 1 ! Fake allocation.
      grd%diff_c_emi_adj_emi => grd%diff_c_emi_nrm_emi ! Adjoint is the same as normal.
    endif

    ! For derivatives.
    if (flag_derivatives) then

      allocate(grd%prefactor_zip(nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,grd,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(grd%reldiff_t_tau(nstrhalf),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(nst,progression,flag_derivatives,grd,stat)
        return
      endif
      progression = progression + 1 ! }}}

    else

      progression = progression + 2 ! Fake allocations.

    endif

    ! Set up double Gaussian quadrature.
    ! the sign of grd%mu is chosen that they are always positive.
    ! The solar and viewing angles are not sampled, because they
    ! are not part of the double-Gaussian grid.
    call double_gauss_grid(nstrhalf,grd,stat)
    ! No need for a progression of the error, because the
    ! double Gauss grid can only raise warnings, no errors.

    ! Check grid for possible instable numbers that occur if different streams differ
    ! only by a small amount. Preferrably, the internal streams are far enough away so that
    ! solar and viewing geometries can fit in between properly.
    do istr = 1,nstrhalf-1
      if ((grd%mu(istr)-grd%mu(istr+1)) .lt. minimum_mudiff) then
        ! An unstable grid overwrites the warning for a vulnerable grid.
        stat = stat - and(stat,warningflag_vulnerablegrid)
        stat = or(stat,warningflag_instablegrid)
      else if ((grd%mu(istr)-grd%mu(istr+1)) .lt. 5.0 * minimum_mudiff) then
        if (and(stat,warningflag_instablegrid) .eq. 0) stat = or(stat,warningflag_vulnerablegrid)
      endif
    enddo
    ! Check if last stream is far enough from zero.
    if (grd%mu(nstrhalf) .lt. minimum_mudiff) then
      stat = stat - and(stat,warningflag_vulnerablegrid)
      stat = or(stat,warningflag_instablegrid)
    else if (grd%mu(nstrhalf) .lt. 5.0 * minimum_mudiff) then
      if (and(stat,warningflag_instablegrid) .eq. 0) stat = or(stat,warningflag_vulnerablegrid)
    endif

    ! Calculate integration prefactors.

    ! The rules for prefactors for a photon that start at the sun and
    ! ends in the instrument, with zero or more directions in between.
    ! - 1/(4*pi) for the source.
    ! - weight/2 for all directions except the sun and the instrument.
    ! - 1/abs(mu) for all direction except the sun (including the instrument).

    ! Implementation.
    ! - Each scatter event takes care of the weight of the source direction
    !   and for dividing by the cosine of the destination direction.
    ! - The source factor 1/(4*pi) is used for the scatter event from the
    !   sun as substitute for the weight/2.
    ! - The response does not take responsibility over the division by the
    !   instrument direction cosine. That is moved to the degeneracy array,
    !   to make all prefactors geometry-independent.

    ! Direct scattering: Scatter from source to instrument.
    grd%prefactor_dir = 0.25 / pi

    ! Source: Scatter from source to stream.
    grd%prefactor_src = 0.25 / (pi*grd%mu)

    ! Response: Scatter from stream to instrument.
    grd%prefactor_resp = 0.5*grd%wt

    ! Scattering: From stream 1 to stream 2.
    ! The source is the left dimension. The destinatino is the right dimension.
    do istr = 1,nstrhalf
      grd%prefactor_scat(:,istr) = 0.5*grd%wt / grd%mu(istr)
    enddo

    ! Adding intermediate direction: Insert stream.
    grd%prefactor_intermediate = 0.5*grd%wt / grd%mu

    ! Stuff for derivatives. A lot of these variables are used as intermediate
    ! results for the forward run as well.

    ! Prefactor for zipping forward and adjoint.
    ! The future of the photon is stored in the adjoint as if it is a source
    ! and some scatters. Therefore, the division by mu was done in both the
    ! last scatter event of the adjoint and in the last scatter event of the
    ! forward, while the 0.5*wt_i has neither been done by the adjoint nor by
    ! the forward. Moreover, both the forward and the adjoint did a source
    ! function, with 1/(4*pi).
    ! Correcting this results in 4*pi * 0.5*wt_i * mu_i = 2*pi*wt_i*mu_i.
    if (flag_derivatives) grd%prefactor_zip = 2.0*pi*grd%wt*grd%mu

    ! Calculate effective direction cosines for two angles.
    ! Compare these to the effective direction cosines in the geometry
    ! module, where the solar and viewing angles are involved. Here,
    ! only grid streams are involved, so these arrays do not depend on
    ! the geometry.
    do istr = 1,nstrhalf
      grd%effmu_pp(:,istr) = (grd%mu*grd%mu(istr)) / (grd%mu + grd%mu(istr))
      ! Dodging division by zero: The line should be as follows:
      ! grd%effmu_mp(:,istr) = (grd%mu*grd%mu(istr)) / (grd%mu - grd%mu(istr)),
      ! but that is a division by zero for diagonal terms.
      ! The code is designed such that the diagonal terms are not touched, so a NaN
      ! would not affect the result. But if the division by zero is caught as an error,
      ! it will not work.
      grd%effmu_mp(1:istr-1,istr) = (grd%mu(1:istr-1)*grd%mu(istr)) / (grd%mu(1:istr-1) - grd%mu(istr))
      grd%effmu_mp(istr+1:nstrhalf,istr) = (grd%mu(istr+1:nstrhalf)*grd%mu(istr)) / (grd%mu(istr+1:nstrhalf) - grd%mu(istr))
      grd%effmu_mp(istr,istr) = 1.0 ! A safe one (not zero) to reduce chances of further.
      ! zero divisions. This value will never be used in the code, so it is possible to
      ! put any value here. Except for zero, because someone is going to divide by this value.
    enddo

    ! Relative derivatives of T with respect to tau, only needed for derivatives.
    if (flag_derivatives) grd%reldiff_t_tau = -1.0 / grd%mu

    ! Generalized spherical functions.
    ! First dimension: Legendre number.
    ! Second dimension: Stream (internal), do not confuse the 0 with the sun. This 0 belongs to the type of the generalized spherical function.
    ! Third dimension: Fourier number.
    do m = 0,nleg
      call calculate_grid_generalized_spherical_functions(nst,nstrhalf,nleg,m,grd)
    enddo

    ! Derivatives with respect to bi-directional reflection function (BDRF)
    ! and fluorescent emission. Will be used for forward run as well.

    ! Surface reflection.
    ! C = BDRF * 4*mu_src*mu_dest
    do istr = 1,nstrhalf
      grd%diff_c_bck_srf_bdrf(:,istr) = 4.0 * grd%mu * grd%mu(istr)
    enddo

    ! Emissivity field.
    ! C = emission * 4*pi*mu.
    do istr = 1,nstrhalf
      grd%diff_c_emi_nrm_emi(:,istr) = 4.0*pi * grd%mu(istr) ! Indifferent of Stokes parameter.
    enddo
    ! With V, we have to copy this array for the adjoint, but flip V.
    ! Wihtout V, we should do nothing with the adjoint, because it is
    ! the same and the array just points to the same data as the normal
    ! array.
    if (fullpol_nst(nst)) then
      do istr = 1,nstrhalf
        grd%diff_c_emi_adj_emi(:,istr) = grd%diff_c_emi_nrm_emi(:,istr) * flip_v
      enddo
    endif

  end subroutine grid_init ! }}}

  !> Constructs the double Gaussian stream grid.
  !!
  !! Such a grid is normalized for Legendre polynomials up to \f$ N_s-1 \f$, which
  !! is the maximum possible. This normalization implies that when integrating
  !! a Legendre polynomial, the result is zero, except for Legendre polynomial
  !! number 0, for which the integral is 2 (when integrating over both downard
  !! and upward streams). Besides the stream direction cosines, the integration
  !! weights are constructed as well. These take the role of \f$ d\mu \f$ when
  !! integrating over \f$ \mu \f$.
  subroutine double_gauss_grid(nstrhalf,grd,stat) ! {{{

    use lintran_constants_module, only: pi, dg_max_iterations, dg_tolerance, warningflag_dgnonconvergence
    use lintran_types_module, only: lintran_grid_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    type(lintran_grid_class), intent(inout) :: grd !< Internal grid structure.
    integer, intent(inout) :: stat !< Error code.

    ! Reference: Davis, P.J. and P. Rabinowitz,
    ! Methods of Numerical Integration,
    ! Academic Press, New York, pp. 87, 1975.

    ! method: Compute the abscissae as roots of the Legendre
    ! polynomial P-sub-M using a cubically convergent
    ! refinement of Newton's method. Compute the
    ! weights from EQ. 2.7.3.8 of Davis/Rabinowitz. Note
    ! that Newton's method can very easily diverge; only a
    ! very good initial guess can guarantee convergence.
    ! The initial guess used here has never led to divergence
    ! even for M up to 1000.

    ! This routine is not protected against numeric instable divisions, but
    ! it is probably impossible to trigger an instability before the double-Gauss
    ! grid fails to converge. Moreover, this routine is called at initialization, so
    ! the crash will at least not happen at an unexpected moment.

    ! Old and new solutions.
    real :: x ! Cosine of the angle.
    real :: xnew ! New cosine of the angle.

    ! Help variables.
    real :: t ! An angle as intermediate solution for the first guess of x.
    real :: inv_s2 ! One divided by sine angle squared.

    ! Iterators.
    integer :: ileg ! Legendre polynomial number.
    integer :: istr ! Over streams.
    integer :: iter ! For iterative convergence.

    ! Dimension sizes.
    integer :: nstrquarter ! Because it is a double Gaussian, it is first
    ! calculated at half nstrhalf (is nstrquarter) streams.

    ! Legendre polynomial value at given angle cosine for recursive
    ! definition.
    real :: leg ! Current.
    real :: leg_m1 ! Previous.
    real :: leg_m2 ! Before previous.

    ! Derivatives of Legendre polynomials.
    real :: dleg ! First derivative.
    real :: d2leg ! Second derivative.

    ! For a two-stream method, the result is easy.
    if (nstrhalf .eq. 1) then
      grd%mu = (/0.5/)
      grd%wt = (/1.0/)
      return
    endif

    ! If nstrhalf is not equal to one, a rather complicated calculation
    ! is executed. First, the first half of the streams are calculated.
    nstrquarter = nstrhalf / 2

    do istr = 1,nstrquarter

      ! Initial guess for istr-th root of Legendre polynomial,
      ! from Davis/Rabinowitz (2.7.3.3a).
      t = (4*istr - 1) * pi / (4*nstrhalf + 2)
      x = cos(t + float(nstrhalf-1) / (8 * nstrhalf**3) / tan(t))

      ! Upward recurrence for Legendre polynomials.
      do iter = 1,dg_max_iterations

        ! Initialize as if ileg=1, Legendre polynomials 0 and 1.
        leg_m1 = 1.0 ! Legendre polynomial 0.
        leg = x ! Legendre polynomial 1.
        do ileg = 2,nstrhalf
          ! Move back previous polynomials.
          leg_m2 = leg_m1
          leg_m1 = leg
          ! Calculate new Legendre polynomial.
          leg = (float(2*ileg-1) * x * leg_m1 - float(ileg-1) * leg_m2) / float(ileg)
        enddo
        ! Now, leg is Legendre polynomial nstrhalf and leg_m1 is Legendre
        ! polynomial nstrhalf-1.

        ! Newton method.
        inv_s2 = 1.0 / (1.0 - x**2.0)

        ! First derivative.
        dleg = nstrhalf * (leg_m1 - x * leg) * inv_s2
        ! Second derivative with the Legendre equation.
        d2leg = (2.0 * x * dleg - (nstrhalf * (nstrhalf+1)) * leg) * inv_s2

        ! Calculate new x by minimizing the value of the nstrhalf'th
        ! Legendre polynomial, using two orders of derivatives.
        xnew = x - (leg / dleg) * (1.0 + (leg / dleg) * d2leg / (2.0 * dleg))

        ! Check for convergence.
        if (abs(xnew-x) .gt. dg_tolerance) then
          x = xnew ! And continue to next iteration.
        else
          exit ! Convergnece: Break the iterative loop.
        endif

      enddo

      ! Check if convergence has been achieved by checking the iteration number.
      if (iter .gt. dg_max_iterations) then
        stat = or(stat,warningflag_dgnonconvergence)
      endif

      ! Iteration finished, calculate weights, abscissae for (-1,1).
      grd%mu(istr) = x
      grd%wt(istr) = 2.0 / (inv_s2 * (nstrhalf * leg_m1)**2.0)

      ! Copy value for negative x. Note that the domain of the grid will later
      ! be transformed so that the ultimate angle cosines (mu) are all positive.
      grd%mu(nstrhalf+1-istr) = -grd%mu(istr)
      grd%wt(nstrhalf+1-istr) = grd%wt(istr)

    enddo

    ! Set middle angle and weight for rules of odd order.
    if (mod(nstrhalf,2) .eq. 1) then
      grd%mu(nstrquarter + 1) = 0.0
      grd%wt(nstrquarter + 1) = 2.0 / product((/(float(istr)/float(istr-1),istr=3,nstrhalf,2)/))**2.0
    endif

    ! Convert the grid from domain [-1..1] to [0..1]. This is where the
    ! Gaussian grid becomes double.
    do istr = 1,nstrhalf
      grd%mu(istr) = 0.5 * grd%mu(istr) + 0.5
      grd%wt(istr) = 0.5 * grd%wt(istr)
    enddo

  end subroutine double_gauss_grid ! }}}

  !> Calculates the generalized spherical harmonics for just the internal streams.
  !!
  !! In this function, generalized spherical harmonics for internal stream angles are calculated.
  !! These are required to calculate any phase matrix element, as well as their derivatives.
  subroutine calculate_grid_generalized_spherical_functions(nst,nstrhalf,nleg,m,grd) ! {{{

    use lintran_constants_module, only: pol_nst
    use lintran_types_module, only: lintran_grid_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: m !< Fourier number.
    type(lintran_grid_class), intent(inout) :: grd !< Internal grid structure.

    integer :: ileg ! Legendre number.
    integer :: ibin ! Inplied do-loop iterator for binomial coefficient.
    real :: prefactor ! Common prefactor for initial Legendre polynomial.

    ! Building blocks of the p and m polynomials, who have their own recursive function.
    real, dimension(nstrhalf) :: temppolyp2
    real, dimension(nstrhalf) :: temppolym2

    ! Recursive definition of the p and m building blocks, save history of the last
    ! two ones. The 'm' after the underscore stands for minus in the sense of 'l'-values,
    ! so that is the previous polynomial. The 'm' (or 'p') before the underscore stands
    ! for the type of function.
    real, dimension(nstrhalf) :: temppolyp2_m1
    real, dimension(nstrhalf) :: temppolyp2_m2
    real, dimension(nstrhalf) :: temppolym2_m1
    real, dimension(nstrhalf) :: temppolym2_m2

    ! Calculate the common prefactor.
    ! The prefactor is 2^(-m) * sqrt(2m over m), where the former is a
    ! power and the latter is the root of a binomial coefficient.
    prefactor = 2.0**(-m)*sqrt(product((/(float(m+ibin)/float(ibin),ibin=1,m)/)))

    ! Set everything at zero for l < m.
    grd%gsf_0(:,0:m-1,m) = 0.0
    if (pol_nst(nst)) then ! Only for vector.
      ! Write a zero for Legendre coefficients 0 and 1, because that is always the
      ! case for gsf_p and gsf_m. In the module, Legendre numbers 0 and 1 are not
      ! calculated, so they are left uninitialized.
      grd%gsf_p(:,0:1,m) = 0.0
      grd%gsf_m(:,0:1,m) = 0.0
      ! If not already set to zero, write a zero for l < m.
      grd%gsf_p(:,2:m-1,m) = 0.0
      grd%gsf_m(:,2:m-1,m) = 0.0
    endif

    ! Set up the polynomials so that recursion can do the rest.

    ! Polynomial 0, set up ileg=m and ileg=m+1.
    grd%gsf_0(:,m,m) = prefactor * sqrt(1.0-grd%mu**2.0)**m
    if (m .lt. nleg) grd%gsf_0(:,m+1,m) = sqrt(float(2*m+1))*grd%mu*grd%gsf_0(:,m,m)

    ! Full recursion for higher Legendre numbers.
    do ileg = m+2,nleg
      grd%gsf_0(:,ileg,m) = (float(2*ileg-1)*grd%mu*grd%gsf_0(:,ileg-1,m) - sqrt(float((ileg-1)**2-m**2))*grd%gsf_0(:,ileg-2,m))/sqrt(float(ileg**2-m**2))
    enddo

    if (pol_nst(nst)) then ! Only for vector.

      ! Polynomial p2 and m2, the lowest possible ileg is 2, or higher if m is higher
      ! than 2.

      ! This lowest ileg will be put in the _m2 polynomial and one higher in _m1.

      select case (m)
        case (0)
          temppolyp2_m2 = -0.25 * sqrt(6.0) * (1.0-grd%mu**2.0)
          temppolym2_m2 = -0.25 * sqrt(6.0) * (1.0-grd%mu**2.0)
        case (1)
          temppolyp2_m2 = 0.5 * sqrt(1.0-grd%mu**2.0) * (1.0+grd%mu)
          temppolym2_m2 = -0.5 * sqrt(1.0-grd%mu**2.0) * (1.0-grd%mu)
        case (2)
          temppolyp2_m2 = -0.25 * (1.0+grd%mu)**2.0
          temppolym2_m2 = -0.25 * (1.0-grd%mu)**2.0
        case default
          temppolyp2_m2 = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-grd%mu**2.0)**(m-2) * (1.0+grd%mu)**2.0
          temppolym2_m2 = -prefactor * sqrt(float(m*(m-1)) / float((m+1)*(m+2))) * sqrt(1.0-grd%mu**2.0)**(m-2) * (1.0-grd%mu)**2.0
      endselect

      ! The ileg is set to two, except for the `case default', where ileg is set to m.
      ileg = max(m,2)

      ! Write results.
      grd%gsf_p(:,ileg,m) = 0.5 * (temppolym2_m2 + temppolyp2_m2)
      grd%gsf_m(:,ileg,m) = 0.5 * (temppolym2_m2 - temppolyp2_m2)

      if (m .ne. nleg) then

        ! Acquire next polynomials, for ileg+1. Polynomial _m2 is now the previous one
        ! and the next-to-previous one does not exist yet.
        temppolyp2_m1 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*grd%mu-2.0*m) * temppolyp2_m2) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2_m1 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*grd%mu+2.0*m) * temppolym2_m2) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))

        ! Write results for the next Legendre number.
        grd%gsf_p(:,ileg+1,m) = 0.5 * (temppolym2_m1 + temppolyp2_m1)
        grd%gsf_m(:,ileg+1,m) = 0.5 * (temppolym2_m1 - temppolyp2_m1)

      endif

      ! And process all ilegs that are higher. Note that the iterator ileg is
      ! one lower than the actual Legendre number, that is how the recursion is
      ! defined and the formulas are complicated.
      do ileg = max(m,2)+1,nleg-1

        temppolyp2 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*grd%mu-2.0*m) * temppolyp2_m1 - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolyp2_m2) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))
        temppolym2 = ((2.0*ileg+1.0)*(ileg*(ileg+1.0)*grd%mu+2.0*m) * temppolym2_m1 - (ileg+1.0)*sqrt((ileg**2.0-4.0)*(ileg**2.0-m**2.0)) * temppolym2_m2) / (ileg*sqrt(((ileg+1.0)**2.0-4.0)*((ileg+1.0)**2.0-m**2.0)))

        ! Write results. The ileg is still one lower than what it should be.
        grd%gsf_p(:,ileg+1,m) = 0.5 * (temppolym2 + temppolyp2)
        grd%gsf_m(:,ileg+1,m) = 0.5 * (temppolym2 - temppolyp2)

        ! Progress previous results for next recursive iteration.
        temppolyp2_m2 = temppolyp2_m1
        temppolym2_m2 = temppolym2_m1
        temppolyp2_m1 = temppolyp2
        temppolym2_m1 = temppolym2

      enddo

    endif

  end subroutine calculate_grid_generalized_spherical_functions ! }}}

  !> Cleans up the grid.
  !!
  !! All arrays that are allocated in grid_init are cleaned up.
  subroutine grid_close(nst,flag_derivatives,grd,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_grid_class), intent(inout) :: grd !< Internal grid structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      deallocate(grd%reldiff_t_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(grd%prefactor_zip,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (fullpol_nst(nst)) then
      deallocate(grd%diff_c_emi_adj_emi_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(grd%diff_c_emi_nrm_emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%diff_c_bck_srf_bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then
      deallocate(grd%gsf_m,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(grd%gsf_p,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(grd%gsf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%effmu_mp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%effmu_pp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%prefactor_intermediate,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%prefactor_scat,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%prefactor_resp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%prefactor_src,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%wt,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(grd%mu,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine grid_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(nst,progression,flag_derivatives,grd,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_grid_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_grid_class), intent(inout) :: grd !< Internal grid structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (flag_derivatives) then
      if (progression .ge. 16) deallocate(grd%reldiff_t_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 15) deallocate(grd%prefactor_zip,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (fullpol_nst(nst)) then
      if (progression .ge. 14) deallocate(grd%diff_c_emi_adj_emi_target,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 13) deallocate(grd%diff_c_emi_nrm_emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(grd%diff_c_bck_srf_bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (pol_nst(nst)) then
      if (progression .ge. 11) deallocate(grd%gsf_m,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 10) deallocate(grd%gsf_p,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 9) deallocate(grd%gsf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(grd%effmu_mp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(grd%effmu_pp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(grd%prefactor_intermediate,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(grd%prefactor_scat,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(grd%prefactor_resp,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(grd%prefactor_src,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(grd%wt,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(grd%mu,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail ! }}}

end module lintran_grid_module
