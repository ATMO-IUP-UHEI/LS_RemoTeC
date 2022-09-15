!> \file lintran_pixel_module.f90
!! Module file for the optical properties of the atmosphere.

!> Module for the optical properties of the atmosphere.
!!
!! This modules handles everything that does depend on the optical properties of
!! the atmosphere. Those are different for each wavelength, so where the grid
!! and geometry can be recycled for different wavelengths, everything that is
!! handled in this module does change. Unlike lintran_tlc_module, the results
!! of this module are still physical terms, such as optical depth. They are not
!! yet converted into radiative transfer.
module lintran_pixel_module

  implicit none

contains

  !> Constructs space for the atmosphere.
  !!
  !! The atmosphere contains the optical properties in a pre-processed,
  !! but they are still optical properties. Also, some chain rule factors
  !! are constructed. Pixel here, stands for one calculation, or one
  !! wavelength pixel.
  subroutine pixel_init(nst,nscat,nstrhalf,nleg,nlay,ngeo,flag_derivatives,pix,stat) ! {{{

    use lintran_constants_module, only: errorflag_allocation
    use lintran_types_module, only: lintran_pixel_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nscat !< Number of independent matrix elements in phase matrix.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, intent(in) :: flag_derivatives !< Flag for possibility of calculating derivatives.
    type(lintran_pixel_class), intent(out), target :: pix !< Internal pixel structure.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr
    integer :: progression

    progression = 0

    ! Optical properties of the atmosphere, possibly pre-processed.
    allocate(pix%tau(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%ssa(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%coefs(nscat,0:nleg,nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%phase_ssg(nscat,nlay,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%bdrf_ssg(nst,nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%bdrf_0(nst,nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%bdrf_v(nst,nst,nstrhalf,0:nleg,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%bdrf(nst,nst,nstrhalf,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%emi(nst,nstrhalf,0:nleg),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%emi_ssg(nst,ngeo),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    allocate(pix%planck_curve(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}

    ! Information for either-or-not taking high Legendre numbers into account.
    allocate(pix%nleg_lay(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}
    ! Two logicals for restricting surface reflection and emission to Fourier number 0 are scalar.

    ! Layer-splitting.
    allocate(pix%nsplit(nlay),stat=ierr) ! {{{
    if (ierr .ne. 0) then
      stat = or(stat,errorflag_allocation)
      call fail(progression,flag_derivatives,pix,stat)
      return
    endif
    progression = progression + 1 ! }}}

    if (flag_derivatives) then

      ! Chain rules for delta-M.
      allocate(pix%chain_ssa_ssa(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_ssa_tau(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_tau_tau(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_ph_ph(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      ! Differentiating with respect to the forward peak.
      allocate(pix%chain_f_ssa(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_f_tau(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_f_ph(nscat,0:nleg,nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_f_ph_elem(nscat,nlay,ngeo),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}

      ! Set pointers that will always remain here. These pointers
      ! are only there to call the same field with a different name.
      pix%chain_zipscat_c => pix%tau ! From scattering zip-integral to C-parameter.
      pix%chain_zipscat_tau => pix%ssa ! From scattering zip-integral to tau.

      ! Chain rules for converting tau/ssa coordintates to taus/taua.
      allocate(pix%chain_taus_ssa(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}
      allocate(pix%chain_taua_ssa(nlay),stat=ierr) ! {{{
      if (ierr .ne. 0) then
        stat = or(stat,errorflag_allocation)
        call fail(progression,flag_derivatives,pix,stat)
        return
      endif
      progression = progression + 1 ! }}}

    else

      progression = progression + 10 ! Fake allocations.

    endif

  end subroutine pixel_init ! }}}

  !> Constructs the atmosphere and arrays for differentiation chain rules.
  !!
  !! The atmosphere is in principle equal to input, but unused additional
  !! Legendre polynomials or such are cut off. Moreover, with &delta;-M, the
  !! optical properties are adapted so that the forward peak is cut off the
  !! phase function. And the optical depths are converted from absorption and
  !! scattering optical depths to extinction optical depth and single-scattering
  !! albedo.
  subroutine pixel_set(nst,nscat,nstrhalf,nleg,nlay,ngeo,execution,flag_derivatives,atm,set,pix) ! {{{

    use lintran_constants_module, only: lowlimit_tau, iscat_a1, iscat_b1, iscat_a2, iscat_a3, iscat_b2, iscat_a4, pol_nst, fullpol_nst
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_pixel_class

    implicit none

    ! Input and output.
    integer, intent(in) :: nst !< Number of Stokes parameters.
    integer, intent(in) :: nscat !< Number of independent matrix elements in phase matrix.
    integer, intent(in) :: nstrhalf !< Half the number of streams.
    integer, intent(in) :: nleg !< Highest Legendre number.
    integer, intent(in) :: nlay !< Number of atmospheric layers.
    integer, intent(in) :: ngeo !< Number of viewing geometries.
    logical, dimension(3), intent(in) :: execution !< Flags for executing single, double and multi-scattering, in that order.
    logical, intent(in) :: flag_derivatives !< Flag for calculating derivatives.
    type(lintran_atmosphere), intent(in) :: atm !< Atmosphere provided by the user.
    type(lintran_settings), intent(in) :: set !< Calculation settings provided by the user.
    type(lintran_pixel_class), intent(inout) :: pix !< Internal pixel structure.

    ! Forward peak (for delta-M).
    real, dimension(nlay) :: forward_peak

    ! Optical depth and single-scattering albedo before delta-M.
    real, dimension(nlay) :: tau_predm
    real, dimension(nlay) :: ssa_predm

    ! Iterators.
    integer :: ileg ! Over Legendre coefficients.
    integer :: ilay ! Over atmospheric layers.
    integer :: igeo ! Over viewing geometries.

    ! Apply delta-M correction. See the following papers:
    ! The Delta-M Method: Rapid Yet Accurate Radiative Flux
    ! Calculations to Strongly Asymmetric Phase Functions
    ! W. J. Wiscombe
    ! J. Atmos. Sci.
    ! For a more basic explanation for a two-stream example, read
    ! The delta-Eddington approximation for radiative flux transfer
    ! J. Joseph, W. J. Wiscombe and J. A. Weiman.
    ! J. Atmos. Sci.

    ! The phase function contains elements from 0 to the number of
    ! streams, actually one too many. The last element is for the delta-m
    ! method, where a forward peak is split off from the phase function.
    ! This means that, instead of setting all higher Legendre elements
    ! to zero, they are all set equal to the last number. As no
    ! information about the high Legendre elements is there, this
    ! approximation must be at least as good as truncation (so zero is
    ! at least equally illogical as the last Legendre coefficient). Here,
    ! the last Legendre coefficient is actually the first Legendre
    ! coefficient that would be truncated without the Delta-M method.

    ! In these comments, where we name Legendre coefficients, we mean their
    ! low definition, (with division by 2l+1 or without multiplication by 2l+1).
    ! Though, in the input in the atmosphere structure, the Legendre coefficients
    ! are defined in their high forms.

    ! With the delta-m method, the phase function is split into a
    ! forward delta-peak and an m-Legendre-coefficient phase function.
    ! P = f*delta + (1-f) * M.
    ! in such a way that Legendre elements 0 to m-1 are preserved, and even
    ! Legendre element m is preserved for the upper left corner I to I.
    ! Then, Legendre coefficients m to infinity for all diagonal matrix elements
    ! are set the same, thus equal to element m for I to I.
    ! Then f is equal to Legendre coefficient m in her low definition.
    ! Then, phase function m has coefficients only from 0 to m-1 and
    ! are defined as
    ! X* = (X-f) / (1-f), where X is the original Legendre coefficient in
    ! the low form (and the result X* is also in her low form).
    ! In the high form, we will muliply left and right sides by 2l+1.
    ! X*h = (Xh-(2l+1)f) / (1-f).
    ! The forward-peak of f can be implemented by removing a fraction
    ! f from the scattering cross section.
    ! For off-diagonal terms, the delta contribution is zero, because the
    ! delta is the identity matrix. Those elements will get the transformation
    ! X*h = Xh / (1-f).
    ! The same will happen for the single-scattering phase function.

    ! Calculate f by turning Legendre coefficient m (for I to I) to the low form.
    ! Note that nleg is m-1, so 2*m+1 is 2*nleg+3.

    ! Translate taus and taua to tau and ssa before delta-M.
    tau_predm = atm%taus + atm%taua
    ! Limits to protect from singularies, especially for pseudo-spherical consitions and
    ! for derivatives. Effectively, a layer with too low tau will have an added contribution
    ! of absorption up to lowlimit_tau. Note that atm%taua is not used after this line.
    where (tau_predm .lt. lowlimit_tau)
      tau_predm = lowlimit_tau
    endwhere
    ssa_predm = atm%taus / tau_predm

    ! Save number of Legendre coefficients per layer for usage after delta-M application.
    ! When bothering about delta-M, the original atm%nleg_lay should be used.
    if (execution(2) .or. execution(3)) pix%nleg_lay = min(atm%nleg_lay,nleg) ! Limits to maximum possible Legendre coefficient.

    if (set%deltam) then

      ! Calculate forward peak. It is the low definition of Legendre coefficient nleg+1
      ! Using the original nleg_lay, because we need a Legendre coefficient that will no
      ! longer exist after application of delta-M.
      where (atm%nleg_lay .gt. nleg)
        forward_peak = atm%coefs(iscat_a1,nleg+1,:) / float(2*nleg+3)
      elsewhere
        forward_peak = 0.0
      endwhere

      ! Calculate new phase function by using the formula for the high
      ! forms: X*h = (Xh-(2l+1)f) / (1-f) for the diagonal terms and
      ! X*h = Xh / (1-f) for the non-diagonal terms.
      if (execution(2) .or. execution(3)) then
        do ileg = 0,nleg
          where (pix%nleg_lay .ge. ileg)
            pix%coefs(iscat_a1,ileg,:) = (atm%coefs(iscat_a1,ileg,:) - (float(2*ileg+1)*forward_peak)) / (1.0 - forward_peak)
          endwhere
          if (pol_nst(nst)) then ! Involving Q and U.
            ! Same where structure.
            where (pix%nleg_lay .ge. ileg)
              pix%coefs(iscat_b1,ileg,:) = atm%coefs(iscat_b1,ileg,:) / (1.0 - forward_peak) ! Off-diagonal.
              pix%coefs(iscat_a2,ileg,:) = (atm%coefs(iscat_a2,ileg,:) - (float(2*ileg+1)*forward_peak)) / (1.0 - forward_peak)
              pix%coefs(iscat_a3,ileg,:) = (atm%coefs(iscat_a3,ileg,:) - (float(2*ileg+1)*forward_peak)) / (1.0 - forward_peak)
            endwhere
          endif
          if (fullpol_nst(nst)) then ! Involving V.
            ! Same where structure.
            where (pix%nleg_lay .ge. ileg)
              pix%coefs(iscat_b2,ileg,:) = atm%coefs(iscat_b2,ileg,:) / (1.0 - forward_peak) ! Off-diagonal.
              pix%coefs(iscat_a4,ileg,:) = (atm%coefs(iscat_a4,ileg,:) - (float(2*ileg+1)*forward_peak)) / (1.0 - forward_peak)
            endwhere
          endif
        enddo
      endif

      ! For single-scattering geometry, it is not that important to peel off
      ! the forward peak, because it will not be used. The dividsion by
      ! 1-f is important, because the angle-resolved phase function works
      ! with the same (reduced) scatteriung coefficients. All independent matrix
      ! elements are handled the same.
      if (execution(1)) then
        do igeo = 1,ngeo
          do ilay = 1,nlay
            pix%phase_ssg(:,ilay,igeo) = atm%phase_ssg(1:nscat,ilay,igeo) / (1.0-forward_peak(ilay))
          enddo
        enddo
      endif

      ! Problems may occur at very low solar and viewing angle and very low
      ! aximuthal difference. Then, single-scattering geometry is actually
      ! quite a forward scattering event.

      ! Now use f to adapt the cross sections and the single-scattering albedo.
      pix%tau = tau_predm - forward_peak*atm%taus ! atm%taus = tau_predm * ssa_predm.
      pix%ssa = (1.0-forward_peak)*ssa_predm / (1.0 - forward_peak*ssa_predm)

      ! Fill chain rules. Note that we are using the tau and ssa before
      ! delta-M. Otherwise, it will become a draconian expression.
      if (flag_derivatives) then

        pix%chain_ssa_ssa = (1.0-forward_peak) / (1.0-forward_peak*ssa_predm)**2.0
        pix%chain_ssa_tau = -forward_peak*tau_predm
        pix%chain_tau_tau = 1.0 - forward_peak*ssa_predm

        ! Chain rules for the phase function.
        pix%chain_ph_ph = 1.0 / (1.0-forward_peak)

        ! Some irrelevant numbers can be multiplied with 0.0, so they are indeed irrelevant.
        ! But uninitialized numbers are dangerous, even if they will be multiplied with zero.
        ! Therefore, this zero-initialization is needed, even outside the if-execution clause.
        pix%chain_f_ph = 0.0
        pix%chain_f_ph_elem = 0.0

        if (execution(2) .or. execution(3)) then
          do ileg = 0,nleg
            where (pix%nleg_lay .ge. ileg)
              pix%chain_f_ph(iscat_a1,ileg,:) = (atm%coefs(iscat_a1,ileg,:) - float(2*ileg+1)) / (1.0-forward_peak)**2.0 / float(2*nleg+3)
            endwhere
            if (pol_nst(nst)) then ! Involving Q and U.
              ! Same where structure.
              where (pix%nleg_lay .ge. ileg)
                pix%chain_f_ph(iscat_b1,ileg,:) = atm%coefs(iscat_b1,ileg,:) / (1.0-forward_peak)**2.0 / float(2*nleg+3) ! Off-diagonal.
                pix%chain_f_ph(iscat_a2,ileg,:) = (atm%coefs(iscat_a2,ileg,:) - float(2*ileg+1)) / (1.0-forward_peak)**2.0 / float(2*nleg+3)
                pix%chain_f_ph(iscat_a3,ileg,:) = (atm%coefs(iscat_a3,ileg,:) - float(2*ileg+1)) / (1.0-forward_peak)**2.0 / float(2*nleg+3)
              endwhere
            endif
            if (fullpol_nst(nst)) then ! Involving V.
              ! Same where structure.
              where (pix%nleg_lay .ge. ileg)
                pix%chain_f_ph(iscat_b2,ileg,:) = atm%coefs(iscat_b2,ileg,:) / (1.0-forward_peak)**2.0 / float(2*nleg+3) ! Off-diagonal.
                pix%chain_f_ph(iscat_a4,ileg,:) = (atm%coefs(iscat_a4,ileg,:) - float(2*ileg+1)) / (1.0-forward_peak)**2.0 / float(2*nleg+3)
              endwhere
            endif
          enddo
        endif
        if (execution(1)) then
          do igeo = 1,ngeo
            do ilay = 1,nlay
              pix%chain_f_ph_elem(:,ilay,igeo) = atm%phase_ssg(1:nscat,ilay,igeo) / (1.0-forward_peak(ilay))**2.0 / float(2*nleg+3)
            enddo
          enddo
        endif

        ! Fill chain rules for differentiating with respect to the forward peak
        ! Note the division by 2*nleg+3 for the conversion between low and high
        ! definition of Legendre coefficient nleg+1.
        pix%chain_f_ssa = (ssa_predm*(ssa_predm-1.0)) / (1.0 - forward_peak*ssa_predm)**2.0 / float(2*nleg+3)
        pix%chain_f_tau = -atm%taus / float(2*nleg+3) ! atm%taus = tau_predm * ssa_predm.

      endif

    else

      ! Ignore the a possible phase coefficient nleg+1.
      if (execution(2) .or. execution(3)) then
        do ilay = 1,nlay
          pix%coefs(:,0:pix%nleg_lay(ilay),ilay) = atm%coefs(1:nscat,0:pix%nleg_lay(ilay),ilay)
        enddo
      endif
      ! Leave tau, ssa and angle-resolved phase function the same.
      pix%tau = tau_predm
      pix%ssa = ssa_predm
      if (execution(1)) pix%phase_ssg = atm%phase_ssg(1:nscat,:,1:ngeo)

    endif

    ! Conversion coordinates tau/ssa to taus/taua, always needed for derivatives
    ! independent of delta-M. These chain rules are applied when delta-M is fully
    ! closed.
    if (flag_derivatives) then
      ! Due to the limit of tau approaching zero, we will not use atm%taua.
      pix%chain_taus_ssa = (tau_predm - atm%taus) / tau_predm**2.0
      pix%chain_taua_ssa = -atm%taus / tau_predm**2.0
    endif

    ! First copy flags for restriction to Fourier number 0.
    pix%bdrf_only_0 = atm%bdrf_only_0
    pix%emi_only_0 = atm%emi_only_0

    ! Copy information that needs no pre-processing.
    pix%sun = atm%sun
    if (execution(3)) then
      if (pix%bdrf_only_0) then
        pix%bdrf(:,:,:,:,0) = atm%bdrf(1:nst,1:nst,:,:,0)
        ! Transpose the matrix, so that source is before destination.
        pix%bdrf(:,:,:,:,0) = reshape(pix%bdrf(:,:,:,:,0) , (/nst,nst,nstrhalf,nstrhalf/), order = (/2,1,4,3/))
      else
        pix%bdrf = atm%bdrf(1:nst,1:nst,:,:,0:nleg)
        ! Transpose the matrix, so that source is before destination.
        pix%bdrf = reshape(pix%bdrf , (/nst,nst,nstrhalf,nstrhalf,nleg+1/), order = (/2,1,4,3,5/))
      endif
    endif
    if (execution(2) .or. execution(3)) then
      if (pix%bdrf_only_0) then
        pix%bdrf_0(:,:,:,0) = atm%bdrf_0(1:nst,1:nst,:,0)
        pix%bdrf_v(:,:,:,0,:) = atm%bdrf_v(1:nst,1:nst,:,0,1:ngeo)
        ! Transpose the matrices, so that source is before destination.
        pix%bdrf_0(:,:,:,0) = reshape(pix%bdrf_0(:,:,:,0) , (/nst,nst,nstrhalf/), order = (/2,1,3/))
        pix%bdrf_v(:,:,:,0,:) = reshape(pix%bdrf_v(:,:,:,0,:) , (/nst,nst,nstrhalf,ngeo/), order = (/2,1,3,4/))
      else
        pix%bdrf_0 = atm%bdrf_0(1:nst,1:nst,:,0:nleg)
        pix%bdrf_v = atm%bdrf_v(1:nst,1:nst,:,0:nleg,1:ngeo)
        ! Transpose the matrices, so that source is before destination.
        pix%bdrf_0 = reshape(pix%bdrf_0 , (/nst,nst,nstrhalf,nleg+1/), order = (/2,1,3,4/))
        pix%bdrf_v = reshape(pix%bdrf_v , (/nst,nst,nstrhalf,nleg+1,ngeo/), order = (/2,1,3,4,5/))
      endif

      if (pix%emi_only_0) then
        pix%emi(:,:,0) = atm%emi(1:nst,:,0)
      else
        pix%emi = atm%emi(1:nst,:,0:nleg)
      endif
    endif
    if (execution(1)) then
      pix%bdrf_ssg = atm%bdrf_ssg(1:nst,1:nst,1:ngeo)
      ! Transpose the matrix, so that source is before destination.
      pix%bdrf_ssg = reshape(pix%bdrf_ssg , (/nst,nst,ngeo/), order = (/2,1,3/))
      pix%emi_ssg = atm%emi_ssg(1:nst,1:ngeo)
    endif
    ! Planck is needed for any execution: 1, 2 or 3. But a flag is needed to turn it on.
    pix%thermal_emission = atm%thermal_emission
    if (pix%thermal_emission) pix%planck_curve = atm%planck_curve

    ! Calculate layer-splitting.
    if (execution(3)) then
      pix%nsplit = max(1,ceiling(pix%tau*pix%ssa/set%taus_split),ceiling(pix%tau*(1.0-pix%ssa)/set%taua_split)) ! This is after delta-M correction.
      ! The atmosphere below a thick cloud does not matter. The layers below such a
      ! thick cloud will not be split and therfore will not cause much additional
      ! computational burden.
      do ilay = nlay-1,1,-1
        if (sum(pix%tau(1:ilay)) .lt. set%tautot_max) then
          ! Atmosphere down to layer ilay+1 must be split, from ilay+2, it need not be split.
          pix%nsplit(ilay+2:nlay) = 1 ! This will do nothing in most cases, when ilay+2 is larger than nlay.
          exit ! Breaks the loop, because the border of tautot_max has been found.
        endif
      enddo
    else
      ! To avoid problems with, for example, stack allocations, we initialize the number of
      ! split layers to what it actually is when there is no split environment.
      pix%nsplit = 1
    endif

  end subroutine pixel_set ! }}}

  !> Cleans up the pixel.
  !!
  !! All arrays that are allocated in pixel_init are cleaned up.
  subroutine pixel_close(flag_derivatives,pix,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_pixel_class

    implicit none

    ! Input and output.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_pixel_class), intent(inout) :: pix !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    if (flag_derivatives) then
      deallocate(pix%chain_taua_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_taus_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_f_ph_elem,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_f_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_f_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_f_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_ph_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_tau_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_ssa_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      deallocate(pix%chain_ssa_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    deallocate(pix%nsplit,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%planck_curve,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    deallocate(pix%tau,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine pixel_close ! }}}

  !> Cleans up garbage that has been created when an error occurs during a routine.
  !!
  !! This is meant to clean up the rubbish that is created before the error occured.
  !! The idea is that the allocation status is the same as before calling the entire
  !! routine where the error occurred.
  subroutine fail(progression,flag_derivatives,pix,stat) ! {{{

    use lintran_constants_module, only: errorflag_deallocation
    use lintran_types_module, only: lintran_pixel_class

    implicit none

    ! Input and output.
    integer, intent(in) :: progression !< Indication of how much garbage has been created.
    logical, intent(in) :: flag_derivatives !< Flag tor whether Lintran was initialized for derivatives.
    type(lintran_pixel_class), intent(inout) :: pix !< Internal pixel structure to be deconstructed.
    integer, intent(inout) :: stat !< Error code.

    ! Error handling.
    integer :: ierr

    ierr = 0 ! Because deallocations happen in if-clauses.

    if (flag_derivatives) then
      if (progression .ge. 23) deallocate(pix%chain_taua_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 22) deallocate(pix%chain_taus_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 21) deallocate(pix%chain_f_ph_elem,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 20) deallocate(pix%chain_f_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 19) deallocate(pix%chain_f_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 18) deallocate(pix%chain_f_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 17) deallocate(pix%chain_ph_ph,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 16) deallocate(pix%chain_tau_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 15) deallocate(pix%chain_ssa_tau,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
      if (progression .ge. 14) deallocate(pix%chain_ssa_ssa,stat=ierr)
      if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    endif
    if (progression .ge. 13) deallocate(pix%nsplit,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 12) deallocate(pix%nleg_lay,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 11) deallocate(pix%planck_curve,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 10) deallocate(pix%emi_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 9) deallocate(pix%emi,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 8) deallocate(pix%bdrf,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 7) deallocate(pix%bdrf_v,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 6) deallocate(pix%bdrf_0,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 5) deallocate(pix%bdrf_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 4) deallocate(pix%phase_ssg,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 3) deallocate(pix%coefs,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 2) deallocate(pix%ssa,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)
    if (progression .ge. 1) deallocate(pix%tau,stat=ierr)
    if (ierr .ne. 0) stat = or(stat,errorflag_deallocation)

  end subroutine fail ! }}}

end module lintran_pixel_module
