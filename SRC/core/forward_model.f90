!------------------------------------------------------------------------------
module forward_model_module
   use header_module
   use spectral_response_module, only: instrument_response, response_internal, spectral_response_stored, spectral_response_custom, &
                                       spectral_response_lbl
   use rad_trans_module, only: &
      derivatives, rad_trans_intf, &
      aero, set_altdis, get_aerosol_properties_lognormal, &
      Mie_lut, cirrus_table, read_aerosol_netcdf, read_cirrus_netcdf, &
      window_ini, settings_flags, file_paths, read_settings, read_win_xsdb
   implicit none
   private

   !*** Public types
   public :: derivatives, aero, Mie_lut, cirrus_table, instrument_response, spectral_response_stored, &
             window_ini, window_spectrum, settings_flags, file_paths

   !*** Public procedures
  public :: forward_model_hi, forward_model_lo, read_aerosol_netcdf, read_cirrus_netcdf, set_altdis, read_settings, read_win_xsdb, &
             get_aerosol_properties_lognormal
   !------------------------------------------------------------------------------
   !> @brief Molecular absorbers
   !> @details Target absorbers are retrieved, global absorbers are not
   !------------------------------------------------------------------------------
   type, public :: absorbers
      integer :: ntype_target                               !< Number of target absorbers
      integer :: ntype_global                               !< Number of interfering absorbers
      integer, dimension(:), allocatable :: type_x_target   !< Target absorbers (dim: ntype_target)
      integer, dimension(:), allocatable :: type_x_global   !< Interfering absorbers (dim: ntype_global)
   end type absorbers

   !------------------------------------------------------------------------------
   !> Retrieval window structure, different for each ground pixel
   !------------------------------------------------------------------------------
   type :: window_spectrum
      !> Geometry
      real(double) :: sza, iza, phi, observer_height

      !*** Derivatives
      !> Lo-res reflectance/radiance derivatives wrt absorber partial column (Dim: nwave_lo)
      real(double), dimension(:, :, :), allocatable :: derivatives_lo
      !> Lo-res reflectance/radiance derivatives wrt temperature offset (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: derivT_lo
      !> Lo-res reflectance/radiance derivatives wrt temperature offset (Dim: nwave_lo,nlay)
      real(double), dimension(:, :), allocatable :: derivP_lo
      !> Lo-res reflectanc/radiancee derivatives wrt 0th order albedo (Dim: nwave_lo,albflag)
      real(double), dimension(:, :), allocatable :: deriv_albedo_lo
      !> Lo-res reflectance/radiance derivatives wrt 0th order intensity offset (Dim: nwave_lo, IOffFlag)
      real(double), dimension(:, :), allocatable :: deriv_IOff_lo
      !> Lo-res reflectance/radiance derivatives wrt 0th shift  of spectrum(Dim: nwave_lo)
      real(double), dimension(:), allocatable :: deriv_specshift0_lo
      !> Lo-res reflectanc/radiancee derivatives wrt 1st shift of spectrum (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: deriv_specshift1_lo
      !> Lo-res reflectance/radiance derivatives wrt 2nd shift of spectrum (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: deriv_specshift2_lo
      !> Lo-res reflectanc/radiancee derivatives wrt oth shift of cross sections (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: deriv_sunshift0_lo
      !> Lo-res reflectance/radiance derivatives wrt aerosol parameters (Dim: nwave_lo,naer)
      real(double), dimension(:, :), allocatable :: deriv_aerosol_lo
      !> Lo-res reflectance/radiance derivatives wrt fluorescence (Dim: nwave_lo, 2)
      real(double), dimension(:, :), allocatable :: deriv_Fs_lo

      !*** Absorber properties:
      !> Partial columns of absorbers, initial guess (dim: natm, ntype)
      real(double), dimension(:, :), allocatable ::dv_x
      !> Partial columns of absorbers, iteratively updated by inversion (dim: ntype,natm)
      real(double), dimension(:, :), allocatable ::x_molec

      !> additive constants to the spectrum:, values plus their derivatives (=1.0)
      real(double), dimension(:), allocatable :: Ioff
      !> Additive term to the spectrum representing fluorescence emission
      real(double), dimension(2) :: Fs

      !*** Spectral variables
      !> Low-res wavelength grid (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: wavelength_lo
      !> Measured spectrum (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: spectrum
      !> Measured spectrum variance (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: spectrum_cov
      !> Lo-res measured sun spectrum (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: sun_spectrum_sat_lo
      !> Measured reflectance (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: reflectance_meas
      !> Measured reflectance variance (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: reflectance_meas_cov
      !> Shifted lo-res wavelength grid (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: wavelength_lo_new
      !> Lo-res modelled spectrum (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: reflectance_lo
      !> Lo-res sun reference (Dim: nwave_lo)
      real(double), dimension(:), allocatable :: sun_spectrum_ref_lo
      !> Shifted hi-res wavelength grid
      real(double), dimension(:), allocatable :: wavelength_hi_new
      !> Shifted hi-res solar spectrum
      real(double), dimension(:), allocatable :: sun_spectrum_ref_hi
      !> hi-res solar spectrum from satellite measurement (eg from deconvolution)
      real(double), dimension(:), allocatable :: sun_spectrum_sat_hi

!*** Multiplicative weighting factors for stokes coefficients provided by NASA L1b
      real(double), dimension(:), allocatable :: measurement_stokesc

      !> Albedo, nth order (albflag)
      real(double), dimension(:), allocatable :: albedo

      !*** Instrument line shape:
      !> delta(wavelength) grid of instrument line shape per wavelength pixel (dim: nils)
      real(double), dimension(:), allocatable :: ilswave
      !> instrument line function per wavelength pixel (dim: nils)
      real(double), dimension(:), allocatable :: ilsfunction
      !> upper limit index on hi-res wavelength grid (dim: nwave_lo)
      integer, dimension(:), allocatable :: ie_store
      !> lower limit index on hi-res wavelength grid (dim: nwave_lo)
      integer, dimension(:), allocatable :: is_store
      !> Array for instrument response (dim: nwave_lo, nils)  (only needed for wavelength dependent function)
      real(double), dimension(:, :), allocatable :: resp_store

      !*** Aerosol and cirri:
      !> Total optical depth at window center
      real(double) :: ot
      !> Cirrus optical depth at window center
      real(double) :: cot

      !> Number of ILS points
      integer :: nils

      !> Number of lo-res wavelengths
      integer :: nwave_lo

   end type window_spectrum

contains
   !------------------------------------------------------------------------------
   !> @brief Compute high-resolution model reflectance and derivatives
   !------------------------------------------------------------------------------
   subroutine forward_model_hi(aero_lut, cirrus_lut, &
                               flag, glintflag, &
                               sza, iza, phi, wspeed, &
                               iwin, &
                               nrt, &
                               nlay, &
                               absorb, &
                               atm_xs, &
                               dvair, &
                               dvair_old, &
                               vmr_h2o, &
                               play_old, &
                               win_ini, &
                               win, &
                               wavelength_hi, &
                               reflectance_hi, &
                               aerosol, &
                               minaotflag, &
                               mincotflag, &
                               nder, &
                               deriv_hi, &
                               ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
      !*** Input
      type(Mie_lut), intent(in) :: aero_lut
      type(cirrus_table), intent(in) :: cirrus_lut
      type(settings_flags), intent(in) :: flag
      integer, intent(in) :: glintflag, iwin, nrt, nlay
      type(absorbers), intent(in) :: absorb
      real(double), intent(in) :: sza, iza, phi, wspeed
      type(atmosphere), intent(in) :: atm_xs
      real(double), dimension(:), intent(in) :: dvair
      real(double), dimension(:), intent(in) :: dvair_old
      real(double), dimension(:), intent(in) :: vmr_h2o
      real(double), dimension(:), intent(in) :: play_old
      type(window_ini), intent(in) :: win_ini
      real(double), dimension(:), intent(in) :: wavelength_hi
      integer, intent(in) :: minaotflag, mincotflag
      integer, dimension(:), intent(in) :: nder
      type(aero), dimension(:), intent(in) :: aerosol
      logical, intent(in) :: deriv_flag
      !   real(double), dimension(:), intent(in) :: albedo
      !   real(double), dimension(:,:), intent(in) :: x_molec
      !*** Output
      !   real(double), intent(out) :: ot, cot
      !*** Input/output
      type(window_spectrum), intent(inout) :: win
      type(derivatives), intent(inout) :: deriv_hi
      !*** Output
      integer, intent(out) :: ExitXSFlag, MaxOTFlag, ierr
      real(double), dimension(:, :), intent(out) :: reflectance_hi
      !*** local variables
      integer :: i, j, k, l, m, n, off, l1, l2, n1, n2, n3, n4, natm, ntype_aer, naer
      real(double) :: w
      real(double), dimension(:, :), allocatable :: rint_fine      ! Modelled reflectance (nwave_hi,nstokes)
      type(derivatives) :: deriv_rt
      real(double), dimension(:), allocatable ::taua_tot
      real(double), parameter :: lambda0 = 755d0, A1 = 1.445d0, A2 = 0.868d0, lambda1 = 736.8d0, &
                                 lambda2 = 685.2d0, sigma1 = 21.2d0, sigma2 = 9.55d0, lambda_o2lo = 759.2, lambda_o2up = 770.8
      real(double) :: continuum
      integer :: nsc
      !---------------------------------------------------------------------

      !*** Initialize
      ierr = 0
      natm = atm_xs%n

      !      phi = DABS(iaz-saz) ! 180.D0 - DABS(DBLE(iaz-saz))
      !***Call radiative transfer
      call rad_trans_intf( &
         aero_lut, cirrus_lut, &
         flag, glintflag, iwin, &
         sza, iza, phi, wspeed, &
         aerosol, &
         minaotflag, mincotflag, &
         nrt, nder, atm_xs, dvair, dvair_old, vmr_h2o, play_old, &
         win_ini, &
         win%albedo, win%x_molec, &
         wavelength_hi, &
         rint_fine, deriv_rt, taua_tot, win%ot, win%cot, ExitXSFlag, MaxOTFlag, ierr, deriv_flag)

      if (ExitXSFlag .ne. 0 .or. MaxOtFlag .ne. 0 .or. ierr .ne. 0) return

      !***Modelled spectrum
      reflectance_hi = rint_fine

      !*** Fluorescence emission
      if (win_ini%Fsflag > 0) then
         !*** Add fluorescence
         do k = 1, win_ini%Fsflag
            reflectance_hi(:, 1) = reflectance_hi(:, 1) + &
                                   win%Fs(k)*(wavelength_hi(:) - lambda0)**(k - 1)* &
                                   (A1*exp(-(wavelength_hi(:) - lambda1)**2.d0/(2.d0*sigma1**2)) + &
                                    A2*exp(-(wavelength_hi(:) - lambda2)**2.d0/(2.d0*sigma2**2)))* &
                                   exp(-taua_tot(:)/cos(dble(iza)/180.*pi))/win%sun_spectrum_ref_hi(:)
         end do
      elseif (win_ini%Fsflag == -1) then
         !*** Add fluorescence
         do k = 1, abs(win_ini%Fsflag)
            reflectance_hi(:, 1) = reflectance_hi(:, 1) + &
                                   win%Fs(k)*(wavelength_hi(:) - lambda0)**(k - 1)* &
                                   (A1*exp(-(wavelength_hi(:) - lambda1)**2.d0/(2.d0*sigma1**2)) + &
                                    A2*exp(-(wavelength_hi(:) - lambda2)**2.d0/(2.d0*sigma2**2)))* &
                                   exp(-taua_tot(:)/cos(dble(iza)/180.*pi))/win%sun_spectrum_ref_hi(:)
         end do
      end if

      !***Fluorescence derivative
      deriv_hi%Fs = 0.D0
      if (win_ini%Fsflag > 0) then
         do k = 1, win_ini%Fsflag
            deriv_hi%Fs(:, k) = 1.d0*(wavelength_hi(:) - lambda0)**(k - 1)* &
                                (A1*exp(-(wavelength_hi(:) - lambda1)**2.d0/(2.d0*sigma1**2)) + &
                                 A2*exp(-(wavelength_hi(:) - lambda2)**2.d0/(2.d0*sigma2**2)))* &
                                exp(-taua_tot(:)/cos(dble(iza)/180.*pi))/win%sun_spectrum_ref_hi(:)
            do l = 1, win_ini%nwave_hi
               if (wavelength_hi(l) > lambda_o2lo .and. wavelength_hi(l) < lambda_o2up) then
                  deriv_hi%Fs(l, k) = 0.d0
               end if
            end do
         end do
      end if

      !***Absorber derivatives
      deriv_hi%densmol = 0.D0
      l = nrt/nlay
      l1 = natm/nrt
      l2 = natm/nlay
      do j = 1, absorb%ntype_target
         do i = 1, win_ini%ntype
            if (win_ini%type_x(i) == absorb%type_x_target(j)) then
               do n = 0, nrt - 1
                  m = int(n/l) + 1
                  n1 = n*l1 + 1
                  n2 = (n + 1)*l1
                  n3 = (m - 1)*l2 + 1
                  n4 = m*l2
                  if (sum(win%x_molec(n3:n4, i)) > 1.d-12) then
                     w = sum(win%x_molec(n1:n2, i))/sum(win%x_molec(n3:n4, i))
                  else
                     w = nlay/nrt
                  end if
!                deriv_hi%densmol(:, m, i, 1) = deriv_hi%densmol(:, m, i, 1) + deriv_rt%densmol(:, n+1, i, 1)*w
                  deriv_hi%densmol(:, m, i, :) = deriv_hi%densmol(:, m, i, :) + deriv_rt%densmol(:, n + 1, i, :)*w
               end do
            end if
         end do
      end do

      do j = 1, absorb%ntype_global
         do i = 1, win_ini%ntype
            if (win_ini%type_x(i) == absorb%type_x_global(j)) then
               do n = 0, nrt - 1
                  n1 = n*l1 + 1
                  n2 = (n + 1)*l1
                  if (sum(win%x_molec(:, i)) > 1.d-12) then
                     w = sum(win%x_molec(n1:n2, i))/sum(win%x_molec(:, i))
                  else
                     w = 1/nrt
                  end if
!                deriv_hi%densmol(:, 1, i, 1) = deriv_hi%densmol(:, 1, i, 1) + (deriv_rt%densmol(:, n+1, i, 1))*w
                  deriv_hi%densmol(:, 1, i, :) = deriv_hi%densmol(:, 1, i, :) + (deriv_rt%densmol(:, n + 1, i, :))*w
               end do
               forall (l=2:nlay) deriv_hi%densmol(:, l, i, :) = deriv_hi%densmol(:, 1, i, :)
            end if
         end do
      end do

      !***Albedo derivatives
      do k = 1, win_ini%albflag
!       deriv_hi%alb(:, k) = (deriv_rt%alb(:, 1))* &
         forall (nsc=1:nstokes) deriv_hi%alb(:, k, nsc) = (deriv_rt%alb(:, 1, nsc))* &
                                                          ((wavelength_hi(:) - wavelength_hi(1))**(k - 1))
      end do

      ntype_aer = size(aerosol)
      naer = 0
      do n = 1, ntype_aer
         naer = naer + sum(aerosol(n)%AerosolFlags)
      end do
      if (naer > 0) then
         !*** Aerosol derivatives
         deriv_hi%aerosol = 0.D0
         off = 0
         do k = 1, ntype_aer
            if (aerosol(k)%AerosolFlags(1) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 1, k, :) ! reff / power
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(2) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 2, k, :) ! veff
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(3) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 3, k, :) ! real refr index
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(4) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 4, k, :) ! im. refr. index
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(5) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 5, k, :) ! aer_col
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(6) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 6, k, :) ! shapefrac
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(7) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 7, k, :) !  height param.1 (central height)
               off = off + 1
            end if
            if (aerosol(k)%AerosolFlags(8) == 1) then
               deriv_hi%aerosol(:, off + 1, 1, :) = deriv_rt%aerosol(:, 8, k, :) !  height param.2 (width)
               off = off + 1
            end if
         end do
      end if !naer>0

      !***Temperature derivative
      deriv_hi%T = 0.D0
      if (flag%temp == 1) then
         deriv_hi%T = deriv_rt%T
      end if

      !***Pressure/O2 column derivatives
      deriv_hi%P = 0.D0
      if (flag%O2 == 1) then
         l = nrt/nlay
         l1 = natm/nrt
         l2 = natm/nlay
         do j = 1, absorb%ntype_target
            if (absorb%type_x_target(j) == 7) then
               do n = 0, nrt - 1
                  m = int(n/l) + 1
                  n1 = n*l1 + 1
                  n2 = (n + 1)*l1
                  n3 = (m - 1)*l2 + 1
                  n4 = m*l2
                  w = sum(dvair(n1:n2))/sum(dvair(n3:n4))
                  deriv_hi%P(:, m, :) = deriv_hi%P(:, m, :) + (deriv_rt%P(:, n + 1, :))*w
               end do
            end if
         end do
         do j = 1, absorb%ntype_global
            if (absorb%type_x_global(j) == 7) then
               do n = 0, nrt - 1
                  n1 = n*l1 + 1
                  n2 = (n + 1)*l1
                  w = sum(dvair(n1:n2))/sum(dvair(:))
                  deriv_hi%P(:, 1, :) = deriv_hi%P(:, 1, :) + (deriv_rt%P(:, n + 1, :))*w
               end do
               forall (l=2:nlay) deriv_hi%P(:, l, :) = deriv_hi%P(:, 1, :)
            end if
         end do
         do i = 1, win_ini%ntype
            if (win_ini%type_x(i) == 7) then
               deriv_hi%P(:, :, :) = deriv_hi%P(:, :, :) + deriv_hi%densmol(:, :, i, :)
            end if
         end do
      end if

   end subroutine forward_model_hi

   !------------------------------------------------------------------------------
   !> @brief Compute low-resolution model reflectance and derivatives, convolved with ILS
   !------------------------------------------------------------------------------
   subroutine forward_model_lo( &
      flag, &
      naer, nlay, absorb, &
      response, &
      win_ini, reflectance_hi, radiance_lo, win, deriv_hi, &
      ierr)
      !*** Input
      type(settings_flags), intent(in) :: flag
      integer, intent(in) :: naer, nlay
      type(absorbers), intent(in) :: absorb
      type(instrument_response), intent(in) :: response
      type(window_ini), intent(in) :: win_ini
      real(double), dimension(:, :), intent(in) :: reflectance_hi
      !*** Input/output
      type(derivatives), intent(inout) :: deriv_hi
      type(window_spectrum), intent(inout) :: win
      !*** Output
      real(double), dimension(:), intent(out) :: radiance_lo
      integer, intent(out) :: ierr
      !*** local variables
      integer :: i, j, k, l, nwave_hi, nwave_lo, nils
      real(double) :: dwave, s1, s2, s3
      real(double), dimension(:), allocatable :: wavelength_lo_per, wavelength_hi_per
      real(double), dimension(:), allocatable :: radiance_lo_per, radiance_hi_per
      real(double), dimension(:), allocatable :: radiance_hi, sun_spectrum_hi_per, sun_spectrum_lo_per
      integer, dimension(:), allocatable :: is_temp, ie_temp
      real(double), dimension(:, :), allocatable :: resp_temp
      character(stringlen) :: message
      !----------------------------------------------------------------

      !*** Initialize
      ierr = 0

      !*** Get dimensions
      nwave_hi = size(reflectance_hi(:, 1))
      nwave_lo = size(radiance_lo)
      nils = size(win%resp_store(1, :))
      !*** Allocate local allocatable arrays
      allocate (wavelength_lo_per(nwave_lo), &
                wavelength_hi_per(nwave_hi), &
                radiance_lo_per(nwave_lo), &
                radiance_hi_per(nwave_hi), &
                radiance_hi(nwave_hi), &
                sun_spectrum_hi_per(nwave_hi), &
                sun_spectrum_lo_per(nwave_lo), &
                stat=ierr)
      if (nstokes > 1) then
         !***Assign multiplicative stokes coefficients
         s1 = win%measurement_stokesc(1)
         s2 = win%measurement_stokesc(2)
         s3 = win%measurement_stokesc(3)
      else
         s1 = 1.d0
         s2 = 0.d0
         s3 = 0.d0
      end if
      if (ierr .ne. 0) then
         write (message, *) 'FORWARD_MODEL_LO: memory allocation error'
         ierr = ierr_all
         goto 999
      end if
      !*** Spectral response is calculated for radiance (not reflectance)
      if (nstokes == 1) then
         radiance_hi = reflectance_hi(:, 1)*win%sun_spectrum_ref_hi
      else
         radiance_hi = (s1*reflectance_hi(:, 1) + s2*reflectance_hi(:, 2) + s3*reflectance_hi(:, 3))*win%sun_spectrum_ref_hi
      end if

      !*** Calculate low-resolution spectrum
      if (flag%ils == 1) then
         allocate (resp_temp(nwave_lo, nils), &
                   is_temp(nwave_lo), &
                   ie_temp(nwave_lo), &
                   stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'FORWARD_MODEL_LO: memory allocation error'
            ierr = ierr_all
            goto 999
         end if
         !*** Recalculate and store response function if need be
         if (abs(win_ini%spsh0flag) + abs(win_ini%spsh1flag) + abs(win_ini%spsh2flag) /= 0) then
            call response_internal( &
               win_ini%wavelength_hi, &   ! high-resolution wavelength grid of model
               win%wavelength_lo_new, &   ! low-resolution wavelength grid of measurement
               response%ils_dwave, &      ! delta(wavelength) grid of instrument line shape
               response%resp_store, &     ! instrument line function
               win%resp_store, &          ! array for storing instrument response
               win%ie_store, &            ! array for storing instrument response
               win%is_store, &            ! array for storing instrument response
               ierr)                      ! error identifier
            if (ierr .ne. 0) return
         end if
         !*** Use stored response function
         call spectral_response_stored( &
            win%resp_store, &
            win%ie_store, &
            win%is_store, &
            radiance_hi, &
            radiance_lo)
      elseif (flag%ils == 2) then
         call spectral_response_custom( &
            win_ini%wavelength_hi, radiance_hi, win_ini%nwave_hi, &
            win%wavelength_lo_new, radiance_lo, nwave_lo, &
            response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
         if (ierr .ne. 0) return
      end if

      !HH: moved from profile_inversion to here
      !*** Convolve model solar spectrum by shifted instrument response function if need be:
      if (abs(win_ini%spsh0flag) + abs(win_ini%spsh1flag) + abs(win_ini%spsh2flag) + abs(win_ini%sunsh0flag) /= 0) then
         if (flag%ils == 1) then
            call spectral_response_stored( &
               win%resp_store, &
               win%ie_store, &
               win%is_store, &
               win%sun_spectrum_ref_hi, &
               win%sun_spectrum_ref_lo)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
               win%wavelength_lo_new, win%sun_spectrum_ref_lo, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
      end if

      !*** Shift+stretch/squeeze of the spectrum wrt sun reference
      !*** Shift
      if (win_ini%spsh0flag == 1) then      !wavelength shift [nm]
         dwave = 1.d-9
         wavelength_lo_per = win%wavelength_lo_new + dwave
         if (flag%ils == 1) then
            !*** Recalculate response function, but store it only temporarily
            call response_internal( &
               win_ini%wavelength_hi, &
               wavelength_lo_per, &
               response%ils_dwave, &
               response%resp_store, &
               resp_temp, &
               ie_temp, &
               is_temp, &
               ierr)
            if (ierr .ne. 0) return
            !*** Use temporarily response function
            call spectral_response_stored( &
               resp_temp, &
               ie_temp, &
               is_temp, &
               radiance_hi, &
               radiance_lo_per)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, radiance_hi, win_ini%nwave_hi, &
               wavelength_lo_per, radiance_lo_per, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
         do k = 1, win%nwave_lo
            win%deriv_specshift0_lo(k) = (radiance_lo_per(k) - radiance_lo(k))/dwave
         end do
      end if

      !*** Stretch/squeeze 1st order
      if (win_ini%spsh1flag == 1) then
         dwave = 1.d-11
         do k = 1, nwave_lo
            wavelength_lo_per(k) = win%wavelength_lo_new(k) + &
                                   dwave*(win%wavelength_lo_new(k) - win%wavelength_lo_new(1))
         end do
         if (flag%ils == 1) then
            !*** Recalculate response function, but store it only temporarily
            call response_internal( &
               win_ini%wavelength_hi, &
               wavelength_lo_per, &
               response%ils_dwave, &
               response%resp_store, &
               resp_temp, &
               ie_temp, &
               is_temp, &
               ierr)
            if (ierr .ne. 0) return
            !*** Use temporarily stored response function
            call spectral_response_stored( &
               resp_temp, &
               ie_temp, &
               is_temp, &
               radiance_hi, &
               radiance_lo_per)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, radiance_hi, win_ini%nwave_hi, &
               win%wavelength_lo_new, radiance_lo_per, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
         do k = 1, nwave_lo
            win%deriv_specshift1_lo(k) = (radiance_lo_per(k) - radiance_lo(k))/dwave
         end do
      end if

      !*** Stretch/squeeze 2nd order
      if (win_ini%spsh2flag == 1) then
         dwave = 1.d-13
         do k = 1, win%nwave_lo
            wavelength_lo_per(k) = win%wavelength_lo_new(k) + &
                                   dwave*(win%wavelength_lo_new(k) - win%wavelength_lo_new(1))**2
         end do
         if (flag%ils == 1) then
            !*** Recalculate response function, but store it only temporarily
            call response_internal( &
               win_ini%wavelength_hi, &
               wavelength_lo_per, &
               response%ils_dwave, &
               response%resp_store, &
               resp_temp, &
               ie_temp, &
               is_temp, &
               ierr)
            if (ierr .ne. 0) return
            !*** Use temporarily stored response function
            call spectral_response_stored( &
               resp_temp, &
               ie_temp, &
               is_temp, &
               radiance_hi, &
               radiance_lo_per)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, radiance_hi, win_ini%nwave_hi, &
               win%wavelength_lo_new, radiance_lo_per, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
         do k = 1, win%nwave_lo
            win%deriv_specshift2_lo(k) = (radiance_lo_per(k) - radiance_lo(k))/dwave
         end do
      end if

      !*** Shift of the sun wrt cross-sections
      if (win_ini%sunsh0flag == 1) then      !wavelength shift
         dwave = 1.d-10
         wavelength_hi_per = win%wavelength_hi_new + dwave
         !*** Interpolate solar spectrum on shifted wavelength grid
         call spline_interpol(wavelength_hi_per, win%sun_spectrum_sat_hi, nwave_hi, &
                              win_ini%wavelength_hi, sun_spectrum_hi_per, nwave_hi, ierr)
         if (ierr .ne. 0) then
            write (message, *) 'FORWARD_MODEL_LO.SPLINE_INTERPOL.SPLINT: bad input'
            ierr = ierr_intrpl
            goto 999
         end if
         !*** Calculate radiance from reflectance with shifted solar spectrum
         if (nstokes == 1) then
            radiance_hi_per = reflectance_hi(:, 1)*sun_spectrum_hi_per
         else
            radiance_hi_per = (s1*reflectance_hi(:, 1) + s2*reflectance_hi(:, 2) + s3*reflectance_hi(:, 3))*sun_spectrum_hi_per
         end if

         if (flag%ils == 1) then
            !*** Convolve radiance
            call spectral_response_stored( &
               win%resp_store, &
               win%ie_store, &
               win%is_store, &
               radiance_hi_per, &
               radiance_lo_per)
            !*** Convolve shifted solar spectrum by instrument response function:
            call spectral_response_stored( &
               win%resp_store, &
               win%ie_store, &
               win%is_store, &
               sun_spectrum_hi_per, &
               sun_spectrum_lo_per)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, radiance_hi_per, win_ini%nwave_hi, &
               win%wavelength_lo_new, radiance_lo_per, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
            call spectral_response_custom( &
               win_ini%wavelength_hi, sun_spectrum_hi_per, win_ini%nwave_hi, &
               wavelength_lo_per, sun_spectrum_lo_per, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
         if (flag%fit == 1) then
            do k = 1, nwave_lo
               win%deriv_sunshift0_lo(k) = &
                  (radiance_lo_per(k)/sun_spectrum_lo_per(k) - radiance_lo(k)/win%sun_spectrum_ref_lo(k))/dwave
            end do
         elseif (flag%fit == 2) then
            do k = 1, nwave_lo
               win%deriv_sunshift0_lo(k) = (radiance_lo_per(k) - radiance_lo(k))/dwave
            end do
         end if
      end if

      !*** Convolve hi-res derivatives
      do j = 1, absorb%ntype_target
         do i = 1, win_ini%ntype
            if (win_ini%type_x(i) == absorb%type_x_target(j)) then
               do k = 1, nlay
                  if (nstokes > 1) then
                     deriv_hi%densmol(:, k, i, 1) = (s1*deriv_hi%densmol(:, k, i, 1) + s2*deriv_hi%densmol(:, k, i, 2) + &
                                                     s3*deriv_hi%densmol(:, k, i, 3))
                  end if
                  if (flag%ils == 1) then
                     call spectral_response_stored( &
                        win%resp_store, &
                        win%ie_store, &
                        win%is_store, &
                        deriv_hi%densmol(:, k, i, 1)*win%sun_spectrum_ref_hi, &
                        win%derivatives_lo(:, k, i))
                  elseif (flag%ils == 2) then
                     call spectral_response_custom( &
                        win_ini%wavelength_hi, deriv_hi%densmol(:, k, i, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                        win%wavelength_lo_new, win%derivatives_lo(:, k, i), nwave_lo, &
                        response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
                     if (ierr .ne. 0) return
                  end if
               end do
            end if
         end do
      end do
      do j = 1, absorb%ntype_global
         do i = 1, win_ini%ntype
            if (win_ini%type_x(i) == absorb%type_x_global(j)) then
               if (nstokes > 1) then
                  deriv_hi%densmol(:, 1, i, 1) = (s1*deriv_hi%densmol(:, 1, i, 1) + s2*deriv_hi%densmol(:, 1, i, 2) + &
                                                  s3*deriv_hi%densmol(:, 1, i, 3))
               end if
               if (flag%ils == 1) then
                  call spectral_response_stored( &
                     win%resp_store, &
                     win%ie_store, &
                     win%is_store, &
                     deriv_hi%densmol(:, 1, i, 1)*win%sun_spectrum_ref_hi, &
                     win%derivatives_lo(:, 1, i))
               elseif (flag%ils == 2) then
                  call spectral_response_custom( &
                     win_ini%wavelength_hi, deriv_hi%densmol(:, 1, i, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                     win%wavelength_lo_new, win%derivatives_lo(:, 1, i), nwave_lo, &
                     response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
                  if (ierr .ne. 0) return
               end if
               forall (l=2:nlay) win%derivatives_lo(:, l, i) = win%derivatives_lo(:, 1, i)

            end if
         end do
      end do

      if (flag%O2 == 1) then
         do k = 1, nlay
            if (nstokes > 1) then
               deriv_hi%P(:, k, 1) = (s1*deriv_hi%P(:, k, 1) + s2*deriv_hi%P(:, k, 2) + &
                                      s3*deriv_hi%P(:, k, 3))
            end if

            if (flag%ils == 1) then
               call spectral_response_stored( &
                  win%resp_store, &
                  win%ie_store, &
                  win%is_store, &
                  deriv_hi%P(:, k, 1)*win%sun_spectrum_ref_hi, &
                  win%derivP_lo(:, k))
            elseif (flag%ils == 2) then
               call spectral_response_custom( &
                  win_ini%wavelength_hi, deriv_hi%P(:, k, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                  win%wavelength_lo_new, win%derivP_lo(:, k), nwave_lo, &
                  response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
               if (ierr .ne. 0) return
            end if
         end do
      end if

      if (flag%temp == 1) then
         if (nstokes > 1) then
            deriv_hi%T(:, 1) = (s1*deriv_hi%T(:, 1) + s2*deriv_hi%T(:, 2) + &
                                s3*deriv_hi%T(:, 3))
         end if
         if (flag%ils == 1) then
            call spectral_response_stored( &
               win%resp_store, &
               win%ie_store, &
               win%is_store, &
               deriv_hi%T(:, 1)*win%sun_spectrum_ref_hi, &
               win%derivT_lo)
         elseif (flag%ils == 2) then
            call spectral_response_custom( &
               win_ini%wavelength_hi, deriv_hi%T(:, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
               win%wavelength_lo_new, win%derivT_lo, nwave_lo, &
               response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
            if (ierr .ne. 0) return
         end if
      end if

      if (naer .gt. 0) then
         do k = 1, naer
            if (nstokes > 1) then
               deriv_hi%aerosol(:, k, 1, 1) = (s1*deriv_hi%aerosol(:, k, 1, 1) + s2*deriv_hi%aerosol(:, k, 1, 2) + &
                                               s3*deriv_hi%aerosol(:, k, 1, 3))
            end if
            if (flag%ils == 1) then
               call spectral_response_stored( &
                  win%resp_store, &
                  win%ie_store, &
                  win%is_store, &
                  deriv_hi%aerosol(:, k, 1, 1)*win%sun_spectrum_ref_hi, &
                  win%deriv_aerosol_lo(:, k))
            elseif (flag%ils == 2) then
               call spectral_response_custom( &
                  win_ini%wavelength_hi, deriv_hi%aerosol(:, k, 1, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                  win%wavelength_lo_new, win%deriv_aerosol_lo(:, k), nwave_lo, &
                  response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
               if (ierr .ne. 0) return
            end if

         end do
      end if

      if (win_ini%albflag > 0) then
         do k = 1, win_ini%albflag
            if (nstokes > 1) then
               deriv_hi%alb(:, k, 1) = (s1*deriv_hi%alb(:, k, 1) + s2*deriv_hi%alb(:, k, 2) + s3*deriv_hi%alb(:, k, 3))
            end if
            if (flag%ils == 1) then
               call spectral_response_stored( &
                  win%resp_store, &
                  win%ie_store, &
                  win%is_store, &
                  deriv_hi%alb(:, k, 1)*win%sun_spectrum_ref_hi, &
                  win%deriv_albedo_lo(:, k))
            elseif (flag%ils == 2) then
               call spectral_response_custom( &
                  win_ini%wavelength_hi, deriv_hi%alb(:, k, 1)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                  win%wavelength_lo_new, win%deriv_albedo_lo(:, k), nwave_lo, &
                  response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
               if (ierr .ne. 0) return
            end if
         end do
      end if

      !***Intensity offset - Additive (win_ini%IOffFlag>0)  or Multiplicative (win_ini%IOffFlag< 0)
      if (win_ini%IOffFlag > 0) then  !Additive (win_ini%IOffFlag>0)
         do k = 1, win_ini%IOffFlag
            win%deriv_Ioff_lo(:, k) = ((win%wavelength_lo_new(:) - &
                                        win%wavelength_lo_new(1))**(k - 1))
            do i = 1, nwave_lo
               radiance_lo(i) = radiance_lo(i) + &
                                win%IOff(k)*((win%wavelength_lo_new(i) - win%wavelength_lo_new(1))**(k - 1))
            end do
         end do
      elseif (win_ini%IOffFlag == -1) then !Multiplicative (win_ini%IOffFlag<0)
         do k = 1, abs(win_ini%IOffFlag)
            win%deriv_Ioff_lo(:, k) = radiance_lo(i)
            do i = 1, nwave_lo
               radiance_lo(i) = radiance_lo(i)*win%IOff(k)
            end do
         end do
      elseif (win_ini%IOffFlag == -2) then
         k = 1
         win%deriv_Ioff_lo(:, k) = radiance_lo(i)
         do i = 1, nwave_lo
            radiance_lo(i) = radiance_lo(i)*win%IOff(k)
         end do
         k = 2
         win%deriv_Ioff_lo(:, k) = ((win%wavelength_lo_new(:) - win%wavelength_lo_new(1))**(k - 1))
         do i = 1, nwave_lo
            radiance_lo(i) = radiance_lo(i) + win%IOff(k)*((win%wavelength_lo_new(i) - win%wavelength_lo_new(1))**(k - 1))
         end do

      end if

      if (win_ini%Fsflag > 0) then
         do k = 1, win_ini%Fsflag
            if (nstokes > 1) then
               deriv_hi%Fs(:, k) = s1*deriv_hi%Fs(:, k)
            end if
            if (flag%ils == 1) then
               call spectral_response_stored( &
                  win%resp_store, &
                  win%ie_store, &
                  win%is_store, &
                  deriv_hi%Fs(:, k)*win%sun_spectrum_ref_hi, &
                  win%deriv_Fs_lo(:, k))
            elseif (flag%ils == 2) then
               call spectral_response_custom( &
                  win_ini%wavelength_hi, deriv_hi%Fs(:, k)*win%sun_spectrum_ref_hi, win_ini%nwave_hi, &
                  win%wavelength_lo_new, win%deriv_Fs_lo(:, k), nwave_lo, &
                  response%ils_dwave(1, :), response%resp_store(1, :), response%nils, ierr)
               if (ierr .ne. 0) return
            end if
         end do
      end if

      if (flag%fit == 1) then !fit reflectance
         do l = 1, nwave_lo
            radiance_lo(l) = radiance_lo(l)/win%sun_spectrum_ref_lo(l)
            if (abs(win_ini%spsh0flag) == 1) win%deriv_specshift0_lo(l) = win%deriv_specshift0_lo(l)/win%sun_spectrum_ref_lo(l)
            if (abs(win_ini%spsh1flag) == 1) win%deriv_specshift1_lo(l) = win%deriv_specshift1_lo(l)/win%sun_spectrum_ref_lo(l)
            if (abs(win_ini%spsh2flag) == 1) win%deriv_specshift2_lo(l) = win%deriv_specshift2_lo(l)/win%sun_spectrum_ref_lo(l)
            if (flag%temp == 1) win%derivT_lo(l) = win%derivT_lo(l)/win%sun_spectrum_ref_lo(l)
            if (win_ini%albflag > 0) win%deriv_albedo_lo(l, :) = win%deriv_albedo_lo(l, :)/win%sun_spectrum_ref_lo(l)
            if (naer .gt. 0) win%deriv_aerosol_lo(l, :) = win%deriv_aerosol_lo(l, :)/win%sun_spectrum_ref_lo(l)
            if (win_ini%IOffFlag .ne. 0) win%deriv_Ioff_lo(l, :) = win%deriv_Ioff_lo(l, :)/win%sun_spectrum_ref_lo(l)
            if (win_ini%Fsflag > 0) win%deriv_Fs_lo(l, :) = win%deriv_Fs_lo(l, :)/win%sun_spectrum_ref_lo(l)
         end do

         do i = 1, win_ini%ntype
            do k = 1, nlay
               do l = 1, nwave_lo
                  win%derivatives_lo(l, k, i) = win%derivatives_lo(l, k, i)/win%sun_spectrum_ref_lo(l)
               end do
            end do
         end do

         if (flag%O2 == 1) then
            do k = 1, nlay
               do l = 1, nwave_lo
                  win%derivP_lo(l, k) = win%derivP_lo(l, k)/win%sun_spectrum_ref_lo(l)
               end do
            end do
         end if

      end if

      deallocate (wavelength_lo_per, wavelength_hi_per, &
                  radiance_lo_per, radiance_hi_per, &
                  radiance_hi, sun_spectrum_hi_per, sun_spectrum_lo_per, &
                  stat=ierr)
      if (ierr .ne. 0) then
         write (message, *) 'FORWARD_MODEL_LO: memory deallocation error'
         ierr = ierr_deall
         goto 999
      end if

      return

999   continue
      call stopretrieval(message)

   end subroutine forward_model_lo
!------------------------------------------------------------------------------

end module forward_model_module
