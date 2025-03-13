module spec_interface_create_module
   use header_module
   use read_settings_module, only: window_ini
   use forward_model_module, only: window_spectrum
   use atmosphere_internal_module, only: atmospheric_scenario
   use spectrum_internal_module, only: spectrum, instrument_interface
   use spectral_response_module, only: instrument_response, response_internal, spectral_response_stored
   use netcdf
   use auxiliary_routines_module, only: check
   implicit none
   private
   !*** types
   public :: instrument_response, spectrum

   !*** Procedures
   public ::  synthetic_interface, instrument_interface, output_l1b, output_l1b_nc, output_lut_nc_js, output_l1b_nc_js, output_l1b_nc_ls, synthetic_interface_init,synthetic_interface_close
   private :: spectral_response_create_gauss

   integer, private :: ncid_spec, ncid_isrf
   integer, dimension(:), allocatable :: fixed_nrow
   type(instrument_response), dimension(:, :), allocatable :: response_row

   type :: row_data
      real(double), dimension(:, :), allocatable :: wavelength
      integer, dimension(:, :), allocatable:: ie_store, is_store
      real(double), dimension(:, :, :), allocatable :: resp_store
      real(double), dimension(:, :), allocatable :: sun_spectrum_sat_lo
   end type row_data

   type(row_data), dimension(:), allocatable :: row

contains

   !------------------------------------------------------------------------------

   subroutine synthetic_interface_init(spectral_file, isrf_file)
      character(len=*), intent(in) :: spectral_file, isrf_file
      integer :: ierr

    call check(nf90_open("/nfs/TROPOMI/users/haili/S5P_orbits/INPUT/L1B/2015/"//trim(spectral_file)//".nc", nf90_nowrite, ncid_spec), ierr)

      call check(nf90_open(trim(isrf_file), nf90_nowrite, ncid_isrf), ierr)

   end subroutine synthetic_interface_init

   !------------------------------------------------------------------------------
   subroutine synthetic_interface_close(ierr)
      integer, intent(out) :: ierr

      call check(NF90_CLOSE(ncid_isrf), ierr)
      call check(NF90_CLOSE(ncid_spec), ierr)
   end subroutine synthetic_interface_close

   !------------------------------------------------------------------------------
   !> @details Set up measurement wavelength grid and ILS (inline calculation of Gaussian ILS)
  !! Get low-resolution solar spectrum
   !------------------------------------------------------------------------------
   subroutine synthetic_interface(ils_flag, spectral_flag, ilspath, win_ini, win)
      integer, intent(in) :: ils_flag, spectral_flag
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(window_spectrum), dimension(:), allocatable, intent(out) :: win
      character(len=*), intent(in) :: ilspath
      !*** local
      integer :: i, l, j, k, n, nils, nwin, ierr, nwave
      type(instrument_response), dimension(:), allocatable :: response, response_in
      type(instrument_response), dimension(:, :), allocatable :: response_tmp
      type(spectrum), dimension(:), allocatable :: measurement
      character*2 :: ch
      character(stringlen) :: filename, band(2), file(2)
      integer :: ncid, grpid1, grpid2, grpid3, dimid, varid, nrow, irow

      nwin = size(win_ini)
      if (allocated(win)) deallocate (win)
      allocate (win(nwin))
      allocate (measurement(nwin), response_in(nwin), response(nwin))
      allocate (fixed_nrow(nwin))
      allocate (row(nwin))

      !*** Set up measurement wavelength grid
      if (spectral_flag == 1) then
         !*** Calculate spectral grid from namelist configuration settings
         do n = 1, nwin
            !*** Set low-resolution wavelength grid
            win(n)%nwave_lo = (win_ini(n)%wave_stop - win_ini(n)%wave_start)/ &
                              (win_ini(n)%fwhm/win_ini(n)%samp) + 1
            allocate (win(n)%wavelength_lo(win(n)%nwave_lo))
            do k = 1, win(n)%nwave_lo
               win(n)%wavelength_lo(k) = win_ini(n)%wave_start + (k - 1)*(win_ini(n)%fwhm/win_ini(n)%samp)
            end do
            allocate (measurement(n)%wavelength(win(n)%nwave_lo))
            measurement(n)%nwave = win(n)%nwave_lo
            measurement(n)%wavelength = win(n)%wavelength_lo
            fixed_nrow(:) = 1
            allocate (row(n)%wavelength(fixed_nrow(n), win(n)%nwave_lo))
            do irow = 1, fixed_nrow(n)
               row(n)%wavelength(irow, :) = measurement(n)%wavelength(:)
            end do
         end do
      end if

      call get_isrf(win_ini, measurement, ils_flag, ilspath, response, ierr)

      do n = 1, nwin
         if (ils_flag > 0) then
            !*** Store instrument response in data type "win"
            nils = size(response(n)%resp_store(1, :))
            if (allocated(win(n)%resp_store)) deallocate (win(n)%resp_store)
            if (allocated(win(n)%is_store)) deallocate (win(n)%is_store)
            if (allocated(win(n)%ie_store)) deallocate (win(n)%ie_store)
            allocate (win(n)%resp_store(win(n)%nwave_lo, nils))
            allocate (win(n)%is_store(win(n)%nwave_lo))
            allocate (win(n)%ie_store(win(n)%nwave_lo))
            !*** Convert from measurement wavelength grid to internal grid
            call response_internal( &
               win_ini(n)%wavelength_hi, & ! high-resolution wavelength grid of model
               win(n)%wavelength_lo, &       ! low-resolution wavelength grid of measurement
               response(n)%ils_dwave, &      ! delta(wavelength) grid of instrument line shape
               response(n)%resp_store, &     ! instrument line function
               win(n)%resp_store, &          ! array for storing instrument response
               win(n)%ie_store, &            ! array for storing instrument response
               win(n)%is_store, &            ! array for storing instrument response
               ierr)

            !*** Convolve solar spectrum by instrument response function:
            allocate (win(n)%sun_spectrum_sat_lo(win(n)%nwave_lo))
            call spectral_response_stored( &
               win(n)%resp_store, &
               win(n)%ie_store, &
               win(n)%is_store, &
               win_ini(n)%sun_spectrum_ref_hi, &
               win(n)%sun_spectrum_sat_lo)
         end if
      end do !n

      return

999   call stopretrieval('SYNTHETIC_INTERFACE: fatal error')
      return

   end subroutine synthetic_interface

   !------------------------------------------------------------------------------

   subroutine get_isrf(win_ini, measurement, ils_flag, ilspath, response, ierr)
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      integer, intent(in) :: ils_flag
      character(len=*), intent(in) :: ilspath
      type(instrument_response), dimension(:), allocatable, intent(out) :: response
      integer, intent(out) :: ierr
      !*** local
      character*2 :: ch
      character(stringlen) :: filename
      integer :: n, nwin

      nwin = size(win_ini)
      allocate(response(nwin))

      if (ils_flag == 1) then
         !*** Inline calculation of Gaussian ISRF
         call calculate_isrf(win_ini, measurement, response)
         !*** Write ISRF to custom file
         do n = 1, nwin
            write (ch, '(i2.2)') n
            filename = trim(ilspath)//'isrf_'//ch//'.nc'
            call write_isrf_netcdf(trim(filename), response(n))
         end do
      else if (ils_flag == 2) then
         ! read ISRF from custom netCDF file
         call read_custom_isrf(win_ini, measurement, trim(ilspath), response, ierr)
         call interpolate_custom_isrf(win_ini, measurement, response)
         if (ierr .ne. 0) then
            write (*, *) "failed to read custom isrf"
            stop
         end if
      else
         write (*, *) "invalid ils_flag"
         stop
      end if
   end subroutine get_isrf

   !------------------------------------------------------------------------------

   subroutine write_isrf_netcdf(filename, response)
      character(len=*), intent(in) :: filename
      type(instrument_response), intent(in) :: response
      !*** local
      integer :: ncid, ierr, dimid1, dimid2, varid1, varid2, varid3, dimids(2), stat, i, j
      real(double), dimension(:, :), allocatable :: response_netcdf
      real(double), dimension(:), allocatable :: dw, wavelength

      ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
      ! overwrite this file, if it already exists.
      call check(nf90_create(trim(filename), nf90_clobber, ncid), stat)

      ! Define the dimensions. NetCDF will hand back an ID for each.
      ! Dimension size: number of wavelength differences
      call check(nf90_def_dim(ncid, "dwl", response%nils, dimid1), stat)
      ! Dimension size: number of measured wavelengths
      call check(nf90_def_dim(ncid, "wl_i", response%nwave, dimid2), stat)

      ! Define the variables
      call check(nf90_def_var(ncid, "Wavelength_differences", NF90_double, dimid1, varid1), stat)
      call check(nf90_def_var(ncid, "Measured_wavelengths", NF90_double, dimid2, varid2), stat)

      dimids = (/dimid1, dimid2/)
      ! Define the variable.
      call check(nf90_def_var(ncid, "Response", nf90_double, dimids, varid3), stat)

      ! End define mode. This tells netCDF we are done defining metadata.
      call check(nf90_enddef(ncid), stat)

      ! Allocate data fields of the ISRF
      allocate (dw(response%nils))
      allocate (wavelength(response%nwave), stat=ierr) ! Measured wavelengths at which the ISRF is defined ! {{{
      allocate (response_netcdf(response%nils, response%nwave), stat=ierr) ! Response function representative at different measured wavelengths as function of the wavelength difference. !

      do i = 1, response%nwave
         wavelength(i) = response%wavelength(i)
         do j = 1, response%nils
            dw(j) = response%ils_dwave(i, j)
            response_netcdf(j, i) = response%resp_store(i, j)
         end do
      end do

      ! Write the data to the file
      call check(nf90_put_var(ncid, varid1, dw), stat)

      call check(nf90_put_var(ncid, varid2, wavelength), stat)

      call check(nf90_put_var(ncid, varid3, response_netcdf), stat)

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      ! Mark success
      stat = 0

   end subroutine write_isrf_netcdf

   !------------------------------------------------------------------------------
   !> @details We compute a Gaussian instrument response function and store it
  !! in the data type "response" to be passed to the retrieval algorithm
  !! This routine is not used in the prototype, where we shall
  !! receive the instrument response in a tbd form.
   !------------------------------------------------------------------------------
   subroutine calculate_isrf(win_ini, measurement, response)
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      type(instrument_response), dimension(:), allocatable, intent(out) :: response
      !*** Local variables
      integer :: i, k, l, n, nils, nband
      real(double) :: fwhm
      real(double), dimension(:), allocatable :: ilswave                 ! ILS wavelength grid
      real(double), dimension(:), allocatable :: ilsfunction         ! ILS
      !---------------------------------------------------------

      allocate (response(size(measurement)))
      nband = size(response)
      do n = 1, nband
         !*** Set low-resolution wavelength grid
         response(n)%nwave = measurement(n)%nwave
         allocate (response(n)%wavelength(response(n)%nwave))
         do k = 1, response(n)%nwave
            response(n)%wavelength(k) = measurement(n)%wavelength(k)
         end do
      end do

      do n = 1, nband
         fwhm = win_ini(n)%fwhm
         nils = int(fwhm*win_ini(n)%wvbd/win_ini(n)%reso)
         if (modulo(nils, 2) == 0) nils = nils + 1   ! make sure number of Gaussian ILS points are uneven
         allocate (ilswave(nils), &
                   ilsfunction(nils))
         call spectral_response_create_gauss( &
            fwhm, win_ini(n)%reso, nils, ilswave, ilsfunction)

         !*** Store in an array(nils) for each measured wavelength as a function of
         !*** wavelength difference
         response(n)%nils = nils
         if (allocated(response(n)%resp_store)) deallocate (response(n)%resp_store)
         if (allocated(response(n)%ils_dwave)) deallocate (response(n)%ils_dwave)
         allocate (response(n)%resp_store(response(n)%nwave, nils))
         allocate (response(n)%ils_dwave(response(n)%nwave, nils))
         do l = 1, response(n)%nwave
            do i = 1, nils
               response(n)%resp_store(l, i) = ilsfunction(i)
               response(n)%ils_dwave(l, i) = ilswave(i)
            end do
         end do
         deallocate (ilswave, ilsfunction)
      end do !loop n

   end subroutine calculate_isrf

   !------------------------------------------------------------------------------
   !> @details Compute Gaussian instrument line shape (gauss) as a function of
  !! wavelength differences (wgauss)
   !------------------------------------------------------------------------------
   subroutine spectral_response_create_gauss( &
      fwhm, dlambda, ngauss, wgauss, gauss)
      integer, intent(IN) :: ngauss
      real(double), intent(IN) :: fwhm, dlambda
      real(double), intent(OUT), dimension(ngauss) :: wgauss, gauss
      !*** local variables
      integer :: k
      real(double) :: gaussint, center, b
      real(double), parameter :: sqrtln2 = 0.8325546

      b = fwhm/sqrtln2/2.
      gaussint = 0.
      gauss = 0.
      center = (int(ngauss/2) + 1)*dlambda
      do k = 1, ngauss
         wgauss(k) = k*dlambda - center
         gauss(k) = DEXP(-(wgauss(k)/b)*(wgauss(k)/b))
         !         gaussint = gaussint + gauss(k)
         gaussint = gaussint + gauss(k)*dlambda  !HH: normalization
      end do
      gauss = gauss/gaussint

   end subroutine spectral_response_create_gauss

   !------------------------------------------------------------------------------
   !> @details read custom isrf from netcdf file
   !------------------------------------------------------------------------------
   subroutine read_custom_isrf(win_ini, measurement, filename, response, ierr)
      type(spectrum), dimension(:), intent(in) :: measurement
      type(window_ini), dimension(:), intent(in) :: win_ini
      character(len=*), intent(in) :: filename
      type(instrument_response), dimension(:), allocatable, intent(out) :: response
      integer, intent(out) :: ierr
      ! local
      integer :: ncid, id, i, j
      real(double), dimension(:, :), allocatable :: resp
      real(double), dimension(:), allocatable :: dw, wavelength, ils_dwave
      real(double) :: ils_dwave_max
      integer :: grpid(3), band, nband, wave, nwave, varid, nils, current_band, win, nwin
      integer :: ierr_band_found  ! number of data bands in file that surrounded fit window (debug)

      nwin = size(win_ini)
      allocate(response(nwin))

      ! read ISRF from netCDF file
      call check(nf90_open(trim(filename), nf90_nowrite, ncid), ierr)

      do win = 1, nwin
         ! get correct band from file. It's wavelength range has to surround
         ! the internal wavelength grid. The netcdf group containing this band needs to be used.
         call check(nf90_inq_grps(ncid, nband, grpid), ierr)
         if (ierr .ne. 0) return

         ierr_band_found = 0

         do band = 1, nband
            ! check if this band contains a wavelength grid that surrounds the fit window
            call check(nf90_inq_dimid(grpid(band), "channel", varid), ierr)
            call check(nf90_inquire_dimension(grpid(band), varid, len=nwave), ierr)

            ! Get wavelengths of the current band and write them into dummy variable wavelength
            if(allocated(wavelength)) deallocate(wavelength)
            allocate(wavelength(nwave))
            call check(nf90_inq_varid(grpid(band), "wavelength_center", varid), ierr)
            call check(nf90_get_var(grpid(band), varid, wavelength), ierr)

            ! check if current band surrounds current fit window. If not, go to the next band
            ! wavelength that needs to be surrounded: measurement(win)%wavelength
            if (wavelength(1) > minval(measurement(win)%wavelength) .or. wavelength(nwave) < maxval(measurement(win)%wavelength)) then
               if (band == nband) then
                  if (ierr_band_found == 0) then
                     print*, "ERROR IN READ_CUSTOM_ISRF: No bands surrounded fit window."
                  else
                     print*, "ERROR IN READ_CUSTOM_ISRF: ", ierr_band_found, " band(s) surrounded fit window but none had sufficiently large ils_dwave grid."
                  end if
               end if
               cycle
            else
               ierr_band_found = ierr_band_found + 1
            end if

            ! check if this band contains a wavelength offset grid that is sufficiently large
            call check(nf90_inq_dimid(grpid(band), "d_channel", varid), ierr)
            call check(nf90_inquire_dimension(grpid(band), varid, len=nils), ierr)

            ! get wavelength offsets of the current band and write them into dummy variable ils_dwave
            if (allocated(ils_dwave)) deallocate(ils_dwave)
            allocate(ils_dwave(nils))
            call check(nf90_inq_varid(grpid(band), "wavelength_offset", varid), ierr)
            call check(nf90_get_var(grpid(band), varid, ils_dwave), ierr)

            ! check if current band has sufficiently large ils_dwave. If not, go to the next band
            ! wavelength offsets that have to be surrounded: win_ini(win)%fwhm * win_ini(win)%wvbd
            ils_dwave_max = win_ini(win)%fwhm * win_ini(win)%wvbd
            if (ils_dwave(1) > -ils_dwave_max .or. ils_dwave(nils) < ils_dwave_max) then
               cycle
            end if

            current_band = band
            exit
         end do  ! loop over band

         ! correct band found to be current_band
         ! we have nwave, nils, wavelength, and ils_dwave, write those into response
         ! also get the correct response

         ! wavelengths on which ils is defined
         response(win)%nwave = nwave
         allocate(response(win)%wavelength(nwave))
         response(win)%wavelength = wavelength

         ! wavelength offsets for which ils is defined
         response(win)%nils = nils
         allocate(response(win)%ils_dwave(nwave, nils))
         do wave = 1, nwave
            response(win)%ils_dwave(wave, :) = ils_dwave
         end do  ! loop over wave

         ! response of ils
         allocate(resp(nils, nwave))
         call check(nf90_inq_varid(grpid(current_band), "response", varid), ierr)
         call check(nf90_get_var(grpid(current_band), varid, resp), ierr)

         allocate(response(win)%resp_store(nwave, nils))
         do wave = 1, nwave
            response(win)%resp_store(wave, :) = resp(:, wave)
         end do  ! loop over wave
      end do ! loop over win

      ! close nc file
      call check(nf90_close(ncid), ierr)

      ! Mark success
      ierr = 0
      return
   end subroutine read_custom_isrf

   subroutine interpolate_custom_isrf(win_ini, measurement, response)
      !*** in
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      !*** inout
      type(instrument_response), dimension(:), intent(inout) :: response
      !*** local
      type(instrument_response), dimension(:), allocatable :: response_source  ! input
      type(instrument_response), dimension(:), allocatable :: response_tmp  ! first interpolation
      type(instrument_response), dimension(:), allocatable :: response_target  ! second interpolation
      integer :: win, nwin
      integer :: wave, ils
      integer :: nwave_source, nils_source
      integer :: nwave_tmp, nils_tmp
      integer :: nwave_target, nils_target
      integer :: ierr

      nwin = size(win_ini)

      response_source = response
      allocate(response_tmp(nwin))
      allocate(response_target(nwin))

      do win = 1, nwin
         ! source dimensions
         nwave_source = response_source(win)%nwave
         nils_source = response_source(win)%nils

         ! target dimensions after second interpolation (over wavelength offsets)
         nwave_target = measurement(win)%nwave
         nils_target = int(win_ini(win)%fwhm*win_ini(win)%wvbd/win_ini(win)%reso)
         if (modulo(nils_target, 2) == 0) nils_target = nils_target + 1  ! make sure number of Gaussian ILS points are uneven
         ! prepare target
         response_target(win)%nwave = nwave_target
         response_target(win)%nils = nils_target
         allocate(response_target(win)%wavelength(nwave_target))
         allocate(response_target(win)%ils_dwave(nwave_target, nils_target))
         allocate(response_target(win)%resp_store(nwave_target, nils_target))

         ! dimensions after first interpolation (over wavelength centers)
         nwave_tmp = nwave_target
         nils_tmp = nils_source
         ! prepare tmp
         response_tmp(win)%nwave = nwave_tmp
         response_tmp(win)%nils = nils_tmp
         allocate(response_tmp(win)%wavelength(nwave_tmp))
         allocate(response_tmp(win)%ils_dwave(nwave_tmp, nils_tmp))
         allocate(response_tmp(win)%resp_store(nwave_tmp, nils_tmp))

         ! get correct wavelength center grid
         response_tmp(win)%wavelength = measurement(win)%wavelength
         response_target(win)%wavelength = measurement(win)%wavelength

         ! for all wavelength centers, get correct dwave grid
         do ils = 1, nils_target
            response_target(win)%ils_dwave(:, ils) = ils * win_ini(win)%reso - int(nils_target/2 + 1) * win_ini(win)%reso
         end do  ! loop over ils

         ! first interpolation:
         ! for all wavelength offsets, interpolate onto the correct wavelength center grid
         do ils = 1, response_source(win)%nils
            call spline_interpol( &
               response_source(win)%wavelength, response_source(win)%ils_dwave(:, ils), response_source(win)%nwave, &
               response_tmp(win)%wavelength, response_tmp(win)%ils_dwave(:, ils), response_tmp(win)%nwave, &
               ierr &
            )
            call spline_interpol( &
               response_source(win)%wavelength, response_source(win)%resp_store(:, ils), response_source(win)%nwave, &
               response_tmp(win)%wavelength, response_tmp(win)%resp_store(:, ils), response_tmp(win)%nwave, &
               ierr &
            )
         end do  ! loop over ils

         ! second interpolation:
         ! for all wavelength centers, interpolate onto the correct wavelength offset grid
         do wave = 1, response_tmp(win)%nwave
            call spline_interpol( &
               response_tmp(win)%ils_dwave(wave, :), response_tmp(win)%resp_store(wave, :), response_tmp(win)%nils, &
               response_target(win)%ils_dwave(wave, :), response_target(win)%resp_store(wave, :), response_target(win)%nils, &
               ierr &
            )
         end do  ! loop over wave
      end do  ! loop over win

      response = response_target
   end subroutine interpolate_custom_isrf

   !------------------------------------------------------------------------------
   !> @details Write simulated L1B data to ascii file
   !------------------------------------------------------------------------------
   subroutine output_l1b(measurement, win_ini, meta, spectrum_file)
      !*** input
      type(spectrum), dimension(:), intent(in) :: measurement
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(metadata), intent(in) :: meta
      character(len=*), intent(in) :: spectrum_file
      !***  local
      integer :: io, nwin, n, k

      open (newunit(io), file=spectrum_file, form='formatted')
      !*** Write instrument info
      write (io, '(7I6)') meta%time(1), meta%time(2), &
         meta%time(3), meta%time(4), meta%time(5), meta%time(6), meta%time(7)
      do k = 1, 5
         write (io, '(2ES23.15)') meta%lat(k), meta%lon(k)
      end do
      write (io, '(ES23.15)') meta%surface_elevation
      write (io, '(4ES23.15)') meta%sza, meta%iza, meta%saz, meta%iaz
      nwin = size(measurement)
      !*** Write synthetic spectrum
      write (io, *) nwin
      do n = 1, nwin
         write (io, *) ''
         write (io, *) ''
         write (io, *) measurement(n)%nwave
         write (io, *) win_ini(n)%fwhm
         write (io, *) ''
         write (io, *) ''
         do k = 1, measurement(n)%nwave
            write (io, '(7ES23.15)') measurement(n)%wavelength(k), &   !wavelength [nm]
               measurement(n)%radiance(k), &       !spectrum with noise
               measurement(n)%radiance_noise(k), &  !noise
               measurement(n)%radiance_error(k), &  !error
               measurement(n)%irradiance(k), &
               measurement(n)%irradiance_noise(k), &
               measurement(n)%irradiance_error(k)
         end do
      end do

      close (io)

   end subroutine output_l1b
   !------------------------------------------------------------------------------
   !> @details Write simulated L1B data to netcdf file. Append spectra if file already exists.
   !------------------------------------------------------------------------------
   subroutine output_l1b_nc(measurement, win_ini, win, meta, spectrum_file, ipixel)
      !*** input
      type(spectrum), dimension(:), intent(in) :: measurement
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(window_spectrum), dimension(:), intent(in):: win
      type(metadata), intent(in) :: meta
      character(len=*), intent(in) :: spectrum_file
      integer, intent(in) :: ipixel
      !***  local
      integer :: ncid, ierr, nwin, n, k, stat, ngroup, dim_nobs, dim_wave, lastindex, dim_ntime, nwave
      integer :: sza_id, vza_id, saz_id, vaz_id, lat_id, lon_id, time_id, pixel_id
    integer, dimension(:), allocatable :: grpid, rad_id, radnoise_id, raderror_id, irrad_id, irradnoise_id, irraderror_id, wave_id, ot_id, cot_id, alb_id
      integer :: dimids(2), start(2)
      character*1 :: ch
      logical :: exst
      character(stringlen) :: group_name
      !date and time
      integer :: hr, min
      character(len=2) :: hour, minute, second
      character(len=8) :: date
      integer :: values(8)
      character(len=20) :: creation_date, validity_start, validity_stop, time

      stat = 0

      nwin = size(measurement)
      allocate (grpid(nwin), &
                wave_id(nwin), &
                rad_id(nwin), &
                radnoise_id(nwin), &
                raderror_id(nwin), &
                irrad_id(nwin), &
                irradnoise_id(nwin), &
                irraderror_id(nwin), &
                ot_id(nwin), &
                cot_id(nwin), &
                alb_id(nwin), stat=ierr)

      inquire (FILE=trim(spectrum_file), EXIST=exst)
      if (.not. exst) then
         !*** Create the netCDF file.
         call check(nf90_create(trim(spectrum_file), nf90_netcdf4, ncid), stat)

         !*** get creation date
         call date_and_time(DATE=date, VALUES=values)
         write (hour, "(i2.2)") values(5)
         write (minute, "(i2.2)") values(6)
         write (second, "(i2.2)") values(7)
         ! local time
         time = hour//minute//second
         creation_date = trim(date)//'T'//trim(time)

         !*** Define global attributes
         call check(nf90_put_att(ncid, nf90_global, "title", "Simulated spectra for CSS"), stat)
         call check(nf90_put_att(ncid, nf90_global, "institution", "DLR / SRON"), stat)
         call check(nf90_put_att(ncid, nf90_global, "source", "RemoTeC simulation package"), stat)
         call check(nf90_put_att(ncid, nf90_global, "date_created", creation_date), stat)

         !*** Define the dimensions. NetCDF will hand back an ID for each.

         !*** Unlimited dimension
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dim_nobs), stat)

!!$       !*** Geometry
!!$       call check( nf90_def_var(ncid, "solar_zenith_angle", NF90_double, dim_nobs, sza_id), stat)
!!$       call check( nf90_def_var(ncid, "viewing_zenith_angle", NF90_double, dim_nobs, vza_id), stat)
!!$       call check( nf90_def_var(ncid, "solar_azimuth_angle", NF90_double, dim_nobs, saz_id), stat)
!!$       call check( nf90_def_var(ncid, "viewing_azimuth_angle", NF90_double, dim_nobs, vaz_id), stat)
!!$       !*** Geodata
!!$       call check( nf90_def_var(ncid, "latitude", NF90_double, dim_nobs, lat_id), stat)
!!$       call check( nf90_def_var(ncid, "longitude", NF90_double, dim_nobs, lon_id), stat)
!!$       !*** Time
!!$       call check(nf90_def_dim(ncid, "ntime", 6, dim_ntime), stat)
!!$       dimids =  (/dim_ntime, dim_nobs/)
!!$       call check( nf90_def_var(ncid, "time", NF90_int, dimids, time_id), stat)

         !*** Pixel identifier
         call check(nf90_def_var(ncid, "pixelID", NF90_int, dim_nobs, pixel_id), stat)

         do n = 1, nwin
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** spectral
            call check(nf90_def_dim(grpid(n), "nwave", measurement(n)%nwave, dim_wave), stat)

            !*** Define the variables
            dimids = (/dim_wave, dim_nobs/)
            call check(nf90_def_var(grpid(n), "wavelength", NF90_double, dim_wave, wave_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance", NF90_double, dimids, rad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_noise", NF90_double, dimids, radnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_error", NF90_double, dimids, raderror_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance", NF90_double, dimids, irrad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_noise", NF90_double, dimids, irradnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_error", NF90_double, dimids, irraderror_id(n)), stat)

            call check(nf90_def_var(grpid(n), "aerosol_optical_thickness", NF90_double, dim_nobs, ot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "cirrus_optical_thickness", NF90_double, dim_nobs, cot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "surface_albedo", NF90_double, dim_nobs, alb_id(n)), stat)

            !*** Define attributes
            call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
            call check(nf90_put_att(grpid(n), rad_id(n), "unit", "mol s-1 m-2 nm-1"), stat)
            call check(nf90_put_att(grpid(n), radnoise_id(n), "unit", "mol s-1 m-2 nm-1"), stat)
            call check(nf90_put_att(grpid(n), raderror_id(n), "unit", "mol s-1 m-2 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irrad_id(n), "unit", "mol s-1 m-2 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irradnoise_id(n), "unit", "mol s-1 m-2 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irraderror_id(n), "unit", "mol s-1 m-2 nm-1"), stat)

            call check(nf90_put_att(grpid(n), ot_id(n), "description", "aerosol optical thickness used in simulation"), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "description", "cirrus optical thickness used in simulation"), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "description", "surface albedo used in simulation"), stat)

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)
         end do

         lastindex = 1
      else
         !*** Open the netCDF file and append
         call check(NF90_OPEN(trim(spectrum_file), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** get modification date
         call date_and_time(DATE=date, VALUES=values)
         write (hour, "(i2.2)") values(5)
         write (minute, "(i2.2)") values(6)
         write (second, "(i2.2)") values(7)
         ! local time
         time = hour//minute//second
         creation_date = trim(date)//'T'//trim(time)
         call check(nf90_put_att(ncid, nf90_global, "history", "Last modified on "//creation_date), stat)
         call check(NF90_INQ_VARID(ncid, "pixelID", pixel_id), stat)
         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            !*** Get variable ID's
            call check(NF90_INQ_DIMID(grpid(n), "nwave", dim_wave), stat)
            call check(NF90_INQUIRE_DIMENSION(grpid(n), dim_wave, len=nwave), stat)
            if (nwave .ne. measurement(n)%nwave) then
               call writelog("OUTPUT_L1B_NC: nwave is not equal to number of spectral channels", 8)
            end if

            call check(NF90_INQ_VARID(grpid(n), "radiance", rad_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "radiance_noise", radnoise_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "radiance_error", raderror_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance", irrad_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance_noise", irradnoise_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance_error", irraderror_id(n)), stat)

            call check(NF90_INQ_VARID(grpid(n), "aerosol_optical_thickness", ot_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "cirrus_optical_thickness", cot_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "surface_albedo", alb_id(n)), stat)
         end do

         !*** Get number of spectra already in file
         call check(NF90_INQ_DIMID(ncid, "nobs", dim_nobs), stat)
         call check(NF90_INQUIRE_DIMENSION(ncid, dim_nobs, len=lastindex), stat)
         lastindex = lastindex + 1
      end if

      ! Write the data to the file
      call check(nf90_put_var(ncid, pixel_id, ipixel, start=(/lastindex/)), stat)

      do n = 1, nwin
         start = (/1, lastindex/)
         call check(nf90_put_var(grpid(n), rad_id(n), measurement(n)%radiance, start=start), stat)
         call check(nf90_put_var(grpid(n), radnoise_id(n), measurement(n)%radiance_noise, start=start), stat)
         call check(nf90_put_var(grpid(n), raderror_id(n), measurement(n)%radiance_error, start=start), stat)
         call check(nf90_put_var(grpid(n), irrad_id(n), measurement(n)%irradiance, start=start), stat)
         call check(nf90_put_var(grpid(n), irradnoise_id(n), measurement(n)%irradiance_noise, start=start), stat)
         call check(nf90_put_var(grpid(n), irraderror_id(n), measurement(n)%irradiance_error, start=start), stat)

         call check(nf90_put_var(grpid(n), ot_id(n), win(n)%ot, start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), cot_id(n), win(n)%cot, start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), alb_id(n), win(n)%albedo(1), start=(/lastindex/)), stat)
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (stat .ne. 0) then
         call writelog('OUTPUT_L1B_NC: Error in writing spectra to netCDF file', 8)
      end if

   end subroutine output_l1b_nc
   !------------------------------------------------------------------------------

   !------------------------------------------------------------------------------
   !> @details Write simulated L1B data to netcdf file. Append spectra if file already exists.
   !------------------------------------------------------------------------------
   subroutine output_l1b_nc_js(measurement, meta, spectrum_file, index_info)
      !*** input
      type(spectrum), dimension(:), intent(in) :: measurement
      type(metadata), intent(in) :: meta
      character(len=*), intent(in) :: spectrum_file, index_info
      !***  local
      integer :: ncid, ierr, nwin, stat, ngroup, dimid_nobs, dimid_wave, lastindex, dimid_time, nwave
      integer :: sza_id, vza_id, saz_id, vaz_id, lat_id, lon_id, time_id, pixel_id, elev_id, x_id, y_id
      integer, dimension(:), allocatable :: grpid, rad_id, radnoise_id, raderror_id, irrad_id, irradnoise_id, irraderror_id, wave_id
      integer :: start(2), dimids_spec(2)
      integer :: i, n, sx, sy
      character*1 :: ch
      logical :: exst
      character(stringlen) :: group_name

      stat = 0

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      nwin = size(measurement)
      allocate (grpid(nwin), &
                wave_id(nwin), &
                rad_id(nwin), &
                radnoise_id(nwin), &
                raderror_id(nwin), &
                irrad_id(nwin), &
                irradnoise_id(nwin), &
                irraderror_id(nwin), stat=ierr)

      inquire (FILE=trim(spectrum_file), EXIST=exst)

      if (.not. exst) then

         !*** Create the netCDF file.
         call check(nf90_create(trim(spectrum_file), nf90_netcdf4, ncid), stat)

         !*** Unlimited dimension for number of observations
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dimid_nobs), stat)

         !*** Indexdata
         call check(nf90_def_var(ncid, "x", NF90_int, dimid_nobs, x_id), stat)
         call check(nf90_def_var(ncid, "y", NF90_int, dimid_nobs, y_id), stat)

         !*** Geometry
         call check(nf90_def_var(ncid, "sza", NF90_double, dimid_nobs, sza_id), stat)
         call check(nf90_def_var(ncid, "vza", NF90_double, dimid_nobs, vza_id), stat)
         call check(nf90_def_var(ncid, "saa", NF90_double, dimid_nobs, saz_id), stat)
         call check(nf90_def_var(ncid, "vaa", NF90_double, dimid_nobs, vaz_id), stat)

         call check(nf90_put_att(ncid, sza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, sza_id, "description", "Solar Zenith Angle"), stat)
         call check(nf90_put_att(ncid, vza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, vza_id, "description", "Viewing Zenith Angle"), stat)
         call check(nf90_put_att(ncid, saz_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, saz_id, "description", "Solar Azimuth Angle"), stat)
         call check(nf90_put_att(ncid, vaz_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, vaz_id, "description", "Viewing Azimuth Angle"), stat)

         !*** Geodata
         call check(nf90_def_var(ncid, "latitude", NF90_double, dimid_nobs, lat_id), stat)
         call check(nf90_def_var(ncid, "longitude", NF90_double, dimid_nobs, lon_id), stat)
         call check(nf90_def_var(ncid, "elevation", NF90_double, dimid_nobs, elev_id), stat)

         call check(nf90_put_att(ncid, lat_id, "unit", "degrees_north"), stat)
         call check(nf90_put_att(ncid, lat_id, "description", "Latitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lon_id, "unit", "degrees_east"), stat)
         call check(nf90_put_att(ncid, lon_id, "description", "Longitude at pixel center"), stat)
         call check(nf90_put_att(ncid, elev_id, "unit", "m"), stat)
         call check(nf90_put_att(ncid, elev_id, "description", "Surface altitude at pixel center"), stat)

         !*** Timedata
         call check(nf90_def_dim(ncid, "ntime", 6, dimid_time), stat)
         call check(nf90_def_var(ncid, "time", NF90_int, dimid_time, time_id), stat)
 call check(nf90_put_var(ncid, time_id, [meta%time(1), meta%time(2), meta%time(3), meta%time(4), meta%time(5), meta%time(6)]), stat)
         call check(nf90_put_att(ncid, time_id, "description", "Date and time as [YYYY,MM,DD,HOUR,MIN,SEC]"), stat)

         do n = 1, nwin
            nwave = measurement(n)%nwave
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** Spectral data
            call check(nf90_def_dim(grpid(n), "nwave", nwave, dimid_wave), stat)

            !*** Define the variables
            dimids_spec = (/dimid_wave, dimid_nobs/)

            call check(nf90_def_var(grpid(n), "wavelength", NF90_double, dimid_wave, wave_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance", NF90_double, dimids_spec, rad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_noise", NF90_double, dimids_spec, radnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_error", NF90_double, dimids_spec, raderror_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance", NF90_double, dimids_spec, irrad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_noise", NF90_double, dimids_spec, irradnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_error", NF90_double, dimids_spec, irraderror_id(n)), stat)

            !*** Define attributes
            call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
            call check(nf90_put_att(grpid(n), rad_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), radnoise_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), raderror_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irrad_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irradnoise_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irraderror_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)
         end do

         lastindex = 1

         !*** GET VARIABLE IDs
      else

         !*** Open the netCDF file and append
         call check(NF90_OPEN(trim(spectrum_file), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** Indexdata
         call check(nf90_inq_varid(ncid, "x", x_id), stat)
         call check(nf90_inq_varid(ncid, "y", y_id), stat)

         !*** Geometry
         call check(nf90_inq_varid(ncid, "sza", sza_id), stat)
         call check(nf90_inq_varid(ncid, "vza", vza_id), stat)
         call check(nf90_inq_varid(ncid, "saa", saz_id), stat)
         call check(nf90_inq_varid(ncid, "vaa", vaz_id), stat)

         !*** Geodata
         call check(nf90_inq_varid(ncid, "latitude", lat_id), stat)
         call check(nf90_inq_varid(ncid, "longitude", lon_id), stat)
         call check(nf90_inq_varid(ncid, "elevation", elev_id), stat)

         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            !*** Get variable ID's
            call check(NF90_INQ_DIMID(grpid(n), "nwave", dimid_wave), stat)
            call check(NF90_INQUIRE_DIMENSION(grpid(n), dimid_wave, len=nwave), stat)
            if (nwave .ne. measurement(n)%nwave) then
               call writelog("OUTPUT_L1B_NC_JS: nwave is not equal to number of spectral channels", 8)
            end if

            !*** Spectral data
            call check(NF90_INQ_VARID(grpid(n), "radiance", rad_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "radiance_noise", radnoise_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "radiance_error", raderror_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance", irrad_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance_noise", irradnoise_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "irradiance_error", irraderror_id(n)), stat)
         end do

         !*** Get number of spectra already in file
         call check(NF90_INQ_DIMID(ncid, "nobs", dimid_nobs), stat)
         call check(NF90_INQUIRE_DIMENSION(ncid, dimid_nobs, len=lastindex), stat)
         lastindex = lastindex + 1

      end if

      !*** WRITE DATA TO FILE

      !*** Indexdata
      call check(nf90_put_var(ncid, x_id, sx, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, y_id, sy, start=(/lastindex/)), stat)

      !*** Geometry
      call check(nf90_put_var(ncid, sza_id, meta%sza, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, vza_id, meta%iza, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, saz_id, meta%saz, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, vaz_id, meta%iaz, start=(/lastindex/)), stat)

      !*** Geodata
      call check(nf90_put_var(ncid, lat_id, meta%lat(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, lon_id, meta%lon(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, elev_id, meta%surface_elevation, start=(/lastindex/)), stat)

      !*** Spectral data
      start = (/1, lastindex/)
      do n = 1, nwin
         call check(nf90_put_var(grpid(n), rad_id(n), measurement(n)%radiance, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), radnoise_id(n), measurement(n)%radiance_noise, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), raderror_id(n), measurement(n)%radiance_error, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irrad_id(n), measurement(n)%irradiance, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irradnoise_id(n), measurement(n)%irradiance_noise, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irraderror_id(n), measurement(n)%irradiance_error, start=(/1, lastindex/)), stat)
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (stat .ne. 0) then
         call writelog('OUTPUT_L1B_NC: Error in writing spectra to netCDF file', 8)
      end if

   end subroutine output_l1b_nc_js
   !------------------------------------------------------------------------------



   !------------------------------------------------------------------------------
   !> @details Write simulated L1B data to netcdf file. Append spectra if file already exists.
   !------------------------------------------------------------------------------
   subroutine output_l1b_nc_ls(measurement, meta, spectrum_file, index_info)
      !*** input
      type(spectrum), dimension(:), intent(in) :: measurement
      type(metadata), intent(in) :: meta
      character(len=*), intent(in) :: spectrum_file, index_info
      !***  local
      integer :: ncid, ierr, nwin, stat, ngroup, dimid_nobs, dimid_wave, lastindex, dimid_time, nwave
      integer :: sza_id, vza_id, saa_id, vaa_id, lat_id, lon_id, time_id, pixel_id, surface_elevation_id, x_id, y_id
      integer, dimension(:), allocatable :: grpid, rad_id, radnoise_id, raderror_id, irrad_id, irradnoise_id, irraderror_id, wave_id
      integer :: start(2), dimids_spec(2)
      integer :: i, n, sx, sy
      character*1 :: ch
      logical :: exst
      character(stringlen) :: group_name

      stat = 0

      ! Extract index in x-dimension to be written
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be written
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      nwin = size(measurement)
      allocate (grpid(nwin), &
                wave_id(nwin), &
                rad_id(nwin), &
                radnoise_id(nwin), &
                raderror_id(nwin), &
                irrad_id(nwin), &
                irradnoise_id(nwin), &
                irraderror_id(nwin), stat=ierr)

      inquire (FILE=trim(spectrum_file), EXIST=exst)

      if (.not. exst) then

         !*** Create the netCDF file
         call check(nf90_create(trim(spectrum_file), nf90_netcdf4, ncid), stat)

         !*** Unlimited dimension for number of observations
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dimid_nobs), stat)

         !*** Indexdata
         call check(nf90_def_var(ncid, "x", nf90_int, dimid_nobs, x_id), stat)
         call check(nf90_def_var(ncid, "y", nf90_int, dimid_nobs, y_id), stat)

         !*** Geometry
         call check(nf90_def_var(ncid, "sza", nf90_double, dimid_nobs, sza_id), stat)
         call check(nf90_def_var(ncid, "vza", nf90_double, dimid_nobs, vza_id), stat)
         call check(nf90_def_var(ncid, "saa", nf90_double, dimid_nobs, saa_id), stat)
         call check(nf90_def_var(ncid, "vaa", nf90_double, dimid_nobs, vaa_id), stat)

         call check(nf90_put_att(ncid, sza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, sza_id, "description", "Solar Zenith Angle"), stat)
         call check(nf90_put_att(ncid, vza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, vza_id, "description", "Viewing Zenith Angle"), stat)
         call check(nf90_put_att(ncid, saa_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, saa_id, "description", "Solar Azimuth Angle"), stat)
         call check(nf90_put_att(ncid, vaa_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, vaa_id, "description", "Viewing Azimuth Angle"), stat)

         !*** Geodata
         call check(nf90_def_var(ncid, "latitude", nf90_double, dimid_nobs, lat_id), stat)
         call check(nf90_def_var(ncid, "longitude", nf90_double, dimid_nobs, lon_id), stat)
         call check(nf90_def_var(ncid, "surface_elevation", nf90_double, dimid_nobs, surface_elevation_id), stat)

         call check(nf90_put_att(ncid, lat_id, "unit", "degrees_north"), stat)
         call check(nf90_put_att(ncid, lat_id, "description", "Latitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lon_id, "unit", "degrees_east"), stat)
         call check(nf90_put_att(ncid, lon_id, "description", "Longitude at pixel center"), stat)
         call check(nf90_put_att(ncid, surface_elevation_id, "unit", "m"), stat)
         call check(nf90_put_att(ncid, surface_elevation_id, "description", "Surface altitude at pixel center"), stat)

         !*** Timedata
         call check(nf90_def_dim(ncid, "ntime", 6, dimid_time), stat)
         call check(nf90_def_var(ncid, "time", nf90_int, dimid_time, time_id), stat)
         call check(nf90_put_var(ncid, time_id, [meta%time(1), meta%time(2), meta%time(3), meta%time(4), meta%time(5), meta%time(6)]), stat)
         call check(nf90_put_att(ncid, time_id, "description", "Date and time as [YYYY,MM,DD,HOUR,MIN,SEC]"), stat)

         do n = 1, nwin
            nwave = measurement(n)%nwave
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** Spectral data
            call check(nf90_def_dim(grpid(n), "nwave", nwave, dimid_wave), stat)

            !*** Define the variables
            dimids_spec = (/dimid_wave, dimid_nobs/)

            call check(nf90_def_var(grpid(n), "wavelength", nf90_double, dimid_wave, wave_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance", nf90_double, dimids_spec, rad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_noise", nf90_double, dimids_spec, radnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "radiance_error", nf90_double, dimids_spec, raderror_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance", nf90_double, dimids_spec, irrad_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_noise", nf90_double, dimids_spec, irradnoise_id(n)), stat)
            call check(nf90_def_var(grpid(n), "irradiance_error", nf90_double, dimids_spec, irraderror_id(n)), stat)

            !*** Define attributes
            call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
            call check(nf90_put_att(grpid(n), rad_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), radnoise_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), raderror_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irrad_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irradnoise_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
            call check(nf90_put_att(grpid(n), irraderror_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)
         end do

         lastindex = 1

         !*** GET VARIABLE IDs
      else

         !*** Open the netCDF file and append
         call check(nf90_open(trim(spectrum_file), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** Indexdata
         call check(nf90_inq_varid(ncid, "x", x_id), stat)
         call check(nf90_inq_varid(ncid, "y", y_id), stat)

         !*** Geometry
         call check(nf90_inq_varid(ncid, "sza", sza_id), stat)
         call check(nf90_inq_varid(ncid, "vza", vza_id), stat)
         call check(nf90_inq_varid(ncid, "saa", saa_id), stat)
         call check(nf90_inq_varid(ncid, "vaa", vaa_id), stat)

         !*** Geodata
         call check(nf90_inq_varid(ncid, "latitude", lat_id), stat)
         call check(nf90_inq_varid(ncid, "longitude", lon_id), stat)
         call check(nf90_inq_varid(ncid, "surface_elevation", surface_elevation_id), stat)

         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            !*** Get variable ID's
            call check(nf90_inq_dimid(grpid(n), "nwave", dimid_wave), stat)
            call check(nf90_inquire_dimension(grpid(n), dimid_wave, len=nwave), stat)
            if (nwave .ne. measurement(n)%nwave) then
               call writelog("OUTPUT_L1B_NC_LS: nwave is not equal to number of spectral channels", 8)
            end if

            !*** Spectral data
            call check(nf90_inq_varid(grpid(n), "radiance", rad_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "radiance_noise", radnoise_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "radiance_error", raderror_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "irradiance", irrad_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "irradiance_noise", irradnoise_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "irradiance_error", irraderror_id(n)), stat)
         end do

         !*** Get number of spectra already in file
         call check(nf90_inq_dimid(ncid, "nobs", dimid_nobs), stat)
         call check(nf90_inquire_dimension(ncid, dimid_nobs, len=lastindex), stat)
         lastindex = lastindex + 1

      end if

      !*** WRITE DATA TO FILE

      !*** Indexdata
      call check(nf90_put_var(ncid, x_id, sx, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, y_id, sy, start=(/lastindex/)), stat)

      !*** Geometry
      call check(nf90_put_var(ncid, sza_id, meta%sza, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, vza_id, meta%iza, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, saa_id, meta%saz, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, vaa_id, meta%iaz, start=(/lastindex/)), stat)

      !*** Geodata
      call check(nf90_put_var(ncid, lat_id, meta%lat(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, lon_id, meta%lon(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, surface_elevation_id, meta%surface_elevation, start=(/lastindex/)), stat)

      !*** Spectral data
      start = (/1, lastindex/)
      do n = 1, nwin
         call check(nf90_put_var(grpid(n), rad_id(n), measurement(n)%radiance, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), radnoise_id(n), measurement(n)%radiance_noise, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), raderror_id(n), measurement(n)%radiance_error, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irrad_id(n), measurement(n)%irradiance, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irradnoise_id(n), measurement(n)%irradiance_noise, start=(/1, lastindex/)), stat)
         call check(nf90_put_var(grpid(n), irraderror_id(n), measurement(n)%irradiance_error, start=(/1, lastindex/)), stat)
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (stat .ne. 0) then
         call writelog('OUTPUT_L1B_NC_LS: Error in writing spectra to netCDF file', 8)
      end if

   end subroutine output_l1b_nc_ls
   !------------------------------------------------------------------------------



   subroutine output_lut_nc_js(measurement_hi, atm_scenario, win, meta, xco2, co2_scaling, lut_file)
      !*** input
      type(spectrum), dimension(:), intent(in) :: measurement_hi
      type(atmospheric_scenario), intent(in) :: atm_scenario
      type(window_spectrum), dimension(:), intent(in):: win
      type(metadata), intent(in) :: meta
      real(double), intent(in) :: xco2, co2_scaling
      character(len=*), intent(in) :: lut_file
      !***  local
      integer :: ncid, ierr, natm, nwin, stat, ngroup, dimid_z, dimid_wave, nwave
      integer :: sza_id, z_id, p_id, t_id, h2o_id, co2_id, xco2_id, co2sca_id
      integer, dimension(:), allocatable :: grpid, rad_id, wave_id, alb_id
      integer :: start(2), dimids_spec(2)
      integer :: i, n
      character*1 :: ch
      logical :: exst
      character(stringlen) :: group_name

      stat = 0

      natm = size(atm_scenario%z)
      nwin = size(measurement_hi)

      allocate (grpid(nwin), &
                wave_id(nwin), &
                rad_id(nwin), &
                alb_id(nwin), stat=ierr)

      ! atm_scenario%co2 = atm_scenario%co2/co2_scaling    ! Re-scale meteo co2 column to the original input values

      !*** Create the netCDF file.
      call check(nf90_create(trim(lut_file), nf90_netcdf4, ncid), stat)

      !*** Number of atmospheric layers
      call check(nf90_def_dim(ncid, "z", natm, dimid_z), stat)

      !*** Geometry
      call check(nf90_def_var(ncid, "sza", NF90_double, sza_id), stat)
      call check(nf90_put_att(ncid, sza_id, "unit", "degrees"), stat)
      call check(nf90_put_att(ncid, sza_id, "description", "Solar Zenith Angle"), stat)
      call check(nf90_put_var(ncid, sza_id, meta%sza), stat)

      !*** Meteo
      call check(nf90_def_var(ncid, "height_lev", NF90_double, dimid_z, z_id), stat)
      call check(nf90_put_att(ncid, z_id, "unit", "m"), stat)
      call check(nf90_put_att(ncid, z_id, "description", "Geopotential height at layer boundaries"), stat)
      call check(nf90_put_var(ncid, z_id, atm_scenario%z), stat)

      call check(nf90_def_var(ncid, "press_lev", NF90_double, dimid_z, p_id), stat)
      call check(nf90_put_att(ncid, p_id, "unit", "hPa"), stat)
      call check(nf90_put_att(ncid, p_id, "description", "Pressure at layer boundaries"), stat)
      call check(nf90_put_var(ncid, p_id, atm_scenario%p), stat)

      call check(nf90_def_var(ncid, "temp_lev", NF90_double, dimid_z, t_id), stat)
      call check(nf90_put_att(ncid, t_id, "unit", "K"), stat)
      call check(nf90_put_att(ncid, t_id, "description", "Temperature at layer boundaries"), stat)
      call check(nf90_put_var(ncid, t_id, atm_scenario%t), stat)

      call check(nf90_def_var(ncid, "h2o_lev", NF90_double, dimid_z, h2o_id), stat)
      call check(nf90_put_att(ncid, h2o_id, "unit", "mol/mol"), stat)
      call check(nf90_put_att(ncid, h2o_id, "description", "Volume mixing ratio of H2O at layer boundaries"), stat)
      call check(nf90_put_var(ncid, h2o_id, atm_scenario%h2o), stat)

      call check(nf90_def_var(ncid, "co2_lev", NF90_double, dimid_z, co2_id), stat)
      call check(nf90_put_att(ncid, co2_id, "unit", "mol/mol"), stat)
      call check(nf90_put_att(ncid, co2_id, "description", "Volume mixing ratio of CO2 at layer boundaries"), stat)
      call check(nf90_put_var(ncid, co2_id, atm_scenario%co2/co2_scaling), stat)

      call check(nf90_def_var(ncid, "xco2", NF90_double, xco2_id), stat)
      call check(nf90_put_att(ncid, xco2_id, "unit", "ppm"), stat)
      call check(nf90_put_att(ncid, xco2_id, "description", "Column-averad dry-air mole fraaction of CO2"), stat)
      call check(nf90_put_var(ncid, xco2_id, xco2), stat)

      call check(nf90_def_var(ncid, "co2_scaling", NF90_double, co2sca_id), stat)
      call check(nf90_put_att(ncid, co2sca_id, "unit", "-"), stat)
    call check( nf90_put_att(ncid, co2sca_id, "description", "Scaling factor applied to the co2-profile ('co2_lev') in order to get the specified XCO2 ('xco2')"), stat)
      call check(nf90_put_var(ncid, co2sca_id, co2_scaling), stat)

      do n = 1, nwin
         nwave = measurement_hi(n)%nwave
         !*** Create a group for each spectral window
         write (ch, '(i1.1)') n
         group_name = 'BAND'//ch
         call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

         !*** Surface albedo
         call check(nf90_def_var(grpid(n), "albedo", NF90_double, alb_id(n)), stat)
         call check(nf90_put_att(grpid(n), alb_id(n), "unit", "-"), stat)
         call check(nf90_put_att(grpid(n), alb_id(n), "description", "Surface albedo"), stat)
         call check(nf90_put_var(grpid(n), alb_id(n), win(n)%albedo(1)), stat)

         !*** Spectrally resolved data
         call check(nf90_def_dim(grpid(n), "wavelength", nwave, dimid_wave), stat)

         !*** Wavelength grid
         call check(nf90_def_var(grpid(n), "wavelength", NF90_double, dimid_wave, wave_id(n)), stat)
         call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
         call check(nf90_put_att(grpid(n), wave_id(n), "description", "High-resolution wavelength grid"), stat)
         call check(nf90_put_var(grpid(n), wave_id(n), measurement_hi(n)%wavelength), stat)

         !*** Spectral radiances
         call check(nf90_def_var(grpid(n), "spectral_radiance", NF90_double, dimid_wave, rad_id(n)), stat)
         call check(nf90_put_att(grpid(n), rad_id(n), "unit", "photons s-1 cm-2 sr-1 nm-1"), stat)
         call check(nf90_put_att(grpid(n), rad_id(n), "description", "High-resolution spectral radiances at telescope"), stat)
         call check(nf90_put_var(grpid(n), rad_id(n), measurement_hi(n)%radiance), stat)

      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (stat .ne. 0) then
         call writelog('OUTPUT_LUT_NC_JS: Error in writing spectra to netCDF file', 8)
      end if

   end subroutine output_lut_nc_js

end module spec_interface_create_module

