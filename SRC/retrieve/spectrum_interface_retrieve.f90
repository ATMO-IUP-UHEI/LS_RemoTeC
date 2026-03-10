module spectrum_interface_module
   use header_module
   use read_settings_module, only: window_ini
   use read_errors_module, only: instrument_errors
   use spectrum_internal_module, only: spectrum, instrument_interface
   use spectral_response_module, only: instrument_response, response_internal
   use auxiliary_routines_module, only: check
   use netcdf
   implicit none
   private
!*** types
   public :: spectrum, instrument_response

!*** Procedures
   public :: read_l1b, get_isrf_interpolated, instrument_interface
   private :: get_ils_wavelengths, get_ils_offsets, get_ils_response, get_ils_response_internal, ils_gauss, get_ils_response_from_file, read_ils_file

contains

   subroutine read_l1b(infile, outputflag, measurement, meta, ierr, synthetic_input_flag, line_number, observer_location, win_ini, instr_errors)
      !** Input
      character(len=*), intent(in) :: infile
      integer, intent(in) :: outputflag
      integer, intent(in) :: synthetic_input_flag
      integer, intent(in) :: observer_location
      type(window_ini), dimension(:), intent(in), optional :: win_ini
      type(instrument_errors), dimension(:), intent(in), optional :: instr_errors
      !*** Output
      type(spectrum), dimension(:), allocatable, intent(out) :: measurement
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
      integer, intent(out) :: line_number
      !*** local variables
      integer :: nstokes_l1b, nst
      real(double), dimension(4) :: s = (/1.d0, 0.D0, 0.D0, 0.d0/)
      real(double) :: lambda_shifted, continuum
      integer :: k, l, i, n, nwave, imid
      integer :: win, nwin, band, nband
      integer, dimension(:), allocatable :: pixelid
      integer :: ncid, grpid(3), varid, dimid_lat, dimid_lon, dimid_time, dimid_wave
      integer :: sx, sy, start1d(1), start2d(2), start3d(3)
      integer :: time_id, sza_id, vza_id, saa_id, vaa_id, observer_altitude_id, lon_id, lat_id! , surface_elevation_id
      real(double), dimension(:), allocatable :: wavelength
      integer :: min_index, max_index
      character(stringlen) :: spectrum_file, index_info, message
      real(double), dimension(:), allocatable :: var

      !-------------------------------------------------------------------

      !*** Set nstokes of measurement same to nstokes of model here
      nstokes_l1b = nstokes

      i = INDEX(infile, '.nc')
      index_info = infile(i + 4:)
      spectrum_file = trim(infile(:i + 2))

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      line_number = sx

      if (outputflag >= 2) then
         call writelog('*** Start of READ_L1B ***', 1)
         call writelog('L1B file: '//trim(spectrum_file), 1)
      end if

      !*** Read spectrum
      call check(nf90_open(trim(spectrum_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('READ_L1B: error opening file: '//trim(spectrum_file), 6)
         return
      end if

      ! location of data in netcdf file
      start1d = (/sy/)
      start2d = (/sx, sy/)
      start3d = (/1, sx, sy/)

      !*** Timedata
      call netcdf_get_var(ncid, "time", meta%seconds_since_reference, start1d)

      !*** Geodata
      call netcdf_get_var(ncid, "latitude", meta%lat(1), start2d)
      call netcdf_get_var(ncid, "longitude", meta%lon(1), start2d)
      ! For now, the center coordinates are used as corner coordinates as well
      meta%lat(:) = meta%lat(1)
      meta%lon(:) = meta%lon(1)

      !*** Geometry
      call netcdf_get_var(ncid, "solar_zenith_angle", meta%sza, start2d)
      call netcdf_get_var(ncid, "viewing_zenith_angle", meta%iza, start2d)
      call netcdf_get_var(ncid, "solar_azimuth_angle", meta%saz, start2d)
      call netcdf_get_var(ncid, "viewing_azimuth_angle", meta%iaz, start2d)
      if (observer_location .eq. 1) then
         call netcdf_get_var(ncid, "observer_altitude", meta%observer_height, start2d)
      end if
      ! calculate relative azimuth angle
      meta%phi = dabs(meta%iaz - meta%saz)

      ! Set up measurement
      ! Get nwin from settings file
      nwin = size(win_ini)

      allocate (measurement(nwin), stat=ierr)
      if (ierr .ne. 0) return

      measurement(:)%sza = meta%sza
      measurement(:)%iza = meta%iza
      measurement(:)%phi = meta%phi
      
      if (observer_location .eq. 1) then
         measurement(:)%observer_height = meta%observer_height
      end if

      call check(nf90_inq_grps(ncid, nband, grpid), ierr)
      if (ierr .ne. 0) return

      do win = 1, nwin
         ! For each fit window get the wavelength and radiance from the correct band in the data
         do band = 1, nband
            ! Check if this band contains a wavelength grid that surrounds the fit window

            ! Get number of spectral points
            call check(nf90_inq_dimid(grpid(band), "channel", dimid_wave), ierr)
            if (ierr .ne. 0) return
            call check(nf90_inquire_dimension(grpid(band), dimid_wave, len=nwave), ierr)
            if (ierr .ne. 0) return

            ! Get wavelengths of the current band and write them into dummy variable wavelength
            if (allocated(wavelength)) deallocate(wavelength)
            allocate(wavelength(nwave))
            call netcdf_get_vector_var(grpid(band), "wavelength", wavelength, start=(/1/))

            ! Check if current band surrounds current fit window. If not, go to the next band.
            if (wavelength(1) > win_ini(win)%wave_start .or. wavelength(nwave) < win_ini(win)%wave_stop) then
               if (band == nband) then
                  print*, "ERROR IN READ_L1B: No bands surround fit window."
               end if
               cycle
            end if

            ! Correct band found. Get the spectral information.
            measurement(win)%nwave = nwave

            ! Deallocate / allocate
            if (allocated(measurement(win)%wavelength)) then
               if (size(measurement(win)%wavelength) /= measurement(win)%nwave) then
                  deallocate (measurement(win)%wavelength)
                  deallocate (measurement(win)%radiance)
                  deallocate (measurement(win)%radiance_noise)
                  if (synthetic_input_flag == 1) then
                     deallocate (measurement(win)%radiance_error)
                  end if
                  deallocate (measurement(win)%mask)
                  deallocate (measurement(win)%measurement_stokesc)
               end if
            end if
            if (.not. allocated(measurement(win)%wavelength)) then
               allocate(measurement(win)%wavelength(measurement(win)%nwave))
               allocate(measurement(win)%radiance(measurement(win)%nwave))
               allocate(measurement(win)%radiance_noise(measurement(win)%nwave))
               if (synthetic_input_flag == 1) then
                  allocate(measurement(win)%radiance_error(measurement(win)%nwave))
               end if
               allocate(measurement(win)%mask(measurement(win)%nwave))
            end if
            if (.not. allocated(measurement(win)%measurement_stokesc) .and. nstokes_l1b > 1) then
               allocate (measurement(win)%measurement_stokesc(nstokes_l1b))
            end if

            !*** Get spectrum
            measurement(win)%wavelength = wavelength
            call netcdf_get_vector_var(grpid(band), "radiance", measurement(win)%radiance, start3d)
            call netcdf_get_vector_var(grpid(band), "radiance_noise", measurement(win)%radiance_noise, start3d)
            if (synthetic_input_flag == 1) then
               call netcdf_get_vector_var(grpid(band), "radiance_error", measurement(win)%radiance_error, start3d)
            end if

            measurement(win)%mask = 0
            if (nstokes_l1b > 1) then
               do nst = 1, nstokes_l1b
                  measurement(win)%measurement_stokesc(nst) = s(nst)
               end do
            end if

            ! Cut down spectrum to necessary wavelength range
            min_index = minloc(wavelength, dim=1, mask=wavelength>=win_ini(win)%wave_start)
            max_index = maxloc(wavelength, dim=1, mask=wavelength<=win_ini(win)%wave_stop)
            nwave = max_index - min_index + 1

            measurement(win)%nwave = nwave

            if (allocated(var)) then
               deallocate(var)
            end if
            allocate(var(nwave))

            var = measurement(win)%wavelength(min_index:max_index)
            deallocate(measurement(win)%wavelength)
            allocate(measurement(win)%wavelength(nwave))
            measurement(win)%wavelength = var

            var = measurement(win)%radiance(min_index:max_index)
            deallocate(measurement(win)%radiance)
            allocate(measurement(win)%radiance(nwave))
            measurement(win)%radiance = var

            var = measurement(win)%radiance_noise(min_index:max_index)
            deallocate(measurement(win)%radiance_noise)
            allocate(measurement(win)%radiance_noise(nwave))
            measurement(win)%radiance_noise = var

            if (synthetic_input_flag == 1) then
               var = measurement(win)%radiance_error(min_index:max_index)
               deallocate(measurement(win)%radiance_error)
               allocate(measurement(win)%radiance_error(nwave))
               measurement(win)%radiance_error = var
            end if

            exit
         end do
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_L1B: Error opening/reading spectrum_file '//trim(spectrum_file)
            call writelog(message, 6)
         end if
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_L1B***', 1)
      end if

   end subroutine read_l1b



   subroutine netcdf_get_var(ncid, varname, values, start)
      ! input
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      integer, dimension(:), intent(in) :: start
      ! output
      real(double), intent(out) :: values
      ! local variables
      integer :: varid, ierr_inq, ierr_get

      call check(nf90_inq_varid(ncid, varname, varid), ierr_inq)
      call check(nf90_get_var(ncid, varid, values, start=start), ierr_get)

      if (ierr_inq /= 0 .or. ierr_get /= 0) then
         print*, "spectrum_interface_retrieve.f90: ERROR reading variable ", varname
         return
      end if
   end subroutine netcdf_get_var



   subroutine netcdf_get_vector_var(ncid, varname, values, start)
      ! input
      integer, intent(in) :: ncid
      character(len=*), intent(in) :: varname
      integer, dimension(:), intent(in) :: start
      ! output
      real(double), dimension(:), intent(out) :: values
      ! local variables
      integer :: varid, ierr_inq, ierr_get

      call check(nf90_inq_varid(ncid, varname, varid), ierr_inq)
      call check(nf90_get_var(ncid, varid, values, start=start), ierr_get)

      if (ierr_inq /= 0 .or. ierr_get /= 0) then
         print*, "spectrum_interface_retrieve.f90: ERROR reading variable ", varname
         return
      end if
   end subroutine netcdf_get_vector_var



!------------------------------------------------------------------------------
!> @details This routine gets the ISRF on the appropriate spectral grids
!! For each measured spectral pixel (lo-reso grid with nwave_lo points),
!! the ISRF is given on the model spectral grid (hi-reso grid with nils points)
!! Thus, we have to interpolate twice when reading the S5P-format file.
!! For testing purposes, one can also choose to calculate a Gaussian ISRF inline or read a custom ISRF file
!------------------------------------------------------------------------------
   subroutine get_isrf_interpolated(flag_ilscalc, flag_output, filename, fixed_nrow, win_ini, measurement, response, ierr)
      integer, intent(in) :: flag_ilscalc, flag_output
      character(len=*) :: filename
      integer, dimension(:), intent(in) :: fixed_nrow
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      !type(instrument_response), dimension(:, :), intent(inout) :: response_out
      type(instrument_response), dimension(:, :), allocatable, intent(out) :: response
      integer, intent(out) :: ierr
      !*** local
      integer :: i, j, k, l, n, nils, nwave, nwin, nband
      integer :: row, win, ils, wave
      ! type(instrument_response) :: response
      ! type(instrument_response), dimension(:), allocatable :: response_tmp
      ! type(instrument_response), dimension(:, :), allocatable :: response_in
      integer :: nrow
      real(double), dimension(:), allocatable :: ils_dwave
      integer :: io
      character*2 :: ch
      character(stringlen) :: message

      ! Initialize instrument_response
      ! assumption: nrow is the same for all windows.
      ! if this is ever changed, change this allocation and all do loops in the below subroutines
      nwin = size(win_ini)
      nrow = fixed_nrow(1)
      allocate(response(nrow, nwin))

      ! wavelengths on which ils are defined
      call get_ils_wavelengths(response, measurement, nwin, nrow)

      ! wavelength offsets for which the ils is defined
      call get_ils_offsets(response, win_ini, nwin, nrow)

      ! response of the ils
      call get_ils_response(response, win_ini, flag_ilscalc, filename, nwin, nrow)

      ! write used ils for each window into a debug output file
      if (flag_output >= 3) then
         do win = 1, size(win_ini)
            write (ch, "(i2.2)") win
            open(newunit(io), file="CONTRL_OUT/used_ils_"//ch//".dat")
            write(io, *) "left wavelength / nm"
            write(io, *) response(win, 1)%wavelength(1)
            write(io, *) "right wavelength / nm"
            write(io, *) response(win, 1)%wavelength(response(1, 1)%nwave)
            write(io, *) "nwave"
            write(io, *) response(win, 1)%nwave
            write(io, *) "ils dwave / nm, left ils, right ils"
            do ils = 1, response(win, 1)%nils
               write(io, *) response(win, 1)%ils_dwave(1, ils), response(win, 1)%resp_store(1, ils), response(win, 1)%resp_store(response(win, 1)%nwave, ils)
            enddo
            close(io)
         enddo
      end if

      print*, "TODO LS: Implement ierr in the subroutines"
      if (ierr .ne. 0) goto 999

      ierr = 0
      return

999   call writelog('GET_ISRF_INTERPOLATED: Fatal error during preparation of ISRF', 8)

   end subroutine get_isrf_interpolated


   subroutine get_ils_wavelengths(response, measurement, nwin, nrow)
      integer, intent(in) :: nwin
      integer, intent(in) :: nrow
      type(spectrum), dimension(:), intent(in) :: measurement
      type(instrument_response), dimension(:, :), allocatable, intent(inout) :: response
      !*** local variables
      integer :: win, row

      do win = 1, nwin
         do row = 1, nrow
            response(row, win)%nwave = measurement(win)%nwave
            response(row, win)%wavelength = measurement(win)%wavelength
         end do ! loop over row
      end do ! loop over win
   end subroutine get_ils_wavelengths


   subroutine get_ils_offsets(response, win_ini, nwin, nrow)
      integer, intent(in) :: nwin
      integer, intent(in) :: nrow
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(instrument_response), dimension(:, :), allocatable, intent(inout) :: response
      !*** local variables
      integer :: win, row, nils, ils, nwave, wave
      real(double), dimension(:), allocatable :: ils_dwave

      do win = 1, nwin
         do row = 1, nrow
            nils = win_ini(win)%fwhm * win_ini(win)%wvbd / win_ini(win)%reso
            ! make sure nils is odd to have a single center point
            if (modulo(nils, 2) == 0) then
               nils = nils + 1
            end if
            response(row, win)%nils = nils
            !
            ! get offsets ils_dwave for each wavelength wave on which ils is defined
            ! ils_dwave is the same for all wavelengths, therefore no loop over wave
            if(allocated(ils_dwave)) deallocate(ils_dwave)
            allocate(ils_dwave(nils))
            do ils = 1, nils
               ! convert index offset into wavelength offset
               ils_dwave(ils) = (ils - (int(nils/2) + 1)) * win_ini(win)%reso
            end do ! loop over ils
            nwave = response(row, win)%nwave
            allocate(response(row, win)%ils_dwave(nwave, nils))
            do wave = 1, nwave
               response(row, win)%ils_dwave(wave, :) = ils_dwave
            end do ! loop over wave
         end do ! loop over row
      end do ! loop over win
   end subroutine get_ils_offsets


   subroutine get_ils_response(response, win_ini, flag_ilscalc, filename, nwin, nrow)
      integer, intent(in) :: flag_ilscalc, nwin, nrow
      type(instrument_response), dimension(:, :), allocatable, intent(inout) :: response
      type(window_ini), dimension(:), intent(in) :: win_ini
      character(len=*) :: filename

      if (flag_ilscalc == 1) then ! Calculate Gaussian ILS
         call get_ils_response_internal(response, win_ini, nwin, nrow)
      elseif (flag_ilscalc == 2) then ! Read ILS from file
         call get_ils_response_from_file(response, filename, nwin, nrow)
      else
         print*, "invalid flag ilscalc"
      end if
   end subroutine get_ils_response


   subroutine get_ils_response_internal(response, win_ini, nwin, nrow)
      type(instrument_response), dimension(:, :), allocatable, intent(inout) :: response
      type(window_ini), dimension(:), intent(in) :: win_ini
      integer, intent(in) :: nwin, nrow
      !*** local variables
      integer :: win, row, wave
      real(double), dimension(:), allocatable :: resp

      do win = 1, nwin
         do row = 1, nrow
            ! assume that ils_dwave grid is the same for all wavelengths, therefore
            ! no loop over wave. Use wave index 1 for all wavelengths.
            if(allocated(resp)) deallocate(resp)
            allocate(resp(response(win, row)%nils))
            call ils_gauss(response(win, row)%ils_dwave(1, :), win_ini(win)%reso, win_ini(win)%fwhm, resp)
            allocate(response(win, row)%resp_store(response(win, row)%nwave, response(win, row)%nils))
            do wave = 1, response(win, row)%nwave
               response(win, row)%resp_store(wave, :) = resp
            end do ! loop over wave
         end do ! loop over row
      end do ! loop over win
   end subroutine get_ils_response_internal


   subroutine ils_gauss(ils_dwave, reso, fwhm, resp)
      real(double), dimension(:), intent(in) :: ils_dwave
      real(double), intent(in) :: fwhm, reso
      real(double), dimension(:), intent(out) :: resp
      !*** local variables
      integer :: ils
      real(double) :: sigma, sqrt2pi, norm

      sigma = fwhm * 0.42466090
      sqrt2pi = 2.50662827
      norm = 0

      do ils = 1, size(ils_dwave)
         resp(ils) = 1 / (sigma * sqrt2pi) * dexp(-0.5 * (ils_dwave(ils)/sigma)**2)
         norm = norm + resp(ils)*reso
      end do ! loop over ils

      !*** normalize
      resp = resp / norm
   end subroutine ils_gauss


   subroutine get_ils_response_from_file(response, filename, nwin, nrow)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nwin, nrow
      type(instrument_response), dimension(:, :), allocatable, intent(inout) :: response
      !*** local variables
      type(instrument_response) :: response_from_file, response_tmp
      integer :: nwave, wave, nils, ils
      integer :: nwave_from_file, nils_from_file
      integer :: nwave_tmp, nils_tmp
      integer :: win, row
      integer :: ierr

      do win = 1, nwin
         do row = 1, nrow
            call read_ils_file(filename, response(win, row), response_from_file)

            ! dimensions from file
            nwave_from_file = response_from_file%nwave
            nils_from_file = response_from_file%nils

            ! target dimensions
            nwave = response(win, row)%nwave
            nils = response(win, row)%nils

            ! dimensions after first interpolation step
            nwave_tmp = nwave
            nils_tmp = nils_from_file

            ! interpolate onto lores wavelength grid
            if(allocated(response_tmp%wavelength)) deallocate(response_tmp%wavelength)
            if(allocated(response_tmp%ils_dwave)) deallocate(response_tmp%ils_dwave)
            if(allocated(response_tmp%resp_store)) deallocate(response_tmp%resp_store)
            allocate(response_tmp%wavelength(nwave_tmp))
            allocate(response_tmp%ils_dwave(nwave_tmp, nils_tmp))
            allocate(response_tmp%resp_store(nwave_tmp, nils_tmp))

            response_tmp%wavelength = response(win, row)%wavelength
            response_tmp%nwave = nwave_tmp
            response_tmp%nils = nils_tmp

            do ils = 1, nils_from_file
               call spline_interpol( &
                  response_from_file%wavelength, response_from_file%ils_dwave(:, ils), nwave_from_file, &
                  response_tmp%wavelength, response_tmp%ils_dwave(:, ils), nwave_tmp, &
                  ierr &
               )
               call spline_interpol( &
                  response_from_file%wavelength, response_from_file%resp_store(:, ils), nwave_from_file, &
                  response_tmp%wavelength, response_tmp%resp_store(:, ils), nwave_tmp, &
                  ierr & 
               )
            end do ! loop over ils

            ! interpolate onto hires ils_dwave grid
            allocate(response(win, row)%resp_store(nwave, nils))

            do wave = 1, nwave_tmp
               call spline_interpol( &
                  response_tmp%ils_dwave(wave, :), response_tmp%resp_store(wave, :), nils_tmp, &
                  response(win, row)%ils_dwave(wave, :), response(win, row)%resp_store(wave, :), nils, &
                  ierr &
               )
            end do ! loop over wave
         end do ! loop over row
      end do ! loop over win
   end subroutine get_ils_response_from_file


   subroutine read_ils_file(filename, response, response_from_file)
      character(len=*), intent(in) :: filename
      type(instrument_response), intent(in) :: response
      type(instrument_response), intent(out) :: response_from_file
      !*** local variables
      real(double), dimension(:), allocatable :: wavelength, ils_dwave
      real(double), dimension(:, :), allocatable :: resp
      integer :: nwave, wave, nils, band, nband, current_band
      integer :: ncid, grpid(3), varid
      integer :: ierr
      integer :: ierr_band_found ! number of data bands in file that surrounded fit window (debug)

      call check(nf90_open(trim(filename), nf90_nowrite, ncid), ierr)

      ! get correct band from file. It's wavelength range has to surround
      ! response%wavelength. The netcdf group containing this band needs to be used
      call check(nf90_inq_grps(ncid, nband, grpid), ierr)
      if (ierr .ne. 0) return

      ierr_band_found = 0

      do band = 1, nband
         ! check if this band contains a wavelength grid that surrounds the fit window
         call check(nf90_inq_dimid(grpid(band), "channel", varid), ierr)
         call check(nf90_inquire_dimension(grpid(band), varid, len=nwave), ierr)

         ! Get wavelengths of the current band and write them into dummy variable wavelength
         if (allocated(wavelength)) deallocate(wavelength)
         allocate(wavelength(nwave))
         call check(nf90_inq_varid(grpid(band), "wavelength_center", varid), ierr)
         call check(nf90_get_var(grpid(band), varid, wavelength), ierr)

         ! check if current band surrounds current fit window. If not, go to the next band
         if (wavelength(1) > response%wavelength(1) .or. wavelength(nwave) < response%wavelength(response%nwave)) then
            if (band == nband) then
               if (ierr_band_found == 0) then
                  print*, "ERROR IN READ_ILS_FILE: No bands surrounded fit window."
               else
                  print*, "ERROR IN READ_ILS_FILE: ", ierr_band_found, " band(s) surrounded fit window but none had sufficiently large ils_dwave grid."
               end if
            end if
            cycle
         else
            ierr_band_found = ierr_band_found + 1
         end if

         ! check if this band contains a wavelength offset grid that is sufficiently large
         call check(nf90_inq_dimid(grpid(band), "d_channel", varid), ierr)
         call check(nf90_inquire_dimension(grpid(band), varid, len=nils), ierr)

         ! Get wavelength offsets of the current band and write them into dummy variable ils_dwave
         if (allocated(ils_dwave)) deallocate(ils_dwave)
         allocate(ils_dwave(nils))
         call check(nf90_inq_varid(grpid(band), "wavelength_offset", varid), ierr)
         call check(nf90_get_var(grpid(band), varid, ils_dwave), ierr)

         ! check if current band has sufficiently large ils_dwave. If not, go to the next band
         if (ils_dwave(1) > maxval(response%ils_dwave(:, 1)) .or. ils_dwave(nils) < minval(response%ils_dwave(:, response%nils))) then
            cycle
         end if

         current_band = band
         exit
      end do

      ! correct band found to be current_band
      ! we have nwave, nils, wavelength, and ils_dwave, write those into response_from_file
      ! also get the correct response

      ! wavelengths on which ils is defined
      response_from_file%nwave = nwave
      allocate(response_from_file%wavelength(nwave))
      response_from_file%wavelength = wavelength

      ! wavelength offsets for which ils is defined
      allocate(response_from_file%ils_dwave(nwave, nils))
      response_from_file%nils = nils
      do wave = 1, nwave
         response_from_file%ils_dwave(wave, :) = ils_dwave
      end do ! loop over wave

      ! response of ils
      allocate(resp(nils, nwave))
      call check(nf90_inq_varid(grpid(current_band), "response", varid), ierr)
      call check(nf90_get_var(grpid(current_band), varid, resp), ierr)

      allocate(response_from_file%resp_store(nwave, nils))
      do wave = 1, nwave
         response_from_file%resp_store(wave, :) = resp(:, wave)
      end do ! loop over wave

      call check(nf90_close(ncid), ierr)

      ! Mark success
      ierr = 0
      return
   end subroutine read_ils_file
end module spectrum_interface_module
