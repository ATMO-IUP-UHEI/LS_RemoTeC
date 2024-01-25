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
   private :: spectral_response_create_gauss, read_isrf_retrieve, calculate_isrf

contains

   subroutine read_l1b(infile, outputflag, measurement, meta, ierr, synthetic_input_flag, observer_location, win_ini, instr_errors)
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
      real(double) :: min_req_wavelength, max_req_wavelength
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

            ! Check if current band surrounds current fit window. If not, go to the next band
            ! This check needs to take into account the wave boundary offset in multiples of the fwhm
            min_req_wavelength = win_ini(win)%wave_start - win_ini(win)%fwhm * win_ini(win)%wvbd
            max_req_wavelength = win_ini(win)%wave_stop + win_ini(win)%fwhm * win_ini(win)%wvbd
            if (.not. (wavelength(1) <= min_req_wavelength .and. max_req_wavelength <= wavelength(nwave))) then
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
            min_index = maxloc(wavelength, dim=1, mask=wavelength<=min_req_wavelength)
            max_index = minloc(wavelength, dim=1, mask=wavelength>=max_req_wavelength)
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
   subroutine get_isrf_interpolated(flag, filename, fixed_nrow, win_ini, measurement, response_out, ierr)
      integer, intent(in) :: flag
      character(len=*) :: filename
      integer, dimension(:), intent(in) :: fixed_nrow
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      type(instrument_response), dimension(:, :), intent(inout) :: response_out
      integer, intent(out) :: ierr
      !*** local
      integer :: i, j, k, l, n, nils, nwave, nband
      type(instrument_response) :: response
      type(instrument_response), dimension(:), allocatable :: response_tmp
      type(instrument_response), dimension(:, :), allocatable :: response_in
      character*2 :: ch
      character(stringlen) :: message

      nband = size(measurement)
      if (flag == 1) then !Calculate Gaussian ISRF

         call calculate_isrf(win_ini, measurement, response_tmp)
         do n = 1, nband
            do i = 1, fixed_nrow(n)
               response_out(i, n) = response_tmp(n)
            end do
         end do

      elseif (flag == 2) then !Read ISRF from custom NetCDF file
         do n = 1, nband
            write (ch, '(i2.2)') n
            call read_isrf_retrieve(trim(filename)//'isrf_'//ch//'.nc', response, ierr)
            if (ierr .ne. 0) goto 999

            do l = 1, size(win_ini)
               !*** Check if wavelength is in fit window
               if (win_ini(l)%wave_start .ge. minval(measurement(n)%wavelength) .and. &
                   win_ini(l)%wave_stop .le. maxval(measurement(n)%wavelength)) then
                  !*** Interpolate ISRF on user-defined wavelength grid
                  nils = 2*int(int(win_ini(l)%fwhm*win_ini(l)%wvbd/win_ini(l)%reso)/2) + 1
                  do i = 1, fixed_nrow(n)
                     response_out(i, n)%nils = nils
                     response_out(i, n)%nwave = response%nwave
                     allocate (response_out(i, n)%ils_dwave(response%nwave, nils), &
                               response_out(i, n)%resp_store(response%nwave, nils), &
                               stat=ierr)
                     if (ierr .ne. 0) goto 999
                     forall (j=1:nils) &
                        response_out(i, n)%ils_dwave(:, j) = -win_ini(l)%reso*int(nils/2) + dble(j - 1)*win_ini(l)%reso

                     do k = 1, response%nwave
                        call spline_interpol(response%ils_dwave(k, :), response%resp_store(k, :), response%nils, &
                                             response_out(i, n)%ils_dwave(k, :), response_out(i, n)%resp_store(k, :), nils, ierr)
!!$                 call linterp(response%ils_dwave(k,:), response%resp_store(k,:),response%nils,&
!!$                      response_out(i,n)%ils_dwave(k,:), response_out(i,n)%resp_store(k,:), nils, ierr)
                        if (ierr .ne. 0) goto 999
                     end do
                  end do !loop over i
               end if
            end do !loop over l
         end do !loop over n

      else

         write (message, '(a)') 'get_isrf_interpolated: Option not implemented.'
         call writelog(message, 6)

      end if

!!$     do n = 1, nband
!!$        write(ch,'(i2.2)') n
!!$        open(50, file='isrf_'//ch//'.dat')
!!$        do i = 1, response_out(1,n)%nils
!!$           write(50,*) response_out(1,n)%ils_dwave(1,i), response_out(1,n)%resp_store(1,i)
!!$        enddo
!!$        close(50)
!!$     enddo

      ierr = 0
      return
999   call writelog('GET_ISRF_INTERPOLATED: Fatal error during preparation of ISRF', 8)

   end subroutine get_isrf_interpolated

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
         do l = 1, size(win_ini)
            if (win_ini(l)%wave_start .le. maxval(measurement(n)%wavelength) .and. &
                win_ini(l)%wave_stop .ge. minval(measurement(n)%wavelength)) then
               fwhm = win_ini(l)%fwhm
               nils = int(fwhm*win_ini(l)%wvbd/win_ini(l)%reso)
               if (modulo(nils, 2) == 0) nils = nils + 1   ! make sure number of Gaussian ILS points are uneven
               allocate (ilswave(nils), &
                         ilsfunction(nils))
               call spectral_response_create_gauss( &
                  fwhm, win_ini(l)%reso, nils, ilswave, ilsfunction)

               !*** Store in an array(nils) for each measured wavelength as a function of
               !*** wavelength difference
               response(n)%nils = nils
               if (allocated(response(n)%resp_store)) deallocate (response(n)%resp_store)
               if (allocated(response(n)%ils_dwave)) deallocate (response(n)%ils_dwave)
               allocate (response(n)%resp_store(response(n)%nwave, nils))
               allocate (response(n)%ils_dwave(response(n)%nwave, nils))
               do k = 1, response(n)%nwave
                  do i = 1, nils
                     response(n)%resp_store(k, i) = ilsfunction(i)
                     response(n)%ils_dwave(k, i) = ilswave(i)
                  end do
               end do
               deallocate (ilswave, ilsfunction)
            end if
         end do !loop over l
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
   subroutine read_isrf_retrieve(filename, response, ierr)
      character(len=*), intent(in) :: filename
      type(instrument_response), intent(out) :: response
      integer, intent(out) :: ierr
      !*** local
      integer :: ncid, id, i, j
      real(double), dimension(:, :), allocatable :: response_netcdf
      real(double), dimension(:), allocatable :: dw, wavelength

      ! Read the ISRF from the NetCDF file
      call check(nf90_open(trim(filename), nf90_nowrite, ncid), ierr)

      ! Dimension size: number of wavelength differences
      call check(nf90_inq_dimid(ncid, "dwl", id), ierr)
      call check(nf90_inquire_dimension(ncid, id, len=response%nils), ierr)

      ! Dimension size: number of measured wavelengths
      call check(nf90_inq_dimid(ncid, "wl_i", id), ierr)
      call check(nf90_inquire_dimension(ncid, id, len=response%nwave), ierr)

      ! Allocate data fields of the ISRF
      allocate (dw(response%nils))
      allocate (response%ils_dwave(response%nwave, response%nils), stat=ierr) ! Wavelength differences for which the ISRF is defined ! {{{
      allocate (wavelength(response%nwave), stat=ierr) ! Measured wavelengths at which the ISRF is defined ! {{{
      allocate (response%wavelength(response%nwave), stat=ierr)
      allocate (response_netcdf(response%nils, response%nwave), stat=ierr) ! Response function representative at different measured wavelengths as function of the wavelength difference. !

      ! Fill fields with data from NetCDF file

      ! Dimension domain: Wavelength differences
      call check(nf90_inq_varid(ncid, "Wavelength_differences", id), ierr)
      call check(nf90_get_var(ncid, id, dw), ierr)

      ! Dimension domain: Measured wavelengths
      call check(nf90_inq_varid(ncid, "Measured_wavelengths", id), ierr)
      call check(nf90_get_var(ncid, id, wavelength), ierr)

      ! Data field: Response
      call check(nf90_inq_varid(ncid, "Response", id), ierr)
      call check(nf90_get_var(ncid, id, response_netcdf), ierr)

      allocate (response%resp_store(response%nwave, response%nils))
      do i = 1, response%nwave
         response%wavelength(i) = wavelength(i)
         do j = 1, response%nils
            response%ils_dwave(i, j) = dw(j)
            response%resp_store(i, j) = response_netcdf(j, i)
         end do
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)

      ! Mark success
      ierr = 0
      return

   end subroutine read_isrf_retrieve

!------------------------------------------------------------------------------
end module spectrum_interface_module
