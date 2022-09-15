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
   public :: read_spectrum, read_l1b_nc_js, read_l1b, get_isrf_interpolated, instrument_interface
   private :: spectral_response_create_gauss, read_isrf_sim, calculate_isrf

contains

!------------------------------------------------------------------------------
!> @details Read in synthetic spectrum
!------------------------------------------------------------------------------
   subroutine read_spectrum(infile, ipixel, outputflag, measurement, meta, ierr)
      !** Input
      character(len=*), intent(in) :: infile
      integer, intent(in) :: ipixel, outputflag
      !*** Output
      type(spectrum), dimension(:), allocatable, intent(out) :: measurement
      type(metadata), intent(in) :: meta
      integer, intent(out) :: ierr
      !*** local variables
      integer :: i, n, nwin, nwave, npixel, nobs
      integer, dimension(:), allocatable :: pixelid
      integer :: ncid, varid, dimid, start(1), count(1), start_prof(2), count_prof(2)
      integer, dimension(2) :: grpid
      real(double), dimension(:), allocatable :: var
      character(stringlen) :: message
      !-------------------------------------------------------------------

      if (outputflag >= 2) then
         call writelog('*** Start of READ_SPECTRUM ***', 1)
         write (message, *) ipixel
         call writelog('L1B file: '//trim(infile)//' PixelID: '//trim(message), 1)
      end if
      !*** Read spectrum
      call check(nf90_open(trim(infile)//'_sim.nc', nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('READ_SPECTRUM: error opening file: '//trim(infile)//'_sim.nc', 6)
         return
      end if

      !*** Get number of observations
      call check(NF90_INQ_DIMID(ncid, "nobs", dimid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_INQUIRE_DIMENSION(ncid, dimid, len=nobs), ierr)
      allocate (pixelid(nobs))
      call check(NF90_INQ_VARID(ncid, "pixelID", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, pixelid), ierr)
      if (ierr .ne. 0) return

      do i = 1, nobs
         if (pixelid(i) == ipixel) then
            npixel = i
            exit
         end if
      end do
      if (npixel > nobs) then
         call writelog("READ_SPECTRUM: pixelid not in spectrumfile", 6)
         ierr = ierr_l1b
         return
      end if
      start = [npixel]
      count = [1]
      start_prof = [1, npixel]

      !*** Get number of bands
      call check(nf90_inq_grps(ncid, nwin, grpid), ierr)
      allocate (measurement(nwin), stat=ierr)
      do n = 1, nwin

         measurement(n)%sza = meta%sza
         measurement(n)%iza = meta%iza
         measurement(n)%phi = dabs(meta%iaz - meta%saz)

         !*** Get number of spectral points
         call check(NF90_INQ_DIMID(grpid(n), "nwave", dimid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_INQUIRE_DIMENSION(grpid(n), dimid, len=nwave), ierr)
         if (ierr .ne. 0) return
         measurement(n)%nwave = nwave
         count_prof = [nwave, 1]

         !*** Deallocate/allocate
         if (.not. allocated(measurement(n)%wavelength)) then
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         elseif (size(measurement(n)%wavelength) .ne. measurement(n)%nwave) then
            deallocate (measurement(n)%wavelength)
            deallocate (measurement(n)%radiance)
            deallocate (measurement(n)%radiance_noise)
            deallocate (measurement(n)%radiance_error)
            deallocate (measurement(n)%irradiance)
            deallocate (measurement(n)%irradiance_noise)
            deallocate (measurement(n)%irradiance_error)
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         end if

         if (allocated(var)) deallocate (var)
         allocate (var(nwave))
         !*** Get wavenumbers
         call check(NF90_INQ_VARID(grpid(n), "wavelength", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%wavelength = var

         !*** Get radiance
         call check(NF90_INQ_VARID(grpid(n), "radiance", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance = var

         !*** Get radiance_noise
         call check(NF90_INQ_VARID(grpid(n), "radiance_noise", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance_noise = var

         !*** Get radiance_error
         call check(NF90_INQ_VARID(grpid(n), "radiance_error", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance_error = var

         !*** Get irradiance
         call check(NF90_INQ_VARID(grpid(n), "irradiance", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance = var

         !*** Get irradiance_noise
         call check(NF90_INQ_VARID(grpid(n), "irradiance_noise", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance_noise = var

         !*** Get irradiance_error
         call check(NF90_INQ_VARID(grpid(n), "irradiance_error", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start_prof, count=count_prof), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance_error = var

      end do

      if (outputflag >= 2) then
         call writelog('*** End of READ_SPECTRUM ***', 1)
      end if

   end subroutine read_spectrum

!------------------------------------------------------------------------------
!>
!------------------------------------------------------------------------------
   subroutine read_l1b_nc_js(infile, outputflag, measurement, meta, ierr, win_ini, instr_errors)
      !** Input
      character(len=*), intent(in) :: infile
      integer, intent(in) :: outputflag
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
      integer :: k, l, i, n, nwin, nwave, imid
      integer, dimension(:), allocatable :: pixelid
      integer :: ncid, grpid(3), varid, dimid_lat, dimid_lon, dimid_time, dimid_wave
      integer :: sx, sy, start2d(2), start3d(3)
      integer :: time_id, sza_id, vza_id, saz_id, vaz_id, lon_id, lat_id, elev_id
      real(double), dimension(:), allocatable :: var
      character(stringlen) :: spectrum_file, index_info, message

      !-------------------------------------------------------------------

      !*** Set nstokes of measurement same to nstokes of model here
      nstokes_l1b = nstokes

      i = INDEX(infile, '.nc')
      index_info = infile(i + 4:)
      spectrum_file = trim(infile(:i + 2))

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 3), '(I3)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 3), '(I3)') sy

      if (outputflag >= 2) then
         call writelog('*** Start of READ_L1B_NC_JS ***', 1)
         call writelog('L1B file: '//trim(spectrum_file), 1)
      end if

      !*** Read spectrum
      call check(nf90_open(trim(spectrum_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('READ_L1B_NC_JS: error opening file: '//trim(spectrum_file), 6)
         return
      end if

      start2d = (/sx, sy/)

      !*** Geometry
      call check(nf90_inq_varid(ncid, "sza", sza_id), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inq_varid(ncid, "vza", vza_id), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inq_varid(ncid, "saa", saz_id), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inq_varid(ncid, "vaa", vaz_id), ierr)
      if (ierr .ne. 0) return

      call check(nf90_get_var(ncid, sza_id, meta%sza, start=start2d), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, vza_id, meta%iza, start=start2d), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, saz_id, meta%saz, start=start2d), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, vaz_id, meta%iaz, start=start2d), ierr)
      if (ierr .ne. 0) return
      meta%phi = dabs(meta%iaz - meta%saz)

      !*** Geodata
      call check(nf90_inq_varid(ncid, "latitude", lat_id), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inq_varid(ncid, "longitude", lon_id), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inq_varid(ncid, "elevation", elev_id), ierr)
      if (ierr .ne. 0) return

      call check(nf90_get_var(ncid, lat_id, meta%lat(1), start=start2d), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, lon_id, meta%lon(1), start=start2d), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, elev_id, meta%surface_elevation, start=start2d), ierr)
      if (ierr .ne. 0) return

      !*** Timedata
      call check(NF90_INQ_VARID(ncid, "time", time_id), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, time_id, meta%time(1:6)), ierr)
      if (ierr .ne. 0) return

      !*** For now the center coordinates are used as corner coordinates as well
      meta%lon(:) = meta%lon(1)
      meta%lat(:) = meta%lat(1)

      call check(nf90_inq_grps(ncid, nwin, grpid), ierr)
      if (ierr .ne. 0) return

      allocate (measurement(nwin), stat=ierr)
      if (ierr .ne. 0) return

      measurement(:)%sza = meta%sza
      measurement(:)%iza = meta%iza
      measurement(:)%phi = meta%phi

      start3d = (/1, sx, sy/)

      do n = 1, nwin

         !*** Get number of spectral points
         call check(NF90_INQ_DIMID(grpid(n), "nwave", dimid_wave), ierr)
         if (ierr .ne. 0) return
         call check(NF90_INQUIRE_DIMENSION(grpid(n), dimid_wave, len=nwave), ierr)
         if (ierr .ne. 0) return
         measurement(n)%nwave = nwave

         !*** Deallocate/allocate
         if (.not. allocated(measurement(n)%measurement_stokesc) .and. nstokes_l1b > 1) then
            allocate (measurement(n)%measurement_stokesc(nstokes_l1b))
         end if
         if (.not. allocated(measurement(n)%wavelength)) then
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%mask(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         elseif (size(measurement(n)%wavelength) .ne. measurement(n)%nwave) then
            deallocate (measurement(n)%wavelength)
            deallocate (measurement(n)%radiance)
            deallocate (measurement(n)%mask)
            deallocate (measurement(n)%radiance_noise)
            deallocate (measurement(n)%radiance_error)
            deallocate (measurement(n)%irradiance)
            deallocate (measurement(n)%irradiance_noise)
            deallocate (measurement(n)%irradiance_error)
            deallocate (measurement(n)%measurement_stokesc)
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%mask(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         end if

         if (allocated(var)) deallocate (var)
         allocate (var(nwave))
         !*** Get wavelengths
         call check(NF90_INQ_VARID(grpid(n), "wavelength", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%wavelength = var

         !*** Get radiance
         call check(NF90_INQ_VARID(grpid(n), "radiance", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance = var

         !*** Get radiance_noise
         call check(NF90_INQ_VARID(grpid(n), "radiance_noise", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance_noise = var

         !*** Get radiance_error
         call check(NF90_INQ_VARID(grpid(n), "radiance_error", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%radiance_error = var

         !*** Get irradiance
         call check(NF90_INQ_VARID(grpid(n), "irradiance", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance = var

         !*** Get irradiance_noise
         call check(NF90_INQ_VARID(grpid(n), "irradiance_noise", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance_noise = var

         !*** Get irradiance_error
         call check(NF90_INQ_VARID(grpid(n), "irradiance_error", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start3d), ierr)
         if (ierr .ne. 0) return
         measurement(n)%irradiance_error = var

         measurement(n)%mask = 0
         if (nstokes_l1b > 1) then
            do nst = 1, nstokes_l1b
               measurement(n)%measurement_stokesc(nst) = s(nst)
            end do
         end if

         if (present(instr_errors)) then
            do l = 1, size(win_ini)

               if (win_ini(l)%wave_start .le. maxval(measurement(n)%wavelength) .and. &
                   win_ini(l)%wave_stop .ge. maxval(measurement(n)%wavelength)) then
                  !*** Add instrument errors specified in INI/errors.in
                  continuum = maxval(measurement(n)%radiance(:))
                  imid = measurement(n)%nwave/2
                  !*** instrument response function !HH: this should be moved to win_ini%fwhm
                  !            measurement(n)%fwhm =  measurement(n)%fwhm*(1.d0+0.01d0*instr_errors(n)%isrf)
                  do k = 1, measurement(n)%nwave
                     !*** radiometric offset calibration error (additive constant)
                    measurement(n)%radiance(k) = max(0.d0, measurement(n)%radiance(k) + 0.01d0*instr_errors(n)%rad_offset*continuum)

                     !*** shift/squeeze in Earth spectrum wavelength grid
                     measurement(n)%wavelength(k) = measurement(n)%wavelength(k) + instr_errors(n)%earth_shift &
                                   + instr_errors(n)%earth_squeeze1*(measurement(n)%wavelength(k) - measurement(n)%wavelength(imid))

                     !*** radiometric gain calibration error (scaling factor)
                     measurement(n)%radiance(k) = measurement(n)%radiance(k)*(1.d0 + 0.01d0*instr_errors(n)%rad_gain)
                  end do
               end if
            end do
         end if

      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_L1B_NC_JS: Error opening/reading spectrum_file '//trim(spectrum_file)
            call writelog(message, 6)
         end if
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_L1B_NC_JS ***', 1)
      end if

   end subroutine read_l1b_nc_js

!------------------------------------------------------------------------------
!  Read in synthetic spectrum in L1B format/units
!------------------------------------------------------------------------------
   subroutine read_l1b(spectrum_file, outputflag, measurement, meta, ierr, win_ini, instr_errors)
      !** Input
      character(len=*), intent(in) :: spectrum_file
      type(window_ini), dimension(:), intent(in), optional :: win_ini
      type(instrument_errors), dimension(:), intent(in), optional :: instr_errors
      integer, intent(in) :: outputflag
      !*** Output
      type(spectrum), dimension(:), allocatable, intent(out) :: measurement
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
      !*** local variables
      integer :: nstokes_l1b
      real(double), dimension(4) :: s = (/1.d0, 0.D0, 0.D0, 0.d0/)
      integer :: k, l, n, nwin_in, io, imid, nst
      real(double) :: lambda_shifted, continuum, sza, iza, saz, iaz, fwhm, time(6)
      character(stringlen) :: message
      logical :: iopen
      !-------------------------------------------------------------------

      if (outputflag >= 2) then
         call writelog('*** Start of READ_L1B ***', 1)
         call writelog('L1B file: '//trim(spectrum_file), 1)
      end if

      !*** Set nstokes of measurement same to nstokes of model here
      nstokes_l1b = nstokes
      !$OMP critical
      open (newunit=io, FILE=trim(spectrum_file), FORM='FORMATTED', status='old', action='read', iostat=ierr)
      !$OMP end critical
      if (ierr .ne. 0) goto 100
      !*** Get instrument info
      read (io, *, iostat=ierr, err=100) meta%time
      do k = 1, 5
         read (io, *, iostat=ierr, err=100) meta%lat(k), meta%lon(k) !1: latitude and longitude of pixel corners
         !2 = right,upper corner
      end do                          !3 = right,lower corner etc.
      read (io, *, iostat=ierr, err=100) meta%surface_elevation
      read (io, *, iostat=ierr, err=100) sza, iza, saz, iaz
      read (io, *, iostat=ierr, err=100) nwin_in
      allocate (measurement(nwin_in), stat=ierr)
      measurement(:)%sza = sza
      measurement(:)%iza = iza
      measurement(:)%phi = dabs(iaz - saz)
      !*** Get spectral measurements
      do n = 1, nwin_in
         read (io, *, iostat=ierr, err=100) measurement(n)%nwave
         read (io, *, iostat=ierr, err=100) fwhm

         if (.not. allocated(measurement(n)%measurement_stokesc) .and. nstokes_l1b > 1) then
            allocate (measurement(n)%measurement_stokesc(nstokes_l1b))
         end if
         if (.not. allocated(measurement(n)%wavelength)) then
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%mask(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         elseif (size(measurement(n)%wavelength) .ne. measurement(n)%nwave) then
            deallocate (measurement(n)%wavelength)
            deallocate (measurement(n)%radiance)
            deallocate (measurement(n)%mask)
            deallocate (measurement(n)%radiance_noise)
            deallocate (measurement(n)%radiance_error)
            deallocate (measurement(n)%irradiance)
            deallocate (measurement(n)%irradiance_noise)
            deallocate (measurement(n)%irradiance_error)
            deallocate (measurement(n)%measurement_stokesc)
            allocate ( &
               measurement(n)%wavelength(measurement(n)%nwave), &
               measurement(n)%radiance(measurement(n)%nwave), &
               measurement(n)%mask(measurement(n)%nwave), &
               measurement(n)%radiance_noise(measurement(n)%nwave), &
               measurement(n)%radiance_error(measurement(n)%nwave), &
               measurement(n)%irradiance(measurement(n)%nwave), &
               measurement(n)%irradiance_noise(measurement(n)%nwave), &
               measurement(n)%irradiance_error(measurement(n)%nwave))
         end if
         do k = 1, measurement(n)%nwave
            read (io, *, iostat=ierr, err=100) measurement(n)%wavelength(k), &
               measurement(n)%radiance(k), &
               measurement(n)%radiance_noise(k), &
               measurement(n)%radiance_error(k), &
               measurement(n)%irradiance(k), &
               measurement(n)%irradiance_noise(k), &
               measurement(n)%irradiance_error(k)

         end do

         measurement(n)%mask = 0
         if (nstokes_l1b > 1) then
            do nst = 1, nstokes_l1b
               measurement(n)%measurement_stokesc(nst) = s(nst)
            end do
         end if
         if (present(instr_errors)) then
            do l = 1, size(win_ini)

               if (win_ini(l)%wave_start .le. maxval(measurement(n)%wavelength) .and. &
                   win_ini(l)%wave_stop .ge. maxval(measurement(n)%wavelength)) then
                  !*** Add instrument errors specified in INI/errors.in
                  continuum = maxval(measurement(n)%radiance(:))
                  imid = measurement(n)%nwave/2
                  !*** instrument response function !HH: this should be moved to win_ini%fwhm
                  !            measurement(n)%fwhm =  measurement(n)%fwhm*(1.d0+0.01d0*instr_errors(n)%isrf)
                  do k = 1, measurement(n)%nwave
                     !*** radiometric offset calibration error (additive constant)
                    measurement(n)%radiance(k) = max(0.d0, measurement(n)%radiance(k) + 0.01d0*instr_errors(n)%rad_offset*continuum)

                     !*** shift/squeeze in Earth spectrum wavelength grid
                     measurement(n)%wavelength(k) = measurement(n)%wavelength(k) + instr_errors(n)%earth_shift &
                                   + instr_errors(n)%earth_squeeze1*(measurement(n)%wavelength(k) - measurement(n)%wavelength(imid))

                     !*** radiometric gain calibration error (scaling factor)
                     measurement(n)%radiance(k) = measurement(n)%radiance(k)*(1.d0 + 0.01d0*instr_errors(n)%rad_gain)
                  end do
               end if
            end do
         end if

      end do
      close (io)

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_L1B: Error opening/reading spectrum_file '//trim(spectrum_file)
            call writelog(message, 6)
         end if
         inquire (unit=io, opened=iopen)
         if (iopen) close (io)
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_L1B ***', 1)
      end if

   end subroutine read_l1b

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
            call read_isrf_sim(trim(filename)//'isrf_'//ch//'.nc', response, ierr)
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
   subroutine read_isrf_sim(filename, response, ierr)
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

   end subroutine read_isrf_sim

!------------------------------------------------------------------------------
end module spectrum_interface_module
