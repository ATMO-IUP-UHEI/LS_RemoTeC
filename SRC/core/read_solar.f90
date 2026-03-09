module read_solar_module
   use header_module
   use read_settings_module, only: window_ini
   use auxiliary_routines_module, only: check
   implicit none
   private

   public :: read_solar, interpolate_solar_spectrum

contains
   !------------------------------------------------------------------------------
   !> @brief Read reference irradiance from standardized solar file
   !> @param[in]  sun_file       file for reference irradiance
   !> @param[out] sun_input      datatype for reference irradiance
   !> @param[out]  ierr          error identifier: 0=normal, 1=error opening file, 2=error in reading, 3=allocation error
   subroutine read_solar(sun_file, sun_input, ierr)
      use netcdf
      !*** input
      character(len=*), intent(in) :: sun_file
      !*** output
      type(sun_spectrum), intent(out) :: sun_input
      integer, intent(out) :: ierr
      !*** local
      integer :: i, ncid, dimid, varid, nwave
      real(double), dimension(:), allocatable :: wavelength, irradiance

      print*, "Start of READ_SOLAR"

      !*** Open NetCDF
      call check(nf90_open(trim(sun_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) return

      !*** Get dimension info, allocate variables
      call check(nf90_inq_dimid(ncid, "channel", dimid), ierr)
      call check(nf90_inquire_dimension(ncid, dimid, len=nwave), ierr)
      allocate(&
         wavelength(nwave),&
         irradiance(nwave),&
         sun_input%wavelength(nwave),&
         sun_input%irradiance(nwave),&
         stat=ierr)
      if (ierr .ne. 0) then
         ierr = ierr_all
         call stopretrieval("READ_SOLAR: memory allocation error")
      end if

      sun_input%nwave = nwave

      !*** Get data
      call check(nf90_inq_varid(ncid, "wavelength", varid), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, varid, wavelength), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inq_varid(ncid, "solar_irradiance", varid), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, varid, irradiance), ierr)
      if (ierr .ne. 0) return

      !*** Close NetCDF
      call check(nf90_close(ncid) ,ierr)

      !** Copy tmp variables into sun_input
      do i = 1, nwave
         sun_input%wavelength(i) = wavelength(i)
         sun_input%irradiance(i) = irradiance(i)
      end do

      print*, "End of READ_SOLAR"

      ierr = 0
      return
   end subroutine

   !------------------------------------------------------------------------------

   subroutine interpolate_solar_spectrum(sun_input, win_ini, ierr)
      type(sun_spectrum), intent(in) :: sun_input
      type(window_ini), dimension(:), intent(inout) :: win_ini
      integer, intent(out) :: ierr
      !*** local
      integer :: nwin, n
      character(199) :: message

      !*** Initialize error identifier
      ierr = 0

      !*** Interpolate hi-res sun to hi-res wavelength grid
      !*** Start loop over all windows
      nwin = size(win_ini)
      do n = 1, nwin
         allocate (win_ini(n)%sun_spectrum_ref_hi(win_ini(n)%nwave_hi), stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'INTERPOLATE_SOLAR_SPECTRUM: memory allocation error'
            ierr = ierr_all
            goto 999
         end if

         call spline_interpol(sun_input%wavelength, sun_input%irradiance, sun_input%nwave, &
                              win_ini(n)%wavelength_hi, win_ini(n)%sun_spectrum_ref_hi(:), win_ini(n)%nwave_hi, ierr)
         if (ierr .ne. 0) then
            write (message, *) 'INTERPOLATE_SOLAR_SPECTRUML.SPLINE_INTERPOL.SPLINT: bad input'
            ierr = ierr_intrpl
            goto 999
         end if
      end do

      return
999   continue
      call stopretrieval(message)

   end subroutine interpolate_solar_spectrum

   !------------------------------------------------------------------------------

end module read_solar_module
