! modified by Andre Galli, February 25, 2011
!        CHANGES: Introduced spectral_response_lbl as default since it is 0.1 sec faster per retrieval
!        than the FFT-based subroutine spectral_response_custom
! modified by Andre Galli, October 4, 2011
!        CHANGES: Switched to updated merged solar line list by Toon [2010, personal comm.], which also allows to
!        observe only a fraction of the solar disk
!Otto: Should we use the approach of van Deelen et al, Appl. Opt., 46, 243-252, 2007?

module solar_model_create_module
   use header_module
   use read_settings_module, only: window_ini
   use auxiliary_routines_module, only: check
   implicit none
   private

   public :: sun_spectrum, read_sun_netcdf, read_sun_tsis1_hsrs, interpolate_solar_spectrum !, solar_model, read_sun_ascii

contains
   !------------------------------------------------------------------------------
   !*** Read in solar line list and compute modelled solar spectrum
   subroutine solar_model(solar_dir, outputflag, win_ini)
      implicit none
      !*** Input
      character(199), intent(in) :: solar_dir
      integer, intent(in) :: outputflag
      !*** Input/output
      type(window_ini), dimension(:), allocatable, intent(inout) :: win_ini
      !*** local variables
      integer :: i, k, n, io, nwin
      integer :: iostat, nlines, intdum
      integer, dimension(1) :: kstart, kstop
      real(double) :: x, d4, y2, x2, ss, sunfrac, ff, sdc, sdi, wdc, wdi, ddc, ddi
      real(double), dimension(:), allocatable :: position, strength, wwidth, dwidth
      real(double), dimension(:), allocatable :: sun_spectrum_hi
      !-----------------------------------------------------------------

      if (outputflag >= 2) then
         write (*, '(a)')
         write (*, '(a)') '*****************************************'
         write (*, '(a)') '***** Reading solar spectrum'
         write (*, '(a)') '*****'
      end if

      ! Solar line input from Foster, Aaron and Toon [JPL, 2010, personal communication]
      !     i=0
      !     OPEN(14,FILE=TRIM(solar_dir)//'solar_di_20071105.101', action = 'read', status = 'old')!
      !     DO
      !        READ(14,*,IOSTAT=iostat)
      !        IF(iostat<0)EXIT
      !        i=i+1
      !     ENDDO
      !     CLOSE(14)
      !     nlines=i
      !
      !     ALLOCATE(position(nlines),strength(nlines),wwidth(nlines),dwidth(nlines))
      !
      !     OPEN(14,FILE=TRIM(solar_dir)//'solar_di_20071105.101',FORM='FORMATTED',action = 'read', status = 'old')
      !     DO k=1,nlines
      !        READ(14,'(i3,f12.6,2e10.3,f5.4,f5.4,f10.4,f8.4,1x,a36)')&
      !               intdum,position(k),strength(k),wwidth(k),dwidth(k),&
      !               doubledum,doubledum,doubledum
      !     ENDDO
      !     CLOSE(14)

      !*** Use the merged solar line list by Toon [2010, private comm] that is identical to the one used by TCCON:
      !*** Note: For backscattered radiation, the sun fraction always is 1.0
      sunfrac = 1.0
      if (sunfrac .gt. 1.0) then
         write (*, *) 'Warning: solar_spectrum: frac > 1', sunfrac
         sunfrac = 1.
      end if
      ff = sunfrac**2

      !*** First determine number of lines
      i = 0
      open (newunit(io), FILE=trim(solar_dir)//'solar_merged_20110401.108', FORM='FORMATTED', &
            action='read', status='old')
      do
         read (io, *, IOSTAT=iostat)
         if (iostat < 0) exit
         i = i + 1
      end do
      close (io)
      nlines = i
      allocate (position(nlines), strength(nlines), wwidth(nlines), dwidth(nlines))
      !*** Now actually read in the file
      open (newunit(io), FILE=trim(solar_dir)//'solar_merged_20110401.108', FORM='FORMATTED', &
            action='read', status='old')
      do k = 1, nlines
         read (io, *) intdum, position(k), sdc, sdi, wdc, wdi, ddc, ddi
         strength(k) = (1 - ff)*sdc + ff*sdi
         wwidth(k) = (1 - ff)*wdc + ff*wdi
         dwidth(k) = (1 - ff)*ddc + ff*ddi
      end do
      close (io)

      !*** Calculate solar spectrum
      !*** Start loop over all windows
      nwin = size(win_ini)
      do n = 1, nwin
         if (allocated(sun_spectrum_hi)) deallocate (sun_spectrum_hi)
         allocate (sun_spectrum_hi(win_ini(n)%nwave_hi))
         kstart = minloc(DABS(position - win_ini(n)%wavelength_hi(1)))
         kstart = max(kstart, 1)
         kstop = minloc(DABS(position - win_ini(n)%wavelength_hi(win_ini(n)%nwave_hi)))
         sun_spectrum_hi = 0.D0
         do k = kstart(1), kstop(1)
            d4 = (dwidth(k)**2 + (5.D-6*position(k)*dsqrt(sunfrac))**2)**2 ! Doppler width + solar rotation term
            y2 = wwidth(k)**2
            ss = strength(k)
            do i = 1, win_ini(n)%nwave_hi
               x = win_ini(n)%wavelength_hi(i) - position(k)
               x2 = x*x
               sun_spectrum_hi(i) = sun_spectrum_hi(i) + &
                                    ss*DEXP(-x2/DSQRT(d4 + x2*y2*(1.D0 + DABS(x/(wwidth(k) + 0.07)))))
            end do
         end do
         sun_spectrum_hi = DEXP(-sun_spectrum_hi)
         allocate (win_ini(n)%sun_spectrum_ref_hi(win_ini(n)%nwave_hi))
         win_ini(n)%sun_spectrum_ref_hi = sun_spectrum_hi

         !**** End loop over all windows
      end do

      if (outputflag >= 2) then
         write (*, '(a)') '***** Done.'
         write (*, '(a)') '*****************************************'
      end if

   end subroutine solar_model
   !------------------------------------------------------------------------------
   !> @brief Read reference irradiance from KNMI's NetCDF file
   !> @param[in]  sun_file      file for reference irradiance
   !> @param[out] sun_input     datatype for reference irradiance
   !> @param[out] ierr          error identifier: 0=normal, 1=error opening file, 2=error in reading, 3=allocation error
   subroutine read_sun_netcdf(sun_file, sun_input, ierr)
      use netcdf
      !*** Input
      character(len=*), intent(in) :: sun_file
      !*** output
      type(sun_spectrum), intent(out) :: sun_input
      integer, intent(out) :: ierr
      !*** local
      integer :: i, ncid, dimID, varID, nwave
      real(double), dimension(:), allocatable ::  wavelength, irradiance

      !*** Open NetCDF LUT
      call check(nf90_open(trim(sun_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_DIMID(ncid, "wavelength_hr", dimID), ierr)
      call check(nf90_inquire_dimension(ncid, dimid, len=nwave), ierr)
      allocate ( &
         wavelength(nwave), &
         irradiance(nwave), &
         sun_input%wavelength(nwave), &
         sun_input%irradiance(nwave), &
         stat=ierr)

      if (ierr .ne. 0) then
         ierr = ierr_all
         call stopretrieval('READ_SUN_NETCDF: memory allocation error')
         return
      end if
      sun_input%nwave = nwave

      call check(NF90_INQ_VARID(ncid, "wavelength_hr", varID), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varID, wavelength), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "irradiance_flux", varID), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varID, irradiance), ierr)
      if (ierr .ne. 0) return

      call check(nf90_close(ncid), ierr)

      do i = 1, nwave
         sun_input%wavelength(i) = wavelength(i)
         sun_input%irradiance(i) = irradiance(i)
      end do

      ierr = 0
      return

   end subroutine read_sun_netcdf

   !------------------------------------------------------------------------------
   !> @brief Read reference irradiance from TSIS-1 HSRS NetCDF file (https://doi.org/10.1029/2020GL091709)
   !> @param[in]  sun_file       file for reference irradiance
   !> @param[out] sun_input      datatype for reference irradiance
   !> @param[out]  ierr          error identifier: 0=normal, 1=error opening file, 2=error in reading, 3=allocation error
   subroutine read_sun_tsis1_hsrs(sun_file, sun_input, ierr)
      use netcdf
      !*** Input
      character(len=*), intent(in) :: sun_file
      !*** Output
      type(sun_spectrum), intent(out) :: sun_input
      integer, intent(out) :: ierr
      !*** Local
      integer :: i, ncid, dimid, varid, nwave
      real(double), dimension(:), allocatable :: wavelength, irradiance

      !*** Open NetCDF LUT
      call check(nf90_open(trim(sun_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inq_dimid(ncid, "wavelength", dimid), ierr)
      call check(nf90_inquire_dimension(ncid, dimid, len=nwave), ierr)
      allocate(&
         wavelength(nwave),&
         irradiance(nwave),&
         sun_input%wavelength(nwave),&
         sun_input%irradiance(nwave),&
         stat=ierr)
      if (ierr .ne. 0) then
         ierr = ierr_all
         call stopretrieval("READ_SUN_TSIS1_HSRS: memory allocation error")
      end if

      sun_input%nwave = nwave

      call check(nf90_inq_varid(ncid, "Vacuum Wavelength", varid), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, varid, wavelength), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inq_varid(ncid, "SSI", varid), ierr)
      if (ierr .ne. 0) return
      call check(nf90_get_var(ncid, varid, irradiance), ierr)
      if (ierr .ne. 0) return

      call check(nf90_close(ncid) ,ierr)

      ! Convert irradiance from W m-2 nm-1 to photons s-1 cm-2 nm-1
      irradiance = irradiance * wavelength * 1d-9 / h_planck / c_light * 1d-4

      do i = 1, nwave
         sun_input%wavelength(i) = wavelength(i)
         sun_input%irradiance(i) = irradiance(i)
      end do

      ierr = 0
      return
   end subroutine

   !------------------------------------------------------------------------------
   !> @brief Read reference irradiance from KNMI's ascii file
   !> @param[in]  sun_file      file for reference irradiance
   !> @param[out] sun_input     datatype for reference irradiance
   !> @param[out] ierr          error identifier: 0=normal, 1=error opening file, 2=error in reading, 3=allocation error
   subroutine read_sun_ascii(sun_file, sun_input, ierr)
      !*** Input
      character(len=*), intent(in) :: sun_file
      !*** output
      type(sun_spectrum), intent(out) :: sun_input
      integer, intent(out) :: ierr
      !*** local
      integer :: i, io, nlines, hlines
      real(double) :: wavelength, sun_mw !in units [mW/m^2/nm]

      !*** Read high resolution sun spectrum
      !*** First, get number of lines
      open (newunit(io), FILE=trim(sun_file), FORM='FORMATTED', action='read', STATUS='OLD', iostat=ierr, err=100)

      read (io, *, err=101) hlines
      do i = 1, hlines - 1
         read (io, *, err=101)
      end do
      i = 0
      do
         read (io, *, IOSTAT=ierr, err=101)
         if (ierr < 0) exit
         i = i + 1
      end do
      close (io)
      nlines = i
      allocate (sun_input%wavelength(nlines), sun_input%irradiance(nlines), stat=ierr)
      if (ierr .ne. 0) then
         ierr = ierr_all
         call stopretrieval('READ_SUN_ASCII: memory allocation error')
         return
      end if
      sun_input%nwave = nlines

      !*** Now, read solar spectrum
      open (newunit(io), FILE=trim(sun_file), FORM='FORMATTED', action='read', STATUS='OLD', iostat=ierr, err=100)
      do i = 1, hlines
         read (io, *, err=101)
      end do
      do i = 1, nlines
         read (io, *, iostat=ierr, err=101) wavelength, sun_mw, sun_input%irradiance(i)
         sun_input%wavelength(i) = wavelength
      end do
      close (io)

      ierr = 0
      return

100   continue
      ierr = ierr_open
      call stopretrieval('READ_SUN_ASCII: error opening file')
      return

101   continue
      ierr = ierr_read
      call stopretrieval('READ_SUN_SPECTRUM: error reading file')
      close (io)
      return

   end subroutine read_sun_ascii

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

end module solar_model_create_module
