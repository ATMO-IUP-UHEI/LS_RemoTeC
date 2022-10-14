!------------------------------------------------------------------------------
!> Read input atmosphere, interpolate onto model grids
!------------------------------------------------------------------------------
module atmosphere_interface_module
   use header_module
   use read_settings_module, only: window_ini
   use auxiliary_routines_module
   use atmosphere_internal_module, only: atmospheric_scenario
   use read_miprep_module
   use netcdf
   implicit none
   private

!*** public types/procedures
   public :: atmospheric_scenario
   public :: read_atmosphere, read_atmosphere_nc_js, read_ecmwf, read_aux

contains
!-----------------------------------------------------------------------------
!> @details Read in synthetic atmospheric scenario on pressure+height grid
!-----------------------------------------------------------------------------
   subroutine read_atmosphere(meteo_file, outputflag, atm_scenario, meta, ierr, meteo_errors)
!*** Input
      character(*), intent(in) :: meteo_file
      real(double), dimension(4), intent(in), optional :: meteo_errors
      integer, intent(in) :: outputflag
!*** Output
      type(atmospheric_scenario), intent(out) :: atm_scenario
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
!*** Local variables
      integer :: ninput, k, io
      real(double) :: sza, iza, saz, iaz, phi, lat(5), lon(5)
      integer :: year, month, day, hour, min
      real(double) :: sec
      character(stringlen) :: message
      logical :: iopen
!-----------------------------------------------------

      if (outputflag >= 2) then
         call writelog('*** Start of READ_ATMOSPHERE ***', 1)
      end if
      !$OMP critical
      open (newunit=io, FILE=meteo_file, FORM='FORMATTED', status='old', action='read', iostat=ierr)
      !$OMP end critical
      if (ierr .ne. 0) goto 100
      read (io, *, iostat=ierr, err=100) ninput

      if (allocated(atm_scenario%z)) deallocate (atm_scenario%z)
      if (allocated(atm_scenario%p)) deallocate (atm_scenario%p)
      if (allocated(atm_scenario%t)) deallocate (atm_scenario%t)
      if (allocated(atm_scenario%h2o)) deallocate (atm_scenario%h2o)
      if (allocated(atm_scenario%co2)) deallocate (atm_scenario%co2)
      if (allocated(atm_scenario%ch4)) deallocate (atm_scenario%ch4)
      if (allocated(atm_scenario%co)) deallocate (atm_scenario%co)
      allocate ( &
         atm_scenario%z(ninput + 1), &
         atm_scenario%p(ninput + 1), &
         atm_scenario%t(ninput + 1), &
         atm_scenario%h2o(ninput + 1), &
         atm_scenario%co2(ninput + 1), &
         atm_scenario%ch4(ninput + 1), &
         atm_scenario%co(ninput + 1))
      do k = 1, ninput + 1
         read (io, *, iostat=ierr, err=100) &
            atm_scenario%z(k), &
            atm_scenario%p(k), &
            atm_scenario%t(k), &
            atm_scenario%h2o(k), &
            atm_scenario%co2(k), &
            atm_scenario%ch4(k), &
            atm_scenario%co(k)
      end do
      read (io, *, iostat=ierr, err=100) year, &
         month, &
         day, &
         hour, &
         min, &
         sec
      do k = 1, 5
         read (io, *, iostat=ierr, err=100) lat(k), lon(k) !1: latitude and longitude of pixel corners
         !2 = right,upper corner
      end do                          !3 = right,lower corner etc.
      read (io, *, iostat=ierr, err=100) atm_scenario%surface_elevation
      read (io, *, iostat=ierr, err=100) sza, iza, saz, iaz
      close (io)

      phi = dabs(iaz - saz)
      atm_scenario%surface_wspeed = 0.d0
      meta%time(1) = year
      meta%time(2) = month
      meta%time(3) = day
      meta%time(4) = hour
      meta%time(5) = min
      meta%time(6) = int(sec)
      meta%lat = lat
      meta%lon = lon
      meta%surface_elevation = atm_scenario%surface_elevation

      if (present(meteo_errors)) then
         !*** Introduce error in pressure profile
         atm_scenario%p = (1.d0 + 0.01*meteo_errors(3))*atm_scenario%p
         atm_scenario%t(:) = atm_scenario%t(:) + meteo_errors(4)
      end if

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_ATMOSPHERE: Error opening/reading meteo_file '//trim(meteo_file)
            call writelog(message, 6)
         end if
         inquire (unit=io, opened=iopen)
         if (iopen) close (io)
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_ATMOSPHERE ***', 1)
      end if

   end subroutine read_atmosphere

   subroutine read_atmosphere_nc_js(infile, outputflag, atm_scenario, meta, ierr, meteo_errors)
      !*** Input
      character(*), intent(in) :: infile
      real(double), dimension(4), intent(in), optional :: meteo_errors
      integer, intent(in) :: outputflag
      !*** Output
      type(atmospheric_scenario), intent(out) :: atm_scenario
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
      !*** Local variables
      integer :: i, k, n, nz
      integer :: ncid, varid, dimid_lat, dimid_lon, dimid_time, dimid_z
      integer :: sx, sy, start2d(2), start3d(3)
      integer :: time_id, lon_id, lat_id, elev_id, z_id
      real(double), dimension(:), allocatable :: var
      character(stringlen) :: meteo_file, index_info, message
!-----------------------------------------------------

      i = INDEX(infile, '.nc')
      index_info = infile(i + 4:)
      meteo_file = trim(infile(:i + 2))

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 3), '(I3)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 3), '(I3)') sy

      if (outputflag >= 2) then
         call writelog('*** Start of READ_ATMOSPHERE_NC_JS ***', 1)
         call writelog('ATM file: '//trim(meteo_file), 1)
      end if

      !*** Read meteodata
      call check(nf90_open(trim(meteo_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('READ_ATMOSPHERE_NC_JS: error opening file: '//trim(meteo_file), 6)
         return
      end if

      start2d = (/sx, sy/)

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

      !*** For now the center corrdinates are used as corner coordinates as well
      meta%lon(:) = meta%lon(1)
      meta%lat(:) = meta%lat(1)

      atm_scenario%surface_elevation = meta%surface_elevation
      atm_scenario%surface_wspeed = 0.d0

      start3d = (/1, sx, sy/)

      !*** Get number of vertical levels
      call check(NF90_INQ_DIMID(ncid, "nz", dimid_z), ierr)
      if (ierr .ne. 0) return
      call check(NF90_INQUIRE_DIMENSION(ncid, dimid_z, len=nz), ierr)
      if (ierr .ne. 0) return

      !*** Deallocate/allocate
      if (allocated(atm_scenario%z)) deallocate (atm_scenario%z)
      if (allocated(atm_scenario%p)) deallocate (atm_scenario%p)
      if (allocated(atm_scenario%t)) deallocate (atm_scenario%t)
      if (allocated(atm_scenario%h2o)) deallocate (atm_scenario%h2o)
      if (allocated(atm_scenario%co2)) deallocate (atm_scenario%co2)
      if (allocated(atm_scenario%ch4)) deallocate (atm_scenario%ch4)
      if (allocated(atm_scenario%co)) deallocate (atm_scenario%co)
      allocate ( &
         atm_scenario%z(nz), &
         atm_scenario%p(nz), &
         atm_scenario%t(nz), &
         atm_scenario%h2o(nz), &
         atm_scenario%co2(nz), &
         atm_scenario%ch4(nz), &
         atm_scenario%co(nz))

      if (allocated(var)) deallocate (var)
      allocate (var(nz))

      !*** Get height
      call check(NF90_INQ_VARID(ncid, "z", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%z = var

      !*** Get pressure
      call check(NF90_INQ_VARID(ncid, "p", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%p = var

      !*** Get temperature
      call check(NF90_INQ_VARID(ncid, "t", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%t = var

      !*** Get water vapor
      call check(NF90_INQ_VARID(ncid, "h2o", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%h2o = var

      !*** Get co2
      call check(NF90_INQ_VARID(ncid, "co2", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%co2 = var

      !*** Get ch4
      call check(NF90_INQ_VARID(ncid, "ch4", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%ch4 = var

      !*** Get co
      call check(NF90_INQ_VARID(ncid, "co", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, var, start=start3d), ierr)
      if (ierr .ne. 0) return
      atm_scenario%co = var

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)

      if (present(meteo_errors)) then
         !*** Introduce error in pressure profile
         atm_scenario%p = (1.d0 + 0.01*meteo_errors(3))*atm_scenario%p
         atm_scenario%t(:) = atm_scenario%t(:) + meteo_errors(4)
      end if

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_ATMOSPHERE_NC_JS: Error opening/reading meteo_file '//trim(meteo_file)
            call writelog(message, 6)
         end if
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_ATMOSPHERE_NC_JS ***', 1)
      end if

   end subroutine read_atmosphere_nc_js

!------------------------------------------------------------------------------
!> @details Read in synthetic atmospheric scenario on ECMWF pressure grid
!! Convert hybrid coefficients to pressure
!! Convert specific humidity to VMR of water
!! Get height profile from ideal gas law
!------------------------------------------------------------------------------
   subroutine read_ecmwf(meteo_file, outputflag, atm_scenario, ierr)
      !*** Input
      character(*), intent(in) :: meteo_file
      integer, intent(in) :: outputflag
      !*** Output
      type(atmospheric_scenario), intent(out) :: atm_scenario
      integer, intent(out) :: ierr
      !*** Local variables
      character(stringlen) :: message
      integer :: i, nlev, io
      real(double) :: q, hyai, hybi, hyam, hybm, psurf, surface_elevation
      logical :: iopen

      if (outputflag >= 2) then
         call writelog('*** Start of READ_ECMWF ***', 1)
      end if
      !$OMP critical
      open (newunit=io, FILE=meteo_file, FORM='FORMATTED', status='old', action='read', iostat=ierr)
      !$OMP end critical
      if (ierr .ne. 0) goto 100
      read (io, *, iostat=ierr, err=100) nlev
      if (.not. allocated(atm_scenario%p)) then
         allocate (atm_scenario%p(nlev + 1))
         allocate (atm_scenario%z(nlev + 1))
         allocate (atm_scenario%t(nlev + 1))
         allocate (atm_scenario%h2o(nlev + 1))
         allocate (atm_scenario%ch4(nlev + 1))
         allocate (atm_scenario%co(nlev + 1))
      elseif (size(atm_scenario%p) .ne. nlev + 1) then
         deallocate (atm_scenario%p)
         deallocate (atm_scenario%z)
         deallocate (atm_scenario%t)
         deallocate (atm_scenario%h2o)
         deallocate (atm_scenario%ch4)
         deallocate (atm_scenario%co)
         allocate (atm_scenario%p(nlev + 1))
         allocate (atm_scenario%z(nlev + 1))
         allocate (atm_scenario%t(nlev + 1))
         allocate (atm_scenario%h2o(nlev + 1))
         allocate (atm_scenario%ch4(nlev + 1))
         allocate (atm_scenario%co(nlev + 1))
      end if
      read (io, *, iostat=ierr, err=100) psurf
      read (io, *, iostat=ierr, err=100) surface_elevation
      surface_elevation = surface_elevation/grav   !convert [m**2 s**-2] to [m]
      do i = 1, nlev
         read (io, *, iostat=ierr, err=100) hyai, hybi, hyam, hybm, atm_scenario%t(i), q, atm_scenario%ch4(i), atm_scenario%co(i)
         !*** Convert hybrid coefficients to pressure [from Pa to hPa]
         atm_scenario%p(i) = (hyam + hybm*psurf)/100.
         !*** Convert specific humidity to dry VMR of water
         !*** Specific humidity s  is mass of water per mass of HUMID air.
         !*** The mass mixing ratio s* is mass of water per mass of DRY air: s* = s / (1 - s).
         !*** The volume mixing ratio is given by: (s*) * Mair / Mwater = s / (1 - s) * Mair / Mwater.
         atm_scenario%h2o(i) = q/(1.0d0 - q)*1.60855
      end do
      close (io)

      !*** Get height from barometric formula, i.e. ideal gas law
      !*** from surface to top of atmopshere
      atm_scenario%z(nlev) = surface_elevation - &
                             rg*atm_scenario%t(nlev)/(air_m*grav)*log(atm_scenario%p(nlev)/(psurf/100.))
      do i = 1, nlev - 1
         atm_scenario%z(nlev - i) = atm_scenario%z(nlev - i + 1) - &
                                 rg*atm_scenario%t(nlev - i)/(air_m*grav)*log(atm_scenario%p(nlev - i)/atm_scenario%p(nlev - i + 1))
      end do
      atm_scenario%surface_elevation = surface_elevation
      atm_scenario%surface_pressure = psurf/100.   !convert [Pa] to [hPa]

      !*** Get surface values
      atm_scenario%p(nlev + 1) = psurf/100
      atm_scenario%z(nlev + 1) = atm_scenario%surface_elevation
      atm_scenario%h2o(nlev + 1) = atm_scenario%h2o(nlev)
      atm_scenario%t(nlev + 1) = atm_scenario%t(nlev) + 0.0065*(atm_scenario%z(nlev) - atm_scenario%z(nlev + 1))
      atm_scenario%surface_wspeed = 0.d0

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_ECMWF: Error opening/reading meteo_file '//trim(meteo_file)
            call writelog(message, 6)
         end if
         inquire (unit=io, opened=iopen)
         if (iopen) close (io)
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_ECMWF ***', 1)
      end if

   end subroutine read_ecmwf

!------------------------------------------------------------------------------
   subroutine read_aux(infile, ipixel, outputflag, atm_scenario, meta, ierr)
      !*** Input
      character(*), intent(in) :: infile
      integer, intent(in) :: ipixel, outputflag
      !*** Output
      type(atmospheric_scenario), intent(out) :: atm_scenario
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
      !*** Local variables
      character(stringlen) :: message
      integer :: ninput
      real(double), dimension(:), allocatable :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable :: colev_in     ! VMR CO

      if (outputflag >= 2) then
         call writelog('*** Start of READ_AUX ***', 1)
      end if

      call read_miprep(infile, ipixel, &
                       heightlev_in, presslev_in, templev_in, &
                       h2olev_in, co2lev_in, ch4lev_in, colev_in, &
                       meta, atm_scenario%surface_pressure, atm_scenario%surface_elevation, atm_scenario%surface_wspeed, ierr)
      if (ierr .ne. 0) then
         call writelog("READ_AUX: Error in read_miprep", 6)
         return
      end if

      !*** Deallocate
      if (allocated(atm_scenario%z)) deallocate (atm_scenario%z)
      if (allocated(atm_scenario%p)) deallocate (atm_scenario%p)
      if (allocated(atm_scenario%t)) deallocate (atm_scenario%t)
      if (allocated(atm_scenario%h2o)) deallocate (atm_scenario%h2o)
      if (allocated(atm_scenario%co2)) deallocate (atm_scenario%co2)
      if (allocated(atm_scenario%ch4)) deallocate (atm_scenario%ch4)
      if (allocated(atm_scenario%co)) deallocate (atm_scenario%co)
      ninput = size(heightlev_in)
      allocate ( &
         atm_scenario%z(ninput), &
         atm_scenario%p(ninput), &
         atm_scenario%t(ninput), &
         atm_scenario%h2o(ninput), &
         atm_scenario%co2(ninput), &
         atm_scenario%ch4(ninput), &
         atm_scenario%co(ninput))
      atm_scenario%z = heightlev_in
      atm_scenario%p = presslev_in
      atm_scenario%t = templev_in
      atm_scenario%h2o = h2olev_in
      atm_scenario%co2 = co2lev_in
      atm_scenario%ch4 = ch4lev_in
      atm_scenario%co = colev_in

      if (outputflag >= 2) then
         call writelog('*** End of READ_AUX ***', 1)
      end if

   end subroutine read_aux
!------------------------------------------------------------------------------

end module atmosphere_interface_module
