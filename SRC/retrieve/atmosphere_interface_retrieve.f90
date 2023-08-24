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
   public :: read_atm

contains

   subroutine read_atm(infile, outputflag, atm_scenario, meta, ierr, meteo_errors)
      !*** Input
      character(len=*), intent(in) :: infile
      integer, intent(in) :: outputflag
      real(double), dimension(4), intent(in), optional :: meteo_errors
      !*** Output
      type(atmospheric_scenario), intent(out) :: atm_scenario
      type(metadata), intent(out) :: meta
      integer, intent(out) :: ierr
      !*** Local variables
      integer :: i, k, n, nlevel
      integer :: ncid, varid, dimid_lat, dimid_lon, dimid_time, dimid_z
      integer :: sx, sy, start1d(1), start2d(2), start3d(3)
      integer :: time_id, lon_id, lat_id, elev_id, z_id
      character(stringlen) :: meteo_file, index_info, message
!-----------------------------------------------------

      i = INDEX(infile, '.nc')
      index_info = infile(i + 4:)
      meteo_file = trim(infile(:i + 2))

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      if (outputflag >= 2) then
         call writelog('*** Start of READ_ATM ***', 1)
         call writelog('ATM file: '//trim(meteo_file), 1)
      end if

      !*** Read meteodata
      call check(nf90_open(trim(meteo_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('READ_ATMOSPHERE: error opening file: '//trim(meteo_file), 6)
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

      !*** Get number of vertical levels
      call check(nf90_inq_dimid(ncid, "level", dimid_z), ierr)
      if (ierr .ne. 0) return
      call check(nf90_inquire_dimension(ncid, dimid_z, len=nlevel), ierr)
      if (ierr .ne. 0) return

      !*** Deallocate/allocate
      if (allocated(atm_scenario%p)) deallocate (atm_scenario%p)
      if (allocated(atm_scenario%t)) deallocate (atm_scenario%t)
      if (allocated(atm_scenario%z)) deallocate (atm_scenario%z)
      if (allocated(atm_scenario%h2o)) deallocate (atm_scenario%h2o)
      if (allocated(atm_scenario%co2)) deallocate (atm_scenario%co2)
      if (allocated(atm_scenario%ch4)) deallocate (atm_scenario%ch4)
      if (allocated(atm_scenario%co)) deallocate (atm_scenario%co)
      allocate ( &
         atm_scenario%z(nlevel), &
         atm_scenario%p(nlevel), &
         atm_scenario%t(nlevel), &
         atm_scenario%h2o(nlevel), &
         atm_scenario%co2(nlevel), &
         atm_scenario%ch4(nlevel), &
         atm_scenario%co(nlevel))

      call netcdf_get_vector_var(ncid, "pressure", atm_scenario%p, start3d)
      ! internally, RemoTeC works with hPa, pressure provided in Pa
      atm_scenario%p = atm_scenario%p / 100
      call netcdf_get_vector_var(ncid, "temperature", atm_scenario%t, start3d)
      call netcdf_get_vector_var(ncid, "geometric_altitude", atm_scenario%z, start3d)

      ! get trace gases only if they are in the netcdf file
      call netcdf_get_vector_var(ncid, "h2o", atm_scenario%h2o, start3d)
      if (nf90_inq_varid(ncid, "co2", i) /= nf90_enotvar) then
         call netcdf_get_vector_var(ncid, "co2", atm_scenario%co2, start3d)
      end if
      if (nf90_inq_varid(ncid, "ch4", i) /= nf90_enotvar) then
         call netcdf_get_vector_var(ncid, "ch4", atm_scenario%ch4, start3d)
      end if
      if (nf90_inq_varid(ncid, "co", i) /= nf90_enotvar) then
         call netcdf_get_vector_var(ncid, "co", atm_scenario%co, start3d)
      end if

      meta%surface_elevation = atm_scenario%z(nlevel)
      atm_scenario%surface_pressure = atm_scenario%p(nlevel)
      atm_scenario%surface_elevation = atm_scenario%z(nlevel)

      atm_scenario%surface_wspeed = 0.d0

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)

      if (present(meteo_errors)) then
         !*** Introduce error in pressure profile
         atm_scenario%p = (1.d0 + 0.01*meteo_errors(3))*atm_scenario%p
         atm_scenario%t(:) = atm_scenario%t(:) + meteo_errors(4)
      end if

100   if (ierr .ne. 0) then
         if (outputflag >= 2) then
            write (message, '(a)') 'READ_ATM: Error opening/reading meteo_file '//trim(meteo_file)
            call writelog(message, 6)
         end if
         return
      end if

      if (outputflag >= 2) then
         call writelog('*** End of READ_ATM ***', 1)
      end if

   end subroutine read_atm



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



end module atmosphere_interface_module
