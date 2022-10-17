!------------------------------------------------------------------------------
!> @details This module contains routines that are only relevant for
!! synthetic measurements
module synthetic_input_module
   use header_module
   use read_errors_module, only: instrument_errors, read_errors
   use spectrum_internal_module, only: spectrum
   use forward_model_module, only: absorbers, window_ini, window_spectrum, aero
   use atmosphere_interface_module, only: atmospheric_scenario, read_atmosphere, read_atmosphere_nc_js
   use atmosphere_internal_module, only: atmosphere_interpolate, altitude_grid
   use retrieval_module, only: get_absorbers
   use auxiliary_routines_module, only: check
   use netcdf
   implicit none
   private

!*** Type
   public :: instrument_errors

!*** Procedures
  public :: get_scenario_info, get_scenario_info_nc_js, get_synthetic_data, shift_solar_spectrum, read_errors, get_synthetic_data_nc

!------------------------------------------------------------------------------
!*** Only for synthetic spectra:
   type, public :: synthetic_data
!> True state vector (dim: nstate)
      real(double), dimension(:), allocatable :: x_true
!> Instrument noise contribution (dim: nwave_lo)
      real(double), dimension(:), allocatable :: y_noise
!> Air column (true, retrieved)
      real(double) :: air_col_true
   end type synthetic_data

contains
!***********************************************************************
!
!    --->  SUBROUTINE COMMENTED OUT BY JS ON 2018-09-04  <---
!
!    get_scenario_info does probably not read in the data correctly.
!    No information on how many sets of albedo, OT, COT etc.
!
!    subroutine use_real_cirrus(meteo_file, aerosol_ini)
!    character(len=*), intent(in) :: meteo_file
!    type(aero), dimension(:), intent(inout) :: aerosol_ini
! !*** Local variables
!    integer :: k, ntype_aer

!    real(double) :: sza, iza, phi, alt1, alt2, reff, tilt, frac
!    real(double), dimension(:), allocatable :: alb, ot, cot

!  !*** Get comparison information from atmospheric scenario
!    call get_scenario_info(meteo_file, &
!         sza, iza, phi, &
!         alb, ot, cot, &
!         alt1, alt2, reff, tilt, frac)

!       ntype_aer = size(aerosol_ini)
!       do k = 1, ntype_aer
!          if (aerosol_ini(k)%CirrusFlag==1 ) then
!             if (aerosol_ini(k)%aeralt1 < -13000.d0) aerosol_ini(k)%aeralt1 = alt1
!             if (aerosol_ini(k)%aeralt2 < 0.d0 ) aerosol_ini(k)%aeralt2 = alt2
!             if (aerosol_ini(k)%shapefrac < 0.d0 )aerosol_ini(k)%shapefrac = frac
!             if (aerosol_ini(k)%reff < 0.d0) aerosol_ini(k)%reff = reff
!             if (aerosol_ini(k)%tilt_angle < 0.d0) aerosol_ini(k)%tilt_angle = tilt
!          endif
!       enddo

!    end subroutine
!
!            --->           !!!!!         <---
!------------------------------------------------------------------------------

   subroutine get_scenario_info(meteo_file, &
                                sza, iza, phi, &
                                wave_start, wave_stop, &
                                alb, ot, cot, &
                                alt1, alt2, reff, tilt, frac)

      character(len=*), intent(in) :: meteo_file

      real(double), intent(out) :: sza, iza, phi, alt1, alt2, reff, tilt, frac
      real(double), intent(out), dimension(:) :: wave_start, wave_stop, alb, ot, cot

      character(len=len(meteo_file)) :: data_file
      real(double) :: saz, iaz
      integer :: ninput, k, io, i, nwin_in

      i = index(meteo_file, '/', BACK=.TRUE.)
      if (meteo_file(i + 1:i + 6) .EQ. "METEO_") then
         data_file = trim(meteo_file(1:i))//'ATM_'//trim(meteo_file(i + 7:))
         print *, "WARNING: no information about date, coordinates and angles in "//trim(meteo_file)
         print *, "         -> Reading data from "//trim(data_file)//" intead!"
      else
         data_file = meteo_file
      end if

      !$OMP critical
      open (newunit(io), FILE=trim(data_file), FORM='FORMATTED')
      !$OMP end critical
      read (io, *) ninput  ! Read number of data levels
      do k = 1, ninput + 9 ! Skip reading atmpsheric data, date, time, and coordinates
         read (io, *)
      end do

      read (io, *) sza, iza, saz, iaz ! Read geometry
      phi = abs(iaz - saz)

      read (io, *) ! Skip empy line
      read (io, *) nwin_in
      read (io, *), wave_start, wave_stop, alb, ot, cot, alt1, alt2, reff, tilt, frac ! Read albedo, OT, COT and aerosol/cirrusx(?) properties

      close (io)

   end subroutine get_scenario_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

   subroutine get_scenario_info_nc_js(infile, alb, ot, cot)

      character(len=*), intent(in) :: infile
      real(double), intent(out), dimension(:) :: alb, ot, cot
      !*** local variables
      integer :: k, l, i, n, nwin, ierr
      integer :: ncid, grpid(3), varid
      integer :: sx, sy, start2d(2)
      character(stringlen) :: meteo_file, index_info, message
      real(double) :: var

      i = INDEX(infile, '.nc')
      index_info = infile(i + 4:)
      meteo_file = trim(infile(:i + 2))

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      !*** Read meteo-file
      call check(nf90_open(trim(meteo_file), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         call writelog('GET_SCENARIO_INFO_NC_JS: error opening file: '//trim(meteo_file), 6)
         return
      end if

      call check(nf90_inq_grps(ncid, nwin, grpid), ierr)
      if (ierr .ne. 0) return

      start2d = (/sx, sy/)

      do n = 1, nwin
         !*** Get surface_albedo
         call check(NF90_INQ_VARID(grpid(n), "alb_inp", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start2d), ierr)
         if (ierr .ne. 0) return
         if (var .LT. 1e-10) then
            var = 0.d0
         end if
         alb(n) = var

         !*** Get aerosol_optical_thickness
         call check(NF90_INQ_VARID(grpid(n), "ot_inp", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start2d), ierr)
         if (ierr .ne. 0) return
         if (var .LT. 1e-10) then
            var = 0.d0
         end if
         ot(n) = var

         !*** Get surface_albedo
         call check(NF90_INQ_VARID(grpid(n), "cot_inp", varid), ierr)
         if (ierr .ne. 0) return
         call check(NF90_GET_VAR(grpid(n), varid, var, start=start2d), ierr)
         if (ierr .ne. 0) return
         if (var .LT. 1e-10) then
            var = 0.d0
         end if
         cot(n) = var

      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), ierr)
      if (ierr .ne. 0) return

   end subroutine get_scenario_info_nc_js

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!*** Get true CH4 profile and noise for comparison from NetCDF files
!------------------------------------------------------------------------------
   subroutine get_synthetic_data_nc( &
      win_ini, atm_scenario, aerosol_ini, &
      measurement, meta, &
      grid, Tflag, fitflag, syn_output)
!*** input
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(aero), dimension(:), intent(in) :: aerosol_ini
      type(spectrum), dimension(:), intent(in) :: measurement
      type(metadata), intent(in) :: meta
      type(altitude_grid), intent(in) :: grid
      integer, intent(in) :: Tflag, fitflag
!*** input/output
      type(atmospheric_scenario), intent(in) :: atm_scenario
!*** output
      type(synthetic_data), intent(inout) :: syn_output
!*** local
      integer :: i, k, n, nwin, naux, nstate, naer, ntype_aer, ierr
      real(double):: true_col, false_col, norm
      real(double) ::  z_tropopause, z_bl
      real(double), dimension(:), allocatable :: dvair
      real(double), dimension(4) :: errors_true
      type(window_spectrum), dimension(:), allocatable :: win_true
      type(absorbers) :: absorb
      type(atmosphere) :: atm_xs, atm_rt, atm_retr

      character(stringlen) :: meteo_file
      integer, parameter :: outputflag = 0
!------------------------------------------------------------------------------
!*** Allocate local win
      nwin = size(win_ini)
      allocate (win_true(nwin))

!*** Interpolate input meteo data (to get win)
      call atmosphere_interpolate( &
         atm_scenario, &
         grid, &
         meta%surface_elevation, &
         meta%lat(1), &
         win_ini, &
         outputflag, &
         win_true, &
         atm_xs, &
         atm_rt, &
         atm_retr, &
         dvair, &
         z_tropopause, &
         z_bl, &
         ierr)

!*** Get true air column:
      syn_output%air_col_true = sum(dvair)

!*** Put absorber in type
      call get_absorbers(win_ini, absorb, ierr)

!*** Calculate state vector dimensions
      naux = absorb%ntype_global
      do n = 1, nwin
         naux = naux + win_ini(n)%albflag + win_ini(n)%IoffFlag &
                + abs(win_ini(n)%spsh0flag) + abs(win_ini(n)%spsh1flag) + abs(win_ini(n)%spsh2flag) + abs(win_ini(n)%sunsh0flag)
      end do
      ntype_aer = size(aerosol_ini)
      naer = 0
      do n = 1, ntype_aer
         naer = naer + sum(aerosol_ini(n)%AerosolFlags)
      end do
      naux = naux + TFlag + naer
      nstate = grid%nlay*absorb%ntype_target + naux
      if (allocated(syn_output%x_true)) deallocate (syn_output%x_true)
      if (allocated(syn_output%y_noise)) deallocate (syn_output%y_noise)

      allocate (syn_output%x_true(nstate))
      do n = 1, nwin
         win_true(n)%nwave_lo = 0
         do i = 1, size(measurement)
            do k = 1, measurement(i)%nwave
               if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                   measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop) then
                  win_true(n)%nwave_lo = win_true(n)%nwave_lo + 1
               end if
            end do
         end do
      end do
      allocate (syn_output%y_noise(sum(win_true(:)%nwave_lo)))

      call get_xtrue(atm_xs%n, grid%nlay, fitflag, &
                     absorb, measurement, &
                     syn_output%x_true, syn_output%y_noise, win_ini, win_true)

   end subroutine get_synthetic_data_nc
!------------------------------------------------------------------------------
!*** Get true CH4 profile and noise for comparison
!*** Use perturbed prior meteo data (not truth)
!------------------------------------------------------------------------------
   subroutine get_synthetic_data( &
      win_ini, atm_scenario, meteo_errors, aerosol_ini, &
      spec_dir, atm_dir, filename, measurement, &
      grid, Tflag, fitflag, atmflag, syn_output, ierr)
!*** input
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(aero), dimension(:), intent(in) :: aerosol_ini
      character(len=*) :: spec_dir, atm_dir, filename
      type(spectrum), dimension(:), intent(in) :: measurement
      type(altitude_grid), intent(in) :: grid
      integer, intent(in) :: Tflag, fitflag, atmflag
      real(double), dimension(:), intent(in) :: meteo_errors
!*** input/output
      type(atmospheric_scenario), intent(inout) :: atm_scenario
!*** output
      integer, intent(out) :: ierr
      type(synthetic_data), intent(inout) :: syn_output
!*** local
      integer :: i, k, n, nwin, naux, nstate, naer, ntype_aer
      real(double):: true_col, false_col, norm
      real(double) ::  z_tropopause, z_bl
      real(double), dimension(:), allocatable :: dvair
      real(double), dimension(4) :: errors_true
      type(window_spectrum), dimension(:), allocatable :: win_true, win
      type(atmospheric_scenario) :: scenario_true
      type(absorbers) :: absorb
      type(atmosphere) :: atm_xs, atm_rt, atm_retr
      type(metadata) :: meta
      character(stringlen) :: meteo_file
      integer, parameter :: outputflag = 0
!------------------------------------------------------------------------------
!*** Allocate local win
      nwin = size(win_ini)
      if (allocated(win_true)) deallocate (win_true)
      if (allocated(win)) deallocate (win)
      allocate (win_true(nwin), win(nwin))

!*** No errors when computing true CH4 profile
      errors_true = 0.d0

!*** File where true CH4 profile is stored
      meteo_file = trim(spec_dir)//trim(filename)

!*** Optionally for comparison purposes:
!*** for retrieval with cirrus, use the same parameters (reff, aeralt1, aeralt2, shapefrac)
!*** as used for creating the synthetic spectrum (only if parameters are negative in retrieval.ini)
!    call use_real_cirrus(meteo_file, aerosol_ini)

      if (atmflag .EQ. 4) then
         call read_atmosphere_nc_js(meteo_file, outputflag, scenario_true, meta, ierr, errors_true)
      else
         call read_atmosphere(meteo_file, outputflag, scenario_true, meta, ierr, errors_true)
      end if
      if (ierr .ne. 0) return

!*** Interpolate input meteo data (to get win)
      call atmosphere_interpolate( &
         atm_scenario, &
         grid, &
         meta%surface_elevation, &
         meta%lat(1), &
         win_ini, &
         outputflag, &
         win, &
         atm_xs, &
         atm_rt, &
         atm_retr, &
         dvair, &
         z_tropopause, &
         z_bl, &
         ierr)
      if (ierr .ne. 0) return

!*** Interpolate true meteo data (to get win_true)
      call atmosphere_interpolate( &
         scenario_true, &
         grid, &
         meta%surface_elevation, &
         meta%lat(1), &
         win_ini, &
         outputflag, &
         win_true, &
         atm_xs, &
         atm_rt, &
         atm_retr, &
         dvair, &
         z_tropopause, &
         z_bl, &
         ierr)
      if (ierr .ne. 0) return
!*** Get true air column:
      syn_output%air_col_true = sum(dvair)

      !*** Perturb CH4/H2O/T profiles
      if (spec_dir .ne. atm_dir) then
         call writelog('GET_SYNTHETIC_DATA: meteo data is different from truth for error sensitivity study', 4)
         !*** First normalize CH4 profile so that total column is equal to truth
         !*** (but profile is not if latitudinal mean is used),
         !*** then introduce scaling error
         do n = 1, nwin
            do i = 1, win_ini(n)%ntype
               if (win_ini(n)%type_x(i) == 6) then
                  false_col = sum(win(n)%dv_x(:, i))
                  true_col = sum(win_true(n)%dv_x(:, i))
                  norm = true_col/false_col
                  atm_scenario%ch4 = atm_scenario%ch4*norm*(1.d0 + 0.01*meteo_errors(1))
               end if
            end do
         end do

         !*** First normalize H2O profile so that total column is equal to truth
         !*** (but profile is not if latitudinal mean is used),
         !*** then introduce scaling error
         do n = 1, nwin
            do i = 1, win_ini(n)%ntype
               if (win_ini(n)%type_x(i) == 1) then
                  false_col = sum(win(n)%dv_x(:, i))
                  true_col = sum(win_true(n)%dv_x(:, i))
                  norm = true_col/false_col
                  atm_scenario%h2o = atm_scenario%h2o*norm*(1.d0 + 0.01*meteo_errors(2))
               end if
            end do
         end do

         !*** First shift Temperature profile so that Tsurf is equal to truth
         !*** (but profile is not if latitudinal mean is used),
         !*** then introduce temperature offset
         n = size(scenario_true%t)
         norm = scenario_true%t(n) - atm_scenario%t(n)
         atm_scenario%t(:) = atm_scenario%t(:) + norm + meteo_errors(4)
      end if

!*** Put absorber in type
      call get_absorbers(win_ini, absorb, ierr)

!*** Calculate state vector dimensions
      naux = absorb%ntype_global
      do n = 1, nwin
         naux = naux + win_ini(n)%albflag + win_ini(n)%IoffFlag &
                + abs(win_ini(n)%spsh0flag) + abs(win_ini(n)%spsh1flag) + abs(win_ini(n)%spsh2flag) + abs(win_ini(n)%sunsh0flag)
      end do
      ntype_aer = size(aerosol_ini)
      naer = 0
      do n = 1, ntype_aer
         naer = naer + sum(aerosol_ini(n)%AerosolFlags)
      end do
      naux = naux + TFlag + naer
      nstate = grid%nlay*absorb%ntype_target + naux
      if (allocated(syn_output%x_true)) deallocate (syn_output%x_true)
      if (allocated(syn_output%y_noise)) deallocate (syn_output%y_noise)
!      if(allocated(syn_output%y_cov)) deallocate(syn_output%y_cov)
      allocate (syn_output%x_true(nstate))
      do n = 1, nwin
         win_true(n)%nwave_lo = 0
         do i = 1, size(measurement)
            do k = 1, measurement(i)%nwave
               if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                   measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop) then
                  win_true(n)%nwave_lo = win_true(n)%nwave_lo + 1
               end if
            end do
         end do
      end do
      allocate (syn_output%y_noise(sum(win_true(:)%nwave_lo)))
!      allocate(syn_output%y_cov(sum(win_true(:)%nwave_lo)))

      call get_xtrue(atm_xs%n, grid%nlay, fitflag, &
                     absorb, measurement, &
                     syn_output%x_true, syn_output%y_noise, win_ini, win_true)

   end subroutine get_synthetic_data
!------------------------------------------------------------------------------
!*** Get true state vector (x_true) and noise (y_noise)
   subroutine get_xtrue(natm, &
                        nlay, &
                        fitflag, &
                        absorb, &
                        measurement, &
                        x_true, &
                        y_noise, &
                        win_ini, &
                        win)
!*** input
      integer, intent(in) :: natm, nlay, fitflag
      type(absorbers), intent(in) :: absorb
      type(spectrum), dimension(:), intent(in) :: measurement
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(window_spectrum), dimension(:), intent(in) :: win
!*** output
      real(double), dimension(:), intent(out) :: x_true     ! True state vector (Dim: nstate)
      real(double), dimension(:), intent(out) :: y_noise  ! Instrument noise contribution (Dim: nwave_lo)
!*** local variables
      integer :: i, j, k, l, n, off, nwin
!--------------------------------------------------------------------------

      nwin = size(win_ini)
      l = natm/nlay
!*** Write true state vector
!*** Restrict x_true only to trace gases (no albedo etc)
      x_true = 0.D0

!*** Target absorber profiles
      off = 0
      do j = 1, absorb%ntype_target
         n = 1
         do while (n .le. nwin)
            i = 1
            do while (i .le. win_ini(n)%ntype)
               if (win_ini(n)%type_x(i) == absorb%type_x_target(j)) then
                  do k = 1, nlay
                     x_true(k + off) = sum(win(n)%dv_x((k - 1)*l + 1:k*l, i))
                  end do
                  i = win_ini(n)%ntype
                  n = nwin
               end if
               i = i + 1
            end do
            n = n + 1
         end do
         off = off + nlay
      end do

!*** Interfering absorber columns
      off = nlay*absorb%ntype_target
      do j = 1, absorb%ntype_global
         n = 1
         do while (n .le. nwin)
            i = 1
            do while (i .le. win_ini(n)%ntype)
               if (win_ini(n)%type_x(i) == absorb%type_x_global(j)) then
                  x_true(off + j) = sum(win(n)%dv_x(:, i))
                  i = win_ini(n)%ntype
                  n = nwin
               end if
               i = i + 1
            end do
            n = n + 1
         end do
      end do
      off = nlay*absorb%ntype_target + absorb%ntype_global + 1

!*** Instrument noise y_noise, only relevant for synthetic spectra
      off = 1
      do n = 1, nwin
         do i = 1, size(measurement)
            do k = 1, measurement(i)%nwave
               if (measurement(i)%wavelength(k) .ge. win_ini(n)%wave_start .and. &
                   measurement(i)%wavelength(k) .le. win_ini(n)%wave_stop) then
                  if (fitflag == 1) then !fitting reflectance
                     y_noise(off) = measurement(i)%radiance_error(k)/measurement(i)%irradiance(k)
                  elseif (fitflag == 2) then !fitting radiance
                     y_noise(off) = measurement(i)%radiance_error(k)
                  end if
                  off = off + 1
               end if
            end do
         end do
      end do

   end subroutine get_xtrue

!------------------------------------------------------------------------------

   subroutine shift_solar_spectrum(dlambda, wavelength, sun_spectrum)
!*** input
      real(double), intent(in) :: dlambda
      real(double), dimension(:), intent(in) ::  wavelength
!*** input/output
      real(double), dimension(:), intent(inout) :: sun_spectrum
!*** local
      integer :: k, nwave
      real(double), dimension(:), allocatable ::  wavelength_per, sun_spectrum_per

      nwave = size(wavelength)
      allocate (wavelength_per(nwave), sun_spectrum_per(nwave))
      do k = 1, nwave
         wavelength_per(k) = wavelength(k) + dlambda
      end do

      call spline_interpol(wavelength_per, sun_spectrum, nwave, &
                           wavelength, sun_spectrum_per, nwave)

      sun_spectrum = sun_spectrum_per

   end subroutine shift_solar_spectrum
!------------------------------------------------------------------------------

end module synthetic_input_module
