module diagnostics_module
   use header_module
   use aerosol_input_module, only: aero
   use synthetic_input_module, only: synthetic_data, get_scenario_info, get_scenario_info_nc_js
   use retrieval_module, only: retrieval_data, window_ini
   use netcdf
   use spectrum_interface_module
   use auxiliary_routines_module, only: check
   implicit none

!Errorflags
!           iter .LE. maxiter -> OK, run converged
!           iter = 90         -> MaxOT (currently set to tau = 1) has been exceeded during retrieval
!!          iter = 91         -> aerosol boundary hit
!           iter = 92         -> maxiter has been exceeded
!           iter = 93         -> maxiter has been exceeded and lambda>0 OR lambda >1d5
!           iter = 95         -> failure in Cross-section table
!           iter = 96         -> chi2 > 1.d10 or NaN
!           iter = 97         -> SVDFlag < 0), i.e., singular value decomposition failed
!           iter = 98         -> degree of freedom = NaN
!           iter = 99         -> boundary hit

contains

!------------------------------------------------------------------------------

   subroutine diagnostics_retrieve( &
      output_dir, runid, &
      debug_output_flag, &
      spectrum_file, meteo_file, &
      win_ini, &
      syn_output, &
      retrieval_output, &
      meta, &
      nlay, nwin_in, outputflag)

      !*** arguments
      character(len=*), intent(in) :: output_dir
      integer, intent(in) :: runid
      integer, intent(in) :: debug_output_flag, nlay, nwin_in, outputflag
      character(len=*), intent(in) :: spectrum_file, meteo_file
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(synthetic_data), intent(in) :: syn_output
      type(retrieval_data), intent(inout) :: retrieval_output
      type(metadata), intent(in) :: meta
      !*** local variables
      integer :: k, l, i, j, i1, i2, length, now(3), io, nstate, nwin, ntype_target, ntype_target_in = 0
      character(299) :: filename
      character(6):: runidstring
      character(8):: lengthid
      character(20) :: format
      real(double) :: scale  !res, noise_term, cov_term, VMR_true, VMR_bias
      real(double) :: sza, iza, phi, aer_alt1, aer_alt2, aer_reff, aer_tilt, aer_frac
      real(double), dimension(nwin_in) :: wave_start_in, wave_stop_in, alb_in, ot_in, cot_in
      real(double), dimension(size(retrieval_output%type_x_target)) :: x, x_err
      real(double), dimension(11) :: x_in
      integer, dimension(11) :: x_in_flag
      logical :: exst
      character(299) :: formatstring
      character(99) :: s_nwin, s_nwin_in, s_ntype_target, s_ntype_target_in
      character(99), dimension(size(retrieval_output%type_x_target)) :: x_name
      character(99), dimension(11) :: x_in_name
      real(double) :: x_air
      real(double) :: psf
      !---------------------------------------------------------------------------------

      x_in_flag(:) = 0

      x_in_name(1) = "H2O"
      x_in_name(2) = "CO2"
      x_in_name(3) = "O3"
      x_in_name(4) = "N2O"
      x_in_name(5) = "CO"
      x_in_name(6) = "CH4"
      x_in_name(7) = "O2"
      x_in_name(8) = "NO"
      x_in_name(9) = "SO2"
      x_in_name(10) = "NO2"
      x_in_name(11) = "NH3"

      ! MEASURE TIME NOW AND TOTAL ELAPSED TIME AT THE END OF THE RUN:
      call itime(now)

      write (runidstring, '(I6.6)') runid
      filename = trim(output_dir)//'selected_output_'//runidstring//'.dat'

      if (outputflag >= 2) then
         write (*, '(a)')
         write (*, '(a)') '*****************************************'
         write (*, '(a)') '***** Diagnostic output'
         write (*, '(a)') '*****'
         write (*, *) trim(filename)
         write (*, *) ''
      end if

      nwin = size(win_ini)
      ntype_target = size(retrieval_output%type_x_target)
      nstate = size(retrieval_output%x_state(:, 1))

      !*** DRY AIR MOLE FRACTION
      x_air = sum(retrieval_output%dvair)
      do j = 1, ntype_target

         if (abs(retrieval_output%type_x_target(j)) == 1) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) == 2) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) == 5) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) == 6) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) == 7) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) == 11) then
            scale = 1.E9
            k = 11
         elseif (abs(retrieval_output%type_x_target(j)) >= 100 .and. abs(retrieval_output%type_x_target(j)) <= 199) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) >= 200 .and. abs(retrieval_output%type_x_target(j)) <= 299) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) >= 500 .and. abs(retrieval_output%type_x_target(j)) <= 599) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) >= 600 .and. abs(retrieval_output%type_x_target(j)) <= 699) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) >= 700 .and. abs(retrieval_output%type_x_target(j)) <= 799) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) >= 1100 .and. abs(retrieval_output%type_x_target(j)) <= 1199) then
            scale = 1.E9
            k = 11
         end if

         i1 = 1 + (nlay*(j - 1))
         i2 = nlay*j

         x(j) = sum(retrieval_output%x_state(i1:i2, retrieval_output%iter))/sum(retrieval_output%dvair)*scale
         x_err(j) = dsqrt(sum(retrieval_output%s_state(i1:i2, i1:i2)))/sum(retrieval_output%dvair)*scale ! What is this, how is it calculated?

         l = len(trim(retrieval_output%x_state_name(i1)))
         x_name(j) = retrieval_output%x_state_name(i1) (l - 7:l) ! extract the last seven cahracters (=TYPEXXXX)

         if (x_in_flag(k) .NE. 1) then
            x_in_flag(k) = 1
            x_in(k) = sum(syn_output%x_true(i1:i2))/syn_output%air_col_true*scale
         end if
          !!TODO: Possibly include a check that current and exeisting values are identical if x_in_flag(k)==1

      end do

      !*** Surface pressure
      psf = retrieval_output%p(nlay + 1)

      do k = 1, 11
         if (x_in_flag(k) .EQ. 1) ntype_target_in = ntype_target_in + 1
      end do

      ! ---------------------------------------------------------- !

      !*** FOR CONVERGED RETRIEVALS ERROR_ID GIVES THE NUMBER OF ITERATIONS
      if (retrieval_output%error_ID == 0) retrieval_output%error_ID = retrieval_output%iter

      ! ---------------------------------------------------------- !

      ! LOAD REFERENCE (INPUT) DATA
      call get_scenario_info(meteo_file, &
                             sza, iza, phi, &
                             wave_start_in, wave_stop_in, &
                             alb_in, ot_in, cot_in, &
                             aer_alt1, aer_alt2, aer_reff, aer_tilt, aer_frac)

      ! ---------------------------------------------------------- !

      ! WRITE HEADER TO FILE (IF NO FILE EXISTS)

      inquire (FILE=trim(filename), EXIST=exst)
      if (.not. exst) then
         open (newunit(io), FILE=trim(filename), STATUS='NEW')

         write (io, '(A)') '# CHAR  SPECTRUM_FILE'
         write (io, '(A)') '# INT   YEAR'
         write (io, '(A)') '# INT   MONTH'
         write (io, '(A)') '# INT   DAY'
         write (io, '(A)') '# FLOAT LAT'
         write (io, '(A)') '# FLOAT LON'
         write (io, '(A)') '# FLOAT SZA'
         write (io, '(A)') '# FLOAT AIRMASS'
         write (io, '(A)') '# INT   FLAG_CONVERGENCE'
         write (io, '(A)') '# INT   NUMBER_OF_ITERATIONS'
         write (io, '(A)') '# INT   ERROR_ID'

         write (io, '(A)') '# FLOAT CHI2'
         write (io, '(A)') '# FLOAT DFS'
         write (io, '(A)') '# FLOAT DFS_TARGET'
         write (io, '(A)') '# FLOAT DFS_SCAT'

         write (io, '(A)') '# INT   NWIN'

         do k = 1, nwin
            write (io, '(A,I2.2)') '# FLOAT WAVE_START_BAND', k
            write (io, '(A,I2.2)') '# FLOAT WAVE_STOP_BAND', k
            write (io, '(A,I2.2)') '# FLOAT ALBEDO_BAND', k
            write (io, '(A,I2.2)') '# FLOAT OT_BAND', k
            write (io, '(A,I2.2)') '# FLOAT COT_BAND', k
         end do

         do k = 1, ntype_target
            write (io, '(A,A)') '# FLOAT X_', trim(x_name(k))
         end do

         do k = 1, ntype_target
            write (io, '(A,A,A)') '# FLOAT X_', trim(x_name(k)), '_ERR'
         end do

         ! REFERENCE (INPUT) DATA
         write (io, '(A)') '# INT   NWIN_IN'

         do k = 1, nwin_in
            write (io, '(A,I2.2)') '# FLOAT WAVE_START_IN_BAND', k
            write (io, '(A,I2.2)') '# FLOAT WAVE_STOP_IN_BAND', k
            write (io, '(A,I2.2)') '# FLOAT ALBEDO_IN_BAND', k
            write (io, '(A,I2.2)') '# FLOAT OT_IN_BAND', k
            write (io, '(A,I2.2)') '# FLOAT COT_IN_BAND', k
         end do

         do k = 1, 11
            if (x_in_flag(k) .EQ. 1) then
               write (io, '(A,A,A)') '# FLOAT X_', trim(x_in_name(k)), '_IN'
            end if
         end do

         close (io)
      end if

      ! ---------------------------------------------------------- !

      ! Define format string
      write (s_nwin, "(I3)"), nwin
      write (s_nwin_in, "(I3)"), nwin_in
      write (s_ntype_target, "(I3)"), ntype_target
      write (s_ntype_target_in, "(I3)"), ntype_target_in

      formatstring = '(A,x,'                                                       !! Spectrum_file
      formatstring = trim(formatstring)//'3(I5,x),'                                !! YYYY, MM, DD
      formatstring = trim(formatstring)//'5(1pE16.8E3,x),'                         !! LAT, LON, SZA, AIRMASS, PSF
      formatstring = trim(formatstring)//'3(I5,x),'                                !! FLAG_CONVERGENCE, NUMBER_OF_ITERATIONS, ERROR_ID
      formatstring = trim(formatstring)//'4(1pE16.8E3,x),'                         !! CHI2, DFS, DFS_TARGET, DFS_SCAT
      formatstring = trim(formatstring)//'1(I5,x),'                                !! NWIN

      formatstring = trim(formatstring)//trim(s_nwin)//'(1pE16.8E3,x),'            !! First wavelength (per band)
      formatstring = trim(formatstring)//trim(s_nwin)//'(1pE16.8E3,x),'            !! Last wavelength (per band)
      formatstring = trim(formatstring)//trim(s_nwin)//'(1pE16.8E3,x),'            !! Retrieved albedo (per band)
      formatstring = trim(formatstring)//trim(s_nwin)//'(1pE16.8E3,x),'            !! Retrieved OT (per band)
      formatstring = trim(formatstring)//trim(s_nwin)//'(1pE16.8E3,x),'            !! Retrieved COT (per band)
      formatstring = trim(formatstring)//trim(s_ntype_target)//'(1pE16.8E3,x),'    !! Retrieved dry air mole fraction (per target absorbers)
      formatstring = trim(formatstring)//trim(s_ntype_target)//'(1pE16.8E3,x),'    !! Retrieved dry air mole fraction error (per target absorbers)

      formatstring = trim(formatstring)//'1(I5,x),'                                !! NWIN_IN

      formatstring = trim(formatstring)//trim(s_nwin_in)//'(1pE16.8E3,x),'         !! Input first wavelength (per band)
      formatstring = trim(formatstring)//trim(s_nwin_in)//'(1pE16.8E3,x),'         !! Input last wavelength (per band)
      formatstring = trim(formatstring)//trim(s_nwin_in)//'(1pE16.8E3,x),'         !! Input albedo (per band)
      formatstring = trim(formatstring)//trim(s_nwin_in)//'(1pE16.8E3,x),'         !! Input OT (per band)
      formatstring = trim(formatstring)//trim(s_nwin_in)//'(1pE16.8E3,x),'         !! Input COT (per band)
      formatstring = trim(formatstring)//trim(s_ntype_target_in)//'(1pE16.8E3,x)'  !! Input dry air mole fraction (per target absorbers)

      formatstring = trim(formatstring)//')'

      ! ---------------------------------------------------------- !

      ! Write selected output to file

      open (newunit(io), FILE=trim(filename), POSITION='APPEND')

      write (io, formatstring) &
         trim(spectrum_file), &
         meta%time(1), &
         meta%time(2), &
         meta%time(3), &
         meta%lat(1), &
         meta%lon(1), &
         sza, &
         x_air, &
         psf, &
         retrieval_output%convergence, &
         retrieval_output%iter, &
         retrieval_output%error_ID, &
         retrieval_output%chi2(retrieval_output%iter)/(sum(retrieval_output%ny) - retrieval_output%dfs), &
         retrieval_output%dfs, &
         retrieval_output%dfs_target, &
         retrieval_output%dfs_scat, &
         nwin, &
            (win_ini(i)%wave_start,win_ini(i)%wave_stop, retrieval_output%albedo(i), retrieval_output%ot(i), retrieval_output%cot(i), i=1,nwin), &
         x(:), &
         x_err(:), &
         nwin_in, &
         (wave_start_in(i), wave_stop_in(i), alb_in(i), ot_in(i), cot_in(i), i=1, nwin_in), &
         x_in(pack([(i, i=1, size(x_in_flag))], x_in_flag == 1))

      close (io)

      if (outputflag >= 2) then
         write (*, '(a)') '***** Done.'
         write (*, 80) now
80       format('***** Time at the end: ', i2.2, ':', i2.2, ':', i2.2)
         write (*, '(a)') '*****************************************'
      end if

   end subroutine diagnostics_retrieve

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

   subroutine diagnostics_retrieve_nc_js( &
      output_dir, runid, &
      meteo_file, &
      win_ini, &
      syn_output, &
      retrieval_output, &
      meta, &
      measurement, &
      nlay, outputflag)
      !*** arguments
      character(len=*), intent(in) :: output_dir
      integer, intent(in) :: runid
      character(len=*), intent(in) :: meteo_file
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(synthetic_data), intent(in) :: syn_output
      type(retrieval_data), intent(inout) :: retrieval_output
      type(metadata), intent(in) :: meta
      type(spectrum), dimension(:), intent(in) :: measurement
      integer, intent(in) :: nlay, outputflag
      !*** local variables
      integer :: nwin
      integer :: k, l, n, i, j, i1, i2, length, now(3), io, nstate, ntype_target, ntype_target_in
      character(299) :: ncfile
      character(6):: runidstring
      real(double) :: scale, chi2
      real(double), dimension(:), allocatable :: x, x_err, wave_start_in, wave_stop_in, alb_in, ot_in, cot_in
      real(double), dimension(11) :: x_in
      integer, dimension(11) :: x_in_flag
      character(99), dimension(11) :: x_in_name, x_in_unit
      logical :: exst
      character(99) :: unit
      character(99), dimension(:), allocatable :: x_name, x_unit
      integer :: sx, sy
      integer :: ncid, ierr, stat, ngroup, dimid_nobs, dimid_wave, lastindex, dimid_time, nwave
      integer :: sza_id, vza_id, xair_id, psf_id, lat_id, lon_id, time_id, x_id, y_id
      integer :: cf_id, chi_id, dfs_id, dfss_id, eid_id, it_id
      integer, dimension(:), allocatable :: grpid, wave_id, ot_id, cot_id, alb_id, otin_id, cotin_id, albin_id, tc_id, dfst_id, tcerr_id, tcin_id
      character*1 :: ch
      character(stringlen) :: group_name, index_info
      real(double) :: x_air
      real(double) :: psf
      !-----------------------------------------------------------

      nwin = size(win_ini)

      if (nwin .NE. size(measurement)) then
         call writelog('DIAGNOSTICS_SIM_NC_JS: expected same number of windows for measurement and retrieval', 8)
         return
      end if

      !---------------------------------------------------------------------------------
      x_in_flag(:) = 0
      x_in_name(1) = "h2o"
      x_in_name(2) = "co2"
      x_in_name(3) = "o3"
      x_in_name(4) = "n2o"
      x_in_name(5) = "co"
      x_in_name(6) = "ch4"
      x_in_name(7) = "o2"
      x_in_name(8) = "no"
      x_in_name(9) = "so2"
      x_in_name(10) = "no2"
      x_in_name(11) = "nh3"

      !---------------------------------------------------------------------------------

      ! MEASURE TIME NOW AND TOTAL ELAPSED TIME AT THE END OF THE RUN:
      call itime(now)

      i = INDEX(meteo_file, '.nc')
      index_info = meteo_file(i + 4:)

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      write (runidstring, '(I6.6)') runid
      ncfile = trim(output_dir)//"RTC_OUT_"//runidstring//'.nc'

      if (outputflag >= 2) then
         write (*, '(a)')
         write (*, '(a)') '*****************************************'
         write (*, '(a)') '***** Diagnostic output'
         write (*, '(a)') '*****'
         write (*, *) trim(ncfile)
         write (*, *) ''
      end if

      ntype_target = size(retrieval_output%type_x_target)
      nstate = size(retrieval_output%x_state(:, 1))

      if (allocated(x)) deallocate (x)
      if (allocated(x_err)) deallocate (x_err)
      if (allocated(x_name)) deallocate (x_name)
      if (allocated(x_unit)) deallocate (x_unit)
      allocate (x(ntype_target), &
                x_err(ntype_target), &
                x_name(ntype_target), &
                x_unit(ntype_target))

      !*** DRY AIR MOLE FRACTION
      x_air = sum(retrieval_output%dvair)
      do j = 1, ntype_target

         if (abs(retrieval_output%type_x_target(j)) == 1) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) == 2) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) == 5) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) == 6) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) == 7) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) == 11) then
            scale = 1.E9
            k = 11
         elseif (abs(retrieval_output%type_x_target(j)) >= 100 .and. abs(retrieval_output%type_x_target(j)) <= 199) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) >= 200 .and. abs(retrieval_output%type_x_target(j)) <= 299) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) >= 500 .and. abs(retrieval_output%type_x_target(j)) <= 599) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) >= 600 .and. abs(retrieval_output%type_x_target(j)) <= 699) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) >= 700 .and. abs(retrieval_output%type_x_target(j)) <= 799) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) >= 1100 .and. abs(retrieval_output%type_x_target(j)) <= 1199) then
            scale = 1.E9
            k = 11
         end if

         !*** surface pressure
         psf = retrieval_output%p(nlay + 1)

         if (scale == 1.) then
            x_unit(j) = '1'
         elseif (scale == 1.E2) then
            x_unit(j) = 'percent'
         elseif (scale == 1.E3) then
            x_unit(j) = 'permil'
         elseif (scale == 1.E6) then
            x_unit(j) = 'ppm'
         elseif (scale == 1.E9) then
            x_unit(j) = 'ppb'
         else
            x_unit(j) = 'unknown'
         end if

         i1 = 1 + (nlay*(j - 1))
         i2 = nlay*j

         x(j) = sum(retrieval_output%x_state(i1:i2, retrieval_output%iter))/sum(retrieval_output%dvair)*scale
         x_err(j) = dsqrt(sum(retrieval_output%s_state(i1:i2, i1:i2)))/sum(retrieval_output%dvair)*scale

         l = len(trim(retrieval_output%x_state_name(i1)))
         x_name(j) = retrieval_output%x_state_name(i1) (l - 7:l) ! extract the last seven cahracters (=TYPEXXXX)

         if (x_in_flag(k) .NE. 1) then
            x_in_flag(k) = 1
            x_in(k) = sum(syn_output%x_true(i1:i2))/syn_output%air_col_true*scale
            x_in_unit(k) = x_unit(j)
         end if

      end do

      ntype_target_in = 0
      do k = 1, 11
         if (x_in_flag(k) .EQ. 1) ntype_target_in = ntype_target_in + 1
      end do

      !*** Convert 'x_name' to lowercase

      do i = 1, ntype_target
         do j = 1, len(x_name(i))
            k = ichar(x_name(i) (j:j))
            if (k >= 65 .and. k < 90) x_name(i) (j:j) = char(k + 32)
         end do
      end do

      ! ---------------------------------------------------------- !

      !*** FOR CONVERGED RETRIEVALS ERROR_ID GIVES THE NUMBER OF ITERATIONS
      if (retrieval_output%error_ID == 0) retrieval_output%error_ID = retrieval_output%iter

      !*** Compute chi2
      chi2 = retrieval_output%chi2(retrieval_output%iter)/(sum(retrieval_output%ny) - retrieval_output%dfs)

      !*** Get reference (input) data
      if (allocated(alb_in)) deallocate (alb_in)
      if (allocated(ot_in)) deallocate (ot_in)
      if (allocated(cot_in)) deallocate (cot_in)
      allocate (alb_in(nwin), &
                ot_in(nwin), &
                cot_in(nwin))
      call get_scenario_info_nc_js(meteo_file, alb_in, ot_in, cot_in)

      allocate (grpid(nwin), &
                wave_id(nwin), &
                ot_id(nwin), &
                cot_id(nwin), &
                alb_id(nwin), &
                otin_id(nwin), &
                cotin_id(nwin), &
                albin_id(nwin), &
                tc_id(ntype_target), &
                dfst_id(ntype_target), &
                tcerr_id(ntype_target), &
                tcin_id(ntype_target_in))

      inquire (FILE=trim(ncfile), EXIST=exst)
      if (.not. exst) then

         !*** Create the netCDF file.
         call check(nf90_create(trim(ncfile), nf90_netcdf4, ncid), stat)

         !*** Unlimited dimension for number of observations
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dimid_nobs), stat)

         !*** Indexdata
         call check(nf90_def_var(ncid, "x", NF90_int, dimid_nobs, x_id), stat)
         call check(nf90_def_var(ncid, "y", NF90_int, dimid_nobs, y_id), stat)

         !*** Geometry
         call check(nf90_def_var(ncid, "sza", NF90_float, dimid_nobs, sza_id), stat)
         call check(nf90_def_var(ncid, "vza", NF90_float, dimid_nobs, vza_id), stat)

         call check(nf90_put_att(ncid, sza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, sza_id, "description", "Solar Zenith Angle"), stat)
         call check(nf90_put_att(ncid, vza_id, "unit", "degrees"), stat)
         call check(nf90_put_att(ncid, vza_id, "description", "Viewing Zenith Angle"), stat)

         !*** Airmass
         call check(nf90_def_var(ncid, "airmass", NF90_float, dimid_nobs, xair_id), stat)
         call check(nf90_put_att(ncid, xair_id, "unit", "molec.cm-2"), stat)
         call check(nf90_put_att(ncid, xair_id, "description", "Airmass vertically integrated"), stat)

         !*** Surface Pressure
         call check(nf90_def_var(ncid, "surface_pressure", NF90_float, dimid_nobs, psf_id), stat)
         call check(nf90_put_att(ncid, psf_id, "unit", "hPa"), stat)
         call check(nf90_put_att(ncid, psf_id, "description", "Surface pressure"), stat)

         !*** Geodata
         call check(nf90_def_var(ncid, "latitude", NF90_float, dimid_nobs, lat_id), stat)
         call check(nf90_def_var(ncid, "longitude", NF90_float, dimid_nobs, lon_id), stat)

         call check(nf90_put_att(ncid, lat_id, "unit", "degrees_north"), stat)
         call check(nf90_put_att(ncid, lat_id, "description", "Latitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lon_id, "unit", "degrees_east"), stat)
         call check(nf90_put_att(ncid, lon_id, "description", "Longitude at pixel center"), stat)

         !*** Timedata
         call check(nf90_def_dim(ncid, "ntime", 6, dimid_time), stat)
         call check(nf90_def_var(ncid, "time", NF90_int, dimid_time, time_id), stat)
 call check(nf90_put_var(ncid, time_id, [meta%time(1), meta%time(2), meta%time(3), meta%time(4), meta%time(5), meta%time(6)]), stat)
         call check(nf90_put_att(ncid, time_id, "description", "Date and time as [YYYY,MM,DD,HOUR,MIN,SEC]"), stat)

         !*** Quality flags
         call check(nf90_def_var(ncid, "convergence", NF90_int, dimid_nobs, cf_id), stat)
         call check(nf90_def_var(ncid, "iter", NF90_int, dimid_nobs, it_id), stat)
         call check(nf90_def_var(ncid, "error_id", NF90_int, dimid_nobs, eid_id), stat)

         call check(nf90_put_att(ncid, cf_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, cf_id, "description", "Binary flag for retrieval convergence"), stat)
         call check(nf90_put_att(ncid, it_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, it_id, "description", "Number of iterations needed for convergence"), stat)
         call check(nf90_put_att(ncid, eid_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, eid_id, "description", "Error-ID for retrieval"), stat)

         !*** Quality measures
         call check(nf90_def_var(ncid, "chi2", NF90_float, dimid_nobs, chi_id), stat)
         call check(nf90_def_var(ncid, "dfs", NF90_float, dimid_nobs, dfs_id), stat)
         call check(nf90_def_var(ncid, "dfs_scat", NF90_float, dimid_nobs, dfss_id), stat)

         call check(nf90_put_att(ncid, chi_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, chi_id, "description", "chi-squared quality measure of retrieval"), stat)
         call check(nf90_put_att(ncid, dfs_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, dfs_id, "description", "Degrees of freedom for retrieval"), stat)
         call check(nf90_put_att(ncid, dfss_id, "unit", "-"), stat)
         call check(nf90_put_att(ncid, dfss_id, "description", "Degrees of freedom for scattering retrieval"), stat)

         !*** Retrieved gas concentrations and corresponding retrievals errors and corresponding degrees of freedom
         do n = 1, ntype_target
            call check(nf90_def_var(ncid, trim('x_')//trim(x_name(n)), NF90_float, dimid_nobs, tc_id(n)), stat)
            call check(nf90_def_var(ncid, trim('x_')//trim(x_name(n))//trim('_err'), NF90_float, dimid_nobs, tcerr_id(n)), stat)
            call check(nf90_def_var(ncid, trim('dfst_')//trim(x_name(n)), nf90_float, dimid_nobs, dfst_id(n)), stat)

            call check(nf90_put_att(ncid, tc_id(n), "unit", trim(x_unit(n))), stat)
            call check(nf90_put_att(ncid, tc_id(n), "description", "Retrieved column-averaged dry-air mole fraction of target gas"), stat)
            call check(nf90_put_att(ncid, tcerr_id(n), "unit", trim(x_unit(n))), stat)
            call check(nf90_put_att(ncid, tcerr_id(n), "description", "Random noise error of retrieved column-averaged dry-air mole fraction of target gas"), stat)
            call check(nf90_put_att(ncid, dfst_id(n), "unit", "-"), stat)
            call check(nf90_put_att(ncid, dfst_id(n), "description", "Degrees of freedom for target retrieval"), stat)
         end do

         !*** Input gas concentrations
         n = 1
         do k = 1, 11
            if (x_in_flag(k) .EQ. 1) then
              call check(nf90_def_var(ncid, trim('x_')//trim(x_in_name(k))//trim('_inp'), NF90_float, dimid_nobs, tcin_id(n)), stat)

               call check(nf90_put_att(ncid, tcin_id(n), "unit", trim(x_in_unit(k))), stat)
                call check(nf90_put_att(ncid, tcin_id(n), "description", "Column-averaged dry-air mole fraction of target gas used to simulate spectrum"), stat)

               n = n + 1
            end if
         end do

         do n = 1, nwin
            nwave = measurement(n)%nwave
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** Spectral data
            call check(nf90_def_dim(grpid(n), "nwave", nwave, dimid_wave), stat)

            !*** Define the variables
            call check(nf90_def_var(grpid(n), "wavelength", NF90_float, dimid_wave, wave_id(n)), stat)
            call check(nf90_def_var(grpid(n), "ot", NF90_float, dimid_nobs, ot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "cot", NF90_float, dimid_nobs, cot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "alb", NF90_float, dimid_nobs, alb_id(n)), stat)
            !call check(nf90_def_var(grpid(n), "ot_inp", NF90_float, dimid_nobs, otin_id(n)), stat)
            !call check(nf90_def_var(grpid(n), "cot_inp", NF90_float, dimid_nobs, cotin_id(n)), stat)
            !call check(nf90_def_var(grpid(n), "alb_inp", NF90_float, dimid_nobs, albin_id(n)), stat)

            call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
            call check(nf90_put_att(grpid(n), wave_id(n), "description", "Wavelength grid of spectrum"), stat)

            call check(nf90_put_att(grpid(n), ot_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), ot_id(n), "description", "Retrieved total optical thickness"), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "description", "Retrieved cirrus optical thickness"), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "description", "Retrieved surface albedo"), stat)

            !call check(nf90_put_att(grpid(n), otin_id(n), "unit", "-"), stat)
            !call check(nf90_put_att(grpid(n), otin_id(n), "description", "Total optical thickness used to simulate spectrum"), stat)
            !call check(nf90_put_att(grpid(n), cotin_id(n), "unit", "-"), stat)
            !call check(nf90_put_att(grpid(n), cotin_id(n), "description", "Cirrus optical thickness used to simulate spectrum"), stat)
            !call check(nf90_put_att(grpid(n), albin_id(n), "unit", "-"), stat)
            !call check(nf90_put_att(grpid(n), albin_id(n), "description", "Surface albedo used to simulate spectrum"), stat)

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)

         end do

         lastindex = 1

         !*** GET VARIABLE IDs
      else

         !*** Open the netCDF file and append
         call check(NF90_OPEN(trim(ncfile), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** Indexdata
         call check(nf90_inq_varid(ncid, "x", x_id), stat)
         call check(nf90_inq_varid(ncid, "y", y_id), stat)

         !*** Geometry
         call check(nf90_inq_varid(ncid, "sza", sza_id), stat)
         call check(nf90_inq_varid(ncid, "vza", vza_id), stat)

         !*** Airmass
         call check(nf90_inq_varid(ncid, "airmass", xair_id), stat)

         !*** Surface Pressure
         call check(nf90_inq_varid(ncid, "surface_pressure", psf_id), stat)

         !*** Geodata
         call check(nf90_inq_varid(ncid, "latitude", lat_id), stat)
         call check(nf90_inq_varid(ncid, "longitude", lon_id), stat)

         !*** Quality flags
         call check(nf90_inq_varid(ncid, "convergence", cf_id), stat)
         call check(nf90_inq_varid(ncid, "iter", it_id), stat)
         call check(nf90_inq_varid(ncid, "error_id", eid_id), stat)

         !*** Quality measures
         call check(nf90_inq_varid(ncid, "chi2", chi_id), stat)
         call check(nf90_inq_varid(ncid, "dfs", dfs_id), stat)
         call check(nf90_inq_varid(ncid, "dfs_scat", dfss_id), stat)

         !*** Retrieved gas concentrations and corresponding retrievals errors
         do n = 1, ntype_target
            call check(nf90_inq_varid(ncid, trim('x_')//trim(x_name(n)), tc_id(n)), stat)
            call check(nf90_inq_varid(ncid, trim('x_')//trim(x_name(n))//trim('_err'), tcerr_id(n)), stat)
            call check(nf90_inq_varid(ncid, trim("dfst_")//trim(x_name(n)), dfst_id(n)), stat)
         end do

         !*** Input gas concentrations
         n = 1
         do k = 1, 11
            if (x_in_flag(k) .EQ. 1) then
               call check(nf90_inq_varid(ncid, trim('x_')//trim(x_in_name(k))//trim('_inp'), tcin_id(n)), stat)
               n = n + 1
            end if
         end do

         !*** Spectrally resolved data
         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            call check(nf90_inq_varid(grpid(n), "ot", ot_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "cot", cot_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "alb", alb_id(n)), stat)
            !call check(nf90_inq_varid(grpid(n), "ot_inp", otin_id(n)), stat)
            !call check(nf90_inq_varid(grpid(n), "cot_inp", cotin_id(n)), stat)
            !call check(nf90_inq_varid(grpid(n), "alb_inp", albin_id(n)), stat)
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

      !*** Airmass
      call check(nf90_put_var(ncid, xair_id, x_air, start=(/lastindex/)), stat)

      !*** Surface Pressure
      call check(nf90_put_var(ncid, psf_id, psf, start=(/lastindex/)), stat)

      !*** Geodata
      call check(nf90_put_var(ncid, lat_id, meta%lat(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, lon_id, meta%lon(1), start=(/lastindex/)), stat)

      !*** Quality flags
      call check(nf90_put_var(ncid, cf_id, retrieval_output%convergence, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, it_id, retrieval_output%iter, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, eid_id, retrieval_output%error_ID, start=(/lastindex/)), stat)

      !*** Quality measures
      call check(nf90_put_var(ncid, chi_id, chi2, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, dfs_id, retrieval_output%dfs, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, dfss_id, retrieval_output%dfs_scat, start=(/lastindex/)), stat)

      !*** Retrieved gas concentrations and corresponding retrievals errors and corresponding degrees of freedom
      do n = 1, ntype_target
         call check(nf90_put_var(ncid, tc_id(n), x(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(ncid, tcerr_id(n), x_err(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(ncid, dfst_id(n), retrieval_output%dfs_target(n), start=(/lastindex/)), stat)
      end do

      !*** Input gas concentrations
      n = 1
      do k = 1, 11
         if (x_in_flag(k) .EQ. 1) then
            call check(nf90_put_var(ncid, tcin_id(n), x_in(k), start=(/lastindex/)), stat)
         end if
      end do

      !*** Spectral data
      do n = 1, nwin
         call check(nf90_put_var(grpid(n), ot_id(n), retrieval_output%ot(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), cot_id(n), retrieval_output%cot(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), alb_id(n), retrieval_output%albedo(n), start=(/lastindex/)), stat)
         !call check(nf90_put_var(grpid(n), otin_id(n), ot_in, start=(/lastindex/)), stat)
         !call check(nf90_put_var(grpid(n), cotin_id(n), cot_in, start=(/lastindex/)), stat)
         !call check(nf90_put_var(grpid(n), albin_id(n), alb_in, start=(/lastindex/)), stat)
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (outputflag >= 2) then
         write (*, '(a)') '***** Done.'
         write (*, 80) now
80       format('***** Time at the end: ', i2.2, ':', i2.2, ':', i2.2)
         write (*, '(a)') '*****************************************'
      end if

   end subroutine diagnostics_retrieve_nc_js

   subroutine diagnostics_retrieve_nc_ls( &
      output_dir, runid, &
      meteo_file, &
      win_ini, &
      synthetic_input_flag, &
      syn_output, &
      retrieval_output, &
      meta, &
      measurement, &
      nlay, outputflag)
      !*** arguments
      character(len=*), intent(in) :: output_dir
      integer, intent(in) :: runid
      character(len=*), intent(in) :: meteo_file
      type(window_ini), dimension(:), intent(in) :: win_ini
      integer, intent(in) :: synthetic_input_flag
      type(synthetic_data), intent(in) :: syn_output
      type(retrieval_data), intent(inout) :: retrieval_output
      type(metadata), intent(in) :: meta
      type(spectrum), dimension(:), intent(in) :: measurement
      integer, intent(in) :: nlay, outputflag
      !*** local variables
      integer :: nwin
      integer :: k, l, n, i, j, i1, i2, length, now(3), io, nstate, ntype_target, ntype_target_in
      character(299) :: ncfile
      character(6):: runidstring
      real(double) :: scale, chi2
      real(double), dimension(:), allocatable :: x, x_err, wave_start_in, wave_stop_in, alb_in, ot_in, cot_in
      real(double), dimension(11) :: x_in
      integer, dimension(11) :: x_in_flag
      character(99), dimension(11) :: x_in_name, x_in_unit
      logical :: exst
      character(99) :: unit
      character(99), dimension(:), allocatable :: x_name, x_unit
      integer :: sx, sy
      integer :: ncid, ierr, stat, ngroup, dimid_nobs, dimid_wave, lastindex, dimid_time, nwave
      integer :: sza_id, vza_id, xair_id, psf_id, lat_id, lon_id, time_id, x_id, y_id
      integer :: cf_id, chi_id, dfs_id, dfss_id, eid_id, it_id
      integer, dimension(:), allocatable :: grpid, wave_id, ot_id, cot_id, alb_id, otin_id, cotin_id, albin_id, tc_id, dfst_id, tcerr_id, tcin_id
      character*1 :: ch
      character(stringlen) :: group_name, index_info
      real(double) :: x_air
      real(double) :: psf
      !-----------------------------------------------------------

      nwin = size(win_ini)

      if (nwin .NE. size(measurement)) then
         call writelog('DIAGNOSTICS_SIM_NC_LS: expected same number of windows for measurement and retrieval', 8)
         return
      end if

      !---------------------------------------------------------------------------------
      x_in_flag(:) = 0
      x_in_name(1) = "h2o"
      x_in_name(2) = "co2"
      x_in_name(3) = "o3"
      x_in_name(4) = "n2o"
      x_in_name(5) = "co"
      x_in_name(6) = "ch4"
      x_in_name(7) = "o2"
      x_in_name(8) = "no"
      x_in_name(9) = "so2"
      x_in_name(10) = "no2"
      x_in_name(11) = "nh3"

      !---------------------------------------------------------------------------------

      ! MEASURE TIME NOW AND TOTAL ELAPSED TIME AT THE END OF THE RUN:
      call itime(now)

      i = INDEX(meteo_file, '.nc')
      index_info = meteo_file(i + 4:)

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 6), '(I6)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 6), '(I6)') sy

      write (runidstring, '(I6.6)') runid
      ncfile = trim(output_dir)//"RTC_OUT_"//runidstring//'.nc'

      if (outputflag >= 2) then
         write (*, '(a)')
         write (*, '(a)') '*****************************************'
         write (*, '(a)') '***** Diagnostic output'
         write (*, '(a)') '*****'
         write (*, *) trim(ncfile)
         write (*, *) ''
      end if

      ntype_target = size(retrieval_output%type_x_target)
      nstate = size(retrieval_output%x_state(:, 1))

      if (allocated(x)) deallocate (x)
      if (allocated(x_err)) deallocate (x_err)
      if (allocated(x_name)) deallocate (x_name)
      if (allocated(x_unit)) deallocate (x_unit)
      allocate (x(ntype_target), &
                x_err(ntype_target), &
                x_name(ntype_target), &
                x_unit(ntype_target))

      !*** DRY AIR MOLE FRACTION
      x_air = sum(retrieval_output%dvair)
      do j = 1, ntype_target

         if (abs(retrieval_output%type_x_target(j)) == 1) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) == 2) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) == 5) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) == 6) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) == 7) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) == 11) then
            scale = 1.E9
            k = 11
         elseif (abs(retrieval_output%type_x_target(j)) >= 100 .and. abs(retrieval_output%type_x_target(j)) <= 199) then
            scale = 1.
            k = 1
         elseif (abs(retrieval_output%type_x_target(j)) >= 200 .and. abs(retrieval_output%type_x_target(j)) <= 299) then
            scale = 1.E6
            k = 2
         elseif (abs(retrieval_output%type_x_target(j)) >= 500 .and. abs(retrieval_output%type_x_target(j)) <= 599) then
            scale = 1.E6
            k = 5
         elseif (abs(retrieval_output%type_x_target(j)) >= 600 .and. abs(retrieval_output%type_x_target(j)) <= 699) then
            scale = 1.E6
            k = 6
         elseif (abs(retrieval_output%type_x_target(j)) >= 700 .and. abs(retrieval_output%type_x_target(j)) <= 799) then
            scale = 1.
            k = 7
         elseif (abs(retrieval_output%type_x_target(j)) >= 1100 .and. abs(retrieval_output%type_x_target(j)) <= 1199) then
            scale = 1.E9
            k = 11
         end if

         !*** surface pressure
         psf = retrieval_output%p(nlay + 1)

         if (scale == 1.) then
            x_unit(j) = '1'
         elseif (scale == 1.E2) then
            x_unit(j) = 'percent'
         elseif (scale == 1.E3) then
            x_unit(j) = 'permil'
         elseif (scale == 1.E6) then
            x_unit(j) = 'ppm'
         elseif (scale == 1.E9) then
            x_unit(j) = 'ppb'
         else
            x_unit(j) = 'unknown'
         end if

         i1 = 1 + (nlay*(j - 1))
         i2 = nlay*j

         x(j) = sum(retrieval_output%x_state(i1:i2, retrieval_output%iter))/sum(retrieval_output%dvair)*scale
         x_err(j) = dsqrt(sum(retrieval_output%s_state(i1:i2, i1:i2)))/sum(retrieval_output%dvair)*scale

         l = len(trim(retrieval_output%x_state_name(i1)))
         x_name(j) = retrieval_output%x_state_name(i1) (l - 7:l) ! extract the last seven characters (=TYPEXXXX)

         if (synthetic_input_flag .eq. 1) then
            if (x_in_flag(k) .NE. 1) then
               x_in_flag(k) = 1
               x_in(k) = sum(syn_output%x_true(i1:i2))/syn_output%air_col_true*scale
               x_in_unit(k) = x_unit(j)
            end if
         end if

      end do
      
      if (synthetic_input_flag .eq. 1) then
         ntype_target_in = 0
         do k = 1, 11
            if (x_in_flag(k) .EQ. 1) ntype_target_in = ntype_target_in + 1
         end do
      end if

      !*** Convert 'x_name' to lowercase

      do i = 1, ntype_target
         do j = 1, len(x_name(i))
            k = ichar(x_name(i) (j:j))
            if (k >= 65 .and. k < 90) x_name(i) (j:j) = char(k + 32)
         end do
      end do

      ! ---------------------------------------------------------- !

      !*** FOR CONVERGED RETRIEVALS ERROR_ID GIVES THE NUMBER OF ITERATIONS
      if (retrieval_output%error_ID == 0) retrieval_output%error_ID = retrieval_output%iter

      !*** Compute chi2
      chi2 = retrieval_output%chi2(retrieval_output%iter)/(sum(retrieval_output%ny) - retrieval_output%dfs)

      !*** Get reference (input) data
      if (allocated(alb_in)) deallocate (alb_in)
      if (allocated(ot_in)) deallocate (ot_in)
      if (allocated(cot_in)) deallocate (cot_in)
      allocate (alb_in(nwin), &
                ot_in(nwin), &
                cot_in(nwin))
      call get_scenario_info_nc_js(meteo_file, alb_in, ot_in, cot_in)

      allocate (grpid(nwin), &
                wave_id(nwin), &
                ot_id(nwin), &
                cot_id(nwin), &
                alb_id(nwin), &
                otin_id(nwin), &
                cotin_id(nwin), &
                albin_id(nwin), &
                tc_id(ntype_target), &
                dfst_id(ntype_target), &
                tcerr_id(ntype_target), &
                tcin_id(ntype_target_in))

      inquire (FILE=trim(ncfile), EXIST=exst)
      if (.not. exst) then

         !*** Create the netCDF file.
         call check(nf90_create(trim(ncfile), nf90_netcdf4, ncid), stat)

         !*** Unlimited dimension for number of observations
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dimid_nobs), stat)

         !*** Indexdata
         call check(nf90_def_var(ncid, "x", nf90_int, dimid_nobs, x_id), stat)
         call check(nf90_def_var(ncid, "y", nf90_int, dimid_nobs, y_id), stat)

         !*** Timedata
         call check(nf90_def_var(ncid, "time", nf90_int, dimid_time, time_id), stat)
         call check(nf90_put_var(ncid, time_id, meta%seconds_since_reference), stat)
         call check(nf90_put_att(ncid, time_id, "long_name", "seconds since reference"), stat)
         call check(nf90_put_att(ncid, time_id, "units", "YYYY-MM-DDThh:mm:ssZ"), stat)

         !*** Geodata
         call check(nf90_def_var(ncid, "latitude", nf90_float, dimid_nobs, lat_id), stat)
         call check(nf90_put_att(ncid, lat_id, "long_name", "Latitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lat_id, "units", "degrees north"), stat)

         call check(nf90_def_var(ncid, "longitude", nf90_float, dimid_nobs, lon_id), stat)
         call check(nf90_put_att(ncid, lon_id, "long_name", "Longitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lon_id, "units", "degrees east"), stat)

         !*** Geometry
         call check(nf90_def_var(ncid, "solar_zenith_angle", nf90_float, dimid_nobs, sza_id), stat)
         call check(nf90_put_att(ncid, sza_id, "long_name", "Solar zenith angle"), stat)
         call check(nf90_put_att(ncid, sza_id, "units", "degrees"), stat)

         call check(nf90_def_var(ncid, "viewing_zenith_angle", nf90_float, dimid_nobs, vza_id), stat)
         call check(nf90_put_att(ncid, vza_id, "long_name", "Viewing zenith angle"), stat)
         call check(nf90_put_att(ncid, vza_id, "units", "degrees"), stat)

         !*** Airmass
         call check(nf90_def_var(ncid, "airmass", nf90_float, dimid_nobs, xair_id), stat)
         call check(nf90_put_att(ncid, xair_id, "long_name", "Vertically integrated airmass"), stat)
         call check(nf90_put_att(ncid, xair_id, "units", "molecules cm-2"), stat)

         !*** Surface Pressure
         call check(nf90_def_var(ncid, "surface_pressure", nf90_float, dimid_nobs, psf_id), stat)
         call check(nf90_put_att(ncid, psf_id, "long_name", "Surface pressure"), stat)
         call check(nf90_put_att(ncid, psf_id, "units", "hPa"), stat)

         !*** Quality flags
         call check(nf90_def_var(ncid, "convergence", nf90_int, dimid_nobs, cf_id), stat)
         call check(nf90_put_att(ncid, cf_id, "long_name", "Binary flag for retrieval convergence"), stat)

         call check(nf90_def_var(ncid, "iter", nf90_int, dimid_nobs, it_id), stat)
         call check(nf90_put_att(ncid, it_id, "long_name", "Number of iterations needed for convergence"), stat)

         call check(nf90_def_var(ncid, "error_id", nf90_int, dimid_nobs, eid_id), stat)
         call check(nf90_put_att(ncid, eid_id, "long_name", "Error-ID for retrieval"), stat)

         !*** Quality measures
         call check(nf90_def_var(ncid, "chi2", nf90_float, dimid_nobs, chi_id), stat)
         call check(nf90_put_att(ncid, chi_id, "long_name", "Chi-squared quality measure of retrieval"), stat)

         call check(nf90_def_var(ncid, "dfs", nf90_float, dimid_nobs, dfs_id), stat)
         call check(nf90_put_att(ncid, dfs_id, "long_name", "Degrees of freedom for retrieval"), stat)

         call check(nf90_def_var(ncid, "dfs_scat", nf90_float, dimid_nobs, dfss_id), stat)
         call check(nf90_put_att(ncid, dfss_id, "long_name", "Degrees of freedom for scattering retrieval"), stat)

         !*** Retrieved gas concentrations and corresponding retrievals errors and corresponding degrees of freedom
         do n = 1, ntype_target
            call check(nf90_def_var(ncid, trim('x_')//trim(x_name(n)), nf90_float, dimid_nobs, tc_id(n)), stat)
            call check(nf90_put_att(ncid, tc_id(n), "long_name", "Retrieved column-averaged dry-air mole fraction of target gas"), stat)
            call check(nf90_put_att(ncid, tc_id(n), "units", trim(x_unit(n))), stat)

            call check(nf90_def_var(ncid, trim('x_')//trim(x_name(n))//trim('_err'), nf90_float, dimid_nobs, tcerr_id(n)), stat)
            call check(nf90_put_att(ncid, tcerr_id(n), "long_name", "Random noise error of retrieved column-averaged dry-air mole fraction of target gas"), stat)
            call check(nf90_put_att(ncid, tcerr_id(n), "units", trim(x_unit(n))), stat)

            call check(nf90_def_var(ncid, trim('dfst_')//trim(x_name(n)), nf90_float, dimid_nobs, dfst_id(n)), stat)
            call check(nf90_put_att(ncid, dfst_id(n), "long_name", "Degrees of freedom for target retrieval"), stat)
         end do

         if (synthetic_input_flag .eq. 1) then
            !*** Input gas concentrations
            n = 1
            do k = 1, 11
               if (x_in_flag(k) .EQ. 1) then
                  call check(nf90_def_var(ncid, trim('x_')//trim(x_in_name(k))//trim('_inp'), nf90_float, dimid_nobs, tcin_id(n)), stat)
                  call check(nf90_put_att(ncid, tcin_id(n), "long_name", "Column-averaged dry-air mole fraction of target gas used to simulate spectrum"), stat)
                  call check(nf90_put_att(ncid, tcin_id(n), "units", trim(x_in_unit(k))), stat)
                  n = n + 1
               end if
            end do
         end if

         do n = 1, nwin
            nwave = measurement(n)%nwave
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** Spectral data
            call check(nf90_def_dim(grpid(n), "nwave", nwave, dimid_wave), stat)

            !*** Define the variables
            call check(nf90_def_var(grpid(n), "wavelength", nf90_float, dimid_wave, wave_id(n)), stat)
            call check(nf90_put_att(grpid(n), wave_id(n), "long_name", "Wavelength grid of spectrum"), stat)
            call check(nf90_put_att(grpid(n), wave_id(n), "units", "nm"), stat)

            call check(nf90_def_var(grpid(n), "ot", nf90_float, dimid_nobs, ot_id(n)), stat)
            call check(nf90_put_att(grpid(n), ot_id(n), "long_name", "Retrieved total optical thickness"), stat)

            call check(nf90_def_var(grpid(n), "cot", nf90_float, dimid_nobs, cot_id(n)), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "long_name", "Retrieved cirrus optical thickness"), stat)

            call check(nf90_def_var(grpid(n), "alb", nf90_float, dimid_nobs, alb_id(n)), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "long_name", "Retrieved surface albedo"), stat)

            if (synthetic_input_flag .eq. 1) then
               call check(nf90_def_var(grpid(n), "ot_inp", nf90_float, dimid_nobs, otin_id(n)), stat)
               call check(nf90_put_att(grpid(n), otin_id(n), "long_name", "Total optical thickness used to simulate spectrum"), stat)

               call check(nf90_def_var(grpid(n), "cot_inp", nf90_float, dimid_nobs, cotin_id(n)), stat)
               call check(nf90_put_att(grpid(n), cotin_id(n), "long_name", "Cirrus optical thickness used to simulate spectrum"), stat)

               call check(nf90_def_var(grpid(n), "alb_inp", nf90_float, dimid_nobs, albin_id(n)), stat)
               call check(nf90_put_att(grpid(n), albin_id(n), "long_name", "Surface albedo used to simulate spectrum"), stat)
            end if

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)

         end do

         lastindex = 1

         !*** GET VARIABLE IDs
      else

         !*** Open the netCDF file and append
         call check(nf90_open(trim(ncfile), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** Indexdata
         call check(nf90_inq_varid(ncid, "x", x_id), stat)
         call check(nf90_inq_varid(ncid, "y", y_id), stat)

         !*** Geometry
         call check(nf90_inq_varid(ncid, "solar_zenith_angle", sza_id), stat)
         call check(nf90_inq_varid(ncid, "viewing_zenith_angle", vza_id), stat)

         !*** Airmass
         call check(nf90_inq_varid(ncid, "airmass", xair_id), stat)

         !*** Surface Pressure
         call check(nf90_inq_varid(ncid, "surface_pressure", psf_id), stat)

         !*** Geodata
         call check(nf90_inq_varid(ncid, "latitude", lat_id), stat)
         call check(nf90_inq_varid(ncid, "longitude", lon_id), stat)

         !*** Quality flags
         call check(nf90_inq_varid(ncid, "convergence", cf_id), stat)
         call check(nf90_inq_varid(ncid, "iter", it_id), stat)
         call check(nf90_inq_varid(ncid, "error_id", eid_id), stat)

         !*** Quality measures
         call check(nf90_inq_varid(ncid, "chi2", chi_id), stat)
         call check(nf90_inq_varid(ncid, "dfs", dfs_id), stat)
         call check(nf90_inq_varid(ncid, "dfs_scat", dfss_id), stat)

         !*** Retrieved gas concentrations and corresponding retrievals errors
         do n = 1, ntype_target
            call check(nf90_inq_varid(ncid, trim('x_')//trim(x_name(n)), tc_id(n)), stat)
            call check(nf90_inq_varid(ncid, trim('x_')//trim(x_name(n))//trim('_err'), tcerr_id(n)), stat)
            call check(nf90_inq_varid(ncid, trim("dfst_")//trim(x_name(n)), dfst_id(n)), stat)
         end do

         if (synthetic_input_flag .eq. 1) then
            !*** Input gas concentrations
            n = 1
            do k = 1, 11
               if (x_in_flag(k) .EQ. 1) then
                  call check(nf90_inq_varid(ncid, trim('x_')//trim(x_in_name(k))//trim('_inp'), tcin_id(n)), stat)
                  n = n + 1
               end if
            end do
         end if

         !*** Spectrally resolved data
         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            call check(nf90_inq_varid(grpid(n), "ot", ot_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "cot", cot_id(n)), stat)
            call check(nf90_inq_varid(grpid(n), "alb", alb_id(n)), stat)
            if (synthetic_input_flag .eq. 1) then
               call check(nf90_inq_varid(grpid(n), "ot_inp", otin_id(n)), stat)
               call check(nf90_inq_varid(grpid(n), "cot_inp", cotin_id(n)), stat)
               call check(nf90_inq_varid(grpid(n), "alb_inp", albin_id(n)), stat)
            end if
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

      !*** Geodata
      call check(nf90_put_var(ncid, lat_id, meta%lat(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, lon_id, meta%lon(1), start=(/lastindex/)), stat)

      !*** Geometry
      call check(nf90_put_var(ncid, sza_id, meta%sza, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, vza_id, meta%iza, start=(/lastindex/)), stat)

      !*** Airmass
      call check(nf90_put_var(ncid, xair_id, x_air, start=(/lastindex/)), stat)

      !*** Surface Pressure
      call check(nf90_put_var(ncid, psf_id, psf, start=(/lastindex/)), stat)

      !*** Quality flags
      call check(nf90_put_var(ncid, cf_id, retrieval_output%convergence, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, it_id, retrieval_output%iter, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, eid_id, retrieval_output%error_ID, start=(/lastindex/)), stat)

      !*** Quality measures
      call check(nf90_put_var(ncid, chi_id, chi2, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, dfs_id, retrieval_output%dfs, start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, dfss_id, retrieval_output%dfs_scat, start=(/lastindex/)), stat)

      !*** Retrieved gas concentrations and corresponding retrievals errors and corresponding degrees of freedom
      do n = 1, ntype_target
         call check(nf90_put_var(ncid, tc_id(n), x(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(ncid, tcerr_id(n), x_err(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(ncid, dfst_id(n), retrieval_output%dfs_target(n), start=(/lastindex/)), stat)
      end do

      if (synthetic_input_flag .eq. 1) then
         !*** Input gas concentrations
         n = 1
         do k = 1, 11
            if (x_in_flag(k) .EQ. 1) then
               call check(nf90_put_var(ncid, tcin_id(n), x_in(k), start=(/lastindex/)), stat)
            end if
         end do
      end if

      !*** Spectral data
      do n = 1, nwin
         call check(nf90_put_var(grpid(n), ot_id(n), retrieval_output%ot(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), cot_id(n), retrieval_output%cot(n), start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), alb_id(n), retrieval_output%albedo(n), start=(/lastindex/)), stat)
         if (synthetic_input_flag .eq. 1) then
            call check(nf90_put_var(grpid(n), otin_id(n), ot_in, start=(/lastindex/)), stat)
            call check(nf90_put_var(grpid(n), cotin_id(n), cot_in, start=(/lastindex/)), stat)
            call check(nf90_put_var(grpid(n), albin_id(n), alb_in, start=(/lastindex/)), stat)
         end if
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (outputflag >= 2) then
         write (*, '(a)') '***** Done.'
         write (*, 80) now
80       format('***** Time at the end: ', i2.2, ':', i2.2, ':', i2.2)
         write (*, '(a)') '*****************************************'
      end if

   end subroutine diagnostics_retrieve_nc_ls


end module diagnostics_module
