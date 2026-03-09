program main
   use header_module
   use read_settings_module
   use read_synsettings_module
   use read_solar_module, only: read_solar, interpolate_solar_spectrum
   use atmo_interface_create_module, only: atmosphere_input, output_atm, output_meteo, output_atm_nc_js
   use atmosphere_internal_module, only: atmospheric_scenario, atmosphere_interpolate
   use calculate_syn_spectrum_module, only: calculate_syn_spectrum, spectrum, &
                                            Mie_lut, read_aerosol_netcdf, cirrus_table, &
                              read_settings, window_ini, settings_flags, file_paths, read_win_xsdb, aero, window_spectrum, absorbers
   use optic_cirrus_module, only: read_cirrus_netcdf
   use spec_interface_create_module, only: synthetic_interface_init, synthetic_interface_close, synthetic_interface, &
                                           output_l1b_nc, output_l1b, output_l1b_nc_js, output_l1b_nc_ls, output_lut_nc_js
   use retrieval_module, only: get_absorbers
   use read_miprep_module, only: open_miprep, close_miprep
   implicit none
   integer :: iargc
   integer :: i, j, n, nfile, stat, io
   character(stringlen) :: atm_name, atm_file, spectrum_file, meteo_file, lut_file, fname_nc              ! file names
   character(stringlen), dimension(:), allocatable :: atm
   character(stringlen) :: runpath
   type(sun_spectrum) :: sun_input  ! input reference solar spectrum
   !*** Variables for retrieval algorithm
   type(spectrum), dimension(:), allocatable :: simulation
   type(spectrum), dimension(:), allocatable :: simulation_hi
   type(atmospheric_scenario) :: atm_scenario
   type(metadata) :: meta
   type(altitude_grid) :: grid
   type(settings_flags) :: flag
   type(file_paths) :: path
   type(window_ini), dimension(:), allocatable :: win_ini                       ! window input in retrieval.ini
   !   type(sun_spectrum) :: sun_input                                              ! input solar spectrum
   type(Mie_lut) :: aero_lut                                                    ! Aerosol scattering LUT
   type(cirrus_table) :: cirrus_lut                                             ! Cirrus scattering LUT
   type(window_spectrum), dimension(:), allocatable :: win
   type(absorbers) :: abs
   type(atmosphere) :: atm_rt, atm_xs, atm_retr
   type(aero), dimension(:), allocatable :: aerosol_ini, aerosol
   type(syn_instrument_noise_settings) :: synsettings
   real(double), dimension(:), allocatable :: dvair
   real(double) :: z_tropopause, z_bl, starttime, stoptime, xco2, xco2_in, co2_scaling
   integer :: exitflag, ierr, nrow(2)
   character(stringlen) :: settings_file, aerosolfile, synsettings_file
   character(6):: runidstring
   integer, dimension(:), allocatable :: pixelid                          ! Identifier for groundpixel
   logical :: exst
   !*** Command line input
   character(stringlen) :: arg, filename, index_info                                        ! inputfile with calculation settings
   integer :: runid                                                       ! Identifier for run
   integer :: j1, j2
   !------------------------------------------------------------------------

   !*** Read in command line arguments RUNID, atm_name, and runpath
   if (iargc() .eq. 2) then
      call getarg(1, arg)
      read (arg, *) runid
      if (runid .lt. 0 .or. runid .gt. 999999) call stopretrieval('RUNID must be between 0 and 999999')
      call getarg(2, atm_file)
      runpath = './'
   elseif (iargc() .eq. 3) then
      call getarg(1, arg)
      read (arg, *) runid
      call getarg(2, atm_file)
      call getarg(3, runpath)
   else
      call stopretrieval('Give correct number of arguments')
   end if

   write (runidstring, '(I6.6)') runid
   settings_file = trim(runpath)//'INI/settings_RTC_create_'//runidstring//'.nml'
   call read_settings(settings_file, &
                      win_ini, &
                      aerosol_ini, &
                      grid, &
                      flag, &
                      path, &
                      ierr)

   synsettings_file = trim(runpath)//'INI/syn_create_'//runidstring//'.nml'
   call read_synsettings(synsettings_file, flag%atm, flag%output, size(win_ini), meta, synsettings, xco2, ierr)
   if (ierr .ne. 0) call stopretrieval("MAIN: error reading synsettings")

   !*** Get nfile = total number of input spectra to be retrieved
   open (newunit(io), file=trim(atm_file), action='read')
   nfile = 1
   read (io, *, iostat=stat)
   do while (stat .ne. -1)
      nfile = nfile + 1
      read (io, *, iostat=stat)
   end do
   close (io)
   nfile = nfile - 1
   !*** Get filenames of input spectra to be retrieved
   allocate (atm(nfile), pixelid(nfile))
   open (newunit(io), file=trim(atm_file), action='read')
   if (flag%atm == 4) then
      do i = 1, nfile
         read (io, *) atm(i), pixelid(i)
      end do
      call open_miprep(trim(path%meteo)//atm(1), ierr)
      if (ierr .ne. 0) call stopretrieval("MAIN: error opening MIPrep files")
   else
      do i = 1, nfile
         read (io, '(A)') atm(i)
      end do
   end if
   close (io)

   !*** Read absorption cross sections from database into memory
   call read_win_xsdb(flag, win_ini, ierr)

   !*** Read in Tables with Mie / T-matrix scattering properties
   if (flag%scat == 2 .or. flag%scat == 3 .or. flag%scat == 5 .or. flag%scat == 6) then
      call read_aerosol_netcdf(path%mie, flag%output, aero_lut, ierr)
   else
      do j = 1, size(aerosol_ini)
         if (aerosol_ini(j)%CirrusFlag == 0) then
            call read_aerosol_netcdf(path%mie, flag%output, aero_lut, ierr)
            exit
         end if
      end do
   end if
   if (ierr .ne. 0) call stopretrieval("MAIN: error reading aerosol table")

   !*** Read in cirrus table
   if (flag%scat == 2 .or. flag%scat == 4 .or. flag%scat == 5 .or. flag%scat == 6) then
      call read_cirrus_netcdf(path%cirrus, 30.d0, flag%output, cirrus_lut, ierr)
   else
      do j = 1, size(aerosol_ini)
         if (aerosol_ini(j)%CirrusFlag == 1) then
            call read_cirrus_netcdf(path%cirrus, aerosol_ini(j)%tilt_angle, flag%output, cirrus_lut, ierr)
            exit
         end if
      end do
   end if
   if (ierr .ne. 0) call stopretrieval("MAIN: error reading cirrus table")

   !*** Read solar irradiance spectrum
   call read_solar(path%sun, sun_input, ierr)
   !*** Interpolate irradiance to internal wavelength grid
   call interpolate_solar_spectrum(sun_input, win_ini, ierr)

   !*** Set up measurement wavelength grid, ILS and convoluted irradiance
   if (flag%ilscalc .eq. 0) then  ! from S5P file
      call synthetic_interface_init(atm(1), path%ils)
   end if
   call synthetic_interface(flag%ilscalc, flag%inv, path%ils, win_ini, win)

   !*** For measuring CPU time
   call cpu_time(starttime)

   do i = 1, nfile

      !** Get atmosphere file name
      atm_name = atm(i)

      !*** Read atmospheric input. Here we need the observation name obs_name for the first time:
      if (flag%scat == 6) then
         call atmosphere_input( &
            trim(path%meteo)//atm_name, &
            pixelid(i), &
            flag%atm, atm_scenario, &
            meta, &
            ierr, &
            aerosolfile, &
            aerosol_ini)
         if (ierr .ne. 0) cycle
      elseif (flag%scat == 5) then
         call atmosphere_input( &
            trim(path%meteo)//atm_name, &
            pixelid(i), &
            flag%atm, atm_scenario, &
            meta, &
            ierr, &
            aerosolfile)
         if (ierr .ne. 0) cycle
      else
         call atmosphere_input( &
            trim(path%meteo)//atm_name, &
            pixelid(i), &
            flag%atm, atm_scenario, &
            meta, &
            ierr)
         aerosolfile = trim(path%meteo)//atm_name
         if (ierr .ne. 0) cycle
      end if

      !*** If generating LUT with high-resolution spectra (flag%atm==7): scale CO2
      !*** column to yield the XCO2 specified in syn_create_??????.nml
      if (flag%atm == 7) then

         !*** First interpolate meteo data to radiative transfer grid and retrieval grid such that input XCO2 can be computed
         call atmosphere_interpolate( &
            atm_scenario, &
            grid, &
            meta%surface_elevation, &
            meta%lat(1), &
            win_ini, &
            flag%output, &
            win, &
            atm_xs, &
            atm_rt, &
            atm_retr, &
            dvair, &
            z_tropopause, &
            z_bl, &
            ierr)
         if (ierr .ne. 0) cycle

         !*** Compute actual XCO2 from the input meteo-file
         do n = 1, size(win_ini)
            do j = 1, win_ini(n)%ntype
               if (win_ini(n)%type_x(j) == 2) then
                  xco2_in = sum(win(1)%dv_x(:, 1))/sum(dvair)
               end if
            end do
         end do

         co2_scaling = xco2/xco2_in*1.E-6   ! Scaling factor needed to scale the meteo co2 column in order to get specified XCO2
         atm_scenario%co2 = atm_scenario%co2*co2_scaling    ! Scale meteo co2 column accordingly
      end if

      !*** Interpolate meteo data to radiative transfer grid and retrieval grid
      call atmosphere_interpolate( &
         atm_scenario, &
         grid, &
         meta%surface_elevation, &
         meta%lat(1), &
         win_ini, &
         flag%output, &
         win, &
         atm_xs, &
         atm_rt, &
         atm_retr, &
         dvair, &
         z_tropopause, &
         z_bl, &
         ierr)
      if (ierr .ne. 0) cycle

      if (flag%output > 0) then
         !*** Print dry air volume mixing ratios
         write (*, *) "INPUT DRY AIR VOLUME MIXING RATIOS OF ABSORBERS:"
         do n = 1, size(win_ini)
            do j = 1, win_ini(n)%ntype
               write (*, *) "***ABSORBER TYPE", win_ini(n)%type_x(j), sum(win(n)%dv_x(:, j))/sum(dvair)
            end do
         end do
      end if

      !*** Put absorber in type
      call get_absorbers(win_ini, abs, ierr)

      !*** Calculate synthetic spectrum
      call calculate_syn_spectrum( &
         aero_lut, cirrus_lut, &
         path%output, runid, &
         aerosolfile, runpath, &
         win_ini, aerosol_ini, &
         grid%nlay, &
         flag, &
         abs, atm_rt, atm_xs, dvair, &
         z_tropopause, z_bl, &
         simulation, simulation_hi, win, aerosol, meta, synsettings, atm_scenario%surface_wspeed, exitflag)

      if (exitflag == 0) then
         if (flag%atm == 4) then !*** netcdf output
            !*** Write synthetic spectrum to output file
            spectrum_file = trim(path%spectrum)//trim(atm(i))//'_sim_'//runidstring//'.nc'
            call output_l1b_nc(simulation, win_ini, win, meta, spectrum_file, pixelid(i))

            if (flag%output > 0) then
               !*** ascii output for testing
               write (filename, '(i6.6)') pixelid(i) !*** groundpixel identifier
               spectrum_file = trim(path%spectrum)//'L1B_'//trim(filename)
               call output_l1b(simulation, win_ini, meta, spectrum_file)
               !*** Write meteo data to output file
               meteo_file = trim(path%spectrum)//'ATM_'//trim(filename)
               call output_atm( &
                  atm_scenario, meta, simulation(1)%sza, &
                  simulation(1)%iza, simulation(1)%phi, win_ini, win, aerosol, meteo_file)
               meteo_file = trim(path%spectrum)//'METEO_'//trim(filename)
               call output_meteo(atm_scenario, win_ini, win, aerosol, meteo_file)
            end if

         elseif (flag%atm == 5 .or. flag%atm == 6) then !*** netcdf output for Indianapolis data
            j1 = INDEX(atm(i), '.nc')
            index_info = atm(i) (j1 + 4:)

            !*** Write synthetic spectrum to output file
            spectrum_file = trim(path%spectrum)//"L1B_"//runidstring//'.nc'
            !call output_l1b_nc_js(simulation, meta, spectrum_file, index_info)
            call output_l1b_nc_ls(simulation, meta, spectrum_file, index_info)

            !*** Write meteo data to output file
            meteo_file = trim(path%spectrum)//"ATM_"//runidstring//'.nc'
            call output_atm_nc_js(atm_scenario, simulation, win, meta, meteo_file, index_info)

            ! !*** Write ascii output for testing purposes
            ! if (flag%output > 0) then
            !    !*** Write synthetic spectrum to output ascii file
            !    !*** spectrum_file = trim(path%spectrum)//"L1B_"//runidstring//'.dat'
            !    spectrum_file = trim(path%spectrum)//'L1B_'//trim(atm(i))
            !    call output_l1b(simulation, win_ini, meta, spectrum_file)

            !    !*** Write meteo data to output ascii file
            !    !*** meteo_file = trim(path%spectrum)//"ATM_"//runidstring//'.dat'
            !    meteo_file = trim(path%spectrum)//'ATM_'//trim(atm(i))
            !    call output_atm(atm_scenario, meta, simulation(1)%sza, &
            !         simulation(1)%iza, simulation(1)%phi, win_ini, win, aerosol, meteo_file)
            ! endif

         elseif (flag%atm == 7) then !*** netcdf output for LUT

            !*** Write high-resolution synthetic spectrum and meteo-data to output file
            lut_file = trim(path%spectrum)//"LUT_"//runidstring//'.nc'
            call output_lut_nc_js(simulation_hi, atm_scenario, win, meta, xco2, co2_scaling, lut_file)

         else                  !*** ascii output
            if (flag%atm == 2) then
               filename = trim(atm_name(29:38))//trim(atm_name(42:len(trim(atm_name)) - 4))
            else
               filename = trim(atm(i))
            end if
            !*** Write synthetic spectrum to output file
            spectrum_file = trim(path%spectrum)//'L1B_'//trim(filename)
            call output_l1b(simulation, win_ini, meta, spectrum_file)
            !*** Write meteo data to output file
            !*** ATM-file
            meteo_file = trim(path%spectrum)//'ATM_'//trim(filename)
            call output_atm( &
               atm_scenario, meta, simulation(1)%sza, &
               simulation(1)%iza, simulation(1)%phi, win_ini, win, aerosol, meteo_file)
            !*** METEO-file
            meteo_file = trim(path%spectrum)//'METEO_'//trim(filename)
            call output_meteo(atm_scenario, win_ini, win, aerosol, meteo_file)
         end if
      elseif (flag%output > 0) then
         print *, 'MAIN: no synthetic measurement created, errorflag = ', exitflag
      end if

      if (flag%output > 0) then
         print *, 'infile:'
         print *, trim(path%meteo)//trim(atm_name)
         print *
         print *, 'Spectrum ', i, ' of ', nfile
         write (*, '(a)') '***** Done.'
         write (*, '(a)') '*****************************************'
      end if

      !*** deallocate
      if (allocated(atm_scenario%z)) deallocate (atm_scenario%z)
      if (allocated(atm_scenario%p)) deallocate (atm_scenario%p)
      if (allocated(atm_scenario%t)) deallocate (atm_scenario%t)
      if (allocated(atm_scenario%h2o)) deallocate (atm_scenario%h2o)
      if (allocated(atm_scenario%co2)) deallocate (atm_scenario%co2)
      if (allocated(atm_scenario%ch4)) deallocate (atm_scenario%ch4)
      if (allocated(atm_scenario%co)) deallocate (atm_scenario%co)
      if (allocated(simulation)) deallocate (simulation)

   end do !i=1,nfile

   if (flag%atm == 4) then
      call close_miprep(atm(1), ierr)
      if (ierr .ne. 0) call stopretrieval("MAIN: error closing MIPrep files")
   end if

   if (flag%ilscalc .le. 0) then
      call synthetic_interface_close(ierr)
   end if

   !*** Mean CPU time per simulated spectra
   call cpu_time(stoptime)
   write (6, '(a, f8.3, a )') ' Mean cpu time per simulated spectra = ', (stoptime - starttime)/nfile, ' s'
   !*** Number of spectra
   write (6, '(a, i6 )') ' Number of synthetic spectra = ', nfile

end program main

