module calculate_syn_spectrum_module
   use header_module
   use read_synsettings_module
   use optic_input_module, only: aero_opt
   use optic_cirrus_module, only: cirrus_table, optic_cirrus_xs, optic_cirrus
   use OpticM_module
   use forward_model_noscat_module
   use forward_model_module
   use spec_interface_create_module
   use auxiliary_routines_module, only: check
   implicit none

!> @instrument noise levels
   type :: instrument_noise_levels
      real(double), dimension(:), allocatable :: SNR
      real(double), dimension(:), allocatable :: signal_shotnoise
      real(double) :: back_shotnoise
      real(double) :: dark_shotnoise
      real(double) :: readout_noise
      real(double) :: quantization_noise
   end type instrument_noise_levels

contains
!*******************************************************************************
!> @details Define geometry, SNR and microphysical properties of scatterers.
!! Compute synthetic spectrum
   subroutine calculate_syn_spectrum( &
      aero_lut, cirrus_lut, output_dir, runid, infile, runpath, win_ini, &
      aerosol_ini, nlay, flag, absorb, atm_rt, atm_xs, dvair, z_tropopause, &
      z_bl, simulation, simulation_hi, win, aerosol, meta, synsettings, &
      wspeed, errorflag)

      use netcdf

      type(Mie_lut), intent(in) :: aero_lut
      type(cirrus_table), intent(in) :: cirrus_lut
      character(len=*), intent(in) :: output_dir, infile, runpath
      type(window_ini), dimension(:), intent(inout) :: win_ini
      type(aero), dimension(:), intent(in) :: aerosol_ini
      type(spectrum), dimension(:), allocatable, intent(out) :: simulation
      type(spectrum), dimension(:), allocatable, intent(out) :: simulation_hi
      integer, intent(in) :: runid, nlay
      type(syn_instrument_noise_settings), intent(in) :: synsettings
      type(settings_flags), intent(in) :: flag
      type(atmosphere), intent(in) :: atm_rt
      type(atmosphere), intent(inout) :: atm_xs
      real(double), dimension(:), intent(inout) :: dvair
      type(absorbers), intent(in) :: absorb
      real(double), intent(in) :: z_tropopause, z_bl, wspeed
      type(window_spectrum), dimension(:), allocatable, intent(inout) :: win
      type(aero), dimension(:), allocatable, intent(out) :: aerosol
      type(metadata), intent(inout) :: meta
      integer, intent(out) :: errorflag

      !*** local
      integer :: i, j, k, l, n, errini, io, nrt, ierr, glintflag
      logical :: nanfound
      real(double), dimension(:), allocatable :: noise_rad, noise_irrad
      integer, dimension(:), allocatable :: nder_dum, indx
      real(double), dimension(:), allocatable :: nder_dbl
      real(double), dimension(3) :: rdn
      character(6):: runidstring
      real(double) :: rdummy1, rdummy2, rdummy3, rdummy4, rdummy5, rdummy6, &
                      rdummy7, rdummy8, rdummy9, rdummy10, rdummy11, rdummy12
      real(double) :: maxdummy, mindummy, albedo_nir, albedo_swir
      real(double), dimension(:), allocatable :: dvair_old  ! Partial air column
      real(double), dimension(:), allocatable :: play_old   ! Pressure, layer center
      real(double), dimension(:), allocatable :: vmr_h2o! water volumn mixing ratio
      integer ::  ntype_aer, nwin
      integer :: exitXSflag, maxotflag, minaotflag, mincotflag, maxd
      integer, parameter :: O2flag = 0, Tflag = 0
      real(double) :: cloud_fraction
      real(double), dimension(:, :), allocatable :: reflectance_hi, reflectance_hi_tmp
      type(derivatives) :: deriv_hi
      integer, dimension(:), allocatable :: nder
      real(double), dimension(:), allocatable :: radiance_hi
      real(double), dimension(:), allocatable :: radiance_lo
      real(double), dimension(:, :), allocatable :: albedo
      real(double), dimension(:, :, :), allocatable :: derivatives_hi     ! Hi-res reflectance derivatives wrt absorber partial columns
      real(double), dimension(:), allocatable :: derivT_hi              ! Hi-res reflectance derivatives wrt T profile
      real(double), dimension(:, :), allocatable :: derivP_hi            ! Hi-res reflectance derivatives wrt p profile
      real(double), dimension(:, :), allocatable :: deriv_albedo_hi      ! Hi-res reflectance derivatives wrt albedo
      character*2 :: ch
      logical :: exst, fluorescence
      logical, parameter :: deriv_flag = .false.
      type(aero), dimension(:), allocatable :: aerosol_tmp
      real(double):: s1, s2, s3

      type(instrument_noise_levels) :: instrument_noise
      integer :: nSB
      real(double) :: radiance_lo_cont, SB_lower, SB_upper, dSB
      real(double), dimension(:), allocatable :: SB   ! SB - Scene Brightness
      integer :: sx, sy, ncid, varid
      real(double) :: albedo_in(1)
      character(stringlen) :: fname_nc, info

      real(double), dimension(:), allocatable :: residual_image_lo

      real(double), dimension(:), allocatable :: sitf_signal
      real(double), dimension(:), allocatable :: sitf_response

      !***********************************************************************

      if (flag%output >= 2) then
         call writelog('*** Start of calculate_syn_spectrum ***', 1)
      end if

      if (nstokes == 3) then
         s1 = 1.
         s2 = 0.
         s3 = 0.
      else
         s1 = 1.d0
         s2 = 0.d0
         s3 = 0.d0
      end if

      !*** Initialize
      ierr = 0
      errorflag = 0
      maxotflag = 0
      nwin = size(win_ini)
      nrt = atm_rt%n
      call rdn01(3, rdn)

      glintflag = meta%oceanglint

      !*** Inquire if the aerosolfile exists
      inquire (FILE=trim(infile), EXIST=exst)

      !*** Get aerosol/cirrus optical properties
      if (flag%scat > 1 .and. flag%scat < 6 .and. exst .and. meta%landflag .ne. 1) then ! aerosol/cirrus optical properties from model and measurements
         call get_aerosol_properties_model( &
            aero_lut, cirrus_lut, infile, win_ini, atm_rt, flag%scat, aerosol)
      else if (flag%scat > 1 .and. flag%scat < 6) then
         call get_aerosol_properties_lognormal( &  ! aerosol parameters from namelist input file, these are ocean pixels
            aero_lut, cirrus_lut, &
            flag%output, glintflag, flag%glintscat, &
            nwin, z_tropopause, z_bl, &
            atm_rt, aerosol_ini, aerosol, ierr)
      else if (flag%scat == 6 .and. exst .and. meta%landflag .ne. 1) then ! aerosol properties from model, cloud properties from namelist
         call get_aerosol_properties_model( &
            aero_lut, cirrus_lut, infile, win_ini, atm_rt, flag%scat, aerosol_tmp)
         ntype_aer = size(aerosol_tmp)
         allocate (aerosol(ntype_aer + 1))
         aerosol(1:ntype_aer) = aerosol_tmp(1:ntype_aer)
         call get_aerosol_properties_lognormal( &
            aero_lut, cirrus_lut, &
            flag%output, glintflag, flag%glintscat, &
            nwin, z_tropopause, z_bl, &
            atm_rt, aerosol_ini, aerosol_tmp, ierr)
         aerosol(ntype_aer + 1) = aerosol_tmp(2)     ! second aerosoltype in namelist file describes clouds
         cloud_fraction = aerosol(ntype_aer + 1)%shapefrac
         aerosol(ntype_aer + 1)%shapefrac = 1.d0
      else if (flag%scat == 1 .or. flag%scat == 6) then ! For flag%scat==6, these are ocean pixels
         call get_aerosol_properties_lognormal( &  ! aerosol parameters from namelist input file
            aero_lut, cirrus_lut, &
            flag%output, glintflag, flag%glintscat, &
            nwin, z_tropopause, z_bl, &
            atm_rt, aerosol_ini, aerosol, ierr)
         if (flag%scat == 6) then         ! last aerosoltype is cloud, shapefrac is actually cloud fraction
            ntype_aer = size(aerosol)
            cloud_fraction = aerosol(ntype_aer)%shapefrac
            aerosol(ntype_aer)%shapefrac = 1.d0
         end if
      else if (flag%scat == 0) then
         allocate (aerosol(0))
      else
         if (flag%output > 1) then
            call writelog("CALCULATE_SYN_SPECTRUM: flag%scat not valid", 6)
         end if
         errorflag = -1
         goto 101
      end if
      ntype_aer = size(aerosol)

      !*** Get albedo
      allocate (albedo(nwin, 3))
      albedo = 0.D0
      if (flag%atm == 2 .or. (flag%atm == 4 .and. exst .and. meta%landflag .ne. 1)) then ! ECHAM input
         !*** Open atmosferic input files (*.out):
         !*** Get surface albedo from satellite data
         open (newunit(io), FILE=trim(infile), action='read', status='old')
         do i = 1, 5
            read (io, *)
         end do
         read (io, *) rdummy1, rdummy2, rdummy3, rdummy4, rdummy5, rdummy6, &
            rdummy7, rdummy8, rdummy9, rdummy10, rdummy11, rdummy12
         close (io)

         albedo_nir = rdummy2 ! Albedo for O2A-band: use later for fluoroscence determination
         albedo_swir = rdummy7*0.7 ! Albedo for 2.4um-band: use later for fluoroscence determination

         do n = 1, nwin
            if (win_ini(n)%albflag == 0) then ! surface albedo from satellite data
               ! The values that wave_start and wave_stop are compared
               ! to were originally given in wavenumbers [cm-1]. Since this
               ! version of RemoTeC works exclusively with wavelengths [nm]
               ! a conversion of the hardcoded parameters is performed here
               ! (actual content of wave_start and wave_stop is not touched
               ! though).
               if (win_ini(n)%wave_stop < 1.d7/20000.) then
                  albedo(n, 1) = rdummy3
               else if (win_ini(n)%wave_stop < 1.d7/16000. .and. win_ini(n)%wave_start >= 1.d7/20000.) then
                  albedo(n, 1) = rdummy4
               else if (win_ini(n)%wave_stop < 1.d7/13400. .and. win_ini(n)%wave_start >= 1.d7/16600.) then
                  albedo(n, 1) = rdummy1
               else if (win_ini(n)%wave_stop <= 1.d7/12700. .and. win_ini(n)%wave_start >= 1.d7/13400.) then
                  albedo(n, 1) = rdummy2
               else if (win_ini(n)%wave_stop <= 1.d7/5800. .and. win_ini(n)%wave_start >= 1.d7/6500.) then
                  albedo(n, 1) = rdummy6
               else if (win_ini(n)%wave_stop <= 1.d7/5000. .and. win_ini(n)%wave_start >= 1.d7/5200.) then
                  albedo(n, 1) = rdummy7
               else if (win_ini(n)%wave_stop <= 1.d7/4700. .and. win_ini(n)%wave_start >= 1.d7/5200.) then
                  albedo(n, 1) = rdummy7
               else if (win_ini(n)%wave_stop <= 1.d7/4100. .and. win_ini(n)%wave_start >= 1.d7/4500.) then
                  albedo(n, 1) = 0.7*rdummy7
                  call writelog('CALCULATE_SYN_SPECTRUM: used proxy &
                      &albedo@2.3mu = 0.7*albedo@2.06mu because of lacking SCIAMACHY data', 5)
               else
                  if (flag%output > 2) then
                     call writelog('CALCULATE_SYN_SPECTRUM: surface albedo not defined for this window', 6)
                  end if
                  errorflag = 1
                  goto 101
               end if
            else ! surface albedo from input file settings**.in
               do i = 1, max(1, win_ini(n)%albflag)
                  albedo(n, i) = win_ini(n)%albedo(i)
               end do
            end if
         end do
      else if (flag%atm == 5 .or. flag%atm == 6) then ! Read albedo from Indianapolis-file (same for older Gaussian plume based files and new ICON based files)
         ! Extract name of file with high-resolution input data
         i = INDEX(infile, '.nc')
         fname_nc = infile(1:i + 2)
         info = infile(i + 3:)

         ! Extract index in x-dimension to be read
         i = INDEX(info, 'X')
         read (info(i + 1:i + 3), '(I3)') sx

         ! Extract index in y-dimension to be read
         i = INDEX(info, 'Y')
         read (info(i + 1:i + 3), '(I3)') sy

         do n = 1, nwin
            if (win_ini(n)%albflag == 0) then   !surface albedo from satellite data
               !*** Open NetCDF file
               call check(nf90_open(trim(fname_nc), nf90_nowrite, ncid), ierr)
               if (ierr .ne. 0) return

               if (win_ini(n)%wave_stop <= 1.d7/4700. .and. win_ini(n)%wave_start >= 1.d7/5200.) then ! 2000nm
                  ! Get variable-id for albedo data at 2000nm
                  call check(NF90_INQ_VARID(ncid, "albedo_2000nm", varid), ierr)
                  if (ierr .ne. 0) return
               else
                  if (flag%output > 2) then
                     call writelog('CALCULATE_SYN_SPECTRUM: high-resolution surface albedo data not defined for this window', 6)
                  end if

                  call check(nf90_close(ncid), ierr)
                  if (ierr .ne. 0) return
                  errorflag = 1
                  goto 101
               end if

               ! Read albedo data from netcdf-file
               call check(NF90_GET_VAR(ncid, varid, albedo_in, [sx, sy], [1, 1]), ierr)
               if (ierr .ne. 0) return
               albedo(n, 1) = albedo_in(1)

               ! Close netcdf-file
               call check(nf90_close(ncid), ierr)
               if (ierr .ne. 0) return
            else ! surface albedo from input file settings**.in
               do i = 1, max(1, win_ini(n)%albflag)
                  albedo(n, i) = win_ini(n)%albedo(i)
               end do
            end if
         end do
      else ! no satellite input
         do n = 1, nwin
            do i = 1, max(1, win_ini(n)%albflag)
               albedo(n, i) = win_ini(n)%albedo(i)
            end do
         end do
      end if

      !*** JS: For the G1 1600nm setup, negative albedo values are present leading to invalid operations
      !*** (square-root of negative numbers) and hence job abortion. Add those cases to faulty spectra
      if (minval(albedo) < 0.0) then
         if (flag%output > 2) then
            call writelog('CALCULATE_SYN_SPECTRUM: negative surface albedo for this window', 6)
         end if
         errorflag = 8
         goto 101
      end if
      !************************************************************************************************

      !*** Turn Fluorescence flag on for land scenes with vegetation (albedo_NIR/albedo_SWIR) > 5)
      fluorescence = .false.

      if (exst .and. meta%landflag .ne. 1) then
         if (albedo_nir/albedo_swir > 5.d0 .and. win_ini(1)%Fsflag > 0) then
            fluorescence = .true.
         end if
      end if

      call get_geometry_sunsync(meta%time, meta%lat(1), meta%altitude, &
                                meta%szaflag, meta%sza, meta%iza, meta%saz, meta%iaz, meta%phi, errorflag)

      if (errorflag .ne. 0 .or. meta%sza > 90.d0) then
         errorflag = 2
         if (flag%output > 1) then
            call writelog('CALCULATE_SYN_SPECTRUM: no sun-synchronous orbit for given latitude or sza>90', 6)
         end if
         goto 101
      end if

      !**********   Reading/computing input data finished. Now calculate spectrum. ****************************
      !*** Flags used for retrieval, set to zero here
      do k = 1, ntype_aer
         aerosol(k)%AerosolFlags = 0
      end do
      !*** Determine layers with non-negligible aerosol optical depth
      if (allocated(nder_dum)) deallocate (nder_dum)
      allocate (nder_dum(sum(aerosol(:)%maxd)))
      j = 0
      do k = 1, ntype_aer
         do i = 1, aerosol(k)%maxd
            j = j + 1
            nder_dum(j) = aerosol(k)%nder(i)
         end do
      end do
      n = 0
      do i = 1, j
         l = nder_dum(i)
         do k = i + 1, j
            if (nder_dum(k) == l .and. nder_dum(k) /= 0) then
               nder_dum(k) = 0
               n = n + 1
            end if
         end do
      end do
      maxd = j - n

      if (allocated(nder)) deallocate (nder)
      if (allocated(nder_dbl)) deallocate (nder_dbl)
      if (allocated(indx)) deallocate (indx)
      allocate (nder(maxd))
      allocate (nder_dbl(maxd))
      allocate (indx(maxd))

      nder = pack(nder_dum, nder_dum /= 0)
      nder_dbl = dble(nder)
      call sort(nder_dbl, maxd, indx)
      nder = int(nder_dbl)

      !*** [JS] Store high-resolution spectra for generation of LUT
      if (flag%atm == 7) then
         if (allocated(simulation_hi)) deallocate (simulation_hi)
         allocate (simulation_hi(nwin))
      end if

      if (allocated(simulation)) deallocate (simulation)
      allocate (simulation(nwin))
      allocate (play_old(atm_xs%n), dvair_old(atm_xs%n), vmr_h2o(atm_xs%n))
      !*** Partial air column
      dvair_old = dvair
      play_old = atm_xs%p

      !*** Get vmr of H2O in xs-layers
      vmr_h2o = 0.d0
      do n = 1, nwin
         do i = 1, win_ini(n)%ntype
            if (abs(win_ini(n)%xsdb(i)%species) == 1 .or. (abs(win_ini(n)%xsdb(i)%species) >= 100 .and. abs(win_ini(n)%xsdb(i)%species) <= 199)) then
               vmr_h2o(:) = win(n)%dv_x(:, i)/dvair(:) ! Use initial VMR
               goto 102
            end if
         end do
      end do

102   continue
      !*** Calculate synthetic spectrum
      do n = 1, nwin
         if (allocated(reflectance_hi)) deallocate (reflectance_hi)
         if (allocated(reflectance_hi_tmp)) deallocate (reflectance_hi_tmp)
         if (allocated(noise_rad)) deallocate (noise_rad)
         if (allocated(noise_irrad)) deallocate (noise_irrad)
         if (allocated(win(n)%reflectance_lo)) deallocate (win(n)%reflectance_lo)
         if (allocated(win(n)%x_molec)) deallocate (win(n)%x_molec)
         if (allocated(win(n)%albedo)) deallocate (win(n)%albedo)
         allocate (simulation(n)%radiance(win(n)%nwave_lo), &
                   simulation(n)%radiance_noise(win(n)%nwave_lo), &
                   simulation(n)%radiance_error(win(n)%nwave_lo), &
                   simulation(n)%irradiance(win(n)%nwave_lo), &
                   simulation(n)%irradiance_noise(win(n)%nwave_lo), &
                   simulation(n)%irradiance_error(win(n)%nwave_lo), &
                   simulation(n)%wavelength(win(n)%nwave_lo), &
                   reflectance_hi(win_ini(n)%nwave_hi, nstokes), &
                   reflectance_hi_tmp(win_ini(n)%nwave_hi, nstokes), &
                   win(n)%reflectance_lo(win(n)%nwave_lo), &
                   noise_rad(win(n)%nwave_lo), &
                   noise_irrad(win(n)%nwave_lo), &
                   win(n)%x_molec(atm_xs%n, win_ini(n)%ntype), &
                   win(n)%albedo(max(1, win_ini(n)%albflag)))
         win(n)%x_molec = win(n)%dv_x
         win(n)%albedo = 0.D0
         do i = 1, max(1, win_ini(n)%albflag)
            win(n)%albedo(i) = albedo(n, i)
         end do

         if (flag%atm == 7) then
            allocate (simulation_hi(n)%radiance(win_ini(n)%nwave_hi), &
                      simulation_hi(n)%wavelength(win_ini(n)%nwave_hi))
         end if

         !*** calculate high resolution and derivatives
         ! derivatives not needed, to remove)
         allocate (deriv_hi%densmol(win_ini(n)%nwave_hi, nlay, win_ini(n)%ntype, nstokes))
         allocate (deriv_hi%alb(win_ini(n)%nwave_hi, max(1, win_ini(n)%albflag), nstokes))
         allocate (deriv_hi%T(win_ini(n)%nwave_hi, nstokes))
         allocate (deriv_hi%P(win_ini(n)%nwave_hi, nlay, nstokes))
         allocate (deriv_hi%aerosol(win_ini(n)%nwave_hi, 3, 1, nstokes))
         allocate (deriv_hi%Fs(win_ini(n)%nwave_hi, max(0, win_ini(n)%Fsflag)))

         simulation(:)%sza = meta%sza
         simulation(:)%iza = meta%iza
         simulation(:)%phi = meta%phi

         if (flag%scat == 0) then
            allocate (derivatives_hi(win_ini(n)%nwave_hi, nlay, win_ini(n)%ntype))
            allocate (deriv_albedo_hi(win_ini(n)%nwave_hi, win_ini(n)%albflag))
            allocate (derivT_hi(win_ini(n)%nwave_hi))
            allocate (derivP_hi(win_ini(n)%nwave_hi, nlay))
            call forward_model_hi_noscat( &
               flag%xs, flag%O2, flag%temp, glintflag, &
               meta%sza, meta%iza, meta%phi, wspeed, n, &
               absorb, atm_xs, dvair, dvair_old, vmr_h2o, play_old, &
               win_ini, win, reflectance_hi, &
               derivatives_hi, deriv_albedo_hi, derivT_hi, derivP_hi, ExitXSFlag, ierr)
            deallocate (derivatives_hi)
            deallocate (derivT_hi)
            deallocate (derivP_hi)
            deallocate (deriv_albedo_hi)
         else
            win(n)%Fs = 0.d0
            if (.not. allocated(win(n)%sun_spectrum_ref_hi)) then
               allocate (win(n)%sun_spectrum_ref_hi(win_ini(n)%nwave_hi))
            end if
            win(n)%sun_spectrum_ref_hi = win_ini(n)%sun_spectrum_ref_hi

            !*** For scenes with vegetation:
            if (win_ini(n)%Fsflag == 1 .and. fluorescence) then
               win(n)%Fs(1) = 3.798d12 ! in photons s-1 cm-2 sr-1 nm-1 = 10 W m-2 sr-1 um-1 ! TODO: check FS
            else if (win_ini(n)%Fsflag == 2 .and. fluorescence) then
               win(n)%Fs(1) = 0.01d0   ! in fraction of continuum
            end if

            call forward_model_hi(aero_lut, cirrus_lut, &
                                  flag, glintflag, &
                                  meta%sza, meta%iza, meta%phi, wspeed, &
                                  n, &
                                  atm_rt%n, &
                                  nlay, &
                                  absorb, &
                                  atm_xs, &
                                  dvair, &
                                  dvair_old, &
                                  vmr_h2o, &
                                  play_old, &
                                  win_ini(n), &
                                  win(n), &
                                  win_ini(n)%wavelength_hi, &
                                  reflectance_hi, &
                                  aerosol, &
                                  minaotflag, &
                                  mincotflag, &
                                  nder, &
                                  deriv_hi, &
                                  ExitXSFlag, &
                                  MaxOTFlag, &
                                  ierr, &
                                  deriv_flag)

            if (maxotflag == 1) then
               call writelog('CALCULATE_SYN_SPECTRUM: AOT+COT above threshold', 6)
               errorflag = 4
               goto 101
            end if

            if (exitxsflag == 1) then
               call writelog('CALCULATE_SYN_SPECTRUM: error in XSDB', 6)
               errorflag = 3
               goto 101
            end if

            !***
!!$          write(ch,'(i2.2)')n
!!$          open(newunit(io),FILE=trim(output_dir)//'spectrum_cloud_'//ch//'.dat')
!!$          write(io,'(A)')'#Wavelength/nm         Reflectance'
!!$          do k = 1, win_ini(n)%nwave_hi
!!$             write(io,'(2(1pE23.15E3,x))') &
!!$                  win_ini(n)%wavelength_hi(k),&
!!$                  reflectance_hi(k)
!!$          enddo
!!$          close(io)
            !***

            !*** Independent pixel approximation for partially clouded pixels
            if (flag%scat == 6 .and. cloud_fraction > 1.d-6 .and. cloud_fraction < 1.d0 - 1d-6) then
               call forward_model_hi(aero_lut, cirrus_lut, &
                                     flag, glintflag, &
                                     meta%sza, meta%iza, meta%phi, wspeed, &
                                     n, &
                                     atm_rt%n, &
                                     nlay, &
                                     absorb, &
                                     atm_xs, &
                                     dvair, &
                                     dvair_old, &
                                     vmr_h2o, &
                                     play_old, &
                                     win_ini(n), &
                                     win(n), &
                                     win_ini(n)%wavelength_hi, &
                                     reflectance_hi_tmp, &
                                     aerosol(1:ntype_aer - 1), &   !*** minus clouds
                                     minaotflag, &
                                     mincotflag, &
                                     nder, &
                                     deriv_hi, &
                                     ExitXSFlag, &
                                     MaxOTFlag, &
                                     ierr, &
                                     deriv_flag)

               if (maxotflag == 1) then
                  call writelog('CALCULATE_SYN_SPECTRUM: AOT+COT above threshold', 6)
                  errorflag = 4
                  goto 101
               end if

               if (exitxsflag == 1) then
                  call writelog('CALCULATE_SYN_SPECTRUM: error in XSDB', 6)
                  errorflag = 3
                  goto 101
               end if

               reflectance_hi = cloud_fraction*reflectance_hi + (1.d0 - cloud_fraction)*reflectance_hi_tmp

!!$             open(newunit(io),FILE=trim(output_dir)//'spectrum_nocloud_'//ch//'.dat')
!!$             write(io,'(A)')'#Wavelength/cm-1         Reflectance'
!!$             do k = 1, win_ini(n)%nwave_hi
!!$                write(io,'(2(1pE23.15E3,x))') &
!!$                     win_ini(n)%wavelength_hi(k),&
!!$                     reflectance_hi_tmp(k)
!!$             enddo
!!$             close(io)
            end if ! clouds
         end if ! flag%scat != 0

         deallocate (deriv_hi%densmol)
         deallocate (deriv_hi%alb)
         deallocate (deriv_hi%T)
         deallocate (deriv_hi%P)
         deallocate (deriv_hi%aerosol)
         deallocate (deriv_hi%Fs)

         !*** Convolve hi-res reflectance spectrum -> lo-res reflectance spectrum
         if (allocated(radiance_hi)) deallocate (radiance_hi)
         if (allocated(radiance_lo)) deallocate (radiance_lo)
         allocate (radiance_hi(size(reflectance_hi(:, 1))))
         allocate (radiance_lo(win(n)%nwave_lo))

         !*** Spectral response is calculated for radiance (not reflectance)
         !radiance_hi = reflectance_hi(:,1)*win_ini(n)%sun_spectrum_ref_hi
         if (nstokes == 1) then
            radiance_hi = reflectance_hi(:, 1)*win_ini(n)%sun_spectrum_ref_hi
         else
          radiance_hi = (reflectance_hi(:, 1)*s1 + reflectance_hi(:, 2)*s2 + reflectance_hi(:, 3)*s3)*win_ini(n)%sun_spectrum_ref_hi
         end if

         !*** CONTROL OUTPUT
         if (flag%output == 3) then
            write (ch, '(i2.2)') n
            open (newunit(io), FILE=trim(output_dir)//'spectrum_hi_'//ch//'.dat')
            if (nstokes == 1) then
               write (io, '(A)') '#Wavelength/nm           Reflectance             Radiance                Irradiance'
               do k = 1, win_ini(n)%nwave_hi
                  write (io, '(6(1pE23.15E3,x))') &
                     win_ini(n)%wavelength_hi(k), &
                     reflectance_hi(k, 1), radiance_hi(k), win_ini(n)%sun_spectrum_ref_hi(k)
               end do
            else
                write(io, '(A)')'#Wavelength/nm           Reflectance1            Reflectance2            Reflectance3            Radiance                Irradiance'
               do k = 1, win_ini(n)%nwave_hi
                  write (io, '(6(1pE23.15E3,x))') &
                     win_ini(n)%wavelength_hi(k), &
                 reflectance_hi(k, 1), reflectance_hi(k, 2), reflectance_hi(k, 3), radiance_hi(k), win_ini(n)%sun_spectrum_ref_hi(k)
               end do
            end if
            close (io)
         end if

         !*** Use stored response function
         call spectral_response_stored( &
            win(n)%resp_store, &
            win(n)%ie_store, &
            win(n)%is_store, &
            radiance_hi, &
            radiance_lo)

         !*** straylight offset
         if (synsettings%straylight_flag == 1) then ! add constant straylight offset
            radiance_lo = radiance_lo + synsettings%straylight_offset*maxval(radiance_lo)
         end if

         !*** residual image
         if (synsettings%residual_image_flag == 1) then ! add residual image spectrum
            call read_residual_image(win(n)%wavelength_lo, win(n)%nwave_lo, radiance_lo, synsettings%residual_image_path, synsettings%residual_image_intensity, residual_image_lo)
            radiance_lo = radiance_lo + residual_image_lo
         end if

         !*** signal intensity transfer function
         !*** here, the sitf is in a form where it maps from spectral radiances to spectral radiances
         if (synsettings%sitf_flag == 1) then ! modify low resolution spectrum using signal intensity transfer function
            call read_signal_intensity_transfer_function(synsettings%sitf_path, sitf_signal, sitf_response)
            call apply_signal_intensity_transfer_function(radiance_lo, sitf_signal, sitf_response)
         end if

         !*** CALCULATE NOISE

         !*** noise on irradiance
         if (synsettings%noiseflag_irrad(n) == -1) then ! no noise
            noise_irrad = 0.d0
         else if (synsettings%noiseflag_irrad(n) == 0) then ! shot noise
            do k = 1, win(n)%nwave_lo
               noise_irrad(k) = 1./synsettings%snr_irrad(n)*sqrt(maxval(win(n)%sun_spectrum_sat_lo)/win(n)%sun_spectrum_sat_lo(k))
            end do
         else
            call stopretrieval('CALCULATE_SYN_SPECTRUM: No noise model implemented for irradiance')
         end if

         !*** noise on radiance
         if (synsettings%noiseflag_rad(n) == -1) then ! no noise
            noise_rad = 0.d0
         else if (synsettings%noiseflag_rad(n) == 0) then ! shot noise scaled wrt. reference scene (reference albdeo and SZA from syn_create_RunID.nml)
            !*** simple shot noise model:
            !*** Represents the noise (in relative units -> divide by 'radiance_lo(k)') of the reflected solar radiance at given reference albedo and SZA
            !*** It is assumed that the noise scales with the sqrt(radiance) if the actual albedo and SZA are not equal to the corresponding reference values
            do k = 1, win(n)%nwave_lo
               noise_rad(k) = 1./synsettings%snr_rad(n)* &
                              sqrt(synsettings%ref_alb/pi*cos(synsettings%ref_sza/180.*pi)* &
                                   maxval(win(n)%sun_spectrum_sat_lo)/radiance_lo(k))
            end do
         else if (synsettings%noiseflag_rad(n) == 1) then ! CO2IMAGE noise model
            if (allocated(instrument_noise%SNR)) deallocate (instrument_noise%SNR)
            if (allocated(instrument_noise%signal_shotnoise)) deallocate (instrument_noise%signal_shotnoise)
            allocate (instrument_noise%SNR(win(n)%nwave_lo), instrument_noise%signal_shotnoise(win(n)%nwave_lo))

          call co2image_noise_model(win_ini(n), synsettings, radiance_lo, win(n)%sun_spectrum_sat_lo, flag%output, instrument_noise)
            noise_rad = 1/instrument_noise%SNR(:)

            if (flag%output >= 2) then
               !*** Write noise levels to file
               write (ch, '(i2.2)') n
               open (newunit(io), FILE=trim(output_dir)//'instrument_noise_levels_'//ch//'.dat')
                write(io, '(A)')'# Background shot noise [# electrons]    Dark current shot noise [# electrons]    Readout noise [# electrons]    Quantization noise [# electrons]'
               write (io, '(A)') '# Wavelength [nm]    Signal shot noise [# electrons]    SNR [-]'
                write(io, '(4(1pE23.15E3,x))') instrument_noise%back_shotnoise, instrument_noise%dark_shotnoise, instrument_noise%readout_noise, instrument_noise%quantization_noise
               do k = 1, win(n)%nwave_lo
              write (io, '(3(1pE23.15E3,x))') win(n)%wavelength_lo(k), instrument_noise%signal_shotnoise(k), instrument_noise%SNR(k)
               end do
               close (io)
            end if
         else if (synsettings%noiseflag_rad(n) == 2) then ! CO2M noise model
            do k = 1, win(n)%nwave_lo
                noise_rad(k) = DSQRT(synsettings%snr_rad_a(n) * radiance_lo(k) + synsettings%snr_rad_b(n)**2) / (synsettings%snr_rad_a(n) * radiance_lo(k))
            end do
         else
            call writelog('CALCULATE_SYN_SPECTRUM: no noise model implemented for the given noise flag', 6)
            errorflag = 6
            goto 101
         end if

         !*** Add noise to the lo-res radiance spectrum
         errini = INT(rdn(2)*10000)
         simulation(n)%radiance_error = radiance_lo   !noise-free spectrum
         simulation(n)%radiance = radiance_lo
         call noise(-errini, simulation(n)%radiance, win(n)%nwave_lo, noise_rad)
         simulation(n)%radiance_error = simulation(n)%radiance - simulation(n)%radiance_error
         !*** add noise to lo-res solar irradiance
         errini = INT(rdn(3)*10000)
         simulation(n)%irradiance_error = win(n)%sun_spectrum_sat_lo  ! noise free irradiance
         simulation(n)%irradiance = win(n)%sun_spectrum_sat_lo
         call noise(-errini, simulation(n)%irradiance, win(n)%nwave_lo, noise_irrad)
         simulation(n)%irradiance_error = simulation(n)%irradiance - simulation(n)%irradiance_error
         !*** Calculate noise in units of radiance/irradiance
         do k = 1, win(n)%nwave_lo
            simulation(n)%radiance_noise(k) = noise_rad(k)*simulation(n)%radiance(k)
            simulation(n)%irradiance_noise(k) = noise_irrad(k)*simulation(n)%irradiance(k)
         end do
         simulation(n)%nwave = win(n)%nwave_lo
         simulation(n)%wavelength = win(n)%wavelength_lo

         !*** Store high-resolution spectra for LUT
         if (flag%atm == 7) then
            simulation_hi(n)%nwave = win_ini(n)%nwave_hi
            simulation_hi(n)%wavelength = win_ini(n)%wavelength_hi
            simulation_hi(n)%radiance = radiance_hi
         end if

         if (flag%output >= 2 .AND. synsettings%noiseflag_rad(n) == 1 .AND. flag%scat == 0) then
            !*** JS:  Write additional information about noise levels and continuum SNR as a function of scence brightness:
            !***      Take the low-resolution refelcted radiance and scale it by the scene's albedo, SZA and pi, such that it
            !***      represents the incoming solar irradiance, but corrected for the two-way atmospheric absorption. Then scale
            !***      by a range of different scene brightness conditions (ALB*cos(SZA)/pi) and calculate the corresponding SNR and noise levels.
            !***      Only allow for non-scattering retrievals for now, not sure if 'radiance_lo' is as easily scalable with a scattering atmosphere

            !*** Find continuum (maximum) radiance and 'remove' the effect of albedo, SZA and spread across the solid angle pi
            radiance_lo_cont = maxval(radiance_lo)
            radiance_lo_cont = radiance_lo_cont/win(n)%albedo(1)/dcos(meta%sza/180.*pi)*pi

            !*** Create a range of scene brightness values
            SB_lower = 0.001 ! Lower limit of scece brightness
            SB_upper = 0.30  ! Upper limit of scece brightness
            dSB = 0.001      ! Scene brightness step size
            nSB = (SB_upper - SB_lower)/dSB + 1

            if (allocated(SB)) deallocate (SB)
            allocate (SB(nSB))

            do k = 1, nSB
               SB(k) = SB_lower + (k - 1)*dSB
            end do

            if (allocated(instrument_noise%SNR)) deallocate (instrument_noise%SNR)
            if (allocated(instrument_noise%signal_shotnoise)) deallocate (instrument_noise%signal_shotnoise)
            allocate (instrument_noise%SNR(nSB), instrument_noise%signal_shotnoise(nSB))

            !*** Calucluate SNR and individual noise levels using instrument noise model
            call co2image_noise_model(win_ini(n), synsettings, radiance_lo_cont * SB(:), SB(:) * 0. + maxval(win(n)%sun_spectrum_sat_lo), flag%output, instrument_noise)

            !*** Write data to ascii file
            write (ch, '(i2.2)') n
            open (newunit(io), FILE=trim(output_dir)//'instrument_continuum_noise_levels_vs_scene_brightness_'//ch//'.dat')
            write(io, '(A)')'# Background shot noise [# electrons]    Dark current shot noise [# electrons]    Readout noise [# electrons]    Quantization noise [# electrons]'
            write (io, '(A)') '# Scene Brightness [-]    Signal shot noise [# electrons]    SNR [-]'
            write(io, '(4(1pE23.15E3,x))') instrument_noise%back_shotnoise, instrument_noise%dark_shotnoise, instrument_noise%readout_noise, instrument_noise%quantization_noise
            do k = 1, nSB
               write (io, '(3(1pE23.15E3,x))') SB(k), instrument_noise%signal_shotnoise(k), instrument_noise%SNR(k)
            end do
            close (io)
         end if
      end do ! n=1, nwin

101   continue
      !*** Write faulty or missing synthetic spectrum to file:
      !*** (1) surface albedo not defined for this window
      !*** (2) no sun-synchronous orbit
      !*** (3) error in XSDB
      !*** (4) MAXOT reached
      !*** (5) no noise model for given SNR
      !*** (6) no noise model for spectral range
      !*** (7) reflectances outside the 0..1 interval
      !*** (8) negative albedo for this window
      write (runidstring, '(I6.6)') runid
      open (newunit(io), FILE=trim(output_dir)//'faulty_spectra_'//runidstring//'.dat', ACCESS='APPEND')

      if (errorflag .ne. 0) then
         write (io, '(a, i3)') trim(infile), errorflag
      else
         maxdummy = 1.0d0
         mindummy = 1d-6
         nanfound = .false.
         do n = 1, nwin
            if (maxval(simulation(n)%radiance/simulation(n)%irradiance) .gt. maxdummy) &
               maxdummy = maxval(simulation(n)%radiance/simulation(n)%irradiance)
            if (minval(simulation(n)%radiance/simulation(n)%irradiance) .lt. mindummy) &
               mindummy = minval(simulation(n)%radiance/simulation(n)%irradiance)
            do l = 1, win(n)%nwave_lo
               if (simulation(n)%radiance(l) .ne. simulation(n)%radiance(l)) nanfound = .true.
            end do
         end do
         if (maxdummy .gt. 1.0d0 .or. mindummy .le. 0.0d0 .or. nanfound) then
            errorflag = 7
            call writelog('in if loop '//trim(output_dir)//'faulty_spectra_'//runidstring//'.dat', 1)
            write (io, '(a, a3, 2E14.5)') trim(infile), '7', maxdummy, mindummy
         end if
      end if

      close (io)

      if (flag%output >= 2) then
         call writelog('*** End of calculate_syn_spectrum ***', 1)
      end if
   end subroutine calculate_syn_spectrum

!*******************************************************************************

   subroutine orbit_sim_sunsync(year, month, day, hour, minu, sec, lat, altitude, sza, ierr)
      implicit none
      integer, intent(in) :: year, month, day, hour, minu, sec
      real(double), intent(in) :: lat, altitude
      !*** output
      integer, intent(out) :: ierr
      real(double), intent(out) :: sza
      !*** local variables
      integer :: j, julianday, LMT, LMT1, LMT2
      real(double) ::  LTAN, degree, kh, lat_rad, n, l_sza, g_sza, lambda_sza, epsilon_sza, y_sza, x_sza
      real(double) :: alpha_sza, delta_sza, ex, i_sza, DT, ha, sza1, sza2
      real(double), parameter :: radius = 6367.d0   !Earth radius in km

      ierr = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MORE ACCURATE CALCULATION FOR SZA AT 13:30 FOR 16th DAY OF MONTH:
      ! Julian day number at 12.00 UT:
      if (month .gt. 2) then
         j = 0
      else
         j = 1
      end if

    julianday = 1721089 + day + floor(-0.75*floor((year - j)/100.0)) + floor(365.25*(year - j)) + floor(367.0*(j+(dble(month)-2.0)/12.0))
      !*** For a sun-synchronous orbit, given Julian day number JD, latitude (deg), local mean time of the ascending node LTAN (hr)
      !*** and altitude (km), produce the solar zenith angle (deg) and LMT crossing time (hr) at the given latitude on the day side
      !*** (nadir viewing). Day side is where the solar zenith angle is smallest.
      LTAN = dble(hour) + dble(minu)/60.d0 + dble(sec)/3600.d0   ! time in hr
      degree = 180.0d0/PI
      kh = 10.10949
      lat_rad = lat/degree       ! latitude in radians
      n = julianday - 2.451545d6 ! days elapsed since January 1, 2000, 12:00 UT
      l_sza = mod(280.46 + 0.9856474*n, 360.0)
      g_sza = mod(357.528 + 0.9856003*n, 360.0)/degree
      lambda_sza = (l_sza + 1.915*sin(g_sza) + 0.02*sin(2.0*g_sza))/degree
      epsilon_sza = (23.439 - n/2500000.)/degree
      y_sza = cos(epsilon_sza)*sin(lambda_sza)
      x_sza = cos(lambda_sza)
      alpha_sza = mod(atan2(y_sza, x_sza)*Degree, 360.0)
      delta_sza = asin(sin(epsilon_sza)*sin(lambda_sza))
      ex = 4.0*(l_sza - alpha_sza)
      i_sza = acos(-(1.0 + altitude/radius)**3.5/kh)
      if (abs(tan(lat_rad)/tan(i_sza)) .ge. 1.d0) then ! no sun-synchronous orbit for this latitude
         ierr = 1
         return
      end if
      DT = asin(tan(lat_rad)/tan(i_sza))/15.*Degree
      LMT1 = mod(LTAN + DT, 24.0)
      LMT2 = mod(LTAN - 12.0 - DT, 24.0)
      ha = (LMT1 - 12.0 + ex/60.0)*15.0/Degree
      sza1 = acos(cos(ha)*cos(delta_sza)*cos(lat_rad) + sin(delta_sza)*sin(lat_rad))
      ha = (LMT2 - 12.0 + ex/60.0)*15.0/Degree
      sza2 = acos(cos(ha)*cos(delta_sza)*cos(lat_rad) + sin(delta_sza)*sin(lat_rad))
      if (sza1 < sza2) then
         SZA = sza1*degree
         LMT = LMT1
      else
         SZA = sza2*degree
         LMT = LMT2
      end if
   end subroutine orbit_sim_sunsync

!-------------------------------------------------------------------------------

   subroutine get_aerosol_properties_model(aero_lut, cirrus_lut, infile, win_ini, atm_rt, scatflag, aerosol)
      type(Mie_lut), intent(in) :: aero_lut
      type(cirrus_table), intent(in) :: cirrus_lut
      character(len=*), intent(in) :: infile
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(atmosphere), intent(in) :: atm_rt
      integer, intent(in) :: scatflag
      type(aero), dimension(:), allocatable, intent(out) :: aerosol
      !*** local
      integer :: i, j, k, n, io, nwin, ntype_aer, nrt
      real(double) :: dummy, sigma
      real(double), dimension(7) :: dummy7, reff, veff
      real(double), dimension(56) :: rm_input_array, fim_input_array
      real(double) :: aot_echam, aot_modis_mean, aot_modis_max, aot_modis_min, aot_modis_median
      real(double) :: cot_calipso_mean, cot_calipso_std, cot_calipso_median, cot_calipso_max, cot_calipso_min
      real(double) :: cth_calipso_mean, cgd_calipso_mean
      real(double), dimension(2) :: rdn
      real(double), dimension(npar_mie) :: aerosol_pars_mode
      real(double) :: sum1
      real(double), dimension(19, 7) :: height_input_array, m_height, r_height
      real(double), dimension(19) :: psyn, altdis_dummy
      logical :: Trunc_Flag
      real(double) :: csca, cabs, rlambda, f

      !*** Dimensions
      nwin = size(win_ini)
      nrt = atm_rt%n
      if (scatflag == 2 .or. scatflag == 5 .or. scatflag == 6) then ! aerosol+cirrus
         ntype_aer = 8
      else if (scatflag == 3) then ! only aerosol
         ntype_aer = 7
      else if (scatflag == 4) then ! only cirrus
         ntype_aer = 1
      end if

      !*** Allocate memory
      if (allocated(aerosol)) deallocate (aerosol)
      allocate (aerosol(ntype_aer))
      do k = 1, ntype_aer
         if (allocated(aerosol(k)%rm)) deallocate (aerosol(k)%rm)
         if (allocated(aerosol(k)%fim)) deallocate (aerosol(k)%fim)
         allocate (aerosol(k)%rm(2*nwin))
         allocate (aerosol(k)%fim(2*nwin))
      end do

      if (scatflag == 2 .or. scatflag == 3 .or. scatflag == 5 .or. scatflag == 6) then
         !*** Open atmosferic input files (*.out):
         !*** For aerosol type k and window i, get column average values (weighted by mass)
         !*** for microphysical properties (ECHAM5):
         !*** aerosol(k)%reff, aerosol(k)%veff, aerosol(k)%rm(i), aerosol(k)%fim(i),
         !*** aerosol(k)%aer_col (= total number of particles)
         open (newunit(io), FILE=trim(infile), action='read', status='old')
         do k = 1, 11
            read (io, *)
         end do
         read (io, *) dummy, dummy, dummy, dummy, dummy, aerosol(1:7)%aer_col, dummy7, &
            reff, veff, rm_input_array, fim_input_array
         read (io, *)
         do j = 1, 19
            read (io, *) psyn(j), dummy, dummy, dummy, dummy, dummy, dummy, dummy, dummy, &
               height_input_array(j, :), m_height(j, :), r_height(j, :)
         end do
         close (io)
         do k = 1, 7
            do n = 1, nwin
               if (win_ini(n)%wave_stop <= 1.d7/12700. .and. win_ini(n)%wave_start >= 1.d7/13700.) then
                  i = 2
               else if (win_ini(n)%wave_stop <= 1.d7/5800. .and. win_ini(n)%wave_start >= 1.d7/6500.) then
                  i = 4
               else if (win_ini(n)%wave_stop <= 1.d7/5000. .and. win_ini(n)%wave_start >= 1.d7/5200.) then
                  i = 6
               else if (win_ini(n)%wave_stop <= 1.d7/4700. .and. win_ini(n)%wave_start >= 1.d7/5200.) then
                  i = 7
               else if (win_ini(n)%wave_stop <= 1.d7/4100. .and. win_ini(n)%wave_start >= 1.d7/4500.) then
                  i = 8
               else
                  call stopretrieval('CALCULATE_SYN_SPECTRUM: &
                      &synthetic aerosol input not defined for this window')
               end if
               aerosol(k)%rm(2*n - 1:2*n) = rm_input_array(i + (k - 1)*8) ! Refractive indexrindex constant within a spectral band
               aerosol(k)%fim(2*n - 1:2*n) = fim_input_array(i + (k - 1)*8) ! Note that in OpticM_module the absolute value is taken
            end do
         end do
      end if

      !*** get rid of faulty negative ECHAM5-HAM input
      do i = 1, 7
         reff(i) = max(0.d0, reff(i))
         veff(i) = max(0.d0, veff(i))
         do j = 1, 19
            r_height(j, i) = max(0.d0, r_height(j, i))
            height_input_array(j, i) = max(0.d0, height_input_array(j, i))
            m_height(j, i) = max(0.d0, m_height(j, i))
            if (r_height(j, i) == 0.d0) then
               height_input_array(j, i) = 0.d0
               m_height(j, i) = 0.d0
            end if
         end do
      end do

      !*** column average weighted by mass !number
      do i = 1, 7
         reff(i) = 0.d0
         do j = 1, 19
            reff(i) = reff(i) + m_height(j, i)*r_height(j, i) ! height_input_array(j, i)*r_height(j, i)
         end do
         if (sum(m_height(:, i)) > 1.d-20) then
            reff(i) = reff(i)/sum(m_height(:, i)) ! sum(height_input_array(:, i))
         else
            reff(i) = 0.d0
         end if
      end do

      !*** Convert ECHAM size parameters to reff, veff
      do i = 1, 7
         sigma = log(veff(i))
         aerosol(i)%veff = exp(sigma*sigma) - 1.d0
         aerosol(i)%reff = reff(i)*(1.d0 + aerosol(i)%veff)**2.5d0
      end do

      !*** Open atmosferic input files (*.out):
      !*** Get optical thickness (ECHAM5, MODIS, MERIS, SCIAMACHY and CALIPSO)
      open (newunit(io), FILE=trim(infile), action='read', status='old')
      do i = 1, 7
         read (io, *)
      end do
      read (io, *) aot_modis_mean, aot_modis_max, aot_modis_min, aot_modis_median
      read (io, *)
      read (io, *) cot_calipso_mean, cot_calipso_std, cot_calipso_median, cot_calipso_max, &
         cot_calipso_min, cth_calipso_mean, cgd_calipso_mean
      !*** convert into meters:
      cth_calipso_mean = 1d3*cth_calipso_mean
      cgd_calipso_mean = 1d3*cgd_calipso_mean
      read (io, *)
      read (io, *) aot_echam
      close (io)
      !*** If no MODIS data are available we have to rely on ECHAM5 values. Modify them to get a typical average AOT:
      if (aot_echam .gt. 0.2d0) then
         aot_echam = aot_echam*0.5D0
      end if
      if (aot_echam .lt. 0.1d0) then
         aot_echam = aot_echam*3.0D0
      end if

      !*** Failsafe checks: aot_echam must be larger than 0.0d0
      if (aot_echam .le. 0.0d0) aot_echam = 1.D-10
      if (aot_modis_mean .le. 0.0d0) aot_modis_mean = aot_echam
      if (aot_modis_median .le. 0.0d0) aot_modis_median = aot_echam
      if (aot_modis_max .le. 0.0d0) aot_modis_max = aot_echam
      if (aot_modis_min .le. 0.0d0) aot_modis_min = aot_echam
      if (cot_calipso_mean .le. 0.0d0) cot_calipso_mean = 1.D-10
      if (cot_calipso_median .le. 0.0d0) cot_calipso_median = 1.D-10
      if (cot_calipso_min .le. 0.0d0) cot_calipso_max = 1.D-10
      if (cot_calipso_max .le. 0.0d0) cot_calipso_min = 1.D-10
      if (cot_calipso_std .le. 0.0d0) cot_calipso_std = 1.D0
      if (cth_calipso_mean .le. 0.0d0) cth_calipso_mean = 0.D0
      if (cgd_calipso_mean .le. 0.0d0) cgd_calipso_mean = 1.D0

      if (scatflag == 2 .or. scatflag == 3 .or. scatflag == 5 .or. scatflag == 6) then
         !*** aerosol
         do k = 1, 7
            aerosol(k)%CirrusFlag = 0
            aerosol(k)%id = 1           ! size distribution flag (lognormal)
            aerosol(k)%shapefrac = 1.d0 !fraction of spherical particles (for an aerosol mixture of spheroids and spheres)
            aerosol(k)%tau_ref = aot_modis_median
         end do
      end if

      if (scatflag == 2 .or. scatflag == 4 .or. scatflag == 5 .or. scatflag == 6) then
         !*** cirrus
         aerosol(ntype_aer)%CirrusFlag = 1
         aerosol(ntype_aer)%altid = 2
         aerosol(ntype_aer)%aeralt1 = cth_calipso_mean - cgd_calipso_mean/2.
         aerosol(ntype_aer)%aeralt2 = cgd_calipso_mean
         aerosol(ntype_aer)%id = 0
         aerosol(ntype_aer)%tau_ref = cot_calipso_median
         call rdn01(2, rdn)
         !aerosol(ntype_aer)%shapefrac = rdn(1) !0.5d0  ! Cannot be random if we want co compare
         aerosol(ntype_aer)%shapefrac = 0.5d0  ! Cannot be random if we want co compare
         !*** Tentative size distribution, Heymsfield et al., JAS, 1985
         !aerosol(ntype_aer)%reff = rdn(2)*2.3+2.2 !3.35d0 ! Cannot be random if we want co compare
         aerosol(ntype_aer)%reff = 0.5*2.3 + 2.2 !3.35d0 ! Cannot be random if we want co compare
         aerosol(ntype_aer)%tilt_angle = 30.d0 ! tilted_angle !specified in retrieval_syn.ini
      end if

      !*** Now explicitly calculate aerosol/cirrus properties:
      sum1 = 0.d0
      do k = 1, ntype_aer
         !*** optical properties @ O2 A-band
         rlambda = 0.765
         if (aerosol(k)%CirrusFlag == 0) then
            aerosol_pars_mode(1) = aerosol(k)%reff
            aerosol_pars_mode(2) = aerosol(k)%veff
            aerosol_pars_mode(3) = aerosol(k)%rm(1)
            aerosol_pars_mode(4) = aerosol(k)%fim(1)
            !*** NOTE: aer_col is the unscaled total particle number.
            !*** It will be scaled to the correct total particle number by using tau_ref at 765nm, see end of loop
            aerosol_pars_mode(5) = aerosol(k)%aer_col
            aerosol_pars_mode(6) = aerosol(k)%shapefrac

            call Modes_calc_xs( &
               aero_lut, &
               aerosol(k)%id, &
               aerosol_pars_mode(1:6), &
               rlambda, &
               csca, &
               cabs)
            aerosol(k)%sig_ref = csca + cabs
            !*** Aerosol height distribution
            if (allocated(aerosol(k)%alt_dis)) deallocate (aerosol(k)%alt_dis)
            if (allocated(aerosol(k)%dalt_daer1)) deallocate (aerosol(k)%dalt_daer1)
            if (allocated(aerosol(k)%dalt_daer2)) deallocate (aerosol(k)%dalt_daer2)
            allocate (aerosol(k)%alt_dis(nrt), &
                      aerosol(k)%dalt_daer1(nrt), &
                      aerosol(k)%dalt_daer2(nrt))

            !*** Open atmospheric input files (*.out):
            !*** Get pressure and particle number for 19 height layers and interpolate (ECHAM5)
            do j = 1, 19
               if (sum(height_input_array(:, k)) /= 0) then
                  altdis_dummy(j) = height_input_array(j, k)/sum(height_input_array(:, k))
               end if
            end do
            call linterp(DLOG(psyn), altdis_dummy, 19, &
                         DLOG(atm_rt%p), aerosol(k)%alt_dis, nrt)
            if (sum(aerosol(k)%alt_dis) /= 0) then
               aerosol(k)%alt_dis = aerosol(k)%alt_dis/sum(aerosol(k)%alt_dis)
            end if
            aerosol(k)%maxd = count(aerosol(k)%alt_dis(:) > maxval(aerosol(k)%alt_dis)*nder_cut)
            if (allocated(aerosol(k)%nder)) deallocate (aerosol(k)%nder)
            allocate (aerosol(k)%nder(aerosol(k)%maxd))
            j = 0
            do i = 1, nrt
               if (aerosol(k)%alt_dis(i) > maxval(aerosol(k)%alt_dis)*nder_cut) then
                  j = j + 1
                  aerosol(k)%nder(j) = i
               else
                  aerosol(k)%alt_dis(i) = 0.D0
                  aerosol(k)%dalt_daer1(i) = 0.D0
                  aerosol(k)%dalt_daer2(i) = 0.D0
               end if
            end do
            !*** sum1 needed to scale particle number
            sum1 = sum1 + aerosol(k)%aer_col*aerosol(k)%sig_ref
         else if (aerosol(k)%CirrusFlag == 1) then
            !*** COT needs to be scaled using extinction cross-sections without delta-approximation:
            Trunc_Flag = .false.
            call OPTIC_CIRRUS_XS( &
               cirrus_lut, &
               Trunc_Flag, &
               aerosol(k)%tilt_angle, &
               aerosol(k)%shapefrac, &
               aerosol(k)%reff, &
               rlambda, &
               csca, &
               cabs, &
               f)
            aerosol(k)%sig_ref = csca + cabs
            aerosol(k)%aer_col = aerosol(ntype_aer)%tau_ref/aerosol(ntype_aer)%sig_ref
            !*** Cirrus height distribution
            if (allocated(aerosol(k)%alt_dis)) deallocate (aerosol(k)%alt_dis)
            if (allocated(aerosol(k)%dalt_daer1)) deallocate (aerosol(k)%dalt_daer1)
            if (allocated(aerosol(k)%dalt_daer2)) deallocate (aerosol(k)%dalt_daer2)
            allocate (aerosol(k)%alt_dis(nrt), &
                      aerosol(k)%dalt_daer1(nrt), &
                      aerosol(k)%dalt_daer2(nrt))
            call set_altdis(atm_rt, &
                            aerosol(k)%altid, &
                            aerosol(k)%aeralt1, &
                            aerosol(k)%aeralt2, &
                            aerosol(k)%alt_dis, &
                            aerosol(k)%dalt_daer1, &
                            aerosol(k)%dalt_daer2)
            aerosol(k)%maxd = &
               count(aerosol(k)%alt_dis(:) > maxval(aerosol(k)%alt_dis)*nder_cut)
            if (allocated(aerosol(k)%nder)) deallocate (aerosol(k)%nder)
            allocate (aerosol(k)%nder(aerosol(k)%maxd))
            j = 0
            do i = 1, nrt
               if (aerosol(k)%alt_dis(i) > maxval(aerosol(k)%alt_dis)*nder_cut) then
                  j = j + 1
                  aerosol(k)%nder(j) = i
               else
                  aerosol(k)%alt_dis(i) = 0.D0
                  aerosol(k)%dalt_daer1(i) = 0.D0
                  aerosol(k)%dalt_daer2(i) = 0.D0
               end if
            end do
         end if ! CirrusFlag
      end do !k = 1, ntype_aer

      !*** Scale the total aerosol particle number
      do k = 1, ntype_aer
         if (aerosol(k)%CirrusFlag == 0) then
            aerosol(k)%aer_col = aerosol(k)%aer_col*aerosol(k)%tau_ref/sum1
            aerosol(k)%tau_ref = aerosol(k)%aer_col*aerosol(k)%sig_ref
         end if
      end do
   end subroutine get_aerosol_properties_model

!-------------------------------------------------------------------------------
!> Read settings for creating synthetic measurements
!-------------------------------------------------------------------------------
   subroutine get_geometry_sunsync(time, lat, altitude, szaflag, sza, iza, iaz, saz, phi, ierr)
      integer, intent(in) :: time(7)
      real(double), intent(in) :: lat, altitude!, lon
      integer, intent(in) :: szaflag
      !*** input/output
      real(double), intent(inout) :: sza, iza, iaz, saz, phi
      !*** local
      integer :: io, ierr, year, month, day, hour, min, sec
      real(double) :: sunpos

      year = time(1)
      month = time(2)
      day = time(3)
      hour = time(4)
      min = time(5)
      sec = time(6)

      if (szaflag == 2) then
         !*** Get geometry from sun-synchronous orbit calculation for given time and location
         call orbit_sim_sunsync(year, month, day, hour, min, sec, lat, altitude, sza, ierr)
         if (ierr .ne. 0) then
            ierr = 1
            call writelog('CALCULATE_SYN_SPECTRUM: no sun-synchronous orbit for given latitude', 8)
            return
         end if
         phi = abs(iaz - saz) ! Definition of azimuth angle for RemoTeC
         !*** Quick estimate of geometry
      else if (szaflag == 3) then
         !*** Quick estimate for solar zenith angle:
         sunpos = 23.439*SIN(2*Pi/365*(month*30 + day - 113.))
         sza = abs(lat - sunpos)
         phi = abs(iaz - saz)
      end if

      !*** If sza .ne. 2 or 3, it is assumed that the desired SZA has already been provided though the synsettings-file!
   end subroutine get_geometry_sunsync

!------------------------------------------------------------------------------

   subroutine co2image_noise_model(win_ini, synsettings, earth_radiance, solar_irradiance, flag_out, noise_levels)
      !*** input
      type(window_ini), intent(in) :: win_ini
      type(syn_instrument_noise_settings), intent(in) :: synsettings
      real(double), dimension(:), intent(in) :: earth_radiance    ! Backscattered Earth radiance
      real(double), dimension(:), intent(in) :: solar_irradiance  ! Solar irradiance
      integer, intent(in) :: flag_out

      !*** output
      type(instrument_noise_levels), intent(out) :: noise_levels

      !*** local
      real(double) :: F_num, eta, q_e, t_int, pixel_pitch, dc, binning_factor
      real(double) :: wave_mid, fwhm, lambda_lower, lambda_upper, band_width, dispersion, signal_etendue
      integer :: i, npix_det, io
      real(double) :: szamin, albmax
real(double), dimension(size(earth_radiance)) :: E_radiance, S_irradiance, Emax_radiance, Ne_signal, Ne_signal_max, signal_shotnoise
      integer :: nlambda
      real(double) :: T_bg, lambda_det_min, lambda_det_max, lambda_bg, dlambda
      real(double), dimension(:), allocatable :: BB_lambda_cryo_full_spectrum, L_lambda_cryo_full_spectrum
      real(double) :: detector_etendue, Integrated_flux_cryo_numerical, Ne_bg, bg_shotnoise
      real(double) :: Ne_dark, dark_shotnoise
      real(double) :: quantization_noise
      real(double) :: ROIC_noise, well_capacity, N_readout
      real(double), dimension(size(earth_radiance)) :: SNR

      !*** constants
      real(double), parameter :: hbar = 1.055d-34 ! [Js]
      real(double), parameter :: kB = 1.38d-23    ! [J/K]
      real(double), parameter :: c = 2.998d8      ! [m/s]
      real(double), parameter :: eV = 1.602E-19; !

      fwhm = win_ini%fwhm                                        ! spectral resolution [nm]

      lambda_lower = win_ini%wave_start                          ! lower end wavelength [nm]
      lambda_upper = win_ini%wave_stop                           ! upper end wavelength [nm]
      band_width = (lambda_upper - lambda_lower)               ! spectral band width [nm]
      npix_det = int(band_width/fwhm*win_ini%samp)        ! Number of detector pixels in the spectral (N/S) dimension [pix]
      dispersion = band_width/npix_det                       ! Wavelength interval covered by *one* detector pixel (spectral dimension) [nm/pix]

      !*** Convert Reflected solar spectral radiance (earth_radiance) and
      !*** incoming solar spectral irradiance at TOA (solar_irradiance) to
      !*** SI units (from 'cm-2' to 'm-2')
      E_radiance(:) = earth_radiance(:)*1.d4          ! Reflected solar radiation in spectral radiances [photons/s/m^2/sr/nm]
      S_irradiance(:) = solar_irradiance(:)*1.d4      ! Incoming solar radiation at TOA in spectral irradiances photons/s/m^2/nm]

      !*** Calculate the maximum spectral radiance that is expected to hit the telescope.
      !*** This is given by the incoming solar irradiance at TOA corrected for the smallest
      !*** expected SZA and the maximum expected surface albedo i.e. ignoring any interaction
      !*** with the atmosphere. This controls how often data have to be read out in order to
      !*** avoid pixel saturation.
      szamin = 0.
      albmax = 0.7
      Emax_radiance(:) = S_irradiance*albmax*cos(szamin*pi/180.)/pi     ! [photons/s/cm^2/sr/nm] -  Maximum spectral radiance expected to hit the telescope

      !*** INSTRUMENT PARAMETERS

      !*** 1) Optics
      f_num = synsettings%f_number ! Ratio between telescope diameter and focal length [-]
      eta = synsettings%optics_efficiency ! [-]

      !*** 2) Detector
      pixel_pitch = synsettings%pixel_pitch ! Width of (quadratic) detector pixel [m] (detector pixel area = pixel_pitch**2)
      Q_E = synsettings%quantum_efficiency ! Quantum efficiency: [e-/photons]

      !*** 3) Sampling settings
      t_int = synsettings%integration_time ! [s] Integration/exposure time
      binning_factor = synsettings%binning_factor ! [-] Binning factor for detectors spatial dimension (East-West).
      ! *** TDB: BINNING FOR SPATIAL NORTH-SOUTH DIMENSION IS CONTROLLED THROUGH THE INTEGRATION TIME? ***

      !*** 4) Signal etendue
      ! Note: To obtain the flux, we need to integrate the spectral radiance over
      !       the entendue. The etendue can be calculated in terms of either
      !           i) Solid angle and area of a _single detector pixel_, where the solid angle can be espressed in terms of f-number, or
      !          ii) Solid angle and area of the _telescope_, where the solid angle is derived from the ground sampling area and orbit altitude.
      !       See the lines below for details:
      signal_etendue = pixel_pitch**2*pi/4*1/f_num**2 ! Etendue of a single detector pixel [sr m^2]
      !       = A_detpix * pi/4 * d_telescope**2 / focal_length**2
      !       = A_detpix * A_telescope / focal_length**2
      !       = -> [magnification formula] ->
      !       = A_telescope * ground_sampling_area / altitude**2

      !*** CALCULATE ACTUAL AND EXPECTED MAXIMUM SIGNAL CHARGE FROM ONE GROUND PIXEL ON ONE DETECTOR PIXEL
      Ne_signal(:) = E_radiance(:)*signal_etendue*eta*Q_E*dispersion*t_int*binning_factor ! Number of electron charges per *effective* detector pixel generated by signal from one ground pixel [e-/pix]
      Ne_signal_max(:) = Emax_radiance(:)*signal_etendue*eta*Q_E*dispersion*t_int*binning_factor ! Number of electron charges per *effective* detector pixel generated by expected maximum signal from one ground pixel [e-/pix]

      !*** CALCULATE PER DETECTOR PIXEL ELECTRON CHARGES FROM DIFFERENT NOISE SOURCES

      !*** 1) Signal photon shot noise from ground pixel
      signal_shotnoise = sqrt(Ne_signal(:))

      !*** 2) Signal photon shot noise from background thrmal radiation
      !***    Reference: emission from black-body sources at detector temperature
      !***    TBD: Background radiation/noise from optics itself? Carsten (DLR-OS):
      !***         "There will be some thermal background radiation from the optics
      !***         itself which could be larger than the other background if the
      !***         temperature of the optics is too high"

      !*** 2a) Calculate integrated background radiance from background
      T_bg = synsettings%background_temperature ! Operating temperature of the detector [K]
      lambda_det_min = synsettings%spectral_range(1) ! Lower cutoff wavelength of detector [m]
      lambda_det_max = synsettings%spectral_range(2) ! Upper cutoff wavelength of detector [m]

      nlambda = 1000
      dlambda = (lambda_det_max - lambda_det_min)/(nlambda - 1)

      if (allocated(BB_lambda_cryo_full_spectrum)) deallocate (BB_lambda_cryo_full_spectrum)
      if (allocated(L_lambda_cryo_full_spectrum)) deallocate (L_lambda_cryo_full_spectrum)
      allocate (BB_lambda_cryo_full_spectrum(nlambda))
      allocate (L_lambda_cryo_full_spectrum(nlambda))

      do i = 1, nlambda
         lambda_bg = lambda_det_min + (i - 1)*dlambda
         BB_lambda_cryo_full_spectrum(i) = hbar*pi*(2.*c)**2./lambda_bg**5./(dexp(hbar*2*pi*c/(lambda_bg*kB*T_bg)) - 1)*pi ! BB spectral irradiance [W/m^2/m]
         L_lambda_cryo_full_spectrum(i) = BB_lambda_cryo_full_spectrum(i)/(hbar*2*pi*c/lambda_bg)/pi                         ! Spectral radiance in SI units [photons/s/m^2/sr/m]
      end do

      Integrated_flux_cryo_numerical = sum(L_lambda_cryo_full_spectrum)*dlambda ! Radiance from background incident on detector [photons/s/m^2/sr]

      !*** 2b) Calculate background signal charge
      !        The solid angle of a detector pixel in terms background signal is given by a full hemisphere (=2*pi)
      !        Hence, the detector etendue in terms of background signal is given by A_detpix * pi
      detector_etendue = pixel_pitch**2*pi   ! etendue seen by detector w.r.t. background signal [sr m^2]
      Ne_bg = Integrated_flux_cryo_numerical*detector_etendue*Q_E*t_int*binning_factor ! Number of electron charges generated per *effective* detector pixel through background radiation [e-/pix]

      !*** 2c) Calulate background photon shot noise
      bg_shotnoise = sqrt(Ne_bg)

      !*** 3) Dark current shot noise

      !*** 3a) Determine number of electrons generatured per detector pixel due to dark current - Detector dependent
      dc = synsettings%dark_current         ! Per-pixel dark current @ 150K (?) for proposed detector [e-/s/pix]
      Ne_dark = dc*1e-15*t_int*1./eV*binning_factor ! Number of electrons generated per *effective* detector pixel by dark current [e-/pix]

      !*** 3b) Calulate dark current shot noise
      dark_shotnoise = sqrt(Ne_dark)

      !*** 4) Readout (integrated circuit) noise

      !*** 4a) Determine number of readout electrons generated per pixel for the given
      !***     exposure time and the detectors maximal well fill capacity
      !*** Preliminary guestimate from *CHROMA 1280 data sheet* with full well = 0.7M
      ROIC_noise = 80.           ! [# electrons/pix] ROIC Noise
      well_capacity = 0.7d6      ! [# electrons/pix] maximal well fill capacity of single pixel
      N_readout = int(maxval(Ne_signal_max(:)/(binning_factor*well_capacity))) + 1 ! [dimless] number of read-outs required to avoid detector overflow
      !***TBD: IMPACT OF BINNING ON N_READOUT AND WELL_CAPACITY***
      !***TBD: DON'T I HAVE TO ADD THE NOISE AS WELL WHEN CALCULATING TOTAL SIGNAL WRT TO FULL WELL/LIMIT FOR READOUT
      ROIC_noise = sqrt(N_readout)*ROIC_noise

      !*** 4b) Preliminary estimate from DLR-OS
      ROIC_noise = synsettings%readout_noise

      !*** 5) Quantization noise
      quantization_noise = synsettings%quantization_noise

      !*** CALCULATE SIGNAL-TO-NOISE RATIO AND RELATIVE NOISE
    SNR(:) = Ne_signal(:) / sqrt(signal_shotnoise(:)**2 + bg_shotnoise**2 + dark_shotnoise**2 + binning_factor*ROIC_noise**2 + quantization_noise**2)

      noise_levels%SNR = SNR
      noise_levels%signal_shotnoise = signal_shotnoise
      noise_levels%back_shotnoise = bg_shotnoise
      noise_levels%dark_shotnoise = dark_shotnoise
      noise_levels%readout_noise = ROIC_noise
      noise_levels%quantization_noise = quantization_noise

      if (flag_out >= 2) then
         print *, '... Radiance (min,max):...........', minval(E_radiance), maxval(E_radiance)
         print *, '... Signal (min,max):.............', minval(Ne_signal), maxval(Ne_signal)
         print *, '... SNR (min,max):................', minval(SNR), maxval(SNR)
         print *, '... signal_shotnoise (min,max):...', minval(signal_shotnoise), maxval(signal_shotnoise)
         print *, '... background_shotnoise:.........', bg_shotnoise
         print *, '... dark_shotnoise:...............', dark_shotnoise
         print *, '... ROIC_noise:...................', ROIC_noise
         print *, '... quantization_noise:...........', quantization_noise
         print *, ' '
         print *, '... Maximum signal [e-]:..........', maxval(Ne_signal_max)
         print *, '... Background signal [e-]:.......', Ne_bg
         print *, '... Dark current signal [e-]:.....', Ne_dark
         print *, '... '
      end if
   end subroutine co2image_noise_model

!-------------------------------------------------------------------------------

   subroutine get_spectrum_hi()
   end subroutine get_spectrum_hi

!-------------------------------------------------------------------------------

   subroutine read_residual_image( &
      wavelength_lo, &
      nwave_lo, &
      radiance_lo, &
      residual_image_path, &
      residual_image_intensity, &
      residual_image_lo)
      use netcdf
      !*** input
      real(double), dimension(:), intent(in) :: wavelength_lo
      integer, intent(in) :: nwave_lo
      real(double), dimension(:), intent(in) :: radiance_lo
      character(stringlen), intent(in) :: residual_image_path
      real(double), intent(in) :: residual_image_intensity
      !*** output
      real(double), dimension(:), allocatable, intent(out) :: residual_image_lo
      !*** local
      integer :: k, ierr, ncid, varid
      real(double), dimension(:), allocatable :: wavelength_read
      real(double), dimension(:), allocatable :: radiance_read

      integer :: nDimensions, nVariables, nAttributes, unlimitedDimId, formatNum

      !*** open netcdf file and read wavelength and radiance for residual image
      call check(NF90_OPEN(trim(residual_image_path), nf90_nowrite, ncid), ierr)

      call check(NF90_INQ_VARID(ncid, 'wavelength', varid), ierr)
      allocate (wavelength_read(nwave_lo))
      call check(NF90_GET_VAR(ncid, varid, wavelength_read, [1]), ierr)

      call check(NF90_INQ_VARID(ncid, 'radiance', varid), ierr)
      allocate (radiance_read(nwave_lo))
      call check(NF90_GET_VAR(ncid, varid, radiance_read, [1]), ierr)

      call check(NF90_CLOSE(ncid), ierr)

      !*** check if wavelengths are correct
      if (DABS(wavelength_lo(1) - wavelength_read(1)) > 0.0001 &
          .or. DABS(wavelength_lo(nwave_lo) - wavelength_read(nwave_lo)) > 0.0001) then
         print *, 'READ RESIDUAL IMAGE: residual image not defined for this wavelength'
         return
      end if

      if (allocated(residual_image_lo)) deallocate (residual_image_lo)
      allocate (residual_image_lo(nwave_lo))

      residual_image_lo = radiance_read
      residual_image_lo = residual_image_lo/maxval(radiance_read)*maxval(radiance_lo)
      residual_image_lo = residual_image_lo*residual_image_intensity
   end subroutine read_residual_image

!-------------------------------------------------------------------------------

   subroutine read_signal_intensity_transfer_function( &
      sitf_path, sitf_signal, sitf_response)
      use netcdf
      !*** input
      character(stringlen), intent(in) :: sitf_path
      !*** output
      real(double), dimension(:), allocatable, intent(out) :: sitf_signal
      real(double), dimension(:), allocatable, intent(out) :: sitf_response
      !*** local
      integer :: ncid, varid, ierr, dimid
      integer :: n_sitf ! number of elements in sitf array

      !*** open netcdf file
      call check(NF90_OPEN(trim(sitf_path), nf90_nowrite, ncid), ierr)

      !*** first we need dimension so we know how long signal and response are
      !*** going to be
      call check(NF90_INQ_DIMID(ncid, 'n_sitf', dimid), ierr)
      call check(NF90_INQUIRE_DIMENSION(ncid, dimid, len=n_sitf), ierr)

      !*** read signal and response
      call check(NF90_INQ_VARID(ncid, 'signal', varid), ierr)
      allocate (sitf_signal(n_sitf))
      call check(NF90_GET_VAR(ncid, varid, sitf_signal, [1]), ierr)

      call check(NF90_INQ_VARID(ncid, 'response', varid), ierr)
      allocate (sitf_response(n_sitf))
      call check(NF90_GET_VAR(ncid, varid, sitf_response, [1]), ierr)

      !*** close netcdf file
      call check(NF90_CLOSE(ncid), ierr)
   end subroutine read_signal_intensity_transfer_function

!-------------------------------------------------------------------------------

   subroutine apply_signal_intensity_transfer_function( &
      radiance_lo, sitf_signal, sitf_response)
      !*** input
      real(double), dimension(:), intent(in) :: sitf_signal
      real(double), dimension(:), intent(in) :: sitf_response
      !*** input/output
      real(double), dimension(:), intent(inout) :: radiance_lo
      !*** local
      integer :: n_wave, n_sitf
      integer :: i_wave, i_sitf
      real(double) :: new_radiance

      n_wave = SIZE(radiance_lo)
      n_sitf = SIZE(sitf_signal)

      do i_wave = 1, n_wave ! for each element in the radiance array
         ! optimization: sort radiance_lo and return indices to get it back into
         ! order again after applying sitf. in that case you only have to step
         ! through the sitf array once because you know that the radiance_lo
         ! array is monotoneously increasing.
         do i_sitf = 1, n_sitf
            if (radiance_lo(i_wave) .lt. sitf_signal(i_sitf)) then
               EXIT
            end if
         end do

         call linear_interpolation( &
            radiance_lo(i_wave), &
            new_radiance, &
            sitf_signal(i_sitf - 1), &
            sitf_signal(i_sitf), &
            sitf_response(i_sitf - 1), &
            sitf_response(i_sitf))

         radiance_lo(i_wave) = new_radiance
      end do

   end subroutine apply_signal_intensity_transfer_function

!-------------------------------------------------------------------------------

   subroutine linear_interpolation(x, y, x0, x1, y0, y1)
      real(double), intent(in) :: x, x0, x1, y0, y1
      real(double), intent(out) :: y

      y = (y0*(x1 - x) + y1*(x - x0))/(x1 - x0)
   end subroutine

!-------------------------------------------------------------------------------

end module calculate_syn_spectrum_module
