module wrapper_retrieve_module
   use header_module
   use atmosphere_interface_module, only: atmospheric_scenario, read_atm
   use retrieval_module, only: retrieval_data, retrieval, aero, &
                               Mie_lut, cirrus_table, read_aerosol_netcdf, read_cirrus_netcdf, &
                               window_ini, settings_flags, file_paths, altitude_grid, read_settings, read_win_xsdb
   use synthetic_input_module, only: instrument_errors, read_errors, synthetic_data, get_synthetic_data, get_synthetic_data_nc
  use spectrum_interface_module, only: spectrum, read_l1b, instrument_response, get_isrf_interpolated
   use solar_model_module, only: sun_spectrum, read_sun_netcdf, read_sun_tsis1_hsrs, interpolate_solar_spectrum
   use diagnostics_module, only: diagnostics_retrieve, diagnostics_retrieve_nc_js, diagnostics_retrieve_nc_ls

   implicit none
   private
   public :: init_shared, free_shared, init_pixel, free_pixel, retrieve_wrapper, write_output

!------------------------------------------------------------------------------
!> Shared data: same for each ground pixel
!------------------------------------------------------------------------------
   type, public :: shared_data
      type(window_ini), dimension(:), allocatable :: win_ini                  ! window input in setting.nml
      type(Mie_lut) :: aero_lut                                               ! Aerosol LUT
      type(cirrus_table) :: cirrus_lut                                        ! Cirrus LUT
      real(double), dimension(4) :: meteo_errors                              ! assumed errors in meteo data
      type(instrument_errors), dimension(:), allocatable :: instr_errors      ! assumed instrument errors
      type(aero), dimension(:), allocatable :: aerosol_ini                    ! Initial guess for aerosol parameters
      type(settings_flags) :: flag                                            ! retrieval settings fixed_flags
      type(altitude_grid) :: grid                                             ! atmospheric grids for retrieval, RT and XS
      type(file_paths) :: path                                                ! paths to input files
      type(instrument_response), dimension(:, :), allocatable :: response      ! Iinstrument Spectral Response Function
   end type shared_data

!------------------------------------------------------------------------------
!> Thread-specific data: different for each ground pixel
!------------------------------------------------------------------------------
   type, public :: pixel_data
      type(spectrum), dimension(:), allocatable :: measurement
      type(atmospheric_scenario) :: atm_scenario
      type(metadata) :: meta
      character(stringlen) :: spectrum_file, meteo_file, filename
      integer :: ipixel
   end type pixel_data

!------------------------------------------------------------------------------
!> Output data: different for each ground pixel
!------------------------------------------------------------------------------
   type, public :: output_data
      type(retrieval_data) :: retrieval_output
      type(synthetic_data) :: syn_output

   end type output_data

contains
!------------------------------------------------------------------------------
!> Initialize shared input data for retrieval algorithm
!------------------------------------------------------------------------------
   subroutine init_shared(fixed, runpath, first_atm, runid, ierr)
      type(shared_data), pointer, intent(inout) :: fixed
      character(stringlen), intent(in) :: runpath, first_atm
      integer, intent(in) :: runid
      integer, intent(out) :: ierr
      !*** local
      type(spectrum), dimension(:), allocatable :: measurement
      type(metadata) :: meta
      type(sun_spectrum) :: sun_input  ! input reference solar spectrum
      character(stringlen) :: settings_file, spectrum_file
      character(6):: runidstring
      integer :: i, k, n
      integer, dimension(:), allocatable :: nrow

      allocate (fixed)
      !*** Read in retrieval settings
      write (runidstring, '(I6.6)') runid
      settings_file = trim(runpath)//'INI/settings_RTC_retrieve_'//runidstring//'.nml'
      call read_settings( &
         settings_file, &
         fixed%win_ini, &
         fixed%aerosol_ini, &
         fixed%grid, &
         fixed%flag, &
         fixed%path, &
         ierr)
      if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_settings.', 8)

      !*** Read in errors for sensitivity analysis
      allocate (fixed%instr_errors(size(fixed%win_ini)))
      call read_errors(runpath, runid, size(fixed%win_ini), fixed%meteo_errors, fixed%instr_errors)

      !*** Apply ISRF error
      do n = 1, size(fixed%win_ini)
         fixed%win_ini(n)%fwhm = fixed%win_ini(n)%fwhm*(1.d0 + 0.01d0*fixed%instr_errors(n)%isrf)
         !*** Calculate hi-resolution wavelength grid
         fixed%win_ini(n)%nwave_hi = int((fixed%win_ini(n)%wave_stop + fixed%win_ini(n)%fwhm*fixed%win_ini(n)%wvbd - &
                               fixed%win_ini(n)%wave_start + fixed%win_ini(n)%fwhm*fixed%win_ini(n)%wvbd)/fixed%win_ini(n)%reso) + 1
         deallocate (fixed%win_ini(n)%wavelength_hi)
         allocate (fixed%win_ini(n)%wavelength_hi(fixed%win_ini(n)%nwave_hi))

         forall (k=1:fixed%win_ini(n)%nwave_hi) fixed%win_ini(n)%wavelength_hi(k) = &
            dble(int(1./fixed%win_ini(n)%reso*(fixed%win_ini(n)%wave_start - fixed%win_ini(n)%fwhm*fixed%win_ini(n)%wvbd)))* &
            fixed%win_ini(n)%reso + dble(k - 1)*fixed%win_ini(n)%reso
      end do

      !*** Read reference irradiance
      if (fixed%flag%solar == 0) then
         call read_sun_netcdf(fixed%path%sun, sun_input, ierr)
         if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_sun_netcdf.', 8)
      else if (fixed%flag%solar == 1) then
         call read_sun_tsis1_hsrs(fixed%path%sun, sun_input, ierr)
         if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_sun_tsis1_hsrs.', 8)
      end if

      !*** Interpolate irradiance to internal wavelength grid
      call interpolate_solar_spectrum(sun_input, fixed%win_ini, ierr)
      if (ierr .ne. 0) call writelog('INIT_SHARED: Error in interpolate_solar_spectrum.', 8)

      !*** Read in Tables with Mie / T-matrix scattering properties
      do i = 1, size(fixed%aerosol_ini)
         if (fixed%aerosol_ini(i)%Cirrusflag == 0) then
            call read_aerosol_netcdf(fixed%path%mie, fixed%flag%output, fixed%aero_lut, ierr)
            if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_aerosol_netcdf.', 8)
            exit
         end if
      end do

      !*** Read in cirrus tables
      do i = 1, size(fixed%aerosol_ini)
         if (fixed%aerosol_ini(i)%Cirrusflag == 1) then
            call read_cirrus_netcdf(fixed%path%cirrus, fixed%aerosol_ini(i)%tilt_angle, fixed%flag%output, fixed%cirrus_lut, ierr)
            !  call read_cirrus_ascii(cirruspath, aerosol_ini(i)%tilt_angle, cirrus_lut)
            if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_cirrus_netcdf.', 8)
            exit
         end if
      end do

      !*** Read absorption cross sections from database into memory
      call read_win_xsdb(fixed%flag, fixed%win_ini, ierr)
      if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_win_xsdb.', 8)

      !*** Get measured spectral grid
      if (fixed%flag%atm == 4) then
         call read_l1b(trim(fixed%path%spectrum)//'L1B_'//first_atm, fixed%flag%output, measurement, meta, ierr, fixed%flag%synthetic_input, fixed%flag%observer_location, fixed%win_ini)
         if (ierr .ne. 0) call writelog('INIT_SHARED: Error in read_l1b.', 8)
      else
         print*, "ERROR: READING ATMOSPHERE WITH FLAG ", fixed%flag%atm, " NOT SUPPORTED ANYMORE"
      end if

      !*** Get Instrument Spectral Response Function on appropiate spectral grids:
      !*** measured spectral grid (lo-reso) and model spectral grid (hi_reso)
      !*** Note that the measured spectral grid needs to be set through datastructure measurement(nband)
      allocate (fixed%response(1, size(measurement)), stat=ierr)
      allocate (nrow(size(measurement)), stat=ierr)
      nrow = 1 !ISRF is identical for all detector rows
      call get_isrf_interpolated( &
         fixed%flag%ilscalc, &
         trim(fixed%path%ils), &
         nrow, &
         fixed%win_ini, &
         measurement, &
         fixed%response, &
         ierr)
      if (ierr .ne. 0) call writelog('INIT_SHARED: Error in get_isrf_interpolated.', 8)

      return

   end subroutine init_shared

!------------------------------------------------------------------------------
!> Deallocate memory of shared data
!------------------------------------------------------------------------------
   subroutine free_shared(fixed)
      type(shared_data), pointer, intent(inout) :: fixed

      deallocate (fixed)

   end subroutine free_shared

!------------------------------------------------------------------------------
!> Prepare ground-pixel dependent input data for retrieval algorithm
!------------------------------------------------------------------------------
   subroutine init_pixel(fixed, varying, output, atm, ipixel, ierr)
      type(shared_data), pointer, intent(in) :: fixed
      type(pixel_data), pointer :: varying
      type(output_data), pointer :: output
      character(stringlen), intent(in):: atm
      integer, intent(in) :: ipixel
      integer, intent(out) :: ierr
      !*** local
      character(6) :: name

      allocate (varying)
      allocate (output)

      !*** pixel identifier
      varying%ipixel = ipixel

      !*** Is the current pixel over ocean-sunglint?
      varying%meta%oceanglint = 0

      !*** a priori value for fluorescence
      varying%meta%FS = 0.d0

      !*** Read data
      !*** Get atmosphere file name from atm(i)
      varying%filename = trim(atm)

      !*** meteo data
      if (fixed%flag%atm == 4) then
         varying%meteo_file = trim(fixed%path%meteo)//'ATM_'//trim(varying%filename)
         call read_atm(varying%meteo_file, fixed%flag%output, varying%atm_scenario, varying%meta, ierr, fixed%meteo_errors)
         if (ierr .ne. 0) return
      else
         print*, "ERROR: READING ATMOSPHERE WITH FLAG ", fixed%flag%atm, " NOT SUPPORTED ANYMORE"
      end if

      !*** spectrum
      if (fixed%flag%atm == 4) then
         varying%spectrum_file = trim(fixed%path%spectrum)//'L1B_'//trim(varying%filename)
         call read_l1b(varying%spectrum_file, fixed%flag%output, varying%measurement, varying%meta, ierr, fixed%flag%synthetic_input, fixed%flag%observer_location, fixed%win_ini, fixed%instr_errors)
         if (ierr .ne. 0) return
      else
         print*, "ERROR: READING ATMOSPHERE WITH FLAG ", fixed%flag%atm, " NOT SUPPORTED ANYMORE"
      end if

      !*******************************************************************
      !*** Only for synthetic spectra
      !*******************************************************************

      if (fixed%flag%synthetic_input == 0) then
         varying%filename = "ATM_"//varying%filename
      else if (fixed%flag%synthetic_input == 1) then
         if (fixed%flag%atm == 3) then
            write (name, '(I6.6)') ipixel
            varying%filename = 'ATM_'//name
        !    call get_synthetic_data( &
        !         fixed%win_ini, varying%atm_scenario, fixed%meteo_errors, fixed%aerosol_ini, &
        !         fixed%path%spectrum, fixed%path%meteo, varying%filename, varying%measurement, &
        !         fixed%grid, fixed%flag%temp, fixed%flag%fit, &
        !         output%syn_output)
            call get_synthetic_data_nc( &
               fixed%win_ini, varying%atm_scenario, fixed%aerosol_ini, &
               varying%measurement, varying%meta, &
               fixed%grid, fixed%flag%temp, fixed%flag%fit, &
               output%syn_output)
         else
            varying%filename = 'ATM_'//varying%filename
            call get_synthetic_data( &
               fixed%win_ini, varying%atm_scenario, fixed%meteo_errors, fixed%aerosol_ini, &
               fixed%path%spectrum, fixed%path%meteo, varying%filename, varying%measurement, &
               fixed%grid, fixed%flag%temp, fixed%flag%fit, fixed%flag%atm, &
               output%syn_output, ierr)
            if (ierr .ne. 0) return
         end if
      end if

   end subroutine init_pixel

!------------------------------------------------------------------------------
!> Call retrieval algorithm
!------------------------------------------------------------------------------
   subroutine retrieve_wrapper(fixed, varying, output, ierr)
      implicit none
      type(shared_data), pointer, intent(in) :: fixed
      type(pixel_data), pointer, intent(inout) :: varying
      type(output_data), pointer, intent(inout) :: output
      integer, intent(out) :: ierr

      call retrieval(varying%measurement, varying%meta, fixed%response(1, :), varying%atm_scenario, &
                     fixed%grid, &
                     fixed%flag, &
                     fixed%win_ini, fixed%aerosol_ini, &
                     fixed%aero_lut, fixed%cirrus_lut, &
                     output%retrieval_output, ierr)

   end subroutine retrieve_wrapper

!------------------------------------------------------------------------------
!> Deallocate memory of thread-specific data
!------------------------------------------------------------------------------
   subroutine free_pixel(varyingData, outputData)
      implicit none
      type(pixel_data), pointer :: varyingData
      type(output_data), pointer :: outputData

      deallocate (outputData)
      deallocate (varyingData)

   end subroutine free_pixel

!------------------------------------------------------------------------------
!> Write output data to file
!------------------------------------------------------------------------------
   subroutine write_output(runId, fixedData, varyingData, atmflag, outputData)
      integer, intent(in) :: runId
      type(shared_data), intent(in), pointer :: fixedData
      type(pixel_data), intent(in), pointer :: varyingData
      integer, intent(in) :: atmflag
      type(output_data), intent(inout), pointer ::outputData
      character(stringlen) :: meteo_file

      if (atmflag .EQ. 4) then
         meteo_file = trim(fixedData%path%spectrum)//trim(varyingData%filename)
         ! call diagnostics_retrieve_nc_js(fixedData%path%output, runid, &
         !                            meteo_file, &
         !                            fixedData%win_ini, &
         !                            outputData%syn_output, outputData%retrieval_output, varyingData%meta, &
         !                            varyingData%measurement, fixedData%grid%nlay, fixedData%flag%output)
         call diagnostics_retrieve_nc_ls(fixedData%path%output, runid, &
                                    meteo_file, &
                                    fixedData%win_ini, &
                                    fixedData%flag%synthetic_input, outputData%syn_output, &
                                    outputData%retrieval_output, varyingData%meta, &
                                    varyingData%measurement, fixedData%grid%nlay, fixedData%flag%output)
      else
         !    meteo_file = trim(fixedData%path%meteo)//trim(varyingData%filename)
         meteo_file = trim(fixedData%path%spectrum)//trim(varyingData%filename)
         call diagnostics_retrieve( &
            fixedData%path%output, runid, fixedData%flag%output, &
            varyingData%spectrum_file, meteo_file, &
            fixedData%win_ini, &
            outputData%syn_output, outputData%retrieval_output, varyingData%meta, &
            fixedData%grid%nlay, size(varyingData%measurement), fixedData%flag%output)
      end if

   end subroutine write_output
!------------------------------------------------------------------------------
end module wrapper_retrieve_module
