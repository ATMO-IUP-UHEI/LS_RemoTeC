module read_synsettings_module
   use header_module
   implicit none

   !> @instrument specifications
   type :: syn_instrument_noise_settings

      ! Flags indicating the type of noise type/model
      integer, dimension(:), allocatable :: noiseflag_rad
      integer, dimension(:), allocatable :: noiseflag_irrad

      ! Reference SNR, albedo and SZA following a simple shotnoise model
      real(double), dimension(:), allocatable :: snr_rad
      real(double), dimension(:), allocatable :: snr_irrad
      real(double) :: ref_alb
      real(double) :: ref_sza

      ! Instrument design parameters used to compute corresponding instrument SNR for CO2IMAGE
      real(double) :: f_number
      real(double) :: optics_efficiency
      real(double) :: pixel_pitch
      real(double) :: quantum_efficiency
      real(double) :: spectral_range(2)
      real(double) :: background_temperature
      real(double) :: dark_current
      real(double) :: readout_noise
      real(double) :: quantization_noise
      real(double) :: integration_time
      real(double) :: binning_factor

      ! Constants for estimating the instrument SNR of CO2M (using the parameterization: SNR = a*L/sqrt(a*L + b**2))
      real(double), dimension(:), allocatable :: snr_rad_a
      real(double), dimension(:), allocatable :: snr_rad_b

      ! Straylight offset. This is a constant offset relative to the maximum intensity of the low resolution spectrum
      integer :: straylight_flag
      real(double) :: straylight_offset

      ! Residual image. This is a spectrum with an intensity relative to the maximum intensity of the low resolution spectrum, that
      ! is added to it.
      integer :: residual_image_flag
      character(stringlen) :: residual_image_path
      real(double) :: residual_image_intensity

      ! Signal intensity transfer function. Detector might have nonlinear response to incoming radiation. A SITF is read which
      ! is used to modify the low resolution spectrum. The intensity of each wavelength is input into the SITF and a new modified
      ! spectrum is output.
      integer :: sitf_flag
      character(stringlen) :: sitf_path

   end type syn_instrument_noise_settings

contains

   subroutine read_synsettings(synsettings_file, atmflag, outputflag, nwin, meta, synsettings, xco2, ierr)
      !*** Input
      integer, intent(in) :: atmflag, outputflag, nwin
      character(len=*), intent(in) :: synsettings_file
      type(metadata), intent(inout) :: meta
      !*** Output
      type(syn_instrument_noise_settings), intent(out) :: synsettings
      real(double), intent(out) :: xco2
      integer, intent(out) :: ierr

      !*** local
      character(stringlen) :: message
      integer, parameter :: nwin_max = 10
      integer :: ierr2, io
      integer :: year, month, day, hour, min, sec
      integer, dimension(nwin_max) :: flag_rad, flag_irrad
      real(double), dimension(nwin_max) :: snr_rad, snr_irrad
      real(double), dimension(nwin_max) :: snr_rad_a, snr_rad_b
      real(double) :: ref_alb, ref_sza
      real(double) :: f_number, optics_efficiency, pixel_pitch, quantum_efficiency, spectral_range(2), background_temperature, dark_current, readout_noise, quantization_noise, binning_factor, integration_time
      integer :: szaflag
      real(double) :: sza, iza, iaz, saz, phi, altitude
      integer :: straylight_flag
      real(double) :: straylight_offset
      integer :: residual_image_flag
      character(stringlen) :: residual_image_path
      real(double) :: residual_image_intensity
      integer :: sitf_flag
      character(stringlen) :: sitf_path

      namelist /meta_time/ year, month, day, hour, min, sec
      namelist /xatm/ xco2
      namelist /noise_flag/ flag_rad, flag_irrad
      namelist /shotnoise_settings/ snr_rad, snr_irrad, ref_alb, ref_sza
      namelist /co2image_noise_settings/ f_number, optics_efficiency, pixel_pitch, quantum_efficiency, spectral_range, background_temperature, dark_current, readout_noise, quantization_noise, binning_factor, integration_time
      namelist /co2m_noise_settings/ snr_rad_a, snr_rad_b
      namelist /satellite/ szaflag, sza, iza, iaz, saz, altitude
      namelist /straylight/ straylight_flag, straylight_offset
      namelist /residual_image/ residual_image_flag, residual_image_path, residual_image_intensity
      namelist /signal_intensity_transfer_function/ sitf_flag, sitf_path

      !*** Default values
      year = 0
      month = 0
      day = 0
      hour = 0
      min = 0
      sec = 0

      xco2 = 0.d0

      flag_rad = -1
      flag_irrad = -1

      snr_rad = 0.d0
      snr_irrad = 0.d0
      ref_alb = 0.d0
      ref_sza = 0.d0

      f_number = 0.d0
      optics_efficiency = 0.d0
      pixel_pitch = 0.d0
      quantum_efficiency = 0.d0
      spectral_range = 0.d0
      background_temperature = 0.d0
      dark_current = 0.d0
      readout_noise = 0.d0
      quantization_noise = 0.d0
      binning_factor = 0.d0
      integration_time = 0.d0

      snr_rad_a = 0.d0
      snr_rad_b = 0.d0

      szaflag = 0
      sza = -999.
      iza = 0.
      iaz = 0.
      saz = 0.
      phi = 0.
      altitude = 0.0

      straylight_flag = 0
      straylight_offset = 0.d0

      residual_image_flag = 0
      residual_image_path = 'no input'
      residual_image_intensity = 0.d0

      sitf_flag = 0
      sitf_path = 'no input'

      !*** open namelist file
      open (newunit(io), file=trim(synsettings_file), status='OLD', iostat=ierr)
      if (ierr == 0) then
         !*** read namelists
         read (io, nml=meta_time, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading metatime data')
         rewind (io)

         if (atmflag == 7) then
            read (io, nml=xatm, iostat=ierr2)
            if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading xatm data')
            rewind (io)
         end if

         read (io, nml=satellite, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading satellite data')
         rewind (io)

         read (io, nml=straylight, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading straylight data')
         rewind (io)

         read (io, nml=residual_image, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading residual image data')
         rewind (io)

         read (io, nml=signal_intensity_transfer_function, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading signal intensity transfer function data')
         rewind (io)

         read (io, nml=noise_flag, iostat=ierr2)
         if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading noise-flag data')
         rewind (io)

         if (ANY(flag_rad == 0) .or. ANY(flag_irrad == 0)) then
            read (io, nml=shotnoise_settings, iostat=ierr2)
            if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading shotnoise-settings data')
            rewind (io)
         end if

         if (ANY(flag_rad == 1) .or. ANY(flag_irrad == 1)) then
            read (io, nml=co2image_noise_settings, iostat=ierr2)
            if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading CO2IMAGE instrument noise settings data')
            rewind (io)
         end if

         if (ANY(flag_rad == 2) .or. ANY(flag_irrad == 2)) then
            read (io, nml=co2m_noise_settings, iostat=ierr2)
            if (ierr2 .ne. 0) call stopretrieval('READ_SYNSETTINGS: error reading CO2M instrument noise settings data')
            rewind (io)
         end if

         close (io)
      else
         call stopretrieval('READ_SYNSETTINGS: error reading synsettings info')
      end if

      meta%time(1) = year
      meta%time(2) = month
      meta%time(3) = day
      meta%time(4) = hour
      meta%time(5) = min
      meta%time(6) = sec
      meta%time(7) = 0.d0

      meta%szaflag = szaflag
      meta%sza = sza
      meta%iza = iza
      meta%saz = saz
      meta%iaz = iaz
      meta%phi = phi
      meta%altitude = altitude

      synsettings%noiseflag_rad = flag_rad(1:nwin)
      synsettings%noiseflag_irrad = flag_irrad(1:nwin)

      synsettings%snr_rad = snr_rad(1:nwin)
      synsettings%snr_irrad = snr_irrad(1:nwin)
      synsettings%ref_alb = ref_alb
      synsettings%ref_sza = ref_sza

      synsettings%f_number = f_number
      synsettings%optics_efficiency = optics_efficiency
      synsettings%pixel_pitch = pixel_pitch*1e-6
      synsettings%quantum_efficiency = quantum_efficiency
      synsettings%spectral_range = spectral_range*1e-9
      synsettings%background_temperature = background_temperature
      synsettings%dark_current = dark_current
      synsettings%readout_noise = readout_noise
      synsettings%quantization_noise = quantization_noise
      synsettings%binning_factor = binning_factor
      synsettings%integration_time = integration_time*1e-3

      synsettings%snr_rad_a = snr_rad_a(1:nwin)
      synsettings%snr_rad_b = snr_rad_b(1:nwin)

      synsettings%straylight_flag = straylight_flag
      synsettings%straylight_offset = straylight_offset

      synsettings%residual_image_flag = residual_image_flag
      synsettings%residual_image_path = residual_image_path
      synsettings%residual_image_intensity = residual_image_intensity

      synsettings%sitf_flag = sitf_flag
      synsettings%sitf_path = sitf_path

      if (outputflag >= 2) then
         write (6, nml=meta_time)
         write (6, nml=satellite)
         write (6, nml=noise_flag)

         if (ANY(flag_rad == 0) .or. ANY(flag_irrad == 0)) then
            write (6, nml=shotnoise_settings)
         end if

         if (ANY(flag_rad == 1) .or. ANY(flag_irrad == 1)) then
            write (6, nml=co2image_noise_settings)
         end if

         if (ANY(flag_rad == 2) .or. ANY(flag_irrad == 2)) then
            write (6, nml=co2m_noise_settings)
         end if

      end if

      if (outputflag > 1) then
         write (message, '(a)') '*** End of READ_SYNSETTINGS ***'
         call writelog(message, 1)
      end if

      ierr = 0 ! Mark success
      return

   end subroutine read_synsettings

!------------------------------------------------------------------------------
end module read_synsettings_module
