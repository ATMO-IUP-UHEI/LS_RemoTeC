module read_settings_module
   use header_module
   use aerosol_input_module, only: aero
   use read_xsdb_module, only: cross_section_db, read_xsdb, read_xsdb_netcdf, interpolate_xsdb, tri_conv_xsdb, rect_conv_xsdb, read_xsdb_frankenberg, read_xsdb_butz, read_xsdb_absco 
   implicit none
   private

!*** types
   public :: window_ini, settings_flags, file_paths, regularization_class, altitude_grid, filter_thresholds

!*** procedures
   public :: read_settings, read_win_xsdb

!> @brief Window-dependent retrieval settings
   type :: window_ini 
!*** Solar spectrum
      real(double), dimension(:), allocatable :: sun_spectrum_ref_hi   !< Hi-res sun reference (Dim: nwave_hi)
      real(double), dimension(:), allocatable :: wavelength_hi   !< hi-res wavelength grid

      real(double) :: wave_start                        !< Window lower wavelength boundary
      real(double) :: wave_stop                         !< Window upper wavelength boundary
      real(double) :: fwhm                              !< FWHM
      real(double) :: fwhm_lores                        !< FWHM for low resolution test (ABUTZ)
      real(double) :: reso                              !< Hi-res wavelength sampling
      real(double) :: samp                              !< oversampling ratio (only relevant for creating synthetic measurement)	
      real(double) :: wvbd                              !< wavelength boundary offset (in multiples of fwhm) 
      real(double) :: albedo(3)                         !< surface albedo (only relevant for creating synthetic measurement)

!*** K-Binning
      real(double) :: ntau, ntau_2nd_max                !< K-binning number of bins, 1st, 2nd species

!*** Flags
      integer :: albflag                                !< Fit nth order albedo
      integer :: spsh0flag                              !< Fit 0th order spectral shift of reflectance wrt sun
      integer :: spsh1flag                              !< Fit 1st order spectral shift of reflectance wrt sun
      integer :: spsh2flag                              !< Fit 2nd order spectral shift of reflectance wrt sun
      integer :: sunsh0flag                             !< Fit 0th order spectral shift of cross sections wrt sun
      integer :: IoffFlag                               !< Fit additive constant to spectrum in this window
      integer :: Fsflag                                 !< Fit fluorescence emission in this window
		       
!*** wavelength grid of model (high-resolution)
     integer :: nwave_hi                                        !< Number of hi-res wavelengths

!*** Absorber properties
      integer :: ntype                                                !< Number of absorbers in this window 
      integer, dimension(:), allocatable :: type_x_flag               !< Flag for target absorber
      integer, dimension(:), allocatable :: type_x                    !< Hitran index of absorber
      integer, dimension(:), allocatable :: type_xsdb                 !< format of xsdb
      character(stringlen), dimension(:), allocatable :: xsdb_path          !< Path of the database file        

      type(cross_section_db), dimension(:), allocatable :: xsdb       !< Cross section database (Dim: ntype)
      real(double), dimension(:), allocatable :: xs_scaling           !< Scalingfactor for cross-sections

   end type window_ini

   type :: filter_thresholds
      real(double) :: T(11)
      real(double) :: BT(9)
      real(double) :: snr_swir, snr_nir, sza, vza, surface_roughness, elevation_diff, minimum_pixels, spectral_offset, chi2
   end type filter_thresholds

   type :: regularization_class
      real(double) :: weight_target(10)
      real(double) :: weight_aerosol
   end type regularization_class
   
   type :: settings_flags
      integer ::  atm, scat, rtm, inv, XS, xs_preprocess, sun, O2, temp, Fs, ils, ilscalc, observer_location, solar, fit, oceanglint, glintscat, output, coreg
      type(regularization_class) :: reg
   end type settings_flags

   type :: file_paths
      character(stringlen) :: spectrum, meteo, ils, sun, mie, cirrus, output, coreg
   end type file_paths

   type :: altitude_grid
      integer :: flag, nlay, nrt, natm 
   end type altitude_grid


   
  contains 

!------------------------------------------------------------------------------  
!> @brief Read retrieval settings from namelist input file
!> @param[in]  settings_file      filename of namelist inputfile with retrieval settings
!> @param[out] win_ini            datatype for window-dependent settings
!> @param[out] aerosol_ini        datatype for aerosol settings and apriori values
!> @param[out] grid               datatype for atmospheric grids of retrieval, RT, and XS
!> @param[out] flag               datatype for retrieval settings flags
!> @param[out] path               datatype for paths to inout files (LUTs, etc) 
!> @param[out] ierr               error identifier: 0=normal, 1=error opening file, 2=error in namelist, 3=error in settings 
!> @param[out] threshhold         datatype for filter thresholds (cloudfilters, geometry, etc)
!------------------------------------------------------------------------------
    subroutine read_settings( &
         settings_file, &
         win_ini, &
         aerosol_ini, &
         grid, &
         flag, &
         path, &
         ierr, &
         threshold, &
         warning)
      !*** Input
      character(len=*), intent(in) :: settings_file
      !*** Output
      type(window_ini), dimension(:), allocatable, intent(out) :: win_ini
      type(aero), dimension(:), allocatable, intent(out) :: aerosol_ini
      type(altitude_grid), intent(out) :: grid
      type(settings_flags), intent(out) :: flag
      type(file_paths), intent(out) :: path
      type(filter_thresholds), optional, intent(out) :: threshold
      type(filter_thresholds), optional, intent(out) :: warning
      integer, intent(out) :: ierr
      !*** local variables
      !*** Maximum dimensions
      integer, parameter :: nwin_max=10
      integer, parameter :: nabs_max=10
      integer, parameter :: naer_max=10
      integer :: i, k, n, io, nwin, ntype_aer
      character(stringlen) :: message
      !*** window
      real(double) :: wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, albedo(3)
      integer :: ntau, ntau_2nd_max
      integer :: albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      !*** aerosol
      integer :: CirrusFlag, height_dist, size_dist, &
           amount_fit, alt1_fit, alt2_fit, &
           rm_fit, fim_fit, shapefrac_fit, reff_fit, veff_fit
      real(double) :: tau_ref, alt1, alt2, shapefrac, reff, veff, tilt_angle
      !*** namelist dummies
      real(double), dimension(2*nwin_max) :: rm, fim
      integer, dimension(nwin_max) :: ntype_abs
      integer, dimension(nabs_max)  :: type_x, type_x_flag, type_xsdb
      character(stringlen), dimension(nabs_max) :: xsdb_path 
      real(double), dimension(nabs_max) :: xs_scaling

      !*** namelists
      namelist /prefilter/ threshold
      namelist /postfilter/ warning
      namelist /flags/  flag
      namelist /paths/ path
      namelist /alt_grid/  grid
      namelist /dim/ nwin, ntype_aer
      namelist /nabs/ ntype_abs
      namelist /window1/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window2/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window3/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window4/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window5/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window6/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window7/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path,  xs_scaling,&
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window8/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window9/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling,&
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /window10/ wave_start, wave_stop, reso, fwhm, fwhm_lores, samp, wvbd, ntau, ntau_2nd_max, &
           type_x, type_x_flag, type_xsdb, xsdb_path, xs_scaling, &
           albedo, albflag, IoffFlag, spsh0flag, spsh1flag, spsh2flag, sunsh0flag
      namelist /aerosol1/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol2/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol3/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol4/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol5/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol6/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol7/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol8/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol9/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle
      namelist /aerosol10/  CirrusFlag,  amount_fit, tau_ref, height_dist,&
           alt1_fit,  alt1, alt2_fit, alt2, rm_fit, rm, fim_fit, fim, size_dist, shapefrac_fit, shapefrac, &
           reff_fit, reff, veff_fit, veff, tilt_angle

      !*** Default settings (if not specified in namelist file)
      !*** setting flags
      flag%atm = 0
      flag%scat = -999
      flag%rtm = 1
      flag%inv = 0
      flag%output = 0
      flag%ils = 1
      flag%ilscalc = 1
      flag%observer_location = 0
      flag%solar = 0
      flag%fit = 2
      flag%xs_preprocess = 1
      flag%sun = 3
      albflag = -999
      ntype_aer = 0
      !*** Debug read cross section flag
      flag%XS = 0
      !*** fit flags
      flag%O2 = 0
      flag%temp = 0
      flag%Fs = 0
      flag%oceanglint = 1
      flag%glintscat = 0
      IoffFlag = 0
      spsh0flag = 0
      spsh1flag = 0
      spsh2flag = 0
      sunsh0flag = 0
      !*** Paths
      path%spectrum = 'no input'
      path%meteo = 'no input'
      path%ils = 'no input'
      path%sun = 'no input'
      path%mie = 'no input'
      path%cirrus = 'no input'
      path%output = 'no input'
      path%coreg = 'no input'
      xsdb_path = 'no input'
      !*** Fill values
      ntype_abs = -999
      type_x = -999
      type_x_flag = -999
      type_xsdb = -999
      rm =  - 999.
      fim = -999.
      albedo = 0.d0
      !*** Aerosol fit flags
      amount_fit = 0
      alt1_fit = 0 
      alt2_fit = 0
      rm_fit = 0
      fim_fit =  0
      reff_fit = 0
      veff_fit = 0
      shapefrac_fit = 0
      !*** Ad hoc regularization, default for GOSAT en OCO-2
      flag%reg%weight_target = 1.
      flag%reg%weight_target(1) = 500.
      flag%reg%weight_target(2) = 100.
      flag%reg%weight_aerosol = 20.
      

      !*** open namelist file
      open(newunit(io), file=settings_file, status='OLD', iostat=ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         write(message, *) 'READ_SETTINGS: error opening namelist input file: ', trim(settings_file), ', iostat = ', ierr
         goto 999
      endif

      !*** read namelists

      !*** Thresholds for prefiltering
      if (present(threshold)) then
         read(io, nml=prefilter, iostat=ierr)
         if(ierr .ne. 0) then
            ierr = ierr_read
            write(message, *)'READ_SETTINGS: error in namelist /prefilter/'
            goto 999
         endif
         rewind(io)
      endif

      !*** Thresholds for postfiltering
      if (present(warning)) then
         read(io, nml=postfilter, iostat=ierr)
         if(ierr .ne. 0) then
            ierr = ierr_read
            write(message, *)'READ_SETTINGS: error in namelist /postfilter/'
            goto 999
         endif
         rewind(io)
      endif

      !*** Directories of input/output
      read(io, nml=paths, iostat=ierr)
      if(ierr .ne. 0) then 
         ierr = ierr_read
         write(message, *) 'READ_SETTINGS: error in namelist /paths/'
         goto 999
      endif
      rewind(io)

      !*** STANDARD ALTITUDE GRID (X, RT, XS)
      read(io, nml=alt_grid, iostat=ierr)
      if(ierr .ne. 0) then
         ierr = ierr_read
         write(message, *) 'READ_SETTINGS: error in namelist /alt_grid/'
         goto 999
      endif
      rewind(io)

      !*** Retrieval flags
      read(io, nml=flags, iostat=ierr)
      if(ierr .ne. 0) then
         ierr = ierr_read
         write(message, *) 'READ_SETTINGS: error in namelist /flags/'
         goto 999
      endif
      rewind(io)

      !*** Number of windows and aerosol types
      read(io, nml=dim, iostat=ierr)
      if(ierr .ne. 0) then
         ierr = ierr_read
         write(message, *) 'READ_SETTINGS: error in namelist /dim/'
         goto 999
      endif
      rewind(io)

      if(nwin<0 ) call stopretrieval( 'READ_SETTINGS: min #windows = 1 ') 
      if(nwin>nwin_max ) call stopretrieval( 'READ_SETTINGS: max #windows exceeded ') 
      if(ntype_aer>naer_max ) call stopretrieval( 'READ_SETTINGS: max #aerosol exceeded') 

      !*** Number of absorbers in each window
      read(io, nml=nabs, iostat=ierr)
      if(ierr .ne. 0) then
         ierr = ierr_read
         write(message, *) 'READ_SETTINGS: error in namelist /nabs/'
         goto 999
      endif
      rewind(io)

      if (flag%output > 1) then
         !*** write namelist input to screen
         write(message,'(a)')  '*** Start of READ_SETTINGS ***'  
         call writelog(message, 1)
         write(6, nml=paths)      
         write(6, nml=alt_grid)
         write(6, nml=flags)
         write(6, nml=dim)
         write(6, nml=nabs)     
      endif

      !*** Window-dependent retrieval settings
      if(allocated(win_ini)) deallocate(win_ini)
      allocate(win_ini(nwin))
      do n = 1, nwin
         xs_scaling = 1.d0    
         fwhm_lores = 0.d0     
         !*** Read window-dependent retrieval settings     
         if (n == 1) then 
            read(io, nml=window1, iostat=ierr)
            if (flag%output > 1) write(6, nml=window1)   
         elseif (n == 2) then
            read(io, nml=window2, iostat=ierr)
            if (flag%output > 1) write(6, nml=window2) 
         elseif (n == 3) then
            read(io, nml=window3, iostat=ierr)
            if (flag%output > 1) write(6, nml=window3) 
         elseif (n == 4) then
            read(io, nml=window4, iostat=ierr)
            if (flag%output > 1) write(6, nml=window4)
         elseif (n == 5) then
            read(io, nml=window5, iostat=ierr)
            if (flag%output > 1) write(6, nml=window5) 
         elseif (n == 6) then
            read(io, nml=window6, iostat=ierr)
            if (flag%output > 1) write(6, nml=window6) 
         elseif (n == 7) then
            read(io, nml=window7, iostat=ierr)
            if (flag%output > 1) write(6, nml=window7) 
         elseif (n == 8) then
            read(io, nml=window8, iostat=ierr)
            if (flag%output > 1) write(6, nml=window8) 
         elseif (n == 9) then
            read(io, nml=window9, iostat=ierr)
            if (flag%output > 1) write(6, nml=window9) 
         elseif (n == 10) then 
            read(io, nml=window10, iostat=ierr)
            if (flag%output > 1)  write(6, nml=window10) 
         endif

         if(ierr .ne. 0) then
            ierr = ierr_read
            write(message, *) 'READ_SETTINGS: error in namelist /window/'
            goto 999
         endif
         rewind(io)    
         !*** Put in type win_ini
         win_ini(n)%wave_start = wave_start
         win_ini(n)%wave_stop = wave_stop
         win_ini(n)%reso = reso
         win_ini(n)%fwhm = fwhm 
         win_ini(n)%fwhm_lores = fwhm_lores 
         win_ini(n)%samp = samp
         win_ini(n)%wvbd = wvbd
         win_ini(n)%ntype = ntype_abs(n)
         win_ini(n)%ntau = ntau
         win_ini(n)%ntau_2nd_max = ntau_2nd_max
         if(win_ini(n)%ntype==1) then
            win_ini(n)%ntau_2nd_max=1
         endif
         if (allocated(win_ini(n)%type_x)) deallocate(win_ini(n)%type_x)
         if (allocated(win_ini(n)%type_x_flag)) deallocate(win_ini(n)%type_x_flag)
         if (allocated(win_ini(n)%type_xsdb)) deallocate(win_ini(n)%type_xsdb)
         if (allocated(win_ini(n)%xsdb_path)) deallocate(win_ini(n)%xsdb_path)   
         if (allocated(win_ini(n)%xs_scaling)) deallocate(win_ini(n)%xs_scaling) 
         allocate(win_ini(n)%type_x(win_ini(n)%ntype))
         allocate(win_ini(n)%type_x_flag(win_ini(n)%ntype))
         allocate(win_ini(n)%type_xsdb(win_ini(n)%ntype))
         allocate(win_ini(n)%xsdb_path(win_ini(n)%ntype))
         allocate(win_ini(n)%xs_scaling(win_ini(n)%ntype))
         do k = 1, win_ini(n)%ntype
            win_ini(n)%type_x(k) = type_x(k)
            win_ini(n)%type_x_flag(k) = type_x_flag(k)
            win_ini(n)%type_xsdb(k) = type_xsdb(k) 
            win_ini(n)%xsdb_path(k) = trim(xsdb_path(k))
            win_ini(n)%xs_scaling(k) = xs_scaling(k)
         enddo
         win_ini(n)%albflag = albflag 
         win_ini(n)%albedo = albedo
         win_ini(n)%IoffFlag = IoffFlag 
         win_ini(n)%spsh0flag =spsh0flag 
         win_ini(n)%spsh1flag = spsh1flag
         win_ini(n)%spsh2flag = spsh2flag
         win_ini(n)%sunsh0flag = sunsh0flag 
         if (win_ini(n)%wave_start > 1.d7/725.d0 .and. win_ini(n)%wave_stop < 1.d7/775.d0) then ! values were originally given as
                                                                                                ! wavenumbers [cm-1], that's why
                                                                                                ! the 1.d7 is here.
            win_ini(n)%Fsflag = flag%Fs   ! O2A-band
         else
            win_ini(n)%Fsflag = 0
         endif
      enddo

      !*** Aerosol retrieval settings and initialization values
      if(allocated(aerosol_ini)) deallocate(aerosol_ini)
      allocate(aerosol_ini(ntype_aer), stat=ierr)
      if (ierr.ne.0) then
         ierr = ierr_all
         write(message,*) 'READ_SETTINGS: memory allocation error'
         goto 999
      endif
      do i = 1 , ntype_aer
         allocate(aerosol_ini(i)%rm(2*nwin),aerosol_ini(i)%fim(2*nwin), stat=ierr)
         if (ierr.ne.0) then
            ierr = ierr_all
            write(message,*) 'READ_SETTINGS: memory allocation error'
            goto 999
         endif
         if (i == 1) then
            read(io, nml=aerosol1, iostat=ierr)     
            if (flag%output > 1) write(6, nml=aerosol1)   
         elseif (i == 2) then
            read(io, nml=aerosol2, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol2)
         elseif (i == 3) then
            read(io, nml=aerosol3, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol3)        
         elseif (i == 4) then
            read(io, nml=aerosol4, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol4)
         elseif (i == 5) then
            read(io, nml=aerosol5, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol5)
         elseif (i == 6) then
            read(io, nml=aerosol6, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol6)
         elseif (i == 7) then
            read(io, nml=aerosol7, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol7)
         elseif (i == 8) then
            read(io, nml=aerosol8, iostat=ierr) 
            if (flag%output > 1)  write(6, nml=aerosol8)   
         elseif (i == 9) then
            read(io, nml=aerosol9, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol9)      
         elseif (i == 10) then
            read(io, nml=aerosol10, iostat=ierr) 
            if (flag%output > 1) write(6, nml=aerosol10)
         endif

         if(ierr .ne. 0) then
            ierr = ierr_read
            write(message, *) 'READ_SETTINGS: error in namelist /aerosol/'
            goto 999
         endif
         !*** Put in type aerosol_ini       
         aerosol_ini(i)%CirrusFlag = CirrusFlag
         aerosol_ini(i)%AerosolFlags(5) = amount_fit
         aerosol_ini(i)%tau_ref = tau_ref
         aerosol_ini(i)%altid = height_dist
         aerosol_ini(i)%AerosolFlags(7) = alt1_fit
         aerosol_ini(i)%aeralt1 = alt1
         aerosol_ini(i)%AerosolFlags(8) = alt2_fit
         aerosol_ini(i)%aeralt2 = alt2
         aerosol_ini(i)%AerosolFlags(3) = rm_fit
         if (rm(2*nwin)==-999.) then !*** Refractive index constant within a spectral band            
            do n = 1, nwin  
               aerosol_ini(i)%rm(2*n-1:2*n) = rm(n)
            enddo
         else !*** 
            aerosol_ini(i)%rm(1:2*nwin) = rm(1:2*nwin)
         endif   
         aerosol_ini(i)%AerosolFlags(4) = fim_fit
         if (fim(2*nwin)==-999.) then !*** Refractive index constant within a spectral band
            do n = 1, nwin  
               aerosol_ini(i)%fim(2*n-1:2*n) = fim(n)
            enddo
         else           
            aerosol_ini(i)%fim(1:2*nwin) = fim(1:2*nwin)
         endif       
         aerosol_ini(i)%id = size_dist
         aerosol_ini(i)%AerosolFlags(6)  = shapefrac_fit 
         aerosol_ini(i)%shapefrac = shapefrac
         aerosol_ini(i)%AerosolFlags(1) = reff_fit 
         aerosol_ini(i)%reff = reff 
         aerosol_ini(i)%AerosolFlags(2) = veff_fit 
         aerosol_ini(i)%veff = veff 
         aerosol_ini(i)%tilt_angle = tilt_angle
      enddo
      close(io, iostat=ierr)
      if (ierr .ne. 0) then
         ierr = ierr_open
         write(message, *) 'READ_SETTINGS: error closing namelist input file'
         goto 999
      endif

      !*** Get dependent retrieval settings
      flag%O2 = 0
      do n = 1, nwin
         do k = 1, win_ini(n)%ntype
            if(win_ini(n)%type_x(k)==7) flag%O2=1
         enddo
      enddo

      !*** Calculate hi-resolution wavelength grid     
      do n = 1, nwin
         win_ini(n)%nwave_hi = int((win_ini(n)%wave_stop + win_ini(n)%fwhm*win_ini(n)%wvbd-&
              win_ini(n)%wave_start + win_ini(n)%fwhm*win_ini(n)%wvbd) / win_ini(n)%reso) + 1 

         allocate(win_ini(n)%wavelength_hi(win_ini(n)%nwave_hi), stat=ierr)
         if (ierr.ne.0) then
            ierr = ierr_all
            write(message,*) 'READ_SETTINGS: memory allocation error'
            goto 999
         endif
         forall(k=1:win_ini(n)%nwave_hi) win_ini(n)%wavelength_hi(k) =&
              dble(int(1./win_ini(n)%reso*(win_ini(n)%wave_start - win_ini(n)%fwhm*win_ini(n)%wvbd)))*&
              win_ini(n)%reso + dble(k-1)*win_ini(n)%reso
      enddo

      call check_settings(win_ini, aerosol_ini, flag, ierr)
      if (ierr .ne. 0) then 
         write(message,*) 'READ_SETTINGS: invalid retrieval settings'
         goto 999
      endif

      if (flag%output > 1) then
         write(message,'(a)')  '*** End of READ_SETTINGS ***'  
         call writelog(message, 1)
      endif

      ierr = 0 ! Mark success
      return

999   continue
      call stopretrieval(message)
      return

    end subroutine read_settings

!------------------------------------------------------------------------------  
!> @brief Check validity of the retrieval settings
!------------------------------------------------------------------------------  
    subroutine check_settings(win_ini, aerosol_ini, flag, ierr)
      type(aero), dimension(:), intent(inout) :: aerosol_ini
      type(window_ini), dimension(:), intent(in) :: win_ini   
      type(settings_flags), intent(in) :: flag
      integer, intent(out) :: ierr
      !*** local
      integer :: i, k, n, ntype_aer, nwin
      character(199) :: message

      nwin = size(win_ini)
      ntype_aer = size(aerosol_ini)

      do i = 1, ntype_aer
         if(aerosol_ini(i)%altid==1)then	       
            if(aerosol_ini(i)%AerosolFlags(3)/=0)then
               ierr = ierr_var
               write(message,*) 'CHECK_SETTINGS: AEROSOL BOUNDARY LAYER HEIGHT CANNOT BE FITTED.'
               goto 999
            endif
         elseif(aerosol_ini(i)%altid==2)then
            if(aerosol_ini(i)%AerosolFlags(3)/=0 .or. aerosol_ini(i)%AerosolFlags(4)/=0)then
               ierr = ierr_var
               write(message,*) 'CHECK_SETTINGDS: AEROSOL BOX HEIGHT DISTRIBUTION PARAMETERS CANNOT BE FITTED.'
               goto 999
            endif
         elseif(aerosol_ini(i)%altid==6)then
            if(aerosol_ini(i)%AerosolFlags(3)/=0 .or. aerosol_ini(i)%AerosolFlags(4)/=0)then
               write(message,*) 'CHECK_SETTINGS: AEROSOL DOUBLE BOX HEIGHT DISTRIBUTION PARAMETERS CANNOT BE FITTED.'
               ierr = ierr_var
               goto 999
            endif
         elseif(aerosol_ini(i)%altid<0 .or. aerosol_ini(i)%altid>6)then	    
            write(message,*) 'CHECK_SETTINGS: AEROSOL HEIGHT DISTRIBUTION ID NOT IMPLEMENTED.'
            ierr = ierr_var
            goto 999  
         endif

         if(aerosol_ini(i)%CirrusFlag==0)then
	    if(aerosol_ini(i)%id==0)then
	       if(aerosol_ini(i)%AerosolFlags(2)/=0)then
                  write(message,*) 'CHECK_SETTINGS: AEROSOL POWER LOW SIZE DISTRIBUTION HAS 1 MODE ONLY.'
                  ierr = ierr_var
                  goto 999
	       endif
	       if(aerosol_ini(i)%veff<0.1)then
                  write(message,*) 'CHECK_SETTINGS: AEROSOL POWER LOW SIZE DISTRIBUTION REQUIRES UPPER BOUNDARY > 0.1 MICRON.'
                  ierr = ierr_var
                  goto 999
	       endif
            elseif(aerosol_ini(i)%id==2)then	   
               if(aerosol_ini(i)%AerosolFlags(2)/=0)then
                  ierr = ierr_var
                  write(message,*) 'CHECK_SETTINGS: AEROSOL GAMMA SIZE DISTRIBUTION HAS 1 MOD ONLY.'
                  goto 999
               endif
            elseif(aerosol_ini(i)%id<0 .or. aerosol_ini(i)%id>2)then	 
               ierr = ierr_var 
               write(message,*) 'CHECK_SETTINGS: AEROSOL SIZE DISTRIBUTION ID NOT IMPLEMENTED.'
               goto 999   
	    endif
         elseif(aerosol_ini(i)%CirrusFlag==1)then	   
	    if(aerosol_ini(i)%id==0)then
	       if(aerosol_ini(i)%AerosolFlags(2)/=0)then
                  ierr = ierr_var
                  write(message,*) 'CHECK_SETTINGS: CIRRUS POWER LOW SIZE DISTRIBUTION HAS 1 MOD ONLY.'
                  goto 999
	       endif
	    elseif(aerosol_ini(i)%id/=0)then
               ierr = ierr_var
               write(message,*) 'CHECK_SETTINGS: ONLY POWER LAW SIZE DISTRIBUTION IMPLEMENTED FOR CIRRUS.'
               goto 999
	    endif
            if (aerosol_ini(i)%AerosolFlags(6)/=0) then
               ierr = ierr_var
               write(message,*) 'CHECK_SETTINGS: SHAPEFRAC CANNOT BE FITTED FOR CIRRUS'
               goto 999
            endif
	 else	
            ierr = ierr_var    
            write(message,*) 'CHECK_SETTINGS: EITHER CIRRUS OR AEROSOL.'
            goto 999
	 endif
      enddo

      !****Check validity of the retrieval flags   
      if(flag%atm<0 .or. flag%atm > 7 ) then      ! JS: increased the number of valid atm-flags
         ierr = ierr_var 
         write(message,*) 'CHECK_SETTINGS: which atmosphere?'
         goto 999
      endif

      !*** Scalar RT
!!$      nstokes = 1
!!$      nper = 1
!!$      nstokes_der = 1      
      if(flag%scat<0 .or. flag%scat >6) then
         ierr = ierr_var 
         write(message,*) 'CHECK_SETTINGS: scatflag not valid'
         goto 999
      elseif (flag%scat == 0) then
         write(message,'(a)') 'CHECK_SETTINGS: aerosol set to zero for non-scattering RTM '
         call writelog(message, 5)
         do i = 1, ntype_aer
            aerosol_ini(i)%AerosolFlags(:) = 0
            aerosol_ini(i)%tau_ref = 0.d0
         enddo
!!$      elseif(flag%scat>=2) then !Vector RT
!!$         nstokes = 3
!!$         nper = 4
!!$         nstokes_der = 3         
      endif
      
      if(flag%XS/=-1 .and. flag%XS/=0 .and. flag%XS /=1) then
         ierr = ierr_var 
         write(message,*) 'CHECK_SETTINGS: READ XS FROM DATABASE OR REUSE FROM LAST RETRIEVAL?'
         goto 999
      endif
      if(flag%inv < -6 .or. flag%inv>3) then
         ierr = ierr_var 
         write(message,*) 'CHECK_SETTINGS: which inversion method?'
         goto 999       
      endif
      if(flag%temp/=0 .and. flag%temp /=1) then
         ierr = ierr_var 
         write(message,*) 'CHECK_SETTINGS: FIT T OR NOT?'
         goto 999
      endif

      do n = 1, ntype_aer
         do i = 1, size(aerosol_ini(i)%AerosolFlags)
            if(aerosol_ini(n)%AerosolFlags(i)/=0 .and. aerosol_ini(n)%AerosolFlags(i)/=1) then
               ierr = ierr_var 
               write(message,*) 'CHECK_SETTINGS: AEROSOL PARAMETER FLAGS NOT VALID.'
               goto 999
            endif
         enddo
      enddo

      do n=1,nwin
         if( (win_ini(n)%spsh0flag/=0 .and. abs(win_ini(n)%spsh0flag) /=1) .or.&
              (win_ini(n)%spsh1flag/=0 .and. abs(win_ini(n)%spsh1flag) /=1) .or.&
              (win_ini(n)%spsh2flag/=0 .and. win_ini(n)%spsh2flag /=1) .or.&
              (win_ini(n)%sunsh0flag/=0 .and. abs(win_ini(n)%sunsh0flag) /=1) ) then
            ierr = ierr_var 
            write(message,*) 'CHECK_SETTINGS: FIT PARAMETER FLAGS NOT VALID.'
            goto 999
         endif
      enddo

      do n=1,nwin
         do k=1,win_ini(n)%ntype
            if(win_ini(n)%type_x_flag(k)/=0 .and. win_ini(n)%type_x_flag(k) /=1) then
               ierr = ierr_var 
               write(message,*) 'CHECK_SETTINGS: TARGET ABSORBER FLAG NOT VALID.'
               goto 999
            endif
         enddo
      enddo

      do n=1, nwin
         if(win_ini(n)%albflag<0 .or. win_ini(n)%albflag>3) then
             ierr = ierr_var 
            write(message,*) 'CHECK_SETTINGS: albflag not valid'    
            goto 999
         endif
         if(win_ini(n)%albedo(1) < 0.d0) then
            ierr = ierr_var 
            write(message,*) 'CHECK_SETTINGS: surface albedo not valid' 
            goto 999
         endif
      enddo

      ierr = 0 ! Mark success
      return

999   continue
      call stopretrieval(message)
      return

    end subroutine check_settings

!------------------------------------------------------------------------------
!> @brief Read cross-section lookup tables and put them in win_ini  
!------------------------------------------------------------------------------
    subroutine read_win_xsdb(flag, win_ini, ierr)
      !*** input
      type(settings_flags), intent(in) :: flag
      !*** Input/output
      type(window_ini), dimension(:), intent(inout) :: win_ini
      integer, intent(out) :: ierr
      !*** local variables
      type(cross_section_db) :: xsdb
      integer :: i, j, k, n, nwin
      character(stringlen) :: message
      !--------------------------------------------------------------------

      if(flag%output >=2) then
         write(message,'(a)') '*** Start reading cross-section databases ***'
         call writelog(message, 1)
      endif

      nwin = size(win_ini)
      !*** Start loop over all windows     
      do n = 1, nwin
!!$         if (allocated(win_ini(n)%xsdb)) deallocate(win_ini(n)%xsdb, stat=ierr)
!!$         if (ierr .ne. 0) then
!!$            write(message,*) 'READ_WIN_XSDB: memory deallocation error'
!!$            ierr = ierr_deall
!!$            goto 999
!!$         endif
         allocate(win_ini(n)%xsdb(win_ini(n)%ntype), stat=ierr)
         if (ierr .ne. 0) then
            write(message,*) 'READ_WIN_XSDB: memory allocation error'
            ierr = ierr_all
            goto 999
         endif
         !*** Start loop over all species 
         if (flag%xs<1)then

            do j = 1, win_ini(n)%ntype
               if(flag%output >=2 ) then
       	          write(message,'(a, i4)') 'molecule ', win_ini(n)%type_x(j)
                  call writelog(message, 1)
        	  write(message,'(a)')  'database '//trim(win_ini(n)%xsdb_path(j))
        	  call writelog(message, 1)
               endif

               xsdb%species=win_ini(n)%type_x(j)
               if (win_ini(n)%type_xsdb(j)==2) then    
        	  call read_xsdb_netcdf(win_ini(n)%xsdb_path(j), xsdb, ierr)
               elseif (win_ini(n)%type_xsdb(j)==3) then    
        	  call read_xsdb(win_ini(n)%xsdb_path(j), xsdb, ierr)
               elseif(win_ini(n)%type_xsdb(j)==0) then
                  call read_xsdb_frankenberg(win_ini(n)%xsdb_path(j), xsdb)
               elseif(win_ini(n)%type_xsdb(j)==1) then
        	  call read_xsdb_butz(win_ini(n)%xsdb_path(j), xsdb)           
               elseif(win_ini(n)%type_xsdb(j)==4) then
        	  call read_xsdb_absco(win_ini(n)%xsdb_path(j), xsdb)           
               endif
               if (ierr .ne. 0) then
                  write(message,*) 'READ_WIN_XSDB: error in reading XSDB LUT'                 
                  goto 999
               endif

               if (flag%xs_preprocess==1) then
                  call interpolate_xsdb(win_ini(n)%wavelength_hi, xsdb, ierr)
               elseif (flag%xs_preprocess==2) then
                  call tri_conv_xsdb(win_ini(n)%wavelength_hi, xsdb, ierr)
               elseif (flag%xs_preprocess==3) then
                  call rect_conv_xsdb(win_ini(n)%wavelength_hi, xsdb, ierr)
               endif
               if (ierr .ne. 0) then
                  write(message,*) 'READ_WIN_XSDB: error in interpolating XSDB LUT'                 
                  goto 999
               endif

               !*** Write temporary cross section database into window structure              
               win_ini(n)%xsdb(j)%species=xsdb%species
               win_ini(n)%xsdb(j)%nwv=xsdb%nwv

!!$               if(allocated(win_ini(n)%xsdb(j)%wv)) deallocate(win_ini(n)%xsdb(j)%wv)
!!$               if(allocated(win_ini(n)%xsdb(j)%p)) deallocate(win_ini(n)%xsdb(j)%p)
               allocate(win_ini(n)%xsdb(j)%wv(xsdb%nwv), &
                    win_ini(n)%xsdb(j)%p(size(xsdb%p)), win_ini(n)%xsdb(j)%vmr(xsdb%nvmr), stat=ierr)
               if (ierr .ne. 0) then
                  write(message,*) 'READ_WIN_XSDB: memory allocation error'
                  ierr = ierr_all
                  goto 999
               endif               
               win_ini(n)%xsdb(j)%wv=xsdb%wv
               win_ini(n)%xsdb(j)%nvmr = xsdb%nvmr            
               win_ini(n)%xsdb(j)%vmr = xsdb%vmr
               do i = 1, size(xsdb%p)

                  win_ini(n)%xsdb(j)%p(i)%p=xsdb%p(i)%p

!!$                  if(allocated(win_ini(n)%xsdb(j)%p(i)%T))deallocate(win_ini(n)%xsdb(j)%p(i)%T)
!!$                  if(allocated(win_ini(n)%xsdb(j)%p(i)%xs))deallocate(win_ini(n)%xsdb(j)%p(i)%xs)
                  allocate(win_ini(n)%xsdb(j)%p(i)%T(size(xsdb%p(i)%T)), &
                       win_ini(n)%xsdb(j)%p(i)%xs(xsdb%nwv, xsdb%nvmr, size(xsdb%p(i)%T)), stat=ierr)
                  if (ierr .ne. 0) then
                     write(message,*) 'READ_WIN_XSDB: memory allocation error'
                     ierr = ierr_all
                     goto 999
                  endif
                  do k = 1, size(xsdb%p(i)%T)
                     win_ini(n)%xsdb(j)%p(i)%T(k)=xsdb%p(i)%T(k)
                     ! HH: apply scaling of cross-sections:
                     win_ini(n)%xsdb(j)%p(i)%xs(:, :, k)=xsdb%p(i)%xs(:, :, k)*win_ini(n)%xs_scaling(j) 
                  enddo

                  deallocate(xsdb%p(i)%xs, xsdb%p(i)%T, stat=ierr)
                  if (ierr .ne. 0) then
                     write(message,*) 'READ_WIN_XSDB: memory deallocation error'
                     ierr = ierr_deall
                     goto 999
                  endif
               enddo

               deallocate(xsdb%p, xsdb%wv, stat=ierr)
               if (ierr .ne. 0) then
                  write(message,*) 'READ_WIN_XSDB: memory deallocation error'
                  ierr = ierr_deall
                  goto 999
               endif
               !*** End loop over all species      
            enddo
         endif
         !*** End loop over all windows      
      enddo

      if(flag%output >= 2) then
         write(message,'(a)')'*** End of reading cross-section databases ***'
         call writelog(message, 1)
      endif


      ierr = 0 ! Mark success
      return

999   continue
      call stopretrieval(message)
      return

    end subroutine read_win_xsdb
!------------------------------------------------------------------------------
end module read_settings_module
