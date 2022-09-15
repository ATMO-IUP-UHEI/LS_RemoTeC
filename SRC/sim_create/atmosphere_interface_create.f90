module atmo_interface_create_module
   use header_module
   use read_settings_module, only: window_ini
   use atmosphere_internal_module, only: atmospheric_scenario
   use spectrum_internal_module, only: spectrum
   use forward_model_module, only: aero, window_spectrum
   use auxiliary_routines_module, only: check
   use read_miprep_module
   use netcdf
   implicit none
   private

!*** procedures
   public :: atmosphere_input, output_atm, output_meteo, output_atm_nc_js

contains
!------------------------------------------------------------------------------
!> Read atmospheric inputfile and put in type
   subroutine atmosphere_input(infile, ipixel, atmflag, atm_input, meta, ierr, aerosolfile, aerosol_ini)
      character(*), intent(in) :: infile
      integer, intent(in) :: atmflag
      type(atmospheric_scenario), intent(out) :: atm_input
      type(metadata), intent(out) :: meta
      character(stringlen), optional, intent(out) :: aerosolfile
      integer, intent(out) :: ierr
      integer, optional, intent(in) :: ipixel
      type(aero), dimension(:), optional, intent(inout) :: aerosol_ini  !*** put cloud in datatype
      !*** local
      integer :: ninput
      real(double), dimension(:), allocatable :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable :: colev_in     ! VMR CO
      real(double) :: delta_angle, lat, lon
      !------------------------------------------------------------------------------

!!$       write(*,'(a)')
!!$       write(*,'(a)')'*****************************************'
!!$       write(*,'(a)')'***** Atmosphere interface'
!!$       write(*,'(a)')'*****'

      if (atmflag == 0) then
         call read_atm_custom( &
            infile, &
            heightlev_in, presslev_in, templev_in, &
            h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            lat, lon, &
            meta%time, atm_input%surface_elevation)
      elseif (atmflag == 1) then
         call read_atm_afgl( &
            infile, &
            heightlev_in, presslev_in, templev_in, &
            h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            lat, lon, &
            meta%time, atm_input%surface_elevation)
      elseif (atmflag == 2) then
         call read_atm_echam( &
            infile, &
            heightlev_in, presslev_in, templev_in, &
            h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            lat, lon, &
            meta%time, atm_input%surface_elevation)
      elseif (atmflag == 3) then
         call read_atm_sicor(infile, &
                             heightlev_in, presslev_in, templev_in, &
                             h2olev_in, co2lev_in, ch4lev_in, colev_in, &
                             lat, lon, &
                             meta%time, atm_input%surface_elevation)
      elseif (atmflag == 4) then
         if (present(aerosol_ini)) then
            call read_miprep(infile, ipixel, &
                             heightlev_in, presslev_in, templev_in, &
                             h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            meta, atm_input%surface_pressure, atm_input%surface_elevation, atm_input%surface_wspeed, ierr, aerosolfile, aerosol_ini)
         elseif (present(aerosolfile)) then
            call read_miprep(infile, ipixel, &
                             heightlev_in, presslev_in, templev_in, &
                             h2olev_in, co2lev_in, ch4lev_in, colev_in, &
                         meta, atm_input%surface_pressure, atm_input%surface_elevation, atm_input%surface_wspeed, ierr, aerosolfile)
         else
            call read_miprep(infile, ipixel, &
                             heightlev_in, presslev_in, templev_in, &
                             h2olev_in, co2lev_in, ch4lev_in, colev_in, &
                             meta, atm_input%surface_pressure, atm_input%surface_elevation, atm_input%surface_wspeed, ierr)
         end if
      elseif (atmflag == 5) then  ! Read atmospheric data for Indianapolis with Gaussian plume dispersion
         call read_atm_gp( &
            infile, &
            heightlev_in, presslev_in, templev_in, &
            h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            lat, lon, &
            meta%time, atm_input%surface_elevation)
      elseif (atmflag == 6 .or. atmflag == 7) then   ! Read atmospheric data for Indianapolis with LES plume dispersion
         call read_atm_icon( &
            infile, &
            heightlev_in, presslev_in, templev_in, &
            h2olev_in, co2lev_in, ch4lev_in, colev_in, &
            lat, lon, &
            meta%time, atm_input%surface_elevation)
      else
         call stopretrieval("ATMOSPHERE_INPUT: not a valid atmflag")
         return
      end if

      if (atmflag .NE. 4) then
         !*** Get latitude and longitude of pixel center and corners
         delta_angle = 0.d0                      ! footprint: TBD
         meta%lat(1) = lat                       ! center
         meta%lon(1) = lon                       ! center
         meta%lat(2) = meta%lat(1) + delta_angle ! right, upper corner
         meta%lon(2) = meta%lon(1) + delta_angle
         meta%lat(3) = meta%lat(1) - delta_angle ! right, lower corner
         meta%lon(3) = meta%lon(1) + delta_angle
         meta%lat(4) = meta%lat(1) - delta_angle ! left, lower corner
         meta%lon(4) = meta%lon(1) - delta_angle
         meta%lat(5) = meta%lat(1) + delta_angle ! left, upper corner
         meta%lon(5) = meta%lon(1) - delta_angle
         meta%surface_elevation = atm_input%surface_elevation
         meta%surface_elevation_stdv = 0.d0
         meta%landflag = 0
         meta%glintflag = 0
         meta%oceanglint = 0
         meta%Fs = 0.d0
         atm_input%surface_wspeed = 0.d0
         if (present(aerosolfile)) aerosolfile = infile
      end if

      !*** Deallocate
      if (allocated(atm_input%z)) deallocate (atm_input%z)
      if (allocated(atm_input%p)) deallocate (atm_input%p)
      if (allocated(atm_input%t)) deallocate (atm_input%t)
      if (allocated(atm_input%h2o)) deallocate (atm_input%h2o)
      if (allocated(atm_input%co2)) deallocate (atm_input%co2)
      if (allocated(atm_input%ch4)) deallocate (atm_input%ch4)
      if (allocated(atm_input%co)) deallocate (atm_input%co)
      ninput = size(heightlev_in)
      allocate ( &
         atm_input%z(ninput), &
         atm_input%p(ninput), &
         atm_input%t(ninput), &
         atm_input%h2o(ninput), &
         atm_input%co2(ninput), &
         atm_input%ch4(ninput), &
         atm_input%co(ninput))
      atm_input%z = heightlev_in
      atm_input%p = presslev_in
      atm_input%t = templev_in
      atm_input%h2o = h2olev_in
      atm_input%co2 = co2lev_in
      atm_input%ch4 = ch4lev_in
      atm_input%co = colev_in

!!$       write(*,'(a)')'***** Done.'
!!$       write(*,'(a)')'*****************************************'

   end subroutine atmosphere_input
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   subroutine read_atm_custom( &
      infile, &
      heightlev_in, presslev_in, templev_in, &
      h2olev_in, co2lev_in, ch4lev_in, colev_in, &
      lat, lon, &
      time, surface_elevation)
!*** input
      character(len=*), intent(in) :: infile
!*** output
      real(double), dimension(:), allocatable, intent(out) :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable, intent(out) :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable, intent(out) :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable, intent(out) :: colev_in     ! VMR CO
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
!*** local
      integer :: ninput, m, n, stat
      integer :: i, k, io, ierr
      real(double), dimension(:), allocatable :: heightlay_in ! Height at pressure-center of the layers [m]
      real(double), dimension(:), allocatable :: presslay_in  ! Pressure at the pressure-center of the layers [hPa]
      real(double), dimension(:), allocatable :: templay_in   ! Temperature at the pressure-center of the layers [K]
      real(double), dimension(:), allocatable :: h2olay_in    ! VMR H2O
      real(double), dimension(:), allocatable :: co2lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: ch4lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: colay_in     ! VMR CO
      real(double), dimension(:), allocatable :: air_in       ! Dry air column above levels [molec cm^-3]
      character(1) :: line

      open (newunit(io), FILE=trim(infile))
      n = 0
      m = 0
      read (io, *, iostat=stat) line
      if (line == '#') m = m + 1
      do while (stat .ne. -1)
         n = n + 1
         read (io, *, iostat=stat) line
         if (line == '#') m = m + 1
      end do
      close (io)
      ninput = n - m - 1

      allocate (heightlev_in(ninput + 1), presslev_in(ninput + 1), templev_in(ninput + 1), &
                heightlay_in(ninput), presslay_in(ninput), templay_in(ninput), &
                h2olev_in(ninput + 1), co2lev_in(ninput + 1), ch4lev_in(ninput + 1), colev_in(ninput + 1), &
                h2olay_in(ninput), co2lay_in(ninput), ch4lay_in(ninput), colay_in(ninput), air_in(ninput + 1))

      open (newunit(io), FILE=trim(infile))
      do i = 1, m
         read (io, *)
      end do
      do k = ninput + 1, 1, -1
         read (io, *, iostat=ierr) heightlev_in(k), presslev_in(k), templev_in(k), &
            h2olev_in(k), co2lev_in(k), ch4lev_in(k), colev_in(k)
      end do
      close (io)
      heightlev_in = heightlev_in*1000

      do k = 1, ninput
         presslay_in(k) = (presslev_in(k + 1) + presslev_in(k))/2.
      end do

      call spline_interpol(DLOG(presslev_in), heightlev_in, ninput + 1, &
                           DLOG(presslay_in), heightlay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), templev_in, ninput + 1, &
                           DLOG(presslay_in), templay_in, ninput, ierr)

      call linterp(presslev_in, h2olev_in, ninput + 1, &
                   presslay_in, h2olay_in, ninput)

      call linterp(presslev_in, co2lev_in, ninput + 1, &
                   presslay_in, co2lay_in, ninput)

      call linterp(presslev_in, ch4lev_in, ninput + 1, &
                   presslay_in, ch4lay_in, ninput)

      call linterp(presslev_in, colev_in, ninput + 1, &
                   presslay_in, colay_in, ninput)

      lat = 45.d0
      lon = 0.d0
      time = 0
      time(2) = 7  ! July
      surface_elevation = heightlev_in(ninput + 1)

      deallocate (air_in)

   end subroutine read_atm_custom
!------------------------------------------------------------------------------
   subroutine read_atm_afgl( &
      infile, &
      heightlev_in, presslev_in, templev_in, &
      h2olev_in, co2lev_in, ch4lev_in, colev_in, &
      lat, lon, &
      time, surface_elevation)
!*** input
      character(len=*), intent(in) :: infile
!*** output
      real(double), dimension(:), allocatable, intent(out) :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable, intent(out) :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable, intent(out) :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable, intent(out) :: colev_in     ! VMR CO
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
!*** local
      integer, parameter :: ninput = 49
      integer :: k, io, ierr
      real(double), dimension(:), allocatable :: heightlay_in ! Height at pressure-center of the layers [m]
      real(double), dimension(:), allocatable :: presslay_in  ! Pressure at the pressure-center of the layers [hPa]
      real(double), dimension(:), allocatable :: templay_in   ! Temperature at the pressure-center of the layers [K]
      real(double), dimension(:), allocatable :: h2olay_in    ! VMR H2O
      real(double), dimension(:), allocatable :: co2lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: ch4lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: colay_in     ! VMR CO
      real(double), dimension(:), allocatable :: air_in       ! Dry air column above levels [molec cm^-3]
      real(double) :: dummy

      allocate (heightlev_in(ninput + 1), presslev_in(ninput + 1), templev_in(ninput + 1), &
                heightlay_in(ninput), presslay_in(ninput), templay_in(ninput), &
                h2olev_in(ninput + 1), co2lev_in(ninput + 1), ch4lev_in(ninput + 1), colev_in(ninput + 1), &
                h2olay_in(ninput), co2lay_in(ninput), ch4lay_in(ninput), colay_in(ninput), air_in(ninput + 1))

      open (newunit(io), FILE=trim(infile))
      read (io, *)
      read (io, *)
      do k = 1, ninput + 1
         read (io, *) heightlev_in(k), presslev_in(k), templev_in(k), air_in(k), dummy, dummy, h2olev_in(k), co2lev_in(k)
      end do
      close (io)
      heightlev_in = heightlev_in*1000
      do k = 1, ninput + 1
         h2olev_in(k) = h2olev_in(k)/air_in(k)
         co2lev_in(k) = co2lev_in(k)/air_in(k)
         ch4lev_in(k) = 1.75D-6
         colev_in(k) = 5.D-8
      end do

      do k = 1, ninput
         presslay_in(k) = (presslev_in(k + 1) + presslev_in(k))/2.
      end do

      call spline_interpol(DLOG(presslev_in), heightlev_in, ninput + 1, &
                           DLOG(presslay_in), heightlay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), templev_in, ninput + 1, &
                           DLOG(presslay_in), templay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), h2olev_in, ninput + 1, &
                           DLOG(presslay_in), h2olay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), co2lev_in, ninput + 1, &
                           DLOG(presslay_in), co2lay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), ch4lev_in, ninput + 1, &
                           DLOG(presslay_in), ch4lay_in, ninput, ierr)

      call spline_interpol(DLOG(presslev_in), colev_in, ninput + 1, &
                           DLOG(presslay_in), colay_in, ninput, ierr)

      lat = 0.d0
      lon = 0.d0
      time = 0
      surface_elevation = heightlev_in(ninput + 1)
      deallocate (air_in)

   end subroutine read_atm_afgl

!------------------------------------------------------------------------------
   subroutine read_atm_echam( &
      infile, &
      heightlev_in, presslev_in, templev_in, &
      h2olev_in, co2lev_in, ch4lev_in, colev_in, &
      lat, lon, &
      time, surface_elevation)
!*** input
      character(*), intent(in) :: infile
!*** output
      real(double), dimension(:), allocatable, intent(out) :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable, intent(out) :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable, intent(out) :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable, intent(out) :: colev_in     ! VMR CO
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
!*** local
      real(double), dimension(:), allocatable :: heightlay_in ! Height at pressure-center of the layers [m]
      real(double), dimension(:), allocatable :: presslay_in  ! Pressure at the pressure-center of the layers [hPa]
      real(double), dimension(:), allocatable :: templay_in   ! Temperature at the pressure-enter of the layers [K]
      real(double), dimension(:), allocatable :: h2olay_in    ! VMR H2O
      real(double), dimension(:), allocatable :: co2lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: ch4lay_in    ! VMR CO2
      real(double), dimension(:), allocatable :: colay_in     ! VMR CO2
      real(double), dimension(:), allocatable :: vair_in      ! Dry air column above levels [molec cm^-2]
      real(double), dimension(:), allocatable :: air_in       ! Dry air column above levels [molec cm^-3]
      integer :: i, k, io, year, month, day, hour, ierr
      logical :: lexist
      character(199) :: filename
      real(double) :: dummy, CT_co2_mean, TM4_ch4_mean, tc
      real(double), dimension(:), allocatable :: CT_p, CT_co2
      real(double), dimension(:), allocatable :: TM4_p, TM4_ch4, TM4_co
      integer, parameter :: ninput = 19, n_CT = 34, n_TM4 = 34
!------------------------------------------------------------------------------
      if (allocated(heightlev_in)) deallocate (heightlev_in)
      if (allocated(heightlay_in)) deallocate (heightlay_in)
      if (allocated(presslev_in)) deallocate (presslev_in)
      if (allocated(presslay_in)) deallocate (presslay_in)
      if (allocated(templev_in)) deallocate (templev_in)
      if (allocated(templay_in)) deallocate (templay_in)
      if (allocated(h2olev_in)) deallocate (h2olev_in)
      if (allocated(h2olay_in)) deallocate (h2olay_in)
      if (allocated(co2lev_in)) deallocate (co2lev_in)
      if (allocated(co2lay_in)) deallocate (co2lay_in)
      if (allocated(ch4lev_in)) deallocate (ch4lev_in)
      if (allocated(ch4lay_in)) deallocate (ch4lay_in)
      if (allocated(colev_in)) deallocate (colev_in)
      if (allocated(colay_in)) deallocate (colay_in)
      if (allocated(vair_in)) deallocate (vair_in)
      if (allocated(air_in)) deallocate (air_in)

      allocate (heightlev_in(ninput + 1), presslev_in(ninput + 1), templev_in(ninput + 1), &
                heightlay_in(ninput), presslay_in(ninput), templay_in(ninput), &
                h2olev_in(ninput + 1), co2lev_in(ninput + 1), ch4lev_in(ninput + 1), colev_in(ninput + 1), &
                h2olay_in(ninput), co2lay_in(ninput), ch4lay_in(ninput), colay_in(ninput))

      allocate (CT_p(n_CT), CT_co2(n_CT), TM4_p(n_TM4), TM4_ch4(n_TM4), TM4_co(n_TM4))

!*** Read ECHAM5+MODIS+CT+TM4 file
      lexist = .false.
      filename = trim(infile)
      inquire (file=filename, exist=lexist)
      if (lexist .eqv. .false.) then
         print *, 'WARNING: ATMOSPHERIC INPUT FILE NOT FOUND', filename
         stop
      end if
      open (newunit(io), FILE=trim(filename), action='read', status='old')
      do i = 1, 3
         read (io, *)
      end do
      read (io, *) year, month, day, hour, lat, lon
      do k = 1, 7
         read (io, *)
      end do
      read (io, *) dummy, presslev_in(ninput + 1), heightlev_in(ninput + 1)
      read (io, *)
      do k = 1, ninput
         read (io, *) presslay_in(k), heightlay_in(k), templay_in(k), dummy, h2olay_in(k)
      end do
      read (io, *)
      do k = 1, n_CT
         read (io, *) CT_p(k), CT_co2(k)
      end do
      read (io, *)
      do k = 1, n_TM4
         read (io, *) TM4_p(k), TM4_ch4(k), TM4_co(k)
      end do
      close (io)

      if (time(1) .eq. 0) then
         time(1) = year
      end if

      if (time(2) .eq. 0) then
         time(2) = month
      end if

      if (time(3) .eq. 0) then
         time(3) = day
      end if

!       time(4) = hour
!       time(5) = 0
!       time(6) = 0

!*** Scale tropospheric (p>250hPa) variability of true CO2 (CT on 1x1 deg -> Gosat 5km radius footprint)
      i = 0
      CT_co2_mean = 0.D0
      do k = 1, n_CT
         if (CT_p(k) > 250.) then
            CT_co2_mean = CT_co2_mean + CT_co2(k)
            i = i + 1
         end if
      end do
      CT_co2_mean = CT_co2_mean/i
      do k = 1, n_CT
         if (CT_p(k) > 250.) then
            CT_co2(k) = CT_co2_mean + (CT_co2(k) - CT_co2_mean)*3.*(CT_p(k)/CT_p(n_CT))
         end if
      end do

!*** Scale tropospheric (p>250hPa) variability of true CH4 (TM4 on 3x2 deg -> Gosat 5km radius footprint)
      i = 0
      TM4_ch4_mean = 0.D0
      do k = 1, n_TM4
         if (TM4_p(k) > 250.) then
            TM4_ch4_mean = TM4_ch4_mean + TM4_ch4(k)
            i = i + 1
         end if
      end do
      TM4_ch4_mean = TM4_ch4_mean/i

      do k = 1, n_TM4
         if (TM4_p(k) > 250.) then
            TM4_ch4(k) = TM4_ch4_mean + (TM4_ch4(k) - TM4_ch4_mean)*5.*(TM4_p(k)/TM4_p(n_TM4))
         end if
      end do

!*** Calculate water vapor mixing ratio
      do k = 1, ninput
         tc = (templay_in(k) - 273.15d0)
!*** Relative humidity -> H2O partial pressure by A.L. Buck, in Journal of Applied Meteorology, 1981:
         h2olay_in(k) = h2olay_in(k)/presslay_in(k)* &
                        (1.00072d0 + presslay_in(k)*(3.2d-6 + 5.9d-10*tc*tc))* &
                        6.1121d0*DEXP(tc*(18.729d0 - 0.0044d0*tc)/(tc + 257.87))
      end do

!*** Calculate pressure at layer boundaries
      do k = ninput, 1, -1
         presslev_in(k) = 2*presslay_in(k) - presslev_in(k + 1)
      end do
!*** constrain pressure grid to height = 1mbar:
      presslev_in(1) = 1.0d0

!*** Interpolate to layer centers (in terms of pressure)
      call linterp(CT_p, CT_co2, n_CT, &
                   presslay_in, co2lay_in, ninput)
      call linterp(TM4_p, TM4_ch4, n_TM4, &
                   presslay_in, ch4lay_in, ninput)
      call linterp(TM4_p, TM4_co, n_TM4, &
                   presslay_in, colay_in, ninput)

!*** Interpolate  to layer boundaries
      call spline_interpol(DLOG(presslay_in), heightlay_in, ninput, &
                           DLOG(presslev_in), heightlev_in, ninput + 1, ierr)
      call spline_interpol(DLOG(presslay_in), templay_in, ninput, &
                           DLOG(presslev_in), templev_in, ninput + 1, ierr)
      call linterp(presslay_in, h2olay_in, ninput, &
                   presslev_in, h2olev_in, ninput + 1)
      call linterp(CT_p, CT_co2, n_CT, &
                   presslev_in, co2lev_in, ninput + 1)
      call linterp(TM4_p, TM4_ch4, n_TM4, &
                   presslev_in, ch4lev_in, ninput + 1)
      call linterp(TM4_p, TM4_co, n_TM4, &
                   presslev_in, colev_in, ninput + 1)

      surface_elevation = heightlev_in(ninput + 1)

      deallocate (CT_p, CT_co2, TM4_p, TM4_ch4, TM4_co)

   end subroutine read_atm_echam

   subroutine read_atm_gp( &
      infile, &
      height_lev, press_lev, temp_lev, &
      h2o_lev, co2_lev, ch4_lev, co_lev, &
      lat, lon, &
      time, surface_elevation)
      use netcdf

      !*** input
      character(*), intent(in) :: infile
      !*** output
      real(double), dimension(:), allocatable, intent(out) :: height_lev ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: press_lev  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: temp_lev   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2o_lev    ! VMR H2O at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: co2_lev    ! VMR CO2 at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: ch4_lev    ! VMR CH4 at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: co_lev     ! VMR CO at layer boundaries [mol/mol]
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
      !*** local
      integer :: i, k, io, ierr, sx, sy, year, month, day
      integer :: ncid, dimid, varid
      integer :: nlay_gp, nlay_bg, nlev_bg
      real(double) :: dx_gp, dy_gp, dz_gp
      real(double) :: surface_elevation_in(1), lat_in(1), lon_in(1)
      character(stringlen) :: fname, info

      real(double), dimension(:), allocatable :: co2_lay_gp          ! CO2 concentration in Gaussian plume model layer
      real(double), dimension(:), allocatable :: height_lay_gp       ! Height at layer center in Gaussian plume model
      real(double), dimension(:), allocatable :: co2_lay             ! CO2 mole fraction in background data model layer
      real(double), dimension(:), allocatable :: airmass_lay         ! Airmass in background model layer
      real(double), dimension(:), allocatable :: temp_lay            ! Temperature at layer center in background data model
      real(double), dimension(:), allocatable :: h2o_lay             ! Specific humidity in background model layer
      real(double), dimension(:), allocatable :: press_lay           ! Pressure at layer center in background model
      real(double), dimension(:), allocatable :: mass_dair, mol_dair, mass_dair2  ! mass and amount of moles of dry air

      real(double) :: molmass_dair = 28.9645  ! [g/mol] ! Molar mass of dry air
      real(double) :: molmass_h2o = 18.01528 ! [g/mol] ! Molar mass of H2O
      real(double) :: molmass_co2 = 44.0095  ! [g/mol] ! Molar mass of CO2

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      ! Extract name of file with high-resolution input data
      i = INDEX(infile, '.nc')
      fname = infile(1:i + 2)
      info = infile(i + 4:)

      ! Read year, month and time from filename
      i = INDEX(fname, 'RTC_INP_INDY_')
      read (fname(i + 16:i + 19), '(I4)') year
      read (fname(i + 20:i + 21), '(I2)') month
      read (fname(i + 22:i + 23), '(I2)') day

      ! Extract index in x-dimension to be read
      i = INDEX(info, 'X')
      read (info(i + 1:i + 3), '(I3)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(info, 'Y')
      read (info(i + 1:i + 3), '(I3)') sy

      if (time(1) .eq. 0) then
         time(1) = year
      end if

      if (time(2) .eq. 0) then
         time(2) = month
      end if

      if (time(3) .eq. 0) then
         time(3) = day
      end if

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      !***  -->  START: READ DATA FROM NETCDF-FILE  <--  ***!

      !*** Open NetCDF file
      call check(nf90_open(trim(fname), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) return

      !*** Read dimensions
      ! zlay_gp
      call check(nf90_inq_dimid(ncid, "zlay_gp", dimid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inquire_dimension(ncid, dimid, len=nlay_gp), ierr)
      if (ierr .ne. 0) return

      ! zlay_ct
      call check(nf90_inq_dimid(ncid, "zlay_ct", dimid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inquire_dimension(ncid, dimid, len=nlay_bg), ierr)
      if (ierr .ne. 0) return

      ! zlev_ct
      call check(nf90_inq_dimid(ncid, "zlev_ct", dimid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inquire_dimension(ncid, dimid, len=nlev_bg), ierr)
      if (ierr .ne. 0) return

      !*** Allocate sufficient memory to read in vertical profiles
      if (allocated(co2_lay_gp)) deallocate (co2_lay_gp)
      if (allocated(height_lay_gp)) deallocate (height_lay_gp)

      if (allocated(co2_lay)) deallocate (co2_lay)
      if (allocated(airmass_lay)) deallocate (airmass_lay)
      if (allocated(temp_lay)) deallocate (temp_lay)
      if (allocated(h2o_lay)) deallocate (h2o_lay)

      if (allocated(press_lev)) deallocate (press_lev)
      if (allocated(height_lev)) deallocate (height_lev)

      allocate (co2_lay_gp(nlay_gp), height_lay_gp(nlay_gp), &
                co2_lay(nlay_bg), airmass_lay(nlay_bg), temp_lay(nlay_bg), h2o_lay(nlay_bg), &
                press_lev(nlev_bg), height_lev(nlev_bg))

      !*** 3D variables [zlay_gp,lon,lat]
      ! Read CO2 concentration from Gaussian plume field [g/m3]
      call check(NF90_INQ_VARID(ncid, "co2_gp", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, co2_lay_gp, [1, sx, sy], [nlay_gp, 1, 1]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_ATT(ncid, varid, 'alongwind_spacing', dx_gp), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_ATT(ncid, varid, 'acrosswind_spacing', dy_gp), ierr)
      if (ierr .ne. 0) return

      !*** 2D variables [lat,lon]
      ! Read latitude [degrees_north]
      call check(NF90_INQ_VARID(ncid, "latitude", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, lat_in, [sx, sy], [1, 1]), ierr)
      if (ierr .ne. 0) return
      lat = lat_in(1)

      ! Read longitude [degrees_east]
      call check(NF90_INQ_VARID(ncid, "longitude", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, lon_in, [sx, sy], [1, 1]), ierr)
      if (ierr .ne. 0) return
      lon = lon_in(1)

      !*** 1D variables [zlay_gp]
      ! Read vertical dimension/grid of Gaussian plume model [m]
      call check(NF90_INQ_VARID(ncid, "z_gp", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, height_lay_gp), ierr)
      if (ierr .ne. 0) return

      !*** 1D variables [zlay_bg]
      ! Read background CO2 mole fraction [micromol/mol]
      call check(NF90_INQ_VARID(ncid, "co2_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, co2_lay), ierr)
      if (ierr .ne. 0) return
      co2_lay = co2_lay/1.0e6  ! Unit conversion: micromol/mol -> mol/mol

      ! Read background air mass [Kg/m2]
      call check(NF90_INQ_VARID(ncid, "airmass_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, airmass_lay), ierr)
      if (ierr .ne. 0) return
      airmass_lay = airmass_lay*1000.  ! Unit conversion: Kg/m2 -> g/m2

      ! Read background temperature [K]
      call check(NF90_INQ_VARID(ncid, "t_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, temp_lay), ierr)
      if (ierr .ne. 0) return

      ! Read background specific humidity [kg/kg]
      call check(NF90_INQ_VARID(ncid, "q_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, h2o_lay), ierr)
      if (ierr .ne. 0) return

      !*** 1D variables [zlev_bg]
      ! Read background temperature [K]
      call check(NF90_INQ_VARID(ncid, "p_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, press_lev), ierr)
      if (ierr .ne. 0) return

      ! Read background geopotential height [m]
      call check(NF90_INQ_VARID(ncid, "gph_bg", varid), ierr)
      if (ierr .ne. 0) return

      call check(NF90_GET_VAR(ncid, varid, height_lev), ierr)
      if (ierr .ne. 0) return

      call check(nf90_close(ncid), ierr)
      if (ierr .ne. 0) return

      !*** Flip data along the veritcal dimension such that the arrays start at the highest altitude
      if (height_lay_gp(1) .lt. height_lay_gp(2)) then
         co2_lay_gp = co2_lay_gp(nlay_gp:1:-1)
         height_lay_gp = height_lay_gp(nlay_gp:1:-1)
      end if

      if (height_lev(1) .lt. height_lev(2)) then
         co2_lay = co2_lay(nlay_bg:1:-1)
         airmass_lay = airmass_lay(nlay_bg:1:-1)
         temp_lay = temp_lay(nlay_bg:1:-1)
         h2o_lay = h2o_lay(nlay_bg:1:-1)

         press_lev = press_lev(nlev_bg:1:-1)
         height_lev = height_lev(nlev_bg:1:-1)
      end if

      surface_elevation = height_lev(nlev_bg)
      height_lay_gp = height_lay_gp + surface_elevation   ! Assume that the plume model's vertical grid starts at the land surface (i.e. at the surface elevation)

      !*** Calulate pressure at background model layer centers (w.r.t. pressure)
      if (allocated(press_lay)) deallocate (press_lay)
      allocate (press_lay(nlay_bg))

      do k = 1, nlay_bg
         press_lay(k) = (press_lev(k) + press_lev(k + 1))/2.
      end do

      !***  -->  END: READ DATA FROM NETCDF-FILE  <--  ***!

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      !*** --> START: CALCULATE THE AMOUNT OF THE ATMOSPHERIC CONSTITUENTS (IN MOLE) IN EACH MODEL GRIDBOX AND MERGE DATA <-- ***!

      !*** BACKGROUND DATA
      ! Allocate sufficient memory
      if (allocated(mass_dair)) deallocate (mass_dair)
      if (allocated(mol_dair)) deallocate (mol_dair)
      allocate (mass_dair(nlay_bg), mol_dair(nlay_bg))

      airmass_lay = airmass_lay*dx_gp*dy_gp ! [g]   Calculate the mass of background (moist) air in each gridbox given a high-resolution spatial grid and the coarser background model vertical grid
      h2o_lay = airmass_lay*h2o_lay     ! [g]   Calculate the corresponding mass of water vapour (airmass x specific humidity)
      mass_dair = airmass_lay - h2o_lay     ! [g]   Calculate the corresponding mass of dry air by subtracting the mass of water vapor from the total air mass
      mol_dair = mass_dair/molmass_dair  ! [mol] Calculate the corresponding number of moles dry air using the mass and molar mass of dry air (~29 g/mol)
      co2_lay = co2_lay*mol_dair        ! [mol] Calculate the corresponding number of moles CO2
      h2o_lay = h2o_lay/molmass_h2o     ! [mol] Calculate the corresponding number of moles water vapor using the mass molar mass of water (~18 g/mol)

      !*** CO2 PLUME DATA
      dz_gp = dabs(height_lay_gp(1) - height_lay_gp(2)) ! [m]   Calulate vertical spacing in Gaussian plume model
      co2_lay_gp = co2_lay_gp*dx_gp*dy_gp*dz_gp            ! [g]   Calculate the mass of plume CO2 in each Gaussian model gridbox
      co2_lay_gp = co2_lay_gp/molmass_co2                  ! [mol] Calculate the number of moles plume CO2 in each Gaussian model gridbox

      !*** MERGE BACKGROUND AND GAUSSIAN PLUME DATA
      ! Find the corresponding background layer by comparing the layer upper boundary height
      ! of the background model and the layer center height from the Gaussian plume model

      do k = 1, nlay_gp
         i = minloc(height_lev(1:nlev_bg - 1) - height_lay_gp(k), 1, ((height_lev(1:nlev_bg - 1) - height_lay_gp(k)) .ge. 0.0))  ! Skip last element in 'height_lev' in order to have the upper gridbox boundaries only
         co2_lay(i) = co2_lay(i) + co2_lay_gp(k)
         mol_dair(i) = mol_dair(i) + co2_lay_gp(k)
      end do

      !*** CALCULATE MOLE FRACTIONS / MIXING RATIOS W.R.T. DRY AIR (mol/mol)
      co2_lay = co2_lay/mol_dair
      h2o_lay = h2o_lay/mol_dair

      !*** --> END: CALCULATE THE AMOUNT OF THE ATMOSPHERIC CONSTITUENTS (IN MOLE) IN EACH MODEL GRIDBOX AND MERGE DATA <-- ***!

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      !*** --> START: INTERPOLATE DATA ONTO BACKGROUND DATA PRESSURE LEVEL BOUNDARY GRID <-- ***!

      !*** Allocate sufficent memory to variables
      if (allocated(temp_lev)) deallocate (temp_lev)
      if (allocated(h2o_lev)) deallocate (h2o_lev)
      if (allocated(co2_lev)) deallocate (co2_lev)
      if (allocated(ch4_lev)) deallocate (ch4_lev)
      if (allocated(co_lev)) deallocate (co_lev)

      allocate (temp_lev(nlev_bg), h2o_lev(nlev_bg), co2_lev(nlev_bg), ch4_lev(nlev_bg), co_lev(nlev_bg))

      !*** Interpolate
      call spline_interpol(DLOG(press_lay), temp_lay, nlay_bg, &
                           DLOG(press_lev), temp_lev, nlev_bg, ierr)

      call linterp(press_lay, h2o_lay, nlay_bg, &
                   press_lev, h2o_lev, nlev_bg)

      call linterp(press_lay, co2_lay, nlay_bg, &
                   press_lev, co2_lev, nlev_bg)

      deallocate (press_lay, temp_lay, &
                  h2o_lay, co2_lay)

      ch4_lev = 0.d0
      co_lev = 0.d0

      !*** --> END: INTERPOLATE DATA ONTO BACKGROUND DATA PRESSURE LEVEL BOUNDARY GRID <-- ***!

      !*** Convert pressure from Pa to hPa/mbar
      press_lev = press_lev/100.

   end subroutine read_atm_gp

   subroutine read_atm_icon( &
      infile, &
      height_lev, press_lev, temp_lev, &
      h2o_lev, co2_lev, ch4_lev, co_lev, &
      lat, lon, &
      time, surface_elevation)
      use netcdf

      !*** input
      character(*), intent(in) :: infile
      !*** output
      real(double), dimension(:), allocatable, intent(out) :: height_lev ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: press_lev  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: temp_lev   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2o_lev    ! VMR H2O at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: co2_lev    ! VMR CO2 at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: ch4_lev    ! VMR CH4 at layer boundaries [mol/mol]
      real(double), dimension(:), allocatable, intent(out) :: co_lev     ! VMR CO at layer boundaries [mol/mol]
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
      !*** local
      integer :: i, k, ierr, sx, sy
      integer :: ncid, dimid, varid
      integer :: nlev, ntime
      integer, dimension(:), allocatable :: datetime_lt
      character(stringlen) :: fname, info

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      ! Extract name of file with high-resolution input data
      i = INDEX(infile, '.nc')
      fname = infile(1:i + 2)

      info = infile(i + 4:)

      ! Extract index in x-dimension to be read
      i = INDEX(info, 'X')
      read (info(i + 1:i + 3), '(I3)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(info, 'Y')
      read (info(i + 1:i + 3), '(I3)') sy

      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !
      ! ----------------------------------------------------------------------------- !

      !***  -->  START: READ DATA FROM NETCDF-FILE  <--  ***!

      !*** Open NetCDF file
      call check(nf90_open(trim(fname), nf90_nowrite, ncid), ierr)
      if (ierr .ne. 0) return

      !*** Read dimensions
      !*** Vertical profiles
      call check(nf90_inq_dimid(ncid, "lev", dimid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inquire_dimension(ncid, dimid, len=nlev), ierr)
      if (ierr .ne. 0) return

      !*** Datetime array
      call check(nf90_inq_dimid(ncid, "time", dimid), ierr)
      if (ierr .ne. 0) return

      call check(nf90_inquire_dimension(ncid, dimid, len=ntime), ierr)
      if (ierr .ne. 0) return

      !*** Allocate sufficient memory to read in vertical profiles
      if (allocated(datetime_lt)) deallocate (datetime_lt)
      if (allocated(press_lev)) deallocate (press_lev)
      if (allocated(height_lev)) deallocate (height_lev)
      if (allocated(temp_lev)) deallocate (temp_lev)
      if (allocated(h2o_lev)) deallocate (h2o_lev)
      if (allocated(co2_lev)) deallocate (co2_lev)
      if (allocated(height_lev)) deallocate (height_lev)

      allocate (datetime_lt(ntime))
      allocate (press_lev(nlev), height_lev(nlev), temp_lev(nlev), h2o_lev(nlev), co2_lev(nlev))

      !*** UTC datetime array [year,month,day,hour,min,second]
      call check(NF90_INQ_VARID(ncid, "datetime_lt", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, datetime_lt), ierr)
      if (ierr .ne. 0) return

      !*** Latitude and longitude from 1D arrays [lat|lon]
      call check(NF90_INQ_VARID(ncid, "latitude", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, lat, [sy]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "longitude", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, lon, [sx]), ierr)
      if (ierr .ne. 0) return

      !*** Vertical profiles from 3D data [lev,lon,lat]
      call check(NF90_INQ_VARID(ncid, "pressure", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, press_lev, [sx, sy, 1], [1, 1, nlev]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "height", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, height_lev, [sx, sy, 1], [1, 1, nlev]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "temperature", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, temp_lev, [sx, sy, 1], [1, 1, nlev]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "h2o", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, h2o_lev, [sx, sy, 1], [1, 1, nlev]), ierr)
      if (ierr .ne. 0) return

      call check(NF90_INQ_VARID(ncid, "co2", varid), ierr)
      if (ierr .ne. 0) return
      call check(NF90_GET_VAR(ncid, varid, co2_lev, [sx, sy, 1], [1, 1, nlev]), ierr)
      if (ierr .ne. 0) return

      call check(nf90_close(ncid), ierr)
      if (ierr .ne. 0) return

      !***  -->  END: READ DATA FROM NETCDF-FILE  <--  ***!

      !*** Add datetime information to output variable ('time') if this has not been defined in input namelist-file (syn_create_XXXXXX.nml)
      do k = 1, ntime
         if (time(k) .eq. 0) then
            time(k) = datetime_lt(k)
         end if
      end do

      !*** surface elevation = height of lowermost model level
      surface_elevation = height_lev(nlev)

      !*** Add ch4 and co
      if (allocated(ch4_lev)) deallocate (ch4_lev)
      if (allocated(co_lev)) deallocate (co_lev)

      allocate (ch4_lev(nlev), co_lev(nlev))
      ch4_lev = 0.d0
      co_lev = 0.d0

   end subroutine read_atm_icon

!------------------------------------------------------------------------------
!> @details Read in synthetic atmospheric scenario on ECMWF pressure grid
!! Convert hybrid coefficients to pressure
!! Convert specific humidity to VMR of water
!! Get height profile from ideal gas law
!------------------------------------------------------------------------------
   subroutine read_atm_sicor( &
      infile, &
      heightlev_in, presslev_in, templev_in, &
      h2olev_in, co2lev_in, ch4lev_in, colev_in, &
      lat, lon, &
      time, surface_elevation)
!*** input
      character(len=*), intent(in) :: infile
!*** output
      real(double), dimension(:), allocatable, intent(out) :: heightlev_in ! Height at layer boundaries [m]
      real(double), dimension(:), allocatable, intent(out) :: presslev_in  ! Pressure at layer boundaries [hPa]
      real(double), dimension(:), allocatable, intent(out) :: templev_in   ! Temperature at layer boundaries [K]
      real(double), dimension(:), allocatable, intent(out) :: h2olev_in    ! VMR H2O
      real(double), dimension(:), allocatable, intent(out) :: co2lev_in    ! VMR CO2
      real(double), dimension(:), allocatable, intent(out) :: ch4lev_in    ! VMR CH4
      real(double), dimension(:), allocatable, intent(out) :: colev_in     ! VMR CO
      integer, intent(out) :: time(6)
      real(double), intent(out) :: surface_elevation, lat, lon
!*** Local variables
      integer :: i, nlev, io, ierr
      real(double), dimension(:), allocatable :: q, hyam, hybm
      real(double) :: psurf
      logical :: iopen

      !$OMP critical
      open (newunit(io), FILE=trim(infile), FORM='FORMATTED', status='old', iostat=ierr)
      !$OMP end critical
      if (ierr .ne. 0) goto 100

      read (io, *)
      read (io, *, iostat=ierr, err=100) surface_elevation
      read (io, *)
      read (io, *, iostat=ierr, err=100) psurf
      read (io, *)
      read (io, *, iostat=ierr, err=100) nlev
      allocate (heightlev_in(nlev), presslev_in(nlev), templev_in(nlev), &
                h2olev_in(nlev), co2lev_in(nlev), ch4lev_in(nlev), colev_in(nlev))
      allocate (hyam(nlev))
      allocate (hybm(nlev))
      allocate (q(nlev))

      !*** Read atmospheric data
      do i = 1, 5
         read (io, *, iostat=ierr, err=100)
      end do
      read (io, *, iostat=ierr, err=100) hyam
      read (io, *, iostat=ierr, err=100)
      read (io, *, iostat=ierr, err=100) hybm
      read (io, *, iostat=ierr, err=100)
      read (io, *, iostat=ierr, err=100) templev_in
      read (io, *, iostat=ierr, err=100)
      read (io, *, iostat=ierr, err=100) q
      read (io, *, iostat=ierr, err=100)
      read (io, *, iostat=ierr, err=100) colev_in
      read (io, *, iostat=ierr, err=100)
      read (io, *, iostat=ierr, err=100) ch4lev_in
      close (io)

      do i = 1, nlev
!*** Convert hybrid coefficients to pressure [from Pa to hPa]
         presslev_in(i) = (hyam(i) + hybm(i)*psurf)/100.
!*** Convert specific humidity to dry VMR of water
!*** Specific humidity s  is mass of water per mass of HUMID air.
!*** The mass mixing ratio s* is mass of water per mass of DRY air: s* = s / (1 - s).
!*** The volume mixing ratio is given by: (s*) * Mair / Mwater = s / (1 - s) * Mair / Mwater.
         h2olev_in(i) = q(i)/(1.0d0 - q(i))*1.60855
      end do

!*** Get height from barometric formula, i.e. ideal gas law
!*** from surface to top of atmopshere
      heightlev_in(nlev) = surface_elevation - &
                           rg*templev_in(nlev)/(air_m*grav)*log(presslev_in(nlev)/(psurf/100.))
      do i = 1, nlev - 1
         heightlev_in(nlev - i) = heightlev_in(nlev - i + 1) - &
                                  rg*templev_in(nlev - i)/(air_m*grav)*log(presslev_in(nlev - i)/presslev_in(nlev - i + 1))
      end do

      lat = 99.d0
      lon = 99.d0
      time = 0

100   if (ierr .ne. 0) then
         inquire (unit=io, opened=iopen)
         if (iopen) close (io)
         return
      end if

   end subroutine read_atm_sicor

!------------------------------------------------------------------------------
! Write old format meteo file
!------------------------------------------------------------------------------
   subroutine output_atm(atm_scenario, meta, sza, iza, phi, win_ini, win, aerosol, meteo_file)
      type(atmospheric_scenario), intent(in) :: atm_scenario
      type(metadata) :: meta
      character(len=*), intent(in) :: meteo_file
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(window_spectrum), dimension(:), intent(in):: win
      type(aero), dimension(:), intent(in) :: aerosol
      real(double), intent(in) :: sza, iza, phi
!*** local
      integer :: io, i, k, ninput, ntype_aer

!*** Write meteo data to output file
      open (newunit(io), file=meteo_file, form='formatted')
      ninput = size(atm_scenario%z) - 1
      write (io, *) ninput
      do k = 1, ninput + 1
         write (io, '(7ES23.15)') &
            atm_scenario%z(k), &
            atm_scenario%p(k), &
            atm_scenario%t(k), &
            atm_scenario%h2o(k), &
            atm_scenario%co2(k), &
            atm_scenario%ch4(k), &
            atm_scenario%co(k)
      end do
      write (io, *) ''
!*** Write instrument info
      write (io, '(7I6)') meta%time(1), meta%time(2), &
         meta%time(3), meta%time(4), meta%time(5), meta%time(6), meta%time(7)
      do k = 1, 5
         write (io, '(2ES23.15)') meta%lat(k), meta%lon(k)
      end do
      write (io, '(ES23.15)') meta%surface_elevation
      write (io, '(4ES23.15)') sza, iza, phi, 0.d0 !pixel_info%saz, pixel_info%iaz
      write (io, *) ''
!*** Write out some extra information for comparison
      ntype_aer = size(aerosol)

      write (io, '(I3)') size(win)

      if (ntype_aer == 0) then
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), &
            (0.d0, i=1, size(win)), & ! No OT - set strictly to zero to avoid bad representation in ascii-file
            (0.d0, i=1, size(win)), & ! No COT -        ----- || -----
            0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      else if (ntype_aer > 1) then
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), win(:)%ot, win(:)%cot, &
            aerosol(ntype_aer)%aeralt1, aerosol(ntype_aer)%aeralt2, &
            aerosol(ntype_aer)%reff, aerosol(ntype_aer)%tilt_angle, aerosol(ntype_aer)%shapefrac
      else
         !! J.S: Not sure if this is correct. If ntype_aer=1, either OT or COT will be zero!?
         !! This is incorrectly printed to the file as X.XXX-314 i.e. without the 'E'
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), win(:)%ot, win(:)%cot, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      end if
      close (io)

   end subroutine output_atm

!------------------------------------------------------------------------------
! Write old format meteo file in nc-format
!------------------------------------------------------------------------------

   subroutine output_atm_nc_js(atm_scenario, measurement, win, meta, meteo_file, index_info)
      type(atmospheric_scenario), intent(in) :: atm_scenario
      type(spectrum), dimension(:), intent(in) :: measurement
      type(window_spectrum), dimension(:), intent(in):: win
      type(metadata) :: meta
      character(len=*), intent(in) :: meteo_file, index_info
      !*** local
      integer :: i, n, sx, sy
      logical :: exst
      integer :: dimids_vert(2)
      integer :: z_id, p_id, t_id, h2o_id, co2_id, ch4_id, co_id, time_id, lat_id, lon_id, x_id, y_id, elev_id
      integer, dimension(:), allocatable :: grpid, wave_id, ot_id, cot_id, alb_id
      integer :: ncid, stat, dimid_nobs, dimid_time, dimid_z, dimid_wave, nwave, nwin, ngroup
      integer :: natm, lastindex
      character*1 :: ch
      character(stringlen) :: group_name

      ! Extract index in x-dimension to be read
      i = INDEX(index_info, 'X')
      read (index_info(i + 1:i + 3), '(I3)') sx

      ! Extract index in y-dimension to be read
      i = INDEX(index_info, 'Y')
      read (index_info(i + 1:i + 3), '(I3)') sy

      natm = size(atm_scenario%z)

      nwin = size(measurement)
      allocate (grpid(nwin), &
                wave_id(nwin), &
                ot_id(nwin), &
                cot_id(nwin), &
                alb_id(nwin), stat=stat)

      inquire (FILE=trim(meteo_file), EXIST=exst)

      if (.not. exst) then

         !*** Create the netCDF file.
         call check(nf90_create(trim(meteo_file), nf90_netcdf4, ncid), stat)

         !*** Unlimited dimension for number of observations
         call check(nf90_def_dim(ncid, "nobs", nf90_unlimited, dimid_nobs), stat)
         call check(nf90_def_dim(ncid, "nz", natm, dimid_z), stat)

         !*** Indexdata
         call check(nf90_def_var(ncid, "x", NF90_int, dimid_nobs, x_id), stat)
         call check(nf90_def_var(ncid, "y", NF90_int, dimid_nobs, y_id), stat)

         !*** Geodata
         call check(nf90_def_var(ncid, "latitude", NF90_double, dimid_nobs, lat_id), stat)
         call check(nf90_def_var(ncid, "longitude", NF90_double, dimid_nobs, lon_id), stat)
         call check(nf90_def_var(ncid, "elevation", NF90_double, dimid_nobs, elev_id), stat)

         call check(nf90_put_att(ncid, lat_id, "unit", "degrees_north"), stat)
         call check(nf90_put_att(ncid, lat_id, "description", "Latitude at pixel center"), stat)
         call check(nf90_put_att(ncid, lon_id, "unit", "degrees_east"), stat)
         call check(nf90_put_att(ncid, lon_id, "description", "Longitude at pixel center"), stat)
         call check(nf90_put_att(ncid, elev_id, "unit", "m"), stat)
         call check(nf90_put_att(ncid, elev_id, "description", "Surface altitude at pixel center"), stat)

         !*** Timedata
         call check(nf90_def_dim(ncid, "ntime", 6, dimid_time), stat)
         call check(nf90_def_var(ncid, "time", NF90_int, dimid_time, time_id), stat)
 call check(nf90_put_var(ncid, time_id, [meta%time(1), meta%time(2), meta%time(3), meta%time(4), meta%time(5), meta%time(6)]), stat)
         call check(nf90_put_att(ncid, time_id, "description", "Date and time as [YYYY,MM,DD,HOUR,MIN,SEC]"), stat)

         !*** Meteodata
         dimids_vert = (/dimid_z, dimid_nobs/)

         call check(nf90_def_var(ncid, "z", NF90_double, dimids_vert, z_id), stat)
         call check(nf90_def_var(ncid, "p", NF90_double, dimids_vert, p_id), stat)
         call check(nf90_def_var(ncid, "t", NF90_double, dimids_vert, t_id), stat)
         call check(nf90_def_var(ncid, "h2o", NF90_double, dimids_vert, h2o_id), stat)
         call check(nf90_def_var(ncid, "co2", NF90_double, dimids_vert, co2_id), stat)
         call check(nf90_def_var(ncid, "ch4", NF90_double, dimids_vert, ch4_id), stat)
         call check(nf90_def_var(ncid, "co", NF90_double, dimids_vert, co_id), stat)

         do n = 1, nwin
            nwave = measurement(n)%nwave
            !*** Create a group for each spectral window
            write (ch, '(i1.1)') n
            group_name = 'BAND'//ch
            call check(nf90_def_grp(ncid, trim(group_name), grpid(n)), stat)

            !*** Spectral data
            call check(nf90_def_dim(grpid(n), "nwave", nwave, dimid_wave), stat)

            !*** Define the variables
            call check(nf90_def_var(grpid(n), "wavelength", NF90_double, dimid_wave, wave_id(n)), stat)
            call check(nf90_def_var(grpid(n), "ot_inp", NF90_double, dimid_nobs, ot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "cot_inp", NF90_double, dimid_nobs, cot_id(n)), stat)
            call check(nf90_def_var(grpid(n), "alb_inp", NF90_double, dimid_nobs, alb_id(n)), stat)

            !*** Define attributes
            call check(nf90_put_att(grpid(n), wave_id(n), "unit", "nm"), stat)
            call check(nf90_put_att(grpid(n), wave_id(n), "description", "Wavelength grid of spectrum"), stat)
            call check(nf90_put_att(grpid(n), ot_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), ot_id(n), "description", "Total optical thickness used to simulate spectrum"), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), cot_id(n), "description", "Cirrus optical thickness used to simulate spectrum"), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "unit", "-"), stat)
            call check(nf90_put_att(grpid(n), alb_id(n), "description", "Surface albedo used to simulate spectrum"), stat)

            !*** write spectral grid
            call check(nf90_put_var(grpid(n), wave_id(n), measurement(n)%wavelength), stat)
         end do

         lastindex = 1

      else

         !*** Open the netCDF file and append
         call check(NF90_OPEN(trim(meteo_file), nf90_write, ncid), stat)
         call check(nf90_redef(ncid), stat)

         !*** Indexdata
         call check(nf90_inq_varid(ncid, "x", x_id), stat)
         call check(nf90_inq_varid(ncid, "y", y_id), stat)

         !*** Geodata
         call check(nf90_inq_varid(ncid, "latitude", lat_id), stat)
         call check(nf90_inq_varid(ncid, "longitude", lon_id), stat)
         call check(nf90_inq_varid(ncid, "elevation", elev_id), stat)

         !*** Meteodata
         call check(nf90_inq_varid(ncid, "z", z_id), stat)
         call check(nf90_inq_varid(ncid, "p", p_id), stat)
         call check(nf90_inq_varid(ncid, "t", t_id), stat)
         call check(nf90_inq_varid(ncid, "h2o", h2o_id), stat)
         call check(nf90_inq_varid(ncid, "co2", co2_id), stat)
         call check(nf90_inq_varid(ncid, "ch4", ch4_id), stat)
         call check(nf90_inq_varid(ncid, "co", co_id), stat)

         !*** Get group ID's
         call check(nf90_inq_grps(ncid, ngroup, grpid), stat)

         do n = 1, ngroup
            !*** Get variable ID's
            call check(NF90_INQ_DIMID(grpid(n), "nwave", dimid_wave), stat)
            call check(NF90_INQUIRE_DIMENSION(grpid(n), dimid_wave, len=nwave), stat)
            if (nwave .ne. measurement(n)%nwave) then
               call writelog("OUTPUT_ATM_NC_JS: nwave is not equal to number of spectral channels", 8)
            end if

            !*** Spectral data
            call check(NF90_INQ_VARID(grpid(n), "ot_inp", ot_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "cot_inp", cot_id(n)), stat)
            call check(NF90_INQ_VARID(grpid(n), "alb_inp", alb_id(n)), stat)
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

      !*** Geodata
      call check(nf90_put_var(ncid, lat_id, meta%lat(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, lon_id, meta%lon(1), start=(/lastindex/)), stat)
      call check(nf90_put_var(ncid, elev_id, meta%surface_elevation, start=(/lastindex/)), stat)

      !*** Meteodata
      call check(nf90_put_var(ncid, z_id, atm_scenario%z, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, p_id, atm_scenario%p, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, t_id, atm_scenario%t, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, h2o_id, atm_scenario%h2o, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, co2_id, atm_scenario%co2, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, ch4_id, atm_scenario%ch4, start=(/1, lastindex/)), stat)
      call check(nf90_put_var(ncid, co_id, atm_scenario%co, start=(/1, lastindex/)), stat)

      !*** Spectral data
      do n = 1, nwin
         call check(nf90_put_var(grpid(n), ot_id(n), win(n)%ot, start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), cot_id(n), win(n)%cot, start=(/lastindex/)), stat)
         call check(nf90_put_var(grpid(n), alb_id(n), win(n)%albedo(1), start=(/lastindex/)), stat)
      end do

      ! Close NetCDF file
      call check(nf90_close(ncid), stat)

      if (stat .ne. 0) then
         call writelog('OUTPUT_ATM_NC_JS: Error in writing meteo-data to netCDF file', 8)
      end if

   end subroutine output_atm_nc_js

!------------------------------------------------------------------------------
! Convert meteo data to ECMWF grid and write new format meteo file
!------------------------------------------------------------------------------
   subroutine output_meteo(atm_input, win_ini, win, aerosol, meteo_file)
!*** input
      type(atmospheric_scenario), intent(in) :: atm_input
      type(window_ini), dimension(:), intent(in) :: win_ini
      type(window_spectrum), dimension(:), intent(in):: win
      type(aero), dimension(:), intent(in) :: aerosol
      character(len=*), intent(in) :: meteo_file
!*** local
      integer :: i, imin_psurf, i1, i2, io, n, ntype_aer, ierr
      integer, parameter :: nlev = 91
      real(double), dimension(nlev) :: hyam, hybm
      real(double), dimension(nlev + 1) :: hyai, hybi
      real(double), dimension(:), allocatable:: q
      real(double), dimension(nlev) :: p_ecmwf, t_ecmwf, q_ecmwf, ch4_ecmwf, co_ecmwf
      real(double) :: psurf ! surface pressure
      real(double), parameter :: g = 9.80665  !gravitational
      data hyam/1.00002002716064, 2.99043607711792, 5.68400907516479, &
         10.1477527618408, 17.1609659194946, 27.683235168457, 42.8497295379639, &
         63.9571285247803, 92.4416084289551, 129.850791931152, 177.811737060547, &
         237.996978759766, 312.09049987793, 401.755142211914, 508.602508544922, &
         634.166290283203, 779.879577636719, 947.05615234375, 1136.8759765625, &
         1350.37469482422, 1588.43676757812, 1851.7919921875, 2141.01495361328, &
         2456.52697753906, 2798.60034179688, 3167.36401367188, 3562.81091308594, &
         3984.80627441406, 4433.09643554688, 4907.31811523438, 5407.00805664062, &
         5931.49780273438, 6479.783203125, 7050.59838867188, 7642.19799804688, &
         8253.77514648438, 8886.46484375, 9540.93310546875, 10216.2211914062, &
         10910.6831054688, 11622.5732421875, 12348.2797851562, 13083.5615234375, &
         13822.6176757812, 14557.34765625, 15280.3696289062, 15983.8071289062, &
         16660.091796875, 17301.9521484375, 17902.1552734375, 18453.9990234375, &
         18950.7568359375, 19386.029296875, 19753.6552734375, 20047.595703125, &
         20262.1552734375, 20391.537109375, 20430.1884765625, 20372.615234375, &
         20213.021484375, 19946.1943359375, 19567.06640625, 19073.798828125, &
         18470.0595703125, 17763.4462890625, 16965.0908203125, 16089.076171875, &
         15149.6059570312, 14159.4326171875, 13130.8017578125, 12075.775390625, &
         11007.2387695312, 9938.2666015625, 8880.73779296875, 7845.70874023438, &
         6844.54272460938, 5888.36181640625, 4986.50927734375, 4146.84008789062, &
         3376.8056640625, 2683.1748046875, 2070.52862548828, 1541.25537109375, &
         1096.42483520508, 735.753845214844, 456.543258666992, 249.407897949219, &
         108.125881195068, 30.3919818401337, 3.28939390194137, 0.00158000004012138/
      data hybm/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.36200000611097e-07, 7.09200001836052e-06, &
         3.42894009008887e-05, 9.30156493268441e-05, 0.000205124451895244, &
         0.000413634450524114, 0.000774259242461994, 0.00135060487082228, &
         0.00223289732821286, 0.00351588393095881, 0.00529460771940649, &
         0.00767857860773802, 0.0107716266065836, 0.0146839204244316, &
         0.0195241114124656, 0.0253994967788458, 0.032418629154563, &
         0.040686521679163, 0.050310181453824, 0.0613952092826366, &
         0.0740467458963394, 0.088370680809021, 0.104471534490585, &
         0.122456915676594, 0.142434172332287, 0.164512306451797, &
         0.188805609941483, 0.215417504310608, 0.244434677064419, &
         0.275773942470551, 0.309161394834518, 0.344265982508659, &
         0.380703672766685, 0.418055549263954, 0.455961376428604, &
         0.494148075580597, 0.53236910700798, 0.570387810468674, &
         0.607938021421432, 0.644746243953705, 0.680578589439392, &
         0.715223699808121, 0.74845165014267, 0.780032128095627, &
         0.809785097837448, 0.837567925453186, 0.863234400749207, &
         0.886642813682556, 0.907709091901779, 0.926403999328613, &
         0.942715436220169, 0.95664045214653, 0.968236565589905, &
         0.977852076292038, 0.985695540904999, 0.991678565740585, &
         0.995917141437531, 0.998815059661865/
      data hyai/0, 2.00004005432129, 3.98083209991455, 7.38718605041504, &
         12.9083194732666, 21.4136123657227, 33.9528579711914, 51.7466011047363, &
         76.1676559448242, 108.715560913086, 150.986022949219, 204.637451171875, &
         271.356506347656, 352.824493408203, 450.685791015625, 566.519226074219, &
         701.813354492188, 857.94580078125, 1036.16650390625, 1237.58544921875, &
         1463.16394042969, 1713.70959472656, 1989.87438964844, 2292.15551757812, &
         2620.8984375, 2976.30224609375, 3358.42578125, 3767.19604492188, &
         4202.41650390625, 4663.7763671875, 5150.85986328125, 5663.15625, &
         6199.83935546875, 6759.72705078125, 7341.4697265625, 7942.92626953125, &
         8564.6240234375, 9208.3056640625, 9873.560546875, 10558.8818359375, &
         11262.484375, 11982.662109375, 12713.8974609375, 13453.2255859375, &
         14192.009765625, 14922.685546875, 15638.0537109375, 16329.560546875, &
         16990.623046875, 17613.28125, 18191.029296875, 18716.96875, &
         19184.544921875, 19587.513671875, 19919.796875, 20175.39453125, &
         20348.916015625, 20434.158203125, 20426.21875, 20319.01171875, &
         20107.03125, 19785.357421875, 19348.775390625, 18798.822265625, &
         18141.296875, 17385.595703125, 16544.5859375, 15633.56640625, &
         14665.6455078125, 13653.2197265625, 12608.3837890625, 11543.1669921875, &
         10471.310546875, 9405.22265625, 8356.2529296875, 7335.16455078125, &
         6353.9208984375, 5422.802734375, 4550.2158203125, 3743.46435546875, &
         3010.14697265625, 2356.20263671875, 1784.85461425781, 1297.65612792969, &
         895.193542480469, 576.314147949219, 336.772369384766, 162.043426513672, &
         54.2083358764648, 6.57562780380249, 0.00316000008024275, 0/
      data hybi/0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.72400001222195e-07, &
         1.39116000354989e-05, 5.46672017662786e-05, 0.00013136409688741, &
         0.000278884806903079, 0.000548384094145149, 0.00100013439077884, &
         0.00170107535086572, 0.00276471930555999, 0.00426704855635762, &
         0.00632216688245535, 0.00903499033302069, 0.0125082628801465, &
         0.0168595779687166, 0.0221886448562145, 0.0286103487014771, &
         0.0362269096076488, 0.0451461337506771, 0.055474229156971, &
         0.0673161894083023, 0.0807773023843765, 0.0959640592336655, &
         0.112979009747505, 0.131934821605682, 0.152933523058891, &
         0.176091089844704, 0.201520130038261, 0.229314878582954, &
         0.259554475545883, 0.291993409395218, 0.326329380273819, 0.3622025847435, &
         0.399204760789871, 0.436906337738037, 0.475016415119171, &
         0.513279736042023, 0.551458477973938, 0.589317142963409, &
         0.626558899879456, 0.662933588027954, 0.69822359085083, &
         0.732223808765411, 0.764679491519928, 0.795384764671326, &
         0.824185431003571, 0.850950419902802, 0.875518381595612, &
         0.897767245769501, 0.917650938034058, 0.935157060623169, &
         0.950273811817169, 0.963007092475891, 0.973466038703918, &
         0.982238113880157, 0.98915296792984, 0.994204163551331, 0.99763011932373, 1/

!*** Find surface pressure for the input surface elevation
      imin_psurf = minval(minloc(DABS(atm_input%surface_elevation - atm_input%z)))
      if (atm_input%z(imin_psurf) .le. atm_input%surface_elevation) i1 = imin_psurf
      if (atm_input%z(imin_psurf) .le. atm_input%surface_elevation) i1 = imin_psurf + 1
      if (imin_psurf .eq. size(atm_input%z)) i1 = imin_psurf
      i2 = i1 - 1
      psurf = log(atm_input%p(imin_psurf)) + &
              (log(atm_input%p(i2)) - log(atm_input%p(i1)))/(atm_input%z(i2) - atm_input%z(i1)) &
              *(atm_input%surface_elevation - atm_input%z(imin_psurf))
      psurf = DEXP(psurf)

      n = size(atm_input%h2o(:))
      allocate (q(n))
!*** Get specific humidity from VMR of water
      do i = 1, n
         q(i) = atm_input%h2o(i)/(1.60855 - atm_input%h2o(i)) ! Mair/Mwater=1.60855
      end do

!*** pressure at mid-levels
      do i = 1, nlev
         p_ecmwf(i) = hyam(i)/100.+hybm(i)*psurf
      end do

!*** interpolate VMR of CO and CH4, specific humidity and temperature to ECMWF pressure-grid
      call spline_interpol(DLOG(atm_input%p), atm_input%t, n, &
                           DLOG(p_ecmwf), t_ecmwf, nlev, ierr)
      call linterp(atm_input%p, q, n, &
                   p_ecmwf, q_ecmwf, nlev)
      call linterp(atm_input%p, atm_input%ch4, n, &
                   p_ecmwf, ch4_ecmwf, nlev)
      call linterp(atm_input%p, atm_input%co, n, &
                   p_ecmwf, co_ecmwf, nlev)

      open (newunit(io), file=meteo_file, form='formatted')
      write (io, *) nlev                             ! number of layers
      write (io, *) psurf*100.d0                     ! [Pa]
      write (io, *) atm_input%surface_elevation*grav ! geopotenial at surface
      write (io, *) ''
      write (io, *) ''
      do i = 1, nlev
         write (io, '(8ES23.15)') &
            hyai(i), &
            hybi(i), &
            hyam(i), &
            hybm(i), &
            t_ecmwf(i), &
            q_ecmwf(i), &
            ch4_ecmwf(i), &
            co_ecmwf(i)
      end do
      write (io, '(2ES23.15)') hyai(nlev + 1), hybi(nlev + 1) !interface has 1 point more
      write (io, *) ''
      write (io, *) ''
      write (io, *) ''
      write (io, *) ''
      write (io, *) ''
      write (io, *) ''
!*** Write out some extra information for comparison
      ntype_aer = size(aerosol)

      write (io, '(I3)') size(win)

      if (ntype_aer == 0) then
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), &
            (0.d0, i=1, size(win)), & ! No OT - set strictly to zero to avoid bad representation in ascii-file
            (0.d0, i=1, size(win)), & ! No COT -        ----- || -----
            0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      else if (ntype_aer > 1) then
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), win(:)%ot, win(:)%cot, &
            aerosol(ntype_aer)%aeralt1, aerosol(ntype_aer)%aeralt2, &
            aerosol(ntype_aer)%reff, aerosol(ntype_aer)%tilt_angle, aerosol(ntype_aer)%shapefrac
      else
         !! J.S: Not sure if this is correct. If ntype_aer=1, either OT or COT will be zero!?
         !! This is incorrectly printed to the file as X.XXX-314 i.e. without the 'E'
         write (io, '(100ES15.7)') win_ini(:)%wave_start, win_ini(:)%wave_stop, &
            (win(i)%albedo(1), i=1, size(win)), win(:)%ot, win(:)%cot, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0
      end if
      close (io)

   end subroutine output_meteo

!------------------------------------------------------------------------------

end module atmo_interface_create_module
