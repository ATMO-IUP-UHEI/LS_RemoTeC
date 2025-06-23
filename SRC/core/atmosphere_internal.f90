!------------------------------------------------------------------------------
!> Read input atmosphere, interpolate onto model grids
!------------------------------------------------------------------------------
module atmosphere_internal_module
   use header_module
   use read_settings_module, only: window_ini, altitude_grid
   use forward_model_module, only: window_spectrum
   implicit none
   private

   !*** public types/procedures
   public :: atmospheric_scenario, altitude_grid
   public :: atmosphere_interpolate

   !------------------------------------------------------------------------------
   !> @brief Atmospheric input information and a priori profiles
   !> @todo move viewing geometry info to type spectrum
   !------------------------------------------------------------------------------
   type :: atmospheric_scenario
      !      private
      real(double), dimension(:), allocatable :: z      !< height profile [m]
      real(double), dimension(:), allocatable :: p      !< pressure profile [hPa]
      real(double), dimension(:), allocatable :: t      !< temperature profile [K]
      real(double), dimension(:), allocatable :: h2o    !< VMR H2O [1]
      real(double), dimension(:), allocatable :: co2    !< VMR CO2 [1]
      real(double), dimension(:), allocatable :: ch4    !< VMR CH4 [1]
      real(double), dimension(:), allocatable :: co     !< VMR CO  [1]
      real(double), dimension(:), allocatable :: n2o    !< VMR N2O [1]
      real(double), dimension(:), allocatable :: hcl    !< VMR HCl [1]
      real(double), dimension(:), allocatable :: hf     !< VMR HF  [1]
      real(double) :: surface_pressure                  !< surface pressure [hPa]
      real(double) :: surface_elevation                 !< surface elevaltion [m]
      real(double) :: surface_wspeed                    !< surface wind speed [m/s]
   end type atmospheric_scenario
   !------------------------------------------------------------------------------

contains
   !------------------------------------------------------------------------------
   !> @details Interpolate the input atmosphere to higher resolution pressure grids:
  !! nlay = number of retrieval layers
  !! nrt = number of layers used in RT code
  !! natm = number of layers for which absorption cross sections are calculated
  !! Calculate initial partial columns of absorbers (win%dv_x) from atmospheric input
   !------------------------------------------------------------------------------
   subroutine atmosphere_interpolate( &
      atm_in, &
      grid, &
      surface_elevation, &
      lat, &
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
      !*** input
      type(altitude_grid), intent(in) :: grid
      real(double), intent(in) ::  surface_elevation, lat
      type(atmospheric_scenario), intent(in) :: atm_in
      type(window_ini), dimension(:), intent(in) :: win_ini
      integer, intent(in) :: outputflag
      !***in/output
      type(window_spectrum), dimension(:), intent(inout) :: win
      !*** output
      type(atmosphere), intent(out) :: atm_retr
      type(atmosphere), intent(out) :: atm_rt
      type(atmosphere), intent(out) :: atm_xs
      real(double), dimension(:), allocatable, intent(out) :: dvair   ! Partial air column, subject to change in O2 retrieval (Dim: natm)
      real(double), intent(out) :: z_tropopause, z_bl
      integer, intent(out) :: ierr
      !*** local variables
      integer, parameter :: gas_units = 1   !1=VMR, 2=number density
      integer :: i, k, n, imin_psurf, i1, i2, nwin, io
      integer :: ninput
      integer, dimension(1) :: imin, imax, zmin
      integer :: nlay, nrt, natm, nstart(1), ip
      real(double) :: psurf, ptop, dp, dz, zsc
      real(double), dimension(:), allocatable :: vair, vh2o, vco2, vch4, vco, vhdo, vn2o, vhcl, vhf
      real(double), dimension(:), allocatable :: dvh2o, dvco2, dvo2, dvch4, dvco, dvhdo, dvn2o, dvhcl, dvhf 
      real(double), dimension(:), allocatable :: dT_dz, dh2o_dz
      real(double) :: rE, gE0
      real(double), dimension(:), allocatable :: gE
      real(double), dimension(:), allocatable :: zlev_atm   ! Layer boundaries (Dim: natm+1)
      real(double), dimension(:), allocatable :: plev_atm   ! Pressure, layer boundaries (Dim: natm+1)
      real(double), dimension(:), allocatable :: tlev_atm   ! Temperature, layer boundaries (Dim: natm+1)
      real(double), dimension(:), allocatable :: zlev_rt    ! Layer boundaries (Dim: nrt+1)
      real(double), dimension(:), allocatable :: plev_rt    ! Pressure, layer boundaries (Dim: nrt+1)
      real(double), dimension(:), allocatable :: tlev_rt    ! Temperature, layer boundaries (Dim: nrt+1)
      real(double), dimension(:), allocatable :: dz_atm     ! Layer thickness (Dim: natm)
      real(double), dimension(:), allocatable :: zlay_atm   ! Layer center (Dim: natm)
      type(atmospheric_scenario) :: atm_input
      character(stringlen) :: message
      !------------------------------------------------------------------------------
      if (outputflag >= 2) then
         write (message, '(a)') '*** Start ATMOSPHERE_INTERPOLATE ***'
         call writelog(message, 1)
      end if

      !*** Initialize error identifier
      ierr = 0
      !*** Dimensions
      nlay = grid%nlay
      nrt  = nlay*grid%nrt
      natm = nrt*grid%natm
      atm_retr%n = nlay
      atm_rt%n   = nrt
      atm_xs%n   = natm

      !*** local allocatables
      allocate (plev_atm(natm + 1), zlev_atm(natm + 1), tlev_atm(natm + 1), &
                plev_rt(nrt + 1), zlev_rt(nrt + 1), tlev_rt(nrt + 1), &
                zlay_atm(natm), dz_atm(natm), dT_dz(natm), dh2o_dz(natm), &
                vair(natm + 1), vh2o(natm), vco2(natm), vch4(natm), vco(natm), vhdo(natm), &
                vn2o(natm), vhcl(natm), vhf(natm), &
                dvh2o(natm), dvco2(natm), dvo2(natm), dvch4(natm), dvco(natm), dvhdo(natm), &
                dvn2o(natm), dvhcl(natm), dvhf(natm), &
                gE(natm + 1), &
                stat=ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
         ierr = ierr_all
         goto 999
      end if

      !*** arguments
      if (.not. allocated(atm_rt%p)) then
         allocate (atm_rt%p(nrt), atm_rt%z(nrt), atm_rt%dz(nrt), atm_rt%t(nrt), &
                   atm_xs%p(natm), atm_xs%t(natm), &
                   atm_retr%p(nlay + 1), atm_retr%z(nlay + 1), atm_retr%t(nlay + 1), &
                   stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
            ierr = ierr_all
            goto 999
         end if
      end if
      if (allocated(dvair)) deallocate (dvair)
      allocate (dvair(natm), stat=ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
         ierr = ierr_all
         goto 999
      end if

      !*** Correct for DEM surface elevation
      n = size(atm_in%z)
      dz = atm_in%z(n) - surface_elevation

      if (dz .ge. 0.) then  !ECMWF height > DEM height
       allocate(atm_input%p(n),  atm_input%t(n),  atm_input%z(n), atm_input%h2o(n),&
                atm_input%ch4(n),atm_input%co(n), atm_input%co2(n),&
                atm_input%n2o(n),atm_input%hcl(n),atm_input%hf(n), stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
            ierr = ierr_all
            goto 999
         end if
         zsc = atm_in%t(n)*rg/(air_m*grav)
         atm_input%p(n) = atm_in%p(n)*exp(dz/zsc)
         atm_input%p(1:n - 1) = atm_in%p(1:n - 1)
         atm_input%z(n) = surface_elevation
         atm_input%z(1:n - 1) = atm_in%z(1:n - 1)
         atm_input%t(n) = atm_in%t(n) + 0.0065*dz
         atm_input%t(1:n - 1) = atm_in%t(1:n - 1)
         atm_input%h2o  = atm_in%h2o
         atm_input%ch4  = atm_in%ch4
         atm_input%co   = atm_in%co
         atm_input%co2  = atm_in%co2
         atm_input%n2o  = atm_in%n2o
         atm_input%hcl  = atm_in%hcl
         atm_input%hf   = atm_in%hf   
      else               !ECMWF height < DEM height
         nstart = minloc(abs(atm_in%z - surface_elevation))
         ip = nstart(1)

         allocate (atm_input%p(ip),  atm_input%t(ip),  atm_input%z(ip), atm_input%h2o(ip), &
                   atm_input%ch4(ip),atm_input%co(ip), atm_input%co2(ip),&
                   atm_input%n2o(ip),atm_input%hcl(ip),atm_input%hf(ip), stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
            ierr = ierr_all
            goto 999
         end if
         dz = atm_in%z(ip) - surface_elevation
         zsc = atm_in%t(ip)*rg/(air_m*grav)
         atm_input%p(ip) = atm_in%p(ip)*exp(dz/zsc)
         atm_input%p(1:ip - 1) = atm_in%p(1:ip - 1)
         atm_input%z(ip) = surface_elevation
         atm_input%z(1:ip - 1) = atm_in%z(1:ip - 1)
         atm_input%t(ip) = atm_in%t(ip) + 0.0065*dz
         atm_input%t(1:ip - 1) = atm_in%t(1:ip - 1)
         atm_input%h2o(1:ip) = atm_in%h2o(1:ip)
         atm_input%ch4(1:ip) = atm_in%ch4(1:ip)
         atm_input%co(1:ip) = atm_in%co(1:ip)
         atm_input%co2(1:ip) = atm_in%co2(1:ip)
         atm_input%n2o(1:ip) = atm_in%n2o(1:ip)
         atm_input%hcl(1:ip) = atm_in%hcl(1:ip)
         atm_input%hf(1:ip) = atm_in%hf(1:ip)
      end if

      ninput = size(atm_input%p) - 1   !number of layers
      psurf = atm_input%p(ninput + 1)

      !*** Get pressure at layer boundaries for retrieval grid with nlay layers
      ptop = max(0.1d0, atm_input%p(1))     !*** constrain pressure grid to height = 0.1 mbar
      dp = (psurf - ptop)/nlay
      do k = 1, nlay
         atm_retr%p(k) = ptop + (k - 1)*dp
      end do
      atm_retr%p(nlay + 1) = psurf

      !*** Get pressure at layer boundaries for radiative transfer grid with nrt layers
      do i = 1, nlay
         dp = (atm_retr%p(i + 1) - atm_retr%p(i))/(nrt/nlay)
         do k = 1, nrt/nlay
            plev_rt((i - 1)*(nrt/nlay) + k) = atm_retr%p(i) + (k - 1)*dp
         end do
      end do
      plev_rt(nrt + 1) = psurf
      !*** Get pressure at layer centers for grid with nrt layers
      do k = 1, nrt
         atm_rt%p(k) = plev_rt(k)/2.+plev_rt(k + 1)/2.
      end do

      !*** Get pressure at layer boundaries for  grid with natm layers
      do i = 1, nlay
         dp = (atm_retr%p(i + 1) - atm_retr%p(i))/(natm/nlay)
         do k = 1, natm/nlay
            plev_atm((i - 1)*(natm/nlay) + k) = atm_retr%p(i) + (k - 1)*dp
         end do
      end do
      plev_atm(natm + 1) = psurf
      !*** Get pressure at layer centers for grid with natm layers
      do k = 1, natm
         atm_xs%p(k) = plev_atm(k)/2.+plev_atm(k + 1)/2.
      end do

      !*** Interpolate atmospheric input onto cross-section/radiative/retrieval grid
      !*** Cross-section grid at layer boundaries
      call spline_interpol(DLOG(atm_input%p), atm_input%z, ninput + 1, &
                           DLOG(plev_atm), zlev_atm, natm + 1, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      atm_xs%z = zlev_atm
      call spline_interpol(DLOG(atm_input%p), atm_input%t, ninput + 1, &
                           DLOG(plev_atm), tlev_atm, natm + 1, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      !*** Cross-section grid at layer centers
      call spline_interpol(DLOG(atm_input%p), atm_input%z, ninput + 1, &
                           DLOG(atm_xs%p), zlay_atm, natm, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      call spline_interpol(DLOG(atm_input%p), atm_input%t, ninput + 1, &
                           DLOG(atm_xs%p), atm_xs%t, natm, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      call linterp(atm_input%p, atm_input%h2o, ninput + 1, &
                   atm_xs%p, vh2o, natm)
      call linterp(atm_input%p, atm_input%co2, ninput + 1, &
                      atm_xs%p, vco2, natm)
      call linterp(atm_input%p, atm_input%ch4, ninput + 1, &
                   atm_xs%p, vch4, natm)
      call linterp(atm_input%p, atm_input%co, ninput + 1, &
                   atm_xs%p, vco, natm)
      call linterp(atm_input%p, atm_input%n2o, ninput + 1, &
                   atm_xs%p, vn2o, natm)
      call linterp(atm_input%p, atm_input%hcl, ninput + 1, &
                   atm_xs%p, vhcl, natm)
      call linterp(atm_input%p, atm_input%hf, ninput + 1, &
                   atm_xs%p, vhf, natm)
      
      !*** Radiative transfer grid at layer boundaries
      call spline_interpol(DLOG(atm_input%p), atm_input%z, ninput + 1, &
                           DLOG(plev_rt), zlev_rt, nrt + 1, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      call spline_interpol(DLOG(atm_input%p), atm_input%t, ninput + 1, &
                           DLOG(plev_rt), tlev_rt, nrt + 1, ierr)
      if (ierr .ne. 0) return
      !*** Radiative transfer grid at layer centers
      call spline_interpol(DLOG(atm_input%p), atm_input%z, ninput + 1, &
                           DLOG(atm_rt%p), atm_rt%z, nrt, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      call spline_interpol(DLOG(atm_input%p), atm_input%t, ninput + 1, &
                           DLOG(atm_rt%p), atm_rt%t, nrt, ierr)
      !*** Retrieval grid ar layer boundaries
      call spline_interpol(DLOG(atm_input%p), atm_input%z, ninput + 1, &
                           DLOG(atm_retr%p), atm_retr%z, nlay + 1, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if
      call spline_interpol(DLOG(atm_input%p), atm_input%t, ninput + 1, &
                           DLOG(atm_retr%p), atm_retr%t, nlay + 1, ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERNAL.SPLINE_INTERPOL.SPLINT: bad input'
         ierr = ierr_intrpl
         goto 999
      end if

      !*** Get layer thickness
      do k = 1, natm
         dz_atm(k) = DABS(zlev_atm(k) - zlev_atm(k + 1))
      end do
      do k = 1, nrt
         atm_rt%dz(k) = DABS(zlev_rt(k) - zlev_rt(k + 1))
      end do

      !*** If the heavy water isotope is retrieved, we assume that HDO = H2O profile times the natural default
      !*** abundance of 2d-3, which is intrinsic to HITRAN spectroscopy.
      vhdo = vh2o

      !*** Effective gravity
      !*** International gravity formula (Moritz, 1980):
      gE0 = 9.780327*(1.D0 + 5.3024D-3*sin(lat/180.D0*Pi)*sin(lat/180.D0*Pi) - 5.8D-6*sin(2*lat/180.D0*Pi)*sin(2*lat/180.D0*Pi))
      ! This seems buggy, rE is larger at poles than equator...
      !    rE = 6378137.D0/(1.006803-0.006706*sin(lat/180.D0*Pi)*sin(lat/180.D0*Pi))
      ! This is formula from Wikipedia:
      rE = 6.378137d6*sqrt( (1.d0-1.334395d-2*sin(lat/180.D0*Pi)*sin(lat/180.D0*Pi) )/(1.d0 - 6.694384d-3*sin(lat/180.D0*Pi)*sin(lat/180.D0*Pi)))
      gE = 9.80
      !*** Iterative formula taking into account dependence of geopotential height on gravity
      do i = 1, 5
         do k = 1, natm + 1
            gE(k) = gE0*(rE**2/(rE + zlev_atm(k)*gE(k)/9.80)**2)
         end do
      end do

      !*** Calculate partial columns
      if (gas_units == 1) then   ! Convert VMR to partial column
         !*** See Wunch et al., AMT, 2010 for calculation of the dry air column
         !*** Note: This is only correct if vh2o (and all other species) is the DRY VMR
         !*** If, on the other hand, we have the moist VMR for water, then vh2o_dry = vh2o_moist/(1.0 - vh2o_moist)
         do k = 1, natm
            dvair(k) = (plev_atm(k + 1) - plev_atm(k))*avoga*1.D-2/ &
                       (air_m* &
                        (gE(k) + gE(k + 1))/2.* &
                        (1.+vh2o(k)/1.60855)) ! Mair/Mwater=1.60855
         end do
         dvair(1) = dvair(1) + plev_atm(1)*avoga*1.D-2/ &
                    (air_m*gE(1)*(1.+vh2o(1)/1.60855)) !Mair/Mwater=1.60855
         do k = 1, natm
            dvh2o(k) = dvair(k)*vh2o(k)
            dvco2(k) = dvair(k)*vco2(k)
            dvch4(k) = dvair(k)*vch4(k)
            dvco(k)  = dvair(k)*vco(k)
            dvhdo(k) = dvair(k)*vhdo(k)
            dvn2o(k) = dvair(k)*vn2o(k)
            dvhcl(k) = dvair(k)*vhcl(k)
            dvhf(k)  = dvair(k)*vhf(k)
         end do
      elseif (gas_units == 2) then   ! Convert number density to partial column
         do k = 1, natm
            dvh2o(k) = dz_atm(k)*1.d2*vh2o(k)
            dvco2(k) = dz_atm(k)*1.d2*vco2(k)
            dvch4(k) = dz_atm(k)*1.d2*vch4(k)
            dvco(k) = dz_atm(k)*1.d2*vco(k)
            dvhdo(k) = dz_atm(k)*1.d2*vhdo(k)
            dvn2o(k) = dz_atm(k)*1.d2*vn2o(k)
            dvhcl(k) = dz_atm(k)*1.d2*vhcl(k)
            dvhf(k) = dz_atm(k)*1.d2*vhf(k)
         end do
         !*** Same as above, but use partial column of H2O instead of VMR
         do k = 1, natm
            dvair(k) = (plev_atm(k + 1) - plev_atm(k))*avoga*1.D-2/ &
                       (air_m*(gE(k) + gE(k + 1))/2.) - dvh2o(k)/1.60855
         end do
         dvair(1) = dvair(1) + plev_atm(1)*avoga*1.D-2/ &
                    (air_m*gE(1)) - dvh2o(1)/1.60855
      end if

      !*** Partical column of O2 is derived from air column
      dvo2 = dvair*relo2

      nwin = size(win_ini)
      !*** Define initial guess profiles for fitted absorbers
      do n = 1, nwin
         if (allocated(win(n)%dv_x)) deallocate (win(n)%dv_x)
         allocate (win(n)%dv_x(natm, win_ini(n)%ntype), stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'ATMOSPHERE_INTERPOLATE: memory allocation error'
            ierr = ierr_all
            goto 999
         end if
         do i = 1, win_ini(n)%ntype
            if (abs(win_ini(n)%type_x(i)) == 1 .or. &
                (abs(win_ini(n)%type_x(i)) >= 100 .and. abs(win_ini(n)%type_x(i)) < 199)) then
               do k = 1, natm
                  if (dvh2o(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE H2O VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = dvh2o(k)
                  end if
               end do
            else if (abs(win_ini(n)%type_x(i)) == 2 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 200 .and. abs(win_ini(n)%type_x(i)) < 299)) then
               do k = 1, natm
                  if (dvco2(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE CO2 VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = dvco2(k)
                  end if
               end do
            else if (abs(win_ini(n)%type_x(i)) == 5 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 500 .and. abs(win_ini(n)%type_x(i)) < 599)) then
               do k = 1, natm
                  if (dvco(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE CO VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                      win(n)%dv_x(k, i) = dvco2(k)/4000. ! Approximate CO profile from CO2
                     call writelog('WARNING: APRIORI CO VMR set to climatological value.', 5)
                 end if
               end do
            else if (abs(win_ini(n)%type_x(i)) == 6 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 600 .and. abs(win_ini(n)%type_x(i)) < 699)) then
               do k = 1, natm
                  if (dvch4(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE CH4 VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = dvch4(k)
                  end if
               end do
            else if (abs(win_ini(n)%type_x(i)) == 7 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 700 .and. abs(win_ini(n)%type_x(i)) < 799)) then
               win(n)%dv_x(:, i) = dvo2
            else if (abs(win_ini(n)%type_x(i)) == 181) then
               do k = 1, natm
                  if (dvhdo(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE HDO VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = dvhdo(k)
                  end if
               end do
            else if (abs(win_ini(n)%type_x(i)) == 4 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 400 .and. abs(win_ini(n)%type_x(i)) < 499)) then
               do k = 1, natm
                  if (dvn2o(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE N2O VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = dvch4(k)/6. ! Approximate N2O profile as 1/6 of CH4 profile
                  end if
               end do
           else if (abs(win_ini(n)%type_x(i)) == 14 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 1400 .and. abs(win_ini(n)%type_x(i)) < 1499)) then
               do k = 1, natm
                  if (dvhf(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE HF VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = 5.E16/natm ! Typical volcanic HF column, Butz et al., AMT, 2017, https://doi.org/10.5194/amt-10-1-2017
                  end if
               end do
          else if (abs(win_ini(n)%type_x(i)) == 15 .or. &
                     (abs(win_ini(n)%type_x(i)) >= 1500 .and. abs(win_ini(n)%type_x(i)) < 1599)) then
               do k = 1, natm
                  if (dvhcl(k) < 0.d0) then
                     call writelog('WARNING: NEGATIVE HCl VMR, set to zero', 5)
                     win(n)%dv_x(k, i) = 0.d0
                  else
                     win(n)%dv_x(k, i) = 2.E17/natm ! Typical volcanic HCl column, Butz et al., AMT, 2017, https://doi.org/10.5194/amt-10-1-2017
                  end if
               end do
            else
               if (outputflag >= 2) then
                  write (message, *) 'NO INITIAL GUESS FOR SPECIES #', win_ini(n)%type_x(i)
                  call writelog(message, 5)
               end if
               win(n)%dv_x(:, i) = 0.
            end if
         end do   ! Close loop over ntype absober types
      end do   ! Close loop over nwin windows

      !*** Height of tropopause is defined as layer with minimum temperature between 5 and 19 km
      imax = minloc(DABS(zlay_atm - 1.9D4))
      imin = minloc(DABS(zlay_atm - 5.D3))
      zmin = minloc(atm_xs%t(imax(1):imin(1)))
      z_tropopause = zlay_atm(imax(1) + zmin(1))
      !*** Define height of boundary layer as layer with maximum temperature gradient between
      !*** 500 m and 5 km above surface elevation
      imax = minloc(DABS(zlay_atm - (zlay_atm(natm) + 5.D3)))
      imin = minloc(DABS(zlay_atm - (zlay_atm(natm) + 5.D2)))
      do k = 1, natm
         dT_dz(k) = (tlev_atm(k + 1) - tlev_atm(k))/(zlev_atm(k + 1) - zlev_atm(k))
      end do
      !*** Or define height of boundary layer as layer with min H2O gradient? from GOSAT code
      dh2o_dz(1) = 0.D0
      DO k = 2, natm
         dh2o_dz(k) = (vh2o(k) - vh2o(k - 1))/(zlay_atm(k) - zlay_atm(k - 1))
      END DO
      !      zmin = maxloc(DABS(dT_dz(imax(1): imin(1))))
      zmin = minloc(dh2o_dz(imax(1):imin(1)))
      z_bl = zlay_atm(imax(1) + zmin(1) - 1)

      !*** Extra output
      if (outputflag >= 3) then
         open (newunit(io), FILE='./CONTRL_OUT/atm_lev_natm.dat')
         do k = 1, natm + 1
            write (io, FMT='(3(ES13.6,X))') zlev_atm(k), plev_atm(k), tlev_atm(k)
         end do
         close (io)
         open (newunit(io), FILE='./CONTRL_OUT/atm_lev_nrt.dat')
         do k = 1, nrt + 1
            write (io, FMT='(3(ES13.6,X))') zlev_rt(k), plev_rt(k), tlev_rt(k)
         end do
         close (io)
         open (newunit(io), FILE='./CONTRL_OUT/atm_lay_natm.dat')
         do k = 1, natm
            write (io, FMT='(10(ES13.6,X))') zlay_atm(k), atm_xs%p(k), atm_xs%t(k), & ! play_atm(k),tlay_atm(k),&
               dvair(k), dvh2o(k), dvco2(k), dvch4(k), &
               dvco(k), dh2o_dz(k), dT_dz(k)
         end do
         close (io)
         open (newunit(io), FILE='./CONTRL_OUT/atm_lay_nrt.dat')
         do k = 1, nrt
            write (io, FMT='(3(ES13.6,X))') atm_rt%z(k), atm_rt%p(k), atm_rt%t(k) !zlay_rt(k),play_rt(k),tlay_rt(k)
         end do
         close (io)
         open (newunit(io), FILE='./CONTRL_OUT/atm_lev_nlay.dat')
         do k = 1, nlay + 1
            write (io, FMT='(3(ES13.6,X))') atm_retr%z(k), atm_retr%p(k), atm_retr%t(k) !zlay_rt(k),play_rt(k),tlay_rt(k)
         end do
         close (io)
         open (newunit(io), FILE='./CONTRL_OUT/atm_input.dat')
         do k = 1, size(atm_input%z)
            write (io, FMT='(3(ES13.6,X))') atm_input%z(k), atm_input%p(k), atm_input%t(k)
         end do
         close (io)

      end if

      deallocate (vair, vh2o, vco2, vch4, vco, dvh2o, dvco2, dvo2, dvch4, dvco, vhdo, dvhdo, &
                  gE, zlev_atm, plev_atm, tlev_atm, zlev_rt, plev_rt, tlev_rt, dz_atm, &
                  zlay_atm, dT_dz, stat=ierr)
      if (ierr .ne. 0) then
         write (message, *) 'ATMOSPHERE_INTERPOLATE: memory deallocation error'
         ierr = ierr_deall
         goto 999
      end if

      if (outputflag >= 2) then
         write (message, '(a)') '*** End of ATMOSPHERE_INTERPOLATE ***'
         call writelog(message, 1)
      end if

      return

999   continue
      call stopretrieval(message)

   end subroutine atmosphere_interpolate
!------------------------------------------------------------------------------

end module atmosphere_internal_module
