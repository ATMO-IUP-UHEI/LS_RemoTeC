!------------------------------------------------------------------------------
!> Calculate initial aerosol information
!------------------------------------------------------------------------------
module aerosol_input_module
   use header_module
   use optic_cirrus_module, only: cirrus_table, optic_cirrus_xs
   use OpticM_module, only: mie_lut, modes_calc_xs
   implicit none
   private

   !*** Public types/procedures
   public :: aero
   public :: get_aerosol_properties_lognormal, set_altdis

   type :: aero
      real(double), dimension(:), allocatable :: alt_dis    ! Aerosol height distribution (Gaussian) (nrt)
      real(double), dimension(:), allocatable :: dalt_daer1 ! Derivative of aerosol height distribution wrt center(nrt)
      real(double), dimension(:), allocatable :: dalt_daer2 ! Derivative of aerosol height distribution wrt width(nrt)
      real(double) :: tilt_angle
      real(double) :: tau_ref                  ! Initial aot
      real(double) :: sig_ref                  ! Initial aerosol extinction cross section
      real(double) :: aer_col                  ! Aerosol column fine mode/power law/gamma
      real(double) :: shapefrac                ! Aerosol mode fraction
      real(double) :: aeralt1                       ! Aerosol height center
      real(double) :: aeralt2                       ! Aerosol height width
      real(double) :: reff                                    ! Effective radius fine mode/ power of power law distribution
      real(double) :: veff                                    ! Effective width fine mode
      real(double), dimension(:), allocatable :: rm    ! Real refractive index (nwin)
      real(double), dimension(:), allocatable :: fim   ! Imaginary refractive index (nwin)

      integer, dimension(:), allocatable :: nder       ! Indices of aerosol layers with significant contribution
      integer, dimension(8) :: AerosolFlags    ! Fit aerosol parameters no:0, yes:1:
      ! (1)r_eff/power mode1, (2)v_eff/max rad mode1,
      ! (3)real refr, (4)im refr,
      ! (5)column mode1, (6) shapefrac
      ! (7)height param1, (8)height param2
      integer :: naer                          ! Number of aerosol parameters to be retrieved
      integer :: CirrusFlag                    ! Cirrus aerosol, no: 0, yes: 1
      integer :: id                            ! Type of aerosol size distribution (0: power law, 1: lognormal, 2: gamma)
      integer :: altid                                    ! Aerosol height distribution type
      integer :: maxd                          ! Number of indices of aerosol layers with significant contribution
   end type aero

contains
   !------------------------------------------------------------------------------
   !> @brief Calculate the initial guess for aerosol parameters
   !> @details Some of the parameters are simply set to the values given in retrieval.ini
  !! (stored in aerosol_ini). The optical properties have to be calculated by calls
  !! to optical_cirrus_xs and modes_calc_xs. The height distribution is calculated in
  !! set_altdis.
   !> @param[in] aero_lut
   !> @param[in] cirrus_lut
   !> @param[in] nwin
   !> @param[in] z_tropopause
   !> @param[in] z_bl
   !> @param[in] atm_rt
   !> @param[in] aerosol_ini
   !> @param[out] aerosol
   !------------------------------------------------------------------------------
   subroutine get_aerosol_properties_lognormal( &
      aero_lut, cirrus_lut, &
      outputflag, glintflag, glintscat, nwin, &
      z_tropopause, z_bl, &
      atm_rt, &
      aerosol_ini, &
      aerosol, &
      ierr)
      !*** input
      type(Mie_lut), intent(in) :: aero_lut
      type(cirrus_table), intent(in) :: cirrus_lut
      integer, intent(in) :: outputflag, glintflag, glintscat, nwin
      real(double), intent(in) :: z_tropopause, z_bl
      type(atmosphere), intent(in) :: atm_rt
      type(aero), dimension(:), intent(in) :: aerosol_ini
      !*** output
      type(aero), dimension(:), allocatable, intent(out) :: aerosol
      integer, intent(out) :: ierr
      !*** local variables
      logical, parameter :: Trunc_Flag = .false.       ! HH: used to be true, but I think it was a bug!
      integer :: i, j, k, ntype_aer, nrt
      real(double) :: csca, cabs, rlambda, f
      real(double), dimension(npar_mie) :: aerosol_pars_mode
      character(stringlen) :: message
      !--------------------------------------------------------------------------------
      !*** Initialize
      ierr = 0
      ntype_aer = size(aerosol_ini)
      nrt = atm_rt%n

      if (allocated(aerosol)) deallocate (aerosol)
      allocate (aerosol(ntype_aer), stat=ierr)
      if (ierr .ne. 0) then
         write (message, *) 'GET_AEROSOL_PROPERTIES_LOGNORMAL: memory allocation problem'
         ierr = ierr_all
         goto 999
      end if

      do k = 1, ntype_aer
         allocate (aerosol(k)%rm(2*nwin), &        !<++++++++++++++++++++++++ Changed dim. to 2*nwin for spectral dependance
                   aerosol(k)%fim(2*nwin), &             !<++++++++++++++++++++++++ Changed dim. to 2*nwin for spectral dependance
                   stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'GET_AEROSOL_PROPERTIES_LOGNORMAL: memory allocation problem'
            ierr = ierr_all
            goto 999
         end if
         aerosol(k)%CirrusFlag = aerosol_ini(k)%CirrusFlag
         aerosol(k)%AerosolFlags = aerosol_ini(k)%AerosolFlags
         aerosol(k)%tau_ref = aerosol_ini(k)%tau_ref
         aerosol(k)%altid = aerosol_ini(k)%altid
         aerosol(k)%aeralt1 = aerosol_ini(k)%aeralt1
         if (aerosol_ini(k)%aeralt1 < 0 .and. aerosol_ini(k)%CirrusFlag .ne. 1) then
            aerosol(k)%aeralt1 = z_bl
         elseif (aerosol_ini(k)%aeralt1 < 0 .and. aerosol_ini(k)%CirrusFlag == 1) then
            aerosol(k)%aeralt1 = z_tropopause
         end if
         aerosol(k)%aeralt2 = aerosol_ini(k)%aeralt2
         aerosol(k)%rm = aerosol_ini(k)%rm
         aerosol(k)%fim = aerosol_ini(k)%fim
         aerosol(k)%id = aerosol_ini(k)%id
         aerosol(k)%shapefrac = aerosol_ini(k)%shapefrac
         aerosol(k)%reff = aerosol_ini(k)%reff
         aerosol(k)%veff = aerosol_ini(k)%veff
         aerosol(k)%tilt_angle = aerosol_ini(k)%tilt_angle

         if (glintflag == 1 .and. glintscat == 0) then
            aerosol(k)%tau_ref = 0.d0
            aerosol(k)%aeralt1 = 2.5D4
            if (outputflag >= 2) call writelog('GET_AEROSOL_PROPERTIES_LOGNORMAL: no aerosols for glint', 5)
         end if

         rlambda = 0.765
         !*** Calculate aerosol optical properties
         if (aerosol(k)%CirrusFlag == 1) then
            call optic_cirrus_xs( &
               cirrus_lut, &
               Trunc_Flag, &
               aerosol(k)%tilt_angle, &
               aerosol(k)%shapefrac, &
               aerosol(k)%reff, &
               rlambda, &
               csca, &
               cabs, &
               f)
         else

            aerosol_pars_mode(1) = aerosol(k)%reff
            aerosol_pars_mode(2) = aerosol(k)%veff
            aerosol_pars_mode(3) = aerosol(k)%rm(1)
            aerosol_pars_mode(4) = aerosol(k)%fim(1)
            aerosol_pars_mode(5) = aerosol(k)%aer_col
            aerosol_pars_mode(6) = aerosol(k)%shapefrac  !fraction of spheres
            call modes_calc_xs( &
               aero_lut, &
               aerosol(k)%id, &
               aerosol_pars_mode(1:6), &
               rlambda, &
               csca, &
               cabs)

         end if
         aerosol(k)%sig_ref = csca + cabs
         aerosol(k)%aer_col = aerosol(k)%tau_ref/aerosol(k)%sig_ref

!!$       if(allocated(aerosol(k)%alt_dis)) deallocate(aerosol(k)%alt_dis)
!!$       if(allocated(aerosol(k)%dalt_daer1)) deallocate(aerosol(k)%dalt_daer1)
!!$       if(allocated(aerosol(k)%dalt_daer2)) deallocate(aerosol(k)%dalt_daer2)
!!$       if(allocated(aerosol(k)%nder)) deallocate(aerosol(k)%nder)
         allocate (aerosol(k)%alt_dis(nrt), &
                   aerosol(k)%dalt_daer1(nrt), &
                   aerosol(k)%dalt_daer2(nrt), &
                   stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'GET_AEROSOL_PROPERTIES_LOGNORMAL: memory allocation problem'
            ierr = ierr_all
            goto 999
         end if
         call set_altdis(atm_rt, &
                         aerosol(k)%altid, &
                         aerosol(k)%aeralt1, &
                         aerosol(k)%aeralt2, &
                         aerosol(k)%alt_dis, &
                         aerosol(k)%dalt_daer1, &
                         aerosol(k)%dalt_daer2)
         aerosol(k)%maxd = count(aerosol(k)%alt_dis(:) > maxval(aerosol(k)%alt_dis)*nder_cut)
         allocate (aerosol(k)%nder(aerosol(k)%maxd), stat=ierr)
         if (ierr .ne. 0) then
            write (message, *) 'GET_AEROSOL_PROPERTIES_LOGNORMAL: memory allocation problem'
            ierr = ierr_all
            goto 999
         end if
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
      end do !end loop over aerosol types

      return

999   continue
      call stopretrieval(message)

   end subroutine get_aerosol_properties_lognormal

!------------------------------------------------------------------------------
!> @brief Set the height distribution for aerosol particles
!> @details Choose from 9 types of height distributions:
!! (0) Gauss, (1) exponential, (2) box, (3) Lorentz, (4) Gauss(1-par),
!! (5) double-Gauss, (6) double box, (7) one layer, (8) Gauss over restricted altitude range
!> @param[in] atm_rt            type containing atmospheric grid for radiative transfer calculations
!> @param[in] id                Identifier for type of height distribution type
!> @param[in] par1              height parameter 1 for  height distribution
!> @param[in] par2              height parameter 2 for height distribution
!> @param[out] alt(nrt)      height distribution
!> @param[out] daltd1(nrt)   derivative of alt to par1
!> @param[out] daltd(nrt)    derivative of alt to par2
!------------------------------------------------------------------------------
   subroutine set_altdis(atm_rt, id, par1, par2, alt, daltd1, daltd2)
!*** arguments
      type(atmosphere), intent(in) :: atm_rt
      integer, intent(in) :: id
      real(double), intent(in) :: par1, par2
      real(double), dimension(:), intent(out) ::  alt, daltd1, daltd2
!*** local  variables
      integer :: k, ilay, nrt
      integer, dimension(1) :: kmin
      real(double) :: GAUSS
      real(double) :: LORENTZ
      real(double) :: x, a, b, dx, norm, s1, s2, s3
      real(double) :: zcenter, zcenter1, zcenter2, zfwhm, zfwhm1, zfwhm2, zfwhm_per
      real(double), dimension(atm_rt%n) :: alt_per
      real(double), dimension(atm_rt%n) :: dalt_dzcenter, dalt_dzfwhm

!*** Function statement
      GAUSS(x, a) = DEXP(-(x*x)/(a*a/4.D0/DLOG(2.D0)))
      LORENTZ(x, a) = 1./pi*a/2./(x*x + a*a/4.)

      nrt = atm_rt%n

      if (id == 0) then
!*** Gauss(center,width)
         b = 1.D0/4.D0/DLOG(2.D0)
         zcenter = par1
         zfwhm = par2
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         alt = 0.D0
         do k = 1, nrt
            dx = atm_rt%z(k) - zcenter
            alt(k) = GAUSS(dx, zfwhm)*atm_rt%dz(k)
         end do
         norm = sum(alt)
         s1 = 0.D0
         s2 = 0.D0
         do k = 1, nrt
            s1 = s1 + (atm_rt%z(k) - zcenter)*alt(k)
            s2 = s2 + (atm_rt%z(k) - zcenter)**2*alt(k)
         end do
         do k = 1, nrt
            dalt_dzcenter(k) = 2.D0/b/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)*alt(k)/norm - &
                               2.D0/b/zfwhm/zfwhm*alt(k)*s1/norm/norm
            daltd1(k) = dalt_dzcenter(k)
            dalt_dzfwhm(k) = 2.D0/b/zfwhm/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)**2*alt(k)/norm - &
                             2.D0/b/zfwhm/zfwhm/zfwhm*alt(k)*s2/norm/norm
            daltd2(k) = dalt_dzfwhm(k)
            alt(k) = alt(k)/norm
         end do
      else if (id == 1) then
!*** Exponential(Boundary Layer height, decay)
         zcenter = par1
         zfwhm = par2
         alt = 0.D0
         do k = nrt, 1, -1
            if (atm_rt%z(k) .lt. zcenter) then
               alt(k) = 1.D0
            else
               alt(k) = alt(k + 1)*(atm_rt%p(k)/atm_rt%p(k + 1))**zfwhm
            end if
         end do
         do k = 1, nrt
            alt(k) = alt(k)*atm_rt%dz(k)
         end do
         alt = alt/sum(alt)
         zfwhm_per = zfwhm + 0.005*zfwhm
         alt_per = 0.D0
         do k = nrt, 1, -1
            if (atm_rt%z(k) .lt. zcenter) then
               alt_per(k) = 1.D0
            else
               alt_per(k) = alt_per(k + 1)*(atm_rt%p(k)/atm_rt%p(k + 1))**zfwhm_per
            end if
         end do
         do k = 1, nrt
            alt_per(k) = alt_per(k)*atm_rt%dz(k)
         end do
         alt_per = alt_per/sum(alt_per)

         do k = 1, nrt
            dalt_dzfwhm(k) = (alt_per(k) - alt(k))/(zfwhm_per - zfwhm)
         end do
         daltd1 = 0.D0
         daltd2 = dalt_dzfwhm
      else if (id == 2) then
!*** Box(center,width)
         zcenter = par1
         zfwhm = par2
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         alt = 0.D0
         do k = 1, nrt
            if (atm_rt%z(k) .le. zcenter + zfwhm/2. .and. atm_rt%z(k) .ge. zcenter - zfwhm/2.) then
               alt(k) = 1.D0
            end if
         end do
         kmin = minloc(DABS(atm_rt%z - zcenter))
         alt(kmin(1)) = 1.D0
         norm = sum(alt)
         alt = alt/norm
         daltd1 = 0.D0
         daltd2 = 0.D0
      elseif (id == 3) then
!*** Lorentz(center,width)
         zcenter = par1
         zfwhm = par2
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         alt = 0.D0
         do k = 1, nrt
            dx = atm_rt%z(k) - zcenter
            alt(k) = LORENTZ(dx, zfwhm)*atm_rt%dz(k)
         end do
         norm = sum(alt)
         s1 = 0.D0
         s2 = 0.D0
         do k = 1, nrt
            s1 = s1 + 2*(atm_rt%z(k) - zcenter)/((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + (zfwhm/2.)*(zfwhm/2.))*alt(k)
            s2 = s2 + 1./pi*(0.5*((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + zfwhm*zfwhm/2./2.) - zfwhm*zfwhm/2./2)/ &
                 (((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + zfwhm*zfwhm/2./2.)**2)*atm_rt%dz(k)
         end do
         do k = 1, nrt
            s3 = 2*(atm_rt%z(k) - zcenter)/((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + (zfwhm/2.)*(zfwhm/2.))*alt(k)
            dalt_dzcenter(k) = s3/norm - &
                               alt(k)*s1/norm/norm
            s3 = 1./pi*(0.5*((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + zfwhm*zfwhm/2./2.) - zfwhm*zfwhm/2./2.)/ &
                 (((atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter) + zfwhm*zfwhm/2./2.)**2)*atm_rt%dz(k)
            dalt_dzfwhm(k) = s3/norm - &
                             alt(k)*s2/norm/norm
         end do
         alt = alt/norm
         daltd1 = dalt_dzcenter
         daltd2 = dalt_dzfwhm
      elseif (id == 4) then
!*** Gauss(center,width(center))
         b = 1.D0/4.D0/DLOG(2.D0)
         zcenter = par1
         zfwhm = par2*GAUSS(zcenter - par2, 2*par2) + 500.D0
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         alt = 0.D0
         do k = 1, nrt
            dx = atm_rt%z(k) - zcenter
            alt(k) = GAUSS(dx, zfwhm)*atm_rt%dz(k)
         end do
         norm = sum(alt)
         s1 = 0.D0
         s2 = 0.D0
         do k = 1, nrt
            s1 = s1 + alt(k)* &
                 (2.*b*(atm_rt%z(k) - zcenter)*zfwhm*zfwhm + &
                  2.*b*(atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter)*zfwhm* &
                  par2*GAUSS(zcenter - par2, 2*par2)*(-1./2./b/par2/par2*(zcenter - par2)))/ &
                 (b*b*zfwhm*zfwhm*zfwhm*zfwhm)
            s2 = s2 + (atm_rt%z(k) - zcenter)**2*alt(k)
         end do
         do k = 1, nrt
            s3 = alt(k)* &
                 (2.*b*(atm_rt%z(k) - zcenter)*zfwhm*zfwhm + &
                  2.*b*(atm_rt%z(k) - zcenter)*(atm_rt%z(k) - zcenter)*zfwhm* &
                  par2*GAUSS(zcenter - par2, 2*par2)*(-1./2./b/par2/par2*(zcenter - par2)))/ &
                 (b*b*zfwhm*zfwhm*zfwhm*zfwhm)
            dalt_dzcenter(k) = s3/norm - &
                               alt(k)*s1/norm/norm
            dalt_dzfwhm(k) = 2.D0/b/zfwhm/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)**2*alt(k)/norm - &
                             2.D0/b/zfwhm/zfwhm/zfwhm*alt(k)*s2/norm/norm* &
                             GAUSS(zcenter - par2, 2*par2)*(1 + par2* &
                                           (-2.*b*(zcenter - par2)*2.*2.*par2*par2 - (zcenter - par2)*(zcenter - par2)*8.*b*par2)/ &
                                                            (4.*b*par2*par2)**2)
         end do
         alt = alt/norm
         daltd1 = dalt_dzcenter
         daltd2 = dalt_dzfwhm

      elseif (id == 5) then
!*** Double-Gauss(center1,center2)), width=2km
         b = 1.D0/4.D0/DLOG(2.D0)
         zcenter = par1
         zcenter2 = par2
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         if (atm_rt%z(1) .lt. zcenter2) zcenter2 = atm_rt%z(1)
         zfwhm = 2000.D0
         alt = 0.D0
         do k = 1, nrt
            dx = atm_rt%z(k) - zcenter
            alt(k) = GAUSS(dx, zfwhm)*atm_rt%dz(k)
            dx = atm_rt%z(k) - zcenter2
            alt(k) = alt(k) + GAUSS(dx, zfwhm)*atm_rt%dz(k)
         end do
         norm = sum(alt)
         s1 = 0.D0
         s2 = 0.D0
         do k = 1, nrt
            s1 = s1 + (atm_rt%z(k) - zcenter)*GAUSS(atm_rt%z(k) - zcenter, zfwhm)*atm_rt%dz(k)
            s2 = s2 + (atm_rt%z(k) - zcenter2)*GAUSS(atm_rt%z(k) - zcenter2, zfwhm)*atm_rt%dz(k)
         end do
         do k = 1, nrt
            dalt_dzcenter(k) = 2.D0/b/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)*GAUSS(atm_rt%z(k) - zcenter, zfwhm)*atm_rt%dz(k)/norm - &
                               2.D0/b/zfwhm/zfwhm*alt(k)*s1/norm/norm
            dalt_dzfwhm(k) = 2.D0/b/zfwhm/zfwhm*(atm_rt%z(k) - zcenter2)*GAUSS(atm_rt%z(k) - zcenter2, zfwhm)*atm_rt%dz(k)/norm - &
                             2.D0/b/zfwhm/zfwhm*alt(k)*s2/norm/norm
         end do
         alt = alt/norm
         daltd1 = dalt_dzcenter
         daltd2 = dalt_dzfwhm
      else if (id == 6) then
!*** Double Box(center1,center2) random width
!         CALL rdn01(2,rdn)
         zcenter1 = par1
         zfwhm1 = 1000. !rdn(1)*2000.  NOT random
         zcenter2 = par2
         zfwhm2 = 1500. !rdn(2)*2000.  NOT random

         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)

         alt = 0.D0
         do k = 1, nrt
            if (atm_rt%z(k) .le. zcenter1 + zfwhm1/2. .and. atm_rt%z(k) .ge. zcenter1 - zfwhm1/2.) then
               alt(k) = 1.D0
            end if
            if (atm_rt%z(k) .le. zcenter2 + zfwhm2/2. .and. atm_rt%z(k) .ge. zcenter2 - zfwhm2/2.) then
               alt(k) = 1.D0
            end if
         end do
         kmin = minloc(DABS(atm_rt%z - zcenter1))
         alt(kmin(1)) = 1.D0
         kmin = minloc(DABS(atm_rt%z - zcenter2))
         alt(kmin(1)) = 1.D0
         norm = sum(alt)
         alt = alt/norm
         daltd1 = 0.D0
         daltd2 = 0.D0
      else if (id == 7) then
!*** only aerosol in one rt layer
         zcenter = par1
         ilay = minval(minloc(DABS(zcenter - atm_rt%z)))
         alt = 0.D0
         alt(ilay) = 1.D0
         daltd1 = 0.D0
         daltd2 = 0.D0
      elseif (id == 8) then
!*** Gauss(center,width(center)) only at 5 different height layers
         b = 1.D0/4.D0/DLOG(2.D0)
         zcenter = par1
         zfwhm = par2
         if (atm_rt%z(1) .lt. zcenter) zcenter = atm_rt%z(1)
         ilay = minval(minloc(DABS(zcenter - atm_rt%z)))
         alt = 0.D0
         do k = max(1, ilay - 2), min(ilay + 2, nrt)
            dx = atm_rt%z(k) - zcenter
            alt(k) = GAUSS(dx, zfwhm)*atm_rt%dz(k)
         end do
         norm = sum(alt)
         s1 = 0.D0
         s2 = 0.D0
         do k = 1, nrt
            s1 = s1 + (atm_rt%z(k) - zcenter)*alt(k)
            s2 = s2 + (atm_rt%z(k) - zcenter)**2*alt(k)
         end do
         do k = 1, nrt
            dalt_dzcenter(k) = 2.D0/b/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)*alt(k)/norm - &
                               2.D0/b/zfwhm/zfwhm*alt(k)*s1/norm/norm
            daltd1(k) = dalt_dzcenter(k)
            dalt_dzfwhm(k) = 2.D0/b/zfwhm/zfwhm/zfwhm*(atm_rt%z(k) - zcenter)**2*alt(k)/norm - &
                             2.D0/b/zfwhm/zfwhm/zfwhm*alt(k)*s2/norm/norm
            daltd2(k) = dalt_dzfwhm(k)
            alt(k) = alt(k)/norm
         end do
      end if

   end subroutine set_altdis

!------------------------------------------------------------------------------

end module aerosol_input_module
