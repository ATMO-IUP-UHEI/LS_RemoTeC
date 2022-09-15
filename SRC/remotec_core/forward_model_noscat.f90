!****************************************************************************************
!*** 2011/06/20, A.Butz: - SUBROUTINE FORWARD_MODEL_LO_NOSCAT, definition of forward model
!***                       vector changed from reflectance to radiance
!***                       units, i.e. division by solar spectrum removed.
!*** 2011/06/20, A.Butz: - SUBROUTINE FORWARD_MODEL_LO_NOSCAT, cross section shift changed
!***                       to shift of the solar spectrum
!****************************************************************************************
module forward_model_noscat_module
   use header_module
   use aerosol_input_module, only: aero
   use forward_model_module, only: absorbers, window_ini, window_spectrum
   use read_xsdb_module
   use spectral_response_module
   use auxiliary_routines_module, only: U0_KASTEN_AND_YOUNG
   use ocean_fresnel_module, only: ocean_ss
   implicit none
   private

   public :: forward_model_hi_noscat

contains
!------------------------------------------------------------------------------
!> @brief Compute high-resolution model reflectance and derivatives
!! (non-scattering)
!------------------------------------------------------------------------------
   subroutine forward_model_hi_noscat( &
      XSFlag, O2Flag, Tflag, glintflag, &
      sza, iza, phi, sfwind, iwin, &
      absorb, &
      atm_xs, &
      dvair, &
      dvair_old, &
      vmr_h2o, &
      play_old, &
      win_ini, &
      win, &
      reflectance_hi, &
      deriv_dummy, &
      deriv_albedo_dummy, &
      derivT_dummy, &
      derivP_dummy, &
      ExitXSFlag, &
      ierr)
      !*** Input
      real(double), intent(in) :: sza, iza, phi, sfwind
      type(absorbers), intent(in) :: absorb
      integer, intent(in) :: iwin, XSFlag, O2Flag, TFlag, glintflag
      type(atmosphere), intent(in) :: atm_xs
      real(double), dimension(:), intent(in) :: dvair
      real(double), dimension(:), intent(in) :: dvair_old
      real(double), dimension(:), intent(in) :: play_old
      type(window_ini), dimension(:), intent(in) :: win_ini
      !*** Input/output
      type(window_spectrum), dimension(:), intent(inout) :: win
      !*** Output
      real(double), dimension(:, :), intent(out) :: reflectance_hi
      real(double), dimension(:, :, :), intent(out) :: deriv_dummy ! Reflectance derivatives wrt. absorber vmr
      real(double), dimension(:, :), intent(out) :: deriv_albedo_dummy! Reflectance derivatives wrt. albedo 0th order
      real(double), dimension(:), intent(out) :: derivT_dummy! Reflectance derivatives wrt. T
      real(double), dimension(:, :), intent(out) :: derivP_dummy! Reflectance derivatives wrt. P
      integer, intent(out) :: ExitXSFlag, ierr
      !*** local variables
      integer :: i, j, k, l, m, n, n1, n2, natm
      real(double) :: dP, dT, u0, w
      real(double), dimension(atm_xs%n) :: vmr_h2o
      real(double), dimension(win_ini(iwin)%nwave_hi) :: albedo_array   ! Integrated optical depth
      real(double), dimension(win_ini(iwin)%nwave_hi) :: tau_int   ! Integrated optical depth
      real(double), dimension(win_ini(iwin)%nwave_hi, atm_xs%n) :: dtau_P   ! Integrated optical depth
      real(double), dimension(win_ini(iwin)%nwave_hi) :: tau_int_Tper   ! Integrated optical depth
      real(double), dimension(win_ini(iwin)%nwave_hi) :: reflectance_hi_Tper   ! Transmittance perturbed by dT
      real(double), dimension(win_ini(iwin)%nwave_hi, atm_xs%n) :: amf ! Enhancement factor for light path (air mass factor)
      real(double), dimension(atm_xs%n) :: tlay_Tper ! T profile perturbed by dT
      real(double), dimension(atm_xs%n) :: play_Pper ! p profile perturbed by dP
      integer :: nlay
      real(double), dimension(:, :, :), allocatable :: cross_section      ! Molecular absorption cross sections (Dim: nwave_hi,natm,ntype)
      real(double), dimension(:, :, :), allocatable :: cross_section_Tper ! Molecular absorption cross sections, perturbed by infinitesimal DT (Dim: nwave_hi,natm,ntype)
      real(double), dimension(:, :, :), allocatable :: cross_section_Pper
      real(double), dimension(atm_xs%n) :: pcor
      character*2 :: ch
      real(double) :: uv
      real(double), dimension(:), allocatable :: bdrf_ss
      !------------------------------------------------------------------------------

      !*** Initialize error identifiers
      ierr = 0
      ExitXSFlag = 0

      !*** Initialize some local variables
      natm = atm_xs%n
      dT = 1.D0
      dP = 1.D0 + 1.D-4
      tlay_Tper = atm_xs%t + dT
      play_Pper = atm_xs%p*dP

      !*** Get water vmr in xs-layers
!    do i = 1,win_ini(iwin)%ntype
!       if(abs(win_ini(iwin)%xsdb(i)%species)==1 .or. (abs(win_ini(iwin)%xsdb(i)%species)>=100 .and. abs(win_ini(iwin)%xsdb(i)%species)<=199))then
!          vmr_h2o(:) =   win(iwin)%x_molec(:,i)/dvair! caculate vmr of water per layer
!          exit
!       endif
!    enddo

      !***Look-up absorption cross section [cm^2]
      allocate (cross_section(win_ini(iwin)%nwave_hi, natm, win_ini(iwin)%ntype))

      if (XSFlag > 0) then

         write (ch, '(i2.2)') iwin
         open (50, FILE='./CONTRL_OUT/cross_section_'//ch//'.dat')
         !         write(*,*) 'Opened absorption cross-section file: ', './CONTRL_OUT/cross_section_'//ch//'.dat', ' in remotec_core/forward_model_noscat/forward_model_hi_noscat (line 94)'
         do i = 1, win_ini(iwin)%ntype
            do k = 1, win_ini(iwin)%nwave_hi
               read (50, *) cross_section(k, :, i)
            end do
         end do
         close (50)

      else

         do i = 1, win_ini(iwin)%ntype
            pcor = 1.
            !*** This is ugly: Self broadening of H2O is non-negligible in some spectral bands.
            !*** The Lorentz self broadeing half-width gamma_self is roughly 5*gamma_air with scatter.
            !*** Since the XS lookup-table does not include a self-broadening dimension, we scale the
            !*** atmospheric pressure when reading the H2O cross sections by a factor [1+4*vmr_H2O].
            !*** This assumes that gamma_self/gamma_air=5 and that Doppler-broadening effects
            !*** is negligeable in the height ranges where the vmr_H2O is high [cf. C.Frankenberg, R. Scheepmaker]
            if (abs(win_ini(iwin)%xsdb(i)%species) == 1 .or. &
     (abs(win_ini(iwin)%xsdb(i)%species)>=100 .and. abs(win_ini(iwin)%xsdb(i)%species)<=199).and. win_ini(iwin)%type_xsdb(1)/=4)then
               pcor = 1.+4.*win(iwin)%x_molec(:, i)/dvair
               do k = 1, natm
                  if (pcor(k)*atm_xs%p(k) > 1100.) pcor(k) = 1100./atm_xs%p(k)
               end do
            else
               pcor = 1.
            end if

            do n = 1, natm
               call get_xs( &
                  win_ini(iwin)%xsdb(i), &
                  vmr_h2o(n), &
                  atm_xs%p(n)*pcor(n), &
                  atm_xs%t(n), &
                  cross_section(:, n, i), &
                  ExitXSFlag, ierr)
               if (ExitXSFlag .ne. 0 .or. ierr .ne. 0) return
            end do
         end do

         if (TFlag == 1) then
            allocate (cross_section_Tper(win_ini(iwin)%nwave_hi, natm, win_ini(iwin)%ntype))
            do i = 1, win_ini(iwin)%ntype
               do n = 1, natm
                  call get_xs( &
                     win_ini(iwin)%xsdb(i), &
                     vmr_h2o(n), &
                     atm_xs%p(n), &
                     tlay_Tper(n), &
                     cross_section_Tper(:, n, i), &
                     ExitXSFlag, ierr)
                  if (ExitXSFlag .ne. 0 .or. ierr .ne. 0) return
               end do
            end do
         end if

         if (O2Flag == 1) then
            allocate (cross_section_Pper(win_ini(iwin)%nwave_hi, natm, win_ini(iwin)%ntype))
            do i = 1, win_ini(iwin)%ntype
               do n = 1, natm
                  call get_xs( &
                     win_ini(iwin)%xsdb(i), &
                     vmr_h2o(n), &
                     play_Pper(n), &
                     atm_xs%t(n), &
                     cross_section_Pper(:, n, i), &
                     ExitXSFlag, ierr)
                  if (ExitXSFlag .ne. 0 .or. ierr .ne. 0) return
               end do
            end do
         end if

         if (XSFlag < 0) then

            write (ch, '(i2.2)') iwin
            open (50, file='CONTRL_OUT/cross_section_'//ch//'.dat')
            !              write(*,*) 'Opened output cross-section file: ', './CONTRL_OUT/cross_section_'//ch//'.dat', ' in remotec_core/forward_model_noscat/forward_model_hi_noscat (line 166)'
            do i = 1, win_ini(iwin)%ntype
               do k = 1, win_ini(iwin)%nwave_hi
                  write (50, '(300(E17.10,X))') cross_section(k, :, i)
               end do
            end do
            close (50)

         end if
      end if

      !***Calculate air mass factor
      call u0_kasten_and_young(dble(sza), u0)
      !     u0 = cos(DBLE(sza)/180.*Pi)
      amf = 1./cos(dble(iza)/180.*pi) + 1./u0

      !***Calculate reflectance spectrum
      tau_int = 0.
      do n = 1, natm
         do i = 1, win_ini(iwin)%ntype
            tau_int(:) = tau_int(:) + cross_section(:, n, i)*win(iwin)%x_molec(n, i)*amf(:, n)
         end do
      end do
      !***If ocean-glint: bidirectional reflection distr. functions calculated from the ocean model
      if (glintflag == 1) then
         uv = DCOS(dble(iza)/180.*pi)
         allocate (bdrf_ss(nstokes))
         call ocean_ss(u0, uv, phi, sfwind, bdrf_ss)
         win(iwin)%albedo(1) = win(iwin)%albedo(1) + bdrf_ss(1)
         deallocate (bdrf_ss)
      end if

      albedo_array = 0.D0
      do i = 1, win_ini(iwin)%albflag
         albedo_array(:) = albedo_array(:) + &
                           (win_ini(iwin)%wavelength_hi(:) - win_ini(iwin)%wavelength_hi(1)) &
                           **(i - 1)*win(iwin)%albedo(i)
      end do
      !*** HH: added albflag=0 (albedo from satellite data) for creating synthetic measurements
      if (win_ini(iwin)%albflag == 0) albedo_array(:) = win(iwin)%albedo(1)

      reflectance_hi = 0.D0
!      reflectance_hi(:) = albedo_array(:)*u0/pi*exp(-tau_int(:))
      reflectance_hi(:, 1) = albedo_array(:)*u0/pi*dexp(-tau_int(:))

      !***Calculate derivatives of reflectance wrt. absorber vmr
      nlay = size(deriv_dummy(1, :, 1))
      deriv_dummy = 0.D0
      l = natm/nlay
      do j = 1, absorb%ntype_target
         do i = 1, win_ini(iwin)%ntype
            if (win_ini(iwin)%type_x(i) == absorb%type_x_target(j)) then
               do n = 0, natm - 1
                  m = int(n/l) + 1
                  n1 = (m - 1)*l + 1
                  n2 = m*l
                  if (sum(win(iwin)%x_molec(n1:n2, i)) > 1.d-12) then
                     w = win(iwin)%x_molec(n + 1, i)/sum(win(iwin)%x_molec(n1:n2, i))
                  else
                     w = 1/l
                  end if
                  deriv_dummy(:, m, i) = deriv_dummy(:, m, i) - &
                                         reflectance_hi(:, 1)*amf(:, n + 1)*cross_section(:, n + 1, i)*w
!                       reflectance_hi(:)*amf(:,n+1)*cross_section(:,n+1,i)*w
               end do
            end if
         end do
      end do

      do j = 1, absorb%ntype_global
         do i = 1, win_ini(iwin)%ntype
            if (win_ini(iwin)%type_x(i) == absorb%type_x_global(j)) then
               do n = 0, natm - 1
                  if (sum(win(iwin)%x_molec(:, i)) > 1.d-12) then
                     w = win(iwin)%x_molec(n + 1, i)/sum(win(iwin)%x_molec(:, i))
                  else
                     w = 1/natm
                  end if
                  deriv_dummy(:, 1, i) = deriv_dummy(:, 1, i) - &
                                         reflectance_hi(:, 1)*amf(:, n + 1)*cross_section(:, n + 1, i)*w
!                       reflectance_hi(:)*amf(:,n+1)*cross_section(:,n+1,i)*w
               end do
               forall (l=2:nlay) deriv_dummy(:, l, i) = deriv_dummy(:, 1, i)
            end if
         end do
      end do

      !***Calculate derivatives with respect to albedo
      deriv_albedo_dummy = 0.D0
      do i = 1, win_ini(iwin)%albflag
         deriv_albedo_dummy(:, i) = u0/pi* &
                                    exp(-tau_int(:))*((win_ini(iwin)%wavelength_hi(:) - &
                                                       win_ini(iwin)%wavelength_hi(1))**(i - 1))
      end do

      !***Calculate derivatives wrt temperature
      derivT_dummy = 0.D0
      if (TFlag == 1) then
         tau_int_Tper = 0.D0
         do n = 1, natm
            do i = 1, win_ini(iwin)%ntype
               tau_int_Tper(:) = tau_int_Tper(:) + cross_section_Tper(:, n, i)*win(iwin)%x_molec(n, i)*amf(:, n)
            end do
         end do
         reflectance_hi_Tper(:) = albedo_array*u0/pi*exp(-tau_int_Tper(:))
!         derivT_dummy(:) = (reflectance_hi_Tper - reflectance_hi(:))/dT
         derivT_dummy(:) = (reflectance_hi_Tper - reflectance_hi(:, 1))/dT
      end if

      !***Calculate derivatives wrt pressure/O2 column
      if (O2Flag == 1) then
         derivP_dummy = 0.D0
         l = natm/nlay
         dtau_P = 0.D0
         do n = 1, natm
            do i = 1, win_ini(iwin)%ntype
               dtau_P(:, n) = dtau_P(:, n) + &
                              (cross_section_Pper(:, n, i) - cross_section(:, n, i))/ &
                              (play_Pper(n) - atm_xs%p(n))* &
                              win(iwin)%x_molec(n, i)*amf(:, n)*play_old(n)
            end do
         end do
         do j = 1, absorb%ntype_target
            if (absorb%type_x_target(j) == 7) then
               do n = 0, natm - 1
                  m = int(n/l) + 1
                  n1 = (m - 1)*l + 1
                  n2 = m*l
                  w = dvair(n + 1)/sum(dvair(n1:n2))
                  derivP_dummy(:, m) = derivP_dummy(:, m) - &
                                       albedo_array(:)*u0/pi* &
                                       exp(-tau_int(:))*dtau_P(:, n + 1)/(dvair_old(n + 1)*RELO2)*w
               end do
            end if
            exit   ! target absorber O2 found, no need to continue do loop
         end do
         do j = 1, absorb%ntype_global
            if (absorb%type_x_global(j) == 7) then
               do n = 1, natm
                  w = dvair(n)/sum(dvair(:))
                  derivP_dummy(:, 1) = derivP_dummy(:, 1) - &
                                       albedo_array(:)*u0/pi* &
                                       exp(-tau_int(:))*dtau_P(:, n)/(dvair_old(n)*RELO2)*w
               end do
               forall (l=2:nlay) derivP_dummy(:, l) = derivP_dummy(:, 1)
            end if
            exit   ! global absorber O2 found, no need to continue do loop
         end do
         do i = 1, win_ini(iwin)%ntype
            if (win_ini(iwin)%type_x(i) == 7) then
               derivP_dummy(:, :) = derivP_dummy(:, :) + deriv_dummy(:, :, i)
            end if
            exit   ! absorber O2 found, no need to continue do loop
         end do
      end if !O2Flag

   end subroutine forward_model_hi_noscat

!------------------------------------------------------------------------------
end module forward_model_noscat_module
