!------------------------------------------------------------------------------
!> @brief Radiative transfer model
!> @details The optical properties are calculated based on LUTs of 
!! Dubovik et al, 2004 for mixture of spheroids and spheres.
!! Phase function matrix for single scattering is also directly extracted from 
!! LUT instead of being calculated from expansion coefficients. 
!! A separate routine is called if only cross sections
!! are needed for aerosol/cirrus and not the phase function (matrix). 
!------------------------------------------------------------------------------
module rad_trans_module
  use header_module
  use auxiliary_routines_module
  use read_settings_module, only: window_ini, settings_flags, file_paths, read_settings, read_win_xsdb
  use optic_input_module, only: Mie_lut, cirrus_table, read_aerosol_netcdf, read_cirrus_netcdf, &
       aero, set_altdis, get_aerosol_properties_lognormal, aero_opt, &
       optic_all_der, optic_all_msc, optic_all_ssc
  use optic_molec_module, only: calculate_optic_mol_prop
  use tau_grid_module, only: interpolate_lbl_2species, interpolate_lbl_1species, set_tau_grid
  use rad_intf_module
  use ocean_fresnel_module, only: ocean_ss, ocean_fou, ocean_fou_lintranV2 
  use header_gsd
  use lintran_interface
  implicit none
  private

  !*** types
  public :: Mie_lut, cirrus_table, aero, derivatives, window_ini, settings_flags, file_paths

  !*** procedures
  public :: rad_trans_intf, read_aerosol_netcdf, read_cirrus_netcdf, set_altdis, get_aerosol_properties_lognormal, &
       read_settings, read_win_xsdb
  private :: rad_trans_efficient, combine_der, combine_der_ss, nakajima_correction, &
       u0_kasten_and_young
  !------------------------------------------------------------------------------
  type :: derivatives
     !> Derivatives wrt. absorber partial column (dim: nwave_hi, nrt, ntype, nstokes)
     real(double), dimension(:,:,:,:), allocatable :: densmol 
     !> Derivatives wrt. temperature offset (dim: nwave_hi, nstokes)
     real(double), dimension(:,:), allocatable :: T          
     !> Derivatives wrt. oxygen column (dim: nwave_hi, nrt, nstokes)
     real(double), dimension(:,:,:), allocatable :: p         
     !> Derivatives wrt. albedo (nwave_hi,order,nstokes)
     real(double), dimension(:,:,:), allocatable :: alb         
     !> Derivatives of intensity wrt. aerosol parameters (dim: nwave_hi, npar, ntype_aer, nstokes)
     real(double), dimension(:,:,:,:), allocatable :: aerosol 
     !> Derivatives wrt. fluorescence paramaters (dim: nwave_hi, 2)
     real(double), dimension(:,:), allocatable :: Fs 
  end type derivatives

contains
  !------------------------------------------------------------------------------
  !> @brief Interface to radiative transfer model
  !------------------------------------------------------------------------------
  subroutine rad_trans_intf( & ! {{{
       aero_lut, cirrus_lut, &
       flag, glintflag, &
       iband, sza, theta, phi, wspeed, &
       aerosol, minaotflag, mincotflag, nrt, nder, &
       atm_xs, dvair, dvair_old, vmr_h2o, play_old, &
       win_ini, albedo, x_molec, wavelength, &
       rint_fine, deriv_rt, taua_tot, ot, cot, ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
    !*** Input
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    type(settings_flags), intent(in) :: flag
    integer, intent(in) :: glintflag, iband
    real(double), intent(in) :: sza, theta, phi, wspeed
    type(aero), dimension(:), intent(in) :: aerosol   
    integer, intent(in) :: minaotflag, mincotflag, nrt
    integer, dimension(:), intent(in) :: nder 
    type(atmosphere), intent(in) :: atm_xs
    real(double), dimension(:), intent(in) :: dvair
    real(double), dimension(:), intent(in) :: dvair_old
    real(double), dimension(:), intent(in) :: vmr_h2o    
    real(double), dimension(:), intent(in) :: play_old 
    type(window_ini), intent(in) :: win_ini
    real(double), dimension(:), intent(in) :: wavelength
    real(double), dimension(:), intent(in) :: albedo
    real(double), dimension(:,:), intent(in) :: x_molec    
    real(double), dimension(:), allocatable :: taua_tot  
    logical, intent(in) :: deriv_flag
    !*** Output
    real(double), intent(out) :: ot, cot
    integer, intent(out) :: ExitXSFlag, MaxOTFlag, ierr
    real(double), dimension(:,:), allocatable, intent(out) :: rint_fine
    type(derivatives), intent(out) :: deriv_rt
    !*** local variables
    integer :: i, l, nmol, ntype_aer, nwave
    real(double), dimension(:), allocatable :: albedo_array 
    real(double) :: starttime, endtime
    character(stringlen) :: message

    !------------------------------------------------------------------------------

    !*** Initialize error identifier 
    ierr = 0

    !*** Transfer window specific variables to RTM variables
    nwave = win_ini%nwave_hi
    nmol  = win_ini%ntype

    ntype_aer = size(aerosol)

    if(allocated(rint_fine)) deallocate(rint_fine)  
    if(allocated(deriv_rt%densmol)) deallocate(deriv_rt%densmol)
    if(allocated(deriv_rt%alb)) deallocate(deriv_rt%alb)
    if(allocated(deriv_rt%T)) deallocate(deriv_rt%T)
    if(allocated(deriv_rt%P)) deallocate(deriv_rt%P)
    if(allocated(deriv_rt%aerosol)) deallocate(deriv_rt%aerosol)
    allocate( &
         rint_fine(nwave, nstokes), &
         deriv_rt%densmol(nwave, nrt, nmol, nstokes), &
         deriv_rt%alb(nwave, 1, nstokes), & 
         deriv_rt%T(nwave, nstokes), &
         deriv_rt%P(nwave, nrt, nstokes), &
         deriv_rt%aerosol(nwave, npar, ntype_aer, nstokes), &
         albedo_array(nwave), stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_INTF: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    if(win_ini%albflag > 0) then
       do i = 1, nwave
          albedo_array(i) = 0.D0
          do l = 1, win_ini%albflag
             albedo_array(i) = albedo_array(i) + albedo(l)*((wavelength(i)-wavelength(1))**(l-1))
          enddo
       enddo
    elseif(win_ini%albflag == 0) then
       do i = 1, nwave 
          albedo_array(i) = albedo(1)
       enddo
    endif

    if (flag%output >= 2 ) call cpu_time(starttime)

    !*** Call RTM efficient/line-by-line
    if (win_ini%ntau > 0) then
       call rad_trans_efficient( & 
            aero_lut, cirrus_lut, &
            flag, glintflag, &
            iband, &
            wavelength, &
            albedo_array, &
            sza, &
            theta, &
            phi, &
            wspeed, &
            aerosol, &
            minaotflag, mincotflag, &
            nrt, &
            nder, &
            atm_xs, &
            dvair, &
            dvair_old, &
            vmr_h2o, &
            play_old, &
            win_ini, &
            x_molec, &
            rint_fine, &
            deriv_rt, &
            taua_tot, &
            ot, cot,  &
            ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
    else
       call rad_trans_lbl( & 
            aero_lut, cirrus_lut, &
            flag, glintflag, &
            iband, &
            wavelength, &
            albedo_array, &
            sza, &
            theta, &
            phi, &
            wspeed, &
            aerosol, &
            minaotflag, mincotflag, &
            nrt, &
            nder, &
            atm_xs, &
            dvair, &
            dvair_old, &
            vmr_h2o, &
            play_old, &
            win_ini, &
            x_molec, &
            rint_fine, &
            deriv_rt, &
            taua_tot, &
            ot, cot,  &
            ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
    endif

    if (flag%output >= 2 ) then
       call cpu_time(endtime)
       write(message,'(a, ES12.3)') 'Time needed for RTM:', endtime-starttime
       call writelog(message, 3)
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine rad_trans_intf ! }}}

  !------------------------------------------------------------------------------    
  !> Efficient calculation of radiative transfer using k-binning for multiple scattering,
  !! line by line calculation for single scattering
  !------------------------------------------------------------------------------
  subroutine rad_trans_efficient( & ! {{{
       aero_lut, cirrus_lut, &            
       flag, glintflag, &
       iband, wavelength, albedo_array, sza, theta, phi, wspeed, &
       aerosol, minaotflag, mincotflag, nrt, nder, atm_xs, &
       dvair, dvair_old, vmr_h2o, play_old, &
       win_ini, x_molec, &
       rint_fine, deriv_rt, taua_tot, ot, cot, &
       ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
    !*** Input
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    type(settings_flags), intent(in) :: flag
    integer, intent(in) :: glintflag, iband
    real(double), dimension(:), intent(in) :: wavelength 
    real(double), dimension(:), intent(in) :: albedo_array 
    real(double), intent(in) :: sza, theta, phi, wspeed
    type(aero), dimension(:), intent(in) :: aerosol 
    integer, intent(in) :: minaotflag, mincotflag, nrt
    integer, dimension(:), intent(in) :: nder     
    type(atmosphere), intent(in) :: atm_xs
    real(double), dimension(:), intent(in) :: dvair
    real(double), dimension(:), intent(in) :: dvair_old
    real(double), dimension(:), intent(in) :: vmr_h2o   
    real(double), dimension(:), intent(in) :: play_old 
    type(window_ini), intent(in) :: win_ini
    real(double), dimension(:,:), intent(in) :: x_molec  
    logical, intent(in) :: deriv_flag
    !*** Input/output
    real(double), dimension(:,:), intent(inout) :: rint_fine
    type(derivatives), intent(inout) :: deriv_rt
    !*** Output
    integer, intent(out) :: ExitXSFlag, MaxOTFlag, ierr
    real(double), intent(out) :: ot, cot
    !*** local variables      
    integer, dimension(nrt) :: ncoefs   
    real(double), dimension(nstokes) :: dalb, rint, dalb_ss, rint_ss
    real(double), dimension(nstokes, nrt) :: dtaua, dtaus, dtaus_ss, dtaua_ss
    real(double), dimension(nstokes, npar) :: k_tmp_phase, k_tmp_taua, k_tmp_taus 
    real(double), dimension(:, :, :, :), allocatable:: dphase 
    real(double), dimension(:,:), allocatable :: rint_fine_ss
    real(double), dimension(:,:), allocatable :: rint_fine_ms
    real(double), dimension(:,:), allocatable :: rint_fine_temp
    real(double), dimension(:,:,:), allocatable :: rint_tau
    real(double), dimension(:,:,:,:,:), allocatable :: kmat_aerosol_phase  		
    real(double), dimension(:,:,:,:,:), allocatable :: kmat_aerosol_taua 	
    real(double), dimension(:,:,:,:,:), allocatable :: kmat_aerosol_taus 	
    real(double), dimension(:,:,:,:), allocatable :: kmat_densmol_ss
    real(double), dimension(:,:), allocatable :: kmat_alb_ss
    real(double), dimension(:,:,:), allocatable :: kmat_alb_tau
    real(double), dimension(:,:,:,:), allocatable :: dtaua_tau
    real(double), dimension(:,:,:,:), allocatable :: kmat_taus_tau
    real(double), dimension(:,:,:), allocatable :: kmat_taus_fine
    real(double), dimension(:,:,:), allocatable :: kmat_taua_fine
    real(double), dimension(:,:,:,:), allocatable :: kmat_aerosol_ss
    real(double), dimension(:, :), allocatable:: dz_ss
    real(double) :: bdrf_ms(nstokes, nstokes, mxhalf+2, mxhalf+2, 0:maxstr-1)
    real(double) :: bdrf_ss(nstokes), bdrf_ss_cst(nstokes)
    real(double), dimension(:), allocatable :: taua_tot
    real(double), dimension(:), allocatable :: taus_tot
    real(double), dimension(:,:), allocatable :: taua_tot_species       
    integer :: i, j, k, l, n, imid, i1, i2, iwave, i_st, &
         itau, itau_2nd,  &
         imol, nmol, &
         imol_1st, imol_2nd, ntau, ntau_2nd_max, ntype_aer, nwave
    real(double) :: wave_mid, albedo_av, u0, uv, scat_angle 
    type(aero_opt), dimension(:), allocatable :: aerosol_optic
    real(double), dimension(:,:), allocatable :: taua_aer
    real(double), dimension(:,:), allocatable :: taus_aer
    real(double), dimension(nrt) :: taua 
    real(double), dimension(nrt) :: taus
    ! Phase matrix coefficients for aerosol scattering
    real(double), dimension(nstokes,nstokes,0:maxleg, nrt) :: plmom_aer
    integer, dimension(nrt) :: maxcoefs       ! Number of coefficients for phase matrix (nrt)
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt) :: plmom_in  !Phase matrix coefficients
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt) :: plmom     ! Phase matrix coefficients after Nakajima correction
    real(double), dimension(nstokes, nrt) :: z_ss          ! Z matrix for single scattering (nstokes, nrt)
    real(double), dimension(4, 4, nrt, size(wavelength)) :: zss_lin! Z matrix for Lintran
    real(double), dimension(:,:,:), allocatable :: z_ss_int    ! Z matrix for single scattering, interpolated for wavelength (nstokes, nrt, nwave_hi)
    !*** K-binning
    integer, dimension(:), allocatable :: ntau_2nd   
    real(double), dimension(:), allocatable :: tau_grid
    real(double), dimension(:,:), allocatable :: tau_grid_2nd
    real(double), dimension(:,:,:), allocatable :: tau_grid_arr 
    !*** Optical properties of molecules
    real(double), dimension(:,:), allocatable :: taua_mol, taus_mol
    real(double), dimension(:,:,:), allocatable ::  taua_mol_species ! Molecular absorption optical thickness per type (nwave_hi, nrt, ntype)
    real(double), dimension(:,:,:), allocatable :: cross_mol         ! Molecular absorption cross section (nwave_hi, nrt, ntype) 
    real(double), dimension(:), allocatable :: csray                 ! Rayleigh scattering cross section (nwave_hi)
    real(double), dimension(:,:), allocatable :: dtaua_T             ! Derivative of absorption optical thickness wrt temperature (nwave_hi, nrt) 
    real(double), dimension(:,:), allocatable :: dtaua_P             ! Derivative of absorption optical thickness wrt pressure (nwave_hi, nrt, ntype) 
    real(double), dimension(nstokes, nstokes, 0:2, nrt) :: plmom_ray ! Phase matrix coefficients for Rayleigh scattering
    integer :: maxd
    type(phase_mat_type),dimension(0:maxstr) :: pm_in
    type(gauss_quad_type) :: gauss_quad
    type(gsf_type) :: gsf
    character(stringlen) :: message
    real(double), dimension(:,:,:,:), allocatable :: tmp_daerosol_phase, & 
         tmp_daerosol_taua, &  
         tmp_daerosol_taus 
    !*** The Lintran input instances
    logical :: cloudflag ! Flag for using cloudy settings
    type(lintran_atmosphere) :: atm
    type(lintran_derivatives) :: drv
    type(lintran_settings) :: set
    type(lintran_class) :: lintran_instance
    integer :: lintran_stat    
    logical, dimension(3) :: execution ! Flags for single, double and multi-scattering, in that order.
    real(double), dimension(nstokes,1) :: rintms, rintss ! Main output
    real(double) :: surf_emi
    real(double), parameter :: degrees = acos(0.0) / 90.0
    real, parameter :: SMALL = 1.D-60
    real(double), dimension(nstokes) :: rint_per
    !*** Indices to convert scattering matrix to array following LINTRAN2 convention
    integer, dimension(4) :: index_ist = (/1,1,2,3 /)
    integer, dimension(4) :: index_jst = (/1,2,2,3 /)   
    !--------------------------------------------------------------------------------

    !*** Initialize error identifier 
    ierr = 0
    MaxOTFlag = 0

    !**** Initialization and Allocation
    uv = DCOS(theta/180.*pi)
    call u0_kasten_and_young(sza, u0)
    scat_angle = DACOS(&
         -u0*uv + &
         DSQRT(1.D0-uv**2)* &
         DSQRT(1.D0-u0**2)* &
         DCOS(phi/180.D0*PI))      
    ntau = win_ini%ntau
    ntau_2nd_max = win_ini%ntau_2nd_max
    nwave = size(wavelength)
    nmol = win_ini%ntype   ! number of absorbers      
    imid = (nwave+1)/2
    wave_mid = wavelength(imid)
    ntype_aer = size(aerosol)
    maxd = size(nder)
    allocate(&
         dz_ss(nstokes, maxd), &
         dphase(nstokes, 0:maxleg, nper, maxd), &
         kmat_densmol_ss(nwave, nrt, nmol, nstokes),&
         kmat_alb_ss(nwave, nstokes),&
         rint_tau(ntau, ntau_2nd_max, nstokes),&
         rint_fine_ms(nwave, nstokes),&
         rint_fine_ss(nwave, nstokes),&
         rint_fine_temp(nwave, nstokes),&
         dtaua_tau(ntau, ntau_2nd_max, nrt, nstokes),&
         kmat_alb_tau(ntau, ntau_2nd_max, nstokes),&
         kmat_taus_tau(ntau, ntau_2nd_max, nrt, nstokes),&
         kmat_taus_fine(nwave, nrt, nstokes),&
         kmat_taua_fine(nwave, nrt, nstokes),&
         kmat_aerosol_ss(nwave, npar, ntype_aer, nstokes),&
         kmat_aerosol_phase(ntype_aer, ntau, ntau_2nd_max, npar, nstokes),& 	
         kmat_aerosol_taua(ntype_aer, ntau, ntau_2nd_max, npar, nstokes),&	
         kmat_aerosol_taus(ntype_aer, ntau, ntau_2nd_max, npar, nstokes), &
         tmp_daerosol_phase(nwave, npar, ntype_aer, nstokes),& 	
         tmp_daerosol_taua(nwave, npar, ntype_aer, nstokes),&	
         tmp_daerosol_taus(nwave, npar, ntype_aer, nstokes), &
         stat=ierr)	

    kmat_aerosol_ss = 0.d0
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    !**** Select taua_tot binning points for multiple scattering calculations
    if(allocated(taua_tot)) deallocate(taua_tot)
    if(allocated(taus_tot)) deallocate(taus_tot)
    if(allocated(taua_tot_species)) deallocate(taua_tot_species)
    allocate(&
         taua_tot(nwave), &
         taus_tot(nwave),&
         taua_tot_species(nwave,win_ini%ntype), stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    !*** Calculate molecular optical quantities
    call calculate_optic_mol_prop( &
         iband, flag%xs, flag%O2, flag%temp, &  !Input
         wavelength, &
         nrt, &
         atm_xs, &
         win_ini%xsdb, &
         win_ini%type_xsdb, &
         x_molec, &
         dvair, &
         dvair_old, &
         vmr_h2o, &
         play_old, &
         taua_mol, &                            !Output
         taus_mol, &
         taua_mol_species, &
         cross_mol, &
         csray, &
         dtaua_T, &
         dtaua_P, &
         plmom_ray, &
         ExitXSFlag, ierr)
    if (ExitXSFlag .ne. 0 .or. ierr .ne. 0) return

    !*** Calculate aerosol optical quantities, only for middle of absorption band
    allocate(aerosol_optic(ntype_aer), &
         taus_aer(nwave, nrt),&
         taua_aer(nwave, nrt), stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    !*** calculate aerosol scattering and absorption optical thickness and phase function for multiple scattering 
    call optic_all_msc( &
         aero_lut, cirrus_lut, &
         iband, wavelength, &
         aerosol, nder, nrt, aerosol_optic, &
         taua_aer, taus_aer, ot, cot, &
         plmom_aer, maxcoefs, ierr)

    ot = ot + cot ! total optical thickness

    if (flag%output >= 2) then   
       write(message,'(a,2ES12.3,a,ES12.3)') 'total OT and COT: ', ot, cot , ' at wavelength ', wave_mid
       call writelog(message, 3)
    endif

    !*** Error checks
    !*** Check for NAN values:
    if(ot .ne. ot)then
       MaxOTFlag = 1
       return
    endif
    !*** Check for too high OTs
    if(ot > 5.d0 .and. .false.) then ! Hack: Turn off total optical depth check
       MaxOTFlag = 1
       return
    endif

    ! Flag for using cloudy settings for Lintran
    cloudflag = ot .gt. 3.0
    !*** Interface aerosol and molecular quantities to RTM, only for middle of absorption band
    !*** taus,taua: total scattering/absorption optical thickness (nmb. of layers)
    !*** dtaus_all,dtaua_all: derivative of total scattering/absorption optical thickness wrt aerosol parameters
    !*** plmom_in, plmom_ss: phase matrix coefficients
    !*** dphase_all: derivatives of coefficients of aerosol scattering matrix wrt aerosol parameters

    call optic_all_der( &
         imid, &
         aerosol_optic, aerosol, nder, nrt, maxcoefs, &
         taua_aer, taus_aer, plmom_aer, &
         taua_mol, taus_mol, plmom_ray, &
         taua, taus, plmom_in, &
         ierr)

    if ( flag%rtm==1) then
       !*** Initilization of RTM
       call perturbation_init(&
            u0,&
            theta,&
            gauss_quad,&
            gsf, ierr)
       if (ierr .ne. 0) return

       !**** Nakajima correction for strongly forward peaked phase functions
       call nakajima_correction(plmom_in, plmom, nder, maxd, nrt)
       !***Nakajima switch: off = (plmom=plmom_in), on = (plmom!=plmom_in)
       !$$         plmom = plmom_in
       !**** Calculate Fourier coefficients of phase matrix
       forall(k=1:nrt) ncoefs(k) = min(maxcoefs(k), maxstr-1)
       !*** calculate fourier coefficients of phase matrix -- Nakajima/Delta-M 
       call pm_coeff(&
            nrt, &
            nstokes,&
            ncoefs,&
            plmom,&
            gsf,&
            pm_in, &
            ierr)
       if (ierr .ne. 0) return
    elseif(flag%rtm==2 .or.flag%rtm==3) then 
       call lintran_init(nstokes, MAXSTR, nrt, 1, .true., lintran_instance, lintran_stat)
       ! Convert the azimuthal difference to the convention of de Haan et al. (1987).
       ! That is what Lintran 2 desires.
       call lintran_provide(real(u0),(/real(uv)/),(/-real(phi*degrees)/),.false.,lintran_instance, lintran_stat)
    endif

    do l = 1, nwave 
       !  cut k-binning grid at max value tatot
       taua_tot(l) = min(sum(taua_mol(l,:)),tatot)
       taus_tot(l) = sum(taus_mol(l,:)) + sum(taus_aer(l,:))
    enddo

    do imol = 1, nmol     
       do l = 1, nwave 
          !  cut k-binning grid at max value tatot 
          taua_tot_species(l,imol) =  min(sum(taua_mol_species(l,:,imol)),tatot)
       enddo
    enddo

    if(nmol.eq.1 .or. ntau_2nd_max .eq. 1) then
       call set_tau_grid(nwave, nrt, ntau, taus, taua_tot, taua_mol, tau_grid, tau_grid_arr, ierr)
       if (ierr .ne. 0) return
       allocate(ntau_2nd(ntau), stat=ierr)
       if (ierr .ne. 0) then
          write(message,*) 'RAD_TRANS_EFFICIENT: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       ntau_2nd = 1
       !*** If number of absorbers eq. 2: 2 grids in taua_tot are needed
       !*** imol_1st should be target absorber, in the retrieval module absorbers are ordered as given in 'INI/retrieval.ini'
    elseif(nmol.eq.2) then
       imol_1st = 1
       imol_2nd = 2
       call set_tau_grid( & 
            nwave, nrt, imol_1st, imol_2nd, ntau, ntau_2nd_max, taus, taua_tot_species, taua_mol, taua_mol_species, &
            ntau_2nd, tau_grid, tau_grid_2nd, tau_grid_arr, ierr)     
       !***If number of absorbers > 2: attribute indices > 2 to 2nd absorber    
    elseif(nmol.gt.2) then
       imol_1st = 1
       imol_2nd = 2
       do imol = 3, nmol
          do l = 1, nwave 
             taua_tot_species(l, 2) = taua_tot_species(l, 2) + sum(taua_mol_species(l, :, imol))
             taua_mol_species(l,:,2) = taua_mol_species(l,:,2) + taua_mol_species(l, :, imol)
             !IF(taua_tot_species(l,2)>100.D0)taua_tot_species(l,2)=100.D0
          enddo
       enddo
       call set_tau_grid(&
            nwave, nrt, imol_1st, imol_2nd, ntau, ntau_2nd_max, taus, taua_tot_species, taua_mol, taua_mol_species, &
            ntau_2nd, tau_grid, tau_grid_2nd, tau_grid_arr, ierr)
    endif


    !**** Albedo is set to the average albedo of the window    
    surf_emi = 0D0  
    albedo_av = sum(albedo_array)/dble(nwave)

    !*** OCEAN glint
    bdrf_ms = 0.D0
    bdrf_ss = 0.D0

    if (glintflag == 1) then
       !*** Bidirectional reflection distr. functions calculated from the ocean model   
       if ( flag%rtm == 1)then
          call ocean_fou(u0, uv, wspeed, gauss_quad, bdrf_ms)
       elseif (flag%rtm==2 .or. flag%rtm==3) then   
          call ocean_fou_lintranV2(u0, uv, wspeed, lintran_instance%grd%mu, bdrf_ms)  
       endif
       call ocean_ss(u0, uv, phi, wspeed, bdrf_ss_cst)
       bdrf_ss(1) =  albedo_av + bdrf_ss_cst(1)
    else
       bdrf_ms(1,1,:,:,0) = albedo_av
       bdrf_ss(1) =  albedo_av 
    endif

    !**** Interface rad/trans variables to lintranV2 internals
    if ( flag%rtm==2 .or. flag%rtm==3) then
       call lintran_allocate( nrt, maxd, atm, drv, set, ierr)      
       call lintranv2_assign(nstokes, nrt, maxd, nder, cloudflag, maxcoefs, plmom_in, bdrf_ms, bdrf_ss, &
            zss_lin(:,:,:,1), surf_emi, taua, taus, atm, drv, set ) ! for multiple scattering zss_lin does nothing	
       execution = (/.false.,.true.,.true./) !(/ssc, dsc, msc/) -> multiple scattering solution
    endif

    !**** MULTIPLE SCATTERING

    !**** Loop over all selected tau values      
    do itau = 1, ntau   
       do itau_2nd = 1, ntau_2nd(itau)            
          do i = 1, nrt
             taua(i) = tau_grid_arr(itau, itau_2nd, i) + taua_aer(imid, i)
             taus(i) = taus_mol(imid, i) + taus_aer(imid, i)
          enddo
          !******* Call forward-adjoint perturbation RTM for multiple scattering
          !******* output @ itau,itau_2nd: rint=reflectance,
          !******* dtaua=deriv. wrt. absorption optical depth,
          !******* dtaus=deriv. wrt. scattering optical depth, 
          !******* dphase=deriv. wrt. phase matrix parameters,
          !******* dalb=deriv. wrt. albedo  
          rint = 0.D0
          dtaua = 0.D0
          dtaus = 0.D0
          dphase = 0.D0
          dalb = 0.D0

          if ( flag%rtm==1 )then

             call fwd_adj_perturbation(&
                  nrt, &
                  taua,&
                  taus,&
                  plmom_in, &
                  pm_in,&
                  gauss_quad,&
                  gsf,&
                  nder,&
                  maxd,&
                  maxcoefs,&
                  bdrf_ms, &
                  u0,&
                  theta,&
                  phi,&
                  rint, &!output
                  dtaua, &!output
                  dtaus,&!output
                  dphase,&!output
                  dalb, &!output
                  ierr)
             if (ierr .ne. 0) return

          elseif ( flag%rtm==2 )then
             atm%taua = taua
             atm%taus = taus
             set%deltam = .true.
             call lintran_calculate(execution, atm, set, lintran_instance, rintms, drv,ierr)
             rint(1:nstokes) = rintms(1:nstokes,1)
             call lintranv2_return( nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua, dtaus, dphase, dalb, drv )
          elseif ( flag%rtm==3 )then
             atm%taua = taua
             atm%taus = taus 
             call lintran_nakajima_multiscattering(atm, set, lintran_instance, rintms, drv, lintran_stat)
             rint(1:nstokes) = rintms(1:nstokes,1)
             call lintranv2_return( nstokes, nrt, maxd, dz_ss,dalb_ss, dtaua, dtaus, dphase, dalb, drv )
          endif


          !*** Take the logarithm of ms reflectance/normalize Stokes components for interpolation
          !*** Stokes vector
          !            if (rint(1)<0.d0 ) pause
          rint_tau(itau,itau_2nd,1) = log(rint(1))
          do i_st = 2, nstokes
             rint_tau(itau, itau_2nd, i_st) = rint(i_st)/rint(1)
          enddo

          !*** Adjust derivatives for logarithm of ms reflectance/normalization of Stokes components

          !*** Absorption optical thickness
          do j = 1, nrt
             dtaua_tau(itau, itau_2nd, j, 1) = (dtaua(1,j))/(rint(1))
             do i_st = 2, nstokes
                dtaua_tau(itau, itau_2nd, j,  i_st) = &
                     dtaua(i_st,j)/rint(1) - &
                     rint(i_st)/rint(1)**2 * dtaua(1, j)
             enddo
          enddo

          !*** Aerosols
          !*** Combine derivatives of ms reflectance wrt. tau with derivatives of tau wrt. parameters
          do i = 1, ntype_aer
             if( ((aerosol(i)%CirrusFlag==1 .and. MinCOTFlag==0) .or. &
                  (aerosol(i)%CirrusFlag .ne. 1 .and. MinAOTFlag==0)) .and. maxd > 0 )then
                i1 = nder(1)
                i2 = nder(maxd)
                call combine_der(&
                     dtaua(:,i1:i2),&
                     dtaus(:,i1:i2),&
                     dphase,&
                     aerosol_optic(i)%dtaua_all,&
                     aerosol_optic(i)%dtaus_all, &
                     aerosol_optic(i)%dphase_all,&
                     maxcoefs,&
                     k_tmp_phase,&
                     k_tmp_taua,&
                     k_tmp_taus,&
                     nder,&
                     maxd,&
                     aerosol(i)%aerosolflags,&
                     npar)
             else
                forall (i_st=1:nstokes, j=1:npar)
                   k_tmp_phase(i_st, j) = 0.d0
                   k_tmp_taua(i_st, j) = 0.d0
                   k_tmp_taus(i_st, j) = 0.d0
                end forall
             endif
             !*** Convert from dI/dx to dlogI/dx
             do j = 1, npar
                kmat_aerosol_phase(i, itau, itau_2nd, j, 1) = (k_tmp_phase(1, j))/(rint(1))
                kmat_aerosol_taua(i, itau, itau_2nd, j, 1)  = (k_tmp_taua(1, j))/(rint(1))
                kmat_aerosol_taus(i, itau, itau_2nd, j, 1)  = (k_tmp_taus(1, j))/(rint(1))
                do i_st = 2, nstokes
                   kmat_aerosol_phase(i, itau, itau_2nd, j, i_st) = &
                        k_tmp_phase(i_st,j)/rint(1) - &
                        rint(i_st)/rint(1)**2 * k_tmp_phase(1,j)
                   kmat_aerosol_taua(i, itau, itau_2nd, j, i_st) = &
                        k_tmp_taua(i_st,j)/rint(1) - &
                        rint(i_st)/rint(1)**2 * k_tmp_taua(1,j)
                   kmat_aerosol_taus(i, itau, itau_2nd, j, i_st) = &
                        k_tmp_taus(i_st,j)/rint(1) - &
                        rint(i_st)/rint(1)**2 * k_tmp_taus(1,j)
                enddo
             enddo ! loop over npar
          enddo ! loop over ntype_aer

          !** Albedo 
          kmat_alb_tau(itau,itau_2nd,1) = (dalb(1))/(rint(1))
          do i_st = 2,nstokes
             kmat_alb_tau(itau,itau_2nd,i_st) = &
                  dalb(i_st)/rint(1) - &
                  rint(i_st)/rint(1)**2 * dalb(1)
          enddo

          !** Scattering optical thickness
          do k = 1, nrt
             kmat_taus_tau(itau, itau_2nd, k, 1) = (dtaus(1,k))/(rint(1))
             do i_st = 2, nstokes
                kmat_taus_tau(itau, itau_2nd, k, i_st) = &
                     dtaus(i_st,k)/rint(1) - &
                     rint(i_st)/rint(1)**2 * dtaus(1, k)
             enddo
          enddo

       enddo ! loop over ntau_2nd
    enddo ! loop over ntau
    !**** End loop over all selected tau values  

    !**** Interpolate ms reflectance and derivatives from tau-subgrid to wavelength grid    

    !** If only one absorber
    if(nmol.eq.1 .or. ntau_2nd_max .eq. 1) then
       call interpolate_lbl_1species(&
            nwave, &
            nrt, &
            tau_grid, &
            tau_grid_arr, &
            aerosol, &
            ntau,&
            taua_tot,&
            taua_mol,&
            rint_tau,&
            rint_fine_temp,&
            dtaua_tau,&
            kmat_taua_fine,&
            kmat_aerosol_phase,& 
            kmat_aerosol_taua,&  
            kmat_aerosol_taus,&
            kmat_alb_tau,&
            kmat_taus_tau,&
            kmat_taus_fine, &
            tmp_daerosol_phase, &
            tmp_daerosol_taua, &
            tmp_daerosol_taus, &
            deriv_rt%alb(:,1,:))   
       !** If two different absorbers
    elseif(nmol.ge.2) then
       call interpolate_lbl_2species(&
            nwave, &
            nrt, &
            ntau_2nd, &
            tau_grid, &
            tau_grid_2nd, &
            tau_grid_arr, &
            aerosol, &
            ntau,&
            imol_1st,&
            imol_2nd,&
            taua_tot_species,&
            taua_tot,&
            taua_mol,&
            rint_tau,&
            rint_fine_temp,&
            dtaua_tau,&
            kmat_taua_fine,&
            kmat_aerosol_phase,& 
            kmat_aerosol_taua,& 
            kmat_aerosol_taus,&
            kmat_alb_tau,&
            kmat_taus_tau,&
            kmat_taus_fine, &
            tmp_daerosol_phase, &
            tmp_daerosol_taua, &
            tmp_daerosol_taus, &
            deriv_rt%alb(:,1,:))
    endif

    !**** Linearly correct the ms reflectance for the spectral dependence
    !**** of scattering optical thickness and albedo
    !**** So far, scattering has only been considered at taus('imid'), albedo as average albedo  
    forall(l = 1:nwave, i = 1:nstokes)
       rint_fine_ms(l,i) = rint_fine_temp(l,i) + &
            sum(kmat_taus_fine(l, :, i)*(taus_mol(l,:)+taus_aer(l,:)-taus(:))) + &
            sum(kmat_taua_fine(l, :, i)*(taua_aer(l,:)-taua_aer(imid,:))) + &
            deriv_rt%alb(l,1,i)*(albedo_array(l) - albedo_av)
       !            deriv_rt%alb(l, i)*(albedo_array(l) - albedo_av)
    end forall

    if (deriv_flag) then
       !*** Linearly correct for spectral dependence of derivatives towards scattering and absorption cross sections
       do i = 1,nstokes
          do j = 1, ntype_aer
             do k = 1, npar_mie ! Only Mie-parameters are linearly interpolated to correct for spectral dependence of dsigma/dMIEparam 
                if(aerosol(j)%AerosolFlags(k).ne.0)then
                   do l = 1,nwave

                      tmp_daerosol_taua(l,k,j,i) = tmp_daerosol_taua(l,k,j,i) * &
                           ( 1D0 + aerosol_optic(j)%dcabs_dwave(k)*( wavelength(l) - wavelength(imid) ) )
                      tmp_daerosol_taus(l,k,j,i) = tmp_daerosol_taus(l,k,j,i) * &
                           ( 1D0 + aerosol_optic(j)%dcsca_dwave(k)*( wavelength(l) - wavelength(imid) ) )

                   enddo
                endif
             enddo
          enddo
       enddo
    endif

    !**** Undo logarithm of ms reflectance/normalization of Stokes coefficients       
    do l = 1, nwave
       rint_fine_ms(l, 1) = exp(rint_fine_ms(l, 1))
       rint_fine_temp(l, 1) = exp(rint_fine_temp(l, 1))
       do i_st = 2, nstokes
          rint_fine_ms(l, i_st) = rint_fine_temp(l, i_st)*rint_fine_ms(l, 1)
          !          rint_fine_ms(l, i_st) = rint_fine_ms(l, i_st)*rint_fine_ms(l, 1)
          !rint_fine_ms(l, i_st) = rint_fine_ms(l, i_st)*rint_fine_temp(l, 1)
       enddo
    enddo

    !*** Get derivatives: undo logarithm of ms reflectance/normalization of Stokes components
    !*** First Stokes parameter
    if (deriv_flag) then
       do k = 1, nrt
          do l = 1, nwave
             kmat_taua_fine(l, k, 1) = kmat_taua_fine(l, k, 1)*rint_fine_temp(l,1)
             kmat_taus_fine(l, k, 1) = kmat_taus_fine(l, k, 1)*rint_fine_temp(l,1)
             !             kmat_taua_fine(l, k, 1) = kmat_taua_fine(l, k, 1)*rint_fine_ms(l,1)
             !             kmat_taus_fine(l, k, 1) = kmat_taus_fine(l, k, 1)*rint_fine_ms(l,1)
          enddo
       enddo
       do l = 1, nwave
          deriv_rt%alb(l, 1, 1) = deriv_rt%alb(l, 1, 1)*rint_fine_temp(l,1)
          !          deriv_rt%alb(l, 1, 1) = deriv_rt%alb(l, 1, 1)*rint_fine_ms(l,1)
       enddo
       do i = 1, ntype_aer
          do k = 1, npar
             do l = 1, nwave
                tmp_daerosol_phase(l, k, i, 1) = tmp_daerosol_phase(l, k, i, 1)*rint_fine_temp(l,1)
                tmp_daerosol_taua(l, k, i, 1)  = tmp_daerosol_taua(l, k, i, 1) *rint_fine_temp(l,1)
                tmp_daerosol_taus(l, k, i, 1)  = tmp_daerosol_taus(l, k, i, 1) *rint_fine_temp(l,1)
                !                tmp_daerosol_phase(l, k, i, 1) = tmp_daerosol_phase(l, k, i, 1)*rint_fine_ms(l,1)
                !                tmp_daerosol_taua(l, k, i, 1)  = tmp_daerosol_taua(l, k, i, 1) *rint_fine_ms(l,1)
                !                tmp_daerosol_taus(l, k, i, 1)  = tmp_daerosol_taus(l, k, i, 1) *rint_fine_ms(l,1)

             enddo
          enddo
       enddo

       !*** Other Stokes parameters
       do i = 2, nstokes
          do l = 1, nwave
             !SIGN(MAX(ABS(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) !!!avoid zero dividing 
             kmat_taua_fine(l, :, i) = &
                  kmat_taua_fine(l, :, i)*rint_fine_temp(l,1)  * &
                  rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                  kmat_taua_fine(l, :, 1)*rint_fine_temp(l,i)
             kmat_taus_fine(l, :, i) = &! 
                  kmat_taus_fine(l, :, i)*rint_fine_temp(l,1)  * &
                  rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                  kmat_taus_fine(l, :, 1)*rint_fine_temp(l,i)

             deriv_rt%alb(l, 1, i) = &
                  deriv_rt%alb(l, 1, i)*rint_fine_ms(l, 1)  * &
                  rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                  deriv_rt%alb(l, 1, 1)*rint_fine_temp(l, i)

             do k = 1, ntype_aer
                tmp_daerosol_phase(l, :, k, i) = &					
                     tmp_daerosol_phase(l, :, k, i)*rint_fine_temp(l,1)  * &
                     rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                     tmp_daerosol_phase(l, :, k, 1)*rint_fine_temp(l,i)

                tmp_daerosol_taua(l, :, k, i) = &				
                     tmp_daerosol_taua(l, :, k, i)*rint_fine_temp(l,1)  * &
                     rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                     tmp_daerosol_taua(l, :, k, 1)*rint_fine_temp(l,i)

                tmp_daerosol_taus(l, :, k, i) = &				
                     tmp_daerosol_taus(l, :, k, i)*rint_fine_temp(l,1)  * &
                     rint_fine_temp(l,i)/sign(max(abs(rint_fine_temp(l,i)), SMALL), rint_fine_temp(l,i)) + &
                     tmp_daerosol_taus(l, :, k, 1)*rint_fine_temp(l,i)
             enddo

          enddo
       enddo

       !*** Combine derivatives
       deriv_rt%aerosol = tmp_daerosol_phase + tmp_daerosol_taua +  tmp_daerosol_taus

    endif

    deallocate(&
         rint_fine_temp, &
         tmp_daerosol_phase, & 
         tmp_daerosol_taua, &
         tmp_daerosol_taus, &
         stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif


    !*** Absorbers
    if (deriv_flag) then
       do imol = 1, nmol   
          !*** For oxygen the derivative of reflectance wrt. to the oxygen subcolumn
          !*** is calculated assuming constant mixing ratio (but variable air column)
          if(win_ini%type_x(imol)==7) then
             forall(l=1:nwave, i=1:nstokes, k=1:nrt)
                deriv_rt%densmol(l, k, imol, i) = kmat_taua_fine(l, k, i)*cross_mol(l,k,imol) + &
                     kmat_taus_fine(l, k, i)*csray(l)/relo2
             end forall
             !*** For absorbers other than oxygen the derivative of reflectance wrt. to the absorber subcolumn
             !*** is calculated assuming constant air column (but variable mixing ratio)
          else
             forall(l=1:nwave, i=1:nstokes, k=1:nrt)
                deriv_rt%densmol(l, k, imol, i) = kmat_taua_fine(l, k, i)*cross_mol(l, k, imol)
             end forall
          endif
       enddo

       !** Temperature offset
       if(flag%temp==1) then
          do i = 1, nstokes
             do l = 1, nwave
                deriv_rt%T(l, i) = sum(kmat_taua_fine(l, :, i)*dtaua_T(l,:))
             enddo
          enddo
       endif

       !** O2 column/pressure
       if(flag%O2==1) then
          do i = 1, nstokes
             do k = 1, nrt
                do l = 1, nwave
                   deriv_rt%P(l, k, i) = kmat_taua_fine(l, k, i)*dtaua_P(l, k)
                enddo
             enddo
          enddo
       endif
    endif

    !*******************************************************************************************
    !**** SINGLE SCATTERING
    !**** Loop over all wavelength (line-by-line) 
    if ( flag%rtm==2 .or. flag%rtm==3 )then
       execution = (/.true.,.false.,.false./)  !(/ssc, dsc, msc/) -> single scattering solution
    endif

    if(allocated(z_ss_int)) deallocate(z_ss_int)
    allocate(z_ss_int(4, nrt, nwave), stat=ierr)
    z_ss_int = 0.d0
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    call optic_all_ssc( &
         aero_lut, cirrus_lut, &
         flag%rtm, &! For polarization, rotation to the local meridian plane is needed only for rtm%flag==1
         iband, wavelength, scat_angle, &
         u0, uv, phi, aerosol, nder, nrt, &
         aerosol_optic, taus_aer, taus_mol, z_ss_int, &
         zss_lin, &! z matrix for Lintran
         ierr)

    do iwave = 1, nwave
       if (flag%rtm == 3) then
          if ( glintflag==1) then
             atm%bdrf_ssg(1,1,1) =  albedo_array(iwave) + bdrf_ss_cst(1)
          else
             atm%bdrf_ssg(1,1,1) = albedo_array(iwave)       
          endif
          do l = 1, nrt
             atm%taua(l) = taua_mol(iwave, l) + taua_aer(iwave, l)
             atm%taus(l) = taus_mol(iwave, l) + taus_aer(iwave, l)
             do n = 1, nper
                atm%phase_ssg(n, l, 1) = zss_lin(index_ist(n), index_jst(n), l, iwave)
             enddo
          enddo
       else
          if ( glintflag==1) then
             bdrf_ss(1) =  albedo_array(iwave) + bdrf_ss_cst(1)
          else
             bdrf_ss(1) = albedo_array(iwave)       
          endif
          do l = 1, nrt
             taua(l) = taua_mol(iwave, l) + taua_aer(iwave, l)
             taus(l) = taus_mol(iwave, l) + taus_aer(iwave, l)
             z_ss(1:nstokes, l) = z_ss_int(1:nstokes, l, iwave)
          enddo
       endif


       !*** Call single-scattering RTM
       !*** output @ iwave: rint_ss=ss reflectance,
       !*** dtaua_ss=deriv. wrt. absorption optical depth,
       !*** dtaus_ss=deriv. wrt. scattering optical depth, 
       !*** dalb_ss=deriv. wrt. albedo 
       if ( flag%rtm==1 )then
          call SINGLE_SCAT(&
               nrt, &
               taua, &
               taus, &
               bdrf_ss, &
               u0, &
               uv, &
               z_ss, &
               rint_ss, &!output
               dtaus_ss, &!output
               dtaua_ss, &!output
               dz_ss, &!output
               dalb_ss, &!output
               nder, &
               maxd)

       elseif (flag%rtm==2) then  
          call lintranv2_assign( nstokes, nrt, maxd, nder, cloudflag, maxcoefs, plmom_in, bdrf_ms, bdrf_ss, &
               zss_lin(:,:,:,iwave), surf_emi, taua, taus, atm, drv, set )!!! for single scattering zss_lin matters
          set%deltam = .true. 
          if (deriv_flag) then
             call lintran_calculate(execution, atm,set, lintran_instance, rintss, drv, ierr)
             call lintranv2_return( nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua_ss, dtaus_ss, dphase, dalb, drv )
          else
             call lintran_calculate(execution, atm,set, lintran_instance, rintss, ierr)
          endif
          rint_ss(1:nstokes) = rintss(1:nstokes,1)
       elseif( flag%rtm==3 ) then
          if (deriv_flag) then
             if (nstokes==1) then
                call remotec_single_scat_with_derivatives_nst1(atm, set, lintran_instance, rintss, drv) 
             elseif (nstokes==3) then
                call remotec_single_scat_with_derivatives_nst3(atm, set, lintran_instance, rintss, drv)
             endif
          else
             if (nstokes==1) then
                call remotec_single_scat_without_derivatives_nst1(atm, lintran_instance, rintss)
             elseif(nstokes==3) then 
                call remotec_single_scat_without_derivatives_nst3(atm, lintran_instance, rintss)
             endif
          endif
          rint_ss(1:nstokes) = rintss(:,1)
          call lintranv2_return_ss(nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua_ss, dtaus_ss, drv )           
       endif

       if (deriv_flag) then

          !**** Combine ss derivatives of reflectance wrt. tau with derivatives of tau wrt. parameters
          !*** Aerosol parameters
          do i = 1, ntype_aer	    
             if( ((aerosol(i)%CirrusFlag==1 .and. MinCOTFlag==0) .or. (aerosol(i)%CirrusFlag .ne. 1 .and. MinAOTFlag==0)) &
                  .and. maxd>0 ) then
                i1 = nder(1)
                i2 = nder(maxd) 
                if (iwave<imid) then

                   aerosol_optic(i)%dzss_all = aerosol_optic(i)%dzss_all_1 + ( aerosol_optic(i)%dzss_all_mid - aerosol_optic(i)%dzss_all_1 )/ (wavelength(imid)-wavelength(1))*(wavelength(iwave)-wavelength(1))
                else
                   aerosol_optic(i)%dzss_all = aerosol_optic(i)%dzss_all_mid + ( aerosol_optic(i)%dzss_all_nwave - aerosol_optic(i)%dzss_all_mid )/ (wavelength(nwave)-wavelength(imid))*(wavelength(iwave)-wavelength(imid))

                endif
                call combine_der_ss(&
                     dtaua_ss(:,i1:i2),&
                     dtaus_ss(:,i1:i2),&
                     dz_ss,&
                     aerosol_optic(i)%dtaua_all,&
                     aerosol_optic(i)%dtaus_all, &
                     aerosol_optic(i)%dzss_all,&
                     k_tmp_phase,&
                     k_tmp_taua,&
                     k_tmp_taus,&
                     nder,&
                     maxd,&
                     aerosol(i)%aerosolflags,&
                     npar)
             else
                forall (i_st=1:nstokes, k = 1:npar)
                   k_tmp_phase(i_st, k) = 0.0d0
                   k_tmp_taua(i_st, k)  = 0.0d0
                   k_tmp_taus(i_st, k)  = 0.0d0
                end forall
             endif
             do k = 1, npar
                do i_st = 1, nstokes
                   if  (k .le. npar_mie) then
                      !*** Linearly correct for spectral dependance off dcsca/dx_aer and dcabs/dx_aer:
                      kmat_aerosol_ss(iwave, k, i, i_st) =  k_tmp_phase(i_st, k) + &
                           k_tmp_taua(i_st, k)*( 1D0 + aerosol_optic(i)%dcabs_dwave(k)* &
                           ( wavelength(iwave) - wavelength(imid) ) ) + & 
                           k_tmp_taus(i_st, k)*( 1D0 + aerosol_optic(i)%dcsca_dwave(k)* &
                           ( wavelength(iwave) - wavelength(imid) ) )
                   else
                      kmat_aerosol_ss(iwave, k, i, i_st) =  k_tmp_phase(i_st, k) + &  !<=========================== D.S
                           k_tmp_taua(i_st, k) / sum(aerosol_optic(i)%taua_aer(imid,:))* sum(aerosol_optic(i)%taua_aer(iwave,:)) + &   !spectral dependence of derivatives to height parameters
                           k_tmp_taus(i_st, k) / sum(aerosol_optic(i)%taus_aer(imid,:))* sum(aerosol_optic(i)%taus_aer(iwave,:))       

                   endif
                enddo
             enddo
          enddo

          !*** Albedo      
          kmat_alb_ss(iwave,1:nstokes) = dalb_ss(1:nstokes)

          !*** Absorbers
          do imol = 1, nmol            
             !*** For oxygen the derivative of reflectance wrt. to the oxygen subcolumn
             !*** is calculated assuming constant mixing ratio (but variable air column)
             if(win_ini%type_x(imol)==7) then
                do k = 1, nrt
                   kmat_densmol_ss(iwave, k, imol, 1:nstokes) = (&
                        dtaua_ss(1:nstokes, k)*cross_mol(iwave,k,imol) + &
                        dtaus_ss(1:nstokes, k)*csray(iwave)/relo2) 
                enddo
                !*** For absorbers other than oxygen the derivative of reflectance wrt. to the absorber subcolumn
                !*** is calculated assuming constant air column (but variable mixing ratio)
             else  
                do k = 1, nrt
                   kmat_densmol_ss(iwave, k, imol, 1:nstokes) = &
                        dtaua_ss(1:nstokes,k)*cross_mol(iwave, k, imol)          
                enddo
             endif
          enddo
          !*** Temperature offset
          if(flag%temp==1) then
             do i = 1, nstokes
                deriv_rt%T(iwave,i) = deriv_rt%T(iwave, i) + sum(dtaua_ss(i,:)*dtaua_T(iwave,:))
             enddo
          endif
          !*** O2 column/pressure
          if(flag%O2==1) then
             do i = 1, nstokes
                do k = 1, nrt
                   deriv_rt%P(iwave,k,i) = deriv_rt%P(iwave,k,i) + dtaua_ss(i,k)*dtaua_P(iwave,k)
                enddo
             enddo
          endif
          rint_fine_ss(iwave,1:nstokes) = rint_ss(1:nstokes)
       else
          rint_fine(iwave,1:nstokes) = rint_fine_ms(iwave,1:nstokes) + rint_ss(1:nstokes)
       endif

    enddo
    !**** End loop over all wavelengths

    !**** Clean up all lintran instances and internals
    if ( flag%rtm==2 .or. flag%rtm==3 ) then
       call lintran_deallocate(atm, set, drv, ierr)
       call lintran_close(lintran_instance,ierr)
    endif

    ! Quickly return to the derivatives if-clause

    !*** The final hi-res spectrum, including single scattering, multiple scattering plus an optional constant, is rint_fine:
    if (deriv_flag) then
       do i_st = 1, nstokes
          rint_fine(:, i_st) = rint_fine_ms(:, i_st) + rint_fine_ss(:, i_st)
          deriv_rt%alb(:, 1, i_st) = deriv_rt%alb(:, 1, i_st) + kmat_alb_ss(:, i_st)
          do imol = 1,nmol           
             do j = 1, nrt
                deriv_rt%densmol(:, j, imol, i_st) = deriv_rt%densmol(:, j, imol, i_st) +&
                     kmat_densmol_ss(:, j, imol, i_st)
             enddo
          enddo
          do n = 1, ntype_aer
             do k = 1, npar
                deriv_rt%aerosol(:, k, n, i_st) = deriv_rt%aerosol(:, k, n, i_st) + kmat_aerosol_ss(:, k, n, i_st)
             enddo
          enddo
       enddo

    endif

    deallocate(aerosol_optic, stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_EFFICIENT: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine rad_trans_efficient ! }}}

  !------------------------------------------------------------------------------  
  !> Time-consuming radiative transfer model. 
  !! Line-by-line calculation for multiple and single scattering
  !------------------------------------------------------------------------------  
  subroutine rad_trans_lbl( & ! {{{
       aero_lut, cirrus_lut, &            
       flag, glintflag, &
       iband, wavelength, albedo_array, sza, theta, phi, wspeed, &
       aerosol, minaotflag, mincotflag, nrt, nder, atm_xs, &
       dvair, dvair_old, vmr_h2o, play_old, &
       win_ini, x_molec, &
       rint_fine, deriv_rt, taua_tot, ot, cot, ExitXSFlag, MaxOTFlag, ierr, deriv_flag)
    !*** Input
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    type(settings_flags), intent(in) :: flag
    integer, intent(in) :: glintflag, iband
    real(double), dimension(:), intent(in) :: wavelength 
    real(double), dimension(:), intent(in) :: albedo_array 
    real(double), intent(in) :: sza, theta, phi, wspeed
    type(aero), dimension(:), intent(in) :: aerosol 
    integer, intent(in) :: minaotflag, mincotflag, nrt
    integer, dimension(:), intent(in) :: nder     
    type(atmosphere), intent(in) :: atm_xs
    real(double), dimension(:), intent(in) :: dvair
    real(double), dimension(:), intent(in) :: dvair_old
    real(double), dimension(:), intent(in) :: vmr_h2o   
    real(double), dimension(:), intent(in) :: play_old 
    type(window_ini), intent(in) :: win_ini
    real(double), dimension(:,:), intent(in) :: x_molec  
    logical, intent(in) :: deriv_flag
    !*** Input/output
    real(double), dimension(:,:), intent(inout) :: rint_fine
    type(derivatives), intent(inout) :: deriv_rt
    !*** Output
    integer, intent(out) :: ExitXSFlag, MaxOTFlag, ierr
    real(double), intent(out) :: ot, cot
    real(double), dimension(:), allocatable, intent(out) :: taua_tot
    !*** local variables        
    integer, dimension(nrt) :: ncoefs   
    real(double), dimension(nstokes) :: dalb, rint, dalb_ss, rint_ss, dalb_tmp
    real(double), dimension(nstokes, nrt) :: dtaua, dtaus, dtaus_ss, dtaua_ss
    real(double), dimension(nstokes, npar) :: k_tmp, k_tmp_ss
    real(double), dimension(nstokes, npar) :: k_tmp_phase, k_tmp_taua, k_tmp_taus
    real(double), dimension(:, :, :, :), allocatable:: dphase, dphase_tmp
    real(double), dimension(:,:), allocatable :: rint_fine_ss
    real(double), dimension(:,:), allocatable :: rint_fine_ms
    real(double), dimension(:,:,:), allocatable :: rint_tau
    real(double), dimension(:,:,:,:), allocatable :: kmat_densmol_ss
    real(double), dimension(:,:), allocatable :: kmat_alb_ss
    real(double), dimension(:,:,:), allocatable :: kmat_alb_tau
    real(double), dimension(:,:,:,:), allocatable :: dtaua_tau
    real(double), dimension(:,:,:,:), allocatable :: kmat_taus_tau
    real(double), dimension(:,:,:), allocatable :: kmat_taus_fine
    real(double), dimension(:,:,:), allocatable :: kmat_taua_fine
    real(double), dimension(:,:,:,:), allocatable :: kmat_aerosol_ss
    real(double), dimension(:, :), allocatable:: dz_ss
    real(double) :: bdrf_ms(nstokes, nstokes, mxhalf+2, mxhalf+2, 0:maxstr-1)
    real(double) :: bdrf_ss(nstokes), bdrf_ss_cst(nstokes)     
    integer :: i, j, k, l, imid, i1, i2, iwave, i_st, &
         imol, nmol, &
         ntau, ntau_2nd_max, ntype_aer, nwave
    real(double) :: wave_mid, albedo, u0, uv, scat_angle 
    type(aero_opt), dimension(:), allocatable :: aerosol_optic
    real(double), dimension(:,:), allocatable :: taua_aer
    real(double), dimension(:,:), allocatable :: taus_aer
    real(double), dimension(nrt) :: taua 
    real(double), dimension(nrt) :: taus
    ! Phase matrix coefficients for aerosol scattering
    real(double), dimension(nstokes,nstokes,0:maxleg,nrt) :: plmom_aer
    integer, dimension(nrt) :: maxcoefs       ! Number of coefficients for phase matrix (nrt)
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt) :: plmom_in  !Phase matrix coefficients
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt) :: plmom     ! Phase matrix coefficients after Nakajima correction
    real(double), dimension(nstokes, nrt) :: z_ss          ! Z matrix for single scattering (nstokes,nrt)
    real(double), dimension(4, 4, nrt, size(wavelength)) :: zss_lin! Z matrix for Lintran
    real(double), dimension(:,:,:), allocatable :: z_ss_int    ! Z matrix for single scattering, interpolated for wavelength (nstokes,nrt,nwave_hi)
    !*** Optical properties of molecules
    real(double), dimension(:,:), allocatable :: taua_mol, taus_mol
    real(double), dimension(:,:,:), allocatable ::  taua_mol_species ! Molecular absorption optical thickness per type (nwave_hi,nrt,ntype)
    real(double), dimension(:,:,:), allocatable :: cross_mol         ! Molecular absorption cross section (nwave_hi,nrt,ntype) 
    real(double), dimension(:), allocatable :: csray                 ! Rayleigh scattering cross section (nwave_hi)
    real(double), dimension(:,:), allocatable :: dtaua_T             ! Derivative of absorption optical thickness wrt temperature (nwave_hi,nrt) 
    real(double), dimension(:,:), allocatable :: dtaua_P             ! Derivative of absorption optical thickness wrt pressure (nwave_hi,nrt,ntype) 
    real(double), dimension(nstokes,nstokes,0:2,nrt) :: plmom_ray    ! Phase matrix coefficients for Rayleigh scattering
    integer :: maxd
    type(phase_mat_type),dimension(0:maxstr) :: pm_in
    type(gauss_quad_type) :: gauss_quad
    type(gsf_type) :: gsf
    character(stringlen) :: message
    !*** The Lintran input instances
    logical :: cloudflag ! Flag for using cloudy settings
    type(lintran_atmosphere) :: atm
    type(lintran_derivatives) :: drv
    type(lintran_settings) :: set
    type(lintran_class) :: lintran_instance
    logical, dimension(3) :: execution ! Flags for single, double and multi-scattering, in that order.
    real(double), dimension(nstokes,1) :: rintms, rintss! Main output LINTRAN
    real(double), dimension(nstokes,1) :: rint_lintran  ! Main output LINTRAN
    real(double) :: surf_emi
    real(double), parameter :: degrees = acos(0.0) / 90.0
    !--------------------------------------------------------------------------------

    !*** Initialize error identifier 
    ierr = 0
    MaxOTFlag = 0

    !**** Initialization and Allocation
    surf_emi = 0D0
    uv = DCOS(theta/180.*pi) 
    call u0_kasten_and_young(sza, u0)
    scat_angle = DACOS(&
         -u0*uv + &
         DSQRT(1.D0-uv**2)* &
         DSQRT(1.D0-u0**2)* &
         DCOS(phi/180.D0*PI))      
    ntau = win_ini%ntau
    ntau_2nd_max = win_ini%ntau_2nd_max
    nwave = size(wavelength)      !win_ini%nwave_hi
    nmol = win_ini%ntype   ! number of absorbers      
    imid = (nwave+1)/2
    wave_mid = wavelength(imid)
    ntype_aer = size(aerosol)
    maxd = size(nder)
    allocate(&
         dz_ss(nstokes, maxd), &
         dphase(nstokes, 0:maxleg, nper, maxd), &
         dphase_tmp(nstokes, 0:maxleg, nper, maxd), &
         kmat_densmol_ss(nwave, nrt, nmol, nstokes),&
         kmat_alb_ss(nwave, nstokes),&
         rint_tau(ntau, ntau_2nd_max, nstokes),&
         rint_fine_ms(nwave, nstokes),&
         rint_fine_ss(nwave, nstokes),&
         dtaua_tau(ntau, ntau_2nd_max, nrt, nstokes),&
         kmat_alb_tau(ntau, ntau_2nd_max, nstokes),&
         kmat_taus_tau(ntau, ntau_2nd_max, nrt, nstokes),&
         kmat_taus_fine(nwave, nrt, nstokes),&
         kmat_taua_fine(nwave, nrt, nstokes),&
         kmat_aerosol_ss(nwave, npar, ntype_aer, nstokes), &
         taua_tot(nwave), & 
         stat=ierr) 
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_LBL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    !*** Calculate molecular optical quantities
    call calculate_optic_mol_prop( &
         iband, flag%xs, flag%O2, flag%temp, &  !Input
         wavelength, &
         nrt, &
         atm_xs, &
         win_ini%xsdb, &
         win_ini%type_xsdb, &
         x_molec, &
         dvair, &
         dvair_old, &
         vmr_h2o, &
         play_old, &
         taua_mol, &                            !Output
         taus_mol, &
         taua_mol_species, &
         cross_mol, &
         csray, &
         dtaua_T, &
         dtaua_P, &
         plmom_ray, &
         ExitXSFlag, ierr)
    if (ExitXSFlag .ne. 0 .or. ierr .ne. 0) return


    ! Get taua_tot for fluorescence calculations:
    do l = 1, nwave 
       !  cut k-binning grid at max value tatot
       taua_tot(l) = min(sum(taua_mol(l,:)),tatot)
    enddo

    !*** Calculate aerosol optical quantities, only for middle of absorption band
    allocate(aerosol_optic(ntype_aer), & 
         taus_aer(nwave, nrt),&
         taua_aer(nwave, nrt), &
         stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_LBL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    !*** Calculate aerosol scattering and absorption optical thicknes and phase function for multiple scattering 
    call optic_all_msc( &
         aero_lut, cirrus_lut, &
         iband, wavelength, &
         aerosol, nder, nrt, aerosol_optic, &
         taua_aer, taus_aer, ot, cot, &
         plmom_aer, maxcoefs, ierr)
    if (ierr .ne. 0) return

    ot = ot + cot ! total optical thickness

    if (flag%output >= 2) then   
       write(message,'(a,2ES12.3,a,ES12.3)') 'total OT and COT: ', ot, cot , ' at wavelength ', wave_mid
       call writelog(message, 3)
    endif

    !*** Error checks
    !*** Check for NAN values:
    if(ot .ne. ot)then
       MaxOTFlag = 1
       return
    endif
    !*** Check for too high OTs
    if(ot > 1.d0 .and. .false.) then ! Hack: Turn off total optical depth check
       MaxOTFlag = 1
       return
    endif

    ! Flag for using cloudy settings for Lintran
    cloudflag = ot .gt. 3.0

    !*** Interface aerosol and molecular quantities to RTM, only for middle of absorption band
    !*** taus,taua: total scattering/absorption optical thickness (nmb. of layers)
    !*** dtaus_all,dtaua_all: derivative of total scattering/absorption optical thickness wrt aerosol parameters
    !*** plmom_in, plmom_ss: phase matrix coefficients
    !*** dphase_all: derivatives of coefficients of aerosol scattering matrix wrt aerosol parameters

    call optic_all_der( &
         imid, &
         aerosol_optic, aerosol, nder, nrt, maxcoefs, &
         taua_aer, taus_aer, plmom_aer, &
         taua_mol, taus_mol, plmom_ray, &
         taua, taus, plmom_in, &
         ierr)
    if (ierr .ne. 0) return

    !*** Initilization of RTM
    if ( flag%rtm==1 )then
       call perturbation_init(&
            u0,&
            theta,&
            gauss_quad,&
            gsf, &
            ierr)
       if (ierr .ne. 0) return         

       !**** Nakajima correction for strongly forward peaked phase functions
       call nakajima_correction(plmom_in, plmom, nder, maxd, nrt)
       !***Nakajima switch: off = (plmom=plmom_in), on = (plmom!=plmom_in)
       !   plmom = plmom_in

       !**** Calculate Fourier coefficients of phase matrix
       forall(k=1:nrt) ncoefs(k) = min(maxcoefs(k), maxstr-1)

       ! Phase coefficients are fixed to wavelength imid. That is an approximation, not
       ! exactly line by line. It is okay as long as the phase function does not change
       ! much in the window.
       call pm_coeff(&
            nrt, &
            nstokes,&
            ncoefs,&
            plmom,&
            gsf,&
            pm_in, &
            ierr)
       if (ierr .ne. 0) return
    elseif ( flag%rtm==2 .or. flag%rtm==3 ) then
       call lintran_init(nstokes, MAXSTR, nrt, 1, .true., lintran_instance, ierr)
       ! Convert the azimuthal difference to the convention of de Haan et al. (1987).
       ! That is what Lintran 2 desires.
       call lintran_provide(real(u0),(/real(uv)/),(/-real(phi*degrees)/),.false.,lintran_instance, ierr)       
    endif

    !*** Calculate single scattering phase matrix
    if(allocated(z_ss_int)) deallocate(z_ss_int)
    allocate(z_ss_int(nstokes,nrt,nwave), stat=ierr)
    if (ierr .ne. 0) then
       write(message,*) 'RAD_TRANS_LBL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif
    call optic_all_ssc( &
         aero_lut, cirrus_lut, flag%rtm, &
         iband, wavelength, scat_angle, &
         u0, uv, phi, aerosol, nder, nrt, &
         aerosol_optic, taus_aer, taus_mol, z_ss_int, &
         zss_lin, &
         ierr)
    if (ierr .ne. 0) return

    bdrf_ss = 0.d0
    bdrf_ms = 0.D0
    !*** OCEAN glint
    if (glintflag==1) then
       !*** Bidirectional reflection distr. functions calculated from the ocean model   
       if ( flag%rtm==1 )then
          call ocean_fou(u0, uv, wspeed, gauss_quad, bdrf_ms)
       elseif ( flag%rtm==2 .or. flag%rtm==3 )then
          call ocean_fou_lintranV2(u0, uv, wspeed, lintran_instance%grd%mu, bdrf_ms)
       endif
       call ocean_ss(u0, uv, phi, wspeed, bdrf_ss_cst)
    endif

    if ( flag%rtm==2 .or. flag%rtm==3) then
       call lintran_allocate(nrt, maxd, atm, drv, set, ierr)   
       if (ierr .ne. 0) return 
    endif

    !**** Loop over all wavelength (line-by-line) 
    do iwave = 1, nwave
       albedo = albedo_array(iwave)
       if (glintflag==1) then
          bdrf_ss(1) =  albedo_array(iwave) + bdrf_ss_cst(1)
       else
          bdrf_ss(1) = albedo_array(iwave)
          bdrf_ms(1,1,:,:,0) = albedo_array(iwave)     
       endif
       do l = 1, nrt
          taua(l) = taua_mol(iwave, l) + taua_aer(iwave, l)
          taus(l) = taus_mol(iwave, l) + taus_aer(iwave, l)
          z_ss(1:nstokes, l) = z_ss_int(1:nstokes, l, iwave)
       enddo


       !***** Interface aerosol and molecular quantities to RTM
       call optic_all_der( &
            iwave, &
            aerosol_optic, aerosol, nder, nrt, maxcoefs, &
            taua_aer, taus_aer, plmom_aer, &
            taua_mol, taus_mol, plmom_ray, &
            taua, taus, plmom_in, &
            ierr)

       !***** Call forward-adjoint perturbation RTM for multiple scattering
       !***** output @ iwave: rint=ms reflectance,
       !***** dtaua=deriv. wrt. absorption optical depth,
       !***** dtaus=deriv. wrt. scattering optical depth, 
       !***** dphase=deriv. wrt. phase matrix parameters,
       !***** dalb=deriv. wrt. albedo 
       if ( flag%rtm==1 ) then
          call fwd_adj_perturbation(&
               nrt, &
               taua,&
               taus,&
               plmom_in, &
               pm_in,&
               gauss_quad,&
               gsf,&
               nder,&
               maxd,&
               maxcoefs,&
               bdrf_ms, &
               u0,&
               theta,&
               phi,&
               rint, &
               dtaua, &
               dtaus,&
               dphase,&
               dalb, &
               ierr)
          if (ierr.ne.0) return
          !*** Call single-scattering RTM
          !*** output @ iwave: rint_ss=ss reflectance,
          !*** dtaua_ss=deriv. wrt. absorption optical depth,
          !*** dtaus_ss=deriv. wrt. scattering optical depth, 
          !*** dalb_ss=deriv. wrt. albedo 
          call SINGLE_SCAT(&
               nrt, &
               taua, &
               taus, &
               bdrf_ss, &
               u0, &
               uv, &
               z_ss, &
               rint_ss, &
               dtaus_ss, &
               dtaua_ss, &
               dz_ss, &
               dalb_ss, &
               nder, &
               maxd)
          !*** Add single and multiple scattering contribution to yield the modelled reflectance
          rint_fine(iwave,1:nstokes) = rint(1:nstokes)+rint_ss(1:nstokes)
       elseif ( flag%rtm==2 ) then
          call lintranv2_assign(nstokes, nrt, maxd, nder, cloudflag, maxcoefs, plmom_in, bdrf_ms, bdrf_ss, &
               zss_lin(:,:,:,iwave), surf_emi, taua, taus, atm, drv, set ) 

          set%deltam = .true. 
          execution = (/.true.,.true.,.true./) ! full solution
          if (deriv_flag) then
             call lintran_calculate(execution, atm, set, lintran_instance, rint_lintran, drv, ierr)
             call lintranv2_return( nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua, dtaus, dphase, dalb, drv )
          else
             call lintran_calculate(execution, atm, set, lintran_instance, rint_lintran, ierr)
          endif
          rint_fine(iwave,1:nstokes) = rint_lintran(1:nstokes,1)
       elseif ( flag%rtm==3 ) then
          call lintranv2_assign(nstokes, nrt, maxd, nder, cloudflag, maxcoefs, plmom_in, bdrf_ms, bdrf_ss, &
               zss_lin(:,:,:,iwave), surf_emi, taua, taus, atm, drv, set)
          if (deriv_flag) then
             call lintran_nakajima_multiscattering(atm, set, lintran_instance, rintms, drv, ierr)
             call lintranv2_return( nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua, dtaus, dphase, dalb, drv )  
             if (nstokes==1) then
                call remotec_single_scat_with_derivatives_nst1(atm, set, lintran_instance, rintss, drv)
             elseif(nstokes==3) then   
                call remotec_single_scat_with_derivatives_nst3(atm, set, lintran_instance, rintss, drv)
             endif
             call lintranv2_return( nstokes, nrt, maxd, dz_ss, dalb_ss, dtaua_ss, dtaus_ss, dphase_tmp, dalb_tmp, drv )
          else
             call lintran_nakajima_multiscattering(atm, set, lintran_instance, rintms, ierr)
             if (nstokes==1) then
                call remotec_single_scat_without_derivatives_nst1(atm, lintran_instance, rintss)                
             elseif(nstokes==3) then   
                call remotec_single_scat_without_derivatives_nst3(atm, lintran_instance, rintss)
             endif
          endif
          rint_fine(iwave,1:nstokes) = rintms(1:nstokes,1) + rintss(1:nstokes,1)
       endif

       !**** Combine ss derivatives of reflectance wrt. tau with derivatives of tau wrt. parameters
       if (deriv_flag) then
          !*** Aerosol parameters
          do i = 1, ntype_aer	    
             i1 = nder(1)
             i2 = nder(maxd)
             if((aerosol(i)%CirrusFlag==1 .and. MinCOTFlag==0) .or. (aerosol(i)%CirrusFlag .ne. 1 .and. MinAOTFlag==0)) then
                if (iwave<imid) then

                   aerosol_optic(i)%dzss_all = aerosol_optic(i)%dzss_all_1 + ( aerosol_optic(i)%dzss_all_mid - aerosol_optic(i)%dzss_all_1 )/ (wavelength(imid)-wavelength(1))*(wavelength(iwave)-wavelength(1))
                else
                   aerosol_optic(i)%dzss_all = aerosol_optic(i)%dzss_all_mid + ( aerosol_optic(i)%dzss_all_nwave - aerosol_optic(i)%dzss_all_mid )/ (wavelength(nwave)-wavelength(imid))*(wavelength(iwave)-wavelength(imid))

                endif

                call combine_der(&
                     dtaua(:, i1:i2),&
                     dtaus(:, i1:i2) ,&
                     dphase,&
                     aerosol_optic(i)%dtaua_all,&
                     aerosol_optic(i)%dtaus_all, &
                     aerosol_optic(i)%dphase_all,&
                     maxcoefs,&
                     k_tmp_phase,&
                     k_tmp_taua,&
                     k_tmp_taus,&
                     nder,&
                     maxd,&
                     aerosol(i)%aerosolflags,&
                     npar)
                forall (i_st=1:nstokes, k = 1:npar_mie)
                   !*** Linearly correct for spectral dependance of dcsca/dx_aer and dcabs/dx_aer:
                   k_tmp(i_st, k) =  k_tmp_phase(i_st, k) + &
                        k_tmp_taua(i_st, k)*( 1D0 + aerosol_optic(i)%dcabs_dwave(k)* &
                        ( wavelength(iwave) - wave_mid ) ) + & 
                        k_tmp_taus(i_st, k)*( 1D0 + aerosol_optic(i)%dcsca_dwave(k)* &
                        ( wavelength(iwave) - wave_mid ) )
                end forall
                forall (i_st=1:nstokes, k = npar_mie+1:npar)
                   !k_tmp(i_st, k) =  k_tmp_phase(i_st, k) +  k_tmp_taua(i_st, k) +  k_tmp_taus(i_st, k)
                   k_tmp(i_st, k) =  k_tmp_phase(i_st, k) + &  
                        k_tmp_taua(i_st, k) / sum(aerosol_optic(i)%taua_aer(imid,:))* sum(aerosol_optic(i)%taua_aer(iwave,:)) + &   !spectral dependence of derivatives to height parameters
                        k_tmp_taus(i_st, k) / sum(aerosol_optic(i)%taus_aer(imid,:))* sum(aerosol_optic(i)%taus_aer(iwave,:))       


                end forall
                call combine_der_ss(&
                     dtaua_ss(:,i1:i2),&
                     dtaus_ss(:,i1:i2),&
                     dz_ss,&
                     aerosol_optic(i)%dtaua_all,&
                     aerosol_optic(i)%dtaus_all, &
                     aerosol_optic(i)%dzss_all,&
                     k_tmp_phase,&
                     k_tmp_taua,&
                     k_tmp_taus,&
                     nder,&
                     maxd,&
                     aerosol(i)%aerosolflags,&
                     npar)
                forall (i_st=1:nstokes, k = 1:npar_mie)
                   !*** Linearly correct for spectral dependance off dcsca/dx_aer and dcabs/dx_aer:
                   k_tmp_ss(i_st, k) =  k_tmp_phase(i_st, k) + &
                        k_tmp_taua(i_st, k)*( 1D0 + aerosol_optic(i)%dcabs_dwave(k)* &
                        ( wavelength(iwave) - wave_mid ) ) + & 
                        k_tmp_taus(i_st, k)*( 1D0 + aerosol_optic(i)%dcsca_dwave(k)* &
                        ( wavelength(iwave) - wave_mid ) )
                end forall
                forall (i_st=1:nstokes, k = npar_mie+1:npar)
                   ! k_tmp_ss(i_st, k) =  k_tmp_phase(i_st, k) +  k_tmp_taua(i_st, k) +  k_tmp_taus(i_st, k)
                   k_tmp_ss(i_st, k) =      k_tmp_phase(i_st, k) + k_tmp_taua(i_st, k) / sum(aerosol_optic(i)%taua_aer(imid,:))* sum(aerosol_optic(i)%taua_aer(iwave,:)) + &   !spectral dependence of derivatives to height parameters
                        k_tmp_taus(i_st, k) / sum(aerosol_optic(i)%taus_aer(imid,:))* sum(aerosol_optic(i)%taus_aer(iwave,:))
                end forall

             else
                forall (i_st=1:nstokes, j=1:npar)
                   k_tmp(i_st, j) = 0.d0
                   k_tmp_ss(i_st, j) = 0.d0
                end forall
             endif
             forall (i_st=1:nstokes, j=1:npar)
                deriv_rt%aerosol(iwave, j, i, i_st) = k_tmp(i_st, j) +  k_tmp_ss(i_st, j)
             end forall
          enddo ! ntype_aer

          !*** Albedo      
          deriv_rt%alb(iwave, 1, 1:nstokes) = dalb(1:nstokes) + dalb_ss(1:nstokes)

          !*** Absorbers
          do imol = 1, nmol            
             !*** For oxygen the derivative of reflectance wrt. to the oxygen subcolumn
             !*** is calculated assuming constant mixing ratio (but variable air column)
             if(win_ini%type_x(imol)==7) then
                do k = 1, nrt
                   deriv_rt%densmol(iwave, k, imol, 1:nstokes)= &
                        (dtaua(1:nstokes,k) + dtaua_ss(1:nstokes, k))*cross_mol(iwave,k,imol) + &
                        (dtaus(1:nstokes, k) + dtaus_ss(1:nstokes, k)*csray(iwave)/relo2) 
                enddo
                !*** For absorbers other than oxygen the derivative of reflectance wrt. to the absorber subcolumn
                !*** is calculated assuming constant air column (but variable mixing ratio)
             else
                do i_st = 1, nstokes
                   do k = 1, nrt
                      deriv_rt%densmol(iwave, k, imol, i_st) = &
                           (dtaua(i_st, k) + dtaua_ss(i_st,k))*cross_mol(iwave, k, imol)          
                   enddo
                enddo
             endif
          enddo
          !*** Temperature offset
          if(flag%temp==1) then
             do i = 1, nstokes
                deriv_rt%T(iwave,i) =  sum((dtaua(i,:)+dtaua_ss(i,:))*dtaua_T(iwave,:))
             enddo
          endif
          !*** O2 column/pressure
          if(flag%O2==1) then
             do i = 1, nstokes
                do k = 1, nrt
                   deriv_rt%P(iwave,k,i) = (dtaua(i,k) + dtaua_ss(i,k))*dtaua_P(iwave,k)
                enddo
             enddo
          endif

       endif

    enddo
    !**** End loop over all wavelengths

    !**** Clean up all lintran instances and internals
    if ( flag%rtm ==2 .or. flag%rtm==3) then
       call lintran_deallocate(atm, set, drv, ierr)
       call lintran_close(lintran_instance,ierr)
    endif

    deallocate(aerosol_optic, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'RAD_TRAN_LBL: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine rad_trans_lbl ! }}}

  !------------------------------------------------------------------------------    
  !> Combine derivatives of rt-model with derivatives of Mie calculations  
  !------------------------------------------------------------------------------
  subroutine COMBINE_DER( & ! {{{
       DTAUA,      DTAUS,     DPHASE, &
       DTAUA_ALL,  DTAUS_ALL, DPHASE_ALL, &
       NC_TRUN, K_WAVE_phase, K_WAVE_taua, K_WAVE_taus, NDER,   MAXD, &
       aerosolflags, NDIM_DER)
    !*** Input
    real(double), intent(in) :: &
         DTAUS_ALL(MAXD,NDIM_DER), &
         DTAUA_ALL(MAXD,NDIM_DER), &
         DPHASE_ALL(NPER,0:MAXLEG,MAXD,NDIM_DER)
    real(double), intent(in) ::  &        !derivatives
         DTAUA(NSTOKES,MAXD), &               !w.r.t. OMEGA
         DTAUS(NSTOKES,MAXD), &               !w.r.t. OMEGA
         DPHASE(NSTOKES,0:MAXLEG,NPER,MAXD)  !w.r.t. PHASE
    integer, dimension(ndim_der), intent(in) :: aerosolflags 
    integer, intent(in) ::  MAXD, ndim_der
    integer, dimension(:), intent(in) ::  NC_TRUN, NDER
    !*** Output
    real(double), intent(out) :: K_WAVE_phase(NSTOKES,NDIM_DER)
    real(double), intent(out) :: K_WAVE_taua(NSTOKES,NDIM_DER)
    real(double), intent(out) :: K_WAVE_taus(NSTOKES,NDIM_DER)
    !*** Local variables
    integer :: i_st, i, k, ik, l, n
    real(double) :: term1, term2, term3
    !------------------------------------------------------------------------------
    K_WAVE_phase = 0.d0
    K_WAVE_taua  = 0.d0
    K_WAVE_taus  = 0.d0

    do I_ST = 1, NSTOKES
       do I = 1, NDIM_DER
          if(aerosolflags(i).ne.0) then
             TERM1 = 0.
             TERM2 = 0.
             TERM3 = 0.
             do K = 1, MAXD
                IK = NDER(K)
                do L = 0, NC_TRUN(IK)
                   do N = 1, NPER
                      TERM1 = TERM1 + &
                           DPHASE(I_ST,L,N,K)* &
                           DPHASE_ALL(N,L,K,I)
                   enddo
                enddo
                TERM2 = TERM2 + DTAUA(I_ST,K)*DTAUA_ALL(K,I)
                TERM3 = TERM3 + DTAUS(I_ST,K)*DTAUS_ALL(K,I)   
             enddo   
             K_WAVE_phase(I_ST,I) = TERM1 !K-matrix dI/daer_params: contributions through dI/dphase function only. 
             K_WAVE_taua(I_ST,I) =  TERM2 !K-matrix dI/daer_params: contributions through dI/dtaua only. 
             K_WAVE_taus(I_ST,I) =  TERM3 !K-matrix dI/daer_params: contributions through dI/dtaus only.      
          endif
       enddo
    enddo
    return
  end subroutine combine_der ! }}}

  !------------------------------------------------------------------------------    
  !> Combine derivatives of rt-model with derivatives of Mie calculations  
  !! for single scattering
  !------------------------------------------------------------------------------
  subroutine COMBINE_DER_SS( & ! {{{
       DTAUA_SS,      DTAUS_SS,     DZ_SS, &
       DTAUA_ALL,     DTAUS_ALL,    DZSS_ALL, &
       K_WAVE_phase, K_WAVE_taua, K_WAVE_taus, NDER,         MAXD, &
       aerosolflags, NDIM_DER)
    !*** Input
    integer, intent(in) :: ndim_der
    real(double), intent(in) :: &
         DTAUS_ALL(MAXD,NDIM_DER), &
         DTAUA_ALL(MAXD,NDIM_DER), &
         DZSS_ALL(nstokes,maxd,ndim_der)
    real(double), intent(in) :: &          !derivatives
         DTAUA_SS(NSTOKES,MAXD), &             !w.r.t. OMEGA
         DTAUS_SS(NSTOKES,MAXD), &             !w.r.t. OMEGA
         DZ_SS(NSTOKES,MAXD)                  !w.r.t. PHASE
    integer,dimension(ndim_der), intent(in) :: aerosolflags
    integer, intent(in) ::  MAXD,  NDER(MAXD)
    !*** Output
    real(double), intent(out) :: K_WAVE_phase(NSTOKES,NDIM_DER)
    real(double), intent(out) :: K_WAVE_taua(NSTOKES,NDIM_DER)
    real(double), intent(out) :: K_WAVE_taus(NSTOKES,NDIM_DER)
    !*** Local variables
    integer :: i_st, i, k, ik
    real(double) :: term1, term2, term3
    !-----------------------------------------------------------------------
    K_WAVE_phase = 0.d0
    K_WAVE_taua  = 0.d0
    K_WAVE_taus  = 0.d0

    do I_ST = 1, NSTOKES
       do I = 1, NDIM_DER
          if(aerosolflags(i).ne.0) then
             TERM1 = 0.d0
             TERM2 = 0.d0
             TERM3 = 0.d0
             do K = 1, MAXD
                IK = NDER(K)
                TERM1 = TERM1 + DZ_SS(I_ST,K)*DZSS_ALL(I_ST,K,I)
                TERM2 = TERM2 + DTAUA_SS(I_ST,K)*DTAUA_ALL(K,I)
                TERM3 = TERM3 + DTAUS_SS(I_ST,K)*DTAUS_ALL(K,I)
             enddo
             K_WAVE_phase(I_ST,I) = TERM1
             K_WAVE_taua(I_ST,I)  = TERM2
             K_WAVE_taus(I_ST,I)  = TERM3
          endif
       enddo
    enddo

    return
  end subroutine combine_der_ss ! }}}

  !-----------------------------------------------------------------------
  !> Nakajima correction for strongly forward peaked phase functions
  subroutine nakajima_correction(phase_in, phase, nder, maxd, nrt)
    !*** input
    integer, intent(in) :: maxd, nder(maxd), nrt
    real(double), intent(in) :: phase_in(nstokes, nstokes, 0:maxleg, nrt)
    !*** ouput
    real(double), intent(out) :: phase(nstokes, nstokes, 0:maxleg, nrt)
    !*** local variables
    integer :: ik, k, l, i_st ! Iterators.
    real(double) :: f ! Forward peak.
    !-------------------------------------------------------------------------

    ! Start with the original phase coefficients. Layers outside the aerosol layers
    ! will not get a Nakajima correction. Those are layers with just Rayleigh scattering
    ! and the highest nonzero Legendre coefficient in those layers is 2.
    phase = phase_in

    ! The plmom for the scattering layer will be overwritten.
    do ik = 1,maxd

       ! Go to the right layer.
       k = nder(ik)
       ! Derive the forward peak.
       f = phase_in(1,1,maxstr,k)/(2.*dble(maxstr) + 1.d0)

       do l = 0,maxstr-1

          ! Diagonal terms.
          do i_st = 1,nstokes
             phase(i_st,i_st,l,k) = phase_in(i_st,i_st,l,k)/(1.d0-f)-(f*(2.d0*dble(l)+1.d0))/(1.d0-f)
          enddo

          ! Off-diagonal terms (HH: this notation is needed to commpile with debug flags)
          do i_st = 1, nstokes-1, 2
             phase(i_st,i_st+1,l,k) = phase_in(i_st,i_st+1,l,k) / (1.d0-f)
             if (i_st==1) phase(i_st+1,i_st,l,k) = phase(i_st,i_st+1,l,k)
             if (i_st==3)  phase(i_st+1,i_st,l,k) = -phase(i_st,i_st+1,l,k)
          enddo
!!$          ! Off-diagonal terms.(Original notation)
!!$          ! If there is polarization we have the elements I-Q, which are equal to
!!$          ! elements Q-I.
!!$          if (nstokes .ge. 3) then ! nstokes .eq. 2 is nonsense.
!!$             phase(2,1,l,k) = phase_in(2,1,l,k) / (1.d0-f)
!!$             phase(1,2,l,k) = phase(2,1,l,k)
!!$          endif
!!$
!!$          ! Elements (3,4) and (4,3) in the highly unlikely case that V-polarization is
!!$          ! taken into account.
!!$          if (nstokes .eq. 4) then
!!$             phase(3,4,l,k) = phase_in(3,4,l,k) / (1.d0-f)
!!$             phase(4,3,l,k) = -phase(3,4,l,k)
!!$          endif

       enddo

    enddo

    return

  end subroutine nakajima_correction
  !-----------------------------------------------------------------------
end module rad_trans_module
