!------------------------------------------------------------------------------
!> @todo error handling
!> @todo quality flags
!------------------------------------------------------------------------------
module profile_inversion_module
  use header_module
  use forward_model_module, only: &
       absorbers, forward_model_hi, forward_model_lo,  &
       derivatives, &
       Mie_lut, cirrus_table, read_aerosol_netcdf, read_cirrus_netcdf, &
       aero, set_altdis, get_aerosol_properties_lognormal, &
       instrument_response, spectral_response_stored, &
       window_ini, settings_flags, file_paths, read_settings, read_win_xsdb, window_spectrum
  use forward_model_noscat_module
  use pt_regularization_module
  use adhoc_inversion_module
  implicit none
  private

  !*** Types
  public :: retrieval_data, absorbers, aero, Mie_lut, cirrus_table, window_ini, settings_flags, file_paths, window_spectrum

  !*** Procedures
  public :: profile_inversion, read_aerosol_netcdf, read_cirrus_netcdf, get_aerosol_properties_lognormal
  public :: read_settings, read_win_xsdb
  private :: init_state_vector, update_forward_model, concat_windows

  !------------------------------------------------------------------------------
  !> Output data of retrieval algorithm
  type :: retrieval_data
     real(double), dimension(:,:), allocatable :: x_state
     real(double), dimension(:), allocatable :: dvair
     real(double), dimension(:,:), allocatable :: ak
     real(double), dimension(:), allocatable :: x_apr
     real(double), dimension(:, :), allocatable :: cf
     real(double), dimension(:,:), allocatable :: s_state

     real(double), dimension(:), allocatable :: chi2_window    ! chi2 per window
     real(double), dimension(:), allocatable :: ot
     real(double), dimension(:), allocatable :: cot
     real(double), dimension(:), allocatable :: albedo
     real(double), dimension(:), allocatable :: albedo_err
     real(double), dimension(:), allocatable :: spectral_shift ! Spectral shift [nm]
     real(double) :: Fs                                        ! Fluorescence intensity

     real(double), dimension(:), allocatable :: spectrum_mod      ! Modeled relfectance/radiance
     real(double), dimension(:), allocatable :: spectrum_meas     ! Measured reflectance/radiance
     real(double), dimension(:), allocatable :: spectrum_meas_cov ! Covariance of measured reflectance/radiance
     real(double), dimension(:), allocatable :: wavelength        ! Wavelength for reflectance/radiance

     integer, dimension(:), allocatable :: type_x_target       ! Target absorbers
     character*25, dimension(:), allocatable :: x_state_name   ! State vector element identifier

     integer, dimension(:), allocatable :: ny   !Number of spectral points used for retrieval per band

     real(double) :: dfs
     real(double), dimension(:), allocatable :: dfs_target
     real(double) :: dfs_scat
     real(double), dimension(:), allocatable :: chi2
     real(double), dimension(:), allocatable :: lambda
     real(double), dimension(:), allocatable :: vza, raa ! viewing zenith angle, relative azimuth angle

     real(double) :: air_col_old
     real(double), dimension(:), allocatable :: p, z, t
     real(double), dimension(:,:), allocatable :: s_apr
     real(double) :: rms

     integer :: convergence
     integer :: iter
     integer :: error_id

     character(stringlen) :: iterflag

  end type retrieval_data

contains
  !------------------------------------------------------------------------------
  !> @details Iterative solver
  !! Initialize state vector, call forward model to compute reflectance and derivatives,
  !! do inversion, update state vector, repeat
  !------------------------------------------------------------------------------
  subroutine profile_inversion( &
       aero_lut, cirrus_lut,  &
       flag, glintflag, &
       Fs_apr, &
       wspeed, &
       nlay, &
       absorb, atm_rt, atm_xs, dvair, &
       response, &
       win_ini, win, aerosol, retrieval_output, ierr, line_number)
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    type(settings_flags), intent(in) :: flag
    integer, intent(in) :: glintflag, nlay
    real(double), intent(in) :: FS_apr, wspeed
    type(absorbers), intent(in) :: absorb
    type(atmosphere), intent(in) :: atm_rt
    type(instrument_response), dimension(:), intent(in) :: response
    type(window_ini), dimension(:), intent(in) :: win_ini
      integer, intent(in) :: line_number  ! hack
    !*** Input/output
    type(atmosphere), intent(inout) :: atm_xs
    real(double), dimension(:), intent(inout) :: dvair      ! Partial air column, subject to change in O2 retrieval (Dim: natm)
    type(window_spectrum), dimension(:), intent(inout) :: win
    type(aero), dimension(:), intent(inout) :: aerosol
    !*** Output
    type(retrieval_data), intent(out) :: retrieval_output
    integer, intent(out) :: ierr
    !*** local variables
    integer :: ExitXSFlag, natm, nwin, ntype_aer
    integer :: i, j, l, n, SVDflag, nlsq, i1, i2
    integer :: aer_red, reduction, reduce_i, reduce_j
    integer :: ExitFlag
    integer, parameter :: maxiter = 30, miniter = 5
    integer :: iter, convergence, slowconvflag
    real(double) :: residual
    real(double) :: state_stop
    real(double) :: lambda
    real(double), parameter ::  minaot = 1.D-20, mincot = 1.D-20   ! HH: from GOSAT (used to be 1.d-3)
    real(double) :: chi2old
    real(double) :: chi2min
    real(double) :: dfs_min
    integer, dimension(:),allocatable :: red_positions
    real(double), dimension(:), allocatable :: x_state_old           ! State vector previous it
    real(double), dimension(:), allocatable :: x_state_min           ! State vector minimum chi2
    real(double), dimension(:), allocatable :: upperx                ! State vector upper boundaries
    real(double), dimension(:), allocatable :: lowerx                ! State vector lower boundaries
    real(double), dimension(:,:), allocatable :: ak_min              ! Averaging kernel with minimum chi2
    real(double), dimension(:,:),allocatable :: call_derivatives_lo  ! Substitution array for derivatives_lo
    real(double), dimension(:),allocatable :: call_x_state           ! Substitution array for x_state
    real(double), dimension(:,:),allocatable :: call_s_state         ! Substitution array for s_state
    real(double), dimension(:),allocatable :: call_x_apr             ! Substitution array for x_apr
    real(double), dimension(:,:),allocatable :: call_ak              ! Substitution array for ak
    real(double), dimension(:,:),allocatable :: call_cf              ! Substitution array for cf
    real(double), dimension(:),allocatable :: call_upperx            ! Substitution array for upperx
    real(double), dimension(:),allocatable :: call_lowerx            ! Substitution array for lowerx
    real(double), dimension(:,:), allocatable :: reflectance_hi
    real(double), dimension(:, :), allocatable :: derivatives_lo     ! Modelled derivatives of the Log of the reflectance concatenated over all windows
    real(double), dimension(:,:), allocatable :: derivatives_lo_old  ! Modelled derivatives of the Log of previous it reflectance concatenated over all windows
    real(double), dimension(sum(win(:)%nwave_lo)) :: regpix          ! Array of pixels to be regularized
    integer, dimension(:), allocatable :: regskill                   ! Array of skill-IDs for regularization
    real(double), dimension(:), allocatable :: ymeas, ymod, ycov, ymod_old, ycov_unscaled
    real(double) :: degfreedom                            ! Degrees of freedom
    real(double) :: dfs_scat                              ! Degrees of freedom for scattering parameters
    real(double), dimension(:),allocatable :: dfs_target  ! Degrees of freedom for target vertical profiles
    type(derivatives), dimension(:), allocatable :: deriv_hi
    integer, dimension(:), allocatable :: nder
    real(double) :: chi2, rms
    real(double) :: chi2_final, covmax
    real(double) :: s1, s2, s3 				!rrae: Multiplicative stokes coefficients
    real(double), dimension(:), allocatable :: covrm
    integer :: error_ID, nwave_lo
    real(double), dimension(:), allocatable :: x_state    ! State vector to be retrieved (Dim: nstate)
    real(double), dimension(:), allocatable :: x_apr      ! Apriori/initial guess vector (Dim: nstate)
    real(double), dimension(:,:), allocatable :: s_state  ! State vector covariance matrix (Dim: nstate,nstate)
    real(double), dimension(:,:), allocatable :: s_apr
    real(double), dimension(:,:), allocatable :: s_y      !Diagonal measurement covariance matrx (Dim: nwave_lo,nwave_lo)
    real(double), dimension(:,:), allocatable :: ak       ! Averaging kernel matrix (Dim: nstate,nstate)
    real(double), dimension(:,:), allocatable :: cf       ! Contribution function matrix (Dim: nstate,nwave_lo)
    real(double), dimension(:,:), allocatable :: cf_min   ! Contribution function matrix minimum chi2
    real(double), dimension(:), allocatable :: play_old   ! Pressure, layer center (Dim: natm)
    real(double), dimension(:), allocatable :: tlay_old   ! Temperature, layer center (Dim: natm)
    real(double), dimension(:), allocatable :: dvair_old  ! Partial air column (Dim: natm)
    real(double), dimension(:), allocatable ::  vmr_h2o
    integer :: MinAOTFlag       ! No aerosol in retrieval if low aerosol optical thickness
    integer :: MinCOTFlag       ! No cirri in retrieval if low cirrus optical thickness
    integer:: maxotflag
    integer :: boundary_flag
    integer :: naux, nstate, naer, off, k, io
    character(stringlen) :: message
    character*25, dimension(:), allocatable :: x_state_name  ! State vector element identifier
    character*2 :: ch
    character(stringlen) :: iterflag                        ! Iteration anomaly flag
    !---------------------------------------------------------------------------------------------------------
    if(flag%output >= 2) then
       write(message, '(a)') '*** Start of profile_inversion ***'
       call writelog(message, 1)
    endif

    !*** Initialize error identifier
    ierr = 0
    maxotflag = 0

    !*** Get dimensions
    nwin = size(win_ini)
    ntype_aer = size(aerosol)
    nwave_lo = sum(win(:)%nwave_lo)
    naux = absorb%ntype_global
    natm = atm_xs%n

    do n = 1, nwin
       naux = naux + win_ini(n)%albflag + abs(win_ini(n)%IOffFlag) + max(0,win_ini(n)%Fsflag) &
            + abs(win_ini(n)%spsh0flag) + abs(win_ini(n)%spsh1flag) + abs(win_ini(n)%spsh2flag) &
            + abs(win_ini(n)%sunsh0flag)
    enddo
    naer = 0
    do n = 1, ntype_aer
       naer = naer + sum(aerosol(n)%AerosolFlags)
    enddo
    naux = naux +  flag%temp + naer
    nstate = nlay*absorb%ntype_target + naux

    !*** Allocate
    allocate(x_state_old(nstate), &
         x_state_min(nstate), &
         upperx(nstate), &
         lowerx(nstate), &
         ak_min(nstate,nstate), &
         x_state(nstate),&
         x_state_name(nstate),&
         x_apr(nstate), &
         s_state(nstate,nstate), &
         s_apr(nstate,nstate), &
         s_y(nwave_lo,nwave_lo), &
         ak(nstate,nstate), &
         cf(nstate,nwave_lo), &
         cf_min(nstate,nwave_lo), &
         derivatives_lo(nwave_lo,nstate), &
         derivatives_lo_old(nwave_lo,nstate), &
         regskill(nstate), &
         ymeas(nwave_lo),&
         ymod(nwave_lo), &
         ymod_old(nwave_lo), &
         ycov(nwave_lo), &
         ycov_unscaled(nwave_lo), &
         dfs_target(absorb%ntype_target),&
         play_old(natm), tlay_old(natm), dvair_old(natm), vmr_h2o(natm), &
         deriv_hi(nwin), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'PROFILE_INVERSION: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    if (.not. allocated(retrieval_output%lambda)) then
       allocate(retrieval_output%lambda(0:maxiter), &
            retrieval_output%x_state(nstate,0:maxiter), &
            retrieval_output%chi2(0:maxiter), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
    endif

    !*** Allocate arrays with fixed dimensions during a run
    if (.not. allocated(retrieval_output%ak)) then
       allocate(retrieval_output%ak(nstate, nstate), &
            retrieval_output%s_state(nstate, nstate), &
            retrieval_output%x_apr(nstate), &
            retrieval_output%dvair(nlay), &
            retrieval_output%dfs_target(absorb%ntype_target), &
            retrieval_output%x_state_name(nstate), &
            retrieval_output%type_x_target(absorb%ntype_target), &
            retrieval_output%chi2_window(nwin), &
            retrieval_output%ny(nwin), &
            retrieval_output%s_apr(nstate,nstate), &
            retrieval_output%ot(nwin), &
            retrieval_output%cot(nwin), &
            retrieval_output%albedo(nwin), &
            retrieval_output%albedo_err(nwin), &
            retrieval_output%spectral_shift(nwin), &
            retrieval_output%vza(nwin), &
            retrieval_output%raa(nwin), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
    endif

    !*** Allocate arrays with possibly variable dimensions during a run
    if (.not. allocated(retrieval_output%cf)) then
       allocate(retrieval_output%cf(nstate, nwave_lo), &
            retrieval_output%spectrum_mod(nwave_lo), &
            retrieval_output%spectrum_meas(nwave_lo), &
            retrieval_output%spectrum_meas_cov(nwave_lo), &
            retrieval_output%wavelength(nwave_lo), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
    elseif(size(retrieval_output%cf(1,:)).ne. nwave_lo) then
       deallocate(retrieval_output%cf, &
            retrieval_output%spectrum_mod, &
            retrieval_output%spectrum_meas, &
            retrieval_output%spectrum_meas_cov, &
            retrieval_output%wavelength, &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory deallocation error'
          ierr = ierr_deall
          goto 999
       endif
       allocate(retrieval_output%cf(nstate, nwave_lo), &
            retrieval_output%spectrum_mod(nwave_lo), &
            retrieval_output%spectrum_meas(nwave_lo), &
            retrieval_output%spectrum_meas_cov(nwave_lo), &
            retrieval_output%wavelength(nwave_lo), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
    endif

    retrieval_output%x_state = nf_fill_double
    retrieval_output%lambda = nf_fill_double
    retrieval_output%chi2 = nf_fill_double

    !*** Initialize variables
    chi2 = 0.d0
    rms = 0.d0
    tlay_old = atm_xs%t
    play_old = atm_xs%p
    dvair_old = dvair
    error_ID = 0 ! a priori no error
    nlsq = 0     ! do nlsq unconstrained least square steps before doing constrained PT-steps
    degfreedom = 0.
    boundary_flag = 0
    state_stop = 1.D-4
    derivatives_lo_old = 0.D0
    chi2old = INF
    chi2min = INF
    chi2_final = INF
    chi2 = INF/10.d0
    ExitFlag = 0
    MinAOTFlag = 0
    MinCOTFlag = 0
    SVDFlag = 0
    slowconvflag = 0
    dfs_target=0.D0
    dfs_scat=0.D0

    !*** Instrument covariance (unscaled)
    off = 0
    if (flag%fit == 1) then !*** fit reflectance
       do n = 1, nwin
          do k = 1, win(n)%nwave_lo
             ycov_unscaled(k+off) = win(n)%spectrum_cov(k)/ &
                  (win(n)%sun_spectrum_sat_lo(k)*win(n)%sun_spectrum_sat_lo(k))
          enddo
          off = off + win(n)%nwave_lo
       enddo
    elseif (flag%fit == 2) then !*** fit radiance
       do n = 1, nwin
          do k = 1, win(n)%nwave_lo
             ycov_unscaled(k+off) = win(n)%spectrum_cov(k)
          enddo
          off = off + win(n)%nwave_lo
       enddo
    endif

    do n = 1, nwin
       !*** Allocate arrays with fixed dimensions during a run
       if( .not. allocated(win(n)%x_molec)) then
          allocate(win(n)%x_molec(natm,win_ini(n)%ntype), &
               win(n)%albedo(max(1,maxval(win_ini(:)%albflag))), &
               win(n)%wavelength_hi_new(win_ini(n)%nwave_hi), &
               win(n)%sun_spectrum_ref_hi(win_ini(n)%nwave_hi), &
               win(n)%IOff(abs(win_ini(n)%IOffFlag)), &
               stat = ierr)
          if (ierr .ne. 0) then
             write(message, *) 'PROFILE_INVERSION: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       endif
       !*** Allocate arrays with possibly variable dimensions during a run
       if( .not. allocated(win(n)%wavelength_lo_new)) then
          allocate(win(n)%wavelength_lo_new(win(n)%nwave_lo), &
               win(n)%reflectance_lo(win(n)%nwave_lo), &
               win(n)%reflectance_meas(win(n)%nwave_lo), &
               win(n)%reflectance_meas_cov(win(n)%nwave_lo), &
               win(n)%derivatives_lo(win(n)%nwave_lo, nlay, win_ini(n)%ntype), &
               win(n)%derivT_lo(win(n)%nwave_lo), &
               win(n)%derivP_lo(win(n)%nwave_lo,nlay), &
               win(n)%deriv_albedo_lo(win(n)%nwave_lo,win_ini(n)%albflag), &
               win(n)%deriv_specshift0_lo(win(n)%nwave_lo), &
               win(n)%deriv_specshift1_lo(win(n)%nwave_lo), &
               win(n)%deriv_specshift2_lo(win(n)%nwave_lo), &
               win(n)%deriv_sunshift0_lo(win(n)%nwave_lo), &
               win(n)%deriv_aerosol_lo(win(n)%nwave_lo,naer), &
               win(n)%deriv_Ioff_lo(win(n)%nwave_lo,abs(win_ini(n)%IOffFlag)), &
               win(n)%deriv_Fs_lo(win(n)%nwave_lo, max(0,win_ini(n)%Fsflag)), &
               stat = ierr)
          if (ierr .ne. 0) then
             write(message, *) 'PROFILE_INVERSION: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       elseif(size(win(n)%wavelength_lo_new) .ne. win(n)%nwave_lo) then
          deallocate(win(n)%wavelength_lo_new, &
               win(n)%reflectance_lo, &
               win(n)%reflectance_meas, &
               win(n)%reflectance_meas_cov, &
               win(n)%derivatives_lo, &
               win(n)%derivT_lo, &
               win(n)%derivP_lo, &
               win(n)%deriv_albedo_lo, &
               win(n)%deriv_specshift0_lo, &
               win(n)%deriv_specshift1_lo, &
               win(n)%deriv_specshift2_lo, &
               win(n)%deriv_sunshift0_lo, &
               win(n)%deriv_aerosol_lo, &
               win(n)%deriv_Ioff_lo, &
               win(n)%deriv_Fs_lo, &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'PROFILE_INVERSION: memory deallocation error'
             ierr = ierr_deall
             goto 999
          endif
          allocate(win(n)%wavelength_lo_new(win(n)%nwave_lo), &
               win(n)%reflectance_lo(win(n)%nwave_lo), &
               win(n)%reflectance_meas(win(n)%nwave_lo), &
               win(n)%reflectance_meas_cov(win(n)%nwave_lo), &
               win(n)%derivatives_lo(win(n)%nwave_lo, nlay, win_ini(n)%ntype), &
               win(n)%derivT_lo(win(n)%nwave_lo), &
               win(n)%derivP_lo(win(n)%nwave_lo,nlay), &
               win(n)%deriv_albedo_lo(win(n)%nwave_lo,win_ini(n)%albflag), &
               win(n)%deriv_specshift0_lo(win(n)%nwave_lo), &
               win(n)%deriv_specshift1_lo(win(n)%nwave_lo), &
               win(n)%deriv_specshift2_lo(win(n)%nwave_lo), &
               win(n)%deriv_sunshift0_lo(win(n)%nwave_lo), &
               win(n)%deriv_aerosol_lo(win(n)%nwave_lo,naer), &
               win(n)%deriv_Ioff_lo(win(n)%nwave_lo,abs(win_ini(n)%IOffFlag)), &
               win(n)%deriv_Fs_lo(win(n)%nwave_lo, 2), &
               stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'PROFILE_INVERSION: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
       endif

       if (flag%fit==1) then
          win(n)%reflectance_meas = win(n)%spectrum/win(n)%sun_spectrum_sat_lo
          win(n)%reflectance_meas_cov = win(n)%spectrum_cov/(win(n)%sun_spectrum_sat_lo**2)
       elseif (flag%fit==2) then
          win(n)%reflectance_meas = win(n)%spectrum
          win(n)%reflectance_meas_cov = win(n)%spectrum_cov
       endif
       win(n)%wavelength_lo_new = win(n)%wavelength_lo
       win(n)%wavelength_hi_new = win_ini(n)%wavelength_hi
       win(n)%sun_spectrum_ref_hi = win(n)%sun_spectrum_sat_hi
       win(n)%derivatives_lo = 0.D0
       win(n)%deriv_albedo_lo = 0.D0
       win(n)%deriv_specshift0_lo = 0.D0
       win(n)%deriv_specshift1_lo = 0.D0
       win(n)%deriv_specshift2_lo = 0.D0
       win(n)%deriv_sunshift0_lo = 0.D0
       win(n)%deriv_aerosol_lo = 0.D0
       if(win_ini(n)%IOffFlag .ne. 0) then
          win(n)%IOff = 0.d0
          if(win_ini(n)%IOffFlag < 0) win(n)%IOff(1) = 1.d0
          win(n)%deriv_IOff_lo = 0.d0
       endif
       if(win_ini(n)%Fsflag .ne. 0) then
          win(n)%Fs(:) = 0.d0
          if(win_ini(n)%Fsflag<0)win(n)%Fs(1) = Fs_apr
          win(n)%deriv_Fs_lo = 0.d0
       endif
       if (nstokes>1) then
          !***Assign multiplicative stokes coefficients
          s1=win(n)%measurement_stokesc(1)
          s2=win(n)%measurement_stokesc(2)
          s3=win(n)%measurement_stokesc(3)
       else
          s1 = 1.d0
       endif
       !*** Maximum of the measured reflectance is the initial guess for albedo
       win(n)%albedo = 0.D0
       win(n)%albedo(1) = (1/s1)*DABS(maxval(win(n)%spectrum)/maxval(win(n)%sun_spectrum_ref_lo))/cos(win(n)%sza/180.*pi)*pi
       if (glintflag>=1) win(n)%albedo(1)=0.D0

       !*** Absorber number densities: x_molec is updated by the iteration, dv_x remains the initial guess
       win(n)%x_molec = win(n)%dv_x
    enddo ! close loop over windows

    !*** Initialize state vector
    call init_state_vector(flag%temp, glintflag, natm, nlay, absorb, win_ini, win, aerosol, nder, &
         regskill, upperx, lowerx, x_state, x_state_name, s_state, ierr)
    if (ierr .ne. 0) return

    x_state_old = x_state*(1.0 + 2.0*state_stop) !set to different value than x_state to avoid convergence in 1st step
    x_apr = x_state ! prior is set to the first guess here
    s_apr = s_state
    x_state_min = 0.D0
    dfs_min = 0.D0
    ak = 0.0d0
    do j = 1, nstate
       ak(j,j) = 1.0d0
    enddo



    !*** Begin iteration
    iter = 0
    iterflag='NOMINAL'
    do
       iter = iter+1
       if (iter < nlsq+1) lambda = 0.d0   ! no step-size reduction for unconstrained least square steps
       if (iter == nlsq+1)lambda = 10.d0  ! starting value for constrained PT-steps

       !*** Get vmr of H2O in xs-layers
       vmr_h2o = 0.d0
       do n = 1, nwin
          do i = 1, win_ini(n)%ntype
             if(abs(win_ini(n)%xsdb(i)%species)==1 .or. (abs(win_ini(n)%xsdb(i)%species)>=100 .and. abs(win_ini(n)%xsdb(i)%species)<=199))then
                vmr_h2o(:) =  win(n)%x_molec(:,i)/dvair(:) ! caculate vmr of water per layer
!!$                vmr_h2o(:) =  win(n)%dv_x(:,i)/dvair(:) ! Use initial VMR
                goto 101
             endif
          enddo
       enddo
101    continue

       !*** Forward model of the reflectance spectrum and its derivatives
       do n = 1, nwin
          if (allocated(reflectance_hi)) deallocate(reflectance_hi)
          allocate(reflectance_hi(win_ini(n)%nwave_hi, nstokes), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'PROFILE_INVERSION: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
          if(.not. allocated(deriv_hi(n)%densmol)) then
             allocate(deriv_hi(n)%densmol(win_ini(n)%nwave_hi, nlay, win_ini(n)%ntype, nstokes), &
                  deriv_hi(n)%alb(win_ini(n)%nwave_hi, win_ini(n)%albflag, nstokes), &
                  deriv_hi(n)%T(win_ini(n)%nwave_hi, nstokes), &
                  deriv_hi(n)%P(win_ini(n)%nwave_hi, nlay, nstokes), &
                  deriv_hi(n)%aerosol(win_ini(n)%nwave_hi, naer, 1, nstokes), &
                  deriv_hi(n)%Fs(win_ini(n)%nwave_hi,max(0,win_ini(n)%Fsflag)), &
                  stat = ierr)
             if (ierr .ne. 0) then
                write(message, *) 'PROFILE_INVERSION: memory allocation error'
                ierr = ierr_all
                goto 999
             endif
             deriv_hi(n)%densmol = 0.D0
             deriv_hi(n)%alb     = 0.D0
             deriv_hi(n)%T       = 0.D0
             deriv_hi(n)%P       = 0.D0
             deriv_hi(n)%aerosol = 0.D0
             deriv_hi(n)%Fs      = 0.D0
          endif
          if (flag%scat > 0) then
             !*** calculate high resolution and derivatives
             call forward_model_hi(aero_lut, cirrus_lut, &
                  flag, glintflag, &
                  win(n)%sza, win(n)%iza, win(n)%phi, wspeed, &
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
                  deriv_hi(n), &
                  ExitXSFlag, MaxOTFlag, ierr, .true.)
             if(ierr.ne.0)then
                ExitFlag = 1
                error_id = 94
                iterflag='RTM_ER'
                goto 100
             endif
          elseif (flag%scat == 0) then
             call forward_model_hi_noscat( &
                  flag%xs, flag%O2, flag%temp, glintflag, &
                  win(n)%sza, win(n)%iza, win(n)%phi, flag%observer_location, win(n)%observer_height, wspeed, n, &
                  absorb, atm_xs, dvair, dvair_old, vmr_h2o, play_old, &
                  win_ini, win, reflectance_hi, &
                  deriv_hi(n)%densmol(:,:,:,1), deriv_hi(n)%alb(:,:,1), deriv_hi(n)%T(:,1), deriv_hi(n)%P(:,:,1), &
                  ExitXSFlag, ierr)
             if (ierr .ne. 0) return
          endif
          if(ExitXSFlag==1)then
             ExitFlag = 1
             error_id = 95
             iterflag='XSDB_EX'
             goto 100
          endif
          !*** calculate low resolution spectrum (convolution with ILS) and derivatives
          call forward_model_lo( &
               flag, &
               naer, nlay, absorb, &
               response(n), &
               win_ini(n), &
               reflectance_hi, win(n)%reflectance_lo, &
               win(n), &
               deriv_hi(n), &
               ierr)
          if (ierr .ne. 0) return

          !*** CONTROL OUTPUT
          if(flag%output >= 3)then
             write(ch,'(i2.2)')n
             open(newunit(io),FILE='./CONTRL_OUT/spectrum_lores_'//ch//'.dat')
             write(io,'(A)')'# Wavelength / nm'
             write(io,'(A)')'# Radiance measured'
             write(io,'(A)')'# Radiance noise'
             write(io,'(A)')'# Radiance modelled'
             write(io,'(A)')'# Measured - modelled '
             write(io,'(A)')'# Solar irradiance'
             do k = 1, win(n)%nwave_lo
                write(io,'(100(1pE16.8E3,x))') &
                     win(n)%wavelength_lo(k),&
                     win(n)%reflectance_meas(k),&
                     sqrt(win(n)%reflectance_meas_cov(k)),&
                     win(n)%reflectance_lo(k),&
                     win(n)%reflectance_meas(k)-win(n)%reflectance_lo(k),&
                     win(n)%sun_spectrum_ref_lo(k)
             enddo
             close(io)
          endif

          if(flag%output >= 3)then
             write(ch,'(i2.2)')n
             open(newunit(io),FILE='./CONTRL_OUT/spectrum_hires_'//ch//'.dat')
             write(io,'(A)')'# Wavelength / nm'
             write(io,'(A)')'# Radiance modelled'
             do k = 1, win_ini(n)%nwave_hi
                write(io,'(100(1pE16.8E3,x))') &
                     win_ini(n)%wavelength_hi(k),&
                     reflectance_hi(k, :)*win(n)%sun_spectrum_ref_hi(k)
             enddo
             close(io)
          endif

       enddo ! close loop over windows

       !*** Check if the max OT condition has anywhere been reached in the previous rad_trans calls (rad_trans_intf).
       !*** If yes, skip this retrieval:
       if(MaxOTFlag==1) then
          ExitFlag = 1
          error_id = 90
          iterflag='AOT_MAX'
          goto 100
       endif

       !*** Concatenate all retrieval windows for calling the inverse method
       call concat_windows(flag%temp, absorb, win_ini, win, ymeas, ymod, ycov, derivatives_lo, regpix)

       if(flag%output >= 3) then
          write(ch,'(i2.2)') iter
          open(newunit(io), FILE='./CONTRL_OUT/kmat_coarse_'//ch//'.out')
          do l = 1, nwave_lo
             write(io,'(100(E28.16e3,X))') ymeas(l),&
                  ycov(l),&
                  ymod(l),&
                  derivatives_lo(l,1:nstate)
          enddo
          close(io)
       endif

       !*** STOP Iteration ?
        residual = DABS(sum(x_state(1:nlay*absorb%ntype_target))/sum(x_state_old(1:nlay*absorb%ntype_target))-1.)
       !residual = abs(sum(x_state(:))/sum(x_state_old(:))-1.)
       if (iter > nlsq+1) then
          state_stop = 5.D-1*DSQRT(sum(s_state(1:nlay*absorb%ntype_target,1:nlay*absorb%ntype_target)))/sum(x_apr(1:nlay*absorb%ntype_target))
          !state_stop = sqrt(sum(s_state(:,:)))/sum(x_apr(:))
       endif

       !*** Only consider convergence if Levenberg-Marquard constraint is zero
       if(iter>nlsq+1) then
          if(residual<state_stop .and. lambda==0.D0  .and. iter>=miniter) then !nominal end of iteration
             !*** Check boundary hit for last iteration
             if(Boundary_Flag==1)then
                error_id = 91
                iterflag='AERO_BD'
             elseif(Boundary_Flag==2)then
                error_id = 99
                iterflag='XALL_BD'
             endif
             ExitFlag = 1
          elseif(iter .ge. maxiter .or. chi2>INF .or. lambda > 1.D5) then
             if(iter .ge. maxiter)then
                error_id = 92
                iterflag='ITER_MAX'
             endif
             if(lambda > 0.)then
                error_id = 93
                iterflag='ITER_LAMBDA'
             endif
             ExitFlag=1
          endif
       endif

       !*** CHI2, housekeeping Levenberg-Marquardt
       chi2 = 0.d0
       rms = 0.d0
       do l = 1, nwave_lo
          chi2 = chi2 + (ymeas(l)-ymod(l))**2/ycov_unscaled(l)
          rms = rms + (ymeas(l) - ymod(l))**2
       enddo
       !*** Get reduced chi2
       chi2 = chi2/(nwave_lo-degfreedom)
       rms = sqrt(rms/nwave_lo)
       if(chi2 > 1d10 .or. chi2 /=chi2) then
          error_id = 96 ! Chi2 too high. Abort retrieval immediately:
          iterflag='ITER_CHI2'
          exitflag = 1
       endif

       !*** Write screen/log output
       if(flag%output >=2) then
          write(message,'(X,A,7X,I2.2)') 'It#:',iter
          call writelog (message, 3)
          write(message,'(X,A,3X,1pE13.6)')'lambda:',lambda
          call writelog(message,3)
          write(message,'(X,A,1X,1(1pE13.6,X))')'CHI2_RED:',chi2
          call writelog (message, 3)
          write(message,'(X,A,6X,1(1pE13.6,X))')'RMS:',rms
          call writelog (message, 3)
          write(message,'(X,A,6X,10(1pE13.6,X))')'DFS:',degfreedom, dfs_target(:), dfs_scat
          call writelog (message, 3)
          do j = 1, absorb%ntype_target
             i1 = 1+nlay*(j-1)
             i2 = nlay*j
             write(message,'(X,A,X,I4.4,X,2(X,1pE13.6))')'abundance of target:', &
                  absorb%type_x_target(j), &
                  sum(x_state(i1:i2))/sum(dvair),DSQRT(sum(s_state(i1:i2,i1:i2)))/sum(dvair)
             call writelog (message, 3)
          enddo
          !*** state vector
          do k=1,nlay*absorb%ntype_target
             write(message,'(X,I2.2,X,A,X,I4.4,A,3(4X,1pE13.6))')k,'TAR',&
                  absorb%type_x_target(int((k-1)/nlay)+1),': ',&
                  x_state(k),&
                  sqrt(s_state(k,k)),&
                  x_state(k)/x_apr(k)
             call writelog (message, 3)
          enddo
          off=nlay*absorb%ntype_target
          do k=1,absorb%ntype_global
             write(message,FMT='(X,I2.2,X,A,X,I4.4,A,3(4X,1pE13.6))')off+k,'ABS',&
                  absorb%type_x_global(k),': ',&
                  x_state(off+k),&
                  sqrt(s_state(off+k,off+k)),&
                  x_state(off+k)/x_apr(off+k)
             call writelog (message, 3)
          enddo
          off=nlay*absorb%ntype_target+absorb%ntype_global+1
          if(flag%temp==1) then
             write(message, FMT='(X,I2.2,X,A,3(4X,1pE13.6))')off,'TEMP:     ',&
                  x_state(off), sqrt(s_state(off,off))
             off=off+1
             call writelog (message, 3)
          endif
          do n=1,nwin
             if(win_ini(n)%albflag>0) then
                do k=1,win_ini(n)%albflag
                   write(message,FMT='(X,I2.2,X,A,I2,X,A,3X,1pE13.6,4X,1pE13.6,4X,1pE13.6)')&
                        off,'Win',n,'Alb: ',&
                        x_state(off),&
                        sqrt(s_state(off,off))
                   off=off+1
                   call writelog (message, 3)
                enddo
             endif
             if(win_ini(n)%IOffFlag .ne. 0) then
                do k=1,abs(win_ini(n)%IOffFlag)
                   write(message,FMT='(X,I2.2,X,A,I2,X,A,3X,1pE13.6,4X,1pE13.6,4X,1pE13.6)')&
                        off,'Win',n,'IOff:',&
                        x_state(off),&
                        sqrt(s_state(off,off))
                   off=off+1
                   call writelog (message, 3)
                enddo
             endif
             if(win_ini(n)%Fsflag>0) then
                do k=1,win_ini(n)%Fsflag
                   write(message,FMT='(X,I2.2,X,A,I2,X,A,3X,1pE13.6,4X,1pE13.6,4X,1pE13.6)')&
                        off,'Win',n,'SIF: ',&
                        x_state(off),&
                        sqrt(s_state(off,off))
                   off=off+1
                   call writelog (message, 3)
                enddo
             endif
             if(abs(win_ini(n)%spsh0flag)==1) then
                write(message,FMT='(X,I2.2,X,A,I2,X,A,2X,1pE13.6,4X,1pE13.6)')&
                     off,'Win',n, 'SpSh0:',&
                     x_state(off),&
                     sqrt(s_state(off,off))
                off=off+1
                call writelog (message, 3)
             endif
             if(abs(win_ini(n)%spsh1flag)==1) then
                write(message,FMT='(X,I2.2,X,A,I2,X,A,2X,1pE13.6,4X,1pE13.6)')&
                     off,'Win',n, 'SpSh1:',&
                     x_state(off),&
                     sqrt(s_state(off,off))
                off=off+1
             endif
             if(win_ini(n)%spsh2flag==1) then
                write(message,FMT='(X,I2.2,X,A,I2,X,A,2X,1pE13.6,4X,1pE13.6)')&
                     off,'Win',n, 'SpSh2:',&
                     x_state(off) ,&
                     sqrt(s_state(off,off))
                off=off+1
                call writelog (message, 3)
             endif
             if(abs(win_ini(n)%sunsh0flag)==1) then
                write(message,FMT='(X,I2.2,X,A,I2,X,A,X,1pE13.6,4X,1pE13.6)')&
                     off,'Win',n, 'SunSh0:',&
                     x_state(off),&
                     sqrt(s_state(off,off))
                off=off+1
                call writelog (message, 3)
             endif
          enddo
          if(naer .gt. 0) then
             do k=1,naer
                write(message,FMT='(X,I2.2,X,A,I2,X,A,2X,2(4X,1pE13.6))')&
                     off,'AERO',k,':',&
                     x_state(off),&
                     sqrt(s_state(off,off))

                off=off+1
                call writelog (message, 3)
             enddo
          endif
       endif ! extra output

       !***If one aerosol parameter only, no LM
       if(naer<2) then
          lambda = 0.D0
       endif

       !***If retrieved aerosol OT below threshold, no aerosols
       !*** Since there is at most one cirrus parameter left this also means that lambda = 0
       if(ntype_aer == 1 .and. sum(win(:)%ot) - sum(win(:)%cot) < minaot) then
          MinAOTFlag = 1
          lambda = 0.D0
       endif
       !*** The same does not hold true if aerosol parameters > 1:
       if(ntype_aer > 1)then
          if(sum(win(:)%cot) < mincot)then
             MinCOTFlag=1
             if(naer<3)lambda=0.D0 !if only 1 aerosol parameter was fitted
          endif
          if(sum(win(:)%ot) - sum(win(:)%cot) < minaot)then
             MinAOTFlag=1
             lambda=0.D0
          endif
       endif
       !*** HH: from GOSAT
       !***If slow convergence, reduce state vector dimension
       if(lambda==0.D0 .and. iter > maxiter-5) then
          SlowConvFlag=1
       endif

       !*** LM convergence test
       if(lambda > 0.)then
          if(chi2 < chi2old*1.1) then
             if(iter > nlsq + 1) then ! contrained PT-steps
                lambda = lambda/4.
             endif
             chi2old = chi2
             x_state_old = x_state
             derivatives_lo_old = derivatives_lo
             ymod_old = ymod
          else   ! if chi2 increases with more than 1.1, than take previous iteration and redo with lambda*2.5
             x_state = x_state_old
             derivatives_lo = derivatives_lo_old
             ymod = ymod_old
             lambda = lambda*2.5
          endif
          if(lambda < 0.2) then
             lambda = 0.D0
          endif
       else
          chi2old = chi2
          x_state_old = x_state
          derivatives_lo_old = derivatives_lo
          ymod_old = ymod
       endif

       !*** Take the iteration with minimun chi^2 as solution
       if(chi2 < chi2min) then
          chi2min = chi2
          x_state_min = x_state
          ak_min = ak
          dfs_min = degfreedom
          cf_min = cf
       endif

       !*** Get iteration-dependent output
       retrieval_output%x_state(:,iter-1) = x_state(:)
       retrieval_output%chi2(iter-1) = chi2*(nwave_lo-degfreedom)  !unreduced chi2
       retrieval_output%lambda(iter) = lambda

       !*** Exit iterative solution
100    if(ExitFlag == 1) exit

       !*** Solve for state vector
       !*** Reduce vectors in call to pt_inversion if aerosol parameters have been switched off:
       !------------------------------------------------------------------
       ! Otto: We should think of doing this different. I do not like
       ! exluding elements from the state vector within the iteration.
       ! Better to regularize such that the boundaries are not hit.
       !-------------------------------------------------------------------
       allocate(red_positions(nstate), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       red_positions = 0
       aer_red = 0

       if (ntype_aer<2) then !Assume that only aerosols are to be fitted:
          if(iter <= nlsq .or. MinAOTFlag .eq. 1 .or. slowconvflag==1 )then
             !*** throw out aerosol parameters from state vector
             aer_red = naer
             red_positions(nstate-naer+1:nstate) = 1
          endif
       else
          !*** If cirrus and aerosols are to be fitted
          !*** Assume 2 retrievable aerosol parameters and 1 cirrus
          if(iter <= nlsq .or. (MinAOTFlag .eq. 1 .and. MinCOTFlag .eq. 1))then
             aer_red = naer
             red_positions(nstate-naer+1:nstate) = 1
          endif
          if(MinCOTFlag .eq. 1 .and. MinAOTFlag .eq. 0) then
             aer_red = 1
             red_positions(nstate) = 1
          endif
          if(MinAOTFlag .eq. 1 .and. MinCOTFlag .eq. 0 .and. naer == 2) then
             aer_red = 1
             red_positions(nstate-1) = 1
          endif
          if(MinAOTFlag .eq. 1 .and. MinCOTFlag .eq. 0 .and. naer == 3) then
             aer_red = 2
             red_positions(nstate-naer+1:nstate-1) = 1
          endif
       endif

       reduction = aer_red

       allocate(call_derivatives_lo(nwave_lo,nstate-reduction), &
            call_x_state(nstate-reduction), &
            call_s_state(nstate-reduction,nstate-reduction), &
            call_x_apr(nstate-reduction), &
            call_ak(nstate-reduction,nstate-reduction), &
            call_cf(nstate-reduction,nwave_lo), &
            call_upperx(nstate-reduction), &
            call_lowerx(nstate-reduction), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory allocation error'
          ierr = ierr_all
          goto 999
       endif

       if (reduction > 0) then
          reduce_i = 0
          reduce_j = 0
          do i = 1,nstate-reduction
             if(red_positions(i+reduce_i) == 1 )then
                reduce_i = reduce_i + 1
             endif
             call_derivatives_lo(:,i) = derivatives_lo(:,i+reduce_i)
             call_x_state(i) = x_state(i+reduce_i)
             call_cf(i,:) = cf(i+reduce_i,:)
             call_upperx(i) = upperx(i+reduce_i)
             call_lowerx(i) = lowerx(i+reduce_i)
             call_x_apr(i) = x_apr(i+reduce_i)
             do j = 1,nstate-reduction
                if(red_positions(j+reduce_j) == 1 )then
                   reduce_j = reduce_j + 1
                endif
                call_s_state(i,j) = s_state(i+reduce_i,j+reduce_j)
                call_ak(i,j) = ak(i+reduce_i,j+reduce_j)
             enddo
             reduce_j = 0
          enddo
       else !No reduction required, the arrays for the call to the pt_inversion subroutine are identical
          call_derivatives_lo = derivatives_lo
          call_x_state = x_state
          call_s_state = s_state
          call_ak = ak
          call_cf = cf
          call_upperx = upperx
          call_lowerx = lowerx
          call_x_apr = x_apr
       endif

       !*** Now call the actual inversion subroutine
       if(flag%inv==0) then
          call tsvd_inversion(&
               absorb%ntype_target, naer-aer_red,  &
               nwave_lo, nstate-reduction, nlay, regskill,&
               call_derivatives_lo, ymeas, ymod, ycov,&
               call_x_state, call_s_state,&
               call_x_apr,&
               call_ak, call_cf, lambda, degfreedom, dfs_target, dfs_scat, &
               call_upperx, call_lowerx, SVDFlag, Boundary_Flag)
       elseif(flag%inv==1) then
          call pt_inversion(&
               absorb%ntype_target, naer-aer_red,  &
               nwave_lo, nstate-reduction,regpix, regskill,&
               call_derivatives_lo, ymeas, ymod, ycov,&
               call_x_state, call_s_state,&
               call_x_apr,&
               call_ak, call_cf, lambda, degfreedom, dfs_target, dfs_scat, &
               call_upperx, call_lowerx, SVDFlag, Boundary_Flag)
       elseif(flag%inv==2) then
          call column_inversion(&
               absorb%ntype_target, naer-aer_red,  &
               nwave_lo, nstate-reduction,regpix, regskill,&
               call_derivatives_lo, ymeas, ymod, ycov,&
               call_x_state, call_s_state,&
               call_x_apr,&
               call_ak, call_cf, lambda, degfreedom, dfs_target, dfs_scat, &
               call_upperx, call_lowerx, SVDFlag, Boundary_Flag)
       elseif(flag%inv==3)then
          call adhoc_inversion(&
               absorb%ntype_target, naer-aer_red, &
               nwave_lo, nstate-reduction, nlay, regskill, flag%reg,&
               call_derivatives_lo, ymeas, ymod, ycov,&
               call_x_state, call_s_state,&
               call_x_apr, &
               call_ak, call_cf, lambda, degfreedom, dfs_target, dfs_scat,&
               call_upperx, call_lowerx, SVDFlag,Boundary_Flag)
       elseif(flag%inv==4 .or. flag%inv==5) then
          call adhoc_inversion_matrix(&
               absorb%ntype_target, naer-aer_red, &
               nwave_lo, nstate-reduction, nlay, regskill, flag%reg,&
               call_derivatives_lo, ymeas, ymod, ycov,&
               call_x_state, call_s_state,&
               call_x_apr, &
               call_ak, call_cf, lambda, degfreedom, dfs_target, dfs_scat,&
               call_upperx, call_lowerx, SVDFlag,Boundary_Flag, line_number,&
               flag%inv)
       endif

       ! print*, "State Vector"
       ! do l = 1, 3
       !    print*, call_x_state(l)
       ! end do

       ! print*, "State Vector Covariance Matrix"
       ! do l = 1, 3
       !    print*, call_s_state(1:4, l)
       ! end do

       ! print*, "Gain Function"
       ! do l = 1, 3
       !    print*, call_cf(1:4, l)
       ! end do

       ! print*, "Averaging Kernel"
       ! do l = 1, 3
       !    print*, call_ak(1:4, l)
       ! end do
       ! stop

       !*** Rewrite the reduced arrays used for the call to pt_inversion into the full arrays needed in profile_inversion
       if (reduction > 0) then
          reduce_i = 0
          reduce_j = 0
          do i = 1, nstate
             if(red_positions(i) == 1 ) then
                reduce_i = reduce_i + 1
             else
                derivatives_lo(:,i) = call_derivatives_lo(:,i-reduce_i)
                x_state(i) = call_x_state(i-reduce_i)
                cf(i,:) = call_cf(i-reduce_i,:)
                upperx(i) = call_upperx(i-reduce_i)
                lowerx(i) = call_lowerx(i-reduce_i)
                !                x_apr(i) = call_x_apr(i-reduce_i)
                do j = 1,nstate
                   if(red_positions(j) == 1 )then
                      reduce_j = reduce_j + 1
                   else
                      s_state(i,j) = call_s_state(i-reduce_i,j-reduce_j)
                      ak(i,j) = call_ak(i-reduce_i,j-reduce_j)
                   endif
                enddo
                reduce_j = 0
             endif
          enddo
       else !No reduction required, the arrays for the call to the pt_inversion subroutine are identical
          derivatives_lo = call_derivatives_lo
          x_state = call_x_state
          s_state = call_s_state
          ak = call_ak
          cf = call_cf
          upperx = call_upperx
          lowerx = call_lowerx
          !          x_apr = call_x_apr
       endif

       deallocate(call_derivatives_lo, &
            call_x_state, &
            call_s_state, &
            call_x_apr, &
            call_ak, &
            call_cf, &
            call_upperx, &
            call_lowerx, &
            red_positions, &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'PROFILE_INVERSION: memory deallocation error'
          ierr = ierr_deall
          goto 999
       endif

       !*** Error checks:
       if(SVDFlag<0)then
          error_id = 97
          iterflag='INV_SVD'
          exit
       endif
       if(degfreedom .ne. degfreedom)then
          error_id = 98
          iterflag='INV_NAN'
          exit
       endif
       if(Boundary_Flag == 2)then
          error_id = 99
          iterflag = 'XALL_BD'
          exit
       endif

       !*** Generate input for next iteration
       call update_forward_model( &
            flag%temp, glintflag, nlay, win_ini, &
            absorb, atm_rt, atm_xs, &
            dvair, dvair_old, play_old, tlay_old, &
            aerosol, &
            minaotflag, mincotflag, nder, win, &
            upperx, lowerx, x_apr, x_state, &
            ierr)
       if (ierr .ne. 0) return

    enddo ! close loop over iterations


    !HH: This is commented in the GOSAT code
    !      x_state = x_state_min
    !      ak = ak_min
    !      degfreedom = dfs_min
    !      cf = cf_min
    call update_forward_model( &
         flag%temp, glintflag, nlay, win_ini, &
         absorb, atm_rt, atm_xs, &
         dvair, dvair_old,  play_old, tlay_old, &
         aerosol, &
         minaotflag, mincotflag, nder, win, &
         upperx, lowerx, x_apr, x_state, &
         ierr)
    if (ierr .ne. 0) return


    !*** Convergence: 1=yes, 0=no
    convergence = 0
    if (error_ID .eq. 0 .and. iter < maxiter) then
       convergence = 1
    endif


    !****************************************************
    !*** Put relevant variables in output type
    !****************************************************

    !*** last iteration
    retrieval_output%x_state(:, iter) = x_state(:)
    retrieval_output%chi2(iter) = chi2*(nwave_lo-degfreedom) !unreduced chi2

    !*** averaging kernel
    retrieval_output%ak = ak

    !*** gain matrix
    retrieval_output%cf = cf

    !*** unscaled covariance
    s_y = 0.d0
    do i = 1, nwave_lo
       s_y(i,i) = ycov_unscaled(i)
    enddo
    retrieval_output%s_state = matmul(cf, matmul(s_y, transpose(cf)))

    !*** A priori state vector
    retrieval_output%x_apr = x_apr

    !*** air columns at retrieval layers
    l = natm/nlay
    do i = 1, nlay
       retrieval_output%dvair(i) = sum(dvair((i-1)*l+1:i*l))
    enddo

    !*** Fluorescence
    retrieval_output%Fs =  win(1)%Fs(1)

    !*** Debug info
    retrieval_output%dfs = degfreedom
    retrieval_output%dfs_target = dfs_target
    retrieval_output%dfs_scat = dfs_scat
    retrieval_output%iter = iter
    retrieval_output%iterflag = iterflag
    retrieval_output%error_id = error_id
    retrieval_output%convergence = convergence

    retrieval_output%x_state_name = x_state_name
    retrieval_output%type_x_target = absorb%type_x_target

    !*** Window dependent output
    off = 0
    do n = 1, nwin
       retrieval_output%vza(n) = win(n)%iza
       retrieval_output%raa(n) = win(n)%phi
       retrieval_output%chi2_window(n)=0.
       do l = 1, win(n)%nwave_lo
          retrieval_output%chi2_window(n) = retrieval_output%chi2_window(n) + &
               (win(n)%reflectance_meas(l) - win(n)%reflectance_lo(l))**2/ &
               ycov_unscaled(l+off)
       enddo
       off = off + win(n)%nwave_lo
       retrieval_output%ot(n) = win(n)%ot
       retrieval_output%cot(n) = win(n)%cot
       retrieval_output%albedo(n) = win(n)%albedo(1)
       retrieval_output%ny(n) = win(n)%nwave_lo
       if(win_ini(n)%spsh0flag==1) then        !wavelength shift [nm]
          retrieval_output%spectral_shift(n) = win(n)%wavelength_lo_new(1) - win(n)%wavelength_lo(1)
       endif
    enddo

    !*** Precision of retrieved albedo
    off = nlay*absorb%ntype_target + absorb%ntype_global+ flag%temp + 1
    do n = 1, nwin
       retrieval_output%albedo_err(n) = sqrt(s_state(off, off))
       off = off + win_ini(n)%albflag + abs(win_ini(n)%IOffFlag) + max(0,win_ini(n)%Fsflag) &
            + abs(win_ini(n)%spsh0flag) + abs(win_ini(n)%spsh1flag) + abs(win_ini(n)%spsh2flag) &
            + abs(win_ini(n)%sunsh0flag)
    enddo
    retrieval_output%s_apr = s_apr
    retrieval_output%air_col_old = sum(dvair_old)
    retrieval_output%rms = rms

    !** Measured and modeled spectra
    retrieval_output%spectrum_mod = ymod
    retrieval_output%spectrum_meas = ymeas
    retrieval_output%spectrum_meas_cov = ycov_unscaled
    off = 0
    do n = 1, nwin
       do k = 1, win(n)%nwave_lo
          retrieval_output%wavelength(off+k) = win(n)%wavelength_lo(k)
       enddo
       off = off + win(n)%nwave_lo
    enddo


    if(flag%output >=2 ) then
       write(message, *) '*** End of profile_inversion ***'
       call writelog(message, 1)
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine profile_inversion

  !------------------------------------------------------------------------------
  !> @brief Initialize state vector
  !> @details This subroutine puts all relevant parameters (target absorbers, interfering
  !! absorbers, aerosol, albedo, auxilary, etc) into one state vector
  ! input: nwin: number of spectral windows considered
  !        win%x_molec : atmospheric trace gas profiles
  !        type_x_target: array with target absorber indices (Hitran)
  !        ntype_target: number of target absorbers considered
  !        win%type_x: the Hitran absorber indices per window
  !        win%ntype: number of absorbers per window
  !        ntype_global: number of interfering absorbers
  !        type_x_global: array with interfering absorber indices
  !        Tflag: flag indicating temperature fit or not
  !        win%albflag: order of polynomial for albedo spectral dependence
  !        win%spsh0flag
  !        win%spsh1flag
  !        win%spsh2flag
  !        win%sunsh0flag
  !        win%IOffFlag
  !        win%Fsflag
  !        The structure aerosol
  ! Output
  !       x_state: initial state vector for retrieval
  !       upperx: array with upper boundaries of state vector elements
  !       lowerx: array with lower boundaries of state vector elements
  !       regskill: array with flags for all state vector elements indicating
  !                 how much weight they get in regularization:
  !                 > 1 : apply regularization (profile)
  !                 = 1 : apply regularization (aerosol parameters)
  !                 = 0 : do not constrain by regularization (small weight
  !                       in side constraint)
  !      s_state: ad hoc prior covariance matrix for state vector
  !               (in case of optimal estimation)
  !------------------------------------------------------------------------------
  subroutine init_state_vector(&
       Tflag, glintflag, &
       natm, nlay, absorb, &
       win_ini, win, aerosol, nder, &
       regskill, upperx, lowerx, x_state, x_state_name, s_state, &
       ierr)
    !*** Input
    integer, intent(in) :: natm, nlay, TFlag, glintflag
    type(absorbers), intent(in) :: absorb
    type(window_ini), dimension(:), intent(in) :: win_ini
    !*** Input/output
    type(window_spectrum), dimension(:), intent(inout) :: win
    type(aero), dimension(:), intent(inout) :: aerosol
    !*** Output
    integer, dimension(:), allocatable, intent(out) :: nder
    real(double), dimension(:), intent(out) :: x_state
    character*25, dimension(:), intent(out) :: x_state_name  ! State vector element identifier
    integer, dimension(size(x_state)), intent(out) :: regskill
    real(double), dimension(size(x_state)), intent(out) :: upperx
    real(double), dimension(size(x_state)), intent(out) :: lowerx
    real(double), dimension(:,:), intent(out) :: s_state
    integer, intent(out) :: ierr
    !*** Local variables
    integer :: i, j, k, l, n, off, nwin, ntype_aer, maxd
    integer, dimension(:), allocatable :: nder_dum, index              ! Aerosol height distribution
    real(double), dimension(:),allocatable :: nder_dbl                 ! Aerosol height distribution (dbl dummy)
    character(stringlen) :: message
    !------------------------------------------------------------
    !*** Initialize
    ierr = 0
    upperx = INF
    lowerx = -INF
    regskill = 0
    nwin = size(win_ini)
    l = natm/nlay
    !*** the index "off" indicates the latest position in the state vector that
    !*** has been filled after each 'set' of parameters.
    off = 0

    !*** First, put target absorber vertical profiles into the state vector
    !*** loop over all target absorbers
    do j = 1, absorb%ntype_target
       !** loop over all windows
       n = 1
       do while(n .le. nwin)
          !*** loop over all absorbers in the window till target absorber is found.
          !*** When the target absorber is put into the state vector, i is set
          !*** to win_ini(n)%ntype+1, so then jump to next target absorber (index j)
          i = 1
          do while(i .le. win_ini(n)%ntype)
             !*** check if the absorber win_ini(n)%type_x(i) is the target
             !*** absorber type_x_target(j)
             if(win_ini(n)%type_x(i)==absorb%type_x_target(j)) then
                !*** If yes, put the absorber profile into the state vector
                !*** Here, nlay sublayers are put into the state vectors
                !*** while the atmosphere has nlay*nrt*natm layers (see read_ini)
                !*** Here layer k for the state vector contains layers
                !*** (k-1)*l+1:k*l from win(n)%x_molec on the full atmosphere grid.
                !*** Remember that:   l=natm/nlay
                do k= 1, nlay
                   x_state(k+off) = sum(win(n)%x_molec((k-1)*l+1:k*l,i))
                   regskill(k+off)=1+j
                   ! regskill is a flag that determines the weight in the regularization
                   upperx(k+off) = INF
                   lowerx(k+off) = NULL
                   write(x_state_name(k+off),'(A,I2.2,A,I4.4)')'DX',k,'_TYPE',win_ini(n)%type_x(i)
                enddo
                ! The target absorber has been found, so set i and n to their
                ! maxvalues. Below 1 is added so that the while loops for i
                ! and n are completed so that then the next target absorber
                ! (index j) is considered
                i = win_ini(n)%ntype
                n = nwin
             endif
             i = i + 1
          enddo
          n = n + 1
       enddo ! end loop over windows
       ! Each target absorber 'occupies' nlay entries in the state vector
       ! so for the next target absorber, jump nlay entries further
       off = off + nlay
    enddo ! end loop over all target absorbers

    !*** Interfering absorber columns
    ! Same way as target absorber, but only the total column (no profile)
    ! is fitted
    off = nlay*absorb%ntype_target
    do j = 1, absorb%ntype_global
       n = 1
       do while(n .le. nwin)
          i = 1
          do while(i .le. win_ini(n)%ntype)
             if(win_ini(n)%type_x(i)==absorb%type_x_global(j)) then
                x_state(off+j) = sum(win(n)%x_molec(:,i)) ! whole column
                regskill(off+j)=0
                upperx(off+j) = INF
                lowerx(off+j) = NULL
                write(x_state_name(off+j),'(A,I4.4)')'X_TYPE',win_ini(n)%type_x(i)
                i = win_ini(n)%ntype
                n = nwin
             endif
             i = i + 1
          enddo
          n = n + 1
       enddo
    enddo

    !*** Temperature scaling
    off = nlay*absorb%ntype_target + absorb%ntype_global+1
    if(TFlag==1) then
       x_state(off) = 0.d0
       regskill(off) = 0
       write(x_state_name(off),'(A)')'TOFFSET'
       off=off+1
    endif

    !*** Auxiliary parameters (albedo,shifts)
    do n = 1, nwin
       if(win_ini(n)%albflag>0) then
          do k = 1, win_ini(n)%albflag
             x_state(off) = win(n)%albedo(k)
             regskill(off) = 0
             if(k==1)lowerx(off) = 0.D0
             if(k==1 .and. glintflag==1) then
                lowerx(off) = -0.08D0
             endif
             write(x_state_name(off),'(A,I2.2,A,I2.2)')'ALB_WIN',n,'_ORDER',k-1
             off = off+1
          enddo
       endif
       if(win_ini(n)%IOffFlag .ne. 0) then
          do k = 1, abs(win_ini(n)%IOffflag)
             x_state(off) = win(n)%IOff(k)
             regskill(off)=0
             write(x_state_name(off),'(A,I2.2,A,I2.2)')'IOFF_WIN',n,'_ORDER',k-1
             off=off+1
          enddo
       endif
       if(win_ini(n)%Fsflag>0) then
          do k = 1, win_ini(n)%Fsflag
             x_state(off) = win(n)%Fs(k)
             regskill(off)=0
             write(x_state_name(off),'(A,I2.2,A,I2.2)')'Fs_WIN',n,'_ORDER',k-1
             off=off+1
          enddo
       endif
       if(abs(win_ini(n)%spsh0flag)==1) then
          x_state(off) = 0.d0
          write(x_state_name(off),'(A,I2.2,A)')'SHIFT_WIN',n,'_ORDER00'
          regskill(off)=0
          off=off+1
       endif
       if(abs(win_ini(n)%spsh1flag)==1) then
          x_state(off) = 0.d0
          write(x_state_name(off),'(A,I2.2,A)')'SHIFT_WIN',n,'_ORDER01'
          regskill(off) = 0
          off=off+1
       endif
       if(abs(win_ini(n)%spsh2flag)==1) then
          x_state(off) = 0.d0
          write(x_state_name(off),'(A,I2.2,A)')'SHIFT_WIN',n,'_ORDER02'
          regskill(off)=0
          off=off+1
       endif
       if(abs(win_ini(n)%sunsh0flag)==1) then
          x_state(off) = 0.d0
          write(x_state_name(off),'(A,I2.2,A)')'SUNSHIFT_WIN',n,'_ORDER00'
          regskill(off) = 0
          off=off+1
       endif
    enddo

    !*** Aerosol parameters
    ntype_aer = size(aerosol)
    do k=1,ntype_aer
       if(aerosol(k)%AerosolFlags(1) == 1) then
          x_state(off) = aerosol(k)%reff
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'SIZE1'
          regskill(off)=1
          upperx(off) = aerosol(k)%reff*1.2
          lowerx(off) = aerosol(k)%reff/1.2
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(2) == 1) then
          x_state(off) = aerosol(k)%veff
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'SIZE2'
          regskill(off)=1
          upperx(off) = aerosol(k)%veff*1.2
          lowerx(off) = aerosol(k)%veff/1.2
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(3) == 1) then
          !DO i=1,nwin
          !x_state(nlay+off) = aerosol(1)%rm(i)
          x_state(off) = aerosol(k)%rm(1)
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'RM'
          regskill(off)=1
          upperx(off) = 1.7D0
          lowerx(off) = 1.3D0
          off=off+1
          !ENDDO
       endif
       if(aerosol(k)%AerosolFlags(4) == 1) then
          !DO i=1,nwin
          !x_state(nlay+off) = aerosol(1)%fim(i)
          x_state(off) = aerosol(k)%fim(1)
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'IM'
          regskill(off)=1
          upperx(off) = -0.0005d0
          lowerx(off) = -0.02d0
          off=off+1
          !ENDDO
       endif
       if(aerosol(k)%AerosolFlags(5) == 1) then
          ! total amount of aerosols
          if(aerosol(k)%aer_col > 0.0d0) then
             x_state(off) = aerosol(k)%aer_col
          else
             x_state(off) = NULL
          endif
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'N'
          regskill(off)=1
          upperx(off) = x_state(off)*10.d0
          lowerx(off) = x_state(off)/10.d0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(6) == 1) then
          ! fraction of spherical/plate particles for aerosol/cirrus
          x_state(off) = aerosol(k)%shapefrac
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'FRAC'
          regskill(off)=1
          upperx(off) = 1.D0
          lowerx(off) = 0.D0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(7) == 1) then
          ! central height of aerosol alt. dis.
          x_state(off) = aerosol(k)%aeralt1
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'HEIGHT1'
          regskill(off)=1
          upperx(off) = min(x_state(off) + 5.d3, 3.D4)
          lowerx(off) = max(x_state(off) - 5.d3, -1.d4)
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(8) == 1) then
          ! width of aerosol alt. dis.
          x_state(off) = aerosol(k)%aeralt2
          write(x_state_name(off),'(A,I2.2,A)')'AER',k,'HEIGHT2'
          regskill(off) = 1
          upperx(off) = 3.D4
          lowerx(off) = 1.D2
          off=off+1
       endif
    enddo

    ! (ad hoc) prior covariance matrix in case of Optimal Estimation
    s_state = 0.
    do k = 1, nlay*absorb%ntype_target
       s_state(k,k)=(.2*x_state(k))**2
    enddo

    off = nlay*absorb%ntype_target
    do k = 1, absorb%ntype_global
       s_state(off+k,off+k)=(1D4*x_state(off+k))**2
    enddo

    !*** T scaling
    off = nlay*absorb%ntype_target + absorb%ntype_global + 1
    if(TFlag==1) then
       s_state(off,off) =(1D4*x_state(off))**2
       off=off+1
    endif

    !*** Auxiliary parameters
    do n=1,nwin
       if(win_ini(n)%albflag>0) then
          do k=1,win_ini(n)%albflag
             s_state(off,off) = (1D4*win(n)%albedo(k))**2
             off=off+1
          enddo
       endif
       if(win_ini(n)%IOffFlag .ne. 0) then
          do k = 1, abs(win_ini(n)%IOffFlag)
             s_state(off,off) = (1D4)**2
             off=off+1
          enddo
       endif
       if(win_ini(n)%Fsflag>0) then
          do k = 1, win_ini(n)%Fsflag
             s_state(off,off) = (1D4)**2
             off=off+1
          enddo
       endif
       if(abs(win_ini(n)%spsh0flag)==1) then
          s_state(off,off) = (1D4)**2
          off=off+1
       endif
       if(abs(win_ini(n)%spsh1flag)==1) then
          s_state(off,off) = (1D4)**2
          off=off+1
       endif
       if(abs(win_ini(n)%spsh2flag)==1) then
          s_state(off,off) = (1D4)**2
          off=off+1
       endif
       if(abs(win_ini(n)%sunsh0flag)==1) then
          s_state(off,off) = (1D4)**2
          off=off+1
       endif
    enddo

    !*** Aerosol parameters
    do k = 1, ntype_aer
       if(aerosol(k)%AerosolFlags(1) == 1) then
          s_state(off,off)=(aerosol(k)%reff*1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(2) == 1) then
          s_state(off,off) = (aerosol(k)%veff*1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(3) == 1) then
          s_state(off,off) = (1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(4) == 1) then
          s_state(off,off) = (1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(5) == 1) then
          if (aerosol(k)%aer_col  > 0.0d0) then
             s_state(off,off) = ((aerosol(k)%aer_col) *1.D4)**2
          else
             s_state(off,off) = 0.0d0
          endif
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(6) == 1) then
          s_state(off,off) = (aerosol(k)%shapefrac*1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(7) == 1) then
          s_state(off,off) = (aerosol(k)%aeralt1*1.D4)**2
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(8) == 1) then
          s_state(off,off) = (aerosol(k)%aeralt2*1.D4)**2
          off = off + 1
       endif
    enddo

    !*** Identify layers with non-negligible aerosol
    allocate(nder_dum(sum(aerosol(:)%maxd)), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'INIT_STATE_VECTOR: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    j = 0
    do k = 1, ntype_aer
       do i = 1, aerosol(k)%maxd
          j = j + 1
          nder_dum(j) = aerosol(k)%nder(i)
       enddo
    enddo
    n = 0
    do i = 1, j
       l = nder_dum(i)
       do k = i+1, j
          if(nder_dum(k)==l .and. nder_dum(k)/=0)then
             nder_dum(k) = 0
             n = n + 1
          endif
       enddo
    enddo
    maxd = j - n

    if(allocated(nder)) deallocate(nder)
    allocate(nder(maxd), &
         nder_dbl(maxd), &
         index(maxd), &
         stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'INIT_STATE_VECTOR: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    nder = pack(nder_dum, nder_dum/=0)
    nder_dbl = dble(nder)
    call sort(nder_dbl, maxd, index)
    nder = int(nder_dbl)

    deallocate(nder_dum,  nder_dbl, index, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'INIT_STATE_VECTOR: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine init_state_vector
     !------------------------------------------------------------------------------
     !> @brief Concatenate the different windows to one array
     !------------------------------------------------------------------------------
  subroutine concat_windows(Tflag, absorb, win_ini, win, ymeas, ymod, ycov, kmod, regpix)
    !*** Input
    integer, intent(in) :: Tflag
    type(absorbers), intent(in) :: absorb
    type(window_ini), dimension(:), intent(in) :: win_ini
    type(window_spectrum), dimension(:), intent(in) :: win
    !*** Output
    real(double), dimension(:), intent(out) :: ymeas  ! Measured reflectance concatenated over all windows
    real(double), dimension(:), intent(out) :: ycov   ! Covariance of measured reflectance concatenated over all windows
    real(double), dimension(:), intent(out) :: ymod   ! Modelled reflectance concatenated over all windows
    real(double), dimension(:,:), intent(out) :: kmod ! Modelled derivatives of reflectance concatenated over all windows
    real(double), dimension(sum(win(:)%nwave_lo)), intent(out) :: regpix  ! Array of pixels to be regularized
    !*** Local variables
    integer ::i, j, k, l, n, nlay, nwin, naer, off, off1, off2
    !-----------------------------------------------------------------------------
    !*** Number of retrieval windows
    nwin = size(win_ini)
    !*** Number of atmosphere layers
    nlay =  size(win(1)%derivP_lo(1, :))

    !*** Concatenate
    off = 0
    ymeas = 0.D0
    ycov = 0.D0
    ymod = 0.D0

    do n = 1, nwin
       do k = 1, win(n)%nwave_lo
          ymeas(k+off) = win(n)%reflectance_meas(k)
          ymod(k+off) = win(n)%reflectance_lo(k)
          ycov(k+off) = win(n)%reflectance_meas_cov(k)
       enddo
       off = off + win(n)%nwave_lo
    enddo

    !*** Element 1 to nlay of state vector: CO2 vmrs
    kmod = 0.D0
    regpix = 0
    off = 0
    do j = 1, absorb%ntype_target
       off1 = 0
       do n = 1, nwin
          do i = 1, win_ini(n)%ntype
             if(win_ini(n)%type_x(i)==absorb%type_x_target(j) .and. absorb%type_x_target(j)/=7) then
                do k = 1, nlay
                   do l = 1, win(n)%nwave_lo
                      kmod(l+off1, off+k) = win(n)%derivatives_lo(l, k, i)
                      regpix(l+off1) = 1
                   enddo
                enddo
             elseif(absorb%type_x_target(j)==7) then
                do k = 1, nlay
                   do l = 1, win(n)%nwave_lo
                      kmod(l+off1,off+k) = win(n)%derivP_lo(l, k)
                      regpix(l+off1) = 1
                   enddo
                enddo
             endif
          enddo
          off1 = off1 + win(n)%nwave_lo
       enddo
       off = off + nlay
    enddo

    off = nlay*absorb%ntype_target
    do j = 1, absorb%ntype_global
       off1 = 0
       do n = 1, nwin
          do i = 1, win_ini(n)%ntype
             if(win_ini(n)%type_x(i)==absorb%type_x_global(j) .and. absorb%type_x_global(j)/=7) then
                do l=1,win(n)%nwave_lo
                   kmod(l+off1,off+j) = win(n)%derivatives_lo(l,1,i)
                   regpix(l+off1) = 1
                enddo
             elseif(absorb%type_x_global(j)==7) then
                do l = 1, win(n)%nwave_lo
                   kmod(l+off1,off+j) = win(n)%derivP_lo(l,1)
                   regpix(l+off1) = 1
                enddo
             endif
          enddo
          off1 = off1 + win(n)%nwave_lo
       enddo
    enddo

    !*** Temperature scaling
    off = nlay*absorb%ntype_target + absorb%ntype_global+1
    if(TFlag==1) then
       off1 = 0
       do n = 1, nwin
          do l = 1, win(n)%nwave_lo
             kmod(l+off1,off) = win(n)%derivT_lo(l)
          enddo
          off1 = off1 + win(n)%nwave_lo
       enddo
       off = off + 1
    endif

    !*** Element nlay+ ... of the state vector
    off2 = 0
    do n = 1, nwin
       if(win_ini(n)%albflag>0) then
          do k=1,win_ini(n)%albflag
             off1=0
             do l=1,win(n)%nwave_lo
                kmod(l+off1+off2, off) = win(n)%deriv_albedo_lo(l,k)
             enddo
             off1 = off1 + win(n)%nwave_lo
             off = off + 1
          enddo
       endif
       if(win_ini(n)%IOffFlag .ne. 0) then
          do k = 1, abs(win_ini(n)%IOffFlag)
             off1 = 0
             do l = 1, win(n)%nwave_lo
                kmod(l+off1+off2, off) = win(n)%deriv_Ioff_lo(l,k)
             enddo
             off1=off1+win(n)%nwave_lo
             off=off+1
          enddo
       endif
       if(win_ini(n)%Fsflag>0) then
          do k = 1, win_ini(n)%Fsflag
             off1 = 0
             do l = 1, win(n)%nwave_lo
                kmod(l+off1+off2, off) = win(n)%deriv_Fs_lo(l,k)
             enddo
             off1=off1+win(n)%nwave_lo
             off=off+1
          enddo
       endif
       if(abs(win_ini(n)%spsh0flag)==1) then
          off1 = 0
          do l = 1, win(n)%nwave_lo
             kmod(l+off1+off2,off)=win(n)%deriv_specshift0_lo(l)
          enddo
          off1 = off1+win(n)%nwave_lo
          off = off+1
       endif
       if(abs(win_ini(n)%spsh1flag)==1) then
          off1=0
          do l=1,win(n)%nwave_lo
             kmod(l+off1+off2,off)=win(n)%deriv_specshift1_lo(l)
          enddo
          off1 = off1+win(n)%nwave_lo
          off = off + 1
       endif
       if(abs(win_ini(n)%spsh2flag)==1) then
          off1=0
          do l=1,win(n)%nwave_lo
             kmod(l+off1+off2,off)=win(n)%deriv_specshift2_lo(l)
          enddo
          off1=off1+win(n)%nwave_lo
          off=off+1
       endif
       if(abs(win_ini(n)%sunsh0flag)==1) then
          off1 = 0
          do l = 1, win(n)%nwave_lo
             kmod(l+off1+off2,off) = win(n)%deriv_sunshift0_lo(l)
          enddo
          off1=off1+win(n)%nwave_lo
          off = off + 1
       endif
       off2 = off2 + win(n)%nwave_lo
    enddo !loop over windows

    !*** Aerosol parameters
    naer = size(win(1)%deriv_aerosol_lo(1,:))
    do k = 1, naer
       off1 = 0
       do n = 1, nwin
          do l = 1, win(n)%nwave_lo
             kmod(l+off1, off) = win(n)%deriv_aerosol_lo(l,k)
             regpix(l+off1) = 1
          enddo
          off1 = off1 + win(n)%nwave_lo
       enddo
       off = off + 1
    enddo

  end subroutine concat_windows
  !------------------------------------------------------------------------------
  !> @brief Update variables and state vector boundaries for next iteration
  !------------------------------------------------------------------------------
  subroutine update_forward_model( &
       Tflag, glintflag, nlay, &
       win_ini, absorb, &
       atm_rt, atm_xs, &
       dvair, dvair_old, play_old, tlay_old, &
       aerosol, minaotflag, mincotflag, nder, &
       win, upperx, lowerx, x_apr, x_state, &
       ierr)
    !*** Input
    integer, intent(in) :: Tflag, glintflag, nlay
    type(window_ini), dimension(:), intent(in) :: win_ini
    type(absorbers) :: absorb
    type(atmosphere), intent(in) :: atm_rt
    real(double), dimension(:), intent(in) :: x_apr
    real(double), dimension(:), intent(inout) :: dvair      ! Partial air column, subject to change in O2 retrieval (Dim: natm)
    real(double), dimension(:), intent(in) :: dvair_old  ! Partial air column (Dim: natm)
    real(double), dimension(:), intent(in):: play_old
    real(double), dimension(:), intent(in):: tlay_old
    real(double), dimension(:), intent(inout) :: x_state
    integer, intent(in) :: minaotflag, mincotflag
    !*** Input/output
    type(atmosphere), intent(inout) :: atm_xs
    type(window_spectrum), dimension(:), intent(inout) :: win
    type(aero), dimension(:), intent(inout) :: aerosol
    integer, dimension(:), allocatable, intent(out) :: nder
    !*** Output
    real(double), dimension(size(x_state)), intent(out) :: upperx
    real(double), dimension(size(x_state)), intent(out) :: lowerx
    integer, intent(out) :: ierr
    !*** Local variables
    integer :: i, j, k, l, m, n, off, natm, nwin, ntype_aer, maxd
    integer, dimension(:), allocatable :: nder_dum, index                ! Aerosol height distribution
    real(double), dimension(:), allocatable :: nder_dbl                 ! Aerosol height distribution (dbl dummy)
    real(double) :: wavelength_lo_new
    character(stringlen) :: message

    !*** Initialize
    ierr = 0
    natm = atm_xs%n
    upperx = INF
    lowerx = -INF
    nwin = size(win_ini)
    !*** Target absorber vertical profiles
    l = natm/nlay
    off = 0
    do j = 1, absorb%ntype_target
       do n = 1, nwin
          do i = 1, win_ini(n)%ntype
             if(win_ini(n)%type_x(i)==absorb%type_x_target(j)) then
                do k = 0, natm-1
                   m = int(k/l) + 1
                   win(n)%x_molec(k+1,i)=x_state(m+off)/x_apr(m+off)*win(n)%dv_x(k+1,i)
                   upperx(m+off) = INF
                   lowerx(m+off) = NULL
                   if(absorb%type_x_target(j)==7) then
                      dvair(k+1)=x_state(m+off)/x_apr(m+off)*dvair_old(k+1)
                      atm_xs%p(k+1)=x_state(m+off)/x_apr(m+off)*play_old(k+1)
                   endif
                enddo
             endif
          enddo
       enddo
       off = off + nlay
    enddo

    !*** Interfering absorber columns
    off = nlay*absorb%ntype_target
    do j = 1, absorb%ntype_global
       do n = 1, nwin
          do i = 1, win_ini(n)%ntype
             if(win_ini(n)%type_x(i)==absorb%type_x_global(j)) then
                do k = 1, natm
                   win(n)%x_molec(k,i)=x_state(off+j)/x_apr(off+j)*win(n)%dv_x(k,i)
                enddo
                if(absorb%type_x_global(j)==7) then
                   do k = 1, natm
                      dvair(k)=x_state(off+j)/x_apr(off+j)*dvair_old(k)
                      atm_xs%p(k) = x_state(off+j)/x_apr(off+j)*play_old(k)
                   enddo
                endif
                upperx(off+j) = INF
                lowerx(off+j) = NULL
             endif
          enddo
       enddo
    enddo

    !*** Temperature scaling
    off = nlay*absorb%ntype_target + absorb%ntype_global + 1
    if(TFlag==1) then
       atm_xs%t = x_state(off) + tlay_old
       off = off + 1
    endif

    do n = 1, nwin

       !*** Albedo, shift + stretch/squeeze
       if(win_ini(n)%albflag>0) then
          do k=1,win_ini(n)%albflag
             win(n)%albedo(k)=x_state(off)
             if(k==1)lowerx(off) = 0.D0
             if(k==1 .and. glintflag==1) then
                lowerx(off) = -0.08D0
             endif
             off=off+1
          enddo
       endif

       !*** Intensity offset
       if(win_ini(n)%IOffFlag .ne. 0) then
          do k = 1, abs(win_ini(n)%IOffFlag)
             win(n)%IOff(k)=x_state(off)
             off=off+1
          enddo
       endif

       !*** Fluorescence emission
       if(win_ini(n)%Fsflag>0) then
          do k = 1, win_ini(n)%Fsflag
             win(n)%Fs(k)=x_state(off)
             off=off+1
          enddo
       endif

       !*** Spectral shift
       do k = 1, win(n)%nwave_lo
          win(n)%wavelength_lo_new(k) = win(n)%wavelength_lo(k)
       enddo
       if(win_ini(n)%spsh0flag==1) then       !wavelength shift [nm]
          do k=1,win(n)%nwave_lo
             win(n)%wavelength_lo_new(k) = win(n)%wavelength_lo_new(k) + &
                  x_state(off)
          enddo
          off = off + 1
       endif
       if(win_ini(n)%spsh1flag==1) then
          do k=1,win(n)%nwave_lo
             win(n)%wavelength_lo_new(k) = win(n)%wavelength_lo_new(k) + &
                  x_state(off)*(win(n)%wavelength_lo(k)-win(n)%wavelength_lo(1))
          enddo
          off = off + 1
       endif
       if (win_ini(n)%spsh2flag==1)then
          do k=1,win(n)%nwave_lo
             win(n)%wavelength_lo_new(k) = win(n)%wavelength_lo_new(k) + &
                  x_state(off)*(win(n)%wavelength_lo(k)-win(n)%wavelength_lo(1))**2
          enddo
          off = off + 1
       endif
       if(win_ini(n)%sunsh0flag==1) then       !wavelength shift
          do k=1,win_ini(n)%nwave_hi
             win(n)%wavelength_hi_new(k) = win_ini(n)%wavelength_hi(k) + x_state(off)
          enddo
          off = off + 1
          !*** Interpolate shifted solar spectrum back to high-resolution model wavelength grid
          call spline_interpol(win(n)%wavelength_hi_new, win(n)%sun_spectrum_sat_hi, win_ini(n)%nwave_hi,&
               win_ini(n)%wavelength_hi, win(n)%sun_spectrum_ref_hi, win_ini(n)%nwave_hi, ierr)
          if (ierr .ne. 0) then
             write(message,*) 'PROFILE_INVERSION.SPLINE_INTERPOL.SPLINT: bad input'
             ierr = ierr_intrpl
             goto 999
          endif
          !*** Convolve shifted solar spectrum by instrument response function:
          ! moved to forward_model_lo
          !            call spectral_response_stored( &
          !                 win(n)%resp_store, &
          !                 win(n)%ie_store, &
          !                 win(n)%is_store, &
          !                 win(n)%sun_spectrum_ref_hi, &
          !                 win(n)%sun_spectrum_ref_lo)
       endif
    enddo ! loop over nwin

    !*** Aerosol parameters
    ntype_aer = size(aerosol)
    do k = 1, ntype_aer
       if(aerosol(k)%AerosolFlags(1)==1) then
          aerosol(k)%reff = x_state(off)
          upperx(off) = aerosol(k)%reff*1.2d0
          lowerx(off) = aerosol(k)%reff/1.2d0
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(2)==1) then
          aerosol(k)%veff = x_state(off)
          upperx(off) = aerosol(k)%veff*1.2d0
          lowerx(off) = aerosol(k)%veff/1.2d0
          off = off + 1
       endif
       if(aerosol(k)%AerosolFlags(3)==1) then
          aerosol(k)%rm(:) = x_state(off)
          upperx(off) = 1.7D0
          lowerx(off) = 1.3D0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(4)==1) then
          aerosol(k)%fim(:)=x_state(off)
          upperx(off) = -0.0005d0
          lowerx(off) = -0.02d0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(5)==1 .and. aerosol(k)%CirrusFlag .ne. 1) then
          if(MinAOTFlag .eq. 1) x_state(off) = 0.d0
          aerosol(k)%aer_col = x_state(off)
          upperx(off) = x_state(off)*10.d0
          lowerx(off) = x_state(off)/10.d0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(5)==1 .and. aerosol(k)%CirrusFlag .eq. 1) then
          if(MinCOTFlag .eq. 1) x_state(off) = 0.d0
          aerosol(k)%aer_col =  x_state(off)
          upperx(off) = x_state(off)*10.d0
          lowerx(off) = x_state(off)/10.d0
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(6)==1) then
          aerosol(k)%shapefrac = x_state(off)
          upperx(off) = 1.d0
          lowerx(off) = 0.d0
          off=off+1
       endif

       if(aerosol(k)%AerosolFlags(7)==1) then
          aerosol(k)%aeralt1 = x_state(off)
          upperx(off) = min(x_state(off) + 5.d3, 3.D4)
          lowerx(off) = max(x_state(off) - 5.d3, -1.d4)
          off=off+1
       endif
       if(aerosol(k)%AerosolFlags(8)==1) then
          aerosol(k)%aeralt2 = x_state(off)
          upperx(off) = 3.d4
          lowerx(off) = 1.d2
          off=off+1
       endif
       call set_altdis(atm_rt, &
            aerosol(k)%altid, &
            aerosol(k)%aeralt1, &
            aerosol(k)%aeralt2, &
            aerosol(k)%alt_dis, &
            aerosol(k)%dalt_daer1, &
            aerosol(k)%dalt_daer2)

       !*** Reset nder, array of layer indices with significant aerosol contribution
       aerosol(k)%maxd=count(aerosol(k)%alt_dis(:)>maxval(aerosol(k)%alt_dis)*nder_cut)

       if(allocated(aerosol(k)%nder))deallocate(aerosol(k)%nder)
       allocate(aerosol(k)%nder(aerosol(k)%maxd), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'UPDATE_FORWARD_MODEL: memory allocation error'
          ierr = ierr_all
          goto 999
       endif

       j = 0
       do i = 1, size(aerosol(k)%alt_dis)
          if(aerosol(k)%alt_dis(i)>maxval(aerosol(k)%alt_dis)*nder_cut) then
             j = j + 1
             aerosol(k)%nder(j) = i
          else
             aerosol(k)%alt_dis(i) = 0.D0
             aerosol(k)%dalt_daer1(i) = 0.D0
             aerosol(k)%dalt_daer2(i) = 0.D0
          endif
       enddo
    enddo !k = 1, ntype_aer

    allocate(nder_dum(sum(aerosol(:)%maxd)), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'UPDATE_FORWARD_MODEL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    j = 0
    do k = 1, ntype_aer
       do i = 1, aerosol(k)%maxd
          j = j + 1
          nder_dum(j) = aerosol(k)%nder(i)
       enddo
    enddo

    n = 0
    do i = 1, j
       l = nder_dum(i)
       do k = i+1, j
          if(nder_dum(k)==l .and. nder_dum(k)/=0)then
             nder_dum(k) = 0
             n = n + 1
          endif
       enddo
    enddo
    maxd = j - n

    if(allocated(nder)) deallocate(nder)
    allocate(nder(maxd), &
         nder_dbl(maxd), &
         index(maxd), &
         stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'UPDATE_FORWARD_MODEL: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    nder = pack(nder_dum, nder_dum/=0)
    nder_dbl = dble(nder)
    call sort(nder_dbl, maxd, index)
    nder = int(nder_dbl)

    deallocate(nder_dum, nder_dbl, index, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'UPDATE_FORWARD_MODEL: memory deallocation error'
       ierr = ierr_deall
       goto 999
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine update_forward_model

  !------------------------------------------------------------------------------
end module profile_inversion_module
