!------------------------------------------------------------------------------
!> @brief compute optical properties of aerosol and cirrus
!> @todo  Possibility to speed further up: in OpticM_module_28_10 use separate LUT
!! for coefficients instead of calculating them from the phase matrix
!------------------------------------------------------------------------------
module optic_input_module
  use header_module
  use aerosol_input_module, only: aero, set_altdis, get_aerosol_properties_lognormal
  use optic_cirrus_module, only: cirrus_table, read_cirrus_netcdf, optic_cirrus, &
       der_par_mode, optic_cirrus_ssc, der_par_mode_ssc
  use OpticM_module, only: mie_lut, read_aerosol_netcdf, modes_calc_ssc, modes_calc
  implicit none
  private

  !*** types
  public :: Mie_lut, cirrus_table, aero, aero_opt

  !*** procedures
  public :: read_aerosol_netcdf,  read_cirrus_netcdf, set_altdis, get_aerosol_properties_lognormal
  public :: optic_all_der, optic_all_msc, optic_all_ssc
  private :: optic_phase_ssc

  !*** Indices to convert scattering matrix to array  
  integer, dimension(4) :: index_ist = (/1,2,3,1 /)
  integer, dimension(4) :: index_jst = (/1,2,3,2 /)

  !------------------------------------------------------------------------------
  !> @brief Optical properties of aerosol
  type :: aero_opt
     !> Aerosol absorption optical thickness (dim: nwave, nrt)             
     real(double), dimension(:,:), allocatable :: taua_aer
     !> Aerosol scattering optical thickness (dim: nwave, nrt)
     real(double), dimension(:,:), allocatable :: taus_aer
     !> Single scattering phase matrix for aerosol scattering (dim: nstokes, nstokes, nrt)   
     real(double), dimension(:, :, :), allocatable :: zss_aer
     !> Derivative of individual aerosol-type scattering optical thickness wrt Mie aerosol parameters 
     real(double), dimension(npar_mie) :: dtaus_aer
     !> Derivative of individual aerosol-type absorption optical thickness wrt Mie aerosol parameters 
     real(double), dimension(npar_mie) :: dtaua_aer
     !> Derivatives of single scattering phase matrix wrt Mie aerosol parameters
     real(double), dimension(nstokes, nstokes, npar_mie) :: d_zss_aer  
     !> Derivative of combined aerosol scattering optical thickness wrt Mie aerosol parameters (dim: nrt, npar_mie)  
     real(double), dimension(:, :), allocatable :: dtaus_tot_aer  
     !> Derivative of combined aerosol absorption optical thickness wrt Mie aerosol parameters (dim: nrt, npar_mie)
     real(double), dimension(:, :), allocatable :: dtaua_tot_aer  
     !> Derivatives of coefficients of individual aerosol-type scattering matrix 
     !! wrt Mie aerosol parameters
     real(double), dimension(nstokes, nstokes, 0:maxleg, npar_mie) :: dphase_aer 
     !> Derivative of total scattering optical thickness wrt all aerosol parameters (dim: maxd, npar)
     real(double), dimension(:,:), allocatable :: dtaus_all
     !> Derivative of total absorption optical thickness wrt all aerosol parameters (dim: maxd, npar)      
     real(double), dimension(:,:), allocatable :: dtaua_all
     !> Derivatives of coefficients of combined aerosol scattering matrix
     !!  wrt Mie aerosol parameters (dim: nstokes, nstokes, maxd, npar_mie) 
     real(double), dimension(:,:,:,:), allocatable :: d_zss_tot_aer
     !> Derivatives of coefficients of combined aerosol scattering matrix 
     !! wrt Mie aerosol parameters (dim: nper, 0:maxleg, maxd, npar_mie)
     real(double), dimension(:,:,:,:), allocatable :: dphase_tot_aer 
     !> Derivatives of coefficients of total scattering matrix 
     !! wrt aerosol parameters (dim: nper, 0:maxleg, maxd, npar))     	  
     real(double), dimension(:,:,:,:), allocatable :: dphase_all
     !> Derivatives of coefficients of total single scattering Z matrix 
     !! wrt aerosol parameters (dim: nstokes, maxd, npar)
     real(double), dimension(:,:,:), allocatable :: dzss_all
!> Derivatives of coefficients of total single scattering Z matrix 
!! wrt aerosol parameters (dim: nstokes, maxd, npar) - at band begin and end (for linear interpolation in rad_trans_intf.f90), 
      real(double), dimension(:,:,:), allocatable :: dzss_all_1	
      real(double), dimension(:,:,:), allocatable :: dzss_all_mid	      
      real(double), dimension(:,:,:), allocatable :: dzss_all_nwave	
     	       
     !> first order specttral correction factors for Mie-derivatives Right now, these ase used in rad_tans_intf.f90
     real(double), dimension( npar_mie) :: dcabs_dwave 
     real(double), dimension( npar_mie) :: dcsca_dwave 
  end type aero_opt

  !------------------------------------------------------------------------------

contains
  !------------------------------------------------------------------------------  
  !> @details Calculate aerosol scattering and absorption optical thickness for 
  !! different altitude layers as a function of wavelength. 
  !! Calculate phase function for multiple scattering at middle of band.
  !------------------------------------------------------------------------------  
  subroutine optic_all_msc( &
       aero_lut, cirrus_lut, &
       iband, wavelength, &
       aerosol, nder, nrt, aerosol_optic, &
       taua_aer, taus_aer, ot, cot, &
       plmom_aer, maxcoefs, ierr)
    !*** Input
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    integer, intent(in) :: iband, nrt
    real(double), dimension(:), intent(in) :: wavelength
    type(aero), dimension(:), intent(in) :: aerosol
    integer, dimension(:), intent(in) :: nder
!!$   real(double), dimension(:,:), intent(in) :: taus_aer
    !*** Input/output
    type(aero_opt), dimension(:), intent(inout) :: aerosol_optic  
    !*** Output
    real(double), dimension(:,:),  intent(out) :: taua_aer  ! Aerosol absorption optical thickness (nrt)
    real(double), dimension(:,:),  intent(out) :: taus_aer  ! Aerosol scattering optical thickness (nrt)
    real(double), intent(out) :: ot   ! aerosol+cirrus optical thickness 
    real(double), intent(out) :: cot  ! Cirrus optical thickness    
    real(double), dimension(nstokes,nstokes,0:maxleg,nrt), intent(out) :: plmom_aer  ! Phase matrix coefficients
    integer, dimension(:), intent(out) :: maxcoefs         
    integer, intent(out) :: ierr
    !*** Local variables
    ! Use cirrus extinction and scattering cross-sections with delta approximation
    logical, parameter :: Trunc_Flag = .true. 
    integer :: nangle, i, j, k, l, n, ik, i_st, j_st, iper, ntype_aer
    real(double) :: rlambda, csca, cabs, f
    ! plmom refer both to the expansion coefficients
    real(double), dimension(nstokes, nstokes, 0:maxleg) :: plmom_aer_lay    
    real(double), dimension(nfull) :: dcsca_full, dcabs_full
    real(double), dimension(nstokes, nstokes, 0:maxleg, nfull) :: dcoefs_full
    real(double), dimension(dim_x) :: rnmb
    real(double), dimension(dim_x, 6) :: der_rnmb     
    real(double), dimension(npar_mie) :: aerosol_pars_mode
    integer :: maxd
    real(double), dimension(npar_mie) :: dcsca_aer, dcabs_aer    
    real(double), dimension(nrt, size(aerosol)) :: maxcoefs_aer  
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt, size(aerosol)) :: plmom_aer_local
    integer :: nwave, iwave, iw, imid
    integer, dimension(3) :: index_wave 
    real(double), dimension(3, npar_mie) :: dcsca_aer_ref	
    real(double), dimension(3, npar_mie) :: dcabs_aer_ref
    character(stringlen) :: message
    !----------------------------------------------------------------------------------------     

    !*** Initialize
    ierr = 0
    ntype_aer = size(aerosol)
    maxd = size(nder)
    nwave = size(wavelength)
    imid = (nwave+1)/2
    index_wave = (/1, nwave, imid/)
    ot = 0.d0
    cot = 0.d0
    plmom_aer_local = 0.d0

    do i = 1, ntype_aer
       if(allocated(aerosol_optic(i)%taus_aer))deallocate(aerosol_optic(i)%taus_aer)
       if(allocated(aerosol_optic(i)%taua_aer))deallocate(aerosol_optic(i)%taua_aer)
       allocate(&
            aerosol_optic(i)%taus_aer(nwave,nrt),&
            aerosol_optic(i)%taua_aer(nwave,nrt), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_ALL_MSC: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       do k = 1, nrt
          do j = 1, nwave
             aerosol_optic(i)%taus_aer(j,k)=0.D0
             aerosol_optic(i)%taua_aer(j,k)=0.D0
          enddo
       enddo
       !*** Calculate optical properties of individual aerosol types
       do iw = 1,3
          iwave = index_wave(iw)
          rlambda = wavelength(iwave) * 1.d-3 ! wavelength [um]
          if(aerosol(i)%CirrusFlag==1) then   !cirrus           
             call optic_cirrus( &
                  cirrus_lut, &
                  Trunc_Flag, &
                  aerosol(i)%tilt_angle, &
                  aerosol(i)%shapefrac, &
                  aerosol(i)%reff, &
                  rlambda, & 
                  csca, &
                  cabs, &
                  f, &
                  plmom_aer_lay, &
                  dcsca_full, &
                  dcabs_full, &
                  dcoefs_full, &
                  nangle, &
                  rnmb, &
                  der_rnmb) 
             call der_par_mode(&
                  dcoefs_full,&
                  dcsca_full,&
                  dcabs_full,&
                  nangle,&
                  der_rnmb,&
                  aerosol(i)%aer_col,&    
                  csca,&
                  cabs,&
                  aerosol_optic(i)%dphase_aer,&
                  dcsca_aer, &
                  dcabs_aer)
          elseif(aerosol(i)%CirrusFlag .ne. 1 )then ! aerosol	
             aerosol_pars_mode(1) = aerosol(i)%reff
             aerosol_pars_mode(2) = aerosol(i)%veff
             aerosol_pars_mode(3) = aerosol(i)%rm(2*iband-1) + (aerosol(i)%rm(2*iband)-aerosol(i)%rm(2*iband-1))/(nwave-1) *(iwave-1)
             aerosol_pars_mode(4) = aerosol(i)%fim(2*iband-1) + (aerosol(i)%fim(2*iband)-aerosol(i)%fim(2*iband-1))/(nwave-1) *(iwave-1)
             aerosol_pars_mode(5) = aerosol(i)%aer_col
             aerosol_pars_mode(6) = aerosol(i)%shapefrac !fraction of spheres           
             call modes_calc(&
                  aero_lut, &
                  aerosol(i)%id,&
                  aerosol_pars_mode,&
                  rlambda,&
                  csca, &
                  cabs, &
                  dcsca_aer,&
                  dcabs_aer,&
                  plmom_aer_lay,&
                  aerosol_optic(i)%dphase_aer, &
                  ierr)
             if (ierr .ne. 0) return      
             nangle = maxleg
          endif ! aerosol or cirrus

          forall(i_st=1:nstokes, j_st=1:nstokes, k=0:maxleg, j=1:nrt)
             plmom_aer_local(i_st, j_st, k, j, i) = plmom_aer_lay(i_st, j_st, k)
          end forall

          do j = 1, maxd
             k = nder(j)
             aerosol_optic(i)%taus_aer(iwave,k) = &
                  (aerosol(i)%aer_col)*csca*aerosol(i)%alt_dis(k)               
             aerosol_optic(i)%taua_aer(iwave,k) = &
                  (aerosol(i)%aer_col)*cabs*aerosol(i)%alt_dis(k)               
          enddo

          do k = 1, npar_mie
             dcsca_aer_ref(iw, k) = dcsca_aer(k)
             dcabs_aer_ref(iw, k) = dcabs_aer(k)
          enddo

          !*** Reference calculation in middle of band 
          !*** Linear spectral dependence to be added (through interpolation) in rad_trans_intf.f90
          if( iwave .eq. imid) then 
             !*** Calculate COT and AOT at middle of the band (assuming it's the average of the edges)
             if(aerosol(i)%CirrusFlag == 1) then   !cirrus              
                cot = cot + sum(aerosol_optic(i)%taus_aer(iwave,:))/(1.-f) + sum(aerosol_optic(i)%taua_aer(iwave,:))
             else
                ot = ot + sum(aerosol_optic(i)%taus_aer(iwave,:)) + sum(aerosol_optic(i)%taua_aer(iwave,:))
             endif
             aerosol_optic(i)%dtaus_aer = dcsca_aer * (aerosol(i)%aer_col) 
             aerosol_optic(i)%dtaua_aer = dcabs_aer * (aerosol(i)%aer_col) 
             maxcoefs_aer(:, i) = 2 ! rayleigh
             forall(k=1:maxd) maxcoefs_aer(nder(k),i) = min(nangle,maxleg)          
             maxcoefs = maxcoefs_aer(:,1)
             if(ntype_aer > 1 ) then
                do j = 1, maxd
                   k = nder(j)
                   if(maxcoefs_aer(k,i) > maxcoefs(k)) maxcoefs(k) = maxcoefs_aer(k, i)
                enddo
             endif
          endif
       enddo !loop over three spectral calc. points (iw)

       !*** Needed for linearly correcting spectral dependence of derivatives towards scattering and absorption cross sections
       do k = 1, npar_mie ! Only Mie-parameters are linearly interpolated to correct for spectral dependence of dsigma/dMIEparam 
          if (dcabs_aer_ref(3,k) .ne. 0.d0) then
             aerosol_optic(i)%dcabs_dwave(k) = (dcabs_aer_ref(2,k) - dcabs_aer_ref(1,k) )/ &
                  (wavelength(nwave) - wavelength(1))/dcabs_aer_ref(3,k)
          else
             aerosol_optic(i)%dcabs_dwave(k) = 0.d0
          endif
          if( dcsca_aer_ref(3,k) .ne. 0.d0)then
             aerosol_optic(i)%dcsca_dwave(k) = (dcsca_aer_ref(2,k) - dcsca_aer_ref(1,k) )/ &
                  ( wavelength(nwave) - wavelength(1))/dcsca_aer_ref(3,k)
          else
             aerosol_optic(i)%dcsca_dwave(k) = 0.d0
          endif
       enddo

    enddo ! end loop over aerosoltype

    !*** linearly interpolate within band
    do i = 1, ntype_aer
       do j = 1, maxd
          k = nder(j)
          do l = 2, imid-1
             aerosol_optic(i)%taus_aer(l,k) = aerosol_optic(i)%taus_aer(1,k) + &
                  (aerosol_optic(i)%taus_aer(imid,k)-aerosol_optic(i)%taus_aer(1,k))&
                  /(wavelength(imid)-wavelength(1)) * (wavelength(l)-wavelength(1))
             aerosol_optic(i)%taua_aer(l,k) = aerosol_optic(i)%taua_aer(1,k) + &
                  (aerosol_optic(i)%taua_aer(imid,k)-aerosol_optic(i)%taua_aer(1,k))&
                  /(wavelength(imid)-wavelength(1)) * (wavelength(l)-wavelength(1))
          enddo
          do l = imid+1, nwave-1
             aerosol_optic(i)%taus_aer(l,k) = aerosol_optic(i)%taus_aer(imid,k) + &
                  (aerosol_optic(i)%taus_aer(nwave,k)-aerosol_optic(i)%taus_aer(imid,k))&
                  /(wavelength(nwave)-wavelength(imid)) * (wavelength(l)-wavelength(imid))
             aerosol_optic(i)%taua_aer(l,k) = aerosol_optic(i)%taua_aer(imid,k) + &
                  (aerosol_optic(i)%taua_aer(nwave,k)-aerosol_optic(i)%taua_aer(imid,k))&
                  /(wavelength(nwave)-wavelength(imid)) * (wavelength(l)-wavelength(imid))            
          enddo
       enddo
    enddo

    if(ntype_aer == 1) then 
       do k = 1, nrt
          do l = 1, nwave  
             taus_aer(l, k) = aerosol_optic(1)%taus_aer(l, k)
             taua_aer(l, k) = aerosol_optic(1)%taua_aer(l, k)
          enddo
       enddo
    elseif(ntype_aer > 1) then
       do k = 1, nrt     
          do l = 1, nwave
             taus_aer(l, k) = 0.D0
             taua_aer(l, k) = 0.D0
          enddo
       enddo
       do i = 1, ntype_aer
          do j = 1, maxd
             k = nder(j)
             do l = 1, nwave
                taus_aer(l, k) = taus_aer(l, k) + aerosol_optic(i)%taus_aer(l, k)
                taua_aer(l, k) = taua_aer(l, k) + aerosol_optic(i)%taua_aer(l, k)
             enddo
          enddo
       enddo
    endif

    !*** Calculate derivatives of individual aerosol-type properties wrt. composite aerosol properties
    plmom_aer = 0.d0      
    do i = 1, ntype_aer
       if(allocated(aerosol_optic(i)%dphase_tot_aer)) deallocate(aerosol_optic(i)%dphase_tot_aer)
       if(allocated(aerosol_optic(i)%dtaus_tot_aer)) deallocate(aerosol_optic(i)%dtaus_tot_aer)
       allocate(aerosol_optic(i)%dphase_tot_aer(nper, 0:maxleg, maxd, npar_mie), &
            aerosol_optic(i)%dtaus_tot_aer(nrt, npar_mie), &
            aerosol_optic(i)%dtaua_tot_aer(nrt, npar_mie), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_ALL_MSC: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       aerosol_optic(i)%dtaus_tot_aer = 0.D0
       aerosol_optic(i)%dtaua_tot_aer = 0.D0
       aerosol_optic(i)%dphase_tot_aer = 0.D0
       do l = 1, npar_mie
          do k = 1, maxd
             ik = nder(k)
             aerosol_optic(i)%dtaus_tot_aer(ik,l) = aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik)
             aerosol_optic(i)%dtaua_tot_aer(ik,l) = aerosol_optic(i)%dtaua_aer(l)*aerosol(i)%alt_dis(ik)
          enddo
       enddo

       do iper = 1, nper
          i_st = index_ist(iper)
          j_st = index_jst(iper)
          if (i_st .le. nstokes .and. j_st .le. nstokes) then
             do k = 1, maxd
                ik = nder(k)
                do l = 1, npar_mie
                   do n = 0, maxcoefs(ik)
                      if(taus_aer(iwave,ik)/=0.D0)then
                         aerosol_optic(i)%dphase_tot_aer(iper, n, k, l) = &
                              (aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik)* &
                              plmom_aer_local(i_st, j_st, n, ik, i) + &
                              aerosol_optic(i)%taus_aer(iwave, ik)*aerosol_optic(i)%dphase_aer(i_st, j_st, n, l))/taus_aer(iwave, ik)
                         !*** The derivatives of the phase function coefficients of a given mode depend on the other modes
                         do j = 1, ntype_aer
                            aerosol_optic(i)%dphase_tot_aer(iper, n, k, l) = &
                                 aerosol_optic(i)%dphase_tot_aer(iper, n, k, l) - &
                                 (aerosol_optic(j)%taus_aer(iwave, ik)*plmom_aer_local(i_st, j_st, n, ik, j)* &
                                 aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik))/taus_aer(iwave,ik)/taus_aer(iwave, ik)
                         enddo
                      else
                         aerosol_optic(i)%dphase_tot_aer(iper, n, k, l) = 0.D0
                      endif
                   enddo ! n
                enddo ! l
             enddo  ! k 
          endif
       enddo ! iper

       !*** Calculate phase matrix coefficients for single aerosol scattering (plmom_aer)
       !*** Combine individual aerosol type properties to composite aerosol properties by weighting with optical thickness 
       do j_st = 1, nstokes
          do i_st = 1, nstokes
             do k = 1, maxd
                ik = nder(k)
                do n = 0, maxcoefs(ik)
                   if(taus_aer(iwave, ik) /= 0.D0) then
                      plmom_aer(i_st, j_st, n, ik) = plmom_aer(i_st, j_st, n, ik) + & 
                           (aerosol_optic(i)%taus_aer(iwave,ik)*plmom_aer_local(i_st, j_st, n, ik, i)) &
                           /taus_aer(iwave,ik)
                   else
                      plmom_aer(i_st, j_st, n, ik) = 0.D0
                   endif
                enddo
             enddo
          enddo
       enddo

    enddo ! end loop over aerosol types

    return

999 continue
    call stopretrieval(message) 

  end subroutine optic_all_msc
  !------------------------------------------------------------------------------
  !> @brief Caculates the single scattering matrix.
  !> @details It is assumed that optic_aot_lognormal is already called, so that the array taus_aer is filled
  ! The output of this subroutine is:
  ! aerosol_optic%dtaus
  ! aerosol_optic%dtaua
  ! aerosol_optic%d_zss_aer
  ! aerosol_optic%zss_aer
  !  --> after combination of modes we get
  ! zss_aer, DIMENSION(nstokes,nstokes,nrt)
  ! aerosol_optic(i)%dtaus_tot_aer
  ! aerosol_optic(i)%dtaua_tot_aer
  ! aerosol_optic(i)%d_zss_tot_aer
  !------------------------------------------------------------------------------
  subroutine optic_phase_ssc( &
       aero_lut, &
       cirrus_lut, &
       iband, &
       iwave, &
       wavelength, &
       scat_angle, &
       aerosol, &
       nder, &
       nrt, &
       aerosol_optic, &
       taus_aer, &
       zss_aer, &
       ierr)
    !*** Input
    type(Mie_lut), intent(in) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    integer, intent(in) :: iband, iwave, nrt
    real(double), dimension(:), intent(in) :: wavelength
    real(double), intent(in) :: scat_angle
    type(aero), dimension(:), intent(in) :: aerosol 
    integer, dimension(:), intent(in) :: nder  
    real(double), dimension(:,:),intent(in) :: taus_aer
    !*** Input/output
    type(aero_opt), dimension(:), intent(inout) :: aerosol_optic    
    !*** Output
    real(double), dimension(:,:,:), intent(out) :: zss_aer   
    integer, intent(out) :: ierr
    !*** local variables
    integer :: i, j, k, l, ik, i_st, j_st, ntype_aer
    real(double) :: rlambda, csca, cabs      
    real(double), dimension(nstokes, nstokes, nfull) :: der_z_ss_full      
    real(double), dimension(nfull) :: dcsca_full, dcabs_full
    real(double), dimension(dim_x) :: rnmb
    real(double), dimension(dim_x, 6) :: der_rnmb      
    real(double), dimension(nstokes, nstokes) :: zss_lay
    real(double), dimension(nstokes, nstokes, npar_mie) :: d_zss_lay    
    real(double), dimension(npar_mie) :: aerosol_pars_mode
    integer :: maxd, nwave
    real(double), dimension(npar_mie) :: dcsca_aer, dcabs_aer 
    character(stringlen) :: message   
    !-----------------------------------------------------------------------------------

    ierr = 0  
    rlambda = wavelength(iwave) * 1.d-3 ! wavelength [um]
    nwave=size(wavelength)
    !*** Calculate optical properties of individual aerosol types
    maxd = size(nder)
    ntype_aer = size(aerosol)
    zss_aer = 0.D0
    do i = 1, ntype_aer
       allocate(aerosol_optic(i)%zss_aer(nstokes, nstokes, nrt), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_PHASE_SSC: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       if(aerosol(i)%CirrusFlag==1) then ! cirrus
          call optic_cirrus_ssc(&
               cirrus_lut, &
               aerosol(i)%tilt_angle,&
               aerosol(i)%shapefrac,&
               aerosol(i)%reff,&
               rlambda,& 
               scat_angle,&
               csca,&
               cabs,&
               zss_lay,&
               dcsca_full,&
               dcabs_full,&
               der_z_ss_full,&
               rnmb,&
               der_rnmb)  
          call der_par_mode_ssc(&
               der_z_ss_full,&
               dcsca_full,&
               dcabs_full,&
               der_rnmb,&
               aerosol(i)%aer_col,&   
               csca,&
               cabs,&
               d_zss_lay,&
               dcsca_aer, &
               dcabs_aer)
       elseif(aerosol(i)%CirrusFlag .ne. 1 )then ! now for aerosols (spheres or ellipsoids)
          aerosol_pars_mode(1) = aerosol(i)%reff
          aerosol_pars_mode(2) = aerosol(i)%veff
          aerosol_pars_mode(3) = aerosol(i)%rm(2*iband-1) + (aerosol(i)%rm(2*iband)-aerosol(i)%rm(2*iband-1))/(nwave-1) *(iwave-1)
          aerosol_pars_mode(4) = aerosol(i)%fim(2*iband-1) + (aerosol(i)%fim(2*iband)-aerosol(i)%fim(2*iband-1))/(nwave-1) *(iwave-1)
          aerosol_pars_mode(5) = aerosol(i)%aer_col
          aerosol_pars_mode(6) = aerosol(i)%shapefrac  !fraction of spheres
          call modes_calc_ssc(& 
               aero_lut, &
               aerosol(i)%id,&
               aerosol_pars_mode,&
               rlambda,&
               scat_angle, &
               dcsca_aer, &
               dcabs_aer,&
               zss_lay,&
               d_zss_lay)
       endif ! aerosol or cirrus
       !*** fill the corresponding arrays  
       aerosol_optic(i)%dtaus_aer = dcsca_aer*aerosol(i)%aer_col
       aerosol_optic(i)%dtaua_aer = dcabs_aer*aerosol(i)%aer_col
       aerosol_optic(i)%d_zss_aer = d_zss_lay
       do ik = 1, maxd
          k = nder(ik)
          aerosol_optic(i)%zss_aer(:,:,k) = zss_lay(:,:)
       enddo
    enddo ! end loop over aerosoltype

    !*** Calculate derivatives of individual aerosol-type properties wrt. composite aerosol properties
    do i = 1, ntype_aer
       if(allocated(aerosol_optic(i)%d_zss_tot_aer)) deallocate(aerosol_optic(i)%d_zss_tot_aer)
       allocate(aerosol_optic(i)%d_zss_tot_aer(4, 4, maxd, npar_mie), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_PHASE_SSC: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       !         aerosol_optic(i)%dtaus_tot_aer = 0.D0
       !         aerosol_optic(i)%dtaua_tot_aer = 0.D0
       aerosol_optic(i)%d_zss_tot_aer = 0.D0
       !         do k= 1, maxd
       !            ik = nder(k)
       !            do l = 1, npar_mie
       !               aerosol_optic(i)%dtaus_tot_aer(ik, l) = aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik)
       !               aerosol_optic(i)%dtaua_tot_aer(ik, l) = aerosol_optic(i)%dtaua_aer(l)*aerosol(i)%alt_dis(ik)
       !            enddo
       !         enddo
       do i_st = 1, nstokes
          do j_st = 1, nstokes         
             do k = 1, maxd
                ik = nder(k)
                do l = 1, npar_mie
                   if(taus_aer(iwave,ik)/=0.D0)then
                      aerosol_optic(i)%d_zss_tot_aer(i_st, j_st, k, l) = &
                           (aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik)* &
                           aerosol_optic(i)%zss_aer(i_st, j_st, ik) + &
                           aerosol_optic(i)%taus_aer(iwave, ik)* &
                           aerosol_optic(i)%d_zss_aer(i_st, j_st, l))/taus_aer(iwave, ik)
                      !***The derivatives of the phase function coefficients of a give mode depend on the other modes
                      do j = 1, ntype_aer                 
                         aerosol_optic(i)%d_zss_tot_aer(i_st, j_st, k, l) = &
                              aerosol_optic(i)%d_zss_tot_aer(i_st, j_st, k, l) - &
                              (aerosol_optic(j)%taus_aer(iwave, ik)*aerosol_optic(j)%zss_aer(i_st, j_st, ik)*&
                              aerosol_optic(i)%dtaus_aer(l)*aerosol(i)%alt_dis(ik)) &
                              /taus_aer(iwave, ik)/taus_aer(iwave, ik)
                      enddo
                   else
                      aerosol_optic(i)%d_zss_tot_aer(i_st, j_st, k, l) = 0.D0
                   endif
                enddo ! l
             enddo ! k
          enddo ! j_st
       enddo ! i_st

       !*** Combine individual aerosol type properties to composite aerosol properties by weighting with optical thickness
       do i_st = 1, nstokes
          do j_st = 1, nstokes
             do k = 1, maxd
                ik = nder(k)
                if(taus_aer(iwave,ik)/=0.D0)then
                   zss_aer(i_st, j_st, ik) = zss_aer(i_st, j_st, ik) + & 
                        (aerosol_optic(i)%taus_aer(iwave,ik)*aerosol_optic(i)%zss_aer(i_st, j_st, ik))/taus_aer(iwave,ik)
                else
                   zss_aer(i_st, j_st, ik) = 0.D0
                endif
             enddo
          enddo
       enddo

    enddo ! end loop over aerosoltype

    do i = 1, ntype_aer
       deallocate(aerosol_optic(i)%zss_aer, stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_PHASE_SSC: memory deallocation error'
          ierr = ierr_deall
          goto 999
       endif
    enddo

    return

999 continue
    call stopretrieval(message)

  end subroutine optic_phase_ssc
  !------------------------------------------------------------------------------  
  !> @details Combine aerosol and molecular optical quantities and
  !! Put aerosol and molecular optical quantities into type aerosol_optic
  !------------------------------------------------------------------------------
  subroutine optic_all_der( &
       iwave, &
       aerosol_optic, aerosol, nder, nrt, maxcoefs, &
       taua_aer, taus_aer, plmom_aer, &
       taua_mol, taus_mol, plmom_ray, &
       taua, taus, plmom_in, &
       ierr)
    integer, intent(in) :: iwave, nrt                               ! which wavelength
    type(aero), dimension(:), intent(in) :: aerosol
    integer, dimension(:), intent(in) :: nder   
    real(double), dimension(:,:,0:,:), intent(in) :: plmom_aer ! Phase matrix coefficients for aerosol scattering  
    integer, dimension(:), intent(in) :: maxcoefs 
    real(double), dimension(:,:), intent(in) :: taua_aer       ! aerosol absorption optical thickness (nwave_hi,nrt)
    real(double), dimension(:,:), intent(in) :: taus_aer       ! aerosol scattering optical thickness (nwave_hi,nrt)
    real(double), dimension(:,:), intent(in) :: taua_mol       ! Molecular absorption optical thickness (nwave_hi,nrt)
    real(double), dimension(:,:), intent(in) :: taus_mol       ! Molecular scattering optical thickness (nwave_hi,nrt)
    real(double), dimension(:,:,0:,:), intent(in) :: plmom_ray ! Phase matrix coefficients for Rayleigh scattering  
    !*** Input/ output
    type(aero_opt), dimension(:), intent(inout) :: aerosol_optic
    !*** Output
    real(double), dimension(:), intent(out) :: taua
    real(double), dimension(:), intent(out) :: taus
    real(double), dimension(nstokes, nstokes, 0:maxleg, nrt), intent(out) :: plmom_in ! Phase matrix coefficients 
    integer, intent(out) :: ierr
    !*** Local variables
    integer :: i, j, k, l, n, i_st, j_st, iper, ik, ntype_aer, maxd
    real(double), dimension(nrt) :: taus_mol_wave, taua_mol_wave
    character(stringlen) :: message
    !------------------------------------------------------------------------------

    ierr = 0
    maxd = size(nder)

    do i = 1, nrt
       taus_mol_wave(i) = taus_mol(iwave,i)
       taua_mol_wave(i) = taua_mol(iwave,i)
       taus(i) = taus_aer(iwave,i) + taus_mol_wave(i)
       taua(i) = taua_aer(iwave,i) + taua_mol_wave(i)
    enddo

    ntype_aer = size(aerosol)
    do i = 1, ntype_aer    
       if(allocated(aerosol_optic(i)%dtaus_all)) deallocate(aerosol_optic(i)%dtaus_all)
       if(allocated(aerosol_optic(i)%dtaua_all)) deallocate(aerosol_optic(i)%dtaua_all)
       allocate(&
            aerosol_optic(i)%dtaus_all(maxd, npar),&
            aerosol_optic(i)%dtaua_all(maxd, npar), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_ALL_DER: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       aerosol_optic(i)%dtaus_all = 0.D0
       aerosol_optic(i)%dtaua_all = 0.D0

       do l = 1, npar_mie
          do ik = 1, maxd
             k = nder(ik)
             aerosol_optic(i)%dtaus_all(ik, l) = aerosol_optic(i)%dtaus_tot_aer(k, l)
             aerosol_optic(i)%dtaua_all(ik, l) = aerosol_optic(i)%dtaua_tot_aer(k, l)
          enddo
       enddo

       do ik = 1, maxd
          k = nder(ik)
          aerosol_optic(i)%dtaus_all(ik, npar_mie+1) = &
               sum(aerosol_optic(i)%taus_aer(iwave,:))*aerosol(i)%dalt_daer1(k)
          aerosol_optic(i)%dtaua_all(ik, npar_mie+1) = &
               sum(aerosol_optic(i)%taua_aer(iwave,:))*aerosol(i)%dalt_daer1(k)
          aerosol_optic(i)%dtaus_all(ik, npar_mie+2) = &
               sum(aerosol_optic(i)%taus_aer(iwave,:))*aerosol(i)%dalt_daer2(k)
          aerosol_optic(i)%dtaua_all(ik, npar_mie+2) = &
               sum(aerosol_optic(i)%taua_aer(iwave,:))*aerosol(i)%dalt_daer2(k)
       enddo
    enddo ! end loop over aerosol types

    !*** Add Rayleigh scattering contribution to composite phase matrix by weighting with optical thickness   
    plmom_in = 0.0d0
    do n = 0, 2
       do i_st = 1, nstokes
          do j_st = 1, nstokes
             do k = 1, nrt
                if(taus(k) == 0.D0) then
                   plmom_in(i_st,j_st,n,k) = 0.D0
                else 
                   plmom_in(i_st,j_st,n,k) = & 
                        (taus_mol_wave(k)*plmom_ray(i_st,j_st,n,k) + &
                        taus_aer(iwave,k)*plmom_aer(i_st,j_st,n,k)) /&
                        (taus(k))
                endif
             enddo
          enddo
       enddo
    enddo

    do i = 1, ntype_aer       
       if(allocated(aerosol_optic(i)%dphase_all))deallocate(aerosol_optic(i)%dphase_all)
       allocate(aerosol_optic(i)%dphase_all(nper,0:maxleg,maxd,npar), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'OPTIC_ALL_DER: memory allocation error'
          ierr = ierr_all
          goto 999
       endif
       aerosol_optic(i)%dphase_all = 0.D0
       do iper = 1, nper
          i_st = index_ist(iper)
          j_st = index_jst(iper)
          if (i_st .le. nstokes .and. j_st .le. nstokes) then
             do ik = 1, maxd
                k = nder(ik)
                do n = 0, 2
                   do l = 1, npar_mie               
                      if(taus(k)==0.D0)then
                         aerosol_optic(i)%dphase_all(iper,n,ik,l) = 0.D0
                      else 
                         aerosol_optic(i)%dphase_all(iper,n,ik,l) = &
                              (aerosol_optic(i)%dtaus_all(ik,l)*plmom_aer(i_st,j_st,n,k) + &
        		      aerosol_optic(i)%dphase_tot_aer(iper,n,ik,l)*taus_aer(iwave,k))/(taus(k)) + &
                              (taus_aer(iwave,k)*plmom_aer(i_st,j_st,n,k) + &
        		      taus_mol_wave(k)*plmom_ray(i_st,j_st,n,k))* &
        		      aerosol_optic(i)%dtaus_all(ik,l)* &
                              (-(taus(k))**(-2))
                      endif
                   enddo
                   do l = npar_mie+1, npar_mie+2
                      if(taus(k)==0.D0)then
                         aerosol_optic(i)%dphase_all(iper,n,ik,l) = 0.D0
                      else 
                         aerosol_optic(i)%dphase_all(iper,n,ik,l) = &
                              (aerosol_optic(i)%dtaus_all(ik,l)*plmom_aer(i_st,j_st,n,k) /&
                              (taus(k)))&
                              +(taus_aer(iwave,k)*plmom_aer(i_st,j_st,n,k) +&
                              taus_mol_wave(k)*plmom_ray(i_st,j_st,n,k))*&
                              aerosol_optic(i)%dtaus_all(ik,l)*&
                              (-(taus(k))**(-2))
                      endif
                   enddo  ! l             
                enddo ! n
             enddo ! ik
          endif
       enddo ! iper

       do j_st = 1, nstokes
          do i_st = 1, nstokes
             do j = 1, maxd
                k = nder(j)
                do n = 3, maxcoefs(k)
                   if(taus(k)==0.D0) then
                      plmom_in(i_st, j_st, n, k) = 0.D0
                   else 
                      plmom_in(i_st, j_st, n, k) = &
                           (taus_aer(iwave,k)*plmom_aer(i_st, j_st, n, k))/(taus(k))
                   endif
                enddo
             enddo
          enddo
       enddo

       do iper = 1, nper
          i_st = index_ist(iper)
          j_st = index_jst(iper)
          if (i_st .le. nstokes .and. j_st .le. nstokes) then 
             do ik = 1, maxd
                k = nder(ik)
                do n = 3, maxstr
                   do l = 1, npar_mie
                      if(taus(k)==0.D0) then
                         aerosol_optic(i)%dphase_all(iper, n, ik, l) = 0.D0
                      else 
                         aerosol_optic(i)%dphase_all(iper ,n, ik ,l) = & 
                              (aerosol_optic(i)%dtaus_all(ik,l)*plmom_aer(i_st, j_st, n, k) + &
                              aerosol_optic(i)%dphase_tot_aer(iper, n, ik, l)*taus_aer(iwave, k))/ &
                              (taus(k))&
                              +(taus_aer(iwave,k)*plmom_aer(i_st, j_st, n, k))*&
                              aerosol_optic(i)%dtaus_all(ik, l)*&
                              (-(taus(k))**(-2))
                      endif
                   enddo
                   do l = npar_mie+1, npar_mie+2
                      if(taus(k)==0.D0) then
                         aerosol_optic(i)%dphase_all(iper, n, ik, l) = 0.D0
                      else 
                         aerosol_optic(i)%dphase_all(iper, n, ik, l) = &
                              (aerosol_optic(i)%dtaus_all(ik,l)*plmom_aer(i_st, j_st, n, k) / &
                              (taus(k))) &
                              +(taus_aer(iwave,k)*plmom_aer(i_st, j_st, n, k))* &
                              aerosol_optic(i)%dtaus_all(ik, l)* &
                              (-(taus(k))**(-2))
                      endif
                   enddo ! l
                enddo ! n
             enddo ! ik
          endif
       enddo ! iper

    enddo ! ntype_aer

    return

999 continue
    call stopretrieval(message)

  end subroutine optic_all_der
  !------------------------------------------------------------------------------
  !> @details Compute Z matrix for single scattering, and derivatives of coefficients 
  !! of Z matrix wrt aerosol parameters (nstokes, maxd, npar)
  !------------------------------------------------------------------------------
  subroutine optic_all_ssc( &
       aero_lut, &
       cirrus_lut, &
       rtmflag, &
       iband, wavelength, &
       scat_angle, u0, uv, phi_in, &
       aerosol, nder, nrt, aerosol_optic, &
       taus_aer, taus_mol, z_ss_int, &
       zss_all_lin, &
       ierr)   
    !*** Input
    type(Mie_lut) :: aero_lut
    type(cirrus_table), intent(in) :: cirrus_lut
    integer, intent(in) :: rtmflag, iband, nrt
    real(double), dimension(:), intent(in) :: wavelength
    real(double), intent(in) :: scat_angle, u0, uv, phi_in
    type(aero), dimension(:), intent(in) :: aerosol
    integer, dimension(:), intent(in) :: nder   
    real(double), dimension(:,:), intent(in) :: taus_aer
    real(double), dimension(:,:), intent(in) :: taus_mol 
    !*** Input/output
    type(aero_opt), dimension(:), intent(inout) :: aerosol_optic
    !*** Output
    real(double), dimension(:,:,:), intent(out) :: z_ss_int  
    real(double), dimension(:, :, :, :), intent(out) :: zss_all_lin
    integer, intent(out) :: ierr
    !*** Local variables
    integer :: i, j, k, i_st, j_st, ik, iwave, imid, iw, ix, ntype_aer, nwave
    real(double), dimension(4, 4) :: zray
    integer, dimension(3) :: index_wave
    real(double), dimension(4, 4) :: F_rot1, F_rot2, Zmat
    real(double) :: u_in, u_out, cos_scat, cos_i1, cos_i2, ang_i1, ang_i2, C0, wave, sin_scat, phi
    real(double), dimension(4, 4, nrt) :: zss_all
    real(double), dimension(:, :, :, :), allocatable :: d_zss_all
    real(double) :: depol 
    real(double), dimension(nstokes, nstokes, nrt) :: zss_aer  
    integer :: maxd
    character(stringlen) :: message
    !------------------------------------------------------------------------------

    ierr = 0
    nwave = size(wavelength)
    maxd = size(nder) 
    allocate (d_zss_all(4, 4, maxd, npar), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'OPTIC_ALL_SSC: memory allocation error'
       ierr = ierr_all
       goto 999
    endif
    d_zss_all = 0.d0
    zss_all = 0.d0
    !*** calculate single scat. properties at 3 spectral points (2 edges and middle of band)
    imid = (nwave+1)/2 ! imid
    index_wave(1) = 1
    index_wave(2) = imid
    index_wave(3) = nwave

    if(nstokes.gt.1 .and. rtmflag==1) then
       !*** prepare rotation matrices for rotation to local meridian plane
       u_out = uv
       u_in = -u0
       cos_scat = DCOS(scat_angle)
       phi = phi_in/ 180.D0*PI
       cos_i1 = 1.D0
       cos_i2 = 1.D0
       if(cos_scat**2.ne.1)then       
          cos_i1 = ( &
               ( u_out*DSQRT(1.D0-u_in**2) - &
               u_in *DSQRT(1.D0-u_out**2)*DCOS(PHI) ) / &
               (DSQRT(1.D0-cos_scat**2)) )       
          cos_i2 = ( &
               ( u_in*DSQRT(1.D0-u_out**2) - &
               u_out *DSQRT(1.D0-u_in**2)*DCOS(PHI) ) / &
               (DSQRT(1.D0-cos_scat**2)) )                   
          cos_i1 = min(cos_i1,1.D0)
          cos_i2 = min(cos_i2,1.D0)      
          cos_i1 = max(cos_i1,-1.D0)
          cos_i2 = max(cos_i2,-1.D0)
       endif
       ang_i1 = DACOS(cos_i1)
       ang_i2 = DACOS(cos_i2)
       if(phi.lt.0.)then
          ang_i1 = -ang_i1
          ang_i2 = -ang_i2
       endif
       F_rot1 = 0.D0              !initialize         
       F_rot1(1,1) = 1.D0
       F_rot1(2,2) = DCOS(2*ang_i2)
       F_rot1(3,3) = DCOS(2*ang_i2)
       F_rot1(2,3) = -DSIN(2*ang_i2)
       F_rot1(3,2) = DSIN(2*ang_i2)

       F_rot2 = 0.D0              !initialize         
       F_rot2(1,1) = 1.D0
       F_rot2(2,2) = DCOS(2*ang_i1)
       F_rot2(3,3) = DCOS(2*ang_i1)
       F_rot2(2,3) = -DSIN(2*ang_i1)
       F_rot2(3,2) = DSIN(2*ang_i1)
    endif

    ntype_aer = size(aerosol_optic)
    do i = 1, ntype_aer
       if(allocated(aerosol_optic(i)%dzss_all)) deallocate(aerosol_optic(i)%dzss_all)
       if(allocated(aerosol_optic(i)%dzss_all_1)) deallocate(aerosol_optic(i)%dzss_all_1)
       if(allocated(aerosol_optic(i)%dzss_all_mid)) deallocate(aerosol_optic(i)%dzss_all_mid)
       if(allocated(aerosol_optic(i)%dzss_all_nwave)) deallocate(aerosol_optic(i)%dzss_all_nwave)
       allocate(aerosol_optic(i)%dzss_all(nstokes, maxd, npar), & 
            aerosol_optic(i)%dzss_all_1(1:nstokes, 1:maxd, 1:npar), &
            aerosol_optic(i)%dzss_all_mid(1:nstokes, 1:maxd, 1:npar), &
            aerosol_optic(i)%dzss_all_nwave(1:nstokes, 1:maxd, 1:npar), &
            stat = ierr)
       if (ierr .ne. 0) then
          write(message,*) 'OPTIC_ALL_SSC: memory alloaction error'
          ierr = ierr_all
          goto 999
       endif
       aerosol_optic(i)%dzss_all = 0.D0
       aerosol_optic(i)%dzss_all_1 = 0.D0
       aerosol_optic(i)%dzss_all_mid = 0.D0
       aerosol_optic(i)%dzss_all_nwave = 0.D0
    enddo

    cos_scat = DCOS(scat_angle) 
    sin_scat = DSIN(scat_angle)
    wave = wavelength(imid)
    c0 = (2. - 2.D0*depol(wave))/(2.D0 + depol(wave))

    zray = 0.D0
    zray(1,1) = c0*0.75D0*(1.D0+cos_scat*cos_scat) +(1.D0-C0)
    zray(1,2) = c0*(-0.75D0)*sin_scat*sin_scat
    zray(2,1) = zray(1,2)
    zray(2,2) = c0*0.75D0*(1.D0+cos_scat*cos_scat)
    zray(3,3) = c0*1.5D0*cos_scat

    do iw = 1, 3
       iwave = index_wave(iw)
       call optic_phase_ssc( &
            aero_lut, cirrus_lut, &
            iband, iwave, wavelength, scat_angle, &
            aerosol, nder, nrt, aerosol_optic, taus_aer, zss_aer, &
            ierr)

       do k = 1, nrt
          !*** combine aerosol and Rayleigh phase matrix into zss_all
          do i_st = 1, nstokes
             do j_st = 1, nstokes
                zss_all(i_st,j_st,k) = &
                     (taus_mol(iwave,k)*zray(i_st,j_st) + &
                     taus_aer(iwave,k)*zss_aer(i_st,j_st,k)) / &
                     (taus_mol(iwave,k)+taus_aer(iwave,k))
    		zss_all_lin(i_st, j_st,k,iwave) = zss_all(i_st,j_st,k)
             enddo
          enddo
       enddo


       !*** For scalar write element 1,1 of the scattering matrix to the array z_ss_int
       !*** For polarization, rotation to the local meridian plane is needed if rtmflag>1.
       if (nstokes.eq.1) then       
          do k = 1, nrt
             z_ss_int(1, k, iwave) = zss_all(1, 1, k)
          enddo
       elseif(rtmflag==1) then
          !*** rotation to local meridian plane          
          do k = 1, nrt
             Zmat(:,:) = matmul(F_rot1,(matmul(zss_all(:,:,k), F_rot2)))
             do i_st = 1, nstokes
                z_ss_int(i_st, k, iwave) = Zmat(i_st, 1)
             enddo
          enddo
       elseif (rtmflag .gt. 1) then            
          do k = 1, nrt
             z_ss_int(1, k, iwave) = zss_all(1, 1, k)           
             z_ss_int(2, k, iwave) = zss_all(2, 2, k)
             z_ss_int(3, k, iwave) = zss_all(3, 3, k)		
             z_ss_int(4, k, iwave) = zss_all(2, 1, k)					
          enddo
       endif



       !*** now derivatives, only for middle wavelength
       !       if(iwave.eq.imid) then
       do i = 1, ntype_aer
          do ik = 1, maxd
             k = nder(ik)
             do i_st = 1, nstokes
                do j_st = 1, nstokes
                   do ix = 1, npar_mie
                      d_zss_all(i_st, j_st, ik, ix) = &
                           (aerosol_optic(i)%dtaus_all(ik,ix)*zss_aer(i_st,j_st,k) + &
                           aerosol_optic(i)%d_zss_tot_aer(i_st,j_st,ik,ix)*taus_aer(iwave,k)) / &
                           (taus_mol(iwave,k)+taus_aer(iwave,k)) + &
                           (taus_aer(iwave,k)*zss_aer(i_st,j_st,k) +&
                           taus_mol(iwave,k)*zray(i_st,j_st)) * &
                           aerosol_optic(i)%dtaus_all(ik,ix)*&
                           (-((taus_mol(iwave,k)+taus_aer(iwave,k))**(-2)))
                   enddo
                   do ix = npar_mie+1, npar_mie+2
                      d_zss_all(i_st, j_st, ik, ix) = &
                           (aerosol_optic(i)%dtaus_all(ik,ix)*zss_aer(i_st,j_st,k) / &
                           (taus_mol(iwave,k)+taus_aer(iwave,k)))  + &
                           (taus_aer(iwave,k)*zss_aer(i_st,j_st,k) + &
                           taus_mol(iwave,k)*zray(i_st,j_st)) * &
                           aerosol_optic(i)%dtaus_all(ik,ix)* &
                           (-((taus_mol(iwave,k)+taus_aer(iwave,k))**(-2)))
                   enddo
                enddo
             enddo
          enddo
          
          !*** rotation to local mer. plane for derivatives if polarization and rtmflag>1
          if(nstokes.eq.1) then
             do ik = 1, maxd
                do ix = 1, npar
                   if(iwave.eq.imid) then
                      aerosol_optic(i)%dzss_all_mid(1, ik, ix) = d_zss_all(1, 1, ik, ix)
                   elseif(iwave.eq.1) then 
                      aerosol_optic(i)%dzss_all_1(1,ik,ix) = d_zss_all(1, 1, ik, ix)
                   elseif(iwave.eq.nwave) then 
                      aerosol_optic(i)%dzss_all_nwave(1,ik,ix) = d_zss_all(1, 1, ik, ix)
                   endif
                enddo
             enddo
          elseif(rtmflag==1)then
             do ik = 1, maxd
                do ix = 1,npar
                   Zmat(:,:) = matmul(F_rot1,(matmul(d_zss_all(:,:,ik,ix),F_rot2)))
                   do i_st = 1, nstokes
                      if(iwave.eq.imid) then
                         aerosol_optic(i)%dzss_all_mid(i_st,ik,ix) = Zmat(i_st,1)
                      elseif(iwave.eq.1) then 
                         aerosol_optic(i)%dzss_all_1(i_st,ik,ix) = Zmat(i_st,1)
                      elseif(iwave.eq.nwave) then 
                         aerosol_optic(i)%dzss_all_nwave(i_st,ik,ix) = Zmat(i_st,1)
                      endif
                   enddo
                enddo
             enddo
          elseif (rtmflag>1) then
             do ik = 1, maxd
                do ix = 1, npar
                   if(iwave.eq.imid) then
                      aerosol_optic(i)%dzss_all_mid(1, ik, ix) = d_zss_all(1, 1, ik, ix)
                      aerosol_optic(i)%dzss_all_mid(2, ik, ix) = d_zss_all(1,2,ik, ix)
                      aerosol_optic(i)%dzss_all_mid(3, ik, ix) = d_zss_all(1,2,ik, ix)		                           
                   elseif(iwave.eq.1) then 
                      aerosol_optic(i)%dzss_all_1(1,ik,ix) = d_zss_all(1, 1, ik, ix)
                      aerosol_optic(i)%dzss_all_1(2, ik, ix) = d_zss_all(1,2,ik, ix)   
                      aerosol_optic(i)%dzss_all_1(3, ik, ix) = d_zss_all(1,2,ik, ix)		                         
                   elseif(iwave.eq.nwave) then 
                      aerosol_optic(i)%dzss_all_nwave(1,ik,ix) = d_zss_all(1, 1, ik, ix)
                      aerosol_optic(i)%dzss_all_nwave(2, ik, ix) = d_zss_all(1,2,ik, ix)
                      aerosol_optic(i)%dzss_all_nwave(3, ik, ix) = d_zss_all(1,2,ik, ix)		      
                   endif
                enddo
             enddo
          endif

       enddo ! i = 1, ntype_aer
       !       endif ! iwave
    enddo ! iw = 1, 3

    !*** Only edges and middle of band are filled, linearly interpolate in between
    do j = 1, nrt
       do i_st = 1, nstokes
          do iwave = 2, imid-1
             z_ss_int(i_st, j, iwave) = z_ss_int(i_st, j, 1) + &
                  (z_ss_int(i_st, j, imid) - z_ss_int(i_st, j, 1))/&
                  (wavelength(imid)-wavelength(1)) * &
                  (wavelength(iwave)-wavelength(1))
          enddo
          do iwave = imid+1, nwave-1
             z_ss_int(i_st,j,iwave) = z_ss_int(i_st,j,imid) + &
                  (z_ss_int(i_st,j,nwave) - z_ss_int(i_st,j,imid))/&
                  (wavelength(nwave)-wavelength(imid)) * &
                  (wavelength(iwave)-wavelength(imid))
          enddo
       enddo
    enddo

    !*** Only edges and middle of band are filled, linearly interpolate in between
    !*** prepare for lintran single scattering
    do j = 1, nrt
       do i_st = 1, nstokes
          do j_st = 1, nstokes
             do iwave = 2, imid-1
                zss_all_lin(i_st, j_st, j, iwave) = zss_all_lin(i_st, j_st, j, 1) + &
                     (zss_all_lin(i_st, j_st, j, imid) - zss_all_lin(i_st, j_st, j, 1))/&
                     (wavelength(imid)-wavelength(1)) * &
                     (wavelength(iwave)-wavelength(1))
             enddo
             do iwave = imid+1, nwave-1
                zss_all_lin(i_st, j_st, j,iwave) = zss_all_lin(i_st, j_st, j,imid) + &
                     (zss_all_lin(i_st, j_st, j,nwave) - zss_all_lin(i_st, j_st,j,imid))/&
                     (wavelength(nwave)-wavelength(imid)) * &
                     (wavelength(iwave)-wavelength(imid))
             enddo
          enddo
       enddo
    enddo
    return

999 continue
    call stopretrieval(message)

  end subroutine optic_all_ssc
  !------------------------------------------------------------------------------

end module optic_input_module
