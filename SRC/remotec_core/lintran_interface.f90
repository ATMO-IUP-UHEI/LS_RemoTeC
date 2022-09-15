module lintran_interface
  use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class
  use lintran_module, only: lintran_init, lintran_provide, lintran_calculate, lintran_close, lintran_nakajima_multiscattering
  use lintran_remotec_module
  use header_module
  implicit none
  !private
contains  
  !------------------------------------------------------------------------------
  subroutine lintran_allocate( nrt, & ! {{{
       maxd, &
       atm, &
       drv, &
       set, &
       ierr)
    !*** input
    integer, intent(in) ::  nrt, maxd
    !in-/output
    type(lintran_atmosphere), intent(inout) :: atm
    type(lintran_derivatives), intent(inout) :: drv 
    type(lintran_settings), intent(inout) :: set
    !*** output
    integer, intent(out) :: ierr
    !*** local 
    integer :: nlay_deriv_base 
    integer :: nlay_deriv_ph
    integer :: nleg_deriv 

    !*** Initialize error identifier
    ierr = 0

    nlay_deriv_base = nrt !maxd
    nlay_deriv_ph = maxd!nrt		
    nleg_deriv = maxleg
    allocate(atm%taus(nrt), &
         atm%taua(nrt), &
                  atm%coefs(nper, 0:maxleg, nrt), &
         !atm%coefs(nper, 0:maxstr, nrt), &
         atm%phase_ssg(nper, nrt,1), &
         atm%bdrf(nstokes, nstokes, MXHALF, MXHALF, 0:maxstr-1), &
         atm%bdrf_0(nstokes, nstokes, MXHALF, 0:maxstr-1), &
         atm%bdrf_v(nstokes, nstokes, MXHALF, 0:maxstr-1,1), &
         atm%bdrf_ssg(nstokes, nstokes, 1), &
         atm%emi(nstokes, MXHALF, 0:maxstr-1), &
         atm%emi_ssg(nstokes, 1), &
         atm%nleg_lay(nrt), &
         set%ilay_deriv_base(nlay_deriv_base), &
         set%ilay_deriv_ph(nlay_deriv_ph), &
         set%execute_stokes(nstokes), &               
         set%differentiate_stokes(nstokes), &               
         drv%drint_taus(nlay_deriv_base, nstokes, 1), &
         drv%drint_taua(nlay_deriv_base, nstokes, 1), &
         drv%drint_coefs(nper, 0:nleg_deriv, nlay_deriv_ph, nstokes, 1), &
         drv%drint_phase_ssg(nper, nlay_deriv_ph, nstokes, 1), &
         !drv%drint_bdrf(nstokes, nstokes, MXHALF, MXHALF, 0:maxleg, nstokes, 1), &
         !drv%drint_bdrf_0(nstokes, nstokes, MXHALF, 0:maxleg, nstokes, 1), &
         !drv%drint_bdrf_v(nstokes, nstokes, MXHALF, 0:maxleg, nstokes, 1), &
         drv%drint_bdrf(nstokes, nstokes, MXHALF, MXHALF, 0:maxstr-1, nstokes, 1), &
         drv%drint_bdrf_0(nstokes, nstokes, MXHALF, 0:maxstr-1, nstokes, 1), &
         drv%drint_bdrf_v(nstokes, nstokes, MXHALF, 0:maxstr-1, nstokes, 1), &        
         drv%drint_bdrf_ssg(nstokes, nstokes, nstokes, 1), &
         drv%drint_emi(nstokes, MXHALF, 0:maxstr-1, nstokes, 1), &
         drv%drint_emi_ssg(nstokes, nstokes, 1), &
         stat = ierr)
    if (ierr .ne. 0) then
       ierr = ierr_all
       call stopretrieval('LINTRAN_ALLOCATE: memory allocation error')
    endif

    set%nleg_deriv =  nleg_deriv

    return

  end subroutine lintran_allocate ! }}}

  !------------------------------------------------------------------------------

  subroutine lintranv2_assign(nst, & ! {{{
       nrt, &
       maxd, &
       nder, &
       cloudflag, &
       maxcoefs, &
       plmom_in, &
       bdrf_ms, &
       bdrf_ss, &
       !       z_ss, & 
       zss_all, & 
       surf_emi, &
       taua, & 
       taus, &
       atm, &
       drv, &
       set )
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class
    use header_module
    use lintran_header
    !input
    integer, intent(in) :: nst, nrt, maxd
    integer, dimension(maxd), intent(in) :: nder
    logical, intent(in) :: cloudflag
    integer, dimension(nrt), intent(in) :: maxcoefs
    real(double), intent(in) :: surf_emi
    real(double), dimension(nst), intent(in) :: bdrf_ss
    !    real(double), dimension(nst,nrt), intent(in) :: z_ss
    real(double), dimension(4, 4, nrt) :: zss_all ! Z matrix for Lintran
    real(double), dimension(nrt), intent(in) ::taua, taus
    real(double), dimension(nst, nst, 0:maxleg, nrt), intent(in) :: plmom_in
    real(double), dimension(nst, nst, mxhalf+2, mxhalf+2, 0:maxstr-1), intent(in) :: bdrf_ms
    !in-/output
    type(lintran_atmosphere), intent(inout) :: atm
    type(lintran_derivatives), intent(inout) :: drv 
    type(lintran_settings), intent(inout) :: set
    !*** lical
    integer :: ilay, i

    !*** Indices to convert scattering matrix to array  
    integer, dimension(4) :: index_ist = (/1,1,2,3/)
    integer, dimension(4) :: index_jst = (/1,2,2,3 /)

    !========================================================================================
    ! Populate the Lintran atmosphere. For now, we know that we have a sun.
    atm%sun = 1.0
    atm%thermal_emission = .false.
    atm%bdrf_only_0 = .true. !glintflag .ne. 1 ! Ingore higher Fourier modes for surface reflection.
    atm%emi_only_0 = .true. ! Ingore higher Fourier modes for fluorescence.    

    atm%taua = taua
    atm%taus = taus 
    atm%coefs = 0.D0
    do i = 1, nper
       atm%coefs(i, 0:maxleg, :) = plmom_in(index_ist(i), index_jst(i), 0:maxleg, :) 
       atm%phase_ssg(i, 1:nrt, 1) = zss_all(index_ist(i), index_jst(i), 1:nrt)
    enddo

    atm%bdrf = 0D0
    atm%bdrf_0 = 0D0
    atm%bdrf_v = 0D0
    atm%bdrf_ssg = 0D0
    atm%bdrf(1:nst, 1:nst, 1:mxhalf, 1:mxhalf, 0:MAXSTR-1) = BDRF_MS(1:nst, 1:nst, 1:MXHALF, 1:MXHALF, 0:MAXSTR-1)
    atm%bdrf_0(1:nst, 1:nst, 1:mxhalf, 0:MAXSTR-1) = BDRF_MS(1:nst, 1:nst, 1:mxhalf, MXHALF+1, 0:MAXSTR-1)
    atm%bdrf_v(1:nst, 1:nst, 1:mxhalf, 0:MAXSTR-1,1) = BDRF_MS(1:nst, 1:nst, MXHALF+2, 1:mxhalf, 0:MAXSTR-1)

    atm%bdrf_ssg(1,1,1) = bdrf_ss(1)
    atm%emi = 0.0
    atm%emi(1, :, 0) = surf_emi

    atm%emi_ssg = 0D0
    atm%emi_ssg(1, 1) = surf_emi

    atm%nleg_lay = min(maxcoefs,maxstr) ! Not maxstr-1, because the number maxstr will be used to apply Delta-M internally.
    atm%bdrf_only_0 = .true. ! Ingore higher Fourier modes for surface reflection.
    atm%emi_only_0 = .true. ! Ingore higher Fourier modes for fluorescence.

    set%execute_stokes = .true.
    set%differentiate_stokes = .false.
    set%differentiate_stokes(1) = .true.
    set%taua_split = taua_split 
    set%tautot_max = tautot_max 
    set%fourier_tolerance = fourier_tolerance
    set%gs_tolerance = gs_tolerance
    set%gs_maxiter = gs_maxiter

    if (cloudflag) then
       set%taus_split = taus_split_cld
       set%split_double = split_double_cld
       set%interpolation = interpolation_cld
       set%solver = solver_cld
    else
       set%taus_split = taus_split_clr
       set%split_double = split_double_clr
       set%interpolation = interpolation_clr
       set%solver = solver_clr
    endif

    set%nleg_deriv =  maxleg
    set%ilay_deriv_base = (/(ilay, ilay=1,nrt)/)
    set%ilay_deriv_ph = nder

    drv%drint_taus = 0D0
    drv%drint_taua = 0D0
    drv%drint_coefs = 0D0
    drv%drint_phase_ssg = 0D0
    drv%drint_bdrf = 0D0
    drv%drint_bdrf_0 = 0D0
    drv%drint_bdrf_v = 0D0
    drv%drint_bdrf_ssg = 0D0
    drv%drint_emi = 0D0
    drv%drint_emi_ssg = 0D0

  end subroutine lintranv2_assign ! }}}

  !------------------------------------------------------------------------------
  subroutine lintran_deallocate(atm, set, drv, ierr) ! {{{
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class
    !in-/output
    type(lintran_atmosphere), intent(inout) :: atm
    type(lintran_derivatives), intent(inout) :: drv
    type(lintran_settings), intent(inout) :: set
    !*** output
    integer, intent(out) :: ierr
 
    !*** Initialize error identifier 
    ierr = 0
 
    deallocate(atm%taus, &
         atm%taua, &
         atm%coefs, &
         atm%phase_ssg, &
         atm%bdrf, &
         atm%bdrf_0, &
         atm%bdrf_v, &
         atm%bdrf_ssg, &
         atm%emi, &
         atm%emi_ssg, &
         atm%nleg_lay, &
         set%ilay_deriv_base, &
         set%ilay_deriv_ph, &
         set%execute_stokes, &
         set%differentiate_stokes, &
         drv%drint_taus, &
         drv%drint_taua, &
         drv%drint_coefs, &
         drv%drint_phase_ssg, &
         drv%drint_bdrf, &
         drv%drint_bdrf_0, &
         drv%drint_bdrf_v, &
         drv%drint_bdrf_ssg, &
         drv%drint_emi, &
         drv%drint_emi_ssg, &
         stat=ierr)
    if (ierr .ne. 0) then
       ierr = 4
       call stopretrieval('LINTRAN_DEALLOCATE: memory deallocation error')
    endif 
   

  end subroutine lintran_deallocate ! }}}
!------------------------------------------------------------------------------
  subroutine lintranv2_return( nst, nrt, maxd, dz_ss, dalb_ss, dtaua, dtaus, dphase, dalb, drv ) ! {{{
    use header_module
    use lintran_types_module, only: lintran_atmosphere, lintran_settings, lintran_derivatives, lintran_class
    !in-/output
    type(lintran_derivatives), intent(inout) :: drv
    integer, intent(in) :: nst, nrt, maxd            
    real(double),  dimension(nst, maxd), intent(inout) :: dz_ss
    real(double),  dimension(nst), intent(inout) :: dalb_ss
    real(double),  dimension(nst,nrt), intent(inout) :: dtaua, dtaus        
    real(double),  dimension(nst,0:maxleg,nper,maxd), intent(inout) :: dphase  
    real(double),  dimension(nst), intent(inout) :: dalb
    integer :: i, k, l, n 
    !========================================================================================


    dalb_ss(1:nst)        = drv%drint_bdrf_ssg(1,1,1:nst,1)

    do n =1,nst
       dz_ss(n,1:maxd)       = drv%drint_phase_ssg(min(n,2),1:maxd,n,1)
    enddo
    
    do n = 1,nst
       dtaua(n,1:nrt)    	= drv%drint_taua(1:nrt,n,1)
       dtaus(n,1:nrt)    	= drv%drint_taus(1:nrt,n,1)
    enddo

    do n = 1, nper
       do k= 1, maxd
          do l = 0, maxleg
             do i = 1, nst
                dphase(i,l,n,k)   = drv%drint_coefs(n, l, k, i, 1)
             enddo
          enddo
       enddo
    enddo

    do i = 1, nst
       dalb(i) = sum(drv%drint_bdrf(1,1,:,:,0,i,1)) +&
            sum(drv%drint_bdrf_0(1,1,:,0,i,1)) +&
            sum(drv%drint_bdrf_v(1,1,:,0,i,1)) +&
            drv%drint_bdrf_ssg(1,1,i,1)
    enddo

  end subroutine lintranv2_return ! }}}
  !------------------------------------------------------------------------------
  subroutine lintranv2_return_ss( nst, nrt, maxd, dz_ss, dalb_ss, dtaua, dtaus, drv ) ! {{{
    !in-/output
    type(lintran_derivatives), intent(inout) :: drv
    integer, intent(in) :: nst, nrt, maxd

    ! real,  dimension(nst,nrt), intent(inout) :: dtaua_ss, dtaus_ss
    integer :: n
    real(double),  dimension(nst, maxd), intent(inout) :: dz_ss
    real(double),  dimension(nst), intent(inout) :: dalb_ss
    real(double),  dimension(nst,nrt), intent(inout) :: dtaua, dtaus        
    !========================================================================================


    do n = 1, nst
       dalb_ss(n)        = drv%drint_bdrf_ssg(1,1,n,1)
       dz_ss(n,1:maxd)       = drv%drint_phase_ssg(min(n,2),1:maxd,n,1)
       dtaua(n,1:nrt)    	= drv%drint_taua(1:nrt,n,1)
       dtaus(n,1:nrt)    	= drv%drint_taus(1:nrt,n,1)
    enddo

  end subroutine lintranv2_return_ss ! }}}


!----------------------------------------------------------------------

end module lintran_interface  
