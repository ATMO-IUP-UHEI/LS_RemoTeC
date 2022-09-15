module tau_grid_module
  use header_module
  implicit none
  private

  !*** Create generic interface (function overloading)
  interface set_tau_grid
     module procedure set_tau_grid_1species
     module procedure set_tau_grid_2species
  end interface set_tau_grid

  public :: interpolate_lbl_2species, interpolate_lbl_1species, set_tau_grid

  !*** Make specific procedures private
  private :: set_tau_grid_2species, set_tau_grid_1species, &
       interpolate_tau_2nd, interpolate_tau_2nd_pt, interpolate_tau_1st

contains
  !*****************************************************************

  subroutine set_tau_grid_2species( &
       nwave, &
       nrt, &
       imol_1st, &
       imol_2nd, &
       ntau, &
       ntau_2nd_max, &
       taus, & 
       taua_tot_species, &
       taua_mol, &
       taua_mol_species, &       
       ntau_2nd, tau_grid, tau_grid_2nd, tau_grid_arr, &
       ierr)
    !*** Input
    integer, intent(in) :: nwave, nrt, imol_1st, imol_2nd, ntau, ntau_2nd_max
    real(double), dimension(:), intent(in) :: taus
    real(double), dimension(:,:),intent(in) :: taua_mol     
    real(double), dimension(:,:), intent(in) :: taua_tot_species
    real(double), dimension(:,:,:), intent(in) ::  taua_mol_species 
    !*** Output
    integer, dimension(:), allocatable, intent(out) :: ntau_2nd   
    real(double), dimension(:), allocatable, intent(out) :: tau_grid
    real(double), dimension(:,:), allocatable, intent(out) :: tau_grid_2nd
    real(double), dimension(:,:,:), allocatable, intent(out) :: tau_grid_arr
    integer, intent(out) :: ierr
    !*** local variables
    real(double), dimension(:), allocatable:: tau_grid_ref, tau_grid_ref_2nd
    integer :: ndim_ratio, ndim_ratio_half, iloc_median, ic, k, l, m,  i1, i2,&
         itau, itau_min, j, j1, j2, i, &
         ndim_int_2nd, imin_2nd, l2, k1, k2, lmin, lplus    
    real(double) :: tmin, tmax, x1, x2, y1, y2,&
         step, taua_min, taua_max, ratio_tau_median, tmin_2nd
    real(double), dimension(:), allocatable :: ratio_tau, taua_tot_sort
    real(double), dimension(:), allocatable :: taua_2nd_dum, taua_2nd_dum_sort
    real(double), dimension(:,:), allocatable :: taua_2nd_dum_arr
    integer, dimension(:), allocatable :: index_sort_dum
    real(double), dimension(ntau_2nd_max) :: tau_grid_2nd_dum
    integer, dimension(:), allocatable :: index_ratio, index_sort, index_sort_wave
    integer, dimension(1) :: imin
    integer :: ndim_inf, errorflag
    character(stringlen) :: message 
    !-----------------------------------------------------------------------

    !*** Initialize error identifier
    ierr = 0

    if(allocated(tau_grid_arr)) deallocate(tau_grid_arr)
    if(allocated(tau_grid)) deallocate(tau_grid)
    if(allocated(tau_grid_2nd)) deallocate(tau_grid_2nd)
    if(allocated(tau_grid_ref)) deallocate(tau_grid_ref)
    if(allocated(tau_grid_ref_2nd)) deallocate(tau_grid_ref_2nd)
    if(allocated(taua_tot_sort)) deallocate(taua_tot_sort)
    if(allocated(index_sort_wave)) deallocate(index_sort_wave)
    if(allocated(ntau_2nd)) deallocate(ntau_2nd)

    allocate(tau_grid_arr(ntau, ntau_2nd_max, nrt), &
         tau_grid(ntau), tau_grid_ref(ntau), &
         tau_grid_ref_2nd(ntau_2nd_max),tau_grid_2nd(ntau, ntau_2nd_max), &
         taua_tot_sort(nwave), index_sort_wave(nwave), &
         ntau_2nd(ntau), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'SET_TAU_GRID_2SPECIES: memory allocation error'
       ierr = ierr_all
       goto 999 
    endif

    tmin_2nd = minval(taua_tot_species(:, imol_2nd))
    imin_2nd = minval(minloc(taua_tot_species(:, imol_2nd)))

    ic = 0
    do k = 1, nrt
       ic = ic+1
       if(sum(taus(1:k)).gt.0.5*sum(taus(:)))exit
    enddo
    i1 = ic
    i2 = nrt
    call indexx(nwave,taua_tot_species(:,imol_1st), index_sort_wave, errorflag)
    if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6)
    forall(l=1:nwave) taua_tot_sort(l) = taua_tot_species(index_sort_wave(l),imol_1st)
    taua_max = taua_tot_sort(nwave)      
    tmax = taua_max
    tmin = minval(taua_tot_species(:,imol_1st))
    step = (DLOG(tmax)-DLOG(tmin))/dble(ntau-1)

    do l = 1, ntau
       tau_grid_ref(l) = DEXP(DLOG(tmin)+step*dble(l-1))
    enddo
    do itau = 1,ntau
       itau_min = minval(minloc(DABS(tau_grid_ref(itau)-taua_tot_species(:,imol_1st))))
       tau_grid(itau) = taua_tot_species(itau_min,imol_1st)
    enddo

    do l = 1, ntau
       lmin = max(l-2,1)
       lplus = min(l+2,ntau)
       x1 = tau_grid(l) - 1.0*(tau_grid(l)-tau_grid(lmin))
       x2 = tau_grid(l) + 1.0*(tau_grid(lplus)-tau_grid(l))       
       j1 = minval(minloc(DABS(x1-taua_tot_sort))) 
       j2 = minval(minloc(DABS(x2-taua_tot_sort))) 

       !*** now look for taua of other species corresponding to taua 1st species between x1 and x2
       ic = 0
       do j = j1, j2
          if(taua_tot_species(index_sort_wave(j),imol_2nd).gt.1.D-5)then
             ic = ic + 1
          endif
       enddo
       ndim_int_2nd = ic
       if(ndim_int_2nd.eq.0)then
          ntau_2nd(l) = 1
          tau_grid_2nd(l,1) = tmin_2nd
          tau_grid_arr(l,1,:) = taua_mol_species(imin_2nd,:,imol_2nd)
          cycle
       endif
       if(allocated(taua_2nd_dum))deallocate(taua_2nd_dum)
       if(allocated(taua_2nd_dum_arr))deallocate(taua_2nd_dum_arr)
       if(allocated(taua_2nd_dum_sort))deallocate(taua_2nd_dum_sort)
       if(allocated(index_sort_dum))deallocate(index_sort_dum)
       allocate(taua_2nd_dum(ndim_int_2nd), &
            taua_2nd_dum_arr(ndim_int_2nd,nrt), &
            taua_2nd_dum_sort(ndim_int_2nd), &
            index_sort_dum(ndim_int_2nd), &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'SET_TAU_GRID_2SPECIES: memory allocation error'
          ierr = ierr_all
          goto 999 
       endif

       ic = 0
       do j = j1, j2
          if(taua_tot_species(index_sort_wave(j),imol_2nd).gt.1.D-5)then
             ic = ic+1
             taua_2nd_dum(ic) = taua_tot_species(index_sort_wave(j),imol_2nd)
             taua_2nd_dum_arr(ic,:) = taua_mol(index_sort_wave(j),:)
          endif
       enddo
       call indexx(ndim_int_2nd,taua_2nd_dum,index_sort_dum, errorflag)
       if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6)
       forall(i=1:ndim_int_2nd)taua_2nd_dum_sort(i)=taua_2nd_dum(index_sort_dum(i))
       taua_max = maxval(taua_2nd_dum)
       taua_min = minval(taua_2nd_dum)       
       if(ndim_int_2nd.eq.1)then
          ntau_2nd(l) = 1
          tau_grid_2nd(l,1) = taua_2nd_dum(1)
          tau_grid_arr(l,1,:) = taua_2nd_dum_arr(1,:)
          cycle
       endif
       tmax = taua_max
       tmin = taua_min         
       step = (DLOG(tmax)-DLOG(tmin))/dble(ntau_2nd_max-1)
       do itau=1,ntau_2nd_max
          tau_grid_ref_2nd(itau) = DEXP(DLOG(tmin)+step*dble(itau-1))
       enddo
       do itau = 1,ntau_2nd_max
          itau_min = minval(minloc(DABS(tau_grid_ref_2nd(itau)-taua_2nd_dum)))
          tau_grid_2nd_dum(itau) =taua_2nd_dum(itau_min) 
       enddo

       !*** check for double values in tau_grid
       tau_grid_2nd(l,1) = tau_grid_2nd_dum(1)
       ic = 1
       do itau = 2, ntau_2nd_max
          ! two values are considered the same if they differ less than 1.D-5
          if(dabs(tau_grid_2nd_dum(itau)-tau_grid_2nd(l,ic)).gt.1.D-5)then
             ic = ic+1
             tau_grid_2nd(l,ic) = tau_grid_2nd_dum(itau)
          endif
       enddo
       ntau_2nd(l) = ic
       imin = maxloc(taua_2nd_dum)
       tau_grid_arr(l,ntau_2nd(l),:) = taua_2nd_dum_arr(imin(1),:)           
       imin = minloc(taua_2nd_dum)
       tau_grid_arr(l,1,:) = taua_2nd_dum_arr(imin(1),:)
       call indexx(ndim_int_2nd,taua_2nd_dum,index_sort_dum, errorflag)
       if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6) 
       forall(i=1:ndim_int_2nd)taua_2nd_dum_sort(i)=taua_2nd_dum(index_sort_dum(i))

       do l2 = 1,ntau_2nd(l)        
          lmin = max(l2-1,1)
          lplus = min(l2+1,ntau_2nd(l))
          y1 = tau_grid_2nd(l,l2)- 0.49*(tau_grid_2nd(l,l2)-tau_grid_2nd(l,lmin))
          y2 = tau_grid_2nd(l,l2) + 0.49*(tau_grid_2nd(l,lplus)-tau_grid_2nd(l,l2))          
          k1 = minval(minloc(DABS(y1-taua_2nd_dum_sort))) 
          k2 = minval(minloc(DABS(y2-taua_2nd_dum_sort))) 
          ndim_ratio = max(k2-k1+1,1)
          allocate(ratio_tau(ndim_ratio),index_ratio(ndim_ratio),index_sort(ndim_ratio), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'SET_TAU_GRID_2SPECIES: memory allocation error'
             ierr = ierr_all
             goto 999 
          endif

          ic = 0
          do k = k1, k2
             ic = ic+1
             m=index_sort_dum(k)
             if (i1 > 1) then !!HH: avoid segmentation fault for i1=1
                ratio_tau(ic) = sum(taua_2nd_dum_arr(m,1:i1-1))/sum(taua_2nd_dum_arr(m,i1:i2))
             else
                ratio_tau(ic) = 0.d0
             endif
             index_ratio(ic) = m
          enddo
          call indexx(ndim_ratio,ratio_tau,index_sort, errorflag)
          if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6)
          ndim_inf = size(pack(ratio_tau,ratio_tau.gt.1.D5))
          ndim_ratio_half = max((ndim_ratio-ndim_inf)/2,1)
          ratio_tau_median = ratio_tau(index_sort(ndim_ratio_half))                   
          iloc_median = index_ratio(index_sort(ndim_ratio_half))          
          tau_grid_arr(l,l2,:) = (tau_grid(l)+tau_grid_2nd(l,l2))/sum(taua_2nd_dum_arr(iloc_median,:)) * &
               taua_2nd_dum_arr(iloc_median,:)
          deallocate(ratio_tau, index_ratio, index_sort, stat=ierr)   
          if (ierr .ne. 0) then
             write(message, *) 'SET_TAU_GRID_2SPECIES: memory deallocation error'
             ierr = ierr_deall
             goto 999 
          endif
       enddo
       deallocate(&
            taua_2nd_dum,&
            taua_2nd_dum_arr,&
            taua_2nd_dum_sort,&
            index_sort_dum, &
            stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'SET_TAU_GRID_2SPECIES: memory deallocation error'
          ierr = ierr_deall
          goto 999 
       endif
    enddo

    deallocate(taua_tot_sort, index_sort_wave, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'SET_TAU_GRID_2SPECIES: memory deallocation error'
       ierr = ierr_deall
       goto 999 
    endif

    return

999 continue
    call stopretrieval(message)

  end subroutine set_tau_grid_2species

!***************************************************************************************

  subroutine set_tau_grid_1species(&
       nwave, & 
       nrt, &
       ntau, & 
       taus, &
       taua_tot, &
       taua_mol, &
       tau_grid, &
       tau_grid_arr, & 
       ierr)
    !*** Input
    integer, intent(in) :: nwave, ntau, nrt
    real(double), dimension(:), intent(in) :: taus 
    real(double), dimension(:), intent(in) :: taua_tot
    real(double), dimension(:,:), intent(in) :: taua_mol 
    !*** Output 
    real(double), dimension(:), allocatable, intent(out) :: tau_grid
    real(double), dimension(:,:,:), allocatable, intent(out) :: tau_grid_arr
    integer, intent(out) :: ierr
    !*** local variables 
    real(double), dimension(:), allocatable:: tau_grid_ref
    integer :: ndim_ratio, ndim_ratio_half, iloc_median, ic, k, l, i1, i2,&
         itau, itau_min, j, j1, j2, lmin, lplus, ntau_2nd_max    
    real(double) :: tmin, tmax, x1, x2, &
         step, taua_max, ratio_tau_median
    real(double), dimension(:), allocatable :: ratio_tau, taua_tot_sort
    integer, dimension(:), allocatable :: index_ratio, index_sort, index_sort_wave
    integer :: ndim_inf, errorflag
    character(stringlen) :: message
    !-----------------------------------------------------------------------

    !*** Initialize error identifier
    ierr = 0

    if(allocated(tau_grid_arr)) deallocate(tau_grid_arr)
    if(allocated(tau_grid)) deallocate(tau_grid)
    if(allocated(tau_grid_ref)) deallocate(tau_grid_ref)
    if(allocated(taua_tot_sort)) deallocate(taua_tot_sort)
    if(allocated(index_sort_wave)) deallocate(index_sort_wave)

    ntau_2nd_max = 1

    allocate(tau_grid_arr(ntau, ntau_2nd_max, nrt), &
         tau_grid(ntau), tau_grid_ref(ntau), &
         taua_tot_sort(nwave), index_sort_wave(nwave), stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'SET_TAU_GRID_1SPECIES: memory allocation error'
       ierr = ierr_all
       goto 999 
    endif

    ic = 0
    do k = 1, nrt
       if(sum(taus(1:k)).gt.0.5*sum(taus(:))) exit
       ic = ic + 1
    enddo

    i1 = ic
    i2 = nrt

    call indexx(nwave, taua_tot, index_sort_wave, errorflag)
    if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6)
    forall(l=1:nwave) taua_tot_sort(l) = taua_tot(index_sort_wave(l))

    taua_max = taua_tot_sort(nwave)     
    tmax = taua_max
!    tmin = minval(taua_tot)  !HH: Error for very low or negative taua_tot)
    tmin = max(minval(taua_tot), 1.d-5)

    step = (DLOG(tmax)-DLOG(tmin))/dble(ntau-1)
    do l = 1, ntau
       tau_grid_ref(l) = DEXP(DLOG(tmin) + step*dble(l-1))
    enddo

    do itau = 1,ntau
       itau_min = minval(minloc(DABS(tau_grid_ref(itau)-taua_tot(:))))
       tau_grid(itau) = taua_tot(itau_min)
    enddo

    do l = 1, ntau
       lmin = max(l-1,1)
       lplus = min(l+1,ntau)
       x1 = tau_grid(l) - 0.49*(tau_grid(l)-tau_grid(lmin))
       x2 = tau_grid(l) + 0.49*(tau_grid(lplus)-tau_grid(l))       
       j1 = minval(minloc(DABS(x1-taua_tot_sort))) 
       j2 = minval(minloc(DABS(x2-taua_tot_sort))) 
       ndim_ratio = max(j2-j1+1, 1)
       allocate(ratio_tau(ndim_ratio),index_ratio(ndim_ratio),index_sort(ndim_ratio), stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'SET_TAU_GRID_1SPECIES: memory allocation error'
          ierr = ierr_all
          goto 999 
       endif
       ic = 0
       do j = j1, j2
          ic = ic + 1
          k = index_sort_wave(j)
          if (i1 > 1) then !HH: avoid segmentation fault for i1=1
             ratio_tau(ic) = sum(taua_mol(k,1:i1-1))/sum(taua_mol(k,i1:i2))
          else
             ratio_tau(ic) = 0.d0
          endif
          index_ratio(ic) = k
       enddo
       call indexx(ndim_ratio, ratio_tau, index_sort, errorflag)
       if (errorflag .ne. 0) call writelog('INDEXX: NSTACK too small', 6)
       ndim_inf = size(pack(ratio_tau,ratio_tau.gt.1.D5))
       ndim_ratio_half = max((ndim_ratio-ndim_inf)/2,1)
       ratio_tau_median = ratio_tau(index_sort(ndim_ratio_half))
       iloc_median = index_ratio(index_sort(ndim_ratio_half))       
       tau_grid_arr(l,1,:) = tau_grid(l)/sum(taua_mol(iloc_median,:)) * &
            taua_mol(iloc_median,:)     
       deallocate(ratio_tau, index_ratio, index_sort, stat=ierr)
       if (ierr .ne. 0) then
          write(message, *) 'SET_TAU_GRID_1SPECIES: memory deallocation error'
          ierr = ierr_deall
          goto 999 
       endif
    enddo

    deallocate(taua_tot_sort, index_sort_wave, stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'SET_TAU_GRID_1SPECIES: memory deallocation error'
       ierr = ierr_deall
       goto 999 
    endif


    return

999 continue
    call stopretrieval(message)


  end subroutine set_tau_grid_1species
!****************************************************************************************
!*** Interpolate from tau-grid to wavelength grid for two absorbers 
  subroutine interpolate_lbl_2species(&
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
       rint_fine_ms,&
       dtaua_tau,&
       kmat_taua_fine,&
!       kmat_aerosol_tau,&
       kmat_aerosol_phase,&
       kmat_aerosol_taua,&  
       kmat_aerosol_taus,& 
       kmat_alb_tau,&
       kmat_taus_tau,&
       kmat_taus_fine, &
!       daerosol, &
       daerosol_phase, &
       daerosol_taua, &
       daerosol_taus, &
       dalb) 
   use aerosol_input_module, only: aero
   use auxiliary_routines_module
!*** Input
   integer, intent(in) :: nwave, nrt
   type(aero), dimension(:), intent(in) :: aerosol 
   integer, dimension(:), intent(in) :: ntau_2nd 
   real(double), dimension(:), intent(in) :: tau_grid
   real(double), dimension(:,:), intent(in) :: tau_grid_2nd
   real(double), dimension(:,:,:), intent(in) :: tau_grid_arr
   integer, intent(in) :: ntau, imol_1st, imol_2nd    
   real(double), dimension(:,:), intent(in) :: taua_tot_species
   real(double), dimension(:), intent(in) :: taua_tot
   real(double), dimension(:,:), intent(in) :: taua_mol
   real(double), dimension(:,:,:), intent(in) :: rint_tau
   real(double), dimension(:,:,:,:), intent(in) :: dtaua_tau, kmat_taus_tau
!   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_tau
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_phase
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_taua
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_taus 
   real(double), dimension(:,:,:), intent(in) :: kmat_alb_tau
!*** output  
!   real(double), dimension(:,:,:,:), intent(out) :: daerosol

   real(double), dimension(:,:,:,:), intent(out) :: daerosol_phase !<=========================== D.S
   real(double), dimension(:,:,:,:), intent(out) :: daerosol_taua  !<=========================== D.S
   real(double), dimension(:,:,:,:), intent(out) :: daerosol_taus  !<=========================== D.S

   real(double), dimension(:,:), intent(out) :: dalb  
   real(double), dimension(:,:), intent(out) :: rint_fine_ms
   real(double), dimension(:,:,:), intent(out) :: kmat_taus_fine, kmat_taua_fine
!*** Local variables  
   integer :: i0, i1, i2, l, i, k, ic, i_st, ntype_aer
   integer, dimension(1) :: imin
   integer, parameter :: npow = 3
   real(double) :: rint_i0, rint_i1, rint_i2
   integer :: ipow
   real(double), dimension(3,npow) :: kmat_poly
   real(double), dimension(3) :: y_poly, x_poly  
!---------------------------------------------------------------------  
      ntype_aer = size(aerosol)

!*** Loop over all wavelength  
      do l = 1, nwave
!*** Polynomial interpolation for small and intermediate optical densities:
!*** Look for 2nd order polynomial in taua_tot_1st_species with polynomial coefficents x
!*** that describes best the modelled quantity (reflectance,derivatives) at three adjacent points y,
!*** i.e. y0 = x0 + x1*tau + x2*tau^2, y1 = ..., y2 = ...
!*** i.e. y=K*x where K is dyj/dxi=tau_j^i.
!*** Then use polynomial to calculate reflectance at intermediate wavelength.	  
  
!*** Choose grid points from tau_grid that are closest to taua_tot at the considered wavelength
         imin = minloc(DABS(taua_tot_species(l,imol_1st) - tau_grid(:)))
         ic = max(imin(1), 2)         
         i1 = min(ic, ntau-1)  !** center point
         i0 = i1 - 1           !** center point - 1
         i2 = i1 + 1           !** center point + 1
!*** Construct K = tau_j^i
         forall(ipow=1:npow) kmat_poly(1,ipow) = tau_grid(i0)**dble(ipow-1)
         forall(ipow=1:npow) kmat_poly(2,ipow) = tau_grid(i1)**dble(ipow-1)
         forall(ipow=1:npow) kmat_poly(3,ipow) = tau_grid(i2)**dble(ipow-1)

!*** now first interpolate to the correct taua_tot_species(l,imol_2nd) of the second
!*** absorber for the 3 grid points i0, i1, and i2 corresponding to the grid
!*** of the 1st absorber 
         do i_st = 1, nstokes_der
!*** Interpolate reflectance
            call interpolate_tau_2nd_pt(&
               nrt, &
               ntau_2nd, &
               tau_grid, &
               tau_grid_2nd, &
               tau_grid_arr, & 
               taua_tot_species(l, imol_2nd),&
               taua_tot(l),&
               taua_mol(l, :),&
               i0,&
               rint_tau(i0, :, i_st),&
               dtaua_tau(:, :, :, i_st),&
               rint_i0)
            call interpolate_tau_2nd_pt(&
               nrt, &
               ntau_2nd, &
               tau_grid, &
               tau_grid_2nd, &
               tau_grid_arr, & 
               taua_tot_species(l, imol_2nd),&
               taua_tot(l),&
               taua_mol(l, :),&
               i1,&
               rint_tau(i1, :, i_st),&
               dtaua_tau(:, :, :, i_st),&
               rint_i1)          
            call interpolate_tau_2nd_pt(&
               nrt, &
               ntau_2nd, &
               tau_grid, &
               tau_grid_2nd, &
               tau_grid_arr, & 
               taua_tot_species(l, imol_2nd),&
               taua_tot(l),&
               taua_mol(l, :),&
               i2,&
               rint_tau(i2, :, i_st),&
               dtaua_tau(:, :, :, i_st),&
               rint_i2)

!*** Interpolate by 2nd order polynomial in taua_tot_1st_species 
!*** to intermediate wavelength, making use of rint_i0, rint_i1, and rint_i2  
            call interpolate_tau_1st( &
               tau_grid, &
               i0, rint_i0, &
               i1, rint_i1, &
               i2, rint_i2, &
               taua_tot_species(l, imol_1st), &
               rint_fine_ms(l, i_st))
          
!*** Interpolate derivatives wrt aerosol parameters         
            do i = 1, ntype_aer             
               do k = 1, npar
                  if(aerosol(i)%AerosolFlags(k).ne.0) then                
!!$!** Linearly interpolate in taua_tot_2nd_species to i0,i1,i2  
!!$!*** Interpolate K-matrix containing SUMMED contribs stemming from dI/dphase * dphase/dx_aer, dI/taua * dtaua/dx_aer and dI/dtaus * dtaus/dx_aer 
!!$                     call interpolate_tau_2nd(&
!!$                        ntau_2nd, &
!!$                        tau_grid_2nd, &
!!$                        taua_tot_species(l, imol_2nd), &
!!$                        i0, &
!!$                        kmat_aerosol_tau(i, i0, :, k, i_st), & 
!!$                        y_poly(1))
!!$                     call interpolate_tau_2nd(&
!!$                        ntau_2nd, &
!!$                        tau_grid_2nd, &
!!$                        taua_tot_species(l, imol_2nd), &
!!$                        i1, &
!!$                        kmat_aerosol_tau(i, i1, :, k, i_st), &
!!$                        y_poly(2))
!!$                     call interpolate_tau_2nd(&
!!$                        ntau_2nd, &
!!$                        tau_grid_2nd, &
!!$                        taua_tot_species(l, imol_2nd), &
!!$                        i2, &
!!$                        kmat_aerosol_tau(i, i2, :, k, i_st), &
!!$                        y_poly(3))                
!!$!** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                     
!!$                     x_poly = &
!!$                     matmul(&
!!$                     INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
!!$                     matmul(transpose(kmat_poly),y_poly))                
!!$!** Interpolate by 2nd order polynomial in taua_tot_1st_species                     
!!$                     daerosol(l, k, i, i_st) = &
!!$                     x_poly(1) + &
!!$                     x_poly(2)*taua_tot_species(l, imol_1st) + &
!!$                     x_poly(3)*taua_tot_species(l, imol_1st)**2            

!*** Interpolate K-matrix containing ONLY contribs stemming from dI/dphase * dphase/dx_aer, 
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i0, &
                        kmat_aerosol_phase(i, i0, :, k, i_st), & 
                        y_poly(1))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i1, &
                        kmat_aerosol_phase(i, i1, :, k, i_st), &
                        y_poly(2))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i2, &
                        kmat_aerosol_phase(i, i2, :, k, i_st), &
                        y_poly(3))                
!** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                     
                     x_poly = &
                     matmul(&
                     INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                     matmul(transpose(kmat_poly),y_poly))                
!** Interpolate by 2nd order polynomial in taua_tot_1st_species                     
                     daerosol_phase(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot_species(l, imol_1st) + &
                     x_poly(3)*taua_tot_species(l, imol_1st)**2


!*** Interpolate K-matrix containing ONLY contribs stemming from  dI/taua * dtaua/dx_aer 
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i0, &
                        kmat_aerosol_taua(i, i0, :, k, i_st), & 
                        y_poly(1))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i1, &
                        kmat_aerosol_taua(i, i1, :, k, i_st), &
                        y_poly(2))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i2, &
                        kmat_aerosol_taua(i, i2, :, k, i_st), &
                        y_poly(3))                
!** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                     
                     x_poly = &
                     matmul(&
                     INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                     matmul(transpose(kmat_poly),y_poly))                
!** Interpolate by 2nd order polynomial in taua_tot_1st_species                     
                     daerosol_taua(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot_species(l, imol_1st) + &
                     x_poly(3)*taua_tot_species(l, imol_1st)**2


!*** Interpolate K-matrix containing ONLY contribs stemming from  dI/dtaus * dtaus/dx_aer 
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i0, &
                        kmat_aerosol_taus(i, i0, :, k, i_st), & 
                        y_poly(1))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i1, &
                        kmat_aerosol_taus(i, i1, :, k, i_st), &
                        y_poly(2))
                     call interpolate_tau_2nd(&
                        ntau_2nd, &
                        tau_grid_2nd, &
                        taua_tot_species(l, imol_2nd), &
                        i2, &
                        kmat_aerosol_taus(i, i2, :, k, i_st), &
                        y_poly(3))                
!** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                     
                     x_poly = &
                     matmul(&
                     INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                     matmul(transpose(kmat_poly),y_poly))                
!** Interpolate by 2nd order polynomial in taua_tot_1st_species                     
                     daerosol_taus(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot_species(l, imol_1st) + &
                     x_poly(3)*taua_tot_species(l, imol_1st)**2

                  endif
               enddo		   
            enddo ! end loop over aerosol types
  
!*** Interpolate derivatives wrt absorption optical density
            do k = 1, nrt             
!*** Linearly interpolate in taua_tot_2nd_species to i0,i1,i2  
               call interpolate_tau_2nd(&
                    ntau_2nd, &
                    tau_grid_2nd, &
                    taua_tot_species(l,imol_2nd), &
                    i0, dtaua_tau(i0, :, k, i_st), &
                    y_poly(1))
               call interpolate_tau_2nd(&
                     ntau_2nd, &
                     tau_grid_2nd, &
                     taua_tot_species(l,imol_2nd), &
                     i1, &
                     dtaua_tau(i1, :, k, i_st), &
                     y_poly(2))
               call interpolate_tau_2nd(&
                    ntau_2nd, &
                    tau_grid_2nd,  &
                    taua_tot_species(l,imol_2nd), &
                    i2, &
                    dtaua_tau(i2, :, k, i_st), &
                    y_poly(3))  
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                          
               x_poly = &
                  matmul(&
                  INVERSE_3BY3(matmul(transpose(kmat_poly), kmat_poly)),&
                  matmul(transpose(kmat_poly), y_poly))             
!*** Interpolate by 2nd order polynomial in taua_tot_1st_species       
             
               kmat_taua_fine(l, k, i_st) = &
                  x_poly(1) + &
                  x_poly(2)*taua_tot_species(l,imol_1st) + &
                  x_poly(3)*taua_tot_species(l,imol_1st)**2 
            enddo
          
!*** Interpolate derivatives wrt albedo 
!*** Linearly interpolate in taua_tot_2nd_species to i0,i1,i2  
            call interpolate_tau_2nd(&
                 ntau_2nd, &
                 tau_grid_2nd, &
                 taua_tot_species(l,imol_2nd), &
                 i0, &
                 kmat_alb_tau(i0,:,i_st), &
                 y_poly(1))
            call interpolate_tau_2nd(&
                 ntau_2nd, &
                 tau_grid_2nd, &
                 taua_tot_species(l,imol_2nd), &
                 i1, &
                 kmat_alb_tau(i1,:,i_st), &
                 y_poly(2))
            call interpolate_tau_2nd(&
                 ntau_2nd, &
                 tau_grid_2nd, &
                 taua_tot_species(l,imol_2nd), &
                 i2, &
                 kmat_alb_tau(i2,:,i_st), &
                 y_poly(3))          
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                          
            x_poly = &
               matmul(&
               INVERSE_3BY3(matmul(transpose(kmat_poly), kmat_poly)),&
               matmul(transpose(kmat_poly), y_poly)) 
!*** Interpolate by 2nd order polynomial in taua_tot_1st_species                   
            dalb(l,i_st) = &
               x_poly(1) + &
               x_poly(2)*taua_tot_species(l,imol_1st) + &
               x_poly(3)*taua_tot_species(l,imol_1st)**2 

!*** Interpolate derivatives wrt taus  
            do k = 1, nrt           
!*** Linearly interpolate in taua_tot_2nd_species to i0,i1,i2  
               call interpolate_tau_2nd(&
                    ntau_2nd, &
                    tau_grid_2nd, &
                    taua_tot_species(l, imol_2nd), &
                    i0, &
                    kmat_taus_tau(i0, :, k, i_st), &
                    y_poly(1))
               call interpolate_tau_2nd(&
                    ntau_2nd, &
                    tau_grid_2nd, &
                    taua_tot_species(l, imol_2nd), &
                    i1, &
                    kmat_taus_tau(i1, :, k, i_st), &
                    y_poly(2))
               call interpolate_tau_2nd(&
                    ntau_2nd, &
                    tau_grid_2nd,   taua_tot_species(l, imol_2nd), &
                    i2, &
                    kmat_taus_tau(i2, :, k, i_st), &
                    y_poly(3))             
!** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                          
               x_poly = &
                  matmul(&
                  INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                  matmul(transpose(kmat_poly),y_poly))         
!*** Interpolate by 2nd order polynomial in taua_tot_1st_species                   
               kmat_taus_fine(l, k, i_st) = &
                  x_poly(1) + &
                  x_poly(2)*taua_tot_species(l,imol_1st) + &
                  x_poly(3)*taua_tot_species(l,imol_1st)**2 
             
            enddo          
         enddo
      enddo

!**** End loop over all wavelengths
  
  end subroutine interpolate_lbl_2species

!*****************************************************************
!*** Interpolate from tau-grid to wavelength grid for one absorber   
   subroutine interpolate_lbl_1species(&
       nwave, &
       nrt, &
       tau_grid, &
       tau_grid_arr, &
       aerosol, &
       ntau,&
       taua_tot,&
       taua_mol,&
       rint_tau,&
       rint_fine_ms,&
       dtaua_tau,&
       kmat_taua_fine,&
!       kmat_aerosol_tau,&
       kmat_aerosol_phase,& 
       kmat_aerosol_taua,& 
       kmat_aerosol_taus,&  
       kmat_alb_tau,&
       kmat_taus_tau,&
       kmat_taus_fine, &
!       daerosol, &
       daerosol_phase, &
       daerosol_taua, &
       daerosol_taus, &
       dalb)
   use aerosol_input_module, only: aero
   use auxiliary_routines_module
!*** Input
   integer, intent(in) :: nwave, nrt
   type(aero), dimension(:), intent(in) :: aerosol 
   real(double), dimension(:), intent(in) :: tau_grid
   real(double), dimension(:,:,:), intent(in) :: tau_grid_arr
   integer, intent(in) :: ntau  
   real(double), dimension(:), intent(in) :: taua_tot
   real(double), dimension(:,:), intent(in) :: taua_mol
   real(double), dimension(:,:,:), intent(in) :: rint_tau
   real(double), dimension(:,:,:,:), intent(in) :: dtaua_tau, kmat_taus_tau
!   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_tau
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_phase
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_taua 
   real(double), dimension(:,:,:,:,:), intent(in) :: kmat_aerosol_taus
   real(double), dimension(:,:,:), intent(in) :: kmat_alb_tau
!*** Output
!   real(double), dimension(:,:,:,:), intent(out) :: daerosol
   real(double), dimension(:,:,:,:), intent(out) :: daerosol_phase !<=========================== D.S
   real(double), dimension(:,:,:,:), intent(out) :: daerosol_taua  !<=========================== D.S
   real(double), dimension(:,:,:,:), intent(out) :: daerosol_taus  !<=========================== D.S
   real(double), dimension(:,:), intent(out) :: dalb
!*** Output
   real(double), dimension(:,:,:), intent(out) :: kmat_taus_fine, kmat_taua_fine
   real(double), dimension(:,:), intent(out) :: rint_fine_ms
!*** Local variables 
    real(double), dimension(:),allocatable :: taua_mol_scale
    integer :: i0, i1, i2, l, i, k, ic, i_st, ntype_aer  
    integer, dimension(1) :: imin
    integer, parameter :: npow = 3  
    real(double) :: scale
    integer :: ipow
    real(double), dimension(3,npow) :: kmat_poly
    real(double), dimension(3) :: y_poly, x_poly
!---------------------------------------------------------------  
      if(allocated(taua_mol_scale)) deallocate(taua_mol_scale)
      allocate(taua_mol_scale(nrt))
    
     ntype_aer = size(aerosol)
!*** Loop over all wavelengths
      do l = 1, nwave

!*** Polynomial interpolation for small and intermediate optical densities:  
!*** Look for 2nd order polynomial in taua_tot with polynomial coefficents x
!*** that describes best the modelled quantity (reflectance,derivatives) at three adjacent points y,
!*** i.e. y0 = x0 + x1*tau + x2*tau^2, y1 = ..., y2 = ...
!*** i.e. y=K*x where K is dyj/dxi=tau_j^i.
!*** Then use polynomial to calculate reflectance at intermediate wavelength.	  
	  
!*** Choose grid points from tau_grid that are closest to taua_tot at the considered wavelength
         imin = minloc(DABS(taua_tot(l)-tau_grid))
         ic = max(imin(1),2)        
         i1 = min(ic,ntau-1) !** center point
         i0 = i1-1           !** center point - 1
         i2 = i1+1           !** center point + 1
          
!*** Construct K = tau_j^i
         forall(ipow=1:npow) kmat_poly(1,ipow) = tau_grid(i0)**dble(ipow-1)
         forall(ipow=1:npow) kmat_poly(2,ipow) = tau_grid(i1)**dble(ipow-1)
         forall(ipow=1:npow) kmat_poly(3,ipow) = tau_grid(i2)**dble(ipow-1)
	  	  
         do i_st = 1, nstokes_der
!** Interpolate reflectance
	  
!*** Construct y = reflectance, linearly corrected for different absorption optical depth at considered wavelength
!*** compared to binned taua_tot         
            scale = tau_grid(i0)/taua_tot(l)
            do k =1, nrt
               taua_mol_scale(k) = scale*taua_mol(l,k)
            enddo
          
            y_poly(1) = rint_tau(i0,1,i_st) + &
                     sum(dtaua_tau(i0, 1, :, i_st) *(taua_mol_scale-tau_grid_arr(i0,1,:))) 
          
            scale = tau_grid(i1)/taua_tot(l)
            do k =1 ,nrt
               taua_mol_scale(k) = scale*taua_mol(l,k)
            enddo
          
            y_poly(2) = rint_tau(i1,1,i_st) + &
                     sum(dtaua_tau(i1, 1, :, i_st) *(taua_mol_scale-tau_grid_arr(i1,1,:)))
             
            scale = tau_grid(i2)/taua_tot(l)
            do k =1 ,nrt
               taua_mol_scale(k) = scale*taua_mol(l,k)
            enddo
             
            y_poly(3) = rint_tau(i2,1,i_st) + &
                     sum(dtaua_tau(i2, 1, :, i_st) *(taua_mol_scale-tau_grid_arr(i2,1,:)))
                   
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                
            x_poly = matmul(&
                  INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                  matmul(transpose(kmat_poly),y_poly))
          
!** Interpolate by 2nd order polynomial in taua_tot                    
            rint_fine_ms(l,i_st) = &
               x_poly(1) + &
               x_poly(2)*taua_tot(l) + &
               x_poly(3)*taua_tot(l)**2 
 
!*** Interpolate derivatives wrt aerosol parameters             
            do i = 1, ntype_aer           
               do k = 1, npar 
                  if(aerosol(i)%AerosolFlags(k).ne.0)then
!!$!*** Construct y = aerosol derivatives		
!!$!*** Interpolate K-matrix containing SUMMED contribs stemming from dI/dphase * dphase/dx_aer, dI/taua * dtaua/dx_aer and dI/dtaus * dtaus/dx_aer 
!!$                     y_poly(1) = kmat_aerosol_tau(i, i0, 1, k, i_st)
!!$                     y_poly(2) = kmat_aerosol_tau(i, i1, 1, k, i_st)
!!$                     y_poly(3) = kmat_aerosol_tau(i, i2, 1, k, i_st)
!!$                  
!!$!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                    
!!$                     x_poly = matmul(&
!!$                           INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
!!$                           matmul(transpose(kmat_poly),y_poly))
!!$
!!$!*** Interpolate by 2nd order polynomial in taua_tot                   
!!$                     daerosol(l, k, i, i_st) = &
!!$                     x_poly(1) + &
!!$                     x_poly(2)*taua_tot(l) + &
!!$                     x_poly(3)*taua_tot(l)**2 		

!*** Interpolate K-matrix containing ONLY contribs stemming from dI/dphase * dphase/dx_aer
                     y_poly(1) = kmat_aerosol_phase(i, i0, 1, k, i_st)
                     y_poly(2) = kmat_aerosol_phase(i, i1, 1, k, i_st)
                     y_poly(3) = kmat_aerosol_phase(i, i2, 1, k, i_st)
                  
		    !*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                    
                     x_poly = matmul(&
                           INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                           matmul(transpose(kmat_poly),y_poly))

		    !*** Interpolate by 2nd order polynomial in taua_tot                   
                     daerosol_phase(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot(l) + &
                     x_poly(3)*taua_tot(l)**2 	

!*** Interpolate K-matrix containing ONLY contribs stemming from dI/taua * dtaua/dx_aer  
                     y_poly(1) = kmat_aerosol_taua(i, i0, 1, k, i_st)
                     y_poly(2) = kmat_aerosol_taua(i, i1, 1, k, i_st)
                     y_poly(3) = kmat_aerosol_taua(i, i2, 1, k, i_st)
                  
		    !*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                    
                     x_poly = matmul(&
                           INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                           matmul(transpose(kmat_poly),y_poly))

		    !*** Interpolate by 2nd order polynomial in taua_tot                   
                     daerosol_taua(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot(l) + &
                     x_poly(3)*taua_tot(l)**2 	


!*** Interpolate K-matrix containing ONLY contribs stemming from  dI/dtaus * dtaus/dx_aer 
                     y_poly(1) = kmat_aerosol_taus(i, i0, 1, k, i_st)
                     y_poly(2) = kmat_aerosol_taus(i, i1, 1, k, i_st)
                     y_poly(3) = kmat_aerosol_taus(i, i2, 1, k, i_st)
                  
	    !*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                    
                     x_poly = matmul(&
                           INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                           matmul(transpose(kmat_poly),y_poly))

	    !*** Interpolate by 2nd order polynomial in taua_tot                   
                     daerosol_taus(l, k, i, i_st) = &
                     x_poly(1) + &
                     x_poly(2)*taua_tot(l) + &
                     x_poly(3)*taua_tot(l)**2 	


                  endif
               enddo		
            enddo
             
!*** Interpolate derivatives wrt absorption optical density  
            do k=1,nrt                
!*** Construct y = absorption optical density derivatives		
               y_poly(1) = dtaua_tau(i0, 1, k, i_st)
               y_poly(2) = dtaua_tau(i1, 1, k, i_st)
               y_poly(3) = dtaua_tau(i2, 1, k, i_st)
             
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                       
               x_poly = &
                  matmul(&
                  INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                  matmul(transpose(kmat_poly),y_poly))
                
!*** Interpolate by 2nd order polynomial in taua_tot                       
               kmat_taua_fine(l, k, i_st) = &
                  x_poly(1) + &
                  x_poly(2)*taua_tot(l) + &
                  x_poly(3)*taua_tot(l)**2 
            enddo
  
!*** Interpolate derivatives wrt albedo
!*** Construct y = albedo derivatives		
            y_poly(1) = kmat_alb_tau(i0,1,i_st)
            y_poly(2) = kmat_alb_tau(i1,1,i_st)
            y_poly(3) = kmat_alb_tau(i2,1,i_st)
  
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                       
            x_poly = &
               matmul(&
               INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
               matmul(transpose(kmat_poly),y_poly))
             
!*** Interpolate by 2nd order polynomial in taua_tot                      
            dalb(l, i_st) = &
               x_poly(1) + &
               x_poly(2)*taua_tot(l) + &
               x_poly(3)*taua_tot(l)**2 
  
!*** Interpolate derivatives wrt scattering optical density             
            do k = 1, nrt
                
!*** Construct y = scattering optical density derivatives            
               y_poly(1) = kmat_taus_tau(i0, 1, k, i_st)
               y_poly(2) = kmat_taus_tau(i1, 1, k, i_st)
               y_poly(3) = kmat_taus_tau(i2, 1, k, i_st)
             
!*** Least squares solution of y=K*x: x=(KT*K)^-1*KT*y                       
               x_poly = &
                  matmul(&
                  INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
                  matmul(transpose(kmat_poly),y_poly))
  
!*** Interpolate by 2nd order polynomial in taua_tot                       
               kmat_taus_fine(l, k, i_st) = &
                  x_poly(1) + &
                  x_poly(2)*taua_tot(l) + &
                  x_poly(3)*taua_tot(l)**2 
            enddo
         enddo
      enddo
!*** End loop over all wavelengths
    
      deallocate(taua_mol_scale)

   end subroutine interpolate_lbl_1species

!************************************************************************************  
!*** Linear interpolation on grid of 2nd species without correction for vertical profile 
!*** (used for Jacobian) 
   subroutine interpolate_tau_2nd(ntau_2nd, tau_grid_2nd, taua_tot_2nd, ix, rint_tau, rint_ix)
   integer, dimension(:), intent(IN) :: ntau_2nd  
   real(double), dimension(:,:), intent(in) :: tau_grid_2nd
   real(double), intent(in) :: taua_tot_2nd  
   integer, intent(in) :: ix
   real(double), dimension(:), intent(in) :: rint_tau  
   real(double), intent(out) :: rint_ix
!*** local variables  
   integer :: j0, j1  
   real(double) :: rint_j0, rint_j1 
!-----------------------------------------------------------------------------------
      if(ntau_2nd(ix).eq.1) then
         j0 = 1
         rint_ix = rint_tau(j0)
         goto 99
      endif  
      j0 = minval(minloc(DABS(taua_tot_2nd-tau_grid_2nd(ix,:)),&
           taua_tot_2nd.ge.tau_grid_2nd(ix,:)))    
      j0 = min(j0, ntau_2nd(ix)-1)
      j0 = max(j0, 1)
      j1 = j0 + 1   
      rint_j0 = rint_tau(j0) 
      rint_j1 = rint_tau(j1)    
      rint_ix = rint_j0 + (rint_j1-rint_j0)/ &
         (tau_grid_2nd(ix,j1)-tau_grid_2nd(ix,j0)) * &
         (taua_tot_2nd-tau_grid_2nd(ix,j0))
  
  99  continue  
  
   end subroutine interpolate_tau_2nd

!************************************************************************************************* 
!*** Interpolate to grid of 2nd absorber, for 1 point on the 1st absorber grid. 
!*** Correct for vertical profile in linear approximation
!*** (Hasekamp and Butz, JGR, 2008) 
   subroutine interpolate_tau_2nd_pt(&
       nrt, &
       ntau_2nd, &
       tau_grid, &
       tau_grid_2nd, &
       tau_grid_arr, &
       taua_tot_2nd,&
       taua_tot_all,&
       taua_arr,&
       ix,&
       rint_tau,&
       dtaua_tau,&
       rint_ix)
    use auxiliary_routines_module
!*** Input
   integer, dimension(:), intent(in) :: ntau_2nd  
   real(double), dimension(:), intent(in) :: tau_grid
   real(double), dimension(:,:), intent(in) :: tau_grid_2nd
   real(double), dimension(:,:,:), intent(in) :: tau_grid_arr
   real(double), intent(in) :: taua_tot_2nd, taua_tot_all
   real(double), dimension(:), intent(in) :: taua_arr
   real(double), dimension(:) , intent(in):: rint_tau
   real(double), dimension(:,:,:), intent(in) :: dtaua_tau
   integer, intent(in) :: nrt, ix
!*** Output
   real(double), intent(out) :: rint_ix
!*** Local variables
   integer :: j0, j1, j2, k, l
   real(double) :: scale
   real(double) :: delta_taua(nrt)
   integer, parameter :: npow = 3
   integer :: ipow  
   real(double), dimension(3,npow) :: kmat_poly
   real(double), dimension(3) :: y_poly, x_poly
!-----------------------------------------------------------------------------------
   
      k = 0
      if(ntau_2nd(ix).eq.1) then
         j0 = 1
         scale = (tau_grid(ix)+tau_grid_2nd(ix,j0))/taua_tot_all
         do l = 1,nrt
            delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j0,l))
         enddo
         rint_ix = rint_tau(j0) + &
            sum(dtaua_tau(ix,j0,:) * delta_taua)

         goto 999
      endif

!   j1=MINVAL(MINLOC(DABS(taua_tot_2nd-tau_grid_2nd(ix,:)),taua_tot_2nd.GE.tau_grid_2nd(ix,:)))
      j1 = minval(minloc(DABS(taua_tot_2nd-tau_grid_2nd(ix,:))))

      if(j1.eq.ntau_2nd(ix).and.taua_tot_2nd.ge.tau_grid_2nd(ix,j1)) then
         scale = (tau_grid(ix)+tau_grid_2nd(ix,j1))/taua_tot_all
!      scale = (tau_grid(ix))/taua_tot_all
         do l = 1,nrt
            delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j1,l))
         enddo
         rint_ix = rint_tau(j1) + &
            sum(dtaua_tau(ix,j1,:) * delta_taua)
         y_poly(2) = rint_tau(j1) + &
            sum(dtaua_tau(ix,j1,:) * delta_taua)
         j0 = j1 - 1
         scale = (tau_grid(ix)+tau_grid_2nd(ix,j0))/taua_tot_all
         do l = 1,nrt
           delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j0,l))
         enddo
         y_poly(1) = rint_tau(j0) + &
            sum(dtaua_tau(ix,j0,:) * delta_taua)
         rint_ix = y_poly(1) + &
            (y_poly(2)-y_poly(1))/(tau_grid_2nd(ix,j1)-tau_grid_2nd(ix,j0))*&
            (taua_tot_2nd-tau_grid_2nd(ix,j0)) 
        return
      endif

      if(j1.eq.1.and.taua_tot_2nd.le.tau_grid_2nd(ix,j1)) then
         scale = (tau_grid(ix)+tau_grid_2nd(ix,j1))/taua_tot_all
!      scale = (tau_grid(ix))/taua_tot_all
         do l = 1,nrt
            delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j1,l))
         enddo
         y_poly(1) = rint_tau(j1) + &
            sum(dtaua_tau(ix,j1,:) * delta_taua)
         j2 = j1 + 1
         scale = (tau_grid(ix)+tau_grid_2nd(ix,j2))/taua_tot_all
         do l = 1,nrt
            delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j2,l))
         enddo
         y_poly(2) = rint_tau(j2) + &
            sum(dtaua_tau(ix,j2,:) * delta_taua)
         rint_ix = y_poly(1) + &
            (y_poly(2)-y_poly(1))/(tau_grid_2nd(ix,j2)-tau_grid_2nd(ix,j1))*&
            (taua_tot_2nd-tau_grid_2nd(ix,j1))
         return
      endif
      j1 = max(j1,2)
      j1 = min(j1,ntau_2nd(ix)-1)
      j0 = j1 - 1
      j2 = j1 + 1  
      scale = (tau_grid(ix)+tau_grid_2nd(ix,j0))/taua_tot_all
 !  scale = (tau_grid(ix))/taua_tot_all
      do l = 1,nrt
         delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j0,l))
      enddo
      y_poly(1) = rint_tau(j0) + &
          sum(dtaua_tau(ix,j0,:) * delta_taua)  
      scale = (tau_grid(ix)+tau_grid_2nd(ix,j1))/taua_tot_all
 !  scale = (tau_grid(ix))/taua_tot_all
      do l = 1,nrt
         delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j1,l))
      enddo  
      y_poly(2) = rint_tau(j1) + &
          sum(dtaua_tau(ix,j1,:) * delta_taua)

      scale = (tau_grid(ix)+tau_grid_2nd(ix,j2))/taua_tot_all
 !  scale = (tau_grid(ix))/taua_tot_all
      do l = 1,nrt
         delta_taua(l) = scale*taua_arr(l) - (tau_grid_arr(ix,j2,l))
      enddo 
      y_poly(3) = rint_tau(j2) + &
          sum(dtaua_tau(ix,j2,:) * delta_taua)        
      forall(ipow=1:npow) kmat_poly(1,ipow) = tau_grid_2nd(ix,j0)**dble(ipow-1)
      forall(ipow=1:npow) kmat_poly(2,ipow) = tau_grid_2nd(ix,j1)**dble(ipow-1)
      forall(ipow=1:npow) kmat_poly(3,ipow) = tau_grid_2nd(ix,j2)**dble(ipow-1)     
      x_poly = &
          matmul(&
          INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)),&
          matmul(transpose(kmat_poly),y_poly))  
      rint_ix = &
          x_poly(1) + &
          x_poly(2)*taua_tot_2nd + &
          x_poly(3)*taua_tot_2nd**2 

    999 continue
   
   end subroutine interpolate_tau_2nd_pt

!***********************************************************************************

   subroutine interpolate_tau_1st( &
              tau_grid, &
              i0, &
              rint_i0, &
              i1, &
              rint_i1, &
              i2, &
              rint_i2, &
              taua_tot_1st, &
              rint_int)  
   use auxiliary_routines_module
   real(double), dimension(:), intent(in) :: tau_grid
   integer, intent(in) :: i0, i1, i2
   real(double), intent(in) :: rint_i0, rint_i1, rint_i2
   real(double), intent(in) :: taua_tot_1st 
   real(double), intent(out) :: rint_int 
!*** local variables
   integer, parameter :: npow = 3
   integer :: ipow
   real(double), dimension(3, npow) :: kmat_poly
   real(double), dimension(3) :: y_poly, x_poly
!-------------------------------------------------------    
      forall(ipow=1:npow) kmat_poly(1,ipow) = tau_grid(i0)**dble(ipow-1)
      forall(ipow=1:npow) kmat_poly(2,ipow) = tau_grid(i1)**dble(ipow-1)
      forall(ipow=1:npow) kmat_poly(3,ipow) = tau_grid(i2)**dble(ipow-1)
                                 
      y_poly(1) = rint_i0 
      y_poly(2) = rint_i1 
      y_poly(3) = rint_i2 
        
      x_poly = matmul( &
               INVERSE_3BY3(matmul(transpose(kmat_poly),kmat_poly)), &
               matmul(transpose(kmat_poly),y_poly))
      
      rint_int = &
         x_poly(1) + &
         x_poly(2)*taua_tot_1st + &
         x_poly(3)*taua_tot_1st**2 
    
   end subroutine interpolate_tau_1st

!***********************************************************************

end module tau_grid_module
