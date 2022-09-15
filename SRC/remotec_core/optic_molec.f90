module optic_molec_module
  use header_module
  use read_xsdb_module, only: cross_section_db, get_xs
  implicit none
  private

  !*** procedures
  public :: calculate_optic_mol_prop
  private :: rayleigh, rayleigh_phase_function

  !*** types
  public :: cross_section_db

contains

  !****************************************************************************
  !*** Calculate optical properties of molecules (Rayleigh scattering)
  !*** Previously the cross-sections were only calculated once and stored in the win structure, 
  !*** unless the cross-sections, oxygen total column or temperature offset were fitted.
  !*** For simplicity, we now calculate the cross-sections at every call to this routine.
  subroutine calculate_optic_mol_prop( &
       iband, XSFlag, O2Flag, TFlag,  &
       wavelength, &
       nrt, &
       atm_xs, &
       xsdb, &
       type_xsdb, &
       x_molec, &
       dvair, &
       dvair_old, &
       vmr_h2o, &
       play_old, & 
       taua_mol, &
       taus_mol, &
       taua_mol_species, &
       cross_mol, &
       csray, &
       dtaua_T, &
       dtaua_P, &
       plmom_ray, &
       ExitXSFlag, &
       ierr)
    !*** Input
    integer, intent(in) :: iband, XSFlag, O2Flag, TFlag, nrt
    real(double), dimension(:), intent(in) :: wavelength     ! wavelength grid (Dim: nwave
    type(atmosphere), intent(in) :: atm_xs                   ! Atmosphere grid for cross-sections
    real(double), dimension(:,:), intent(in) :: x_molec      ! Partial columns of absorbers, (Dim: ntype,natm)
    type(cross_section_db), dimension(:), intent(in) :: xsdb ! cross-section database
    integer, dimension(:)  :: type_xsdb
    real(double), dimension(:), intent(in) :: dvair          ! Partial air column, subject to change in O2 retrieval (Dim: natm)
    real(double), dimension(:), intent(in) :: dvair_old      ! Partial air column (Dim: natm)
    real(double), dimension(:), intent(in) :: play_old       ! Pressure, layer center (Dim: natm)   
    !*** Output 
    real(double), dimension(:,:), allocatable, intent(out) :: taua_mol, taus_mol   
    real(double), dimension(:,:,:), allocatable, intent(out) ::  taua_mol_species, cross_mol        
    real(double), dimension(:), allocatable, intent(out) :: csray                  
    real(double), dimension(:,:), allocatable, intent(out) :: dtaua_T, dtaua_P                                 
    real(double), dimension(:,:,:,:), intent(out) :: plmom_ray  
    integer, intent(out) :: ExitXSFlag, ierr
    !*** local variables
    integer :: i, j, k, l, m, imid, imol, natm, nwave, ntype    
    real(double) ::  dT, dP, wave_mid
    real(double), dimension(atm_xs%n) :: vmr_h2o    
    real(double), dimension(:,:), allocatable :: norm
    real(double), dimension(atm_xs%n) :: tlay_Tper
    real(double), dimension(atm_xs%n) :: play_Pper
    real(double), dimension(:,:,:), allocatable :: cross_section      ! Molecular absorption cross sections (Dim: nwave_hi,natm,ntype)
    real(double), dimension(:,:,:), allocatable :: cross_section_Tper ! Molecular absorption cross sections, perturbed by infinitesimal DT
    real(double), dimension(:,:,:), allocatable :: cross_section_Pper ! Molecular absorption cross sections, perturbed by infinitesimal DP
    real(double), dimension(atm_xs%n) :: pcor
    integer     :: io
    character*2 :: ch
    character(199) :: message
    !*** absco
    real(double) ::  dvmr
    real(double), dimension(atm_xs%n) :: vmr_h2o_per    
    real(double), dimension(:,:,:), allocatable :: cross_section_h2o_per ! Molecular absorption cross sections, perturbed by h2o vmr
    !------------------------------------------------------ 
    !*** Initialize error identifiers
    ierr = 0
    ExitXSFlag = 0

    !*** Set dimensions
    nwave = size(wavelength)
    natm = atm_xs%n
    ntype = size(xsdb)
    
    if(TFlag==1) then
       dT = 1.D0
       do i = 1, natm
          tlay_Tper(i) = atm_xs%t(i) + dT
       enddo
    endif
    if(O2Flag==1) then
       dP = 1.001
       do i = 1, natm
          play_Pper(i) = atm_xs%p(i)*dP
       enddo
    endif

!*** ABSCO
    if(type_xsdb(1)==4)then
      dvmr = 0.001! which is reasonable small d?
      do i = 1, natm
         vmr_h2o_per(i) = vmr_h2o(i) + dvmr
      enddo
    endif

    allocate(cross_section(nwave, natm, ntype), cross_section_h2o_per(nwave, natm, ntype),stat=ierr)
    cross_section = 0.D0
    cross_section_h2o_per = 0.D0       
    if (ierr .ne. 0) then
       write(message, *) 'CALCULATE_OPTIC_MOL_PROP: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    if(XSFlag>0) THEN
       write(ch,'(i2.2)')iband
       open(newunit(io),FILE='./CONTRL_OUT/cross_section_'//ch//'.dat')
       do i=1,ntype
          do k=1,nwave
             read(io,*)cross_section(k,:,i)
          enddo
       enddo
       close(io)

    else
       do i = 1, ntype 
          pcor=1.
          !*** This is ugly: Self broadening of H2O is non-negligible in some spectral bands. 
          !*** The Lorentz self broadening half-width gamma_self is roughly 5*gamma_air with scatter.
          !*** Since the XS lookup-table does not include a self-broadening dimension, we scale the
          !*** atmospheric pressure when reading the H2O cross sections by a factor [1+4*vmr_H2O]. 
          !*** This assumes that gamma_self/gamma_air=5 and that Doppler-broadening effects
          !*** is negligeable in the height ranges where the vmr_H2O is high [cf. C.Frankenberg, R. Scheepmaker]
          if(abs(xsdb(i)%species)==1 .or. (abs(xsdb(i)%species)>=100 .and. abs(xsdb(i)%species)<=199) .and. type_xsdb(1)/=4)then
             pcor = 1. + 4. * x_molec(:,i)/dvair
             do k = 1, natm
                if(pcor(k)*atm_xs%p(k)>1100.)pcor(k)=1100./atm_xs%p(k)
             enddo
          else
             pcor=1.
          endif

          do k = 1, natm
             call get_xs( &
                  xsdb(i),&                 
                  vmr_h2o(k), &
                  atm_xs%p(k)*pcor(k), &
                  atm_xs%t(k), &
                  cross_section(1:nwave,k,i), &
                  ExitXSFlag, ierr)
             if (ExitXSFlag==1 .or. ierr .ne. 0) return
          enddo

!*** influence of h2o self bordening ABSCO
          if((abs(xsdb(i)%species)==1 .or. (abs(xsdb(i)%species)>=100 .and. abs(xsdb(i)%species)<=199)) .and. type_xsdb(1)==4)then
            do k = 1, natm
               call get_xs( &
                    xsdb(i),&                 
                    vmr_h2o_per(k), &
                    atm_xs%p(k)*pcor(k), &
                    atm_xs%t(k), &
                    cross_section_h2o_per(1:nwave,k,i), &
                    ExitXSFlag, ierr)
               if (ExitXSFlag==1 .or. ierr .ne. 0) return
            enddo

            do k = 1, natm!*** tau + N*dtau/dN
                    cross_section(1:nwave,k,i) = cross_section(1:nwave,k,i) + &
	          			      vmr_h2o(k)*(cross_section_h2o_per(1:nwave,k,i)-cross_section(1:nwave,k,i))/dvmr 
	    enddo
          endif



       enddo

       if(TFlag==1) then
          allocate(cross_section_Tper(nwave, natm, ntype), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'CALCULATE_OPTIC_MOL_PROP: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
          do i = 1, ntype 
             if(abs(xsdb(i)%species)==1 .or. (abs(xsdb(i)%species)>=100 .and. abs(xsdb(i)%species)<=199).and. type_xsdb(1)/=4)then
                pcor = 1. + 4. * x_molec(:,i)/dvair
                do k = 1, natm
                   if(pcor(k)*atm_xs%p(k)>1100.)pcor(k)=1100./atm_xs%p(k)
                enddo
             else
                pcor=1.
             endif
             do k = 1, natm
                call get_xs( &
                     xsdb(i),&                      
                     vmr_h2o(k), &                    
                     atm_xs%p(k)*pcor(k), &
                     tlay_Tper(k), &
                     cross_section_Tper(1:nwave,k,i), &
                     ExitXSFlag, ierr)
                if (ExitXSFlag==1 .or. ierr .ne. 0) return
             enddo
          enddo
       endif
       if(O2Flag==1) then
          allocate(cross_section_Pper(nwave, natm, ntype), stat=ierr)
          if (ierr .ne. 0) then
             write(message, *) 'CALCULATE_OPTIC_MOL_PROP: memory allocation error'
             ierr = ierr_all
             goto 999
          endif
          do i = 1, ntype 
             if(abs(xsdb(i)%species)==1 .or. (abs(xsdb(i)%species)>=100 .and. abs(xsdb(i)%species)<=199).and. type_xsdb(1)/=4)then
                pcor = 1. + 4. * x_molec(:,i)/dvair
                do k = 1, natm
                   if(play_Pper(k)*pcor(k)>1100.)pcor(k)=1100./play_Pper(k)
                enddo
             else
                pcor=1.
             endif
             do k = 1, natm
                call get_xs( &
                     xsdb(i),& 
                     vmr_h2o(k), &
                     play_Pper(k)*pcor(k), &
                     atm_xs%t(k), &
                     cross_section_Pper(1:nwave,k,i), &
                     ExitXSFlag, ierr)
                if (ExitXSFlag==1 .or. ierr.ne.0 ) return
             enddo
          enddo
       endif

       if(XSFlag<0)then
          write(ch,'(i2.2)')iband
          open(newunit(io),file='CONTRL_OUT/cross_section_'//ch//'.dat')
          do i=1,ntype
             do k=1,nwave
                write(io,'(300(E17.10,X))')cross_section(k,:,i)
             enddo
          enddo
          close(io)
       endif
    endif
    !      write(*,*) 'Opened output file: ', 'CONTRL_OUT/cross_section_'//ch//'.dat', ' in remotec_core/optic_molec/calculate_optic_mol_prop (line 169)'

    if(allocated(taua_mol)) deallocate(taua_mol)
    if(allocated(taua_mol_species)) deallocate(taua_mol_species)
    if(allocated(taus_mol)) deallocate(taus_mol)
    if(allocated(cross_mol)) deallocate(cross_mol)
    if(allocated(dtaua_P)) deallocate(dtaua_P)
    if(allocated(dtaua_T)) deallocate(dtaua_T)
    if(allocated(norm)) deallocate(norm)
    if(allocated(csray)) deallocate(csray)     
    allocate(&
         taua_mol(nwave,nrt), &
         taus_mol(nwave,nrt), &
         taua_mol_species(nwave,nrt,ntype), &
         cross_mol(nwave,nrt, ntype), & 
         dtaua_P(nwave,nrt), &
         dtaua_T(nwave,nrt), &
         norm(nrt, ntype),&
         csray(nwave), &
         stat=ierr)
    if (ierr .ne. 0) then
       write(message, *) 'CALCULATE_OPTIC_MOL_PROP: memory allocation error'
       ierr = ierr_all
       goto 999
    endif

    call rayleigh(wavelength, csray, nwave)
    imid = (nwave+1)/2
    wave_mid = wavelength(imid)
    call rayleigh_phase_function(wave_mid, plmom_ray, nrt)      
    m = natm/nrt
    do i = 1,nwave
       do j = 1,nrt
          do k = 1,ntype
             taua_mol_species(i,j,k) = 0.D0
             cross_mol(i,j,k) = 0.D0
          enddo
          dtaua_T(i,j) = 0.D0
          dtaua_P(i,j) = 0.D0
       enddo
    enddo
    do j = 1,nrt
       do k = 1, ntype
          norm(j,k) = 0.D0
       enddo
    enddo

!!!       do n = 1, nwin
!!!          do i = 1, win_ini(n)%ntype
!!!             if(abs(win_ini(n)%xsdb(i)%species)==1 .or. (abs(win_ini(n)%xsdb(i)%species)>=100 .and. abs(win_ini(n)%xsdb(i)%species)<=199))then       
!!!                vmr_h2o(:) =  win(n)%x_molec(:,i)/dvair(:) ! caculate vmr of water per layer
!!!!!$                vmr_h2o(:) =  win(n)%dv_x(:,i)/dvair(:) ! Use initial VMR
!!!                goto 101
!!!             endif
!!!          enddo               
!!!       enddo

    do imol = 1, ntype         
       do k = 0, natm-1
          j = int(k/m) + 1
          taua_mol_species(:,j,imol) = taua_mol_species(:,j,imol) + &
               cross_section(:,k+1,imol)*x_molec(k+1,imol)
          cross_mol(:,j,imol) = cross_mol(:,j,imol) + &
               cross_section(:,k+1,imol)*x_molec(k+1,imol)    
          norm(j,imol) = norm(j,imol) + &
               x_molec(k+1,imol)
          if(TFlag==1) dtaua_T(:,j) = dtaua_T(:,j) + &
               (cross_section_Tper(:,k+1,imol) - cross_section(:,k+1,imol))/&
               (tlay_tper(k+1) - atm_xs%t(k+1))*x_molec(k+1,imol)   
          if(O2Flag==1) dtaua_P(:,j) = dtaua_P(:,j) + &
               (cross_section_Pper(:,k+1,imol) - cross_section(:,k+1,imol))/&
               (play_Pper(k+1) - atm_xs%p(k+1))*x_molec(k+1,imol)*play_old(k+1)/dvair_old(k+1)/relo2/m    
       enddo
       do k = 1, nrt
          if(norm(k,imol)>0)then
             do i = 1,nwave
                cross_mol(i,k,imol) =  cross_mol(i,k,imol)/norm(k,imol)  
             enddo
          else
             do i = 1, nwave
                cross_mol(i,k,imol) = 0.D0
             enddo
          endif
       enddo
    enddo

    do i = 1,nwave
       do j = 1,nrt
          taus_mol(i,j)=0.D0
          taua_mol(i,j)=0.D0
       enddo
    enddo
    do k = 0, natm-1
       j = int(k/m)+1
       do i = 1,nwave 
          taus_mol(i,j) = taus_mol(i,j) + csray(i)*dvair(k+1)
       enddo
    enddo

    ! Change by Otto (11-2-2013): Do not cut off optical depth profiles
    ! at max value tatot but only do this for k-binning grid 
    ! (see rad_trans_intf.f90)
    !      do l = 1, nwave
    !         do k = 1, nrt
    !            if(sum(taua_mol_species(l,1:k,:))>tatot) then
    !               taua_mol_species(l,k+1:nrt,:) = 0.D0
    !               taua_mol_species(l,1:k,:) = taua_mol_species(l,1:k,:)*tatot / sum(taua_mol_species(l,1:k,:))
    !               exit
    !           endif 
    !         enddo
    !      enddo

    do l = 1, nwave
       do k = 1, nrt
          taua_mol(l,k) = sum(taua_mol_species(l,k,:))
       enddo
    enddo

    return

999 continue
    call stopretrieval(message)

  end subroutine calculate_optic_mol_prop

!**********************************************************************
   subroutine RAYLEIGH(WAVELENGTH, CROSS_RAY, NWAVE)
!*** Input
   integer, intent(in) :: nwave
   real(double), intent(in) :: WAVELENGTH(NWAVE)
!*** Output
   real(double), intent(out) :: CROSS_RAY(NWAVE)
!*** local variables
   integer :: l
   real(double) :: wmu, xp
!--------------------------------------------------------------------
      do L = 1, NWAVE
         Wmu           = WAVELENGTH(L) *1.d-3         !change units mn -> um
	 XP            = 0.389*Wmu+0.09426/Wmu-0.3228
         CROSS_RAY(L)  = 4.02d-28/Wmu**(4.+XP)
      enddo
      
      return
      end subroutine RAYLEIGH

!************************************************************************
!*** scattering and absorption optical depths per layer for the       *
!*** given GOME wavelength grid using effective cross sections        *
!************************************************************************
      subroutine RAYLEIGH_PHASE_FUNCTION(WAVE, PLMOM_RAY, nrt)   
        integer, intent(in) :: nrt   
        real(double), intent(in) :: wave
        !*** Output  
        real(double),  intent(out) :: PLMOM_RAY(NSTOKES, NSTOKES, 0:2, NRT)
        !*** Local variables
        real(double) :: plmom(3,3,0:2,nrt)   
        integer :: k, l, i_st, j_st
        real(double) :: c0, depol
        !--------------------------------------------------------------
        ! phase matrix
        C0 = (2. - 2.*DEPOL(WAVE))/(2. + DEPOL(WAVE))
        do k = 1, nrt
           do l = 0, 2
              do I_ST = 1,NSTOKES
                 do J_ST = 1,NSTOKES
                    PLMOM(J_ST, I_ST, L, K) = 0.
                 enddo
              enddo
           enddo
        enddo

        do k = 1, nrt
           PLMOM(1, 1, 0, K) = 1.
           PLMOM(1, 1, 2, K) = 0.5*C0
           !PLMOM(2, 1, 0, K) = -0.5	!rrae !JadB: Trijf.
           !PLMOM(2, 1, 2, K) = 0.5	!rrae !JadB: Trijf.
           PLMOM(2, 1, 2, K) = sqrt(1.5) * C0 !  0.5	!rrae
           !PLMOM(1, 2, 0, K) = -0.5	!rrae !JadB: Trijf.
           ! PLMOM(1, 2, 2, K) = 0.5	!rrae !JadB: Trijf.
           PLMOM(1, 2, 2, K) = sqrt(1.5) * C0 !  0.5	!rrae
           !PLMOM(2, 2, 0, K) = 1.	!rrae !JadB: Trijf.
           !PLMOM(2, 2, 2, K) = 0.5*C0	!rrae !JadB: Trijf.
           PLMOM(2, 2, 2, K) = 3.0*C0	!rrae
           ! PLMOM(3, 3, 1, K) = 3./2	!rrae !JadB: Trijf.
        enddo

        PLMOM_RAY(1:nstokes, 1:nstokes, 0:2, 1:nrt) = plmom(1:nstokes, 1:nstokes, 0:2, 1:nrt)     
        return
      end subroutine RAYLEIGH_PHASE_FUNCTION

!**************************************************************************************************
end module optic_molec_module
