!------------------------------------------------------------------------------ 
!------------------------------------------------------------------------------ 
module rad_intf_module
  use header_gsd
  use gautool_acs_module
  use gaus_pol_module
  use viewing_module
  use single_scat_module
  use perturbation_module
  implicit none 

contains
  !------------------------------------------------------------------------------ 
  !> initialization of the radiation transsport model  
  !! 1. gaussian quadratures 
  !! 2. spherical functions
  !------------------------------------------------------------------------------ 
  subroutine PERTURBATION_INIT(&                                     
       U0,&
       THETA,&
       gauss_quad,&
       gsf,&
       ierr)
    !*** input
    real(double), intent(inout) :: u0
    real(double), intent(in):: theta
    type(gauss_quad_type), intent(out) :: gauss_quad
    !*** output
    type(gsf_type), intent(out) :: gsf
    integer, intent(out) :: ierr

    !*** local variables
    real(double), parameter :: EPS = 1.D-4
    integer :: i
    real(double) :: uv
    !------------------------------------------------------------------------------ 
    !*** Initialize error identifier
    ierr = 0
    !*** Fix Gaussian quadrature, viewing and solar zenith and azimuthal
    !*** cosines

    call DOUBLE_GAUSS(gauss_quad, ierr)            !gaussian quadrature
    if (ierr .ne. 0) return

    UV = DCOS(THETA/180.*PI) 

    !*** check that U0 and UV are not identical with the gaussian points
    do I = 1,MAXSTR
       if(abs(U0-gauss_quad%DG_MU(I)).lt.EPS) U0 = U0 - EPS
       if(abs(UV-gauss_quad%DG_MU(I)).lt.EPS) UV = UV - EPS
    enddo

    !*** the generalized spherical functions
    call SET_PLM(&
         UV,&
         U0,&
         gauss_quad,&
         gsf)

    return

  end subroutine PERTURBATION_INIT

  !------------------------------------------------------------------------------

  subroutine FWD_ADJ_PERTURBATION(&
       nrt, &
       TAUA_IN,&    
       TAUS_IN, & 
       PHASE_IN,& 
       pm_in,&
       gauss_quad,&
       gsf,&
       NDER,& 
       MAXD,  & 
       MAXCOEFS,& 
       BDRF_MS,&    
       U0,&   
       THETA, & 
       PHI,     & 
       RINT_MSC, & 
       DTAUA,&  
       DTAUS,&    
       DPHASE,&    
       DALB,&
       ierr)
    !*** input     
    type(phase_mat_type),intent(in), dimension(0:maxstr) :: pm_in
    type(gauss_quad_type), intent(in) :: gauss_quad
    type(gsf_type), intent(in) :: gsf
    integer, intent(IN) :: nrt, MAXD              !coresponding layer index
    integer, dimension(maxd), intent(IN):: nder !coresponding layer index
    real(double), dimension(nrt), intent(in) :: taua_in, taus_in
    real(double), dimension(NSTOKES, NSTOKES, 0:MAXLEG, nrt), intent(in) :: PHASE_IN
    integer, intent(in), dimension(nrt) :: MAXCOEFS
    real(double), dimension(NSTOKES, NSTOKES, MXHALF+2, MXHALF+2, 0:MAXSTR-1), &
         intent(in) :: BDRF_MS
    real(double), intent(IN) :: u0, theta, phi
    !*** output
    real(double), intent(OUT), dimension(nstokes) :: RINT_MSC   ! intensity vector at TOA 
    real(double), intent(OUT), dimension(NSTOKES,nrt) :: dtaua, dtaus         !dI/dTAU, dI/dW0
    real(double), intent(OUT), dimension(NSTOKES,0:MAXLEG,NPER,MAXD) :: dphase   !dI/dPHASE
    real(double), intent(OUT), dimension(nstokes) :: dalb                        !dI/dALBEDO
    integer, intent(out) :: ierr
    !*** local variables
    real(double), dimension(nrt) :: taua, taus, taua_nk, taus_nk
    integer :: i_st, i, k, ifwd_adj, ik, imom, isrc_fwd, nmat, isrc_adj, iper, l, n, n1, n2
    real(double) :: tau, tau_in, uv, w0, eps1
    real(double), dimension(nstokes) :: rint_pre !intensity vector at TOA of previous iteration 
    integer :: NF  !index of Fourier component
    logical :: SW_SURF              !switch for surface reflectance     
    !*** filter functions for auto. cell spliting
    integer, dimension(MXK) :: nfilt1
    integer, dimension(nrt) :: nfilt2
    integer :: MNK                  !number of internal layers
    integer, dimension(nrt) :: NCOEFS       !number of relevant phase fct coeff.
    real(double), dimension(mxk) :: tautot !total optical depth
    real(double), dimension(0:mxk) ::eintuv !exp(-TAUINT(K)/UV) 
    real(double), dimension(0:mxk) ::eintu0 !attenuation factors EXP(-TAUINT(K)/U0)
    real(double), dimension(nrt) :: omega !single scattering albedo
    real(double), dimension(mxk) :: expu0  !attenuation factors EXP(-tau(k)/U0)
    real(double), dimension(mxk) :: expuv !attenuation factors EXP(-tau(k)/UV)
    real(double), dimension(MXHALF, MXK) :: expmu !attenuation factors EXP(-tau(k)/mu(i))
    !*** forward and adjoint intensity field
    real(double), dimension(NSTOKES, MXHALF, 0:MXK) :: RFWD_UP, RFWD_DN, RADJ_UP, RADJ_DN
    real(double), dimension(NSTOKES, 0:MXK) :: RADJ_VW, RFWD_VW
    !*** For Nakajima(1988) MS correction
    real(double), dimension(nstokes) ::  RSS_DM, RSS_NK 
    real(double), dimension(NSTOKES, nrt) ::&
         DTAUT_DM_FOU,&
         DTAUS_DM_FOU,&
         DTAUT_NK_FOU,&
         DTAUS_NK_FOU
    real(double), dimension(NSTOKES, 0:MAXLEG, NPER, MAXD) :: &
         DPHASE_DM_FOU,&
         DPHASE_NK_FOU
    real(double), dimension(nstokes) :: DALB_DM,DALB_NK
    real(double), dimension(NSTOKES, nrt) ::&
         DTAUT_DM,&
         DTAUS_DM,&
         DTAUT_NK,&
         DTAUS_NK
    real(double), dimension(NSTOKES, 0:MAXLEG, NPER, MAXD) :: DPHASE_DM, DPHASE_NK
    real(double), dimension(NSTOKES, nrt) :: dtaut
    real(double), dimension(0:maxleg) :: COSMPHI, SINMPHI, delfac
    real(double), dimension(4) :: DELPL, DELMI
    real(double), dimension(4, 0:MAXLEG) :: bplus, bmin
    real(double), dimension(nrt) :: f, g ! for Nakajima correction
    !*** determines trunction of Gauss-Seidel iteration loop and Fourier loop
    real(double), parameter ::  EPST = 1.D-4   
    !*** determines truncation of phase matrix expansion
    real(double), parameter :: EPSP = 1.D-6      
    !------------------------------------------------------------------------------
    !*** Initialize error identifier
    ierr = 0
    !*** setup for source index indication  which Stokes component of the
    !*** source is non zero. Should be always 1 for the forward mode.
    ISRC_FWD = 1

    !*** calculate the diagonal matrices 1-DELTA and 1+DELTA  
    DELMI(1) = 0.0D0
    DELMI(2) = 0.0D0
    DELMI(3) = 2.0D0
    DELMI(4) = 2.0D0
    DELPL(1) = 2.0D0
    DELPL(2) = 2.0D0
    DELPL(3) = 0.0D0
    DELPL(4) = 0.0D0

    UV = DCOS(THETA/180.*PI) 

    !*** calculate the basis vectors B+ and B- of the Fourier expansion
    DTAUT = 0.D0
    do NF = 0, MAXLEG
       COSMPHI(NF) = DCOS(dble(NF)*(PHI/180.*PI))
       SINMPHI(NF) = DSIN(dble(NF)*(PHI/180.*PI))
       BPLUS(1,NF) =  COSMPHI(NF)
       BPLUS(2,NF) =  COSMPHI(NF)
       BPLUS(3,NF) =  SINMPHI(NF)
       BPLUS(4,NF) =  SINMPHI(NF)
       BMIN(1,NF)  = -SINMPHI(NF)
       BMIN(2,NF)  = -SINMPHI(NF)
       BMIN(3,NF)  =  COSMPHI(NF)
       BMIN(4,NF)  =  COSMPHI(NF)
       DELFAC(NF) = 2.d0               !factor 0.5*(2-delta(m,0)) 
       if(NF.eq.0)DELFAC(NF) = 1.D0
    enddo

    do I_ST = 1, NSTOKES
       RINT_PRE(I_ST) = 0.D0
       RINT_MSC(I_ST) = 0.D0
    enddo

    do K = 1, nrt
       do I_ST = 1,NSTOKES
          DTAUT(I_ST,K)    = 0.D0
          DTAUS(I_ST,K)    = 0.D0
          DTAUT_DM(I_ST,K) = 0.D0
          DTAUS_DM(I_ST,K) = 0.D0
          DTAUT_NK(I_ST,K) = 0.D0
          DTAUS_NK(I_ST,K) = 0.D0
       enddo
    enddo

    do K = 1,MAXD
       do IPER = 1,NPER
          do I = 0, MAXLEG  
             do I_ST = 1,NSTOKES
                DPHASE(I_ST,I,IPER,K)    = 0.D0
                DPHASE_DM(I_ST,I,IPER,K) = 0.D0
                DPHASE_NK(I_ST,I,IPER,K) = 0.D0
             enddo
          enddo
       enddo
    enddo

    !*** for the calculation of the diffuse field use only MAXSTR-1 
    !*** coefficients of the phase function. Higher moments are only
    !*** taken into account in the single scattering field. Due to that, 
    !*** the phase matrix does not to be renormalized. Seems to work fine for 
    !*** at least 16 streams. 
    do K = 1, nrt
       NCOEFS(K) = min(MAXCOEFS(K), MAXSTR-1)
    enddo

    do K = 1, nrt
       TAUS(K) = TAUS_IN(K)
       TAUA(K) = TAUA_IN(K)
       TAUS_NK(K) = TAUS_IN(K)
       TAUA_NK(K) = TAUA_IN(K)
    enddo

    do IK = 1,MAXD
       K=NDER(IK)
       TAU_IN = TAUS_IN(K) + TAUA_IN(K)
       TAU = TAU_IN
       W0 = TAUS_IN(K)/TAU_IN
       !*** Nakajima switch: off = (f(k)=0.), on = (f(k)!=0)
       !         f(k) = 0.
       f(k) = PHASE_IN(1,1,MAXSTR,K)/(2.*dble(MAXSTR) + 1.D0)
       g(k) = f(k)* W0
       TAU = TAU_IN - f(k)*TAUS_IN(K)
       TAUS(K) = (1.-f(k))*TAUS_IN(K)
       TAUA(K) = TAU-TAUS(K)
       TAUS_NK(K) = (1.-f(k))*TAUS_IN(K)
       TAUA_NK(K) = TAU_IN - TAUS_NK(K)
    enddo

    call ATTEN_FAC(&   
         nrt, &  
         TAUA, &  
         TAUS,  &
         U0,    &  
         UV, &
         gauss_quad,&
         OMEGA, &
         TAUTOT,   &
         EXPMU,  &
         EXPU0,&
         EXPUV, &
         EINTU0, &    
         EINTUV, &
         MNK,&
         NFILT1, & 
         NFILT2,  &  
         SW_SURF, &
         ierr)
    if (ierr .ne. 0) return
    

    !*** Fourier loop         
    do NF = 0, MAXSTR-1       
       NMAT = NSTOKES
       if(NF.eq.0) NMAT=min(2,NSTOKES)

       IFWD_ADJ  = 1      !flag for forward calculation
       ISRC_FWD  = 1      !Stokes index for forward source 

       call GAUSS_SEIDEL(&
            nrt, &
            NF,&
            MNK,&  
            NFILT1,&   
            NMAT,&
            IFWD_ADJ,& 
            ISRC_FWD,& 
            EPST, &
            U0, &
            SW_SURF,&
            pm_in(nf),&
            gauss_quad,&
            BDRF_MS, &
            OMEGA, &    
            EXPMU,&
            EXPU0, &
            EINTU0, &
            RFWD_UP,&   
            RFWD_DN,&
            ncoefs)

       call VIEWING( &
            gauss_quad,&
            pm_in(nf),&
            NF,&
            MNK,&     
            NFILT1,&
            NMAT,&  
            IFWD_ADJ,& 
            SW_SURF,& 
            BDRF_MS,& 
            OMEGA,& 
            EXPUV,& 
            RFWD_UP, &
            RFWD_DN,&
            RFWD_VW)

       call SINGLE_SCAT_FOU_DER(&
            gsf,&
            pm_in(nf),&
            NF,&      
            UV,&        
            U0,&
            nrt, &    
            TAUS, &
            TAUA,&
            BDRF_MS,&
            RSS_DM, & 
            DTAUT_DM_FOU,&
            DTAUS_DM_FOU,& 
            DPHASE_DM_FOU,&
            DALB_DM, & 
            NCOEFS,&     
            NDER, & 
            MAXD)

       call SINGLE_SCAT_FOU_DER(&
            gsf,&
            pm_in(nf),&
            NF,&
            UV,& 
            U0,& 
            nrt, &         
            TAUS_NK, &
            TAUA_NK,&
            BDRF_MS,&    
            RSS_NK, &   
            DTAUT_NK_FOU, &
            DTAUS_NK_FOU,&
            DPHASE_NK_FOU,&   
            DALB_NK,& 
            NCOEFS,&   
            NDER,&        
            MAXD)

       do I_ST = 1, NSTOKES
          RINT_PRE(I_ST) = RINT_MSC(I_ST)
          RINT_MSC(I_ST) =  RINT_MSC(I_ST) + &
               BPLUS(I_ST,NF)* DELFAC(NF)* &
               (RFWD_VW(I_ST,0)+ RSS_DM(I_ST) - RSS_NK(I_ST))
       enddo

       IFWD_ADJ = 2      !flag for adjoint calculation
       do ISRC_ADJ = 1, 1!Stokes index for adjoint source
          call GAUSS_SEIDEL(&
               nrt, &
               NF,&
               MNK,&
               NFILT1,&
               NMAT,&
               IFWD_ADJ,& 
               ISRC_ADJ,& 
               EPST,&   
               UV,& 
               SW_SURF,&
               pm_in(nf),&
               gauss_quad,&
               BDRF_MS,& 
               OMEGA, &  
               EXPMU,&
               EXPUV,&  
               EINTUV,&   
               RADJ_UP,&  
               RADJ_DN,&
               ncoefs)

          call VIEWING(&
               gauss_quad,&
               pm_in(nf),&
               NF,&
               MNK, &
               NFILT1,&   
               NMAT,&
               IFWD_ADJ,&
               SW_SURF,& 
               BDRF_MS,&
               OMEGA, &
               EXPU0,& 
               RADJ_UP,& 
               RADJ_DN,&
               RADJ_VW) 

          call PERT_INTEGRALS(&
               gauss_quad,&
               gsf,&
               pm_in(nf),&
               ISRC_ADJ,&
               ISRC_FWD,&
               NF,&
               NMAT,&
               NDER,&
               MAXD, &
               nrt, &
               NCOEFS,& 
               NFILT2,&
               MNK,&
               COSMPHI,& 
               SINMPHI,&   
               TAUTOT,&
               EXPMU,& 
               EXPU0, &
               EINTU0,&  
               EXPUV, &
               EINTUV, &
               OMEGA,&    
               RFWD_UP,&
               RFWD_DN,&
               RFWD_VW, &
               RADJ_UP,&  
               RADJ_DN,& 
               RADJ_VW,& 
               U0, &
               UV,  &    
               SW_SURF,&
               DTAUT, &
               DTAUS, & 
               DPHASE,& 
               DALB,&
               ierr)
       enddo

       if(NF.eq.0) then
          do I_ST = 1, NSTOKES
             DALB(I_ST) = DALB(I_ST)+DELFAC(NF)*COSMPHI(NF)*(&
                  DALB_DM(I_ST)-DALB_NK(I_ST))
          enddo
       endif

       do IK = 1, MAXD
          K = NDER(IK)
          do I_ST = 1, min(2, nstokes)
             DTAUT_DM(I_ST,K) = DTAUT_DM(I_ST,K) +&
                  DELFAC(NF)*COSMPHI(NF)*&
                  (DTAUT_DM_FOU(I_ST,K))
             DTAUS_DM(I_ST,K) = DTAUS_DM(I_ST,K) +&
                  DELFAC(NF)*COSMPHI(NF)* &
                  (DTAUS_DM_FOU(I_ST,K))
             DTAUT_NK(I_ST,K) = DTAUT_NK(I_ST,K) + &
                  DELFAC(NF)*COSMPHI(NF)* &
                  (DTAUT_NK_FOU(I_ST,K))
             DTAUS_NK(I_ST,K) = DTAUS_NK(I_ST,K) + &
                  DELFAC(NF)*COSMPHI(NF)* &
                  (DTAUS_NK_FOU(I_ST,K))
          enddo

          do I_ST = 3, nstokes
             DTAUT_DM(I_ST,K) = DTAUT_DM(I_ST,K) +&
                  DELFAC(NF)*SINMPHI(NF)*&
                  (DTAUT_DM_FOU(I_ST,K))
             DTAUS_DM(I_ST,K) = DTAUS_DM(I_ST,K) +&
                  DELFAC(NF)*SINMPHI(NF)* &
                  (DTAUS_DM_FOU(I_ST,K))
             DTAUT_NK(I_ST,K) = DTAUT_NK(I_ST,K) + &
                  DELFAC(NF)*SINMPHI(NF)* &
                  (DTAUT_NK_FOU(I_ST,K))
             DTAUS_NK(I_ST,K) = DTAUS_NK(I_ST,K) + &
                  DELFAC(NF)*SINMPHI(NF)* &
                  (DTAUS_NK_FOU(I_ST,K))
          enddo

          do L = NF, NCOEFS(K)
             DPHASE_DM(1,L,1,IK) = DPHASE_DM(1,L,1,IK) + &
                  DELFAC(NF)*COSMPHI(NF)* &
                  (DPHASE_DM_FOU(1,L,1,IK))
             do n=2,nstokes
                if (n==2) DPHASE_DM(n,L,nper,IK) = DPHASE_DM(n,L,nper,IK) + &
                     DELFAC(NF)*COSMPHI(NF)*&
                     (DPHASE_DM_FOU(n,L,nper,IK))

                if (n==3) DPHASE_DM(n,L,nper,IK) = DPHASE_DM(n,L,nper,IK) +&
                     DELFAC(NF)*SINMPHI(NF)*&
                     (DPHASE_DM_FOU(n,L,nper,IK))
             enddo

             DPHASE_NK(1,L,1,IK) = DPHASE_NK(1,L,1,IK) + &
                  DELFAC(NF)*COSMPHI(NF)* &
                  (DPHASE_NK_FOU(1,L,1,IK))

             do n=2,nstokes
                if (n==2) DPHASE_NK(n,L,nper,IK) = DPHASE_NK(n,L,nper,IK) +&
                     DELFAC(NF)*COSMPHI(NF)*&
                     (DPHASE_NK_FOU(n,L,nper,IK))

                if (n==3) DPHASE_NK(n,L,nper,IK) = DPHASE_NK(n,L,nper,IK) +&
                     DELFAC(NF)*SINMPHI(NF)*&
                     (DPHASE_NK_FOU(n,L,nper,IK))
             enddo
          enddo

       enddo

       if(NF.ge.2)then
          if(RINT_MSC(1).ne.0)then
             EPS1 = abs(RINT_PRE(1)-RINT_MSC(1))/RINT_MSC(1)
          else
             EPS1 = 0.
          endif
          if(EPS1.le.EPST) goto 111
       endif
    enddo

111 continue

    !*** Transform Nakajima-scaled values to unscaled values
    do I_ST = 1, NSTOKES
       do IK = 1, MAXD
          K = NDER(IK)
          TAU_IN = TAUS_IN(K) + TAUA_IN(K)
          W0 = TAUS_IN(K)/TAU_IN
          do L=0,MAXSTR-1

             !HH: Changes made to compile without out-of-bounds errors for nstokes=nper=1 
             DPHASE(I_ST,MAXSTR,1,IK) = DPHASE(I_ST,MAXSTR,1,IK) +&
                  DPHASE(I_ST,L,1,IK)/(2.*DBLE(MAXSTR)+1.) *&
                  (-(2.*DBLE(L)+1)/(1.-f(k))    + &
                  ((PHASE_IN(1,1,L,K) -(2.*DBLE(L)+1)*f(k))&
                  /(1-f(k))**2 ))

             DPHASE_NK(I_ST,MAXSTR,1,IK) = &
                  DPHASE_NK(I_ST,MAXSTR,1,IK) +&
                  DPHASE_NK(I_ST,L,1,IK)/(2.*DBLE(MAXSTR)+1.) *&
                  (-(2.*DBLE(L)+1)/(1.-f(k))    + &
                  ((PHASE_IN(1,1,L,K) -(2.*DBLE(L)+1)*f(k))&
                  /(1-f(k))**2 ))

             DPHASE_DM(I_ST,MAXSTR,1,IK) = &
                  DPHASE_DM(I_ST,MAXSTR,1,IK) +&
                  DPHASE_DM(I_ST,L,1,IK)/(2.*DBLE(MAXSTR)+1.) *&
                  (-(2.*DBLE(L)+1)/(1.-f(k))    + &
                  ((PHASE_IN(1,1,L,K) -(2.*DBLE(L)+1)*f(k))&
                  /(1-f(k))**2 ))

             do n=2,nper
                if (n<4) then
                   n1=n
                   n2=n
                else
                   n1=2
                   n2=1
                endif
                DPHASE_NK(I_ST,MAXSTR,1,IK) = &
                     DPHASE_NK(I_ST,MAXSTR,1,IK) +&
                     DPHASE_NK(I_ST,L,n,IK)/(2.*DBLE(MAXSTR)+1.) * (&
                     -PHASE_IN(n1,n2,L,K)*W0/(1.-f(k)) +&
                     PHASE_IN(n1,n2,L,K)*(1.-g(k))/(1.-f(k))**2)

                DPHASE_DM(I_ST,MAXSTR,1,IK) = &
                     DPHASE_DM(I_ST,MAXSTR,1,IK) +&
                     DPHASE_DM(I_ST,L,n,IK)/(2.*DBLE(MAXSTR)+1.) * (&
                     -PHASE_IN(n1,n2,L,K)*W0/(1.-f(k)) +&
                     PHASE_IN(n1,n2,L,K)*(1.-g(k))/(1.-f(k))**2)

                DPHASE(I_ST,MAXSTR,1,IK) = DPHASE(I_ST,MAXSTR,1,IK) +&
                     DPHASE(I_ST,L,n,IK)/(2.*DBLE(MAXSTR)+1.) * (&
                     -PHASE_IN(n1,n2,L,K)*W0/(1.-f(k)) +&
                     PHASE_IN(n1,n2,L,K)*(1.-g(k))/(1.-f(k))**2)
             enddo

          enddo

          DPHASE(I_ST,MAXSTR,1,IK) = DPHASE(I_ST,MAXSTR,1,IK) &
               +DTAUT(I_ST,K)*(-TAUS_IN(K))/(2.*dble(MAXSTR)+1.D0) &
               +DTAUS(I_ST,K)*(-TAUS_IN(K))/(2.*dble(MAXSTR)+1.D0)
          DPHASE_DM(I_ST,MAXSTR,1,IK) = DPHASE_DM(I_ST,MAXSTR,1,IK)&
               +DTAUT_DM(I_ST,K)*(-TAUS_IN(K))/(2.*dble(MAXSTR)+1.D0)&
               +DTAUS_DM(I_ST,K)*(-TAUS_IN(K))/(2.*dble(MAXSTR)+1.D0)
          DPHASE_NK(I_ST,MAXSTR,1,IK) = DPHASE_NK(I_ST,MAXSTR,1,IK)&
               +DTAUS_NK(I_ST,K)*(-TAUS_IN(K))/(2.*dble(MAXSTR)+1.D0)

          DTAUS(I_ST,K) = DTAUS(I_ST,K)*(1.-f(k)) + DTAUT(I_ST,K)*(-f(k))  
          DTAUS_DM(I_ST,K) = DTAUS_DM(I_ST,K)*(1.-f(k)) + DTAUT_DM(I_ST,K)*(-f(k))  
          DTAUS_NK(I_ST,K) = DTAUS_NK(I_ST,K)*(1.-f(k)) 

          DTAUT_DM(I_ST,K) =  DTAUT_DM(I_ST,K)
          DTAUT_NK(I_ST,K) =  DTAUT_NK(I_ST,K)

          do n = 2, nper
             if (n<4) then
                n1=n
                n2=n
             else
                n1=2
                n2=1
             endif
             DO L=0,MAXSTR-1               
                DTAUS(I_ST,K) = DTAUS(I_ST,K) + &                   
                     (-f(k))/(TAU_IN*(1.-f(k)))*(&
                     DPHASE(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K))

                DTAUS_DM(I_ST,K) = DTAUS_DM(I_ST,K) +&
                     (-f(k))/(TAU_IN*(1.-f(k)))*&
                     DPHASE_DM(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K) 

                DTAUS_NK(I_ST,K) = DTAUS_NK(I_ST,K) +&
                     (-f(k))/(TAU_IN*(1.-f(k)))*&
                     DPHASE_NK(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K) 

                DTAUT(I_ST,K) = DTAUT(I_ST,K) +&
                     f(k)*TAUS_IN(K)/(TAU_IN**2*(1.-f(k)))*&
                     DPHASE(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K) 

                DTAUT_DM(I_ST,K) = DTAUT_DM(I_ST,K) +&
                     f(k)*TAUS_IN(K)/(TAU_IN**2*(1.-f(k)))*&
                     DPHASE_DM(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K) 

                DTAUT_NK(I_ST,K) = DTAUT_NK(I_ST,K) +&
                     f(k)*TAUS_IN(K)/(TAU_IN**2*(1.-f(k)))*&
                     DPHASE_NK(I_ST,L,n,IK)*PHASE_IN(n1,n2,L,K)                   

             enddo
          ENDDO


          do IMOM = 0, NCOEFS(K)
             do IPER = 1, 1
                DPHASE(I_ST,IMOM,IPER,IK) = &
                     DPHASE(I_ST,IMOM,IPER,IK)& 
                     /(1.-f(k))
                DPHASE_DM(I_ST,IMOM,IPER,IK) = &
                     DPHASE_DM(I_ST,IMOM,IPER,IK)& 
                     /(1.-f(k))
                DPHASE_NK(I_ST,IMOM,IPER,IK) = &
                     DPHASE_NK(I_ST,IMOM,IPER,IK)& 
                     /(1.-f(k))
             enddo
             do IPER = 2, NPER
                DPHASE(I_ST,IMOM,IPER,IK) = &
                     DPHASE(I_ST,IMOM,IPER,IK) &
                     *(1.D0-g(k))/(1.D0-f(k))
                DPHASE_DM(I_ST,IMOM,IPER,IK) = &
                     DPHASE_DM(I_ST,IMOM,IPER,IK) &
                     *(1.D0-g(k))/(1.D0-f(k))
                DPHASE_NK(I_ST,IMOM,IPER,IK) = &
                     DPHASE_NK(I_ST,IMOM,IPER,IK) &
                     *(1.D0-g(k))/(1.D0-f(k))
             enddo
          enddo

          do IPER = 1, NPER
             do L = 0, NCOEFS(K)+1
                DPHASE(I_ST, L, IPER, IK) = DPHASE(I_ST, L, IPER, IK) + &
                     DPHASE_DM(I_ST, L, IPER, IK) - &
                     DPHASE_NK(I_ST, L, IPER, IK)
             enddo
          enddo
       enddo
    enddo

    do I_ST = 1, NSTOKES
       DTAUS(I_ST,:) = DTAUS(I_ST,:) + DTAUS_DM(I_ST,:)&
            -DTAUS_NK(I_ST,:) 
       DTAUT(I_ST,:) = DTAUT(I_ST,:) + DTAUT_DM(I_ST,:)&
            -DTAUT_NK(I_ST,:)    
    enddo

    !*** Do transformation from se
    !*** dI/dtautot with tausca = const, dI / dtausca with tautot = const
    !*** to the derivative set
    !*** dI/dtauabs with tausca = const, dI / dtausca with tauabs = const
    do K = 1, nrt    
       do I_ST = 1, NSTOKES
          DTAUA(I_ST,K) =  DTAUT(I_ST,K)
          DTAUS(I_ST,K) =  DTAUS(I_ST,K) + DTAUT(I_ST,K)
       enddo
    enddo



    return

  end subroutine FWD_ADJ_PERTURBATION

  !------------------------------------------------------------------------------   

end module rad_intf_module
