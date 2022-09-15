module single_scat_module
   use header_gsd
   implicit none 

   contains
!------------------------------------------------------------------------------   
     subroutine SINGLE_SCAT(&
          nrt, &
          TAUA,&
          TAUS,&
          BDRF_SS,&
          U0,&
          UV,&
          Z_SS,&
          RINT_SS,&
          DTAUS_SS,&
          DTAUA_SS,&
          DZ_SS,&
          DALB_SS,&
          NDER,& 
          MAXD)
       !*** input
       integer, intent(in) :: nrt
       real(double), intent(IN), dimension(NSTOKES, NRT) :: z_ss
       real(double), intent(IN), dimension(NRT) :: taua, taus
       integer, intent(IN) :: maxd
       integer, dimension(maxd), intent(in) :: nder
       real(double), intent(in) :: U0, UV
       real(double), dimension(nstokes), intent(in) :: bdrf_ss
       !*** output
       real(double), dimension(nstokes), intent(out) :: RINT_SS
       !intensity vector at top of atmosphere
       !in single scattering approx.
       real(double), dimension(NSTOKES,NRT), intent(out) :: &
            DTAUS_SS, &            !d I_SS / d tau_scat
            DTAUA_SS               !d I_SS  / d tau_abs
       real(double), dimension(NSTOKES, MAXD), intent(out) :: dz_ss
       real(double), dimension(nstokes), intent(out) :: dalb_ss
       !*** local variables      
       real(double), parameter :: EPS = 1.D-10
       real(double), dimension(NSTOKES,NRT) :: DTAUT_SS  !d I_SS/d tau_tot
       real(double), dimension(NSTOKES,NRT) :: &
            DTA,& !d I_SS / d tau_tot
            DW0   !d I_SS / d w0
       integer :: i_st, ik, k, klay
       real(double) :: DTK, dummy, explay, expu0, expuv, tau_above, taut, omega
       !-----------------------------------------------------------------------
       !***  initialization
       RINT_SS = 0.d0            !array
       DTA = 0.D0                !array
       DW0 = 0.D0                !array
       !------------------------------------------------------------------
       !     single scatering contribution of model layers
       !------------------------------------------------------------------

       TAU_ABOVE = 0.d0
       EXPU0     = 1.D0
       EXPUV     = 1.D0

       KLAY = 1
       do K = 1, NRT
          TAUT  = TAUA(K) + TAUS(K)
          OMEGA =  min(1.D0-EPS,TAUS(K)/TAUT)
          EXPLAY = DEXP(-TAUT*(1.D0/UV + 1.D0/U0))
          DUMMY =  U0/(4.D0*PI*(U0+UV))*EXPU0*EXPUV
          do I_ST = 1, NSTOKES
             !     contribution to Stokes parameter per layer,
             !     integation over layers
             RINT_SS(I_ST) = RINT_SS(I_ST)+&
                  OMEGA*Z_SS(I_ST,K)*(1.D0-EXPLAY)*DUMMY

             !     RINT_SS(I_ST) = Z_SS(I_ST)
             !     1. d I / d W0 for tau_tot = const

             DW0(I_ST,K) = DUMMY*Z_SS(I_ST,K)*(1.D0-EXPLAY)

             !     2. d I / d tau_tot for w0 = const

             !     if scattering point in layer K
             DTA(I_ST,K) = OMEGA*DUMMY*Z_SS(I_ST,K)* &
                  EXPLAY*(1.D0/UV + 1.D0/U0)

             !     if scattering point lies below layer K
             DTK = -OMEGA*DUMMY*Z_SS(I_ST,K)*(1.D0-EXPLAY)*&
                  (1.D0/UV + 1.D0/U0)         
             do IK = 1,K-1
                DTA(I_ST,IK) = DTA(I_ST,IK) + DTK
             enddo
          enddo

          if(KLAY .le. MAXD)then
             if(K.eq.NDER(KLAY))then              
                !           derivatives
                DUMMY = DUMMY*(1.D0-EXPLAY)*OMEGA
                DZ_SS(:,klay) = DUMMY
                KLAY = KLAY+1
             endif
          endif

          TAU_ABOVE = TAU_ABOVE+TAUT
          EXPU0  = DEXP(-TAU_ABOVE/U0)
          EXPUV  = DEXP(-TAU_ABOVE/UV)

       enddo

       !-------------------------------------------------------------------
       !     Reflection at the surface or
       !     Fresnel reflection over ocean
       !-------------------------------------------------------------------

       DUMMY =  EXPU0*EXPUV*U0/PI

       do I_ST = 1,NSTOKES
          RINT_SS(I_ST) = RINT_SS(I_ST) + BDRF_SS(I_ST)*DUMMY 
       enddo

       do IK = 1,NRT
          do I_ST = 1,NSTOKES
             DTA(I_ST,IK) = DTA(I_ST,IK) - &
                  BDRF_SS(I_ST)*DUMMY*(1.D0/UV + 1.D0/U0)
          enddo
       enddo

       DALB_SS(1) = DUMMY
       do I_ST = 2, NSTOKES
          DALB_SS(I_ST) = 0.D0
       enddo
       !-------------------------------------------------------------------
       !     Do transformation to d I / d tau_scat for tau = const
       !     and d I / d tau_tot for tau_scat = const
       !-------------------------------------------------------------------

       do K = 1, NRT

          TAUT = TAUA(K) + TAUS(K)
          OMEGA =  min(1.D0-EPS,TAUS(K)/TAUT)
          do I_ST = 1,NSTOKES

             DTAUS_SS(I_ST,K) =  DW0(I_ST,K)/TAUT
             DTAUT_SS(I_ST,K) =  DTA(I_ST,K) - DW0(I_ST,K)*OMEGA/TAUT


          enddo

       enddo

       do K = 1,NRT
          do I_ST = 1,NSTOKES
             DTAUA_SS(I_ST,K) =  DTAUT_SS(I_ST,K)
             DTAUS_SS(I_ST,K) =  DTAUS_SS(I_ST,K) + DTAUT_SS(I_ST,K)
          enddo
       enddo




       return
     end subroutine SINGLE_SCAT
!------------------------------------------------------------------------------  
     subroutine SINGLE_SCAT_FOU_DER(&
          gsf,&
          pm,&
          NF, &
          UV, &   
          U0,&   
          nrt, &
          TAUSCA, &
          TAUABS,&
          BDRF_MS,& 
          RSS_FOU,& 
          DTAU_FOU,&
          DTAUS_FOU,&
          DPHASE_FOU,&
          DALBSS_FOU, &
          NCOEFS,&  
          NDER,&    
          MAXD)
       !*** input
       type(gsf_type), intent(in) :: gsf
       type(phase_mat_type), intent(in) :: pm
       integer, intent(in) :: nrt
       real(double), dimension(nrt), intent(in) :: tausca, tauabs
       real(double), dimension(NSTOKES, NSTOKES, MXHALF+2, MXHALF+2, 0:MAXSTR-1), &
            intent(in) :: bdrf_ms
       integer, intent(in) :: nf
       real(double), intent(in) :: uv, u0
       integer, dimension(nrt), intent(in) :: ncoefs
       integer, intent(in) :: maxd
       integer, dimension(maxd), intent(in) :: nder
       !*** output
       real(double), dimension(nstokes), intent(out) :: rss_fou
       real(double), dimension(NSTOKES, NRT), intent(out) :: dtau_fou, dtaus_fou
       real(double), dimension(NSTOKES, 0:MAXLEG, NPER, MAXD), intent(out) :: dphase_fou
       real(double), dimension(NSTOKES), intent(out) :: dalbss_fou
       !*** local variables
       integer :: i, i_st, k, j_out, ik, i_in, l, NSTR, klay, n
       real(double) :: TAULAY, TAU_ABOVE, OMEGA, EXPU0, EXPLAY, DER2 , DER1, EXPUV, FLX_IN
       real(double), parameter :: eps = 1.D-10     
       real(double), dimension(nstokes) :: rvw, dzss
       !--------------------------------------------------------------------

       NSTR = MXHALF+1

       I_IN = MAXSTR+1
       J_OUT = MAXSTR+2

       do I_ST = 1,NSTOKES
          RSS_FOU(I_ST) = 0.D0
          DALBSS_FOU(I_ST)= 0.D0
       enddo


       do K = 1, MAXD
          do I = 1, NPER
             do L = 0, NCOEFS(K)
                do I_ST = 1, NSTOKES
                   DPHASE_FOU(I_ST, L, I, K) = 0.D0
                enddo
             enddo
          enddo
       enddo

       KLAY = 1
       TAU_ABOVE = 0.D0
       TAULAY = 0.D0
       do K=1,NRT
          TAU_ABOVE = TAU_ABOVE + TAULAY
          TAULAY = TAUSCA(K) + TAUABS(K)
          OMEGA = min(1.D0-EPS,TAUSCA(K)/TAULAY)
          FLX_IN = DEXP(-TAU_ABOVE/U0)
          EXPU0 =  DEXP(-TAU_ABOVE/U0)
          EXPUV = DEXP(-TAU_ABOVE/UV)
          EXPLAY = DEXP(-TAULAY*(1.D0/UV + 1.D0/U0))

          do I_ST = 1, NSTOKES
             RVW(I_ST) = &
                  (OMEGA*U0*EXPU0)/(4.D0*PI*(U0+UV)) * &
                  pm%PM_UP_DN(I_ST,1,NSTR,NSTR,K)* &
                  (1.D0-EXPLAY)

             RSS_FOU(I_ST) = RSS_FOU(I_ST) +RVW(I_ST)*EXPUV

             DZSS(I_ST) = (OMEGA*U0*EXPU0)/(4.D0*PI*(U0+UV)) *&
                  (1.D0-EXPLAY)*EXPUV

             DTAUS_FOU(I_ST,K) = ((U0*EXPU0)/(4.D0*PI*(U0+UV)) *&
                  pm%PM_UP_DN(I_ST,1,NSTR,NSTR,K)* &
                  (1.D0-EXPLAY))*EXPUV/TAULAY

             DTAU_FOU(I_ST,K) =  &
                  ((OMEGA*U0*EXPU0)/(4.D0*PI*(U0+UV)) * &
                  pm%PM_UP_DN(I_ST,1,NSTR,NSTR,K) * &
                  (1.D0/UV + 1.D0/U0)*EXPLAY*EXPUV)

             DER1 = -(OMEGA*U0)/(4.D0*PI*(U0+UV)) * &
                  pm%PM_UP_DN(I_ST,1,NSTR,NSTR,K)* &
                  (1.D0-EXPLAY) * EXPU0/U0

             DER2 =  -EXPUV/UV

             do IK = 1,K-1
                DTAU_FOU(I_ST,IK) = DTAU_FOU(I_ST,IK) + &
                     RVW(I_ST)*DER2 +EXPUV*DER1
             enddo

          enddo

          if(KLAY .le. MAXD) then
             if(K.eq.NDER(KLAY))then

                do L = NF,NCOEFS(K)
                   DPHASE_FOU(1,L,1,KLAY) = DZSS(1)* &
                        gsf%PLM_0(NF,L,J_OUT) * gsf%PLM_0(NF,L,I_IN)               
                   do n = 2, nstokes   
                      if (n==2) DPHASE_FOU(n,L,nper,KLAY) = DZSS(n)* &
                           gsf%PLM_P(NF,L,J_OUT) * gsf%PLM_0(NF,L,I_IN)

                      if (n==3) DPHASE_FOU(n,L,nper,KLAY) = DZSS(n)* &
                           gsf%PLM_M(NF,L,J_OUT) * gsf%PLM_0(NF,L,I_IN)
                   enddo
                enddo

                KLAY = KLAY+1
             endif
          endif

       enddo

       TAU_ABOVE = TAU_ABOVE+TAULAY
       EXPU0 =  DEXP(-TAU_ABOVE/U0)
       EXPUV = DEXP(-TAU_ABOVE/UV)
       DER1 = -EXPU0/U0
       DER2 = -EXPUV/UV

       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       !c     Check!

       if(NF.eq.0)then
          DALBSS_FOU(1) = U0/PI*EXPU0*EXPUV
       endif

       do I_ST = 1,NSTOKES
          RSS_FOU(I_ST) = RSS_FOU(I_ST) +&
                                !     $        BDRF_MS(I_ST,1,MXHALF+1,MXHALF+1,NF)*U0/PI*EXPU0*EXPUV
                                !     Changed on 2010-03-24
               BDRF_MS(I_ST,1,MXHALF+1,MXHALF+2,NF)*U0/PI*EXPU0*EXPUV

          do IK = 1,NRT
             DTAU_FOU(I_ST,IK) = DTAU_FOU(I_ST,IK) +  &
                                !     $           BDRF_MS(I_ST,1,MXHALF+1,MXHALF+1,NF)
                                !     Changed on 2010-03-24
                  BDRF_MS(I_ST,1,MXHALF+1,MXHALF+2,NF) &
                  *U0/PI*(EXPU0*DER2 + EXPUV*DER1)
          enddo
       enddo

       !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


       do K = 1, NRT
          do I_ST = 1, NSTOKES
             DTAU_FOU(I_ST,K) = DTAU_FOU(I_ST,K) +  &
                  DTAUS_FOU(I_ST,K)*(-TAUSCA(K))/(TAUSCA(K)+TAUABS(K))
          enddo
       enddo

       return

     end subroutine SINGLE_SCAT_FOU_DER

!-----------------------------------------------------------------------------

end module single_scat_module
