module perturbation_module
  use header_gsd
  use header_module, only: stopretrieval, writelog
  implicit none

contains
  !**************************************************************************
  subroutine PERT_INTEGRALS(&
       gauss_quad,&
       gsf,&
       pm,&
       ISRC_ADJ,&
       ISRC_FWD,&
       NF,&
       NMAT,&
       NDER,& 
       MAXD,&
       nrt, &
       NCOEFS,& 
       NFILT2, &
       MNK,&
       COSMPHI,&
       SINMPHI,&
       TAUTOT,& 
       EXPMU,&
       EXPU0,&
       EINTU0,&
       EXPUV,&
       EINTUV,&
       OMEGA,&
       RFWD_UP,&
       RFWD_DN,&
       RFWD_VW,&
       RADJ_UP,&
       RADJ_DN,&
       RADJ_VW,&
       U0,&
       UV,&
       SW_SURF,&
       DTAUT, &
       DTAUS,&
       DPHASE,&
       DALB,&
       ierr) 
    !**************************************************************************

    !     Modified by Andre Galli, September 16, 2010
    !     CHANGES: Introduced a failsafe check after encountering few spectra that ran into negative indices:
    !             IF(ISRC_ADJ .LE. 0 .OR. K .LE. 0 .OR. IC .LE. 0 .OR. K .LE. 0 .OR. NFILT2(K).LE. 0) THEN RETURN
    logical, parameter :: PERT_ALB = .TRUE.
    logical, parameter :: PERT_ABS = .TRUE.
    logical, parameter :: PERT_SCA = .TRUE.
    logical, parameter :: PERT_PHA = .TRUE.
    !     inputs
    type(gauss_quad_type), intent(IN) :: gauss_quad
    type(gsf_type), intent(IN) :: gsf
    type(phase_mat_type), intent(IN) :: pm
    integer ::  nrt, &
         ISRC_ADJ,&            !adjoint source index indication
         NF, &                 !Fourier index
         NMAT,&                !matrix dimension for vector RT
         NCOEFS(NRT),&      !number of Fourier terms
         NFILT2(NRT),&      !filter function 2 for auto. cell splitting
         MNK,&                 !number of internal layers
         MAXD,&                !Number of layers
         NDER(MAXD)           !coresponding layer index
    logical :: & 
         SW_SURF              !switch for ground reflection
    real(double) ::  &         
         COSMPHI(0:MAXLEG),&   !cos(M*PHI)
         SINMPHI(0:MAXLEG)    !sin(M*PHI)
    real(double) :: &         !forward and adjoint radiance field
         RFWD_UP(NSTOKES,MXHALF,0:MXK),&     
         RFWD_DN(NSTOKES,MXHALF,0:MXK),&     
         RFWD_VW(NSTOKES,0:MXK),&
         RADJ_UP(NSTOKES,MXHALF,0:MXK),&     
         RADJ_DN(NSTOKES,MXHALF,0:MXK),&
         RADJ_VW(NSTOKES,0:MXK)
    real(double) :: &         !attenuation factors 
         EXPMU(MXHALF,MXK),&   !exp(-tau(k)/mu(i))
         EXPUV(MXK),  &        !exp(-tau(k)/U0)
         EXPU0(MXK), &         !exp(-tau(k)/UV)
         EINTU0(0:MXK),&       !exp(-TAUINT(K)/U0)
         EINTUV(0:MXK)        !exp(-TAUINT(K)/UV) where TAUINT is the
    !                                integrated optical depth
    real(double) :: &
         TAUTOT(MXK), &        !total optical depth
         OMEGA(NRT)        !single scattering albedo
    real(double) :: &
         U0,     &            !cosine of solar zenith angle
         UV                  !cosine of solar zenith angle
    !     output
    real(double)  ::  &                 !derivatives
         DALB(NSTOKES), &               ! d I / d Albedo
         DTAUT(NSTOKES,NRT),&        ! d I / d tau
         DTAUS(NSTOKES,NRT),&        ! d I / d W0
         DPHASE(NSTOKES,0:MAXLEG,NPER,MAXD) !d  I / d PHASE
    integer, intent(out) :: ierr  

    !     internals
    integer :: ISRC_FWD, i_st, ic, imom, i1, i2, klay, kk, ki, k, i, j, j_st, iper
    real(double) ::  FDIR_FWD, FDIR_ADJ, FDIF_ADJ, DUMMY, e1, e2, e3, rdum, sum1, one_exp
    real(double) :: OM_D0M, FDIF_FWD, EXPUU, DELFWD, dum1, AADJ_DN,   AFWD_UP,  DUM3  
    real(double) :: uu, TWO_D0M, sum2, sum3, sum4, ONE_D0M, XI_ADJ_DN, XI_ADJ_UP,  XI_ADJ_VW 
    real(double) ::    XI_FWD_VW , XI_FWD_UP ,    XI_FWD_DN,  DELADJ, dum2, AADJ_UP, AFWD_DN 
    real(double) :: &
         PHA_UP_DN(NSTOKES,NSTOKES,MXHALF+1,MXHALF+1),&
         PHA_UP_UP(NSTOKES,NSTOKES,MXHALF+1,MXHALF  ),&
         PHA_DN_DN(NSTOKES,NSTOKES,MXHALF,  MXHALF+1),&
         PHA_DN_UP(NSTOKES,NSTOKES,MXHALF,  MXHALF  )

    !     dummies

    real(double) :: &
         SUM11(NSTOKES, NSTOKES, MXHALF, MXHALF),&
         SUM12(NSTOKES, NSTOKES, MXHALF, MXHALF),&
         SUM13(NSTOKES, NSTOKES, MXHALF, MXHALF),&
         SUM14(NSTOKES, NSTOKES, MXHALF, MXHALF),&
         SUM21(NSTOKES, MXHALF),&
         SUM22(NSTOKES, MXHALF),&
         SUM31(NSTOKES, MXHALF),&
         SUM32(NSTOKES, MXHALF)

    real(double) :: DELFAC(0:MAXLEG)   !factor (2-delta(m,0)) 
    real(double), parameter :: ONE_4PI = 1/(4.*PI)   ! 1/(4 PI) 

    real(double) :: &
         BFWD(4),&
         BADJ(4),&
         BMIXED(4)

    real(double) :: &
         FAC_E3(NSTOKES)

    real(double),target :: &
         PPER(NSTOKES,NSTOKES,NPER)

    real(double),pointer ::  PPER_tmp(:,:,:)



    !     precalculated single scattering terms ZETA

    real(double) :: &
         ZETA_FWD_UP(MXHALF+1),&
         ZETA_FWD_DN(MXHALF),&
         ZETA_ADJ_UP(MXHALF+1),&
         ZETA_ADJ_DN(MXHALF)

    !     functions related to U0, UV and DG_MU(I)

    real(double) :: &
         U1M(MXHALF),  U1P(MXHALF),  U2M(MXHALF),  U2P(MXHALF),&
         U3M(MXHALF),  U3P(MXHALF),  U4M(MXHALF),  U4P(MXHALF),&
         U5(MXHALF),   U6(MXHALF),   U7,           U8,&
         U9

    real(double) :: &
         K1_FWD,  K2_FWD,&
         K1_ADJ,  K2_ADJ

    real(double) :: &
         LA1, LA2, LA3, LA4,&
         LA5, LA6, LA7, LA8

    real(double) :: &
         EXPU0UV(0:MXK)       !EXP(-TAUINT*(U0+UV)/(UV*U0))

    !-------------------------------------------------------------------------


    DELFAC(NF) = 2.D0         
    if(NF.eq.0)DELFAC(NF)=1.D0

    !     --------------------------------------------------------------------
    !     derivative with respect to Lambertian surface albedo
    !     --------------------------------------------------------------------
    if(PERT_ALB)then


       if(NF.eq.0) then

          if(SW_SURF) then

             FDIF_FWD = 0.
             FDIF_ADJ = 0.
             do J=1,MXHALF
                FDIF_FWD = FDIF_FWD + RFWD_DN(1,J,MNK)*gauss_quad%DG_MU(J)*gauss_quad%DG_WT(J)
                FDIF_ADJ = FDIF_ADJ + RADJ_DN(1,J,MNK)*gauss_quad%DG_WT(J)*gauss_quad%DG_MU(J)
             enddo
             FDIF_FWD =  2*FDIF_FWD 
             FDIF_ADJ =  2*FDIF_ADJ/UV

             FDIR_FWD = U0/PI*EINTU0(MNK)
             if(ISRC_ADJ.eq.1) then
                FDIR_ADJ = EINTUV(MNK)/PI
             else
                FDIR_ADJ = 0.D0
             endif

             DALB(ISRC_ADJ) = PI*(FDIF_FWD*FDIF_ADJ + &
                  FDIF_FWD*FDIR_ADJ + FDIR_FWD*FDIF_ADJ)

          else
             DALB(ISRC_ADJ) = 0.D0
          endif

       endif

    endif

    !     --------------------------------------------------------------------
    !     derivative with respect to absorption optical depth
    !     --------------------------------------------------------------------

    if(PERT_ABS) then

       !     some abbreviations

       BFWD(1) =  COSMPHI(NF)
       BFWD(2) =  COSMPHI(NF)
       BFWD(3) = -SINMPHI(NF)
       BFWD(4) = -SINMPHI(NF)

       if(ISRC_ADJ.eq.1.or.ISRC_ADJ.eq.2)then

          BMIXED(1) =  COSMPHI(NF)
          BMIXED(2) =  COSMPHI(NF)
          BMIXED(3) = -COSMPHI(NF)
          BMIXED(4) = -COSMPHI(NF)

          BADJ(1)   =  COSMPHI(NF)
          BADJ(2)   =  COSMPHI(NF)
          BADJ(3)   =  SINMPHI(NF)
          BADJ(4)   =  SINMPHI(NF)

       else

          BMIXED(1) =  SINMPHI(NF)
          BMIXED(2) =  SINMPHI(NF)
          BMIXED(3) = -SINMPHI(NF)
          BMIXED(4) = -SINMPHI(NF)

          BADJ(1)   =  SINMPHI(NF)
          BADJ(2)   =  SINMPHI(NF)
          BADJ(3)   =  COSMPHI(NF)
          BADJ(4)   =  COSMPHI(NF)

       endif

       if(NF.eq.0) then           
          ONE_D0M   = 2.*PI
          TWO_D0M   = 1.
       else
          ONE_D0M   = 4.*PI
          TWO_D0M   = 2.
       endif

       do I = 1,MXHALF
          U1M(I) = 1./(UV-gauss_quad%DG_MU(I))
          U1P(I) = 1./(UV+gauss_quad%DG_MU(I))
          U2M(I) = U0/(U0-gauss_quad%DG_MU(I))
          U2P(I) = U0/(U0+gauss_quad%DG_MU(I))
          U3M(I) = U2M(I)*gauss_quad%DG_MU(I)
          U3P(I) = U2P(I)*gauss_quad%DG_MU(I)
          U4M(I) = U1M(I)*gauss_quad%DG_MU(I)*UV
          U4P(I) = U1P(I)*gauss_quad%DG_MU(I)*UV
          U5(I)  = U1M(I)*UV**2
          U6(I)  = U0**2/(U0-gauss_quad%DG_MU(I))
       enddo

       U7 = 1./(U0+UV)    
       U8 = U7*U0
       U9 = U8*UV

       !     integrated optical depth

       EXPU0UV(0)  = 1.D0
       do K = 0,MNK
          EXPU0UV(K) = EINTUV(K)*EINTU0(K)             
       enddo

       do K = 1,nrt

          if(NFILT2(K).ne.0)then
             DUMMY  = 0

             IC = 0
             do KK = 1,K-1
                IC = IC + NFILT2(KK)
             enddo

             do KI = 1,NFILT2(K)

                IC = IC + 1

                OM_D0M  = OMEGA(K)/(4*PI)*TWO_D0M

                do I_ST = 1,NMAT

                   DELADJ = 0.
                   DELFWD = 0.
                   if(I_ST.eq.1)    DELADJ = 1.
                   if(I_ST.eq.ISRC_ADJ) DELFWD = 1. 

                   !              -------------------------------------------------------
                   !              multiple scattering terms XI(mu_i) in layer k
                   !              -------------------------------------------------------

                   DUM1=0.D0
                   DUM2=0.D0

                   XI_FWD_VW = 0.D0
                   XI_ADJ_VW = 0.D0

                   do I = 1,MXHALF
                      do J_ST = 1,NMAT

                         AFWD_UP = RFWD_UP(J_ST,I,IC)+RFWD_UP(J_ST,I,IC-1)
                         AFWD_DN = RFWD_DN(J_ST,I,IC)+RFWD_DN(J_ST,I,IC-1)
                         AADJ_UP = RADJ_UP(J_ST,I,IC)+RADJ_UP(J_ST,I,IC-1)
                         AADJ_DN = RADJ_DN(J_ST,I,IC)+RADJ_DN(J_ST,I,IC-1)

                         XI_FWD_VW = XI_FWD_VW + gauss_quad%DG_WT(I)/2 *(&
                              pm%PM_UP_UP(I_ST,J_ST,MXHALF+1,I,K)*AFWD_UP +&
                              pm%PM_UP_DN(I_ST,J_ST,MXHALF+1,I,K)*AFWD_DN)
                         XI_ADJ_VW = XI_ADJ_VW + gauss_quad%DG_WT(I)/2 *(&
                              pm%PM_UP_UP(I_ST,J_ST,MXHALF+2,I,K)*AADJ_UP +&
                              pm%PM_UP_DN(I_ST,J_ST,MXHALF+2,I,K)*AADJ_DN)

                      enddo
                   enddo

                   XI_FWD_VW = OMEGA(K)/2.*XI_FWD_VW
                   XI_ADJ_VW = OMEGA(K)/2.*XI_ADJ_VW

                   !              ------------------------------------------------------
                   !              single scattering terms Zeta(mu_i) in layer k
                   !              ------------------------------------------------------

                   do J = 1,MXHALF
                      ZETA_FWD_UP(J) = ONE_4PI*OMEGA(K)*U2P(J)*&
                           pm%PM_UP_DN(I_ST,ISRC_FWD,J,MXHALF+1,K)
                      ZETA_FWD_DN(J) = ONE_4PI*OMEGA(K)*U2M(J)*&
                           pm%PM_DN_DN(I_ST,ISRC_FWD,J,MXHALF+1,K)
                      ZETA_ADJ_UP(J) = ONE_4PI*OMEGA(K)*U1P(J)*&
                           pm%PM_UP_DN(I_ST,ISRC_ADJ,J,MXHALF+2,K)
                      ZETA_ADJ_DN(J) = ONE_4PI*OMEGA(K)*U1M(J)*&
                           pm%PM_DN_DN(I_ST,ISRC_ADJ,J,MXHALF+2,K)
                   enddo

                   ZETA_ADJ_UP(MXHALF+1)=ONE_4PI*OMEGA(K)*U7*&
                        pm%PM_UP_DN(I_ST,ISRC_ADJ,MXHALF+2,MXHALF+2,K)
                   ZETA_FWD_UP(MXHALF+1)=ONE_4PI*OMEGA(K)*U8*&
                        pm%PM_UP_DN(I_ST,ISRC_FWD,MXHALF+1,MXHALF+1,K)

                   !              direct beam contributions

                   K1_FWD = RFWD_VW(I_ST,IC)*EINTUV(IC)*TAUTOT(IC)/UV 
                   K2_FWD = XI_FWD_VW*EINTUV(IC-1)*(UV-EXPUV(IC)*&
                        (UV+TAUTOT(IC)))/UV  
                   K1_ADJ = RADJ_VW(I_ST,IC)*EINTU0(IC)* TAUTOT(IC)/UV 
                   K2_ADJ = XI_ADJ_VW*EINTU0(IC-1)*(U0-EXPU0(IC)*&
                        (U0+TAUTOT(IC)))/UV

                   DUMMY = DUMMY + TWO_D0M*(&
                        (K1_ADJ+K2_ADJ)*DELADJ*BADJ(I_ST) +&
                        (K1_FWD+K2_FWD)*DELFWD*BFWD(I_ST) )

                   !              diffuse contributions

                   do J=1,MXHALF

                      UU     = gauss_quad%DG_MU(J)
                      EXPUU  = EXPMU(J,IC)
                      ONE_EXP= 1.D0-EXPUU

                      XI_FWD_UP = 0.
                      XI_FWD_DN = 0.
                      XI_ADJ_UP = 0.
                      XI_ADJ_DN = 0.

                      do I = 1,MXHALF
                         do J_ST = 1,NMAT

                            AFWD_UP=RFWD_UP(J_ST,I,IC-1)+RFWD_UP(J_ST,I,IC) 
                            AFWD_DN=RFWD_DN(J_ST,I,IC-1)+RFWD_DN(J_ST,I,IC)
                            AADJ_UP=RADJ_UP(J_ST,I,IC-1)+RADJ_UP(J_ST,I,IC) 
                            AADJ_DN=RADJ_DN(J_ST,I,IC-1)+RADJ_DN(J_ST,I,IC)

                            XI_FWD_DN = XI_FWD_DN + gauss_quad%DG_WT(I)/2.*(&
                                 pm%PM_DN_UP(I_ST,J_ST,J,I,K)*AFWD_UP +&
                                 pm%PM_DN_DN(I_ST,J_ST,J,I,K)*AFWD_DN)
                            XI_FWD_UP = XI_FWD_UP + gauss_quad%DG_WT(I)/2.*(&
                                 pm%PM_UP_UP(I_ST,J_ST,J,I,K)*AFWD_UP +&
                                 pm%PM_UP_DN(I_ST,J_ST,J,I,K)*AFWD_DN)
                            XI_ADJ_DN = XI_ADJ_DN + gauss_quad%DG_WT(I)/2.*(&
                                 pm%PM_DN_UP(I_ST,J_ST,J,I,K)*AADJ_UP +&
                                 pm%PM_DN_DN(I_ST,J_ST,J,I,K)*AADJ_DN)
                            XI_ADJ_UP = XI_ADJ_UP + gauss_quad%DG_WT(I)/2.*(&
                                 pm%PM_UP_UP(I_ST,J_ST,J,I,K)*AADJ_UP +&
                                 pm%PM_UP_DN(I_ST,J_ST,J,I,K)*AADJ_DN)

                         enddo
                      enddo

                      XI_FWD_DN = OMEGA(K)/2.*XI_FWD_DN
                      XI_FWD_UP = OMEGA(K)/2.*XI_FWD_UP
                      XI_ADJ_DN = OMEGA(K)/2.*XI_ADJ_DN
                      XI_ADJ_UP = OMEGA(K)/2.*XI_ADJ_UP

                      LA1 =(RFWD_UP(I_ST,J,IC)  *RADJ_DN(I_ST,J,IC-1)/UV+ &
                           RFWD_DN(I_ST,J,IC-1)*RADJ_UP(I_ST,J,IC)/UV )*&
                           EXPUU*TAUTOT(IC)
                      LA2 =(RFWD_UP(I_ST,J,IC)  *XI_ADJ_DN+ &
                           RADJ_UP(I_ST,J,IC)  *XI_FWD_DN+&
                           RFWD_DN(I_ST,J,IC-1)*XI_ADJ_UP+&
                           RADJ_DN(I_ST,J,IC-1)*XI_FWD_UP)*&
                           (UU-EXPUU*(UU+TAUTOT(IC)))/UV
                      LA3 = RFWD_UP(I_ST,J,IC)* ZETA_ADJ_DN(J)*&
                           EINTUV(IC-1)*(U4M(J)*EXPUV(IC)-&
                           EXPUU*(U4M(J)+TAUTOT(IC))) +&
                           RADJ_UP(I_ST,J,IC)/UV*ZETA_FWD_DN(J)*&
                           EINTU0(IC-1)*(U3M(J)*EXPU0(IC)-&
                           EXPUU*(U3M(J)+TAUTOT(IC)))
                      LA4 =(XI_ADJ_DN*XI_FWD_UP + &
                           XI_ADJ_UP*XI_FWD_DN )/UV *&
                           ((TAUTOT(IC)-2*UU)+(TAUTOT(IC)+2*UU)*EXPUU)
                      LA5 = ZETA_ADJ_DN(J)*XI_FWD_UP* &
                           EINTUV(IC-1)*(UV-UU-(UV+U4M(J))*&
                           EXPUV(IC)+(UU+U4M(J)+TAUTOT(IC))*EXPUU) +&
                           ZETA_FWD_DN(J)*XI_ADJ_UP*&
                           EINTU0(IC-1)*(U0-UU-(U0+U3M(J))*&
                           EXPU0(IC)+(UU+U3M(J)+TAUTOT(IC))*EXPUU)/UV
                      LA6 = RADJ_DN(I_ST,J,IC-1)/UV*ZETA_FWD_UP(J)*&
                           EINTU0(IC-1)*(U3P(J)-EXPUU*EXPU0(IC)*&
                           (U3P(J) + TAUTOT(IC))) + &
                           RFWD_DN(I_ST,J,IC-1)*ZETA_ADJ_UP(J)*&
                           EINTUV(IC-1)*(U4P(J)-EXPUU*EXPUV(IC)*&
                           (U4P(J) + TAUTOT(IC)))
                      LA7 = ZETA_FWD_UP(J)*XI_ADJ_DN* &
                           EINTU0(IC-1)*(U0-U3P(J)-(U0+UU)*EXPU0(IC)+&
                           (UU+U3P(J)+TAUTOT(IC))*EXPUU*EXPU0(IC))/UV +&
                           ZETA_ADJ_UP(J)*XI_FWD_DN*&
                           EINTUV(IC-1)*(UV-U4P(J)-(UV+UU)*EXPUV(IC)+&
                           (UU+U4P(J)+TAUTOT(IC))*EXPUU*EXPUV(IC))
                      LA8 = ZETA_FWD_UP(J)*ZETA_ADJ_DN(J)*&
                           ((U9-U3P(J))*EXPU0UV(IC-1)-(U9+U4M(J))*&
                           EXPU0UV(IC)+(TAUTOT(IC)+U3P(J)+U4M(J))*EXPUU*&
                           EINTU0(IC)*EINTUV(IC-1)) +&
                           ZETA_FWD_DN(J)*ZETA_ADJ_UP(J)*&
                           ((U9-U4P(J))*EXPU0UV(IC-1)-(U9+U3M(J))*&
                           EXPU0UV(IC) +(TAUTOT(IC)+U4P(J)+U3M(J))*EXPUU*&
                           EINTUV(IC)*EINTU0(IC-1))

                      DUMMY = DUMMY + gauss_quad%DG_WT(J)*ONE_D0M*BMIXED(I_ST)*&
                           (LA1+LA2+LA3+LA4+LA5+LA6+LA7+LA8)

                   enddo                                  !stream index J

                enddo                                 !I_ST loop

             enddo                                !internal layer loop KI
             ! Failsafe:
             if(ISRC_ADJ .le. 0 .or. K .le. 0 .or. IC .le. 0 .or. K .le. 0&
                  .or. NFILT2(K).le. 0) then
                ierr = ierr_conv
                !          Pathologic_Flag = 1
                !          write(*,*) 'pathology in PERT_INTEGRALS'
                call writelog('PERT_INTEGRALS: pathology', 6)
                !            goto 301 !Jump to return flag
                RETURN
             endif

             if(ISRC_ADJ.le.2)then
                DTAUT(ISRC_ADJ,K) = DTAUT(ISRC_ADJ,K) -&
                     DUMMY/(TAUTOT(IC)*DFLOAT(NFILT2(K)))
             else
                DTAUT(ISRC_ADJ,K) = DTAUT(ISRC_ADJ,K) +&
                     DUMMY/(TAUTOT(IC)*DFLOAT(NFILT2(K)))
             endif

          endif
       enddo                                  !loop over layer index K

    endif

    !     --------------------------------------------------------------------
    !     derivatives with respect to TAUSCA
    !     --------------------------------------------------------------------

    if(PERT_SCA) then

       FAC_E3(1:nstokes) =  1.D0
       if (nstokes==3) then
          FAC_E3(nstokes) = -FAC_E3(nstokes)
       endif


       do K=1,NRT

          if(NFILT2(K).ne.0)then

             IC = 0
             do KK = 1,K-1
                IC = IC + NFILT2(KK)
             enddo

             RDUM = 0

             do KI = 1,NFILT2(K)

                IC = IC + 1

                E1 = 0
                E2 = 0
                E3 = 0


                DUM1 = PI*TAUTOT(IC)*0.5/UV
                DUM2 = 0.25*EINTUV(IC-1)*(1.D0-EXPUV(IC))
                DUM3 = 0.25*U0/UV*&
                     EINTU0(IC-1)*(1.D0-EXPU0(IC))

                do I1 = 1,MXHALF

                   do I2 = 1,MXHALF

                      SUM1 = 0
                      SUM2 = 0
                      SUM3 = 0
                      SUM4 = 0

                      do J_ST = 1,NMAT
                         do I_ST = 1,NMAT

                            SUM1 = SUM1 + FAC_E3(I_ST)*&
                                 (RADJ_DN(I_ST,I1,IC)  *RFWD_UP(J_ST,I2,IC)+&
                                 RADJ_DN(I_ST,I1,IC-1)*RFWD_UP(J_ST,I2,IC-1))*&
                                 pm%PM_UP_UP(I_ST,J_ST,I1,I2,K)
                            SUM2 = SUM2 + FAC_E3(I_ST)*&
                                 (RADJ_UP(I_ST,I1,IC)  *RFWD_DN(J_ST,I2,IC)+&
                                 RADJ_UP(I_ST,I1,IC-1)*RFWD_DN(J_ST,I2,IC-1))*&
                                 pm%PM_DN_DN(I_ST,J_ST,I1,I2,K)
                            SUM3 = SUM3 + FAC_E3(I_ST)*&
                                 (RADJ_UP(I_ST,I1,IC)  *RFWD_UP(J_ST,I2,IC)+&
                                 RADJ_UP(I_ST,I1,IC-1)*RFWD_UP(J_ST,I2,IC-1))*&
                                 pm%PM_DN_UP(I_ST,J_ST,I1,I2,K)
                            SUM4 = SUM4 + FAC_E3(I_ST)*&
                                 (RADJ_DN(I_ST,I1,IC)  *RFWD_DN(J_ST,I2,IC)+&
                                 RADJ_DN(I_ST,I1,IC-1)*RFWD_DN(J_ST,I2,IC-1))*&
                                 pm%PM_UP_DN(I_ST,J_ST,I1,I2,K)

                         enddo
                      enddo

                      E1 = E1 + &
                           DUM1*gauss_quad%DG_WT(I1)*gauss_quad%DG_WT(I2)*(SUM1+SUM2+SUM3+SUM4)

                   enddo

                   SUM1 = 0.D0
                   SUM2 = 0.D0

                   do I_ST = 1,NMAT

                      SUM1 = SUM1 + FAC_E3(I_ST)*&
                           (RADJ_UP(I_ST,I1,IC)+RADJ_UP(I_ST,I1,IC-1))*&
                           pm%PM_DN_DN(I_ST,ISRC_FWD,I1,MXHALF+1,K)
                      SUM2 = SUM2 +  FAC_E3(I_ST)*&
                           (RADJ_DN(I_ST,I1,IC)+RADJ_DN(I_ST,I1,IC-1))*&
                           pm%PM_UP_DN(I_ST,ISRC_FWD,I1,MXHALF+1,K)

                   enddo

                   E3   = E3 + DUM3*gauss_quad%DG_WT(I1)*(SUM1 + SUM2)

                   SUM1 = 0.D0
                   SUM2 = 0.D0

                   do I_ST = 1,NMAT

                      SUM1 = SUM1 + &
                           pm%PM_UP_UP(ISRC_ADJ,I_ST,MXHALF+1,I1,K)*&
                           (RFWD_UP(I_ST,I1,IC)+RFWD_UP(I_ST,I1,IC-1))
                      SUM2 = SUM2 + &
                           pm%PM_UP_DN(ISRC_ADJ,I_ST,MXHALF+1,I1,K)*&
                           (RFWD_DN(I_ST,I1,IC)+RFWD_DN(I_ST,I1,IC-1))

                   enddo

                   E2   = E2 + DUM2*gauss_quad%DG_WT(I1)*(SUM1 + SUM2)

                enddo

                if(ISRC_ADJ.le.2)then
                   RDUM = RDUM + DELFAC(NF)*(E1+E2+E3)
                else
                   RDUM = RDUM + DELFAC(NF)*(E2-E3-E1)
                endif

             enddo

             if(ISRC_ADJ.le.2)then
                DTAUS(ISRC_ADJ,K) = DTAUS(ISRC_ADJ,K) +& 
                     RDUM*COSMPHI(NF)/(TAUTOT(IC)*DFLOAT(NFILT2(K)))
             else
                DTAUS(ISRC_ADJ,K) = DTAUS(ISRC_ADJ,K) + &
                     RDUM*SINMPHI(NF)/(TAUTOT(IC)*DFLOAT(NFILT2(K)))
             endif

             !         IF(IC.GE.MNK)GOTO 200
          endif
       enddo

       ! 200  CONTINUE

    endif

    !     --------------------------------------------------------------------
    !     derivatives with respect to element a(IMOM,IPER) of the scattering 
    !     phase matrix
    !     --------------------------------------------------------------------

    if(PERT_PHA) then

       !     ---------------------------------------------------
       !     calculation of the perturbation integrals
       !     ---------------------------------------------------
       pper_tmp => pper
       do IPER = 1,NPER         
          do j_st = 1,NSTOKES
             do i_st = 1,NSTOKES
                PPER_tmp(I_ST,J_ST,IPER) = 0.D0
             enddo
          enddo
       enddo

       PPER_tmp(1,1,1) = 1.D0
       if (nstokes > 1) then
          PPER_tmp(2,2,2) = 1.D0
          PPER_tmp(3,3,3) = 1.D0
          PPER_tmp(1,2,4) = 1.D0
          PPER_tmp(2,1,4) = 1.D0
       endif

!       pper(1:nstokes, 1:nstokes, 1:nper) = pper_tmp(1:nstokes, 1:nstokes, 1:nper)

       do KLAY = 1,MAXD

          K    = NDER(KLAY)
          if(NFILT2(K).ne.0)then

             IC = 0
             do KK = 1,K-1
                IC = IC + NFILT2(KK)
             enddo

             RDUM = 0

             do I2 = 1,MXHALF
                do I1 = 1,MXHALF
                   do J_ST = 1,NMAT
                      do I_ST = 1,NMAT
                         SUM11(I_ST,J_ST,I1,I2) = 0.D0
                         SUM12(I_ST,J_ST,I1,I2) = 0.D0
                         SUM13(I_ST,J_ST,I1,I2) = 0.D0
                         SUM14(I_ST,J_ST,I1,I2) = 0.D0
                      enddo
                   enddo
                enddo
             enddo

             do I1 = 1,MXHALF
                do I_ST = 1,NMAT
                   SUM21(I_ST,I1) = 0.D0
                   SUM22(I_ST,I1) = 0.D0
                   SUM31(I_ST,I1) = 0.D0
                   SUM32(I_ST,I1) = 0.D0
                enddo
             enddo

             do KI = 1,NFILT2(K)

                IC = IC + 1

                DUM1 = TAUTOT(IC)
                DUM2 = EINTUV(IC-1)*(1.D0-EXPUV(IC))
                DUM3 = EINTU0(IC-1)*(1.D0-EXPU0(IC))

                do I1 = 1,MXHALF
                   do I2 = 1,MXHALF

                      do J_ST = 1,NMAT
                         do I_ST = 1,NMAT

                            SUM11(I_ST,J_ST,I1,I2) =  &
                                 SUM11(I_ST,J_ST,I1,I2) + DUM1*FAC_E3(I_ST)*&
                                 (RADJ_DN(I_ST,I1,IC)  *RFWD_UP(J_ST,I2,IC)+&
                                 RADJ_DN(I_ST,I1,IC-1)*RFWD_UP(J_ST,I2,IC-1))
                            SUM12(I_ST,J_ST,I1,I2) = &  
                                 SUM12(I_ST,J_ST,I1,I2) + DUM1*FAC_E3(I_ST)*&
                                 (RADJ_UP(I_ST,I1,IC)  *RFWD_DN(J_ST,I2,IC)+&
                                 RADJ_UP(I_ST,I1,IC-1)*RFWD_DN(J_ST,I2,IC-1))
                            SUM13(I_ST,J_ST,I1,I2) =   &
                                 SUM13(I_ST,J_ST,I1,I2) + DUM1*FAC_E3(I_ST)*&
                                 (RADJ_UP(I_ST,I1,IC)  *RFWD_UP(J_ST,I2,IC)+&
                                 RADJ_UP(I_ST,I1,IC-1)*RFWD_UP(J_ST,I2,IC-1))
                            SUM14(I_ST,J_ST,I1,I2) =   &
                                 SUM14(I_ST,J_ST,I1,I2) + DUM1*FAC_E3(I_ST)*&
                                 (RADJ_DN(I_ST,I1,IC)  *RFWD_DN(J_ST,I2,IC)+&
                                 RADJ_DN(I_ST,I1,IC-1)*RFWD_DN(J_ST,I2,IC-1))

                         enddo
                      enddo

                   enddo

                   do I_ST = 1,NMAT

                      SUM21(I_ST,I1) = SUM21(I_ST,I1) + DUM2*&
                           (RFWD_UP(I_ST,I1,IC)+RFWD_UP(I_ST,I1,IC-1))
                      SUM22(I_ST,I1) = SUM22(I_ST,I1) + DUM2*&
                           (RFWD_DN(I_ST,I1,IC)+RFWD_DN(I_ST,I1,IC-1))

                   enddo

                   do I_ST = 1,NMAT

                      SUM31(I_ST,I1) = SUM31(I_ST,I1) + DUM3*FAC_E3(I_ST)*&
                           (RADJ_UP(I_ST,I1,IC)+RADJ_UP(I_ST,I1,IC-1))
                      SUM32(I_ST,I1) = SUM32(I_ST,I1) + DUM3*FAC_E3(I_ST)*&
                           (RADJ_DN(I_ST,I1,IC)+RADJ_DN(I_ST,I1,IC-1))

                   enddo


                enddo               !loop over I1

             enddo                  !loop over KI

             do IPER = 1,NPER         

                do IMOM = NF,NCOEFS(K)

                   call PHI_COEFF_PHASE(&
                        gsf,&
                        NF,&
                        IMOM,&
                        NMAT,&
                        IPER,&   
                        PPER,&
                        PHA_UP_DN,&
                        PHA_UP_UP,&
                        PHA_DN_DN,&
                        PHA_DN_UP)

                   DUM1 = PI*OMEGA(K)*0.5/UV
                   DUM2 = 0.25*OMEGA(K)
                   DUM3 = 0.25*OMEGA(K)*U0/UV

                   E1 = 0.
                   E2 = 0.
                   E3 = 0.

                   do I1 = 1,MXHALF
                      do I2 = 1,MXHALF

                         SUM1 = 0
                         SUM2 = 0
                         SUM3 = 0
                         SUM4 = 0

                         do J_ST = 1,NMAT
                            do I_ST = 1,NMAT

                               SUM1 = SUM1 + SUM11(I_ST,J_ST,I1,I2)*&
                                    PHA_UP_UP(I_ST,J_ST,I1,I2)
                               SUM2 = SUM2 + SUM12(I_ST,J_ST,I1,I2)*&
                                    PHA_DN_DN(I_ST,J_ST,I1,I2)
                               SUM3 = SUM3 + SUM13(I_ST,J_ST,I1,I2)*&
                                    PHA_DN_UP(I_ST,J_ST,I1,I2)
                               SUM4 = SUM4 + SUM14(I_ST,J_ST,I1,I2)*&
                                    PHA_UP_DN(I_ST,J_ST,I1,I2)

                            enddo
                         enddo

                         E1 = E1 + &
                              DUM1*gauss_quad%DG_WT(I1)*gauss_quad%DG_WT(I2)*(SUM1+SUM2+SUM3+SUM4)

                      enddo

                      SUM1 = 0.D0
                      SUM2 = 0.D0

                      do I_ST = 1,NMAT

                         SUM1 = SUM1 + SUM21(I_ST,I1)*&
                              PHA_UP_UP(ISRC_ADJ,I_ST,MXHALF+1,I1)
                         SUM2 = SUM2 + SUM22(I_ST,I1)*&
                              PHA_UP_DN(ISRC_ADJ,I_ST,MXHALF+1,I1)

                      enddo

                      E2 = E2 + DUM2*gauss_quad%DG_WT(I1)*(SUM1 + SUM2)

                      SUM1 = 0.D0
                      SUM2 = 0.D0

                      do I_ST = 1,NMAT

                         SUM1 = SUM1 + SUM31(I_ST,I1)*&
                              PHA_DN_DN(I_ST,ISRC_FWD,I1,MXHALF+1)
                         SUM2 = SUM2 + SUM32(I_ST,I1)*&
                              PHA_UP_DN(I_ST,ISRC_FWD,I1,MXHALF+1)

                      enddo

                      E3 = E3 + DUM3*gauss_quad%DG_WT(I1)*(SUM1 + SUM2)

                   enddo

                   if(ISRC_ADJ.le.2)then
                      RDUM = DELFAC(NF)*(E1+E2+E3)
                   else
                      RDUM = DELFAC(NF)*(E2-E3-E1)
                   endif

                   if(ISRC_ADJ.le.2)then
                      DPHASE(ISRC_ADJ,IMOM,IPER,KLAY) = &
                           DPHASE(ISRC_ADJ,IMOM,IPER,KLAY) + &
                           RDUM*COSMPHI(NF)
                   else
                      DPHASE(ISRC_ADJ,IMOM,IPER,KLAY) = &
                           DPHASE(ISRC_ADJ,IMOM,IPER,KLAY) + &
                           RDUM*SINMPHI(NF)
                   endif

                enddo               !loop over IMOM

             enddo                  !loop over IPER

             !         IF(IC.GE.MNK)GOTO 300

          endif
       enddo                     !loop over KLAY/K 

       ! 300  continue

    endif

301 return
  end subroutine PERT_INTEGRALS

  !*********************************************************************
  subroutine PM_COEFF(&
       nrt, &
       NMAT,  &
       NCOEFS, &
       PLMOM,&
       gsf,&
       pm_in,&
       ierr)
    use header_gsd
    !*                                                                   *
    !*     calculate the Fourier coefficients of the phase matrix        *
    !*     P_IN_OUT =  (-1)^m sum_l=nf^MAXCOEFS P^l_m(-out) S^l P^l_m(-IN) *
    !*     see eg. Eq 3.18 Otto's thesis                                 *
    !*********************************************************************
    !     input
    integer, intent(in) :: nrt, nmat
    real(double) :: &
         PLMOM(NSTOKES,NSTOKES,0:MAXLEG,NRT)
    integer :: &
         NCOEFS(NRT)
    type(gsf_type), intent(IN) :: gsf
    !     output
    integer, intent(out) :: ierr
    type(phase_mat_type), intent(OUT),dimension(0:maxstr) :: pm_in
    !     Fourier coefficients of the phase matrix 
    !     notation PM_DIRout_DIRin(STOKESout,STOKESin,STREAMout,STREAMin,LAYER)
    !     internals
    !     matrix Zm = P^l_m(-out) S^l P^l_m(-in)
    real(double) :: &
         Zm(NSTOKES,NSTOKES)
    integer :: i, i_in, i_st, j_out, k, l, nf, j, j_st
    real(double) :: diff

    !-------------------------------------------------------------------------
    !     
    if (nrt > maxlay) then
       ierr = ierr_var
       call stopretrieval('PM_COEFF: nrt > maxlay')  
       return
    endif

    !     U is -mu from the paper so that P^l_m(-mu) = P^l_m(U)
    !     so U(1:MXHALF) describes downward streams,  U(1+MXHALF:2*MXHALF) the
    !     corresponding upward streams,  U(MAXSTR + 1) the solar beam and
    !     U(MAXSTR + 2) the viewing geometry. U(MAXSTR+3) and U(MAXSTR+4) are
    !     needed for the adjoint formulation.

    !     -------------------------------------------------------------------
    !     initialize for L = NF
    !     -------------------------------------------------------------------

    do NF = 0,MAXSTR-1


       L = NF
       do K = 1,NRT

          do I = 1,MXHALF+2
             do J = 1,MXHALF+2

                !           ----------------------------------------------------------
                !           PM_UP_DN
                !           ----------------------------------------------------------

                I_IN  = I                              !+mu_i
                J_OUT = MXHALF+J                       !-mu_j

                if(I.eq.MXHALF+1) I_IN  = MAXSTR + 1   !-mu_0
                if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   ! mu_v
                if(I.eq.MXHALF+2) I_IN  = MAXSTR + 3   !-mu_v
                if(J.eq.MXHALF+2) J_OUT = MAXSTR + 4   ! mu_o

                call CALC_Zm( nrt, &
                     I_IN,      J_OUT,   K,      L,  &
                     NF,        PLMOM,   gsf,    Zm)

                do J_ST = 1,NMAT
                   do I_ST = 1,NMAT
                      pm_in(NF)%PM_UP_DN(I_ST,J_ST,J,I,K)=Zm(I_ST,J_ST)
                   enddo
                enddo

             enddo
          enddo

          !        ------------------------------------------------------
          !        PM_UP_UP
          !        ------------------------------------------------------

          do I = 1,MXHALF
             do J = 1,MXHALF+2

                I_IN  = MXHALF+I                       ! mu_i
                J_OUT = MXHALF+J                       ! mu_j
                if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   ! mu_v
                if(J.eq.MXHALF+2) J_OUT = MAXSTR + 4   ! mu_o

                call CALC_Zm( nrt, &
                     I_IN,        J_OUT,       K,    L, &
                     NF,          PLMOM,       gsf,  Zm)

                do J_ST = 1,NMAT
                   do I_ST = 1,NMAT
                      pm_in(NF)%PM_UP_UP(I_ST,J_ST,J,I,K)=Zm(I_ST,J_ST)
                   enddo
                enddo

             enddo
          enddo

          !        ---------------------------------------------------------------
          !        PM_DN_DN          -mu_o ->- mu_i
          !        ---------------------------------------------------------------

          I_IN  = MAXSTR+1       !-mu_0

          do I = 1,MXHALF

             J_OUT = I           !-mu_i

             call CALC_Zm( nrt, &
                  I_IN,        J_OUT,   K,      L, &
                  NF,          PLMOM,   gsf,    Zm)

             do J_ST = 1,NMAT
                do I_ST = 1,NMAT
                   pm_in(NF)%PM_DN_DN(I_ST,J_ST,I,MXHALF+1,K) = Zm(I_ST,J_ST)
                enddo
             enddo


          enddo

          I_IN  = MAXSTR+3       !-mu_v

          do I = 1,MXHALF

             J_OUT = I           !-mu_i

             call CALC_Zm( nrt, &
                  I_IN,        J_OUT,   K,      L, &
                  NF,          PLMOM,   gsf,    Zm)

             do J_ST = 1,NMAT
                do I_ST = 1,NMAT
                   pm_in(NF)%PM_DN_DN(I_ST,J_ST,I,MXHALF+2,K) = Zm(I_ST,J_ST)
                enddo
             enddo

          enddo

       enddo

       !     ---------------------------------------------------------------
       !     summation for L = NF+1,MAXCOEFS
       !     ---------------------------------------------------------------

       do K = 1,NRT
          do L = NF+1,NCOEFS(K)


             DIFF = 0.

             do I = 1,MXHALF+2
                do J = 1,MXHALF+2

                   !        -------------------------------------------------------------
                   !        PM_UP_DN
                   !        -------------------------------------------------------------

                   I_IN  = I                              !+mu_i
                   J_OUT = MXHALF+J                       !-mu_j

                   if(I.eq.MXHALF+1) I_IN  = MAXSTR + 1   !-mu_o
                   if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   ! mu_v
                   if(I.eq.MXHALF+2) I_IN  = MAXSTR + 3   !-mu_v
                   if(J.eq.MXHALF+2) J_OUT = MAXSTR + 4   ! mu_o

                   call CALC_Zm(nrt, &
                        I_IN,      J_OUT,   K,      L,  & 
                        NF,        PLMOM,   gsf,    Zm)

                   do J_ST = 1,NMAT
                      do I_ST = 1,NMAT
                         pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K) = &
                              pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K) + Zm(I_ST,J_ST)
                      enddo
                   enddo

                enddo
             enddo

             !        ------------------------------------------------------
             !        PM_UP_UP
             !        ------------------------------------------------------

             do I = 1,MXHALF
                do J = 1,MXHALF+2

                   I_IN  = MXHALF+I                       !-mu_i
                   J_OUT = MXHALF+J                       !-mu_j
                   if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   !-mu_v
                   if(J.eq.MXHALF+2) J_OUT = MAXSTR + 4   !-mu_o

                   call CALC_Zm(nrt, &
                        I_IN,        J_OUT,       K,    L, &
                        NF,          PLMOM,       gsf,  Zm)

                   do J_ST = 1,NMAT
                      do I_ST = 1,NMAT
                         pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K) = &
                              pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K) + Zm(I_ST,J_ST) 
                      enddo
                   enddo

                enddo
             enddo

             !        ---------------------------------------------------------------
             !        PM_DN_DN          - mu_o -> - mu_i
             !        ---------------------------------------------------------------

             I_IN = MAXSTR+1

             do I = 1,MXHALF

                J_OUT = I

                call CALC_Zm( nrt, &
                     I_IN,        J_OUT,   K,     L, &
                     NF,          PLMOM,   gsf,   Zm)

                do J_ST = 1,NMAT
                   do I_ST = 1,NMAT
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,I,MXHALF+1,K) = &
                           pm_in(nf)%PM_DN_DN(I_ST,J_ST,I,MXHALF+1,K) + &
                           Zm(I_ST,J_ST)
                   enddo
                enddo


             enddo

             I_IN = MAXSTR+3

             do I = 1,MXHALF

                J_OUT = I

                call CALC_Zm(nrt, I_IN, J_OUT, K, L, NF, PLMOM, gsf,Zm)

                do J_ST = 1,NMAT
                   do I_ST = 1,NMAT
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,I,MXHALF+2,K) = &
                           pm_in(nf)%PM_DN_DN(I_ST,J_ST,I,MXHALF+2,K) + &
                           Zm(I_ST,J_ST)
                   enddo
                enddo

             enddo

          enddo

       enddo

       !------------------------------------------------------------------
       !     use symmetry relations to define the missing element of 
       !     PM_DN_DN and PM_UP_DN
       !------------------------------------------------------------------

       do K = 1,NRT
          do I = 1,MXHALF
             do J = 1,MXHALF

                do J_ST = 1, min(2,nstokes)
                   do I_ST = 1, min(2, nstokes)                     
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,J,I,K) = &
                           pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K)
                      pm_in(nf)%PM_DN_UP(I_ST,J_ST,J,I,K) = &
                           pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K)
                   enddo
                enddo
                if (nstokes>1) then
                do I_ST = 3,NMAT
                   do J_ST = 3,NMAT
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,J,I,K) = &
                           pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K)
                      pm_in(nf)%PM_DN_UP(I_ST,J_ST,J,I,K) = &
                           pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K)
                   enddo
                enddo


                DO J_ST = 1,2
                   DO I_ST = 3,NMAT
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,J,I,K) = &
                           -pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K)
                      pm_in(nf)%PM_DN_UP(I_ST,J_ST,J,I,K) =&
                           -pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K)
                   ENDDO
                ENDDO

                DO J_ST = 3,NMAT
                   DO I_ST = 1,2
                      pm_in(nf)%PM_DN_DN(I_ST,J_ST,J,I,K) =&
                           -pm_in(nf)%PM_UP_UP(I_ST,J_ST,J,I,K)
                      pm_in(nf)%PM_DN_UP(I_ST,J_ST,J,I,K) =&
                           -pm_in(nf)%PM_UP_DN(I_ST,J_ST,J,I,K)
                   ENDDO
                ENDDO
                endif


             enddo
          enddo
       enddo

    enddo

    return
  end subroutine PM_COEFF
  !*********************************************************************
  !> calculate the Fourier coefficients of the phase matrix    
  !! P_IN_OUT =  (-1)^m sum_l=nf^NCOEFS P^l_m(-out) S^l P^l_m(-IN) 
  !! see eg. Eq 3.18 Otto's thesis  
  subroutine PHI_COEFF_PHASE(&
       gsf,&
       NF,& 
       IMOM,&
       NMAT,&
       IPER,&
       PPER,&
       P_UP_DN,&
       P_UP_UP,&
       P_DN_DN,&
       P_DN_UP)
    !*** input
    integer, intent(in) :: imom, nf, nmat, iper
    type(gsf_type), intent(IN) :: gsf
    real(double), intent(in) :: PPER(NSTOKES, NSTOKES, NPER)
    !*** output
    !c     coefficients of the phase matrix 
    !c     notation P_DIRout_DIRin(STOKESout,STOKESin,STREAMout,STREAMin,LAYER)
    real(double), intent(out) :: &
         P_UP_DN(NSTOKES,NSTOKES,MXHALF+1,MXHALF+1),&
         P_UP_UP(NSTOKES,NSTOKES,MXHALF+1,MXHALF  ),&
         P_DN_DN(NSTOKES,NSTOKES,MXHALF,  MXHALF+1),&
         P_DN_UP(NSTOKES,NSTOKES,MXHALF,  MXHALF  )
    !     internals
    !     matrix Zm = P^l_m(-out) S^l P^l_m(-in)
    real(double) :: Zm(NSTOKES, NSTOKES)
    integer :: i_st, j, j_st, l, j_out, i_in, i
    !c-------------------------------------------------------------------------
    !c     
    !c     U is -mu from the paper so that P^l_m(-mu) = P^l_m(U)
    !c     so U(1:MXHALF) describes downward streams,  U(1+MXHALF:2*MXHALF) the
    !c     corresponding upward streams,  U(MAXSTR + 1) the solar beam and
    !c     U(MAXSTR + 2) the viewing geometry. U(MAXSTR+3) and U(MAXSTR+4) are
    !c     needed for the adjoint formulation.

    L = IMOM
    do I = 1, MXHALF+1
       do J = 1, MXHALF+1
          !c        ----------------------------------------------------------
          !c        P_UP_DN
          !c        ----------------------------------------------------------
          I_IN  = I                              !+mu_i
          J_OUT = MXHALF+J                       !-mu_j
          if(I.eq.MXHALF+1) I_IN  = MAXSTR + 1   !-mu_0
          if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   ! mu_v
          call CALC_Zm_PHA(&
               I_IN, &
               J_OUT, &
               IMOM, &
               IPER, &
               NF, &   
               PPER, &  
               gsf, &
               Zm)         
          do I_ST = 1, NMAT
             do J_ST = 1, NMAT
                P_UP_DN(I_ST,J_ST,J,I) = Zm(I_ST,J_ST)
             enddo
          enddo
       enddo
    enddo

    !c     ------------------------------------------------------
    !c     P_UP_UP
    !c     ------------------------------------------------------

    do I = 1,MXHALF
       do J = 1,MXHALF+1
          I_IN  = MXHALF+I                       ! mu_i
          J_OUT = MXHALF+J                       ! mu_j
          if(J.eq.MXHALF+1) J_OUT = MAXSTR + 2   ! mu_v
          call CALC_Zm_PHA(&
               I_IN, &
               J_OUT, &
               IMOM, &
               IPER,&
               NF, &
               PPER, &
               gsf, &
               Zm)
          do I_ST = 1, NMAT
             do J_ST = 1, NMAT
                P_UP_UP(I_ST,J_ST,J,I) = Zm(I_ST,J_ST)
             enddo
          enddo
       enddo
    enddo

    !c     ---------------------------------------------------------------
    !c     P_DN_DN          -mu_o ->- mu_i
    !c     ---------------------------------------------------------------

    I_IN  = MAXSTR+1       !-mu_0

    do I = 1,MXHALF
       J_OUT = I           !-mu_i
       call CALC_Zm_PHA(&
            I_IN,        J_OUT,   IMOM,  IPER,& 
            NF,          PPER,    gsf,   Zm)      

       do I_ST = 1, NMAT
          do J_ST = 1,NMAT
             P_DN_DN(I_ST,J_ST,I,MXHALF+1) = Zm(I_ST,J_ST)
          enddo
       enddo
    enddo

    !c------------------------------------------------------------------
    !c     use symmetry relations to define the missing element of 
    !c     P_DN_DN and P_UP_DN
    !c------------------------------------------------------------------
    do I = 1,MXHALF
       do J = 1, MXHALF 
          do J_ST = 1, NMAT           !abutz: was I_ST = 1, 1
             do I_ST = 1, NMAT        !abutz: was J_ST = 1, 1             
                P_DN_DN(I_ST,J_ST,J,I) = P_UP_UP(I_ST,J_ST,J,I)
                P_DN_UP(I_ST,J_ST,J,I) = P_UP_DN(I_ST,J_ST,J,I)
             enddo
          enddo
       enddo
    enddo

    return
  end subroutine PHI_COEFF_PHASE
  !---------------------------------------------------------------------------
  subroutine CALC_Zm_PHA(&
       I_IN, &
       J_OUT, &
       L, &
       IPER, &
       NF, &
       PPER_in, &
       gsf, &
       Zm_out)      
    !*** input
    type(gsf_type), intent(IN) :: gsf
    integer, intent(in) :: i_in, l, j_out, nf, iper
    real(double), target, intent(in) ::  PPER_in(NSTOKES, NSTOKES, NPER)
    !*** output
    real(double), target, intent(out) :: Zm_out(NSTOKES, NSTOKES)
    !***
    real(double),pointer :: PPER(:,:,:), Zm(:,:)    
    !---------------------------------------------------------------------------

    pper => pper_in
    zm => zm_out

    Zm(1,1) = gsf%PLM_0(NF,L,J_OUT) *PPER(1,1,iper)* gsf%PLM_0(NF,L,I_IN)    
    if (nstokes>1) then
       Zm(2,1) =  gsf%PLM_P(NF,L,J_OUT) *PPER(2,1,iper)*  gsf%PLM_0(NF,L,I_IN)
       Zm(3,1) =  gsf%PLM_M(NF,L,J_OUT) *PPER(2,1,iper)*  gsf%PLM_0(NF,L,I_IN)

       Zm(1,2) =  gsf%PLM_0(NF,L,J_OUT) *PPER(1,2,iper)*  gsf%PLM_P(NF,L,I_IN)
       Zm(2,2) =  gsf%PLM_P(NF,L,J_OUT) *PPER(2,2,iper)*  gsf%PLM_P(NF,L,I_IN)+&
            gsf%PLM_M(NF,L,J_OUT) *PPER(3,3,iper)*  gsf%PLM_M(NF,L,I_IN)
       Zm(3,2) =  gsf%PLM_M(NF,L,J_OUT) *PPER(2,2,iper)*  gsf%PLM_P(NF,L,I_IN)+&
            gsf%PLM_P(NF,L,J_OUT) *PPER(3,3,iper)*  gsf%PLM_M(NF,L,I_IN)

       Zm(1,3) =  gsf%PLM_0(NF,L,J_OUT) *PPER(1,2,iper)*  gsf%PLM_M(NF,L,I_IN)
       Zm(2,3) =  gsf%PLM_P(NF,L,J_OUT) *PPER(2,2,iper)*  gsf%PLM_M(NF,L,I_IN)+&
            gsf%PLM_M(NF,L,J_OUT) *PPER(3,3,iper)*  gsf%PLM_P(NF,L,I_IN)
       Zm(3,3) =  gsf%PLM_M(NF,L,J_OUT) *PPER(2,2,iper)*  gsf%PLM_M(NF,L,I_IN)+&
            gsf%PLM_P(NF,L,J_OUT) *PPER(3,3,iper)*  gsf%PLM_P(NF,L,I_IN)
    endif
    

!    zm_out(1:nstokes,1:nstokes) = zm(1:nstokes,1:nstokes)

    return
  end subroutine CALC_Zm_PHA

  !----------------------------------------------------------------------------

  subroutine CALC_Zm(&
       nrt, &
       I_IN, &
       J_OUT,&
       K, &
       L, &  
       NF, &
       PLMOM_in, &
       gsf, &
       Zm_out)
    !*** input
    integer, intent(in) :: nrt, i_in, j_out, k, l, nf 
    real(double), target, intent(in) :: &
         PLMOM_in(NSTOKES,NSTOKES,0:MAXLEG,NRT)
    type(gsf_type), intent(in) :: gsf
    !*** output
    real(double), target, intent(out) :: Zm_out(NSTOKES,NSTOKES)
    !***
    integer :: n1, n2
    real(double), pointer :: PLMOM(:,:,:,:), Zm(:,:)
    !---------------------------------------------------------------------------

    plmom => plmom_in
    zm => zm_out


    Zm(1,1) = gsf%PLM_0(NF,L,J_OUT) *PLMOM(1,1,l,k)* gsf%PLM_0(NF,L,I_IN)
    if (nstokes>1) then
       Zm(2,1) = gsf%PLM_P(NF,L,J_OUT) *PLMOM(2,1,l,k)* gsf%PLM_0(NF,L,I_IN)
       Zm(3,1) = gsf%PLM_M(NF,L,J_OUT) *PLMOM(2,1,l,k)* gsf%PLM_0(NF,L,I_IN)

       Zm(1,2) = gsf%PLM_0(NF,L,J_OUT) *PLMOM(1,2,l,k)* gsf%PLM_P(NF,L,I_IN)
       Zm(2,2) = gsf%PLM_P(NF,L,J_OUT) *PLMOM(2,2,l,k)* gsf%PLM_P(NF,L,I_IN)+&
            gsf%PLM_M(NF,L,J_OUT) *PLMOM(3,3,l,k)* gsf%PLM_M(NF,L,I_IN)
       Zm(3,2) = gsf%PLM_M(NF,L,J_OUT) *PLMOM(2,2,l,k)* gsf%PLM_P(NF,L,I_IN)+&
            gsf%PLM_P(NF,L,J_OUT) *PLMOM(3,3,l,k)* gsf%PLM_M(NF,L,I_IN)

       Zm(1,3) = gsf%PLM_0(NF,L,J_OUT) *PLMOM(1,2,l,k)* gsf%PLM_M(NF,L,I_IN)
       Zm(2,3) = gsf%PLM_P(NF,L,J_OUT) *PLMOM(2,2,l,k)* gsf%PLM_M(NF,L,I_IN)+&
            gsf%PLM_M(NF,L,J_OUT) *PLMOM(3,3,l,k)* gsf%PLM_P(NF,L,I_IN)
       Zm(3,3) = gsf%PLM_M(NF,L,J_OUT) *PLMOM(2,2,l,k)* gsf%PLM_M(NF,L,I_IN)+&
            gsf%PLM_P(NF,L,J_OUT) *PLMOM(3,3,l,k)* gsf%PLM_P(NF,L,I_IN)
    endif


    return

  end subroutine CALC_Zm

  !----------------------------------------------------------------------------

end module perturbation_module
