module gaus_mat_module
   use header_gsd
   implicit none
   contains 

!***************************************************************************
   subroutine GAUITR_MAT(&
              nrt, &
           NF,&
           MNK,&
           IFWD_ADJ,&  
           I_SRC,&
           NMAT, &
           KFILT, &  
           EPST,&   
           UBEAM,&
           FLXBEAM,&
           SW_SURF,&
           PM,&
           gauss_quad,&
           BDRF_MS, & 
           OMEGA, &
           EXPMU,& 
           EXPINT,&  
           RSS_UP,&
           RSS_DN,& 
           RAD_UP, &
           RAD_DN, &  
           NCOEFS)

!     inputs
      integer :: nrt, &
          NF,&                       !index of Fourier component
          MNK,&                      !number of internal layers
          KFILT(MXK)
      type(phase_mat_type),intent(IN) :: pm
      type(gauss_quad_type), intent(IN) :: gauss_quad
      integer, dimension(nrt) :: ncoefs
      real(double) :: &
          EPST                      !determines trunction of Gauss-Seidel 
                                     !iteration series
      real(double) :: &
          UBEAM,&                     !cosine of solar zenith angle 
          BDRF_MS(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,0:MAXSTR-1),& !BDRF, Fourier decomposition
          FLXBEAM,&                   !direct downward flux at the surface
          OMEGA(NRT), &            !single scattering albedo
          EXPMU(MXHALF,MXK),&         !attenuation fact. EXP(-tau(k)/UU(i))
          EXPINT(0:MXK),&             !attenuation fact. EXP(-TAUINT(K)/UBEAM)
          RSS_UP(NSTOKES,MXHALF,MXK),&!single scattering contribution (upward)
          RSS_DN(NSTOKES,MXHALF,MXK) !(downward)          
! internals
      integer :: i, i_mat, i_mat_start, ic, ik_form, IFWD_ADJ, I_SRC, j_mat, NMAT, MAXMAT, mdb
      integer :: i_st, ia, ia_start, ik, j, j_st, ja, k, ja_start, itr
      real(double) :: DIFF, DUMMY, FACT,  FLX0, MDB2, MDB3, RAVG_UP, RAVG_DN
      real(double) :: &
          GAMFAC2_DN((MNK+1)*MXHALF*NMAT),&
          GAMFAC2_UP((MNK+1)*MXHALF*NMAT),&
          GAMSS_UP((MNK+1)*MXHALF*NMAT),&
          GAMSS_DN((MNK+1)*MXHALF*NMAT),&
          GAMMU((MNK+1)*MXHALF*NMAT),&
          GAMFAC((MNK+1)*MXHALF*NMAT)

      real(double) :: &
          A_UP(2*MXHALF*NSTOKES,NRT*MXHALF*NSTOKES),&
          A_DN(2*MXHALF*NSTOKES,NRT*MXHALF*NSTOKES)
!     out- and inputs
      real(double) :: &
          RAD_UP(NSTOKES,MXHALF,0:MXK),&     !Fourier comp. of upward intensity  
          RAD_DN(NSTOKES,MXHALF,0:MXK)      !Fourier comp. of downward intensity
!     multiple scattering source function at the layer boundaries beside
!     a factor 0.5*OMEGA(K) that has to multiply later.
      real(double) :: &
          FAC2_DN(NSTOKES,MXHALF,MXK),&
          FAC2_UP(NSTOKES,MXHALF,MXK)
      real(double) :: &
          GAMMA_UP((MNK+1)*MXHALF*NMAT),&
          GAMMA_DN((MNK+1)*MXHALF*NMAT)
      real(double) :: &
          RCON_TRA(MXHALF), &   !reflected and transmited component from
          RCON_RFL(MXHALF)     !previous iteration step (conv. check)
      real(double), parameter :: EPS = 1.D-20
      real(double) :: &              !precalculated coefficients
          FAC(MXHALF,MXK), &          !which are not affected by the iteration
          FLXFAC(MXHALF)
      logical :: SW_SURF

!------------------------------------------------------------------------

      MDB = MXHALF*NMAT
      MDB2 = 2*MDB
      MDB3 = 3*MDB
      MAXMAT = (MNK+1)*MDB
      
!     make vector gamma

      do I = 1,MXHALF
         do K = 1,MNK
            IK = KFILT(K)
            FAC(I,K) = 0.5 *  OMEGA(IK) * (1-EXPMU(I,K))
         enddo
         FLXFAC(I) = 2.D0*gauss_quad%DG_MU(I)*gauss_quad%DG_WT(I)
      enddo
      FLX0 = FLXBEAM/PI

      do K=1,MNK
         IK = KFILT(K)
         do I = 1,MXHALF
            DUMMY = (1.D0 - EXPMU(I,K))*OMEGA(IK)*gauss_quad%DG_WT(I)/4.
            do I_ST = 1,NMAT
               FAC2_DN(I_ST,I,K) = &
                   1./(1-pm%PM_DN_DN(I_ST,I_ST,I,I,IK)*DUMMY)
               FAC2_UP(I_ST,I,K) = &
                   1./(1-pm%PM_UP_UP(I_ST,I_ST,I,I,IK)*DUMMY)
            enddo
         enddo
      enddo


      IC = 1
      do K=0,MNK
      do I=1,MXHALF
      do I_ST=1,NMAT
         GAMMA_DN(IC) = RAD_DN(I_ST,I,K)
         GAMMA_UP(IC) = RAD_UP(I_ST,I,K)
         IC = IC+1
      enddo
      enddo
      enddo

      IC = MDB+1
      do K=1,MNK
      do I=1,MXHALF
      do I_ST = 1,NMAT
         GAMFAC2_DN(IC) = FAC2_DN(I_ST,I,K)
         GAMFAC2_UP(IC) = FAC2_UP(I_ST,I,K)
         GAMSS_UP(IC)   = RSS_UP(I_ST,I,K)
         GAMSS_DN(IC)   = RSS_DN(I_ST,I,K)
         GAMMU(IC)      = EXPMU(I,K)
         GAMFAC(IC)     = FAC(I,K)
         IC = IC+1
      enddo
      enddo
      enddo

      IA = 0
      JA_START = 0
      do K = 1,NRT
      do I = 1,MXHALF
      do I_ST = 1,NMAT
         IA = IA +1
         JA = JA_START
         do J=1,MXHALF
            FACT = gauss_quad%DG_WT(J)/2.
            do J_ST = 1,NMAT
               JA = JA+1
               A_UP(JA,IA)     = FACT*pm%PM_UP_DN(I_ST,J_ST,I,J,K)
               A_UP(JA+MDB,IA) = FACT*pm%PM_UP_UP(I_ST,J_ST,I,J,K)
               A_DN(JA,IA)     = FACT*pm%PM_DN_DN(I_ST,J_ST,I,J,K)
               A_DN(JA+MDB,IA) = FACT*pm%PM_DN_UP(I_ST,J_ST,I,J,K)
            enddo
         enddo
      enddo
      enddo
      enddo

!     store in RCONV intensities at top and bottom of the atmosphere 
!     for convergence checking   

      do I=1,MXHALF       
         RCON_TRA(I) = RAD_DN(1,I,MNK) 
         RCON_RFL(I) = RAD_UP(1,I,0)
      enddo
  
      do ITR=1,MAXITR

         IA = MDB
         JA_START = MDB

         I_MAT_START = 0
         IK_FORM =1
         do K = 1,MNK
            IK = KFILT(K)
            
            I_MAT_START = I_MAT_START + (IK-IK_FORM)*MDB
            I_MAT = I_MAT_START
            do I = 1,MDB
               IA = IA+1
               I_MAT = I_MAT+1
               GAMMA_DN(IA) = 0.
               DUMMY = 0.

               JA = JA_START
               J_MAT = 0

               if(NF.le.ncoefs(ik))then
!               IF(ncoefs(ik).GT.2)THEN

               do J=1,MDB
                  JA = JA+1
                  J_MAT = J_MAT+1

                  RAVG_UP = GAMMA_UP(JA-MDB) + GAMMA_UP(JA) 
                  RAVG_DN = GAMMA_DN(JA-MDB) + GAMMA_DN(JA)
                  
                  DUMMY = DUMMY +&
                      A_DN(J_MAT,I_MAT)*RAVG_DN +&
                      A_DN(J_MAT+MDB,I_MAT)*RAVG_UP
                  
               enddo
               endif


               GAMMA_DN(IA) = GAMFAC2_DN(IA)*&
                   (GAMMA_DN(IA-MDB) * GAMMU(IA) +&
                   DUMMY * GAMFAC(IA) + GAMSS_DN(IA) )

            enddo

            JA_START = JA_START+MDB
            IK_FORM = IK
            
         enddo

!-------ground reflection --------------------------------------------

         IA_START = MAXMAT-MDB+1
         do IA = IA_START,MAXMAT
            GAMMA_UP(IA) = 0.
         enddo
         
         IA = MAXMAT-MDB

         
         if(SW_SURF)then  

            do I=1,MXHALF
            do I_ST = 1,NMAT
               IA = IA+1
               JA = MAXMAT-MDB+1
                  
               GAMMA_UP(IA) = UBEAM*EXPINT(MNK)/PI* &
                             BDRF_MS(I_ST,I_SRC,MXHALF+IFWD_ADJ,I,NF)

               do J=1,MXHALF
               do J_ST = 1,NMAT
                  GAMMA_UP(IA) = GAMMA_UP(IA) + GAMMA_DN(JA) * &
                      (2.*gauss_quad%DG_WT(J)*gauss_quad%DG_MU(J))*BDRF_MS(I_ST,J_ST,J,I,NF)
                  JA = JA+1
               enddo
               enddo

            enddo
            enddo

         endif 

!----------------------------------------------------------------------------

         IA = MAXMAT-MDB2
         JA_START = MAXMAT-MDB2
         
         I_MAT_START = NRT*NMAT*MXHALF-MDB
!c         IK_FORM = KFILT(MNK-1)   !bug changed 5.11.03
         IK_FORM = KFILT(MNK)

         do K=MNK-1,0,-1
          
            IK = KFILT(K+1)
         
!     calculate  averaged source function for layer K

            I_MAT_START = I_MAT_START - (IK_FORM-IK)*MDB

            I_MAT = I_MAT_START
            do I = 1,MDB
               IA = IA+1
               I_MAT = I_MAT+1
               GAMMA_UP(IA) = 0.
               DUMMY = 0.

               JA = JA_START
               J_MAT = 0

               
               if(nf.le.ncoefs(ik))then
!               IF(ncoefs(ik).GT.2)THEN
               do J=1,MDB
                  J_MAT = J_MAT +1
                  JA = JA+1
                  
                  RAVG_UP = GAMMA_UP(JA+MDB) + GAMMA_UP(JA)
                  
                  RAVG_DN =  GAMMA_DN(JA+MDB) + GAMMA_DN(JA)

                  DUMMY = DUMMY +&
                      A_UP(J_MAT,I_MAT)*RAVG_DN +&
                      A_UP(J_MAT+MDB,I_MAT)*RAVG_UP
                  
               enddo
               endif

               GAMMA_UP(IA) = &
                   GAMFAC2_UP(IA+MDB)*&
                   (GAMMA_UP(IA+MDB) * GAMMU(IA+MDB)    +&
                   DUMMY * GAMFAC(IA+MDB)      +&
                   GAMSS_UP(IA+MDB) )
                  
            enddo
            
            IA = IA-MDB2
            JA_START = JA_START-MDB
            
            IK_FORM = IK

         enddo

!-------check for convergence -------------------------------------------
         
         IA = MAXMAT-MDB+1
         JA = 1
         do I=1,MXHALF                               
            
            if(abs(RCON_TRA(I)).gt.EPS)then
               DIFF=abs( (GAMMA_DN(IA)-RCON_TRA(I))/RCON_TRA(I) ) 
               if(DIFF.gt.EPST) goto 10                 
            endif
            IA = IA + NMAT

            if(abs(RCON_RFL(I)).gt.EPS) then   
               DIFF=abs( (GAMMA_UP(JA)-RCON_RFL(I))/RCON_RFL(I) )      
               if(DIFF.gt.EPST) goto 10                 
            endif
            JA = JA+NMAT
            
         enddo
         
         goto 100               !jump out the iteration loop
         

 10      IA = MAXMAT-MDB+1
         JA=1
         do I=1,MXHALF          !reset RCON arrays
            RCON_TRA(I) = GAMMA_DN(IA)             
            RCON_RFL(I) = GAMMA_UP(JA)
            IA = IA+NMAT
            JA = JA+NMAT
                  
         enddo


      enddo

!     convergence has failed            

      ! Remove write-statement for operational RemoTeC version.
      !if(ITR.eq.MAXITR) then
      !   write(*,*)'convergence criterion is not satisfied after',&
      !          MAXITR,'  iterations.'       
      !endif

 100  continue
      
      IC = 1
      do K=0,MNK
      do I=1,MXHALF
      do I_ST=1,NMAT
         RAD_DN(I_ST,I,K) = GAMMA_DN(IC) 
         RAD_UP(I_ST,I,K) = GAMMA_UP(IC)
         IC = IC+1
      enddo
      enddo
      enddo
         
      return
    end subroutine GAUITR_MAT

!**************************************************************************
      subroutine SSC_LAYER_POL(&
           nrt, &
           gauss_quad,&
           pm,&
           IFWD_ADJ,&
           MNK,&
           OMEGA,&
           FLXBEAM,&
           UBEAM,&
           EXPMU,&
           EXPBEAM,&  
           EXPINT, &
           RSS_UP,&
           RSS_DN,&
           I_SRC,&
           NFILT1)
        use header_gsd
!*                                                                        *
!*     integration of single scattering source function over              *
!*     model layers                                                       *
!**************************************************************************                          
      implicit none 
!     inputs      
      integer, intent(in) :: nrt
      type(phase_mat_type),intent(IN) :: pm
      type(gauss_quad_type), intent(IN) :: gauss_quad
      integer :: NFILT1(MXK)
      real(double) :: &
          OMEGA(NRT),  &           !single scattering albedo
          UBEAM,  &                   !cosine of beam zenith angle
          EXPMU(MXHALF,MXK),  &       !attenuation factors EXP(-tau(k)/mu(i))
          EXPBEAM(MXK), &             !attenuation factors EXP(-tau(k)/UBEAM)
          EXPINT(0:MXK)             !exp(-TAUINT(K)/UBEAM)  
!     outputs

      real(double) ::  &
          RSS_UP(NSTOKES,MXHALF,MXK),  &      !single scattering source (upward)
          RSS_DN(NSTOKES,MXHALF,MXK)        !(downward)      
      real(double) ::  &
          FLXBEAM                           !downward flux of the beam at the surface
!     internals
      integer :: ifwd_adj, mnk, i_src, i, i_st, ic, ik, k, nstr
      real(double) ::  &
          FRACMIN(MXHALF),  &                 !UBEAM/(UBEAM-DG_MU)
          FRACPLU(MXHALF)                   !UBEAM/(UBEAM+DG_MU)
      real(double) :: one_4pi
!----------------------------------------------------------------
      IC = 0

      NSTR = MXHALF + IFWD_ADJ

      ONE_4PI=1./(4.*PI)  

      do I = 1,MXHALF
         FRACMIN(I) = ONE_4PI*UBEAM/(UBEAM - gauss_quad%DG_MU(I))
         FRACPLU(I) = ONE_4PI*UBEAM/(UBEAM + gauss_quad%DG_MU(I))
      enddo

      do I_ST = 1,NSTOKES
         do K = 1,MNK
            IK = NFILT1(K)
            do I = 1,MXHALF
               RSS_DN(I_ST,I,K) = OMEGA(IK)*&
                   pm%PM_DN_DN(I_ST,I_SRC,I,NSTR,IK)*FRACMIN(I)*&
                   EXPINT(K-1)*(EXPBEAM(K) -EXPMU(I,K))
               RSS_UP(I_ST,I,K) = OMEGA(IK)* &
                   pm%PM_UP_DN(I_ST,I_SRC,I,NSTR,IK)*FRACPLU(I)* &
                  EXPINT(K-1)*(1-EXPBEAM(K)*EXPMU(I,K))
            enddo
         enddo
      enddo

      FLXBEAM=0.
      if(I_SRC.eq.1) FLXBEAM = UBEAM*EXPINT(MNK)

      return                                                                  
    end subroutine SSC_LAYER_POL
                                                             
!**************************************************************************
      subroutine GAUSSEID_INIT_SSC_POL(&
          NF,        MNK,     IFWD_ADJ, I_SRC, &
          EXPMU,     FLXDIR,  SW_SURF,  BDRF_MS, &
          RSS_UP,    RSS_DN,  RAD_UP,   RAD_DN)
        use header_gsd
!*                                                                        *
!*     calculate the initialization of the Gauss-Seidel iteration         *
!*     using the single scattering approximation                          *
!**************************************************************************
      implicit none 
!     inputs
      real(double) ::  &
          EXPMU(MXHALF,MXK), &          !attenuation factor
          FLXDIR, &                     !direct downward flux at the surface
          BDRF_MS(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,0:MAXSTR-1), & 
!                                      !BDRF, Fourier decomposition
          RSS_UP(NSTOKES,MXHALF,MXK), & !single scattering contribution (upward)
          RSS_DN(NSTOKES,MXHALF,MXK)  !(downward)      

      logical ::  &
           SW_SURF                     !switch for surface reflectance

      real(double) ::  &
          RAD_UP(NSTOKES,MXHALF,0:MXK), & !Fourier comp. of upward intensity  
          RAD_DN(NSTOKES,MXHALF,0:MXK)  !Fourier comp. of downward intensity
      integer :: nf, mnk, ifwd_adj, i_src, i, k, i_st

!----------------------------------------------------------------

!     boundary condition for incident intensity at TOA

      do I=1,MXHALF           
         do I_ST = 1,NSTOKES
            RAD_DN(I_ST,I,0) = 0.0           
         enddo
      enddo

      do K=1,MNK               
         do I=1,MXHALF 
            do I_ST = 1,NSTOKES
               RAD_DN(I_ST,I,K) = RAD_DN(I_ST,I,K-1)*EXPMU(I,K) + &
                   RSS_DN(I_ST,I,K) 
            enddo
         enddo
      enddo

!     boundary condition using BDRF

      if(SW_SURF)then

         do I_ST = 1,NSTOKES
         do I=1,MXHALF
            RAD_UP(I_ST,I,MNK) = 1./PI * FLXDIR *& 
                BDRF_MS(I_ST,I_SRC,MXHALF+IFWD_ADJ,I,NF)
         enddo
         enddo

      else

         do I_ST = 1,NSTOKES
         do I = 1,MXHALF
            RAD_UP(I_ST,I,MNK) = 0.
         enddo
         enddo

      endif

!     add up single scattering contribution

      do K= MNK,1,-1
         do I_ST = 1,NSTOKES
            do I=1,MXHALF 
               RAD_UP(I_ST,I,K-1) = RAD_UP(I_ST,I,K)*EXPMU(I,K) + &
                   RSS_UP(I_ST,I,K)
            enddo
         enddo
      enddo

      return                                                                  
    end subroutine GAUSSEID_INIT_SSC_POL


end module gaus_mat_module
