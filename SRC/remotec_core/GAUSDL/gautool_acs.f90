module gautool_acs_module

contains
  !************************************************************************
  !*** calculates the generalized spherical functions needed for the   
  !*** expansian of the phase matrix. Follows 'De Haan et al. 1987'    
  subroutine SET_PLM(&
       UVIEW,&
       U0,&
       gauss_quad,&
       gsf)
    use header_gsd
    implicit none
    !*** Input
    type(gauss_quad_type), intent(IN) :: gauss_quad
    real(double), intent(in) :: uview, u0
    !*** Output
    type(gsf_type), intent(OUT) :: gsf
    !*** Local variables
    integer :: i
    integer, parameter :: MAXM = MAXSTR-1, MAXI = 2
    real(double), dimension(MAXSTR+4) :: u    
    !---------------------------------------------------------------------------
    !     U is -mu from the paper so that P^l_m(-mu) = P^l_m(U)
    !     so U(1:MXHALF) describes downward streams,  U(1+MXHALF:2*MXHALF) the
    !     corresponding upward streams,  U(MAXSTR + 1) the solar beam and
    !     U(MAXSTR + 2) the viewing geometry. U(MAXSTR+3) and U(MAXSTR+4) are
    !     needed for the adjoint formulation.

    do I = 1, MXHALF
       U(I)        = -gauss_quad%DG_MU(I)
       U(MXHALF+I) =  gauss_quad%DG_MU(I)
    enddo

    U(MAXSTR + 1) = -U0                       
    U(MAXSTR + 2) =  UVIEW
    U(MAXSTR + 3) = -UVIEW
    U(MAXSTR + 4) =  U0

    call GENERALIZED_SPHERICAL_FCT(&
         MAXM,  MAXSTR+4, U, &
         gsf%PLM_0, gsf%PLM_M, gsf%PLM_P)

    return

  end subroutine SET_PLM

  !************************************************************************
  !*** calculates the generalized spherical functions needed for the
  !*** expansian of the phase matrix. Follows 'De Haan et al. 1987'  
  subroutine GENERALIZED_SPHERICAL_FCT(&
       MAXM, MAXI, U, &
       PLM0, PLMM, PLMP)
    use header_gsd
    implicit none
    integer, intent(in) :: maxm, maxi  
    real(double), dimension(maxi), intent(in) :: u
    real(double), dimension(0:MAXM,0:MAXM,MAXI), intent(out) :: PLM0, PLMM, PLMP
    !*** Local variables
    integer :: i, l, m, n
    real(double) :: bin2mm, twom, binfac, rootu, binf2, fac1p, bin2m2, urootm, fac1m, fac2, fac3
    real(double), dimension(0:MAXM,0:MAXM) :: plmp2, plmm2
    !-----------------------------------------------------------------------
    !*** start loop over angles
    do I = 1,MAXI

       !*** initialization
       do M=0,MAXM
          do L=0,MAXM
             PLM0(M,L,I) = 0.
             PLMP2(M,L)  = 0.
             PLMM2(M,L)  = 0.
          enddo
       enddo

       do M = 0, MAXM
          !*** Compute the binomial factor 2m above m needed in Eq.(77)    
          !*** of de Haan et al. 87  
          BIN2MM = 1.D0
          do N=1,M
             BIN2MM = BIN2MM*dble(N+M)/dble(N)
          enddo
          TWOM   = 2.D0**(-M)
          BINFAC = TWOM*DSQRT(BIN2MM) 
          !            BINFAC = DSQRT(FAC(2*M)/(FAC(M)*FAC(M)))
          !            BINFAC = 2.D0**(-M)*BINFAC
          ROOTU  = DSQRT(DABS(1.D0-U(I)**2))

          !*** start with the generalized spherical function P^m_{m,0}
          !*** Eq. (77) and ignore factor i
          if (M.eq.0)then
             PLM0(M,M,I) = BINFAC 
          else
             PLM0(M,M,I) = BINFAC*ROOTU**M
          endif

          !*** use recurrence relatation Eq. 81 to calculate higher
          !*** P^l_{m,0}, with l = M+1,MAXLEG 
          do L = M,MAXM-1
             if(L.ne.0)then
                PLM0(M,L+1,I) = dble(2.*L+1)*U(I)*PLM0(M,L,I) -&
                     DSQRT(dble(L**2-M**2))*PLM0(M,L-1,I)
             else
                PLM0(M,L+1,I) = dble(2.*L+1)*U(I)*PLM0(M,L,I) 
             endif
             PLM0(M,L+1,I) = PLM0(M,L+1,I)/DSQRT(dble((L+1)**2-M**2))
          enddo
       enddo

       !*** generalized spherical function P^2_{m,+2} and  P^2_{m,-2}
       !*** for m = 0,1,2
       PLMP2(0,2) = -1./4.*DSQRT(6.D0) * (1.-U(I)**2)
       PLMM2(0,2) =  PLMP2(0,2)

       PLMP2(1,2) =  1./2.*DSQRT(1.-U(I)**2) * (1.+U(I))
       PLMM2(1,2) = -1./2.*DSQRT(1.-U(I)**2) * (1.-U(I))

       PLMP2(2,2) = -1./4.*(1.D0 + U(I))**2
       PLMM2(2,2) = -1./4.*(1.D0 - U(I))**2

       do M = 3, MAXM
          !*** Compute the binomial factor 2m above m-2 needed in 
          !*** Eq.(80) of de Haan et al. 87 
          BIN2MM = 1.D0
          do N=1,M
             BIN2MM = BIN2MM*(dble(N+M)/dble(N))
          enddo
          TWOM   = 2.D0**(-M)
          BINFAC = TWOM*DSQRT(BIN2MM)
          BIN2M2 = BIN2MM*dble(M)*dble(M-1)/(dble(M+1)*&
               dble(M+2))
          BINF2  = -TWOM*DSQRT(BIN2M2)            
          ROOTU = DSQRT(DABS(1.D0-U(I)**2))
          UROOTM = ROOTU**(M-2)
          BINF2  = -TWOM*DSQRT(BIN2M2) 
          !            TWOM   = 2.D0**(-M)
          !            BINF2  = -TWOM*DSQRT(FAC(2*M)/(FAC(M+2)*FAC(M-2)))
          !            ROOTU  = DSQRT(DABS(1.D0-U(I)**2))
          !            UROOTM = ROOTU**(M-2)
          PLMM2(M,M) = BINF2*UROOTM*(1.-U(I))*(1.-U(I))
          PLMP2(M,M) = BINF2*UROOTM*(1.+U(I))*(1.+U(I))           
       enddo

       do M = 0, MAXM
          do L = 0, MAXM-1               
             if(L.ge.max(M,2))then              
                FAC1M = (2.D0*L+1D0)*(L*(L+1D0)*U(I)-2.*M)
                FAC1P = (2.D0*L+1D0)*(L*(L+1D0)*U(I)+2.*M)
                FAC2  = (L+1)*DSQRT((L**2D0-4D0)*(L**2D0-M**2))
                FAC3  = 1D0/(L*DSQRT(((L+1D0)**2-4D0)*((L+1D0)**2-M**2)))
                PLMP2(M,L+1) = FAC3*(FAC1M*PLMP2(M,L)-FAC2*PLMP2(M,L-1))
                PLMM2(M,L+1) = FAC3*(FAC1P*PLMM2(M,L)-FAC2*PLMM2(M,L-1))               
             endif
          enddo
       enddo

       !*** finally, the general spherical function  P^2_{m,+} and  P^2_{m,-}
       do M = 0,MAXM
          do L = 0,MAXM
             PLMP(M,L,I) = 1./2.*(PLMM2(M,L)+PLMP2(M,L))
             PLMM(M,L,I) = 1./2.*(PLMM2(M,L)-PLMP2(M,L))
          enddo
       enddo

    enddo ! end loop over angles

    return

  end subroutine GENERALIZED_SPHERICAL_FCT

  !****************************************************************************** 
  !*     Set up double gaussian quadrature:
  !*     The sign of DG_MU is chosen that DG_MU and -DG_MU denote the downward
  !*     and upward direction, respectively. Thus, U0 is positive, and the 
  !*     direction of the direct beam is (U0,0) refering to the principle
  !*     plane.   
  subroutine DOUBLE_GAUSS(gauss_quad, ierr)
    use header_gsd         
    implicit none                 
    type(gauss_quad_type), intent(OUT) :: gauss_quad
    integer, intent(out) :: ierr
    !*** Local variables
    integer :: i
    real(double), dimension(maxstr) :: tmp_mu,tmp_wt
    !----------------------------------------------------------------------------
    ierr = 0 
    if((MAXSTR+1)/2 .ne. MAXSTR/2 .or. MAXSTR/2 .lt. 1 ) then
       ierr = ierr_var
       call stopretrieval('DOUBLE_GAUSS: Bad value of MAXSTR! Choose an even integer for MAXSTR.')
       return
    endif

    call DGAUSS(MAXSTR/2, TMP_MU, TMP_WT) !only positive values are calculated

    do I = 1, MAXSTR/2                     !downward direction
       gauss_quad%DG_MU(MAXSTR/2-I+1) =  TMP_MU(I) !resort MAXSTR/2   -> 1,
       !                                                 MAXSTR/2-1 -> 2, ... 
       gauss_quad%DG_WT(MAXSTR/2-I+1) =  TMP_WT(I)
    enddo
    do I = 1, MAXSTR/2                     !upward values
       gauss_quad%DG_MU(MAXSTR+1-I)   = -gauss_quad%DG_MU(I)
       gauss_quad%DG_WT(MAXSTR+1-I)   =  gauss_quad%DG_WT(I)
    enddo

    return   

  end subroutine DOUBLE_GAUSS

  !***************************************************************************
  !*     Compute weights and abscissae for ordinary gaussian quadrature      *
  !*     (no weight function inside integral) on the interval (0,1)          *
  !*                                                                         *
  !*     input :    M                      order of quadrature rule          *
  !*                                                                         *
  !*     output :   GMU(I)  I = 1 TO M     array of abscissae                *
  !*                GWT(I)  I = 1 TO M     array of weights                  *
  !*                                                                         *
  !*     reference:  Davis, P.J. and P. Rabinowitz,                          *
  !*                 Methods of Numerical Integration,                       *
  !*                 Academic Press, New York, pp. 87, 1975.                 *
  !*                                                                         *
  !*     method:  Compute the abscissae as roots of the Legendre             *
  !*              polynomial P-sub-M using a cubically convergent            *
  !*              refinement of Newton's method.  Compute the                *
  !*              weights from EQ. 2.7.3.8 of Davis/Rabinowitz.  Note        *
  !*              that Newton's method can very easily diverge; only a       *
  !*              very good initial guess can guarantee convergence.         *
  !*              The initial guess used here has never led to divergence    *
  !*              even for M up to 1000.                                     *
  !*                                                                         *
  !*     accuracy:  at least 13 significant digits                           *
  !*                                                                         *
  !*     internal variables:                                                 *
  !*                                                                         *
  !*    ITER      : number of Newton Method iterations                       *
  !*    MAXIT     : maximum allowed iterations of Newton Method              *
  !*    PM2,PM1,P : 3 successive Legendre polynomials                        *
  !*    PPR       : derivative of Legendre polynomial                        *
  !*    P2PRI     : 2nd derivative of Legendre polynomial                    *
  !*    TOL       : convergence criterion for Legendre poly root iteration   *
  !*    X,XI      : successive iterates in cubically-convergent version      *
  !*                of Newtons Method (seeking roots of Legendre poly.)      *
  subroutine  DGAUSS( M, GMU, GWT)
    use header_gsd
    implicit none
    !*** input
    integer, intent(in) :: m
    !*** output
    real(double), intent(out) :: GMU(*), GWT(*)
    !*** Local variables 
    integer :: k, nn
    real(double) :: CONA, T
    integer, parameter :: MAXIT = 1000
    integer :: ITER, LIM, NP1	
    real(double) :: EN, NNP1, P, PM1, PM2, PPR, P2PRI, PROD, TMP, X, XI   
    real(double), parameter :: ONE = 1.D0, TWO = 2.D0, EPSILON = 1.D-14
    !-----------------------------------------------------------------------

!!$    if ( M.lt.1 )  call stopretrieval('DGAUSS: Bad value of M')

    if ( M.eq.1 )  then
       GMU( 1 ) = 0.5
       GWT( 1 ) = 1.0
       return
    end if

    EN   = M
    NP1  = M + 1
    NNP1 = M * NP1
    CONA = FLOAT( M-1 ) / ( 8 * M**3 )
    LIM  = M / 2
    do  K = 1, LIM
       !*** initial guess for k-th root of Legendre polynomial, 
       !*** from Davis/Rabinowitz (2.7.3.3a)
       T = ( 4*K - 1 ) * PI / ( 4*M + 2 )
       X = cos ( T + CONA / DTAN( T ) )
       ITER = 0

       !*** upward recurrence for Legendre polynomials
10     ITER = ITER + 1
       PM2 = ONE
       PM1 = X
       do NN = 2, M
          P   = ( ( 2*NN - 1 ) * X * PM1 - ( NN-1 ) * PM2 ) / NN
          PM2 = PM1
          PM1 = P
       enddo

       !*** Newton method
       TMP   = ONE / ( ONE - X**2 )
       PPR   = EN * ( PM2 - X * P ) * TMP
       P2PRI = ( TWO * X * PPR - NNP1 * P ) * TMP
       XI    = X - ( P / PPR ) * ( ONE +( P / PPR ) * P2PRI / ( TWO * PPR ) )

       !*** check for convergence
       if ( DABS(XI-X) .gt. EPSILON) then
          if( ITER.gt.MAXIT ) then
             call stopretrieval('DGAUSS: MAX ITERATION COUNT')
          endif
          X = XI
          GO TO 10
       endif

       !*** iteration finished--calculate weights, abscissae for (-1,1)
       GMU( K ) = - X
       GWT( K ) = TWO / ( TMP * ( EN * PM2 )**2 )
       GMU( NP1 - K ) = - GMU( K )
       GWT( NP1 - K ) =   GWT( K )         
    enddo

    !*** set middle abscissa and weight for rules of odd order
    if ( mod( M,2 ) .ne. 0 )  then
       GMU( LIM + 1 ) = 0.0
       PROD = ONE
       do  K = 3, M, 2
          PROD = PROD * K / ( K-1 )
       enddo
       GWT( LIM + 1 ) = TWO / PROD**2
    endif

    !*** convert from (-1,1) to (0,1) and resort
    do K = 1, M
       GMU(K) = 0.5 * GMU( K ) + 0.5
       GWT(K) = 0.5 * GWT( K )
    enddo

  end subroutine  DGAUSS

!*****************************************************************************
!*** attenuation factors at various atmospheric layers
  subroutine ATTEN_FAC( &
       nrt, &
       TAUA,&
       TAUS,&
       U0,&
       UV,&
       gauss_quad,&
       OMEGA,&
       TAUTOT,& 
       EXPMU,&
       EXPU0,&
       EXPUV,&
       EINTU0,&
       EINTUV,&
       MNK,&
       NFILT1,&
       NFILT2,&
       SW_SURF, &
       ierr)
    use header_gsd         
    !*** Input
    integer, intent(in) :: nrt
    real(double), intent(in) :: &
         TAUA(NRT), &         !absorption optical depth
         TAUS(NRT), &         !scattering optical depth
         U0, &                   !cosine od SZA
         UV                      !cosine of viewing zenith angle
    type(gauss_quad_type), intent(IN) :: gauss_quad
    !*** Output
    real(double), intent(out) :: &
         TAUTOT(MXK),    &        !optical depths of the internal layers
         OMEGA(NRT),  &        !single scatering albedo
         EXPMU(MXHALF,MXK),  &    !attenuation factors EXP(-tau(k)/mu(i))
         EXPUV(MXK), &            !attenuation factors EXP(-tau(k)/U0)
         EXPU0(MXK),  &           !attenuation factors EXP(-tau(k)/UV)
         EINTU0(0:MXK), &         !exp(-TAUINT(K)/U0)
         EINTUV(0:MXK)            !exp(-TAUINT(K)/UV) where TAUINT is the integrated optical depth
    integer, intent(out) ::  &
         MNK,   &               !number of internal layers
         NFILT2(NRT),  &        !filter for auto. cell splitting
         NFILT1(MXK), &
         ierr                   ! error identifier 
    !*** Local variables
    integer :: i, ik, k, mslay, malay, mlay
    logical :: SW_SURF            !switch for ground reflection 
    real(double) :: tauint
    real(double), parameter :: EPS = 1.D-10
    real(double), parameter :: tslim = 0.02
    real(double), parameter :: talim = 0.4
    !----------------------------------------------------------------------------

    ierr = 0

    SW_SURF = .false.
    IK = 0
    TAUINT = 0      
    do K = 1, NRT
       NFILT2(K) = 0
    enddo
    do K = 1, NRT !***ABUTZ: K=1,LAYBOT IF NEC.
       MSLAY = aint(TAUS(K)/TSLIM) + 1
       MALAY = aint(TAUA(K)/TALIM) + 1
       MLAY = max(MSLAY,MALAY)
       NFILT2(K) = MLAY
       do I = 1,MLAY
          if (IK .ge. MXK) then
             ierr = 1
             call writelog("ATTEN_FAC: Number of sublayers exceeds array dimension", 6)
             return
          endif
          IK = IK + 1
          TAUTOT(IK) = (TAUS(K) + TAUA(K))/dble(MLAY)
          TAUINT = TAUINT + (TAUA(K))/dble(MLAY)!HH: TAUTOT(IK) !Check if abs OT is too large (not scat OT)
          NFILT1(IK) = K
       enddo
       if(TAUINT.gt.TATOT) goto 100
    enddo
    SW_SURF = .true.
100 continue
    MNK = IK
    !*** HH, do not go over array boundaries
    if (MNK > MXK) then
       ierr = 1
       call writelog("ATTEN_FAC: Number of sublayers exceeds array dimension", 6)
       return
    endif

    do I = 1, MXHALF                         
       do K = 1, MNK                     
          EXPMU(I,K) = exp(-TAUTOT(K)/gauss_quad%DG_MU(I))
       enddo
    enddo
    EINTU0(0) = 1.D0
    EINTUV(0) = 1.D0
    do K=1,MNK   
       EXPU0(K) = exp(-TAUTOT(K)/U0)
       EXPUV(K) = exp(-TAUTOT(K)/abs(UV))
       EINTU0(K)= EINTU0(K-1)*EXPU0(K)   
       EINTUV(K)= EINTUV(K-1)*EXPUV(K)   
    enddo
    do K = 1,NRT
       OMEGA(K) = min(1.D0-EPS,TAUS(K)/(TAUS(K) + TAUA(K)))
    enddo

    return            

  end subroutine ATTEN_FAC
!*********************************************************************
end module gautool_acs_module
