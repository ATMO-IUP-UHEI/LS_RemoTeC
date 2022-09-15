  module ocean_fresnel_module
    use header_module
    use header_gsd
    implicit none

    contains
 !---------------------------------------------------------------------
 
      subroutine ocean_ss(&
           u0,&
           uv,&
           phi,&
           wspeed,&
           roc_ss)
        real(double), intent(in) :: u0, uv, phi, wspeed
        real(double), dimension(:), intent(out) :: roc_ss
!*** local
        real(double) :: sigma2, phi_in, phi_out, u0_in
!        complex*16 :: CN1,CN2
        real(double), dimension(4,4) :: R        
     
        ROC_SS = 0.d0
!        CN1=(1D0,0D0)
!        CN2=(1.33D0, 0.0D0)

!     Cox & Munk: SIGMA2 = 0.5 * (0.003 + 0.00512 *Windspeed[m/s])

        SIGMA2 = 0.5 * (0.003 + 0.00512 *wspeed)
        
        PHI_IN =  PHI /180.D0*PI
           
        PHI_OUT = 0.D0
           
        u0_in = u0  
        call fresnel_reflection_model(&
              u0_in,&
              uv,&
              phi_in,&
              1.33D0,&
              DSQRT(sigma2),&
              R)

        ROC_SS(1:nstokes) = R(1:nstokes,1)

        return
	
    end subroutine ocean_ss
    
!**********************************************************************

    subroutine ocean_fou(&
         u0,&
         uv,&
         wspeed,&
         gauss_quad,&         
         roc_ms)
      !*** input      
      real(double), intent(in) :: u0, uv, wspeed
      type(gauss_quad_type), intent(in) :: gauss_quad
      !*** output
      real(double), dimension(:,:,:,:,0:), intent(out) :: roc_ms
      !*** local
      real(double) :: sigma2, phi_in,phi_out,start_phi,end_phi,cosmphi,sinmphi,u_in,u_out
      !     complex*16 :: CN1,CN2
      real(double), dimension(4,4) :: R
      integer :: maxvza      
      integer, parameter :: nphi = 500     
      integer :: nf,i,j,np,i_st,j_st

      !      COMMON /DG_QUAD/&          !double gauss quadrature
      !           DG_MU(MAXSTR),&       !quadrature points
      !           DG_WT(MAXSTR)        !quadrature weights

      real(double), dimension(nphi) :: xphi, wphi

      real(double), dimension(:), allocatable :: u


      !---------------------------------------------------------------------

      maxvza = 1
      allocate(u(mxhalf+maxvza+1))

!!$      CN1=(1D0,0D0)
!!$      CN2=(1.33D0, 0.0D0)

      !     Cox & Munk: SIGMA2 = 0.5 * (0.003 + 0.00512 *Windspeed[m/s])

      SIGMA2 = 0.5 * (0.003 + 0.00512 * wspeed)

      PHI_OUT = 0.D0
      START_PHI = 0.D0
      END_PHI = 2.D0*PI 

      call gauleg_otto(NPHI,NPHI,START_PHI,END_PHI,XPHI,WPHI)

      do I=1,MXHALF
         U(I) = GAUSS_QUAD%DG_MU(I)
      enddo
      U(MXHALF+1) = U0
      U(MXHALF+2:MXHALF+2) = UV


      ROC_MS = 0.D0
      do I=1,MXHALF+MAXVZA+1             !incoming angles
         do J=1,MXHALF+MAXVZA+1          !outgoing angles

            U_IN = U(I)
            U_OUT = U(J)

            do NP = 1,NPHI
               PHI_IN = XPHI(NP)

               call fresnel_reflection_model(&
                    u_in,&
                    u_out,&
                    phi_in,&
                    1.33D0,&
                    DSQRT(sigma2),&
                    R)

               do NF=0,MAXSTR-1

                  COSMPHI = DCOS(dble(NF)*PHI_IN)
                  SINMPHI = DSIN(dble(NF)*PHI_IN)

                  do i_st=1, min(nstokes,3)
                     do j_st = 1, min(nstokes,3)
                        if (i_st<3) then
                           ROC_MS(i_st,j_st,I,J,NF) = ROC_MS(i_st,j_st,I,J,NF) +&
                                R(i_st,j_st)*WPHI(NP)/(2.D0*PI) * (COSMPHI-SINMPHI)
                        else
                           ROC_MS(i_st,j_st,I,J,NF) = ROC_MS(i_st,j_st,I,J,NF) +&
                                R(i_st,j_st)*WPHI(NP)/(2.D0*PI) * (COSMPHI+SINMPHI)
                        endif
                     enddo
                  enddo

               enddo
            enddo
         enddo
      enddo

      deallocate(U)

      return
    end subroutine ocean_fou

!**********************************************************************

    subroutine ocean_fou_lintranV2(&
         u0,&
         uv,&
         wspeed,&
         quad_mu,&         
         roc_ms)
      !*** input      
      real(double), intent(in) :: u0, uv, wspeed
      real(double), dimension(mxhalf), intent(in) :: quad_mu
      !*** output
      real(double), dimension(:,:,:,:,0:), intent(out) :: roc_ms
      !*** local
      real(double) :: sigma2, phi_in,phi_out,start_phi,end_phi,cosmphi,sinmphi,u_in,u_out
      !     complex*16 :: CN1,CN2
      real(double), dimension(4,4) :: R
      integer :: maxvza      
      integer, parameter :: nphi = 500     
      integer :: nf,i,j,np,i_st,j_st

      !      COMMON /DG_QUAD/&          !double gauss quadrature
      !           DG_MU(MAXSTR),&       !quadrature points
      !           DG_WT(MAXSTR)        !quadrature weights

      real(double), dimension(nphi) :: xphi, wphi

      real(double), dimension(:), allocatable :: u


      !---------------------------------------------------------------------

      maxvza = 1
      allocate(u(mxhalf+maxvza+1))

!!$      CN1=(1D0,0D0)
!!$      CN2=(1.33D0, 0.0D0)

      !     Cox & Munk: SIGMA2 = 0.5 * (0.003 + 0.00512 *Windspeed[m/s])

      SIGMA2 = 0.5 * (0.003 + 0.00512 * wspeed)

      PHI_OUT = 0.D0
      START_PHI = 0.D0
      END_PHI = 2.D0*PI 

      call gauleg_otto(NPHI,NPHI,START_PHI,END_PHI,XPHI,WPHI)

      do I=1,MXHALF
         U(I) = quad_mu(I)
      enddo
      U(MXHALF+1) = U0
      U(MXHALF+2:MXHALF+2) = UV


      ROC_MS = 0.D0
      do I=1,MXHALF+MAXVZA+1             !incoming angles
         do J=1,MXHALF+MAXVZA+1          !outgoing angles

            U_IN = U(I)
            U_OUT = U(J)

            do NP = 1,NPHI
               PHI_IN = XPHI(NP)

               call fresnel_reflection_model(&
                    u_in,&
                    u_out,&
                    phi_in,&
                    1.33D0,&
                    DSQRT(sigma2),&
                    R)

               do NF=0,MAXSTR-1

                  COSMPHI = DCOS(dble(NF)*PHI_IN)
                  SINMPHI = DSIN(dble(NF)*PHI_IN)

                  !                 ROC_MS(1,1,I,J,NF) = ROC_MS(1,1,I,J,NF) +&
                  !                      R(1,1)*WPHI(NP)/(2.D0*PI) * (COSMPHI-SINMPHI)

                  do i_st=1, min(nstokes,3)
                     do j_st = 1, min(nstokes,3)
                        if (i_st<3) then
                           ROC_MS(i_st,j_st,I,J,NF) = ROC_MS(i_st,j_st,I,J,NF) +&
                                R(i_st,j_st)*WPHI(NP)/(2.D0*PI) * (COSMPHI-SINMPHI)
                        else
                           ROC_MS(i_st,j_st,I,J,NF) = ROC_MS(i_st,j_st,I,J,NF) +&
                                R(i_st,j_st)*WPHI(NP)/(2.D0*PI) * (COSMPHI+SINMPHI)
                        endif
                     enddo
                  enddo



               enddo
            enddo
         enddo
      enddo

      deallocate(U)

      return
    end subroutine ocean_fou_lintranV2

!****************************************************************************

  subroutine fresnel_reflection_model(&
       mu_in,&
       mu_out,&
       phi_in,&
       Rm_surf,&
       sigma,&
       R)
!-------------------------------------------------------------------------!
! NOTE: Rm_surf is the refractive index of the medium to which radiation  !
! is transmitted divided by the refractive index of the medium from       !
! which the radiation originates.                                         !
!--------------------------------------------------------------------------
!*** input
    real(double), intent(inout) :: mu_in
    real(double), intent(in) :: mu_out, phi_in, Rm_surf, sigma
!*** output
    real(double), dimension(4,4), intent(out) :: R
!*** local
    real(double), dimension(4,4) ::   Rmat,  R_Fresnel, F_rot1, F_rot2
    real(double) :: u_in,u_out, phi
    real(double), parameter :: epsilon=1.D-6
    real(double) :: &
         cos_scat, scat_angle, theta_in, theta_trans,&
         alpha, eta, gamma_re, cos_i1, cos_i2, ang_i1,ang_i2,&
         mu_n, pdf, v_in, v_out, F_v_in, F_v_out, shadow,&
         fac1, sign
!----------------------------------------------------------



    sign = 1.D0
    if(DSIN(phi_in).lt.0)sign=-1.D0

    PHI = sign*(DACOS(DCOS(phi_in)))


    if(mu_in.eq.mu_out)mu_in = mu_in-epsilon

    u_in = -mu_in
    u_out = mu_out


    cos_scat =  &
         u_in*u_out + &
         DSQRT((1.D0-u_out**2)) * DSQRT((1.D0-u_in**2)) * &
         DCOS(phi)
    
    scat_angle = DACOS(cos_scat)

    theta_in = 0.5D0*(PI-scat_angle)
      
    theta_trans = DASIN(DSIN(theta_in)/Rm_surf)
      


!----------------------------------------------------
! facet orientation

    mu_n = (mu_in+mu_out)/(2.D0*DCOS(theta_in))

    fac1 = 1./(4.*mu_in*mu_out*mu_n)

    pdf = (1.D0/(2.*sigma**2*mu_n**3))*DEXP(-((1-mu_n**2)/(2.*sigma**2*mu_n**2)))

!    PRINT*,'pdf',pdf,fac1

!-----------------------------------------------------
! shadowing

    v_in = mu_in / (sigma*(DSQRT(1-mu_in**2)))
    v_out = mu_out / (sigma*(DSQRT(1-mu_out**2)))
    
    F_v_in = &
         DEXP(-v_in**2)/ DSQRT(PI*v_in**2) - DERFC(v_in)

    F_v_in = 0.5D0*F_v_in

    F_v_out = &
         DEXP(-v_out**2)/ DSQRT(PI*v_out**2) - DERFC(v_out)

    F_v_out = 0.5D0*F_v_out

    shadow = 1.D0 / (1.D0+F_v_in+F_v_out)

!-------------------------------------------------------------------
   

    alpha = & 
         0.5D0*( DTAN(theta_in-theta_trans) /  &
         DTAN(theta_in+theta_trans) )**2

    
    eta = &
         0.5D0*( DSIN(theta_in-theta_trans) / &
         DSIN(theta_in+theta_trans) )**2



    gamma_re = &
         -( DTAN(theta_in-theta_trans) * &
         DSIN(theta_in-theta_trans) ) / &
         ( DTAN(theta_in+theta_trans) * &
         DSIN(theta_in+theta_trans) )

    R_Fresnel = 0.D0          !initialize
    
    R_Fresnel(1,1) = alpha + eta
    R_Fresnel(2,2) = alpha + eta
    R_Fresnel(3,3) = gamma_re
    R_Fresnel(4,4) = gamma_re
    R_Fresnel(2,1) = alpha - eta
    R_Fresnel(1,2) = alpha - eta
    

! Transformation matrix
    
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
    




    Rmat = matmul(F_rot1,(matmul(R_Fresnel,F_rot2)))



    R = shadow*fac1*pdf*Rmat


  end subroutine fresnel_reflection_model



!---------------------------------------------------------------


!!$  SUBROUTINE fresnel_ocean_model_wind_direction(&
!!$       mu_in,&
!!$       mu_out,&
!!$       phi_in,&
!!$       Rm_surf,&
!!$       wspeed,&
!!$       wind_direction_in,&
!!$       R)
!!$    IMPLICIT NONE
!!$!-------------------------------------------------------------------------!
!!$! NOTE: Rm_surf is the refractive index of the medium to which radiation  !
!!$! is transmitted divided by the refractive index of the medium from       !
!!$! which the radiation originates.                                         !
!!$!--------------------------------------------------------------------------
!!$
!!$    REAL(DOUBLE) :: PI
!!$
!!$    REAL(KIND=8), DIMENSION(4,4) :: & 
!!$         Rmat, R_Fresnel, F_rot1, F_rot2, R
!!$
!!$    REAL(KIND=8) :: wspeed, wind_direction_in
!!$
!!$    REAL(KIND=8) :: mu_in,mu_out,u_in,u_out, phi,phi_in,Rm_surf, sigma, wind_direction
!!$
!!$    REAL(KIND=8), PARAMETER :: epsilon=1.D-6
!!$
!!$    REAL(KIND=8) :: Z_x, Z_y, Zx_prime, Zy_prime, sigma_x, sigma_y, Zeta_squared, eta_squared,theta_s,theta_v
!!$
!!$    REAL(KIND=8) :: &
!!$         cos_scat, scat_angle, theta_in, theta_trans,&
!!$         alpha, eta, gamma_re, cos_i1, cos_i2, ang_i1,ang_i2,&
!!$         mu_n, pdf, v_in, v_out, F_v_in, F_v_out, shadow,&
!!$         factor,fac1,fac2,cos2chi, sign,Beta, alpha_prime,omega,cos_alpha
!!$!----------------------------------------------------------
!!$    PI = 2*DASIN(1.D0)
!!$
!!$
!!$    wind_direction = DACOS(wind_direction_in)
!!$
!!$    sign = 1.D0
!!$    IF(DSIN(phi_in).LT.0)sign=-1.D0
!!$    PHI = sign*(DACOS(DCOS(phi_in)))
!!$
!!$
!!$
!!$
!!$    IF(mu_in.EQ.mu_out)mu_in = mu_in-epsilon
!!$
!!$    u_in = -mu_in
!!$    u_out = mu_out
!!$
!!$
!!$    cos_scat =  &
!!$         u_in*u_out + &
!!$         DSQRT((1.D0-u_out**2)) * DSQRT((1.D0-u_in**2)) * &
!!$         DCOS(phi)
!!$    
!!$    scat_angle = DACOS(cos_scat)
!!$
!!$    theta_in = 0.5D0*(PI-scat_angle)
!!$
!!$    omega =  0.5D0*(PI-scat_angle)
!!$ 
!!$
!!$     
!!$    theta_trans = DASIN(DSIN(theta_in)/Rm_surf)
!!$
!!$!---------------------------------------------------------------------
!!$
!!$    theta_v = DACOS(mu_out)
!!$    theta_s = DACOS(mu_in)
!!$
!!$    Z_x = -DSIN(theta_v)*DSIN(pi+phi) / (DCOS(theta_s)+DCOS(theta_v))
!!$
!!$    Z_y = DSIN(theta_s) + DSIN(theta_v)*DCOS(pi+phi) / (DCOS(theta_s)+DCOS(theta_v))
!!$
!!$
!!$    Zx_prime = DCOS(pi+wind_direction)*Z_x + DSIN(pi+wind_direction)*Z_y
!!$    ZY_prime = -DSIN(pi+wind_direction)*Z_x + DCOS(pi+wind_direction)*Z_y
!!$!----------------------------------------------------------------------
!!$
!!$!    Beta = DACOS((mu_in+mu_out)/(2.*DCOS(theta_in)))
!!$!    alpha = DACOS((DSIN(DACOS(mu_out))*DCOS(phi)-DSIN(DACOS(mu_in)))/(2.*DCOS(theta_in)))
!!$
!!$    Beta = DACOS((mu_in+mu_out)/(2.*DCOS(omega)))
!!$
!!$    cos_alpha = ((DSIN(theta_v)*DCOS(phi)-DSIN(theta_s))/(2.*DCOS(omega)*DSIN(Beta)))
!!$    IF(cos_alpha.GT.1.D0)cos_alpha = 1.D0
!!$    IF(cos_alpha.LT.1.D0)cos_alpha = -1.D0
!!$
!!$    alpha = DABS(DACOS(cos_alpha))
!!$
!!$    IF(phi.LT.0.)alpha = -alpha
!!$ 
!!$
!!$
!!$    alpha_prime = ((alpha-wind_direction))
!!$
!!$!    SIGMA = DSQRT((0.003 + 0.00512 *wspeed))    
!!$!    sigma_x = 1./DSQRT(2.D0) *sigma
!!$!    sigma_y = sigma_x
!!$
!!$
!!$    sigma_x = DSQRT((0.003+0.00192*wspeed))
!!$    sigma_y = DSQRT((0.00316*wspeed))
!!$    sigma = DSQRT(sigma_x**2+sigma_y**2)
!!$
!!$    Zeta_squared = (DSIN(alpha_prime)*DTAN(Beta)/sigma_x)**2
!!$    eta_squared = (DCOS(alpha_prime)*DTAN(Beta)/sigma_y)**2
!!$
!!$
!!$!    Zeta_squared = (Zx_prime / sigma_x)**2
!!$!    eta_squared = (Zy_prime / sigma_y)**2
!!$
!!$!----------------------------------------------------
!!$! facet orientation
!!$
!!$    mu_n = (mu_in+mu_out)/(2.D0*DCOS(theta_in))
!!$
!!$!    PRINT*,alpha,cos_alpha,wind_direction
!!$
!!$
!!$    fac1 = 1./(4.*mu_in*mu_out*mu_n)
!!$
!!$!    pdf = (1.D0/(2.*sigma**2*mu_n**3))*DEXP(-((1-mu_n**2)/(2.*sigma**2*mu_n**2)))
!!$!    pdf = (1.D0/(2.*sigma_x*sigma_y*mu_n**3))*DEXP(-((1-mu_n**2)/(2.*sigma**2*mu_n**2)))
!!$    
!!$    pdf = (1.D0/(2.*sigma_x*sigma_y*mu_n**3))*DEXP(-(Zeta_squared+eta_squared)/2.D0)
!!$
!!$
!!$!    PRINT*,'pdf',pdf,eta_squared,zeta_squared
!!$!-----------------------------------------------------
!!$! shadowing
!!$
!!$    v_in = mu_in / (sigma*(DSQRT(1-mu_in**2)))
!!$    v_out = mu_out / (sigma*(DSQRT(1-mu_out**2)))
!!$    
!!$    F_v_in = &
!!$         DEXP(-v_in**2)/ DSQRT(PI*v_in**2) - DERFC(v_in)
!!$
!!$    F_v_in = 0.5D0*F_v_in
!!$
!!$    F_v_out = &
!!$         DEXP(-v_out**2)/ DSQRT(PI*v_out**2) - DERFC(v_out)
!!$
!!$    F_v_out = 0.5D0*F_v_out
!!$
!!$    shadow = 1.D0 / (1.D0+F_v_in+F_v_out)
!!$
!!$!-------------------------------------------------------------------
!!$   
!!$
!!$    alpha = & 
!!$         0.5D0*( DTAN(theta_in-theta_trans) /  &
!!$         DTAN(theta_in+theta_trans) )**2
!!$
!!$    
!!$    eta = &
!!$         0.5D0*( DSIN(theta_in-theta_trans) / &
!!$         DSIN(theta_in+theta_trans) )**2
!!$
!!$
!!$
!!$    gamma_re = &
!!$         -( DTAN(theta_in-theta_trans) * &
!!$         DSIN(theta_in-theta_trans) ) / &
!!$         ( DTAN(theta_in+theta_trans) * &
!!$         DSIN(theta_in+theta_trans) )
!!$
!!$    R_Fresnel = 0.D0          !initialize
!!$    
!!$    R_Fresnel(1,1) = alpha + eta
!!$    R_Fresnel(2,2) = alpha + eta
!!$    R_Fresnel(3,3) = gamma_re
!!$    R_Fresnel(4,4) = gamma_re
!!$    R_Fresnel(2,1) = alpha - eta
!!$    R_Fresnel(1,2) = alpha - eta
!!$    
!!$
!!$! Transformation matrix
!!$    
!!$    cos_i1 = 1.D0
!!$    cos_i2 = 1.D0
!!$
!!$    IF(cos_scat**2.NE.1)THEN
!!$       
!!$       cos_i1 = ( &
!!$            ( u_out*DSQRT(1.D0-u_in**2) - &
!!$            u_in *DSQRT(1.D0-u_out**2)*DCOS(PHI) ) / &
!!$            (DSQRT(1.D0-cos_scat**2)) )
!!$       
!!$       cos_i2 = ( &
!!$            ( u_in*DSQRT(1.D0-u_out**2) - &
!!$            u_out *DSQRT(1.D0-u_in**2)*DCOS(PHI) ) / &
!!$            (DSQRT(1.D0-cos_scat**2)) )
!!$       
!!$       cos_i1 = MIN(cos_i1,1.D0)
!!$       cos_i2 = MIN(cos_i2,1.D0)
!!$       
!!$       cos_i1 = MAX(cos_i1,-1.D0)
!!$       cos_i2 = MAX(cos_i2,-1.D0)
!!$         
!!$         
!!$
!!$    ENDIF
!!$
!!$    ang_i1 = DACOS(cos_i1)
!!$    ang_i2 = DACOS(cos_i2)
!!$
!!$
!!$    IF(phi.LT.0.)THEN
!!$       ang_i1 = -ang_i1
!!$       ang_i2 = -ang_i2
!!$    ENDIF
!!$
!!$    
!!$
!!$    F_rot1 = 0.D0              !initialize
!!$
!!$    F_rot1(1,1) = 1.D0
!!$    F_rot1(2,2) = DCOS(2*ang_i2)
!!$    F_rot1(3,3) = DCOS(2*ang_i2)
!!$    F_rot1(2,3) = -DSIN(2*ang_i2)
!!$    F_rot1(3,2) = DSIN(2*ang_i2)
!!$    
!!$    F_rot2 = 0.D0              !initialize
!!$
!!$    F_rot2(1,1) = 1.D0
!!$    F_rot2(2,2) = DCOS(2*ang_i1)
!!$    F_rot2(3,3) = DCOS(2*ang_i1)
!!$    F_rot2(2,3) = -DSIN(2*ang_i1)
!!$    F_rot2(3,2) = DSIN(2*ang_i1)
!!$    
!!$
!!$
!!$
!!$
!!$    Rmat = MATMUL(F_rot1,(MATMUL(R_Fresnel,F_rot2)))
!!$
!!$
!!$
!!$    R = shadow*fac1*pdf*Rmat
!!$
!!$
!!$  END SUBROUTINE fresnel_ocean_model_wind_direction

!---------------------------------------------------------------------

      subroutine gauleg_otto(ndim,ngauss,a,b,x,w)

!**********************************************************************
!* Given the lower and upper limits of integration a and b, 
!* and given the number of Gauss-Legendre points ngauss,
!* this routine returns through array x the abscissas and through
!* array w the weights of the Gauss-Legendre quadrature formula.
!* Eps is the desired accuracy of the abscissas.
!* This routine is documented further in:
!*   W.H. Press et al. 'Numerical Recipes' Cambridge Univ. Pr. (1987)
!*   page 125 ISBN 0-521-30811-9                                   
!**********************************************************************
!      IMPLICIT DOUBLE PRECISION (a-h,o-z)
!*** input
   integer, intent(in) :: ndim, ngauss
   real(double), intent(in) :: a, b
!*** output
   real(double), intent(out) ::  x(ndim),w(ndim)
!*** local
   integer :: i, j, m 
   real(double), parameter :: eps= 1.d-14
   double precision ::  xm,xl,z,p1,p2,p3,pp,z1!,pi

!**********************************************************************
!      pi=4.D0*datan(1.D0)
      m=(ngauss+1)/2
      xm=0.5D0*(a+b)
      xl=0.5D0*(b-a)

      do 12 i=1,m
         z= dcos(pi*(dble(i)-0.25D0)/(dble(ngauss)+0.5D0))
1        continue
            p1=1.D0
            p2=0.D0
            do j=1,ngauss
               p3= p2
               p2= p1
               p1=((dble(2*j)-1.d0)*z*p2-(dble(j)-1.d0)*p3)/dble(j)
            enddo
            pp=ngauss*(z*p1-p2)/(z*z-1.d0)
            z1= z
            z= z1-p1/pp
          if (dabs(z-z1).gt.eps) goto 1
          x(i)= xm-xl*z
          x(ngauss+1-i)= xm+xl*z
          w(i)=2.D0*xl/((1.D0-z*z)*pp*pp)
          w(ngauss+1-i)= w(i)
12    continue
!**********************************************************************
      return
      end subroutine gauleg_otto

end module ocean_fresnel_module
