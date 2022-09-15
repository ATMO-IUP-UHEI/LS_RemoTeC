module gaus_pol_module
   use header_gsd
   use gaus_mat_module
   implicit none
   
   contains
!------------------------------------------------------------------------------ 
!> using Gauss Seidel iteration scheme for solving M I = C           
!------------------------------------------------------------------------------  
   subroutine GAUSS_SEIDEL(&
              nrt, &
              NF,&
              MNK,& 
              NFILT1, &
              NMAT,& 
              IFWD_ADJ,&
              ISRC,&   
              EPST,&   
              UBEAM,& 
              SW_SURF, &
              pm,&
              gauss_quad,&
              BDRF_MS,& 
              OMEGA,&   
              EXPMU,&
              EXPBEAM,& 
              EXPINT,& 
              RAD_UP,&  
              RAD_DN,&
              ncoefs)
!*** input
   integer, intent(in) ::  nrt, NF                        !index of fourier component
   integer, dimension(nrt), intent(in) :: ncoefs
   type(phase_mat_type), intent(in) :: pm
   type(gauss_quad_type), intent(in) :: gauss_quad
   logical, intent(in) :: SW_SURF    !switch for surface reflectance
   integer, intent(in) :: &
       NFILT1(MXK), &                !filter function for auto. cell splitting
       MNK                           !number of internal layers
   real(double), intent(in) :: EPST  !determines trunction of Gauss-Seidel iteration series
   real(double), intent(in) :: &
       UBEAM,&                       !cosine of beam zenith angle
       OMEGA(nrt),&               !single scattering albedo
       EXPMU(MXHALF,MXK), &          !attenuation factors EXP(-TAU(K)/MU(I))
       EXPBEAM(MXK),&                !attenuation factors EXP(-TAU(K)/UBEAM)
       EXPINT(0:MXK)                 !attenuation factors EXP(-TAUINT(K)/UBEAM)
   real(double), intent(in) :: &
      BDRF_MS(NSTOKES, NSTOKES, MXHALF+2, MXHALF+2, 0:MAXSTR-1) 
!c                                   !reflection matrix, Fourier decomposition
!C     OUTPUT
!C     INTENSITY FIELD
   real(double), intent(out) :: &
      RAD_UP(NSTOKES,MXHALF,0:MXK),& !Fourier comp. of upward intensity  
      RAD_DN(NSTOKES,MXHALF,0:MXK)   !Fourier comp. of downward intensity
!*** local variables
   real(double) :: &                 !single scattering solution
      RSS_UP(NSTOKES,MXHALF,MXK), &  !upward direction
      RSS_DN(NSTOKES,MXHALF,MXK)     !downward direction
   real(double) :: FLXBEAM           !direct downward flux at the surface
   integer :: nmat, ifwd_adj, isrc
!------------------------------------------------------------------------------ 
!*** initialization of the gauss-seidel iteration by single scattering
      call SSC_LAYER_POL(&
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
           EXPINT,&
           RSS_UP,&
           RSS_DN,&
           ISRC,& 
           NFILT1)

      call GAUSSEID_INIT_SSC_POL(&
          NF,        MNK,       IFWD_ADJ,  ISRC,&
          EXPMU,     FLXBEAM,   SW_SURF,   BDRF_MS,&
          RSS_UP,    RSS_DN,    RAD_UP,    RAD_DN)

!*** matrix inversion using Gauss-Seidel approach
      call GAUITR_MAT(&
           nrt, &
           NF,&
           MNK, &  
           IFWD_ADJ,& 
           ISRC,&
           NMAT,&
           NFILT1,& 
           EPST,&   
           UBEAM,&
           FLXBEAM, &
           SW_SURF,& 
           pm,&
           gauss_quad,&
           BDRF_MS,&
           OMEGA,&
           EXPMU,&
           EXPINT, &
           RSS_UP,& 
           RSS_DN, &
           RAD_UP,&
           RAD_DN,&  
           ncoefs)

      return
    end subroutine GAUSS_SEIDEL

!------------------------------------------------------------------------------ 

end module gaus_pol_module
