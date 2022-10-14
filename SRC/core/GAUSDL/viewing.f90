module viewing_module
   use header_gsd
   implicit none


contains
!------------------------------------------------------------------------------ 
      subroutine VIEWING(&
           gauss_quad,&
           pm,&
           NF,&
           MNK,&  
           NFILT1,&  
           NMAT, &
           IFWD_ADJ,&
           SW_SURF,&
           BDRF_MS,&
           OMEGA,& 
           EXPVIEW,&
           RAD_UP,&
           RAD_DN,&
           RAD_VW)

!     inputs
      type(gauss_quad_type), intent(IN) :: gauss_quad
      type(phase_mat_type), intent(IN) :: pm
      integer :: &
          NF, &                          !Fourier index
          IFWD_ADJ, &                    !switch forward/adjoint (1,2)
          NFILT1(MXK),&
          MNK                           !number of internal layers
      logical :: &
          SW_SURF                       !switch for surface reflectance
      real(double) :: &
          BDRF_MS(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,0:MAXSTR-1) 
!c                                        !BDRF, Fourier decomposition
      real(double), dimension(:) :: OMEGA                 !single scattering albedo
      real(double) :: &
          RAD_UP(NSTOKES,MXHALF,0:MXK),& !diffuse intensity field (downward)
          RAD_DN(NSTOKES,MXHALF,0:MXK)  !diffuse intensity field (upward)
      real(double) :: &
          EXPVIEW(MXK)                  !atten. factors EXP(-tau(k)/UVIEW)
!c     output
      real(double) :: &
          RAD_VW(NSTOKES,0:MXK) !interpol. intensity field in viewing
                                 !direction
!c     internals
      real(double) :: &
          RAVG_UP(NSTOKES),&
          RAVG_DN(NSTOKES)
      real(double) ::  &            
          SRCMS_VW(NSTOKES,MXK) !multiple scattering source function 
      integer :: nmat, i, j_stream, k, nstr, i_st, ik, j_st
      real(double) :: dummy, fac, fact, one_4pi

!-----------------------------------------------------------------------
      ONE_4PI=1./(4.*PI) 

      NSTR = MXHALF + IFWD_ADJ

!     source function in viewing direction

!     initialization of SRCMS_VW

      do K = 1,MNK
         IK = NFILT1(K)

         do I_ST =1,NMAT
            RAVG_UP(I_ST) = RAD_UP(I_ST,1,K)+RAD_UP(I_ST,1,K-1)
            RAVG_DN(I_ST) = RAD_DN(I_ST,1,K)+RAD_DN(I_ST,1,K-1)
         enddo

         FACT = gauss_quad%DG_WT(1)/2
         do I_ST =1,NMAT
            SRCMS_VW(I_ST,K)=0.
            do J_ST = 1,NMAT
               SRCMS_VW(I_ST,K) = SRCMS_VW(I_ST,K) + FACT *&
                   (pm%PM_UP_UP(I_ST,J_ST,NSTR,1,IK)*RAVG_UP(J_ST) + &
                    pm%PM_UP_DN(I_ST,J_ST,NSTR,1,IK)*RAVG_DN(J_ST))
            enddo
         enddo

!        summation over I = 2,MAXHALF

         do I=2,MXHALF                     
               
            FACT = gauss_quad%DG_WT(I)/2

            do I_ST = 1,NMAT
               RAVG_UP(I_ST) = RAD_UP(I_ST,I,K)+RAD_UP(I_ST,I,K-1)
               RAVG_DN(I_ST) = RAD_DN(I_ST,I,K)+RAD_DN(I_ST,I,K-1)
            enddo

            do I_ST = 1,NMAT
               DUMMY = 0.
               do J_ST = 1,NMAT
                  DUMMY  = DUMMY + FACT*&
                   (pm%PM_UP_UP(I_ST,J_ST,NSTR,I,IK)*RAVG_UP(J_ST)+ &
                    pm%PM_UP_DN(I_ST,J_ST,NSTR,I,IK)*RAVG_DN(J_ST))
               enddo
               SRCMS_VW(I_ST,K) = SRCMS_VW(I_ST,K) + DUMMY
            enddo

         enddo

      enddo


!c---- Lambertian ground reflection of diffuse dowmward radiation ------

      do I_ST = 1,NMAT
         RAD_VW(I_ST,MNK)= 0.D0
      enddo

      J_STREAM = MXHALF+2                       !for forward simulation
      if(IFWD_ADJ.eq.2)J_STREAM = MXHALF+1      !for adjoint simulation

      if(SW_SURF) then

         do I_ST = 1,NMAT
         do I    = 1,MXHALF

            DUMMY = 2.D0*gauss_quad%DG_WT(I)*gauss_quad%DG_MU(I)
            do J_ST = 1,NMAT
               RAD_VW(I_ST,MNK) = RAD_VW(I_ST,MNK) + DUMMY *&
                   BDRF_MS(I_ST,J_ST,I,J_STREAM,NF)*RAD_DN(J_ST,I,MNK)
            enddo

         enddo
         enddo

      endif

!c---- upward direction ------------------------------------------------

      do I_ST = 1,NMAT
      do K=MNK,1,-1      

         IK = NFILT1(K)

         FAC = 0.5*OMEGA(IK)*(1 - EXPVIEW(K))

         RAD_VW(I_ST,K-1)=&
             RAD_VW(I_ST,K)   * EXPVIEW(K)&  
          +  SRCMS_VW(I_ST,K) * FAC         
      enddo
      enddo

      return
    end subroutine VIEWING

end module viewing_module
