!------------------------------------------------------------------------------
!> In this module types (structures) are defined used in the subroutines to
!! solve the RT equation using the Gauss-Seidel method 
!------------------------------------------------------------------------------
module header_gsd
   use header_module
   implicit none
 
   integer, parameter :: MXK = 5000  !< max number of layers (incl sub layers) in the model atmosphere
   integer, parameter :: MAXITR = 25 !< max numbers of gauss seidel iterations
   integer, parameter :: maxlay = 72 !< max number layers in model atmopshere (nrt) 

!> contains the scattering phase matrix
   type phase_mat_type
!> outgoing upward, incoming downward       
     real(double), dimension(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,MAXLAY) :: PM_UP_DN 
!> outgoing upward, incoming upward
     real(double), dimension(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,MAXLAY) :: PM_UP_UP
!> outgoing downward, incoming downward
     real(double), dimension(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,MAXLAY) :: PM_DN_DN
!> outgoing downward, incoming upward
     real(double), dimension(NSTOKES,NSTOKES,MXHALF+2,MXHALF+2,MAXLAY) :: PM_DN_UP  
   end type phase_mat_type

   type gauss_quad_type
      real(double), dimension(maxstr) :: dg_mu  !< cosine of quadrature angles
      real(double), dimension(maxstr) :: dg_wt  !< quadrature weights
   end type gauss_quad_type

!> generalized spherical functions
   type gsf_type
      real(double), dimension(0:MAXSTR-1,0:MAXSTR-1,MAXSTR+4) :: plm_0 
      real(double), dimension(0:MAXSTR-1,0:MAXSTR-1,MAXSTR+4) :: plm_m
      real(double), dimension(0:MAXSTR-1,0:MAXSTR-1,MAXSTR+4) :: plm_p
   end type gsf_type

!------------------------------------------------------------------------------

end module header_gsd
