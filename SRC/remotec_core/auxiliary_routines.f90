module auxiliary_routines_module
   use header_module
   implicit none
   private

   public :: mat_inverse, inverse_lu, INVERSE_3BY3, intrpl_spline, U0_KASTEN_AND_YOUNG, stdv, check, linear_interpol
   private :: inf_trans, E01BAF, E02BBF

contains

!------------------------------------------------------------------------------
!     F. Kasten and T. Young, revised optical air mass tabels and
!     approximation formula (1989) Apllied Optics Vol 28, No. 22 p. 4735
!     and J. Lenoble, atmospheric radiative transfer (1993), p. 236
!------------------------------------------------------------------------------
   subroutine U0_KASTEN_AND_YOUNG(SZA, u0)
      real(double), intent(in) :: sza
      real(double), intent(out) :: u0
!*** local variables
      real(double) :: a, b, c, gamma
!------------------------------------------------------------------------------
      A = 0.50572
      B = 6.07995
      C = 1.63640

      GAMMA = 90 - SZA
      U0 = sin(GAMMA*pi/180.) + A*(GAMMA + B)**(-C)
      U0 = min(1.D0, U0)

      return
   end subroutine U0_KASTEN_AND_YOUNG

!------------------------------------------------------------------------------

!> Matrix inversion though singular value decomposition
   function MAT_INVERSE(a, ny, nx)
      integer i, j, iostat
      integer :: nx, ny
      real(double), dimension(nx, ny) :: mat_inverse
      real(double), dimension(nx, ny) :: u_tr
      real(double), dimension(ny, nx) :: a
      real(double), dimension(ny, nx) :: u
      real(double), dimension(nx) :: w
      real(double), dimension(nx, nx) :: v, w_inv
      real(double) :: wmax, wmin, min_wmin, cond_nmb
!------------------------------------------------------------------------------
      if (nx .gt. ny) then
         do i = nx + 1, ny
            do j = 1, nx
               a(i, j) = 0.
            end do
         end do
      end if

      do i = 1, ny
         do j = 1, nx
            u(i, j) = a(i, j)
         end do
      end do
      call dsvdcmp(u, ny, nx, ny, nx, w, v, iostat)
      if (iostat < 0) call writelog('MAT_INVERSE: singular matrix in SVD', 8)
      wmax = 0.
      do j = 1, nx
         if (w(j) .gt. wmax) wmax = w(j)
      end do
      wmin = wmax
      do j = 1, nx
         if (w(j) .lt. wmin) wmin = w(j)
      end do
      cond_nmb = wmax/wmin
      min_wmin = wmax*1.0D-20

      do j = 1, nx
         if (w(j) .lt. min_wmin) w(j) = 0.
      end do

      do i = 1, nx
         do j = 1, nx
            w_inv(i, j) = 0.D0
            if (w(i) .ne. 0.) w_inv(i, i) = 1/w(i)
         end do
      end do
      u_tr = transpose(u)
      mat_inverse = matmul(v, matmul(w_inv, u_tr))
      return

   end function MAT_INVERSE

!------------------------------------------------------------------------------
!> Find inverse of square matrix through LU decomposition
   subroutine inverse_lu(a, m, inverse, ierr)
      integer, intent(in) :: m
      real(double), dimension(m, m), intent(in) :: a
      real(double), dimension(m, m), intent(out) :: inverse
      integer, intent(out) :: ierr
!*** local
      integer, dimension(m) :: indx
      real(double) :: d, a_lu(m, m)
      integer :: i, j
!---------------------------------------------------------------------------

!*** set up identity matrix
      do i = 1, m
         do j = 1, m
            inverse(i, j) = 0.d0
         end do
         inverse(i, i) = 1.d0
      end do

      a_lu = a
      call ludcmp(a_lu, m, m, indx, d, ierr)
!HH: added ierr to identify singular matrices
!Previously, this resulted in a pause statement, now the run can continue
      if (ierr == -1) then
!This will ensure non-convergence for current retrieval
         call writelog('INVERSE_LU: Singular matrix in LU decomposition', 8)
         inverse(1:m, 1:m) = 0.d0
         return
      end if

      do J = 1, M
         call lubksb(a_lu, m, m, indx, inverse(:, j))
      end do
      return

   end subroutine inverse_lu

!--------------------------------------------------------------------------------
!> Find inverse of 3x3 matrix
   function INVERSE_3BY3(A)
      real(double), dimension(3, 3) :: INVERSE_3BY3 !output
      real(double) :: A(3, 3)          ! input
      real(double) :: Adet

      INVERSE_3BY3(1, 1) = A(2, 2)*A(3, 3) - A(3, 2)*A(2, 3)

      INVERSE_3BY3(1, 2) = A(3, 2)*A(1, 3) - A(1, 2)*A(3, 3)

      INVERSE_3BY3(1, 3) = A(1, 2)*A(2, 3) - A(1, 3)*A(2, 2)

      Adet = INVERSE_3BY3(1, 1)*A(1, 1) + INVERSE_3BY3(1, 2)*A(2, 1) + INVERSE_3BY3(1, 3)*A(3, 1)

      INVERSE_3BY3(1, 1) = INVERSE_3BY3(1, 1)/Adet

      INVERSE_3BY3(1, 2) = INVERSE_3BY3(1, 2)/Adet

      INVERSE_3BY3(1, 3) = INVERSE_3BY3(1, 3)/Adet

      INVERSE_3BY3(2, 1) = (A(2, 3)*A(3, 1) - A(2, 1)*A(3, 3))/Adet

      INVERSE_3BY3(2, 2) = (A(1, 1)*A(3, 3) - A(3, 1)*A(1, 3))/Adet

      INVERSE_3BY3(2, 3) = (A(2, 1)*A(1, 3) - A(1, 1)*A(2, 3))/Adet

      INVERSE_3BY3(3, 1) = (A(2, 1)*A(3, 2) - A(2, 2)*A(3, 1))/Adet

      INVERSE_3BY3(3, 2) = (A(3, 1)*A(1, 2) - A(1, 1)*A(3, 2))/Adet

      INVERSE_3BY3(3, 3) = (A(1, 1)*A(2, 2) - A(1, 2)*A(2, 1))/Adet

      return     ! end of invert3x3 matrix code

   end function INVERSE_3BY3
!--------------------------------------------------------------------------------

   function INF_TRANS(flag, lobd, upbd, x)

      real(double)  :: INF_TRANS

      integer :: flag
      real(double)  :: lobd, upbd
      real(double)  :: x
      real(double) :: x_out

      !*** Transform x with boundaries [lobd,upbd] into [-INF,INF]
      if (flag == 1) then

         x_out = tan(Pi/2.d0*(2.d0*(x - lobd)/(upbd - lobd) - 1.d0))

         !*** Transform x with boundaries [-INF,INF] into [lobd,upbd]
      elseif (flag == -1) then

         x_out = (upbd - lobd)/Pi*atan(x) + upbd/2.+lobd/2.

         !*** Derivative of x with boundaries [lobd,upbd] wrt x with boundaries [-INF,INF]
      elseif (flag == 0) then

         x_out = (upbd - lobd)/Pi/(1 + tan(Pi/2.d0*(2.d0*(x - lobd)/(upbd - lobd) - 1.d0))**2)
      end if

      INF_TRANS = x_out

   end function INF_TRANS
!*************************************************************************
! Interpolation based on NAG Library routines
   subroutine intrpl_spline(MM, XX, YY, XARG, YFIT, key, KS, CS, ifail)
      !*** input
      integer, intent(IN)    :: MM
      real(double), intent(IN)    :: XX(MM), YY(MM), XARG
      !*** input/output
      integer, intent(INOUT) :: key
      real(double), intent(INOUT) :: KS(MM + 4), CS(MM + 4)
      !*** output
      real(double), intent(OUT)   :: YFIT
      integer, intent(out) :: ifail
      !*** local
      integer   :: LCK, LWRK
      real(double) ::  WRK(6*MM + 16)

      if (key .eq. 0) then
         LCK = MM + 4
         LWRK = 6*MM + 16
         !          IFAIL=0   ! hard failure (not thread-safe)
         IFAIL = 1   ! soft failure
         call E01BAF(MM, XX, YY, KS, CS, LCK, WRK, LWRK, IFAIL)
         if (ifail .ne. 0) then
            call stopretrieval('INTRPL_SLPINE: failure in NAG spline interpolation routine')
            return
         end if
         key = 1
      end if

      !          IFAIL=0
      IFAIL = 1   ! soft failure
      call E02BBF(MM + 4, KS, CS, XARG, YFIT, IFAIL)
      if (ifail .ne. 0) call stopretrieval('INTRPL_SPLINE: failure in NAG spline interpolation routine')

      return
   end subroutine intrpl_spline

!*******************SPLINE**INTEREPOLATION****ROUTINES**********************************************
   subroutine E01BAF(M, X, Y, K, C, LCK, WRK, LWRK, IFAIL)
!     MARK 8 RELEASE. NAG COPYRIGHT 1979.
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!
!     ******************************************************
!
!     NPL ALGORITHMS LIBRARY ROUTINE SP3INT
!
!     CREATED 16/5/79.                        RELEASE 00/00
!
!     AUTHORS ... GERALD T. ANTHONY, MAURICE G.COX
!                 J.GEOFFREY HAYES AND MICHAEL A. SINGER.
!     NATIONAL PHYSICAL LABORATORY, TEDDINGTON,
!     MIDDLESEX TW11 OLW, ENGLAND
!
!     ******************************************************
!
!     E01BAF.  AN ALGORITHM, WITH CHECKS, TO DETERMINE THE
!     COEFFICIENTS IN THE B-SPLINE REPRESENTATION OF A CUBIC
!     SPLINE WHICH INTERPOLATES (PASSES EXACTLY THROUGH) A
!     GIVEN SET OF POINTS.
!
!     INPUT PARAMETERS
!        M        THE NUMBER OF DISTINCT POINTS WHICH THE
!                    SPLINE IS TO INTERPOLATE.
!                    (M MUST BE AT LEAST 4.)
!        X        ARRAY CONTAINING THE DISTINCT VALUES OF THE
!                    INDEPENDENT VARIABLE. NB X(I) MUST BE
!                    STRICTLY GREATER THAN X(J) WHENEVER I IS
!                    STRICTLY GREATER THAN J.
!        Y        ARRAY CONTAINING THE VALUES OF THE DEPENDENT
!                    VARIABLE.
!        LCK      THE SMALLER OF THE ACTUALLY DECLARED DIMENSIONS
!                    OF K AND C. MUST BE AT LEAST M + 4.
!
!     OUTPUT PARAMETERS
!        K        ON SUCCESSFUL EXIT, K CONTAINS THE KNOTS
!                    SET UP BY THE ROUTINE. IF THE SPLINE IS
!                    TO BE EVALUATED (BY NPL ROUTINE E02BEF,
!                    FOR EXAMPLE) THE ARRAY K MUST NOT BE
!                    ALTERED BEFORE CALLING THAT ROUTINE.
!       C        ON SUCCESSFUL EXIT, C CONTAINS THE B-SPLINE
!                    COEFFICIENTS OF THE INTERPOLATING SPLINE.
!                    THESE ARE ALSO REQUIRED BY THE EVALUATING
!                    ROUTINE E02BEF.
!        IFAIL    FAILURE INDICATOR
!                    0 - SUCCESSFUL TERMINATION.
!                    1 - ONE OF THE FOLLOWING CONDITIONS HAS
!                        BEEN VIOLATED -
!                        M AT LEAST 4
!                        LK AT LEAST M + 4
!                        LWORK AT LEAST 6 * M + 16
!                    2 - THE VALUES OF THE INDEPENDENT VARIABLE
!                        ARE DISORDERED. IN OTHER WORDS, THE
!                        CONDITION MENTIONED UNDER X IS NOT
!                        SATISFIED.
!
!     WORKSPACE (AND ASSOCIATED DIMENSION) PARAMETERS
!        WRK     WORKSPACE ARRAY, OF LENGTH LWRK.
!        LWRK    ACTUAL DECLARED DIMENSION OF WRK.
!                    MUST BE AT LEAST 6 * M + 16.
!
!     .. Parameters ..
      character*6 SRNAME
      parameter(SRNAME='E01BAF')
!     .. Scalar Arguments ..
      integer IFAIL, LCK, LWRK, M
!     .. Array Arguments ..
      real(double) C(LCK), K(LCK), WRK(LWRK), X(M), Y(M)
!     .. Local Scalars ..
      real(double) SS
      integer I, IERROR, M1, M2
!     .. Local Arrays ..
!      character*1       P01REC(1)
!     .. External Functions ..
!      integer           P01ABF
!      external          P01ABF
!     .. External Subroutines ..
!      EXTERNAL          E02BAF
!     .. Data statements ..
      real(double), parameter :: ONE = 1.0D+0
!     .. Executable Statements ..
      IERROR = 1
!
!     TESTS FOR ADEQUACY OF ARRAY LENGTHS AND THAT M IS GREATER
!     THAN 4.
!
      if (LWRK .lt. 6*M + 16 .or. M .lt. 4) GO TO 80
      if (LCK .lt. M + 4) GO TO 80
!
!     TESTS FOR THE CORRECT ORDERING OF THE X(I)
!
      IERROR = 2
      do 20 I = 2, M
         if (X(I) .le. X(I - 1)) GO TO 80
20       continue

!     INITIALISE THE ARRAY OF KNOTS AND THE ARRAY OF WEIGHTS

         WRK(1) = ONE
         WRK(2) = ONE
         WRK(3) = ONE
         WRK(4) = ONE
         if (M .eq. 4) GO TO 60
         do 40 I = 5, M
            K(I) = X(I - 2)
            WRK(I) = ONE
40          continue
60          M1 = M + 1
            M2 = M1 + M
!
!     CALL THE SPLINE FITTING ROUTINE
!
            IERROR = 0
            call E02BAF(M, M + 4, X, Y, WRK, K, WRK(M1), WRK(M2), C, SS, IERROR)
!
!     ALL THE TESTS PERFORMED BY E02BAF ARE REDUNDANT
!     BECAUSE OF THE ABOVE TESTS AND ASSIGNMENTS, AND SO
!     IERROR = 0 AFTER THIS CALL.

80          IFAIL = ierror !P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
            return
!
!     END OF E01BAF.
!
            end subroutine E01BAF
!*****************************************************************************
            subroutine E02BAF(M, NCAP7, X, Y, W, K, WORK1, WORK2, C, SS, IFAIL)
!     NAG COPYRIGHT 1975
!     MARK 5 RELEASE
!     MARK 6 REVISED  IER-84
!     MARK 8 RE-ISSUE. IER-224 (APR 1980).
!     MARK 9A REVISED. IER-356 (NOV 1981)
!     MARK 11.5(F77) REVISED. (SEPT 1985.)
!
!     NAG LIBRARY SUBROUTINE  E02BAF
!
!     E02BAF  COMPUTES A WEIGHTED LEAST-SQUARES APPROXIMATION
!     TO AN ARBITRARY SET OF DATA POINTS BY A CUBIC SPLINE
!     WITH KNOTS PRESCRIBED BY THE USER.  CUBIC SPLINE
!     INTERPOLATION CAN ALSO BE CARRIED OUT.
!
!     COX-DE BOOR METHOD FOR EVALUATING B-SPLINES WITH
!     ADAPTATION OF GENTLEMAN*S PLANE ROTATION SCHEME FOR
!     SOLVING OVER-DETERMINED LINEAR SYSTEMS.
!
!     USES NAG LIBRARY ROUTINE  P01AAF.
!
!     STARTED - 1973.
!     COMPLETED - 1976.
!     AUTHOR - MGC AND JGH.
!
!     REDESIGNED TO USE CLASSICAL GIVENS ROTATIONS IN
!     ORDER TO AVOID THE OCCASIONAL UNDERFLOW (AND HENCE
!     OVERFLOW) PROBLEMS EXPERIENCED BY GENTLEMAN*S 3-
!     MULTIPLICATION PLANE ROTATION SCHEME
!
!     WORK1  AND  WORK2  ARE WORKSPACE AREAS.
!     WORK1(R)  CONTAINS THE VALUE OF THE  R TH  DISTINCT DATA
!     ABSCISSA AND, SUBSEQUENTLY, FOR  R = 1, 2, 3, 4,  THE
!     VALUES OF THE NON-ZERO B-SPLINES FOR EACH SUCCESSIVE
!     ABSCISSA VALUE.
!     WORK2(L, J)  CONTAINS, FOR  L = 1, 2, 3, 4,  THE VALUE OF
!     THE  J TH  ELEMENT IN THE  L TH  DIAGONAL OF THE
!     UPPER TRIANGULAR MATRIX OF BANDWIDTH  4  IN THE
!     TRIANGULAR SYSTEM DEFINING THE B-SPLINE COEFFICIENTS.
!
!     .. Parameters ..
               character*6 SRNAME
               parameter(SRNAME='E02BAF')
!     .. Scalar Arguments ..
               real(double) SS
               integer IFAIL, M, NCAP7
!     .. Array Arguments ..
               real(double) C(NCAP7), K(NCAP7), W(M), WORK1(M), &
                  WORK2(4, NCAP7), X(M), Y(M)
!     .. Local Scalars ..
               real(double) ACOL, AROW, CCOL, COSINE, CROW, D, D4, D5, D6, &
                  D7, D8, D9, DPRIME, E2, E3, E4, E5, K0, K1, K2, &
                  K3, K4, K5, K6, N1, N2, N3, RELEMT, S, SIGMA, &
                  SINE, WI, XI
               integer I, IERROR, IPLUSJ, IU, J, JOLD, JPLUSL, JREV, L, &
                  L4, LPLUS1, LPLUSU, NCAP, NCAP3, NCAPM1, R
!     .. Local Arrays ..
!      character*1       P01REC(1)
!     .. External Functions ..
!      integer           P01ABF
!      external          P01ABF
!     .. Intrinsic Functions ..
               intrinsic ABS, SQRT
!     .. Executable Statements ..
               IERROR = 4
!     CHECK THAT THE VALUES OF  M  AND  NCAP7  ARE REASONABLE
               if (NCAP7 .lt. 8 .or. M .lt. NCAP7 - 4) GO TO 420
               NCAP = NCAP7 - 7
               NCAPM1 = NCAP - 1
               NCAP3 = NCAP + 3

!     IN ORDER TO DEFINE THE FULL B-SPLINE BASIS, AUGMENT THE
!     PRESCRIBED INTERIOR KNOTS BY KNOTS OF MULTIPLICITY FOUR
!     AT EACH END OF THE DATA RANGE.

               do 20 J = 1, 4
                  I = NCAP3 + J
                  K(J) = X(1)
                  K(I) = X(M)
20                continue

!     TEST THE VALIDITY OF THE DATA.

!     CHECK THAT THE KNOTS ARE ORDERED AND ARE INTERIOR
!     TO THE DATA INTERVAL.

                  IERROR = 1
                  if (K(5) .le. X(1) .or. K(NCAP3) .ge. X(M)) GO TO 420
                  do 40 J = 4, NCAP3
                     if (K(J) .gt. K(J + 1)) GO TO 420
40                   continue

!     CHECK THAT THE WEIGHTS ARE STRICTLY POSITIVE.

                     IERROR = 2
                     do 60 I = 1, M
                        if (W(I) .le. 0.0D0) GO TO 420
60                      continue

!     CHECK THAT THE DATA ABSCISSAE ARE ORDERED, THEN FORM THE
!     ARRAY  WORK1  FROM THE ARRAY  X.  THE ARRAY  WORK1  CONTAINS
!     THE
!     SET OF DISTINCT DATA ABSCISSAE.

                        IERROR = 3
                        WORK1(1) = X(1)
                        J = 2
                        do 80 I = 2, M
                           if (X(I) .lt. WORK1(J - 1)) GO TO 420
                           if (X(I) .eq. WORK1(J - 1)) GO TO 80
                           WORK1(J) = X(I)
                           J = J + 1
80                         continue
                           R = J - 1

!     CHECK THAT THERE ARE SUFFICIENT DISTINCT DATA ABSCISSAE FOR
!     THE PRESCRIBED NUMBER OF KNOTS.

                           IERROR = 4
                           if (R .lt. NCAP3) GO TO 420

!     CHECK THE FIRST  S  AND THE LAST  S  SCHOENBERG-WHITNEY
!     CONDITIONS ( S = MIN(NCAP - 1, 4) ).

                           IERROR = 5
                           do 100 J = 1, 4
                              if (J .ge. NCAP) GO TO 160
                              I = NCAP3 - J + 1
                              L = R - J + 1
                              if (WORK1(J) .ge. K(J + 4) .or. K(I) .ge. WORK1(L)) GO TO 420
100                           continue

!     CHECK ALL THE REMAINING SCHOENBERG-WHITNEY CONDITIONS.

                              if (NCAP .le. 5) GO TO 160
                              R = R - 4
                              I = 4
                              do 140 J = 5, NCAPM1
                                 K0 = K(J + 4)
                                 K4 = K(J)
120                              I = I + 1
                                 if (WORK1(I) .le. K4) GO TO 120
                                 if (I .gt. R .or. WORK1(I) .ge. K0) GO TO 420
140                              continue

!     INITIALISE A BAND TRIANGULAR SYSTEM (I.E. A
!     MATRIX AND A RIGHT HAND SIDE) TO ZERO. THE
!     PROCESSING OF EACH DATA POINT IN TURN RESULTS
!     IN AN UPDATING OF THIS SYSTEM. THE SUBSEQUENT
!     SOLUTION OF THE RESULTING BAND TRIANGULAR SYSTEM
!     YIELDS THE COEFFICIENTS OF THE B-SPLINES.

160                              do 200 I = 1, NCAP3
                                    do 180 L = 1, 4
                                       WORK2(L, I) = 0.0D0
180                                    continue
                                       C(I) = 0.0D0
200                                    continue
                                       SIGMA = 0.0D0
                                       J = 0
                                       JOLD = 0
                                       do 340 I = 1, M

!        FOR THE DATA POINT  (X(I), Y(I))  DETERMINE AN INTERVAL
!        K(J + 3) .LE. X .LT. K(J + 4)  CONTAINING  X(I).  (IN THE
!        CASE  J + 4 .EQ. NCAP  THE SECOND EQUALITY IS RELAXED TO
!        INCLUDE
!        EQUALITY).

                                          WI = W(I)
                                          XI = X(I)
220                                       if (XI .lt. K(J + 4) .or. J .gt. NCAPM1) GO TO 240
                                          J = J + 1
                                          GO TO 220
240                                       if (J .eq. JOLD) GO TO 260

!        SET CERTAIN CONSTANTS RELATING TO THE INTERVAL
!        K(J + 3) .LE. X .LE. K(J + 4).

                                          K1 = K(J + 1)
                                          K2 = K(J + 2)
                                          K3 = K(J + 3)
                                          K4 = K(J + 4)
                                          K5 = K(J + 5)
                                          K6 = K(J + 6)
                                          D4 = 1.0D0/(K4 - K1)
                                          D5 = 1.0D0/(K5 - K2)
                                          D6 = 1.0D0/(K6 - K3)
                                          D7 = 1.0D0/(K4 - K2)
                                          D8 = 1.0D0/(K5 - K3)
                                          D9 = 1.0D0/(K4 - K3)
                                          JOLD = J

!        COMPUTE AND STORE IN  WORK1(L) (L = 1, 2, 3, 4)  THE VALUES
!        OF
!        THE FOUR NORMALIZED CUBIC B-SPLINES WHICH ARE NON-ZERO AT
!        X=X(I).

260                                       E5 = K5 - XI
                                          E4 = K4 - XI
                                          E3 = XI - K3
                                          E2 = XI - K2
                                          N1 = WI*D9
                                          N2 = E3*N1*D8
                                          N1 = E4*N1*D7
                                          N3 = E3*N2*D6
                                          N2 = (E2*N1 + E5*N2)*D5
                                          N1 = E4*N1*D4
                                          WORK1(4) = E3*N3
                                          WORK1(3) = E2*N2 + (K6 - XI)*N3
                                          WORK1(2) = (XI - K1)*N1 + E5*N2
                                          WORK1(1) = E4*N1
                                          CROW = Y(I)*WI

!        ROTATE THIS ROW INTO THE BAND TRIANGULAR SYSTEM USING PLANE
!        ROTATIONS.

                                          do 320 LPLUS1 = 1, 4
                                             L = LPLUS1 - 1
                                             RELEMT = WORK1(LPLUS1)
                                             if (RELEMT .eq. 0.0D0) GO TO 320
                                             JPLUSL = J + L
                                             L4 = 4 - L
                                             D = WORK2(1, JPLUSL)
                                             if (abs(RELEMT) .ge. D) DPRIME = abs(RELEMT) &
                                                                              *sqrt(1.0D0 + (D/RELEMT)**2)
                                             if (abs(RELEMT) .lt. D) DPRIME = D*sqrt(1.0D0 + (RELEMT/D)**2)
                                             WORK2(1, JPLUSL) = DPRIME
                                             COSINE = D/DPRIME
                                             SINE = RELEMT/DPRIME
                                             if (L4 .lt. 2) GO TO 300
                                             do 280 IU = 2, L4
                                                LPLUSU = L + IU
                                                ACOL = WORK2(IU, JPLUSL)
                                                AROW = WORK1(LPLUSU)
                                                WORK2(IU, JPLUSL) = COSINE*ACOL + SINE*AROW
                                                WORK1(LPLUSU) = COSINE*AROW - SINE*ACOL
280                                             continue
300                                             CCOL = C(JPLUSL)
                                                C(JPLUSL) = COSINE*CCOL + SINE*CROW
                                                CROW = COSINE*CROW - SINE*CCOL
320                                             continue
                                                SIGMA = SIGMA + CROW**2
340                                             continue
                                                SS = SIGMA

!     SOLVE THE BAND TRIANGULAR SYSTEM FOR THE B-SPLINE
!     COEFFICIENTS. IF A DIAGONAL ELEMENT IS ZERO, AND HENCE
!     THE TRIANGULAR SYSTEM IS SINGULAR, THE IMPLICATION IS
!     THAT THE SCHOENBERG-WHITNEY CONDITIONS ARE ONLY JUST
!     SATISFIED. THUS IT IS APPROPRIATE TO EXIT IN THIS
!     CASE WITH THE SAME VALUE  (IFAIL=5)  OF THE ERROR
!     INDICATOR.

                                                L = -1
                                                do 400 JREV = 1, NCAP3
                                                   J = NCAP3 - JREV + 1
                                                   D = WORK2(1, J)
                                                   if (D .eq. 0.0D0) GO TO 420
                                                   if (L .lt. 3) L = L + 1
                                                   S = C(J)
                                                   if (L .eq. 0) GO TO 380
                                                   do 360 I = 1, L
                                                      IPLUSJ = I + J
                                                      S = S - WORK2(I + 1, J)*C(IPLUSJ)
360                                                   continue
380                                                   C(J) = S/D
400                                                   continue
                                                      IERROR = 0
420                                                   if (IERROR) 440, 460, 440
440                                                   IFAIL = ierror !P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
                                                      return
460                                                   IFAIL = 0
                                                      return
                                                      end subroutine E02BAF
!*********************************************************************************
                                                      subroutine E02BBF(NCAP7, K, C, X, S, IFAIL)
!     NAG LIBRARY SUBROUTINE  E02BBF

!     E02BBF  EVALUATES A CUBIC SPLINE FROM ITS
!     B-SPLINE REPRESENTATION.

!     DE BOOR*S METHOD OF CONVEX COMBINATIONS.

!     USES NAG LIBRARY ROUTINE  P01AAF.

!     STARTED - 1973.
!     COMPLETED - 1976.
!     AUTHOR - MGC AND JGH.

!     NAG COPYRIGHT 1975
!     MARK 5 RELEASE
!     MARK 7 REVISED IER-141 (DEC 1978)
!     MARK 11.5(F77) REVISED. (SEPT 1985.)

!     .. Parameters ..
                                                         character*6 SRNAME
                                                         parameter(SRNAME='E02BBF')
!     .. Scalar Arguments ..
                                                         real(double) S, X
                                                         integer IFAIL, NCAP7
!     .. Array Arguments ..
                                                         real(double) C(NCAP7), K(NCAP7)
!     .. Local Scalars ..
                                                         real(double) C1, C2, C3, E2, E3, E4, E5, K1, K2, K3, K4, K5, K6
                                                         integer IERROR, J, J1, L
!     .. Local Arrays ..
!      character*1       P01REC(1)
!     .. External Functions ..
!      integer           P01ABF
!      external          P01ABF
!     .. Executable Statements ..
                                                         IERROR = 0
                                                         if (NCAP7 .ge. 8) GO TO 20
                                                         IERROR = 2
                                                         GO TO 120
20                                                       if (X .ge. K(4) .and. X .le. K(NCAP7 - 3)) GO TO 40
                                                         IERROR = 1
                                                         S = 0.0D0
                                                         GO TO 120

!     DETERMINE  J  SUCH THAT  K(J + 3) .LE. X .LE. K(J + 4).

40                                                       J1 = 0
                                                         J = NCAP7 - 7
60                                                       L = (J1 + J)/2
                                                         if (J - J1 .le. 1) GO TO 100
                                                         if (X .ge. K(L + 4)) GO TO 80
                                                         J = L
                                                         GO TO 60
80                                                       J1 = L
                                                         GO TO 60

!     USE THE METHOD OF CONVEX COMBINATIONS TO COMPUTE  S(X).

100                                                      K1 = K(J + 1)
                                                         K2 = K(J + 2)
                                                         K3 = K(J + 3)
                                                         K4 = K(J + 4)
                                                         K5 = K(J + 5)
                                                         K6 = K(J + 6)
                                                         E2 = X - K2
                                                         E3 = X - K3
                                                         E4 = K4 - X
                                                         E5 = K5 - X
                                                         C2 = C(J + 1)
                                                         C3 = C(J + 2)
                                                         C1 = ((X - K1)*C2 + E4*C(J))/(K4 - K1)
                                                         C2 = (E2*C3 + E5*C2)/(K5 - K2)
                                                         C3 = (E3*C(J + 3) + (K6 - X)*C3)/(K6 - K3)
                                                         C1 = (E2*C2 + E4*C1)/(K4 - K2)
                                                         C2 = (E3*C3 + E5*C2)/(K5 - K3)
                                                         S = (E3*C2 + E4*C1)/(K4 - K3)
120                                                      if (IERROR) 140, 160, 140
140                                                      IFAIL = ierror !P01ABF(IFAIL,IERROR,SRNAME,0,P01REC)
                                                         return
160                                                      IFAIL = 0
                                                         return
                                                      end subroutine E02BBF

                                                      !-----------------------------------------------------------------------------
                                                      !> Mass interpolation routine.
  !!
  !! This routine interpolates a given function represented as an x-y table on
  !! a number of points, returning an array of interpolated values. Both arrays of the
  !! independent (x) values of the function and the points on which the interpolation
  !! has to be performed should be strictly monotonous. However, they can be either
  !! ascending or descending. And it is also possible to have one ascending array
  !! and the other one descending.
                                                      subroutine linear_interpol(sz1, sz2, z1, z2, f1, f2) ! {{{

                                                         implicit none

                                                         ! Input and output.
                                                         integer, intent(in) :: sz1 !< Size of the x-y table of the function to be interpolated.
                                                         integer, intent(in) :: sz2 !< Size of the array of interpolated values.
                                                         real, dimension(sz1), intent(in) :: z1 !< Independent (x) values of the function to be interpolated (should be monotonous).
                                                         real, dimension(sz2), intent(in) :: z2 !< Independent (x) values on which the interpolation should take place (should be monotonous).
                                                         real, dimension(sz1), intent(in) :: f1 !< Dependent (y) values of the function to be interpolated.
                                                         real, dimension(sz2), intent(out) :: f2 !< Interpolation values (y).

                                                         real, dimension(sz1) :: x1
                                                         real, dimension(sz1) :: y1
                                                         real, dimension(sz2) :: x2
                                                         real, dimension(sz2) :: y2

                                                         real :: a
                                                         real :: b

                                                         integer :: k1
                                                         integer :: k2
                                                         integer :: m1
                                                         integer :: m2

                                                         ! z1,f1 and z2,f2 may go from top to bottom, x1,y1 and x2,y2
                                                         ! go from bottom to top.

                                                         ! z1,f1 and z2,f2 must be either ascending or descending, not
                                                         ! a combination of both.

                                                         ! Defines ascending arrays x1 and y1. These are either copies of
                                                         ! z1 and f1, or their reverses, depending on the order of z1.
                                                         if (z1(1) .gt. z1(sz1)) then
                                                            do k1 = 1, sz1
                                                               x1(k1) = z1(sz1 + 1 - k1)
                                                               y1(k1) = f1(sz1 + 1 - k1)
                                                            end do
                                                         else
                                                            x1 = z1
                                                            y1 = f1
                                                         end if

                                                         ! The same for x2 with respect of z2, only there is no y2 or f2 yet.
                                                         ! The array x2 will be ascending. Later on, the interpolated array
                                                         ! may be reversed, depending on whether x2 is a reverse of z2 or not.
                                                         if (z2(1) .gt. z2(sz2)) then
                                                            do k2 = 1, sz2
                                                               x2(k2) = z2(sz2 + 1 - k2)
                                                            end do
                                                         else
                                                            x2 = z2
                                                         end if

                                                         ! Interpolation indices. Initialize at lower index at 1, which means
                                                         ! that the value should be between 1 and 2. These values may increase
                                                         ! later on. Even if the value on which is interpolated is outside
                                                         ! domain x1, it will be interpreted as if it is between the last two
                                                         ! indices at that side.
                                                         m1 = 1

                                                         do k2 = 1, sz2
                                                            ! Update interpolation index m2 if necessary, but not further than
                                                            ! sz1. As x2 is increasing, the interpolation indices are expected
                                                            ! to increase steadily. At least, they never need to decrease.
                                                            do m2 = m1 + 1, sz1 - 1 ! So that sz1 is the end result if x1 is very large.
                                                               if (x2(k2) .le. x1(m2)) exit
                                                            end do
                                                            ! Lower index is always one lower than the upper index.
                                                            m1 = m2 - 1

                                                            ! Linear relationship.
                                                            a = (y1(m1) - y1(m2))/(x1(m1) - x1(m2))
                                                            b = y1(m1) - a*x1(m1)

                                                            y2(k2) = a*x2(k2) + b
                                                         end do

                                                         ! Reverse the result if the target domain was reversed at the start.
                                                         if (z2(1) .gt. z2(sz2)) then
                                                            do k2 = 1, sz2
                                                               f2(k2) = y2(sz2 + 1 - k2)
                                                            end do
                                                         else
                                                            f2 = y2
                                                         end if

                                                      end subroutine linear_interpol ! }}}
                                                      !-----------------------------------------------------------------------------

                                                      function stdv(y, n)
                                                         integer, intent(in) :: n
                                                         real(double), dimension(n), intent(in) :: y
                                                         real(double) :: stdv
                                                         !*** local
                                                         integer :: i
                                                         real(double) :: mean, variance

                                                         mean = 0.0 ! compute mean
                                                         do i = 1, n
                                                            mean = mean + y(i)
                                                         end do
                                                         mean = mean/n
                                                         Variance = 0.0 ! compute variance
                                                         do i = 1, n
                                                            variance = variance + (y(i) - mean)**2
                                                         end do
                                                         variance = variance/n
                                                         stdv = sqrt(variance) ! compute standard deviation
                                                         return
                                                      end function stdv

                                                      !-----------------------------------------------------------------------------
                                                      !> Check calls to netCDF functions
                                                      subroutine check(status, ierr)
                                                         use NETCDF
                                                         integer, intent(in) :: status
                                                         integer, intent(out) :: ierr

                                                         if (status /= nf90_noerr) then
                                                            ierr = ierr_read
                                                            call writelog('CHECK: '//trim(nf90_strerror(status)), 6)
                                                         else
                                                            ierr = 0
                                                         end if

                                                      end subroutine check
                                                      !-----------------------------------------------------------------------------

                                                      end module auxiliary_routines_module
