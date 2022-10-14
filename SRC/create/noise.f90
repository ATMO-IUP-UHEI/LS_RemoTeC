
      subroutine NOISE(NINIT, Y, NDIM, ERROR)

!     put noise (gaussian) on the array Y(NDIM). The array ERROR specifies
!     the relative standard deviation of the random error, i.e. 0.01
!     means 1 %. NINIT is an integer used for initialization, and can
!     be any negative integer.

!     Modified by Andre Galli, June 4, 2010
!     CHANGES: If array(i) + noise(i) < 0 then array(i) + nose(i) = array (i)

!     Modified by Andre Galli, January 20, 2011
!     CHANGES: Extracted function ran1.f to a separate file

!     Modified by Andre Galli, March 16, 2011
!     CHANGES: If odd number of dimension, call RAND one more time to have exactly the same random pattern for every retrieval
!              if always the same random seed NINIT is used.

!--------------------------------------------------------------------------
         implicit none
         integer, intent(in) :: ninit, ndim
         double precision, intent(in) :: error(ndim)
         double precision, intent(inout) :: Y(ndim)
         integer :: CHECK, i, l
         double precision :: YT(NDIM), brange, erange, dint, rand, err
!-------------- variables to check error deviation ------------------------
         integer, parameter :: nint = 50
         double precision :: BB(NINT), EB(NINT), count(NINT)
         real :: GASDEV
         data BRANGE/0.97/, ERANGE/1.03/
!--------------------------------------------------------------------------
         CHECK = 0
         do L = 1, NDIM
            RAND = dble(GASDEV(NINIT))
            YT(L) = Y(L)*(1.0d0 + RAND*ERROR(L))
            if (YT(L) .le. 0.0d0) then
               YT(L) = Y(L) ! avoid negative values
            end if
         end do

         if (mod(ndim, 2) .eq. 1) then
            RAND = dble(GASDEV(NINIT))
         end if

         if (CHECK .eq. 1) then

!     sort the error in different bins to check the distribution

            DINT = (ERANGE - BRANGE)/real(NINT) !bin size

            do I = 1, NINT
               BB(I) = BRANGE + (I - 1)*DINT
               EB(I) = BRANGE + I*DINT
               count(I) = 0
            end do

            do L = 1, NDIM
               ERR = YT(L)/Y(L)
               do I = 1, NINT
                  if (BB(I) .le. ERR .and. ERR .lt. EB(I)) &
                     count(I) = count(I) + 1
               end do
            end do

            open (1, FILE='err_dist.dat', STATUS='unknown')
            write (1, '(8F8.4)') (0.5*(BB(I) + EB(I)), I=1, NINT)
            write (1, '(8I8)') COUNT

            close (1)

         end if

         do L = 1, NDIM
            Y(L) = YT(L)
         end do

! If an additive constant Deltay is to be added to the synthetic spectrum:

!      DO L = 1,NDIM
!         RAND  = DBLE(GASDEV(NINIT-L))
!         YT(L) = RAND * 0.01 * deltay + deltay
!      ENDDO
!
!       if (MOD(ndim,2) .EQ. 1) then
!         RAND  = DBLE(GASDEV(NINIT-L))
!       endif
!
!      DO L = 1,NDIM
!         Y(L) = Y(L) + YT(L)
!      ENDDO

         return
      end subroutine NOISE
!--------------------------------------------------------------------------

      subroutine RDN01(n, rdn)
         implicit none
!     Draw random number rdn between 0 and 1.
!--------------------------------------------------------------------------

         integer :: init, n, k
         integer*4 now(3)
         integer, dimension(8) :: values
         double precision rdn(n)
         real ran1

         !**********************************************************************
         !*** JS [29/05/2020]: init is now calculated as a sum of, hour, minute,
         !*** second *and* *millisecond*. Before the latter was included, the
         !*** radnomization did not work properly at the Mistral cluster since
         !*** blocks of scenes in the same job (lst-file) got the same random
         !*** number because the computations were initiated within the same
         !*** second.

         ! call itime(now)
         ! init=-abs(now(1)+now(2)+now(3))

         call date_and_time(VALUES=values)
         init = -abs(values(5) + values(6) + values(7) + values(8))

         !**********************************************************************

         do k = 1, n
            rdn(k) = dble(ran1(init))
         end do

         return
      end subroutine RDN01
!----------------------------------------------------------------------------

      function gasdev(idum)
         implicit none
         integer idum
         real gasdev
!    USES ran1
         integer iset
         real fac, gset, rsq, v1, v2, ran1
         save iset, gset
         data iset/0/
         if (iset .eq. 0) then
1           v1 = 2.*ran1(idum) - 1.
            v2 = 2.*ran1(idum) - 1.
            rsq = v1**2 + v2**2
            if (rsq .ge. 1. .or. rsq .eq. 0.) goto 1
            fac = sqrt(-2.*log(rsq)/rsq)
            gset = v1*fac
            gasdev = v2*fac
            iset = 1
         else
            gasdev = gset
            iset = 0
         end if
         return
      end function gasdev

!---------------------------------------------------------------------
      FUNCTION ran1(idum)
         implicit none
         INTEGER idum, IA, IM, IQ, IR, NTAB, NDIV
         REAL ran1, AM, EPS, RNMX
         PARAMETER(IA=16807, IM=2147483647, AM=1./IM, IQ=127773, IR=2836, &
                   NTAB=32, NDIV=1 + (IM - 1)/NTAB, EPS=1.2e-7, RNMX=1.-EPS)
         INTEGER j, k, iv(NTAB), iy
         SAVE iv, iy
         DATA iv/NTAB*0/, iy/0/
         if (idum .le. 0 .or. iy .eq. 0) then
            idum = max(-idum, 1)
            do 11 j = NTAB + 8, 1, -1
               k = idum/IQ
               idum = IA*(idum - k*IQ) - IR*k
               if (idum .lt. 0) idum = idum + IM
               if (j .le. NTAB) iv(j) = idum
11             continue
               iy = iv(1)
               end if
               k = idum/IQ
               idum = IA*(idum - k*IQ) - IR*k
               if (idum .lt. 0) idum = idum + IM
               j = 1 + iy/NDIV
               iy = iv(j)
               iv(j) = idum
               ran1 = min(AM*iy, RNMX)
               return
            END
!  (C) Copr. 1986-92 Numerical Recipes Software Y_'i&.
