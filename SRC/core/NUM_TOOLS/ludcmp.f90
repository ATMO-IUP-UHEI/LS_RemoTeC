      subroutine ludcmp(a, n, np, indx, d, iostat)
      implicit none
      include '../INCLUDES/precision.inc'
      integer, intent(in) :: n, np
      integer, intent(out) :: indx(n), iostat
      real(double), intent(inout) :: a(np,np)
      real(double), intent(out) :: d
!*** local variables
      integer, parameter :: NMAX=2000
      real(double), parameter :: TINY=1.0e-20
      integer :: i, imax, j, k
      real(double) :: aamax, dum, sum, vv(NMAX)

      iostat = 0
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
          !pause 'singular matrix in ludcmp'
          iostat = -1
          return
        endif
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      end subroutine ludcmp
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
