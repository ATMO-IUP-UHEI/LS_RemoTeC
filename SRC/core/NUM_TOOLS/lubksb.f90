      subroutine lubksb(a, n, np, indx, b)
      implicit none
      include '../INCLUDES/precision.inc' 
      integer, intent(in) ::  n, np, indx(n)
      real(double), intent(in) :: a(np,np)
      real(double), intent(out) :: b(n)
!*** local variables
      integer :: i, ii, j, ll
      real(double) :: sum

      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      end
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
