!*************************************************************************
!*** Natural cubic spline: zero second derivatives at boundaries
subroutine spline_interpol(x, y, n, xnew, ynew, nnew, ierr)
  implicit none 
  include '../INCLUDES/precision.inc'
  integer, intent(in) :: n, nnew
  real(double), dimension(n), intent(in) :: x, y
  real(double), dimension(nnew), intent(in) :: xnew
  real(double), dimension(nnew), intent(out) :: ynew
  integer, intent(out) :: ierr
  !*** local variables
  integer :: k
  real(double), dimension(n) :: deriv2

  call spline(x, y, n, 1.1D30, 1.1D30, deriv2, ierr)
  if (ierr .ne. 0) return
  do k = 1, nnew
     call splint(x, y, deriv2, n, xnew(k), ynew(k), ierr)
     if (ierr .ne. 0) return
  enddo

end subroutine spline_interpol

!********************************************************************
!*** Cubic spline with specified first derivatives at boundaries
subroutine spline_int(x, y, n, xnew, ynew, nnew, ierr)
  implicit none 
  include '../INCLUDES/precision.inc'
  integer, intent(in) :: n, nnew
  real(double), dimension(n), intent(in) :: x, y
  real(double), dimension(nnew), intent(in) :: xnew
  real(double), dimension(nnew), intent(out) :: ynew
  integer, intent(out) :: ierr
  !*** local variables
  integer :: k
  real(double) :: yp1, yp2
  real(double), dimension(n) :: deriv2    

  ierr = 0
  yp1 = (y(2) - y(1))/(x(2) - x(1))      ! first derivative at point 1 
  yp2 = (y(n)- y(n-1))/(x(n) - x(n-1))   ! first derivative at point n
  call spline(x, y, n, yp1, yp2, deriv2, ierr)
  if (ierr .ne. 0) return
  do k = 1, nnew
     call splint(x, y, deriv2, n, xnew(k), ynew(k), ierr)
     if (ierr .ne. 0) return
  enddo

end subroutine spline_int

!************************************************************************
subroutine spline(x, y, n, yp1, ypn, y2, ierr)
  implicit none
  include '../INCLUDES/precision.inc'
  integer, intent(in) :: n !,NMAX
  real(double), intent(in) :: yp1, ypn, x(n), y(n)
  real(double), intent(out) :: y2(n)
  integer, intent(out) :: ierr
  !*** local variables
  !      PARAMETER (NMAX=4000000)  !this will cause stack overflow with openMP
  integer i, k
  real(double) p, qn, sig, un
  real(double), dimension(:), allocatable :: u

  !      if (NMAX.lt.n) stop 'Increase NMAX in subroutine spline'
  ierr = 0

  allocate(u(n), stat=ierr)
  if (ierr .ne. 0) return

  if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do 11 i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
          1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11   continue
     if (ypn.gt..99e30) then
        qn=0.
        un=0.
     else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
     endif
     y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
     do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12      continue
        deallocate(u, stat=ierr)    

        return

      end subroutine spline
      !  (C) Copr. 1986-92 Numerical Recipes Software &1%+%+5)+.

      !******************************************************************
      subroutine splint(xa, ya, y2a, n, x, y, ierr)
        include '../INCLUDES/precision.inc'
        integer, intent(in) :: n
        real(double), intent(in) :: x, xa(n), y2a(n), ya(n)
        real(double), intent(out) :: y
        integer, intent(out) :: ierr
        !*** local variables
        integer k, khi, klo
        real(double) a, b, h

        ierr = 0
        klo=1
        khi=n
1       if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if(xa(k).gt.x)then
              khi=k
           else
              klo=k
           endif
           goto 1
        endif
        h=xa(khi)-xa(klo)
        if (h.eq.0.) then
!           	     print*,'bad xa input in splint'
           ierr = 1
           return
        endif
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        return

      end subroutine splint
      !  (C) Copr. 1986-92 Numerical Recipes Software &1%+%+5)+.
