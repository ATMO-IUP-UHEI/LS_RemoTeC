
      subroutine LINTERP (X1, Y1, N1, X, Y, N)
      implicit none
      include '../INCLUDES/precision.inc'
      integer, intent(in) ::  n1, n 
      real(double), intent(in) :: x1(n1), y1(n1), x(n)
      real(double), intent(out) :: y(n)
!*** local variables 
      integer :: i, j, j0

      j0 = 2
      do i = 1, n
         do j = j0, n1
            if (x(i).lt.x1(j)) then
               y(i) = y1(j-1) + (x(i)-x1(j-1))/(x1(j)-x1(j-1))* &
                      (y1(j)-y1(j-1))
               j0 = j
               goto 10
            endif
         enddo
         j = n1
         y(i) = y1(j-1) + (x(i)-x1(j-1))/(x1(j)-x1(j-1))* &
                (y1(j)-y1(j-1))
         j0 = n1 + 1
 10      continue
      enddo
      return
      end

