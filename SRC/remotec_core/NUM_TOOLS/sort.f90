   subroutine sort(x_in, nx, index_out)
   implicit none
   include '../INCLUDES/precision.inc'
   integer, intent(in) :: nx
   real(double), dimension(nx), intent(inout) :: x_in
   integer, dimension(nx), intent(out) :: index_out
!*** local variables
   real(double), dimension(nx) :: x_dum
   real(double) :: first_value, min_value
   integer :: l, iswap

      x_dum = x_in
      do l=1,nx
         first_value = x_dum(l)
         iswap = l-1+minval(minloc(x_dum(l:nx)))
         min_value = minval(x_dum(l:nx)) 
         x_dum(l) = min_value
         index_out(l) = minval(minloc(DABS(x_in-min_value)))
         x_dum(iswap) = first_value
      enddo
      x_in = x_dum

   end subroutine sort
