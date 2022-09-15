      subroutine dtwofft(data1, data2, fft1, fft2, n)
      implicit none
      include '../INCLUDES/precision.inc'
      integer, intent(in) :: n
      real(double), intent(in) :: data1(n),data2(n)
      complex(double), intent(out) :: fft1(n), fft2(n)
!     USES four1
      integer j, n2
      complex(double) :: h1, h2, c1, c2

      c1=cmplx(0.5,0.0)
      c2=cmplx(0.0,-0.5)
      do 11 j=1,n
        fft1(j)=cmplx(data1(j),data2(j))
11    continue
      call dfour1(fft1,n,1)
      fft2(1)=cmplx(dimag(fft1(1)),0.0)
      fft1(1)=cmplx(dble(fft1(1)),0.0)
      n2=n+2
      do 12 j=2,n/2+1
        h1=c1*(fft1(j)+conjg(fft1(n2-j)))
        h2=c2*(fft1(j)-conjg(fft1(n2-j)))
        fft1(j)=h1
        fft1(n2-j)=conjg(h1)
        fft2(j)=h2
        fft2(n2-j)=conjg(h2)
12    continue
      return
      end subroutine dtwofft
!  (C) Copr. 1986-92 Numerical Recipes Software Y_'i&.
