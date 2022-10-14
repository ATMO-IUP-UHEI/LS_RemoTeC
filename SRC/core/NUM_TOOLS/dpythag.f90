     function dpythag(a, b)
     include '../INCLUDES/precision.inc'
     real(double), intent(in) :: a, b
     real(double) :: dpythag
     real(double) :: absa, absb

      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        dpythag=absa*sqrt(1.0d0+(absb/absa)**2)
      else
        if(absb.eq.0.0d0)then
          dpythag=0.0d0
        else
          dpythag=absb*sqrt(1.0d0+(absa/absb)**2)
        endif
      endif
      return
      end function dpythag
!  (C) Copr. 1986-92 Numerical Recipes Software Y_'i&.
