      function qpythag(a,b)
      include '../INCLUDES/precision.inc'
      real(quad) :: a,b,qpythag
      real(quad) ::  absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        qpythag=absa*sqrt(1.0q0+(absb/absa)**2)
      else
        if(absb.eq.0.0q0)then
          qpythag=0.0q0
        else
          qpythag=absb*sqrt(1.0q0+(absa/absb)**2)
        endif
      endif
      return
      end function qpythag
!  (C) Copr. 1986-92 Numerical Recipes Software Y_'i&.
