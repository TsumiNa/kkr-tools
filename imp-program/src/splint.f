      Subroutine splint(xa,ya,y2a,n,x,y)
      Implicit None
c----------------------------------------------------------------------
c     This routine returns a cubic-spline interpolated value y at the
c     point x.
c
c     Input: arrays x and y of length n contain a tabulated function
c            y_i = f(x_i), x_1 < x_2 < ... < x_n.
c            array y2a is the output from spline routine
c            x is the point, where the interpolated value is wanted
c     Output: y is the interpolated value, i.e., y = f(x)
c
c     See: Press et al., Numerical Recipes (in Fortran),
c          (Cambridge University Press, 1989), pp. 88-89.
c
c                                     T. Korhonen, August 1997
c----------------------------------------------------------------------
c
c     Scalar arguments:
      Integer n
      Double Precision x, y
c
c     Array arguments:
      Double Precision xa(n), ya(n), y2a(n)
c
c     Local scalars:
      Integer klo, khi, k
      Double Precision h, a, b
c
      klo = 1
      khi = n
 1    Continue
      If ( (khi-klo) .gt. 1 ) Then
        k = (khi+klo)/2
        If ( xa(k) .gt. x ) Then
          khi = k
        Else
          klo = k
        End If
        Goto 1
      End If
c
      h = xa(khi) - xa(klo)
      If ( h .eq. 0.d0 ) Then
        Write(6,*) 'SPLINT: Bad xa input'
        Stop
      End If
c
      a = ( xa(khi) - x ) / h
      b = ( x - xa(klo) ) / h
      y = a*ya(klo) + b*ya(khi) +
     $     ( (a**3 - a)*y2a(klo) + (b**3 - b)*y2a(khi) )*(h**2) / 6.d0
c
      Return
      End
