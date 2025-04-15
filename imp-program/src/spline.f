      Subroutine spline(x,y,n,yp1,ypn,y2)
      Implicit None
c----------------------------------------------------------------------
c     This subroutine is used to calculate the first derivatives
c     needed in the spline interpolation routine splint. The
c     calculated 1st derivatives are stored in the array y2.
c     See: Press et al., Numerical Recipes (in Fortran),
c          (Cambridge University Press, 1989), pp. 88-89.
c
c     Input: arrays x and y of length n contain a tabulated function
c            y_i = f(x_i), x_1 < x_2 < ... < x_n.
c            yp1, ypn 1st derivatives at the end points, if the
c            value(s) is(are) greater than 1e30, a natural boundary
c            condition is used (i.e., the corresponding 2nd derivative
c            is put to zero).
c     Output: array y2 of length n containing the 2nd derivative of
c             the interpolating function  at points x_i.
c
c                                     T. Korhonen, August 1997
c----------------------------------------------------------------------
c
c     Parameters:
      Integer irid,irikd
      Parameter (irid=435,irikd=435)
      Integer irimad
      Parameter ( irimad = max0(irid,irikd) )
c
c     Scalar arguments:
      Integer n
      Double Precision yp1, ypn
c
c     Array arguments:
      Double Precision x(n), y(n), y2(n)
c
c     Local arrays:
      Double Precision u(irimad)
c
c     Local scalars:
      Integer i, k
      Double Precision sig, p, qn, un
c
      If ( irimad .lt. n ) Stop 'SPLINE: increase irid'
c
      If ( yp1 .gt. 0.99e30 ) Then
c     The lower boundary condition is set to be 'natural'
        y2(1) = 0.d0
        u(1)  = 0.d0
      Else
c     The lower boundary condition is set to have a specified
c     first derivative (yp1)
        y2(1) = -0.5d0
        u(1)  = ( 3.d0/(x(2)-x(1)) )*( (y(2)-y(1))/(x(2)-x(1)) - yp1 )
      End If
c
      Do i = 2, n-1
        sig   = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p     = sig*y2(i-1) + 2.d0
        y2(i) = (sig - 1.d0)/p
        u(i)  = ( 6.d0*( (y(i+1)-y(i))/(x(i+1)-x(i)) - (y(i)-y(i-1))
     $       /(x(i)-x(i-1)) ) / (x(i+1)-x(i-1)) - sig*u(i-1) ) / p
      End Do
c
      If ( ypn .gt. 0.99e30 ) Then
c     The upper boundary condition is set to be 'natural'
        qn = 0.d0
        un = 0.d0
      Else
c     The upper boundary condition is set to have a specified
c     first derivative (ypn)
        qn = 0.5d0
        un = ( 3.d0 / (x(n)-x(n-1)) ) * ( ypn - (y(n)-y(n-1)) /
     $       (x(n)-x(n-1)) )
      End If
c
      y2(n) = ( un - qn*u(n-1)) / ( qn*y2(n-1) + 1.d0 )
c
      Do k = n-1, 1, -1
        y2(k) = y2(k)*y2(k+1) + u(k)
      End Do
c
      Return
      End
