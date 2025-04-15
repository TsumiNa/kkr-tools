
      Double Precision Function fun5(x,alpha)
      Implicit None
c
c     Calculates the value of the derivative of Fermi function
c
      Double Precision x, alpha
c     Lorenzial (not good!!)
c     fun5(x,alpha) = alpha*alpha/(x*x+alpha*alpha)

c     Fermi function, better localized (used !!)
      If ( Dabs(x/alpha) .lt. 80.d0 ) Then
c        fun5 = alpha*( Exp(-x/alpha) + 2.d0 + Exp(x/alpha) )
c        fun5 = 1.d0 / fun5
        fun5 = (1.d0/alpha)*exp(-x/alpha) / ( exp(-x/alpha) + 1 )**2
      Else
        fun5 = 0.d0
      End If
c
      End
