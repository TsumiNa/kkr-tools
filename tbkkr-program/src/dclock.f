C ************************************************************************
      DOUBLE PRECISION FUNCTION DCLOCK()
      implicit none
C ************************************************************************
C     .. External Functions ..
c     INTEGER MCLOCK
c     EXTERNAL MCLOCK
C     ..
c     DCLOCK = MCLOCK()/100.
c
c     Next is for Cray UNICOS
c      Real*8 tnow
c      External second
c      call second(tnow)
c      dclock = tnow
c
c     Next is for DEC Alphas
      Real etime, tarry(2)
      ! External etime
      dclock = Dble (ETIME(tarry))
c
c     Next is for no timing
c      DCLOCK = 0.
c
      RETURN
      END
