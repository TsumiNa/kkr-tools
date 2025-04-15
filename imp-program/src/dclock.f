      REAL*8 FUNCTION DCLOCK()
      Implicit None
c
c     Next is for IBM RS6000 workstations
c      Integer mclock
c      External mclock
c      dclock = mclock()/100.
c
c     Next is for Cray UNICOS
c      Real*8 tnow
c      External second
c      call second(tnow)
c      dclock = tnow
c
c     Next is for DEC Alphas
      Real etime, tarry(2)
      External etime
      dclock = Dble (etime(tarry))
c
c     Next is for no timing
c      DCLOCK = 0.
c
      RETURN
      END
