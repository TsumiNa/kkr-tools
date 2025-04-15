C ************************************************************************
      SUBROUTINE ESAVE(E1,N1,IP)
      implicit none
c ************************************************************************
c     p.zahn, april 96
c ------------------------------------------------------------------------
      include 'inc.fi'
c     .. arguments
      INTEGER N1,IP
      DOUBLE PRECISION E1
c     .. locals
      DOUBLE PRECISION E,EE,EEN
      INTEGER OPT,N
c     ..arrays in common
      INTEGER     IUP,IDO,NUM
      DOUBLE PRECISION WUP,WDO
      COMMON /EIGENV/ WUP,WDO,NUM,IUP,IDO

      DOUBLE PRECISION EE1(0:M2D,KPOIBZ),EE2(0:M2D,KPOIBZ),
     +                 WW1(0:M2D,KPOIBZ),WW2(0:M2D,KPOIBZ)
      INTEGER II1(0:M2D,KPOIBZ),II2(0:M2D,KPOIBZ)
      COMMON /EE/ EE1,EE2,WW1,WW2,II1,II2 ! index arrays for eigenvalue
                                          ! analysis

      INTEGER NBASIS(KPOIBZ),NMIN(KPOIBZ),NMAX(KPOIBZ)
      COMMON /NBASIS / NBASIS,NMIN,NMAX

      DOUBLE PRECISION TH
      DATA  TH   / 5.0D2 /

      SAVE
c ------------------------------------------------------------------------
      N = N1 - NMIN(IP)
      IF (N.LT.0 .OR. N.GT.M2D) THEN
        write(6,*) 'IP,N1,NMIN(IP),M2D :',ip,n1,nmin(ip),m2d
        STOP 'ESAVE'
      END IF
      E = E1
      EE =  EE1(N,IP)
      EEN = -EE2(N+1,IP)
      IF (EEN.LT.-TH) EEN = TH
      IF ((((EE.LT.TH).AND.(E.GT.EE)) .OR.
     +     (EE.GT.TH)) .AND.
     +     (E.LE.EEN) ) THEN
        EE1(N,IP) = E
        WW1(N,IP) = WUP
        II1(N,IP) = IUP
      END IF

      E = -E1
      EE =  EE2(N,IP)
      EEN = -EE1(N-1,IP)
      IF (N.EQ.0 .OR. EEN.LT.-TH ) EEN = TH
      IF ((((EE.LT.TH).AND.(E.GT.EE)) .OR.
     +     (EE.GT.TH)) .AND.
     +     (E.LE.EEN) ) THEN
        EE2(N,IP) = E
        WW2(N,IP) = WDO
        II2(N,IP) = IDO
      END IF

c      write(6,FMT='(2i6,2f12.6)') N,IP,EE1(N,IP),EE2(N,IP)

      RETURN
      END 
