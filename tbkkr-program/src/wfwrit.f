c ************************************************************************
      SUBROUTINE WFWRIT(PNS,QNS,ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ,IOWFCT,INS)
      implicit none
c ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER IRMD,IRNSD,LMAXD,LMX
c      PARAMETER (IRMD=1484,IRNSD=508,LMAXD=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=(LMAXD+1)**2)
      INTEGER IRLLD
      PARAMETER (IRLLD= (IRNSD+1)*LMMAXD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DET
      INTEGER INS,IOWFCT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               FZ(IRMD,0:LMAXD),PNS(IRLLD),PZ(IRMD,0:LMAXD),
     +               QNS(IRLLD),QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD)
C     ..
c
c---> store wavefunctions and matrices
c
      IF (INS.NE.0) WRITE (IOWFCT) PNS,QNS
      WRITE (IOWFCT) ALPHA,DET,AR,CR,PZ,QZ,FZ,SZ
      RETURN
      END                           ! WFWRIT
