      SUBROUTINE WFREAD(PNS,QNS,IOWFCT)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NSWR
      PARAMETER (NSWR=2)
      INTEGER IRNSD,LMAXD,LMX
      PARAMETER (irnsd=508,lmaxd=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER IRLLD
      PARAMETER (IRLLD= (IRNSD+1)*LMMAXD*LMMAXD*NSWR)
C     ..
C     .. Scalar Arguments ..
      INTEGER IOWFCT
C     ..
C     .. Array Arguments ..
c
c---> store regular wavefunctions and matrices on the buffer memory
c
      COMPLEX*16 PNS(IRLLD),QNS(IRLLD)
C     ..
      READ (IOWFCT) PNS,QNS
      RETURN

      END
