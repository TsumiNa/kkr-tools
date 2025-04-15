      SUBROUTINE TMWRIT(TMATLL,ITMAT)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER LMAXD,LMSQ
      PARAMETER (lmaxd=4,LMSQ= (LMAXD+1)**4)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT
C     ..
C     .. Array Arguments ..
      COMPLEX*16 TMATLL(LMSQ)
C     ..
      WRITE (ITMAT) TMATLL
      RETURN

      END
