c ************************************************************************
      SUBROUTINE TMWRIT(TMATLL,ITMAT)
      implicit none
c ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      INTEGER LMSQ
      PARAMETER (LMSQ= (LMAXD+1)**4)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX TMATLL(LMSQ)
C     ..
      WRITE (ITMAT) TMATLL
      RETURN

      END
