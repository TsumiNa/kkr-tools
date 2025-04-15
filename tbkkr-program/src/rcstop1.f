c ************************************************************************
      SUBROUTINE RCSTOP1(C)
      implicit none
c ************************************************************************
C     .. Scalar Arguments ..
      CHARACTER*8 C
C     ..
C     .. Local Scalars ..
      INTEGER I
      DOUBLE COMPLEX CC(1)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSEND,FORFLUSH
C     ..
      I = 0
      CC(1) = DCMPLX(0.d0,0.d0)
      PRINT *,C
      CALL FORFLUSH(CC)
      CALL CSEND(1,CC,4,-1,0)
      STOP
      END
