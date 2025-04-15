      SUBROUTINE RCSTOP(C)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      CHARACTER*8 C
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
C     .. External Subroutines ..
      EXTERNAL CSEND,FORFLUSH
C     ..
      I = 0
      PRINT *,C
      CALL FORFLUSH(6)
      CALL CSEND(1,I,4,-1,0)
      STOP
      END
