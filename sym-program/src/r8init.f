      SUBROUTINE R8INIT(R,N,RR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      REAL*8 RR
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL*8 R(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        R(I) = RR
   10 CONTINUE
      RETURN
      END
