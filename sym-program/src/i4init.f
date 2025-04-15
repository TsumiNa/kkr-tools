      SUBROUTINE I4INIT(J,N,JR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER JR,N
C     ..
C     .. Array Arguments ..
      INTEGER J(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        J(I) = JR
   10 CONTINUE
      RETURN
      END
