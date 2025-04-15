      SUBROUTINE C8INIT(R,N,RR)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      COMPLEX*16 RR
      INTEGER N
C     ..
C     .. Array Arguments ..
      COMPLEX*16 R(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        R(I) = RR
   10 CONTINUE
      RETURN
      END
