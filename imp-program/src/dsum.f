      REAL*8 FUNCTION DSUM(N,V,IV)
      IMPLICIT NONE
CODER REAL*4  V(*),DSUM
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      REAL*8 V(*)
C     ..
C     .. Local Scalars ..
      REAL*8 VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = 0D0
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      DSUM = VSUM
      END
