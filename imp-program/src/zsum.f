      COMPLEX*16 FUNCTION ZSUM(N,V,IV)
      IMPLICIT NONE
CODER COMPLEX*8 V(*),ZSUM
C     .. Scalar Arguments ..
      INTEGER IV,N
C     ..
C     .. Array Arguments ..
      COMPLEX*16 V(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 VSUM
      INTEGER I,IBOT,ITOP
C     ..
      IF (IV.GE.0) THEN
        IBOT = 1
        ITOP = 1 + (N-1)*IV

      ELSE
        IBOT = 1 - (N-1)*IV
        ITOP = 1
      END IF

      VSUM = (0D0,0D0)
      DO 10 I = IBOT,ITOP,IV
        VSUM = VSUM + V(I)
   10 CONTINUE
      ZSUM = VSUM
      END
