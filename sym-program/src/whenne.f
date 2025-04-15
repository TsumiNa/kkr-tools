C-----------------------------------------------------------------------
      SUBROUTINE WHENNE(N,V,IV,T,IN,L)
      IMPLICIT NONE
C     REAL*8  V(*),T
C     INTEGER V(*),T  FUER REAL*4,INTEGER UND LOGICAL
C     INTEGER*2 V(*),T
C     .. Scalar Arguments ..
      INTEGER IV,L,N,T
C     ..
C     .. Array Arguments ..
      INTEGER IN(N),V(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,J
C     ..
      J = 1
      L = 0
      IF (IV.LT.0) J = 1 + (N-1)* (-IV)
      DO 10 I = 1,N
        IF (V(J).NE.T) THEN
          L = L + 1
          IN(L) = I
        END IF
        J = J + IV
   10 CONTINUE
      END
