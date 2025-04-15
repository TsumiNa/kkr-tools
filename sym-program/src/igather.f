      SUBROUTINE IGATHER(N,A,B,INDEX)
      IMPLICIT NONE
C     REAL*8 A(N),B(*)
C     REAL*4 A(N),B(*)   FUER REAL*4, INTEGER*4, LOGICAL
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER A(N),B(*),INDEX(N)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        A(I) = B(INDEX(I))
   10 CONTINUE
      END
