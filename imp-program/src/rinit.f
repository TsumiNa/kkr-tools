      SUBROUTINE RINIT(N,A)
      IMPLICIT NONE
c-----------------------------------------------------------------------
C     INITIALIZE THE FIRST N VALUES OF A REAL ARRAY A WITH ZERO
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL*8 A(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        A(I) = 0.0D0
   10 CONTINUE
      END
