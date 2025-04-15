      SUBROUTINE CINIT(N,A)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     initialize the first n values of a complex array a with zero
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      COMPLEX*16 A(*)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        A(I) = 0.0D0
   10 CONTINUE
      END
