C-----------------------------------------------------------------------
      SUBROUTINE ISCATTER(N,A,INDEX,B)
      IMPLICIT NONE
C     REAL*8 A(*),B(N)
C     REAL*4 A(N),B(*)   FUER REAL*4, INTEGER*4, LOGICAL
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      INTEGER A(*),B(*),INDEX(N)
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        A(INDEX(I)) = B(I)
c	write(6,*) I,index(i),a(i),b(i)
   10 CONTINUE
      END
