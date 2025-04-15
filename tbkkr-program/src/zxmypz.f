c ************************************************************************
      SUBROUTINE ZXMYPZ(N,X,Y,Z)
      implicit none
c ************************************************************************
C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX X(*),Y(*),Z(*) 
C     ..
C     .. Local Scalars ..
      INTEGER I
C     ..
      DO 10 I = 1,N
        Z(I) = Z(I) + X(I)*Y(I)
 10   CONTINUE
      RETURN
c
      END
