c ************************************************************************
      DOUBLE PRECISION FUNCTION ZNORM2(N,A,INC)
      implicit none
c ************************************************************************
      DOUBLE COMPLEX A(*)
      INTEGER N,INC
c
      DOUBLE PRECISION NORM
      INTEGER I
      INTRINSIC DCONJG
c ------------------------------------------------------------------------
      NORM = 0.d0
      DO 10 I = 1,N,INC
        NORM = NORM + A(I)*DCONJG(A(I))
 10   END DO
c
      ZNORM2 = DSQRT(NORM)
c
      RETURN
      END
