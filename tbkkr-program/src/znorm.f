c ************************************************************************
      SUBROUTINE ZNORM(CC,MAX,VAL)
      implicit none
c ************************************************************************
c    
c     normalize the vector CC (with MAX double complex elements) 
c     to a length VAL
c
c     p.zahn, april 96
c ------------------------------------------------------------------------
c     .. arguments
      DOUBLE COMPLEX CC(*)
      DOUBLE PRECISION VAL
      INTEGER MAX
c     .. locals
      INTEGER I
      DOUBLE COMPLEX V
c     .. external
      INTRINSIC DSQRT,ABS
c ------------------------------------------------------------------------
      V = (0.0D0,0.0D0)
      DO 10 I=1,MAX
        V = V + ZABS(CC(I)*CC(I))
 10   END DO
      V = ZSQRT(V/VAL)*CC(1)/ABS(CC(1))
      DO 20 I=1,MAX
        CC(I) = CC(I) / V
 20   END DO
      RETURN
      END
