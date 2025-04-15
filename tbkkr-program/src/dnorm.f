c ************************************************************************
      SUBROUTINE DNORM(CC,MAX,VAL)
      implicit none
c ************************************************************************
c    
c     normalize the vector CC (with MAX double precision elements) 
c     to a length VAL
c
c     p.zahn, april 96
c ------------------------------------------------------------------------
c     .. arguments
      DOUBLE PRECISION CC(*),VAL
      INTEGER MAX
c     .. locals
      INTEGER I
      DOUBLE PRECISION V
c     .. external
      INTRINSIC DSQRT
c ------------------------------------------------------------------------
      V = 0.0D0
      DO 10 I=1,MAX
        V = V + CC(I)*CC(I)
 10   END DO
      V = DSQRT(V)/VAL
      DO 20 I=1,MAX
        CC(I) = CC(I) / V
 20   END DO
      RETURN
      END
