C ************************************************************************
      LOGICAL FUNCTION VTEST(a,b)
      implicit none
C ************************************************************************
      double precision a(*),b(*)
      INTRINSIC ABS
      VTEST=.FALSE.
      IF (ABS(A(1)-B(1))+ABS(A(2)-B(2))+ABS(A(3)-B(3)).LE.1.0D-4)
     +     VTEST=.TRUE.
      RETURN
      END
