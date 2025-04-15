C ************************************************************************
      SUBROUTINE VEQ(a,b)
      implicit none
C ************************************************************************
      double precision a(*),b(*)
      integer i
      do 1 i=1,3
        b(i)=a(i)
 1    continue
      return
      END
