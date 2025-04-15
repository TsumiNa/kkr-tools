C ************************************************************************
      SUBROUTINE VSUB(a,b,c)
      implicit none
C ************************************************************************
      double precision a(*),b(*),c(*)
      integer i
      do 1 i=1,3
        c(i)=a(i)-b(i)
 1    continue
      return
      END
