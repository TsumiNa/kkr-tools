C ************************************************************************
      SUBROUTINE VMUL(a,b,c)
      implicit none
C ************************************************************************
      double precision a(*),b,c(*)
      integer i
      do 1 i=1,3
        c(i)=b*a(i)
 1    continue
      return
      END
