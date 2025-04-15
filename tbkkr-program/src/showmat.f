c ************************************************************************
      subroutine showmat(str,a,DIM,I,J)
      implicit none
c ************************************************************************
c     output of matrix elements for test purpose
c
c ------------------------------------------------------------------------
      integer dim,I,J
      double complex a(DIM,*)
      character*4 str
c ------------------------------------------------------------------------
c
      write(6,FMT='(a4,1p,2i6,2d20.10)') str,I,J,a(I,J)
c
      return
      end
