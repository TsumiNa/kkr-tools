


      subroutine tripmt(u,b,ust,ndi1,ndi2,ndim)
c======================
c
c vectorized routine for triple product of rectangular matrices
c
c      implicit real*8(a-h,o-z)
      implicit none
      include 'inc.fi'

c
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ=(LMAXD+1)**2)
      INTEGER NNDIM
      PARAMETER (NNDIM=LMAXSQ*NPRINCD)
c
       integer ndim,ndi1,ndi2
      double complex cone,czero
      double complex u(ndim,ndim),ust(ndim,ndim),b(ndim,ndim)
      double complex c(nndim,nndim),temp(nndim,nndim)
c
      data cone/(1.d0,0.d0)/
      data czero/(0.d0,0.d0)/
c      external zgemul
      external zgemm



      call zgemm('N','N',ndim,ndim,ndim,cone,b(1,1),ndim,
     +     ust(1,1),ndim,czero,temp(1,1),ndim)
      call zgemm('N','N',ndim,ndim,ndim,cone,u(1,1),ndim,
     +     temp(1,1),ndim,czero,b(1,1),ndim)

      return
      end
