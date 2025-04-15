      subroutine surfgf(ml,m0,mr,x,kdim,ndim,itermax,errmax,ichck)
c============================================================
c
c solve surface green's function: f(x)=ml*(m0-x)**(-1)*mr
c method: decimation technique
c NEW VERSION (speeded up) by V.Bellini (march,1999)
c
c input:  ml,m0,mr - complex rectangular matrices
c         kdim     - actual dimension of matrices
c         ndim     - declared dimension of matrices
c output: x        - result, matrix of same type as before
c=============================================================
c
c
       implicit none

      include 'inc.fi'

c
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ=(LMAXD+1)**2)
      INTEGER NNDIM
      PARAMETER (NNDIM=LMAXSQ*NPRINCD)
c
      double complex cone,czero
      double complex detl

      double complex ml(nndim,nndim),m0(nndim,nndim),mr(nndim,nndim),
     &               x(nndim,nndim)
      double complex alfa(nndim,nndim),beta(nndim,nndim),
     &               eps(nndim,nndim),y1(nndim,nndim),y2(nndim,nndim)
      double complex tempin(nndim,nndim),tempout(nndim,nndim)
      double complex cunit(nndim,nndim),aa(nndim,nndim),bb(nndim,nndim)
      double complex cc(nndim,nndim)
      double precision sum,err,xre,xim,errmax
      integer i,j,kdim,iter,itermax,ichck,ndim
      integer info,ipvt(nndim),n
      data cone/(1.d0,0.d0)/
      data czero/(0.d0,0.d0)/
      external cinit,zcopy,zgetrf,zgetrs,zaxpy
c
      if(ndim.gt.nndim) stop 'increase nndim!!'

      call cinit(kdim*kdim,cunit)
      do n = 1,kdim
        cunit(n,n) = cone
      enddo

      call zcopy(kdim*kdim,m0(1,1),1,eps(1,1),1)
      call zcopy(kdim*kdim,ml(1,1),1,alfa(1,1),1)
      call zcopy(kdim*kdim,mr(1,1),1,beta(1,1),1)
      call zcopy(kdim*kdim,m0(1,1),1,x(1,1),1)

      iter=1
    1 continue

      call zcopy(kdim*kdim,eps(1,1),1,y1(1,1),1)
      call zcopy(kdim*kdim,y1(1,1),1,tempin(1,1),1)
      call zgetrf(kdim,kdim,tempin(1,1),kdim,ipvt,info)

c     aa = eps^-1 * alfa
      call zcopy(kdim*kdim,alfa(1,1),1,tempout(1,1),1)
      call zgetrs('N',kdim,kdim,tempin(1,1),kdim,ipvt,
     +     tempout(1,1),kdim,info)
      call zcopy(kdim*kdim,tempout(1,1),1,aa(1,1),1)

c     bb = eps^-1 * beta

      call zcopy(kdim*kdim,beta(1,1),1,tempout(1,1),1)
      call zgetrs('N',kdim,kdim,tempin(1,1),kdim,ipvt,
     +     tempout(1,1),kdim,info)
      call zcopy(kdim*kdim,tempout(1,1),1,bb(1,1),1)

c     alfa_new = alfa * aa

      call zgemm('N','N',kdim,kdim,kdim,cone,alfa(1,1),kdim,
     +     aa(1,1),kdim,czero,y1(1,1),kdim)

c     beta_new = beta * bb

      call zgemm('N','N',kdim,kdim,kdim,cone,beta(1,1),kdim,
     +     bb(1,1),kdim,czero,y2(1,1),kdim)

c     cc = - alfa * bb

      call zgemm('N','N',kdim,kdim,kdim,-cone,alfa(1,1),kdim,
     +     bb(1,1),kdim,czero,cc(1,1),kdim)

c     x_new = x + cc

      call zaxpy(kdim*kdim,cone,cc(1,1),1,x(1,1),1)

c     cc = eps + cc

      call zaxpy(kdim*kdim,cone,cc(1,1),1,eps(1,1),1)

c     eps_new = cc - beta * aa

      call zgemm('N','N',kdim,kdim,kdim,-cone,beta(1,1),kdim,
     +     aa(1,1),kdim,cone,eps(1,1),kdim)

      call zcopy(kdim*kdim,y1(1,1),1,alfa(1,1),1)
      call zcopy(kdim*kdim,y2(1,1),1,beta(1,1),1)

      sum=0.d0
      do i=1,kdim
      do j=1,kdim
        xre=dreal(alfa(i,j))
        xim=dimag(alfa(i,j))
        sum=sum+xre*xre+xim*xim
      end do
      end do

      err=dsqrt(sum)
      if(err.lt.errmax.or.iter.gt.itermax) goto 2
      iter=iter+1
      goto 1

    2 continue

      call zcopy(kdim*kdim,x(1,1),1,tempin(1,1),1)
      call zcopy(kdim*kdim,cunit(1,1),1,tempout(1,1),1)
      call zgetrf(kdim,kdim,tempin(1,1),kdim,ipvt,info)
      call zgetrs('N',kdim,kdim,tempin(1,1),kdim,ipvt,
     +     tempout(1,1),kdim,info)
      call zcopy(kdim*kdim,tempout(1,1),1,x(1,1),1)

      call tripmt(ml,x,mr,kdim,kdim,ndim)



      if(iter.gt.itermax) then
        write(6,'('' itermax too small !  iter='',i3)') iter
      end if
      if(ichck.eq.0) return
c       write(6,'('' Surfgf:  iter='',i4,''  error='',d12.7)') iter,err
c      write(6,'(/'' X matrix'')')
c      call outmat(x,kdim,kdim,ndim,6)
c      write(6,*)
c
      return
      end
