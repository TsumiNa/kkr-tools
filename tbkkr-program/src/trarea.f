                                                           
cdeck
c subroutine trarea
c
c     from complex to real  (differenciated spherical harmonics)
c
      SUBROUTINE TRAREA(a,b,lmax)
      implicit REAL*8 (a-h,o-z) 
c     INTEGER LPOTD,LMMAXD,IJD,lmax
c      parameter(LPOTD=8,LMMAXD=(LPOTD+1)**2)
      include 'inc.fi'
      INTEGER LLMAXD,lmax
      parameter(LMMAXD=(LPOTD+1)**2)
      INTEGER NN
      parameter(NN=2*(LPOTD+1))
      COMPLEX a(LMMAXD),ci
      REAL*8 b(LMMAXD)
      data rtwo/1.414213562373d0/
      data pi,ci/3.14159265359d0,(0.d0,1.d0)/


C
C    calculate real the spherical harmonics derivetived
C
         i=0
         do 60 l=0,lmax
         i=i+l +1
         b(i)=REAL(a(i)) 
c        write(6,9000) a(i),b(i)
 9000    format(1x,' a=',4f10.5,' b=',f10.5)
         sgm=-1.e0
         do 70 m=1,l
         b(i-m)= REAL( ci*( a(i-m)-conjg(a(i-m)) ) )/rtwo
         b(i+m)= sgm*REAL((a(i+m)+conjg(a(i+m))))/rtwo
c        write(6,9000) a(i-m),conjg(a(i-m)),b(i-m)
         sgm=-sgm
   70 continue
      i=i+l
   60 continue
      return
      end
