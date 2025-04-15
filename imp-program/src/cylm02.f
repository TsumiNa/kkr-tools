c                                                                       
cdeck
C     subroutine cylm02(lmax,np,cosx,fai)
      subroutine cylm02(lmax,np,cosx,fai)
c.....------------------------------------------------------------------
c     preparation of cylm0(=ylm(ip,i)), cylmt1(=dylm/dtheta),
c     cylmt2(=d2ylm/dt2), 
c     cylmf1, cylmf2 are for fai.
c     cylmtf=d2ylm/dfdt
c     i=1,2,....,(lmax+1)**2                  
c.....------------------------------------------------------------------
      implicit real*8 (a-h,o-z)   
  	include 'inc.fi'                                      
c       INTEGER LPOTD
c       PARAMETER(LPOTD=8)
      INTEGER N,IJD,LMMAXD
C     PARAMETER (N=2* (LPOTD+1),IJD=2*N**2,LMMAXD=(LPOTD+1)**2)
      PARAMETER (N=2* (LPOTD+1),IJD=434,LMMAXD=(LPOTD+1)**2)



c.....------------------------------------------------------------------
      complex cylm0(LMMAXD)
      complex cylmt1(LMMAXD),cylmt2(LMMAXD)
      complex cylmf1(LMMAXD),cylmf2(LMMAXD),cylmtf(LMMAXD)
      complex ep1f,ep2f,em1f,em2f,ci
      real*8 bb1(lmmaxd)
      integer llmax      
c.....------------------------------------------------------------------
      common/cylmdnp/thet(IJD),ylm(IJD,LMMAXD),dylmt1(IJD,LMMAXD),
     &  dylmt2(IJD,LMMAXD),dylmf1(IJD,LMMAXD),dylmf2(IJD,LMMAXD),
     &  dylmtf(IJD,LMMAXD)
c.....------------------------------------------------------------------
      integer np
      real*8 cosx(*),fai(*)
c.....
C     common/der2/np,cosx,fai                                
      real*8 yl(21)                                                  
c.....------------------------------------------------------------------
c                                                                       
      ci=cmplx(0.d0,1.d0)
      one=1
      pi=4*atan(one)                                                
      llmax=(lmax+1)**2

        nph=np/2

      write(6,*) ' in cylm02 np=',np        
      do 11 ip=1,np                                                     

c
c
c     if(ip.gt.nph) cosx(ip)=-cosx(ip)
c    
       thet(ip)=acos(cosx(ip))
       fi=fai(ip)
       di=2*fai(ip)
       ep1f=cmplx(cos(fi),sin(fi))
       em1f=conjg(ep1f)
       ep2f=cmplx(cos(di),sin(di))
       em2f=conjg(ep2f)

       do 21 l=0,lmax                                                  
c
          call spher(yl,l,cosx(ip))                                     
        do 20 m=-l,l                                                    
          mm=l+m+1                                                      
          i=(l+1)**2-l+m                                                
          aaa=m*fai(ip)
          ccc=cos(aaa)
          sss=sin(aaa)
          cylm0(i) = yl(mm) * dcmplx(ccc,sss)
c         write(6,9010) i,cylm0(i)
c9010  format(1x,' i cylm0',i5,2f10.5)
   20   continue 

         do 22 m=-l,l                                                   
          i=(l+1)**2-l+m                                                
          cylmt1(i)=0.
          cylmt2(i)=0.
          cylmtf(i)=0.
   22    continue 

        do 23 m=-l,l                                                    
          i=(l+1)**2-l+m                                                

          lmm1m=l-m-1
          lmm=l-m
          lmm1=l-m+1
          lmm2=l-m+2
          lm1m=l+m-1
          lm=l+m
          lm1=l+m+1
          lm2=l+m+2

          cylmt2(i)=cylmt2(i)-(lmm*lm1+lmm1*lm)/4.*cylm0(i)

          if(m+2.le.l) cylmt2(i)=cylmt2(i)+
     &      sqrt(float(lmm1m*lmm*lm1*lm2))/4*cylm0(i+2)*em2f

          if(m+1.le.l) cylmt1(i)=cylmt1(i)+
     &        sqrt(float(lmm*lm1))/2*cylm0(i+1)*em1f

          if(m-1.ge.-l) cylmt1(i)=cylmt1(i)-
     &        sqrt(float(lm*lmm1))/2*cylm0(i-1)*ep1f

          if(m-2.ge.-l) cylmt2(i)=cylmt2(i)+
     &      sqrt(float(lmm1*lmm2*lm1m*lm))/4*cylm0(i-2)*ep2f

   23   continue 

         do 24 m=-l,l                                                   
          i=(l+1)**2-l+m                                                
          cylmf1(i)=ci*m*cylm0(i)
          cylmf2(i)=-m*m*cylm0(i)
          cylmtf(i)=ci*m*cylmt1(i)
   24    continue 

   21    continue 
c  
c        calculate real spherical harmonics differenciated
c
c
c        write(6,9005) (cylm0(i),i=1,5)
 9005 format(1x,' cylm0',4f10.5)
         call trarea(cylm0,bb1,lmax) 

         DO 31 M=1,llmax
   31    ylm(ip,m)=bb1(m)
c
c        if(ip.eq.5) write(6,9006) (ylm(ip,i),i=1,5)
c9006 format(1x,' ylm',10f10.5)
c
c
         call trarea(cylmt1,bb1,lmax)
         DO 32 M=1,llmax
   32    dylmt1(ip,m)=bb1(m)
C
         call trarea(cylmt2,bb1,lmax)
         DO 33 M=1,llmax
   33    dylmt2(ip,m)=bb1(m)
c
         call trarea(cylmf1,bb1,lmax)
         DO 34 M=1,llmax
   34    dylmf1(ip,m)=bb1(m)
C
         call trarea(cylmf2,bb1,lmax)
         DO 35 m=1,llmax
   35    dylmf2(ip,m)=bb1(m)
C
         call trarea(cylmtf,bb1,lmax)
         DO 36 m=1,llmax
   36    dylmtf(ip,m)=bb1(m)
C
         if(ip.le.5.or.ip.eq.np) then
C
         write(6,9008) thet(ip),cosx(ip),fai(ip)
 9008    format(1x, ' thta cosx fai',6f12.8)
         write(6,9007) thet(ip),cosx(ip),fai(ip),
     &         (ylm(ip,m),dylmt1(ip,m),dylmt2(ip,m),
     &         dylmf1(ip,m),dylmf2(ip,m),dylmtf(ip,m),m=1,5)
 9007    format(1x,' thet cosx fai ylm  ',12f8.3)
C        
        end if
C
   11 continue                                                          
C        stop  ' cylm02'
         return                                                         
         end                                                        
