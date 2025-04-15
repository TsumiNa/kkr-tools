
cdeck
c subroutine spher  ====*====3====*====4====*====5====*====6====*====7  
c                                                                       
c      spherical haminics exept the facter exp(i*m*phi)                 
c                                                                       
c      m=-l to l , for given l.                                         
c      x=cos(theta)                                                     
c                                                                       
c --*----1----*----2----*----3----*----4----*----5----*----6----*----7  
c                                                                       
      SUBROUTINE SPHER(ylm,l,x)                                         
c                                                                       
      IMPLICIT REAL*8 (a-h,o-z) 
      include 'inc.fi'
C     INTEGER LPOTD,LMMAXD 
c     PARAMETER (LPOTD=8,LMD0=2*LPOTD+1)                          
c                                                                       
      REAL*8 ylm(2*l+1)                                              
c     REAL*8 ylm(LMD0)                                              
c                                                                       
      pi = 4.0*DATAN(1.0d0) 
      IF (L.GT.LPOTD) THEN  
      WRITE(6,*) 'ERROR IN SPHER : ',L,LPOTD 
      STOP 
      END IF                                               
c                                                                       
c                                                                       
      ovr1=ABS(x)-1.d0                                                  
      if(ovr1.gt.0.1d-12) then                                          
        write(6,200) x                                                  
  200   format(//3x,'==invalid argument for spher; x=',e24.16,' ==')    
        stop                                                            
      else if(abs(ovr1).lt.1.d-10) then                                 
        if(x.gt.0.0) then                                               
          fac=1.0                                                       
        else                                                            
          fac=(-1)**l                                                   
        end if                                                          
        l2=2*l+1                                                        
        do 10 i = 1,l2                                                  
   10   ylm(i) = 0.0                                                    
        ylm(l+1) = SQRT(DFLOAT(l2)/(4.0*pi))*fac                        
        return                                                          
      end if                                                            
c                                                                       
c l<0                                                                   
      if(l.lt.0) then                                                   
        write(6,*) ' === l=',l,' < 0  : in sub.spher. ==='              
        stop '=== stop in sub.spher. (l<0) ==='                         
c l=0                                                                   
      else if(l.eq.0) then                                              
        ylm(1) = SQRT(1.0/(4.0*pi))                                     
c l=1                                                                   
      else if(l.eq.1) then                                              
        fac = SQRT(3.0/(4.0*pi))                                        
        ylm(1) = fac*SQRT((1.0-x*x)/2.0)                                
        ylm(2) = fac*x                                                  
        ylm(3) = -ylm(1)                                                
c l>1                                                                   
      else                                                              
        ylm(1) = 1.0                                                    
        ylm(2) = x                                                      
        do 20 i = 2,l                                                   
   20   ylm(i+1) = ((2*i-1)*x*ylm(i)-(i-1)*ylm(i-1))/i                  
        fac = 1.0d0/SQRT(1.0d0- x*x)                                    
        do 30 m = 1,l                                                   
          lm = l + m                                                    
          ylm(lm+1)=fac*(-(l-m+1)*x*ylm(lm)+(lm-1)*ylm(l))              
          if (m.lt.l) then                                              
            nn = m +1                                                   
            do 32 i = nn,l                                              
            ii = l-i + nn                                               
   32       ylm(ii)=fac*(-(ii-m)*x*ylm(ii)+(ii+m-2)*ylm(ii-1))          
          end if                                                        
   30   continue                                                        
        fac = SQRT((2*l+1)/(4.0d0*pi))                                  
        ylm(l+1) = fac * ylm(l+1)                                       
        do 40 m = 1,l                                                   
          fac = -fac/SQRT(DFLOAT((l+m)*(l-m+1)))                        
          lm = l + 1 + m                                                
          ln = l + 1 - m                                                
          qq = ylm(lm)                                                  
          ylm(lm) = fac * qq                                            
          ylm(ln) = ABS(fac) * qq                                       
   40   continue                                                        
      end if                                                            
c                                                                       
      return                                                            
      end    
