
      SUBROUTINE EXCH91(d,s,u,v,exl,exg,vxl,vxg)                          
c.....----------------------------------------------------------------- 
c     gga91 exchange for a spin-unpolarized electronic system           
c.....----------------------------------------------------------------- 
c     input d: density                                                  
c           s:  abs(grad d)/(2*kf*d)                                    
c           u:  (grad d)*grad(abs(grad d))/(d**2 * (2*kf)**3)           
c           v: (laplacian d)/(d*(2*kf)**2)                              
c     output:  exchange energy (ex) and potential (vx) in ry.           
                                                                        
c       kf=cbrt(3*pai**2*d).                                            
c.....----------------------------------------------------------------- 
      IMPLICIT REAL*8 (a-h,o-z)                                 
c.....----------------------------------------------------------------- 
      data a1,a2,a3,a4/0.19645d0,0.27430d0,0.15084d0,100.d0/            
      data ax,a,b1/-0.7385588d0,7.7956d0,0.004d0/                       
      data thrd,thrd4/0.333333333333d0,1.33333333333d0/                 
c.....----------------------------------------------------------------- 
      fac = ax*d**thrd                                                  
      s2 = s*s                                                          
      s3 = s2*s                                                         
      s4 = s2*s2                                                        
      p0 = 1.d0/SQRT(1.d0+a*a*s2)                                      
      p1 = LOG(a*s+1.d0/p0)                                            
      p2 = EXP(-a4*s2)                                                 
      p3 = 1.d0/(1.d0+a1*s*p1+b1*s4)                                    
      p4 = 1.d0+a1*s*p1+(a2-a3*p2)*s2                                   
      f = p3*p4                                                         
      ex = fac*f                                                        
c  local exchange exl                                                   
      exl = fac                                                         
      exg=ex-exl                                                        
                                                                        
c  energy done. now the potential:                                      
      p5 = b1*s2-(a2-a3*p2)                                             
      p6 = a1*s*(p1+a*s*p0)                                             
      p7 = 2.d0*(a2-a3*p2)+2.d0*a3*a4*s2*p2-4.d0*b1*s2*f                
      fs = p3*(p3*p5*p6+p7)                                             
      p8 = 2.d0*s*(b1-a3*a4*p2)                                         
      p9 = a1*p1+a*a1*s*p0*(3.d0-a*a*s2*p0*p0)                          
      p10 = 4.d0*a3*a4*s*p2*(2.d0-a4*s2)-8.d0*b1*s*f-4.d0*b1*s3*fs      
      p11 = -p3*p3*(a1*p1+a*a1*s*p0+4.d0*b1*s3)                         
      fss = p3*p3*(p5*p9+p6*p8)+2.d0*p3*p5*p6*p11+p3*p10+p7*p11         
      vx = fac*(thrd4*f-(u-thrd4*s3)*fss-v*fs)                          
c  local exchange vxl:                                                  
      vxl = fac*thrd4                                                   
      vxg=vx-vxl                                                        
      if(ipr.eq.0) then
        write(6,'(/'' p5,p6,p7,fs,p8,p9,p10,p11,fss,vxl,g='',12d11.5)')
     &  p5,p6,p7,fs,p8,p9,p10,p11,fss,vxl,vxg
        ipr=1
      endif
                                                                        
c in ry and energy density.                                             
      exl=exl*2.*d                                                      
      exg=exg*2.*d                                                      
      vxl=vxl*2.                                                        
      vxg=vxg*2.                                                        
                                                                        
      return                                                            
      end    
