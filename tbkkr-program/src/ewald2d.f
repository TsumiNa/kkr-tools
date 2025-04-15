      subroutine ewald2d(LPOT,ALAT,VEC1,VEC2,RM2,NRMAX,NSHLR,NSR,GN2,
     &                   NGMAX,NSHLG,NSG,SUM2D,SUMDIR,VOL)
c 10.01.2000 
      implicit none
c-----------------------------------------------------------------------
c
c     calculation of lattice sums for l .le. 2*lpot :
c
c                      ylm( q(i) - q(j) + rm )
c           sum      ===========================
c                    | q(i) - q(j) + rm |**(l+1)
c
c            - summed over all 2D lattice vectors rm  -
c
c     ylm       : real spherical harmic to given l,m
c     
c     The sum is done different in the plane (qi-qj)z = 0
c     and out of the plane. In plane an Ewald procedure similar
c     to the 3d is used and we perform 2 sums (real and reciprocal)
c     the l= 2,4 m=0 terms are calculated with a different method
c    
c
c     The l=0 term is calculated with a extra factror sqrt(4*pi) this 
c     is for transparency reasons (so that the correction terms *r=0,g=0* 
c     can be followed in the program) 
c     Literature : lm = (0,0), (1,0) terms :PRB 40, 12164 (1989)
c                                           PRB 49, 2721 (1994)  
c                                           PRB 47, 16525 (1993)
c                                           Zimman , p.39-40
c         l=2,4 (m=0) terms are done with recursive diferentiation
c                                                    v. 16.8.99
c         The l multipols are treated using the expansion
c         for a complex plane wave. 
c         eq.(33) , M. Weinert, J. Math Phys. 22, 2439 (1981)
c                                                    
c     Final version : 11.01.2000   (No direct sum neaded everything
c                                   is done with Ewald method)                    
c----------------------------------------------------------------------
c     .. parameters ..
      include 'inc.fi'
      integer l2maxd,l2mmxd
      parameter (l2maxd=2*LPOTD,l2mmxd= (l2maxd+1)**2)
c     ..
c     .. scalar arguments ..
      double precision alat,vol
      integer lpot,natom,ngmax,nrmax,nshlg,nshlr
c     ..
c     .. array arguments ..
      double precision  gn2(2,*),rm2(2,*),vec1(3),vec2(3),
     &                  sum2d(l2mmxd),SUMDIR(l2mmxd)
      integer nsg(*),nsr(*)
c     ..
c     .. local scalars ..
      double complex bfac,cfac,ci,bfac1,bf,simag
      double precision alpha,beta,bound,dq1,dq2,dq3,dqdotg,expbsq,fpi,
     +                 g1,g2,g3,ga,dqr,factor,ga0,dot1,con,ga1,ga2,
     +                 lamda,pi,r,r1,r2,r3,rfac,s,sgm,fun,stest0,r0
      double precision rp,tpi
      integer i,i1,i2,it,l,lm,lmax,lmmax,m,nge,ngs,nre,nrs,nstart
      integer im,ir,is
      logical TEST
c     ..
c     .. local arrays ..
      double complex stest(l2mmxd),s0(l2mmxd)
      double complex stestnew(l2mmxd)
      double complex apref,apref1,expon,aprefpp,aprefmm
      double precision gr(0:4),ylm(l2mmxd),gi(0:4),
     &                 pref0(0:l2maxd),g(0:l2maxd)
      double precision dfac(0:2*l2maxd+1),facl(0:l2maxd+1)
      double precision erfcex,signrz,PL0(0:l2maxd) ! Legendre Plm(0)
c     ..
c     .. external subroutines ..
      external fplaner,fplaneg
c     ..
c     .. intrinsic functions ..
      intrinsic abs,aimag,atan,exp,real,sqrt
c     ..
c     .. data statements ..
      data ci/ (0.0d0,1.0d0)/,bound/1.0d-9/
c     ..
      pi    = 4.0d0*datan(1.0d0)
      fpi   = 4.0d0*pi
      tpi   = 2.0d0*pi
c Factorial 
      dfac(0) = 1
      do l=1,2*l2maxd+1
       dfac(l) = dfac(l-1)*l
      end do
      do l=0,l2maxd
       pref0(l) = 0.d0
      end do
      LMAX = 2*LPOT
      LMMAX = (lmax+1)**2
c
      pref0(2) = sqrt(5.d0/pi)/2.d0/2.d0
      pref0(4) = 3.d0*sqrt(9.d0/pi)/16.d0/9.d0
c      pref0(6) = sqrt(13.d0/pi)/32.d0/13.d0
c
c---> choose proper splitting parameter
c
      lamda = sqrt(pi)/alat
c
      dq1 = (VEC2(1) - VEC1(1))*ALAT ! SCALE WITH ALAT
      dq2 = (VEC2(2) - VEC1(2))*ALAT
      dq3 = (VEC2(3) - VEC1(3))*ALAT 
c Initialize
      do lm = 1,lmmax
         stest(lm) = 0.0d0         
         stestnew(lm) = 0.d0
      end do
c
c---> Add correction if rz = 0
c
      if (dabs(dq3).lt.1.d-6) then
         stest(1) = stest(1) - 2.d0*lamda/sqrt(pi) 
     &                       - 2.d0*sqrt(pi)/lamda/vol
         stestnew(1) = stestnew(1) - 2.d0*lamda/sqrt(pi) 
     &                       - 2.d0*sqrt(pi)/lamda/vol

         if ((dq1*dq1+dq2*dq2).gt.1.d-6) then
            stest(1) = stest(1) + 2.d0*lamda/sqrt(pi)
            stestnew(1) = stestnew(1) + 2.d0*lamda/sqrt(pi)
         end if
      else
c     
c---> Add correction if rz<> 0
c   
         stest(1) = stest(1) - dabs(dq3)*fpi/2.d0/vol
         stest(3) = stest(3) - dabs(dq3)/dq3*sqrt(3.d0*fpi)/2.d0/vol ! -d/dz
c the correction for higher l vanishes...
      end if

      if (dabs(dq3).lt.1.d-6) THEN 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                 I N          P L A N E      M = 0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc            
c
c---> Real space sum  
c
c         WRITE(6,*) 'IN PLANE *********'       

         do ir = 1,NRMAX        
            r1 = dq1 - rm2(1,ir) 
            r2 = dq2 - rm2(2,ir) 
            r3 = 0.d0
            r = sqrt(r1*r1+r2*r2)

            if (r.gt.1.d-8) then
               alpha = lamda*r 
               call fplaner(alpha,gr,l2maxd,r)
               do l=0,4
                  lm = l*(l+1) + 1 ! m =0 
                  stest(lm) = stest(lm) +  gr(l)
               end do
c                 
                  call ymy(r1,r2,r3,r0,ylm,l2maxd)
                  call gamfc(alpha,g,l2maxd,r)
                  ylm(1) = 1.d0  ! just definition matter

                  do l = 0,lmax
                     rfac = g(l)/sqrt(pi)
                     do m = -l,l
                        lm = l* (l+1) + m + 1
                        stestnew(lm) = stestnew(lm)+ylm(lm)*rfac
                     end do
                  end do 

               if (ir.eq.(NRMAX-NSR(NSHLR))) THEN 
c     keep the value before the last shell to test convergence
                  do l=0,lmax
                     do m = -l,l
                     lm = l*(l+1) + m + 1 
                     s0(lm) = stestnew(lm)
                     end do
                  end do 
               end if   
            end if  ! r <> 0

         end do                 ! ir loop    
c     
c---  > Check convergence 
c  
         s = 0.d0
         do l=0,lmax
            do m=-l,l
            lm = l*(l+1) + m + 1 
            stest0 = abs(s0(lm)-stestnew(lm)) 
            if (s.lt.stest0) s = stest0
            end do
         end do            
         IF (s.gt.bound) write (6,fmt=9000) abs(s)
         IF (TEST('electro ')) THEN
         do l=0,lmax
            do m = -l,l
               lm = l*(l+1) + m + 1 
               if (abs(stestnew(lm)).gt.1.d-9) then
                  write(6,2000) lm,stestnew(lm)
               end if
            end do
         end do            
 2000    format('real sum',I4,6D14.7)
        END IF
c     
c---  > Sum in reciprocal lattice
c     

         con = fpi/2.d0/vol

         do im = 1,ngmax         
            g1 = gn2(1,im)
            g2 = gn2(2,im)
            g3 = 0.d0
            ga = sqrt(g1*g1+g2*g2)
c -------------------------------------------------
            dot1 = dq1*g1+dq2*g2
            CALL FPLANEG(lamda,GI,pref0,l2maxd,ga,vol)            
            simag = exp(CI*dot1)
            do l=0,4
               lm = l* (l+1) + 1
               stest(lm) = stest(lm) + gi(l)*simag
            end do   
c --------------------------------------------------
c 
               if (ga.gt.1.d-6) then      
                  call ymy(g1,g2,g3,ga,ylm,l2maxd)
c     
                  beta = ga/lamda
                  expbsq = erfcex(beta/2.d0)
c     
                  bfac = con*simag*expbsq
                  stestnew(1) = stestnew(1) + bfac/ga
                  
                  do  l = 0,lmax    
                     IF (L.NE.0) THEN
                        do  m = -l,l
                           lm = l* (l+1) + m + 1
                           stestnew(lm) = stestnew(lm)+
     &                                         ylm(lm)*bfac*ga**(l-1)
                        end do
                     END IF
                     bfac = bfac/ci/real(2*l+1)
                  end do       
               end if

               
               if (im.eq.(NGMAX-NSG(NSHLG))) THEN 
c     keep the value before the last shell to test convergence
                  do lm=1,lmmax
                     s0(lm) = stestnew(lm)
                  end do 
               end if
               
            end do              ! imaginary loop
c     
c---  > test convergence
c     
         do lm=1,lmmax
            stest0 = abs(s0(lm)-stestnew(lm)) 
            if (s.lt.stest0) s = stest0
         end do
c
C Correction due to r=0 term only for DRn = 0
c
         if ((dq1*dq1+dq2*dq2).LT.1.D-6) THEN ! if statement Added on 28.11.1999   
            do l=2,4,2
               pref0(l) = pref0(l)*dfac(l)/dfac(l/2)/(l+1)
            end do
            
            i = 1
            do l=2,4,2
               i = i + 1
               lm = l*(l+1)+1 
               stest(lm) = stest(lm) 
     &              + (-1)**i*2.d0/sqrt(pi)*pref0(l)*lamda**(l+1)
            end do
         END IF         ! Added on 28.11.1999
c--------------------------------------------------------
 
         do l=2,4,2
            lm = l*(l+1)+1
            STESTNEW(lm) = STEST(lm)
         end do
c
c     end of correction
c
c--------------------------------------------------------       
         IF (TEST('electro ')) THEN
         do l=0,lmax
            do m=-l,l
            lm = l*(l+1) + m+ 1 
            if (abs(stestnew(lm)).gt.1.d-9) then
               write(6,3000) lm,stestnew(lm)
            end if
            end do
         end do 
 3000    format('Im sum',I4,6D14.7)
         END IF

c * * * * * * * * * * * * * * * * * * * * * * * * * * *
c     Fill up the array with the direct sum results
c        do l=1,2*LPOT             
c !           do m=-l,l
c !                 lm = l*(l+1)+m+1
c !                 if (lm.eq.7.or.lm.eq.21) then  ! these two are found differently
c !                    stest(lm) = stest(lm) 
c !                    else
c !                    stest(lm) = sumdir(lm)     
c !                 end if
c !           end do 
c !        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc         
         IF (s.gt.bound) write (6,fmt=9000) abs(s)
c
      do lm=1,lmmax
         stest(lm) = stestnew(lm)
      end do     
      else                      ! dabs(dq3) IN PLANE ENDED 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     O U T      O F      T H E     P L A N E
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c---  > Sum in reciprocal space   
c     
         do i = 2,ngmax  
c
c   Exclude the origin all terms vanish for g=0 
c   except the l = 0 componets which are treated 
c   separetly. (look at the begining of the sub)
c     
            signrz = dq3/dabs(dq3)
            g1 = gn2(1,i)
            g2 = gn2(2,i)                  
            ga = sqrt(g1*g1+g2*g2)
      
            expbsq = exp(-ga*dabs(dq3))
            dqdotg = (dq1*g1 + dq2*g2)
c
c     In case rz < 0 then multiply by (-1)**(l-m)
c     !
c     ! M. Weinert J. Math Phys. 22, 2439 (1981) formula 33
c     ! compare also formula A9 
c     
            do l= 0,lmax
c     m = 0
               apref1 = sqrt((2.d0*l+1.d0)/fpi)/dfac(l)
               if (l.eq.0) apref1 =1.d0 ! This is for consistency with amn 
               apref = apref1*tpi/vol*ga**l*(-signrz)**l 
      
               lm = l* (l+1) + 1
               stest(lm) = stest(lm) + apref*exp(CI*dqdotg)*expbsq/ga
c     m <> 0     
c
               expon = 1.d0
      
               do m = 1,l 
                  apref1 = 
     &            sqrt( (2.d0*l+1.d0)/2./fpi/dfac(l+m)/dfac(l-m) )
                  apref = apref1*tpi/vol*
     &                 ga**l*(-1)**l* signrz**(l-m) 
c      
c Go from the <usual> Jackson Ylm to the ones we use
c
                  expon = (g1+CI*g2)/ga*expon    ! exp(i*m*fi)
c    
                  aprefpp = ((-CI)**m/expon + expon/(CI)**m )  
                  aprefmm = ((-CI)**m/expon - expon/(CI)**m )*CI
c     m > 0   
                  lm = l* (l+1) + m + 1
                  stest(lm) = stest(lm) + 
     &                 apref*aprefpp*exp(CI*dqdotg)*expbsq/ga
c     m < 0                  
                  lm = l* (l+1) - m + 1
                  stest(lm) = stest(lm) + 
     &                 apref*aprefmm*exp(CI*dqdotg)*expbsq/ga
      
               end do
            end do
            
ccccccccccccccccccccccccccccccccccccccccccccccccccc
            if (i.eq.(NGMAX-NSG(NSHLG))) THEN 
c     keep the value before the last shell to test convergence
               do lm=1,lmmax
                  s0(lm) = stest(lm)
               end do 
            end if
c     
         end do                 ! i loop (rec lattice sum) 
c         
c---  > test convergence
c     
         s = 0.d0
         do lm=2,lmmax
            stest0 = abs(s0(lm)-stest(lm)) 
            if (s.lt.stest0) s = stest0
         end do  
         IF (s.gt.bound) write (6,fmt=9000) abs(s)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      end if                    ! Both cases finished
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      do lm = 1,LMMAX
         if (ABS(DIMAG(STEST(LM))).GT.BOUND) GOTO 120 
         sum2d(lm) = real(stest(lm))
         stest(lm) = 0.0d0
      end do
c     
      !write(6,*) '>>>> EWALD 2D'
c This is for consistency with the madelcoef sub
       sum2d(1) = sum2d(1)/sqrt(fpi)
c
      IF (TEST('electro ')) THEN
      do lm = 1,lmmax
         if (abs(sum2d(lm)).gt.1.0d-9) write (6,9010) 
     &        lm,sum2d(lm)
      end do
      END IF
c     
      return
      
 120  stop ' imaginary contribution to real 2d lattice sum '
      
 9000 format (1x,' convergence of sum2d : ',d10.5,
     +     ' is less than 1.0d-8 - use more lattice vectors ')
 9010 format('EWALD 2d ',I4,2D18.9)
      end
