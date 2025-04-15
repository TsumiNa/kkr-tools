       SUBROUTINE MADELCOEF(I1,I2,LPOT,A,B,SMAT,CLEB,ICLEB,
     &                      IEND)
       implicit none
       include 'inc.fi'
       integer lmpotd,l3d,lm3d
       parameter (lmpotd= (lpotd+1)**2,l3d=2*lpotd,lm3d= (l3d+1)**2)
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
       integer ncleb1
       parameter (ncleb1=lm3d*lmpotd)
      
       INTEGER LPOT,IEND 
       double precision cleb(ncleb1),a(lmpotd,lmpotd),b(lmpotd),
     &                  smat(NAEZD,NAEZD,LM3D),dfac(0:lpotd,0:lpotd)
       integer icleb(ncleb1,3),loflm(lm3d)
       integer lmpot,lm1,lm2,i,lm3,l1,l2,l3,l,m,i1,i2
       double precision pi,fpi
c    
       pi = 4.0d0*datan(1.0d0)
       fpi = 4.0d0*pi
       lmpot = (lpot+1)**2
      i = 1
c
c---> determine the l-value for given lm
c
      do 10 l = 0,2*lpot
         do 20 m = -l,l
            loflm(i) = l
            i = i + 1
   20    continue
   10 continue      
c
c---> calculate:                             (2*(l+l')-1)!!
c                 dfac(l,l') = 4pi**2 *  ----------------------
c                                        (2*l+1)!! * (2*l'+1)!!
c
      dfac(0,0) = fpi*fpi
      do 30 l1 = 1,lpot
         dfac(l1,0) = dfac(l1-1,0)*real(2*l1-1)/real(2*l1+1)
         dfac(0,l1) = dfac(l1,0)
         do 40 l2 = 1,l1
            dfac(l1,l2) = dfac(l1,l2-1)*real(2* (l1+l2)-1)/real(2*l2+1)
            dfac(l2,l1) = dfac(l1,l2)
   40    continue
   30 continue       
c
c---> initialize
c
       
               do 130 lm1 = 1,lmpot
                  b(lm1) = 0.0d0
  130          continue
               do 140 lm1 = 1,lmpot
                  do 150 lm2 = 1,lmpot
                     a(lm1,lm2) = 0.0d0
  150             continue
  140          continue
c
c---> calculate b(lm1)
c
               do 160 lm1 = 1,lmpot
                  l1 = loflm(lm1)
                  b(lm1) = b(lm1) -
     +                           2.0d0*fpi/real(2*l1+1)*smat(i1,i2,lm1)
  160          continue
c
c---> calculate a(i1,i2,lm1,lm2)
c
               do 170 i = 1,iend
                  lm1 = icleb(i,1)
                  lm2 = icleb(i,2)
                  lm3 = icleb(i,3)
                  l1 = loflm(lm1)
                  l2 = loflm(lm2)
c
c---> this loop has to be calculated only for l1+l2=l3
c
                  a(lm1,lm2) = a(lm1,lm2) +
     +                               2.0d0*dfac(l1,l2)*smat(i1,i2,lm3)*
     +                               cleb(i)
  170          continue
         END
