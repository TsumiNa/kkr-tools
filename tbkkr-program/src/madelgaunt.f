      
      
      
      
      
       SUBROUTINE MADELGAUNT(LPOT,YR,W,CLEB,ICLEB,IEND)
       implicit none
       include 'inc.fi'
       INTEGER LASSLD
       PARAMETER (LASSLD=4*LMAXD)
       integer lmpotd,l3d,lm3d
       parameter (lmpotd= (lpotd+1)**2,l3d=2*lpotd,lm3d= (l3d+1)**2)
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
       integer ncleb1
       parameter (ncleb1=lm3d*lmpotd)
      
       INTEGER LPOT 
       double precision w(LASSLD),yr(LASSLD,0:LASSLD,0:LASSLD)
       double precision cleb(ncleb1)
       integer icleb(ncleb1,3)
       double precision  clecg,factor,s
       integer i,i1,i2,iend,j,l,l1,l2,l3,lm1,lm2,lm3,lmpot,m,m1,m1a,m1s,
     &        m2,m2a,m2s,m3,m3a,m3s,n 
c     ..
c     .. intrinsic functions ..
      intrinsic abs,atan,real,sign                    
c

      lmpot = (lpot+1)**2

c
c---> set up of the gaunt coefficients with an index field
c     recognize that they are needed here only for l3=l1+l2
c
      if (2*lpot.gt.LASSLD) then
       write(6,*) 'Dim ERROR in MADELGAUNT'
       STOP
      end if
c
      i = 1
      do 50 l1 = 0,lpot
         do 60 l2 = 0,lpot
            l3 = l1 + l2
            do 70 m1 = -l1,l1
               do 80 m2 = -l2,l2
                  do 90 m3 = -l3,l3
                     m1s = sign(1,m1)
                     m2s = sign(1,m2)
                     m3s = sign(1,m3)
c
                     if (m1s*m2s*m3s.ge.0) then
c
                        m1a = abs(m1)
                        m2a = abs(m2)
                        m3a = abs(m3)
c
                        factor = 0.0d0
c
                        if (m1a+m2a.eq.m3a) factor = factor +
     +                      real(3*m3s+sign(1,-m3))/8.0d0
                        if (m1a-m2a.eq.m3a) factor = factor +
     +                      real(m1s)/4.0d0
                        if (m2a-m1a.eq.m3a) factor = factor +
     +                      real(m2s)/4.0d0
c
                        if (factor.ne.0.0d0) then
c
                           if (m1s*m2s.ne.1 .or. m2s*m3s.ne.1 .or.
     +                         m1s*m3s.ne.1) factor = -factor
c
                           s = 0.0d0
                           do 100 j = 1,LASSLD
                              s = s + w(j)*yr(j,l1,m1a)*yr(j,l2,m2a)*
     +                            yr(j,l3,m3a)
  100                      continue
                           clecg = s*factor
                           if (abs(clecg).gt.1.d-10) then
                              cleb(i) = clecg
                              icleb(i,1) = l1* (l1+1) + m1 + 1
                              icleb(i,2) = l2* (l2+1) + m2 + 1
                              icleb(i,3) = l3* (l3+1) + m3 + 1
                              i = i + 1
                           end if
 
                        end if
 
                     end if
 
   90             continue
 
   80          continue
   70       continue
   60    continue
   50 continue
      iend = i - 1
      if (ncleb1.lt.iend) then
         write (6,fmt='(i10)') iend,ncleb
         stop ' Dim stop in MADELGAUNT '
 
      END IF
      END
