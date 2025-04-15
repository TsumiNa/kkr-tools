     
CDECK  GRADR 
      SUBROUTINE GRADR(nspin,ist1,mesh,dx,drdi,drdi2,ro,zta,
     &  drr,ddrr,drru,ddrru,
     &  rou) 
C.....-----------------------------------------------------------------
C     evaluates d(ro)/dr,d{d(ro)/dr}/dr.
C     drr=d(ro)/dr, ddrr=d(drr)/dr.
C     coded by T.Asada. Feb.1994.
C.....-----------------------------------------------------------------
      IMPLICIT REAL*8(a-h,o-z)
C.....-----------------------------------------------------------------
c      INTEGER IRMD
c      PARAMETER(IRMD=484)
      include 'inc.fi'
      common/cxcf/igl,igh,imj,ihb,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,
     &  iex,xlf
      common/cndvpt/ndvpt
c.....------------------------------------------------------------------
      REAL*8 drdi(IRMD),drdi2(IRMD),ro(IRMD),zta(IRMD),
     &                  rou(IRMD)
      REAL*8 drr(IRMD),ddrr(IRMD),drru(IRMD),ddrru(IRMD)
c.....-----------------------------------------------------------------
C     double function
c      double precision f131,f132,f133,f141,f142,f143,f144
c      double precision fl61,fl62,fl63,fl64,fl65,fl66 
c      double precision fl51,fl52,fl53,fl54,fl55
c      double precision f231,f232,f233,f241,f242,f243,f244
c      double precision f251,f252,f253,f254,f255
c      double precision f261,f262,f263,f264,f265,f266  

      data ndvpt/6/
      data igl,igh,imj,ibh,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,iex,
     &     xlf/0,0,0,0,0,1,0,0,1,0,0,1,0,0,0.00/
c  statement functions:

c.....three point formula for the 1st deriv.
      f131(f0,f1,f2,d)=(-3*f0+4*f1  -f2)/(2*d)
      f132(g1,f0,f1,d)=(-1*g1-0*f0  +f1)/(2*d)
      f133(g2,g1,f0,d)=(   g2-4*g1+3*f0)/(2*d)
      
c.....four point formula for the 1st deriv.
      f141(f0,f1,f2,f3,d)=(-11*f0+18*f1 -9*f2 +2*f3)/(6*d)
      f142(g1,f0,f1,f2,d)=( -2*g1 -3*f0 +6*f1   -f2)/(6*d)
      f143(g2,g1,f0,f1,d)=(    g2 -6*g1 +3*f0 +2*f1)/(6*d)
      f144(g3,g2,g1,f0,d)=( -2*g3 +9*g2-18*g1+11*f0)/(6*d)

c.....five point formula for the 1st deriv.
      f151(f0,f1,f2,f3,f4,d)=(-50*f0+96*f1-72*f2+32*f3 -6*f4)/(24*d)
      f152(g1,f0,f1,f2,f3,d)=( -6*g1-20*f0+36*f1-12*f2 +2*f3)/(24*d)
      f153(g2,g1,f0,f1,f2,d)=(  2*g2-16*g1- 0*f0+16*f1 -2*f2)/(24*d)
      f154(g3,g2,g1,f0,f1,d)=( -2*g3+12*g2-36*g1+20*f0 +6*f1)/(24*d)
      f155(g4,g3,g2,g1,f0,d)=(  6*g4-32*g3+72*g2-96*g1+50*f0)/(24*d)

c.....six point formula for the 1st deriv.
      f161(f0,f1,f2,f3,f4,f5,d)=(-274*f0+600*f1-600*f2+400*f3-150*f4
     &                            +24*f5)/(120*d)
      f162(g1,f0,f1,f2,f3,f4,d)=( -24*g1-130*f0+240*f1-120*f2 +40*f3
     &                             -6*f4)/(120*d)
      f163(g2,g1,f0,f1,f2,f3,d)=(   6*g2 -60*g1 -40*f0+120*f1 -30*f2 
     &                             +4*f3)/(120*d)
      f164(g3,g2,g1,f0,f1,f2,d)=(  -4*g3 +30*g2-120*g1 +40*f0 +60*f1 
     &                             -6*f2)/(120*d)
      f165(g4,g3,g2,g1,f0,f1,d)=(   6*g4 -40*g3+120*g2-240*g1+130*f0 
     &                            +24*f1)/(120*d)
      f166(g5,g4,g3,g2,g1,f0,d)=( -24*g5+150*g4-400*g3+600*g2-600*g1
     &                           +274*f0)/(120*d)

c.....three point formula for the 2nd deriv.
      f231(f0,f1,f2,d)=(f0-2*f1+f2)/(d*d)
      f232(g1,f0,f1,d)=(g1-2*f0+f1)/(d*d)
      f233(g2,g1,f0,d)=(g2-2*g1+f0)/(d*d)

c.....four point formula for the 2nd deriv.
      f241(f0,f1,f2,f3,d)=( 6*f0-15*f1+12*f2-3*f3)/(3*d*d)
      f242(g1,f0,f1,f2,d)=( 3*g1 -6*f0 +3*f1+0*f2)/(3*d*d)
      f243(g2,g1,f0,f1,d)=( 0*g2 +3*g1 -6*f0+3*f1)/(3*d*d)
      f244(g3,g2,g1,f0,d)=(-3*g3 +2*g2+15*g1+6*f0)/(3*d*d)
   
c.....five point formula for the 2nd deriv.
      f251(f0,f1,f2,f3,f4,d)=(35*f0-104*f1+114*f2 -56*f3+11*f4)/(12*d*d)
      f252(g1,f0,f1,f2,f3,d)=(11*g1 -20*f0  +6*f1  +4*f2   -f3)/(12*d*d)
      f253(g2,g1,f0,f1,f2,d)=(  -g2 +16*g1 -30*f0 +16*f1   -f2)/(12*d*d)
      f254(g3,g2,g1,f0,f1,d)=(  -g3  +4*g2  +6*g1 -20*f0+11*f1)/(12*d*d)
      f255(g4,g3,g2,g1,f0,d)=(11*g4 -56*g3+114*g2-104*g1+35*f0)/(12*d*d)

c.....six point formula for the 2nd deriv.
      f261(f0,f1,f2,f3,f4,f5,d)=(225*f0-770*f1+1070*f2 -780*f3
     &                                         +305*f4  -50*f5)/(60*d*d)
      f262(g1,f0,f1,f2,f3,f4,d)=( 50*g1 -75*f0  -20*f1  +70*f2
     &                                          -30*f3   +5*f4)/(60*d*d)
      f263(g2,g1,f0,f1,f2,f3,d)=( -5*g2 +80*g1 -150*f0  +80*f1  
     &                                           -5*f2   +0*f3)/(60*d*d)
      f264(g3,g2,g1,f0,f1,f2,d)=(  0*g3  -5*g2  +80*g1 -150*f0
     &                                          +80*f1   -5*f2)/(60*d*d)
      f265(g4,g3,g2,g1,f0,f1,d)=(  5*g4 -30*g3  +70*g2  -20*g1 
     &                                          -75*f0  +50*f1)/(60*d*d)
      f266(g5,g4,g3,g2,g1,f0,d)=(-50*g5+305*g4 -780*g3+1070*g2  
     &                                         -770*g1 +225*f0)/(60*d*d)
      
c.....-----------------------------------------------------------------
      data iwr/0/
      if(iwr.eq.1)
     +  write(6,'(/''  igl,igh,imj,ihb,ica,icg,ivn,ipw,ipg,'',
     &  ''ivg,ip9,igd,ixlf,iex,xlf=''14i2,f10.4)') igl,igh,imj,ihb,ica,
     & icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,iex,xlf
      iwr=0

      ist=ist1
c     write(6,*) 'ndvpt ist mesh dx drdi2' ,ndvpt,ist,mesh,dx,
c    &            drdi2(ist)

      if(ndvpt.lt.3 .or. ndvpt.gt.6) then
        write(6,1265) ndvpt
 1265   format(/' ndvpt should be ge.4 .or. le.6. ndvpt=',i3)
        stop18
      endif
c.....
c.....ro: total(core+val)(up+down) charge density.

      do 50 i=ist,mesh
   50 rou(i)=ro(i)*(zta(i)+1.)/2.
c.....
      if(igd.le.0) then

        do 41 i=ist,mesh
         drr(i)=0.
         ddrr(i)=0.
         drru(i)=0.
         ddrru(i)=0.
   41   continue
        go to 200

      endif

      i1=ist
      i2=ist+1
      i3=ist+2
      i4=ist+3
      i5=ist+4
      i6=ist+5

c.....drr:d(ro)/dr, ddrr=d(d(ro)/dr)/dr
cc.... drru,ddrru: for up   spin,
c.....

      if(nspin.eq.1) go to 100
c.....
       if(ndvpt.eq.3) then

         drx1=  f131(ro(i1), ro(i2), ro(i3),dx)
         drxu1= f131(rou(i1),rou(i2),rou(i3),dx)
         drxx1= f231(ro(i1), ro(i2), ro(i3),dx)
         drxxu1=f231(rou(i1),rou(i2),rou(i3),dx)

       elseif(ndvpt.eq.4) then

         drx1=  f141(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxu1= f141(rou(i1),rou(i2),rou(i3),rou(i4),dx)
         drxx1= f241(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxxu1=f241(rou(i1),rou(i2),rou(i3),rou(i4),dx)
         drx2=  f142(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxu2= f142(rou(i1),rou(i2),rou(i3),rou(i4),dx)
         drxx2= f242(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxxu2=f242(rou(i1),rou(i2),rou(i3),rou(i4),dx)

       elseif(ndvpt.eq.5) then

         drx1=  f151(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxu1= f151(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
         drxx1= f251(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxxu1=f251(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
         drx2=  f152(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxu2= f152(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)
         drxx2= f252(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxxu2=f252(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),dx)

       elseif(ndvpt.eq.6) then

         drx1=  f161(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxu1= f161(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
         drxx1= f261(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxxu1=f261(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
         drx2=  f162(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxu2= f162(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
         drxx2= f262(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxxu2=f262(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
         drx3=  f163(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxu3= f163(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)
         drxx3= f263(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxxu3=f263(rou(i1),rou(i2),rou(i3),rou(i4),rou(i5),rou(i6),dx)

       endif

          drr(i1)=drx1/drdi(i1)
          ddrr(i1)=(drxx1-drx1*drdi2(i1))/drdi(i1)**2
          drru(i1)=drxu1/drdi(i1)
          ddrru(i1)=(drxxu1-drxu1*drdi2(i1))/drdi(i1)**2

        if(ndvpt.gt.3) then

          drr(i2)=drx2/drdi(i2)
          ddrr(i2)=(drxx2-drx2*drdi2(i2))/drdi(i2)**2
          drru(i2)=drxu2/drdi(i2)
          ddrru(i2)=(drxxu2-drxu2*drdi2(i2))/drdi(i2)**2

          if(ndvpt.eq.6) then
            drr(i3)=drx3/drdi(i3)
            ddrr(i3)=(drxx3-drx3*drdi2(i3))/drdi(i3)**2
            drru(i3)=drxu3/drdi(i3)
            ddrru(i3)=(drxxu3-drxu3*drdi2(i3))/drdi(i3)**2
          endif

        endif

        nred=DFLOAT(ndvpt)/2+.1

        do 52 j=nred+ist,mesh-nred

          if(ndvpt.eq.3) then

            drx=  f132(ro(j-1), ro(j), ro(j+1),dx)
            drxu= f132(rou(j-1),rou(j),rou(j+1),dx)
            drxx= f232(ro(j-1), ro(j), ro(j+1),dx)
            drxxu=f232(rou(j-1),rou(j),rou(j+1),dx)

          elseif(ndvpt.eq.4) then

            drx=  f142(ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxu= f142(rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
            drxx= f242(ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxxu=f242(rou(j-1),rou(j),rou(j+1),rou(j+2),dx)

          elseif(ndvpt.eq.5) then

            drx=  f153(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxu= f153(rou(j-2),rou(j-1),rou(j),rou(j+1),rou(j+2),dx)
            drxx= f253(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxxu=f253(rou(j-2),rou(j-1),rou(j),rou(j+1),rou(j+2),dx)

          elseif(ndvpt.eq.6) then

            drx=  f164(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1),
     &                                                   ro(j+2),dx)
            drxu= f164(rou(j-3),rou(j-2),rou(j-1),rou(j),rou(j+1),
     &                                                   rou(j+2),dx)
            drxx= f264(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1),
     &                                                   ro(j+2),dx)
            drxxu=f264(rou(j-3),rou(j-2),rou(j-1),rou(j),rou(j+1),
     &                                                   rou(j+2),dx)

          endif

            drr(j)=drx/drdi(j)
            ddrr(j)=(drxx-drx*drdi2(j))/drdi(j)**2
            drru(j)=drxu/drdi(j)
            ddrru(j)=(drxxu-drxu*drdi2(j))/drdi(j)**2

   52   continue
c.....
        if(ndvpt.eq.3) then

          drx0=  f133(ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxu0= f133(rou(mesh-2),rou(mesh-1),rou(mesh),dx)
          drxx0= f233(ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxxu0=f233(rou(mesh-2),rou(mesh-1),rou(mesh),dx)

        elseif(ndvpt.eq.4) then

          drx1=  f143(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxu1= f143(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
          drxx1= f243(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxxu1=f243(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
          drx0=  f144(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxu0= f144(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)
          drxx0= f244(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxxu0=f244(rou(mesh-3),rou(mesh-2),rou(mesh-1),rou(mesh),dx)

        elseif(ndvpt.eq.5) then

          drx1= f154(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)
          drxu1=f154(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1),
     &                                                rou(mesh),dx)
          drxx1=f254(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)
          drxxu1=f254(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1),
     &                                                 rou(mesh),dx)
          drx0=f155(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                               ro(mesh),dx)
          drxu0=f155(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1),
     &                                                rou(mesh),dx)
          drxx0=f255(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)
          drxxu0=f255(rou(mesh-4),rou(mesh-3),rou(mesh-2),rou(mesh-1),
     &                                                 rou(mesh),dx)

        elseif(ndvpt.eq.6) then

          drx2=f164(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxu2=f164(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                     rou(mesh-1),rou(mesh),dx)
          drxx2=f264(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)
          drxxu2=f264(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                      rou(mesh-1),rou(mesh),dx)

          drx1=f165(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxu1=f165(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                     rou(mesh-1),rou(mesh),dx)
          drxx1=f265(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)
          drxxu1=f265(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                      rou(mesh-1),rou(mesh),dx)

          drx0=f166(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxu0=f166(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                     rou(mesh-1),rou(mesh),dx)
          drxx0=f266(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)
          drxxu0=f266(rou(mesh-5),rou(mesh-4),rou(mesh-3),rou(mesh-2),
     &                                      rou(mesh-1),rou(mesh),dx)


        endif

       if(ndvpt.gt.3) then

         if(ndvpt.eq.6) then
           drr(mesh-2)=drx2/drdi(mesh-2)
            drru(mesh-2)=drxu2/drdi(mesh-2)
           ddrr(mesh-2)=(drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
           ddrru(mesh-2)=(drxxu2-drxu2*drdi2(mesh-2))/drdi(mesh-2)**2
         endif

          drr(mesh-1)=drx1/drdi(mesh-1)
          drru(mesh-1)=drxu1/drdi(mesh-1)
          ddrr(mesh-1)=(drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2
          ddrru(mesh-1)=(drxxu1-drxu1*drdi2(mesh-1))/drdi(mesh-1)**2

       endif

          drr(mesh)=drx0/drdi(mesh)
          drru(mesh)=drxu0/drdi(mesh)
          ddrr(mesh)=(drxx0-drx0*drdi2(mesh))/drdi(mesh)**2
          ddrru(mesh)=(drxxu0-drxu0*drdi2(mesh))/drdi(mesh)**2

      go to 200

  100 continue

c.....
       if(ndvpt.eq.3) then

         drx1=  f131(ro(i1), ro(i2), ro(i3),dx)
         drxx1= f231(ro(i1), ro(i2), ro(i3),dx)

       elseif(ndvpt.eq.4) then

         drx1=  f141(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxx1= f241(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drx2=  f142(ro(i1), ro(i2), ro(i3), ro(i4),dx)
         drxx2= f242(ro(i1), ro(i2), ro(i3), ro(i4),dx)

       elseif(ndvpt.eq.5) then

         drx1=  f151(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxx1= f251(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drx2=  f152(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)
         drxx2= f252(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5),dx)

       elseif(ndvpt.eq.6) then

         drx1=  f161(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxx1= f261(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drx2=  f162(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxx2= f262(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drx3=  f163(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)
         drxx3= f263(ro(i1), ro(i2), ro(i3), ro(i4), ro(i5), ro(i6),dx)

       endif

          drr(i1)=drx1/drdi(i1)
          ddrr(i1)=(drxx1-drx1*drdi2(i1))/drdi(i1)**2

        if(ndvpt.gt.3) then

          drr(i2)=drx2/drdi(i2)
          ddrr(i2)=(drxx2-drx2*drdi2(i2))/drdi(i2)**2

          if(ndvpt.eq.6) then
            drr(i3)=drx3/drdi(i3)
            ddrr(i3)=(drxx3-drx3*drdi2(i3))/drdi(i3)**2
          endif

        endif

        nred=DFLOAT(ndvpt)/2+.1

        if(mesh-nred .le. ist) then
          write(6,'(/'' mesh-nred.lt.ist. mesh,nred,ist=''3i4)') mesh,
     &      nred,ist
          stop13
        endif

        do 53 j=nred+ist,mesh-nred

          if(ndvpt.eq.3) then

            drx=  f132(ro(j-1), ro(j), ro(j+1),dx)
            drxx= f232(ro(j-1), ro(j), ro(j+1),dx)

          elseif(ndvpt.eq.4) then

            drx=  f142(ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxx= f242(ro(j-1), ro(j), ro(j+1), ro(j+2),dx)

          elseif(ndvpt.eq.5) then

            drx=  f153(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2),dx)
            drxx= f253(ro(j-2), ro(j-1), ro(j), ro(j+1), ro(j+2),dx)

          elseif(ndvpt.eq.6) then

            drx=  f164(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1),
     &                                                   ro(j+2),dx)
            drxx= f264(ro(j-3), ro(j-2), ro(j-1), ro(j), ro(j+1),
     &                                                   ro(j+2),dx)

          endif

            drr(j)=drx/drdi(j)
            ddrr(j)=(drxx-drx*drdi2(j))/drdi(j)**2
c           write(6,9000) j,drr(j)
c9000       format(1x,' j drr(j)',i5,e15.5)
   53   continue
c.....
        if(ndvpt.eq.3) then

          drx0=  f133(ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxx0= f233(ro(mesh-2), ro(mesh-1), ro(mesh),dx)

        elseif(ndvpt.eq.4) then

          drx1=  f143(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxx1= f243(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drx0=  f144(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)
          drxx0= f244(ro(mesh-3), ro(mesh-2), ro(mesh-1), ro(mesh),dx)

        elseif(ndvpt.eq.5) then

          drx1= f154(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)
          drxx1=f254(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)
          drx0=f155(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                               ro(mesh),dx)
          drxx0=f255(ro(mesh-4), ro(mesh-3), ro(mesh-2), ro(mesh-1),
     &                                                ro(mesh),dx)

        elseif(ndvpt.eq.6) then

          drx2=f164(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxx2=f264(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)

          drx1=f165(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxx1=f265(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)

          drx0=f166(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                    ro(mesh-1), ro(mesh),dx)
          drxx0=f266(ro(mesh-5), ro(mesh-4), ro(mesh-3), ro(mesh-2),
     &                                     ro(mesh-1), ro(mesh),dx)


        endif

       if(ndvpt.gt.3) then

         if(ndvpt.eq.6) then
           drr(mesh-2)=drx2/drdi(mesh-2)
           ddrr(mesh-2)=(drxx2-drx2*drdi2(mesh-2))/drdi(mesh-2)**2
         endif

          drr(mesh-1)=drx1/drdi(mesh-1)
          ddrr(mesh-1)=(drxx1-drx1*drdi2(mesh-1))/drdi(mesh-1)**2

       endif

          drr(mesh)=drx0/drdi(mesh)
          ddrr(mesh)=(drxx0-drx0*drdi2(mesh))/drdi(mesh)**2

  200 continue
C
C
C
C      write(6,8000) nspin,ist1,mesh,dx  
C8000 format(1x,' nspin ist1 mesh dx',3i5,2d20.10)
C     write(6,8001) (ro(kk),drr(kk),ddrr(kk),
c    &  drdi(kk),drdi2(kk), kk=ist1,mesh,20)
c8001 format(1x,' ro drr ddrr drdi drdi2',5f12.5)
      return
      end
