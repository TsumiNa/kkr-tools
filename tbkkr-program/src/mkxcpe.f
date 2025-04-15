*DECK MKXCPE
c
      subroutine mkxcpe(nspin,ir,np,l1max,rv,rl,vxcp,excp)
c.....------------------------------------------------------------------
      implicit real*8(a-h,o-z)                                          
      INTEGER IRMD,LPOTD,NSPIND
      parameter(IRMD=884,LPOTD=8,NSPIND=2)
      INTEGER LMPOTD,N,IJD
C     parameter (LMPOTD=(LPOTD+1)**2,N=2*(LPOTD+1),IJD=2*N**2)
      parameter (LMPOTD=(LPOTD+1)**2,N=2*(LPOTD+1),IJD=434)
C     INTEGER L3D,LM3D,NLCEB
C     parameter (L3D=2*LPOTD,LM3D=(L3D+1)**2,NCLEB=LM3D*LMPOTD)

c.....------------------------------------------------------------------
c.....------------------------------------------------------------------
      common/cylmdnp/thet(IJD),ylm(IJD,LMPOTD),dylmt1(IJD,LMPOTD),
     &  dylmt2(IJD,LMPOTD),dylmf1(IJD,LMPOTD),dylmf2(IJD,LMPOTD),
     &  dylmtf(IJD,LMPOTD)
C
      common/cgradr/drrs,ddrrs,
     &        drrus,ddrrus

      real*8 drrs(IRMD,LMPOTD),ddrrs(IRMD,LMPOTD),
     &       drrus(IRMD,LMPOTD),ddrrus(IRMD,LMPOTD)
      real*8 vxcp(ijd,nspind),excp(ijd)
      real*8 rl(LMPOTD,nspind)
c.....------------------------------------------------------------------
      real*8 ry(IJD),dry(IJD),ddry(IJD),
     &  rdt1(IJD),rdt2(IJD),rdf1(IJD),
     &  rdf2(IJD),rdtf(IJD),drdt(IJD),
     &  drdf(IJD),
     &  ryu(IJD),dryu(IJD),ddryu(IJD),
     &  rdt1u(IJD),rdt2u(IJD),rdf1u(IJD),
     &  rdf2u(IJD),rdtfu(IJD),drdtu(IJD),
     &  drdfu(IJD),
     &  ryd(IJD),dryd(IJD),ddryd(IJD),
     &  rdt1d(IJD),rdt2d(IJD),rdf1d(IJD),
     &  rdf2d(IJD),rdtfd(IJD),drdtd(IJD),
     &  drdfd(IJD)

      common /lebedev/ICHECK(IJD)


       real*8 rv
      data rdspr,zero,zero1/9.0,0.e0,1.e-12/
c.....------------------------------------------------------------------
c     rl: charge=sumlm(rl*ylm)
c     ry=sumlm(ro*ylm), dry=sumlm(drr*ylm), ddry=sumlm(ddrr*ylm), 
cc    rdt1=sumlm(ro*dylmt1), rdt2=sumlm(ro*dylmt2), ...
cc    rdf1=sumlm(ro*dylmf1), rdf2=sumlm(ro*dylmf2), ...
cc    rdtf=sumlm(ro*dylmtf), rdf2=sumlm(ro*dylmf2), ...
cc    drdt=sumlm(drr*dylmt1),drdf=sumlm(drr*dylmf1),

c     agr: abs(grad(ro)), g2r: laplacian(ro), 
cc    gggr: grad(ro)*grad(agr),
cc    grgru,d: grad(ro)*grad(rou),for rod., gzgr: grad(zeta)*grad(ro). 

c     dagrr,-t,-f: d(agr)/dr, d(agr)/dth, d(agr)/dfi.
c.....------------------------------------------------------------------
c     if(meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or. 
c    &   mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD) then
c       write(6,'(/'' meshx.ne.IRMD .or. lLMPOTD.ne.LMPOTD .or. '',
c    &    ''mesh.gt.IRMD  .or. l1max.gt.LMPOTD .or. np.gt.IJD.''/
c    &    '' meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD='',10i4)'
c    &    ) meshx,IRMD,lLMPOTD,LMPOTD,mesh,IRMD,l1max,LMPOTD,np,IJD
c       stop14
c     endif
check  ist=mesh
      llmax=l1max*l1max
      lmax=l1max-1
c     lmax2=lmax*2
c     llmax2=(lmax2+1)**2
c     lmax3=lmax*1
c     llmax3=(lmax3+1)**2 
c     
C     write(6,*) ' r',rv
c
c
c
c     if(ir.eq.375)
c    & write(6,9030) (ii,drrs(ii,1),ddrrs(ii,1),drrus(ii,1),
c    &               ddrrus(ii,1),ii=ir,ir)
c9030 format(1x,' ist drrs  ddrrs drrus ddrrus',i5,4e12.5)
c


      do 20 ip=1,np
       
       ry(ip)=0.
       dry(ip)=0.
       ddry(ip)=0.
       rdt1(ip)=0.
       rdt2(ip)=0.
       rdf1(ip)=0.
       rdf2(ip)=0.
       rdtf(ip)=0.
       drdt(ip)=0.
       drdf(ip)=0.

       ryu(ip)=0.
       dryu(ip)=0.
       ddryu(ip)=0.
       rdt1u(ip)=0.
       rdt2u(ip)=0.
       rdf1u(ip)=0.
       rdf2u(ip)=0.
       rdtfu(ip)=0.
       drdtu(ip)=0.
       drdfu(ip)=0.

       ryd(ip)=0.
       dryd(ip)=0.
       ddryd(ip)=0.
       rdt1d(ip)=0.
       rdt2d(ip)=0.
       rdf1d(ip)=0.
       rdf2d(ip)=0.
       rdtfd(ip)=0.
       drdtd(ip)=0.
       drdfd(ip)=0.

   20  continue
C     write(6,*) ' in mkxcpe ' 
C     write(6,'(/'' nspin,mesh,np,l1max='',4i5)') nspin,mesh,np,l1max

      lm=0

      do 30 l1=1,l1max

        ll=l1-1

      

      do 50 im=-ll,ll

        lm=lm+1


         ro=rl(lm,1)*2.e0
         rou=ro/2.
         rod=rou
c
          if(nspin.ne.1) then
c
           ro=rl(lm,1)+rl(lm,2)
           rou=rl(lm,2)
           rod=rl(lm,1)
c        write(6,9001) ro,rou,rod
 9001    format(1x,' ro rou ro ',3e15.5)
          endif 
            drr=drrs(ir,lm)
            ddrr=ddrrs(ir,lm)
           drru=drrus(ir,lm)
           ddrru=ddrrus(ir,lm)
         drrd=drr-drru
        ddrrd=ddrr-ddrru


          do 60 ip=1,np

            rylm=ylm(ip,lm)
            dt1=dylmt1(ip,lm)
            dt2=dylmt2(ip,lm)
            df1=dylmf1(ip,lm)
            df2=dylmf2(ip,lm)
            dtf=dylmtf(ip,lm)

c        IF(ip.eq.5) write(6,9000) rylm,dt1,dt2,df1,df2,dtf
c9000    format(1x,' rylm dt1 dt2 df1 df2 dtf ',6e12.4)

            ry(ip)=ry(ip)+ro*rylm
            dry(ip)=dry(ip)+drr*rylm
            ddry(ip)=ddry(ip)+ddrr*rylm

            ryu(ip)=ryu(ip)+rou*rylm
            dryu(ip)=dryu(ip)+drru*rylm
            ddryu(ip)=ddryu(ip)+ddrru*rylm

            ryd(ip)=ryd(ip)+rod*rylm
            dryd(ip)=dryd(ip)+drrd*rylm
            ddryd(ip)=ddryd(ip)+ddrrd*rylm

              rdt1(ip)=rdt1(ip)+ro*dt1
              rdt2(ip)=rdt2(ip)+ro*dt2
              rdf1(ip)=rdf1(ip)+ro*df1
              rdf2(ip)=rdf2(ip)+ro*df2
              rdtf(ip)=rdtf(ip)+ro*dtf
              drdt(ip)=drdt(ip)+drr*dt1
              drdf(ip)=drdf(ip)+drr*df1

              rdt1u(ip)=rdt1u(ip)+rou*dt1
              rdt2u(ip)=rdt2u(ip)+rou*dt2
              rdf1u(ip)=rdf1u(ip)+rou*df1
              rdf2u(ip)=rdf2u(ip)+rou*df2
              rdtfu(ip)=rdtfu(ip)+rou*dtf
              drdtu(ip)=drdtu(ip)+drru*dt1
              drdfu(ip)=drdfu(ip)+drru*df1

              rdt1d(ip)=rdt1d(ip)+rod*dt1
              rdt2d(ip)=rdt2d(ip)+rod*dt2
              rdf1d(ip)=rdf1d(ip)+rod*df1
              rdf2d(ip)=rdf2d(ip)+rod*df2
              rdtfd(ip)=rdtfd(ip)+rod*dtf
              drdtd(ip)=drdtd(ip)+drrd*dt1
              drdfd(ip)=drdfd(ip)+drrd*df1

c             write(6,9907) ro,rou,rod
c9907         format(1x,' ro rou rod',3e12.6) 


   60     continue
  50      continue
 30       continue

C        write(6,*) ' after 30 continue '
c


      do 80 ip=1,np
       sint1=sin(thet(ip))
       if(abs(sint1).le.1.E-10) go to 80
       sint2=sint1**2
       tant1=tan(thet(ip))

c  999
c     write(6,'(/'' ip,sint1,2,tnt1='',i3,3f12.7)') ip,sint1,sint2,tant1


       rv2=rv**2
       rv3=rv**3


       rvsin1=rv*sint1
       rvsin2=rv2*sint2 

       grr=dry(ip)
       grt=rdt1(ip)/rv
       grf=rdf1(ip)/rvsin1
       ry2=ry(ip)**2

       agr=sqrt(grr**2+grt**2+grf**2)
C      if(ip.le.6) write(6,9955) ip,grr,grt,grf,rv,thet(ip)
C9955  format(1x,' ip grr grt grf rv   thet',I5,6E15.3)

        dagrr=( dry(ip)*ddry(ip)*rv3+
     &          rdt1(ip)*(drdt(ip)*rv-rdt1(ip))+
     &          rdf1(ip)*(drdf(ip)*rv-rdf1(ip))/sint2
     &        )/agr/rv3

        dagrt=( dry(ip)*drdt(ip)*rv2+
     &          rdt1(ip)*rdt2(ip)+
     &          rdf1(ip)*(-rdf1(ip)/tant1+rdtf(ip))/sint2
     &        )/(agr*rv3)

        dagrf=( dry(ip)*drdf(ip)*rv2+
     &          rdt1(ip)*rdtf(ip)+
     &          rdf1(ip)*rdf2(ip)/sint2
     &        )/(agr*rv3*sint1)


        dzdr=( (dryu(ip)-dryd(ip))*ry(ip)-
     &         (ryu(ip)-ryd(ip))*dry(ip)
     &       )/ry2

        dzdtr=( (rdt1u(ip)-rdt1d(ip))*ry(ip)-
     &          (ryu(ip)-ryd(ip))*rdt1(ip)
     &        )/ry2/rv

        dzdfs=( (rdf1u(ip)-rdf1d(ip))*ry(ip)-
     &           (ryu(ip)-ryd(ip))*rdf1(ip)
     &         )/ry2/rvsin1

       g2r=ddry(ip)+2.*dry(ip)/rv + 
     &   (rdt2(ip)+rdt1(ip)/tant1+rdf2(ip)/sint2)/rv2

       gggr=grr*dagrr+grt*dagrt+grf*dagrf 

       gzgr=dzdr*grr+dzdtr*grt+dzdfs*grf      

c       
       chg=ry(ip) 
       spi=ryu(ip)-ryd(ip)
       chg=max(1.0e-12,chg)
       smag=sign(1.0e0,real(spi))
       spi=smag*min(chg-1.0e-12,abs(spi))
       zta=spi/chg
c

       grru=dryu(ip)
       grtu=rdt1u(ip)/rv
       grfu=rdf1u(ip)/rvsin1

       agru=sqrt(grru**2+grtu**2+grfu**2)

        dagrru=( dryu(ip)*ddryu(ip)*rv3+
     &           rdt1u(ip)*(drdtu(ip)*rv-rdt1u(ip))+
     &           rdf1u(ip)*(drdfu(ip)*rv-rdf1u(ip))/sint2
     &         )/agru/rv3

        dagrtu=( dryu(ip)*drdtu(ip)*rv2+
     &           rdt1u(ip)*rdt2u(ip)+
     &           rdf1u(ip)*(-rdf1u(ip)/tant1+rdtfu(ip))/sint2
     &         )/(agru*rv3)

        dagrfu=( dryu(ip)*drdfu(ip)*rv2+
     &           rdt1u(ip)*rdtfu(ip)+
     &           rdf1u(ip)*rdf2u(ip)/sint2
     &         )/(agru*rv3*sint1)



       g2ru=ddryu(ip)+2.*dryu(ip)/rv + 
     &   (rdt2u(ip)+rdt1u(ip)/tant1+rdf2u(ip)/sint2)/rv2

       gggru=grru*dagrru+grtu*dagrtu+grfu*dagrfu 

       grgru=grr*grru+grt*grtu+grf*grfu


        grrd=dryd(ip)
        grtd=rdt1d(ip)/rv
        grfd=rdf1d(ip)/rvsin1

       agrd=sqrt(grrd**2+grtd**2+grfd**2)

        dagrrd=( dryd(ip)*ddryd(ip)*rv3+
     &           rdt1d(ip)*(drdtd(ip)*rv-rdt1d(ip))+
     &           rdf1d(ip)*(drdfd(ip)*rv-rdf1d(ip))/sint2
     &         )/agrd/rv3

        dagrtd=( dryd(ip)*drdtd(ip)*rv2+
     &           rdt1d(ip)*rdt2d(ip)+
     &           rdf1d(ip)*(-rdf1d(ip)/tant1+rdtfd(ip))/sint2
     &         )/(agrd*rv3)

        dagrfd=( dryd(ip)*drdfd(ip)*rv2+
     &           rdt1d(ip)*rdtfd(ip)+
     &           rdf1d(ip)*rdf2d(ip)/sint2
     &         )/(agrd*rv3*sint1)



       g2rd=ddryd(ip)+2.*dryd(ip)/rv + 
     &   (rdt2d(ip)+rdt1d(ip)/tant1+rdf2d(ip)/sint2)/rv2

       gggrd=grrd*dagrrd+grtd*dagrtd+grfd*dagrfd 

       grgrd=grr*grrd+grt*grtd+grf*grfd


      idspr=0
      if(rv.gt.rdspr) idspr=1
      
c
       if(ip.eq.5.and.ir.eq.2)
     & write(6,9911) agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,
     &              gggrd,grgru,grgrd,gzgr
c
c
c
c for debug
C     write(6,*) ' before gxcpt '
       call gxcpt(idspr,chg,zta,agr,agru,
     &  agrd, g2r,g2ru,g2rd,
     &  gggr,gggru,gggrd, grgru,grgrd,gzgr,
     &  vxcp(ip,2),vxcp(ip,1),excp(ip),
     &  vxl1,vxl2,vcl1,vcl2,xedl,
     &  cedl,
     &  vxg1,vxg2,vcg1,vcg2,xedg,
     &  cedg)
c
      GO TO 82
c
   81 vxcp(ip,1) = vxcp(ip-1,1)
      vxcp(ip,2) = vxcp(ip-1,2)
      excp(ip)   = excp(ip-1)
c
   82 continue
c
      if(ir.eq.2.and.ip.le.1) then
      write(6,9912) ir,ip,ry(ip),zta
 9912 format(1x,' ir ip ry zta',2i5,2e15.6)
c
      write(6,9911) agr,agru,agrd,g2r,g2ru,g2rd,gggr,gggru,
     &              gggrd,grgru,grgrd,gzgr
 9911 format(1x,' agr  ',6E15.6)

      write(6,7777) vxcp(ip,1),vxcp(ip,2),
     &              vxl1,vxl2,vcl1,vcl2
      write(6,7778) vxg1,vxg2,vcg1,vcg2
 7777 format(1x,'vxcp(1,2) vxl(1,2) vcl(1,2) ',6E15.6)
 7778 format(1x,'vxg(1,2) vcg(1,2)  (asada)  ',4E15.6)

      end if
   80 continue
c
c    interpolation
       do 85 ip=1,np
       IF(icheck(ip).eq.0) go to 85
       jj=icheck(ip)
       vxcp(ip,1)= vxcp(jj,1)
       vxcp(ip,2)= vxcp(jj,2)
       excp(ip)= excp(jj)
C      if(ip.ne.5) go to 85
C     write(6,7799) ip,jj,vxcp(ip,1),vxcp(ip,2),excp(ip)
C7799 format(1x,' interpolation',2I5,3D12.5)
   85 continue

c
c
c     if(ir.eq.2) stop 'mkxcpe'
c
c  check
c
c     iir=10
c     llm=4
c     write(6,9998) iir,llm
c9998  format(1x,' iir llm',10i5)
c     write(6,9999) thet(iir),ylm(iir,llm),dylmtf(iir,llm)
c    &   ,dylmf2(iir,llm),dylmf1(iir,llm),dylmt2(iir,llm),
c    &    dylmt1(iir,llm)
c9999 format(1x,' thet ylm dylmtf',6f10.5) 
c     write(6,9997) drrs(iir,llm),ddrrs(iir,llm)
c9997 format(1x,' drrs ddrrs',10f10.5)
c    
C     stop99
c
      return
      end 
