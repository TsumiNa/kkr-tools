c
c for debug
c
      SUBROUTINE GXCPT(idspr,ro,zta,  agr,agru,agrd,  g2r,g2ru,g2rd,
     &  gggr,gggru,gggrd,  grgru,grgrd,gzgr,  xcptu,xcptd,xced,
     &  vxlu,vxld,vclu,vcld,xedl,cedl,
     &  vxgu,vxgd,vcgu,vcgd,xedg,cedg )

c.....-----------------------------------------------------------------
c.....gxcp: exchange-correlation potential in ry. also total-energy.
c.....-----------------------------------------------------------------
      IMPLICIT REAL*8 (a-h,o-z)
c     common/cxcf/igl,igh,imj,ibh,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,
c    &  iex,xlf
c     common/ctrns7/hugeo,huges,hugef,dspr,rdspr,idspr
c.....-----------------------------------------------------------------
      COMMON/gga/igga
      fncf(x)=(1.d0+x*x*x)*LOG(1.d0+1.d0/x)+x/2.d0-x*x-0.333333333d0
      fncecl(r,g,b1,b2)=g/(1.d0+b1*SQRT(r)+b2*r)
      fncvcl(ce,r,b1,b2)=ce*(1.d0+1.16666667d0*b1*SQRT(r)+1.33333333d0*
     &                 b2*r)/(1.d0+b1*SQRT(r)+b2*r)
      fncecs(r,a,b,c,d)=a*LOG(r)+b+c*r*LOG(r)+d*r
      fncvcs(r,a,b,c,d)=a*LOG(r)+(b-a/3.d0)+0.666666667d0*c*r*LOG(r)+
     &                      (2.d0*d-c)*r/3.d0
      ffz(zta)=1.923661051d0*((1.d0+zta)**1.3333333333d0+
     &                        (1.d0-zta)**1.3333333333d0-2.d0)
      fdfdz(zta)=2.564881401d0*((1.d0+zta)**.333333333333d0-(1.d0-zta)**
     & .333333333333d0)
      fvq(b,c)=SQRT(4.d0*c-b**2)
      fvnec(a,x,xl,x0,xl0,b,q)=a*(LOG(x*x/xl)+2.d0*b/q*DATAN(q/(2.d0*x+
     & b))-b*x0/xl0*(LOG((x-x0)**2/xl)+2.d0*(b+2.d0*x0)/q*DATAN(q/(2.d0*
     & x+b))))
      fbet(fdd0,ecf,ecp,alc)=fdd0*(ecf-ecp)/alc-1.d0
      fdedr(ro,x,a,x0,xl,xl0,xld,b,q)=-x/(6.d0*ro)*a*((2.d0*xl-x*xld)/
     & (x*xl)-b*(4.d0/(xld**2+q**2)+x0/xl0*((2.d0*xl-(x-x0)*xld)/
     & ((x-x0)*xl)-4.d0*(b+2.d0*x0)/(xld**2+q**2))))
c.....-----------------------------------------------------------------
c.....Perdew-zunger parametrization of Ceperley-Alder. g,a,b,c,d in ry.
      data gp,gf,b1p,b1f,b2p,b2f,cp,cf,dp,df/-.2846d0,-.1686d0,
     &  1.0529d0,1.3981d0,  0.3334d0,0.2611d0,  0.0040d0,0.0014d0,
     &  -.0232d0,-.0096d0/
      data ap,bp,af,bf/0.0622d0,-.096d0,  0.0311d0,-0.0538d0/
c.....for vwn. a1,x01,b1,c1 for para(zta=0), -2 for zta=1, -3 for alfac.
      data a1,x01,b1,c1/.0621814d0,-.10498d0,3.72744d0,12.9352d0/
      data a2,x02,b2,c2/.0310907d0,-.32500d0,7.06042d0,18.0578d0/
      data a3,x03,b3,c3/-.03377373d0,-.0047584d0,1.13107d0,13.0045d0/
c.....fdd0: fz''(o).
      data fdd0/1.70992093d0/
c.....-----------------------------------------------------------------
      data hugeo,huges,hugef,dspr/625.,1.d+6,50.,1.d-4/
      data igl,igh,imj,ibh,ica,icg,ivn,ipw,ipg,ivg,ip9,igd,ixlf,iex,xlf
     &     /0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00/
C
C     PW91   ip9=1  igd=1
C
C     PW86      imj=1  igd=1
C

      ip9=1
      igd=1
C
C
      pi=ACOS(-1.d0)
      sml=1.d-12
      if(zta.gt. 1.d0-sml) zta= 1.d0-sml
c.....
c     vxlu,vxld,vxgu,vxgd: exchange potential in ry.(local,grad),(up,dw)
c     vclu,vcld,vcgu,vcgd: correl. potential in ry.(local,grad),(up,dw)
c     xedl,xedg: exchange energy density (local,grad.exp.) in ry.
c     cedl,cedg: exchange energy density (local,grad.expnd.) in ry.
c.....
      vxlu=0.0
      vclu=0.0
      vxld=0.0
      vcld=0.0
      xedl=0.0
      cedl=0.0
      vxgu=0.0
      vcgu=0.0
      vxgd=0.0
      vcgd=0.0
      xedg=0.0
      cedg=0.0
c.....
      if(ro.lt.sml) go to 200
c.....
      c13=1./3.
      c23=2./3.
      c32=3./2.
      c43=4./3.
      c53=5./3.
      c76=7./6.
      c113=11./3.
c.....ca=2.**(-.33333333)
      ca=0.793700526d0
c.....alf=-3*(3/4*pai)**(1/3).
      alf=-1.861051473d0
c.....
      ro2=ro*ro
      ro13=ro**c13
      ro43=ro**c43
      ro83=ro43**2
      ro76=ro**c76
      ro113=ro**c113
c.....
      rou=ro*(1.+zta)/2.d0
      rou3=rou**3
      rou13=rou**c13
      rou23=rou**c23
      rou43=rou**c43
c.....
      rod=ro-rou
      rod3=rod**3
      rod13=rod**c13
      rod23=rod**c23
      rod43=rod**c43
c.....
c     gr2=drr*drr
c     gr2u=drru**2
c     drrd=drr-drru
c     gr2d=drrd**2
c     ddrrd=ddrr-ddrru
c.....
      fz=ffz(zta)
      dfdz=fdfdz(zta)
c.....
c.....gz,gz2,gz3: for Wang-Perdew ssf.
      gz=((1.+zta)**c23+(1.-zta)**c23)/2.
      gz2=gz**2
      gz3=gz**3
c.....
      zta3=zta**3
      zta4=zta**4
      zt13p=(1.+zta)**c13
      zt13m=(1.-zta)**c13
c.....
      rs=0.620350491d0/ro13
      rs2=rs*rs
      rs3=rs*rs2
c.....
c.....xedl: exchange-energy-density in ry.
      xedl=alf*(rou43+rod43)
c.....
c.....exchange-potential, vxp,vxlu,vxld: v-exchange-(para,up,dw).
      vxlu=c43*alf*rou13
      vxld=c43*alf*rod13

c.....
      if(iex.eq.1) go to 100
c.....

c.....xlfa.
      if(ixlf.ne.0) then
        xlf=2./3.
        vclu=(xlf*c32-1.)*vxlu
        vcld=(xlf*c32-1.)*vxld
        cedl=(xlf*c32-1.)*xedl
        go to 100
      endif
c.....
c.....Gunnarson-Lundqvist.(p.r.b13('76),4274,eqs(54)~(56).) or
cc....  g-l but with beta by hedin-lundq.(j.phys.c.4('71),2064))
      if((igl.ne.0).or.(igh.ne.0)) then
        xp=rs/11.4d0
        xf=rs/15.9d0
        cep=-0.0666d0*fncf(xp)
        cef=-0.0406d0*fncf(xf)
        ce=cep+(cef-cep)*fz
        cedl=ce*ro
c.....
        if(igl.ne.0) beta=1.+0.0545d0*rs*log(1.+11.4d0/rs)
        if(igh.ne.0) beta=1.+.03683d0*rs*log(1.+21./rs)
        dlta=1.-0.036d0*rs+1.36d0*rs/(1.+10.*rs)
c.....
        vxp=c43*alf*ca*ro13
        vclu=vxp*(beta+dlta/3.*zta/(1.+0.297d0*zta))-vxlu
        vcld=vxp*(beta-dlta/3.*zta/(1.-0.297d0*zta))-vxld
c.....
        go to 100
      endif
c.....
c.....Hedin-von Barth. (j.phys.c.5('72),1629) or moruzzi-janak-williams.
      if((ibh.ne.0).or.(imj.ne.0)) then
          if(ibh.ne.0) then
            crp=30.
            crf=75.
            ccp=0.0504d0
            ccf=0.0254d0
          elseif(imj.ne.0) then
            crp=21.
            crf=52.916684d0
            ccp=0.045d0
            ccf=0.0225d0
c           write(6,*) 'MJW'
          endif
        xp=rs/crp
        xf=rs/crf
        cep=-ccp*fncf(xp)
        cef=-ccf*fncf(xf)
        ce=cep+(cef-cep)*fz
        cedl=ce*ro
c       vclu,vcld: v-correlation-(up,dw). potential.(ry)
        rnc=c43*ca/(1.-ca)*(cef-cep)
        vcp=-ccp*log(1.+crp/rs)
        brs=vcp-rnc
        vclu=rnc*zt13p+brs
        vcld=rnc*zt13m+brs
c.....
        go to 100
      endif
c.....
c.....Ceperley-Alder.(paramtrzd by Perdew-zunger.(p.r.23('81),5048)).
      if(ica.ne.0) then
c.....
        if(rs.ge.1.d0) then
          cep=fncecl(rs,gp,b1p,b2p)
          cef=fncecl(rs,gf,b1f,b2f)
          vcp=fncvcl(cep,rs,b1p,b2p)
          vcf=fncvcl(cef,rs,b1f,b2f)
        else
          cep=fncecs(rs,ap,bp,cp,dp)
          cef=fncecs(rs,af,bf,cf,df)
          vcp=fncvcs(rs,ap,bp,cp,dp)
          vcf=fncvcs(rs,af,bf,cf,df)
        endif
c.....
        ce=cep+(cef-cep)*fz
        cedl=ce*ro
c.....
c.....
        vcl2=(cef-cep)*dfdz
        vcl1=vcp+(vcf-vcp)*fz-vcl2*zta
          vclu=vcl1+vcl2
          vcld=vcl1-vcl2
c.....
        go to 100
      endif
c.....
c.....Ceperley-Alder.with Wang-Perdew spin-scaling-factor.
      if(icg.ne.0) then
c.....
        if(rs.ge.1.d0) then
          cep=fncecl(rs,gp,b1p,b2p)
          vcp=fncvcl(cep,rs,b1p,b2p)
        else
          cep=fncecs(rs,ap,bp,cp,dp)
          vcp=fncvcs(rs,ap,bp,cp,dp)
        endif
c.....
        ce=cep*gz3
        cedl=ce*ro
c.....
        cgz=cep*gz2*(1./zt13p-1./zt13m)
        vcl1=vcp*gz3-cgz*zta
        vclu=vcp*gz3+cgz
        vcld=vcp*gz3-cgz
c.....
        go to 100
      endif
c.....
c.....Vosko-Wilk-Nusair. Phys.Rev..22,3812,'80.
      if(ivn.ne.0) then
c.....
c.....xl:x-large. xld:d(xl)/dx. xl0:x-large for x=x0.
        xs=sqrt(rs)
          xl1=xs**2+b1*xs+c1
          xl2=xs**2+b2*xs+c2
          xl3=xs**2+b3*xs+c3
           xld1=2.*xs+b1
           xld2=2.*xs+b2
           xld3=2.*xs+b3
          xl01=x01**2+b1*x01+c1
          xl02=x02**2+b2*x02+c2
          xl03=x03**2+b3*x03+c3
            q1=fvq(b1,c1)
            q2=fvq(b2,c2)
            q3=fvq(b3,c3)
       ecp=fvnec(a1,xs,xl1,x01,xl01,b1,q1)
       ecf=fvnec(a2,xs,xl2,x02,xl02,b2,q2)
       alc=fvnec(a3,xs,xl3,x03,xl03,b3,q3)
         beta=fbet(fdd0,ecf,ecp,alc)
         bz41=1.+beta*zta4
c.....
        ce=ecp+alc*fz/fdd0*bz41
        cedl=ce*ro
c.....
c.....alc: alfac.
c.....decdrp,decdrf: d(ec)/dro-para(zta=0), -(zta=1).
c.....dacdr: d(alc)/dro.
c.....dbdr: d(beta)/dro.
        decdrp=fdedr(ro,xs,a1,x01,xl1,xl01,xld1,b1,q1)
        decdrf=fdedr(ro,xs,a2,x02,xl2,xl02,xld2,b2,q2)
        dacdr =fdedr(ro,xs,a3,x03,xl3,xl03,xld3,b3,q3)
c.....
        dbdr=fdd0*((decdrf-decdrp)*alc-(ecf-ecp)*dacdr)/alc**2
         vcl1=ce+ro*(decdrp+(dacdr*fz*bz41+alc*fz*dbdr*zta4)/fdd0)
         vcl2=2.d0*alc/(fdd0*ro)*(dfdz*bz41+4.*fz*beta*zta3)
          vclu=vcl1+vcl2*rod
          vcld=vcl1+vcl2*(-rou)
c.....
        go to 100
      endif
c.....
      if(ip9.eq.1) then
        go to 100
      endif
c.....
  100 continue

c.....gradient expansion.
c.....
      if(igd.le.0) go to 200
c       write(6,*)  '  GGA '
c.....
        gr2=agr**2
        gr2u=agru**2
        gr2d=agrd**2

        c56=5./6.
        c115=1./15.
        c1415=14./15.
        c2915=29./15.
        c2q23=2.**c23
        c83=8./3.
c.....
c.....  dsprs: divergence-suppress-factor.
c       if((log(dspr)+2.*log(agr)-c83*log(ro)).gt.8.0) go to 200
          dsprs=1.
          if(idspr.eq.1) dsprs=exp(-dspr*gr2/ro**c83)
c.....
c     agr,agru,agrd: abs(grad(rho)), for all, up, and down.
cc    gr2,gr2u,gr2d: grad(rho_all)**2, grad(rho_up)**2, grad(rho_d)**2.
c     g2r,g2ru,g2rd: laplacian rho_all, _up and _down.
c     gggru,-d: grad(rho)*grad(abs(grad(rho))) for all,up and down.
c     grgru,-d: grad(rho_all)*grad(rhor_up) and for down.

c       g2r=ddrr+2.*drr/rv
c.....
        rou53=rou**c53
c.....
c.....  edrru: d(abs(d(rou)/dr))/dr, edrrd for down.
c       edrru=ddrru
c       if(drru.lt.0.) edrru=-ddrru
c.....
c       agr,agbru,-d: abs(grad(rho)),for rou, rod.
c       gggru,-d: grad(rho)*grad(abs(grad(rho))) for up and down.
c.....  su:at ro=2*rou. 1/(2(3*pai**2)**(1/3))*|grad(rou)|/rou**(4/3).
        su=0.128278244d0*agru/rou43
        if(su.gt.huges) go to 200
c       g2ru=ddrru+2.*drru/rv
        tu=.016455307d0*g2ru/rou53
        uu=0.002110857d0*gggru/rou3

      if(ip9.ne.1) then

        su2=su*su
        su3=su*su2
        su4=su2*su2
        su6=su2*su4
c.....
        f1u=1.d0+1.296d0*su2+14.d0*su4+.2d0*su6
        f2u=2.592d0+56.d0*su2+1.2d0*su4
        f3u=112.d0*su+4.8d0*su3
c.....
c.....  fu: fgga(su) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
c.....  sidfu: su**(-1)*d(fu)/d(su)).
c.....  dsdfu: d(sidfu)/d(su).
c.....  xedgu; exchange energy density xe at ro=2*rou.(16) of p.w.
cc....      xedgu=ax*rou**(4/3)*(fu-1). ax=2**(4/3)*1.47711(ry).

        fu=f1u**c115
        sidfu=c115*f1u**(-c1415)*f2u
        dsdfu=c115*f1u**(-c2915)*(-c1415*su*f2u**2+f1u*f3u)
c.....
        xedgu=-3.722102942d0*(fu-1.d0)*rou43
c.....
        vxgu=dsprs*alf*rou13*(c43*(fu-1.)-tu*sidfu-(uu-c43*su3)*dsdfu)

      else

        dbrou=rou*2.

        call exch91(dbrou,su,uu,tu,xedlu,xedgu,vxlu,vxgu)

         xedl=xedlu/2.d0

      endif

c.....
c.....bxu,bxd,bx: grad-coeff. for exchange.

      bxu=xedgu/gr2u*rou43
c.....
        rod53=rod**c53
c       edrrd=ddrrd
c       if(drrd.lt.0.) edrrd=-ddrrd

        sd=0.128278244d0*agrd/rod43
        if(sd.gt.huges) go to 200

c       g2rd=ddrrd+2.*drrd/rv

        td=.016455307d0*g2rd/rod53
        ud=0.002110857d0*gggrd/rod3

      if(ip9.ne.1) then

        sd2=sd*sd
        sd3=sd*sd2
        sd4=sd2*sd2
        sd6=sd2*sd4
c.....
        f1d=1.d0+1.296d0*sd2+14.d0*sd4+.2d0*sd6
        f2d=2.592d0+56.d0*sd2+1.2d0*sd4
        f3d=112.d0*sd+4.8d0*sd3
c.....
c.....  fd: fgga(sd) eq.(20) of Perdew-Wang.(Phys.Rev..b33,8800,'86.)
c.....  sidfd: sd**(-1)*d(fd)/d(sd)).
c.....  dsdfd: d(sidfd)/d(sd).
c.....  xedgd; exchange energy density xe at ro=2*rod.(16) of p.w.
cc....      xedgd=ax*rod**(4/3)*(fd-1). ax=2**(4/3)*1.47711(ry).

        fd=f1d**c115
        sidfd=c115*f1d**(-c1415)*f2d
        dsdfd=c115*f1d**(-c2915)*(-c1415*sd*f2d**2+f1d*f3d)
c.....
        xedgd=-3.722102942d0*(fd-1.d0)*rod43
c.....
        vxgd=dsprs*alf*rod13*(c43*(fd-1.)-td*sidfd-(ud-c43*sd3)*dsdfd)

      else

        dbrod=rod*2.

        call exch91(dbrod,sd,ud,td,xedld,xedgd,vxld,vxgd)

        xedl=xedl+xedld/2.d0

      endif

      bxd=xedgd/gr2d*rod43
c.....

      xedg=dsprs*(xedgu+xedgd)/2.d0

      bx=(bxu+bxd)/2.d0
 
      if(iex.eq.1) go to 200

c.....
c.... cro: c(n) of (6),Phys.Rev..b33,8822('86). in ry.
c.... dcdr: d(cro)/d(ro).
c.....0.001625816=1.745*f(=0.11)*cro(rs=0).
   
      if(ip9.ne.1) then

       crr1=.005136d0+.046532d0*rs+1.4778d-5*rs2
       crr2=1.+8.723d0*rs+.472d0*rs2+.07389d0*rs3
       cro=.003334d0+crr1/crr2
       dcdr=((.046532d0+2.9556d-5*rs)*crr2-crr1*(8.723d0+.944d0*rs+
     &      .22167d0*
     &      rs2))/crr2/crr2*(-rs/ro/3.)
c.....
        fai=0.001625816d0/cro*agr/ro76
        if(fai.gt.hugef) go to 200
        fai2=fai*fai
        expfai=exp(-fai)
c.....
c.....
        if(ipg.eq.0) then

          dd=0.707106781d0*SQRT((1.+zta)**c53+(1.-zta)**c53)
c.....    ssfc: spin-scaling-factor for gradient correlation energy.
          ssfc=1./dd
          crdc=c56/(ro113*dd**2)*c2q23
          vc45u=-crdc*(rou23-rod23)*((1.-fai)*rod*gr2-(2.-fai)*
     &          ro*grgrd)
          vc45d=-crdc*(rod23-rou23)*((1.-fai)*rou*gr2-(2.-fai)*
     &          ro*grgru)

        elseif(ipg.eq.1) then

          ssfc=gz
          crdc=c2q23/(3.*gz*ro83)
          vc45u=crdc*(1./rou13-1./rod13)*((1.-fai)*rod*gr2-(2.-fai)*
     &          ro*grgrd)
          vc45d=crdc*(1./rod13-1./rou13)*((1.-fai)*rou*gr2-(2.-fai)*
     &          ro*grgru)

        elseif(ivg.eq.1) then

         write(6,'(/'' non-spher modification not completed for vg'')')
         stop30

          if(ivn.eq.0) then
            write(6,223) ivn,ivg
  223       format(/' ivn should be 1 for ivg=1. ivn,ivg=',2i5/)
            stop16
          endif
c.....
          dfdz=fdfdz(zta)
          vz=(1.+alc/ecp*fz/fdd0*bz41)**c13
c.....
          ssfc=vz
c.....
c.....    dvdru,dvdrd: d(vz)/drou,-d.
          ef3vi=1.d0/(ecp*fdd0*3.d0*vz**2)
          dvdr1=(dacdr*bz41-alc/ecp*decdrp*bz41+alc*dbdr*zta4)*
     &        fz*ef3vi
          dvdr2=2.d0*(dfdz*bz41+4.d0*fz*beta*zta3)*alc/ro2*ef3vi
          dvdru=dvdr1+dvdr2*rod
          dvdrd=dvdr1-dvdr2*rou
c.....
          vc45u=((1.-fai)*gr2*dvdru-(2.-fai)*grgrd*(dvdru-dvdrd))/
     &        (vz*ro)
          vc45d=((1.-fai)*gr2*dvdrd-(2.-fai)*grgru*(dvdrd-dvdru))/
     &          (vz*ro)

        endif

c.....  cedg: correlation-energy-density from grad.expansion.
c.....  bcr: grad-coeff. for correlation.
        bcr=ssfc*expfai*cro
        cedg=dsprs*bcr*gr2/ro43
c.....
c.....  vccf:v-correlation-coeff.
        vccf=-ssfc*expfai*cro/ro13
        vc13=(2.-fai)*g2r/ro-(c43-c113*fai+c76*fai2)*gr2/ro2+
     &    fai*(fai-3.)*gggr/agr/ro
c    &    fai*(fai-3.)*ddrr/ro
        vc6=-gr2/ro*(fai2-fai-1.)/cro*dcdr
c.....
        vcgu=dsprs*vccf*(vc13+vc6+vc45u)
c.....
        vcgd=dsprs*vccf*(vc13+vc6+vc45d)

      else

c       PW91

        call corlsd(rs,zta,ec,vclu,vcld,ecrs,eczta,alfc)

          vclu=vclu*2.
          vcld=vcld*2.
          cedl=ec*2.*ro

          fk=1.91915829d0/rs
          sk=SQRT(4.*fk/pi)
          tksg=2.*sk*gz
          tc=agr/(ro*tksg)
c           gagr: d(ABS(d(ro)/dr))/dr.       
c           gagr=ddrr
c           if(drr.lt.0.) gagr=-ddrr
          uc=gggr/(ro2*tksg**3)
c         uc=drr*gagr/(ro2*tksg**3)
          vc=g2r/(ro*tksg**2)
          wc=gzgr/(ro*tksg**2)
c         wc=drr*dzr/(ro*tksg**2)

        call cpw91(fk,sk,gz,ec,ecrs,eczta,rs,zta,tc,uc,vc,wc,
     &    cedg,vcgu,vcgd)

          vcgu=vcgu*2.
          vcgd=vcgd*2.
          cedg=cedg*ro*2.*dsprs

          bcr=cedg/gr2*ro43

      endif
c.....
  200 continue
      
      xcptu=vxlu+vclu+vxgu+vcgu
      xcptd=vxld+vcld+vxgd+vcgd
check
c     ro is small
c
      xced=0.0
      if(ro.gt.sml) xced=(xedl+cedl+xedg+cedg)/ro
        
c       write(6,'(/'' vxlu,vxld,vclu,vcld,xedl,cedl ro='',7f11.5)') vxlu,
c    &    vxld,vclu,vcld,xedl,cedl,ro
c       write(6,'(/'' vxgu,vxgd,vcgu,vcgd,xedg,cedg='',6f12.7)') vxgu,
c    &    vxgd,vcgu,vcgd,xedg,cedg

      return
      end
