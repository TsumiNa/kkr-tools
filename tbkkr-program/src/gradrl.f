      SUBROUTINE GRADRL(nspin,mesh,l1max,dx,rhol,
     &  rv,drdi,ipan,ipand5,ircut)
c.....------------------------------------------------------------------
c     gradient of rl with rl defined by charge density=sum(rl*ylm).
c     mesh,l1max: max of mesh and l+1.
c     IRMD,LMPOTD: maxima of corresponding dimension parameters.
c     drrl=d(rl)/dr, ddrrl=d(drrl)/dr, drrul=d(rl-up)/dr, 
c     ztal: zeta for each l-component necessary to get down-components.
c.....------------------------------------------------------------------
       implicit REAL*8(a-h,o-z)                                
c.....------------------------------------------------------------------
c       INTEGER IRMD,IRNSD,LPOTD,NSPIND
c       PARAMETER(IRMD=484,IRNSD=208,LPOTD=8,NSPIND=2)
       include 'inc.fi'
       INTEGER LMPOTD,ipand5
       PARAMETER (LMPOTD= (LPOTD+1)**2)
      REAL*8 rhol(IRMD,NSPIND,LMPOTD),
     &           rv(IRMD),drdi(IRMD)
c
c
c
      REAL*8  drrl(IRMD,LMPOTD),ddrrl(IRMD,LMPOTD),
     &           drrul(IRMD,LMPOTD),
     &           ddrrul(IRMD,LMPOTD)
      REAL*8 rl1(IRMD),ztal(IRMD),rl1udm(IRMD)
      INTEGER  ircut(0:ipand)
      REAL*8 drdi2(IRMD) 
       common/cgradr/drrl,ddrrl,drrul,ddrrul
c      common/cgradr/drrl(IRMD,LMPOTD),ddrrl(IRMD,LMPOTD),
c    &               drrul(IRMD,LMPOTD),ddrrul(IRMD,LMPOTD)
c    
       data zero,zero1/0.e0,1.0d-12/
c.....------------------------------------------------------------------
      pi=ACOS(-1.d0)
      s4=SQRT(4.*pi)
      llmax=l1max*l1max
      write(6,9018) l1max,mesh,nspin,ipan
 9018 format(1x,' l1max=',i5,' mesh=',i5, 'nspi=',i5,
     & ' ipan=',i5)
C     write(6,9121) rv(424),drdi(424)
C9121 format(1x,' rv drdi',5d20.10)
C     write(6,9019) (ircut(ii),ii=0,ipan)
C9019 format(1x,' ircut',10i5)
c     write(6,9120) dx,rhol(mesh,1,1)
C9120 format(1x,' dx rhol =',2d20.10)
c
      if(l1max.gt.LMPOTD) then
	write(6,'(/'' l1max.gt.LMPOTD. l1max,LMPOTD='',2i3)') 
     &  l1max,LMPOTD
	stop20
      endif

      
      do 60 ip=1,ipan
      ist=ircut(ip-1)+1
      ien=ircut(ip)
      write(6,9050) ip,ist,ien
 9050 format(1x,'  ip ist ien',3i5)
      if(ip.eq.1) then
      do 61 ir=ist,ien
      drdi2(ir)=dx
   61 continue
      else
      do 62 ir=ist,ien
      drdi2(ir)=zero
   62 continue
      end if
   60 continue
c
c
      do 10 I1=1,llmax  


        if(nspin.eq.1) go to 7777 

        do 15 ir=2,mesh
         r2=rv(ir)*rv(ir) 
c        rl1(ir)=EXP(rv(ir))
c        if(nspin.eq.1) then
c        rl1(ir)=rhol(ir,1,I1)*2.e0
c        ztal(ir)=zero
c        else
c          rl1(ir)=EXP(rv(ir))
           chgden=rhol(ir,1,I1)+rhol(ir,2,I1)
           spiden=rhol(ir,2,I1)-rhol(ir,1,I1)
           if(ABS(chgden).ge.zero1) then
           rl1(ir)=chgden
           ztal(ir)=spiden/chgden

           else
           rl1(ir)=zero 
           ztal(ir)=zero
 
          end if
   15 continue   

       go to 7778


 7777   continue
        do 16 ir=2,mesh
         r2=rv(ir)*rv(ir) 
           rl1(ir)=rhol(ir,1,I1)+rhol(ir,2,I1) 
           ztal(ir)=zero
C     check
C
C        rl1(ir)=EXP(rv(ir))
C        ztal(ir)=zero
C
C 
   16 continue   

 7778 continue  
c
C        write(6,9010) I1,rv(ir),EXP(rv(ir)),rhol(ir,1,I1),
c    &                 rhol(ir,2,I1)
c9010    format(1x, ' I1 rv exp rhol',I5,4e15.6)   
c
        rl1(1)=rl1(2)
        ztal(1)=ztal(2)
c
C       write(6,*) 'mesh',mesh
C       write(6,9000) (ir,rv(ir),rl1(ir),ztal(ir),
C    &     drdi(ir),drdi2(ir),ir=1,mesh,20)
C9000   format(1x,' ir rv rl1 ztal drdi drdi2',i5,5F15.5)
c

        do 20 ip=1,ipan
        ist=ircut(ip-1)+1
        ien=ircut(ip)
c
C      write(6,9005) ip,ist,ien,I1
C9005  format(1x,/,' ip=',i5,' ist=',i5, ' ien=',i5,' I1=',i5)
C       write(6,*) ' before gradr '
        call gradr(nspin,ist,ien,1.d0,drdi,drdi2,rl1,ztal,
     &    drrl(1,I1),ddrrl(1,I1),drrul(1,I1),ddrrul(1,I1),rl1udm)
c
       if(ip.eq.1) then
       do 21 ir=1,4
       drrl(ir,I1)=drrl(5,I1)
       ddrrl(ir,I1)=ddrrl(5,I1)
       drrul(ir,I1)=drrul(5,I1)
       ddrrul(ir,I1)=ddrrul(5,I1)
   21  continue
      end if
C
      if(nspin.eq.1) then
      do 17 ir=ist,ien
      drrul(ir,I1)=drrl(ir,I1)/2.e0
      ddrrul(ir,I1)=ddrrl(ir,I1)/2.e0
   17 continue
      end if
C

C       if(I1.eq.1.or.I1.eq.4) then
C       write(6,9001) (ir,rv(ir),rl1(ir),drrl(ir,I1),ddrrl(ir,I1),
C    &    drrul(ir,I1),ddrrul(ir,I1),ir=ist,ien,20)
C9001 format(1x,' ir rv rl1 drrl ddrrl drrul ddrrul',i5,6F12.5)
C        end if

C     stop77
   20 continue

   10 continue
C
c     iir=10
c     llm=4
c     write(6,9996) iir,llm
c9996 format(1x,' iir llm',10i5)
c     write(6,9997) drrl(iir,llm),ddrrl(iir,llm),
c    &              drrul(iir,llm),ddrrul(iir,llm)
c9997 format(1x,' drrl ddrrl drrul ddrrul',5f10.5)
C
      return
      end 
