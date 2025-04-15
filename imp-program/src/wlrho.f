      Subroutine wlrho(ifile,natyp,nspin,z,alat,rmt,rmtnew,rws,ititle,r,
     $                 drdi,rho2ns,rhoc,irws,a,b,txc,kxc,irns,lpot,irc,
     $                 kshape,efermi,vbc,ecore,lcore,ncore,igga)
      Implicit None
c-----------------------------------------------------------------------
c      this subroutine stores in 'ifile' the charge density
c
c        (see to subroutine start , where most of the arrays are
c         described)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      Integer natypd,nspind
      Parameter (NATYPD=38,NSPIND=2)
      Integer lpotd
      Parameter (lpotd=8)
      Integer irmd
      Parameter (irmd=1484)
      Integer lmpotd
      Parameter (lmpotd= (lpotd+1)**2)
C     ..
C     .. Scalar Arguments ..
      Double precision alat
      Integer ifile,igga,kshape,kxc,lpot,natyp,nspin
C     ..
C     .. Array Arguments ..
      Double precision a(*),b(*),drdi(irmd,*),ecore(20,*),efermi(2),
     $                 r(irmd,*),rho2ns(irmd,lmpotd,natypd,nspind),
     $                 rhoc(irmd,*),rmt(*),rmtnew(*),rws(*),vbc(2),z(*)
      Integer irc(*),irns(*),irws(*),ititle(20,*),lcore(20,*),ncore(*)
      Character*24 txc(*)
C     ..
C     .. Local Scalars ..
      Double precision a1,b1,rhobound,rmax,rmt1,rmtnw1,rv,sum,z1
      Integer i,icore,ih,ip,ir,irmin,irns1,is,isave,j,lm,lmnr,lmpot,
     $        ncore1,nr
C     ..
C     .. Local Arrays ..
      Double precision ecore1(20),vm2za(irmd),u(irmd),vmcore(irmd)
      Integer lcore1(20)
C     ..
C     .. Intrinsic Functions ..
      Intrinsic sqrt
C     ..
      isave = 1
c
      lmpot = (lpot+1)* (lpot+1)

      Do 20 ih = 1,natyp
        Do 10 is = 1,nspin
          rhobound = 1.e-16
          If ( is .eq. 2 ) rhobound = rhobound / 100.d0

          ip = nspin* (ih-1) + is
          rmt1 = rmt(ih)
          rmtnw1 = rmtnew(ih)
          z1 = z(ih)
          rmax = rws(ih)

          If (kshape.eq.0) Then
            nr = irws(ih)
          Else
            nr = irc(ih)
          End if

          irns1 = irns(ih)
          irmin = nr - irns1
          a1 = a(ih)
          b1 = b(ih)
          ncore1 = ncore(ip)
c
          Do j = 2,nr
c
c     Store the 'true' charge density (rho2ns is rho*r*r)
            vm2za(j) = rho2ns(j,1,ih,is)/ (r(j,ih)*r(j,ih))
            vmcore(j) = rhoc(j,ip)/ (r(j,ih)*r(j,ih))
          End do
c     Next is 'bad' extrapolation
          vm2za(1) = vm2za(2)
          vmcore(1) = vmcore(2)
c
          If (ncore1.ge.1) Then
            Do j = 1,ncore1
              lcore1(j) = lcore(j,ip)
              ecore1(j) = ecore(j,ip)
            End do
          End if
c
          If (igga.eq.0) Then
            Write (ifile,FMT=9000) (ititle(i,ip),i=1,7),txc(kxc+1)
          Else
            Write (ifile,FMT=9000) (ititle(i,ip),i=1,7),txc(igga+3)
          End if
          Write (ifile,FMT=9010) rmt1,alat,rmtnw1
          Write (ifile,FMT=9020) z1,rmax,efermi(is),vbc(is)
          Write (ifile,FMT=9030) nr,a1,b1,ncore1
          If (ncore1.ge.1) Write (ifile,FMT=9040) (lcore1(icore),
     $        ecore1(icore),icore=1,ncore1)
c
c---> store the full charge density, but the non spherical contribution
c     only from irns1 up to irws1 ;
c
          Write (ifile,FMT=9061) nr,irns1,lmpot,isave,is
          Write (ifile,FMT=9070) (vmcore(ir),ir=1,nr)
          If ( is .eq. 1 ) Write (ifile,FMT=9062) nr,irns1,lmpot,isave
          If ( is .eq. 2 ) Write (ifile,FMT=9063) nr,irns1,lmpot,isave
          Write (ifile,FMT=9070) (vm2za(ir),ir=1,nr)
c
          lmnr = 1
          Do lm = 2,lmpot
c            Do ir = irmin,nr
c     rho2ns is calculated everywhere, but potential is adsumed to be
c     spherically symmetric if r < r(irmin)
            Do ir = 2,nr
              u(ir) = rho2ns(ir,lm,ih,is)/(r(ir,ih)*r(ir,ih))
            End Do
c     Next is 'bad' extrapolation
            u(1) = u(2)
            sum = 0.0d0
c            Do ir = irmin,nr
            Do ir = 2,nr
              rv = u(ir)*r(ir,ih)
              sum = sum + rv*rv*drdi(ir,ih)
            End do

            If (sqrt(sum).gt.rhobound) Then
              lmnr = lmnr + 1
              Write (ifile,FMT=9060) lm
c              Write (ifile,FMT=9070) (u(ir),ir=irmin,nr)
              Write (ifile,FMT=9070) (u(ir),ir=1,nr)
            End if

          End do
c
c---> write a one to mark the end
c
          If (lmnr.lt.lmpot) Write (ifile,FMT=9060) isave

c     is-loop (spins)
   10   Continue
c     ih-loop (atoms)
   20 Continue



 9000 Format (7a4,6x,'  exc:',a24,3x,a10)
 9010 Format (3f12.8)
 9020 Format (f10.5,/,f10.5,2f15.10)
 9030 Format (i5,/,2d15.8,/,i2)
 9040 Format (i5,1p,d15.6)
 9050 Format (1p,2d15.6,1p,d15.8)
 9060 Format (10i5)
 9061 Format (4i5,1x,' core density for ispin:',i5)
 9062 Format (4i5,1x,' total charge density ')
 9063 Format (4i5,1x,' spin density ')
 9070 Format (1p,4d20.13)
      End
