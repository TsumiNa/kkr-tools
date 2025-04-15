      SUBROUTINE VXCGGA(EXC,KTE,KXC,LMAX,NSPIN,NSTART,NEND,RHO2NS,V,R,
     +                  DRDI,A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM,
     +                  IMAXSH,IFUNM,THETAS,YR,WTYR,IJEND,WORK,LMSP)
      implicit none
c-----------------------------------------------------------------------
c     add the exchange-correlation-potential to the given potential
c     and if total energies should be calculated (kte=1) the exchange-
c     correlation-energies are calculated .
c     use as input the charge density times r**2 (rho2ns(...,1)) and
c     in the spin-polarized case (nspin=2) the spin density times r**2
c     (rho2ns(...,2)) .
c     the density times 4 pi is generated at an angular mesh .
c     the exchange-correlation potential and the exchange-correlation
c     energy are calculated at those mesh points with a subroutine .
c     in the non-spin-polarized case the "spin-density" is
c     set equal zero .
c     after that the exchange-correlation potential and in the case of
c     total energies (kte=1) the exchange-correlation energy are
c     expanded into spherical harmonics .
c     the ex.-cor. potential is added to the given potential .
c     the expansion into spherical harmonics uses the orthogonality
c     of these harmonics . - therefore a gauss-legendre integration
c     for "theta" and a gauss-tschebyscheff integration for "phi"
c     is used .
c     all needed values for the angular mesh and angular integration
c     are generate in the subroutine sphere .
c
c     the ex.-cor. potential is extrapolated to the origin only
c     for the lm=1 value .
c
c                               b.drittler   june 1987
c
c     modified for shape functions
c                                       b. drittler oct. 1989
c     simplified and modified for Paragon X/PS
c                                       R. Zeller Nov. 1993
c-----------------------------------------------------------------------
C     .. Parameters ..
c     INTEGER NATYPD,NSPIND
c     PARAMETER (natypd=1,NSPIND=2)
c     INTEGER IRMD,LPOTD
c     PARAMETER (IRMD=484,LPOTD=8)
c     INTEGER NFUND,IRID,NGSHD
c     PARAMETER (NFUND=131,IRID=635,NGSHD=13079)
c     INTEGER IPAND
c     PARAMETER (IPAND=10)
      include 'inc.fi'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
c     INTEGER IJD
c     PARAMETER (IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,KSHAPE,KTE,KXC,LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      REAL*8 DRDI(IRMD,*),EXC(0:LPOTD,*),GSH(*),R(IRMD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*),V(IRMD,LMPOTD,*),
     +                 WORK,WTYR(IJD,*),YR(IJD,*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
C     INTEGER IGGA
C     ..
C     .. Scalars in Common ..
      LOGICAL strco1
      CHARACTER*5 LMACH
C     ..
C     .. Common blocks ..
      COMMON /CMACH/strco1,LMACH
C     ..
C     .. Local Scalars ..
      REAL*8 ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3
      INTEGER IATYP,ICELL,IFUN,IJ,IPAN1,IPOT,IR,IRC1,IRH,IRS1,IS,ISPIN,
     +        J,L,LM,LM2,LMMAX,M
C     ..
C     .. Local Arrays ..
      REAL*8 ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJD),
     +                 FPRHO(IJD,2),VXC(IJD,2),VXCR(2:3,NSPIND),
     +                 COSX(IJD),FAI(IJD)
C     ..
      REAL*8 A(NATYPD),DX,R1,R2,CHGden,SPIden,rpoint,
     &          rhol(irmd,nspind,lmpotd),rholm(lmpotd,nspind)
      INTEGER MESH,l1max,nspin2

C     .. External Functions ..
      REAL*8 DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,SIMP3,SIMPK,gradrl,mkxcpe
C     ..
      REAL*8 ST1,ET1
      REAL*8 DCLOCK
      EXTERNAL DCLOCK


C     .. Intrinsic Functions ..
      INTRINSIC DATAN,MAX1
C     ..
      REAL*8 zero,zero1
      data zero,zero1/0.d0,1.d-12/

      write(6,*) ' GGA CALCULATION '
      FPI = 16.0D0*DATAN(1.0D0)
      LMMAX = (LMAX+1)* (LMAX+1)
c
c---> loop over given representive atoms
c
      DO 220 IATYP = NSTART,NEND
        IF (KSHAPE.NE.0) THEN
          IPAN1 = IPAN(IATYP)
          ICELL = NTCELL(IATYP)
          IRC1 = IRCUT(IPAN1,IATYP)
          IRS1 = IRCUT(1,IATYP)
c        write(6,*) 'in vxc',ipan1,icell,irc1,irs1
        ELSE

          IRC1 = IRWS(IATYP)
          IRS1 = IRC1
          ipan1=1   
        END IF

        write(6,*) 'ipan1=',ipan1
        DO 10 ISPIN = 1,NSPIN
          VXCR(2,ISPIN) = 0.0D0
          VXCR(3,ISPIN) = 0.0D0
   10   CONTINUE
c
c---> initialize for ex.-cor. energy
c
        IF (KTE.EQ.1) THEN
          DO 30 L = 0,LMAX
            EXC(L,IATYP) = 0.0D0
            DO 20 IR = 1,IRC1
              ER(IR,L) = 0.0D0
   20       CONTINUE
   30     CONTINUE
c
          DO 50 LM = 1,LMMAX
            DO 40 IR = 1,IRC1
              ESTOR(IR,LM) = 0.0D0
   40       CONTINUE
   50     CONTINUE
        END IF
C
        l1max=lmax+1
        mesh=irws(iatyp)
        dx=a(iatyp)
C
        if(nspin.eq.2) then
        do 41 lm=1,lmmax
        do 42 ir=2,mesh
        r1=r(ir,iatyp) 
        r2=r1*r1
        chgden=rho2ns(ir,lm,iatyp,1)/r2
        spiden=rho2ns(ir,lm,iatyp,2)/r2
        if(ABS(chgden).le.zero1) chgden=zero
        if(ABS(spiden).le.zero1) spiden=zero
        rhol(ir,2,lm)=(chgden+spiden)/2.d0
        rhol(ir,1,lm)=(chgden-spiden)/2.d0
   42   continue
c
c       extrapolate
c
        rhol(1,1,lm)=rhol(2,1,lm)
        rhol(1,2,lm)=rhol(2,2,lm)
   41   continue
C
        else
C
        Do 43 lm=1,lmmax
        Do 44 ir=2,mesh
        r1=r(ir,iatyp)
        r2=r1*r1
c
        chgden=rho2ns(ir,lm,iatyp,1)/r2
        if(ABS(chgden).le.zero1) chgden=zero
        rhol(ir,1,lm)=chgden/2.d0
        rhol(ir,2,lm)=chgden/2.d0
   44   continue
c
c       extrapolate
        rhol(1,1,lm)=rhol(2,1,lm)
        rhol(1,2,lm)=rhol(2,2,lm)
   43   continue
        end if 
c
        ST1= DCLOCK()
C
        write(6,*) ' before gradrl'
        write(6,9111) iatyp,(ircut(ir,iatyp),ir=0,4)
 9111   format(1x,' iatyp ircut',10i5)
c       write(6,9112) dx,rhol(424,1,1),r(mesh,iatyp)
c    &             ,drdi(424,iatyp)
c9112   format(1x,' dx rholi r drdi=',4d15.10)
        write(6,*) ' iapn1 ipand',ipan1,ipand

        call gradrl(nspin,mesh,l1max,dx,rhol,r(1,iatyp),
     &      drdi(1,iatyp),ipan1,ipand,ircut(0,iatyp))
        write(6,*) ' after gradrl '
C
        ET1= DCLOCK()
        WRITE(6,*) ' TIME IN GRADRL',ET1-ST1 


c
c---> loop over radial mesh
c
        ST1= DCLOCK()

        DO 150 IR = 2,IRC1
        rpoint=r(ir,iatyp)

c
c---> generate the densities on an angular mesh
c
c         DO 70 IS = 1,2
c           DO 60 IJ = 1,IJEND
c             FPRHO(IJ,IS) = 0.D0
c  60       CONTINUE
c  70     CONTINUE
c
c          FPIPR2 = FPI/R(IR,IATYP)**2
c          DO 100 ISPIN = 1,NSPIN
c            DO 80 LM = 1,LMMAX
c              CALL DAXPY(IJEND,RHO2NS(IR,LM,IATYP,ISPIN)*FPIPR2,
c     +                   YR(1,LM),1,FPRHO(1,ISPIN),1)
c   80       CONTINUE
C
C---> REMOVE NEGATIVE VALUES OF SMALL DENSITIES
C
C            DO 90 IJ = 1,IJEND
C              FPRHO(IJ,ISPIN) = DMAX1(FPRHO(IJ,ISPIN),1.D-12)
C   90       CONTINUE
C  100     CONTINUE
c
c---> calculate the ex.-cor. potential
c
         nspin2=2

         do 71 ispin=1,nspin2
         do 73 lm=1,lmmax
         rholm(lm,ispin)=rhol(ir,ispin,lm)
   73    continue
   71    continue
c
c    only for spin-polarized
c
         write(6,*) ' before mkxcpe '
         call mkxcpe(nspin2,ir,ijend,l1max,rpoint,rholm,
     &   vxc,excij)
c
c
c
c
c---> expand the ex.-cor. potential into spherical harmonics ,
c       using the orthogonality
c   
          DO 120 ISPIN = 1,NSPIN
c
c---> determine the corresponding potential number
c
            IPOT = NSPIN* (IATYP-1) + ISPIN
            DO 110 LM = 1,LMMAX
              VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
              V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC
c
c---> store the ex.-c. potential of ir=2 and =3 for the extrapolation
c
              IF (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,
     +            ISPIN) = VLMXC
  110       CONTINUE
  120     CONTINUE
c
c---> file er in case of total energies
c
          IF (KTE.EQ.1) THEN
c
c---> expand ex.-cor. energy into spherical harmonics
c       using the orthogonality
c
            DO 140 L = 0,LMAX
              DO 130 M = -L,L
                LM = L*L + L + M + 1
                ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)
c
c---> multiply the lm-component of the ex.-cor. energy with the same
c     lm-component of the charge density times r**2 and sum over lm
c     this corresponds to a integration over the angular .
c
                IF ((KSHAPE.NE.0) .AND. (IR.GT.IRS1)) THEN
                  ESTOR(IR,LM) = ELMXC

                ELSE

                  ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,IATYP,1)*ELMXC
                END IF

  130         CONTINUE

  140       CONTINUE

          END IF

  150   CONTINUE
         ET1= DCLOCK()
         WRITE(6,*) ' TIME IN VXC ', ET1-ST1
c
c---> integrate er in case of total energies to get exc
c
        IF (KTE.EQ.1) THEN
          IF (KSHAPE.EQ.0) THEN
            DO 160 L = 0,LMAX
              CALL SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI(1,IATYP))
  160       CONTINUE

          ELSE

            DO 200 L = 0,LMAX
              DO 190 M = -L,L
                LM = L*L + L + M + 1
c
c---> convolute with shape function
c
                DO 180 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
                  LM2 = ILM(J,2)
                  IF (LMSP(ICELL,ILM(J,3)).GT.0) THEN
                    IFUN = IFUNM(ICELL,ILM(J,3))
                    DO 170 IR = IRS1 + 1,IRC1
                      IRH = IR - IRS1
                      ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,IATYP,1)*
     +                           GSH(J)*THETAS(IRH,IFUN,ICELL)*
     +                           ESTOR(IR,LM2)
  170               CONTINUE
                  END IF
  180           CONTINUE
  190         CONTINUE
              CALL SIMPK(ER(1,L),EXC(L,IATYP),IPAN1,IRCUT(0,IATYP),
     +                   DRDI(1,IATYP))
  200       CONTINUE
          END IF

        END IF
c
c---> extrapolate ex.-cor potential to the origin only for lm=1
c
        DO 210 ISPIN = 1,NSPIN
          IPOT = NSPIN* (IATYP-1) + ISPIN
c
          VXC2 = VXCR(2,ISPIN)
          VXC3 = VXCR(3,ISPIN)
          VXC1 = VXC2 - R(2,IATYP)* (VXC3-VXC2)/ (R(3,IATYP)-R(2,IATYP))
c
          V(1,1,IPOT) = V(1,1,IPOT) + VXC1
  210   CONTINUE
  220 CONTINUE
c
      END
