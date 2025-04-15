      SUBROUTINE COREL(IPF,ITC,NITMAX,KHYP,KCORE,IPR,IRNUMX,IP,IRM,
     +                 NSPIN,RHOC,NATREF,IMT,KSHAPE,ECORE,LCORE,NCORE,
     +                 CFG)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     this subroutine is responsible for relaxation of core-charge dens
c     it calls int=integration of scalarrelativistic eqation and find o
c     the appropriate eigenvalues
c     as long we don't have relativistic greensfunctions use kcore inst
c     kcore and also the statements kcor=kcore,if(ip.eq.1)kcor=2,just t
c     sure we calculate the host relativistically.
c     lmxc = lmaxcore = (0,1,2,...), argon core : lmxc = 1, krypton : l
c     kfg = configuration of core states f.e. argon core: 3300=3s,3p,0d
c                                           krypton core: 4430=4s,4p,3d
c                                             xenon core: 5540=5s,5p,4d
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IRMD
      PARAMETER (irmd=1484)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IP,IPF,IPR,IRM,IRNUMX,ITC,KCORE,KHYP,KSHAPE,NATREF,NITMAX,
     +        NSPIN
C     ..
C     .. Array Arguments ..
      REAL*8 ECORE(20,NPOTD),RHOC(IRMD,*)
      INTEGER CFG(4,NATYPD),IMT(*),LCORE(20,NPOTD),NCORE(NPOTD)
C     ..
C     .. Scalars in Common ..
      INTEGER LMAXP1
C     ..
C     .. Arrays in Common ..
      REAL*8 A1(NATYPD),B1(NATYPD),DRDI(IRMD,NATYPD),
     +                 HYPSUM(10,NPOTD),R(IRMD,NATYPD),RHYPF(150,NPOTD),
     +                 RWS(NATYPD),RWSM1(NATYPD),V(IRMD,NPOTD),
     +                 Z1(NATYPD)
      INTEGER IRT(NATYPD),IRWS(NATYPD),KFG(4,NATYPD),LMXC(NATYPD),
     +        NCMAX(NATYPD)
C     ..
C     .. Local Scalars ..
      REAL*8 A,B,DIFF,E,E1,E2,EDIFF,EI,FOURPI,GAUSS,QC,R0,
     +                 R2RHO1,R2RHO2,RLP1,RMAX,RR,SLOPE,SUM,TOL,VALUE,Z,
     +                 ZERO
      INTEGER I,IARRAY,ID,IDD,IIR,IK,IL,IN,INUC,INUCP1,IPOT,IR,IRIPE,
     +        IRIPST,IRM2,IRR,IS,ISY,IU,IUU,K,KCOR,L,LMP1,LP1,NC,NMAX,
     +        NN,NONSRA,NR,NREM
      LOGICAL VLNC
C     ..
C     .. Local Arrays ..
      REAL*8 DRH(4),F(IRMD,2),G(IRMD,2),RHO(IRMD,2),WGT(2)
      INTEGER NRE(2)
      CHARACTER*4 SPN(2),TEXT(5)
      CHARACTER*8 TEXTAT(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL BREIT,INTCOR,SIMP3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX,REAL
C     ..
C     .. Common blocks ..
      COMMON /DLOG87/V,RWS,RWSM1,DRDI,R,Z1,A1,B1,IRWS,IRT,LMAXP1
      COMMON /NUCDAT/RHYPF,HYPSUM,KFG,LMXC,NCMAX
C     ..
C     .. Save statement ..
      SAVE SPN,TEXT,TEXTAT,ZERO,FOURPI,GAUSS,/DLOG87/,/NUCDAT/
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','s   ','p   ','d   ','f   ','g   '/
      DATA TEXTAT/'    host','impurity'/
      DATA ZERO,FOURPI,GAUSS/0.0D0,12.56637D0,524.D0/
C     ..
c
      ISY = 1
      IF (IP.GT.NATREF) ISY = 2
      IK = IP - (ISY-1)*NATREF
      VLNC = .false.
      VALUE = 1.D-8
      SLOPE = -1.D-8
      E2 = 50.0D0
      IF (KCORE.LT.3) THEN
        KCOR = KCORE
        IF (IP.LE.NATREF) KCOR = 2

      ELSE
        KCOR = 1
      END IF

      NONSRA = KCOR - 1
c
      ID = NSPIN* (IP-1) + 1
      IU = ID + NSPIN - 1
      IDD = 1
      IUU = NSPIN
c
      IF (ITC.EQ.1) THEN
        DO 10 IL = 1,4
          KFG(IL,IP) = 0
          CFG(IL,IP) = IL
   10   CONTINUE
        DO 20 IL = 1,NCORE(IU)
          IF (LCORE(IL,IU).EQ.0) KFG(1,IP) = KFG(1,IP) + 1
          IF (LCORE(IL,IU).EQ.1) KFG(2,IP) = KFG(2,IP) + 1
          IF (LCORE(IL,IU).EQ.2) KFG(3,IP) = KFG(3,IP) + 1
          IF (LCORE(IL,IU).EQ.3) KFG(4,IP) = KFG(4,IP) + 1
   20   CONTINUE
        IF (KFG(2,IP).NE.0) KFG(2,IP) = KFG(2,IP) + 1
        IF (KFG(3,IP).NE.0) KFG(3,IP) = KFG(3,IP) + 2
        IF (KFG(4,IP).NE.0) KFG(4,IP) = KFG(4,IP) + 3
        IF (KFG(2,IP).NE.0) LMXC(IP) = 1
        IF (KFG(3,IP).NE.0) LMXC(IP) = 2
        IF (KFG(4,IP).NE.0) LMXC(IP) = 3
c       DO 30 IL = 1,NCORE(IU)
        DO 30 IL = 1,4
          CFG(IL,IP) = MAX(CFG(IL,IP),KFG(IL,IP)+1)
   30   CONTINUE
        WRITE (6,FMT=9000) LMXC(IP), (KFG(IL,IP),IL=1,4)
        WRITE (6,FMT=9000) LMXC(IP), (CFG(IL,IP),IL=1,4)
      END IF
c
c---> attention renormalization confined to muffin tin sphere in case
c---> of shape-corrected calculation
      IF (KSHAPE.NE.0) THEN
        NR = IMT(IP)
        RMAX = R(NR,IP)
        IRM2 = NR

      ELSE


        NR = IRWS(IP)
        RMAX = RWS(IP)
        IRM2 = IRM
      END IF

      A = A1(IP)
      B = B1(IP)
      Z = Z1(IP)
      TOL = 1.0D-12* (Z*Z+1.D0)
      LMP1 = LMXC(IP) + 1
      NC = 0
      INUC = -IRNUMX
c
      DO 40 IR = 1,IRM2
        RHOC(IR,ID) = ZERO
        RHOC(IR,IU) = ZERO
        RHO(IR,IDD) = ZERO
        RHO(IR,IUU) = ZERO
   40 CONTINUE
      DO 50 IR = 1,IRNUMX
        HYPSUM(IR,ID) = ZERO
        HYPSUM(IR,IU) = ZERO
        RHYPF(IR,ID) = ZERO
        RHYPF(IR,IU) = ZERO
   50 CONTINUE
c
      IF (IPR.EQ.0) WRITE (IPF,FMT=9010) IK,TEXTAT(ISY)
c
c    begin loop over l
c
      DO 160 LP1 = 1,LMP1
        L = LP1 - 1
        RLP1 = REAL(LP1)
        E1 = (-5.D0- ((Z+1.D0)/RLP1)**2)*1.5D0 - 50.D0
        NMAX = KFG(LP1,IP)
        IF (NMAX.NE.0) THEN
c
c    core states
c
          DO 150 IN = LP1,NMAX
            NN = IN - LP1
            NC = NC + 1
            INUC = INUC + IRNUMX
            DO 60 IS = 1,NSPIN
              I = NSPIN* (IP-1) + IS
              E = ECORE(NC,I)
              EI = ECORE(NC,I)
              IF (IPR.NE.0) WRITE (IPF,FMT=9020) IN,TEXT(LP1),NN,
     +            SPN(IS),IK,TEXTAT(ISY),E
              CALL INTCOR(E1,E2,RHO(1,IS),G(1,IS),F(1,IS),V(1,I),VALUE,
     +                    SLOPE,L,NN,E,SUM,NRE(IS),VLNC,A,B,Z,RMAX,NR,
     +                    TOL,IRM2,IPR,NITMAX,NONSRA)
              EDIFF = E - EI
              ECORE(NC,I) = E
              WGT(IS) = REAL(L+L+1)/SUM
              IF (IPR.NE.0) WRITE (IPF,FMT=9030) EI,EDIFF,E
              If (ipr .ne. 0 ) Then
                Write (6,fmt='(4x,a,2d14.6)') 'Density at boundary',
     $                rho(nr,is)*wgt(is)
              End If
   60       CONTINUE
c
            IF (NSPIN.EQ.1) THEN
              NREM = NRE(IUU)

            ELSE
              NREM = MAX(NRE(IUU),NRE(IDD))
            END IF
c
            IF (KHYP.NE.0) THEN
c
              INUCP1 = INUC + 1
              RHYPF(INUCP1,ID) = ZERO
              RHYPF(INUCP1,IU) = ZERO
c
              DO 70 IR = 2,IRNUMX
                IRR = IR + INUC
                RR = R(IR,IP)
                RR = FOURPI*RR*RR
                R2RHO1 = WGT(IDD)*RHO(IR,IDD)
                R2RHO2 = WGT(IUU)*RHO(IR,IUU)
                SUM = (R2RHO1+R2RHO2)/RR
c
c                    diff=gauss*(rho2ns2rho2nso1)/rr
c
                CALL BREIT(DIFF,NONSRA,IR,IUU,NREM,L,FOURPI,GAUSS,
     +                     R2RHO1,R2RHO2,RR,WGT,R(1,IP),DRDI(1,IP),G,F)
                RHYPF(IRR,IU) = DIFF
                RHYPF(IRR,ID) = SUM
   70         CONTINUE
c
              IRIPST = 2
              IRIPE = 4
              DO 120 IS = 1,NSPIN
                I = NSPIN* (IP-1) + IS
                R0 = R(1,IP)
                DO 80 IR = 1,IRIPE
                  IRR = INUC + IR
                  DRH(IR) = RHYPF(IRR,I)
   80           CONTINUE
                DO 100 K = 1,2
                  DO 90 IR = IRIPST,IRIPE - K
                    IIR = IRIPST + IRIPE - IR
                    DRH(IIR) = (DRH(IIR)-DRH(IIR-1))/
     +                         (R(IIR,IP)-R(IIR-K,IP))
   90             CONTINUE
  100           CONTINUE
                IR = 1
                IRR = INUC + IR
                RHYPF(IRR,I) = DRH(IRIPE)
                DO 110 IR = IRIPST,IRIPE - 1
                  IIR = IRIPST + IRIPE - IR - 1
                  RHYPF(IRR,I) = RHYPF(IRR,I)* (R0-R(IIR,IP)) + DRH(IIR)
  110           CONTINUE
  120         CONTINUE
              IF (L.EQ.0) THEN
                DO 130 IR = 1,IRNUMX
                  IRR = IR + INUC
                  HYPSUM(IR,ID) = HYPSUM(IR,ID) + RHYPF(IRR,ID)
                  IF (NSPIN.EQ.2) HYPSUM(IR,IU) = HYPSUM(IR,IU) +
     +                RHYPF(IRR,IU)
  130           CONTINUE
              END IF

            END IF
c
c---> sum up contributions to total core charge
c
            DO 140 IR = 2,NREM
              RHOC(IR,ID) = RHOC(IR,ID) + RHO(IR,IDD)*WGT(IDD)
              RHOC(IR,IU) = RHOC(IR,IU) + RHO(IR,IUU)*WGT(IUU)
              RHO(IR,IDD) = ZERO
              RHO(IR,IUU) = ZERO
  140       CONTINUE
  150     CONTINUE
        END IF

  160 CONTINUE

      NCMAX(IP) = NC
      IARRAY = NC*IRNUMX
      IF (IARRAY.GT.150 .OR. IRNUMX.GT.10) WRITE (IPF,FMT=9040)
      IF (IARRAY.GT.150 .OR. IRNUMX.GT.10) THEN
        STOP 85

      ELSE


        DO 170 IS = 1,NSPIN
          IPOT = NSPIN* (IP-1) + IS
c
c---> integrate core density to get core charge
c
          CALL SIMP3(RHOC(1,IPOT),QC,1,NR,DRDI(1,IP))

          WRITE (IPF,FMT=9050) Z1(IP),QC
  170   CONTINUE
        Write (6,fmt='(4x,a,2d14.6)') 'Core density at boundary',
     $       (rhoc(nr,nspin*(ip-1)+is), is=1,nspin)

      END IF



c

 9000 FORMAT (1x,i1,1x,4i1)
 9010 FORMAT (1x,5 ('*'),' core-relaxation for ',i3,'th ',a8,'-cell',
     +       '  spin=',a4,'   l = ',a4,' was done ',5 ('*'))
 9020 FORMAT (1x,90 ('*'),/,'  n = ',i1,'  l = ',a4,'   nnode = ',i1,
     +       '  spin=',a4,i5,'th ',a8,'-cell','    einput = ',1p,d16.8)
 9030 FORMAT (1x,'  einput =',1p,d16.8,'   eout - ein =',1p,d16.8,
     +       '   eoutput = ',1p,d16.8)
 9040 FORMAT (1x,
     +' space of arrays rhypf or hypsum to small stop in
     +  subroutine corel')
 9050 FORMAT (1x,/,4x,'nuclear charge  ',f10.6,9x,'core charge =   ',
     +       f10.6)
      END
