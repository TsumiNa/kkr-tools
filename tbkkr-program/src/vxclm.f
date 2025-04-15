c 13.10.95 ***************************************************************
      SUBROUTINE VXCLM(EXC,KTE,KXC,LMAX,NSPIN,NSTART,NEND,RHO2NS,V,R,
     +                 DRDI,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM,
     +                 IMAXSH,IFUNM,THETAS,YR,WTYR,IJEND,WORK,LMSP)
      implicit none
c ************************************************************************
c     add the exchange-correlation-potential to the given potential
c     and if total energies should be calculated (kte=1) the exchange-
c     correlation-energies are calculated .
c     use as input the charge density times r**2 (rho2ns(...,1)) and
c     in the spin-polarized case (nspin=2) the spin density times r**2
c     (rho2ns(...,2)) .
c     the density times 4 pi is generated at an angular mesh .
c     the exchange-correlation potential and the exchange-correlation
c     energy are calculated at those mesh points with a subroutine .
c     in the paramagnetic case the "spin-density" is set equal zero .
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
      include 'inc.fi'
c      INTEGER NATYPD,NSPIND
c      PARAMETER (NATYPD=1,NSPIND=2)
c      INTEGER IRMD,LPOTD
c      PARAMETER (IRMD=1484,LPOTD=8)
c      INTEGER NFUND,IRID,NGSHD
c      PARAMETER (NFUND=24,IRID=435,NGSHD=3079)
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
c      INTEGER IJD
c      PARAMETER (IJD=8*LMPOTD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,KSHAPE,KTE,KXC,LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),EXC(0:LPOTD,*),GSH(*),R(IRMD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*),V(IRMD,LMPOTD,*),
c p.zahn, 10.3.99
c     +                 WORK(IRMD,LMPOTD),WTYR(IJD,*),YR(IJD,*)
     +                 WORK,WTYR(IJD,*),YR(IJD,*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Scalars in Common ..
      LOGICAL strco1
      CHARACTER*5 LMACH
C     ..
C     .. Common blocks ..
      COMMON /CMACH/strco1,LMACH
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3
      INTEGER IATYP,ICELL,IFUN,IJ,IN,IPAN1,IPOT,IQ,IQ1,IR,IRC1,IRH,IRS1,
     +        IS,ISPIN,J,L,LM,LM2,LMMAX,M,WLEN
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJD),
     +                 FPRHO(IJD,2),VXC(IJD,2),VXCR(2:3,NSPIND)
      INTEGER WLEN1(140),WLEN2(140)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DDOT
      INTEGER MAPBLOCK,MYNODE
C     ,NUMNODES
      EXTERNAL DDOT,MAPBLOCK,MYNODE
C     ,NUMNODES
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,GCOLX,GDSUM,SIMP3,SIMPK,VOSKO,VXCSPO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DMAX1
C     ..
      FPI = 16.0D0*DATAN(1.0D0)
      LMMAX = (LMAX+1)* (LMAX+1)
c
c---> loop over given representive atoms
c
      DO 270 IATYP = NSTART,NEND

        IF (KSHAPE.NE.0) THEN
          IPAN1 = IPAN(IATYP)
          ICELL = NTCELL(IATYP)
          IRC1 = IRCUT(IPAN1,IATYP)
          IRS1 = IRCUT(1,IATYP)
        ELSE
          IRC1 = IRWS(IATYP)
          IRS1 = IRC1
        END IF

        DO 10 ISPIN = 1,NSPIN
          VXCR(2,ISPIN) = 0.0D0
          VXCR(3,ISPIN) = 0.0D0
   10   CONTINUE
c
c--->   initialize for ex.-cor. energy
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
c
c--->   loop over radial mesh
c
        IF (LMACH.EQ.'INTEL') THEN
          DO 60 IN = 1,140
            WLEN = 0
            WLEN2(IN) = 4
   60     CONTINUE
        END IF

        DO 160 IR = 2,IRC1

          IF (LMACH.EQ.'INTEL') THEN
            IN = MYNODE()
c           IF (IN.NE.MAPBLOCK(IR,2,IRC1,1,0,NUMNODES()-1)) THEN
c             GO TO 160

c           ELSE
c             WLEN = WLEN + 8
c           END IF

          END IF

c
c--->     generate the densities on an angular mesh
c
          DO 80 IS = 1,2
            DO 70 IJ = 1,IJEND
              FPRHO(IJ,IS) = 0.D0
   70       CONTINUE
   80     CONTINUE

          FPIPR2 = FPI/R(IR,IATYP)**2
          DO 110 ISPIN = 1,NSPIN
            DO 90 LM = 1,LMMAX
              CALL DAXPY(IJEND,RHO2NS(IR,LM,IATYP,ISPIN)*FPIPR2,
     +                   YR(1,LM),1,FPRHO(1,ISPIN),1)
   90       CONTINUE
C
C--->       remove negative values of small densities
C
c            DO 100 IJ = 1,IJEND
c              FPRHO(IJ,ISPIN) = DMAX1(FPRHO(IJ,ISPIN),1.D-12)
c  100       CONTINUE
  110     CONTINUE
c
c--->     calculate the ex.-cor. potential
c
          IF (KXC.LE.1) THEN
            CALL VXCSPO(EXCIJ,FPRHO,VXC,KXC,IJEND,IJD)
          ELSE
            CALL VOSKO(EXCIJ,FPRHO,VXC,IJEND,IJD)
          END IF
c
c--->     expand the ex.-cor. potential into spherical harmonics ,
c         using the orthogonality
c
          DO 130 ISPIN = 1,NSPIN
c
c--->       determine the corresponding potential number
c
            IPOT = NSPIN* (IATYP-1) + ISPIN

            DO 120 LM = 1,LMMAX
              VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
              V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC
c
c--->         store the ex.-c. potential of 
c             ir=2 and =3 for the extrapolation
c
              IF (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR,
     +            ISPIN) = VLMXC
  120       CONTINUE
  130     CONTINUE
c
c--->     file er in case of total energies
c
          IF (KTE.EQ.1) THEN
c
c--->     expand ex.-cor. energy into spherical harmonics
c         using the orthogonality
c
            DO 150 L = 0,LMAX
              DO 140 M = -L,L
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

  140         CONTINUE

  150       CONTINUE

          END IF

  160   CONTINUE

        IF (LMACH.EQ.'INTEL') THEN
          STOP 'LMACH = ''INTEL'' in VXCLM'
c         p. zahn, 10.3.99
c          CALL GCOLX(WLEN,WLEN2,WLEN1)
c          IN = MYNODE()
c          IQ = 2
c          DO 170 IQ1 = 1,IN
c            IQ = IQ + WLEN1(IQ1)/8
c  170     CONTINUE
c          WRITE (6,FMT=*) IRC1,IN,IQ,WLEN1(IN+1)
c          DO 200 ISPIN = 1,NSPIN
c            IPOT = NSPIN* (IATYP-1) + ISPIN
c            DO 190 LM = 1,LMMAX
c              DO 180 IR = 2,IRC1
c                WORK(IR,LM) = V(IR,LM,IPOT)
c  180         CONTINUE
c              CALL GCOLX(WORK(IQ,LM),WLEN1,V(2,LM,IPOT))
c  190       CONTINUE
c  200     CONTINUE
c          CALL GDSUM(ESTOR,IRMD*LMPOTD,WORK)
c          CALL GDSUM(ER,IRMD* (LPOTD+1),WORK)
c          CALL GDSUM(VXCR,2*NSPIND,WORK)
        END IF                      ! (LMACH.EQ.'INTEL')

c
c--->   integrate er in case of total energies to get exc
c
        IF (KTE.EQ.1) THEN

          IF (KSHAPE.EQ.0) THEN

            DO 210 L = 0,LMAX
              CALL SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI(1,IATYP))
  210       CONTINUE

          ELSE                      ! (KSHAPE.EQ.0)

            DO 250 L = 0,LMAX
              DO 240 M = -L,L
                LM = L*L + L + M + 1
c
c--->           convolute with shape function
c
                DO 230 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
                  LM2 = ILM(J,2)
                  IF (LMSP(ICELL,ILM(J,3)).GT.0) THEN
                    IFUN = IFUNM(ICELL,ILM(J,3))
                    DO 220 IR = IRS1 + 1,IRC1
                      IRH = IR - IRS1
                      ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,IATYP,1)*
     +                           GSH(J)*THETAS(IRH,IFUN,ICELL)*
     +                           ESTOR(IR,LM2)
  220               CONTINUE
                  END IF
  230           CONTINUE
  240         CONTINUE
              CALL SIMPK(ER(1,L),EXC(L,IATYP),IPAN1,IRCUT(0,IATYP),
     +                   DRDI(1,IATYP))
  250       CONTINUE

          END IF                    ! (KSHAPE.EQ.0)

        END IF                      ! (KTE.EQ.1)
c
c--->   extrapolate ex.-cor potential to the origin only for lm=1
c
        DO 260 ISPIN = 1,NSPIN
          IPOT = NSPIN* (IATYP-1) + ISPIN
c
          VXC2 = VXCR(2,ISPIN)
          VXC3 = VXCR(3,ISPIN)
          VXC1 = VXC2 - R(2,IATYP)* (VXC3-VXC2)/ (R(3,IATYP)-R(2,IATYP))
c
          V(1,1,IPOT) = V(1,1,IPOT) + VXC1
  260   CONTINUE

  270 CONTINUE                      ! IATYP = NSTART,NEND
c
      RETURN
c
      END
