      SUBROUTINE CPLXWF(IE,E,KVREL,NVSP,ISPIN,NATREF,ALPHA,NSPIN,NT,
     +                  IPAN,IRCUT,MASS,C,DROR,RS,S,PZ,FZ,QZ,SZ,TMAT,
     +                  VM2Z,RWS,RWSM1,DRDI,R,Z,A,B,IRWS,IRT,LMAXP1,
     $                  irmd)
      Implicit None
c-----------------------------------------------------------------------
c  on output: the common block crdfun contains the regular radial wave-
c             functions (in array pz) and the non-regular ones corres-
c             ponding to the hankel-functions in free space (in array qz
c             it also contains the t-matrices in array tmat.
c
c             the generalized phase shifts can be calculated by using
c             a wronski relation :
c
c                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0
c
c             where hl is the free hankel function and rl the regular
c             solution . using the analytical behaviour of rl at the
c             origin : rl = alphal * r**(l+1)  ; r->0
c             therefore the generalized phase shifts can be calcu-
c             lated directly with the renormalization alphal .
c
c
c     attention:    this version is slightly modified .
c                   to get proper t - matrices the t - matrix has to
c                   be determined at the mt radius in case of mt
c                   calculation and at the ws sphere in case of ws
c                   calculation . the array irt contains irmt in
c                   case of mt calculation and irws in case of ws
c                   calculation .
c
c                   the inwards integration starts now at irt1
c                   since the discontinuity of the potential at
c                   irws or irmt causes problems for the inwards
c                   integration .
c                                           b.drittler nov.1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NSPIND=2)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
c      PARAMETER (irmd=1484,lmaxd=4,LMX=LMAXD+1)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 E
      REAL*8 C
      INTEGER IE,ISPIN,KVREL,LMAXP1,NATREF,NSPIN,NT,NVSP
      Integer irmd
C     ..
C     .. Array Arguments ..
      COMPLEX*16 ALPHA(0:LMAXD,*),FZ(IRMD,0:LMAXD,NATYPD),
     +     MASS(IRMD,NATYPD),PZ(IRMD,0:LMAXD,NATYPD),
     +     QZ(IRMD,0:LMAXD,NATYPD),SZ(IRMD,0:LMAXD,NATYPD),
     +     TMAT(0:LMAXD,NATYPD)
      REAL*8 A(NATYPD),B(NATYPD),DRDI(IRMD,NATYPD),
     +     DROR(IRMD,NATYPD),R(IRMD,NATYPD),
     +     RS(IRMD,0:LMAXD,NATYPD),RWS(NATYPD),
     +     RWSM1(NATYPD),S(0:LMAXD,NATYPD),VM2Z(IRMD,NPOTD),
     +                 Z(NATYPD)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRT(NATYPD),IRWS(NATYPD)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 ALPHAL,ARG,BL,CZERO,EK,EKLFAC,FAC,HL,PN,QF,SLOPE,
     +     TL,TLSQEZ,VALUE,W,X,Y
      REAL*8 R1,RIRC,RSIRC,S1
      INTEGER I,I1,IH,IPE,IRC1,IRWS1,KVSRA,L,LMAX,N,NEND,NSTART
      LOGICAL LCALL
C     ..
C     .. Local Arrays ..
      COMPLEX*16 ALPHHO(0:LMAXD,IEMXD,NTREFD,NSPIND),BESSJW(0:LMX),
     +     BESSYW(0:LMX),DLOGDP(0:LMAXD),HAMF(IRMD,0:LMAXD),
     +     HANKWS(0:LMX),TMATHO(0:LMAXD,IEMXD,NTREFD,NSPIND)
C     ..
C     .. External Subroutines ..
      EXTERNAL BESSEL,IRWSOL,REGSOL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDSQRT,REAL
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA LCALL/.false./
C     ..

      LMAX = LMAXP1 - 1
      IF (LMAXP1.GT.LMX) THEN
        STOP 'cplxwf'

      ELSE


        IF (NT.EQ.1) THEN
          NSTART = 1
          NEND = NVSP/NSPIN

        ELSE
          NSTART = NATREF + 1
          NEND = NVSP/NSPIN
c
c---> restore host alpha - and t - matrix
c
          DO 20 IH = 1,NATREF
            DO 10 L = 0,LMAX
              TMAT(L,IH) = TMATHO(L,IE,IH,ISPIN)
              ALPHA(L,IH) = ALPHHO(L,IE,IH,ISPIN)
   10       CONTINUE
   20     CONTINUE
        END IF

        DO 110 IH = NSTART,NEND
          I1 = NSPIN* (IH-1) + ISPIN
          IF ((KVREL.EQ.1.AND.IH.GT.NATREF) .OR. KVREL.EQ.2) THEN
            KVSRA = 1
            EK = CDSQRT(E+E*E/ (C*C))

          ELSE
            KVSRA = 0
            EK = CDSQRT(E)
          END IF

          IPE = IPAN(IH)
          IRC1 = IRCUT(IPE,IH)
          IRWS1 = IRWS(IH)
          RIRC = R(IRC1,IH)
          ARG = RIRC*EK
          CALL BESSEL(BESSJW,BESSYW,HANKWS,ARG,LMX,LMAXP1,.true.,.true.,
     +                .true.,LCALL)
c
c---> calculate regular wavefunctions
c
          CALL REGSOL(C,E,KVSRA,LMAX,DLOGDP,FZ(1,0,IH),HAMF,MASS(1,IH),
     +                PZ(1,0,IH),DRDI(1,IH),DROR(1,IH),R(1,IH),S(0,IH),
     +                VM2Z(1,I1),Z(IH),IPAN(IH),IRCUT(0,IH),irmd)
c
          EKLFAC = EK
c
          DO 40 L = 0,LMAX
            S1 = S(L,IH)
            RSIRC = RS(IRC1,L,IH)
            EKLFAC = EKLFAC/EK*REAL(2*L+1)
c
c---> determine t - matrix
c
            PN = PZ(IRC1,L,IH)*RSIRC
            N = L + 1
            QF = REAL(L)/RIRC
            HL = HANKWS(L)
            BL = BESSJW(L)
            X = QF*HL - EK*HANKWS(N)
            Y = QF*BL - EK*BESSJW(N)
            W = DLOGDP(L)
            TLSQEZ = (BL*W-Y)/ (X-HL*W)
            TL = TLSQEZ/EK
            TMAT(L,IH) = TL
c
c---> determine the renormalization
c
            ALPHAL = (BL+HL*TLSQEZ)*RIRC/PN
c
c---> determine the alpha matrix
c
            ALPHA(L,IH) = ALPHAL*EKLFAC
c
c---> store host properties
c
            IF (IH.LE.NATREF) THEN
              TMATHO(L,IE,IH,ISPIN) = TMAT(L,IH)
              ALPHHO(L,IE,IH,ISPIN) = ALPHA(L,IH)
            END IF

            DO 30 I = 2,IRC1
              PZ(I,L,IH) = PZ(I,L,IH)*ALPHAL
              FZ(I,L,IH) = FZ(I,L,IH)*ALPHAL
   30       CONTINUE
c
            VALUE = HL*RIRC*RSIRC
            SLOPE = REAL(L+1)*HL - RIRC*EK*HANKWS(L+1)
            SLOPE = (SLOPE*RSIRC+S1/RIRC*VALUE)
            QZ(IRC1,L,IH) = VALUE
            SZ(IRC1,L,IH) = (SLOPE*RIRC- (S1+1.0D0)*VALUE)/
     +                      MASS(IRC1,IH)*DROR(IRC1,IH)
   40     CONTINUE
c
c---> calculate irregular wavefunctions
c
          CALL IRWSOL(EK,LMAX,FZ(1,0,IH),HAMF,MASS(1,IH),PZ(1,0,IH),
     +                QZ(1,0,IH),SZ(1,0,IH),DROR(1,IH),S(0,IH),IPAN(IH),
     +                IRCUT(0,IH),irmd)

          DO 70 L = 0,LMAX
            IF (KVSRA.EQ.1) THEN
              DO 50 I = 2,IRC1
                PZ(I,L,IH) = PZ(I,L,IH)*RS(I,L,IH)
                QZ(I,L,IH) = QZ(I,L,IH)/RS(I,L,IH)
                FZ(I,L,IH) = FZ(I,L,IH)*RS(I,L,IH)/C
                SZ(I,L,IH) = SZ(I,L,IH)/RS(I,L,IH)/C
   50         CONTINUE

            ELSE
              DO 60 I = 2,IRC1
                PZ(I,L,IH) = PZ(I,L,IH)*RS(I,L,IH)
                QZ(I,L,IH) = QZ(I,L,IH)/RS(I,L,IH)
                FZ(I,L,IH) = CZERO
                SZ(I,L,IH) = CZERO
   60         CONTINUE
            END IF
   70     CONTINUE
c
c---> in case of mt-calculation : calculate the wavefunctions between
c             mt-radius and ws-radius analytically
c
          IF (IRC1.LT.IRWS1) THEN
            DO 100 I = IRC1 + 1,IRWS1
              R1 = R(I,IH)
              ARG = R1*EK
c
              CALL BESSEL(BESSJW,BESSYW,HANKWS,ARG,LMX,LMAXP1,.true.,
     +                    .true.,.true.,LCALL)
c
              DO 80 L = 0,LMAX
                TLSQEZ = TMAT(L,IH)*EK
                PZ(I,L,IH) = (BESSJW(L)+TLSQEZ*HANKWS(L))*R1
                QZ(I,L,IH) = HANKWS(L)*R1
                FZ(I,L,IH) = CZERO
                SZ(I,L,IH) = CZERO
   80         CONTINUE
c
c---> calculate small component in case of sra
c
              IF (KVSRA.EQ.1) THEN
                FAC = C/ (E+C*C)
                DO 90 L = 0,LMAX
                  N = L + 1
                  TLSQEZ = TMAT(L,IH)*EK
                  FZ(I,L,IH) = (REAL(L)* (BESSJW(L)+TLSQEZ*HANKWS(L))-
     +                         ARG* (BESSJW(N)+TLSQEZ*HANKWS(N)))*FAC
                  SZ(I,L,IH) = (REAL(L)*HANKWS(L)-ARG*HANKWS(N))*FAC
   90           CONTINUE
c
              END IF

  100       CONTINUE

          END IF

  110   CONTINUE

      END IF

      END
