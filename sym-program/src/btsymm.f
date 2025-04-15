      SUBROUTINE BTSYMM(LMAXSQ,ND,RIMP,NTYPIM,NTYPHO,NIMP,JSHELL,NSHELL,
     +                  NREP,IGROUP,LBTSYM,CN,STR,IMAX,IMIN,ISTR,MAXSTR,
     +                  MI,MINSTR,MSTR,N0ATOM)
      IMPLICIT NONE
      
C     .. Parameters ..
      INTEGER LMAXD,NIMPD,NSHELD
      PARAMETER (LMAXD=4,NIMPD=4290,NSHELD=4020)
      INTEGER NATMX,LMMAXD
      PARAMETER (NATMX=NIMPD,LMMAXD= (LMAXD+1)**2)
      INTEGER NUMCOL,JGMN,NGR,JSHLL
      PARAMETER (NUMCOL=3,JGMN=556,NGR=51577,JSHLL=NSHELD)
      INTEGER NCOEFF
      PARAMETER (NCOEFF=130000)
      INTEGER NSEC,NREPMX,NBASIS,NISTR,NSTR,NMSTR
      PARAMETER (NSEC=700,NREPMX=10,NBASIS=3000,
     +          NISTR=2600000,NSTR=2600000,
     +          NMSTR=1800000)
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,NSEC,NUMCOL),RIMP(3,NIMPD),STR(NSTR)
      INTEGER IMAX(NSEC,NATMX,NUMCOL),IMIN(NSEC,NATMX,NUMCOL),
     +        ISTR(NISTR),JSHELL(5,NIMPD),MAXSTR(NMSTR),
     +        MI(NBASIS,NSEC,NUMCOL),MINSTR(NMSTR),MSTR(NMSTR),
     +        N0ATOM(NBASIS,NSEC,NUMCOL),ND(48,3,3),NTYPHO(NIMPD),
     +        NTYPIM(NIMPD)
C     ..
C     .. Local Scalars ..
      REAL*8 GBCF,GFF
      INTEGER I,I1,I2,IC,ICG,IE,IG,IIFL,IIFM,IIFX,IS,ITFEND,J,KGMN,KKQ,
     +        LMAX,M,MILEN,N,N1,N2,NA,NATOM,NCOL,NDM,NK,NLAST,NN,NP,
     +        NREPQ,NSH,NSV,NSW,NT0A,NTXA,NV
C     ..
C     .. Local Arrays ..
      REAL*8 GCOEFF(NCOEFF),GHOST(NSEC,NSEC),
     +                 GM(LMMAXD,LMMAXD,JSHLL),GMN(LMMAXD,LMMAXD,JGMN),
     +                 GR(NGR),SIGNLN(NSEC,NREPMX)
      INTEGER GSTRC(LMMAXD,LMMAXD,JSHLL),I1F(NCOEFF),I2F(NCOEFF),
     +        I3F(NCOEFF),I4F(NCOEFF),I5F(NCOEFF),LN(NSEC,NREPMX),
     +        LOFLM(25),N0(NATMX),NAF(JGMN),NBF(JGMN),NDG(NREPMX),
     +        NDIM(NREPMX),NEQ(NATMX),NFM(10),NLEQ(NATMX),
     +        NMS(NSEC,NUMCOL),NSHF(100),NTERMS(NATMX),NUMBER(NATMX)
      LOGICAL ISZ(JSHLL)
C     ..
C     .. External Subroutines ..
      EXTERNAL G1MN,I4INIT,MATSYM,MICHK,R8INIT,SYMBCK,SYMRCL,SYMSTO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,NIMP,NREP,NSHELL
      LOGICAL LBTSYM
      CHARACTER*5 IGROUP
C     ..
C     .. Data statements ..
      DATA LOFLM/1,3*2,5*3,7*4,9*5/
C     ..
      MILEN = 1
      CALL R8INIT(CN,NBASIS*NSEC*NUMCOL,0.0D0)
      CALL R8INIT(GR,NGR,0.0D0)
      CALL I4INIT(MI,NBASIS*NSEC*NUMCOL,0)
      CALL I4INIT(LN,NSEC*NREPMX,0)
      CALL I4INIT(NMS,NSEC*NUMCOL,0)
      REWIND 19
      READ (19,FMT=9010) NREPQ,NATOM,LMAX,NSH, (NSHF(I),I=1,NSH)
      WRITE (35,FMT=9030) NREP,NATOM,LMAX,NSH, (NSHF(I),I=1,NSH)
      WRITE (6,FMT=9020) NREP,NATOM,LMAX
      WRITE (35,FMT=9000) ((RIMP(I,M),I=1,3),NTYPHO(M),M,NTYPIM(M),M=1,
     +  NIMP)
      CALL R8INIT(GHOST,NSEC*NSEC,0.0D0)
      NEQ(1) = 0
      NSW = NTYPIM(1)
      NSV = 0
      DO 20 N = 2,NATOM
        IF (NSW.EQ.NTYPIM(N)) GO TO 10
        NEQ(N) = 0
        NSV = N
        NSW = NTYPIM(N)
        GO TO 20
   10   NEQ(N) = NSV
   20 CONTINUE
      WRITE (6,FMT=9040) NSTR,NISTR,NMSTR,NBASIS,NSEC,NUMCOL,NBASIS,
     +  NSEC,NUMCOL,NSEC,NATMX,NUMCOL,NSEC,NATMX,NUMCOL,NSEC,NSEC
C
      CALL SYMSTO(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,N0ATOM,IMIN,
     +            IMAX,NATOM,NSTR,NISTR,NMSTR,NSEC,NUMCOL,NBASIS,NREP,
     +            NMS,LN,LMAX,NLEQ,N0,NTERMS,NATMX,NDIM,NDG)
      WRITE (35,FMT=9050) (NDIM(NP),NP=1,NREP)
      WRITE (35,FMT=9050) (NDG(NP),NP=1,NREP)
      WRITE (6,FMT=9060) NATOM
      DO 30 N = 1,NATOM
        NUMBER(N) = 0
   30 CONTINUE
      DO 60 N = 1,NATOM
        IF (NEQ(N).NE.0) GO TO 40
        NUMBER(N) = 1
        GO TO 50
   40   NUMBER(NEQ(N)) = NUMBER(NEQ(N)) + 1
   50   CONTINUE
   60 CONTINUE
      MILEN = 1
C----------------------------------------------------------------------
C     BACK-SYMMETRIZATION
C----------------------------------------------------------------------
      IIFL = 1
      IIFX = 1
      IIFM = 1
      DO 90 NP = 1,NREP
        CALL SYMRCL(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,IMIN,IMAX,
     +              IIFL,IIFX,IIFM,NATOM,LN(1,NP),NBASIS,NSEC,NUMCOL,NP,
     +              NLEQ,N0,NTERMS,NATMX,NDIM,NDG)
        NCOL = NDG(NP)/2
        DO 80 NA = 1,NATOM
          IF (NEQ(NA).NE.0) GO TO 80
          IF (NLEQ(NA).EQ.0) GO TO 80
          KKQ = NLEQ(NA)
          NT0A = N0(KKQ)
          NTXA = NT0A + NTERMS(KKQ) - 1
          DO 70 NK = 1,NCOL
            CALL MICHK(CN(1,1,NK),MI(1,1,NK),IMIN(1,NA,NK),
     +                 IMAX(1,NA,NK),NBASIS,NT0A,NTXA,MILEN)
   70     CONTINUE
   80   CONTINUE
   90 CONTINUE
      IIFL = 1
      IIFX = 1
      IIFM = 1
      CALL R8INIT(GMN,NATMX*LMMAXD*LMMAXD,0.0D0)
      DO 180 NP = 1,NREP
        IC = 0
        CALL R8INIT(GHOST,NSEC*NSEC,0.0D0)
        CALL SYMRCL(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,IMIN,IMAX,
     +              IIFL,IIFX,IIFM,NATOM,LN(1,NP),NBASIS,NSEC,NUMCOL,NP,
     +              NLEQ,N0,NTERMS,NATMX,NDIM,NDG)
        NDM = NDIM(NP)
        NCOL = NDG(NP)/2
        DO 170 NN = 1,NTYPIM(NATOM)
          NSW = 1
          DO 160 NA = 1,NATOM
            IF (NSW.EQ.0 .OR. NN.NE.NTYPIM(NA)) GO TO 160
            NSW = 0
            IF (NEQ(NA).NE.0) GO TO 160
            IF (NLEQ(NA).EQ.0) GO TO 160
            KKQ = NLEQ(NA)
            NT0A = N0(KKQ)
            NTXA = NT0A + NTERMS(KKQ) - 1
            DO 150 N1 = NT0A,NTXA
              DO 140 N2 = N1,NTXA
                DO 100 NK = 1,NCOL
                  CALL SYMBCK(CN(1,1,NK),MI(1,1,NK),IMIN(1,NA,NK),
     +                        IMAX(1,NA,NK),GHOST,GMN(1,1,NA),LMMAXD,
     +                        NSEC,NBASIS,NT0A,NTXA,0,N1,N2)
  100           CONTINUE
                DO 130 I1 = 1,LMAXSQ
                  DO 110 I2 = I1,LMAXSQ
C**** THE PURPOSE OF THE NEXT STATEMENTS (DIFFERENT FROM THOSE
C**** OF THE ORIGINAL VERSION OF THIS PROGRAM),IS TO REDUCE THE
C**** DIMENSION OF THE ARRAY CONTAINING THE BACKSYMMETRIZED
C**** GREEN'S FUNCTIONS (SEE IMPURITY PROGRAM  ARRAY DTB)
                    GBCF = (GMN(I1,I2,NA)+GMN(I2,I1,NA))*0.5D0
                    IF (DABS(GBCF).LT.1.D-6) GO TO 110
                    IC = IC + 1
                    I1F(IC) = I1
                    I2F(IC) = I2
                    I3F(IC) = N1
                    I4F(IC) = N2
                    I5F(IC) = NN
                    GCOEFF(IC) = GBCF
  110             CONTINUE
                  GMN(I1,I1,NA) = 0.0D0
                  DO 120 I2 = I1 + 1,LMAXSQ
                    GMN(I1,I2,NA) = 0.0D0
                    GMN(I2,I1,NA) = 0.0D0
  120             CONTINUE
  130           CONTINUE
  140         CONTINUE
  150       CONTINUE
  160     CONTINUE
  170   CONTINUE
        ICG = IC
        WRITE (33) NP,ICG, (GCOEFF(J),J=1,ICG), (I1F(J),J=1,ICG),
     +    (I2F(J),J=1,ICG), (I5F(J),J=1,ICG), (I3F(J),J=1,ICG),
     +    (I4F(J),J=1,ICG)
        WRITE (6,FMT=9080) ICG
  180 CONTINUE
      CLOSE (33)
      IF (.NOT.LBTSYM) THEN
        WRITE (6,FMT=*) 'tu coefficients are not calculated'
        RETURN
      END IF
C----------------------------------------------------------------------
C     SYMMETRIZATION FOR THE SITE-DIAGONAL TLL' MATRIX
C----------------------------------------------------------------------
      NV = NATOM
      NTYPIM(NATOM) = NTYPIM(NATOM)
      DO 190 IS = 1,NSHELL
        ITFEND = 0
        WRITE (35,FMT='(A5,43X,I3)') IGROUP,ITFEND
  190 CONTINUE
      IG = 0
      DO 220 IS = 1,NSHELL
        DO 210 J = 1,LMAXSQ
          DO 200 I = 1,LMAXSQ
            IG = IG + 1
            GSTRC(I,J,IS) = IG
  200     CONTINUE
  210   CONTINUE
  220 CONTINUE
      NLAST = IG
      IC = 0
      CALL R8INIT(GMN,LMMAXD*LMMAXD*JGMN,0.0D0)
      DO 270 IE = 1,NLAST
        CALL G1MN(GMN,GM,GSTRC,LMMAXD,LMMAXD,NATOM,IE,NAF,NBF,KGMN,ISZ,
     +            GR,NSHELL,NTYPIM)
        IF (KGMN.EQ.0) GO TO 270
        IIFL = 1
        IIFX = 1
        IIFM = 1
        DO 260 NP = 1,NREP
          CALL SYMRCL(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,IMIN,IMAX,
     +                IIFL,IIFX,IIFM,NATOM,LN(1,NP),NBASIS,NSEC,NUMCOL,
     +                NP,NLEQ,N0,NTERMS,NATMX,NDIM,NDG)
          CALL R8INIT(GHOST,NSEC*NSEC,0.0D0)
          NCOL = NDG(NP)/2
          DO 230 NN = 1,NDIM(NP)
            SIGNLN(NN,NP) = (-1)**LN(NN,NP)
  230     CONTINUE
          CALL MATSYM(CN,MI,IMIN,IMAX,NUMCOL,GHOST,GMN,LMMAXD,NSEC,
     +                NDIM(NP),LN(1,NP),SIGNLN(1,NP),NEQ,NBASIS,NCOL,
     +                NLEQ,NUMBER,N0,NTERMS,NAF,NBF,KGMN,NATMX)
          NDM = NDIM(NP)
          DO 250 I2 = 1,NDM
            DO 240 I1 = I2,NDM
              GFF = GHOST(I1,I2)
              IF (DABS(GFF).LT.1.D-6) GO TO 240
              IC = IC + 1
              GCOEFF(IC) = GFF
              I1F(IC) = I1
              I2F(IC) = I2
              I3F(IC) = NP
              I4F(IC) = IE
  240       CONTINUE
  250     CONTINUE
  260   CONTINUE
        IF (IC.GT.NCOEFF) GO TO 300
  270 CONTINUE
      WRITE (35,FMT=9095) IC
      WRITE (35,FMT=9100) (GCOEFF(J),J=1,IC)
      WRITE (35,FMT=9090) (I1F(J),J=1,IC)
      WRITE (35,FMT=9090) (I2F(J),J=1,IC)
      WRITE (35,FMT=9090) (I3F(J),J=1,IC)
      WRITE (35,FMT=9090) (I4F(J),J=1,IC)
      DO 280 NP = 1,NREP
        WRITE (35,FMT=9090) (LN(J,NP),J=1,NDIM(NP))
  280 CONTINUE
      WRITE (6,FMT=9070) IC
  290 CONTINUE
      RETURN
  300 WRITE (6,FMT=9110) IC,NCOEFF
      STOP 55
C9010 FORMAT (A5,4X,L1,4X,L1,4X,L1,2I5)
C9020 FORMAT (A50)
C9040 FORMAT (' LMAX=',I1,' ON INPUT AND LMX=',I1,' IN DIMENSION',
C    +       ' STATEMENTS ARE INCONSISTENT ')
 9000 FORMAT (3F10.6,I2,I4,I5)
 9010 FORMAT (8I5)
 9020 FORMAT (I5,'  REPRESENTATIONS CONTRIBUTE ',/,I5,'  CLUSTER ATOMS',
     +       /,' LMAX=',I2,'  IS USED')
 9030 FORMAT (11I5)
 9040 FORMAT (1H ,33H******************************** ,/,
     +       ' THE FOLLOWING ARRAYS HAVE BEEN ALLOCATED',/,'  STR(',I7,
     +       ')  ISTR(',I7,')  MSTR(',I7,')',/,'  CN(',I3,',',I4,',',I1,
     +       ') MI(',I3,',',I4,',',I1,')',/,' IMIN(',I4,',',I4,',',I1,
     +       ')  IMAX(',I4,',',I4,',',I1,')',/,' GHOST(',I4,',',I4,')')
 9050 FORMAT (10I5)
 9060 FORMAT (1H0,30X,18HNUMBER OF CENTERS=,I4)
 9070 FORMAT (I9,' tu COEFFICIENTS WRITTEN')
 9080 FORMAT (I9,' back COEFFICIENTS WRITTEN')
 9095 FORMAT (I6)
 9090 FORMAT (16I5)
 9100 FORMAT (1P,5D16.9)
 9110 FORMAT ('  ERROR STOP IN MAIN: DIMENSIONS TOO SMALL',/,
     &'  IC is:',I8,10x,'NCOEFF is:',I8 ,/,
     & 'NCOEFF must be larger than IC')
      END
