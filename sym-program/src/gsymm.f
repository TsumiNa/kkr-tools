      SUBROUTINE GSYMM(LMAXSQ,NEZ,DROT,IROTMN,ISHELL,NIMP,NSHELL,MSHELL,
     +                 IJ,NHSPIN,CN,STR,IMAX,IMIN,ISTR,MAXSTR,MI,MINSTR,
     +                 MSTR,N0ATOM,ihandle)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NIMPD,NSHELD
      PARAMETER (NIMPD=4290,NSHELD=4020)
      INTEGER NSEC,NREPMX,NBASIS,NISTR,NSTR,NMSTR
      PARAMETER (NSEC=700,NREPMX=10,NBASIS=3000
     +          ,NISTR=2600000,NSTR=2600000,
     +          NMSTR=1800000)
      INTEGER LMAXD,NRATOM
      PARAMETER (LMAXD=4,NRATOM=4290)
      INTEGER LMMAXD,NGD
      PARAMETER (LMMAXD= (LMAXD+1)**2,NGD=LMMAXD*NIMPD)
      INTEGER NUMCOL
      PARAMETER (NUMCOL=3)
      COMPLEX*16 CZERO
      PARAMETER (CZERO= (0.D0,0.D0))
C     ..
C     .. Local Scalars ..
      COMPLEX*16 DF,EZ,ezk
      COMPLEX DF80,EZ80
      INTEGER I,I1,I2,IE,IIFL,IIFM,IIFX,LM1,LM2,LMAX,M,MS,NATM,NATOM,NM,
     +        NN,NP,NREP,NS,NSH,NSHDUM,ilin,ndimnp,nlin
      LOGICAL LSTART
C     ..
C     .. Local Arrays ..
      COMPLEX*16 GLLNEW(LMMAXD,LMMAXD,NSHELD),
     +               GM(LMMAXD,LMMAXD,NSHELD),GMN(NGD,LMMAXD),
     +               GSYM(NSEC,NSEC)
      COMPLEX*16 GTEMP(NSEC,NSEC)
      Complex*8 glin(nsec*nsec)
      INTEGER LN(NSEC,NREPMX),N0(NRATOM),NDG(NREPMX),NDIM(NREPMX),
     +        NFM(10),NLEQ(NRATOM),NMS(NSEC,NUMCOL),NTERMS(NRATOM)
C     ..
C     .. External Subroutines ..
      EXTERNAL C8INIT,G0MN,GMATSYM,I4INIT,R8INIT,SYMRCL,SYMSTO
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,NSEC,NUMCOL),DROT(48,LMMAXD,LMMAXD),STR(NSTR)
      INTEGER IJ(NSHELD),IMAX(NSEC,NRATOM,NUMCOL),
     +        IMIN(NSEC,NRATOM,NUMCOL),IROTMN(NIMPD,0:NIMPD),
     +        ISHELL(NIMPD,0:NIMPD),ISTR(NISTR),MAXSTR(NMSTR),
     +        MI(NBASIS,NSEC,NUMCOL),MINSTR(NMSTR),MSTR(NMSTR),
     +        N0ATOM(NBASIS,NSEC,NUMCOL)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,MSHELL,NEZ,NHSPIN,NIMP,NSHELL,ihandle
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Data statements ..
      DATA LSTART/.TRUE./
C     ..
      IF (LSTART) THEN
        CALL R8INIT(CN,NBASIS*NSEC*NUMCOL,0.0D0)
        CALL I4INIT(MI,NBASIS*NSEC*NUMCOL,0)
        CALL I4INIT(LN,NSEC*NREPMX,0)
        CALL I4INIT(NMS,NSEC*NUMCOL,0)
        READ (19,FMT=9050) NREP,NATM,LMAX,NSH, (NSHDUM,I=1,NSH)
        NATOM = NIMP
        WRITE (6,FMT=9070) NREP,NATM,LMAX
c     IF (NREP.GT.NREPMX .OR. NATM.NE.NATOM .OR. LMAXD.NE.LMAX) STOP 22
        WRITE (6,FMT=9080) NSTR,NISTR,NMSTR,NBASIS,NSEC,NUMCOL,NBASIS,
     +    NSEC,NUMCOL,NSEC,NRATOM,NUMCOL,NSEC,NRATOM,NUMCOL,NSEC,NSEC
C
        CALL SYMSTO(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,N0ATOM,IMIN,
     +              IMAX,NATOM,NSTR,NISTR,NMSTR,NSEC,NUMCOL,NBASIS,NREP,
     +              NMS,LN,LMAXD,NLEQ,N0,NTERMS,NRATOM,NDIM,NDG)
        WRITE (6,FMT=9090) NATOM
        LSTART = .FALSE.
      END IF
      CALL C8INIT(GM,NGD*LMMAXD,CZERO)
      DO 130 IE = 1,NEZ
        READ (15) EZ,DF
        WRITE (6,FMT=*) IE,EZ,DF,MSHELL
        DO 10 MS = 1,MSHELL

      READ (15) ((GLLNEW(LM2,LM1,MS),LM2=1,LMAXSQ),LM1=1,LMAXSQ)
c	write(6,*) ' ms gllnew',ms,GLLNEW(1,1,1)
c	stop
   10   CONTINUE
c       write(6,*) ' after 10'
        DO 60 NS = 1,NSHELL
          IF (IJ(NS).NE.0) THEN
            DO 30 LM1 = 1,LMAXSQ
              DO 20 LM2 = 1,LMAXSQ
                GM(LM2,LM1,NS) = GLLNEW(LM2,LM1,IJ(NS))
   20         CONTINUE
   30       CONTINUE
          ELSE
            DO 50 LM1 = 1,LMAXSQ
              DO 40 LM2 = 1,LMAXSQ
                GM(LM2,LM1,NS) = (0.0D0,0.0D0)
   40         CONTINUE
   50       CONTINUE
          END IF

   60   CONTINUE
c        write(6,*) ' after 60'
C----------------------------------------------------------------------
C     SYMMETRIZATION FOR THE GREEN FUNCTION MATRIX
C----------------------------------------------------------------------
        IIFL = 1
        IIFX = 1
        IIFM = 1
        DO 120 NP = 1,NREP
c	write(6,*) 'before symrcl'
          CALL SYMRCL(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,IMIN,IMAX,
     +                IIFL,IIFX,IIFM,NATOM,LN,NBASIS,NSEC,NUMCOL,NP,
     +                NLEQ,N0,NTERMS,NRATOM,NDIM,NDG)
          CALL C8INIT(GSYM,NSEC*NSEC,CZERO)
          DO 70 M = 1,NATOM
            IF (NLEQ(M).EQ.0) GO TO 70
            CALL G0MN(GMN,GM,NIMP,LMAXSQ,M,DROT,IROTMN,ISHELL)
            CALL GMATSYM(CN,MI,IMIN,IMAX,GSYM,M,GMN,LMMAXD,NSEC,NGD,
     +                   NBASIS,NLEQ,N0,NTERMS,NATOM)
   70     CONTINUE
          DO 90 NN = 2,NDIM(NP)
            DO 80 NM = 1,NN - 1
              GSYM(NM,NN) = GSYM(NN,NM)
   80       CONTINUE
   90     CONTINUE
c
c---> store struc. host green's functions on disc
c
          DO 110 I1 = 1,NDIM(NP)
            DO 100 I2 = 1,I1
              GTEMP(I2,I1) = GSYM(I2,I1)
  100       CONTINUE
  110     CONTINUE
          EZ80 = EZ
          DF80 = DF
       WRITE (80) EZ,SQRT(EZ),DF,((GTEMP(I2,I1),I2=1,I1),I1=1,
     +      NDIM(NP))
c$$$          WRITE (80) EZ80,SQRT(EZ80),DF80,
c$$$     +      ((GTEMP(I2,I1),I2=1,I1),I1=1,NDIM(NP))
            ilin = 1
            ndimnp = ndim(np)
            Do i1 = 1,ndimnp
              Do i2 = 1,i1
                glin(ilin) = gtemp(i2,i1)
                ilin = ilin + 1
              End Do
            End Do
            nlin = ndimnp* (ndimnp+1)/2
            If ((ilin-1).ne.nlin) Stop 'READ FXDR'
            ezk = sqrt(ez)
c           Call fxdrdbl(ihandle,ez,2)
c      write(ihandle) ez
c            Call fxdrdbl(ihandle,ezk,2)
c      write(ihandle) ezk
c            Call fxdrdbl(ihandle,df,2)
c      write(ihandle) df
c            Call fxdrrl(ihandle,glin,2*nlin)
c       write(ihandle) (glin(i1),i1=1,nlin)
  120   CONTINUE
c      write(6,*) ' after 120'
  130 CONTINUE

      RETURN
 9000 FORMAT (I5)
 9010 FORMAT (1P,6D12.5)
 9020 FORMAT (I5)
 9030 FORMAT (1P,6D12.5)
 9040 FORMAT (3F10.6,I1,I4,I5)
 9050 FORMAT (8I5)
 9060 FORMAT (' LMAX=',I1,' ON INPUT AND LMX=',I1,' IN DIMENSION',
     +       ' STATEMENTS ARE INCONSISTENT ')
 9070 FORMAT (I5,'  REPRESENTATIONS CONTRIBUTE ',/,I5,'  CLUSTER ATOMS',
     +       /,' LMAX=',I2,'  IS USED')
 9080 FORMAT (1H ,33H******************************** ,/,
     +       ' THE FOLLOWING ARRAYS HAVE BEEN ALLOCATED',/,'  STR(',I7,
     +       ')  ISTR(',I7,')  MSTR(',I7,')',/,'  CN(',I5,',',I4,',',I1,
     +       ') MI(',I5,',',I4,',',I1,')',/,' IMIN(',I4,',',I4,',',I1,
     +       ')  IMAX(',I4,',',I4,',',I1,')',/,' GSYM(',I4,',',I4,')')
 9090 FORMAT (1H0,30X,18HNUMBER OF CENTERS=,I4)
 9100 FORMAT (1P,5D16.9)
      END
