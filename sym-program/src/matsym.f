      SUBROUTINE MATSYM(CN,MI,IMIN,IMAX,NUMCOL,GSYM,GMN,LMMAXD,NSEC,
     +                  NDIMNP,LN,SIGNLN,NEQ,NBASIS,NCOL,NLEQ,NUMBER,N0,
     +                  NTERMS,NAF,NBF,KGMN,NATMX)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER KGMN,LMMAXD,NATMX,NBASIS,NCOL,NDIMNP,NSEC,NUMCOL
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,NSEC,*),GMN(LMMAXD,LMMAXD,*),
     +                 GSYM(NSEC,NSEC),SIGNLN(*)
      INTEGER IMAX(NSEC,NATMX,NUMCOL),IMIN(NSEC,NATMX,NUMCOL),LN(*),
     +        MI(NBASIS,NSEC,*),N0(*),NAF(*),NBF(*),NEQ(*),NLEQ(*),
     +        NTERMS(*),NUMBER(*)
C     ..
C     .. Local Scalars ..
      REAL*8 SA,SB
      INTEGER IGMN3,KKQ,NA,NB,NK,NM,NN,NT0A,NT0B,NT1A,NTXA,NTXB
C     ..
C     .. External Subroutines ..
      EXTERNAL MATSBM,MATSBP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
      DO 40 IGMN3 = 1,KGMN
        NA = NAF(IGMN3)
        NB = NBF(IGMN3)
        IF (NEQ(NB).NE.0 .AND. NEQ(NA).NE.0 .OR.
     +      NLEQ(NA)*NLEQ(NB).EQ.0) GO TO 40
        IF (NEQ(NB).NE.0) GO TO 20
        KKQ = NLEQ(NA)
        NT0A = N0(KKQ)
        NTXA = NT0A + NTERMS(KKQ) - 1
        KKQ = NLEQ(NB)
        NT0B = N0(KKQ)
        NTXB = NT0B + NTERMS(KKQ) - 1
        NT1A = MAX0(NT0A,NT0B)
        DO 10 NK = 1,NCOL
          CALL MATSBP(CN(1,1,NK),MI(1,1,NK),IMIN(1,NA,NK),IMAX(1,NA,NK),
     +                IMIN(1,NB,NK),IMAX(1,NB,NK),NTXA,GSYM,
     +                GMN(1,1,IGMN3),LMMAXD,NSEC,NBASIS,NUMBER(NB),NCOL,
     +                NT0B,NTXB,NT1A)
   10   CONTINUE
   20   CONTINUE
        IF (NEQ(NA).NE.0) GO TO 40
        IF (NA.EQ.NB) GO TO 40
        NB = NAF(IGMN3)
        NA = NBF(IGMN3)
        KKQ = NLEQ(NA)
        NT0A = N0(KKQ)
        NTXA = NT0A + NTERMS(KKQ) - 1
        KKQ = NLEQ(NB)
        NT0B = N0(KKQ)
        NTXB = NT0B + NTERMS(KKQ) - 1
        NT1A = MAX0(NT0A,NT0B)
        DO 30 NK = 1,NCOL
          CALL MATSBM(CN(1,1,NK),MI(1,1,NK),IMIN(1,NA,NK),IMAX(1,NA,NK),
     +                IMIN(1,NB,NK),IMAX(1,NB,NK),NTXA,GSYM,
     +                GMN(1,1,IGMN3),LMMAXD,NSEC,NBASIS,NUMBER(NB),NCOL,
     +                NT0B,NTXB,NT1A)
   30   CONTINUE
   40 CONTINUE
      IF (NDIMNP.EQ.1) GO TO 70
      DO 60 NN = 2,NDIMNP
        SA = SIGNLN(NN)
        DO 50 NM = 1,NN - 1
          SB = SIGNLN(NM)*SA
          GSYM(NM,NN) = GSYM(NN,NM)*SB
   50   CONTINUE
   60 CONTINUE
   70 RETURN
      END
