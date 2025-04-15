      SUBROUTINE GMATSYM(CNK,LNK,IMIN,IMAX,GSYM,NA,GMNLL,LMMAXD,NSEC,
     +                   NGMN,NBASIS,NLEQ,N0,NTERMS,NATOM)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER LMMAXD,NA,NATOM,NBASIS,NGMN,NSEC
C     ..
C     .. Array Arguments ..
      COMPLEX*16 GMNLL(NGMN,LMMAXD),GSYM(NSEC,NSEC)
      REAL*8 CNK(NBASIS,NSEC)
      INTEGER IMAX(NSEC,*),IMIN(NSEC,*),LNK(NBASIS,NSEC),N0(*),NLEQ(*),
     +        NTERMS(*)
C     ..
C     .. Local Scalars ..
      REAL*8 CI,CJ
      INTEGER I,J,KKQ,LAMA,NB,NLMB,NM,NN,NT0A,NT0B,NTXA,NTXB
C     ..
      IF (NLEQ(NA).EQ.0) GO TO 60
      KKQ = NLEQ(NA)
      NT0A = N0(KKQ)
      NTXA = NT0A + NTERMS(KKQ) - 1
      DO 50 NN = NT0A,NTXA
        DO 40 I = IMIN(NN,NA),IMAX(NN,NA)
          CI = CNK(I,NN)
          LAMA = LNK(I,NN)
          DO 30 NB = 1,NATOM
            IF (NLEQ(NB).EQ.0) GO TO 30
            KKQ = NLEQ(NB)
            NT0B = N0(KKQ)
            NTXB = NT0B + NTERMS(KKQ) - 1
            DO 20 NM = NT0B,NTXB
              IF (NM.GT.NN) GO TO 20
              DO 10 J = IMIN(NM,NB),IMAX(NM,NB)
                CJ = CNK(J,NM)
                NLMB = LNK(J,NM) + (NB-1)*LMMAXD
                GSYM(NN,NM) = GSYM(NN,NM) + GMNLL(NLMB,LAMA)*CI*CJ
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
   60 CONTINUE
      RETURN
      END
