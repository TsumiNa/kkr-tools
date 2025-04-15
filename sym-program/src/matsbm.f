      SUBROUTINE MATSBM(CNK,LNK,JMINA,JMAXA,JMINB,JMAXB,NTXA,GSYM,GMN,
     +                  LMMAXD,NSEC,NBASIS,NUMBER,NCOL,NT0B,NTXB,NT1A)
      IMPLICIT NONE
C     .. Scalar Arguments ..
      INTEGER LMMAXD,NBASIS,NCOL,NSEC,NT0B,NT1A,NTXA,NTXB,NUMBER
C     ..
C     .. Array Arguments ..
      REAL*8 CNK(NBASIS,*),GMN(LMMAXD,LMMAXD),GSYM(NSEC,NSEC)
      INTEGER JMAXA(*),JMAXB(*),JMINA(*),JMINB(*),LNK(NBASIS,*)
C     ..
C     .. Local Scalars ..
      REAL*8 CI,CI1,CI2,CJ,CJ1,CJ2,CTOL
      INTEGER I,IMAXA,IMAXB,IMINA,IMINB,J,LAMA,LAMA1,LAMA2,LBMB,LBMB1,
     +        LBMB2,NM,NN,NTXBA
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,MIN0
C     ..
      CTOL = 1.D-8
      DO 180 NN = NT1A,NTXA
        IMINA = JMINA(NN)
        IMAXA = JMAXA(NN)
        NTXBA = MIN0(NN,NTXB)
        IF (IMINA.EQ.IMAXA) GO TO 80
        IF (IMINA.EQ.IMAXA-1) GO TO 130
        DO 70 NM = NT0B,NTXBA
          IMINB = JMINB(NM)
          IMAXB = JMAXB(NM)
          IF (IMINB.EQ.IMAXB) GO TO 30
          IF (IMINB.EQ.IMAXB-1) GO TO 50
          DO 20 J = IMINB,IMAXB
            CJ = CNK(J,NM)*NUMBER/NCOL
            IF (DABS(CJ).LT.CTOL) GO TO 20
            LBMB = LNK(J,NM)
            DO 10 I = IMINA,IMAXA
              CI = CNK(I,NN)
              LAMA = LNK(I,NN)
              GSYM(NN,NM) = GSYM(NN,NM) + GMN(LBMB,LAMA)*CI*CJ
   10       CONTINUE
   20     CONTINUE
          GO TO 70
   30     CONTINUE
          J = IMINB
          CJ = CNK(J,NM)*NUMBER/NCOL
          IF (DABS(CJ).LT.CTOL) GO TO 70
          LBMB = LNK(J,NM)
          DO 40 I = IMINA,IMAXA
            CI = CNK(I,NN)
            LAMA = LNK(I,NN)
            GSYM(NN,NM) = GSYM(NN,NM) + GMN(LBMB,LAMA)*CI*CJ
   40     CONTINUE
          GO TO 70
   50     CONTINUE
          J = IMINB
          CJ1 = CNK(J,NM)*NUMBER/NCOL
          LBMB1 = LNK(J,NM)
          J = IMINB + 1
          CJ2 = CNK(J,NM)*NUMBER/NCOL
          IF (DABS(CJ1)+DABS(CJ2).LT.CTOL) GO TO 70
          LBMB2 = LNK(J,NM)
          DO 60 I = IMINA,IMAXA
            CI = CNK(I,NN)
            LAMA = LNK(I,NN)
            GSYM(NN,NM) = GSYM(NN,NM) + (GMN(LBMB1,LAMA)*CJ1+
     +                    GMN(LBMB2,LAMA)*CJ2)*CI
   60     CONTINUE
   70   CONTINUE
        GO TO 180
   80   CONTINUE
        I = IMINA
        CI = CNK(I,NN)
        IF (DABS(CI).LT.CTOL) GO TO 180
        LAMA = LNK(I,NN)
        DO 120 NM = NT0B,NTXBA
          IMINB = JMINB(NM)
          IMAXB = JMAXB(NM)
          IF (IMINB.EQ.IMAXB) GO TO 100
          IF (IMINB.EQ.IMAXB-1) GO TO 110
          DO 90 J = IMINB,IMAXB
            CJ = CNK(J,NM)*NUMBER/NCOL
            LBMB = LNK(J,NM)
            GSYM(NN,NM) = GSYM(NN,NM) + GMN(LBMB,LAMA)*CI*CJ
   90     CONTINUE
          GO TO 120
  100     CONTINUE
          J = IMINB
          CJ = CNK(J,NM)*NUMBER/NCOL
          LBMB = LNK(J,NM)
          GSYM(NN,NM) = GSYM(NN,NM) + GMN(LBMB,LAMA)*CI*CJ
          GO TO 120
  110     CONTINUE
          J = IMINB
          CJ1 = CNK(J,NM)*NUMBER/NCOL
          LBMB1 = LNK(J,NM)
          J = IMINB + 1
          CJ2 = CNK(J,NM)*NUMBER/NCOL
          LBMB2 = LNK(J,NM)
          GSYM(NN,NM) = GSYM(NN,NM) + CI*
     +                  (GMN(LBMB1,LAMA)*CJ1+GMN(LBMB2,LAMA)*CJ2)
  120   CONTINUE
        GO TO 180
  130   CONTINUE
        I = IMINA
        CI1 = CNK(I,NN)
        LAMA1 = LNK(I,NN)
        I = IMINA + 1
        CI2 = CNK(I,NN)
        IF (DABS(CI1)+DABS(CI2).LT.CTOL) GO TO 180
        LAMA2 = LNK(I,NN)
        DO 170 NM = NT0B,NTXBA
          IMINB = JMINB(NM)
          IMAXB = JMAXB(NM)
          IF (IMINB.EQ.IMAXB) GO TO 150
          IF (IMINB.EQ.IMAXB-1) GO TO 160
          DO 140 J = IMINB,IMAXB
            CJ = CNK(J,NM)*NUMBER/NCOL
            LBMB = LNK(J,NM)
            GSYM(NN,NM) = GSYM(NN,NM) + (GMN(LBMB,LAMA1)*CI1+
     +                    GMN(LBMB,LAMA2)*CI2)*CJ
  140     CONTINUE
          GO TO 170
  150     CONTINUE
          J = IMINB
          CJ = CNK(J,NM)*NUMBER/NCOL
          LBMB = LNK(J,NM)
          GSYM(NN,NM) = GSYM(NN,NM) + (GMN(LBMB,LAMA1)*CI1+
     +                  GMN(LBMB,LAMA2)*CI2)*CJ
          GO TO 170
  160     CONTINUE
          J = IMINB
          CJ1 = CNK(J,NM)*NUMBER/NCOL
          LBMB1 = LNK(J,NM)
          J = IMINB + 1
          CJ2 = CNK(J,NM)*NUMBER/NCOL
          LBMB2 = LNK(J,NM)
          GSYM(NN,NM) = GSYM(NN,NM) + (GMN(LBMB1,LAMA1)*CI1+
     +                  GMN(LBMB1,LAMA2)*CI2)*CJ1 +
     +                  (GMN(LBMB2,LAMA1)*CI1+GMN(LBMB2,LAMA2)*CI2)*CJ2
  170   CONTINUE
  180 CONTINUE
      RETURN
      END
