      SUBROUTINE SYMBCK(CN,MI,IMIN,IMAX,A,GMN,LMMAXD,NSEC,NBASIS,NT0A,
     +                  NTXA,ISW,N1,N2)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C SIMILAR TO SUBROUTINE SYMBAC, BUT OPTIMIZED FOR DIAGONAL SITES
C FOR ISW=0 IT IS ASSUMED THAT A(N1,N2)=1, A(N2,N1)=1 AND THAT ALL OTHER
C ELEMENTS OF A( , ) ARE SET TO ZERO. SEE CALLING PROGRAM!
C     .. Scalar Arguments ..
      INTEGER ISW,LMMAXD,N1,N2,NBASIS,NSEC,NT0A,NTXA
C     ..
C     .. Array Arguments ..
      REAL*8 A(NSEC,*),CN(NBASIS,*),GMN(LMMAXD,LMMAXD)
      INTEGER IMAX(*),IMIN(*),MI(NBASIS,*)
C     ..
C     .. Local Scalars ..
      REAL*8 CI,CJ
      INTEGER I,IMAXB,IMINB,J,LAMA,LBMB,NM,NN
C     ..
      IF (ISW.EQ.0) THEN
        NN = N1
        DO 20 I = IMIN(NN),IMAX(NN)
          CI = CN(I,NN)
          IF (CI.EQ.0.D0) GO TO 20
          LAMA = MI(I,NN)
C     LEFT INDEX OF GMN
          NM = N2
          IMINB = IMIN(NM)
          IMAXB = IMAX(NM)
          DO 10 J = IMINB,IMAXB
            CJ = CN(J,NM)
            LBMB = MI(J,NM)
C     RIGHT INDEX OF GMN
            GMN(LAMA,LBMB) = GMN(LAMA,LBMB) + CI*CJ
   10     CONTINUE
   20   CONTINUE
        IF (N1.EQ.N2) GO TO 50
        NN = N2
        DO 40 I = IMIN(NN),IMAX(NN)
          CI = CN(I,NN)
          IF (CI.EQ.0.D0) GO TO 40
          LAMA = MI(I,NN)
C     LEFT INDEX OF GMN
          NM = N1
          IMINB = IMIN(NM)
          IMAXB = IMAX(NM)
          DO 30 J = IMINB,IMAXB
            CJ = CN(J,NM)
            LBMB = MI(J,NM)
C     RIGHT INDEX OF GMN
            GMN(LAMA,LBMB) = GMN(LAMA,LBMB) + CI*CJ
   30     CONTINUE
   40   CONTINUE
   50   CONTINUE
      ELSE
        DO 100 NN = NT0A,NTXA
          DO 90 I = IMIN(NN),IMAX(NN)
            CI = CN(I,NN)
            IF (CI.EQ.0.D0) GO TO 90
            LAMA = MI(I,NN)
C     LEFT INDEX OF GMN
            DO 70 NM = NT0A,NTXA
              IMINB = IMIN(NM)
              IMAXB = IMAX(NM)
              IF (IMINB.EQ.IMAXB) GO TO 70
              IF (IMINB+1.EQ.IMAXB) THEN
                J = IMINB + 1
                CJ = CN(J,NM)
                LBMB = MI(J,NM)
C     RIGHT INDEX OF GMN
                GMN(LAMA,LBMB) = GMN(LAMA,LBMB) + CI*A(NN,NM)*CJ
              ELSE
                DO 60 J = IMINB + 1,IMAXB
                  CJ = CN(J,NM)
                  LBMB = MI(J,NM)
C     RIGHT INDEX OF GMN
                  GMN(LAMA,LBMB) = GMN(LAMA,LBMB) + CI*A(NN,NM)*CJ
   60           CONTINUE
              END IF
   70       CONTINUE
            DO 80 NM = NT0A,NTXA
              J = IMIN(NM)
              CJ = CN(J,NM)
              LBMB = MI(J,NM)
C     RIGHT INDEX OF GMN
              GMN(LAMA,LBMB) = GMN(LAMA,LBMB) + CI*A(NN,NM)*CJ
   80       CONTINUE
   90     CONTINUE
  100   CONTINUE
      END IF
      RETURN
      END
