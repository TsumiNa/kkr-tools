      SUBROUTINE MICHK(CN,MI,IMIN,IMAX,NBASIS,NT0A,NTXA,MILEN)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C SIMILAR TO SUBROUTINE SYMBAC, BUT OPTIMIZED FOR DIAGONAL SITES
C     .. Scalar Arguments ..
      INTEGER MILEN,NBASIS,NT0A,NTXA
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,*)
      INTEGER IMAX(*),IMIN(*),MI(NBASIS,*)
C     ..
C     .. Local Scalars ..
      REAL*8 CI
      INTEGER I,IMAXB,IMINB,J,M,MIMAX,NM,NN,NVAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
      DO 50 NN = NT0A,NTXA
        DO 40 I = IMIN(NN),IMAX(NN)
          CI = CN(I,NN)
          IF (CI.EQ.0.D0) GO TO 40
C     LEFT INDEX OF GMN
          DO 30 NM = NT0A,NTXA
            IMINB = IMIN(NM)
            IMAXB = IMAX(NM)
            IF (IMINB.EQ.IMAXB) GO TO 30
C     WRITE(6,*) '121,',(MI(J,NM),J=IMINB,IMAXB)
            MIMAX = 49
            DO 20 M = 1,MIMAX
              NVAL = 0
              MILEN = MAX0(MILEN,IMAXB-IMINB+1)
              DO 10 J = IMINB,IMAXB
                IF (MI(J,NM).EQ.M) NVAL = NVAL + 1
   10         CONTINUE
              IF (NVAL.GT.1) WRITE (6,FMT=*) M,' OCCURS ',NVAL,' TIMES'
              IF (MILEN.GT.64) WRITE (6,FMT=*) MILEN,
     +            ' > 64 IS TOO LONG FOR SYMBCK'
              IF (NVAL.GT.1) STOP
              IF (MILEN.GT.64) STOP
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      RETURN
      END
