      SUBROUTINE SYTMAT(LMAX,TMATLL,IDND,ND,YR,WTYR,RIJ,IJEND)
      IMPLICIT NONE
C     .. Parameters ..

      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER LPOTD
      PARAMETER (lpotd=8)
      INTEGER LMPOTD,IJD
      PARAMETER (LMPOTD= (LPOTD+1)**2,IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,LMAX
C     ..
C     .. Array Arguments ..
      COMPLEX*16 TMATLL(LMMAXD,LMMAXD)
      REAL*8 RIJ(IJD,3),WTYR(IJD,LMPOTD),YR(IJD,LMPOTD)
      INTEGER IDND(6),ND(48,3,3)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO
      INTEGER II,IN,LM1,LM2,LM3,LMMAX
      LOGICAL LJ
C     ..
C     .. Local Arrays ..
      COMPLEX*16 AQ(LMMAXD,LMMAXD),TMSAVE(LMMAXD,LMMAXD)
      REAL*8 DROT(LMMAXD,LMMAXD,40:48)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,ROTTNS,ZAXPY,ZCOPY,ZSCAL
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA LJ/.true./
      DATA CZERO/ (0.D0,0.D0)/
C     ..

      LMMAX = (LMAX+1)* (LMAX+1)
      IF (LJ) THEN
c
c---> only once: determine rotation matrix for t  - matrices
c
        CALL ROTTNS(DROT(1,1,40),40,LMAX,ND,YR,WTYR,RIJ,IJEND)
        CALL ROTTNS(DROT(1,1,41),41,LMAX,ND,YR,WTYR,RIJ,IJEND)
        CALL ROTTNS(DROT(1,1,42),42,LMAX,ND,YR,WTYR,RIJ,IJEND)
        CALL ROTTNS(DROT(1,1,46),46,LMAX,ND,YR,WTYR,RIJ,IJEND)
        CALL ROTTNS(DROT(1,1,47),47,LMAX,ND,YR,WTYR,RIJ,IJEND)
        CALL ROTTNS(DROT(1,1,48),48,LMAX,ND,YR,WTYR,RIJ,IJEND)
c
        LJ = .false.

      END IF

      DO 70 II = 1,6
        IN = IDND(II)
        IF (IN.EQ.1) GO TO 70
        CALL CINIT(LMMAXD**2,AQ)
        CALL ZSCAL(LMMAXD**2, (0.5D0,0.0D0),TMATLL,1)
        CALL ZCOPY(LMMAXD**2,TMATLL,1,TMSAVE,1)
c
        DO 30 LM2 = 1,LMMAX
          DO 20 LM3 = 1,LMMAX
            DO 10 LM1 = 1,LMMAXD
              AQ(LM1,LM2) = AQ(LM1,LM2) +
     +                      DROT(LM1,LM3,IN)*TMATLL(LM3,LM2)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
c
c---> initialize
c
        CALL CINIT(LMMAXD**2,TMATLL)
c
        DO 60 LM2 = 1,LMMAX
          DO 50 LM3 = 1,LMMAX
            DO 40 LM1 = 1,LMMAXD
              TMATLL(LM1,LM2) = TMATLL(LM1,LM2) +
     +                          AQ(LM1,LM3)*DROT(LM2,LM3,IN)
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
        CALL ZAXPY(LMMAXD**2, (1.D0,0.D0),TMSAVE,1,TMATLL,1)
   70 CONTINUE
      RETURN

      END
