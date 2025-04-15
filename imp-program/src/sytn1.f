      SUBROUTINE SYTN1(DTMTRX,INS,LMAX,NATREF,NSTART,NEND,TMATLL,IRMIN,
     +                 NQ,NDIM,NREP,TLLMAT)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     the delta t - matrix is symmetrized
c     (transformed into the irreducible representation)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NSEC,NREPD
      PARAMETER (nsec=689,NREPD=4)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER ICJD
      PARAMETER (ICJD=93139)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER INS,LMAX,NATREF,NEND,NQ,NREP,NSTART
      LOGICAL TLLMAT
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DTMTRX(NSEC,NSEC),TMATLL(*)
      INTEGER IRMIN(*),NDIM(*)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I2,ICJ,IRMIN1,IT,J,LA,LB,LMMAX,NP
      LOGICAL LJ
C     ..
C     .. Local Arrays ..
      REAL*8 JCOEFF(ICJD)
      INTEGER I1J(ICJD),I2J(ICJD),IG(ICJD),LN(NSEC,NREPD),NREPJ(ICJD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA LJ/.true./
C     ..

      LMMAX = (LMAX+1)* (LMAX+1)
c
c
c---> in the first time : read symmetrization coeffients
c
      IF (LJ) THEN
        IRMIN1 = 0
        DO 10 I1 = NSTART,NEND
          IF (IRMIN1.LT.IRMIN(I1)) IRMIN1 = IRMIN(I1)
   10   CONTINUE
c
        READ (35,FMT=9001) ICJ
c
c---> check dimensions
c
        WRITE (6,FMT=9020) ICJ,ICJD
        IF (ICJ.GT.ICJD) STOP
c
        READ (35,FMT=9010) (JCOEFF(J),J=1,ICJ)
        READ (35,FMT=9000) (I1J(J),J=1,ICJ)
        READ (35,FMT=9000) (I2J(J),J=1,ICJ)
        READ (35,FMT=9000) (NREPJ(J),J=1,ICJ)
        READ (35,FMT=9000) (IG(J),J=1,ICJ)
        DO 20 NP = 1,NREP
          READ (35,FMT=9000) (LN(J,NP),J=1,NDIM(NP))
   20   CONTINUE
        LJ = .false.

      END IF

c
c---> initialize dtmtrx
c
      CALL CINIT(NSEC*NSEC,DTMTRX)
c
      DO 30 J = 1,ICJ
        I1 = I1J(J)
        I2 = I2J(J)
        NP = NREPJ(J)
        IF (NP.NE.NQ) GO TO 30
        IT = NATREF*LMMAXD**2 + IG(J)
        DTMTRX(I1,I2) = DTMTRX(I1,I2) + JCOEFF(J)*TMATLL(IT)
   30 CONTINUE
      IF (TLLMAT) THEN
        DO 50 I1 = 2,NDIM(NQ)
          DO 40 I2 = 1,I1 - 1
            DTMTRX(I2,I1) = DTMTRX(I1,I2)
   40     CONTINUE
   50   CONTINUE
      ELSE
        DO 70 I1 = 2,NDIM(NQ)
          LA = LN(I1,NQ)
          DO 60 I2 = 1,I1 - 1
            LB = LN(I2,NQ)
c           DTMTRX(I2,I1) = DTMTRX(I1,I2)* (-1)** (LA-LB)
            DTMTRX(I2,I1) = DTMTRX(I1,I2)* (-1)** (LA+LB)
   60     CONTINUE
   70   CONTINUE
      END IF
c
      RETURN


 9000 FORMAT (16i5)
 9001 FORMAT (16i8)
 9010 FORMAT (5d16.9)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'icj    : (',i5,',',i5,')',/)
      END
