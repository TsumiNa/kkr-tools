      SUBROUTINE WFINT0(CDER,DDER,QZLM,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     +                  NSRA,LKONV,irmd,irmind)
      Implicit None
c-----------------------------------------------------------------------
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c        (This subroutine is used in zeroth order Born approximation,
c         otherwise subroutine WFINT must be used)
c      R. Zeller      Aug. 1994
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER LMAXD
      PARAMETER (lmaxd=4)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IRC1,IRMIN,NSRA,LKONV,irmd,irmind
C     ..
C     .. Array Arguments ..
      COMPLEX*16 CDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2)
      REAL*8 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Local Arrays ..
      REAL*8 vzwll(LMMAXD,LMMAXD)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 V1
      INTEGER IR,LM1,LM2,LMMKONV
C     ..
c
      IF (LKONV .NE. LMAXD) THEN
        LMMKONV=(LKONV+1)**2
        DO IR=IRMIN,IRC1
          VZWLL=0.0D0
          DO LM1=1,LMMKONV
            DO LM2=1,LMMKONV
              VZWLL(LM1,LM2)=VNSPLL(LM1,LM2,IR)
            END DO
          END DO
          DO LM1=1,LMMAXD
            DO LM2=1,LMMAXD
              VNSPLL(LM1,LM2,IR)=VZWLL(LM1,LM2)
            END DO
          END DO
        END DO
      END IF
c
      DO 50 IR = IRMIN,IRC1
        DO 20 LM2 = 1,LMMAXD
          DO 10 LM1 = 1,LMMAXD
            V1 = VNSPLL(LM1,LM2,IR)*QZLM(LM2,IR,1)
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1)*V1
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1)*V1
   10     CONTINUE
   20   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 40 LM2 = 1,LMMAXD
            DO 30 LM1 = 1,LMMAXD
              V1 = VNSPLL(LM1,LM2,IR)*QZLM(LM2,IR,2)
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) + QZEKDR(LM1,IR,2)*V1
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) + PZEKDR(LM1,IR,2)*V1
   30       CONTINUE
   40     CONTINUE
        END IF

   50 CONTINUE
      END
