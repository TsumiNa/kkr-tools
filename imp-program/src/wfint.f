      SUBROUTINE WFINT(QNS,CDER,DDER,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     +                 NSRA,LKONV,irmd,irmind)
      Implicit None
c-----------------------------------------------------------------------
c      determines the integrands CDER, DDER or ADER, BDER in the
c        integral equations for the non-spherical wavefunctions from
c        the non-spherical contributions of the potential vinsPLL.
c
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
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2)
      REAL*8 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Local Scalars ..
      INTEGER IR,LM1,LM2,LMMKONV
C     ..
C     .. External Subroutines ..
      EXTERNAL DGEMM
C     ..
C     .. Local Arrays ..
      REAL*8 QNSI(LMMAXD,LMMAXD),QNSR(LMMAXD,LMMAXD),
     +     VTQNSI(LMMAXD,LMMAXD),VTQNSR(LMMAXD,LMMAXD),
     $     vzwll(LMMAXD,LMMAXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX,DIMAG,DREAL
C     ..
c
      IF (LKONV .NE. LMAXD) THEN
         LMMKONV=(LKONV+1)**2
         VZWLL=0.0D0
         DO IR=IRMIN,IRC1
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
      DO 90 IR = IRMIN,IRC1
        DO 20 LM2 = 1,LMMAXD
          DO 10 LM1 = 1,LMMAXD
            QNSR(LM1,LM2) = DREAL(QNS(LM1,LM2,IR,1))
            QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,1))
   10     CONTINUE
   20   CONTINUE
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
        CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +             LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
        DO 40 LM1 = 1,LMMAXD
          DO 30 LM2 = 1,LMMAXD
            CDER(LM1,LM2,IR) = QZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
            DDER(LM1,LM2,IR) = PZEKDR(LM1,IR,1)*
     +                         DCMPLX(VTQNSR(LM1,LM2),VTQNSI(LM1,LM2))
   30     CONTINUE
   40   CONTINUE
        IF (NSRA.EQ.2) THEN
          DO 60 LM2 = 1,LMMAXD
            DO 50 LM1 = 1,LMMAXD
              QNSR(LM1,LM2) = DREAL(QNS(LM1,LM2,IR,2))
              QNSI(LM1,LM2) = DIMAG(QNS(LM1,LM2,IR,2))
   50       CONTINUE
   60     CONTINUE
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSR,LMMAXD,0.D0,VTQNSR,LMMAXD)
          CALL DGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,1.D0,VNSPLL(1,1,IR),
     +               LMMAXD,QNSI,LMMAXD,0.D0,VTQNSI,LMMAXD)
          DO 80 LM2 = 1,LMMAXD
            DO 70 LM1 = 1,LMMAXD
              CDER(LM1,LM2,IR) = CDER(LM1,LM2,IR) +
     +                           QZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
              DDER(LM1,LM2,IR) = DDER(LM1,LM2,IR) +
     +                           PZEKDR(LM1,IR,2)*DCMPLX(VTQNSR(LM1,
     +                           LM2),VTQNSI(LM1,LM2))
   70       CONTINUE
   80     CONTINUE
        END IF

   90 CONTINUE
      END
