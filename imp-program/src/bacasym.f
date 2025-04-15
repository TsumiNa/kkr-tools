
      SUBROUTINE BACASYM(DTB,GMAT,LMAXSQ,NREP,NP,KSYMMAT,
     +                   SOCCUP,QOCCUP,ISPIN)
      IMPLICIT NONE
CASYM--------------------------------------------
CASYM    T H I S   I S  F O R  A S Y M    
CASYM--------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD)
      INTEGER NSEC,NREPD,NSPIND
      PARAMETER (nsec=689,NREPD=4,NSPIND=2)
      INTEGER LMAXD,LMX,LMMAXD
      PARAMETER (lmaxd=4,LMX=LMAXD+1,LMMAXD=LMX**2)
      INTEGER JCOEFF
      PARAMETER (JCOEFF=187261)
C     .. 
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,NP,NREP,KSYMMAT
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DTB(LMMAXD,LMMAXD,*),GMAT(NSEC,NSEC)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO
      INTEGER I1,I2,I3,ICG,J,J1,J2,NR,ISPIN
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Local Arrays ..
      REAL*8 GCOENP(JCOEFF),QOCCUP(NSPIND,NREPD)
      INTEGER I1NP(JCOEFF),I2NP(JCOEFF),I3NP(JCOEFF),IANP(JCOEFF),
     +        IBNP(JCOEFF),SOCCUP(NSPIND,NREPD)
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..

      IF (NP.EQ.1) THEN
        REWIND 71
        DO 30 I3 = 1,NTPERD
          DO 20 I1 = 1,LMAXSQ
            DO 10 I2 = 1,LMAXSQ
              DTB(I1,I2,I3) = CZERO
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
      END IF
      READ (71) NR,ICG, (GCOENP(J),J=1,ICG), (I1NP(J),J=1,ICG),
     +  (I2NP(J),J=1,ICG), (I3NP(J),J=1,ICG), (IANP(J),J=1,ICG),
     +  (IBNP(J),J=1,ICG)
      IF (ICG.GT.JCOEFF) STOP 'JCOEFF - BACASYM'
      IF (KSYMMAT .EQ. NP) THEN

        IF (SOCCUP(ISPIN,NP).EQ.1) THEN
        DO 40 J = 1,ICG
        I1 = I1NP(J)
        I2 = I2NP(J)
        I3 = I3NP(J)
        J1 = IANP(J)
        J2 = IBNP(J)
        DTB(I1,I2,I3) = DTB(I1,I2,I3) + QOCCUP(ISPIN,NP)*GCOENP(J)
     +                  *GMAT(J1,J2)
   40   CONTINUE
        ELSE
        DO  J = 1,ICG
          I1 = I1NP(J)
          I2 = I2NP(J)
          I3 = I3NP(J)
          J1 = IANP(J)
          J2 = IBNP(J)
          DTB(I1,I2,I3) = DTB(I1,I2,I3) + GCOENP(J)*GMAT(J1,J2)
	END DO
        END IF  

      END IF
      IF (NP.EQ.NREP) THEN
        DO 70 I3 = 1,NTPERD
          DO 60 I2 = 1,LMAXSQ
            DO 50 I1 = 1,I2 - 1
              DTB(I2,I1,I3) = DTB(I1,I2,I3)
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
      END IF

      END
