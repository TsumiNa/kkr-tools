      SUBROUTINE WRITLOWL(E,EK,DF,GHOST,NDIMNP,IGHWRIT,P4TO3,
     +                    NDIMLOW)
      IMPLICIT NONE
CO--->WRITE OUT TRANSFORMED GREENS-FUNCTION, AFTER
CO--->TRANSFORMATION OR FIRST DYSON STEP IN CASE OF
CO--->LATTICE RELAXATION

C     .. Parameters ..
      INTEGER NSEC
      PARAMETER (nsec=689)
C     ..
c
c---> read sym. host green's functions from disc
c
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E,EK
      INTEGER NDIMNP,IGHWRIT,IHANDLE,NDIMLOW
      INTEGER P4TO3(NSEC)
C     ..
C     .. Array Arguments ..
      COMPLEX*16 GHOST(NSEC,NSEC)
      COMPLEX*8 GLOW(NSEC,NSEC)
CFXDR COMPLEX*8 GLIN(NSEC*NSEC)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I2,ILIN,NLIN,IN1,IN2
C     ..
CFXDR      NLIN = NDIMNP*(NDIMNP+1)/2
      DO I1=1,NDIMNP
         DO I2=1,NDIMNP
            IN1=P4TO3(I1)
            IN2=P4TO3(I2)
            IF ((IN1 .GT. NDIMLOW) .OR. (IN2 .GT. NDIMLOW))
     +      STOP 'NDIMLOW'
            IF ((IN1 .EQ. 0) .OR. (IN2 .EQ. 0)) GOTO 500
            GLOW(IN1,IN2)=GHOST(I1,I2)
 500        CONTINUE
         END DO
      END DO
      WRITE (IGHWRIT) E,EK,DF, ((GLOW(I2,I1),I2=1,I1),I1=1,NDIMLOW)
      END
