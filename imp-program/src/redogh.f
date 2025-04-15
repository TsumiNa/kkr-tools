      SUBROUTINE REDOGH(E,EK,DF,GHOST,NDIMNP,GTEMP,IHANDLE)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NSEC
      PARAMETER (nsec=689)
C     ..
c
c---> read sym. host green's functions from disc
c
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E,EK
      INTEGER NDIMNP,IHANDLE
C     ..
C     .. Array Arguments ..
      COMPLEX*16 GHOST(NSEC,NSEC)
      COMPLEX*16 GTEMP(NSEC,NSEC),GLIN(NSEC*NSEC)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I2,ILIN,NLIN
C     ..
      NLIN = NDIMNP*(NDIMNP+1)/2
      READ (IHANDLE) E,EK,DF, ((GTEMP(I2,I1),I2=1,I1),I1=1,NDIMNP)
      DO 20 I1 = 1,NDIMNP
        DO 10 I2 = 1,I1
          GHOST(I2,I1) = GTEMP(I2,I1)
          GHOST(I1,I2) = GTEMP(I2,I1)
   10   CONTINUE
   20 CONTINUE
      ILIN = 1
      END
