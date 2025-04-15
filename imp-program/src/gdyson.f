      SUBROUTINE GDYSON(DTMTRX,GMAT,DET,NDIMNP,LJUSTD)
      IMPLICIT NONE
c
c---> solve the dyson equation to get disturbed green's functions
c
C     .. Parameters ..
      INTEGER NSEC
      PARAMETER (nsec=689)
      COMPLEX*16 CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
C     ..
C     .. Local Scalars ..
      INTEGER I,IB,INFO
C     ..
C     .. Local Arrays ..
      COMPLEX*16 BSYM(NSEC,NSEC),GSYM(NSEC,NSEC)
      INTEGER IPVT(NSEC)
C     ..
C     .. External Subroutines ..
      EXTERNAL ZCOPY,ZGEMM,ZGETRF,ZGETRS
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DET
      INTEGER NDIMNP
      LOGICAL LJUSTD
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DTMTRX(NSEC,NSEC),GMAT(NSEC,NSEC)
C     ..
      CALL ZGEMM('N','N',NDIMNP,NDIMNP,NDIMNP,-CONE,GMAT,NSEC,DTMTRX,
     +           NSEC,CZERO,GSYM,NSEC)
      CALL ZCOPY(NSEC*NDIMNP,GMAT,1,BSYM,1)
c
      DO 10 I = 1,NDIMNP
        GSYM(I,I) = CONE + GSYM(I,I)
   10 CONTINUE
c
c---> solve the system of linear equations and calculate determinant
c
      CALL ZGETRF(NDIMNP,NDIMNP,GSYM,NSEC,IPVT,INFO)
      DET = CONE
      DO 20 I = 1,NDIMNP
        IF (IPVT(I).NE.I) DET = -DET
        DET = GSYM(I,I)*DET
   20 CONTINUE
      IF (LJUSTD) RETURN
      CALL ZGETRS('N',NDIMNP,NDIMNP,GSYM,NSEC,IPVT,BSYM,NSEC,INFO)
      DO 40 IB = 1,NDIMNP
        DO 30 I = 1,NDIMNP
          GMAT(I,IB) = BSYM(I,IB)
   30   CONTINUE

   40 CONTINUE

      RETURN

      END
