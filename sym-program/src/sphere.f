      SUBROUTINE SPHERE(LMAX,YR,WTYR,RIJ,IJEND,IJD)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     generate an angular mesh and spherical harmonics at those
c     mesh points. For an angular integration the weights are ge-
c     rated .
c
c     R. Zeller      Feb. 1996
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IJD,IJEND,LMAX
C     ..
C     .. Local Scalars ..
      REAL*8 PI,R,R1,R2,R3
      INTEGER IJ,LM1,IHAND
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Array Arguments ..
      REAL*8 RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*),COSX(IJD),FAI(IJD)
      REAL*8 WEI
C     ..
C     .. Local Arrays ..
      REAL*8 W(1000),Y(1000)
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
C     ..
C      CALL FXDROPN('lebedev.fxdr',
C     +       'DECODE',IHAND)
C      CALL FXDRINT(IHAND,IJEND,1)
      OPEN (89,FILE='lebedev',FORM='formatted')
      READ (89,*) IJEND
      WRITE(6,*) IJEND       
      IF (IJEND.GT.IJD) STOP 'SPHERE'
      IF (IJEND.GT.1000) STOP 'SPHERE'
c
c
      DO 30 IJ = 1,IJEND
c      CALL FXDRDBL(IHAND,R1,1)
c      CALL FXDRDBL(IHAND,R2,1)
c      CALL FXDRDBL(IHAND,R3,1)
c      CALL FXDRDBL(IHAND,WEI,1)
c      W(IJ)=WEI
c      CALL FXDRDBL(IHAND,COSX(IJ),1)
c      CALL FXDRDBL(IHAND,FAI(IJ),1)
        READ (89,9001) R1,R2,R3,W(IJ)
        RIJ(IJ,1) = R1
        RIJ(IJ,2) = R2
        RIJ(IJ,3) = R3
        CALL YMY(R1,R2,R3,R,Y,LMAX)
        DO 10 LM1 = 1, (LMAX+1)**2
          YR(IJ,LM1) = Y(LM1)
   10   CONTINUE
c
c---> multiply the spherical harmonics with the weights
c
        DO 20 LM1 = 1, (LMAX+1)**2
          WTYR(IJ,LM1) = YR(IJ,LM1)*W(IJ)*PI*4.D0
   20   CONTINUE
   30 CONTINUE
c      CALL FXDRCLS(IHAND)
 9001 format(4D20.12)
      END
