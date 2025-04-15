      SUBROUTINE ROTTNS(C,K,LMAX,ND,YR,WTYR,RIJ,IJEND)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     rottmt computes the integrals of y(l,m)*rotated{y(l,m')} up
c     to lmax . this matrix is needed for the rotation of the
c     non spherical t - matrices
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER LMAXD,LMX,LPOTD
      PARAMETER (lmaxd=4,LMX=LMAXD+1,lpotd=8)
      INTEGER LMPOTD,IJD
      PARAMETER (LMPOTD= (LPOTD+1)**2,IJD=434)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX*LMX)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,K,LMAX
C     ..
C     .. Array Arguments ..
      REAL*8 C(LMMAXD,LMMAXD),RIJ(IJD,3),WTYR(IJD,LMPOTD),
     +                 YR(IJD,LMPOTD)
      INTEGER ND(48,3,3)
C     ..
C     .. Local Scalars ..
      REAL*8 ROTR,ROTR1,ROTR2,ROTR3
      INTEGER IJ,L,LM1,LM2,LMMAX,M1,M2
C     ..
C     .. Local Arrays ..
      REAL*8 YRROT(LMMAXD,IJD)
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      LMMAX = (LMAX+1)* (LMAX+1)
c---> generate the rotated spherical harmonics on the angular mesh
c
      DO 10 IJ = 1,IJEND
        ROTR1 = ND(K,1,1)*RIJ(IJ,1) + ND(K,1,2)*RIJ(IJ,2) +
     +          ND(K,1,3)*RIJ(IJ,3)
        ROTR2 = ND(K,2,1)*RIJ(IJ,1) + ND(K,2,2)*RIJ(IJ,2) +
     +          ND(K,2,3)*RIJ(IJ,3)
        ROTR3 = ND(K,3,1)*RIJ(IJ,1) + ND(K,3,2)*RIJ(IJ,2) +
     +          ND(K,3,3)*RIJ(IJ,3)
        CALL YMY(ROTR1,ROTR2,ROTR3,ROTR,YRROT(1,IJ),LMAX)
   10 CONTINUE
c
      DO 30 LM2 = 1,LMMAX
        DO 20 LM1 = 1,LMMAX
          C(LM1,LM2) = 0.0D0
   20   CONTINUE
   30 CONTINUE
c
c---> now integrate
c
      DO 70 L = 0,LMAX
        DO 60 M1 = -L,L
          DO 50 M2 = -L,L
            LM1 = L* (L+1) + M1 + 1
            LM2 = L* (L+1) + M2 + 1
            DO 40 IJ = 1,IJEND
              C(LM1,LM2) = C(LM1,LM2) + WTYR(IJ,LM1)*YRROT(LM2,IJ)
   40       CONTINUE
            IF (ABS(C(LM1,LM2)).LE.1.0D-10) C(LM1,LM2) = 0.0D0
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE

      END
