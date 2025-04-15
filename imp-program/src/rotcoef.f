      SUBROUTINE ROTCOEF(C,IATOM,IOPER,LMAX,ND,YR,WTYR,RIJ,IJEND)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     rotcoef computes the integrals of y(l,m)*rotated{y(l,m')} up
c     to lmax . these value are needed for calculating the charge den-
c     sity of a given atom -iatom- when the charge density of the
c     representive is known .
c     the spherical harmonics and the spherical harmonics times the
c     weights generated on the unit sphere , the points and the
c     weights for the gauss-legendre integration are calculated in
c     the subroutine sphere and stored in the common block rsphere .
c     therefore only the rotated spherical harmonics have to calcu-
c     lated .
c     attention : for a better vectorization combined indices
c                 are introduced .
c
c                                  b.drittler    may 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER LPOTD
      PARAMETER (lpotd=8)
      INTEGER LMPOTD,N,IJD
c     PARAMETER (LMPOTD= (LPOTD+1)**2,N=2* (LPOTD+1),IJD=2*N**2)
      PARAMETER (LMPOTD= (LPOTD+1)**2,N=2* (LPOTD+1),IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER IATOM,IJEND,LMAX
C     ..
C     .. Array Arguments ..
      REAL*8 C(0:LPOTD,-LPOTD:LPOTD,-LPOTD:LPOTD),RIJ(IJD,3),
     +       WTYR(IJD,LMPOTD),YR(IJD,LMPOTD)
      INTEGER IOPER(*),ND(48,3,3)
C     ..
C     .. Local Scalars ..
      REAL*8 ROTR,ROTR1,ROTR2,ROTR3
      INTEGER IJ,IL,K,LM1,LM2,M1,M2
C     ..
C     .. Local Arrays ..
      REAL*8 YRROT(LMPOTD,IJD)
C     ..
C     .. External Subroutines ..
      EXTERNAL YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
c
c---> determine the number of rotation matrix which transforms
c     iatom into its representive .
c
      K = IOPER(IATOM)
c
c---> generate the rotated spherical harmonics on the angular mesh
c      (mesh points are generated in deck sphere)
c
      DO 10 IJ = 1,IJEND
        ROTR1 = ND(K,1,1)*RIJ(IJ,1) + ND(K,1,2)*RIJ(IJ,2) +
     +          ND(K,1,3)*RIJ(IJ,3)
        ROTR2 = ND(K,2,1)*RIJ(IJ,1) + ND(K,2,2)*RIJ(IJ,2) +
     +          ND(K,2,3)*RIJ(IJ,3)
        ROTR3 = ND(K,3,1)*RIJ(IJ,1) + ND(K,3,2)*RIJ(IJ,2) +
     +          ND(K,3,3)*RIJ(IJ,3)
c
        CALL YMY(ROTR1,ROTR2,ROTR3,ROTR,YRROT(1,IJ),LMAX)
   10 CONTINUE
c
c
c---> now integrate
c
      DO 50 IL = 0,LMAX
        DO 40 M1 = -IL,IL
          DO 30 M2 = -IL,IL
            C(IL,M1,M2) = 0.D0
            LM1 = IL* (IL+1) + M1 + 1
            LM2 = IL* (IL+1) + M2 + 1
            DO 20 IJ = 1,IJEND
              C(IL,M1,M2) = C(IL,M1,M2) + WTYR(IJ,LM1)*YRROT(LM2,IJ)
   20       CONTINUE
            IF (ABS(C(IL,M1,M2)).LE.1.0d-10) C(IL,M1,M2) = 0.0D0
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE

      END
