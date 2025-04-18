
C ************************************************************************
      SUBROUTINE GREFSY(GTMAT,GMAT,NDIM)
      implicit none
C ************************************************************************
C
C---> SOLVE THE DYSON EQUATION TO GET REFERENCE GREEN FUNCTION
C
C     .. PARAMETERS ..
      INCLUDE 'inc.fi'
      INCLUDE 'inc.cls'
C
      INTEGER LMAX,NATOMD
      PARAMETER (LMAX=LMAXD,NATOMD=NACLSD)
      INTEGER LMAXSQ,NGD
      PARAMETER (LMAXSQ= (LMAX+1)**2,NGD=LMAXSQ*NATOMD)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. LOCAL SCALARS ..
      INTEGER I,INFO
C     ..
C     .. LOCAL ARRAYS ..
      INTEGER IPVT(NGD)
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL ZGETRF,ZGETRS
C     ..
C     .. SAVE STATEMENT ..
      SAVE
C     ..
C     .. SCALAR ARGUMENTS ..
      INTEGER NDIM
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX GMAT(NGD,LMAXSQ),GTMAT(NGD,NGD)
C     ..
C
      DO 10 I = 1,NDIM
        GTMAT(I,I) = CONE + GTMAT(I,I) ! GTMAT= 1 - G * T
   10 CONTINUE
C
C---> SOLVE THE SYSTEM OF LINEAR EQUATIONS
C
      CALL ZGETRF(NDIM,NDIM,GTMAT,NGD,IPVT,INFO)
      CALL ZGETRS('N',NDIM,LMAXSQ,GTMAT,NGD,IPVT,GMAT,NGD,INFO)
      RETURN

      END
