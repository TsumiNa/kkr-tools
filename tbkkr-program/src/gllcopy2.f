c 03.06.97 ***************************************************************
      SUBROUTINE GLLCOPY2(GLLKE0,G,D1,D2)
      implicit none
c ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER NDIM
      PARAMETER (NDIM=LMAXSQ*NPRINCD)
c
c     .. arguments
      DOUBLE COMPLEX GLLKE0(NDIM,*),
     +               G(LMAXSQ,*)
      INTEGER D1,D2
C
c     .. Local
      INTEGER D1LM,D2LM,LM1,LM2
c     .. external
      EXTERNAL RCSTOP
c ------------------------------------------------------------------------
      IF (D1.LT.1 .OR. D1.GT.NPRINCD .OR. 
     +     D2.LT.1 .OR. D2.GT.NPRINCD ) THEN
        WRITE(6,*) 'D1, D2 : ',D1,D2
        CALL RCSTOP('GLLCOPY2')
      END IF

      D1LM = (D1-1)*LMAXSQ
      D2LM = (D2-1)*LMAXSQ
      
      DO 10 LM1 = 1,LMAXSQ
        DO 20 LM2 = 1,LMAXSQ
          G(LM1,LM2) = GLLKE0(D2LM+LM2,D1LM+LM1)
 20     END DO
 10   END DO

      RETURN
      END
