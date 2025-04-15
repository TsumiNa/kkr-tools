c 20.07.96 ***************************************************************
      SUBROUTINE GLLCOPY1(GLLKE0,G,D1,D2)
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
      EXTERNAL RCSTOP,ZCOPY
c ------------------------------------------------------------------------
      IF (D1.LT.1 .OR. D1.GT.NPRINCD .OR. 
     +     D2.LT.1 .OR. D2.GT.NPRINCD ) THEN
        WRITE(6,*) 'D1, D2 : ',D1,D2
        CALL RCSTOP('GLLCOPY1')
      END IF

      D1LM = (D1-1)*LMAXSQ+1
      D2LM = (D2-1)*LMAXSQ
      
      IF (NPRINCD.GT.1) THEN
        DO 10 LM2 = 1,LMAXSQ
          CALL ZCOPY(LMAXSQ,GLLKE0(D1LM,D2LM+LM2),1,G(1,LM2),1)
 10     END DO
      ELSE
        CALL ZCOPY(LMAXSQ*LMAXSQ,GLLKE0(1,1),1,G(1,1),1)
      END IF

      RETURN
      END
