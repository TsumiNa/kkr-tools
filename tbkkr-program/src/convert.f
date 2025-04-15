C ************************************************************************
      SUBROUTINE CONVERT(PL1,PL2,BLOCK,GIN)
C ************************************************************************
      implicit none

c     .. parameters ..
      include 'inc.fi'

      INTEGER LMAXSQ
      PARAMETER (LMAXSQ=(LMAXD+1)**2)
      INTEGER ALMD,NDIM
      PARAMETER (ALMD= NAEZD*LMAXSQ,NDIM = NPRINCD*LMAXSQ)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
      DOUBLE COMPLEX GIN(ALMD,ALMD),BLOCK(NDIM,NDIM)
      INTEGER IP1,II1,IL1,LDI1,IP2,II2,IL2,LDI2,LM1,LM2
      INTEGER PL1,PL2

      DO IP1 = 1,NPRINCD
      DO IP2 = 1,NPRINCD
         II1 = (PL1-1)*NPRINCD+IP1
         II2 = (PL2-1)*NPRINCD+IP2
         DO LM1 = 1,LMAXSQ
         DO LM2 = 1,LMAXSQ
            LDI1 = LMAXSQ*(IP1-1)+LM1
            IL1 = LMAXSQ*(II1-1)+LM1
            LDI2 = LMAXSQ*(IP2-1)+LM2
            IL2 = LMAXSQ*(II2-1)+LM2
            GIN(IL1,IL2) = BLOCK(LDI1,LDI2)
         ENDDO
         ENDDO
      ENDDO
      ENDDO


      RETURN

      END
