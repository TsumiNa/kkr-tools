c ************************************************************************
      SUBROUTINE EINIT(MAX)
      implicit none
c ************************************************************************
c     p.zahn, april 96
c ------------------------------------------------------------------------
      include 'inc.fi'
      INTEGER MAX

c     .. locals
      INTEGER N,OPT,IP,LM1,LM2

c     .. arrays in common
      DOUBLE PRECISION EE1(0:M2D,KPOIBZ),EE2(0:M2D,KPOIBZ),
     +                 WW1(0:M2D,KPOIBZ),WW2(0:M2D,KPOIBZ)
      INTEGER II1(0:M2D,KPOIBZ),II2(0:M2D,KPOIBZ)
      COMMON /EE/ EE1,EE2,WW1,WW2,II1,II2 ! index arrays for eigenvalue
                                          ! analysis

      DOUBLE PRECISION ZERO,TH
      DATA ZERO / 0.0D0 /
     +     TH   / 1.0D3 /

      SAVE
c ------------------------------------------------------------------------
      DO 10 N = 0,M2D
        DO 20 IP=1,MAX
          EE1(N,IP)  = TH
          EE2(N,IP)  = TH
          WW1(N,IP) = ZERO
          WW2(N,IP) = ZERO
          II1(N,IP) = 0
          II2(N,IP) = 0
 20     END DO
 10   END DO

      RETURN
      END 
