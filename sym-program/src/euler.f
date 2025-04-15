      SUBROUTINE EULER (KT,
     1          ALPHA,BETA,GAMMA,OMEGA,X,Y,Z)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N48=48)
C-----------------------------------------------------------------------
C THE SUBROUTINES CALCULATES THE EULER ANGLES ALPHA, BETA, GAMMA FOR THE
C SYMMETRY OPERATIONS DEFINED BY KT.
C (SEE MESSIAH: QUANTUM MECHANICS)
C-----------------------------------------------------------------------
      DIMENSION ALPHA(N48),BETA(N48),GAMMA(N48)
      PI=4.D0*DATAN(1.D0)
      OMEGA=DMOD(OMEGA,2.0D0)*PI
      ALPHA(KT)=0.0D0
      BETA(KT)=0.0D0
      GAMMA(KT)=0.0D0
      IF (OMEGA .EQ. 0.0D0) RETURN
      R=DSQRT(X*X+Y*Y+Z*Z)
      IF (R .EQ. 0.0D0) GO TO 30
      X=X/R
      Y=Y/R
      Z=Z/R
      P=DSQRT(X*X+Y*Y)
      IF (P .EQ. 0.0D0) GO TO 10
      IF (Z .EQ. 0.0D0) GO TO 20
      C=DCOS(0.5D0*OMEGA)
      S=DSIN(0.5D0*OMEGA)
      STH=P*S
      CTH=DSQRT(1.0D0-STH**2)
      SC=STH*CTH
      IF (DABS(SC) .LT. 1.D-3) GO TO 30
      CA=(C*Y+X*Z*S)*S/SC
      CG=(C*Y-X*Z*S)*S/SC
      SA=(-C*X+Y*Z*S)*S/SC
      SG=( C*X+Y*Z*S)*S/SC
      IF(DABS(SA).LE.1.D-5) GO TO 4
      ALPHA(KT)=2.0D0*DATAN2(SA,1.0D0+CA)/PI
4     ST=DABS(STH)
      IF(ST.LE.1.D0) GO TO 5
      IF(ST.GT.1.0001D0) GO TO 30
      IF(STH.LT.0.D0) STH=-1.D0
      IF(STH.GT.0.D0) STH=1.D0
5     BETA(KT)=2.0D0*DASIN(STH)/PI
      IF (DABS(SG) .LE.1.D-5) GO TO 6
      GAMMA(KT)=2.0D0*DATAN2(SG,1.0D0+CG)/PI
6     OMEGA=OMEGA/PI
      IF(DABS(SA).LE.1.D-5.AND.DABS(CA+1.).LE.1.D-5) ALPHA(KT)=1.D0
      IF(DABS(SG).LE.1.D-5.AND.DABS(CG+1.).LE.1.D-5) GAMMA(KT)=1.D0
      X=X*R
      Y=Y*R
      Z=Z*R
      RETURN
10    GAMMA(KT)=Z*OMEGA/PI
      OMEGA=OMEGA/PI
      Z=Z*R
      RETURN
20    CA=Y
      CG=Y
      SA=-X
      SG=X
      IF(DABS(SA).LE.1.D-5) GO TO 21
      ALPHA(KT)=2.0D0*DATAN2(SA,1.0D0+CA)/PI
21    BETA(KT)=OMEGA/PI
      IF(DABS(SG).LE.1.D-5) GO TO 22
      GAMMA(KT)=2.0D0*DATAN2(SG,1.0D0+CG)/PI
22    OMEGA=OMEGA/PI
      IF(DABS(SA).LE.1.D-5.AND.DABS(CA+1.).LE.1.D-5) ALPHA(KT)=1.D0
      IF(DABS(SG).LE.1.D-5.AND.DABS(CG+1.).LE.1.D-5) GAMMA(KT)=1.D0
      X=X*R
      Y=Y*R
      RETURN
30    WRITE (14,40) R,STH,CTH
40    FORMAT('  FAULTY  INPUT  DATA  ',3D14.7)
      STOP
      END
