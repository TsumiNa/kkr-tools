C **********************************************************************
      SUBROUTINE BESHAN(HL,JL,NL,Z,LMIN,LMAX)
      implicit none
C **********************************************************************
C  CALCULATES SPHERICAL BESSEL, HANKEL AND NEUMANN FUNCTIONS
C  FOR THE ORDERS 0 .LE. L .LE. LMAX.
C  FOR ARGUMENTS Z .LT. 1 THE TAYLOR EXPANSIONS OF JL AND NL ARE USED.
C  FOR ARGUMENTS Z .GE. 1 THE EXPLICIT EXPRESSIONS FOR HL(+), HL(-) ARE
C  USED.
C     .. Parameters ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX Z
      INTEGER LMAX,LMIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX TERMJ,TERMN,Z2,ZJ,ZN
      DOUBLE PRECISION RL,RN,RNM
      INTEGER L,M,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDABS,CDEXP
C     ..
      IF (CDABS(Z).LT.1.D0) THEN
        ZJ = 1.D0
        ZN = 1.D0
        Z2 = Z*Z
        DO 20 L = 0,LMAX
          RL = L + L
          TERMJ = -0.5D0/ (RL+3.D0)*Z2
          TERMN = 0.5D0/ (RL-1.D0)*Z2
          JL(L) = 1.D0
          NL(L) = 1.D0
          N = 1
          DO 10 WHILE ((CDABS(TERMJ) +CDABS(TERMN)).GE.1.0D-20)
             N = N + 1 
             JL(L) = JL(L) + TERMJ
             NL(L) = NL(L) + TERMN
             RN = N + N
             TERMJ = -TERMJ/ (RL+RN+1.D0)/RN*Z2
             TERMN = TERMN/ (RL-RN+1.D0)/RN*Z2
 10       CONTINUE
          JL(L) = JL(L)*ZJ
          NL(L) = -NL(L)*ZN/Z
          HL(L) = JL(L) + NL(L)*CI

          ZJ = ZJ*Z/ (RL+3.D0)
          ZN = ZN/Z* (RL+1.D0)
   20   CONTINUE

      ELSE
        DO 40 L = 0,LMAX
          HL(L) = 0.D0
          NL(L) = 0.D0
          RNM = 1.D0
          DO 30 M = 0,L
            HL(L) = HL(L) + RNM/ (-CI* (Z+Z))**M
            NL(L) = NL(L) + RNM/ (CI* (Z+Z))**M
            RNM = RNM* (L*L+L-M*M-M)/ (M+1.D0)
   30     CONTINUE
          HL(L) = HL(L)* (-CI)**L*CDEXP(CI*Z)/ (CI*Z)   ! HL +
          NL(L) = NL(L)*CI**L*CDEXP(-CI*Z)/ (-CI*Z)     ! HL -
          JL(L) = (HL(L)+NL(L))*0.5D0
          NL(L) = (HL(L)-JL(L))/CI
   40   CONTINUE
      END IF

      RETURN

      END
