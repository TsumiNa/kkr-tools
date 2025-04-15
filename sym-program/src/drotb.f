      FUNCTION DROTB(L,MP,M,BETA)
C-----------------------------------------------------------------------
C  FUNCTION TO CALCULATE
C
C   L
C  D   (BETA)={(L-M)!(L+M)!(L-M')!(L+M')!}**(1/2)
C   M'M
C          (-1)**K*(COS(BETA/2))**(2L+M-M'-2K)*(-SIN(BETA/2))**(M'-M+2K)
C    *SUM  -------------------------------------------------------------
C        K                (L+M-K)!(L-M'-K)!K!(K+M'-M)!
C
C WITH MAX(M-M',0).LE.K.LE.MIN(L-M',L+M)
C
C AS GIVEN IN EQ.(4.13) OF "ELEMENTARY THEORY OF ANGULAR MOMENTA"
C BY M.E. ROSE, JOHN WILEY AND SONS, NEW YORK.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C     D COEFFICIENT.    EQ.4.13, ROSE.
      DIMENSION NF(4)
      EQUIVALENCE (N1,NF(1)),(N2,NF(2)),(N3,NF(3)),(N4,NF(4))
      DATA L0,M0,MP0/-1,2*1/
      IF(L.NE.L0) GO TO 2
      IF((IABS(M ).EQ.IABS(M0).AND.IABS(MP).EQ.IABS(MP0)).OR.
     1   (IABS(MP).EQ.IABS(M0).AND.IABS(M ).EQ.IABS(MP0))) GO TO 1
    2 FF=1.D0
      IF(IABS(M).LE.L.AND.IABS(MP).LE.L) GO TO 3
      WRITE(14,100) L,M,MP
  100 FORMAT('     L=',I5,'    M=',I5,'    MP =',I5)
      RETURN
    3 N1   =L+M
      N2   =L-M
      N3   =L+MP
      N4   =L-MP
      L0=L
      M0=M
      MP0=MP
      DO 4 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 4
      DO 5 I=1,NN
    5 FF=FF*I
    4 CONTINUE
      FF=DSQRT(FF)
    1 BETA2=BETA/2.D0
      COSB=DCOS(BETA2)
      SINB=-DSIN(BETA2)
      IF(DABS(COSB).LT.1.D-4) GO TO 9
      IF(DABS(SINB).LT.1.D-4) GO TO 11
      KMAX=MIN0(L-MP,L+M)
      KMIN=MAX0(M-MP,0)
      TERM=COSB**(2*L+M-MP-2*KMIN)*SINB**(MP-M+2*KMIN)*FF
      GO TO 12
    9 LTRM=L
      TERM=FF
      IF(SINB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
      GO TO 14
   11 LTRM=0
      TERM=FF
      IF(COSB.LT.0.D0.AND.MOD(MP-M,2).NE.0) TERM=-TERM
   14 KMAX=M-MP
      IF(MOD(KMAX,2).NE.0) GO TO 13
      KMAX=LTRM+KMAX/2
      IF(KMAX.LT.MAX0(M-MP,0)) GO TO 13
      IF(KMAX.GT.MIN0(L-MP,L+M)) GO TO 13
      KMIN=KMAX
   12 IF(MOD(KMIN,2).NE.0) TERM=-TERM
      N1   =L-MP-KMIN
      N2   =L+M-KMIN
      N3   =KMIN+MP-M
      N4   =KMIN
      DO 6 N=1,4
      NN=NF(N)
      IF(NN.EQ.0) GO TO 6
      DO 7 I=1,NN
    7 TERM=TERM/I
    6 CONTINUE
      DROTB=TERM
      IF(KMIN.EQ.KMAX) RETURN
      KMIN=KMIN+1
      COSB=COSB**2
      SINB=SINB**2
      N3=N3   +1
      DO 8 K=KMIN,KMAX
      TERM=-N1*N2*TERM*SINB/(COSB*K*N3)
      DROTB=DROTB+TERM
      N1=N1-1
      N2=N2-1
    8 N3=N3+1
      RETURN
   13 DROTB=0.D0
      RETURN
      END
