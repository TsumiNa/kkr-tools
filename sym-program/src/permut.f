      SUBROUTINE PERMUT(R,KT,
     1          ALPHA,BETA,GAMMA,
     2          MT,MTAU,INVER)
C-----------------------------------------------------------------------
C THE SUBROUTINE CALCULATES THE PERMUTATION OF THE SITES IN THE
C CONSIDERED SHELL FOR THE SYMMETRY OPERATION DEFINED BY KT.
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N48=48)
      DIMENSION R(3,N48)
      DIMENSION ALPHA(N48),BETA(N48),GAMMA(N48)
      DIMENSION MT(N48,N48)
      DIMENSION A(3,3),S(3)
      CA=DCOS(ALPHA(KT))
      SA=DSIN(ALPHA(KT))
      CB=DCOS(BETA(KT))
      SB=DSIN(BETA(KT))
      CG=DCOS(GAMMA(KT))
      SG=DSIN(GAMMA(KT))
      A(1,1)=CA*CB*CG-SA*SG
      A(2,1)=SA*CB*CG+CA*SG
      A(3,1)=-SB*CG
      A(1,2)=-CA*CB*SG-SA*CG
      A(2,2)=-SA*CB*SG+CA*CG
      A(3,2)=SB*SG
      A(1,3)=CA*SB
      A(2,3)=SA*SB
      A(3,3)=CB
      DO 40 M=1,MTAU
      DO 10 I=1,3
      S(I)=0.0D0
      DO 10 J=1,3
10    S(I)=S(I)+A(I,J)*R(J,M)
      IF (KT .LE. INVER) GO TO 20
      S(1)=-S(1)
      S(2)=-S(2)
      S(3)=-S(3)
20    DO 30 MP=1,MTAU
      RS=DABS(R(1,MP)-S(1))+DABS(R(2,MP)-S(2))+DABS(R(3,MP)-S(3))
      IF (RS .LT. 0.01D0) GO TO 40
30    CONTINUE
      WRITE (14,50) KT,M,S(1),S(2),S(3)
      STOP
40    MT(KT,M)=MP
      RETURN
50    FORMAT (20X,'STOP IN PERMUTE',3X,2I5,3D14.7)
      END
