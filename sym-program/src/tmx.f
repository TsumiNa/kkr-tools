      SUBROUTINE TMX (IREP,
     1          ALPHA,BETA,GAMMA,IDEN
     2          ,KTAU,INVER)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (N48=48)
      DIMENSION ALPHA(N48),BETA(N48),GAMMA(N48),IDEN(N48)
      DIMENSION RA(5,5),RB(5,5),RG(5,5),T(5,5,N48)
      IF (IREP .EQ. 0) RETURN
      DO 10 KT=1,KTAU
      S=1.0D0
      IF (KT .GT. INVER) S=-1.0D0
      CA=DCOS(ALPHA(KT))
      SA=DSIN(ALPHA(KT))
      CB=DCOS(BETA(KT))
      SB=DSIN(BETA(KT))
      CG=DCOS(GAMMA(KT))
      SG=DSIN(GAMMA(KT))
      T(1,1,KT)=(CA*CB*CG-SA*SG)*S
      T(1,2,KT)=(SA*CB*CG+CA*SG)*S
      T(1,3,KT)=-SB*CG*S
      T(2,1,KT)=(-CA*CB*SG-SA*CG)*S
      T(2,2,KT)=(-SA*CB*SG+CA*CG)*S
      T(2,3,KT)=SB*SG*S
      T(3,1,KT)=CA*SB*S
      T(3,2,KT)=SA*SB*S
10    T(3,3,KT)=CB*S
      DO 30 KT=1,KTAU,6
      KTM=KT-1
      WRITE(14,100) (IDEN(KTM+K),K=1,6)
      DO 20 J=1,3
20    WRITE (14,110) ((T(I,J,KTM+K),I=1,3),K=1,6)
30    CONTINUE
      IF (IREP .EQ. 1) RETURN
      IF (IREP .EQ. -1) GO TO 90
      RT=DSQRT(3.0D0)
      DO 60 KT=1,KTAU
      DO 40 I=1,5
      DO 40 J=1,5
      RA(I,J)=0.0D0
      RB(I,J)=0.0D0
      RG(I,J)=0.0D0
40    T(I,J,KT)=0.0D0
      CA=DCOS(ALPHA(KT))
      SA=DSIN(ALPHA(KT))
      CB=DCOS(BETA(KT))
      SB=DSIN(BETA(KT))
      CG=DCOS(GAMMA(KT))
      SG=DSIN(GAMMA(KT))
      RA(1,1)=1.0D0
      RA(2,2)=2.0D0*CA*CA-1.0D0
      RA(2,3)=2.0D0*CA*SA
      RA(3,2)=-RA(2,3)
      RA(3,3)=RA(2,2)
      RA(4,4)=CA
      RA(4,5)=SA
      RA(5,4)=-SA
      RA(5,5)=CA
      RG(1,1)=1.0D0
      RG(2,2)=2.0D0*CG*CG-1.0D0
      RG(2,3)=2.0D0*CG*SG
      RG(3,2)=-RG(2,3)
      RG(3,3)=RG(2,2)
      RG(4,4)=CG
      RG(4,5)=SG
      RG(5,4)=-SG
      RG(5,5)=CG
      RB(2,2)=0.5D0*CB*CB+0.5D0
      RB(1,1)=3.0D0*RB(2,2)-2.0D0
      RB(4,4)=4.0D0*RB(2,2)-3.0D0
      RB(3,3)=CB
      RB(3,5)=-SB
      RB(5,3)=SB
      RB(4,2)=CB*SB
      RB(2,4)=-RB(4,2)
      RB(1,4)=RB(4,2)*RT
      RB(4,1)=-RB(1,4)
      RB(5,5)=CB
      RB(1,2)=0.5D0*RT*SB*SB
      RB(2,1)=RB(1,2)
      DO 50 I=1,5
      DO 50 J=1,5
      DO 50 K=1,5
      DO 50 M=1,5
50    T(I,J,KT)=T(I,J,KT)+RG(I,K)*RB(K,M)*RA(M,J)
60    CONTINUE
      DO 80 KT=1,KTAU,3
      KTM=KT-1
      WRITE (14,120) (IDEN(KTM+K),K=1,3)
      DO 70 J=1,5
70    WRITE (14,130) ((T(I,J,KTM+K),I=1,5),K=1,3)
80    CONTINUE
      IF (IREP .EQ. 2) RETURN
90    WRITE (14,140) IREP
      STOP
100   FORMAT (//2X,5(7X,A4,7X),7X,A4)
110   FORMAT (2H  ,6(1H*,F6.3,2F6.3))
120   FORMAT (//2X,2(13X,A4,13X),13X,A4)
130   FORMAT (2H  ,3(1H*,F6.3,4F6.3 ))
140   FORMAT (//50X,'  STOP  IN  TMX   IREP= ',I6)
      END
