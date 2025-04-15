c ************************************************************************
      SUBROUTINE THCORE(E,R,T,F)
      implicit none
c ************************************************************************
      include 'inc.fi'
      include 'inc.cls'
      INTEGER LMMAXD,LMAX
      PARAMETER (LMMAXD= (LMAXD+1)**2,LMAX=LMAXD)

      DOUBLE COMPLEX E,T(LMMAXD,LMMAXD)
      DOUBLE PRECISION F,R

c
      DOUBLE PRECISION 
     +     A,B,DELTA,PI,RWS,RMT
      DOUBLE COMPLEX 
     +     HL(0:LMAX),
     +     JL(0:LMAX),
     +     NL(0:LMAX),
     +     Z,S,TI

      INTEGER I,LM,L,M

      DOUBLE COMPLEX EXPIDL
      EXTERNAL EXPIDL

      INTRINSIC DREAL,DSQRT,SQRT

      DOUBLE COMPLEX CI
      PARAMETER (CI = (0.D0,1.D0))
      PARAMETER (PI = 3.14159265358979312D0)
      INTEGER LF(144)
      DATA LF /0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11 /
c ------------------------------------------------------------------------
      CALL CINIT(LMMAXD*LMMAXD,T)

      Z = F*R*SQRT(E)   
      CALL BESHAN(HL,JL,NL,Z,0,LMAX)

      LM = 1
      DO 91 L = 0,LMAX
        TI = SQRT(E)*(-NL(L)/JL(L)+CI)
        DO 92 M = -L,L
          T(LM,LM) = 1.d0/TI
          LM = LM + 1
 92     END DO
 91   END DO

      END
