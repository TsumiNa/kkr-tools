      SUBROUTINE JMTRX(X1,X2,X3,E,LMAX,JMAT,LCALL)
      implicit none
C ****************************************************************
C * This subroutine calculates the transformation matrix Ull'(E) *
C * for the Green's function in the shifted position. According  *
C * N. Stefanou et al PRB 36 6372 (1987) eq.(6).                 *
C * Lmax is the maximum l cutoff for the matrix Ull':=JMAT       *
C * Since the Gaunt coefficients are used, subroutines gaunt and *
C * Gaunt2 have to be called first to set up the common block    *
C * GAUNTC.                                                      *
C ****************************************************************
C
C     .. Parameters ..
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
      INTEGER LMX,LPOTD,LMAXD
      PARAMETER (lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER LMMAXD,LM2D,LMPOTD
      PARAMETER (LMMAXD= (LMAXD+1)**2,LM2D= (2*LMAXD+1)**2,
     +          LMPOTD= (LPOTD+1)**2)
      INTEGER TWOLMAX,LMMAX2D
      PARAMETER (TWOLMAX=2*LMAXD,LMMAX2D=(TWOLMAX+1)**2)
      INTEGER N,LASSLD
      PARAMETER (N=4*LMAXD,LASSLD=N)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
C
      INTEGER LMAX,IA,L1,L2,L3,LM1,LM2,LM3
      INTEGER LMAXSQ,I,LM,L
      COMPLEX*16 EI,ARG,ZSUM,CZERO,E,CLLL1
      REAL*8 SABS,FOURP,ZERO,X1,X2,X3,C0LL,X123
      LOGICAL  LCALL
C
C     .. Array arguments ..
C
      REAL*8 YLM(LMMAX2D)
      COMPLEX*16 BJ(0:TWOLMAX),H(0:TWOLMAX),Y(0:TWOLMAX)
      COMPLEX*16 JMAT(LMMAXD,LMMAXD)
C
C     .. Intrinsic functions ..
C
      INTRINSIC DSQRT,DATAN,CDSQRT
C
C     .. External subroutines ..
C
      EXTERNAL YMY,BESSEL1
C
C     .. Scalars in Common ..
C
      INTEGER IEND
C     ..
C     .. Arrays in Common ..
C
      REAL*8 CLEB(NCLEB,2)
      INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(LM2D)
C
      COMMON /GAUNTC/CLEB,LOFLM,ICLEB,IEND,JEND
      DATA EI/(0.0D0,1.0D0)/
      DATA CZERO/(0.0D0,0.0D0)/
      FOURP=16.D0*DATAN(1.D0)
      C0LL=1.D0/SQRT(FOURP)
      LMAXSQ=(LMAX+1)**2
C
      X123=SQRT(X1*X1+X2*X2+X3*X3)
      IF (X123.LT.1.D-10) GOTO 99
      CALL YMY(X1,X2,X3,SABS,YLM,2*LMAX)
      ARG=SQRT(E)*SABS
      CALL BESSEL1(BJ,Y,H,ARG,TWOLMAX,2*LMAX,.TRUE.,.TRUE.,.TRUE.
     &           ,LCALL)
      DO 2 LM=1,LMAXSQ
      DO 2 LM1=1,LMAXSQ
  2   JMAT(LM,LM1)=CZERO
C
      DO 1 IA=1,IEND
      LM1=ICLEB(IA,1)
      LM2=ICLEB(IA,2)
      LM3=ICLEB(IA,3)
      L1=LOFLM(LM1)
      L2=LOFLM(LM2)
      L3=LOFLM(LM3)
C
c calculate jmat for lm1.ge.lm2
C
      CLLL1=FOURP*CLEB(IA,1)*(EI)**(L1+L3-L2)
      JMAT(LM1,LM2)=JMAT(LM1,LM2)+CLLL1*BJ(L3)*YLM(LM3)
   1  CONTINUE
C
c add the l=0 component in the diagonal elements.
C
      DO 6 LM1=1,LMAXSQ
  6   JMAT(LM1,LM1)=JMAT(LM1,LM1)+FOURP*C0LL*BJ(0)*YLM(1)
c
C   Create the rest of the matrix by symmetry relation.
c
      DO 5 LM1=1,LMAXSQ
      DO 5 LM2=1,LM1-1
      L1=LOFLM(LM1)
      L2=LOFLM(LM2)
 5    JMAT(LM2,LM1)=JMAT(LM1,LM2)*(-1.D0)**(L1+L2)
      RETURN
 99   CONTINUE
      DO 4 LM=1,LMAXSQ
      DO 3 LM1=1,LMAXSQ
  3   JMAT(LM,LM1)=CZERO
  4   JMAT(LM,LM)=CMPLX(1.0D0,0.0D0)
      RETURN
      END
