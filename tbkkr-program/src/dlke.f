c ************************************************************************
      SUBROUTINE DLKE(GLLKE,ADUM,BZKP,IE)
      implicit none
c ************************************************************************
C     .. Parameters ..
      INTEGER NATOMD
      PARAMETER (NATOMD=79)
      INTEGER LMAX
      PARAMETER (LMAX=4)
      INTEGER LMAXSQ,NLM
      PARAMETER (LMAXSQ= (LMAX+1)**2,NLM=LMAXSQ*NATOMD)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ADUM
      INTEGER IE
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GLLKE(LMAXSQ,LMAXSQ)
      DOUBLE PRECISION BZKP(3)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,CONVPU
      INTEGER I,IELAST,IESAVE,IG,LM2,M,ML,N,N1,NATOM,NL
      LOGICAL LSTART
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX ARGX(NATOMD),ARGY(NATOMD),ARGZ(NATOMD),
     +               EIKR(NATOMD),GINP(NLM,LMAXSQ),GMN(NLM,LMAXSQ)
      DOUBLE PRECISION RM(3,NATOMD)
      INTEGER LF(144)
C     ..
C     .. External Subroutines ..
      EXTERNAL ZAXPY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC EXP
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
c      DATA LF/0,3*1,5*2,7*3,9*4/
      DATA LSTART/.TRUE./
C     ..
      IF (LSTART) THEN
        IESAVE = 0
        OPEN (69,FILE='gdos.input',status='old',FORM='unformatted')
        READ (69) ALAT,NATOM,IELAST
C       NATOM = 1
        WRITE (6,FMT=*) ALAT,NATOM,IELAST
        READ (69) ((RM(I,M),I=1,3),M=1,NATOM)
        WRITE (6,FMT=9000) ((RM(I,M),I=1,3),M=1,NATOM)
 9000   format(3f12.4)
        CONVPU = ALAT/2.D0/ 3.14159265358979312d0
        DO 10 M = 1,NATOM
          ARGX(M) = -CI/CONVPU*RM(1,M)
          ARGY(M) = -CI/CONVPU*RM(2,M)
          ARGZ(M) = -CI/CONVPU*RM(3,M)
   10   CONTINUE
        LSTART = .FALSE.
      END IF
      IF (IE.NE.IESAVE) THEN
        READ (69) ((GINP(N,M),M=1,LMAXSQ),N=1,LMAXSQ*NATOM)
c        WRITE(6,9010) (GINP(2+N*LMAXSQ-LMAXSQ,2),N=1,NATOM)
 9010 format(2f12.6)
        DO 40 M = 1,LMAXSQ
          DO 30 N = 1,LMAXSQ
            DO 20 N1 = 1,NATOM
              NL = N1*LMAXSQ - LMAXSQ + N
              ML = N1*LMAXSQ - LMAXSQ + M
              GMN(NL,M) = GINP(NL,M)
C             GMN(NL,M) = ((-1)** (LF(M)+LF(N))*GINP(ML,N)+
C    +                      GINP(NL,M))*0.5D0
C             GMN(NL,M) = ((-1)** (LF(M)+LF(N))*GINP(ML,N)+
C    +                      GINP(NL,M))*0.5D0
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
      END IF
      IESAVE = IE
      DO 60 M = 1,LMAXSQ
        DO 50 N = 1,LMAXSQ
          GLLKE(N,M) = CZERO
   50   CONTINUE
   60 CONTINUE

      DO 70 M = 1,NATOM
        EIKR(M) = EXP(BZKP(1)*ARGX(M)+BZKP(2)*ARGY(M)+BZKP(3)*ARGZ(M))
     +            *CONVPU
C    +            *CONVPU*6.76D0/2.D0/3.14159265358979312D0
   70 CONTINUE
      DO 90 M = 1,NATOM
        IG = 1 + (M-1)*LMAXSQ
        DO 80 LM2 = 1,LMAXSQ
          CALL ZAXPY(LMAXSQ,EIKR(M),GMN(IG,LM2),1,GLLKE(1,LM2),1)
   80   CONTINUE
   90 CONTINUE

      RETURN

      END
