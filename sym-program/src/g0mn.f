      SUBROUTINE G0MN(GMN,GM,NIMP,LMAXSQ,M,DROT,IROTMN,ISHELL)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER LMAXD,NIMPD
      PARAMETER (LMAXD=4,NIMPD=4290)
      INTEGER LMMAXD,NGD
      PARAMETER (LMMAXD= (LMAXD+1)**2,NGD=LMMAXD*NIMPD)
      COMPLEX*16 CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,M,NIMP
C     ..
C     .. Array Arguments ..
      COMPLEX*16 GM(LMMAXD,LMMAXD,*),GMN(NGD,LMMAXD)
      REAL*8 DROT(48,LMMAXD,LMMAXD)
      INTEGER IROTMN(NIMPD,0:NIMPD),ISHELL(NIMPD,0:NIMPD)
C     ..
C     .. Local Scalars ..
      INTEGER I,IROT,IS,J,N,NLM1
C     ..
C     .. Local Arrays ..
      COMPLEX*16 DGD(LMMAXD,LMMAXD),DLL(LMMAXD,LMMAXD),
     +               GD(LMMAXD,LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL ZGEMM
C     ..
      DO 50 N = 1,NIMP
        IS = ISHELL(M,N)
        IROT = IROTMN(M,N)
C----  PRODUCE DLL AND ROTATE GM INTO GMN
        DO 20 J = 1,LMAXSQ
          DO 10 I = 1,LMAXSQ
            DLL(I,J) = DROT(IROT,I,J)
   10     CONTINUE
   20   CONTINUE
        CALL ZGEMM('T','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,DLL,LMMAXD,
     +             GM(1,1,IS),LMMAXD,CZERO,GD,LMMAXD)
        CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,GD,LMMAXD,DLL,
     +             LMMAXD,CZERO,DGD,LMMAXD)
        DO 40 J = 1,LMAXSQ
          NLM1 = J + LMMAXD* (N-1)
          DO 30 I = 1,LMAXSQ
            GMN(NLM1,I) = DGD(J,I)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
      RETURN
      END
