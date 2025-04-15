      SUBROUTINE SYTN0(AR,ADET,DETALF,INS,ITMAT,IREF,KTE,LMAX,NATPER,
     +                 NATREF,NSTART,NEND,TMATLL,NSHELL)
      IMPLICIT NONE
c
c     the delta t - matrix is calculated
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      COMPLEX*16 CONE
      PARAMETER (CONE= (1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      INTEGER INS,ITMAT,KTE,LMAX,NATPER,NATREF,NEND,NSTART
C     ..
C     .. Array Arguments ..
      COMPLEX*16 ADET(LMMAXD,LMMAXD),AR(LMMAXD,LMMAXD,*),
     +               DETALF(NATYPD),TMATLL(LMMAXD,LMMAXD,*)
      INTEGER IREF(*),NSHELL(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 DET
      INTEGER I,I1,IA,IHOST,INFO,LM1,LM2,LMMAX
C     ..
C     .. Local Arrays ..
      INTEGER IPVT(LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,TMREAD,ZCOPY,ZGETRF
C     ..
C     .. Save statement ..
      SAVE
C     ..

      LMMAX = (LMAX+1)* (LMAX+1)
c
c
c---> in the first time : read symmetrization coeffients
c


      DO 10 I1 = 1,NATREF
        CALL TMREAD(TMATLL(1,1,I1),ITMAT)
   10 CONTINUE

c
c---> loop over representive atoms
c
      DO 50 I1 = NSTART,NEND

        IA = I1 - NATREF
c
        IHOST = IREF(I1-NATREF)
        DO 30 LM2 = 1,LMMAX
          DO 20 LM1 = 1,LMMAXD
            TMATLL(LM1,LM2,I1) = TMATLL(LM1,LM2,I1) -
     +                           TMATLL(LM1,LM2,IHOST)
   20     CONTINUE
   30   CONTINUE
        DET = CONE
        IF (KTE.EQ.1 .AND. INS.GT.0) THEN
          CALL ZCOPY(LMMAXD**2,AR(1,1,I1),1,ADET,1)
          CALL ZGETRF(LMMAXD,LMMAXD,ADET,LMMAXD,IPVT,INFO)
          DO 40 I = 1,LMMAXD
            IF (IPVT(I).NE.I) DET = -DET
            DET = ADET(I,I)*DET
   40     CONTINUE
        END IF
        DETALF(I1) = DET
   50 CONTINUE
c
c---> set t matrices of cutted shells equal zero
c
      DO 60 I1 = NEND + 1,NATREF + NATPER
        CALL CINIT(LMMAXD**2,TMATLL(1,1,I1))
   60 CONTINUE

c
      RETURN

      END
