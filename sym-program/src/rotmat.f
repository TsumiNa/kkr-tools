      SUBROUTINE ROTMAT(LMAX,ND,DROT,YR,WTYR,RIJ,IJEND)
      IMPLICIT NONE
C     .. Parameters ..

      INTEGER LMAXD
      PARAMETER (LMAXD=4)
      INTEGER LMMAXD,IJD
      PARAMETER (LMMAXD= (LMAXD+1)**2,IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,LMAX
C     ..
C     .. Array Arguments ..
      REAL*8 DROT(48,LMMAXD,LMMAXD),RIJ(IJD,3),
     +                 WTYR(IJD,LMMAXD),YR(IJD,LMMAXD)
      INTEGER ND(48,3,3)
C     ..
C     .. Local Scalars ..
      INTEGER K,LM1,LM2,LMMAX
C     ..
C     .. Local Arrays ..
      REAL*8 C(LMMAXD,LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL ROTTNS
C     ..
C     .. Save statement ..
      SAVE
C     ..

      LMMAX = (LMAX+1)* (LMAX+1)
      DO K = 1,48
        CALL ROTTNS(C,K,LMAX,ND,YR,WTYR,RIJ,IJEND)
        DO LM2 = 1,LMMAXD
          DO LM1 = 1,LMMAXD
            DROT(K,LM1,LM2) = C(LM2,LM1)
          END DO
        END DO
      END DO
c     DO LM2 = 1,LMMAX
c       WRITE (6,FMT=9000) LM2
c       WRITE (6,FMT=9000) (LM1,DROT(1,LM1,LM2),LM1=1,LMMAX)
c     END DO
      RETURN

 9000 FORMAT (7 (I3,F8.4))
      END
