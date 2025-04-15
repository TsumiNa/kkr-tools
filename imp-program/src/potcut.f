      SUBROUTINE POTCUT(IMT1,IRC1,INS,LMPOT,R,VM2Z,vspsme,VINS,Z1)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     set potential equal zero between muffin tin sphere and
c       outer sphere
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER IRMD,IRNSD
      PARAMETER (irmd=1484,irnsd=508)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      REAL*8 Z1
      INTEGER IMT1,INS,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      REAL*8 R(*),VINS(IRMIND:IRMD,*),VM2Z(*),vspsme(*)
C     ..
C     .. Local Scalars ..
      INTEGER IR,IST,LM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
      DO 10 IR = IMT1 + 1,IRC1
        VM2Z(IR) = 2.0D0*Z1/R(IR)
        vspsme(ir) = 2.0d0*z1/r(ir)
   10 CONTINUE
c
      IF (INS.GE.1) THEN
        IST = MAX(IRMIND,IMT1+1)
        DO 30 IR = IST,IRC1
          DO 20 LM = 2,LMPOT
            VINS(IR,LM) = 0.0D0
   20     CONTINUE
   30   CONTINUE
      END IF

      END
