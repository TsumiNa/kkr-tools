      SUBROUTINE G1MN(GMN,GM,GSTRC,LMMAXD,LMAXSQ,NIMP,IE,NAF,NBF,KGMN,
     +                ISZ,GR,NSHELL,NTYPIM)
      IMPLICIT NONE
C----------------------------------------------------------------------
C     THIS SUBROUTINE PRODUCES A SYMBOLIC REPRESENTATION FOR THE
C     STRUCTURAL GREEN FUNCTION MATRIX. THE SYMBOLIC REPRESENTATION
C     FROM 'G0STR'IS USED.
C----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IE,KGMN,LMAXSQ,LMMAXD,NIMP,NSHELL
C     ..
C     .. Array Arguments ..
      REAL*8 GM(LMMAXD,LMMAXD,*),GMN(LMMAXD,LMMAXD,*),GR(*)
      INTEGER GSTRC(LMMAXD,LMMAXD,*),NAF(*),NBF(*),NTYPIM(*)
      LOGICAL ISZ(*)
C     ..
C     .. Local Scalars ..
      REAL*8 SUM
      INTEGER I,IGMN3,IJABS,IJIS,IS,ISM,IV,J,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,MAX0
C     ..
      GR(IE) = 1.D0
C-----------------------------------------------------------------------
C     CONVERTING INTEGER 'GSTRC' TO REAL 'GM'
C-----------------------------------------------------------------------
      DO 50 IS = 1,NSHELL
        DO 20 J = 1,LMAXSQ
          DO 10 I = 1,LMAXSQ
            IJIS = GSTRC(I,J,IS)
            IJABS = MAX0(IJIS,-IJIS,1)
            GM(I,J,IS) = (IJIS/IJABS)*GR(IJABS)
   10     CONTINUE
   20   CONTINUE
        ISZ(IS) = .FALSE.
        SUM = 0.0D0
        DO 40 J = 1,LMAXSQ
          DO 30 I = 1,LMAXSQ
            SUM = SUM + DABS(GM(I,J,IS))
   30     CONTINUE
   40   CONTINUE
        IF (SUM.GT.1.D-6) ISZ(IS) = .TRUE.
   50 CONTINUE
      GR(IE) = 0.D0

      IV = 0
      IGMN3 = 0
      DO 80 N = 1,NIMP
        IV = IV + 1
        ISM = NTYPIM(IV)
        IF (.NOT.ISZ(ISM)) GO TO 80
        IGMN3 = IGMN3 + 1
        NAF(IGMN3) = N
        NBF(IGMN3) = N
        IF (IGMN3.EQ.1) THEN
          DO 70 I = 1,LMAXSQ
            DO 60 J = 1,LMAXSQ
              GMN(J,I,IGMN3) = GM(J,I,ISM)
   60       CONTINUE
   70     CONTINUE

        ELSE
          GO TO 90

        END IF

   80 CONTINUE
   90 CONTINUE
      KGMN = IGMN3
      RETURN

      END
