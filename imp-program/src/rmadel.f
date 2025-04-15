      SUBROUTINE RMADEL(AVMAD,BVMAD,ALAT,NATYP,LMPOT)
      IMPLICIT NONE
c
c---> read coeffients for madelung potential
c
C     .. Parameters ..
      INTEGER NATYPD,NTREFD
      PARAMETER (NATYPD=38,NTREFD=1)
      INTEGER LPOTD
      PARAMETER (lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALAT
      INTEGER LMPOT,NATYP
C     ..
C     .. Array Arguments ..
      REAL*8 AVMAD(NTREFD,NTREFD,LMPOTD,LMPOTD),
     +                 BVMAD(NTREFD,NTREFD,LMPOTD)
C     ..
C     .. Local Scalars ..
      REAL*8 VALAT
      INTEGER I1,I2,IJ1,IJ2,LM,LM2,LM2POT
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      READ (40,FMT='(i5,d20.12)') LM2POT,VALAT

      IF (LM2POT.LT.LMPOT .OR. LM2POT.GT.LMPOTD .OR.
     +    ABS(ALAT-VALAT).GT.1.0D-6) THEN
        WRITE (6,FMT=*) '  lmpot : ',LMPOT,LM2POT,LMPOTD
        WRITE (6,FMT=*) '  alat  : ',ALAT,VALAT
        CALL RCSTOP('vmad    ')

      END IF

      DO 20 IJ1 = 1,NATYP
        DO 10 IJ2 = 1,NATYP
          READ (40,FMT='(2i5)') I1,I2
          READ (40,FMT='(4d20.12)') (BVMAD(I1,I2,LM),LM=1,LM2POT)
   10   CONTINUE
   20 CONTINUE
c
      DO 40 IJ1 = 1,NATYP
        DO 30 IJ2 = 1,NATYP
          READ (40,FMT='(2i5)') I1,I2
          READ (40,FMT='(4d20.12)') ((AVMAD(I1,I2,LM,LM2),LM=1,LM2POT),
     +      LM2=1,LM2POT)
   30   CONTINUE
   40 CONTINUE
      RETURN

      END
