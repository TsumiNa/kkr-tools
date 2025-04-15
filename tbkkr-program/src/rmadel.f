c 17.10.95 ***************************************************************
      SUBROUTINE RMADEL(AVMAD,BVMAD,ALAT,NATYP,LMPOT)
      implicit none
c ************************************************************************
c
c  read coeffients for madelung potential
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD,NTREFD
c      PARAMETER (NATYPD=1,NTREFD=0)
c      INTEGER LPOTD
c      PARAMETER (LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER LMPOT,NATYP
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVMAD(NVMADD,NVMADD,LMPOTD,*),
     +                 BVMAD(NVMADD,NVMADD,*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION VALAT
      INTEGER I1,I2,IJ1,IJ2,LM,LM2,LM2POT
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      write(6,*) 'RMADEL : Read Madelung Potential Coefficients'

      READ (40,FMT='(i5,d20.12)') LM2POT,VALAT

      IF (LM2POT.LT.LMPOT .OR. LM2POT.GT.LMPOTD .OR.
     +    ABS(ALAT-VALAT).GT.1.0D-6) THEN
        WRITE (6,FMT=*) '  lmpot : ',LMPOT,LM2POT,LMPOTD
        WRITE (6,FMT=*) '  alat  : ',ALAT,VALAT
        CALL RCSTOP('VMAD    ')
      END IF

      DO 20 IJ1 = 1,NATYP
        DO 10 IJ2 = 1,NATYP
          READ (40,FMT='(2i5)') I1,I2
          IF (I1.NE.IJ1 .OR. I2.NE.IJ2) STOP 'rmadel - BVMAD'
          READ (40,FMT='(4d20.12)') (BVMAD(I1,I2,LM),LM=1,LM2POT)
   10   CONTINUE
   20 CONTINUE
c
      DO 40 IJ1 = 1,NATYP
        DO 30 IJ2 = 1,NATYP
          READ (40,FMT='(2i5)') I1,I2
          IF (I1.NE.IJ1 .OR. I2.NE.IJ2) STOP 'rmadel - AVMAD'
          READ (40,FMT='(4d20.12)') ((AVMAD(I1,I2,LM,LM2),LM=1,LM2POT),
     +      LM2=1,LM2POT)
   30   CONTINUE
   40 CONTINUE
      rewind(40)
       write(6,*) 'File madel was closed by rmadel'
      RETURN

      END
