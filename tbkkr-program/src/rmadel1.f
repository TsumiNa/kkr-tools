c 30.06.97 ***************************************************************
      SUBROUTINE RMADEL1(AVMAD,BVMAD,ALAT,NATYP,LMPOT,
     +                  IMAD,JMAD,NVMAD)
      implicit none
c ************************************************************************
c
c  read coeffients for madelung potential
c
c modified for use of madelung 'shells'   p.zahn, june 97
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
      INTEGER LMPOT,NATYP,NVMAD
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVMAD(NVMADD,LMPOTD,*),
     +                 BVMAD(NVMADD,*)
      INTEGER
     +     IMAD(NAEZD,*),           ! madelung constant for 
                                    ! position pair(i,j)
     +     JMAD(NVMADD,*)           ! position indeces for madelung constant 
C     ..
C     .. Locals ..
      DOUBLE PRECISION BV,VALAT,RFPI
      DOUBLE PRECISION AV1(LMPOTD,LMPOTD),
     +                 BV1(LMPOTD)
      INTEGER I1,I2,II,IJ1,IJ2,IOS,LM,LM2,LM2POT,NV,I,IM,
     +     INMAD(NVMADD)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL RCSTOP,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
c ------------------------------------------------------------------------
      write(6,*) '>>> RMADEL : Read Madelung Potential Coefficients'

      RFPI=DSQRT(16.D0*DATAN(1.0D0))

      NV=0
      READ (40,FMT='(i5,d20.12,I6)',ERR=900) LM2POT,VALAT,NV
 900  CONTINUE
      write(6,*) 'NV :',nv

      IF (LM2POT.LT.LMPOT .OR. LM2POT.GT.LMPOTD .OR.
     +    ABS(ALAT-VALAT).GT.1.0D-6) THEN
        WRITE (6,FMT=*) '  lmpot : ',LMPOT,LM2POT,LMPOTD
        WRITE (6,FMT=*) '  alat  : ',ALAT,VALAT
        CALL RCSTOP('VMAD   1')
      END IF
      IF (NV.GT.0 .AND. NV.LT.NVMAD) THEN
        WRITE (6,FMT=*) '  NVMAD : ',NVMAD,NV
        CALL RCSTOP('VMAD   2')
      END IF

      IF (NV.EQ.0) THEN
c
c ---> old madelung constants file
c
        write(6,*) 'OLD FORMAT ------------',NVMAD
        DO 101 I=1,NVMAD
          INMAD(I) = 0
 101    END DO
        write(6,*) 'BVMAD'
        DO 20 IJ1 = 1,NATYP
          DO 10 IJ2 = 1,NATYP
            READ (40,FMT='(2i5)') I1,I2
            IF (I1.NE.IJ1 .OR. I2.NE.IJ2) STOP 'rmadel - BVMAD'
            READ (40,FMT='(4d20.12)') (BV1(LM),LM=1,LM2POT)
            IM = IMAD(I1,I2)
            IF (INMAD(IM).EQ.0) THEN
              INMAD(IM) = 1
              DO 90 LM=1,LM2POT
                BVMAD(IM,LM) = BV1(LM)
 90           END DO
              IF (TEST('IMAD    ')) 
     +             WRITE(6,FMT='(2I4,I8)') I1,I2,IM
            END IF
 10       CONTINUE
 20     CONTINUE
c
        DO 102 I=1,NVMAD
          INMAD(I) = 0
 102    END DO
        write(6,*) 'AVMAD'
        DO 40 IJ1 = 1,NATYP
          DO 30 IJ2 = 1,NATYP
            READ (40,FMT='(2i5)') I1,I2
            IF (I1.NE.IJ1 .OR. I2.NE.IJ2) STOP 'rmadel - AVMAD'
            READ (40,FMT='(4d20.12)') ((AV1(LM,LM2),LM=1,LM2POT),
     +           LM2=1,LM2POT)
            IM = IMAD(I1,I2)
            IF (INMAD(IM).EQ.0) THEN
              INMAD(IM) = 1
              DO 70 LM=1,LM2POT
                DO 80 LM2=1,LM2POT
                  AVMAD(IM,LM,LM2) = AV1(LM,LM2)
 80             END DO
 70           END DO
              IF (TEST('IMAD    ')) 
     +             WRITE(6,FMT='(2I4,I8)') I1,I2,IM
            END IF
 30       CONTINUE
          II = 1
          DO 105 I=1,NVMAD
            IF (INMAD(I).EQ.0) THEN
              II = 0
              end if
 105      END DO
          IF (II.EQ.1) THEN
            write(6,*) 'All AVMAD read. I1 = ',I1,'.'
            RETURN
          END IF
 40     CONTINUE

      ELSE                          ! (NV.EQ.0)

        write(6,*) 'NEW FORMAT ++++++++++++'
        DO 104 I = 1,NVMAD
          INMAD(I) = 0
 104    END DO

        write(6,*) 'BVMAD'
        DO 100 II = 1,NV
          READ (40,FMT='(2i5)') I1,I2
          READ (40,FMT='(4d20.12)') (BV1(LM),LM=1,LM2POT)
          IM = IMAD(I1,I2)
          IF (INMAD(IM).EQ.0) THEN
            INMAD(IM) = 1
            DO 110 LM=1,LM2POT
              BVMAD(IM,LM) = BV1(LM)
 110        END DO
            IF (TEST('IMAD    ')) 
     +           WRITE(6,FMT='(2I4,I8)') I1,I2,IM
          END IF
 100    CONTINUE                    ! II = 1,NV
c
        DO 103 I = 1,NVMAD
          IF (INMAD(I).EQ.0) THEN
            write(6,*) 'SHELL',I,'NOT FOUND IN FILE 40'
            CALL RCSTOP('VMAD   3')
          END IF
          INMAD(I) = 0
 103    END DO

        write(6,*) 'AVMAD'
        DO 120 II = 1,NV
          READ (40,FMT='(2i5)') I1,I2
          READ (40,FMT='(4d20.12)') ((AV1(LM,LM2),LM=1,LM2POT),
     +         LM2=1,LM2POT)
          IM = IMAD(I1,I2)
          IF (INMAD(IM).EQ.0) THEN
            INMAD(IM) = 1
            DO 130 LM=1,LM2POT
              DO 140 LM2=1,LM2POT
                AVMAD(IM,LM,LM2) = AV1(LM,LM2)
 140          END DO
 130        END DO
            IF (TEST('IMAD    ')) 
     +           WRITE(6,FMT='(2I4,I8)') I1,I2,IM
          END IF
 120    CONTINUE                    ! II = 1,NV
c
c ---> shift BVMAD(i,1) and AVMAD(i,1,1)
c
        BV = 0.D0
        DO 106 I = 1,NATYP
          BV = BV + BVMAD(I,1)/DBLE(NATYP)
 106    END DO

        DO 107 I = 1,NVMAD
          BVMAD(I,1)   = BVMAD(I,1) - BV
          AVMAD(I,1,1) = AVMAD(I,1,1) + BV*RFPI
 107    END DO


      END IF                        ! (NV.EQ.0)

      write(6,*) '<<< RMADEL'

      RETURN

      END
