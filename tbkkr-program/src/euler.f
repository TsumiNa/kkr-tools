c 20.09.96 ***************************************************************
      SUBROUTINE EULER(ND,ALPHA,BETA,GAMMA)
      implicit none
c ************************************************************************
C  CALCULATE THE EULER ANGLES FOR THE ROTATION MATRIX ND
c
c  modified from OHROT
c
c  by p. zahn, sept. 96
c ------------------------------------------------------------------------
c     .. arguments
      DOUBLE PRECISION ND(3,3),ALPHA,BETA,GAMMA
C     .. locals
      DOUBLE PRECISION D(3,3)
      DOUBLE PRECISION ALPHAK,BETAK,GAMMAK,DIFF,
     +                 SINA,SINB,SINC,
     +                 COSA,COSB,COSC,
     +                 PI
      INTEGER I,J
C
C     .. Scalars in Common ..
      LOGICAL strco1
      CHARACTER*5 LMACH
C
C     .. Common blocks ..
      COMMON /CMACH/strco1,LMACH
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL RCSTOP,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DABS,DACOS,DASIN,DATAN,DCOS,DSIN,SIN

      DOUBLE PRECISION SMALL
      DATA SMALL / 1.D-10 /
C-----------------------------------------------------------------------
      PI = 4.D0*DATAN(1.D0)
      BETA = DACOS(ND(3,3))
      IF (DABS(ND(3,3)-1.d0) .LT. SMALL ) THEN
        GAMMA = 0.0D0
        ALPHA = DACOS(ND(1,1))
        IF (DABS(DSIN(ALPHA)-ND(2,1)).GT.SMALL) 
     +       ALPHA = - ALPHA 
        
      ELSE
        IF (DABS(ND(3,3)+1.d0) .LT. SMALL) THEN
          GAMMA = 0.0D0
          ALPHA = -DACOS(-ND(1,1))
          IF (DABS(-DSIN(ALPHA)-ND(2,1)).GT.SMALL) 
     +         ALPHA = - ALPHA
          
        ELSE
          ALPHA = DASIN(ND(2,3)/SIN(BETA))
          IF (DABS(DCOS(ALPHA)*DSIN(BETA)-ND(1,3)).GT.SMALL) 
     +         ALPHA = PI - ALPHA
          GAMMA = DASIN(ND(3,2)/SIN(BETA))
          IF (DABS(-DCOS(GAMMA)*DSIN(BETA)-ND(3,1)).GT.SMALL) 
     +         GAMMA = PI - GAMMA
        END IF
        
      END IF
      
C     WRITE(6,291) K,ALPHA(K)/DATAN(1.D0),
C     1             K,BETA(K)/DATAN(1.D0),K,GAMMA(K)/DATAN(1.D0)
C     291  FORMAT(6X,'ALPHA(',I2,')=',F15.9,'*DATAN(1.D0)',/,
C     1       6X,'BETA(',I2,') =',F15.9,'*DATAN(1.D0)',/,
C     2       6X,'GAMMA(',I2,')=',F15.9,'*DATAN(1.D0)')
C-----------------------------------------------------------------------
C     BELOW THE ROTATION MATRICES ARE RECALCULATED FROM THE EULER ANGLES.
C-----------------------------------------------------------------------
      COSA = DCOS(ALPHA)
      COSB = DCOS(BETA)
      COSC = DCOS(GAMMA)
      SINA = DSIN(ALPHA)
      SINB = DSIN(BETA)
      SINC = DSIN(GAMMA)
      D(1,1) = COSA*COSB*COSC - SINA*SINC
      D(2,1) = SINA*COSB*COSC + COSA*SINC
      D(3,1) = -SINB*COSC
      D(1,2) = -COSA*COSB*SINC - SINA*COSC
      D(2,2) = -SINA*COSB*SINC + COSA*COSC
      D(3,2) = SINB*SINC
      D(1,3) = COSA*SINB
      D(2,3) = SINA*SINB
      D(3,3) = COSB
C-----------------------------------------------------------------------
C  THE EQUIVALENCE OF THE RECALCULATED AND THE ORIGINAL ROTATION
C  MATRICES IS CHECKED BELOW
C-----------------------------------------------------------------------
      DIFF = 0.0D0
      DO 230 I = 1,3
        DO 220 J = 1,3
          DIFF = DIFF + DABS(ND(I,J)-D(I,J))
 220    CONTINUE
 230  CONTINUE
      IF (DABS(DIFF).GT.1.D-4) THEN
        WRITE (6,FMT=9000)
        WRITE (6,FMT=9010) ALPHA/PI,BETA/PI,GAMMA/PI
        WRITE (6,FMT=9020) ((ND(I,J),J=1,3), (D(I,J),J=1,3),I=1,3)
        CALL RCSTOP('EULER   ')
        
      ELSE
        IF (ALPHA.LT.0.0D0) ALPHA = ALPHA + 8.D0*DATAN(1.D0)
        IF (GAMMA.LT.0.0D0) GAMMA = GAMMA + 8.D0*DATAN(1.D0)
        ALPHAK = ALPHA/4.D0/DATAN(1.D0)
        BETAK = BETA/4.D0/DATAN(1.D0)
        GAMMAK = GAMMA/4.D0/DATAN(1.D0)

        IF (TEST('ND      ')) THEN
          WRITE(6,FMT='(3F12.4)') ALPHAK,BETAK,GAMMAK
          WRITE(6,*) 
          WRITE(6,FMT='(3F12.4)') ((ND(I,J),J=1,3),I=1,3)
        END IF
      END IF                     

      RETURN

 9000 FORMAT (' ERROR IN ROTATION MATRICES')
 9010 FORMAT (4X,'ALPHA(K)=',F4.2,'*PI',4X,' BETA(K)=',F4.2,'*PI',4X,
     +       'GAMMA(K)=',F4.2,'*PI')
 9020 FORMAT (3 (4X,3F6.2,4X,3F6.2,/))

      END
