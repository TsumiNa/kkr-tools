c ************************************************************************
      SUBROUTINE EREST(IP,EA,EE,II,NUM1,NUM2,W1,W2,IVAL1,IVAL2)
      implicit none
c ************************************************************************
c     p.zahn, april 96
c ------------------------------------------------------------------------
      include 'inc.fi'
      INTEGER II,IP,IVAL1,IVAL2,NUM1,NUM2
      DOUBLE PRECISION EA,EE,W1,W2

c     .. locals
      LOGICAL LTEST
      INTEGER I,J,MAX

c     .. arrays in common
      DOUBLE PRECISION EE1(0:M2D,KPOIBZ),EE2(0:M2D,KPOIBZ),
     +                 WW1(0:M2D,KPOIBZ),WW2(0:M2D,KPOIBZ)
      INTEGER II1(0:M2D,KPOIBZ),II2(0:M2D,KPOIBZ)
      COMMON /EE/ EE1,EE2,WW1,WW2,II1,II2 ! index arrays for eigenvalue
                                          ! analysis

      INTEGER NBASIS(KPOIBZ),NMIN(KPOIBZ),NMAX(KPOIBZ)
      COMMON /NBASIS / NBASIS,NMIN,NMAX

      DOUBLE PRECISION ZERO,TH
      DATA ZERO / 0.0D0 /
     +     TH   / 5.0D2 /

      SAVE
c ------------------------------------------------------------------------
      I = II - NMIN(IP)
      MAX = NMAX(IP) - NMIN(IP)

      IF (MAX.GT.M2D) THEN
        write(6,*) 'IP,NMIN,NMAX :',IP,NMIN(IP),NMAX(IP)
        STOP 'EREST'
      END IF

c      write(6,*) 'EREST, ii,i ',II,I
c      write(6,*) 'EREST, nmin,nmax ',NMIN(IP),NMAX(IP)

      LTEST = .TRUE.
      DO 10 J = I-1,0,-1
        IF (EE1(J,IP).LT.TH .AND. LTEST) THEN
          EA    = EE1(J,IP)
          W1    = WW1(J,IP)
          IVAL1 = II1(J,IP)
          NUM1  = J + NMIN(IP)
          LTEST = .FALSE.
        END IF
 10   END DO

      LTEST = .TRUE.
      DO 20 J = I,MAX
        IF (EE2(J,IP).LT.TH .AND. LTEST) THEN
          EE    = -EE2(J,IP)
          W2    = WW2(J,IP)
          IVAL2 = II2(J,IP)
          NUM2  = J + NMIN(IP)
          LTEST = .FALSE.
        END IF
 20   END DO

      RETURN
      END 
