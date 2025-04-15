c ************************************************************************
      SUBROUTINE EREST1(IP,E1,IO)
      implicit none
c ************************************************************************
c     p.zahn, april 96
c ------------------------------------------------------------------------
      include 'inc.fi'
      INTEGER IP,IO
      DOUBLE PRECISION E1

c     .. locals
      INTEGER I,MAX,N
      DOUBLE PRECISION E,E0
c
c     .. arrays in common
      DOUBLE PRECISION EE1(0:M2D,KPOIBZ),EE2(0:M2D,KPOIBZ),
     +                 WW1(0:M2D,KPOIBZ),WW2(0:M2D,KPOIBZ)
      INTEGER II1(0:M2D,KPOIBZ),II2(0:M2D,KPOIBZ)
      COMMON /EE/ EE1,EE2,WW1,WW2,II1,II2 ! index arrays for eigenvalue
                                          ! analysis

      INTEGER NBASIS(KPOIBZ),NMIN(KPOIBZ),NMAX(KPOIBZ)
      COMMON /NBASIS / NBASIS,NMIN,NMAX

      DOUBLE PRECISION ZERO,TH
      integer J
      
      DATA ZERO / 0.0D0 /
     +     TH   / 5.0D2 /

      SAVE
c ------------------------------------------------------------------------
      MAX = NMAX(IP) - NMIN(IP)
      IO = 0

      IF (MAX.GT.M2D .OR. MAX.LT.0) THEN
        write(6,*) 'IP,NMIN,NMAX :',IP,NMIN(IP),NMAX(IP)
        STOP 'EREST1'
      END IF

      E = E1
      DO 10 J = MAX,0,-1
        E0 = EE1(J,IP)
        IF (E0.LT.TH .AND. E0.GE.E) N = J
 10   END DO

      E = -E1
      E0 = EE2(N,IP)
      IF (E0.LT.TH .AND. E0.GE.E) IO = 1

      RETURN
      END 
