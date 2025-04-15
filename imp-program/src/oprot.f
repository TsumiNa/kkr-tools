      SUBROUTINE OPROT(IOPER,NATPER,ND,NSHELL,RM)
      IMPLICIT NONE
c
c---> attention changed version : call rotcub is now in the main program
c
C     .. Scalar Arguments ..
      INTEGER NATPER
C     ..
C     .. Array Arguments ..
      REAL*8 RM(3,*)
      INTEGER IOPER(*),ND(48,3,3),NSHELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 ARG,X,Y,Z
      INTEGER K,N,NIR,NMAX,NMIN,NS,NST
      INTEGER C4V(9),NC4V,C3V(6),NC3V,NTD,TD(24)
      INTEGER NC2V,C2V(4) 
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
c changed on 20/11/1996 ny nikos
c this is a test to work for C4V only
      DATA C4V/1,12,10,18,40,45,46,48,41/
      DATA NC4V/9/
CO
CO    Changed to work for C3v
CO
      DATA C3V/1,2,3,46,47,48/
      DATA NC3V/6/
      DATA C2V    /1,18,45,46/
      DATA NC2V   /4/
      DATA TD /1,2,3,4,5,6,7,8,9,16,17,18,34,35,36,37,
     +             38,39,43,44,45,46,47,48/
      
      DATA NTD/24/
CO                                 ------->C4V
c      WRITE(6,*) ' SUB OPROT , CHANGED FOR C4V '
c 
      NIR = 1
      DO 40 NS = 1,NATPER
        NST = NIR
        NIR = NIR + NSHELL(NS)
        NMIN = NST
        NMAX = NST + NSHELL(NS) - 1
        DO 30 N = NMIN,NMAX
c
          DO 10 K = 1,48
c this was added
cCO             ------->NTD
c          DO 10 K1 = 1,NC2V
cCO     -------->TD
c            K = C2V(K1)
c
cCO             ------->NC4V
c          DO 10 K1 = 1,NC4V
cCO     -------->C4V        
c            K = C4V(K1)
c  this was added
            X = ND(K,1,1)*RM(1,N) + ND(K,1,2)*RM(2,N) +
     +          ND(K,1,3)*RM(3,N)
            Y = ND(K,2,1)*RM(1,N) + ND(K,2,2)*RM(2,N) +
     +          ND(K,2,3)*RM(3,N)
            Z = ND(K,3,1)*RM(1,N) + ND(K,3,2)*RM(2,N) +
     +          ND(K,3,3)*RM(3,N)
            ARG = ABS(X-RM(1,NST)) + ABS(Y-RM(2,NST)) + ABS(Z-RM(3,NST))
            IF (ARG.LE.0.01D0) GO TO 20
   10     CONTINUE
          GO TO 50
c
c---> determine the operation k of the oh-group which transforms the
c     atom n into its representive shell atom . ( active rotation of
c     the lattice vector !)
c
   20     IOPER(N) = K
   30   CONTINUE
   40 CONTINUE
      RETURN

   50 WRITE (6,FMT=9000)
      STOP


 9000 FORMAT (/,33x,'stop 2 in subroutine oprot')
      END
