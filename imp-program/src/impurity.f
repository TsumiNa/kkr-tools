      PROGRAM IMPURITY
        IMPLICIT NONE
c   use  imslf90
C     .. Local Scalars ..
      INTEGER I,IHANDLE
C     ..
C     .. Local Arrays ..
      CHARACTER*60 FILNAM(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL FPIMPU
C     ..
      DO 10 I = 1,7
        READ (5,FMT=9000) FILNAM(I)
   10 CONTINUE
c     OPEN (14,FILE=FILNAM(1),FORM='unformatted')
c
      OPEN (12,FILE=FILNAM(2),FORM='formatted')
      OPEN (89,FILE=FILNAM(3),FORM='unformatted')
C   OPEN (89,FILE=FILNAM(3),FORM='formatted')
      OPEN (71,FILE=FILNAM(4),FORM='unformatted')
      OPEN (35,FILE=FILNAM(5),FORM='formatted')
      OPEN (19,FILE=FILNAM(6),FORM='formatted')
      OPEN (40,FILE=FILNAM(7),FORM='formatted')
      OPEN (50,FILE='wfct',FORM='unformatted')
      OPEN (51,FILE='wfct1',FORM='unformatted')
      OPEN (80,FILE=FILNAM(1),FORM='unformatted',access='stream')
      OPEN (20,FILE='broy',FORM='unformatted')
c     CALL FXDRINI
      IHANDLE = 80
      CALL FPIMPU(IHANDLE,FILNAM)
      STOP
 9000 FORMAT (A50)
      END
