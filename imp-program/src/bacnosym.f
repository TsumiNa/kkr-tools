      SUBROUTINE BACNOSYM(DTB,GMAT,LMAXSQ,NREP,NP)
      IMPLICIT NONE
C
C    This is for the case of no symettry
c
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD)
      INTEGER NSEC
      PARAMETER (nsec=689)
      INTEGER LMAXD,LMX,LMMAXD
      PARAMETER (lmaxd=4,LMX=LMAXD+1,LMMAXD=LMX**2)
      INTEGER JCOEFF
      PARAMETER (JCOEFF=187261)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,NP,NREP
      integer ntest
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DTB(LMMAXD,LMMAXD,*),GMAT(NSEC,NSEC)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO
      INTEGER I1,I2,I3,ICG,J,J1,J2,NR
CO--------ryu-------
      integer t_IN
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..

      ntest = ntperd*lmaxsq
      if (ntest.ne.nsec) then
       write(6,*) 'Look at backnosym nsec <> ntperd*lmaxsq ',nsec,
     &             ntest
       stop 'backnosym'
      end if
        DO 70 I3 = 1,NTPERD
          t_IN = LMAXSQ*(I3-1)
          DO 60 I2 = 1,LMAXSQ
            DO 50 I1 = 1,LMAXSQ
              DTB(I2,I1,I3) = GMAT(t_IN+I2,t_IN+I1)
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE

      END
