      SUBROUTINE CUTSHELL(ICUT,IREF,IRWS,LPOT,NSPIN,NATREF,NATYP,RHO2NS)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     if less shells should be taken into account this subroutine
c     has to be used .
c     the charge density and the spin density (in the case of spin
c     polarization) of the representive atoms which should be not
c     taken into account are replaced by those of the host atoms .
c
c                               b.drittler   sept. 1987
c-----------------------------------------------------------------------
c
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,LPOTD
      PARAMETER (irmd=1484,LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER ICUT,LPOT,NATREF,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      REAL*8 RHO2NS(IRMD,LMPOTD,NATYPD,*)
      INTEGER IREF(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IATYP,IRWS1,LM,LMPOT,NREF,NSTART
C     ..
C     .. Save statement ..
      SAVE
C     ..
      LMPOT = (LPOT+1)**2
      NSTART = NATYP - ICUT + 1
c
      IF (NSTART.LE.NATYP) THEN
        WRITE (6,FMT=9000) ICUT
c---> loop over reference atoms
        DO 40 IATYP = NSTART,NATYP
          IRWS1 = IRWS(IATYP)
          NREF = IREF(IATYP-NATREF)
          DO 30 LM = 1,LMPOT
            DO 10 I = 1,IRWS1
              RHO2NS(I,LM,IATYP,1) = RHO2NS(I,LM,NREF,1)
              RHO2NS(I,LM,IATYP,NSPIN) = RHO2NS(I,LM,NREF,NSPIN)
   10       CONTINUE

            DO 20 I = IRWS1 + 1,IRMD
              RHO2NS(I,LM,IATYP,1) = 0.0D0
              RHO2NS(I,LM,IATYP,NSPIN) = 0.0D0
   20       CONTINUE
   30     CONTINUE
   40   CONTINUE
      END IF


 9000 FORMAT (13x,'the last',i2,'shells of the cluster are cutted off')
      END
