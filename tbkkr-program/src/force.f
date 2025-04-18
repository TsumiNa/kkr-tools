      SUBROUTINE FORCE(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,DRDI,
     +                 IRWS)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m
c     from a given non spherical charge density at the nucleus site r
c     with core correction (coulomb contribution)
 
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMD,*),
     +       RHOC(IRMD,*),V(IRMD,LMPOTD,*)
      INTEGER IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DV,FAC,PI,RWS,VINT1
      INTEGER I,IATYP,IPOT,IRWS1,ISPIN,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION FLMH(-1:1,NATYPD),V1(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
c
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DSQRT
C     ..
      PI = 4.D0*DATAN(1.D0)
      FAC = DSQRT((4.0D0*PI)/3.0D0)
      IF (LMAX.LT.1) THEN
         WRITE (6,FMT=9000)
         STOP
 
      END IF
c
c---> loop over rep. atoms
c
      DO 10 IATYP = NSTART,NEND
c
c
         IRWS1 = IRWS(IATYP)
         RWS = R(IRWS1,IATYP)
c
c
 
         DO 20 M = -1,1
            LM = 2 + M + 1
c
c---> initialize v1
c
            DO 30 I = 1,IRWS1
               V1(I) = 0.0D0
   30       CONTINUE
c
            DO 40 ISPIN = 1,NSPIN
c
c---> determine the right potential numbers
c
               IPOT = NSPIN* (IATYP-1) + ISPIN
c
c---> determine the derivative of the potential using a 5-point formular
c
               DV = (-3.0D0*V(1,LM,IPOT)-10.0D0*V(2,LM,IPOT)+
     +            18.0D0*V(3,LM,IPOT)-6.0D0*V(4,LM,IPOT)+V(5,LM,IPOT))/
     +              (12.0D0*DRDI(2,IATYP))
c
               V1(2) = RHOC(2,IPOT)* (2.0D0*V(2,LM,IPOT)/R(2,IATYP)+DV)/
     +                 (4.0D0*PI) + V1(2)
c
               DO 50 I = 3,IRWS1 - 2
c
                  DV = (V(I-2,LM,IPOT)-V(I+2,LM,IPOT)+
     +                 8.0D0* (V(I+1,LM,IPOT)-V(I-1,LM,IPOT)))/
     +                 (12.0D0*DRDI(I,IATYP))
c
                  V1(I) = RHOC(I,IPOT)* (2.0D0*V(I,LM,IPOT)/R(I,IATYP)+
     +                    DV)/ (4.0D0*PI) + V1(I)
   50          CONTINUE
c
               DV = (-V(IRWS1-4,LM,IPOT)+6.0D0*V(IRWS1-3,LM,IPOT)-
     +              18.0D0*V(IRWS1-2,LM,IPOT)+10.0D0*V(IRWS1-1,LM,IPOT)+
     +             3.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1-1,IATYP))
               V1(IRWS1-1) = RHOC(IRWS1-1,IPOT)*
     +                       (2.0D0*V(IRWS1-1,LM,IPOT)/R(IRWS1-1,IATYP)+
     +                       DV)/ (4.0D0*PI) + V1(IRWS1-1)
c
               DV = (3.0D0*V(IRWS1-4,LM,IPOT)-16.0D0*V(IRWS1-3,LM,IPOT)+
     +              36.0D0*V(IRWS1-2,LM,IPOT)-48.0D0*V(IRWS1-1,LM,IPOT)+
     +              25.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1,IATYP))
c
               V1(IRWS1) = RHOC(IRWS1,IPOT)*
     +                     (2.0D0*V(IRWS1,LM,IPOT)/R(IRWS1,IATYP)+DV)/
     +                     (4.0D0*PI) + V1(IRWS1)
   40       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
            FLMH(M,IATYP) = FAC*FLM(M,IATYP)
            FLMC(M,IATYP) = -FAC*VINT1
            FLM(M,IATYP) = FLMH(M,IATYP) + FLMC(M,IATYP)
c
 
   20    CONTINUE
c
c
   10 CONTINUE
c
 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least !')
 
      END
