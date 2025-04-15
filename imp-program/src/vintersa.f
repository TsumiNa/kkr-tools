      SUBROUTINE VINTERSA(AMAT,BMAT,CMOM,IREF,LMAX,NATPER,NATREF,NSPIN,
     +     NSTART,NEND,V,Z,R,IRWS,IRCUT,IPAN,KSHAPE,
     +     CMINST,AVMAD,BVMAD,irm,ntim)
      Implicit None
C-----------------------------------------------------------------------
C     CALCULATE THE INTERCELL-POTENTIALS AND ADD THESE TO THE POTEN-
C     TIAL V  (IN THE SPIN-POLARIZED CASE FOR EACH SPIN-DIRECTION
C     THE INTERCELL-POTENTIAL IS THE SAME . )
C     IT USES THE STRUCTURE DEPENDENT MATRICES AMAT AND BMAT WHICH
C     ARE CALCULATE ONCE IN THE SUBROUTINE AMN .
C     THE CHARGE-MOMENTS ARE CALCULATED IN THE SUBROUTINE VINTRA2 ,
C     THEREFORE VINTRA2 HAS TO BE CALLED FIRST .
C     THE INTERCELL-POTENTIAL IS EXPANDED INTO SPHERICAL HARMONICS .
C     THE LM-TERM OF THE INTERCELL-POTENTIAL V OF THE REPRESENTIVE
C     ATOM I IS GIVEN BY
C
C      V(R,LM,I) =  (-R)**L * {AMAT(I1,I2,LM,L'M')*CMOM(I2,L'M')
C                                               +BMAT(I1,I2,LM)*Z(I2)}
C
C     SUMMED OVER I2 (ALL SHELLS) AND L'M' .    (I1=I-NATREF)
C             (SEE NOTES BY B.DRITTLER)
C
C     IN CASE OF SHAPE CORRECTION THE MADELUNG POTENTIAL OF THE HOST
C        IS TAKEN INTO ACCOUNT . IN ALL OTHER CASE THE MADELUNG POTEN-
C        TIAL OF THE HOST IS SET TO BE ZERO .
C        AS ACTUAL VALUES FOR Z AND CMOM THE DIFFERENCES BETWEEN THE
C        VALUES OF THE  REPRESENTIVE ATOMS AND THOSE OF THE REFERENCES
C        ARE USED .
C
C     ATTENTION : THE FIRST INDEX OF CMOM (MOMENT OF THE CHARGE
C                 DENSITY - CALCULATED IN VINTR2) AND OF Z (NUCLEAR
C                 CHARGE OF THE ATOMS IN THE SHELL) IS IN THE PROGRAM
C                 DIFFERENT DEFINED , THERE ONE HAS TO USE :
C                     CMOM(NATREF+I2,L'M')  AND  Z(NATREF+I2)
C
C                               B.DRITTLER   JUNE 1987
C-----------------------------------------------------------------------
C     .. PARAMETERS ..
      INTEGER NTREFD,NTPERD,NATOMD
c next line changed by T.Korhonen Apr 97
c      PARAMETER(NTREFD=1,NTPERD=19,NATOMD=79)
      PARAMETER(NTREFD=1,NATOMD=102)
      INTEGER NSPIND
      PARAMETER(NSPIND=2)
      INTEGER LPOTD
      PARAMETER(lpotd=8)
      INTEGER IPAND
      PARAMETER(IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NATYPD,NPOTD
c next line changed by T.Korhonen Apr 97
c      PARAMETER (NATYPD=NTPERD+NTREFD,NPOTD=NSPIND*NATYPD)
      PARAMETER (NATYPD=38,NTPERD=NATYPD-NTREFD,NPOTD=NSPIND*NATYPD)
C     ..
C     .. SCALAR ARGUMENTS ..
      INTEGER KSHAPE,LMAX,NATPER,NATREF,NEND,NSPIN,NSTART,NREF2,
     $     irm, ntim
C     ..
C     .. ARRAY ARGUMENTS ..
      REAL*8 AMAT(NTPERD,NTPERD,LMPOTD,LMPOTD),
     +     AVMAD(NTREFD,NTREFD,LMPOTD,LMPOTD),
     +     BMAT(NTPERD,NTPERD,LMPOTD),BVMAD(NTREFD,NTREFD,LMPOTD),
     +     CMINST(LMPOTD,NATYPD),CMOM(LMPOTD,NATYPD),
     +     R(IRM,NATYPD),V(IRM,LMPOTD,NPOTD),Z(NATYPD)
      INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD),IREF(NTPERD),
     +        IRWS(NATYPD)
C     ..
C     .. LOCAL SCALARS ..
      REAL*8 AC,PI,SUM
      INTEGER I,I1,I2,IATYP,IATYP2,IPOT,IRS1,ISPIN,L,LM,LM2,LMMAX,M,NREF
C     ..
C     .. LOCAL ARRAYS ..
      REAL*8 ACH(LMPOTD,NTREFD,2)
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. SAVE STATEMENT ..
      SAVE ACH
C     ..
      PI = 4.0*DATAN(1.0D0)
      LMMAX = (LMAX+1)* (LMAX+1)

       DO I1=1,25
          DO I2=1,4
             WRITE (53,*) i1,i2,CMOM(I1,I2)
          END DO
       END DO
c      STOP
      DO 10 IATYP = NSTART,NEND
C---> DETERMINE SHELL INDEX
         I1 = IATYP - NATREF
         IF (KSHAPE.NE.0) THEN
      write(6,*) nstart,nend
            IRS1 = IRCUT(IPAN(IATYP),IATYP)
         ELSE
            IRS1 = IRWS(IATYP)
         END IF
C
         DO 20 L = 0,LMAX
            DO 30 M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0E0
	       IF (KSHAPE .eq. 0) then
               IF (I1.LT.1) THEN
                        DO 53 I2 = 1,NATREF
                           AC = AC + BVMAD(IATYP,I2,LM)*Z(I2)
                           DO 63 LM2 = 1,LMMAX
                              SUM = CMOM(LM2,I2) 
                              AC = AC + AVMAD(IATYP,I2,LM,LM2)*SUM
   63                      CONTINUE
   53                   CONTINUE
                     ACH(LM,IATYP,ntim) = AC
ccc                     write(6,*) 'aa',lm,iatyp, ach(lm,iatyp,ntim)
               END IF
               end if
               IF (KSHAPE.GT.0) THEN
                  IF (I1.LT.1) THEN
C
C---> MADELUNG POTENTIAL - ONLY IN CASE SHAPE CORRECTION
C
                     IF (NATREF.EQ.1) THEN
C---> LM = 1 COMPONENT DISAPPEARS IF THERE IS ONLY ONE HOST ATOM
                        DO 40 LM2 = 2,LMMAX
C---> TAKE MOMENTS OF MT SPHERE AND INTERSTIAL
                           SUM = CMOM(LM2,1) + CMINST(LM2,1)
                           AC = AC + AVMAD(IATYP,1,LM,LM2)*SUM
   40                   CONTINUE
                     ELSE
                        DO 50 I2 = 1,NATREF
                           AC = AC + BVMAD(IATYP,I2,LM)*Z(I2)
                           DO 60 LM2 = 1,LMMAX
                              SUM = CMOM(LM2,I2) + CMINST(LM2,I2)
                              AC = AC + AVMAD(IATYP,I2,LM,LM2)*SUM
   60                      CONTINUE
   50                   CONTINUE
                     END IF
                     ACH(LM,IATYP,ntim) = AC
                  ELSE
C
C---> INTERCELL POTENTIAL IN CASE OF SHAPE CORRECTIONS
C
                     NREF2=IREF(I1)
                     DO 70 I2 = 1,NATPER
C---> DETERMINE THE REFERENCE AND REPRESENTIVE ATOM INDEX
                        NREF = IREF(I2)
                        IATYP2 = NATREF + I2
                        DO 80 LM2 = 1,LMMAX
C---> TAKE MOMENTS OF MT SPHERE AND INTERSTIAL
                           SUM = CMOM(LM2,IATYP2) - CMOM(LM2,NREF) +
     +                           CMINST(LM2,IATYP2) - CMINST(LM2,NREF)
                           AC = AC + AMAT(I1,I2,LM,LM2)*SUM
   80                   CONTINUE
                        AC = AC + BMAT(I1,I2,LM)* (Z(IATYP2)-Z(NREF))
   70                CONTINUE
                     AC = AC + ACH(LM,NREF2,ntim)
                  END IF
C
C---> INTERCELL POTENTIAL WITHOUT MADELUNG POTENTIAL IN CASE OF A.S.A.
C
               ELSE IF (I1.GE.1) THEN
                     NREF2 = IREF(I1)
                  DO 90 I2 = 1,NATPER
C---> DETERMINE THE REFERENCE AND REPRESENTIVE ATOM INDEX
                     NREF = IREF(I2)
                     IATYP2 = NATREF + I2
                     DO 100 LM2 = 1,LMMAX
                        SUM = CMOM(LM2,IATYP2) - CMOM(LM2,NREF)
                        AC = AC + AMAT(I1,I2,LM,LM2)*SUM
  100                CONTINUE
                     AC = AC + BMAT(I1,I2,LM)* (Z(IATYP2)-Z(NREF))
   90             CONTINUE
                     AC = AC + ACH(LM,NREF2,ntim)
               END IF
               IF (LM.EQ.1) WRITE (6,FMT=9000) I1, (AC/SQRT(4.E0*PI))
C
C---> ADD TO V THE INTERCELL-POTENTIAL
C
               DO 110 ISPIN = 1,NSPIN
C
C---> DETERMINE THE RIGHT POTENTIAL NUMBER
C
                  IPOT = NSPIN* (IATYP-1) + ISPIN
C
C---> IN THE CASE OF L=0 : R(1)**L IS NOT DEFINED
C
                  IF (L.EQ.0) V(1,1,IPOT) = V(1,1,IPOT) + AC
                  DO 120 I = 2,IRS1
                     V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IATYP))**L*AC
  120             CONTINUE
C
C     WRITE(6,1100) IPOT,LM,IRS1,IATYP,NSPIN
 1100 FORMAT(1X,' IPOT=',I5,' LM=',I5,' IRS1=',I5,' IATYPE=',I5,' NSPIN=
     &',I5)
C
C     IF(LM.EQ.1) WRITE(6,1200) (V(I,LM,IPOT),I=1,IRS1)
 1200 FORMAT(1X,12D10.4)
C
  110          CONTINUE
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
C
 9000 FORMAT (1X,'SPHERICALLY AVERAGED INTERCELL-POTENTIAL FOR SHELL',
     +       I2,' :',1P,D14.6)
      END
