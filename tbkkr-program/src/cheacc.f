c ************************************************************************
      SUBROUTINE CHEACC(VISP,VINS,V,ALPHA,NATPS,NATYP,NSPIN,INS,IRMIN,
     +                  IRC,IPF,LMPOT,IE) 
      implicit none
c ************************************************************************
c      this program uses chebycheff acceleration mixing formular to
c      speed up convergence, according to paper h.akai and p.h.dederichs,
c      j.phys.c, 18(1985), pp. 2455-2460
c
c                          stefan bluegel , kfa , may 1987
c
c      modified for general potential  , b. drittler , aug. 1988
c
c      modified p.zahn, may 1996
c-----------------------------------------------------------------------
c
C     .. Parameters ..
      INCLUDE 'inc.fi'
      INTEGER NTPERD
      PARAMETER (NTPERD=NATYPD-NTREFD)
C      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
C      PARAMETER (NATYPD=1,NTREFD=0,NTPERD=NATYPD-NTREFD,NSPIND=2)
C      INTEGER IRMD,IRNSD,LPOTD
C      PARAMETER (IRMD=1484,IRNSD=508,LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NTIRD
      PARAMETER (NTIRD= ( IRMD*NTPERD + 
     +                    (IRNSD+1)*(LMPOTD-1)*NATYPD )*NSPIND)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA
      INTEGER IE,INS,IPF,LMPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION V(IRMD,LMPOTD,*),VINS(IRMIND:IRMD,LMPOTD,*),
     +                 VISP(IRMD,*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BETAC,ETA,F1,F2,RMIXIV
      INTEGER I,IH,IP,IR,IRC1,IRMIN1,IS,ISAVE,JI,LM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Save statement ..
      SAVE ISAVE
C     ..
C     .. Arrays in Common ..
      DOUBLE PRECISION DVM(NTIRD),VM(NTIRD),VM1(NTIRD)
C     ..
C     .. Common blocks ..
      COMMON /CHEBY/DVM,VM,VM1
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Data statements ..
      DATA ISAVE/0/
C     ..
c--->
c     initialize vm1
c---<
      IF (ISAVE.EQ.0) THEN
        DO 10 I = 1,NTIRD
          VM1(I) = 0.0D0
   10   CONTINUE
      END IF
c
c---->  map potentials of each mt-sphere in one single vector
c
      RMIXIV = 1.0D0/ALPHA
      JI = 0
      DO 60 IS = 1,NSPIN
        DO 50 IH = NATPS,NATYP
          IP = NSPIN* (IH-1) + IS

          IRC1 = IRC(IH)
          DO 20 IR = 1,IRC1
            JI = JI + 1
            VM(JI) = VISP(IR,IP)
            DVM(JI) = RMIXIV* (V(IR,1,IP)-VM(JI))
   20     CONTINUE
          IF (INS.GT.0 .AND. LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IH)
            DO 40 LM = 2,LMPOT
              DO 30 IR = IRMIN1,IRC1
                JI = JI + 1
                VM(JI) = VINS(IR,LM,IP)
                DVM(JI) = RMIXIV* (V(IR,LM,IP)-VM(JI))
   30         CONTINUE
   40       CONTINUE
          END IF

   50   CONTINUE
   60 CONTINUE
c
      IF (JI.GT.NTIRD) THEN
        CALL RCSTOP('cheacc  ')

      ELSE                          ! (JI.GT.NTIRD)
c
c-----> determine mixing parameter betac
c
c       the factor 0.9 is an empirical one
c       concerning the smallest eigenvector
c
c        ETA = 1.D0 - ALPHA*0.9D0
        ETA = 1.D0 - ALPHA*4.0D0
        BETAC = (1.D0/ETA-SQRT(1.D0/ETA/ETA-1.D0))**2

        IF (ISAVE.NE.0 .AND. IE.EQ.ISAVE+1) THEN
c        IF (ISAVE.NE.0) THEN
          WRITE (IPF,FMT=9000)
        ELSE
          BETAC = 0.0D0
        END IF
c
c-----> perform mixing
c
        F1 = 1.D0 + BETAC
        F2 = ALPHA*F1
        DO 70 IR = 1,JI
          DVM(IR) = F1*VM(IR) + F2*DVM(IR) - BETAC*VM1(IR)
   70   CONTINUE
        DO 80 IR = 1,JI
          VM1(IR) = VM(IR)
   80   CONTINUE
c
c-----> map back to each muffin tin cell
c
        JI = 0
        DO 130 IS = 1,NSPIN
          DO 120 IH = NATPS,NATYP
            IP = NSPIN* (IH-1) + IS

            IRC1 = IRC(IH)
            DO 90 IR = 1,IRC1
              JI = JI + 1
              V(IR,1,IP) = DVM(JI)
   90       CONTINUE
            IF (INS.GT.0 .AND. LMPOT.GT.1) THEN
              IRC1 = IRC(IH)
              IRMIN1 = IRMIN(IH)
              DO 110 LM = 2,LMPOT
                DO 100 IR = IRMIN1,IRC1
                  JI = JI + 1
                  V(IR,LM,IP) = DVM(JI)
  100           CONTINUE
  110         CONTINUE
            END IF

  120     CONTINUE
  130   CONTINUE
        ISAVE = IE
      END IF                        ! (JI.GT.NTIRD)

      RETURN

 9000 FORMAT ('       chebychev   used    !')
      END
