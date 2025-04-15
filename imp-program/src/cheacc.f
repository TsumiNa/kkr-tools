      SUBROUTINE CHEACC(VISP,VINS,V,vspsme,vspsmo,ALPHA,NATPS,NATYP,
     $     NSPIN,INS,IRMIN,IRC,IPF,LMPOT,lsmear)
      Implicit None
c-----------------------------------------------------------------------
c      this program uses chebycheff acceleration mixing formular to
c      speed up convergence, according to paper h.akai and p.h.dederichs
c
c                          stefan bluegel , kfa , may 1987
c
c      modified for general potential  , b. drittler , aug. 1988
c
c-----------------------------------------------------------------------
c
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER IRMD,IRNSD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NTIRD
c     now also the smeared spherical potential is mixed
c      PARAMETER (NTIRD= (IRMD+ (IRNSD+1)* (LMPOTD-1))*NSPIND*NTPERD)
      PARAMETER (NTIRD= (2*IRMD+ (IRNSD+1)* (LMPOTD-1))*NSPIND*NTPERD)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALPHA
      INTEGER INS,IPF,LMPOT,NATPS,NATYP,NSPIN,lsmear
C     ..
C     .. Array Arguments ..
      REAL*8 V(IRMD,LMPOTD,*),VINS(IRMIND:IRMD,LMPOTD,*),
     +     VISP(IRMD,*),vspsmo(irmd,*),vspsme(irmd,*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      REAL*8 BETAC,ETA,F1,F2,RMIXIV
      INTEGER I,IH,IP,IR,IRC1,IRMIN1,IS,ISAVE,JI,LM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Save statement ..
      SAVE ISAVE
C     ..
C     .. Arrays in Common ..
      REAL*8 DVM(NTIRD),VM(NTIRD),VM1(NTIRD)
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

c     next for smeared spherical potential
          If ( lsmear .gt. 0 ) Then
            Do ir = 1, irc1
              ji = ji + 1
              vm(ji) = vspsme(ir,ip)
              dvm(ji) = rmixiv* (vspsmo(ir,ip)-vm(ji))
            End Do
          End If

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

      ELSE
        IF (ISAVE.NE.0) WRITE (IPF,FMT=9000)
c
c-----> determine mixing parameter betac
c
*       the factor 0.9 is an empirical one
*       concerning the smallest eigenvector
        ETA = 1.D0 - ALPHA*0.9D0
        BETAC = (1.D0/ETA-SQRT(1.D0/ETA/ETA-1.D0))**2
        IF (ISAVE.EQ.0) BETAC = 0.0D0
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

c     next for smeared spherical potential
            If ( lsmear .gt. 0 ) Then
              Do ir = 1, irc1
                ji = ji + 1
                vspsmo(ir,ip) = dvm(ji)
              End Do
            End If

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
        ISAVE = ISAVE + 1
      END IF


 9000 FORMAT ('       chebychev   used    !')
      END
