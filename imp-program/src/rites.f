      SUBROUTINE RITES(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,
     +     ITITLE,R,DRDI,VM2Z,vspsme,IRWS,A,B,TXC,KXC,INS,IRNS,
     +     LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +     LCORE,NCORE,IGGA,lsmear,alpha)
      Implicit None
c-----------------------------------------------------------------------
c      this subroutine stores in 'ifile' the necessary results
c      (potentials e.t.c.) to start self-consistency iterations
c
c      modified for the full potential case - if ins .gt. 0 there
c       is written a different potential card
c       if the sum of absolute values of an lm component of vins (non
c       spher. potential) is less than the given rms error qbound this
c       component will not be stored .
c
c        (see to subroutine start , where most of the arrays are
c         described)
c
c                            modified by b. drittler  aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER IRMD,IRNSD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALAT,QBOUND,alpha
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN,IGGA,lsmear
C     ..
C     .. Array Arguments ..
      REAL*8 A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI(2),
     +     R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),vspsme(irmd,*),
     +     VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*24 TXC(*)
C     ..
C     .. Local Scalars ..
      REAL*8 A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1
      INTEGER I,ICORE,IH,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,LMPOT,
     +        NCORE1,NR
C     ..
C     .. Local Arrays ..
      REAL*8 DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD),
     $     vm2sme(irmd)
      INTEGER LCORE1(20)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
      ISAVE = 1
c
      LMPOT = (LPOT+1)* (LPOT+1)
      DO 60 IH = 1,NATYP
        DO 50 IS = 1,NSPIN
          IF (IS.EQ.NSPIN) THEN
            SIGN = 1.0D0

          ELSE
            SIGN = -1.0D0
          END IF
          IP = NSPIN* (IH-1) + IS

          RMT1 = RMT(IH)
          RMTNW1 = RMTNEW(IH)
          Z1 = Z(IH)
          RMAX = RWS(IH)
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS(IH)

          ELSE
            NR = IRC(IH)
          END IF

          IRNS1 = IRNS(IH)
          IRMIN = NR - IRNS1
          A1 = A(IH)
          B1 = B(IH)
          NCORE1 = NCORE(IP)
c
          DO 10 J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
c
c---> store only lm=1 component of the potential
c
            VM2ZA(J) = VM2Z(J,IP)
            vm2sme(J) = vspsme(J,IP)
   10     CONTINUE
c
          IF (NCORE1.GE.1) THEN
c
            DO 20 J = 1,NCORE1
              LCORE1(J) = LCORE(J,IP)
              ECORE1(J) = ECORE(J,IP)
   20       CONTINUE
          END IF
c
c
          IF (IGGA.EQ.0) THEN
          WRITE (IFILE,FMT=9000) (ITITLE(I,IP),I=1,7),TXC(KXC+1)
          ELSE
          WRITE (IFILE,FMT=9000) (ITITLE(I,IP),I=1,7),TXC(IGGA+3)
          END IF
          WRITE (IFILE,FMT=9010) RMT1,ALAT,RMTNW1
          WRITE (IFILE,FMT=9020) Z1,RMAX,EFERMI(IS),VBC(IS)
          WRITE (IFILE,FMT=9030) NR,A1,B1,NCORE1
          IF (NCORE1.GE.1) WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +        ECORE1(ICORE),ICORE=1,NCORE1)
c
          IF (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) THEN
c
c---> store only the spherically averaged potential (in mt or as - case)
c        this is done always for the host
c
            WRITE (IFILE,FMT=9050) (RA(IR),DRADI(IR),VM2ZA(IR),IR=1,NR)

          ELSE
c
c---> store the full potential , but the non spherical contribution
c         only from irns1 up to irws1 ;
c     remember that the lm = 1 contribution is multiplied
c         by a factor 1/sqrt(4 pi)
c
            WRITE (IFILE,FMT=9060) NR,IRNS1,LMPOT,ISAVE
            WRITE (IFILE,FMT=9070) (VM2ZA(IR),IR=1,NR)
c
c     lsmear=3: no smeared potential for host
            If ( (lsmear.eq.3) .and. (ih.ge.natps) ) Then
              Write (6,*) 'Smeared spherical potential saved on disk,',
     $             ' atom i = ',ih
              Write (ifile,fmt=9080) nr,irns1,lmpot,isave,alpha
              Write (ifile,fmt=9070) (vm2sme(ir),ir=1,nr)
            End If
c
c     lsmear=4: host+clu have smeared potential
            If ( lsmear .eq. 4 ) Then
              Write (6,*) 'Smeared spherical potential saved on disk,',
     $             ' atom i ',ih
              Write (ifile,fmt=9080) nr,irns1,lmpot,isave,alpha
              Write (ifile,fmt=9070) (vm2sme(ir),ir=1,nr)
            End If

            IF (LPOT.GT.0) THEN
              LMNR = 1
              DO 40 LM = 2,LMPOT
                SUM = 0.0D0
                DO 30 IR = IRMIN,NR
                  RV = VINS(IR,LM,IP)*RA(IR)
                  SUM = SUM + RV*RV*DRADI(IR)
   30           CONTINUE

                IF (SQRT(SUM).GT.QBOUND) THEN
                  LMNR = LMNR + 1
                  WRITE (IFILE,FMT=9060) LM
                  WRITE (IFILE,FMT=9070) (VINS(IR,LM,IP),IR=IRMIN,NR)
                END IF

   40         CONTINUE
c
c---> write a one to mark the end
c
              IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
            END IF

          END IF

   50   CONTINUE
   60 CONTINUE



 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i5,/,2d15.8,/,i2)
 9040 FORMAT (i5,1p,d15.6)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
 9080 FORMAT (4i5,' alpha=',1f10.7)
      END
