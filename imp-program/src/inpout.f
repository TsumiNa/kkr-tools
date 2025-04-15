      SUBROUTINE INPOUT(FILNAM,IRESIST,IFILE,IPE,IWTMAT,NSPIN,IRM,INS,
     +                  ICST,KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC,KTE,
     +                  KPRE,KSPH,KEFG,KVMAD,KF,IGGA,NATREF,NATPER,ICUT,
     +                  KSHAPE,KMOLD,IREF,NTCELL,IRNS,ITCLST,IMIX,IEF,
     +                  STRMIX,FCM,QBOUND,BRYMIX,EF,HFIELD,VBC,VCONSTH,
     +                  KLATR,KSYMM,LKONV,KSYMMAT,KESYM,OCCUP,ASARY,RM,
     +                  RM1,KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD,
     +                  K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST,IOBROY,
     +                  IOBROG,IOFILE,STEPL,DELTAMAX,T_DEBYE,DIMASS,
     +                  EPS_MD,KATDYN,SOCCUP,KOCCUP,QOCCUPH)

C     If mdyncs.f ends succesfully this subroutine
C     writes out a new inputcard "- input_kkr".
C     This inputcard can be used to make new shapes and
C     to interpolate the potentials.
C     For detailed description of the parameters have a
C     look at reainp.f.
C                                     good luck
C                                     Holger Hoehler , 1999

      IMPLICIT NONE
C     .. Parameters ..
c
      INTEGER NATYPD,NTREFD,NATOMD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD,
     +          NSPIND=2)
      INTEGER NREPD
      PARAMETER (NREPD=4)
C     ..
C     .. Local Scalars ..
      REAL*8 DELTAMAX,EPS_MD,FPI,PI,RFPI,STEPL,T_DEBYE
      INTEGER I,I1,I2,IOBROG,IOBROY,IOFILE,IOWFCT,ITDMD,KMDYN,KPRI,
     +        KTEST,K_ADPT_DT,K_EPS_MD,K_PRPFX,K_SET_F_M,N,NATOM,NFATOM,
     +        NFATOM_DYN,NREP
C     ..
C     .. Local Arrays ..
      REAL*8 DIMASS(NATOMD)
      INTEGER KATDYN(NTPERD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
CASYM
C     .. Scalar Arguments ..
      REAL*8 BRYMIX,FCM,HFIELD,QBOUND,STRMIX,VCONSTH
      INTEGER ICST,ICUT,IEF,IFILE,IGGA,IMIX,INS,IPE,IRESIST,IRM,ITCLST,
     +        IWTMAT,KCOR,KEFG,KESYM,KF,KFSR,KHFELD,KHYP,KLATR,KMOLD,
     +        KPRE,KSHAPE,KSPH,KSYMM,KSYMMAT,KTE,KVMAD,KVREL,KWS,KXC,
     +        LKONV,NATPER,NATREF,NSPIN,KOCCUP
C     ..
C     .. Array Arguments ..
      REAL*8 EF(NSPIND),RM(3,NATOMD),RM1(3,NATOMD),VBC(2)
      INTEGER ASARY(NREPD,NSPIND),IREF(NTPERD),IRNS(NATYPD),
     +        NTCELL(NATYPD),OCCUP(NATYPD),SOCCUP(NSPIND,NREPD),
     +        QOCCUPH(NSPIND,NREPD)
      CHARACTER*60 FILNAM(7)
C     ..

c
c------------ definition of input parameter -----------
c
      PI = 4.D0*DATAN(1.D0)
      FPI = 4.0D0*PI
      RFPI = SQRT(FPI)
      IOWFCT = 50
      IGGA = 0
      LKONV = 0
CASYM
CASYM
      OPEN (700,FILE='input_kkr')
      FILNAM(2) = 'intpot_in'
      FILNAM(6) = 'shapefun_in'
      DO I = 1,7
          WRITE (700,FMT=9320) FILNAM(I)
      END DO
      WRITE (700,FMT=9120) IRESIST,IFILE,IPE,IWTMAT
      WRITE (700,FMT=9120) NSPIN,IRM,INS,ICST
      WRITE (700,FMT=9120) KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC
      WRITE (700,FMT=9120) KTE,KPRE,KSPH,KEFG,KVMAD,KF,IGGA
      WRITE (700,FMT=9120) NATREF,NATPER,ICUT,KSHAPE,KMOLD
      WRITE (700,FMT=9120) (IREF(I),I=1,NATPER)
      WRITE (700,FMT=9120) (NTCELL(I),I=1,NATREF+NATPER)
      IF (KMOLD.EQ.1) THEN
          WRITE (700,FMT=9120) KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD
          WRITE (700,FMT=9120) K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST
          WRITE (700,FMT=9120) IOBROY,IOBROG,IOFILE
          WRITE (700,FMT=9280) STEPL,DELTAMAX,T_DEBYE,DIMASS(1)
          WRITE (700,FMT=9310) EPS_MD
          WRITE (700,FMT=9120) (KATDYN(I),I=1,NATPER)
      END IF
      IF (INS.GT.0) WRITE (700,FMT=9120) (IRNS(I),I=1,NATREF+NATPER)
      WRITE (700,FMT=9120) ITCLST,IMIX,IEF
      WRITE (700,FMT=9330) STRMIX,FCM,QBOUND
      WRITE (700,FMT=9090) BRYMIX,EF(1),EF(2)
      WRITE (700,FMT=9090) HFIELD,VBC(1),VCONSTH
      WRITE (700,FMT=9120) KLATR,KSYMM,LKONV,KSYMMAT,KESYM,KOCCUP
CBZ-INTEGRATION CZ
      IF (KLATR.EQ.10) THEN
          WRITE (700,FMT=9120) (OCCUP(I),I=1,NATPER)
      END IF
CASYM ----> This is read in, because NREP is used in next READ
      IF (KSYMMAT.NE.0) THEN
          DO I2 = 1,NSPIN
              DO I1 = 1,NREPD
                  IF (ASARY(I1,I2).gt.0) ASARY(I1,I2)=1
              END DO
              WRITE (700,FMT=9120) (ASARY(I1,I2),I1=1,NREPD)
          END DO
      END IF
      IF (KOCCUP.EQ.1.AND.KSYMMAT.NE.0) THEN
          DO I2 = 1,NSPIN
              WRITE (700,FMT=9120) (QOCCUPH(I2,I1),I1=1,NREPD)
          END DO
      END IF
      IF (KLATR.GT.0) THEN
          DO N = 1,NATOMD
              WRITE (700,FMT=9280) (RM(I,N),I=1,3), (RM1(I,N),I=1,3)
          END DO

      END IF


c
      RETURN
 9090 FORMAT (3f10.6)
 9120 FORMAT (8i5)
 9121 FORMAT (8f5.0)
 9280 FORMAT (6F10.5)
 9300 FORMAT (i5)
 9310 FORMAT (4E10.5)
 9320 FORMAT (A50)
 9330 FORMAT (2f10.6,1E10.1)
      END
