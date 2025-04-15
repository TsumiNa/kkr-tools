      SUBROUTINE REAINP(EF,RM,IREF,IRNS,NDG,NDIM,NSHELL,NTCELL,TXC,
     +                  BRYMIX,FCM,FPI,HFIELD,PI,QBOUND,RFPI,STRMIX,VBC,
     +                  VCONSTH,ICST,ICUT,IEF,IFILE,IWTMAT,IMIX,INS,
     +                  IOWFCT,IPE,IPF,IPFE,IRESIST,IRM,ITCLST,KCOR,
     +                  KEFG,KF,KFSR,KHFELD,KHYP,KPRE,KSHAPE,KSPH,KTE,
     +                  KVMAD,KVREL,KWS,KXC,LMAX,LMAXP1,LMAXSQ,LMPOT,
     +                  LPOT,NATOM,NATPER,NATPS,NATREF,NATYP,NPTPS,
     +                  NPTREF,NREP,NSPIN,NSPPOT,KLATR,RELAX,NIMP,RM1,
     +                  DRM,KSYMM,IGGA,LKONV,LSMEAR,OCCUP,KSYMMAT,KESYM,
     +                  ASARY,SOCCUP,QOCCUP,KMOLD,KMDYN,K_PRPFX,
     +                  NFATOM_DYN,NFATOM,K_EPS_MD,K_ADPT_DT,K_SET_F_M,
     +                  IOBROY,IOBROG,IOFILE,KPRI,KTEST,ITDMD,STEPL,
     +                  DELTAMAX,T_DEBYE,EPS_MD,DIMASS,KSHAPEH,KATDYN,
     +                  IFILEH,KOCCUP,QOCCUPH)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NATOMD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD,
     +          NSPIND=2)
c      INTEGER NLSTD
ccc      PARAMETER (NLSTD=11686)
      INTEGER NSEC,NREPD
      PARAMETER (nsec=689,NREPD=4)
      INTEGER IRMD,IRNSD,LMAXD,LPOTD,LMX
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER IRID,IRIKD
      PARAMETER (irid=435,irikd=435)
c
      INTEGER IRNSKD
      PARAMETER (IRNSKD=IRIKD+ (IRNSD-IRID))
C     ..
C     .. Local Scalars ..
      INTEGER I,LMXTW1,M,NKWS,NKXC,NP,NPTPER,NQ,NT,N1,N
      CHARACTER*43 TSHAPE
C     ..
C     .. Local Arrays ..
      INTEGER NATM(NATOMD),NHOST(NATOMD),NIMP(NATOMD)
      CHARACTER*4 TSPIN(2)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,MIN,SQRT
C     ..
CASYM
      INTEGER KSYMMAT,KESYM,ASARY(NREPD,NSPIND),I1,I2,KMOLD
C     .. Scalar Arguments ..
      REAL*8 BRYMIX,FCM,FPI,HFIELD,PI,QBOUND,RFPI,STRMIX,VCONSTH
      REAL*8 STEPL,DELTAMAX,T_DEBYE,EPS_MD,DIMASS(NATOMD)
      INTEGER ICST,ICUT,IEF,IFILE,IMIX,INS,IOWFCT,IPE,IPF,IPFE,IRESIST,
     +        IRM,ITCLST,IWTMAT,KCOR,KEFG,KF,KFSR,KHFELD,KHYP,KPRE,
     +        KSHAPE,KSPH,KTE,KVMAD,KVREL,KWS,KXC,LMAX,LMAXP1,LMAXSQ,
     +        LMPOT,LPOT,NATOM,NATPER,NATPS,NATREF,NATYP,NPTPS,NPTREF,
     +        NREP,NSPIN,NSPPOT,KLATR,KSYMM,IGGA,LKONV,LSMEAR,KOCCUP,
     +        KSHAPEH,IFILEH

C     ..
C     .. Array Arguments ..
      REAL*8 EF(NSPIND),RM(3,NATOMD),VBC(2),RELAX(NTPERD)
      REAL*8 RM1(3,NATOMD),RM2(3,NATOMD),DRM0(3,NATOMD)
      REAL*8 DRM(3,NTPERD),QOCCUP(NSPIND,NREPD)
      INTEGER IREF(NTPERD),IRNS(NATYPD),NDG(NREPD),NDIM(NREPD),
     +        NSHELL(NTPERD),NTCELL(NATYPD),SOCCUP(NSPIND,NREPD),
     +        QOCCUPH(NSPIND,NREPD)
      REAL*8 NDGH(NREPD),QOCCUPH2(NSPIND,NREPD)
      INTEGER OCCUP(NATYPD)
      INTEGER KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,K_EPS_MD,K_ADPT_DT,
     +        K_SET_F_M,IOBROY,IOBROG,IOFILE,KPRI,KTEST,ITDMD,
     +        KATDYN(NTPERD)
      CHARACTER*24 TXC(5)
C     ..
C     .. Data statements ..
      DATA TSPIN/'non-','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/' non relativistic calculation              ',
     +     ' s.r.a. calculation for cluster            ',
     +     ' s.r.a. calculation                        '/
      DATA TKCOR/' frozen core approximation                 ',
     +     ' core relaxation s.r.a.  for cluster       ',
     +     ' core relaxation nonsra                    ',
     +     ' core relaxation                           '/
      DATA TINS/' spherical averaged input potential        ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/
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
      KSYMMAT = 0
      KESYM = 0
      ASARY = 0
      KOCCUP = 0
CASYM
      READ (5,FMT=9120) IRESIST,IFILE,IPE,IWTMAT
      READ (5,FMT=9120) NSPIN,IRM,INS,ICST
      READ (5,FMT=9120) KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC
      READ (5,FMT=9120) KTE,KPRE,KSPH,KEFG,KVMAD,KF,IGGA
      READ (5,FMT=9120) NATREF,NATPER,ICUT,KSHAPE,KMOLD
      READ (5,FMT=9120) (IREF(I),I=1,NATPER)
c      WRITE(6,FMT=9120) (IREF(I),I=1,NATPER)
      READ (5,FMT=9120) (NTCELL(I),I=1,NATREF+NATPER)
c If one wants to calculate the relaxed positions
c set kmold = 1 !
c In the next lines some keys and parameters are read
c in for the mdyncs.f
c For further informations one should have a look at
c this subroutine.
c
c               edited by Holger Hoehler
c
      KSHAPEH = KSHAPE
      IFILEH = IFILE
      IF (KMOLD.EQ.1) THEN
          WRITE (6,FMT=*)
          WRITE (6,FMT=*) 'NEW POSITIONS ARE CALCULATED WITH MDYNC.F !'
          READ (5,FMT=9120) KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD
          READ (5,FMT=9120) K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST
          READ (5,FMT=9120) IOBROY,IOBROG,IOFILE
          READ (5,FMT=9280) STEPL,DELTAMAX,T_DEBYE,DIMASS(1)
          READ (5,FMT=9310) EPS_MD
          READ (5,FMT=9120) (KATDYN(I),I=1,NATPER)
          DO I = 2,NATOMD
              DIMASS(I) = DIMASS(1)
          END DO
          WRITE (6,FMT=*) 'KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD'
          WRITE (6,FMT=9120) KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD
          WRITE (6,FMT=*) 'K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST'
          WRITE (6,FMT=9120) K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST
          WRITE (6,FMT=*) 'IOBROY,IOBROG,IOFILE'
          WRITE (6,FMT=9120) IOBROY,IOBROG,IOFILE
          WRITE (6,FMT=*) 'STEPL,DELTAMAX,T_DEBYE,DIMASS(1)'
          WRITE (6,FMT=9280) STEPL,DELTAMAX,T_DEBYE,DIMASS(1)
          WRITE (6,FMT=*) 'EPS_MD'
          WRITE (6,FMT=9310) EPS_MD
      END IF
c
c     end kmold
c
      IF (INS.GT.0) READ (5,FMT=9120) (IRNS(I),I=1,NATREF+NATPER)
      READ (5,FMT=9120) ITCLST,IMIX,IEF
      READ (5,FMT=9090) STRMIX,FCM,QBOUND
      READ (5,FMT=9090) BRYMIX,EF(1),EF(2)
      READ (5,FMT=9090) HFIELD,VBC(1),VCONSTH
      READ (5,FMT=9120) KLATR,KSYMM,LKONV,KSYMMAT,KESYM,KOCCUP
CBZ-INTEGRATION CZ
      IF (KLATR.EQ.10) THEN
          WRITE (*,FMT=*) 'GREEN FUNCTION FOR RELAXATION DONE BY'
          WRITE (*,FMT=*) 'BZ-INTEGRATION --- GTRANS HAS TO BE READ IN!'
          READ (5,FMT=9120) (OCCUP(I),I=1,NATPER)
      END IF
c
c Print the input from the job card
c
      WRITE (6,FMT=*) 'IRESIST,IFILE,IPE,IWTMAT'
      WRITE (6,FMT=9120) IRESIST,IFILE,IPE,IWTMAT
      WRITE (6,FMT=*) 'NSPIN,IRM,  INS,ICST'
      WRITE (6,FMT=9120) NSPIN,IRM,INS,ICST
      WRITE (6,FMT=*) 'KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC'
      WRITE (6,FMT=9120) KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC
      WRITE (6,FMT=*) 'KTE,  KPRE,KSPH,KEFG,KVMAD,KF ,IGGA'
      WRITE (6,FMT=9120) KTE,KPRE,KSPH,KEFG,KVMAD,KF,IGGA
      WRITE (6,FMT=*) 'NATREF,NATPER,ICUT,KSHAPE'
      WRITE (6,FMT=9120) NATREF,NATPER,ICUT,KSHAPE
      WRITE (6,FMT=*) '(IREF(I),I=1,NATPER)'
      WRITE (6,FMT=9120) (IREF(I),I=1,NATPER)
      WRITE (6,FMT=*) '(NTCELL(I),I=1,NATREF+NATPER)'
      WRITE (6,FMT=9120) (NTCELL(I),I=1,NATREF+NATPER)
      IF (INS.GT.0) WRITE (6,FMT=*) '(IRNS(I),I=1,NATREF+NATPER)'
      IF (INS.GT.0) WRITE (6,FMT=9120) (IRNS(I),I=1,NATREF+NATPER)
      WRITE (6,FMT=*) 'ITCLST,IMIX,IEF'
      WRITE (6,FMT=9120) ITCLST,IMIX,IEF
      WRITE (6,FMT=*) 'STRMIX,      FCM,     QBOUND'
      WRITE (6,FMT=9090) STRMIX,FCM,QBOUND
      WRITE (6,FMT=*) 'BRYMIX,     EF(1),    EF(2)'
      WRITE (6,FMT=9090) BRYMIX,EF(1),EF(2)
      WRITE (6,FMT=*) 'HFIELD,VBC(1),VCONSTH'
      WRITE (6,FMT=9090) HFIELD,VBC(1),VCONSTH
      WRITE (6,FMT=*) 'KLATR,KSYMM,LKONV,KSYMMAT,KESYM,KOCCUP'
      WRITE (6,FMT=9120) KLATR,KSYMM,LKONV,KSYMMAT,KESYM,KOCCUP
      IF (KLATR.EQ.10) THEN
          WRITE (6,FMT=*) 'OCCUPATION:'
          WRITE (6,FMT=9120) (OCCUP(I),I=1,NATPER)
      END IF
c
CASYM ----> This is read in, because NREP is used in next READ
      READ (35,FMT=9020) NREP,NATOM,LMAX,NATPER, (NSHELL(I),I=1,NATPER)

      DO I1 = 1,NSPIND
          DO I2 = 1,NREPD
              SOCCUP(I1,I2) = 0
              QOCCUPH(I1,I2) = 1
          END DO
      END DO

      IF (KSYMMAT.NE.0) THEN
          DO I2 = 1,NSPIN
              READ (5,FMT=9120) (ASARY(I1,I2),I1=1,NREP)
              DO I1 = 1,NREP
                  IF (KOCCUP.EQ.1) THEN
                      SOCCUP(I2,I1) = ASARY(I1,I2)
                  END IF
                  IF (ASARY(I1,I2).NE.0) ASARY(I1,I2) = I1
              END DO
              WRITE (6,FMT=*) 'ASARY'
              WRITE (6,FMT=9120) (ASARY(I1,I2),I1=1,NREP)
              WRITE (6,FMT=*) 'SOCCUP'
              WRITE (6,FMT=9120) (SOCCUP(I2,I1),I1=1,NREP)
          END DO
      END IF

      IF (KOCCUP.EQ.1 .AND. KSYMMAT.NE.0) THEN
          WRITE (6,FMT=*) 'QOCCUPOLD(NSPIN,NREPD)'
          DO I2 = 1,NSPIN
              READ (5,FMT=9120) (QOCCUPH(I2,I1),I1=1,NREP)
              WRITE (6,FMT=9120) (QOCCUPH(I2,I1),I1=1,NREP)
          END DO
      END IF

      IF (LKONV.EQ.0) LKONV = LMAXD
      IF (LKONV.NE.LMAXD) WRITE (6,FMT=*) 'LKONV ',LKONV,' IS USED'
c
      WRITE (6,FMT=*)
      WRITE (6,FMT=*)
c
c
      VBC(2) = VBC(1)
      IF (IGGA.EQ.1) WRITE (6,FMT=*) ' GGA CALCULATION WITH PW91'
      IF (IGGA.EQ.2) WRITE (6,FMT=*) ' GGA CALCULATION WITH PW86'
c
      IF (KSHAPE.NE.0) KWS = 2
      IF (KSHAPE.LE.2) LSMEAR = 0
      IF (KSHAPE.EQ.3) LSMEAR = 1
      IF (KSHAPE.EQ.4) LSMEAR = 2
      IF (KSHAPE.EQ.5) LSMEAR = 3
      IF (KSHAPE.GE.6) LSMEAR = 4
      IF (KSHAPE.GE.3) KSHAPE = 2
c     lsmear=0: host+clu with kinks, smepot is not read in
c     lsmear=1: host with kinks, clu eq.mesh, smepot is not read in
c     lsmear=2: host+clu eq.mesh, smepot is read in for host
c     lsmear=3: host with kinks, clu eq.mesh, smepot is read in
c     lsmear=4: host+clu eq.mesh, smepot is read in
c
      NPTREF = NSPIN*NATREF
      NPTPER = NSPIN*NATPER
      NSPPOT = NPTREF + NPTPER
      NATYP = NATREF + NATPER
      NATPS = NATREF + 1
      NPTPS = NPTREF + 1
c
      IPF = 6
      IPFE = IPF + 3
      NKWS = KWS + 1
      NKXC = KXC + 1
c
CASYM READ (35,FMT=9020) NREP,NATOM,LMAX,NATPER, (NSHELL(I),I=1,NATPER)
      WRITE (6,FMT=*) ' NREP,NATOM,LMAX,NATPER'
      WRITE (6,FMT=9020) NREP,NATOM,LMAX,NATPER
      WRITE (6,FMT=*) ' NSHELL(I),I=1,NATPER'
      WRITE (6,FMT=9020) (NSHELL(I),I=1,NATPER)
      READ (35,FMT=9010) ((RM(I,M),I=1,3),NHOST(M),NATM(M),NIMP(M),M=1,
     +  NATOM)
      READ (35,FMT=9040) (NDIM(NP),NP=1,NREP)
      READ (35,FMT=9040) (NDG(NP),NP=1,NREP)
      WRITE (6,FMT=9050)
      WRITE (6,FMT=9060) (NATM(M), (RM(I,M),I=1,3),NHOST(M),NIMP(M),M=1,
     +  NATOM)
C    In case of lattice relaxation read new positions - old positions
C    of the atoms (new positions first)
C ...
      IF (KLATR.GT.0) THEN
          DO N = 1,NATOM
              READ (5,FMT=*) M,(RM2(I,N),I=1,3), (RM1(I,N),I=1,3)
CT          READ(5,FMT=9280) (RM2(I,N),I=1,3),(RM1(I,N),I=1,3)
c         WRITE(6,FMT=9280) (RM2(I,N),I=1,3),(RM1(I,N),I=1,3)
          END DO
          WRITE (6,FMT=9290)
          WRITE (6,FMT=9060) (NATM(M), (RM2(I,M),I=1,3),NHOST(M),
     +      NIMP(M),M=1,NATOM)
          DO N = 1,NATOM
              DO I = 1,3
                  IF (ABS(RM(I,N)-RM1(I,N)).GT.
     +                1.D-6) STOP ' Inconsistent lattice '
                  RM(I,N) = RM2(I,N)
                  DRM0(I,N) = RM(I,N) - RM1(I,N)
              END DO
          END DO
          N1 = 0
          DO N = 1,NATPER
              DO I = 1,3
                  DRM(I,N) = DRM0(I,N1+1)
              END DO
              N1 = N1 + NSHELL(N)
          END DO
      END IF
      DO 20 NP = 1,NREP
          NDGH(NP) = NDG(NP)
          NDG(NP) = NDG(NP)/2
          WRITE (6,FMT=9030) NP,NDIM(NP),NDIM(NP),NDG(NP)

          IF (NDIM(NP).GT.NSEC) STOP 17

   20 CONTINUE

      IF (KOCCUP.EQ.1) THEN
          WRITE (6,FMT=*) 'QOCCUPNEW=QOCCUPOLD*NSPIN/NDG'
          DO I1 = 1,NSPIN
              DO I2 = 1,NREPD
                  QOCCUPH2(I1,I2)=QOCCUPH(I1,I2)
                  QOCCUP(I1,I2) = NSPIN*QOCCUPH2(I1,I2)/NDGH(I2)
                  IF (QOCCUP(I1,I2).GT.1.0) THEN
                      WRITE (6,FMT=*) 'QOCCUP in reainp is wrong!'
                  END IF
              END DO
              WRITE (6,FMT=9280) (QOCCUP(I1,I2),I2=1,NREPD)          
          END DO
      END IF




      WRITE (6,FMT=9210) VBC(1)
      WRITE (6,FMT=9270) EF(1),EF(2)
c
c
      IF (QBOUND.LT.1.D-15) QBOUND = 1.D-4
      IF (IMIX.EQ.3 .OR. IMIX.EQ.4) THEN
          WRITE (6,FMT=9170) (IMIX-2)
          IF (FCM.NE.1.0D0) WRITE (6,FMT=9190) FCM
      END IF

      IF (IMIX.EQ.5) WRITE (6,FMT=9180)
      WRITE (6,FMT=9150) QBOUND
      WRITE (6,FMT=9110) STRMIX,FCM
      WRITE (6,FMT=9200) BRYMIX
c
      LMAXP1 = LMAX + 1
      LPOT = MIN(2*LMAX,LPOTD)
      LMXTW1 = LMAX + LMAX + 1
      LMAXSQ = LMAXP1*LMAXP1
      LMPOT = (LPOT+1)* (LPOT+1)
c
      WRITE (6,FMT=9000) LMAXP1,LMX,NATREF,NTREFD,NATPER,NTPERD,NREP,
     +  NREPD,IRM,IRMD,NSPIN,NSPIND,NATOM,NATOMD

      IF (LMAXP1.GT.LMX .OR. NATREF.GT.NTREFD .OR. NATPER.GT.NTPERD .OR.
     +    NREP.GT.NREPD .OR. IRM.GT.IRMD .OR. NSPIN.GT.NSPIND .OR.
     +    NATOM.GT.NATOMD) STOP 18

      IF (INS.GT.0) THEN
          WRITE (6,FMT=9230)
          DO 30 I = 1,NATYP
              WRITE (6,FMT=9240) I,IRNS(I),IRNSKD

              IF (IRNS(I).GT.IRNSKD) STOP

   30     CONTINUE
          IF (LMAXP1.NE.LMX) THEN
              WRITE (6,FMT=9220)

              STOP

          END IF

      END IF

      IF (KHFELD.EQ.1) WRITE (6,FMT=9070) HFIELD
      WRITE (6,FMT=9100) TSPIN(NSPIN)
      WRITE (6,FMT=9260) TVREL(KVREL)
      WRITE (6,FMT=9260) TKCOR(KCOR)
      IF (KSHAPE.EQ.0) THEN
          WRITE (6,FMT=9130) TKWS(NKWS)

      ELSE
          WRITE (6,FMT=9260) TSHAPE
      END IF
      IF (IGGA.EQ.0) THEN
          WRITE (6,FMT=9160) TXC(KXC+1)
      ELSE
          WRITE (6,FMT=9160) TXC(IGGA+3)
      END IF
      IF (INS.GT.0) WRITE (6,FMT=9250) TINS(INS),ICST
      WRITE (6,FMT=9140)
c
      RETURN

 9000 FORMAT (/,13x,'check of dimension-data consistency',/,13x,
     +       35 ('-'),/,20x,'lmax+1 : (',i6,',',i6,')',/,20x,
     +       'natref : (',i6,',',i6,')',/,20x,'natper : (',i6,',',i6,
     +       ')',/,20x,'nrep   : (',i6,',',i6,')',/,20x,'irm    : (',i6,
     +       ',',i6,')',/,20x,'nspin  : (',i6,',',i6,')',/,20x,
     +       'natom  : (',i6,',',i6,')',/)
 9010 FORMAT (3f10.6,i3,i4,i5)
 9020 FORMAT (11i5)
 9030 FORMAT (' for represention no.',i3,' matrix dimensions are:',i3,
     +       2H x,i3,'  degeneracy: ',i2)
 9040 FORMAT (16i5)
 9050 FORMAT (' atomic coordinates and types')
 9060 FORMAT (' atom:',i3,2x,3f10.6,' host type: ',i2,
     +       ' impurity type: ',i2)
 9070 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9080 FORMAT (3x,4i1)
 9090 FORMAT (3f10.6)
 9100 FORMAT (1x,a4,'spin polarized calculation')
 9110 FORMAT (' mixing factor used:',f15.9)
 9120 FORMAT (8i5)
 9121 FORMAT (8f5.0)
 9130 FORMAT (1x,' calculation with',a8,'-potential')
 9140 FORMAT (1x,79 ('*'))
 9150 FORMAT (' convergence quality required:  ',1p,d12.2)
 9160 FORMAT (1x,a24,'exchange-correlation potential')
 9170 FORMAT (' broyden"s method # :',i1,' used')
 9180 FORMAT (' generalized Anderson scheme used')
 9190 FORMAT (/,' spin mixing should not be enhanced ==> fcm=1.0 ',/,1x,
     +       'but actual fcm : ',f5.2)
 9200 FORMAT (' optional parameter for broyden-update :',f15.9)
 9210 FORMAT (1x,'constant shift of the potentials :',f10.6)
 9220 FORMAT (1x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmx has to be set equal lmax + 1 ! ')
 9230 FORMAT (' full potential calulation ',
     +       '- cut off of non spherical potential')
 9240 FORMAT (' representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9250 FORMAT (1x,a43,/,1x,' using',i3,'-th. born approximation ')
 9260 FORMAT (1x,a43)
 9270 FORMAT (1X,'LDOS FROM EF  EF(1)=',F10.6,'EF(2)=',F10.6)
 9280 FORMAT (6F10.5)
 9290 FORMAT (' The lattice is relaxed, new atomic positions ',/)
 9300 FORMAT (i5)
 9310 FORMAT (4E10.5)
      END
