      SUBROUTINE ETOTALsc(DNEINT,DSELOC,ECOU,EFERMI,EGR,EPOTIN,ESPC,
     $     ESPV,EXC,IREF,KPRE,KFSR,LMAX,LPOT,NATPER,NATREF,
     +     NSPIN,NSHELL,dslsc1,dslsc2)
      Implicit None
c-----------------------------------------------------------------------
c     calculate the total energy of the cluster .
c     gather all energy-parts which are calculated in different
c     subroutines .
c     since the program uses group theory only shell-indices
c     are used instead of atom-indices .
c
c                               b.drittler   may 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD
      PARAMETER (NATYPD=38,NTREFD=1)
      INTEGER LMAXD,LMX,LPOTD
      PARAMETER (lmaxd=4,LMX=LMAXD+1,lpotd=8)
      INTEGER LD
      PARAMETER (LD=LMX-1)
C     ..
C     .. Scalar Arguments ..
      REAL*8 DNEINT,DSELOC,EFERMI, dslsc1, dslsc2
      INTEGER KFSR,KPRE,LMAX,LPOT,NATPER,NATREF,NSPIN
C     ..
C     .. Array Arguments ..
      REAL*8 ECOU(0:LPOTD,*),EGR(*),EPOTIN(*),
     +                 ESPC(0:LD,NATYPD,*),ESPV(0:LD,NATYPD,*),
     +                 EXC(0:LPOTD,*)
      INTEGER IREF(*),NSHELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 EDIFF,ET,ETOT
      INTEGER I1,IATYP,IPRE,IS,ISPIN,L,NREF
C     ..
C     .. Local Arrays ..
      REAL*8 EREF(NTREFD)
      CHARACTER*4 TEXTL(0:5)
      CHARACTER*13 TEXTS(3)
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h ='/
      DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
C     ..
c
      IF (KPRE.EQ.1) THEN
        IPRE = 6

      ELSE

        IPRE = 56
      END IF

      ETOT = 0.0D0
      EDIFF = 0.0D0

      WRITE (6,FMT=9000)
      WRITE (6,FMT=9010) EFERMI
      WRITE (IPRE,FMT=9210)
      WRITE (IPRE,FMT=9170)
      WRITE (IPRE,FMT=9210)
c
c---> loop over host atoms
c
      DO 40 IATYP = 1,NATREF
        WRITE (IPRE,FMT=9020) IATYP
        ET = 0.0D0
        IS = 0
        IF (NSPIN.EQ.1) IS = IS + 2
        DO 20 ISPIN = 1,NSPIN
          IS = IS + 1
          WRITE (IPRE,FMT=9040) TEXTS(IS)
          WRITE (IPRE,FMT=9050) (TEXTL(L),ESPC(L,IATYP,ISPIN),L=0,LMAX)
          WRITE (IPRE,FMT=9060) (TEXTL(L),ESPV(L,IATYP,ISPIN),L=0,LMAX)
          DO 10 L = 0,LMAX
            ET = ET + ESPC(L,IATYP,ISPIN) + ESPV(L,IATYP,ISPIN)
   10     CONTINUE
   20   CONTINUE
c
        WRITE (IPRE,FMT=9070) ET
        WRITE (IPRE,FMT=9080) (L,ECOU(L,IATYP),L=0,LPOT)
        WRITE (IPRE,FMT=9090) (L,EXC(L,IATYP),L=0,LPOT)
        DO 30 L = 0,LPOT
          ET = ET + ECOU(L,IATYP) + EXC(L,IATYP)
   30   CONTINUE
        ET = ET + EPOTIN(IATYP) + EGR(IATYP)
        WRITE (IPRE,FMT=9100) ET
        EREF(IATYP) = ET
   40 CONTINUE
c
      WRITE (IPRE,FMT=9210)
      WRITE (IPRE,FMT=9180)
      WRITE (IPRE,FMT=9210)
c
c---> loop over shells (first index)
c
      DO 80 I1 = 1,NATPER
        NREF = IREF(I1)
        IATYP = I1 + NATREF
        WRITE (IPRE,FMT=9030) I1
        ET = 0.0D0
        IS = 0
        IF (NSPIN.EQ.1) IS = IS + 2
        DO 60 ISPIN = 1,NSPIN
          IS = IS + 1
          WRITE (IPRE,FMT=9040) TEXTS(IS)
          WRITE (IPRE,FMT=9050) (TEXTL(L),ESPC(L,IATYP,ISPIN),L=0,LMAX)
          WRITE (IPRE,FMT=9060) (TEXTL(L),ESPV(L,IATYP,ISPIN),L=0,LMAX)
          DO 50 L = 0,LMAX
            ET = ET + ESPC(L,IATYP,ISPIN) + ESPV(L,IATYP,ISPIN)
   50     CONTINUE
   60   CONTINUE
c
        WRITE (IPRE,FMT=9070) ET
        WRITE (IPRE,FMT=9080) (L,ECOU(L,IATYP),L=0,LPOT)
        WRITE (IPRE,FMT=9090) (L,EXC(L,IATYP),L=0,LPOT)
        DO 70 L = 0,LPOT
          ET = ET + ECOU(L,IATYP) + EXC(L,IATYP)
   70   CONTINUE
        ET = ET + EPOTIN(IATYP) + EGR(IATYP)
        WRITE (IPRE,FMT=9110) ET
        WRITE (IPRE,FMT=9160) ET - EREF(NREF),
     +    NSHELL(I1)* (ET-EREF(NREF))
        ETOT = ETOT + NSHELL(I1)*ET
        EDIFF = EDIFF + NSHELL(I1)* (ET-EREF(NREF))
   80 CONTINUE
      WRITE (IPRE,FMT=9210)
      WRITE (6,FMT=9120) ETOT
      WRITE (6,FMT=9220) EDIFF
      IF (KFSR.EQ.1) THEN
        WRITE (IPRE,FMT=9210)
        WRITE (6,FMT=9190) DSELOC
        WRITE (6,FMT=9200) DNEINT
        WRITE (6,FMT=9130) DNEINT - DSELOC
        Write (6,fmt=9300) dslsc1,dslsc2
        Write (6,fmt=9310) dslsc2 - dslsc1
        WRITE (IPRE,FMT=9210)
        WRITE (6,FMT=9140) ETOT + DNEINT - DSELOC
        WRITE (6,FMT=9230) EDIFF + DNEINT - DSELOC
      END IF

      WRITE (6,FMT=9150)



 9000 FORMAT (/,1x,33 ('-'),' total energies ',33 ('-'),/)
 9010 FORMAT (/,13x,'fermi energy respectively the electrostatic zero ',
     +       'in ryd : ',f10.6)
 9020 FORMAT (3x,'total energy of the ',i3,'-th. host atom :')
 9030 FORMAT (3x,'contribution of the ',i3,'-th. shell :')
 9040 FORMAT (5x,'single particle energies ',a13)
 9050 FORMAT (7x,'  core   contribution : ',6 (a4,f15.8))
 9060 FORMAT (7x,'valence  contribution : ',6 (a4,f15.8))
 9070 FORMAT (7x,'total contribution of the single particle energies :',
     +       1x,f15.8)
 9080 FORMAT (7x,'coulomb  contribution : ',4 (i3,1x,f15.8),/,
     +       4 (i3,1x,f15.8))
 9090 FORMAT (7x,'ex.-cor. contribution : ',4 (i3,1x,f15.8),/,
     +       4 (i3,1x,f15.8))
 9100 FORMAT (/,3x,'total contribution of the host atom : ',f15.8,/)
 9110 FORMAT (/,3x,'total contribution of the representive atom : ',
     +       f15.8,/)
 9120 FORMAT (1x,/,3x,'total energy in ryd. : ',f17.8)
 9130 FORMAT (7x,'correction of the valence single particle energies',
     +       ' with lloyd"s formula : ',f17.8)
 9140 FORMAT (/,3x,'total energy calculated with lloyd"s formula ',
     +       'in ryd. : ',f17.8)
 9150 FORMAT (1x,82 ('-'))
 9160 FORMAT (1x,/,5x,
     +       'energy difference to host atom per atom  in ryd :',f15.8,
     +       /,5x,'energy difference to host atom per shell in ryd :',
     +       f15.8)
 9170 FORMAT (1x,/,1x,' total energy of the host atoms ',/,1x)
 9180 FORMAT (1x,/,1x,' total energy of the disturbed atoms ',/,1x)
 9190 FORMAT (1x,/,9x,' change of single particle energies ',
     +       '(local integrated) :',f15.8)
 9200 FORMAT (9x,' change of single particle energies ',
     +       '(lloyd"s formula)  :',f15.8)
 9210 FORMAT (1x,'>')
 9220 FORMAT (1x,3x,'total energy difference in ryd. : ',f17.8)
 9230 FORMAT (3x,'energy difference calculated with lloyd"s formula ',
     +       'in ryd. : ',f17.8)
 9300 FORMAT (9x,'change of sc s.p.e. (old, scaled) : ',2f15.8)
 9310 FORMAT (9x,'correction to sc s.p.e by scscaling: ',f17.8)
      END
