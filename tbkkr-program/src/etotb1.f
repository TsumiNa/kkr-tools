c 13.04.2000**************************************************************
      SUBROUTINE ETOTB1(ECOU,EFERMI,EPOTIN,ESPC,ESPV,EXC,KPRE,LMAX,LPOT,
     +                 NSPIN,NATYP,NSHELL,espcor,espco2,npntsc)
      implicit none
c ************************************************************************
c     calculate the total energy of the cluster .
c     gather all energy-parts which are calculated in different
c     subroutines .
c     since the program uses group theory only shell-indices
c     are used instead of atom-indices .
c
c                               b.drittler   may 1987
c
c     modified for supercells with nshell(i) atoms of type i in the
c     unit cell
c                               p.zahn       oct. 95
c
c     takes semicore correction into account
c                               h.hoehler    april 2000
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
C      INTEGER NATYPD
C      PARAMETER (NATYPD=1)
C      INTEGER LMAXD,LPOTD
C      PARAMETER (LMAXD=4,LPOTD=8)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EFERMI,espcor,espco2
      INTEGER KPRE,LMAX,LPOT,NATYP,NSPIN,npntsc
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ECOU(0:LPOTD,*),EPOTIN(*),ESPC(0:LMAXD,NATYPD,*),
     +                 ESPV(0:LMAXD,NATYPD,*),EXC(0:LPOTD,*)
      INTEGER NSHELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ECOUS,EXCS,EDC,ET,ETOT,EFCTOR
      INTEGER IATYP,IS,ISPIN,L
C     ..
C     .. Local Arrays ..
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*13 TEXTS(3)
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
C     ..
c     .. externals
      LOGICAL TEST
      EXTERNAL TEST
c
c ------------------------------------------------------------------------
      EFCTOR = 1.0D0/13.6058D0
C      IF (KPRE.EQ.1) THEN

      ETOT = 0.0D0

      IF (KPRE.EQ.1) THEN
        WRITE (6,FMT=9000)
        WRITE (6,FMT=9010) EFERMI
        WRITE (6,FMT=9110)
        WRITE (6,FMT=9110)

c     
c---> loop over host atoms
c     
      DO 40 IATYP = 1,NATYP

        WRITE (6,FMT=9020) IATYP
        
        EDC = 0.0D0
        ET  = 0.0D0
        ECOUS=0.0D0
        EXCS=0.0D0

        IS  = 0
        IF (NSPIN.EQ.1) IS = IS + 2

        DO 20 ISPIN = 1,NSPIN
          IS = IS + 1
          
          WRITE (6,FMT=9030) TEXTS(IS)
          WRITE (6,FMT=9040) (TEXTL(L),ESPC(L,IATYP,ISPIN),L=0,LMAX)
          WRITE (6,FMT=9050) (TEXTL(L),ESPV(L,IATYP,ISPIN),L=0,LMAX)

          
          DO 10 L = 0,LMAX
            ET = ET + ESPC(L,IATYP,ISPIN) + ESPV(L,IATYP,ISPIN)
 10       CONTINUE
 20     CONTINUE

        WRITE (6,FMT=9060) ET
        WRITE (6,FMT=9070) (L,ECOU(L,IATYP),L=0,LPOT)
        WRITE (6,FMT=9080) (L,EXC(L,IATYP),L=0,LPOT)

        DO 31 L = 0,LPOT
          ET = ET + ECOU(L,IATYP) + EXC(L,IATYP)
 31     CONTINUE
          ET = ET + EPOTIN(IATYP)
        IF (NATYP.GT.1) WRITE (6,FMT=9090) ET
          ETOT = ETOT + ET
c
c
c --->  sum up Coulomb and Ex.-Corel. contribution
c
        DO 30 L = 0,LPOT
          ECOUS  = ECOUS  + ECOU(L,IATYP)
          EXCS   = EXCS + EXC(L,IATYP)
 30     CONTINUE
C
  
        WRITE (6,FMT=9060) ET
        WRITE (6,FMT=9070) (L,ECOU(L,IATYP),L=0,LPOT)
        WRITE (6,FMT=9150) ECOUS
        WRITE (6,FMT=9080) (L,EXC(L,IATYP),L=0,LPOT)
        WRITE (6,FMT=9140) EXCS
        WRITE (6,FMT=9120) EPOTIN(IATYP)      
        
        IF (TEST('NoMadel ')) GOTO 90

        ET  = ET  + ECOUS + EXCS
        EDC = EDC + ECOUS + EXCS
        
        ET  = ET  + EPOTIN(IATYP)
        EDC = EDC + EPOTIN(IATYP)

        IF (KPRE.EQ.1) WRITE (6,FMT=9130) EDC

 90     IF (NATYP.GT.1 .OR. NSHELL(IATYP).GT.1) 
     +       WRITE (6,FMT=9090) ET

        ETOT = ETOT + ET*DBLE(NSHELL(IATYP))

 40   CONTINUE                      ! IATYP = 1,NATYP

      WRITE (6,FMT=9100) ETOT,ETOT/EFCTOR

        If ( npntsc .gt. 0) Then
          Write (6,fmt=9121) etot+espcor+espco2
          Write (6,FMT='(3x,a,2f15.8)')
     $         'semicore corrections to energy : ',
     $         espcor,espco2
        End If

      ELSE

        ETOT = 0.0D0
        WRITE (6,FMT=9000)
        WRITE (6,FMT=9010) EFERMI
        WRITE (6,FMT=9110)
        WRITE (6,FMT=9110)
c
c---> loop over host atoms
c
        DO 80 IATYP = 1,NATYP
          ET = 0.0D0
          IS = 0
          IF (NSPIN.EQ.1) IS = IS + 2
          DO 60 ISPIN = 1,NSPIN
            IS = IS + 1

          WRITE (6,FMT=9030) TEXTS(IS)
          WRITE (6,FMT=9040) (TEXTL(L),ESPC(L,IATYP,ISPIN),L=0,LMAX)
          WRITE (6,FMT=9050) (TEXTL(L),ESPV(L,IATYP,ISPIN),L=0,LMAX)

            DO 50 L = 0,LMAX
              ET = ET + ESPC(L,IATYP,ISPIN) + ESPV(L,IATYP,ISPIN)
   50       CONTINUE
   60     CONTINUE

          WRITE (6,FMT=9060) ET,ETOT/EFCTOR
          WRITE (6,FMT=9070) (L,ECOU(L,IATYP),L=0,LPOT)
          WRITE (6,FMT=9080) (L,EXC(L,IATYP),L=0,LPOT)

          DO 70 L = 0,LPOT
            ET = ET + ECOU(L,IATYP) + EXC(L,IATYP)
   70     CONTINUE

          ET = ET + EPOTIN(IATYP)
          IF (NATYP.GT.1) WRITE (6,FMT=9090) ET
          ETOT = ETOT + ET

   80   CONTINUE
        WRITE (6,FMT=9100) ETOT,ETOT/EFCTOR
        If ( npntsc .gt. 0 ) Then
          Write (6,fmt=9121) etot+espcor+espco2
          Write (6,FMT='(3x,a,2f15.8)')
     $         'semicore corrections to energy : ',
     $         espcor,espco2
        End If
      END IF

c ------------------------------------------------------------------------
C      ELSE                          ! (KPRE.EQ.1)
C
C        ETOT = 0.0D0
c
c---> loop over host atoms
c
C        DO 80 IA = 1,NAEZ
C          INV   = EQINV(IA)
C          IATYP = KAOEZ(INV)
C          ET = 0.0D0
C          IS = 0
C          IF (NSPIN.EQ.1) IS = IS + 2
C          DO 60 ISPIN = 1,NSPIN
C            IS = IS + 1
C            DO 50 L = 0,LMAX
C              ET = ET + ESPC(L,IATYP,ISPIN) + ESPV(L,IATYP,ISPIN)
C   50       CONTINUE
C   60     CONTINUE
C          DO 70 L = 0,LPOT
C            ET = ET + ECOU(L,IATYP) + EXC(L,IATYP)
C   70     CONTINUE
C          ET = ET + EPOTIN(IATYP)
C          IF (NATYP.GT.1) WRITE (6,FMT=9090) ET
C          ETOT = ETOT + ET
C   80   CONTINUE
C        WRITE (6,FMT=9100) ETOT,ETOT/EFCTOR
C
C      END IF                        ! (KPRE.EQ.1)
c ------------------------------------------------------------------------

      RETURN

 9000 FORMAT (/,1x,30 ('-'),' total energies ',30 ('-'))
 9010 FORMAT (/,13x,'fermi energy in ryd : ',f10.6)
 9020 FORMAT (3x,'total energy of the ',i3,'-th. atom :')
 9030 FORMAT (5x,'single particle energies ',a13)
 9040 FORMAT (7x,'  core   contribution : ',2(a4,f15.8),/,
     +     (31x,2(a4,f15.8)))
 9050 FORMAT (7x,'valence  contribution : ',2(a4,f15.8),/,
     +     (31x,2(a4,f15.8)))
 9060 FORMAT (7x,
     +       'total contribution of the single particle energies :',
     +       1x,f15.8)
 9070 FORMAT (7x,'coulomb  contribution : ',2(i3,1x,f15.8),/,
     +     (31x,2(i3,1x,f15.8)))
 9080 FORMAT (7x,'ex.-cor. contribution : ',2(i3,1x,f15.8),/,
     +     (31x,2(i3,1x,f15.8)))
 9090 FORMAT (/,3x,'total contribution of the atom : ',f15.8,/)
 9100 FORMAT (1x,/,3x,'total energy in ryd. : ',f17.8,/
     +             3x,'                 eV  : ',F17.8)
 9110 FORMAT (1x,'>')
 9120 FORMAT (7x,'eff. pot. contribution     : ',f15.8)
 9121 FORMAT (3x,'total energy in ryd. : ',f17.8,1x,
     $       '(with semicore correction)',/)
 9130 FORMAT (7x,
     +       'total double counting contribution                 :',
     +       1x,f15.8)
 9140 FORMAT (7x,'tot. ex.-cor. contribution : ',f15.8)
 9150 FORMAT (7x,'tot. coublomb contribution : ',f15.8)

      END                           ! ETOTB1
