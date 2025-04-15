C ************************************************************************
      SUBROUTINE EIGEN1(HEAD,ITITLE,M2,MMIN,MMAX,QBOUND,
     +     EB,E1,E2,TK,NPNT2,RBASIS,
     +     IOPSYS,IOPREF,IPF,IPFE,IPE)
      implicit none
C ************************************************************************
c
c     this is the driver for the eigen value search
c
c     to save storage the wavefunctions are not stored for
c     all atoms in the unit cell
c
c      p. zahn, febr. 96
C ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
c      include 'inc.lay'

      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LM2D
      PARAMETER (LM2D = (2*LMAXD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER NREFPOTD
      PARAMETER (NREFPOTD=NSPIND*NREFD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER ALM,NDIM
      PARAMETER (ALM = NDIMGK*LMAXSQ,NDIM = NPRINCD*LMAXSQ)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
      INTEGER NAUX,NAUXSYM
      PARAMETER (NAUX=2*ALM**2+5*ALM,NAUXSYM=ALM)
      DOUBLE COMPLEX CONE,CZERO,CONEM,CI
      PARAMETER (
     +     CONE  = ( 1.0D0,0.0D0),
     +     CZERO = ( 0.0D0,0.0D0),
     +     CONEM = (-1.0D0,0.0D0),
     +     CI    = ( 0.0D0,1.0D0))
      DOUBLE PRECISION PI,TPI,ZERO,ONE
      PARAMETER (
     +     PI    = 3.14159265358979312D0,
     +     TPI   = 2.0d0 * PI,
     +     ZERO  = 0.0D0,
     +     ONE   = 1.0D0)
C     ..
C     .. Scalar Arguments ..
      INTEGER 
     +     IPF,IPFE,IPE,
     +     IOPSYS,IOPREF,           ! files for input of potentials
     +     M2,MMIN,MMAX,            ! maximum number of bands
     +     NPNT2                    ! energy intervalls in (E1,E2)

      DOUBLE PRECISION EB,E1,E2,TK,QBOUND
C     ..
C     .. Array Arguments ..
      INTEGER ITITLE(20,*)
      DOUBLE PRECISION RBASIS(3,*)  ! position of atoms in the unit cell
                                    ! in units of bravais vectors
      CHARACTER*1          HEAD(80)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  
     +     ALATNEW,AKSCALE,BKSCALE,CKSCALE,
     +     CONVPU,CSQR,
     +     D,DE,DOS,
     +     EA,ED,EDIFF,EE,EFCTOR,EFERMI,EIMAG,EM,ETK,E2IN,
     +     FAC,FACI,
     +     KB,
     +     POLBOUND,
     +     TIMEE,TIMEI,TIMES,
     +     W1,W2,DCLOCK,TIMET
      DOUBLE COMPLEX 
     +     E
      INTEGER 
     +     I,IA,I1,IC,IE,II,IL,ILM1,IMIN,IMAX,INV,INFO,ITOT,
     +     IOPT,IP,IRF,ISPIN,IVAL1,IVAL2,
     +     J,JLM1,JLM2,JRF,
     +     LAYMAX,
     +     LM,LM1,LM2,
     +     N,N1,NPNT,NPOINTS,NPOLES,NSHELL,NUE,NUEL,
     +     NUM0,NUM01,NUM02,NUM1,NUM2,
     +     RF,IO,N2
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION 
     +     A1(NREFD),
     +     B1(NREFD),               ! constants for exponential r mesh
     +     ECORE(20,NPOTD),         ! core states
     +     ECORE1(20,NREFPOTD),     ! core states
     +     RMTNEW(NATYPD),RMTNEW1(NREFD), ! adapted MT radius
     +     RWS1(NREFD),             ! Wigner Seitz radius
     +     RMT1(NREFD),             ! Muffin-Tin-radius
     +     THETAS1(IRID,NFUND,NREFD),
     +     Q(6),QQ(6),
     +     QKP(6,KPOIBZ),

     +     VF(3),
     +     EQ(0:NAEZD*LMMAXD)

      DOUBLE COMPLEX 
     +     EZ(NPNTD)

      INTEGER 
     +     IFUNM(NATYPD,LMXSPD),
     +     IPVT(LMAXSQ),
     +     IRC(NATYPD),
     +     LCORE(20,NPOTD),         ! angular momentum of core states
     +     LLMSP(NATYPD,NFUND),     ! lm=(l,m) of 'nfund'th nonvanishing
                                    ! component of non-spherical pot.
     +     LMSP(NATYPD,LMXSPD),     ! 0,1 : non/-vanishing lm=(l,m) component
                                    ! of non-spherical potential
     +     NCORE(NPOTD),            ! number of core states
     +     NFU(NATYPD),
c
c ---> for reference system
c
     +     IFUNM1(NREFD,LMXSPD),
     +     IMT1(NREFD),
     +     IRC1(NREFD),
c     +     IRNS1(NREFD),
     +     LCORE1(20,NREFPOTD),         ! angular momentum of core states
     +     LLMSP1(NREFD,NFUND),     ! lm=(l,m) of 'nfund'th nonvanishing
                                    ! component of non-spherical pot.
     +     LMSP1(NREFD,LMXSPD),     ! 0,1 : non/-vanishing lm=(l,m) component
                                    ! of non-spherical potential
     +     NCORE1(NREFPOTD),            ! number of core states
     +     NFU1(NREFD),

     +     ISTEP(NAEZD*LMMAXD)      ! number of iteration steps for 
                                    ! energy eigenvalues

      LOGICAL 
     +     FOUND,LSTART,LSYM,LTEST, 
     +     TEST,OPT
      LOGICAL SELECT(ALM)

c     .. external statement
      EXTERNAL 
     +     EINIT,EREST,EREST1,ESAVE,
     +     EIGWFCT,EIGISO,
     +     OPT,RINIT,TEST,
     +     DELTAMAT,MINEIG,NPOLES,STARTB1
      INTRINSIC DABS,DIMAG,DMIN1,DREAL,DSIGN,LOG,MAX,ABS
c     .. data
      INTEGER LF(144)
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
     +     KB       / 0.6333659D-5/
c
c     .. arrays in common
      INTEGER     IUP,IDO,NUM
      DOUBLE PRECISION WUP,WDO
      COMMON /EIGENV/ WUP,WDO,NUM,IUP,IDO

      INTEGER NDIFF(NSPIND,20),NZERO(NSPIND)
      DOUBLE PRECISION EZERO(NSPIND,20),EBOT
      COMMON /EZERO / EBOT,EZERO,NZERO,NDIFF

      DOUBLE PRECISION EE1(0:M2D,KPOIBZ),EE2(0:M2D,KPOIBZ),
     +                 WW1(0:M2D,KPOIBZ),WW2(0:M2D,KPOIBZ)
      INTEGER II1(0:M2D,KPOIBZ),II2(0:M2D,KPOIBZ)
      COMMON /EE/ EE1,EE2,WW1,WW2,II1,II2 ! index arrays for eigenvalue
                                          ! analysis

      INTEGER NBASIS(KPOIBZ),NMIN(KPOIBZ),NMAX(KPOIBZ)
      COMMON /NBASIS / NBASIS,NMIN,NMAX

      DOUBLE COMPLEX CL(NAEZD*LMMAXD)
      COMMON /CL/ CL

c ------------------------------------------------------------------------
      INTEGER
     +     IPAN(NATYPD),
     +     IRMIN(NATYPD),
     +     IRCUT(0:IPAND,NATYPD)
      INTEGER
     +     IRWS1(NREFD),            ! r point at WS radius
     +     IPAN1(NREFD),
     +     IRMIN1(NREFD),
     +     IRCUT1(0:IPAND,NREFD)

      DOUBLE PRECISION 
     +     R1(IRMD,NREFD),
     +     RS1(IRMD,0:LMAXD,NREFD),
     +     S1(0:LMAXD,NREFD),
     +     DRDI1(IRMD,NREFD),
     +     DROR1(IRMD,NREFD)
      DOUBLE PRECISION 
     +     Z1(NREFD),
     +     VINS1(IRMIND:IRMD,LMPOTD,NREFPOTD),
     +     VISP1(IRMD,NREFPOTD)        ! spherical input potential

      COMMON / REFPARA /
     +         R1,RS1,S1,DRDI1,DROR1,
     +         Z1,VINS1,VISP1,
     +         IPAN,IRMIN,IRCUT,
     +         IRWS1,IPAN1,IRMIN1,IRCUT1
c ------------------------------------------------------------------------
      DOUBLE PRECISION 
     +     ALATC,C,VBC,
     +     A(NATYPD),B(NATYPD),     ! contants for exponential r mesh
     +     R(IRMD,NATYPD),          ! r mesh
     +     RR(3,0:NRD),
     +     RS(IRMD,0:LMAXD,NATYPD),
     +     S(0:LMAXD,NATYPD),
     +     DRDI(IRMD,NATYPD),
     +     DROR(IRMD,NATYPD),
     +     THETAS(IRID,NFUND,NCELLD),
     +     VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     +     VISP(IRMD,NPOTD),         ! spherical input potential
     +     Z(NATYPD)
      DOUBLE COMPLEX 
     +     DELTALL(LMMAXD,LMMAXD,NAEZD),
     +     GINP(NACLSD*LMMAXD,LMMAXD,NCLSD),     ! cluster GF(ref syst.)
     +     TMATLL(LMMAXD,LMMAXD,NAEZD),
     +     TMATLL1(LMMAXD,LMMAXD,NATYPD),
     +     TREFLL(LMMAXD,LMMAXD,NREFD)
      INTEGER
     +     INS,INSREF,
     +     KVREL,
     +     KSHAPE,
     +     KWS,                     ! 0 (MT), 1(ASA)
     +     KMT,
     +     KHFELD,                  ! 0,1: no / yes external magnetic field
     +     KSCOEF
      INTEGER
     +     LOFLM(LM2D),
     +     ICLEB(NCLEB,4),
     +     IEND,
     +     LMAX,LMMAX,LPOT,
     +     NSPIN,NREF,NATYP,
     +     NAEZ,                    ! number of atoms in unit cell
     +     NINEQ,
     +     NSHELL1(0:NSHELD),
     +     NZ,                      ! number of atoms at centers of inversion
     +     NLAYER,
     +     IMT(NATYPD),
     +     IRWS(NATYPD),            ! r point at WS radius
     +     IRNS(NATYPD),
     +     NTCELLR(NREFD),
     +     NTCELL(NATYPD)
      DOUBLE PRECISION 
     +     HFIELD,
     +     VCONST,
     +     MTFAC(NATYPD),
     +     CLEB(NCLEB,2)
      DOUBLE PRECISION 
     +     RWS(NATYPD),                  ! Wigner Seitz radius
     +     RMT(NATYPD),              ! Muffin-Tin-radius
     +     RCLS(3,NACLSD,NCLSD)    ! real space position of atom in cluster
      INTEGER
     +     REFPOT(NATYPD+NEMBD),
     +     NCLS,                    ! number of cluster
     +     NACLS(NCLSD),
     +     CLS(NATYPD),
     +     EZOA(NACLSD,NAEZD),
     +     ATOM(NACLSD,NAEZD),
     +     EQINV(NAEZD),
     +     KAOEZ(NAEZD+NEMBD),
     +     INIPOL(NATYPD),          ! initial spin polarisation
     +     LATT,                   
     +     ICC,
     +     ICST                     ! Born approximation

      COMMON / SYSPARA /
     +         TMATLL,TMATLL1,TREFLL,DELTALL,GINP,
     +         ALATC,C,
     +         A,B,R, RR, RS, S, DRDI, DROR,
     +         THETAS,VINS,VISP,Z,
     +         HFIELD,VCONST,MTFAC,CLEB,
     +         RWS,RMT,RCLS,
     +         INS,INSREF,KSHAPE,KVREL,KWS,KMT,KSCOEF,KHFELD,
     +         LOFLM,ICLEB,IEND,
     +         LMAX,LMMAX,LPOT,NSPIN,
     +         NREF,NATYP,NAEZ,NINEQ,NSHELL1,NZ,NLAYER,
     +         IMT,IRWS,IRNS,NTCELL,NTCELLR,
     +         REFPOT,NCLS,NACLS,
     +         CLS,EZOA,ATOM,EQINV,KAOEZ,INIPOL,LATT,ICC,ICST
c ------------------------------------------------------------------------
      IF(TEST('flow    '))
     +     WRITE(6,*) '>>> EIGEN1: Eigenvalue determination'
c ------------------------------------------------------------------------
      NSHELL = NSHELL1(0)
c ------------------------------------------------------------------------
      write(6,9010) IPF,IPFE,IPE
      write(6,9020) IOPREF,IOPSYS
      write(6,9030) LMAX,LMMAX,LPOT,NSPIN
      write(6,9040) NREF,NATYP,NAEZ,NINEQ,NSHELL,NZ,NLAYER,NPNT2
      write(6,9050) INS,KSHAPE,KVREL,KSCOEF,KHFELD
      write(6,9060) KAOEZ(NAEZ),REFPOT(NATYP)
 9010 FORMAT(' IPF IPFE IPE  :',3i4)
 9020 FORMAT(' IOPREF IOPSYS :',2i4)
 9030 FORMAT(' LMAX LMMAX LPOT NSPIN :',4i4)
 9040 FORMAT(' NREF NATYP NAEZ NINEQ NSHELL NZ NLAYER NPNT2 :',8i4)
 9050 FORMAT(' INS INSREF KSHAPE KVREL KSCOEF KHFELD        :',6i4)
 9060 FORMAT(' KAOEZ(NAEZ) REFPOT(NATYP) :',2i4)
c ------------------------------------------------------------------------
      EBOT = EB
      LSTART = .TRUE.
      CONVPU = ALATC / TPI
      POLBOUND = QBOUND/20.0d0
      ITOT = 0
      E2IN = 0.0d0
      LSYM = .TRUE.
      IF (OPT('COMPLEX ')) LSYM = .FALSE.
c
c ---> test of minimum and maximum band indices
c
      IF (MMIN.LE.0) MMIN = 1
      IF (MMAX.LE.0) MMAX = LMMAX*NAEZ
      IF (MMAX.LT.MMIN) MMAX = MMIN
      IF (MMAX-MMIN+1.GT.M2D) THEN
        write(6,*) 'Please, change the parameter m2d in',
     *   ' inc.fi to',MMAX-MMIN+1
        stop 'in EIGEN'
      END IF
c ------------------------------------------------------------------------
c
c ---> input of potential cards for reference system () 
c      and real system ()
c
      call STARTB1(IOPREF,IPF,IPFE,IPE,KVREL,KWS,0,LMAX,
     +     1,NREF,
     +     ALATNEW,RMTNEW1(1),RMT1(1),      ! KMT,MTFAC,
     +     ITITLE,HFIELD,IMT1(1),IRC1(1),0.0d0,
     +     INSREF,IRNS(1),LPOT,NSPIN,VINS1(IRMIND,1,1),
     +     IRMIN1(1),KSHAPE,NTCELLR(1),
     +     IRCUT1,IPAN1,THETAS1,IFUNM1,
     +     NFU1(1),LLMSP1,LMSP1,E2IN,
     +     VBC,C,DROR1,RS1,S1,VISP1,RWS1,
     +     ECORE1,LCORE1,NCORE1,DRDI1,
     +     R1,Z1(1),A1(1),B1(1),IRWS1(1),INIPOL,1)
c
      write(6,*) 'VISP REF'
c      
      call STARTB1(IOPSYS,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +     1,NATYP,
     +     ALATNEW,RMTNEW(1),RMT(1),      ! KMT,MTFAC,
     +     ITITLE,HFIELD,IMT(1),IRC(1),0.0D0,
     +     INS,IRNS(1),LPOT,NSPIN,VINS(IRMIND,1,1),
     +     IRMIN(1),KSHAPE,NTCELL(1),
     +     IRCUT,IPAN,THETAS,IFUNM,
     +     NFU(1),LLMSP,LMSP,E2IN,
     +     VBC,C,DROR,RS,S,VISP,RWS,
     +     ECORE,LCORE,NCORE,DRDI,
     +     R,Z(1),A(1),B(1),IRWS(1),INIPOL,1)

      write(6,*) 'VISP SYS'
      write(6,FMT='('' E2IN    : '',f15.6)') E2IN
c      write(6,fmt='(8f9.3)') (visp(i,1),i=1,irws(1))
c ------------------------------------------------------------------------
      write(6,*) 'MMIN,MMAX :',MMIN,MMAX
c
c ---> correct order of EBOT, E1 and E2
c
      write(6,FMT='('' EBOT,E1,E2 :'',3f12.6)') EBOT,E1,E2
      IF (.NOT.(EBOT.LE.E1 .AND. E1.LE.E2)) 
     +     STOP 'EIGEN1 : order EB,E1,E2'
c ------------------------------------------------------------------------
c
c ---> create complex E-mesh on the real axis between 
c      (EBOT,0) and (E2,0)
c
      IF (NPNT2.GT.NPNTD) then
        write(6,*) 'Please, change the parameter npntd in',
     *       ' inc.fi to a value greater equal ',npnt2
        CALL RCSTOP('NPNTD   ')
      ENDIF
c
      EIMAG = 0.d-5 ! ZERO
      IF (NPNT2.gt.1) THEN
        DO I = 1,NPNT2
          EZ(I)= CMPLX(EBOT,EIMAG) 
     +         + CMPLX(E2-EBOT,EIMAG)*(I-1)/(NPNT2-1)
        END DO
      END IF
c ------------------------------------------------------------------------
      IF (OPT('wfct    ')) THEN
        CALL EIGWFCT(NSPIN,NAEZ,NINEQ,LMAX,LSYM,RBASIS,KAOEZ,QBOUND)
        GOTO 999
      END IF                        ! (OPT('wfct    '))
c ------------------------------------------------------------------------
      DO 120 ISPIN = 1,NSPIN
        NZERO(ISPIN) = 0
 120  END DO
c ------------------------------------------------------------------------
c
c ---> open file 'kpoints' with mesh points
c
      OPEN(30, FILE='kpoints',STATUS='OLD')
      READ(30,*) NPOINTS,AKSCALE,BKSCALE,CKSCALE,FAC

c      IF (.not.OPT('use fac ')) fac = 1.0d0

      WRITE(6,9005) NPOINTS,AKSCALE,BKSCALE,CKSCALE,FAC
      
      AKSCALE = TPI/AKSCALE ! /FAC

      IF (NPOINTS.GT.KPOIBZ) THEN
        write(6,*) 'Please change the parameter kpoibz in inc.fi ',
     +       'to a value greater than ',npoints
        STOP
      END IF
c ------------------------------------------------------------------------
      IF (TEST('e(k)    ')) THEN
c
c --->  eigenvalues at E-mesh-points
c
        TIMES= DCLOCK()
        DO 60 IP=1,NPOINTS
          write(6,*) 'IP:',ip
c
          II = 3
          IF (OPT('COMPLEX ')) II = 6
          READ(30,9002) (QQ(I),I=1,II)
c
          DO 431 I = 1,II
            Q(I) = QQ(I)/AKSCALE
 431      END DO
c
          IF (OPT('COMPLEX ')) THEN
            WRITE (6,FMT=9006) (QQ(I),I=1,3),(Q(I),I=1,3),
     +           (QQ(I),I=4,II)
          ELSE
            WRITE (6,FMT=9008) (QQ(I),I=1,3),(Q(I),I=1,3)
          END IF
c
          DO 105 ISPIN=1,NSPIN
            write(6,*) 'ISPIN :',ISPIN
            DO 100 I = 1,NPNT2
              CALL EIGEN0(ISPIN,Q,EZ(I),LSYM,0)
              if (.not.test('fix sym ') ) goto 220
              CALL EIGEN0(ISPIN,Q,EZ(I),.not.LSYM,0)
              
 220          IF (I.EQ.1) NBASIS(IP) = NUM

c              NUM = NUM + NPOLES(DREAL(EZ(I)),ISPIN)-NBASIS(IP)
c              write(6,FMT='(''E , N :'',f13.6,i6)') dreal(ez(i)),num

              write(6,FMT='(''E , N :'',f13.6,i6,1p,2d12.2)') 
     +             dreal(ez(i)),num,WDO,WUP

 100        END DO                  ! I = 1,NPNT2
 105      END DO                    ! ISPIN=1,NSPIN
 60     END DO                      ! IP=1,NPOINTS
        write(6,FMT='('' time in E(K)   :'',f12.2,'' sec.'' )') 
     +       DCLOCK()-TIMES
        GOTO 999
      END IF                        ! (TEST('e(k)    '))
c ------------------------------------------------------------------------
      II = 3
      IF (OPT('COMPLEX ')) II = 6
      DO 40 IP=1,NPOINTS
c
c --->  read and normalize all k-points
c
        READ(30,9002) (QQ(I),I=1,II)
c
        DO 430 I = 1,II
          QKP(I,IP) = QQ(I)/AKSCALE
 430    END DO

        IF ( TEST('k-net    ') .OR. 
     +       (IP.LE.5) .OR. 
     +       (IP.GT.NPOINTS-5)     ) THEN
          IF (OPT('COMPLEX ')) THEN
            WRITE (6,FMT=9007) 
     +           IP,(QQ(I),I=1,3),(QKP(I,IP),I=1,3),(QQ(I),I=4,II)
          ELSE
            WRITE (6,FMT=9009) IP,(QQ(I),I=1,3),(QKP(I,IP),I=1,3)
          END IF
        END IF

 40   END DO
c ------------------------------------------------------------------------
      IF (OPT('iso surf')) THEN
        IF (E2IN.NE.ZERO .AND. .NOT. TEST('fix E2  ')) E2 = E2IN

        CALL EIGISO(HEAD,NSPIN,NAEZ,NINEQ,LMAX,LSYM,E2,
     +       VCONST,QKP,NPOINTS,QBOUND,AKSCALE)
        GOTO 999
      END IF                        ! (OPT('iso surf'))
c ------------------------------------------------------------------------
c
c ---> calculate EIGENVALUES of the bands MMIN to MMAX at given K-POINTS
c
      EFCTOR = 1.D0
      IF (TEST('EV      ')) EFCTOR = 13.6058D0
c
      OPEN(31, FILE='evonkdown',STATUS='UNKNOWN')
      OPEN(32, FILE='evonkup',STATUS='UNKNOWN')

      WRITE(31,102) NPOINTS
      WRITE(32,102) NPOINTS
 102  FORMAT('!! number of k-points :',I6)

      TIMET = DCLOCK()
      DO 50 ISPIN = 1,NSPIN
        TIMES = DCLOCK()
        WRITE(6,*) 'ISPIN :',ISPIN
        CALL EINIT(NPOINTS)
c
        DO 250 IP=1,NPOINTS
          NBASIS(IP) = 0
 250    END DO
c
        DO 260 IP=1,NPOINTS
c
c --->    determine NMIN
c
          IF (DABS(E1-EBOT).GT.1.0D-10) THEN
            E = CMPLX(E1,ZERO)  
            CALL EIGEN0(ISPIN,QKP(1,IP),E,LSYM,0)
c            NMIN(IP) = NUM + NPOLES(E1,ISPIN) - NBASIS(IP)
            NMIN(IP) = NUM
          ELSE
            NMIN(IP) = 0
          END IF
          CALL ESAVE(E1,NMIN(IP),IP)
 260    END DO

        DO 270 IP=1,NPOINTS
c
c --->    determine NMAX
c
          E = CMPLX(E2,ZERO)
          CALL EIGEN0(ISPIN,QKP(1,IP),E,LSYM,0)
c
c          NMAX(IP) = NUM + NPOLES(E2,ISPIN) - NBASIS(IP)
          NMAX(IP) = NUM
          CALL ESAVE(E2,NMAX(IP),IP)

          IF (TEST('nbasis  '))
     +         write(6,FMT='('' IP,NBASIS,NMIN,NMAX :'',I7,3I5 )')
     +         IP,NBASIS(IP),NMIN(IP),NMAX(IP)

 270    END DO
c
c --->  test maximum number of bands in energy intervall (E1,E2)
c       at all k-points
c
        LTEST = .FALSE.
        I = 0
        DO 420 IP=1,NPOINTS
          IF (NMAX(IP)-NMIN(IP) .GT. M2D) THEN
            LTEST = .TRUE.
            I = MAX(I,NMAX(IP)-NMIN(IP))
          END IF
 420    END DO
        IF (LTEST) THEN
          write(6,*) 'Please, change the parameter m2d in',
     *         ' inc.fi to',i
          STOP 'M2D in eigen1'
        END IF
c
c --->  E-SCAN
c
        write(6,*) 'E-SCAN :'
        D = 5*NACLS(1)/NAEZ
c
c --->  not optimized E-SCAN at 2**N E-points
c
c        N = 2*1.4427*LOG(1.4427*D*NPOINTS*(MMAX-MMIN+1)/(D+NPOINTS))+1
c
c --->  optimized E-SCAN at 2**N E-points
c

c        N = 1.4427*NPOINTS*(MMAX-MMIN+1)*(D-1)/D
c        N = MIN(N,I)
 
        I = 1.4427*LOG(1.d0*(E2-E1)/QBOUND)
        N = MAX(I,1)

        IF (TEST('fix NE  ')) THEN
          N = NPNT2
          write(6,*) 'NE for opt. E-scan is changed to ',N
        END IF
c
c --->  create E-MESH for E-SCAN
c
        NPNT = 0
        II = 1
        DE = E2-E1
        DO 400 I = 1,N
          DE = DE/2.0D0
          DO 410 J = 0,II-1
            NPNT = NPNT+1
            IF (NPNT.LE.NPNTD) THEN
              EZ(NPNT) = CMPLX(E1+DE*DBLE(1+2*J),ZERO)
              IF (DABS(DREAL(EZ(NPNT))).LT.1.0D-4) 
     +             EZ(NPNT) = CMPLX(1.0D-4,DIMAG(EZ(NPNT)))
              IF (TEST('e-scan1 ')) 
     +             WRITE (6,FMT='(I6,f14.6)') NPNT,DREAL(EZ(NPNT))
            END IF
 410      END DO
          II = 2*II
 400    END DO

        IF (NPNT.GT.NPNTD) then
          write(6,*) 'Please, change the parameter npntd in',
     +         ' inc.fi to a value greater equal ',npnt
          CALL RCSTOP('NPNTD   ')
        END IF

        write(6,*) 'Points for E-SCAN used : ',NPNT,N
        
        TIMEE = DCLOCK()
        ITOT = 0
        IE = 0

        DO 280 I = 1,NPNT
        EE = DREAL(EZ(I))
C
        IF (TEST('e-scan  ')) THEN 
          write(6,FMT='('' IE,EE :'',I8,F15.6 )') I,EE
        ELSE
          IF (MOD(I,NPNT/80+1).EQ.0) write(6,FMT='(''*'',$)') 
        END IF
C
        ltest=.false.
        DO 290 IP = 1,NPOINTS
          CALL EREST1(IP,EE,IO)
          IF (TEST('simpl E ')) IO = 0 ! simple E-scan
          IF (IO.LT.1) THEN
            ltest=.true.
            ITOT = ITOT + 1
            CALL EIGEN0(ISPIN,QKP(1,IP),EZ(I),LSYM,0)
            N = NUM + NPOLES(EE,ISPIN) - NBASIS(IP)
            CALL ESAVE(EE,N,IP)
            IF (TEST('e-scan  ')) WRITE(6,*) 'IP,N  :',IP,N
          END IF                    ! (IO.LT.1)
 290    END DO                      ! IP = 1,NPOINTS
        IF (LTEST) IE = IE + 1
c
 280    END DO                        ! I = 1,NPNT
c
        write(6,FMT='('' ITOT in E-SCAN :'',I12)') ITOT
        write(6,FMT='('' IE   in E-SCAN :'',I12)') IE
        ITOT = 0
c
        write(6,FMT='('' time in E-SCAN :'',f12.2,'' sec.'' )') 
     +       DCLOCK()-TIMEE
c
        TIMEI = DCLOCK()
c ------------------------------------------------------------------------
c
c --->  looop over k-mesh points
c
 8000   DO 300 IP = 1,NPOINTS
          write(6,*) 'IP :',IP
c
c --->    init step counter
c
          DO 140 I=1,NAEZD*LMMAXD
            ISTEP(I) = 0
 140      END DO
c
c --->    init of lowest and highest band in intervall E1, E2
c         at every k-point
c
          IMIN = MMIN
          IMAX = MMAX
          IF (MMIN.LE.NMIN(IP)) THEN
            DO I=MMIN,NMIN(IP)
              EQ(I) = E1
            END DO
            IMIN = NMIN(IP) + 1
          END IF

          IF (MMAX.GT.NMAX(IP)) THEN
            DO I=NMAX(IP)+1,MMAX
              EQ(I) = E2
            END DO
            IMAX = NMAX(IP)
          END IF

          IF (TEST('ip-scan ')) write(6,*) 'IMIN,IMAX :',IMIN,IMAX
c
          I = IMIN
c ------------------------------------------------------------------------
c
c --->  loop for eigenvalue search
c
        DO 10 WHILE (I.LE.IMAX)
c
          IF (TEST('ip-scan ')) THEN
            write(6,*) 'I :',I
          ELSE
            WRITE(6,FMT='(''I :'',I7,$)') I
          END IF
c
c --->    determine energy intervall for incremental searching
c
          CALL EREST(IP,EA,EE,I,NUM1,NUM2,W1,W2,IVAL1,IVAL2)
          IF (TEST('ip-scan ')) 
     +           write(6,fmt='(''ea,ee :'',2f14.9)') ea,ee

          N2 = NUM2
          FAC = 1.0d0
          IF ((ABS(EA)+ABS(EE)).LT. 0.1D0 ) FAC = 3.D-2
          DO 20 WHILE (EE-EA .gt. FAC*QBOUND)

            EM = (EA + EE)/2.0D0
            IF (DABS(EM) .LT. 1.0d-9) EM = (EA+EM)/2.0D0
            E = CMPLX(EM,ZERO)
            CALL EIGEN0(ISPIN,QKP(1,IP),E,LSYM,0)
            
            N = NUM + NPOLES(EM,ISPIN) - NBASIS(IP)

            CALL ESAVE(EM,N,IP)

            IF (TEST('ip-scan ')) 
     +           write(6,fmt='(''em,n :'',f14.9,i6)') em,n

            IF (N.GE.I) THEN
              EE    = EM
              N2    = N
              W2    = WDO
              IVAL2 = IDO
            ELSE
              EA    = EM
              W1    = WUP
              IVAL1 = IUP
            END IF
            
            ISTEP(I) = ISTEP(I) +1 ! step counter
            ITOT  = ITOT + 1
 20       END DO

          IF (TEST('ip-scan ')) THEN
            write(6,FMT='('' EA,EE :'',2f14.9 )') ea,ee
            write(6,FMT='('' VAL,W :'',2I6,2d14.3)') ival1,ival2,w1,w2
          END IF
  
          IF (IVAL1*IVAL2.GT.0) THEN
            EM = (W1*EE-W2*EA)/(W1-W2)
          ELSE
            EM = (EA+EE)/2.0D0
c            IF (QBOUND.GT. 1.0D-6) STOP 'Decrease QBOUND in in5 !!'
          END IF
c
          IF (N2.GT.IMAX) N2=IMAX
c
          DO 30 J = I,N2
            EQ(J) = EM
 30       END DO
          write(6,FMT='(f16.10)') EQ(I)
          IF (N2.gt.I) write(6,FMT='((10x,f16.10))') (EQ(J),J = I+1,N2)
c 
          I = N2 + 1
 10     END DO                      ! WHILE (I.LE.IMAX)

        IF (TEST('ip-scan ')) 
     +       write(6,FMT='(i5,F18.6,i8)')
     +       (i,EQ(I),ISTEP(I),I=MMIN,MMAX)

        IF (TEST('sub ef  ')) THEN
          DO I = MMIN,MMAX
            EQ(I) = EQ(I)-E2IN
          END DO
        END IF

        IF (TEST('larg_acc')) THEN
          write(30+ISPIN,9004)  (EQ(I)*EFCTOR,I=MMIN,MMAX)
        ELSE
          write(30+ISPIN,9003)  (EQ(I)*EFCTOR,I=MMIN,MMAX)
        END IF
        
 300    END DO                        ! IP=1,NPOINTS

        write(6,FMT='('' time in IP    loop :'',f12.2,'' sec.'' )') 
     +       DCLOCK()-TIMEI
        write(6,FMT='('' time in ISPIN loop :'',f12.2,'' sec.'' )') 
     +       DCLOCK()-TIMES

 50   END DO                        ! ISPIN = 1,NSPIN
c
      write(6,*) 
     +     'Total number of steps ',ITOT,
     +     ' for ',NPOINTS,' kpoints and ',MMAX-MMIN+1,' bands.'
c
      close(30)
      close(31)
      close(32)

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 999  IF (TEST('flow    ')) WRITE(6,*) '<<< EIGEN1'

      RETURN

 9002 FORMAT((10X,3F10.6))
 9003 FORMAT(20F9.5)
 9004 FORMAT(10F16.12)
 9005 FORMAT('npoints,akscale,bkscale,ckscale,fac : ',I7,4f8.3)
 9006 FORMAT(  10X,3F10.6,5x,3F10.6,/,10X,3F10.6)
 9007 FORMAT(i6,4X,3F10.6,5x,3F10.6,/,10X,3F10.6)
 9008 FORMAT(  10X,3F10.6,5x,3F10.6)
 9009 FORMAT(i6,4X,3F10.6,5x,3F10.6)

      END
