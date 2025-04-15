      SUBROUTINE FPIMPU(IHANDLE,FILNAM)
c	include 'cxml_include.f90'
      IMPLICIT NONE
c                -------------------------------------
c               |                                     |
c               |   i m p u r i t y - p r o g r a m   |
c               |                                     |
c               ---------------------------------------
c        this program is used to perform selfconsistent,all-electron,
c        scalarrelativistic calculations of the electronic structure
c        of a defect-cluster of many shells  perfectly embedded in a
c        ordered alloy - taking into account the exact shape of each
c        cell
c
c
c   original version : cluster of 13 atoms, complex energy
c                                                   by r.zeller 1982
c   further developments : core-relaxation,scalarrelativistic,
c                          exponential mesh,high accuracy for radial-
c                          wavefunction,hyperfinefield and isomershift
c                                                   by s.bluegel 1984
c                          extension to f-electrons in greenowsfunction
c                                                   by r.zeller  1984
c                          use of symmetrized greensfunction
c                                                   by r.zeller  1986
c                          extension to different reference atoms,more
c                          shells,dynamic-programming and vectorization
c                          single-site boundary condition for core,
c                          hyperfine fields and isomershifts
c                                                   by s.bluegel 1986
c                          broyden update method implemented
c                                                   by s.bluegel 1987
c                          non-spherical potentials and total energy
c                          calculation,lloyd's formular
c                                                  by b.drittler 1987
c                          correct boundary conditions for the wave-
c                          functions,software reconstruction
c                                                  by b.drittler 1988
c                          non-spherical input potentials (nonsra and
c                          s.r.a.),electric field gradients
c                                                  by b.drittler 1989
c                          exact cell treatment
c                                  by b.drittler and u. klemradt 1989
c
c                          Reconstruction of the code using blas
c                          routines
c                                                 by R. Zeller   1994
c
c                          Implementation of lattice relaxations
c                                            by N. Papanikolaou  1996
c
c                          Extended for more than one host atom
c
c                          Implementation of rigid shift of potential
c                          and void methods for lattice relaxations,
c                          (see: klatr key)
c                                                  by A. Settels 1997
c
c                          Calculation for special representations
c                          possible.
c                          This can be used to calculate different
c                          charged defect levels in the gap.
c                          Charges and Energies are tested.
C                          DOS for special representations.
c
c----------------------------------------------------------------------
c                   ****** legende *****
c    cfg       : configuration of valence states 4,4,3,0= 4s,4p,3d,0f,
c                compare with kfg in subroutine corel
c
c   khfeld=(0,1)  key for applying an external magnetic  hfield
c                 on the impurity : 0 means no field
c                                   1 means hfield is read in & applied
c   kvrel=(0,1,2) key valence electron : 0 means calculating valence
c                                          electrons by schroedingereq
c                                        1 means calculating valence
c                                          electrons by sra for cluster
c                                        2 means calculating valence
c                                          electrons by sra
c   kcor=(0,1,2,3) :key of calculating core-relaxation
c                                        0 means no
c                                        1 means yes by sra for cluster
c                                        2 means yes by schroedingereq.
c                                        3 means yes by sra
c   khyp=(0,1) :key of calculating hyperfinefields + isomerieshifts +
c                                  nuclear-spin-lattice relaxation-time
c                                        0 means no
c                                        1 means yes
c   kws =(0,1,2) :  key for doing ws or mt-calculation
c                                        0 means mt for all atoms
c                                        1 means ws for the cluster and
c                                                  mt for the host
c                                        2 means ws for cluster and host
c   kfsr=(0,1) :  key of calculating friedel's screenning rule
c                                        0 means no
c                                        1 means yes
c   kxc=(0,1,2) : key for using different exchange-correlation potential
c                      0 means moruzzi,janak,williams exchange (default)
c                      1 means von barth,hedin exchange
c                      2 means vosko,wilk,nusair exchange
c   kte=(0,1)   : key for calculating total energies
c                                        0 means no calculation
c                                        1 means calculation
c   kefg=(0,1)  : key for calculating electric field gradients
c                                        0 means no calculation
c                                        1 means calculation
c   kf=(0,1)    : key for calculating forces
c                                        0 means no calculation
c                                        1 means calculation
c   ksph=(0,1)  : key for calculating change of specific heat
c                                        0 means no calculation
c                                        1 means calculation
c   igsym=(0,1) : key for using symmetrized host green's functions
c   kshape=(0,1,2) : key for exact cell treatment (shape corrections)
c                     0 means always atomic sphere approximation (asa)
c                     1 means start with asa , but end with exact cell
c                             treatment
c                     2 means exact cell treatment
c                     3/4 not working in this version
c                     5 use shape with kinks for host, equally spaced
c                       shapes for the cluster
c                     6 use equally spaced shapes everywhere
c   ins=(0,1,2,3) : key for non spherical input potential
c                     0 means always spherical averaged input potential
c                     1 means start with spherical averaged input poten-
c                             tial , but end with non sph. input poten-
c                             tial for the cluster
c                     2 means spherical averaged input potential for
c                             host,non sph. input potential for cluster
c                     3 means non spherical input potential for host
c                             and cluster
c   icst=(0,..,n) : key for number of born approximations
c                                  0 means zero-th. born approximation
c                                  1 means first born approximation
c                                  n means n-th. born approximation
c   iprcor=(0,1,2): key for printinformation after using core-relaxation
c                                        0 means no information
c                                        1 means short information
c                                        2 means large information
c   ipe=(0,1): key for printing the information of the last iteration
c              of an different file
c                                        0 means no
c                                        1 means yes
c   icut =(0,1,2,...) : key for cutting off some shells
c                                        0 no cut off
c                                        1 after the first shell
c                                        2 after the second shell
c                                        i after the i.th   shell
c   imix=(0,..,5) : key for using different mixing schemes
c                                        0 means straight mixing
c                                        1 means chebycheff used after
c                                                19th. iteration
c                                        2 means chebycheff used right
c                                                away
c                                        3 broyden's first method used
c                                        4 broyden's second method used
c                                        5 generalized andersen method
c   klatr=(0,1)   : key for lattice relaxation
c                   present version works only for full potential
c                                   0 means no lattice relaxation
c                                   1 means lattice relaxation
c                                      (in this case read new lattice
c                                       coordinates (look reainp))
c             1: 'old method'   _   ~o   ~o     ~o _
c                               G = G  + G (t - t )G
c
c             2: rigid shift    _o  ~o   ~o  o  ~o _o
c                of pots        G = G  + G (t - t )G
c                               _   _o   _o      o _
c                               G = G  + G (t - t )G
c
c             3: void      v   o    o      o  v       ~v     v T
c                         G = G  + G (0 - t )G   and  G = U G U
c
c                         _   ~v   ~v       _
c                         G = G  + G (t - 0)G
c
c   ksymm=(0,1..??): key for symmetrization procedure
c                                        0 means no symmetry is used
c                                        1 means symmetrization coefs
c                                          from R. Zeller for bulk
c
c
c                                            ver. 10.6.1996
c
c   ksymmat 0/1     :   0 - nothing happens
c                       1 select special symmetry
c                       representations to calculate.
C
C   >>>>> USED ONLY IF KSYMMAT=1 >>>>>
c   kesym   1-ielast:   This is the first energy point
c                       using this selection for special
c                       representations.
c   ASARY           :   is only read in in case of ksymmat=1
c                       the input is between klatr...  and
c                       the relaxation vectors.
c                       Only the calculated spin directions are
c                       read in.
C                       0 - NEGLECT THIS REPRESENTATION
C                       1 - USE THIS     REPRESENTATION
C   EXAMPLE for Oh SYMMETRY ==> 10 representations
C   and NSPIN=2:
C   ...
C   1    1    0    1    1               klatr,ksymm,lkonv,ksymmat,kesym
c   1    1    1    1    1    1    1    1    ASARY FOR SPIN DOWN 1..NREP
C   1    1                                  ASARY FOR SPIN DOWN 1..NREP
c   1    1    1    1    1    1    1    1    ASARY FOR SPIN UP   1..NREP
C   1    1                                  ASARY FOR SPIN UP   1..NREP
c   .00000    .00000    .00000    .00000    .00000    .00000
C   ...
C   HERE: THIS EXAMPLE SAME RESULT WITH KSYMMAT=0
C   <<<<< USED ONLY IF KSYMMAT=1 <<<<<
c
c--------------------------------------------------------------------
c     Modified by T. Korhonen,          Sep. 1997
c     See the README file for details.
c
c     This version uses equally spaced shape functions for the
c     cluster and the host (kshape=6) or host is calculated
c     using the old shapes with kinks (kshape=5).
c
c     Input shapes: old shapes with kinks. Number of radial mesh
c              points may be large (parameter IRIKD). Host
c              shapes should have IRID points, not IRIKD,
c              if kshape=5 (i.e., shapes with kinks are
c              used for host).
c
c     Input potentials: cluster potential in equally spaced mesh
c              in shape area, host equal mesh or mesh
c              with kinks (kshape=5/6). If host potential
c              is in a mesh with kinks (kshape=5), then
c              there should not be smeared spherical
c              potential for host. Otherwise, one should
c              have smeared spherical potential on the
c              disk also.
c
c     IRMD, IRID:   number of points in the equally spaced mesh
c     IRMKD, IRIKD: max. number of points in the original shapes,
c                   which are read in
c     *** IRMKD .ge. IRMD and IRIKD .ge. IRID works now ***
c
c     Now uses spline fitting to get the shapes in the equally spaced
c     mesh from the shapes in the mesh with kinks (l>0). For the l=0
c     component: use smearing with a tiny alpha and correct the volume
c     to the input shape volume. Smeared shapes need in the calculation
c     of the smeared spherical potential are calculated only for the
c     values l <= lmpot.
c---------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NATOMD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD,
     +            NSPIND=2)
      INTEGER IEMXD,NSEC,NREPD
      PARAMETER (iemxd=150,nsec=689,NREPD=4)
      INTEGER IRMD,IRNSD,LMAXD,LPOTD,LMX
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER NFUND,IRID,NGSHD
      PARAMETER (NFUND=289,irid=435,NGSHD=54287)
c
      INTEGER IRMKD,IRIKD
      PARAMETER (irmkd=1484,irikd=435)
c
      INTEGER NCELLD,IPAND
      PARAMETER (NCELLD=20,IPAND=80)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IJD,LASSLD
      PARAMETER (IJD=434,LASSLD=4*LMAXD)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMX-1)**2)
      COMPLEX*16 CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
      REAL*8 ONE,ONEM
      PARAMETER (ONE=1.D0,ONEM=-1.D0)
C     ..
C     .. Scalars in Common ..
      REAL*8 C
      INTEGER LKONV,KGV,NSTOP
      INTEGER IEND,LMAXP1
      LOGICAL STRCO1
      CHARACTER*5 LMACH
C     ..
C     .. Arrays in Common ..
      COMPLEX*16 FZ(IRMD,LMX,NATYPD),MASS(IRMD,NATYPD),
     +           PZ(IRMD,LMX,NATYPD),QZ(IRMD,LMX,NATYPD),
     +           SZ(IRMD,LMX,NATYPD),TMAT(LMX,NATYPD)
      REAL*8 A(NATYPD),B(NATYPD),CLEB(NCLEB,2),DRDI(IRMKD,NATYPD),
     +       R(IRMKD,NATYPD),RWS(NATYPD),RWSM1(NATYPD),
     +       S(0:LMAXD,NATYPD),VISP(IRMD,NPOTD),Z(NATYPD)
      INTEGER ICLEB(NCLEB,4),IRT(NATYPD),IRWS(NATYPD),
     +        JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(LM2D)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 DET,DF,E,EK,DETH,DFLLOYD
      REAL*8 ALAT,ARG,BRYMIX,DNEINT,DSELOC,EFERMI,ETIME,FCM,DNEIHO,FPI,
     +       HFIELD,MIXING,PI,QBOUND,RFCTOR,RFPI,RMSAVM,RMSAVQ,STIME,
     +       STRMIX,VBCC,VCONSTH,XXX,YYY,ZERO,ZZZ,ETIME2,STIME2,
     +       ETIMETOT,STIMETOT,TTIMDY,TTIMRH,TTIMES,TTIMSY,TTIMRE,
     +       TTIMCP,TTIMTM,ALFSME,X,Y,ATIM
      REAL*8 STEPL,DELTAMAX,T_DEBYE,EPS_MD,DIMASS(NATOMD),STEPLN
      INTEGER I,I1,I1DISK,I2,IA,IATYP,ICST,ICUT,IDIM,IE,IEF,IELAST,
     +        IFILE,II,IJEND,IKST,IMIX,INFO,INS,IOWFCT,IP,IPE,IPF,IPFE,
     +        IPFIL,IPOT,IPRF,IR,IRC1,IRESIST,IRM,IRMIN1,IROT,IS,ISPIN,
     +        IT,ITC,ITCLST,IWTMAT,J,K,KCOR,KEFG,KF,KFSR,KHFELD,KHYP,
     +        KHYPO,KOSZIL,KPRE,KSHAPE,KSPH,KTE,KVMAD,KVREL,KWS,KXC,L,
     +        LM,LM1,LM2,LMAX,LMAXSQ,LMPOT,LPOT,M,MR,N,NATOM,NATPER,
     +        NATPS,NATREF,NATYP,NEND,NHSPIN,NIR,NP,NPTPS,NPTREF,NREP,
     +        NSPIN,NSPPOT,NSRA,NSTART,KLATR,INDX,L1,L2,KSYMM,IGGA,
     +        LSMEAR,NPTS,NPTS1,IPR,INDX_REF,IREF_LAST
      INTEGER IHANDLE,IHOST,IHANDLE2,IHANDLE3,KOCCUP
      LOGICAL SURFACE
CKGV--->INTEGER DEFINITION
      INTEGER NREPLTL,NATLTL,LHIGH,LLOW,NDIMHIGH,NDIMLOW,NICHEK,
     +        P4TO3(NSEC)
      LOGICAL NONSPH,WFDISK,FXDRKEY
C     ..
C     .. Local Arrays ..
      COMPLEX*16 ALPHAHOST(LMX,NATYPD,IEMXD,NSPIND)
      COMPLEX*16 ADET(LMMAXD,LMMAXD),ALPHA(LMX,NATYPD),
     +           ALPHA1(LMX,NATYPD),AR(LMMAXD,LMMAXD,NATYPD),
     +           CR(LMMAXD,LMMAXD,NATYPD),DEN(IEMXD,0: (LMX-1),NATYPD,
     +           NSPIND),DETALF(NATYPD,IEMXD,NSPIND),
     +           REFHOST(LMMAXD,LMMAXD,NTREFD),
     +           DTB(LMMAXD,LMMAXD,NTPERD),DTMTRX(NSEC,NSEC),EKL(LMX),
     +           GJHOST(LMX,2,IEMXD,NTREFD),GJM(LMX,NATYPD),
     +           GMAT(NSEC,NSEC),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,
     +           NATYPD),QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,NATYPD),
     +           TMATLL(LMMAXD,LMMAXD,NATYPD),
     +           TOMTOR(LMMAXD,LMMAXD,NATYPD)
c     use next lines, if WFDISK = true
c     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,1),
c     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,1),
      COMPLEX*16 ADETH(LMMAXD,LMMAXD),DETALFH(NATYPD,IEMXD,NSPIND)
      COMPLEX*16 DETJJM1(NATYPD),TMATREL(LMMAXD,LMMAXD,NATYPD),
     +           TMATHOST(LMMAXD,LMMAXD,NATYPD)
      REAL*8 AMATI(NTPERD,NTPERD,LMPOTD,LMPOTD),
     +       BMATI(NTPERD,NTPERD,LMPOTD),RM1(3,NATOMD),DX(3),
     +       RELAX(NTPERD),DRM(3,NTPERD),QOCCUP(NSPIND,NREPD)
      REAL*8 AMAT(NTPERD,NTPERD,LMPOTD,LMPOTD),ATWGHT(NATYPD),
     +       AVMAD(NTREFD,NTREFD,LMPOTD,LMPOTD),
     +       BMAT(NTPERD,NTPERD,LMPOTD),BVMAD(NTREFD,NTREFD,LMPOTD),
     +       CMINST(LMPOTD,NATYPD),CMOM(LMPOTD,NATYPD),ECORE(20,NPOTD),
     +       CMAMV(LMPOTD,NATYPD),ONSITE(LMPOTD,NATYPD),
     +       ECOU(0:LPOTD,NATYPD),EF(2),EGR(NATYPD),EPOTIN(NATYPD),
     +       ESPC(0: (LMX-1),NATYPD,NSPIND),
     +       ESPV(0: (LMX-1),NATYPD,NSPIND),EXC(0:LPOTD,NATYPD),
     +       FLM(-1:1,NATYPD),FLMC(-1:1,NATYPD),GSH(NGSHD),
     +       IEN(IEMXD,NSPIND),IENI(IEMXD,NSPIND),
     +       RHO2NS(IRMKD,LMPOTD,NATYPD,NSPIND),RHOC(IRMD,NPOTD),
     +       RIJ(IJD,3),RM(3,NATOMD),RMT(NATYPD),RMTNEW(NATYPD),
     +       SUMNS(NATYPD),THETAS(IRIKD,NFUND,NCELLD),
     +       THESME(IRID,NFUND,NCELLD),VASUM(0:LMAXD,NATYPD),
     +       THEEQM(IRID,NFUND,NCELLD),VSPSME(IRMD,NPOTD),
     +       VSPSMO(IRMD,NPOTD),VBC(2),VINS(IRMIND:IRMD,LMPOTD,NPOTD),
     +       VONS(IRMKD,LMPOTD,NPOTD),WG(LASSLD),
     +       WORK(IRMKD,LMPOTD,NATYPD),WTYR(IJD,LMPOTD),
     +       WORKEQ(IRMD,LMPOTD,NATYPD),RHOCK(IRMKD,NPOTD),
     +       YR(IJD,LMPOTD),YRG(LASSLD,0:LASSLD,0:LASSLD),
     +       F1XYZ(3,NATOMD),F_OLD(3*NATOMD),TAUXYZ(3,NATOMD)
      REAL*8 TAU(3,NATOMD),TAUPRO(3,NATOMD)
      COMPLEX*16 GTEMP(NSEC,NSEC)
      REAL*8 EFMTZ
CASYM
      REAL*8 GL
      COMPLEX*16 EKLASYM(LMX,NATYPD)
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,NTPERD,NSPIND)
      INTEGER KSYMMAT,KESYM,NSYMMAT,ASARY(NREPD,NSPIND)
      INTEGER CFG(4,NATYPD),IDND(6,NATYPD),IFUNM(NATYPD,LMXSPD),
     +        ILM(NGSHD,3),IMAXSH(0:LMPOTD),IMT(NATYPD),IOPER(NATOMD),
     +        IPAN(NATYPD),IPVT(LMMAXD),IRC(NATYPD),
     +        IRCUT(0:IPAND,NATYPD),IREF(NTPERD),IRMIN(NATYPD),
     +        IRNS(NATYPD),ITITLE(20,NPOTD),LCORE(20,NPOTD),
     +        LLMSP(NATYPD,NFUND),LME(9),LMS(9),LMSP(NATYPD,LMXSPD),
     +        NCORE(NPOTD),ND(48,3,3),NDG(NREPD),NDIM(NREPD),
     +        NFU(NATYPD),NLST(2),NSHELL(NTPERD),NTCELL(NATYPD),KMOLD,
     +        KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,K_EPS_MD,K_ADPT_DT,
     +        K_SET_F_M,IOBROY,IOBROG,IOFILE,KPRI,KTEST,ITDMD,MIT,ITER,
     +        IFATOM(NATOMD),ICNT,IFILEH,IATDYN(NATOMD),KSHAPEH,
     +        KATDYN(NTPERD)
      INTEGER NIMP(NATOMD),SOCCUP(NSPIND,NREPD),QOCCUPH(NSPIND,NREPD)
c     Next are used for the equally spaced radial mesh in the
c     shape function region (used in the solution of the radial
c     equations
      DOUBLE PRECISION REQ(IRMD,NATYPD),DRDIEQ(IRMD,NATYPD),
     +                 DROREQ(IRMD,NATYPD),RSEQ(IRMD,0:LMAXD,NATYPD),
     +                 R2NSEQ(IRMD,LMPOTD,NATYPD,NSPIND),R2NS2D(IRID),
     +                 DROR(IRMKD,NATYPD),RS(IRMKD,0:LMAXD,NATYPD)
      INTEGER IRWSEQ(NATYPD),IRMIEQ(NATYPD),IRCUEQ(0:IPAND,NATYPD),
     +        IPANEQ(NATYPD),IRCEQ(NATYPD)
C
c     Next are needed for Timo's test
      DOUBLE PRECISION EXCEQ(0:LPOTD,NATYPD),VONSEQ(IRMD,LMPOTD,NPOTD),
     +                 CMOMEQ(LMPOTD,NATYPD),CMSTEQ(LMPOTD,NATYPD)
C     ..
CO---> SAVE T-MATRIX FOR SHIFTED POSITION
      INTEGER IOTMAT,IGHWRIT,NOKTE,NOTREAD
C___BZ-INTEGRATION
      INTEGER OCCUP(NATYPD)

      CHARACTER*13 TEXTS(3)
      CHARACTER*24 TXC(5)
      CHARACTER*60 FILNAM(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL AMN,BACSYM,BRYDBM,CHEACC,CONVOL,COREL,CPLXWF,CUTSHELL,
     +         ECOULOM,EFGRAD,EPOTINS,ESINGPC,ESINGPV,
     +         ETOTAL,FORCE,FORCEH,FORCXC,GAUNT,GAUNT2,GDYSON,HYPISO,
     +         LLOYD,MFIELD,MIXSTR,OPROT,REAINP,REDOGH,RESIST,RHOLM,
     +         RHONS,RHOTOT,RINIT,RITES,RMADEL,ROTCUB,SHAPE,SPHERE,
     +         START,SYTMAT,SYTN0,SYTN1,TMATNS,TMREAD,TMWRIT,VBOUND,
     +         VINTERS,VINTRAS,VXCLM,WKSPH,WLDOS,RLXLAT,
     +         VINTREL,SPLINE,SPLINT,MTZSME
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DIMAG,MAX,MAX0,MIN0,MOD,REAL
C     ..
C     .. Common blocks ..
      COMMON /CMACH/STRCO1,LMACH
      COMMON /CRDFUN/PZ,FZ,QZ,SZ,TMAT
      COMMON /DLOG87/VISP,RWS,RWSM1,DRDIEQ,REQ,Z,A,B,IRWSEQ,IRT,LMAXP1
      COMMON /GAUNTC/CLEB,LOFLM,ICLEB,IEND,JEND
      COMMON /MTSOL/MASS,C,S
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. External Functions ..
      REAL*8 DCLOCK
      EXTERNAL DCLOCK
C     ..
C     .. Data statements ..
CO--->CHANNEL FOR T-MATRIX
      DATA IOTMAT/91/
      DATA NOTREAD/0/
      DATA NOKTE/0/
      DATA ICUT/0/
      DATA ZERO/0.0D0/
      DATA TXC/' Moruzzi,Janak,Williams ',' von Barth,Hedin        ',
     +     ' Vosko,Wilk,Nusair      ',' GGA PW91               ',
     +     ' GGA PW86               '/
      DATA TEXTS/' spin down   ',' spin  up    ',' paramagnetic'/
      DATA LMS,LME/1,2,5,10,17,26,37,50,65,1,4,9,16,25,36,49,64,81/

c beretesta
      DATA SURFACE/.false./
c bereteste
C
c     Next are used for timing
CFXDR-----> IN CASE OF NOT USING FXDR-FORMAT:SET OUTPUT CHANNEL
      FXDRKEY = .false.
      write(6,*) ' ihandle fxdrkey',ihandle,fxdrkey
      IGHWRIT = 81
      TTIMRE = 0.d0
      TTIMCP = 0.d0
      TTIMTM = 0.d0
      TTIMSY = 0.d0
      TTIMDY = 0.d0
      TTIMRH = 0.d0
      TTIMES = 0.d0
      STIMETOT = DCLOCK()
      STIME = DCLOCK()
c
      LSMEAR = 0
cholger
      ITER = 0
      MIT = 1
C     ..
      CALL REAINP(EF,RM,IREF,IRNS,NDG,NDIM,NSHELL,NTCELL,TXC,BRYMIX,FCM,
     +            FPI,HFIELD,PI,QBOUND,RFPI,STRMIX,VBC,VCONSTH,ICST,
     +            ICUT,IEF,IFILE,IWTMAT,IMIX,INS,IOWFCT,IPE,IPF,IPFE,
     +            IRESIST,IRM,ITCLST,KCOR,KEFG,KF,KFSR,KHFELD,KHYP,KPRE,
     +            KSHAPE,KSPH,KTE,KVMAD,KVREL,KWS,KXC,LMAX,LMAXP1,
     +            LMAXSQ,LMPOT,LPOT,NATOM,NATPER,NATPS,NATREF,NATYP,
     +            NPTPS,NPTREF,NREP,NSPIN,NSPPOT,KLATR,RELAX,NIMP,RM1,
     +            DRM,KSYMM,IGGA,LKONV,LSMEAR,OCCUP,KSYMMAT,KESYM,ASARY,
     +            SOCCUP,QOCCUP,KMOLD,KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,
     +            K_EPS_MD,K_ADPT_DT,K_SET_F_M,IOBROY,IOBROG,IOFILE,
     +            KPRI,KTEST,ITDMD,STEPL,DELTAMAX,T_DEBYE,EPS_MD,DIMASS,
     +            KSHAPEH,KATDYN,IFILEH,KOCCUP,QOCCUPH)
CASYM-------------------------------------->
CASYM -----> Attention: This is for symmetry only
      NSYMMAT = KSYMMAT
      IF (KSYMMAT.GT.0) THEN
          IF (KSYMM.EQ.0) STOP 'KSYMMAT'
          DTB1 = (0.0D0,0.0D0)
          DO ISPIN = 1,NSPIN
              IF (NOT(FXDRKEY)) THEN
             write(6,*) ' ihandle ',ihandle
                  OPEN (IHANDLE,FILE='green_UNIT',FORM='UNFORMATTED')
                  READ (IHANDLE) I1,I2,EFMTZ,EFMTZ
c	write(6,*) ' i1 i2 efmtz efmtz',i1,i2,efmtz
              ELSE
c                  CALL FXDROPN('green_UNIT','DECODE',IHANDLE)
c                  CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                  CALL FXDRDBL(IHANDLE,EFMTZ,1)
              END IF
              I1 = 0
              I2 = 0
              WRITE (6,FMT=*) NREP,LMAXSQ,KSYMMAT
              DO NP = 1,NREP
                  IF (NOT(FXDRKEY)) THEN
                      CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,IHANDLE)
                  ELSE
c                     CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),GTEMP,
c    +                               IHANDLE)
                  END IF
                  NSYMMAT = ASARY(NP,ISPIN)
                  IF (KSYMMAT.EQ. (NREP+1)) NSYMMAT = NP
                  CALL BACASYM(DTB1(1,1,1,ISPIN),GMAT,LMAXSQ,NREP,NP,
     +                         NSYMMAT,SOCCUP,QOCCUP,ISPIN)
                  NSYMMAT = ASARY(NP,ISPIN)
              END DO
              IF (NOT(FXDRKEY)) THEN
                  CLOSE (IHANDLE)
              ELSE
c                  CALL FXDRCLS(IHANDLE)
              END IF
              DO I2 = 1,NTPERD
                  WRITE (*,FMT=*) 'NTPERD:',NTPERD
                  DO I1 = 1,LMAXSQ
                      WRITE (*,FMT=*) DTB1(I1,I1,I2,ISPIN)
                  END DO
              END DO
              DTB = (0.0D0,0.0D0)
          END DO
      END IF
CASYM -----> END of ASYM-extra
c
      IF (LSMEAR.EQ.0 .AND. (IRMD.NE.IRMKD)) THEN
          WRITE (6,FMT=*) '*** ERROR: lsmear is 0 and irmk .ne. irmkd'
          WRITE (6,FMT=*) '           irmd ',IRMD,', irmkd ',IRMKD
          STOP 'Stop in FPIMPU: irmk .ne. irmkd'
      END IF
c
      KGV = 0
      NSTOP = 0
      IF (KLATR.EQ.4) THEN
          KLATR = 3
C-------->WRITE OUT TILDE(G-VOID)
          KGV = 1
          NSTOP = 1
          WRITE (*,FMT=*) 'WRITE OUT GTRANS AND STOP - Gvoid'
      ELSE IF (KLATR.EQ.5) THEN
          KLATR = 3
C-------->READ ALSO TILDE(G-VOID) IN
          KGV = 2
          WRITE (*,FMT=*) 'READ IN GTRANS, HAS TO BE COPIED! - Gvoid'
      ELSE IF (KLATR.EQ.6) THEN
          KLATR = 2
C-------->WRITE OUT TILDE(G-VOID)
          KGV = 1
          NSTOP = 1
          WRITE (*,FMT=*) 'WRITE OUT GTRANS AND STOP - Rigid Shift'
      ELSE IF (KLATR.EQ.7) THEN
          KLATR = 2
C-------->READ ALSO TILDE(G-VOID) IN
          KGV = 2
          WRITE (*,FMT=*)
     +      'READ IN GTRANS, HAS TO BE COPIED!- Rigid Shift '
      ELSE IF (KLATR.EQ.8) THEN
          WRITE (*,FMT=*) 'WRITE OUT GTRANS AND ONE ITER,TEST'
          KLATR = 3
C-------->READ ALSO TILDE(G-VOID) IN
          KGV = 1
      ELSE IF (KLATR.EQ.9) THEN
          KLATR = 3
C-------->THIS IS ONLY FOR TEST PURPOSE LKVONV=LMAXD
C-------->ATTENTION
          WRITE (*,FMT=*) 'ATTENTION: ONLY FOR TEST PURPOSE!!!'
          KGV = 0
      ELSE IF (KLATR.EQ.10) THEN
          WRITE (*,FMT=*) 'GTRANS HAS TO BE READ IN FROM BZ-INT !!!'
          WRITE (*,FMT=*) 'OCCUPATION HAS TO BE SET FOR T0-MATRIX'
          WRITE (*,FMT=*) 'CHANGE JOBCARD - OCCUP IS READ IN'
      END IF
      IF (KLATR.EQ.2) THEN
          WRITE (*,FMT=*)
     +      'RIGID SHIFT OF POTENTIALS USED FOR RELAXATION'
      ELSE IF (KLATR.EQ.3) THEN
          WRITE (*,FMT=*) 'G-VOID-METHOD USED FOR RELAXATION'
      END IF
      IF (LKONV.NE.LMAXD) THEN
          WRITE (*,FMT=*) 'GREENS FUNCTION HAS LOWER ANGULAR MOMENTUM!'
          WRITE (*,FMT=*) 'CALCULATION UP TO :',LKONV,'!'
          IF (KGV.EQ.0) THEN
              WRITE (*,FMT=*) 'THIS CASE IS NOT IMPLEMENTED'
              WRITE (*,FMT=*) 'WRITE OUT GREENS FUNCTION OR READ IT IN'
              WRITE (*,FMT=*) 'USE VOID METHOD'
              STOP 'KGV'
          END IF
      END IF
c     flushes the buffered output
      CALL FLUSH(6)
c
      CALL ROTCUB(ND)

      DO N = 1,NATPER
          DO II = 1,6
              IDND(II,N) = 1
          END DO
          DO II = 1,7
              READ (35,FMT=9000) IROT
              IF (IROT.EQ.0) THEN
                  GO TO 30
              ELSE
                  IDND(II,N) = IROT
              END IF
          END DO
c
   30     CONTINUE
      END DO

      CALL OPROT(IOPER,NATPER,ND,NSHELL,RM)
      CALL SPHERE(LPOT,YR,WTYR,RIJ,IJEND,IJD)
      CALL GAUNT2(WG,YRG)
      CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND)
      DO 80 NP = 1,NPOTD
          DO 70 I1 = 1,LMPOTD
              DO 60 I2 = IRMIND,IRMD
                  VINS(I2,I1,NP) = ZERO
   60         CONTINUE
   70     CONTINUE
   80 CONTINUE
c     flushes the buffered output
      CALL FLUSH(6)
c
      ETIME = DCLOCK()
      WRITE (6,FMT=*) 'time for initialization ',ETIME - STIME
c
c
c------- end of array set up and definition of input parameter ---------
c
c                     -----------------------------------------
c                   |   do loop over selfconsitency iterations  |
c                     _________________________________________
c
      DO 440 ITC = 1,ITCLST
          NSTART = NATREF + 1
          IF (ITC.EQ.1) NSTART = 1
c
          IF (ITC.EQ.ITCLST .AND. ITCLST.GT.0 .AND. IPE.EQ.1) IPF = IPFE
          IKST = MIN0(ITC-1,1)*NATREF + 1
          KOSZIL = MAX0(-1,-MOD(ITC,5)) + 1
          IF (ITC.EQ.1 .OR. ITC.EQ.ITCLST) KOSZIL = 1
          KHYPO = KOSZIL*KHYP
          IF (ITC.NE.1) IFILE = 0
c
c     11.9.97 Start is ok for different number of points in shapes
          CALL START(IFILE,IPF,IPE,IPFE,KVREL,KWS,LMAX,NATYP,ALAT,RMT,
     +               RMTNEW,ITITLE,IMT,RFCTOR,NATREF,VCONSTH,INS,IRNS,
     +               LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,IRCUT,IPAN,
     +               THETAS,THEEQM,THESME,ALFSME,VSPSME,LSMEAR,IFUNM,
     +               NFU,LLMSP,LMSP,ECORE,LCORE,NCORE,R,DRDI,DROREQ,
     +               RSEQ,IRWS,IRMIEQ,IRCUEQ,IPANEQ,DROR,RS,VBC)


          IF (ITC.EQ.1) THEN

              IF (LSMEAR.GE.1) WRITE (6,FMT=*
     +            ) 'Smearing of spherical potential, alpha= ',ALFSME
              IF ((LSMEAR.EQ.1) .OR. (LSMEAR.EQ.2)) WRITE (6,
     +            FMT=*) 'Use non-smeared potential in 1st iteration'
              IF ((LSMEAR.EQ.3) .OR. (LSMEAR.EQ.4)) WRITE (6,
     +            FMT=*) 'Smeared spherical potential read in'
              IF ((LSMEAR.EQ.1) .OR. (LSMEAR.EQ.3)) WRITE (6,
     +            FMT=*) 'Cluster in equally spaced mesh, host ',
     +            'in original mesh'
              IF ((LSMEAR.EQ.2) .OR. (LSMEAR.EQ.4)) WRITE (6,
     +            FMT=*) 'Equally spaced mesh for host + cluster'
          END IF

c     Smeared shapes now in the array thesme and smeared spherical
c     potential in vspsme (lsmear=3,4). If lsmear=1,2, then vspsme is
c     equal to the non-smeared spherical potential visp.
c     Only in the first iteration can lsmear be equal to 1,2,
c     If host is smeared (lsmear=2,4), then smeared host potential
c     is read in and stored to vspsme.
c
c     Original shapes are in array thetas, shapes interpolated to
c     an equally spaced mesh are in array theeqm and the smeared
c     equally spaced shapes in array thesme ( lsmear=1,2,3 or 4)
          IF (IMIX.GE.3) THEN
              MIXING = BRYMIX
              FCM = 1.0D0

          ELSE


              MIXING = STRMIX
          END IF
c
c     flushes the buffered output
          CALL FLUSH(6)
c
          IF (ITC.EQ.1) THEN
c#######################################################################
c                   this part runs only one time during first iteration
c
              STIME = DCLOCK()
c
              IF (KSHAPE.NE.0) THEN
c
c---> set up of gaunt coefficients
c
                  CALL SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,
     +                       YRG)
c
c---> read coeffients for madelung potential
c
                  IF (.NOT.SURFACE) CALL RMADEL(AVMAD,BVMAD,ALAT,NATREF,
     +                                   LMPOT)

              END IF

              IF (KSHAPE.EQ.0 .AND. NATREF.GT.1) THEN
c
c--> Reading Madelung in case of ASA for more atoms per unit cell
c
                  IF (.NOT.SURFACE) CALL RMADEL(AVMAD,BVMAD,ALAT,NATREF,
     +                                   LMPOT)
              END IF

              DO 90 IATYP = 1,NATYP
                  IRC(IATYP) = IRCUT(IPAN(IATYP),IATYP)
                  IRCEQ(IATYP) = IRCUEQ(IPANEQ(IATYP),IATYP)
   90         CONTINUE
c
              IF (KLATR.EQ.0) THEN
                  CALL AMN(AMAT,ALAT,BMAT,IOPER,LPOT,NATPER,ND,NSHELL,
     +                     RM,WG,YRG,YR,WTYR,RIJ,IJEND)
              ELSE
                  DO I = 1,NATPER
                      DO II = 1,3
                          DRM(II,I) = DRM(II,I)*ALAT
                      END DO
                  END DO
                  CALL AMN(AMAT,ALAT,BMAT,IOPER,LPOT,NATPER,ND,NSHELL,
     +                     RM,WG,YRG,YR,WTYR,RIJ,IJEND)
                  CALL AMN(AMATI,ALAT,BMATI,IOPER,LPOT,NATPER,ND,NSHELL,
     +                     RM1,WG,YRG,YR,WTYR,RIJ,IJEND)
              END IF
c
              IF (KCOR.EQ.0) THEN
c
c---> calculate core densities only at the first iteration
c
                  DO 100 IP = 1,NSPPOT
                      IPFIL = IPF
                      IF (IP.LE.NATREF) IPFIL = IPFE
c     ipr=0 : do not write state dependent information
c     ipr=1 : write something
c     ipr=2 : write all (for debugging)
                      IPR = 0
                      IF (ITC.EQ.1) IPR = 1
c
                      CALL COREL(IPFIL,ITC,40,KHYPO,2,IPR,10,IP,IRM,
     +                           NSPIN,RHOC,NATREF,IMT,KSHAPE,ECORE,
     +                           LCORE,NCORE,CFG)
c              CALL COREL(IPFIL,ITC,40,KHYPO,2,0,10,IP,IRM,NSPIN,RHOC,
c     +                   NATREF,IMT,KSHAPE,ECORE,LCORE,NCORE,CFG)
  100             CONTINUE
              END IF
c
              ETIME = DCLOCK()
              WRITE (6,FMT=*) 'time for 1st iter set up ',ETIME - STIME
c

c     flushes the buffered output
              CALL FLUSH(6)
c
c       itc.eq.1, i.e., only for the first iteration
          END IF
c
c
          IF (KCOR.NE.0) THEN
c     this part is only used for core-relaxation
c
              DO 110 IP = IKST,NATYP
                  IPFIL = IPF
                  IF (IP.LE.NATREF) IPFIL = IPFE
c     ipr=0 : do not write state dependent information
c     ipr=1 : write something
c     ipr=2 : write all (for debugging)
                  IPR = 0
                  IF (ITC.EQ.1) IPR = 1
c
                  CALL COREL(IPFIL,ITC,40,KHYPO,KCOR,IPR,10,IP,IRM,
     +                       NSPIN,RHOC,NATREF,IMT,KSHAPE,ECORE,LCORE,
     +                       NCORE,CFG)

  110         CONTINUE
c
          END IF
c                     ----------------------------------
c                   |   do loop over spins and energies  |
c                     __________________________________
c
          IF (ITC.EQ.1) THEN
c------ host calculation only in first iteration-------
c
              STIME = DCLOCK()
c                     ----------------------------------
              NHSPIN = 1
c new
              IF (NOT(FXDRKEY)) THEN
                  write(6,*) ' ihandle ',ihandle
                  OPEN(IHANDLE,FILE='green',FORM='UNFORMATTED',
     +                                  access='stream')
              ELSE
c                  CALL FXDROPN('green','DECODE',IHANDLE)
              END IF

              DO 200 ISPIN = 1,NSPIN
                  CALL RINIT(IRMKD*LMPOTD*NATREF,RHO2NS(1,1,1,ISPIN))
                  CALL RINIT(IRMD*LMPOTD*NATREF,R2NSEQ(1,1,1,ISPIN))
                  CALL RINIT(NATREF,SUMNS)
                  CALL RINIT((LMAXD+1)*NATREF,VASUM)
                  CALL RINIT((LMAXD+1)*NATREF,ESPV(0,1,ISPIN))

                  IF (ISPIN.EQ.2 .AND. NHSPIN.EQ.1) THEN

                      IF (NOT(FXDRKEY)) THEN
                          CLOSE (IHANDLE)
                          OPEN (IHANDLE,FILE='green',FORM='UNFORMATTED',
     +                     access='stream')
                      ELSE
c                          CALL FXDRCLS(IHANDLE)
c                          CALL FXDROPN('green','DECODE',IHANDLE)
                      END IF

                  END IF
ccccccccccccccccc
                  IF (NOT(FXDRKEY)) THEN
                  write(6,*) ' ihandle  (2)',ihandle
                    READ(IHANDLE) IELAST,nhspin,EFMTZ,EFMTZ
                   write(6,*) ielast,nhspin,efmtz,efmtz
                  ELSE
c                      CALL FXDRINT(IHANDLE,IELAST,1)
c                      CALL FXDRINT(IHANDLE,NHSPIN,1)
c                      CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                      CALL FXDRDBL(IHANDLE,EFMTZ,1)
                  END IF
ccccccccccccccccc
                  write (6,fmt=*) 'NERGIES ',ielast,nhspin,efmtz

                  IF (ISPIN.EQ.1) THEN
                      REWIND IOWFCT + 1
                      REWIND IOWFCT + 2
                  END IF
                  IF (ISPIN.NE.2 .OR. NHSPIN.NE.1) THEN

                      WRITE (6,FMT='(a14,1p,d10.3)') ' fermi-energy ',
     +                  EFMTZ
                      NLST(ISPIN) = IELAST
                      IF (NSPIN.EQ.2 .AND. NHSPIN.EQ.
     +                    1) NLST(NSPIN) = IELAST
                      IF (IELAST.GT.IEMXD) THEN
                          STOP 19

                      END IF
                  END IF
c
                  DO 190 IE = 1,IELAST
                      DO 180 NP = 1,NREP
c
c---> read struc. host green's fct from disc
c
                          IF (NOT(FXDRKEY)) THEN
                              CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,
     +                                    IHANDLE)
                          ELSE
C                              CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),
C    +                                       GTEMP,IHANDLE)
                          END IF

                WRITE(6,*) IE,E,DF
                          IEN(IE,ISPIN) = REAL(E)
                          IENI(IE,ISPIN) = DIMAG(E)
                          IF (NSPIN.EQ.2 .AND. NHSPIN.EQ.1) THEN
                              IEN(IE,NSPIN) = REAL(E)
                              IENI(IE,NSPIN) = DIMAG(E)
                          END IF

                          DF = DF/REAL(NSPIN)
                          IF (NP.EQ.1) THEN
                              CALL CPLXWF(IE,E,KVREL,NPTREF,ISPIN,
     +                                    NATREF,ALPHA,NSPIN,1,IPANEQ,
     +                                    IRCUEQ,MASS,C,DROREQ,RSEQ,S,
     +                                    PZ,FZ,QZ,SZ,TMAT,VSPSME,RWS,
     +                                    RWSM1,DRDIEQ,REQ,Z,A,B,IRWSEQ,
     +                                    IRT,LMAXP1,IRMD)
CO
CO------> SET UP HOST ALPHA-MATRIX FOR LLOYD'S FORMULA
CO

                              IF ((KLATR.EQ.2) .OR. (KLATR.EQ.3)) THEN
                                  DO I1 = 1,NATREF
                                      DO L1 = 1,LMAX + 1
                                          ALPHAHOST(L1,I1,IE,
     +                                      ISPIN) = ALPHA(L1,I1)
                                      END DO
                                  END DO
                              END IF
c
c
c---> calculate non spherical wavefunctions and t - matrix
c
                              IF (KVREL.GT.1) THEN
                                  NSRA = 2

                              ELSE

                                  NSRA = 1
                              END IF

                              REWIND IOWFCT

c     wfdisk=true : write WFs to the disk (this should be changed
c     in the routines fpimpu, tmatns, and rhons)
c                  WFDISK = .true.
                              WFDISK = .false.

                              IF (INS.LE.2) THEN
                                  NONSPH = .TRUE.
                              ELSE
                                  NONSPH = .FALSE.
                              END IF
                              DO 130 I1 = 1,NATREF
                                  IPOT = NSPIN* (I1-1) + ISPIN
                                  I1DISK = I1
                                  IF (WFDISK) I1DISK = 1

                                  CALL TMATNS(AR(1,1,I1),CR(1,1,I1),
     +                                        DRDIEQ(1,I1),E,ICST,
     +                                        IOWFCT,LMAX,PZ(1,1,I1),
     +                                        QZ(1,1,I1),FZ(1,1,I1),
     +                                        SZ(1,1,I1),PNS(1,1,IRMIND,
     +                                        1,I1DISK),QNS(1,1,IRMIND,
     +                                        1,I1DISK),TMATLL(1,1,I1),
     +                                        VINS(IRMIND,1,IPOT),
     +                                        VISP(1,IPOT),
     +                                        VSPSME(1,IPOT),IRWSEQ(I1),
     +                                        IPANEQ(I1),IRCUEQ(0,I1),
     +                                        IRMIEQ(I1),NSRA,C,CLEB,
     +                                        ICLEB,IEND,LOFLM,
     +                                        TMAT(1,I1),NONSPH,LKONV)

c                   CALL SYTMAT(LMAX,TMATLL(1,1,I1),IDND(1,I1),ND,YR,
c    +                          WTYR,RIJ,IJEND)

                                  CALL TMWRIT(TMATLL(1,1,I1),IOWFCT+1)



c
c   In the case of lattice ralaxation t-host is transformed, here it is
c   saved for later use (look at RLXLAT)
c
                                  IF (KLATR.GT.0) CALL TMWRIT(TMATLL(1,
     +                                1,I1),IOWFCT+2)
c
CASYM-----------------> SAVE HOST MATRICES IN CASE OF SYMMETRY
                                  IF (KSYMMAT.GT.0) THEN
                                      CALL TMWRIT(AR(1,1,I1),92)
                                  END IF
CASYM----------------->END ASYM
c
                                  DET = CONE
                                  IF (KTE.EQ.1 .AND. INS.GE.3) THEN
c
c---> solve the system of linear equations and calculate determinant
c
                                      CALL ZCOPY(LMMAXD**2,AR(1,1,I1),1,
     +                                           ADET,1)
                                      CALL ZGETRF(LMMAXD,LMMAXD,ADET,
     +                                            LMMAXD,IPVT,INFO)
                                      DO 120 I = 1,LMMAXD
                                          IF (IPVT(I).NE.I) DET = -DET
                                          DET = ADET(I,I)*DET
  120                                 CONTINUE
                                  END IF

                                  DETALF(I1,IE,ISPIN) = DET
                                  DETALFH(I1,IE,ISPIN) = DET
                                  DETH = DET
  130                         CONTINUE

                          END IF

                          IF (KSYMM.GT.0) THEN

                              CALL BACSYM(DTB,GMAT,LMAXSQ,NREP,NP)
                          ELSE

                              CALL BACNOSYM(DTB,GMAT,LMAXSQ,NREP,NP)

                          END IF
                          IF (NP.EQ.NREP) THEN
c
c The mapping of the Green's function is done
c for the host atoms
c The first impurity always disturbs the first reference potential
c
                   indx = 1
                   indx_ref =1
                   iref_last=1

                  DO I1=2,NATPER
                     if (indx_ref .eq. natref) exit
                    indx = indx + 1
                     if (iref(i1) .ne. iref_last) then
                        indx_ref = indx_ref + 1
                        iref_last = iref(i1)
                        dtb(:,:,indx_ref) = dtb(:,:,indx)
                     end if
                  end do
c
c------ end mapping of greensfunction for host calculation----
c
c$$$c
c$$$c THIS is for more atoms per unit cell.
c$$$c The mapping of the Green's function is not done correctly
c$$$c for the host atoms (A. Settels)
c$$$CO ----> THIS WORKS ONLY FOR SYMMETRY
c$$$CO ----> THEN INDX HAS TO BE THE NUMBER OF THE
c$$$CO ----> REPRESENTATIVE ATOM
c$$$c
C                              IF (NATREF.GT.1) THEN
C                                  DO I1 = 1,NATREF
C                                      IF (I1.EQ.1) THEN
C                                          INDX = 1
C                                      ELSE IF (I1.EQ.2) THEN
C                                          INDX = 4
C                                      ELSE IF (I1.EQ.3) THEN
C                                          INDX = 6
C                                      ELSE IF (I1.EQ.4) THEN
C                                          INDX = 2
C                                      END IF
C                                      DO L1 = 1,LMAXSQ
C                                          DO L2 = 1,LMAXSQ
C                                              REFHOST(L1,L2,
C     +                                          I1) = DTB(L1,L2,INDX)
C                                          END DO
C                                      END DO
C                                  END DO
C
C                                  DO I1 = 1,NATREF
C                                      DO L1 = 1,LMAXSQ
C                                          DO L2 = 1,LMAXSQ
C                                              DTB(L1,L2,I1)
C     +                                          = REFHOST(L1,L2,I1)
C                                          END DO
C                                      END DO
c$$$c                      WRITE (*,*) DTB(1,1,I1)
C                                  END DO
c$$$c
c$$$c More atoms per unit cell
c$$$c
C                              END IF
c
c---> set up of arrays for the hyperfine field calculation
c
                              IF (KHYPO.NE.0) THEN
                                  DO 140 L = 1,LMAXP1
                                      EKL(L) = EK*REAL(L+L-1)
  140                             CONTINUE
                                  DO 170 L = 1,LMAXP1
                                      DO 160 M = 1,NATREF
                                          GJM(L,M) = CZERO
                                          LM1 = LMS(L)
                                          LM2 = LME(L)
CASYM------------------>This is for host atoms !!!!
                                          EKLASYM(L,M) = EKL(L)
                                          DO 150 LM = LM1,LM2
                                              GJM(L,M) = GJM(L,M) +
     +                                                   DTB(LM,LM,M)
  150                                     CONTINUE
                                          GJHOST(L,ISPIN,IE,M) = GJM(L,
     +                                      M)
  160                                 CONTINUE
  170                             CONTINUE
                                  CALL HYPISO(IPF,1,KVREL,KCOR,NATREF,
     +                                        10,IRWSEQ,NATREF,ISPIN,
     +                                        NSPIN,IELAST,LMAXP1,CFG,
     +                                        REQ,DRDIEQ,RHOC,IE,DF,EKL,
     +                                        GJM,EKLASYM)

                              END IF



                              IF (INS.GE.3) THEN
c
c---> in case of non-spherical host potential
c
                                  CALL RHONS(DEN,DF,DRDIEQ,DTB,E,IE,
     +                                       IELAST,IOWFCT,ISPIN,IRMIEQ,
     +                                       IRWSEQ,LMAX,NATREF,NSPIN,1,
     +                                       NATREF,R2NSEQ,IPANEQ,
     +                                       IRCUEQ,THEEQM,NTCELL,IFUNM,
     +                                       LMSP,KVREL,QNS,PNS,AR,CR,C,
     +                                       PZ,FZ,QZ,SZ,CLEB(1,1),
     +                                       ICLEB,JEND,IEND,KESYM,
     +                                       DTB1(1,1,1,ISPIN),KSYMMAT)
                              ELSE
c
c---> in case of spherical host potential
c
                                  CALL RHOLM(DEN,DF,DTB,E,IE,IELAST,IPF,
     +                                       ISPIN,IRWSEQ,KVREL,LMAX,
     +                                       LMAXSQ,NSPIN,NATREF,1,
     +                                       NATREF,R2NSEQ,DRDIEQ,
     +                                       IPANEQ,IRCUEQ,THEEQM,
     +                                       NTCELL,IFUNM,LMSP,SUMNS,
     +                                       VASUM,C,PZ,FZ,QZ,SZ,
     +                                       CLEB(1,1),ICLEB,IEND,JEND,
     +                                       KSYMMAT,KESYM,
     +                                       DTB1(1,1,1,ISPIN))

                              END IF
c
                              IF (KTE.EQ.1) CALL ESINGPV(DEN,DF,E,ESPV,
     +                                           IE,IELAST,LMAX,ISPIN,
     +                                           IREF,NSHELL,1,NATREF,
     +                                           NATREF,DSELOC)
c

                          END IF

  180                 CONTINUE
c             nrep loop ends

  190             CONTINUE
c           energy loop ends

  200         CONTINUE

              IF (NOT(FXDRKEY)) THEN
                  CLOSE (IHANDLE)
              ELSE
c                  CALL FXDRCLS(IHANDLE)
              END IF

c         spin loop ends


              DO ISPIN = 1,NSPIN
                  DO I1 = 1,NATREF
                      DO L = 1,LMPOTD
c
                          IPOT = NSPIN* (I1-1) + ISPIN
                          DO I = 1,IMT(I1)
                              RHO2NS(I,L,I1,ISPIN) = R2NSEQ(I,L,I1,
     +                          ISPIN)
                              RHOCK(I,IPOT) = RHOC(I,IPOT)
                          END DO
                          IF (KSHAPE.GT.0 .AND.
     +                        (IMT(I1).NE.IRCUEQ(1,I1)))
     +                        STOP 'IMT .ne. IMTEQ, Stop in Main 1'
c     lsmear=2,4: host has equal mesh
c     lsmear=1,3: host has original mesh with kinks
                          IF ((LSMEAR.EQ.4) .OR. (LSMEAR.EQ.2)) THEN
                              ATIM = 1.0d35
                              II = IMT(I1) + 1
                              NPTS = IRMD - IMT(I1)
                              IF (NPTS.NE.IRID) STOP
     +                            ' npts .ne. irid in msgf'
c
                              CALL SPLINE(REQ(II,I1),
     +                                    R2NSEQ(II,L,I1,ISPIN),NPTS,
     +                                    ATIM,ATIM,R2NS2D)
                              DO I = IMT(I1) + 1,IRWS(I1)
                                  X = R(I,I1)
                                  CALL SPLINT(REQ(II,I1),
     +                                        R2NSEQ(II,L,I1,ISPIN),
     +                                        R2NS2D,NPTS,X,Y)
                                  RHO2NS(I,L,I1,ISPIN) = Y
                              END DO
c
                              CALL SPLINE(REQ(II,I1),RHOC(II,IPOT),NPTS,
     +                                    ATIM,ATIM,R2NS2D)
                              DO I = IMT(I1) + 1,IRWS(I1)
                                  X = R(I,I1)
                                  CALL SPLINT(REQ(II,I1),RHOC(II,IPOT),
     +                                        R2NS2D,NPTS,X,Y)
                                  RHOCK(I,IPOT) = Y
                              END DO
c
                          ELSE
                              IF (IRWS(I1).NE.IRWSEQ(I1))
     +                            STOP 'MAIN: irws .ne. irwseq'
                              DO I = IMT(I1) + 1,IRWS(I1)
                                  RHO2NS(I,L,I1,ISPIN) = R2NSEQ(I,L,I1,
     +                              ISPIN)
                                  RHOCK(I,IPOT) = RHOC(I,IPOT)
                              END DO
                          END IF
c
                      END DO
                  END DO
              END DO

              IF (NSPIN.EQ.2) THEN
                  IDIM = IRMKD*LMPOTD*NATREF
                  CALL DCOPY(IDIM,RHO2NS(1,1,1,NSPIN),1,WORK,1)
                  CALL DAXPY(IDIM,ONEM,RHO2NS(1,1,1,1),1,
     +                       RHO2NS(1,1,1,NSPIN),1)
                  CALL DAXPY(IDIM,ONE,WORK,1,RHO2NS(1,1,1,1),1)

                  IDIM = IRMD*LMPOTD*NATREF
                  CALL DCOPY(IDIM,R2NSEQ(1,1,1,NSPIN),1,WORKEQ,1)
                  CALL DAXPY(IDIM,ONEM,R2NSEQ(1,1,1,1),1,
     +                       R2NSEQ(1,1,1,NSPIN),1)
                  CALL DAXPY(IDIM,ONE,WORKEQ,1,R2NSEQ(1,1,1,1),1)
              END IF
c
              ETIME = DCLOCK()
              WRITE (6,FMT=*) 'time for unperturbed atoms',ETIME - STIME
c
c

c     flushes the buffered output
              CALL FLUSH(6)
c

          END IF
c       host atom(s) only in the first iteration  -IF ends here

          STIME = DCLOCK()
c                     ----------------------------------

          IF ((KLATR.EQ.0) .OR. (ITC.EQ.1)) THEN

              IF (NOT(FXDRKEY)) THEN
                  OPEN (IHANDLE,FILE='green',FORM='UNFORMATTED',
     +            access='stream')
              ELSE
c                  CALL FXDROPN('green','DECODE',IHANDLE)
              END IF

          END IF
          IF ((KLATR.GT.0) .AND. (ITC.NE.1)) THEN

              IF (NOT(FXDRKEY)) THEN
                  OPEN (IGHWRIT,FILE='green_TRANS',FORM='UNFORMATTED',
     +             access='stream')
              ELSE
c                  CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
              END IF

              IHANDLE = IGHWRIT

          END IF

          IF ((ITC.EQ.1) .AND. (KLATR.GT.0)) THEN
C---------RELAXATIONS--------------------------------------------------
c
              STIME2 = DCLOCK()
c
              DNEIHO = 0.0E0
C------> LOOP OVER SPINS
              DO ISPIN = 1,NSPIN
C
                  IF (ISPIN.EQ.2 .AND. NHSPIN.EQ.1) THEN

                      IF (NOT(FXDRKEY)) THEN
                          CLOSE (IHANDLE)
                          OPEN (IHANDLE,FILE='green',FORM='UNFORMATTED',
     +                     access='stream')
                      ELSE
c                          CALL FXDRCLS(IHANDLE)
c                          CALL FXDROPN('green','DECODE',IHANDLE)
                      END IF

                      IF (NOT(FXDRKEY)) THEN
                          IF (KGV.EQ.2) then
                             OPEN (IGHWRIT,FILE='green_TRANS',
     +                         FORM='UNFORMATTED',access='stream')
                          END IF
                          CLOSE (IGHWRIT)
                          OPEN (IGHWRIT,FILE='green_TRANS',
     +                         FORM='UNFORMATTED',access='stream')
                      ELSE
                          IF (KGV.EQ.2) then
c                          CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
                          END IF
c                          CALL FXDRCLS(IGHWRIT)
c                          CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
                      END IF

                  END IF
c

                  IF (ISPIN.EQ.1) THEN
                      REWIND IOWFCT + 1
                      REWIND IOWFCT + 2
                      IF ((KLATR.GT.0) .AND. (ITC.GT.1)) THEN
                          REWIND IOTMAT
                      END IF
                  END IF
                  IF (NOT(FXDRKEY)) THEN
                      READ (IHANDLE) IELAST,NHSPIN,EFMTZ,EFMTZ
                  ELSE
c                      CALL FXDRINT(IHANDLE,IELAST,1)
c                      CALL FXDRINT(IHANDLE,NHSPIN,1)
c                      CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                      CALL FXDRDBL(IHANDLE,EFMTZ,1)
                  END IF
C-----> Loop over energies and representations
                  DO IE = 1,IELAST
                      DO NP = 1,NREP

                          IF (NOT(FXDRKEY)) THEN
                              CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,
     +                                    IHANDLE)
                          ELSE
c                              CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),
c     +                                       GTEMP,IHANDLE)
                          END IF
                          DFLLOYD = DF/REAL(NSPIN)
                          IF (KLATR.EQ.1) THEN
CO
CO
CO----> KLATR == 1 -----------------------------------------------
CO
CO--------------------> THIS METHOD IS DONE
CO                      _ ~    ~     ~  ~
CO                      G=G0 + G0 (T-T0)G0
CO
CO---------> DO THIS U-TRAFO ONLY IN THE FIRST ITERATION
CO
CO---------> TRANSFORM G-MATRIX AND T-MATRIX

                              CALL RLXLAT(GMAT,TMATREL,E,LMAX,NATPER,
     +                                    NATREF,NATYP,DRM,NDIM(NP),NP,
     +                                    IOWFCT+2,INS,IRMIN,IREF,NDIM,
     +                                    NREP,KSYMM,KLATR,NOTREAD,
     +                                    TMATHOST)
CO---------> SAVE NEW GREENS FUNCTION IN FXDR - FORMAT
                              IF ((IE.EQ.1) .AND. (NP.EQ.1)) THEN
                                  IF (ISPIN.EQ.1) THEN
                                      IF (NOT(FXDRKEY)) THEN
                                          OPEN (IGHWRIT,
     +                                         FILE='green_TRANS',
     +                                         FORM='UNFORMATTED',
     +                                         access='stream')
                                      ELSE
c                                          CALL FXDROPN('green_TRANS',
c     +                                                 'ENCODE',IGHWRIT)
                                      END IF
                                  END IF
                                  IF (NOT(FXDRKEY)) THEN
                                      WRITE (IGHWRIT) IELAST,NHSPIN,
     +                                  EFMTZ
                                  ELSE
c                                      CALL FXDRINT(IGHWRIT,IELAST,1)
c                                      CALL FXDRINT(IGHWRIT,NHSPIN,1)
c                                      CALL FXDRDBL(IGHWRIT,EFMTZ,1)
c                                      CALL FXDRDBL(IGHWRIT,EFMTZ,1)
                                  END IF
                              END IF
                              IF (NOT(FXDRKEY)) THEN
                                  CALL WRITGTR(E,EK,DF,GMAT,NDIM(NP),
     +                                         IGHWRIT)
                              ELSE
c                                 CALL WRITGTR_FX(E,EK,DF,GMAT,NDIM(NP),
c     +                                            IGHWRIT)
                              END IF
                              IF (NP.EQ.1) THEN
                                  DO I1 = 1,NATYP
                                      CALL TMWRIT(TMATREL(1,1,I1),
     +                                            IOTMAT)
                                  END DO
                              END IF
CO
CO
CO----> END KLATR == 1 -------------------------------------------
CO
CO
                          ELSE IF (KLATR.EQ.2) THEN
CO
CO----> KLATR == 2     -------------------------------------------
CO
CO----> RIGID SHIFT OF POTENTIALS
CO
CO                                        _
CO----> ONE NON SELFCONSISTENT STEP: GET  G0
CO
CO------------------> THIS METHOD IS DONE
CO                     _  ~  ~     ~  _
CO                     G0=G0+G0(T0-T0)G0
CO
CO----> ATTENTION:TMATHOST IS RELAXED T0 !!! TMATHOST:T0^~
CO                TMATREL IS ONLY HOST - T MATRIX !!!!: T0
CO----> THIS IS DONE FOR THE IMPURITY-CALCULATION
CO
                              NEND = NATYP - ICUT
CO
CO----> TRANSFORM G-MATRIX AND T-MATRIX
CO
                              CALL RLXLAT(GMAT,TMATHOST,E,LMAX,NATPER,
     +                                    NATREF,NATYP,DRM,NDIM(NP),NP,
     +                                    IOWFCT+2,INS,IRMIN,IREF,NDIM,
     +                                    NREP,KSYMM,KLATR,NOTREAD,
     +                                    TMATREL)
CO
CO-----> COPY HOST T-MATRIX ---- TOMTOR T0 Minus T0 Relaxed
CO
                              DO I1 = 1,NATYP
                                  DO L1 = 1,LMAXSQ
                                      DO L2 = 1,LMAXSQ
                                          TOMTOR(L1,L2,I1) = TMATREL(L1,
     +                                      L2,I1)
                                      END DO
                                  END DO
                              END DO
CO
CO-----> CALCULATE T0-T0^~
CO
                              NOKTE = 3
CKGV--------------> DONE ONLY IF GTRANS IS NOT READ IN
                              IF (KGV.NE.2) THEN
                                  CALL SYTN0REL(AR,ADET,
     +                                          DETALFH(1,IE,ISPIN),INS,
     +                                          IOWFCT+1,IREF,NOKTE,
     +                                          LMAX,NATPER,NATREF,
     +                                          NATPS,NEND,TOMTOR,
     +                                          NSHELL,TMATHOST)

CO-----> DYSON STEP
                                  IF (KSYMM.GT.0) THEN
                                      CALL SYTN1(DTMTRX,INS,LMAX,NATREF,
     +                                           NATPS,NEND,TOMTOR,
     +                                           IRMIN,NP,NDIM,NREP,
     +                                           .TRUE.)
                                  ELSE
c----------------->this is for no symmetry
                                      CALL SYNOTN1(DTMTRX,INS,LMAX,
     +                                             NATREF,NATPS,NEND,
     +                                             TOMTOR,IRMIN,NP,NDIM,
     +                                             NREP,.TRUE.,NATPER)
                                  END IF
CO--------> DYSON STEP
                                  CALL GDYSON(DTMTRX,GMAT,DETH,NDIM(NP),
     +                                        .FALSE.)
CO--------> COPY ALPHA-MATRIX (IN CASE OF HOST) FOR LLOYD'S FORMULA
                              END IF
CO---------------->END KGV
                              IF (KFSR.EQ.1) THEN
                                  DO I1 = 1,NATREF
                                      DO L1 = 1,LMAX + 1
                                          ALPHA1(L1,I1) = ALPHAHOST(L1,
     +                                      I1,IE,ISPIN)
                                      END DO
                                  END DO
                                  DO I1 = NATPS,NEND
                                      IHOST = IREF(I1-NATREF)
                                      DO L1 = 1,LMAX + 1
                                          ALPHA1(L1,I1) = ALPHAHOST(L1,
     +                                      IHOST,IE,ISPIN)
                                      END DO
                                      DETALFH(I1,IE,ISPIN)
     +                                  = DETALFH(IHOST,IE,ISPIN)
                                  END DO
CO---------> CALCULATE LLOYDS FORMULA FOR NEW REFERENCE SYSTEM
                                  CALL LLOYD(ALPHA1,DETALFH,IREF,NDG,
     +                                       NSHELL,TEXTS,DETH,DFLLOYD,
     +                                       E,DNEIHO,IE,IELAST,ISPIN,
     +                                       NATPS,NATREF,NEND,NP,NREP,
     +                                       NSPIN,LCORE,NCORE,KSYMMAT,
     +                                       KESYM,DTB1(1,1,1,ISPIN))
                              END IF
CO
CO-----> END OF DYSON STEP
CO                                           _
CO----> SAVE GREEN FUNCTION AFTER DYSON STEP G0
CO
CKGV-------------->DONE ONLY IF GTRANS IS NOT READ IN
                              IF (KGV.NE.2) THEN
                                  IF ((IE.EQ.1) .AND. (NP.EQ.1)) THEN
                                      IF (ISPIN.EQ.1) THEN
                                          IF (NOT(FXDRKEY)) THEN
                                              OPEN (IGHWRIT,
     +                                             FILE='green_TRANS',
     +                                             FORM='UNFORMATTED',
     +                                             access='stream')
                                          ELSE
c                                              CALL FXDROPN('green_TRANS'
c     +                                             ,'ENCODE',IGHWRIT)
                                          END IF
                                      END IF
                                      IF (NOT(FXDRKEY)) THEN
                                          WRITE (IGHWRIT) IELAST,
     +                                           NHSPIN, EFMTZ
                                      ELSE
c                                         CALL FXDRINT(IGHWRIT,IELAST,1)
c                                         CALL FXDRINT(IGHWRIT,NHSPIN,1)
c                                         CALL FXDRDBL(IGHWRIT,EFMTZ,1)
c                                         CALL FXDRDBL(IGHWRIT,EFMTZ,1)
                                      END IF
                                  END IF
CO-----> SAVE GREEN's FUNCTION
                                  IF (NOT(FXDRKEY)) THEN
                                      CALL WRITGTR(E,EK,DF,GMAT,
     +                                             NDIM(NP),IGHWRIT)
                                  ELSE
c                                      CALL WRITGTR_FX(E,EK,DF,GMAT,
c     +                                                NDIM(NP),IGHWRIT)
                                  END IF
                              END IF
CKGV------------->END KGV
                              IF (NP.EQ.1) THEN
CO-----> SAVE T-MATRIX FOR FIRST IRREDUCIBLE REPRESENTATION OF G
                                  DO I1 = 1,NATYP
                                      CALL TMWRIT(TMATREL(1,1,I1),
     +                                            IOTMAT)
                                  END DO
                              END IF
CO
CO
CO----> END KLATR == 2 -------------------------------------------
CO
CO
                          ELSE IF (KLATR.EQ.3) THEN
CO
CO----> KLATR == 3     -------------------------------------------
CO
CO----> G-VOID-METHOD
CO
CO
CO----> ONE NON SELFCONSISTENT STEP: GET  Gv
CO
CO------------------> THIS METHOD IS DONE
CO                                      ~        T
CO              Gv=G0+G0 (0-t0) Gv THEN Gv=U Gv U
CO
CO
                              NEND = NATYP - ICUT
CO
CO----> TRANSFORM G-MATRIX AND T-MATRIX
CO
CO----> READ ONLY T-MATRIX:
                              CALL RLXLAT(GMAT,TMATHOST,E,LMAX,NATPER,
     +                                    NATREF,NATYP,DRM,NDIM(NP),NP,
     +                                    IOWFCT+2,INS,IRMIN,IREF,NDIM,
     +                                    NREP,KSYMM,KLATR,NOTREAD,
     +                                    TMATREL)
CO
CO-----> TMATREL IS ONLY HOST MATRIX
CO
CO-----> COPY HOST T-MATRIX ---- TOMTOR 0 Minus T0 = -T0
CO
                              DO I1 = 1,NATYP
                                  DO L1 = 1,LMAXSQ
                                      DO L2 = 1,LMAXSQ
                                          TOMTOR(L1,L2,I1) = (0.0D0,
     +                                      0.0D0)
                                      END DO
                                  END DO
                              END DO

CO
CO-----> CALCULATE T0-T0^~
CO
                              NOKTE = 3
CKGV--------------> DONE ONLY IF GTRANS IS NOT READ IN
                              IF (KGV.NE.2) THEN
                                  CALL SYTN0REL(AR,ADET,
     +                                          DETALFH(1,IE,ISPIN),INS,
     +                                          IOWFCT+1,IREF,NOKTE,
     +                                          LMAX,NATPER,NATREF,
     +                                          NATPS,NEND,TOMTOR,
     +                                          NSHELL,TMATHOST)

CO-----> DYSON STEP
                                  IF (KSYMM.GT.0) THEN
                                      CALL SYTN1(DTMTRX,INS,LMAX,NATREF,
     +                                           NATPS,NEND,TOMTOR,
     +                                           IRMIN,NP,NDIM,NREP,
     +                                           .TRUE.)
                                  ELSE
c----------------->this is for no symmetry
                                      CALL SYNOTN1(DTMTRX,INS,LMAX,
     +                                             NATREF,NATPS,NEND,
     +                                             TOMTOR,IRMIN,NP,NDIM,
     +                                             NREP,.TRUE.,NATPER)
                                  END IF
CO--------> DYSON STEP
                                  CALL GDYSON(DTMTRX,GMAT,DETH,NDIM(NP),
     +                                        .FALSE.)
CO                                      ~
CO--------> TRANSFORM ONLY Gv !!!! GET: Gv
CO
CO                          NOTREAD=7
CO                 CALL RLXLAT(GMAT,TMATHOST,E,LMAX,NATPER,
CO     +                      NATREF,NATYP,DRM,NDIM(NP),NP,
CO     +                      IOWFCT+2,INS,IRMIN,IREF,NDIM,NREP,
CO     +                      KSYMM,KLATR,NOTREAD,TMATREL)
CO
CO--------> COPY ALPHA-MATRIX (IN CASE OF HOST) FOR LLOYD'S FORMULA
CO
                              END IF
CO---------------->END KGV
                              IF (KFSR.EQ.1) THEN
                                  DO I1 = 1,NATREF
                                      DO L1 = 1,LMAX + 1
                                          ALPHA1(L1,I1) = ALPHAHOST(L1,
     +                                      I1,IE,ISPIN)
                                      END DO
                                  END DO
                                  DO I1 = NATPS,NEND
                                      IHOST = IREF(I1-NATREF)
                                      DO L1 = 1,LMAX + 1
                                          ALPHA1(L1,I1) = ALPHAHOST(L1,
     +                                      IHOST,IE,ISPIN)
                                      END DO
                                      DETALFH(I1,IE,ISPIN)
     +                                  = DETALFH(IHOST,IE,ISPIN)
                                  END DO
CO---------> CALCULATE LLOYDS FORMULA FOR NEW REFERENCE SYSTEM
                                  CALL LLOYD(ALPHA1,DETALFH,IREF,NDG,
     +                                       NSHELL,TEXTS,DETH,DFLLOYD,
     +                                       E,DNEIHO,IE,IELAST,ISPIN,
     +                                       NATPS,NATREF,NEND,NP,NREP,
     +                                       NSPIN,LCORE,NCORE,KSYMMAT,
     +                                       KESYM,DTB1(1,1,1,ISPIN))
                              END IF
CO
CO-----> END OF DYSON STEP
CO                                           ~
CO----> SAVE GREEN FUNCTION AFTER DYSON STEP Gv
CO
CKGV-------------->DONE ONLY IF GTRANS IS NOT READ IN
                              IF (KGV.NE.2) THEN
                                  IF ((IE.EQ.1) .AND. (NP.EQ.1)) THEN
                                      IF (ISPIN.EQ.1) THEN
                                          IF (NOT(FXDRKEY)) THEN
                                              OPEN (IGHWRIT,
     +                                             FILE='green_TRANS',
     +                                             FORM='UNFORMATTED',
     +                                             access='stream')
                                          ELSE
c                                             CALL FXDROPN('green_TRANS'
c    +                                             ,'ENCODE',IGHWRIT)
                                          END IF
c                    WRITE (IGHWRIT) IELAST,NHSPIN,EFMTZ
                                      END IF
                                      IF (NOT(FXDRKEY)) THEN
                                          WRITE (IGHWRIT) IELAST,NHSPIN,
     +                                                    EFMTZ
                                      ELSE
c                                        CALL FXDRINT(IGHWRIT,IELAST,1)
c                                        CALL FXDRINT(IGHWRIT,NHSPIN,1)
c                                        CALL FXDRDBL(IGHWRIT,EFMTZ,1)
c                                        CALL FXDRDBL(IGHWRIT,EFMTZ,1)
                                      END IF
                                  END IF
CO-----> SAVE GREEN's FUNCTION
                                  IF (NOT(FXDRKEY)) THEN
                                      CALL WRITGTR(E,EK,DF,GMAT,
     +                                             NDIM(NP),IGHWRIT)
                                  ELSE
c                                      CALL WRITGTR_FX(E,EK,DF,GMAT,
c    +                                                NDIM(NP),IGHWRIT)
                                  END IF
                              END IF
CKGV------------->END KGV
CO
CO
CO----> END KLATR == 3 -------------------------------------------
CO
CO
                          END IF
c              klatr -IF ends (klatr = 1 or 2 or 3)

                      END DO
                  END DO
c          energy and nrep loops end


              IF (ISPIN .EQ. NSPIN) THEN

                  REWIND IOTMAT

                  IF (NOT(FXDRKEY)) THEN
                      CLOSE (IHANDLE)
                  ELSE
c                     CALL FXDRCLS(IHANDLE)
                  END IF


                  IF (KGV.NE.2) THEN

                      IF (NOT(FXDRKEY)) THEN
                          CLOSE (IGHWRIT)
                      ELSE
c                         CALL FXDRCLS(IGHWRIT)
                      END IF

                  END IF

                  IHANDLE = IGHWRIT

              END IF

              END DO
c        spin loop ends
c
              ETIME2 = DCLOCK()
              WRITE (6,FMT=*) 'time for relaxation set up ',
     +          ETIME2 - STIME2
c
          END IF
c      itc = 1 and klatr > 0   -IF ends

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
          IF (KGV.EQ.1) THEN
              DO ISPIN = 1,NHSPIN
                  OPEN (84,FILE='LOW_DATA',FORM='UNFORMATTED')
                  READ (84,FMT='(4I5)') NREPLTL,NATLTL,LHIGH,LLOW
                  IF (ISPIN.EQ.1) THEN
                      IF (NOT(FXDRKEY)) THEN
                          OPEN (85,FILE='gtrans_LOW',FORM='UNFORMATTED')
                          OPEN (86,FILE='ghost_LOW',FORM='UNFORMATTED')
                      ELSE
                          IHANDLE2 = 85
                          IHANDLE3 = 86
c                         CALL FXDROPN('gtrans_LOW','ENCODE',IHANDLE2)
c                         CALL FXDROPN('ghost_LOW','ENCODE',IHANDLE3)
                      END IF
C-------> READ HOST
                      IHANDLE = 80
                      IGHWRIT = 81
                      IF (NOT(FXDRKEY)) THEN
                          OPEN (IHANDLE,FILE='green',
     +                          FORM='UNFORMATTED',access='stream')
                          OPEN (IGHWRIT,FILE='green_TRANS',
     +                          FORM='UNFORMATTED',access='stream')
                          READ (IHANDLE) IELAST,NHSPIN,EFMTZ
                          READ (IGHWRIT) IELAST,NHSPIN,EFMTZ
                      ELSE
c                         CALL FXDROPN('green','DECODE',IHANDLE)
c                         CALL FXDRINT(IHANDLE,IELAST,1)
c                         CALL FXDRINT(IHANDLE,NHSPIN,1)
c                         CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                         CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                         CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
c                         CALL FXDRINT(IGHWRIT,IELAST,1)
c                         CALL FXDRINT(IGHWRIT,NHSPIN,1)
c                         CALL FXDRDBL(IGHWRIT,EFMTZ,1)
c                         CALL FXDRDBL(IGHWRIT,EFMTZ,1)
                      END IF
                  END IF
                  IF (NOT(FXDRKEY)) THEN
                      WRITE (85) IELAST,NHSPIN,EFMTZ
                      WRITE (86) IELAST,NHSPIN,EFMTZ
                  ELSE
c                     CALL FXDRINT(IHANDLE2,IELAST,1)
c                     CALL FXDRINT(IHANDLE2,NHSPIN,1)
c                     CALL FXDRDBL(IHANDLE2,EFMTZ,1)
c                     CALL FXDRDBL(IHANDLE2,EFMTZ,1)
c                     CALL FXDRINT(IHANDLE3,IELAST,1)
c                     CALL FXDRINT(IHANDLE3,NHSPIN,1)
c                     CALL FXDRDBL(IHANDLE3,EFMTZ,1)
c                     CALL FXDRDBL(IHANDLE3,EFMTZ,1)
                  END IF
                  DO IE = 1,IELAST
                      REWIND (84)
                      READ (84,FMT='(4I5)') NREPLTL,NATLTL,LHIGH,LLOW
                      WRITE (*,FMT='(4I5)') NREPLTL,NATLTL,LHIGH,LLOW
                      IF (LHIGH.NE.LMAX) STOP 'LHIGH'
                      IF (LLOW.NE.LKONV) STOP 'LKONV - LLOW'
                      IF (NREP.NE.NREPLTL) STOP 'NREP - NREPLTL'
                      DO NP = 1,NREP
                          READ (84,FMT='(4I5)') NDIMHIGH,NDIMLOW
                          WRITE (*,FMT='(4I5)') NDIMHIGH,NDIMLOW
                          IF (NDIMHIGH.NE.NDIM(NP)) THEN
                              WRITE (*,FMT=*) NDIMHIGH,NDIM(NP)
                              STOP 'NSTOP'
                          END IF
                          DO I1 = 1,NDIMHIGH
                              READ (84,FMT='(2I5)') NICHEK,P4TO3(I1)
                              WRITE (*,FMT='(2I5)') NICHEK,P4TO3(I1)
                              IF (NICHEK.NE.I1) STOP 'NICHEK'
                          END DO
                          IF (NOT(FXDRKEY)) THEN
C--------------->READ  L_HIGH-HOST-GREEN-FUNCTION
                              CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,
     +                                    IHANDLE)
C--------------->WRITE L_HIGH-HOST-GREEN-FUNCTION
                              CALL WRITLOWL(E,EK,DF,GMAT,NDIM(NP),86,
     +                                      P4TO3,NDIMLOW)
C--------------->READ  L_LOW-HOST-GREEN-FUNCTION
                              CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,
     +                                    IGHWRIT)
C--------------->WRITE L_LOW-HOST-GREEN-FUNCTION
                              CALL WRITLOWL(E,EK,DF,GMAT,NDIM(NP),85,
     +                                      P4TO3,NDIMLOW)
                          ELSE
C--------------->READ  L_HIGH-HOST-GREEN-FUNCTION
c                             CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),
c    +                                       GTEMP,IHANDLE)
C--------------->WRITE L_HIGH-HOST-GREEN-FUNCTION
c                             CALL WRITLOWL_FX(E,EK,DF,GMAT,NDIM(NP),
c    +                                         IHANDLE3,P4TO3,NDIMLOW)
C--------------->READ  L_LOW-HOST-GREEN-FUNCTION
c                             CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),
c    +                                       GTEMP,IGHWRIT)
C--------------->WRITE L_LOW-HOST-GREEN-FUNCTION
c                             CALL WRITLOWL_FX(E,EK,DF,GMAT,NDIM(NP),
c    +                                         IHANDLE2,P4TO3,NDIMLOW)
                          END IF

                      END DO
                  END DO
                  CLOSE (84)
              END DO
              IF (NOT(FXDRKEY)) THEN
                  CLOSE (IHANDLE)
                  CLOSE (IGHWRIT)
                  CLOSE (85)
                  CLOSE (86)
              ELSE
c                 CALL FXDRCLS(IHANDLE)
c                 CALL FXDRCLS(IGHWRIT)
c                 CALL FXDRCLS(IHANDLE2)
c                 CALL FXDRCLS(IHANDLE3)
              END IF
          END IF

c-----> Spin loop starts
          IF (NSTOP.EQ.1) STOP 'GTRANS'

          IF (KLATR.EQ.0) THEN
              IF (NOT(FXDRKEY)) THEN
                  OPEN(IHANDLE,FILE='green',FORM='UNFORMATTED',
     +             access='stream')
              ELSE
c                 CALL FXDROPN('green','DECODE',IHANDLE)
              END IF
          ELSE
              IF (NOT(FXDRKEY)) THEN
                  OPEN(IGHWRIT,FILE='green_TRANS',FORM='UNFORMATTED',
     +             access='stream')
              ELSE
c                 CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
              END IF
              IHANDLE = IGHWRIT
          END IF

          DO 300 ISPIN = 1,NSPIN
              CALL RINIT(IRMKD*LMPOTD*NATPER,RHO2NS(1,1,NATPS,ISPIN))
              CALL RINIT(IRMD*LMPOTD*NATPER,R2NSEQ(1,1,NATPS,ISPIN))
              CALL RINIT(NATPER,SUMNS(NATPS))
              CALL RINIT((LMAXD+1)*NATPER,VASUM(0,NATPS))
              CALL RINIT((LMAXD+1)*NATPER,ESPV(0,NATPS,ISPIN))
              IF (ISPIN.EQ.2 .AND. NHSPIN.EQ.1) THEN
                  IF (KLATR.EQ.0) THEN
                      IF (NOT(FXDRKEY)) THEN
                          CLOSE (IHANDLE)
                          OPEN (IHANDLE,FILE='green',FORM='UNFORMATTED',
     +                     access='stream')
                      ELSE
c                         CALL FXDRCLS(IHANDLE)
c                         CALL FXDROPN('green','DECODE',IHANDLE)
                      END IF
                  ELSE
                      IF (NOT(FXDRKEY)) THEN
                          CLOSE (IGHWRIT)
                          OPEN (IGHWRIT,FILE='green_TRANS',
     +                         FORM='UNFORMATTED',access='stream')
                      ELSE
c                         CALL FXDRCLS(IGHWRIT)
c                         CALL FXDROPN('green_TRANS','DECODE',IGHWRIT)
                      END IF
                      IHANDLE = IGHWRIT
                  END IF
              END IF
              IF (ISPIN.EQ.1) THEN
                  REWIND IOWFCT + 1
                  REWIND IOWFCT + 2
CASYM------>92 is ar
                  REWIND (92)
                  IF ((KLATR.GT.0) .AND. (ITC.GT.1)) THEN
                      REWIND IOTMAT
                  END IF
              END IF
              IF (NOT(FXDRKEY)) THEN
                  READ (IHANDLE) IELAST,NHSPIN,EFMTZ
              ELSE
c                 CALL FXDRINT(IHANDLE,IELAST,1)
c                 CALL FXDRINT(IHANDLE,NHSPIN,1)
c                 CALL FXDRDBL(IHANDLE,EFMTZ,1)
c                 CALL FXDRDBL(IHANDLE,EFMTZ,1)
              END IF
c------> Energy and nrep loops start
              DO 290 IE = 1,IELAST
                  DO 220 NP = 1,NREP
CO----> IN CASE OF LATTICE RELAXATION READ T-MATRIX AGAIN
                      IF ((KLATR.EQ.1) .OR. (KLATR.EQ.2)) THEN
                          IF (NP.EQ.1) THEN
                              DO I1 = 1,NATYP
                                  CALL TMREAD(TMATREL(1,1,I1),IOTMAT)
                              END DO
                          END IF
                      ELSE IF (KLATR.EQ.3) THEN
                          DO I1 = NATPS,NATYP
                              DO L1 = 1,LMAXSQ
                                  DO L2 = 1,LMAXSQ
                                      TMATREL(L1,L2,I1) = (0.0D0,0.0D0)
                                  END DO
                              END DO
                          END DO
                      END IF
CO----> READ GREEN'S FUNCTION
                      STIME2 = DCLOCK()
                      IF (NOT(FXDRKEY)) THEN
                          CALL REDOGH(E,EK,DF,GMAT,NDIM(NP),GTEMP,
     +                                IHANDLE)
                      ELSE
c                         CALL REDOGH_FX(E,EK,DF,GMAT,NDIM(NP),GTEMP,
c    +                                   IHANDLE)
                      END IF
                      ETIME2 = DCLOCK()
cz
cz  Attention:  occup(NATOMD)  or  occup(NTPERD)  ?
                      IF (KLATR.EQ.10) THEN
                          IF (NP.EQ.1) THEN
                              DO I1 = 1,NATREF
                                  CALL TMREAD(TMATLL(1,1,I1),IOWFCT+2)
                              END DO
                          END IF
                          DO I1 = NATREF + 1,NATYP
                              IF (OCCUP(I1-NATREF).EQ.0) THEN
                                  DO L1 = 1,LMAXSQ
                                      DO L2 = 1,LMAXSQ
                                          TMATREL(L1,L2,I1) = (0.0d0,
     +                                      0.0d0)
                                      END DO
                                  END DO
                              ELSE
                                  IHOST = IREF(I1-NATREF)
                                  DO L1 = 1,LMAXSQ
                                      DO L2 = 1,LMAXSQ
                                          TMATREL(L1,L2,I1) = TMATLL(L1,
     +                                      L2,IHOST)
                                      END DO
                                  END DO
                              END IF
                          END DO

                      END IF
                      TTIMRE = TTIMRE + (ETIME2-STIME2)
                      DF = DF/REAL(NSPIN)
                      IF (NP.EQ.1) THEN
c
                          NEND = NATYP - ICUT
c
c---> apply magnetic field if khfeld equal 1
c
                          IF (KHFELD.EQ.1 .AND. NSPIN.EQ.2) THEN
                              CALL MFIELD(THEEQM,THESME,VISP,VSPSME,
     +                                    IPANEQ,IRCUEQ,NTCELL,HFIELD,
     +                                    RFPI,ISPIN,KSHAPE,NATPS,NEND,
     +                                    NSPIN)
                          END IF
                          STIME2 = DCLOCK()
                          CALL CPLXWF(IE,E,KVREL,NSPPOT,ISPIN,NATREF,
     +                                ALPHA,NSPIN,2,IPANEQ,IRCUEQ,MASS,
     +                                C,DROREQ,RSEQ,S,PZ,FZ,QZ,SZ,TMAT,
     +                                VSPSME,RWS,RWSM1,DRDIEQ,REQ,Z,A,B,
     +                                IRWSEQ,IRT,LMAXP1,IRMD)

                          ETIME2 = DCLOCK()
                          TTIMCP = TTIMCP + (ETIME2-STIME2)
c
c
c---> calculate wavefunctions and t - matrix
c
                          IF (KVREL.GE.1) THEN
                              NSRA = 2

                          ELSE

                              NSRA = 1
                          END IF

                          REWIND IOWFCT
                          IF (INS.EQ.0) THEN
                              NONSPH = .TRUE.
                          ELSE
                              NONSPH = .FALSE.
                          END IF
                          STIME2 = DCLOCK()
                          DO 210 I1 = NATPS,NEND
                              IPOT = NSPIN* (I1-1) + ISPIN
                              I1DISK = I1
                              IF (WFDISK) I1DISK = 1

                              CALL TMATNS(AR(1,1,I1),CR(1,1,I1),
     +                                    DRDIEQ(1,I1),E,ICST,IOWFCT,
     +                                    LMAX,PZ(1,1,I1),QZ(1,1,I1),
     +                                    FZ(1,1,I1),SZ(1,1,I1),
     +                                    PNS(1,1,IRMIND,1,I1DISK),
     +                                    QNS(1,1,IRMIND,1,I1DISK),
     +                                    TMATLL(1,1,I1),
     +                                    VINS(IRMIND,1,IPOT),
     +                                    VISP(1,IPOT),VSPSME(1,IPOT),
     +                                    IRWSEQ(I1),IPANEQ(I1),
     +                                    IRCUEQ(0,I1),IRMIEQ(I1),NSRA,
     +                                    C,CLEB,ICLEB,IEND,LOFLM,
     +                                    TMAT(1,I1),NONSPH,LKONV)

c                 CALL SYTMAT(LMAX,TMATLL(1,1,I1),IDND(1,I1-NATREF),ND,
c    +                        YR,WTYR,RIJ,IJEND)
                              IF (IWTMAT.EQ.1) CALL TMWRIT(TMATLL(1,1,
     +                            I1),IOWFCT+2)
  210                     CONTINUE
                          ETIME2 = DCLOCK()
                          TTIMTM = TTIMTM + (ETIME2-STIME2)
CASYM -------->READ AR-MATRICES in CASE OF SYMMETRY
                          IF (KTE.EQ.1) THEN
                              IF (KSYMMAT.GT.0) THEN
                                  DO I1 = 1,NATREF
                                      CALL TMREAD(AR(1,1,I1),92)
                                  END DO
                              END IF
                          END IF
CASYM -------->END ASYM
                          IF (KLATR.EQ.0) THEN

                              CALL SYTN0(AR,ADET,DETALF(1,IE,ISPIN),INS,
     +                                   IOWFCT+1,IREF,KTE,LMAX,NATPER,
     +                                   NATREF,NATPS,NEND,TMATLL,
     +                                   NSHELL)
                          ELSE
c
c In case of lattice relaxation calculate dt - matrix
c
                              CALL SYTN0REL(AR,ADET,DETALF(1,IE,ISPIN),
     +                                      INS,IOWFCT+1,IREF,KTE,LMAX,
     +                                      NATPER,NATREF,NATPS,NEND,
     +                                      TMATLL,NSHELL,TMATREL)
                          END IF
c
c---> subtract magnetic field
c
                          IF (KHFELD.EQ.1 .AND. NSPIN.EQ.2) THEN
                              CALL MFIELD(THEEQM,THESME,VISP,VSPSME,
     +                                    IPANEQ,IRCUEQ,NTCELL,-HFIELD,
     +                                    RFPI,ISPIN,KSHAPE,NATPS,NEND,
     +                                    NSPIN)
                          END IF
                      END IF
c            above is done only for the first representation

                      STIME2 = DCLOCK()
                      IF (KSYMM.GT.0) THEN
                          CALL SYTN1(DTMTRX,INS,LMAX,NATREF,NATPS,NEND,
     +                               TMATLL,IRMIN,NP,NDIM,NREP,.TRUE.)

                      ELSE
c this is for no symmetry
                          CALL SYNOTN1(DTMTRX,INS,LMAX,NATREF,NATPS,
     +                                 NEND,TMATLL,IRMIN,NP,NDIM,NREP,
     +                                 .TRUE.,NATPER)
                      END IF
                      ETIME2 = DCLOCK()
                      TTIMSY = TTIMSY + (ETIME2-STIME2)

                      STIME2 = DCLOCK()
                      CALL GDYSON(DTMTRX,GMAT,DET,NDIM(NP),.FALSE.)
                      IF (KFSR.EQ.1) THEN
CASYM--------begin---this is for fp-energies (friedel sum rule)--
                          IF ((KSYMMAT.GT.0) .AND. (IE.GE.KESYM)) THEN
                              CALL SYMLLOYD(DETALF,AR,NATREF,NATPER,
     +                                      LMAXSQ,LMAX,NATPS,NEND,NREP,
     +                                      INS,NDIM,IREF,IRMIN,ASARY,
     +                                      ISPIN,IE,NDG,NSHELL)
                              IF (ASARY(NP,ISPIN).NE.NP) DET = (1.0D0,
     +                            0.0D0)
                          END IF
CASYM--------end-------------------------------------------------
                          CALL LLOYD(ALPHA,DETALF,IREF,NDG,NSHELL,TEXTS,
     +                               DET,DF,E,DNEINT,IE,IELAST,ISPIN,
     +                               NATPS,NATREF,NEND,NP,NREP,NSPIN,
     +                               LCORE,NCORE,KSYMMAT,KESYM,
     +                               DTB1(1,1,1,ISPIN))
                      END IF
                      ETIME2 = DCLOCK()
                      TTIMDY = TTIMDY + (ETIME2-STIME2)

                      IF (IRESIST.GT.0 .AND. IE.EQ.IELAST .AND.
     +                    NATREF.EQ.1 .AND. ITC.EQ.1) CALL RESIST(ALAT,
     +                    ISPIN,LMAX,LMAXSQ,N,NATOM,NATYP,NREP,NSPIN,
choshino
c    +                    GMAT,TMATLL,RM,ITITLE,NDIM,NDIM(NP),NSHELL,NP)
     +                    GMAT,TMATLL,RM,ITITLE,NDIM,NDIM(NP),NSHELL)
c---> back symmetrize green's function
c
c
                      STIME2 = DCLOCK()
                      IF (KSYMM.GT.0) THEN
                          IF (KSYMMAT.EQ.0) THEN
                              CALL BACSYM(DTB,GMAT,LMAXSQ,NREP,NP)
                          ELSE IF ((KSYMMAT.GT.0) .AND.
     +                             (IE.GE.KESYM)) THEN
                              NSYMMAT = ASARY(NP,ISPIN)
                              IF (KSYMMAT.EQ. (NREP+1)) NSYMMAT = NP
                              CALL BACASYM(DTB,GMAT,LMAXSQ,NREP,NP,
     +                                     NSYMMAT,SOCCUP,QOCCUP,ISPIN)
                              NSYMMAT = ASARY(NP,ISPIN)
                          ELSE
                              CALL BACSYM(DTB,GMAT,LMAXSQ,NREP,NP)
                          END IF
                      ELSE
c this is for no symmetry
                          CALL BACNOSYM(DTB,GMAT,LMAXSQ,NREP,NP)
                      END IF
                      ETIME2 = DCLOCK()
                      TTIMSY = TTIMSY + (ETIME2-STIME2)
  220             CONTINUE
c           nrep loop ends

c this is new
                  IF (IWTMAT.EQ.1) THEN
                      MR = 0
                      DO 230 M = NATPS,NATYP
                          MR = MR + 1
                          CALL TMWRIT(DTB(1,1,MR),IOWFCT+3)
  230                 CONTINUE
                  ELSE IF (IWTMAT.EQ.2) THEN
                      MR = 0
                      DO 240 M = NATPS,NATYP
                          MR = MR + 1
                          CALL TMREAD(DTB(1,1,MR),IOWFCT+3)
  240                 CONTINUE
                  END IF
ccccc
c
                  IF (KHYPO.NE.0) THEN
                      DO 280 L = 1,LMAXP1
                          EKL(L) = EK*REAL(L+L-1)
                          DO 250 M = 1,NATREF
                              GJM(L,M) = GJHOST(L,ISPIN,IE,M)
  250                     CONTINUE
                          MR = 0
                          DO 270 M = NATPS,NATYP
                              MR = MR + 1
                              GJM(L,M) = CZERO
                              LM1 = LMS(L)
                              LM2 = LME(L)
                              GL = 0.0D0
                              DO 260 LM = LM1,LM2
                                  GJM(L,M) = GJM(L,M) + DTB(LM,LM,MR)
                                  GL = GL + REAL(DTB1(LM,LM,MR,ISPIN))
  260                         CONTINUE
CASYM-------------begin
                              IF ((I1.GT.NATREF) .AND.
     +                            (KSYMMAT.GT.0) .AND.
     +                            (IE.GE.KESYM)) THEN
                                  EKLASYM(L,M) = EK*GL
                              ELSE
                                  EKLASYM(L,M) = EKL(L)
                              END IF
CASYM-------------end
  270                     CONTINUE
  280                 CONTINUE
c
                      CALL HYPISO(IPF,NATREF+1,KVREL,KCOR,NATREF,10,
     +                            IRWSEQ,NATYP,ISPIN,NSPIN,IELAST,
     +                            LMAXP1,CFG,REQ,DRDIEQ,RHOC,IE,DF,EKL,
     +                            GJM,EKLASYM)
CASYM------------------------>
                  END IF
c
                  STIME2 = DCLOCK()
                  IF (INS.GE.1) THEN
c
c---> in case of non-spherical input potential
c
                      CALL RHONS(DEN,DF,DRDIEQ,DTB,E,IE,IELAST,IOWFCT,
     +                           ISPIN,IRMIEQ,IRWSEQ,LMAX,NATREF,NSPIN,
     +                           NATPS,NEND,R2NSEQ,IPANEQ,IRCUEQ,THEEQM,
     +                           NTCELL,IFUNM,LMSP,KVREL,QNS,PNS,AR,CR,
     +                           C,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,JEND,
     +                           IEND,KESYM,DTB1(1,1,1,ISPIN),KSYMMAT)

                  ELSE

c
c---> in case of spherical input potential
c
                      CALL RHOLM(DEN,DF,DTB,E,IE,IELAST,IPF,ISPIN,
     +                           IRWSEQ,KVREL,LMAX,LMAXSQ,NSPIN,NATREF,
     +                           NATPS,NEND,R2NSEQ,DRDIEQ,IPANEQ,IRCUEQ,
     +                           THEEQM,NTCELL,IFUNM,LMSP,SUMNS,VASUM,C,
     +                           PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,IEND,JEND,
     +                           KSYMMAT,KESYM,DTB1(1,1,1,ISPIN))

                  END IF
                  ETIME2 = DCLOCK()
                  TTIMRH = TTIMRH + (ETIME2-STIME2)
c
                  IF (KTE.EQ.1) THEN
                      STIME2 = DCLOCK()
                      CALL ESINGPV(DEN,DF,E,ESPV,IE,IELAST,LMAX,ISPIN,
     +                             IREF,NSHELL,NATPS,NEND,NATREF,DSELOC)
                      ETIME2 = DCLOCK()
                      TTIMES = TTIMES + (ETIME2-STIME2)
                  END IF

  290         CONTINUE
c         energy loop ends

  300     CONTINUE

          IF ((KLATR.EQ.0) .OR. (ITC.EQ.1)) THEN
              IF (NOT(FXDRKEY)) THEN
                  CLOSE (IHANDLE)
              ELSE
c                 CALL FXDRCLS(IHANDLE)
              END IF
          ELSE
              IF (NOT(FXDRKEY)) THEN
                  CLOSE (IGHWRIT)
              ELSE
c                 CALL FXDRCLS(IGHWRIT)
              END IF
          END IF
c       spin loop ends

          DO ISPIN = 1,NSPIN
c          Do i1 = nstart, nend
c     Host part is done already
              DO I1 = NATREF + 1,NEND
                  DO L = 1,LMPOTD
c
                      IPOT = NSPIN* (I1-1) + ISPIN
                      DO I = 1,IMT(I1)
                          RHO2NS(I,L,I1,ISPIN) = R2NSEQ(I,L,I1,ISPIN)
                          RHOCK(I,IPOT) = RHOC(I,IPOT)
                      END DO
                      IF (KSHAPE.GT.0 .AND. (IMT(I1).NE.IRCUEQ(1,I1)))
     +                    STOP 'IMT .ne. IMTEQ, Stop in Main 2'
c     lsmear>0: two shape function sets
                      IF (LSMEAR.GE.1) THEN
                          ATIM = 1.0d35
                          II = IMT(I1) + 1
                          NPTS = IRMD - IMT(I1)
                          IF (NPTS.NE.IRID) STOP
     +                        ' npts .ne. irid in msgf'
c
                          CALL SPLINE(REQ(II,I1),R2NSEQ(II,L,I1,ISPIN),
     +                                NPTS,ATIM,ATIM,R2NS2D)
                          DO I = IMT(I1) + 1,IRWS(I1)
                              X = R(I,I1)
                              CALL SPLINT(REQ(II,I1),
     +                                    R2NSEQ(II,L,I1,ISPIN),R2NS2D,
     +                                    NPTS,X,Y)
                              RHO2NS(I,L,I1,ISPIN) = Y
                          END DO
c
                          CALL SPLINE(REQ(II,I1),RHOC(II,IPOT),NPTS,
     +                                ATIM,ATIM,R2NS2D)
                          DO I = IMT(I1) + 1,IRWS(I1)
                              X = R(I,I1)
                              CALL SPLINT(REQ(II,I1),RHOC(II,IPOT),
     +                                    R2NS2D,NPTS,X,Y)
                              RHOCK(I,IPOT) = Y
                          END DO
c
                      ELSE
c     lsmear=0 case: only one shape function set
                          IF (IRWS(I1).NE.IRWSEQ(I1))
     +                        STOP 'MAIN: irws .ne. irwseq'
                          DO I = IMT(I1) + 1,IRWS(I1)
                              RHO2NS(I,L,I1,ISPIN) = R2NSEQ(I,L,I1,
     +                          ISPIN)
                              RHOCK(I,IPOT) = RHOC(I,IPOT)
                          END DO
                      END IF
c
                  END DO
              END DO
          END DO


          IF (NSPIN.EQ.2) THEN
              IDIM = IRMKD*LMPOTD*NATPER
              CALL DCOPY(IDIM,RHO2NS(1,1,NATPS,NSPIN),1,WORK,1)
              CALL DAXPY(IDIM,ONEM,RHO2NS(1,1,NATPS,1),1,
     +                   RHO2NS(1,1,NATPS,NSPIN),1)
              CALL DAXPY(IDIM,ONE,WORK,1,RHO2NS(1,1,NATPS,1),1)
c
              IDIM = IRMD*LMPOTD*NATPER
              CALL DCOPY(IDIM,R2NSEQ(1,1,NATPS,NSPIN),1,WORKEQ,1)
              CALL DAXPY(IDIM,ONEM,R2NSEQ(1,1,NATPS,1),1,
     +                   R2NSEQ(1,1,NATPS,NSPIN),1)
              CALL DAXPY(IDIM,ONE,WORKEQ,1,R2NSEQ(1,1,NATPS,1),1)
          END IF

          ETIME = DCLOCK()
          WRITE (6,FMT=*) 'time for perturbed atoms  ',ETIME - STIME

          STIME = DCLOCK()

c
c                 ------------------------------------------
c                |  end of do loop over spins and energies  |
c                 __________________________________________
c
c---> only during the first iteration calculate the densities and
c     the potentials for the reference atoms
c
          IF (ITC.EQ.1) CALL VBOUND(EFERMI,EGR,IEN,KTE,NLST,1,NATYP,
     +                              NSPIN,VBC,VBCC,Z)
c
c---> determine last reference atom which is taken into account
c
          IF (ICUT.GT.0) CALL CUTSHELL(ICUT,IREF,IRWSEQ,LPOT,NSPIN,
     +                                 NATREF,NATYP,R2NSEQ)

          IF (ICUT.GT.0) CALL CUTSHELL(ICUT,IREF,IRWS,LPOT,NSPIN,NATREF,
     +                                 NATYP,RHO2NS)
c
c---> determine total charge density expanded in spherical harmonics
c     determine charge density expansion coefficients cdeco by means of
c     a least square fit in case of a shape corrected calculation,
c     Adds the core charge to array rho2ns, does not do anything
c     else to rho2ns.

          IF (LSMEAR.GT.0) THEN
              WRITE (6,FMT=*)
     +          'RHOTOT: using shapes on radial mesh with kinks'
              CALL RHOTOT(ITC,IPF,IREF,NATREF,NATPER,NSHELL,NSPIN,
     +                    NSTART,NEND,RHO2NS,RHOCK,Z,DRDI,IRWS,IRCUT,
     +                    LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,
     +                    IRMKD,IRIKD,2)

c     Calculates intracell-potential part VONS and moments CMOM
              CALL VINTRAS(CMOM,CMINST,LPOT,NSPIN,NSTART,NATYP,RHO2NS,
     +                     VONS,R,DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,
     +                     ILM,IFUNM,IMAXSH,GSH,THETAS,LMSP,IRMKD,IRIKD)
              WRITE (6,FMT=*)
     +          'RHOTOT: using shapes on equally spaced mesh'
          END IF

c
c     Adds the core charge to array rho2ns, does not do anything
c     else to rho2ns.
          CALL RHOTOT(ITC,IPF,IREF,NATREF,NATPER,NSHELL,NSPIN,NSTART,
     +                NEND,R2NSEQ,RHOC,Z,DRDIEQ,IRWSEQ,IRCUEQ,LPOT,NFU,
     +                LLMSP,THEEQM,NTCELL,KSHAPE,IPANEQ,IRMD,IRID,1)

c     Calculates potential VONS and moments CMOM
          CALL VINTRAS(CMOMEQ,CMSTEQ,LPOT,NSPIN,NSTART,NATYP,R2NSEQ,
     +                 VONSEQ,REQ,DRDIEQ,IRWSEQ,IRCUEQ,IPANEQ,KSHAPE,
     +                 NTCELL,ILM,IFUNM,IMAXSH,GSH,THEEQM,LMSP,IRMD,
     +                 IRID)
c
          IF (LSMEAR.EQ.0) THEN
              DO I1 = NSTART,NATYPD
ccc          Do i1 = 1,natypd
                  DO L = 1,LMPOTD
                      CALL DCOPY(IRMD,R2NSEQ(1,L,I1,1),1,
     +                           RHO2NS(1,L,I1,1),1)
                      IF (NSPIN.EQ.2) CALL DCOPY(IRMD,
     +                                     R2NSEQ(1,L,I1,NSPIN),1,
     +                                     RHO2NS(1,L,I1,NSPIN),1)
                      CMOM(L,I1) = CMOMEQ(L,I1)
                      CMINST(L,I1) = CMSTEQ(L,I1)
                  END DO
              END DO
              IDIM = IRMD*LMPOTD*NATYP*NSPIN
              CALL DCOPY(IDIM,VONSEQ,1,VONS,1)
          END IF
c
c---> only for the disturbed atoms the intercell potential is not zero
c
          IF (KLATR.EQ.0) THEN
              IF (SURFACE) THEN
                  IF (LSMEAR.GT.0) THEN
c                     CALL VINTERS_SURFACE(AMAT,BMAT,CMOM,IREF,LPOT,
c    +                                     NATPER,NATREF,NSPIN,NSTART,
c    +                                     NEND,VONS,Z,R,IRWS,IRCUT,
c    +                                     IPAN,KSHAPE,CMINST,IRMKD,2)
                  END IF
c                 CALL VINTERS_SURFACE(AMAT,BMAT,CMOMEQ,IREF,LPOT,
c    +                                 NATPER,NATREF,NSPIN,NSTART,NEND,
c    +                                 VONSEQ,Z,REQ,IRWSEQ,IRCUEQ,
c    +                                 IPANEQ,KSHAPE,CMSTEQ,IRMD,1)
              ELSE
                  IF (NATREF.EQ.1) THEN
c     This routine changes only the potential VONS
                      IF (LSMEAR.GT.0) THEN
                          CALL VINTERS(AMAT,BMAT,CMOM,IREF,LPOT,NATPER,
     +                                 NATREF,NSPIN,NSTART,NEND,VONS,Z,
     +                                 R,IRWS,IRCUT,IPAN,KSHAPE,CMINST,
     +                                 AVMAD,BVMAD,IRMKD,2)
                      END IF
                      CALL VINTERS(AMAT,BMAT,CMOMEQ,IREF,LPOT,NATPER,
     +                             NATREF,NSPIN,NSTART,NEND,VONSEQ,Z,
     +                             REQ,IRWSEQ,IRCUEQ,IPANEQ,KSHAPE,
     +                             CMSTEQ,AVMAD,BVMAD,IRMD,1)
                  ELSE
c this is for more atoms per unit shell
                      IF (LSMEAR.GT.0) THEN
                          CALL VINTERSA(AMAT,BMAT,CMOM,IREF,LPOT,NATPER,
     +                                  NATREF,NSPIN,NSTART,NEND,VONS,Z,
     +                                  R,IRWS,IRCUT,IPAN,KSHAPE,CMINST,
     +                                  AVMAD,BVMAD,IRMKD,2)
                      END IF
                      CALL VINTERSA(AMAT,BMAT,CMOMEQ,IREF,LPOT,NATPER,
     +                              NATREF,NSPIN,NSTART,NEND,VONSEQ,Z,
     +                              REQ,IRWSEQ,IRCUEQ,IPANEQ,KSHAPE,
     +                              CMSTEQ,AVMAD,BVMAD,IRMD,1)
                  END IF
              END IF

          ELSE
              IF (NATREF.EQ.1) THEN
c     This routine changes only the potential VONS
                  IF (LSMEAR.GT.0) THEN
                      CALL VINTREL(AMAT,BMAT,CMOM,IREF,LPOT,NATPER,
     +                             NATREF,NSPIN,NSTART,NEND,VONS,Z,R,
     +                             IRWS,IRCUT,IPAN,KSHAPE,CMINST,AVMAD,
     +                             BVMAD,AMATI,BMATI,DRM,WG,YRG,IRMKD,2)
                  END IF
                  CALL VINTREL(AMAT,BMAT,CMOMEQ,IREF,LPOT,NATPER,NATREF,
     +                         NSPIN,NSTART,NEND,VONSEQ,Z,REQ,IRWSEQ,
     +                         IRCUEQ,IPANEQ,KSHAPE,CMSTEQ,AVMAD,BVMAD,
     +                         AMATI,BMATI,DRM,WG,YRG,IRMD,1)
              ELSE
c     this is for more atoms per unit shell
                  IF (LSMEAR.GT.0) THEN
                      CALL VINTRELA(AMAT,BMAT,CMOM,IREF,LPOT,NATPER,
     +                              NATREF,NSPIN,NSTART,NEND,VONS,Z,R,
     +                              IRWS,IRCUT,IPAN,KSHAPE,CMINST,AVMAD,
     +                              BVMAD,AMATI,BMATI,DRM,WG,YRG,IRMKD,
     +                              2)
                  END IF
                  CALL VINTRELA(AMAT,BMAT,CMOMEQ,IREF,LPOT,NATPER,
     +                          NATREF,NSPIN,NSTART,NEND,VONSEQ,Z,REQ,
     +                          IRWSEQ,IRCUEQ,IPANEQ,KSHAPE,CMSTEQ,
     +                          AVMAD,BVMAD,AMATI,BMATI,DRM,WG,YRG,IRMD,
     +                          1)
              END IF
          END IF

          IF (LSMEAR.EQ.0) THEN
              IDIM = IRMD*LMPOTD*NPOTD
              CALL DCOPY(IDIM,VONSEQ,1,VONS,1)
          END IF


          IF (KHYP.EQ.1.and.NSPIN.EQ.2) THEN
              CALL CMAM(CMAMV,LPOT,NSPIN,NSTART,NATYP,RHO2NS,
     +                  R,DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,
     +                  ILM,IFUNM,IMAXSH,GSH,THETAS,LMSP,IRMKD,
     +                  IRIKD,ONSITE)

              CALL HYPDIP(AMAT,CMAMV,LPOT,NSPIN,NATPS,NATYP,RHO2NS,R,
     +                    DRDI,IRWS,NATREF,NATPER,NSHELL,IOPER,ND,
     +                    YR,WTYR,RIJ,IJEND,RM,ONSITE)


          END IF

          IF (KEFG.EQ.1) THEN
              IF (KSHAPE.EQ.0) THEN
                  CALL EFGRAD(CMOM,LPOT,NSPIN,NATPS,NATYP,RHO2NS,VONS,R,
     +                        DRDI,IRWS,NATREF)

              ELSE


                  CALL EFGRAD(CMOM,LPOT,NSPIN,NATPS,NATYP,RHO2NS,VONS,R,
     +                        DRDI,IMT,NATREF)
              END IF

          END IF

          IF (KF.EQ.1) THEN
              IF (KSHAPE.EQ.0) THEN
                  CALL FORCEH(CMOM,FLM,LPOT,NSPIN,NATPS,NATYP,RHO2NS,
     +                        VONS,R,DRDI,IRWS,Z)

                  CALL FORCE(FLM,FLMC,LPOT,NSPIN,NATPS,NATYP,RHOCK,VONS,
     +                       R,DRDI,IRWS)

              ELSE


                  CALL FORCEH(CMOM,FLM,LPOT,NSPIN,NATPS,NATYP,RHO2NS,
     +                        VONS,R,DRDI,IMT,Z)
                  CALL FORCE(FLM,FLMC,LPOT,NSPIN,NATPS,NATYP,RHOCK,VONS,
     +                       R,DRDI,IMT)
              END IF

          END IF
c
          IF (KTE.EQ.1) THEN
              CALL ESINGPC(ESPC,IEN,NLST,NSPIN,NSTART,NEND,ECORE,LCORE,
     +                     NCORE)
c
c     Calculates energy integral  Int( dr*n(r)*V_eff(r) )
              CALL EPOTINS(EPOTIN,NSPIN,NSTART,NEND,R2NSEQ,VISP,REQ,
     +                     DRDIEQ,INS,IRMIEQ,IRWSEQ,LPOT,VINS,IRCUEQ,
     +                     IPANEQ,Z)
              CALL ECOULOM(CMOM,ECOU,LPOT,NSPIN,NSTART,NEND,RHO2NS,VONS,
     +                     Z,R,DRDI,IRWS,KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,
     +                     IFUNM,ILM,NTCELL,GSH,THETAS,VBC,LMSP)
          END IF
c
C
C
C
C    IGGA=1 GGA-CALCULATION
C    PW91 for IGGA=1
C    PW86 for IGGA=1
C
C
          IF (IGGA.EQ.0) THEN

c     VONS is used only on output, the xc-part of the potential
c     is added to it.
              CALL VXCLM(EXCEQ,KTE,KXC,LPOT,NSPIN,NSTART,NEND,R2NSEQ,
     +                   VONSEQ,REQ,DRDIEQ,IRWSEQ,IRCUEQ,IPANEQ,NTCELL,
     +                   KSHAPE,GSH,ILM,IMAXSH,IFUNM,THEEQM,YR,WTYR,
     +                   IJEND,WORKEQ,LMSP,IRMD,IRID)
              IF (LSMEAR.GT.0) THEN
                  CALL VXCLM(EXC,KTE,KXC,LPOT,NSPIN,NSTART,NEND,RHO2NS,
     +                       VONS,R,DRDI,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,
     +                       GSH,ILM,IMAXSH,IFUNM,THETAS,YR,WTYR,IJEND,
     +                       WORK,LMSP,IRMKD,IRIKD)
              END IF
          ELSE

              IF (IGGA.GT.3) STOP 'GGA'
              CALL VXCGGA(EXCEQ,KTE,KXC,LPOT,NSPIN,NSTART,NEND,R2NSEQ,
     +                    VONSEQ,REQ,DRDIEQ,A,IRWSEQ,IRCUEQ,IPANEQ,
     +                    NTCELL,KSHAPE,GSH,ILM,IMAXSH,IFUNM,THEEQM,YR,
     +                    WTYR,IJEND,WORKEQ,LMSP)
              IF (LSMEAR.GT.0) THEN
                  CALL VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,NSTART,NEND,RHO2NS,
     +                        VONS,R,DRDI,A,IRWS,IRCUT,IPAN,NTCELL,
     +                        KSHAPE,GSH,ILM,IMAXSH,IFUNM,THETAS,YR,
     +                        WTYR,IJEND,WORK,LMSP)
              END IF

          END IF
          IF (LSMEAR.EQ.0) THEN
              IDIM = IRMD*LMPOTD*NATYP*NSPIN
              CALL DCOPY(IDIM,VONSEQ,1,VONS,1)
              IDIM = (LPOTD+1)*NATYP
              CALL DCOPY(IDIM,EXCEQ(0,1),1,EXC(0,1),1)
          END IF
c
          IF (KF.EQ.1) THEN
              IF (KSHAPE.EQ.0) THEN
                  CALL FORCXC(FLM,FLMC,LPOT,NSPIN,NATPS,NATYP,RHOCK,
     +                        VONS,R,ALAT,RM,NSHELL,DRDI,IRWS,NATREF,ND,
     +                        IOPER,F1XYZ)

              ELSE



                  CALL FORCXC(FLM,FLMC,LPOT,NSPIN,NATPS,NATYP,RHOCK,
     +                        VONS,R,ALAT,RM,NSHELL,DRDI,IMT,NATREF,ND,
     +                        IOPER,F1XYZ)
              END IF

          END IF
c
c---> add vbc
c
          DO 330 IS = 1,NSPIN
              DO 320 IATYP = NSTART,NEND
                  IPOT = NSPIN* (IATYP-1) + IS
                  DO 310 IR = 1,IRCEQ(IATYP)
                      VONSEQ(IR,1,IPOT) = VONSEQ(IR,1,IPOT) +
     +                                    VBC(IS)*RFPI
  310             CONTINUE
  320         CONTINUE
  330     CONTINUE

c
c---> fill the outer potentials with the host ones in case of cut off
c
          IF (ICUT.GT.0) THEN
              DO 360 IS = 1,NSPIN
                  DO 350 NP = NEND + 1,NATYP
                      IPOT = NSPIN* (NP-1) + IS
                      IPRF = NSPIN* (IREF(NP-NATREF)-1) + IS
                      DO 340 I = 1,IRWSEQ(NP)
                          VONSEQ(I,1,IPOT) = VISP(I,IPRF)
                          VSPSMO(I,IPOT) = VSPSME(I,IPRF)
  340                 CONTINUE
  350             CONTINUE
  360         CONTINUE

          END IF
c
c---> convolute potential with shape function for next iteration
c
          IF (KSHAPE.NE.0) THEN
              DO 380 IS = 1,NSPIN
                  DO 370 I1 = NATPS,NEND
                      IPOT = NSPIN* (I1-1) + IS
c     now vspsmo contains smeared shape convoluted spherical
c     part of the potential
                      CALL CONVOL(IRCUEQ(1,I1),IRCEQ(I1),NTCELL(I1),
     +                            IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     +                            THEEQM,THESME,Z(I1),RFPI,REQ(1,I1),
     +                            VONSEQ(1,1,IPOT),VSPSMO(1,IPOT),LMSP)
  370             CONTINUE
  380         CONTINUE
          END IF

c     Next makes vspsmo for the host atom (is the same as the
c     smeared input potential)
          IF ((LSMEAR.GE.1) .AND. (ITC.EQ.1)) THEN
c          nend = natyp - icut
              DO I1 = 1,NATREF
                  DO IS = 1,NSPIN
                      IPOT = NSPIN* (I1-1) + IS
                      DO IR = 1,IRCEQ(I1)
                          VSPSMO(IR,IPOT) = VSPSME(IR,IPOT)
                      END DO
                  END DO
              END DO
          END IF


c     Make sure, that  Integral{ V(r)*4pi*r^2 } is equal to
c     Integral{ Vsme(r)*4pi*r^2}. This is needed in order to
c     have the correct MT-zero (VBC)
          CALL MTZSME(LMPOT,NATYP,NSPIN,VONSEQ,REQ,DRDIEQ,IMT,IRCUEQ,
     +                IPANEQ,NTCELL,LMSP,IRWSEQ,VSPSMO,LSMEAR,NSTART)

c     VONSEQ(i,1,*) = sqrt(4pi)*V_spherical
c
c
c     flushes the buffered output
          CALL FLUSH(6)
c

          CALL MIXSTR(RMSAVQ,RMSAVM,INS,LPOT,LMPOT,NATREF,NSHELL,NATPS,
     +                NEND,NSPIN,ITC,RFPI,FPI,IPF,MIXING,FCM,IRCEQ,
     +                IRMIEQ,REQ,DRDIEQ,VONSEQ,VISP,VINS,VSPSMO,VSPSME,
     +                LSMEAR)

c     Only in the first iteration can lsmear be equal to one, after
c     that the smeared spherical potential is already stored in the
c     potential file
          IF (LSMEAR.EQ.1) LSMEAR = 3
          IF (LSMEAR.EQ.2) LSMEAR = 4

          IF (MAX(RMSAVQ,RMSAVM).GE.QBOUND) THEN
c
c----> chebycheff mixing scheme
c
              IF (IMIX.EQ.2 .OR. (IMIX.EQ.1.AND.ITC.GE.19))
     +            CALL CHEACC(VISP,VINS,VONSEQ,VSPSME,VSPSMO,MIXING,
     +            NATPS,NATYP,NSPIN,INS,IRMIEQ,IRCEQ,IPF,LMPOT,LSMEAR)
c
c----> broyden updating schemes
c
              IF (IMIX.GE.3) THEN
                  DO 390 IA = NATPS,NATYP
                      IS = IA - NATREF
                      IRC1 = IRCEQ(IA)
                      ATWGHT(IA) = REAL(NSHELL(IS))/REAL(NATOM)
  390             CONTINUE
                  CALL BRYDBM(VISP,VONSEQ,VINS,VSPSME,VSPSMO,INS,LMPOT,
     +                        REQ,DRDIEQ,MIXING,ATWGHT,IRCEQ,IRMIEQ,
     +                        NSPIN,NATPS,NATYP,40,IMIX,20,IPF,LSMEAR)
              END IF
c     flushes the buffered output
              CALL FLUSH(6)
c
c
c----> reset to start new iteration
c
              DO 430 I = NPTPS,NSPPOT
                  IF (NSPIN.EQ.2) THEN
                      IT = (I+1)/2

                  ELSE


                      IT = I
                  END IF

                  IRC1 = IRCEQ(IT)
                  IF (LSMEAR.GT.0) THEN
                      DO 400 J = 1,IRC1
                          VISP(J,I) = VONSEQ(J,1,I)
                          VSPSME(J,I) = VSPSMO(J,I)
  400                 CONTINUE
                  ELSE
                      DO J = 1,IRC1
                          VISP(J,I) = VONSEQ(J,1,I)
                          VSPSME(J,I) = VISP(J,I)
                      END DO
                  END IF

                  IF (INS.NE.0 .AND. LPOT.GT.0) THEN

                      IRMIN1 = IRMIEQ(IT)
                      DO 420 LM = 2,LMPOT
                          DO 410 J = IRMIN1,IRC1
                              VINS(J,LM,I) = VONSEQ(J,LM,I)
  410                     CONTINUE
  420                 CONTINUE

                  END IF

  430         CONTINUE



          END IF

          REWIND 11

c
          CALL RITES(11,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,ITITLE,
     +               REQ,DRDIEQ,VISP,VSPSME,IRWSEQ,A,B,TXC,KXC,INS,IRNS,
     +               LPOT,VINS,QBOUND,IRCEQ,KSHAPE,EF,VBC,ECORE,LCORE,
     +               NCORE,IGGA,LSMEAR,ALFSME)
          CLOSE (11)

          OPEN (11)
c
          IF ((KLATR.EQ.2) .OR. (KLATR.EQ.3)) THEN
              DNEINT = DNEINT + DNEIHO
          END IF
c
c     Use correct int( n(r)V_eff(r) dr ):
c     This should be calculated using the same potential as in the
c     solution of radial equation (i.e., in this case the potential
c     on the equally spaced mesh) so replace epotin by epotin.
c
          IF (KTE.EQ.1) CALL ETOTAL(DNEINT,DSELOC,ECOU,EFERMI,EGR,
     +                              EPOTIN,ESPC,ESPV,EXC,IREF,KPRE,KFSR,
     +                              LMAX,LPOT,NATPER-ICUT,NATREF,NSPIN,
     +                              NSHELL)
c
c
          ETIME = DCLOCK()
          WRITE (6,FMT=*) 'time for mixing etc. ',ETIME - STIME
c
c     flushes the buffered output
          CALL FLUSH(6)
c
          IF (MAX(RMSAVQ,RMSAVM).LT.QBOUND) GO TO 450
  440 CONTINUE
c     loop over iterations ends


c
c              ------------------------------------------------
c             |  end of do loop over selfconsitency iterations  |
c              ________________________________________________
c
  450 CONTINUE
c
cholger status is a test condition

      IF (KMOLD.EQ.1 .AND. KF.EQ.1 .AND. KLATR.GT.0 .AND.
     +    MAX(RMSAVQ,RMSAVM).LT.QBOUND) THEN

          DO I = 1,NATOM
              DO J = 1,3
                  TAUXYZ(J,I) = RM(J,I)
                  TAU(J,I) = RM1(J,I)
                  TAUPRO(J,I) = RM(J,I)
              END DO
          END DO

          ICNT = 1
          DO I = 1,NATPER
              DO J = 1,NSHELL(I)
                  IFATOM(ICNT) = ICNT
                  IF (KATDYN(I).EQ.1) THEN
                      IATDYN(ICNT) = ICNT
                  ELSE IF (KATDYN(I).EQ.0) THEN
                      IATDYN(ICNT) = -1
                  END IF
                  ICNT = ICNT + 1
              END DO
          END DO

          OPEN (600,FILE='lastiter')
          READ (600,FMT=9120) ITER,MIT,STEPLN
          WRITE (6,FMT=*) ITER,MIT,STEPLN
          ITER = ITER + 1

          CALL MDYNCS(IOBROY,IOBROG,IOFILE,KMDYN,KPRI,KTEST,K_PRPFX,
     +                TAUPRO,IATDYN,NFATOM_DYN,ITER,IFATOM,ITDMD,STEPL,
     +                DELTAMAX,T_DEBYE,NFATOM,DIMASS,EPS_MD,K_EPS_MD,
     +                K_ADPT_DT,K_SET_F_M,TAU,TAUXYZ,F1XYZ,MIT,F_OLD,
     +                STEPLN)

          CALL INPOUT(FILNAM,IRESIST,IFILEH,IPE,IWTMAT,NSPIN,IRM,INS,
     +                ICST,KCOR,KVREL,KWS,KHYP,KHFELD,KFSR,KXC,KTE,KPRE,
     +                KSPH,KEFG,KVMAD,KF,IGGA,NATREF,NATPER,ICUT,
     +                KSHAPEH,KMOLD,IREF,NTCELL,IRNS,ITCLST,IMIX,IEF,
     +                STRMIX,FCM,QBOUND,BRYMIX,EF,HFIELD,VBC,VCONSTH,
     +                KLATR,KSYMM,LKONV,KSYMMAT,KESYM,OCCUP,ASARY,
     +                TAUXYZ,RM1,KMDYN,K_PRPFX,NFATOM_DYN,NFATOM,ITDMD,
     +                K_EPS_MD,K_ADPT_DT,K_SET_F_M,KPRI,KTEST,IOBROY,
     +                IOBROG,IOFILE,STEPL,DELTAMAX,T_DEBYE,DIMASS,
     +                EPS_MD,KATDYN,SOCCUP,KOCCUP,QOCCUPH)

          REWIND (600)
          WRITE (600,FMT=9120) ITER,MIT,STEPLN

      ELSE IF (KMOLD.EQ.1 .AND. (KF.NE.1.OR.KLATR.EQ.0) .AND.
     +         RMSAVM.GE.QBOUND) THEN

          WRITE (6,FMT=*)
     +  'mdyncs: electronic structure is not converged or you have
     +   forgotten to set kf eq 1 and klatr gt 0 !'

      END IF

      IF (KSPH.EQ.1) THEN
c
c---> calculate specific heat and change induced by defects
c
          CALL WKSPH(DEN,IREF,NLST,NSHELL,TEXTS,PI,LMAX,NATPS,NATREF,
     +               NEND,NSPIN)
      END IF
c
c
c---> store densities of state
c
      IF (IEF.GT.0) THEN

c     Write unconvoluted charge density (=rho2ns/r^2)
c        Rewind 49
          OPEN (49)
          CALL WLRHO(49,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,RWS,ITITLE,REQ,
     +               DRDIEQ,R2NSEQ,RHOC,IRWSEQ,A,B,TXC,KXC,IRNS,LPOT,
     +               IRCEQ,KSHAPE,EF,VBC,ECORE,LCORE,NCORE,IGGA)
          CLOSE (49)
c
c     Write LDOS
          CALL WLDOS(DEN,EF,IEN,ITITLE,PI,IA,IELAST,LMAX,NATYP,NPTPS,
     +               NSPIN)
      END IF
c
      ETIMETOT = DCLOCK()
      WRITE (6,FMT=*) 'total time used ',ETIMETOT - STIMETOT
c
      WRITE (6,FMT=*) 'Times in different parts:'
      WRITE (6,FMT=*) ' REDOCH ',TTIMRE
      WRITE (6,FMT=*) ' CPLXWF ',TTIMCP
      WRITE (6,FMT=*) ' TMATNS ',TTIMTM
      WRITE (6,FMT=*) ' SYMMET ',TTIMSY
      WRITE (6,FMT=*) ' DYSON  ',TTIMDY
      WRITE (6,FMT=*) ' RHO    ',TTIMRH
      WRITE (6,FMT=*) ' ENERGY ',TTIMES
c
      RETURN
  460 WRITE (6,FMT=*) 'rotation matrix',IDND(II,N),'is not o.k.'
      STOP 'rotation'


 9000 FORMAT (48X,I3)
 9120 FORMAT (2i5,1f10.4)
 9320 FORMAT (A50)
      END
