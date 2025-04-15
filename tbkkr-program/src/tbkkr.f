c main.f (for self consistent iterations)
c 23.04.96 ***************************************************************
C                            MAIN PROGRAM
c 1.6.99  Reconstructed...
c ************************************************************************
C     
C 
      PROGRAM TBKKR
c	include 'cxml_include.f90'
      IMPLICIT NONE
      include 'inc.fi' 
      include 'inc.cls'
C 
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMMAXD= (LMAXD+1)**2,LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER NREFPOTD
      PARAMETER (NREFPOTD=NSPIND*NREFD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LASSLD
      PARAMETER (LASSLD=4*LMAXD)
      INTEGER LM2D
      PARAMETER (LM2D = (2*LMAXD+1)**2)
      INTEGER NCPAIRD
      PARAMETER(NCPAIRD=10)
c    parameter QDOSKDIM is used for q-dos plots (set to 1 to save space)
      INTEGER QDOSKDIM
      PARAMETER (QDOSKDIM=2000)
      DOUBLE PRECISION ONE,ONEM,TWO
      PARAMETER (ONE = 1.D0,ONEM = -1.D0,TWO=2.D0)
      DOUBLE COMPLEX CONE,CONEM,CZERO,CI
      PARAMETER (CONE  = ( 1.0D0,0.0D0),
     +           CONEM = (-1.0D0,0.0D0),
     +           CZERO = ( 0.0D0,0.0D0),
     +           CI    = ( 0.0D0,1.0D0))
C     ..
C     .. Local Scalars ..
c
      LOGICAL       
     +     COMPLX,
     +     LINIPOL,LSTART,LRHOSYM,
     +     OPT,
     +     TEST 
C
      CHARACTER*1          HEAD(80)
      CHARACTER*8          DAT,VERS
      CHARACTER*26         STARTDATE
      CHARACTER*40         I12,I13,I19,I25,I40
      CHARACTER*9          I16,I17,I18,I21,I22,I23
C
C     .. DOUBLE PRECISION ....
C
      DOUBLE PRECISION 
     +     ALATC,BLATC,CLATC,       ! lattice constants (in a.u.)
     +     C,                       ! speed of light
     +     ABASIS,BBASIS,CBASIS,    ! scaling factors for rbasis
     +     DOSTOT(0:LMAXD),         ! total DOS for output
     +     DOS,                     ! integrated DOS for output
     +     FAC,FACL,                ! factor (+/-1.0d0)
     +     GAMMA,                   ! lattice distortion
     +     LATSP,                   ! lattice distortion
     +     E3,
     +     E1,E2IN,E2,              ! energies needed in EMESHT
     +     EFERMI,                  ! Fermi energy (for scaling DOS)
     +     HFIELD,                  ! external magnetic field, for 
                                    ! initial potential shift in 
                                    ! spin polarised case
     +     STIME,STIME0,STIME1,     ! real time
     +     STIMEK,STIMEREF,sumdf(NATYPD),espcor,factor,espco2
      DOUBLE PRECISION 
     +     ETIME,ETIME0,ETIME1,
     +     VOLUC,                   ! volume of unit cell in 
                                    ! units of alat**3
     +     TK,                      ! temperature
     +     KB,                      ! Boltzmann constant
     +     VCONST,                  ! potential shift
c
     +     ALAT,ALATNEW,EIN2,ESHIFT,FCM,MIX,MIXING,MIXING0,QBOUND,
     +     CHRGNT,DENEF,E2SHIFT,EFCTOR,
     +     CHRGOLD,DENOLD,DENI,
     +     F,FF,FPI,PI,RFPI,RTEMP,RMSAVM,RMSAVQ,RMSAV0,RMSAV1,SIGN,X,
     +     RFOURIER,                ! cutting radius for Fourier-transf.
     +     ebotsc,eupsc,            ! energy in Ryd for begining\end of semic. contur
     +     eimsc,                   ! goes to this energy in the imag part of semic. contur
     +     ebotgp,eupgp,            ! energy in Ryd for begining\end of gap contur
     +     eimgp                    ! goes to this energy in the imag part of gap contur
C
C     .. DOUBLE PRECISION ARRAYS ....
C
      DOUBLE PRECISION
     +     ATWGHT(NATYPD),
     +     BRAVAIS(3,3),            ! bravais lattice vectors
     +     RECBV(3,3),              ! reciprocal basis vectors 
     +     MTFAC(NATYPD),           ! scaling factor for radius MT
c
     +     A(NATYPD),B(NATYPD),     ! contants for exponential r mesh
     +     CLEB(NCLEB,2),           ! GAUNT coefficients (GAUNT)
     +     DRDI(IRMD,NATYPD),       ! derivative dr/di
     +     DROR(IRMD,NATYPD),       ! logarithmic derivative (1/r)*(dr/di)
     +     ECORE(20,NPOTD),         ! core states
     +     R(IRMD,NATYPD),          ! radial r mesh (in units a Bohr)
     +     RBASIS(3,NAEZD+NEMBD),   ! position of atoms in the unit cell
                                    ! in units of bravais vectors
     +     RBASIS1(3,NAEZD+NEMBD),  ! pos. of atoms in atomic units
     +     RCLS(3,NACLSD,NCLSD),    ! real space position of atom in cluster
     +     RCLS1(3,NACLSD),         ! only for test use
     +     RMT(NATYPD),             ! Muffin-Tin-radius
     +     RMTNEW(NATYPD),          ! adapted MT radius
     +     RR(3,0:NRD),             ! set of real space vectors (in a.u.)
     +     RS(IRMD,0:LMAXD,NATYPD), ! r mesh for relat. calc.
     +     RWS(NATYPD),             ! Wigner Seitz radius
     +     S(0:LMAXD,NATYPD)
      DOUBLE PRECISION
     +     THETAS(IRID,NFUND,NCELLD), ! shape function
                                      !         ( 0 outer space
                                      ! THETA = (
                                      !         ( 1 inside WS cell
                                      ! in spherical harmonics expansion
     +     VBC(2),                  ! potential constants
     +     VINS(IRMIND:IRMD,        ! nonsperical potential
     +                 LMPOTD,NSPOTD),       
     +     VISP(IRMD,NPOTD),        ! spherical input potential
     +     VONS(IRMD,LMPOTD,NPOTD), ! output potential
     +     WG(LASSLD),              ! integr. weights for Legendre pol.
                                    ! (GAUNT2)
     +     YRG(LASSLD,0:LASSLD,
     +                0:LASSLD),    ! spherical harmonics (GAUNT2)
     +     Z(NATYPD)               ! nuclues charge
c
      DOUBLE PRECISION
     +     CMINST(LMPOTD,NATYPD),   ! charge moment of interstitial
     +     CMOM(LMPOTD,NATYPD),     ! LM moment of total charge
     +     ECOU(0:LPOTD,NATYPD),    ! Coulomb energy
     +     EOLD(2),ENEW(2),         ! old, new FERMI energy
     +     EPOTIN(NATYPD),          ! energy of input potential (EPOTINB)
     +     ESPC(0:LMAXD,NATYPD,NSPIND), ! energy single particle core
     +     ESPV(0:LMAXD,NATYPD,NSPIND), ! energy single particle valence
     $                 espsc1(0:lmaxd,natypd,nspind),
     $                 espsc2(0:lmaxd,natypd,nspind), 
     +     EXC(0:LPOTD,NATYPD),     ! E(XC)
     +     GSH(NGSHD),
     +     R2NEF(IRMD,LMPOTD,NATYPD,NSPIND), ! rho at FERMI energy
     +     RHO2NS(IRMD,LMPOTD,NATYPD,NSPIND), ! radial density
c                     nspin=1 : (*,*,*,1)  radial charge density
c                     nspin=2 : (*,*,*,1)  rho(2) + rho(1)  -> charge
c                               (*,*,*,2)  rho(2) - rho(1)  -> mag. moment
     +     RHOC(IRMD,NPOTD),        ! core charge density
     +     RIJ(IJD,3),
     +     SUMNS(NATYPD),           ! sum of valence charge
     +     VASUM(0:LMAXD,NATYPD),   ! sum of l-dep. valence charge
c    p. zahn, 10.3.99,
c     +     WORK(IRMD,LMPOTD,NATYPD),
     +     WORK,
     +     WTYR(IJD,LMPOTD),
     +     YR(IJD,LMPOTD)
      DOUBLE PRECISION C00(LMPOTD)   ! added 29.10.99
      double precision FLM(-1:1,NATYPD),         ! FORCES    added 20.04.00
     +                 FLMC(-1:1,NATYPD)
c
       DOUBLE PRECISION
     +     MIXSAV(200),MIXW(200),MIXV(200),MIXC(3),XX,FM ! mixing optimization 
c                                                          by p.z. ??
c
c   new madelung variables
       INTEGER NGMAX,NRMAX,NSHLG,NSHLR
       INTEGER NSG(ISHLD),NSR(ISHLD)
       DOUBLE PRECISION GN(4,NMAXD),RM(4,NMAXD),
     &                  MADELSMAT(NAEZD,NAEZD,LMXSPD)
       DOUBLE PRECISION VOLUME0,RCUTZ,RCUTXY,TEMP(3)

 
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DOUBLE PRECISION SMALL
      DOUBLE COMPLEX DELTA,T
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     
C     .. INTEGER ....
C
      integer
     +     INTERVX,INTERVY,INTERVZ, ! number of intervals in x-,y-,z-direction
                                    ! for k-net in IB of the BZ
     +     ICC,                     ! center of cluster for output of GF
c                                   ! for DOS calc. icc = 0 
     +     ICLS,                    ! counter for cluster
     +     ICST,                    ! number of Born approximation
     +     IELAST,                  ! number of complex energy values
     +     IEND,                    ! number of nonzero gaunt coeffizients
     +     IFILE,                   ! unit specifier for potential card
     +     IINFO,IINFO1,            ! info's about r-mesh in CALRMT
     +     IN,                      ! I/O channel for system potential
     +     SINN,SOUT,RIN,ROUT,      ! in/out of GF, WF, T-matrices
     +     INS,                     ! 0 (MT), 1(ASA), 2(Full Potential)
     +     INSREF,                  ! INS for reference pot. (usual 0)
     +     IOWFCT,IOWREF,           ! unit specifier for wave functions
     +     IOPREF,IOPSYS,           !                    potentials
     +     IOTREF,IOTSYS,           !                    T matrices
     +     IOGREF,                  !                    GF(ref sys.)
     +     IPE,IPF,IPFE,            ! not real used, IPFE should be 0
     +     IGF,IMIX,IPOTOU,IPRCOR,
     +     IRM,IRNUMX,
     +     KSCOEF,                  ! 0,1: read shell structure from
c                                   !      file 25
     +     ISHIFT,ITCCOR,ITCLST,ITDBRY,
     +     KHFELD,                  ! 0,1: no / yes external magnetic field
     +     KSHAPE,                  ! exact treatment of WS cell
     +     KVREL,                   ! 0,1 : non / scalar relat. calculation
     +     KWS,                     ! 0 (MT), 1(ASA)
     +     KMT,                     ! scaling of RMT with MTFAC
c                                   ! 0: RMT from crystal structure 
c                                   ! 1: RMT from crystal structure
c                                   !    scaled by RMTFAC
c                                   ! 2: RMT = RMTFAC
c                                   ! 3: RMT from ref. pot. card
     +     KCOR,KEFG,KFROZN,
     +     KHYP,KPRE,KTE,KVMAD,KXC,KFORCE
      INTEGER
     +     LATT,                    ! 0: scu  1: fcc  2: bcc
c                                   ! 3: tet  4: bct  5: hex
c                                   !                 8: ort
     +     LMAX,                    ! maximum l component in 
c                                   ! wave function expansion
     +     LPOT,                    ! maximum l component in 
c                                   ! potential expansion
     +     LMMAX,LMPOT,MD,LTSP,
     +     MAXMESH,SAVEMESH,        ! number of k-mesh sets
     +     MMIN,MMAX,               ! min. max. eigenvalue
     +     M2,                      ! maximum number of bands
c
     +     NMESH,                   ! number of k set 
     +     NAEZ,                    ! number of atoms in unit cell
     +     NATYP,                   ! number of kinds of atoms in unit cell
     +     NCLS,                    ! number of reference clusters
     +     NEMB,                    ! number of 'embedding' positions
     +     NEMBZ,                   ! inequiv. 'embedding' positions 
     +     NINEQ,                   ! number of ineq. positions in  unit cell
     +     NLAYER,                  ! number of principal layer
     +     NPNT1,NPNT2,NPNT3,       ! number of E points (EMESHT)
     +     NPOL,                    ! number of Matsubara Pols (EMESHT)
     +     NR,                      ! number of real space vectors rr
     +     NREF,                    ! number of diff. ref. potentials
     +     NSPIN,                   ! counter for spin directions
     +     NSTEPS,                  ! number of iterations
c     +     NVMAD,                   ! number of madelung shells
     +     NZ,                      ! number of atoms at centers of inversion
     +     STOP_MARK,               ! dimension tests
c 
     +     I,J,IK,                  ! counters
     +     IPOT,ISPIN,
     +     I1,IC,IE,II,INV,             
     +     IA,IATYP,IDIM,IH,IJEND,IP,IPARA,
     +     IR,IRC1,IRMIN1,IT,ITC,IWRIT,
     +     L,LM,LM1,LM2,
     +     N,N1,N2,NITMIX,NMIX,
     +     RF,IHANDLE,IHANDLE1,NKDOS,
     +     npntsc,                    ! number of energy points in real part  (sc)
     +     npntim,                    ! number of energy points in imag. part (sc)
     +     npntgp,                    ! number of energy points in real part  (gap)
     +     npntimgp,                   ! number of energy points in imag. part (gap)
     +     igga                       ! number of pw91-gga (=1) 
C     
C     .. INTEGER ARRAYS ....
C
      INTEGER
     +     EQINV(NAEZD),            ! site equiv. by invers. symmetry
     +     INIPOL(NATYPD),          ! initial spin polarisation
     +     IXIPOL(NATYPD),          ! constraint of spin pol.
     +     KAOEZ(NAEZD+NEMBD),      ! kind of atom at site in elem. cell
c
c     ..  cluster arrays
c
     +     ATOM(NACLSD,NAEZD),     ! atom at site in cluster
     +     CLS(NATYPD),             ! cluster around atom
     +     COCLS(NCLSD),            ! center of cluster (kaez )
     +     NACLS(NCLSD),            ! number of atoms in cluster
     +     EZOA(NACLSD,NAEZD),      ! EZ of atom at site in cluster
c
     +     ICLEB(NCLEB,4),          ! pointer array
     +     IFUNM(NATYPD,LMXSPD),
     +     IMT(NATYPD),             ! r point at MT radius
     +     IPAN(NATYPD),            ! number of panels in non-MT-region
     +     IRC(NATYPD),             ! r point for potential cutting
     +     IRCUT(0:IPAND,NATYPD),   ! r points of panel borders
     +     IRMIN(NATYPD),           ! max r for spherical treatment
     +     IRNS(NATYPD),            ! number r points for non spher. treatm.
     +     IRWS(NATYPD),            ! r point at WS radius
c     +     IMAD(NAEZD,NAEZD),       ! madelung constant for 
                                    ! position pair(i,j)
     +     ITITLE(20,NPOTD)         ! title line in potential card
      INTEGER
c     +     JMAD(NVMADD,2*NVMADD+1), ! position indices for madelung constant 
     +     JEND(LMPOTD,             ! pointer array for icleb()
     +             0:LMAXD,0:LMAXD),
     +     KFG(4,NATYPD),
     +     LCORE(20,NPOTD),         ! angular momentum of core states
     +     LLMSP(NATYPD,NFUND),     ! lm=(l,m) of 'nfund'th nonvanishing
                                    ! component of non-spherical pot.
     +     LMSP(NATYPD,LMXSPD),     ! 0,1 : non/-vanishing lm=(l,m) component
                                    ! of non-spherical potential
     +     LOFLM(LM2D),             ! l of lm=(l,m) (GAUNT)
     +     LMXC(NATYPD),
     +     NCORE(NPOTD),            ! number of core states
     +     NFU(NATYPD),
     +     NSHELL(0:NSHELD),        ! number of atoms in shell
c                                   ! nshell(0) = number of shells
     +     NTCELL(NATYPD),          ! index for WS cell 
     +     NTCELLR(NREFD),          ! index for WS cell of ref. pot.
     +     REFPOT(NATYPD+NEMBD),    ! ref. pot. card  at position
c     
     +     ILM(NGSHD,3),
     +     IMAXSH(0:LMPOTD),ATOMIMP(NATOMIMPD)
      INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR,IEGFOUT ! for conductivity
C                                                                   ! calculation
C
C     .. DOUBLE COMPLEX ....
c     
      DOUBLE COMPLEX DF,E,EK,CFAC,CFACL
c
C
C     .. DOUBLE COMPLEX ARRAYS ....
c     
      DOUBLE COMPLEX 
c
     +     ALPHA(0:LMAXD,NATYPD),
     +     AR(LMMAXD,LMMAXD,NATYPD),
     +     CR(LMMAXD,LMMAXD,NATYPD),
     +     DEN(IEMXD,0:LMAXD,NATYPD,NSPIND), ! density of states
     +     DETALF(NATYPD,IEMXD,NSPIND),
     +     DEZ(IEMXD),              ! weights for complex energy integration
     +     EZ(IEMXD),               ! energy points for compl.integration
     +     FZ(IRMD,0:LMAXD,NATYPD),
     +     GMATLL(LMMAXD,LMMAXD,NSHELD), ! diag elements of G matrix(system)
     +     PNS(LMMAXD,LMMAXD,       ! non-sph. eigen states of single pot
     +                 IRMIND:IRMD,2,NATYPD),
     +     PZ(IRMD,0:LMAXD,NATYPD),
     +     QZ(IRMD,0:LMAXD,NATYPD),
     +     SZ(IRMD,0:LMAXD,NATYPD),
     +     QNS(LMMAXD,LMMAXD,       ! non-sph. eigen states of single pot
     +                 IRMIND:IRMD,2,NATYPD), 
     +     GINP(NACLSD*LMMAXD,LMMAXD,NCLSD), ! cluster GF(ref syst.)
     +     WEZ(IEMXD)               ! modified weights for energy integration
                                    ! = -2.0/pi*DEZ
c added 17.07.2000 to simulate finite basis set
      double complex ARREF(LMMAXD,LMMAXD),CRREF(LMMAXD,LMMAXD),
     &             PZREF(IRMD,0:LMAXD),QZREF(IRMD,0:LMAXD),
     &             SZREF(IRMD,0:LMAXD),FZREF(IRMD,0:LMAXD),
     &             PNSREF(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     &             QNSREF(LMMAXD,LMMAXD,IRMIND:IRMD,2)


c
      DOUBLE COMPLEX 
     +     TLL(LMMAXD,LMMAXD),
     +     TMATLL(LMMAXD,LMMAXD,NAEZD),
     +     TMATLL1(LMMAXD,LMMAXD,NATYPD),
     +     TREFLL(LMMAXD,LMMAXD,NREFD),
     +     DELTALL(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX LEFTTINVLL(LMMAXD,LMMAXD,NEMBD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD)
      DOUBLE COMPLEX GSQDOS(LMMAXD,LMMAXD,NAEZD,QDOSKDIM)  ! 19.4.2000 
      DOUBLE PRECISION QDOSKP(3,QDOSKDIM)     ! 19.4.2000
      INTEGER NQDOSKP           ! 19.4.2000
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Lines with a !new1 in the end are added after 
c    20.9.99
      INTEGER NLAY,                               !new1 
     +        NLBASIS,NRBASIS,                    !new1 
     +        NLEFT,NRIGHT,IER                    !new1
      DOUBLE PRECISION                            !new1 
     +     ZPERLEFT(3),ZPERIGHT(3),               !new1
     +     TLEFT(3,NEMBD),TRIGHT(3,NEMBD)         !new1
      DOUBLE PRECISION CMOMHOST(LMPOTD,NEMBD)     !new1
      double precision tx,ty,tz
c     zvec(3,nlay)       : coordinates of the layers
c     nlay               : number of layers
c     nlbasis            : number of basis layers of left 
c                          host (repeated units)
c     nleft              : number of repeated basis 
c                          for left host to get converged  
c                          electrostatic potentials
c
c     zperleft(3)        : vector to define how to repeat
c                          the basis of the left host 
c     tleft(3,nlbasis)   : vectors of the basis for the 
c                          left host.
c     nrbasis            : number of basis layers of right 
c                          host (repeated units)
c     nright             : number of repeated basis 
c                          for right host to get converged  
c                          electrostatic potentials
c
c     zperight(3)        : vector to define how to repeat
c                          the basis of the right host 
c     tright(3,nlbasis)  : vectors of the basis for the 
c                          right host.
c
c     cmomhol(lmpot,nlbasis) : charge moments of each atom
c                              of the left host
c     
c     cmomhor(lmpot,nrbasis) : charge moments of each atom
c                              of the right host
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      DOUBLE COMPLEX BESSJW1(0:LMAXD+1),BESSYW1(0:LMAXD+1),
     &               HANKWS1(0:LMAXD+1),BESSJW2(0:LMAXD+1),
     &               BESSYW2(0:LMAXD+1),HANKWS2(0:LMAXD+1)
      double precision vrep,rmt1
      double complex a1,b1,tmatanal(0:LMAXD)
      logical lcall

ccccccccccccccccccc
      LOGICAL STRCO1,LINTERFACE,LCARTESIAN,VACFLAG(2)        !new1
      CHARACTER*5 LMACH
      CHARACTER*50 HOME
      CHARACTER*80 UIO
c
c
C     .. Common blocks ..
c
c      COMMON /CLL/CLEB,ICLEB,LOFLM,IEND
      COMMON /CMACH/STRCO1,LMACH
c      COMMON /PARA/AVMAD,BVMAD,ALAT,E1,EIN2,TK,ESHIFT,FCM,HFIELD,MIXING,
c     +       QBOUND,VBC,TXC,KFG,LMXC,IRNS,NTCELL,ICST,IFILE,IGF,IMIX,
c     +       INS,IPE,IPF,IPFE,IPOTOU,IPRCOR,IRM,IRNUMX,ISHIFT,ITCCOR,
c     +       ITCLST,ITDBRY,KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
c     +       KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,MD,NATYP,NSPIN,
c     +       NPOL,NPNT1,NPNT2,NPNT3
      COMMON /RFOUR/ RFOURIER
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
     +         NREF,NATYP,NAEZ,NINEQ,NSHELL,NZ,NLAYER,
     +         IMT,IRWS,IRNS,NTCELL,NTCELLR,
     +         REFPOT,NCLS,NACLS,
     +         CLS,EZOA,ATOM,EQINV,KAOEZ,INIPOL,LATT,ICC,ICST
c
c     .. EXTERNAL STATEMENT
c
      EXTERNAL   
     +     CLSGEN99,DELTAMAT,
     +     RTEST,ROPT,TEST,OPT,
     +     RINPUT1,TMATRX,
     +     TMWRIT,WFWRIT,
c
     +     BRYDBM,CHEACC,CINIT,CONVOL,CORELB,CRECV,CSEND,
c    +     DAXPY,
c    +     DCOPY,
     +     ECOUB,EMESHT,EPOTINB,ESPCB,ESPVB,ETOTB1,FORFLUSH,
     +     GAUNT,GAUNT2,GDSUM,KLOOPZ1,MIXSTR,MTZERO,RCSTOP,RHOLMB,
     +     RHONSB,RHOTOTB,RINIT,RITES,RMADEL,SHAPE,SPHERE,
     +     STARTB1,TMREAD,VINTRAS,VMADEL,VXCLM,WFREAD,
c
     +     cylm02,trarea,vxcgga,gradrl,gradr,mkxcpe,gxcpt,cpw91,
     +     corlsd,exch91,gcor91
c
C     .. Intrinsic Functions ..
      INTRINSIC ACOS,DABS,DATAN,DIMAG,DMAX1,
     +          DMIN1,DSIGN,DSQRT,LOG,MAX,REAL,SIN,ZABS,ZEXP
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DCLOCK
      DOUBLE COMPLEX EXPIDL
      INTEGER MAPBLOCK,MYNODE,NUMNODES
      EXTERNAL DCLOCK,EXPIDL,MAPBLOCK,MYNODE,NUMNODES
	real*8 rx(1000),ry(1000),rz(1000)
	integer numdiff,typ1(1000),typ2(1000),ix,iy,iz,
     +numat1(1000),numat2(1000)
	common/sym_file/rx,ry,rz,numdiff,typ1,typ2
     +,NUMAT1,NUMAT2
C     ..
      CHARACTER*24 TXC(3)
C     ..
C     .. Data statements ..
c
      DATA ATWGHT   / NATYPD*1.D0  /,
     +     DAT      / '        '   /,
     +     VERS     / '01.06.99'   /,
     +     SMALL    / 1.0D-8       /,
     +     KB       / 0.6333659D-5 /
c
      DATA I16 /'tmat.ref '/,
     +     I17 /'gll.ref  '/,
     +     I18 /'wfct.ref '/,
     +     I21 /'tmat.sys '/,
     +     I22 /'gll.sys  '/
c
C     .. Save statement ..
c
      SAVE
c
c ---------------------------------------------------CONSTANTS -----------
      PI = 4.0D0*DATAN(1.0D0)
      FPI = 4.0D0*PI
      RFPI = DSQRT(FPI)
      EFERMI = 0.0d0
      RMSAV0 = 1.0d10
      NITMIX = 0
c ------------------------------------------------------------------------
      STIME0 = DCLOCK()
c ------------------------------------------------------------------------
c   sym_file   file29
      WRITE(6,*) ' write sym_file'
      OPEN(29,FILE='sym_file',form='formatted')
      rewind 29
	read(29,*) numdiff
	WRITE(6,*) ' NUMDIFF=',NUMDIFF
	if(numdiff.gt.1000) then
	write(6,*) ' numdiff is too large',numdiff
	endif
	do 2111 ix=1,numdiff
	read(29,*) IY,NUMAT1(IX),NUMAT2(IX),RX(IX),RY(IX),
     +RZ(IX),TYP1(IX),TYP2(IX)
	IF(IY.GT.NUMDIFF) STOP
	WRITE(6,2112) NUMAT1(IX),NUMAT2(IX),RX(IX),
     +RY(IX),RZ(IX),TYP1(IX),TYP2(IX)
 2111 CONTINUE
 2112 FORMAT(1X,' NUMAT1-2,rx,ry,rz,typ1-2',2I3,3F10.5,2I5)

c
choshino  file 91 _____  green
      open(91,FILE=' green',FORM='UNFORMATTED')
	rewind 91
c     >>>>  TEST <<<<<
c      write(6,*) 'lmaxd,m2d,irmd ',
c     +     lmaxd,m2d,irmd
c      write(6,*) 'ispard,ifulld,islabd,natypd,naezd,nzd',
c     +     ispard,ifulld,islabd,natypd,naezd,nzd
c      stop

c ----------------------------------------------------- INPUT ------------
       CALL RINPUT99(ALAT,RBASIS,ABASIS,BBASIS,CBASIS,
     &           ALATC,BLATC,CLATC,LATT,CLS,NCLS,
     &           E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,ESHIFT,
     &           ITCLST,NSTEPS,IMIX,MIXING,QBOUND,FCM,ITDBRY,
     &           IRNS,NTCELL,NAEZ,NEMB,KAOEZ,EQINV,IRM,Z,
     &           NINEQ,NREF,NTCELLR,
     &           ICST,IFILE,IGF,INS,INSREF,IPE,IPF,IPFE,
     &           KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     &           KFG,KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT, 
     &           NATYP,NSPIN,
     &           LMXC,TXC,KSCOEF,ICC,REFPOT,
     &           IPOTOU,IPRCOR,IRNUMX,ISHIFT,ITCCOR,
     &           MD,INTERVX,INTERVY,INTERVZ,
     &           HFIELD,COMPLX,
     &           KMT,MTFAC,VBC,VCONST,LINIPOL,INIPOL,IXIPOL,LRHOSYM,
     &           MMIN,MMAX,SINN,SOUT,RIN,ROUT,M2,I12,I13,I19,I25,I40,
     &           NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,    
     &           TLEFT,TRIGHT,LINTERFACE,RCUTZ,RCUTXY,EBOTSC,EUPSC,     ! new holger
     &           EIMSC,NPNTSC,NPNTIM,EBOTGP,EUPGP,EIMGP,NPNTGP,NPNTIMGP)! new holger     
c
c     write(6,*) 'cls',(cls(i1),i1=1,ncls)
c
       if(kxc.lt.0) then
       write(6,*) ' pw91-gga calculations (kexcor=kxc < 0) '
       else
       write(6,*) ' lsda calculations (kxc >= 0)'
       endif
c
c
c
       MIXING0 = MIXING
c
       CALL TESTDIM(nspin,naez,nemb,natyp,lmax,irm,ins,insref,
     &                    nref,M2,IRNS,ncls,nlayer,LATT)     
c
c      OPEN (13,FILE=I13,status='old',FORM='formatted')
      IF (INS.GT.0) OPEN (19,FILE=I19,status='old',FORM='formatted')
      OPEN (25,FILE=I25,status='unknown',FORM='formatted')
      OPEN (20,FORM='unformatted',STATUS='unknown')
      IF (.NOT.TEST('NoMadel ')) 
     +     OPEN (40,FILE=I40,STATUS='old',FORM='formatted')
      OPEN (61,FORM='unformatted',STATUS='unknown')    ! corrected 17.10.1999
c      OPEN (61,FORM='unformatted',STATUS='scratch')

C *********************************************Input-End ********

C ------------------------------------------------------------------------

c ------------------------------------------------------------------------
c
c ---> dependencies of input parameters
c
      IF (TEST('c/a     ')) CLATC = CLATC*ALATC
      ALAT = ALATC
      MD = intervx

      IF (ICC.NE.0 .AND. ITCLST.GT.1) THEN
        ITCLST =1
        write(6,*) 'CLUSTER CALCULATION : ITCLST is set to 1.'
      END IF

      IF (ICC.NE.0 .AND. IGF.EQ.0) THEN
        IGF = 1
        write(6,*) 'CLUSTER CALCULATION : ',
     +             'Output of Greens Function. IGF is set to 1.'
      END IF

      IF (ICC.EQ.0 .AND. IGF.NE.0) THEN
        ICC = -1
        write(6,*) 
     +       'Output of Cluster Greens Function. ICC is set to -1.'
      END IF

      IF (ICC.EQ.0 .AND. KSCOEF.NE.0) THEN
        KSCOEF = 0
        write(6,*) 'ELECTRON DENSITY calculation : '
        write(6,*) 'No input of Shell Structure. KSCOEF is set to 0.'
      END IF

      IF(ICC.GT.NINEQ) THEN
        write(6,*) 'The choice of ICC (=',ICC,') is not possible.'
        write(6,*) 'ICC is changed to the equivalent atom ',
     +       EQINV(ICC),'.'
        ICC = EQINV(ICC)
      ENDIF

c
c ---> switch on external splitting field if INIPOL is given
c
      IF (LINIPOL) KHFELD = 1
      IF (KHFELD.NE.0 .AND. IFILE.NE.13 ) THEN
        KHFELD = 0
        write(6,*) 'Input of a previous output potential : '
        write(6,*) 'KHFELD is changed to ',KHFELD,'.'
      END IF
c
c ---> cut Fourier Transformation at a given radius ( for test of the method)
c
      RFOURIER = 0.0D0
      IF (OPT('cut Four')) THEN
        RFOURIER = ESHIFT
        write(6,9300) RFOURIER
        write(6,2100) 
      END IF

 9300 format(' CUT cluster for FOURIER transformation at radius',
     +     f12.5,' .')

      IF(IPE.GT.0)THEN
        IPFE = 9
        open(IPFE,file='ipfe',status='unknown')
      ENDIF

      IF(OPT('rigid-ef').OR.OPT('DECIMATE'))THEN
        ISHIFT = 2
        write(6,*) 'Rigid Fermi Energy, ISHIFT is set to ',
     +       ISHIFT,'.'
      END IF

c
c ---> determination of properties at Fermi level
c
      IF (OPT('GF-EF   ')) THEN
        ITCLST =1
        write(6,*) 'CLUSTER CALCULATION at EF : ITCLST is set to 1.'
        IGF = 1
        write(6,*) 'CLUSTER CALCULATION : ',
     +             'Output of Greens Function. IGF is set to 1.'
        ICC = -1
        write(6,*) 
     +       'Output of Cluster Greens Function. ICC is set to -1.'
        IF (NPOL.GT.0) NPOL = 0
        IF (NPOL.LT.0) NPNT1 = 0
        IF (NPOL.LT.0) NPNT3 = 0
        NPNT2 = 1
      END IF
      IF (OPT('DOS-EF  ')) THEN
        ITCLST =1
        NPOL = 0
        NPNT2 = 1
      END IF
      IF (OPT('deci-out')) THEN
        ITCLST =1
        WRITE(6,*) ' *************************************'
        WRITE(6,*) ' *************************************'
        WRITE(6,*) ' ***    WRITING DECIMATION FILE  *****'
        WRITE(6,*) ' *************************************'
        WRITE(6,*) ' *************************************'
      END IF
c
c ---> for optimization of MIXING parameter simple mixing whitout
c      chebychev-mixing is used
c
      IF (TEST('opt mix ')) THEN
        IMIX = 0
        write(6,*) 'IMIX is set to ',IMIX
      END IF
c
c ---> end of dimension and input data check
c
      write(6,2100) 
c
c ------------------------------------------------------------------------
      LTSP=0   ! not used !!
      GAMMA=0.0D0 ! not used !!
c
c--->   determines the direct and reciprocal space basis vectors
c       and a number of direct space vectors RR(3,0:NRD)
c
c Changes in LATTIX start here. A new subroutine, LATTIX99, is introduced,
c that reads the bravais vectors directly from the input card and returns
c the same things as LATTIX: BRAVAIS, RECBV, RR, NR.
c BRAVAIS in units of au/alat
c
c
         CALL LATTIX99(ALATC,BRAVAIS,RECBV,RR,NR)
c
c---->  normalization of basis vectors
c
         CALL SCALEVEC(RBASIS,ABASIS,BBASIS,CBASIS,NLBASIS,
     &     NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT,             ! new1
     &     TLEFT,TRIGHT,LINTERFACE,NAEZ,NEMB,BRAVAIS,KAOEZ)
c
c NOW !!! RBASIS are the basis vectors in units of au/alat in (xyz) 
c         reference
c
c On the inputcard RBASIS are the basis vactors in the
c units of the bravais vectors and au/alat !
c
c --------------------------------------------------------------
c
      IF (OPT('CONDUCT ')) THEN
c
c   In case of Conductivity calculation 
c  
c  
      CALL CONDINPUT(IATCONDL,IATCONDR,NCPAIRD,NCONDPAIR,NAEZ,RBASIS,
     &                      ALAT,IEGFOUT) 
      IGF = 1 
      END IF 
c
c---> Calculates the real and reciprocal lattice vectors  
c
      IF (.NOT.(LINTERFACE)) THEN
      CALL LATTICE3D(ALAT,BRAVAIS,NGMAX,NRMAX,NSHLG,NSHLR,NSG,NSR,
     &              GN,RM,NAEZ,VOLUME0)
c
c---> Calculates Ewald sums MADELSMAT to be used in VMADELBLK 
c
      CALL STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,
     &            GN,RM,RBASIS,MADELSMAT,VOLUME0)
      ELSE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Calculates 2d lattice                                        !new1
      CALL LATTICE2D(ALAT,BRAVAIS,NGMAX,NRMAX,NSHLG,NSHLR,NSG, !new1
     &                NSR,GN,RM,NAEZ,VOLUME0)                  !new1
      END IF 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     
c--->     determines the configuration of the clusters in the lattice
c     write(6,*) 'cls',(cls(i1),i1=1,ncls)

      CALL CLSGEN99 (NATYP,NAEZ,NEMB,RR,NR,RBASIS,
     &               KAOEZ,Z,CLS,NCLS,NINEQ,EQINV,
     &               NACLS,ATOM,EZOA, 
     &               NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,
     &               ZPERIGHT,TLEFT,TRIGHT,
     &               RCLS,RWS,KMT,RMT,MTFAC,RCUTZ,RCUTXY,LINTERFACE,
     &               ALAT)           ! 12.03.2000
c
c     write(6,*) 'cls',(cls(i1),i1=1,ncls)
c
       IF (OPT('full inv')) THEN
      WRITE(6,*) '>>>>>>>>>>> ================ <<<<<<<<<<<'
      WRITE(6,*) '>>    F U L L   I N V E R S I O N     <<'
      WRITE(6,*) '>>>>>>>>>>> ================ <<<<<<<<<<<'
      END IF
c     
c      SUBROUTINE MSGF
c      -------------------------------------------------------
c     |                                                       |
c     |   k k r - b a n d s t r u c t u r e - p r o g r a m   |
c     |                                                       |
c      -------------------------------------------------------
c        this program is used to perform selfconsistent,all-electron,
c        scalarrelativistic calculations of the electronic structure
c        of an ideal crystal - taking into account the exact shape of
c        each cell
c
c----------------------------------------------------------------------
c
      LMACH='AIX  '
      IWRIT = 2
c      IOWFCT = 62
      IOPREF = 12
      IOPSYS = 13
      IOTREF = 16
      IOGREF = 17
      IOWREF = 18
      IOTSYS = 21
      IOWFCT = 23
      OPEN(IOPREF,FILE=I12,STATUS='OLD',FORM='FORMATTED')
      IF (IFILE.EQ.IOPSYS) 
     +     OPEN(IFILE,FILE=I13,STATUS='OLD',FORM='FORMATTED')
      IF (IFILE.EQ.11) THEN
        OPEN(11,STATUS='OLD',FORM='FORMATTED')
      ELSE
        OPEN(11,STATUS='UNKNOWN',FORM='FORMATTED')
      END IF
c
      NSHELL(0) = 1
      LSTART = .true.

c        
       CALL GAUNT2(WG,YRG)            
       CALL GAUNT(LMAX,LPOT,WG,YRG,CLEB,LOFLM,ICLEB,IEND,JEND) 
c
c     open(89,file='lebedev',status='old')
      CALL SPHERE(LPOT,YR,WTYR,RIJ,IJEND,IJD)
c     close(89)      
c     
c ------------------------------------------------------------------------
c
      IF (OPT('EigenV  ') .OR.
     +    OPT('wfct    ') .OR.
     +    OPT('iso surf')) THEN
c
c --->  calcalate eigenvalues and eigenvectors
c
        call EIGEN1(HEAD,ITITLE,M2,MMIN,MMAX,QBOUND,
     +       E1,E3,E2,TK,NPNT2,RBASIS,
     +       IFILE,IOPREF,IPF,IPFE,IPE)
        
        GOTO 999
      END IF
c
c           -----------------------------------------
c         |   do loop over selfconsitency iterations  |
c           -----------------------------------------
c
      DO 250 ITC = 1,ITCLST         ! iteration steps

        write(6,*) '>>> NUMBER OF ITERATION : ',itc

        STIME0 = DCLOCK()
        
        IF (ITC.EQ.ITCLST .AND. ITCLST.GT.0 .AND. IPE.EQ.1) IPF = IPFE

        IN = 11
        IF (ITC.EQ.1) IN = IFILE 
        
        REWIND (IN)
        E2IN = E2

        IINFO = 0
        IF (ITC.EQ.1 .OR. TEST('pot head')) IINFO = 1
c
c --->  read real potential for core relaxation,
c       in case of shape corrections (KSHAPE.GT.0) also THETAS, LMSP,..
c
        call STARTB1(IN,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +       1,NATYP,
     +       ALATNEW,RMTNEW,RMT,    ! KMT,MTFAC,
     +       ITITLE,HFIELD,IMT,IRC,VCONST,
     +       INS,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +       IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,
     +       VBC,C,DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,
     +       R,Z,A,B,IRWS,INIPOL,IINFO)
c       constant potential shift in first iteration
        VCONST = 0.D0
c Z3
c
c ---> setup of GAUNT coefficients C(l,m;l',m';l'',m'') for all 
c      nonvanishing (l'',m'')-components of the shape functions THETAS
c
      IF (KSHAPE.NE.0) 
     +     CALL SHAPE(LPOT,NATYP,GSH,ILM,IMAXSH,LMSP,NTCELL,WG,YRG)
c
c ------------------------------------------------------------------------
        IF (NPOL.EQ.0) EFERMI = E2IN
        IF (OPT('GF-EF   ')) THEN
          E1 = E2IN
          write(6,FMT=9110) E1
        END IF
        IF (OPT('DOS-EF  ')) THEN
          E1 = E2IN
          write(6,FMT=9120) E1
        END IF
c ------------------------------------------------------------------------
c
c --->  test E2IN and EMU (EMESHT)
c
        IF (DABS(E2IN-E2).GT.1d-10.AND.NPOL.NE.0) THEN
          write(6,*) 'EFERMI(file = ',in,' ) .NE. E2(emesht)',
     +         ' .AND. NPOL.NE.0 <<<<<<'
          write(6,FMT='('' E2IN = '',f14.6,''  E2 = '',f14.6)') 
     +         E2IN,E2
c
c         DECIMATE change
          E2 = E2IN ! Do this always Changed 5.11.1999
c
          IF (ISHIFT.LT.2 .AND. (ITC.EQ.1)) THEN
            E2 = E2IN
            write(6,FMT='('' FERMI ENERGY = '',f12.5)') E2

c     
c --->      doesn't read TREFLL, GREFLL, TMATLL1 and WF's from files
c     
            SINN = 0
            RIN = 0
            write(6,FMT='('' SINN, RIN     = '',2I6)') SINN,RIN
          END IF
            
        END IF                      ! (DABS(E2IN-E2).GT.1d-10.AND.NPOL.NE.0)
c        
c ------------------------------------------------------------------------
        IF (KFROZN.NE.0 .OR. ITC.EQ.1) THEN
c
c --->    if no frozen core is used or in the first iteration
c         the core states are relaxed and the potential is saved
c         to file 11, 
c
          DO 310 IP = 1,NATYP
            CALL CORELB(IPF,ITCCOR,KHYP,KCOR-1,IPRCOR,IRNUMX,IP,IRM,
     +                  NSPIN,RHOC,IMT,KSHAPE,VISP,RWS,ECORE,DRDI,R,Z,A,
     +                  B,IRWS,LMXC,KFG)
 310      CONTINUE

          REWIND (11)
          CALL RITES(11,1,NATYP,NSPIN,Z,ALATC,RMT,RMTNEW,RWS,
     +               ITITLE,R,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,
     +               LPOT,VINS,QBOUND,IRC,KSHAPE,E2IN,VBC,ECORE,
     +               LCORE,NCORE)
          REWIND (11)
          IF (OPT('GENPOT  ')) THEN
             rewind(3)
             CALL GENERALPOT (3,1,NATYP,NSPIN,Z,ALATC,RMT,RMTNEW,RWS,
     +            ITITLE,R,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,
     +            LPOT,VINS,QBOUND,IRC,KSHAPE,E2IN,VBC,ECORE,
     +            LCORE,NCORE)      
             rewind(3)
          END IF
c          CLOSE (11)
c          OPEN (11,status='old')
          IFILE = 11
          KHFELD = 0

        END IF                      ! (KFROZN.NE.0 .OR. ITC.EQ.1)
c ------------------------------------------------------------------------
c
c --->  creates a mesh in the complex energy plane
c
c implemented semicore contur
c 13.4.2000 Holger Hoehler

c implemented gap contur
c 4.10.2000 Holger Hoehler

c change here 13.4.2000

        
        CALL EMESH_SEMI(EZ,DEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,
     $                  npntsc,ebotsc,EUPSC,npntim,eimsc,EBOTGP,EUPGP,
     $                  EIMGP,NPNTGP,NPNTIMGP)


c ------------------------------------------------------------------------
c
c --->  E-mesh for special purposes
c
        IF (TEST('IE=1    ')) THEN
          IELAST=1
          write(6,*) 'IELAST is changed to 1. <<<<<<<<<<<<<<<<<<<<'
        ENDIF
c
c ------------------------------------------------------------------------
        WRITE (6,FMT=9200) IELAST
        WRITE (6,FMT=9210) TK
        WRITE (6,FMT=9220) NPOL,NPNT1,NPNT2,NPNT3

        IF (IELAST.GT.IEMXD) then
          write(6,*) 'Please, change the parameter iemxd in',
     *         ' inc.fi to a value greater equal ',ielast
          CALL RCSTOP('IEMXD   ')
        ENDIF
c
c
c ---> scales the weight WEZ for the energy integration and determine
c      number of necessary k-mesh sets
c
        DO IE = 1,IELAST
          IF (TEST('e-net   ')) 
     +         WRITE (6,FMT='(2(f14.5,f10.5))') EZ(IE),DEZ(IE)
          WEZ(IE) = -2.d0/PI*DEZ(IE)
        END DO
c ---------------------------------------------------29.10.99
c
c Write Header of Decimation file   ----------
c ------------------------------------------------------------------------
c 
        IF (OPT('deci-out').AND.(ITC.EQ.1)) THEN
c     write label of file
           OPEN(37,file='decifile',status='unknown')
           WRITE(37,*) 'INVERSE T-MATRIX AND CMOMS'
           WRITE(37,1100)
           WRITE(37,1110) ALAT,NSPIN,NAEZ,LMMAX,INS
           WRITE(37,1120) BRAVAIS
           WRITE(37,1125)
           DO IH=1,NAEZ
              WRITE(37,1130) (RBASIS(I,NAEZ),I=1,3)
           END DO
           WRITE(37,1140) E2,TK
           WRITE(37,1150) NPNT1,NPNT2,NPNT3,NPOL
        END IF 
 1100 FORMAT(' Vectors in lattice constant units'/
     &       '                                 ')
 1110 FORMAT('ALAT=',F9.6,' NSPIN=',I2,'  NAEZ=',I3,' LMMAX=',I3,
     &       ' INS=',I1)
 1120 FORMAT('BRAVAIS '/3F8.4/3F8.4/3F8.4)
 1125 FORMAT('RBASIS')
 1130 FORMAT(3F8.4)
 1140 FORMAT('EF=',F10.6,' TEMP=',F10.4,' Kelvin')
 1150 FORMAT('N1=',I3,' N2=',I3,' N3=',I3,' NPOL=',I3)
      IF (OPT('DECIMATE')) THEN

c     IENERGY = 0 READ ONLY THE BANNER
         
         CALL DECIMAREAD(EZ,TK,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,
     &        LEFTTINVLL,RIGHTTINVLL,vacflag,0,
     &        NLBASIS,NRBASIS,NAEZ,KAOEZ)    
         
      END IF
c     ------------------------------------------------------29.10.99
        MAXMESH = 1

c changed for semicore 13.04.2000 holger hoehler
c changed for gap contour 5.10.2000 holger hoehler

        IF (ITC.EQ.1 .AND. .NOT.TEST('fix mesh')) THEN
          DO IE = 1+NPNTSC+2*NPNTIM,IELAST-NPNTGP-2*NPNTIMGP
            N=1+log(DIMAG(EZ(IE))/
     +          DIMAG(EZ(IELAST-NPNTGP-2*NPNTIMGP)))/log(2.0d0)
            MAXMESH=MAX(MAXMESH,N)
          END DO
          SAVEMESH=MAXMESH
          write(6,*) 'MAXMESH : ',maxmesh
        END IF                      ! (ITC.EQ.1)
c ------------------------------------------------------------------------
c
c --->  init of radial meshs for densities
c
        CALL RINIT(IRMD*LMPOTD*NATYPD*NSPIND,RHO2NS)
        CALL RINIT(IRMD*LMPOTD*NATYPD*NSPIND,R2NEF)
        CALL CINIT(IEMXD* (LMAXD+1)*NATYPD*NSPIND,DEN)
c ------------------------------------------------------------------------
c
c --->  files for in/output of the reference and the real system
c
        IF (ITC.EQ.1.AND.RIN.NE.0) THEN
          OPEN(IOTREF,file=I16,status='OLD',FORM='unformatted')
          OPEN(IOGREF,file=I17,status='OLD',FORM='unformatted')
        ELSE
          IF (ROUT.GE.IOTREF) 
     +         OPEN(IOTREF,file=I16,status='unknown',FORM='unformatted')
          IF (ROUT.GE.IOGREF)
     +         OPEN(IOGREF,file=I17,status='unknown',FORM='unformatted')
C          IF (ROUT.GE.IOWREF)
C     +         OPEN(IOWREF,file=I18,status='unknown',FORM='unformatted')
        END IF
        IF (ITC.GT.1.AND.RIN.NE.0) THEN
          REWIND(IOTREF)
          REWIND(IOGREF)
        END IF

        IF (ITC.EQ.1.AND.SINN.NE.0) THEN
          OPEN(IOTSYS,file=I21,status='OLD',FORM='unformatted')
          OPEN(IOWFCT,file=I23,status='OLD',FORM='unformatted')
        ELSE
          IF (SOUT.GE.IOTSYS) 
     +         OPEN(IOTSYS,file=I21,status='unknown',FORM='unformatted')
          IF (SOUT.GE.IOWFCT)
     +         OPEN(IOWFCT,file=I23,status='unknown',FORM='unformatted')
        END IF
c ------------------------------------------------------------------------
c c4
        STIME = DCLOCK()
c
c      ------------------------------------------
c     |   BEGIN do loop over spins and energies  |
c      ------------------------------------------
c
c
c     added to write out the Green Function in the correct way. 
c    
        
        IF (IGF.NE.0) THEN                        ! 23.2.2000  
           write(6,*) ' green ihandle',ihandle
c          CALL FXDROPN('green','ENCODE',IHANDLE) !
        IF (OPT('CONDUCT ')) THEN
c          CALL FXDROPN('green-k-res','ENCODE',IHANDLE1)
c           write(6,*) 'OPENING green-k-res ********************' 
           END IF
        ENDIF                                     !
cccccccccccccccccccc


        CALL RINIT((LMAXD+1)*NATYPD*NSPIND,ESPV)
c
c next lines are added at 13.04.00 by Holger Hoehler for semicore

        factor=1.0
        CALL RINIT(NATYPD,SUMNS)
        CALL RINIT( (lmaxd+1)*natypd*nspind, espsc1)
        CALL RINIT( (lmaxd+1)*natypd*nspind, espsc2)

        DO 370 ISPIN = 1,NSPIN
c a2
          IF (TEST('UP      ') .AND. ISPIN.EQ.1) GOTO 370
          IF (TEST('DOWN    ') .AND. ISPIN.EQ.2) GOTO 370
c
          write(6,*) 'ISPIN = ',ispin
c
          CALL RINIT(NATYPD,SUMNS)
          CALL RINIT((LMAXD+1)*NATYPD,VASUM)
c
c
          IF (IGF.NE.0) THEN    ! added 23.2.2000      
             
c            CALL FXDRINT(IHANDLE,IELAST,1)
choshino
             write(91) IELAST
c            CALL FXDRINT(IHANDLE,NSPIN,1)
choshino
             write(91) NSPIN
c            CALL FXDRDBL(IHANDLE,TK*KB,1)
choshino
             write(91) TK*KB
c            CALL FXDRDBL(IHANDLE,E2,1) 
choshino
             write(91) E2     
choshino gf_indata
      if(ispin.eq.1) then
	write(23) ielast,nspin,numdiff,lmaxd
	do 2119 ix=1,numdiff
	write(23) rx(ix),ry(ix),rz(ix),typ1(ix),typ2(ix)
 2119 continue
      endif
c
             IF (OPT('CONDUCT ')) THEN
c            CALL FXDRINT(IHANDLE1,NSPIN,1)
c            CALL FXDRDBL(IHANDLE1,E2,1)
             END IF
          ENDIF    !  IF (IGF.NE.0)    ! up to here 23.2.2000 
c ------------------------------------------------------------------------
          IF (NPOL.EQ.0 .OR. TEST('DOS     ') 
     &       .OR. OPT('QDOS    ') .OR. OPT('QDOSEF  ') ) THEN
c
c--->       store densities of state
c
            EFCTOR = 1.0D0
            IF (TEST('EV      ')) EFCTOR = 13.6058D0
            N1 = 1
            N2 = NATYP
            IF (OPT('DOS-EF  ')) THEN
              N1 = 0
              N2 = 0
            END IF
            DO 270 I1 = N1,N2
              IPOT = NSPIN* (I1-1) + ISPIN
              WRITE (70+I1,FMT=9000) (ITITLE(IA,IPOT),IA=1,19)
              WRITE (70+I1,FMT=9005) I1
              WRITE (70+I1,FMT=9010) ISPIN,IELAST,E1,E2,EFERMI,EFCTOR
              WRITE (70+I1,FMT=9012) EFERMI
              WRITE (70+I1,FMT=9015) TK,PI*KB*TK,LATT,ALATC,
     +             INTERVX,INTERVY,INTERVZ,NACLS(1)
 270        CONTINUE                ! I1 = 1,NATYP
            
            DO I=0,LMAX
              DOSTOT(I) =0.0D0
            END DO
          END IF                    ! (NPOL.EQ.0 .OR. TEST('DOS     '))
c ------------------------------------------------------------------------

          DO 360 IE = 1,IELAST
c ------------------------------------------------------------------------
            IF(TEST('ie      ')) 
     +           write(6,9170) IE,EZ(IE)
            IF(TEST('ONE-E   ').AND.IE.gt.3) GOTO 360  
            
            E = EZ(IE)
            DF = WEZ(IE)/REAL(NSPIN)

            NMESH = 1
            IF (.NOT.TEST('fix mesh'))
c changed for gap contour 5.10.2000 holger hoehler 
     +           NMESH = 1+
     +           log(DIMAG(EZ(IE))/
     +           DIMAG(EZ(IELAST-NPNTGP-2*NPNTIMGP)))/log(2.0d0)

c added 13.4.2000 Holger Hoehler

            write(60,*) NMESH,IE
            IF(NMESH.gt.SAVEMESH) NMESH=SAVEMESH   
           

            IF ((ITC.EQ.1 .OR. OPT('rigid-ef')) .AND.
     +           RIN.NE.0) THEN
c ------------------------------------------------------------------------
c
c--->         read t - matrices and GF of reference system
c
              DO 440 I1 = 1,NREF
                IF(TEST('flow    ')) write(6,*) 'IREFPOT = ',i1
                CALL TMREAD(TREFLL(1,1,I1),IOTREF)

                IF (TEST('tmat    ')) THEN
                  WRITE(6,*) 'TREFLL (',I1,' ) :'
                  DO LM1 = 1,12
                    WRITE(6,*) TREFLL(LM1,LM1,I1)
                  END DO
                  write(6,*) 
                END IF

 440          CONTINUE

              DO 450 I1=1,NCLS
                IF(TEST('flow    ')) write(6,*) 'ICLS = ',i1
                CALL GREFREAD(NACLS(I1),GINP(1,1,I1),IOGREF)

                IF (TEST('grefll  ')) THEN
                  write(6,*) 'GINP ( 2, 2 ',I1,' ):'
                  write(6,9050) (GINP((N-1)*LMMAXD+2,2,I1),N=1,NACLSD)
                END IF

 450          CONTINUE

            ELSE                    ! (ITC.EQ.1.AND.RIN.NE.0)
c ------------------------------------------------------------------------
c     
c        -----------------------------------------------
c       |                 B E G I N                     |
c       |    t(ll'), wavefunctions and G(ll')           |
c       |    of the reference system (on clusters)      |
c        -----------------------------------------------
c     
              IF(TEST('flow    ')) 
     +             write(6,*) ' >>> calculate t matrices and G(ll)',
     +             ' of reference system (on clusters)'
              
              STIMEREF = DCLOCK()
              REWIND (IOPREF)
              
              IINFO = 0
              IF (IE.EQ.1 .AND. (ITC.EQ.1 .OR. TEST('pot head'))) 
     +             IINFO = 1
              IINFO1 = IINFO
              IF (OPT('T_hcore ')) IINFO1 = 0
c
c --->        read potential card of the reference system
c
              call STARTB1(IOPREF,IPF,IPFE,IPE,KVREL,KWS,0,LMAX,
     +             1,NREF,
     +             ALATNEW,RMTNEW,RMT, ! KMT,MTFAC,
     +             ITITLE,HFIELD,IMT,IRC,VCONST,
     +             INSREF,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELLR,
     +             IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,
     +             VBC,C,DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,
     +             R,Z,A,B,IRWS,INIPOL,IINFO1)
                
              E = EZ(IE)
                
              IF(TEST('flow    ')) write(6,*) 'rmt(1) : ',rmt(1)

c     
c--->         calculate wavefunctions and t - matrix of ref. system
c     
              DO 20 I1 = 1,NREF
                IPOT = NSPIN* (I1-1) + ISPIN

                IF (OPT('T_hcore ')) THEN 
                  F = 1.D0
                  IF (TEST('spec RMT')) F = MTFAC(I1)
                  IF (IINFO.NE.0) write(6,FMT=9150) F,RMT(I1)
                  CALL THCORE(E,RMT(I1),TREFLL(1,1,I1),F)
                ELSE
                IF (OPT('FINBASIS')) THEN
c added on 17.07.2000
                  IF (NREF.NE.1) STOP 'finbasis'
                  CALL TMATRX(ALPHA(0,I1),DETALF(I1,IE,ISPIN),
     +               ARREF(1,1),
     +               CRREF(1,1),DRDI(1,I1),E,ICST,INSREF,LMAX,
     +               PZREF(1,0),QZREF(1,0),FZREF(1,0),SZREF(1,0),
     +               PNSREF(1,1,IRMIND,1),QNSREF(1,1,IRMIND,1),
     +               TREFLL(1,1,I1),R(1,I1),VINS(IRMIND,1,IPOT),
     +               VISP(1,IPOT),Z(I1),IRWS(I1),IPAN(I1),
     +               IRCUT(0,I1),IRMIN(I1),KVREL,C,DROR(1,I1),
     +               RS(1,0,I1),S(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND)
c
                ELSE
                  CALL TMATRX(ALPHA(0,I1),DETALF(I1,IE,ISPIN),
     +               AR(1,1,I1),
     +               CR(1,1,I1),DRDI(1,I1),E,ICST,INSREF,LMAX,
     +               PZ(1,0,I1),QZ(1,0,I1),FZ(1,0,I1),SZ(1,0,I1),
     +               PNS(1,1,IRMIND,1,I1),QNS(1,1,IRMIND,1,I1),
     +               TREFLL(1,1,I1),R(1,I1),VINS(IRMIND,1,IPOT),
     +               VISP(1,IPOT),Z(I1),IRWS(I1),IPAN(I1),
     +               IRCUT(0,I1),IRMIN(I1),KVREL,C,DROR(1,I1),
     +               RS(1,0,I1),S(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND)
                END IF
                END IF

c ------------------------------------------------------------------------
                IF (TEST('trefll  ')) THEN
                  WRITE(6,*) 'TREFLL (',I1,' ) :'
c                  DO LM1 = 1,lmmax
c                    write(14,FMT='(1p,2(d14.4,d12.4))') 
c     +                   TREFLL(LM1,LM1,I1),1.d0/TREFLL(LM1,LM1,I1)
c                  END DO
c                  DO L = 0,lmax
c                    LM1=(L+1)*(L+1)
c                    write(14,FMT='(1p,2(d14.4,d12.4))') 
c     +                   TREFLL(LM1,LM1,I1),1.d0/TREFLL(LM1,LM1,I1)
c                  END DO
c                  DO LM1 = 1,lMmaxd
c                    write(6,FMT='(1p,2(d14.4,d12.4))') 
c     +                   TREFLL(LM1,LM1,I1),1.d0/TREFLL(LM1,LM1,I1)
c                  END DO
c                  write(6,*) 
c
c tmat(l) =   aj(l+1,aR)j(l,bR) - bj(l,aR)j(l+1,bR)
c           - ------------------------------------
c             j(l,bR)h(l+1,aR)a - bh(l,aR)j(l+1,bR)
c
c     a = sqrt(E),  b = sqrt(E-Vo)
c
c
       VREP = 4.D0
       
       LCALL=.true.
       rmt1 = 2.25383111d0
       A1 = SQRT(E)*RMT1
       B1 = SQRT(E - VREP)*RMT1
       CALL BESSEL(BESSJW1,BESSYW1,HANKWS1,A1,LMAXD+1,LMAX+1,.true.,
     &      .true.,.true.,LCALL)
       CALL BESSEL(BESSJW2,BESSYW2,HANKWS2,B1,LMAXD+1,LMAX+1,.true.,
     &      .true.,.true.,LCALL)
       DO L=0,LMAX
          a1 = SQRT(E)*       BESSJW1(L+1)* BESSJW2(L  ) - 
     &         SQRT(E - VREP)*BESSJW1(L  )* BESSJW2(L+1)

          b1 = SQRT(E)*       HANKWS1(L+1)* BESSJW2(L  ) -
     &         SQRT(E - VREP)*HANKWS1(L  )* BESSJW2(L+1)
          
          TMATANAL(L) =  - 1.d0/SQRT(E)*a1/b1
          
       end do
       do l=0,lmax
          lm1 = l*(l+1) + 1  
          write(6,FMT='(1p,4(d12.4,d12.4))')
     +         TREFLL(LM1,LM1,I1),tmatanal(l),1.d0/TREFLL(LM1,LM1,I1)
       end do
       do l=0,lmax
          do i=-l,l
             lm1= L*(L+1) +i + 1
             TREFLL(LM1,LM1,I1) = tmatanal(l)
          end do
       end do
c     
c                  
c                  stop 'option : trefll'
                END IF
c ------------------------------------------------------------------------
c ------------------------------------------------------------------------
                IF (TEST('Tphase  ')) THEN
                  CALL DELTAMAT(TREFLL(1,1,I1),DELTALL(1,1,I1),E)
                  write(14,FMT='(f14.5,$)') DREAL(E)
                  write(6,FMT='(f14.5,$)') DREAL(E)
                  DO L = 0,lmax
                    LM1=(L+1)*(L+1)
                    DELTA=1.D0/PI*ZLOG(DELTALL(LM1,LM1,I1))
                    write(14,FMT='(1p,d14.5,$)') DIMAG(DELTA)
                    write(6,FMT='(1p,d14.5,$)') DIMAG(DELTA)
                  END DO
                  write(14,*)
                  write(6,*)
                  goto 360
                END IF
c ------------------------------------------------------------------------

                IF(TEST('flow    ')) 
     +               write(6,*) 'tll(ref),  i1 = ',i1
                if (ROUT.ge.IOTREF)
     +               CALL TMWRIT(TREFLL(1,1,I1),IOTREF)

 20           CONTINUE              ! I1 = 1,NREF

              IF(TEST('flow    ')) write(6,*) 't-mat(Ref) o.k.'
             
C      write(6,*) ' cls=',(cls(i1),i1=1,ncls)
              DO 90 ICLS=1,NCLS
                
c
c --->          search for atom at center
c     
                I1 = 1
                IC = 0
                DO WHILE (IC.EQ.0 .AND. I1.LE.NINEQ)
C          write(6,*)  ic,i1,nineq,icls,ncls,cls(i1)
                  IF (CLS(I1).EQ.ICLS) IC = I1
                  I1 = I1 + 1
                END DO
c         write(6,*)  ' after ',ic,i1,nineq,icls,ncls,cls(i1)
                IF (IC.EQ.0) STOP 'Error in CLS(*) array in main'
                IF(TEST('flow    ')) 
     +               write(6,*) 'CLUSTER ',ICLS,' at ATOM ',IC
c
c --->          calculate GF of reference system on clusters
c
                IF (ROUT.GE.IOGREF) THEN
c
c --->            write GF to file IOGREF
c
                  WRITE (IOGREF) NACLS(ICLS)
                  CALL GLL95(LSTART,E,CLEB(1,2),ICLEB,LOFLM,IEND,
     +                 TREFLL,ATOM(1,IC),KAOEZ,REFPOT,
     +                 RCLS(1,1,ICLS),NACLS(ICLS),ALATC,IOGREF,
     +                 GINP(1,1,ICLS))
                ELSE
                  CALL GLL95(LSTART,E,CLEB(1,2),ICLEB,LOFLM,IEND,
     +                 TREFLL,ATOM(1,IC),KAOEZ,REFPOT,
     +                 RCLS(1,1,ICLS),NACLS(ICLS),ALATC,0,
     +                 GINP(1,1,ICLS))
                !  do i=1,16
                !   do j=1,16
                !  write(6,*) 'ref g',i,j,ginp(5*16+i,j,1) 
                !  end do
                ! end do
                END IF

 90           CONTINUE              ! ICLS=1,NCLS
c             write(6,*) ' after do 90' 
              IF(TEST('flow    ')) 
     +             write(6,*)  'G(n,lm,n,lm) (Ref) o.k.'
c     
c            -----------------------------------------------
c           |    t(ll'), G(ll') and wavefunctions           |
c           |    of the reference system (on clusters)      |
c           |                 E N D                         |
c            -----------------------------------------------
c     
              IF (IE.EQ.1 .AND. ITC.EQ.1)
     +        WRITE (6,FMT=9080) 
     +        DCLOCK()-STIMEREF
c ------------------------------------------------------------------------
            END IF                  ! (ITC.EQ.1.AND.RIN.NE.0)

c ------------------------------------------------------------------------
            IF (OPT('DOS REF ')) THEN
c
c --->        only valid for one type of atom in unit cell
c
              DO LM1=1,LMMAXD
                DO LM2=1,LMMAXD
                  GMATLL(LM1,LM2,1) =  GINP(LM1,LM2,1)
                END DO
              END DO
c
c --->        jump over k-Loop
c
              GOTO 990
            END IF
c ------------------------------------------------------------------------

            IF (ITC.EQ.1.AND.SINN.NE.0) THEN
c ------------------------------------------------------------------------
C B1
c
c--->       read t - matrices and WF of real system
c           write(6,*)  ' before do 350'
c
              DO 350 I1 = 1,NATYP
                IF(TEST('flow    ')) write(6,*) 'I1 = ',i1
                IPOT = NSPIN* (I1-1) + ISPIN
                CALL TMREAD(TMATLL1(1,1,I1),IOTSYS)
                
                IF (TEST('tmat    ')) THEN
                  WRITE(6,*) 'TMATLL1 (',I1,' ) :'
                  DO LM1 = 1,6
                    WRITE(6,*) TMATLL1(LM1,LM1,I1)
                  END DO
                  write(6,*) 
                END IF

                CALL WFREAD(PNS(1,1,IRMIND,1,I1),QNS(1,1,IRMIND,1,I1),
     +               ALPHA(0,I1),DETALF(I1,IE,ISPIN),AR(1,1,I1),
     +               CR(1,1,I1),PZ(1,0,I1),QZ(1,0,I1),FZ(1,0,I1),
     +               SZ(1,0,I1),IOWFCT,INS)
 350          CONTINUE
c ------------------------------------------------------------------------
            ELSE                    ! (ITC.EQ.1.AND.SINN.NE.0)
c ------------------------------------------------------------------------
C B2
c     
c            -----------------------------------------------
c           |                 B E G I N                     |
c           |           t(ll') and wavefunctions            |
c           |              of the real system               |
c            -----------------------------------------------
c     
              IF(TEST('flow    ')) 
     +             write(6,*) ' >>> calculate t matrices and wfct',
     +             ' of real system'
              
              IN = 11
              IF (ITC.EQ.1) IN = IFILE 
              REWIND (IN)

              IINFO = 0
c              IF (IE.EQ.1 .AND. (ITC.EQ.1 .OR. TEST('pot head')))
              IF (IE.EQ.1 .AND. TEST('pot head') )
     +             IINFO = 1
              
c             write(6,*) ' before startb1(3)',in
              call STARTB1(IN,IPF,IPFE,IPE,KVREL,KWS,KHFELD,LMAX,
     +             1,NATYP,
     +             ALATNEW,RMTNEW,RMT, ! KMT,MTFAC,
     +             ITITLE,HFIELD,IMT,IRC,VCONST,
     +             INS,IRNS,LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,
     +             IRCUT,IPAN,THETAS,IFUNM,NFU,LLMSP,LMSP,E2IN,
     +             VBC,C,DROR,RS,S,VISP,RWS,ECORE,LCORE,NCORE,DRDI,
     +             R,Z,A,B,IRWS,INIPOL,IINFO)
                   
C     B3
              E = EZ(IE)
c     
c--->         calculate wavefunctions and t - matrix of real system
c  
c     write(6,*) ' before 120'
              DO 120 I1 = 1,NATYP
                
                 
                IPOT = NSPIN* (I1-1) + ISPIN
                CALL TMATRX(ALPHA(0,I1),DETALF(I1,IE,ISPIN),
     +               AR(1,1,I1),
     +               CR(1,1,I1),DRDI(1,I1),E,ICST,INS,LMAX,
     +               PZ(1,0,I1),QZ(1,0,I1),FZ(1,0,I1),SZ(1,0,I1),
     +               PNS(1,1,IRMIND,1,I1),QNS(1,1,IRMIND,1,I1),
     +               TMATLL1(1,1,I1),R(1,I1),VINS(IRMIND,1,IPOT),
     +               VISP(1,IPOT),Z(I1),IRWS(I1),IPAN(I1),
     +               IRCUT(0,I1),IRMIN(I1),KVREL,C,DROR(1,I1),
     +               RS(1,0,I1),S(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND)
                
c               write(6,*) ' after tmatrx' 
c
c added on 17.07.2000
c
              IF (OPT('FINBASIS')) THEN
            !  if (ins.gt.0) then
            !     write(6,*) 'full potential version not ready !'
            !     STOP
             ! end if
c
c Map the l=3 components of wf and t-mat of the repulsive potential to  
c the real potential 
c
              write(6,*)
              write(6,*) ' SCREENING L= 4 COMPONENTS OF REAL SYSTEM '
              write(6,*)
              write(6,*) ' %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-' 
              do lm1=17,lmmaxd
                  do lm2=17,lmmaxd
                    ar(lm1,lm2,i1) = arref(lm1,lm2)
                    cr(lm1,lm2,i1) = crref(lm1,lm2)
                    tmatll1(lm1,lm2,i1) = trefll(lm1,lm2,1)
                    if (lm1.eq.lm2) 
     &          tmatll1(lm1,lm2,i1) = trefll(lm1,lm2,1) + 1.d-10
                    do i=IRMIND,IRMD
                     pns(lm1,lm2,i,1,i1) = pnsref(lm1,lm2,i,1)
                     qns(lm1,lm2,i,1,i1) = qnsref(lm1,lm2,i,1)
                    end do
                  end do
              end do
              do l=4,lmaxd
                 do i=1,irmd
                    pz(i,l,i1) = pzref(i,l)
                    qz(i,l,i1) = qzref(i,l)
                 end do
              end do  
                   
              END IF
c
                
c ------------------------------------------------------------------------
                IF (TEST('tmat    ')) THEN
                  WRITE(6,*) 'TMATLL1 (',I1,' ) :'
                  do lm1=1,lmmaxd
                    WRITE(6,FMT='(1p,(2d18.8,4x,2d18.8))')
     +                   TMATLL1(LM1,LM1,I1),1.d0/TMATLL1(LM1,LM1,I1)
                  end do
                  write(6,*) 
                END IF

                IF (TEST('setTzero')) THEN
                  do lm1=1,lmmaxd
                    TMATLL1(LM1,LM1,I1) = czero
                  end do
                END IF
c ------------------------------------------------------------------------

                IF(TEST('flow    ')) 
     +               write(6,*) 'tll(sys),  i1 = ',i1
                
                IF (SOUT.GE.IOTSYS) 
     +               CALL TMWRIT(TMATLL1(1,1,I1),IOTSYS)

                IF(TEST('flow    ')) 
     +               write(6,*) 'wfct(sys), i1 = ',i1

                IF (SOUT.ge.IOWFCT)
     +               CALL WFWRIT(
     +                 PNS(1,1,IRMIND,1,I1),QNS(1,1,IRMIND,1,I1),
     +                 ALPHA(0,I1),DETALF(I1,IE,ISPIN),AR(1,1,I1),
     +                 CR(1,1,I1),PZ(1,0,I1),QZ(1,0,I1),FZ(1,0,I1),
     +                 SZ(1,0,I1),
     +                 IOWFCT,
     +                 INS)

 120          CONTINUE              ! I1 = 1,NATYP

              IF(TEST('flow    ')) 
     +             write(6,*) 't-mat and WF of real system o.k.'
c
c          -----------------------------------------------
c         |    t(ll') and wavefunctions                   |
c         |    of the real system                         |
c         |                 E N D                         |
c          -----------------------------------------------
c     
c ------------------------------------------------------------------------
            END IF                  ! (ITC.EQ.1.AND.SINN.NE.0)

            IF (TEST('LLOYD   ')) THEN
              DO 460 I1=1,NATYP
                WRITE (74,FMT=*) 'LLOYD ALPHA',
     +               (-DIMAG(LOG(ALPHA(L,I1))),L=0,LMAXD),
     +               -DIMAG(LOG(DETALF(I1,IE,ISPIN))),'FOR IE=',IE
 460          CONTINUE
            END IF                  ! (TEST('LLOYD   '))

c ------------------------------------------------------------------------
c
c --->      determine tmat(sys) - tmat(ref) = TMATLL1 - TMATLL
c
            DO 430 I1 = 1,NAEZ
              
              INV = EQINV(I1)
              
              CFAC = CONE
              IF (I1.NE.INV) CFAC = CONEM

              INV = KAOEZ(INV)
              RF  = REFPOT(INV)

c              WRITE(6,*) 'CFAC( ',I1,' ) = ',CFAC

              DO 340 LM1 = 1,LMMAXD
                DO 330 LM2 = 1,LMMAXD
                  CFACL = CFAC**(LOFLM(LM1)+LOFLM(LM2))
                  TMATLL(LM2,LM1,I1) = 
     +                 CFACL*(TMATLL1(LM2,LM1,INV) - TREFLL(LM2,LM1,RF))
 330            CONTINUE
 340          CONTINUE
              
              IF (TEST('tmat    ')) THEN
                WRITE(6,*) 'DELTA_TMATLL (',I1,' )'
                WRITE(6,fmt='(1p,2d18.8)') 
     +               (TMATLL(LM1,LM1,I1),LM1=1,LMMAX)
                write(6,*) 
              END IF
                
 430        CONTINUE                ! I1 = 1,NAEZ
c ------------------------------------------------------------------------
c
c --->     Fourier transformation of GF(ref. system)
c          solution of DYSON equation
c          BZ integration
c
            STIMEK = DCLOCK()
c
            IF (OPT('DECIMATE')) THEN
c 
c we are in an   ie,nspin loop                        -1
c In the case of decimation read left and right host t  -matrices
c for each energy and spin . From files 37,38
c  
             CALL DECIMAREAD(EZ,TK,NPNT1,NPNT2,NPNT3,NPOL,ISPIN,
     &                      LEFTTINVLL,RIGHTTINVLL,vacflag,IE,
     &                      NLBASIS,NRBASIS,NAEZ,KAOEZ)               
 
            END IF
c
c Write out Cluster GF        added 23.2.2000
c
            IF (IGF.NE.0) THEN
               
               IF (KVREL.GE.1) THEN
                  EK = SQRT(E+E*E/(C*C))
               ELSE
                  EK = SQRT(E)
               END IF
c              CALL FXDRDBL(IHANDLE,E,2)
choshino
               write(91) E
	write(6,*) ' e=',e
c              CALL FXDRDBL(IHANDLE,EK,2)
choshino       
                write(91) EK
	WRITE(6,*) ' ek=',ek
c              CALL FXDRDBL(IHANDLE,DF*DBLE(NSPIN),2)
choshino
                write(91) DF*DBLE(NSPIN)
c    gf-indata
           write(23) E,DF*DBLE(NSPIN)
      WRITE(6,*) ' df*DBLE(NSPIN)',DF*DBLE(NSPIN)
               IF (OPT('CONDUCT ')) THEN 
                  IF (IEGFOUT.GT.IELAST) THEN
                    WRITE(6,*) 'I cannot write energy :',iegfout
                    write(6,*) 'Program is stoping '
                    STOP
                  END IF
                  IF (IE.EQ.IEGFOUT) THEN
                     WRITE(6,*) 'WROTE K-RES ENERGIES'
c                    CALL FXDRINT(IHANDLE1,IE,1)
c                    CALL FXDRINT(IHANDLE1,ISPIN,1)
c                    CALL FXDRDBL(IHANDLE1,E,2)
c                    CALL FXDRDBL(IHANDLE1,EK,2)
c                    CALL FXDRDBL(IHANDLE1,DF*DBLE(NSPIN),2)
                     write(6,*) ie,ispin,e,ek,DF*DBLE(NSPIN)
                  END IF
               END IF
              
            ENDIF               ! added 23.2.2000
ccc
            IF (OPT('CONDUCT ').AND.(IE.NE.IEGFOUT).AND.(IE.NE.1))
     &                                       GOTO 990
            IF (OPT('QDOSEF  ')) THEN
               CALL IoInput('IEGFOUT   ' ,UIO,1,7,IER)
               READ (UNIT=UIO,FMT=*) IEGFOUT
               WRITE(6,*) 'Writing QDOS for energy ....',IEGFOUT
               IF ((IE.NE.IEGFOUT).AND.(IE.NE.1))   GOTO 990 
            END IF 
c
c           write(6,*) ' before kloopz1 ?'  
            CALL KLOOPZ1(E,DF,GMATLL,TMATLL,DELTALL,
     +           INS,ALAT,LMAX,LMMAX,MD,E2,NSPIN,MAXMESH,NMESH,
     +           IE,IELAST,LSTART,IGF,CLEB(1,2),ICLEB,LOFLM,IEND,
     +           KSCOEF,NSHELL,INTERVY,INTERVZ,
     +           BLATC/ALATC,CLATC/ALATC,NAEZ,NATYP,
     +           CLS,EQINV,NACLS,RR,
     +           RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
     +           BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
     &           LEFTTINVLL,RIGHTTINVLL,vacflag,NLBASIS,NRBASIS,
     &           IHANDLE,IHANDLE1,ATOMIMP,IATCONDL,IATCONDR,NCONDPAIR,
     &           IEGFOUT,GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM)
c
            IF (IE.EQ.1 .AND. ITC.EQ.1)
     +         WRITE (6,FMT=9070) DCLOCK()-STIMEK
            
c ------------------------------------------------------------------------

 990        CONTINUE     ! k-loop jump ends here
            IF (OPT('QDOSEF  ').AND.(IE.NE.IEGFOUT).AND.(IE.NE.1))
     &                                       GOTO 360
            
            NKDOS = 1
            IF (OPT('QDOS    ').or.OPT('QDOSEF  ')) NKDOS = NQDOSKP

            DO 888 IK = 1, NKDOS    ! k-lool for k-dos calculation
               IF (OPT('QDOS    ').or.OPT('QDOSEF  ')) THEN
                  if (nqdoskp.gt.qdoskdim) then
                    write(6,*) 'Increase qdoskdim to....', nqdoskp
                     STOP
                  end if
                  IF (INS.GT.0) THEN
c     
c     --->  init of radial meshs for densities
c     
                     CALL RINIT(IRMD*LMPOTD*NATYPD*NSPIND,RHO2NS)
                     CALL RINIT(IRMD*LMPOTD*NATYPD*NSPIND,R2NEF)
                     CALL CINIT(IEMXD* (LMAXD+1)*NATYPD*NSPIND,DEN)
                     
                  END IF
                  DO IATYP=1,NATYP                 
                     DO LM1=1,LMMAX
                        DO LM2=1,LMMAX
                       GMATLL(LM1,LM2,IATYP) = GSQDOS(LM1,LM2,IATYP,IK) 
                        END DO
                     END DO
                  END DO
               
               END IF 
                IF (OPT('CONDUCT ').AND.(IE.EQ.IEGFOUT)) THEN
c
c     Write out the wavefunction for conductance calculation
c     
          CALL CONDUCTWF(E,EK,IE,ISPIN,IRMIN,IRWS,IRNS,LMAX,NSPIN,NATYP,
     &               AR,PNS,PZ,R,KVREL,IATCONDL,IATCONDR,NCONDPAIR,INS)
                END IF
                
               IF (INS.GE.1) THEN
c     
c---  >         in case of non-spherical host potential
c     
              
              CALL RHONSB(DEN,DF,DRDI,GMATLL,E,IE,IELAST,ISPIN,IRMIN,
     +             IRWS,LMAX,NSPIN,1,NATYP,RHO2NS,R2NEF,IPAN,
     +             IRCUT,THETAS,NTCELL,IFUNM,LMSP,KVREL,AR,CR,
     +             PNS,QNS,C,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,IEND,
     +             JEND,DENEF)
              
           ELSE                 ! (INS.GE.1)
c     
c--->         in case of spherical host potential
c
c ------------------------------------------------------------------------
              IF (TEST('gmatll  ')) THEN
                do IATYP=1,NATYP
                  write(6,FMT='(''GMATLL  atom '',I4)') 
     +                 IATYP
                  do lm1=1,lmmax
                    do LM=1,lmmax
c                      if (abs(GMATLL(LM1,LM,IATYP)).gt.1.d-12)
                          write(6,FMT='( 1p,2I4,d22.12,d22.12 )') 
     +                     lm1,lm,GMATLL(LM1,LM,IATYP)
                    end do
                  end do
                end do
              END IF
c ------------------------------------------------------------------------
              CALL RHOLMB(DEN,DF,GMATLL,E,IE,IELAST,IPF,ISPIN,IRWS,
     +             KVREL,LMAX,NSPIN,NATYP,RHO2NS,R2NEF,DRDI,IPAN,
     +             IRCUT,THETAS,NTCELL,IFUNM,LMSP,SUMNS,VASUM,
     +             DENEF,C,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB,IEND,JEND,
     +             NSHELL)


            END IF                  ! (INS.GE.1)

c     Integrate DOS over semicore energy contour
c     change 13.4.2000 Holger Hoehler

                   
             CALL SEMICORE(NATYP,LMAX,SUMDF,IE,npntsc,npntim,
     +                    ISPIN,IPAN,IRCUT,NTCELL,rho2ns,NFU,
     +                    llmsp,THETAS,EZ,DRDI,ielast,DF,DEN,
     +                    factor,espco2)

CSEMI            IF (KTE.EQ.1) CALL ESPVB(DEN,DF,E,ESPV,IE,LMAX,ISPIN,NATYP)
            
c ------------------------------------------------------------------------
            IF (NPOL.EQ.0 .OR. TEST('DOS     ') 
     &          .OR. OPT('QDOS    ') .OR. OPT('QDOSEF  ') ) THEN
c
c--->         store densities of state
c
              DF = WEZ(IE)/REAL(NSPIN)
              SIGN = 1.0D0
              IF (ISPIN.NE.NSPIN) SIGN = -1.0D0
              EFCTOR = 1.0D0
              IF (TEST('EV      ')) EFCTOR = 13.6058D0
              
              DO 280 I1 = 1,NATYP
                IPOT = NSPIN* (I1-1) + ISPIN
c
c --->          summation of l ( and E )
c
                DOS = 0.0D0
                DO L=0,LMAX
                  DOS = DOS -1.0D0/PI*DIMAG(DEN(IE,L,I1,ISPIN))
                  DOSTOT(L) = DOSTOT(L) + DIMAG(DF*DEN(IE,L,I1,ISPIN))
                END DO
c     
                E3 = REAL(EZ(IE))*EFCTOR
                IF (TEST('sub ef  ')) E3 = E3 - EFERMI*EFCTOR 
c     
                
                IF (OPT('QDOS    ').or.OPT('QDOSEF  ')) THEN
                  WRITE (70+I1,FMT=9020)E3,QDOSKP(1,IK),QDOSKP(2,IK),
     &                  DOS*SIGN/EFCTOR
                ELSE            ! OPT('QDOS    ')
                   IF (OPT('DOS-EF  ')) THEN
                      WRITE (70,FMT=9021) I1,
     +                (-1.0D0/PI*DIMAG(DEN(IE,L,I1,ISPIN))*sign/EFCTOR,
     +                     L=0,LMAX),
     +                     DOS*SIGN/EFCTOR
                   ELSE         ! (OPT('DOS-EF  '))
                      WRITE (70+I1,FMT=9020) E3,
     +                (-1.0D0/PI*DIMAG(DEN(IE,L,I1,ISPIN))*sign/EFCTOR,
     +                     L=0,LMAX),
     +                     DOS*SIGN/EFCTOR
                   END IF       ! (OPT('DOS-EF  '))
                END IF          ! OPT('QDOS    ')
                
C     WRITE (50+I1,FMT=9020) REAL(EZ(IE)),
C     +             (-1.0D0/PI*REAL(DEN(IE,L,I1,ISPIN)),L=0,LMAX)

                IF (IE.EQ.IELAST ) THEN
                IF ((.NOT.(OPT('QDOS    '))).AND.
     &              (.NOT.(OPT('QDOSEF  ')))) THEN
                  WRITE (70+I1,FMT=9025) (DOSTOT(L)/EFCTOR,L=0,LMAX)
                  IF (ISPIN.NE.NSPIN) WRITE(70+I1,FMT=9026)
                END IF
        END IF      
 280          CONTINUE              ! I1 = 1,NATYP
              IF (ISPIN.NE.NSPIN .AND. OPT('DOS-EF  ')) 
     +             WRITE(70,FMT=9026)
            END IF                  ! (NPOL.EQ.0 .OR. TEST('DOS     '))
            IF (OPT('QDOS    ')) THEN
               IF (IE.EQ.1.AND.(ISPIN.EQ.1)) THEN
               WRITE (70+NATYP+2,FMT=9020)  QDOSKP(1,IK), QDOSKP(2,IK)
               END IF
            END IF
 888        CONTINUE     ! q-dos loop terminates here
c ------------------------------------------------------------------------
            IF (OPT('QDOS    ')) THEN
               IF (ISPIN.EQ.1) THEN
               WRITE (70+NATYP+1,FMT=9020) E3
               END IF
            END IF
 360      CONTINUE                  ! IE = 1,IELAST

c     Calculate single particle energies.
c     change 13.4.2000 Holger Hoehler
c begin

        If ( kte.eq.1 ) Then
          If ( lmach.eq.'INTEL' .or. lmach.eq.'PARAG' ) Then
            Stop 'No parallel machine in this version of ELOOP'
          ElseIf (lmach.eq.'C64  '.or.lmach.eq.'Schne'.or.
     +            lmach.eq.'CASIO') Then
            Stop 'Not used any more -- to slow, 
     +           it is faster to do it by hand'
          End If
          If ( ispin .eq. 1 ) espcor = 0.d0

c     Calculate semicore contribution to the single particle energies
          Do ie = 1,  (npntsc + 2*npntim)
            e  =  ez(ie)
            df = wez(ie)/real(nspin)
            Call espscb(den,df,e,espsc1,ie,lmax,ispin,natyp,factor)
            Call espscb(den,df,e,espsc2,ie,lmax,ispin,natyp,1.d0)
          End Do
c     Calculate valence contribution to the single particle energies
          Do ie =  (npntsc + 2*npntim) + 1, ielast
            e  =  ez(ie)
            df = wez(ie)/real(nspin)
            Call espvb(den,df,e,espv,ie,lmax,ispin,natyp)
          End Do

          Do i1 = 1, natyp
            Do lm1 = 0, lmax
              espv(lm1,i1,ispin) = espv(lm1,i1,ispin) +
     +             espsc2(lm1,i1,ispin)
c     espcor is the correction to the semicore energies
c     which should be added to the total single particle energy
              espcor = espcor +
     $             ( espsc1(lm1,i1,ispin) - espsc2(lm1,i1,ispin) )
            End Do
          End Do
        End If
c end

 370    CONTINUE                    !  ISPIN = 1,NSPIN
        IF (OPT('QDOS    ').or.OPT('QDOSEF  ')) THEN
           do i1=1,natyp
              CLOSE(70+i1)
           end do
           WRITE(6,*) 'Q-DOS Calculation finished, program is STOPING'
           STOP
        END IF
        
c
c    ------------------------------------------
c   |  END of do loop over spins and energies  |
c    ------------------------------------------
c
        ETIME = DCLOCK()

C C2
        IF (ROUT.GE.IOTREF) CLOSE(IOTREF)
        IF (ROUT.GE.IOGREF) CLOSE(IOGREF)
C        IF (ROUT.GE.IOWREF) CLOSE(IOWREF)
        IF (SOUT.GE.IOTSYS) CLOSE(IOTSYS)
        IF (SOUT.GE.IOWFCT) CLOSE(IOWFCT)
        


        IF (NSPIN.EQ.2) THEN
          IDIM = IRMD*LMPOTD*NATYPD
C          CALL DCOPY(IDIM,RHO2NS(1,1,1,NSPIN),1,WORK,1)
          CALL DAXPY(IDIM,ONEM,RHO2NS(1,1,1,1),1,RHO2NS(1,1,1,NSPIN),1)
          CALL DSCAL(IDIM,TWO,RHO2NS(1,1,1,1),1)
          CALL DAXPY(IDIM,ONE,RHO2NS(1,1,1,NSPIN),1,RHO2NS(1,1,1,1),1)
C
C          CALL DCOPY(IDIM,R2NEF(1,1,1,NSPIN),1,WORK,1)
          CALL DAXPY(IDIM,ONEM,R2NEF(1,1,1,1),1,R2NEF(1,1,1,NSPIN),1)
          CALL DSCAL(IDIM,TWO,R2NEF(1,1,1,1),1)
          CALL DAXPY(IDIM,ONE,R2NEF(1,1,1,NSPIN),1,R2NEF(1,1,1,1),1)
        END IF

c
c--->   determine total charge density expanded in spherical harmonics
c
        IF(TEST('flow    ')) write(6,*) '>>> RHOTOTB'
        
        CALL RHOTOTB(IPF,NATYP,NSPIN,RHO2NS,RHOC,Z,DRDI,IRWS,IRCUT,
     +       LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,CHRGNT,
     +       ITC,NSHELL)

        IF(TEST('flow    ')) write(6,*) '<<< RHOTOTB'
c ------------------------------------------------------------------------
        IF (TEST('rhons   ')) then
          OPEN(75,file='rho2ns',status='unknown')
          DO IATYP=1,NATYP
            DO LM=1,9
              write(75,9180) 
     +             IATYP,LM,
     +             (R(I,IATYP),RHO2NS(I,LM,IATYP,1)*R(I,IATYP)**2,
     +             I = 2,IRWS(IATYP))
            END DO
          END DO
          CLOSE(75)
        END IF                      ! (TEST('rhons   '))
c ------------------------------------------------------------------------
c
c --->  determine new Fermi level due to valence charge up to 
c       old Fermi level E2 and density of states DENEF
c


        IF ( ITC .GT. 1 .AND.
     +       CHRGNT*CHRGOLD .LT. 0.D0 .AND.
     +       ABS(CHRGNT) .GT. 5.D-2) THEN
          E2SHIFT = CHRGNT/(CHRGNT-CHRGOLD)*(E2-EOLD(1))
        ELSE
          E2SHIFT = CHRGNT/DENEF
        END IF

        E2SHIFT = DMIN1(DABS(E2SHIFT),0.05D0)*DSIGN(1.0D0,E2SHIFT)
        EOLD(1) = E2
        EOLD(2) = E2
        DENOLD = DENEF
        CHRGOLD = CHRGNT

        IF (ISHIFT.LT.2) E2 = E2 - E2SHIFT

        WRITE (6,FMT=9030) ITC,E2,DENEF/REAL(NSPIN*NAEZ)

        DF = 2.0D0/PI*E2SHIFT/REAL(NSPIN)

        DO 420 ISPIN = 1,NSPIN

          IF (KTE.EQ.1) THEN
            
            DO 380 I1 = 1,NATYP
              ESPV(0,I1,ISPIN) = ESPV(0,I1,ISPIN) -
     +                           EOLD(ISPIN)*CHRGNT/REAL(NSPIN*NAEZ)
 380        CONTINUE
          END IF                    ! (KTE.EQ.1)
c
c--->     get correct density
c
          IF (.NOT.(OPT('DECIMATE'))) THEN  ! NEW for decimation!!!!
          DO 410 I1 = 1,NATYP
            DO 400 LM = 1,LMPOT
              DO 390 I = 1,IRC(I1)
                RHO2NS(I,LM,I1,ISPIN) = RHO2NS(I,LM,I1,ISPIN) +
     +                                  REAL(DF)*R2NEF(I,LM,I1,ISPIN)
 390          CONTINUE
 400        CONTINUE
 410      CONTINUE
          END IF

 420    CONTINUE                    ! ISPIN = 1,NSPIN
     
c
c--->   potential part
c
        IF (LRHOSYM)
     +       CALL RHOSYMM(LMPOT,NSPIN,1,NATYP,RHO2NS,IXIPOL,
     +                    IRWS,IRCUT,IPAN,KSHAPE)

       CALL VINTRAS(CMOM,CMINST,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,R,DRDI,
     +            IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM,IFUNM,IMAXSH,GSH,
     +            THETAS,LMSP)

        IF (.NOT.TEST('NoMadel ')) THEN
!     +  CALL VMADEL1(AVMAD0,BVMAD0,CMOM,CMINST,LPOT,NSPIN,NATYP,VONS,
!     +       Z,R,IRWS,IRCUT,IPAN,KSHAPE,NAEZ,EQINV,KAOEZ,
!     +      IMAD,JMAD,NVMAD)
C
c     +   CALL VMADEL(AVMAD,BVMAD,CMOM,CMINST,LPOT,NSPIN,NATYP,VONS,
c     +       Z,R,IRWS,IRCUT,IPAN,KSHAPE)
c new
c                                    ---------------- 29.10.99
c Write the CMOMS to a file...
c
           IF (OPT('deci-out').AND.(ITC.EQ.1)) THEN
              WRITE(37,1080) NAEZ,LMPOT
              DO IH=1,NAEZ
                 WRITE(37,*) IH
                 DO LM=1,LMPOT
                    C00(LM) = CMOM(LM,IH)
                    IF (INS.NE.0) C00(LM) = C00(LM)+CMINST(LM,IH)
                 END DO
                 C00(1) = C00(1) - Z(IH)/rfpi
                 WRITE(37,1090) (C00(lm),LM=1,LMPOT)
                 
              END DO  
 1080         FORMAT('CMOMC',2I6)
 1090         FORMAT(4D22.14)
              
           END IF               ! OPT('deci-out')
c   Read the host CMOM from a file !
           IF (OPT('DECIMATE').AND.(ITC.EQ.1)) THEN
c I read this only once and store , since I read all the 
c t-matrices the file is in the correct place.
c DO NOT forget to rewind it before the iteration
c   
             CALL CMOMSREAD(NLBASIS,NRBASIS,NAEZ,CMOMHOST,VACFLAG,KAOEZ)
              
           END IF               ! OPT('DECIMATE')
c                                         ---------------- 29.10.99
           if (.NOT.(LINTERFACE)) THEN  
              
              CALL VMADELBLK(CMOM,CMINST,MADELSMAT,LPOT,NSPIN,NATYP,
     &             VONS,Z,R,IRWS,IRCUT,IPAN,KSHAPE,YRG,WG)
           else
c     
              CALL VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NAEZ,VONS, ! new1
     &             Z,R,IRWS,IRCUT,IPAN,KSHAPE,ALAT, ! new1
     &             BRAVAIS,NGMAX,NRMAX,NSHLG,NSHLR, ! new1
     &             NSG,NSR,GN,RM,VOLUME0, ! new1
     &             RBASIS,NLBASIS,NLEFT, ZPERLEFT,TLEFT, ! new1
     &             NRBASIS,NRIGHT,ZPERIGHT,TRIGHT, ! new1
     &             CMOMHOST,CHRGNT,    ! new1
     &             YRG,WG,ICC,ATOMIMP)      ! new1
           end if               ! LATTICE3d
        END IF
c     end new
     
c     
c     Force calculation 18.5.2000
c
        KFORCE = 0
        IF (INS.GT.0) THEN      
        CALL IoInput('KFORCE   ',UIO,1,7,IER)
                        READ (UNIT=UIO,FMT=*) KFORCE
        if (Ier.gt.0) then 
         write(6,*) 'Please include :  KFORCE= 0   in inputcard'
        end if
        END IF
        IF (KFORCE.EQ.1) THEN
           IF (INS.EQ.0) THEN
              CALL FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,
     &                     R,DRDI,IRWS,Z)
              CALL FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     &                     DRDI,IRWS)
             ELSE
                CALL FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,
     +                     R,DRDI,IMT,Z)
                CALL FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                     DRDI,IMT)
              END IF
          END IF
c
c Force Calculation stops here look after vxclm
c
c
        IF (KTE.EQ.1) THEN
          CALL ESPCB(ESPC,NSPIN,NATYP,ECORE,LCORE,NCORE)

          CALL EPOTINB(EPOTIN,NSPIN,NATYP,RHO2NS,VISP,R,DRDI,INS,IRMIN,
     +                 IRWS,LPOT,VINS,IRCUT,IPAN,Z)

          CALL ECOUB(CMOM,ECOU,LPOT,NSPIN,NATYP,RHO2NS,VONS,Z,R,DRDI,
     +               IRWS,KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,IFUNM,ILM,
     +               NTCELL,GSH,THETAS)
        END IF                      ! (KTE.EQ.1)
c
c
c       igga=1  kxc < 0
        if(kxc.ge.0) then
        write(6,*) ' lsda'
c
c       lsda
c
        CALL VXCLM(EXC,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,R,DRDI,
     +             IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM,IMAXSH,IFUNM,
     +             THETAS,YR,WTYR,IJEND,WORK,LMSP)
        else
c
c          pw91-gga (kexcor=kxc < -1)
        write(6,*)  ' pw91-gga'
c
        CALL VXCGGA(EXC,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,R,DRDI,
     +      A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM,IMAXSH,IFUNM,
     +             THETAS,YR,WTYR,IJEND,WORK,LMSP)
        endif
c
c Force calculation start here 18.5.2000
c
        IF (KFORCE.EQ.1) THEN
          IF (KSHAPE.EQ.0) THEN
            CALL FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                  ALAT,DRDI,IRWS,0)
  
          ELSE
            CALL FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                  ALAT,DRDI,IMT,0)
          END IF
  
        END IF
c
c Force calculation ends
c
        CALL MTZERO(LMPOT,NATYP,NSPIN,VONS,VBC,Z,R,DRDI,IMT,IRCUT,IPAN,
     +              NTCELL,LMSP,IFUNM,THETAS,IRWS,E2SHIFT,ISHIFT,
     +              NSHELL,LINTERFACE)
c
c--->   convolute potential with shape function for next iteration
c
        IF (KSHAPE.NE.0) THEN
          DO 200 ISPIN = 1,NSPIN
            DO 190 I1 = 1,NATYP
              IPOT = NSPIN* (I1-1) + ISPIN
              CALL CONVOL(IRCUT(1,I1),IRC(I1),NTCELL(I1),IMAXSH(LMPOT),
     +                    ILM,IFUNM,LMPOT,GSH,THETAS,Z(I1),RFPI,R(1,I1),
     +                    VONS(1,1,IPOT),LMSP)
  190       CONTINUE
  200     CONTINUE
        END IF                      ! (KSHAPE.NE.0)
c
c--->   final construction of the potentials (straight mixing)
c
        MIX = MIXING
        IF (TEST('alt mix ')) MIX = MIXING/REAL(1+MOD(ITC,2))
        IF (TEST('spec mix')) 
     +       MIX = MIXING/
     +       (1.0D0 + 1.0D+3 * ABS(CHRGNT)/REAL(NAEZ*NSPIN))

        CALL MIXSTR(RMSAVQ,RMSAVM,INS,LPOT,LMPOT,0,NSHELL,
     +              1,NATYP,NSPIN,
     +              ITC,RFPI,FPI,IPF,
     +              MIX,
     +              FCM,IRC,IRMIN,R,DRDI,VONS,
     +              VISP,VINS)

        WRITE(6,FMT=9160) ITC,MIX

        IF (ITC.EQ.1) RMSAV0 = 1.0d2*MAX(RMSAVQ,RMSAVM)
        IF (MAX(RMSAVQ,RMSAVM).GE.QBOUND) THEN
c
c --->    potential mixing procedures
c
          IF ( IMIX.EQ.1 
     +         .AND. ABS(CHRGNT)/REAL(NAEZ*NSPIN).LE.1.D-3 
     +         .AND. .NOT. (TEST('spec mix').OR.TEST('alt mix ')) )
c
c --->      switch to Chebycheff acceleration scheme
c
     +         IMIX = 2
          
c
c---->    Chebycheff mixing scheme
c
          IF ( IMIX.EQ.2 )
     +         CALL CHEACC(VISP,VINS,VONS,MIX,1,NATYP,NSPIN,INS,IRMIN,
     +                     IRC,IPF,LMPOT,ITC)
c
c---->    Broyden or Andersen updating schemes
c
          IF (IMIX.GE.3) THEN

            CALL BRYDBM(VISP,VONS,VINS,INS,LMPOT,R,DRDI,MIX,ATWGHT,
     +           IRC,IRMIN,NSPIN,1,NATYP,ITDBRY,IMIX,20,IPF)

c ------------------------------------------------------------------------
c            IF ((LMACH.EQ.'INTEL') .AND. MYNODE().NE.0) THEN
c              CALL CRECV(4,VONS,IRMD*LMPOTD*NPOTD*8)
c            ELSE
c              CALL BRYDBM(VISP,VONS,VINS,INS,LMPOT,R,DRDI,MIX,ATWGHT,
c     +                    IRC,IRMIN,NSPIN,1,NATYP,ITDBRY,IMIX,20,IPF)
c              IF ((LMACH.EQ.'INTEL') .AND. MYNODE().EQ.0) THEN
c                CALL CSEND(4,VONS,IRMD*LMPOTD*NPOTD*8,-1,0)
c              END IF
c            END IF
c ------------------------------------------------------------------------

          END IF                    ! (IMIX.GE.3)
c
c---->    reset to start new iteration
c
          DO 230 I = 1,NSPIN*NATYP

            IF (NSPIN.EQ.2) THEN
              IT = (I+1)/2
            ELSE
              IT = I
            END IF

            IRC1 = IRC(IT)
            CALL DCOPY(IRC1,VONS(1,1,I),1,VISP(1,I),1)

C              DO 440 J = 1,IRC1
C                 visp(J,I) = vons(J,1,I)
C 440          CONTINUE

            IF (INS.NE.0 .AND. LPOT.GT.0) THEN
              IRMIN1 = IRMIN(IT)
              DO 220 LM = 2,LMPOT
                DO 210 J = IRMIN1,IRC1
                  VINS(J,LM,I) = VONS(J,LM,I)
  210           CONTINUE
  220         CONTINUE
            END IF                  ! (INS.NE.0 .AND. LPOT.GT.0)

  230     CONTINUE                  ! I = 1,NSPIN*NATYP

        END IF                      ! (MAX(RMSAVQ,RMSAVM).GE.QBOUND)

        REWIND 11
        IF (IPOTOU.GT.0) THEN
          IF ((LMACH.EQ.'INTEL') .AND. MYNODE().NE.0) GO TO 240
          ENEW(1) = E2
          ENEW(2) = E2
          IF (OPT('rigid-ef').OR.OPT('DECIMATE')) THEN
            ENEW(1) = EOLD(1)
            ENEW(2) = EOLD(2)
          END IF
          CALL RITES(11,1,NATYP,NSPIN,Z,ALATC,RMT,RMTNEW,RWS,ITITLE,R,
     +               DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,LPOT,VINS,
     +               QBOUND,IRC,KSHAPE,ENEW,VBC,ECORE,LCORE,NCORE)
          CLOSE (11)
          OPEN (11,status='old')
        END IF
         IF (OPT('GENPOT  ')) THEN
         rewind (3)
      CALL GENERALPOT (3,1,NATYP,NSPIN,Z,ALATC,RMT,RMTNEW,RWS,
     +                 ITITLE,R,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE)      
         CLOSE(3)
         open(3,status='old')
      END IF

  240   CONTINUE
c
        IF (KTE.EQ.1 .AND. ICC.EQ.0) THEN
          CALL ETOTB1(ECOU,E2,EPOTIN,ESPC,ESPV,EXC,KPRE,LMAX,LPOT,NSPIN,
     +               NATYP,NSHELL(1),espcor,espco2,npntsc)
        END IF

        ETIME0 = DCLOCK()
      
        WRITE (6,FMT=9090) 
     +       ETIME - STIME
        WRITE (6,FMT=9100)
     +       ITC,ETIME0 - STIME0

        IF (MAX(RMSAVQ,RMSAVM).LT.QBOUND) THEN
          write(6,*) 'ITERATION FINISHED +++'
          GO TO 260
        END IF

        IF (MAX(RMSAVQ,RMSAVM).GT.RMSAV0) THEN
          write(6,*) 'ITERATION DIVERGED ---'
          GO TO 260
        END IF

c ------------------------------------------------------------------------
        NMIX = 6
        IF (TEST('opt mix ') .AND. ITC.LE.NMIX*ITDBRY+1) THEN
          IF (ITC.GT.2) THEN
            NITMIX = NITMIX + 1
            MIXV(NITMIX) = MAX(RMSAVQ,RMSAVM)/RMSAV1
            MIXSAV(NITMIX) = MIX
            MIXW(NITMIX) = 1.D0
            write(6,FMT='('' OPT MIX'',3F15.8)') 
     +           MIXSAV(NITMIX),MIXV(NITMIX),MIXW(NITMIX)
          END IF
c
c --->    change MIXING for determination of optimal mixing parameter
c
          IF (ITC.EQ.MAX(1,ITC/ITDBRY)*ITDBRY+1) 
     +         MIXING = MIXING*MAX(HFIELD,1.d-4)
          IF (ITC.EQ.MAX(1,ITC/ITDBRY)*ITDBRY+2) 
     +         NITMIX = NITMIX - 1
          IF (ITC.EQ.NMIX*ITDBRY+1) THEN
            write(6,9140) (I,MIXSAV(i),MIXV(i),MIXW(i),i=1,NITMIX)
            CALL REG2(NITMIX,MIXSAV,MIXV,MIXW,MIXC,XX,FM)
            write(6,FMT='('' XX,FM         :'',1p,2d10.2)') XX,1.d0-FM
            write(6,FMT='('' C(1),C(2),C(3):'',1p,3d10.2)') 
     +           MIXC(1),MIXC(2),MIXC(3)
            MIXING = MIXING0
            IF (XX .GT. 1.D-2*MIXING0 .AND.
     +          XX .LT. 1.D+2*MIXING0 .AND.
     +          FM .LT. 1.d0 .AND. FM .GT. 0.d0 ) THEN
              MIXING = XX
            ELSE
              write(6,*) 'Optimization of MIXING factor not successful.'
              write(6,*) 'mix(opt), conv(opt) :',XX,1.d0-FM
              STOP 'OPT MIX'
            END IF
            write(6,FMT=9130) MIXING
            IMIX = 1
          END IF
          RMSAV1 = MAX(RMSAVQ,RMSAVM)
        END IF
c ------------------------------------------------------------------------
      CLOSE(37)
      CLOSE(38) 
  250 CONTINUE                      ! ITC = 1,ITCLST
c
c         -------------------------------------------------
c        |  END of do loop over selfconsitency iterations  |
c         -------------------------------------------------
c
 260  CONTINUE                      ! jump mark
C
C ---> close k-mesh file
C
      CLOSE(11)
      CLOSE(52)
C
c ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2025 FORMAT((3F15.8,I6))
 9000 FORMAT ('#',19a4)
 9005 FORMAT ('# I1    :',I8)
 9010 FORMAT ('# ISPIN :',I8,  '   IELAST :',I5,/,
     +        '# E1,E2 :',2f12.5,' EFERMI :',f12.5,'   EFCTR',f10.6)
 9012 FORMAT ('# FERMI :',f12.5)
 9015 FORMAT ('# TK    =',f8.1,'   Kelvin =',3p,f8.3,' mRyd',0p,/
     +        '# LATT  :',I8,  '   ALAT   :',f12.5,/,
     +        '# INTERV X,Y,Z  :',3I5,/,
     +        '# NACLS :',I8)
 9020 FORMAT (1p,7d11.3)
 9021 FORMAT (1p,I6,7d11.3)
 9025 FORMAT ('# Integrated DOS ',1p,d10.3,7d11.3)
 9026 FORMAT ('&')
 9030 FORMAT (I4,' E FERMI ',F12.6,10x,'DEN OF ST ',F12.6)
 9040 format ('ie,e =',I4,'(',2F10.4,' )')
 9050 format (2f12.6)
 9060 format (1p,d14.3,d11.3)

 9070 format (' TIME IN KLOOPZ1           : ',f9.2)
 9080 format (' TIME IN GLL95             : ',f9.2)
 9090 format (' time in loop 370          : ',f9.2)
 9100 format (' TIME IN ITERATION     ',i3,' : ',f9.2)

 9110 FORMAT(' Determination of Cluster GF at EF =',f12.6)
 9120 FORMAT(' Determination of DOS  at EF =',f12.6)
 9130 FORMAT(' OPT MIX, MIXING set to ',1p,d12.3)
 9140 FORMAT(i5,3f12.6)
 9150 FORMAT(' HARDCORE ref. sys. : R_CORE=',f12.4,'*',f12.4)
 9160 FORMAT(I4,' Mixing Factor used : ',1P,D12.2)
 9170 FORMAT(' IE    = ',I3,', E = (',1p,d11.3,' ,',d11.3,')')
 9180 FORMAT('& RHO2NS, atom ',I6,' LM =',I6,/,(f15.6,d15.4))
 9200 FORMAT(' The calculations are done with',I4,' points.')
 9210 FORMAT(' The temperature used is',F14.6,' K.')
 9220 FORMAT(I4,' poles and',I4,' and',I4,' and',I4,
     +     ' points on the contours are used.')
 999  continue

      write(6,*) 'End of KKR'

      END
