C ************************************************************************
      SUBROUTINE EIGEN0(ISPIN,Q,E,LSYMIN,IO)
      implicit none
C ************************************************************************
c
c     determination of eigenvalues of the kkr-matrix
c
c     input : IO = 0 count negative EV of KKR matrix
c                  1                   of TINVLL
c                  2 WF coefficients at given q-vector and energy
c
c      p. zahn, febr. 96
C ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
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
      INTEGER ALM,NLM,NDIM
      PARAMETER (ALM  = NDIMGK*LMAXSQ,
     +           NDIM = NPRINCD*LMAXSQ,
     +           NLM  = NAEZD*LMAXSQ)
c      INTEGER NCLEB
c      PARAMETER (NCLEB = LMPOTD*LMMAXD)
      INTEGER NAUX,NAUXSYM
      PARAMETER (NAUX    = MAX(2*ALM**2+5*ALM,2*NDIM**2+5*NDIM),
     +           NAUXSYM = 4*ALM)
      INTEGER NAUXD
      PARAMETER(NAUXD = NDIM*NDIM*8)
      DOUBLE COMPLEX CONE,CZERO,CONEM,CI,CTWO
      PARAMETER (
     +     CONE  = ( 1.0D0,0.0D0),
     +     CTWO  = ( 2.0D0,0.0D0),
     +     CZERO = ( 0.0D0,0.0D0),
     +     CONEM = (-1.0D0,0.0D0),
     +     CI    = ( 0.0D0,1.0D0))
      DOUBLE PRECISION PI,TPI,ZERO
      PARAMETER (
     +     PI    = 3.14159265358979312D0,
     +     TPI   = 2.0d0 * PI,
     +     ZERO  = 0.0D0 )
c ------------------------------------------------------------------------
C     ..
C     .. Scalar Arguments ..
      INTEGER 
     +     IO,                      ! option   0 : eigenvalues KKR matrix
                                    !          1 : eigenvalues TINVLL
                                    !          2 : eigenvector of KKR matrix
     +     ISPIN

      DOUBLE COMPLEX E
      LOGICAL LSYMIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  Q(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION  CONVPU,DUMMY,QDECI,WMIN,RTEMP,FTEMP
      DOUBLE COMPLEX 
     +     CFAC,CFACL,CTEMP,CTEMP1,C1,
     +     DET,
     +     ELAST,
     +     FACI,
     +     GTEMP,
     +     T
      INTEGER
     +     BASIS,
     +     DIM,
     +     I,IDUMMY,I1,IE,ILM1,IMAX,IMIN,IN,INFO,INV,
     +     IOPT,IPOT,IRF,ISPINLAST,
     +     J,JLM1,JLM2,JN,JRF,
     +     L,LM,LM1,LM2,LMJ,
     +     NL,NUMT,
     +     POS,POST,
     +     RF,ICLS,II,IC,ILM,N,IILM1
      LOGICAL LT
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX 
     +     ALPHA(0:LMAXD),
     +     AR(LMMAXD,LMMAXD),
     +     CR(LMMAXD,LMMAXD),
     +     DETALF(NATYPD),
     +     DETALF1(NREFD),
     +     FZ(IRMD,0:LMAXD),
     +     PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +     PZ(IRMD,0:LMAXD),
     +     QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +     QZ(IRMD,0:LMAXD),
     +     SZ(IRMD,0:LMAXD)
c
      DOUBLE COMPLEX 
     +     AUX(NAUX),
     +     AUXD(NAUXD),
c     +     AUXN(NAUXN),
     +     CC(NLM),
     +     CCSYM(ALM,ALM),
     +     DMAT(NDIM,NDIM),
     +     EMAT(NDIM,NDIM),
     +     GLLKES(ALM,ALM),
c     +     GLLKE2S(NDIM,NDIM,NLAYERD),
     +     GLLKESYM(ALM*(ALM+1)/2),
     +     TAU(NDIM,NDIM),
     +     TINVLL(LMAXSQ,LMAXSQ,NAEZD),
     +     TLL(LMAXSQ,LMAXSQ),
     +     W(NLM),
     +     ZTAU(MAX(ALM,NDIM))
c
      DOUBLE PRECISION 
     +     AUXSYM(NAUXSYM),
c     +     CCSYM(ALM,ALM),
c     +     GLLKESYM(ALM*(ALM+1)/2),
     +     WSYM(NLM)
      INTEGER 
     +     IPIV(NDIM),
     +     IPVT(LMAXSQ),
     +     IPVT1(NLM),
     +     IND(NLM),
     +     IFAIL(MAX(ALM,NDIM))

      DOUBLE COMPLEX EXPIDL
      LOGICAL 
     +     LSTART,LSYM,
     +     OPT,
     +     TEST
      LOGICAL SELECT(NLM)

c     .. external statement
      EXTERNAL 
     +     OPT,TEST,ZGEEVY,
     +     EXPIDL,DELTAMAT,PZNORM,TMATRX,ZSORT,ZNORM
      INTRINSIC DABS,DCONJG,DIMAG,DREAL,LOG,ABS

c     .. arrays in common
c ------------------------------------------------------------------------
      INTEGER     IUP,IDO,NUM
      DOUBLE PRECISION WUP,WDO
      COMMON /EIGENV/ WUP,WDO,NUM,IUP,IDO
c ------------------------------------------------------------------------
      DOUBLE COMPLEX CL(NAEZD*LMMAXD)
      COMMON /CL/ CL
c ------------------------------------------------------------------------
      DOUBLE COMPLEX
     +     PZSQ(LMMAXD,NATYPD)
      COMMON /PZSQ/ PZSQ
c ------------------------------------------------------------------------
      DOUBLE COMPLEX 
     +     GLLKE(ALM,ALM),
     +     GLLKE1(NDIM,NDIM,NLAYERD),
     +     GLLKE2(NDIM,NDIM,NLAYERD),
     +     GLLKE3(NDIM,NDIM,NLAYERD)
      COMMON / GLLKE / GLLKE,GLLKE1,GLLKE2,GLLKE3,CTEMP,CTEMP1
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
     +     ALATC,C,
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
     +     NSHELL(0:NSHELD),
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
     +         NREF,NATYP,NAEZ,NINEQ,NSHELL,NZ,NLAYER,
     +         IMT,IRWS,IRNS,NTCELL,NTCELLR,
     +         REFPOT,NCLS,NACLS,
     +         CLS,EZOA,ATOM,EQINV,KAOEZ,INIPOL,LATT,ICC,ICST
c ------------------------------------------------------------------------
c     .. data
      INTEGER LF(144)
      DOUBLE PRECISION BOUND 
      DATA LF /0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/
     +     BOUND  / 1.0d-10 /
     +     LSTART / .TRUE. /
c ------------------------------------------------------------------------
      IF (DABS(DREAL(E)).LT.BOUND) THEN
        write(6,*) 'E :',E
        RETURN
      END IF

      LSYM = LSYMIN
      IF (.NOT. (TEST('fix sym ') .OR. OPT('COMPLEX ')) ) 
     +     LSYM = (DABS(DIMAG(E)).LT.BOUND)
      if (test('fix true') ) LSYM = .TRUE.
      if (test('fix fals') ) LSYM = .FALSE.

      CONVPU = ALATC / TPI
      NUM  = 0
      NUMT = 0
      IUP  = 1
      IDO  = 1
      WUP = ZERO
      WDO = ZERO
      LSTART =  LSTART .OR. 
     +         (IO.EQ.1) .OR.
     +         (ISPIN.NE.ISPINLAST) .OR. 
     +         (ABS(E-ELAST).GT.BOUND)
c ------------------------------------------------------------------------
      IF (LSTART) THEN
c
c --->  t-matrices and wave functions only in case of changed energy
c       or spin
c
        IF (TEST('eigen0   ')) WRITE(6,*) 't matrices and WF''s'
c
c --->  t-matrices for reference potential
c
        DO 20 I1 = 1,NREF
          IPOT = NSPIN* (I1-1) + ISPIN
          
          IF (OPT('T_hcore ')) THEN
            FTEMP = 1.D0
            IF (TEST('spec RMT')) FTEMP = MTFAC(I1)
            IF (TEST('flow    ')) 
     +           write(6,FMT='(''R_MT_FAC,R_MT    '',f12.6,f12.4)') 
     +           FTEMP,RMT(I1)
            CALL THCORE(E,RMT(I1),TREFLL(1,1,I1),FTEMP)
          ELSE                      ! OPT('T_hcore ')
            CALL TMATRX(ALPHA(0),DETALF1(I1),
     +         AR(1,1),
     +         CR(1,1),DRDI1(1,I1),E,ICST,INSREF,LMAX,
     +         PZ(1,0),QZ(1,0),FZ(1,0),SZ(1,0),
     +         PNS(1,1,IRMIND,1),QNS(1,1,IRMIND,1),
     +         TREFLL(1,1,I1),R1(1,I1),VINS1(IRMIND,1,IPOT),
     +         VISP1(1,IPOT),Z1(I1),IRWS1(I1),IPAN1(I1),
     +         IRCUT1(0,I1),IRMIN1(I1),KVREL,C,DROR1(1,I1),
     +         RS1(1,0,I1),S1(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND)
          END IF                    ! OPT('T_hcore ')

          CALL DELTAMAT(TREFLL(1,1,I1),DELTALL(1,1,I1),E)

          IF (TEST('tmat    ')) THEN
            write(6,*) 'TREFLL(',i1,')'
            write(6,fmt='(1p,2(d14.3,d11.3))') 
     +           (TREFLL(LM1,LM1,I1),
     +           cone/TREFLL(LM1,LM1,I1),LM1=1,LMMAXD)
          END IF

 20     END DO

        IF (IO.NE.1) THEN
c
c --->    calculate GF of reference system on clusters
c
          DO 50 ICLS=1,NCLS
c
c --->      search for atom at center
c
            I1 = 1
            IC = 0
            DO WHILE (IC.EQ.0 .AND. I1.LE.NINEQ)
              IF (CLS(I1).EQ.ICLS) IC = I1
              I1 = I1 + 1
            END DO
            IF (IC.EQ.0) STOP 'Error in CLS(*) array in EIGEN0'
c
            CALL GLL95(.TRUE.,E,CLEB(1,2),ICLEB,LOFLM,IEND,
     +           TREFLL,ATOM(1,IC),KAOEZ,REFPOT,
     +           RCLS(1,1,ICLS),NACLS(ICLS),ALATC,0,
     +           GINP(1,1,ICLS))
 50       CONTINUE                  ! ICLS=1,NCLS
        END IF                      ! (IO.NE.1)
c
c --->  t-matrix and WF of real system
c
        DO 30 I1 = 1,NATYP
          IPOT = NSPIN* (I1-1) + ISPIN
          
          IF (TEST('HARDCORE')) THEN
c
c --->      calculate eigenvalues of hard sphere system
c           (only for test)
c
            FTEMP = 1.D0
            IF (TEST('spec RMT')) FTEMP = HFIELD
            IF (TEST('flow    ')) 
     +           write(6,FMT='(''R_MT_FAC,R_MT_SYS'',f12.6,f12.4)') 
     +           FTEMP,RMT(I1)
            CALL THCORE(E,RMT(I1),TMATLL1(1,1,I1),FTEMP)
          ELSE                     
            CALL TMATRX(ALPHA(0),DETALF(I1),
     +         AR(1,1),
     +         CR(1,1),DRDI(1,I1),E,ICST,INS,LMAX,
     +         PZ(1,0),QZ(1,0),FZ(1,0),SZ(1,0),
     +         PNS(1,1,IRMIND,1),QNS(1,1,IRMIND,1),
     +         TMATLL1(1,1,I1),R(1,I1),VINS(IRMIND,1,IPOT),
     +         VISP(1,IPOT),Z(I1),IRWS(I1),IPAN(I1),
     +         IRCUT(0,I1),IRMIN(I1),KVREL,C,DROR(1,I1),
     +         RS(1,0,I1),S(0,I1),CLEB(1,1),LOFLM,ICLEB,IEND)
c
c --->      calculate norm of regular radial wave functions PZ,FZ
c
            CALL PZNORM(PZSQ(1,I1),E,6,ISPIN,IRWS(I1),KVREL,
     +         LMAX,NSPIN,DRDI(1,I1),IPAN(I1),IRCUT(0,I1),THETAS,
     +         NTCELL(I1),C,PZ(1,0),FZ(1,0))
c
            IF (TEST('pznorm  ')) THEN
              WRITE(6,*) 'Norm of radial wave function :',I1
              WRITE(6,FMT='(1p,2d12.4,d16.4)') 
     +             (PZSQ(L*L,I1),ABS(PZSQ(L*L,I1)),L=1,LMAX+1)
            END IF
          END IF
c ------------------------------------------------------------------------
          IF (TEST('setTzero')) THEN
            DO LM=1,LMMAXD
              TMATLL1(LM,LM,I1) = CZERO
            END DO
          END IF
c
          IF (TEST('tmatll1 ')) THEN
            write(6,*) 'TMATLL1(',i1,')'
            write(6,fmt='(1p,2(d14.3,d11.3),0p,f14.6)') 
     +           (TMATLL1(L*L,L*L,I1),
     +           cone/TMATLL1(L*L,L*L,I1),
     +       DIMAG(1.D0/PI*ZLOG(EXPIDL(TMATLL1(L*L,L*L,I1),E,LF(L*L)))),   
     +           L=1,LMAXD+1)
          END IF
c ------------------------------------------------------------------------

 30     END DO                      ! I1 = 1,NATYP

c ------------------------------------------------------------------------
c
c --->  determine delta_tmat = TMATLL1 - TREFLL
c
        DO 40 I1 = 1,NAEZ
              
          CFAC = CONE
          IF (I1.GT.NINEQ) CFAC = CONEM
          
          INV = KAOEZ(I1)
          RF  = REFPOT(INV)

          DO 340 LM1 = 1,LMMAXD
            DO 330 LM2 = 1,LMMAXD
              CFACL = CFAC**(LF(LM1)+LF(LM2))
c              write(6,FMT='(1p,2d20.10)') cfacl
              TMATLL(LM2,LM1,I1) = 
     +             CFACL*(TMATLL1(LM2,LM1,INV) - TREFLL(LM2,LM1,RF))
 330        CONTINUE
 340      CONTINUE
              
          IF (TEST('tmat    ')) THEN
            write(6,*) 'TMATLL(',i1,') = delta_tmat'
            write(6,fmt='(1p,(d14.3,d11.3))') 
     +           (TMATLL(LM1,LM1,I1),LM1=1,LMMAXD)
          END IF
c
c --->    transformation to hermitian KKR matrix
c
          IF (INSREF.EQ.0 .AND. INS.EQ.0) THEN
            DO 470 LM1 = 1,LMMAXD
              TMATLL(LM1,LM1,I1) = TMATLL(LM1,LM1,I1)
     +             /DELTALL(LM1,LM1,RF)/DELTALL(LM1,LM1,RF)
 470        END DO
          ELSE
            DO 471 LM1 = 1,LMMAXD
              DO 472 LM2 = 1,LMMAXD
                TMATLL(LM1,LM2,I1) = TMATLL(LM1,LM2,I1)
     +               /DELTALL(LM1,LM1,RF)/DELTALL(LM2,LM2,RF)
 472          END DO
 471        END DO
          END IF
c
c --->    init TINVLL
c
          CALL CINIT(LMAXSQ*LMAXSQ,TINVLL(1,1,I1))
c
c --->    TINVLL = TMATLL**(-1)
c
          IF (INSREF.EQ.0 .AND. INS.EQ.0) THEN
c            
            DO 350 LM1 = 1,LMMAX
              TINVLL(LM1,LM1,I1) = CONVPU/TMATLL(LM1,LM1,I1)
 350        CONTINUE
c
c --->      negativ eigenvalues of TINVLL (!! only correct for INS.LT.1 )
c
            DO 170 LM1 = 1,LMMAXD
              IF (dreal(TINVLL(LM1,LM1,I1)).lt.zero) NUMT = NUMT + 1
 170        END DO

          ELSE                      ! (INS.LT.1)
            
            DO 360 LM1 = 1,LMMAX
              TINVLL(LM1,LM1,I1) = CONVPU
 360        CONTINUE
c            
            CALL ZCOPY(LMMAX*LMMAX,TMATLL(1,1,I1),1,TLL,1)
c
c --->      invert t-matrix and count negative eigenvalues (real part)
c
            CALL ZGETRF(LMMAX,LMMAX,TLL,LMAXSQ,IPVT,INFO)
            DO 171 LM = 1,LMMAXD
              RTEMP = 1.d0
              IF (IPVT(LM).NE.LM) RTEMP = - 1.d0
              IF (DREAL(1.D0/TLL(LM,LM))*RTEMP .LT. ZERO ) 
     +             NUMT = NUMT + 1
 171        END DO
            CALL ZGETRS('N',LMMAX,LMMAX,TLL,LMAXSQ,
     +           IPVT,TINVLL(1,1,I1),LMAXSQ,INFO)
          END IF                    ! (INS.LT.1)
c
          IF (TEST('tmat    ')) THEN
            write(6,*) 'TINVLL(',i1,')'
            write(6,fmt='(1p,(d14.3,d11.3))') 
     +           (TINVLL(LM1,LM1,I1),LM1=1,LMMAXD)
          END IF
          
 40     END DO                      ! I1 = 1,NAEZ

        IF (IO.EQ.1) THEN
c
c --->    only counts the negative EV of TINVLL
c
          NUM = NUMT
          IF (TEST('eigen0   '))
     +         WRITE(6,FMT='(''E,N :'',F12.4,I6)') DREAL(E),NUM
          RETURN
        END IF

      END IF                        ! (LSTART)

c ************************************************************************
      IF (OPT('full inv')) THEN
c ************************************************************************
c
c --->  Fourier transformation of GINP
c
        CALL DLKE0(GLLKE,ALATC,NAEZ,CLS,EQINV,NACLS,
     +       RR,EZOA,ATOM,Q,IDUMMY,KAOEZ,RCLS,GINP)
c
c --->  GLLKE = (i)^LF(LM1)*DELTALL*GLLKE*DELTALL*(i)^(-LF(LM2))
c
        FACI = CI
        DO 90 I = 1,NAEZ
          IRF = REFPOT(KAOEZ(I))
          DO 80 J = 1,NAEZ  
            JRF = REFPOT(KAOEZ(J))
            DO 70 LM1 = 1,LMAXSQ
              ILM1=LMAXSQ*(I-1) + LM1
              DO 60 LM2 = 1,LMAXSQ
                JLM2=LMAXSQ*(J-1) + LM2
c
                GLLKE(ILM1,JLM2) =  
     +               GLLKE(ILM1,JLM2)
     +               *(FACI**(LF(LM1)-LF(LM2)))
     +               *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
c                
 60           END DO
 70         END DO
 80       END DO
 90     END DO
c
c --->  KKR matrix GLLKE =  GLLKE - TINVLL
c       Loop over diagonal elements
c
        DO 130 J = 1,NAEZ
          LMJ = LMAXSQ*(J-1)
          DO 120 LM2 = 1,LMAXSQ
            JLM2 = LMJ + LM2
            DO 110 LM1 = 1,LMAXSQ
              JLM1 = LMJ + LM1
              GLLKE(JLM1,JLM2) = GLLKE(JLM1,JLM2)-TINVLL(LM1,LM2,J)
 110        CONTINUE
 120      CONTINUE
 130    CONTINUE
c
c --->  eigenvalues of GLLKE
c
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (.not.LSYM) THEN
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c --->    eigenvalues of general complex matrix
c
c
c --->    save GLLKE for eigen vector determination
c
          IF (IO.GT.1) 
     +         CALL ZCOPY(ALM*ALM,GLLKE,1,GLLKES,1)

c         CALL ZGEEVY(0,GLLKE,ALM,W,CC,ALM,1,SELECT,ALM,AUX,NAUX,
c    +         ZTAU,IFAIL)
c          CALL ZGEEV(0,GLLKE,ALM,W,CC,ALM,SELECT,ALM,AUX,NAUX)

          IF (IO.LT.2) THEN
c
c --->      sorting of complex eigenvalues due to real part
c
            CALL ZSORT(W,IND,ALM,POS)
c
c --->      number of eigenvalues with real part smaller then zero
c
            DET = CZERO
            DO 160 I = 1,ALM
              DET = DET + ZLOG(W(I))
              ii = ind(i)
              IF (dreal(w(II)).lt.zero) NUM = NUM + 1
              IF (TEST('www     ')) 
     +             write(6,FMT=1010) 
     +             ii,w(ii),dabs(dimag(w(ii))/dreal(w(ii)))
 160        END DO
c
            IF (TEST('www     ')) write(6,9040) DET
c  
            WDO = DREAL(W(IND(NUM)))
            WUP = DREAL(W(IND(NUM+1)))

          ELSE                      ! (IO.LT.2)
c
c --->      calculation of EIGEN VECTOR
c
c
c --->      initialization of array SELECT for ZGEEV
c
            DO 190 I = 1,ALM
              SELECT(I) = .FALSE.
 190        END DO
c
c --->      selection of the eigenvalue nearest to zero
c
            WMIN = ABS(W(1))
            IMIN = 1
            DO 180 I = 2,ALM
              IF (ABS(w(I)).lt.WMIN) THEN
                IMIN = I
                WMIN = ABS(W(I))
              END IF
 180        END DO
            SELECT(IMIN) = .TRUE.

            NUM = IMIN
            WUP = WMIN
            WDO = DREAL(W(IMIN))

c ------------------------------------------------------------------------
            if (test('w(imin) ')) 
     +           write(6,FMT=1020) W(IMIN),IMIN
c ------------------------------------------------------------------------
c
c --->      calculate eigenvalue and eigenvector selected by the
c           true element of logical array SELECT
c
c           CALL ZGEEVY(2,GLLKES,ALM,W,CC,NLM,1,SELECT,ALM,AUX,NAUX,
c    +           ZTAU,IFAIL)
c            CALL ZGEEV(2,GLLKES,ALM,W,CC,NLM,SELECT,ALM,AUX,NAUX)

c ------------------------------------------------------------------------
            if (test('w(imin) ')) write(6,FMT=1020) W(IMIN)
            IF (TEST('cc      ')) THEN
              write(6,*) 'eigenvector CC/cc(1) :'
              write(6,FMT=1000) (LM,CC(LM)/cc(1),lm=1,naez*lmmaxd)
            END IF
c ------------------------------------------------------------------------

            IF (INS.GE.1) STOP 'INS.GE.1 in EIGEN0 !!'
c
c --->      transformation of eigenvector of GLLKE to the 
c           eigenvector of the KKR matrix
c
            DO 200 I1 = 1,NAEZ
              INV = KAOEZ(I1)       ! real potential at site I1
              RF  = REFPOT(INV)     ! reference potential at site I1
              DO 210 LM1 = 1,LMMAXD
                ILM = (I1-1)*LMMAXD + LM1
                CL(ILM) = CC(ILM)
     +               *CI**LF(LM1)
     +               *TINVLL(LM1,LM1,I1)
     +               *SQRT(PZSQ(LM1,INV))
     +               /DELTALL(LM1,LM1,RF)
 210          END DO                ! LM1 = 1,LMMAXD
 200        END DO                  ! I1 = 1,NAEZ

            CALL ZNORM(CL,ALM,1.0d0) ! normalization of CL

          END IF                    ! (IO.LT.2)

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ELSE                        ! (.not.LSYM)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c --->    complex hermitian matrix GLLKE
c         save in upper-packed storage mode in GLLKESYM         
c
          J = 0
          DO 140 II = 1,NAEZ*LMAXSQ
            DO 150 I = 1,II
              J = J+1
              GLLKESYM(J)=(GLLKE(I,II)+DCONJG(GLLKE(II,I)))/CTWO
 150        END DO
 140      END DO
          
          IF (IO.LT.2) THEN

c            CALL ZHPEV(20,GLLKESYM,WSYM,DUMMY,1,ALM,AUXSYM,NAUXSYM)
            CALL ZHPEVY(20,GLLKESYM,WSYM,DUMMY,1,ALM,AUXSYM,NAUXSYM)
c
c --->      counts number of eigenvalues smaller then zero 
c           (wsym(*) are in increasing order)
c
            IF (TEST('www     ')) then 
              DET = CZERO
              write(6,FMT=1030) (ii,wsym(ii),ii=1,alm)
              DO 950 I=1,ALM
                DET = DET + LOG(WSYM(I))
 950          END DO
              write(6,9040) DREAL(DET)
            END IF
c
            WDO = ZERO
            WUP = WSYM(1)
            I   = 1
            NUM = 0
            DO 940 WHILE (WSYM(I).LT.ZERO .AND. I.LE.ALM)
              WDO = WSYM(I)
              NUM = I
              I   = I+1
              WUP = WSYM(I)
 940        END DO
            
          ELSE                      ! (IO.LT.2)
c
c --->      calculate eigenvalues and all eigenvectors of GLLKE (io=2)
c           in symmetrized form
c
c            CALL ZHPEV(21,GLLKESYM,WSYM,CCSYM,ALM,ALM,AUXSYM,NAUXSYM)
            CALL ZHPEVY(21,GLLKESYM,WSYM,CCSYM,ALM,ALM,AUXSYM,NAUXSYM)
c
c --->      determine smallest eigenvalue
c
            WMIN = DABS(WSYM(1))
            IMIN = 1
            DO 390 I = 2,ALM
              IF (DABS(WSYM(I)).lt.WMIN) THEN
                IMIN = I
                WMIN = DABS(WSYM(I))
              END IF
 390        END DO

            NUM = IMIN
            WUP = WMIN
            WDO = WSYM(IMIN)

c ------------------------------------------------------------------------
            if (test('w(imin) ')) 
     +           write(6,FMT='('' WSYM(IMIN) :'',1P,D12.4,I6)') 
     +           WSYM(IMIN),IMIN
            IF (TEST('cc      ')) THEN
              write(6,*) 'eigenvector CC/cc(1) :'
              write(6,FMT=1000) 
     +             (LM,CCSYM(LM,IMIN)/CCSYM(1,IMIN),lm=1,naez*lmmaxd)
            END IF
c ------------------------------------------------------------------------
c
c --->      back transformation of eigenvector
c
            DO 410 I1 = 1,NAEZ
              INV = KAOEZ(I1)      ! real potential at site I1
              RF  = REFPOT(INV)     ! reference potential at site I1
c              write(6,*) 'I1,INV,RF',i1,inv,rf
              DO 400 LM1 = 1,LMMAXD
                ILM = (I1-1)*LMMAXD + LM1
                CL(ILM) = CCSYM(ILM,IMIN)
     +               *CI**LF(LM1)
     +               *TINVLL(LM1,LM1,I1)
     +               *ZSQRT(PZSQ(LM1,INV))
     +               /DELTALL(LM1,LM1,RF)
 400          END DO                ! LM1 = 1,LMMAXD
 410        END DO                  ! I1 = 1,NAEZ

            CALL ZNORM(CL,ALM,1.0d0)  ! normalization of CL

          END IF                    ! (IO.LT.2)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        END IF                      ! (.not.LSYM)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c ************************************************************************
      ELSE                          ! (OPT('full inv'))
c ************************************************************************
c
c ---> TB formalism with O(N) matrix inversion
c
c ------------------------------------------------------------------------
        CALL CINIT(NDIM*NDIM*NLAYER,GLLKE1)
        CALL CINIT(NDIM*NDIM*NLAYER,GLLKE2)
        CALL CINIT(NDIM*NDIM*NLAYER,GLLKE3)
c
c --->  GLLKE1,GLLKE2,GLLKE3 = fouriertransformed GF of reference system
c
        CALL DLKEB(GLLKE1,GLLKE2,GLLKE3,ALATC,NAEZ,NZ,CLS,EQINV,
     +       NACLS,RR,EZOA,ATOM,Q,IE,KAOEZ,RCLS,GINP)

        CALL CONVGLL(NLAYER,NDIM,REFPOT,KAOEZ,DELTALL,LF)
c ------------------------------------------------------------------------
      IF (TEST('gllke   ')) THEN
        DIM = 32
        write(6,*) 'TEST GLLKE2'
        DO J=1,NLAYER
          DO LM1=1,DIM
            DO LM2=MAX(1,LM1-1),MIN(LM1+1,DIM)
C            DO LM2=1,DIM
              IF (ABS(GLLKE2(LM1,LM2,J)).GT.BOUND)      
     +             write(6,FMT='(3I6,1p,2d19.9)') 
     +             J,LM1,LM2,GLLKE2(LM1,LM2,J)
            END DO
          END DO
        END DO
      END IF                        ! (TEST('gllke   '))
c ------------------------------------------------------------------------
c
c --->  construct KKR matrix GLLKE2 =  GLLKE2 - TINVLL
c
        DO 520 J=1,NAEZ
c     
c --->    position in GLLKE2
c     
          I=(J-1)/NPRINCD+1
          II=(J-(I-1)*NPRINCD-1)*LMAXSQ
          DO 510 LM1 = 1,LMAXSQ
            IILM1 = II + LM1
            DO 500 LM2 = 1,LMAXSQ
              GLLKE2(IILM1,II+LM2,I) = 
     +             GLLKE2(IILM1,II+LM2,I) - TINVLL(LM1,LM2,J)
 500        CONTINUE
 510      CONTINUE
          
 520    CONTINUE
c ------------------------------------------------------------------------
        IF (OPT('DECI    ')) THEN
c
c --->    DECIMATION TECHNIQUE FOR SEMIINFINITE CRYSTAL
c
          QDECI = 1.d-6
          IOPT  = 3
          IMAX  = 20
          LT    = .false.
c          IF (ZABS(E).LT.1.d-1 .AND. IMAX.GE.0) IMAX = 3*IMAX
c
c         LEFT
c
          CALL DECI(GLLKE3(1,1,1),GLLKE2(1,1,1),
     +              GLLKE1(1,1,1),TAU,
     +              QDECI,NDIM,AUXD,NAUXD,IPIV,IOPT,IMAX,INFO)
          if (info.lt.0 .or. LT) 
     +      write(6,9030) 'L',info,(dreal(auxd(i)),i=1,2),info,q(1),q(2)
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,GLLKE3(1,1,1),NDIM,
     +         TAU,NDIM,CZERO,DMAT,NDIM)
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,DMAT,NDIM,
     +         GLLKE1(1,1,1),NDIM,CONE,GLLKE2(1,1,1),NDIM)
c
c         RIGHT
c
          CALL DECI(GLLKE1(1,1,NLAYER-1),GLLKE2(1,1,NLAYER),
     +              GLLKE3(1,1,NLAYER-1),TAU,
     +              QDECI,NDIM,AUXD,NAUXD,IPIV,IOPT,IMAX,INFO)
          if (info.lt.0 .or. LT) 
     +     write(6,9030) 'R',info,(dreal(auxd(i)),i=1,2),info,q(1),q(2)
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,GLLKE1(1,1,NLAYER-1),
     +         NDIM,TAU,NDIM,CZERO,DMAT,NDIM)
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,DMAT,NDIM,
     +         GLLKE3(1,1,NLAYER-1),NDIM,CONE,GLLKE2(1,1,NLAYER),NDIM)
          
          CALL CINIT(NDIM*NDIM,GLLKE1(1,1,NLAYER))
          CALL CINIT(NDIM*NDIM,GLLKE3(1,1,NLAYER))
        END IF                      ! (OPT('DECI    '))
c ------------------------------------------------------------------------
        IF (LSYM) THEN
c
c --->    symmetrize KKR matrix (= complex hermitian matrix)
c
          DO 220 NL = 1,NLAYER
            DO 230 LM1 = 1,NDIM
              DO 240 LM2 = 1,LM1
                GLLKE2(LM1,LM2,NL) =
     +               (GLLKE2(LM1,LM2,NL) + 
     +               DCONJG(GLLKE2(LM2,LM1,NL)))/CTWO
                GLLKE2(LM2,LM1,NL) =
     +               DCONJG(GLLKE2(LM1,LM2,NL))
 240          END DO
              DO 250 LM2 = 1,NDIM
                GLLKE1(LM1,LM2,NL) =
     +               (GLLKE1(LM1,LM2,NL) + 
     +               DCONJG(GLLKE3(LM2,LM1,NL)))/CTWO
                GLLKE3(LM2,LM1,NL) =
     +               DCONJG(GLLKE1(LM1,LM2,NL))
 250          END DO
 230        END DO
 220      END DO                    ! NL = 1,NLAYER
        END IF                      ! (LSYM)
c ------------------------------------------------------------------------
c
c --->  O(N) ALGORITHM FOR MATRIX INVERSION
c
        IOPT = 1
        IF (OPT('SLAB    ')) IOPT = 11          ! count neg. eigenvalues
c
        IF (IO.GT.1 .AND. IOPT.EQ.1) IOPT = 2   ! eigenvector determin.
c
        CALL INVERT(NLAYER,NDIM,GLLKE1,GLLKE2,GLLKE3,IOPT,
     +       CTEMP,CTEMP1)
c
c --->  IOPT = 1,11 : GLLKE2 contains the elements of the diagonal matrix 
c                     D**(-1) of the factorized KKR matrix 
c                     g - t**-1 = L**-1 D**-1 U**-1
c       IOPT = 2    : GLLKE2 contains D**(-1)
c                     GLLKE3          C
c
c ------------------------------------------------------------------------
c
c --->    eigenvalues of diagonal blocks
c
          DO 525 NL=1,NLAYER
c
c --->      eigenvalues of GLLKE2(NL)
c
            CALL ZCOPY(NDIM*NDIM,GLLKE2(1,1,NL),1,DMAT(1,1),1)
c            CALL ZGEEV(0,DMAT(1,1),NDIM,W((NL-1)*NDIM+1),
c     +           CC,NLM,SELECT,NDIM,AUX,NAUX)
c           CALL ZGEEVY(0,DMAT(1,1),NDIM,W((NL-1)*NDIM+1),
c    +           CC,NLM,1,SELECT,NDIM,AUX,NAUX,ZTAU,IFAIL)
 525      END DO                    ! NL=1,NLAYER

          IF (IO.LT.2) THEN
c
c --->      sorting of complex eigenvalues due to real part
c
            CALL ZSORT(W,IND,NLAYER*NDIM,POS)
c
c --->      number of eigenvalues with real part smaller then zero
c
            DET = CZERO
            DO 526 I = 1,NLAYER*NDIM
              DET = DET + ZLOG(W(I))
              ii=ind(i)
              IF (dreal(w(I)).lt.zero) NUM = NUM + 1
              IF (TEST('www     ')) 
     +             write(6,FMT=1010) 
     +             ii,w(ii),dabs(dimag(w(ii))/dreal(w(ii)))
 526        END DO
c
            IF (TEST('www     ')) write(6,9040) DET
c  
            WDO = DREAL(W(IND(NUM)))
            WUP = DREAL(W(IND(NUM+1)))

          ELSE                      ! (IO.LT.2)
c
c --->      calculation of eigen vector
c
c
c --->      selection of the eigenvalue nearest to zero
c
            WMIN = ABS(W(1))
            IMIN = 1
            IF (DREAL(W(1)).LT.ZERO) NUM = NUM + 1
            DO 528 I = 2,NLM
              IF (DREAL(W(I)).LT.ZERO) NUM = NUM + 1
              IF (ABS(W(I)).LT.WMIN) THEN
                IMIN = I
                WMIN = ABS(W(I))
              END IF
 528        END DO

            WDO = DREAL(W(IMIN))
            IF (WDO.GT.ZERO) NUM = NUM + 1

            POS = (IMIN-1)/NDIM +1
            IMIN = IMIN - (POS-1)*NDIM
c
c --->      initialization of array SELECT for ZGEEV
c
            DO 527 I = 1,NLM
              SELECT(I) = .FALSE.
 527        END DO
            SELECT(IMIN) = .TRUE.

            IF (POS.NE.NLAYER) THEN
              write(6,*) 'PROBLEM in ',
     +           'eigenvector determination. << !!'
            END IF
c ------------------------------------------------------------------------
            if (test('w(imin) ')) 
     +           write(6,
     +           FMT='(''POS,IMIN,W(IMIN) :'',2I6,1p,2d12.2 )') 
     +           POS,IMIN,W(IMIN+(POS-1)*NDIM)
c ------------------------------------------------------------------------
            CALL CINIT(NLM,CC)
c
c --->      calculate eigenvalue and eigenvector selected by the
c           true element of logical array SELECT
c
            CALL ZCOPY(NDIM*NDIM,GLLKE2(1,1,POS),1,DMAT(1,1),1)
c            CALL ZGEEV(2,DMAT(1,1),NDIM,W((POS-1)*NDIM+1),
c     +           CC((POS-1)*NDIM+1),NDIM,SELECT,NDIM,AUX,NAUX)
c           CALL ZGEEVY(2,DMAT(1,1),NDIM,W((POS-1)*NDIM+1),
c    +           CC((POS-1)*NDIM+1),NDIM,1,SELECT,NDIM,AUX,NAUX,
c    +           ZTAU,IFAIL)

            IF (INS.GE.1) STOP 'INS.GE.1 in EIGEN0 !!'

c ------------------------------------------------------------------------
            IF (TEST('cc      ')) THEN
              write(6,*) 'eigenvector cc :'
              write(6,FMT=1000) (lm,CC(LM),lm=1,naez*lmmaxd)
            END IF
c ------------------------------------------------------------------------
c
c --->      transformation of eigenvector of GLLKE2 to the 
c           eigenvector of the KKR matrix
c
            IF (NLAYER.GT.1) THEN
              DO 530 N = NLAYER-1,1,-1
                IF (N.GE.POS) GOTO 530
                BASIS=(N-1)*NDIM
c
c --->          EMAT = D(N) = GLLKE2(N)**(-1)
c
                CALL ZCOPY(NDIM*NDIM,GLLKE2(1,1,N),1,DMAT,1)
                CALL ZGETRF(NDIM,NDIM,DMAT,NDIM,IPVT1,INFO)
                CALL CINIT(NDIM*NDIM,EMAT)
                DO LM1=1,NDIM
                  EMAT(LM1,LM1) = CONE
                END DO
                CALL ZGETRS('N',NDIM,NDIM,DMAT,NDIM,
     +               IPVT1,EMAT,NDIM,INFO)
c
c --->          DMAT = D(N)*C(N) (see factorization in routine INVERT)
c
                CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,EMAT(1,1),
     +               NDIM,GLLKE3(1,1,N),NDIM,CZERO,DMAT(1,1),NDIM)
                
                DO 540 LM1=1,NDIM
                  C1 = CZERO
                  DO 550 LM2=1,NDIM  
                    C1 = C1 + DMAT(LM1,LM2)*CC((NLAYER-1)*NDIM+LM2)
 550              END DO
                  CC(BASIS+LM1)= CC(BASIS+LM1) - C1
 540            END DO
          
                IF (N.LT.NLAYER-1) THEN
c     
c --->            DMAT = D(N)*M(N,N+1) (see factorization in routine INVERT)
c
                  CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,EMAT(1,1),
     +                 NDIM,GLLKE1(1,1,N),NDIM,CZERO,DMAT(1,1),NDIM)

                  DO 560 LM1=1,NDIM
                    C1 = CZERO
                    DO 570 LM2=1,NDIM  
                      C1 = C1 + DMAT(LM1,LM2)*CC(BASIS+NDIM+LM2)
 570                END DO
                    CC(BASIS+LM1)= CC(BASIS+LM1) - C1
 560              END DO
                END IF              ! (N.LT.NLAYER-1)
                
 530          END DO                ! N = NLAYER-1,1-1
            END IF                  ! (NLAYER.GT.1)

c ------------------------------------------------------------------------
            IF (TEST('cc      ')) THEN
              write(6,*) 'eigenvector CC :'
              write(6,FMT=1000) (lm,CC(LM),lm=1,naez*lmmaxd)
            END IF
c ------------------------------------------------------------------------

            DO 580 I1 = 1,NAEZ
              INV = KAOEZ(I1)       ! real potential at site I1
              RF  = REFPOT(INV)     ! reference potential at site I1
              DO 590 LM1 = 1,LMMAXD
                ILM = (I1-1)*LMMAXD + LM1
                CL(ILM) = CC(ILM)
     +               *CI**LF(LM1)
     +               *TINVLL(LM1,LM1,I1)
     +               *SQRT(PZSQ(LM1,INV))
     +               /DELTALL(LM1,LM1,RF)
 590          END DO                ! LM1 = 1,LMMAXD
 580        END DO                  ! I1 = 1,NAEZ

            CALL ZNORM(CL,NLM,1.0d0) ! normalization of CL

          END IF                    ! (IO.LT.2)
c ************************************************************************
      END IF                        ! (OPT('full inv'))
c ************************************************************************

      IF (NUM.EQ.0)   IDO = 0
      IF (NUM.EQ.NLM) IUP = 0
      
      LSTART    = .FALSE.
      ELAST     = E
      ISPINLAST = ISPIN
c     ------------------------------------------------------------------------
      if (test('eigen0   '))
     +     write(6,FMT='(f12.6,2I4,1p,2d15.3,L8)')
     +     DREAL(E),NUM,NUMT,wdo,wup,lsym

 999  RETURN

 1000 FORMAT(1p,I6,2d14.3)
 1010 FORMAT(1p,I6,2d14.4,d14.1)
 1020 FORMAT(' W(IMIN) :',1P,D12.4,I6)
 1030 FORMAT(1p,I6,d14.4)
 9030 FORMAT(A1,4x,':',I6,2d12.3,i6,2f10.4)
 9040 FORMAT(' DET :',1p,2d20.10)
      END
