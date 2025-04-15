c 13.10.95 ***************************************************************
      SUBROUTINE RHONSB(DEN,DF,DRDI,GMAT,E,IE,IELAST,ISPIN,IRMIN,IRWS,
     +                  LMAX,NSPIN,NSTART,NEND,RHO2NS,R2NEF,IPAN,IRCUT,
     +                  THETAS,NTCELL,IFUNM,LMSP,KSRA,AR,CR,PNS,QNS,C,
     +                  PZ,FZ,QZ,SZ,CLEB,ICLEB,IEND,JEND,DENEF)
      implicit none
c ************************************************************************
c
c     this is the driver for the valence charge density and density
c        of states in the case of a non spherically symmetric input
c        potential - calling subroutine rhoin and rhooub
c
c      the full non spherical wavefcts are calculated only for the
c       first irns1 points inwards the wigner seitz sphere (see e.g.
c       subroutine irwns) ; the charge density for that region is
c       calculated by subroutine rhooub .
c
c      inside these irns1 points  the non spher. wavefunctions are
c       approximated by r - independent matrices times the muffin
c       tin wave functions ;the charge density for that region is
c       calculated by the subroutine rhoin .
c
c     general information :
c
c     calculates in the paramagnetic case (nspin=1) :
c         the valence charge density times r**2 from the greensfunction
c     calculates in the spin-polarized case (nspin=2) :
c         the valence charge density times r**2 and the valence spin
c         density times r**2 from the greensfunction ,
c         ( convention spin density :=
c                            density(spin up)-density(spin down) )
c
c     therefore an implicit energy-spin integration is done :
c        this subroutine is called for each energy and spin value
c        and n(r,e) times df (the energy weight) is calculated .
c      in the paramagnetic case this is added only to rho2ns(...,1)
c         which contains the charge density .
c      in the spin-polarized case this is added to rho2ns(...,1),
c         which contains the charge density , and spin dependent
c         added to or subtracted from rho2ns(...,2) , which con-
c         tains in that case the spin density .
c        when the loop is finished rho2ns(...,1)  contains the
c        charge density times r**2 and in the spin-polarized case
c        rho2ns(...,2) contains the spin density times r**2 .
c
c     the charge density is developed in spherical harmonics :
c
c             rho(r) =   { rho(lm,r) * y(r,lm) }       (summed over lm)
c
c          rho(lm,r) =   { do rho(r) * y(r,lm)         (integrated over
c                                                           unit sphere)
c     in the case of spin-polarization :
c       the spin density is developed in spherical harmonics :
c
c            sden(r) =   { sden(lm,r) * y(r,lm) }      (summed over lm)
c
c         sden(lm,r) =   { do sden(r) * y(r,lm)        (integrated over
c                                                           unit sphere)
c     n(r,e) is developed in
c
c        n(r,e) = { y(r,l'm') * n(l'm',lm,r,e) * y(r,lm) }
c
c     therefore a faltung of n(l'm',lm,r,e) with the gaunt coeffients
c     has to be used to calculate the lm-contribution of the charge
c     density .
c
c
c     calculate the valence density of states , in the spin-polarized
c      case spin dependent .
c     recognize that the density of states is always complex also in
c      the case of "real-energy-integation" (ief>0) since in that case
c      the energy integration is done parallel to the real energy axis
c      but not on the real energy axis .
c     in the last energy-spin loop the l-contribution of the valence
c      charge is calculated .
c
c     modified for the use of shape functions
c
c     attention : irmin + 3 has to be less then imt
c                 if shape functions are used
c
c     changed for band structure code
c
c                               b.drittler   oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
c      INTEGER IEMXD
c      PARAMETER (IEMXD=193)
c      INTEGER IRMD,IRNSD,LMAXD,LPOTD
c      PARAMETER (IRMD=1484,IRNSD=508,LMAXD=4,LPOTD=8)
c      INTEGER NFUND,IRID
c      PARAMETER (NFUND=24,IRID=435)
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
      INTEGER LMPOTD,LMMAXD
      PARAMETER (LMPOTD= (LPOTD+1)**2,LMMAXD= (LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,E
      DOUBLE PRECISION C,DENEF
      INTEGER IE,IELAST,IEND,ISPIN,KSRA,LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Save statement ..
      SAVE ZSUM,ZSUMNS,SPN,TEXT
 
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD,*),
     +               CR(LMMAXD,LMMAXD,*),
     +               DEN(IEMXD,0:LMAXD,NATYPD,*),
     +               FZ(IRMD,0:LMAXD,*),
     +               GMAT(LMMAXD,LMMAXD,*),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,*),
     +               PZ(IRMD,0:LMAXD,*),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,*),
     +               QZ(IRMD,0:LMAXD,*),
     +               SZ(IRMD,0:LMAXD,*)
      DOUBLE PRECISION CLEB(*),DRDI(IRMD,*),
     +                 R2NEF(IRMD,LMPOTD,NATYPD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),IPAN(*),IRCUT(0:IPAND,*),
     +        IRMIN(*),IRWS(*),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DENNS,EK,V1
      DOUBLE PRECISION PI,SUM
      INTEGER I,I1,ICELL,IMT1,IPAN1,IRC1,IRMIN1,IS,L,LM,LMMAX,M,NSRA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD)
      DOUBLE PRECISION ZSUM(0:LMAXD,NATYPD),ZSUMNS(NATYPD)
      INTEGER IRCUTM(0:IPAND)
      CHARACTER*4 SPN(2),TEXT(7)
      LOGICAL TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL CSIMPK,RHOINB,RHOOUTB,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DIMAG,DMAX1,REAL,SQRT
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','  s=','  p=','  d=','  f=','  g=',
     +     '  h=','  i='/
C     ..
      PI = 4.0D0*DATAN(1.0D0)

      LMMAX = (LMAX+1)* (LMAX+1)
c
      IF (KSRA.GE.1) THEN
        NSRA = 2
        EK = SQRT(E+E*E/ (C*C))

      ELSE

        NSRA = 1
        EK = SQRT(E)
      END IF
c
c---> set up efac(lm) = sqrt(e))**l/(2l - 1)!!
c
      EFAC(1) = 1.0D0
      V1 = 1.0D0
      DO 20 L = 1,LMAX
        V1 = V1*EK/REAL(2*L-1)
        DO 10 M = -L,L
          LM = L* (L+1) + M + 1
          EFAC(LM) = V1
   10   CONTINUE
   20 CONTINUE
c     
c
c---> loop over representive atoms
c
      DO 50 I1 = NSTART,NEND
        IRMIN1 = IRMIN(I1)
        IPAN1 = IPAN(I1)
        IF (IPAN1.EQ.1) THEN
          IRC1 = IRWS(I1)
          IRCUTM(0) = 0
          IRCUTM(1) = IRC1

        ELSE
          ICELL = NTCELL(I1)
          IMT1 = IRCUT(1,I1)
          IRC1 = IRCUT(IPAN(I1),I1)
          DO 30 I = 0,IPAN1
            IRCUTM(I) = IRCUT(I,I1)
   30     CONTINUE
        END IF
c
c---> calculate charge density from irmin up to irc1 with the full
c        non-spherical wavefunctions
c
        CALL RHOOUTB(CDEN,DF,GMAT(1,1,I1),EK,I1,IE,ISPIN,IRMIN1,IRC1,
     +               LMAX,LMMAX,PNS(1,1,IRMIND,1,I1),
     +               QNS(1,1,IRMIND,1,I1),IELAST,RHO2NS,R2NEF,THETAS,
     +               ICELL,IFUNM,IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,
     +               IEND)
c
c---> calculate charge density from the origin up to irmin with
c        nearly r - independent matrices
c
        CALL RHOINB(AR(1,1,I1),CDEN,CR(1,1,I1),DF,GMAT(1,1,I1),EK,I1,IE,
     +              IELAST,ISPIN,IRMIN1,IRC1,LMAX,LMMAX,NSPIN,-1,RHO2NS,
     +              R2NEF,NSRA,EFAC,PZ,FZ,QZ,SZ,CLEB,ICLEB,JEND,IEND)
      
c
c---> calculate complex density of states
c
        DO 40 L = 0,LMAX
          IF (IE.EQ.1) ZSUM(L,I1) = 0.0D0
c
c---> call integration subroutine
c
          CALL CSIMPK(CDEN(1,L),DEN(IE,L,I1,ISPIN),IPAN1,IRCUTM,
     +                DRDI(1,I1))
          ZSUM(L,I1) = ZSUM(L,I1) + DIMAG(DF*DEN(IE,L,I1,ISPIN))
c       write(6,*) l,i1,zsum(l,i1)
   40   CONTINUE
c       stop ' rhonsb'

        IF (IPAN1.GT.1) THEN
          IF (IE.EQ.1) ZSUMNS(I1) = 0.0D0
          CALL CSIMPK(CDENNS,DENNS,IPAN1,IRCUTM,DRDI(1,I1))
          DEN(IE,0,I1,ISPIN) = DEN(IE,0,I1,ISPIN) + DENNS
          ZSUMNS(I1) = ZSUMNS(I1) + DIMAG(DF*DENNS)
        END IF

   50 CONTINUE


c
      DENEF = 0.0D0
      IF (IE.EQ.IELAST .AND. IELAST.GT.1) THEN
        DO 80 IS = 1,NSPIN
          DO 70 I1 = 1,NEND
            DO 60 L = 0,LMAX
              DENEF = DENEF - 2.0D0*DIMAG(DEN(IELAST,L,I1,IS))/PI/
     +                REAL(NSPIN)       !*NEND) corrected 14.10.99
   60       CONTINUE
   70     CONTINUE
   80   CONTINUE
c        DENEF = DMAX1(DENEF,1.0D0)
        IF (NSPIN.EQ.2) THEN
          WRITE (6,FMT=9000) SPN(ISPIN)

        ELSE

          WRITE (6,FMT=9010)
        END IF

        DO 100 I1 = NSTART,NEND
          SUM = 0.0D0
          DO 90 L = 0,LMAX
            SUM = SUM + ZSUM(L,I1)
   90     CONTINUE
          WRITE (6,FMT=9020) SUM, (TEXT(L+1),ZSUM(L,I1),L=0,LMAX)
          IF (IPAN1.GT.1) WRITE (6,FMT=9030) SUM + ZSUMNS(I1),
     +        ZSUMNS(I1)
  100   CONTINUE
      END IF
       
      RETURN

 9000 FORMAT (2x,/,2x,'ws-cell valence-charges for spin ',a4)
 9010 FORMAT (2x,/,2x,'ws-cell valence-charges')
 9020 FORMAT (2x,'sum=',f11.7,5 (a4,f11.7))
 9030 FORMAT (2x,'sum=',f11.7,2x,'sumns=',f11.7)

      END
