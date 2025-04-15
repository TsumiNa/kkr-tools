      SUBROUTINE RHONS(DEN,DF,DRDI,GMAT,E,IE,IELAST,IOWFCT,ISPIN,IRMIN,
     +                 IRWS,LMAX,NATREF,NSPIN,NSTART,NEND,RHO2NS,IPAN,
     +                 IRCUT,THETAS,NTCELL,IFUNM,LMSP,KSRA,QNS,PNS,AR,
     +                 CR,C,PZ,FZ,QZ,SZ,CLEB,ICLEB,JEND,IEND,
     +                 KESYM,DTB1,KSYMMAT)
      IMPLICIT NONE
C____________ASYM___CHANGED
c-----------------------------------------------------------------------
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
c     attention : in this subroutine the arrays pns and qns ( con-
c                 taining the non spherical wavefunctions ) are linear
c                 arrays to read faster from buffer memory !
c                 this feature is only used if WFDISK=.TRUE.
c
c             (see notes by b.drittler)
c
c                               b.drittler   aug. 1988
c
c     modified for the use of shape functions
c
c     attention : irmin + 3 has to be less then imt
c                 if shape functions are used
c
c                               b.drittler   july 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NSRAD
      PARAMETER (NATYPD=38,NSRAD=2)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER IRMD,IRNSD,LMAXD,LMX,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,LMX=LMAXD+1,lpotd=8)
      INTEGER NFUND,IRID
      PARAMETER (NFUND=289,irid=435)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD,LMMAXD
      PARAMETER (LMPOTD= (LPOTD+1)**2,LMMAXD=LMX**2)
      INTEGER IRLLD
      PARAMETER (IRLLD= (IRNSD+1)*LMMAXD*LMMAXD*NSRAD)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E
      REAL*8 C
      INTEGER IE,IELAST,IEND,IOWFCT,ISPIN,KSRA,LMAX,NATREF,NEND,NSPIN,
     +        NSTART
C     ..
C     .. Array Arguments ..
      COMPLEX*16 AR(LMMAXD,LMMAXD,*),CR(LMMAXD,LMMAXD,*),
     +               DEN(IEMXD,0:LMAXD,NATYPD,*),
     +               FZ(IRMD,0:LMAXD,NATYPD),GMAT(LMMAXD,LMMAXD,*),
     +               PNS(IRLLD,NATYPD),PZ(IRMD,0:LMAXD,NATYPD),
     +               QNS(IRLLD,NATYPD),
     +               QZ(IRMD,0:LMAXD,NATYPD),SZ(IRMD,0:LMAXD,NATYPD)
c     use next line, if WFDISK is true
c     +               PNS(IRLLD,1),PZ(IRMD,0:LMAXD,NATYPD),QNS(IRLLD,1),
      REAL*8 CLEB(*),DRDI(IRMD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),THETAS(IRID,NFUND,*)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),IPAN(*),IRCUT(0:IPAND,*),
     +        IRMIN(*),IRWS(*),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 DENNS,EK,V1
      REAL*8 SUM
      INTEGER I,I1,I1DISK,ICELL,IMT1,IPAN1,IRC1,IRMIN1,L,LM,LMMAX,M,NSRA
      LOGICAL WFDISK
C     ..
C     .. Local Arrays ..
      COMPLEX*16 CDEN(IRMD,0:LMAXD),CDENNS(IRMD),EFAC(LMMAXD)
      REAL*8 ZSUM(0:LMAXD,NATYPD),ZSUMNS(NATYPD)
      INTEGER IRCUTM(0:IPAND)
      CHARACTER*4 SPN(2),TEXT(7)
C     ..
CASYM
      INTEGER KESYM,KSYMMAT
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,*)
CASYM end

C     .. External Subroutines ..
      EXTERNAL CSIMPK,RHOIN,RHOOUT,WFREAD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDSQRT,DIMAG,REAL
C     ..
C     .. Save statement ..
      SAVE ZSUM,ZSUMNS,SPN,TEXT
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','  s=','  p=','  d=','  f=','  g=',
     +     '  h=','  i='/
C     ..

c     wfdisk=true : write WFs to the disk (this should be changed
c     in the routines fpimpu, tmatns, and rhons)
c      WFDISK = .true.
      WFDISK = .false.

      LMMAX = (LMAX+1)* (LMAX+1)
c
      IF ((KSRA.EQ.1.AND.NSTART.GT.NATREF) .OR. KSRA.GT.1) THEN
        NSRA = 2
        EK = CDSQRT(E+E*E/ (C*C))

      ELSE

        NSRA = 1
        EK = CDSQRT(E)
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
      REWIND IOWFCT
c
c---> loop over representive atoms
c
      DO 50 I1 = NSTART,NEND
        I1DISK = I1
        IF (WFDISK) I1DISK = 1
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
c       write(6,*) 'rhoout ',icell,imt1
c
c---> read wavefunctions and matrices from the buffer memory
c
        IF (WFDISK) CALL WFREAD(PNS,QNS,IOWFCT)

c
c---> calculate charge density from irmin up to irc1 with the full
c        non-spherical wavefunctions
c
CASYM------->changed !
        CALL RHOOUT(CDEN,DF,GMAT,EK,I1,IE,ISPIN,IRMIN1,IRC1,LMAX,LMMAX,
     +              NSPIN,NATREF,PNS(1,I1DISK),QNS(1,I1DISK),RHO2NS,
     +              THETAS,ICELL,IFUNM,IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,
     +              ICLEB,IEND,KESYM,DTB1,KSYMMAT)
c
c---> calculate charge density from the origin up to irmin with
c        nearly r - independent matrices
c
        CALL RHOIN(AR(1,1,I1),CDEN,CR(1,1,I1),DF,GMAT,EK,I1,IE,IELAST,
     +             ISPIN,IRMIN1,IRC1,LMAX,LMMAX,NSPIN,NATREF,RHO2NS,
     +             NSRA,EFAC,PZ,FZ,QZ,SZ,CLEB,ICLEB,JEND,IEND,
     +             KESYM,DTB1,KSYMMAT)
CASYM------->end changed!
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
   40   CONTINUE

        IF (IPAN1.GT.1) THEN
          IF (IE.EQ.1) ZSUMNS(I1) = 0.0D0
          CALL CSIMPK(CDENNS,DENNS,IPAN1,IRCUTM,DRDI(1,I1))
          DEN(IE,0,I1,ISPIN) = DEN(IE,0,I1,ISPIN) + DENNS
          ZSUMNS(I1) = ZSUMNS(I1) + DIMAG(DF*DENNS)
        END IF

   50 CONTINUE

c
      IF (IE.EQ.IELAST) THEN
        IF (NSPIN.EQ.2) THEN
          WRITE (6,FMT=9000) SPN(ISPIN)

        ELSE

          WRITE (6,FMT=9010)
        END IF

        DO 70 I1 = NSTART,NEND
          SUM = 0.0D0
          DO 60 L = 0,LMAX
            SUM = SUM + ZSUM(L,I1)
   60     CONTINUE
          WRITE (6,FMT=9020) SUM, (TEXT(L+1),ZSUM(L,I1),L=0,LMAX)
          IF (IPAN1.GT.1) WRITE (6,FMT=9030) SUM + ZSUMNS(I1),
     +        ZSUMNS(I1)
   70   CONTINUE
      END IF
c
      RETURN
c
   80 CONTINUE
      WRITE (6,FMT=9040)
      STOP


c
c

 9000 FORMAT (2x,/,2x,'ws-cell valence-charges for spin ',a4)
 9010 FORMAT (2x,/,2x,'ws-cell valence-charges')
 9020 FORMAT (2x,'sum=',f10.6,5 (a4,f14.10))
 9030 FORMAT (2x,'sum=',f10.6,2x,'sumns=',f14.10)
c 9020 FORMAT (2x,'sum=',f10.6,5 (a4,f10.6))
c 9030 FORMAT (2x,'sum=',f10.6,2x,'sumns=',f10.6)
 9040 FORMAT (13x,'error reading of the wcfts on the buffer memory !')
      END
