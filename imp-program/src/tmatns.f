      SUBROUTINE TMATNS(AR,CR,DRDI,E,ICST,IOWFCT,LMAX,PZ,QZ,FZ,SZ,PNS,
     +     QNS,TMATLL,VINS,vm2z,vspsme,IRWS,IPAN,IRCUT,IRMIN,NSRA,C,
     +     CLEB,ICLEB,IEND,LOFLM,TMAT,NONSPH,LKONV)
      Implicit None
c-----------------------------------------------------------------------
c
c     this is the driver for the non spherically sym. wavefunctions
c        and t - matrices .
c        test calculations showed that it is suffient to integrate
c        only 30 - 50 points inwards for the irregular solution ,
c        getting with that the proper starting matrices (especially
c        a good alpha matrix) for the outwards integration back to
c        the wigner seitz or muffin tin sphere .
c        inside this 30 - 50 points - given in the array irns - the
c        non spher. wavefunctions are approximated with the nearly
c        r - independent matrices in the following way :
c
c           the regular one (ir < irmin = irc - irns) :
c
c              pns(ir,lm1,lm2) = pz(ir,l1) * alfmat(lm1,lm2)
c
c          where pz is the regular wavefct of the spherically symmetric
c          part of the potential and alfmat the alpha matrix .
c
c
c           the irregular one (ir < irmin = irc - irns) :
c
c              qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)
c
c                                    + qz(ir,l1) * dr(lm1,lm2)
c
c          where pz is the regular and qz is the irregular
c          wavefct of the spherically symmetric part of the
c          potential and cr , dr the matrices calculated
c          at the point irmin = irc - irns .
c
c     to save storage the non spherical wavefunctions are stored only
c        from the sphere boundary up to irc - irns
c
c     after that the delta t - matrix is calculated and symmetrized
c        (transformed into the irreduzible representation)
c
c     attention : in this subroutine the arrays pns and qns ( con-
c                 taining the non spherical wavefunctions ) are linear
c                 arrays to write faster on buffer memory !
c                 this feature is only used if WFDISK=.TRUE.
c
c
c             (see notes by b.drittler)
c
c                               b.drittler   aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER IRMD,IRNSD,LMAXD,LPOTD,LMX
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER IRID
      PARAMETER (irid=435)
      INTEGER IPAND
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      PARAMETER (IPAND=80)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 E
      REAL*8 C
      INTEGER ICST,IEND,IOWFCT,IPAN,IRMIN,IRWS,LMAX,NSRA,LKONV
      LOGICAL NONSPH
C     ..
C     .. Array Arguments ..
      COMPLEX*16 AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +     FZ(IRMD,0:LMAXD),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +     PZ(IRMD,0:LMAXD),QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +     QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),TMAT(0:LMAXD),
     +     TMATLL(LMMAXD,LMMAXD)
      REAL*8 CLEB(NCLEB,2),DRDI(IRMD),VINS(IRMIND:IRMD,LMPOTD),
     $     vm2z(irmd),vspsme(irmd)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),LOFLM(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO,EK
      INTEGER IRC1,LM1,LMMAX
      INTEGER LMMKONV
      LOGICAL WFDISK
C     ..
C     .. Local Arrays ..
      COMPLEX*16 CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),
     +     DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD),
     +     EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +     PZLM(LMMAXD,IRMIND:IRMD,2),
     +     QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +     QZLM(LMMAXD,IRMIND:IRMD,2)
      REAL*8 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,IRWNS,REGNS,VLLNS,WFTSCA,WFWRIT
C     ..
C     .. Save statement ..
      SAVE CZERO
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDSQRT
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.D0,0.D0)/
C     ..
      LMMAX = (LMAX+1)* (LMAX+1)
      IF (LKONV.NE.LMAXD) THEN 
         LMMKONV=(LKONV+1)*(LKONV+1)
      ELSE
         LMMKONV=LMMAX
      END IF
      CALL CINIT(LMMAXD*LMMAXD,TMATLL)
      IF (NONSPH) THEN
        DO 10 LM1 = 1,LMMAX
          TMATLL(LM1,LM1) = TMAT(LOFLM(LM1))
   10   CONTINUE
      ELSE

c     wfdisk=true : write WFs to the disk (this should be changed
c     in the routines fpimpu, tmatns, and rhons)
c        WFDISK = .true.
        WFDISK = .false.

        IF (NSRA.GT.1) THEN
          EK = CDSQRT(E+E*E/ (C*C))

        ELSE

          EK = CDSQRT(E)
        END IF
c
c
c---> loop over representive atoms
c
c
        IRC1 = IRCUT(IPAN)
c
c---> determine the lm,lm' dependent potential
c
c
c     In a new version you should calculate the difference 
c     between actual l=0 pot and smeared l=0 pot and put
c     the difference in the V_LL's
        CALL VLLNS(IRMIN,IRC1,LMMAX,VNSPLL,VINS,vm2z,vspsme,
     $       CLEB,ICLEB,IEND,irmd,irmind)
c
c---> get wfts of same magnitude by scaling with efac
c
        CALL WFTSCA(DRDI,EFAC,LMAX,PZ,QZ,IRMIN,IRWS,IPAN,IRCUT,FZ,SZ,
     +       NSRA,PZLM,QZLM,PZEKDR,QZEKDR,EK,LOFLM,irmd,irmind)
c
c---> determine the irregular non sph. wavefunction
c
        CALL IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +       NSRA,PZLM,QZLM,PZEKDR,QZEKDR,QNS(1,1,IRMIND,1),CMAT,
     +       QNS(1,1,IRMIND,2),DMAT,LKONV,irmd,irmind)
c
c---> determine the regular non sph. wavefunction
c
        CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +       PZLM,QZLM,PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +       PNS(1,1,IRMIND,2),DMAT,NSRA,LKONV,irmd,irmind)
c

        DO 20 LM1 = 1,LMMKONV
          TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1))
   20   CONTINUE
c
c---> store regular wavefunctions and matrices on the buffer memory
c
        IF (WFDISK) CALL WFWRIT(PNS,QNS,IOWFCT)
      END IF
      
      RETURN

c
c---> stop here in case of i/o error
c
   30 CONTINUE
      WRITE (6,FMT=9000)
      STOP


 9000 FORMAT (13x,'error writing of the wcfts on the buffer memory !')
      END
