
C ************************************************************************
      SUBROUTINE TMATRX(ALPHA,DET,AR,CR,DRDI,E,ICST,INS,LMAX,PZ,QZ,FZ,
     +                  SZ,PNS,QNS,TMATLL,R,VINS,VM2Z,Z,IRWS,IPAN,IRCUT,
     +                  IRMIN,KSRA,C,DROR,RS,S,CLEB,LOFLM,ICLEB,IEND)
      implicit none
C ************************************************************************
c
c     this is the driver for the wavefunctions and t - matrices
c
c     in case of non spherical input potential the non spher.
c     wavefunctions are approximated inside a given sphere
c     with the nearly r - independent matrices :
c
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
c             (see notes by b.drittler)
c
c             changed for band structure code
c
c                               b.drittler   nov. 1989
c                               valerio 14.6.99
C ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
      DOUBLE COMPLEX CONE
      PARAMETER (CONE= (1.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DET,E
      DOUBLE PRECISION C,Z
      INTEGER ICST,IEND,INS,IPAN,IRMIN,IRWS,KSRA,LMAX
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),CR(LMMAXD,LMMAXD),
     +               FZ(IRMD,0:LMAXD),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               PZ(IRMD,0:LMAXD),QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +               QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
     +               TMATLL(LMMAXD,LMMAXD)
      DOUBLE PRECISION CLEB(*),DRDI(IRMD),DROR(IRMD),R(IRMD),
     +                 RS(IRMD,0:LMAXD),S(0:LMAXD),
     +                 VINS(IRMIND:IRMD,LMPOTD),VM2Z(IRMD)
      INTEGER ICLEB(NCLEB,4),IRCUT(0:IPAND),LOFLM(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX EK
      INTEGER I,INFO,IRC1,LM1,LMMAX,NSRA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX CMAT(LMMAXD,LMMAXD,IRMIND:IRMD),
     +               DMAT(LMMAXD,LMMAXD,IRMIND:IRMD),DR(LMMAXD,LMMAXD),
     +               EFAC(LMMAXD),PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               PZLM(LMMAXD,IRMIND:IRMD,2),
     +               QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +               QZLM(LMMAXD,IRMIND:IRMD,2),TMAT(0:LMAXD)
      DOUBLE PRECISION VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IPVT(LMMAXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,CPLXWB,IRWNS,REGNS,VLLNS,WFTSCA,ZGETRF
C     ..
      LMMAX = (LMAX+1)* (LMAX+1)

      IF (KSRA.GT.1) THEN
         NSRA = 2

      ELSE

         NSRA = 1
      END IF

c     this was giving a value of 3 for KSRA=2 which was WRONG!
c      NSRA = KSRA + 1

      CALL CINIT(LMMAXD*LMMAXD,TMATLL)
c
c---> determine wavefunctions and t matrix for spherical averaged pot.
c
      CALL CPLXWB(PZ,FZ,QZ,SZ,TMAT,ALPHA,E,EK,KSRA,IPAN,IRCUT,IRWS,VM2Z,
     +            DRDI,R,Z,LMAX,C,DROR,RS,S)
      DET = CONE
c
      IF (INS.EQ.0) THEN
        DO 10 LM1 = 1,LMMAX
          TMATLL(LM1,LM1) = TMAT(LOFLM(LM1))
c          write(113,*) LM1,TMATLL(LM1,LM1)
   10   CONTINUE
        
      ELSE
c
c---> non spherical input potential
c
        IRC1 = IRCUT(IPAN)
c
c---> determine the lm,lm' dependent potential
c
        CALL VLLNS(IRMIN,IRC1,LMMAX,VNSPLL,VINS,CLEB,ICLEB,IEND)
c
c---> get wfts of same magnitude by scaling with efac
c
        CALL WFTSCA(DRDI,EFAC,LMAX,PZ,QZ,IRMIN,IRWS,IPAN,IRCUT,FZ,SZ,
     +              NSRA,PZLM,QZLM,PZEKDR,QZEKDR,EK,LOFLM)
c
c---> determine the irregular non sph. wavefunction
c
        CALL IRWNS(CR,DR,EFAC,QNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +             NSRA,PZLM,QZLM,PZEKDR,QZEKDR,QNS(1,1,IRMIND,1),CMAT,
     +             QNS(1,1,IRMIND,2),DMAT)
c
c---> determine the regular non sph. wavefunction
c
        CALL REGNS(AR,TMATLL,EFAC,PNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +             PZLM,QZLM,PZEKDR,QZEKDR,EK,PNS(1,1,IRMIND,1),CMAT,
     +             PNS(1,1,IRMIND,2),DMAT,NSRA)
c
        DO 20 LM1 = 1,LMMAX
          TMATLL(LM1,LM1) = TMATLL(LM1,LM1) + TMAT(LOFLM(LM1))
   20   CONTINUE
c
c---> calculate determinant of alpha matrix
c
        CALL ZGETRF(LMMAXD,LMMAXD,AR,LMMAXD,IPVT,INFO)
        DO 30 I = 1,LMMAXD
          IF (IPVT(I).NE.I) DET = -DET
          DET = AR(I,I)*DET
   30   CONTINUE


      END IF

      END
