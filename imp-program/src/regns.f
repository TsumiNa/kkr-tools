
      SUBROUTINE REGNS(AR,BR,EFAC,PNS,VNSPLL,ICST,IRMIN,IRWS,IPAN,IRCUT,
     +                 PZLM,QZLM,PZEKDR,QZEKDR,EK,ADER,AMAT,BDER,BMAT,
     +                 NSRA,LKONV,irmd,irmind)
      Implicit None
c-----------------------------------------------------------------------
c     determines the regular non spherical wavefunctions , the
c       alpha matrix and the t - matrix in the n-th. born appro-
c       ximation ( n given by input parameter icst )
c
c
c     using the wave functions pz and qz ( regular and irregular
c       solution ) of the spherically averaged potential , the
c       regular wavefunction pns is determined by
c
c           pns(ir,lm1,lm2) = ar(ir,lm1,lm2)*pz(ir,l1)
c
c                                   + br(ir,lm1,lm2)*qz(ir,l1)
c
c      the matrices ar and br are determined by integral equations
c        containing pns and only the non spherical contributions of
c        the potential , stored in vinspll . these integral equations
c        are  solved iteratively with born approximation up to given n.
c
c     the original way of writing the cr and dr matrices in the equa-
c        tions above caused numerical troubles . therefore here are used
c        rescaled ar and br matrices :
c
c              ~
c              ar(ir,lm1,lm2) = sqrt(e)**(l1-l2)
c
c                             * ar(ir,lm1,lm2)*((2*l2-1)!!/(2*l1-1)!!)
c
c              ~
c              br(ir,lm1,lm2) = sqrt(e)**(-l1-l2)
c
c                             * br(ir,lm1,lm2)/((2*l1-1)!!*(2*l2-1)!!)
c
c
c
c     for lloyd's formular is only the determinant of the alpha -
c        matrix is needed which is identical with the determinant
c        of the rescaled ar - matrix at the innerst point .
c
c     the non spherical t - matrix is the br matrix at r(irc)
c
c     numerical tests showed that it is sufficient to integrate only
c        ca. 70 points inwards to get reliable results . then the
c        rescaled ar and br matrices are nearly r independent . the
c        first point of the inwards integration is irmin
c
c     total energies are converged for first born approximation
c
c     modified for the use of shape functions
c
c                              (see notes by b.drittler)
c
c                                b.drittler   mar.  1989
c-----------------------------------------------------------------------
c     modified by R. Zeller      Aug. 1994
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER LMAXD
      PARAMETER (lmaxd=4)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      COMPLEX*16 CONE
      PARAMETER (CONE= (1.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 EK
      INTEGER ICST,IPAN,IRMIN,IRWS,NSRA,LKONV,irmd,irmind
C     ..
C     .. Array Arguments ..
      COMPLEX*16 ADER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +     AMAT(LMMAXD,LMMAXD,IRMIND:IRMD),AR(LMMAXD,LMMAXD),
     +     BDER(LMMAXD,LMMAXD,IRMIND:IRMD),
     +     BMAT(LMMAXD,LMMAXD,IRMIND:IRMD),BR(LMMAXD,LMMAXD),
     +     EFAC(*),PNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     +     PZEKDR(LMMAXD,IRMIND:IRMD,2),
     +     PZLM(LMMAXD,IRMIND:IRMD,2),
     +     QZEKDR(LMMAXD,IRMIND:IRMD,2),
     +     QZLM(LMMAXD,IRMIND:IRMD,2)
      REAL*8 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 EFAC1,EFAC2
      INTEGER I,IPNS,IR,IRC1,J,LM1,LM2
C     ..
C     .. External Subroutines ..
      EXTERNAL CSINWD,CSOUT,WFINT,WFINT0
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0
C     ..
      IRC1 = IRCUT(IPAN)
      DO 70 I = 0,ICST
c---> set up integrands for i-th born approximation
        IF (I.EQ.0) THEN
          CALL WFINT0(ADER,BDER,PZLM,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     +                NSRA,LKONV,irmd,irmind)

        ELSE
          CALL WFINT(PNS,ADER,BDER,QZEKDR,PZEKDR,VNSPLL,IRMIN,IRC1,
     $                NSRA,LKONV,irmd,irmind)
        END IF
c---> call integration subroutines
        CALL CSINWD(ADER,AMAT,IRMIN,IPAN,IRCUT)
        CALL CSOUT(BDER,BMAT,IRMIN,IPAN,IRCUT)
        DO 20 IR = IRMIN,IRC1
          DO 10 LM2 = 1,LMMAXD
            AMAT(LM2,LM2,IR) = CONE + AMAT(LM2,LM2,IR)
   10     CONTINUE
   20   CONTINUE
c---> calculate non sph. wft. in i-th born approximation
        DO 60 J = 1,NSRA
          DO 50 IR = IRMIN,IRC1
            DO 40 LM1 = 1,LMMAXD
              DO 30 LM2 = 1,LMMAXD
                PNS(LM1,LM2,IR,J) = (AMAT(LM1,LM2,IR)*PZLM(LM1,IR,J)+
     +                              BMAT(LM1,LM2,IR)*QZLM(LM1,IR,J))
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
   70 CONTINUE
      DO 90 LM2 = 1,LMMAXD
        EFAC2 = EFAC(LM2)
c---> store alpha and t - matrix
        DO 80 LM1 = 1,LMMAXD
          EFAC1 = EFAC(LM1)
          AR(LM1,LM2) = AMAT(LM1,LM2,IRMIN)
c---> t-matrix
          BR(LM1,LM2) = BMAT(LM1,LM2,IRC1)*EFAC1*EFAC2/EK
   80   CONTINUE
   90 CONTINUE
c---> in case of muffin tin calculation : fill pns
      DO 130 J = 1,NSRA
        DO 120 IR = IRC1 + 1,IRWS
          DO 110 LM1 = 1,LMMAXD
            DO 100 LM2 = 1,LMMAXD
              PNS(LM1,LM2,IR,J) = AMAT(LM1,LM2,IRC1)*PZLM(LM1,IR,J) +
     +                            BMAT(LM1,LM2,IRC1)*QZLM(LM1,IR,J)
  100       CONTINUE
  110     CONTINUE
  120   CONTINUE
  130 CONTINUE
c---> rescale with efac
      IPNS = MAX0(IRC1,IRWS)
      DO 170 J = 1,NSRA
        DO 160 LM2 = 1,LMMAXD
          EFAC2 = EFAC(LM2)
          DO 150 IR = IRMIN,IPNS
            DO 140 LM1 = 1,LMMAXD
              PNS(LM1,LM2,IR,J) = PNS(LM1,LM2,IR,J)*EFAC2
  140       CONTINUE
  150     CONTINUE
  160   CONTINUE
  170 CONTINUE
      END
