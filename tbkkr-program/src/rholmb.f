c 13.10.95 ***************************************************************
      SUBROUTINE RHOLMB(DEN,DF,GMAT,E,IE,IELAST,IPF,ISPIN,IRWS,KSRA,
     +                  LMAX,NSPIN,NEND,RHO2NS,R2NEF,DRDI,IPAN,IRCUT,
     +                  THETAS,NTCELL,IFUNM,LMSP,SUMNS,VASUM,DENEF,C,PZ,
     +                  FZ,QZ,SZ,CLEB,ICLEB,IEND,JEND,NSHELL)
      implicit none
c ************************************************************************
c     calculate in the paramagnetic case (nspin=1) :
c         the valence charge density times r**2 from the greensfunction
c     calculate in the spin-polarized case (nspin=2) :
c         the valence charge density times r**2 and the valence spin
c         density times r**2 from the greensfunction ,
c         ( convention spin density :=
c                            density(spin up)-density(spin down) )
c     calculate the valence density of states , in the spin-polarized
c      case spin dependent ; splitted into its l-contributions .
c
c     in this subroutine an implicit energy-spin integration is  done :
c        this subroutine is called for each energy and spin value
c        and n(r,e) times df (the energy weight) is calculated .
c
c     recognize that the density of states is always complex also in
c      the case of "real-energy-integation" (ief>0) since in that case
c      the energy integration is done parallel to the real energy axis
c      but not on the real energy axis .
c      in the paramagnetic case only rho2ns(irmd,lmxtsq,natypd,1)
c      is used containing  the charge density times r**2 .
c      in the spin-polarized case rho2ns(...,1) contains the charge
c      density times r**2 and rho2ns(...,2) the spin density times
c      r**2 .
c
c     the charge density is expanded in spherical harmonics :
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
c             (see notes by b.drittler)
c
c     attention : the gaunt coeffients are stored in an index array
c                 (see subroutine gaunt)
c                 the structure part of the greens-function (gmat) is
c                 symmetric in its lm-indices . the gaunt coeffients
c                 are symmetric too (since the are calculated for
c                 real spherical harmonics) . that is why the lm2-
c                 loop only goes up to lm1 and the summands are
c                 multiplied by a factor of 2 in the case of lm1
c                 not equal to lm2 .
c
c     changed for use of shapes and modified for bandstructure code
c
c                                   b.drittler nov. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
c      INTEGER IEMXD
c      PARAMETER (IEMXD=193)
c      INTEGER IRMD,LMAXD,LPOTD
c      PARAMETER (IRMD=1484,LMAXD=4,LPOTD=8)
c      INTEGER NFUND,IRID
c      PARAMETER (NFUND=24,IRID=435)
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,E
      DOUBLE PRECISION C,DENEF
      INTEGER IE,IELAST,IEND,IPF,ISPIN,KSRA,LMAX,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(IEMXD,0:LMAXD,NATYPD,*),
     +               FZ(IRMD,0:LMAXD,*),
     +               GMAT(LMMAXD,LMMAXD,*),
     +               PZ(IRMD,0:LMAXD,*),
     +               QZ(IRMD,0:LMAXD,*),
     +               SZ(IRMD,0:LMAXD,*)
      DOUBLE PRECISION CLEB(*),DRDI(IRMD,*),
     +                 R2NEF(IRMD,LMPOTD,NATYPD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 SUMNS(*),
     +                 THETAS(IRID,NFUND,*),
     +                 VASUM(0:LMAXD,*)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),IPAN(*),IRCUT(0:IPAND,*),
     +        IRWS(*),JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(NATYPD,*),
     +        NSHELL(0:NSHELD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CZERO,EK,EKL,FFZ,GMATL,PPZ,V1
      DOUBLE PRECISION C0LL,FACSYM,PI,SUM
      INTEGER I,I1,ICELL,IFUN,IMT1,IPAN1,IRC1,IS,J,J0,J1,L,L1,L2,LM,LM1,
     +        LM2,LM3,LM3MAX,M
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DENR(IRMD),WR(IRMD,0:LMAXD,0:LMAXD)
      INTEGER IRCUTM(0:IPAND)
      CHARACTER*4 SPN(2),TEXT(7)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CSIMPK,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DIMAG,DMAX1,DREAL,DSQRT,REAL,SQRT,ZSQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','  s=','  p=','  d=','  f=','  g=',
     +     '  h=','  i='/
      DATA CZERO/ (0.0D0,0.0D0)/
c ------------------------------------------------------------------------
      stop ' rholmb'
c
      PI = 4.0D0*DATAN(1.0D0)
      C0LL = 1.0D0/DSQRT(4.0D0*PI)
c
c---> loop over reference atoms
c
      DO 230 I1 = 1,NEND
        ICELL = NTCELL(I1)
        IPAN1 = IPAN(I1)
        IF (KSRA.GE.1) THEN
          EK = ZSQRT(E+E*E/ (C*C))

        ELSE
          EK = ZSQRT(E)
        END IF

        IF (IPAN1.EQ.1) THEN
          IRC1 = IRWS(I1)
          IRCUTM(0) = 0
          IRCUTM(1) = IRC1

        ELSE
          IMT1 = IRCUT(1,I1)
          IRC1 = IRCUT(IPAN(I1),I1)
          DO 10 I = 0,IPAN1
            IRCUTM(I) = IRCUT(I,I1)
   10     CONTINUE
        END IF

        LM3MAX = ICLEB(IEND,3)
        ! Added on 28.3.2000 to switch off the higher L from the charge
        IF (TEST('CHARGE0 '))  LM3MAX = 0 ! Only the spherical charge
c
c---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)
c
        IF (KSRA.GE.1) THEN
          DO 40 L1 = 0,LMAX
            DO 30 L2 = 0,L1
              DO 20 I = 2,IRC1
                WR(I,L1,L2) = PZ(I,L1,I1)*PZ(I,L2,I1) +
     +                        FZ(I,L1,I1)*FZ(I,L2,I1)
   20         CONTINUE
   30       CONTINUE
   40     CONTINUE

        ELSE                        ! (KSRA.GE.1)

          DO 70 L1 = 0,LMAX
            DO 60 L2 = 0,L1
              DO 50 I = 2,IRC1
                WR(I,L1,L2) = PZ(I,L1,I1)*PZ(I,L2,I1)
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE

        END IF                      ! (KSRA.GE.1)
c
c---> first calculate only the spherically symmetric contribution
c
        LM = 0
        DO 140 L = 0,LMAX
          GMATL = CZERO
          EKL = EK*REAL(2*L+1)
          DO 80 M = -L,L
            LM = LM + 1
            GMATL = GMATL + GMAT(LM,LM,I1)
   80     CONTINUE
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
          DENR(1) = CZERO
          IF (KSRA.GE.1) THEN
            IF (IE.NE.IELAST) THEN
              DO 90 I = 2,IRC1
                PPZ = PZ(I,L,I1)
                FFZ = FZ(I,L,I1)
                DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1)) +
     +                    FFZ* (GMATL*FFZ+EKL*SZ(I,L,I1))
                RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                                 C0LL*DIMAG(DF*DENR(I))
   90         CONTINUE

            ELSE
              DO 100 I = 2,IRC1
                PPZ = PZ(I,L,I1)
                FFZ = FZ(I,L,I1)
                DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1)) +
     +                    FFZ* (GMATL*FFZ+EKL*SZ(I,L,I1))
                RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                                 C0LL*DIMAG(DF*DENR(I))
                R2NEF(I,1,I1,ISPIN) = R2NEF(I,1,I1,ISPIN) +
     +                                C0LL*DIMAG(DENR(I))
  100         CONTINUE
            END IF

          ELSE                      ! (KSRA.GE.1)

            IF (IE.NE.IELAST) THEN
              DO 110 I = 2,IRC1
                PPZ = PZ(I,L,I1)
                DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1))
                RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                                 C0LL*DIMAG(DF*DENR(I))
  110         CONTINUE

            ELSE
              DO 120 I = 2,IRC1
                PPZ = PZ(I,L,I1)
                DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1))
                RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                                 C0LL*DIMAG(DF*DENR(I))
                R2NEF(I,1,I1,ISPIN) = R2NEF(I,1,I1,ISPIN) +
     +                                C0LL*DIMAG(DENR(I))
  120         CONTINUE
            END IF

          END IF                    ! (KSRA.GE.1)
c
c
c
          IF (IPAN1.NE.1) THEN
            DO 130 I = IMT1 + 1,IRC1
              DENR(I) = DENR(I)*THETAS(I-IMT1,1,ICELL)*C0LL
  130       CONTINUE
          END IF
c
c--->     calculate density of states
c
          CALL CSIMPK(DENR,DEN(IE,L,I1,ISPIN),IPAN(I1),IRCUTM,
     +                DRDI(1,I1))
          VASUM(L,I1) = VASUM(L,I1) + DIMAG(DF*DEN(IE,L,I1,ISPIN))

  140   CONTINUE                    ! L = 0,LMAX
c
c---> calculate the non spherically symmetric contribution
c        to speed up the pointer jend generated in gaunt is used
c        remember that the wavefunctions are l and not lm dependent
c
        J0 = 1
c
        DO 150 I = 1,IRC1
          DENR(I) = 0.0D0
  150   CONTINUE
        DO 220 LM3 = 2,LM3MAX
          DO 210 L1 = 0,LMAX
            DO 200 L2 = 0,L1
c
              J1 = JEND(LM3,L1,L2)
c
              IF (J1.NE.0) THEN
c
                GMATL = CZERO
c
c---> sum over m1,m2 for fixed lm3,l1,l2
c
                DO 160 J = J0,J1
                  FACSYM = 2.0D0
                  LM1 = ICLEB(J,1)
                  LM2 = ICLEB(J,2)
                  IF (LM1.EQ.LM2) FACSYM = 1.0D0

                  GMATL = GMATL + FACSYM*CLEB(J)*GMAT(LM2,LM1,I1)

c                  IF (TEST('GMAT    ') .AND.
c     +                 ABS((GMAT(LM2,LM1,I1)-GMAT(LM1,LM2,I1))/
c     +                 (GMAT(LM2,LM1,I1)+GMAT(LM1,LM2,I1))) .GT. 1.D-8
c     +                 .AND.
c     +                 ABS(GMAT(LM2,LM1,I1)+GMAT(LM1,LM2,I1))
c     +                 .GT. 1.D-14)
c     +                 write(6,FMT='(3i4,1p,4d15.6)')
c     +                 I1,LM1,LM2,GMAT(LM2,LM1,I1),GMAT(LM1,LM2,I1)

  160           CONTINUE
c
                J0 = J1 + 1
c
                IF (IE.NE.IELAST) THEN
                  DO 170 I = 2,IRC1
                    RHO2NS(I,LM3,I1,ISPIN) = RHO2NS(I,LM3,I1,ISPIN) +
     +                                       DIMAG(DF*GMATL*WR(I,L1,L2))
  170             CONTINUE

                ELSE                ! (IE.NE.IELAST)
                  DO 180 I = 2,IRC1
                    RHO2NS(I,LM3,I1,ISPIN) = RHO2NS(I,LM3,I1,ISPIN) +
     +                                       DIMAG(DF*GMATL*WR(I,L1,L2))
                    R2NEF(I,LM3,I1,ISPIN) = R2NEF(I,LM3,I1,ISPIN) +
     +                                      DIMAG(GMATL*WR(I,L1,L2))
  180             CONTINUE

                END IF              ! (IE.NE.IELAST)

                IF (IPAN1.NE.1 .AND. LMSP(ICELL,LM3).GT.0) THEN
                  IFUN = IFUNM(ICELL,LM3)
                  DO 190 I = IMT1 + 1,IRC1
                    DENR(I) = DENR(I) + GMATL*WR(I,L1,L2)*
     +                        THETAS(I-IMT1,IFUN,ICELL)
  190             CONTINUE

                END IF

              END IF                ! (J1.NE.0)

  200       CONTINUE

  210     CONTINUE

  220   CONTINUE

        IF (IPAN1.NE.1) THEN
c
c---> calculate non-sph. contribution - not added to density of state !
c
          CALL CSIMPK(DENR,V1,IPAN(I1),IRCUTM,DRDI(1,I1))
c
          SUMNS(I1) = SUMNS(I1) + DIMAG(DF*V1)
        END IF

  230 CONTINUE                      ! I1 = 1,NEND
c
      DENEF = 1.D-14

      IF (IE.EQ.IELAST .AND. IELAST.GT.1) THEN
        DO 260 IS = 1,NSPIN
          DO 250 I1 = 1,NEND
            DO 240 L = 0,LMAX
              DENEF = DENEF - 2.0D0*DIMAG(DEN(IELAST,L,I1,IS))/PI*
     +                REAL(NSHELL(I1))/REAL(NSPIN)
  240       CONTINUE
  250     CONTINUE
  260   CONTINUE

c        DENEF = DMAX1(DENEF,1.0D0)

        IF (NSPIN.NE.2) THEN
          IF (IPAN1.GT.1) THEN
            WRITE (IPF,FMT=9040)
          ELSE
            WRITE (IPF,FMT=9010)
          END IF
        ELSE IF (IPAN1.GT.1) THEN
          WRITE (IPF,FMT=9030) SPN(ISPIN)
        ELSE
          WRITE (IPF,FMT=9000) SPN(ISPIN)
        END IF

        DO 280 I1 = 1,NEND
          SUM = 0.0D0
          DO 270 L = 0,LMAX
            SUM = SUM + VASUM(L,I1)
  270     CONTINUE
          WRITE (IPF,FMT=9020) I1,SUM, (TEXT(L+1),VASUM(L,I1),L=0,LMAX)
          stop 'rholmb'
          IF (IPAN1.GT.1) WRITE (IPF,FMT=9050) SUM + SUMNS(I1),
     +        SUMNS(I1)
  280   CONTINUE
      END IF                        ! (IE.EQ.IELAST .AND. IELAST.GT.1)

      RETURN

 9000 FORMAT (2x,/,2x,'ws-sphere valence-charges for spin ',a4)
 9010 FORMAT (2x,/,2x,'ws-sphere valence-charges')
 9020 FORMAT (1x,I3,' su1=',f10.6,5 (a4,f10.6))
 9030 FORMAT (2x,/,2x,'ws-cell valence-charges for spin ',a4)
 9040 FORMAT (2x,/,2x,'ws-cell valence-charges')
 9050 FORMAT (2x,2x,' sum=',f10.6,' non spherical contribution=',f10.6)

      END
