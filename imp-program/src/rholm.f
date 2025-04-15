      SUBROUTINE RHOLM(DEN,DF,GMAT,E,IE,IELAST,IPF,ISPIN,IRWS,KVREL,
     +                 LMAX,LMAXSQ,NSPIN,NATREF,NSTART,NEND,RHO2NS,DRDI,
     +                 IPAN,IRCUT,THETAS,NTCELL,IFUNM,LMSP,SUMNS,VASUM,
     +                 C,PZ,FZ,QZ,SZ,CLEB,ICLEB,IEND,JEND,
     +                 KSYMMAT,KESYM,DTB1)
      IMPLICIT NONE
c-----------------------------------------------------------------------
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
c                 symmetric in its lm-indices , therefore only one
c                 half of the matrix is calculated in the subroutine
c                 for the back-symmetrisation . the gaunt coeffients
c                 are symmetric too (since the are calculated for
c                 real spherical harmonics) . that is why the lm2-
c                 loop only goes up to lm1 and the summands are
c                 multiplied by a factor of 2 in the case of lm1
c                 not equal to lm2 .
c
c                               b.drittler   may 1987
c                                   changed  dec 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER IRMD,LMAXD,LMX,LPOTD
      PARAMETER (irmd=1484,lmaxd=4,LMX=LMAXD+1,lpotd=8)
      INTEGER NFUND,IRID
      PARAMETER (NFUND=289,irid=435)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
CASYM
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,*)
      INTEGER KSYMMAT,KESYM
      REAL*8 GL
CASYM
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E
      REAL*8 C
      INTEGER IE,IELAST,IEND,IPF,ISPIN,KVREL,LMAX,LMAXSQ,NATREF,NEND,
     +        NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DEN(IEMXD,0:LMAXD,NATYPD,*),
     +               FZ(IRMD,0:LMAXD,NATYPD),GMAT(LMMAXD,LMMAXD,*),
     +               PZ(IRMD,0:LMAXD,NATYPD),QZ(IRMD,0:LMAXD,NATYPD),
     +               SZ(IRMD,0:LMAXD,NATYPD)
      REAL*8 CLEB(*),DRDI(IRMD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),SUMNS(NATYPD),
     +                 THETAS(IRID,NFUND,*),VASUM(0:LMAXD,NATYPD)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),IPAN(*),IRCUT(0:IPAND,*),
     +        IRWS(*),JEND(LMPOTD,0:LMAXD,0:LMAXD),LMSP(NATYPD,*),
     +        NTCELL(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO,EK,EKL,FFZ,GMATL,PPZ,V1
      REAL*8 C0LL,FACSYM,PI,SUM
      INTEGER I,I1,ICELL,IFUN,IMT1,IPAN1,IRC1,J,J0,J1,KSRA,L,L1,L2,LM3,
     +        LM3MAX,LN,LN1,LN2,LNE,LNS,N
C     ..
C     .. Local Arrays ..
      COMPLEX*16 DENR(IRMD),WR(IRMD,0:LMAXD,0:LMAXD)
      INTEGER IRCUTM(0:IPAND)
      CHARACTER*4 SPN(2),TEXT(7)
C     ..
C     .. External Subroutines ..
      EXTERNAL CSIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDSQRT,DATAN,DIMAG,DSQRT,REAL
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA SPN,TEXT/'down','up  ','  s=','  p=','  d=','  f=','  g=',
     +     '  h=','  i='/
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
      PI = 4.0D0*DATAN(1.0D0)
      C0LL = 1.0D0/DSQRT(4.0D0*PI)
c

      LM3MAX = ICLEB(IEND,3)
      N = 0
c
c---> loop over reference atoms
c
      DO 200 I1 = NSTART,NEND
        ICELL = NTCELL(I1)
        IPAN1 = IPAN(I1)
        IF ((KVREL.EQ.1.AND.I1.GT.NATREF) .OR. KVREL.EQ.2) THEN
          KSRA = 1
          EK = CDSQRT(E+E*E/ (C*C))

        ELSE

          KSRA = 0
          EK = CDSQRT(E)
        END IF

        IF (IPAN1.EQ.1) THEN
          IRC1 = IRWS(I1)
          IRCUTM(0) = 0
          IRCUTM(1) = IRC1

        ELSE
          IMT1 = IRCUT(1,I1)
c         write(6,*) 'i1,imt1,icell.irc1' , i1,imt1,icell,irc1
          IRC1 = IRCUT(IPAN(I1),I1)
          DO 10 I = 0,IPAN1
            IRCUTM(I) = IRCUT(I,I1)
   10     CONTINUE
        END IF

        LN = LMAXSQ*N
c
c---> set up of wr(ir,l1,l2) = pz(ir,l1)*pz(ir,l2)
c
        IF (KSRA.EQ.1) THEN
          DO 40 L1 = 0,LMAX
            DO 30 L2 = 0,L1
              DO 20 I = 2,IRC1
                WR(I,L1,L2) = PZ(I,L1,I1)*PZ(I,L2,I1) +
     +                        FZ(I,L1,I1)*FZ(I,L2,I1)
   20         CONTINUE
   30       CONTINUE
   40     CONTINUE

        ELSE

          DO 70 L1 = 0,LMAX
            DO 60 L2 = 0,L1
              DO 50 I = 2,IRC1
                WR(I,L1,L2) = PZ(I,L1,I1)*PZ(I,L2,I1)
   50         CONTINUE
   60       CONTINUE
   70     CONTINUE

        END IF
c
c---> first calculate only the spherically symmetric contribution
c
CASYM-------------------------
        DO 120 L = 0,LMAX
          GL=0.0D0
          GMATL = CZERO
          LNS = L*L + 1
          LNE = LNS + 2*L
          EKL = EK*REAL(2*L+1)
          DO 80 LN1 = LNS,LNE
            GMATL = GMATL + GMAT(LN1,LN1,N+1)
CASYM------>
            GL=GL+REAL(DTB1(LN1,LN1,N+1))
   80     CONTINUE
        IF (I1 .GT. NATREF) THEN
        IF (KSYMMAT .GT. 0) THEN
        IF (IE  .GE. KESYM) THEN
CTESTASYM           write (*,*) 'calculation in rhoin:',' l:',l,' gl:',gl
           EKL=EK*GL
        END IF
        END IF
        END IF
CASYM-------------------------
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
          DENR(1) = CZERO
          IF (KSRA.EQ.1) THEN
            DO 90 I = 2,IRC1
              PPZ = PZ(I,L,I1)
              FFZ = FZ(I,L,I1)
              DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1)) +
     +                  FFZ* (GMATL*FFZ+EKL*SZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*DENR(I))
   90       CONTINUE

          ELSE

            DO 100 I = 2,IRC1
              PPZ = PZ(I,L,I1)
              DENR(I) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*DENR(I))
  100       CONTINUE
          END IF
c
c
c
          IF (IPAN1.NE.1) THEN
            DO 110 I = IMT1 + 1,IRC1
              DENR(I) = DENR(I)*THETAS(I-IMT1,1,ICELL)*C0LL
  110       CONTINUE
          END IF
c
c---> calculate density of states
c
          CALL CSIMPK(DENR,DEN(IE,L,I1,ISPIN),IPAN(I1),IRCUTM,
     +                DRDI(1,I1))
          VASUM(L,I1) = VASUM(L,I1) + DIMAG(DF*DEN(IE,L,I1,ISPIN))
  120   CONTINUE
c
c---> calculate the non spherically symmetric contribution
c        to speed up the pointer jend generated in gaunt is used
c        remember that the wavefunctions are l and not lm dependent
c
        J0 = 1
c
        DO 130 I = 1,IRC1
          DENR(I) = 0.0D0
  130   CONTINUE
        DO 190 LM3 = 2,LM3MAX
          DO 180 L1 = 0,LMAX
            DO 170 L2 = 0,L1
c
              J1 = JEND(LM3,L1,L2)
c
              IF (J1.NE.0) THEN
c
                GMATL = CZERO
c
c---> sum over m1,m2 for fixed lm3,l1,l2
c
                DO 140 J = J0,J1
                  FACSYM = 2.0D0
                  LN1 = ICLEB(J,1)
                  LN2 = ICLEB(J,2)
                  IF (LN1.EQ.LN2) FACSYM = 1.0D0
                  GMATL = GMATL + FACSYM*CLEB(J)*DF*GMAT(LN2,LN1,N+1)
  140           CONTINUE
c
                J0 = J1 + 1
c
                DO 150 I = 2,IRC1
                  RHO2NS(I,LM3,I1,ISPIN) = RHO2NS(I,LM3,I1,ISPIN) +
     +                                     DIMAG(GMATL*WR(I,L1,L2))
  150           CONTINUE

                IF (IPAN1.NE.1 .AND. LMSP(ICELL,LM3).GT.0) THEN
                  IFUN = IFUNM(ICELL,LM3)
C ???? V1 = GMATL  or next statement ????
                  V1 = GMATL/DF
                  DO 160 I = IMT1 + 1,IRC1
                    DENR(I) = DENR(I) + V1*WR(I,L1,L2)*
     +                        THETAS(I-IMT1,IFUN,ICELL)
  160             CONTINUE

                END IF

              END IF

  170       CONTINUE

  180     CONTINUE

  190   CONTINUE

        IF (IPAN1.NE.1) THEN
c
c---> calculate non-sph. contribution - not added to density of state !
c
          CALL CSIMPK(DENR,V1,IPAN(I1),IRCUTM,DRDI(1,I1))
c
          SUMNS(I1) = SUMNS(I1) + DIMAG(DF*V1)
        END IF

        N = N + 1
  200 CONTINUE
c
c
      IF (IE.EQ.IELAST) THEN
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

        DO 220 I1 = NSTART,NEND
          SUM = 0.0D0
          DO 210 L = 0,LMAX
            SUM = SUM + VASUM(L,I1)
  210     CONTINUE
          WRITE (IPF,FMT=9020) SUM, (TEXT(L+1),VASUM(L,I1),L=0,LMAX)
          IF (IPAN1.GT.1) WRITE (IPF,FMT=9050) SUM + SUMNS(I1),
     +        SUMNS(I1)
  220   CONTINUE
      END IF

 9000 FORMAT (2x,/,2x,'ws-cell valence-charges for spin ',a4)
 9010 FORMAT (2x,/,2x,'ws-cell valence-charges')
 9020 FORMAT (2x,'sum=',f10.6,5 (a4,f10.6))
 9030 FORMAT (2x,/,2x,'ws-cell valence-charges for spin ',a4)
 9040 FORMAT (2x,/,2x,'ws-cell valence-charges')
 9050 FORMAT (2x,'sum=',f10.6,' non spherical contribution=',f10.6)
      END
