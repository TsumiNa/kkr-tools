c 18.10.95 ***************************************************************
      SUBROUTINE RHOINB(AR,CDEN,CR,DF,GMAT,EK,I1,IE,IELAST,ISPIN,IRMIN,
     +                  IRC1,LMAX,LMMAX,NSPIN,NATREF,RHO2NS,R2NEF,NSRA,
     +                  EFAC,PZ,FZ,QZ,SZ,CLEB,ICLEB,JEND,IEND)
      implicit none
c ************************************************************************
c
c     calculates the charge density inside r(irmin) in case
c      of a non spherical input potential .
c
c     fills the array cden for the complex density of states
c
c      the non spher. wavefunctions are approximated in that region
c       in the following way :
c
c           the regular one (ir < irmin = irws-irns) :
c
c              pns(ir,lm1,lm2) = pz(ir,l1) * ar(lm1,lm2)
c
c          where pz is the regular wavefct of the spherically symmetric
c          part of the potential and ar the alpha matrix .
c          (see subroutine regns)
c
c
c           the irregular one (ir < irmin) :
c
c              qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)
c
c                                    + qz(ir,l1) * dr(lm1,lm2)
c
c          where pz is the regular and qz is the irregular
c          wavefct of the spherically symmetric part of the
c          potential and cr , dr the matrices calculated
c          at the point irmin .  (see subroutine irwns)
c
c     to safe cpu time first all the matrices which are not r -
c      dependent are multiplied and summed up - taking into account
c      that the muffin tin wavefunctions pz and pz are only l dependent
c      in the final step the charge density is obtained by multiply-
c      ing this with the muffin tin wavefunctions .
c
c     attention : therefore the gaunt coeffients which are used here
c                 are ordered in a special way !   (see subroutine
c                 gaunt)
c
c                 remember that the matrices ar,cr,dr are rescaled !
c                 (see subroutines irwns and regns)
c
c                 arrays rho2ns and cden are initialize in subroutine
c                 rhoout .
c
c
c     the structured part of the greens-function (gmat) is symmetric in
c       its lm-indices , therefore only one half of the matrix is
c       calculated in the subroutine for the back-symmetrisation .
c       the gaunt coeffients are symmetric too (since the are calculated
c       using the real spherical harmonics) . that is why the lm2- and
c       the lm02- loops are only only going up to lm1 or lm01 and the
c       summands are multiplied by a factor of 2 in the case of lm1 .ne.
c       lm2 or lm01 .ne. lm02 .
c
c             (see notes by b.drittler)
c
c                               b.drittler   aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
c      INTEGER IRMD,LMAXD,LPOTD
c      PARAMETER (IRMD=1484,LMAXD=4,LPOTD=8)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
c      INTEGER NCLEB
c      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER I1,IE,IELAST,IEND,IRC1,IRMIN,ISPIN,LMAX,LMMAX,NATREF,
     +        NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD),CDEN(IRMD,0:LMAXD),
     +               CR(LMMAXD,LMMAXD),EFAC(*),FZ(IRMD,0:LMAXD,*),
     +               GMAT(LMMAXD,LMMAXD),PZ(IRMD,0:LMAXD,*),
     +               QZ(IRMD,0:LMAXD,*),SZ(IRMD,0:LMAXD,*)
      DOUBLE PRECISION CLEB(*),R2NEF(IRMD,LMPOTD,NATYPD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*)
      INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CZERO,EFAC1,EFAC2,EKL,FFZ,GMATL,PPZ,V1,V2
      DOUBLE PRECISION C0LL
      INTEGER I,IR,J,J0,J1,L,L1,L2,LM1,LM2,LM3,LM3MAX,LN,LN2,LN3,M,N
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX VR(LMMAXD,LMMAXD),WF(IRMD,0:LMAXD,0:LMAXD),
     +               WR(LMMAXD,LMMAXD)
C     ..
C     .. External Functions ..
      DOUBLE COMPLEX ZDOTU
      EXTERNAL ZDOTU
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DIMAG,DSQRT,REAL
C     ..
C     .. Save statement ..
      SAVE CZERO
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
C     C0LL = 1/sqrt(4*pi)
      C0LL = 1.0d0/DSQRT(16.0D0*DATAN(1.0D0))
c
      IF (I1.GT.NATREF) THEN
        N = I1 - (NATREF+1)

      ELSE

        N = I1 - 1
      END IF

      IF (NATREF.EQ.-1) N = 0

      LN = LMMAX*N
      LM3MAX = ICLEB(IEND,3)
c
c---> set up array wr(lm1,lm2)
c        use first vr
c
      DO 20 LM2 = 1,LMMAX
        LN2 = LN + LM2
        V2 = EFAC(LM2)*EFAC(LM2)*GMAT(LN2,LN2)
        DO 10 LM1 = 1,LMMAXD
          VR(LM1,LM2) = EK*CR(LM1,LM2) + V2*AR(LM1,LM2)
   10   CONTINUE
   20 CONTINUE
c
c---> using symmetry of structural green's function
c
      DO 50 LM2 = 2,LMMAX
        LN2 = LN + LM2
        EFAC2 = 2.0D0*EFAC(LM2)
        DO 40 LM3 = 1,LM2 - 1
          LN3 = LN + LM3
          V1 = EFAC2*GMAT(LN3,LN2)*EFAC(LM3)
          DO 30 LM1 = 1,LMMAXD
            VR(LM1,LM2) = VR(LM1,LM2) + V1*AR(LM1,LM3)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
c
      DO 70 LM1 = 1,LMMAX
        EFAC1 = EFAC(LM1)
        WR(LM1,LM1) = ZDOTU(LMMAX,AR(LM1,1),LMMAXD,VR(LM1,1),LMMAXD)/
     +                (EFAC1*EFAC1)
        DO 60 LM2 = 1,LM1 - 1
c
c---> using symmetry of gaunt coeffients
c
          EFAC2 = EFAC(LM2)
          WR(LM1,LM2) = (ZDOTU(LMMAX,AR(LM1,1),LMMAXD,VR(LM2,1),LMMAXD)+
     +                  ZDOTU(LMMAX,AR(LM2,1),LMMAXD,VR(LM1,1),LMMAXD))/
     +                  (EFAC1*EFAC2)
   60   CONTINUE
   70 CONTINUE
c
c---> set up array wf(l1,l2) = pz(l1)*pz(l2)
c
      IF (NSRA.EQ.2) THEN
        DO 100 L1 = 0,LMAX
          DO 90 L2 = 0,L1
            DO 80 IR = 2,IRMIN
              WF(IR,L1,L2) = PZ(IR,L1,I1)*PZ(IR,L2,I1) +
     +                       FZ(IR,L1,I1)*FZ(IR,L2,I1)
   80       CONTINUE
   90     CONTINUE
  100   CONTINUE

      ELSE
        DO 130 L1 = 0,LMAX
          DO 120 L2 = 0,L1
            DO 110 IR = 2,IRMIN
              WF(IR,L1,L2) = PZ(IR,L1,I1)*PZ(IR,L2,I1)
  110       CONTINUE
  120     CONTINUE
  130   CONTINUE
      END IF
c
c---> first calculate only the spherically symmetric contribution
c     remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
      DO 190 L = 0,LMAX
        GMATL = CZERO
        EKL = EK*REAL(2*L+1)
        DO 140 M = -L,L
          LM1 = L* (L+1) + M + 1
          GMATL = GMATL + WR(LM1,LM1)
  140   CONTINUE
c
        IF (IE.EQ.IELAST) THEN
          IF (NSRA.EQ.2) THEN
            DO 150 I = 2,IRMIN
              PPZ = PZ(I,L,I1)
              FFZ = FZ(I,L,I1)
              CDEN(I,L) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1)) +
     +                    FFZ* (GMATL*FFZ+EKL*SZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*CDEN(I,L))
              R2NEF(I,1,I1,ISPIN) = R2NEF(I,1,I1,ISPIN) +
     +                              C0LL*DIMAG(CDEN(I,L))
  150       CONTINUE

          ELSE
            DO 160 I = 2,IRMIN
              PPZ = PZ(I,L,I1)
              CDEN(I,L) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*CDEN(I,L))
              R2NEF(I,1,I1,ISPIN) = R2NEF(I,1,I1,ISPIN) +
     +                              C0LL*DIMAG(CDEN(I,L))
  160       CONTINUE
          END IF

        ELSE
          IF (NSRA.EQ.2) THEN
            DO 170 I = 2,IRMIN
              PPZ = PZ(I,L,I1)
              FFZ = FZ(I,L,I1)
              CDEN(I,L) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1)) +
     +                    FFZ* (GMATL*FFZ+EKL*SZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*CDEN(I,L))
  170       CONTINUE

          ELSE
            DO 180 I = 2,IRMIN
              PPZ = PZ(I,L,I1)
              CDEN(I,L) = PPZ* (GMATL*PPZ+EKL*QZ(I,L,I1))
              RHO2NS(I,1,I1,ISPIN) = RHO2NS(I,1,I1,ISPIN) +
     +                               C0LL*DIMAG(DF*CDEN(I,L))
  180       CONTINUE
          END IF

        END IF

  190 CONTINUE
c
c---> calculate the non spherically symmetric contribution
c        to speed up the pointer jend generated in gaunt is used
c        remember that the wavefunctions are l and not lm dependent
c
      J0 = 1
c
      DO 250 LM3 = 2,LM3MAX
        DO 240 L1 = 0,LMAX
          DO 230 L2 = 0,L1
c
            J1 = JEND(LM3,L1,L2)
c
            IF (J1.NE.0) THEN
c
              GMATL = CZERO
c
c---> sum over m1,m2 for fixed lm3,l1,l2
c
              DO 200 J = J0,J1
                LM1 = ICLEB(J,1)
                LM2 = ICLEB(J,2)
                GMATL = GMATL + CLEB(J)*WR(LM1,LM2)
  200         CONTINUE

c
              J0 = J1 + 1
c
              IF (IE.EQ.IELAST) THEN
                DO 210 I = 2,IRMIN
                  RHO2NS(I,LM3,I1,ISPIN) = RHO2NS(I,LM3,I1,ISPIN) +
     +                                     DIMAG(DF*GMATL*WF(I,L1,L2))
                  R2NEF(I,LM3,I1,ISPIN) = R2NEF(I,LM3,I1,ISPIN) +
     +                                    DIMAG(GMATL*WF(I,L1,L2))
  210           CONTINUE

              ELSE
                GMATL = DF*GMATL
                DO 220 I = 2,IRMIN
                  RHO2NS(I,LM3,I1,ISPIN) = RHO2NS(I,LM3,I1,ISPIN) +
     +                                     DIMAG(GMATL*WF(I,L1,L2))
  220           CONTINUE
              END IF

            END IF

  230     CONTINUE

  240   CONTINUE

  250 CONTINUE

      END
