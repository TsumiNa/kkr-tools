      SUBROUTINE VLLNS(IRMIN,IRC,LMMAX,VNSPLL,VINS,vm2z,vspsme,
     $     CLEB,ICLEB,IEND,irmd,irmind)
      Implicit None
c-----------------------------------------------------------------------
c     to determine the non - spherical wavefunctions the potential
c         has to be lm1 and lm2 dependent . the potential is stored
c         only as lm dependent , therefore a transformation in the
c         following way has to be done :
c
c
c        vinspll(r,lm1,lm2)   =   {  c(lm1,lm2,lm3) *vinss(r,lm3)  }
c
c                                  (summed over lm3 at the right site )
c
c        where c(lm1,lm2,lm3) are the gaunt coeffients .
c
c
c             (see notes by b.drittler)
c
c     attention : the gaunt coeffients are stored in an index array
c                  only for lm1.gt.lm2
c                 (see subroutine gaunt)
c
c     attention : here only the non spherical contributions of the
c                 input potential are needed . this mean vins(ir,lm=1)
c                 is not taken into account
c
c                               b.drittler   july 1988
c
c     modified by T. Korhonen    Jun. 1997
c       Now adds the spherical component of the potential
c       difference to VLLNS: vm2z-vspsme, where vm2z is the
c       total spherical potential and vspsme is the spherical
c       potential used in the radial Schr. equations
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER LMAXD,LPOTD
      PARAMETER (lmaxd=4,lpotd=8)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NCLEB
      PARAMETER (NCLEB=LMPOTD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND,IRC,IRMIN,LMMAX,irmd,irmind
C     ..
C     .. Array Arguments ..
      REAL*8 CLEB(NCLEB,2),VINS(IRMIND:IRMD,LMPOTD),
     +     VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD),vm2z(irmd),
     $     vspsme(irmd)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,J,LM1,LM2,LM3
      Double Precision fpi,rfpi
C     ..

      fpi  = 4.d0*4.d0*Datan(1.d0)
      rfpi = Dsqrt(fpi)

      DO 30 LM1 = 1,LMMAX
        DO 20 LM2 = 1,LM1
          DO 10 IR = IRMIN,IRC
            VNSPLL(LM1,LM2,IR) = 0.0D0
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
c
      DO 50 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        DO 40 I = IRMIN,IRC
          VNSPLL(LM1,LM2,I) = VNSPLL(LM1,LM2,I) +CLEB(J,1)*VINS(I,LM3)
   40   CONTINUE
   50 CONTINUE
c
c---> use symmetry of the gaunt coef.
c
      DO 80 LM1 = 1,LMMAX
        DO 70 LM2 = 1,LM1 - 1
          DO 60 I = IRMIN,IRC
            VNSPLL(LM2,LM1,I) = VNSPLL(LM1,LM2,I)
   60     CONTINUE
   70   CONTINUE
   80 CONTINUE

      Do lm1 = 1, lmmax
        Do i = irmin, irc
c     Remember, that vdelta is already defined to be the
c     difference between potentials, see subroutine intpot
ccccc          vnspll(lm1,lm1,i) = vnspll(lm1,lm1,i) + vdelta(i)
          vnspll(lm1,lm1,i) = vnspll(lm1,lm1,i) +
     $         (vm2z(i)-vspsme(i))
        End Do
      End Do


      END
