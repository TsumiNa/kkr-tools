C ************************************************************************
      SUBROUTINE VLLNS(IRMIN,IRC,LMMAX,VNSPLL,VINS,CLEB,ICLEB,IEND)
      implicit none
C ************************************************************************
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
C     ..
C     .. Scalar Arguments ..
      INTEGER IEND,IRC,IRMIN,LMMAX
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(*),VINS(IRMIND:IRMD,LMPOTD),
     +                 VNSPLL(LMMAXD,LMMAXD,IRMIND:IRMD)
      INTEGER ICLEB(NCLEB,4)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,J,LM1,LM2,LM3
C     ..
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
          VNSPLL(LM1,LM2,I) = VNSPLL(LM1,LM2,I) + CLEB(J)*VINS(I,LM3)
c     write(6,*) lm1,lm2,vnspll(lm1,lm2,i)
   40   CONTINUE
c     stop 'vllns'
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

      END
