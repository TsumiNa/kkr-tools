

c 29.04.96 ***************************************************************
      SUBROUTINE PZNORM(PZSQ,E,IPF,ISPIN,IRWS,KSRA,
     +                  LMAX,NSPIN,DRDI,IPAN,IRCUT,THETAS,
     +                  NTCELL,C,PZ,FZ)
      implicit none
c ************************************************************************
c     calculate the norm of the regular solution of the radial
c     Schroedinger equation 
c        in nonrelativistic case PZ
c        in scalar relat. case   PZ,FZ
c
c     p. zahn, april 1996
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX   E
      DOUBLE PRECISION C
      INTEGER          IPF,IPAN,IRWS,ISPIN,KSRA,LMAX,NSPIN,NTCELL
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX   FZ(IRMD,0:LMAXD),
     +                 PZSQ(*),
     +                 PZ(IRMD,0:LMAXD)
      DOUBLE PRECISION DRDI(*),
     +                 THETAS(IRID,NFUND,*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   CZERO,W
      DOUBLE PRECISION C0LL,PI
      INTEGER          I,I1,ICELL,IMT1,IRC1,L,L1,LM
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX   WR(IRMD,0:LMAXD)
      INTEGER          IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CSIMPK,TEST
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DSQRT,ZABS
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
      PI = 4.0D0*DATAN(1.0D0)
      C0LL = 1.0D0/DSQRT(4.0D0*PI)

      IF (IPAN.EQ.1) THEN
        IRC1 = IRWS
        IRCUTM(0) = 0
        IRCUTM(1) = IRC1
      ELSE
        IMT1 = IRCUT(1)
        IRC1 = IRCUT(IPAN)
        DO 10 I = 0,IPAN
          IRCUTM(I) = IRCUT(I)
 10     CONTINUE
      END IF

c
c---> set up of wr(ir,l1) = pz(ir,l1)*pz(ir,l1)
c
      IF (KSRA.GE.1) THEN
        
        DO 40 L1 = 0,LMAX
          WR(1,L1) = CZERO
          DO 20 I = 2,IRC1
            WR(I,L1) = PZ(I,L1)*PZ(I,L1) +
     +                 FZ(I,L1)*FZ(I,L1)
 20       CONTINUE
 40     CONTINUE
        
      ELSE                          ! (KSRA.GE.1)
        
        DO 70 L1 = 0,LMAX
          WR(1,L1) = CZERO
          DO 50 I = 2,IRC1
            WR(I,L1) = PZ(I,L1)*PZ(I,L1)
c            WR(I,L1) = ZABS(WR(I,L1))
 50       CONTINUE
 70     CONTINUE
        
      END IF                        ! (KSRA.GE.1)
c     
      LM = 0
      DO 140 L = 0,LMAX
c     
        IF (IPAN.NE.1) THEN
          DO 130 I = IMT1 + 1,IRC1
            WR(I,L) = WR(I,L)*THETAS(I-IMT1,1,NTCELL)*C0LL
 130      CONTINUE
        END IF
c     
c--- >  calculate length of PZ*PZ + FZ*FZ
c     
        CALL CSIMPK(WR(1,L),W,IPAN,IRCUTM,DRDI(1))

        DO 150 L1 = -L,L
          LM = LM + 1
          PZSQ(LM) = W
          IF (TEST('testpzsq')) 
     +         write(6,FMT='(2i5,f12.4)') 
     +         l,l1,datan(dimag(w)/dreal(w))
 150    END DO
 140  CONTINUE                      ! L = 0,LMAX

      RETURN

      END
