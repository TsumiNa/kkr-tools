c 13.10.95 ***************************************************************
      SUBROUTINE ESPCB(ESPC,NSPIN,NATYP,ECORE,LCORE,NCORE)
      implicit none
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c                 since input potential and single particle energies
c                 are using muffin tin zero as zero the energy shift
c                 is cancelled in the kinetic energy contribution !
c
c     calculate the core contribution of the single particle energies
c     l and spin dependent .
c     attention : here are the results of the subroutine corel (stored
c                 in the common block core) used .
c                                        (see notes by b.drittler)
c
c                 modified for bandstructure code
c                               b.drittler   jan 1990
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
c      INTEGER LMAXD
c      PARAMETER (LMAXD=4)
C     ..
C     .. Scalar Arguments ..
      INTEGER NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION ECORE(20,*),ESPC(0:LMAXD,NATYPD,*)
      INTEGER LCORE(20,*),NCORE(*)
C     ..
C     .. Local Scalars ..
      INTEGER I1,IPOT,IS,L,N
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
c
c---> loop over reference atoms
c
      DO 40 I1 = 1,NATYP
        DO 30 IS = 1,NSPIN
c
c---> initialize espc
c
          DO 10 L = 0,LMAXD
            ESPC(L,I1,IS) = 0.0D0
   10     CONTINUE
c
c---> determine correct potential indices
c
          IPOT = NSPIN* (I1-1) + IS
c
c---> loop over all core states
c
          DO 20 N = 1,NCORE(IPOT)
            L = LCORE(N,IPOT)
            ESPC(L,I1,IS) = ESPC(L,I1,IS) +
     +                      ECORE(N,IPOT)*REAL(2*L+1)*REAL(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      END
