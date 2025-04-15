      SUBROUTINE ESINGPC(ESPC,IEN,NLST,NSPIN,NSTART,NEND,ECORE,LCORE,
     +                   NCORE)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculate the core contribution of the single particle energies
c     l and spin dependent .
c     attention : here are the results of the subroutine corel (stored
c                 in the common block core) used .
c                                        (see notes by b.drittler)
c
c                               b.drittler   may 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER LD
      PARAMETER (LD=LMX-1)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Scalar Arguments ..
      INTEGER NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      REAL*8 ECORE(20,NPOTD),ESPC(0:LD,NATYPD,*),IEN(IEMXD,*)
      INTEGER LCORE(20,NPOTD),NCORE(NPOTD),NLST(*)
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
      DO 40 I1 = NSTART,NEND
        DO 30 IS = 1,NSPIN
c
c---> initialize espc
c
          DO 10 L = 0,LD
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
     +                      (ECORE(N,IPOT)-IEN(NLST(IS),IS))*
     +                      REAL(2*L+1)*REAL(3-NSPIN)
   20     CONTINUE
   30   CONTINUE
   40 CONTINUE
      END
