c ************************************************************************
      INTEGER FUNCTION NPOLES(E,ISPIN)
      implicit none
c ************************************************************************
c     p.zahn, april 96
c ------------------------------------------------------------------------
C     .. arguments
      include 'inc.fi'
      DOUBLE PRECISION E
      INTEGER ISPIN
c     .. local
      INTEGER I
c     .. common block
      INTEGER NDIFF(NSPIND,20),NZERO(NSPIND)
      DOUBLE PRECISION EZERO(NSPIND,20),EBOT
      COMMON /EZERO / EBOT,EZERO,NZERO,NDIFF

c ------------------------------------------------------------------------
      NPOLES = 0
      RETURN
      END
