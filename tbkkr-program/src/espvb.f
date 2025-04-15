c 13.10.95 **************************************************************
      SUBROUTINE ESPVB(DEN,DF,E,ESPV,IE,LMAX,ISPIN,NATYP)
      implicit none
c ************************************************************************
c
c     attention : energy zero ---> electro static zero
c
c                 since input potential and single particle energies
c                 are using muffin tin zero as zero the energy shift
c                 is cancelled in the kinetic energy contribution !
c
c     calculate the valence contribution of the single particle energies
c       l and spin dependent .
c     calculate the change of the valance single particle energies with
c       the local densities of states.
c     recognize that the density of states is always complex (see sub-
c      routine denst) .
c      this includes an implicit energy integration (see deck rholm) .
c     attention : in the case of complex energy integration the single
c                 single particle energies of a representive atom i
c                 is given by :
c
c                                ef
c          esp(i,l,spin) = imag( s dz z  n(z,i,l,spin)  )
c                                eb
c
c          where ef is the fermi energy ,eb the lowest band energy and
c          n(z,i,l,spin) the complex density of states .
c                                        (see notes by b.drittler)
c
c                 modified for bandstructure code
c                               b.drittler   jan 1990
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD,IEMXD,LMAXD
c      PARAMETER (NATYPD=1,IEMXD=128,LMAXD=4)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,E
      INTEGER IE,ISPIN,LMAX,NATYP
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DEN(IEMXD,0:LMAXD,NATYPD,*)
      DOUBLE PRECISION ESPV(0:LMAXD,NATYPD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I1,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG
C     ..

      DO 20 I1 = 1,NATYP
c
c---> implicit energy integration
c
        DO 10 L = 0,LMAX
          ESPV(L,I1,ISPIN) = ESPV(L,I1,ISPIN) +
     +                       DIMAG(E*DEN(IE,L,I1,ISPIN)*DF)
   10   CONTINUE
   20 CONTINUE
      END
