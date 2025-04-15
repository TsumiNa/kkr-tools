      SUBROUTINE ESINGPV(DEN,DF,E,ESPV,IE,IELAST,LMAX,ISPIN,IREF,NSHELL,
     +                   NSTART,NEND,NATREF,DSELOC)
      Implicit None
c-----------------------------------------------------------------------
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
c          esp(i,l,spin) = imag( { dz  (z - ef) n(z,i,l,spin)  )
c                                eb
c
c          where ef is the fermi energy ,eb the lowest band energy and
c          n(z,i,l,spin) the complex density of states .
c                                        (see notes by b.drittler)
c
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
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E
      REAL*8 DSELOC
      INTEGER IE,IELAST,ISPIN,LMAX,NATREF,NEND,NSTART
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DEN(IEMXD,0:LD,NATYPD,*)
      REAL*8 ESPV(0:LD,NATYPD,*)
      INTEGER IREF(*),NSHELL(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO
      INTEGER I1,L,NREF
C     ..
C     .. Local Arrays ..
      COMPLEX*16 N(0:LD,NATYPD,NSPIND)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG,REAL
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..
c
c--->initialize dseloc
c
      IF (IE.EQ.1 .AND. ISPIN.EQ.1) DSELOC = 0.0D0
c
c---> loop over reference atoms
c
      DO 30 I1 = NSTART,NEND
        DO 10 L = 0,LMAX
c
c---> initalize espv and n for the first energy-spin loop
c
          IF (IE.EQ.1) THEN
            ESPV(L,I1,ISPIN) = 0.0D0
            N(L,I1,ISPIN) = CZERO
          END IF
c
c---> implicit energy integration : n contains in the last energy-spin
c     loop the integrated density of states , espv the integral of the
c     energy times the density of states .
c
          N(L,I1,ISPIN) = N(L,I1,ISPIN) + DEN(IE,L,I1,ISPIN)*DF
          ESPV(L,I1,ISPIN) = ESPV(L,I1,ISPIN) +
     +                       DIMAG(E*DEN(IE,L,I1,ISPIN)*DF)
          IF (IE.EQ.IELAST) THEN
c
c---> gather the parts in the last energy loop
c     BE CAREFULL!!! THIS IS FOR THE TEMPERATURE!!!
            ESPV(L,I1,ISPIN) = ESPV(L,I1,ISPIN) -
     +                         REAL(E)*DIMAG(N(L,I1,ISPIN))
c           write(6,'(13x,2e12.6)') e
c           write(6,'(13x,i3,1x,e12.6)') l,espv(l,i1,ispin)
          END IF
 
   10   CONTINUE
c
c---> calculate change of valance single particle energies
c
        IF (IE.EQ.IELAST .AND. I1.GT.NATREF) THEN
          NREF = IREF(I1-NATREF)
          DO 20 L = 0,LMAX
            DSELOC = DSELOC + NSHELL(I1-NATREF)*
     +               (ESPV(L,I1,ISPIN)-ESPV(L,NREF,ISPIN))
   20     CONTINUE
        END IF
 
   30 CONTINUE
      END
