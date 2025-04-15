      SUBROUTINE LLOYDG(DNE,DNEINT,DF,E,IE,IELAST,NSPIN,NDG,NP,NREP,DET)
      IMPLICIT NONE
c----------------------------------------------------------------------
c     calculates the change of the density of states resulting from the
c     back-scattering contribution via Lloyd's formula for each
c     representation.  ln det ( 1 - g dt )
c                                           b.drittler oct. 1987
c     modified  by  R. Zeller      Sept. 1994
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NREPD
      PARAMETER (NREPD=4)
      COMPLEX*16 CZERO,EI
      PARAMETER (CZERO= (0.0D0,0.0D0),EI= (0.0D0,1.0D0))
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DET,DF,DNE,E
      REAL*8 DNEINT
      INTEGER IE,IELAST,NDG,NP,NREP,NSPIN
C     ..
C     .. Local Scalars ..
      COMPLEX*16 DNENEW,DNEOLD,EONE,ETWO,FEST
      REAL*8 COR,FAC
C     ..
C     .. Local Arrays ..
      COMPLEX*16 FONE(NREPD),FTWO(NREPD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,CDLOG,DATAN,DIMAG
C     ..
C     .. Save statement ..
      SAVE
C     ..
      FAC = 1.D0/ (8.D0*DATAN(1.D0))
c
c--> initialize energy points and values to zero
c
      IF (IE.EQ.1) THEN
        DNEOLD = CZERO
        EONE = CZERO
        ETWO = CZERO
        FONE(NP) = CZERO
        FTWO(NP) = CZERO
      END IF

      DNENEW = -CDLOG(DET)*FAC
C
c---> correct jumps caused by phase shifts
c
      IF (IE.LE.2) THEN
        COR = ANINT(DIMAG(DNENEW-DNEOLD))

      ELSE
c
c---> two point lagrange interpolation
c
        FEST = (FTWO(NP)* (E-EONE)-FONE(NP)* (E-ETWO))/ (ETWO-EONE)
        COR = ANINT(DIMAG(DNENEW-FEST))
      END IF
c
      DNEOLD = DNENEW - EI*COR
c

      FONE(NP) = FTWO(NP)
      FTWO(NP) = DNEOLD
      DNE = DNEOLD*4.D0/NSPIN
      DNEINT = DNEINT + NDG*DIMAG(DF*DNE)/ (FAC*4.D0)*NSPIN

c
c---> update energy points
c
      IF (NP.EQ.NREP) THEN
        EONE = ETWO
        ETWO = E
      END IF

      END
