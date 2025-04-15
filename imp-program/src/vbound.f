      SUBROUTINE VBOUND(EFERMI,EGR,IEN,KTE,NLST,NSTART,NEND,NSPIN,VBC,
     +                  VBCC,Z)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     determine the shift of the potentials .
c     the potentials calculated from the charge density use the
c     electrostatic zero as zero . in this program the green's
c     functions are related to the muffin tin zero . remember that
c     therefore the fermi energy is different for each spin direc-
c     tion - in the spin polarized case .
c     the calculated potentials have to be shifted to that energy
c     scale .
c     the constant potential shift vbcc is given by the input file .
c     it has to be determined self-consistently . in the spin-polarized
c     case the difference between the fermi energies (spin down - spin
c     up) is calculated and half of that difference is added to get
c     the shift of the spin down potential and subtracted to get the
c     one of the spin up potential .
c
c     in the case total energy calculation the energy zero is set
c     to be equal the electrostatic zero .
c     efermi is the fermi energy of the cristall respectively the
c     electrostatic zero .
c     calculate atomic charge times fermi energy to get the total
c     energies in the sense of the grand canonical ensemble .
c             (see notes by b.drittler)
c
c                               b.drittler   june 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER IEMXD
      PARAMETER (iemxd=150)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EFERMI,VBCC
      INTEGER KTE,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION EGR(*),IEN(IEMXD,*),VBC(*),Z(*)
      INTEGER NLST(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DEFERM
      INTEGER I1
C     ..
      VBCC = (VBC(1)+VBC(2))*0.5D0
      IF (NSPIN.EQ.2) THEN
        DEFERM = IEN(NLST(1),1) - IEN(NLST(2),2)
c
c---> all spin down potentials have to be shifted by vbc(1)
c
        VBC(1) = VBCC + 0.5D0*DEFERM
c
c---> all spin up potentials have to be shifted by vbc(2)
c
        VBC(2) = VBCC - 0.5D0*DEFERM

      ELSE
c
c---> paramagnetic case : all potentials have to be shifted by vbc(1)
c
        VBC(1) = VBCC
        VBC(2) = VBCC
      END IF
c
c---> determine value of the fermi-energy related to the
c               electrostatic zero
c
      EFERMI = IEN(NLST(1),1) - VBC(1)
c
c---> in case of total energy calculation calculate z times fermi energy
c
      IF (KTE.EQ.1) THEN
c
c---> loop over reference atoms
c
        DO 10 I1 = NSTART,NEND
          EGR(I1) = Z(I1)*EFERMI
   10   CONTINUE
      END IF

      END
