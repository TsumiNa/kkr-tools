c ************************************************************************
      INTEGER FUNCTION IHALFSPACE(I2,NATOM,NP)
      implicit none
c ************************************************************************
c     returns an index for the region in an interlayer system with 2D
c     periodicity
c     NATOM number of atoms in supercell
c     NP    number of atoms in principal layer
c           each half space (left/right) is represented 
c           by one princip. layer
c     I2    considered atomic index
c-----------------------------------------------------------------------
      INTEGER I2,NATOM,NP
c-----------------------------------------------------------------------
      IF ( I2.LE.NP ) THEN
c       left
        IHALFSPACE = 1
      ELSE IF ( I2.GT.NATOM-NP ) THEN
c       right
        IHALFSPACE = 2
      ELSE
c       interlayer region
        IHALFSPACE = 0
      END IF

      RETURN
      END
