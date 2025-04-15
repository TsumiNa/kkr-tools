c ************************************************************************
      INTEGER FUNCTION IOBEN(R)
      implicit none
c ************************************************************************
c
c                             --   --
c     Calculates the function |  r  |  (next upper or equal integer)
c                             |     |
c
c     Descrition of input parameters:
c
c       r : real number to look for
c
c                                           Rudolf Berrendorf, July 1992
c-----------------------------------------------------------------------
c
C     .. Scalar Arguments ..
      REAL R
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC FLOAT,INT
C     ..
      IF (FLOAT(INT(R)).EQ.R) THEN
        IOBEN = INT(R)

      ELSE
        IOBEN = INT(R+1.0)
      END IF

      RETURN

      END
