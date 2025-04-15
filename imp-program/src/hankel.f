      SUBROUTINE HANKEL(H,L,ARG)
      IMPLICIT NONE
c  this subroutine uses the explicit formulas for the hankel
c  functions. for higher l-values these formulas may lead to
c  loss of significant figures.
C     .. Scalar Arguments ..
      COMPLEX*16 ARG
      INTEGER L
C     ..
C     .. Array Arguments ..
      COMPLEX*16 H(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 A1,A2,A3,A4,I
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC CDEXP
C     ..
C     .. Save statement ..
      SAVE I
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Data statements ..
      DATA I/ (0.D0,1.D0)/
C     ..

      H(1) = -CDEXP(ARG*I)/ARG
      IF (L.NE.1) THEN
        A1 = (1.D0,0.D0) - ARG*I
        H(2) = H(1)*A1/ARG
        IF (L.NE.2) THEN
          A1 = 3.D0*A1
          A2 = ARG*ARG
          H(3) = H(1)* (A1-A2)/A2
          IF (L.NE.3) THEN
            A1 = 5.D0*A1
            A3 = A2*ARG*I
            A4 = A2*ARG
            A2 = 6.D0*A2
            H(4) = H(1)* (A1-A2+A3)/A4
            IF (L.NE.4) THEN
              A1 = 7.D0*A1
              A2 = 7.5D0*A2
              A3 = 10.D0*A3
              A4 = A4*ARG
              H(5) = H(1)* (A1-A2+A3+A4)/A4
              IF (L.NE.5) THEN
                H(6) = (9.0D0,0.0D0)*H(5)/ARG - H(4)
                IF (L.NE.6) THEN
                  WRITE (6,FMT=9000) L
                  CALL RCSTOP('34      ')

                END IF

              END IF

            END IF

          END IF

        END IF

      END IF

      RETURN


 9000 FORMAT (2x,' hankel :  l=',i2,' is too large')
      END
