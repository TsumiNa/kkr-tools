c ************************************************************************
      SUBROUTINE GREFREAD(NATOM,GINP,ITMAT)
      implicit none
c ************************************************************************
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER ITMAT,NATOM
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GINP(LMAXSQ*NACLSD,*)
c     .. local scalars
      INTEGER M,N
c      .. external statement
      LOGICAL TEST
      EXTERNAL RCSTOP,TEST
      SAVE
C     ..
      READ (ITMAT) N
      IF (N.NE.NATOM) THEN 
        write(6,*) 'natom .ne. n',natom,n
        CALL RCSTOP('GREFREAD')
      END IF
      READ (ITMAT) ((GINP(N,M),M=1,LMAXSQ),N=1,LMAXSQ*NACLSD)
      RETURN

      END
