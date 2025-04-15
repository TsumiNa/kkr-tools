c ************************************************************************
      SUBROUTINE ZHPEVY(IOPT,AP,W,Z,LDZ,N,AUX,NAUX)
      implicit none
c ************************************************************************
c  computes eigenvalues and, optionally, all eigenvectors of 
c  a hermitian matrix in upper/lower starage mode AP by using the 
c  LAPACK routine ZHPEV
c
c  (for compatibilty with LAPACK)
c
c  p.zahn, jan. 99
c ------------------------------------------------------------------------
c     .. scalar arguments
      INTEGER IOPT,LDZ,N,NAUX
c     .. array arguments 
      DOUBLE COMPLEX AP(*),AUX(*),Z(LDZ,*)
      DOUBLE PRECISION W(*)
C     .. locals
      INTEGER, PARAMETER :: NZHPEV = 100000
      DOUBLE PRECISION DAUX(NZHPEV)
      CHARACTER*1 JOBZ,UPLO
      INTEGER INFO
c ------------------------------------------------------------------------
      IF (4*N.GT. NZHPEV) THEN
        write(6,*) 'Increase the parameter NZHPEV in routine ZHPEVY'
        STOP 'ZHPEVY - 1'
      END IF
c
      IF (IOPT.EQ.0) THEN
        JOBZ = 'N'
        UPLO = 'L'
      ELSE IF (IOPT.EQ.1) THEN
        JOBZ = 'V'
        UPLO = 'L'
      ELSE IF (IOPT.EQ.20) THEN
        JOBZ = 'N'
        UPLO = 'U'
      ELSE IF (IOPT.EQ.21) THEN
        JOBZ = 'V'
        UPLO = 'U'
      ELSE 
        write(6,*) 'Parameter IOPT not in range [0,1,20,21].'
        STOP 'ZHPEVY -1'
      END IF
c
      CALL ZHPEV(JOBZ,UPLO,N,AP,W,Z,LDZ,AUX,DAUX,INFO)
c
      RETURN
      END
