c ************************************************************************
      SUBROUTINE ZGEEVY(IOPT,A,LDA,W,Z,LDZ,MM,SELECT,N,
     +                  AUX,NAUX,TAU,IFAIL)
      implicit none
c ************************************************************************
c  computes eigenvalues and (no or selected or all) eigenvectors of 
c  a general matrix A by using the LAPACK routines
c  ZGEHRD   transform A to an upper Hessenberg form
c  ZHSEQR   calculate only eigenvalues
c  ZHSEIN   calculate eigenvalues and selected or all eigenvectors
c  ZUNMHR   transform eigenvectors ( AC=l*C, H=Q'AQ, Hc=l*c, C=Qc)
c
c  p.zahn, jan. 99
c ------------------------------------------------------------------------
c     .. scalar arguments
      INTEGER IOPT,LDA,LDZ,MM,N,NAUX
c     .. array arguments 
      INTEGER IFAIL(*)
      DOUBLE COMPLEX A(LDA,*),AUX(*),TAU(*),W(*),Z(LDZ,*)
      LOGICAL SELECT(*)
C     .. LOCAL
      Integer ,PARAMETER :: NZGEEV = 100000
      DOUBLE PRECISION DAUX(NZGEEV)
      DOUBLE COMPLEX ZDUMMY
      INTEGER I,IDUMMY,INFO,M
c ------------------------------------------------------------------------
      IF (4*N .GT. NZGEEV) THEN
        write(6,*) 'Increase the parameter NZGEEV in routine ZHPEVY'
        STOP 'ZGEEVY - 1'
      END IF
c
      IF (IOPT.EQ.0) THEN
C       only eigenvalues
C       A -> Hessenberg form
        CALL ZGEHRD(N,1,N,A,LDA,TAU,AUX,NAUX,INFO)
c       all eigenvalues
        CALL ZHSEQR('E','N',N,1,N,A,LDA,W,ZDUMMY,1,AUX,NAUX,INFO)
      ELSE                          ! (IOPT.EQ.0)
        IF (IOPT.EQ.1) THEN
C         eigenvalues and all right eigenvectors
          CALL ZGEEV('N','V',N,A,LDA,W,ZDUMMY,1,Z,LDZ, 
     +         AUX,NAUX,DAUX,INFO)
        ELSE IF (IOPT.EQ.2) THEN
C         selected right eigenvectors
c         eigenvalues have to be provided in array W
          M = 0
          DO 20 I=1,N
            IF (SELECT(I) .eqv. .TRUE.) M = M + 1
 20       END DO
C         A -> Hessenberg form
          CALL ZGEHRD(N,1,N,A,LDA,TAU,AUX,NAUX,INFO)
c         selected eigenvectors
          CALL ZHSEIN('R','Q','N',SELECT,N,A,LDA,W,
     +         ZDUMMY,1,Z,LDZ,MM,M,AUX,DAUX,IDUMMY,IFAIL,INFO)

c          write(6,*) 'ZHSEIN: IFAIL,INFO',IFAIL,INFO

          CALL ZUNMHR('L','N',N,M,1,N,A,LDA,TAU,Z,LDZ,AUX,NAUX,INFO)

c          write(6,*) 'ZUNMHR: INFO',INFO

        END IF

      END IF                        ! (IOPT.EQ.0)
      RETURN
      END
