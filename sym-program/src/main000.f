      PROGRAM STRUCTURE_USE
      IMPLICIT NONE
c     include 'cxml_include.f90'
C     .. Parameters ..
      INTEGER LMAXD,NBASD,NIMPD
      PARAMETER (LMAXD=4,NBASD=32,NIMPD=4290)
C     ..
C     .. Local Scalars ..
      REAL*8 ALAT,GMAX,RMAX
      INTEGER I,INVER,IPUN,KTAU,LMAX,LMAXSQ,LPOT,N,NBASIS,NC,NCOUNT,
     +        NIMP,NIN,NREP,NSHELL,NTYP,ihandle
      LOGICAL LBTSYM,LGSYMM,LHARMY
      CHARACTER*5 IGROUP,LGROUP,LSTRUCT,SKIPLINE
      CHARACTER*8 AUNITS
      CHARACTER*26 GSYMOUT
C     ..
C     .. Local Arrays ..
      REAL*8 RBASIS(3,NBASD),RIMP(3,NIMPD),RLATT(3,3)
      INTEGER NTYPHO(NIMPD),NTYPIM(NIMPD)
C     ..
C     .. External Subroutines ..
      EXTERNAL GSTRUCT,HARMONY
C     ..
c      CALL FXDRINI
      OPEN (9,FILE='gf_input',FORM='formatted')
      OPEN (10,FILE='structure_input',FORM='formatted')
      OPEN (11,FILE='lattice_input',FORM='formatted')
      OPEN (12,FILE='harmony_input',FORM='formatted')
      OPEN (15,FILE='gf_indata',FORM='unformatted')
      READ (10,FMT='(5X,I2,3(8X,L1))') LMAX,LHARMY,LGSYMM,LBTSYM
      READ (10,FMT='(A26)') GSYMOUT

      ihandle = 80
c     Call fxdropn('greenout_fxdr','ENCODE',ihandle)
      OPEN (80,FILE='greenout_win',FORM='unformatted',access='stream')
      READ (10,FMT=9020) LSTRUCT
      WRITE (9,FMT=9030) LSTRUCT
      READ (10,FMT=9040) LGROUP
      READ (10,FMT=*) ALAT,AUNITS
      READ (10,FMT=*) (RLATT(I,1),I=1,3)
      READ (10,FMT=*) (RLATT(I,2),I=1,3)
      READ (10,FMT=*) (RLATT(I,3),I=1,3)
      READ (10,FMT=*) NBASIS
      DO 10 N = 1,NBASIS
        READ (10,FMT=*) (RBASIS(I,N),I=1,3)
   10 CONTINUE
      READ (10,FMT=9020) SKIPLINE
      READ (10,FMT=9040) IGROUP
      READ (10,FMT=*) NIMP
      DO 20 N = 1,NIMP
        READ (10,FMT=*) (RIMP(I,N),I=1,3),NIN,NTYPIM(N),NTYPHO(N)
   20 CONTINUE
C--->  tests for input data follow
      IF (AUNITS.EQ.'ANGSTROM') THEN
        ALAT = ALAT/0.529177
        AUNITS = 'BOHR    '
      END IF
      IF (AUNITS.NE.'BOHR    ') THEN
        WRITE (6,FMT=*) AUNITS,
     +    ': Wrong unit for lattice structure given'
        STOP
      END IF
      IF (IGROUP.NE.'E    ' .AND. IGROUP.NE.'C2v  ' .AND.
     +    IGROUP.NE.'C4v  ' .AND. IGROUP.NE.'C3v  ' .AND.
     +    IGROUP.NE.'Td   ' .AND. IGROUP.NE.'Oh   ' .AND.
     +    IGROUP.NE.'D2d  '. AND. IGROUP.NE.'D4h  ' .and.
     +    IGROUP.NE.'D2h  ') THEN
        WRITE (6,FMT=*) IGROUP,
     +    ': This point group symmetry is not yet allowed'
        STOP
      END IF
      IF (LGROUP.NE.'Td   ' .AND. LGROUP.NE.'Oh   '.
     +    AND. LGROUP.NE.'C3v   '. AND. IGROUP.NE.'D4h  ') 
     +  THEN
        WRITE (6,FMT=*) LGROUP,
     +    ': This lattice structure is not provided'
c        STOP
      END IF
C--->  output for lattice.fortran
      IF (LGROUP.EQ.'Td   ') THEN
        RMAX = 5.0
        GMAX = 80.0
      END IF
      LMAXSQ = (LMAX+1)**2
      LPOT = LMAX*2
      WRITE (11,FMT=9050) LPOT,ALAT,NBASIS,RMAX,GMAX
      DO 30 N = 1,3
        WRITE (11,FMT=9060) (RLATT(I,N),I=1,3)
   30 CONTINUE
      DO 40 N = 1,NBASIS
        WRITE (11,FMT=9070) (RBASIS(I,N),I=1,3)
   40 CONTINUE
      CLOSE (11)
C--->  output for harmony.fortran
      IF (IGROUP.EQ.'Oh   ') THEN
        NREP = 10
        KTAU = 48
        INVER = 24
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'Td   ') THEN
        NREP = 5
        KTAU = 24
        INVER = 12
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'C4v  ') THEN
        NREP = 5
        KTAU = 8
        INVER = 4
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'C3v  ') THEN
        NREP = 3
        KTAU = 6
        INVER = 3
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'C2v  ') THEN
        NREP = 4
        KTAU = 4
        INVER = 2
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'D2h  ') THEN
        NREP = 8
        KTAU = 8
        INVER = 4
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'E    ') THEN
        NREP = 1
        KTAU = 1
        INVER = 1
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'D2d  ') THEN
        NREP = 5
        KTAU = 8
        INVER = 4
        IPUN = 3
      END IF
      IF (IGROUP.EQ.'D4h  ') THEN
        NREP = 10
        KTAU = 16
        INVER = 8
        IPUN = 3
      END IF
      NTYP = NTYPIM(1)
      NCOUNT = 1
      NSHELL = 0
      WRITE (12,FMT=9000) IGROUP
      WRITE (12,FMT=9010) NTYPIM(NIMP),LMAX,NREP
      DO 70 N = 2,NIMP
        IF (NTYPIM(N).NE.NTYP) THEN
          NSHELL = NSHELL + 1
          WRITE (12,FMT=9080) KTAU,NCOUNT,INVER,IPUN
          DO 50 NC = N - NCOUNT,N - 1
            WRITE (12,FMT=9090) (RIMP(I,NC),I=1,3),NTYPHO(NC),NC,
     +        NTYPIM(NC)
   50     CONTINUE
          NTYP = NTYPIM(N)
          NCOUNT = 1
        ELSE
          NCOUNT = NCOUNT + 1
        END IF
        IF (N.EQ.NIMP) THEN
          WRITE (12,FMT=9080) KTAU,NCOUNT,INVER,IPUN
          DO 60 NC = N + 1 - NCOUNT,N
            WRITE (12,FMT=9090) (RIMP(I,NC),I=1,3),NTYPHO(NC),NC,
     +        NTYPIM(NC)
   60     CONTINUE

        END IF
   70 CONTINUE
      REWIND 12
      write(6,*) ' before harmony'
      IF (LHARMY) CALL HARMONY
      write(6,*) ' after harmony'
      IF (LMAX.GT.LMAXD) STOP 'LMAX'
      REWIND 19
c     CLOSE (12,STATUS='delete')
      CLOSE (12)
      OPEN (33,FILE='back.coeff',FORM='unformatted')
      OPEN (35,FILE='tu.coeff',FORM='formatted')
      write(6,*) ' before gstruct'
      CALL GSTRUCT(RIMP,NTYPIM,NTYPHO,NIMP,IGROUP,LGROUP,NREP,LMAXSQ,
     +             LGSYMM,LBTSYM,ihandle)
      write(6,*) ' after gstruct'
      CLOSE (19)
      STOP
 9000 FORMAT ('pointop/',A5)
 9010 FORMAT (3I10)
 9020 FORMAT (15X,A5)
 9030 FORMAT (A5)
 9040 FORMAT (A5)
 9050 FORMAT (I4,46X,'lpot=2*lmax',/,F10.5,I6,2F10.5,14X,
     +       'alat,nbasis,rmax,gmax')
 9060 FORMAT (3F10.6,20X,'primitive vectors')
 9070 FORMAT (3F10.6,20X,'basis vectors')
 9080 FORMAT (5X,3I5,25X,I5)
 9090 FORMAT (3F10.6,I2,I4,I5)
      END
