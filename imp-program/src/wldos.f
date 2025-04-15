      SUBROUTINE WLDOS(DEN,EF,IEN,ITITLE,PI,IA,IELAST,LMAX,NATYP,NPTPS,
     +                 NSPIN)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER LMAXD
      PARAMETER (lmaxd=4)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Local Scalars ..
      REAL*8 EFCTOR,EFEV,IAN
      INTEGER I1,IE,IPOT,IS,L
      CHARACTER*8 DOSFL1
      CHARACTER*9 DOSFL0
      CHARACTER*10 DOSFL
C     ..
C     .. Local Arrays ..
      REAL*8 WAB(NPOTD,IEMXD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG
C     ..
C     .. Scalar Arguments ..
      REAL*8 PI
      INTEGER IA,IELAST,LMAX,NATYP,NPTPS,NSPIN
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DEN(IEMXD,0:LMAXD,NATYPD,NSPIND)
      REAL*8 EF(NSPIND),IEN(IEMXD,NSPIND)
      INTEGER ITITLE(20,NPOTD)
C     ..
      EFCTOR = 1.0D0/13.6058D0
      DOSFL0 = 'dos.ipot0'
      DOSFL1 = 'dos.ipot'
c
      DO 50 IS = 1,NSPIN
        DO 40 I1 = 1,NATYP
          IPOT = NSPIN* (I1-1) + IS
          REWIND 47
          IF (IPOT.LT.10) WRITE (47,FMT='(A9,I1)') DOSFL0,IPOT
          IF (IPOT.GE.10) WRITE (47,FMT='(A8,I2)') DOSFL1,IPOT
          REWIND 47
          READ (47,FMT='(A10)') DOSFL
          OPEN (48,FILE=DOSFL,FORM='formatted')
          EFEV = IEN(IELAST,IS)/EFCTOR
          WRITE (48,FMT=9000) (ITITLE(IA,IPOT),IA=1,20)
          WRITE (48,FMT=9010) NSPIN,IELAST,EFEV,EFCTOR
          DO 30 IE = 1,IELAST
            IAN = IEN(IE,IS) - EF(IS)
            WAB(IPOT,IE) = 0.0D0
            WRITE (48,FMT=9020) IAN, (-1.0D0/PI*DIMAG(DEN(IE,L,I1,IS)),
     +        L=0,LMAX)
            DO 10 L = 0,LMAX
              WAB(IPOT,IE) = WAB(IPOT,IE) -
     +                       1.0D0/PI*DIMAG(DEN(IE,L,I1,IS))
   10       CONTINUE
   20       CONTINUE
   30     CONTINUE
          CLOSE (48)
   40   CONTINUE
   50 CONTINUE
c
      WRITE (6,FMT=9030) EF(1),EF(2)
      OPEN (48,FILE='dos.ipot00',FORM='formatted')
      WRITE (48,FMT=9000) (ITITLE(IA,NPTPS),IA=1,20)
      WRITE (48,FMT=9010) NSPIN,IELAST,IEN(IELAST,1),EFCTOR
      DO 70 IS = 1,NSPIN
        DO 60 IE = 1,IELAST
          IAN = IEN(IE,IS) - EF(IS)
          WRITE (48,FMT=9020) IAN, (WAB(NSPIN* (I1-1)+IS,IE),I1=1,
     +      NATYP)
   60   CONTINUE
   70 CONTINUE
      CLOSE (48)
      RETURN


 9000 FORMAT (20a4)
 9010 FORMAT (2i5,2f10.6)
 9020 FORMAT (1p,7D11.3)
 9030 FORMAT (33X,'LDOS FROM EF EF(1)=',F10.6,'EF(2)=',F10.6)
      END
