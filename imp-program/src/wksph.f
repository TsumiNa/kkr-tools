      SUBROUTINE WKSPH(DEN,IREF,NLST,NSHELL,TEXTS,PI,LMAX,NATPS,NATREF,
     +                 NEND,NSPIN)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER IEMXD
      PARAMETER (iemxd=150)
      INTEGER LMAXD
      PARAMETER (lmaxd=4)
C     ..
C     .. Local Scalars ..
      REAL*8 DSHC,DST,DSTC,DSTREL,DSTTOT
      INTEGER I1,IE,IL,IR,IS,L
C     ..
C     .. Local Arrays ..
      REAL*8 DD(NATYPD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG,REAL
C     ..
C     .. Scalar Arguments ..
      REAL*8 PI
      INTEGER LMAX,NATPS,NATREF,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DEN(IEMXD,0:LMAXD,NATYPD,NSPIND)
      INTEGER IREF(NTPERD),NLST(2),NSHELL(NTPERD)
      CHARACTER*13 TEXTS(3)
C     ..
      WRITE (6,FMT=9010)
      DO 30 I1 = 1,NATREF
        DSTTOT = 0.0D0
        DO 20 IS = 1,NSPIN
          IL = 2* (2-NSPIN) + IS
          IE = NLST(IS)
          DST = 0.0D0
          DO 10 L = 0,LMAX
            DST = DST - DIMAG(DEN(IE,L,I1,IS))
   10     CONTINUE
          DST = DST*REAL(NSPIN)/ (2.0D0*PI)/13.6056981D0
          DSTTOT = DSTTOT + DST
          WRITE (6,FMT=9020) I1,DST,TEXTS(IL)
   20   CONTINUE
        WRITE (6,FMT=9030) DSTTOT*2.357D0
        DD(I1) = DSTTOT
   30 CONTINUE
c
      DSTC = 0.0D0
      DSHC = 0.0D0
      DO 60 I1 = NATPS,NEND
        IR = IREF(I1-NATREF)
        DSTTOT = 0.0D0
        DO 50 IS = 1,NSPIN
          IL = 2* (2-NSPIN) + IS
          IE = NLST(IS)
          DST = 0.0D0
          DO 40 L = 0,LMAX
            DST = DST - DIMAG(DEN(IE,L,I1,IS))
   40     CONTINUE
          DST = DST*REAL(NSPIN)/ (2.0D0*PI)/13.6056981D0
          DSTTOT = DSTTOT + DST
          WRITE (6,FMT=9040) I1 - NATREF,DST,TEXTS(IL)
   50   CONTINUE
        DSTREL = DSTTOT/DD(IR) - 1.0D0
        WRITE (6,FMT=9050) DSTREL
        DSTC = DSTC + DSTREL*NSHELL(I1-NATREF)
        DSHC = DSHC + DD(IR)*NSHELL(I1-NATREF)
   60 CONTINUE
      WRITE (6,FMT=9060) DSTC
      WRITE (6,FMT=9000)
      RETURN


 9000 FORMAT (1x,82 ('-'))
 9010 FORMAT (/,1x,33 ('-'),' specific heat coefficients ',33 ('-'),/,
     +       1x,'>')
 9020 FORMAT (3x,i3,'-th. host atom : ldos at the fermi',' energy :',
     +       f12.6,' states / ev  for ',a13)
 9030 FORMAT (7x,'specific heat coefficient:',f12.6,' mj/mole/d*g**2',
     +       3 (/,1x, ('>')))
 9040 FORMAT (3x,i3,'-th. rep. atom : ldos at the fermi',' energy :',
     +       f12.6,' states / ev  for ',a13)
 9050 FORMAT (7x,'relative change of the ldos at the fermi energy :',
     +       f12.6,/,1x,'>')
 9060 FORMAT (1x,'>',/,3x,'relative change of the ldos at the fermi ',
     +       'energy in the cluster :',f12.6,' local summation ')
 9070 FORMAT (33X,'LDOS FROM EF EF(1)=',F10.6,'EF(2)=',F10.6)
      END
