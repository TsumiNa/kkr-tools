      SUBROUTINE GSTRUCT(RIMP,NTYPIM,NTYPHO,NIMP,IGROUP,LGROUP,NREP,
     +                   LMAXSQ,LGSYMM,LBTSYM,ihandle)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NIMPD,NSHELD
      PARAMETER (NIMPD=4290,NSHELD=4020)
      INTEGER NATMX
      PARAMETER (NATMX=NIMPD)
      INTEGER NUMCOL
      PARAMETER (NUMCOL=3)
      INTEGER NSEC,NBASIS,NISTR,NSTR,NMSTR
      PARAMETER (NSEC=700,NBASIS=3000,NISTR=2600000,NSTR=2600000,
     +          NMSTR=1800000)
      INTEGER LMAXD
      PARAMETER (LMAXD=4)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER IJD
      PARAMETER (IJD=434)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAXSQ,NIMP,NREP,ihandle
      LOGICAL LBTSYM,LGSYMM
      CHARACTER*5 IGROUP,LGROUP
C     ..
C     .. Array Arguments ..
      REAL*8 RIMP(3,NIMPD)
      INTEGER NTYPHO(NIMPD),NTYPIM(NIMPD)
C     ..
C     .. Local Arrays ..
      REAL*8 CN(NBASIS,NSEC,NUMCOL),DROT(48,LMMAXD,LMMAXD),
     +                 RIJ(IJD,3),STR(NSTR),WTYR(IJD,LMMAXD),
     +                 YR(IJD,LMMAXD)
      INTEGER IJ(NSHELD),IMAX(NSEC,NATMX,NUMCOL),
     +        IMIN(NSEC,NATMX,NUMCOL),IROIND(48),IROTC2V(4),IROTC3V(6),
     +        IROTC4V(8),IROTD2D(8),IROTD(24),IROTMN(NIMPD,0:NIMPD),
     +        ISHELL(NIMPD,0:NIMPD),ISTR(NISTR),JSHELL(5,NSHELD),
     +        JTHELL(5,NIMPD),KSHELL(5,NSHELD),MAXSTR(NMSTR),
     +        MI(NBASIS,NSEC,NUMCOL),MINSTR(NMSTR),MSTR(NMSTR),
     +        N0ATOM(NBASIS,NSEC,NUMCOL),ND(48,3,3),IROTD4H(16)
      integer irotd2h(8)
C     ..
C     .. External Subroutines ..
      EXTERNAL BTSYMM,GSYMM,ROTCUB,ROTMAT,SPHERE
C     ..
C     .. Local Scalars ..
      REAL*8 E2,R1,R2,R3
      INTEGER I,IJEND,IJMAX,IROT,ISPIN,J,J1,J1ROT,J2,J2ROT,J3,J3ROT,J4,
     +        J5,JABS,JROT,LMAXIN,M,MS,MSHELL,N,NEZ,NHSPIN,NIROT,NLROT,
     +        NS,NSHELL,NTHELL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC IABS,MAX
C     ..
C     .. Data statements ..
      DATA IROIND/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
     +     22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,
     +     42,43,44,45,46,47,48/
      DATA IROTD/1,2,3,4,5,6,7,8,9,16,17,18,34,35,36,37,38,39,43,44,45,
     +     46,47,48/
      DATA IROTC4V/1,10,12,18,40,41,45,46/
      DATA IROTC3V/1,2,3,46,47,48/
      DATA IROTC2V/1,21,42,46/
      DATA IROTD2h/1,18,21,22,25,42,45,46/
      DATA IROTD2D/1,16,17,18,34,36,45,46/
      DATA IROTD4H/1,10,12,16,17,18,21,22,25,34,36,40,41,42,45,46/
      !DATA IROTE/1/
C     ..
      NLROT = 48
      IF (LGROUP.EQ.'Td   ') THEN
        NLROT = 24
        DO 10 JROT = 1,NLROT
          IROIND(JROT) = IROTD(JROT)
   10   CONTINUE
      END IF
      IF (LGROUP.EQ.'C3v   ') THEN
        NLROT = 6
        DO 11 JROT = 1,NLROT
          IROIND(JROT) = IROTC3V(JROT)
 11    CONTINUE
      END IF
      IF (LGROUP.EQ.'D4h   ') THEN
        NLROT = 16
        DO 12 JROT = 1,NLROT
          IROIND(JROT) = IROTD4H(JROT)
 12    CONTINUE
      END IF
      CALL ROTCUB(ND)
      CALL SPHERE(LMAXD,YR,WTYR,RIJ,IJEND,IJD)
      CALL ROTMAT(LMAXD,ND,DROT,YR,WTYR,RIJ,IJEND)
      NSHELL = 0
      DO 60 M = 1,NIMP
        DO 50 N = 1,NIMP
          J1 = 1000000* (RIMP(1,N)-RIMP(1,M))
          J2 = 1000000* (RIMP(2,N)-RIMP(2,M))
          J3 = 1000000* (RIMP(3,N)-RIMP(3,M))
          J4 = NTYPHO(N)
          J5 = NTYPHO(M)
c--> rotate and compare with previous shells
          DO 30 JROT = 1,NLROT
            IROT = IROIND(JROT)
            J1ROT = ND(IROT,1,1)*J1 + ND(IROT,1,2)*J2 + ND(IROT,1,3)*J3
            J2ROT = ND(IROT,2,1)*J1 + ND(IROT,2,2)*J2 + ND(IROT,2,3)*J3
            J3ROT = ND(IROT,3,1)*J1 + ND(IROT,3,2)*J2 + ND(IROT,3,3)*J3
            DO 20 NS = 1,NSHELL
              JABS = IABS(J1ROT-JSHELL(1,NS)) +
     +               IABS(J2ROT-JSHELL(2,NS)) +
     +               IABS(J3ROT-JSHELL(3,NS)) + IABS(J4-JSHELL(4,NS)) +
     +               IABS(J5-JSHELL(5,NS))
              IF (JABS.EQ.0) THEN
                ISHELL(M,N) = NS
                IROTMN(M,N) = IROT
                GO TO 40
              END IF
   20       CONTINUE
   30     CONTINUE
          NSHELL = NSHELL + 1
          ISHELL(M,N) = NSHELL
          IROTMN(M,N) = 1
          JSHELL(1,NSHELL) = J1
          JSHELL(2,NSHELL) = J2
          JSHELL(3,NSHELL) = J3
          JSHELL(4,NSHELL) = J4
          JSHELL(5,NSHELL) = J5
          R1 = J1*0.000001D0
          R2 = J2*0.000001D0
          R3 = J3*0.000001D0
          WRITE (6,FMT=9000) NSHELL,N,M,R1,R2,R3,J4,J5
          WRITE (9,FMT=9000) NSHELL,N,M,R1,R2,R3,J4,J5
   40     CONTINUE
C         WRITE(6,*) M,N,ISHELL(M,N),IROTMN(M,N)
   50   CONTINUE
   60 CONTINUE
      N = 0
      WRITE (9,FMT=9000) N
      WRITE (6,FMT=*) 'impurity information'
      WRITE (9,FMT=*) 'impurity information'
      DO 70 N = 1,NIMP
        WRITE (6,FMT=9020) N, (RIMP(J,N),J=1,3),NTYPHO(N),NTYPIM(N)
        WRITE (9,FMT=9020) N, (RIMP(J,N),J=1,3),NTYPHO(N),NTYPIM(N)
   70 CONTINUE
      CLOSE (9)
      IF (LGSYMM) THEN
        READ (15) NEZ,NHSPIN,MSHELL,LMAXIN
        DO N = 1,MSHELL
          READ (15) R1,R2,R3,J4,J5
      write(6,1222) R1,r2,r3,J4,J5
 1222 format(1x,' gfstruct numdiff',3f10.3,2I5)
          KSHELL(1,N) = R1*1000000
          KSHELL(2,N) = R2*1000000
          KSHELL(3,N) = R3*1000000
          KSHELL(4,N) = J4
          KSHELL(5,N) = J5
        END DO
        WRITE (6,FMT=*) NEZ,NHSPIN,MSHELL,LMAXIN
        IJMAX = 0
        DO 90 NS = 1,NSHELL
          IJ(NS) = 0
          DO 80 MS = 1,MSHELL
            IF (KSHELL(1,MS).EQ.JSHELL(1,NS) .AND.
     +          KSHELL(2,MS).EQ.JSHELL(2,NS) .AND.
     +          KSHELL(3,MS).EQ.JSHELL(3,NS) .AND.
     +          KSHELL(4,MS).EQ.JSHELL(4,NS) .AND.
     +          KSHELL(5,MS).EQ.JSHELL(5,NS)) THEN
              IJ(NS) = MS
              WRITE (6,FMT=*) 'MS,KSHELL,NS,JSHELL',MS,
     +          (KSHELL(I,MS),I=1,5),NS, (JSHELL(I,NS),I=1,5)
            END IF
   80     CONTINUE
          IJMAX = MAX(IJMAX,IJ(NS))
          WRITE (6,FMT=*) 'IJ(NS),NS,IJMAX',IJ(NS),NS,IJMAX
          IF (IJ(NS).EQ.0) THEN
            WRITE (6,FMT=*) ' ERROR IN GSTRUCT'
C         STOP '57'
          END IF
   90   CONTINUE
        DO 100 ISPIN = 1,NHSPIN
           WRITE (80) NEZ,NHSPIN,E2,E2
c          WRITE (81,*) NEZ,NHSPIN,E2,E2
c          Call fxdrint(ihandle,nez,1)
c       write(ihandle) nez
c          Call fxdrint(ihandle,nhspin,1)
c       write(ihandle) nhspin
c          Call fxdrdbl(ihandle,e2,1)
choshino  e2=0.d0  ?????
        e2=0.D0
c      write(ihandle) e2
c          Call fxdrdbl(ihandle,e2,1)
c       write(ihandle) e2
       write(6,*) 'before gsymm'
          CALL GSYMM(LMAXSQ,NEZ,DROT,IROTMN,ISHELL,NIMP,NSHELL,MSHELL,
     +               IJ,NHSPIN,CN,STR,IMAX,IMIN,ISTR,MAXSTR,MI,MINSTR,
     +               MSTR,N0ATOM,ihandle)
      write(6,*) 'after gsymm'
  100   CONTINUE
      END IF
      CLOSE (80)
      NIROT = 48
      IF (IGROUP.EQ.'E    ') THEN
      IROIND(1) = 1
      ENDIF
      IF (IGROUP.EQ.'Td   ') THEN
        NIROT = 24
        DO 110 JROT = 1,NIROT
          IROIND(JROT) = IROTD(JROT)
  110   CONTINUE
      END IF
      IF (IGROUP.EQ.'C4v  ') THEN
        NIROT = 8
        DO 120 JROT = 1,NIROT
          IROIND(JROT) = IROTC4V(JROT)
  120   CONTINUE
      END IF
      IF (IGROUP.EQ.'C3v  ') THEN
        NIROT = 6
        DO 130 JROT = 1,NIROT
          IROIND(JROT) = IROTC3V(JROT)
  130   CONTINUE
      END IF
      IF (LGROUP.EQ.'D4h   ') THEN
        NIROT = 16
        DO 142 JROT = 1,NLROT
          IROIND(JROT) = IROTD4H(JROT)
 142    CONTINUE
      END IF
      IF (IGROUP.EQ.'C2v  ') THEN
        NIROT = 4
        DO 140 JROT = 1,NIROT
          IROIND(JROT) = IROTC2V(JROT)
  140   CONTINUE
      END IF
      IF (IGROUP.EQ.'D2h  ') THEN
        NIROT = 8
        DO 143 JROT = 1,NIROT
          IROIND(JROT) = IROTD2h(JROT)
  143   CONTINUE
      write(6,1111) NIROT
      write(6,1112) (IROIND(JROT),JROT=1,NIROT)
 1111  format(1x,' nirot',i5)
 1112  format(1x,' iroind',8I5)
      END IF
      IF (IGROUP.EQ.'D2d  ') THEN
        NIROT = 8
        DO 141 JROT = 1,NIROT
          IROIND(JROT) = IROTD2D(JROT)
  141   CONTINUE
      END IF
      NTHELL = 0
      DO 180 M = 1,NIMP
        J1 = 1000000*RIMP(1,M)
        J2 = 1000000*RIMP(2,M)
        J3 = 1000000*RIMP(3,M)
        J4 = NTYPIM(M)
        J5 = NTYPIM(M)
c--> rotate and compare with previous shells
        DO 160 JROT = 1,NIROT
          IROT = IROIND(JROT)
          J1ROT = ND(IROT,1,1)*J1 + ND(IROT,1,2)*J2 + ND(IROT,1,3)*J3
          J2ROT = ND(IROT,2,1)*J1 + ND(IROT,2,2)*J2 + ND(IROT,2,3)*J3
          J3ROT = ND(IROT,3,1)*J1 + ND(IROT,3,2)*J2 + ND(IROT,3,3)*J3
          DO 150 NS = 1,NTHELL
            JABS = IABS(J1ROT-JTHELL(1,NS)) + IABS(J2ROT-JTHELL(2,NS)) +
     +             IABS(J3ROT-JTHELL(3,NS)) + IABS(J4-JTHELL(4,NS)) +
     +             IABS(J5-JTHELL(5,NS))
            IF (JABS.EQ.0) THEN
              ISHELL(M,0) = NS
              IROTMN(M,0) = IROT
              GO TO 170
            END IF
  150     CONTINUE
  160   CONTINUE
        NTHELL = NTHELL + 1
        ISHELL(M,0) = NTHELL
        IROTMN(M,0) = 1
        JTHELL(1,NTHELL) = J1
        JTHELL(2,NTHELL) = J2
        JTHELL(3,NTHELL) = J3
        JTHELL(4,NTHELL) = J4
        JTHELL(5,NTHELL) = J5
        R1 = J1*0.000001D0
        R2 = J2*0.000001D0
        R3 = J3*0.000001D0
        WRITE (6,FMT=9010) NTHELL,M,R1,R2,R3,J4,J5
  170   CONTINUE
  180 CONTINUE
      CALL BTSYMM(LMAXSQ,ND,RIMP,NTYPIM,NTYPHO,NIMP,JTHELL,NTHELL,NREP,
     +            IGROUP,LBTSYM,CN,STR,IMAX,IMIN,ISTR,MAXSTR,MI,MINSTR,
     +            MSTR,N0ATOM)
      RETURN
 9000 FORMAT (3I5,3F10.6,2I5)
 9010 FORMAT (2I5,3F10.6,2I5)
 9020 FORMAT (I5,3F10.6,2I5)
      END
