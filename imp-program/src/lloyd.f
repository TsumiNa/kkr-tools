      SUBROUTINE LLOYD(ALPHA,DETALF,IREF,NDG,NSHELL,TEXTS,DET,DF,E,
     +                 DNEINT,IE,IELAST,ISPIN,NATPS,NATREF,NEND,NP,NREP,
     +                 NSPIN,LCORE,NCORE,KSYMMAT,KESYM,DTB1)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER IEMXD,NREPD
      PARAMETER (iemxd=150,NREPD=4)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
C     ..
CASYM LLOYD
      INTEGER KSYMMAT,KESYM
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,*)
C     ..
C     .. Local Scalars ..
      REAL*8 DNMOM,DNTOT,ZERO
      INTEGER I1,IS,L,NR
C     ..
C     .. External Subroutines ..
      EXTERNAL LLOYDA,LLOYDG
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DET,DF,E
      REAL*8 DNEINT
      INTEGER IE,IELAST,ISPIN,NATPS,NATREF,NEND,NP,NREP,NSPIN
C     ..
C     .. Array Arguments ..
      COMPLEX*16 ALPHA(LMX,NATYPD),DETALF(NATYPD,IEMXD,NSPIND)
      INTEGER IREF(NTPERD),LCORE(20,NPOTD),NCORE(NPOTD),NDG(NREPD),
     +        NSHELL(NTPERD)
      CHARACTER*13 TEXTS(3)
C     ..
C     .. Local Arrays ..
      COMPLEX*16 DNE(NSPIND,NREPD),DNEALF(0:LMX,NSPIND,NATYPD)
C     ..
C     .. Data statements ..
      DATA ZERO/0.0D0/
C     ..
      IF (NP.EQ.1) THEN
        IF (IE.EQ.1 .AND. ISPIN.EQ.1) DNEINT = ZERO
c
        IF (IE.EQ.IELAST) THEN
          IF (ISPIN.EQ.1) THEN
            DNTOT = ZERO
            DNMOM = ZERO
          END IF
c
          IF (ISPIN.EQ.2) DNMOM = -DNMOM
        END IF
      END IF

      CALL LLOYDG(DNE(ISPIN,NP),DNEINT,DF,E,IE,IELAST,NSPIN,NDG(NP),NP,
     +            NREP,DET)
      IF (NP.EQ.1) THEN
        DO 20 I1 = NATPS,NEND
          CALL LLOYDA(DNEALF(0,ISPIN,I1),I1,DNEINT,DF,E,ALPHA,IE,IELAST,
     +                NSPIN,DETALF(1,IE,ISPIN),IREF,NATREF,NATPS,NEND,
     +                NSHELL,ISPIN,LCORE,NCORE,KSYMMAT,KESYM,DTB1)
          IF (IE.EQ.IELAST) THEN
            DO 10 L = 0,LMAXD + 1
              DNTOT = DNTOT + DIMAG(DNEALF(L,ISPIN,I1))
              DNMOM = DNMOM + DIMAG(DNEALF(L,ISPIN,I1))
   10       CONTINUE
          END IF
   20   CONTINUE
      END IF
      IF (IE.EQ.IELAST) THEN
        DNTOT = DNTOT + NDG(NP)*DIMAG(DNE(ISPIN,NP))
        DNMOM = DNMOM + NDG(NP)*DIMAG(DNE(ISPIN,NP))
      END IF
c
c---> print results of friedel's sum rule
c
      IF (IE.EQ.IELAST .AND. ISPIN.EQ.NSPIN .AND. NP.EQ.NREP) THEN
        WRITE (6,FMT=9000)
        DO 40 IS = 1,NSPIN
          WRITE (6,FMT=9010) TEXTS(IS),E,
     +      (DIMAG(NDG(NR)*DNE(IS,NR)),NR=1,NREP)
          DO 30 I1 = NATPS,NEND
            WRITE (6,FMT=9020) I1 - NATREF,
     +        DIMAG(DNEALF(LMAXD+1,IS,I1)),
     +        (L,DIMAG(DNEALF(L,IS,I1)),L=0,LMAXD)
   30     CONTINUE
   40   CONTINUE
        WRITE (6,FMT=9030) DNTOT
        IF (NSPIN.EQ.2) WRITE (6,FMT=9040) DNMOM
        WRITE (6,FMT=9050)
      END IF
c
c
      RETURN

 9000 FORMAT (/,1x,33 ('-'),' friedel"s sum rule ',33 ('-'))
 9010 FORMAT (3x,'for ',a13,' at the fermi energy ',
     +       'on the gf scale e = :',f10.8,1x,f10.8,/,5x,
     +       'symmetry decomposition : ',/,1x,10 (1x,f12.6))
 9020 FORMAT (3x,'impurity shell no.',i3,' non spherical contribution',
     +       ' of alpha matrix:',f12.6,/,3x,'spherical contribution',
     +       ' of alpha matrix:',/,3x,5 ('   l= ',i1,f12.6))
 9030 FORMAT (/,3x,' total change of charge in the cluster :',f15.8)
 9040 FORMAT (3x,' total change of moment in the cluster :',f15.8)
 9050 FORMAT (1x,82 ('-'))
      END
