      SUBROUTINE LLOYDA(DNEALF,I1,DNEINT,DF,E,ALPHA,IE,IELAST,NSPIN,
     +                  DETALF,IREF,NATREF,NATPS,NEND,NSHELL,ISPIN,
     +                  LCORE,NCORE,KSYMMAT,KESYM,DTB1)
c----------------------------------------------------------------------
c     calculates the change of the density of states resulting from the
c     alpha matrices for each atom. (similar to LLOYDG)
c     R. Zeller      Sept. 1994
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      COMPLEX*16 CZERO,EI
      PARAMETER (CZERO= (0.0D0,0.0D0),EI= (0.0D0,1.0D0))
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
C     ..
CASYM LLOYD
      INTEGER KSYMMAT,KESYM,LM1,M
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,*)
      REAL*8 GL
C     .. Scalar Arguments ..
      COMPLEX*16 DF,E
      REAL*8 DNEINT
      INTEGER I1,IE,IELAST,ISPIN,NATPS,NATREF,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      COMPLEX*16 ALPHA(0:LMAXD,NATYPD),DETALF(NATYPD),DNEALF(0:LMX)
      INTEGER IREF(NTPERD),LCORE(20,NPOTD),NCORE(NPOTD),NSHELL(NTPERD)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CRIMAG,EONE,ETWO
      REAL*8 FAC,RLL1
      INTEGER IA,IHOST,IPOTH,IPOTI,L,N
C     ..
C     .. Local Arrays ..
      COMPLEX*16 DN3NEW(0:LMX,NATYPD),DN3OLD(0:LMX,NATYPD),
     +               FEST3(0:LMX,NATYPD),FONE(0:LMX,NATYPD),
     +               FTWO(0:LMX,NATYPD)
      REAL*8 COR3(0:LMX,NATYPD)
      INTEGER LLCORE(0:LMAXD,NATYPD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,CDLOG,DATAN,DIMAG,MOD
C     ..
C     .. Save statement ..
      SAVE
C     ..
      FAC = 1.D0/ (8.D0*DATAN(1.D0))
      IA = I1 - NATREF
      IHOST = IREF(I1-NATREF)
c
c---  determine correction for Lloyd's formula due to core states
c
      IPOTI = NSPIN* (I1-1) + ISPIN
      IPOTH = NSPIN* (IHOST-1) + ISPIN
      DO 10 L = 0,LMAXD
        LLCORE(L,I1) = 0
   10 CONTINUE
      DO 20 N = 1,NCORE(IPOTI)
        L = LCORE(N,IPOTI)
        LLCORE(L,I1) = LLCORE(L,I1) + 1
   20 CONTINUE
      DO 30 N = 1,NCORE(IPOTH)
        L = LCORE(N,IPOTH)
        LLCORE(L,I1) = LLCORE(L,I1) - 1
   30 CONTINUE
c
c--> initialize energy points and values to zero
c
      IF (IE.EQ.1) THEN
        IF (I1.EQ.NATPS) THEN
          EONE = CZERO
          ETWO = CZERO
        END IF
        DO 40 L = 0,LMAXD + 1
          DN3OLD(L,I1) = CZERO
          FONE(L,I1) = CZERO
          FTWO(L,I1) = CZERO
   40   CONTINUE
      END IF
CASYM--------------------------------------------------
      DO 50 L = 0,LMAXD
         GL=0.0D0
         DO M=-L,L
            LM1 = L* (L+1) + M + 1
            IF ((I1-NATREF) .LT. 1) STOP 'CASYM-LLOYDA'
            GL=GL+REAL(DTB1(LM1,LM1,I1-NATREF))
         END DO
         GL=GL/(REAL(2*L+1))
         IF (KSYMMAT .EQ. 0) GL=1.0D0
         IF ((KSYMMAT.GT.0).AND.(IE.LT.KESYM)) GL=1.0D0 
        CRIMAG = EI*MOD(LLCORE(L,I1),2)/2.D0
CASYM---------------->-->
        DN3NEW(L,I1) =GL*(CDLOG(ALPHA(L,I1)/ALPHA(L,IHOST))*FAC 
     +                - CRIMAG)
CASYM------------------------->
COLD        DN3NEW(L,I1) = CDLOG(ALPHA(L,I1)/ALPHA(L,IHOST))*FAC - CRIMAG
   50 CONTINUE
CASYM--------------------------------------------------
COLD      DO 50 L = 0,LMAXD
COLD        CRIMAG = EI*MOD(LLCORE(L,I1),2)/2.D0
COLD        DN3NEW(L,I1) = CDLOG(ALPHA(L,I1)/ALPHA(L,IHOST))*FAC - CRIMAG
COLD   50 CONTINUE
      DN3NEW(LMAXD+1,I1) = CDLOG(DETALF(I1)/DETALF(IHOST))*FAC
c
c---> correct jumps caused by non spherical contribution to alpha matrix
c
      IF (IE.LE.2) THEN
        DO 60 L = 0,LMAXD + 1
          COR3(L,I1) = ANINT(DIMAG(DN3NEW(L,I1)-DN3OLD(L,I1)))
   60   CONTINUE

      ELSE
c
c---> two point lagrange interpolation
c
        DO 70 L = 0,LMAXD + 1
          FEST3(L,I1) = (FTWO(L,I1)* (E-EONE)-FONE(L,I1)* (E-ETWO))/
     +                  (ETWO-EONE)
          COR3(L,I1) = ANINT(DIMAG(DN3NEW(L,I1)-FEST3(L,I1)))
   70   CONTINUE
      END IF
      DO 80 L = 0,LMAXD + 1
        RLL1 = (L+L+1)*4.D0/NSPIN
        IF (L.EQ.LMAXD+1) RLL1 = 4.D0/NSPIN
        DN3OLD(L,I1) = DN3NEW(L,I1) - EI*COR3(L,I1)
        FONE(L,I1) = FTWO(L,I1)
        FTWO(L,I1) = DN3OLD(L,I1)
        DNEALF(L) = DN3OLD(L,I1)*RLL1*NSHELL(IA)
        DNEINT = DNEINT + DIMAG(DF*DNEALF(L))/ (FAC*4.D0)*NSPIN
   80 CONTINUE
c
c---> update energy points
c
      IF (I1.EQ.NEND) THEN
        EONE = ETWO
        ETWO = E
      END IF

      END
