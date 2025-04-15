      SUBROUTINE HYPISO(IPF,IKST,KVREL,KCOR,NATREF,IRNUMX,IRWS,NPOT,
     +                  ISPIN,NSPIN,NLAST,LMAXP1,CFG,R,DRDI,RHOC,IE,DF,
     +                  EKL,GJM,EKLASYM)
      IMPLICIT NONE
c
c this subroutine calculates hyperfinefields and isomerieshifts
c of host,impurity and next-neighbour atom depending on quantum
c number n,l
c the common block nucleus data (nucdat) is evaluated in subroutine
c corel.rhypf contains the hyperfinefields and isomerieshift of every
c corelevel and hypsum is summed up over all states and is so forth
c the spindensity of the up to tenth nearest point at the nucleus.
c remember (see also sub corel) hypsum = spindensity = radialspin
c desity /4/pi/r/r     (/ means dividet)
c rho2ns and r2sum correspond to rhypf and hypsum in valance case
c sumrrr (summ/r/r r=reduced)was introduced to recalculate only the
c impurity and neighbor fields and not host again. it's necessary!!!
c
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IRMD,LMAXD,LMX
      PARAMETER (irmd=1484,lmaxd=4,LMX=LMAXD+1)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DF

      INTEGER IE,IKST,IPF,IRNUMX,ISPIN,KCOR,KVREL,LMAXP1,NATREF,NLAST,
     +        NPOT,NSPIN
C     ..
C     .. Array Arguments ..
CASYM--------
      COMPLEX*16 EKLASYM(LMX,NATYPD)
CASYM--------
      COMPLEX*16 EKL(*),GJM(LMX,*)
      REAL*8 DRDI(IRMD,*),R(IRMD,*),RHOC(IRMD,*)
      INTEGER CFG(4,*),IRWS(*)
C     ..
C     .. Arrays in Common ..
      REAL*8 HYPSUM(10,NPOTD),RHORR(10,LMX,NATYPD),
     +                 RHYPF(150,NPOTD),SUMRR(10,NATYPD)
      INTEGER KFG(4,NATYPD),LMXC(NATYPD),NCMAX(NATYPD)
C     ..
C     .. Local Scalars ..
      REAL*8 DR23,DR24,DR34,FOUR,GAUSS,PI,R2,R3,R4,RHOCD,
     +                 RHOCU,RL2,RL3,RL4,RR,SIGNSP,ZERO
      INTEGER ID,IK,IL,INC1,INCM,IP,IR,IS,ISY,IU,IUD,L,LMCP1,LP1,N,NMAX
C     ..
C     .. Local Arrays ..
      REAL*8 R2SUM(10,NPOTD),RHO2NS(10,LMX,NPOTD),SUMRRR(10)
      CHARACTER*4 TEXT(5)
      CHARACTER*8 TEXTAT(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL HYPNRL,HYPREL
C     ..
C     .. Common blocks ..
      COMMON /NUCDAT/RHYPF,HYPSUM,KFG,LMXC,NCMAX
      COMMON /T1T/RHORR,SUMRR
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA TEXT/'=s  ','=p  ','=d  ','=f  ','=g  '/
      DATA TEXTAT/'    host','impurity'/
      DATA PI,FOUR,ZERO/3.14159265359D0,4.D0,0.D0/
      DATA GAUSS/524.0D0/
C     ..
c
      SIGNSP = -1.0D0
      IF (ISPIN.EQ.2) SIGNSP = +1.0D0
c
c---> paramagnetic case
c
      IF (NSPIN.EQ.1) SIGNSP = 0.0D0
      IF (IE.LE.1) THEN
        IF (ISPIN.NE.2) THEN
          DO 30 IP = IKST,NPOT
            IU = NSPIN*IP
            ID = IU - NSPIN + 1
            DO 20 IR = 1,IRNUMX
              R2SUM(IR,ID) = ZERO
              R2SUM(IR,IU) = ZERO
              SUMRR(IR,IP) = ZERO
              DO 10 L = 1,LMAXP1
                RHO2NS(IR,L,ID) = ZERO
                RHO2NS(IR,L,IU) = ZERO
                RHORR(IR,L,IP) = ZERO
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
        END IF

      END IF
c
c
      IF (KVREL.GE.1) THEN
        CALL HYPREL(IKST,NPOT,NATREF,IRWS,IRNUMX,LMAXP1,NSPIN,ISPIN,
     +              SIGNSP,GJM,EKL,IE,NLAST,DF,R,DRDI,RHO2NS,R2SUM,
     +              EKLASYM)

      ELSE

        CALL HYPNRL(IKST,NPOT,IRNUMX,LMAXP1,SIGNSP,GJM,EKL,DF,R,RHO2NS,
     +              R2SUM,NSPIN,EKLASYM)
      END IF
c - - - - -  lagrange extrapolation at origin
      DO 50 IP = IKST,NPOT
        IU = NSPIN*IP
        ID = IU - NSPIN + 1
        R2 = R(2,IP)
        R3 = R(3,IP)
        R4 = R(4,IP)
        DR23 = R2 - R3
        DR24 = R2 - R4
        DR34 = R3 - R4
        RL2 = (R3*R4)/ (DR23*DR24)
        RL3 = - (R2*R4)/ (DR23*DR34)
        RL4 = (R2*R3)/ (DR24*DR34)
c
        R2SUM(1,ID) = R2SUM(2,ID)*RL2 + R2SUM(3,ID)*RL3 +
     +                R2SUM(4,ID)*RL4
        R2SUM(1,IU) = R2SUM(2,IU)*RL2 + R2SUM(3,IU)*RL3 +
     +                R2SUM(4,IU)*RL4
        DO 40 L = 1,LMAXP1
          RHO2NS(1,L,ID) = RHO2NS(2,L,ID)*RL2 + RHO2NS(3,L,ID)*RL3 +
     +                     RHO2NS(4,L,ID)*RL4
          RHO2NS(1,L,IU) = RHO2NS(2,L,IU)*RL2 + RHO2NS(3,L,IU)*RL3 +
     +                     RHO2NS(4,L,IU)*RL4
   40   CONTINUE
   50 CONTINUE
c
      IF (IE.EQ.NLAST) THEN
        IF (ISPIN.EQ.NSPIN) THEN
c
          IF (KCOR.EQ.0) THEN
            DO 70 IP = IKST,NPOT
              IU = NSPIN*IP
              ID = IU - NSPIN + 1
              DO 60 IR = 2,IRNUMX
                RR = R(IR,IP)
                RR = FOUR*PI*RR*RR
                RHOCD = RHOC(IR,ID)
                RHOCU = RHOC(IR,IU)
                HYPSUM(IR,IU) = 0.5D0*GAUSS* (RHOCU-RHOCD)/RR
                HYPSUM(IR,ID) = 0.5D0* (RHOCU+RHOCD)/RR
   60         CONTINUE
   70       CONTINUE
          END IF
c
          DO 130 IS = 1,NSPIN
c
            IF (IS.EQ.1) WRITE (IPF,FMT=9000)
            IF (IS.EQ.2) WRITE (IPF,FMT=9010)
c
            DO 120 IP = IKST,NPOT
              ISY = 1
              IF (IP.GT.NATREF) ISY = 2
              IK = IP - (ISY-1)*NATREF
              IUD = NSPIN* (IP-1) + IS
              LMCP1 = LMXC(IP) + 1
              INCM = ZERO
              WRITE (IPF,FMT=9020) IK,TEXTAT(ISY)
c
c                                         core states
c
              WRITE (IPF,FMT=9030)
              WRITE (IPF,FMT=9040) (R(IR,IP),IR=1,IRNUMX)
              IF (KCOR.NE.0) THEN
                DO 90 LP1 = 1,LMCP1
                  NMAX = KFG(LP1,IP)
c
                  IF (NMAX.NE.0) THEN
                    DO 80 N = LP1,NMAX
                      INC1 = INCM + 1
                      INCM = INCM + IRNUMX
chyp                      WRITE (IPF,FMT=9050) N,TEXT(LP1),
chyp     +                  (RHYPF(IR,IUD),IR=INC1,INCM)
                      WRITE (IPF,FMT=9050) TEXT(LP1),
     +                  (RHYPF(IR,IUD),IR=INC1,INCM)
   80               CONTINUE
                  END IF

   90           CONTINUE

              END IF

              WRITE (IPF,FMT=9060) (HYPSUM(IR,IUD),IR=1,IRNUMX)
c
c                                          valence states
c
              WRITE (IPF,FMT=9070)
              DO 100 IL = 1,LMAXP1
chyp                N = CFG(IL,IP)
chyp                WRITE (IPF,FMT=9050) N,TEXT(IL),
chyp     +            (RHO2NS(IR,IL,IUD),IR=1,IRNUMX)
                WRITE (IPF,FMT=9050) TEXT(IL),
     +            (RHO2NS(IR,IL,IUD),IR=1,IRNUMX)
  100         CONTINUE
              WRITE (IPF,FMT=9060) (R2SUM(IR,IUD),IR=1,IRNUMX)
              DO 110 IR = 1,IRNUMX
                SUMRRR(IR) = R2SUM(IR,IUD) + HYPSUM(IR,IUD)
  110         CONTINUE
              WRITE (IPF,FMT=9060) (SUMRRR(IR),IR=1,IRNUMX)
  120       CONTINUE
  130     CONTINUE
c
          IF (NSPIN.EQ.1) WRITE (6,FMT=9080)
        END IF

      END IF




 9000 FORMAT (/,1x,37 ('$'),' isomer shifts ',37 ('$'))
 9010 FORMAT (/,1x,37 ('$'),' hyperfine fields ',36 ('$'))
 9020 FORMAT (1x,/,35x,i3,'th ',a8,'-cell')
 9030 FORMAT (1x,10 ('#'),' core states ',12 ('#'))
 9040 FORMAT (1x,'radius   :',1p,5d16.8, (/,11x,1p,5d16.8))
c 9050 FORMAT (1x,' n=',i1,' l',a4,5f16.6, (/,11x,5f16.6))
 9060 FORMAT (1x,'total    :',5f16.6, (/,11x,5f16.6))
 9070 FORMAT (1x,10 ('#'),' valance states ',10 ('#'))
 9080 FORMAT (1x,20 ('$'),' no spin polarisation - no hyperfinefields ',
     +       20 ('$'))
 9050 FORMAT (1x,' l',a4,5f16.6, (/,11x,5f16.6)) 
      END
