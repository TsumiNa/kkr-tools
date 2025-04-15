      SUBROUTINE HYPREL(IKST,NPOT,NATREF,IRWS,IRNUMX,LMAXP1,NSPIN,ISPIN,
     +                  SIGNSP,GJM,EKL,IE,NLAST,DF,R,DRDI,RHO2NS,R2SUM,
     +                  EKLASYM)
      IMPLICIT NONE
c***********************************************************************
c     this subroutine calculates hyperfinefields and isomerieshifts
c     and nuclear-spin lattice relaxation time (t1t) for
c     valance states in case of relativistic treatment
c     charge-contactterm and breit-term is calculated
c     gtf = g*f   g = large component of sra   f = small component of sr
c     gtf/rr = integrand of breitformula
c***********************************************************************
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,LMAXD,LMX
      PARAMETER (irmd=1484,lmaxd=4,LMX=LMAXD+1)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DF
      REAL*8 SIGNSP
      INTEGER IE,IKST,IRNUMX,ISPIN,LMAXP1,NATREF,NLAST,NPOT,NSPIN
C     ..
C     .. Array Arguments ..
CASYM-------->
      COMPLEX*16 EKLASYM(LMX,NATYPD)
CASYM-------------------------------
      COMPLEX*16 EKL(*),GJM(LMX,*)
      REAL*8 DRDI(IRMD,*),R(IRMD,*),R2SUM(10,*),
     +                 RHO2NS(10,LMX,*)
      INTEGER IRWS(*)
C     ..
C     .. Arrays in Common ..
      COMPLEX*16 FZ(IRMD,LMX,NATYPD),PZ(IRMD,LMX,NATYPD),
     +               QZ(IRMD,LMX,NATYPD),SZ(IRMD,LMX,NATYPD),
     +               TMAT(LMX,NATYPD)
      REAL*8 RHOTRR(10,LMX,NATYPD),SUMTRR(10,NATYPD)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 FFZ,PPZ,QQZ,SSZ
      REAL*8 AL,C,FOUR,GAUSS,GTFRR,PG1,PG1RR,PG2,PG2RR,PI,RL,
     +                 RR,S1,S2,SUM1,SUM2,ZERO
      INTEGER ID,IH,IL,IP,IR,IRWS1,IS,IU,JR,JRP1,L,NATPS
C     ..
C     .. Local Arrays ..
      REAL*8 DRDI1(IRMD),GTF(IRMD,LMX,NATYPD),GTFR(IRMD),
     +                 GTFRT1(IRMD),GTFT1T(IRMD,LMX,NATYPD),
     +                 GTFT2T(IRMD),R1(IRMD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG,MOD,REAL
C     ..
C     .. Common blocks ..
      COMMON /CRDFUN/PZ,FZ,QZ,SZ,TMAT
      COMMON /T1T/RHOTRR,SUMTRR
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA PI,FOUR,ZERO/3.14159265359D0,4.D0,0.D0/
      DATA GAUSS,C/524.0D0,274.074D0/
C     ..
c
      NATPS = NATREF + 1
      IF (IE.LE.1) THEN
        IF (ISPIN.NE.2) THEN
          DO 30 IP = IKST,NPOT
            IRWS1 = IRWS(IP)
            DO 20 IL = 1,LMAXP1
              DO 10 IR = 1,IRWS1
                GTF(IR,IL,IP) = ZERO
                GTFT1T(IR,IL,IP) = ZERO
   10         CONTINUE
   20       CONTINUE
   30     CONTINUE
        END IF

      END IF
c
c-----------------------------------------------------------------------
c                     this step is only used as long we don't have relat
c                     green'sfunctions,in the other case we can ignore t
c                     and change 'do natps,npot' to 'do ikst,npot'
chyp-------------->1
chyp  IF (IKST.LE.1) THEN
      IF (IKST.LE.-100) THEN
        DO 60 IS = 1,NATREF
          IU = NSPIN*IS
          ID = IU - NSPIN + 1
          DO 50 IR = 2,IRNUMX
            S1 = ZERO
            S2 = ZERO
            RR = R(IR,IS)
            RR = FOUR*PI*RR*RR
            L = 1
            PPZ = PZ(IR,L,IS)
            QQZ = QZ(IR,L,IS)
            PG1 = DIMAG(PPZ*DF* (PPZ*GJM(L,IS)+EKL(L)*QQZ))
            PG1RR = PG1/RR
            PG2RR = SIGNSP*GAUSS*PG1RR
            S1 = S1 + PG1RR
            S2 = S2 + PG2RR
            RHO2NS(IR,L,ID) = RHO2NS(IR,L,ID) + PG1RR
            RHO2NS(IR,L,IU) = RHO2NS(IR,L,IU) + PG2RR
            RHOTRR(IR,L,IS) = PG2RR
c
            DO 40 L = 2,LMAXP1
              PPZ = PZ(IR,L,IS)
              QQZ = QZ(IR,L,IS)
              PG1 = DIMAG(PPZ*DF* (PPZ*GJM(L,IS)+EKL(L)*QQZ))
              PG1RR = PG1/RR
              PG2RR = SIGNSP*GAUSS*PG1RR
              RHO2NS(IR,L,ID) = RHO2NS(IR,L,ID) + PG1RR
              RHO2NS(IR,L,IU) = RHO2NS(IR,L,IU) + PG2RR
              RHOTRR(IR,L,IS) = PG2RR
   40       CONTINUE
            R2SUM(IR,ID) = R2SUM(IR,ID) + S1
            R2SUM(IR,IU) = R2SUM(IR,IU) + S2
            SUMTRR(IR,IS) = S2
   50     CONTINUE
   60   CONTINUE
      END IF
c-----------------------------------------------------------------------
chyp      DO 140 IP = NATPS,NPOT
      DO 140 IP = IKST,NPOT
        IU = NSPIN*IP
        ID = IU - NSPIN + 1
        IRWS1 = IRWS(IP)
        DO 130 L = 1,LMAXP1
CASYM
           EKL(L)=EKLASYM(L,IP)
CASYM
          DO 70 IR = 2,IRNUMX
            RR = R(IR,IP)
            RR = FOUR*PI*RR*RR
            PPZ = PZ(IR,L,IP)
            FFZ = FZ(IR,L,IP)
            QQZ = QZ(IR,L,IP)
            SSZ = SZ(IR,L,IP)
            PG1 = DIMAG(PPZ*DF* (PPZ*GJM(L,IP)+EKL(L)*QQZ))
            PG2 = DIMAG(FFZ*DF* (FFZ*GJM(L,IP)+EKL(L)*SSZ))
            S1 = (PG1+PG2)/RR
            RHO2NS(IR,L,ID) = RHO2NS(IR,L,ID) + S1
            IF (L.EQ.1) R2SUM(IR,ID) = R2SUM(IR,ID) + S1
   70     CONTINUE
          DO 80 IR = 2,IRWS1
            PPZ = PZ(IR,L,IP)
            FFZ = FZ(IR,L,IP)
            SSZ = SZ(IR,L,IP)
            GTFRR = SIGNSP*DIMAG(PPZ*DF* (FFZ*GJM(L,IP)+EKL(L)*SSZ))
            GTF(IR,L,IP) = GTF(IR,L,IP) + GTFRR
            GTFT2T(IR) = GTFRR
   80     CONTINUE
          IF (IE.EQ.NLAST) THEN
            DO 90 IR = 2,IRWS1
              GTFT1T(IR,L,IP) = GTFT1T(IR,L,IP) + GTFT2T(IR)
   90       CONTINUE
            IF (ISPIN.EQ.NSPIN) THEN
c
              RL = REAL(L-1)
              AL = (3.0D0*RL*RL-RL-1.0D0)/ (4.0D0*RL*RL-1)
c -------- al is a weight depending on l which comes from integral of sp
c -------- harmonics (a*a+b*b)*ylm*ylm l is real l ( l=0,1,2,3....)
c
              DO 100 IR = 2,IRWS1
                GTFR(IR) = GTF(IR,L,IP)
                GTFRT1(IR) = GTFT1T(IR,L,IP)
                DRDI1(IR) = DRDI(IR,IP)
                R1(IR) = R(IR,IP)
  100         CONTINUE
              DO 120 JR = 2,IRNUMX
                JRP1 = JR + 1
                SUM1 = ZERO
                SUM2 = ZERO
c------------------------------ simpson integration of breit-formula
                DO 110 IR = JRP1,IRWS1
                  RR = R1(IR)*R1(IR)
                  RR = DRDI1(IR)/RR
                  IH = 1 + MOD(JRP1+1+IR,2)
                  SUM1 = SUM1 + REAL(IH)*RR*GTFR(IR)
                  SUM2 = SUM2 + REAL(IH)*RR*GTFRT1(IR)
  110           CONTINUE
                RR = DRDI1(JR)/R1(JR)/R1(JR)
                SUM1 = SUM1 + SUM1 + RR*GTFR(JR)
                SUM2 = SUM2 + SUM2 + RR*GTFRT1(JR)
                SUM1 = -2.0D0*C*GAUSS/FOUR/PI*AL*SUM1/3.0D0
                SUM2 = -2.0D0*C*GAUSS/FOUR/PI*AL*SUM2/3.0D0
                RHO2NS(JR,L,IU) = SUM1
                RHOTRR(JR,L,IP) = SUM2
                IF (L.EQ.1) R2SUM(JR,IU) = R2SUM(JR,IU) + SUM1
                SUMTRR(JR,IP) = SUMTRR(JR,IP) + SUM2
  120         CONTINUE
            END IF

          END IF

  130   CONTINUE

  140 CONTINUE
c
      END
