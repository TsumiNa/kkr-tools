      SUBROUTINE HYPNRL(IKST,NPOT,IRNUMX,LMAXP1,SIGNSP,GJM,EKL,DF,R,
     +                  RHO2NS,R2SUM,NSPIN,EKLASYM)
      IMPLICIT NONE
c***********************************************************************
c     this subroutine calculates hyperfinefields and isomerieshifts
c     and nuclear-spin lattice relaxation time (t1t) for
c     valance states in case of nonrelativistic treatment
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
      INTEGER IKST,IRNUMX,LMAXP1,NPOT,NSPIN
C     ..
C     .. Array Arguments ..
CASYM------->
      COMPLEX*16 EKLASYM(LMX,NATYPD)
C-----------
      COMPLEX*16 EKL(*),GJM(LMX,*)
      REAL*8 R(IRMD,*),R2SUM(10,*),RHO2NS(10,LMX,*)
C     ..
C     .. Arrays in Common ..
      COMPLEX*16 FZ(IRMD,LMX,NATYPD),PZ(IRMD,LMX,NATYPD),
     +               QZ(IRMD,LMX,NATYPD),SZ(IRMD,LMX,NATYPD),
     +               TMAT(LMX,NATYPD)
      REAL*8 RHOTRR(10,LMX,NATYPD),SUMTRR(10,NATYPD)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 PPZ,QQZ
      REAL*8 FOUR,GAUSS,PG1,PG1RR,PG2RR,PI,RR,S1,S2,ZERO
      INTEGER ID,IR,IS,IU,L
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG
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
      DATA GAUSS/524.0D0/
C     ..
c
      DO 30 IS = IKST,NPOT
        IU = NSPIN*IS
        ID = IU - NSPIN + 1
        DO 20 IR = 2,IRNUMX
          S1 = ZERO
          S2 = ZERO
          RR = R(IR,IS)
          RR = FOUR*PI*RR*RR
          L = 1
          EKL(L)=EKLASYM(L,IS)
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
          DO 10 L = 2,LMAXP1
            EKL(L)=EKLASYM(L,IS)  
            PPZ = PZ(IR,L,IS)
            QQZ = QZ(IR,L,IS)
            PG1 = DIMAG(PPZ*DF* (PPZ*GJM(L,IS)+EKL(L)*QQZ))
            PG1RR = PG1/RR
            PG2RR = SIGNSP*GAUSS*PG1RR
            RHO2NS(IR,L,ID) = RHO2NS(IR,L,ID) + PG1RR
            RHO2NS(IR,L,IU) = RHO2NS(IR,L,IU) + PG2RR
            RHOTRR(IR,L,IS) = PG2RR
   10     CONTINUE
          R2SUM(IR,ID) = R2SUM(IR,ID) + S1
          R2SUM(IR,IU) = R2SUM(IR,IU) + S2
          SUMTRR(IR,IS) = S2
   20   CONTINUE
   30 CONTINUE
c
      END
