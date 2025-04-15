      SUBROUTINE AMN(A,ALAT,B,IOPER,LMAX,NATPER,ND,NSHELL,RM,W,YR,YRA,
     +               WTYR,RIJ,IJEND)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculate the structure dependent matrices a and b which are used
c     for the determination of the nonspherical intercell-potential
c     since the program uses group theory only the potential of
c     one atom per shell has to be calculated ,i.e. shell-indices
c     can be used instead of atom-indices .
c     the intercell-potential is expanded into spherical harmonics .
c     the lm-term of the intercell-potential v of the shell i1
c     is given by
c
c      v(r,lm,i1) =  (-r)**l * {a(i1,i2,lm,l'm')*cmom(i2,l'm')
c                                                  +b(i1,i2,lm)*z(i2)}
c
c     summed over i2 (all shells) and l'm' .
c             (see notes by b.drittler)
c     attention : the first index of cmom (moment of the charge
c                 density - calculated in vintr2) and of z (nuclear
c                 charge of the atoms in the shell) is in the program
c                 different defined , there one has to use :
c                     cmom(natref+i2,l'm')  and  z(natref+i2)
c
c                 gaunt2 has to be called bevor to set up the common
c                 block assleg
c
c                               b.drittler   may 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
c
c---> attention : ncleb is an empirical factor - it has to be optimized
c
      INTEGER NATYPD,NTREFD,NTPERD
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD)
      INTEGER LMAXD,LMX,LPOTD
      PARAMETER (lmaxd=4,LMX=LMAXD+1,lpotd=8)
      INTEGER LMMAXD,L3D,LM3D
      PARAMETER (LMMAXD= (LPOTD+1)**2,L3D=2*LPOTD,LM3D= (L3D+1)**2)
      INTEGER N,LASSLD
      PARAMETER (N=4* (LMX-1),LASSLD=N)
      INTEGER LMPOTD,IJD
      PARAMETER (LMPOTD=(LPOTD+1)**2,IJD=434)
      INTEGER NCLEB
      PARAMETER (NCLEB=LM3D*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALAT
      INTEGER LMAX,NATPER
C     ..
C     .. Array Arguments ..
      REAL*8 A(NTPERD,NTPERD,LMMAXD,*),B(NTPERD,NTPERD,*),
     +       RIJ(IJD,3),RM(3,*),W(N),WTYR(IJD,LMPOTD),
     +       YR(N,0:LASSLD,0:LASSLD),YRA(IJD,LMPOTD)
      INTEGER IOPER(*),ND(48,3,3),NSHELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 CLECG,EPI,FACTOR,FPI,PI,R,R1,R2,R3,S
      INTEGER I,I1,I2,IATOM,IEND,IL,ILM,IMAX,IMIN,J,L,L1,L2,L3,L3MAX,
     +        LM1,LM2,LM3,LM4,LMMAX,LX,LY,M,M1,M1A,M1S,M2,M2A,M2S,M3,
     +        M3A,M3S,M4,NATOM,NI1R,NI2R
C     ..
C     .. Local Arrays ..
      REAL*8 C(0:LPOTD,-LPOTD:LPOTD,-LPOTD:LPOTD),
     +                 CLEB(NCLEB,2),DFAC(0:LPOTD,0:LPOTD),Y(LM3D)
      INTEGER ICLEB(NCLEB,4),LOFLM(LM3D)
C     ..
C     .. Statement Functions ..
      INTEGER MOFLM
      integer IJEND
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
C     .. External Subroutines ..
      EXTERNAL ROTCOEF,YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,SIGN
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
C     ..
C     .. Statement Function definitions ..
c
c---> calculate the m-value for given lm and l
c
      MOFLM(ILM,IL) = ILM - IL**2 - IL - 1
C     ..

c
c---> determine the l-value for given lm
c
      I = 1
      DO 20 L = 0,L3D
        DO 10 M = -L,L
          LOFLM(I) = L
          I = I + 1
   10   CONTINUE
   20 CONTINUE
c
      FPI = 4.0D0*PI
      EPI = 8.0D0*PI
      L3MAX = 2*LMAX
      LMMAX = (LMAX+1)**2
c
c--->calculate:                  (2*(l+l')-1)!!
c                dfac(l,l')= ----------------------
c                            (2*l+1)!! * (2*l'+1)!!
c
      DFAC(0,0) = 1.D0
      DO 40 LX = 1,LMAX
        DFAC(LX,0) = DFAC(LX-1,0)*REAL(2*LX-1)/REAL(2*LX+1)
        DFAC(0,LX) = DFAC(LX,0)
        DO 30 LY = 1,LX
          DFAC(LX,LY) = DFAC(LX,LY-1)*REAL(2* (LX+LY)-1)/REAL(2*LY+1)
          DFAC(LY,LX) = DFAC(LX,LY)
   30   CONTINUE
   40 CONTINUE
c
c---> set up of the gaunt coefficients with an index field
c     recognize that they are needed here only for l3=l1+l2
c
      I = 1
      DO 100 L1 = 0,LMAX
        DO 90 L2 = 0,LMAX
          L3 = L1 + L2
          DO 80 M1 = -L1,L1
            DO 70 M2 = -L2,L2
              DO 60 M3 = -L3,L3
                M1S = SIGN(1,M1)
                M2S = SIGN(1,M2)
                M3S = SIGN(1,M3)
c
                IF (M1S*M2S*M3S.GE.0) THEN
c
                  M1A = ABS(M1)
                  M2A = ABS(M2)
                  M3A = ABS(M3)
c
                  FACTOR = 0.0D0
c
                  IF (M1A+M2A.EQ.M3A) FACTOR = FACTOR +
     +                REAL(3*M3S+SIGN(1,-M3))/8.0D0
                  IF (M1A-M2A.EQ.M3A) FACTOR = FACTOR + REAL(M1S)/4.0D0
                  IF (M2A-M1A.EQ.M3A) FACTOR = FACTOR + REAL(M2S)/4.0D0
c
                  IF (FACTOR.NE.0.0D0) THEN
c
                    IF (M1S*M2S.NE.1 .OR. M2S*M3S.NE.1 .OR.
     +                  M1S*M3S.NE.1) FACTOR = -FACTOR
c
                    S = 0.0D0
                    DO 50 J = 1,N
                      S = S + W(J)*YR(J,L1,M1A)*YR(J,L2,M2A)*
     +                    YR(J,L3,M3A)
   50               CONTINUE
                    CLECG = S*FACTOR
                    IF (ABS(CLECG).GT.1.D-10) THEN
                      CLEB(I,1) = CLECG
                      ICLEB(I,1) = L1* (L1+1) + M1 + 1
                      ICLEB(I,2) = L2* (L2+1) + M2 + 1
                      ICLEB(I,3) = L3* (L3+1) + M3 + 1
                      I = I + 1
                    END IF

                  END IF

                END IF

   60         CONTINUE
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE
      IEND = I - 1
      IF (NCLEB.LT.IEND) THEN
        STOP 13

      ELSE
        WRITE (6,FMT='(i10)') IEND
c
c---> initialize a(i1,i2,lm1,lm2) and b(i1,i2,lm1)
c
        DO 140 I2 = 1,NATPER
          DO 130 I1 = 1,NATPER
            DO 120 LM1 = 1,LMMAX
              B(I1,I2,LM1) = 0.0D0
              DO 110 LM2 = 1,LMMAX
                A(I1,I2,LM1,LM2) = 0.0D0
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
  140   CONTINUE
c
c---> calculate a(i1,i2,lm1,lm2) and b(i1,i2,lm1)
c

        NI2R = 1
c
c---> loop over shells (second index)
c
        DO 200 I2 = 1,NATPER
          IMIN = NI2R
          NI2R = NI2R + NSHELL(I2)
          IMAX = NI2R - 1
c
c---> sum over all atoms of the shell i2
c
          DO 190 IATOM = IMIN,IMAX
c
            CALL ROTCOEF(C,IATOM,IOPER,LMAX,ND,YRA,WTYR,RIJ,IJEND)
c
            NI1R = 1
c
c---> loop over shells (first index)
c
            DO 180 I1 = 1,NATPER
              NATOM = NI1R
              NI1R = NI1R + NSHELL(I1)
c
c---> sum only over different natoms
c
              IF (NATOM.NE.IATOM) THEN
                R1 = RM(1,NATOM) - RM(1,IATOM)
                R2 = RM(2,NATOM) - RM(2,IATOM)
                R3 = RM(3,NATOM) - RM(3,IATOM)
c
                CALL YMY(R1,R2,R3,R,Y,L3MAX)
c
                R = R*ALAT
c
c---> calculate b(i1,i2,lm1)
c
                DO 150 LM1 = 1,LMMAX
                  L1 = LOFLM(LM1)
                  B(I1,I2,LM1) = B(I1,I2,LM1) -
     +                           EPI/REAL(2*L1+1)*Y(LM1)/ (R** (L1+1))
  150           CONTINUE
c
c---> calculate a(i1,i2,lm1,lm2)
c
                DO 170 I = 1,IEND
                  LM1 = ICLEB(I,1)
                  LM4 = ICLEB(I,2)
                  LM3 = ICLEB(I,3)
                  L1 = LOFLM(LM1)
                  L2 = LOFLM(LM4)
                  M4 = MOFLM(LM4,L2)
c
c---> this loop has to be calculated only for l1+l2=l3
c
                  DO 160 M2 = -L2,L2
                    LM2 = L2* (L2+1) + M2 + 1
                    A(I1,I2,LM1,LM2) = A(I1,I2,LM1,LM2) +
     +                                 EPI*FPI*DFAC(L1,L2)*Y(LM3)*
     +                                 CLEB(I,1)*C(L2,M4,M2)/
     +                                 (R** (L1+L2+1))
  160             CONTINUE
  170           CONTINUE
              END IF

  180       CONTINUE
  190     CONTINUE
  200   CONTINUE
      END IF

      END
