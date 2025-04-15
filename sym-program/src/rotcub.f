      SUBROUTINE ROTCUB(ND)
      IMPLICIT NONE
c-----------------------------------------------------------------
c-----------------------------------------------------------------
c  program  to  generate  rotation  matrices  of  cubic  group  Oh
c-----------------------------------------------------------------
c-----------------------------------------------------------------
C     .. Array Arguments ..
      INTEGER ND(48,3,3)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,J,K,K1,KI,L,LK,N1,N2,N3
C     ..
C     .. Local Arrays ..
      INTEGER NDET(48),NF1(24),NF2(24),NR(48,3,3),NRDET(48),NRTR(48),
     +        NRTR2(48),NTR(48),NTR2(48),NTR2T(5),NTRT(5)
C     ..
C     .. Save statement ..
      SAVE NTRT,NTR2T
C     ..
C     .. Data statements ..
      DATA NTRT/3,0,1,-1,-1/,NTR2T/3,0,1,3,1/
C     ..
c-----------------------------------------------------------------------
c    the next statements produce the 48 rotation matrices 'nd' for the
c    three dimensional vector space
c-----------------------------------------------------------------------
      DO 30 I = 1,3
        DO 20 J = 1,3
          DO 10 IR = 1,48
            ND(IR,I,J) = 0
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
      ND(1,1,1) = 1
      ND(1,2,2) = 1
      ND(1,3,3) = 1
      ND(2,2,1) = 1
      ND(2,1,2) = 1
      ND(2,3,3) = 1
      ND(3,3,1) = 1
      ND(3,1,3) = 1
      ND(3,2,2) = 1
      ND(4,3,2) = 1
      ND(4,2,3) = 1
      ND(4,1,1) = 1
      ND(5,1,2) = 1
      ND(5,2,3) = 1
      ND(5,3,1) = 1
      ND(6,2,1) = 1
      ND(6,3,2) = 1
      ND(6,1,3) = 1
      DO 80 K = 1,6
        DO 60 L = 1,3
          LK = K + 6*L
          DO 50 I = 1,3
            DO 40 J = 1,3
              ND(LK,I,J) = ND(K,I,J)
   40       CONTINUE
   50     CONTINUE
   60   CONTINUE
        DO 70 I = 1,3
          ND(K+6,I,1) = -ND(K+6,I,1)
          ND(K+12,I,2) = -ND(K+12,I,2)
          ND(K+18,I,3) = -ND(K+18,I,3)
   70   CONTINUE
   80 CONTINUE
c---------------------------------------------------------------
c     24 matrices are produced above. the other matrices are now
c     produced by inversion symmetry.
c---------------------------------------------------------------
      DO 110 K = 1,24
        DO 100 I = 1,3
          DO 90 J = 1,3
            ND(K+24,I,J) = -ND(K,I,J)
   90     CONTINUE
  100   CONTINUE
  110 CONTINUE
c------------------------------------------------------------------
c     below det(nd),trace(nd), and sum of squared diagonal elements
c     of 'nd' are calculated
c------------------------------------------------------------------
      DO 120 K = 1,48
        N1 = ND(K,2,2)*ND(K,3,3) - ND(K,2,3)*ND(K,3,2)
        N2 = ND(K,2,1)*ND(K,3,3) - ND(K,2,3)*ND(K,3,1)
        N3 = ND(K,2,1)*ND(K,3,2) - ND(K,3,1)*ND(K,2,2)
        NDET(K) = N1*ND(K,1,1) - N2*ND(K,1,2) + N3*ND(K,1,3)
        NTR(K) = ND(K,1,1) + ND(K,2,2) + ND(K,3,3)
        NTR2(K) = ND(K,1,1)**2 + ND(K,2,2)**2 + ND(K,3,3)**2
  120 CONTINUE
      K1 = 0
c---------------------------------------------------------------------
c     the next statements order the rotation matrices according to the
c     10 classes of the group 'oh' as:
c     e  8c3  6c4  3c42  6c2  i  8s6  6s4  3sigmah  6sigmad
c     the information in ndet,ntr,ntr2 is used
c----------------------------------------------------------------------
      DO 130 K = 1,48
        IF (NDET(K).NE.-1) THEN
          K1 = K1 + 1
          NF1(K1) = K
        END IF

  130 CONTINUE
      K1 = 0
      DO 150 L = 1,5
        DO 140 K = 1,24
          KI = NF1(K)
          IF (NTR(KI).EQ.NTRT(L)) THEN
            IF (NTR2(KI).EQ.NTR2T(L)) THEN
              K1 = K1 + 1
              NF2(K1) = NF1(K)
            END IF

          END IF

  140   CONTINUE
  150 CONTINUE
      DO 180 I = 1,3
        DO 170 J = 1,3
          DO 160 K = 1,24
            KI = NF2(K)
            NR(K,I,J) = ND(KI,I,J)
            NRDET(K) = NDET(KI)
            NRTR(K) = NTR(KI)
            NRTR2(K) = NTR2(KI)
            NR(K+24,I,J) = -ND(KI,I,J)
            NRDET(K+24) = -NDET(KI)
            NRTR(K+24) = -NTR(KI)
            NRTR2(K+24) = NTR2(KI)
  160     CONTINUE
  170   CONTINUE
  180 CONTINUE
      DO 210 K = 1,48
        DO 200 I = 1,3
          DO 190 J = 1,3
            ND(K,I,J) = NR(K,I,J)
  190     CONTINUE
  200   CONTINUE
  210 CONTINUE

      END
