      SUBROUTINE GRULE(N,X,W)
      IMPLICIT NONE
c
c***********************************************************************
c
c     determines the (n+1)/2 nonnegative points x(i) and
c     the corresponding weights w(i) of the n-point
c     gauss-legendre integration rule, normalized to the
c     interval [-1,1]. the x(i) appear in descending order.
c
c     this routine is from 'methods of numerical integration',
c     p.j. davis and p. rabinowitz, page 369.
c
c***********************************************************************
c

C     .. Scalar Arguments ..
      INTEGER N
C     ..
C     .. Array Arguments ..
      REAL*8 W(*),X(*)
C     ..
C     .. Local Scalars ..
      REAL*8 D1,D2PN,D3PN,D4PN,DEN,DP,DPN,E1,FX,H,P,PI,PK,
     +                 PKM1,PKP1,T,T1,U,V,X0
      INTEGER I,IT,K,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC COS,DATAN
C     ..
      PI = 4.D0*DATAN(1.D0)
      M = (N+1)/2
      E1 = N* (N+1)
      DO 30 I = 1,M
        T = (4*I-1)*PI/ (4*N+2)
        X0 = (1.D0- (1.D0-1.D0/N)/ (8.D0*N*N))*COS(T)
c
c--->    iterate on the value  (m.w. jan. 1982)
c
        DO 20 IT = 1,2
          PKM1 = 1.D0
          PK = X0
          DO 10 K = 2,N
            T1 = X0*PK
            PKP1 = T1 - PKM1 - (T1-PKM1)/K + T1
            PKM1 = PK
            PK = PKP1
   10     CONTINUE
          DEN = 1.D0 - X0*X0
          D1 = N* (PKM1-X0*PK)
          DPN = D1/DEN
          D2PN = (2.D0*X0*DPN-E1*PK)/DEN
          D3PN = (4.D0*X0*D2PN+ (2.D0-E1)*DPN)/DEN
          D4PN = (6.D0*X0*D3PN+ (6.D0-E1)*D2PN)/DEN
          U = PK/DPN
          V = D2PN/DPN
          H = -U* (1.D0+.5D0*U* (V+U* (V*V-U*D3PN/ (3.D0*DPN))))
          P = PK + H* (DPN+.5D0*H* (D2PN+H/3.D0* (D3PN+.25D0*H*D4PN)))
          DP = DPN + H* (D2PN+.5D0*H* (D3PN+H*D4PN/3.D0))
          H = H - P/DP
          X0 = X0 + H
   20   CONTINUE
        X(I) = X0
        FX = D1 - H*E1* (PK+.5D0*H* (DPN+H/3.D0* (D2PN+.25D0*H* (D3PN+
     +       .2D0*H*D4PN))))
        W(I) = 2.D0* (1.D0-X(I)*X(I))/ (FX*FX)
   30 CONTINUE
      IF (M+M.GT.N) X(M) = 0.D0
      END
