      SUBROUTINE CSIMPK(CF,CFINT,IPAN,IRCUT,DRDI)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     this subroutine does an integration up to rcut of an
c     complex function cf with an extended 3-point-simpson :
c
c                             rcut
c                      cfint = { cf(r') dr'
c                              0
c
c     modified for functions with kinks - at each kink the
c     integration is restarted .
c
c     attention : input cf is destroyed !
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER IPAND
      PARAMETER (IPAND=80)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 CFINT
      INTEGER IPAN
C     ..
C     .. Array Arguments ..
      COMPLEX*16 CF(*)
      REAL*8 DRDI(*)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      REAL*8 A1,A2
      INTEGER I,IEN,IP,IST,N
C     ..
C     .. External Functions ..
      COMPLEX*16 ZSUM
      EXTERNAL ZSUM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD
C     ..
      A1 = 4.0D0/3.0D0
      A2 = 2.0D0/3.0D0
      CFINT = 0.0D0
c
      DO 20 IP = 1,IPAN
c
c---> loop over kinks
c
        IST = IRCUT(IP-1) + 1
        IEN = IRCUT(IP)
c
        DO 10 I = IST,IEN
          CF(I) = CF(I)*DRDI(I)
   10   CONTINUE
c
        IF (MOD(IEN-IST,2).EQ.0) THEN
          CFINT = CFINT + (CF(IST)-CF(IEN))/3.0D0
          IST = IST + 1
          N = (IEN-IST+1)/2

        ELSE
c---> four point lagrange integration for the first step
          CFINT = CFINT + (9.0D0*CF(IST)+19.0D0*CF(IST+1)-
     +            5.0D0*CF(IST+2)+CF(IST+3))/24.0D0 +
     +            (CF(IST+1)-CF(IEN))/3.0D0
          IST = IST + 2
          N = (IEN-IST+1)/2
        END IF
c
c---> calculate with an extended 3-point-simpson
c
        CFINT = CFINT + A1*ZSUM(N,CF(IST),2) + A2*ZSUM(N,CF(IST+1),2)
   20 CONTINUE
c
      END
