C **********************************************************************
      SUBROUTINE BESHAN1(HL,JL,NL,Z,LMAX)
      implicit none
C **********************************************************************
C  CALCULATES BESSEL, HANKEL AND NEUMANN FUNCTIONS OF INTEGER ORDER
C  FOR THE ORDERS 0 .LE. L .LE. LMAX.
C  THE TAYLOR EXPANSIONS OF JL AND NL ARE USED.
C  FROM M.ABRAMOWITSCH, I.A. STEGUN, 'HANDBOOK OF MATH. FUNCTIONS'   
c ------------------------------------------------------------------------
C
C     .. PARAMETERS ..
      DOUBLE COMPLEX CI
      PARAMETER (CI= (0.0D0,1.0D0))
C     ..
C     .. SCALAR ARGUMENTS ..
      DOUBLE COMPLEX Z
      INTEGER LMAX
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE COMPLEX HL(0:LMAX),JL(0:LMAX),NL(0:LMAX)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE COMPLEX TERMJ,TERMN,TERMN1,Z2,ZJ,ZN
      DOUBLE PRECISION L0,LN,LLN,RL,RN,RNM,PI
      INTEGER L,M,N
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC CDABS,CDEXP
C     ..
      PI = 4.0D0 *DATAN(1.0D0)
      ZJ = 1.D0
      ZN = -2.0D0*(CDLOG(Z/2.0D0) +0.5772156649D0)
      Z2 = Z*Z
      DO 20 L = 0,LMAX
         RL = L + L
         L0 = L
         TERMJ = ZJ
         TERMN1= ZN
         JL(L) = ZJ
         IF (L.EQ.0) THEN 
            NL(L) = ZJ*ZN
            TERMN = 1.0D0 /ZJ
         ELSE 
            NL(L) = 1.0D0 /ZJ /L0 + ZJ*ZN
            TERMN = 1.0D0 /ZJ /L0
         END IF   
         N=0
         DO 10 WHILE ((CDABS(TERMJ)+CDABS(TERMJ*TERMN1)).GE.1.0D-20)
            N  = N + 1
            RN = N + N
            LN = N + L
            IF (N.LT.L) THEN
               LLN = L - N
               TERMN = 0.5D0 *TERMN /RN /LLN *Z2 
               NL(L) = NL(L) + TERMN
            END IF   
C     *********************************************************
C     
C     TERMJ= (-Z*Z/4)**N L! /N! /(N+L)!
C     
C                    Z              1         1
C     TERMN1 = -2(LN(-)-EXP(-1)) + (- + ... + -) +
C                    2              1         L
C     
C                1        1     1        1
C             +( - + .. + - )+( - + .. + - )
C               L+1      L+N    1        N
C     
C              WRITE(6,9000) TERMJ
C 9000         FORMAT(2F35.20)
C     *********************************************************
            TERMJ  = - 0.5D0 *TERMJ /RN /LN *Z2 
            TERMN1 = TERMN1 + 2.0D0 /RN + 1.0D0 / LN
            JL(L) = JL(L) + TERMJ
            NL(L) = NL(L) + TERMN1 * TERMJ
 10      CONTINUE                       ! WHILE (ABS() > ..)
         NL(L) = -NL(L)/PI
         HL(L) = JL(L) + NL(L)*CI
         
         ZJ = ZJ *Z /(RL+2.D0)
         ZN = ZN +2.0D0 /(RL+2.0D0)
C        WRITE(6,9001) L,ZJ
C 9001   FORMAT(I10,2F35.20)
 20   CONTINUE                          !  L = 0,LMAX
      
      RETURN

      END
