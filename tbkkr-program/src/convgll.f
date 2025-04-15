c ************************************************************************
      SUBROUTINE CONVGLL(NLAYER,DIM,REFPOT,KAOEZ,DELTALL,LF)
      implicit none
c ************************************************************************
c     converts the GLLKE1, GLLKE2, and GLLKE3 matrices that G - t**(-1)
c     is complex hermitian
c     p.zahn, may 96
c ------------------------------------------------------------------------
      include 'inc.fi'
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER ALM,NLM,NDIM
      PARAMETER (ALM  = NDIMGK*LMMAXD,
     +           NDIM = NPRINCD*LMMAXD,
     +           NLM  = NAEZD*LMMAXD)
c
c     .. arguments
      INTEGER NLAYER,DIM
      INTEGER KAOEZ(*),
     +        LF(*),
     +        REFPOT(*)
      DOUBLE COMPLEX 
     +     DELTALL(LMMAXD,LMMAXD,*)
C
c     .. locals
      INTEGER I,II,IRF,IILM1,LM1,LM2,
     +        J,JJ,JRF,JJLM2,
     +        N,NP,NNP
      DOUBLE COMPLEX C1,C2,FACI
      LOGICAL LSTART
      DOUBLE COMPLEX FAC(LMMAXD,LMMAXD)
      save
c
c     .. arrays in common
c
      DOUBLE COMPLEX 
     +     GLLKE(ALM,ALM),
     +     GLLKE1(NDIM,NDIM,NLAYERD),
     +     GLLKE2(NDIM,NDIM,NLAYERD),
     +     GLLKE3(NDIM,NDIM,NLAYERD)
      COMMON / GLLKE / GLLKE,GLLKE1,GLLKE2,GLLKE3,C1,C2

      DATA LSTART / .TRUE. /
c ------------------------------------------------------------------------
      IF (DIM. NE. NDIM) STOP 'CONVGLL'

      NP = NPRINCD
      IF (LSTART) THEN
c
c --->  init of fac(lmmaxd,lmmaxd)
c
        FACI = (0.0D0,1.0D0)
        DO 1 LM1 = 1,LMMAXD
          DO 2 LM2 = 1,LMMAXD
            FAC(LM1,LM2) = FACI**(LF(LM1)-LF(LM2))
 2        END DO
 1      END DO
        LSTART = .FALSE.
      END IF

c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (NLAYER .GT. 1) THEN
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DO 10 N = 1,NLAYER-1
          NNP = (N-1)*NP
          DO 20 II = 1,NP
            IILM1 = (II-1)*LMMAXD
            DO 100 JJ = 1,NP
              JJLM2 = (JJ-1)*LMMAXD
c
c --->        GLLKE2
c
              I   = II + NNP
              J   = JJ + NNP
              IRF = REFPOT(KAOEZ(I))
              JRF = REFPOT(KAOEZ(J))
              DO 40 LM1=1,LMMAXD
                DO 50 LM2=1,LMMAXD
                  GLLKE2(IILM1+LM1,JJLM2+LM2,N) =
     +                 GLLKE2(IILM1+LM1,JJLM2+LM2,N)
     +                 *FAC(LM1,LM2)
     +                 *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 50             END DO              ! LM2=1,LMMAXD
 40           END DO                ! LM1=1,LMMAXD
c     
c --->        GLLKE3
c
              I   = II + NNP
              J   = JJ + NNP + NP
              IRF = REFPOT(KAOEZ(I))
              JRF = REFPOT(KAOEZ(J))
              DO 70 LM1=1,LMMAXD
                DO 80 LM2=1,LMMAXD
                  GLLKE3(IILM1+LM1,JJLM2+LM2,N) =
     +                 GLLKE3(IILM1+LM1,JJLM2+LM2,N)
     +                 *FAC(LM1,LM2)
     +                 *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 80             END DO              ! LM2=1,LMMAXD
 70           END DO                ! LM1=1,LMMAXD
c     
c --->        GLLKE1
c
              I   = II + NNP + NP
              J   = JJ + NNP
              IRF = REFPOT(KAOEZ(I))
              JRF = REFPOT(KAOEZ(J))
              DO 110 LM1=1,LMMAXD
                DO 120 LM2=1,LMMAXD
                  GLLKE1(IILM1+LM1,JJLM2+LM2,N) =
     +                 GLLKE1(IILM1+LM1,JJLM2+LM2,N)
     +                 *FAC(LM1,LM2)
     +                 *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 120            END DO              ! LM2=1,LMMAXD
 110          END DO                ! LM1=1,LMMAXD

 100        END DO                  ! JJ=1,NP
 20       END DO                    ! II=1,NP
          
 10     END DO                      ! N = 1,NLAYER-1
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END IF                        ! (NLAYER .GT. 1)
c ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      N = NLAYER
      NNP = (N-1)*NP

      DO 200 II=1,NP
        IILM1 = (II-1)*LMMAXD

        DO 210 JJ=1,NP
          JJLM2 = (JJ-1)*LMMAXD
c
c --->    GLLKE2
c
          I   = II + NNP
          J   = JJ + NNP
          IRF = REFPOT(KAOEZ(I))
          JRF = REFPOT(KAOEZ(J))
          DO 220 LM1=1,LMMAXD
            DO 230 LM2=1,LMMAXD
              GLLKE2(IILM1+LM1,JJLM2+LM2,N) =
     +             GLLKE2(IILM1+LM1,JJLM2+LM2,N)
     +             *FAC(LM1,LM2)
     +             *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 230        END DO                  ! LM2=1,LMMAXD
 220      END DO                    ! LM1=1,LMMAXD
c     
c --->    GLLKE3
c
          I   = II + NNP
          J   = JJ
          IRF = REFPOT(KAOEZ(I))
          JRF = REFPOT(KAOEZ(J))
          DO 250 LM1=1,LMMAXD
            DO 260 LM2=1,LMMAXD
              GLLKE3(IILM1+LM1,JJLM2+LM2,N) =
     +             GLLKE3(IILM1+LM1,JJLM2+LM2,N)
     +             *FAC(LM1,LM2)
     +             *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 260        END DO                  ! LM2=1,LMMAXD
 250      END DO                    ! LM1=1,LMMAXD
c     
c --->    GLLKE1
c     
          I   = II
          J   = JJ + NNP
          IRF = REFPOT(KAOEZ(I))
          JRF = REFPOT(KAOEZ(J))
          DO 290 LM1=1,LMMAXD
            DO 300 LM2=1,LMMAXD
              GLLKE1(IILM1+LM1,JJLM2+LM2,N) =
     +             GLLKE1(IILM1+LM1,JJLM2+LM2,N)
     +             *FAC(LM1,LM2)
     +             *DELTALL(LM1,LM1,IRF)*DELTALL(LM2,LM2,JRF)
 300        END DO                  ! LM2=1,LMMAXD
 290      END DO                    ! LM1=1,LMMAXD
          
 210    END DO                      ! JJ=1,NP
 200  END DO                        ! II=1,NP

      RETURN
      END
