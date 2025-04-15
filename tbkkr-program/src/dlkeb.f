c 24.11.95 ***************************************************************
      SUBROUTINE DLKEB(GLLKE1,GLLKE2,GLLKE3,ALAT,NAEZ,NZ,CLS,EQINV,
     +                  NACLS,RR,EZOA,ATOM,BZKP,IE,KAOEZ,RCLS,GINP)
      implicit none
c ************************************************************************
c
c     Fourier transformation of the reference system Greens function
c     GLLKE is splitted into GLLKE1, GLLKE2, GLLKE3 due to linear
c     algorithm for matrix inversion
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
c      INTEGER NATOMD
c      PARAMETER (NATOMD=79)
      INTEGER LMAX
c      PARAMETER (LMAX=4)
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM,NDIM
      PARAMETER (ALM=LMAXSQ*NAEZD,NDIM=LMAXSQ*NPRINCD)
      DOUBLE COMPLEX CZERO,CI
      PARAMETER (CZERO= (0.0D0,0.0D0),CI= (0.0D0,1.0D0))
      DOUBLE COMPLEX CONE,CONEM
      PARAMETER (CONE= (1.0D0,0.0D0), CONEM= (-1.0D0,0.0D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IE,NAEZ
C     ..
C     .. Array Arguments ..
      INTEGER 
     +     ATOM(NACLSD,*),
     +     CLS(*),
     +     EZOA(NACLSD,*),
     +     EQINV(*),
     +     KAOEZ(*),
     +     NACLS(*)
c
      DOUBLE COMPLEX 
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +     GLLKE1(NDIM,NDIM,*),
     +     GLLKE2(NDIM,NDIM,*),
     +     GLLKE3(NDIM,NDIM,*)
c
      DOUBLE PRECISION 
     +     BZKP(*),
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,NZ,
     +     DI,DILM,DJ,DJLM,DJLM1,
     +     IC,IM,INV,
     +     J,JJ,JM,JLM,JLM1,
     +     LM1,LM2,M,NLAYER,N1,N2
      DOUBLE COMPLEX FAC
      LOGICAL OPT,TEST
c     ..
c     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE0(ALM,LMAXSQ)
      DOUBLE PRECISION KP(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,DLKE1,GLLCOPY,OPT,TEST,ZAXPY
C     ..
C     .. Save statement ..
      SAVE
C     ..
c ------------------------------------------------------------------------
c      write(6,*) '>>> DLKEB : Fourier-transforms the ',
c     +           'GF of reference system'

      NLAYER=NAEZ/NPRINCD

      DO 20 I=1,NAEZ

        FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around inv = eqinv(i)
c
        IF (I.NE.EQINV(I)) FAC = CONEM

        KP(1) = BZKP(1)
        KP(2) = BZKP(2)
        KP(3) = BZKP(3)
        IF (OPT('COMPLEX ')) THEN
          KP(4) = BZKP(4)
          KP(5) = BZKP(5)
          KP(6) = BZKP(6)
        END IF

c        write(6,*) 'fac :',fac

        IC  = CLS(KAOEZ(I))

        II   = (I-1)/NPRINCD + 1
        DI   = I - (II-1)*NPRINCD

        CALL DLKE1(GLLKE0,ALAT,NACLS,RR,
     +       EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +       RCLS(1,1,IC))

c        lm=3*lmaxsq-1
c        write(6,*) gllke0(lm+16,1)
c        write(6,*) gllke0(lm+17,16)

        DO 100 J = 1,NAEZ

          JJ   = (J-1)/NPRINCD + 1
          DJ   = J - (JJ-1)*NPRINCD

          IF (II.EQ.JJ) THEN 
c
c --->      main diagonal
c
            CALL GLLCOPY(GLLKE0,GLLKE2(1,1,II),J,DJ,DI)
          END IF

          IF (II.EQ.JJ+1) THEN 
c
c --->      upper diagonal
c
            CALL GLLCOPY(GLLKE0,GLLKE1(1,1,JJ),J,DJ,DI)
          END IF

          IF (II+1.EQ.JJ) THEN 
c
c --->      lower diagonal
c
            CALL GLLCOPY(GLLKE0,GLLKE3(1,1,II),J,DJ,DI)
          END IF

          IF (II.EQ.NLAYER .AND. JJ.EQ.1) THEN 
c
c --->      upper corner element
c
            CALL GLLCOPY(GLLKE0,GLLKE3(1,1,II),J,DJ,DI)
          END IF

          IF (JJ.EQ.NLAYER .AND. II.EQ.1) THEN 
c
c --->      lower corner element
c
            CALL GLLCOPY(GLLKE0,GLLKE1(1,1,JJ),J,DJ,DI)
          END IF
 100    END DO
c ------------------------------------------------------------------------

        IF (OPT('symG(k) ')) THEN
          write(6,*) 'SYMMETRIZATION in DLKEB NOT POSSIBLE.'
          STOP
        END IF                      ! (OPT('symG(k) '))

 20   CONTINUE                      ! I=1,NAEZ

      RETURN

      END                           ! DLKEB
