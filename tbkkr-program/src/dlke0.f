

c 04.10.95 ***************************************************************
      SUBROUTINE DLKE0(GLLKE,ALAT,NAEZ,CLS,EQINV,NACLS,
     +                 RR,EZOA,ATOM,BZKP,IE,KAOEZ,RCLS,GINP)
c ************************************************************************
      implicit none
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
      INTEGER ALM,CLM
      PARAMETER (ALM=LMAXSQ*NDIMGK,CLM=LMAXSQ*NACLSD)
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
     +     GLLKE(ALM,*)
c
      DOUBLE PRECISION 
     +     BZKP(*),
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IC,IM,INV,J,JM,M,N,JN
      DOUBLE COMPLEX FAC
      LOGICAL OPT
c     ..
c     .. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(ALM,LMAXSQ)
      DOUBLE PRECISION KP(6)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,DLKE1,OPT,ZAXPY
C     ..
C     .. Save statement ..
      SAVE
C     ..
c      write(6,*) '>>> DLKE0 : Fourier-transforms the ',
c     +           'GF of reference system'
c ------------------------------------------------------------------------

      CALL CINIT(ALM*ALM,GLLKE(1,1))

      DO 20 I=1,NAEZ

        FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around eqinv(i)
c
        IF ( I.NE.EQINV(I) ) FAC = CONEM

        KP(1) = BZKP(1)
        KP(2) = BZKP(2)
        KP(3) = BZKP(3)
        IF (OPT('COMPLEX ')) THEN
          KP(4) = BZKP(4)
          KP(5) = BZKP(5)
          KP(6) = BZKP(6)
        END IF

        IC  = CLS(KAOEZ(I))

c        write(6,*) 'fac :',fac
c        write(6,*) 'ic :',ic

        CALL DLKE1(GLLKE1,ALAT,NACLS,RR,
     +       EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +       RCLS(1,1,IC))   ! Changed on 22.03.2000
                               ! Correction for fourier trans.
      
        DO 140 M=1,LMAXSQ
          IM=(I-1)*LMAXSQ+M
          DO 150 JN=1,LMAXSQ*NAEZ
            GLLKE(JN,IM) = GLLKE(JN,IM)+ GLLKE1(JN,M)
 150      CONTINUE
 140    CONTINUE
        


 20   CONTINUE                      ! I=1,NAEZ

c ------------------------------------------------------------------------
      IF (OPT('symG(k) ')) THEN
c     
c -->   symmetrization
c     
        DO 90 I=1,NAEZ
            
          FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around inv = eqinv(i)
c
          IF (I.NE.EQINV(I)) FAC = CONEM
            
          KP(1) = -BZKP(1)
          KP(2) = -BZKP(2)
          KP(3) = -BZKP(3)
          IF (OPT('COMPLEX ')) THEN
            KP(4) = -BZKP(4)
            KP(5) = -BZKP(5)
            KP(6) = -BZKP(6)
          END IF
          
          IC  = CLS(KAOEZ(I))
          
          CALL DLKE1(GLLKE1,ALAT,NACLS,RR,
     +         EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +         RCLS(1,1,IC))

          DO 120 J=1,NAEZ
            DO 110 M=1,LMAXSQ
              IM=(I-1)*LMAXSQ+M
              DO 100 N=1,LMAXSQ
                JN=(J-1)*LMAXSQ+N
                GLLKE(IM,JN) = (GLLKE(IM,JN)+ GLLKE1(JN,M))/2.0D0
 100          CONTINUE
 110        CONTINUE
 120      CONTINUE
          
 90     CONTINUE                    ! I=1,NAEZ

      END IF                        ! (OPT('symG(k) '))
c ------------------------------------------------------------------------

      RETURN

      END                           ! DLKE0
