c 27.11.98 ***************************************************************
      SUBROUTINE DLKESP(GSP,INGSP,NUMGSP,ALAT,NAEZ,NZ,CLS,EQINV,
     +                  NACLS,RR,EZOA,ATOM,BZKP,IE,KAOEZ,RCLS,GINP)
      implicit none
c ************************************************************************
c
c     Fourier transformation of the reference system Greens function
c     GLLKE is stored in linear array GSP due to SPARSE MATRIX
c     algorithm for matrix inversion
c     for dimension of GSP,INGSP,NUMGSP see routine SP2
c 
c     p. zahn
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM,NDIM
      PARAMETER (ALM=LMAXSQ*NAEZD,NDIM=LMAXSQ*NPRINCD)
      DOUBLE COMPLEX CZERO,CI,CONE,CONEM
      PARAMETER (
     +     CZERO = ( 0.D0,0.D0),
     +     CI    = ( 0.D0,1.D0),
     +     CONE  = ( 1.D0,0.D0), 
     +     CONEM = -CONE  )
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER IE,NAEZ,NUMGSP,NZ
C     ..
C     .. Array Arguments ..
      INTEGER 
     +     ATOM(NACLSD,*),
     +     CLS(*),
     +     EZOA(NACLSD,*),
     +     EQINV(*),
     +     KAOEZ(*),
     +     INGSP(NLSPD,*),
     +     NACLS(*)
c
      DOUBLE COMPLEX 
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +     GSP(NDIM,NDIM,*)
c
      DOUBLE PRECISION 
     +     BZKP(*),
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*)
C     ..
C     .. Local Scalars ..
      INTEGER I,II,
     +     DI,DJ,
     +     J,JJ,
     &     N,
     &     IC 
      DOUBLE COMPLEX FAC
      LOGICAL LGSP,OPT,TEST
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
      DATA LGSP /.TRUE./
      
      
c ------------------------------------------------------------------------
      IF (TEST('flow    ')) write(6,*) '>>> DLKESP'

      IF (LGSP) THEN
c
c --->  determination of index array INGSP of nonzero blocks of
c       GLLKE which are stored in linear block array GSP
c
        NUMGSP = 0
        DO 10 II = 1,NLSPD
          DO 11 JJ = 1,NLSPD
            INGSP(II,JJ) = 0
 11       END DO
 10     END DO          
C
        DO 30 J = 1,NAEZ
          JJ   = (J-1)/NPRINCD + 1
          IF (TEST('flow    '))       
     +         write(6,*) J,NACLS(CLS(KAOEZ(J)))
          DO 40 N = 1,NACLS(CLS(KAOEZ(J)))
            I  = ATOM(N,J)
            IF (TEST('flow    ')) write(6,*) '  ',N,I
            IF (I.GT.0) THEN
c
c --->        avoid artifical couplings in systems with 
c             reduced dimensionality (SLAB, WIRE)
c
              II = (I-1)/NPRINCD + 1
              IF (INGSP(II,JJ).EQ.0) THEN
                NUMGSP = NUMGSP + 1
                IF (NUMGSP.GT.NAUXSPD) THEN
                  write(6,*) 'Increase the parameter NAUXSPD in inc.fi.'
                  STOP '1 - DLKESP'
                END IF
                INGSP(II,JJ) = NUMGSP
              END IF
            END IF
 40       END DO
 30     END DO
c
        IF (TEST('INGSP   ')) THEN
          write(6,*) 'INGSP(ii,jj) :'
          DO 50 II = 1,NLSPD
            write(6,9000) (INGSP(II,JJ),JJ=1,NLSPD)
 50       END DO
        END IF
        IF (TEST('ingsp   ')) THEN
          write(6,*) 'ingsp(ii,jj) :'
          write(6,9009) (II,II=1,NLSPD)
          DO 60 II = 1,NLSPD
            write(6,9010) II
            DO 61 JJ = 1,NLSPD
              IF (INGSP(II,JJ).GT.0) THEN 
                WRITE(6,9011) 
              ELSE 
                WRITE(6,9012) 
              END IF
 61         END DO
            write(6,*) 
 60       END DO
        END IF
c
        write(6,9030) NUMGSP
        LGSP = .FALSE.
      END IF

      DO 20 I = 1,NAEZ

        FAC = CONE
c
c -->     cluster around i is invers symmetric 
c         to cluster around inv = eqinv(i)
c
        IF (I.NE.EQINV(I)) FAC = CONEM

        KP(1)= BZKP(1)
        KP(2)= BZKP(2)
        KP(3)= BZKP(3)
        IF (OPT('COMPLEX ')) THEN
          KP(4) = BZKP(4)
          KP(5) = BZKP(5)
          KP(6) = BZKP(6)
        END IF

        IC  = CLS(KAOEZ(I))

        II   = (I-1)/NPRINCD + 1
        DI   = I - (II-1)*NPRINCD

        CALL DLKE1(GLLKE0,ALAT,NACLS,RR,
     +       EZOA(1,I),ATOM(1,I),KP,IE,IC,FAC,GINP(1,1,IC),
     +       RCLS(1,1,IC))

        DO 100 J = 1,NAEZ
          JJ   = (J-1)/NPRINCD + 1
          DJ   = J - (JJ-1)*NPRINCD
          IF (INGSP(JJ,II).NE.0) 
     +         CALL GLLCOPY(GLLKE0,GSP(1,1,INGSP(JJ,II)),J,DJ,DI)
 100    END DO
c ------------------------------------------------------------------------
        IF (OPT('symG(k) ')) THEN
          write(6,*) 'SYMMETRIZATION in DLKESP NOT POSSIBLE.'
          STOP
        END IF                      ! (OPT('symG(k) '))

 20   CONTINUE                      ! I=1,NAEZ

      RETURN
 9000 FORMAT(20I4)
 9009 FORMAT(4x,200I4)
 9010 FORMAT(I4,$)
 9011 FORMAT('   x',$)
 9012 FORMAT('    ',$)
 9030 FORMAT(' NUMGSP =',I6)
      END                           ! DLKEB
