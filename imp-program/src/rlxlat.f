      SUBROUTINE RLXLAT(GMAT,TMATLL1,E,LMAX,NATPER,NATREF,
     &                  NATYP,DRM,NDIMNP,NP,ITMAT,INS,IRMIN,IREF,
     &                  NDIM,NREP,KSYMM,KLATR,NOTREAD,TMATHOST)
      implicit none
C  ******************************************************************
C  * This subroutine transforms the Green function to the shifted
C  * positions in the case of lattice relaxations
C  * This subroutine is called for each representation.
C  * On output:
C  *                        
C  *                ~o      
C  *      GMAT  ==  G  
C  *                   
C  *
C  *                  ~o 
C  *      TMATLL  =   t  
C  *
C  *                                             PRB. 36,6372,(1987)         
C  ******************************************************************
C     .. Parameters ..
      INTEGER NTREFD,NTPERD,NATOMD,NATYPD
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD)
      INTEGER NSEC,NREPD
      PARAMETER (nsec=689,NREPD=4)
      INTEGER LMX,LPOTD,LMAXD
      PARAMETER (lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
c   
c ... Array arguments ...
c
      INTEGER NDIMNP,NREP,KLATR,IOTMAT
      INTEGER IRMIN(NATYPD),IREF(NTPERD),NDIM(NREPD)
      REAL*8 DRM(3,NTPERD)
      COMPLEX*16 GMAT(NSEC,NSEC),GMATSH(NSEC,NSEC)
C
C ... Scalar arguments...
C
      INTEGER  LMAX,NATPER,NATREF,NATYP,ITMAT,INS,KSYMM
      COMPLEX*16 E
C
C ... Local arrays...
C
      COMPLEX*16 JMR(NSEC,NSEC),JMAT(LMMAXD,LMMAXD,NATYPD)
      COMPLEX*16 AR(NSEC,NSEC)
      COMPLEX*16 TMATHOST(LMMAXD,LMMAXD,NATYPD)
      COMPLEX*16 TMATLL(LMMAXD,LMMAXD,NATYPD)
      COMPLEX*16 TMATLL1(LMMAXD,LMMAXD,NATYPD)
      REAL*8 X,Y,Z
c
c ... Scalar arguments..
c
      INTEGER I,I1,NP,I2,IC,IA,LM1,LM2,LMMAX,I0,NOTREAD
      INTEGER IB,INFO,IHOST,l1,l2,m1,m2
      COMPLEX*16 CZERO,CONE
C ... Logical variables
      LOGICAL LCALL
c
c ... External routines...
c
      EXTERNAL CINIT,JMTRX,ZCOPY,SYTN1
C
C ... Save Statement ...
C 
      SAVE
C
      IF(NSEC.LT.LMMAXD) STOP 'rlxlat'
      CZERO = CMPLX(0.D0,0.D0)
      CONE = CMPLX(1.D0,0.D0)
      LCALL = .false.
      LMMAX = (LMAX+1)*(LMAX+1)
      
      IF (NP.EQ.1) THEN
CO        IF (NOTREAD .NE. 7) THEN
        DO 11 I1 = 1,NATREF
         CALL TMREAD(TMATLL(1,1,I1),ITMAT)
 11     CONTINUE
CO
CO----> IN CASE OF RIGID SHIFT SAVE HOST T-MATRIX
CO
        IF (KLATR .EQ. 2) THEN
           
        DO I1=1,NATREF
          DO L1=1,LMMAX
             DO L2=1,LMMAX
                TMATHOST(L1,L2,I1)=TMATLL(L1,L2,I1)
             END DO
          END DO
        END DO
        DO I1=NATREF+1,NATYP
           IHOST=IREF(I1-NATREF)
          DO L1=1,LMMAX
             DO L2=1,LMMAX
                TMATHOST(L1,L2,I1)=TMATLL(L1,L2,IHOST)
             END DO
          END DO
        END DO        
        END IF
CO        END IF
CO        IF (NOTREAD .EQ. 6) RETURN
CO
CO
CO
c 
c In this subroutine tmatll is the host t-matrix but transformed
c in the new coordinate system!
CO---->IN CASE OF G-VOID METHOD: DON'T TRANSFORM T-MATRIX
         DO 20 I1 = NATREF+1,NATYP
           IA = I1 - NATREF
           IHOST = IREF(I1 - NATREF)
           X = DRM(1,IA)
           Y = DRM(2,IA)
           Z = DRM(3,IA)
           CALL JMTRX(X,Y,Z,E,LMAX,JMAT(1,1,I1),LCALL)
           CALL CINIT(NSEC**2,AR)
           CALL CINIT(NSEC**2,GMATSH)
             DO 150 I0 = 1,LMMAX
               DO 130 I2 = 1,LMMAX
                  DO 140 IC = 1,LMMAX
                     AR(I0,I2) = AR(I0,I2) +
     +                         JMAT(I0,IC,I1)*TMATLL(IC,I2,IHOST)
 140              CONTINUE
 130           CONTINUE
 150         CONTINUE
            DO 180 I0 = 1,LMMAX
               DO 160 I2 = 1,LMMAX
                  DO 170 IC = 1,LMMAX
                    GMATSH(I0,I2) = GMATSH(I0,I2) +
     +                                  AR(I0,IC)*JMAT(I2,IC,I1)
 170              CONTINUE
 160           CONTINUE
 180        CONTINUE
            DO 10 LM1=1,LMMAX
              DO 10 LM2=1,LMMAX
 10             TMATLL1(LM1,LM2,I1) = GMATSH(LM1,LM2)
 20      CONTINUE
       END IF
      IF (KSYMM.GT.0) THEN
      CALL SYTN1(JMR,INS,LMAX,NATREF,1,NATPER,JMAT,IRMIN,NP,NDIM,NREP,
     +           .FALSE.) 
      ELSE
      CALL SYNOTN1(JMR,INS,LMAX,NATREF,1,NATPER,JMAT,IRMIN,NP,NDIM,NREP,
     +           .FALSE.,NATPER) 
      END IF 
      CALL CINIT(NSEC**2,AR)
      CALL CINIT(NSEC**2,GMATSH)
c
        DO 50 I1 = 1,NDIMNP
           DO 30 I2 = 1,NDIMNP
              DO 40 IC = 1,NDIMNP
                 AR(I1,I2) = AR(I1,I2) +
     +                       JMR(I1,IC)*GMAT(IC,I2)
  40          CONTINUE
  30       CONTINUE
  50    CONTINUE
          DO 80 I1 = 1,NDIMNP
             DO 60 I2 = 1,NDIMNP
                DO 70 IC = 1,NDIMNP
                  GMATSH(I1,I2) = GMATSH(I1,I2) +
     +                                  AR(I1,IC)*JMR(I2,IC)
  70            CONTINUE
  60         CONTINUE
  80      CONTINUE
c
       DO 200 IB = 1,NDIMNP
         DO 210 I = 1,NDIMNP
           GMAT(I,IB)   = GMATSH(I,IB)
210      CONTINUE
200    CONTINUE
CO       IF (NOTREAD .EQ. 6 ) STOP 'NOTREAD'
      END
