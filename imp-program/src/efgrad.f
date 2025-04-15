      SUBROUTINE EFGRAD(CMOM,LMAX,NSPIN,NSTART,NEND,RHO2NS,V,R,DRDI,
     +                  IRWS,NATREF)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     calculates the electric field gradients from a given non -
c     spherical  charge density . the lattice summation is done
c     implicitly by using the coulomb potential at the sphere
c     boundary - suggested by m. weinert (private com. 1984)
c     the different components of the electric field gradients
c     stored in the following way :
c
c              efg(..,1) : intra atomic contribution
c              efg(..,2) : extra atomic contribution
c              efg(..,3) : total contribution
c
c                               b.drittler   sep. 1988
c                               changed June 1996
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMKD,LPOTD
      PARAMETER (irmkd=1484,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,NATREF,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      REAL*8 CMOM(LMPOTD,*),DRDI(IRMKD,*),R(IRMKD,*),
     +     RHO2NS(IRMKD,LMPOTD,NATYPD,*),V(IRMKD,LMPOTD,*)
      INTEGER IRWS(*)
C     ..
C     .. Local Scalars ..
      REAL*8 ETA,FAC,FACETA,PI,QX,QXX,QY,QYY,QZ,QZZ,RNORM,RWS,
     +     TEMP1,TEMP2,VINT1
      INTEGER I,IATYP,IPOT,IRWS1,J,J1,J2,J3,K,LM,M,N
C     ..
C     .. Local Arrays ..
      REAL*8 EFG(3,3,3),EFGLM(-2:2,3),EVAL(3,3),EVEC(3,3),
     +     EVNO(3,3),V1(IRMKD)
      CHARACTER*15 TEXT(3)
C     ..
      double precision ::  PARACd
C     .. External Functions ..
      INTEGER IDAMAX,IDAMIN
      EXTERNAL IDAMAX,IDAMIN
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
C     ..
C     .. Save statement ..
      SAVE PI,TEXT
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
      DATA TEXT/'intra atomic : ','extra atomic : ','       total : '/
C     ..
c
      FAC = SQRT(5.0D0/ (4.0D0*PI))
      IF (LMAX.LT.2) THEN
        WRITE (6,FMT=9000)
        STOP

      ELSE
c
        WRITE (6,FMT=9020)
        WRITE (6,FMT=9010)
        WRITE (6,FMT=9020)
c
        DO 100 IATYP = NSTART,NEND
c
          WRITE (6,FMT=9040) IATYP - NATREF
          WRITE (6,FMT=9060)
c
          IRWS1 = IRWS(IATYP)
          RWS = R(IRWS1,IATYP)
c
c---> determine the right potential numbers
c
          IPOT = NSPIN* (IATYP-1) + 1

          DO 20 M = -2,2
            LM = 6 + M + 1
c
            V1(1) = 0.0D0
            DO 10 I = 2,IRWS1
              V1(I) = RHO2NS(I,LM,IATYP,1)* (R(I,IATYP)** (-3))
   10       CONTINUE
c
c---> integrate with simpson subroutine
c
            CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
            EFGLM(M,1) = 8.0D0*PI/5.0D0*VINT1
c
c---> use coulomb potential to determine extra atomic contribution
c
            EFGLM(M,2) = (V(IRWS1,LM,IPOT)-
     +                   8.0D0*PI/5.0D0*CMOM(LM,IATYP)/ (RWS**3))/
     +                   (RWS**2)
            EFGLM(M,3) = EFGLM(M,1) + EFGLM(M,2)
   20     CONTINUE
c
c---> convert the efgs into cartesian representation
c
          DO 30 I = 1,3
c---> xx contribution
            EFG(1,1,I) = (EFGLM(2,I)*SQRT(3.0D0)-EFGLM(0,I))*FAC
c---> yy contribution
            EFG(2,2,I) = - (EFGLM(2,I)*SQRT(3.0D0)+EFGLM(0,I))*FAC
c---> zz contribution
            EFG(3,3,I) = EFGLM(0,I)*FAC*2.0D0
c---> xy contribution
            EFG(1,2,I) = EFGLM(-2,I)*FAC*SQRT(3.0D0)
            EFG(2,1,I) = EFG(1,2,I)
c---> xz contribution
            EFG(1,3,I) = EFGLM(1,I)*FAC*SQRT(3.0D0)
            EFG(3,1,I) = EFG(1,3,I)
c---> yz contribution
            EFG(2,3,I) = EFGLM(-1,I)*FAC*SQRT(3.0D0)
            EFG(3,2,I) = EFG(2,3,I)
   30     CONTINUE
c
c---> calculate eigenvalues of the tensor
c
C          CALL DEVCSF(3,EFG(1,1,3),3,EVAL(1,3),EVEC,3)
c
c---> sort : qzz > qyy > qxx
c
          J1 = IDAMAX(3,EVAL(1,3),1)
          QZZ = EVAL(J1,3)
          J3 = IDAMIN(3,EVAL(1,3),1)
          QXX = EVAL(J3,3)
          J2 = 6 - J1 - J3
          QYY = EVAL(J2,3)
c
          ETA = 0.0D0
          IF (ABS(QZZ).GT.1.0D-7) ETA = (QXX-QYY)/QZZ
          FACETA = SQRT(1.0D0+ETA*ETA/3.0D0)
c
          DO 50 N = 1,3
            RNORM = SQRT(EVEC(1,N)**2+EVEC(2,N)**2+EVEC(3,N)**2)
            DO 40 M = 1,3
              EVNO(M,N) = EVEC(M,N)/RNORM
   40       CONTINUE
   50     CONTINUE
c
          DO 90 I = 1,2
            DO 80 M = 1,3
              TEMP1 = 0.0D0
              DO 70 K = 1,3
                TEMP2 = 0.0D0
                DO 60 J = 1,3
                  TEMP2 = TEMP2 + EFG(K,J,I)*EVNO(J,M)
   60           CONTINUE
                TEMP1 = TEMP1 + EVNO(K,M)*TEMP2
   70         CONTINUE
              EVAL(M,I) = TEMP1
   80       CONTINUE
            WRITE (6,FMT=9070) TEXT(I),EVAL(J1,I),EVAL(J2,I),EVAL(J3,I)
   90     CONTINUE
c
          WRITE (6,FMT=9070) TEXT(3),QZZ,QYY,QXX
          QZ = 6.748335D0*QZZ
          QY = 6.748335D0*QYY
          QX = 6.748335D0*QXX
          WRITE (6,FMT=9080) TEXT(I),QZ,QY,QX
c---> nuclear quadr. moment 0.83 for Cd !!!!!!!
           PARACd = 0.83d0
c--->  This depends on the nucleus!!!!!!!!!
          QZ = paracd*1.174825D5*ABS(QZZ)*FACETA
          WRITE (6,FMT=9090) QZ
c---> nuclear quadr. moment 0.15
c          QZ = 0.15D0*1.174825D5*ABS(QZZ)*FACETA
c          WRITE (6,FMT=9090) QZ
          WRITE (6,FMT=9050) ETA
          WRITE (6,FMT=9100) (EVEC(J,J1),J=1,3)
          WRITE (6,FMT=9110) (EVEC(J,J2),J=1,3)
          WRITE (6,FMT=9120) (EVEC(J,J3),J=1,3)
  100   CONTINUE
c
        WRITE (6,FMT=9020)
        WRITE (6,FMT=9030)
        WRITE (6,FMT=9020)
      END IF

c

 9000 FORMAT (13x,'error stop in subroutine efgrad :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=2 at least !')
 9010 FORMAT (1x,33 ('-'),' electric field gradients ',33 ('-'),/,34x,
     +       ' in units ry/(a(bohr)**2) ')
 9020 FORMAT (1x,'>')
 9030 FORMAT (1x,82 ('-'))
 9040 FORMAT (3x,i5,'th shell')
 9050 FORMAT (30x,' eta : ',f10.6)
 9060 FORMAT (1x,'>',/,25x,' qzz : ',6x,' qyy : ',6x,' qxx : ',/,1x,'>')
 9070 FORMAT (7x,a15,f10.6,3x,f10.6,3x,f10.6,3x,
     +       'in units ry/(a(bohr)**2)')
 9080 FORMAT (7x,a15,f10.6,3x,f10.6,3x,f10.6,3x,
     +       'in units of 10**24 cm**-3')
 9090 FORMAT (22x,f14.4,' khz ')
 9100 FORMAT (1x,/,9x,'eigenvector of qzz : ',f10.6,2 (',',f10.6))
 9110 FORMAT (9x,'eigenvector of qyy : ',f10.6,2 (',',f10.6))
 9120 FORMAT (9x,'eigenvector of qxx : ',f10.6,2 (',',f10.6),/,1x,
     +       82 ('.'))
      END
