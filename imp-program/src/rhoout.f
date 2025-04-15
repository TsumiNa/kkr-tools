      SUBROUTINE RHOOUT(CDEN,DF,GMAT,EK,I1,IE,ISPIN,IRMIN,IRC1,LMAX,
     +                  LMMAX,NSPIN,NATREF,PNS,QNS,RHO2NS,THETAS,ICELL,
     +                  IFUNM,IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,
     +                  IEND,KESYM,DTB1,KSYMMAT)
CASYM------------------------>
c-----------------------------------------------------------------------
c
c     calculates the charge density from r(irmin) to r(irc)
c      in case of a non spherical input potential .
c
c     fills the array cden for the complex density of states
c
c     attention : the gaunt coeffients are stored in index array
c                   (see subroutine gaunt)
c
c     the structured part of the greens-function (gmat) is symmetric in
c       its lm-indices , therefore only one half of the matrix is
c       calculated in the subroutine for the back-symmetrisation .
c       the gaunt coeffients are symmetric too (since the are calculated
c       using the real spherical harmonics) . that is why the lm2- and
c       the lm02- loops are only only going up to lm1 or lm01 and the
c       summands are multiplied by a factor of 2 in the case of lm1 .ne.
c       lm2 or lm01 .ne. lm02 .
c
c             (see notes by b.drittler)
c
c                               b.drittler   aug. 1988
c-----------------------------------------------------------------------
C     .. Parameters ..
      IMPLICIT NONE
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,IRNSD,LMAXD,LPOTD,LMX
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,lpotd=8,LMX=LMAXD+1)
      INTEGER NFUND,IRID
      PARAMETER (NFUND=289,irid=435)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND,IRLMD
      PARAMETER (IRMIND=IRMD-IRNSD,IRLMD= (IRNSD+1)*LMMAXD)
      INTEGER NCLEB,IRLLD
      PARAMETER (NCLEB=LMPOTD*LMMAXD,IRLLD=IRLMD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      COMPLEX*16 DF,EK
      INTEGER I1,ICELL,IE,IEND,IMT1,IPAN1,IRC1,IRMIN,ISPIN,LMAX,LMMAX,
     +        NATREF,NSPIN,NSRA
C     ..
C     .. Array Arguments ..
      COMPLEX*16 CDEN(IRMD,0:*),CDENNS(*),GMAT(LMMAXD,LMMAXD,*),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,*)
      COMPLEX*16 QHELP(LMMAXD,LMMAXD)
      REAL*8 CLEB(*),RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 THETAS(IRID,NFUND,*)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),LMSP(NATYPD,*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CLTDF,CONE,CZERO
      REAL*8 C0LL
      INTEGER I,IFUN,IR,IRLN,J,L1,LM1,LM2,LM3,LN,M1,N
C     ..
CASYM---------------------------------- 
      INTEGER KESYM,KSYMMAT
      COMPLEX*16 DTB1(LMMAXD,LMMAXD,*)
CASYM----------------------------------
C     .. Local Arrays ..
      COMPLEX*16 WR(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,DIMAG,DSQRT
C     ..
C     .. Save statement ..
      SAVE CZERO,CONE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
C     ..
c
      IRLN = IRLLD*NSRA
C     C0LL = 1/sqrt(4*pi)
      C0LL = 1.0d0/DSQRT(16.0D0*DATAN(1.0D0))
c


      IF (I1.GT.NATREF) THEN
        N = I1 - (NATREF+1)

      ELSE
        N = I1 - 1
      END IF

      LN = LMMAX*N
c
c---> initialize array for complex charge density
c
      DO 20 L1 = 0,LMAX
        DO 10 I = 1,IRMD
          CDEN(I,L1) = CZERO
   10   CONTINUE
   20 CONTINUE
CASYM---->
      IF (I1 .GT. NATREF) THEN
      IF (KSYMMAT .GT. 0) THEN
      IF (IE .GE. KESYM)  THEN
CTESTASYM         write (*,*) 'rhoout:irmin:',IRMIN,IRC1,' natom:',n
      DO IR=IRMIN + 1,IRC1
        QHELP(:,:)=QNS(:,:,IR,1)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,DTB1(1,1,N+1),LMMAXD,
     +             QHELP,LMMAXD,CZERO,QNS(1,1,IR,1),LMMAXD)        
      END DO
      END IF
      END IF
      END IF
CASYM---->
CASYM---->
      IF (NSRA .EQ. 2)    THEN
      IF (I1 .GT. NATREF) THEN
      IF (KSYMMAT .GT. 0) THEN
      IF (IE .GE. KESYM)  THEN
CTESTASYM         write (*,*) 'rhoout:irmin:',IRMIN,IRC1,' natom:',n
      DO IR=IRMIN + 1,IRC1
        QHELP(:,:)=QNS(:,:,IR,2)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,DTB1(1,1,N+1),LMMAXD,
     +             QHELP,LMMAXD,CZERO,QNS(1,1,IR,2),LMMAXD)        
      END DO
      END IF
      END IF
      END IF
      END IF
CASYM---->
C------------------------------------------------------------------
c
c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
c                                      summed over lm3
c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
c                                               summed over lm3
      DO 50 IR = IRMIN + 1,IRC1
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
     +             GMAT(1,1,N+1),LMMAXD,EK,QNS(1,1,IR,1),LMMAXD)
        CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
     +             QNS(1,1,IR,1),LMMAXD,CZERO,WR(1,1,IR),LMMAXD)
        IF (NSRA.EQ.2) THEN
          CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
     +               LMMAXD,GMAT(1,1,N+1),LMMAXD,EK,QNS(1,1,IR,2),
     +               LMMAXD)
          CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
     +               LMMAXD,QNS(1,1,IR,2),LMMAXD,CONE,WR(1,1,IR),LMMAXD)
        END IF

        DO 40 LM1 = 1,LMMAX
          DO 30 LM2 = 1,LM1 - 1
            WR(LM1,LM2,IR) = WR(LM1,LM2,IR) + WR(LM2,LM1,IR)
   30     CONTINUE
   40   CONTINUE
   50 CONTINUE
c
c---> first calculate only the spherically symmetric contribution
c
      DO 100 L1 = 0,LMAX
        DO 70 M1 = -L1,L1
          LM1 = L1* (L1+1) + M1 + 1
          DO 60 IR = IRMIN + 1,IRC1
c
c---> fill array for complex density of states
c
            CDEN(IR,L1) = CDEN(IR,L1) + WR(LM1,LM1,IR)
   60     CONTINUE
   70   CONTINUE
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
        DO 80 IR = IRMIN + 1,IRC1
          RHO2NS(IR,1,I1,ISPIN) = RHO2NS(IR,1,I1,ISPIN) +
     +                            C0LL*DIMAG(CDEN(IR,L1)*DF)
   80   CONTINUE
c
c
c
        IF (IPAN1.GT.1) THEN
          DO 90 I = IMT1 + 1,IRC1
            CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1,ICELL)*C0LL
   90     CONTINUE
        END IF

  100 CONTINUE
c
      IF (IPAN1.GT.1) THEN
        DO 110 I = 1,IRC1
          CDENNS(I) = 0.0D0
  110   CONTINUE
      END IF

      DO 140 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        CLTDF = DF*CLEB(J)
c
c---> calculate the non spherically symmetric contribution
c
        DO 120 IR = IRMIN + 1,IRC1
          RHO2NS(IR,LM3,I1,ISPIN) = RHO2NS(IR,LM3,I1,ISPIN) +
     +                              DIMAG(CLTDF*WR(LM1,LM2,IR))
  120   CONTINUE
c
        IF (IPAN1.GT.1 .AND. LMSP(ICELL,LM3).GT.0) THEN
c       IF (IPAN1.GT.1) THEN
          IFUN = IFUNM(ICELL,LM3)
          DO 130 I = IMT1 + 1,IRC1
            CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
     +                  THETAS(I-IMT1,IFUN,ICELL)
  130     CONTINUE

        END IF

  140 CONTINUE


      END
