c 18.10.95 ***************************************************************
      SUBROUTINE RHOOUTB(CDEN,DF,GMAT,EK,I1,IE,ISPIN,IRMIN,IRC1,LMAX,
     +                   LMMAX,PNS,QNS,IELAST,RHO2NS,R2NEF,THETAS,ICELL,
     +                   IFUNM,IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,
     +                   IEND)
      implicit none
c ************************************************************************
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
      include 'inc.fi'
c     INTEGER NATYPD
c     PARAMETER (NATYPD=1)
c     INTEGER IRMD,IRNSD,LMAXD,LPOTD
c     PARAMETER (IRMD=1484,IRNSD=508,LMAXD=4,LPOTD=8)
c     INTEGER NFUND,IRID
c     PARAMETER (NFUND=24,IRID=435)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND,IRLMD
      PARAMETER (IRMIND=IRMD-IRNSD,IRLMD= (IRNSD+1)*LMMAXD)
c      INTEGER NCLEB,IRLLD
c      PARAMETER (NCLEB=LMPOTD*LMMAXD,IRLLD=IRLMD*LMMAXD)
      INTEGER IRLLD
      PARAMETER (IRLLD=IRLMD*LMMAXD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER I1,ICELL,IE,IELAST,IEND,IMT1,IPAN1,IRC1,IRMIN,ISPIN,LMAX,
     +        LMMAX,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(*),GMAT(LMMAXD,*),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,*),
     &               QNS1(LMMAXD,LMMAXD)
      DOUBLE PRECISION CLEB(*),R2NEF(IRMD,LMPOTD,NATYPD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),THETAS(IRID,NFUND,*)
      INTEGER ICLEB(NCLEB,4),IFUNM(NATYPD,*),LMSP(NATYPD,*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CLT,CONE,CZERO
      DOUBLE PRECISION C0LL
      INTEGER I,IFUN,IR,IRLN,J,L1,LM1,LM2,LM3,M1
      LOGICAL OPT
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
C     .. Arrays in Common ..
      DOUBLE COMPLEX WR(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. Common blocks ..
      COMMON /ZELLER/WR
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
c---> initialize array for complex charge density
c
      DO 20 L1 = 0,LMAX
        DO 10 I = 1,IRMD
          CDEN(I,L1) = CZERO
   10   CONTINUE
   20 CONTINUE
      IF (OPT('QDOS    ')) THEN
c
c Initialize WR ( stored in common)
c
         DO LM1 = 1,LMMAXD
            do lm2= 1,lmmaxd
               DO I = irmind,IRMD
                  wr(lm1,lm2,i) = CZERO
               end do 
            end do
         end do         
      end if
c
c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
c                                      summed over lm3
c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
c                                               summed over lm3
      DO 50 IR = IRMIN + 1,IRC1
c
c        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
c     +             GMAT,LMMAXD,EK,QNS(1,1,IR,1),LMMAXD)
c        CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
c     +             QNS(1,1,IR,1),LMMAXD,CZERO,WR(1,1,IR),LMMAXD)
c
c This change was done for the QDOS
c
c        
        DO LM1 = 1,LMMAXD
            do lm2= 1,lmmaxd
              qns1(lm1,lm2) = qns(lm1,lm2,ir,1)
c       write(6,*) qns1(lm1,lm2)
c       STOP 'QNS1'
            end do 
        end do
c       write(6,*) qns1(1,1),pns(1,1,ir,1),gmat(1,1)
        CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
     +             GMAT,LMMAXD,EK,QNS1(1,1),LMMAXD)
c       write(6,*) qns1(1,1),pns(1,1,1,1),gmat(1,1)
c       stop ' zgemm'

        CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,1),LMMAXD,
     +             QNS1(1,1),LMMAXD,CZERO,WR(1,1,IR),LMMAXD)

        IF (NSRA.EQ.2) THEN

          DO LM1 = 1,LMMAXD
            do lm2= 1,lmmaxd
              qns1(lm1,lm2) = qns(lm1,lm2,ir,2)
            end do 
        end do
c
c This change was done for the q-dos 30.5.2000
c the problem is that qns is overwriten so on multiple calls
c for the same energy but different k-points the q-dos was not
c correct
c
          CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
     +               LMMAXD,GMAT,LMMAXD,EK,QNS1(1,1),LMMAXD)
          CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
     +               LMMAXD,QNS1(1,1),LMMAXD,CONE,WR(1,1,IR),LMMAXD)

c          CALL ZGEMM('N','N',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
c     +               LMMAXD,GMAT,LMMAXD,EK,QNS(1,1,IR,2),LMMAXD)
c          CALL ZGEMM('N','T',LMMAX,LMMAX,LMMAX,CONE,PNS(1,1,IR,2),
c     +               LMMAXD,QNS(1,1,IR,2),LMMAXD,CONE,WR(1,1,IR),LMMAXD)
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
      DO 110 L1 = 0,LMAX
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
        IF (IE.EQ.IELAST) THEN
          DO 80 IR = IRMIN + 1,IRC1
            R2NEF(IR,1,I1,ISPIN) = R2NEF(IR,1,I1,ISPIN) +
     +                             C0LL*DIMAG(CDEN(IR,L1))
            RHO2NS(IR,1,I1,ISPIN) = RHO2NS(IR,1,I1,ISPIN) +
     +                              C0LL*DIMAG(CDEN(IR,L1)*DF)
   80     CONTINUE

        ELSE
          DO 90 IR = IRMIN + 1,IRC1
            RHO2NS(IR,1,I1,ISPIN) = RHO2NS(IR,1,I1,ISPIN) +
     +                              C0LL*DIMAG(CDEN(IR,L1)*DF)
c               write(6,*) ir,1,i1,ispin,rho2ns(ir,1,i1,ispin)
c        stop ' rhooutb'
   90     CONTINUE
        END IF
c
c
c
        IF (IPAN1.GT.1) THEN
          DO 100 I = IMT1 + 1,IRC1
            CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1,ICELL)*C0LL
  100     CONTINUE
        END IF

  110 CONTINUE
c
      IF (IPAN1.GT.1) THEN
        DO 120 I = 1,IRC1
          CDENNS(I) = 0.0D0
  120   CONTINUE
      END IF

      DO 160 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        CLT = CLEB(J)
c
c---> calculate the non spherically symmetric contribution
c
        IF (IE.EQ.IELAST) THEN
          DO 130 IR = IRMIN + 1,IRC1
            RHO2NS(IR,LM3,I1,ISPIN) = RHO2NS(IR,LM3,I1,ISPIN) +
     +                                DIMAG(CLT*DF*WR(LM1,LM2,IR))
            R2NEF(IR,LM3,I1,ISPIN) = R2NEF(IR,LM3,I1,ISPIN) +
     +                               DIMAG(CLT*WR(LM1,LM2,IR))
  130     CONTINUE

        ELSE
          DO 140 IR = IRMIN + 1,IRC1
            RHO2NS(IR,LM3,I1,ISPIN) = RHO2NS(IR,LM3,I1,ISPIN) +
     +                                DIMAG(CLT*DF*WR(LM1,LM2,IR))
  140     CONTINUE
        END IF
c
        IF (IPAN1.GT.1 .AND. LMSP(ICELL,LM3).GT.0) THEN
          IFUN = IFUNM(ICELL,LM3)
          DO 150 I = IMT1 + 1,IRC1
            CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
     +                  THETAS(I-IMT1,IFUN,ICELL)
  150     CONTINUE

        END IF

  160 CONTINUE

      END
