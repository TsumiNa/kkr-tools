c ************************************************************************
      SUBROUTINE TAUTOGQDOS(GSQ,TINVLL,DSYMLL,NSHELL,RFCTOR,IGF,
     +                   TAUVBZ,NSYMAT,NSH1,NSH2,RATOM,NQDOSKP)
      implicit none
c ************************************************************************
c
c     GLL0 = GMATLL in subroutine GLL96
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.D0,0.D0))
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX TAUVBZ
      DOUBLE PRECISION RFCTOR
      INTEGER IGF,NSYMAT,NSHELL,NQDOSKP
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX 
     +     GSQ(LMAXSQ,LMAXSQ,NAEZD,*),
     +     TINVLL(LMAXSQ,LMAXSQ,*),
     +     DSYMLL(LMAXSQ,LMAXSQ,*)
      DOUBLE PRECISION RATOM(3,*)
      INTEGER NSH1(*),NSH2(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IND,IU,LM1,LM2,NSLL,NS,J,IK
      LOGICAL LDIA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX GLL(LMAXSQ,LMAXSQ),
     +               TPG(LMAXSQ,LMAXSQ),
     +               XC(LMAXSQ,LMAXSQ)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CINIT,ZCOPY,ZSCAL,ZGEMM,TEST
      INTRINSIC ABS,ZABS
c
C     .. Save statement ..
      SAVE
c ------------------------------------------------------------------------
C ..
      IF (NAEZD.LT.NSHELL) THEN 
      WRITE(6,*) ' Error stop in tautogqdos.f '
      STOP
      END IF
      DO IK = 1,NQDOSKP

         DO 70 NS = 1,NSHELL
            I = NSH1(NS)
            J = NSH2(NS)
c     
        LDIA = 
     +       (DABS(RATOM(1,NS)**2+
     +            RATOM(2,NS)**2+
     +            RATOM(3,NS)**2  ) .LT. 1.0D-6)
c
            CALL ZCOPY(LMAXSQ*LMAXSQ,GSQ(1,1,NS,IK),1,GLL,1)
            CALL ZSCAL(LMAXSQ*LMAXSQ,TAUVBZ,GLL,1)

c
c --->  XC = TINVLL(I) * GLL
c
        CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TINVLL(1,1,I),
     +       LMAXSQ,GLL,LMAXSQ,CZERO,XC,LMAXSQ)
        
        IF (LDIA) THEN
c     
c --->    GLL = -TINVLL - TINVLL * GLL* TINVLL
c     
          CALL ZCOPY(LMAXSQ**2,TINVLL(1,1,I),1,GLL,1)
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +         TINVLL(1,1,I),LMAXSQ,-CONE,GLL,LMAXSQ)

        ELSE                        ! (LDIA)
c     
c --->    GLL =  - TINVLL(I) * GLL* TINVLL(J)
c
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +         TINVLL(1,1,J),LMAXSQ,CZERO,GLL,LMAXSQ)
          
        END IF                      ! (LDIA)
c
c --->  GLL0 = GLL/RFCTOR
c
        DO 30 LM1 = 1,LMAXSQ
          DO 20 LM2 = 1,LMAXSQ
            GSQ(LM2,LM1,NS,IK) = GLL(LM2,LM1)/RFCTOR
 20       CONTINUE
 30     CONTINUE

 70   CONTINUE                      ! NS = 1,NSHELL
      END DO ! Loop in ik
      RETURN
 9000 format(2f18.10)

      END
