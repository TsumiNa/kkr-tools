      Subroutine convol(imt1,irc1,icell,imaxsh,ilm,ifunm,lmpot,gsh,
     +     thetas,thesme,z,rfpi,r,vons,vspsmo,lmsp)
      Implicit None
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,LPOTD
      PARAMETER (irmd=1484,lpotd=8)
      INTEGER NFUND,IRID,NGSHD
      PARAMETER (NFUND=289,irid=435,NGSHD=54287)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      REAL*8 RFPI,Z
      INTEGER ICELL,IMAXSH,IMT1,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      REAL*8 GSH(*),R(*),THETAS(IRID,NFUND,*),VONS(IRMD,*),
     $     vspsmo(irmd),thesme(irid,nfund,*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),LMSP(NATYPD,*)
C     ..
C     .. Local Scalars ..
      REAL*8 ZZOR
      INTEGER I,IFUN,IR,IRH,LM,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      REAL*8 VSTORE(IRID,LMPOTD), vstsme(irid,lmpotd)
C     ..
      DO 20 LM = 1,LMPOT
        DO 10 IR = 1,IRC1 - IMT1
          VSTORE(IR,LM) = 0.0D0
          vstsme(ir,lm) = 0.0d0
   10   CONTINUE
   20 CONTINUE

      DO 30 IR = IMT1 + 1,IRC1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VONS(IR,1) - ZZOR
   30 CONTINUE

      DO 50 I = 1,IMAXSH
        LM1 = ILM(I,1)
        LM2 = ILM(I,2)
        LM3 = ILM(I,3)
        IF (LMSP(ICELL,LM3).GT.0) THEN
          IFUN = IFUNM(ICELL,LM3)
          DO 40 IR = IMT1 + 1,IRC1
            IRH = IR - IMT1
            VSTORE(IRH,LM1) = VSTORE(IRH,LM1) +
     +           GSH(I)*VONS(IR,LM2)*THETAS(IRH,IFUN,ICELL)
            vstsme(irh,lm1) = vstsme(irh,lm1) +
     +           gsh(i)*vons(ir,lm2)*thesme(irh,ifun,icell)
   40     CONTINUE
        END IF
   50 CONTINUE

      DO 60 IR = IMT1 + 1,IRC1
        IRH = IR - IMT1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VSTORE(IRH,1) + ZZOR
        vspsmo(ir) = (vstsme(irh,1) + zzor) /rfpi

   60 CONTINUE

c     Copy the part inside the MT sphere
      Do ir = 1,imt1
        vspsmo(ir) = vons(ir,1)/rfpi
      End Do

      DO 80 LM = 2,LMPOT
        DO 70 IR = IMT1 + 1,IRC1
          IRH = IR - IMT1
          VONS(IR,LM) = VSTORE(IRH,LM)
   70   CONTINUE
   80 CONTINUE
C     ..


      END
