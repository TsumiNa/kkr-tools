      SUBROUTINE VINTREL(AMAT,BMAT,CMOM,IREF,LMAX,NATPER,NATREF,NSPIN,
     +     NSTART,NEND,V,Z,R,IRWS,IRCUT,IPAN,KSHAPE,
     +     CMINST,AVMAD,BVMAD,AMATI,BMATI,DRM,WG,YRG,irm,ntim)
      implicit none
c-----------------------------------------------------------------------
c     This subroutine is used for the case of lattice relaxations.
c
c     Calculate the intercell-potentials and add these to the poten-
c     tial v  (in the spin-polarized case for each spin-direction
c     the intercell-potential is the same . )
c     It uses the structure dependent matrices amat and bmat which
c     correspond to the new, shifted lattice as well as the amati and
c     bmati which correspond to the ideal crystall lattice, and are
c     calculated in subroutine amn .
c     In the case of lattice relaxations always shape corrections are
c     used, this is because of the big overlap of the spheres in the
c     MT or ASA case. Therefore the madelung potential is always taken
c     into account
c
c  The algorithm:1.Obtain the madelung potential of the ideal host.
c  ------------- 2.Subtract the intercell potential of the host
c                  for the cluster, obtain Vout.
c                3.Transform the potential of the outer region Vout
c                  to the new centers (look SHFTVOUT).
c                4.Add the transformed outer potential to the cluster
c                  potential expanded around the new centers.
c
c     the charge-moments are calculated in the subroutine vintra2 ,
c     therefore vintra2 has to be called first .
c
c     attention : the first index of cmom (moment of the charge
c                 density - calculated in vintr2) and of z (nuclear
c                 charge of the atoms in the shell) is in the program
c                 different defined , there one has to use :
c                     cmom(natref+i2,l'm')  and  z(natref+i2)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NTREFD,NTPERD,NATYPD,NATOMD
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD)
      INTEGER LMAXD,LPOTD
      PARAMETER (lmaxd=4,lpotd=8)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER LASSLD
      PARAMETER (LASSLD=4*LMAXD)
C     ..
C     .. Scalar Arguments ..
      INTEGER LMAX,NATPER,NATREF,NEND,NSPIN,NSTART,KSHAPE, irm,ntim
C     ..
C     .. Array Arguments ..
      REAL*8  AMAT(NTPERD,NTPERD,LMPOTD,*),
     +      AVMAD(NTREFD,NTREFD,LMPOTD,*),
     +      BMAT(NTPERD,NTPERD,*),BVMAD(NTREFD,NTREFD,*),
     +      CMINST(LMPOTD,*),CMOM(LMPOTD,*),R(IRM,*),V(IRM,LMPOTD,*),
     +      Z(*),DRM(3,NTPERD)
      REAL*8 AMATI(NTPERD,NTPERD,LMPOTD,*),
     +                 BMATI(NTPERD,NTPERD,*),WG(LASSLD),
     +                 YRG(LASSLD,0:LASSLD,0:LASSLD)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IREF(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      REAL*8  AC,PI,SUM
      INTEGER I,I1,I2,IATYP,IATYP2,IPOT,IRS1,ISPIN,L,LM,LM2,LMMAX,M,
     +          NREF
C     ..
C     .. Local Arrays ..
      REAL*8 ACH(LMPOTD,NTREFD,2),VOUT(LMPOTD),VOUT1(LMPOTD),
     +                 SN(3)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE ACH
C     ..
      PI = 4.D0*DATAN(1.D0)
      LMMAX = (LMAX+1)* (LMAX+1)
      DO 10 IATYP = NSTART,NEND
         I1 = IATYP - NATREF
         IF (KSHAPE.NE.0) THEN
           IRS1 = IRCUT(IPAN(IATYP),IATYP)
         ELSE
           IRS1 = IRWS(IATYP)
         END IF
         IF (KSHAPE.GT.0) THEN
            IF (I1.LT.1) THEN
              DO 20 L = 0,LMAX
              DO 30 M = -L,L
              LM = L*L + L + M + 1
                AC = 0.0D0
               IF (NATREF.EQ.1) THEN
c
c---> lm = 1 component disappears if there is only one host atom
c
                 DO 40 LM2 = 2,LMMAX
                   SUM = CMOM(LM2,1) + CMINST(LM2,1)
                   AC = AC + AVMAD(IATYP,1,LM,LM2)*SUM
   40            CONTINUE
                ELSE
                 DO 50 I2 = 1,NATREF
                   AC = AC + BVMAD(IATYP,I2,LM)*Z(I2)
                    DO 60 LM2 = 1,LMMAX
                     SUM = CMOM(LM2,I2) + CMINST(LM2,I2)
                     AC = AC + AVMAD(IATYP,I2,LM,LM2)*SUM
   60               CONTINUE
   50            CONTINUE
                END IF
                ACH(LM,IATYP,ntim) = AC
                VOUT(LM)=AC
   30          CONTINUE
   20          CONTINUE
             ELSE
C ******* end of ideal host madelung potential calculation *********
               DO 120 L = 0,LMAX
               DO 130 M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0D0
               DO 70 I2 = 1,NATPER
               NREF = IREF(I2)
               DO 80 LM2 = 1,LMMAX
               SUM = - CMOM(LM2,NREF) - CMINST(LM2,NREF)
               AC = AC + AMATI(I1,I2,LM,LM2)*SUM
   80          CONTINUE
               AC = AC + BMATI(I1,I2,LM)* (-Z(NREF))
   70          CONTINUE
               AC = AC + ACH(LM,NREF,ntim)
               VOUT(LM) = AC
  130          CONTINUE
  120          CONTINUE
C
c  the vout = vmad - vclust (host) is calculated,shift it to new pos.
C
               DO 150 I=1,3
               SN(I)=DRM(I,I1)
  150          CONTINUE
               CALL SHFTVOUT(VOUT,VOUT1,SN,LMAX,WG,YRG)
C
               DO 220 L = 0,LMAX
               DO 230 M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0D0
               DO 170 I2 = 1,NATPER
               IATYP2=NATREF+ I2
               DO 180 LM2 = 1,LMMAX
               SUM =  CMOM(LM2,IATYP2) + CMINST(LM2,IATYP2)
               AC = AC + AMAT(I1,I2,LM,LM2)*SUM
  180          CONTINUE
               AC = AC + BMAT(I1,I2,LM)*Z(IATYP2)
  170          CONTINUE
               VOUT(LM) = AC + VOUT1(LM)
  230          CONTINUE
  220          CONTINUE
              END IF
        

        ELSE IF (I1.GE.1) THEN
             DO 420 L=0,LMAX
             DO 430 M = -L,L
             LM = L*L + L + M + 1
              AC = 0.0D0
              DO 370 I2=1,NATPER
              NREF = IREF(I2)
              IATYP2 = NATREF + I2
c AMATI : first index is in new coordinates, second index
c         in old coordinates
c AMAT  : Both indeces in new coordinates
c 
              DO 360 LM2 = 1,LMMAX
              AC = AC - AMATI(I1,I2,LM,LM2)*CMOM(LM2,NREF)
360           CONTINUE
              AC = AC - BMATI(I1,I2,LM)*Z(NREF)
370           CONTINUE
              VOUT(LM) = AC 
430           CONTINUE
420           CONTINUE
               DO 450 I=1,3
               SN(I)= DRM(I,I1)
  450          CONTINUE
               CALL SHFTVOUT(VOUT,VOUT1,SN,LMAX,WG,YRG)
               DO 520 L = 0,LMAX
               DO 530 M = -L,L
               LM = L*L + L + M + 1
               AC = 0.0D0
               DO 570 I2 = 1,NATPER
               IATYP2=NATREF+ I2
               DO 580 LM2 = 1,LMMAX
               SUM =  CMOM(LM2,IATYP2)
               AC = AC + AMAT(I1,I2,LM,LM2)*SUM
  580          CONTINUE
               AC = AC + BMAT(I1,I2,LM)*Z(IATYP2)
  570          CONTINUE
               VOUT(LM) = AC + VOUT1(LM)
  530          CONTINUE
  520          CONTINUE

         END IF


        DO 320 L = 0,LMAX
        DO 330 M = -L,L
        LM = L*L + L + M + 1

        IF (LM.EQ.1) WRITE (6,FMT=9000) I1, (VOUT(LM)/SQRT(4.D0*PI))
c
c---> add to v the intercell-potential
c
        DO 240 ISPIN = 1,NSPIN
c
c---> determine the right potential number
c
        IPOT = NSPIN* (IATYP-1) + ISPIN
c
c---> in the case of l=0 : r(1)**l is not defined
c
         IF (L.EQ.0) V(1,1,IPOT) = V(1,1,IPOT) + VOUT(LM)
         DO 250 I = 2,IRS1
         V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IATYP))**L*VOUT(LM)
  250    CONTINUE
  240    CONTINUE
  330    CONTINUE
  320    CONTINUE

   10 CONTINUE

 9000 FORMAT (1x,'spherically averaged intercell-potential for shell',
     +       i2,' :',1p,d14.6)

      END
