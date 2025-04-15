      SUBROUTINE EPOTINS(EPOTIN,NSPIN,NSTART,NEND,RHO2NS,VM2Z,R,DRDI,
     +                   INS,IRMIN,IRWS,LPOT,VINS,IRCUT,IPAN,Z)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     attention : redefinition of energy zero ---> now muffin tin zero
c
c     calculate the energy of the input potential which is spherically
c     averaged .
c     the energy for the representive atom i is given by
c
c                               rws
c       epotin(i) = - sqrt(4 pi) {  dr' vm2z(r',i)*rho2ns(r',1,i)
c                                0
c
c     in case of non spherical input potential one has to add
c
c                 rirt
c            {  -  {  dr' vins(r',lm,i)rho2nso(r',1,lm,1)   }
c                 rmin
c                                        (summed over lm)
c
c     remember : the non spherical part of the input potential is
c                different from zero only between r(irmin) and r(irt)
c
c             (see notes by b.drittler)
c
c     attention: vm2z is the spherically averaged input potential ,
c                vins contains the non spherical contribution of the
c                potential and rho2ns(...,1) is the  real charge density
c                times r**2. vins and rho2ns are expanded into spherical
c                harmonics. (see deck rholm or rhons)
c
c
c     attention : in case of shape corrections it is important that
c                 muffin tin zero is used as zero of the energy scale .
c                 the shift vbc between mt zero and the electro static
c                 zero is treated a external potential - see deck
c                 ecoulom !
c
c     remember :  in case of shape corrections  the contribution of
c                 the nuclear potential - 2*Z/r has to be explicitly
c                 taken into account between muffin tin sphere and
c                 circum scribed sphere .
c                 only within the muffin tin sphere this term is
c                 analytically cancelled wtih the contribution of
c                 the coulomb potential - see deck ecoulom
c
c
c                 modified for non spherical potential and shape correc-
c                  tions
c
c                               b.drittler   oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,IRNSD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lpotd=8)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      INTEGER INS,LPOT,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      REAL*8 DRDI(IRMD,*),EPOTIN(*),R(IRMD,*),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRMIN(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      REAL*8 PI,R2RHOD,R2RHOU,RFPI,TEMP,ZZOR
      INTEGER I,IATYP,IC,IPAN1,IPOTD,IPOTU,IRC1,IRMIN1,IRS1,L1,LM,M1
C     ..
C     .. Local Arrays ..
      REAL*8 ENS(0:LPOTD,NATYPD),ER(IRMD)
      INTEGER IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
      PI = 4.D0*DATAN(1.D0)
      RFPI = SQRT(4.0D0*PI)

      DO 90 IATYP = NSTART,NEND
        IPAN1 = IPAN(IATYP)
        IRC1 = IRCUT(IPAN1,IATYP)
        IF (IPAN1.GT.1) THEN
          IRS1 = IRCUT(1,IATYP)

        ELSE
          IRS1 = IRWS(IATYP)
        END IF
c
        IF (NSPIN.EQ.1) THEN
          IPOTU = IATYP
          IPOTD = IATYP

        ELSE
          IPOTU = 2*IATYP - 1
          IPOTD = 2*IATYP
        END IF

        DO 10 I = 1,IRS1
c
c---> calculate charge density times input potential
c
          R2RHOU = (RHO2NS(I,1,IATYP,1)-RHO2NS(I,1,IATYP,NSPIN))/2.0D0
          R2RHOD = (RHO2NS(I,1,IATYP,1)+RHO2NS(I,1,IATYP,NSPIN))/2.0D0
          ER(I) = - (R2RHOU*VM2Z(I,IPOTU)+R2RHOD*VM2Z(I,IPOTD))*RFPI
   10   CONTINUE
c
c--->  remember the form of vm2z between mt sphere and rirc
c
        IF (IPAN1.GT.1) THEN
          DO 20 I = IRS1 + 1,IRC1
            R2RHOU = (RHO2NS(I,1,IATYP,1)-RHO2NS(I,1,IATYP,NSPIN))/2.0D0
            R2RHOD = (RHO2NS(I,1,IATYP,1)+RHO2NS(I,1,IATYP,NSPIN))/2.0D0
            ZZOR = 2.0D0*Z(IATYP)/R(I,IATYP)
            ER(I) = - (R2RHOU* (VM2Z(I,IPOTU)-ZZOR)+
     +              R2RHOD* (VM2Z(I,IPOTD)-ZZOR))*RFPI
   20     CONTINUE
        END IF
c
c---> now integrate er to get epotin
c
        IF (IPAN1.GT.1) THEN
          CALL SIMPK(ER,TEMP,IPAN(IATYP),IRCUT(0,IATYP),DRDI(1,IATYP))

        ELSE

          CALL SIMP3(ER,TEMP,1,IRS1,DRDI(1,IATYP))
        END IF
c
        EPOTIN(IATYP) = TEMP
        ENS(0,IATYP) = TEMP
c
c---> add non spher. contribution in case of non spher. input potential
c
        DO 30 L1 = 1,LPOT
          ENS(L1,IATYP) = 0.0D0
   30   CONTINUE
c
        IF (INS.NE.0) THEN
          IRMIN1 = IRMIN(IATYP)
          IF (IRMIN1.LE.IRS1) THEN
            IRCUTM(0) = IRMIN1 - 1
            DO 40 IC = 1,IPAN1
              IRCUTM(IC) = IRCUT(IC,IATYP)
   40       CONTINUE
c
            DO 80 L1 = 1,LPOT
              DO 50 I = 1,IRMD
                ER(I) = 0.0D0
   50         CONTINUE
              DO 70 M1 = -L1,L1
                LM = L1* (L1+1) + M1 + 1
                DO 60 I = IRMIN1,IRC1
c
c---> calculate charge density times potential
c
                  R2RHOU = (RHO2NS(I,LM,IATYP,1)-
     +                     RHO2NS(I,LM,IATYP,NSPIN))/2.0D0
                  R2RHOD = (RHO2NS(I,LM,IATYP,1)+
     +                     RHO2NS(I,LM,IATYP,NSPIN))/2.0D0
                  ER(I) = ER(I) - R2RHOU*VINS(I,LM,IPOTU) -
     +                    R2RHOD*VINS(I,LM,IPOTD)
   60           CONTINUE
   70         CONTINUE
              CALL SIMPK(ER,TEMP,IPAN1,IRCUTM,DRDI(1,IATYP))
c
              EPOTIN(IATYP) = EPOTIN(IATYP) + TEMP
              ENS(L1,IATYP) = TEMP
   80       CONTINUE

          END IF

        END IF

   90 CONTINUE

      IF (1.GT.2) THEN
        DO 100 IATYP = NSTART,NEND
          WRITE (6,FMT='(4(i3,1x,f15.8))') (L1,ENS(L1,IATYP),L1=0,LPOT),
     +      IATYP,EPOTIN(IATYP)
  100   CONTINUE
      END IF


      END
