      SUBROUTINE RHOTOT(ITC,IPF,IREF,NATREF,NATPER,NSHELL,NSPIN,NSTART,
     +                  NEND,RHO2NS,RHOC,Z,DRDI,IRWS,IRCUT,LPOT,NFU,
     +                  LLMSP,THETAS,NTCELL,KSHAPE,IPAN,irm,iri,ntim)
      Implicit None
c-----------------------------------------------------------------------
c     add core and valence density expanded in spherical harmonics
c         ( convention see subroutine rholm )
c     in the paramagnetic case (nspin=1) the core valence charge times
c         r**2 is add to the valence charge density times r**2
c         then only rho2ns(irmd,lmxtsq,natypd,1) is used .
c     in the spin-polarized case (nspin=2) the spin-splitted core
c         charge density times r**2 is converted into core charge
c         density times r**2 and core spin density times r**2 .
c         then these parts are added to corresponding parts of
c         the valence densities times r**2 , that are rho2ns(...,1)
c         which contains the charge density  and rho2ns(...,2) which
c         contains in that case the spin density .
c             (see notes by b.drittler)
c
c     attention : the core density is spherically averaged and multi-
c                 plied by 4 pi. therefore the core density is only
c                 added to l=0 part .
c
c                               b.drittler   june 1987
c      changed for more atoms per unit shell         11.6.1996
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IRMD,IRMKD,LPOTD
      PARAMETER (irmd=1484,irmkd=1484,lpotd=8)
      INTEGER NFUND
      PARAMETER (NFUND=289)
      Integer irmaxd
      Parameter ( irmaxd=max0(irmd,irmkd) )
c
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IPF,ITC,KSHAPE,LPOT,NATPER,NATREF,NEND,NSPIN,NSTART,
     $     irm,iri,ntim
C     ..
C     .. Array Arguments ..
      REAL*8 DRDI(IRM,*),RHO2NS(IRM,LMPOTD,NATYPD,*),
     +                 RHOC(IRM,*),THETAS(IRI,NFUND,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IREF(*),IRWS(*),LLMSP(NATYPD,*),
     +        NFU(*),NSHELL(*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 CHRCHG,DIFF,FACTOR,RFPI,SUM
      INTEGER I,IATYP,ICELL,IFUN,IP,IPAN1,IPOTD,IPOTU,IR,IRC1,IRS1,
     +        ISPIN,IT,LM,LMPOT,NATYP,NEND1
C     ..
C     .. Local Arrays ..
      REAL*8 C(NATYPD,NSPIND,2),RHO(irmaxd)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE C
C     ..
      RFPI = SQRT(16.0D0*DATAN(1.0D0))
      LMPOT = (LPOT+1)**2
      NATYP = NATREF + NATPER
c
      IF (ITC.EQ.1) THEN
        NEND1 = NATYP

      ELSE

        NEND1 = NEND
      END IF
c
c---> loop over atoms
c
      DO 70 IATYP = NSTART,NEND1
c
c---> determine the right potential numbers for rhoc
c
        IF (NSPIN.EQ.2) THEN
          IPOTD = 2*IATYP - 1
          IPOTU = 2*IATYP
          FACTOR = 1.0D0

        ELSE

          IPOTD = IATYP
          IPOTU = IATYP
          FACTOR = 0.5D0
        END IF

        IF (KSHAPE.NE.0) THEN
          IPAN1 = IPAN(IATYP)
          IRS1 = IRCUT(1,IATYP)
          IRC1 = IRCUT(IPAN1,IATYP)

        ELSE

          IRS1 = IRWS(IATYP)
          IRC1 = IRS1
        END IF

        DO 10 I = 2,IRS1
c
c---> convert core density
c
          SUM = (RHOC(I,IPOTD)+RHOC(I,IPOTU))*FACTOR/RFPI
          DIFF = (RHOC(I,IPOTU)-RHOC(I,IPOTD))/RFPI
c
c---> add this to the lm=1 component of rho2ns
c
          RHO2NS(I,1,IATYP,1) = RHO2NS(I,1,IATYP,1) + SUM
          RHO2NS(I,1,IATYP,NSPIN) = RHO2NS(I,1,IATYP,NSPIN) + DIFF
   10   CONTINUE
c
c---> calculate  charge and moment of the cell
c
        DO 60 ISPIN = 1,NSPIN
c
          IF (KSHAPE.EQ.0) THEN
c
c---> integrate over wigner seitz sphere - no shape correction
c
            CALL SIMP3(RHO2NS(1,1,IATYP,ISPIN),SUM,1,IRS1,DRDI(1,IATYP))
c
c---> the result has to be multiplied by sqrt(4 pi)
c       (4 pi for integration over angle and 1/sqrt(4 pi) for
c       the spherical harmonic y(l=0))
c
            SUM = SUM*RFPI

          ELSE

c
c---> convolute charge density with shape function to get the
c      charge in the exact cell - if kshape .gt. 0
c
            ICELL = NTCELL(IATYP)

            DO 20 I = 1,IRS1
              RHO(I) = RHO2NS(I,1,IATYP,ISPIN)*RFPI
   20       CONTINUE

            DO 30 I = IRS1 + 1,IRC1
              RHO(I) = 0.0D0
   30       CONTINUE

            DO 50 IFUN = 1,NFU(ICELL)
              LM = LLMSP(ICELL,IFUN)
              IF (LM.LE.LMPOT) THEN
                DO 40 I = IRS1 + 1,IRC1
                  RHO(I) = RHO(I) + RHO2NS(I,LM,IATYP,ISPIN)*
     +                     THETAS(I-IRS1,IFUN,ICELL)

   40           CONTINUE
              END IF

   50       CONTINUE

c
c---> integrate over circum scribed sphere
c
            CALL SIMPK(RHO,SUM,IPAN1,IRCUT(0,IATYP),DRDI(1,IATYP))
          END IF


          C(IATYP,ISPIN,ntim) = SUM
          IF (ISPIN.NE.1) THEN

            IF (KSHAPE.NE.0) THEN
              WRITE (IPF,FMT=9010) SUM

            ELSE
              WRITE (IPF,FMT=9050) SUM
            END IF

          ELSE IF (KSHAPE.NE.0) THEN
            WRITE (IPF,FMT=9000) SUM

          ELSE
            WRITE (IPF,FMT=9040) SUM
          END IF

   60   CONTINUE

   70 CONTINUE


c
c---> print information of the cutted shells
c     remember: cutted shells contain host data
      DO 80 IATYP = NEND1 + 1,NATYP
        IP = IATYP - NATREF

        IF (NSPIN.NE.1) THEN

          IF (KSHAPE.NE.0) THEN
            WRITE (IPF,FMT=9000) C(IREF(IP),1,ntim)
            WRITE (IPF,FMT=9010) C(IREF(IP),NSPIN,ntim)

          ELSE
            WRITE (IPF,FMT=9040) C(IREF(IP),1,ntim)
            WRITE (IPF,FMT=9050) C(IREF(IP),NSPIN,ntim)
          END IF

        ELSE IF (KSHAPE.NE.0) THEN
          WRITE (IPF,FMT=9000) C(IREF(IP),1,ntim)

        ELSE
          WRITE (IPF,FMT=9040) C(IREF(IP),1,ntim)
        END IF

   80 CONTINUE

      CHRCHG = 0.0D0
      DO 90 IP = 1,NATPER
c this is new
        IR=IREF(IP)
c this is new
        IT = IP + NATREF
c       CHRCHG = CHRCHG + NSHELL(IP)* (C(IT,1)-Z(IT))
        CHRCHG = CHRCHG + NSHELL(IP)* ((C(IT,1,ntim)-
     +    Z(IT))-(C(IR,1,ntim)-Z(IR)))
   90 CONTINUE
      WRITE (IPF,FMT=9020) CHRCHG
      IF (NSPIN.EQ.2) THEN
        CHRCHG = 0.0D0
        DO 100 IP = 1,NATPER
          IT = IP + NATREF
          IR = IREF(IP)
          CHRCHG = CHRCHG + NSHELL(IP)* (C(IT,2,ntim)-
     $         C(IR,2,ntim))
  100   CONTINUE
        WRITE (IPF,FMT=9030) CHRCHG
      END IF



 9000 FORMAT (4x,' charge in wigner seitz cell =',f10.6)
 9010 FORMAT (4x,' moment in wigner seitz cell =',f10.6)
 9020 FORMAT ('  ******   charge neutrality in cluster = ',f10.6)
 9030 FORMAT ('  ******   change of moment in cluster  = ',f10.6)
 9040 FORMAT (4x,' charge in wigner seitz sphere =',f10.6)
 9050 FORMAT (4x,' moment in wigner seitz sphere =',f10.6)
      END
