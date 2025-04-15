c 13.10.95 ***************************************************************
      SUBROUTINE RHOTOTB(IPF,NATYP,NSPIN,RHO2NS,RHOC,Z,DRDI,IRWS,IRCUT,
     +                   LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,
     +                   CHRGNT,ITC,NSHELL)
      implicit none
c ************************************************************************
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
c                               b.drittler   nov. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER NATYPD,NSPIND
c      PARAMETER (NATYPD=1,NSPIND=2)
c      INTEGER IRMD,LPOTD
c      PARAMETER (IRMD=1484,LPOTD=8)
c      INTEGER NFUND,IRID
c      PARAMETER (NFUND=24,IRID=435)
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION CHRGNT
      INTEGER ITC,IPF,KSHAPE,LPOT,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION DRDI(IRMD,*),RHO2NS(IRMD,LMPOTD,NATYPD,*),
     +                 RHOC(IRMD,*),THETAS(IRID,NFUND,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRWS(*),LLMSP(NATYPD,*),NFU(*),
     +        NSHELL(0:NSHELD),NTCELL(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION DIFF,FACTOR,RFPI,SUM,TOTMOM
      INTEGER I,I1,IATYP,ICELL,IFUN,IPAN1,IPOTD,IPOTU,IRC1,IRS1,ISPIN,
     +        LM,LMPOT
C     ..
C     .. Local Arrays ..
      INTEGER LMXCD0,KFGD0(4),CLSD0,REFPOTD0,NTCELLD0,IRNSD0,IER,J,IR
      DOUBLE PRECISION MTFACD0,SCFACTOR(NATYPD),FACTLCUT1,FACTLCUT,
     &                 ZD0
      CHARACTER*80 UIO
      DOUBLE PRECISION C(NATYPD,NSPIND),RHO(IRMD)
c
c
      LOGICAL OPT
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
      RFPI = SQRT(16.0D0*DATAN(1.0D0))
      LMPOT = (LPOT+1)**2
c
c---> loop over atoms
c
      IF (OPT('SCALE-L ')) THEN
      DO I=1,NATYP
         CALL IoInput('ATOMINFO  ',UIO,I+1,7,IER)
                           READ (UNIT=UIO,FMT=*)    ZD0,
     +                        LMXCD0,
     +                       (KFGD0(J),J=1,4),
     +                        CLSD0,
     +                        REFPOTD0,
     +                        NTCELLD0,
     +                        MTFACD0,
     +                        IRNSD0,SCFACTOR(I)
      END DO
      END IF 

      DO 70 IATYP = 1,NATYP
c
c--->   determine the right potential numbers for rhoc
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
c
c Scaling of the valence charge density to correct for the l-cutoff
c
c                            ***   3.4.2000  ***
c
        IF (OPT('SCALE-L ')) THEN
           DO ISPIN=1,NSPIN  
              FACTLCUT = SCFACTOR(IATYP)/NSPIN
              FACTLCUT1 = FACTLCUT
              DO IR=2,IRS1
        IF (DABS(SCFACTOR(IATYP)).LT.0.5D0) 
     &         FACTLCUT1 = 1.d0+FACTLCUT/RFPI/RHO2NS(IR,1,IATYP,ISPIN)
        IF (KSHAPE.GT.0) 
     &              FACTLCUT1 = 1.d0+FACTLCUT/RHO2NS(IR,1,IATYP,ISPIN)

                    RHO2NS(IR,1,IATYP,ISPIN) = 
     &                   RHO2NS(IR,1,IATYP,ISPIN)*FACTLCUT1
                 END DO
              END DO
              WRITE(6,*) 'Charge Scaled by : ',FACTLCUT
           END IF
c     
c
c End of scaling            ***   3.4.2000  ***
c

        DO 10 I = 2,IRS1
c
c--->     convert core density
c
          SUM = (RHOC(I,IPOTD)+RHOC(I,IPOTU))*FACTOR/RFPI
          DIFF = (RHOC(I,IPOTU)-RHOC(I,IPOTD))/RFPI
c
c--->     add this to the lm=1 component of rho2ns
c
          RHO2NS(I,1,IATYP,1) = RHO2NS(I,1,IATYP,1) + SUM
          RHO2NS(I,1,IATYP,NSPIN) = RHO2NS(I,1,IATYP,NSPIN) + DIFF

   10   CONTINUE                    ! I = 2,IRS1
c
c--->   calculate  charge and moment of the cell
c
        DO 60 ISPIN = 1,NSPIN
c
          IF (KSHAPE.EQ.0) THEN
c
c--->       integrate over wigner seitz sphere - no shape correction
c
            CALL SIMP3(RHO2NS(1,1,IATYP,ISPIN),SUM,1,IRS1,DRDI(1,IATYP))
c
c--->       the result has to be multiplied by sqrt(4 pi)
c           (4 pi for integration over angle and 1/sqrt(4 pi) for
c           the spherical harmonic y(l=0))
c
            SUM = SUM*RFPI

          ELSE                      ! (KSHAPE.EQ.0)
c
c--->       convolute charge density with shape function to get the
c           charge in the exact cell - if kshape .gt. 0
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
c--->       integrate over circum scribed sphere
c
            CALL SIMPK(RHO,SUM,IPAN1,IRCUT(0,IATYP),DRDI(1,IATYP))

          END IF                    ! (KSHAPE.EQ.0)

          C(IATYP,ISPIN) = SUM

          IF (ISPIN.NE.1) THEN

            IF (KSHAPE.NE.0) THEN
              WRITE (IPF,FMT=9010) IATYP,SUM
            ELSE
              WRITE (IPF,FMT=9050) IATYP,SUM
            END IF

          ELSE                      ! (ISPIN.NE.1)

            IF (KSHAPE.NE.0) THEN
              WRITE (IPF,FMT=9000) IATYP,SUM
            ELSE
              WRITE (IPF,FMT=9040) IATYP,SUM
            END IF

          END IF                    ! (ISPIN.NE.1)

   60   CONTINUE                    ! ISPIN = 1,NSPIN

   70 CONTINUE                      ! IATYP = 1,NATYP


      CHRGNT = 0.0D0
      DO 80 I1 = 1,NATYP
        CHRGNT = CHRGNT + REAL(NSHELL(I1))*(C(I1,1) - Z(I1))
   80 CONTINUE
      WRITE (IPF,FMT=9020) ITC,CHRGNT

      IF (NSPIN.EQ.2) THEN
        TOTMOM = 0.0D0
        DO 90 I1 = 1,NATYP
          TOTMOM = TOTMOM + REAL(NSHELL(I1))*C(I1,NSPIN)
   90   CONTINUE
        WRITE (IPF,FMT=9030) ITC,TOTMOM
      END IF

      RETURN

 9000 FORMAT (I4,' charge in wigner seitz cell =',f10.6)
 9010 FORMAT (I4,' moment in wigner seitz cell =',f10.6)
 9020 FORMAT (I4,' ******   charge neutrality in unit cell = ',f12.6)
 9030 FORMAT (I4,' ******   total mag. moment in unit cell = ',f12.6)
 9040 FORMAT (I4,' charge in wigner seitz sphere =',f10.6)
 9050 FORMAT (I4,' moment in wigner seitz sphere =',f10.6)

      END
