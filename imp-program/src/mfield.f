      SUBROUTINE MFIELD(THETAS,thesme,VISP,vspsme,IPAN,IRCUT,
     $     NTCELL,HFIELD,RFPI,ISPIN,
     +     KSHAPE,NATPS,NEND,NSPIN)
      Implicit None
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IRMD
      PARAMETER (irmd=1484)
      INTEGER NFUND,IRID
      PARAMETER (NFUND=289,irid=435)
      INTEGER NCELLD,IPAND
      PARAMETER (NCELLD=20,IPAND=80)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
C     ..
C     .. Local Scalars ..
      INTEGER I1,ICELL,IMT1,IP,J
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC REAL
C     ..
C     .. Scalar Arguments ..
      REAL*8 HFIELD,RFPI
      INTEGER ISPIN,KSHAPE,NATPS,NEND,NSPIN
C     ..
C     .. Array Arguments ..
      REAL*8 THETAS(IRID,NFUND,NCELLD),VISP(IRMD,NPOTD),
     $     thesme(irid,nfund,ncelld),vspsme(irmd,npotd)
      INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD),NTCELL(NATYPD)
C     ..
      DO 30 I1 = NATPS,NEND
        IP = NSPIN* (I1-1) + ISPIN
        DO 10 J = 1,IRCUT(1,I1)
          VISP(J,IP) = VISP(J,IP) + REAL(2*ISPIN-3)*HFIELD
          vspsme(j,ip) = vspsme(j,ip) + real(2*ispin-3)*hfield
   10   CONTINUE

        IF (KSHAPE.GE.1) THEN
          ICELL = NTCELL(I1)
          IMT1 = IRCUT(1,I1)
          DO 20 J = IMT1 + 1,IRCUT(IPAN(I1),I1)
            VISP(J,IP) = VISP(J,IP) + REAL(2*ISPIN-3)*HFIELD*
     +                   THETAS(J-IMT1,1,ICELL)/RFPI
            vspsme(j,ip) = vspsme(j,ip) + real(2*ispin-3)*hfield*
     +                   thesme(j-imt1,1,icell)/rfpi
   20     CONTINUE
        END IF

   30 CONTINUE

      RETURN
      END
