      SUBROUTINE VINTERS(AMAT,BMAT,CMOM,IREF,LMAX,NATPER,NATREF,NSPIN,
     +     NSTART,NEND,V,Z,R,IRWS,IRCUT,IPAN,KSHAPE,
     +     CMINST,AVMAD,BVMAD,irm,ntim)
      Implicit None
c-----------------------------------------------------------------------
c     calculate the intercell-potentials and add these to the poten-
c     tial v  (in the spin-polarized case for each spin-direction
c     the intercell-potential is the same . )
c     it uses the structure dependent matrices amat and bmat which
c     are calculate once in the subroutine amn .
c     the charge-moments are calculated in the subroutine vintra2 ,
c     therefore vintra2 has to be called first .
c     the intercell-potential is expanded into spherical harmonics .
c     the lm-term of the intercell-potential v of the representive
c     atom i is given by
c
c      v(r,lm,i) =  (-r)**l * {amat(i1,i2,lm,l'm')*cmom(i2,l'm')
c                                               +bmat(i1,i2,lm)*z(i2)}
c
c     summed over i2 (all shells) and l'm' .    (i1=i-natref)
c             (see notes by b.drittler)
c
c     in case of shape correction the madelung potential of the host
c        is taken into account . in all other case the madelung poten-
c        tial of the host is set to be zero .
c        as actual values for z and cmom the differences between the
c        values of the  representive atoms and those of the references
c        are used .
c
c     attention : the first index of cmom (moment of the charge
c                 density - calculated in vintr2) and of z (nuclear
c                 charge of the atoms in the shell) is in the program
c                 different defined , there one has to use :
c                     cmom(natref+i2,l'm')  and  z(natref+i2)
c
c                               b.drittler   june 1987
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD)
      INTEGER LPOTD
      PARAMETER (lpotd=8)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMAX,NATPER,NATREF,NEND,NSPIN,NSTART,
     $     irm, ntim
C     ..
C     .. Array Arguments ..
      REAL*8 AMAT(NTPERD,NTPERD,LMPOTD,*),
     +                 AVMAD(NTREFD,NTREFD,LMPOTD,*),
     +                 BMAT(NTPERD,NTPERD,*),BVMAD(NTREFD,NTREFD,*),
     +                 CMINST(LMPOTD,*),CMOM(LMPOTD,*),R(IRM,*),
     +                 V(IRM,LMPOTD,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IREF(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      REAL*8 AC,PI,SUM
      INTEGER I,I1,I2,IATYP,IATYP2,IPOT,IRS1,ISPIN,L,LM,LM2,LMMAX,M,NREF
C     ..
C     .. Local Arrays ..
      REAL*8 ACH(LMPOTD,NTREFD,2)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
C     .. Save statement ..
      SAVE ACH
C     ..
      PI = 4.D0*DATAN(1.D0)
      LMMAX = (LMAX+1)* (LMAX+1)
      DO 120 IATYP = NSTART,NEND
c---> determine shell index
        I1 = IATYP - NATREF
        IF (KSHAPE.NE.0) THEN
          IRS1 = IRCUT(IPAN(IATYP),IATYP)

        ELSE
          IRS1 = IRWS(IATYP)
        END IF
c
        DO 110 L = 0,LMAX
          DO 100 M = -L,L
            LM = L*L + L + M + 1
            AC = 0.0D0

            IF (KSHAPE.GT.0) THEN
              IF (I1.LT.1) THEN
c
c---> madelung potential - only in case shape correction
c
                IF (NATREF.EQ.1) THEN
c---> lm = 1 component disappears if there is only one host atom
                  DO 10 LM2 = 2,LMMAX
c---> take moments of mt sphere and interstial
                    SUM = CMOM(LM2,1) + CMINST(LM2,1)
                    AC = AC + AVMAD(IATYP,1,LM,LM2)*SUM
   10             CONTINUE

                ELSE

                  DO 30 I2 = 1,NATREF
                    AC = AC + BVMAD(IATYP,I2,LM)*Z(I2)
                    DO 20 LM2 = 1,LMMAX
                      SUM = CMOM(LM2,I2) + CMINST(LM2,I2)
                      AC = AC + AVMAD(IATYP,I2,LM,LM2)*SUM
   20               CONTINUE
   30             CONTINUE
                END IF

                ACH(LM,IATYP,ntim) = AC

              ELSE

c
c---> intercell potential in case of shape corrections
c
                DO 50 I2 = 1,NATPER
c---> determine the reference and representive atom index
                  NREF = IREF(I2)
                  IATYP2 = NATREF + I2
                  DO 40 LM2 = 1,LMMAX
c---> take moments of mt sphere and interstial
                    SUM = CMOM(LM2,IATYP2) - CMOM(LM2,NREF) +
     +                    CMINST(LM2,IATYP2) - CMINST(LM2,NREF)
                    AC = AC + AMAT(I1,I2,LM,LM2)*SUM
   40             CONTINUE
                  AC = AC + BMAT(I1,I2,LM)* (Z(IATYP2)-Z(NREF))
   50           CONTINUE
                AC = AC + ACH(LM,NREF,ntim)

              END IF
c
c---> intercell potential without madelung potential in case of a.s.a.
c
            ELSE IF (I1.GE.1) THEN
              DO 70 I2 = 1,NATPER
c---> determine the reference and representive atom index
                NREF = IREF(I2)
                IATYP2 = NATREF + I2
                DO 60 LM2 = 1,LMMAX
                  SUM = CMOM(LM2,IATYP2) - CMOM(LM2,NREF)
                  AC = AC + AMAT(I1,I2,LM,LM2)*SUM
   60           CONTINUE
                AC = AC + BMAT(I1,I2,LM)* (Z(IATYP2)-Z(NREF))
   70         CONTINUE

            END IF

            IF (LM.EQ.1) WRITE (6,FMT=9000) I1, (AC/SQRT(4.D0*PI))
c
c---> add to v the intercell-potential
c
            DO 90 ISPIN = 1,NSPIN
c
c---> determine the right potential number
c
              IPOT = NSPIN* (IATYP-1) + ISPIN
c
c---> in the case of l=0 : r(1)**l is not defined
c
              IF (L.EQ.0) V(1,1,IPOT) = V(1,1,IPOT) + AC
              DO 80 I = 2,IRS1
                V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,IATYP))**L*AC
   80         CONTINUE
   90       CONTINUE
  100     CONTINUE
  110   CONTINUE
  120 CONTINUE



 9000 FORMAT (1x,'spherically averaged intercell-potential for shell',
     +       i2,' :',1p,d14.6)
      END
