      SUBROUTINE CMAM(CMAMV,LMAX,NSPIN,NSTART,NEND,RHO2NS,R,
     +     DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM,IFUNM,
     +     IMAXSH,GSH,THETAS,LMSP,irm,iri,ONSITE)
      Implicit None
c-----------------------------------------------------------------------
c     This routine calculates some integrals for the contribution
c     of the dipolar Hyperfine fields.
c     the lm contribution of the charge moment of the representive
c     atom i is given by
c
c                             rcut
c              cmamv(lm,i) =    s dr' r'** l rho2ns(r',lm,i,2)
c                              0
c           
c                             rcut
c              onsite(lm,i) =    s dr' r'** (-3) rho2ns(r',lm,i,2)
c                              0 
c
c             (see notes by b.drittler and h. hoehler)
c
c              rcut is muffin tin or wigner seitz sphere radius,
c              depending on kshape turned on or off
c
c     attention : rho2ns(...,1) is the real charge density times r**2
c                 developed into spherical harmonics . (see deck rholm)
c
c                               b.drittler   may 1987
c                               h.hoehler    nov 1999
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMD,IRMKD,LPOTD
      PARAMETER (irmd=1484,irmkd=1484,lpotd=8)
      INTEGER NFUND,NGSHD
      PARAMETER (NFUND=289,NGSHD=54287)
      Integer irmaxd
      Parameter ( irmaxd=max0(irmkd,irmd) )
c
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMAX,NEND,NSPIN,NSTART, irm, iri,IRWIG
C     ..
C     .. Array Arguments ..
      REAL*8 CMAMV(LMPOTD,*),DRDI(IRM,*),ONSITE(LMPOTD,*),
     +                 GSH(*),R(IRM,*),RHO2NS(IRM,LMPOTD,NATYPD,*),
     +                 THETAS(IRI,NFUND,*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +        IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 FAC,PI,RL
      INTEGER I,IATYP,ICELL,IEND,IFUN,IPOT,IRC1,IRS1,ISTART,J,L,LM,LM2,
     +        LM3,M
C     ..
C     .. Local Arrays ..
      REAL*8 V1(IRMaxD),V2(IRMaxD),VINT1(IRMaxD),VINT2(IRMaxD)
      INTEGER IRCUTM(0:IPAND)
C     ..
C     .. External Subroutines ..
      EXTERNAL SINWK,SOUTK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,REAL
      SAVE
C     ..
      PI = 4.D0*DATAN(1.D0)
c
      DO 110 IATYP = NSTART,NEND
        IF (KSHAPE.NE.0) THEN
          IRS1 = IRCUT(1,IATYP)
          IRC1 = IRCUT(IPAN(IATYP),IATYP)
          ICELL = NTCELL(IATYP)

          DO 10 I = 0,IPAN(IATYP)
            IRCUTM(I) = IRCUT(I,IATYP)
   10     CONTINUE

        ELSE

          IRS1 = IRWS(IATYP)
          IRC1 = IRS1
          IRCUTM(0) = IRCUT(0,IATYP)
          IRCUTM(1) = IRC1
        END IF

        DO 100 L = 0,LMAX
          FAC = 8.0D0*PI/REAL(2*L+1)
          DO 90 M = -L,L
            LM = L*L + L + M + 1
c
c---> set up of the integrands v1 and v2
c
            V1(1) = 0.0D0
            V2(1) = 0.0D0
            DO 20 I = 2,IRS1
              RL = R(I,IATYP)**L
              V1(I) = RHO2NS(I,LM,IATYP,2)*RL*DRDI(I,IATYP)
              V2(I) = RHO2NS(I,LM,IATYP,2)/R(I,IATYP)/RL*DRDI(I,IATYP)
   20       CONTINUE
c
c---> convolute charge density of interstial with shape function
c        if kshape.gt.0
c
            IF (KSHAPE.NE.0) THEN
              DO 30 I = IRS1 + 1,IRC1
                V1(I) = 0.0D0
   30         CONTINUE
              ISTART = IMAXSH(LM-1) + 1
              IEND = IMAXSH(LM)
              DO 50 J = ISTART,IEND
                LM2 = ILM(J,2)
                LM3 = ILM(J,3)
                IF (LMSP(ICELL,LM3).GT.0) THEN
                  IFUN = IFUNM(ICELL,LM3)
                  DO 40 I = IRS1 + 1,IRC1
                    V1(I) = V1(I) + GSH(J)*RHO2NS(I,LM2,IATYP,2)*
     +                      THETAS(I-IRS1,IFUN,ICELL)
   40             CONTINUE
                END IF
   50         CONTINUE

              DO 60 I = IRS1 + 1,IRC1
                RL = R(I,IATYP)**L
                V2(I) = V1(I)/R(I,IATYP)/RL*DRDI(I,IATYP)
                V1(I) = V1(I)*RL*DRDI(I,IATYP)*DRDI(I,IATYP)
   60         CONTINUE
            END IF
c
c---> now integrate v1 and v2
c
            IRWIG=IRWS(IATYP) 
            CALL SOUTK(V1,VINT1,IPAN(IATYP),IRCUTM)
            CALL SOUTK(V2,VINT2,IPAN(IATYP),IRCUTM) 

c
c---> gather all parts
c

c
c---> store charge moment - in case of kshape.gt.0 this is the moment
c      of the charge in the muffin tin sphere
c
            CMAMV(LM,IATYP) = VINT1(IRS1)            
            ONSITE(LM,IATYP) = VINT2(IRS1)*8.0D+00*PI/5.0D+00
            IF(KSHAPE.NE.0) THEN
              CMAMV(LM,IATYP) = VINT1(IRC1)            
              ONSITE(LM,IATYP) = VINT2(IRC1)*8.0D+00*PI/5.0D+00
            END IF 
c
c---> store charge moment of interstial in case of kshape.gt.0
c
c


   90     CONTINUE

  100   CONTINUE

  110 CONTINUE


      END
