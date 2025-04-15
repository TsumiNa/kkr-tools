      SUBROUTINE ECOULOM(CMOM,ECOU,LMAX,NSPIN,NSTART,NEND,RHO2NS,VM2Z,Z,
     +     R,DRDI,IRWS,KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,
     +     IFUNM,ILM,NTCELL,GSH,THETAS,VBC,LMSP)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     attention : new energy zero ---> muffin tin zero
c
c     calculate the electrostatic potential-energies without the
c     electron-nuclear interaction in the cell itself .
c     the energy of the representive atom i is given by
c
c                          rc
c      ecou(i) =  1/2 (  {  s dr' vm2z(r',lm,i)*rho2ns(r',lm,i,1) }
c                           0
c
c                                       -  z(i) * vmad ( ri )     )
c
c                          rc
c                   +       s dr' vbc*rho2ns(r',lm,i)
c                           0
c
c                                         ( {..} = summed over lm )
c             (see notes by b.drittler)
c     vm2z is the coulomb potential of the atom without the nuclear
c             potential of the atom
c     rho2ns(...,1) is the real charge density times r**2
c
c      both developed into spherical harmonics . (see deck rholm)
c
c     z    is the nuclear charge of the atom
c
c     vmad ( ri ) is a generalized madelung potential
c                 = 1/sqrt(4 pi) * vm2z(irws,1,is)
c                         - sqrt(4 pi) * 2 * cmom(1,ipot) / rws
c
c                                        ( <..> = spherical averaged )
c
c     vbc is the constant shift of the potential from electrostatic
c         zero to muffin tin zero - this has to be treated as external
c         potential ! in case of spin polarisation vbc is spin - de-
c         pendend therefore it has to be integrated with the correct
c         spin density
c
c     attention : this subroutine has to be called bevor the
c                 exchange correlation potential is added to
c                 the potential vm2z .
c                 the energy calculated here is splitted into
c                 l-dependent parts to see the l -convergency .
c
c     attention : in case of shape corrections the contribution of
c                 the coulomb potential the of the nucleus is
c                 analytically cancelled only in the muffin tin sphere
c                 in the interstial region it has to be taken into
c                 account ! see deck epotins
c
c                               b.drittler   oct. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER IRMKD,LPOTD
      PARAMETER (irmkd=1484,lpotd=8)
      INTEGER NFUND,IRIKD,NGSHD
      PARAMETER (NFUND=289,irikd=435,NGSHD=54287)
      INTEGER IPAND
      PARAMETER (IPAND=80)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,KVMAD,LMAX,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      REAL*8 CMOM(LMPOTD,*),DRDI(IRMKD,*),ECOU(0:LPOTD,*),
     +     GSH(*),R(IRMKD,*),RHO2NS(IRMKD,LMPOTD,NATYPD,*),
     +     THETAS(IRIKD,NFUND,*),VBC(*),VM2Z(IRMKD,LMPOTD,*),
     +     Z(*)
      INTEGER IFUNM(NATYPD,*),ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
     +     IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 PI,RFPI,RHOSP,SIGN,VM,VMAD
      INTEGER I,IATYP,ICELL,IFUN,IPAN1,IPOT,IR,IRC1,IRH,IRS1,ISPIN,J,L,
     +     LM,LM2,M
C     ..
C     .. Local Arrays ..
      REAL*8 ER(IRMKD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMP3,SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
      PI = 4.D0*DATAN(1.D0)
      RFPI = SQRT(4.0D0*PI)
c
      DO 110 IATYP = NSTART,NEND
        IF (KSHAPE.NE.0) THEN
          IPAN1 = IPAN(IATYP)
          ICELL = NTCELL(IATYP)
          IRS1 = IRCUT(1,IATYP)
          IRC1 = IRCUT(IPAN1,IATYP)

        ELSE

          IRS1 = IRWS(IATYP)
          IRC1 = IRS1
        END IF
c
c---> determine the right potential numbers - the coulomb potential
c      is not spin dependend
c
        IPOT = IATYP*NSPIN

        DO 90 L = 0,LMAX
          DO 10 I = 1,IRC1
            ER(I) = 0.0D0
   10     CONTINUE
          DO 80 ISPIN = 1,NSPIN
            IF (ISPIN.EQ.NSPIN) THEN
              SIGN = 1.0D0

            ELSE
              SIGN = -1.0D0
            END IF

            DO 70 M = -L,L
              LM = L*L + L + M + 1
              IF (LM.EQ.1) THEN
c
c---> treat shift vbc as external potential
c
                DO 20 I = 1,IRS1
                  RHOSP = (RHO2NS(I,LM,IATYP,1)+
     +                    SIGN*RHO2NS(I,LM,IATYP,NSPIN))/2.0D0
                  ER(I) = ER(I) + RHOSP* (VM2Z(I,LM,IPOT)/2.0D0+
     +                    VBC(ISPIN)*RFPI)
   20           CONTINUE

              ELSE

                DO 30 I = 1,IRS1
                  RHOSP = (RHO2NS(I,LM,IATYP,1)+
     +                    SIGN*RHO2NS(I,LM,IATYP,NSPIN))/2.0D0
                  ER(I) = ER(I) + RHOSP*VM2Z(I,LM,IPOT)/2.0D0
   30           CONTINUE
              END IF

              IF (KSHAPE.NE.0) THEN
c
c---> convolute with shape function
c
                DO 60 J = IMAXSH(LM-1) + 1,IMAXSH(LM)
                  LM2 = ILM(J,2)
                    IFUN = IFUNM(ICELL,ILM(J,3))
                    IF (LM2.EQ.1) THEN
                      DO 40 IR = IRS1 + 1,IRC1
                        IRH = IR - IRS1
                        RHOSP = (RHO2NS(IR,LM,IATYP,1)+
     +                          SIGN*RHO2NS(IR,LM,IATYP,NSPIN))/2.0D0
c
c---> remember that in the interstial -2z/r has to be taken into account
c
                        IF (LMSP(ICELL,ILM(J,3)).GT.0) THEN
                        ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                           THETAS(IRH,IFUN,ICELL)*
     +                           (VM2Z(IR,1,IPOT)/2.0D0+
     +                           (VBC(ISPIN)-Z(IATYP)/R(IR,IATYP))*RFPI)
                        END IF
   40                 CONTINUE

                    ELSE

                      DO 50 IR = IRS1 + 1,IRC1
                        IRH = IR - IRS1
                        RHOSP = (RHO2NS(IR,LM,IATYP,1)+
     +                          SIGN*RHO2NS(IR,LM,IATYP,NSPIN))/2.0D0
                        IF (LMSP(ICELL,ILM(J,3)).GT.0) THEN
                        ER(IR) = ER(IR) + RHOSP*GSH(J)*
     +                           THETAS(IRH,IFUN,ICELL)*
     +                           VM2Z(IR,LM2,IPOT)/2.0D0
                        END IF
   50                 CONTINUE
                    END IF
   60           CONTINUE

              END IF

   70       CONTINUE

   80     CONTINUE
c
c---> now integrate
c
          IF (KSHAPE.EQ.0) THEN
            CALL SIMP3(ER,ECOU(L,IATYP),1,IRS1,DRDI(1,IATYP))

          ELSE

            CALL SIMPK(ER,ECOU(L,IATYP),IPAN1,IRCUT(0,IATYP),
     +                 DRDI(1,IATYP))
          END IF

   90   CONTINUE

c
c---> calculate the madelung potential
c
        VMAD = VM2Z(IRS1,1,IPOT)/RFPI -
     +         RFPI*2.0D0*CMOM(1,IATYP)/R(IRS1,IATYP)
c
c---> add to ecou
c
        ECOU(0,IATYP) = ECOU(0,IATYP) - Z(IATYP)*VMAD/2.0D0
c
c---> option to calculate full generalized madelung potential
c                                  rc
c     vm(rn) = vmad +2*sqrt(4*pi)* s  dr*r*rho(lm=1,r)
c                                  0
        IF (KVMAD.EQ.1) THEN
          ER(1) = 0.0D0
          DO 100 I = 2,IRS1
            ER(I) = RHO2NS(I,1,IATYP,1)/R(I,IATYP)
  100     CONTINUE
          CALL SIMP3(ER,VM,1,IRS1,DRDI(1,IATYP))
          VM = 2.0D0*RFPI*VM + VMAD
c
c  atom nr. iatyp is the iatyp-th atom on the potential cards
c  e. g., in binary alloys iatyp=1 and iatyp=2 refer to host
c
          WRITE (6,FMT=9000) IATYP,VM

        END IF

  110 CONTINUE



 9000 FORMAT (13x,'full generalized madelung pot. for atom',1x,i3,1x,
     +       ': ',D14.8)
      END
