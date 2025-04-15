      SUBROUTINE FORCXC(FLM,FLMC,LMAX,NSPIN,NSTART,NEND,RHOC,V,R,ALAT,
     +                RM,NSHELL,DRDI,IRWS,NATREF,ND,IOPER,f1xyz)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     calculates the force on nucleus m
c     from a given non spherical charge density at the nucleus site r
c     with core correction(exchange contribution)
c
c     Added the MT boundary correction (needed if rhocore(r_mt) .ne. 0)
c                                   A.Settels and T.Korhonen
c
c     Added the rotation of forces for no symmetry calculation with
c     mdyn.f 
c                                   Holger Hoehler
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NATOMD
      PARAMETER (NATYPD=38,NATOMD=102)
      INTEGER IRMKD,LPOTD
      PARAMETER (irmkd=1484,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALAT
      INTEGER LMAX,NATREF,NEND,NSPIN,NSTART,icount,NMIN,NMAX,
     +        k,I2
C     ..
C     .. Array Arguments ..
      REAL*8 DRDI(IRMKD,*),FLM(-1:1,*),FLMC(-1:1,*),R(IRMKD,*),
     +                 RHOC(IRMKD,*),RM(3,*),V(IRMKD,LMPOTD,*)
      INTEGER IRWS(*),NSHELL(*)
C     ..
C     .. Local Scalars ..
      REAL*8 DV,DVOL,FAC,PI,RWS,TRP,VINT1,VOL
      INTEGER I,IATYP,IPER,IPOT,IREP,IRWS1,ISPIN,J,LM,M
C     ..
C     .. Local Arrays ..
      REAL*8 F(3,NATYPD),FLMH(-1:1,NATYPD),FLMXC(-1:1,NATYPD),
     +                 P(NATYPD),V1(IRMKD),DYN(3),
     +                 DYN0(NATYPD,3),fpart(-1:1)
      REAL*8 f1xyz(3,NATOMD)

C$$$     +                 DYN0(NATYPD,3),fpart(-1:1,natypd)
      INTEGER ND(48,3,3),IOPER(NATOMD)
      integer KATOM,IAT,I1,jj
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP,SIMP3
C     ..
C     .. Save statement ..
      SAVE PI
C     ..
c
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
      PI = 4.0D0*DATAN(1.0D0)
      FAC = SQRT((4.0D0*PI)/3.0D0)
      TRP = 0.0D0
      icount=1
       DO J=1,3
        DYN(J) = 0.d0
       END DO
      IF (LMAX.LT.1) THEN
        WRITE (6,FMT=9000)
        CALL RCSTOP('32      ')

      END IF
c
      WRITE (6,FMT=9040)
      WRITE (6,FMT=9030)
      WRITE (6,FMT=9040)
c
        IREP = 1
choshino       IREP=0
      DO 60 IATYP = NSTART,NEND
c
        fpart = 0.d0
        IPER = IATYP - NATREF
        P(IPER) = 0.0D0
        DO J=1,3
          DYN0(IPER,J) = 0.d0
        END DO
        WRITE (6,FMT=9050) IPER
c
        IRWS1 = IRWS(IATYP)
        RWS = R(IRWS1,IATYP)
        VOL = 0.25D0*ALAT**3
c
c---> determine the right potential numbers
c

        DO 40 M = -1,1
          LM = 2 + M + 1
c
          DO 10 I = 1,IRWS1
            V1(I) = 0.0D0
   10     CONTINUE
c
          DO 30 ISPIN = 1,NSPIN
c
            IPOT = NSPIN* (IATYP-1) + ISPIN
c
c
            DV = (-3.0D0*V(1,LM,IPOT)-10.0D0*V(2,LM,IPOT)+
     +           18.0D0*V(3,LM,IPOT)-6.0D0*V(4,LM,IPOT)+V(5,LM,IPOT))/
     +           (12.0D0*DRDI(2,IATYP))
c
            V1(2) = RHOC(2,IPOT)* (2.0D0*V(2,LM,IPOT)/R(2,IATYP)+DV)/
     +              (4.0D0*PI) + V1(2)
c
            DO 20 I = 3,IRWS1 - 2
c
              DV = (V(I-2,LM,IPOT)-V(I+2,LM,IPOT)+
     +             8.0D0* (V(I+1,LM,IPOT)-V(I-1,LM,IPOT)))/
     +             (12.0D0*DRDI(I,IATYP))
c
              V1(I) = RHOC(I,IPOT)* (2.0D0*V(I,LM,IPOT)/R(I,IATYP)+DV)/
     +                (4.0D0*PI) + V1(I)
   20       CONTINUE
c
            DV = (-V(IRWS1-4,LM,IPOT)+6.0D0*V(IRWS1-3,LM,IPOT)-
     +           18.0D0*V(IRWS1-2,LM,IPOT)+10.0D0*V(IRWS1-1,LM,IPOT)+
     +           3.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1-1,IATYP))
            V1(IRWS1-1) = RHOC(IRWS1-1,IPOT)*
     +                    (2.0D0*V(IRWS1-1,LM,IPOT)/R(IRWS1-1,IATYP)+
     +                    DV)/ (4.0D0*PI) + V1(IRWS1-1)
c
            DV = (3.0D0*V(IRWS1-4,LM,IPOT)-16.0D0*V(IRWS1-3,LM,IPOT)+
     +           36.0D0*V(IRWS1-2,LM,IPOT)-48.0D0*V(IRWS1-1,LM,IPOT)+
     +           25.0D0*V(IRWS1,LM,IPOT))/ (12.0D0*DRDI(IRWS1,IATYP))
c
            V1(IRWS1) = RHOC(IRWS1,IPOT)*
     +                  (2.0D0*V(IRWS1,LM,IPOT)/R(IRWS1,IATYP)+DV)/
     +                  (4.0D0*PI) + V1(IRWS1)
c
c     Correction to the  force comming from the partial integration:
c     If core charge is not zero at the MT boundary, one has to add
c     the following part coming from the partial integration.
c
          fpart(m) = fpart(m) +
     +                  fac/(4.d0*pi)*rhoc(irws1,ipot)*v(irws1,lm,ipot)
c
c
   30     CONTINUE
c
c---> integrate with simpson subroutine
c
          CALL SIMP3(V1,VINT1,1,IRWS1,DRDI(1,IATYP))
c
          FLMH(M,IATYP) = FLM(M,IATYP) - FLMC(M,IATYP)
          FLMXC(M,IATYP) = -FAC*VINT1 - FLMC(M,IATYP)
          FLM(M,IATYP) = FLM(M,IATYP) + FLMXC(M,IATYP)
c

   40   CONTINUE
c
        WRITE (6,FMT=9060) FLMH(1,IATYP),FLMC(1,IATYP),FLMXC(1,IATYP),
     +    FLM(1,IATYP)
        WRITE (6,FMT=9070) FLMH(-1,IATYP),FLMC(-1,IATYP),
     +    FLMXC(-1,IATYP),FLM(-1,IATYP)
        WRITE (6,FMT=9080) FLMH(0,IATYP),FLMC(0,IATYP),FLMXC(0,IATYP),
     +    FLM(0,IATYP)
c
c
c     Core density correction
        Write (6,*)  'Core density correction at MT boundary'
        Write (6,*)  '( It is not added to the total force)'
        Write (6,fmt='(47x,a,d12.6)') 'fccx=',fpart( 1)
        Write (6,fmt='(47x,a,d12.6)') 'fccy=',fpart(-1)
        Write (6,fmt='(47x,a,d12.6)') 'fccz=',fpart( 0)
        Write (6,*)
c
c
        F(1,IATYP) = FLM(1,IATYP)
        F(2,IATYP) = FLM(-1,IATYP)
        F(3,IATYP) = FLM(0,IATYP)
c
c
        DO 50 J = 1,3
          P(IPER) = P(IPER) + RM(J,IREP)*NSHELL(IPER)*F(J,IATYP)*ALAT
   50   CONTINUE

c This part is written for test conditions of mdyncs.f
c It rotates the forces for the representive atoms to
c the forces belonging to the other shell atoms.
c
c unit f1xyz : ryd/a_Bohr 
c
c edited by Holger Hoehler
        
        Do 400 I2 = 1,NSHELL(IPER)
        k=ioper(icount)
        DO 500 I = 1,3
           do 600 J = 1,3
            f1xyz(I,icount)=f1xyz(I,icount)+ND(K,J,I)*F(J,IATYP)
 600       continue 
 500   continue
       icount = icount +1
 400   continue

c rotation ends

        TRP = TRP + P(IPER)
cccccccccccccccccccccccccccccccccccc
         DO 80 IAT=1,NSHELL(IPER)
           KATOM= IOPER(IREP+IAT-1)
      write(6,*) ' iat nst nen',iatyp,nstart,nend
      WRITE(6,*) 'FORCE',IREP+IAT-1,KATOM,IOPER(IREP+IAT-1)   
     +,nshell(iper)
           DO I1=1,3
                  DO JJ=1,3
                  DYN0(IPER,I1) = DYN0(IPER,I1) +
     &            F(JJ,IATYP)* ND(KATOM,JJ,I1)
                  END DO
           END DO

   80    CONTINUE 
           WRITE (6,FMT=9120) (DYN0(IPER,J),J=1,3)
           DO JJ=1,3
                  DYN(JJ) = DYN(JJ) + DYN0(IPER,JJ)
           END DO
ccccccccccccccccccccccccccccccccccc
        IREP = IREP + NSHELL(IPER)
c
        WRITE (6,FMT=9090) P(IPER)
c
   60 CONTINUE
c
      DVOL = TRP/ (3.0D0*VOL)
c
      WRITE (6,FMT=9040)
      WRITE (6,FMT=9010)
      WRITE (6,FMT=9040)
      WRITE (6,FMT=9100) DVOL
      WRITE (6,FMT=9040)
      WRITE (6,FMT=9020)
      WRITE (6,FMT=9110) (DYN(J),J=1,3)
      WRITE (6,FMT=9040)



 9000 FORMAT (13x,'error stop in subroutine force :',
     +       ' the charge density has to contain non spherical',
     +       ' contributions up to l=1 at least !')
 9010 FORMAT (1x,33 ('-'),' volume change ',33 ('-'),/,34x,
     +       ' in units Ry/(a(Bohr)**3 ')
 9020 FORMAT (1x,81 ('-'))
 9030 FORMAT (1x,33 ('-'),' force on the nucleus ',33 ('-'),/,34x,
     +       ' in units Ry/(a(Bohr) ')
 9040 FORMAT (1x,'>')
 9050 FORMAT (3x,i5,'th shell')
 9060 FORMAT (7x,'fhx=',d12.6,2x,'fcx=',d12.6,2x,'fxcx=',d12.6,2x,'fx=',
     +       d12.6,' Ry/(a(Bohr))')
 9070 FORMAT (7x,'fhy=',d12.6,2x,'fcy=',d12.6,2x,'fxcy=',d12.6,2x,'fy=',
     +       d12.6,' Ry/(a(Bohr))')
 9080 FORMAT (7x,'fhz=',d12.6,2x,'fcz=',d12.6,2x,'fxcz=',d12.6,2x,'fz=',
     +       d12.6,' Ry/(a(Bohr))')
 9090 FORMAT (10x,'contribution to the trace of the dipol force tensor:'
     +       ,3x,d12.6,' Ry')
 9100 FORMAT (7x,' volume change dvol/vol=',2x,d12.6,' Ry/(a(Bohr))**3',
     +       /,7x,'( notice: has to be divided',
     +       ' by the bulk modulus of the host)')
 9110 FORMAT (5x,'Sum rule in case of phonon calculation  ',3d12.6)
 9120 FORMAT (3x,'Shell contribution to sum rule ',3d12.6)
 9900 FORMAT (i5,3d14.6,3F10.6)
      END
