        SUBROUTINE CONDUCTWF(E,EK,IE,ISPIN,IRMIN,IRWS,IRNS,LMAX,NSPIN,
     &          NATYP,AR,PNS,PZ,R,KSRA,IATCONDL,IATCONDR,NCONDPAIR,INS)
c
c *******************************************************
c * This subroutine writes out the regular wavefunctions
c * to be used for the conductance formula.
c *                                           3.03.2000
c *******************************************************
c
      implicit none
      include 'inc.fi'
      INTEGER LMMAXD
      PARAMETER(LMMAXD=(LMAXD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      DOUBLE COMPLEX E,EK
      INTEGER I,IE,ISPIN,NSPIN,NATYP,LMAX,IRMIN(*),KSRA,INS
      INTEGER IATCONDL(*),IATCONDR(*),NCONDPAIR,IRNS(*),IRWS(*)
      DOUBLE PRECISION DRDI(IRMD,NATYPD),
     &                 R(IRMD,NATYPD)
      DOUBLE COMPLEX AR(LMMAXD,LMMAXD,NATYPD),
     &               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2,NATYPD),
     &               PZ(IRMD,0:LMAXD,NATYPD),
     &               WF(IRMD),EFAC(LMMAXD)
      INTEGER NS,NS1,IATOM,L1,M1,LM1,IR,LM2,LMMAX,IRWS1,IRMIN1,
     &        IRNS1,iflag
      DOUBLE PRECISION SMIN,RWS
      DOUBLE COMPLEX CONE,V1
      character*7 text
      data text/'ATOM  '/
      CONE = CMPLX(1.d0,0.d0)
c
c
c The sub is called for each energy and spin
c this works only for non relativistic WF.
c scallar relativistic option not implemented!
c
c    
      IF (KSRA.GT.0) STOP 'CONDUCTWF: SRA VERSION NOT AVAILIABLE' 
      LMMAX = (LMAX+1)**2
c
c---> set up array efac : efac(lm) = sqrt(e)**l/(2l - 1)!!
c
c      EFAC(1) = CONE
c      V1 = CONE
c      DO 20 L = 1,LMAX
c        V1 = V1*EK/REAL(2*L-1)
c        DO 10 M = -L,L
c          LM = L* (L+1) + M + 1
c          EFAC(LM) = V1
c   10   CONTINUE
c   20 CONTINUE

      WRITE(66,*) NCONDPAIR,E
      DO NS1 = 1,NCONDPAIR
         DO NS = 1,2
            IF (NS.EQ.1) IATOM =  IATCONDL(NS1)
            IF (NS.EQ.2) IATOM =  IATCONDR(NS1)
            IRWS1 = IRWS(IATOM)
            IRMIN1 = IRMIN(IATOM)
            IRNS1 = IRNS(IATOM)
            SMIN = R(IRMIN1,IATOM)
            RWS  = R(IRWS1,IATOM)
            WRITE(66,*)TEXT,IATOM
            WRITE(66,*)LMAX,IRWS1,IRMIN1,IRNS1,INS

            IF (INS.GT.0) THEN
               WRITE(66,1000) (R(I,IATOM),I=1,IRWS1)
               DO L1=0,LMAX
                  DO IR=1,IRMIN1-1
                     WF(IR) = PZ(IR,L1,IATOM)
                  END DO
                  write(66,*) L1
                  WRITE(66,1000) (WF(I),I=1,IRMIN1-1)
                  DO M1=-L1,L1
                     LM1= L1*(L1+1) + M1 + 1
                     DO LM2=1,LMMAX
                        write(66,*) lm1,lm2,AR(LM1,LM2,IATOM)
                        if (abs(AR(LM1,LM2,IATOM)).gt.1.d-14) then
                           WRITE(66,1000) 
     &                         (PNS(LM1,LM2,IRMIN1+I,1,IATOM),I=0,IRNS1)
                        end if
                     END DO
                  END DO
               END DO
            ELSE
c Spherical wf
            WRITE(66,1000) (R(I,IATOM),I=1,IRWS1)
            DO L1=0,LMAX
                  DO IR=1,IRWS1
                     WF(IR) = PZ(IR,L1,IATOM)
                  END DO
                  write(66,*) L1
                  WRITE(66,1000) (WF(I),I=1,IRWS1)
            END DO
            END IF
c -----------------------------------            
            END DO
         END DO
 999  format(A7,I5)
 1000 FORMAT(4D19.12)
      END 
