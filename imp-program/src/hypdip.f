      SUBROUTINE HYPDIP(AMAT,CMAMV,LMAX,NSPIN,NSTART,NEND,RHO2NS,R,DRDI,
     +                  IRWS,NATREF,NATPER,NSHELL,IOPER,ND,YR,WTYR,RIJ,
     +                  IJEND,RM,ONSITE)
      IMPLICIT NONE
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NATOMD,NTPERD
      PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD)
      INTEGER IJD,NSPIND
      PARAMETER (IJD=434,NSPIND=2)
      INTEGER IRMKD,LPOTD
      PARAMETER (irmkd=1484,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IJEND,LMAX,NATPER,NATREF,NEND,NSPIN,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AMAT(NTPERD,NTPERD,LMPOTD,LMPOTD),
     +                 CMAMV(LMPOTD,NATYPD),DRDI(IRMKD,NATYPD),
     +                 R(IRMKD,NATYPD),RHO2NS(IRMKD,LMPOTD,NATYPD,
     +                 NSPIND),RIJ(IJD,3),RM(3,NATOMD),WTYR(IJD,LMPOTD),
     +                 YR(IJD,LMPOTD),ONSITE(LMPOTD,NATYPD)
      INTEGER IRWS(NATYPD),ND(48,3,3),NSHELL(NTPERD),IOPER(NATOMD)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A0,BR,FAC,HYPFAC,PI,UB,VINT1
      INTEGER I,I1,I2,IA,IATOM,IATYP,IATYP2,IMAX,IMIN,IRWS1,J,L,LM,LM2,
     +        LMMAX,M,M1,M2,NI1R
C     ..
C     .. Local Arrays ..

      DOUBLE PRECISION B(3,3),C(0:LPOTD,-LPOTD:LPOTD,-LPOTD:LPOTD),
     +                 HYPDP(-LMAX:LMAX,NATYPD,3),V1(IRMKD),VROT(3),
     +                 YLM(9)
      REAL*8 TESTV
C     ..
C     .. External Functions ..
C     ..
C     .. External Subroutines ..
      EXTERNAL ROTCOEF,SIMP3,YMY
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
C     ..
C     .. Save statement ..
C     ..
C     .. Data statements ..

      DATA A0/0.529177e-08/,UB/9.2741e-21/
      DATA PI/3.14159265359D0/
C     ..
c
      HYPFAC = UB/ (A0**3)*1.0e-3
      FAC = SQRT(5.0D0/ (4.0D0*PI))
      LMMAX = (LMAX+1)* (LMAX+1)
      IF (LMAX.LT.2) THEN

          STOP 'hypdip'

      ELSE

          DO 10 IATYP = NSTART,NEND
              IRWS1 = IRWS(IATYP)
              DO 20 L = 0,LMAX
                  DO 30 M = -L,L
                      LM = L**2 + L + M + 1
 
                      IF(L.EQ.2) then

                          HYPDP(M,IATYP,1) = ONSITE(LM,IATYP)
                          HYPDP(M,IATYP,2) = HYPDP(M,IATYP,1)

                      END IF
  30             CONTINUE
  20         CONTINUE
  10     CONTINUE


          DO 60 IATYP = NSTART,NEND
              I1 = IATYP - NATREF
              L = 2
              DO 70 M = -L,L
                  LM = 6 + M + 1

                  DO 80 I2 = 1,NATPER
                      IATYP2 = NATREF + I2

                      DO 90 LM2 = 1,LMMAX
                          HYPDP(M,IATYP,1) = HYPDP(M,IATYP,1) +
     +                                   AMAT(I1,I2,LM,LM2)*
     +                                   CMAMV(LM2,IATYP2)
                          
  90                 CONTINUE
  80             CONTINUE

                  HYPDP(M,IATYP,3) = 0.0
                  DO 95 LM2 = 1,LMMAX
                      HYPDP(M,IATYP,3) = HYPDP(M,IATYP,3) +
     +                               AMAT(I1,1,LM,LM2)*CMAMV(LM2,NSTART)
  95             CONTINUE
                  
  70         CONTINUE
  60     CONTINUE


          IF (NSPIN.EQ.2) THEN

              B(1,1) = 0.0
              B(1,2) = 0.0
              B(1,3) = 1.0
              B(2,1) = 1.0
              B(2,2) = 1.0
              B(2,3) = 1.0
              B(3,1) = 1.0
              B(3,2) = 1.0
              B(3,3) = 0.0

              DO 101 I = 1,3

                  WRITE (6,FMT=9000) (B(I,J),J=1,3)

                  CALL YMY(B(I,1),B(I,2),B(I,3),BR,YLM,2)

                  NI1R = 1

                  DO 110 I1 = 1,NATPER
                      WRITE (6,FMT=9010) I1
                      IA = I1 + NATREF
                      IMIN = NI1R
                      NI1R = NI1R + NSHELL(I1)
                      IMAX = NI1R - 1

                      DO 120 IATOM = IMIN,IMAX


                          CALL ROTCOEF(C,IATOM,IOPER,2,ND,YR,WTYR,RIJ,
     +                                 IJEND)

                          DO 125 J = 1,3
                          VROT(J) = 0.0
                          DO 130 M1 = -2,2
                              DO 140 M2 = -2,2
                          VROT(J) = VROT(J) + YLM(6+M1+1)*C(2,M1,M2)*
     +                                   HYPDP(M2,IA,J)*HYPFAC
 140                         CONTINUE
 130                     CONTINUE
 125                     CONTINUE 
                          WRITE (6,FMT=9020) IATOM, (RM(J,IATOM),J=1,3),
     +                      (VROT(J),J=1,3)
 120                 CONTINUE
 110             CONTINUE
 101         CONTINUE

          ELSE
              STOP 'Hypdip:nspin'
          END IF
      END IF

 9000 FORMAT (1x,'>',/,7x,'dipol contribution to hyperfine field ',/,9x,
     +       ' orientation of h field : ',3f6.3,/,1x,'>')
 9010 FORMAT (/,9x,' shell ',i3,/,31x,'dipol contribution in units ',
     +       ' of kGauss ',/)
 9020 FORMAT (13x,' natom : ',i3,/,31x,' rm : ',3f6.3,'total : ',f10.5,
     +        3x,' on site : ',f10.5,3x,' impurity : ',f10.5)
 9040 FORMAT (3x,i5,'th shell')

      END
