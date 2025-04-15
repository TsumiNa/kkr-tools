      SUBROUTINE MIXSTR(RMSAVQ,RMSAVM,INS,LPOT,LMPOT,NATREF,NSHELL,
     +     NSTART,NEND,NSPIN,ITC,RFPI,FPI,IPF,MIXING,FCM,
     +     IRC,IRMIN,R,DRDI,VONS,VISP,VINS,vspsmo,vspsme,lsmear)
      Implicit None
c-----------------------------------------------------------------------
c
c     straight-mixing scheme
c     with determination of root-mean-square error
c
c
c                            R. Zeller   Jan. 1993
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER IRMD,IRNSD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lpotd=8)
      INTEGER LMPOTD,NPOTD,IRMIND
      PARAMETER (LMPOTD= (LPOTD+1)**2,NPOTD=NSPIND*NATYPD,
     +          IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      REAL*8 FCM,FPI,MIXING,RFPI,RMSAVM,RMSAVQ
      INTEGER INS,IPF,ITC,LMPOT,LPOT,NATREF,NEND,NSPIN,NSTART,lsmear
C     ..
C     .. Array Arguments ..
      REAL*8 DRDI(IRMD,NATYPD),R(IRMD,NATYPD),
     +     VINS(IRMIND:IRMD,LMPOTD,NPOTD),VISP(IRMD,NPOTD),
     +     VONS(IRMD,LMPOTD,NPOTD),vspsmo(irmd,npotd),
     $     vspsme(irmd,npotd)
      INTEGER IRC(NATYPD),IRMIN(NATYPD),NSHELL(NTPERD)
C     ..
C     .. Local Scalars ..
      REAL*8 FAC,RMSERM,RMSERQ,VMN,VNM,VNP,VOLDM,VOLDP,VPN
      INTEGER I,IH,IHP1,IRC1,IRMIN1,J,LM,NATOM,NP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MOD,REAL,SQRT
C     ..
      RMSAVQ = 0.0D0
      RMSAVM = 0.0D0
c
c---> final construction of the potentials
c     attention : the spherical averaged potential is the lm=1
c                     component of vons times sqrt(4 pi).
c
c     first mixing scheme : straight mixing
c---> determination of the root mean sqare error
c
      NATOM = 0
      DO 40 NP = NSTART,NEND
        I = NP - NATREF
        NATOM = NATOM + NSHELL(I)
        IF (NSPIN.EQ.2) THEN
          IH = 2*NP - 1
          IHP1 = IH + 1
        ELSE


          IH = NP
          IHP1 = IH
        END IF

        IRC1 = IRC(NP)
        RMSERQ = 0.0D0
        RMSERM = 0.0D0
        FAC = 0.5D0/RFPI
c     
        DO 10 J = 1,IRC1
          VNP = FAC* (VONS(J,1,IH)+VONS(J,1,IHP1))
          VNM = FAC* (VONS(J,1,IH)-VONS(J,1,IHP1))
          VOLDP = 0.5D0* (VISP(J,IH)+VISP(J,IHP1))
          VOLDM = 0.5D0* (VISP(J,IH)-VISP(J,IHP1))
          RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J,NP)*R(J,NP)*
     +             DRDI(J,NP)* (VNP-VOLDP)* (VNP-VOLDP)
          RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J,NP)*R(J,NP)*
     +             DRDI(J,NP)* (VNM-VOLDM)* (VNM-VOLDM)
          VPN = VOLDP + MIXING* (VNP-VOLDP)
          VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
          VONS(J,1,IHP1) = VPN - VMN
          VONS(J,1,IH) = VPN + VMN
   10   CONTINUE

        If ( lsmear .ge. 3 ) Then
          Do j = 1,irc1
            vnp   = 0.5d0* (vspsmo(j,ih)+vspsmo(j,ihp1))
            vnm   = 0.5d0* (vspsmo(j,ih)-vspsmo(j,ihp1))
            voldp = 0.5d0* (vspsme(j,ih)+vspsme(j,ihp1))
            voldm = 0.5d0* (vspsme(j,ih)-vspsme(j,ihp1))
            vpn   = voldp + mixing* (vnp-voldp)
            vmn   = voldm + fcm*mixing* (vnm-voldm)
            vspsmo(j,ihp1) = vpn - vmn
            vspsmo(j,ih)   = vpn + vmn
          End Do
        End If
          
        If ( (lsmear .eq. 1) .or. (lsmear .eq. 2) ) Then
          Do j = 1,irc1
            vspsme(j,ihp1) = vspsmo(j,ihp1)
            vspsme(j,ih)   = vspsmo(j,ih)
          End Do
        End If
          

        RMSERQ = RMSERQ/ (R(IRC1,NP)**3)
        RMSERM = RMSERM/ (R(IRC1,NP)**3)
        RMSAVQ = RMSAVQ + RMSERQ*NSHELL(I)
        RMSAVM = RMSAVM + RMSERM*NSHELL(I)
        IF (NSPIN.EQ.2) THEN
          WRITE (IPF,FMT=9000) I,SQRT(RMSERQ),SQRT(RMSERM)
        ELSE


          WRITE (IPF,FMT=9020) I,SQRT(RMSERQ)
        END IF

        IF (INS.NE.0 .AND. LPOT.GT.0) THEN

          RMSERQ = 0.0D0
          RMSERM = 0.0D0
          IRMIN1 = IRMIN(NP)
          DO 30 LM = 2,LMPOT
            DO 20 J = IRMIN1,IRC1
              VNP = 0.5D0* (VONS(J,LM,IH)+VONS(J,LM,IHP1))
              VNM = 0.5D0* (VONS(J,LM,IH)-VONS(J,LM,IHP1))
              VOLDP = 0.5D0* (VINS(J,LM,IH)+VINS(J,LM,IHP1))
              VOLDM = 0.5D0* (VINS(J,LM,IH)-VINS(J,LM,IHP1))
              RMSERQ = RMSERQ + 2.0D0*REAL(1+MOD(J,2))*R(J,NP)*R(J,NP)*
     +                 DRDI(J,NP)* (VNP-VOLDP)* (VNP-VOLDP)
              RMSERM = RMSERM + 2.0D0*REAL(1+MOD(J,2))*R(J,NP)*R(J,NP)*
     +                 DRDI(J,NP)* (VNM-VOLDM)* (VNM-VOLDM)
              VPN = VOLDP + MIXING* (VNP-VOLDP)
              VMN = VOLDM + FCM*MIXING* (VNM-VOLDM)
              VONS(J,LM,IHP1) = VPN - VMN
              VONS(J,LM,IH) = VPN + VMN
   20       CONTINUE
   30     CONTINUE
          RMSERQ = RMSERQ/ (R(IRC1,NP)**3)/FPI
          RMSERM = RMSERM/ (R(IRC1,NP)**3)/FPI
          RMSAVQ = RMSAVQ + RMSERQ*NSHELL(I)
          RMSAVM = RMSAVM + RMSERM*NSHELL(I)
          IF (NSPIN.EQ.2) THEN
            WRITE (IPF,FMT=9010) I,SQRT(RMSERQ),SQRT(RMSERM)
          ELSE


            WRITE (IPF,FMT=9030) I,SQRT(RMSERQ)
          END IF

        END IF
   40 CONTINUE


      RMSAVQ = SQRT(RMSAVQ/REAL(NATOM))
      RMSAVM = SQRT(RMSAVM/REAL(NATOM))
      IF (NSPIN.EQ.2) THEN
        WRITE (IPF,FMT=9040) ITC,RMSAVQ,RMSAVM
      ELSE


        WRITE (IPF,FMT=9050) ITC,RMSAVQ
      END IF

 9000 FORMAT (5x,' rms-error for cell',i3,1x,':','v+ + v- = ',1p,d11.4,
     +       2x,',',2x,'v+ - v- = ',1p,d11.4)
 9010 FORMAT (5x,' rms-error non spherical contribution for cell ',i3,
     +       1x,':','v+ + v- = ',1p,d11.4,02x,',',2x,'v+ - v- = ',1p,
     +       d11.4)
 9020 FORMAT (5x,' rms-error for cell',i3,1x,':','v+ + v- = ',1p,d11.4)
 9030 FORMAT (5x,' rms-error non spherical contribution for cell ',i3,
     +       1x,':','v+ + v- = ',1p,d11.4)
 9040 FORMAT (10x,'iteration no.:',i4,' average rms-error : v+ + v- = ',
     +       1p,d11.4,/,48x,' v+ - v- = ',1p,d11.4)
 9050 FORMAT (10x,'iteration no.:',i4,' average rms-error : v+ + v- = ',
     +       1p,d11.4)
      END
