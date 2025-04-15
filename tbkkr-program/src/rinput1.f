c ************************************************************************
      SUBROUTINE RINPUT1(ALAT,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3,ESHIFT,
     +                  FCM,
     +                  HFIELD,MIXING,QBOUND,VCONST,KFG,LMXC,IRNS,
     +                  NTCELL,
     +                  ICST,IFILE,IGF,IMIX,INS,INSREF,IPE,IPF,IPFE,
     +                  IPOTOU,
     +                  IPRCOR,IRM,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,
     +                  KCOR,KEFG,KFROZN,KHFELD,KHYP,KPRE,KSHAPE,KTE,
     +                  KVMAD,KVREL,KWS,KXC,LMAX,LMMAX,LMPOT,LPOT,MD,
     +                  NATYP,NSPIN,NATYPD,NSPIND,IEMXD,IRMD,IRNSD,
     +                  LMAXD,LPOTD,TXC,KSCOEF,ICC)
      implicit none
c ************************************************************************
C     .. Local Arrays ..
      CHARACTER*4 TSPIN(2)
      CHARACTER*8 TKWS(3)
      CHARACTER*43 TINS(0:3),TKCOR(0:3),TVREL(0:2)
C     ..
C     .. External Subroutines ..
      EXTERNAL RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MIN
C     ..
C     .. Array Arguments ..
      INTEGER IRNS(*),KFG(4,*),LMXC(*),NTCELL(*)
      CHARACTER*24 TXC(3)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,E1,E2,ESHIFT,FCM,HFIELD,MIXING,QBOUND,TK,
     +       VCONST
      INTEGER ICC,ICST,IEMXD,IFILE,IGF,IMIX,INS,INSREF,
     +        IPE,IPF,IPFE,IPOTOU,IPRCOR,
     +        IRM,IRMD,IRNSD,IRNUMX,ISHIFT,ITCCOR,ITCLST,ITDBRY,KCOR,
     +        KEFG,KFROZN,KHFELD,KHYP,KPRE,KSCOEF,KSHAPE,KTE,KVMAD,
     +        KVREL,KWS,KXC,LMAX,LMAXD,LMMAX,LMPOT,LPOT,LPOTD,MD,
     +        NATYP,NATYPD,NPNT1,NPNT2,NPNT3,NPOL,NSPIN,NSPIND
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION BRYMIX,STRMIX
      INTEGER I,IL,IP
      CHARACTER*43 TSHAPE
C     ..
C     .. Data statements ..
      DATA TSPIN/'non-','    '/
      DATA TSHAPE/' exact cell treatment (shape correction)  '/
      DATA TVREL/
     +     ' non relativistic calculation              ',
     +     ' s.r.a. calculation                        ',
     +     ' s.r.a. calculation                        '/
      DATA TKCOR/
     +     ' frozen core approximation                 ',
     +     ' core relaxation s.r.a.                    ',
     +     ' core relaxation nonsra                    ',
     +     ' core relaxation                           '/
      DATA TINS/' spherical averaged input potential        ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential for cluster ',
     +     ' non spherical input potential             '/
      DATA TKWS/' full mt','   ws   ',' full ws'/
C     ..
c
c------------ array set up and definition of input parameter -----------
c
      TXC(1) = ' Morruzi,Janak,Williams '
      TXC(2) = ' von Barth,Hedin        '
      TXC(3) = ' Vosko,Wilk,Nusair      '
c
c---> read input
c
      read (7,FMT=9180) LMAX
      read (7,FMT=9190) E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3
      read (7,FMT=9060) IRNUMX,ITCCOR,IPRCOR
      read (7,FMT=9200) IFILE,IPE,ISHIFT,ESHIFT
      read (7,FMT=9060) KSHAPE,IRM,INS,ICST,INSREF
      read (7,FMT=9060) KCOR,KVREL,KWS,KHYP,KHFELD,KXC
      KFROZN = KCOR
      IF (KCOR.EQ.0) KCOR = 2
      read (7,FMT=9060) KTE,KPRE,KEFG,KVMAD,KSCOEF
c      read (7,FMT=9060) NATYP,KSHAPE
c      read (7,FMT=9060) (NTCELL(I),I=1,NATYP)
c      read (7,FMT=9060) (IRNS(I),I=1,NATYP)
      read (7,FMT=9060) IMIX,IPOTOU,IGF,ICC
      read (7,FMT=9060) ITDBRY
      read (7,FMT=9040) STRMIX,FCM,QBOUND
      read (7,FMT=9040) BRYMIX
      read (7,FMT=9040) HFIELD,VCONST
c ------------------------------------------------------------------------
      WRITE (6,9210) LMAX
      WRITE (6,9301)
      WRITE (6,9220) E1,E2,TK
      WRITE (6,9302)
      WRITE (6,9230) NPOL,NPNT1,NPNT2,NPNT3
      WRITE (6,9304)
      WRITE (6,9240) IRNUMX,ITCCOR,IPRCOR
      WRITE (6,9303)
      WRITE (6,9250) IFILE,IPE,ISHIFT,ESHIFT
      WRITE (6,9305)
      WRITE (6,9260) KSHAPE,IRM,INS,ICST,INSREF
      WRITE (6,9309)
      WRITE (6,9270) KCOR,KVREL,KWS,KHYP,KHFELD,KXC
      WRITE (6,9306)
      WRITE (6,9330) KTE,KPRE,KEFG,KVMAD,KSCOEF
      WRITE (6,9309)
      WRITE (6,9290) IMIX,IPOTOU,IGF,ICC
      WRITE (6,9304)
      WRITE (6,9300) ITDBRY
      WRITE (6,9307)
      WRITE (6,9310) STRMIX,FCM,QBOUND
      WRITE (6,9302)
      WRITE (6,9320) BRYMIX
      WRITE (6,9308)
      WRITE (6,9280) HFIELD,VCONST
c ------------------------------------------------------------------------
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9240 FORMAT (' IRNUMX ITCCOR IPRCOR'/,3i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT ESHIFT'/,3i7,f12.6)
 9260 FORMAT (' KSHAPE    IRM    INS   ICST INSREF'/,5i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHYP KHFELD    KXC'/,6i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,
     +        ' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX IPOTOU    IGF    ICC'/,4i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KEFG  KVMAD KSCOEF'/,5i7)
 9301 format(   3(1H-),1H+  ,75(1H-))
 9302 format( 3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
c ------------------------------------------------------------------------
c      DO 10 IP = 1,NATYP
c        READ (5,FMT=9010) LMXC(IP), (KFG(IL,IP),IL=1,4)
c        WRITE (6,FMT=9010) LMXC(IP), (KFG(IL,IP),IL=1,4)
c   10 CONTINUE
c
      IF (KSHAPE.NE.0) KWS = 2
c
      IPF = 6
      IPFE = IPF + 3
c
      IF (QBOUND.LT.1.D-15) QBOUND = 1.D-4
      IF (IMIX.GT.2) THEN
        FCM = 1.0D0
        MIXING = BRYMIX
      ELSE
        MIXING = STRMIX
      END IF
c
      IF (IMIX.GE.6) WRITE (6,FMT=9110) (IMIX-5),ITDBRY - 1
c
      WRITE (6,FMT=9090) MIXING,QBOUND
c
      LMMAX = (LMAX+1)**2
      LPOT  = MIN(2*LMAX,LPOTD)
      LMPOT = (LPOT+1)* (LPOT+1)
c
      WRITE (6,FMT=9020) LMAX,LMAXD,NATYP,NATYPD,IRM,IRMD,NSPIN,NSPIND

c      IF (LMAX.GT.LMAXD .OR. NATYP.GT.NATYPD .OR. IRM.GT.IRMD .OR.
c     +    NSPIN.GT.NSPIND) CALL RCSTOP('18      ')

      IF (INS.GT.0) THEN
        WRITE (6,FMT=9130)
        WRITE (6,FMT=9140)
        DO 20 I = 1,NATYP
          WRITE (6,FMT=9150) I,IRNS(I),IRNSD

          IF (IRNS(I).GT.IRNSD) CALL RCSTOP('19      ')

   20   CONTINUE

        IF (LMAX.NE.LMAXD) THEN
          WRITE (6,FMT=9120)

          CALL RCSTOP('20      ')

        END IF

      END IF


      WRITE (6,FMT=9130)
c
c
c
c
      IF (KHFELD.EQ.1) WRITE (6,FMT=9030) HFIELD
      WRITE (6,FMT=9050) TSPIN(NSPIN)
      WRITE (6,FMT=9170) TVREL(KVREL)
      WRITE (6,FMT=9170) TKCOR(KFROZN)
      IF (KSHAPE.EQ.0) THEN
        WRITE (6,FMT=9070) TKWS(KWS+1)

      ELSE
        WRITE (6,FMT=9170) TSHAPE
      END IF

      WRITE (6,FMT=9100) TXC(KXC+1)
      IF (INS.GT.0) WRITE (6,FMT=9160) TINS(INS),ICST
      WRITE (6,FMT=9080)

      RETURN


 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,
     +       35 ('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,
     +       'natyp  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,
     +       ')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',
     +       f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,
     +        ' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,
     +       ' is used up to iteration-      ',/,20x,'depth :',i3,
     +       '  then jacobian is fixed and potential      ',/,20x,
     +       'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts ',
     +       'the parameter lmaxd has to be set equal lmax ! ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation ',
     +       '- cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
      END                           ! RINPUT1
