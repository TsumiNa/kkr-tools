      SUBROUTINE BRYDBM(VISP,V,VINS,vspsme,vspsmo,INS,LMPOT,R,DRDI,
     $     ALPHA,ATWGHT,IRC,
     +     IRMIN,NSPIN,NATPS,NATYP,ITDEPT,IMIX,IOBROY,IPF,lsmear)
      Implicit None
c*********************************************************************
c     imix :
c       3      broyden's            f i r s t  m e t h o d
c       4      broyden's          s e c o n d  m e t h o d
c       5      anderson's     g e n e r a l i z e d   m e t h o d
c
c     implemented here according to notes of s.b.
c     broyden's iteration scheme following the papers of :
c     srivastava, j. phys. , 17 (1984) , pp l317
c     c.g. broyden in math.comput., 19 , pp 577, 1965
c     c.g. broyden in ibid, 21 ,pp 368 ,1967
c     the method has been generalized to include a metric.the
c     definition of the necessary inner products are similar to the
c     discription given in the notes of m.weinert.the algorithm
c     discribed in the paper srivastava  has been simplified
c     ( see notes of s.b.)
c     the files ui,vi are stored on high speed ssd memory.
c     broyden's update treats charge and spin on the same footing
c                  s. bluegel , kfa , may 1987
c     the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
c     been generalized and reformulated as an improvement of broyden's
c     second method. successive linesearch is replaced by successive
c     search on hyperplanes. ( see notes of s.b. )
c                  s. bluegel , issp , july 1989
c     the files ui,vi are stored on high speed ssd memory.
c
c     modified for non spherical potential
c                  b. drittler , aug. 1988
c*********************************************************************
C     .. Parameters ..
      INTEGER NATYPD,NTREFD,NTPERD,NSPIND
      PARAMETER (NATYPD=38,NTREFD=1,NTPERD=NATYPD-NTREFD,NSPIND=2)
      INTEGER IRMD,IRNSD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NTIRD
c     now also the smeared spherical potential may be mixed
c      PARAMETER (NTIRD= (IRMD+ (IRNSD+1)* (LMPOTD-1))*NSPIND*NTPERD)
      Parameter (ntird= (2*irmd + (irnsd+1)* (lmpotd-1))*nspind*ntperd)
      INTEGER ITDTHD
      PARAMETER (ITDTHD=40)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALPHA
      INTEGER IMIX,INS,IOBROY,IPF,ITDEPT,LMPOT,NATPS,NATYP,NSPIN,
     $     lsmear
C     ..
C     .. Array Arguments ..
      REAL*8 ATWGHT(*),DRDI(IRMD,*),R(IRMD,*),
     +     V(IRMD,LMPOTD,*),VINS(IRMIND:IRMD,LMPOTD,*),
     +     VISP(IRMD,*),vspsmo(irmd,*),vspsme(irmd,*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      REAL*8 CMM,ONE,RMIXIV,SMNORM,VMDENO,VMNORM,VOLINV,ZERO
      INTEGER IA,IJ,IMAP,IR,IRC1,IRMIN1,ISP,IT,LM,MIT
C     ..
C     .. External Functions ..
      REAL*8 DDOT
      INTEGER MYNODE
      EXTERNAL DDOT
	external MYNODE
C     ..
C     .. External Subroutines ..
      EXTERNAL BRYSH1,BRYSH2,BRYSH3,CREAD,CWRITE,DAXPY,DSCAL,RCSTOP
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
C     .. Save statement ..
c     SAVE MIT,ZERO,ONE,WIT
      SAVE
C     ..
C     .. Local Arrays ..
      REAL*8 AM(2:ITDTHD-1),BM(2:ITDTHD-1),FM(NTIRD),
     +                 FM1(NTIRD),G(NTIRD),SM(NTIRD),SM1(NTIRD),
     +                 VI3(NTIRD),WIT(2:200)
C     ..
C     .. Scalars in Common ..
      LOGICAL STRCO1
      CHARACTER*5 LMACH
C     ..
C     .. Common blocks ..
      COMMON /CHEBY/UI2,UI3,VI2
      COMMON /CMACH/STRCO1,LMACH
C     ..
C     .. Arrays in Common ..
      REAL*8 UI2(NTIRD),UI3(NTIRD),VI2(NTIRD)
C     ..
C     .. Data statements ..
      DATA MIT/1/,ZERO,ONE/0.0D0,1.0D0/
C     ..

      IF (ITDEPT.GT.ITDTHD .OR. ITDTHD.GT.200) CALL RCSTOP('ITDBRY  ')

      IF (IMIX.LE.2 .OR. IMIX.GT.5) CALL RCSTOP('IMIXD   ')

      IF (MIT.GT.ITDEPT) MIT = 1
      IF (IMIX.EQ.3) WRITE (IPF,FMT='('' broyden"s 1st method used '')')
      IF (IMIX.EQ.4) WRITE (IPF,FMT='('' broyden"s 2nd method used '')')
      IF (IMIX.EQ.5) WRITE (IPF,FMT=
     +    '('' generalized anderson method used '')')
c
      RMIXIV = ONE/ALPHA
c
c---->  the following block is activated only one iteration before
c        broyden iteration scheme is used
c---->  set up of : sm1 = rho(1) ; fm1=fm[1]=f(rho(1)) - rho(1) ;
c                   metric  g := r*r*drdi
c---->  map data of all muffin-tin spheres into one single vector
c
      IF (MIT.EQ.1) THEN
c
c
c     order of potentials: v0,v0sme,v1,v2,...
c
        CALL BRYSH3(SM1,VISP,VINS,vspsme,INS,IRMIN,IRC,NATPS,NATYP,
     $       NSPIN,IMAP,LMPOT,lsmear)
        CALL BRYSH1(FM1,V,vspsmo,INS,IRMIN,IRC,NATPS,NATYP,NSPIN,
     $       IMAP,LMPOT,lsmear)

        IF (IMAP.GT.NTIRD) CALL RCSTOP('NIRDBRY ')

        DO 10 IJ = 1,IMAP
          FM1(IJ) = RMIXIV* (FM1(IJ)-SM1(IJ))
   10   CONTINUE
c
        IJ = 0
        DO 60 ISP = 1,NSPIN
          DO 50 IA = NATPS,NATYP
            IRC1 = IRC(IA)
            VOLINV = 3.0D0/ (R(IRC1,IA)**3)
            DO 20 IR = 1,IRC1
              IJ = IJ + 1
              G(IJ) = ATWGHT(IA)*VOLINV*R(IR,IA)*R(IR,IA)*DRDI(IR,IA)
   20       CONTINUE
c
c     next for smeared spherical potential
            If ( lsmear .gt. 0 ) Then
              Do  ir = 1, irc1
                ij = ij + 1
                g(ij) = atwght(ia)*volinv*r(ir,ia)*r(ir,ia)*drdi(ir,ia)
              End Do
            End If
c
            IF (INS.GE.1 .AND. LMPOT.GT.1) THEN
              IRMIN1 = IRMIN(IA)
              DO 40 LM = 2,LMPOT
                DO 30 IR = IRMIN1,IRC1
                  IJ = IJ + 1
                  G(IJ) = ATWGHT(IA)*VOLINV*R(IR,IA)*R(IR,IA)*
     +                    DRDI(IR,IA)
   30           CONTINUE
   40         CONTINUE
            END IF

   50     CONTINUE

   60   CONTINUE
c
        MIT = MIT + 1

      ELSE

c
c----> map rho(m) of all mt-spheres into one single vector
c
        CALL BRYSH3(SM,VISP,VINS,vspsme,INS,IRMIN,IRC,NATPS,NATYP,
     $       NSPIN,IMAP,LMPOT,lsmear)
c
c----> map f[m] = f(m) - rho(m) = f(rho(m)) - rho(m) of all mt-spheres
c      into one single vector
c
        CALL BRYSH1(FM,V,vspsmo,INS,IRMIN,IRC,NATPS,NATYP,NSPIN,
     $       IMAP,LMPOT,lsmear)
        DO 70 IJ = 1,IMAP
          FM(IJ) = RMIXIV* (FM(IJ)-SM(IJ))
   70   CONTINUE
c
c----> calculate  sm = rho(m) - rho(m-1)
c----> calculate dfm = f[m] - f[m-1]
c
        DO 80 IJ = 1,IMAP
          SM1(IJ) = SM(IJ) - SM1(IJ)
          FM1(IJ) = FM(IJ) - FM1(IJ)
   80   CONTINUE
c
c----> loop to generate u[m] = u(ij,mit)
c
        DO 90 IJ = 1,IMAP
          UI3(IJ) = ALPHA*FM1(IJ) + SM1(IJ)
   90   CONTINUE
        REWIND IOBROY
        DO 100 IT = 2,MIT - 1
          IF (LMACH.NE.'INTEL') THEN
            READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP)

          ELSE
            CALL CREAD(IOBROY,UI2,IMAP*8)
            CALL CREAD(IOBROY,VI2,IMAP*8)
          END IF

          AM(IT) = DDOT(IMAP,FM1,1,VI2,1)
          CALL DAXPY(IMAP,-AM(IT),UI2,1,UI3,1)
  100   CONTINUE
c
c----> print amj = the importance of the history of ui
c
        WRITE (IPF,FMT='(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +    MIT - 1, (AM(IT),IT=2,MIT-1)
c
c
        IF (IMIX.EQ.3) THEN
c-------->     b r o y d e n ' s   f i r s t   m e t h o d
c
c
c----> calculate dsmnorm
c
          SMNORM = ZERO
          DO 110 IJ = 1,IMAP
            SMNORM = SMNORM + SM1(IJ)*G(IJ)*SM1(IJ)
  110     CONTINUE
c
c----> convolute dsm with the metric g
c
          DO 120 IJ = 1,IMAP
            SM1(IJ) = G(IJ)*SM1(IJ)
  120     CONTINUE
c
c----> loop to generate v[m] = v(ij,mit)
c
          DO 130 IJ = 1,IMAP
            VI3(IJ) = ALPHA*SM1(IJ)
  130     CONTINUE
          REWIND IOBROY
          DO 140 IT = 2,MIT - 1
            IF (LMACH.NE.'INTEL') THEN
              READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP)

            ELSE
              CALL CREAD(IOBROY,UI2,IMAP*8)
              CALL CREAD(IOBROY,VI2,IMAP*8)
            END IF

            BM(IT) = DDOT(IMAP,SM1,1,UI2,1)
            CALL DAXPY(IMAP,-BM(IT),VI2,1,VI3,1)
  140     CONTINUE
c
c----> complete the evaluation of v[m]
c
          VMDENO = DDOT(IMAP,SM1,1,UI3,1) - SMNORM

          IF (ABS(VMDENO).LT.1D-70) CALL RCSTOP('BRY1SN  ')

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)
c
c----> print bmj = the importance of the history of vi
c
          WRITE (IPF,FMT='(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))'
     +      ) MIT - 1, (BM(IT),IT=2,MIT-1)
c
        ELSE IF (IMIX.EQ.4) THEN
c-------->     b r o y d e n ' s   s e c o n d    m e t h o d
c
c----> calculate v[m] ; convoluted with the metric g
c
          DO 150 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  150     CONTINUE
c
c----> calculate #vm# and normalize v[m]
c
          VMNORM = DDOT(IMAP,VI3,1,FM1,1)
          CALL DSCAL(IMAP,ONE/VMNORM,VI3,1)
c
        ELSE IF (IMIX.EQ.5) THEN
c-------->     g e n e r a l i z e d   a n d e r s o n   m e t h o d
c
c----> calculate v[m] ; convoluted with the metric g
c
          DO 160 IJ = 1,IMAP
            VI3(IJ) = G(IJ)*FM1(IJ)
  160     CONTINUE
          REWIND IOBROY
          DO 170 IT = 2,MIT - 1
            IF (LMACH.NE.'INTEL') THEN
              READ (IOBROY) (UI2(IJ),IJ=1,IMAP), (VI2(IJ),IJ=1,IMAP)

            ELSE
              CALL CREAD(IOBROY,UI2,IMAP*8)
              CALL CREAD(IOBROY,VI2,IMAP*8)
            END IF

            CALL DAXPY(IMAP,-AM(IT)*WIT(IT),VI2,1,VI3,1)
  170     CONTINUE
c
c----> complete the evaluation of v[m]
c
          VMDENO = DDOT(IMAP,FM1,1,VI3,1)

          IF (ABS(VMDENO).LT.1D-70) CALL RCSTOP('BRY1SN  ')

          CALL DSCAL(IMAP,ONE/VMDENO,VI3,1)
c
c----> save wit(mit) for next iteration
c
          WIT(MIT) = VMDENO
c
        END IF
c
c----> write u3(ij) and v3(ij) on disk
c
        IF (LMACH.NE.'INTEL') THEN
          WRITE (IOBROY) (UI3(IJ),IJ=1,IMAP), (VI3(IJ),IJ=1,IMAP)

        ELSE
choshino
c	   IF (MYNODE().EQ.0) THEN
c          CALL CWRITE(IOBROY,UI3,IMAP*8)
c          CALL CWRITE(IOBROY,VI3,IMAP*8)
        END IF
c
c----> update f[m-1] = f[m]  ; rho(m) = rho(m-1)
c
        DO 180 IJ = 1,IMAP
          FM1(IJ) = FM(IJ)
          SM1(IJ) = SM(IJ)
  180   CONTINUE
c
c----> calculate cmm
c
        CMM = DDOT(IMAP,FM,1,VI3,1)
C           WRITE (IPF,FMT='(5X,'' CMM = '',1P,D12.4)') CMM
c
c----> update rho(m+1)
c
        CALL DAXPY(IMAP,ONE-CMM,UI3,1,SM,1)
c
c----> map solution back into each mt-sphere
c
        CALL BRYSH2(SM,V,vspsmo,INS,IRMIN,IRC,NATPS,NATYP,NSPIN,
     $       IMAP,LMPOT,lsmear)
c
        MIT = MIT + 1
        RETURN

  190   CALL RCSTOP('broy10  ')
  200   CALL RCSTOP('broy11  ')
  210   CALL RCSTOP('broy12  ')
  220   CALL RCSTOP('broy13  ')

      END IF


      END
