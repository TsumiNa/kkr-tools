
      SUBROUTINE START(IFILE,IPF,IPE,IPFE,KVREL,KWS,LMAX,NEND,ALAT,RMT,
     +     RMTNEW,ITITLE,IMT,RFCTOR,NATREF,VCONSTH,INS,IRNS,
     +     LPOT,NSPIN,VINS,IRMIN,KSHAPE,NTCELL,IRCUT,IPAN,
     +     THETAS,theeqm,thesme,alpha,vspsme,lsmear,
     $     IFUNM,NFU,LLMSP,LMSP,ECORE,LCORE,NCORE,
     $     r,drdi,droreq,rseq,irws,irmieq,ircueq,ipaneq,
     $     dror,rs,vbc)
      Implicit None
c-----------------------------------------------------------------------
c   reads the input potentials
c
c    units :       ry - units for energy
c                  the lattice constant and all other lengths
c                                           given in bohr units
c                  the planck constant h/2pi=1
c                  the electron charge e=sqrt(2)
c                  the electron mass m=1/2
c                  the speed of light c = 2/alpha = 274.0720442
c                      with the fine structure constant alpha
c
c    remember that the input potentials do not include the electro-
c             static contribution of the nucleus of the cell itself
c             this has to be added explicitly !
c
c   as input is used: lmax=maximum angular momentum
c                    nend=number of different atoms
c
c
c     in case of shape corrections this routine  reads from unit 19
c     a suitable radial  mesh 'xrn',its derivate 'drn' and the shape
c     functions 'thetas' .          thus, the region from the muffin
c     tin to the circumscribed  sphere radii is divided  into  'npan'
c     pannels, each one containing 'nm(ipan)' points in order to take
c     care of the  discontinuities of the shape-function  derivative.
c     at the output one obtains :
c            llmsp (icell,ifun)       = integer array giving the com-
c                                       posite  index  lm=l*(l+1)+m+1
c                                       of the ifun-th shape function
c            lmsp  (icell,lm)         = (0,1)  if the lm-th component
c                                       is vanishing or not
c            nfu   (icell)            = number  of   shape   function
c                                       components in cell 'icell'
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD,NSPIND
      PARAMETER (NATYPD=38,NSPIND=2)
      INTEGER IRMD,IRNSD,LMAXD,LPOTD
      PARAMETER (irmd=1484,irnsd=508,lmaxd=4,lpotd=8)
      INTEGER NFUND,IRID
      PARAMETER (NFUND=289,irid=435)
      Integer    ndense
      Parameter ( NDENSE=1001 )
c
      Integer irmkd, irikd
      Parameter (irmkd=1484,irikd=435)
c
      INTEGER NCELLD,IPAND
      PARAMETER (NCELLD=20,IPAND=80)
      INTEGER NPOTD
      PARAMETER (NPOTD=NSPIND*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALAT,RFCTOR,VCONSTH,alpha
      INTEGER IFILE,INS,IPE,IPF,IPFE,KSHAPE,KVREL,KWS,LMAX,LPOT,NATREF,
     +     NEND,NSPIN,lsmear
C     ..  
C     .. Array Arguments ..
      REAL*8 ECORE(20,NPOTD),RMT(*),RMTNEW(*),
     +     THETAS(IRIKD,NFUND,*),VINS(IRMIND:IRMD,LMPOTD,*),
     $     thesme(irid,nfund,*),vspsme(irmd,*),
     $     theeqm(irid,nfund,*),vbc(2)
      INTEGER IFUNM(NATYPD,*),IMT(*),IPAN(*),IRCUT(0:IPAND,*),IRMIN(*),
     +     IRNS(*),ITITLE(20,*),LCORE(20,NPOTD),LLMSP(NATYPD,*),
     +     LMSP(NATYPD,*),NCORE(NPOTD),NFU(*),NTCELL(*)
      Integer irwseq(natypd), irmieq(natypd), ircueq(0:ipand,natypd),
     $     ipaneq(natypd)
C     ..
C     .. Scalars in Common ..
      REAL*8 C
      INTEGER LMAXP1
C     ..
C     .. Arrays in Common ..
      COMPLEX*16 MASS(IRMD,NATYPD)
      REAL*8 A(NATYPD),B(NATYPD),DRDI(IRMKD,NATYPD),
     +     R(IRMKD,NATYPD),RWSeq(NATYPD),RWSM1(NATYPD),
     $     S(0:LMAXD,NATYPD),VM2Z(IRMD,NPOTD),Z(NATYPD)
      INTEGER IRT(NATYPD),IRWS(NATYPD)
C     ..
c     Next are used for the equally spaced radial mesh in the
c     shape function region (used in the solution of the radial
c     equations
      Double Precision req(irmd,natypd), drdieq(irmd,natypd),
     $     droreq(irmd,natypd), rseq(irmd,0:lmaxd,natypd),
     $     dror(irmkd,natypd), rs(irmkd,0:lmaxd,natypd)
C     ..
C     .. Local Scalars ..
      REAL*8 A1,B1,EA,PI,S1,Z1, mu, alpha1
      INTEGER I,IA,ICELL,ICORE,IFUN,IH,IMT1,IMTM1,IPAN1,IR,IRC1,IRI,
     +     IRMINM,IRMINP,IRNS1P,IRT1P,IRWS1,ISAVE,ISPIN,ISUM,J,L,LM,
     +     LM1,LMPOT,LMPOTP,N,NCELL,NFUN,NR,idum1,idum2,idum3,idum4,
     $     nreq, ist, ialku, ien, npts, lask, ip, ipec
      Logical lvol, lvol1, lhost
      Double Precision reqeka, deltar, dist, mu1, mu2, alphamin,
     $     x, atim, value, volu, efermi, vbcin
C     ..
C     .. Local Arrays ..
      REAL*8 DRN(IRIKD,NCELLD),SCALE(NCELLD),U(IRMD),
     +     XRN(IRIKD,NCELLD),ZOLD(NATYPD), xrn1(irid),
     $     xdrn1(irid), rws(natypd), thet2d(irikd)
      INTEGER MESHN(NCELLD),NM(IPAND,NCELLD),NPAN(NCELLD),nmeq(2)
C     ..
C     .. External Subroutines ..
      EXTERNAL CALRMT,POTCUT, smear2, splint, spline, simpk
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ANINT,DATAN,EXP,LOG,MAX,MOD,REAL,SQRT,dsqrt
C     ..
C     .. Common blocks ..
      COMMON /DLOG87/VM2Z,RWS,RWSM1,DRDIEQ,REQ,Z,A,B,IRWSEQ,IRT,LMAXP1
      COMMON /MTSOL/MASS,C,S
C     ..
C     .. Save statement ..
      SAVE
C     ..
      PI = 4.D0*DATAN(1.D0)
      mu = 0.d0
      zold = 0.d0
c
c---> set speed of light
c
      C = 274.0720442D0
c
c---> read radial mesh information of the shape functions and
c     shape functions in the first iteration - if needed
c
      IF ((KSHAPE.NE.0) .AND. (IFILE.NE.0)) THEN
        READ (19,FMT=9000) NCELL
        WRITE (6,FMT=*) '  ncell : ',NCELL,NCELLD
        IF (NCELL.GT.NCELLD) STOP 'start'

        READ (19,FMT=9010) (SCALE(ICELL),ICELL=1,NCELL)
        DO 30 ICELL = 1,NCELL

c     Check, if the shape is used for the host atoms and set the
c     flag 'lhost'.
          lhost = .false.
          Do i = 1, natref
            If ( icell .eq. ntcell(i) ) lhost = .true.
          End Do

          READ (19,FMT=9000) NPAN(ICELL),MESHN(ICELL)
          IF (NPAN(ICELL)+1.GT.IPAND) THEN
            WRITE (6,FMT=*) 'INCREASE IPAND ',NPAN(ICELL) + 1,IPAND
            STOP
          END IF

          If ( irikd .lt. irid ) Then
            Write (6,*) 'The mesh with kinks should have more (or the',
     $           'same number) of points'
            Write (6,*) 'than the equally spaced mesh used in radial ',
     $           'equation'
            Write (6,*) ' irid = ',irid,'   irikd = ',irikd
            Stop
          End If
          If ( (lsmear .ge. 1) .and. (meshn(icell) .lt. irid) ) Then
            Write (6,*) 'START: icell, meshn, irid ',
     $           icell,meshn(icell),irid
            Stop
          End If

c     Next is working only, if ntcell_host < ntcell_cluster
c          If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
c     $         (icell .lt. ntcell(natref)) .and.
c     $         (meshn(icell) .gt. irid) ) Then
c     Next works (perhaps) in all cases
          If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
     $         (lhost) .and.
     $         (meshn(icell) .gt. irid) ) Then
            Write (6,*) 'START: icell, meshn, irid ',
     $           icell,meshn(icell),irid
            Stop
          End If
          If ( (lsmear .eq. 3 .or. lsmear .eq. 1)
     $         .and. ( .not. lhost ) .and.
     $         (meshn(icell) .gt. irikd) ) Then
            Write (6,*) 'START: icell, meshn, irikd ',
     $           icell,meshn(icell),irikd
            Stop
          End If
          If ( (lsmear .eq. 4 .or. lsmear .eq. 2) .and.
     $         (meshn(icell) .gt. irikd) ) Then
            Write (6,*) 'START: icell, meshn, irikd ',
     $           icell,meshn(icell),irikd
            Stop
          End If

          READ (19,FMT=9000) (NM(IPAN1,ICELL),IPAN1=2,NPAN(ICELL)+1)
          READ (19,FMT=9010) (XRN(IR,ICELL),DRN(IR,ICELL),IR=1,
     +      MESHN(ICELL))

          READ (19,FMT=9000) NFU(ICELL)
          NFUN = NFU(ICELL)
          WRITE (6,FMT=*) '  nfun  : ',NFUN,NFUND
          IF (NFUN.GT.NFUND) STOP 'start'

          DO 10 LM = 1,LMXSPD
            LMSP(ICELL,LM) = 0
   10     CONTINUE

          DO 20 IFUN = 1,NFUN
            READ (19,FMT=9000) LM
            LLMSP(ICELL,IFUN) = LM
            LMSP(ICELL,LM) = 1
            IFUNM(ICELL,LM) = IFUN
            READ (19,FMT=9010) (THETAS(N,IFUN,ICELL),N=1,MESHN(ICELL))
c     Make equally spaced shapes theeqm and smear these shapes
c     also to array thesme
            If ( lsmear .ge. 1 ) Then

              If ( (icell .eq. 1) .and. (ifun .eq. 1) ) Then

                Write (6,*) ' Choose smearing parameter alpha using',
     $               ' the first shape read in'
                alphamin = 2.5d0*( xrn(meshn(1),1) - xrn(1,1) )
     $               / dble(ndense)
c     Next with 135 points gives same alpha as 8.0 with 235 points
c     Next is used with 135 points
                alpha = 4.58d0*( xrn(meshn(1),1) - xrn(1,1) )
     $               / (irid-1)
c     Next is used, when only 45 points in the shapes
c                alpha = 3.58d0*( xrn(meshn(1),1) - xrn(1,1) )
c     $               / (irid-1)
                alpha = Max(alphamin,alpha)
              End If
c
c     Calculate the volume of the input shape (l=0 component)
              If ( ifun .eq. 1 ) Then
                Do ir = 1, meshn(icell)
                  thet2d(ir) = 4.d0*pi*thetas(ir,ifun,icell)*
     $                 xrn(ir,icell)**2
                End Do
                ircut(0,icell) = 0
                isum = 0
                Do ir = 1, npan(icell)
                  isum = isum + nm(ir+1,icell)
                  ircut(ir,icell) = isum
                End Do
              
                volu = 0.d0
                Call simpk( thet2d, volu, npan(icell), ircut(0,icell),
     $               drn(1,icell) )
                volu = volu + (4.d0*pi/3.d0)*Dsqrt(4.d0*pi)*
     $               xrn(1,icell)**3
                volu = volu/Dsqrt(4.d0*pi)
                Write (6,*) '  Volume of the input shape ',icell,
     $               ' is ',volu
              End If

c              If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
c     $             icell .eq. 1 ) Then
c     Next works, if ntcell_host < ntcell_cluster
c              If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
c     $             icell .le. ntcell(natref) ) Then

              If ( (lsmear .eq. 3 .or. lsmear .eq. 1)
     $             .and. lhost ) Then
                Do ir = 1, irid
                  theeqm(ir,ifun,icell) = thetas(ir,ifun,icell)
                End Do
                If ( lm .le. lmpotd ) Then
                  Do ir = 1, irid
                    thesme(ir,ifun,icell) = thetas(ir,ifun,icell)
                  End Do
                End If

              Else

                If ( ifun .eq. 1 ) Then
                  dist   = xrn(meshn(icell),icell) - xrn(1,icell)
                  deltar = dist / Dble( irid - 1 )
                  reqeka = xrn(1,icell)
                  Do ir = 1, irid
                    xrn1(ir)  = reqeka + (ir-1)*deltar
                    xdrn1(ir) = deltar
                  End Do
                  xrn1(irid) = xrn(meshn(icell),icell)
                  alpha1 = 2.5d0*dist / Dble(ndense)
                End If
                lvol1 = .false.
                If ( lm .eq. 1 ) lvol1 = .true.
                If ( lvol1 ) mu1 = 0.d0
                If ( lvol1 ) mu2 = 0.d0

c$$$                If ( .true. ) Then
                If ( lm .eq. 1 ) Then

                  Call smear2( thetas(1,ifun,icell),
     $                 theeqm(1,ifun,icell),xrn(1,icell), xrn1,
     $                 xdrn1, nm(1,icell),irikd, irid, meshn(icell),
     $                 npan(icell), alpha1, mu1, lvol1, volu )
c     Next is a direct interpolation from the 'kinky' shapes to
c     the an equally spaced shapes
                Else

                  ist   = 1 
                  ialku = 1
                  Do ip = 2, npan(icell) + 1
                    atim = 1.0d35
                    ien  = ist + nm(ip,icell) - 1
                    npts = nm(ip,icell)
                    If ( npts .ne. ien-ist+1 ) Stop 'ntps error'
                    Call spline( xrn(ist,icell),
     $                   thetas(ist,ifun,icell),
     $                   npts, a, a, thet2d(ist) )

                    lask = 0
                    Do ir = ialku, irid
                      x = xrn1(ir)
                      If ( x .gt. xrn(ien,icell) ) Goto 100
c                      If ( (x .eq. xrn(ien,icell) ) .and.
c     $                     (x .eq. xrn1(ir-1) ) ) Then
c                        Goto 100
c                      End If
                      lask = lask + 1
                      Call splint( xrn(ist,icell),
     $                     thetas(ist,ifun,icell),
     $                     thet2d(ist), npts, x, value )
                      theeqm(ir,ifun,icell) = value
                    End Do
 100                Continue
                    ialku = ialku + lask
                    ist   = ist + nm(ip,icell)
                  End Do
                End If

                nmeq(1)=0
                nmeq(2)=irid
                If ( lm .le. lmpotd ) Call smear2( theeqm(1,ifun,icell),
     $               thesme(1,ifun,icell), xrn1, xrn1, xdrn1,
     $               nmeq, irid, irid, irid,
     $               1 , alpha, mu2, lvol1, volu )

              End If

            Else
c     lsmear = 0 part: use only the input shapes
              Do ir = 1, irid
                theeqm(ir,ifun,icell) = thetas(ir,ifun,icell)
              End Do
              If ( lm .le. lmpotd ) Then
                Do ir = 1, irid
                  thesme(ir,ifun,icell) = thetas(ir,ifun,icell)
                End Do
              End If

c     lsmear .ge. 1 if ends
            End If

c     ifun loop ends
   20     CONTINUE

c     icell loop ends
   30   CONTINUE
c     
      END IF

      LMPOT = (LPOT+1)* (LPOT+1)
      LMAXP1 = LMAX + 1
      DO 150 IH = 1,NEND
        DO 140 ISPIN = 1,NSPIN
          I = NSPIN* (IH-1) + ISPIN

          IF (IFILE.NE.0) THEN
            IRCUT(0,IH) = 0
            ircueq(0,ih) = 0
            IF (KSHAPE.NE.0) THEN
              ICELL      = NTCELL(IH)

              IPAN(IH)   = 1 + NPAN(ICELL)
              ipaneq(ih) = 2
              If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
     $             ih .le. natref ) ipaneq(ih) = ipan(ih)
              If ( lsmear .eq. 0 ) ipaneq(ih) = ipan(ih)

            ELSE

              IPAN(IH) = 1
              ipaneq(ih) = 1
            END IF
c
c---> read title of potential card
c
            READ (IFILE,FMT=9020) (ITITLE(IA,I),IA=1,20)
            WRITE (6,FMT=9070) (ITITLE(IA,I),IA=1,20)

c
c---> read muffin-tin radius , lattice constant and new muffin radius
c      (new mt radius is adapted to the given radial mesh)
c
            READ (IFILE,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
c$$$            write (6,FMT=9030) RMT(IH),ALAT,RMTNEW(IH)
c
c---> read nuclear charge , maybe a shift , lmax of the core states ,
c     wigner seitz radius , last mesh point of the exponential mesh
c     (in case of ws input-potential: last mesh point corresponds
c     to ws-radius, in case of shape-corrected input-potential
c     last mesh point of the exponential mesh corresponds to
c     mt-radius/nevertheless this point is always in the array
c     irws(ih)),number of points for the radial non-muffin-tin
c     mesh  needed for shape functions, the constants a and b
c     for the radial exponential mesh : r(i) = b*(exp(a*(i-1))-1)
c     the no. of different core states and some other stuff
c
            If ( (ih .eq. 1) .and. (ispin .eq. 1) ) Then
              READ (IFILE,FMT=9041) Z(IH),ZOLD(IH),RWS(IH),
     $             efermi,vbcin,
     $             IRWSeq(IH),A(IH),B(IH),NCORE(I)
              If ( dabs(vbcin - vbc(1)) .gt. 1.d-6 ) Then
                Write (6,*) '*** Warning: different VBC in the input',
     $               ' card and in the potential'
                Write (6,*) '*** Using VBC from the input card ',vbc(1)
              Else
                Write (6,*) '*** Using VBC from the potential ',vbcin
                vbc(1) = vbcin
                vbc(2) = vbcin
              End If
            Else
              READ (IFILE,FMT=9040) Z(IH),ZOLD(IH),RWS(IH),IRWSeq(IH),
     +             A(IH),B(IH),NCORE(I)
            End If
c$$$            write (6,FMT=9040) Z(IH),ZOLD(IH),RWS(IH),IRWSeq(IH),
c$$$     +        A(IH),B(IH),NCORE(I)
c
c     ASA case:
            If ( kshape .eq. 0 ) Then

              irws(ih) = irwseq(ih)
            Else

c              irws(ih)   = (irmd-irid) + meshn(icell)
              irws(ih) = irwseq(ih) - irid + meshn(icell)
              If ( lsmear .eq. 0 ) irws(ih) = irwseq(ih)
              If ( (lsmear .eq. 3 .or. lsmear .eq. 1) .and.
     $             ih .le. natref ) irws(ih) = irwseq(ih)

            End If

            NR = IRWS(IH)
            NReq = IRWSeq(IH)
c
c---> read the different core states : l and energy
c
            IF (NCORE(I).GE.1) READ (IFILE,FMT=9060) (LCORE(ICORE,I),
     +          ECORE(ICORE,I),ICORE=1,NCORE(I))
c$$$            IF (NCORE(I).GE.1) write (6,FMT=9060) (LCORE(ICORE,I),
c$$$     +          ECORE(ICORE,I),ICORE=1,NCORE(I))
c
            IF ((INS.LE.1.OR.IH.LE.NATREF) .AND. INS.LE.2) THEN
c
c--->read radial mesh points, its derivative, the spherically averaged
c      charge density and the input potential without the nuclear pot.
c
              READ (IFILE,FMT=9050) (R(IR,IH),DRDI(IR,IH),VM2Z(IR,I),
     +          IR=1,NReq)

              Do ir = 1,nreq
                vspsme(ir,i) = vm2z(ir,i)
              End Do
            ELSE
c
c---> read full potential - the non spherical contribution from irmin
c      to irt - remember that the lm = 1 contribution is multiplied by
c      1/sqrt(4 pi)
c
              READ (IFILE,FMT=9080) IRT1P,IRNS1P,LMPOTP,ISAVE
              If ( irns1p .ne. irnsd )Then
                Write (6,*) 'ERROR in start: irns1p ne irnsd ',
     $               irns1p,irid
                Stop
              End If

              IRMINP = IRT1P - IRNS1P
              IRMINM = MAX(IRMINP,IRMIND)
              READ (IFILE,FMT=9090) (VM2Z(IR,I),IR=1,NReq)
c
c     lsmear=4: host+clu eq.mesh, smepot is read in
              If ( lsmear .eq. 4 ) Then
                Read (ifile,fmt=9081) idum1,idum2,idum3,idum4,alpha1
                Read (ifile,fmt=9090) (vspsme(ir,i),ir=1,nreq)
                Write (6,*) '*** Smearing read in: alpha is ',alpha1
              End If
c
c     lsmear=3: host with kinks, clu eq.mesh, smepot is read in
              If ( lsmear .eq. 3 ) Then
                If ( ih .gt. natref ) Then
                  Read (ifile,fmt=9081) idum1,idum2,idum3,idum4,alpha1
                  Read (ifile,fmt=9090) (vspsme(ir,i),ir=1,nreq)
                  Write (6,*) '*** Smearing read in: alpha is ',alpha1
                Else
                  Do ir = 1, nreq
                    vspsme(ir,i) = vm2z(ir,i)
                  End Do
                End If
              End If
c
c     lsmear=2: host+clu eq.mesh, smepot is read in for host
              If ( lsmear .eq. 2 ) Then
                If ( ih .le. natref ) Then
                  Read (ifile,fmt=9081) idum1,idum2,idum3,idum4,alpha1
                  Read (ifile,fmt=9090) (vspsme(ir,i),ir=1,nreq)
                  Write (6,*) '*** Smearing read in: alpha is ',alpha1
                Else
                  Do ir = 1, nreq
                    vspsme(ir,i) = vm2z(ir,i)
                  End Do
                End If
              End If
c
c     lsmear=1: host with kinks, clu eq.mesh, smepot is not read in
c     lsmear=0: host+clu with kinks, smepot is not read in
              If ( lsmear .le. 1 ) Then
                  Do ir = 1, nreq
                    vspsme(ir,i) = vm2z(ir,i)
                  End Do
              End If
              If ( lsmear .gt. 4 ) Stop 'START: lsmear '
c
              IF (LMPOTP.GT.1) THEN
                LM1 = 2
                DO 50 LM = 2,LMPOTP
                  IF (LM1.NE.1) THEN
                    IF (ISAVE.EQ.1) THEN
                      READ (IFILE,FMT=9080) LM1

                    ELSE
                      LM1 = LM
                    END IF

                    IF (LM1.GT.1) THEN

                      READ (IFILE,FMT=9090) (U(IR),IR=IRMINP,NReq)

                      IF (LM1.LE.LMPOT) THEN
                        DO 40 IR = IRMINM,NReq
                          VINS(IR,LM1,I) = U(IR)
   40                   CONTINUE
                      END IF

                    END IF

                  END IF

   50           CONTINUE

              END IF

            END IF
c
c---> conversion factor for energy from p.u. to ry
c                      (only for green's functions)
c
            RFCTOR = 0.5D0*ALAT/PI
            IRWS1 = IRWS(IH)
c
c---> redefine new mt-radius in case of shape corrections
c
            IF (KSHAPE.NE.0) THEN
              RMTNEW(IH) = SCALE(ICELL)*ALAT*XRN(1,ICELL)
              IMTM1 = ANINT(LOG(RMTNEW(IH)/B(IH)+1.0D0)/A(IH))
              IMT1 = IMTM1 + 1
c
c---> for proper core treatment imt must be odd
c     shift potential by one mesh point if imt is even
c
              IF (MOD(IMT1,2).EQ.0) THEN
                IMT1 = IMT1 + 1
                DO 60 IR = IMT1,2,-1
                  VM2Z(IR,I) = VM2Z(IR-1,I)
                  vspsme(ir,i) = vspsme(ir-1,i)
   60           CONTINUE
              END IF
c
              IMT(IH) = IMT1
              B(IH) = RMTNEW(IH)/ (EXP(A(IH)*REAL(IMT1-1))-1.0D0)

            Else
              b(ih) = rws(ih)/ (exp(a(ih)*real(irws(ih)-1))-1.0d0)
            END IF
c
c---> generate radial mesh - the format statements are inaccurate
c
            A1 = A(IH)
            B1 = B(IH)
            R(1,IH) = 0.0D0
            DRDI(1,IH) = A1*B1
            req(1,IH)    = 0.0D0
            drdieq(1,IH) = A1*B1
            DO 70 IR = 2,IRWS1
              EA = EXP(A1*REAL(IR-1))
              R(IR,IH) = B1* (EA-1.0D0)
              DRDI(IR,IH) = A1*B1*EA
              DROR(IR,IH) = A1/ (1.0D0-1.0D0/EA)
   70       CONTINUE
c$$$            Do ir = 2, irwseq(ih)
c$$$              ea = exp(a1*real(ir-1))
c$$$              droreq(ir,ih) = a1/ (1.0d0-1.0d0/ea)
c$$$            End Do

            If ( kshape .eq. 0 ) Then
              Do ir = 1, irws(ih)
                req(ir,ih)    =    r(ir,ih)
                drdieq(ir,ih) = drdi(ir,ih)
                droreq(ir,ih) = dror(ir,ih)
              End Do
            Else
              Do ir = 1, imt1
                req(ir,ih)    =    r(ir,ih)
                drdieq(ir,ih) = drdi(ir,ih)
                droreq(ir,ih) = dror(ir,ih)
              End Do
            End If
c
c---> fill cell-type depending mesh points in the non-muffin-tin-region
c
            IF (KSHAPE.NE.0) THEN


c     smeared potential inside the MT sphere should be equal
c     to the non smeared potential:
              If ( lsmear .ge. 2 ) Then
                Do ir = 1, imt1
                  vspsme(ir,i) = vm2z(ir,i)
                End Do
              End If

               Do IRI = 1,MESHN(ICELL)
                IR = IRI + IMT1
                R(IR,IH) = SCALE(ICELL)*ALAT*XRN(IRI,ICELL)
                DRDI(IR,IH) = SCALE(ICELL)*ALAT*DRN(IRI,ICELL)
                DROR(IR,IH) = DRDI(IR,IH)/R(IR,IH)
              End Do
c     Make an equally spaced mesh in the shape reqion. This mesh
c     is used in the solution of the radial equations.
              dist   = scale(icell)*alat*
     $             ( xrn(meshn(icell),icell) - xrn(1,icell) )
              deltar = dist / ( irid - 1 )
              reqeka = scale(icell)*alat*xrn(1,icell)
              Do iri = 1, irid
                ir = iri + imt1
                req(ir,ih)    = reqeka + (iri-1)*deltar
                drdieq(ir,ih) = deltar
                droreq(ir,ih) = deltar/req(ir,ih)
              End Do
              If ( Dabs( req(irmd,ih) - r(meshn(icell)+imt1,ih) )
     $             .gt. 1.0d-10 ) Then
                Write (6,*) 'STARTB: Error, radial mesh is not ok'
                Write (6,*) 'irmd, meshn ', irmd, meshn(icell)
                Write (6,*) req(irmd,ih), r(meshn(icell)+imt1,ih)
                Stop 'STARTB'
              End If
              req(irmd,ih) = r(meshn(icell)+imt1,ih)
c
              If ( (lsmear .eq. 1 .or. lsmear .eq. 3) .and.
     $             ih .le. natref) Then
                If ( irid .ne. meshn(icell) ) Stop 'START: 21'
                Do iri = 1, irid
                  ir = iri + imt1
                  req(ir,ih)    = r(ir,ih)
                  drdieq(ir,ih) = drdi(ir,ih)
                  droreq(ir,ih) = drdieq(ir,ih)/req(ir,ih)
                End Do
              End If
c
              If ( lsmear .eq. 0 ) Then
                If ( irid .ne. meshn(icell) ) Stop 'START: 22'
                Do iri = 1, irid
                  ir = iri + imt1
                  req(ir,ih)    = r(ir,ih)
                  drdieq(ir,ih) = drdi(ir,ih)
                  droreq(ir,ih) = drdieq(ir,ih)/req(ir,ih)
                End Do
              End If
c
c     kshape .ne. 0 if ends
            End If
c
            RWSM1(IH) = R(IRWS1-1,IH)
            RWS(IH)   = R(IRWS1,IH)
            RWSeq(IH) = Req(IRWSeq(ih),IH)

c
c---> kshape.eq.0 : calculate new rmt adapted to exp. mesh
c
            ipec = mod(ipe,2)

            CALL CALRMT(IPF,IPFE,IPEc,IMT(IH),Z(IH),RMT(IH),RWS(IH),
     +           RMTNEW(IH),ALAT,DRDI(1,IH),A(IH),B(IH),
     $           IRWS(IH),R(1,IH),IFILE,KSHAPE)

            CALL CALRMT(IPF,IPFE,IPEc,IMT(IH),Z(IH),RMT(IH),RWSeq(IH),
     +           RMTNEW(IH),ALAT,DRDIeq(1,IH),A(IH),B(IH),
     $           IRWSeq(ih),Req(1,IH),IFILE,KSHAPE)
c
            IF (KSHAPE.GT.0) THEN
              IRCUT(1,IH) = IMT(IH)
              ircueq(1,ih) = imt(ih)
              ircueq(2,ih) = irmd
              ISUM = IMT(IH)
              DO 90 IPAN1 = 2,IPAN(IH)
                ISUM = ISUM + NM(IPAN1,ICELL)
                IRCUT(IPAN1,IH) = ISUM
                If ( (lsmear .eq. 1 .or. lsmear .eq. 3) .and.
     $               ih .le. natref) ircueq(ipan1,ih) = isum
                If ( lsmear .eq. 0 ) ircueq(ipan1,ih) = isum
   90         CONTINUE
              NR = ISUM
              NReq = irwseq(ih)

            ELSE

              NR   = IRWS(IH)
              nreq = irwseq(ih)
              IF ((KWS.EQ.1.AND.IH.GT.NATREF) .OR. KWS.EQ.2) THEN
                IRCUT(1,IH) = IRWS1
                ircueq(1,ih) = irmd

              ELSE
                IRCUT(1,IH) = IMT(IH)
                ircueq(1,ih) = imt(ih)
              END IF

            END IF
c
c---> fill array irmin in case of full potential
c
            IF (INS.NE.0) Then
              IRMIN(IH) = NR - IRNS(IH)
              irmieq(ih) = irmind
              If ( (lsmear .eq. 1 .or. lsmear .eq. 3) .and.
     $             ih .le. natref ) irmieq(ih) = irmin(ih)
              If ( lsmear .eq. 0 )  irmieq(ih) = irmin(ih)
            End If
c
c---> generate arrays for the calculation of the wave functions
c
            Z1 = Z(IH)
            DO 110 L = 0,LMAX
              IF ((KVREL.EQ.1.AND.IH.GT.NATREF) .OR. KVREL.EQ.2) THEN
                S1 = SQRT(REAL(L*L+L+1)-4.0D0*Z1*Z1/ (C*C))
                IF (Z1.EQ.0.0D0) S1 = REAL(L)

              ELSE

                S1 = REAL(L)
              END IF

              S(L,IH) = S1
              rseq(1,L,IH) = 0.0D0
              RS(1,L,IH) = 0.0D0
              Do ir = 2, irmd
                RS(IR,L,IH) = R(IR,IH)**S1
                rseq(ir,l,ih) = req(ir,ih)**s1
              End Do
  110       CONTINUE
c
c---> cut input potential at rmt if given only at exponential mesh
c
            IF (KSHAPE.EQ.1) THEN
              Stop 'START: no pot cutting in this version (yet)'
              IMT1 = IMT(IH)
              IRC1 = IRCUT(IPAN(IH),IH)
              CALL POTCUT(IMT1,IRC1,INS,LMPOT,R(1,IH),VM2Z(1,I),
     +             vspsme(1,i),VINS(IRMIND,1,I),Z(IH))
            END IF
c
c--->  first iteration : shift all potentials (only for test purpose)
c
            DO 120 J = 1,NReq
              VM2Z(J,I) = VM2Z(J,I) + VCONSTH
              vspsme(j,i) = vspsme(j,i) + vconsth
  120       CONTINUE
c
c--->  convert potential from zold to z (generate impurity with z)
c
            IF (ZOLD(IH).NE.0.0D0) THEN
              Stop 'zold not in this version'
c$$$              DO 130 J = 2,NReq
c$$$                VM2Z(J,I) = VM2Z(J,I)+2.0D0*(Z(IH)-ZOLD(IH))/Req(J,IH)
c$$$                vspsme(j,i) = vspsme(j,i) +
c$$$     $               2.0d0* (z(ih)-zold(ih))/req(j,ih)
c$$$  130         CONTINUE
            END IF

          END IF

          IF (KSHAPE.EQ.0 .AND. ((KWS.EQ.1.AND.IH.LE.NATREF).OR.
     +        KWS.EQ.0)) THEN
c
c---> in case of a mt calculation cut potential at mt radius
c
            IMT1 = IMT(IH)
            IRWS1 = IRWS(IH)
            CALL POTCUT(IMT1,IRWS1,INS,LMPOT,R(1,IH),VM2Z(1,I),
     +           vspsme(1,i),VINS(IRMIND,1,I),Z(IH))
          END IF

  140   CONTINUE

  150 CONTINUE

 9000 FORMAT (16i5)
 9010 FORMAT (4D20.12)
 9020 FORMAT (20a4)
 9030 FORMAT (3f12.8)
 9040 FORMAT (2F10.5,/,F10.5,/,I5,/,2D15.8,/,I2)
 9041 FORMAT (2F10.5,/,F10.5,2f15.10/,I5,/,2D15.8,/,I2)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9060 FORMAT (i5,1p,d15.6)
 9070 FORMAT (10x,20a4)
 9080 FORMAT (10i5)
 9081 FORMAT (4i5,7x,1f10.7)
 9090 FORMAT (1p,4D20.13)
      END
