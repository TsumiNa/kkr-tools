      Subroutine smear2( thetas, thesme, xrn, xrn1, xdrn1, nm, irikd,
     $     irid, meshn, npan, alpha, mu, lvol, volu )
      Implicit None
c     ******************************************************************
c     * Interpolates the input shape to a dense mesh and after
c     * that does the smearing. 
c     * The radial meshes of the input and smeared shapes should
c     * be equal (the other part of the program needs this, not
c     * this subroutine, if it is run as a separate program).
c     * The radial mesh of smeared shapes is adsumed to be
c     * equally spaced with IRID points.
c     *                                 5.6.1997 T. Korhonen
c     ******************************************************************
c     Passed variables:
c     irid    : dimension of the shapes array
c     meshn   : number of points in the old shapes 
c     npan    : number of panels in the old shapes
c     thetas  : old shapes
c     thesme  : smeared shapes (output)
c     xrn     : radial mesh of the old shapes
c     alpha   : smearing parameter used if the Fermi function
c     mu      : shift of the Fermi function (chemical pot)
c               used to get the correct l=0 volume 
c     lvol    : iterate the correct volume 
c     volu    : volume of the WS cell (bcc=4,fcc=2,...)
      Integer irikd, meshn, npan, nm(npan+1), irid
      Double Precision thetas(irikd), thesme(irid), xrn(irikd), mu,
     $     alpha, xrn1(irid), xdrn1(irid), volu
      Logical lvol
c     
c     ndense:  the number of points in the mesh where
c     input shapes are interpolated before smearing
c     number of points inside and outside of the range
c     of shape function radial mesh. Needed in order
c     to calculate smearing integral correctly.
      Integer    ndense
      Parameter ( NDENSE=1001 )

c     Local arrays
      Double Precision v1(ndense), rdense(ndense+1)
c     Local scalars
      Double Precision drdens, r, r01, w, sum, sum1,
     $     rfpi, pi, fpi, dist, x, fnol, fral, frlo,
     $     ral, rlo, tim1
      Integer j, itno, nhere, ipan, np, itmax, iold, iden, i,
     $     nin

c     Local statement function
      Double Precision fun5
      External fun5
c     Fermi function, better localized (used !!)
c      fun5(x,alpha) = (1.d0/alpha)*exp(-x/alpha)/(exp(-x/alpha)+1)**2
c     
      pi       = 4.d0*Datan( 1.d0 )
      fpi      = 4.d0*pi
      rfpi     = Dsqrt( 4.d0*pi )
c      mu       = 0.d0
c      write(6,*) '**** mu is ',mu
      itmax    = 1
      tim1 = 0.d0
      If (lvol) tim1  = rfpi
      If (lvol) itmax = 10
      nin    = ndense/6
      If ( Mod( nin,2 )    .eq. 0 ) nin    = nin   - 1

      If ( alpha .lt. 0.d0 ) Then
        Write (6,*) '** SMEAR2: alpha lt 0 ',alpha
        Stop
      End If

c     The new mesh for which the old shapes are interpolated before
c     smearing
      dist   =  xrn(meshn) - xrn(1)
      drdens = dist / Dble(ndense-1)
      Do i = 1, ndense
        rdense(i) = xrn(1) + drdens*(i-1)
      End Do
      rdense(ndense) = xrn(meshn)
c     Next is used as a terminator flag
      rdense(ndense+1) = 9.9d12

c     Interpolate the old shapes to a dense mesh
c     Next works every original shape mesh, in the original
c     shapes the boundary point between two panels is given
c     twice.

      v1(1) = thetas(1)
      iold      = 1
      iden      = 1

      Do ipan = 1, npan
        np    = nm(ipan+1)
        nhere = 0
 10     Continue
        If ( rdense(iden+nhere) .le. xrn(iold+np-1) ) Then
          nhere = nhere + 1
c          Write(6,*) 'ipan, nhere ',ipan,nhere
          Goto 10
        End If
        iold = iold + 1

        If ( np .gt. nhere ) Then
          Do j = 1, np - 1
            v1(iden) = ( thetas(iold-1)*(xrn(iold)-rdense(iden)) +
     $           thetas(iold)*(rdense(iden)-xrn(iold-1)) ) /
     $           ( xrn(iold) - xrn(iold-1) )
            If ( xrn(iold) .ge. rdense(iden) ) iden = iden + 1
            iold = iold + 1
          End Do
        Else
          Do j = 1, nhere
            v1(iden) = ( thetas(iold-1)*(xrn(iold)-rdense(iden)) +
     $           thetas(iold)*(rdense(iden)-xrn(iold-1)) ) /
     $           ( xrn(iold) - xrn(iold-1) )
            iden = iden + 1
            If ( xrn(iold) .lt. rdense(iden) ) iold = iold + 1
          End Do
        End If
c        Write(6,*) 'ipan, iold, iden ',ipan,iold,iden
      End Do

      v1(1)      = thetas(1)
      v1(ndense) = thetas(meshn)
      If ( (iden-1) .ne. ndense .or. (iold-1) .ne. meshn ) Then
        Write (6,*) 'SMEAR2: *** iden, iold, ndense, meshn',
     $       iden, iold, ndense, meshn
        Stop
      End If


      fnol   = fun5(0.d0,alpha)
 123  Continue
      If ( fun5(rdense(nin)-rdense(1),alpha)/fnol .gt. 1.e-4 ) Then
c     Next with 45 points
c      If ( fun5(rdense(nin)-rdense(1),alpha)/fnol .gt. 1.e-3 ) Then
        nin = nin + 2
        If ( nin .ge. ndense ) Stop 'SMEAR: nin too big'
        Goto 123
      End If

c     Iterate the volume (if lvol is true)
      Do itno = 1,itmax

c     Loop over new radial mesh points xrn1
        Do i = 1, irid
          sum  = 0.0d0
          sum1 = 0.0d0

          Do j = 1, ndense
            w = 4.d0/3.d0
            If ( Mod(j,2) .ne. 0 )             w = 2.d0/3.d0
            If ( j .eq. 1 .or. j .eq. ndense ) w = 1.d0/3.d0
            r    = rdense(j) - mu - xrn1(i)
            sum  = sum  + w*drdens*fun5(r,alpha)*v1(j)
            sum1 = sum1 + w*drdens*fun5(r,alpha)
          End Do

          Do j = 1, nin
            w = 4.d0/3.d0
            If ( Mod(j,2) .ne. 0 )          w = 2.d0/3.d0
            If ( j .eq. 1 .or. j .eq. nin ) w = 1.d0/3.d0
c     Next points inside the sphere
            r01  = xrn(1) - (j-1)*drdens
            r    = r01 - mu - xrn1(i)
            sum  = sum  + w*drdens*fun5(r,alpha)*v1(1)
            sum1 = sum1 + w*drdens*fun5(r,alpha)
c     Next points outside the sphere
            r01  = xrn(meshn) + (j-1)*drdens
            r    = r01 - mu - xrn1(i)
            sum  = sum  + w*drdens*fun5(r,alpha)*0.d0
            sum1 = sum1 + w*drdens*fun5(r,alpha)
          End Do

          thesme(i) = sum/sum1

          ral   = xrn1(i) - xrn(1)
          rlo   = xrn(meshn) - xrn1(i)
          fral  = fun5(ral,alpha)/fnol
          frlo  = fun5(rlo,alpha)/fnol
          thesme(i) = ( tim1*fral + (1.d0 - fral)*thesme(i) )
          thesme(i) = ( 0.d0*frlo + (1.d0 - frlo)*thesme(i) ) 

        End Do

c     Integrate the volume
        If (lvol) Then
          sum  = 0.0d0
          sum1 = 0.0d0
          Do i = 1, irid
            w    = 4.d0/3.d0
            If ( Mod(i,2) .ne. 0 )            w = 2.d0/3.d0
            If ( i .eq. 1 .or. i .eq. irid ) w = 1.d0/3.d0
            sum1 = sum1 + w*4.d0*pi*thesme(i)*xdrn1(i)*xrn1(i)**2
          End Do
c     Add the volume of the 'inner' sphere
          sum1 = sum1 + (4.d0*pi/3.d0)*rfpi*xrn1(1)**3
          sum1 = sum1 / Dsqrt(4.d0*pi)
c
c     Next line is not used, if one uses the Voronoi construction
c     for the shapes. For 'old' shapes, it was only a test case,
c     so one could comment it out.
c$$$          If ( Dabs(sum1 - volu) .gt. 0.1d0 ) Then
c$$$            Write (6,*) 'Error: Volume should be close to ',volu
c$$$            Stop ' in SMEAR '
c$$$          End If
c     Shift the 'chemical potential' of the Fermi function
          mu = mu + ( 3.d0*sum1 / fpi )**(1.d0/3.d0) -
     $         ( 3.d0*volu / fpi )**(1.d0/3.d0)
        End If

c     End itno loop
      End Do

      Return
      End
