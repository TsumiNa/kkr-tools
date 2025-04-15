      Subroutine mtzsme(lmpot,natyp,nspin,v,r,drdi,imt,ircut,ipan,
     +     ntcell,lmsp,irws,vspsmo,lsmear,nstart)
      Implicit None
c-----------------------------------------------------------------------
c
c     determine muffin tin zero and shift potential to muffin tin zero
c
c     for spin polarized calculations muffin tin zero is related to
c         spin up (majority spin direction)
c
c                                            may,1990  b. drittler
c
c     This routine calculates the difference between the mtz of the
c     unsmeared radial potential and the smeared one. It then 
c     iterates the correct mtz for the smeared one, so finally
c     both mtz's are the same:
c
c     Delta VMTZ =  { (V_o(r) - V~_o(r))*Y_0*4pi*r^2 dr = 
c                  int.
c           
c                      o      o
c                = rfpi* { (V(r) - V~(r))*r^2 dr = 0
c
c                                          June, 1997  T. Korhonen
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      Integer natypd
      Parameter (NATYPD=38)
      Integer irmd,lpotd
      Parameter (irmd=1484,lpotd=8)
      Integer nfund,irid
      Parameter (NFUND=289,irid=435)
      Integer ipand
      Parameter (IPAND=80)
      Integer lmpotd
      Parameter (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      Double Precision eshift,vbc
      Integer ishift,lmpot,natyp,nspin,lsmear,nstart
C     ..
C     .. Array Arguments ..
      Double Precision drdi(irmd,*),r(irmd,*),v(irmd,lmpotd,*),
     $     vspsmo(irmd,*)
      Integer imt(*),ipan(*),ircut(0:ipand,*),irws(*),
     +     lmsp(natypd,*),ntcell(*)
C     ..
C     .. Local Scalars ..
      Double Precision fpi,rfpi,x,alpha,ral,rlo,deltav,fral,frlo,fnol,
     $     vints,vint1,vint2
      Integer icell,ih,imt1,ipan1,ipot,ir,irc1,it,is
C     ..
C     .. Local Arrays ..
      Double Precision v1(irmd),v2(irmd),vs(irmd)
C     ..
C     .. External Subroutines ..
      External simpk
C     ..
C     .. Intrinsic Functions ..
      Intrinsic Datan, Sqrt
C     ..
C     .. Local Statement Function ..
      Double Precision fun1
C     Fermi function (actually its derivative)
      fun1(x,alpha) = (1.d0/alpha)*exp(-x/alpha)/(exp(-x/alpha)+1)**2
C     ..
C     ..
      fpi   = 16.0d0*Datan(1.0d0)
      rfpi  = Sqrt(fpi)

      If ( lsmear .lt. 1 ) Return


c     V(i,1,*) = sqrt(4pi)*V_spherical
c---  > full potential calculation with smearing of rad. potential

      Do is = 1, nspin
        Do ih = nstart, natyp

          ipot  = nspin*(ih-1) + is
          ipan1 = ipan(ih)
          imt1  = imt(ih)
          irc1  = ircut(ipan1,ih)
          icell = ntcell(ih)

          Do ir = 1, irmd
            v1(ir) = 0.0d0
            v2(ir) = 0.0d0
            vs(ir) = 0.0d0
          End Do
          alpha = 5.d0*( r(irc1,ih) - r(imt1+1,ih) )/ (irc1 - imt1)
          fnol  = fun1(0.d0,alpha)
          Do ir = imt1 + 1, irc1
            v2(ir) = rfpi*r(ir,ih)**2
            v1(ir) = (v(ir,1,ipot)/rfpi)*v2(ir)
          End Do
          vint1 = 0.d0
          vint2 = 0.d0
          Call simpk( v2, vint2, ipan1, ircut(0,ih), drdi(1,ih) )
          Call simpk( v1, vint1, ipan1, ircut(0,ih), drdi(1,ih) )
c          Write (6,*) 'MTZSME: ispin,vint1 vint2 ',is,vint1,vint2

          deltav = 0.d0
          Do ir = imt1 + 1, irc1
            vs(ir) = vspsmo(ir,ipot)
          End Do
          Do it = 1, 20
            vints = 0.d0
            Do ir = imt1 + 1, irc1
              vs(ir) = vs(ir)*rfpi*r(ir,ih)**2
            End Do
            Call simpk( vs, vints, ipan1, ircut(0,ih), drdi(1,ih) )
            deltav = deltav + (vint1 - vints)/vint2
c            Write (6,*) 'ispin,vints deltav ',is,vints,deltav
            Do ir = imt1 + 1, irc1
              ral = r(ir,ih) - r(imt1+1,ih)
              rlo = r(irc1,ih) - r(ir,ih)
              fral  = fun1(ral,alpha)/fnol
              frlo  = fun1(rlo,alpha)/fnol
              vs(ir) = vspsmo(ir,ipot) +
     $             1.15d0*deltav
              vs(ir) = ( v(ir,1,ipot)*fral/rfpi +
     $             (1.d0 - fral)*vs(ir) )
              vs(ir) = ( v(ir,1,ipot)*frlo/rfpi +
     $             (1.d0 - frlo)*vs(ir) ) 
            End Do
            
          End Do
          Do ir = imt1 + 1, irc1
            vspsmo(ir,ipot) = vs(ir)
          End Do
          
        End Do

      End Do

      Return
      End
