      SUBROUTINE SEMICORE(NATYP,LMAX,SUMDF,IE,npntsc,npntim,
     +                    ISPIN,IPAN,IRCUT,NTCELL,rho2ns,NFU,
     +                    llmsp,THETAS,EZ,DRDI,ielast,DF,DEN,
     +                    factor,espco2)

      implicit none
      

      include 'inc.fi'

      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
   
      INTEGER I1,LM1,IE,npntsc,npntim,ispin,ir,ipan1,irs1,irc1,
     +        ifun,NATYP,LMAX,icell,ielast
      INTEGER IPAN(NATYPD),IRCUT(0:IPAND,NATYPD),NTCELL(NATYPD),
     +        NFU(NATYPD),LLMSP(NATYPD,NFUND)
      DOUBLE PRECISION PI,FPI,RFPI,espco2,sum,sumtot,sumint,
     +                 factor
      DOUBLE PRECISION SUMDF(natypd),sum1(natypd),rho(irmd),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,NSPIND),
     +                 THETAS(IRID,NFUND,NCELLD),DRDI(IRMD,NATYPD)
      DOUBLE COMPLEX DEN(IEMXD,0:LMAXD,NATYPD,NSPIND),EZ(IEMXD),DF


      If ( ie .eq. 1 ) sumdf  = 0.d0 
      If ( ie .le.  (npntsc + 2*npntim) ) Then
            Do i1 = 1, natyp
              Do lm1 = 0, lmax
                sumdf(i1) = sumdf(i1) + Dimag(den(ie,lm1,i1,ispin)*df)
              End Do
            End Do
      End If

c     Integrate the semicore charge
c     Inside the MT-sphere only l=0 component has charge, the
c     higher l's only redistribute the charge (i.e., l=0 component
c     is the spherical average). 
c     Outside the MT-sphere one should use the shape functions.

      PI = 4.0D0*DATAN(1.0D0)
      FPI = 4.0D0*PI
      RFPI = DSQRT(FPI)
      
      If ( ie .eq. (npntsc + 2*npntim) ) Then
      If ( ispin .eq. 1 ) espco2 = 0.d0
            sum    = 0.d0
            sumtot = 0.d0
            Do i1 = 1, natyp
              sum1(i1) = 0.d0
              ipan1 = ipan(i1)
              irs1  = ircut(1,i1)
              irc1  = ircut(ipan1,i1)
              icell = ntcell(i1)
              Do ir = 1, irs1
                rho(ir) = rho2ns(ir,1,i1,ispin)*rfpi
              End Do
              Do ir = irs1 + 1, irc1
                rho(ir) = 0.d0
              End Do
              Do ifun = 1, nfu(icell)
                lm1 = llmsp(icell,ifun)
                If ( lm1 .le. lmpotd ) Then
                  Do ir = irs1 + 1, irc1
                    rho(ir) = rho(ir) + rho2ns(ir,lm1,i1,ispin)*
     +                   thetas(ir-irs1,ifun,icell)
                  End Do
                End If
              End Do
              Call simpk(rho,sum1(i1),ipan1,ircut(0,i1),drdi(1,i1))
              sum    = sum    + sum1(i1)
              sumtot = sumtot + sumdf(i1)
              
              Write (6,FMT='(a,i4,a,f16.10,a,f16.10,a)')
     $             ' atom ',i1,' semicore charge = ', sum1(i1),
     $             ' (rho) ', sumdf(i1),' (dos)'
c     Loop over atoms
            End Do
            
            sumint = Anint( sum )
            Write (6,FMT='(a,f16.10)') ' total semicore charge = ',sum
            Write (6,FMT='(a,f16.10,a)') ' total semicore charge = ',
     $           sumint, ' used for semicore s.p. energies'

            If ( sumint .gt. 0.1d0 ) Then
              factor = sumint/sum
c                                   Delta N_sc * E_F
              espco2 = espco2 + (sum-sumint)*Dble(ez(ielast))
            End If
            
            If ( Dabs(sumint - sum) .gt. 0.1d0 ) Then
              Write (6,*) '*** Semicore charge is not close to an ',
     $             'integer'
              Write (6,*) '   semicore charge = ',sum
              Stop
            End If
         END IF

          END
