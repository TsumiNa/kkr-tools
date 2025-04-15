



c 13.04.00 ***************************************************************
      SUBROUTINE EMESH_SEMI(EZ,DF,NPNT,EBOT,EMU,TK,NPOL,NPNT1,NPNT2,
     +                      NPNT3,NPNTSC,EBOTSC,EUPSC,NPNTIM,EIMSC,
     +                      EBOTGP,EUPGP,EIMGP,NPNTGP,NPNTIMGP)
      Implicit None
C-----------------------------------------------------------------------
C This subroutine provides the energy mesh in array EZ and the
C appropriate integration weights in array DF.
C
C Poles of the Fermi function C (Matsubara frequencies) and
C a contour in the complex energy are used as described in (????).
C
C The contour consists of three straight lines with
C NPNT1, NPNT2, and NPNT3 integration points and is determined by
C the input arguments: EBOT, EMU, TK, and NPOL.
C
C            TK   = temperature in K
C            EMU  = chemical potential in Ry
C            EBOT = bottom of contour in Ry
C            NPOL = number of Matsubara frequencies
C
C The three lines are defined by:
C
C  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK
C              with NPNT1 integration points (Gauss-Legendre rule)
C
C  2. the line from EBOT+2*NPOL*pi*i*k*TK to EMU+(2*NPOL*pi*i-30)*k*TK
C              with NPNT2 integration points (Gauss-Legendre rule)
C
C  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity
C              with NPNT3 integration points (Gauss-Fermi-Dirac rule)
C
C  The total number of integration points is given by:
C              NPNT=NPNT1+NPNT2+NPNT3+NPOL
C
C  The integration points and weights on three lines are chosen
C  according to Gauss integration rules. Only in third interval
C  the Fermi function matters since exp(x) < 10**(-10) for x < -25.
C
C  There are two special cases determined by NPOL = 0 and NPOL < 0.
C
C  a) NPOL = 0 leads to density-of-states calculations
C  with constant integration weights and equally distributed points
C  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.
C
C  The total number of integration points is given by:
C              NPNT=NPNT2
C
C
C  b) NPOL < 0 is meant for calculations where the Fermi-Dirac function
C  is replaced by a step function with step at EMU. When this option is
C  used no poles of the Fermi-Dirac function are used and the contour
C  consists of the three straight lines:
C
C  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK
C              with NPNT1 integration points (Gauss-Legendre rule)
C
C  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK
C              with NPNT2 integration points (Gauss-Legendre rule)
C
C  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU
C              with NPNT3 integration points (Gauss-Legendre rule)
C
C
C  The total number of integration points is given by:
C              NPNT=NPNT1+NPNT2+NPNT3
C
C  Holger's change:
C     If npntsc > 0 then energy points for semicore are added using
C     Simpson's method. Equally spaced npntsc points are added between
C     ebotsc and 'old' ebot. The total number of points is NPNT+npntsc
C
C  Holger's change:
C     If npntgp > 0 then energy points for gap contour are added using
C     Simpson's method. Equally spaced npntgp points are added between
C     ebotgp and 'old' ebot. The total number of points is NPNT+npntgp
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..

      include 'inc.fi'

      DOUBLE PRECISION EBOT,EMU,TK,ebotsc,eimsc,eupsc,
     +                 ebotgp,eimgp,eupgp       
      INTEGER NPNT,NPNT1,NPNT2,NPNT3,NPOL,npntsc,npntgp
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DF(*),EZ(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DE, eimasc, detot,eimagp
      DOUBLE PRECISION ER,ETK,KB,PI, deltae, enow, dfnow
      INTEGER I, npntup, npntdn, npntim,RYD,npntimgp
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WI(IEMXD),XI(IEMXD),XN(IEMXD),WN(IEMXD)
C     ..
C     .. External Subroutines ..
      EXTERNAL GAUFD,RCSTOP,GAULEGNEW
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
C     ..
C     .. Data statements ..
      DATA PI/3.14159265359D0/
      DATA KB/0.6333659D-5/
      DATA RYD/13.6058D0/
      LOGICAL OPT
C     ..

      write(6,*) '>>> EMESH_SEMI: generates a complex E contour'
      ETK = PI*KB*TK
      IF (NPOL.EQ.0) THEN
        DE = (EMU-EBOT)
        IF (NPNT2.GT.1) DE = DE/(NPNT2-1)
        NPNT = 0
        DO 10 I = 1,NPNT2
          NPNT = NPNT + 1
          ER = EBOT + (I-1)*DE
          EZ(NPNT) = DCMPLX(ER,ETK)
          DF(NPNT) = DE
   10   CONTINUE
        WRITE (6,FMT=9000) NPNT,ETK,ETK*RYD
c 9000   format('density-of-states calculation',/,
c     +       'for',I4,' energy points with broadening',E12.4,'Ry'
 9000   format(' density-of-states calculation',/,
     +         ' for',I4,' energy points with broadening',
     +           3p,f10.3,' mRy = ',f10.3,' meV')

      ELSE IF (NPOL.GT.0) THEN

        NPNT = 0

        If ( npntsc .gt. 0 ) Then
c     Use a semicore contour
c     etk = pi*kT
c     eimasc: imaginary part of the semicore energy loop
c
c                      npntsc
c              ----------------------  eimsc
c              |                    |
c       npntup |                    | npntdn
c       ---------------------------------------------> Re-axis
c           ebotsc                ebot
          eimasc = Dcmplx( 0.d0, eimsc )
          npntup = npntim
          npntdn = npntim


          Write (6,*)
          Write (6,*) 'Semicore contour is used, npntsc, ebotsc: ',
     $         npntsc,ebotsc
          Write (6,*) '                          npntup, npntdn: ',
     $         npntup, npntdn
          If ( mod(npntsc,2) .eq. 0 ) Then
            Write (6,*) '*** input npntsc is even ',npntsc
            npntsc = npntsc + 1
            Write (6,*) '*** npntsc used ',npntsc
          End If
          Write (6,*)
          If ( npntsc .lt. 3 ) Stop 'NPNTSC too small'
          If ( (eupsc - ebotsc) .lt. 1.0d-1 ) Then
            Write (6,*) 'EMETIM: EBOTSC is not good'
            Write (6,*) 'eupsc, ebotsc ',ebot,ebotsc
            Stop 'EBOTSC not good'
          End If
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Semicore going up
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ebotsc --> ebotsc + i*eimsc (npntup points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(0.d0,eimsc,xn,wn,npntup,iemxd)
          Call gaulegnew(0.d0,eimsc,xn,wn,npntup)
          Write (6,*)
     $         '*** Using Gauss-Legendre mesh for going up **'
          Do i = 1, npntup
            npnt = npnt + 1
            ez(npnt) = dcmplx(ebotsc, xn(i))
            write(8,*) ez(npnt),npnt
            df(npnt) = dcmplx(0.d0, wn(i))
          End Do
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Semicore horizontal line
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ebotsc + i*eimsc --> ebupsc + i*eimsc (npntsc points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(ebotsc,eupsc,xn,wn,npntsc,iemxd)
          Write (6,*) '*** Using Gauss-Legendre mesh  ***'
          Do i = 1, npntsc
            npnt = npnt + 1
            ez(npnt) = xn(i) + eimasc
            write(8,*) ez(npnt),npnt
            df(npnt) = dcmplx(wn(i),0.d0)
          End Do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Semicore going down
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     eupsc + i*eimsc --> eupsc (npntdn points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(eimsc,0.d0,xn,wn,npntup,iemxd)
          Write (6,*)
     $         '*** Using Gauss-Legendre mesh for going down ***'
          Do i = 1, npntup
            npnt = npnt + 1
            ez(npnt) = dcmplx(eupsc, xn(i))
            write(8,*) ez(npnt),npnt 
            df(npnt) = dcmplx(0.d0, wn(i))
          End Do

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Semicore ends
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

        Else
c     npntsc = 0, so no semicore
          ebotsc = ebot
        End If      
c     npntsc > 0 -if ends
c
c     Valence contour starts
c
c     Gauss-Legendre points from Numerical Recipes routine
        If ( npntgp .gt. 0 ) emu=ebotgp
        CALL GAULEG(XI,WI,NPNT1)
        DE = NPOL*DCMPLX(0.0D0,ETK)
        DO 20 I = 1,NPNT1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT
          write(8,*) EZ(NPNT),NPNT
          DF(NPNT) = WI(I)*DE
   20   CONTINUE
        CALL GAULEG(XI,WI,NPNT2)
        DE = (EMU-30*KB*TK-EBOT)*0.5D0
        DO 30 I = 1,NPNT2
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
          write(8,*) EZ(NPNT),NPNT
          DF(NPNT) = WI(I)*DE
   30   CONTINUE
        CALL GAUFD(XI,WI,NPNT3)
        DE = 30*KB*TK
        DO 40 I = 1,NPNT3
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
          write(8,*) EZ(NPNT),NPNT
          DF(NPNT) = WI(I)*DE
   40   CONTINUE
        DO 50 I = NPOL,1,-1
          NPNT = NPNT + 1
          EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
          write(8,*) EZ(NPNT),NPNT
          DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
   50   CONTINUE

        If ( npntgp .gt. 0 ) Then
c     Use a gap contour
c     etk = pi*kT
c     eimagp: imaginary part of the semicore energy loop
c
c                      npntgp
c              ----------------------  eimgp
c              |                    |
c       npntup |                    | npntdn
c       ---------------------------------------------> Re-axis
c           ebotgp                ebot
          eimagp = Dcmplx( 0.d0, eimgp )
          npntup = npntimgp
          npntdn = npntimgp
          ebotgp = DREAL(EZ(NPNT))


          Write (6,*)
          Write (6,*) 'GAP contour is used, npntgp, ebotgp: ',
     $         npntgp,ebotgp
          Write (6,*) '                          npntup, npntdn: ',
     $         npntup, npntdn
          If ( mod(npntgp,2) .eq. 0 ) Then
            Write (6,*) '*** input npntgp is even ',npntgp
            npntgp = npntgp + 1
            Write (6,*) '*** npntgp used ',npntgp
          End If
          Write (6,*)
          If ( npntgp .lt. 3 ) Stop 'NPNTGP too small'
          If ( (eupgp - ebotgp) .lt. 1.0d-2 ) Then
            Write (6,*) 'EMETIM: EBOTGP is not good'
            Write (6,*) 'eupgp, ebotgp ',eupgp,ebotgp
            Stop 'EBOTGP not good'
          End If
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     GAP going up
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ebotgp --> ebotgp + i*eimgp (npntup points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(0.d0,eimgp,xn,wn,npntup,iemxd)
          Write (6,*)
     $         '*** Using Gauss-Legendre mesh for going up **'
          Do i = 1, npntup
            npnt = npnt + 1
            ez(npnt) = dcmplx(ebotgp, xn(i))
            write(8,*) ez(npnt),npnt
            df(npnt) = dcmplx(0.d0, wn(i))
          End Do
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     GAP horizontal line
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     ebotgp + i*eimgp --> ebupgp + i*eimgp (npntgp points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(ebotgp,eupgp,xn,wn,npntgp,iemxd)
          Write (6,*) '*** Using Gauss-Legendre mesh  ***'
          Do i = 1, npntgp
            npnt = npnt + 1
            ez(npnt) = xn(i) + eimagp
            write(8,*) ez(npnt),npnt
            df(npnt) = dcmplx(wn(i),0.d0)
          End Do
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     GAP going down
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     eupgp + i*eimgp --> eupgp (npntdn points)
c     Gauss-Legendre points from Numerical Recipes routine
c         Call gaulegnew(eimgp,0.d0,xn,wn,npntup,iemxd)
          Write (6,*)
     $         '*** Using Gauss-Legendre mesh for going down ***'
          Do i = 1, npntup
            npnt = npnt + 1
            ez(npnt) = dcmplx(eupgp, xn(i))
            write(8,*) ez(npnt),npnt 
            df(npnt) = dcmplx(0.d0, wn(i))
          End Do

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     GAP ends
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

        Else
c     npntgp = 0, so no gap
          ebotgp = emu
        End If
c     npntgp > 0 -if ends

      ELSE

        IF (NPNT1.GT.0) CALL GAULEG(XI,WI,NPNT1)
        DE = -NPOL*DCMPLX(0.0D0,ETK)
        NPNT = 0
        DO 60 I = 1,NPNT1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT
          DF(NPNT) = WI(I)*DE
   60   CONTINUE
        CALL GAULEG(XI,WI,NPNT2)
        DE = (EMU-EBOT)*0.5D0
        DO 70 I = 1,NPNT2
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
          IF (OPT('GF-EF   ')) EZ(NPNT) = EMU + NPOL*DCMPLX(0.0D0,ETK)
         DF(NPNT) = WI(I)*DE
   70   CONTINUE
        IF (NPNT3.GT.0) CALL GAULEG(XI,WI,NPNT3)
        DE = -NPOL*DCMPLX(0.0D0,ETK)
        DO 80 I = NPNT3,1,-1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EMU
          DF(NPNT) = -WI(I)*DE
   80   CONTINUE



      END IF

      RETURN

      END
