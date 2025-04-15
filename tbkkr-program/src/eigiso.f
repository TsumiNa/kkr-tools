c ************************************************************************
      SUBROUTINE EIGISO(HEAD,NSPIN,NAEZ,NINEQ,LMAX,LSYM,E2,
     +     VCONST,QKP,NPOINTS,QBOUND,AKSCALE)
      implicit none
c ************************************************************************
c
c ---> calculate section points of an isoenergetic surface (e=EFERMI)
c      with given lines in the reciprocal space
c
c ------------------------------------------------------------------------
      include 'inc.fi'
C
c     .. arguments
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
C
      INTEGER LMAX,NPOINTS,NSPIN,NINEQ,NAEZ
      LOGICAL LSYM
      DOUBLE PRECISION
     +     AKSCALE,
     +     E2,
     +     VCONST,
     +     QBOUND,QKP(6,*)
      CHARACTER*1          HEAD(80)
c ------------------------------------------------------------------------
c     .. locals
      DOUBLE PRECISION  
     +     CSQR,
     +     DE,DOS,
     +     E1,EA,EE,EM,EFERMI,EDIFF,
     +     FAC,FACI,
     +     TIMEI,TIMES,TIMET,
     +     W1,W2
      DOUBLE PRECISION 
     +     EQ(0:NAEZD*LMMAXD),
     +     Q(6)
      DOUBLE COMPLEX E
      INTEGER
     +     ILAY(NAEZD),             ! layer of atom
     +     ISTEP(NAEZD*LMMAXD),
     +     NLAY(NAEZD)               ! number of atoms in layer
      INTEGER
     +     I,I1,I2,IA,IBAND,IE,IL,IP,ILOOP,ISPIN,IVAL1,IVAL2,
     +     J,JEND,
     +     LAYMAX,L,LM,LOOP,
     +     N,NA,NE,NL,NLINES,NUM1,NUM2
C
      DOUBLE PRECISION PI,TPI,ZERO,ONE
      PARAMETER (
     +     PI    = 3.14159265358979312D0,
     +     TPI   = 2.0d0 * PI,
     +     ZERO  = 0.0D0,
     +     ONE   = 1.0D0)
C
      DOUBLE PRECISION DCLOCK
      LOGICAL OPT,TEST
      EXTERNAL DCLOCK, EIGEN0,EINIT,EREST,MINEIG,OPT,RINIT,TEST
      INTRINSIC DMIN1,DSIGN,DABS

      CHARACTER*9 FILEN(3,2)
      DATA FILEN / 'isodown  ','isodown.+','isodown.-',
     +             'isoup    ','isoup.+  ','isoup.-  ' /
C
C     .. COMMONS
C
      INTEGER     IUP,IDO,NUM,ITOT,N2
      DOUBLE PRECISION WUP,WDO
      COMMON /EIGENV/ WUP,WDO,NUM,IUP,IDO

      INTEGER NBASIS(KPOIBZ),NMIN(KPOIBZ),NMAX(KPOIBZ)
      COMMON /NBASIS / NBASIS,NMIN,NMAX
      
c ------------------------------------------------------------------------
      OPEN(29, FILE='lines',STATUS='OLD')

      TIMET = DCLOCK()
      ITOT = 0

      IP = 1
      JEND = 3
      IF (OPT('COMPLEX ')) JEND = 6

      LOOP = 3
      IF (TEST('at ef   ')) LOOP = 1
      IF (LOOP.GT.1 .AND. ABS(VCONST).LT.1.D-10) VCONST=1.D-4
      write(6,FMT=9030) LOOP,VCONST

      DO 10 ILOOP=1,LOOP
        IF (ILOOP.EQ.2) E2 = E2 + VCONST
        IF (ILOOP.EQ.3) E2 = E2 - 2.d0*VCONST
        IF (TEST('no ef-de') .AND. ILOOP.EQ.3) GOTO 10
        IF (TEST('no ef+de') .AND. ILOOP.EQ.2) GOTO 10
        IF (TEST('no ef   ') .AND. ILOOP.EQ.1) GOTO 10

        DO 310 ISPIN = 1,NSPIN

          IF (TEST('UP      ') .AND. ISPIN.EQ.1) GOTO 310
          IF (TEST('DOWN    ') .AND. ISPIN.EQ.2) GOTO 310

          OPEN(31,FILE=FILEN(ILOOP,ISPIN),STATUS='UNKNOWN')
          WRITE(31,103) (HEAD(I),I=1,76)

 103      FORMAT('!! ',76A1)

          TIMES = DCLOCK()

          WRITE(6,*) 'ISPIN :',ISPIN
          REWIND(29)
          READ (29,*) NLINES,N
          IF (N.NE.NPOINTS) THEN
            write(6,*) 'N(',N,').NE.NPOINTS(',NPOINTS,')'
          END IF
          WRITE(31,104) NLINES,NPOINTS,E2
          WRITE(6,104) NLINES,NPOINTS,E2
 104      FORMAT(2I8,1p,d20.10,'   NLINES NPOINTS EFERMI')
C     
          E = CMPLX(E2,ZERO)
c
c --->    loop over lines
c
          DO 360 IL = 1,NLINES
            READ (29,*) NL,IA,IE

c ------------------------------------------------------------------------
            IF (TEST('NFE     ')) THEN
              EA=sqrt(QKP(1,IA)**2+QKP(2,IA)**2+QKP(3,IA)**2)*AKSCALE
              EE=sqrt(QKP(1,IE)**2+QKP(2,IE)**2+QKP(3,IE)**2)*AKSCALE
              DE=sqrt((QKP(1,IA)-QKP(1,IE))**2
     +             +(QKP(2,IA)-QKP(2,IE))**2
     +             +(QKP(3,IA)-QKP(3,IE))**2)*AKSCALE
              if (abs(ea-sqrt(e2)).gt.1.3d0*de .or.
     +            abs(ee-sqrt(e2)).gt.1.3d0*de ) goto 360
            END IF
c ------------------------------------------------------------------------

            IF (TEST('eigiso  ')) write(6,9020) NL
            IF (TEST('eigiso  ')) write(6,*) ' KPA,KPE :',IA,IE
c
c --->      number of eigenvalues at begin and end of the line
c
            CALL EIGEN0(ISPIN,QKP(1,IA),E,LSYM,0)
            NA = NUM                ! + NPOLES(E2,ISPIN)

            CALL EIGEN0(ISPIN,QKP(1,IE),E,LSYM,0)
            NE = NUM                ! + NPOLES(E2,ISPIN)

            IF (TEST('eigiso  ')) write(6,*) '   NA,NE :',NA,NE

            IF (NA.NE.NE) THEN
c
c --->        line is crossed by the isoenergetic plane
c
c     
              DO 350 I=1,NAEZD*LMMAXD
                ISTEP(I) = 0
 350          END DO
              CALL EINIT(IP)
              FACI = ZERO
c
              IF (NA.GT.NE) THEN
c
c --->          interchange begin and end of line
c
                FACI = ONE
                I = IA
                IA = IE
                IE = I
                I = NA
                NA = NE
                NE = I
              END IF

              NMIN(IP) = NA
              CALL ESAVE(ZERO,NA,IP)

              NMAX(IP) = NE
              CALL ESAVE(ONE,NE,IP)

              I = NA+1
              FAC = ONE
              IF (ABS(E2) .LT. 0.1D0) FAC = 5.D-2
c ------------------------------------------------------------------------
c
c --->      loop for section point search
c
              DO 330 WHILE (I.LE.NE)
c                WRITE(6,*) 'I:',I
c
c --->          determine intervall for incremental searching
c
                CALL EREST(IP,EA,EE,I,NUM1,NUM2,W1,W2,IVAL1,IVAL2)
                N2 = NUM2
                
                DO 320 WHILE (EE-EA .gt. FAC*QBOUND)
                
                  EM = (EA + EE)/2.0D0

                  DO 370 J = 1,JEND
                    Q(J) = QKP(J,IA) + EM*(QKP(J,IE)-QKP(J,IA))
 370              END DO
                  
                  CALL EIGEN0(ISPIN,Q,E,LSYM,0)
                  
                  N = NUM           ! + NPOLES(E2,ISPIN)
                  
                  IF (N.GE.NMIN(IP)) CALL ESAVE(EM,N,IP)
                  
                  IF (TEST('eigiso  ')) 
     +                 write(6,fmt='(''em,n :'',f14.9,i6)') em,n
                  
                  IF (N.GE.I) THEN
                    EE    = EM
                    N2    = N
                    W2    = WDO
                    IVAL2 = IDO
                  ELSE
                    EA    = EM
                    W1    = WUP
                    IVAL1 = IUP
                  END IF
                  
                  ISTEP(I) = ISTEP(I) +1 ! step counter
                  ITOT  = ITOT + 1
 320            END DO              ! WHILE (EE-EA .gt. FAC*QBOUND)
              
                IF (TEST('eigiso  ')) THEN
                  write(6,FMT=9008) ea,ee
                  write(6,FMT=9009) ival1,ival2,w1,w2
                END IF
  
                IF (IVAL1*IVAL2.GT.0) THEN
                  EM = (W1*EE-W2*EA)/(W1-W2)
                ELSE
                  EM = (EA+EE)/2.0D0
c                  IF (QBOUND.GT. 1.0D-6) 
c     +                 STOP 'Decrease QBOUND in in5 !!'
                END IF
                
                IF (N2.GT.NE) N2 = NE
                
                DO 380 J = I,N2
                  EQ(J) = EM*(1.0D0-2.0D0*FACI) + FACI
                  write(6,FMT=9007) NL,J,EQ(J)
 380            END DO
                
                I = N2 + 1
                
 330          END DO                ! WHILE (I.LE.NE)
              
              IF (TEST('eigiso  ')) 
     +             write(6,FMT=9006)
     +             (I,EQ(I),ISTEP(I),I=NA+1,NE)
              
              write(31,9005)  NL,NE-NA,(I,EQ(I),I=NA+1,NE)
              
            END IF                  ! (NA.NE.NE)
            
 360      END DO                    ! IL=1,NLINES
          
          write(6,FMT=9011) DCLOCK()-TIMES
          
          CLOSE(31)
          
 310    END DO                      ! ISPIN = 1,NSPIN

 10   END DO                        ! ILOOP=1,LOOP
c     
      write(6,*) 
     +     'Total number of steps ',ITOT,' for ',NLINES,' lines.'
c     
      close(29)
      close(30)
      
      RETURN

 9005 format(I8,I4,/,(I4,f16.12))
 9006 format(i5,F18.6,i8)
 9007 format(' IL,IB :',2I6,F18.10)
 9008 format(' EA,EE :',2f14.9)
 9009 format(' VAL,W :',2I6,2d14.3)
 9010 format(' time in IL    loop :',f12.2,' sec.')
 9011 format(' time in ISPIN loop :',f12.2,' sec.')
 9020 format(' IL    :',I6)
 9030 format(' LOOP,VCONST :',I6,1p,d10.2)
        
      END
