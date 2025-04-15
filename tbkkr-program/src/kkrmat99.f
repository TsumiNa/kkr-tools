c 28.9.99 ***************************************************************
      SUBROUTINE KKRMAT99(BZKP,NOFKS,GS,VOLCUB,E,LF,
     +                 TINVLL,TMATLL,DELTALL,RROT,DSYMLL,
     +                 NSHELL,RFCTOR,IE,IELAST,ALAT,NSYMAT,
     +                 NAEZ,CLS,EQINV,NACLS,RR,EZOA,ATOM,KAOEZ,
     +                 NSH1,NSH2,GINP,ICC,RBASIS,RCLS,FACTINV,
     &                 TINVBUP,TINVBDOWN,vacflag,NLBASIS,NRBASIS,
     &                 TAUVBZ,RATOM,IHANDLE1,IGF,IATCONDL,IATCONDR,
     &                 NCONDPAIR,IEGFOUT,GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM)
      implicit none
c ************************************************************************
c   performs k-space integration,
c   determines scattering path operator (g(k,e)-t**-1)**-1 and
c   Greens function of the real system -> GS(*,*,*,*),
c
c   NEW VERSION 10.99
c   up -> left , down -> right, for decimation 
c ------------------------------------------------------------------------
c     .. parameters ..
      include 'inc.fi'
      include 'inc.cls'
      INTEGER LMAX,NSYMAXD
C      PARAMETER (LMAX=4)
      PARAMETER (LMAX=LMAXD,NSYMAXD=48)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      INTEGER ALM,NDIM
      PARAMETER (ALM = NAEZD*LMAXSQ,NDIM = NPRINCD*LMAXSQ)
      INTEGER NAUX
      PARAMETER(NAUX = NDIM*NDIM*8)
      DOUBLE COMPLEX CI,CZERO,CONE
      PARAMETER (CI=(0.D0,1.D0),CZERO=(0.D0,0.D0),CONE=(1.D0,0.D0))
c     ..
c     .. scalar arguments ..
      DOUBLE COMPLEX E,TAUVBZ
      DOUBLE PRECISION ALAT,RFCTOR
      INTEGER ICC,IE,IELAST,IUMAX,IVMAX,NAEZ,NOFKS,NSHELL,NZ,NSYMAT
      INTEGER NLBASIS,NRBASIS,IGF,IEGFOUT,NKLINES
c     ..
c     .. array arguments ..
      DOUBLE COMPLEX
     +     DELTALL(LMAXSQ,LMAXSQ,*),
     +     DSYMLL(LMAXSQ,LMAXSQ,*),
     +     GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +     GS(LMAXSQ,LMAXSQ,NSYMAXD,*),
     +     TINVLL(LMAXSQ,LMAXSQ,*),
     +     TMATLL(LMAXSQ,LMAXSQ,*),
     +     TINVBUP(LMAXSQ,LMAXSQ,*),
     +     TINVBDOWN(LMAXSQ,LMAXSQ,*)
      DOUBLE COMPLEX GSK(LMAXSQ,LMAXSQ,NSHELD)
      DOUBLE PRECISION
     +     BZKP(3,*),
     +     RROT(48,3,*),
     +     VOLCUB(*),
     +     RBASIS(3,*),             ! position of atoms in the unit cell
                                    ! in units of bravais vectors
     +     RR(3,0:NRD),
     +     RCLS(3,NACLSD,*),RATOM(3,NSHELD)
      INTEGER
     +     ATOM(NACLSD,*),
     +     CLS(*),
     +     EQINV(*),
     +     EZOA(NACLSD,*),
     +     KAOEZ(*),
     +     NACLS(*),
     +     LF(*),
     +     NSH1(*),NSH2(*),
     &     ICOUPLE(NAEZD,NAEZD),
     &     IATCONDL(*),IATCONDR(*),NCONDPAIR
      LOGICAL VACFLAG(2)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE COMPLEX CARG,GTEMP,g1,g2
      INTEGER DI,DJ,DIM,
     +        I,I1,I2,ID,IETA,II,IICC,
     +        ILM,ILM1,IMAX,INCNP,INFO,IO,IOPT,IROUND,ISYM,IU,
     +        J,J1,JJ,JLM,IL1,IL2,IP1,IP2,
     +        K,II1,II2,LDI1,LDI2,
     +        LL,LM,LM1,LM2,LMJ,
     +        NI,NL,NLAYER,NROUND,NS,NSTAU,NUMGSP,NUMTAU,
     +        IDECI,INVMOD,IER,IL,IHOST,IHANDLE1
      INTEGER L1,l2,ip1t,ip2t,m1,m2,ldi1t,ldi2t,ISTEP1,ISTEP2,ILT1,ILT2
      DOUBLE PRECISION QDECI,TPI,VOL,BOUND,FACTINV,test1
      DOUBLE COMPLEX GSQDOS(LMAXSQ,LMAXSQ,NAEZD,*)
      DOUBLE PRECISION QDOSKP(3,*)
      INTEGER NQDOSKP,QDOSKDIM
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX
     +     AUX(NAUX),
     +     ETAIKR(2*NSYMAXD,NSHELD),   ! ETAIKR(LMAXSQ,LMAXSQ,NSYMAXD,NSHELD),
     +     EXARG(2*NSYMAXD),
     +     G(LMAXSQ,LMAXSQ),GTRANS(LMAXSQ,LMAXSQ),
     +     GK(ALM,ALM),
     +     GSP(NDIM,NDIM,NAUXSPD),
     +     TAU(NDIM,NDIM),
     +     TAUSP(NDIM,NDIM,NSHELD),
     +     TESTLL(LMAXSQ,LMAXSQ)
      DOUBLE PRECISION BZKPK(3),KP(3),WEIGHT
      INTEGER IPVT(ALM),IPIV(NDIM),NLM2,NLM1,N1,N2,
     +     INGSP(NLSPD,NLSPD),
     +     INGSP0(NLSPD,NLSPD),
     +     INTAU(NLSPD,NLSPD),
     +     ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),ITERMAX,ICHCK
      DOUBLE PRECISION ERRMAX,FACTL(0:LMAXD,0:LMAXD)
      LOGICAL LGSP,LIO,TEST,OPT,TLAYPOS,LINTERFACE
      CHARACTER*80 UIO
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CINIT,DLKEB,DLKE0,INVERT,OPT,RCSTOP,
     +         TEST,TLAYPOS,
     +         ZGETRF,ZGETRS,ZXMYPZ
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DATAN,EXP
      DATA BOUND/1.0D-10/
C
c     .. arrays in common
C
      DOUBLE COMPLEX
     +     GLLKE(ALM,ALM),
     +     A1(NDIM,NDIM),B1(NDIM,NDIM),C1(NDIM,NDIM),
     +     AN(NDIM,NDIM),BN(NDIM,NDIM),CN(NDIM,NDIM),
     +     X1(NDIM,NDIM),XN(NDIM,NDIM)
c
C     .. SAVE STATEMENT ..
      SAVE
      DATA LGSP,LIO /2*.TRUE./
c     ..
c ------------------------------------------------------------------------
      if(test('flow     '))
     +     write(6,*) '>>> kkrmat1: loop over k-points'
c
      TPI = 8.D0*DATAN(1.D0)    ! = 2*PI 
c
c Testing the variables
c
      IF ((.NOT.OPT('full inv')).AND.(MOD(NAEZ,NPRINCD).NE.0)) THEN
         WRITE(6,*) 'NAEZ = ',NAEZ,' NPRINCD=',NPRINCD
         WRITE(6,*) 'Use **only** Full inversion in this case '
         STOP
      END IF
      
      IF (OPT('DECIMATE')) THEN
         IF (MOD(NPRINCD,NLBASIS).NE.0) THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NLBASIS=',NLBASIS
            STOP
         END IF
         IF ( MOD(NPRINCD,NRBASIS).NE.0)  THEN
            WRITE(6,*) ' Decimation cannot continue '
            WRITE(6,*) 'NPRINCD=',NPRINCD,' NRBASIS=',NRBASIS
            STOP
         END IF
c
      END IF
c
c Added 1.02.2000
 
      do l1=0,lmax
          do l2=0,lmax
             factl(l1,l2) = (-1)**(l1+l2)
          enddo
       enddo
c ******************************
      NLAYER=NAEZ/NPRINCD  
c     
c     Test ok go on
c     
      
      INCNP = 0                 ! flag to increase NPRINC
      
      DO 10 NS = 1,NSHELL
         DO 20 IU = 1,NSYMAXD
            CALL CINIT(LMAXSQ**2,GS(1,1,IU,NS))
 20      CONTINUE
 10   CONTINUE
c     
c ---> use site symmetry cluster GF
c      G(n,n',L,L')(-k) = G(n',n,L',L)(k)
c
C     Set up the Inversion mode     ----------  29.10.99
c     here it has to construct the icheck matrix, which
c     gives the informations of which elements (in the case of the
c     slab or supercell algorithm) has to be calculated.
        
                                ! changes 29.10.99
       IL = 1   
       CALL IoInput('INTERFACE ',UIO,IL,7,IER)
                    READ (UNIT=UIO,FMT=*) LINTERFACE
        INVMOD = 2                    ! Supercell mode                 

        IF (LINTERFACE) INVMOD = 1    ! band diagonal mode

cccccccccccc &&&&&&&&&&&&&&&&&&&&&&&&&& ccccccccccccccccccccc
c     It prepares the icheck matrix which determines which
c     Green's function elements are neaded.
c    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         CALL GFMASK(ICHECK,ICC,INVMOD,NSH1,NSH2,NLAYER,naez,IE,
     &                  NSHELL,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT)

      IF (OPT('full inv'))  INVMOD = 0     ! PERFORMS FULL INVERSION
         
c                                ----------  29.10.99 
      IF (OPT('QDOS    ')) THEN
         CALL QDOSKNET(NOFKS,NKLINES,QDOSKP,NQDOSKP,IE,QDOSKDIM)
         NOFKS = NQDOSKP
      END IF
      IF (OPT('QDOSEF  ')) THEN
         NQDOSKP = NOFKS
         if (NQDOSKP.gt.QDOSKDIM) then
            write(6,*) 'Increase QDOSKDIM to ....',NQDOSKP
            STOP
         end if
         do i=1,NQDOSKP
            do i1=1,3
               QDOSKP(i1,i) = bzkp(i1,i)
            end do
         end do
      END IF
      IF (OPT('CONDUCT ').AND.(IE.EQ.IEGFOUT)) THEN  
c        CALL FXDRINT(IHANDLE1,NCONDPAIR,1)
c        CALL FXDRINT(IHANDLE1,NOFKS,1)
c        CALL FXDRINT(IHANDLE1,LMAXSQ,1)
         WRITE(6,*) ' I write out:',NCONDPAIR,NOFKS,LMAXSQ
      END IF

      DO 300 K = 1,NOFKS                                 ! K-POINT-LOOP

        IF (OPT('G cont  ') .AND. BZKP(2,K) .gt. .2d0) GOTO 300
        IF (OPT('G_2x2   ') .AND.
     +         ABS(BZKP(1,K))+ABS(BZKP(2,K)).gt. .28d0) GOTO 300
        IF (OPT('GG cont ') .AND. BZKP(2,K) .gt. .05d0) GOTO 300
        IF (OPT('X cont  ') .AND. BZKP(2,K) .le. .2d0) GOTO 300

        IF (OPT('QDOS    ')) THEN
           KP(1) = QDOSKP(1,K)
           KP(2) = QDOSKP(2,K)
           KP(3) = QDOSKP(3,K)
        ELSE
           KP(1) = BZKP(1,K)
           KP(2) = BZKP(2,K)
           KP(3) = BZKP(3,K)
        END IF
        DO 40 NS = 1,NSHELL
          I = NSH1(NS)
          J = NSH2(NS)
          DO 50 ISYM  = 1,NSYMAT

              CARG =  (KP(1)*
     +                  (RROT(ISYM,1,NS)-RBASIS(1,J)+RBASIS(1,I) ) +
     +                  KP(2)*
     +                  (RROT(ISYM,2,NS)-RBASIS(2,J)+RBASIS(2,I) ) +
     +                  KP(3)*
     +                  (RROT(ISYM,3,NS)-RBASIS(3,J)+RBASIS(3,I) )
     +                 )*CI*TPI
              IF (OPT('ONEBULK ')) THEN     ! added on 25.02.2000
c
c For the correct phases look at subroutine dlke1 
c
               CARG =  (KP(1)*(RROT(ISYM,1,NS) ) + 
     &                  KP(2)*(RROT(ISYM,2,NS) ) +
     &                  KP(3)*(RROT(ISYM,3,NS) )
     &                 )*CI*TPI
              END IF
              EXARG(ISYM) = VOLCUB(K)*EXP(-CARG)
              EXARG(NSYMAT+ISYM) = VOLCUB(K)*EXP(+CARG)
c
              ETAIKR(ISYM,NS) = EXARG(ISYM)
              IF (FACTINV.GT.0.5d0) THEN
                 ETAIKR(NSYMAT+ISYM,NS) = EXARG(NSYMAT+ISYM)
              END IF

c ------------------------------------------------------------------------
              IF (TEST('EXARG   ')
c     +             .and.k.eq.1
     +             ) THEN
                write(6,FMT='('' EXARG :'',2f12.4)')
     +               EXARG(ISYM)/VOLCUB(K)      !,EXARG(IV+IVMAX)
                if (k.eq.1)
     +               write(6,*)
     +               (RROT(ISYM,1,NS)+RBASIS(1,J)-RBASIS(1,I)),
     +               (RROT(ISYM,2,NS)+RBASIS(2,J)-RBASIS(2,I)),
     +               (RROT(ISYM,3,NS)+RBASIS(3,J)-RBASIS(3,I))
              END IF
c ----------------------------------------------------------------------


 50       CONTINUE                  ! ISYM = 1,NSYMAT
 40     CONTINUE                    ! NS = 1,NSHELL

        BZKPK(1) = KP(1)
        BZKPK(2) = KP(2)
        BZKPK(3) = KP(3)
c
c ---> fourier transformation
c

        CALL DLKE0(GLLKE,ALAT,NAEZ,
     +       CLS,EQINV,NACLS,RR,EZOA,ATOM,BZKPK,IE,KAOEZ,RCLS,GINP)

      if(k.eq.1) then
      DO N1=1,NAEZ
         DO N2=1,NAEZ
         DO LM1=1,LMAXSQ
             DO LM2=1,LMAXSQ
             NLM1 = (N1-1)*LMAXSQ + LM1
             NLM2 = (N2-1)*LMAXSQ + LM2
c             write(101,8000) N1,N2,LM1,LM2,GLLKE(NLM1,NLM2)/
c     +           (ALAT/2.D0/3.14159265358979312D0)
c,GINP(NLM1,LM1,N2)
             END do
         end do
         end do
      end do
 8000 format(4I5,4F15.8)
      end if
        
        IDECI=0
        IF (OPT('DECIMATE')) IDECI=1

        IF (IDECI.EQ.1) THEN    ! Perform Decimation 3.11.99

c     Parameters for the "decimation" technique.
           itermax = 300
           errmax = 1.0D-10
           ichck = 1

           IF (.NOT.(VACFLAG(1))) THEN  ! if we have no vacuum
c
C Get the matrix B1
c
              CALL BACKCONVERT(1,1,B1,GLLKE)

              if(test('nikos    ').and.(k.lt.5)) then
              CALL BACKCONVERT(NLAYER,NLAYER,BN,GLLKE)
                 write(6,*) 'K-POI=',k,nlayer
                 do lm1=1,ndim
                    do lm2=1,ndim
                       test1 = abs (b1(lm1,lm2)-bn(lm1,lm2))
                       if (test1.gt.1.d-8) then
                          write(6,*) lm1,lm2,b1(lm1,lm2),bn(lm1,lm2)
                       end if
                    end do
                 end do 
                 
              end if            ! testing ok b1=bn
c Now Subtract t-mat of left host
              DO IP1=1,NPRINCD
                 IHOST = NLBASIS - MOD(IP1,NLBASIS) ! get right host atom
                 DO LM1 = 1,LMAXSQ
                 DO LM2 = 1,LMAXSQ
                    IL1 = LMAXSQ*(IP1-1)+LM1
                    IL2 = LMAXSQ*(IP1-1)+LM2   
                    B1(IL1,IL2) =(B1(IL1,IL2) - TINVBUP(LM1,LM2,IHOST)*
     &                                                  RFCTOR) 
                 END DO
                 END DO    
              END DO               

           CALL BACKCONVERT(1,2,C1,GLLKE)
           CALL BACKCONVERT(2,1,A1,GLLKE)

c     it performs the 'space decimation' iterative procedure.        
c             CALL SURFGF(A1,B1,C1,
c    +             X1,NDIM,NDIM,ITERMAX,ERRMAX,ICHCK)
c     adds to the matrix GLLKE the elements that couples the
c     interface to the two half-spaces.
              DO IP1 = 1,NPRINCD
                 DO IP2 = 1,NPRINCD
                    II1=IP1
                    II2=IP2
                    DO LM1 = 1,LMAXSQ
                    DO LM2 = 1,LMAXSQ
                       LDI1 = LMAXSQ*(IP1-1)+LM1
                       IL1  = LMAXSQ*(II1-1)+LM1
                       LDI2 = LMAXSQ*(IP2-1)+LM2
                       IL2  = LMAXSQ*(II2-1)+LM2
c                       GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - X1(LDI1,LDI2)
                    ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              
           END IF  ! if we have no vacuum uper (left) host

           IF (.NOT.(VACFLAG(2))) THEN  ! if we have no vacuum 

c     If 'ONEBULK' is activated then it calculates the xn decimated element
c     from the x1 element: this is just in the case of equal bulks on the 
c     two sides of the interface regions!!!!!!!

              IF (.NOT.OPT('ONEBULK ')) THEN    ! Added 1.02.2000

c     
C     Get the matrix BN
c
             CALL BACKCONVERT(NLAYER,NLAYER,BN,GLLKE)
 
c Now Substract t-mat right host  
! Notes : the indexing is easier like that 
             DO IP1=1,NPRINCD
                 IHOST = NRBASIS - MOD(IP1,NRBASIS)  ! get right host atom
                 DO LM1 = 1,LMAXSQ
                 DO LM2 = 1,LMAXSQ
                    IL1 = LMAXSQ*(IP1-1)+LM1
                    IL2 = LMAXSQ*(IP1-1)+LM2   
                    BN(IL1,IL2) =(BN(IL1,IL2)-TINVBDOWN(LM1,LM2,IHOST)*
     &                                                   RFCTOR)
                    !write(6,*)'right',lm1,lm2,TINVBDOWN(LM1,LM2,IHOST) 
                 END DO
                 END DO    
              END DO

           CALL BACKCONVERT(NLAYER,NLAYER-1,AN,GLLKE)
           CALL BACKCONVERT(NLAYER-1,NLAYER,CN,GLLKE)

           if (test('nikos    ').and.(k.lt.5)) then
              write(6,*) 'K-POI THIRD TIME=',k,nlayer
              do lm1=1,ndim
                 do lm2=1,ndim

                  test1 = abs (a1(lm1,lm2)-an(lm1,lm2))
                  if (test1.gt.1.d-8) then
                   write(6,*) 'A1 CN'
                   write(6,*) lm1,lm2,a1(lm1,lm2),cn(lm1,lm2)
                  end if

                  test1 = abs (c1(lm1,lm2)-cn(lm1,lm2))
                  if (test1.gt.1.d-8) then
                   write(6,*) 'C1 AN'
                   write(6,'(4F10.5)') lm1,lm2,c1(lm1,lm2),an(lm1,lm2)
                  end if
                 end do
              end do 
              end if
           
c     it performs the 'space decimation' iterative procedure.
c             CALL SURFGF(CN,BN,AN,
c    +             XN,NDIM,NDIM,ITERMAX,ERRMAX,ICHCK)
c
c
c
           ELSE                 ! .NOT.OPT('ONEBULK ')  1.02.2000
ccccccccc     
                 do ip1 = 1,nprincd
                 do ip2 = 1,nprincd
                    ip1t = (nprincd+1) - ip2
                    ip2t = (nprincd+1) - ip1
                    do l1 = 0,lmax
                    do m1 = -l1,l1
                       lm1 = l1*(l1+1) + m1 + 1
                       do l2 = 0,lmax
                       do m2 = -l2,l2
                          lm2 = l2*(l2+1) + m2 + 1
                          
                          ldi1 = lmaxsq*(ip1-1) + lm1
                          ldi2 = lmaxsq*(ip2-1) + lm2
                          ldi1t = lmaxsq*(ip1t-1) + lm2
                          ldi2t = lmaxsq*(ip2t-1) + lm1                          
c                  xn(ldi1t,ldi2t) = factl(l1,l2)*x1(ldi1,ldi2)
                       enddo
                       enddo
                    enddo
                    enddo                    
                 enddo
                 enddo
                 
             END IF            ! end of OPT('ONEBULK ')
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             Added on 1.02.2000         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c     adds to the matrix GLLKE the elements that couples the
c     interface to the two half-spaces.
              DO IP1 = 1,NPRINCD
                 DO IP2 = 1,NPRINCD
                    II1=(NLAYER-1)*NPRINCD+IP1
                    II2=(NLAYER-1)*NPRINCD+IP2
                    DO LM1 = 1,LMAXSQ
                    DO LM2 = 1,LMAXSQ
                       LDI1 = LMAXSQ*(IP1-1)+LM1
                       IL1  = LMAXSQ*(II1-1)+LM1
                       LDI2 = LMAXSQ*(IP2-1)+LM2
                       IL2  = LMAXSQ*(II2-1)+LM2
c                       GLLKE(IL1,IL2) = GLLKE(IL1,IL2) - XN(LDI1,LDI2)
                    ENDDO
                    ENDDO
                 ENDDO
              ENDDO

           END IF  ! if we have no vacuum lower (down) host

        ENDIF                   !end of the decimation step        

c     It constructs the matrix M=[-(t)^-1 + G^r] and it is stored
c     in the same matrix GLLKE where G^r was stored.

        DO I1=1,NAEZ
           DO LM1=1,LMAXSQ
              DO LM2=1,LMAXSQ
                 IL1=LMAXSQ*(I1-1)+LM1
                 IL2=LMAXSQ*(I1-1)+LM2
                 GLLKE(IL1,IL2)=GLLKE(IL1,IL2)-TINVLL(LM1,LM2,I1)
c c                write(104,FMT=8001) I1,LM1,LM2,GLLKE(IL1,IL2) 
              ENDDO
           ENDDO
           
        ENDDO
 8001   format(3I5,2F15.8)
 
c     it performs the inversion of the matrix M
c     the output is the scattering path operator TAU!
c       do JLM=1,16
c           do ILM=1,16
c              write(6,*) jlm,ilm,GLLKE(JLM,ILM),tinvll(jlm,ilm,1)
c           end dod
c        end do

        CALL INVERSION(GLLKE,INVMOD,ICHECK)

        DO I1=1,NAEZ
           DO LM1=1,LMAXSQ
              DO LM2=1,LMAXSQ
                 IL1=LMAXSQ*(I1-1)+LM1
                 IL2=LMAXSQ*(I1-1)+LM2
c                 write(105,FMT=8001) I1,LM1,LM2,GLLKE(IL1,IL2) 
              ENDDO
           ENDDO
           
        ENDDO
c 8001   format(3I5,2F15.8)
        
c        do ILM=1,ALM
c           do JLM=1,ALM
c              GKTRANS(ILM,JLM) = GLLKE(JLM,ILM)
c           end do
c        end do
        

        DO 120 NS = 1,NSHELL
           I = NSH1(NS)
           J = NSH2(NS)
           ILM = LMAXSQ*(I-1) + 1

           DO 140 LM = 1,LMAXSQ
              JLM = LMAXSQ*(J-1) + LM
              CALL ZCOPY(LMAXSQ,GLLKE(ILM,JLM),1,G(1,LM),1)
c
c     do the same for the transpose
c
c              CALL ZCOPY(LMAXSQ,GKTRANS(ILM,JLM),1,GTRANS(1,LM),1)
           do lm1=1,LMAXSQ
c               GTRANS(LM1,LM) = GKTRANS(ILM+LM1,JLM) 
                GTRANS(LM1,LM) = GLLKE(JLM,ILM+LM1-1)
           end do       
 140       CONTINUE
c     FACTINV = 0.0 if we have the inversion
c     Else it is 1.d0
           DO 110 ISYM = 1,NSYMAT
              do lm1=1,lmaxsq
                 do lm2=1,lmaxsq
                    GS(lm1,lm2,ISYM,NS) =
     &                   ETAIKR(ISYM,NS)*G(lm1,lm2)+
     &           FACTINV*ETAIKR(ISYM+NSYMAT,NS)*GTRANS(lm1,lm2)+ ! the transpose is added
     &                   GS(lm1,lm2,ISYM,NS)
                 end do
              end do
 110       CONTINUE             ! ISYM = 1,NSYMAX
c c        write(6,*) ' after 110 kkrmat99 gs',gs(1,1,1,1)
c    Added for the q-dos  19.4.2000
              IF (OPT('QDOS    ').or.OPT('QDOSEF  ')) THEN
                 IF (NS.GT.NAEZD) THEN
                 WRITE(6,*) ' QDOS Dimension problem !!'
                 STOP
                 END IF
                 do lm1=1,lmaxsq 
                    do lm2=1,lmaxsq
                       GSQDOS(lm1,lm2,NS,K) = G(lm1,lm2)
                    end do
                 end do
c              if (ns.eq.4) then
c              write(6,1001) K, bzkpk(1),GSQDOS(1,1,ns,K),etaikr(1,ns)
c 1001         format(i5,6D16.6)
c              end if
              END IF
ccccc    q-dos
c
c  This is for the case of k-resolved Green's function
c  new added 24.2.2000 
           IF (OPT('ONEBULK ')) THEN ! added on 25.02.2000
c     
c     Just multiply the phase out 
c     
              CARG =  (KP(1)*(RBASIS(1,J)-RBASIS(1,I) ) +
     +                 KP(2)*(RBASIS(2,J)-RBASIS(2,I) ) +
     +                 KP(3)*(RBASIS(3,J)-RBASIS(3,I) )
     +                )*CI*TPI
              CARG = 0.d0  
           ELSE
              CARG = 0.d0
           END IF
           do lm1=1,lmaxsq
             do lm2=1,lmaxsq
           GSK(lm1,lm2,NS) = G(lm1,lm2)*EXP(-CARG)
c           if ((k.lt.2).and.(IE.EQ.IEGFOUT)) then
c           write(6,9980) ns,lm1,lm2,gsk(lm1,lm2,ns)
c           end if
             end do
             end do
 9980        format('LA ',3I5,4D18.6)
 120    CONTINUE                ! NS = 1,NSHELL
        IF (OPT('CONDUCT ').AND.(IE.EQ.IEGFOUT)) THEN      !
           WEIGHT = VOLCUB(K)
           write(6,*) ' K-Point ....',K
           CALL KRESTAU2G(KP,WEIGHT,GSK,TINVLL,NSHELL,NSH1,NSH2,
     &          RFCTOR,RBASIS,IHANDLE1,IATCONDL,IATCONDR,
     &          NCONDPAIR,IEGFOUT,IGF)
        END IF
 300  END DO                    ! K = 1,NOFKS


      IF(TEST('flow    ')) write(6,*) '<<< KKRMAT1'
      IF (INCNP.NE.0 .AND. IE.EQ.IELAST)
     +     WRITE(6,*) 'KKRMAT1: Please increase the parameter ',
     +     'NPRINCD in file inc.fi.'

      RETURN

 9000 FORMAT(20I4)
 9010 FORMAT(3f12.5)
 9011 FORMAT('K = ',i4,3f12.5)
 9030 FORMAT(A1,4x,':',I6,2d12.3,2f10.4)
 9040 FORMAT(' NUMTAU =',I6)
 9050 FORMAT(' NUM_LR =',I6)

      END
