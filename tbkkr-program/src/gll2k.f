
c ************************************************************************
      SUBROUTINE GLL2K(LSTART,NBX,TINVLL,TMATLL,GMATLL,DELTALL,
     +                 EZ,DF,RFCTOR,E2,NSPIN,MAXMESH,NMESH,
     +                 IE,IELAST,IGF,LOFLM,KSCOEF,NSHELL,NBY,NBZ,
     +                 BBYA,CBYA,NAEZ,NATYP,CLS,EQINV,NACLS,RR,
     +                 RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
     +                 BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
     &                 LEFTTINVLL,RIGHTTINVLL,vacflag,NLBASIS,
     &                 NRBASIS,IHANDLE,IHANDLE1,ATOMIMP,
     &                 IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT,
     &                 GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM)
      implicit none
c ************************************************************************
c  input file 25 : determines which green's function elements are
c                  calculated in case of cluster GF determination
c                  (shell structure).
c                  KSCOEF.NE.0                   : input
c                  KSCOEF.EQ.0 .and. IGF.NE.0    : output
c  call KKRMAT1 : k-space integration
c  call TAUTOG1 : calculate structural Greens Function (GMATLL) from 
c                 the scattering path operator (GS).
c-----------------------------------------------------------------------
c     .. parameters ..
      include 'inc.fi'
      include 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER (CONE= (1.D0,0.D0),CZERO= (0.D0,0.D0))
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAX+1)**2)
      DOUBLE PRECISION PI
      PARAMETER (PI= 3.14159265358979312D0)
c     ..
c     .. scalar arguments ..
      DOUBLE COMPLEX DF,EZ          ! E-WEIGHT AND E-POINT
      DOUBLE PRECISION 
     +       BBYA,CBYA,             ! b/a, c/a
     +       E2,                    ! Fermi energy
     +       RFCTOR,                ! ALATA/2/PI
     +       STIME0
      INTEGER ICC,IDGRP,IE,IELAST,IGF,IORIGIN,
     +        IJEND,
     +        KSCOEF,
     +        LATT,LPOT,
     +        MAXMESH,
     +        NAEZ,NATYP,NBY,NBZ,NSPIN,
     +        NMESH,
     +        NZ,                    ! number of atoms at centers of inversion
     +        NSYMAT,NLBASIS,NRBASIS,IHANDLE,IHANDLE1
      INTEGER  IATCONDL(*),IATCONDR(*),NCONDPAIR,IEGFOUT
      LOGICAL LSTART,LIRR,VACFLAG(2)
c     ..
c     .. array arguments ..
      DOUBLE COMPLEX DELTALL(LMAXSQ,LMAXSQ,*),
     +               GMATLL(LMAXSQ,LMAXSQ,*),
     +               GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +               TINVLL(LMAXSQ,LMAXSQ,*),
     +               TMATLL(LMAXSQ,LMAXSQ,*)
      DOUBLE COMPLEX LEFTTINVLL(LMAXSQ,LMAXSQ,*),
     &               RIGHTTINVLL(LMAXSQ,LMAXSQ,*)
      INTEGER ATOM(NACLSD,*),
     +        CLS(*),
     +        EQINV(*),EZOA(NACLSD,*),
     +        KAOEZ(*),
     +        LOFLM(*),
     +        NACLS(*),
     +        NSHELL(0:NSHELD),ATOMIMP(NATOMIMPD)
       DOUBLE PRECISION RBASIS(3,*),RR(3,*),RCLS(3,NACLSD,*),
     +                  BRAVAIS(3,3),RECBV(3,3)
       DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
       DOUBLE COMPLEX GSQDOS(LMAXSQ,LMAXSQ,NAEZD,*)
       DOUBLE PRECISION QDOSKP(3,*)
       INTEGER NQDOSKP,QDOSKDIM
C     ..
c     .. local scalars ..
      DOUBLE COMPLEX CCONST,DFZ,EK,TAUVBZ
      DOUBLE PRECISION ALAT,BSC,CSC,R,STIME,VOLBZ,FACTINV,RMAX
      double precision vol1
      INTEGER I,II,IC,ID,IN,IU,IFILE1,ISHELL,
     +        J,
     +        KS,
     +        L,LM,LM1,LM2,LMAXIN,
     +        M,
     +        N,NBX,NDUM,NHSPIN,NMESHOLD,NOFKS,NS,
     +        POS
      INTEGER NATOMIMP,N1,N2,N3,K
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      CHARACTER*5 STRUCT
      CHARACTER*40 NAME,NEW
      CHARACTER*10 ROTNAME(64)
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX DSYMLL(LMAXSQ,LMAXSQ,NSYMAXD),
     +               GCF(NGFD),
     +               GLL(LMAXSQ,LMAXSQ),
     +               GS(LMAXSQ,LMAXSQ,NSYMAXD,NSHELD),
     +               TLL(LMAXSQ,LMAXSQ)
      DOUBLE PRECISION BZKP(3,KPOIBZ),
     +               RSYMAT(64,3,3),
     +               RATOM(3,NSHELD),
     +               RATOMS(3,NSHELD),
     +               RROT(48,3,NSHELD),
     +               VOLCUB(KPOIBZ)
      INTEGER LF(LMAXSQ),
     +               ISORT(NSHELD),
     +               NSB(2,NSHELD),
     +               NSH1(NSHELD),NSH2(NSHELD),
     +               NXYZ(3),ISYMINDEX(NSYMAXD)
      LOGICAL LSHELL(NSHELD)
C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION DCLOCK
      LOGICAL OPT,TEST
      EXTERNAL DCLOCK,DSORT,SHELLGEN2k,OPT,TEST
      INTRINSIC DATAN,IDINT
C     ..
C     .. EXTERNAL SUBROUTINES ..
C      EXTERNAL BZINTG,BZMESH,KKRMAT,RCSTOP,SYMMGRP,ZAXPY,ZCOPY,ZGEMM
      EXTERNAL BZIRR3D,KKRMAT99,RCSTOP,ROTBRILL,POINTGRP,FINDGROUP,
     +     TAUTOG1,ZAXPY,ZCOPY,ZGEMM
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DSQRT,SQRT
C     ..
C     .. SAVE STATEMENT ..
      SAVE
C     ..
C     ..
      RMAX = 5.d0
      IF(TEST('flow    '))
     + WRITE(6,*) '>>> GLL99: DRIVER FOR BZ INTEGRATION'
C      WRITE(6,*) 'LSTART = ',LSTART
C      WRITE(6,*) 'RFCTOR = ',RFCTOR

      IF (LSTART) THEN
        STIME = DCLOCK()
        WRITE(6,*) 'KSCOEF  = ',KSCOEF
c
c  read the group IDGRP
c
c        ifile1=39
c        IDGRP = 225
c        IORIGIN = 1   
c       CALL READGRP(IDGRP,IORIGIN,RSYMAT,NSYMAT,IFILE1)
c 
c
        CALL POINTGRP(rsymat,rotname)
        CALL FINDGROUP(bravais,recbv,rbasis,rfctor,naez,
     &                 rsymat,rotname,isymindex,nsymat)
c
c
c     find if the inversion is a symmetry op of the real lattice
c     This is going in sub kkrmat99 so that the transpose is
c     added or not !
c
        FACTINV = 1.d0
        do i=1,nsymat
           if (ISYMINDEX(i).eq.25) FACTINV=0.d0
c And in case of 2d              ! Changed on 20.01.2000 
        if (bravais(1,3).eq.0.d0.and.bravais(2,3).eq.0.d0.and.
     &       bravais(3,3).eq.0.d0) THEN
           IF (ISYMINDEX(i).eq.12) FACTINV= 0.d0
        END IF
        end do
c     
c      
c     
        CALL ROTBRILL(DSYMLL,LMAX,NSYMAT,RSYMAT,ISYMINDEX,
     &                LPOT,YR,WTYR,RIJ,IJEND)
c
c Now DSYMLL hold NSYMAT symmetrization matrices
c
        CALL GFSHELLS(ICC,NATOMIMP,NSH1,NSH2,NSHELL,NAEZ,NATYP,
     &                    RBASIS,BRAVAIS,RATOM,RATOMS,RCLSIMP,
     &                    NSYMAT,ISYMINDEX,RSYMAT,RCLS,CLS,
     &                    EQINV,KAOEZ,NACLS,ATOM,ATOMIMP,RFCTOR,
     &                    IATCONDL,IATCONDR,NCONDPAIR)

ccccccc   added on 23.2.2000     
c
        WRITE(6,*) 'ICC     = ',ICC

c
c --->  creates difference vectors RROT for BZ integration in KKRMAT1
c
        DO 107 I=1,NSHELL(0)
          DO 108 IN=1,3
            RATOMS(IN,I) = RATOM(IN,I) 
c     +           - 2.D0*( RBASIS(IN,NSH2(I))-RBASIS(IN,NSH1(I)) )
 108      END DO
 107    END DO
c
        CALL CRTSTAR(RATOMS,NSHELL(0),RSYMAT,NSYMAT,ISYMINDEX,RROT)
c ------------------------------------------------------------------------
        IF (TEST('RROT    ')) THEN
          DO ISHELL=1,NSHELL(0)
            WRITE(6,FMT='(I4)') ISHELL
            WRITE(6,FMT='((I4,3F10.1))') 
     +           (IU,(RROT(IU,I,ISHELL),I=1,3),
     +           IU=1,NSYMAT)   !2*IUMAX*IVMAX)
          END DO
        END IF
c ------------------------------------------------------------------------
        LM = 0
        DO 100 L = 1,LMAX + 1
          DO 90 M = 1,L + L - 1
            LM = LM + 1
            LF(LM) = L
   90     CONTINUE
  100   CONTINUE

        IF (NBY.EQ.0) NBY = NBX/BBYA
        IF (NBZ.EQ.0) NBZ = NBX/CBYA
        IF (OPT('NBZ=0   ') .OR. OPT('SLAB    ')) NBZ = 0
        IF (OPT('WIRE    ')) THEN
          NBX = 0
          NBY = 0
        END IF

        NAME='mesh/kp'

        CALL SNAME(NAME,NEW,LATT)
c	write(6,*) 'name new latt  ',name,new,latt
        CALL SNAME(NEW,NAME,NBX)
c	write(6,*) 'new name nbx   ',new,name,nbx
        CALL SNAME(NAME,NEW,NBY)
c	write(6,*) 'name new nby   ',name,new,nby
        CALL SNAME(NEW,NAME,NBZ)
c	write(6,*) 'new name nbz   ',new,name,nbz

        OPEN(52,FILE=NAME,status='old',ERR=1000)
c  ??????
c        go to 1000
c  ???????
        write(6,*) 'k-mesh file exist : ',name
        DO 242 L=1,MAXMESH
c	write(6,*) 'L MAXMESH',L,MAXMESH
          READ (52,FMT='(I8,3f15.10)') NOFKS,VOLBZ,BSC,CSC
c ------------------------------------------------------------------------
          IF (NOFKS.GT.KPOIBZ) THEN
            write(6,*) 'Error in GLL: (NOFKS.GT.KPOIBZ)'
            write(6,*) 'Please increase the parameter KPOIBZ in ',
     +           'file ''inc.fi''.'
            STOP 'GLL'
          END IF
c ------------------------------------------------------------------------
          READ (52,*) ((BZKP(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
          write(6,FMT=9020) L,NOFKS,VOLBZ
          IF (BSC .NE. 0.D0) THEN
            DO KS = 1,NOFKS
              BZKP(2,KS)=BZKP(2,KS)*BSC/BBYA
              BZKP(3,KS)=BZKP(3,KS)*CSC/CBYA
            END DO
          END IF
          IF (TEST('k-net   ')) THEN
            DO KS = 1,NOFKS
              WRITE(6,9000) (BZKP(I,KS),I=1,3),VOLCUB(KS)
            END DO
          END IF
 242    END DO
        GOTO 1010
       
 1000   write(6,*) 'Create k-mesh, write to file ',name

        DO 230 L=1,MAXMESH
          IF (L.GT.1) THEN
            NBX = (NBX)/1.4
            NBY = (NBY)/1.4
            NBZ = (NBZ)/1.4
          END IF

          LIRR=.TRUE.
          nxyz(1) = nbx
          nxyz(2) = nby
          nxyz(3) = nbz
c	write(6,*) ' L MAXMESH bzirr3d',L,MAXMESH,nbx,nby,nbz
          call BZIRR3D(NOFKS,nxyz,KPOIBZ,BZKP,RECBV,BRAVAIS,RFCTOR,
     &                 VOLCUB,volbz,rsymat,nsymat,ISYMINDEX,LIRR)
c
           write(6,*) 'volbz',volbz 
c
c
          IF (L.EQ.1) OPEN(52,FILE=NAME,status='new')

          WRITE(52,FMT='(I8,3f15.10,/,(3f12.8,d20.10))') 
     +         NOFKS,VOLBZ,BBYA,CBYA,
     +         ((BZKP(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
c ------------------------------------------------------------------------
c
c --->      output of k-mesh created
c
          IF (TEST('k-net   ')) THEN
            DO 200 KS = 1,NOFKS
              WRITE(6,9000) (BZKP(I,KS),I=1,3),VOLCUB(KS) ! /VOLCUB(1)
 200        CONTINUE
          END IF
c ------------------------------------------------------------------------
 230    END DO                      !  L=1,MAXMESH
c ------------------------------------------------------------------------
c
c --->  FOR TEST
c
        IF (TEST('crt mesh')) stop 'in GLL99: ''crt mesh'' '
c ------------------------------------------------------------------------
c
 1010   NMESHOLD = 0
c
C
C --->  ALAT=2*PI*RFCTOR
C
        ALAT = RFCTOR*8.D0*DATAN(1.0D0)       
C
C
        WRITE (6,*) 'IGF    = ',IGF
        IF (IGF.EQ.3) READ (81) NDUM

        WRITE (6,FMT=9070) DCLOCK()-STIME
 9070   format (' TIME IN BZMESH            : ',f9.2)
      END IF                        ! (LSTART)

      LSTART = .FALSE.

      IF (IGF.EQ.3) THEN
        WRITE(6,*) 'READ TINVLL FROM FILE 81'
        DO 150 I=1,NAEZ
          READ (81) EZ,EK,DF, ( (TINVLL(LM1,LM2,I),LM1=1,LMAXSQ),
     +                          LM2=1,LMAXSQ)
 150    CONTINUE
        DF = DF/RFCTOR**2
      END IF

      IF (NMESH.NE.NMESHOLD) THEN
C         WRITE(6,*) '>>> READ K-MESH FROM FILE ''mesh.kp'' :'
        REWIND(52)
        DO 240 L=1,NMESH
          READ (52,FMT='(I8,3f15.10)') NOFKS,VOLBZ,BSC,CSC
          READ (52,*) ((BZKP(ID,I),ID=1,3),VOLCUB(I),I=1,NOFKS)
          IF (BSC.NE.0.D0 .AND. L.EQ.NMESH) THEN
            DO KS = 1,NOFKS
              BZKP(2,KS)=BZKP(2,KS)*BSC/BBYA
              BZKP(3,KS)=BZKP(3,KS)*CSC/CBYA
            END DO
          END IF                    ! BSC.NE.0.D0 .AND. L.EQ.NMESH
 240    END DO
        IF (TEST('k-mesh  ')) THEN
          write(6,*) 'NMESH : ',NMESH
          WRITE(6,*) 'NOFKS : ',NOFKS
        END IF
      END IF                        ! (NMESH.NE.NMESHOLD)

      TAUVBZ = 1.D0/VOLBZ
C ------------------------------------------------------------------------
C
C ---> symmetrize inverted t-matrix with eta-matrices, by (r,r')-symmetry 
c      due to real potential and convert to p.u.
c      changed on 26.4.99 by p.z.
C
      DO 140 I=1,NAEZ
C
        DO 120 LM2 = 1,LMAXSQ
          DO 110 LM1 = 1,LM2
            TINVLL(LM1,LM2,I) = TINVLL(LM1,LM2,I)*RFCTOR
            TINVLL(LM2,LM1,I) = TINVLL(LM1,LM2,I)
 110      CONTINUE
 120    CONTINUE 

C     
C     
C --->  symmetrize inverted t-matrix with dll-matrices
C
        CCONST = 1.d0/DFLOAT(NSYMAT)
        
        DO 130 IU = 1,NSYMAT
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,DSYMLL(1,1,IU),
     +               LMAXSQ,TINVLL(1,1,I),LMAXSQ,CZERO,GLL,LMAXSQ)
          IF (IU.EQ.1) THEN
            CALL ZGEMM('N','T',LMAXSQ,LMAXSQ,LMAXSQ,CCONST,GLL,
     +                 LMAXSQ,DSYMLL(1,1,IU),LMAXSQ,CZERO,TLL,LMAXSQ)
            
          ELSE
            CALL ZGEMM('N','T',LMAXSQ,LMAXSQ,LMAXSQ,CCONST,GLL,
     +                 LMAXSQ,DSYMLL(1,1,IU),LMAXSQ,CONE,TLL,LMAXSQ)
          END IF
          
 130    CONTINUE

        CALL ZCOPY(LMAXSQ**2,TLL,1,TINVLL(1,1,I),1)

 140  CONTINUE                      ! 140 I=1,NAEZ
c       do i=1,naez
c          do lm1=1,lmaxsq
c             write(6,*) 'tmat',i,lm1,tinvll(lm1,lm1,i)
c          end do
c       end do

      EZ = EZ*RFCTOR**2
      DFZ = DF*RFCTOR**2
      STIME0 = DCLOCK()


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       if (TEST('nikos   ')) write(6,*) '>>>before kkrmat99' 
c       write(6,*) ' before kkrmat99 gs',gs(1,1,1,1)
       CALL KKRMAT99(BZKP,NOFKS,GS,VOLCUB,EZ,LF,
     +             TINVLL,TMATLL,DELTALL,RROT,DSYMLL,
     +             NSHELL(0),RFCTOR,IE,IELAST,ALAT,NSYMAT,
     +             NAEZ,CLS,EQINV,NACLS,RR,EZOA,ATOM,KAOEZ,
     +             NSH1,NSH2,GINP,ICC,RBASIS,RCLS,FACTINV,
     &             LEFTTINVLL,RIGHTTINVLL,vacflag,NLBASIS,NRBASIS,
     &             TAUVBZ,RATOM,IHANDLE1,IGF,IATCONDL,IATCONDR,
     &             NCONDPAIR,IEGFOUT,GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM)
c         write(6,*) ' kkrmat99'

       if (TEST('nikos   ')) write(6,*) '>>>after kkrmat99'

      CALL TAUTOG1(GS,TINVLL,DSYMLL,NSHELL(0),RFCTOR,GMATLL,IGF,
     +             TAUVBZ,NSYMAT,NSH1,NSH2,RATOM)
c      write(6,*) ' tautog1 gmatll',gmatll(1,1,1)
      IF (OPT('QDOS    ')) THEN
      CALL TAUTOGQDOS(GSQDOS,TINVLL,DSYMLL,NSHELL(0),RFCTOR,IGF,
     +                TAUVBZ,NSYMAT,NSH1,NSH2,RATOM,NQDOSKP)
      END IF
c
c ------------------------------------------------------------------------
c     it calculates the rest of the G n n' matrix from the 
c     knowledge of the representative pairs (shells) using the
c     real space symmetries (added 23.2.2000)

      IF (ICC.NE.0) THEN

         CALL ROTGLL(GMATLL,RCLSIMP,NATOMIMP,ATOMIMP,RATOM,
     +        NSHELL(0),NSH1,NSH2,DSYMLL,NSYMAT,RSYMAT,ISYMINDEX,
     +        IHANDLE,IGF)

c     changed nshell to nshell(0) in the calling list 4/2/2000


      ENDIF  

      if(test('flow    ')) write(6,*) '<<< GLL99'

      RETURN

 9000 FORMAT(3F12.5,F15.8)
 9010 FORMAT(2F12.6)
 9020 FORMAT(' NMESH : ',I4,'  NOFKS : ',I7,'  VOLBZ :',f14.8)
 9030 FORMAT(I3,I7,I5,F14.3,2F11.3,I8,F8.3)

      END
