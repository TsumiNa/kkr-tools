c ************************************************************************
      SUBROUTINE CLSGEN99(NATYP,NAEZ,NEMB,RR,NR,RBASIS,
     &                   KAOEZ,Z,CLS,NCLS,NINEQ,EQINV,
     &                   NACLS,ATOM,EZOA, 
     &                   NLBASIS,NRBASIS,NLEFT,NRIGHT,ZPERLEFT,ZPERIGHT, ! new1
     &                   TLEFT,TRIGHT,
     &                   RCLS,RWS,KMT,RMT,MTFAC,RCUT,RCUTXY,L2DIM,
     &                   ALAT)
      implicit none
c ************************************************************************
c This subroutine is used to create the clusters around each atom 
c where repulsive potentials will be positioned.
c
c STRATEGY : 
c Calculate the cluster of each atom by the lattice
c parameters avaliable. Sort the atoms in a unique way :big r, big z, big y
c compare the positions with the previous clusters to see if there is 
c a difference. If not keep only previous clusters and make indexing if
c a new cluster is found then check dimensions and continue for the new
c atom.  
c
c
      include 'inc.fi'
      include 'inc.cls'
c
c
c     .. arguments
c
      DOUBLE PRECISION ALAT         ! lattice constant A
      DOUBLE PRECISION RCUT,RCUTXY
      DOUBLE PRECISION
     +     RBASIS(3,*),             ! pos. of basis atoms in EZ
     +     RCLS(3,NACLSD,*),        ! real space position of atom in cluster
     +     RR(3,0:NRD),             ! set of lattice vectors
     +     RWS(*),
     +     RMT(*),
     +     MTFAC(*),                ! factor to scale RMT (not relly used)
     +     Z(*)                     ! nucleus charge
c
      INTEGER
     +     KMT,                     ! scaling of RMT with MTFAC (not used)
c                                   ! 0: RMT from crystal structure 
c                                   ! 1: RMT - " -  scaled by RMTFAC
c                                   ! 2: RMT = RMTFAC
c                                   ! 3: RMT from ref. pot. card
     +     NATYP,                   ! number of sorts of atoms
     +     NAEZ,                    ! number of atoms in EZ
     +     NEMB,                    ! number of embedding postions
     +     NCLS,                    ! number of diff. clusters
     +     NINEQ,                   ! number of nonequivalent atomic 
                                    ! positions in EZ
     +     NR                       ! number of lattice vectors RR
c
      INTEGER
     +     CLS(*),                  ! sort of cluster around atom
     +     KAOEZ(*),                ! sort of atom at position in EZ
     +     EQINV(*),                ! atom position equivalent by 
                                    ! inversional symmetry
     +     !COCLS(*),                ! center of cluster
     +     NACLS(*),                ! number of atoms in cluster
     +     ATOM (NACLSD,*),         ! index to atom in elem/cell at site in cluster
     +     EZOA (NACLSD,*)          ! index to bravais lattice  at site in cluster
c
c     .. locals
c
      INTEGER 
     +     AJ,C,I,J,N1,INUM,ISUM,
     +     NA,NUMBER,N,NC,NPRIN,ITEST1,ITEST,
     +     POS,IA,IN,IB,II,JATOM,ICU,IC,IAT,I0,I1,ICLUSTER
      INTEGER IATOM(NACLSD),IEZOA(NACLSD),
     +     ISORT(NACLSD),ICOUPLMAT(NAEZD,NAEZD)
c
      DOUBLE PRECISION  
     +     R,R1,R2,RABS,RD,T,EPSSHL,
     +     ASC(3),RCLS1(3,NACLSD),
     +     R0(3,20),RG(3,NACLSD),TMP(3),RSORT(NACLSD)
      INTEGER NLAY,                                
     +        NLBASIS,NRBASIS,                    
     +        NLEFT,NRIGHT           
      DOUBLE PRECISION                             
     +        ZPERLEFT(3),ZPERIGHT(3),            
     +        TLEFT(3,*),TRIGHT(3,*)
      DOUBLE PRECISION RCUT2,RCUTXY2,RXY2 
c
      LOGICAL  L2DIM,CLUSTCOMP
c
c
      LOGICAL TEST,OPT
      EXTERNAL DSORT,CLUSTCOMP
      INTRINSIC MIN,SQRT
c
      DATA     EPSSHL   / 1.0D-4 /
c
c ------------------------------------------------------------------------
      write(6,*) '>>> CLSGEN99: generation of cluster coordinates'
c This is generating the clusters which have a distance smaller
c than RCUT and RCUTXY in plane .
c The cluster atoms are ordered with radious and then z>y>x 
c The ordering allows an easy comparison of clusters
c The principal layer for each layer (atom in unit cell) is
c calculated also for each cluster and the maximum number
c is returned. Some dimension tests are also done      

c      LSPHER=.FALSE.

      WRITE(6,*) 'RCUT = ',rcut,' RCUTXY = ',rcutxy
      IF (ABS(rcutxy - rcut).LT.1.D-4) THEN
          WRITE(6,*) 'Spherical Clusters are created'
c          LSPHER=.TRUE.
      END IF 
      IF (TEST('clusters')) THEN
         OPEN(8,FILE='clusters',status='unknown')
         write(8,9005) NAEZ
         write(8,9030) ALAT
         write(8,9010) (Z(KAOEZ(i)),i=1,NAEZ)
         write(8,9020) (KAOEZ(i),i=1,NAEZ)
      END IF
      
      ICLUSTER = 1
      RCUTXY2 = (RCUTXY+EPSSHL)*(RCUTXY+EPSSHL)
      RCUT2   = (RCUT+EPSSHL)*(RCUT+EPSSHL)
            
      if (TEST('clusters'))  write(21,*) NAEZ*40
      if (TEST('clusters'))  write(21,*) 'Layers '
      DO 2 JATOM = 1,NAEZ       ! loop in all atoms or layers
c Printing the layers
c
      if (TEST('clusters')) THEN
      I0 = Z(JATOM)
      do n=0,39
         write(21,9080) I0, (RR(I,N)+RBASIS(I,JATOM),I=1,3)
      end do
      END IF
c     
c End of printing
         CLS(JATOM) = 0   
c     ----------------------------------------------         
         NUMBER = 0             ! counter for atoms in cluster
         DO NA = 1,NAEZ  ! loop in all atoms
         !write(6,*) 'Calculating atom ',NA!,rbasis(i,na),i=1,3)
            DO N=0,NR    ! loop in all bravais vectors    
               DO I=1,3
                  TMP(I) = RR(I,N)+RBASIS(I,NA)-RBASIS(I,JATOM)
               END DO
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2   =  TMP(3)**2 + TMP(1)**2+TMP(2)**2

               !IF (LSPHER) THEN
               !   RXY2 = RXY2 + RZ2
               !   RZ2 = 0.d0
               !END IF

               IF ( (RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2) )  THEN
                  NUMBER = NUMBER + 1
                  ATOM(NUMBER,JATOM) = NA ! store the atom in elem cell
                  EZOA(NUMBER,JATOM) = N ! store the bravais vector
                  DO I=1,3
                     RCLS1(I,NUMBER) = TMP(I)
                  END DO
c                  !write(6,800) number,na,n,(tmp(i),i=1,3),
c     &      sqrt(tmp(1)**2+tmp(2)**2+tmp(3)**2)
               END IF
 800           format(3I5,4F8.4)
            END DO              ! N loop in bravais
            
         END DO                 ! NA loop in NAEZ
  
c     
c     In the case of 2 dimensional case loop in the atoms
c     outside.
c     
         IF (L2DIM) THEN
c     Somehow meshy
c     ATOM gives the kind of atom 
c   
            DO N=0,NR
               DO I=NLEFT,1,-1  ! loop in some layers on left side
                  DO I1=NLBASIS,1,-1 ! loop in representative atoms on left side
                  DO I0=1,3
                  TMP(I0) = RR(I0,N) + TLEFT(I0,i1) + (I-1)*ZPERLEFT(I0)
     &                               - RBASIS(I0,JATOM)
                  END DO
                  if (TEST('clusters').and.n.lt.39.and.jatom.eq.1) THEN
                     if (I.le.2) then
                     write(21,*) 'Cu ', (tmp(i0),i0=1,3)
                     end if
                  END IF  
                  RXY2 =  TMP(1)**2+TMP(2)**2
                  R2  =   TMP(3)**2+TMP(1)**2+TMP(2)**2
                  !IF (LSPHER) THEN
                  !   RXY2 = RXY2 + RZ2
                  !   RZ2 = 0.d0
                  !END IF
                  
                  IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                     
                     NUMBER = NUMBER + 1
                     ATOM(NUMBER,JATOM) = -NAEZ-I1 ! negative values are used in dlke1.f
                     EZOA(NUMBER,JATOM) = N ! I,I1 are negative
                     DO I0=1,3
                        RCLS1(I0,NUMBER) = TMP(I0)
                     END DO
                  END IF
                  
               END DO 
            END DO
c     
c     
            DO I=1,NRIGHT
               DO I1=1,NRBASIS
               DO I0=1,3
                 TMP(I0) = RR(I0,N)+ TRIGHT(I0,i1) + (I-1)*ZPERIGHT(I0)
     &                             - RBASIS(I0,JATOM)  
               END DO
               if (TEST('clusters').and.n.lt.39.and.jatom.eq.1) THEN
                     if (I.le.2) then
                     write(21,*) 'Ag ', (tmp(i0),i0=1,3)
                     end if
                  END IF
               RXY2 =  TMP(1)**2+TMP(2)**2
               R2  =  TMP(3)**2 + TMP(1)**2+TMP(2)**2
               !IF (LSPHER) THEN
               !   RXY2 = RXY2 + RZ2
               !   RZ2 = 0.d0
               !END IF
               IF ((RXY2.LE.RCUTXY2).AND.(R2.LE.RCUT2)) THEN
                  NUMBER = NUMBER + 1
                  ATOM(NUMBER,JATOM) = -NAEZ-NLBASIS-I1 
                  EZOA(NUMBER,JATOM) = N
                  DO I0=1,3
                     RCLS1(I0,NUMBER) = TMP(I0)
                  END DO
               END IF
            END DO 
         END DO
c     
         
      END DO                    ! loop in all bravais lattices
      END IF                 ! L2DIM Interface calculation
c     
c     Now the atom JATOM Has it's cluster first 
c     sort the atoms of the cluster in increasing order. First by distance
c     Then by z then by y    
c     
         IF (NUMBER.GT.NACLSD) THEN 
            write(6,*) 'Please, increase the parameter NACLSD in',
     &              ' inc.cls to a value greater equal ',number,' .'
               STOP 'Dimension error.'
        
         END IF
         do ia=1,number
            rsort(ia) = SQRT(RCLS1(1,ia)**2+
     &                       RCLS1(2,ia)**2+
     &                       RCLS1(3,ia)**2)
            rsort(ia) = 1000.d0*rsort(ia)-
     &                    10.d0*RCLS1(3,ia)-
     &                   0.1d0*RCLS1(2,ia)-
     &                  0.001d0*RCLS1(1,ia) 
         end do      
c     
         CALL DSORT(RSORT,ISORT,NUMBER,POS)
c     Rearange exchange ia with ib
c MAP temporarily to another array         
         do IA=1,NUMBER       
            do I=1,3
               RG(I,IA)    = RCLS1(I,IA)
            end do
            IATOM(IA) = ATOM(IA,JATOM)
            IEZOA(IA) = EZOA(IA,JATOM)
         end do    
c Now use correct order
         do IA =1,NUMBER
            IB = ISORT(IA)
             do I=1,3
               RCLS1(I,IA) = RG(I,IB)
            end do
            ATOM(IA,JATOM) = IATOM(IB)
            EZOA(IA,JATOM) = IEZOA(IB) 
         END DO
c     
c     Now the clusters have a unique sorting and can be compared with
c     each other Check if ICLUSTER was found previously       
c     
         DO ICU = 1,ICLUSTER-1
            
            N1 = NACLS(ICU)
	    write(6,*) ' clsgen99 icu',icu,n1,number
            IF( CLUSTCOMP(RCLS,ICU,N1,RCLS1,NUMBER) ) THEN ! return true if found before
	     CLS(JATOM) = ICU
            END IF
         END DO
         IF (CLS(JATOM).EQ.0) THEN
            IF (ICLUSTER.GT.NCLSD) THEN
               write(6,*) 'Please, increase the parameter NCLSD in',
     &              ' inc.cls to a value greater equal ',ICLUSTER,' .'
               STOP 'Dimension error.' 
            END IF
            CLS(JATOM) = ICLUSTER
            NACLS(ICLUSTER) = NUMBER
            DO IN = 1,NUMBER
               DO II=1,3
               RCLS(II,IN,ICLUSTER) = RCLS1(II,IN)
               END DO
               write(6,800) jatom,atom(in,jatom),ezoa(in,jatom),
     &                  (rcls1(i,in),i=1,3),
     &      sqrt(rcls1(1,in)**2+rcls1(2,in)**2+rcls1(3,in)**2)
            END DO   
            ICLUSTER = ICLUSTER + 1
         END IF 
c ******************************************************
       write(6,*) 'Atom ',JATOM,' has cluster ', CLS(JATOM),
     &            'with ',NUMBER,' sites'
 2    CONTINUE                      ! JATOM = 1,NAEZ
c
c Now all clusters of all atoms are found print out
c and test the results...
c
c
      DO 22 JATOM = 1,NAEZ
c     
c ------------------------------------------------------------------------
          IC = CLS(JATOM)
          NUMBER = NACLS(IC)
          IF (TEST('clusters')) THEN
            write(8,FMT=1030) NUMBER
            write(8,FMT=1030) JATOM,IC
            DO 105 I=1,NUMBER
              r = sqrt(rcls(1,i,ic)**2+rcls(2,i,ic)**2+rcls(3,i,ic)**2)
c              write(8,1040) (RCLS(II,I,IC),II=1,3),
c     +             KAOEZ(ABS(ATOM(I,JATOM))),ATOM(I,JATOM),r
               write(8,1040) (RCLS(II,I,IC)*ALAT,II=1,3)      
 105        END DO
          END IF
c
c  Now print out the coupling matrix
c
c 
          DO IAT = 1,NAEZ
             ICOUPLMAT(JATOM,IAT) = 0
             DO I=1,NUMBER 
                IF (ATOM(I,JATOM).EQ.IAT) THEN
                   ICOUPLMAT(JATOM,IAT) = 1
                END IF
             END DO
          END DO
          !IF (JATOM.EQ.1) WRITE(6,9050) NAEZ,NAEZ

          WRITE(6,9060) JATOM,(ICOUPLMAT(JATOM,IAT),IAT=1,NAEZ)          

 22       END DO ! Do loop in JATOM (second test loop)
c
c Now testing the shape of the dyson equation
c
          ITEST1 = ICOUPLMAT(NAEZ,1)+ICOUPLMAT(1,NAEZ)
          IF (ITEST1.NE.0) THEN 
          WRITE (6,*) ' This is not a banded matrix ' 
          END IF
          nprin = 0
          do i=1,naez-1
             itest = icouplmat(1,i)
             itest1 = icouplmat(1,i+1)
             if (itest.eq.1.and.itest1.eq.0) then
                nprin = i-1
             end if
          end do
          if (nprin.eq.0) nprin = 1
          WRITE(6,*) '***********************************************'
          WRITE(6,*) '********** TESTING THE COUPLING MATRIX ********'
           WRITE(6,*) '***********************************************'
          WRITE(6,9090) NPRIN
          IF (NPRIN.NE.NPRINCD) THEN 
              WRITE(6,*) 'Please change NPRINCD in your inc.fi file'
              WRITE(6,*) 'from ',NPRINCD, 'to ',NPRIN
              WRITE(6,*) ' ******** RESULTS COULD BE WRONG ******** '
           END IF 
c Now check if you can divide the matrix correctly
          IF (MOD(NAEZ,NPRIN).NE.0) THEN
              WRITE(6,*) ' Your matrix cannot be divided in '
              WRITE(6,*) ' Principal layers. Use a number of layers '
              write(6,*) ' which is multiple of ',NPRIN
          END IF 
c                       NPR 
c NL + 2*NPR*NL - 4* sum  {n}
c                      n=1 
          isum = 0
          do i=1,nprin
             isum = isum + i
          end do
          inum = NAEZ + 2*NPRIN*NAEZ - 2*ISUM
c
c
          isum = 0
          do i=1,naez
             do i0=1,naez
               isum = isum + icouplmat(i,i0)         
             end do
          end do

          IF (ISUM.EQ.INUM) THEN 
              WRITE(6,*) ' Your matrix is BAND DIAGONAL' 
          ELSE 
              WRITE(6,*) ' Your matrix is *NOT* BAND DIAGONAL',ISUM,INUM
          END IF  
          WRITE(6,*) 
          WRITE(6,*) ' Sub clsgen99  exiting <<<<<<<<<<<<<'
c ------------------------------------------------------------------------

 1000   format(' cluster around atom     ',10I4/,
     +         ' number atoms in cluster ',10I4)
 1001   format(I4,2I5,3F8.2,I6,4f7.2)
 1002   format(' cocls : naez =',I3,' (x,y,z)= ',3F10.4)
 1010   format(12x,I6,3F10.4)
 1020   format('  Nr  naez kaoez     x       y       z',
     +         '  ezoa  RR(1)  RR(2)  RR(3)      R')
 1030   FORMAT(3I8)
 1040   FORMAT('Cu  ',3D24.16,2I8,F18.12)
 1050   FORMAT(3F12.7,'  scaling factor')
 1060   FORMAT(I4,3F12.7,'  center',/,(I4,3f12.7))
 1070   FORMAT(I4,3F12.7,'  center of gravity')
 1080   FORMAT('contains ',I4,'  atoms.')
 9005   FORMAT(I4)
 9010   FORMAT(('# Z     ',20F4.0))
 9020   FORMAT(('# KAOEZ ',20I4))
 9030   FORMAT(F12.7,6x,'ALAT')
 9040   FORMAT('> cluster ',I4,' at atom ',I4,
     +       ' of type ',I4,'.')
 9050   FORMAT('      **** COUPLING MATRIX ****'I3,' x ',I3)
 9060   FORMAT(I4,1X,200I1)
 9080   FORMAT(I6,3F15.6)
 9090   FORMAT('The Number of layers is each Principal Layer = ',
     &                I5)
      RETURN
      END
