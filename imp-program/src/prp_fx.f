c ======================================================================
*DECK prp_fx
      SUBROUTINE PRP_FX(IOFILE,KPRI,K_PRPFX,K_SET,MBRYD,NATOM,NFATOM,
     +                  IFATOM,TAUPRO,IATDYN,NFATOM_DYN,NBRYD,F1XYZ,
     +                  TAUXYZ,F_M,XKS_M,F_OLD,XKS_MP1,TAU,XKS_I)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c $Id: prp_fx.f,v 1.2 1998/12/09 08:38:50 wikromen Exp $
C
C    This subroutine is used for
c    assignment and reassignment of forces and coordinates as dyn variab
c    It replaces 'set position m'
c            and 'RESET COORDINATES FOR NEXT ITERATION' in mdyncs.f
c
c    assignment and reassignment is governed by two keys:
c
c    k_prpfx : reducing number of dyn atoms;
c              using symmetry or not;
c              projecting to hyperplane or not
c
c    k_set : 1 assign forc and coord to dynamical variables;
c            2 reassign dyn variables to forces and coordinates
c                           Kurt Schroeder, IFF, KFA Juelich, Oct. 1996

C-----------------------------------------------------------------------
c $Log: prp_fx.f,v $
c Revision 1.2  1998/12/09  08:38:50  wikromen
c ID and LOG added; print ID included;
c
c-----------------------------------------------------------------------
C
C----> USE STATMENTS
C

C     IMPLICIT NONE
C
      INTEGER NATOMD
      PARAMETER (NATOMD=102)
      INTEGER MBRYD_MD
      PARAMETER (MBRYD_MD=3*NATOMD)

C---> FILE NUMBER FOR READ AND WRITE
C
      INTEGER IOFILE
C
C---> RUNNING MODE PARAMETERS
C
      INTEGER KPRI
csym  integer ksym
      INTEGER K_PRPFX,K_SET
C
C---> STRUCTURAL INFORMATION
C
cccc  integer natypd
      INTEGER NATOM,NKATOM
      REAL*8 TAUXYZ(3,NATOMD),TAU(3,NATOMD)
      REAL*8 DEL_TAUXYZ(3,NATOMD)
C
C---> FORCE INFORMATION
C
cccc  integer natypd
      INTEGER NFATOM
      INTEGER IFATOM(NATOMD)
      REAL*8 F1XYZ(3,NATOMD)
C
C---> SYMMETRY  INFORMATION
C
cccc  integer msym
csym  integer ksym,nsym
csym  integer ifshell(2,msym,natypd)
csym  real    crot(3,3,msym)
c
c---> quasi newton step variables
c
cccc  integer mbryd_md
      INTEGER MBRYD
      INTEGER NBRYD
      REAL*8 XKS_M(MBRYD),XKS_MP1(MBRYD),XKS_I(MBRYD)
      REAL*8 F_M(MBRYD)
      REAL*8 F_OLD(MBRYD_MD)
      INTEGER NFATOM_DYN
      REAL*8 TAUPRO(3,NATOMD)
      INTEGER IATDYN(NATOMD)

C
C---> LOCAL VARIABLES
C
      INTEGER ICNT,JCNT,IJCNT,IAT,JAT,IDYN,JDYN
      REAL*8 FREST(3,NATOMD),FPAR(NATOMD)
      REAL*8 TREST(3,NATOMD),TPAR(NATOMD)
      REAL*8 TAUHLP(3,NATOMD),TAUHLPN(NATOMD)
      REAL*8 ZERO,MRYD
      PARAMETER (ZERO=0.0d0,MRYD=1.0d3)
C
C---> ABBREVIATIONS
C
c      k_prpfx: key for assigning dynamical variables
c               1  no special condition
c              11  atoms in upper half of slab mobile plus inv images
c              21  symmetry acc to symmetry file (atom # in ifshell)
c               2  atoms out of atom-list iatdyn mobile in z-dir and
c                        all other atoms mobile without restriction
c              12  atoms out of atom-list iatdyn mobile in z-dir and
c                  atoms in upper half of slab mobile plus inv images
c               3  only atoms out of atom-list iatdyn mobile
c              13  only atoms out of atom-list iatdyn mobile plus inv
c                  images
c               4  only first atom mobile in z-direct
c              14  only first atom mobile in z-direct plus inv image
c               6  only atoms out of atom-list iatdyn mobile in restr
c                  geom
c              16  only atoms out of atom-list iatdyn mobile in restr
c                  geom
c                                                     plus inv images
c               7  atoms out of atom-list iatdyn mobile in restr geom
c                  and all other atoms mobile without restriction
c              17  atoms out of atom-list iatdyn mobile in restr geom
c                  and atoms in upper half of slab mobile plus inv
c                  images
c               8  all atoms mobile in restr geometry
c              18  all atoms mobile in restr geometry plus inv images
c             135  special for 110-surf of 3-5-comp semiconductors:
c                        atoms in upper half of slab mobile
c                                              plus yz-plane-mirr images
c
c      k_set  : key for call of prp_fx
c               1  assign coordinates and forces as dynamical variables
c               2  reassign dynamical variables to coordinates and
c                  forces
c      ksym   : key for symmetry; 0: no symmetry; 1: symmetry in
c               sym_file
c      crot   : rotation matrices for symmetry
c      ifshell: atomic # of atom images created by application of crot
c
C      NFATOM : NUMBER OF MOBILE ATOMS
C      TAUXYZ : CARTESIAN COORDINATES OF ALL ATOMS IN A.U.(1:3,1:NATOMD)
C               ARRAY IS CHANGED IN THE SUBROUTINE
C      F1XYZ  : CONTAINS THE X,Y,Z--COMPONENTS OF THE
C               TOTAL HELLMANN--FEYNMAN FORCE (EWALD-FORCE DUE TO ALL
C               BARE IONS AND FORCE DUE TO ELECTRON DISTRIBUTION)
C               OF ALL ATOMS FOR WHICH THE FORCE SHOULD BE CALCULATED.
C               F1XYZ IS CALCULATED IN THE SUBROUTINE FORCE
C               WITHOUT UPDATING ELECTRONIC PARAMETERS
C      taupro : CARTESIAN COORDINATES OF projection direction;
C               hypervector: (x1,y1,z1;....;xN,yN,zN)
C               dynamics considered in plane perpendicular to taupro
C      iatdyn : list of atom-Nrs (corresponding to ifatom)
c               which are considered in the dynamics
C
C----------------------------------------------------------------------
c
c
c ---> PRINT REVISION-NUMBER
c
      LOGICAL FIRST_CALL

      FIRST_CALL = .true.
      IF (FIRST_CALL) THEN
          FIRST_CALL = .false.
          WRITE (IOFILE,FMT=*) '$RCSfile: prp_fx.f,v $ $Revision: 1.2 $'
      END IF
c
cholger
      IJCNT = 0
      DO ICNT = 1,3
         DO JCNT = 1,NFATOM
            IJCNT = IJCNT + 1
            IAT = IFATOM(JCNT)
            XKS_I(IJCNT) = TAU(ICNT,IAT)
         END DO
      END DO
c end
      IF (KPRI.GE.1) THEN
          WRITE (IOFILE,FMT='(/)')
          WRITE (IOFILE,FMT='(''       *<* prp_fx *>* '')')
          WRITE (IOFILE,FMT='(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~ '')')
          WRITE (IOFILE,FMT=2900)
          WRITE (IOFILE,FMT='(/)')
      END IF
c
      IF (K_SET.NE.1 .AND. K_SET.NE.2) THEN
          WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
          STOP 'k_set in prp_fx'
      END IF
c
c
      IF (K_PRPFX.EQ.1) THEN
c     all atoms mobile; no special condition
c
          IF (NFATOM_DYN.NE.NFATOM) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn .ne. nfatom'
              STOP 'nfatom_dyn in prp_fx 1'
          END IF
c
          NBRYD = 3*NFATOM
c
c---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      IAT = IFATOM(JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                  END DO
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END DO
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 1'
          END IF
c
      ELSE IF (K_PRPFX.EQ.2) THEN
          STOP 'prp_fx : NOT TESTET !'
c     atoms specified in atom-list iatdyn mobile in z-dir;
c     all other atom without restriction
c                              no special symmetry
c     nfatom_dyn here has the meaning 'number of atoms mobile in z dir'
c
          IF (NFATOM_DYN.GT.NFATOM) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn .gt. nfatom'
              STOP 'nfatom_dyn in prp_fx 2'
          END IF
c
          NBRYD = 3*NFATOM - 2*NFATOM_DYN
c
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      ICNT = 3
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                          XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      ICNT = 3
                      IJCNT = IJCNT + 1
cccc                  del_tauxyz(icnt,iat) = xks_mp1(ijcnt)
cccc >                                     - tauxyz(icnt,iat)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
cccc                  del_tauxyz(icnt,iat) = xks_mp1(ijcnt)
cccc >                                     - tauxyz(icnt,iat)
                          TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                          F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
                      END DO
                  END IF
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 2'
          END IF
c
      ELSE IF (K_PRPFX.EQ.3) THEN
c     only atoms out of atom-list iatdyn mobile ;
c                              no special symmetry
c
          NBRYD = 3*NFATOM_DYN
          IJCNT = 0
          DO JDYN = 1,NFATOM_DYN
              JCNT = 0
              DO IDYN = 1,NFATOM
                  IF (IFATOM(IDYN).EQ.IATDYN(JDYN)) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no corresponding atom found in ifatom'
cholger                    stop 'iatdyn in prp_fx 1'
              END IF
              IAT = IFATOM(JCNT)
              DO ICNT = 1,3
                  IJCNT = IJCNT + 1
c
c---> set position m
c
                  IF (K_SET.EQ.1) THEN
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
                  ELSE IF (K_SET.EQ.2) THEN
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END IF
              END DO
          END DO
c
      ELSE IF (K_PRPFX.EQ.4) THEN
          STOP 'PRP_fx : NOT TESTET'
c     only first atom (for which force is calc) mobile in z-direct;
c                                         no special symmetry
c
          IF (NFATOM_DYN.NE.1) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,'  nfatom_dyn .ne. 1'
              STOP 'nfatom_dyn in prp_fx 4'
          END IF
c
          NBRYD = NFATOM_DYN
          JCNT = 1
          IAT = IFATOM(JCNT)
          ICNT = 3
c
c---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 1
              F_M(IJCNT) = F1XYZ(ICNT,JCNT)
              XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 1
              TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
              F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
              F1XYZ(ICNT,JCNT) = ZERO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 5'
          END IF
c
      ELSE IF (K_PRPFX.EQ.6) THEN
          STOP 'prp_fx : NOT TESTET'
c     only atoms out of atom-list iatdyn mobile in restricted geometry;
c                                       no special symmetry
c
c
c ---> print out incoming fields (for testing)
c
          WRITE (IOFILE,FMT=*) K_PRPFX,K_SET,' k_prpfx,k_set'
          WRITE (IOFILE,FMT=*) NFATOM_DYN,' nfatom_dyn'
          DO JDYN = 1,NFATOM_DYN
              WRITE (IOFILE,FMT=
     +'(1x,i4,3(d17.10,1x),1(1x,i4),                              ''  ia
     +tdyn,taupro,jdyn'')') IATDYN(JDYN), (TAUPRO(ICNT,JDYN),ICNT=1,3),
     +          JDYN
          END DO
c      end testing
c
          NBRYD = 3*NFATOM_DYN
c
c
c---> print (and sort) projection vector
c
          DO JDYN = 1,NFATOM_DYN
              IAT = IATDYN(JDYN)
              WRITE (IOFILE,FMT=*) IAT,JDYN,' iat=iatdyn(jdyn),jdyn'
              WRITE (IOFILE,FMT='(a35)')
     +          ' hypervector perpend to hyperplane '
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JDYN),ICNT=1,3),JDYN
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JDYN)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JDYN,IAT
c
c---> (1) calculate square norm of projection vector
c
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (2) calculate parallel components
c
          IJCNT = 0
          DO JDYN = 1,NFATOM_DYN
              JCNT = 0
              DO IDYN = 1,NFATOM
                  IF (IFATOM(IDYN).EQ.IATDYN(JDYN)) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no corresponding atom found in ifatom'
                  STOP 'iatdyn in prp_fx 6'
              END IF
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
              DO ICNT = 1,3
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
c
c---> (3) calculate force and coord restricted to hyperplane
c
              DO ICNT = 1,3
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
              DO ICNT = 1,3
                  IJCNT = IJCNT + 1
c
c---> set position m
c
                  IF (K_SET.EQ.1) THEN
                      F_M(IJCNT) = FREST(ICNT,JCNT)
                      XKS_M(IJCNT) = TREST(ICNT,IAT)
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
                  ELSE IF (K_SET.EQ.2) THEN
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TREST(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                   DEL_TAUXYZ(ICNT,IAT)
                      F_OLD(IJCNT) = FREST(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END IF
              END DO
c
          END DO
c
c
      ELSE IF (K_PRPFX.EQ.7) THEN
          STOP 'prp_fx : NOT TESTET'
c     atoms specified in atom-list iatdyn mobile in restr geometry
c     all other atom without restriction
c                              no special symmetry
c     nfatom_dyn here has the meaning 'number of atoms mobile in restr'
c
          IF (NFATOM_DYN.GT.NFATOM) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn .gt. nfatom'
              STOP 'nfatom_dyn in prp_fx 7'
          END IF
c
          NBRYD = 3*NFATOM
c
c---> print (and sort) projection vectors
c
          DO JDYN = 1,NFATOM_DYN
              IAT = IATDYN(JDYN)
c
c ---> find equivalent atom in list ifatom
c
              JCNT = 0
              DO IDYN = 1,NFATOM
                  IF (IFATOM(IDYN).EQ.IAT) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no or wrong corresponding atom found in ifatom'
                  STOP 'iatdyn in prp_fx 7'
              END IF
              WRITE (IOFILE,FMT='(a35)')
     +          ' hypervector perpend to hyperplane '
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JDYN),ICNT=1,3),JDYN
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JDYN)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JDYN,IAT
c
c---> (1) calculate square norm of projection vector
c
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (2) calculate parallel components
c
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
              DO ICNT = 1,3
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
c
c---> (3) calculate force and coord restricted to hyperplane
c
              DO ICNT = 1,3
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
c
          END DO
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (perp to proj vector)
                  IF (IDYN.EQ.1) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = FREST(ICNT,JCNT)
                          XKS_M(IJCNT) = TREST(ICNT,IAT)
                      END DO
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                          XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                      TREST(ICNT,IAT)
                          TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                       DEL_TAUXYZ(ICNT,IAT)
                          F_OLD(IJCNT) = FREST(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
                      END DO
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                      TAUXYZ(ICNT,IAT)
                          TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                          F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
                      END DO
                  END IF
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 7'
          END IF
c
      ELSE IF (K_PRPFX.EQ.8) THEN
          STOP 'prp_fx : NOT TESTET'
c     all atoms mobile in restricted geometry;
c                                       no special symmetry
c
          IF (NFATOM_DYN.NE.NFATOM) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn .ne. nfatom'
              STOP 'nfatom_dyn in prp_fx 8'
          END IF
c
          NBRYD = 3*NFATOM_DYN
c
c---> print (and sort) incomimg projection vector
c
          WRITE (IOFILE,FMT='(a35)')
     +      ' hypervector perpend to hyperplane '
          DO JCNT = 1,NFATOM_DYN
              IAT = IFATOM(JCNT)
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JCNT),ICNT=1,3),JCNT
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JCNT)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JCNT,IAT
c
c---> (1) calculate square norm of projection vector
c
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (2) calculate parallel components
c
          DO JCNT = 1,NFATOM
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
          END DO
          DO ICNT = 1,3
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
          END DO
          DO JCNT = 1,NFATOM
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
          END DO
c
c---> (3) calculate force and coord restricted to hyperplane
c
          DO ICNT = 1,3
              DO JCNT = 1,NFATOM
                  IAT = IFATOM(JCNT)
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
          DO JCNT = 1,NFATOM
              IAT = IFATOM(JCNT)
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
          END DO
c
c---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = FREST(ICNT,JCNT)
                      IAT = IFATOM(JCNT)
                      XKS_M(IJCNT) = TREST(ICNT,IAT)
                  END DO
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TREST(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                   DEL_TAUXYZ(ICNT,IAT)
                      F_OLD(IJCNT) = FREST(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END DO
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 8'
          END IF
c
      ELSE IF (K_PRPFX.EQ.11) THEN
          STOP 'prp_fx : NOT TESTET'
c     atoms in upper half of slab (for which force is calc) are mobile;
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
          IF (NFATOM_DYN.NE.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn.ne.nfatom/2'
              STOP 'nfatom_dyn in prp_fx 11'
          END IF
c
          NBRYD = 3*NFATOM_DYN
c
c---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      IAT = IFATOM(JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                  END DO
              END DO
C
C ---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TAUXYZ(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END DO
              END DO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
C
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IAT = IFATOM(JCNT)
c ---> here we assume, that inverse images have consecutive numbers!!!
                      JAT = IAT + 1
                      TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                   DEL_TAUXYZ(ICNT,IAT)
                  END DO
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 11'
          END IF
c
      ELSE IF (K_PRPFX.EQ.12) THEN
          STOP 'prp_fx : NOT TESTET'
c     atoms specified in atom-list iatdyn mobile in z-dir;
c     all other atoms of upper half slab without restriction
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
c     nfatom_dyn here has the meaning 'number of atoms mobile in z dir'
c
          IF (NFATOM_DYN.GT.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          ' nfatom_dyn .gt. nfatom/2'
              STOP 'nfatom_dyn in prp_fx 12'
          END IF
c
          NBRYD = 3*NFATOM/2 - 2*NFATOM_DYN
c
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      ICNT = 3
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                          XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
c ---> here we assume, that inverse images have consecutive numbers!!!
                  JAT = IAT + 1
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      ICNT = 3
                      IJCNT = IJCNT + 1
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TAUXYZ(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                      TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                   DEL_TAUXYZ(ICNT,IAT)
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                      TAUXYZ(ICNT,IAT)
                          TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                          F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                          TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                       DEL_TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 12'
          END IF
c
      ELSE IF (K_PRPFX.EQ.13) THEN
          STOP 'prp_fx : NOT TESTET'
c     only atoms out of atom-list iatdyn mobile
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
          IF (NFATOM_DYN.GT.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn.gt.nfatom/2'
              STOP 'nfatom_dyn in prp_fx 13'
          END IF
          NBRYD = 3*NFATOM_DYN
          IJCNT = 0
          DO JDYN = 1,NFATOM_DYN
              JCNT = 0
              DO IDYN = 1,NFATOM,2
                  IF (IFATOM(IDYN).EQ.IATDYN(JDYN)) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0 .OR. MOD(JCNT,2).NE.1) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no or wrong corresponding atom found in ifatom'
                  STOP 'iatdyn in prp_fx 13'
              END IF
              IAT = IFATOM(JCNT)
c ---> here we assume, that inverse images have consecutive numbers!!!
              JAT = IAT + 1
              DO ICNT = 1,3
                  IJCNT = IJCNT + 1
c
c---> set position m
c
                  IF (K_SET.EQ.1) THEN
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
                  ELSE IF (K_SET.EQ.2) THEN
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TAUXYZ(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
C
C---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                      TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                   DEL_TAUXYZ(ICNT,IAT)
                  END IF
              END DO
c
          END DO
c
c
      ELSE IF (K_PRPFX.EQ.14) THEN
          STOP 'prp_fx : NOT TESTET'
c     only first atom (for which force is calc) mobile in z-direct
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
          NFATOM_DYN = 1
          NBRYD = 1
          JCNT = 1
          ICNT = 3
          IAT = IFATOM(JCNT)
c
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 1
              F_M(IJCNT) = F1XYZ(ICNT,JCNT)
              XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
C
C ---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 1
              DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) - TAUXYZ(ICNT,IAT)
              TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
              F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
              F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
C
c ---> here we assume, that inverse images have consecutive numbers!!!
              JAT = IAT + 1
              TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) - DEL_TAUXYZ(ICNT,IAT)
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx'
          END IF
c
      ELSE IF (K_PRPFX.EQ.16) THEN
          STOP 'prp_fx : NOT TESTET'
c     only atoms out of atom-list iatdyn mobile in restricted geometry;
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
          IF (NFATOM_DYN.GT.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn.gt.nfatom/2'
              STOP 'nfatom_dyn in prp_fx 16'
          END IF
c
          NBRYD = 3*NFATOM_DYN
c
c
c---> print (and sort) projection vector
c
          DO JDYN = 1,NFATOM_DYN
              IAT = IATDYN(JDYN)
              WRITE (IOFILE,FMT='(a35)')
     +          ' hypervector perpend to hyperplane '
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JDYN),ICNT=1,3),JDYN
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JDYN)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JDYN,IAT
c
c---> (1) calculate square norm of projection vector
c
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (2) calculate parallel components
c
          IJCNT = 0
          DO JDYN = 1,NFATOM_DYN
              JCNT = 0
c ---> find equivalent atom in list ifatom
              DO IDYN = 1,NFATOM,2
                  IF (IFATOM(IDYN).EQ.IATDYN(JDYN)) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0 .OR. MOD(JCNT,2).NE.1) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no or wrong corresponding atom found in ifatom'
                  STOP 'iatdyn in prp_fx 16'
              END IF
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
              DO ICNT = 1,3
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
c
c---> (3) calculate force and coord restricted to hyperplane
c
              DO ICNT = 1,3
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
c ---> here we assume, that inverse images have consecutive numbers!!!
              JAT = IAT + 1
              DO ICNT = 1,3
                  IJCNT = IJCNT + 1
c
c---> set position m
c
                  IF (K_SET.EQ.1) THEN
                      F_M(IJCNT) = FREST(ICNT,JCNT)
                      XKS_M(IJCNT) = TREST(ICNT,IAT)
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
                  ELSE IF (K_SET.EQ.2) THEN
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TREST(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                   DEL_TAUXYZ(ICNT,IAT)
                      F_OLD(IJCNT) = FREST(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                      TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                   DEL_TAUXYZ(ICNT,IAT)
                  END IF
              END DO
c
          END DO
c
c
      ELSE IF (K_PRPFX.EQ.17) THEN
          STOP 'prp_fx : NOT TESTET'
c     atoms specified in atom-list iatdyn mobile in restr geometry
c     all other atoms of upper half slab without restriction
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
c     nfatom_dyn here has the meaning 'number of atoms mobile in restr'
c
          IF (NFATOM_DYN.GT.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          ' nfatom_dyn .gt. nfatom/2'
              STOP 'nfatom_dyn in prp_fx 17'
          END IF
c
          NBRYD = 3*NFATOM/2
c
c---> print (and sort) projection vectors
c
          DO JDYN = 1,NFATOM_DYN
              IAT = IATDYN(JDYN)
c
c ---> find equivalent atom in list ifatom
c
              JCNT = 0
              DO IDYN = 1,NFATOM,2
                  IF (IFATOM(IDYN).EQ.IAT) JCNT = IDYN
              END DO
              IF (JCNT.EQ.0 .OR. MOD(JCNT,2).NE.1) THEN
                  WRITE (IOFILE,FMT=*) JDYN,IATDYN(JDYN),
     +              '  jdyn, iatdyn(jdyn);',
     +              '  no or wrong corresponding atom found in ifatom'
                  STOP 'iatdyn in prp_fx 17'
              END IF
              WRITE (IOFILE,FMT='(a35)')
     +          ' hypervector perpend to hyperplane '
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JDYN),ICNT=1,3),JDYN
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JDYN)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JDYN,IAT
c
c---> (1) calculate square norm of projection vector
c
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (2) calculate parallel components
c
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
              DO ICNT = 1,3
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
c
c---> (3) calculate force and coord restricted to hyperplane
c
              DO ICNT = 1,3
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
c
          END DO
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (perp to proj vector)
                  IF (IDYN.EQ.1) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = FREST(ICNT,JCNT)
                          XKS_M(IJCNT) = TREST(ICNT,IAT)
                      END DO
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                          XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
C
C---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
c ---> here we assume, that inverse images have consecutive numbers!!!
                  JAT = IAT + 1
                  IDYN = 3
                  DO JDYN = 1,NFATOM_DYN
                      IF (IFATOM(JCNT).EQ.IATDYN(JDYN)) IDYN = 1
                  END DO
c (a) for the atoms with restricted dynamics (parallel z dir)
                  IF (IDYN.EQ.1) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                      TREST(ICNT,IAT)
                          TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                       DEL_TAUXYZ(ICNT,IAT)
                          F_OLD(IJCNT) = FREST(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                          TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                       DEL_TAUXYZ(ICNT,IAT)
                      END DO
c (b) for the other atoms with unrestricted dynamics
                  ELSE IF (IDYN.EQ.3) THEN
                      DO ICNT = 1,3
                          IJCNT = IJCNT + 1
                          DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                      TAUXYZ(ICNT,IAT)
                          TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
                          F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                          F1XYZ(ICNT,JCNT) = ZERO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
c ---> here we assume, that inverse images have consecutive numbers!!!
C
                          TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                       DEL_TAUXYZ(ICNT,IAT)
                      END DO
                  END IF
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 17'
          END IF
c
      ELSE IF (K_PRPFX.EQ.18) THEN
          STOP 'prp_fx : NOT TESTET'
c     all atoms of upper half slab mobile in restricted geometry;
c     inversion symmetry (inverse image atoms have consecutive numbers)
c
          IF (NFATOM_DYN.NE.NFATOM/2) THEN
              WRITE (IOFILE,FMT=*) NFATOM_DYN,NFATOM,
     +          '  nfatom_dyn.ne.nfatom/2'
              STOP 'nfatom_dyn in prp_fx 18'
          END IF
c
          NBRYD = 3*NFATOM_DYN
c
c ---> print (and sort) incomimg fields
c
          WRITE (IOFILE,FMT='(a35)')
     +      ' hypervector perpend to hyperplane '
          DO JCNT = 1,NFATOM_DYN
c ---> here we assume, that inverse images have consecutive numbers!!!
              IAT = IFATOM(2*JCNT-1)
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),1(1x,i4),''       taupro'')')
     +          (TAUPRO(ICNT,JCNT),ICNT=1,3),JCNT
              DO ICNT = 1,3
                  TAUHLP(ICNT,IAT) = TAUPRO(ICNT,JCNT)
              END DO
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  tauhlp'')')
     +          (TAUHLP(ICNT,IAT),ICNT=1,3),JCNT,IAT
          END DO
c
c---> project forces and coordinates to plane perpendicular to tauhlp
c
c
c---> (1) calculate square norm of projection vector
c
          DO JCNT = 1,NFATOM,2
              IAT = IFATOM(JCNT)
              TAUHLPN(IAT) = ZERO
              DO ICNT = 1,3
                  TAUHLPN(IAT) = TAUHLPN(IAT) +
     +                           TAUHLP(ICNT,IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
c
c---> (2) calculate parallel components
c
          DO JCNT = 1,NFATOM,2
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = ZERO
              TPAR(IAT) = ZERO
          END DO
          DO ICNT = 1,3
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
                  FPAR(JCNT) = FPAR(JCNT) +
     +                         F1XYZ(ICNT,JCNT)*TAUHLP(ICNT,IAT)
                  TPAR(IAT) = TPAR(IAT) + TAUXYZ(ICNT,IAT)*
     +                        TAUHLP(ICNT,IAT)
              END DO
          END DO
          DO JCNT = 1,NFATOM,2
              IAT = IFATOM(JCNT)
              FPAR(JCNT) = FPAR(JCNT)/TAUHLPN(IAT)
              TPAR(IAT) = TPAR(IAT)/TAUHLPN(IAT)
          END DO
c
c---> (3) calculate force and coord restricted to hyperplane
c
          DO ICNT = 1,3
              DO JCNT = 1,NFATOM,2
                  IAT = IFATOM(JCNT)
                  FREST(ICNT,JCNT) = F1XYZ(ICNT,JCNT) -
     +                               FPAR(JCNT)*TAUHLP(ICNT,IAT)
                  TREST(ICNT,IAT) = TAUXYZ(ICNT,IAT) -
     +                              TPAR(IAT)*TAUHLP(ICNT,IAT)
              END DO
          END DO
          DO JCNT = 1,NFATOM,2
              IAT = IFATOM(JCNT)
              WRITE (IOFILE,FMT=
     +          '(3(d17.10,1x),2(1x,i4),''       frest '')')
     +          (FREST(ICNT,JCNT)*MRYD,ICNT=1,3),JCNT,IAT
              WRITE (IOFILE,FMT='(3(d17.10,1x),2(1x,i4),''  trest '')')
     +          (TREST(ICNT,IAT),ICNT=1,3),JCNT,IAT
          END DO
c
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = FREST(ICNT,JCNT)
                      IAT = IFATOM(JCNT)
                      XKS_M(IJCNT) = TREST(ICNT,IAT)
                  END DO
              END DO
C
C ---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
                      DEL_TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT) -
     +                                       TREST(ICNT,IAT)
                      TAUXYZ(ICNT,IAT) = TAUXYZ(ICNT,IAT) +
     +                                   DEL_TAUXYZ(ICNT,IAT)
                      F_OLD(IJCNT) = FREST(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END DO
              END DO
C
C ---> couple COORDINATES of lower surface atoms to upper surface
C
              IJCNT = 0
              DO ICNT = 1,3
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
c ---> here we assume, that inverse images have consecutive numbers!!!
                      JAT = IAT + 1
                      TAUXYZ(ICNT,JAT) = TAUXYZ(ICNT,JAT) -
     +                                   DEL_TAUXYZ(ICNT,IAT)
                  END DO
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx'
          END IF
c
      ELSE IF (K_PRPFX.EQ.21) THEN
          STOP 'prp_fx : NOT TESTET'
c     all atoms mobile;
c     full symmetry of point group considered
          WRITE (IOFILE,FMT=*) K_PRPFX,'   k_prpfx; full sym not implem'
          STOP 'k_prpfx in prp_fx 21'
c
      ELSE IF (K_PRPFX.EQ.135) THEN
          STOP 'prp_fx : NOT TESTET'
c     all atoms mobile;
c     special symmetry for relaxation of (110)-surf of
c     3-5-semiconductors:
c       mirror plane perpendicular to (long) x-dir (= surf normal)
c       no relaxation in y direction
c ---> here we assume, that symmetry related atoms have consecutive
c      numbers!!!
c
c
c ---> set position m
c
          IF (K_SET.EQ.1) THEN
              NBRYD = NFATOM
              IJCNT = 0
              DO ICNT = 1,3,2
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      F_M(IJCNT) = F1XYZ(ICNT,JCNT)
                      IAT = IFATOM(JCNT)
                      XKS_M(IJCNT) = TAUXYZ(ICNT,IAT)
                  END DO
              END DO
c
C
C ---> RESET COORDINATES FOR NEXT ITERATION
C
          ELSE IF (K_SET.EQ.2) THEN
              IJCNT = 0
              DO ICNT = 1,3,2
                  DO JCNT = 1,NFATOM,2
                      IJCNT = IJCNT + 1
                      IAT = IFATOM(JCNT)
                      TAUXYZ(ICNT,IAT) = XKS_MP1(IJCNT)
c ---> symmetry related atoms have consecutive numbers!!!
                      IF (ICNT.EQ.1) TAUXYZ(ICNT,IAT+1) = -TAUXYZ(ICNT,
     +                    IAT)
                      IF (ICNT.EQ.3) TAUXYZ(ICNT,IAT+1) = TAUXYZ(ICNT,
     +                    IAT)
                      F_OLD(IJCNT) = F1XYZ(ICNT,JCNT)
                      F1XYZ(ICNT,JCNT) = ZERO
                  END DO
              END DO
c
          ELSE
              WRITE (IOFILE,FMT=*) K_SET,'   k_set; wrong choice'
              STOP 'k_set in prp_fx 135'
          END IF
c
      ELSE
          WRITE (IOFILE,FMT=*) K_PRPFX,'   k_prpfx; wrong choice'
          STOP 'k_prpfx in prp_fx else'
      END IF
      RETURN
C
C---> FORMAT STATEMENTS
C
 2900 FORMAT (1x,'assign forces and coord to dynamical variables')
      END
