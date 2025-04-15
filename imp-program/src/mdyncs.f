*DECK mdyncs
      SUBROUTINE MDYNCS(IOBROY,IOBROG,IOFILE,KMDYN,KPRI,KTEST,K_PRPFX,
     +                  TAUPRO,IATDYN,NFATOM_DYN,ITER,IFATOM,ITDMD,
     +                  STEPL,DELTAMAX,T_DEBYE,NFATOM,DIMASS,EPS_MD,
     +                  K_EPS_MD,K_ADPT_DT,K_SET_F_M,TAU,TAUXYZ,F1XYZ,
     +                  MIT,F_OLD,STEPLN)
      IMPLICIT NONE
c---------------------------------------------------------------------
c $Id: mdyncs.f,v 1.2 1998/12/09 08:30:22 wikromen Exp $
C---------------------------------------------------------------------
C     TREE STRUCTURE:
C
C                 MDYNCS--
C                        |----prp_fx
C                        |----verlet-------random
C                        |----grad
C                        |----qnewton------saxpy
C                                      |---sscal
C
C---------------------------------------------------------------------
C
C     This subroutine calculates the new positions and velocities
C     of the mobile particles using
c
C    --  STEEPEST DESCENT (KMDYN = 1)
c    --  broyden methods  (KMDYN = 2-6)
C    --  Davidon-Fletcher-Powell rank 2 quasi-Newton algorithm (KMDYN=7)
C    --  Broyden-Fletcher-Goldfarb-Shanno rank 2 quasi-Newton algorithm
C                                                         (KMDYN = 8)
C
C     Extended to contain quasi_newton methods FOR MINIMUM SEARCH
C
C                                       K. SCHROEDER, IFF, JULY 1992
c     debugged
c                                       B. ENGELS, IFF, FEB 1995
c
c     extended
c        to find estimate for largest displacement which
c        scales with m_e/M_at * (debye_temp)**(-2):
c
c     del_x = alpha * f .lt. deltamax (e.g.0.15 a.u.)
c
c     alpha = Min{300*(m_e/M_at)*(300/debye_temp)**2; deltamax/max_forc}
c
c     The idea was communicated by M.Weinert
c
c     extended (see also qnewton_md.f) to ensure
c                  positive definite Hesse matrix in qnewton_md.f
c                                       K. SCHROEDER, IFF, JULY 1995
c     changes in 'set position m'
c            and 'RESET COORDINATES FOR NEXT ITERATION'
c     assignment and reassignment of forces and coordinates done in
c            subroutine prp_fx.f
c     keys k_prpfx (for symmetry use and # of dof)
c      and k_set (for assign and reassign, resp) introduced
c
c     ktest.ge.4: require local force to converge to max_fxyz/ten
c                                                 of current it_md step
c                            Kurt Schroeder, IFF, KFA Juelich, Oct. 1996

C-----------------------------------------------------------------------
c $Log: mdyncs.f,v $
c Revision 1.2  1998/12/09  08:30:22  wikromen
c  ID and LOG added; print ID included;
c
c-----------------------------------------------------------------------
C
C
C     .. Parameters ..
      INTEGER NATOMD
      PARAMETER (NATOMD=102)
      INTEGER MBRYD_MD
      PARAMETER (MBRYD_MD=3*NATOMD)
      REAL*8 ZERO,FOUR,ONE
      PARAMETER (ZERO=0.0d0,FOUR=4.0d0,ONE=1.0d0)
      REAL*8 TWO,HALF
      PARAMETER (TWO=2.0d0,HALF=0.5d0)
      REAL*8 MRYD
      PARAMETER (MRYD=1.d3)
      REAL*8 THREE_2,FIVE_2
      PARAMETER (THREE_2=3.0d2,FIVE_2=5.0d2)
      REAL*8 TOL_8
      PARAMETER (TOL_8=1.0d-8)
C     ..
C     .. Scalar Arguments ..
      REAL*8 DELTAMAX,EPS_MD,STEPL,STEPLN,T_DEBYE
      INTEGER IOBROG,IOBROY,IOFILE,ITDMD,ITER,KMDYN,KPRI,KTEST,
     +        K_ADPT_DT,K_EPS_MD,K_PRPFX,K_SET_F_M,MIT,NFATOM,NFATOM_DYN
C     ..
C     .. Array Arguments ..
      REAL*8 DIMASS(NATOMD),F1XYZ(3,NATOMD),F_OLD(MBRYD_MD),
     +       TAU(3,NATOMD),TAUPRO(3,NATOMD),TAUXYZ(3,NATOMD)
      INTEGER IATDYN(NATOMD),IFATOM(NATOMD)
C     ..
C     .. Local Scalars ..
      REAL*8 DELTA,MAX_FXYZ,MAX_MASS,MAX_TXYZ,PI,RMS_FXYZ,RMS_TXYZ,TPINV
      INTEGER IAT,ICNT,IFA,IFMAX,IJCNT,IR,ITMAX,I_QN_FAIL,K_SET,K_STOP,
     +        MBRYD,MITDEP,NATOM,NBRYD
      LOGICAL FIRST_CALL
C     ..
C     .. Local Arrays ..
      REAL*8 AM(NATOMD-1),BM(NATOMD-1),FM(MBRYD_MD),FM1(MBRYD_MD),
     +       F_M(MBRYD_MD),G(MBRYD_MD),SM(MBRYD_MD),SM1(MBRYD_MD),
     +       UI(MBRYD_MD,3),VI(MBRYD_MD,3),WIT(NATOMD-1),
     +       XKS_M(MBRYD_MD),XKS_MP1(MBRYD_MD),XKS_I(MBRYD_MD)
C     ..
C     .. External Subroutines ..
      EXTERNAL GRAD,PRP_FX,QNEWTON_MD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,MAX,MOD,SQRT
C     ..
      FIRST_CALL = .true.
      IF (FIRST_CALL) THEN
          FIRST_CALL = .false.
          WRITE (IOFILE,FMT=*) '$RCSfile: mdyncs.f,v $ $Revision: 1.2 $'
      END IF
c
      MBRYD = 3*NATOMD
      MITDEP = ITDMD
C
C---> ABBREVIATIONS
C
C      KMDYN  : RUNNING MODE PARAMETER GOVERNING TYPE OF
C               ALGORITHM USED IN MOLECULAR DYNAMICS
C      KPRI   : RUNNING MODE PARAMETER GOVERNING PRINT INFORMATION
C      KTEST  : RUNNING MODE PARAMETER GOVERNING TEST INFORMATION
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
c                  and  all other atoms mobile without restriction
c              17  atoms out of atom-list iatdyn mobile in restr geom
c                  and atoms in upper half of slab mobile plus inv
c                  images
c               8  all atoms mobile in restr geometry
c              18  all atoms mobile in restr geometry plus inv images
c             135  special for 110-surf of 3-5-comp semiconductors:
c                        atoms in upper half of slab mobile
c                                           plus yz-plane-mirr images
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
C      TAU    : ATOMIC POSITIONS (basis vectors) WITH RESPECT TO BRAVAIS
C               LATTICE; RELATIVE COORDINATES USED. 0<TAU(i,IAT)<1
C               ARRAY IS CHANGED IN THE SUBROUTINE
C      TAUXYZ : CARTESIAN COORDINATES OF ALL ATOMS IN A.U.(1:3,1:NATOMD)
C               ARRAY IS CHANGED IN THE SUBROUTINE
C      F1XYZ  : CONTAINS THE X,Y,Z--COMPONENTS OF THE
C               TOTAL HELLMANN--FEYNMAN FORCE (EWALD-FORCE DUE TO ALL
C               BARE IONS AND FORCE DUE TO ELECTRON DISTRIBUTION)
C               OF ALL ATOMS FOR WHICH THE FORCE SHOULD BE CALCULATED.
C               F1XYZ IS CALCULATED IN THE SUBROUTINE FORCE
C      DIMASS : DYNAMICAL ION MASS IN A.U.
C      ITER   : ACTUAL NUMBER OF TIME STEPS FOR Ionic minimum search
C      IDTION : ACTUAL NUMBER OF TIME STEPS FOR IONIC MOTION
C      DTION  : LENGTH OF TIME STEP IN (Ry)**(-1)= 4.83776874*10-17 sec
C               FOR IONIC MOTION
C      (IST IM HAUPTPROGRAM UMGETAUSCHT WORDEN VON DT--> DTION)
C      ZUSAETZLICH IST DTEL FUER ELECTRONISCHE  FREIHEITSGRADE DEFINIERT
c
c      i_qn_fail : key to indicate negative definite JACOBIAN matrix:
c               i_qn_fail=0: matrix ok
c               i_qn_fail=1: matrix not ok; stop broyden, reset mit=1
C
C----------------------------------------------------------------------
c
      PI = FOUR*ATAN(ONE)
      TPINV = HALF/PI
      WRITE (IOFILE,FMT=*) '---> k_eps_md,eps_md:     ',K_EPS_MD,EPS_MD
      WRITE (IOFILE,FMT=*) '---> k_adpt_dt,k_set_f_m: ',K_ADPT_DT,
     +  K_SET_F_M

      IF (KPRI.GE.4) THEN
          WRITE (IOFILE,FMT='(/)')
          WRITE (IOFILE,FMT='(''       *<* mdyncs *>* '')')
          WRITE (IOFILE,FMT='(3x,'' ~~~~~~~~~~~~~~~~~~~~~~~~ '')')
          WRITE (IOFILE,FMT=2900)
          WRITE (IOFILE,FMT='(/)')
      END IF
c
c---> set key for performance of qnewton methods
c
      IF (ITER.EQ.1) I_QN_FAIL = 0
c
c
c---> set position m   changed for call of subroutine prp_fx
c
      K_SET = 1
c
      CALL PRP_FX(IOFILE,KPRI,K_PRPFX,K_SET,MBRYD,NATOM,NFATOM,IFATOM,
     +            TAUPRO,IATDYN,NFATOM_DYN,NBRYD,F1XYZ,TAUXYZ,F_M,XKS_M,
     +            F_OLD,XKS_MP1,TAU,XKS_I)
C
c
c---> check, if local force is converged
c            to yield meaningful total force to at least two digits
c            if YES: proceed with molecular dynamics
c            if NO : return to electronic selfconsistency
c                    WITHOUT atomic displacement step
c
      MAX_FXYZ = ZERO
      DO IJCNT = 1,NBRYD
          MAX_FXYZ = MAX(MAX_FXYZ,ABS(F_M(IJCNT)))
      END DO
c
      IF (KMDYN.GT.0 .AND. KMDYN.LT.9 .OR.
     +    KMDYN.GT.11 .AND. KMDYN.LT.19) THEN
C
C   ---> QUASI-NEWTON RANK 1 or 2 METHODS (use SYMMETRIC HESSIAN MATRIX)
C                                  (KMDYN = 2-8,11-18)
C        OR       STEEPEST DESCENT (KMDYN = 1  )
C
C   ---> (A) Coordinates changed along the direction of force ;
c        iter=1 : estimate step length with debye frequency
C        KMDYN=1: step length dynam. adjusted to STEPLN for k_adpt_dt=1
C        KMDYN=2-9,12-19: step length kept const,
c                   except when reset mit=1: then new estimate used
c
          IF (ITER.EQ.1 .OR. (KMDYN.EQ.1.AND.STEPLN.LT.TOL_8) .OR.
     +        (MIT.EQ.1.AND. (KMDYN.GT.1.AND.KMDYN.LT.9.OR.
     +        KMDYN.GT.11.AND.KMDYN.LT.19))) THEN
c
c---> first setting of step lenghts stepln
c
c
c ---> calculate maximal mass and force
C
              MAX_MASS = ZERO
              DO IFA = 1,NFATOM
                  IAT = IFATOM(IFA)
                  MAX_MASS = MAX(MAX_MASS,TWO*DIMASS(IAT))
              END DO
              MAX_FXYZ = ZERO
              DO ICNT = 1,NBRYD
                  MAX_FXYZ = MAX(MAX_FXYZ,ABS(F_M(ICNT)))
              END DO
              WRITE (IOFILE,FMT=*) MAX_MASS,' max_mass '
              WRITE (IOFILE,FMT=*) MAX_FXYZ,' max_fxyz '
              WRITE (IOFILE,FMT=*) T_DEBYE,' debye_temp '
c
c ---> following choice for steplength for first steepest descent step
c                   copied from m.weinert's bfgs.f 21.7.95
c
c--->       choose a reasonable first guess for scaling, but
c--->       limit displacement to a maximum deltamax (e.g. 0.25 a.u.)
c--->       (may need to be changed for different systems)
c--->       this choice is based on a Debye temperature of 330K;
c--->       modify as needed (by changing deltamax)
c
              STEPLN = (FIVE_2/MAX_MASS)* ((THREE_2/T_DEBYE)**2)
c
              IF (STEPLN*MAX_FXYZ.GT.DELTAMAX) THEN
                  WRITE (IOFILE,FMT=*) ITER,STEPLN,
     +              ' iter, stepln too large;',
     +              ' reset to ensure del_x .lt. deltamax'
                  STEPLN = DELTAMAX/MAX_FXYZ
                  WRITE (IOFILE,FMT=*) ITER,STEPLN,' iter,new stepln'
              ELSE
                  WRITE (IOFILE,FMT=*) ITER,STEPLN,
     +              ' iter,stepln in expected range'
              END IF
          END IF
C
          CALL GRAD(IOFILE,KMDYN,K_ADPT_DT,ITER,MBRYD,NBRYD,XKS_M,F_M,
     +              F_OLD,STEPLN,XKS_MP1)
c
c   ---> set rms_error
c
          RMS_TXYZ = ZERO
          MAX_TXYZ = ZERO
          DO IJCNT = 1,NBRYD
              DELTA = XKS_MP1(IJCNT) - XKS_M(IJCNT)
              RMS_TXYZ = RMS_TXYZ + DELTA**2
              MAX_TXYZ = MAX(MAX_TXYZ,ABS(DELTA))
          END DO
          RMS_TXYZ = SQRT(RMS_TXYZ)/NBRYD
          WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5)')
     +      'it_md =',ITER,'rms-error:','rms_dtauxy(i) (au)=',RMS_TXYZ
          WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5)')
     +      'it_md =',ITER,'rms-error:','max_dtauxy(i) (au)=',MAX_TXYZ
C
C        for KMDYN = 2-8,12-18
C   ---> additional coordinate changes USING RANK 2 QUASI-NEWTON
C        METHODS, DEPENDING ON KMDYN (SEE ABOVE)
C
          IF (KMDYN.GT.1 .AND. KMDYN.LT.9 .OR.
     +        KMDYN.GT.11 .AND. KMDYN.LT.19) THEN
c
  100         CONTINUE
c
              IF (K_ADPT_DT.EQ.1) THEN
                  WRITE (IOFILE,FMT=*) KMDYN,K_ADPT_DT,
     +              ' kmdyn,k_adpt_dt; this combination not allowed'
                  STOP 'k_adpt_dt in mdyncs'
              END IF
c
              IF (K_SET_F_M.EQ.3) THEN
                  DO IJCNT = 1,NBRYD
                      F_M(IJCNT) = STEPLN*F_M(IJCNT)
                  END DO
              ELSE
                  WRITE (IOFILE,FMT=*) K_SET_F_M,
     +             ' k_set_f_m; this value not allowed, use k_set_f_m=3'
                  STOP 'mdyncs: k_set_f_m'
              END IF
c
c      ---> read and update the iteration number of Broyden iteration:
c           mit
c
              IF (KMDYN.GE.10) THEN
                  KMDYN = MOD(KMDYN,10)
                  REWIND IOBROG
                  READ (IOBROG,END=22,ERR=22) MIT
                  REWIND IOBROG
                  WRITE (IOFILE,FMT=*)
     +             'old broyden jakobian is read in for first iteration'
                  WRITE (IOFILE,FMT=*) 'with mit = ',MIT
                  GO TO 21
   22             CONTINUE
                  WRITE (IOFILE,FMT=*) 'WARNING:'
                  WRITE (IOFILE,FMT=*) 'read in of old broyden jakobian'
     +              ,' failed for first iteration'
                  MIT = 1
   21             CONTINUE
              END IF

c
c      ---> quasi newton update step
c

              CALL QNEWTON_MD(IOBROY,IOBROG,IOFILE,KMDYN,KPRI,KTEST,
     +                        MITDEP,ITDMD,I_QN_FAIL,MIT,MBRYD,NBRYD,
     +                        STEPL,DELTAMAX,XKS_M,F_M,XKS_MP1,G,FM,FM1,
     +                        SM,SM1,UI,VI,WIT,AM,BM,XKS_I,T_DEBYE)
C
c   ---> check scalar product  <del_f,del_x>
c        (a) if <del_f,del_x> > 0: proceed with qnewton
c        (b) if <del_f,del_x> < 0: stop qnewton,
c                                  remove hesse_matrix memory
c                                  proceed with steepest descent step
c        is done inside qnewton!!!!!
c
              IF (I_QN_FAIL.EQ.1) THEN
c
                  WRITE (IOFILE,FMT=*)
     +              'WARNING: <del_f,del_x> > 0: i_qn_fail=1'
                  WRITE (IOFILE,FMT=*) I_QN_FAIL,ITER,ITDMD,MIT,
     +              ' i_qn_fail,iter,itdmd,mit;',
     +        ' STOP qnewton in md, PROCEED with steepest descent step '

                  MIT = 1
                  I_QN_FAIL = 0

c --->          qnewton failed, and mit=1 reset
c                 in this case the displacement del_x
c                 calculated in qnewton_md
c                 with the accumulated Hessian matrix cannot be used.
c --->          Start with new steepest descent step;
c                 limit del_x to deltamax for each component
C
c ---> determine new stepln
c
                  MAX_MASS = ZERO
                  DO IFA = 1,NFATOM
                      IAT = IFATOM(IFA)
                      MAX_MASS = MAX(MAX_MASS,TWO*DIMASS(IAT))
                  END DO
c
c --->   reset forces and calculate stepln
c
                  MAX_FXYZ = ZERO
                  IF (K_SET_F_M.EQ.3) THEN
c remember f_m has been scaled before :f_m(icnt)=stepln*f_m(icnt)
                      DO ICNT = 1,NBRYD
                          F_M(ICNT) = F_M(ICNT)/STEPLN
                          MAX_FXYZ = MAX(MAX_FXYZ,ABS(F_M(ICNT)))
                      END DO
                  ELSE
                      WRITE (IOFILE,FMT=*) K_SET_F_M,
     +             ' k_set_f_m; this value not allowed, use k_set_f_m=3'
                      STOP 'mdyncs: k_set_f_m'
                  END IF
c
                  WRITE (IOFILE,FMT=*) MAX_MASS,' max_mass '
                  WRITE (IOFILE,FMT=*) MAX_FXYZ*MRYD,
     +              ' max_fxyz(mRy/au) '
                  WRITE (IOFILE,FMT=*) T_DEBYE,' debye_temp '

                  STEPLN = (FIVE_2/MAX_MASS)* ((THREE_2/T_DEBYE)**2)

                  IF (STEPLN*MAX_FXYZ.GT.DELTAMAX) THEN
                      WRITE (IOFILE,FMT=*) ITER,STEPLN,
     +                  ' iter, stepln too large;',
     +                  ' reset to ensure del_x .lt. deltamax'
                      STEPLN = DELTAMAX/MAX_FXYZ
                      WRITE (IOFILE,FMT=*) ITER,STEPLN,
     +                  ' iter,new stepln'
                  ELSE
                      WRITE (IOFILE,FMT=*) ITER,STEPLN,
     +                  ' iter,stepln in expected range'
                  END IF
c
c --- > call qnewton again
c
                  GO TO 100
c
              END IF
c
          END IF

C
c   ---> test output of iterated coordinates
c
          IF (KTEST.GE.3) THEN
              WRITE (IOFILE,FMT=
     +'(1x,'' mdyncs '',4x,''tauxyz(1,ir)'',5x,                '' f1xyz(
     +1,ir)'')')
              WRITE (IOFILE,FMT=
     +'(1x,''ir='',i4,2x,'','',1pe14.7,'','',2x,               1pe14.7)'
     +          ) (IFATOM(IR),TAUXYZ(1,IFATOM(IR)),F1XYZ(1,IR),IR=1,
     +          NFATOM)
          END IF
C
C         -----------------------------------------
C         |                                        |
C   --->  |   END OF QUASI-NEWTON METHODS          |
C         |                                        |
C         ------------------------------------------
C
      END IF

C
C---> RESET COORDINATES FOR NEXT ITERATION
C
      K_SET = 2
c
      CALL PRP_FX(IOFILE,KPRI,K_PRPFX,K_SET,MBRYD,NATOM,NFATOM,IFATOM,
     +            TAUPRO,IATDYN,NFATOM_DYN,NBRYD,F1XYZ,TAUXYZ,F_M,XKS_M,
     +            F_OLD,XKS_MP1,TAU,XKS_I)
C
c
c---> set rms_error
c
      RMS_TXYZ = ZERO
      MAX_TXYZ = ZERO
      ITMAX = 0
      RMS_FXYZ = ZERO
      MAX_FXYZ = ZERO
      IFMAX = 0
      DO IJCNT = 1,NBRYD
          DELTA = XKS_MP1(IJCNT) - XKS_M(IJCNT)
          RMS_TXYZ = RMS_TXYZ + DELTA**2
          IF (MAX_TXYZ.LT.ABS(DELTA)) THEN
              ITMAX = IJCNT
              MAX_TXYZ = ABS(DELTA)
              WRITE (IOFILE,FMT=*) ITMAX,MAX_TXYZ,' itmax, max_txyz'
          END IF
          RMS_FXYZ = RMS_FXYZ + F_OLD(IJCNT)**2
          IF (MAX_FXYZ.LT.ABS(F_OLD(IJCNT))) THEN
              IFMAX = IJCNT
              MAX_FXYZ = ABS(F_OLD(IJCNT))
              WRITE (IOFILE,FMT=*) IFMAX,MAX_FXYZ,' ifmax, max_fxyz'
          END IF
      END DO
      RMS_TXYZ = SQRT(RMS_TXYZ)/NBRYD
      RMS_FXYZ = SQRT(RMS_FXYZ)/NBRYD
      WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5)') 'it_md =',
     +  ITER,'rms-error:','rms_dtauxy    (au)=',RMS_TXYZ
      WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5,2x,a4,i5)')
     +  'it_md =',ITER,'rms-error:','max_dtauxy    (au)=',MAX_TXYZ,
     +  'dof=',ITMAX
      WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5)') 'it_md =',
     +  ITER,'rms-error:','rms_fxyz  (mRy/au)=',RMS_FXYZ*MRYD
      WRITE (IOFILE,FMT='(5x,a7,i3,2x,a10,2x,a19,2x,d12.5,2x,a4,i5)')
     +  'it_md =',ITER,'rms-error:','max_fxyz  (mRy/au)=',MAX_FXYZ*MRYD,
     +  'dof=',IFMAX
c
c---> stop program
c
      K_STOP = 0
      IF (K_EPS_MD.EQ.1 .AND. RMS_TXYZ.LT.EPS_MD) THEN
          WRITE (IOFILE,FMT=*) '---> required md-convergence reached'
          WRITE (IOFILE,FMT=*) '     rms_txyz < ',EPS_MD
          K_STOP = 1
      END IF
      IF (K_EPS_MD.EQ.2 .AND. MAX_TXYZ.LT.EPS_MD) THEN
          WRITE (IOFILE,FMT=*) '---> required md-convergence reached'
          WRITE (IOFILE,FMT=*) '     max_txyz < ',EPS_MD
          K_STOP = 1
      END IF
      IF (K_EPS_MD.EQ.3 .AND. RMS_FXYZ*MRYD.LT.EPS_MD) THEN
          WRITE (IOFILE,FMT=*) '---> required md-convergence reached'
          WRITE (IOFILE,FMT=*) '     rms_fxyz < ',EPS_MD
          K_STOP = 1
      END IF
      IF (K_EPS_MD.EQ.4 .AND. MAX_FXYZ*MRYD.LT.EPS_MD) THEN
          WRITE (IOFILE,FMT=*) '---> required md-convergence reached'
          WRITE (IOFILE,FMT=*) '     max_fxyz < ',EPS_MD
          K_STOP = 1
      END IF
c
c
c
      DO IFA = 1,NFATOM
          IAT = IFATOM(IFA)
          WRITE (IOFILE,FMT=650) TAUXYZ(1,IAT),TAUXYZ(2,IAT),
     +      TAUXYZ(3,IAT),TAU(1,IAT),TAU(2,IAT),TAU(3,IAT)
      END DO
c
c---> stop
c
      IF (K_STOP.EQ.1) THEN
          WRITE (6,FMT=*) 'mdyncs: convergence'
          STOP
      END IF

      RETURN
C
C---> FORMAT STATEMENTS
C
  650 FORMAT (6F10.5)
 2900 FORMAT (1x,'calculate new positions and velocities',
     +       ' of the mobile particles ')
      END
