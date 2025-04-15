*DECK grad
      SUBROUTINE GRAD(IOFILE,KMDYN,K_ADPT_DT,ITER,MBRYD,NBRYD,XKS_M,F_M,
     +                F_OLD,STEPLN,XKS_MP1)
c-----------------------------------------------------------------------
c $Id: grad.f,v 1.2 1998/12/07 17:50:39 wikromen Exp $
C***********************************************************************
C      this subroutine performs a simple step
c      in direction f_m of steplengths f_m*stepln
c
C      xks_mp1  =xks_m + f_m*stepln
c
c      stepln is adapted dynamically, if kmdyn=1, k_adpt_dt=1
c
C
C***********************************************************************
c $Log: grad.f,v $
c Revision 1.2  1998/12/07  17:50:39  wikromen
c  ID and LOG added;   print ID included.
c
c-----------------------------------------------------------------------
C
      IMPLICIT NONE
C
C
C---> FILE NUMBER FOR READ AND WRITE
C
      INTEGER IOFILE
C
C---> DYNAMICS  INFORMATION
C
      INTEGER KMDYN,K_ADPT_DT
      INTEGER NBRYD,MBRYD,ITER
      REAL*8 STEPLN
      REAL*8 XKS_M(MBRYD),F_OLD(MBRYD),XKS_MP1(MBRYD),F_M(MBRYD)
C
C---> LOCAL VARIABLES
C
      REAL*8 EONE,TOL,KRIT
      INTEGER ICNT,I_MAX
      REAL*8 ZERO,ONE,TWO
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,TWO=2.0d0)
C
C---> INTRINSIC FUNCTIONS
C
      INTRINSIC EXP,SQRT,MAX,ABS
c----------------------------------------------------------------
C
C---> ABBREVIATIONS
C
C     IOFILE  : FILE NUMBER WHICH CONTAINS THE OUTPUT
C     ITER    : NUMBER OF ITERATIONS FOR UPDATING ionic coordinates
C     STEPLN  : adapted steplength for steepest descent step;
C                                  scaled with maximal force
C
c----------------------------------------------------------------
c
      LOGICAL FIRST_CALL

      FIRST_CALL = .true.
      IF (FIRST_CALL) THEN
          FIRST_CALL = .false.
          WRITE (IOFILE,FMT=*) '$RCSfile: grad.f,v $ $Revision: 1.2 $'
      END IF
c

      EONE = TWO
c
      WRITE (IOFILE,FMT=*) STEPLN,' stepln,transferred from mdyncs'
c
c---> dynamical set of step lenght stepln
c          (actually applied for kmdyn=1, k_adpt_dt=1)
c
      TOL = ZERO
      I_MAX = 1
      DO ICNT = 1,NBRYD
          KRIT = F_M(ICNT)*F_OLD(ICNT)
          KRIT = ABS(KRIT)
          IF (KRIT.GT.TOL) THEN
              TOL = KRIT
              I_MAX = ICNT
          END IF
      END DO

      KRIT = F_M(I_MAX)*F_OLD(I_MAX)
      WRITE (IOFILE,FMT=*) 'k_adpt_dt,i_max,krit:',K_ADPT_DT,I_MAX,KRIT
c
      IF (KMDYN.EQ.1) THEN
c

          IF (K_ADPT_DT.EQ.1 .AND. KRIT.LT.ZERO .AND.
     +        ABS(F_M(I_MAX)).GT.ABS(F_OLD(I_MAX))/TWO) THEN
              STEPLN = STEPLN/EONE
              WRITE (IOFILE,FMT=*) ITER,' iter; stepln lowered:     ',
     +          STEPLN
          ELSE IF (K_ADPT_DT.EQ.1 .AND. KRIT.GT.ZERO .AND.
     +             ABS(F_M(I_MAX)).GT.ABS(F_OLD(I_MAX))/TWO) THEN

              STEPLN = STEPLN*SQRT(EONE)
              WRITE (IOFILE,FMT=*) ITER,' iter; stepln increased:   ',
     +          STEPLN
          ELSE IF (K_ADPT_DT.EQ.1) THEN
              WRITE (IOFILE,FMT=*) ITER,' iter; stepln unchanged:   ',
     +          STEPLN
          END IF
c
      END IF

      IF (K_ADPT_DT.EQ.0 .AND. KRIT.LT.ZERO .AND.
     +    ABS(F_M(I_MAX)).GT.ABS(F_OLD(I_MAX))/TWO) THEN
          WRITE (IOFILE,FMT=*) ITER,' iter; stepln lower guess: ',STEPLN
      ELSE IF (K_ADPT_DT.EQ.0 .AND. KRIT.GT.ZERO .AND.
     +         ABS(F_M(I_MAX)).GT.ABS(F_OLD(I_MAX))/TWO) THEN
          WRITE (IOFILE,FMT=*) ITER,' iter; stepln incr. guess: ',STEPLN
      ELSE IF (K_ADPT_DT.EQ.0) THEN
          WRITE (IOFILE,FMT=*) ITER,' iter; stepln unch. guess: ',STEPLN
      END IF

C
C---> new coordinates
C
      DO ICNT = 1,NBRYD
          XKS_MP1(ICNT) = F_M(ICNT)*STEPLN + XKS_M(ICNT)
      END DO

      RETURN
      END
