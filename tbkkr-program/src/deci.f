c ***********************************************************************
      SUBROUTINE DECI(ML,M0,MR,TAU,
     +                QBOUND,DIM,AUX,NAUX,
     +                IPIV,IOPT,IMAX,IINFO)
      implicit none
c ***********************************************************************
c     on input:
c             ML    M0    MR 
c     LEFT    M_10  M_00  M_01
c     RIGHT   M_01  M_00  M_10
c
c     QBOUND upper bound for norm of alpha matrix
c     dimension of block matrices
c     AUX   auxiliarry array of dimension NAUX (minimum 8*DIM**2)
c     IPIV  integer array of at least dimension DIM
c     IOPT  not used (see DECI1)
c     IMAX  maximum number of iterations until criterion 
c           of QBOUND is reached (IMAX.GT.0)
c
c     on output :
c     TAU        Surface Green Function 
c     IINFO > 0  iteration converged
c                AUX(1) F-norm of alpha
c                AUX(2) F-norm of beta
c     IINFO < 0  iteration not converged,
c ************************************************************************
      INTEGER DIM,NAUX,IMAX,IOPT,IINFO,IPIV(*)
      DOUBLE PRECISION QBOUND
C
      DOUBLE COMPLEX 
     +     ML(DIM,*),
     +     M0(DIM,*),
     +     MR(DIM,*),
     +     TAU(DIM,*),
     +     AUX(DIM,DIM,*)
C
      INTEGER DIMM,IC,ILOOP,LM1,INFO
      DOUBLE PRECISION 
     +     BOUND,
     +     NORMFA,NORMFB
      logical select(1000)
C
      DOUBLE COMPLEX CONE,CZERO
      PARAMETER(CONE=(1.D0,0.D0),CZERO=(0.D0,0.D0))
C
      DOUBLE PRECISION ZNORM2
      EXTERNAL ZGEINV1,ZCOPY,ZNORM2,ZGEMM,ZGEEV
      INTRINSIC MAX
c ************************************************************************
      IF (NAUX.LT.8*DIM*DIM) THEN
        write(6,*) 'Increase NAUX in subroutine DECI.'
        STOP 'DECI - 1'
      END IF
      IF (IMAX.LE.0) THEN
        write(6,*) 'Set parameter IMAX to a nonzero value.'
        STOP 'DECI - 2'
      END IF
      DIMM = DIM*DIM
c
      IC = 0
c ------------------------------------------------------------------------

      CALL CINIT(NAUX,AUX)
c     aux(1) = alpha  ( = - M_10 )
      CALL ZAXPY(DIMM,-CONE,ML,1,AUX(1,1,1),1)
c     aux(2) = beta   ( = - M_01 )
      CALL ZAXPY(DIMM,-CONE,MR,1,AUX(1,1,2),1)
c     aux(3) = eps    ( = M_00 )
      CALL ZCOPY(DIMM,M0,1,AUX(1,1,3),1)
c     aux(4) = eps_s  ( = M_00 )
      CALL ZCOPY(DIMM,M0,1,AUX(1,1,4),1)
c
c ------------------------------------------------------------------------
      NORMFA = ZNORM2(DIMM,AUX(1,1,1),1)
      NORMFB = ZNORM2(DIMM,AUX(1,1,2),1)
c ------------------------------------------------------------------------
c
c ---> iteration loop
c
      BOUND = MIN(NORMFA,NORMFB)
c
      DO 20 WHILE ( BOUND.GT.QBOUND .AND. IC.LT.IMAX  )
c
c           aux(5) = eps**(-1) 
c
        CALL CINIT(DIM*DIM,AUX(1,1,5))
        DO 11 LM1=1,DIM
          AUX(LM1,LM1,5) = CONE
 11     END DO
        CALL ZCOPY(DIM*DIM,AUX(1,1,3),1,AUX(1,1,8),1)
        CALL ZGETRF(DIM,DIM,AUX(1,1,8),DIM,IPIV,INFO)
        CALL ZGETRS('N',DIM,DIM,AUX(1,1,8),DIM,IPIV,AUX(1,1,5),DIM,INFO)
c
c           aux(6) = eps**(-1) * alpha
c
        CALL ZGEMM('N','N',DIM,DIM,DIM,CONE,AUX(1,1,5),DIM,
     +       AUX(1,1,1),DIM,CZERO,AUX(1,1,6),DIM)
c
c           aux(7) = eps**(-1) * beta
c
        call zgemm('n','n',dim,dim,dim,cone,aux(1,1,5),dim,
     +       aux(1,1,2),dim,czero,aux(1,1,7),dim)
c
c           eps_s = eps_s - alpha * eps**(-1) * beta
c
        call zgemm('n','n',dim,dim,dim,-cone,aux(1,1,1),dim,
     +       aux(1,1,7),dim,cone,aux(1,1,4),dim)
c
c           eps = eps - alpha * eps**(-1) * beta
c
        CALL ZGEMM('N','N',DIM,DIM,DIM,-CONE,AUX(1,1,1),DIM,
     +       AUX(1,1,7),DIM,CONE,AUX(1,1,3),DIM)
c     
c           eps = eps - beta * eps**(-1) * alpha
c
        CALL ZGEMM('N','N',DIM,DIM,DIM,-CONE,AUX(1,1,2),DIM,
     +       AUX(1,1,6),DIM,CONE,AUX(1,1,3),DIM)
c
c           alpha = alpha * eps**(-1) * alpha
c
        CALL ZCOPY(DIMM,AUX(1,1,1),1,AUX(1,1,8),1)
        CALL ZGEMM('N','N',DIM,DIM,DIM,CONE,AUX(1,1,8),DIM,
     +       AUX(1,1,6),DIM,CZERO,AUX(1,1,1),DIM)
C
c           beta = beta * eps**(-1) * beta
c
        CALL ZCOPY(DIMM,AUX(1,1,2),1,AUX(1,1,8),1)
        CALL ZGEMM('N','N',DIM,DIM,DIM,CONE,AUX(1,1,8),DIM,
     +       AUX(1,1,7),DIM,CZERO,AUX(1,1,2),DIM)
c
        IC = IC + 1
c
        NORMFA = ZNORM2(DIMM,AUX(1,1,1),1)
        NORMFB = ZNORM2(DIMM,AUX(1,1,2),1)
c
        BOUND = MIN(NORMFA,NORMFB)
c
 20   END DO      ! WHILE (BOUND.GT.QBOUND .AND. ...
c
      IINFO = IC
      AUX(1,1,1) = NORMFA
      AUX(2,1,1) = NORMFB
      IF (BOUND.GT.QBOUND) THEN
        IINFO = -1
        RETURN
      END IF
c
c --->    TAU = eps_s**(-1)
c
      CALL CINIT(DIM*DIM,TAU)
      DO 12 LM1=1,DIM
        TAU(LM1,LM1) = CONE
 12   END DO
      CALL ZGETRF(DIM,DIM,AUX(1,1,4),DIM,IPIV,INFO)
      CALL ZGETRS('N',DIM,DIM,AUX(1,1,4),DIM,IPIV,TAU,DIM,INFO)
c
      RETURN
      END
