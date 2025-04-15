*DECK qnewton_md
      SUBROUTINE QNEWTON_MD(IOBROY,IOBROG,IOFILE,KMIX,KPRI,KTEST,
     +                      MITDEPT,NITDEPT,I_FAIL,ITDEPT,MDIM,NDIM,
     +                      ALPHA,DELTAMAX,XKS_M,F_M,XKS_MP1,G,FM,FM1,
     +                      SM,SM1,UI,VI,WIT,AM,BM,XKS_I,T_DEBYE)
      IMPLICIT NONE

c---------------------------- ----------------------------------------
c $Id: qnewton_md.F,v 1.2 1998/12/09 08:52:31 wikromen Exp $
c*********************************************************************
c     this subroutine performs the quasi newton update step
c     of XKS_M to XKS_MP1
c'    according to broyden's iteration schemes or improved versions
c     depending on the chosen key kmix:
c
c'      2    broyden's    f i r s t  m e t h o d               (ibroy=1)
c'      3    broyden's    s e c o n d  m e t h o d             (ibroy=2)
c'      4    anderson's   g e n e r a l i z e d   m e t h o d  (ibroy=3)
c       5    use ibroy=1 and save jacobian accumulated up to nitdept-1
c            to update xks
c       6    use ibroy=2 and save jacobian accumulated up to nitdept-1
c            to update xks
c       7    davidon, fletcher, powell  r a n k 2  m e t h o d (ibroy=6)
c       8    broyden, fletcher, goldfarb,shanno
c                                       r a n k 2  m e t h o d (ibroy=7)
c
c     methods 1 - 6
c     implemented here according to notes of s.b.
c     broyden iteration scheme following the papers of :
c     srivastava, j. phys. , 17 (1984) , pp l317
c     c.g. broyden in math.comput., 19 , pp 577, 1965
c     c.g. broyden in ibid, 21 ,pp 368 ,1967
c     the method has been generalized to include a metric. the
c     definition of the necessary inner products are similar to the
c     discription given in the notes of m.weinert. the algorithm
c     discribed in the paper of srivastava has been simplified
c     ( see notes of s.b.)
c
c'    broyden's update treats charge and spin on the same footing
c                                         s. bluegel , kfa , may 1987
c
c     the anderson method (d.g. anderson, j. acm 12, 547 (1964)) has
c'    been generalized and reformulated as an improvement of broyden's
c     second method. successive linesearch is replaced by successive
c     search on hyperplanes. ( see notes of s.b. )
c                                         s. bluegel , issp , july 1989
c
c     the files ui,vi are stored on high speed ssd memory with
c                                                      write(iobroy,....
c     the files fm,sm are stored on high speed ssd memory with
c                                                      write(iobrog,....
c
c     modified for potential mixing in real space
c     for use with plane wave representation of wave functions
c                                                k. schroeder, jan. 1992
c
c     methods 7 - 8
c     implemented according to notes of k. schroeder
c     basic formulae can be found in phd-thesis of s.bluegel and in
c     r.fletcher, methods in computation 2nd ed. 1985, chapter 3
c                                                k. schroeder, july 1992
c
c---> for fixpoint surch:  xks=f(xks)
c                          F(xks)=f(xks)-xks
c
c     input:  xks_m
c             f_m = F(xks_m)
c     output: xks_mp1
c
c---> for minimum search:   min E(xks) = ?
c                          f(xks) = - dE/dx
c
c     input:  xks_m
c             f_m =  f(xks_m)
c     output: xks_mp1
c             or i_fail=1, if negative definite {- Hesse matrix}
c                             cannot be guaranteed
c
c---> subroutine developed out of brydbm.f
c     brydbm was driver and update step
c     qnewton performes only update step and needs a driver routine
c                                       b.engels, jan 1995
c
c---> check for sign of scalar product <del_x,del_f> included
c     to ensure positiv definite Hessian matrix.
c     In our code this means NEGATIV definit Jacobian matrix.
c     For <del_x,del_f> > 0: set i_fail=1,
c                            stop qnewton, remove hessian memory,
c                            proceed with steepest descent step
c                                    k.schroeder, july 1995
c*********************************************************************
c $Log: qnewton_md.F,v $
c Revision 1.2  1998/12/09  08:52:31  wikromen
c ID and LOG added; print ID included;
c
c
c----------------------------------------------------------------------
C
C
c---> running mode parameters
c
c
C---> file numbers
c
C
C---> broyden dimensions
C
c
c---> arrays arguments
c

c
c---> local work arrays
c
c
c---> local scalars
c
c
c---> externals
c
c
c---> intrinsics
c

c*********************************************************************
c
c ---> PRINT REVISION-NUMBER
c
C     .. Parameters ..
      REAL*8 ZERO,ONE,TOL_70
      PARAMETER (ZERO=0.0d0,ONE=1.0d0,TOL_70=1.0d-70)
C     ..
C     .. Scalar Arguments ..
      REAL*8 ALPHA,DELTAMAX
      INTEGER IOBROG,IOBROY,IOFILE,ITDEPT,I_FAIL,KMIX,KPRI,KTEST,MDIM,
     +        MITDEPT,NDIM,NITDEPT,T_DEBYE
C     ..
C     .. Array Arguments ..
      REAL*8 AM(MITDEPT-1),BM(MITDEPT-1),FM(MDIM),FM1(MDIM),
     +                 F_M(MDIM),G(MDIM),SM(MDIM),SM1(MDIM),UI(MDIM,3),
     +                 VI(MDIM,3),WIT(MITDEPT-1),XKS_M(MDIM),
     +                 XKS_MP1(MDIM),XKS_I(MDIM)
C     ..
C     .. Local Scalars ..
      REAL*8 CI,CMM,SCAL_P,SMNORM,VMDENO,VMNORM
      INTEGER IBROY,IJ,IT
      LOGICAL FIRST_CALL
C     ..
C     .. External Functions ..
      REAL*8 DDOT
      EXTERNAL DDOT
C     ..
C     .. External Subroutines ..
      EXTERNAL DAXPY,DSCAL
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
      FIRST_CALL = .true.
      IF (FIRST_CALL) THEN
          FIRST_CALL = .false.
          WRITE (IOFILE,FMT=*)
     +      '$RCSfile: qnewton_md.F,v $ $Revision: 1.2 $'
      END IF
c

      I_FAIL = 0

      IF (NITDEPT.GT.MITDEPT) STOP 'qnewton: nitdept.gt.mitdept'
      IF (NDIM.GT.MDIM) STOP 'qnewton: ndim.gt.mdim'

      IBROY = KMIX - 1
      IF (IBROY.LE.0 .OR. IBROY.GT.7) STOP 'ibroyd'
      IF ((IBROY.LE.3.OR.IBROY.GE.7) .AND. ITDEPT.GT.NITDEPT) ITDEPT = 1
      IF ((IBROY.EQ.4.OR.IBROY.EQ.5) .AND.
     +    ITDEPT.GE.NITDEPT) ITDEPT = NITDEPT
      IF (IBROY.EQ.4 .AND. ITDEPT.LT.NITDEPT) IBROY = 1
      IF (IBROY.EQ.5 .AND. ITDEPT.LT.NITDEPT) IBROY = 2
      IF (IBROY.EQ.1) WRITE (IOFILE,FMT=
     +    '('' broyden"s 1st method used '')')
      IF (IBROY.EQ.2) WRITE (IOFILE,FMT=
     +    '('' broyden"s 2nd method used '')')
      IF (IBROY.EQ.3) WRITE (IOFILE,FMT=
     +    '('' generalized anderson method used '')')
      IF (IBROY.EQ.4) WRITE (IOFILE,FMT=
     +    '(''jacobian is fixed after broyden"s 1st method used'')')
      IF (IBROY.EQ.5) WRITE (IOFILE,FMT=
     +    '(''jacobian is fixed after broyden"s 2nd method used'')')
      IF (IBROY.EQ.6) WRITE (IOFILE,FMT=
     +    '(''davidon, fletcher, powell rank 2 method used '')')
      IF (IBROY.EQ.7) WRITE (IOFILE,FMT=
     +   '(''broyden, fletcher, goldfarb, shanno rank 2 method used '')'
     +    )
c
c---> for test purpose: print out nitdept and itdept
c
      WRITE (IOFILE,FMT='(1x,''nitdept= '',i4,'', itdept= '',i4)')
     +  NITDEPT,ITDEPT

      IF (ITDEPT.EQ.1) THEN
c
c   ---> set up : sm1 = xks(1) ; fm1=fm[1]=f(xks(1)) - xks(1) ;
c                 metric  g := 1
c
          DO IJ = 1,NDIM
              SM1(IJ) = XKS_M(IJ)
              FM1(IJ) = F_M(IJ)
              G(IJ) = ONE
              XKS_MP1(IJ) = ALPHA*FM1(IJ) + XKS_M(IJ)
          END DO

c
c   ---> test output of iterated potential
c
          IF (KTEST.GE.3) THEN
              WRITE (IOFILE,FMT=
     +'(1x,/,'' test metric and last iteration fields '',
     +/,14x,'' g '',15x,''sm1'',14x,''fm1'')')
              WRITE (IOFILE,FMT=
     +'(1x,''ij='',i4,2x,1pe14.7,'','',2x,1pe14.7,'','',
     +                   2x,1pe14.7)') (IJ,G(IJ),SM1(IJ),FM1(IJ),IJ=2,
     +          10,2)
          END IF
c
c   ---> store on ssd metric g and last iteration fields sm1 and fm1
c
c
          DO IJ = 2,MITDEPT - 1
              WIT(IJ) = ZERO
              AM(IJ) = ZERO
              BM(IJ) = ZERO
          END DO
          REWIND IOBROG
          WRITE (IOBROG,ERR=2201) ITDEPT + 1
          WRITE (IOBROG,ERR=2202) (G(IJ),IJ=1,NDIM),
     +      (SM1(IJ),IJ=1,NDIM), (FM1(IJ),IJ=1,NDIM)
          WRITE (IOBROG,ERR=2203) (WIT(IJ),IJ=2,MITDEPT-1),
     +      (AM(IJ),IJ=2,MITDEPT-1), (BM(IJ),IJ=2,MITDEPT-1)
c
      ELSE
c
c   ---> read in metric g and last iteration fields sm1 and fm1
c
          REWIND IOBROG
          READ (IOBROG,ERR=2204) ITDEPT
          READ (IOBROG,ERR=2205) (G(IJ),IJ=1,NDIM), (SM1(IJ),IJ=1,NDIM),
     +      (FM1(IJ),IJ=1,NDIM)
          READ (IOBROG,ERR=2206) (WIT(IJ),IJ=2,MITDEPT-1),
     +      (AM(IJ),IJ=2,MITDEPT-1), (BM(IJ),IJ=2,MITDEPT-1)
c
c   ---> set up : sm = xks_m ; fm = f(xks_m) - xks_m ;
c
          DO IJ = 1,NDIM
              SM(IJ) = XKS_M(IJ)
              FM(IJ) = F_M(IJ)
          END DO
c
c   ---> calculate dsm = pot(m) - pot(m-1)
c   ---> calculate dfm = f[m] - f[m-1]
c
          DO IJ = 1,NDIM
              SM1(IJ) = SM(IJ) - SM1(IJ)
              FM1(IJ) = FM(IJ) - FM1(IJ)
              WRITE (110,FMT=*) SM1(IJ),FM1(IJ),SM1(IJ)*FM1(IJ)
          END DO
c
c   ---> check scalar product  <dsm,dfm>
c        We iterate the NEGATIV Hesse matrix:
c        (a) if <dsm,dfm> < 0:     Hesse matrix stays positive definite;
c                                  proceed with qnewton
c        (b) if <dsm,dfm> > 0:     set i_fail=1,
c                                  stop qnewton,
c                                  remove Hesse_matrix memory
c                                  proceed with steepest descent step
c
          SCAL_P = DDOT(NDIM,SM1,-1,FM1,-1)

          WRITE (IOFILE,FMT=*) NITDEPT,ITDEPT,SCAL_P,
     +      ' nitdept,itdept,scal_p in md'
c
          IF (SCAL_P.GT.ZERO) THEN
              I_FAIL = 1
              WRITE (IOFILE,FMT=*)
     +          'WARNING: <del_f,del_x> > 0: i_fail=1'
              WRITE (IOFILE,FMT=*) I_FAIL,NITDEPT,ITDEPT,SCAL_P,
     +          ' i_fail,nitdept,itdept,scal_p'
c    >       ' STOP qnewton in md, PROCEED with steepest descent step '
c           itdept = 1
          END IF
C
c
c   ---> branching for rank 1 methods (ibroy .le. 5)
c   --->           and rank 2 methods (ibroy .ge. 6)
c
          IF (IBROY.LE.5) THEN
c
c      ---> use rank 1 methods to update jacobian
c
c      ---> loop to generate u[m] = u(ij,itdept)
c           (remember: all v(ij,2) are convoluted with the metric g)
c
              DO IJ = 1,NDIM
                  UI(IJ,3) = ALPHA*FM1(IJ) + SM1(IJ)
              END DO
              REWIND IOBROY
              DO IT = 2,ITDEPT - 1
                  READ (IOBROY,ERR=2101) (UI(IJ,2),IJ=1,NDIM),
     +              (VI(IJ,2),IJ=1,NDIM)
                  AM(IT) = DDOT(NDIM,VI(1,2),-1,FM1,-1)
                  CALL DAXPY(NDIM,-AM(IT),UI(1,2),1,UI(1,3),1)
              END DO
c
c      ---> print amj = the importance of the history of ui
c
              WRITE (IOFILE,FMT=
     +          '(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +          ITDEPT - 1, (AM(IT),IT=2,ITDEPT-1)
c
c
              IF (IBROY.EQ.1) THEN
c'        ---> b r o y d e n ' s   f i r s t   m e t h o d
c
c         ---> calculate dsmnorm
c
                  SMNORM = ZERO
                  DO IJ = 1,NDIM
                      SMNORM = SMNORM + SM1(IJ)*G(IJ)*SM1(IJ)
                  END DO
c
c         ---> convolute dsm with the metric g
c
                  DO IJ = 1,NDIM
                      SM1(IJ) = G(IJ)*SM1(IJ)
                  END DO
c
c         ---> loop to generate v[m] = v(ij,itdept) ;
c              convoluted with the metric g
c
                  DO IJ = 1,NDIM
                      VI(IJ,3) = ALPHA*SM1(IJ)
                  END DO
                  REWIND IOBROY
                  DO IT = 2,ITDEPT - 1
                      READ (IOBROY,ERR=2102) (UI(IJ,2),IJ=1,NDIM),
     +                  (VI(IJ,2),IJ=1,NDIM)
                      BM(IT) = DDOT(NDIM,UI(1,2),-1,SM1,-1)
                      CALL DAXPY(NDIM,-BM(IT),VI(1,2),1,VI(1,3),1)
                  END DO
c
c         ---> complete the evaluation of v[m]
c
                  VMDENO = DDOT(NDIM,UI(1,3),-1,SM1,-1) - SMNORM
c
c         ---> test output of vmdeno
c
                  IF (KTEST.GE.3) THEN
                      WRITE (IOFILE,FMT=
     +'(1x,/,'' vmdeno for broydens 1st method (kmix=2)''
     +,e14.7)') VMDENO
                  END IF
                  IF (ABS(VMDENO).LT.TOL_70) STOP 'qnewton: vmdeno'

                  CALL DSCAL(NDIM,ONE/VMDENO,VI(1,3),1)
c
c         ---> print bmj = the importance of the history of vi
c
                  WRITE (IOFILE,FMT=
     +              '(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +              ITDEPT - 1, (BM(IT),IT=2,ITDEPT-1)
c
              ELSE IF (IBROY.EQ.2) THEN
c'        ---> b r o y d e n ' s   s e c o n d    m e t h o d
c
c         ---> calculate v[m] ; convoluted with the metric g
c
                  DO IJ = 1,NDIM
                      VI(IJ,3) = G(IJ)*FM1(IJ)
                  END DO
c
c         ---> calculate #vm# and normalize v[m]
c
                  VMNORM = DDOT(NDIM,VI(1,3),-1,FM1,-1)
c
c         ---> test output of vmnorm
c
                  IF (KTEST.GE.3) THEN
                      WRITE (IOFILE,FMT=
     +'(1x,/,'' vmnorm for broydens 2nd method (kmix=3)''
     +,e14.7)') VMNORM
                  END IF
                  CALL DSCAL(NDIM,ONE/VMNORM,VI(1,3),1)

              ELSE IF (IBROY.EQ.3) THEN
c         ---> g e n e r a l i z e d   a n d e r s o n   m e t h o d
c
c         ---> calculate v[m] ; convoluted with the metric g
c
                  DO IJ = 1,NDIM
                      VI(IJ,3) = G(IJ)*FM1(IJ)
                  END DO
                  REWIND IOBROY
                  DO IT = 2,ITDEPT - 1
                      READ (IOBROY,ERR=2103) (UI(IJ,2),IJ=1,NDIM),
     +                  (VI(IJ,2),IJ=1,NDIM)
                      CALL DAXPY(NDIM,-AM(IT)*WIT(IT),VI(1,2),1,VI(1,3),
     +                           1)
                  END DO
c
c         ---> complete the evaluation of v[m]
c
                  VMDENO = DDOT(NDIM,FM1,-1,VI(1,3),-1)
c
c         ---> test output of vmdeno
c
                  IF (KTEST.GE.3) THEN
                      WRITE (IOFILE,FMT=
     +'(1x,/,'' vmdeno for andersons method (kmix=4)'',
     +e14.7)') VMDENO
                  END IF
                  IF (ABS(VMDENO).LT.TOL_70) STOP 'qnewton: vmdeno'

                  CALL DSCAL(NDIM,ONE/VMDENO,VI(1,3),1)
c
c         ---> save wit(itdept) for next iteration
c
                  WIT(ITDEPT) = VMDENO
c
              END IF
c
c      ---> update f[m-1] <-- f[m]  ; pot(m-1) <-- pot(m)
c      ---> store on ssd metric g and last iteration fields sm and fm
c
              REWIND IOBROG
              WRITE (IOBROG,ERR=2207) ITDEPT + 1
              WRITE (IOBROG,ERR=2208) (G(IJ),IJ=1,NDIM),
     +          (SM(IJ),IJ=1,NDIM), (FM(IJ),IJ=1,NDIM)
              WRITE (IOBROG,ERR=2209) (WIT(IJ),IJ=2,MITDEPT-1),
     +          (AM(IJ),IJ=2,MITDEPT-1), (BM(IJ),IJ=2,MITDEPT-1)
c
c      ---> test output of iterated potential
c
              IF (KTEST.GE.3) THEN
                  WRITE (IOFILE,FMT=
     +'(1x,/,'' test metric and last iteration fields '',
     +/,14x,'' g '',15x,''sm'',14x,''fm'')')
                  WRITE (IOFILE,FMT=
     +'(1x,''ij='',i4,2x,1pe14.7,'','',2x,1pe14.7,'','',
     +2x,1pe14.7)') (IJ,G(IJ),SM(IJ),FM(IJ),IJ=2,10,2)
              END IF
c
c      ---> store on ssd u(ij,itdept) and v(ij,itdept)
c
              IF (IBROY.LE.3) THEN
                  WRITE (IOBROY,ERR=2104) (UI(IJ,3),IJ=1,NDIM),
     +              (VI(IJ,3),IJ=1,NDIM)
              END IF
c
c      ---> calculate cmm
c
              IF (IBROY.LE.3) THEN
                  CMM = DDOT(NDIM,VI(1,3),-1,FM,-1)
                  WRITE (IOFILE,FMT='(5x,'' cmm = '',1p,d12.4)') CMM
              ELSE
                  CMM = ZERO
              END IF
c
c      ---> update pot(m+1)
c
              CALL DAXPY(NDIM,ONE-CMM,UI(1,3),1,SM,1)
c
          ELSE IF (IBROY.GE.6) THEN
c
c      ---> use rank 2 methods to update jacobian
c
              AM(ITDEPT) = DDOT(NDIM,FM1,-1,SM1,-1)
c
c      ---> test output of um-norm
c
              IF (KTEST.GE.3) THEN
                  WRITE (IOFILE,FMT=
     +'(1x,/,'' um-norm for dfp rank 2 method (kmix=7)''
     +,e14.7)') AM(ITDEPT)
              END IF
              IF (ABS(AM(ITDEPT)).LT.TOL_70) STOP 'qnewton: am(itdept)'
              AM(ITDEPT) = ONE/AM(ITDEPT)
c
c      ---> loop to generate u[m] = u(ij,itdept)
c           and part of v[m] = v(ij,itdept)
c
              DO IJ = 1,NDIM
                  UI(IJ,3) = SM1(IJ)
                  VI(IJ,3) = ALPHA*FM1(IJ)
              END DO
c
              IF (IBROY.EQ.6) THEN
c
c         ---> davidon, fletcher, powell rank 2 method
c
c         ---> loop to complete generation of v[m] = v(ij,itdept)
c              in two steps
c
                  REWIND IOBROY
                  DO IT = 2,ITDEPT - 1
                      READ (IOBROY,ERR=2105) (UI(IJ,2),IJ=1,NDIM),
     +                  (VI(IJ,2),IJ=1,NDIM)

                      CI = AM(IT)*DDOT(NDIM,UI(1,2),-1,FM1,-1)
                      CALL DAXPY(NDIM,-CI,UI(1,2),1,VI(1,3),1)

                      CI = BM(IT)*DDOT(NDIM,VI(1,2),-1,FM1,-1)
                      CALL DAXPY(NDIM,-CI,VI(1,2),1,VI(1,3),1)
                  END DO

                  BM(ITDEPT) = DDOT(NDIM,VI(1,3),-1,FM1,-1)
c
c         ---> test output of vm-norm
c
                  IF (KTEST.GE.3) THEN
                      WRITE (IOFILE,FMT=
     +'(1x,/,'' vm-norm for dfp rank 2 method (kmix=7)''
     +,e14.7)') BM(ITDEPT)
                  END IF
                  IF (ABS(BM(ITDEPT)).LT.TOL_70)
     +                STOP 'qnewton: bm(itdept)'
                  BM(ITDEPT) = ONE/BM(ITDEPT)
c
c         ---> print amj )
c         ---> print bmj ) = the importance of the history of vi
c
                  WRITE (IOFILE,FMT=
     +              '(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +              ITDEPT - 1, (AM(IT),IT=2,ITDEPT-1)
                  WRITE (IOFILE,FMT=
     +              '(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +              ITDEPT - 1, (BM(IT),IT=2,ITDEPT-1)
c
              ELSE IF (IBROY.EQ.7) THEN
c
c         ---> broyden, fletcher, goldfarb, shanno rank 2 method
c

c         ---> loop to complete generation of v[m] = v(ij,itdept)
c              in two steps
c
                  REWIND IOBROY
                  DO IT = 2,ITDEPT - 1
                      READ (IOBROY,ERR=2106) (UI(IJ,2),IJ=1,NDIM),
     +                  (VI(IJ,2),IJ=1,NDIM)

                      CI = AM(IT)*DDOT(NDIM,UI(1,2),-1,FM1,-1)
                      CALL DAXPY(NDIM,-CI,VI(1,2),1,VI(1,3),1)

                      CI = CI* (ONE-AM(IT)/BM(IT))
                      CI = CI + AM(IT)*DDOT(NDIM,VI(1,2),-1,FM1,-1)
                      CALL DAXPY(NDIM,-CI,UI(1,2),1,VI(1,3),1)
                  END DO

                  BM(ITDEPT) = DDOT(NDIM,VI(1,3),-1,FM1,-1)
c
c         ---> test output of vm-norm
c
                  IF (KTEST.GE.3) THEN
                      WRITE (IOFILE,FMT=
     +'(1x,/,'' vm-norm for bfgs rank 2 method (kmix=8)''
     +,e14.7)') BM(ITDEPT)
                  END IF
                  IF (ABS(BM(ITDEPT)).LT.TOL_70) STOP 'bm(itdept)'
                  BM(ITDEPT) = ONE/BM(ITDEPT)
c
c         ---> print amj )
c         ---> print bmj ) = the importance of the history of vi
c
                  WRITE (IOFILE,FMT=
     +              '(5x,'' amj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +              ITDEPT - 1, (AM(IT),IT=2,ITDEPT-1)
                  WRITE (IOFILE,FMT=
     +              '(5x,'' bmj , ---> j=2,'',i3,/,(9x,1p,7d10.2))')
     +              ITDEPT - 1, (BM(IT),IT=2,ITDEPT-1)
c
              END IF
c
c      ---> update f[m-1] <-- f[m]  ; pot(m-1) <-- pot(m)
c      ---> store on ssd metric g and last iteration fields sm and fm
c
              REWIND IOBROG
              WRITE (IOBROG,ERR=2210) ITDEPT + 1
              WRITE (IOBROG,ERR=2211) (G(IJ),IJ=1,NDIM),
     +          (SM(IJ),IJ=1,NDIM), (FM(IJ),IJ=1,NDIM)
              WRITE (IOBROG,ERR=2212) (WIT(IJ),IJ=2,MITDEPT-1),
     +          (AM(IJ),IJ=2,MITDEPT-1), (BM(IJ),IJ=2,MITDEPT-1)
c
c      ---> test output of iterated potential
c
              IF (KTEST.GE.3) THEN
                  WRITE (IOFILE,FMT=
     +'(1x,/,'' test metric and last iteration fields '',
     +/,14x,'' g '',15x,''sm'',14x,''fm'')')
                  WRITE (IOFILE,FMT=
     +'(1x,''ij='',i4,2x,1pe14.7,'','',2x,1pe14.7,'','',
     +2x,1pe14.7)') (IJ,G(IJ),SM(IJ),FM(IJ),IJ=2,10,2)
              END IF
c
c      ---> store on ssd u(ij,itdept) and v(ij,itdept)
c
              WRITE (IOBROY,ERR=2107) (UI(IJ,3),IJ=1,NDIM),
     +          (VI(IJ,3),IJ=1,NDIM)
c
c      ---> update pot(m+1)
c
c      ---> calculate dpot(m+1) in two steps
c
              CMM = AM(ITDEPT)*DDOT(NDIM,UI(1,3),-1,FM,-1)
              WRITE (IOFILE,FMT='(5x,'' cmm1 = '',1p,d12.4)') CMM

              IF (IBROY.EQ.6) THEN
                  CALL DAXPY(NDIM,ONE-CMM,UI(1,3),1,SM,1)
                  CMM = BM(ITDEPT)*DDOT(NDIM,VI(1,3),-1,FM,-1)
                  WRITE (IOFILE,FMT='(5x,'' cmm2 = '',1p,d12.4)') CMM
                  CALL DAXPY(NDIM,ONE-CMM,VI(1,3),1,SM,1)
              ELSE IF (IBROY.EQ.7) THEN
                  CALL DAXPY(NDIM,ONE-CMM,VI(1,3),1,SM,1)
                  CMM = CMM* (ONE-AM(ITDEPT)/BM(ITDEPT))
                  WRITE (IOFILE,FMT='(5x,'' cmm2a = '',1p,d12.4)') CMM
                  CMM = CMM + AM(ITDEPT)*DDOT(NDIM,VI(1,3),-1,FM,-1)
                  WRITE (IOFILE,FMT='(5x,'' cmm2b = '',1p,d12.4)') CMM
                  CALL DAXPY(NDIM,ONE-CMM,UI(1,3),1,SM,1)
              ELSE
                  WRITE (IOFILE,FMT=
     +              '(5x,''wrong choice of ibroy ='',i4)') IBROY
                  STOP 'ibroy'
              END IF
c
          END IF
c
c   ---> define new potential
c
          DO IJ = 1,NDIM
              XKS_MP1(IJ) = SM(IJ)
              IF (ABS(XKS_MP1(IJ)-XKS_I(IJ)).GT.DELTAMAX) THEN
                  I_FAIL = 1
                  ITDEPT = 0
                  T_DEBYE=T_DEBYE/2
                  WRITE (6,FMT=*)
     +              'QNEWTON: Positions were out of range!'
              END IF
          END DO
c
      END IF

      ITDEPT = ITDEPT + 1
      RETURN

 2101 STOP 'iobroy 1'
 2102 STOP 'iobroy 2'
 2103 STOP 'iobroy 3'
 2104 STOP 'iobroy 4'
 2105 STOP 'iobroy 5'
 2106 STOP 'iobroy 6'
 2107 STOP 'iobroy 7'
 2201 STOP 'iobrog 1'
 2202 STOP 'iobrog 2'
 2203 STOP 'iobrog 3'
 2204 STOP 'iobrog 4'
 2205 STOP 'iobrog 5'
 2206 STOP 'iobrog 6'
 2207 STOP 'iobrog 7'
 2208 STOP 'iobrog 8'
 2209 STOP 'iobrog 9'
 2210 STOP 'iobrog 10'
 2211 STOP 'iobrog 11'
 2212 STOP 'iobrog 12'

      END
