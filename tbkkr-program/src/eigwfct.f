c ************************************************************************
      SUBROUTINE EIGWFCT(NSPIN,NAEZ,NINEQ,LMAX,LSYM,RBASIS,KAOEZ,
     +                   QBOUND)
      implicit none
c ************************************************************************
c
c --->  output of wavefunction coefficients for a layered structure
c
c ------------------------------------------------------------------------
      include 'inc.fi'
C
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
C
      INTEGER LMAX,NSPIN,NINEQ,NAEZ,
     +     KAOEZ(*)
      LOGICAL LSYM
C
      DOUBLE PRECISION  
     +     AKSCALE,BKSCALE,CKSCALE,
     +     CSQR,
     +     DOS,
     +     E1,EFERMI,EDIFF,ED0,ED1,ED2,EDMAX,
     +     FAC,
     +     QBOUND,
     +     RQ,
     +     W1,W2
      DOUBLE PRECISION  
     +     AR(-7:7),
     +     RBASIS(3,*),             ! position of atoms in the unit cell
                                    ! in units of bravais vectors
     +     VF(3),Q(3),QQ(3),
     +     WFCT(0:LMAXD+1,NAEZD)
      DOUBLE COMPLEX E,IKR,
     +     CL1(-7:7),
     +     EXPIKR(NAEZD)
      INTEGER
     +     ILAY(NAEZD),             ! layer of atom
     +     NLAY(NAEZD)              ! number of atoms in layer
      INTEGER
     +     I,I1,I2,IBAND,II,IILM,IL,IN,INMAX,IPZ,ISPIN,
     +     LAYMAX,L,LM,LMMAX,L2,M,
     +     NP,NPOINTS,NUEL
      CHARACTER*1   HEAD(80)
C
      DOUBLE PRECISION PI,TPI,ZERO,ONE
      DOUBLE COMPLEX CI,CC,CZERO
      PARAMETER (
     +     PI    = 3.14159265358979312D0,
     +     TPI   = 2.0d0 * PI,
     +     ZERO  = 0.0D0,
     +     ONE   = 1.0D0,
     +     CZERO = ( 0.0D0,0.0D0),
     +     CI    = (0.D0,1.D0))
C
      LOGICAL OPT,TEST
      EXTERNAL EIGEN0,MINEIG,OPT,RINIT,TEST
      INTRINSIC DMIN1,DSIGN,DABS,DATAN
C
C     .. COMMONS
C
c ------------------------------------------------------------------------
      DOUBLE COMPLEX CL(NAEZD*LMMAXD)
      COMMON /CL/ CL
c ------------------------------------------------------------------------
      DOUBLE COMPLEX
     +     PZSQ(LMMAXD,NATYPD)
      COMMON /PZSQ/ PZSQ
c ------------------------------------------------------------------------
      INTEGER     IUP,IDO,NUM
      DOUBLE PRECISION WUP,WDO
      COMMON /EIGENV/ WUP,WDO,NUM,IUP,IDO
      
       integer N !for fo function
       integer NUE
       integer J
c ------------------------------------------------------------------------
c     bounds for q refinement for eigenvector determination
      LMMAX =(LMAX+1)*(LMAX+1)
      EDMAX = .05D0                 ! maximum energy shift for eigenvalue 
                                    ! determination
      QBOUND = MIN(QBOUND,1.D-6)    ! maximum value of smallest KKR matrix 
                                    ! eigenvalue
      INMAX = 6                     ! maximum number of refinement steps
c ------------------------------------------------------------------------
c
c --->  input of layer structure
c
      IF (TEST('lay ind ')) THEN
        OPEN(30,FILE='layer.index',STATUS='OLD')
        READ(30,*) (ILAY(I),I=1,NAEZ)
        CLOSE(30)
        LAYMAX = 1
        DO 161 I = 1,NAEZD
          NLAY(I) = 0
 161    END DO
        DO 162 I=1,NAEZ
          LAYMAX=MAX(LAYMAX,ILAY(I))
          NLAY(ILAY(I)) = NLAY(ILAY(I)) + 1
 162    END DO
      ELSE
        DO 163 I=1,NAEZ
          ILAY(I) = I
          NLAY(I) = 1
 163    END DO
        LAYMAX = NAEZ
      END IF                        ! (TEST('lay ind '))
      write(6,FMT='('' ILAY(I) :'',10I4,/,(10x,10I4))') 
     +     (ILAY(I),I=1,NAEZ)
c ------------------------------------------------------------------------
        OPEN(30,FILE='numdat',STATUS='OLD')
        REWIND 30

        OPEN(31,FILE='numdat.wfct',STATUS='UNKNOWN')
        REWIND 31

        IF (TEST('C(L)    ')) THEN
          OPEN(32,FILE='numdat.cl',STATUS='UNKNOWN')
          REWIND 32
        END IF

        DO 90 ISPIN=1,NSPIN

          READ(30,FMT=9020) (HEAD(I),I=1,79)
          READ(30,*) EFERMI,AKSCALE,BKSCALE,CKSCALE,FAC
          READ(30,*) I1,I2
c
c --->    FAC is only used for test purposes
c
c          IF (.not.OPT('use fac ')) fac = 1.0d0

          WRITE(31,FMT=9020) (HEAD(I),I=1,79)
          WRITE(31,FMT=9030) 
     +         EFERMI,AKSCALE,BKSCALE,CKSCALE,FAC
          WRITE(31,FMT=9050) 
     +         I1,I2,LMMAX,NAEZ

          IF (TEST('C(L)    ')) THEN
            WRITE(32,FMT=9020) (HEAD(I),I=1,79)
            WRITE(32,FMT=9030) 
     +           EFERMI,AKSCALE,BKSCALE,CKSCALE,FAC
            WRITE(32,FMT=9050) 
     +           I1,I2,LMMAX,NAEZ
          END IF

          WRITE(6,FMT=9020) 
     +         (HEAD(I),I=1,79)
          WRITE(6,FMT=9030) 
     +         EFERMI,AKSCALE,BKSCALE,CKSCALE,FAC
          WRITE(6,FMT=9050) 
     +         I1,I2,LMMAX,NAEZ

          IF (DABS(EFERMI-E1).gt.1.0d-5) THEN
            write(6,*) 
     +         'WARNING : EFERMI in numdat is not equal E1 in in5.'
            IF (.NOT.OPT('rigid-ef')) THEN
              E1 = EFERMI
              write(6,FMT='('' E1 is set to '',f12.6 )') E1 
            END IF
          END IF

          AKSCALE = TPI/AKSCALE ! /FAC

          write (6,*) 'akscale :',akscale

          DO 230 IBAND =I1,I2

            WRITE(6,*) 'IBAND :',IBAND

            READ(30,*) I,NPOINTS
            NP = NPOINTS
            IF (TEST('UP      ') .AND. ISPIN.EQ.1) NP = 0
            IF (TEST('DOWN    ') .AND. ISPIN.EQ.2) NP = 0

            WRITE(31,FMT='(2I5)') I,NP
            IF (TEST('C(L)    ')) WRITE(32,FMT='(2I5)') I,NP

            IF (NPOINTS.GT.0) THEN
            DO 240 N =1,NPOINTS
              READ(30,FMT='(F13.10,6F12.7)') 
     +             DOS,(VF(I),I=1,3),(QQ(I),I=1,3)

              IF (TEST('UP      ') .AND. ISPIN.EQ.1) GOTO 240
              IF (TEST('DOWN    ') .AND. ISPIN.EQ.2) GOTO 240

              RQ = 1.D0
              IF (TEST('NFE     ')) THEN
                RQ = SQRT(E1)/SQRT(QQ(1)**2+QQ(2)**2+QQ(3)**2)
              END IF

              Q(1)=QQ(1)/AKSCALE*RQ
              Q(2)=QQ(2)/AKSCALE*RQ
              Q(3)=QQ(3)/AKSCALE*RQ

              WRITE (6,FMT='(I5,5X,3F10.6,5x,3F10.6)') 
     +             N,(QQ(I),I=1,3),(Q(I),I=1,3)

              E = CMPLX(E1,ZERO)

              CALL EIGEN0(ISPIN,Q,E,LSYM,0)

              CALL MINEIG(W1)
              ED1 = DMIN1(DABS(W1/2.0D0),0.01d0)*DSIGN(1.0d0,W1)
              IF (.NOT.OPT('full inv')) ED1 = ED1/1.D2      
              write(6,FMT=9110) W1,ED1
              DO I=1,3
                Q(I) = (QQ(I)-VF(I)*ED1)/AKSCALE
              END DO

              E = CMPLX(E1,ZERO)
              CALL EIGEN0(ISPIN,Q,E,LSYM,0)
              CALL MINEIG(W2)

              ED2 = ED1
              ED1 = 0.d0
              
              IN = 0
              DO WHILE (ABS(W2).GT.QBOUND    .AND. 
     +                  IN.LE.INMAX          .AND.
     +                  .NOT.TEST('noqshift')   )
c
c --->          refinement of q-vector due to value of smallest 
c               KKR-matrix-eigenvalue W2
c
                ED0 = ED2
                ED2 = (ED2-ED1)*W1/(W1-W2) + ED1
                IF (ABS(ED2).GT.EDMAX) THEN
                  write(6,FMT='(1P,D12.4,''  %%%%%%%%%%%'')') ED2
                  DOS = .0d0
                  IN = INMAX
                END IF
                ED2 = DMIN1(DABS(ED2),EDMAX)*DSIGN(1.d0,ED2)
                ED1 = ED0
                W1 = W2
                write(6,FMT=9100) W2,ED2
                DO I=1,3
                  Q(I) = (QQ(I)-VF(I)*ED2)/AKSCALE
                END DO
                IF (TEST('qs_out  ')) write(6,FMT=9130) IN,(Q(I),I=1,3)
c
                E = CMPLX(E1,ZERO)
                CALL EIGEN0(ISPIN,Q,E,LSYM,0)
                CALL MINEIG(W2)
c
                IN = IN + 1
              END DO

              ED0 = ED2
              ED2 = (ED2-ED1)*W1/(W1-W2) + ED1
              ED2 = DMIN1(DABS(ED2),EDMAX)*DSIGN(1.d0,ED2)
              IF (TEST('noqshift')) ED2 = 0.d0
              ED1 = ED0
              W1 = W2
              write(6,FMT=9100) W2,ED2
              DO I=1,3
                Q(I) = (QQ(I)-VF(I)*ED2)/AKSCALE
              END DO

              E = CMPLX(E1,ZERO)
              CALL EIGEN0(ISPIN,Q,E,LSYM,2)
              write(6,FMT=9120) WDO

              IF (NUM.NE.IBAND) THEN
                write(6,FMT='(2I6,1P,D12.2,''   <<<<<>>>>>'')') 
     +               IBAND,NUM,WUP
                DOS = 0.d0
              END IF
              IF (ABS(WDO).GT.QBOUND) THEN
                write(6,FMT='(1P,D12.2,''  !!!!!!!!!!!'')') WDO
                DOS = 0.d0
              END IF
c
c ------------------------------------------------------------------------
c
c --->        cancel the BLOCH factor EXPIKR(i)=exp(i*q*rbasis(i))
c
              DO 10 II=1,NAEZ
                IILM=(II-1)*LMMAXD
                IKR=-CI*TPI*( RBASIS(1,II)*Q(1)
     +                     + RBASIS(2,II)*Q(2)
     +                     + RBASIS(3,II)*Q(3) )
                EXPIKR(II)=ZEXP(IKR)
                DO 20 LM=1,LMMAXD
                  CL(IILM+LM)=CL(IILM+LM)*EXPIKR(II)/EXPIKR(1)
 20             END DO
 10           END DO
c ------------------------------------------------------------------------
              DO 430 II = 1,NINEQ
                IF (TEST('csq(L)  ') .OR.
     +              TEST('phase(i)') .OR.
     +              TEST('c(L)    ') )
     +               write(6,FMT='('' atom '',I4)') II

                IF (TEST('csq(L)  ')) THEN
                  DO 431 L=0,LMAX
                    L2=(II-1)*LMMAXD+L*(L+1)+1
                    WRITE(6,FMT=9090)
     +                   L,(ABS(CL(L2+M))**2,M=-L,L)
 431              END DO
                END IF              ! (TEST('csq(L)  '))

                IF (TEST('phase(i)')) 
     +               write(6,FMT=9080) 
     +               ii,dreal(cl((ii-1)*lmmaxd+1)/cl(1))

                IF (TEST('c(L)    ')) THEN
                  DO 432 L = 0,LMAX
                    LM = (II-1)*LMMAXD+L*(L+1)+1
                    DO 434 M = -L,L
                      CL1(M)=CL(LM+M)
                      AR(M)=DREAL(-CI*ZLOG(CL1(M)))/TPI
                      IF (ABS(CL1(M)).LT.1.D-9) AR(M)  = 0.D0
                      IF (ABS(CL1(M)).LT.1.D-9) CL1(M) = CZERO
                      IF (AR(M).LT.-1.d-8) AR(M) = AR(M) + 1.D0 
 434                END DO          ! M = -L,L
                    WRITE(6,FMT=9070) L,(CL1(M),AR(M),M=-L,L)
 432              END DO            ! L = 0,LMAX
                END IF              ! (TEST('c(L)    '))

 430          END DO                ! II = 1,NINEQ
c ------------------------------------------------------------------------
c
c --->        init of wave function coeffizients
c
             CALL RINIT(NAEZD*(LMAXD+2),WFCT)
c
c --->        sum up the l-contributions in the layers defined by
c             NLAY()
c
              DO 180 I=1,NAEZ
                IL= ILAY(I)
                II= (I-1)*LMMAXD+1
                DO 170 L= 0,LMAX
                  NUEL = II + L*(L+1)
                  DO 190 M = -L,L
                    NUE             = NUEL + M
                    CSQR = ABS(CL(NUE)*CL(NUE))
                    WFCT(L,IL)      = WFCT(L,IL) 
     +                   + CSQR/DBLE(NLAY(IL))
                    WFCT(LMAX+1,IL) = WFCT(LMAX+1,IL) 
     +                   + CSQR/DBLE(NLAY(IL))
 190              END DO            ! M = -L,L
 170            END DO              ! L = 0,LMAX
 180          END DO                ! I = 1,NAEZ

              IF (TEST('C(l)    ')) THEN
c
c --->          output of l-dep. (wfct coeff.)**2
c
                IF (TEST('larg_acc')) THEN
                  WRITE(31,FMT=9001) 
     +                 DOS,(VF(I),I=1,3),(QQ(I),I=1,3),
     +                 ((WFCT(I,J),I=0,LMAX+1),J=1,LAYMAX)
                ELSE
                  WRITE(31,FMT=9000) 
     +                 DOS,(VF(I),I=1,3),(QQ(I),I=1,3),
     +                 ((WFCT(I,J),I=0,LMAX+1),J=1,LAYMAX)
                END IF              ! TEST('larg_acc')
              ELSE                  ! TEST('C(l)    ')
c
c --->          only output of sum of (wfct coeff.)**2
c
                IF (TEST('larg_acc')) THEN
                  WRITE(31,FMT=9001) 
     +                 DOS,(VF(I),I=1,3),(QQ(I),I=1,3),
     +                 (WFCT(LMAX+1,J),J=1,LAYMAX)
                ELSE
                  WRITE(31,FMT=9000) 
     +                 DOS,(VF(I),I=1,3),(QQ(I),I=1,3),
     +                 (WFCT(LMAX+1,J),J=1,LAYMAX)
                END IF              ! TEST('larg_acc')
              END IF                ! TEST('C(l)    ')

              IF (TEST('C(L)    ')) THEN
c
c --->            output of (l,m)-dep. wfct coeff.
c
                WRITE(32,FMT=9010) 
     +               ((CL((I-1)*LMMAXD+LM),LM=1,LMMAX),I=1,NAEZ)
              END IF                ! (TEST('C(L)    '))

 240        END DO                  ! N = 1, NPOINTS
            END IF                  ! (NPOINTS.GT.0)
 230      END DO                    ! IBAND =I1,I2
c ------------------------------------------------------------------------
c
c --->    write norm of radial wave functions PZSQ to file
C         (see COMMON PZSQ)
c
          IF (TEST('C(L)    ')) 
     +         WRITE(32,FMT=9060) 
     +         ((PZSQ(L,KAOEZ(I)),L=1,LMMAX),I=1,NAEZ)
          IF (TEST('C(L)    ')) 
     +         WRITE(6,*) (KAOEZ(I),I=1,NAEZ)
c ------------------------------------------------------------------------
 90     END DO                       ! ISPIN=1,NSPIN

        CLOSE(30)
        CLOSE(31)
        IF (TEST('C(L)    ')) CLOSE(32)

        RETURN

 9000   FORMAT(F13.10,6F12.7,100F8.5)
 9001   FORMAT(F13.10,6F12.7,100F12.9)
 9010   FORMAT(1p,5d16.8)
 9020   FORMAT(80A1)
 9030   FORMAT(f10.6,3f10.5,f12.8)
 9040   FORMAT(2I5)
 9050   FORMAT(4I5)
 9060   FORMAT(1P,2D16.8)
 9070   FORMAT(1P,I10,2D14.4,0P,F10.5,/,(10X,1P,2D14.4,0P,F10.5))
 9080   FORMAT(I6,f12.4)
 9090   FORMAT(I10,7F8.5,/,(10X,7F8.5))
 9100   FORMAT(6x,'W2,ED2 :',d12.2,d12.5)
 9110   FORMAT(6x,'W1,ED1 :',d12.2,d12.5)
 9120   FORMAT(6x,'W2     :',d12.2)
 9130   FORMAT(6x,i6,3f14.8)
        END
