      SUBROUTINE DECIMAREAD(EZ,TK,NPTP1,NPTP2,NPTP3,NPOL,ISPIN,
     &                      LEFTTINVLL,RIGHTTINVLL,vacflag,IENERGY,
     &                      NLBASIS,NRBASIS,NAEZ,KAOEZ)
      IMPLICIT NONE
c     *******************************************************
c     * This subroutine reads in the t-matrices of the left
c     * and right host for the decimation method.
c     * 
c     * The t-matrices are writen out in kloopz1
c     * 
c     * The host files contain the CMOMS neaded to create 
c     * the interatomic potential in sub VINTERFACE
c     * This is going to be read in there.
c     * The input files are for this reason not rewinded.
c     * In case of vacuum in one side read in zero points
c     * For the energy mesh and set vacuum flag ??
c     *******************************************************
c     * This sub works in different modes. if IENERGY = 0     
c     * THis means that we are not in the energy loop 
c     * if IENERGY <> 0 then we want to read the energy IENERGY 
c     * which is then returned
c     *******************************************************
      include 'inc.fi' 
      INTEGER LMMAXD
      PARAMETER (LMMAXD=(LMAXD+1)*(LMAXD+1)) 
      integer nptp1,nptp2,nptp3,npol,ispin,IENERGY
      integer NLBASIS,NRBASIS,NAEZ,KAOEZ(NAEZD+NEMBD)
      double precision tk
      double complex ez(IEMXD),
     &               LEFTTINVLL(LMMAXD,LMMAXD,NEMBD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD)
c
c ..Local variables
c
      INTEGER IL,IERROR,ios
      INTEGER npt1L,npt2L,npt3L,npolL
      INTEGER npt1R,npt2R,npt3R,npolR
      CHARACTER*40 fileleft,fileright  
      CHARACTER*80 BANER1,UIO
      LOGICAL VACFLAG(2)
      DOUBLE PRECISION alatL,alatR,TEMPR,TEMPL,e2r,e2l
      INTEGER nspinL,naezL,lmmaxL,nspinR,naezR,lmmaxR,insr,insl
      DOUBLE PRECISION bravaisL(3,3),bravaisR(3,3),rbasisL(3,NEMBD),
     &                 rbasisR(3,NEMBD)
      DOUBLE COMPLEX ER,DFR,EL,DFL
      INTEGER iesp,iel,ier,ie1,ih,lm1,lm2,lmmax,i,is,idum,ne,ih1
c
c ... SAVE Statement
      SAVE
c 
c -----------------------------------------------------  
c     
      IF (IENERGY.LT.1) THEN ! READ just the headers

         WRITE(6,*) '>>>>> Reading the host t-matrices '
         IL=1
         CALL IoInput('DECIFILES ',UIO,IL,7,IERROR)
         READ (UNIT=UIO,FMT='(A40)')  fileleft
         CALL IoInput('DECIFILES ',UIO,IL+1,7,IERROR)
         READ (UNIT=UIO,FMT='(A40)')  fileright
         vacflag(1) = .false.
         vacflag(2) = .false.  
                  
         IF (fileleft(1:7).ne.'vacuum') THEN
            
            OPEN(37,file=fileleft,status='old',IOSTAT=ios)
            IF (IOS.GT.0) THEN
               WRITE(6,*) ' File : ',fileleft
               WRITE(6,*) ' not found '
               write(6,*) ' use   vacuum  as filename in case of vacuum'
               STOP 'READ DECIMATION'
            END IF     
            write(6,*) '>>>>> Reading file : ',fileleft
            READ  (37,2210) BANER1
            READ  (37,2210) BANER1
            READ  (37,2210) BANER1
            READ  (37,2220) alatL,nspinL,naezL,lmmaxL,insL
            WRITE(6,1110) alatL,nspinL,naezL,lmmaxL,insL
            IF (NAEZL.NE.NLBASIS) THEN 
               WRITE(6,*) 'Left Host atoms NAEZL should be',NLBASIS
               STOP 'DECIMAREAD ERROR '
            END IF
            READ  (37,2210) BANER1
            READ  (37,1130) bravaisL
            WRITE (6,1120) BRAVAISL
            READ  (37,*)    BANER1
            WRITE (6,*)  BANER1
            DO IH=1,naezL             
               READ(37,1130) (rbasisL(I,IH),I=1,3)
               WRITE(6,1130) (rbasisL(I,IH),I=1,3)
            END DO
            READ(37,2240) e2L,TempL
            WRITE(6,1140) e2L,TempL 
            READ(37,2250) npt1L,npt2L,npt3L,npolL
            WRITE(6,1150) npt1L,npt2L,npt3L,npolL
c     
            NE = NPT1L+NPT2L+NPT3L+NPOLL 
            IF ((NPT1L.NE.NPTP1).OR.(NPT2L.NE.NPTP2).OR.
     &           (NPT3L.NE.NPTP3).OR.(NPOLL.NE.NPOL).OR.
     &           (abs(templ-tk).GT.1.d-6)) THEN
               WRITE(6,*) ' ERROR During reading Left Host ',FILELEFT
               WRITE(6,*) ' for the Decimation, The energy mesh '
               WRITE(6,*) ' is not Consistent'
               WRITE(6,*) ' *** PROGRAM IS STOPING *** '
               STOP  
            END IF

         ELSE           
            WRITE(6,*) ' Vacuum is used on the left side '
            VACFLAG(1) = .true.
         END IF

c     ------   Now Read Right Host Header    -------
c     
c     
         IF (fileright(1:7).ne.'vacuum') THEN
            
            OPEN(38,file=fileright,status='old',IOSTAT=ios)
            IF (IOS.GT.0) THEN
               WRITE(6,*) ' File : ',fileright
               WRITE(6,*) ' not found '
               write(6,*) ' use   vacuum  as filename in case of vacuum'
               STOP 'READ DECIMATION'
            END IF
            write(6,*) '>>>>> Reading file : ',fileright         
            READ  (38,2210) BANER1
            READ  (38,2210) BANER1
            READ  (38,2210) BANER1
            READ  (38,2220) alatR,nspinR,naezR,lmmaxR,insR
            WRITE(6,1110) alatR,nspinR,naezR,lmmaxR,insR
            IF (NAEZR.NE.NRBASIS) THEN
               WRITE(6,*) 'Right Host atoms NAEZR should be',NRBASIS
               STOP 'DECIMAREAD ERROR '
            END IF
            READ (38,2210) BANER1 
            READ  (38,1130) bravaisR             
            WRITE (6,1120) BRAVAISR
            READ  (38,*)    BANER1
            WRITE (6,*)  BANER1
            DO IH=1,naezR             
               READ(38,1130) (rbasisR(I,IH),I=1,3)
               WRITE(6,1130) (rbasisR(I,IH),I=1,3)
            END DO
            READ(38,2240) e2R,TempR
            WRITE(6,1140) e2R,TempR 
            READ(38,2250) npt1R,npt2R,npt3R,npolR
            WRITE(6,1150) npt1R,npt2R,npt3R,npolR
            NE = NPT1R+NPT2R+NPT3R+NPOLR 
            IF ((NPT1R.NE.NPTP1).OR.(NPT2R.NE.NPTP2).OR.
     &           (NPT3R.NE.NPTP3).OR.(NPOLR.NE.NPOL).OR.
     &           (abs(tempR-tk).GT.1.d-6)) THEN
               WRITE(6,*) ' ERROR During reading Right Host ',FILERIGHT
               WRITE(6,*) ' for the Decimation, The energy mesh '
               WRITE(6,*) ' is not Consistent'
               WRITE(6,*) ' *** PROGRAM IS STOPING *** '
               STOP  
            END IF
         ELSE 
            WRITE(6,*) ' Vacuum is used on the right side '
            VACFLAG(2) = .true.
         END IF  
c
c Headers read in this now returns to the END IF at the end  
c

      ELSE   ! NOW READ T-MATRICES IENERGY > 0  


         IF (.NOT.VACFLAG(1)) THEN ! not vacuum
            IF ((ISPIN.NE.2).OR.(NSPINL.NE.1)) THEN 
c     
C     Now Reading the t-matrices
c     
               READ(37,2260) ieL,el,dfl
               IF (ABS(eL-EZ(IENERGY)).GT.1.D-6) THEN
                  WRITE(6,*)' ERROR During reading Left Host ',FILELEFT
                  WRITE(6,*)' for the Decimation, The energy mesh '
                  WRITE(6,*)' is not Consistent ',IENERGY,EL,EZ(IENERGY)
                  WRITE(6,*)' *** PROGRAM IS STOPING *** '
                  STOP     
               END IF
               write(6,*) 'ISPIN=',ISPIN,INSL,LMMAXL
               DO IH=1,naezL
                  READ(37,2270) IDUM
                  IF (IDUM.NE.IH) STOP ' DECIMAREAD ERROR 1'
                  IH1 = KAOEZ(NAEZ+IH) ! map the t-matrix corectly
                  !write(6,*) 'DEDEDEDEDEDDE',ih1
                  IF (INSL.GT.0) THEN ! read full potential
                     READ(37,2280) ((LEFTTINVLL(LM1,LM2,IH1),
     &                    LM2=1,LMMAXL),LM1=1,LMMAXL)
                  ELSE          ! read spherical
                     DO LM1=1,LMMAXL
                        DO LM2=1,LMMAXL
                           LEFTTINVLL(LM1,LM2,IH1) = 0.D0
                        END DO
                     END DO
                     READ(37,2280) (LEFTTINVLL(LM1,LM1,IH1),
     &                    LM1=1,LMMAXL)
                  END IF           
               END DO           ! Loop in atoms          
c     
c     T- matrices read ok!                 
            END IF              !   ((ISPIN.NE.2).OR.(NSPINL.NE.1)) 
c     
c     In case of spin polarized calculation with paramagnetic host
c     do NOT read the file just keep the previous values 
c     with the SAVE statement.
c     
      
         END IF                 ! if it is vacuum do nothing
c     ==========================================
c     Right Host
c     ==========================================
         
         IF (.NOT.VACFLAG(2)) THEN       ! not vacuum  
            IF ((ISPIN.NE.2).OR.(NSPINR.NE.1)) THEN 
               READ(38,2260) ieR,eR,dfR
               IF (ABS(eR-EZ(IENERGY)).GT.1.D-6) THEN
                  WRITE(6,*) ' ERROR During reading Right Host ',
     &                 FILERIGHT
                  WRITE(6,*) ' for the Decimation, The energy mesh '
                 WRITE(6,*) ' is not Consistent ',IENERGY,ER,EZ(IENERGY)
                  WRITE(6,*) ' *** PROGRAM IS STOPING *** '
                  STOP     
               END IF
               write(6,*) 'ISPIN=',ISPIN,INSL,LMMAXL
               DO IH=1,naezR
                  READ(38,2270) IDUM
                  IF (IDUM.NE.IH) STOP ' DECIMAREAD ERROR 1'
                  IH1 = KAOEZ(NAEZ+NLBASIS+IH)
                  !write(6,*) 'RIGHT lkadlhash',ih1
                  IF (INSR.GT.0) THEN ! read full potential
                     READ(38,2280) ((RIGHTTINVLL(LM1,LM2,IH1),
     &                    LM2=1,LMMAXR),LM1=1,LMMAXR)
                  ELSE          ! read spherical
                     DO LM1=1,LMMAXR
                        DO LM2=1,LMMAXR
                           RIGHTTINVLL(LM1,LM2,IH1) = 0.D0
                        END DO
                     END DO
                     READ(38,2280) (RIGHTTINVLL(LM1,LM1,IH1),
     &                    LM1=1,LMMAXR)
                  END IF
c     
c     In case of spin polarized calculation with paramagnetic host
c     do NOT read the file just keep the previous values 
c     with the SAVE statement.
c     
                  
                  
               END DO           ! Loop in the atoms
            END IF              ! ((ISPIN.NE.2).OR.(NSPINR.NE.1)) 
         END IF                 ! VACFLAG
         
      END IF                    ! IENERGY = 0
      
      
      WRITE(6,*) ' >>>>>>>> Host t-matrices read in '
      RETURN
c     Read in Format
 1110 FORMAT('ALAT=',F9.6,' NSPIN=',I2,'  NAEZ=',I3,' LMMAX=',I3,
     &     ' INS=',I1)
 1120 FORMAT('BRAVAIS '/3F8.4/3F8.4/3F8.4)
 1125 FORMAT('RBASIS')
 1130 FORMAT(3F8.4)
 1140 FORMAT('EF=',F10.6,' TEMP=',F10.4,' Kelvin')
 1150 FORMAT('N1=',I3,' N2=',I3,' N3=',I3,' NPOL=',I3)
c      
 2210 FORMAT(A80)
 2220 FORMAT(5X,F9.6,7X,I2,7X,I3,7X,I3,5X,I1)
 2240 FORMAT(3X,F10.6,6X,F10.4)
 2250 FORMAT(3X,I3,4X,I3,4X,I3,6X,I3)
 2260 FORMAT(7X,I5,4D16.8)
 2270 FORMAT(5X,I3)
 2280 FORMAT(4D22.14)
      
      END          
