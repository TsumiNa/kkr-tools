      SUBROUTINE GENERALPOT(IFILE,NATPS,NATYP,NSPIN,Z,ALAT,RMT,RMTNEW,
     +                 RWS,ITITLE,R,DRDI,VM2Z,IRWS,A,B,TXC,KXC,INS,IRNS,
     +                 LPOT,VINS,QBOUND,IRC,KSHAPE,EFERMI,VBC,ECORE,
     +                 LCORE,NCORE)
c **************************************************
c * The subroutine writes out the potential cards
c * in a standard r-mesh that can be read in and 
c * interpolated to a different r-mesh from subroutine
c * start No shape function information is neaded
c * and all nessecery data are stored in the potential 
c * card.
c *                                      ver. 18.5.2000
c ***************************************************
      implicit none
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ALAT,QBOUND
      INTEGER IFILE,INS,KSHAPE,KXC,LPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*),B(*),DRDI(IRMD,*),ECORE(20,*),EFERMI,
     +                 R(IRMD,*),RMT(*),RMTNEW(*),RWS(*),VBC(2),
     +                 VINS(IRMIND:IRMD,LMPOTD,*),VM2Z(IRMD,*),Z(*)
      INTEGER IRC(*),IRNS(*),IRWS(*),ITITLE(20,*),LCORE(20,*),NCORE(*)
      CHARACTER*24 TXC(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION A1,B1,RMAX,RMT1,RMTNW1,RV,SIGN,SUM,Z1,parsum,
     &                 PARSUMDERIV,R0,RINTER,DR,MAXA,VCON
      INTEGER I,ICORE,IH,INEW,IP,IR,IRMIN,IRNS1,IS,ISAVE,J,LM,LMNR,
     &        LMPOT,NCORE1,NR,NZ1,NR_U,IRMIN_U,IRNS_U,
     &        IMT1,LM1,IRNSTOT
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION DRADI(IRMD),ECORE1(20),RA(IRMD),VM2ZA(IRMD),
     &                 RR_U(IRMD),DRDI_U(IRMD)
      double precision vm2zb(IRMD),vm2z_U(irmd),
     &            vins_u(IRMIND:IRMD,LMPOTD),vinsa(IRMIND:IRMD,LMPOTD),
     &            vinsb(IRMIND:IRMD,LMPOTD)
      INTEGER LCORE1(20)
       CHARACTER*4 ELEMNAME(0:113)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC SQRT
c        1      2      3      4      5      6      7      8      9    
      DATA ELEMNAME/'VAC',
     &  'H   ','He  ','Li  ','Be  ','B   ','C   ','N   ','O   ','F   ',
     &  'Ne  ',
     &  'Na  ','Mg  ','Al  ','Si  ','P   ','S   ','Cl  ','Ar  ','K   ',
     &  'Ca  ',
     &  'Sc  ','Ti  ','V   ','Cr  ','Mn  ','Fe  ','Co  ','Ni  ','Cu  ',
     &  'Zn  ',
     &  'Ga  ','Ge  ','As  ','Se  ','Br  ','Kr  ','Rb  ','Sr  ','Y   ',
     &  'Zr  ',
     &  'Nb  ','Mo  ','Tc  ','Ru  ','Rh  ','Pd  ','Ag  ','Cd  ','In  ',
     &  'Sn  ',
     &  'Sb  ','Te  ','I   ','Xe  ','Cs  ','Ba  ','La  ','Ce  ','Pr  ',
     &  'Nd  ',
     &  'Pm  ','Sm  ','Eu  ','Gd  ','Tb  ','Dy  ','Ho  ','Er  ','Tm  ',
     &  'Yb  ',
     &  'Lu  ','Hf  ','Ta  ','W   ','Re  ','Os  ','Ir  ','Pt  ','Au  ',
     &  'Hg  ',
     &  'Tl  ','Pb  ','Bi  ','Po  ','At  ','Rn  ','Fr  ','Ra  ','Ac  ',
     &  'Th  ',
     &  'Pa  ','U   ','Np  ','Pu  ','Am  ','Cm  ','Bk  ','Cf  ','Es  ',
     &  'Fm  ',
     &  'Md  ','No  ','Lr  ','Rf  ','Db  ','Sg  ','Bh  ','Hs  ','Mt  ',
     &  'Uun ','Uuu ','Uub ','NoE '/

C     ..
      ISAVE = 1
      LMPOT = (LPOT+1)* (LPOT+1)
   
      DO 60 IH = 1,NATYP
        DO 50 IS = 1,NSPIN
          VINSA = 0.d0
          VINSB = 0.D0       
          IP = NSPIN* (IH-1) + IS
          RMT1 = RMT(IH)
          RMTNW1 = RMTNEW(IH)
          Z1 = Z(IH)
          RMAX = RWS(IH)
          IF (KSHAPE.EQ.0) THEN
            NR = IRWS(IH)
            IRNS1 = 0
          ELSE
            NR = IRC(IH)
            IRNS1 = IRNS(IH)
          END IF

          IRMIN = NR - IRNS1
          A1 = A(IH)
          B1 = B(IH)
          NCORE1 = NCORE(IP)
c
          DO 10 J = 1,NR
            RA(J) = R(J,IH)
            DRADI(J) = DRDI(J,IH)
            VM2ZA(J) = VM2Z(J,IP) 
   10     CONTINUE
          do lm1=1,lmpot
            do j=irmind,irmd
            VINSA(j,LM1) = vins(j,lm1,ip)
            end do
          end do 
c
          IF (NCORE1.GE.1) THEN
c
            DO 20 J = 1,NCORE1
              LCORE1(J) = LCORE(J,IP)
              ECORE1(J) = ECORE(J,IP)
   20       CONTINUE
          END IF
c
c  Generate uniform mesh RUNI
c     
          NR_U = NR
          IRNS_U = IRNS1
          IRMIN_U = NR_U
          if (ins.gt.0) IRMIN_U = NR_U - IRNS_U
  
          IF (INS.EQ.0) THEN
             DO I=1,NR_U
                RR_U(I) = RA(I)
                DRDI_U(I) = DRADI(I)
             END DO 
          ELSE
             IMT1 = ANINT(LOG(RMTNW1/B1+1.0D0)/A1) + 1
             DO I=1,IMT1
                RR_U(I) = RA(I)
                DRDI_U(I) = DRADI(I)
             END DO
             RINTER =  RMAX - RMTNW1
             DR = RINTER/FLOAT(NR-IMT1) 
             Do I=1,NR-IMT1
                DRDI_U(IMT1+I) = DR
                RR_U(IMT1+I) = RR_U(IMT1) + DR*FLOAT(I)
             END DO
             CALL DOUBLERAUS1(NR,IRMIN,LMPOT,RA,DRADI,VM2ZA,VINSA)
c     
c     After this sub the arrays are rearanged and nr is not 
c     the same anymore in the case of FP. If ins.eq.0 there is
c     no nead for doubleraus1. IRMIN should remain the same
c  
          END IF
c ----------------------------------------------------------------
c Now the new mesh is generated
c     
c test     
         ! write(6,*) nr_u,imt1,irns_u
         ! do i=1,nr
         !   write(60,*) i,ra(i),VM2ZA(i)  
         ! end do
          
c test
c
          maxa = 1.d35
          CALL SPLINE(IRMD,RA,VM2ZA,NR,maxa,maxa,VM2ZB)
          IF (INS.GT.0) THEN
             DO LM1=1,LMPOT 
                IRNSTOT = NR - IRMIN ! nr has changed irmin is the same
                !write(6,*) ' Testing ',nr,irmin,irnstot,irmind 
                CALL SPLINE(IRMD-IRMIND,RA(IRMIND),
     &                      VINSA(IRMIND,LM1),IRNSTOT,maxa,maxa,
     &                      VINSB(IRMIND,LM1))
             ENDDO              ! LM1
          END IF
c     
c OK with spline 
c
          DO IR = 1,NR_U
             R0 = RR_U(IR)
             CALL SPLINT(RA,VM2ZA,VM2ZB,NR,R0,PARSUM,PARSUMDERIV)
             VM2Z_U(IR) = PARSUM
          END DO
          IF (INS.GT.0) THEN
             !IRNSTOT = NR_U - IRMIN_U 
             DO LM1=1,LMPOT
                DO IR = IRMIN_U,NR_U
                   R0 = RR_U(IR)
                   CALL SPLINT(RA(IRMIND),VINSA(IRMIND,LM1),
     &                         VINSB(IRMIND,LM1),IRNSTOT,R0,
     &                         PARSUM,PARSUMDERIV)
                   VINS_U(IR,LM1) = PARSUM 
                END DO
             end do
          END IF 
          !write(6,*) ' All interpolation ok now write'             
c     --------------------------------------------------------------
          WRITE (IFILE,FMT=8000)
          NZ1 = Z1
          IF (NSPIN.EQ.1) THEN
          WRITE (IFILE,FMT=8010) ELEMNAME(NZ1),Z1
          ELSEIF (IS.EQ.1) THEN 
          WRITE (IFILE,FMT=8012) ELEMNAME(NZ1),Z1
          ELSEIF (IS.EQ.2) THEN 
          WRITE (IFILE,FMT=8011) ELEMNAME(NZ1),Z1
          END IF
          WRITE (IFILE,FMT=8020)
c          write (ifile,*) ALAT,RMAX,RMTNW1,RMT1
          WRITE (IFILE,FMT=8030) ALAT,RMAX,RMTNW1,RMT1
          WRITE (IFILE,FMT=8040) NR_U,IMT1,IRNS1
          WRITE (IFILE,FMT=8050) A1,B1
	    WRITE (IFILE,FMT=8060) EFERMI,VBC(IS)
          WRITE (IFILE,FMT=8070) NCORE1,LMPOT
          IF (NCORE1.GE.1) WRITE (IFILE,FMT=9040) (LCORE1(ICORE),
     +         ECORE1(ICORE),ICORE=1,NCORE1)
          
          IF (INS.EQ.0 .OR. (IH.LT.NATPS.AND.INS.LE.2)) THEN
c     
c---  >       store only the spherically averaged potential 
c     (in mt or as - case) 
c     this is done always for the host
c     
             WRITE (IFILE,FMT=9051) (VM2Z_U(IR),IR=1,NR_U)
          ELSE
c     
c---  >     store the full potential , but the non spherical contribution
c     only from irns1 up to irws1 ;
c     remember that the lm = 1 contribution is multiplied
c     by a factor 1/sqrt(4 pi)
c     
             WRITE (IFILE,FMT=9060) NR_U,IRNS1,LMPOT,ISAVE
             WRITE (IFILE,FMT=9070) (VM2Z_U(IR),IR=1,NR_U)
             IF (LPOT.GT.0) THEN
                LMNR = 1
                DO 40 LM = 2,LMPOT
                   SUM = 0.0D0
                   DO 30 IR = IRMIN,NR_U
                      RV = VINS_U(IR,LM)*RR_U(IR)
                      SUM = SUM + RV*RV*DRADI(IR)
 30                CONTINUE
                   
                   IF (SQRT(SUM).GT.QBOUND) THEN
                      LMNR = LMNR + 1
                    WRITE (IFILE,FMT=9060) LM
                    WRITE (IFILE,FMT=9070) (VINS_U(IR,LM),IR=IRMIN,NR_U)
                   END IF
                   
 40             CONTINUE
c     
c---  >         write a one to mark the end
c     
                IF (LMNR.LT.LMPOT) WRITE (IFILE,FMT=9060) ISAVE
             END IF
             
          END IF
          
 50    CONTINUE
 60   CONTINUE
      
      
 8000 format (' GENERAL POTENTIAL MESH             exc:')
 8010 Format ('#  ',A4,'POTENTIAL             Z = ',F8.3)
 8011 format ('#  ',A4,'POTENTIAL SPIN UP     Z=  ',F8.3)
 8012 format ('#  ',A4,'POTENTIAL SPIN DOWN   Z=  ',F8.3)
 8020 format ('#')
 8030 format (4f12.8, '   # alat, rmax, rmaxlog, rmt')
 8040 format (1p,3I6,31X,'  # IRWS, IMR, IRNS ')
 8050 format (2D15.8,19X,'  # A , B ')
 8060 format (3f12.8,13X,'  # Ef, vbc, vcon ')
 8070 format (1p,2I5,39X,'  # NCORE, LMPOT' )
 9000 FORMAT (7a4,6x,'  exc:',a24,3x,a10)
 9010 FORMAT (3f12.8)
 9020 FORMAT (f10.5,/,f10.5,2f15.10)
 9030 FORMAT (i5,/,2d15.8,/,2i2)
 9040 FORMAT (i5,1p,d20.11)
 9050 FORMAT (1p,2d15.6,1p,d15.8)
 9051 FORMAT (1p,4d20.12)
 9060 FORMAT (10i5)
 9070 FORMAT (1p,4d20.13)
      END
