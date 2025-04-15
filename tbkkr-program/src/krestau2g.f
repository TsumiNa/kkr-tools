c ************************************************************************
      SUBROUTINE KRESTAU2G(KPOI,WEIGHT,GSK,TINVLL,NSHELL,NSH1,NSH2,
     &                     RFCTOR,RBASIS,IHANDLE1,IATCONDL,IATCONDR,
     &                     NCONDPAIR,IEGFOUT,IGF)
      implicit none
c ************************************************************************
c
c     GLL0 = GMATLL in subroutine GLL96
c     This prepares k-resolved Green's function from the tau matrix
c     To be used in the current calculation in the Landauer formalism.
c     24.02.2000
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.D0,0.D0))
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER NCLUSTD,NSIZE
      PARAMETER (NSIZE=NATOMIMPD*LMAXSQ,NCLUSTD=NSIZE*(NSIZE+1)/2)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION RFCTOR,WEIGHT
      INTEGER IGF,NSYMAT,NSHELL,NSH1(*),NSH2(*),IHANDLE1,i1,j1
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX GSK(LMAXSQ,LMAXSQ,*),TINVLL(LMAXSQ,LMAXSQ,*)
	complex, allocatable ::  GCLUST(:)
      DOUBLE PRECISION RBASIS(3,*),KPOI(3)
      INTEGER IATCONDL(*),IATCONDR(*),NCONDPAIR,IEGFOUT
C     ..
C     .. Local Scalars ..
      INTEGER I,IND,IU,LM1,LM2,NSLL,NS,J,LIN,NSHELL0,ns0
      LOGICAL LDIA
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX GLL(LMAXSQ,LMAXSQ),
     +               TPG(LMAXSQ,LMAXSQ),
     +               XC(LMAXSQ,LMAXSQ)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CINIT,ZCOPY,ZSCAL,ZGEMM,TEST
      INTRINSIC ABS,ZABS
c
C     .. Save statement ..
      SAVE
c ------------------------------------------------------------------------
C     ..
c
c  I am writing out k-resolved Green's function for all 
c  shells except the diagonal.
c
c
      !write(6,*) 'in krestau2g',ncondpair,ihandle1
      allocate (GCLUST(nclustd))
c      On Win OS whatever the os is x86 or x64, Windows limits static code and data to 2GB, The solution is to change the arrays from being declared with fixed bounds to being allocatable, and then using ALLOCATE to make them the desired size. 

      
      IF (NCONDPAIR.GT.0) THEN
         DO I=1,3
c           CALL FXDRDBL(IHANDLE1,KPOI(i),1)
         END DO
c        CALL FXDRDBL (IHANDLE1,WEIGHT,1)
      END IF

      DO 70 NS = 1,NCONDPAIR

         I = IATCONDL(NS) 
         J = IATCONDR(NS) 
c         WRITE(6,*) 'WRITINIG PAIR :',I,J 
c
c Now find NSHELL0 
c         
         NSHELL0 = 0
         DO NS0 = 1,NSHELL
           IF ( (I.EQ.NSH1(NS0)) .AND. (J.EQ.NSH2(NS0)) ) THEN
           NSHELL0 = NS0
           EnD IF
         EnD DO
c         WRITE (6,*) 'NOW NSHELL :',nshell0
         !do lm1=1,3
         ! write(6,*) 'in green :',lm1,gsk(lm1,lm1,nshell0)
         !end do
         IF (NSHELL0.EQ.0) THEN
             WRITE(6,*) 'THIS PAIR IS NOT CALCULATED LOOK '
             WRITE(6,*) 'AT SUB GFSHELLS '
             STOP 'KRESTAU2G is STOPPING'
         END IF 
c      
         LDIA = (DABS  ((RBASIS(1,J)-RBASIS(1,I) )**2 +
     +        (RBASIS(2,J)-RBASIS(2,I) )**2 +
     +        (RBASIS(3,J)-RBASIS(3,I) )**2) .LT. 1.D-6 )
         
c     
c     
         CALL ZCOPY(LMAXSQ*LMAXSQ,GSK(1,1,NSHELL0),1,GLL,1)
c     
c     --->  XC = TINVLL(I) * GLL
c                                                                I is before 
         CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TINVLL(1,1,I),
     +        LMAXSQ,GLL,LMAXSQ,CZERO,XC,LMAXSQ)


         IF (LDIA) THEN
c     
c     --->    GLL = -TINVLL - TINVLL * GLL* TINVLL
c     
            CALL ZCOPY(LMAXSQ**2,TINVLL(1,1,I),1,GLL,1)
            CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +           TINVLL(1,1,I),LMAXSQ,-CONE,GLL,LMAXSQ)
          write(6,*) ' Diagonal element of GF in Conductance !!!!'  
         ELSE                   ! (LDIA)
c     
c     --->    GLL =  - TINVLL(I) * GLL* TINVLL(J)
c     
            CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +           TINVLL(1,1,J),LMAXSQ,CZERO,GLL,LMAXSQ)
                  !         J is before
         END IF                 ! (LDIA)
c     
c     --->  GLL0 = GLL/RFCTOR
c     
c         DO 30 LM1 = 1,LMAXSQ
c            DO 20 LM2 = 1,LMAXSQ
c               GLL0(LM2,LM1) = GLL(LM2,LM1)/RFCTOR
c 20         CONTINUE
c 30      CONTINUE
c     
c     Now write out the k-resolved Green's function
c     Map in a linear array and write out
c     
         !do lm1=1,3
         !  write(6,*) 'Green ',lm1,gll(lm1,lm1)/rfctor
         !end do
         LIN = 0
         DO LM2=1,LMAXSQ
            DO LM1=1,LMAXSQ
               LIN = LIN + 1
               GCLUST(LIN) = GLL(LM1,LM2)/RFCTOR
            ENDDO
         ENDDO
         IF (LIN.GT.NCLUSTD) THEN
            WRITE(6,*) ' Error in SUBROUTINE KRESTAU2G '
            WRITE(6,*) 'LIN = ',LIN,' NCLUSTD= ',NCLUSTD
            WRITE(6,*) ' The GCLUST array is overflown !!'
            STOP
         END IF 
c        CALL FXDRINT(IHANDLE1,I,1)
c        CALL FXDRINT(IHANDLE1,J,1)            
c        CALL FXDRRL(IHANDLE1,GCLUST,2*LIN) 
         
 70   CONTINUE                  ! NS = 1,NPAIRS
      deallocate(GCLUST)
      RETURN
 9000 format(2f18.10)

      END
