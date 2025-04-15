       SUBROUTINE CONDINPUT(IATCONDL,IATCONDR,NCPAIRD,NCONDPAIR,NAEZ,
     &                      RBASIS,ALAT,IEGFOUT)
      implicit none
c ****************************************************
c * This sub just decides which atoms contribute to the
c * conductivity.
c ****************************************************
      INTEGER NCPAIRD
      INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR
      INTEGER NAEZ,IEGFOUT
      DOUBLE PRECISION RBASIS(3,*),ALAT
      logical cartesian
      integer il,ier,ns,I,J
      double precision r1,r2
      double precision ZCURR1,ZCURR2,RMTCURR
      character*80 UIO
c
c
         IL = 1   
         CALL IoInput('ZCURR1    ',UIO,IL,7,IER)
         READ (UNIT=UIO,FMT=*) ZCURR1
         CALL IoInput('ZCURR2    ',UIO,IL,7,IER)
         READ (UNIT=UIO,FMT=*) ZCURR2
         CALL IoInput('RMTCURR   ',UIO,IL,7,IER)
         READ (UNIT=UIO,FMT=*) RMTCURR
         CALL IoInput('CARTESIAN ',UIO,IL,7,IER)
         READ (UNIT=UIO,FMT=*) cartesian
         CALL IoInput('IEGFOUT   ' ,UIO,IL,7,IER)
         READ (UNIT=UIO,FMT=*) IEGFOUT
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     Read System setup from the inputcard
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  
         ZCURR1 = ZCURR1*ALAT
         ZCURR2 = ZCURR2*ALAT
         RMTCURR = RMTCURR*ALAT
c     
c     We want to calculate the current between Z1, Z2 planes
c     we must find which atomic polyhedra are cutted by these
c     planes to determine the atomic pairs we have to use.
c     if the planes are outside our structure just write an error!
c     
c     Since we will use 2d-geometry the Z1,Z2 are just in units of
c     alat.
c     The mufin-tin radious RMT0 will be used for this. 
c     RBASIS are in units of alat
         ns = 0
         DO I=1,NAEZ
            R1 = RBASIS(3,I)*ALAT - Zcurr1
            IF (ABS(R1).LT.RMTCURR) THEN
               DO J=1,NAEZ
                  R2 = RBASIS(3,J)*ALAT - Zcurr2
                  IF (ABS(R2).LT.RMTCURR) THEN
                     ns = ns+1
                     Iatcondl(ns) = I
                     Iatcondr(ns) = J
                  end IF
               EnD DO
            END IF
         end do 
         write(6,*) 
         write(6,*) 
         write(6,*) ' CONDUCTANCE CALCULATION ' 
         write(6,*)            
         write(6,*) 'We have found ',ns,' PAIRS'
         IF (NS.LE.0) THEN
             WRITE(6,*) '>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<'
             WRITE(6,*) ' C H E C K   T H E   P L A I N S   '
             WRITE(6,*) '    *NO*  k-res GF output          '
             WRITE(6,*) '  ... PROGRAM  CONTINUES ...       '
             WRITE(6,*) '>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<'        

         END IF
         IF (NS.GT.NCPAIRD) then
         WRITE(6,*) '** Error in  SUBROUTINE CONDINPUT '
         WRITE(6,*) 'change NCPAIR in the main program also '
         WRITE(6,*) 'NCPAIR=',NCPAIRD,' must be bigger than ',ns
         STOP         
         END IF
         NCONDPAIR = NS
         WRITE(6,*) 'The junction extends from '
         write(6,*) '                  Z min = ',RBASIS(3,1)*alat
         write(6,*) '               to Z max = ',RBASIS(3,NAEZ)*alat
         write(6,*) 'The PLANE Z1= ',Zcurr1,' is passing through'
         do i=1,ns
            write (6,*) ' ATOM NO : ',Iatcondl(ns),'with z= ',
     &           RBASIS(3,Iatcondl(ns))
         end do
         write(6,*) 'The PLANE Z2= ',Zcurr2,' is passing through'
         do i=1,ns
            write (6,*) ' ATOM NO : ',IatcondR(ns),'with z= ',
     &           RBASIS(3,IatcondR(ns))
         end do
         write(6,*) 'The k-resolved GF for energy no:',IEGFOUT,
     &              ' will be writen out'
         write(6,*)
         write(6,*)  
         END        
