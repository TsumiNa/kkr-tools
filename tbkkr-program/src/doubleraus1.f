c************************************************************************
      SUBROUTINE DOUBLERAUS1(IRMAX,IRMIN,LMPOT,RR,DRDI,VPOT,VINS)
c Gets rid of the double-points in the radial mesh, i.e. the points
c where RR(I) = RR(I-1). Returns the "new" mesh in the same array,
c and rearranges accordingly the WAVEF defined at the same mesh.
c IRMAX is also altered to the new value.
      implicit none
      INCLUDE 'inc.fi'
      INTEGER LMPOTD
      PARAMETER(LMPOTD=(2*LMAXD+1)**2)
      INTEGER IRMIND 
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NCOUNTMAX
      PARAMETER(NCOUNTMAX=500)
c Input and output:
      INTEGER IRMAX
      double precision RR(IRMD),DRDI(IRMD),VPOT(IRMD),
     &                 VINS(IRMIND:IRMD,LMPOTD)
c Inside:
      INTEGER IR,ICOUNT,NCOUNT,LMPOT,IRMIN
      INTEGER LM1,LM2
      INTEGER IDOUBLE(NCOUNTMAX)

c Find double-points:
      NCOUNT = 0
      DO IR = 2,IRMAX
         IF (ABS(RR(IR)-RR(IR-1)).LT.1.D-15) THEN
            NCOUNT = NCOUNT + 1
            IDOUBLE(NCOUNT) = IR
         ENDIF
      ENDDO

      IF (NCOUNT+1.GT.NCOUNTMAX) 
     &                STOP 'DOUBLERAUS2: Too many double-points.'
      IDOUBLE(NCOUNT+1) = IRMAX+1   ! To be used below.
      
c Rearrange the arrays.
      DO ICOUNT = 1,NCOUNT
         DO IR = IDOUBLE(ICOUNT)-ICOUNT+1,IDOUBLE(ICOUNT+1)-ICOUNT
            RR(IR) = RR(IR+ICOUNT)
            DRDI(IR) = DRDI(IR+ICOUNT)
            VPOT(IR) = VPOT(IR+ICOUNT)
	write(6,*) ' ir icount',IR,icount
         ENDDO
      ENDDO
      IRMAX = IRMAX - NCOUNT
       NCOUNT = 0
       IDOUBLE = 0 
      DO IR = IRMIN,IRMAX
         IF ((RR(IR)-RR(IR-1)).LT.1.D-20) THEN
            NCOUNT = NCOUNT + 1
            IDOUBLE(NCOUNT) = IR
         ENDIF
      ENDDO
c Rearrange the arrays.
      DO ICOUNT = 1,NCOUNT
         DO IR = IDOUBLE(ICOUNT)-ICOUNT+1,IDOUBLE(ICOUNT+1)-ICOUNT
            do lm1=1,lmpot
            VINS(IR,lm1) = VINS(IR+ICOUNT,lm1)
            end do
         ENDDO
      ENDDO      
      RETURN
      END
