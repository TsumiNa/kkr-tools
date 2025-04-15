*DECK SPHERE
CDECK SPHERE
      subroutine SPHERE(LMAX,YR,WTYR,RIJ,IJEND,IJD)
C-----------------------------------------------------------------------
C     GENERATE A ANGULAR MESH AND THE SPHERICAL HARMONICS AT THOSE
C     MESH POINTS . FOR AN ANGULAR INTEGRATION THE WEIGHTS ARE GE-
C     RATED .
C     THE SPHERICAL HARMONICS AND THE SPHERICAL HARMONICS TIMES THE
C     WEIGHTS GENERATED ON THE UNIT SPHERE , THE POINTS AND THE
C     WEIGHTS FOR THE GAUSS-LEGENDRE INTEGRATION ARE STORED IN THE
C     COMMON BLOCK RSPHERE .
C     FOR A BETTER VECTORIZATION COMBINED INDICES ARE INTRODUCED .
C
C                                  B.DRITTLER    MAY 1987
C
C
C     USING GAUSSIAN QUADRATURE AS GIVEN BY
C     M. ABRAMOWITZ AND I.A. STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS,
C     NBS APPLIED MATHEMATICS SERIES 55 (1968), PAGES 887 AND 916
C     M. WEINERT AND E. WIMMER
C     NORTHWESTERN UNIVERSITY MARCH 1980
C
C
C     MODIFIED TO USE CALCULATED POINTS AND WEIGHTS
C     TO MAKE IT DYNAMIC.   (M.W.  JAN. 1982)
C-----------------------------------------------------------------------
C     .. PARAMETERS ..
      implicit real*8(a-h,o-z)
	include 'inc.fi'
c      INTEGER LPOTD
c      PARAMETER(LPOTD=8)
      INTEGER LMMAXD,N1,IJD,IJDM
C     PARAMETER (LMMAXD= (LPOTD+1)**2,N1=2* (LPOTD+1),IJD=2*N1**2)
      PARAMETER (LMMAXD= (LPOTD+1)**2,N1=2* (LPOTD+1),IJDM=434)
C     ..
C     .. SCALAR ARGUMENTS ..
      INTEGER LMAX,LMAX3
C     ..
C     .. SCALARS IN COMMON ..
      INTEGER IJEND
C     ..
C     .. ARRAYS IN COMMON ..
      REAL*8 RIJ(IJD,*),W(N1),WTYR(IJD,*),X(N1),YR(IJD,*)
      REAL*8 Y(LMMAXD),WEI
C     ..
C     .. LOCAL SCALARS ..
      REAL*8 FAC,HW,PI,R,R1,R2,R3,WJ
      INTEGER I,IJ,IS,J,LM1,NN,III
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL GRULE,YMY
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC COS,REAL,SIN,SQRT
C     ..
C     .. COMMON BLOCKS ..
C
C     COMMON /RSPHERE/YR,WTYR,RIJ,X,W,IJEND
      COMMON /LEBEDEV/ICHECK(IJDM)

C
c.....------------------------------------------------------------------
C     dimension cosx(IJD),fai(IJD)
      real*8 cosx(IJDM),fai(IJDM)
c     double dimension cosx(IJDM),fai(IJDM)
c.....
c      common/der2/np,cosx,fai                                

C     ..
C     .. SAVE STATEMENT ..
      SAVE
C     ..
C     .. DATA STATEMENTS ..
      DATA PI/3.14159265359E0/
C     ..
C
c      open(89,FORM='unformatted')
c      rewind 89
      WRITE(6,*) ' FORT.89 in sphee ijend lmax',ijend,lmax
      READ(89) IJEND
      WRITE(6,*) IJEND
      IF(IJEND.GT.IJD) STOP 'SPHERE'
C
      DO 30 IJ=1,IJEND
      READ(89) R1,R2,R3,WEI,COSX(IJ),FAI(IJ)
      WRITE(6,9005) IJEND,IJ,R1,R2,R3,WEI,COSX(IJ),FAI(IJ)
      RIJ(IJ,1)=R1
      RIJ(IJ,2)=R2
      RIJ(IJ,3)=R3
C
      CALL YMY(R1,R2,R3,R,Y,LMAX) 
C
      DO 11 LM1=1,(LMAX+1)**2
   11 YR(IJ,LM1) = Y(LM1)
C
C
      DO 12 LM1=1,(LMAX+1)**2
   12 WTYR(IJ,LM1)=YR(IJ,LM1)*WEI*PI*4.E0
C
   30 CONTINUE
C
      DO 31 IJ=1,IJEND
   31 READ(89) III,ICHECK(IJ)

      DO 32 IJ=1,IJEND
C
C
      IF(ICHECK(IJ).NE.0) write(6,*) ' IJ ICHECK = ',IJ,ICHECK(IJ)
   32  CONTINUE
C     STOP ' SPHERE'
      np=ijend
c     write(6,*) 'cosx'
c     write(6,9011) np
 9011 format(i5) 
c     write(6,9010) (cosx(i),i=1,ijend)
 9010 format(1x,10f10.6)
C
       write(6,9001) IJEND,LMAX
 9001 format(1x,' IJEND LMAX ',10I5) 
        Do 21 i=1,5
        aaaa=acos(cosx(i))
        write(6,9007) aaaa,cosx(i),fai(i)
 9007  format(1x, ' thet cosx fai',10F8.3)
   21  continue
C
      call cylm02(lmax,np,cosx,fai)
c
c
 9000 FORMAT (9D12.6)
 9005 FORMAT (1x,' IJEND IJ R1 R2 R3 W COSX FAI',2I6,10f10.5)
      END
