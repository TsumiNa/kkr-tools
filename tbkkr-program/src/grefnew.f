       SUBROUTINE GREFNEW(LSTART,EZ,CLEB,ICLEB,LOFLM,IEND,
     +                    RATOM,NATOM,ALAT,GREF1)
       implicit none
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
c
      INTEGER LMAX,NATOMD
      PARAMETER (LMAX=LMAXD,NATOMD=NACLSD)
      INTEGER LMAXSQ,NGD
      PARAMETER (LMAXSQ= (LMAX+1)**2,NGD=LMAXSQ*NATOMD)
      DOUBLE COMPLEX CONE,CZERO,CONEM,CTWO
      PARAMETER (CONE  = (1.D0,0.D0),
     +           CTWO  = (2.D0,0.D0),
     +           CZERO = (0.D0,0.D0),
     +           CONEM = (-1.D0,0.D0))
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EZ,EZ1
      INTEGER IEND,NATOM
      DOUBLE PRECISION ALAT
      LOGICAL LSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION CLEB(*),RATOM(3,*)
c
      DOUBLE COMPLEX GREF1(NACLSD*LMAXSQ,LMAXSQ)
      INTEGER ICLEB(NCLEB,*),LOFLM(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,L,LM,LM1,LM2,M,N,N1,N2,NDIM,NLM1,NLM2
      DOUBLE COMPLEX A,B
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX 
     +     GLL(LMAXSQ,LMAXSQ),
     +     GREF(NGD,NGD),
     +     GTREF(NGD,LMAXSQ)
c
      DOUBLE PRECISION RDIFF(3)
c
      INTEGER LF(144)
      DATA LF/0,3*1,5*2,7*3,9*4,11*5,13*6,15*7,17*8,19*9,21*10,23*11/

c
c ---> construct free Green's function
c
      EZ1 = EZ - 4.D0
      DO 90 N1 = 1,NATOM
      DO 80 N2 = 1,NATOM
          DO 30 I = 1,3
            RDIFF(I) = -(RATOM(I,N1) - RATOM(I,N2))*ALAT
 30       CONTINUE
          IF (N1.NE.N2) THEN
            CALL GFREE(RDIFF,EZ1,GLL,CLEB,ICLEB,LOFLM,IEND)
            DO 50 LM2 = 1,LMAXSQ
              NLM2 = (N2-1)*LMAXSQ + LM2
              DO 40 LM1 = 1,LMAXSQ
                NLM1 = (N1-1)*LMAXSQ + LM1
                GREF(NLM1,NLM2) = GLL(LM1,LM2)
 40           CONTINUE
 50         CONTINUE
          ELSE
            DO 70 LM2 = 1,LMAXSQ
              NLM2 = (N2-1)*LMAXSQ + LM2
              DO 60 LM1 = 1,LMAXSQ
                NLM1 = (N1-1)*LMAXSQ + LM1
                GREF(NLM1,NLM2) = CZERO
 60           CONTINUE
 70         CONTINUE
          END IF

 80     CONTINUE
 90   CONTINUE
      DO N1=1,NATOM
         DO LM1=1,LMAXSQ
             DO LM2=1,LMAXSQ
             NLM1 = (N1-1)*LMAXSQ + LM1
             GREF1(NLM1,Lm2) = GREF(NLM1,LM2)
             END do
         end do
      end do 
      RETURN 
      END
