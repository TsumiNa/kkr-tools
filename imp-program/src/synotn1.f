      SUBROUTINE SYNOTN1(DTMTRX,INS,LMAX,NATREF,NSTART,NEND,TMATLL,
     +                IRMIN,NQ,NDIM,NREP,TLLMAT,natper)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c
c     the delta t - matrix is symmetrized ( THIS IS CHANGED!!!)
c     this sub is for the case of no symmetry , it does a simple
c     maping of dtmtrx = tmatll
c     (transformed into the irreducible representation)
c
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER NATYPD
      PARAMETER (NATYPD=38)
      INTEGER NSEC,NREPD
      PARAMETER (nsec=689,NREPD=4)
      INTEGER LMAXD,LMX
      PARAMETER (lmaxd=4,LMX=LMAXD+1)
      INTEGER ICJD
      PARAMETER (ICJD=93139)
      INTEGER LMMAXD
      PARAMETER (LMMAXD=LMX**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER INS,LMAX,NATREF,NEND,NQ,NREP,NSTART
      integer lm1,lm2
      LOGICAL TLLMAT
C     ..
C     .. Array Arguments ..
      COMPLEX*16 DTMTRX(NSEC,NSEC),TMATLL(LMMAXD,LMMAXD,NATYPD)
      INTEGER IRMIN(*),NDIM(*)
C     ..
C     .. Local Scalars ..
      INTEGER I1,I2,ICJ,IRMIN1,IT,J,LA,LB,LMMAX,NP,natper,n1
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT
C     ..
C     ..

      LMMAX = (LMAX+1)* (LMAX+1)

c
c---> initialize dtmtrx
c
      CALL CINIT(NSEC*NSEC,DTMTRX)

         do i1=1,natper
          n1 = i1 + natref
          i2 = i1 -1 
          i2= i2*lmmax
          do lm1 =1,lmmax
            do lm2=1,lmmax
               dtmtrx(i2+lm1,i2+lm2) = tmatll(lm1,lm2,n1)
            end do
          end do
         end do
      RETURN


      END
