
c 13.10.95 ***************************************************************
      SUBROUTINE KLOOPZ1(E,DF,GMATLL,TMATLL,DELTALL,
     +                   INS,ALAT,LMAXDD,LMMAX,MD,E2,NSPIN,
     +                   MAXMESH,NMESH,
     +                   IE,IELAST,LSTART,IGF,CLEB,ICLEB,LOFLM,IEND,
     +                   KSCOEF,NSHELL,INTERVY,INTERVZ,
     +                   BBYA,CBYA,NAEZ,NATYP,
     +                   CLS,EQINV,NACLS,RR,
     +                   RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
     +                   BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
     &                   LEFTTINVLL,RIGHTTINVLL,vacflag,NLBASIS,
     &                   NRBASIS,IHANDLE,IHANDLE1,ATOMIMP,
     &                   IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT,
     &                   GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM)
c ************************************************************************
      implicit none
C     .. Parameters ..
      include 'inc.fi'
      include 'inc.cls'
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
C     ..
      DOUBLE COMPLEX DF,E
      DOUBLE PRECISION ALAT,BBYA,CBYA,E2
      INTEGER ICC,IE,IELAST,IEND,IGF,INS,INTERVY,INTERVZ,
     +        IJEND,
     +        KSCOEF,
     +        LMAXDD,LMMAX,LATT,LPOT,
     +        MAXMESH,
     +        MD,
     +        NAEZ,NATYP,NMESH,NSPIN,NZ,NLBASIS,NRBASIS,IHANDLE,
     &        IHANDLE1,NCONDPAIR,IATCONDL(*),IATCONDR(*),IEGFOUT
      LOGICAL LSTART,VACFLAG(2)
C     ..
C     .. Array Arguments ..
C     ..
      DOUBLE COMPLEX DELTALL(LMAXSQ,LMAXSQ,*),
     +               GMATLL(LMAXSQ,LMAXSQ,*),
     +               GINP(LMAXSQ*NACLSD,LMAXSQ,*),
     +               TMATLL(LMAXSQ,LMAXSQ,*)
      DOUBLE COMPLEX LEFTTINVLL(LMAXSQ,LMAXSQ,*),
     &               RIGHTTINVLL(LMAXSQ,LMAXSQ,*)
      DOUBLE COMPLEX GSQDOS(LMAXSQ,LMAXSQ,NAEZD,*)
      DOUBLE PRECISION QDOSKP(3,*)
      INTEGER NQDOSKP,QDOSKDIM 
      DOUBLE PRECISION CLEB(NCLEB),
     +                 RBASIS(3,*),RR(3,*),RCLS(3,NACLSD,*),
     +                 BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION RIJ(IJD,3),WTYR(IJD,*),YR(IJD,*)
      INTEGER ATOM(NACLSD,*),
     +        CLS(*),
     +        EQINV(*),EZOA(NACLSD,*),
     +        ICLEB(NCLEB,4),
     +        LOFLM(*),
     +        KAOEZ(*),
     +        NACLS(*),
     +        NSHELL(0:NSHELD),ATOMIMP(NATOMIMPD)
C     ..
C     .. Local Scalars ..
C     ..
      DOUBLE COMPLEX CONE,CZERO,EZ
      DOUBLE PRECISION RFCTOR
      INTEGER IH,INFO,LM1,LM2
C     ..
C     .. Local Arrays ..
C     ..
      DOUBLE COMPLEX TINVLL(LMAXSQ,LMAXSQ,NAEZD),
     +               TLL(LMAXSQ,LMAXSQ)
      INTEGER IPVT(LMAXSQ)
C     ..
C     .. External Subroutines ..
      LOGICAL TEST,OPT                               ! opt adeed 29.10.99
      EXTERNAL GLL2k,TEST,ZGETRF,ZGETRS,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE / (1.0D0,0.0D0)/
c ------------------------------------------------------------------------
C
      IF(TEST('flow    ')) write(6,*) 
     +     '>>> KLOOPZ1: invert delta_t and call GLL96'

      EZ = E
C     RFCTOR=A/(2*PI) conversion factor to p.u.
      RFCTOR = ALAT/ (8.D0*DATAN(1.0D0))          ! = ALAT/(2*PI)
c
c ---> loop over all atoms in unit cell
c
      DO 70 IH = 1,NAEZ
c
c --->  init TINVLL
c
        DO 20 LM2 = 1,LMMAX
          DO 10 LM1 = 1,LMMAX
            TINVLL(LM1,LM2,IH) = CZERO
 10       CONTINUE
 20     CONTINUE

        IF (INS.LT.1) THEN

          DO 30 LM1 = 1,LMMAX
            IF(TMATLL(LM1,LM1,IH).NE.0) then
            TINVLL(LM1,LM1,IH) = CONE/TMATLL(LM1,LM1,IH)
            ELSE
             TINVLL(LM1,LM1,IH) =0
            END IF
c            write(112,*) IH,LM1,TMATLL(LM1,LM1,IH)
 30       CONTINUE

c          IF (TEST('tmat    ')) THEN
          IF (TEST('TINVLL  ')) THEN
            WRITE(6,*) 'TINVLL (',IH,' ) :'
            WRITE(6,fmt=*) 
c            WRITE(6,fmt='(1p,(d14.3,d11.3))') 
     +           (TINVLL(LM1,LM1,IH),LM1 = 1,LMMAX)
          END IF

        ELSE                        ! (INS.LT.1)

          DO 40 LM1 = 1,LMMAX
            TINVLL(LM1,LM1,IH) = CONE
 40       CONTINUE

          DO 60 LM1 = 1,LMMAX
            DO 50 LM2 = 1,LMMAX
              TLL(LM2,LM1) = TMATLL(LM2,LM1,IH)
 50         CONTINUE
 60       CONTINUE
c
c --->   invert t-matrix
c
          CALL ZGETRF(LMMAX,LMMAX,TLL,LMAXSQ,IPVT,INFO)
          CALL ZGETRS('N',LMMAX,LMMAX,TLL,LMAXSQ,
     +                IPVT,TINVLL(1,1,IH),LMAXSQ,INFO)

        END IF                      ! (INS.LT.1)

 70   CONTINUE                      ! IH = 1,NAEZ
c THIS PART IS ALL ADDED ON --------------------------    29.10.99
c  The inverse of the t-matrix is calculated, now write
c  this on a file in case of decimation output.
c
      IF(OPT('deci-out')) THEN

         WRITE(37,1160) IE,E,DF
         DO IH=1,NAEZ
            WRITE(37,1170) IH
            IF (INS.NE.0) THEN  ! Write full-pot t-mat
               WRITE(37,1180) ((TINVLL(LM1,LM2,IH),LM2=1,LMMAX),
     &              LM1=1,LMMAX)
            ELSE                ! write spherical tmat
               WRITE(37,1180) (TINVLL(LM1,LM1,IH),LM1=1,LMMAX)
            END IF
         END DO

      END IF 
c ----------------------------------------------------    29.10.99
      write(6,*) ' before gll2k'
      CALL GLL2K(LSTART,MD,TINVLL,TMATLL,GMATLL,DELTALL,
     +     EZ,DF,RFCTOR,E2,NSPIN,MAXMESH,NMESH,
     +     IE,IELAST,IGF,LOFLM,KSCOEF,NSHELL,INTERVY,INTERVZ,
     +     BBYA,CBYA,NAEZ,NATYP,CLS,EQINV,NACLS,RR,
     +     RBASIS,EZOA,ATOM,RCLS,KAOEZ,LATT,ICC,GINP,
     +     BRAVAIS,RECBV,LPOT,YR,WTYR,RIJ,IJEND,
     &     LEFTTINVLL,RIGHTTINVLL,vacflag,NLBASIS,NRBASIS,IHANDLE,
     &     IHANDLE1,ATOMIMP,IATCONDL,IATCONDR,NCONDPAIR,IEGFOUT,
     &     GSQDOS,QDOSKP,NQDOSKP,QDOSKDIM) 
c      WRITE(6,*) 'GMATLL kloopz1',GMATLL(1,1,1) 
      IF(TEST('flow    ')) write(6,*) '<<< KLOOPZ1'
      
      RETURN

 1160 FORMAT('ENERGY ',I5,4D16.8)
 1170 FORMAT('ATOM ',I3)
 1180 FORMAT(4D22.14)
c Read in Format
c 2210  FORMAT(A80//A80)
c 2220  FORMAT(A6,F9.6,A8,I2,A7,I3,A8,I3)
c 2230  FORMAT(3F8.4)
c 2240  FORMAT(A6,F8.4,A5,,F8.4,A7,,F8.4)
c 2250  FORMAT(A4,I3,A5,I3,A5,I3,A7,I3)
c 2260  FORMAT(2I5,4D16.8)
c 2270  FORMAT(A5,I3)
c 2280  FORMAT(4D22.14)
      END
