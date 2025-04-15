      SUBROUTINE SYMSTO(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,N0ATOM,
     +                  IMIN,IMAX,NAT,NSTR,NISTR,NMSTR,NSEC,NUMCOL,
     +                  NBASIS,NREP,NMS,LN,LMAX,NLEQ,N0,NTERMS,NATMX,
     +                  NDIM,NDG)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER LMAX,NAT,NATMX,NBASIS,NISTR,NMSTR,NREP,NSEC,NSTR,NUMCOL
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,NSEC,*),STR(*)
      INTEGER IMAX(NSEC,NATMX,NUMCOL),IMIN(NSEC,NATMX,NUMCOL),ISTR(*),
     +        LN(*),MAXSTR(*),MI(NBASIS,NSEC,*),MINSTR(*),MSTR(*),N0(*),
     +        N0ATOM(NBASIS,NSEC,*),NDG(*),NDIM(*),NFM(*),NLEQ(*),
     +        NMS(NSEC,*),NTERMS(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IFL,IFM,IFX,LAM,MIMAX,MIMIN,MPL,N,NCOL,NDIMNP,NF,NK,NM,
     +        NN,NP
C     ..
C     .. External Subroutines ..
      EXTERNAL I4INIT,IGATHER,SETUP,WHENNE
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX0,MIN0
C     ..
      IFL = 1
      IFX = 1
      IFM = 1
      NF = NSEC*NATMX*NUMCOL
      N = 0
      MIMAX = 0
      MIMIN = 0
      DO 50 NP = 1,NREP
        CALL I4INIT(IMAX,NF,0)
        CALL SETUP(CN,MI,IMIN,IMAX,N0ATOM,NBASIS,NSEC,NUMCOL,NAT,LMAX,
     +             NMS,LN,NLEQ,N0,NTERMS,NP,NATMX,NDIM,NDG)
C     IF(NDIM(NP).NE.NTERMS(1)) GOTO 120
C     IF(NLEQ(1)) 110,120,110
C 110 WRITE(6,30)
C  30 FORMAT(/'   THIS REPRESENTATION GIVES NO CONTRIBUTION, AND IS NOT
C    SSTORED'//)
C     GO TO 311
C 120 N=N+1
        N = N + 1
        NCOL = NDG(NP)/2
        NDIMNP = NDIM(NP)
        DO 30 NN = 1,NDIMNP
          ISTR(IFX) = LN(NN)
          LAM = LN(NN)**2 + LN(NN) + 1
          IFX = IFX + 1
          DO 20 NK = 1,NCOL
            NM = NMS(NN,NK)
            ISTR(IFX) = NM
            IFX = IFX + 1
            DO 10 I = 1,NM
              MPL = MI(I,NN,NK) + LAM
              ISTR(IFX) = MPL
              MIMAX = MAX0(MIMAX,MPL)
              MIMIN = MIN0(MIMIN,MPL)
              IFX = IFX + 1
              STR(IFL) = CN(I,NN,NK)
              IFL = IFL + 1
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
        DO 40 I = 1,NAT
          ISTR(IFX) = N0(I)
          IFX = IFX + 1
          ISTR(IFX) = NTERMS(I)
          IFX = IFX + 1
          ISTR(IFX) = NLEQ(I)
          IFX = IFX + 1
   40   CONTINUE
        CALL WHENNE(NF,IMAX,1,0,MSTR(IFM),NFM(NP))
        CALL IGATHER(NFM(NP),MAXSTR(IFM),IMAX,MSTR(IFM))
        CALL IGATHER(NFM(NP),MINSTR(IFM),IMIN,MSTR(IFM))
        IFM = IFM + NFM(NP)
        IF (IFM.LE.NMSTR .AND. IFX.LE.NISTR .AND.
     +      IFL.LE.NSTR) GO TO 50
        WRITE (6,FMT=9000) IFL,IFX,IFM
        STOP
   50 CONTINUE
      WRITE (6,FMT=9010) MIMIN,MIMAX
      NREP = N
      IF (NREP) 70,60,70
   60 WRITE (6,FMT=9020)
      STOP
   70 WRITE (6,FMT=9000) IFL,IFX,IFM
      RETURN
 9000 FORMAT (' MAXIMUM TEMPORARY STORAGE SUBSCRIPTS ARE ',3I8,/,/)
 9010 FORMAT (' MINIMUM AND MAXIMUM IN ARRAY MI:',2I5)
 9020 FORMAT ('  THERE IS NO CONTRIBUTING REPRESENTATION LEFT')
      END
