      SUBROUTINE SYMRCL(STR,ISTR,MSTR,MINSTR,MAXSTR,NFM,CN,MI,IMIN,IMAX,
     +                  IIFL,IIFX,IIFM,NAT,LN,NBASIS,NSEC,NUMCOL,NP,
     +                  NLEQ,N0,NTERMS,NATMX,NDIM,NDG)
      IMPLICIT NONE
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER IIFL,IIFM,IIFX,NAT,NATMX,NBASIS,NP,NSEC,NUMCOL
C     ..
C     .. Array Arguments ..
      REAL*8 CN(NBASIS,NSEC,*),STR(*)
      INTEGER IMAX(NSEC,NATMX,NUMCOL),IMIN(NSEC,NATMX,NUMCOL),ISTR(*),
     +        LN(*),MAXSTR(*),MI(NBASIS,NSEC,*),MINSTR(*),MSTR(*),N0(*),
     +        NDG(*),NDIM(*),NFM(*),NLEQ(*),NTERMS(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IFL,IFM,IFX,NCOL,NDIMNP,NF,NK,NM,NN
C     ..
C     .. External Subroutines ..
      EXTERNAL ISCATTER
C     ..
      IFL = IIFL
      IFX = IIFX
      IFM = IIFM
      NF = NSEC*NATMX*NUMCOL
      NCOL = NDG(NP)/2
      NDIMNP = NDIM(NP)
      DO 30 NN = 1,NDIMNP
        LN(NN) = ISTR(IFX)
        IFX = IFX + 1
        DO 20 NK = 1,NCOL
          NM = ISTR(IFX)
          IFX = IFX + 1
C THE NEXT LOOP IS REPLACED BY MORE EFFICIENT CODE
C     DO 3 I=1,NM
C     MI(I,NN,NK)=ISTR(IFX)
C     IFX=IFX+1
C     CN(I,NN,NK)=STR(IFL)
C   3 IFL=IFL+1
          DO 10 I = 1,NM
            MI(I,NN,NK) = ISTR(IFX-1+I)
            CN(I,NN,NK) = STR(IFL-1+I)
   10     CONTINUE
          IFX = IFX + NM
          IFL = IFL + NM
   20   CONTINUE
   30 CONTINUE
      DO 70 I = 1,NAT
        N0(I) = ISTR(IFX)
        IFX = IFX + 1
        NTERMS(I) = ISTR(IFX)
        IFX = IFX + 1
        NLEQ(I) = ISTR(IFX)
        IFX = IFX + 1
C THE NEXT TWO LOOPS ARE REPLACED BY MORE EFFICIENT CODE IF NCOL=1,2,3
C     DO 5 NN=1,NDIMNP
C     DO 5 NK=1,NCOL
C     IMIN(NN,I,NK)=ISTR(IFX)
C     IFX=IFX+1
C     IMAX(NN,I,NK)=ISTR(IFX)
C     IFX=IFX+1
C   5 CONTINUE
C     CALL R8INIT(IMAX,NF,0)
C     CALL R8INIT(IMIN,NF,0)
        IF (NCOL.EQ.1) THEN
          DO 40 NN = 1,NDIMNP
            IMIN(NN,I,1) = 0
            IMAX(NN,I,1) = 0
   40     CONTINUE
        END IF
        IF (NCOL.EQ.2) THEN
          DO 50 NN = 1,NDIMNP
            IMIN(NN,I,1) = 0
            IMAX(NN,I,1) = 0
            IMIN(NN,I,2) = 0
            IMAX(NN,I,2) = 0
   50     CONTINUE
        END IF
        IF (NCOL.EQ.3) THEN
          DO 60 NN = 1,NDIMNP
            IMIN(NN,I,1) = 0
            IMAX(NN,I,1) = 0
            IMIN(NN,I,2) = 0
            IMAX(NN,I,2) = 0
            IMIN(NN,I,3) = 0
            IMAX(NN,I,3) = 0
   60     CONTINUE
        END IF
   70 CONTINUE
      CALL ISCATTER(NFM(NP),IMAX,MSTR(IFM),MAXSTR(IFM))
      CALL ISCATTER(NFM(NP),IMIN,MSTR(IFM),MINSTR(IFM))
      IFM = IFM + NFM(NP)
C
      IIFL = IFL
      IIFX = IFX
      IIFM = IFM
      RETURN
      END
