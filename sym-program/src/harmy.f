      SUBROUTINE HARMY(R,IFILE,MTAOS,NATOM,NHOST,NIMP,L1,NBASUM)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (LMAX=48)
      PARAMETER (N48=48)
      PARAMETER (NBASMX=100)
      PARAMETER (LMAXP1=LMAX+1)
      DIMENSION COEF(N48,LMAXP1,2,NBASMX,3)
      DIMENSION ALPHA(N48),BETA(N48),GAMMA(N48),IDEN(N48)
      DIMENSION R(3,N48)
      DIMENSION DPL(LMAXP1,LMAXP1,8),DMN(LMAXP1,LMAXP1,8),DSRY(N48,3)
      DIMENSION COSN(LMAXP1,12),SINN(LMAXP1,12)
      DIMENSION MT(N48,N48),MALPH(N48),MBETA(N48),MGAM(N48)
      LOGICAL CLEAR
      DIMENSION ALPHAO(12),BETAO(8),TITLE(20)
      CHARACTER*4 WORD(2),NGRP,NREP,NGO
      DIMENSION NATOM(1),ICS(2),NM(N48)
      DIMENSION NHOST(1),NIMP(1),NBASUM(1)
      DATA ICS/1,-1/
      DATA WORD/' COS',' SIN'/
      DATA NGRP,NREP/'GRP ','REP '/
      JREP=0
      WRITE(14,*) IFILE
      REWIND IFILE
      PI=4.D0*DATAN(1.D0)
      SQR2=DSQRT(2.D0)
10    READ(12,360) KTAU,MTAU,INVER,NOORTH,IREP,IPUN
      WRITE(14,360) KTAU,MTAU,INVER,NOORTH,IREP,IPUN
      MTAOS=MTAU
      NBETA=0
      NALPH=0
      READ(12,501) ((R(I,M),I=1,3),NHOST(M),NATOM(M),NIMP(M),M=1,MTAU)
 501  FORMAT(3F10.6,I2,I4,I5)
      WRITE (14,500) ((R(I,M),I=1,3),NATOM(M),M=1,MTAU)
      WRITE (14,380)
C---------------------------------------------------------------------
C THE LOOP KT=1,KTAU STORES ALL DIFFERENT EULER ANGLES IN THE ARRAYS
C ALPHAO AND BETAO AND THE CORRESPONDING ADDRESSES IN ARRAYS MALPH,
C MBETA,AND MGAM:
C                ALPHAO(MALPH(KT)) = ALPHA(KT)
C                ALPHAO(MGAM(KT))  = GAMMA(KT)
C                BETAO(MBETA(KT))  = BETA(KT)
C---------------------------------------------------------------------
      REWIND 18
      DO 90 KT=1,KTAU
      READ(18,370) ID,OM,EGA,X,Y,Z
      IDEN(KT)=ID
      OMEGA=OM/EGA
C---------------------------------------------------------------------
      CALL EULER(KT,
     1          ALPHA,BETA,GAMMA,OMEGA,X,Y,Z)
C---------------------------------------------------------------------
      WRITE(14,390) ID,X,Y,Z,OMEGA,ALPHA(KT),BETA(KT),GAMMA(KT)
      ALPHA(KT)=PI*ALPHA(KT)
      BETA(KT)=BETA(KT)*PI
      GAMMA(KT)=PI*GAMMA(KT)
C---------------------------------------------------------------------
      CALL PERMUT(R,KT,
     1          ALPHA,BETA,GAMMA,
     2          MT,MTAU,INVER)
C---------------------------------------------------------------------
      WRITE(14,400) (MT(KT,M),M=1,MTAU)
      IF (NALPH .EQ. 0) GO TO 30
      DO 20 NK=1,NALPH
      IF (DABS(ALPHAO(NK)-ALPHA(KT)) .GT. 1.0D-4) GO TO 20
      MALPH(KT)=NK
      GO TO 40
20    CONTINUE
30    NALPH=NALPH+1
      MALPH(KT)=NALPH
      ALPHAO(NALPH)=ALPHA(KT)
40    DO 50 NK=1,NALPH
      IF (DABS(ALPHAO(NK)-GAMMA(KT)) .GT. 1.0D-4) GO TO 50
      MGAM(KT)=NK
      GO TO 60
50    CONTINUE
      NALPH=NALPH+1
      MGAM(KT)=NALPH
      ALPHAO(NALPH)=GAMMA(KT)
60    IF (NBETA .EQ. 0) GO TO 80
      DO 70 NK=1,NBETA
      IF (DABS(BETAO(NK)-BETA(KT)) .GT. 1.0D-4) GO TO 70
      MBETA(KT)=NK
      GO TO 90
70    CONTINUE
80    NBETA=NBETA+1
      MBETA(KT)=NBETA
      BETAO(NBETA)=BETA(KT)
90    CONTINUE
      IF (NALPH .GT. 12 .OR. NBETA .GT. 8) WRITE (14,410) NALPH,NBETA
      CALL TMX(IREP,
     1          ALPHA,BETA,GAMMA,IDEN
     2         ,KTAU,INVER)
      M=0
100   DO 110 NK=1,NALPH
      COSN(M+1,NK)=DCOS(M*ALPHAO(NK))
110   SINN(M+1,NK)=DSIN(M*ALPHAO(NK))
      M=M+1
      IF (M .LE. L1) GO TO 100
C-----------------------------------------------------------------------
C GO TO NEXT REPRESENTATION
C-----------------------------------------------------------------------
120   READ(18,359) NDIM
      JREP=JREP+1
      NBASOS=0
      REWIND 17
      IF (NDIM .EQ. 0) NDIM=1
      DO 130 NCOL=1,NDIM
      READ(18,420)TITLE
      WRITE (14,430)TITLE
      READ(18,440) (DSRY(KT,NCOL),KT=1,KTAU)
      DO 5441 KT=1,KTAU
      IF(DABS(DSRY(KT,NCOL)-DSQRT(0.75D0)).LT.1.D-6)
     +     DSRY(KT,NCOL)=DSQRT(0.75D0)
      IF(DABS(DSRY(KT,NCOL)+DSQRT(0.75D0)).LT.1.D-6)
     +     DSRY(KT,NCOL)=-DSQRT(0.75D0)
 5441 CONTINUE
130   WRITE(14,450)NCOL,(DSRY(KT,NCOL),KT=1,KTAU)
      NDG=NDIM*2
      NDGOS=NDG
      DO 350 LP1=1,L1+1
C-----------------------------------------------------------------------
C CONSTRUCTION OF LINEAR COMBINATIONS:
C
C                  L                  L
C            FAC*{D   (BETA)+(-1)**M*D    (BETA)}
C                  M'M                M'-M
C
C                  L                  L
C AND        FAC*{D   (BETA)-(-1)**M*D    (BETA)}
C                  M'M                M'-M
C
C TO BE USED IN SUBROUTINE PRJECT IN THE APPLICATION OF THE PROJECTION
C OPERATOR TO THE REAL SPHERICAL HARMONICS
C
C-----------------------------------------------------------------------
      LPX=LP1
      L=LP1-1
      FAC2=1.D0
      M=0
140   FAC1=1.D0
      MP=0
150   CONTINUE
      FAC=FAC1*FAC2/2.D0
      DO 170 NK=1,NBETA
      D1=DROTB(L,MP,M,BETAO(NK))
      D2=DROTB(L,MP,-M,BETAO(NK))
      IF (MOD(M,2) .NE. 0) D2=-D2
      DPL(MP+1,M+1,NK)=(D1+D2)*FAC
      DMN(MP+1,M+1,NK)=(D1-D2)*FAC
      IF (MOD(M+MP,2) .NE. 0) GO TO 160
      DPL(M+1,MP+1,NK)=DPL(MP+1,M+1,NK)
      DMN(M+1,MP+1,NK)=DMN(MP+1,M+1,NK)
      GO TO 170
160   DMN(M+1,MP+1,NK)=-DMN(MP+1,M+1,NK)
      DPL(M+1,MP+1,NK)=-DPL(MP+1,M+1,NK)
170   CONTINUE
      FAC1=SQR2
      MP=1+MP
      IF (MP .LE. M) GO TO 150
      FAC2=SQR2
      M=1+M
      IF (M .LE. L) GO TO 140
      NBAS=1
C---------------------------------------------------------------------
C
C CALCULATION OF      S                      FOR J0=1
C                      LMN,J0,LAMBDA,MU(L)
C---------------------------------------------------------------------
      DO 210 MTO=1,MTAU
      IMAX=1
      M=0
180   CONTINUE
      DO 200 I=1,IMAX
C---------------------------------------------------------------------
      CALL PRJECT(M,I,MTO,NBAS,1,1.0D0,.TRUE.,DSRY,DPL,DMN
     1           ,COSN,SINN,COEF
     2           ,MT,MALPH,MBETA,MGAM,MTAU,KTAU,L,LPX,INVER)
C---------------------------------------------------------------------
      IF (NOORTH .NE. 0) GO TO 190
      CALL ORTHOG(NBAS,COEF,MTAU,LPX)
      GO TO 200
190   NBAS=NBAS+1
200   CONTINUE
      IMAX=2
      M=M+1
      IF(M .LE. L) GO TO 180
210   CONTINUE
      NBAS=NBAS-1
      WRITE (14,510) NDIM,NBAS
      IF(NBAS .EQ. 0) GO TO 340
      IF (NDIM .EQ. 1) GO TO 290
C---------------------------------------------------------------------
C
C CALCULATION OF      S                      FOR J'.NE.J0
C                      LMN,J',LAMBDA,MU(L)
C---------------------------------------------------------------------
      DO 280 NCOL=2,NDIM
      IF(NBAS.GT.NBASMX) STOP 77
      DO 280 NB=1,NBAS
      CLEAR=.TRUE.
      IMAX=1
      DO 240 MTO=1,MTAU
      M=0
220   CONTINUE
      DO 230 I=1,IMAX
      AC=COEF(MTO,M+1,I,NB,1)
      IF (DABS(AC) .LT. 1.0D-4) GO TO 230
C---------------------------------------------------------------------
      CALL PRJECT(M,I,MTO,NB,NCOL,AC,CLEAR,DSRY,DPL,DMN
     1           ,COSN,SINN,COEF
     2           ,MT,MALPH,MBETA,MGAM,MTAU,KTAU,L,LPX,INVER)
C---------------------------------------------------------------------
      CLEAR=.FALSE.
230   CONTINUE
      IMAX=2
      M=1+M
      IF (M .LE. L) GO TO 220
240   CONTINUE
      FAC=NDIM
      FAC=FAC/KTAU
      DO 270 MTO=1,MTAU
      IMAX=1
      DO 260 M1=1,LP1
      DO 250 I=1,IMAX
      COEF(MTO,M1,I,NB,NCOL)=COEF(MTO,M1,I,NB,NCOL)*FAC
250   CONTINUE
260   IMAX=2
270   CONTINUE
280   CONTINUE
290   NDG=NDIM+NDIM
C
      NBASOS=NBASOS+NBAS
C
      IF(NBAS.GT.NBASMX) STOP 77
      DO 330 NB=1,NBAS
      DO 330 NCOL=1,NDIM
      IF (NDIM .NE. 1) WRITE (14,460) L,NCOL
      IF (NDIM .EQ. 1) WRITE (14,470) L
      M=0
300   CONTINUE
      IMAX=2
      IF (M .EQ. 0) IMAX=1
      DO 320 I=1,IMAX
      SUM=0.0D0
      DO 310 MTO=1,MTAU
      IF (DABS(COEF(MTO,M+1,I,NB,NCOL)).LT.1.D-4)
     1     COEF(MTO,M+1,I,NB,NCOL) =0.0D0
310   SUM=SUM+DABS(COEF(MTO,M+1,I,NB,NCOL))
      IF (SUM .LT. 1.D-4) GO TO 320
      WRITE (14,480) M,WORD(I),(COEF(MTO,M+1,I,NB,NCOL)
     1     ,NATOM(MTO),MTO=1,MTAU)
320   CONTINUE
      M=1+M
      IF(M .LE. L) GO TO 300
      IF(NCOL.GT.IPUN) GO TO 330
      NMS=0
      DO 322 MTO=1,MTAU
      NM(MTO)=0
      DO 321 M1=1,LP1
      M=M1-1
      IMAX=1
      IF(M.GT.0) IMAX=2
      DO 321 I=1,IMAX
      IF(DABS(COEF(MTO,M1,I,NB,NCOL)).LT.1.D-4) GO TO 321
      NM(MTO)=NM(MTO)+1
      NMS=NMS+1
321   CONTINUE
      IF(NM(MTO).EQ.0) NMS=NMS+1
322   CONTINUE
      IF (NDIM .NE. 1) WRITE(17,361) L,NMS,NCOL
      IF (NDIM .EQ. 1) WRITE(17,362) L,NMS
      DO 324 MTO=1,MTAU
      IF(NM(MTO).GT.0) GO TO 323
      MX=0
      IC=1
      NA=NATOM(MTO)
      C=0.D0
      WRITE(17,520) MX,IC,NA,C
      GO TO 324
323   DO 325 M1=1,LP1
      M=M1-1
      IMAX=1
      IF(M.GT.0) IMAX=2
      DO 325 I=1,IMAX
      IF(DABS(COEF(MTO,M1,I,NB,NCOL)).LT.1.D-4) GO TO 325
      WRITE(17,520) M,ICS(I),NATOM(MTO),COEF(MTO,M1,I,NB,NCOL)
325   CONTINUE
324   CONTINUE
330   CONTINUE
      GO TO 350
340   WRITE (14,490)L
350   CONTINUE
C---------------------------------------------------------------------
C     READING OS-DATA
C---------------------------------------------------------------------
      REWIND 17
      IF(IPUN.NE.0) WRITE(IFILE,5300) NBASOS,NDGOS,(TITLE(K),K=1,2)
      NDGOS=NDGOS/2
      NNOS=NBASOS*NDGOS
      IF(NBASOS.EQ.0) GO TO 3270
      DO 3260 J=1,NNOS
      IF(NDIM.NE.1) READ(17,361) LOS,NMSOS,NCOLOS
      IF(NDIM.EQ.1) READ(17,362) LOS,NMSOS
      IF(NDIM.EQ.1) NCOLOS=0
      WRITE(IFILE,361) LOS,NMSOS,NCOLOS
      DO 3240 I=1,NMSOS
      READ(17,520) MOS,ICSOS,NATOS,COEFOS
      WRITE(IFILE,520) MOS,ICSOS,NATOS,COEFOS
3240  CONTINUE
3260  CONTINUE
3270  CONTINUE
      NBASUM(JREP)=NBASUM(JREP)+NBASOS
      READ(18,420) NGO
      IF (NGO .EQ. NGRP) GO TO 10
      IF (NGO .EQ. NREP) GO TO 120
      RETURN
359   FORMAT(11I5)
360   FORMAT(5X,3I5,15X,3I5)
361   FORMAT(3I5)
362   FORMAT(2I5)
370   FORMAT (1X,A4,5X,F3.0,3X,F4.3,3F10.7)
380   FORMAT (//' OPERATION',7X,'AXIS OF ROTATION',8X,'ANGLE/PI',7X,'(EU
     LLER ANGLES)/PI',17X,'PERMUTATION')
390   FORMAT(1X,A4,5X,7F10.4,8I5)
400   FORMAT(1H+,80X,8I5/(81X,8I5))
410   FORMAT('  NALPH=',I3,'  NBETA=',I3,'  CHECK DIMENSIONS OF',' DPL,
     DDMN, COSN, SINN, ALPHAO, AND BETAO.')
420   FORMAT(20A4)
430   FORMAT(//20A4)
440   FORMAT(8F10.0)
450   FORMAT(60X,'DSRY FOR COL',I3,/,(1X,8F14.8))
460   FORMAT (/50X,'L=',I3,'   COLUMN',I3)
470   FORMAT(/60X,'L=',I3)
480   FORMAT(' M=',I3,A4,':',8(F12.8,I3)/(11X,8(F12.8,I3)))
490   FORMAT (//40X,'NO  BASIS  FOUND  FOR  L= ',I3/)
500   FORMAT(3F10.6,1X,I4)
510   FORMAT (//40X,'NDIM=',I8,27X,'NBAS=',I8)
520   FORMAT (3I5,F17.14)
5300  FORMAT (I5,I5,6X,8A4)
      END
