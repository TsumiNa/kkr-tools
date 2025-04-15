c ************************************************************************
      SUBROUTINE TAUTOG1(GS,TINVLL,DSYMLL,NSHELL,RFCTOR,GLL0,IGF,
     +                   TAUVBZ,NSYMAT,NSH1,NSH2,RATOM)
      implicit none
c ************************************************************************
c
c     GLL0 = GMATLL in subroutine GLL96
c
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
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX TAUVBZ
      DOUBLE PRECISION RFCTOR
      INTEGER IGF,NSYMAT,NSHELL,i1,j1
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX 
     +     GLL0(LMAXSQ,LMAXSQ,*),
     +     GS(LMAXSQ,LMAXSQ,NSYMAXD,*),
     +     TINVLL(LMAXSQ,LMAXSQ,*),
     +     DSYMLL(LMAXSQ,LMAXSQ,*)
      DOUBLE PRECISION RATOM(3,*)
      INTEGER NSH1(*),NSH2(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IND,IU,LM1,LM2,NSLL,NS,J
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
      DO 70 NS = 1,NSHELL
        I = NSH1(NS)
        J = NSH2(NS)
c
        LDIA = 
     +       (DABS(RATOM(1,NS)**2+
     +            RATOM(2,NS)**2+
     +            RATOM(3,NS)**2  ) .LT. 1.0D-6)
c


       
        DO 10 IU = 1,NSYMAT

         DO LM1=1,LMAXSQ
          DO LM2=1,LMAXSQ
c              write(106,FMT=8001) NS,IU,LM1,LM2,GS(LM1,LM2,IU,NS)
          END DO
        END DO
 8001   format(4I5,2F15.8)
c
c --->    GLL = sum(i=1,iumax)(tauvbz * DLL(i) * GS * DLL(i)^T)
c
          IF (IU.EQ.1) THEN
c
c --->      ull(1) is equal to unity matrix
c
            CALL ZCOPY(LMAXSQ*LMAXSQ,GS(1,1,1,NS),1,GLL,1)
c        write(6,*) ' after zcopy gs',gs(1,1,1,1)
            CALL ZSCAL(LMAXSQ*LMAXSQ,TAUVBZ,GLL,1)
C            CALL ZGEMM('N','T',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TPG,LMAXSQ,
C    +           ULL(1,1,IU),LMAXSQ,CZERO,GLL,LMAXSQ)
          ELSE
c
c --->      tpg = tauvbz * DLL * GS
c                         N
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,TAUVBZ,DSYMLL(1,1,IU),
     +           LMAXSQ,GS(1,1,IU,NS),LMAXSQ,CZERO,TPG,LMAXSQ)
c
c --->    GLL = GLL + TPG * DLL(i)^T
c                           T
            CALL ZGEMM('N','T',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TPG,LMAXSQ,
     +           DSYMLL(1,1,IU),LMAXSQ,CONE,GLL,LMAXSQ)
          END IF
c     
 10     CONTINUE                    ! IU = 1,NSYMAT
c
c --->  XC = TINVLL(I) * GLL
c

         DO LM1=1,LMAXSQ
          DO LM2=1,LMAXSQ
c              write(107,FMT=8000) NS,LM1,LM2,GLL(LM1,LM2)
c              write(6,FMT=8000) NS,LM1,LM2,GLL(LM1,LM2)
c              GLL(LM1,LM2)=GLL(LM1,LM2)*4
          END DO
        END DO
 8000   format('tautog1',3I5,2F15.8)

        CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,TINVLL(1,1,I),
     +       LMAXSQ,GLL,LMAXSQ,CZERO,XC,LMAXSQ)



        IF (LDIA) THEN
c     
c --->    GLL = -TINVLL - TINVLL * GLL* TINVLL
c     
          CALL ZCOPY(LMAXSQ**2,TINVLL(1,1,I),1,GLL,1)
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +         TINVLL(1,1,I),LMAXSQ,-CONE,GLL,LMAXSQ)

        ELSE                        ! (LDIA)
c     
c --->    GLL =  - TINVLL(I) * GLL* TINVLL(J)
c
          CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,-CONE,XC,LMAXSQ,
     +         TINVLL(1,1,J),LMAXSQ,CZERO,GLL,LMAXSQ)
          
        END IF                      ! (LDIA)
c
c --->  GLL0 = GLL/RFCTOR

c
        DO 30 LM1 = 1,LMAXSQ
          DO 20 LM2 = 1,LMAXSQ
            GLL0(LM2,LM1,NS) = GLL(LM2,LM1)/RFCTOR
 20       CONTINUE
 30     CONTINUE
c This was removed on 23.2.2000
c        IF (IGF.NE.0) THEN
c          WRITE (103,*) ((GLL0(LM2,LM1,NS),LM2=1,LMAXSQ),LM1=1,LMAXSQ)
c ------------------------------------------------------------------------
c          IF (TEST('GLL0    ')) THEN
c            write(6,*) 'GLL0 NS :',NS
c            write(6,*) 'LDIA :',ldia
c            DO LM1=1,LMAXSQ
c              DO LM2=LM1,LM1
c                IF ( zabs(GLL0(LM1,LM2,NS)).gt.1.d-10)
c     +               WRITE (6,FMT='(3i4,1p,2d15.6,d17.6,0p,f8.4)') 
c     +               NS,LM1,LM2,GLL0(LM1,LM2,NS),
c     +               zabs(GLL0(LM1,LM2,NS)),
c     +               DREAL(GLL0(LM1,LM2,NS)/zabs(GLL0(LM1,LM2,NS)))
c              END DO
c            END DO
c          END IF
c ------------------------------------------------------------------------
c        END IF

 70   CONTINUE                      ! NS = 1,NSHELL

      RETURN
 9000 format(2f18.10)

      END
