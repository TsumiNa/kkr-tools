
      SUBROUTINE ROTGLL(GMATLL,RCLSIMP,NATOMIMP,
     +           ATOMIMP,RATOM,NSHELL,NSH1,NSH2,
     +           DSYMLL,NSYMAT,ND,ISYMINDEX,IHANDLE,IGF)

      implicit none
c ************************************************************************
c
c     it calculates all the elements of the Green Function of
c     the impurity cluster using the GF calculated for the 
c     representative pairs.
c       _ _
c       n n'                      n n'
c       m m'               T      m m'       
c      G    (E) = SUM    D     * G    (E) * D 
c       L L'      L1 L2   L L1    L1L2      L2 L'
c
c                   _                 _
c                n    n            n'   n'
c     where   D R  = R     and  D R  = R 
c                m    m            m    m
c ------------------------------------------------------------------------

      include 'inc.fi'
      include 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)

      DOUBLE COMPLEX CZERO,CONE
      PARAMETER (CZERO= (0.0D0,0.0D0),CONE= (1.D0,0.D0))
      INTEGER LMAXSQ
      PARAMETER (LMAXSQ= (LMAXD+1)**2)
      INTEGER NSHELL,NSH1(*),NSH2(*),NSYMAT
      INTEGER NCLUSTD,NSIZE
      PARAMETER (NSIZE=NATOMIMPD*LMAXSQ,NCLUSTD=NSIZE*(NSIZE+1)/2)
      DOUBLE COMPLEX GMATLL(LMAXSQ,LMAXSQ,*),
c    +     GLL(LMAXSQ,LMAXSQ,NATOMIMPD,NATOMIMPD)
     +     GLL(LMAXSQ,LMAXSQ),
     +     DSYMLL(LMAXSQ,LMAXSQ,NSYMAXD),
     +     TPG(LMAXSQ,LMAXSQ)
c      COMPLEX*8 GCLUST(NCLUSTD)
	complex, allocatable ::  gclust2(:)
	integer iii,jjj,lmi,lmj
      DOUBLE PRECISION RATOM(3,*),RCLSIMP(3,*),SMALL,R1,
     +     ND(64,3,3),RATOMI(3)
      INTEGER I,II,J,JJ,K,LM1,LM2,NATOMIMP,ID,ISYM,
     +     ATOMIMP(NATOMIMPD),ISYMINDEX(*),IHANDLE,
     +     LIN,IGF
c
c  sym_file  file29    gf-indata  file 23
	real*8 rx(1000),ry(1000),rz(1000)
	integer numdiff,typ1(1000),typ2(1000),ix,iy,iz,
     +numat1(1000),numat2(1000)
	common/sym_file/rx,ry,rz,numdiff,typ1,typ2
     +,NUMAT1,NUMAT2
c
      EXTERNAL CINIT,ZCOPY

      DATA SMALL /  1.0D-10/

      
C      DO 10 I = 1,NATOMIMP
C      DO 10 J = 1,NATOMIMP
      
      allocate (gclust2(nclustd) )    
c      On Win OS whatever the os is x86 or x64, Windows limits static code and data to 2GB, The solution is to change the arrays from being declared with fixed bounds to being allocatable, and then using ALLOCATE to make them the desired size. 
 
      
      DO 10 IX=1,NUMDIFF
	  I=NUMAT1(IX)
	  J=numat2(IX)    
         DO II = 1,NSHELL 
            
            DO ID = 1,NSYMAT
               ISYM = ISYMINDEX(ID)

               
               DO K=1,3
                  RATOMI(K) = ND(ISYM,K,1)*RATOM(1,II) +
     +                 ND(ISYM,K,2)*RATOM(2,II) + 
     +                 ND(ISYM,K,3)*RATOM(3,II)
               ENDDO
               
               IF ( (ATOMIMP(I).EQ.NSH1(II) .AND.
     +              ATOMIMP(J).EQ.NSH2(II)) .OR.
     +              (ATOMIMP(I).EQ.NSH2(II) .AND. 
     +              ATOMIMP(J).EQ.NSH1(II)) ) THEN                  
                  
                  R1 = (RCLSIMP(1,J)-RCLSIMP(1,I)-RATOMI(1))**2  +
     +                 (RCLSIMP(2,J)-RCLSIMP(2,I)-RATOMI(2))**2  +
     +                 (RCLSIMP(3,J)-RCLSIMP(3,I)-RATOMI(3))**2
                  
                  IF (R1.LT.SMALL) THEN                     

c                     write (6,*) 'I',I,'J',J,'shell',II,'rot',ISYM
                     CALL ZGEMM('T','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,
     +                    DSYMLL(1,1,ID),LMAXSQ,GMATLL(1,1,II),
     +                    LMAXSQ,CZERO,TPG,LMAXSQ)   
                     
                     CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,
     +                    TPG,LMAXSQ,DSYMLL(1,1,ID),
c    +                    LMAXSQ,CZERO,GLL(1,1,I,J),LMAXSQ)
     +                    LMAXSQ,CZERO,GLL(1,1),LMAXSQ)
c
c
c
      if(igf.ne.0) then
	write(23) ((gll(iY,iZ),iY=1,lmaxsq),iZ=1,lmaxsq)
	endif


c
       
c	  	    
c
c                     write (6,*) '      i ','      j ',
c     +                            '      ii ','      id'
c                     write (6,*) i,j,ii,id
                     
                     GOTO 10
                     


                  ENDIF
                  
               ENDIF
                                 
            ENDDO
            
         ENDDO  
         
 10   CONTINUE
      go to 9999
 
       DO 11 I = 1,NATOMIMP
       DO 11 J = 1,NATOMIMP
            
         DO II = 1,NSHELL 
            
            DO ID = 1,NSYMAT
               ISYM = ISYMINDEX(ID)

               
               DO K=1,3
                  RATOMI(K) = ND(ISYM,K,1)*RATOM(1,II) +
     +                 ND(ISYM,K,2)*RATOM(2,II) + 
     +                 ND(ISYM,K,3)*RATOM(3,II)
               ENDDO
               
               IF ( (ATOMIMP(I).EQ.NSH1(II) .AND.
     +              ATOMIMP(J).EQ.NSH2(II)) .OR.
     +              (ATOMIMP(I).EQ.NSH2(II) .AND. 
     +              ATOMIMP(J).EQ.NSH1(II)) ) THEN                  
                  
                  R1 = (RCLSIMP(1,J)-RCLSIMP(1,I)-RATOMI(1))**2  +
     +                 (RCLSIMP(2,J)-RCLSIMP(2,I)-RATOMI(2))**2  +
     +                 (RCLSIMP(3,J)-RCLSIMP(3,I)-RATOMI(3))**2
                  
                  IF (R1.LT.SMALL) THEN                     

c                     write (6,*) 'I',I,'J',J,'shell',II,'rot',ISYM
                     CALL ZGEMM('T','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,
     +                    DSYMLL(1,1,ID),LMAXSQ,GMATLL(1,1,II),
     +                    LMAXSQ,CZERO,TPG,LMAXSQ)   
                     
                     CALL ZGEMM('N','N',LMAXSQ,LMAXSQ,LMAXSQ,CONE,
     +                    TPG,LMAXSQ,DSYMLL(1,1,ID),
c    +                    LMAXSQ,CZERO,GLL(1,1,I,J),LMAXSQ)
     +                    LMAXSQ,CZERO,GLL(1,1),LMAXSQ)
c
c
c
c      if(i.eq.1.and.j.eq.1) then
c	write(97,*) i,j,lmaxsq
c	write(97,1113) ((gll(lmi,lmj,i,j),lmi=1,lmj),lmj=1,lmaxsq)
c 1113 format(10f8.3)
c	endif
c	if(i.eq.2.and.j.eq.1) then
c	write(97,*) i,j,lmaxsq
c	write(97,1113) ((gll(lmi,lmj,i,j),lmi=1,lmaxsq),lmj=1,lmaxsq)
c	endif
c
        if(I.lt.J) then
	  do lmj=1,lmaxsq
	  do lmi=1,lmaxsq
	  jjj=(J-1)*lmaxsq+lmj-1
	  jjj=jjj*(jjj+1)/2
	  jjj=jjj+(i-1)*lmaxsq+lmi
	  gclust2(jjj)=gll(lmi,lmj)
	  enddo
	  enddo
	  endif
	  if(j.eq.i) then
	  do lmj=1,lmaxsq
	  do lmi=1,lmj
	  jjj=(J-1)*lmaxsq+lmj-1
	  jjj=jjj*(jjj+1)/2
	  jjj=jjj+(i-1)*lmaxsq+lmi
	  gclust2(jjj)=gll(lmi,lmj)
	  enddo
	  enddo
	  endif
c	  	    
c
c                     write (6,*) '      i ','      j ',
c     +                            '      ii ','      id'
c                     write (6,*) i,j,ii,id
                     
                     GOTO 11
                     


                  ENDIF
                  
               ENDIF
                                 
            ENDDO
            
         ENDDO  
         
 11   CONTINUE
                   
c      DO I=1,NATOMIMP
c         DO J=1,NATOMIMP
c            DO LM1 = 1 ,LMAXSQ
c               DO LM2 = 1 ,LMAXSQ
c                  WRITE (6,FMT='(4i4,2x,2d15.8)') I,J,
c     +                 LM1,LM2,GLL(LM1,LM2,I,J)
c               ENDDO
c            ENDDO
c         ENDDO
c      ENDDO

c      LIN = 0
c      DO  J=1,NATOMIMP
c         DO LM2=1,LMAXSQ
c            DO  I=1,J
c               IF (I.NE.J) THEN
c                  DO LM1=1,LMAXSQ
c                     LIN = LIN + 1
c                    GCLUST(LIN) = GLL(LM1,LM2,I,J)
c     for test purpouse
c                     write (6,*) (i-1)*lmaxsq +lm1,(j-1)*lmaxsq +lm2
c                  ENDDO
c               ELSE
c                  DO LM1=1,LM2
c                     LIN = LIN + 1
c                     GCLUST(LIN) = GLL(LM1,LM2,I,J)
c     for test purpouse
c                     write (6,*) (i-1)*lmaxsq +lm1,(j-1)*lmaxsq +lm2
c                  ENDDO
c               ENDIF
c
c            ENDDO
c      ENDDO



      IF (IGF.NE.0) THEN
c        CALL FXDRRL(IHANDLE,GCLUST,2*LIN)
choshino
      lin=jjj
      if(lin.ne.nclustd) then
	write(6,*) ' rotgll (lin.ne.nclustd)'
	stop
	endif
      write(91) lin,(GCLUST2(lm1),lm1=1,lin)
      write(6,*) 'ihandle lin nclustd=',ihandle,LIN,nclustd
ccc    test1
c	write(99,1111) LIN
c	WRITE(99,1112) (GCLUST(lm1),lm1=1,lin)
c	write(98,1111) jjj
c	write(98,1112) (Gclust2(lm1),lm1=1,jjj)
c1111 format(1x,'lin=',I12)
c 1112 format(10f8.3)
ccc   test1
c	stop
	
c         do lm1=1,lmaxsq
c            do lm2=1,lmaxsq
c               write (6,*) lm1,lm2
c               write (6,*) GLL(lm1,lm2,1,1)
c       	write (6,*) GLL(lm1,lm2,18,18)
c	      enddo
c          enddo
c       	write (6,*) 'GLL(25,25,1,1)',GLL(25,25,1,1)
      ENDIF
      

c      stop
 9999 continue
      deallocate(gclust2)
      RETURN

      END                       ! SUBROUTINE ROTGLL
