       SUBROUTINE CHGLATT(ALAT,RM,RM1,DRM,RELAX,NATPER,NSHELL,NATOM,
     +                    IOPER,NIMP,ND)
       implicit none
C  ****************************************************************
C  *  In case of lattice relaxation this subroutine changes the
C  *  lattice to the new position
C  *  on input RM is the ideal host lattice 
C  *  on output RM is the new shifted lattice and RM1 is the old
C  *  ideal host lattice DRM is the shift.
C  ****************************************************************
C
C  ...Parameter statements
       INTEGER NATYPD,NTREFD,NATOMD,NTPERD,NSPIND
       PARAMETER (NATYPD=38,NTREFD=1,NATOMD=102,NTPERD=NATYPD-NTREFD)
       INTEGER NATPER,NATOM
       INTEGER IOPER(NATOMD),NSHELL(NTPERD),NIMP(NATOMD),ND(48,3,3)
       INTEGER N1,I,N,K,I1
       REAL*8 RM1(3,NATOMD),RM(3,NATOMD),DX(3),
     +                  DRM(3,NTPERD),RELAX(NTPERD)
       REAL*8 ALAT
c
c
       N1=0
       DO 10 N=1,NATPER
          DO 20 I=1,3
             DRM(I,N) = RM(I,N1+1)*RELAX(N)/100.d0
 20       CONTINUE
          N1=NSHELL(N)+N1
 10    CONTINUE
       WRITE(6,9550)
       DO 30 I=1,NATPER
          WRITE(6,9560) I,RELAX(I)
 30    CONTINUE
C
C Using the representative atom shift find the position of the
c  shifted lattice. Notice that the inverse rotation is done.
C
       DO 40 N = 1,NATOM
          K = IOPER(N)
          N1 = NIMP(N)
          DX(1)=ND(K,1,1)*DRM(1,N1) + ND(K,2,1)*DRM(2,N1) +
     +          ND(K,3,1)*DRM(3,N1)
          DX(2)=ND(K,1,2)*DRM(1,N1) + ND(K,2,2)*DRM(2,N1) +
     +          ND(K,3,2)*DRM(3,N1)
          DX(3)=ND(K,1,3)*DRM(1,N1) + ND(K,2,3)*DRM(2,N1) +
     +          ND(K,3,3)*DRM(3,N1)
          DO 50 I=1,3
             RM1(I,N) = RM(I,N)
             RM(I,N) = RM(I,N) + DX(I)
 50       CONTINUE
          WRITE(6,FMT='(3F10.6)') ( RM(I1,N),I1=1,3)
 40    CONTINUE
       DO 60 N=1,NATPER
         DO 70 I=1,3
           DRM(I,N) = DRM(I,N)*ALAT
 70      CONTINUE
 60    CONTINUE
c
9550   FORMAT ('**** THE LATTICE IS RELAXED **** ',//)
9560   FORMAT ('Shell ',i3,' is relaxed by ',f5.2,' %'/)
       END
