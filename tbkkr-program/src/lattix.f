C ************************************************************************
      SUBROUTINE LATTIX(LATT,AA,BBYA,CBYA,LTSP,GAMMA,
     +     VOLUC,BV1,RBV1,
     +     RR,NR)
      implicit none
C ************************************************************************
C LATTIX GENERATES THE REAL SPACE AND RECIPROCAL LATTICES.
C LATT DETERMINES THE BRAVAIS LATTICE, DEFINITONS SEE BELOW.
C BV(I,J) ARE INTERNAL BASIS VECTORS, WITH I=X,Y,Z AND J=A,B,C.
C BV1 ARE THE TRUE(!) BASIS VECTORS, IN UNITS OF A0.
C RECIPROCAL SPACE VECTORS ARE IN UNITS OF 2*PI/A0.
C HERE AA IS THE LATTICE CONSTANT AND BBYA, CBYA ARE IN UNITS OF AA,
C RR ARE THE DIRECT SPACE VECTORS.
C NR+1 IS THE NUMBER OF DIRECT SPACE VECTORS CREATED
C (STRUCTURE DEPENDENT OUTPUT). 
C VOL/8 IS THE VOLUME OF THE UNIT CELL, IN UNITS OF (A0)**3.
C ************************************************************************
C
      INCLUDE 'inc.fi'
C
      INTEGER
     +     NR,                       ! number of real space vectors
     +     LTSP,
     &    LATT         !      LATT = 0  : SIMPLE CUBIC
                            !      LATT = 3  : SIMPLE TETRAGONAL
                            !      LATT = 8  : SIMPLE ORTHORHOMBIC
                            !      LATT = 10 : SIMPLE ORTHORHOMBIC
                            !     see below for detial 

      
      
      integer I   , J        ! for do statement
C
      DOUBLE PRECISION
     +     VOL,VOLUC,
     +     AA,BBYA,CBYA,
     +     LATSP,GAMMA,
     +     A,B,
     +     FACV,PI,TPI
C
      DOUBLE PRECISION
     +     BV(3,3),                 ! INTERNAL BASIS
     +     BV1(3,3),                ! TRUE BASIS VECTORS, IN UNITS AA
     +     RBV(3,3),                ! INTERN RECIPROCAL BASIS 
     +     RBV1(3,3),               ! RECIPROCAL BASIS IN 2*PI/A
     +     RR(3,0:NRD),
     +     ZERO_VEC(3)
C
c
      EXTERNAL CROSPR,SPATPR,VADD,VSUB,VEQ
c
      PARAMETER (PI   = 3.14159265358979312D0)
      PARAMETER (TPI  = 2.0d0 * PI)
c
      double precision SRTW, SRTR
      DATA ZERO_VEC / 0.0D0,0.0D0,0.0d0 /
      DATA SRTW/1.4142135623730950488D0/
      DATA SRTR/1.7320508075688772935D0/
c
c ------------------------------------------------------------------------
      IF (LATT.GT.11) THEN
           WRITE(6,6010) LATT
 6010      FORMAT(' LATTIX: LATT=',I3,' NOT IMPLEMENTED.STOP')
           STOP
      END IF
C
      DO 110 I = 1,3
        DO 111 J = 1,3
          BV(J,I) = 0.D0
 111    END DO
 110  END DO
C
C VOLUC CONTAINS THE UNIT CELL VOLUME IN UNITS OF A0**3.
C
      VOLUC=0.D0
C
C  !! ATTENTION : THERE MAY BE SEVERAL POSSIBILITIES FOR A CERTAIN
C                 VALUE OF LATT.   ONE OF THEM IS CHOSEN ACCORDING
C                 TO THE FOLLOWING GOTO STATEMENT.
C                 FOR LATT=0 THE FOLLOWING GOTO IS NOT EFFECTIVE!
c      GOTO (1110,1210,1010,1210,1510,1610,1710,1010,1910,1920)  ! ori
C     LATT =  1    2    3    4    5    6    7    8    9    10
      GOTO (1110,1210,1010,1210,1510,1510,1710,1010,1910,1010,
c            11
     +     1010 )  
     +     LATT
C
C     LATT = 0  : SIMPLE CUBIC
C     LATT = 3  : SIMPLE TETRAGONAL
C     LATT = 8  : SIMPLE ORTHORHOMBIC
C     LATT = 10 : SIMPLE ORTHORHOMBIC
C                (DEPENDING ON THE ABOVE GOTO)
C
      write(6,*) '>>> LATTIX: simple cubic lattice is created'
 1010 continue 
      if(latt.eq.3) write(6,*) 
     +     '>>> LATTIX: simple tetragonal lattice is created'
      if(latt.eq.8 .or. latt.eq.10) write(6,*) 
     +     '>>> LATTIX: simple orthorhombic lattice is created'
      if(latt.eq.11) write(6,*) 
     +     '>>> LATTIX: lattice for WIRE calculation is created'
      BV(1,1)=2.D0
      BV(2,2)=2.D0*BBYA
      BV(3,3)=2.D0*CBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      GOTO 210
C
C     LATT = 1 : FACE CENTERED CUBIC
C                INCLUDING TRIGONAL DISTORTION, IF LTSP = 1
C
 1110 write(6,*) '>>> lattix: face centered cubic lattice is created'
      BV(2,1)=BBYA
      BV(3,1)=CBYA
      BV(1,2)=1.D0
      BV(3,2)=CBYA
      BV(1,3)=1.D0
      BV(2,3)=BBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      IF (LATT.NE.1.OR.LTSP.NE.1) GOTO 210
C
      WRITE(6,6020) GAMMA
 6020 FORMAT(/' TRIGONAL DISTORTION GAMMA = ',F13.6/,
     &' CHECK UNITS OF GAMMA',/)
      A=DEXP(-GAMMA*0.5D0)*SRTW/SRTR
      B=DEXP(GAMMA)*SRTW/SRTR
      BV(1,1)=A
      BV(2,1)=0.D0
      BV(3,1)=SRTW*B
      BV(1,2)=-A*0.5D0
      BV(2,2)=SRTR*0.5D0*A
      BV(3,2)=SRTW*B
      BV(1,3)=-A*0.5D0
      BV(2,3)=-SRTR*0.5D0*A
      BV(3,3)=SRTW*B
      GOTO 210
C
C     LATT = 2 : BODY CENTERED CUBIC
C                INCLUDING TRIGONAL DISTORTION, IF LTSP = 1
C     LATT = 4 : BODY CENTERED TETRAGONAL
C
 1210 if(latt.eq.2)
     * write(6,*) '>>> lattix: body centered cubic lattice is created'
      if(latt.eq.4)
     *     write(6,*) '>>> lattix: body centered tetragonal',
     *     'lattice is created'
      BV(1,1)=-1.D0
      BV(2,1)=BBYA
      BV(3,1)=CBYA
      BV(1,2)=1.D0
      BV(2,2)=-BBYA
      BV(3,2)=CBYA
      BV(1,3)=1.D0
      BV(2,3)=BBYA
      BV(3,3)=-CBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      IF (LATT.NE.2.OR.LTSP.NE.1) GOTO 210
C
      WRITE(6,6020) GAMMA
      A=DEXP(-GAMMA*0.5D0)*SRTW/SRTR*2.D0
      B=DEXP(GAMMA)*SRTW/SRTR*0.25D0
      BV(1,1)=A
      BV(2,1)=0.D0
      BV(3,1)=SRTW*B
      BV(1,2)=-A*0.5D0
      BV(2,2)=SRTR*0.5D0*A
      BV(3,2)=SRTW*B
      BV(1,3)=-A*0.5D0
      BV(2,3)=-SRTR*0.5D0*A
      BV(3,3)=SRTW*B
      GOTO 210
C
C     LATT = 5 : HEXAGONAL
C
 1510 write(6,*) '>>> lattix: hexagonal lattice is created'
      BV(1,1)=SRTR
      BV(2,1)=-1.D0
      BV(1,2)=SRTR
      BV(2,2)=1.D0
      BV(3,3)=2.D0*CBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      GOTO 210
C
C     LATT = 6 : SPECIAL TRIGONAL LATTICE FOR AF2 FCC
C
c 1610 write(6,*) '>>> lattix: special trigonal lattice for',
c     *           ' af2 fcc is created'
c      BV(1,1)=2.D0
c      BV(1,2)=1.D0
c      BV(2,2)=SRTR
c      BV(1,3)=2.D0*CBYA
c      BV(2,3)=2.D0*CBYA/SRTR
c      BV(3,3)=4.D0*CBYA*SRTW/SRTR
c      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
c      GOTO 210
C
C     LATT = 7 : MONOCLINIC (SPECIAL AXES)
C
 1710 write(6,*) '>>> lattix: monoclinic (special axes)',
     +     ' lattice is created'
      BV(1,1)=2.D0
      BV(2,1)=2.D0*BBYA
      BV(1,2)=-1.D0
      BV(3,2)=CBYA
      BV(2,3)=-BBYA
      BV(3,3)=CBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      GOTO 210
C
C     LATT = x : SIMPLE CUBIC (SPECIAL AXES)
C                (DEPENDING ON THE ABOVE GOTO)
C NOT FINISHED!
c1810 BV(1,1)=2.D0/SRTR
c     BV(2,1)=-SRTW
c     BV(3,1)=SRTW/SRTR
c     BV(1,2)=2.D0/SRTR
c     BV(2,2)=SRTW
c     BV(3,2)=SRTW/SRTR
c     BV(1,3)=-2.D0/SRTR
c     BV(3,3)=2.D0*SRTW/SRTR
c     GOTO 210
C
C     LATT = 9 : BASE CENTERED ORTHORHOMBIC
C                (DEPENDING ON THE ABOVE GOTO)
C
 1910 write(6,*) '>>> lattix: base centered orthorhombic',
     +     ' lattice is created'
      BV(1,1)=1.D0
      BV(2,1)=-BBYA
      BV(1,2)=1.D0
      BV(2,2)=+BBYA
      BV(3,3)=2.D0*CBYA
      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
      GOTO 210
c ------------------------------------------------------------------------
C
C     LATT = 10 : special case
C                (DEPENDING ON THE ABOVE GOTO)
C
c 1920 write(6,*) '>>> lattix: special',
c     +     ' lattice is created'
c      BV(1,1)=2.D0
c      BV(2,2)=2.D0*BBYA
c      BV(3,3)=2.D0*CBYA
c      CALL SPATPR(BV(1,2),BV(1,3),BV(1,1),VOL)
c      GOTO 210
c ------------------------------------------------------------------------
c
c ---> start computing bravais lattices;
c      make correction for atomic units!
c
  210 CONTINUE
      VOLUC=0.125D0*DABS(VOL)*AA**(3.D0)
      write(6,1001) vol,voluc
c
c ---> true basis vectors bv1:
c
      write(6,*) 'true basis vectors'
      A = 0.5D0*AA
      DO 225 J = 1,3
         DO 220 I = 1,3
            BV1(I,J) = BV(I,J)*A
 220     CONTINUE
         WRITE(6,1000) J, BV1(1,J),BV1(2,J),BV1(3,J)
 225  CONTINUE
c
c ---> construct the reciprocal basis vectors:
c
      CALL CROSPR(BV1(1,2),BV1(1,3),RBV(1,1))
      CALL CROSPR(BV1(1,3),BV1(1,1),RBV(1,2))
      CALL CROSPR(BV1(1,1),BV1(1,2),RBV(1,3))
c
c ---> multiply with appropriate prefactors:
c
      write(6,*) 'reciprocal basis vectors (in units of 2*pi/alatc)'
      FACV=TPI/VOLUC
      DO 240 I=1, 3
         DO 230 J=1, 3
            RBV(J,I)=RBV(J,I)*FACV
            RBV1(J,I)=RBV(J,I)*AA/TPI
 230     CONTINUE
c         write(6,1000) i, rbv(1,i),rbv(2,i),rbv(3,i)
         write(6,1000) i, rbv1(1,i),rbv1(2,i),rbv1(3,i)
  240 CONTINUE
c
c ---> now generate the real space lattice vectors:
c
      CALL RRGEN(BV1,AA,LATT,RR,NR)
c
c ---> test on volume unit cell:
c
      IF (VOLUC.LT.1.0D-5) THEN
         WRITE(6,909) LATT, VOLUC
  909    FORMAT(//' STOP. VOLUC IN LATTIX IS ZERO FOR LATTICE TYPE:',
     +   I4,' VOLUC= ',D13.4)
         STOP
      END IF
C
 1000 FORMAT(i4,3F15.10)
 1001 FORMAT( ' vol,voluc: ',2F10.4)
      RETURN
      END
