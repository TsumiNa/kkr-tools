      SUBROUTINE SAVEGH(GHOST,RFCTOR,NDIM,IFILE,NLAST,NSPIN,NLST,IEN,
     +                  IENI,SAVKEY,NHSPIN,ISW80,NREP,GTEMP)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     the green's are symmetrisized in this subroutine if SAVKEY is
c     false .
c     to save storage and speed up the program  the sym. green's
c     function are stored on a temporay disc .
c
c
c             e = energy
c             ek = sqrt(energy)
c             df = weight of energy-intervall
c
c   THIS IS A CHANGED VERSION TO DEAL WITH SEMICORE STATES ALSO
C                                                nikos 7/5/96
c-----------------------------------------------------------------------
C     .. Parameters ..
      INTEGER ICGD,NLSTD
      PARAMETER (ICGD=6000,NLSTD=15686)
      INTEGER IEMXD,NSEC
      PARAMETER (iemxd=150,nsec=689)
C     ..
C     .. Scalar Arguments ..
      REAL*8 RFCTOR
      INTEGER IFILE,ISW80,NHSPIN,NLAST,NREP,NSPIN
      LOGICAL SAVKEY
C     ..
C     .. Array Arguments ..
      COMPLEX*16 GHOST(NSEC,NSEC)
      REAL*8 IEN(IEMXD,*),IENI(IEMXD,*)
      COMPLEX GTEMP(NSEC,NSEC)
      INTEGER NDIM(*),NLST(*)
C     ..
C     .. Local Scalars ..
      COMPLEX*16 CZERO,DF,E,EFMTZ,EK
      INTEGER I,I1,I2,ICG,IE,IELAST,IG,ISPIN,J,NDIMNP,NQ,NR
      integer nsemcore,iecore
      double complex ecore(iemxd),dfcore(iemxd)
      double precision e1,e2
C     ..
C     .. Local Arrays ..
      COMPLEX*16 GZ(NLSTD)
      REAL*8 GCOEFF(ICGD)
      INTEGER I1G(ICGD)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DIMAG,REAL
C     ..
C     .. Save statement ..
      SAVE CZERO
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
C     ..

      CALL CINIT(NLSTD,GZ)
c
      REWIND IFILE
      NHSPIN = 1
      nsemcore = 0
      DO 120 ISPIN = 1,NSPIN
c
        IF (ISPIN.NE.2 .OR. NHSPIN.NE.1) THEN

          READ (IFILE) IELAST,NHSPIN,EFMTZ
          EFMTZ = EFMTZ/ (RFCTOR*RFCTOR)
          WRITE (6,FMT='(a14,1p,2d10.3)') ' fermi-energy ',EFMTZ
             if (nsemcore.eq.1) then
                read(5,1001) iecore,e1,e2
                ielast = ielast+iecore  
1001    format (i5,4d20.12)
             do i = 1,iecore 
              read(5,1002) ecore(i),dfcore(i)
             end do
1002    format (4d20.12)
             end if
c
          NLST(ISPIN) = IELAST
          IF (NSPIN.EQ.2 .AND. NHSPIN.EQ.1) NLST(NSPIN) = IELAST 
c
          IF (IELAST.GT.IEMXD) THEN
            STOP 19

          ELSE

c
            DO 110 IE = 1,IELAST
c
c---> read  struc. host green's functions
c
              if (ie.le.iecore) then 
                e = ecore(ie)
                ek = sqrt(e)
                df = dfcore(ie)
                do i=1,nlast
                gz(i) = czero
                end do
              else
              READ (IFILE) E,EK,DF, (GZ(I),I=1,NLAST)
              end if
              IF (IE.EQ.IELAST) WRITE(6,*) 'LAST ENERGY = ',
     &                                      E/(RFCTOR*RFCTOR)                       
c
c******** energy is converted from pu. to ry   ************************
c
              E = E/ (RFCTOR*RFCTOR)
              EK = EK/RFCTOR
c                 if(kvrel.eq.1)ek=cdsqrt(e+(e/c)**2)
              DF = DF/RFCTOR/RFCTOR
              DO 10 I = 1,NLAST
                GZ(I) = GZ(I)/RFCTOR
   10         CONTINUE
c
c**********************************************************************
c
              write(6,*) 'Energy (Ry) :',e
              IF (.NOT.SAVKEY) THEN
                DO 100 NR = 1,NREP
                  NDIMNP = NDIM(NR)
                  DO 30 I2 = 1,NDIMNP
                    DO 20 I1 = 1,NDIMNP
                      GHOST(I1,I2) = CZERO
   20               CONTINUE
   30             CONTINUE
c
c---> read symmetrisation coeffients
c
                  REWIND 20 + NR
                  READ (20+NR) NQ,NLAST
                  DO 50 IG = 1,NLAST
                    READ (20+NR) ICG, (I1G(J),J=1,ICG),
     +                (GCOEFF(J),J=1,ICG)
                    IF (ICG.GT.ICGD) THEN
                      WRITE (6,FMT=*)
     +                  'ICG > SPECIFIED DIMENSIONS: STOP IN SAVEGH'
                      STOP 'SAVEGH'

                    END IF

                    DO 40 J = 1,ICG
                      I1 = I1G(J)/10000
                      I2 = I1G(J) - I1*10000
                      GHOST(I1,I2) = GHOST(I1,I2) + GCOEFF(J)*GZ(IG)
   40               CONTINUE

   50             CONTINUE
                  DO 70 I1 = 1,NDIMNP
                    DO 60 I2 = I1,NDIMNP
                      GHOST(I2,I1) = GHOST(I1,I2)
   60               CONTINUE
   70             CONTINUE
c
c---> store struc. host green's functions on disc
c
                  IF (ISW80.EQ.0) THEN
                    DO 90 I1 = 1,NDIMNP
                      DO 80 I2 = 1,I1
                        GTEMP(I2,I1) = GHOST(I2,I1)
   80                 CONTINUE
   90               CONTINUE
                    WRITE (80) E,EK,DF, ((GTEMP(I2,I1),I2=1,I1),I1=1,
     +                NDIMNP)
                  END IF
  100           CONTINUE

              END IF

              IEN(IE,ISPIN) = REAL(E)
              IENI(IE,ISPIN) = DIMAG(E)
c
              IF (NSPIN.EQ.2 .AND. NHSPIN.EQ.1) THEN
                IEN(IE,NSPIN) = REAL(E)
                IENI(IE,NSPIN) = DIMAG(E)
              END IF

  110       CONTINUE
          END IF

        END IF

  120 CONTINUE
      RETURN

      END
