c ************************************************************************
      SUBROUTINE RIMAD (IFILE,NAEZ,RBASIS,IMAD,JMAD,NVMAD)
      implicit none
c ************************************************************************
c
      include 'inc.fi'
c
c     .. arguments
c
      DOUBLE PRECISION RBASIS(3,*)  ! pos. of basis atoms in EZ
c
      INTEGER
     +     IFILE,
     +     NAEZ,                    ! number of atoms in EZ
     +     NVMAD                    ! number of madelung shells
c
      INTEGER
     +     IMAD(NAEZD,*),           ! madelung constant for 
                                    ! position pair(i,j)
     +     JMAD(NVMADD,*)           ! position indices for madelung constant 
c
c     .. locals
c
      INTEGER I,J,JJ,JJJ,N
      DOUBLE PRECISION DR,R(3)
      INTRINSIC SQRT
c ------------------------------------------------------------------------
      write(6,*) '>>> RIMAD : read madelung shells from file'
      READ (IFILE,*) N,NVMAD
      write(6,*) 'N,NVMAD ',N,NVMAD
      IF (N.NE.NAEZ) THEN
        write(6,9000) N,NAEZ
        GOTO 900
      END IF
      IF (NVMAD.GT.NVMADD) THEN
        write(6,9010) NVMAD
        STOP '1 - RIMAD'
      END IF
      DO 10 J=1,NAEZ
        READ (IFILE,*) (R(I),I=1,3)
        DR=SQRT((R(1)-RBASIS(1,J))**2+
     +          (R(2)-RBASIS(2,J))**2+
     +          (R(3)-RBASIS(3,J))**2 )
        IF (DR.GT.1.D-5) THEN
          write(6,9020) J
          GOTO 900
        END IF
 10   END DO

      READ (IFILE,*) ((IMAD(I,J),J=1,NAEZ),I=1,NAEZ)
      READ (IFILE,*) (JJJ,JMAD(N,1),(JMAD(N,J+1),J=1,2*JMAD(N,1)),
     +                N=1,NVMAD)

      write(6,*) '<<< RIMAD : o.k.'
      CLOSE(IFILE)
      IFILE = 0
      RETURN

 900  write(6,*) '<<< RIMAD : ERROR OCCUR. <<< <<< <<< <<< <<<'
      CLOSE(IFILE)
      RETURN

 9000 FORMAT('Number of atoms in unit cell not equal (',
     +       I5,',',I5,').')
 9010 FORMAT('Increase parameter NVMADD in inc.fi to ',I6,'.')
 9020 FORMAT('Error in basis vectors at position',I5,'.')

      END
