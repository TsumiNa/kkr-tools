       SUBROUTINE QDOSKNET(NKPOI,NKLINES,QDOSKP,NQDOSKP,IE,QDOSKDIM)
c **********************************************
c * This subroutine produces the k-point mesh
c * along lines in the IBZ This is done for the 
c * q-dos calculation       
c *                                 19.4.2000
c **********************************************
       implicit none
       INTEGER NLDIM
       PARAMETER (NLDIM=5)
       INTEGER NKPOI,QDOSKDIM,NQDOSKP,IE
       CHARACTER*80 UIO
       INTEGER IER,NKLINES,KSUM,I,IK,I1,I2,I3,I4,I5,I6,J,K
       INTEGER NKMESH(NLDIM)
       DOUBLE PRECISION KPOINTS(6,NLDIM),BZKP0(6),DK(6)
       DOUBLE PRECISION QDOSKP(3,*)
c     
c     reading from the inputcard the points defining the symmetry
c     directions for a band structure calculation.
c     For each symmetry direction the points are 2.
       
       CALL IoInput('BNDKPNTS  ',UIO,1,7,IER)
       READ (UNIT=UIO,FMT=*)  NKLINES
       if (nklines.lt.1) then 
            write(6,*) 'Option QDOS without NKLINES > 1 '
            write(6,*) 'Cannot create lines for QDOS '
            write(6,*) 'Program is STOPING '
            STOP
         end if
       IF (NKLINES.GT.0) THEN
          WRITE(6,*) ' Q-DOS ALONG LINES K-NET IS MADE NEW !!! '
          IF (NKLINES.GT.NLDIM) THEN
          WRITE(6,*) 'ERROR STOP in sub QDOSKNET '
          WRITE(6,*) 'Increase papameter NLDIM to : ',NKLINES
          STOP
          END IF
       ELSE
          WRITE(6,*) ' Q-DOS IN FULL IRR. B.Z. !!! '
          KSUM = NKPOI
          NQDOSKP = NKPOI
          IF (KSUM.GT.QDOSKDIM) THEN
          write(6,*) 'Error stop QDOSKNET'
          write(6,*) 'Increase parameter QDOSKDIM in main.f to:',KSUM
          STOP
          END IF
          RETURN
       END IF
       
       DO I=1,NKLINES
          CALL IoInput('KPNTLIST  ',UIO,I,7,IER)
          IF (I.NE.1) THEN
             READ (UNIT=UIO,FMT=*)  (KPOINTS(J,I),J=1,6),
     +            NKMESH(I)
          ELSE
             READ (UNIT=UIO,FMT=*)  (KPOINTS(J,I),J=1,6)
          ENDIF
c     write (6,*) (KPOINTS(J,I),J=1,6)
       ENDDO
       KSUM=0
       DO I=1,NKLINES-1
          DO J=1,3
             DK(J) = KPOINTS(J,I+1) - KPOINTS(J,I)
             IF (NKMESH(I+1).GT.1) DK(J) = DK(J)/(NKMESH(I+1)-1)
          ENDDO
          DO K=1,NKMESH(I+1)
             KSUM = KSUM + 1
             IF (KSUM.GT.QDOSKDIM) THEN
             WRITE(6,*) ' *** Increase Dimensions for q-dos ploting ***'
             WRITE(6,*) ' Change parameter QDOSKDIM in main program '
             WRITE(6,*) ' To a value greater than : ',NQDOSKP
             WRITE(6,*) '   PROGRAM IS STOPING '
             STOP 
             END IF
             DO J=1,3
                QDOSKP(J,KSUM) = KPOINTS(J,I) +  (K-1)*DK(J)
             ENDDO
          ENDDO
          
       ENDDO              
       IF (IE.EQ.1) THEN
          write (6,*) 'K POINTS FOR BAND CALC. ',KSUM
       ENDIF
       NQDOSKP = KSUM
       
       RETURN
       END        
