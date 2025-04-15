      SUBROUTINE CMOMSREAD(NLBASIS,NRBASIS,NAEZ,CMOMHOST,VACFLAG,KAOEZ)
      implicit none
c *****************************************************
c * Reads the CMOMHOST from the decimation file 
c * This sub must be called after the t-matrices for
c * all the energies are read in.
c * It returns the cmomhost array. First the left bulk
c * then the right bulk are indexed.
c * Left is in unit 37 right in unit 38 
c *                                        29.10.99
c *****************************************************
      include 'inc.fi'
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NLBASIS,NRBASIS
      INTEGER KAOEZ(NAEZD+NEMBD)
      DOUBLE PRECISION CMOMHOST(LMPOTD,NEMBD)
      LOGICAL VACFLAG(2)
c     ..Local variables ..
      INTEGER IH,IH1,LMPOTL,LM,LMPOTR,NAEZL,NAEZR,IHL,IHR,NAEZ
      DOUBLE PRECISION C00(LMPOTD)

 
c     I read this only once and store , since I read all the 
c t-matrices the file is in the correct place.
c DO NOT forget to rewind it before the iteration
c   
      !CMOMHOST = 0.d0
      IF (VACFLAG(1)) THEN
c     You have vacuum on this side
         DO IH=1,NLBASIS
            DO LM=1,LMPOTD
               CMOMHOST(LM, IH) = 0.d0   
            END DO 
         END DO 
      ELSE                      ! in case of no vacuum read from file
         READ(37,1180) NAEZL,LMPOTL
         IF (NAEZL.NE.NLBASIS) THEN
            WRITE(6,*) ' DECIMATION ERROR '
            WRITE(6,*) ' NLBASIS=',NLBASIS,' BUT ',NAEZL,
     &           ' ATOMS FOUND IN FILE '
            STOP 'READ DECIMATION CMOMS'
         END IF
         DO IH=1,NAEZL
            READ(37,*) IHL
            IF (IHL.NE.IH) STOP 'DECIFILE ERROR IN CMOMS'
            READ(37,1090) (C00(lm),LM=1,LMPOTL)

            IH1 = KAOEZ(NAEZ+IH)  ! get the correct atom
            write(6,*) 'Host left ',ih1
            DO LM=1,LMPOTL
               CMOMHOST(LM, IH1) = C00(LM)   
            END DO 
         END DO
      END IF     
c     *************** Now the Right Host **************
      IF (VACFLAG(2)) THEN
c     You have vacuum on this side
         DO IH=1,NRBASIS
            DO LM=1,LMPOTL
               CMOMHOST(LM, NLBASIS+IH) = 0.d0   
            END DO 
         END DO 
      ELSE                      ! in case of no vacuum read from file
         READ(38,1180) NAEZR,LMPOTR
         IF (NAEZR.NE.NRBASIS) THEN
            WRITE(6,*) ' DECIMATION ERROR '
            WRITE(6,*) ' NRBASIS=',NRBASIS,' BUT ',NAEZR,
     &           ' ATOMS FOUND IN FILE '
            STOP 'READ DECIMATION CMOMS'
         END IF
         DO IH=1,NAEZR
            READ(38,*) IHR
            IF (IHR.NE.IH) STOP 'DECIFILE ERROR IN CMOMS'
            READ(38,1090) (C00(lm),LM=1,LMPOTR)

            IH1 = KAOEZ(NAEZ+NLBASIS+IH)  ! get the correct atom
            WRITE(6,*) 'filling up right host',naezl+ih1
            DO LM=1,LMPOTR
               CMOMHOST(LM,NLBASIS+IH1) = C00(LM)   
            END DO 
         END DO    
      END IF                    ! both cases ok for right host
c     Now map to the CMOMHOST array the ordering is important
c     Consider this later...
 1090 FORMAT(4D22.14)    
 1180 FORMAT(5X,2I6)
      RETURN
      END
