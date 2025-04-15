       LOGICAL FUNCTION CLUSTCOMP(RCLS,IC1,N1,RCLS1,N2)
      implicit none
c  This function returns true if cluster number ic1
c  is equal to cluster ic2
c  RCLS        clusters coordinates
c  IC1         First cluster
c  N1          Number of atoms in IC1 cluster
c  IC2         Second cluster
c  
c        include 'inc.fi'
        include 'inc.cls'
        DOUBLE PRECISION RCLS(3,NACLSD,*),RCLS1(3,NACLSD)
        INTEGER IC1,N1,IC2,N2
        integer  N,I
        DOUBLE PRECISION R
        CLUSTCOMP = .FALSE.
        IF (N1.EQ.N2) THEN
           R = 0
           DO N=1,N1
              DO I=1,3
              R = R + ( RCLS(I,N,IC1) - RCLS1(I,N))**2 
              END DO
           END DO
        IF (ABS(R).LT.1.D-4) CLUSTCOMP = .TRUE.
        END IF
        END 
