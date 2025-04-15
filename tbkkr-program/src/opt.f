C ***********************************************************************
      LOGICAL FUNCTION OPT(STRING)
      implicit none
C ***********************************************************************
C                                                                      
C     OPT = 'STRING  ' IS CONTAINED IN /OPTC/.                        
C                                                                       
C ------------------------------------------------------------------------
C                                                                      
      COMMON/OPTC/  OPTC(8)                                           
      save  /optc/
C                                                                    
      character*8      STRING   ,OPTC      
      integer I
C                                                                  
C                                                                      
      OPT=.FALSE.                                                     
      DO 1 I=1,8                                                     
        IF(STRING.EQ.OPTC(I)) OPT=.TRUE.
 1    END DO
      RETURN                                                       
      END                                                         
