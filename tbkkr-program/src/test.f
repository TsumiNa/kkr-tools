



C *********************************************************** 17.05.91 **
      LOGICAL FUNCTION TEST(STRING) 
      implicit none
C ***********************************************************************
C                                                                      
C     TEST = 'STRING  ' IS CONTAINED IN /TESTC/.                      
C                                                                    
C ------------------------------------------------------------------------
C                                                                      
      COMMON/TESTC/ TESTC(16)                                         
      save  /testc/
C                                                                    
      integer i
      character*8      STRING   ,TESTC                              
C                                                                       
      TEST=.FALSE.                                                    
      DO 1 I=1,16                                                    
        IF(STRING.EQ.TESTC(I)) TEST=.TRUE.                         
 1    END DO
      RETURN                                                       
      END                                                         
