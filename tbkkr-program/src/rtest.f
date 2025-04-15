C ************************************************************************
      SUBROUTINE RTEST(LUN) 
      implicit none
C ************************************************************************
C                                                                      
C     INPUT OF TEST OPTIONS FROM LUN: 8(A8,2X)/8(A8,2X)             
C                                                                    
C ------------------------------------------------------------------------
C                                                                      
      CHARACTER*8   TESTC   
      integer LUN
C                                                                    
      COMMON/TESTC/ TESTC(16)                                     
      save  /testc/
C                                                                  
C                                                                      
      READ(LUN,1) TESTC                                              
    1 FORMAT(8(A8,2X)/8(A8,2X))                                    
      WRITE(6,2) TESTC                                                 
    2 FORMAT(79('-')/' TEST OPTIONS:'/2(1X,A8,7('//',A8)/)/79('-'))
      RETURN                                                         
      END                                                           
