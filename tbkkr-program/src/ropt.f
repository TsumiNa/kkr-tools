C ***********************************************************************
      SUBROUTINE ROPT(LUN) 
      implicit none
C ***********************************************************************
C                                                                      
C     INPUT OF EXECUTION OPTIONS FROM LUN: 8(A8,2X)                  
C                                                                    
C ------------------------------------------------------------------------
C                                                                      
      CHARACTER*8   OPTC     
      integer LUN
C                                                                    
      COMMON/OPTC/  OPTC(8)                                           
      save  /optc/
C                                                                      
C                                                                      
      READ(LUN,1) OPTC                                               
    1 FORMAT(8(A8,2X))                                              
      WRITE(6,2) OPTC                                              
    2 FORMAT(79('-')/' EXECUTION OPTIONS:'/1X,A8,7('//',A8)/79('-'))
      RETURN                                                     
      END                                                       
