C *********************************************************** 26.07.96 **
      LOGICAL FUNCTION TLAYPOS(NL,I,J)  
      implicit none
C ***********************************************************************
      INTEGER NL,I,J

      TLAYPOS= (I.EQ.J .OR. I.EQ.J+NL .OR. I.EQ.J-NL)

      RETURN                                                       
      END                                                         
