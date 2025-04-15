C ************************************************************************
      SUBROUTINE CROSPR(X,Y,Z)
      implicit none
C ************************************************************************
C     CROSP COMPUTES THE CROSS PRODUCT OF X AND Y RETURNING
C     IT INTO Z.
C ------------------------------------------------------------------------
      DOUBLE PRECISION X(*), Y(*), Z(*)
      Z(1)=X(2)*Y(3)-X(3)*Y(2)
      Z(2)=X(3)*Y(1)-X(1)*Y(3)
      Z(3)=X(1)*Y(2)-X(2)*Y(1)
      RETURN
      END
