      double precision function dnrm2(r)
!-  'norm'-square
! ----------------------------------------------------------------------
!i Inputs:
!i   r     :vector
!o Outputs:
!o   dnrm2 :norm-square
!r Remarks:
! ----------------------------------------------------------------------
      implicit none
! Passed parameters:
      double precision r(3)

      dnrm2=r(1)*r(1)+r(2)*r(2)+r(3)*r(3)

      end
