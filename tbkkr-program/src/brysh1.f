c ************************************************************************
      SUBROUTINE BRYSH1(Y,X,INS,IRMIN,IRC,NATPS,NATYP,NSPIN,IMAP,LMPOT)
      implicit none
c*********************************************************************
c     shifts the density or potential of all mt-cell into one single
c     vector and projects out the coulomb part only.
c                                    s. bluegel , kfa , 1987
c
c ------------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER IRMD,LPOTD
c      PARAMETER (IRMD=1484,LPOTD=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IMAP,INS,LMPOT,NATPS,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION X(IRMD,LMPOTD,*),Y(*)
      INTEGER IRC(*),IRMIN(*)
C     ..
C     .. Local Scalars ..
      INTEGER IA,IP,IR,IRC1,IRMIN1,IS,LM
C     ..
      IMAP = 0
      DO 50 IS = 1,NSPIN
        DO 40 IA = NATPS,NATYP
          IP = NSPIN* (IA-1) + IS
          IRC1 = IRC(IA)
          DO 10 IR = 1,IRC1
            IMAP = IMAP + 1
            Y(IMAP) = X(IR,1,IP)
   10     CONTINUE
c
          IF (INS.GT.0 .AND. LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IA)
            DO 30 LM = 2,LMPOT
              DO 20 IR = IRMIN1,IRC1
                IMAP = IMAP + 1
                Y(IMAP) = X(IR,LM,IP)
   20         CONTINUE
   30       CONTINUE
          END IF
c
   40   CONTINUE
   50 CONTINUE
c
      END
