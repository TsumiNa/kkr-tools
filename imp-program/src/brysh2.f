      SUBROUTINE BRYSH2(Y,X,xsme,INS,IRMIN,IRC,NATPS,NATYP,NSPIN,
     $     IMAP,LMPOT,lsmear)
      Implicit None
c*********************************************************************
c     maps the density or potential back from one single vector into
c     the proper bins of each single mt-cell . the magnetization
c     density is also added in.
c                                    s. bluegel , kfa , 1987
c
c*********************************************************************
C     .. Parameters ..
      INTEGER IRMD,LPOTD
      PARAMETER (irmd=1484,lpotd=8)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IMAP,INS,LMPOT,NATPS,NATYP,NSPIN,lsmear
C     ..
C     .. Array Arguments ..
      REAL*8 X(IRMD,LMPOTD,*),Y(*),xsme(irmd,*)
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
            X(IR,1,IP) = Y(IMAP)
   10     CONTINUE
c
c     next for smeared spherical potential
          If ( lsmear .gt. 0 ) Then
            Do ir = 1, irc1
              imap = imap + 1
              xsme(ir,ip) = y(imap)
            End Do
          End If
c
          IF (INS.GT.0 .AND. LMPOT.GT.1) THEN
            IRMIN1 = IRMIN(IA)
            DO 30 LM = 2,LMPOT
              DO 20 IR = IRMIN1,IRC1
                IMAP = IMAP + 1
                X(IR,LM,IP) = Y(IMAP)
   20         CONTINUE
   30       CONTINUE
          END IF
c
   40   CONTINUE
   50 CONTINUE
c
      END
