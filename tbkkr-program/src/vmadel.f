c 13.10.95 ***************************************************************
      SUBROUTINE VMADEL(AVMAD,BVMAD,CMOM,CMINST,LMAX,NSPIN,NATYP,V,Z,R,
     +                  IRWS,IRCUT,IPAN,KSHAPE)
      implicit none
c ************************************************************************
c     calculate the vmadelung-potentials and add these to the poten-
c     tial v  (in the spin-polarized case for each spin-direction
c     this is the same . )
c     it uses the structure dependent matrices avmad and bvmad which
c     are calculate once in the subroutine amn .
c     the charge-moments are calculated in the subroutine vintra2 ,
c     therefore vintra2 has to be called first .
c     the madelung-potential is expanded into spherical harmonics .
c     the lm-term of the potential v of the atom i is given by
c
c      v(r,lm,i) =  (-r)**l * {avmad(i,i2,lm,l'm')*cmom(i2,l'm')
c                                               +bvmad(i,i2,lm)*z(i2)}
c
c     summed over i2 (all atoms) and l'm' .
c             (see notes by b.drittler)
c
c                               b.drittler   nov. 1989
c-----------------------------------------------------------------------
C     .. Parameters ..
      include 'inc.fi'
c      INTEGER IRMD,LPOTD
c      PARAMETER (IRMD=1484,LPOTD=8)
c      INTEGER IPAND
c      PARAMETER (IPAND=4)
c      INTEGER NATYPD
c      PARAMETER (NATYPD=1)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMAX,NATYP,NSPIN
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION AVMAD(NATYPD,NATYPD,LMPOTD,*),
     +                 BVMAD(NATYPD,NATYPD,*),CMINST(LMPOTD,*),
     +                 CMOM(LMPOTD,*),R(IRMD,*),V(IRMD,LMPOTD,*),Z(*)
      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRWS(*)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION AC,PI
      INTEGER I,I1,I2,IPOT,IRS1,ISPIN,L,LM,LM2,LMMAX,M
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DATAN,SQRT
C     ..
      PI = 4.0D0*DATAN(1.0D0)
      LMMAX = (LMAX+1)* (LMAX+1)
c       DO I1=1,LMMAX
c             WRITE (53,8010) I1,(CMOM(I1,I2),I2=1,4)
c             WRITE (53,8010) I1,(CMINST(I1,I2),I2=1,4)
c       END DO
c8010   format(I3,10F12.8)
c       STOP

      DO 100 I1 = 1,NATYP
        IF (KSHAPE.NE.0) THEN
          IRS1 = IRCUT(IPAN(I1),I1)

        ELSE
          IRS1 = IRWS(I1)
        END IF
c
        DO 90 L = 0,LMAX
          DO 80 M = -L,L
            LM = L*L + L + M + 1
            AC = 0.0D0

            IF (NATYP.EQ.1) THEN
c---> lm = 1 component disappears if there is only one host atom
              DO 10 LM2 = 2,LMMAX
c---> take moments of sphere
                AC = AC + AVMAD(I1,1,LM,LM2)*CMOM(LM2,1)
   10         CONTINUE
              IF (KSHAPE.NE.0) THEN
                DO 20 LM2 = 2,LMMAX
c---> add contribution of interstial in case of shapes
                  AC = AC + AVMAD(I1,1,LM,LM2)*CMINST(LM2,1)
   20           CONTINUE
              END IF

            ELSE

              DO 50 I2 = 1,NATYP
                AC = AC + BVMAD(I1,I2,LM)*Z(I2)
                DO 30 LM2 = 1,LMMAX
c---> take moments of sphere
                  AC = AC + AVMAD(I1,I2,LM,LM2)*CMOM(LM2,I2)
   30           CONTINUE
                IF (KSHAPE.NE.0) THEN
                  DO 40 LM2 = 1,LMMAX
c---> add contribution of interstial in case of shapes
                    AC = AC + AVMAD(I1,I2,LM,LM2)*CMINST(LM2,I2)
   40             CONTINUE
                END IF

   50         CONTINUE
            END IF

            IF (LM.EQ.1) WRITE (6,FMT=9000) I1, (AC/SQRT(4.D0*PI))
c
c---> add to v the intercell-potential
c
            DO 70 ISPIN = 1,NSPIN
c
c---> determine the right potential number
c
              IPOT = NSPIN* (I1-1) + ISPIN
c
c---> in the case of l=0 : r(1)**l is not defined
c
              IF (L.EQ.0) V(1,1,IPOT) = V(1,1,IPOT) + AC
              DO 60 I = 2,IRS1
                V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,I1))**L*AC
   60         CONTINUE
   70       CONTINUE
   80     CONTINUE
   90   CONTINUE
  100 CONTINUE

      RETURN

 9000 FORMAT (1x,'spherically averaged madelung-potential for atom',i3,
     +       ' :',1p,d14.6)
      END
