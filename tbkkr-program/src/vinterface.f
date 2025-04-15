      SUBROUTINE VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NLAYERS,V,
     &                Z,R,IRWS,IRCUT,IPAN,KSHAPE,alat,
     &                bravais,ngmax,nrmax,nshlg,nshlr,nsg,nsr,gn,rm,vol,
     &                ZVEC,NLBASIS,NLEFT, ZPERLEFT,TLEFT,
     &                     NRBASIS,NRIGHT,ZPERIGHT,TRIGHT,
     &                     CMOMHOST,CHRGNT,
     &                     YRG,WG,ICC,ATOMIMP)

      implicit none
c**************************************************
c  This is calculating the intratomic contibution
c  of the potential in the case of an interface taking
c  into account the bulk potential on the 2 sides.
c     zvec(3,nlayers)    : coordinates of the layers
c     nlay               : number of layers
c     nlbasis            : number of basis layers of left
c                          host (repeated units)
c     nleft              : number of repeated basis
c                          for left host to get converged
c                          electrostatic potentials
c
c     zperleft(3)        : vector to define how to repeat
c                          the basis of the left host
c     tleft(3,nlbasis)   : vectors of the basis for the
c                          left host.
c     nrbasis            : number of basis layers of right
c                          host (repeated units)
c     nright             : number of repeated basis
c                          for right host to get converged
c                          electrostatic potentials
c
c     zperight(3)        : vector to define how to repeat
c                          the basis of the right host
c     tright(3,nlbasis)  : vectors of the basis for the
c                          right host.
c
c     cmomhol(lmpot,nlbasis) : charge moments of each atom
c                              of the left host
c
c     cmomhor(lmpot,nrbasis) : charge moments of each atom
c                              of the right host
c --------------------------------------------------
c ==================================================
      include 'inc.fi'
      integer LMPOTD
      parameter (LMPOTD=(LPOTD+1)**2)
      integer lm3d,ncleb1
      parameter (lm3d=(2*lpotd+1)**2)
      parameter (ncleb1=lm3d*lmpotd)
      INTEGER LASSLD
      PARAMETER (LASSLD=4*LMAXD)
c     ..
c     .. scalar arguments ..
      double precision vol,CHRGNT
      integer ngmax,nrmax,nshlg,nshlr,NLAY,NTOTAL,kshape
      integer nleft,nright,nlbasis,nrbasis,nlayers,iofile
c     ..
c     .. array arguments ..

      integer nsg(ishld),nsr(ishld)

      double precision gn2(2,nmaxd),rm2(2,nmaxd),br(3,3)
      double precision gn(4,nmaxd),rm(4,nmaxd),bravais(3,3)
      double precision ZVEC(3,*)
      double precision ZPERLEFT(3),ZPERIGHT(3),TLEFT(3,*)
      double precision TRIGHT(3,*)
      double precision CMINST(LMPOTD,*),CMOM(LMPOTD,*),
     &                 R(IRMD,*),V(IRMD,LMPOTD,*),Z(*),
     &                 CMOMHOST(LMPOTD,*)
      double precision CLEB(NCLEB1),AVMAD(LMPOTD,LMPOTD)
      double precision SUM(LM3D),SUMDIR(LM3D)
      double precision WG(LASSLD),
     +                 YRG(LASSLD,0:LASSLD,
     +                     0:LASSLD)
      INTEGER ICLEB(NCLEB1,3),IEND

      INTEGER IPAN(*),IRCUT(0:IPAND,*),IRWS(*)
      double precision AC(lmpotd),ACTOT(lmpotd),cm(LMPOTD)
      double precision vec1(3),vec2(3),diff,MONOPOL(NAEZD)
      double precision cm1,charge(2)

      integer lpot,ih1,ib,lm2,ilay1,iatom,lm1,i1,i2,nat
      integer n1,n2,lm,l,ilay,lmpot,ispin,nspin,irs1,ipot,i,m,j
      integer icc,iatommin,iatommax,natomimp,atomimp(natomimpd)
      integer nrefl,irefs(NAEZD)
      double precision rclsimp(3,natomimpd),vinters(lmpotd,naezd)
      double precision pi,fpi,alat
      logical OPT,TEST
!----------------------------------------------------------------------------
c       LMPOT = (LPOT+1)**2
c       DO I1=1,NLAYERS
c             READ(53,*) I2,Z(I2)
c             READ(53,8010) (CMOM(I2,I1),I2=1,LMPOT)
c       END DO
c8010   format(10F12.8)
c     Go from scaled units to a.u.
c
       write(6,*) '>>>>>> Vinterface'
       do i1=1,3
       do i=1,3
         br(i,i1) = bravais(i,i1)*alat
       end do
       end do
c
c Now br is in true a.u. units
c
       CALL MADELGAUNT(LPOT,YRG,WG,CLEB,ICLEB,IEND)
       pi = 4.d0*atan(1.d0)
       fpi = 4.d0*pi
       lmpot = (lpot+1)**2
       do i1 = 1,nmaxd
          do i=1,2
            gn2(i,i1) = gn(i,i1)
            rm2(i,i1) = rm(i,i1)
          end do
       end do

c     setup the charges to put in the ghost layers
c     in the case of decimation technique to achieve
c     charge neutrality



      if (OPT('DECIMATE')) then

            charge(1) = -chrgnt/(2.d0*dsqrt(fpi))
            charge(2) = -chrgnt/(2.d0*dsqrt(fpi))

c            write (6,*) 'charge left  = ',charge(1)
c            write (6,*) 'charge right = ',charge(2)

      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        START CALCULATION IN THE LAYERS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       Do ILAY=1,NLAYERS

          if (TEST('electro ')) write(29,*)
     &                          'NOW Calculating LAYER = ',ILAY
          if (kshape.ne.0) then
             irs1 = ircut(ipan(ilay),ilay)
          else
             irs1 = irws(ilay)
          end if

          do lm=1,lmpot
             AC(lm) = 0.d0
          end do
c
c     Make summation in all other layers
c     in three parts
c
c ***************************************************
c  all positions must be scaled with alat to get them
c  correct.
c  This is done in EWALD2D, SLAYDIR
c
c ***************************************************
c
c     1.  Summation in all layers in the slab
c
          do i=1,3
             VEC1(i) = ZVEC(i,ILAY)
          end do

          do ilay1 = 1,nlayers
          if (TEST('electro ')) write(29,*) 'with layer ***** ',ilay1

             IATOM = ILAY1   ! This is used on the cmoms

             do i=1,3
                VEC2(i) = ZVEC(i,ILAY1)
             end do

             SUMDIR = 0.d0
             diff = (vec1(3) - vec2(3))**2
             if (dabs(diff).lt.1.d-6) THEN
             CALL SLAYDIRECT(LPOT,VEC1,VEC2,ALAT,BR,RM2,NRMAX,NSHLR,NSR,
     &               SUMDIR)
             end if
c
c     make ewald sumation in plane and
c     Inverse space sum if rz<>0 (out of plane)
c
          CALL EWALD2D(LPOT,ALAT,VEC1,VEC2,RM2,NRMAX,NSHLR,NSR,
     &         gn2,ngmax,nshlg,nsg,SUM,SUMDIR,Vol)

          CALL MADELCOEFS(LPOT,AVMAD,SUM,CLEB,ICLEB,IEND)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Hoshino is doing (SMONOPOL(I) -SMONOPOL(0))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          IF (TEST('electro ')) THEN
             if (ILAY.eq.1.and.ilay1.le.5) THEN
                do lm1=1,lmpotd
                   if (dabs(SUM(lm1)).gt.1.d-8) then
                      write(6,*) 'SUMDIR',lm1,sum(lm1),sumdir(lm1)
                   end if
                end do
                do lm1=1,lmpotd

                   do lm2=1,lmpotd
                      if (dabs(avmad(lm1,lm2)).gt.1.d-8) then
                         write(29,*) 'lm1, lm2',lm1,lm2,avmad(lm1,lm2)
                      end if
                   end do
                end do
             end if             ! ILAY.eq.1.....
          end if                ! electro
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccc
c Keep the monopole term
cccccccccccccccccccccccccccccccccc
          IF (ILAY1.eq.ILAY) THEN
              MONOPOL(ILAY1) = AVMAD(1,1)
          END IF
          do lm=1,lmpot
             cm(lm) = CMOM(LM,IATOM)
c---  > add contribution of interstial in case of shapes
             if (kshape.ne.0) cm(lm) = cm(lm) + CMINST(LM,IATOM)
          end do

          CM(1) = CM(1) - Z(IATOM)/sqrt(fpi)
c
          Do LM =1,LMPOT
             Do LM2 = 1,LMPOT
                AC(LM) = AC(LM) + AVMAD(LM,LM2)*CM(LM2)
             End Do
          End Do

       end do                   ! ilay1 loop in all interface planes
c     #####################################################
       DO ILAY1=1,NLAYERS

           IATOM = ILAY1           ! before
          cm1 = CMOM(1,IATOM)
          if (kshape.ne.0) cm1 = cm1 + CMINST(1,IATOM)
          cm1 = CM1 - Z(IATOM)/sqrt(fpi)
          AC(1) = AC(1) - MONOPOL(ILAY)*CM1
       END DO

c     Correction look for example at P.Lang charge neutrality
c     is imposed ....
c     ####################################################
c
c     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
c     2.  Summation in the LEFT bulk side
c     vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

       IF (OPT('DECIMATE')) THEN
c
c
         do ih1 = 1,NLEFT
             do ib = 1,NLBASIS
                do i=1,3
                   VEC2(i) = (tleft(i,ib) + (ih1-1)*zperleft(i))
                end do

                IATOM = IB

                SUMDIR = 0.d0
                diff = (vec1(3) - vec2(3))**2 ! changed on 29.11.99
                if (dabs(diff).lt.1.d-6) THEN
             CALL SLAYDIRECT(LPOT,VEC1,VEC2,ALAT,BR,RM2,NRMAX,NSHLR,NSR,
     &                  SUMDIR)
                end if
c
c     make ewald sumation for m= 0 l<5 rz=0 (in plane)and
c     Inverse space sum if rz<>0 (out of plane)
c
                CALL EWALD2D(LPOT,ALAT,VEC1,VEC2,RM2,NRMAX,NSHLR,NSR,
     &               gn2,ngmax,nshlg,nsg,SUM,SUMDIR,Vol)
                CALL MADELCOEFS(LPOT,AVMAD,SUM,CLEB,ICLEB,IEND)
                do lm=1,lmpot
                   cm(lm) = CMOMHOST(LM,IATOM)
                end do
c
                Do Lm =1,LMpot
                   Do LM2 = 1,LMpot
                      AC(LM) = AC(LM) + AVMAD(LM,LM2)*CM(LM2)
                   End Do
                End Do
                                ! Added 19.11.99
                if (OPT('DECIMATE').and.(ih1.eq.1).and.(ib.eq.1)) then
                   AC(1) = AC(1) + (AVMAD(1,1)-MONOPOL(ILAY))*CHARGE(1)
c     write (6,*) 'monopol contr. from the left ghost layer'
c     write (6,*) '-> ',(AVMAD(1,1)-MONOPOL(ILAY))*CHARGE(1)
                endif
                                ! Added 19.11.99

             end do             ! ib loop in right host basis
          end do                ! ih1 loop in layers to get convergence


c
c     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c     3.  Summation in the RIGHT bulk side
c     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
c
          do ih1 = 1,NRIGHT
             do ib = 1,NRBASIS
                do i=1,3
                   VEC2(i) = (tright(i,ib) + (ih1-1)*zperight(i))
                end do

                IATOM = NLBASIS + IB ! corrected 27.10.99

                SUMDIR = 0.d0

c     diff = (vec1(1) - vec2(1))**2 + (vec1(2) - vec2(2))**2 +
c     &            (vec1(3) - vec2(3))**2
                diff =  (vec1(3) - vec2(3))**2 ! changed on 29.11.99
                if (dabs(diff).lt.1.d-6) THEN
             CALL SLAYDIRECT(LPOT,VEC1,VEC2,ALAT,BR,RM2,NRMAX,NSHLR,NSR,
     &                  SUMDIR)
                end if
c
c     make ewald sumation (in plane)and
c     Inverse space sum if rz<>0 (out of plane)
c
                CALL EWALD2D(LPOT,ALAT,VEC1,VEC2,RM2,NRMAX,NSHLR,NSR,
     &               gn2,ngmax,nshlg,nsg,SUM,SUMDIR,Vol)
                CALL MADELCOEFS(LPOT,AVMAD,SUM,CLEB,ICLEB,IEND)
                do lm=1,lmpot
                   cm(lm) = CMOMHOST(LM,IATOM)
c---  > add contribution of interstial in case of shapes
                end do

                Do Lm =1,LMpot
                   Do LM2 = 1,LMpot
                      AC(LM) = AC(LM) + AVMAD(LM,LM2)*CM(LM2)
                   End Do
                End Do
                                ! Added 19.11.99
                if (OPT('DECIMATE').and.(ih1.eq.1).and.(ib.eq.1)) then
                   AC(1) = AC(1) + (AVMAD(1,1)-MONOPOL(ILAY))*CHARGE(2)
c     write (6,*)'monopol contr. from the right ghost layer'
c     write (6,*)'-> ',(AVMAD(1,1)-MONOPOL(ILAY))*CHARGE(2)
                endif
                                ! Added 19.11.99
             end do             ! ib loop in right host basis
          end do                ! ih1 loop in layers to get convergence
c
c
       END IF                   ! (OPT(DECIMATE)
c
c     -----------------------------------------------------------------
       WRITE (6,FMT=9000) ILAY, (AC(1)/SQRT(4.D0*PI)),
     &                             (AC(3)/SQRT(4.D0*PI))
       DO  ISPIN = 1,NSPIN
c
c---  > determine the right potential number
c
          IPOT = NSPIN* (ILAY-1) + ISPIN
c
c---  > in the case of l=0 : r(1)**l is not defined
c

             V(1,1,IPOT) = V(1,1,IPOT) + AC(1)

             do L = 0,LPOT

                do M=-L,L
                   LM = L*L + L + m + 1
                   DO I = 2,IRS1

                   V(I,LM,IPOT) = V(I,LM,IPOT) + (-R(I,ILAY))**L*AC(LM)

                   End do
                end do
             end do
          End Do                ! ISPIN LOOP

          if (icc.gt.0) then
             do L = 0,LPOT
                do M=-L,L
                   LM = L*L + L + m + 1
                   vinters(lm,ilay) = ac(lm)
                enddo
             enddo
          endif

       END DO                   ! ILAY loop
c
c Now Prepare output for Impurity calculation
c
       if (icc.gt.0) then

C     it reads again the informations about the cluster
C     around the impurity.

          REWIND 25
          READ (25,FMT=*) NATOMIMP
          DO I=1,NATOMIMP
             READ (25,FMT=*) (RCLSIMP(J,I),J=1,3)!,ATOMIMP(I)
             IREFS(I) = ATOMIMP(I) +ICC
          ENDDO

          NREFL=1
          DO I = 1,NATOMIMP
             J = ICC + ATOMIMP(I)
             IF (J.NE.IREFS(NREFL)) THEN
                NREFL = NREFL + 1
                IREFS(NREFL) = J
             END IF
          ENDDO


c         call fxdropn('reference_pot','ENCODE',iofile)

          do i=1,nrefl
c            call fxdrdbl(iofile,vinters(1,irefs(i)),lmpotd)

             do lm=1,lmpotd
             write (6,*) lm,irefs(i)
             write (6,*) vinters(lm,irefs(i))
             enddo
          end do

       end if

       RETURN

 9000  FORMAT (1x,'madelung monopole and dipole for layer',i3,
     +      ' :',1p,2d13.6)
 9001  FORMAT (1x,'madelung before host added   for layer',i3,
     +      ' :',1p,2d12.5)
       END
