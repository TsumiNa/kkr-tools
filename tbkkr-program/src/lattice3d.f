



      subroutine lattice3d(alat,bravais,ngmax,nrmax,nshlg,nshlr,nsg,
     +                    nsr,gn,rm,natom,vol)
      implicit none
c-----------------------------------------------------------------------
c
c     generate lattice vectors of direct and reciprocal space from
c     basic translation vectors br
c
c     alat            : lattice constant
c     br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors
c                       *** in a.u. **** 
c     rmax            : maximum radius in real space
c     gmax            : maximum radius in reciprocal space
c     ngmax           : Number of reciprocal lattice vectors
c     gn(4,nmaxd)     : x,y,z,r of reciprocal lattice vectors
c     nrmax           : Number of real lattice vectors
c     rm(4,nmaxd)     : x,y,z,r of real space vectors
c     nshlg           : shells in reciprocal space
c     nshlr           : shells in real space
c     nsg,nsr         : integer arrays, number of atoms in each shell
c     vol             : Elementary cell volume 
c-----------------------------------------------------------------------
 
c     .. parameters ..
      !integer nmaxd,ishld
      !PARAMETER (NMAXD=9500,ISHLD=500)
      include 'inc.fi'
c     ..
c     .. scalar arguments ..
      double precision alat,vol
      integer iprint,ngmax,nrmax,nshlg,nshlr,natom
c     ..
c     .. array arguments ..
      double precision gn(4,nmaxd),rm(4,nmaxd),bravais(3,3)
      integer nsg(ishld),nsr(ishld)
c     ..
c     .. local scalars ..
      double precision a,absgm,absrm,ag,ar,b,c,da,db,gmax,gx,gy,gz,pi,
     +                 rmax,rx,ry,rz,
     +                 vmin
      integer i,i1,i2,k,l,m,n,n1,ng,nr,nsh,nshl,numg,numgh,numr,numrh
      INTEGER IER
c     ..
c     .. local arrays ..
      double precision absg(3),absr(3),bg(3,3),br(3,3),cj(4,nmaxd)
      CHARACTER*80 UIO
c     ..
c     .. logical ..
      LOGICAL TEST
c     ..
c     .. intrinsic functions ..
      intrinsic abs,datan,idint,max,mod,real,sqrt
c     ..
      EXTERNAL TEST,IOINPUT
      pi = 4.0d0*datan(1.0d0)
c
c---> read lattice constant , no. of atoms per unit cell and cutoffs
c
c      read (5,*) rmax,gmax
c     read (5,fmt=9000) alat,natom,rmax,gmax
      CALL IoInput('RMAX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) RMAX
      CALL IoInput('GMAX      ',UIO,1,7,IER)
                      READ (UNIT=UIO,FMT=*) GMAX
      write(6,*) '>>>>>>>> lattice3d '
      write (6,fmt=9010) alat,rmax,gmax
      rmax = rmax*alat
      gmax = gmax/alat
c
c---> basic trans. vectors and basis vectors
c
      iprint = 0
      do i=1,3
         br(1,i) = bravais(1,i)*alat
         br(2,i) = bravais(2,i)*alat
         br(3,i) = bravais(3,i)*alat
      end do
      !if (TEST('latshell')) iprint = 1
      if (TEST('latshell'))  THEN
         write (6,fmt=9020)
         write(6,*) 'normalized bravais in lattice3d'
         do i = 1,3
            write (6,fmt=9040) bravais(1,i),bravais(2,i),bravais(3,i)
         END DO
         write(6,*) 'true bravais in lattice3d'
         DO I=1,3 
            write (6,fmt=9040) br(1,i),br(2,i),br(3,i)
         END DO
         
      END IF                    ! test latshell
c     
c---> generate primitive vectors bg of reciprocal space
c
         do 30 i = 1,3
            i1 = 1 + mod(i,3)
            i2 = 1 + mod(i1,3)
c
c---> cross product
c
            bg(1,i) = br(2,i1)*br(3,i2) - br(2,i2)*br(3,i1)
            bg(2,i) = br(3,i1)*br(1,i2) - br(3,i2)*br(1,i1)
            bg(3,i) = br(1,i1)*br(2,i2) - br(1,i2)*br(2,i1)
   30    continue
c
         vol = abs(br(1,1)*bg(1,1)+br(2,1)*bg(2,1)+br(3,1)*bg(3,1))
c
         if (TEST('latshell'))  write (6,fmt=9050)
         do 40 i = 1,3
            bg(1,i) = bg(1,i)/vol*2.0d0*pi
            bg(2,i) = bg(2,i)/vol*2.0d0*pi
            bg(3,i) = bg(3,i)/vol*2.0d0*pi
            if (TEST('latshell')) 
     &      write (6,fmt=9040) bg(1,i),bg(2,i),bg(3,i)
   40    continue
c
c---> estimate no. of lattice vectors
c
         do 50 i = 1,3
            absr(i) = sqrt(br(1,i)**2+br(2,i)**2+br(3,i)**2)
            absg(i) = sqrt(bg(1,i)**2+bg(2,i)**2+bg(3,i)**2)
   50    continue
         absrm = max(absr(1),absr(2),absr(3))
         absgm = max(absg(1),absg(2),absg(3))
         absrm = 2.0d0*pi/absrm
         absgm = 2.0d0*pi/absgm
         numr = 2* (idint(rmax/absgm)+1) + 1
         numg = 2* (idint(gmax/absrm)+1) + 1
         numrh = numr/2 + 1
         numgh = numg/2 + 1
         if (TEST('latshell')) write (6,fmt=9110) numr,numg
c
c---> generate lattice vectors of real space
c
         if (iprint.gt.0) write (6,fmt=9060)
         nr = 0
         do 60 l = 1,numr
            a = real(l-numrh)
            do 70 m = 1,numr
               b = real(m-numrh)
               do 80 n = 1,numr
                  c = real(n-numrh)
                  rx = a*br(1,1) + b*br(1,2) + c*br(1,3)
                  ry = a*br(2,1) + b*br(2,2) + c*br(2,3)
                  rz = a*br(3,1) + b*br(3,2) + c*br(3,3)
                  ar = sqrt(rx*rx+ry*ry+rz*rz)
                  if (ar.le.rmax) then
                     nr = nr + 1
                     cj(1,nr) = rx
                     cj(2,nr) = ry
                     cj(3,nr) = rz
                     cj(4,nr) = ar
                  end if
 
   80          continue
 
   70       continue
 
   60    continue
c
         if (nr.gt.nmaxd) then
            write(6,*) nr,nmaxd
            stop ' 1 - lattice3d '
 
         else
 
            nrmax = nr
c
c---> sort vectors in order of increasing absolute value
c
            da = 1.d-06
            nsh = 0
            nshl = -1
            do 90 k = 1,nr
               vmin = rmax + 1.0d0
               do 100 n = 1,nr
                  if (cj(4,n)-vmin.lt.0) then
                     vmin = cj(4,n)
                     n1 = n
                  end if
 
  100          continue
 
 
               nshl = nshl + 1
               rm(1,k) = cj(1,n1)
               rm(2,k) = cj(2,n1)
               rm(3,k) = cj(3,n1)
               rm(4,k) = cj(4,n1)
               db = vmin
               if (db.gt.da+1.d-06) then
 
                  nsh = nsh + 1
                  if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                  nsr(nsh) = nshl
                  if (iprint.gt.0) write (6,fmt=9070) k,rm(1,k),rm(2,k),
     +                rm(3,k),db
                  nshl = 0
                  da = db
 
               else if (iprint.gt.0) then
                  write (6,fmt=9070) k,rm(1,k),rm(2,k),rm(3,k),db
               end if
 
               cj(4,n1) = rmax + 1.0d0
   90       continue
            nsh = nsh + 1
            nshl = nshl + 1
            nsr(nsh) = nshl
            if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
            if (nsh.gt.ishld) then
               write(6,*) nsh,ishld
               stop ' 2 - lattice3d '
 
            else
 
               nshlr = nsh
c
c---> generate lattice vectors of reciprocal space
c
               if (iprint.gt.0) write (6,fmt=9090)
               ng = 0
               do 110 l = 1,numg
                  a = real(l-numgh)
                  do 120 m = 1,numg
                     b = real(m-numgh)
                     do 130 n = 1,numg
                        c = real(n-numgh)
                        gx = a*bg(1,1) + b*bg(1,2) + c*bg(1,3)
                        gy = a*bg(2,1) + b*bg(2,2) + c*bg(2,3)
                        gz = a*bg(3,1) + b*bg(3,2) + c*bg(3,3)
                        ag = sqrt(gx*gx+gy*gy+gz*gz)
                        if (ag.le.gmax) then
                           ng = ng + 1
                           cj(1,ng) = gx
                           cj(2,ng) = gy
                           cj(3,ng) = gz
                           cj(4,ng) = ag
                        end if
 
  130                continue
 
  120             continue
 
  110          continue
c
               if (ng.gt.nmaxd) then
                  write(6,*) ng,nmaxd
                  stop ' 3 - lattice3d '
 
               else
 
                  ngmax = ng
c
c---> sort vectors in order of increasing abs. value
c
                  da = 1.d-06
                  nsh = 0
                  nshl = -1
                  do 140 k = 1,ng
                     vmin = gmax + 1.0d0
                     do 150 n = 1,ng
                        if (cj(4,n)-vmin.lt.0) then
                           vmin = cj(4,n)
                           n1 = n
                        end if
 
  150                continue
 
 
                     nshl = nshl + 1
                     gn(1,k) = cj(1,n1)
                     gn(2,k) = cj(2,n1)
                     gn(3,k) = cj(3,n1)
                     gn(4,k) = cj(4,n1)
                     db = vmin
                     if (db.gt.da+1.d-07) then
 
                        nsh = nsh + 1
                        if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                        if (iprint.gt.0) write (6,fmt=9070) k,gn(1,k),
     +                      gn(2,k),gn(3,k),db
                        nsg(nsh) = nshl
                        nshl = 0
                        da = db
 
                     else if (iprint.gt.0) then
                        write (6,fmt=9070) k,gn(1,k),gn(2,k),gn(3,k),db
                     end if
 
                     cj(4,n1) = gmax + 1.0d0
  140             continue
                  nsh = nsh + 1
                  nshl = nshl + 1
                  nsg(nsh) = nshl
                  if (iprint.gt.0) write (6,fmt=9080) nsh,nshl
                  if (nsh.gt.ishld) then
                     write(6,*) nsh,ishld
                     stop ' 4 - lattice3d '
 
                  else
 
                     nshlg = nsh
 
                     write (6,fmt=9100) nrmax,ngmax,vol
 
 
 
                  end if
 
               end if
 
            end if
 
         end if
 

 
 
 9000 format (f10.5,i5,2f10.5)
 9010 format (/,1x,' lattice constant : ',f10.5,/,7x,' rmax : ',f10.5,
     +       ' gmax : ',f10.5)
 9020 format (/,1x,' primitive vectors of direct lattice ',/)
 9030 format (/,1x,i3,' basis vectors per unit cell ',/)
 9040 format (3f10.5)
 9050 format (/,1x,' primitive vectors of reciprocal lattice ',/)
 9060 format (/,1x,' real space lattice vectors ',/)
 9070 format (5x,i5,4f10.5)
 9080 format (1x,' shell no. ',i3,' contains ',i3,' points ',/)
 9090 format (/,1x,' reciprocal space lattice vectors ',/)
 9100 format (/,1x,' no. of lattice vectors   : ',i4,/,1x,
     +       ' no. of rec. lat. vectors : ',i4,/,1x,
     +       ' volume of the unit cell  : ',f10.5)
 9110 format (/,5x,' numr ',i3,' numg ',i3)
 
      end
