c number of type of atoms 
      integer natypd,NATOMD,NTREFD,NTPERD
      parameter (NATYPD=13,NATOMD=201,NTREFD=1,NTPERD=NATYPD-NTREFD)
c highest valence orbital quantum number     
      integer lmaxd,lmx,lpotd
      parameter (LMAXD=4,lmx=lmaxd+1,lpotd=2*lmaxd)
c higheset number of seqular equation
      integer nsec,NREPD
      parameter (nsec=339,NREPD=10)
C     MESH NUMBER OF ENERGY CONTOUR
      INTEGER IEMXD
      PARAMETER (iemxd=128)      
c number of spin
      integer nspind
      parameter (NSPIND = 2)
c     see output: newsym + GREEN-FUNCTION
      integer NLSTD,ICGD,ICJD
      parameter(NLSTD=11686,ICGD=6000,ICJD=93139)
c     radial-mesh
      integer irmd,IRNSD,irmind
      parameter (irmd=884,irnsd=608,irmind=irmd-irnsd)
c     shape
      integer nfund,irid,ngshd
      parameter (nfund=289,irid=535,ngshd=54287)
c


