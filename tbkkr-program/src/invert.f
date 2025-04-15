c ************************************************************************
      SUBROUTINE INVERT(NL,DIM,M1,M2,M3,IOPT,M4,M5)
      implicit none
c ************************************************************************
c     INPUT :
c
c     nl       dimension of tridiagonal matrix 
c     dim      dimension of matrix elements
c     m1,m2,m3 matrix elements of the tridiagonal matrix to invert,
c
c              storage of the diagonal elements
c
c                    M(n     ,n     )  -->  M2(n)
c                    M(n+1   ,n     )  -->  M3(n)
c                    M(n     ,n+1   )  -->  M1(n)
c
c              storage of the corner elements
c
c                    M(1     ,nlayer)  -->  M3(nlayer)
c                    M(nlayer,1)       -->  M1(nlayer)
c
c     ALGORITHM :
c
c     iopt     =  0, 1, 2, 3   supercell algorithm
c                10,11,  (13)   slab algorithm
c
c     FUNCTION :
c
c     iopt     =  0,10      determination of diagonal elements of M**(-1)
c
c                 1,11      factorization M = L**(-1) D**(-1) U**(-1)
c                           to count negative eigenvalues
c
c                 2,12      factorization M = L**(-1) D**(-1) U**(-1)
c                           for eigenvector determination
c
c                 3,(13)    determination of cluster GF elements
c
c     OUTPUT :
c
c     iopt = 0,10    M**(-1)(n,n)   -->  M2(n)      n = 1,nl 
c
c            1,11    D**(-1)(n)     -->  M2(n)      n = 1,nl 
c
c            2,12    D**(-1)(n)     -->  M2(n)      n = 1,nl 
c                    C(n)           -->  M3(n)      n = 1,nl 
c
c            3,(13)  M**(-1)(n  ,n)          -->  M2(n)    n = 1,nl 
c                    M**(-1)(n  ,n+1)=Y(n)   -->  M1(n)    n = 1,nl-1
c                    M**(-1)(n  ,n+2)=YY(n)  -->  M3(n)    n = 1,nl-2
c                    M**(-1)(n+1,n  )=X(n)   -->  M4(n)    n = 1,nl-1
c                    M**(-1)(n+2,n  )=XX(n)  -->  M5(n)    n = 1,nl-2
c
c ------------------------------------------------------------------------
c     .. Parameters ..
      include 'inc.fi'
c
      INTEGER LMMAXD,NDIM
      PARAMETER (LMMAXD=(LMAXD+1)**2,NDIM=LMMAXD*NPRINCD)
      INTEGER NAUX
      PARAMETER (NAUX=2*NDIM**2+5*NDIM)
C
c     .. scalar arguments
      INTEGER DIM,NL,ICC,IOPT
C
c     .. array arguments
      DOUBLE COMPLEX 
     +     M1(NDIM,NDIM,*),
     +     M2(NDIM,NDIM,*),
     +     M3(NDIM,NDIM,*),
     +     M4(NDIM,NDIM,*),
     +     M5(NDIM,NDIM,*)
C
c     .. local scalars
      INTEGER I,J,N,LM,INFO,LM1,LM2,LM3,
     +        SAVC,SAVD,SAVY,SLAB
c
c     .. local arrays
      INTEGER IPVT(NDIM)
      DOUBLE COMPLEX A(NDIM,NDIM),
     +               B(NDIM,NDIM,NLAYERD),
     +               C(NDIM,NDIM,NLAYERD),
     +               D(NDIM,NDIM,NLAYERD),
     +               E(NDIM,NDIM),
     +               F(NDIM,NDIM),
     +               G(NDIM,NDIM),
     +               V(NDIM,NDIM),V1(NDIM,NDIM),V2(NDIM,NDIM),
     +               W(NDIM,NDIM),W1(NDIM,NDIM),W2(NDIM,NDIM),
     +               X(NDIM,NDIM),
     +               XX(NDIM,NDIM),
     +               Y(NDIM,NDIM),
     +               YY(NDIM,NDIM),
     +               CUNIT(NDIM,NDIM),
     +     Z(NDIM,NDIM),
     +     AUX(NAUX)
c ------------------------------------------------------------------------
c     .. Data Statements ..
      DOUBLE COMPLEX CZERO,CONE,BOUND
      PARAMETER (CZERO=(0.0D0,0.0D0),CONE=(1.0D0,0.0D0))
c
c     .. External Subroutines ..
      LOGICAL TEST
      EXTERNAL CINIT,TEST,
     +         ZCOPY,ZGEEV,ZGEMM,
c     +         ZGEMUL,
     +         ZGETRF,ZGETRS
      INTRINSIC ABS,DIMAG,DREAL
      DATA BOUND/1.0D-12/
      save
c
c ------------------------------------------------------------------------
c
c ---> translate IOPT to internal execution option SLAB,INV,SAVC,SAVD,SAVY
c
      IF (IOPT.EQ.13) THEN
c
c --->  take supercell algorithm instead of slab algorithm in case of 
c       nondiagonal element determination 
c
        IOPT = 3                   
      END IF

      SLAB = 0                      ! slab algorithm
      IF (IOPT.GE.10) SLAB = 1      
c
      SAVD = 0                      ! save D**(-1)(n) in M2(n)
      IF (IOPT.EQ.1 .OR. IOPT.EQ.2 
     +     .OR. IOPT.EQ.11 .OR. IOPT.EQ.12) SAVD = 1
c
      SAVC = 0                      ! save C(n) in M3(n)
      IF (IOPT.EQ.2 .OR. IOPT.EQ.12 ) SAVC = 1
C
      SAVY = 0                      ! save Y(n) in M1(n), YY(n) in M3(n)
      IF (IOPT.EQ.3) SAVY = 1
c ------------------------------------------------------------------------
c
c ---> test of dimensions
c
      IF (NL.GT.NLAYERD .OR. DIM.NE.NDIM) THEN
        IF (NL.GT.NLAYERD) THEN
          write(6,*) 'Number of layer (',NL,
     +         ') is greater than defined (',NLAYERD,'). '
          write(6,*) 'Increase the parameter NLAYERD. '
          write(6,*) 
        END IF
        
        IF (DIM.NE.NDIM) 
     +       write(6,*) 'Dimension of matrix elements are not equal.',
     +       ' DIM= ',DIM,', NDIM= ',NDIM,'.'

        STOP 'INVERT   '
      END IF

c ------------------------------------------------------------------------
c
c ---> initialize local arrays
c
      call cinit(ndim*ndim,a)
      call cinit(ndim*ndim*nlayerd,b)
      call cinit(ndim*ndim*nlayerd,c)
      call cinit(ndim*ndim*nlayerd,d)
      call cinit(ndim*ndim,e)
      call cinit(ndim*ndim,f)
      call cinit(ndim*ndim,v)
      call cinit(ndim*ndim,w)
      call cinit(ndim*ndim,cunit)
c ------------------------------------------------------------------------
c
c ---> cunit = complex unity matrix of order NDIM
c
      DO 10 N = 1,NDIM
        CUNIT(N,N) = CONE
 10   CONTINUE
c ************************************************************************
C
      IF (SLAB.EQ.0) THEN
C
c ************************************************************************
c
c ---> ALGORITM FOR SUPERCELL GEOMETRY
c
c ------------------------------------------------------------------------
c
c ---> factorization D ^-1 = (prod L) * M * (prod U)
c      
c      see notes R. Zeller
c
c ------------------------------------------------------------------------
c
c ---> N =1
c
      CALL ZCOPY(NDIM*NDIM,M2(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NDIM*NDIM,CUNIT,1,D(1,1,1),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,1),NDIM,INFO)

      IF (NL.EQ.1) GOTO 980

      CALL ZCOPY(NDIM*NDIM,M2(1,1,NL),1,A(1,1),1)
      CALL ZCOPY(NDIM*NDIM,M1(1,1,NL),1,B(1,1,1),1)
      CALL ZCOPY(NDIM*NDIM,M3(1,1,NL),1,C(1,1,1),1)

c ------------------------------------------------------------------------
c
c ---> 2 <= N < NL-1
c
      IF (NL.EQ.2) GOTO 970

      IF (NL.EQ.3) GOTO 960

      DO 20 N = 2,NL-2
        
c
c ---> E = D(N-1) * C(N-1)
c
c        CALL ZGEMUL(D(1,1,N-1),NDIM,'N',C(1,1,N-1),NDIM,'N',
c     +       E(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +       C(1,1,N-1),NDIM,CZERO,E(1,1),NDIM)
c
c ---> F = D(N-1) * M1(N-1)
c
c        CALL ZGEMUL(D(1,1,N-1),NDIM,'N',M1(1,1,N-1),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +       M1(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)
        
c
c ---> A = A - B(N-1)*D(N-1)*C(N-1)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +       E(1,1),NDIM,CONE,A(1,1),NDIM)
        
c
c ---> B(N) = - B(N-1)*D(N-1)*M1(N-1)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +       F(1,1),NDIM,CZERO,B(1,1,N),NDIM)
        
c
c ---> C(N) = - M3(N-1)*D(N-1)*C(N-1)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +       E(1,1),NDIM,CZERO,C(1,1,N),NDIM)
c
c ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
c
      CALL ZCOPY(NDIM*NDIM,M2(1,1,N),1,E(1,1),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +     F(1,1),NDIM,CONE,E(1,1),NDIM)

      IF (SAVD.GT.0)
     +     CALL ZCOPY(NDIM*NDIM,E(1,1),1,M2(1,1,N),1)

      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,N),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,N),NDIM,INFO)

 20   CONTINUE
c ------------------------------------------------------------------------
c
c ---> N = NL - 1
c
c
 960  N = NL - 1
c
c ---> E = D(N-1) * C(N-1)
c
c        CALL ZGEMUL(D(1,1,N-1),NDIM,'N',C(1,1,N-1),NDIM,'N',
c     +       E(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +       C(1,1,N-1),NDIM,CZERO,E(1,1),NDIM)
c
c ---> F = D(N-1) * M1(N-1)
c
c        CALL ZGEMUL(D(1,1,N-1),NDIM,'N',M1(1,1,N-1),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +       M1(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)
        
c
c ---> A = A - B(N-1)*D(N-1)*C(N-1)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +       E(1,1),NDIM,CONE,A(1,1),NDIM)
        
c
c ---> B(N) = - B(N-1)*D(N-1)*M1(N-1) + M3(N)
c
        CALL ZCOPY(NDIM*NDIM,M3(1,1,N),1,B(1,1,N),1)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,N-1),NDIM,
     +       F(1,1),NDIM,CONE,B(1,1,N),NDIM)
        
c
c ---> C(N) = - M3(N-1)*D(N-1)*C(N-1) + M1(N)
c
        CALL ZCOPY(NDIM*NDIM,M1(1,1,N),1,C(1,1,N),1)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +       E(1,1),NDIM,CONE,C(1,1,N),NDIM)
c
c ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
c
        
      CALL ZCOPY(NDIM*NDIM,M2(1,1,N),1,E(1,1),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +     F(1,1),NDIM,CONE,E(1,1),NDIM)

      IF (SAVD.GT.0)
     +     CALL ZCOPY(NDIM*NDIM,E(1,1),1,M2(1,1,N),1)

      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,N),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,N),NDIM,INFO)

c ------------------------------------------------------------------------
c
c ---> N = NL
c
 970  N = NL
c
c ---> D(NL) = (A - B(NL-1)*D(NL-1)*C(NL-1))^-1
c
c
c ---> E = D(NL-1) * C(NL-1)
c
c      CALL ZGEMUL(D(1,1,NL-1),NDIM,'N',C(1,1,NL-1),NDIM,'N',
c     +     E(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL-1),NDIM,
     +     C(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)
c
c ---> A = A - B(NL-1) * E
c
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,B(1,1,NL-1),NDIM,
     +     E(1,1),NDIM,CONE,A(1,1),NDIM)
c
c ---> D(NL) = (A)^-1
c
      IF (SAVD.GT.0)
     +     CALL ZCOPY(NDIM*NDIM,A(1,1),1,M2(1,1,NL),1)

      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,NL),1)
      CALL ZGETRF(NDIM,NDIM,A(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,A(1,1),NDIM,IPVT,D(1,1,NL),NDIM,INFO)

 980  CONTINUE                      ! jump label for NL=1
c ------------------------------------------------------------------------
c
c --->  END OF FACTORIZATION
c
c ------------------------------------------------------------------------
c
c ---> in case of iopt=1 RETURN
c                      2 save C(n) and RETURN
c
      IF (SAVC.GT.0 .AND. NL.GT.1) THEN
        DO 110 N = 1,NL-1
          CALL ZCOPY(NDIM*NDIM,C(1,1,N),1,M3(1,1,N),1)
 110    END DO
      END IF
c
      IF (SAVC.GT.0 .OR. SAVD.GT.0) RETURN
c ------------------------------------------------------------------------
c
c ---> DETERMINATION OF THE DIAGONAL ELEMENTS
c
c      of the invers of matrix M ( stored in M2(1,1,N) )
c
c
c ---> N = NL
c
      N = NL
c
c ---> M2(N) = D(NL)
c
      CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,M2(1,1,NL),1)

c ------------------------------------------------------------------------
c
c ---> N = NL -1
c
      IF (NL.EQ.1) GOTO 990

      N = NL -1
c
c ---> M2(NL-1)= D(NL-1) + D(NL-1)*C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
c
c
c ---> E = 1 + C(NL-1)*D(NL)*B(NL-1)*D(NL-1)
c
c      CALL ZGEMUL(B(1,1,NL-1),NDIM,'N',D(1,1,NL-1),NDIM,'N',
c     +     E(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,NL-1),NDIM,
     +     D(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)
c      CALL ZGEMUL(D(1,1,NL),NDIM,'N',E(1,1),NDIM,'N',
c     +     F(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL),NDIM,
     +     E(1,1),NDIM,CZERO,F(1,1),NDIM)
c      CALL ZGEMUL(C(1,1,NL-1),NDIM,'N',F(1,1),NDIM,'N',
c     +     E(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,NL-1),NDIM,
     +     F(1,1),NDIM,CZERO,E(1,1),NDIM)
      DO 30 LM = 1,NDIM
        E(LM,LM) = CONE + E(LM,LM)
 30   CONTINUE
c
c ---> M2(NL-1) = D(NL-1) * E
c
c      CALL ZGEMUL(D(1,1,NL-1),NDIM,'N',E(1,1),NDIM,'N',
c     +     M2(1,1,NL-1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL-1),NDIM,
     +     E(1,1),NDIM,CZERO,M2(1,1,NL-1),NDIM)

      IF (NL.EQ.2 .AND. SAVY.EQ.0) GOTO 990

c
c ---> W = - D(NL-1)*C(NL-1)*D(NL)
c
c      CALL ZGEMUL(C(1,1,NL-1),NDIM,'N',D(1,1,NL),NDIM,'N',
c     +     E(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,NL-1),NDIM,
     +     D(1,1,NL),NDIM,CZERO,E(1,1),NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,D(1,1,NL-1),NDIM,
     +     E(1,1),NDIM,CZERO,W(1,1),NDIM)
c
c ---> V = - D(NL)*B(NL-1)*D(NL-1)
c
c      CALL ZGEMUL(B(1,1,NL-1),NDIM,'N',D(1,1,NL-1),NDIM,'N',
c     +     E(1,1),NDIM,NDIM,NDIM,NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,NL-1),NDIM,
     +     D(1,1,NL-1),NDIM,CZERO,E(1,1),NDIM)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,D(1,1,NL),NDIM,
     +     E(1,1),NDIM,CZERO,V(1,1),NDIM)
c
c ---> save elements for cluster GF
c
      CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,V1(1,1),1)
      CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,W1(1,1),1)
c
c ---> save elements for cluster gf determination
c
      IF (SAVY.NE.0) THEN
        CALL ZCOPY(NDIM*NDIM,W(1,1),1,M1(1,1,NL-1),1)
        CALL ZCOPY(NDIM*NDIM,V(1,1),1,M4(1,1,NL-1),1)
      END IF

      IF (NL.EQ.2) GOTO 990
c ------------------------------------------------------------------------
c
c ---> NL-1 > N >=1
c
      DO 40 N = (NL-2),1,(-1)

c ------------------------------------------------------------------------
        IF (SAVY.NE.0) THEN
c
c --->    save elements for cluster gf determination
c
          CALL ZCOPY(NDIM*NDIM,V1(1,1),1,V2(1,1),1)
          CALL ZCOPY(NDIM*NDIM,W1(1,1),1,W2(1,1),1)
          CALL ZCOPY(NDIM*NDIM,V(1,1),1, V1(1,1),1)
          CALL ZCOPY(NDIM*NDIM,W(1,1),1, W1(1,1),1)
        END IF
c ------------------------------------------------------------------------
c
c ---> W = D(N) * [-M1(N)* W - C(N)*D(NL)]
c
c        CALL ZCOPY(NDIM*NDIM,C(1,1,N) ,1,E(1,1),1) ! E = C(1,1,N)
c        CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,F(1,1),1) ! F = D(1,1,NL)
c        CALL ZGEMUL(C(1,1,N),NDIM,'N',D(1,1,NL),NDIM,'N',
c     +       E(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,N),NDIM,
     +       D(1,1,NL),NDIM,CZERO,G(1,1),NDIM)

c        CALL ZCOPY(NDIM*NDIM,M1(1,1,N),1,E(1,1),1) ! E = M1(1,1,N)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M1(1,1,N),NDIM,
     +       W(1,1),NDIM,-CONE,G(1,1),NDIM)

c        CALL ZCOPY(NDIM*NDIM,D(1,1,N),1,E(1,1),1) ! E = D(1,1,N)
c        CALL ZGEMUL(D(1,1,N),NDIM,'N',E(1,1),NDIM,'N',
c     +       W(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N),NDIM,
     +       G(1,1),NDIM,CZERO,W(1,1),NDIM)
c     
c ---> M(2,N) = D(N) + D(N)*[ M1(N)*M2(N+1) + C(N)*V ]*M3(N)*D(N) -
c               W*B(N)*D(N)
c
c     
c ---> F = C(N)*V
c
c        CALL ZGEMUL(C(1,1,N),NDIM,'N',V(1,1),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,N),NDIM,
     +       V(1,1),NDIM,CZERO,F(1,1),NDIM)
c
c        CALL ZCOPY(NDIM*NDIM,M1(1,1,N),1,E(1,1),1) ! E = M1(1,1,N)
c
c ---> F = [M1(N)*M2(N+1) + C(N)*V]
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M1(1,1,N),NDIM,
     +       M2(1,1,N+1),NDIM,CONE,F(1,1),NDIM)
c ------------------------------------------------------------------------
        IF (SAVY.NE.0) THEN
c
c --->    Y = -D(N)*F
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,D(1,1,N),NDIM,
     +         F(1,1),NDIM,CZERO,Y(1,1),NDIM)
        END IF
c ------------------------------------------------------------------------
c
c ---> F = 1 + [ M1(N)*M2(N+1) + C(N)*V ]*M3(N)*D(N)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,F(1,1),NDIM,
     +       M3(1,1,N),NDIM,CZERO,G(1,1),NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,G(1,1),NDIM,
     +       D(1,1,N),NDIM,CZERO,F(1,1),NDIM)
        DO 50 LM = 1,NDIM
          F(LM,LM) = CONE + F(LM,LM)
 50     CONTINUE
c
c ---> M2(N) = W*B(N)*D(N)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,B(1,1,N),NDIM,
     +       D(1,1,N),NDIM,CZERO,E(1,1),NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,W(1,1),NDIM,
     +       E(1,1),NDIM,CZERO,M2(1,1,N),NDIM)
c
c ---> M2(N)= D(N)*F - M2(N)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N),NDIM,
     +       F(1,1),NDIM,-CONE,M2(1,1,N),NDIM)
c
c ---> V = [ -V*M3(N) - D(NL)*B(N) ]*D(N)
c
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,NL),NDIM,
     +       B(1,1,N),NDIM,CZERO,E(1,1),NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,V(1,1),NDIM,
     +       M3(1,1,N),NDIM,-CONE,E(1,1),NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,E(1,1),NDIM,
     +       D(1,1,N),NDIM,CZERO,V(1,1),NDIM)

c ------------------------------------------------------------------------
        IF (SAVY.NE.0) THEN
c
c --->    save elements for cluster GF determination
c
c
c --->    YY = D(N) * [-M1(N)* M1(N+1) - C(N)*V2]
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,C(1,1,N),NDIM,
     +         V2(1,1),NDIM,CZERO,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M1(1,1,N),NDIM,
     +         M1(1,1,N+1),NDIM,-CONE,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N),NDIM,
     +         G(1,1),NDIM,CZERO,YY(1,1),NDIM)

c
c --->    X = [-M2(N+1)* M3(N) - W1* B(N)] * D(N)
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,W1(1,1),NDIM,
     +         B(1,1,N),NDIM,CZERO,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M2(1,1,N+1),NDIM,
     +         M3(1,1,N),NDIM,-CONE,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,G(1,1),NDIM,
     +         D(1,1,N),NDIM,CZERO,X(1,1),NDIM)
c
c --->    XX = [-M4(N+1)* M3(N) - W2* B(N)] * D(N)
c
          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,W2(1,1),NDIM,
     +         B(1,1,N),NDIM,CZERO,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M4(1,1,N+1),NDIM,
     +         M3(1,1,N),NDIM,-CONE,G(1,1),NDIM)

          CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,G(1,1),NDIM,
     +         D(1,1,N),NDIM,CZERO,XX(1,1),NDIM)

          CALL ZCOPY(NDIM*NDIM,Y(1,1) ,1,M1(1,1,N),1) ! M1(N) = Y
          CALL ZCOPY(NDIM*NDIM,YY(1,1),1,M3(1,1,N),1) ! M3(N) = YY
          CALL ZCOPY(NDIM*NDIM,X(1,1) ,1,M4(1,1,N),1) ! M4(N) = X
          CALL ZCOPY(NDIM*NDIM,XX(1,1),1,M5(1,1,N),1) ! M5(N) = XX

        END IF                      ! (SAVY.NE.0)
c ------------------------------------------------------------------------
C
 40   CONTINUE                      ! N = (NL-2),1,(-1)

 990  CONTINUE                      ! jump mark for small NL
c
c ------------------------------------------------------------------------
C     END OF SUPERCELL ALGORITHM
c ************************************************************************
c
      ELSE                          ! (SLAB.EQ.0)
c
c ************************************************************************
c
c ---> ALGORITHM FOR SLAB GEOMETRY
c
c ------------------------------------------------------------------------
C     N = 1
C
      CALL ZCOPY(NDIM*NDIM,M2(1,1,1),1,E(1,1),1)
      CALL ZCOPY(NDIM*NDIM,CUNIT,1,D(1,1,1),1)
      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,1),NDIM,INFO)

      IF (NL.EQ.1) GOTO 90
c ------------------------------------------------------------------------
C     2 <= N <= NL
C
      DO 60 N = 2,NL
c
c ---> F = D(N-1) * M1(N-1)
c
c        CALL ZGEMUL(D(1,1,N-1),NDIM,'N',M1(1,1,N-1),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N-1),NDIM,
     +       M1(1,1,N-1),NDIM,CZERO,F(1,1),NDIM)
        
c
c ---> D(N) = [ M2(N) - M3(N-1)*D(N-1)*M1(N-1) ]^-1
c
        
      CALL ZCOPY(NDIM*NDIM,M2(1,1,N),1,E(1,1),1)
      CALL ZGEMM('N','N',NDIM,NDIM,NDIM,-CONE,M3(1,1,N-1),NDIM,
     +     F(1,1),NDIM,CONE,E(1,1),NDIM)

      IF (SAVD.GT.0)
     +     CALL ZCOPY(NDIM*NDIM,E(1,1),1,M2(1,1,N),1)

      CALL ZGETRF(NDIM,NDIM,E(1,1),NDIM,IPVT,INFO)
      CALL ZCOPY(NDIM*NDIM,CUNIT(1,1),1,D(1,1,N),1)
      CALL ZGETRS('N',NDIM,NDIM,E(1,1),NDIM,IPVT,D(1,1,N),NDIM,INFO)
        
 60   CONTINUE

 90   CONTINUE                      ! jump label for nl=1
c
c ------------------------------------------------------------------------
c  end of factorization
c ------------------------------------------------------------------------
c
c ---> in case of iopt=11 RETURN
c                      12 set m3(n.lt.nl-1)=0 and
c                         m3(nl-1)=m1(nl-1) and RETURN
c
      IF (SAVC.GT.0 .AND. NL.GT.1) THEN
        CALL CINIT(NDIM*NDIM*(NL-1),M3(1,1,1))
        CALL ZCOPY(NDIM*NDIM,M1(1,1,NL-1),1,M3(1,1,NL-1),1)
      END IF
c
      IF (SAVC.GT.0 .OR. SAVD.GT.0) RETURN
c ------------------------------------------------------------------------
c
c ---> determination of the diagonal elements of the invers of M
c      storing in M2(1,1,n)
c
c ------------------------------------------------------------------------
c
c ---> N = NL
c
      N = NL
c
c ---> M2(N) = D(NL)
c
      CALL ZCOPY(NDIM*NDIM,D(1,1,NL),1,M2(1,1,NL),1)

      IF (NL.EQ.1) GOTO 100
c ------------------------------------------------------------------------
c
      DO 70 N = NL-1, 1 ,(-1)
c     
c ---> M(2,N) = D(N) + D(N)*M1(N)*M2(N+1)*M3(N)*D(N) 
c
c     
c ---> F = 1 + M1(N)*M2(N+1)*M3(N)*D(N)
c
c        CALL ZGEMUL(M3(1,1,N),NDIM,'N',D(1,1,N),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M3(1,1,N),NDIM,
     +       D(1,1,N),NDIM,CZERO,F(1,1),NDIM)
c        CALL ZGEMUL(M2(1,1,N+1),NDIM,'N',F(1,1),NDIM,'N',
c     +       E(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M2(1,1,N+1),NDIM,
     +       F(1,1),NDIM,CZERO,E(1,1),NDIM)
c        CALL ZGEMUL(M1(1,1,N),NDIM,'N',E(1,1),NDIM,'N',
c     +       F(1,1),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,M1(1,1,N),NDIM,
     +       E(1,1),NDIM,CZERO,F(1,1),NDIM)
        DO 80 LM = 1,NDIM
          F(LM,LM) = CONE + F(LM,LM)
 80     CONTINUE
c     
c ---> M2(N) = D(N)*F 
c
c        CALL ZGEMUL(D(1,1,N),NDIM,'N',F(1,1),NDIM,'N',
c     +       M2(1,1,N),NDIM,NDIM,NDIM,NDIM)
        CALL ZGEMM('N','N',NDIM,NDIM,NDIM,CONE,D(1,1,N),NDIM,
     +       F(1,1),NDIM,CZERO,M2(1,1,N),NDIM)

 70   CONTINUE

 100  CONTINUE                      ! jump mark for nl=1

      END IF                        ! (SLAB.EQ.0)

      RETURN
      
      END
