cThis subroutine is NOT IDENTICAL to the one found used in betaFIT,
c in order to selectively calculate Dm(r) when using the DS form for
c only m values used, the value of ISTATE which only appears in potfit 

c***********************************************************************
      SUBROUTINE dampF(r,RHOab,NCMM,m,MMAX,i,IDF,WD,DM,DMP,DMPP)

c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 19 August 2009 ------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* RHOab - Effective scaling parameter for systems != H2, calculated
clater by RHOab= 2*RHOa*RHOb/(RHOa+RHOb) Each RHOa is defined as 
c (I_p^A/I_p^H)^.66 where I_p^A is the ionization potential of atom A and I_p^H
c is the ionization potential of atomic hydrogen. If A=B RHOab = RHOa= RHOb    
c* r - the radial distance at which the damping function is being evaluated
c* MMAX - The exponent of the highest order long range term used in the 
c         analysis
c* IDF= 2 yields the Tang-Toennies damping functions [JCP 80,3726(1984)]
c         defined such that  Dm(r)/r^m \propto{r^1} --> 0  as r --> 0
c* IDF= 1 yields the Douketis-Scoles damping functions defined such that
c                Dm(r)/r^m --> {const.}   as r --> 0
c* IDF= 0 yields Douketis-Scoles damping functions defined such that (l=-1)
c                Dm(r)/r^m \propto  r^{-1/2}  as r --> 0
c    all assuming  R  in units Angst. 
c    Exit if  IDF is not one of these values
c----------------------------------------
c                 Upon Output
c  DM(i) - The value of the damping function for the long range term 
c          C_i/r^i
c  DMP(i) - The value of the first derivative of the damping function for 
c           the long range term C_i/r^i
c  DMPP(i) - The value of the second derivative of the damping function 
c            for the long range term C_i/r^i
c-----------------------------------------------------------------------
      INTEGER MMAX,k,m(6,4),n,IDF,l,Mi,NCMM, WD ,FIRST,i
      REAL*8 b,tnf,cnf(-1:0),bnf(-1:0), r,br,XP,SSm,SSm1,SSm2
     1    ,DMP(MMAX),IRM(20),DMPP(MMAX),SM(-3:MMAX+10),mfact,TK,DM(MMAX)
     2  , RHOab,RHOr, ba,bb
c The following values were normalized to the Hydrogen data presented
c  by Kreek and Meath in J. Chem. Phys. Volume 50, Page 2289 (1969).

 
c value for the numerical factor used in the TT form
       DATA  tnf  / 3.13/
c value for the numberical factor b used in the DS form for l=-1,0, respectively
       DATA  bnf / 3.66,3.95/
c value for numerical factor c in the DS form, for l= -1, 0 respectively
       DATA cnf / 0.03,0.025 /
       DATA FIRST/ 1/
       SAVE FIRST, IRM
       
c    2  QM(MMAX),QMP(MMAX),QMPP(MMAX)
c------------------------------------------------------------------------
      IF(IDF.GT.2) THEN
          WRITE(6,600)  IDF
          STOP
          ENDIF
      IF(MMAX.LT.(1-IDF)) THEN
          WRITE(6,602)  MMAX, IDF
          STOP
         END IF
      IF (RHOab.LE.0) THEN
          WRITE (6, 604) RHOab
          STOP
 
          ENDIF
         
c This section  can be removed
c------------------------------------------------------------------------
c  The value of b used has been normalized to the Hydrogen data presented
c  by Kreek and Meath in J. Chem. Phys. Volume 50, Page 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c  potential of atom A. For the heteroatomic diatomic species XY, 
c       b= 2*(b^X * b^Y)/(b^X + b^Y)
c  FACTOR(IDF=-1)= 2.44;  FACTOR(IDF=0)= 2.78;  FACTOR(IDF=+1)= 3.13;  
c------------------------------------------------------------------------

      IF (IDF.EQ.2) THEN
      l = 1 
      b = RHOab*tnf
      br= b*r
      XP= DEXP(-br)
      SM(-3)= 0.d0
      SM(-2)= 0.d0
      SM(-1)= 0.d0
      SM(0)=  1.d0
      TK= 1.d0
      IF(br.GT.0.5d0) THEN
          DO  k= 1,MMAX-1+l
              TK= TK*br/DFLOAT(k)
              SM(k)= SM(k-1)+ TK
              ENDDO
          DO N= 1, MMAX
              DM(N)= 1.d0 - XP*SM(N-1+l)
              DMP(N)= b*XP*(SM(N-1+l) - SM(N-2+l))
              DMPP(N)= -b*b*XP*(SM(N-1+l) - 2.d0*SM(N-2+l)
     1                                                  + SM(N-3+l))
              ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r. However, as stated above, at very small
c  r, the damping function becomes unstable due to numerical noise. Thus,
c  if the value of the damping function falls below a certain value it
c  is re-evaluated by the different functional definition below.
c-----------------------------------------------------------------------
        ELSE
          DO  k= 1,MMAX+4+l
c... NOTE that here SM(k) is the k'th term  (b*r)^k/k!  [not a sum]
              SM(k)= SM(k-1)*br/DFLOAT(k)
              ENDDO
          DO  N= 1, MMAX
              k= N+l
              DM(N)= XP*(SM(k)+ SM(k+1)+ SM(k+2)+ SM(k+3)+ SM(k+4))
              DMP(N)= b*XP*SM(k-1)
              DMPP(N)= b*b*XP*(SM(k-2) - SM(k-1))
              ENDDO
        ENDIF
      ELSEIF ((IDF.EQ.0) .or. (IDF.EQ.1)) THEN
       l  =  IDF -1
      IF (FIRST.EQ.1) THEN
         DO K=1 , 20
              IRM(K) = 1/DSQRT(DFLOAT(K))
         END DO
         FIRST = 0 
      END IF
      b = bnf(l) * RHOab
      br = b*r
      DO N = 1 , NCMM
         K = (2*M(N,i)+l)/2
         XP= DEXP (-br/M(N,i)-cnf(l)*(br)**2*IRM(M(N,i)))
         
      Dm(M(N,i)) = (1- XP)**(K)
         
        IF (WD.EQ.1) THEN
         Dmp =  K*(bnf(l)*RHOr/M(N,i)-2*cnf(l)*b**2*r*IRM(M(N,i)))
     1          *XP * (1-XP)**(K-1)
        END IF
        End DO
      END IF  
      RETURN
  600 FORMAT(/,'ERROR!! Input variable IDF=',i3,'  does not give a valid
     1 damping function.',/)
  602 FORMAT(/' *** dampF  dies because   MMAX=',i3,'   l=',i3)
  604 FORMAT( /,'ERROR!!! Input variable RHOab= '
     1 , F7.4, '  not valid') 
      END
              
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
