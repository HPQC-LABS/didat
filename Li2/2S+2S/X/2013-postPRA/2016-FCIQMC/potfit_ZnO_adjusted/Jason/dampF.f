c***********************************************************************
      SUBROUTINE dampF(r,RHOab,NCMM,m,MMAX,i,IDF,WD,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c---------------------- RJL Version of 19 August 2009 ------------------
c-----------------------------------------------------------------------
c                 Upon Input
c* RHOab - Scaling coeffcient used for systems != H_2
c          RHOab= 2*(RHOa*RHOb)/(RHOa+RHOb) where RHOa = (I_p^A/I_p^H)^0.66
c                where I_p^A is the ionization potential of atom A
c                and I_p^H is the ionization potential of atomic hydrogen
c* r - the radial distance at which the damping function is being evaluated
c* MMAX - The exponent of the highest order long range term used in the 
c         analysis
c* l= 2 yields the Tang-Toennies damping functions [JCP 80,3726(1984)]
c         defined such that  Dm(r)/r^m \propto{r^1} --> 0  as r --> 0
c* l= 1 yields Douketis-Scoles damping function defined such that
c                Dm(r)/r^m --> {const.} \propto  r^0  as r --> 0
c* l=  0 yields modified Douketis-Scoles damping functions defined such that
c                Dm(r)/r^m \propto  r^{-1}  as r --> 0
c    all assuming  R  in units Angst. 
c    Exit if  IDF is not one of these values
c  WD  = Designates if WD.GT.0 the second derivative is calculated 
 
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
     1    ,DMP(MMAX),DMPP(MMAX),SM(-3:MMAX+10),mfact,TK,DM(MMAX)
     2  , RHOab,RHOr,bpm(6), cpirm(6) ,coef,coef2
c------------------------------------------------------------------------
c  The following values for the numerical factors used in both TT and DS
c  were  normalized to the Hydrogen data presented
c  by Kreek and Meath in J. Chem. Phys. Volume 50, Page 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c------------------------------------------------------------------------
       DATA  tnf  / 3.13/
       DATA  bnf / 3.69,3.95/
       DATA cnf / 0.4  , .025 /
       DATA FIRST/ 1/
       SAVE FIRST, bpm, cpirm
       
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
      l = IDF -1
      IF (IDF.EQ.2) THEN
       
      br= RHOab*r* tnf
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

      IF (FIRST.EQ.1) THEN
         DO K=1 , NCMM 
              bpm(k) = bnf(l)*RHOab/M(K,i)
              cpirm(k)= cnf(l) *RHOab**2/DSQRT(DFLOAT(M(K,i)))
         END DO
         FIRST = 0 
      END IF
      DO N = 1 , NCMM
         K = (2*M(N,i)+l)/2
         XP= DEXP (-bpm(N)*r-cpirm(N)*r**2)
              
         
      Dm(M(N,i)) = (1- XP)**(K)
        IF (WD.GT.0) THEN
         Coef2 =  2*cpirm(N)
         Coef= K*(bpm(N)+Coef2*r)
         Dmp(M(N,i)) = K*Coef*XP*(1-XP)**(K-1)
c        Dmpp(M(N,i)) =  Coef2*XP*(1-XP)**(K-1)
c    1                - Coef*Coef2*XP*(1-XP)**(K-1) 
c    2                + K*Coef*Coef2*XP*(1-XP)**(K-2) 
        END IF
      ENDDO   
      END IF  
      RETURN
  600 FORMAT(/,'ERROR!! Input variable IDF=',i3,'  does not give a valid
     1 damping function.',/)
  602 FORMAT(/' *** dampF  dies because   MMAX=',i3,'   l=',i3)
  604 FORMAT( /,'ERROR!!! Input variable RHOa=', F7.4,' or RHOb= '
     1 , F7.4, '  not valid') 
      END
              
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
