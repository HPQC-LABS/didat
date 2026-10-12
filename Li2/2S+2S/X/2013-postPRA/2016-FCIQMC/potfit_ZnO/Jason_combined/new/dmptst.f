      REAL*8 r,b,DM(20),DMP(20),DMPP(20),QM(20),QMP(20),QMPP(20)
      INTEGER i,m,MMAX,IDF

      MMAX= 10
      IDF= 1
    2 r= 1.0d0
      b= 2.0d0
      WRITE(6,600)  IDF
  600 FORMAT(/'   IDF=',i3)
      DO  i= 1,40
          CALL dampF(r,b,MMAX,IDF,DM,DMP,DMPP)
          WRITE(6,601) r,b*r,  DM(6),DM(6)/r**6 ,DM(8),DM(8)/r**8,
     1                               DM(10),DM(10)/r**10
  601 Format(f9.6,f9.6,3(1PD14.6,D16.8))

          IF(r.GT.0.10000001d0) then
              r= r-0.1d0
              goto 10
              endif
          IF(r.GT.0.01000001d0) then
              r= r-0.01d0
              goto 10
              endif
          IF(r.GT.0.00100001d0) then
              r= r-0.001d0
              goto 10
              endif
          IF(r.GT.0.00010001d0) then
              r= r-0.0001d0
              goto 10
              endif
   10     continue
          ENDDO
      IDF= IDF-1
      GOTO 2

      stop
      end


c***********************************************************************
      SUBROUTINE dampF(r,b,MMAX,IDF,DM,DMP,DMPP)
c** Subroutine to generate values 'Dm' and its first `Dmp' and second
c   'Dmpp' derivatives w.r.t. R of the chosen version of the incomplete
c    gamma function damping function, for  m= 1 to MMAX.
c***********************************************************************
c                 Upon Input
c* b - the value of the damping coeffient (method of calculation described
c      later in the subroutine)
c* r - the radial distance at which the damping function is being evaluated
c* MMAX - The exponent of the highest order long range term used in the 
c         analysis
c* IDF= 1 yields the Tang-Toennies damping functions [JCP 80,3726(1984)]
c         defined such that  Dm(r)/r^m \propto{r^1} --> 0  as r --> 0
c* IDF= 0 yields damping functions defined such that
c                Dm(r)/r^m --> {const.} \propto  r^0  as r --> 0
c* IDF= -1 yields damping functions defined such that
c                Dm(r)/r^m \propto  r^{-1}  as r --> 0
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
      INTEGER MMAX,k,m,n,IDF
      REAL*8 b,r,br,XP,SSm,SSm1,SSm2,mfact,TK,DM(MMAX),DMP(MMAX),
     1       DMPP(MMAX),SM(-3:MMAX+10) 
c    2  QM(MMAX),QMP(MMAX),QMPP(MMAX)
c------------------------------------------------------------------------
      IF((IDF.LT.-1).OR.(IDF.GT.1)) THEN
          WRITE(6,600)  IDF
          STOP
          ENDIF
      IF(MMAX.LT.(1-IDF)) THEN
          WRITE(6,602) 
  602 FORMAT(/' *** dampF  dies because   MMAX=',i3,'   IDF=',i3)
          STOP
          ENDIF
c------------------------------------------------------------------------
c  The value of b used has been normalized to the Hydrogen data presented
c  by Kreek and Meath in J. Chem. Phys. Volume 50, Page 2289 (1969).
c  The ratio has been chosen such that  b= FACTOR*(I_p^X / I_p^H)^{2/3}
c  for the homoatomic diatomic species X_2, where I_p^A is the ionization
c  potential of atom A. For the heteroatomic diatomic species XY, 
c       b= 2*(b^X * b^Y)/(b^X + b^Y)
c  FACTOR(IDF=0)= 2.78;  FACTOR(IDF=1)= X.XX;  FACTOR(IDF=-1)= X.XX;  
c------------------------------------------------------------------------
      br= b*r
      XP= DEXP(-br)
      SM(-3)= 0.d0
      SM(-2)= 0.d0
      SM(-1)= 0.d0
      SM(0)=  1.d0
      TK= 1.d0
      IF(br.GT.0.5d0) THEN
          DO  k= 1,MMAX-1+IDF
              TK= TK*br/DFLOAT(k)
              SM(k)= SM(k-1)+ TK
              ENDDO
          DO m= 1, MMAX
              DM(m)= 1.d0 - XP*SM(m-1+IDF)
              DMP(m)= b*XP*(SM(m-1+IDF) - SM(m-2+IDF))
              DMPP(m)= -b*b*XP*(SM(m-1+IDF) - 2.d0*SM(m-2+IDF)
     1                                                  + SM(m-3+IDF))
              ENDDO
c-----------------------------------------------------------------------
c  The above section handles the calculation of the value of the damping
c  function for most values of r. However, as stated above, at very small
c  r, the damping function becomes unstable due to numerical noise. Thus,
c  if the value of the damping function falls below a certain value it
c  is re-evaluated by the different functional definition below.
c-----------------------------------------------------------------------
        ELSE
          DO  k= 1,MMAX+4+IDF
c... NOTE that here SM(k) is the k'th term  (b*r)^k/k!  [not a sum]
              SM(k)= SM(k-1)*br/DFLOAT(k)
              ENDDO
          DO  m= 1, MMAX
              k= m+IDF
              DM(m)= XP*(SM(k)+ SM(k+1)+ SM(k+2)+ SM(k+3)+ SM(k+4))
              DMP(m)= b*XP*SM(k-1)
              DMPP(m)= b*b*XP*(SM(k-2) - SM(k-1))
              ENDDO
        ENDIF
      RETURN
  600 FORMAT(/,'ERROR!! Input variable IDF=',i3,'  does not give a valid
     1 damping function.',/)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

