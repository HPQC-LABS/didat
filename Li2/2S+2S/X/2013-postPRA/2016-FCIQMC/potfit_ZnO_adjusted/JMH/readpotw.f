c**********************************************************************
      SUBROUTINE READPOT(ISTATE)
c**********************************************************************
c** This subroutine will generate one of five possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program. These four families are:
c     - The Generalized Morse Oscillator
c     - The Expanded Morse Oscillator
c     - The Modified Morse Oscillator
c     - The Modified Lennard-Jones Oscillator
c     - The Generalized Potential Energy Function
c     - The Double-Exponential Long-Range Potential
c** This subroutine can also calculate the trial parameters for each of
c   the four families of analytical potential functions from the low
c   order molecular constants: Be, We, WeXe, ae.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 2000  by  Jenning Seto       (ver of 12/10/1997) +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NPOT    is a flag to determine if a new potential is to be read in.
c            NPOT = 0 : Use old potential paramters.
c            NPOT = 1 : Read in new potential paramters.
c
c    ISTATE  is the electronic state being fitted to.
c
      INCLUDE 'arrsizes.h'

      INTEGER ISTATE, I

      REAL*8 RSR
c
c** Type statements & common block for data
c
      REAL*8 FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)

      INTEGER COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)

      CHARACTER*2 NAME(2)

      CHARACTER*2 SLABL(-3:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
     3  ,nBetaIn(NSTATEMX),nPB(NSTATEMX),inAndOut(NSTATEMX)
     4  ,nPq(NSTATEMX),nPu(NSTATEMX),nPL(NSTATEMX)
 
      COMMON /nPwBLK/nPB,nPq,nPu,nPL
      COMMON /nBetaInBLK/nBetaIn
      COMMON /inAndOutBLK/inAndOut

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),Dn(NSTATEMX,numCnMax),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

c** For Expanded Morse Long-Range Potential
      INTEGER numCn(NSTATEMX),nLR(NSTATEMX,numCnMax),Df
      REAL*8 dampBeta(NSTATEMX),Cn(NSTATEMX,numCnMax),Demo(NSTATEMX),
     &  DeOffset,ReOffset,Dnp,Dnpp,ReBar
      COMMON /DELRBLK/numCn,nLR,dampBeta,Cn,ReBar,Df
      COMMON /OffSetBLK/Demo,DeOffset,ReOffset

c** The possible repulsive chemical bonding
      REAL*8 C(NSTATEMX),DVDC(NPNTMX,NSTATEMX)
      INTEGER ifC(NSTATEMX)
      COMMON /repulsiveBLK/C,ifC,DVDC

      REAL*8 cRp(NSTATEMX)
      COMMON /cRpBLK/cRp

c** For Modified Lennard-Jones-Gaussian Potential
      INTEGER ifAmp(NSTATEMX),ifR0(NSTATEMX),ifFWHM(NSTATEMX)
      REAL*8 Amp(NSTATEMX),R0(NSTATEMX),FWHM(NSTATEMX)
      COMMON /MLJG/ifAmp,ifR0,ifFWHM,Amp,R0,FWHM
c
c** Common block that will hold all of the information for the potential
c   and the partial derivative arrays.
c
c  BETA     holds the actual exponent parameters used for the potential
c
c  NDATPT   number of data points used for each potential
c
c  VPOT	holds the potential arrays for each state
c
c** Born-Oppenheimer Correction Parameters.
c
 
      INTEGER PoEpZ(NSTATEMX)

      COMMON /PoEpZBLK/PoEpZ
c
c** Born-Oppenheimer Breakdown Parameters.
c
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NQA(NSTATEMX),NQB(NSTATEMX),
     1  IFUA(0:NABOMX,NSTATEMX),IFUB(0:NABOMX,NSTATEMX),
     2  IFQA(0:NRNAMX,NSTATEMX),IFQB(0:NRNAMX,NSTATEMX)

      REAL*8 UA(0:NABOMX,NSTATEMX),UB(0:NABOMX,NSTATEMX),
     1  QA(0:NRNAMX,NSTATEMX),QB(0:NRNAMX,NSTATEMX)

      COMMON /BPBLK/NUA,NUB,NQA,NQB,IFUA,IFUB,IFQA,IFQB,UA,UB,QA,QB
c
c** Born-Oppenheimer Breakdown Functions.
c     
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1  QAR(NPNTMX,NSTATEMX),QBR(NPNTMX,NSTATEMX),
     2  ZMUA(NISTPMX),ZMUB(NISTPMX),ZMQA(NISTPMX),ZMQB(NISTPMX)

      COMMON /BFBLK/UAR,UBR,QAR,QBR,ZMUA,ZMUB,ZMQA,ZMQB

**
** Lambda doubling function 
**
      INTEGER nLamdCoef(NSTATEMX),ifLamdCoef(0:NRNAMX,NSTATEMX)
      REAL*8 lamdCoef(0:NRNAMX,NSTATEMX),lambdaR(NPNTMX,NSTATEMX)
      
      COMMON /lambdaBLK/nLamdCoef,ifLamdCoef,lamdCoef,lambdaR
c-----------------------------------------------------------------------
c** Now to choose potential model and desired Born-Oppenheimer
c   correction expansion parameter.
c
c  PSEL(s)  choses the type of analytical potential to be fitted to:
c         = 0     : Do forward calculation
c         = +/- 1 : Use the Generalized Morse Potential.
c         = +/- 2 : Use the Expanded Morse Potential.
c         = +/- 3 : Use the Modified Morse Potential.
c         = +/- 4 : Use the Modified Lennard-Jones Potential.
c         = +/- 5 : Use the Generalized Potential Energy Function.
c         = +/- 6 : Use the Expanded Morse Long-Range Potential.
c         = +/- 7 : Use the Modified Lennard-Jones-Gaussian Potential.
c         = +/- 8 : Use the Total Expanded Morse Potential.
c         = +/- 9 : Total Expanded Morse Long-Range Potential
c                   (trail version of DELR).
c         = +/- 10: Use the Double-Exponential Long-Range Potential.
c           if NBETA < 0 then read in parmaters are in fact the
c           molecular constants for that molecule and the trial
c           exponential paramters must be calculated from these
c           values. The order that they are read in are:
c           De, re, We (Y10), WeXe (-Y20), Be (Y01), ae (-Y11).
c
c  PoEpZ(s) specifies type of expansion used for Born-Oppenheimer
c           correction, lambda doubling and beta functions.
c            = 0 : use half of OT expansion (R-RE)/(R+RE)
c            =-1 : use Coxon expansion (R-RE)
c            > 0 : use new variable (R^m-Re^m)/(R^m+Re^m) for BOB
c            = 2 : use new variable (R^m-Re^m)/(R^m+Re^m) for beta
c            = 3 : use new exponent (replace betaR*z) for MLJ
c  inAndOut(s) indicates the #terms of the beta function used in inner wall
c            = 1: beta function uses N_S terms at R < Re 
c                 beta function uses N_L terms at R >= Re  (N_S < N_L)
c            = 0: beta function of repulsive exponential uses N_S terms
c                 beta function of attractive exponential uses N_L terms 
c                 e.g. 'inAndOut(s) = 0, N_S = N_L = 7' is the same as
c                      'inAndOut(s) = 1, N_S = N_L = 7'
c                  but 'inAndOut(s) = 0, N_S = 3, N_L = 7' differ from
c                      'inAndOut(s) = 1, N_S = 3, N_L = 7'  
c  nPB(s)   the value of power m for beta function   
c  nPq(s)   the value of power m for q(r) function   
c  nPu(s)   the value of power m for u(r) function   
c  nPL(s)   the value of power m for Lambda(r) function   
c  C(s)     the coefficient of the possible repulsive chemical bonding
c           set to 0 if is attractive chemical bonding
c  ifC(s)   if C(s) is free or fixed
c  Df       the form of damping function
c  ReBar    Approximate Re value used to calculate fd(R) function
c-----------------------------------------------------------------------
      READ(5,*) PSEL(ISTATE),PoEpZ(ISTATE),inAndOut(ISTATE)
      IF (PoEpZ(ISTATE) .GT. 0) THEN
        READ(5,*) nPB(ISTATE),nPu(ISTATE),nPq(ISTATE),nPL(ISTATE)
      END IF
      IF (ABS(PSEL(ISTATE)) .EQ. 10) THEN
        READ(5,*) C(ISTATE), cRp(ISTATE)
        READ(5,*) ifC(ISTATE)
      END IF
c-----------------------------------------------------------------------
c** For the GPEF, read in the parameters for w(R,a,b,k)
c
c  (A/B/K)GPEF(s) are parameters used to define the Generalized
c                 Potential Energy Function.
c
c                                 k    k      k      k
c                  w(R,k,a,b) = (R - RE )/(a*R + b*RE )
c-----------------------------------------------------------------------
      IF (ABS(PSEL(ISTATE)).EQ.5) THEN
        READ(5,*) AGPEF(ISTATE),BGPEF(ISTATE),KGPEF(ISTATE)
      ELSE IF (ABS(PSEL(ISTATE)).EQ.6 .OR. ABS(PSEL(ISTATE)).GE.9) THEN
        READ(5,*) numCn(ISTATE),dampBeta(ISTATE),Df
        IF (Df .EQ. 4) READ(5,*) ReBar
        READ(5,*) (nLR(ISTATE,i),Cn(ISTATE,i),i=1,numCn(ISTATE))
        IF (Df .EQ. 1) THEN
          WRITE(6,602) dampBeta(ISTATE)
        ELSE IF (Df .EQ. 2) THEN
          WRITE(6,606) dampBeta(ISTATE)
        ELSE IF (Df .EQ. 3 .OR. Df .EQ. 4) THEN
          WRITE(6,607) dampBeta(ISTATE)
        END IF
      ENDIF
  602 FORMAT(/' POTENTIAL ATTRACTION INCLUDES TANG-TOENNIES DAMPING FUNC
     1TION:'/1x,'(J.C.P.80,3726(1984))'/' 1 - EXP(-betaR)·SUM(betaR^K÷K!
     2)'/10X,'WHERE:   betaR = beta·R =',F11.8,'·R')
  606 FORMAT(/' POTENTIAL ATTRACTION INCLUDES DAMPED-DISPERSION DAMPING
     1FUNCTION:'/1x,'(M.P.,52,763(1984))'/' [1 - EXP(-2.1*beta*r/n-0.109
     2*(beta*r)²/n½)]^n'/10X,'WHERE:   beta = ',F11.8)
  607 FORMAT(/' POTENTIAL ATTRACTION INCLUDES DAMPED-DISPERSION DAMPING
     1FUNCTION:'/1x,'(M.P.,52,763(1984))'/' [1 - EXP(-3.97*beta*r/n-0.39
     2*(beta*r)²/n½)]^n'/10X,'WHERE:   beta = ',F11.8)
c-----------------------------------------------------------------------
c** Now to read in the trial dissociation energy and equilibrium
c   radial distance for the state.
c
c  DE(s)   is the dissociation energy for each state.
c
c  RE(s)   is the equilibrium radial distance for each state.
c          Note: If Re < 0 and PSEL < 0, then calculate Re from Be.
c
c  IFDE(s) indicates whether the dissociation energy will be:
c           = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c
c  IFRE(s) indicates whether the equilibrium radial distance will be:
c           = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c-----------------------------------------------------------------------
      READ(5,*) DE(ISTATE),RE(ISTATE)
      IF (ABS(PSEL(ISTATE)).EQ.6 .OR. ABS(PSEL(ISTATE)).GE.9) THEN
        Demo(ISTATE) = DE(ISTATE) 
      END IF
      IF (ABS(PSEL(ISTATE)).EQ.6 .OR. ABS(PSEL(ISTATE)).EQ.9) THEN
        READ(5,*) DeOffset,ReOffset
      END IF
      READ(5,*) IFDE(ISTATE),IFRE(ISTATE)
************************************************************************
* If PSEL = 7, uses the Modified Lennard-Jones-Gaussian Potential model,
* read in the maximum height Amp, the corresponding distance R0 and the
* FWHM of Gaussian function
* Amp(s)  is the amplitude of the Gaussian function
* R0(s)   is the distance where Gaussian function reaches its maximum
* FWHM(s) is full width half maximum
* ifAmp(s) indicates whether the amplitude fixed at read-in value
*          if ifAmp(s) = 1 or determined from fits
*          if ifAmp(s) <= 0
* ifR0(s) = 1  fixed R0(s) at read-in value
*         <= 0 determined from fits
*ifFWHM(s) = 1 fixed FWHM(s) at read-in value
*         <= 0 determined from fits   
************************************************************************
      IF (ABS(PSEL(ISTATE)).EQ.7) THEN
        READ(5,*) Amp(ISTATE),R0(ISTATE),FWHM(ISTATE)
        READ(5,*) ifAmp(ISTATE),ifR0(ISTATE),ifFWHM(ISTATE)
      END IF 
c-----------------------------------------------------------------------
c** Read in trial Born-Oppenheimer potential parameters, or for unknown
c   parameter values, read in specific molecular constants which will
c   be used to generate the first three trial parameters.
c
c  NBETA(s)    is the number of beta coefficients that will be used.
c
c  nBetaIn(s)  is the number of beta coefficients to represent the inner
c              wall of the potential curve.
c
c  BETA(b,s)   contains the parameters that define the Born-Oppenheimer
c              potential functions:
c            PSEL = 1-4, > 6: read in values are beta parameters for the
c                             family of morse potentials.
c            PSEL = 5: read in values are leading coefficients in the
c                      generalized potential energy function.
c            PSEL < 0: read in values are in fact molecular constants:
c                      We, WeXe, Be, ae.
c                      Note: If Be < 0 then calculate Be from Re.
c
c  IFBETA(b,s) indicates whether each BO potential parameter will be:
c               = 1: held fixed at read-in values.
c              <= 0: determined from fits.
c-----------------------------------------------------------------------
      READ(5,*) NBETA(ISTATE), nBetaIn(ISTATE) 
      IF (PSEL(ISTATE).LT.0) THEN
        READ(5,*) (BETA(I,ISTATE),I=0,3)
c
c** Calculate Be if only Re known and Re if only Be known.
c
        IF (BETA(2,ISTATE).LT.0.0d0) THEN
          BETA(2,ISTATE) = 16.85762908D0/ZMASS(3,1)/RE(ISTATE)**2
        ELSEIF (RE(ISTATE).LT.0.0d0) THEN
          RE(ISTATE) = DSQRT(16.85762908D0/ZMASS(3,1)/BETA(2,ISTATE))
        ENDIF
      ELSE
        READ(5,*) (BETA(I,ISTATE),I=0,NBETA(ISTATE))
      ENDIF
      READ(5,*) (IFBETA(I,ISTATE),I=0,NBETA(ISTATE))
c     IF (PoEpZ(ISTATE) .EQ. 3) THEN
c       BETA(0,ISTATE) = 0.0d0
c       IFBETA(0,ISTATE) = 1
c     END IF
c-----------------------------------------------------------------------
c** Read in constraints for the BETA function.
c
c  RSR(s)   determines if BETA(R) is fixed in the short range region.
c            <= 0.0 : Do not fix BETA(R)
c            >  0.0 : Fix BETA(R) at BETA(RSR) for R <= RSR
c
c  CNVAL(s) is the value of the first Cn constant for the MLJ potential.
c            <= 0.0 : Do not fix the Cn constant.
c            >  0.0 : Fix the Cn constant to be CNVAL using the
c                     switching function.
c
c                                            1
c                     fsw(R) =       ------------------
c                                 e^(ALPHA * (R - RX)) + 1
c
c  ALPHA(s) effectively defines the "strength" of switching function.
c
c  RX(s)    defines location at which fsw is one half: fsw(RX) = 0.5
c-----------------------------------------------------------------------
      READ(5,*) RSR, CNVAL(ISTATE), ALPHA(ISTATE), RX(ISTATE)
c-----------------------------------------------------------------------
      IF (RSR.GT.0.0d0) WRITE(6,640) RSR
c
c** Also determine radial position of short range constraint if chosen.
c
      IRSR(ISTATE) = (RSR - RMIN(ISTATE))/RH(ISTATE) + 1.001d0

      IF (ABS(PSEL(ISTATE)).EQ.4.AND.(CNVAL(ISTATE).GT.0.0d0)) THEN
        WRITE(6,645) NCN(ISTATE), CNVAL(ISTATE)
        WRITE(6,650) ALPHA(ISTATE), RX(ISTATE)
      END IF
*      IF ((PSEL(ISTATE).EQ.4).AND.(CNVAL(ISTATE).LT.0.0d0)) 
*       write(6,*) psel(istate), cnval(istate)
  640 FORMAT('Constrain BETA(R) at for values of R less than',F14.5)
  645 FORMAT('Constrain the C(',I2,') constant to be ',F14.5)
  650 FORMAT('Using switching function ',F14.5,2X,F14.5)
c-----------------------------------------------------------------------
c** Read in the Born-Oppenheimer breakdown corrrection terms.
c** First read in Adiabatic correction terms.
c
c  NABO(A/B)(s)    specifies the highest adiabatic correction term to
c                  be read in for atom A/B.
c                  >= 0: read in correction terms ABO(0 - NABO)
c                  <  0: do not read in any adiabatic correction terms.
c
c  ABO(A/B)(a,s)   are the adiabatic correction terms.
c
c  IFABO(A/B)(a,s) indicates whether each adiabatic correction term
c                  will be:
c                   = 1: held fixed at read-in values.
c                  <= 0: determined from fits.
c-----------------------------------------------------------------------
      READ(5,*) NUA(ISTATE)
      IF (NUA(ISTATE).GE.0) THEN
        READ(5,*) (UA(I,ISTATE),I=0,NUA(ISTATE))
        READ(5,*) (IFUA(I,ISTATE),I=0,NUA(ISTATE))
      ENDIF
      READ(5,*) NUB(ISTATE)
      IF (NUB(ISTATE).GE.0) THEN
        READ(5,*) (UB(I,ISTATE),I=0,NUB(ISTATE))
        READ(5,*) (IFUB(I,ISTATE),I=0,NUB(ISTATE))
      ENDIF
c-----------------------------------------------------------------------
c** Read in Rotational Non-Adiabatic correction terms.
c** For only one isotopomer, read only one set of Rotational
c   Non-Adiabatic correction terms.
c
c  NRNA(A/B)(s)    specifies the the highest rotational non-adiabatic
c                  correction term to be read for atom A/B.
c                  >= 0: read in correction terms RNA(0 - NRNA)
c                  <  0: do not read in any rotational correction terms.
c
c  RNA(A/B)(r,s)   are the rotational non-adiabatic correction terms.
c
c  IFRNA(A/B)(r,s) indicates whether each rotational non-adiabatic
c                  correction term will be:
c                   = 1: held fixed at read-in values.
c                  <= 0: determined from fits.
c-----------------------------------------------------------------------
      READ(5,*) NQA(ISTATE)
      IF (NQA(ISTATE).GE.0) THEN
        READ(5,*) (QA(I,ISTATE),I=0,NQA(ISTATE))
        READ(5,*) (IFQA(I,ISTATE),I=0,NQA(ISTATE))
      ENDIF
      READ(5,*) NQB(ISTATE)
      IF (NQB(ISTATE).GE.0) THEN
        READ(5,*) (QB(I,ISTATE),I=0,NQB(ISTATE))
        READ(5,*) (IFQB(I,ISTATE),I=0,NQB(ISTATE))
      ENDIF
c-----------------------------------------------------------------------
c** For fits to a single isotopomer, we know that there are no adiabatic
c   correction terms and only one set of rotational correction terms can
c   determined.
c
      IF (NISTP.EQ.1) THEN
        NUA(ISTATE) = -1
        NUB(ISTATE) = -1
        NQB(ISTATE) = -1
      ENDIF
c
c** For fits to homonuclear diatomics, there is only one adiabatic
c   correction function and one rotational non-adiabatic correction
c   function, so ignore the second set of correction parameters.
c
      IF (AN(1).EQ.AN(2)) THEN
        NUB(ISTATE) = -1
        NQB(ISTATE) = -1
      ENDIF
c-----------------------------------------------------------------------
c** Read in lambda doubling correction terms
      READ(5,*) nLamdCoef(ISTATE)
      IF (nLamdCoef(ISTATE).GE.0) THEN
        READ(5,*) (lamdCoef(I,ISTATE),I=0,nLamdCoef(ISTATE))
        READ(5,*) (ifLamdCoef(I,ISTATE),I=0,nLamdCoef(ISTATE))
      ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
