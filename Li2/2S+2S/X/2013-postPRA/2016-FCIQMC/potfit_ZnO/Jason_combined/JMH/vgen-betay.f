c**********************************************************************
      SUBROUTINE VGEN(ISTATE,Rd,Vd,DATN)
c**********************************************************************
c** This subroutine will generate one of five possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program. These four families are:
c     - The Generalized Morse Oscillator
c     - The Expanded Morse Oscillator
c     - The Modified Morse Oscillator
c     - The Modified Lennard-Jones Oscillator
c     - The Generalized Potential Energy Function
c** This subroutine can also calculate the trial parameters for each of
c   the four families of analytical potential functions from the low
c   order molecular constants: Be, We, WeXe, ae.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 1997  by  Jenning Seto       (ver of 12/10/1997) +++        
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
c    Rd      is the distance, if Rd > 0, calculate that point potential
c            value Vd
c** On entry via common blocks:
c    PSEL   is used to choose the type of potential used.
c      |PSEL| = 0 : Do forward calculation.
c      |PSEL| = 1 : Use the Generalized Morse Potential.
c      |PSEL| = 2 : Use the Expanded Morse Potential.
c      |PSEL| = 3 : Use the Modified Morse Potential.
c      |PSEL| = 4 : Use the Modified Lennard-Jones Potential.
c      |PSEL| = 5 : Use the Generalized Potential Energy Function.
c      |PSEL| = 6 : Use the Expanded Morse Long-Range Potential. 
c      |PSEL| = 7 : Use the Modified Lennard-Jones-Gaussian Potential. 
c      |PSEL| = 8 : Use the Total Expanded Morse Potential.
c      |PSEL| = 9 : Total Expanded Morse Long-Range Potential
c                   (trial version of DELR).
c      |PSEL| = 10: Use the Double-Exponential Long-Range Potential.
c           < 0 : Read in parmaters are in fact the molecular
c                 constants for that molecule and the trial exponential
c                 paramters must be calculated from these values.
c    NBETA  is the number of beta parameters that will be used.
c    NCN    is the power of the dampening term.
c    nPB(s)   the value of power m for beta function
c    nPq(s)   the value of power m for q(r) function
c    nPu(s)   the value of power m for u(r) function
c    nPL(s)   the value of power m for Lambda(r) function
c    YPARM  is the array of parameters used for the potential.
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    BETA   is the array of potential parameters for each state
c           or (if PSEL < 0) selected molecular constants.
c    NDATPT is the number of meshpoints used for the array.
c
c  PoEpZ(s) specifies type of expansion used for Born-Oppenheimer
c           correction, lambda doubling and beta functions.
c            = 0 : use half of OT expansion (R-RE)/(R+RE)
c            =-1 : use Coxon expansion (R-RE)
c            > 0 : use new variable (R^m-Re^m)/(R^m+Re^m) for BOB
c            = 2 : use new variable (R^m-Re^m)/(R^m+Re^m) for beta
c            = 3 : use exponent beta*y and Cn constraint as
c                  beta=beta(sum)*(1-y) + beta(inf)*y
c  inAndOut(s) indicates the terms of the beta function used in inner wall
c            = 1: beta function uses N_S terms at R < Re 
c                 beta function uses N_L terms at R >= Re  (N_S < N_L)
c            = 0: beta function of repulsive exponential uses N_S terms
c                 beta function of attractive exponential uses N_L terms
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    R      is the distance array
c    VPOT   is the potential that is generated.
c    BMORSE is used to contain the beta function.
c    DVDD   is the array of the partial derivatives of V with DE.
c    DVDR   is the array of the partial derivatives of V with RE.
c    DVDB   is the array of the partial derivatives of V with the
c           various beta parameters.
c
c** Temp:
c    BTEMP     is used to represent the sum used for dV/dRe.
c              is used in GPEF for De calculations.
c    BINF      is used to represent the beta infinity value.
c    RTEMP     is used to represent (R-Re)/(R+Re) or R-Re.
c                                         n    -BMORSE*RTEMP
c    XTEMP     is used to represent (Re/R)  * e 
c                                         n    -BINF*RTEMP
c    XITEMP    is used to represent (Re/R)  * e
c    PBTEMP    is used to calculate dV/dBi.
c    PETEMP    is used to calculate dV/dBi.
c    AZERO     is used for the trial exponential calculations.
c    AONE      is used for the trial exponential calculations.
c    ATWO      is used for the trial exponential calculations.
c    DBETA     is used in the MMO trial exponential calculations.
c    AZTEMP    is used in the MMO trial exponential calculations.
c              is used in the GPEF = (a+b)/k
c    AOTEMP    is used in the GPEF = [a(k+1)-b(k-1)]/k
c    ATTEMP    is used in the GPEF = [a^2(k+1)-b^2(k-1)]/k
c    ARTEMP    is used in the GPEF = [a^3(k+1)-b^3(k-1)]/k
c    FSW       is used to represent the MLJ switching function.
c    RX(ISTATE) is the Rx value used in the MLJ switching function.
c
      INCLUDE 'arrsizes.h'

      INTEGER ISTATE, IPV, I, J, NPPS(NSTATEMX),l,step
      COMMON /NPPSBLK/NPPS 

      REAL*8 BTEMP,BINF,RDVAL,RTEMP,RM2,XTEMP,XITEMP,PBTEMP,PETEMP,
     1  AZERO,AONE,ATWO,DBETA(0:2),AZTEMP,AOTEMP,ATTEMP,ARTEMP,FSW,
     2  Xtemp2,Btemp2,betaCoef,BetaInfinity,BMtemp,BMtemp2,Rminimum,RM4
     3  ,PBtemp2,ZTEMP
c
c** Born-Oppenheimer Potential Variables.
c
      LOGICAL emoPrinted,outputPrint

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
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,VP,VPP,EFP,EFPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL      
      
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
     3  ,nBetaIn(NSTATEMX),nPB(NSTATEMX),inAndOut(NSTATEMX)
     4  ,nPq(NSTATEMX),nPu(NSTATEMX),nPL(NSTATEMX)

      COMMON /nPwBLK/nPB,nPq,nPu,nPL
      COMMON /nBetaInBLK/nBetaIn
      COMMON /inAndOutBLK/inAndOut

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),bMorsIn(NPNTMX,NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      REAL*8 R(NPNTMX,NSTATEMX),Uad,Ubd,Qad,Qbd,Ld,Rd,Vd,dVdPk(NPOTMX),
     &       JFCT,JfL
      INTEGER DATN, NBAND, IISTP, efParity, NDatPtS
      COMMON /RBLK/R
      COMMON /dVdPkBLK/dVdPk
      COMMON /OUTPUT/outputPrint

c***********************************************************************
c** For Expanded Morse Long-Range Potential
      INTEGER numCn(NSTATEMX),nLR(NSTATEMX,numCnMax),tempN,Df
      REAL*8 dampBeta(NSTATEMX),Cn(NSTATEMX,numCnMax),Demo(NSTATEMX),
     &  DeOffset,ReOffset,Dn(NSTATEMX,numCnMax,NPNTMX),
     &  Vlr(NPNTMX,NSTATEMX),Dm,Dmp,Dmpp,A(NSTATEMX),B(NSTATEMX),
     &  dAdRe,dBdRe,dVdBtemp,cTemp,fd(NPNTMX,NSTATEMX),VlrB,ReBar

      COMMON /DELRBLK/numCn,nLR,dampBeta,Cn,ReBar,Df
      COMMON /OffSetBLK/Demo,DeOffset,ReOffset

c** The possible repulsive chemical bonding
      REAL*8 C(NSTATEMX),DVDC(NPNTMX,NSTATEMX)
      INTEGER ifC(NSTATEMX)
      COMMON /repulsiveBLK/C,ifC,DVDC

      REAL*8 cRp(NSTATEMX)
      COMMON /cRpBLK/cRp
c***********************************************************************

c** For Modified Lennard-Jones-Gaussian Potential
      INTEGER ifAmp(NSTATEMX),ifR0(NSTATEMX),ifFWHM(NSTATEMX)
      REAL*8 Amp(NSTATEMX),R0(NSTATEMX),FWHM(NSTATEMX)
      COMMON /MLJG/ifAmp,ifR0,ifFWHM,Amp,R0,FWHM

      REAL*8 dVdA(NPNTMX,NSTATEMX),dVdR0(NPNTMX,NSTATEMX)
     &  ,dVdW(NPNTMX,NSTATEMX),Gauss(NPNTMX,NSTATEMX)
c***********************************************************************
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
     3  IFQA(0:NRNAMX,NSTATEMX),IFQB(0:NRNAMX,NSTATEMX)

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

c
c** Partial derivative arrays.
c
      REAL*8 DVDD(NPNTMX,NSTATEMX),DVDR(NPNTMX,NSTATEMX),
     1  DVDB(0:NBETAMX,NPNTMX,NSTATEMX),
     2  DVDUA(0:NABOMX,NPNTMX,NSTATEMX),DVDUB(0:NABOMX,NPNTMX,NSTATEMX),
     3  DUADR(NPNTMX,NSTATEMX),DUBDR(NPNTMX,NSTATEMX),
     4  DVDQA(0:NRNAMX,NPNTMX,NSTATEMX),DVDQB(0:NRNAMX,NPNTMX,NSTATEMX),
     5  DQADR(NPNTMX,NSTATEMX),DQBDR(NPNTMX,NSTATEMX)
     6  ,dVdL(0:NRNAMX,NPNTMX,NSTATEMX),dLdRe(NPNTMX,NSTATEMX)
     7  ,DBDB(0:NBETAMX,NPNTMX,NSTATEMX),DBDR(NPNTMX,NSTATEMX)

      REAL*8 DVBODP(NPNTMX,NPOTMX),DUADP(NPNTMX,NPOTMX),
     1  DUBDP(NPNTMX,NPOTMX),DQADP(NPNTMX,NPOTMX),DQBDP(NPNTMX,NPOTMX)
     2  ,dLdP(NPNTMX,NPOTMX),DBDP(NPNTMX,NPOTMX)

      COMMON /PDVBLK/DVBODP,DUADP,DUBDP,DQADP,DQBDP,DBDP
      COMMON /PdVLBLK/dLdP 

c      REAL*8 VdT
c      LOGICAL out
c      COMMON /VdTBLK/VdT,out
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Initializing variables.
c
c** Normally data point starts from 1
      NDatPtS = 1
c** When calculating only one potential point
      IF (Rd .GT. 0) THEN
        NDATPT(ISTATE) = NDATPT(ISTATE) + 1
        NDatPtS = NDATPT(ISTATE)
      END IF
      Vd = 0.0d0
      Uad = 0.0d0
      Ubd = 0.0d0
      Qad = 0.0d0
      Qbd = 0.0d0
      Ld  = 0.0d0
      PBTEMP = 0.0d0
      PETEMP = 0.0d0
      DO I = 1,NPNTMX
        BMORSE(I,ISTATE) = 0.0d0
        bMorsIn(I,ISTATE) = 0.0d0 
      END DO
c
c** If PSEL < 0, then generate trial parameters from the given molecular
c   constatnts.
c
c   These calculated values for a(0,1,2) are not the same as those in
c   the manual. These values are used to simplify calculations.
c
c                        a0 = (We^2/4*DE*Be)^1/2
c                        a1 = -(We*ae/6*Be^2)-1
c                        a2 = WeXe/Be
c
      IF (PSEL(ISTATE).LT.0) THEN
        AZERO = DSQRT((BETA(0,ISTATE)**2)
     1          /(4.0d0*DE(ISTATE)*BETA(2,ISTATE)))
        AONE = -BETA(0,ISTATE)*BETA(3,ISTATE)
     1         /(6.0d0*BETA(2,ISTATE)**2)-1.0d0
        ATWO = BETA(1,ISTATE)/BETA(2,ISTATE)
      ENDIF

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For Forward Calculation.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (PSEL(ISTATE).EQ.0) THEN
        DO 5 I = 1,NDATPT(ISTATE) 
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          DO 7 J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
    7     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            BMORSE(I,ISTATE) = BMORSE(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
          BTEMP = 0.0d0 
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            DO 9 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
   9        CONTINUE
          ELSE
            DO 11 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  11        CONTINUE
            BTEMP = BTEMP*(RDVAL-RE(ISTATE)) 
          END IF
c
c** Now to calculate the remaining partial derivatives.
c
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) /
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (1.0d0-XTEMP)
          DVDR(I,ISTATE) = -PBTEMP * (BMORSE(I,ISTATE)+PETEMP*BTEMP)
          DO 13 J = 0,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * (RDVAL-RE(ISTATE)) * RTEMP**J
   13    CONTINUE
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
   5   CONTINUE
      ENDIF
******7***** End Forward Calculation ***********************************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Generalized Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-1) THEN
        PSEL(ISTATE) = 1
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/(2.0d0*RE(ISTATE)**2)
        BETA(2,ISTATE) = AZERO*(3.0d0*AONE*(AZERO+AONE)+2.0d0*AZERO**2
     1                    -2.0d0*ATWO)/(6.0d0*RE(ISTATE)**3)
        DO 10 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
   10   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.1) THEN
        DO 20 I = 1,NDATPT(ISTATE) 
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          RTEMP = (RDVAL-RE(ISTATE))
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          DO 30 J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
   30     CONTINUE
          XTEMP = DEXP(-BMORSE(I,ISTATE)*RTEMP)
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = (1d0-XTEMP)**2 - 1.0d0
          BTEMP = BETA(0,ISTATE)
          DO 40 J = 1,NBETA(ISTATE)
            BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J+1) * RTEMP**J 
   40     CONTINUE
          DVDR(I,ISTATE) = -2.0d0*DE(ISTATE)*XTEMP*(1-XTEMP)*BTEMP
          DVDB(0,I,ISTATE) = 2.0d0*DE(ISTATE)*XTEMP*(1-XTEMP)*RTEMP
          DO 50 J = 1,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = DVDB(J-1,I,ISTATE) * RTEMP
   50     CONTINUE
*          DBDR(I,ISTATE) = 0.0d0 
*          DO J = 1,NBETA(ISTATE)
*            DBDR(I,ISTATE) = DBDR(I,ISTATE)
*     1                     - DBLE(J)*BETA(J,ISTATE)*RTEMP**(J-1)
*          END DO
*          DBDB(0,I,ISTATE) = 1.0d0 
*          DO J = 1,NBETA(ISTATE)
*            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
*          END DO
   20   CONTINUE
      ENDIF
******7***** End Generalized Morse Oscillator Potential Function *******

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-2) THEN
        PSEL(ISTATE) = 2
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/RE(ISTATE)
        BETA(2,ISTATE) = AZERO*((2.0d0*AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2-ATWO)*4.0d0/3.0d0)/RE(ISTATE)
        DO 60 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
   60   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.2) THEN
        DO 70 I = 1,NDATPT(ISTATE) 
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          DO 80 J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
   80     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            BMORSE(I,ISTATE) = BMORSE(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
          BTEMP = 0.0d0 
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            DO 90 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  90        CONTINUE
          ELSE
            DO 95 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  95        CONTINUE
            BTEMP = BTEMP*(RDVAL-RE(ISTATE)) 
          END IF
c
c** Now to calculate the remaining partial derivatives.
c
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) /
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (1.0d0-XTEMP)
          DVDR(I,ISTATE) = -PBTEMP * (BMORSE(I,ISTATE)+PETEMP*BTEMP)
          DO 100 J = 0,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * (RDVAL-RE(ISTATE)) * RTEMP**J
  100     CONTINUE
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
*          DBDR(I,ISTATE) = -2.0d0*BETA(1,ISTATE)*RDVAL/
*     1                      (RDVAL+RE(ISTATE))**2
*          DO J = 2,NBETA(ISTATE)
*            DBDR(I,ISTATE) = DBDR(I,ISTATE)
*     1                     - DBLE(J)*BETA(J,ISTATE)*RTEMP**(J-1)
*          END DO
*          DBDB(0,I,ISTATE) = 1.0d0 
*          DO J = 1,NBETA(ISTATE)
*            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
*          END DO
   70   CONTINUE
      ENDIF
******7***** End Expanded Morse Oscillator Potential Function **********

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For Dulick's Modified Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
c** Need recursive algorithm to calculate the exponential parameters.
c
      IF (PSEL(ISTATE).EQ.-3) THEN
        PSEL(ISTATE) = 3
        AZTEMP = AZERO
c
c** Finding corrections to the trial exponential parameters.
c
  102   DBETA(0) = 2.0d0*AZTEMP
        DBETA(1) = 2.0d0*AZTEMP*(AZTEMP+AONE+1.0d0)
        DBETA(2) = AZTEMP*(4.0d0*(AONE+1.0d0)*(AZTEMP+AONE)
     1             +(AZTEMP**2-ATWO)*8.0d0/3.0d0+2.0d0)
c
c** Making corrections to the trial exponential parameters.
c
        DO 104 I=0,2
          BETA(I,ISTATE) = DBETA(I)
  104   CONTINUE
        BINF = BETA(0,ISTATE)+BETA(1,ISTATE)+BETA(2,ISTATE)
        AZTEMP = AZERO * (1.0d0 - DEXP(-BINF))
c
c** Check to see if the relative change is less than 0.1% of the 
c   value of the parmaeter itself.
c
        DO 106 I=0,2
          IF (DABS(DBETA(I)/BETA(I,ISTATE)-1.0d0).GT.1.0d-3) GOTO 102
  106   CONTINUE
        DO 110 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  110   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.3) THEN
        DO 120 I = 1,NDATPT(ISTATE) 
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          BINF = BETA(0,ISTATE)
          DO 130 J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
            BINF = BINF + BETA(J,ISTATE)
  130     CONTINUE
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*RTEMP)
          XITEMP = DEXP(-BINF)
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*((1.0d0-XTEMP)/(1.0d0-XITEMP))**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = ((1.0d0-XTEMP)/(1.0d0-XITEMP))**2 - 1.0d0
          BTEMP = BETA(0,ISTATE)
          DO 140 J = 1,NBETA(ISTATE)
            BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J+1) * RTEMP**J 
  140     CONTINUE
c
c** Now to calculate the remaining partial derivatives.
c
          DVDR(I,ISTATE) = -2.0d0 * DE(ISTATE) * XTEMP * 
     1      (1.0d0-XTEMP) * 2.0d0*BTEMP*RDVAL/(RDVAL+RE(ISTATE)**2)
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP *
     1       (1.0d0-XTEMP) / (1.0d0-XITEMP)**2
          PETEMP = (XITEMP*(1.0d0-XTEMP))/(XTEMP*(1.0d0-XITEMP))
          DO 150 J = 0,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * (RTEMP**(J+1) - PETEMP)
  150     CONTINUE
  120   CONTINUE
      ENDIF
******7***** End Dulick's Modified Morse Potential Function ************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Modified Lennard-Jones Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-4) THEN
        PSEL(ISTATE) = 4
        BETA(0,ISTATE) = 2.0d0*AZERO-2.0d0*DBLE(NCN(ISTATE))
        BETA(1,ISTATE) = 2.0d0*AZERO*(AZERO+AONE+1.0d0)
        BETA(2,ISTATE) = AZERO*(4.0d0*(AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2.0d0-ATWO)*8.0d0/3.0d0+2.0d0)
     2                    -(2.0d0/3.0d0)*DBLE(NCN(ISTATE))
        DO 160 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  160   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.4) THEN
        DO 170 I = 1,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          bMorsIn(I,ISTATE) = BETA(0,ISTATE)
          DO J = 1,nBetaIn(ISTATE)
            bMorsIn(I,ISTATE)=bMorsIn(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
          END DO
          BMORSE(I,ISTATE) = bMorsIn(I,ISTATE)
          DO J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
          END DO
          IF (PoEpZ(ISTATE) .EQ. 3) THEN
            betaCoef = 2.0d0 * RE(ISTATE)**nPB(ISTATE) /
     &                 (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))  
            BETA(NBETA(ISTATE)+1,ISTATE) = DLOG(2.0d0*DE(ISTATE)*
     &                     RE(ISTATE)**NCN(ISTATE)/CNVAL(ISTATE)) 
            BetaInfinity = BETA(NBETA(ISTATE)+1,ISTATE)
            bMorsIn(I,ISTATE) = bMorsIn(I,ISTATE) * betaCoef
            BMORSE(I,ISTATE) = BETA(NBETA(ISTATE)+1,ISTATE) * RTEMP + 
     &        BMORSE(I,ISTATE) * betaCoef
          END IF 
c-----------------------------------------------------------------------
c** If chosen, constrain the Cn value to CNVAL by using a switching
c   function.
c
c             BETAsw(R) = fsw(R) [ BETA(R) - BINF ] + BINF
c
c                fsw(R) = 1 / [e^(ALPHA'*(R - Rx)) + 1]
c
          IF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).GT.0.0d0) THEN
            BINF = DLOG(2.0d0*DE(ISTATE)*RE(ISTATE)**NCN(ISTATE)
     1             /CNVAL(ISTATE))
            BetaInfinity = BINF
            FSW = 1.0d0 / (1.0d0+DEXP(ALPHA(ISTATE)*(RDVAL-RX(ISTATE))))
c
c** Update the generalized morse exponent.
c
            BMORSE(I,ISTATE) = FSW*(BMORSE(I,ISTATE) - BINF) + BINF
*          ENDIF
c**********************************
c
c** If chosen, constrain the Cn value to BVC by fixing the NBVCth BETA
c   parameter BETA(NBVC).
c
          ELSEIF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).LT.0.0d0) THEN
            BINF = BETA(0,ISTATE)
            DO J = 1,NBETA(ISTATE)
              BINF = BINF + BETA(J,ISTATE)
            END DO
            BETA(NBETA(ISTATE)+1,ISTATE) = - BINF +
     1         DLOG(2.0d0*DE(ISTATE)
     2        *RE(ISTATE)**NCN(ISTATE)/DABS(CNVAL(ISTATE)))
            BetaInfinity = BINF + BETA(NBETA(ISTATE)+1,ISTATE)
c
c** Update the generalized morse exponent.
c
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)
     1        +BETA(NBETA(ISTATE)+1,ISTATE)*RTEMP**(NBETA(ISTATE)+1)
          ENDIF
c**********************************

c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            BMORSE(I,ISTATE) = BMORSE(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c** Calculate some temporary variables.
c
          IF (PoEpZ(ISTATE).NE.3) THEN
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     1          DEXP(-BMORSE(I,ISTATE)*RTEMP)
              BMtemp = BMORSE(I,ISTATE)
            ELSE
              IF (RDVAL .LT. RE(ISTATE)) THEN
                XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     &            DEXP(-bMorsIn(I,ISTATE)*RTEMP)  
                BMtemp = bMorsIn(I,ISTATE)
              ELSE
                XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     1            DEXP(-BMORSE(I,ISTATE)*RTEMP) 
                BMtemp = BMORSE(I,ISTATE)
              END IF
            END IF
          ELSE
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) * 
c    &          DEXP(-BMORSE(I,ISTATE))
     &          DEXP(-BMORSE(I,ISTATE)*RTEMP)
            ELSE
              IF (RDVAL .LT. RE(ISTATE)) THEN
                XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
c    &            DEXP(-bMorsIn(I,ISTATE))
     &            DEXP(-bMorsIn(I,ISTATE)*RTEMP)
              ELSE
                XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
c    &            DEXP(-BMORSE(I,ISTATE))
     &            DEXP(-BMORSE(I,ISTATE)*RTEMP)
              END IF
            END IF 
          END IF
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate temporary variables for the partial derivatives.
c
          BTEMP = 0.0d0
          IF (PoEpZ(ISTATE).NE.3) THEN 
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              DO 190 J = 1,NBETA(ISTATE)
                BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J 
  190         CONTINUE
            ELSE
              IF (RDVAL .LT. RE(ISTATE)) THEN
                DO J = 1,nBetaIn(ISTATE)
                  BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
                END DO
              ELSE
                DO J = 1,NBETA(ISTATE)
                  BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
                END DO
              END IF
            END IF  
          ELSE
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              DO J = 1,NBETA(ISTATE)
                BTEMP = BTEMP + BETA(J,ISTATE) * (2.0d0*DBLE(J) *
     &            RE(ISTATE)**nPB(ISTATE) - (RDVAL**nPB(ISTATE) - 
     &            RE(ISTATE)**nPB(ISTATE))) / (RDVAL**nPB(ISTATE) +
     &            RE(ISTATE)**nPB(ISTATE)) * RTEMP**(J-1)
              END DO
            ELSE
              IF (RDVAL .LT. RE(ISTATE)) THEN
                DO J = 1,nBetaIn(ISTATE)
                  BTEMP = BTEMP + BETA(J,ISTATE) * (2.0d0*DBLE(J) *
     &              RE(ISTATE)**nPB(ISTATE) - (RDVAL**nPB(ISTATE) -
     &              RE(ISTATE)**nPB(ISTATE))) / (RDVAL**nPB(ISTATE) +
     &              RE(ISTATE)**nPB(ISTATE)) * RTEMP**(J-1)
                END DO
              ELSE
                DO J = 1,NBETA(ISTATE)
                  BTEMP = BTEMP + BETA(J,ISTATE) * (2.0d0*DBLE(J) *
     &              RE(ISTATE)**nPB(ISTATE) - (RDVAL**nPB(ISTATE) -
     &              RE(ISTATE)**nPB(ISTATE))) / (RDVAL**nPB(ISTATE) +
     &              RE(ISTATE)**nPB(ISTATE)) * RTEMP**(J-1)
                END DO
              END IF
            END IF
          END IF 
c**********************************
          IF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).LT.0.0d0) THEN
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              BTEMP = BTEMP + BETA(NBETA(ISTATE)+1,ISTATE) * 
     &          DBLE(NBETA(ISTATE)+1)*RTEMP**(NBETA(ISTATE)+1)
            ELSE
              IF (RDVAL .GT. RE(ISTATE)) THEN
                BTEMP = BTEMP + BETA(NBETA(ISTATE)+1,ISTATE) *
     &            DBLE(NBETA(ISTATE)+1)*RTEMP**(NBETA(ISTATE)+1)
              END IF
            END IF
          END IF 
c**********************************
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) / 
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (1.0d0-XTEMP)
c-----------------------------------------------------------------------
c** Now to calculate the partial derivatives.
c** Case when switching function constraint is used.
c
          IF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).GT.0.0d0) THEN
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
     1                       + PBTEMP/DE(ISTATE)*(1-FSW)
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1        *(1-RTEMP+RTEMP*FSW) + PETEMP*(BMtemp+FSW*BTEMP))
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * RTEMP**(J+1) * FSW
            END DO
c**********************************
c** Case when polynomial constraint is used.
c
          ELSEIF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).LT.0.0d0) THEN
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
     1                     + PBTEMP/DE(ISTATE)*RTEMP**(NBETA(ISTATE)+2)
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1        *(1-RTEMP**(NBETA(ISTATE)+2))+PETEMP*(BMtemp+BTEMP))
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * (RTEMP**(J+1) - 
     &          RTEMP**(NBETA(ISTATE)+2))
            END DO
c**********************************
c** Case when new exponent function with Binf constraint is used
          ELSEIF (PoEpZ(ISTATE). EQ. 3) THEN
c           DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
     1                      + 2.d0*RTEMP**2 * XTEMP*(1.D0-XTEMP)
c           DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
c    1            + PETEMP * (BETA(NBETA(ISTATE)+1,ISTATE) + BTEMP))
            DBDY=BetaInfinity
            DO J = 0,NBETA(ISTATE)
              DBDY=DBDY+J*BETA(J,ISTATE)*RTEMP**(J-1)*(1.D0-RTEMP)
     1          - BETA(J,ISTATE)*RTEMP**J
            END DO
            DBDY = BMORSE(I,ISTATE) + RTEMP*DBDY
            DVDR(I,ISTATE) = - PBTEMP / RE(ISTATE) 
     1          * (1.D0-RTEMP) * (1.D0+RTEMP) 
     2          * (DBLE(NCN(ISTATE)) + 0.5D0*nPB(ISTATE)*DBDY)
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * betaCoef * RTEMP**(J+1)
            END DO
            DBDR(I,ISTATE)=-PETEMP*(BETA(NBETA(ISTATE)+1,ISTATE)+BTEMP)
            DO J = 0,NBETA(ISTATE)
              DBDB(J,I,ISTATE) = betaCoef * RTEMP**(J+1)
            END DO
c**********************************
c** Case when no constraints are used.
c
          ELSE
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1                       + PETEMP * (BMtemp + BTEMP))
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * RTEMP**(J+1)
            END DO
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
            DBDR(I,ISTATE) = BETA(1,ISTATE)
            DO J = 2,NBETA(ISTATE)
              DBDR(I,ISTATE) = DBDR(I,ISTATE)
     1                       + DBLE(J)*BETA(J,ISTATE)*RTEMP**(J-1)
            END DO
            DBDR(I,ISTATE) = DBDR(I,ISTATE) * (-PETEMP)
            DBDB(0,I,ISTATE) = 1.0d0
            DO J = 1,NBETA(ISTATE)
              DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
            END DO
          END IF
  170   CONTINUE
      END IF
******7***** End Modified Lennard-Jones Potential Function *************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Generalized Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If PSEL = -5, then calculate GPEF parameters from the given
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-5) THEN
        PSEL(ISTATE) = 5
        AZERO = (BETA(1,ISTATE)**2)/(4.0d0*BETA(0,ISTATE))
        AONE = -BETA(3,ISTATE)*BETA(1,ISTATE)
     1         /(6.0d0*BETA(0,ISTATE)**2)-1.0d0
        ATWO = (5.0d0*AONE**2)/4.0d0 
     1       - (2.0d0*BETA(2,ISTATE))/(3.0d0*BETA(0,ISTATE))
        AZTEMP = (AGPEF(ISTATE)+BGPEF(ISTATE))/DBLE(KGPEF(ISTATE))
        AOTEMP = (AGPEF(ISTATE)*(DBLE(KGPEF(ISTATE)+1))
     1           - BGPEF(ISTATE)*(DBLE(KGPEF(ISTATE)-1)))
     2           /DBLE(KGPEF(ISTATE))
        ATTEMP = (AGPEF(ISTATE)**2*(DBLE(KGPEF(ISTATE)+1))
     1           + BGPEF(ISTATE)**2*(DBLE(KGPEF(ISTATE)-1)))
     2           /DBLE(KGPEF(ISTATE))
*        ARTEMP = (AGPEF(ISTATE)**3*(DBLE(KGPEF(ISTATE)+1))
*     1           + BGPEF(ISTATE)**3*(DBLE(KGPEF(ISTATE)-1)))
*     2           /DBLE(KGPEF(ISTATE))
        BETA(0,ISTATE) = AZERO*AZTEMP**2
        BETA(1,ISTATE) = AONE*AZTEMP + AOTEMP
        BETA(2,ISTATE) = ATWO*AZTEMP**2 + 3.0d0*AONE*AZTEMP*AOTEMP/2.0d0
     1     + 7.0d0*AOTEMP**2/12.0d0 + ATTEMP/3.0d0
*        BETA(3,ISTATE) = ATHREE*AZTEMP**3 + 2.0d0*ATWO*AZTEMP**2*AOTEMP
*          + AONE*AZTEMP*(5.0d0*AOTEMP**2/4.0d0 + ATTEMP/2.0d0)
*          + AOTEMP**3/4.0d0 + 5.0d0*AOTEMP*ATTEMP/12.0d0 + ARTEMP/6.0d0
        DO 210 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  210   CONTINUE
      ENDIF
      
      IF (PSEL(ISTATE).EQ.5) THEN
c
c** First, we calculate the Dissociation Energy
c
        IF (AGPEF(ISTATE).NE.0.0d0) THEN
          BTEMP = 1.0d0
          DO 212 I = 1,NBETA(ISTATE)
            BTEMP = BTEMP + BETA(I,ISTATE)/AGPEF(ISTATE)**I
 212      CONTINUE
          BTEMP = BTEMP*BETA(0,ISTATE)/AGPEF(ISTATE)**2
        ENDIF
        IF (IFDE(ISTATE).LE.0) DE(ISTATE) = BTEMP
c
c** If considering a constrained De, then calculate additional
c   constrained parameter g(n+1)
c
        IF (IFDE(ISTATE).EQ.1) THEN
          BETA(NBETA(ISTATE)+1,ISTATE) = DE(ISTATE) / BTEMP
        ENDIF

        DO 215 I = 1,NDATPT(ISTATE) 
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          RTEMP = (RDVAL**KGPEF(ISTATE)-RE(ISTATE)**KGPEF(ISTATE))
     1           / (AGPEF(ISTATE)*RDVAL**KGPEF(ISTATE)
     2             + BGPEF(ISTATE)*RE(ISTATE)**KGPEF(ISTATE))
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = 1.0d0 
          DO 220 J = 1,NBETA(ISTATE)
            VPOT(I,ISTATE) = VPOT(I,ISTATE) + BETA(J,ISTATE)*RTEMP**J 
  220     CONTINUE
          VPOT(I,ISTATE) = VPOT(I,ISTATE)*BETA(0,ISTATE)*RTEMP**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDR(I,ISTATE) = 2.0d0
          DO 225 J = 1,NBETA(ISTATE)
            DVDR(I,ISTATE) = DVDR(I,ISTATE) 
     1                     + DBLE(J+2)*BETA(J,ISTATE)*RTEMP**J
  225     CONTINUE
          DVDR(I,ISTATE) = -DVDR(I,ISTATE) * BETA(0,ISTATE) * RTEMP
     1       * (AGPEF(ISTATE)+BGPEF(ISTATE))*DBLE(KGPEF(ISTATE))
     2       * RDVAL**KGPEF(ISTATE) * RE(ISTATE)**(KGPEF(ISTATE)-1)
     3       / (AGPEF(ISTATE)*RDVAL**KGPEF(ISTATE)
     4       + BGPEF(ISTATE)*RE(ISTATE)**KGPEF(ISTATE))
          DVDB(0,I,ISTATE) = 1.0d0 
          DO 230 J = 1,NBETA(ISTATE)
            DVDB(0,I,ISTATE) = DVDB(0,I,ISTATE)+BETA(J,ISTATE)*RTEMP**J 
  230     CONTINUE
          DVDB(0,I,ISTATE) = DVDB(0,I,ISTATE)*RTEMP**2
          DVDB(1,I,ISTATE) = BETA(0,ISTATE)*RTEMP**3
          DO 235 J = 2,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = DVDB(J-1,I,ISTATE)*RTEMP
  235     CONTINUE
  215   CONTINUE
      ENDIF
******7***** End Generalized Potential Energy Function *****************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Longe-Range Potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-6) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        PSEL(ISTATE) = 6
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/RE(ISTATE)
        BETA(2,ISTATE) = AZERO*((2.0d0*AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2-ATWO)*4.0d0/3.0d0)/RE(ISTATE)
        DO 360 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  360   CONTINUE
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.6) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        DO 370 I = NDatPtS,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL 
          IF (Rd .GT. 0) RDVAL = Rd
          IF (PoEpZ(ISTATE) .LT. 2) THEN 
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE 
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          DO 380 J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  380     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            BMORSE(I,ISTATE) = BMORSE(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
          Vlr(I,ISTATE) = 0.0d0
          DO J = 1,numCn(ISTATE)
            tempN = nLR(ISTATE,J)
            CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
            Dn(ISTATE,J,I) = Dm
*            WRITE(7,800) RDVAL,Dn(ISTATE,J,I),nLR(ISTATE,J)
*  800       FORMAT(F10.4,F16.10,I6)
            Vlr(I,ISTATE) = Vlr(I,ISTATE) + Cn(ISTATE,J)/
     1       RDVAL**nLR(ISTATE,J)*Dm
          END DO
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2 - DeOffset 
     1                     + Vlr(I,ISTATE) + VLIM(ISTATE) - Demo(ISTATE)
          IF (Rd .GT. 0) Vd = VPOT(I,ISTATE)
*          IF (Rd .GT. 0 .AND. out)  VdT = DE(ISTATE)*
*     &              (1.0d0-DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE))))**2
*          WRITE(24,802) RDVAL,Vlr(I,ISTATE),VPOT(I,ISTATE)-Vlr(I,ISTATE) 
*     &                  ,VPOT(I,ISTATE) 
*  802     FORMAT (F10.4,3F16.8)
*          WRITE(23,801) RDVAL,DE(ISTATE)*(1.0d0-XTEMP)**2 - DeOffset +
*     1      Vlr(I,ISTATE), Vlr(I,ISTATE), DE(ISTATE)*(1.0d0-XTEMP)**2-
*     2      DeOffset, DE(ISTATE)*(1.0d0-XTEMP)**2,VPOT(I,ISTATE)
*  801     FORMAT(F10.4,5F16.8)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
          BTEMP = 0.0d0 
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            DO 390 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J 
  390       CONTINUE
          ELSE
            DO 395 J = 1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  395       CONTINUE
            BTEMP = BTEMP*(RDVAL-RE(ISTATE))
          END IF 
c
c** Now to calculate the remaining partial derivatives.
c
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) /
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (1.0d0-XTEMP)
          DVDR(I,ISTATE) = -PBTEMP * (BMORSE(I,ISTATE)+PETEMP*BTEMP)
          DO 400 J = 0,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * (RDVAL-RE(ISTATE)) * RTEMP**J
  400     CONTINUE
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0 
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
  370   CONTINUE
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
******7***** End Expanded Morse Longe-Range Potential Function *********

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Modified Lennard-Jones-Gaussian potential function
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-7) THEN
        PSEL(ISTATE) = 7 
        BETA(0,ISTATE) = 2.0d0*AZERO-2.0d0*DBLE(NCN(ISTATE))
        BETA(1,ISTATE) = 2.0d0*AZERO*(AZERO+AONE+1.0d0)
        BETA(2,ISTATE) = AZERO*(4.0d0*(AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2.0d0-ATWO)*8.0d0/3.0d0+2.0d0)
     2                    -(2.0d0/3.0d0)*DBLE(NCN(ISTATE))
        DO 460 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  460   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.7) THEN
        DO 470 I = 1,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL 
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          BMORSE(I,ISTATE) = BETA(0,ISTATE)
          DO J = 1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
          END DO
c-----------------------------------------------------------------------
c** If chosen, constrain the Cn value to CNVAL by using a switching
c   function.
c
c             BETAsw(R) = fsw(R) [ BETA(R) - BINF ] + BINF
c
c                fsw(R) = 1 / [e^(ALPHA'*(R - Rx)) + 1]
c
          IF (CNVAL(ISTATE).GT.0.0d0) THEN
            BINF = DLOG(2.0d0*DE(ISTATE)*RE(ISTATE)**NCN(ISTATE)
     1             /CNVAL(ISTATE))
            FSW = 1.0d0 / (1.0d0+DEXP(ALPHA(ISTATE)*(RDVAL-RX(ISTATE))))
c
c** Update the generalized morse exponent.
c
            BMORSE(I,ISTATE) = FSW*(BMORSE(I,ISTATE) - BINF) + BINF
*          ENDIF
c**********************************
c
c** If chosen, constrain the Cn value to BVC by fixing the NBVCth BETA
c   parameter BETA(NBVC).
c
          ELSEIF (CNVAL(ISTATE).LT.0.0d0) THEN
            BINF = BETA(0,ISTATE)
            DO J = 1,NBETA(ISTATE)
              BINF = BINF + BETA(J,ISTATE)
            END DO
            BETA(NBETA(ISTATE)+1,ISTATE) = - BINF +
     1         DLOG(2.0d0*DE(ISTATE)
     2        *RE(ISTATE)**NCN(ISTATE)/DABS(CNVAL(ISTATE)))
c
c** Update the generalized morse exponent.
c
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)
     1        +BETA(NBETA(ISTATE)+1,ISTATE)*RTEMP**(NBETA(ISTATE)+1)
          ENDIF
c**********************************

c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            BMORSE(I,ISTATE) = BMORSE(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c** Calculate some temporary variables.
c
          XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     1            DEXP(-BMORSE(I,ISTATE)*RTEMP)
c
c** Now to calculate the actual potential
c
          Gauss(I,ISTATE) = Amp(ISTATE)*DEXP(-(RDVAL-R0(ISTATE))**2
     &                    /FWHM(ISTATE)**2)
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0-XTEMP)**2 + Gauss(I,ISTATE)
     1                   + VLIM(ISTATE) - DE(ISTATE)
c
c** Now to calculate temporary variables for the partial derivatives.
c
          BTEMP = 0.0d0
          DO 490 J = 1,NBETA(ISTATE)
            BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  490     CONTINUE
c**********************************
          IF (CNVAL(ISTATE).LT.0.0d0) BTEMP = BTEMP+
     1       BETA(NBETA(ISTATE)+1,ISTATE)*DBLE(NBETA(ISTATE)+1)
     2       *RTEMP**(NBETA(ISTATE)+1)
c**********************************
          PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (1.0d0-XTEMP)
c-----------------------------------------------------------------------
c** Now to calculate the partial derivatives.
c** Case when switching function constraint is used.
c
          IF (CNVAL(ISTATE).GT.0.0d0) THEN
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
     1                       + PBTEMP/DE(ISTATE)*(1-FSW)
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1       *(1-RTEMP+RTEMP*FSW) + PETEMP*(BMORSE(I,ISTATE)+FSW*BTEMP))
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * RTEMP**(J+1) * FSW
            END DO
c**********************************
c** Case when polynomial constraint is used.
c
          ELSEIF (CNVAL(ISTATE).LT.0.0d0) THEN
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
     1                     + PBTEMP/DE(ISTATE)*RTEMP**(NBETA(ISTATE)+2)
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1   *(1-RTEMP**(NBETA(ISTATE)+2))+PETEMP*(BMORSE(I,ISTATE)+BTEMP))
            DO J = 0,NBETA(ISTATE)
         DVDB(J,I,ISTATE)=PBTEMP*(RTEMP**(J+1)-RTEMP**(NBETA(ISTATE)+2))
            END DO
c**********************************
c
c** Case when no constraints are used.
c
          ELSE
            DVDD(I,ISTATE) = (1.0d0-XTEMP)**2 - 1.0d0
            DVDR(I,ISTATE) = - PBTEMP * (DBLE(NCN(ISTATE))/RE(ISTATE)
     1                       + PETEMP * (BMORSE(I,ISTATE) + BTEMP))
            DO J = 0,NBETA(ISTATE)
              DVDB(J,I,ISTATE) = PBTEMP * RTEMP**(J+1)
            END DO
          END IF
          dVdA(I,ISTATE) = Gauss(I,ISTATE)/Amp(ISTATE)
          dVdR0(I,ISTATE) = 2*(RDVAL-R0(ISTATE))/FWHM(ISTATE)**2 *
     &                      Gauss(I,ISTATE)  
          dVdW(I,ISTATE) = 2*(RDVAL-R0(ISTATE))**2/FWHM(ISTATE)**3 *
     &                      Gauss(I,ISTATE)
  470   CONTINUE
      END IF
******7***** End Modified Lennard-Jones-Gaussian potential function ****

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Total Expanded Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-8) THEN
        PSEL(ISTATE) = 8
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/RE(ISTATE)
        BETA(2,ISTATE) = AZERO*((2.0d0*AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2-ATWO)*4.0d0/3.0d0)/RE(ISTATE)
        DO 513 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  513   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.8) THEN
        DO 603 I = NDatPtS,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (Rd .GT. 0) RDVAL = Rd
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          bMorsIn(I,ISTATE) = BETA(0,ISTATE)
          DO 523 J = 1,nBetaIn(ISTATE)
            bMorsIn(I,ISTATE)=bMorsIn(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J 
  523     CONTINUE
          BMORSE(I,ISTATE) = bMorsIn(I,ISTATE)
          DO 533 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  533     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < RSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            bMorsIn(I,ISTATE) = bMorsIn(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
          Xtemp2= DEXP(-2.0d0*bMorsIn(I,ISTATE)*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0 + Xtemp2 - 2.0d0*XTEMP)
     &                   + VLIM(ISTATE) - DE(ISTATE)
          IF (Rd .GT. 0) Vd = VPOT(I,ISTATE)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = 1.0d0 + Xtemp2 - 2.0d0*XTEMP - 1.0d0
          BTEMP = 0.0d0
          Btemp2= 0.0d0
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            DO 543 J = 1,nBetaIn(ISTATE)
              Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  543       CONTINUE          
            BTEMP = Btemp2
            DO 553 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  553       CONTINUE
          ELSE
            DO 563 J = 1,nBetaIn(ISTATE)
              Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  563       CONTINUE       
            BTEMP = Btemp2
            Btemp2 = Btemp2*(RDVAL-RE(ISTATE))
            DO 573 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  573       CONTINUE
            BTEMP = BTEMP*(RDVAL-RE(ISTATE))
          END IF
c
c** Now to calculate the remaining partial derivatives.
c
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) /
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (RDVAL-RE(ISTATE))
          PBtemp2= 2.0d0 * DE(ISTATE) * Xtemp2* (RDVAL-RE(ISTATE))
          DVDR(I,ISTATE) = 2.0d0*DE(ISTATE)*(Xtemp2*(bMorsIn(I,ISTATE)
     &        + PETEMP*Btemp2) - XTEMP*(BMORSE(I,ISTATE)+PETEMP*BTEMP)) 
          DO 583  J = 0,nBetaIn(ISTATE)
            DVDB(J,I,ISTATE) = (-PBtemp2 + PBTEMP) * RTEMP**J
  583     CONTINUE
          DO 593  J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * RTEMP**J
  593     CONTINUE
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
  603   CONTINUE
      ENDIF
******7***** End Total Expanded Morse Potential Function ***************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Long-Range Total Expanded Morse Oscillator (trial version).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-9) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        PSEL(ISTATE) = 9
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/RE(ISTATE)
        BETA(2,ISTATE) = AZERO*((2.0d0*AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2-ATWO)*4.0d0/3.0d0)/RE(ISTATE)
        DO 613 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  613   CONTINUE
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.9) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        DO 703 I = NDatPtS,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (Rd .GT. 0) RDVAL = Rd
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          bMorsIn(I,ISTATE) = BETA(0,ISTATE)
          DO 623 J = 1,nBetaIn(ISTATE)
            bMorsIn(I,ISTATE)=bMorsIn(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  623     CONTINUE
          BMORSE(I,ISTATE) = bMorsIn(I,ISTATE)
          DO 633 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  633     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < RSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            bMorsIn(I,ISTATE) = bMorsIn(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
          Xtemp2= DEXP(-2.0d0*bMorsIn(I,ISTATE)*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
          Vlr(I,ISTATE) = 0.0d0
          DO J = 1,numCn(ISTATE)
            tempN = nLR(ISTATE,J)
            CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
            Dn(ISTATE,J,I) = Dm
            Vlr(I,ISTATE) = Vlr(I,ISTATE) + Cn(ISTATE,J)/
     &       RDVAL**nLR(ISTATE,J)*Dm
          END DO
          VPOT(I,ISTATE) = DE(ISTATE)*(1.0d0 + Xtemp2 - 2.0d0*XTEMP)
     &        - DeOffset + Vlr(I,ISTATE) + VLIM(ISTATE) - Demo(ISTATE)
          IF (Rd .GT. 0) Vd = VPOT(I,ISTATE)
*          IF (MOD(I-1,10) .EQ. 0) THEN 
*            WRITE (60,111) RDVAL,DE(ISTATE)*Xtemp2,
*     &        -DE(ISTATE)*2.0d0*XTEMP,DE(ISTATE)*(Xtemp2-*2.0d0*XTEMP)
*          END IF 
*  111     FORMAT (F10.4, 3G16.8)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = 1.0d0 + Xtemp2 - 2.0d0*XTEMP - 1.0d0
          BTEMP = 0.0d0
          Btemp2= 0.0d0
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            DO 643 J = 1,nBetaIn(ISTATE)
              Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  643       CONTINUE
            BTEMP = Btemp2
            DO 653 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**J
  653       CONTINUE
          ELSE
            DO 663 J = 1,nBetaIn(ISTATE)
              Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  663       CONTINUE
            BTEMP = Btemp2
            Btemp2 = Btemp2*(RDVAL-RE(ISTATE))
            DO 673 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  673       CONTINUE
            BTEMP = BTEMP*(RDVAL-RE(ISTATE))
          END IF
c
c** Now to calculate the remaining partial derivatives.
c
          IF (PoEpZ(ISTATE) .LT. 2) THEN
            PETEMP = 2.0d0*RDVAL/(RDVAL+RE(ISTATE))
          ELSE
            PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &        * RDVAL**nPB(ISTATE) /
     &        (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          END IF
          PBTEMP = 2.0d0 * DE(ISTATE) * XTEMP * (RDVAL-RE(ISTATE))
          PBtemp2= 2.0d0 * DE(ISTATE) * Xtemp2* (RDVAL-RE(ISTATE))
          DVDR(I,ISTATE) = 2.0d0*DE(ISTATE)*(Xtemp2*(bMorsIn(I,ISTATE)
     &        + PETEMP*Btemp2) - XTEMP*(BMORSE(I,ISTATE)+PETEMP*BTEMP))
          DO 683  J = 0,nBetaIn(ISTATE)
            DVDB(J,I,ISTATE) = (-PBtemp2 + PBTEMP) * RTEMP**J
  683     CONTINUE
          DO 693  J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            DVDB(J,I,ISTATE) = PBTEMP * RTEMP**J
  693     CONTINUE
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
  703   CONTINUE
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
******7* End Total Expanded Morse Long-Range Potential Function (Trial)*

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Double-Exponential Long-Range Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If selected, generate trial exponential parameters from low order
c   molecular constants.
c
      IF (PSEL(ISTATE).EQ.-10) THEN
        PSEL(ISTATE) = 10 
        BETA(0,ISTATE) = AZERO / RE(ISTATE)
        BETA(1,ISTATE) = AZERO*(AZERO+AONE)/RE(ISTATE)
        BETA(2,ISTATE) = AZERO*((2.0d0*AONE+1.0d0)*(AZERO+AONE)
     1                    +(AZERO**2-ATWO)*4.0d0/3.0d0)/RE(ISTATE)
        DO 713 I = 3,NBETA(ISTATE)
          BETA(I,ISTATE) = 0.0d0
  713   CONTINUE
      ENDIF
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.10) THEN
        DO 803 I = NDatPtS,NDATPT(ISTATE)
          RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
          R(I,ISTATE) = RDVAL
          IF (Rd .GT. 0) RDVAL = Rd
          IF (PoEpZ(ISTATE) .LT. 1) THEN
            WRITE (6,715)
  715       FORMAT (/,1X,16('  '),'WARNING!!!'/,1X,39('%%'),/' For Total
     & Expanded Morse Long-Range potential, the variable power (PoEpZ) m
     &ust'/,' be greater than 1, please try again'/,1X,39('%%'))
            GO TO 999
          ELSE
            RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &              (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          END IF
          bMorsIn(I,ISTATE) = BETA(0,ISTATE)
          DO 723 J = 1,nBetaIn(ISTATE)
            bMorsIn(I,ISTATE)=bMorsIn(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  723     CONTINUE
          BMORSE(I,ISTATE) = bMorsIn(I,ISTATE)
          DO 733 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            BMORSE(I,ISTATE) = BMORSE(I,ISTATE)+BETA(J,ISTATE)*RTEMP**J
  733     CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < RSR.
c
          IF (I.LT.IRSR(ISTATE)) THEN
            bMorsIn(I,ISTATE) = bMorsIn(IRSR(ISTATE),ISTATE)
          END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
          IF (inAndOut(ISTATE) .EQ. 0) THEN
            XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
            Xtemp2= DEXP(-2.0d0*bMorsIn(I,ISTATE)*(RDVAL-RE(ISTATE)))
            BMtemp = BMORSE(I,ISTATE)
            BMtemp2 = bMorsIn(I,ISTATE)
          ELSE
            IF (RDVAL .LT. RE(ISTATE)) THEN
              XTEMP = DEXP(-bMorsIn(I,ISTATE)*(RDVAL-RE(ISTATE)))
              Xtemp2= DEXP(-2.0d0*bMorsIn(I,ISTATE)*(RDVAL-RE(ISTATE)))
              BMtemp = bMorsIn(I,ISTATE)
              BMtemp2 = bMorsIn(I,ISTATE)
            ELSE
              XTEMP = DEXP(-BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
              Xtemp2= DEXP(-2.0d0*BMORSE(I,ISTATE)*(RDVAL-RE(ISTATE)))
              BMtemp = BMORSE(I,ISTATE)
              BMtemp2 = BMORSE(I,ISTATE)
            END IF
          END IF
            
c
c** Now to calculate the actual potential, (see notes for the math)
c
          Vlr(I,ISTATE) = 0.0d0
          A(ISTATE) = 0.0d0
          B(ISTATE) = 0.0d0
          dAdRe = 0.0d0
          dBdRe = 0.0d0
          dVdBtemp = 0.0d0
          DO J = 1,numCn(ISTATE)
            tempN = nLR(ISTATE,J)
            CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
            Dn(ISTATE,J,I) = Dm
            Vlr(I,ISTATE) = Vlr(I,ISTATE) + Cn(ISTATE,J)/
     &        RDVAL**nLR(ISTATE,J)*Dm
**$            VlrB = Vlr(I,ISTATE)
            CALL dampF(RE(ISTATE),dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
            A(ISTATE)=A(ISTATE)+Cn(ISTATE,J)/RE(ISTATE)**nLR(ISTATE,J)
     &                *Dm*(1.0d0+Dmp/(BETA(0,ISTATE)*Dm)-nLR(ISTATE,J)
     &                /(BETA(0,ISTATE)*RE(ISTATE)))
            B(ISTATE)=B(ISTATE)+Cn(ISTATE,J)/RE(ISTATE)**nLR(ISTATE,J)
     &                *Dm*(2.0d0+Dmp/(BETA(0,ISTATE)*Dm)-nLR(ISTATE,J)
     &                /(BETA(0,ISTATE)*RE(ISTATE)))
            dAdRe=dAdRe+Cn(ISTATE,J)*((Dmp*RE(ISTATE)-nLR(ISTATE,J)*Dm)
     &        /RE(ISTATE)**(nLR(ISTATE,J)+1)*(1.0d0+Dmp/(BETA(0,ISTATE)
     &        *Dm)-nLR(ISTATE,J)/(BETA(0,ISTATE)*RE(ISTATE))) 
     &        + Dm/RE(ISTATE)**nLR(ISTATE,J)*((Dm*Dmpp-Dmp**2)
     &        /(BETA(0,ISTATE)*Dm**2)+nLR(ISTATE,J)/(BETA(0,ISTATE)
     &        *RE(ISTATE)**2)))
            dBdRe=dBdRe+Cn(ISTATE,J)*((Dmp*RE(ISTATE)-nLR(ISTATE,J)*Dm)
     &        /RE(ISTATE)**(nLR(ISTATE,J)+1)*(2.0d0+Dmp/(BETA(0,ISTATE)
     &        *Dm)-nLR(ISTATE,J)/(BETA(0,ISTATE)*RE(ISTATE)))   
     &        + Dm/RE(ISTATE)**nLR(ISTATE,J)*((Dm*Dmpp-Dmp**2)
     &        /(BETA(0,ISTATE)*Dm**2)+nLR(ISTATE,J)/(BETA(0,ISTATE)
     &        *RE(ISTATE)**2)))
            dVdBtemp=dVdBtemp+Cn(ISTATE,J)*Dm/(RE(ISTATE)**nLR(ISTATE,J)
     &        *BETA(0,ISTATE)**2)*(nLR(ISTATE,J)/RE(ISTATE)-Dmp/Dm)  
          END DO
          fd(I,ISTATE) = 1.0d0
          IF (Df .EQ. 4) THEN
            fd(I,ISTATE) = 1.0d0/(1.0d0 + DEXP(-1.3*dampBeta(ISTATE)*
     &        (RDVAL-ReBar)))
**$      &        (RDVAL-RE(ISTATE))))
            Vlr(I,ISTATE) = fd(I,ISTATE)*Vlr(I,ISTATE)
          END IF
          dAdRe = dAdRe * fd(I,ISTATE)
          dBdRe = dBdRe * fd(I,ISTATE)
          dVdBtemp = dVdBtemp * fd(I,ISTATE)
          cTemp = C(ISTATE)*DEXP(-cRp(ISTATE)*BETA(0,ISTATE) 
     &            * (RDVAL-RE(ISTATE)))
          Vlr(I,ISTATE) = Vlr(I,ISTATE) + cTemp
          A(ISTATE) = fd(I,ISTATE)*A(ISTATE)+DE(ISTATE) 
     &                + (1.0d0-cRp(ISTATE))*C(ISTATE)
          B(ISTATE) = fd(I,ISTATE)*B(ISTATE) + 2.0d0*DE(ISTATE)
     &                + (2.0d0-cRp(ISTATE))*C(ISTATE)
          VPOT(I,ISTATE) = A(ISTATE)*Xtemp2 - B(ISTATE)*XTEMP+DE(ISTATE)
     &                     + Vlr(I,ISTATE) + VLIM(ISTATE) - Demo(ISTATE)
          IF (Rd .GT. 0) Vd = VPOT(I,ISTATE)
*          DEP = A(ISTATE)*Xtemp2 - B(ISTATE)*XTEMP+DE(ISTATE)
*          DELR = DEP + Vlr(I,ISTATE)
*          fd(I,ISTATE) = 1.0d0/(1.0d0 + DEXP(-1.3*dampBeta(ISTATE)*
*     &      (RDVAL-2.9)))
*          IF (MOD(I-1,10) .EQ. 0) THEN 
*            WRITE (44,456) RDVAL,Vlr(I,ISTATE),DEP,DELR,
*     &            fd(I,ISTATE)*Vlr(I,ISTATE)
*          END IF
*  456     FORMAT (4x,5G14.5)
c
c** Now to calculate the partial derivatives
c
          DVDD(I,ISTATE) = 1.0d0 + Xtemp2 - 2.0d0*XTEMP - 1.0d0
          BTEMP = 0.0d0
          Btemp2= 0.0d0
          IF (inAndOut(ISTATE) .EQ. 0) THEN
            DO 763 J = 1,nBetaIn(ISTATE)
              Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  763       CONTINUE
            BTEMP = Btemp2
            DO 773 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
  773       CONTINUE
          ELSE
            IF (RDVAL .LT. RE(ISTATE)) THEN
              DO J = 1,nBetaIn(ISTATE)
                Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
              END DO
              BTEMP = Btemp2
            ELSE
              DO J = 1,NBETA(ISTATE)
                Btemp2= Btemp2+ BETA(J,ISTATE) * DBLE(J) * RTEMP**(J-1)
              END DO
              BTEMP = Btemp2
            END IF
          END IF  
          Btemp2 = Btemp2*(RDVAL-RE(ISTATE))
          BTEMP = BTEMP*(RDVAL-RE(ISTATE))
c
c** Now to calculate the remaining partial derivatives.
c
          PETEMP = 2.0d0 * nPB(ISTATE) * RE(ISTATE)**(nPB(ISTATE)-1)
     &      * RDVAL**nPB(ISTATE) /
     &      (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
          PBTEMP = B(ISTATE) * XTEMP * (RDVAL-RE(ISTATE))
          PBtemp2= 2.0d0 * A(ISTATE) * Xtemp2* (RDVAL-RE(ISTATE))
          DVDR(I,ISTATE)=2.0d0*A(ISTATE)*Xtemp2*(BMtemp2+PETEMP*Btemp2)
     &      + Xtemp2*dAdRe - B(ISTATE)*XTEMP*(BMtemp+PETEMP*BTEMP)
     &      - XTEMP*dBdRe + cRp(ISTATE)*BETA(0,ISTATE)*cTemp
**$          IF (Df .EQ. 4) THEN
**$            DVDR(I,ISTATE) = DVDR(I,ISTATE) - 1.3*dampBeta(ISTATE) * 
**$     &        DEXP(-1.3*dampBeta(ISTATE)*(RDVAL-RE(ISTATE))) * 
**$     &        fd(I,ISTATE)**2 * VlrB 
**$          END IF
          DVDB(0,I,ISTATE) = -PBtemp2 + PBTEMP + dVdBtemp*(Xtemp2-XTEMP)
     &                       - cRp(ISTATE)*(RDVAL - RE(ISTATE)) * cTemp 
          DO 783  J = 1,nBetaIn(ISTATE)
            DVDB(J,I,ISTATE) = (-PBtemp2 + PBTEMP) * RTEMP**J
  783     CONTINUE
          DO 793  J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
            IF (inAndOut(ISTATE) .EQ. 0) THEN
              DVDB(J,I,ISTATE) = PBTEMP * RTEMP**J
            ELSE
              IF (RDVAL .LT. RE(ISTATE)) THEN
                DVDB(J,I,ISTATE) = 0.0d0
              ELSE
                DVDB(J,I,ISTATE) = PBTEMP * RTEMP**J
              END IF
            END IF
  793     CONTINUE
          DVDC(I,ISTATE) = DEXP(-cRp(ISTATE)*BETA(0,ISTATE) * 
     &      (RDVAL-RE(ISTATE))) + (1.0d0-cRp(ISTATE))*Xtemp2 - 
     &      (2.0d0-cRp(ISTATE))*XTEMP    
c *** DBDR and DBDB is used in uncertainty calculation, see fununc.f
          DBDR(I,ISTATE) = -PETEMP * BTEMP / (RDVAL-RE(ISTATE))
          DBDB(0,I,ISTATE) = 1.0d0
          DO J = 1,NBETA(ISTATE)
            DBDB(J,I,ISTATE) = DBDB(J-1,I,ISTATE) * RTEMP
          END DO
  803   CONTINUE
      ENDIF
******7* End Total Expanded Morse Long-Range Potential Function ********

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Born-Oppenheimer breakdown correction terms
c** and the lambda doubling correction terms.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (PSEL(ISTATE) .EQ. 6 .OR. PSEL(ISTATE) .EQ. 9) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
      ENDIF
      DO 599 I = NDatPtS,NDATPT(ISTATE)
        RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
        IF (Rd .GT. 0) RDVAL = Rd
        RM2 = 1/RDVAL**2
        RM4 = 1/RDVAL**4
c
c** Choose Born-Oppenheimer correction expansion parameter.
c
        IF (PoEpZ(ISTATE).LE.0) THEN
          IF (PoEpZ(ISTATE).EQ.0) RTEMP = (RDVAL-RE(ISTATE))
     1                                 /(RDVAL+RE(ISTATE))
          IF (PoEpZ(ISTATE).EQ.-1) RTEMP = (RDVAL-RE(ISTATE))
c
c** Adiabatic correction terms and state dependent partial derivatives.
c
          UAR(I,ISTATE) = 0.0d0
          DUADR(I,ISTATE) = 0.0d0
          DO 250 J = 0,NUA(ISTATE)
            UAR(I,ISTATE) = UAR(I,ISTATE) + UA(J,ISTATE)*RTEMP**J
            DVDUA(J,I,ISTATE) = RTEMP**J
            DUADR(I,ISTATE) = DUADR(I,ISTATE) 
     1                       - DBLE(J)*UA(J,ISTATE)*RTEMP**(J-1)
  250     CONTINUE
          IF (PoEpZ(ISTATE).EQ.0) DUADR(I,ISTATE) = DUADR(I,ISTATE)
     1      *2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          IF (Rd .GT. 0) Uad = UAR(I,ISTATE) 

          UBR(I,ISTATE) = 0.0d0
          DUBDR(I,ISTATE) = 0.0d0
          DO 260 J = 0,NUB(ISTATE)
            UBR(I,ISTATE) = UBR(I,ISTATE) + UB(J,ISTATE)*RTEMP**J
            DVDUB(J,I,ISTATE) = RTEMP**J
            DUBDR(I,ISTATE) = DUBDR(I,ISTATE)
     1                       - DBLE(J)*UB(J,ISTATE)*RTEMP**(J-1)
  260     CONTINUE
          IF (PoEpZ(ISTATE).EQ.0) DUBDR(I,ISTATE) = DUBDR(I,ISTATE) 
     1      *2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          IF (Rd .GT. 0) Ubd = UBR(I,ISTATE)
c
c** Rotational non-adiabatic correction terms and state dependent
c   partial derivatives.
c
          QAR(I,ISTATE) = 0.0d0
          DQADR(I,ISTATE) = 0.0d0
          DO 270 J = 0,NQA(ISTATE)
            QAR(I,ISTATE) = QAR(I,ISTATE) + QA(J,ISTATE) * RTEMP**J
            DVDQA(J,I,ISTATE) = RTEMP**J * RM2
            DQADR(I,ISTATE) = DQADR(I,ISTATE)
     1                       - DBLE(J)*QA(J,ISTATE)*RTEMP**(J-1)
  270     CONTINUE
          IF (PoEpZ(ISTATE).EQ.0) DQADR(I,ISTATE) = DQADR(I,ISTATE) 
     1      *2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          DQADR(I,ISTATE) = DQADR(I,ISTATE) * RM2
          IF (Rd .GT. 0) Qad = QAR(I,ISTATE) 
         
          QBR(I,ISTATE) = 0.0d0
          DQBDR(I,ISTATE) = 0.0d0
          DO 280 J = 0,NQB(ISTATE)
            QBR(I,ISTATE) = QBR(I,ISTATE) + QB(J,ISTATE) * RTEMP**J
            DVDQB(J,I,ISTATE) = RTEMP**J * RM2
            DQBDR(I,ISTATE) = DQBDR(I,ISTATE)
     1                       - DBLE(J)*QB(J,ISTATE)*RTEMP**(J-1)
  280     CONTINUE
          IF (PoEpZ(ISTATE).EQ.0) DQBDR(I,ISTATE) = DQBDR(I,ISTATE) 
     1      *2.0d0*RDVAL/(RDVAL+RE(ISTATE))**2
          DQBDR(I,ISTATE) = DQBDR(I,ISTATE) * RM2
          IF (Rd .GT. 0) Qbd = QBR(I,ISTATE) 
          
c** For lambda doubling
          lambdaR(I,ISTATE) = 0.0d0
          dLdRe(I,ISTATE) = 0.0d0
          DO 300 J = 0,nLamdCoef(ISTATE)
            lambdaR(I,ISTATE) = lambdaR(I,ISTATE) + 
     &        lamdCoef(J,ISTATE)/(RDVAL+RE(ISTATE))**J 
            dVdL(J,I,ISTATE) = 1/(RDVAL+RE(ISTATE))**J * RM4
            dLdRe(I,ISTATE) = dLdRe(I,ISTATE) - DBLE(J)*
     &        lamdCoef(J,istate)/(RDVAL+RE(ISTATE))**(J+1)
  300     CONTINUE
          dLdRe(I,ISTATE) = dLdRe(I,ISTATE) * RM4
          IF (Rd .GT. 0) Ld=lambdaR(I,ISTATE)
******7***************************************************************72
**  for PoEpZ > 0, introduce a new expression for BOB and lambda 
**  doubling functions
******7***************************************************************72
        ELSE IF (PoEpZ(ISTATE).GT.0) THEN
          RTEMP = 1.0d0/(RDVAL**nPu(ISTATE)+RE(ISTATE)**nPu(ISTATE))
          ZTEMP = RDVAL**nPu(ISTATE)-RE(ISTATE)**nPu(ISTATE)
          UAR(I,ISTATE) = UA(0,ISTATE)*2*RE(ISTATE)**nPu(ISTATE)*RTEMP
     &      + UA(NUA(ISTATE),ISTATE)*ZTEMP*RTEMP
          DVDUA(0,I,ISTATE) = 2*RE(ISTATE)**nPu(ISTATE)*RTEMP
          DVDUA(NUA(ISTATE),I,ISTATE) = ZTEMP*RTEMP          
          DUADR(I,ISTATE)=(UA(0,ISTATE)-UA(NUA(ISTATE),ISTATE))*RTEMP**2
          DO 500 J = 1,NUA(ISTATE)-1
            UAR(I,ISTATE) = UAR(I,ISTATE) + UA(J,ISTATE)
     &        *2*RE(ISTATE)**nPu(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
            DVDUA(J,I,ISTATE) = 2*RE(ISTATE)**nPu(ISTATE)*ZTEMP**J
     &        *RTEMP**(J+1)
            DUADR(I,ISTATE) = DUADR(I,ISTATE)+UA(J,ISTATE)*ZTEMP**(J-1) 
     &        *(ZTEMP-2*J*RE(ISTATE)**nPu(ISTATE))*RTEMP**(J+2)
  500     CONTINUE
          DUADR(I,ISTATE) = DUADR(I,ISTATE)*2
     &      *nPu(ISTATE)*RE(ISTATE)**(nPu(ISTATE)-1)*RDVAL**nPu(ISTATE)
          IF (Rd .GT. 0) Uad = UAR(I,ISTATE)

          UBR(I,ISTATE) = UB(0,ISTATE)*2*RE(ISTATE)**nPu(ISTATE)*RTEMP
     &      + UB(NUB(ISTATE),ISTATE)*ZTEMP*RTEMP
          DVDUB(0,I,ISTATE) = 2*RE(ISTATE)**nPu(ISTATE)*RTEMP
          DVDUB(NUB(ISTATE),I,ISTATE) = ZTEMP*RTEMP          
          DUBDR(I,ISTATE)=(UB(0,ISTATE)-UB(NUB(ISTATE),ISTATE))*RTEMP**2
          DO 520 J = 1,NUB(ISTATE)-1
            UBR(I,ISTATE) = UBR(I,ISTATE) + UB(J,ISTATE)
     &        *2*RE(ISTATE)**nPu(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
            DVDUB(J,I,ISTATE) = 2*RE(ISTATE)**nPu(ISTATE)*ZTEMP**J
     &        *RTEMP**(J+1)
            DUBDR(I,ISTATE) = DUBDR(I,ISTATE)+UB(J,ISTATE)*ZTEMP**(J-1) 
     &        *(ZTEMP-2*J*RE(ISTATE)**nPu(ISTATE))*RTEMP**(J+2)
  520     CONTINUE
          DUBDR(I,ISTATE) = DUBDR(I,ISTATE)*2
     &      *nPu(ISTATE)*RE(ISTATE)**(nPu(ISTATE)-1)*RDVAL**nPu(ISTATE)
          IF (Rd .GT. 0) Ubd = UBR(I,ISTATE)
**
** Note: the 1/R^2 factor in the Hamiltonian has been incorporated into 
**       the partial derivatives
**
          RTEMP = 1.0d0/(RDVAL**nPq(ISTATE)+RE(ISTATE)**nPq(ISTATE))
          ZTEMP = RDVAL**nPq(ISTATE)-RE(ISTATE)**nPq(ISTATE)
          QAR(I,ISTATE) = QA(0,ISTATE)*2*RE(ISTATE)**nPq(ISTATE)*RTEMP
          DVDQA(0,I,ISTATE) = 2*RE(ISTATE)**nPq(ISTATE)*RTEMP*RM2
          DQADR(I,ISTATE) = UB(0,ISTATE)*RTEMP**2
          DO 540 J = 1,NQA(ISTATE)
            QAR(I,ISTATE) = QAR(I,ISTATE) + QA(J,ISTATE)
     &        *2*RE(ISTATE)**nPq(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
            DVDQA(J,I,ISTATE) = 2*RE(ISTATE)**nPq(ISTATE)*ZTEMP**J
     &        *RTEMP**(J+1)*RM2
            DQADR(I,ISTATE) = DQADR(I,ISTATE)+QA(J,ISTATE)*ZTEMP**(J-1)
     &        *(ZTEMP-2*J*RE(ISTATE)**nPq(ISTATE))*RTEMP**(J+2)
  540     CONTINUE
          DQADR(I,ISTATE) = DQADR(I,ISTATE)*RM2 * 2
     &      *nPq(ISTATE)*RE(ISTATE)**(nPq(ISTATE)-1)*RDVAL**nPq(ISTATE)
          IF (Rd .GT. 0) Qad = QAR(I,ISTATE)

          QBR(I,ISTATE) = QB(0,ISTATE)*2*RE(ISTATE)**nPq(ISTATE)*RTEMP
          DVDQB(0,I,ISTATE) = 2*RE(ISTATE)**nPq(ISTATE)*RTEMP*RM2
          DQBDR(I,ISTATE) = UB(0,ISTATE)*RTEMP**2
          DO 560 J = 1,NQB(ISTATE)
            QBR(I,ISTATE) = QBR(I,ISTATE) + QB(J,ISTATE)
     &        *2*RE(ISTATE)**nPq(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
            DVDQB(J,I,ISTATE) = 2*RE(ISTATE)**nPq(ISTATE)*ZTEMP**J
     &        *RTEMP**(J+1)*RM2
            DQBDR(I,ISTATE) = DQBDR(I,ISTATE)+QB(J,ISTATE)*ZTEMP**(J-1)
     &        *(ZTEMP-2*J*RE(ISTATE)**nPq(ISTATE))*RTEMP**(J+2)
  560     CONTINUE
          DQBDR(I,ISTATE) = DQBDR(I,ISTATE)*RM2 * 2
     &      *nPq(ISTATE)*RE(ISTATE)**(nPq(ISTATE)-1)*RDVAL**nPq(ISTATE)
          IF (Rd .GT. 0) Qbd = QBR(I,ISTATE)

c** For   lambda doubling
          RTEMP = 1.0d0/(RDVAL**nPL(ISTATE)+RE(ISTATE)**nPL(ISTATE))
          ZTEMP = RDVAL**nPL(ISTATE)-RE(ISTATE)**nPL(ISTATE)
          lambdaR(I,ISTATE)=lamdCoef(0,ISTATE)*2*RE(ISTATE)**nPL(ISTATE)
     &      *RTEMP + lamdCoef(nLamdCoef(ISTATE),ISTATE)*ZTEMP*RTEMP
          dVdL(0,I,ISTATE) = 2*RE(ISTATE)**nPL(ISTATE)*RTEMP * RM4
          dVdL(nLamdCoef(ISTATE),I,ISTATE) = ZTEMP*RTEMP * RM4
          dLdRe(I,ISTATE) = (lamdCoef(0,ISTATE) - 
     &      lamdCoef(nLamdCoef(ISTATE),ISTATE))*RTEMP**2
          DO 590 J = 1,nLamdCoef(ISTATE)-1
            lambdaR(I,ISTATE) = lambdaR(I,ISTATE) + lamdCoef(J,ISTATE)
     &        *2*RE(ISTATE)**nPL(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
            dVdL(J,I,ISTATE) = 2*RE(ISTATE)**nPL(ISTATE)*ZTEMP**J
     &        *RTEMP**(J+1) * RM4
            dLdRe(I,ISTATE) = dLdRe(I,ISTATE) + lamdCoef(J,ISTATE) * 
     &        ZTEMP**(J-1)*(ZTEMP-2*J*RE(ISTATE)**nPL(ISTATE)) * 
     &        RTEMP**(J+2)
  590     CONTINUE
          dLdRe(I,ISTATE) = dLdRe(I,ISTATE)*2*nPL(ISTATE) * 
     &      RE(ISTATE)**(nPL(ISTATE)-1)*RDVAL**nPL(ISTATE) * RM4

************ if the simple y variable is used **************************
c          DO 590 J = 0,nLamdCoef(ISTATE)
c            lambdaR(I,ISTATE) = lambdaR(I,ISTATE) + lamdCoef(J,ISTATE)
c     &                          * (ZTEMP*RTEMP)**J
c            dVdL(J,I,ISTATE) = (ZTEMP*RTEMP)**J * RM4
c            dLdRe(I,ISTATE)=dLdRe(I,ISTATE)-DBLE(2*J)*lamdCoef(J,ISTATE)
c     &       *nPL(ISTATE)*RE(ISTATE)**(nPL(ISTATE)-1)*RDVAL**nPL(ISTATE)              
c     &       *ZTEMP**(J-1)*RTEMP**(J+1)
c  590     CONTINUE
c          dLdRe(I,ISTATE) = dLdRe(I,ISTATE) * RM4
************************************************************************

          IF (Rd .GT. 0) Ld=lambdaR(I,ISTATE)
        END IF
  599 CONTINUE
      IF (PSEL(ISTATE) .EQ. 6 .OR. PSEL(ISTATE) .EQ. 9) THEN
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
******7********************** BLOCK END ******************************72

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Print the potential within two dissociation energies from the
c   asymptotic limit. (values outside of this range can be discarded)
c
c                   VLIM - 2*De < V(I) < VLIM + 2*De
c
c   Also use a radial mesh of approximately 0.05 angs.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (outputPrint) THEN
        IF (PSEL(ISTATE).EQ.4 .AND. CNVAL(ISTATE).NE.0.0d0) THEN
           WRITE(6,905) BetaInfinity
        END IF
        IF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).LT.0.0d0) THEN
          WRITE(6,888) 'Be','ta',NBETA(ISTATE)+1,
     &                                      BETA(NBETA(ISTATE)+1,ISTATE)
        END IF
        IF (PSEL(ISTATE) .GE. 8) THEN 
          WRITE (6,907) A(ISTATE),B(ISTATE),C(ISTATE)
          IF (inAndOut(ISTATE) .EQ. 0) THEN
            WRITE (6,909) nBetaIn(ISTATE),NBETA(ISTATE) 
          END IF
        END IF
        IF (inAndOut(ISTATE) .EQ. 1) THEN
          WRITE (6,911) nBetaIn(ISTATE),NBETA(ISTATE)
        END IF
        IF (ABS(PSEL(ISTATE)).EQ.6 .OR. ABS(PSEL(ISTATE)).GE.9) THEN
          step = 0.05 / RH(ISTATE)
          DO J = 1,numCn(ISTATE)
            WRITE (7,877) nLR(ISTATE,J)
            DO I = 1,NDATPT(ISTATE),step
              RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
              WRITE (7,879) RDVAL, Dn(ISTATE,J,I)
            END DO
          END DO
        END IF 
  877 FORMAT (' n =',I3)
  879 FORMAT (4x,2F10.5)
*        step = 0.05 / RH(ISTATE)
*        DO I = 1,NDATPT(ISTATE),step
*          IF(DABS(VPOT(I,ISTATE)-VLIM(ISTATE)).LT.2.0d0*DE(ISTATE)) THEN
*            RDVAL = RMIN(ISTATE) + RH(ISTATE)*DBLE(I-1)
*            WRITE(10,900) RDVAL,VPOT(I,ISTATE)
*            WRITE(11,900) RDVAL,BMORSE(I,ISTATE)
*            WRITE(12,900) RDVAL,UAR(I,ISTATE)
*            WRITE(13,900) RDVAL,UBR(I,ISTATE)
*            WRITE(14,900) RDVAL,QAR(I,ISTATE)
*            WRITE(15,900) RDVAL,QBR(I,ISTATE)
*          END IF
*        END DO
      END IF
  888 FORMAT (3X,A2,A2,'(',I2,')',1PD19.10,x,/,x,'defined by the beta in
     &finity minus the sum of all the other beta coefficients')
  900 FORMAT(5X,G18.8,5X,G18.8)
  905 FORMAT(1x,78('=')/,/' Beta function limit is ',G18.8,' as R goes t
     &o infinity')
  907 FORMAT(1x,78('=')/,5X,'A = ',G16.8,/,5X,'B = ',G16.8,/,5X,'C = ',
     &G16.8)
  909 FORMAT(1x,'Inner beta function is summed over 0 to',I3,/,
     &      1x,'Outer beta function is summed over 0 to',I3)
  911 FORMAT(1x,'Beta function is summed over 0 to ',I2,' at R < Re',/,
     &       1x,'Beta function is summed over 0 to ',I2,' at R >= Re') 
******7********************** BLOCK END ******************************72

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to determine the relevant partial derivatives for the parameters
c   that will be fitted.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      l = 0
      DO i = 1, ISTATE - 1
        l = l + NPPS(i)
      END DO
      IPV = l
      IF (PSEL(ISTATE).NE.5) THEN
        IF (IFDE(ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = DVDD(J,ISTATE)
            DBDP(J,IPV) = 0.0d0
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      ENDIF
      IF (IFRE(ISTATE).LE.0) THEN
        IPV = IPV + 1
        DO J = NDatPtS,NDATPT(ISTATE)
          DVBODP(J,IPV) = DVDR(J,ISTATE) 
          DBDP(J,IPV) = DBDR(J,ISTATE)
          DUADP(J,IPV) = DUADR(J,ISTATE)
          DUBDP(J,IPV) = DUBDR(J,ISTATE)
          DQADP(J,IPV) = DQADR(J,ISTATE)
          DQBDP(J,IPV) = DQBDR(J,ISTATE)
          dLdP(J,IPV) = dLdRe(J,ISTATE)
        END DO
      ENDIF
      IF (PSEL(ISTATE).EQ.7) THEN
        IF (ifAmp(ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = dVdA(J,ISTATE)
            DBDP(J,IPV) = 0.0d0
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
        IF (ifR0(ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = dVdR0(J,ISTATE)
            DBDP(J,IPV) = 0.0d0
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
        IF (ifFWHM(ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = dVdW(J,ISTATE)
            DBDP(J,IPV) = 0.0d0
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      ENDIF
      DO I = 0,NBETA(ISTATE)
        IF (IFBETA(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = DVDB(I,J,ISTATE)
            DBDP(J,IPV) = DBDB(I,J,ISTATE)
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END DO
      DO I = 0,NUA(ISTATE)
        IF (IFUA(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = 0.0d0
            DBDP(J,IPV) = 0.0d0 
            DUADP(J,IPV) = DVDUA(I,J,ISTATE)
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END DO
      DO I = 0,NUB(ISTATE)
        IF (IFUB(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = 0.0d0
            DBDP(J,IPV) = 0.0d0 
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = DVDUB(I,J,ISTATE)
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END DO
      DO I = 0,NQA(ISTATE)
        IF (IFQA(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = 0.0d0
            DBDP(J,IPV) = 0.0d0 
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = DVDQA(I,J,ISTATE)
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END DO
      DO I = 0,NQB(ISTATE)
        IF (IFQB(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = 0.0d0
            DBDP(J,IPV) = 0.0d0 
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = DVDQB(I,J,ISTATE)
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END DO
      DO I = 0,nLamdCoef(ISTATE)
        IF (ifLamdCoef(I,ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = 0.0d0
            DBDP(J,IPV) = 0.0d0 
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = dVdL(I,J,ISTATE)
          END DO
        ENDIF
      END DO
      IF (PSEL(ISTATE).EQ.10) THEN
        IF (ifC(ISTATE).LE.0) THEN
          IPV = IPV + 1
          DO J = NDatPtS,NDATPT(ISTATE)
            DVBODP(J,IPV) = DVDC(J,ISTATE)
            DBDP(J,IPV) = 0.0d0
            DUADP(J,IPV) = 0.0d0
            DUBDP(J,IPV) = 0.0d0
            DQADP(J,IPV) = 0.0d0
            DQBDP(J,IPV) = 0.0d0
            dLdP(J,IPV) = 0.0d0
          END DO
        ENDIF
      END IF
******7********************** BLOCK END ******************************72

******7***************************************************************72
c** Now calculate the potential value at Rd when distance Rd > 0 and 
c   partial derivative dVdPk at Rd
c   Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabaric BOB correction terms.
c
cV(R) = (2*mu/hbar^2)*(VPOT(R)+U(R))+(J(J+1)-OMEGA^2)(1+q(R))/R^2
c       ±(hbar^2/2*mu)(J(J+1))f(R)/R^4
c     = BFCT*(VPOT+U(R))+(JROT(JROT+1)-IOMEGA**2)(1+q(R))*RH**2/R^2
c       ±(hbar^2/2*mu)(JROT(JROT+1))f(R)*RH**2/R^4
******7***************************************************************72
      IF (Rd .GT. 0) THEN 
        NBAND = IB(DATN)
        IISTP = ISTP(NBAND)
*        WRITE (40,644) IISTP,Rd,RDVAL,Vd,I,NDATPT(ISTATE)
*  644   FORMAT ('IISTP =',I3,' Rd =',G16.8,' RDVAL =',G16.8,
*     &          ' Vd =',G16.8,' I =',I6,' NDATPT =',I6)
        JFCT = 16.85762908d0/ZMASS(3,IISTP)
     &       * DBLE(JPP(DATN)*(JPP(DATN)+1)-IOMEG(ISTATE)*IOMEG(ISTATE))
        JfL=(16.85762908d0/ZMASS(3,IISTP))**2 
     &       * DBLE(JPP(DATN)*(JPP(DATN)+1))
        Vd = Vd + ZMUA(IISTP)*Uad + ZMUB(IISTP)*Ubd      
*        WRITE (40,646) JFCT,Rd,Vd
*  646   FORMAT ('JFCT =',G16.8,' Rd =',G16.8,' Vd =',G16.8)
c** For J > 0, also include the centrifugal distortion and lambda doubling.
        IF (ISTATE .EQ. 1) THEN
          efParity = EFPP(DATN)
        ELSEIF (ISTATE .EQ. 2) THEN
          efParity = EFP(DATN)
        ENDIF
        IF (JPP(DATN).GT.0) THEN
          Vd = Vd + JFCT*(1.0d0+ZMQA(IISTP)*Qad+ZMQB(IISTP)*Qbd)/Rd**2 
     &                                     + efParity * JfL * Ld/Rd**4
        END IF
*        WRITE (40,648) JPP(DATN),efParity,Rd,Vd
*  648   FORMAT ('J =',I3,' efParity =',I3,' Rd =',G16.8,' Vd =',G16.8/)

        J = NDATPT(ISTATE)
        DO IPV = l+1,l+NPPS(ISTATE)
          dVdPk(IPV) = DVBODP(J,IPV) + ZMUA(IISTP)*DUADP(J,IPV)
     &      + ZMUB(IISTP)*DUBDP(J,IPV) + JFCT*(ZMQA(IISTP)*DQADP(J,IPV)
     &      + ZMQB(IISTP)*DQBDP(J,IPV)) + efParity * JfL * dLdP(J,IPV)
        END DO
      END IF
******7********************** BLOCK END ******************************72
      IF (Rd .GT. 0) NDATPT(ISTATE) = NDATPT(ISTATE) - 1
******7***************************************************************72

      RETURN
 999  END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
