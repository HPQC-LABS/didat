c** Version 0.9 dated Aug. 8, 2000
c***********************************************************************
        PROGRAM DSPotFit
c***********************************************************************
c** Program "D(iatomic)S(inglet)Pot(ential)Fit" (DSPotFit) for
c   performing least squares fits of diatomic spectral data to
c   molecular potential functions for multiple electronic states.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++         COPYRIGHT 1997  by  Jenning Seto and R. J. Le Roy        +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c++ Uses least-squares subroutine NLLSSRR written by Dulick & Le Roy +++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
c** This program can perform the following types of calculations:
c (i)  From a set of read-in constants, make predictions for any chosen
c    input data set consisting of diatomic singlet-singlet transitions,
c    and calculate deviations [calc.-obs.]
c (ii)  Fit a data set made up of any combination of MW, IR or
c    electronic vibrational bands, and/or fluorescence series, involving
c    one or more electronic states and one or more isotopomers, to
c    parameters defining the observed levels of each state.
c
c=======================================================================
c** Dimensioning parameters intrinsic to the program are input through
c   'nparams.h'
c
c** Parameters characterizing the problem and governing the fits are
c   read on  channel-5  while the raw data are read on  channel-4 .
c   Principle output goes to  channel-6  while higher channel numbers
c   are used for secondary or more detailed/voluminous output.
c***********************************************************************
      INCLUDE 'arrsizes.h'

      INTEGER NSTATES,PROMC,CNDIFF,OSEL(8,NSTATEMX),DATN,
     1  JTRUNC(NSTATEMX),NDAT(0:NVIBMX,NISTPMX,NSTATEMX),
     2  I,J,ISTATE,IISTP,CHARGE,NOWIDTHS,EFSEL(NSTATEMX)

      COMMON /WIDTHBLK/NOWIDTHS

      REAL*8 UCUTOFF,ZMASE,RMAX(NSTATEMX),DECM(NSTATEMX),XX
      COMMON /RMAXBLK/RMAX
c                                                      
c** Parameters required within NLLSSRR.
c
      INTEGER NPTOT,IROUND,LPRINT,NPPS(NSTATEMX),                  
     1        VMIN(NSTATEMX),VMAX(NSTATEMX),IFTE
      COMMON /NPPSBLK/NPPS
                                                                       
      REAL*8 PV(NPARMX),PU(NPARMX),PS(NPARMX),CM(NPARMX,NPARMX),
     1  TSTPS,TSTPU,DSE,TEVAL(NPARMX)

      COMMON /NLSBLK/NSTATES,VMIN,VMAX,IFTE,TEVAL
c
c** MASSES variables.
c
      CHARACTER*2  CATOM

      INTEGER GELGS(2,NISTPMX),GNS(2,NISTPMX)

      REAL*8 ABUND(2,NISTPMX)
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
      LOGICAL outputPrint

      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  SIFDE(NSTATEMX),SIFRE(NSTATEMX),SIFBETA(0:NBETAMX,NSTATEMX),
     3  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      REAL*8 Rd,Vd,Vsr(NPNTMX,NSTATEMX),bMi,bMo,Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Vsr,Bsr,nPointSR

      REAL*8 Plr(NPNTMX,NSTATEMX),Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Plr,Blr,nPointLR

c** For Double-Exponential Long-Range Potential
      INTEGER numCn(NSTATEMX),nLR(NSTATEMX,numCnMax),tempN,step,Df
      REAL*8 dampBeta(NSTATEMX),Cn(NSTATEMX,numCnMax),Demo(NSTATEMX),
     &  DeOffset,ReOffset,Dn(NSTATEMX,numCnMax,NPNTMX),
     &  Vlr(NPNTMX,NSTATEMX),Dm,Dmp,Dmpp,A(NSTATEMX),B(NSTATEMX),dAdRe,
     &  dBdRe,dVdBtemp,ReBar

      COMMON /DELRBLK/numCn,nLR,dampBeta,Cn,ReBar,Df
      COMMON /OffSetBLK/Demo,DeOffset,ReOffset

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      COMMON /OUTPUT/outputPrint
*      COMMON /POPBLK/IOMEG,NDATPT,NCN,RMIN,RH,VLIM,EPS
c
c** Common block that will hold all of the information for the potential
c   and the partial derivative arrays.
c
c  BETA     holds the actual exponent parameters used for the potential
c
c  NDATPT   number of data points used for each potential
c
c  VPOT	holds the potential arrays for each state
c------------------------------------------------------------------------
c
c** Born-Oppenheimer Breakdown Parameters.
c
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NQA(NSTATEMX),NQB(NSTATEMX),
     1  IFUA(0:NABOMX,NSTATEMX),IFUB(0:NABOMX,NSTATEMX),
     2  IFQA(0:NRNAMX,NSTATEMX),IFQB(0:NRNAMX,NSTATEMX),
     3  SIFUA(0:NABOMX,NSTATEMX),SIFUB(0:NABOMX,NSTATEMX),
     4  SIFQA(0:NRNAMX,NSTATEMX),SIFQB(0:NRNAMX,NSTATEMX)
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
      INTEGER nLamdCoef(NSTATEMX),ifLamdCoef(0:NRNAMX,NSTATEMX),
     1                                  SifLamdCoef(0:NRNAMX,NSTATEMX)
      REAL*8 lamdCoef(0:NRNAMX,NSTATEMX),lambdaR(NPNTMX,NSTATEMX)
      COMMON /lambdaBLK/nLamdCoef,ifLamdCoef,lamdCoef,lambdaR
**************************
      REAL*8 RR(NPNTMX),RM2(NPNTMX),VV(NPNTMX),VLIMT
      INTEGER NCNN

      DATA ZMASE /5.48579903D-04/
c=======================================================================
      SLABL(-2)='WI'
      SLABL(-1)='BV'
      SLABL(0)='FS'
c=======================================================================
      DO 5 I=1,NPARMX
        PV(I) = 0.0d0
        PU(I) = 0.0d0
        PS(I) = 0.0d0
    5 CONTINUE
c=======================================================================
c** Start by reading parameters describing the overall nature of the
c   case and placing chosen restrictions on the data set to be used.
c
c  AN(1) & AN(2) are atomic numbers identifying the atoms forming the
c          molecule.
c
c  CHARGE (+/- integer) is the charge on the molecule (=0 for neutral).
c  If(CHARGE.ne.0) use Watson's(JMS 1980) charge-modified reduced mass.
c
c  NISTP   is the number of isotopomers to be simultaneously considered.
c
c  NSTATES  is the number of electronic states associated with the data
c        set to be analysed:  NSTATES = 1  for fits to IR/MW and/or
c        fluorescence data for a single electronic state, while
c        NSTATES > 1  for multi-state fits.
c        Upper states of fluorescence series NOT included in this count.
c-----------------------------------------------------------------------
      READ(5,*)  AN(1),AN(2),CHARGE,NISTP,NSTATES
c-----------------------------------------------------------------------
*      IF (NISTP.LT.1) THEN
*        WRITE(6,666)
*        NISTP = 1
*      ELSEIF (NISTP.GT.NISTPMX) THEN
*        WRITE(6,666)
*        NISTP = NISTPMX
*      END IF
*      IF (NSTATES.LT.1) THEN
*        WRITE(6,666)
*        NSTATES = 1
*      ELSEIF (NSTATES.GT.NSTATEMX) THEN
*        WRITE(6,666)
*        NSTATES = NSTATEMX
*      END IF
c-----------------------------------------------------------------------
c  UCUTOFF   Neglect any input data with uncertainties > UCUTOFF (cm-1)
c
c  NOWIDTHS >  0  causes the program to ignore any tunneling widths in
c                 the data set.
c           <= 0  causes the program to fit to tunneling widths
c           <  0  use simple version of dWdP, ignoring the partial 
c                 derivative of t_vib which involves k = 1 phase integral
c  IROUND  specifies the level of rounding inside NLLSSRR if:
c          > 0 : requires that Sequential Rounding & Refitting be
c                performed, with each parameter being rounded at the
c                IROUND'th sig. digit of its local uncertainty.
c          <=0 : simply stops after full convergence (without rounding).
c
c  LPRINT  specifies the level of printing inside NLLSSRR if:
c            = 0 : no print except for failed convergence.
c           <  0 : only converged, unrounded parameters, PU & PS's
c           >= 1 : print converged parameters, PU & PS's
c           >= 2 : also print parameter change each rounding step
c           >= 3 : also indicate nature of convergence
c           >= 4 : also print convergence tests on each cycle
c           >= 5 : also parameters changes & uncertainties, each cycle
c
c  PROMC   <> 0  causes a statistical summary of predicted [obs.-calc.]
c                results to be printed (on channel-6).
c           > 0  causes a full listing of [obs.-calc.] values to be
c                written to channel-8.
c
c  CNDIFF  flag to denote use of differences for adiabatic BO breakdown
c          corrections
c          = 0    fit using differences
c          = 1    fit using clamped nuclei
c-----------------------------------------------------------------------
      READ(5,*)  UCUTOFF,NOWIDTHS,IROUND,LPRINT,PROMC,CNDIFF
c-----------------------------------------------------------------------
      I= 999
      WRITE(6,600)
      WRITE(6,601) NISTP
      DO 10 IISTP= 1,NISTP
c
c** Read the mass numbers of the atoms in each of the isotopomers
c
c  MN(i,IISTP)  is the mass number for atom with atomic number AN(i)
c       [NOTE: be sure order of MN values consistent with that of AN's].
c       Choosing it .ne. value for some known isotope if that species
c       causes the average atomic mass to be used.
c-----------------------------------------------------------------------
        READ(5,*) MN(1,IISTP),MN(2,IISTP)
c-----------------------------------------------------------------------
        I= MIN(I,MN(1,IISTP),MN(2,IISTP))
        CALL MASSES(AN(1),MN(1,IISTP),CATOM,GELGS(1,IISTP),
     1              GNS(1,IISTP),ZMASS(1,IISTP),ABUND(1,IISTP))
        IF (IISTP.EQ.1) NAME(1) = CATOM
        CALL MASSES(AN(2),MN(2,IISTP),CATOM,GELGS(2,IISTP),
     1              GNS(2,IISTP),ZMASS(2,IISTP),ABUND(2,IISTP))
        IF (IISTP.EQ.1) NAME(2) = CATOM
        ZMASS(3,IISTP) = (ZMASS(1,IISTP)*ZMASS(2,IISTP))/
     1                     (ZMASS(1,IISTP)+ZMASS(2,IISTP)-CHARGE*ZMASE)
        WRITE(6,602) NAME(1),MN(1,IISTP),NAME(2),
     1       MN(2,IISTP),(ZMASS(J,IISTP),J=1,3)
        XX = ZMASS(3,1)/ZMASS(3,IISTP)
        DO J=0,8
          RMUP(J,IISTP)= XX**J
        END DO
        RSQMU(IISTP) = DSQRT(XX)
        XX= 1.0d0
        DO J=0,NDUNMX
          RSQMUP(J,IISTP) = XX
          XX = XX*RSQMU(IISTP)
        END DO
c
c** Calculate mass scaling terms for the adiabatic (ZMUA, ZMUB) and
c   rotational non-adiabaric (ZMQA, ZMQB) Born-Oppenheimer breakdown
c   functions.
c** For Watson compatiblity.
c
        IF (CNDIFF.EQ.1) THEN
          ZMUA(IISTP) = 1.0D0 / ZMASS(1,IISTP)
          ZMUB(IISTP) = 1.0D0 / ZMASS(2,IISTP)
          ZMQA(IISTP) = 1.0d0 / ZMASS(1,IISTP)
          ZMQB(IISTP) = 1.0d0 / ZMASS(2,IISTP)
c
c** Using mass differences for the adiabatic corrections (ZMUA, ZMUB):
c
c          ZMUA = (ZMASSA(IISTP) - ZMASSA(1)) / ZMASSA(IISTP)
c          ZMUB = (ZMASSB(IISTP) - ZMASSB(1)) / ZMASSB(IISTP)
c
c   and mass ratios for the rotational corrections (ZMQA, ZMQB):
c
c                 ZMQA = ZMASSA(1) / ZMASSA(IISTP)
c                 ZMQB = ZMASSB(1) / ZMASSB(IISTP)
c
        ELSE
          ZMUA(IISTP) = 1.0d0 - ZMASS(1,1)/ZMASS(1,IISTP)
          ZMUB(IISTP) = 1.0d0 - ZMASS(2,1)/ZMASS(2,IISTP)
          ZMQA(IISTP) = ZMASS(1,1) / ZMASS(1,IISTP)
          ZMQB(IISTP) = ZMASS(2,1) / ZMASS(2,IISTP)
        END IF
c
c** For homonuclear diatomics, set the first mass scaling term for each
c   set of correction terms to be the sum of the two original mass
c   scaling factors, and set the second mass term to zero. 
c
        IF (AN(1).EQ.AN(2)) THEN
          ZMUA(IISTP) = ZMUA(IISTP) + ZMUB(IISTP)
          ZMUB(IISTP) = 0.0d0
          ZMQA(IISTP) = ZMQA(IISTP) + ZMQB(IISTP)
          ZMQB(IISTP) = 0.0d0 
        END IF
c
c** When fitting to a single isotopic species using differences, only
c   one set of correction terms (QAR) can be determined. Therefore, set
c   the mass scaling factor ZMQA to be the sum of the originial
c   rotational non-adiabatic mass scaling terms. Also set the remaining
c   mass terms to zero.
c
        IF ((CNDIFF.EQ.0).AND.(NISTP.EQ.1)) THEN
          ZMUA(IISTP) = 0.0d0
          ZMUB(IISTP) = 0.0d0
          ZMQA(IISTP) = ZMQA(IISTP) + ZMQB(IISTP)
          ZMQB(IISTP) = 0.0d0 
        END IF
   10 CONTINUE
      IF (I.EQ.0) WRITE(6,603)
      WRITE(6,*)
      IF (AN(1).EQ.AN(2)) WRITE(6,604)
      IF (CNDIFF.EQ.0) THEN
        WRITE(6,605)
        IF (NISTP.EQ.1) WRITE(6,606)
      END IF
  600 FORMAT(/,1X,39('--'))
  601 FORMAT(2X,'Input data for',I3,'  isotopomer(s)'/2X,16('**')/2X,
     1  '    Isotopomer       Mass of atom-1   Mass of atom-2    Reduced
     2 mass'/ 2X,'----------------- ',3('   --------------'))
  602 FORMAT(2X,A2,'(',I3,') - ',A2,'(',I3,')',3(3X,F14.9))
  603 FORMAT('  Note that   (Mass Number) = 0   causes the average atomi
     1c mass to be used.')
  604 FORMAT('  For homonuclear diatomics, BO correction terms are the s
     1ame for both atoms.'/'  Only the first set of correction terms wil
     2l be used, UA(R) and QA(R),'/'  with a mass scaling factor equal t
     3o the sum of the two individual terms.')
  605 FORMAT('  Using differences for adiabatic BOB mass scaling terms.'
     1/'  Using ratios for rotational BOB mass scaling terms.')
  606 FORMAT('  For fitting to a single isotopic species'/'  only fit to
     1 one set of rotational BOB terms: QA(R)')
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now read in limits that will be used on the data set.
c
c  V(MIN/MAX)(s) Neglect data for electronic state vibrational levels
c                outside the range  VMIN  to  VMAX.
c
c  JTRUNC(s)     data with J > JTRUNC are not included in the fit.
c
c  EFSEL(s)  allows a user to consider data for:
c          * ONLY the e-parity levels of this state, if EFSEL > 0
c          * ONLY the f-parity levels of this state, if EFSEL < 0
c          * BOTH e- and f-parity levels of thsi state, if EFSEL = 0
c-----------------------------------------------------------------------
      DO 20 ISTATE=1,NSTATES 
        READ(5,*) VMIN(ISTATE),VMAX(ISTATE),JTRUNC(ISTATE),EFSEL(ISTATE)
c-----------------------------------------------------------------------
c** Now to read in potential parameters that are not model specific.
c
c  SLABL(s)  is a 2-character alphameric label enclosed in single quotes
c            to identify the electronic state; e.g., 'X0', 'A1', ... etc.
c
c  IOMEG(s)      is the electronic angular momentum of each state.
c
c  NCN(s)   specifies the limiting inverse-power behaviour of the 
c           potential tail V(R) -> VLIM - CNN/R^NCN
c
c  VLIM(s)       is the absolute energy of the potential asymptote
c                (cm-1).
c-----------------------------------------------------------------------
        READ(5,*) SLABL(ISTATE),IOMEG(ISTATE),NCN(ISTATE),VLIM(ISTATE)
c-----------------------------------------------------------------------
c** Read in variables for selecting output. 
c
c  OSEL(n,s) selects the type of output for each generated potential
c       OSEL < 0: Export to file every OSELth point of initial function
c       OSEL > 0: Export to file every OSELth point of final function
c       OSEL(1) : Potential function (V) exported to channel 10.
c       OSEL(2) : Beta function (BMORSE) exported to channel 11.
c       OSEL(3) : Adiabatic correction function (UAR) exported to
c                 channel 12.
c       OSEL(4) : Adiabatic correction function (UBR) exported to
c                 channel 13.
c       OSEL(5) : Rotational non-adiabatic correction function (QAR)
c                 exported to channel 14.
c       OSEL(6) : Rotational non-adiabatic correction function (QBR)
c                 exported to channel 15.
c       OSEL(7) : Lambda doubling correction function (lambdaR)
c                 exported to channel 16.
c       OSEL(8) : Write out partial derivatives to channels 31 to 35.
c-----------------------------------------------------------------------
        READ(5,*) (OSEL(I,ISTATE),I=1,8)
c-----------------------------------------------------------------------
c** Now to read in the range over which the potentials will be defined.
c
c  R(MIN/MAX)  define the range over which the potentials are defined.
c
c  RH(s)       specifies the meshvalue of the potential (distance
c              between two points that are defined on the potential).
c-----------------------------------------------------------------------
        READ(5,*) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
c-----------------------------------------------------------------------
        NDATPT(ISTATE) = (RMAX(ISTATE)-RMIN(ISTATE))/RH(ISTATE)+1.0001d0
        NDATPT(ISTATE) = MIN(NPNTMX,NDATPT(ISTATE))
        RMAX(ISTATE) = RMIN(ISTATE) + RH(ISTATE)*DBLE(NDATPT(ISTATE)-1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL READPOT(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   20 CONTINUE
      WRITE(6,609) NSTATES
      DO ISTATE = 1,NSTATES
        WRITE(6,610) SLABL(ISTATE),RMIN(ISTATE),RMAX(ISTATE),RH(ISTATE)
      END DO
  609 FORMAT(/'  Limits and increment of integration for',I3,'  state(s)
     1'/2X,26('**')/3X,'Electronic State     Minimum Limit    Maximum Li
     2mit      Mesh size'/'  ----------------- ',3('   --------------'))
  610 FORMAT(8X,A3,8X,3(3X,F14.9))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(6,620) NSTATES
      DO ISTATE = 1,NSTATES
        WRITE(6,621) SLABL(ISTATE), DE(ISTATE),
     1               DE(1)-VLIM(1)+VLIM(ISTATE)-DE(ISTATE), VLIM(ISTATE)
      END DO
*  620 FORMAT(/'  Atomic energies and potential asymptotes for',I3,'  sta
*     1te(s)'/2X,57('*')/3X,'Electronic State     Atom-A Limit     Atom-B
*     2 Limit    Pot. Asymptote'/'  ----------------- ',3('   -----------
*     3---'))
  620 FORMAT(/'  Potential energies and asymptotes for ',I3,'  state(s)'
     1      /2X,52('*')/3X,'Electronic State     Dissociation       Te V
     2alue      Pot. Asymptote'/'  ----------------- ',3('   -----------
     3---'))
  621 FORMAT(8X,A3,8X,3(3X,F14.6))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to input experimental data in specified
c   band-by-band, and do bookkeeping to characterize amounts of data or
c   each type.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL READATA(NSTATES,UCUTOFF,JTRUNC,EFSEL,VMIN,VMAX,NDAT,
     &                                                  NOWIDTHS,PROMC)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generating internal NLLSSRR array (PV) from external variables and
c   calculating the number of "free" parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      NPTOT = NPARMX
      CALL NMCHNG(NSTATES,NPTOT,PV,0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine that will generate the potential function and its 
c   partial derivatives for each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO ISTATE=1,NSTATES
        CALL VGEN(ISTATE,-1.0d0,Vd,1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine PREPOT that will generate the potential function 
c   for each state if forward calculation is required.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (PSEL(NSTATES) .EQ. 0) THEN
          DO I = 1, NDATPT(ISTATE)
            RR(I) = RMIN(ISTATE) + DBLE(I) * RH(ISTATE) - RH(ISTATE)
            RM2(I) = 1.0d0 / (RR(I) ** 2)
          END DO
          CALL PREPOT(1,AN(1),AN(2),MN(1,1),MN(2,1),NDATPT(ISTATE),RR,
     &                RM2,VLIMT,VV,NCNN)
          DO I = 1, NDATPT(ISTATE)
            VPOT(I,ISTATE) = VV(I) - DE(ISTATE) + VLIM(ISTATE)
          END DO
          J = 0.05 / RH(ISTATE)
          DO I = 1,NDATPT(ISTATE),J
            WRITE(17,900) RR(I),VPOT(I,ISTATE)
          END DO
        END IF
      END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set the energy convergence criterion to be 1/100th of the smallest
c   experimental uncertainty. 
c
      DO ISTATE=1,NSTATES
        EPS(ISTATE) = UCUTOFF/100.0d0
c
c** Initialize the dissociation energy correlation matrix.
c
        DECM(ISTATE) = 0.0d0
      END DO
      NPTOT = NPTOT + NFSTOT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If fitting to fluorescence data, set the fluorescence band
c   parameters to zero:
c
c                    PV(NPTOT-NFSTOT+1 < I < NPTOT).
c
      IFTE = 0
      DO I = 1,NFSTOT
        TEVAL(I) = 0.0d0
      END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to write out a summary of the initial potential parameters for
c   each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL WRIDATA(1,NSTATES,SLABL,NISTP,NAME,DECM,PU,PS)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine NLLSSRR to calculate converged parameters from trial
c   values and spectroscopic data.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (IROUND.NE.0) WRITE(6,685) ABS(IROUND)
      IF (IROUND.GT.0) WRITE(6,686)
      IF (IROUND.LT.0) WRITE(6,687)
  685 FORMAT(/' Apply "Sequential Rounding & Refitting" at digit-',
     1  i1,' of the (local) parameter')
  686 FORMAT(4x,'uncertainty, selecting remaining parameter with largest
     1 relative uncertainty')
  687 FORMAT(4x,'uncertainty, proceeding sequentially from the LAST para
     1meter to the FIRST.')
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NFSTOT.GT.0) THEN
c** If have fluorescence series ... first fix ALL potential parameters
c  and fit to determine the series origins, and only THEN free potential
c  parameters too.  First, save read-in values of  'IF(fix)' parameters
          DO  ISTATE= 1, NSTATES
              SIFDE(ISTATE)= IFDE(ISTATE)
              SIFRE(ISTATE)= IFRE(ISTATE)
              IFDE(ISTATE)= 1
              IFRE(ISTATE)= 1
            DO  I= 0,NBETA(ISTATE)
                  SIFBETA(I,ISTATE)= IFBETA(I,ISTATE)
                  IFBETA(I,ISTATE)= 1
            ENDDO
            IF(NUA(ISTATE).GT.0) THEN
              DO  I= 0,NUA(ISTATE)
                      SIFUA(I,ISTATE)= IFUA(I,ISTATE)
                      IFUA(I,ISTATE)= 1
              ENDDO
            ENDIF
            IF(NUB(ISTATE).GT.0) THEN
              DO  I= 0,NUB(ISTATE)
                      SIFUB(I,ISTATE)= IFUB(I,ISTATE)
                      IFUB(I,ISTATE)= 1
              ENDDO
            ENDIF
            IF(NQA(ISTATE).GT.0) THEN
              DO  I= 0,NQA(ISTATE)
                      SIFQA(I,ISTATE)= IFQA(I,ISTATE)
                      IFQA(I,ISTATE)= 1
              ENDDO
            ENDIF
            IF(NQB(ISTATE).GT.0) THEN
              DO  I= 0,NQB(ISTATE)
                      SIFQB(I,ISTATE)= IFQB(I,ISTATE)
                      IFQB(I,ISTATE)= 1
              ENDDO
            ENDIF
            IF(nLamdCoef(ISTATE).GE.0) THEN
              DO  I= 0,nLamdCoef(ISTATE)
                      SifLamdCoef(I,ISTATE)= ifLAMdCoef(I,ISTATE)
                      ifLamdCoef(I,ISTATE)= 1
              ENDDO
            ENDIF
          ENDDO
          CALL NMCHNG(NSTATES,NPTOT,PV,0)
          NPTOT = NPTOT + NFSTOT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generating internal NLLSSRR array (PV) from external variables and
c   calculating the number of "free" parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
*c?????????????????????/
*
*c ... now call NLLSSRR to get Fluorescence series origins ....
        CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,IROUND,NFSTOT,LPRINT,
     1                    FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
*          CALL NLLSSRR(COUNTOT,NFSTOT,NPARMX,IROUND,NFSTOT,LPRINT,
*     1                   FREQ,UFREQ,DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c  parameters too.  First, save read-in values of  'IF(fix)' parameters
c  Now ... reset "IF(fix)" parameterd to read-in values ... & proceed ..
          DO  ISTATE= 1, NSTATES
              IFDE(ISTATE)= SIFDE(ISTATE)
              IFRE(ISTATE)= SIFRE(ISTATE)
            DO  I= 0,NBETA(ISTATE)
                  IFBETA(I,ISTATE)= SIFBETA(I,ISTATE)
            ENDDO
            IF(NUA(ISTATE).GT.0) THEN
              DO  I= 0,NUA(ISTATE)
                      IFUA(I,ISTATE)= SIFUA(I,ISTATE)
              ENDDO
            ENDIF
            IF(NUB(ISTATE).GT.0) THEN
              DO  I= 0,NUB(ISTATE)
                      IFUB(I,ISTATE)= SIFUB(I,ISTATE)
              ENDDO
            ENDIF
            IF(NQA(ISTATE).GT.0) THEN
              DO  I= 0,NQA(ISTATE)
                      IFQA(I,ISTATE)= SIFQA(I,ISTATE)
              ENDDO
            ENDIF
            IF(NQB(ISTATE).GT.0) THEN
              DO  I= 0,NQB(ISTATE)
                      IFQB(I,ISTATE)= SIFQB(I,ISTATE)
              ENDDO
            ENDIF
            IF(nLamdCoef(ISTATE).GE.0) THEN
              DO  I= 0,nLamdCoef(ISTATE)
                      ifLamdCoef(I,ISTATE)= SifLAMdCoef(I,ISTATE)
              ENDDO
            ENDIF
          ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generating internal NLLSSRR array (PV) from external variables and
c   calculating the number of "free" parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c?????????????????????/
        DO I = 1,NFSTOT
          TEVAL(I) = PV(I)
        END DO
        CALL NMCHNG(NSTATES,NPTOT,PV,0)
        NPTOT = NPTOT + NFSTOT
        DO I = 1,NFSTOT
          PV(NPTOT-NFSTOT+I) = TEVAL(I)
        END DO
      ENDIF
c
      WRITE(22,688)
  688 FORMAT('   v       ','E',10x,'Bv',11x,'-Dv',
     1  13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
      CALL NLLSSRR(COUNTOT,NPTOT,NPARMX,IROUND,NFSTOT,LPRINT,FREQ,UFREQ,
     1             DFREQ,PV,PU,PS,CM,TSTPS,TSTPU,DSE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Writing out the general information of the fit.
c-----------------------------------------------------------------------
      WRITE(6,690)
      WRITE(6,691) NPTOT,COUNTOT,DSE
c-----------------------------------------------------------------------
c** Writing out the fluorescence band results.
c-----------------------------------------------------------------------
      IF (NFSTOT.GT.0) THEN
        WRITE(6,690)
        WRITE(6,692)
        WRITE(6,693)
        J = NPTOT - NFSTOT
        DO I = 1,NFSTOT
          WRITE(6,694) 'T',VP(FSBAND(I)),VPP(FSBAND(I)),
     &   EFP(IFIRST(FSBAND(I))),ISTP(FSBAND(I)),TEVAL(I),PU(J+I),PS(J+I)
        END DO
      ENDIF
  690 FORMAT(/,1X,39('=='))
  691 FORMAT(/' For fit of',I5,'  parameters to',I5,'  transitions, DSE=
     1',G15.8)
  692 FORMAT(/" The following Fluorescence Series Term Values were deter
     &mined.")
  693 FORMAT(/" T(v', J', p';ISTP)      Final Value:    Uncertainty:  Se
     &nsitivity:")
  694 FORMAT(2X,A1,'(',I3,',',I3,',',SP,I3,SS,';',I2,')',1X,1PD19.10,3X,
     &1PD8.1,6X,1PD8.1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculation of the uncertainties for Te for each potential require 
c   elements from the correlation matrix.
c
      I = 1
      DO ISTATE = 1,NSTATES
        IF ((IFDE(1).EQ.0).AND.(IFDE(ISTATE).EQ.0)) THEN
          DECM(ISTATE) = CM(1,I)
        ELSE
          DECM(ISTATE) = 0.0d0
        END IF
        I = I + NPPS(ISTATE)
      END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine NMCHNG to convert internal NLLSSRR parameter array
c   back into external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL NMCHNG(NSTATES,NPTOT,PV,1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine VGEN to generate the potential function from the
c   final calculated converged parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      outputPrint = .TRUE.
      IF (PSEL(NSTATES) .NE. 0) THEN
        DO ISTATE = 1,NSTATES
          nPointSR(ISTATE) = RMIN(ISTATE)/RH(ISTATE) 
          IF (OSEL(1,ISTATE) .NE. 0) THEN
            IF (RMAX(ISTATE) .GT. 80.0) THEN
              nPointLR(ISTATE) = 0
            ELSE
              nPointLR(ISTATE)=(80.0-RMAX(ISTATE))
     &                                   /(RH(ISTATE)*OSEL(1,ISTATE))
            END IF
          END IF
          DO J = 1,numCn(ISTATE)
            tempN = nLR(ISTATE,J)
            WRITE (7,877) nLR(ISTATE,J)
            step = 0.05 / RH(ISTATE)
            DO I = 1,nPointSR(ISTATE),step
              Rd = 0.1d0 + RH(ISTATE)*DBLE(I-1)
              CALL dampF(Rd,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
              WRITE (7,879) Rd, Dm
            END DO
          END DO
          CALL VGEN(ISTATE,-1.0d0,Vd,1)
          IF (OSEL(1,ISTATE) .NE. 0) THEN
            DO I = 1,nPointSR(ISTATE),OSEL(1,ISTATE)
              Rd = RH(1)*DBLE(I-1) 
              CALL vValue(ISTATE,Rd,Vd,bMi,bMo,0,.FALSE.)
              Vsr(I,ISTATE) = Vd
              Bsr(I,ISTATE) = bMi
            END DO
          END IF 
          IF (OSEL(1,ISTATE) .NE. 0) THEN
            DO I = 1,nPointLR(ISTATE)
              Rd = RMAX(ISTATE) + RH(ISTATE)*DBLE(I*OSEL(1,ISTATE))
              CALL vValue(ISTATE,Rd,Vd,bMi,bMo,0,.FALSE.)
              Plr(I,ISTATE) = Vd
              Blr(I,ISTATE) = bMo
            END DO
          END IF
        END DO
      END IF
  877 FORMAT (' n =',I3)
  879 FORMAT (4x,2F10.5)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out a summary of the converged and fixed
c   values to standard output and a list of the parameters to channel-9.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL WRIDATA(2,NSTATES,SLABL,NISTP,NAME,DECM,PU,PS)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If chosen, output file(s) will be created for the export of the
c   generated functions: V, BMORSE, UAR/UBR, or QAR/QBR and their
c   respective uncertainties.
c** First rewind output files so as to write over debugging information.
c
      DO ISTATE = 1,NSTATES
        IF (OSEL(1,ISTATE).NE.0) REWIND(10)
        IF (OSEL(2,ISTATE).NE.0) REWIND(11)
        IF (OSEL(3,ISTATE).NE.0) REWIND(12)
        IF (OSEL(4,ISTATE).NE.0) REWIND(13)
        IF (OSEL(5,ISTATE).NE.0) REWIND(14)
        IF (OSEL(6,ISTATE).NE.0) REWIND(15)
        IF (OSEL(7,ISTATE).NE.0) REWIND(16)
      END DO
      IF (PSEL(NSTATES) .NE. 0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out the generated functions and their
c   respective uncertainties.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL FUNUNC(NSTATES,OSEL,NPPS,PU,CM,XX)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculate the first five partial derivatives of the potential and
c   print them out to channel-31 to 35.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        CALL DVCALC(NSTATES,OSEL)
      ELSE 
        J = 0.05 / RH(NSTATES)
        DO I = 1,NDATPT(NSTATES),J
          WRITE(17,900) RR(I),VPOT(I,NSTATES)
        END DO
  900   FORMAT(5X,G18.8,5X,G18.8)
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to print out summary of dimensionless standard
c   errors on a band-by-band basis, and (if desired) print [obs.-calc.]
c   values to channel-8.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (PROMC.NE.0) CALL DIFFSTATS(NSTATES,PROMC)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
