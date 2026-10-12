***********************************************************************
      SUBROUTINE DYIDPJ(DATN,NPTOT,YC,PV,PD)
c***********************************************************************
c** This program assumes that the upper state IEP is at a higher energy
c   than the lower state IEPP.
c** This subroutine will fit to the predefined potentials given in VGEN.
c   Will calculate the PD(j) for each datum (DATN) with respect to each
c   of the parameters (PV) varied in the fit. (every parameter for each 
c   potential state)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 1997  by  Jenning Seto       (ver of 12/10/1997) +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    DATN     is the number of the current observable being considered. 
c    NPTOT    is the total number of parameters being varied.
c    PV(i)    is the array of parameters being varied.
c
c** On exit:
c    YC       is the calculated value of the DATNth observable.
c    PV(i)    is the array of parameters being varied.
c    PD(i)    is the partial derivative array (de/dp).
c
      INCLUDE 'arrsizes.h'
      INTEGER DATN,NPTOT,NBAND,IISTP,ISTATE,VIMX,I,J
c
      REAL*8 YC,PV(NPTOT),PD(NPTOT),
     1  EHI,ELO,UPPER(NPARMX),LOWER(NPARMX)
     2  ,EO,DEDPK(NPARMX),Vd
c ** Width parameters
      REAL*8 width,dWdPk(NPARMX)
c
c** Parameters required within NLLSSRR.
c
      INTEGER NSTATES,NPPS(NSTATEMX),VMIN(NSTATEMX),VMAX(NSTATEMX),IFTE
      COMMON /NPPSBLK/NPPS
c
      REAL*8 TEVAL(NPARMX)
      COMMON /NLSBLK/NSTATES,VMIN,VMAX,IFTE,TEVAL
c
c** Vibrational energies and rotational constants.
c
      REAL*8 ZK(0:NVIBMX,0:MORDRMX,NISTPMX,NSTATEMX)
      COMMON /ZKBLK/ZK
c
c** Type statements & common block for data
c
      REAL*8 FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)
c
      INTEGER COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)
c
      CHARACTER*2 NAME(2)
      CHARACTER*2 SLABL(-3:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL

      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),RM2(NPNTMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      INTEGER IPASS

      SAVE IPASS,RM2

      DATA IPASS /0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Need to update partial derivatives and trial energy values for
c   every new itteration (DATN=1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(DATN.EQ.1) THEN
         IPASS = IPASS + 1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to use internal NLLSSRR parameter array PV to update
c   external (physical) variable set.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         CALL NMCHNG(NSTATES,NPTOT,PV,1)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If fitting to fluorescence data, convert fluorescence band data back
c   to external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ((NFSTOT.GT.0).AND.(IFTE.EQ.0)) THEN
            J = NPTOT - NFSTOT
            DO I = 1,NFSTOT
               TEVAL(I) = PV(J+I)
               END DO
            END IF
         DO ISTATE=1,NSTATES
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine to generate the potential function and its partial 
c   derivatives with respect to each fitted parameter for each state.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF (PSEL(ISTATE) .NE. 0) CALL VGEN(ISTATE,-1.0d0,Vd,1)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to generate the molecular constants for each state and
c   isotopomer for generating trial eigenvalues in fit calculations.
c** To take into account the different possible vibrational ranges for
c   different isotopomers, only read in values that can be found up to
c   the reported VMAX from ALF. (And cross your fingers.)
            DO IISTP=1,NISTP
               VIMX = VMAX(ISTATE) * DSQRT(ZMASS(3,IISTP)/ZMASS(3,1))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine INITDD to calculate molecular constants for trial
c   energy calculations as well as the partial derivatives.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
               CALL INITDD(ISTATE,IISTP,VIMX,ZMASS(3,IISTP),RM2,ZK)
               END DO
            END DO
         END IF
c ++++++++++++++++++++++++++++++++++++++++++end of datum-1 set-up ++++++
c** Initialize variables for current datum
      EHI = 0.0D0
      ELO = 0.0D0
      DO I = 1,NPTOT
         PD(I) = 0.0d0
         END DO
      DO I = 1,NPARMX
         UPPER(I) = 0.0d0
         LOWER(I) = 0.0d0
         END DO
      NBAND = IB(DATN) 
      IISTP = ISTP(NBAND)
c
c** Calculate the partial derivatives and BETA values for fitted BETA
c   data.
c
*      IF (IEP(NBAND).LT.0) THEN
*        ISTATE = -IEP(NBAND)
*        RTEMP = (RDPNT(DATN) - RE(ISTATE))/(RDPNT(DATN) + RE(ISTATE))
*        DO I = 1,NBETA(ISTATE)
*          DVDBi = RTEMP**I
*          DVDR = BETA(I,ISTATE)*RTEMP**(I-1)
*        END DO
*        DVDR = - DVDR * RDPNT(DATN) / (RDPNT(DATN) + RE(ISTATE))**2
c
c** Now to determine the partials with respect to the upper and lower
c   state for fluorescence data.
c
      IF (IEP(NBAND).EQ.0) THEN
        J = NPTOT - NFSTOT
        DO I = 1,NFSTOT
          IF (FSBAND(I).EQ.NBAND) THEN
c** PV(J+I) is the energy of the I'th fluorescence band.
            EHI = TEVAL(I)
            IF (IFTE.EQ.0) UPPER(J + I) = 1.0d0
          ENDIF
        END DO
        CALL DEDP(DATN,IEPP(NBAND),IISTP,NPPS,ZMASS(3,IISTP),
     1                      JP(DATN),JPP(DATN),ELO,width,ZK,RM2,LOWER)
c** For rotational constants
      ELSE IF (IEP(NBAND).EQ.-1) THEN
        WRITE(6,*) 'Rotational constants not implemented in DsPotFit'
c** For widths of tunneling predissociation quasibound levels
      ELSE IF (IEP(NBAND).EQ.-2) THEN
        CALL DEDP(DATN,IEPP(NBAND),IISTP,NPPS,ZMASS(3,IISTP),
     1                       JP(DATN),JPP(DATN),EO,width,ZK,RM2,DEDPK)
        IF(PSEL(IEPP(NBAND)).GT.0) 
c ... for forward calculation, use widths from SCHRQ in DEDP
     1      CALL DWDP(DATN,IEPP(NBAND),NPPS,FREQ(DATN),ZMASS(3,IISTP),
     2                        JP(DATN),JPP(DATN),EO,width,DEDPK,dWdPk)
      ELSE IF (IEP(NBAND).EQ.-3) THEN
c** For fitting to scattering length / v_D, do derivatives by
c   symmetric differences between two potentials
c
        CALL DVDDP(DATN,IEPP(NBAND),IISTP,NPPS,ZMASS(3,IISTP),
     1   width,dWdPk)
c 
c** Now to determine the partials with respect to the upper and lower
c   state for microwave, infrared, and visible data.
c
      ELSE
        CALL DEDP(DATN,IEP(NBAND),IISTP,NPPS,ZMASS(3,IISTP),
     1                      VP(NBAND),JP(DATN),EHI,width,ZK,RM2,UPPER)
        CALL DEDP(DATN,IEPP(NBAND),IISTP,NPPS,ZMASS(3,IISTP),
     1                    VPP(NBAND),JPP(DATN),ELO,width,ZK,RM2,LOWER)
      ENDIF
c
c** Collect only those de/dp for parameters being changed.
c
      DO I = 1,NPTOT
        PD(I) = UPPER(I) - LOWER(I)
      END DO
*** For first pass of fit to fluorescence bands, first "fix" potential
*** parameters and determine fluorescence band energies.
*      IF ((NFSTOT.GT.0).AND.(IPASS.EQ.1)) THEN
*        DO I = 1,NPTOT-NFSTOT+1
*          PD(I) = PD(I)*1.0d99
*        END DO
*      ENDIF
c
c** Find the calculated value for the DATNth observable from energy
c   levels.
c
      YC = EHI - ELO
      IF (IEP(NBAND).LE.-2) THEN
         DO I = 1,NPTOT
            PD(I) = dWdPk(i)
            END DO
         YC = width
         END IF 
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE INITDD(ISTATE,IISTP,VIMX,ZMASS,RM2,ZK)
c***********************************************************************
c** This subroutine updates the trial vibrational energies & rotational
c   constants on each iteration
c
c** On entry:   ISTATE  is the states being considered.
c               IISTP   is the isotopomer being considered.
c               VIMX    is the upper vibrational bound for each state.
c               ZMASS   is the molecular reduced mass for this isotopomer.
c** On exit:   RM2      is the array of values of 1/R**2
c              ZK       are the molecular constants for this system.
c
      INCLUDE 'arrsizes.h'
c
      INTEGER ISTATE, IISTP, VIMX, AFLAG, INNR(0:NVIBMX),
     1  VIN, NBEG, NEND, WARN, IWR, LPRWF, I, J
c
      REAL*8 ZMASS,ZK(0:NVIBMX,0:MORDRMX,NISTPMX,NSTATEMX),
     1  GV(0:NVIBMX),RCNST(MORDRMX),
     2  RM2(NPNTMX),FWHM,PMAX,BFCT,BvWN,V(NPNTMX),SWF(NPNTMX)
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)
c
      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR
c
c** Born-Oppenheimer Breakdown Functions.
c     
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1  QAR(NPNTMX,NSTATEMX),QBR(NPNTMX,NSTATEMX),
     2  ZMUA(NISTPMX),ZMUB(NISTPMX),ZMQA(NISTPMX),ZMQB(NISTPMX)
c
      COMMON /BFBLK/UAR,UBR,QAR,QBR,ZMUA,ZMUB,ZMQA,ZMQB
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      VIN = VIMX
      AFLAG = 0
      WARN = 0
*      IWR = 2
      IWR = 0
      LPRWF = 0
c
c** Calculate potential scaling factor BFCT for ALF.
c
      BFCT = (ZMASS/16.85762908D0)*RH(ISTATE)*RH(ISTATE)
c
c** Now generate the BFCT scaled, adiabatically corrected, and
c   rotationally perturbed potentials for each isotope for use in
c   SCHRQ.
c
      DO I=1,NDATPT(ISTATE)
         V(I) = BFCT * (VPOT(I,ISTATE) + ZMUA(IISTP) * UAR(I,ISTATE)
     1                                  + ZMUB(IISTP) * UBR(I,ISTATE))
         END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine ALF that will calculate the vibrational levels and
c   rotational constants for each state and isotopomer.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL ALF(NDATPT(ISTATE),RMIN(ISTATE),RH(ISTATE),V,SWF,
     1 VLIM(ISTATE),VIN,AFLAG,ZMASS,EPS(ISTATE),NCN(ISTATE),GV,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a serious error occured during within ALF, then print out a 
c   warning, record the constants that we have and hope the program
c   doesn't call on the constants for which ALF could not calculate.
c
      IF (AFLAG.LT.0) THEN
          WRITE(6,600) ISTATE,IISTP
          IF (AFLAG.EQ.-1) WRITE(6,601) VIMX, AFLAG
          IF (AFLAG.EQ.-2) WRITE(6,602)
          IF (AFLAG.EQ.-3) WRITE(6,603)
          IF (AFLAG.EQ.-4) WRITE(6,604)
          IF (AFLAG.EQ.-5) WRITE(6,605)
          IF (AFLAG.EQ.-6) WRITE(6,606)
          IF (AFLAG.EQ.-8) THEN
              WRITE(6,608)
        ELSE
          WRITE(6,607) VIN
        END IF
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now to calculate the rotational constants for each vibrational level
c   of each state and isotopomer.
c
c** Generate the 1/R**2 and V(R) arrays for CDJOEL.
c
      BvWN = 16.85762908D0/ZMASS
      DO I = 1,NDATPT(ISTATE)
          RM2(I) = 1/(RMIN(ISTATE)+RH(ISTATE)*DBLE(I-1))**2
          END DO
      DO I = 0,VIN
          AFLAG = 0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to calculate the wavefunction required by
c   CDJOEL to calculate the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(I,AFLAG,GV(I),FWHM,PMAX,VLIM(ISTATE),V,SWF,BFCT,
     1        EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2        NBEG,NEND,INNR(I),IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine CDJOEL to determine the rotational constants.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL CDJOEL(GV(I),NBEG,NEND,BvWN,RH(ISTATE),WARN,V,SWF,RM2,
     1                                                          RCNST)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Store molecular constants in array   ZK(v,J,isotope,state)
c
          ZK(I,0,IISTP,ISTATE) = GV(I)
          DO J = 1,MORDRMX
              ZK(I,J,IISTP,ISTATE) = RCNST(J)
              END DO
          WRITE(22,2201) I,GV(I),(RCNST(J),J=1,MORDRMX)
          END DO
 2201 FORMAT(I4,f12.4,f14.10,6(1PD15.7))
c
c** If all is well, then (without further ado) continue with the
c   calculations.
c
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(/'  *** INITDD ERROR ***',/4X,'For state',I3,
     1'  of isotope',I3,'  a serious error has occured:')
  601 FORMAT(4X,'Highest found level is less than desired (v=',I3,',J='
     1,I3,')')
  602 FORMAT(4X,'The Schrodinger Solver was unable to use the initial tr
     1ial energy.')
  603 FORMAT(4X,'The Schrodinger Solver was unable to use the calculated
     1 trial energy')
  604 FORMAT(4X,'The Automatic Level Finder could not find the first vib
     1rational level.')
  605 FORMAT(4X,'The next calculated trial energy was too low.')
  606 FORMAT(4X,'The next calculated trial energy was too high.')
  607 FORMAT(4X,'Could not find vibrational levels beyond (v=',I3,')')
  608 FORMAT(4X,'A second minimum exists in this potential')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DEDP(DATN,ISTATE,IISTP,NPPS,ZMU,KVLEV,JROT,EO,FWHM,
     1                                                   ZK,RM2,DEDPK)
c***********************************************************************
c** This subroutine calculates  dE/dp from the expectation values of the
c   partial derivatives of an analytical potential with respect to its
c   parameters {p(k)}:    dE/dp(k) = <vj|dV/dp(k)|vj>  stored in DEDPK.
c
c** On entry:
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopomer being considered.
c    NPPS(i) are the number of parameters varied in each state.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    KVLEV   is the vibrational quantum number.
c    JROT    is the rotational quantum number.
c    EO      is the initial trial energy (DE in cm-1).
c    RM2     is the inverse distance array  1/r**2
c
c** On exit:
c    EO       is the final calculated energy (DE in cm-1).
c    DEDPK(i) are the values of the partial derivative dE/dP(k).
c
c** Flags: Use only when debugging.
c    INNER specifies wave function matching (& initiation) conditions.
c          = 0 : Match inward & outward solutions at outermost wave
c                function maximum
c          <>0 : Match at inner edge of classically allowed region.
c          < 0 : uses zero slope inner boundary condition.
c          For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
c            To find symmetric eigenfunctions of a symmetric potential,
c            set INNER < 0  & start integration (set RMIN) at potential
c            mid point.
c    IWR   specifies the level of printing inside SCHRQ
c          <> 0 : print error & warning descriptions.
c          >= 1 : also print final eigenvalues & node count.
c          >= 2 : also show end-of-range wave function amplitudes.
c          >= 3 : print also intermediate trial eigenvalues, etc.
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c
      INCLUDE 'arrsizes.h'
c
      INTEGER DATN, bandN, efParity
c
      INTEGER ISTATE, IISTP, NPPS(NSTATEMX), KVLEV, JROT,
     1   KV, JIN, NBEG, NEND, INNER, INNR(0:NVIBMX), IWR, LPRWF,
     2   IPV, I, K
c
      REAL*8 ZMU, EO, BFCT, JFCT, FWHM, VMAX, DEDPV, RM2(NPNTMX),
     1   GV(0:NVIBMX),Vtotal(NPNTMX),JfctL,JfctTQ, JfctTL, muFct,
     2   SWF(NPNTMX), V(NPNTMX), DVDP(NPNTMX), DEDPK(NPARMX)
c
      COMMON /VBLIK/Vtotal
c
c** Vibrational energies and rotational constants.
c
      REAL*8 ZK(0:NVIBMX,0:MORDRMX,NISTPMX,NSTATEMX)
c 
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)
      INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)
      CHARACTER*2 NAME(2)
c
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
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)
c
      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR
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
c
      COMMON /lambdaBLK/nLamdCoef,ifLamdCoef,lamdCoef,lambdaR
c
c** Partial derivative arrays.
c
      REAL*8 DVBODP(NPNTMX,NPOTMX),DUADP(NPNTMX,NPOTMX),
     1  DUBDP(NPNTMX,NPOTMX),DQADP(NPNTMX,NPOTMX),DQBDP(NPNTMX,NPOTMX)
     2  ,dLdP(NPNTMX,NPOTMX),DBDP(NPNTMX,NPOTMX)
c
      COMMON /PDVBLK/DVBODP,DUADP,DUBDP,DQADP,DQBDP,DBDP
      COMMON /PdVLBLK/dLdP
c
c** Initializing the arrays and variables
c
      DATA IWR/0/,LPRWF/0/
c
c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
c
      BFCT = ZMU * RH(ISTATE) * RH(ISTATE) / 16.85762908d0
c
c** Calculating the trial energy value from the vibrational energy 
c   level spacings.
c** For J > 0, also include the rotational distortion constants.
c
      muFct = 16.85762908d0/ZMU
      JFCT = DBLE(JROT*(JROT+1)-IOMEG(ISTATE)*IOMEG(ISTATE))
*      JfctL = DBLE(JROT*(JROT+1))
      JfctL = muFct * DBLE(JROT*(JROT+1))
      JfctTQ = muFct * JFCT 
      JfctTL = muFct * JfctL
      GV(KVLEV) = ZK(KVLEV,0,IISTP,ISTATE)
      IF (JROT.GT.0) THEN
          DO I=1,MORDRMX
              GV(KVLEV) = GV(KVLEV) + JFCT**I*ZK(KVLEV,I,IISTP,ISTATE)
              END DO
          END IF
c
c** Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabatic BOB correction terms.
c
cV(R) = (2*mu/hbar^2)*(VPOT(R)+U(R))+(J(J+1)-OMEGA^2)(1+q(R))/R^2
c       ±(hbar^2/2*mu)(J(J+1))f(R)/R^4
c     = BFCT*(VPOT+U(R))+(JROT(JROT+1)-IOMEGA**2)(1+q(R))*RH**2/R^2
c       ±(hbar^2/2*mu)(JROT(JROT+1))f(R)*RH**2/R^4
c     = BFCT*(VPOT+ZMUA*UAR+ZMUB*UBR)+JFCT(1+ZMQA*QAR+ZMQB*QBR)/R^2
c       + efParity*JfctL*lambdaR/R^4
c
      DO I= 1,NDATPT(ISTATE)
          V(I) = BFCT* (VPOT(I,ISTATE) + ZMUA(IISTP)* UAR(I,ISTATE)
     1                                 + ZMUB(IISTP)* UBR(I,ISTATE))
          Vtotal(I) = V(I)/BFCT
          END DO
c
c** For J > 0, also include the centrifugal distortion and lambda doubling.
c
      JFCT = JFCT * RH(ISTATE) * RH(ISTATE)
      JfctL = JfctL * RH(ISTATE) * RH(ISTATE)
      bandN = IB(DATN)
      IF(ISTATE .EQ. 1) THEN
          efParity = EFPP(DATN)
        ELSEIF (ISTATE .EQ. 2) THEN   
          efParity = EFP(DATN)
        ENDIF
      IF(JROT.GT.0) THEN
          DO I= 1,NDATPT(ISTATE)
          RM2(I) = 1/(RMIN(ISTATE)+RH(ISTATE)*DBLE(I-1))**2
c         if(jrot.eq.1) write(6,666) i,rmin(istate)+rh(istate)*(i-1),
c    1   rm2(i),qar(i,istate),qbr(i,istate)
c 666 format(i3,6f10.5)
              V(I) = V(I) + JFCT* (1.0d0 + ZMQA(IISTP)* QAR(I,ISTATE)
     1                + ZMQB(IISTP)*QBR(I,ISTATE))*RM2(I) 
     2                + efParity* JfctL* lambdaR(I,ISTATE)* RM2(I)**2 
              Vtotal(I) = V(I)/BFCT
              END DO
          END IF
      KV = KVLEV
      INNER = 0
      EO = GV(KVLEV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1      EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2      NBEG,NEND,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated wavefunction is for the wrong vibrational level,
c   then try finding the level in the inner well (if possible).
c
      IF(KV.NE.KVLEV) THEN
          WRITE(6,600) ISTATE, IISTP
          WRITE(6,610) KVLEV, JROT
          IF(KV.GE.0) THEN
              WRITE(6,620) KV, JROT, GV(KVLEV)
              KV = KVLEV
              INNER = 1
              EO = GV(KVLEV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1              EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                                      NBEG,NEND,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ELSE
              WRITE(6,621)
            END IF
          END IF
c
c** If this fails, then try calling ALF to recalculate the
c   vibrational energy from the centrifugal distorted potential.
c
      IF(KV.NE.KVLEV) THEN
          IF(INNER.EQ.1) THEN
              WRITE(6,611) KVLEV, JROT
              IF (KV.GE.0) THEN
                  WRITE(6,620) KV, JROT, GV(KVLEV)
                ELSE
                  WRITE(6,621)
                END IF
              END IF
          WRITE(6,612)
          GV(KVLEV) = 0.0d0
          KV = KVLEV
          JIN = JROT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL ALF(NDATPT(ISTATE),RMIN(ISTATE),RH(ISTATE),V,SWF,
     1      VLIM(ISTATE),KV,JIN,ZMU,EPS(ISTATE),NCN(ISTATE),GV,INNR,IWR)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          KV = KVLEV
          EO = GV(KVLEV)
c
c** Call SCHRQ again with the new calculated energy.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          CALL SCHRQ(KV,JROT,EO,FWHM,VMAX,VLIM(ISTATE),V,SWF,BFCT,
     1              EPS(ISTATE),RMIN(ISTATE),RH(ISTATE),NDATPT(ISTATE),
     2                                  NBEG,NEND,INNR(KVLEV),IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          END IF
c
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning and skip the calculation of the
c   partial derivatives (hence setting them to zero).
c
      IF(KV.NE.KVLEV) THEN
          IF (INNER.EQ.1) THEN
              WRITE(6,611) KVLEV, JROT
              IF (KV.GE.0) THEN
                  WRITE(6,620) KV, JROT, GV(KVLEV)
                ELSE
                  WRITE(6,621)
                END IF
              END IF
          WRITE(6,613)
          WRITE(6,614) KVLEV, JROT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated wavefunction is for the right vibrational level,
c   then continue with calculation of the partial derivatives by
c   integration using the modified trapezoidal rule.
c
c** First determine which partials need to be changed, increment so
c   that the (fluorescence term values and) lower states are skipped.
c
        ELSE
          K = 0
          DO I = 1,ISTATE - 1
              K = K + NPPS(I)
              END DO
c
c** Now calculate the rotational factor JFCT and wavefunction
c   contribution WF^2.
c
          JFCT = JFCT / BFCT
          SWF(NBEG) = 0.5d0 * SWF(NBEG) * SWF(NBEG)
          SWF(NBEG+1) = 23.0d0/24.0d0 * SWF(NBEG+1) * SWF(NBEG+1)
          DO I = NBEG+2,NEND-2
              SWF(I) = SWF(I) * SWF(I)
              END DO
          SWF(NEND-1) = 23.0d0/24.0d0 * SWF(NEND-1) * SWF(NEND-1)
          SWF(NEND) = 0.5d0 * SWF(NEND) * SWF(NEND)
c
c** For DVDP - partials with respect to the each parameter PV(IPV).
c            - including adiabatic correction terms.
c            - including rotational non-adiabatic correction terms.
c            - including lambda doubling correction terms.
          DO IPV = K+1,K+NPPS(ISTATE)
              DO I = NBEG,NEND
                  DVDP(I) = DVBODP(I,IPV) + ZMUA(IISTP)*DUADP(I,IPV)
     1                                    + ZMUB(IISTP)*DUBDP(I,IPV)
     2     + JFCT * (ZMQA(IISTP)*DQADP(I,IPV)+ZMQB(IISTP)*DQBDP(I,IPV))
     3     + efParity * JfctTL * dLdP(I,IPV)
              END DO
c
c** Beginning integration.
c
          DEDPV = 0.0d0
          DO I = NBEG,NEND
              DEDPV = DEDPV + SWF(I) * DVDP(I)
              END DO
          DEDPK(IPV) = DEDPV * RH(ISTATE)
          END DO
        END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT(/'  *** DEDP ERROR ***',/4X,'For state',I3,
     1'  of isotope',I3,'  a serious error has occured:')
  610 FORMAT(4X,'Unable to find v=',I3,' J=',I3,' level in outer well.')
  611 FORMAT(4X,'Unable to find v=',I3,' J=',I3,' level in inner well.')
  612 FORMAT(4X,'Molecular constants give improper trial energy.')
  613 FORMAT(4X,'ALF unable to determine trial energy for this level.')
  614 FORMAT(4X,'Calculation of the v=',I3,' J=',I3,' level has failed.'
     1)
  620 FORMAT(4X,'SCHRQ has obtained the v=',I3,' J=',I3,' level for E=',
     1F11.3,'.')
  621 FORMAT(4X,'SCHRQ was unable to calculate any level.')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE DWDP(DATN,ISTATE,NPPS,EXPT,ZMU,vb,Jr,EO,width,DEDPK,
     1                                                          dWdPk)
c***********************************************************************
c** This subroutine calculates dW/dP using Pajunen's quadrature method 
c   Ref: JCP. 71(6), 2618, 1979.
c** The values of dW/dP(k) for each parameter P(k) are stored in dWdPk.
c********************RJL's Version of 11 May 2002***********************
c
c** On entry:
c    DATN    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopomer being considered.
c    NPPS(i) are the number of parameters varied in each state.
c    ZMU     is the reduced mass of the diatom in atomic units.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
c    width   is the FWHM level width calculated in SCHRQ
c    DEDPK(i) are the values of the partial derivative dE/dP(k).
c    
c** On exit:
c    width   is the width calculated herein (DE in cm-1).
c    dWdPk(i) are the values of the partial derivative dW/dP(k).
c
c***********************************************************************
      INCLUDE 'arrsizes.h'
*** Variables declaration
*----------------------------------------------------------------------*
      INTEGER I,ISTATE,DATN,NN,NST,NPPS(NSTATEMX),vb,Jr,IPV,IPP,j
c ** Phase integral and partial derivative arrays ....
      REAL*8 Iwell(0:NPARMX,-1:1),Ibarr(0:NPARMX,-1:1), DEDPK(NPARMX),
     1       dWdPk(NPARMX)
      REAL*8 ARG,dARG,ZMU,EO,FCTOR,FWHM,Pi,EXPT,width, R1, R2, R3,
     1  KAPPA,KAPPAcl,EMSC,EPSRJ,dKdEPS, XX,TI,COR,dvdG2pi
      DATA pi/3.141592653589793D0/
*----------------------------------------------------------------------*
      FCTOR = SQRT(ZMU/16.85762908d0)
*** Begin by calling subroutine locateTP to find the three turning 
*   points of the quasi-bound level 
      CALL locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
c *** if some turning points are not found, ignore this datum
      IF((R1.LT.0.d0) .OR. (R2.LT.0.d0) .OR. (R3.LT.0.d0)) THEN
          width = EXPT
          DO IPV = IPP+1, IPP+NPPS(ISTATE)
              dWdPk(IPV) = 0.0d0
              END DO 
          WRITE (6,197) DATN
  197 FORMAT(' <<Energy of Datum ',I5,' is above barrier maximum, henc
     &e eliminated>>')
          RETURN
          END IF
c** Call subroutine  phaseIntegral  to calculate the phase integrals 
c   using Pajunen's quadrature method
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R2,R3,-1,EO,DEDPK,Ibarr)
      CALL PhaseIntegral(DATN,ISTATE,NPPS,R1,R2,0,EO,DEDPK,Iwell)
      IPP = 0
      DO j= 1, ISTATE-1
          IPP= IPP + NPPS(j)
          END DO
c ... save WIDTH calculated in SCHRQ to test vs. new value
      FWHM= width
*** Calculate the width
      EPSRJ= 2.0d0*FCTOR*Ibarr(0,-1)
      KAPPAcl= DEXP(-EPSRJ)
      KAPPA= DSQRT(1.d0 + KAPPAcl) - 1.d0
c ... alternate calculation to give better precision for small TUN0
      IF(KAPPAcl.LT.1.d-5) KAPPA= KAPPAcl*(0.5d0- KAPPAcl*(0.125d0-
     1                                                0.0625d0*KAPPAcl))
      KAPPA= 4.d0* KAPPA /(KAPPA+ 2.d0)
c** Derivative of complex gamma function argument calculated as
c .....  EPSRJ= -2.* PI* EMSC
c  per eq.(6.1.27) in Abramowitz and Stegun.
      EMSC= - EPSRJ/(2.d0*pi)
      NST= DABS(EMSC)*1.d2
      NST= MAX0(NST,4)
      ARG= -1.963510026021423d0
      dARG= 0.d0
      DO  I= 0,NST
          NN= I
          XX= I + 0.5d0
          TI= 1.d0/((XX/EMSC)**2 + 1.d0)
          dARG= dARG+ XX*TI*TI
          TI= TI/XX
          ARG= ARG+TI
          IF(DABS(TI).LT.1.d-10) GO TO 233
          ENDDO
c ... and use integral approximation for tails of summations ...
  233 COR= 0.5d0*(EMSC/(NN + 1.d0))**2
      ARG= ARG + COR - COR**2
      dARG= dARG + EMSC**2 *(COR - 2.d0* COR**2)
      dvdG2pi = FCTOR * (Iwell(0,0) 
     1           + (DLOG(DABS(EMSC)) - ARG)* Ibarr(0,0)/(2.d0*pi) )
      width = KAPPA/dvdG2pi
c??????
c      WRITE (32,320) vb,Jr,EO,width, width/FWHM- 1.d0, FWHM
c      WRITE (33,320) vb,Jr,EO,width, width/FWHM- 1.d0, FWHM
c 320  FORMAT ('Level v=',I3,' J=',I3,'  Ep=',G16.8,'    FWHM(new)=',
c    1  1Pd16.8/ 12x,'Relative Diff.=',d9.2,'    FWHM(SCHRQ)=', d16.8)
c??????
      dKdEPS= DSQRT(1.d0 + KAPPAcl)
      dKdEPS= -4.d0*KAPPAcl/(dKdEPS*(dKdEPS+ 1.d0)**2)
      IF( NPPS(ISTATE) .gt. 0 ) THEN
          DO  IPV= 1, NPPS(ISTATE)
              dWdPk(IPP+IPV) = dKdEPS * FCTOR * Ibarr(IPV,0)/dvdG2pi
     1             - (KAPPA/dvdG2pi**2)* (0.5d0 * FCTOR * Iwell(IPV,1)
     2         + (FCTOR**2 *Ibarr(0,0)/(2.d0*pi* EPSRJ))* Ibarr(IPV,0) 
     3                            * (1.d0 - dARG* 8.d0*pi**2/EPSRJ**2)
     4    - (FCTOR* Ibarr(IPV,1)/(4.d0*pi))* (DLOG(DABS(EMSC)- ARG)) )
              END DO
c??????
c         write(32,321)  (dWdPk(IPP+IPV), IPV= 1, NPPS(ISTATE))
c         write(33,321)  (dWdPk(IPP+IPV), IPV= 1, NPPS(ISTATE))
c         write(33,321)  (dEdPk(IPP+IPV), IPV= 1, NPPS(ISTATE))
c         write(33,*)
c 321 format(1Pd11.3,6d11.3/(11x,6d11.3))
c??????
          ENDIF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c*****7***************************************************************72
      SUBROUTINE DeltaBeta(ISTATE,BetaN,deltaP,Change)
c*****7***************************************************************72
*  Subroutine DeltaBeta change the beta parameters manually
*  If Change = .TRUE.  change the beta parameters
*  If Change = .FALSE. recover beta parameters to its original values
c*****7***************************************************************72
      INCLUDE 'arrsizes.h'
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)
c
      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR      
      REAL*8 deltaP,percent
      INTEGER ISTATE,BetaN
      LOGICAL Change
      percent = 0.00001d0
      IF (Change) THEN
          deltaP = BETA(BetaN,ISTATE) * percent
          BETA(BetaN,ISTATE) = BETA(BetaN,ISTATE)*(1.0d0+percent)
        ELSE IF (.NOT.Change) THEN 
          BETA(BetaN,ISTATE) = BETA(BetaN,ISTATE)/(1.0d0+percent) 
          deltaP = BETA(BetaN,ISTATE) * percent
        END IF
******7***************************************************************72
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

******7***************************************************************72
      SUBROUTINE locateTP(DATN,ISTATE,vb,Jr,EO,R1,R2,R3)
******7***************************************************************72
*  Subroutine locateTP locates the turning points R1, R2 and R3 for    
*  a given energy E.                                                   
******7***************************************************************72
*  To find the R1, R2 and R3 for a given energy E, first obtain trial R
*  values by searching the trial energy on the potential array passed 
*  from the common block from RMIN to RMAX, then compare it with the 
*  given energy. This process is repeated until convergence is achieved.
******7***************************************************************72 
*** On entry
c    DATN    is the experimental data number
c    ISTATE  is the molecular state being considered.
c    vb      is the vibrational quantum number.
c    Jr      is the rotational quantum number.
c    EO      is the energy for this level (DE in cm-1).
*** On exit
*  R1          is the array of inner turning points.                 
*  R2          is the array of second turning points.                
*  R3          is the array of third turning points.                 
**********************************************************************72
*  Re          is the Equilibrium Distance for each state.          
*  R1Up        is the upper trial R1.                              
*  R1Down      is the lower trial R1.                                
*  R2Up        is the upper trial R2.                                
*  R2Down      is the lower trial R2.                                
*  R3Up        is the upper trial R3.                                
*  R3Down      is the lower trial R3.                                
*  RH          is the mesh distance used to locate the trial R.   
******7***************************************************************72
      INCLUDE 'arrsizes.h'
c
c** Type statements & common block for data
c
      REAL*8 FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)
c
      INTEGER COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)
c
      CHARACTER*2 NAME(2), SLABL(-3:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
**************************
      REAL*8 RR(1),RM2(NPNTMX),VV(1),VLIMT
      INTEGER NCNN
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)
c
      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR
c
      REAL*8 R(NPNTMX,NSTATEMX),RMAX(NSTATEMX)
      COMMON /RBLK/R
      COMMON /RMAXBLK/RMAX
c
      INTEGER  ISTATE,DATN,vb,Jr,index
c
      REAL*8 R1,R2,R3,R1Up,R1Down,R2Up,R2Down,R3Up,R3Down,EO,
     &  bMi,bMo
      REAL*8 prevEMinusV,currEMinusV,Rmid,Vmid,Vtotal(NPNTMX)
c
      COMMON /VBLIK/Vtotal
c
******7************** Locating inner turning point R1 ******************
      index = 1
  10  IF (EO .LT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN
          WRITE(6,80) vb,Jr
  80      FORMAT('The inner turning point R1 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R1 = -1.0d0
          GO TO 600
        END IF
        GO TO 10
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R1 = R(index,ISTATE)
        GO TO 22
      END IF
      R1Down = R(index,ISTATE)
      R1Up = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index-1) - Vtotal(index) 
  20  Rmid = (R1Up+R1Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
        ELSE
          RR(1) = Rmid
          RM2(1) = 1.0d0 / (RR(1) ** 2)
          CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,
     1                                                        VV,NCNN)
          Vmid = VV(1) 
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
        END IF
      currEMinusV = EO - Vmid
*      WRITE (30,21) R1Up,Rmid,R1Down,EO,Vmid
*  21  FORMAT ('R1U=',F8.6,' Rmid=',F8.6,' R1Down=',F8.6,' E=',G14.7,
*     &' Vmid=',G14.7)
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R1 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R1Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 20 
      ELSE IF (currEMinusV .EQ. 0) THEN
        R1 = Rmid
      ELSE
        R1Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 20
  22  END IF
      WRITE(30,100) vb,Jr,R1,EO,currEMinusV,prevEMinusV
 100  FORMAT('Level v =',I3,' J =',I3,' R1 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8)
******7************** Locating second turning point R2 *****************
      index = IDNINT((RE(ISTATE) - RMIN(ISTATE))/RH(ISTATE))
  25  IF (EO .GT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN    
          WRITE(6,200) vb,Jr
 200      FORMAT('The second turning point R2 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R2 = -1.0d0
          GO TO 600
        END IF
        GO TO 25
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R2 = R(index,ISTATE)
        GO TO 35
      END IF
      R2Up = R(index,ISTATE)
      R2Down = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index) - Vtotal(index-1)
  30  Rmid = (R2Up+R2Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
        ELSE
          RR(1) = Rmid
          RM2(1) = 1.0d0 / (RR(1) ** 2)
          CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,
     1                                                        VV,NCNN)
          Vmid = VV(1)
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
        END IF
      currEMinusV = EO - Vmid
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R2 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R2Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 30
      ELSE IF (currEMinusV .EQ. 0) THEN
        R2 = Rmid
      ELSE 
        R2Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 30
  35  END IF
      WRITE(30,250) vb,Jr,R2,EO,currEMinusV,prevEMinusV
 250  FORMAT('Level v =',I3,' J =',I3,' R2 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8)
******7************** Locating third turning point R3 ******************
  40  IF (EO .LT. Vtotal(index)) THEN
        index = index + 1
        IF (index .GT. NDATPT(ISTATE)) THEN
          WRITE(6,300) vb,Jr
 300      FORMAT('The third turning point R3 of level v =',I3,
     &           ' J =',I3,' is greater than RMAX')
          R3 = -1.0d0
          GO TO 600
        END IF
        GO TO 40
      END IF
      IF (EO .EQ. Vtotal(index)) THEN
        R3 = R(index,ISTATE)
        GO TO 600
      END IF
      R3Down = R(index,ISTATE)
      R3Up = R(index-1,ISTATE)
      prevEMinusV = Vtotal(index-1) - Vtotal(index)   
  50  Rmid = (R3Up+R3Down)/2.0d0
      IF (PSEL(ISTATE) .NE. 0) THEN
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.FALSE.)
        ELSE
          RR(1) = Rmid
          RM2(1) = 1.0d0 / (RR(1) ** 2)
          CALL PREPOT(0,AN(1),AN(2),MN(1,1),MN(2,1),1,RR,RM2,VLIMT,
     1                                                        VV,NCNN)
          Vmid = VV(1)
          CALL vValue(ISTATE,Rmid,Vmid,bMi,bMo,DATN,.TRUE.)
        END IF
      currEMinusV = EO - Vmid
      IF (DABS(currEMinusV) .GE. DABS(prevEMinusV) .AND. 
     &   ((currEMinusV.GE.0 .AND. prevEMinusV.GE.0) .OR. 
     &   (currEMinusV.LE.0 .AND. prevEMinusV.LE.0))) THEN
        R3 = Rmid
      ELSE IF (currEMinusV .GT. 0) THEN
        R3Down = Rmid
        prevEMinusV = currEMinusV
        GO TO 50
      ELSE IF (currEMinusV .EQ. 0) THEN
        R3 = Rmid
      ELSE 
        R3Up = Rmid
        prevEMinusV = currEMinusV
        GO TO 50
      END IF
      WRITE(30,350) vb,Jr,R3,EO,currEMinusV,prevEMinusV
 350  FORMAT('Level v =',I3,' J =',I3,' R3 =',F8.5,' E =',G16.8,
     &  ' EmVc =', G16.8,' EmVp =', G16.8/)
 600  RETURN
******7***************************************************************72
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c*****7***************************************************************72
      SUBROUTINE PhaseIntegral(DATN,ISTATE,NPPS,Ri,Ro,k1,EO,DEDPK,Int)
c*****7***************************************************************72
c  Subroutine PhaseIntegral calculates phase integrals of the integrand
c  f(R)/[E-V(R)]**(k+1/2), for k=-1,0 or 1  on interval between turning 
c  points  Ri and Ro, where f(R) is either 1 or a derivative of the
c  potential w.r.t. some potential parameter. 
c             Ro                                         N
c      Int = Int dr f(r)/|E-V(r)|^{k+1/2} = ½(Ro - Ri) Sum{Wi*F(Zi)}
c             Ri                                        i=1
c  Evaluate integrals for  k=k1 to 1  where for k=k1 f(r)=1
c                                           for k>k1 f(r)=[dE/dp-dV/dp]
c********************RJL's Version of 11 May 2002***********************
c  Tests against proper Gaussian (25.4.38 & 25.4.40) gave agreement of
c    ca. 10^{-10} to 10^{-8) for I_{-1}^{barr}({1})  and  ranging from 
c    10^{-13} to [rare worst case] 10^{-6} fo I_0^{barr}({deriv})
c    typically 1-^{-6} for I_0^{well}({1}).
c* For k=1 integrals, comparing PP "k=1" vs. "k=3" integration gave
c  typical agreement for I_1^{barr}({deriv.}) of ca. 10^{-6}, but
c  occasionally as bad as 10^{-1} and as good as 10^{-9}
c   [?? check cgce. of turning point search!].
c** I_1^{well}({f})  derivatives always poor cgce. w.r.t. k=1 vs. k=3
c               PP method integration tests ... ???
c* Notice deriv.'s w.r.t. q(r) parameters typically orders of magnitude
c     poorer than others (??)
c
c* On entry:  DATN    is the experimental datum number
c-----------  ISTATE  is the molecular state being considered.
c             NPPS(i) are the number of parameters varied in each state.
c             Ri      is the inner turning point
c             Ro      is the outer turning point
c             EO      is the energy for this level (DE in cm-1).
c             k       is the powers in the phase integral Int
c             n=0  for  f(r)= 1;  n=1  for  f(r)= [dV(r)/dp_k - dE/dp_k]
c             DEDPK(i) are the values of the partial derivative dE/dp_k.
c* On exit:  Int(j,k)    are the phase integral(s)
c=======================================================================
      INCLUDE 'arrsizes.h'
c** Type statements & common block for data
      INTEGER i,IPP,IPV,ISTATE, j,k,k1,kmx,M,DATN, NPPS(NSTATEMX)
c** M  is the number of quadrature mesh points
      PARAMETER (M = 21)
c
      REAL*8 dVdPk(NPOTMX)
      COMMON /dVdPkBLK/dVdPk
c
      REAL*8 Ri,Ro,EO,Fzi,Fzt,Rd,Vd,EMinusV,EMinusVt,
     1  Int(0:NPARMX,-1:1),DEDPK(NPARMX), ARG,Pi
c
      REAL*8 Zp21(M),Wp21k1(M),Wp21k3(M)
ccc  1   ,Zp11(M),Wp11k1(M),Wp11k3(M),Zp15(M),Wp15k1(M),Wp15k3(M)
ccc
ccc   DATA Zp11/0.9659258262890683, 0.8660254037844386,
ccc  1          0.7071067811865475, 0.5000000000000000,
ccc  2          0.2588190451025208, 
ccc  3          0.0000000000000000,
ccc  4         -0.2588190451025208,-0.5000000000000000,
ccc  5         -0.7071067811865475,-0.8660254037844386,
ccc  6         -0.9659258262890683,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp11k1/22.06633397656787,-37.69911184307752,
ccc  1            35.60471674068432,-37.69911184307752,
ccc  2            36.57672889044161,
ccc  3           -37.69911184307752,
ccc  4            36.57672889044161,-37.69911184307752,
ccc  5            35.60471674068432,-37.69911184307752,
ccc  6            22.06633397656787,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp11k3/228.3214545745810, -720.4719152232592,
ccc  1           1139.3509357018984,-1357.1680263507907,
ccc  2           1447.1946273399755,
ccc  3          -1474.4541520848097,
ccc  4           1447.1946273399755,-1357.1680263507907,
ccc  5           1139.3509357018984, -720.4719152232592,
ccc  6            228.3214545745810,
ccc  7          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Zp15/0.9807852804032304, 0.9238795325112868, 
ccc  1          0.8314696123025452, 0.7071067811865475,
ccc  2          0.5555702330196022, 0.3826834323650898,
ccc  3          0.1950903220161283, 
ccc  4          0.0000000000000000,
ccc  5         -0.1950903220161283,-0.3826834323650898,
ccc  6         -0.5555702330196022,-0.7071067811865475,
ccc  7         -0.8314696123025452,-0.9238795325112868,
ccc  8         -0.9807852804032304,
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp15k1/29.62981929591175,-50.26548245743668,
ccc  1            47.72092686124880,-50.26548245743664,
ccc  2            49.12943331558201,-50.26548245743658,
ccc  3            49.44900912828539,
ccc  4           -50.26548245743656,
ccc  5            49.44900912828539,-50.26548245743658,
ccc  6            49.12943331558201,-50.26548245743664,
ccc  7            47.72092686124880,-50.26548245743668,
ccc  8            29.62981929591175, 
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
ccc   DATA Wp15k3/1164.639428963841,-3580.432803129281,
ccc  1            5525.620597073791,-6534.512719466765,
ccc  2            7020.542275282852,-7276.911407677031,
ccc  3            7400.700330802901,
ccc  4           -7439.291403700618,
ccc  5            7400.700330802901,-7276.911407677031,
ccc  6            7020.542275282852,-6534.512719466765,
ccc  7            5525.620597073791,-3580.432803129281,
ccc  8            1164.639428963841,
ccc  9          0.0,0.0,0.0,0.0,0.0,0.0/
ccc
      DATA Zp21/0.9898214418809327, 0.9594929736144974,
     1          0.9096319953545184, 0.8412535328311812,
     2          0.7557495743542583, 0.6548607339452851,
     3          0.5406408174555976, 0.4154150130018864,
     4          0.2817325568414297, 0.1423148382732851,
     5          0.0000000000000000,
     6         -0.1423148382732851,-0.2817325568414297,
     7         -0.4154150130018864,-0.5406408174555976,
     8         -0.6548607339452851,-0.7557495743542583, 
     9         -0.8412535328311812,-0.9096319953545184,
     a         -0.9594929736144974,-0.9898214418809327/
c
      DATA Wp21k1/40.91258980361040,-69.11503837897816,
     1            65.80507790523560,-69.11503837898373,
     2            67.78308420797106,-69.11503837899778,
     3            68.30792716563759,-69.11503837900795,
     4            68.49459295516724,-69.11503837901213,
     5            68.54383971474920,
     6           -69.11503837901213, 68.49459295516724,
     7           -69.11503837900795, 68.30792716563759,
     8           -69.11503837899778, 67.78308420797106,
     9           -69.11503837898373, 65.80507790523560,
     a           -69.11503837897816, 40.91258980361040/
ccc
ccc   DATA Wp21k3/6364.91821744174,-19267.83229098791,
ccc  1           29317.26172550868,-34478.42376106038,
ccc  2           37049.70046340271,-38499.30026119479,
ccc  3           39357.74970048209,-39887.54536375314,
ccc  4           40205.64256060925,-40378.03411697539,
ccc  5           40431.72625305376,
ccc  6          -40378.03411697539, 40205.64256060925,
ccc  7          -39887.54536375314, 39357.74970048209,
ccc  8          -38499.30026119479, 37049.70046340271,
ccc  9          -34478.42376106038, 29317.26172550868,
ccc  a          -19267.83229098791,  6364.91821744174/
c
      DATA Pi/3.141592653589793D0/
ccc   SAVE Pi, Zp21, Wp21k1, Wp21l3
******7***************************************************************72
* In Pajunen method, there is a option to choose different points of 
* weight. To do this, just change the variable names for the Zp* & Wp*
* e.g., in order to use only 15 points
*                   Zi(i) = Zp15(i) and Wi(i) = -Wp15k1(i), 
* and to use only 11 points 
*                   Zi(i) = Zp11(i) and Wi(i) = -Wp11k1(i)
******7***************************************************************72  
      IPP = 0
c ... determine parameter count for previous states
      IF(ISTATE.GT.1) THEN
          DO j = 1, ISTATE-1
              IPP = IPP + NPPS(j)
              END DO
          ENDIF
c*** Zero integral arrays
      kmx= k1
      IF(k1 .eq. -1) kmx= 0
      IF( NPPS(ISTATE) .GT. 0 ) kmx= 1
      DO  k= k1,kmx
          DO  IPV= 0, MAX(0, NPPS(ISTATE))
              Int(IPV,k) = 0.0d0
              ENDDO
          ENDDO
c
c*** Begin quadrature loop for sums over M mesh points
      DO  i= 1,M
c ... first get potential and derivatives at Pajunen k=1 point
          Rd = 0.5d0*(Ro+Ri + (Ro-Ri)*Zp21(i))
          CALL VGEN(ISTATE,Rd,Vd,DATN)
          EMinusV = DABS(EO- Vd)
          DO  k= k1, kmx
c ... loop over posible  k  values
              IF(k .EQ. -1) THEN
c* For k= -1, using Pajunen k=1 quadrature method ......
                  Fzt= DSQRT(EminusV/(1-Zp21(i)**2))* (1-Zp21(i)**2)**2
                ELSEIF(k .EQ. 0) THEN 
c* For k= 0, using Pajunen k=1 quadrature method ......
                  Fzt= DSQRT((1-Zp21(i)**2)/EminusV) * (1-Zp21(i)**2)
c* For k= 1, use Pajunen quadrature method (JCP, 71(6), 2618, 1979)
                ELSEIF(k .EQ. 1) THEN
c* If  k = 1 points are used
                  Fzt= ( DSQRT((1-Zp21(i)**2)/EminusV) )**3
                END IF
c*** Accumulate the phase integrals
              Int(0,k) = Int(0,k) - Wp21k1(i)* Fzt
              IF((k .gt. k1). and. ( NPPS(ISTATE) .gt. 0 )) THEN
                  DO IPV = 1, NPPS(ISTATE)
c ... Loop over free parameters accumulating phase integral
                      Int(IPV,k) = Int(IPV,k) - Wp21k1(i)* Fzt
     1                            * (dVdPk(IPP+IPV)- DEDPK(IPP+IPV))
                      ENDDO
                  END IF
c .... end of loop over  k values
              ENDDO
c ... end of quadrature loop over mesh points ....................
          ENDDO
      DO  k= k1,kmx
          DO  IPV= 0, MAX(0, NPPS(ISTATE))
c*** In Pajunen's method, a factor of (1/2) is incorporated because 
c    the phase integral is a contour integral
              Int(IPV,k) = 0.25d0*(Ro-Ri) * Int(IPV,k)
              ENDDO
          ENDDO
c????????????????????????
c** Test that 'regular' & Pajunen integrals agree for k= -1 & 0
c   and that Pajunen  k=1 & 3 integrals agree for k=1 case
ccc   DO  IPV= 0, MAX(0, NPPS(ISTATE))
ccc       WRITE(32, 320) EO, (Int(IPV,k), k= k1, kmx)
c         WRITE(32, 321)  (Int(IPV,k), k= k1, kmx)
c         Do  k= k1, kmx
c             tst(k)= 0.d0
c             IF(DABS(Ink(IPV,k)) .gt. 0.d0) 
c    1                              tst(k)= Int(IPV,k)/Ink(IPV,k)-1.d0
c             enddo
c         write(32, 322) (tst(k), k= k1,kmx)
ccc       ENDDO
cc320 FORMAT('  At   E=',F12.6,1P3D16.8)
cc321 FORMAT(21x,1P3D16.8)
cc322 FORMAT(21x,3(1PD10.2,6x))
c????????????????????????
      RETURN
      END      
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE vValue(ISTATE,Rd,Vd,bMi,bMo,DATN,Forward)
c**********************************************************************
c** This subroutine will generate potential value Vd given a known 
c   distance variable Rd
c   Only the following types of potential form are included in current
c   version
c     - The Long-Range Expanded Morse Oscillator
c     - Modified Lennard-Jones Oscillator
c     - Total Expanded Morse Oscillator
c     - Long-Range Total Expanded Morse Oscillator (trial version)
c     - Double-Exponential Long-Range Oscillator
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 2001 by Yiye Huang
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c
c    ISTATE  is the electronic state being fitted to.
c
c    Rd      is the distance, if Rd > 0, calculate that point potential
c            value Vd
c    DATN    is the date number
c    Forward if Forward = .TRUE., do forward calculation, Vd is an input  
c** On entry via common blocks:
c    PSEL   is used to choose the type of potential used.
c           |PSEL| = 0 : Do forward calculation.
c           |PSEL| = 1 : Use the Generalized Morse Potential.
c           |PSEL| = 2 : Use the Expanded Morse Potential.
c           |PSEL| = 3 : Use the Modified Morse Potential.
c           |PSEL| = 4 : Use the Modified Lennard-Jones Potential.
c           |PSEL| = 5 : Use the Generalized Potential Energy Function.
c           |PSEL| = 6 : Use the Expanded Morse Long-Range Potential. 
c           |PSEL| = 7 : Use the Modified Lennard-Jones-Gaussian Potential. 
c           |PSEL| = 8 : Use the Total Expanded Morse Potential.
c           |PSEL| = 9 : Total Expanded Morse Long-Range Potential(trail ver.).
c           |PSEL| = 10: Use the Total Expanded Morse Long-Range Potential.
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
c** On exit via common blocks:
c    Vd     is the potential value that is generated.
c
******7*****************************************************************
      INCLUDE 'arrsizes.h'

      INTEGER ISTATE, IPV, I, J, NPPS(NSTATEMX),l
      COMMON /NPPSBLK/NPPS 

      REAL*8 BTEMP,BINF,RDVAL,RTEMP,RM2,XTEMP,XITEMP,PBTEMP,PETEMP,
     1  AZERO,AONE,ATWO,DBETA(0:2),AZTEMP,AOTEMP,ATTEMP,ARTEMP,FSW,
     2  Xtemp2,Btemp2,betaCoef,BetaInfinity,bMi,bMo,A,B, BMtemp,ZTEMP

      LOGICAL emoPrinted,Forward
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

      REAL*8 Uad,Ubd,Qad,Qbd,Ld,Rd,Vd,dVdPk(NPOTMX),JFCT,JfL
      INTEGER DATN, NBAND, IISTP, efParity, NDatPtS

c** For Expanded Morse Long-Range Potential
      INTEGER numCn(NSTATEMX),nLR(NSTATEMX,numCnMax),tempN,Df
      REAL*8 dampBeta(NSTATEMX),Cn(NSTATEMX,numCnMax),Demo(NSTATEMX),
     &  DeOffset,ReOffset,Dn(NSTATEMX,numCnMax,NPNTMX),Vlr,Dm,Dmp,Dmpp,
     &  ReBar,fdr

      COMMON /DELRBLK/numCn,nLR,dampBeta,Cn,ReBar,Df
      COMMON /OffSetBLK/Demo,DeOffset,ReOffset

c** The possible repulsive chemical bonding
      REAL*8 C(NSTATEMX),DVDC(NPNTMX,NSTATEMX)
      INTEGER ifC(NSTATEMX)
      COMMON /repulsiveBLK/C,ifC,DVDC
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
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Initializing variables.
c
      Uad = 0.0d0
      Ubd = 0.0d0
      Qad = 0.0d0
      Qbd = 0.0d0
      Ld  = 0.0d0
      PBTEMP = 0.0d0
      PETEMP = 0.0d0
      bMi = 0.0d0
      bMo = 0.0d0
      RDVAL = Rd
      IF (Forward) THEN
        Vd = Vd + VLIM(ISTATE) - DE(ISTATE)
        GO TO 444
      END IF
      Vd = 0.0d0
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.2) THEN
        IF (PoEpZ(ISTATE) .LT. 2) THEN
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
        ELSE
          RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        END IF
        bMi = BETA(0,ISTATE)
        DO 80 J = 1,NBETA(ISTATE)
          bMi = bMi + BETA(J,ISTATE)*RTEMP**J
   80   CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
        IF (DATN.LT.IRSR(ISTATE)) THEN
          bMi = BMORSE(IRSR(ISTATE),ISTATE)
        END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
        XTEMP = DEXP(-bMi*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
        Vd = DE(ISTATE)*(1.0d0-XTEMP)**2
     1                   + VLIM(ISTATE) - DE(ISTATE)
        bMo = bMi
      ENDIF
******7***** End Expanded Morse Oscillator Potential Function **********

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Modified Lennard-Jones Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.4) THEN
        IF (PoEpZ(ISTATE) .LT. 2) THEN
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
        ELSE
          RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        END IF
        bMi = BETA(0,ISTATE)
        DO J = 1,nBetaIn(ISTATE)
          bMi=bMi+BETA(J,ISTATE)*RTEMP**J
        END DO
        bMo = bMi
        DO J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
          bMo = bMo+BETA(J,ISTATE)*RTEMP**J
        END DO
        IF (PoEpZ(ISTATE) .EQ. 3) THEN
          betaCoef = 2.0d0 * RE(ISTATE)**nPB(ISTATE) /
     &               (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
          BETA(NBETA(ISTATE)+1,ISTATE) = DLOG(2.0d0*DE(ISTATE)*
     &                   RE(ISTATE)**NCN(ISTATE)/CNVAL(ISTATE))
          BetaInfinity = BETA(NBETA(ISTATE)+1,ISTATE)
          bMi = bMi * betaCoef
          bMo = BETA(NBETA(ISTATE)+1,ISTATE) * RTEMP +
     &      bMo * betaCoef
        END IF
c-----------------------------------------------------------------------
c** If chosen, constrain the Cn value to CNVAL by using a switching
c   function.
c
c           BETAsw(R) = fsw(R) [ BETA(R) - BINF ] + BINF
c
c              fsw(R) = 1 / [e^(ALPHA'*(R - Rx)) + 1]
c
        IF (PoEpZ(ISTATE).NE.3 .AND. CNVAL(ISTATE).GT.0.0d0) THEN
          BINF = DLOG(2.0d0*DE(ISTATE)*RE(ISTATE)**NCN(ISTATE)
     1           /CNVAL(ISTATE))
          BetaInfinity = BINF
          FSW = 1.0d0 / (1.0d0+DEXP(ALPHA(ISTATE)*(RDVAL-RX(ISTATE))))
c
c** Update the generalized morse exponent.
c
          bMo = FSW*(bMo - BINF) + BINF
*        ENDIF
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
     1       DLOG(2.0d0*DE(ISTATE)
     2      *RE(ISTATE)**NCN(ISTATE)/DABS(CNVAL(ISTATE)))
          BetaInfinity = BINF + BETA(NBETA(ISTATE)+1,ISTATE)
c
c** Update the generalized morse exponent.
c
          bMo = bMo
     1      +BETA(NBETA(ISTATE)+1,ISTATE)*RTEMP**(NBETA(ISTATE)+1)
        ENDIF
c**********************************

c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < IRSR.
c
        IF (DATN.LT.IRSR(ISTATE)) THEN
          bMo = BMORSE(IRSR(ISTATE),ISTATE)
        END IF
c-----------------------------------------------------------------------
c** Calculate some temporary variables.
c
        IF (PoEpZ(ISTATE).NE.3) THEN
          IF (inAndOut(ISTATE) .EQ. 0) THEN
            XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     1        DEXP(-bMo*RTEMP)
            BMtemp = bMo
          ELSE
            IF (RDVAL .LT. RE(ISTATE)) THEN
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     &          DEXP(-bMi*RTEMP)
              BMtemp = bMi
            ELSE
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
     1          DEXP(-bMo*RTEMP)
              BMtemp = bMo
            END IF
          END IF
        ELSE
          IF (inAndOut(ISTATE) .EQ. 0) THEN
            XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
c    1        DEXP(-bMo)
     &        DEXP(-bMo*RTEMP)
          ELSE
            IF (RDVAL .LT. RE(ISTATE)) THEN
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
c    1          DEXP(-bMi)
     &          DEXP(-bMi*RTEMP)
            ELSE
              XTEMP = (RE(ISTATE)/RDVAL)**NCN(ISTATE) *
c    1          DEXP(-bMo)
     &          DEXP(-bMo*RTEMP)
            END IF
          END IF
        END IF
c
c** Now to calculate the actual potential
c
        Vd = DE(ISTATE)*(1.0d0-XTEMP)**2 + VLIM(ISTATE) - DE(ISTATE)
      END IF
******7***** End Modified Lennard-Jones Potential Function *************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Longe-Range Potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.6) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        IF (PoEpZ(ISTATE) .NE. 2) THEN 
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
        ELSE 
          RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        END IF
        bMi = BETA(0,ISTATE)
        DO 380 J = 1,NBETA(ISTATE)
          bMi = bMi+BETA(J,ISTATE)*RTEMP**J
  380   CONTINUE
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
        XTEMP = DEXP(-bMi*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
        Vlr = 0.0d0
        DO J = 1,numCn(ISTATE)
          tempN = nLR(ISTATE,J)
          CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
          Dn(ISTATE,nLR(ISTATE,J),I) = Dm
          Vlr = Vlr + Cn(ISTATE,J)/RDVAL**nLR(ISTATE,J)*Dm
        END DO
        Vd = DE(ISTATE)*(1.0d0-XTEMP)**2 - DeOffset 
     1                   + Vlr + VLIM(ISTATE) - Demo(ISTATE)
        bMo = bMi
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
******7***** End Expanded Morse Longe-Range Potential Function *********

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Total Expanded Morse Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.8) THEN
        IF (PoEpZ(ISTATE) .LT. 2) THEN
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
        ELSE
          RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        END IF
        bMi = BETA(0,ISTATE)
        DO 523 J = 1,nBetaIn(ISTATE)
          bMi = bMi+BETA(J,ISTATE)*RTEMP**J
  523   CONTINUE
        bMo = bMi
        DO 533 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
          bMo = bMo + BETA(J,ISTATE)*RTEMP**J
  533   CONTINUE
c-----------------------------------------------------------------------
c** If chosen, constrain BETA(R) at BETA(IRSR) for R < RSR.
c
        IF (DATN.LT.IRSR(ISTATE)) THEN
          bMi = bMorsIn(IRSR(ISTATE),ISTATE)
        END IF
c-----------------------------------------------------------------------
c
c** Calculate some temporary variables.
c
        XTEMP = DEXP(-bMo * (RDVAL-RE(ISTATE)))
        Xtemp2= DEXP(-2.0d0*bMi*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
        Vd = DE(ISTATE)*(1.0d0 + Xtemp2 - 2.0d0*XTEMP)
     &                 + VLIM(ISTATE) - DE(ISTATE)
      ENDIF
******7***** End Total Expanded Morse Potential Function ***************

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Long-Range Total Expanded Morse Oscillator (trial version).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.9) THEN
        DE(ISTATE) = DE(ISTATE) + DeOffset
        RE(ISTATE) = RE(ISTATE) + ReOffset
        IF (PoEpZ(ISTATE) .LT. 2) THEN
          RTEMP = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
        ELSE
          RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        END IF
        bMi = BETA(0,ISTATE)
        DO 623 J = 1,nBetaIn(ISTATE)
          bMi = bMi + BETA(J,ISTATE)*RTEMP**J
  623   CONTINUE
        bMo = bMi
        DO 633 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
          bMo = bMo + BETA(J,ISTATE)*RTEMP**J
  633   CONTINUE
c
c** Calculate some temporary variables.
c
        XTEMP = DEXP(-bMo * (RDVAL-RE(ISTATE)))
        Xtemp2= DEXP(-2.0d0*bMi*(RDVAL-RE(ISTATE)))
c
c** Now to calculate the actual potential
c
        Vlr = 0.0d0
        DO J = 1,numCn(ISTATE)
          tempN = nLR(ISTATE,J)
          CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
          Dn(ISTATE,nLR(ISTATE,J),I) = Dm
          Vlr = Vlr + Cn(ISTATE,J)/RDVAL**nLR(ISTATE,J)*Dm
        END DO
        Vd = DE(ISTATE)*(1.0d0 + Xtemp2 - 2.0d0*XTEMP)
     &      - DeOffset + Vlr + VLIM(ISTATE) - Demo(ISTATE)
        DE(ISTATE) = DE(ISTATE) - DeOffset
        RE(ISTATE) = RE(ISTATE) - ReOffset
      ENDIF
******7* End Total Expanded Morse Long-Range Potential Function (Trial)*

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Double-Exponential Long-Range Oscillator.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** Now to calculate the generalized morse exponent
c
      IF (PSEL(ISTATE).EQ.10) THEN
        RTEMP = (RDVAL**nPB(ISTATE)-RE(ISTATE)**nPB(ISTATE))/
     &            (RDVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))
        bMi = BETA(0,ISTATE)
        DO 723 J = 1,nBetaIn(ISTATE)
          bMi = bMi + BETA(J,ISTATE)*RTEMP**J
  723   CONTINUE
        bMo = bMi
        DO 733 J = nBetaIn(ISTATE)+1,NBETA(ISTATE)
          bMo = bMo + BETA(J,ISTATE)*RTEMP**J
  733   CONTINUE
c
c** Calculate some temporary variables.
c
        IF (inAndOut(ISTATE) .EQ. 0) THEN
          XTEMP = DEXP(-bMo * (RDVAL-RE(ISTATE)))
          Xtemp2= DEXP(-2.0d0*bMi*(RDVAL-RE(ISTATE)))
        ELSE
          IF (RDVAL .LT. RE(ISTATE)) THEN
            XTEMP = DEXP(-bMi*(RDVAL-RE(ISTATE)))
            Xtemp2= DEXP(-2.0d0*bMi*(RDVAL-RE(ISTATE)))
          ELSE
            XTEMP = DEXP(-bMo * (RDVAL-RE(ISTATE)))
            Xtemp2= DEXP(-2.0d0 * bMo * (RDVAL-RE(ISTATE)))
          END IF
        END IF

c
c** Now to calculate the actual potential, (see notes for the math)
c
        Vlr = 0.0d0
        A = 0.0d0
        B = 0.0d0
        DO J = 1,numCn(ISTATE)
          tempN = nLR(ISTATE,J)
          CALL dampF(RDVAL,dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
          Dn(ISTATE,nLR(ISTATE,J),I) = Dm
          Vlr = Vlr + Cn(ISTATE,J)/RDVAL**nLR(ISTATE,J)*Dm
          CALL dampF(RE(ISTATE),dampBeta(ISTATE),tempN,Dm,Dmp,Dmpp,Df)
          A = A + Cn(ISTATE,J)/RE(ISTATE)**nLR(ISTATE,J)*Dm*
     &          (1.0d0+Dmp/(BETA(0,ISTATE)*Dm)-nLR(ISTATE,J)
     &          /(BETA(0,ISTATE)*RE(ISTATE)))
          B = B + Cn(ISTATE,J)/RE(ISTATE)**nLR(ISTATE,J)*Dm*
     &          (2.0d0+Dmp/(BETA(0,ISTATE)*Dm)-nLR(ISTATE,J)
     &          /(BETA(0,ISTATE)*RE(ISTATE)))
        END DO
        IF (Df .EQ. 4) THEN
          fdr = 1.0d0/(1.0d0+DEXP(-1.3*dampBeta(ISTATE)*(RDVAL-ReBar)))
          Vlr = fdr * Vlr
        END IF
        Vlr=Vlr+C(ISTATE)*DEXP(-BETA(0,ISTATE)/2.0d0*(RDVAL-RE(ISTATE)))
        A = A + DE(ISTATE)
        B = B + 2.0d0*DE(ISTATE)
        Vd=A*Xtemp2-B*XTEMP+DE(ISTATE)+Vlr+VLIM(ISTATE)-Demo(ISTATE)
      END IF
******7* End Total Expanded Morse Long-Range Potential Function ********

      IF (DATN .EQ. 0) GO TO 999
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Born-Oppenheimer breakdown correction terms
c** and the lambda doubling correction terms.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 444  RM2 = 1/RDVAL**2
c
c** Choose Born-Oppenheimer correction expansion parameter.
c
      IF (PoEpZ(ISTATE).LE.0) THEN
        IF (PoEpZ(ISTATE).EQ.0) RTEMP = (RDVAL-RE(ISTATE))
     1                               /(RDVAL+RE(ISTATE))
        IF (PoEpZ(ISTATE).EQ.-1) RTEMP = (RDVAL-RE(ISTATE))
c
c** Adiabatic correction terms and state dependent partial derivatives.
c
        Uad = 0.0d0
        DO 250 J = 0,NUA(ISTATE)
          Uad = Uad + UA(J,ISTATE)*RTEMP**J
  250   CONTINUE

        Ubd = 0.0d0
        DO 260 J = 0,NUB(ISTATE)
          Ubd = Ubd + UB(J,ISTATE)*RTEMP**J
  260   CONTINUE
c
c** Rotational non-adiabatic correction terms and state dependent
c   partial derivatives.
c
        Qad = 0.0d0
        DO 270 J = 0,NQA(ISTATE)
          Qad = Qad + QA(J,ISTATE) * RTEMP**J
  270   CONTINUE
       
        Qbd = 0.0d0
        DO 280 J = 0,NQB(ISTATE)
          Qbd = Qbd + QB(J,ISTATE) * RTEMP**J
  280   CONTINUE
        
c** For lambda doubling
        Ld = 0.0d0
        DO 300 J = 0,nLamdCoef(ISTATE)
          Ld = Ld + lamdCoef(J,ISTATE)/(RDVAL+RE(ISTATE))**J 
  300   CONTINUE
******7***************************************************************72
**  for PoEpZ > 0, introduce a new expression for BOB and lambda 
**  doubling functions
******7***************************************************************72
      ELSE IF (PoEpZ(ISTATE).GT.0) THEN
        RTEMP = 1.0d0/(RDVAL**nPu(ISTATE)+RE(ISTATE)**nPu(ISTATE))
        ZTEMP = RDVAL**nPu(ISTATE)-RE(ISTATE)**nPu(ISTATE)
        Uad = UA(0,ISTATE)*2*RE(ISTATE)**nPu(ISTATE)*RTEMP
     &    + UA(NUA(ISTATE),ISTATE)*ZTEMP*RTEMP
        DO 500 J = 1,NUA(ISTATE)-1
          Uad = Uad + UA(J,ISTATE)
     &      *2*RE(ISTATE)**nPu(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
  500   CONTINUE

        Ubd = UB(0,ISTATE)*2*RE(ISTATE)**nPu(ISTATE)*RTEMP
     &    + UB(NUB(ISTATE),ISTATE)*ZTEMP*RTEMP
        DO 520 J = 1,NUB(ISTATE)-1
          Ubd = Ubd + UB(J,ISTATE)
     &      *2*RE(ISTATE)**nPu(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
  520   CONTINUE
**
** Note: the 1/R^2 factor in the Hamiltonian has been incorporated into 
**       the partial derivatives
**
        RTEMP = 1.0d0/(RDVAL**nPq(ISTATE)+RE(ISTATE)**nPq(ISTATE))
        ZTEMP = RDVAL**nPq(ISTATE)-RE(ISTATE)**nPq(ISTATE)
        Qad = QA(0,ISTATE)*2*RE(ISTATE)**nPq(ISTATE)*RTEMP
        DO 540 J = 1,NQA(ISTATE)
          Qad = Qad + QA(J,ISTATE)
     &      *2*RE(ISTATE)**nPq(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
  540   CONTINUE

        Qbd = QB(0,ISTATE)*2*RE(ISTATE)**nPq(ISTATE)*RTEMP
        DO 560 J = 1,NQB(ISTATE)
          Qbd = Qbd + QB(J,ISTATE)
     &      *2*RE(ISTATE)**nPq(ISTATE)*RTEMP*(ZTEMP*RTEMP)**J
  560   CONTINUE

c** For   lambda doubling
        RTEMP = 1.0d0/(RDVAL**nPL(ISTATE)+RE(ISTATE)**nPL(ISTATE))
        ZTEMP = RDVAL**nPL(ISTATE)-RE(ISTATE)**nPL(ISTATE)
        Ld = 0.0d0
        DO 590 J = 0,nLamdCoef(ISTATE)
          Ld = Ld + lamdCoef(J,ISTATE) * (ZTEMP*RTEMP)**J
  590   CONTINUE
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
      NBAND = IB(DATN)
      IISTP = ISTP(NBAND)
      JFCT = 16.85762908d0/ZMASS(3,IISTP)
     &     * DBLE(JPP(DATN)*(JPP(DATN)+1)-IOMEG(ISTATE)*IOMEG(ISTATE))
*      JfL=(16.85762908d0/ZMASS(3,IISTP))*DBLE(JPP(DATN)*(JPP(DATN)+1))
      JfL = (16.85762908d0/ZMASS(3,IISTP))**2 * DBLE(JPP(DATN) 
     &    * (JPP(DATN)+1))
      Vd = Vd + ZMUA(IISTP)*Uad + ZMUB(IISTP)*Ubd      
c** For J > 0, also include the centrifugal distortion and lambda doubling.
      IF (ISTATE .EQ. 1) THEN
        efParity = EFPP(DATN)
      ELSEIF (ISTATE .EQ. 2) THEN
        efParity = EFP(DATN)
      ENDIF
      IF (JPP(DATN).GT.0) THEN
        Vd = Vd + JFCT*(1.0d0 + ZMQA(IISTP)*Qad + ZMQB(IISTP)*Qbd)/Rd**2
     &       + efParity * JfL * Ld/Rd**4
      END IF
******7********************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
