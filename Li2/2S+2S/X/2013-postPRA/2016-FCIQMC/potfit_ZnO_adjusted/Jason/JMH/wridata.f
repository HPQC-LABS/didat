c***********************************************************************
      SUBROUTINE WRIDATA(NPASS,NSTATES,SLABL,NISTP,NAME,DECM,PU,PS)
c***********************************************************************
c** Subroutine to print out and correct calculated values for the
c   parameter-expansion fits to standard output.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 1997  by  Jenning Seto       (ver of 07/01/1999) +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NPASS   is the number of times WRIDATA was called.
c    NSTATES is the number of states being considered.
c    NISTP   is the number of isotopomers being considered.
c    NAME    is the name of the molecule.
c    PU      are the parameter uncertainties.
c    PS      are the parameter sensitivities.
c
      INCLUDE 'arrsizes.h'

      CHARACTER*2  NAME(2)

      INTEGER NPASS, NSTATES, NISTP, ISTATE, IPARM, I

      REAL*8 DECM(NSTATEMX), BTEMP, UAT, UBT
c
c** NLLSSRR variables used in output
c
      REAL*8 PU(NPARMX), PS(NPARMX)
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      CHARACTER*2 SLABL(-3:NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      INTEGER PoEpZ(NSTATEMX),
     &  nPB(NSTATEMX),nPq(NSTATEMX),nPu(NSTATEMX),nPL(NSTATEMX)

      COMMON /nPwBLK/nPB,nPq,nPu,nPL

      COMMON /PoEpZBLK/PoEpZ

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
c** Born-Oppenheimer Breakdown Parameters.
c
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NQA(NSTATEMX),NQB(NSTATEMX),
     1  IFUA(0:NABOMX,NSTATEMX),IFUB(0:NABOMX,NSTATEMX),
     3  IFQA(0:NRNAMX,NSTATEMX),IFQB(0:NRNAMX,NSTATEMX)

      REAL*8 UA(0:NABOMX,NSTATEMX),UB(0:NABOMX,NSTATEMX),
     1  QA(0:NRNAMX,NSTATEMX),QB(0:NRNAMX,NSTATEMX)

      COMMON /BPBLK/NUA,NUB,NQA,NQB,IFUA,IFUB,IFQA,IFQB,UA,UB,QA,QB
**
** Lambda doubling function
**
      INTEGER nLamdCoef(NSTATEMX),ifLamdCoef(0:NRNAMX,NSTATEMX)
      REAL*8 lamdCoef(0:NRNAMX,NSTATEMX),lambdaR(NPNTMX,NSTATEMX)

      COMMON /lambdaBLK/nLamdCoef,ifLamdCoef,lamdCoef,lambdaR
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Write out energy convergence criterion (same for all states).
c-----------------------------------------------------------------------
      IF (NPASS.EQ.1) WRITE(6,690) EPS(1)
c-----------------------------------------------------------------------
c** Writing out state specific information.
c-----------------------------------------------------------------------
      IPARM = 0
      DO ISTATE=1,NSTATES
        WRITE(6,600)
        IF (PSEL(ISTATE).EQ.1) THEN
          WRITE(6,601) SLABL(ISTATE),NAME(1),NAME(2)
        ELSEIF (PSEL(ISTATE).EQ.2) THEN
          WRITE(6,602) SLABL(ISTATE),NAME(1),NAME(2)
        ELSEIF (PSEL(ISTATE).EQ.3) THEN
          WRITE(6,603) SLABL(ISTATE),NAME(1),NAME(2)
        ELSEIF (PSEL(ISTATE).EQ.4) THEN
          WRITE(6,604) SLABL(ISTATE),NAME(1),NAME(2)
          WRITE(6,605) NCN(ISTATE)
        ELSEIF (PSEL(ISTATE).EQ.5) THEN
          WRITE(6,606) SLABL(ISTATE),NAME(1),NAME(2)
          WRITE(6,607) AGPEF(ISTATE),BGPEF(ISTATE),KGPEF(ISTATE)
        ELSEIF (PSEL(ISTATE).EQ.6) THEN
          WRITE(6,608) SLABL(ISTATE),NAME(1),NAME(2)
          WRITE (6,545) 
          IF (PoEpZ(ISTATE) .EQ. 2) THEN 
            WRITE(6,555) nPB(ISTATE),nPB(ISTATE),nPB(ISTATE),nPB(ISTATE)
            WRITE(6,557) nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
     &                  ,nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
     &                  ,nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
            WRITE(6,559) nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
     &                  ,nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
     &                  ,nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
          ELSE
            WRITE (6,565)
          END IF 
        ELSEIF (PSEL(ISTATE).EQ.7) THEN
          WRITE(6,609) SLABL(ISTATE),NAME(1),NAME(2)
        ELSEIF (PSEL(ISTATE).EQ.10) THEN
          WRITE(6,613) SLABL(ISTATE),NAME(1),NAME(2)
          WRITE (6,575) cRp(ISTATE)
          IF (PoEpZ(ISTATE) .EQ. 2) THEN
            WRITE(6,555) nPB(ISTATE),nPB(ISTATE),nPB(ISTATE),nPB(ISTATE)
            WRITE(6,557) nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
     &                  ,nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
     &                  ,nPu(ISTATE),nPu(ISTATE),nPu(ISTATE),nPu(ISTATE)
            WRITE(6,559) nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
     &                  ,nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
     &                  ,nPq(ISTATE),nPq(ISTATE),nPq(ISTATE),nPq(ISTATE)
          ELSE
            WRITE (6,565)
          END IF
        ENDIF
        
c
c** Write out potential constraints for each state (if any).
c
        IF (NPASS.EQ.1) WRITE(6,614)
        IF (NPASS.EQ.2) WRITE(6,615)
c-----------------------------------------------------------------------
c** Writing out the Ae information.
c-----------------------------------------------------------------------
        WRITE(6,630) 'Ae',VLIM(ISTATE)
c-----------------------------------------------------------------------
c** Writing out the Te information.
c-----------------------------------------------------------------------
        IF (ISTATE.GT.1) THEN
          UAT = 0.0d0
          UBT = 0.0d0
          IF (IFDE(1).EQ.0) UAT = PU(1)
          IF (IFDE(ISTATE).EQ.0) UBT = PU(IPARM+1)
          UBT = (UAT+UBT)*DSQRT(DECM(ISTATE)+1.0d0)
          UAT = DE(1)-VLIM(1)+VLIM(ISTATE)-DE(ISTATE)
          WRITE(6,620) 'Te',UAT,UBT,0.0d0
        END IF
c-----------------------------------------------------------------------
c** Writing out the De information.
c-----------------------------------------------------------------------
        IF ((PSEL(ISTATE).NE.5).AND.(IFDE(ISTATE).EQ.0)) THEN
          IPARM = IPARM + 1
          IF (DABS(DE(ISTATE)).GT.PU(IPARM)) THEN
            WRITE(6,620) 'De',DE(ISTATE),PU(IPARM),PS(IPARM)
          ELSE
            WRITE(6,621) 'De',DE(ISTATE),PU(IPARM),PS(IPARM)
          ENDIF
        ELSE
          WRITE(6,630) 'De',DE(ISTATE)
        ENDIF
c-----------------------------------------------------------------------
c** Writing out the Re information.
c-----------------------------------------------------------------------
        IF (IFRE(ISTATE).EQ.0) THEN
          IPARM = IPARM + 1
          IF (DABS(RE(ISTATE)).GT.PU(IPARM)) THEN
            WRITE(6,620) 'Re',RE(ISTATE),PU(IPARM),PS(IPARM)
          ELSE
            WRITE(6,621) 'Re',RE(ISTATE),PU(IPARM),PS(IPARM)
          ENDIF
        ELSE
          WRITE(6,630) 'Re',RE(ISTATE)
        ENDIF
c** For Modified Lennard-Jones-Gaussian Potential
        IF (PSEL(ISTATE).EQ.7) THEN
c-----------------------------------------------------------------------
c** Writing out the Gaussian function maximum amplitude.
c-----------------------------------------------------------------------
          IF (ifAmp(ISTATE).LE.0) THEN
            IPARM = IPARM + 1
            IF (DABS(Amp(ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,624) 'Amp',Amp(ISTATE),PU(IPARM),PS(IPARM)
            ELSE
              WRITE(6,625) 'Amp',Amp(ISTATE),PU(IPARM),PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,635) 'Amp',Amp(ISTATE)
          ENDIF
c-----------------------------------------------------------------------
c** Writing out the Gaussian function R0.
c-----------------------------------------------------------------------
          IF (ifR0(ISTATE).LE.0) THEN
            IPARM = IPARM + 1
            IF (DABS(R0(ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,620) 'R0',R0(ISTATE),PU(IPARM),PS(IPARM)
            ELSE
              WRITE(6,621) 'R0',R0(ISTATE),PU(IPARM),PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,630) 'R0',R0(ISTATE)
          ENDIF
c-----------------------------------------------------------------------
c** Writing out the Gaussian function FWHM.
c-----------------------------------------------------------------------
          IF (ifFWHM(ISTATE).LE.0) THEN
            IPARM = IPARM + 1
            IF (DABS(FWHM(ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,626) 'FWHM',FWHM(ISTATE),PU(IPARM),PS(IPARM)
            ELSE
              WRITE(6,627) 'FWHM',FWHM(ISTATE),PU(IPARM),PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,636) 'FWHM',FWHM(ISTATE)
          ENDIF
        ENDIF
c-----------------------------------------------------------------------
c** Writing out the exponent parameter information.
c-----------------------------------------------------------------------
        DO I=0,NBETA(ISTATE)
          IF (IFBETA(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (DABS(BETA(I,ISTATE)).GT.PU(IPARM)) THEN
             WRITE(6,640) 'Be','ta',I,BETA(I,ISTATE),PU(IPARM),PS(IPARM)
            ELSE
             WRITE(6,641) 'Be','ta',I,BETA(I,ISTATE),PU(IPARM),PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,650) 'Be','ta',I,BETA(I,ISTATE)
          ENDIF
        END DO
c
c** Write out Cn constant for the MLJ.
c
        IF (PSEL(ISTATE).EQ.4) THEN
          WRITE(6,650) ' ','C ',NCN(ISTATE),CNVAL(ISTATE)
        END IF
c-----------------------------------------------------------------------
c** Writing out the Adiabatic information for atom A.
c-----------------------------------------------------------------------
        UAT = 0.0D0
        DO I=0,NUA(ISTATE)
          IF (IFUA(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (DABS(UA(I,ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,640) ' u',NAME(1),I,UA(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ELSE
              WRITE(6,641) ' u',NAME(1),I,UA(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,650) ' u',NAME(1),I,UA(I,ISTATE)
          ENDIF
          UAT = UAT + UA(I,ISTATE)
        END DO
c-----------------------------------------------------------------------
c** Writing out the Adiabatic information for atom B.
c-----------------------------------------------------------------------
        UBT = 0.0D0
        DO I=0,NUB(ISTATE)
          IF (IFUB(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (DABS(UB(I,ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,640) ' u',NAME(2),I,UB(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ELSE
              WRITE(6,641) ' u',NAME(2),I,UB(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,650) ' u',NAME(2),I,UB(I,ISTATE)
          ENDIF
          UBT = UBT + UB(I,ISTATE)
        END DO
c-----------------------------------------------------------------------
c** Writing out the Rotational Non-Adiabatic information for atom A.
c-----------------------------------------------------------------------
        DO I=0,NQA(ISTATE)
          IF (IFQA(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (NISTP.EQ.1) THEN
              IF (DABS(QA(I,ISTATE)).GT.PU(IPARM)) THEN
                WRITE(6,640) ' q',' ',I,QA(I,ISTATE),PU(IPARM),PS(IPARM)
              ELSE
                WRITE(6,641) ' q',' ',I,QA(I,ISTATE),PU(IPARM),PS(IPARM)
              ENDIF
            ELSE
              IF (DABS(QA(I,ISTATE)).GT.PU(IPARM)) THEN
                WRITE(6,640) ' q',NAME(1),I,QA(I,ISTATE),PU(IPARM),
     1                       PS(IPARM)
              ELSE
                WRITE(6,641) ' q',NAME(1),I,QA(I,ISTATE),PU(IPARM),
     1                       PS(IPARM)
              ENDIF
            ENDIF
          ELSE
            IF (NISTP.EQ.1) THEN
              WRITE(6,650) ' q',' ',I,QA(I,ISTATE)
            ELSE
              WRITE(6,650) ' q',NAME(1),I,QA(I,ISTATE)
            ENDIF
          ENDIF
        END DO
c-----------------------------------------------------------------------
c** Writing out the Rotational Non-Adiabatic information for atom B.
c-----------------------------------------------------------------------
        DO I=0,NQB(ISTATE)
          IF (IFQB(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (DABS(QB(I,ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,640) ' q',NAME(2),I,QB(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ELSE
              WRITE(6,641) ' q',NAME(2),I,QB(I,ISTATE),PU(IPARM),
     1                     PS(IPARM)
            ENDIF
          ELSE
            WRITE(6,650) ' q',NAME(2),I,QB(I,ISTATE)
          ENDIF
        END DO
c-----------------------------------------------------------------------
c** Writing out lambda doubling coefficients
c-----------------------------------------------------------------------
        DO I=0,nLamdCoef(ISTATE)
          IF (ifLamdCoef(I,ISTATE).EQ.0) THEN
            IPARM = IPARM + 1
            IF (DABS(lamdCoef(I,ISTATE)).GT.PU(IPARM)) THEN
              WRITE(6,642) 'Lambda',I,lamdCoef(I,ISTATE),PU(IPARM),
     &                     PS(IPARM)
            ELSE
              WRITE(6,643) 'Lambda',I,lamdCoef(I,ISTATE),PU(IPARM),
     &                     PS(IPARM)
            ENDIF
          ELSE         
            WRITE(6,651) 'Lambda',I,lamdCoef(I,ISTATE)
          ENDIF
        END DO
      END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Write out the unmodified final parameter values to fort.9.
c
      IF (NPASS.EQ.2) THEN
        DO ISTATE = 1,NSTATES 
          IF (PSEL(ISTATE).EQ.10) THEN
            WRITE(9,*) C(ISTATE)
          END IF
          WRITE(9,*) DE(ISTATE)
          WRITE(9,*) RE(ISTATE)
          IF (PSEL(ISTATE).EQ.7) THEN
            WRITE(9,*) Amp(ISTATE)
            WRITE(9,*) R0(ISTATE)
            WRITE(9,*) FWHM(ISTATE) 
          ENDIF
          DO I = 0,NBETA(ISTATE)
            WRITE(9,*) BETA(I,ISTATE)
          END DO
          DO I = 0,NUA(ISTATE)
            WRITE(9,*) UA(I,ISTATE)
          END DO
          DO I = 0,NUB(ISTATE)
            WRITE(9,*) UB(I,ISTATE)
          END DO
          DO I = 0,NQA(ISTATE)
            WRITE(9,*) QA(I,ISTATE)
          END DO
          DO I = 0,NQB(ISTATE)
            WRITE(9,*) QB(I,ISTATE)
          END DO
          DO I = 0,nLamdCoef(ISTATE)
            WRITE(9,*) lamdCoef(I,ISTATE)
          END DO
        END DO
      ENDIF
      WRITE(6,600)
      RETURN
c-----------------------------------------------------------------------
  545 FORMAT(3x,'V(R) = De[1 - exp(-beta(z)*(R - Re))]^2 + Vlr(R)',/,
     &1x,'where Vlr(R) = SUM(Cn*Dn(R)/R^n) over inverse-power n',/,
     &1x,'in which Dn(R) is the damping function, and',/,
     &3x,'beta(z) = beta0 + beta1*z + ... + betaN*z^N') 
  555 FORMAT(3x,'z = (R^',I1,' - Re^',I1,')/(R^',I1,' + Re^',I1,')')
  557 FORMAT(1x,'Adiabatic BOB is expanded as',/,
     &3x,'UA(R) = u0*(2Re/(R^',I1,'+Re^',I1,')) + uN*((R^',I1,'-Re^'
     &,I1,')/(R^',I1,'+Re^',I1,')) ',/,
     &9x,'+ (2Re/(R^',I1,'+Re^',I1,'))SUM(ui*((R^',I1,'-Re^',I1,')/(R^',
     &I1,'+Re^',I1,'))^i from 1 to N-1',/,
     &1x,'UB(R) is expanded in the same way')
  559 FORMAT(1x,'Non-adiabatic BOB is expanded as',/,
     &3x,'QA(R) = q0*(2Re/(R^',I1,'+Re^',I1,')) + qN*((R^',I1,'-Re^'
     &,I1,')/(R^',I1,'+Re^',I1,')) ',/,
     &9x,'+ (2Re/(R^',I1,'+Re^',I1,'))SUM(qi*((R^',I1,'-Re^',I1,')/(R^',
     &I1,'+Re^',I1,'))^i from 1 to N-1',/,
     &1x,'QB(R) is expanded in the same way')
  565 FORMAT(1x,'z = (R - Re)/(R + Re)')
  575 FORMAT(3x,'V(R) = A*exp(-2*beta(z)*(R - Re)) - B*exp(-beta(z)*(R -
     & Re)) + De + Vlr(R)',/,
     &1x,'where Vlr(R) = C*exp(-',F3.1,'*beta0*(R - Re)) + SUM(Cn*Dn(R)/
     &R^n) over inverse-power n',/,1x,'in which Dn(R) is the damping fun
     &ction, and',/,3x,'beta(z) = beta0 + beta1*z + ... + betaN*z^N')
  600 FORMAT(/,1X,39('=='))
  601 FORMAT(/' The potential used for state',A3,'  of the molecule ',A2
     1,A2,/' is the Generalized Morse Oscillator Model for Diatoms.')
  602 FORMAT(/' The potential used for state',A3,'  of the molecule ',A2
     1,A2,/' is the Expanded Morse Oscillator Model for Diatoms.')
  603 FORMAT(/' The potential used for state',A3,'  of the molecule ',A2
     1,A2,/' is the Modified Morse Oscillator Model for Diatoms.')
  604 FORMAT(/' The potential used for state',A3,'  of the molecule ',A2
     1,A2,/' is the Modified Lennard-Jones Oscillator Model for Diatoms.
     2')
  605 FORMAT(' The potential tail decreases as 1/R^',I1)
  606 FORMAT(' The potential used for state',A3,'  of the molecule ',A2,
     1A2/,' is the Generalized Potential Energy Function for Diatoms.')
  607 FORMAT(' The parameters used are: a =',F6.3,', b =',F6.3,', and k
     1=',I3)
  608 FORMAT(' The potential used for state',A3,'  of the molecule ',A2,
     1A2/,' is the Long-range Expanded Morse potential for Diatoms.')
  609 FORMAT(' The potential used for state',A3,'  of the molecule ',A2,
     1A2/,' is the Modified Lennard-Jones-Gaussian potential for Diatoms
     2.')
  613 FORMAT(' The potential used for state',A3,' of the molecule ',A2,A
     &2/,' is the Double-Exponential Long-Range potential for Diatoms.') 
  614 FORMAT(/'  Parameter:   Initial Value:   Uncertainty:  Sensitivity
     1:')
  615 FORMAT(/'  Parameter:    Final Value:    Uncertainty:  Sensitivity
     1:')
  620 FORMAT(6X,A2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  621 FORMAT(5X,'*',A2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  624 FORMAT(5X,A3,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  625 FORMAT(4X,'*',A2,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  626 FORMAT(4X,A4,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  627 FORMAT(3X,'*',A4,3X,1PD19.10,3X,1PD8.1,6X,1PD8.1)
  630 FORMAT(6X,A2,3X,1PD19.10,7X,'--',12X,'--')
  635 FORMAT(5X,A3,3X,1PD19.10,7X,'--',12X,'--')
  636 FORMAT(4X,A4,3X,1PD19.10,7X,'--',12X,'--')
  640 FORMAT(3X,A2,A2,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  642 FORMAT(1X,A6,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
*  645 FORMAT(3X,A3,'(',I3,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  641 FORMAT(2X,'*',A2,A2,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  643 FORMAT('*',A6,'(',I2,')',1PD19.10,3X,1PD8.1,6X,1PD8.1)
  650 FORMAT(3X,A2,A2,'(',I2,')',1PD19.10,7X,'--',12X,'--')
  651 FORMAT(1X,A6,'(',I2,')',1PD19.10,7X,'--',12X,'--')
*  660 FORMAT('  For isotope ',I2,' the adiabatic correction is:',F19.10,
*     1' wavenumbers.')
  690 FORMAT(/' Setting the Schrodinger Solver energy convergence criter
     1ion to  ',1PD7.1,'  cm-1')
*  691 FORMAT(/' Constrain the C(',I2,') coefficient to be',1PD17.10)
      END
