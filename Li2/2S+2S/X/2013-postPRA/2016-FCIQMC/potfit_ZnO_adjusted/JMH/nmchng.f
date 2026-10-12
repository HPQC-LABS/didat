c***********************************************************************
      SUBROUTINE NMCHNG(NSTATES,NPTOT,PV,FORBAC)
c***********************************************************************
c** This subroutine will convert internal NLLSSRR parameter array back
c   into external (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 1997  by  Jenning Seto       (ver of 12/10/1997) +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    NSTATES is the number of states being considered.
c    NPTOT   is the total number of parameters being varied.
c    NPPS
c    PV(i)   is the NLLSSRR parameter array.
c    FORBAC  is a flag to determine which way the parameters are mapped
c            FORBAC = 0 : Map internal PV to external variables.
c            FORBAC = 1 : Map external varuables to internal PV.
c
      INCLUDE 'arrsizes.h'

      INTEGER NSTATES, NPTOT, NPPS(NSTATEMX), FORBAC, ISTATE, IPV, I
      COMMON /NPPSBLK/NPPS

      REAL*8 PV(NPTOT)
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

c** The possible repulsive chemical bonding
      REAL*8 C(NSTATEMX),DVDC(NPNTMX,NSTATEMX)
      INTEGER ifC(NSTATEMX)
      COMMON /repulsiveBLK/C,ifC,DVDC
c***********************************************************************
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

c** For Modified Lennard-Jones-Gaussian Potential
      INTEGER ifAmp(NSTATEMX),ifR0(NSTATEMX),ifFWHM(NSTATEMX)
      REAL*8 Amp(NSTATEMX),R0(NSTATEMX),FWHM(NSTATEMX)
      COMMON /MLJG/ifAmp,ifR0,ifFWHM,Amp,R0,FWHM

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generating internal NLLSSRR array (PV) from external variables and
c   calculating the number of "free" parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF (FORBAC.EQ.0) THEN
        NPTOT = 0
        DO ISTATE = 1,NSTATES
          NPPS(ISTATE) = 0
          IF (PSEL(ISTATE).NE.5) THEN
            IF (IFDE(ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = DE(ISTATE)
            ENDIF
          ENDIF
          IF (IFRE(ISTATE).LE.0) THEN
            NPPS(ISTATE) = NPPS(ISTATE) + 1
            PV(NPTOT+NPPS(ISTATE)) = RE(ISTATE)
          ENDIF
          IF (PSEL(ISTATE).EQ.7) THEN
            IF (ifAmp(ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = Amp(ISTATE)
            ENDIF
            IF (ifR0(ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = R0(ISTATE)
            ENDIF
            IF (ifFWHM(ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = FWHM(ISTATE)
            ENDIF
          ENDIF
          DO I = 0,NBETA(ISTATE)
            IF (IFBETA(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = BETA(I,ISTATE)
            ENDIF
          END DO
          DO I = 0,NUA(ISTATE)
            IF (IFUA(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = UA(I,ISTATE)
            ENDIF
          END DO
          DO I = 0,NUB(ISTATE)
            IF (IFUB(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = UB(I,ISTATE)
            ENDIF
          END DO
          DO I = 0,NQA(ISTATE)
            IF (IFQA(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = QA(I,ISTATE)
            ENDIF
          END DO
          DO I = 0,NQB(ISTATE)
            IF (IFQB(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = QB(I,ISTATE)
            ENDIF
          END DO
          DO I = 0, nLamdCoef(ISTATE)
            IF (ifLamdCoef(I,ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = lamdCoef(I,ISTATE)
            ENDIF
          END DO
          IF (PSEL(ISTATE).EQ.10) THEN
            IF (ifC(ISTATE).LE.0) THEN
              NPPS(ISTATE) = NPPS(ISTATE) + 1
              PV(NPTOT+NPPS(ISTATE)) = C(ISTATE) 
            END IF
          END IF
          NPTOT = NPTOT + NPPS(ISTATE)
        END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Convert internal NLLSSRR parameter array back into external
c   (logical) variable system.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSEIF (FORBAC.EQ.1) THEN
        IPV = 0
        DO ISTATE=1,NSTATES
          IF (ABS(PSEL(ISTATE)).NE.5) THEN
            IF (IFDE(ISTATE).LE.0) THEN
              IPV = IPV + 1
              DE(ISTATE) = PV(IPV)
            ENDIF
          ENDIF
          IF (IFRE(ISTATE).LE.0) THEN
            IPV = IPV + 1
            RE(ISTATE) = PV(IPV)
          ENDIF
          IF (PSEL(ISTATE).EQ.7) THEN
            IF (ifAmp(ISTATE).LE.0) THEN
              IPV = IPV + 1
              Amp(ISTATE) = PV(IPV)
            ENDIF
            IF (ifR0(ISTATE).LE.0) THEN
              IPV = IPV + 1
              R0(ISTATE)  = PV(IPV)
            ENDIF
            IF (ifFWHM(ISTATE).LE.0) THEN
              IPV = IPV + 1
              FWHM(ISTATE) = PV(IPV)
            ENDIF
          ENDIF
         DO I=0,NBETA(ISTATE)
            IF (IFBETA(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              BETA(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          DO I=0,NUA(ISTATE)
            IF (IFUA(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              UA(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          DO I=0,NUB(ISTATE)
            IF (IFUB(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              UB(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          DO I=0,NQA(ISTATE)
            IF (IFQA(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              QA(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          DO I=0,NQB(ISTATE)
            IF (IFQB(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              QB(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          DO I=0,nLamdCoef(ISTATE)
            IF (ifLamdCoef(I,ISTATE).LE.0) THEN
              IPV = IPV + 1
              lamdCoef(I,ISTATE) = PV(IPV)
            ENDIF
          END DO
          IF (PSEL(ISTATE).EQ.10) THEN
            IF (ifC(ISTATE).LE.0) THEN
              IPV = IPV + 1
              C(ISTATE) = PV(IPV)
            END IF
          END IF
        END DO
      ENDIF
      RETURN
      END
