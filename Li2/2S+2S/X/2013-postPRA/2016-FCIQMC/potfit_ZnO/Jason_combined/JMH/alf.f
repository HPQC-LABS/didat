c** Version 0.9x dated August 15, 2000.
c***********************************************************************
      SUBROUTINE ALF(NDP,RMIN,RH,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,NCN,GV,
     1               INNR,IWR) 
c***********************************************************************
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++   COPYRIGHT 1998 - 1999  by  Jenning Seto and Robert J. Le Roy   +++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Please inform me of any bugs, by phone at: (519)888-4567, ext. 4099 +
c++++++++ by e-mail to: jyseto@uwaterloo.ca , or write me at: ++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Based on the automatic level finding routine found in LEVEL 6.0
c   written by Robert J. Le Roy
c** Uses the Schrodinger solver subroutine SCHRQ.
c
c** On entry:
c    NDP    is the number of datapoints used for the potential.
c    RMIN   is the inner radial distance of the potential (ang).
c    RH     is the meshvalue (ang).
c           NDP, RMIN, and RH define the radial range over which the
c           potential is defined.
c    V(i)   is the scaled input potential (cm-1).
c           The scaling factor BFCT is (2*mu/hbar^2)*RH^2.
c    VLIM   is the potential asymptote (cm-1).
c    KVMAX  is the maximum vibrational level for which we wish to find.
c    AFLAG  is the rotational state of the potential.
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    NCN    is the near dissociation limit radial exponential.
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   returns the highest allowed vibrational quantum number if
c            less than the inputed KVMAX.
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : Subroutine has functioned normally.
c             = -1 : KVMAX larger than number of allowed levels.
c             = -2 : Initial trial energy is unusable.
c             = -3 : Calculated trial energy is unusable.
c             = -4 : Cannot find first vibrational level.
c             = -5 : Calculated trial energy too low.
c             = -6 : Calculated trial energy too high.
c             = -7 : An impossible situation occured.
c             = -8 : Potential found to have a second minimum.
c    GV(v)   contains the vibrational energys for each level.
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
c    MCO   specifies the level of printing of molecular constants.
c          >  0 : print out vibrational energies to channel-21.
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
c    LPRWF specifies option of printing out generated wavefunction
c          > 0 : print wave function every LPRWF-th  point.
c          < 0 : compactly write to channel-7 every |LPRWF|-th wave
c                function value.
c          A lead "card" identifies the level, gives the position of
c          1-st point and radial mesh, & states No. of  points.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The dimensioning parameters must be consistant with the sizes of the
c   arrays used in the calling program.
c
c    NVIBMX  is the maximum number of vibrational levels considered.
c            Note: NVIBMX should be larger than KVMAX.
c
      INTEGER NVIBMX
      PARAMETER (NVIBMX = 400)
c
c** NF counts levels found in automatic search option
c
c** OWL holds the vibrational levels that are contained in the outer
c   well.
c** IWL holds the vibrational levels that are contained in the inner
c   well (if present).
c
      INTEGER NDP,KVMAX,NCN,KV,AFLAG,NF,NBEG,NEND,INNR(0:KVMAX),IWR,
     1  I,IZPE,IVDIF,IVCOR,IQT,IEG,LTRY,TRIALE,AWO,MCO,INNER,LPRWF,JROT,
     2  NPMIN, NPMAX, NIWL,IWL(0:NVIBMX),NOWL,OWL(0:NVIBMX)
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,
     1  GV(0:KVMAX),BV(0:NVIBMX),BVDOUT,BVDIN,AO,VD,
     2  BZ,BFCT,PW,PWI,GAMA,VMIN,VMAX,RE,PMAX,VDMV,VDL,VDU,
     3  VPMIN(10), RPMIN(10), VPMAX(10), RPMAX(10),
     4  LHIE, LLOE,
     5  ZQ, ZH, Z1, Z2
c
      DATA AWO/-1/,MCO/2/,LPRWF/0/
c
      DATA ZQ/0.25D0/,ZH/0.5D0/,Z1/1.D0/,Z2/2.D0/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check if the dimensions are adequate.
c
      IF (KVMAX.GT.NVIBMX) THEN
        WRITE(6,610)
        WRITE(6,613) KVMAX, NVIBMX
        STOP
      END IF
c
c** Initialize level counters for each well.
c
      DO I = 0,KVMAX
        INNR(I) = -1
        IWL(I) = 0
        OWL(I) = 0
      END DO
c
c** Initialize the remaining variables and flags.
c
      NF = 0
      NIWL = 0
      NOWL = 0
      KV = 0
      INNER = 0
      LTRY = 0
      LHIE = 1.0d99
      CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** Store rotational quantum number.
c
      JROT = AFLAG
c
c** Numerical factor  16.85762908  based on 1998 physical constants.
c
      BZ = ZMU/16.85762908d0
      BFCT = BZ*RH*RH
c
c** RMAX is the outer radial distance over which the potential is
c   defined. 
c
      RMAX = RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
c
      NPMIN = 0
      DO I = 2,NDP-1
        IF ((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
          NPMIN = NPMIN + 1
          RPMIN(NPMIN) = RMIN + DBLE(I-1)*RH
          VPMIN(NPMIN) = V(I) / BFCT
          IF (NPMIN.EQ.10) GOTO 100
        END IF
      END DO
c
c** If a minimum cannot be found, then print a warning and exit.
c
  100 IF (NPMIN.EQ.0) THEN
        WRITE(6,610)
        WRITE(6,614)
        STOP
      END IF
c
c** If more than two minima are found, then print a warning and exit.
c
      IF (NPMIN.GT.2) THEN
        WRITE(6,605)
        WRITE(6,615) NPMIN, 'minima'
*        STOP
      END IF
c
c** Locate the potential maxima (if it exists).
c
      NPMAX = 0
      DO I = 2,NDP-1
        IF ((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
          NPMAX = NPMAX + 1
          RPMAX(NPMAX) = RMIN + DBLE(I-1)*RH
          VPMAX(NPMAX) = V(I) / BFCT
          IF (NPMAX.EQ.10) GOTO 150
        END IF
      END DO
c
c** If no maxima were found, then set the maximum to be the value at the
c   end of the range.
c
  150 IF (NPMAX.EQ.0) THEN
        NPMAX = 1
        RPMAX(NPMAX) = RMAX
        VPMAX(NPMAX) = V(NDP) / BFCT
      END IF
c
c** If more than three maxima are found, then print a warning and exit.
c
      IF (NPMAX.GT.3) THEN
        WRITE(6,605)
        WRITE(6,615) NPMAX, 'maxima'
*        STOP
      END IF
c
c** If there is no rotationless barrier to assiciation, then set the
c   final VPMAX to be the value at the end of the range.
c
      IF (RPMAX(NPMAX).LT.RPMIN(NPMIN)) THEN
        NPMAX = NPMAX + 1
        RPMAX(NPMAX) = RMAX
        VPMAX(NPMAX) = V(NDP) / BFCT
      END IF
c
c** If a maxima occurs before a minima, then the potential turns over in
c   short range region and should not be used. Print a warning and exit.
c
      IF (RPMAX(1).LT.RPMIN(1)) THEN
        WRITE(6,610)
        WRITE(6,616) RPMAX(1)
*        STOP
      END IF
c
c** Now find the absolute potential minimum.
c
      VMIN = VPMIN(1)
      RE = RPMIN(1)
      DO I = 2,NPMIN
        IF (VMIN.GT.VPMIN(I)) THEN
          VMIN = VPMIN(I)
          RE = RPMIN(I)
        END IF
      END DO 
c
c** Now find the absolute potential maximum.
c
      VMAX = VPMAX(1)
      DO I = 2,NPMAX
        IF (VMAX.LT.VPMAX(I)) VMAX = VPMAX(I)
      END DO 
c
c** If the absolute potential maximum is lower than the absolute
c   potential minimum, then print out an error statement and quit.
c
      IF (VMAX.LE.VMIN) THEN
        WRITE(6,610)
        WRITE(6,617)
        STOP
      END IF
c
c** Otherwise, print out the results if desired.
c
      IF (AWO.GT.0) THEN
        WRITE(6,650) NPMIN, VMIN
        WRITE(6,651) NPMAX, VMAX
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculate 2*NCN/(NCN - 2) for use when calculating trial energies. 
c
      PW = 20.0d0
      IF ((NCN.GT.0).AND.(NCN.NE.2)) PW = Z2*DBLE(NCN)/(DBLE(NCN)-Z2)
      IF (VMAX.GT.VLIM) PW = Z2
      PWI = Z1/PW
c
c** Use Lennard-Jones estimation of zero point energy to determine the
c   initial trial energy.
c                            _____________________________
c            vD + 0.5 = ao \/ZMU * De * Re^2 / 16.85762908
c
c                 De = A (vD - v)^3 = A (vD + 0.5)^3
c
c              E(v=0) = VMIN + A [(vD + 0.5)^3 - vD^3]
c
c** Choose AO to have a value of 0.25.
c
      AO = ZQ
      VD = AO * DSQRT(BZ*(VMAX-VMIN)) * RE - ZH
      AO = (VMAX-VMIN)*(Z1 - (VD/(VD+ZH))**3)
      EO = VMIN + AO
c
c** If desired, write out energy level information.
c
      IF (MCO.GE.1) THEN
        WRITE(21,2100)
        WRITE(21,2110) RMIN, RMAX, RH, BZ, ZMU
        WRITE(21,2111) EPS
        WRITE(21,2112)
        WRITE(21,2101)
      END IF
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc up to the KVMAXth vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
   10 IF (AWO.GT.0) THEN
        WRITE(6,601)
        IF (INNER.EQ.0) WRITE(6,602)
        IF (INNER.EQ.1) WRITE(6,603)
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1  NBEG,NEND,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If a level higher than NF has already been found and the level found
c   not the correct level, ignore all errors and continue with the 
c   binary search.
c
      IF((KV.LT.NF).AND.(KV.GE.0).AND.((IVCOR.GT.0).OR.(IEG.GT.3))) THEN
        KV = NF + 1
        EO = VMAX + 1.0d0
        IF (IVCOR.GT.0) THEN
          LLOE = (LHIE + LLOE) / 2.0d0
        ELSE
          LLOE = GV(NF-1)
          LHIE = VMAX
          IVCOR = 1
        END IF
        TRIALE = 0
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF ((KV.LT.0).AND.(EO.GT.VMAX).AND.(IVCOR.LT.20)) THEN
        KV = NF + 1
        IF (IVCOR.EQ.0) THEN
          LLOE = GV(NF-1)
          LHIE = VMAX
          IVCOR = 1
        END IF
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.
c   There are three possible situations to consider:
c     EO > VMAX : Trial energy greater than potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the case when the next trial energy is higher than the potential
c   maximum, try one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
c
      IF ((KV.LT.0).AND.(EO.GT.VMAX)) THEN
        IF (LTRY.LT.1) THEN
          LTRY = 1
          KV = 999
          EO = VMAX - 1.0d-2
c
c** If unsuccessful, then print out a warning and exit.
c
        ELSE
          IF (AWO.NE.0) THEN
            WRITE(6,605)
            WRITE(6,606) NF, EO, VMAX
          END IF
          AFLAG = -1
          GOTO 200
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find the first vibrational level (v=0), then double the
c   zero point energy estimation (AO).
c
c                         E(v=0) = VMIN + IQT*AO
c
      ELSEIF ((KV.LT.0).AND.(NF.EQ.0)) THEN
        IF (IQT.GT.1) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,611)
            WRITE(6,620) IQT, EO
          END IF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -2
            GOTO 200
          END IF
        END IF
        IQT = IQT + 1
        EO = VMIN + DBLE(IQT)*AO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If trying to find other vibrational levels (v > 0) then switch to
c   use of differences for estimating spacing.
c
      ELSEIF ((KV.LT.0).AND.(NF.GT.0)) THEN
        IF (IVDIF.GT.0) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,621) NF,IVDIF
          END IF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -3
            GOTO 200
          END IF
        END IF
        IVDIF = 1
        IF (INNER.EQ.0) THEN
          CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
        ELSE
          CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If first level found isn't v=0, try up to 3 times to 'harmonically'
c   estimate improved trial ground state energy.
c
c           E(v=0) = E(v=KV) - (E(v=KV) - VMIN)/(1 + KV/2)
c
      ELSEIF ((KV.GT.0).AND.(NF.EQ.0)) THEN
        IF (IZPE.GT.3) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,611)
            WRITE(6,622) IZPE,GV(0),KV,EO
          END IF
c
c** If differences fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -4
            GOTO 200
          END IF
        END IF
        IZPE = IZPE + 1
        EO = EO - (EO-VMIN)/(Z1+ZH/KV)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the next three cases, KV >= 0 and NF > 0.
c** If the calculated vibrational level is less than the next expected
c   level, then the estimated trial energy is too low.
c** Perhaps the difference in energy between vibrational levels v and
c   v-1 is much greater than the energy between levels v-1 and v-2. 
c
c                 E(v) - E(v-1) >> E(v-1) - E(v-2)
c
c   In which case (most likely a potential with a shelf), try twice to
c   estimate a higher trial energy.
c
c   E(v) = E(v-1) + (1+IEG/2) * (2*(E(v-1)-E(v-2)) - (E(v-2)-E(v-3)))
c
      ELSEIF (KV.LT.NF) THEN
        IF (IEG.GT.3) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,623) NF, KV
          END IF
c
c** If this fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   zero point level.
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -5
            GOTO 200
          END IF
        END IF
        IEG = IEG + 1
c
c** If a second minimum is present, then the next vibrational level may
c   be in the inner well. If so, use the inner well vibrational levels
c   to estimate the next trial energy.
c
        IF (INNER.EQ.0) THEN
          CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
        ELSE
          CALL DTENG(IEG,NF,NIWL,IWL,NVIBMX,VMIN,GV,EO)
        END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is the next expected level, then
c   continue.
c
      ELSEIF (KV.EQ.NF) THEN
        GV(NF) = EO
        INNR(NF) = INNER
        LTRY = 0
        TRIALE = 1
        LLOE = EO
        LHIE = VMAX
        CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c-----------------------------------------------------------------------
c
c** To ease confusion when using a potential with a second minimum, keep
c   track of levels that are in the outer well seperate from the levels
c   in the inner well.
c** First, calculate the rotational constant for this level (Bv).
c
   20   BV(NF) = ZH*((SWF(NBEG)/(RMIN+DBLE(NBEG-1)*RH))**2
     1           + (SWF(NEND)/(RMIN+DBLE(NEND-1)*RH))**2)
        DO I = NBEG+1,NEND-1
          BV(NF) = BV(NF) + (SWF(I)/(RMIN+DBLE(I-1)*RH))**2
        END DO
        BV(NF) = BV(NF)*RH/BZ
c
c** Double check that the calculated level is in fact located in the
c   correct well. This can be done (for v <> 0) by comparing the Bv
c   value of the new level and with the Bv values in each well. If the
c   difference is greater than 1.5 times the difference in the other
c   well, then the calculated level is probably in the wrong well.
c
        IF (NOWL.GT.0) THEN
          BVDOUT = DABS(BV(NF) - BV(OWL(NOWL-1)))
        ELSE
          BVDOUT = 9999.9d0
        END IF
        IF (NIWL.GT.0) THEN
          BVDIN = DABS(BV(NF) - BV(IWL(NIWL-1)))
        ELSE
          BVDIN = 9999.9d0
        END IF
        IF (INNER.EQ.0) THEN
          IF ((NOWL.GT.0).AND.(BVDOUT.GT.(1.5d0*BVDIN))) THEN
            IF (MCO.GE.1) THEN
              WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,BV(NF)
            END IF
            IWL(NIWL) = NF
            NIWL = NIWL + 1
          ELSE
            IF (MCO.GE.1) THEN
              WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,BV(NF)
            END IF
            OWL(NOWL) = NF
            NOWL = NOWL + 1
          END IF
        ELSE
          IF ((NIWL.GT.0).AND.(BVDIN.GT.(1.5d0*BVDOUT))) THEN
            IF (MCO.GE.1) THEN
              WRITE(21,2113) NF,'Outer',NOWL,GV(NF)-VMIN,BV(NF)
            END IF
            OWL(NOWL) = NF
            NOWL = NOWL + 1
          ELSE
            IF (MCO.GE.1) THEN
              WRITE(21,2113) NF,'Inner',NIWL,GV(NF)-VMIN,BV(NF)
            END IF
            IWL(NIWL) = NF
            NIWL = NIWL + 1
          END IF
          INNER = 0
        END IF
c
c** Look for the next uncalculated level.
c
        NF = NF + 1
        IF (INNR(NF).GE.0) THEN
c   next line commented out by JMH after discussion with RJL 24/7/02
c         NF = NF + 1
          IF (NF.LE.KVMAX) GOTO 20
        END IF
c-----------------------------------------------------------------------
c** Now estimate trial energy for next higher vibrational energy level
c   by using the Near-Dissociation Theory result that:
c
c                  (binding energy)**((NCN-2)/(2*NCN))
c
c   is (at least locally) linear in vibrational quantum number.
c
        IF (NF.EQ.1) THEN
          VDMV = ZH/(((VMAX-VMIN)/(VMAX-GV(0)))**PWI - Z1)
        ELSE
          VDMV = Z1/(((VMAX-GV(NF-2))/(VMAX-GV(NF-1)))**PWI - Z1)
        END IF
c
c** If unable to calculate the next trial energy, see if all of the
c   desired levels have been calculated. If not then turn on the warning
c   flag and quit, otherwise print out success message and quit.
c
        IF ((VDMV.LT.Z1).AND.(NCN.GT.2)) THEN
          IF (NF.LE.KVMAX) THEN
            AFLAG = -1
            WRITE(6,640) JROT, NF - 1 + VDMV
          ELSEIF (AWO.GT.0) THEN
            WRITE(6,630) KVMAX
          END IF
          GOTO 200
        END IF
c
c** Now calculate the next trial energy.
c
        EO = VMAX - (VMAX-GV(NF-1))*(Z1-Z1/VDMV)**PW
c
c** However, if the level is above the dissociation limit (for
c   potentials with barriers) then use differences to calculate the
c   next trial energy.
c
        IF (EO.GT.VMAX) CALL DTENG(IEG,NF,NOWL,OWL,NVIBMX,VMIN,GV,EO)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If the calculated vibrational level is higher then the next expected
c   level, then try once to interpolate harmonically for the missed
c   level.
c
c                   E(v) = E(v-1) + (E(KV) - E(v-1)) / 2
c
      ELSEIF (KV.GT.NF) THEN
c
c** Record vibrational level (if haven't already) for posterity.
c
        IF ((INNR(KV).EQ.-1).AND.(EO.LT.VMAX)) THEN
          GV(KV) = EO
          INNR(KV) = INNER
        END IF
        IF (IVCOR.GT.19) THEN
          IF (AWO.NE.0) THEN
            WRITE(6,610)
            WRITE(6,612)
            WRITE(6,624) IVCOR,KV,EO,(NF-1),GV(NF-1)
          END IF
c
c** If interpolation fails, then try changing the wavefunction matching
c   condition (INNER) to see if a possible second minimum contains the
c   missing level. 
c
          IF (INNER.EQ.0) THEN
            INNER = 1
            LLOE = GV(NF-1)
            LHIE = VMAX
            CALL INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c
c** If both attempts fail, then print out warning message and exit the
c   subroutine.
c
          ELSE
            AFLAG = -6
            GOTO 200
          END IF
        END IF
c
c** Use NDE theory to determine the missing level.
c
        IF (IVCOR.EQ.0) THEN
          VDU = (VMAX-EO)**PWI
          VDL = (VMAX-GV(NF-1))**PWI
          EO = VMAX - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
*        IF ((NPMIN.EQ.1).AND.(IVCOR.EQ.0)) THEN
*          VDU = (VPMAX(1)-EO)**PWI
*          VDL = (VPMAX(1)-GV(OWL(NOWL-1)))**PWI
*          EO = VPMAX(1) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
*        ELSEIF ((((INNER.EQ.1).AND.(NIWL.GT.0)).OR.
*     1           ((INNER.EQ.0).AND.(NOWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
*          VDU = (VPMAX(1)-EO)**PWI
*          VDL = (VPMAX(1)-GV(IWL(NIWL-1)))**PWI
*          EO = VPMAX(1) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
*        ELSEIF ((((INNER.EQ.0).AND.(NOWL.GT.0)).OR.
*     1           ((INNER.EQ.1).AND.(NIWL.EQ.0))).AND.(IVCOR.EQ.0)) THEN
*          VDU = (VPMAX(2)-EO)**PWI
*          VDL = (VPMAX(2)-GV(OWL(NOWL-1)))**PWI
*          EO = VPMAX(2) - (VDL + (VDU - VDL) / DBLE(KV - NF + 1))**PW
c
c** Otherwise use a binary search.
c
        ELSE
          IF ((TRIALE.EQ.1).AND.(EO.LT.LHIE)) THEN
            LHIE = EO
            IVCOR = 0
          ELSEIF (TRIALE.EQ.1) THEN
            LHIE = (LHIE + LLOE) / 2.0d0
          END IF
          TRIALE = 1
          EO = (LHIE + LLOE) / 2.0d0
        END IF
        IVCOR = IVCOR + 1
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If an unknown case occurs (quite impossible but don't quote me on
c   it) then write out an error message and exit.
c
      ELSE
        IF (AWO.NE.0) THEN
          WRITE(6,610)
          WRITE(6,666) KV,NF
        END IF
        AFLAG = -7
        GOTO 200
      END IF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set KV to the next uncalculated vibrational level to be found unless
c   looking for the highest vibrational level.
c
      IF (KV.NE.999) KV = NF
c
c** If still haven't found all of the vibrational levels then
c   look for the next vibrational level.
c
      IF ((KV.LE.KVMAX).OR.(KV.EQ.999)) GOTO 10
c
c** Otherwise, print out a message saying that all is well.
c
      IF ((KV.GT.KVMAX).AND.(AWO.GT.0)) WRITE(6,630) KVMAX
c
c** If the potential has levels in a second minimum, then print out a
c   list of those levels to channel-21 if desired.
c
      IF ((NIWL.GT.0).AND.(NOWL.GT.0)) THEN
        IF (MCO.GE.1) WRITE(21,2114) NIWL, NOWL
        IF (AWO.NE.0) THEN
          WRITE(6,605)
          WRITE(6,607)
        END IF
        AFLAG = -8
      END IF
c
c** If an error has occured, then set KVMAX to the quantum number of the
c   highest vibrational level found and print out that quantum number
c   and the energy of that level.
c
  200 IF (AFLAG.LT.0) THEN
        KVMAX = NF - 1
        IF (AWO.NE.0) WRITE(6,626) KVMAX, GV(KVMAX)
      END IF
      IF (MCO.GE.1) THEN
*        WRITE(21,2100)
*        DO I = 0,NF-1
*          IF (INNR(I).EQ.0) WRITE(21,2113) I,'Outer',I,GV(I)-VMIN,BV(I)
*          IF (INNR(I).EQ.1) WRITE(21,2113) I,'Inner',I,GV(I)-VMIN,BV(I)
*        END DO
        WRITE(21,2100)
      END IF
      RETURN
c-----------------------------------------------------------------------
  601 FORMAT(/' Solve by matching inward and outward solutions at')
  602 FORMAT(' the outermost wave function maximum, S(max), where  R = R
     1R(M)')
  603 FORMAT(' the innermost turning point  R1 = R(M)')
  605 FORMAT(/'  *** ALF WARNING ***')
  606 FORMAT(4X,'Next estimated trial energy  E(v=',I3,') =',G15.8/,4X,
     1'is greater than the potential maximum  VMAX =',G15.8)
  607 FORMAT(4X,'Potential found to have a second minimum.')
  610 FORMAT(/'  *** ALF ERROR ***')
  611 FORMAT(4X,'Attempt to find zero point level fails!')
  612 FORMAT(4X,'Attempt to find next higher vibrational level fails!')
  613 FORMAT(4X,'Number of vib levels requested=',i4,' exceeds internal 
     1ALF array dimension  NVIBMX=',i4)
  614 FORMAT(4X,'Unable to find a potential minimum.')
  615 FORMAT(4X,'There are',I3,'  potential ',A6,' in this potential. St
     1op searching after 10.')
  616 FORMAT(4X,'The potential turns over in the short range region at R
     1 = ',G15.8)
  617 FORMAT(4X,'VMAX =',G15.8,' found to be less than VMIN =',G15.8)
  620 FORMAT(4X,'Use of energy ',I1,'0% up the potential well (E =',
     1G15.8,')'/4X,' fails to produce a viable vibrational eigenstate.')
  621 FORMAT(4X,'Use of differences to estimate the energy for the next'
     1/4X,' vibrational level (v=',I3,') failed after',I3,'  attempt.')
  622 FORMAT(4X,'After',I3,' tries to harmonically estimate the zero-poi
     1nt energy,'/4X,' initial trial energy',G15.8,'   had yielded   E(v
     2=',I3,') =',G15.8)
  623 FORMAT(4X,'Expecting to find level (v=',I3,') but found level (v='
     1,I3,')')
  624 FORMAT(4X,'After',I3,' tries, failed to interpolate trial energy b
     1etween'/4X,'E(v=',I3,') =',G15.8,'   and   E(v=',I3,') =',G15.8)
  626 FORMAT(4X,'The highest calculated level is  E(v=',I3,') =',G15.8)
  630 FORMAT(/' ALF successfully finds all vibrational levels up to KVMA
     1X=',I3)
  640 FORMAT(/' ALF finds all  J=',i3,'  vib. levels below  vD=',F7.3,
     1  '  estimated by N-D theory')
  650 FORMAT(/' There were',I3,'  potential minima found with the absolu
     1te minimum'/4X,'VMIN =',G15.8,'  cm-1.')
  651 FORMAT(/' There were',I3,'  potential maxima found with the absolu
     1te maximum'/4X,'VMAX =',G15.8,'  cm-1.')
  666 FORMAT(4X,'Undefined case for automatic search.'/,4X,'Values of KV
     1 =',I3,'  and NF =',I3)
 2100 FORMAT(/1X,39('=='))
 2101 FORMAT(/1X,39('--'))
 2110 FORMAT(/' Limits and increment of integration (in Angstroms):'
     1 /'    RMIN =',F6.3,'    RMAX =',F7.3,'    RH =',F9.6,
     2 //' Generate    BZ =',G19.12,' ((1/cm-1)(1/Angstroms**2))'
     3 /' from ZMU:',F15.11,' (amu)')
 2111 FORMAT(/' Eigenvalue convergence criterion is   EPS =',G11.4,'(cm-
     11)')
 2112 FORMAT(/' Calculating properties of the potential described above.
     1 '/' Use Airy function at 3-rd turning point as outer boundary'
     2 /' condition for quasibound levels.')
 2113 FORMAT(' v=',I3,4X,'v(',A5,')=',I3,4X,'Gv=',F16.9,4X,'Bv=',F16.12)
 2114 FORMAT(/' Found',I4,' level(s) in the inner well and',I4,' level(s
     1) in the outer well.')
      END
c***********************************************************************
      SUBROUTINE INITVAL(IQT,IVDIF,IZPE,IEG,IVCOR)
c***********************************************************************
c** This subroutine reinitializes the condition flags when considering a
c   new case (found next vibrational level or finding level in inner
c   well - INNER = 1).
c
c** On entry and exit:
c    IQT     Case when KV < 0 and NF = 0
c            determines the value used for the initial trial energy.
c    IVDIF   Case when KV < 0 and NF > 0
c            is the flag denoting the use of differences to calculate
c            trial energies.
c    IZPE    Case when KV > 0 and NF = 0
c            is the number of times the zero point energy (v = 0) has
c            been estimated harmonically.
c    IEG     Case when KV < NF and NF > 0
c            are the number of times that a larger trial energy is used
c            to find the next level.
c    IVCOR   Case when KV > NF and NF > 0
c            are the number of times that a smaller trial energy is used
c            to find the next level.
c
      INTEGER IZPE,IVDIF,IVCOR,IQT,IEG
c
      IQT = 1
      IVDIF = 0
      IZPE = 0
      IEG = 0
      IVCOR = 0
      RETURN
      END
c***********************************************************************
      SUBROUTINE DTENG(IEG,NF,NVEL,VEL,NVIBMX,VMIN,GV,EO)
c***********************************************************************
c** This subroutine calculates the next trial energy using differences.
c
c** On entry:
c    IEG     factor by which a larger trial energy should be calculated:
c              NVEL = 2 : Increase correction by increments of 25%
c              NVEL > 2 : Increase correction by increments of 50%
c    NF      is the highest calculated vibrational level.
c    NVEL    is the number of levels found in the potential well.
c    VEL(v)  keeps track of all levels in the potential well.
c    NVIBMX  is the maximum number of vibrational levels (dimension).
c    VMIN    is the absolute value of the potential minimum (cm-1).
c    GV(v)   contains the vibrational energy level spacings
c            and rotational constants for each level (cm-1).
c
c** On exit:
c    EO      is the calculated trial energy.
c
      INTEGER IEG,NF,NVEL,NVIBMX,VEL(0:NVIBMX)
c
      REAL*8 VMIN,GV(0:NVIBMX),EO,ZQ,Z1,Z2
c
      DATA ZQ/0.25D0/,Z1/1.D0/,Z2/2.D0/
c
c** If determining the first (non-zero point energy) level in the well,
c   then use the last determined level in the other well plus a larger
c   than harmonic correction that becomes smaller with each new
c   itteration.
c
c           E(v=0) = E(v=NF-1) + (E(v=NF-1)-VMIN)/(NF-1+IEG/4) 
c
      IF (NVEL.EQ.0) THEN
        EO = GV(NF-1) + (GV(NF-1) - VMIN)/(NF - 1 + ZQ*DBLE(IEG))
c
c** Try to get v = 1 using smaller-than-harmonic spacing.
c
c              E(v=1) = E(v=0) + (1.0+IEG/4)*(E(v=0)-VMIN)
c
      ELSEIF (NVEL.EQ.1) THEN
        EO = GV(VEL(0)) + (Z1+DBLE(IEG)*ZQ)*(GV(VEL(0))-VMIN)
c
c** Try to get v = 2 using a sequentially increasing correction.
c
c             E(v=2) = E(v=1) + (0.8+IEG/4)*(E(v=1)-E(v=0))
c
      ELSEIF (NVEL.EQ.2) THEN
        EO = GV(VEL(1)) + (0.8d0+DBLE(IEG)*ZQ)*(GV(VEL(1))-GV(VEL(0)))
c
c** Try to get v > 2 using a sequentially increasing correction.
c
c        E(v) = E(v-1) + (1.0+IEG)*(2.0*E(v-1)-3.0*E(v-2)+E(v-3))
c
      ELSE
        EO = GV(VEL(NVEL-1)) + (Z1+DBLE(IEG))
     1       *(Z2*GV(VEL(NVEL-1))-3.0d0*GV(VEL(NVEL-2))+GV(VEL(NVEL-3)))
      END IF
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
