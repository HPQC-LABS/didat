c***********************************************************************
      SUBROUTINE ALF(NDP,RMIN,RH,V,SWF,VLIM,KVMAX,AFLAG,ZMU,EPS,GV,
     1                                                INNODE,INNR,IWR)
c***********************************************************************
c   Version 2.0 dated  September 21, 2008. 
c-----------------------------------------------------------------------
c** The subroutine ALF (Automatic vibrational Level Finder) will
c   automatically generate the eigenvalues from the first vibrational
c   level (v=0) to a user specified level (v=KVMAX) or the highest
c   allowed vibrational level of a given smooth single (or double)
c   minimum potential (V). These energies are stored and returned to the
c   calling program in the molecular constants array GV(v=0-KVMAX).
c** For any errors that cannot be resolved within the subroutine, ALF
c   returns AFLAG with a value that defines which error had occured.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++   COPYRIGHT 2008  by  Robert J. Le Roy   ++++++++++++++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c     of it without the express written permission of the authors.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++++++ Please inform me of any bugs, by phone at: (519)888-4051 +++++++
c+++++++++ by e-mail to: leroy@uwaterloo.ca , or by Post at: +++++++++++
c+++ Dept. of Chemistry, Univ. Waterloo, Waterloo, Ontario  N2L 3G1 ++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
c    KVMAX  is v for the highest vibrational level we wish to find.
c    AFLAG  is rot.quantum J for the (centrifugally distorted) potential
c    ZMU    is the reduced mass of the diatom (amu).
c    EPS    is the energy convergence criterion (cm-1).
c    INNODE specifies whether wave fx. initiation @ RMIN starts with a
c        note (normal case: INNODE > 0) or zero slope (when INNODE.le.0)
c    IWR    specifies the level of printing inside SCHRQ
c           <> 0 : print error & warning descriptions.
c           >= 1 : also print final eigenvalues & node count.
c           >= 2 : also show end-of-range wave function amplitudes.
c           >= 3 : print also intermediate trial eigenvalues, etc.
c
c** On exit:
c    KVMAX   is vib.quantum number for the highest vibrational level
c            found (may be less than the input value of KVMAX).
c    AFLAG   returns calculation outcome to calling program.
c            >=  0 : found all levels to v=KVMAX{input} & AFLAG= J 
c             = -1 : KVMAX larger than number of levels found.
c    GV(v)   contains the vibrational energy levels found for v=0-KVMAX
c    INNR(v) labels each level as belonging to the inner (INNR = 1) or
c            outer (INNR = 0) well.
c
c** Flags: Modify only when debugging.
c    AWO   specifies the level of printing inside ALF
c          <> 0 : print error & warning descriptions.
c          >  0 : also print intermediate ALF messages.
c    INNER specifies wave function matching (& initiation) conditions.
c        .le.0 : Match inward & outward solutions at outermost well t.p.
c          > 0 : Match at innermost well inner turning point
c        For most normal cases set INNER = 0,  but ......
c            To find "inner-well-dominated" solutions of an asymmetric
c            double minimum potential, set  INNER > 0.
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
      PARAMETER (NVIBMX= 400)
c
c** NF counts levels found in automatic search option
c
      INTEGER NDP,KVMAX,KV,KVB,KVBB,AFLAG,NF,NBEG,NEND,NBEGG(0:NVIBMX),
     1  NENDD(0:NVIBMX),INNR(0:NVIBMX),ICOR,IWR,
     2  I,LTRY,AWO,INNODE,INNER,
     3  LPRWF,JROT,NPMIN, NPMAX, IPMIN(10),IPMINN
c
      REAL*8 RMIN,RMAX,RH,V(NDP),SWF(NDP),VLIM,EO,ZMU,EPS,
     1  BZ,BFCT,GAMA,VMIN,VMAX,RE,PMAX, ESAV, ZPEHO, DGDV2, BMAX,
     2  GV(0:KVMAX),VPMIN(10),RPMIN(10),VPMAX(10),RPMAX(10)
c
      DATA AWO/1/,LPRWF/0/,KVB/-1/,KVBB/-2/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Check that the array dimensions are adequate.
      IF(KVMAX.GT.NVIBMX) THEN
          WRITE(6,602) KVMAX, NVIBMX
          STOP
          ENDIF
c
c** Initialize the remaining variables and flags.
      NF= 0                                ! NF is label of level being sought
      KV= 0
      INNER= 0
      LTRY= 0
c** Initialize level counters for each well.
      DO  I= 0,KVMAX
          INNR(I)= -1
          ENDDO
c** Store input rotational quantum number.
      JROT= AFLAG
      AFLAG= -1
c
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920d0
      BFCT= BZ*RH*RH
c
c** RMAX is the outermost radial distance at which potential is defined
      RMAX= RMIN + DBLE(NDP-1)*RH
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Locate the potential minima.
      NPMIN= 0
      DO  I= 2,NDP-1
          IF((V(I).LT.V(I-1)).AND.(V(I).LT.V(I+1))) THEN
              NPMIN= NPMIN + 1
              IPMIN(NPMIN)= I
              RPMIN(NPMIN)= RMIN + DBLE(I-1)*RH
              VPMIN(NPMIN)= V(I) / BFCT
              IF(NPMIN.EQ.10) GOTO 10
              ENDIF
          END DO
   10 IF(NPMIN.EQ.0) THEN
          IF(V(2).LE.V(1)) THEN
c** If NO minimum can be found, then print a warning and stop.
              WRITE(6,604) JROT
              KVMAX= -1
              RETURN
              ENDIF
c...  but if potl. alway has positive slope, mesh point 1 is minimum
          NPMIN= 1
          IPMIN(NPMIN)= 1
          VPMIN(NPMIN)= V(1)/BFCT
          RPMIN(NPMIN)= RMIN
          WRITE(6,606) VPMIN(1),RMIN
          ENDIF
c** If more than two minima are found, then print a warning and hope ...
      IF(NPMIN.GT.2) THEN
          WRITE(6,608) NPMIN, 'minima'
          ENDIF
c** Now find the absolute potential minimum.
      VMIN= VPMIN(1)
      IPMINN= IPMIN(1)
      RE= RPMIN(1)
      IF(NPMIN.GE.2) THEN
          DO  I= 2,NPMIN
              IF(VMIN.GT.VPMIN(I)) THEN
                  IPMINN= IPMIN(I)
                  VMIN= VPMIN(I)
                  RE= RPMIN(I)
                  ENDIF
              END DO 
          ENDIF
c
c** Locate any potential maxima (if they exists).
      NPMAX= 0
      DO  I= 2,NDP-1
          IF((V(I).GT.V(I-1)).AND.(V(I).GT.V(I+1))) THEN
              NPMAX= NPMAX + 1
              RPMAX(NPMAX)= RMIN + DBLE(I-1)*RH
              VPMAX(NPMAX)= V(I) / BFCT
              IF(NPMAX.EQ.10) GOTO 20
              ENDIF
          END DO
c** If more than three maxima found, then print a warning and hope ...
   20 IF(NPMAX.GT.3) THEN
          WRITE(6,608) NPMAX, 'maxima'
          ENDIF
      IF(NPMAX.EQ.0) THEN
c** If no maxima were found, then set energy maximum to be the energy 
c   of the potential asymptote  VLIM
          NPMAX= 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= VLIM
          ENDIF
c** If there is no barrier past outermost minimum, then set a final
c   VPMAX to be  VLIM 
      IF(RPMAX(NPMAX).LT.RPMIN(NPMIN)) THEN
          NPMAX= NPMAX + 1
          RPMAX(NPMAX)= RMAX
          VPMAX(NPMAX)= VLIM
          ENDIF
c
c** If innermost maximum lies inside innermost minimum, the potential 
c   turns over in short range region and should not be used, so PRINT  
c   a Warning and STOP.
      IF(RPMAX(1).LT.RPMIN(1)) THEN
          WRITE(6,610) RPMAX(1)
          STOP
          ENDIF
c
c** Now find the highest potential maximum.
      VMAX= VPMAX(1)
      IF(NPMAX.GT.1) THEN
          DO  I= 2,NPMAX
              IF(VMAX.LT.VPMAX(I)) VMAX= VPMAX(I)
              END DO 
          ENDIF
c
c** If the absolute potential maximum is lower than the absolute
c   potential minimum, then print out an error statement and quit.
      IF(VMAX.LE.VMIN) THEN
          WRITE(6,612)
          STOP
          ENDIF
c
c** Otherwise, print out potential extrema count, if desired.
      IF(NPMIN.GT.1) THEN
          WRITE(6,614) NPMIN, (VPMIN(I),I= 1,NPMIN)
          WRITE(6,616) (RPMIN(I), I= 1,NPMIN)
          WRITE(6,618) NPMAX, (VPMAX(I),I= 1,NPMAX)
          WRITE(6,616) (RPMAX(I), I= 1,NPMAX)
          IF(NPMIN.GT.2) THEN
              WRITE(6,620)
              STOP
              ENDIF
          ENDIF
c
c** Set BMAX as barrier height of double-minimum potential
      BMAX= -9.d+09
      IF(NPMIN.GT.1) THEN
          DO  I= 1,NPMAX
              IF((RPMAX(I).GT.RPMIN(1)).AND.(RPMAX(I).LT.RPMIN(2)))
     1            BMAX= VPMAX(I)
              ENDDO
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c*** Use harmonic approximation to estimate zero point energy.
      ZPEHO= DSQRT((V(IPMINN+20)-V(IPMINN))/400.d0)/BFCT
      EO= VMIN + ZPEHO
c
c=========== Begin Actual Eigenvalue Calculation Loop Here =============
c** Compute eigenvalues ... etc. up to the KVMAX'th vibrational level.
c** When attempts to find the next eigenvalue fails, then perhaps the
c   next level is located in a second (inner) well. If so, then the
c   subroutine will set INNER = 1, and attempt to find that level.
c
      ICOR= 0
  100 KVBB= KVB
      KVB= KV
      KV= NF
  110 ESAV= EO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Call subroutine SCHRQ to find eigenvalue EO and eigenfunction SWF(I).
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL SCHRQ(KV,JROT,EO,GAMA,PMAX,VLIM,V,SWF,BFCT,EPS,RMIN,RH,NDP,
     1                               NBEG,NEND,INNODE,INNER,IWR,LPRWF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.LT.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** The SCHRQ error condition is KV < 0.  Allow for 3 cases:
c     EO > VMAX : energy from previous trial above potential maximum
c     NF = 0 : Looking for the first vibrational level (v = 0)
c     NF > 0 : Looking for the other vibrational levels (v > 0)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(EO.GT.VMAX) THEN
c** For the case when the previous trial gave energy above the potential
c   maximum, make one last ditch attempt to find the highest bound level
c   (quasi or otherwise) in the potential.
              IF(LTRY.LT.1) THEN
                  LTRY= 1
                  KV= 999
                  EO= VMAX - 0.01d0
                  GOTO 110
c... if that was unsuccessful, then print out a warning and exit.
                ELSE
                  WRITE(6,622) NF, EO, VMAX
                  KV= NF-1
                  GOTO 200
                ENDIF
              ENDIF
          WRITE(6,624) NF,JROT,ESAV
c.. eigenvalue of -9.9d9 signifies that eigenvalue search failed completely
          KVMAX= NF-1
          EO= -9.9d9
          RETURN
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If calculated vibrational level is the desired level, NF, then ...
c   call SCECOR to calculate dG/dv and predict next higher level
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.EQ.NF) THEN
      NBEGG(KV)= NBEG
      NENDD(KV)= NEND
          GV(NF)= EO
          INNR(NF)= INNER
  120     NF= NF + 1
          IF(NF.LE.KVMAX) THEN
              IF(INNR(NF).GT.0) GOTO 120
c... if the next level was found earlier in overshoot ... 
            ELSE
              IF(AWO.GT.0) WRITE(6,626) JROT,KVMAX
              AFLAG= JROT
              RETURN
            ENDIF
          ICOR= 0
          CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,V,
     1                                                BMAX,VLIM,DGDV2)
          IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value below it
              EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
              ICOR= 20
            ELSE
              IF(DGDV2.LT.0.d0) THEN
                  WRITE(6,628) JROT,EO
                  AFLAG= -1
                  GOTO 200
                  ENDIF
            ENDIF
          LTRY= 0
          KV= NF
          GOTO 100
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(KV.NE.NF) THEN
c*** If last level found is not the desired one ...
          IF(INNR(KV).EQ.-1) THEN
c... Record vibrational level (if haven't already) for posterity.
              GV(KV)= EO
              INNR(KV)= INNER
              ENDIF
          ICOR= ICOR+1
          IF(ICOR.LE.20) THEN
c... Call subroutine using semiclassical methods to estimate correct energy
              CALL SCECOR(KV,NF,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,V,
     1                                                BMAX,VLIM,DGDV2)
              IF(EO.GT.VPMAX(NPMAX)) THEN
c... if estimated energy above highest barrier, set value below it
                  KV= 999
                  EO=  VPMAX(NPMAX) - 0.05d0*DGDV2
                  ENDIF
              GOTO 100
              ENDIF
c** If the calculated wavefunction is still for the wrong vibrational
c   level, then write out a warning return
          WRITE(6,630) NF,JROT
          KVMAX= NF-1
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  200 IF(AFLAG.LT.0) THEN
c** If unable to find all KVMAX+1 levels requested, then return KVMAX as
c  v for the highest vibrational level actually found, and print out the
c  the energy of that level.
          KVMAX= KV
          IF(AWO.NE.0) WRITE(6,632) KVMAX, GV(KVMAX)
          ENDIF
      RETURN
c-----------------------------------------------------------------------
  602 FORMAT(/'  *** ALF ERROR ***'/4X,'Number of vib levels requested='
     1 ,i4,' exceeds internal ALF array dimension  NVIBMX=',i4)
  604 FORMAT(/' *** ALF ERROR ***   Find NO potential minima for   J=',
     1  i4)
  606 FORMAT(/'  ALF  finds onee potential minimum of',1PD15.7,
     1  '  at  R(1)=',0Pf9.6)
  608 FORMAT(/'  *** ALF WARNING ***'/4X,'There are',I3,'  potential ',
     1  A6,' in this potential. Stop searching after 10.')
  610 FORMAT(/'  *** ALF ERROR ***'/ 4X,'The potential turns over in the
     1 short range region at  R= ',G15.8)
  612 FORMAT(/'  *** ALF ERROR ***'/4X,'VMAX =',G15.8,' found to be less
     1 than VMIN =',G15.8)
  614 FORMAT(' Find',I3,'  potential minima:   Vmin=',8F11.3)
  616 FORMAT(19x,'located at   R =',8f11.5)
  618 FORMAT(' Find',I3,'  potential maxima:   Vmax=',8F11.3)
  620 FORMAT(' *** So  STOP !!!!')
  622 FORMAT(/' ALF search finds next estimated trial energy  E(v=',I3,
     1 ')=',G15.8/8X,'lies above potential maximum or asymptote at  VMAX
     2=',G15.8)
  624 FORMAT(/' *** SCHRQ FAILS in ALF when searching for  v=',i3,
     1  ' J=',i3,'   with   EO=',f9.3/5x,'Check range and/or contact R.J
     2. Le Roy [leroy@uwaterloo.ca]')
  626 FORMAT(/' ALF successfully finds all (J=',i3,') vibrational levels
     1 up to   v= KVMAX=',I3)
  628 FORMAT(/' *** ERROR:   at   E(J=',i3,')=',f10.3,"  SCECOR  finds n
     1o outer' turning point")
  630 FORMAT(4x,'ALF fails to find level   v=',i3,', J=',i3)
  632 FORMAT(' ALF finds the highest calculated level is  E(v=',I3,
     1  ')=',G15.8 /)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE SCECOR(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,V,
     1                                                BMAX,VLIM,DGDV2)
c** Subroutine calculates (approximate!) semiclassical estimate of 
c  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential 
c  {V(i),i=1,NDP} (in 'internal units' {cm-1 * BFCT}), and uses those
c  results to estimate energy of level KVLEV
c** If the 'clever' semiclassical procedure fails - try a brute force
c  step-by-step search, using alternately INNER & OUTER well starting
c** BMAX is internal barrier maximum energy for double-well case, 
c   and very large negative number for single-well potential
c** On return, negative DGDV2 signals error!  No phase integrals found
c
      INTEGER I,II,I1,I2,I3,I4,IV1,IV2,INNER,ICOR,JROT,KV,KVB,KVLEV,
     1  KVDIF,NDP,IDIF,BRUTE,IB,IWR,NPMAX
      REAL*8 EO,DE0,RH,BFCT,ARG2,ARG3,EINT,VPH1,VPH2,DGDV1,DGDV2,DGDVM,
     1  DGDVB,EBRUTE,DEBRUTE,DE1,DE2,Y1,Y2,Y3,RT,ANS1,ANS2,XDIF,VLIM,
     2  BMAX,V(NDP)
      SAVE BRUTE,EBRUTE,DEBRUTE,DGDVB
      DATA DGDVB/-1.d0/,KVB/-1/
c
      DGDV2= -1.d0
      EINT= EO*BFCT
      KVDIF= KVLEV- KV
      IF(ICOR.EQ.1) BRUTE= 0
      I3= NDP
c*** For Quasibound levels, first search inward to classically forbidden
      IF(EO.GT.VLIM) THEN
          DO  I= NDP,1,-1
              I3= I
              IF(V(I).GT.EINT) GOTO 8
              ENDDO
          ENDIF
c*** First, search inward for outermost turning point
    8 DO  I= I3,1,-1
          I4= I
          IF(V(I).LT.EINT) GOTO 10
          ENDDO
c*** If never found an 'outer' turning point (e.g., above qbdd. barier)
c  then simply return with negative  DGDV2  as error flag
      RETURN
c... Now collect vibrational phase and its energy deriv. over outer well
   10 Y1= EINT- V(I4+1)
      Y2= EINT- V(I4)
      Y3= EINT- V(I4-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH2= 0.5d0*ARG2 + ANS2/RH
      DGDV2= 0.5d0/ARG2 + ANS1/RH
      DO  I= I4-2,1,-1
c... now, collect (v+1/2) and dv/dG integrals to next turning point ...
          II= I
          IF(V(I).GT.EINT) GO TO 12
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH2= VPH2+ ARG2
          DGDV2= DGDV2+ 1.d0/ARG2
          ENDDO
   12 I3= II+1
      Y1= EINT- V(I3-1)
      Y2= EINT- V(I3)
      Y3= EINT- V(I3+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH2= (VPH2 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV2= DGDV2 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV2= 6.283185308d0/(BFCT*DGDV2)
c*** Next, search for innermost turning point
      DO  I= 1,NDP 
          I1= I
          IF(V(I).LT.EINT) GOTO 20
c... then collect vibrational phase and its energy deriv. over outer well
          ENDDO
c
   20 IF(I1.EQ.1) THEN
          WRITE(6,602) JROT,EO
          STOP
          ENDIF
      IF(I1.GE.I3) THEN
c*** For single-well potential or above barrier of double-well potential
          IF(IWR.GE.2) WRITE(6,600) ICOR,KV,JROT,EO,VPH2-0.5d0,DGDV2
          IF((KV.NE.(KVLEV-1)).AND.(DGDVB.GT.0.d0)) THEN
c... If got wrong level ........
              IF((EO-BMAX).GT.(2.d0*DGDV2)) THEN
c... 'Normal' case: use B-S plot area to estimate correct energy
                  DE0= KVDIF*(DGDV2- 0.5d0*(DGDV2-DGDVB)/DFLOAT(KV-KVB))
                  EO= EO+ DE0 
                  KVB= KV
                  KV= KVLEV
                  RETURN
                ELSE
c... but close to barrier in double-well potential, switch to 'BRUTE'
                  BRUTE=BRUTE+ 1
                  DGDV1= DGDV2
                  XDIF= SIGN(1,KVDIF)  
                  GOTO 54
                ENDIF
              ENDIF
          IF(KV.EQ.0) THEN
c** Normally:  use B-S plot considerations to estimate next level energy
c... use harmonic estimate for v=1 energy
              EO= EO+ DGDV2
            ELSE
c... estimate Delta(G) based on linear Birge-Sponer
              DE0= 0.5d0*(3.d0*DGDV2 - DGDVB)
              IF((2.d0*DGDV2).GT.DGDVB) THEN
c... if linear Birge-Sponer predicts another level, then use it
                  EO= EO+ DE0
                ELSE
c... otherwise, consider Delta{G} extrapolation for next level...
                  DE0= 0.5d0*(3.d0*DSQRT(DGDV2)- DSQRT(DGDVB))
                  IF(2.d0*DE0.GT.DSQRT(DGDV2)) THEN
c... if linear NDT(n=6) predicts another level use its predicted spacing
                      EO= EO+ DE0**2
                    ELSE   
c... but if NDT(n=6) it predicts no more levels, go back to last DGDV
                      EO= EO+ DGDV2
                    ENDIF
                ENDIF 
            ENDIF
          DGDVB= DGDV2
          KVB= KV
          INNER= 0
          RETURN
          ENDIF
c
c*** For a double-well potential, collect vibrational phase and its 
c   energy derivative over the inner well
      Y1= EINT- V(I1-1)
      Y2= EINT- V(I1)
      Y3= EINT- V(I1+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH1= 0.5d0*ARG2 + ANS2/RH
      DGDV1= 0.5d0/ARG2 + ANS1/RH
      DO  I= I1+2,NDP
c... now, collect integral and count nodes outward to next turning point ...
          IF(V(I).GT.EINT) GO TO 22
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH1= VPH1+ ARG2
          DGDV1= DGDV1+ 1.d0/ARG2
          ENDDO
   22 I2= I-1
      Y1= EINT- V(I2+1)
      Y2= EINT- V(I2)
      Y3= EINT- V(I2-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH1= (VPH1 - ARG2 - 0.5d0*ARG3 + ANS2/RH)/3.141592654d0
      DGDV1= DGDV1 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1/RH
      DGDV1= 6.28318531d0/(BFCT*DGDV1)
      DGDVM= DGDV1*DGDV2/(DGDV1+DGDV2)
      IF(KVDIF.EQ.0) THEN
c** If already at level sought, return
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
          RETURN
          ENDIF
c
c** Check whether looking for higher or lower level ...
      IDIF= SIGN(1,KVDIF)
      XDIF= IDIF
      IF((ICOR.GE.6).AND.((IABS(KVDIF).EQ.1).OR.(BRUTE.GT.0))) GOTO 50
c*** 'Conventional' semiclassical search for neared INNER or OUTER well level
      IF(INNER.LE.0) THEN
c... and current energy EO is for an outer-well level ...
          DE2= DGDV2*XDIF
          IV1= INT(VPH1+ 0.5d0)
          DE1= (DFLOAT(IV1) + 0.5d0 - VPH1)*DGDV1*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   30     IF(DABS(DE1).LT.DABS(DE2)) THEN
              INNER= 1
              EO= EO+ DE1
              DE1= DGDV1*XDIF
            ELSE
              INNER= 0
              EO= EO+ DE2
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 30
          ENDIF
      IF(INNER.GT.0) THEN
c... and current energy EO is for an inner-well level ...
          DE1= DGDV1*XDIF
          IV2= INT(VPH2+ 0.5d0)
          DE2= (DFLOAT(IV2) + 0.5d0 - VPH2)*DGDV2*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   40     IF(DABS(DE2).LT.DABS(DE1)) THEN
              INNER= 0
              EO= EO+ DE2
              DE2= DGDV2*XDIF
            ELSE
              INNER= 1
              EO= EO+ DE1
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 40
          ENDIF
   50 BRUTE= BRUTE+ 1
c*** Now .. Brute force search for desired level !
      IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   54 IF(BRUTE.EQ.1) THEN
c... in first brute-force step, use previous energy with opposite INNER
          EBRUTE= EO 
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          DEBRUTE= DMIN1(DGDV1,DGDV2)*XDIF*0.3d0
          RETURN
          ENDIF
      IB= BRUTE/2
c... in subsequent even steps, lower EO by DEBRUTE/10 for same INNER
      IF((IB+IB).EQ.BRUTE) THEN
          EBRUTE= EBRUTE+ DEBRUTE
          EO= EBRUTE
          RETURN
        ELSE
c... in subsequent odd steps, lower repeat previous EO with INNER changed
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          EO= EBRUTE
          RETURN
        ENDIF
c     RETURN
  600 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',f10.2,
     1 '  v(SC)=',F8.3,'  dGdv=',f8.3)
  602 FORMAT(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2 )
  610 FORMAT(' Double well   E(v=',i3,', J=',I3,')=',f9.3,
     1 ':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,
     2 ' (ICOR=',I2,')',8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2 )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

