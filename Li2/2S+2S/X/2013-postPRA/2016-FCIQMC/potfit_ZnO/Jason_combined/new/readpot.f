c**********************************************************************N(
      SUBROUTINE READPOT(ISTATE,SLABL,OSEL)
c**********************************************************************
c** This subroutine will generate one of four possible types of
c   analytical molecular potentials for the direct Hamiltonian fitting
c   program. These four families are:
c     - The Expanded Morse Oscillator (EMO potential
c     - The Morse/Long-Range (MLR) or Morse/Lennard-Jones (MLJ) potential
c     - The Double-Exponential Long-Range (DELR) potential
c     - The Generalized Potential Energy Function (GPEF)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                   Version of 19 August 2009
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being fitted to
c   SLABL   is the two-character label identifying that state
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOB.h'
c-----------------------------------------------------------------------
c** Type statements for input or local variables
      INTEGER I, ISTATE, IISTP, m, MMN, npow, OSEL(NSTATEMX)
      CHARACTER*2 SLABL(-3:NSTATEMX)
      REAL*8 ZMASE, RR(NPNTMX),VV(NPNTMX)
      DATA ZMASE /5.4857990945D-04/
c
c** Set some defaults for parameters not common to all models ...
      IFXDe(ISTATE)= 1
      IFXRe(ISTATE)= 1
      DO  m= 1, NCMMmax
          IFXCm(m,ISTATE)= 1
          ENDDO
c-----------------------------------------------------------------------
c** First choose potential model and select form of BOB representation
c  PSEL(s)  choses the type of analytical potential to be fitted to:
c         = -2 : represent each distinct observed level of this state as
c               an independent term value [an alternative to an 'FS' 
c               treatmentof transitions involving that state]
c         = 0 : Use fixed potential defined by LEVEL's PREPOT routine
c         = 1 : Use the Expanded Morse Oscillator (p)
c         = 2 : Use the Morse Lennard-Jones Potential.
c         = 3 : Use the Double-Exponential Long-Range Potential.
c         = 4 : Use the Generalized Potential Energy Function.
c  VLIM(s)  is the fixed absolute energy of the potential asymptote
c  BOBCN  is a flag to denote reference & scaling for BOB corrections 
c      = 0  using differences as per RJL [JMS 194,189(1999)]
c      = 1  use 'clamped nuclei' limit, m_e1/MASS scaling.
c  OSEL(s)  controls printout of radial function arrays to Ch. 10-16.
c       OSEL > 0: Export to file every OSEL'th point of final function
c=======================================================================
      READ(5,*) PSEL(ISTATE), VLIM(ISTATE), BOBCN(ISTATE), OSEL(ISTATE)
c=======================================================================
c** For term value fits ... no further reads needed!
      IF(PSEL(ISTATE).EQ. -2) THEN
          WRITE(6,604) SLABL(ISTATE)
          NUA(ISTATE)= -1
          NUB(ISTATE)= -1
          NTA(ISTATE)= -1
          NTB(ISTATE)= -1
          RETURN
          ENDIF
c-----------------------------------------------------------------------
c** Now to read in the range and mesh for the numerical integration
c  RMIN/MAX(s)  define the range over which this potential is defined.
c  RH(s)  specifies radial mesh for numerical integration for this state
c=======================================================================
      READ(5,*) RMIN(ISTATE), RMAX(ISTATE), RH(ISTATE)
c=======================================================================
      NDATPT(ISTATE)= (RMAX(ISTATE)-RMIN(ISTATE))/RH(ISTATE)+1.0001d0
      NDATPT(ISTATE)= MIN(NPNTMX-1,NDATPT(ISTATE))
      RMAX(ISTATE)= RMIN(ISTATE) + RH(ISTATE)*DBLE(NDATPT(ISTATE)-1)
      DO  I= 1, NDATPT(ISTATE)
          RR(I)= RMIN(ISTATE)+ DBLE(I-1)*RH(ISTATE)
          RD(I,ISTATE)= RR(I)
          ENDDO
      IF(PSEL(ISTATE).EQ.0) THEN
c-----------------------------------------------------------------------
c** For case of a fixed potential defined by read-in turning points, 
c  subroutine PREPOTT reads those points & generates potential array
c-----------------------------------------------------------------------
          WRITE(6,600)SLABL(ISTATE),RMIN(ISTATE),RMAX(ISTATE),RH(ISTATE)
          CALL PREPOTT(1,AN(1),AN(2),MN(1,1),MN(2,1),NDATPT(ISTATE),
     1                                             VLIM(ISTATE),RR,VV)
          DO  I= 1, NDATPT(ISTATE)
              VPOT(I,ISTATE)= VV(I)
              ENDDO
          NUA(ISTATE)= -1
          NUB(ISTATE)= -1
          NTA(ISTATE)= -1
          NTB(ISTATE)= -1
          NwCFT(ISTATE)= -1
          RETURN
          ENDIF
      IF(PSEL(ISTATE).EQ.2) THEN
c-----------------------------------------------------------------------
c** For MLJ or MLR, read number of inverse-power terms NCMM in long-range
c   form   V(R) -> VLIM - CmVAL(1)/R**MMLR(1) - CmVAL(2)/R**MMLR(2) - ...
c** If IDF=0,1,2   apply 'normalized' incomplete gamma function damping
c     function of order  m  to each Cm/r^m  term.  Otherwise, no damping
c    if IDF < 0
c,  IDF specifies damping form,
c      = 0 , Douketis-Scoles form with l = -1 (recommended) 
c      = 1 , Unmodifiend Douketis-Scoles form (l = 0)
c      = 2 , Tang-Toennies form
c   RHOab is the scaling factor for molecules other than H2 calculated 
c using  RHOab=  2*RHOa*RHOb/(RHOa+RHOb) where RHO for a molecule is equal to
c (I_p^A/I_p^H)^0.66 where I_p^A is the ionization potential of atom A 
c and Ih is the
c ionization potential of atomic hydrogen.
c** For alkali dimer A-state: Aubert-Frecon ULR(r) [PRA 55, 3458 (1997)]
c  is special case with MMLR= {3, 0, 6, 6, 8, 8, ...} 
c=======================================================================
          READ(5,*) NCMM(ISTATE), IDF(ISTATE), RHOab(ISTATE)
     1     
          DO  m= 1,NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE), IFXCm(m,ISTATE)
              ENDDO
          ENDIF
c
      IF(PSEL(ISTATE).EQ.3) THEN
c-----------------------------------------------------------------------
c** For DELR potential, read parameter to define form of `LR' part
c  Switch to Tang-Toennies-type damping function with order s.th.
c      Dm(r)*Cm/r^m --> {constant}  as  r --> 0
c=======================================================================
          READ(5,*) NCMM(ISTATE), IDF(ISTATE), RHOab(ISTATE)
     1         
          
          DO  m= 1, NCMM(ISTATE)
              READ(5,*) MMLR(m,ISTATE), CmVAL(m,ISTATE)
              ENDDO
c=======================================================================
          ENDIF
c-----------------------------------------------------------------------
      IF(PSEL(ISTATE).EQ.4) THEN
c-----------------------------------------------------------------------
c** For GPEF potential, read parameters defining the expansion variable
c                                 p    p      p      p
c                  y(R;k,a,b) = (R - Re )/(a*R + b*Re )
c=======================================================================
          READ(5,*) AGPEF(ISTATE), BGPEF(ISTATE)
c=======================================================================
          RREF(ISTATE)= -1.d0
          ENDIF
      WRITE(6,626) SLABL(ISTATE),RMIN(ISTATE),RMAX(ISTATE),RH(ISTATE)
c
c** Now to read in the trial dissociation energy and equilibrium
c   radial distance for the state.
c  De(s)   is the dissociation energy for each state.
c  Re(s)   is the equilibrium radial distance for each state.
c  IFDe(s) indicates whether the dissociation energy will be:
c        = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c  IFRe(s) indicates whether the equilibrium radial distance will be:
c           = 1: held fixed at read-in values.
c          <= 0: determined from fits.
c=======================================================================
      READ(5,*) DE(ISTATE), IFXDE(ISTATE)
      READ(5,*) RE(ISTATE), IFXRE(ISTATE)
c=======================================================================
      IF(PSEL(ISTATE).EQ.4) IFXDE(ISTATE)= 1
c=======================================================================
c** Read in parameters defining the exponent expansion \phi(r)
c  NSpow(s)  is the order of the phi_i exponent expansion for  R < Re
c  NLpow(s)  is the order of the phi_i exponent expansion for R > Re
c  nPB(s)  is the power  p  for  beta(r)  basic exponent variable
c  nQB(s)  is the power  q  for  beta(r)  exponent expansion variable
c  RREF(s)  defines the reference distance in the potential exponent 
c    expansion variable:  * for  RREF.le.0 , define parameter  RREF = Re
c      * for  RREF.gt.0 , fix parameter  RREF   at its read-in value
c=======================================================================
      READ(5,*) NSpow(ISTATE), NLpow(ISTATE), nPB(ISTATE), nQB(ISTATE),
     1                                                     RREF(ISTATE)
c=======================================================================
      IF((PSEL(ISTATE).EQ.2).OR.(PSEL(ISTATE).EQ.3)) THEN
          MMN= MMLR(NCMM(ISTATE),ISTATE)- MMLR(1,ISTATE)
c         IF((NCMM(ISTATE).GT.1).AND.(MMLR(1,ISTATE).NE.3).AND.
c    1                  ((nPB(ISTATE).LE.MMN).OR.(MMN.LE.0))) THEN
c             WRITE(6,628) nPB(ISTATE),MMN
c             STOP
c             ENDIF
          ENDIF
      IF(NSpow(ISTATE).GE.0) THEN
          npow= MAX(NSpow(ISTATE),NLpow(ISTATE))
          DO  I= 0,npow
c** Read in trial initial trial parameters for exponent \phi(r)
c
c  PHI(b,s)   contains the parameters that define the Born-Oppenheimer
c              potential functions:
c        PSEL .ne. 4 : read-in values are  phi_i  parameters for the
c               exponent expansion of the 'clamped nuclei' potential.
c        PSEL = 4 :    read-in values are leading coefficients in 
c                Surkus' Generalized Potential Energy Function (GPEF).
c  IFXPHI(b,s) indicates whether each exponent potential expansion 
c             coefficient will be:    = 1: held fixed at read-in values.
c                                  .LE. 0: determined from fits.
c** For NSpow .ge.0  exponent i\beta(yp) is 'conventional Huang-type 
c                 constrained polynomial in  yp.
c   For NSpow < 0  exponent is natural spline function defined at the yp 
c                  values ypPHI and fixed to equal ypINF at  ypPHI=1
c=======================================================================
              READ(5,*) PHI(I,ISTATE), IFXPHI(I,ISTATE)
c=======================================================================
              ENDDO
        ELSE
          npow= NLpow(ISTATE)
          DO  I= 1,npow
c-----------------------------------------------------------------------
c** For NSpow < 0  exponent is a natural spline function with values PHI
c       at the yp values ypPHI, and fixed to equal ypINF at  ypPHI=1
c=======================================================================
              READ(5,*) ypBETA(I,ISTATE),PHI(I,ISTATE),IFXPHI(I,ISTATE)
c=======================================================================
              ENDDO
          NLpow(ISTATE)= NLpow(ISTATE)+1
          ypBETA(NLpow(ISTATE),ISTATE)= 1.d0
          IFXPHI(NLpow(ISTATE),ISTATE)= 1
        ENDIF
c=======================================================================
c** Read parameters defining the BOB adiabatic radial functions
c  NUA/NUB(s)  specifies the order of the polynomial in  yp  defining 
c              the adiabatic BOB function for atom A/B
c         if < 0   do not read in any adiabatic BOB function parameters
c  pAD(s)/qAD(s)  are the powers defining the expansion variables
c  RREFadA/B(s) define reference distances in radial expansion variables
c    yp(r;RREFad{A/B}) for adiabatic radial strength function(s) Sad^{A/B}
c   * for  RREFad{A/B} .le.0 , define parameter  RREFad{A/B} = Re
c   * for  RREFad{A/B} > 0   , fix parameter at this read-in value
c  UA/UB(a,s)   are the adiabatic BOB function expansion coefficients
c  IFXU(A/B)(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  UAinf/UBinf  is the limiting asymptotic value of uA(r)/uB(r), as per
c               Theochem paper [internally stored as  UA(NUA+1), etc.]
c  IFXUAinf/IFXUBinf  specifies whether (>0) or not (.le.0)  UAinf/UBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NUA(ISTATE),NUB(ISTATE),pAD(ISTATE),qAD(ISTATE),
     1                                 RREFadA(ISTATE),RREFadB(ISTATE)
      IF(NUA(ISTATE).GE.0) THEN
c... NOTE that  parameters  NUA(ISTATE)+1  are  UAinf & IFXUAinf ...
          NUA(ISTATE)= NUA(ISTATE)+ 1
          DO  I= 0, NUA(ISTATE)
              READ(5,*) UA(I,ISTATE), IFXUA(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UA(NUA(ISTATE),ISTATE)= 0.d0
              IFXUA(NUA(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF
c=======================================================================
      IF(NUB(ISTATE).GE.0) THEN
c... NOTE that  parameters  NUB(ISTATE)+1  are  UBinf & IFXUBinf ...
          NUB(ISTATE)= NUB(ISTATE)+ 1
          DO  I= 0, NUB(ISTATE)
              READ(5,*) UB(I,ISTATE), IFXUB(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              UB(NUB(ISTATE),ISTATE)= 0.d0
              IFXUB(NUB(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF 
c***********************************************************************
c** Read parameters defining the BOB non-adiabatic centrifugal functions
c** If NISTP= 1 , read only one set of non-adiabatic parameters
c
c  NTA/NTB(s)  specifies the order of the polynomial in  yp  defining
c              the non-adiabatic centrifugal BOB functions for atom A/B
c         if < 0   do not read in any non-adiabatic BOB parameters
c  pNA(s)/qNA(s)  are the powers defining the expansion variables
c  TA/TB(a,s)   are the non-adiabatic centrifugal BOB expansion coeffts
c  IFXTA/IFXTB(a,s) indicates whether each expansion coefficient is to be
c           > 0 :  held fixed at read-in value, or
c        .le. 0 :  varied in the fit
c  TAinf/TBinf  is the limiting asymptotic value of qA(r)/qB(r), as per 
c               Theochem paper [internally stored as  TA(NTA+1), etc.]
c  IFXTAinf/IFXTBinf  specifies whether (>0) or not (.le.0)  TAinf/TBinf 
c               is to be held fixed at the read-in value
c=======================================================================
      READ(5,*) NTA(ISTATE),NTB(ISTATE),pNA(ISTATE),qNA(ISTATE)
      IF(NTA(ISTATE).GE.0) THEN
c... NOTE that  parameters  NTA(ISTATE)+1  are  TAinf & IFXTAinf ...
          NTA(ISTATE)= NTA(ISTATE)+ 1
          DO  I= 0, NTA(ISTATE)
              READ(5,*) TA(I,ISTATE), IFXTA(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TA(NTA(ISTATE),ISTATE)= 0.d0
              IFXTA(NTA(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF 
c=======================================================================
      IF(NTB(ISTATE).GE.0) THEN
c... NOTE that  parameters  NTB(ISTATE)+1  are  TBinf & IFXTBinf ...
          NTB(ISTATE)= NTB(ISTATE)+ 1
          DO  I= 0, NTB(ISTATE)
              READ(5,*) TB(I,ISTATE), IFXTB(I,ISTATE)
              ENDDO
c=======================================================================
          IF(BOBCN(ISTATE).GT.0) THEN
              TB(NTB(ISTATE),ISTATE)= 0.d0
              IFXTB(NTB(ISTATE),ISTATE)= 1
              ENDIF
          ENDIF
c
      NwCFT(ISTATE)= -1
      IF((IOMEG(ISTATE).GT.0).OR.(IOMEG(ISTATE).EQ.-1)) THEN
c-----------------------------------------------------------------------
c** If electronic angular momentum not zero for this state, read Lambda 
c       doubling or doublet Sigma radialfunction parameters.
c*  NwCFT(s)  is order of the polynomial representing the radial fx.
c*  nPw(s)  defined nature of radial expansion variable:  
c         yp = [R^p - Re^p]/[R^p + Re^p]    with  p= nPw
c*  efREF(s) defines reference level for the Lambda doubling splitting
c            = -1  treats f level as the reference
c            =  0  treats the mid-point between e and f as reference
c            =  1  treats e level as the reference
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          READ(5,*) NwCFT(ISTATE), nPw(ISTATE), efREF(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          IF(NwCFT(ISTATE).GE.0) THEN
c... NOTE that  parameters  NwCFT(ISTATE)+1  are  wCFTinf & IFXwCFTinf 
              NwCFT(ISTATE)= NwCFT(ISTATE)+ 1
              DO  I= 0, NwCFT(ISTATE)
                  READ(5,*) wCFT(I,ISTATE), IFXwCFT(I,ISTATE)
                  ENDDO
c=======================================================================
              IF(IOMEG(ISTATE).GT.0) THEN
                  IF(EFREF(ISTATE).EQ.-1) WRITE(6,640) SLABL(ISTATE)
                  IF(EFREF(ISTATE).EQ.0) WRITE(6,642) SLABL(ISTATE)
                  IF(EFREF(ISTATE).GE.1) WRITE(6,644) SLABL(ISTATE)
                  ENDIF
              ENDIF
          ENDIF
c
c** Calculate BOB mass scaling factors for the adiabatic (ZMUA, ZMUB) &
c   non-adiabaric centrifugal (ZMTA, ZMTB) BOB functions
      DO  IISTP= 1, NISTP
          IF(BOBCN(ISTATE).GE.1) THEN
c** For Watson/Coxon/Ogilvie-type clamped-nuclei reference species:
              ZMUA(IISTP,ISTATE)= ZMASE/ZMASS(1,IISTP)
              ZMUB(IISTP,ISTATE)= ZMASE/ZMASS(2,IISTP)
              ZMTA(IISTP,ISTATE)= ZMASE/ZMASS(1,IISTP)
              ZMTB(IISTP,ISTATE)= ZMASE/ZMASS(2,IISTP)
            ELSE
c** Using RJL's mass differences for adiabatic corrections (ZMUA, ZMUB):
              ZMUA(IISTP,ISTATE)= 1.0d0 - ZMASS(1,1)/ZMASS(1,IISTP)
              ZMUB(IISTP,ISTATE)= 1.0d0 - ZMASS(2,1)/ZMASS(2,IISTP)
c   and mass ratios for the rotational corrections (ZMTA, ZMTB):
              ZMTA(IISTP,ISTATE)= ZMASS(1,1)/ZMASS(1,IISTP)
              ZMTB(IISTP,ISTATE)= ZMASS(2,1)/ZMASS(2,IISTP)
            END IF
c
c** For homonuclear diatomics, set the first mass scaling term for each
c   set of correction terms to be the sum of the two original mass
c   scaling factors, and set the second mass term to zero.
c
          IF(AN(1).EQ.AN(2)) THEN
              ZMUA(IISTP,ISTATE)= ZMUA(IISTP,ISTATE)+ ZMUB(IISTP,ISTATE)
              ZMUB(IISTP,ISTATE)= 0.0d0
              ZMTA(IISTP,ISTATE)= ZMTA(IISTP,ISTATE)+ ZMTB(IISTP,ISTATE)
              ZMTB(IISTP,ISTATE)= 0.0d0
              END IF
          ENDDO
c-----------------------------------------------------------------------
      RETURN
  600 FORMAT(/'For state ',A2,' use a fixed potential defined by LEVEL s
     1ubroutine PREPOT'/4x,'Integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  604 FORMAT(/' For state ',A2,' represent level energies by independent
     1 term values')
  626 FORMAT(/' For state   ',A2/4x,'integrate from   RMIN=',f5.2,
     1  '   to   RMAX=',f6.2,'   with mesh   RH=',f8.5)
  628 FORMAT(' ** Since  p=',i2,' .LE. [MMLR(NCMM)-MMLR(1)]=',i2,
     1  '  or   [MMLR(NCMM)-MMLR(1)].le.0,   STOP !!')
  640 FORMAT(/' ', A2,' state energies referenced to f-parity levels')
  642 FORMAT(/' ', A2,' state energies referenced to the mid-point betwe
     1en e and f-parity levels')
  644 FORMAT(/' ', A2,' state energies referenced to e-parity levels')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
