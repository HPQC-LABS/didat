c** Utility routine to summarize dimensioning of arrays
c
      INTEGER NISTPMX,NPARMX,NphiMX,NBOBMX,NPOTMX,NDATAMX,
     1  NBANDMX,NBNDMX,NVIBMX,NSTATEMX,NPNTMX,MORDRMX,numCmMax 

c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX =  4)

c*  NPARMX  is the largest number of free parameters allowed for fit
      PARAMETER (NPARMX  = 2500)

c*  NphiMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NphiMX  = 22)

c*  NBOBMX-1  is the highest-order polynomial expansion allowed for the
c      adiabatic or centrifugal Born-Oppenheimer breakdown functions, or 
c      the Lambda-doubling or 2\Sigma splitting radial strength functions
      PARAMETER (NBOBMX  = 10)

c*  NPOTMX  is the largest number of radial function parameters of all 
c    types (potential energy, BOB. etc.) allowed for a given state.
c           NPOTMX >= 2 + (NphiMX+1) + 5*(NBOBMX)
      PARAMETER (NPOTMX = 55)

c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 15000)

c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 2500)

c*  NBNDMX  band-by-band fits allow for up to NBNDMX data in a band
      PARAMETER (NBNDMX  = 20)

c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 2)

c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 120)

c*  NPNTMX  is the largest number of potential data points that can be
c           stored
      PARAMETER (NPNTMX = 60000)

c*  MORDRMX  is the largest number of centrifugal distortion constants
c            which may be considered
      PARAMETER (MORDRMX = 7)

c*  numCmMax is the largest number of Cm terms in the DELR long-range
c            potential
      PARAMETER (numCmMax = 6)

