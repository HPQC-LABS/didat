c** Utility routine to summarize dimensioning of arrays
c
      INTEGER NISTPMX,NPARMX,NBETAMX,NABOMX,NRNAMX,NPOTMX,NDATAMX,
     1  NBANDMX,NBNDMX,NVIBMX,NSTATEMX,MORDRMX,NPNTMX,NDUNMX
     2  ,numCnMax

c*  NISTPMX  is the maximum number of isotopomers allowed for fit
      PARAMETER (NISTPMX = 5)
c*  NPARMX  is the largest number of free parameters allowed for fit
      PARAMETER (NPARMX  = 2000)
c*  NBETAMX  is the largest number of exponent parameters allowed for fit
      PARAMETER (NBETAMX  = 21)
c*  NABOMX  is the largest number of adiabatic Born-Oppenheimer breakdown 
c           parameters per atom allowed for fit
      PARAMETER (NABOMX  = 7)
c*  NRNAMX  is the largest number of rotational non-adiabatic 
c           Born-Oppenheimer breakdown parameters per atom allowed for fit
      PARAMETER (NRNAMX  = 7)
c*  NPOTMX  is the largest number of potential parameters allowed.
c           NPOTMX >= 2 + (NBETAMX+1) + 2*(NABOMX+1) + 2*(NRNAMX+1)
      PARAMETER (NPOTMX = 55)
c*  NDATAMX  is largest No. of individual data which may be considered
      PARAMETER (NDATAMX = 10000)
c*  NBANDMX  is largest No. of bands/series which may be considered
      PARAMETER (NBANDMX = 1500)
c*  NBNDMX  band-by-band fits allow for up to NBNDMX data in a band
      PARAMETER (NBNDMX  = 40)
c*  NSTATEMX  is maximum no. of electronic states which can be
c             simultaneously fitted to
      PARAMETER (NSTATEMX = 2)
c*  NVIBMX  is the maximum number of vibrational levels of a single
c           state for which data are to be considered
      PARAMETER (NVIBMX    = 120)
c*  MORDRMX  is the largest number of centrifugal distortion constants
c            which may be considered
      PARAMETER (MORDRMX   = 7)
c*  NPNTMX  is the largest number of potential data points that can be
c           stored
      PARAMETER (NPNTMX = 38000)
c*  NDUNMX  is the largest number of dunham constants which may be
c           considered
      PARAMETER (NDUNMX = 16)
c*  numCnMax is the largest number of Cn values in the long-range
c            potential
      PARAMETER (numCnMax = 6)
