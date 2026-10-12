c=======================================================================
c** Isotope/isotopologue numbers, masses & BOB mass scaling factors
      INTEGER NISTP,NDUNMX,AN(2),MN(2,NISTPMX)
c** NDUNMX is a dummy parameter reqd. for portability of READATA
      PARAMETER (NDUNMX=0)
c
      REAL*8  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     1 RMUP(0:9,NISTPMX),ZMUA(NISTPMX,NSTATEMX),ZMUB(NISTPMX,NSTATEMX),
     2 ZMTA(NISTPMX,NSTATEMX),ZMTB(NISTPMX,NSTATEMX)
c
      COMMON /BLKISOT/ZMASS,RSQMU,RSQMUP,RMUP,ZMUA,ZMUB,ZMTA,ZMTB,
     1  NISTP,AN,MN
c=======================================================================
