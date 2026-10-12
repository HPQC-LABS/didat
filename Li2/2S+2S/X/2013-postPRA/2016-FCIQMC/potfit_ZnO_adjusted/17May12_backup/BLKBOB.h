c=======================================================================
c** Born-Oppenheimer Breakdown & doubling function parameters.
c**                       March 16 2012
c=======================================================================
      INTEGER NUA(NSTATEMX),NUB(NSTATEMX),NTA(NSTATEMX),NTB(NSTATEMX),
     1  IFXUA(0:NBOBMX,NSTATEMX),IFXUB(0:NBOBMX,NSTATEMX),
     2  IFXTA(0:NBOBMX,NSTATEMX),IFXTB(0:NBOBMX,NSTATEMX),
     3  NwCFT(NSTATEMX),IFXwCFT(0:NBOBMX,NSTATEMX),efREF(NSTATEMX)
c
      REAL*8 RREFadA(NSTATEMX),RREFadB(NSTATEMX),RREFnaA(NSTATEMX),
     1  RREFnaB(NSTATEMX),RREFw(NSTATEMX),UA(0:NBOBMX,NSTATEMX),
     1  UB(0:NBOBMX,NSTATEMX),TA(0:NBOBMX,NSTATEMX),
     2  TB(0:NBOBMX,NSTATEMX),wCFT(0:NBOBMX,NSTATEMX)
c
      COMMON /BLKBOB/RREFadA,RREFadB,RREFnaA,RREFnaB,RREFw,UA,UB,TA,TB,
     1  wCFT,NUA,NUB,NTA,NTB,NwCFT,IFXUA,IFXUB,IFXTA,IFXTB,IFXwCFT,efREF
c=======================================================================
