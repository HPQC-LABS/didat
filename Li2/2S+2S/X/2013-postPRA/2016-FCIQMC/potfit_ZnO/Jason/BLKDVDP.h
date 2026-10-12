c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Partial derivative arrays for fits and uncertainties (fununc)
c
c=======================================================================
      REAL*8 DVtot(HPARMX,NPNTMX),DLDDRe(NPNTMX,NSTATEMX),
     1  DUADRe(NPNTMX,NSTATEMX),DUBDRe(NPNTMX,NSTATEMX),
     2  DTADRe(NPNTMX,NSTATEMX),DTBDRe(NPNTMX,NSTATEMX),
     3  DBDB(0:NphiMX,NPNTMX,NSTATEMX),DBDRe(NPNTMX,NSTATEMX)
      COMMON/BLKDVDP/DVtot,DUADRe,DUBDRe,DTADRe,DTBDRe,DLDDRe,DBDB,DBDRe
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
