c***********************************************************************
      SUBROUTINE DVDDP(DATN,ISTATE,IISTP,NPPS,ZMU,VD,DVDPK)
      IMPLICIT NONE
c***********************************************************************
c** This subroutine calculates  d(vD)/dp 
c
c** On entry:
c    ISTATE  is the molecular state being considered.
c    IISTP   is the isotopomer being considered.
c    NPPS(i) are the number of parameters varied in each state.
c    ZMU     is the reduced mass of the diatom in atomic units.
c
c** On exit:
c    DVDPK(i) are the values of the partial derivative d(vD)/dP(k).
c
c** Flags: Use only when debugging.
c
      INCLUDE 'arrsizes.h'
c
      INTEGER DATN, ISTATE, IISTP, NPPS(NSTATEMX),
     1   KV, IPV, I, K, ISCHR
c
      REAL*8 ZMU, BFCT, VL, ASCAT, VDFRAC, VD1, VD2, VD, DELTA,
     2   SWF(NPNTMX), V(NPNTMX), DVDP(NPNTMX), DVDPK(NPARMX)
c
c 
c** Type statements & common block for data
c
c     REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
c    1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
c    2  RMUP(0:8,NISTPMX)
c     INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
c    1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
c    2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
c    3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
c    4  ILAST(NBANDMX)
c     CHARACTER*2 NAME(2)
c
c     CHARACTER*2 SLABL(-3:NSTATEMX)
c     COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
c    1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
c    2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)
c
      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)
c
      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR
c
c** Born-Oppenheimer Breakdown Functions.
c     
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1  QAR(NPNTMX,NSTATEMX),QBR(NPNTMX,NSTATEMX),
     2  ZMUA(NISTPMX),ZMUB(NISTPMX),ZMQA(NISTPMX),ZMQB(NISTPMX)
      
      COMMON /BFBLK/UAR,UBR,QAR,QBR,ZMUA,ZMUB,ZMQA,ZMQB
c
c** Partial derivative arrays.
c
      REAL*8 DVBODP(NPNTMX,NPOTMX),DUADP(NPNTMX,NPOTMX),
     1  DUBDP(NPNTMX,NPOTMX),DQADP(NPNTMX,NPOTMX),DQBDP(NPNTMX,NPOTMX)
     2  ,dLdP(NPNTMX,NPOTMX),DBDP(NPNTMX,NPOTMX)
c
      COMMON /PDVBLK/DVBODP,DUADP,DUBDP,DQADP,DQBDP,DBDP
      COMMON /PdVLBLK/dLdP
c
      BFCT = ZMU / 16.85762908d0
c
c** Generating potential function including adiabatic BOB terms
c
c     V(R) = (2*mu/hbar^2)*(VPOT(R)+U(R))
c          = BFCT*(VPOT+ZMUA*UAR+ZMUB*UBR)
c
      DO I= 1,NDATPT(ISTATE)
          V(I) = BFCT* (VPOT(I,ISTATE) + ZMUA(IISTP)* UAR(I,ISTATE)
     1                                 + ZMUB(IISTP)* UBR(I,ISTATE))
          END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c** First determine which partials need to be calculated
          K = 0
          DO I = 1,ISTATE - 1
              K = K + NPPS(I)
              END DO
c
c** For DVDP - partials with respect to each parameter PV(IPV).
c            - including adiabatic correction terms.
          DO IPV = K+1,K+NPPS(ISTATE)
              DO I = 1,NDATPT(ISTATE)
                  DVDP(I) = BFCT * (DVBODP(I,IPV) 
     1                                    + ZMUA(IISTP)*DUADP(I,IPV)
     1                                    + ZMUB(IISTP)*DUBDP(I,IPV))
              END DO
c
              DELTA=1.D-3
              DO I = 1,NDATPT(ISTATE)
                  V(I) = V(I) + DELTA*DVDP(I)
              END DO
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              VL = BFCT*VLIM(ISTATE)
              CALL SCATLN(SWF,RH(ISTATE),RMIN(ISTATE),V,NDATPT(ISTATE),
     1          VL,NCN(ISTATE),KV,ASCAT,VDFRAC,ISCHR)
              VD1=DBLE(KV-1)+VDFRAC
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              DO I = 1,NDATPT(ISTATE)
                  V(I) = V(I) - 2.D0*DELTA*DVDP(I)
              END DO
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              VL = BFCT*VLIM(ISTATE)
              CALL SCATLN(SWF,RH(ISTATE),RMIN(ISTATE),V,NDATPT(ISTATE),
     1          VL,NCN(ISTATE),KV,ASCAT,VDFRAC,ISCHR)
              VD2=DBLE(KV-1)+VDFRAC
c
              DO I = 1,NDATPT(ISTATE)
                  V(I) = V(I) + DELTA*DVDP(I)
              END DO
c
              VD=(VD1+VD2)/2.D0
              DVDPK(IPV)=(VD1-VD2)/(2.D0*DELTA)
c             write(6,*) 'For parm ', ipv,', vD =',vd1,vd2,dvdpk(ipv) 
          END DO
C
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
