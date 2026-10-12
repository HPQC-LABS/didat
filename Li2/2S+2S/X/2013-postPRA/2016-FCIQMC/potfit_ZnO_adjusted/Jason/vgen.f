c***********************************************************************
      SUBROUTINE VGEN(ISTATE,RDIST,VDIST,PHIDIST,IDAT)
c***********************************************************************
c** This subroutine will generate one of six possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program. These four families are:
c     - The Extended Morse Oscillator
c     - The Modified Lennard-Jones Oscillator
c     - The Double-Exponential/Long-Range potential
c     - The Generalized Potential Energy Function
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2006-2009  by  R.J. Le Roy, Jenning Seto & Yiye Huang ++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           ----- Version of 22 August 2009 -----
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST: If RDIST > 0, calculate potl & derivs only @ that onee distance
c   -----    * return potential function at that point as VDIST and
c              potential function exponent as PHIDIST
c            * skip partial derivative calculation if  IDAT.le.0
c        * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks:
c PSEL   specifies the type of potential used.
c      PSEL = 0 : Use fixed potential function defined in READPOT.
c      PSEL = 1 : Use the Expanded Morse Oscillator(p) potential.
c      PSEL = 2 : Use the Morse/Lennard-Jones(p) potential.
c      PSEL = 3 : Use the Double-Exponential Long-Range Potential.
c      PSEL = 4 : Use the Generalized Potential Energy Function.
c      PSEL = 5 : Use the [Tiemann Polynomial]                  
c  NSpow(s)  is the order of the  phi(r)  exponent expansion for R.le.Re
c  NLpow(s)  is the order of the  phi(r)  exponent expansion for R > Re
c!!!  The code (currently) assumes  NLpow .ge. NSpow !!!!
c    MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c    nPB(s)  the basic value of power p for the beta(r)  exponent function
c    nQB(s)  the power p for the power series expansion variable in beta(r)
c    pAD(s) & qAD(s) the values of power p for adiabatic u(r) BOB functions
c    nNA(s) & qNA(s) the values of power p for centrifugal q(r) BOB functions
c    nPw(s)  the value of power m for Lambda(r) function
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    PHI    is the array of potential (exponent) expansion parameters 
c    NDATPT is the number of meshpoints used for the array.
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    R      is the distance array
c    VPOT   is the potential that is generated.
c    PHIFX  is used to contain the  phi(r)  function.
c** Internal partial derivative arrays ...
c    DUADRe & DUBDRe are p.derivs of adiabatic fx. w.r.t.  Re
c    DVDQA & DVDQB are p.derivs of non-adiabatic fx. wrt q_A(i) & q_B(i)
c    DTADRe & DTBDRe are p.derivs of non-adiabatic fx. w.r.t.  Re
c    dVdL & dLDDRe are p.derivatives of f_\lambda(r) w.r.t. phi_i & Re
c    DBDB & DBDRe are p.derives of phi(r) w.r.t. \phi_i & Re, respectively
c
c** Temp:
c    BTEMP     is used to represent the sum used for dV/dRe.
c              is used in GPEF for De calculations.
c    BINF      is used to represent the  beta(\infty)  value.
c    YP        is used to represent (R^p-Re^p)/(R^p+Re^p) or R-Re.
c    XTEMP     is used to represent (uLR/uLR_e)* exp{-BINF*RTEMP}
c    PBTEMP    is used to calculate dV/dBi.
c    PETEMP    is used to calculate dV/dBi.
c    AZERO     is used for the trial exponential calculations.
c    AONE      is used for the trial exponential calculations.
c    ATWO      is used for the trial exponential calculations.
c    AZTEMP    is used in the MMO trial exponential calculations.
c              is used in the GPEF = (a+b)/k
c    AOTEMP    is used in the GPEF = [a(k+1)-b(k-1)]/k
c    ATTEMP    is used in the GPEF = [a^2(k+1)-b^2(k-1)]/k
c    ARTEMP    is used in the GPEF = [a^3(k+1)-b^3(k-1)]/k
c    FSW       is used to represent the MLJ switching function.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
c-----------------------------------------------------------------------
c** Define local variables ...
      INTEGER MMAX
      PARAMETER (MMAX=20)
      INTEGER I,J,I1,ISTATE,IPV,IPVSTART,ISTART,ISTOP,LAMB2,m,npow,
     1 nread,IDAT, NBAND, IISTP, NPP
      REAL*8 BTEMP,BINF,RVAL,RTEMP,RM2,XTEMP,PBTEMP,PETEMP,RET,
     1 FSW,Xtemp2,Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,C3VAL,C3bar,C6bar,
     2 C6adj,C9adj,YP,YQ,YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,YQPA,
     3 YQPB,REp,REq,RDp,RDq,DYPDRE,DYQDRE,VAL,DVAL,HReP,HReQ,SL,SLB,
     4 AREF,AREFp,AREFq, RE3,RE6,RE8,T0,T0P,T1,ULRe,Scalc,dLULRedCm(9),
     5 dLULRedRe,dLULRedDe,dULRdDe,dULRdCm(9),RD3,RD6,RD8,DVDD,RDIST,
     6 VDIST,PHIDIST,BFCT,JFCT,JFCTLD,xBETA(NphiMX),rKL(NphiMX,NphiMX),
     7 RETSig,RETPi,RETp,RETm,AREFadA,AREFadB,AREFpADA,AREFpADB,
     8 AREFqADA,AREFqADB,BT,Rinn,Rout,A1,A2,A3,B1,B2,B3,B4,B5
c***********************************************************************
c** Common block for partial derivatives at the one distance RDIST
      REAL*8 dVdPk(HPARMX)
      COMMON /dVdPkBLK/dVdPk
c***********************************************************************
c** Temporary variables for MLR and DELR potentials
      INTEGER MMLRP
      REAL*8 ULR,dAAdRe,dBBdRe,dVdBtemp,CmVALL,tDm,tDmp,tDmpp,Dm(MMAX),
     1  Dmp(MMAX),Dmpp(MMAX)
c***********************************************************************
c** Initializing variables.
      WD(ISTATE)= 0
      REP= RE(ISTATE)**nPB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFp= AREF**nPB(ISTATE)
      AREFq= AREF**nQB(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= NDATPT(ISTATE)
c** When calculating only one potential point
      IF(RDIST.GT.0) THEN
          ISTART= NPNTMX
          ISTOP= NPNTMX
          VDIST= 0.0d0
          PHIDIST= 0.d0
          ENDIF
      PBTEMP= 0.0d0
      PETEMP= 0.0d0
      DO  I= ISTART,ISTOP
          PHIFX(I,ISTATE)= 0.0d0
          UAR(I,ISTATE)= 0.d0
          UBR(I,ISTATE)= 0.d0
          TAR(I,ISTATE)= 0.d0
          TBR(I,ISTATE)= 0.d0
          UAR(I,ISTATE)= 0.d0
          WRAD(I,ISTATE)= 0.d0
          ENDDO
c** Initialize parameter counter for this state ...
      IPVSTART= POTPARI(ISTATE) - 1

      IF(PSEL(ISTATE).EQ.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For forward data simulation using a known input potential 
c?????? Most/all of this section is dummy ?!?!?!?!?!?!?!?!?!/
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
c             RTEMP= (RVAL**nPB(ISTATE)-REP)/(RVAL**nPB(ISTATE)+REP)
c             PHIFX(I,ISTATE)= PHI(0,ISTATE)
c             DO  J= 1,NLpow(ISTATE)
c                 PHIFX(I,ISTATE)=PHIFX(I,ISTATE)+PHI(J,ISTATE)*RTEMP**J
c                 ENDDO
c             XTEMP= DEXP(-PHIFX(I,ISTATE)*(RVAL-RE(ISTATE)))
c** Now to calculate the actual potential
c             VPOT(I,ISTATE)= DE(ISTATE)*(1.0d0-XTEMP)**2
c    1                   + VLIM(ISTATE) - DE(ISTATE)
              VDIST= VPOT(I,ISTATE)
              ENDDO
          ENDIF
******7***** End Forward Calculation ***********************************

      IF(PSEL(ISTATE).EQ.1) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Expanded Morse Oscillator - with  NSpow/NLpow
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First ... calculate the Generalized Morse Oscillator exponent
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**nPB(ISTATE)
              YP= (RDp - AREFp)/(RDp + AREFp)
              VAL= PHI(0,ISTATE) 
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              YPP= 1.d0
              npow= NSpow(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
              DO  J= 1,npow
                  DVAL= DVAL + PHI(J,ISTATE)* DBLE(J)* YPP
                  YPP= YPP*YP
                  VAL= VAL + PHI(J,ISTATE)*YPP 
                  DBDB(J,I,ISTATE)= YPP
                  ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
              DBDRe(I,ISTATE)= 0.d0
              PHIFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*(RVAL-RE(ISTATE)))
c** First calculate the partial derivative w.r.t. DE
              IPV= IPVSTART+ 1
              DVtot(IPV,I)= XTEMP*(XTEMP- 2.d0)
              DVDD= DVtot(IPV,I)
c** Now calculate the actual potential
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= VAL
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
c... now generate remaining partial derivatives
              YPP= 2.0d0 * DE(ISTATE) * XTEMP * (1.d0 - XTEMP) 
              IF(RREF(ISTATE).LE.0.d0) THEN
                  DBDRe(I,ISTATE)= -0.5d0*nPB(ISTATE)*(1.d0-YP**2)
     1                                                *DVAL/RE(ISTATE)
                  VAL= VAL - (RVAL- RE(ISTATE))*DBDRe(I,ISTATE)
                  ENDIF
              IPV= IPV+1
              DVtot(IPV,I)= -YPP*VAL
              YPP= YPP*(RVAL - RE(ISTATE))
              DO  J= 0,npow
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP
                  YPP= YPP*YP
                  ENDDO
              ENDDO
          ENDIF
ccccc Print for testing
cc      write(10,610) (RD(i,ISTATE),vpot(i,istate),PHIFX(i,istate),
cc   1                                        i= 1, NDATPT(ISTATE),50)
ccccc
c******** End preparation of Expanded Morse Potential Function *********

      IF(PSEL(ISTATE).EQ.2) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the  {Morse/Long-Range}_p  potential.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First - define values & derivatives of uLR at Re
          
          IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... for Aubert-Frecon treatment of {C3,C6,C8} for Alkali A-state
              RE3= 1.d0/RE(ISTATE)**3
              RE6= RE3*RE3
c** Include QED retardation function in C3 terms for Li2(A) !!
c  NOTE ... the numerical factor is  2\pi/\lambda  for this case
c  Ignore effect of retardation on the partial derivatives & C6adj & C9adj.
              RET= 9.364128d-4*RE(ISTATE)
              RETSig= DCOS(RET) + (RET)*DSIN(RET)
              RETPi= RETSig - RET**2 *DCOS(RET)
              RETp= RETSig + 0.5d0*RETPi
              RETm= RETSig - 0.5d0*RETPi
              C3VAL= CmVAL(1,ISTATE)
              C3bar= C3VAL/DE(ISTATE)
              C6bar= CmVAL(3,ISTATE)/DE(ISTATE)
              C6adj= CmVAL(3,ISTATE) + C3VAL**2/(4.d0*DE(ISTATE))
              C9adj= 0.5d0*C3VAL*C6adj/DE(ISTATE)
              T1= (C3VAL*RETm + (C6adj - CmVAL(4,ISTATE))*RE3)*RE3/3.d0
              IF(NCMM(ISTATE).GT.4) THEN
c ... extension for Aubert-Frecon Li2(A) {3,0,6,6,8,8} case ...
                  RE8= RE6/RE(ISTATE)**2
                  T1= T1+ (CmVAL(5,ISTATE)- CmVAL(6,ISTATE))*RE8/3.d0
                  ENDIF
              T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
              ULRe= 0.5d0*(-CmVAL(2,ISTATE) + RE3*(C3VAL*RETp
     1    + RE3*(C6adj + CmVAL(4,ISTATE)))) + RE3*RE6*C9adj + 0.5d0*T0
              IF(NCMM(ISTATE).GT.4) ULRe= ULRe
     1                 + 0.5d0*(CmVAL(5,ISTATE) + CmVAL(6,ISTATE))*RE8
c ... now collect derivative terms ...
              T0P= 0.5d0*(9.d0*T1-CmVAL(2,ISTATE))/T0
              dLULRedCm(1)= RE3*((9.d0 + T0P) + RE3*(C3bar*(3.d0+T0P)
     1               + RE3*(6.d0*C6bar + 4.5d0*C3bar**2)))/(12.d0*ULRe)
              dLULRedCm(3)= RE6*((3.d0+T0P)+3.0d0*RE3*C3bar)/(6.d0*ULRe)
              dLULRedCm(4)= (3.d0-T0P)*RE6/(6.d0 * ULRe)
              dLULRedRe= -RE3*(0.25d0*C3VAL*(9.d0+T0P) 
     1        + RE3*(C6adj*(3.d0+T0P) + CmVAL(4,ISTATE)*(3.d0-T0P)
     2                            + RE3*9.d0*C9adj))/(ULRe*RE(ISTATE))
              dLULRedDe= -RE6*((3.d0+T0P)*C3bar**2 
     1           - RE3*C3bar*(12.d0*C6bar + 6.d0*C3bar**2))/(24.d0*ULRe)
              IF(NCMM(ISTATE).GT.4) THEN
                  dLULRedCm(5)= RE8*(3.d0 + T0P)/(6.d0 * ULRe)
                  dLULRedCm(6)= RE8*(3.d0 - T0P)/(6.d0 * ULRe)
                  dLULRedRe= dLULRedRe 
     1                        - RE8*4.d0*(CmVAL(5,ISTATE)*(3.d0 + T0P)
     2          + CmVAL(6,ISTATE)*(3.d0 - T0P))/(3.d0*ULRe*RE(ISTATE))
                  ENDIF
            ELSE
c ... for ordinary NCMM-term MLR uLR(r)
              ULRe= 0.d0
              T1= 0.d0
              IF(IDF(ISTATE).GE.0) CALL dampF(RE(ISTATE),RHOab(ISTATE)
     1  ,NCMM(ISTATE),MMLR,MMLR(NCMM(ISTATE),
     2  ISTATE),ISTATE,IDF(ISTATE),WD(ISTATE),Dm,Dmp,Dmpp)
              DO  m= 1,NCMM(ISTATE)
                  IF(IDF(ISTATE).LT.0.d0) THEN
                      dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
                    ELSE
                      dLULRedCm(m)= Dm(MMLR(m,ISTATE))
     1                                     /RE(ISTATE)**MMLR(m,ISTATE)
                    ENDIF
                  T0= CmVAL(m,ISTATE)*dLULRedCm(m)
                  ULRe= ULRe + T0
                  T1= T1 + MMLR(m,ISTATE)*T0
                  ENDDO
              DO  m= 1,NCMM(ISTATE)
                  dLULRedCm(m)= dLULRedCm(m)/ULRe
                  ENDDO
              dLULRedRe= -T1/(ULRe*RE(ISTATE))
            ENDIF
          BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
          betaINF(ISTATE)= BINF
          IF(NSpow(ISTATE).LT.0) THEN
c*** For Pashov-natural-spline exponent coefficient ...
              nQB(ISTATE)=nPB(ISTATE)
              DO  I= 1,NLpow(ISTATE)
                  xBETA(I)= ypBETA(I,ISTATE)
                  ENDDO
              PHI(NLpow(ISTATE),ISTATE)= BINF
              CALL Lkoef(NLpow(ISTATE),xBETA,rKL,NphiMX)
              ENDIF
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**nPB(ISTATE)
              RDq= RVAL**nQB(ISTATE)
              YPE= (RDp-REP)/(RDp+REP)
              YP= (RDp-AREFp)/(RDp+AREFp)
              YQ= (RDq-AREFq)/(RDq+AREFq)
              YPM= 1.d0 - YP
              DYPDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YP**2)/RE(ISTATE)
              DYQDRE= -0.5d0*nQB(ISTATE)*(1.d0 - YQ**2)/RE(ISTATE)
              YPP= 1.d0
              DVAL= 0.d0
              DBDB(0,I,ISTATE)= 1.0d0
              IF(NSpow(ISTATE).GE.0) THEN
c** For 'conventional' power series or A-F exponent function ...
                  VAL= PHI(0,ISTATE)
                  npow= NSpow(ISTATE)
                  IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
                  DO  J= 1,npow
c... now calculate power series part of the Morse-like exponent
                      DVAL= DVAL + PHI(J,ISTATE)* DBLE(J)* YPP
                      YPP= YPP*YQ
                      VAL= VAL + PHI(J,ISTATE)*YPP
                      DBDB(J,I,ISTATE)= YPM*YPP
                      ENDDO
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
                  DBDRe(I,ISTATE)= -YP*dLULRedRe
                  IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= 
     1                            DBDRe(I,ISTATE)+ (BINF - VAL)*DYPDRE 
     2                                         + (1.d0-YP)*DVAL*DYQDRE
                  VAL= YP*BINF + (1.d0- YP)*VAL
                ELSE
c... now calculate Pashov-spline exponent coefficient & its derivatives
                  VAL= 0.d0
                  npow= NLpow(ISTATE)
                  DO  J= 1,npow
                      DBDB(J,I,ISTATE)= 
     1                               Scalc(YQ,J,npow,xBETA,rKL,NphiMX)
                      VAL= VAL+ DBDB(J,I,ISTATE)*PHI(J,ISTATE)
                      ENDDO
                  DBDRe(I,ISTATE)= -DBDB(npow,I,ISTATE)*dLULRedRe
                ENDIF
              PHIFX(I,ISTATE)= VAL
              XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... generate ULR for Aubert-Frecon Li2(A) {3,0,6,6} type case ...
                  RD3= 1.d0/RVAL**3
                  RD6= RD3*RD3
c** Include QED retardation function in definition of C3 ... as above!!
                  RET= 9.364128d-4*RVAL
                  RETSig= DCOS(RET) + RET*DSIN(RET)
                  RETPi= RETSig - RET**2 *DCOS(RET)
                  RETp= RETSig + 0.5d0*RETPi
                  RETm= RETSig - 0.5d0*RETPi
                  T1= (C3VAL*RETm + (C6adj
     1                                - CmVAL(4,ISTATE))*RD3)*RD3/3.d0
                  IF(NCMM(ISTATE).GT.4) THEN
c ... extension for Li2(A) {3,0,6,6,8,8} case ...
                      RD8= RD6/RVAL**2
                      T1=T1+(CmVAL(5,ISTATE)-CmVAL(6,ISTATE))*RD8/3.d0
                      ENDIF
                  T0= DSQRT((T1- CmVAL(2,ISTATE))**2 + 8.d0*T1**2)
                  ULR= 0.5d0*(-CmVAL(2,ISTATE)+ (C3VAL*RETp
     1                + (C6adj + CmVAL(4,ISTATE))*RD3)*RD3) + 0.5d0*T0
     2                                                 + C9adj*RD3*RD6
                  IF(NCMM(ISTATE).GT.4) ULR= ULR
     1                 + 0.5d0*(CmVAL(5,ISTATE) + CmVAL(6,ISTATE))*RD8
                  T0P= (9.d0*T1-CmVAL(2,ISTATE))/T0
                  dULRdCm(1)= RD3*((9.d0+T0P) + RD3*(C3bar*(3.d0+T0P)
     1                    + RD3*(6.d0*C6bar + 4.5d0*C3bar**2)))/12.d0
                  dULRdCm(2)= 0.d0
                  dULRdCm(3)= RD6*((3.d0+T0P) + 3.d0*RD3*C3bar)/6.d0
                  dULRdCm(4)= RD6*(3.d0-T0P)/6.d0
                  dULRdDe= -RD6*((3.d0+T0P)*C3bar**2
     1                + RD3*C3bar*(12.d0*C6bar + 6.d0*C3bar**2))/24.d0
                  IF(NCMM(ISTATE).GT.4) THEN
                      dULRdCm(5)= RD8*(3.d0+T0P)/6.d0
                      dULRdCm(6)= RD8*(3.d0-T0P)/6.d0
                      ENDIF
                ELSE
c ... for the case of a 'normal' MLR or MLJ function
                  ULR= 0.d0
                  IF(IDF(ISTATE).GE.0.d0) CALL dampF(RVAL,RHOab(ISTATE),
     1  NCMM(ISTATE),MMLR,MMLR(NCMM(ISTATE),ISTATE),ISTATE,
     2  IDF(ISTATE),WD(ISTATE),Dm,Dmp,Dmpp)
                  DO  m= 1,NCMM(ISTATE)
                      IF(IDF(ISTATE).LT.0.d0) THEN
                          dULRdCm(m)= 1.d0/RVAL**MMLR(m,ISTATE)
                        ELSE
                          dULRdCm(m)= Dm(MMLR(m,ISTATE))
     1                                           /RVAL**MMLR(m,ISTATE)
                        ENDIF
                      ULR= ULR + CmVAL(m,ISTATE)*dULRdCm(m)
                      ENDDO
                ENDIF
              XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is asymptote ...
              DVDD= XTEMP*(XTEMP - 2.D0)  
              VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= VAL
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
              YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
              IPV= IPVSTART+2
              IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0)) THEN
c ... derivatives w.r.t. long-range parameters for Aubert-Frecon  uLR
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(1)
     1                                               - dULRdCm(1)/ULR)
                  IPV= IPV+1
c... note that derivative w.r.t. atomic level splitting FIXED
                  DVtot(IPV,I)= 0.d0
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(3)
     1                                               - dULRdCm(3)/ULR)
                  IPV= IPV+1
                  DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(4)
     1                                               - dULRdCm(4)/ULR)
                  IF(NCMM(ISTATE).GT.4) THEN
c... and then derivatives w.r.t. Auber-Frecon C8 values ...
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(5)
     1                                               - dULRdCm(5)/ULR)
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP*((1.d0 - YP*YPE)*dLULRedCm(6)
     1                                               - dULRdCm(6)/ULR)
                      ENDIF
                ELSE
c ... derivative w.r.t. Cm's for ordinary MLR/MLJ case ...
                  DO  m= 1, NCMM(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= -YPP*(dLULRedCm(m)*(YP*YPE- 1.d0)
     1                                               + dULRdCm(m)/ULR)
                      ENDDO
                ENDIF
c... derivative w.r.t. Re  
              DVtot(IPVSTART+2,I)= YPP*(YPE*DBDRe(I,ISTATE)
     1                                       + VAL*DYPDRE + dLULRedRe)
              IF(NSpow(ISTATE).GE.0) THEN
c... derivative w.r.t. De  for 'conventional' power-series exponent
                  DVDD= DVDD + YPP*YP*YPE/DE(ISTATE)
                  IF((NCMM(ISTATE).GE.4).AND.(MMLR(2,ISTATE).EQ.0))
c... derivative w.r.t. De  for Aubert-Frecon 2x2 exponent
     1        DVDD= DVDD+ YPP*((1.d0- YP*YPE)*dLULRedDe - dULRdDe/ULR)
                  YPP= YPP*YPE*(1.d0 - YP)
                  DVtot(IPVSTART+1,I)= DVDD
c... finally ... derivatives w.r.t. exponent expansion coefficients
                  DO  J= 0,npow
                      IPV= IPV+1
                      DVtot(IPV,I)= YPP
                      YPP= YPP*YQ
                      ENDDO
                ELSE
c... for Pashov-spline exponent cases...
                  YPP= YPP*YPE
                  DO  J= 1,NLpow(ISTATE)
                      IPV= IPV+ 1
                      DVtot(IPV,I)= DBDB(J,I,ISTATE)*YPP
                      ENDDO
                  DVtot(IPVSTART+1,I)= DVDD
     1                             + YPP*DBDB(npow,I,ISTATE)/DE(ISTATE)
                ENDIF
              ENDDO
ccccc Print for testing
      rewind(10)
      write(10,610) (RD(i,ISTATE),vpot(i,istate),PHIFX(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
  610 FORMAT(/(f10.4,f15.5,f12.6))
ccccc End of Print for testing
          ENDIF
c********* End for Morse/Lennard-Jones(p) potential function ***********

      IF(PSEL(ISTATE).EQ.3) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Double-Exponential Long-Range Oscillator.
          WD(ISTATE)= 1
          AA(ISTATE)= 0.0d0
          BB(ISTATE)= 0.0d0
          dAAdRe= 0.0d0
          dBBdRe= 0.0d0
          dVdBtemp= 0.0d0
          CALL dampF(RE(ISTATE),RHOab(ISTATE),NCMM(ISTATE), 
     1  MMLR,MMLR(NCMM(ISTATE),ISTATE),ISTATE,IDF(ISTATE),WD(ISTATE),
     2  Dm,Dmp,Dmpp)
c... first, get  AA & BB and their derivatives!
          DO  m= 1,NCMM(ISTATE)
              MMLRP= MMLR(m,ISTATE)
              CmVALL= CmVAL(m,ISTATE)
              tDm= Dm(MMLR(m+1,ISTATE))
              tDmp= Dm(MMLR(m+1,ISTATE))
              tDmpp= Dm(MMLR(m+1,ISTATE))
              AA(ISTATE)= AA(ISTATE)- CmVALL*tDm/RE(ISTATE)**MMLRP
     2                              *(1.0d0 + tDmp/(PHI(0,ISTATE)*tDm)
     3                             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
              BB(ISTATE)= BB(ISTATE)- CmVALL*tDm/RE(ISTATE)**MMLRP
     2                              *(2.0d0 + tDmp/(PHI(0,ISTATE)*tDm)
     3                             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
              dAAdRe= dAAdRe - CmVALL*((tDmp*RE(ISTATE) - MMLRP*tDm)
     2          /RE(ISTATE)**(MMLRP+1)*(1.0d0+tDmp/(PHI(0,ISTATE)*tDm)
     3             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE))) 
     4             + tDm/RE(ISTATE)**MMLRP*((tDm*tDmpp-tDmp**2)
     5 /(PHI(0,ISTATE)*tDm**2) + MMLRP/(PHI(0,ISTATE)*RE(ISTATE)**2)))
              dBBdRe= dBBdRe - CmVALL*((tDmp*RE(ISTATE) - MMLRP*tDm)
     1         /RE(ISTATE)**(MMLRP+ 1)*(2.0d0+tDmp/(PHI(0,ISTATE)*tDm)
     2             - MMLRP/(PHI(0,ISTATE)*RE(ISTATE)))
     3             + tDm/RE(ISTATE)**MMLRP*((tDm*tDmpp-tDmp**2)
     4 /(PHI(0,ISTATE)*tDm**2) + MMLRP/(PHI(0,ISTATE)*RE(ISTATE)**2)))
              dVdBtemp= dVdBtemp- CmVALL*tDm/
     1 (RE(ISTATE)**MMLRP*PHI(0,ISTATE)**2)*(MMLRP/RE(ISTATE)-tDmp/tDm)
              ENDDO
          AA(ISTATE)= AA(ISTATE) + DE(ISTATE) 
          BB(ISTATE)= BB(ISTATE) + 2.0d0*DE(ISTATE)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          DO  I= ISTART,ISTOP
c** First, calculate the generalized Morse exponent
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**NPB(ISTATE)
              YP= (RDp-REP)/(RDp+REP)
              npow= NSpow(ISTATE)
              IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
              PHIFX(I,ISTATE)= PHI(0,ISTATE)
              DO  J= 1,npow
                  PHIFX(I,ISTATE)= PHIFX(I,ISTATE) +
     1                                         PHI(J,ISTATE)*YP**J
                  ENDDO
c** Calculate some temporary variables.
              XTEMP= DEXP(-PHIFX(I,ISTATE)*(RVAL-RE(ISTATE)))
              Xtemp2= XTEMP*XTEMP
              BMtemp= PHIFX(I,ISTATE)
              BMtemp2= PHIFX(I,ISTATE)
c
c** Now to calculate the actual potential and partial derivatives:
c
              ULR= 0.0d0
              CALL dampF(RVAL,RHOab(ISTATE),NCMM(ISTATE),
     1  MMLR, MMLR(NCMM(ISTATE),ISTATE),ISTATE,IDF(ISTATE),WD(ISTATE),
     2  Dm,Dmp,Dmpp)
              DO  m= 1,NCMM(ISTATE)
                  MMLRP= MMLR(m,ISTATE)
                  CmVALL= CmVAL(m,ISTATE)
                  ULR= ULR + CmVALL*Dm(MMLR(m,ISTATE))/RVAL**MMLRP
                  ENDDO
              VPOT(I,ISTATE)= (AA(ISTATE)*XTEMP - BB(ISTATE))*XTEMP 
     1                                            - ULR + VLIM(ISTATE)
              IF(RDIST.GT.0.d0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= PHIFX(I,ISTATE)
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
c
c** Now to calculate the partial derivatives
c
              IPV= IPVSTART + 1
c ... first, derivative of the potential w.r.t. De
              DVtot(IPV,I)= (XTEMP - 2.0d0)*XTEMP
              DVDD= (XTEMP - 2.0d0)*XTEMP 
              BTEMP= 0.0d0
              Btemp2= 0.0d0
              IF(RVAL .LT. RE(ISTATE)) THEN
                  DO  J= 1,NSpow(ISTATE)
                      Btemp2= Btemp2 + PHI(J,ISTATE)*DBLE(J)
     1                                                   *YP**(J-1)
                      ENDDO
                  BTEMP= Btemp2
                ELSE
                  DO  J= 1,NLpow(ISTATE)
                      Btemp2= Btemp2 + PHI(J,ISTATE)*DBLE(J) 
     1                                                   *YP**(J-1)
                      ENDDO
                  BTEMP= Btemp2
                ENDIF
              Btemp2= Btemp2*(RVAL-RE(ISTATE))
              BTEMP= BTEMP*(RVAL-RE(ISTATE))
c
c** Now to calculate the remaining partial derivatives.
c
              PETEMP= 2.0d0*nPB(ISTATE)*RE(ISTATE)**(nPB(ISTATE)-1)
     &                                          * RVAL**nPB(ISTATE) /
     &                 (RVAL**nPB(ISTATE)+RE(ISTATE)**nPB(ISTATE))**2
              PBTEMP= BB(ISTATE) * XTEMP * (RVAL-RE(ISTATE))
              PBtemp2= 2.0d0 * AA(ISTATE) * Xtemp2* (RVAL-RE(ISTATE))
              IPV= IPV + 1
c ... then, derivative of the potential w.r.t. Re
              DVtot(IPV,I)= 2.0d0*AA(ISTATE)*Xtemp2*
     1                       (BMtemp2+PETEMP* Btemp2) + Xtemp2*dAAdRe
     2        - BB(ISTATE)*XTEMP*(BMtemp+PETEMP*BTEMP) - XTEMP*dBBdRe
c ... finally, derivatives of the potential w.r.t. the \phi_i
              IPV= IPV+ 1
              DVtot(IPV,I)= -PBtemp2 + PBTEMP + dVdBtemp*(Xtemp2-XTEMP)
              DO  J= 1,NSpow(ISTATE)
                  IPV= IPV+ 1
                  DVtot(IPV,I)= (-PBtemp2 + PBTEMP) * YP**J
                  ENDDO
              DO  J= NSpow(ISTATE)+1, NLpow(ISTATE)
                  IPV= IPV+ 1
                  IF (RVAL .LT. RE(ISTATE)) THEN
                      DVtot(IPV,I)= 0.0d0
                    ELSE
                      DVtot(IPV,I)= PBTEMP * YP**J
                    ENDIF
                  ENDDO
c *** DBDRe and DBDB is used in uncertainty calculation, see fununc.f
              DBDRe(I,ISTATE)= -PETEMP*BTEMP/(RVAL-RE(ISTATE))
              DBDB(0,I,ISTATE)= 1.0d0
              DO  J= 1, NLpow(ISTATE)
                  DBDB(J,I,ISTATE)= DBDB(J-1,I,ISTATE)*YP
                  ENDDO
              ENDDO
          ENDIF
******7* End Double-Exponential Long-Range Potential Function ********

      IF(PSEL(ISTATE).EQ.4) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Surkus-type Generalized Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** First, we calculate the implied Dissociation Energy
        IF(AGPEF(ISTATE).NE.0.0d0) THEN
            YPP= 1.d0/AGPEF(ISTATE)**2
            VAL= YPP
            DO  I= 1, NLpow(ISTATE)
                YPP= YPP/AGPEF(ISTATE)
                VAL= VAL + PHI(I,ISTATE)*YPP
                ENDDO
            DE(ISTATE)= VAL*PHI(0,ISTATE)
            ENDIF
        DO  I= ISTART, ISTOP
            RVAL= RD(I,ISTATE)
            IF(RDIST.GT.0.d0) RVAL= RDIST
            RDp= RVAL**nPB(ISTATE)
            YP= (RDp-REP)/(AGPEF(ISTATE)*RDp + BGPEF(ISTATE)*REP)
c** Now to calculate the actual potential
            YPP= 1.d0
            VAL= 1.d0
            DVAL= 2.d0
            DO  J= 1, NLpow(ISTATE)
                YPP= YPP*YP
                VAL= VAL + PHI(J,ISTATE)*YPP
                DVAL= DVAL+ (J+2)*PHI(J,ISTATE)*YPP
                ENDDO
            VPOT(I,ISTATE)= VAL*PHI(0,ISTATE)*YP**2 + VLIM(ISTATE) 
            IF(RDIST.GT.0) THEN
                VDIST= VPOT(I,ISTATE)
                PHIDIST= 0.d0
                IF(IDAT.LE.0) GO TO 999
                ENDIF
            DVAL= DVAL*PHI(0,ISTATE)*YP
c** Now to calculate the partial derivatives
            DVDD= 0.d0
            IPV= IPVSTART + 1
c ... derivative of the potential w.r.t. Re
            DVtot(IPV,I)= -DVAL*REP*RDp*(AGPEF(ISTATE)+BGPEF(ISTATE))
     1                  *(nPB(ISTATE)/RE(ISTATE))/(AGPEF(ISTATE)*RDp +
     2                                           BGPEF(ISTATE)*REP)**2
c ... and derivatives w.r.t. the phi_i expansion coefficients ...
            IPV= IPV+ 1
            DVtot(IPV,I)= VAL*YP**2
            IPV= IPV+ 1
            DVtot(IPV,I)= PHI(0,ISTATE)*YP**3
            DO  J= 2, NLpow(ISTATE)
                IPV= IPV+ 1
                DVtot(IPV,I)= DVtot(IPV-1,I)*YP
                ENDDO
            ENDDO
        IF(RDIST.LE.0.d0) VLIM(ISTATE)= VPOT(NDATPT(ISTATE),ISTATE)
c????
        rewind(10)
        write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
  612 FORMAT(/(f10.4,f15.5))
c????
        ENDIF
c*********** End Generalized Potential Energy Function *****************

      IF(PSEL(ISTATE).EQ.5) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For the Tiemann Potential Energy Function.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          BT= PHI(NLPow(ISTATE)+1, ISTATE)
          Rinn= PHI(NLPow(ISTATE)+2, ISTATE)
          Rout= PHI(NLPow(ISTATE)+3, ISTATE)
c** Determine analytic function attaching smoothly to inner wall of 
c  polynomial expansion at  R= Rinn < Rm
          YP= (Rinn - RE(ISTATE))/(Rinn+ BT*RE(ISTATE))
          YPP= 1.d0
          A1= PHI(0,ISTATE) 
          A2= 0.d0
          DO  J= 1, NLpow(ISTATE)
              A2= A2+ J*YPP*PHI(J,ISTATE)
              YPP= YPP*YP
              A1= A1+ YPP*PHI(J,ISTATE)
              ENDDO
          A2= A2*(RE(ISTATE)+ BT*RE(ISTATE))/(Rinn + BT*RE(ISTATE))**2
c* If inward extrapolation is exponential:   A1*exp(-A2*(R-Rinn))
          A2= -A2/A1
c** With long-range tail an NCMM-term inverse-power sum, add 1 exponential
c   term to ensure continuity and smoothness at  Rout
          YP= (Rout - RE(ISTATE))/(Rout+ BT*RE(ISTATE))
          YPP= 1.d0
          B1= PHI(0,ISTATE) + VLIM(ISTATE) - DE(ISTATE)
          B2= 0.d0
          DO  J= 1, NLpow(ISTATE)
              B2= B2+ J*YPP*PHI(J,ISTATE)
              YPP= YPP*YP
              B1= B1+ YPP*PHI(J,ISTATE)
              ENDDO
          B2= B2*(RE(ISTATE)+ BT*RE(ISTATE))/(Rout + BT*RE(ISTATE))**2
          B3= VLIM(ISTATE)
          B4= 0.d0
          DO  J= 1, NCMM(ISTATE)
              B5= CmVAL(J,ISTATE)/Rout**MMLR(J,ISTATE)
              B3= B3- B5
              B4= B4+ MMLR(J,ISTATE)*B5
              ENDDO
          B3= B3- B1
          B4= (B2- B4/Rout)/B3
c ... now generate potential as a Tiemann-type expansion
cc              WRITE(6,698) IPVSTART+ 1
          DO  I= ISTART, ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0) RVAL= RDIST
              IF(RVAL.LE.Rinn) THEN
c ... for exponential inward extrapolation ...
                  VPOT(I,ISTATE)= A1*DEXP(-A2*(RVAL- Rinn)) 
     1                                     + VLIM(ISTATE) - De(ISTATE)
                ELSEIF(RVAL.LE.Rout) THEN
                  YP= (RVAL - RE(ISTATE))/(RVAL + BT*RE(ISTATE))
                  A3= PHI(0,ISTATE) + VLIM(ISTATE) - DE(ISTATE)
                  YPP= 1.d0
                  DO  J= 1,NLpow(ISTATE)
                      YPP= YPP*YP
                      A3= A3+ PHI(J,ISTATE)*YPP
                      ENDDO
                  VPOT(I,ISTATE)= A3
                ELSEIF(RVAL.GT.Rout) THEN
                  A3= VLIM(ISTATE)
                  DO  J= 1, NCMM(ISTATE)
                      A3= A3- CmVAL(J,ISTATE)/RVAL**MMLR(J,ISTATE)
                      ENDDO
                  VPOT(I,ISTATE)= A3 - B3*DEXP(-B4*(RVAL- Rout))
cc     1                                 + VLIM(ISTATE) - De(ISTATE)
                ENDIF
              IF(RDIST.GT.0) THEN
                  VDIST= VPOT(I,ISTATE)
                  PHIDIST= 0.d0
                  IF(IDAT.LE.0) GO TO 999
                  ENDIF
              IPV= IPVSTART+ 1
              DVtot(IPV,I)= -1.d0
              IPV= IPV+1
              DVtot(IPV,I)= 0.d0
              DO  J=1, NCMM(ISTATE)
                  IPV= IPV+1
                  DVtot(IPV,I)= 0.d0
                  ENDDO
              IPV= IPV+ 1
              DVtot(IPV,I)= 1.d0
              DO  J= 1, NLpow(ISTATE)
                 IPV= IPV+ 1
                 DVtot(IPV,I)= DVtot(IPV-1,I)*YP
                 ENDDO
              ENDDO

cc              WRITE(6,699) IPV

      rewind(10)
      write(10,612) (RD(i,ISTATE),vpot(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
          ENDIF
 698        FORMAT('IPV before loop:'I5)
 699        FORMAT('IPV after loop:'I5)
c*********** End Tiemann Potential Energy Function *****************

      IF((NUA(ISTATE).GE.0).OR.(NUB(ISTATE).GT.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'adiabatic' BOB radial potential functions here ...
c        u_A(r) = yp*uA_\infty + [1 - yp]\sum_{i=0,NUA} {uA_i yq^i} 
c   where the  u_\infty  values stored/fitted as  UA(NUA(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REp= RE(ISTATE)**pAD(ISTATE)
          REq= RE(ISTATE)**qAD(ISTATE)
          HReP= 0.5d0*pAD(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qAD(ISTATE)/RE(ISTATE)
          IF(RREFadA(ISTATE).LE.0) AREFadA= RE(ISTATE)
          IF(RREFadB(ISTATE).LE.0) AREFadB= RE(ISTATE)
          IF(RREFadA(ISTATE).GT.0) AREFadA= RREFadA(ISTATE)
          IF(RREFadB(ISTATE).GT.0) AREFadB= RREFadB(ISTATE)
          AREFpADA= AREFadA**pAD(ISTATE)
          AREFqADA= AREFadA**qAD(ISTATE)
          AREFpADB= AREFadB**pAD(ISTATE)
          AREFqADB= AREFadB**qAD(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pAD(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
c ... reset parameter counter ...
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RDp= RVAL**pAD(ISTATE)
              RDq= RVAL**qAD(ISTATE)
              YPA= (RDp - AREFpADA)/(RDp + AREFpADA)
              YQA= (RDq - AREFqADA)/(RDq + AREFqADA)
              YPB= (RDp - AREFpADB)/(RDp + AREFpADB)
              YQB= (RDq - AREFqADB)/(RDq + AREFqADB)
              YPMA= 1.d0 - YPA
              YPMB= 1.d0 - YPB
              IF(BOBCN(ISTATE).GE.1) THEN
c** If  BOBCN > 0  &  p= 1,  assume use of Ogilvie-Tipping vble.
                  IF(pAD(ISTATE).EQ.1) THEN
                      YPA= 2.d0*YPA
                      YPB= 2.d0*YPB
                      ENDIF
                  ENDIF
              IF(NUA(ISTATE).GE.0) THEN
c ... Now ... derivatives of UA w.r.t. expansion coefficients
                  VAL= UA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPMA
                  YQPA= 1.d0
                  IF(NUA(ISTATE).GE.2) THEN
                      DO  J= 1,NUA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQPA*UA(J,ISTATE)
                          YQPA= YQPA*YQA
                          VAL= VAL+ UA(J,ISTATE)*YQPA
                          IPV= IPV+ 1
                          DVtot(IPV,I)= YPMA*YQPA
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPA
                  UAR(I,ISTATE)= VAL*YPMA + YPA*UA(NUA(ISTATE),ISTATE)
                  DUADRe(I,ISTATE)= 0.d0
c ... and derivative of UA w.r.t. Re ... if  RREFadA.LE.0
                  IF(RREFadA(ISTATE).LE.0.d0) 
     1               DUADRe(I,ISTATE)= -HReQ*(1.d0 - YQA**2)*YPMA*DVAL
     2            + HReP*(1.d0 - YPA**2)*(VAL- UA(NUA(ISTATE),ISTATE))
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
c ... Now ... derivatives of UB w.r.t. expansion coefficients
                  VAL= UB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NUA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YPMB
                  YQPB= 1.d0
                  IF(NUB(ISTATE).GE.2) THEN
                      DO  J= 1,NUB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQPB*UB(J,ISTATE)
                          YQPB= YQPB*YQB
                          VAL= VAL+ UB(J,ISTATE)*YQPB
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPMB*YQPB
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YPB
                  UBR(I,ISTATE)= VAL*YPMB + YPB*UB(NUB(ISTATE),ISTATE)
                  DUBDRe(I,ISTATE)= 0.d0
c ... and derivative of UB w.r.t. Re ... if  RREFadB.LE.0
                  IF(RREFadB(ISTATE).LE.0.d0) 
     1               DUBDRe(I,ISTATE)= -HReQ*(1.d0 - YQB**2)*YPMB*DVAL
     2            + HReP*(1.d0 - YPB**2)*(VAL- UB(NUB(ISTATE),ISTATE))
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of adiabatic potential BOB function++++++++++++++

      IF((NTA(ISTATE).GE.0).OR.(NTB(ISTATE).GE.0)) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any 'non-adiabatic' centrifugal BOB functions here ...
c       q_A(r) = yp*qA_\infty + [1 - yp]\sum_{i=0,NTA} {qA_i yq^i}
c   where the  q_\infty  values stored/fitted as  TA(NTA(ISTATE))
c  Incorporate the 1/r^2 factor into the partial derivatives (but not in
c     the g(r) functions themselves, since pre-SCHRQ takes care of that).
c    Need to add  M_A^{(1)}/M_A^{(\alpha)}  factor later too
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          REp= RE(ISTATE)**pNA(ISTATE)
          REq= RE(ISTATE)**qNA(ISTATE)
          HReP= 0.5d0*pNA(ISTATE)/RE(ISTATE)
          HReQ= 0.5d0*qNA(ISTATE)/RE(ISTATE)
          IF((BOBCN(ISTATE).GE.1).AND.(pNA(ISTATE).EQ.0)) THEN
              HReP= 2.d0*HReP
              HReQ= 2.d0*HReQ
              ENDIF
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RM2= 1/RVAL**2
              RDp= RVAL**pNA(ISTATE)
              RDq= RVAL**qNA(ISTATE)
              YP= (RDp - REp)/(RDp + REp)
              YQ= (RDq - REq)/(RDq + REq)
              YPM= 1.d0 - YP
              IF(BOBCN(ISTATE).GE.1) THEN
                  YPM= 1.d0
                  YP= 2.d0*YP
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN 
c ... Now ... derivatives of TA w,r,t, expansion coefficients
                  VAL= TA(0,ISTATE)
                  DVAL= 0.d0
                  IPV= IPVSTART + 1
                  DVtot(IPV,I)= YPM*RM2
                  YQP= 1.d0
                  IF(NTA(ISTATE).GE.2) THEN
                      DO  J= 1,NTA(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TA(J,ISTATE)
                          YQP= YQP*YQ
                          VAL= VAL+ TA(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TAR(I,ISTATE)= VAL*YPM + YP*TA(NTA(ISTATE),ISTATE)
c ... and derivative of TA w.r.t. Re ...
                  DTADRe(I,ISTATE)= (-HReQ*(1.d0 - YQ**2)*YPM*DVAL
     1        + HReP*(1.d0 - YP**2)*(VAL- TA(NTA(ISTATE),ISTATE)))*RM2
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
c ... Now ... derivatives of TB w,r,t, expansion coefficients
                  VAL= TB(0,ISTATE)
                  DVAL= 0.d0
                  IF(NTA(ISTATE).LT.0) THEN
                      IPV= IPVSTART + 1
                    ELSE
                      IPV= IPV + 1
                    ENDIF
                  DVtot(IPV,I)= YPM*RM2
                  YQP= 1.d0
                  IF(NTB(ISTATE).GE.2) THEN
                      DO  J= 1,NTB(ISTATE)-1
                          DVAL= DVAL+ DBLE(J)*YQP*TB(J,ISTATE)
                          YQP= YQP*YQ 
                          VAL= VAL+ TB(J,ISTATE)*YQP
                          IPV= IPV + 1
                          DVtot(IPV,I)= YPM*YQP*RM2
                          ENDDO
                      ENDIF
                  IPV= IPV + 1
                  DVtot(IPV,I)= YP*RM2
                  TBR(I,ISTATE)= VAL*YPM + YP*TB(NTB(ISTATE),ISTATE)
c ... and derivative of TA w.r.t. Re ...
                  DTBDRe(I,ISTATE)= (-HReQ*(1.d0 - YQ**2)*YPM*DVAL
     1        + HReP*(1.d0 - YP**2)*(VAL- TB(NTB(ISTATE),ISTATE)))*RM2
                  ENDIF
              ENDDO
          ENDIF
c++++ END of treatment of non-adiabatic centrifugal BOB function++++++++

      IF(NwCFT(ISTATE).GE.0) THEN
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Treat any Lambda- or 2\Sigma-doubling radial strength functions here
c    representing it as   f(r)= yp*w_inf + [1-yp]*Sum{ w_i * yp^i}
c    where the  w_inf  values stored/fitted as  wCFT(NwCFT(ISTATE))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          LAMB2= 2*IOMEG(ISTATE)
          REP= RE(ISTATE)**nPw(ISTATE)
          HReP= 0.5d0*nPw(ISTATE)/RE(ISTATE)
          IPVSTART= IPV
          DO  I= ISTART,ISTOP
              RVAL= RD(I,ISTATE)
              IF(RDIST.GT.0.d0) RVAL= RDIST
              RMF= 1.d0/RVAL**2
              IF(IOMEG(ISTATE).GT.0)  RMF= RMF**LAMB2
              RDp= RVAL**nPw(ISTATE)
              YP= (RDp - REP)/(RDp + REP)
              YPM= 1.d0 - YP
              DVAL= 0.d0
              YPP= YPM*RMF
              VAL= wCFT(0,ISTATE)*YPP
              IPV= IPVSTART + 1
              DVtot(IPV,I)= YPP
              IF(NwCFT(ISTATE).GE.2) THEN
                  DO  J= 1,NwCFT(ISTATE)-1
                      DVAL= DVAL+ DBLE(J)*YPP*wCFT(J,ISTATE)
                      YPP= YPP*YP
                      IPV= IPV + 1
                      DVtot(IPV,I)= YPP
                      VAL= VAL+ wCFT(J,ISTATE)*YPP
                      ENDDO
                  ENDIF
              IPV= IPV + 1
              DVtot(IPV,I)= YP*RMF
              wRAD(I,ISTATE)= VAL + YP*wCFT(NwCFT(ISTATE),ISTATE)*RMF
              dLDDRe(I,NSTATEMX)= -HReP*(1.d0 - YP**2)*(DVAL - VAL +
     1                                 wCFT(NwCFT(ISTATE),ISTATE)*RMF)
              ENDDO
          ENDIF
c++++ END of treatment of Lambda/2-sigma centrifugal BOB function+++++++

c++++ Test for inner wall inflection above the asymptote, and if it ++++
c++++ occurs, replace inward potential with linear approximation +++++++
      I1= (RE(ISTATE)-RD(1,ISTATE))/(RD(2,ISTATE)-RD(1,ISTATE))
      IF(I1.GT.3) THEN
          SL= 0.d0
          DO  I= I1-2, 1, -1
              SLB= SL
              SL= VPOT(I,ISTATE) - VPOT(I+1,ISTATE)
              IF((SL.LE.SLB).AND.(VPOT(I,ISTATE).GT.VLIM(ISTATE))) THEN
                  DO  J= I,1,-1
                      VPOT(J,ISTATE)= VPOT(I,ISTATE) + (I-J)*SL
                      ENDDO
                  WRITE(6,606) SLABL(ISTATE),RD(I,ISTATE),VPOT(I,ISTATE)
                  GOTO 66
                  ENDIF
              ENDDO
          ENDIF
   66 CONTINUE
  606 FORMAT(9('===')/'!*!* Extrapolate to correct ',A2,' inner-wall inf
     1lection at   R=',f6.4,'   V=',f8.0/9('==='))
c++++++++++++End of Inner Wall Test/Correction code+++++++++++++++++++++
c======================================================================
c** At the one distance RDIST calculate total effective potential VDIST
c  including (!!) centrifugal and Lambda/2Sigma doubling terms,
c  and its partial derivatives w.r.t. Hamiltonian parameters dVdPk.
c** This case only for simulation & fitting of tunneling width data.
c
      IF((RDIST.GT.0).AND.(IDAT.GT.0)) THEN 
          NBAND= IB(IDAT)
          IISTP= ISTP(NBAND)
cccccccc
c         WRITE (40,644) IISTP,RDIST,RVAL,VDIST,I,NDATPT(ISTATE)
c 644 FORMAT ('IISTP =',I3,' RDIST =',G16.8,' RVAL =',G16.8,
c    &          ' VDIST =',G16.8,' I =',I6,' NDATPT =',I6)
cccccccc
          BFCT= 16.85762920d0/(ZMASS(3,IISTP)*RDIST**2)
          JFCT= DBLE(JPP(IDAT)*(JPP(IDAT)+1))
          IF(IOMEG(ISTATE).GT.0) JFCT= JFCT - IOMEG(ISTATE)**2
          IF(IOMEG(ISTATE).EQ.-2) JFCT= JFCT + 2.D0
          JFCT= JFCT*BFCT
c ... First get total effective potential, including BOB terms
          VDIST= VDIST + JFCT
          IF(NUA(ISTATE).GE.0) VDIST= VDIST 
     1                          + ZMUA(IISTP,ISTATE)*UAR(ISTOP,ISTATE)
          IF(NUB(ISTATE).GE.0) VDIST= VDIST 
     1                          + ZMUB(IISTP,ISTATE)*UBR(ISTOP,ISTATE)
          IF(NTA(ISTATE).GE.0) VDIST= VDIST 
     1                     + JFCT*ZMTA(IISTP,ISTATE)*TAR(ISTOP,ISTATE)
          IF(NTB(ISTATE).GE.0) VDIST= VDIST 
     1                     + JFCT*ZMTB(IISTP,ISTATE)*TBR(ISTOP,ISTATE)
          JFCTLD= 0.d0
          IF(IOMEG(ISTATE).NE.0) THEN
              IF(IOMEG(ISTATE).GT.0) THEN
c ... for Lambda doubling case ...
                  JFCTLD= (EFPP(IDAT)-EFREF(ISTATE)) 
     1         *(DBLE(JPP(IDAT)*(JPP(IDAT)+1))*BFCT**2)**IOMEG(ISTATE)
                  ENDIF
              IF(IOMEG(ISTATE).EQ.-1) THEN
c ... for doublet Sigma doubling case ...
                  IF(EFPP(IDAT).GT.0) JFCTLD= 0.5d0*JPP(IDAT)*BFCT
                  IF(EFPP(IDAT).EQ.0) JFCTLD= 0.d0
                  IF(EFPP(IDAT).LT.0) JFCTLD= -0.5d0*(JPP(IDAT)+1)*BFCT
                  ENDIF
              VDIST= VDIST + JFCTLD* WRAD(ISTOP,ISTATE)
             ENDIF
cccccccc
c       WRITE (40,648) JPP(IDAT),EFPP(IDAT),RDIST,VDIST
c 648   FORMAT ('J =',I3,' efPARITY =',I3,' RDIST =',G16.8,' VDIST =',
c    1                                                           G16.8/)
cccccccc
          IF(PSEL(ISTATE).GT.0) THEN
              DO  IPV= 1,TOTPOTPAR
                  dVdPk(IPV)= 0.d0
                  ENDDO
c** Now ... generate requisite partial derivatives.
              DO IPV= POTPARI(ISTATE),POTPARF(ISTATE)
                  dVdPk(IPV)= DVtot(IPV,ISTOP)
                  ENDDO
              IF(NUA(ISTATE).GE.0) THEN
                  DO IPV= UAPARI(ISTATE),UAPARF(ISTATE)
                      dVdPk(IPV)= ZMUA(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
              IF(NUB(ISTATE).GE.0) THEN
                  DO IPV= UBPARI(ISTATE),UBPARF(ISTATE)
                      dVdPk(IPV)= ZMUB(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
              IF(NTA(ISTATE).GE.0) THEN
                  DO IPV= TAPARI(ISTATE),TAPARF(ISTATE)
                     dVdPk(IPV)=JFCT*ZMTA(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                     ENDDO
                  ENDIF
              IF(NTB(ISTATE).GE.0) THEN
                  DO IPV= TBPARI(ISTATE),TBPARF(ISTATE)
                     dVdPk(IPV)=JFCT*ZMTB(IISTP,ISTATE)*DVtot(IPV,ISTOP)
                     ENDDO
                  ENDIF
              IF(NwCFT(ISTATE).GE.0) THEN
                  DO IPV= LDPARI(ISTATE),LDPARF(ISTATE)
                      dVdPk(IPV)= JFCTLD*DVtot(IPV,ISTOP)
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
c*****7********************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
