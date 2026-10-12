c***********************************************************************
      SUBROUTINE DVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.
c         *****  Version left by  Alek  31 May 2010 ****
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
       INTEGER check(0:NPARMX+3),j,counter,nParams,ISTATE,ISTART,
     1 IDAT,i,m,k,kk,jj
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,BTOT(3),
     1 dBTOTdp(NPARMX),CONST(3),jump,b,YVAL,RVAL,
     2 VDIST,EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,Bq2,
     3 INTEGRALS(NPARMX+3),error,BVIR,BETADISTa,dVdR,d2VdR,XTEMP1,
     4 dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU
c
       REAL*8, PARAMETER :: k_boltz=1.3806504D-23, NA=6.02214179D23,
     1  h = 6.62606896D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538782D-27,
     3  Cu= 16.857629D0
c* Common block data for quadrature weights and points
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
c***********************************************************************
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO
      DO J= 1,HPARMX
          dBVIRdP(J)= 0.d0
          ENDDO
      BTemp = k_cm*Temp(IDAT)
c     BTEMP= DE(ISTATE)*5.d0
c     PRINT*, 'temp', BTemp
      BTempInv = 1.d0/BTemp
      BTOT = 0.d0
      dBTOTdP = 0.d0
      Const(1)= -2*pi*0.602214179d0
      Const(2)= Cu*pi*BTempInv**3/(6.d0*ZMU)
      Const(3)= -(pi/6.d0)*(Cu*BTempInv**2/ZMU)**2
      nParams = 3 + NCMM(ISTATE) + NLpow(ISTATE)
c.. this first loop is here to repeat the calculations until the values
c.  converge to within error, or until they are repeated 7 times without
c.  convergence
c     DO k= 1,2
        check = 1
        DO  kk= 0,7
          IF(check(0).EQ.-1) EXIT
          counter= 2.d0**kk
          INTEGRALS= 0.d0
          jump= 2.d0/DBLE(counter)
          b= -1.d0 + jump
c. the first loop is over each Gaussian point
          DO i= 1,8
              YVAL= (b + 1.d0)*XG(i)/2.d0 + (b - 1.d0)/2
c. this next loop is over each subinterval within (-1,1), so that every
c. Gaussian point can be evaluated over each subinterval
              DO m= 1,counter
c                 IF(k.EQ.1) RVAL= Re(ISTATE)*(1.1 + 0.9*YVAL)
c    1                                                  /(1.d0 - YVAL)
c                 IF(k.EQ.2) RVAL= Re(ISTATE)*DSQRT((1.0001d0 
c    1                                 + 0.9999d0*YVAL)/(1.d0 - YVAL))
                  RVAL= Re(ISTATE)*DSQRT((1.0001d0 + 0.9999d0*YVAL)
     1                                                 /(1.d0 - YVAL))
                  VDIST= 0.d0
                  CALL VGENp(ISTATE,RVAL,VDIST,BETADISTa,IDAT,dVdR,
     1                                                          d2VdR)
c. some of the terms that will be used in later calculations of the
c. virial coefficients are constructed here
                  EXP_TERM= DEXP(-VDIST*BTempINV)
                  RINV= 1.d0/RVAL
c                 IF(k.EQ.1) int_fact= Re(ISTATE)*(b + 1.d0)
c    1                                               /(1.d0 - YVAL)**2
c                 IF(k.EQ.2) int_fact= (Re(ISTATE)/(1.d0 - YVAL))**2
c    1                                          *RINV*(b + 1.d0)*0.5d0
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL))**2
     1                                          *RINV*(b + 1.d0)*0.5d0
                  Vsq= VDIST**2
                  VPsq= dVdR**2
                  Rsq= RVAL**2
                  XTEMP= d2VdR**2/10.d0 + VPsq*RINV**2/5.d0
     1         + dVdR**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients are calculated
                  DO J = 1,nParams
c                     jj = ISTART + J - 1
                      jj = ISTART + J
                      IF(check(J).eq.1) THEN
c* the derivative of the classical expression
                          dBcdP(J)= -BTempInv*EXP_TERM
     1                                  *DVtot(jj,NPNTMX)*Rsq*int_fact
c* and of the first quantum correction 
                          XTEMP1 = -VPsq*BTempInv*DVtot(jj,NPNTMX)
     1                                       + 2*dVdR*dVpdP(jj,NPNTMX)
                          dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
  
c..  As the final step these terms all added together in a weighted sum
  
                          INTEGRALS(J) = INTEGRALS(J) + Const(1)
     1                      *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
c                 IF(b.EQ.-0.875d0) WRITE(70+J,664) (YVAL, dBcdP(J),
c    1                                                      dBq1dP(J))
                          ENDIF
                      ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. now the integrands are evaluated at each particular Gaussian point
c. and summed together with the proper weighting
c                 Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
c                 IF(k.EQ.1) Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
c                 IF(K.EQ.2) THEN
c                 Bq1= EXP_TERM*VPsq*Rsq*int_fact
c                 Bq2= EXP_TERM*XTEMP*Rsq*int_fact
c                   ENDIF
c                 IF(k.EQ.1) THEN
                    IF(check(nParams+1).EQ.1) THEN
                      Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
                      INTEGRALS(nParams+1)= INTEGRALS(NParams+1) 
     1                                                 + Bclass*WG(i)
                      ENDIF
c                   ENDIF
c                 IF(k.EQ.2) THEN
                    IF(check(nParams+2).EQ.1) THEN
                      Bq1= EXP_TERM*VPsq*Rsq*int_fact
                      INTEGRALS(nParams+2)= INTEGRALS(nParams+2) 
     1                                                    + Bq1*WG(i)
                      ENDIF
                    IF(check(nParams+3).EQ.1) THEN
                      Bq2= EXP_TERM*XTEMP*Rsq*int_fact
                      INTEGRALS(nParams+3)= INTEGRALS(nParams+3)
     1                                                    + Bq2*WG(i)
                      ENDIF
c                   ENDIF
c. this write statement is to examine the behaviour of the integrands
c                 IF(b.EQ.-0.875d0) THEN
c                      IF (k.EQ.1) WRITE(46,661) (YVAL, TEMP(IDAT),
c    1                                              Bclass, Bq1, Bq2)
c                      IF (k.EQ.2) WRITE(45,661) (YVAL, TEMP(IDAT),
c    1                                              Bclass, Bq1, Bq2)
c                      WRITE(46,661) (YVAL, TEMP(IDAT),
c    1                                              Bclass, Bq1, Bq2)
c                      ENDIF
                  YVAL= YVAL + jump
                  ENDDO
              ENDDO
c. this section checks each integral to see if it has converged,
c. assigning a check value of -1 to those that have, which will stop
c. further individual calculations for them
          IF(counter.GE.2) THEN
              error= 0.000001
              DO i= 1,nParams
                  IF(check(i).EQ.1) THEN
                      IF(abs(INTEGRALS(i) - dBTOTdp(i)).LT.(abs(error
     1                                  *INTEGRALS(i)))) THEN
                         check(i)= -1
c                        IF(k.EQ.1) WRITE(61,662) (i, counter,
c    1                                                     Temp(IDAT))
c                        WRITE(60,662) (i, counter, Temp(IDAT))
c                        IF(k.EQ.2) WRITE(60,662) (i, counter,
c    1                                                     Temp(IDAT))
                         ENDIF
                      dBTOTdp(i)= INTEGRALS(i)
                      INTEGRALS(i)= 0.d0
                      ENDIF
                  ENDDO
              DO i= 1,3
                  IF(check(nParams+i).EQ.1) THEN
                     IF(abs(INTEGRALS(nParams+i)-BTOT(i)).LT.(abs(error
     1                                   *INTEGRALS(nParams+i)))) THEN
                         check(nParams+i)= -1
c                        IF(k.EQ.1) WRITE(61,662) (i, counter,
c    1                                                     Temp(IDAT))
c                        WRITE(61, 662) (i, counter, Temp(IDAT))
c                        IF(k.EQ.2) WRITE(60,662) (i, counter,
c    1                                                     Temp(IDAT))
                         ENDIF
                      BTOT(i)= INTEGRALS(nParams+i)
                      INTEGRALS(nParams+i)= 0.d0
                      ENDIF
                  ENDDO
c. this next section checks if all the integrals have converged. If they
c. have, then check(0) is set to -1 to stop the calculations there
              check(0)= -1
              DO i= 1,nParams+3
                  IF(check(i).EQ.1) THEN
                      check(0)= 1
                      EXIT
                      ENDIF
                  ENDDO
              ELSE
              DO i= 1,3
                  BTOT(i)= INTEGRALS(nParams+i)
                  ENDDO
              DO j= 1,nParams
                  dBTOTdp(j)= INTEGRALS(j)
                  ENDDO
              ENDIF
c         IF(k.EQ.1) WRITE(63,663) (counter, b, Temp(IDAT), BTOT(1:3))
c         IF(k.EQ.2) WRITE(62,663) (counter, b, Temp(IDAT), BTOT(1:3))
          WRITE(62,663) counter, Temp(IDAT), (BTOT(i),i=1,3)
          WRITE(63,663) counter, Temp(IDAT), (dBTOTdP(i),i=1,3)
          ENDDO
c       PRINT *, dBTOTdp(1)
c.. now that the integrals have been evaluated, they can be summed
c.  together to form the second virial coefficient, B
c       ENDDO
      BVIR= Const(1)*BTOT(1) + Const(2)*BTOT(2) + Const(3)*BTOT(3)
      DO J= 1,nParams
          dBVIRdP(J)= dBTOTdp(J)
          ENDDO
c. this is just for testing
      WRITE (44,660) (TEMP(IDAT),BVIR)
  660 FORMAT(G12.5, G18.10, G18.10)
  661 FORMAT(G18.10, F6.0, G18.10, G18.10, G18.10)
  662 FORMAT(I4, I4, ' ',  F5.0)
  663 FORMAT(I4, F5.0, 1PD18.10, 1PD18.10, 1PD18.10)
  664 FORMAT(G18.10, G18.10, G18.10)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12


c***********************************************************************
      SUBROUTINE VGENP(ISTATE,RDIST,VDIST,BETADIST,IDAT,dVdR,d2VdR)
c***********************************************************************
c** This subroutine will generate one of six possible families of
c   analytical molecular potentials for the direct hamiltonian fitting
c   program. This family is:
c     - The Modified Lennard-Jones Oscillator
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++  COPYRIGHT 2006-2009  by  R.J. Le Roy, Jenning Seto & Yiye Huang ++
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c           ----- Version of 22 February 2010 -----
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c   ISTATE  is the electronic state being considered in this CALL.
c   RDIST: If RDIST > 0, calculate potl & derivs only @ that onee distance
c   -----    * return potential function at that point as VDIST and
c              potential function exponent as BETADIST
c            * skip partial derivative calculation if  IDAT.le.0
c        * If RDIST.le.0  calculate partial derivatives at distances
c             given by array RD(i,ISTATE) & return them in array DVtot
c** On entry via common blocks:
c  NSpow(s)  is the order of the  beta(r)  exponent expansion for R.le.Re
c  NLpow(s)  is the order of the  beta(r)  exponent expansion for R > Re
c!!!  The code (currently) assumes  NLpow .ge. NSpow !!!!
c    MMLR(j,s)  are long-range inverse-powers for an MLR or DELR potential
c    nPB(s)  the basic value of power p for the beta(r)  exponent function
c    nQB(s)  the power p for the power series expansion variable in beta(r)
c    pAD(s) & qAD(s) the values of power p for adiabatic u(r) BOB functions
c    nNA(s) & qNA(s) the values of power p for centrifugal q(r) BOB functions
c    nPw(s)  the value of power m for Lambda(r) function
c    DE     is the Dissociation Energy for each state.
c    RE     is the Equilibrium Distance for each state.
c    BETA    is the array of potential (exponent) expansion parameters 
c    NDATPT is the number of meshpoints used for the array.
c-----------------------------------------------------------------------
c** On exit via common blocks:
c    R      is the distance array
c    VPOT   is the potential that is generated.
c    BETAFX  is used to contain the  beta(r)  function.
c** Internal partial derivative arrays ...
c    DUADRe & DUBDRe are p.derivs of adiabatic fx. w.r.t.  Re
c    DVDQA & DVDQB are p.derivs of non-adiabatic fx. wrt q_A(i) & q_B(i)
c    DTADRe & DTBDRe are p.derivs of non-adiabatic fx. w.r.t.  Re
c    dVdL & dLDDRe are p.derivatives of f_\lambda(r) w.r.t. beta_i & Re
c    DBDB & DBDRe are p.derives of beta(r) w.r.t. \beta_i & Re, respectively
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
     1  IDAT, NBAND, IISTP,MMLR1D(NCMMax),KDER
      REAL*8 BTEMP,BINF,RVAL,RTEMP,RM2,XTEMP,PBTEMP,PETEMP,RET,
     1 FSW,Xtemp2,Btemp2,BMtemp,BMtemp2,RMF,PBtemp2,C3VAL,C3bar,C6bar,
     2 C6adj,C9adj,YP,YQ,YPA,YPB,YQA,YQB,YPE,YPM,YPMA,YPMB,YPP,YQP,YQPA,
     3 YQPB,REp,Req,RDp,RDq,DYPDRE,DYQDRE,VAL,DVAL,HReP,HReQ,SL,SLB,
     4 AREF,AREFp,AREFq, RE3,RE6,RE8,T0,T0P,T1,ULRe,Scalc,dLULRedCm(9),
     5 dLULRedRe,dLULRedDe,dULRdDe,dULRdCm(9),RD3,RD6,RD8,DVDD,RDIST,
     6 VDIST,BETADIST,BFCT,JFCT,JFCTLD,RETSig,RETPi,RETp,RETm,AREFadA,
     7 AREFadB,AREFpADA,AREFpADB,AREFqADA,AREFqADB,D2VAL,dYPdR,dYPEdR,
     8 dYQdR,d2YPdR,d2YQdR,d2YPEdR,RINV,dDULRdR,d2DULRdR,dULRdR,d2ULRdR,
     9 DXTEMP,D2XTEMP,dVdR,d2VdR,dLULRdR,YPPP,dBdR,d2BdR,VVAL,YVAL,
     x dULRdRCm(9),dXdP(HPARMX),dXpdP(HPARMX),dLULRdCm(9),DYPEDRE,
     y dVALdRe,dYBdRe,dBpdRe,DYPpDRE,DYPEpdRE,DYQpDRE,dYBpdRe,
     z xBETA(NbetaMX),rKL(NbetaMX,NbetaMX)
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
      REP= RE(ISTATE)**nPB(ISTATE)
      IF(RREF(ISTATE).LE.0) AREF= RE(ISTATE)
      IF(RREF(ISTATE).GT.0) AREF= RREF(ISTATE)
      AREFp= AREF**nPB(ISTATE)
      AREFq= AREF**nQB(ISTATE)
c** Normally data point starts from 1
      ISTART= 1
      ISTOP= NDATPt(istate)
c** When calculating only one potential point
      IF(RDIST.GT.0) THEN
          ISTART= NPNTMX
          ISTOP= NPNTMX
          VDIST= 0.0d0
          BETADIST= 0.d0
          ENDIF
      PBTEMP= 0.0d0
      PETEMP= 0.0d0
      DO  I= ISTART,ISTOP
          BETAFX(I,ISTATE)= 0.0d0
          UAR(I,ISTATE)= 0.d0
          UBR(I,ISTATE)= 0.d0
          TAR(I,ISTATE)= 0.d0
          TBR(I,ISTATE)= 0.d0
          UAR(I,ISTATE)= 0.d0
          WRAD(I,ISTATE)= 0.d0
          ENDDO
c** Initialize parameter counter for this state ...
      IPVSTART= POTPARI(ISTATE) - 1
c** First - define values & derivatives of uLR at Re fro MLR potential
      ULRe= 0.d0
      T1= 0.d0
      IF(rhoAB(ISTATE).GT.0.d0) THEN
c ... save uLR powers in a 1D array
          DO  m= 1, NCMM(ISTATE)
              MMLR1D(m)= MMLR(m,ISTATE)
              ENDDO
          KDER= 1 
          CALL dampF(RE(ISTATE),rhoAB(ISTATE),NCMM(ISTATE),
     1      MMLR1D,IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
          ENDIF
      DO  m= 1,NCMM(ISTATE)
          IF(rhoAB(ISTATE).LE.0.d0) THEN
              dLULRedCm(m)= 1.d0/RE(ISTATE)**MMLR(m,ISTATE)
            ELSE
              dLULRedCm(m)= Dm(m)/RE(ISTATE)**MMLR(m,ISTATE)
            ENDIF
          T0= CmVAL(m,ISTATE)*dLULRedCm(m)
          ULRe= ULRe + T0
          T1= T1 + MMLR(m,ISTATE)*T0
          ENDDO
      dLULRedRe= -T1/(ULRe*RE(ISTATE))
      DO  m= 1,NCMM(ISTATE)
          dLULRedCm(m)= dLULRedCm(m)/ULRe
          IF(rhoAB(ISTATE).GT.0) THEN
              dLULRedRe= dLULRedRe + dLULRedCm(m)*Dmp(m)/Dm(m)
              ENDIF
          ENDDO
      BINF= DLOG(2.0d0*DE(ISTATE)/ULRe)
      betaINF(ISTATE)= BINF
      DO  I= ISTART,ISTOP
          RVAL= RD(I,ISTATE)
          IF(RDIST.GT.0.d0) RVAL= RDIST
          RINV= 1.d0/RVAL
          RDp= RVAL**nPB(ISTATE)
          RDq= RVAL**nQB(ISTATE)
          YPE= (RDp-REP)/(RDp+REP)
          YP= (RDp-AREFp)/(RDp+AREFp)
          YQ= (RDq-AREFq)/(RDq+AREFq)
          YPM= 1.d0 - YP
          DYPDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YP**2)/RE(ISTATE)
          DYQDRE= -0.5d0*nQB(ISTATE)*(1.d0 - YQ**2)/RE(ISTATE)
c----------------------------------------------------------------------
          DYPEDRE= -0.5d0*nPB(ISTATE)*(1.d0 - YPE**2)/RE(ISTATE)
          DYPDR= -DYPDRE*RE(ISTATE)*RINV
          DYPEDR= 0.5d0*nPB(ISTATE)*RINV*(1.d0 - YPE**2)
          DYQDR= -DYQDRE*RE(ISTATE)*RINV
          D2YPDR= -DYPDR*RINV*(1.d0 + nPB(ISTATE)*YP)
          D2YPEDR= -DYPEDR*RINV*(1.d0 + nPB(ISTATE)*YPE)
          D2YQDR= -DYQDR*RINV*(1.d0 + nQB(ISTATE)*YQ)
          DYPpDRE= -nPB(ISTATE)*YP*RINV*DYPDRE
          DYPEpDRE= -nPB(ISTATE)*YPE*RINV*DYPEDRE
          DYQpDRE= -nQB(ISTATE)*YQ*RINV*DYQDRE
          D2VAL= 0.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          YPP= 1.d0
          DVAL= 0.d0
          DBDB(0,I,ISTATE)= 1.0d0
c-------------------------------------------------------------------
          VAL= BETA(0,ISTATE) + YQ*BETA(1,ISTATE)
          DVAL= BETA(1,ISTATE)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          npow= NSpow(ISTATE)
          IF(RVAL.GT.RE(ISTATE)) npow= NLpow(ISTATE)
c-------------------------------------------------------------------
          DO  J= 2,npow
c... now calculate power series part of the Morse-like exponent,along
c    with its radial derivatives
              D2VAL= D2VAL + BETA(J,ISTATE)* DBLE(J)
     1                                            *DBLE(J - 1) *YPP
              YPP= YPP*YQ
              DVAL= DVAL + BETA(J,ISTATE)* DBLE(J)* YPP
              YPPP= YPP* YQ
              VAL= VAL + BETA(J,ISTATE)*YPPP
              DBDB(J,I,ISTATE)= YPM*YPPP
              ENDDO
          YPP= YPPP
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*** DBDB & DBDRe= dBeta/dRe  used in uncertainty calculation in fununc.f
          DBDRe(I,ISTATE)= -YP*dLULRedRe
          dVALdRe= DBDRe(I,ISTATE) + (BINF - VAL)*DYPDRE
     1                                   + (1.d0 - YP)*DVAL*DYQDRE
          IF(RREF(ISTATE).LE.0.d0) DBDRe(I,ISTATE)= dVALdRe
c-----------------------------------------------------------------------
c... now the power series and its radial derivatives are used in the
c    construction of the derivatives  with respect to the parameters
          dBpdRe= DYPpDRE*(BINF - VAL) - DYPDR*dLULRedRe
     1   + (-DYPDR*DYQDRE + (1.d0 - YP)*DYQpDRE - DYPDRE*DYQDR)*DVAL
     2                              + (1.d0 - YP)*DYQDR*DYQDRE*D2VAL
          D2VAL= (BINF - VAL)*D2YPDR - 2.d0*DYPDR*DYQDR*DVAL
     1                   + (1.d0- YP)*(D2YQDR*DVAL + DYQDR**2*D2VAL)
          DVAL= (BINF - VAL)*DYPDR + (1.d0- YP)*DYQDR*DVAL
          VAL= YP*BINF + (1.d0- YP)*VAL
          dBdR= dYPEdR*VAL + YPE*DVAL
          d2BdR= d2YPEdR*VAL + 2.d0*dYPEdR*DVAL + YPE*D2VAL
          dYBdRe= DYPEDRE*VAL + YPE*dVALdRe
          dYBpdRe= VAL*DYPEpDRE + DYPEDRE*DVAL + DYPEDR*dVALdRe
     1                                                  + YPE*dBpdRe
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          BETAFX(I,ISTATE)= VAL
          XTEMP= DEXP(-VAL*YPE)
c** Now begin by generating  uLR(r)
          ULR= 0.d0
c-------------------------------------------------------------------
          dULRdR= 0.d0
          d2ULRdR= 0.d0
          dULRdRCm= 0.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          IF(rhoAB(ISTATE).GT.0.d0) THEN
              KDER= 2
              CALL dampF(RVAL,rhoAB(ISTATE),NCMM(ISTATE),MMLR1D,
     1                     IDF(ISTATE),IDSTT(ISTATE),KDER,Dm,Dmp,Dmpp)
              ENDIF
          DO  m= 1,NCMM(ISTATE)
              IF(rhoAB(ISTATE).LE.0.d0) THEN
c-----------------------------------------------------------------------
                  dULRdCm(m)= 1.d0*RINV**MMLR(m,ISTATE)
                  dULRdRCm(m)= -dULRdCm(m)*RINV*DBLE(MMLR(m,ISTATE))
                  dDULRdR= 0.d0
                  d2DULRdR= 0.d0
                ELSE
                  dULRdCm(m)= Dm(m)*RINV**MMLR(m,ISTATE)
                  dULRdRCm(m)= -dULRdCm(m)*RINV *DBLE(MMLR(m,ISTATE)) 
     2                                   + Dmp(m)*RINV**MMLR(m,ISTATE)
                  dDULRdR= Dmp(m)*RINV**MMLR(m,ISTATE)
                  d2DULRdR= Dmpp(m)*RINV**MMLR(m,ISTATE)
                ENDIF
              ULR= ULR + CmVAL(m,ISTATE)*dULRdCm(m)
              dULRdR= dULRdR + CmVAL(m,ISTATE)*(dDULRdR
     1                         - dULRdCm(m)*RINV*DBLE(MMLR(m,ISTATE)))
              d2ULRdR= d2ULRdR + CmVAL(m,ISTATE)*(d2DULRdR
     1   - 2.d0*dDULRdR*RINV*DBLE(MMLR(m,ISTATE)) + dULRdCm(m)*RINV**2
     2               *DBLE(MMLR(m,ISTATE))*DBLE((MMLR(m,ISTATE) + 1)))
              ENDDO
          dLULRdR= dULRdR/ULR
          DO m= 1,NCMM(ISTATE)
              dLULRdCm(m)= dULRdCm(m)/ULR
              ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          XTEMP= XTEMP*ULR/ULRe
c... note ... reference energy for each state is asymptote ...
          DVDD= XTEMP*(XTEMP - 2.D0)  
          VPOT(I,ISTATE)= DE(ISTATE)*DVDD + VLIM(ISTATE)
          IF(RDIST.GT.0.d0) THEN
              VDIST= VPOT(I,ISTATE)
              BETADIST= VAL
              IF(IDAT.LE.0) GO TO 999
              ENDIF
c-              ENDDO
          YPP= 2.d0*DE(ISTATE)*(1.0d0-XTEMP)*XTEMP
          IPV= IPVSTART+2
c-----------------------------------------------------------------------
c... derivatives w.r.t R
          YVAL= (RVAL/RE(ISTATE) - 1.1d0)/(RVAL/RE(ISTATE) + 0.9d0)
          DXTEMP= XTEMP*(dLULRdR - dBdR)
          D2XTEMP= XTEMP*(dBdR**2 - d2BdR + (d2ULRdR 
     1                                           - 2*dBdR*dULRdR)/ULR)
          VVAL= DE(ISTATE)*(1.d0- XTEMP)**2 - DE(ISTATE)
          dVdR= 2.d0*DE(ISTATE)*DXTEMP*(XTEMP - 1.d0)
          d2VdR= 2.d0*DE(ISTATE)*(DXTEMP**2 + D2XTEMP
     1                                                *(XTEMP - 1.d0))
c *** This is just to write the derivatives for testing
c             IF(RDIST.LT.0) WRITE (40,640) (RVAL,VVAL,dVdR,d2VdR,
c    1                                                           YVAL)
c 640         FORMAT(G12.5, G18.10, G18.10, G18.10, G14.7)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c ... derivative w.r.t. Cm's
          DO  m= 1, NCMM(ISTATE)
              IPV= IPV+ 1
              dXdP(IPV)= XTEMP*(dLULRdCm(m) + (YPE*YP - 1.d0)
     1                                                  *dLULRedCm(m))
              dXpdP(IPV)= DEXP(-VAL*YPE)/ULRe*(dULRdRCm(m)
     1                    - dBdR*dULRdCm(m)) + (DXTEMP*(YPE*YP - 1.d0)
     2                   + XTEMP*(dYPEdR*YP + YPE*dYPdR))*dLULRedCm(m)
              dVpdP(IPV,I)= 2.d0*DE(ISTATE)*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
              DVtot(IPV,I)= -YPP*(dLULRedCm(m)*(YP*YPE- 1.d0)
     1                                               + dULRdCm(m)/ULR)
              ENDDO
c... derivative w.r.t. Re  
          dXdP(IPVSTART+2)= -XTEMP*(dYBdRe + dLULRedRe)
          dXpdP(IPVSTART+2)= -DXTEMP*(dYBdRe + dLULRedRe)
     1                                                 - XTEMP*dYBpdRe
          dVpdP(IPVSTART+2,I)= 2.d0*DE(ISTATE)*(dXdP(IPVSTART+2)
     1                     *DXTEMP + (XTEMP - 1.d0)*dXpdP(IPVSTART+2))
          DVtot(IPVSTART+2,I)= YPP*(dYBdRe + dLULRedRe)
c... derivative w.r.t. De
          dXdP(IPVSTART+1)= -XTEMP*YPE*YQ
          dXpdP(IPVSTART+1)= -(XTEMP*(YPE*DYQDR + DYPEDR*YQ)
     1                                            + YPE*YQ*DXTEMP)
          DVDD= DVDD + YPP*YP*YPE/DE(ISTATE)
          YPP= YPP*YPE*(1.d0 - YP)
          dVpdP(IPVSTART+1,I)= 2.d0*(dXdP(IPVSTART+1)*DXTEMP
     1                             + (XTEMP - 1.d0)*dXpdP(IPVSTART+1)) 
     2                                    + 2.d0*(XTEMP - 1.d0)*DXTEMP
          DVtot(IPVSTART+1,I)= DVDD
c... finally ... derivatives w.r.t. exponent expansion coefficients
          DO  J= 0,npow
              IPV= IPV+1
              dXdP(IPV)= XTEMP*YPE*(1.d0 - YP)*YQ**J
              dXpdP(IPV)= (XTEMP*((1.d0 - YP)*DYPEDR - DYPDR*YQ)
     1                  + YPE*(1.d0 - YP)*DXTEMP)*YQ**J + XTEMP*J*(YPE
     2                                       *(1.d0 - YP))*YQ**(J - 1)
              dVpdP(IPV,I)= 2.d0*DE(ISTATE)*(dXdP(IPV)*DXTEMP
     1                                    + (XTEMP - 1.d0)*dXpdP(IPV))
              DVtot(IPV,I)= YPP
              YPP= YPP*YQ
              ENDDO
          ENDDO
ccccc Print for testing
      rewind(10)
      write(10,610) (RD(i,ISTATE),vpot(i,istate),BETAFX(i,istate),
     1                                        i= 1, NDATPT(ISTATE),50)
  610 FORMAT(/(f10.4,f15.5,f12.6))
ccccc End of Print for testing

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
c*****7********************** BLOCK END ******************************72
  999 RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
