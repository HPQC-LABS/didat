c*******************************************************************
      SUBROUTINE PHSHTD(L,ZE,ADLT,STD,NODX,N,V,XM2,BZ,RH,RMIN,
     1  PHTST,TDTST,IWR,LIFE)
c**** R.J.LE ROY  SUBROUTINE 'PHSHTD' , AS OF 16 NOVEMBER, 1978  ***
c** SUBROUTINE TO CALCULATE THE ABSOLUTE PHASE SHIFT  'ADLT' , AND
c  THE TIME-DELAY 'STD' FOR POTENTIAL 'V(I)' (INCLUDING FACTOR 'BZ')
c  AT ENERGY 'ZE' AND PARTIAL WAVE 'L'.  SOLVE THE RADIAL SCHROEDINGER
c  EQUATION IN THE FORM    D2S/DR2= -(ZE*BZ-V(I)-L*(L+1)*XM2)*S ,
c  WHERE XM2= 1/X**2 AND S IS THE WAVE FUNCTION, USING NUMEROV-S METHOD.
c** INTEGRATE OUTWARD FROM  'RMIN' WITH AN INCREMENT OF  'RH' (BOTH
c  IN UNITS OF  ANGSTROMS ,  FOR UP TO 'N' STEPS.
c** 'NODX' IS NO. OF WAVE FX. NODES BEHIND THE POTENTIAL BARRIER MAX.
c** IF(IWR.GT.0) PRINTS ALL TRIAL PHASE SHIFTS AND NODE POSITIONS.
c  IF(IWR.NE.0) PRINT CONVERGED PHASE SHIFTS AND NODE POSITIONS
c** IF(LIFE.GT.0) COMPUTE COLLISIONAL TIME DELAY. ** IN THIS CASE
c  RADIAL INTEGRATION STOPS WHEN TIME DELAYS AT 3 SUCCESSIVE NODES
c  HAVE RELATIVE DIFFERENCES SMALLER THAN 'TDTST'.
c** IF(LIFE.LE.0) COMPUTE ONLY PHASE SHIFT. ** HERE RADIAL INTEGRATION
c  FINISHED WHEN CORRECTED VALUES AT 3 SUCCESSIVE NODES AGREE TO
c  WITHIN 'PHTST' (RADIANS).
c** CALLS SUBROUTINES 'PREPOT', 'BSELFX' AND 'CORTD'.
c-------------------------------------------------------------------
      INTEGER I,II,IAN1,IAN2,IFIRST,IG,IJ,IK,ILK1,IM,IMN1,IMN2,ISG,IST,
     1        IWR,IQ,ISTATE,IDAT,J,JN,J1,K,L,LD,LE,LH,LI,LIFE,L4,FE,LMX,
     2        LNPT,LPREV,M,N,NCN,NGT,NNL(66),NODE,NODX
      REAL*8 ADLT,ANR2,ARG,BETADIST,BZ,COR1,COR2,COR3,COTA,DLT,DLTX,DLTY
     1  ,DLT1,DLT2,DLT3,DSA1,DSA2,DSI1,DSI2,E,EJ,EMV,GI,GN,GMA,H,H2,ILK,
     2  NDT,PHTST,PI,PSC,RATIN,RG2,RH,RMIN,RMINN,RRM2,R4,SC,SCOR,SB,SBJ,
     3  SBN,SCRT,SGB,SGI,SI,SIK,SIR4,SM,SN,SNRM,SP,STD1,SQE,SQEMV,SRTGI,
     4  SRTGN,STD,STD2,TDTST,TH,X12,XJC,XJN,XLK,XNDT,XND1,XND2,XND3,
     5  XNORM,XNRM,VB,VI,VN,Y1,Y2,Y3,ZE,ZU,V(N),XM2(N),XX(8),WW(8),
     6  VV(8),RR(8),RM2(8)

      COMMON /GQP/XX,WW,NGT

      DATA PI/3.141592653589793D0/,IAN1/-1/,IAN2/-1/,IMN1/-1/,IMN2/-1/,
     1        IFIRST/0/

c** DEFINE AND INITIALIZE VARIOUS CONSTANTS AND PARAMETERS.
      LNPT= -1
      X12= 1.d0/12.D0
c** Parameters for call to VGEN in DPotFit
      IDAT= 0
      ISTATE= 1
c** XNORM  IS A FACTOR USED IN THE NORMALIZATION OF THE WAVE
c  FUNCTION TO A CONSTANT DENSITY OF STATES.
      XNORM= BZ/(PI*2.99792458D10)
      H2= RH*RH
      EJ= L*(L+1)*H2
      ISG= (-1)**L
      E= ZE*BZ*H**2
      SQE= DSQRT(ZE*BZ)
      NODX= 0
      NODE= 0
      STD= 1.d0
      STD2= 2.d0
      DLT2= 0.d0
      DLT3= 0.d0
      COR2= 0.d0
      COR3= 0.d0
      XND2= -.000001
      XND3= -.000001
      RMINN= RMIN-RH
      IF(IWR.GT.0) WRITE(6,613)
c** START TESTING FOR PHASE SHIFT CONVERGENCE AT NODES PAST THE
c  MESH POINT 'ILK' WHERE   E = L*(L+1)/R**2
      XLK= DSQRT(1.d0/(ZE*BZ))
      IF(L.GT.0) XLK= DSQRT(EJ/(ZE*BZ))
      ILK= (XLK-RMINN)/RH+ 0.5d0
      IF(ILK.LT.0) ILK= 1
      IF(ILK.GT.N) GO TO 1002
c** IF HAVE POTENTIAL (SUCH AS EXP-6) WITH NEGATIVE 'V(R)' AT VERY
c  SMALL  R , OR IF (V(I)-E) SO HIGH THAT 'RH' IS (LOCALLY) MUCH TOO
c  LARGE, THEN SHIFT INNER STARTING POINT OUTWARD.
      IST= 1
   10 IST= IST+1
      IF(IST.GE.(N-8)) GO TO 1001
      GI= (V(IST)+ EJ*XM2(IST)-E)
      IF((GI.LE.0.d0).OR.(GI.GT.11.D0)) GO TO 10
      IF((IWR.GT.0).AND.(IST.GT.2)) WRITE(6,601) L,ZE,IST
      GN= (V(IST-1)+ EJ*XM2(IST-1)-E)
      IST= IST+1
      IF(GN.LE.0.d0) GO TO 10
c** INITIALIZE WAVE FUNCTION FOR OUTWARD INTEGRATION.
      IF((RMINN+(IST-2)*RH).LE.0.d0) GO TO 9
c** IF INTEGRATION STARTS AT  R.GT.0 , INITIALIZE WAVE FUNCTION AS
c  THE EXPONENTIALLY INCREASING JWKB SOLUTION.
      SRTGI= DSQRT(GI)
      SRTGN= DSQRT(GN)
      SB= 1.d0
      SI= DSQRT(SRTGN/SRTGI)*DEXP(0.5d0*(SRTGN+SRTGI))
      IF(SI.GE.SB) GO TO 11
      WRITE(6,606) L,ZE,SB,SI
    9 SB= 0.d0
      SI= 1.d0
   11 Y1= SB*(1.d0-GN*X12)
      Y2= SI*(1.d0-GI*X12)
      RATIN= SI
      SM= 0.d0
c** INTEGRATE OUT TO MESH POINT -ILK- OR CLASSICALLY ALLOWED REGION
      DO 20 I= IST,N
          II= I
          SC= SB
          SB= SI
          Y3= Y2+Y2-Y1+GI*SB
          GI= (V(I)+EJ*XM2(I)-E)
          SI= Y3/(1.d0-GI*X12)
          IF(SI.LE.1.D+16) GO TO 12
c** NORMALIZE IN CASE OF OVERFLOW FROM WAVE FUNCTION GROWTH UNDER POTENT
          SP= 1.d0/SI
          IF(DABS(RATIN).GT.1.D-16) RATIN= RATIN*SP
          SC= SC*SP
          SB= SB*SP
          SI= 1.d0
          Y3= Y3*SP
          Y2= Y2*SP
          SM= SM*SP*SP
   12     Y1= Y2
          Y2= Y3
          SM= SM+SI*SI
c** TEST FOR AND COUNT NODES
          IF(SB) 14,16,14
   14     IF((SB*SC).GE.0.d0) GO TO 18
   16     NODE= NODE+1
          VI= V(I)+EJ*XM2(I)
          VB= V(I-1)+EJ*XM2(I-1)
          IF(VI.GE.VB) NODX= NODE
   18     IF((GI.LT.0.d0).AND.(I.GE.ILK)) GO TO 22
   20     CONTINUE
      GO TO 1002
c** AT THIS POINT, NORMALIZE BEFORE PROCEEDING
   22 SP= 1.d0/(DABS(SI)+DABS(SC))
      IF(DABS(RATIN).GT.1.D-16) RATIN= RATIN*SP
      SC= SC*SP
      SB= SB*SP
      SI= SI*SP
      Y1= Y1*SP
      Y2= Y2*SP
      SM= SM*SP*SP
      IF(IWR.GT.0) WRITE(6,611)
c** NUMERICAL INTEGRATION CALCULATING TRIAL PHASE SHIFTS STARTS HERE
      ILK1= II+1
      DO 70 I= ILK1,N
      Y3= Y2+Y2-Y1+GI*SI
      GI= (V(I)+EJ*XM2(I)-E)
      SN= Y3/(1.d0-GI*X12)
c** LOOK FOR WAVE FUNCTION NODE
      IF(SI.NE.0.d0) THEN
          IF((SB*SI).GE.0.d0) GO TO 60
          ENDIF
      NODE= NODE+1
      VN= V(I)+EJ*XM2(I)
      VI= V(I-1)+EJ*XM2(I-1)
      IF(VN.GE.VI) NODX= NODE
c** EVALUATE NODAL PHASE SHIFTS AND TEST FOR CONVERGENCE.
      JN= I
      XJN= RMINN+JN*RH
      XJC= XJN-3.d0*RH
      XND1= XND2
      XND2= XND3
      XND3= XJN+RH*(SN*SB*SC/((SI-SN)*(SI-SB)*(SI-SC))
     1  +2.d0*SN*SI*SC/((SB-SN)*(SB-SI)*(SB-SC))+3.d0*SN*SI*SB/((SC-SN)*
     2  (SC-SI)*(SC-SB)))
      IF((XND3.GT.XJC).AND.(XND3.LT.XJN)) GO TO 30
      WRITE(6,612)  NODE,L,ZE,JN,XND3
      XND3= XJN-RH*(1.d0+SI/(SI-SB))
      IF((XND3.GT.XJC).AND.(XND3.LT.XJN)) GO TO 30
      WRITE(6,610)  NODE,L,ZE,JN,XND3
      GO TO 60
   30 ARG= XND3*SQE
      CALL BSELFX(L,ARG,SBJ,SBN)
      DLTX= DATAN(SBJ/SBN)
      DLT1= DLT2
      DLT2= DLT3
      DLT3= DLTX
c** ADD JWKB PHASE SHIFT CORRECTION DUE TO POTENTIAL TAIL BEYOND 'XND3'
      COR1= COR2
      COR2= COR3
c** START QUADRATURE CALCULATION OF TRUNCATION CORRECTION
c** FIRST GET (DIMENSIONLESS) POTENTIAL AT APPROPRIATE GAUSSIAN POINTS
      DO  IK= 1,NGT
          RR(IK)= XND3/XX(IK)
          RM2(IK)= 1.d0/RR(IK)**2
c
c*** Fix for DPotFit ....NOTE ASSUMES ISTATE= 1
          CALL VGEN(ISTATE,RR(IK),VV(IK),BETADIST,IDAT)
c
          ENDDO
ccc   CALL PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NGT,RR,RM2,VLIM,VV,NCN)
c
c** NOW PERFORM ACTUAL QUADRATURE
      PSC= 0.d0
      DO  IK= 1,NGT
          RRM2= 1.d0/RR(IK)**2
          EMV= E-BZ*VV(IK)-EJ*RRM2
          IF(EMV.LT.0.d0) GO TO 60
          SQEMV= DSQRT(EMV)
          PSC= PSC+WW(IK)*(SQEMV - DSQRT(E-EJ*RRM2))
          ENDDO
      COR3= PSC*0.5d0*XND3
      DLT3= DLT3+COR3
      SNRM= 1.d0/SN
      IF(LIFE.LE.0) THEN
c** When time delay NOT calculated, use a local JWKB fit to the wavefx
c  to estimate phase shift (DLTSC) & amplitude (XNRMSC) - for fun!??
c?!?!  RJL ... check this out sometime ?!?!!?
c         SGI= DSQRT(-GI)
c         SGB= DSQRT((E-V(N-1)-EJ*XM2(N-1)))
c         DLT= (SGI+SGB)*0.5d0
c         SIK= SI*DSQRT(SGI)
c         COTA= (DCOS(DLT)-SB*DSQRT(SGB)/SIK)/DSIN(DLT)
c         ARGX= DATAN(1.d0/COTA)
c         IF(ARGX.LT.0.d0) ARGX= ARGX+PI
c         DLTSC= (ARGX+NODE*PI)-DSQRT(E)*(RMIN+(I-1)*RH)+L*PI*0.5d0
c         XNRMSC= (SQE*RH/SGI)/((COTA**2+1.d0)*SIK**2/RH)
c         IF(IWR.GT.0) 
c    1         WRITE(6,666) NODE,XND3,DLTX,DLT3,COR3,XNRM,XNRMSC,DLTSC
c 666 FORMAT(I6,F15.10,F14.10,F16.10,F14.10,2G14.6,F15.10)
c++++++++++++++++++++++++
          IF(IWR.GT.0) WRITE(6,607) NODE,XND3,DLTX,DLT3,COR3
c** Phase shift convergence test ...
          IF((DABS(DLT2-DLT3).LT.PHTST).AND.(DABS(DLT1-DLT2).LT.PHTST))
     1          GO TO 100
        ELSE
c** CALCULATE APPARENT COLLISIONAL TIME DELAY AT NODE -NODE-.
          IF(E.LE.(EJ*XM2(JN-1))) GO TO 60
          STD1= STD2
          STD2= STD
          RG2= XJN*SQE
          CALL BSELFX(L,RG2,SBJ,SBN)
c** Find JWKB-corrected local phase shift at XJN
          DSA1= DSQRT(E-VN)
          DSA2= DSQRT(E-VI)
          DSI1= DSQRT(E-EJ*XM2(JN))
          DSI2= DSQRT(E-EJ*XM2(JN-1))
          DLTY= DLTX+ 0.5d0*((DSA1-DSI1)*(2.d0*RH-(XJN-XND3)) + 
     1                           (DSA2-DSI2)*(XJN-XND3))*(XJN-XND3)/RH
c** Get JWKB-corrected effective local amplitude factor at XJN
          XNRM= (DSI1/DSA1)*((DCOS(DLTY)*SBJ-DSIN(DLTY)*SBN)*RG2/SN)**2
          ANR2= XNORM/SQE
          STD= (SM+ 0.5d0*SN*SN)*RH*XNRM - 
     1                0.5d0*(XJN- 0.5d0*ISG*DSIN(2.d0*(RG2+DLTY))/SQE)
          CALL CORTD(L,RG2,DLTY,SCOR)
          SCOR= SCOR*ANR2*DSI1/(DSA1*SQE)
          STD= STD*ANR2+SCOR
          SNRM= DSQRT(XNRM)
          IF(IWR.GT.0) WRITE(6,607) NODE,XND3,DLTX,DLT3,COR3,STD,SCOR
c** Time Delay convergence test ...
          IF((DABS(STD2/STD-1.d0).LT.TDTST).AND.
     1                      (DABS(STD1/STD2-1.d0).LT.TDTST)) GO TO 100
          ENDIF
c++++++++++++++++++++++

   60 SC= SB
      SB= SI
      SI= SN
      SM= SM+SN*SN
      Y1= Y2
   70 Y2= Y3
c** End of integration loop which calculates nodal phase shifts.
c* Since nodal phase shift (and amplitude) not converged, use local
c  JWKB approximation to estimate them.
      SGI= DSQRT(-GI)
      SGB= DSQRT(E-V(N-1)-EJ*XM2(N-1))
      DLT= (SGI+SGB)*0.5d0
      SIK= SI*DSQRT(SGI)
      COTA= (SIK/(SB*DSQRT(SGB))-DCOS(DLT))/DSIN(DLT)
      XNRM= 1.d0/((COTA**2+1.d0)*SIK**2)
      WRITE(6,609) DLT3,COR3,DLT2,COR2,DLT1,COR1,XND3,XND2,XND1
c** GET ABSOLUTE PHASE SHIFT USING EQ.(9.2.29) OF ABRAMOWITZ AND STEGUN
  100 ZU= 4.d0*EJ+1.d0
      R4= 4.d0*ARG
      SIR4= 1.d0/R4**2
      TH= ARG- 0.5d0*L*PI+ (0.5d0+ SIR4*((ZU-25.D0)/6.d0+ 
     1   SIR4*((1073.D0+ ZU*(-114.D0+ ZU))/5.d0 + SIR4*(-375733.D0 + 
     2   ZU*(54703.D0+ZU*(-1535.D0 + ZU*5.d0)))/14.D0)))*(ZU- 1.d0)/R4
      XNDT= (TH+DLT3)/PI
      NDT= (XNDT+0.5d0)
      ADLT= DLT3+(NODE-NDT)*PI
      RATIN= RATIN*SNRM
      IF(DABS(RATIN).GT.1.0D-10) WRITE(6,608) L,ZE,RATIN
      IF(IWR.NE.0) THEN
          WRITE(6,603) L,ZE,ADLT,NODE,XND3,COR3
          IF(LIFE.GT.0) THEN
              SCRT= SCOR/STD
              GMA= 2.d0/(STD*2.99792458D10*PI)
              IF(GMA.LT.0.d0) GMA= 0.d0
              WRITE(6,604) STD,GMA,SCRT
              ENDIF
          ENDIF
  999 RETURN
 1001 WRITE(6,602) L,ZE,IST
      GO TO 999
 1002 WRITE(6,605) L,ZE,XLK,ILK
      GO TO 999
  601 FORMAT(' L= ',I3,'    E=',F9.3,'   integration cannot start till p
     1ast mesh point #',I5)
  602 FORMAT(/' *** ERROR *** FOR   L =',I3,'   E =',F9.3/15x,'Integrati
     1on apparently cannot start till mesh point',I5)
  603 FORMAT(/' Phase Sh [L=',I3,', E=',F9.3,']=',F12.8,'   X(Node',I3,
     1  ')=',F7.3,'   Cor=',1PD9.2)
  604 FORMAT(' Time Delay =',1PD14.7,'(sec)   FWHM=',D12.5,
     1  '(cm-1)    CORRN/STD=',D8.1)
  605 FORMAT(/' *** ERROR ***   For   L= ',I3,'   E=',F9.3,'   XLK=',
     1  f7.3,'   ILK=',I8/17x,'Reach end of range before start search f
     2or nodal phase shifts')
  606 FORMAT(/' *** WARNING ***  For  L=',I3,'   E=',F9.3,'   JWKB initi
     1alization gives   SB/SI=',D10.3,'/',D10.3,',  so set   SB=0.')
  607 FORMAT(I5,F13.8,2F15.10,1PD11.3,D16.7,'(',D12.5,')')
  608 FORMAT(/' *** CAUTION    L= ',I3,'   E=',F9.3,'   wave function ma
     1y be unreliable as   RATIN= ',D10.3,',  so make RMIN smaller!')
  609 FORMAT(' *** CAUTION *** NOT CONVERGED BY END OF RANGE'/8x,
     1 'DLT= ',3(0PF11.8,'(',1PD9.1,')')/7X,'@  X= ',3(0PF12.8,10X))
  610 FORMAT(/' INSTABILITY *** in linear interpolation at  NODE= ',
     1  I3,'   FOR  L= ',I3,'   E=',E11.5,'(cm-1)   JN=',I4,'   XND3=',
     2  F10.6)
  611 FORMAT(/'  NODE    X(NODE)      DLT(NODE)      DLT(COR3)     COR3'
     1 ,'       STD(SEC)'/1X,39('=='))
  612 FORMAT(/' INSTABILITY *** in cubic interpolation at   NODE=',i3,
     1  '   For  L= ',I3,'   E=',E11.5,'(cm-1)   JN=',I4,'   XND3=',
     2  F10.6)
  613 FORMAT(1X,39('--'))
      END
c*******************************************************************
      SUBROUTINE CORTD(L,ARG,DLT,SUM)
c** SUBROUTINE TO CORRECT COLLISIONAL TIME DELAY FOR RESIDUAL
c  DIFFERENCES BETWEEN THE APPROPRIATE COMBINATION OF BESSEL FUNCTIONS
c  AND THE PHASE-SHIFTED SINE FUNCTION AT   R .GT. ARG/DSQRT(E)
      INTEGER IFIRST,I,IG,IJ,IM,ISG,J,J1,K,L,LMX,LPREV,M,N,NM1,NM2,
     1        NNL(66)
      REAL*8 LD
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SNL(2211),NNL(66),A(66),B(66),C(66)
      DATA IFIRST/0/
      SAVE LPREV,IFIRST,LMX,NNL,SNL
      SUM= 0.d0
      IF(IFIRST.GT.0) GO TO 20
c** ON VERY FIRST CALL TO -CORTD- SET UP ARRAY OF NORMALIZED COEFFICIENT
c  SNL(J*(J+1)/2+M+1)= N(J,M)/L**M   OF THE RAYLEIGH POWER SERIES
c  EXPANSION FOR SPHERICAL BESSEL FUNCTIONS.  AS IN RJL-RBB PAPER,
c  Z*J(L,Z)= (-1)**L *SUM( SNL(NNL(L)+M+1)*SIN(Z+(L-M)*PI/2)*(L/Z)**M )
c** ARRAYS DIMENSIONED FOR  L  UP TO  LMX= 65 .
      IFIRST= 1
      LPREV= -1
      LMX= 65
      NNL(1)= 2
      SNL(1)= 1.d0
      SNL(2)= 1.d0
      SNL(3)= -1.d0
      IJ= 2
      DO  J= 2,LMX
          IJ= IJ+J
          NNL(J)= IJ
          SNL(IJ)= 1.d0
          J1= J-1
          IM= IJ+1
          ZJ= J
          SNL(IM)= (1.d0-IJ)/ZJ
          IF(J.NE.2) THEN
              RT= (ZJ-1.d0)/ZJ
              F= RT
              DO  M= 2,J1
                  IM= IM+1
                  F= F*RT
                  SNL(IM)= (SNL(IM-J)-SNL(IM-J-1)*(J+M-1)/(ZJ-1.d0))*F
                  ENDDO
              ENDIF
          IM= IM+1
          SNL(IM)= -SNL(IM-1)/ZJ
          ENDDO
   20 IF(L.LE.0) GO TO 99
c** GENERATE NORMALIZED (BY L**(2*M)) VALUES OF COEFFICIENTS A(L,M),
c  B(L,M) &C(L,M) OF EQS.(A8) & (A9) IN RJL-RBB PAPER.
c** IF IMMEDIATELY PRECEEDING CASE INVOLVED SAME VALUE OF -L-,
c  THEN SKIP THIS PART.
      IF(L.EQ.LPREV) GO TO 50
      LPREV= L
      IF(L.GT.LMX) THEN
          WRITE(6,602) L,LMX
          GO TO 99
          ENDIF
      ISG= (-1)**L
      DO  M= 1,L
          ISG= -ISG
          NM1= MIN0(M,L-M+1)
          NM2= MIN0(M,L-M)
          AM= 0.d0
          I= NNL(L)+M
          DO  N= 1,NM1
              AM= AM+SNL(I+N-1)*SNL(I-N)
              ENDDO
          A(M)= ISG*AM/L
          BM= 0.5d0*SNL(I)**2
          CM= BM
          IG= 1
          IF(NM2.GE.1) THEN
              DO  N= 1,NM2
                  D= SNL(I+N)*SNL(I-N)
                  IG= -IG
                  BM= BM+D
                  CM= CM+IG*D
                  ENDDO
              ENDIF
          B(M)= ISG*BM
          C(M)= CM/(2*M-1)
          ENDDO
c** ACTUAL SUM OF RJL-RBB EQ.(A8) PERFORMED HERE.
   50 IF(L.GT.LMX) GO TO 99
      ZI= ARG*(L/ARG)**(L+L)
      CSZ= 0.5d0*ZI*DCOS(2.d0*(ARG+DLT))
      SNZ= 0.5d0*(ZI/ARG)*DSIN(2.d0*(ARG+DLT))
c** NOTE ** AS FOUND BY RJL (JANUARY 1978), ASIDE FROM TRUNCATION
c  EFFECTS, VALUE OF -SUM- INDEPENDENT OF INITIAL VALUE OF C(2*M)= CI
c  (OR IF SUM STARTED AT M= 1, OF ASSUMED S(1) VALUE).
      CI= 0.d0
      SI= CSZ-(L-0.5d0)*CI
      M= L
      SUM= A(M)*SI-B(M)*CI+C(M)*ZI
      DSMX= DABS(SUM)
      IF(L.LE.1) GO TO 99
      FD= 1.d0/L**2
      FN= (ARG/L)**2
      DO  K= 2,L
          M= M-1
          CSZ= CSZ*FN
          SNZ= SNZ*FN
          ZI= ZI*FN
          CI= -SNZ+M*FD*SI
          SI= CSZ-(M-0.5d0)*CI
          DSM= A(M)*SI-B(M)*CI+C(M)*ZI
c** WARNING ** FOR LARGE -L-, THESE -DSM- VALUES MAY BE MUCH LARGER
c  THAN FINAL -SUM-, SO MAY GET CONSIDERABLE LOSS OF ACCURACY DUE TO
c  CANCELLATION OF MANY OF THE SIGNIFICANT DIGITS. ** IN PARTICULAR,
c  NO. OF GOOD DIGITS APPROX. = 15 - LOG10(MAX(DSM)/SUM)  .
c** E.G. FOR L= 50, ARG=75, HAVE ONLY 4 GOOD DIGITS
c**      FOR L= 60, ARG=90, HAVE ONLY 1-2 GOOD DIGITS.
c**      FOR L.GE.63, ARG= 1.5*L, HAVE  NO  GOOD DIGITS.
          ADSM= DABS(DSM)
          IF(ADSM.GT.DSMX) DSMX= ADSM
          SUM= SUM+DSM
          ENDDO
      ADSM= DSMX/DABS(SUM)
      IF(ADSM.GE.1.D+13) THEN
          LD= DLOG10(ADSM)
          WRITE(6,601) L,ARG,LD,DSMX,SUM
          IF(LD.GE.15) SUM= 0.d0
          ENDIF
   99 RETURN
  601 FORMAT(/'  WARNING *** TRUNCATION EFFECTS IN -CORTD- LOSE',I3,
     1 ' SIGNIFICANT DIGITS, SINCE   DSMX = ',D12.5,'   SUM =',D12.5)
  602 FORMAT(/' *** WARNING *** no resudual correction to Time Delays fo
     1r   L=',I3/5x,'  Since required arrays larger than allowed by   LM
     2X=',I3)
      END
c*******************************************************************
      SUBROUTINE BSELFX(L,ARG,SBJ,SBN)
c** SUBROUTINE TO CALCULATE SPHERICAL BESSEL FUNCTIONS OF THE FIRST AND
c   SECOND KINDS,  SBJ  AND  SBN  , RESPECTIVELY, OF ORDER  L  , WITH
c   ARGUMENT  ARG.
      INTEGER L,LH,LE,L4,I,IQ,K
      IMPLICIT REAL*8 (A-H,O-Z)
      LH= L/2
      PM= 1.0D0
      LE= L
      AR1= 1.D0/ARG
      XQ= DSQRT(1.570796326794896*ARG)
      AQ= DCOS(ARG)/XQ
      BQ= DSIN(ARG)/XQ
      XQ= XQ*AR1
   12 CQ= PM
      IF(LE) 14,28,16
   14 QQ= AQ
      RQ= BQ
      GO TO 18
   16 IF(ARG.LT.(1.D0*L)) GO TO 30
      QQ= BQ
      RQ= AQ
   18 DO 20 IQ= 1,L
      PQ= QQ*(CQ+CQ-PM)*AR1-RQ
      CQ= CQ+PM
      RQ= QQ
   20 QQ= PQ
      PQ= PQ*XQ
      IF(PM.LT.0.) GO TO 26
      SBJ= PQ
   22 PM= -1.0D0
      LE= -L
      GO TO 12
   26 SBN= PQ
      K= L-2*LH
      IF(K.EQ.0) SBN= -PQ
      RETURN
   28 SBJ= BQ*XQ
      SBN= -AQ*XQ
      RETURN
c** CALCULATE FX OF 1-ST KIND FOR SMALL ARGUMENT BY ASCENDING POWER
c   EXPANSION.
   30 QQ= 1.0D0
      DO  I= 1,L
          QQ= QQ*ARG/(2*I+1)
          ENDDO
      SBJ= QQ
      L4= 50
      IF(L.GT.12) L4= 4*L
      CQ= (2*L+1)
      QQ= -0.5D0*ARG*ARG
      RQ= 1.0D0
      PQ= RQ
      DO  I= 1,L4
          K= I
          RQ= RQ*QQ/(I*(CQ+2*I))
          PQ= PQ+RQ
          IF(DABS(RQ/PQ).LT.1.D-10) GO TO 36
          CONTINUE
          ENDDO
      RQ= RQ/PQ
      WRITE(6,601) L,ARG,L4,RQ
   36 SBJ= SBJ*PQ
      GO TO 22
  601 FORMAT(/'  ERROR IN   SBJ(',I3,',',E10.5,') . AFTER',I4,' ASCENDIN
     1G TERMS   RATIO = ',E10.5 )
      END
c23456789012345678901234567890123456789012345678901234567890123456789012
