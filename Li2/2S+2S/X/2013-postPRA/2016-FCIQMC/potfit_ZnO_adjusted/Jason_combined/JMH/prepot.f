c***********************************************************************
      SUBROUTINE PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,RR,RM2,VLIM,
     1                                                        VV,NCN)
c** Driver subroutine of package to read parameters and/or generate
c  values of a potential V(I) at the NPP input distances RR(I).
c====================== Version of 20 July 2001 ====================
c**** Subroutine Input:
c----------------------
c  LNPT  is an integer specifying the operational mode:
c      *  LNPT > 0  : for a new case for which all potential-defining
c                     parameters are read in & a description printed
c      *  LNPT.le.0 : if potential points are to be generated in exactly
c                     the same manner as on preceding call, but at
c                     different distances RR(I) (no reads or writes)
c  IAN1 & IAN2 are the atomic numbers and IMN1 & IMN2 the mass numbers
c        of atoms #1 & 2, used (if needed) to specify isotope masses for
c        calculating adiabatic and/or non-adiabatic B-O-B correction fx.
c  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
c        at which potential values  VV(i) (in cm-1) are to be generated
c  RR  (real array) is set of NPP distances where potential calculated
c  RM2 (real array) on input is the (centrifugal) array of  1/RR(i)**2
c----------------------
c**** Subroutine Output:
c----------------------
c  VLIM (cm-1) is the absolute energy at the potential asymptote
c  VV (real array) is the set of function values generated (in cm-1)
c  RM2 values returned may (if appropriate) be modified to include B-O-B
c      corrections to the (centrifugal) potential  1/RR(i)**2
c  NCN is an integer power defining the asymptotically-dominant 
c      inverse-power long-range potential tail:  CNN/R**NCN 
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) ,  or POTGEN ++
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Set maximum array dimension for the input function values to be
c  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600) 
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LPPOT,LWR,
     1  NCN,NLIN,NPP,NROW,NTP,NUSE,NPRS,NPRF
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP),
     1  XI(NTPMX),YI(NTPMX),RWR(20),VWR(20),VWRB(3),D1V(3),D1VB(3),
     2  D2V(3),CNN
c
c** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,LPPOT,NTP,NUSE
      SAVE CNN,VSHIFT,XI,YI
c
      DATA VWRB/3*0.D0/,D1VB/3*0.D0/
      LPPOT= 0
      NPRS= 1
      NPRF= NPP
c
      IF(LNPT.GT.0) THEN
c** If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c** If LPPOT > 0 : at every |LPPOT|-th point, print potential and 
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c** VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
          READ(5,*) NTP, LPPOT, VLIM
c-----------------------------------------------------------------------
          WRITE(6,600) VLIM
          IF(NTP.GT.0) THEN
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c-----------------------------------------------------------------------
              READ(5,*) NUSE, IR2, ILR, NCN, CNN
c-----------------------------------------------------------------------
              IF(NTP.GT.NTPMX) THEN
                  WRITE(6,602) NTP,NTPMX
                  STOP
                  ENDIF
              IF(NUSE.GT.0) WRITE(6,604) NUSE,NTP
              IF(NUSE.LE.0) WRITE(6,606) NTP
              IF(IR2.GT.0) WRITE(6,608)
              IF((ILR.GT.1).AND.(DABS(CNN).GT.0.D0))WRITE(6,610)CNN,NCN
c** Read in turning points to be interpolated over
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (may = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
              READ(5,*) RFACT, EFACT, VSHIFT
              READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
              WRITE(6,612) RFACT,EFACT
              NROW= (NTP+2)/3
              DO  J= 1,NROW
                  IF(EFACT.LE.10.D0) THEN
                      WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
                    ELSE
                      WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
                    ENDIF
                  ENDDO
              WRITE(6,618) VSHIFT
              DO  I= 1,NTP
                  YI(I)= YI(I)*EFACT+ VSHIFT
                  XI(I)= XI(I)*RFACT
                  ENDDO
              IF(IR2.GT.0) THEN
                  DO  I= 1,NTP
                      YI(I)= YI(I)*XI(I)**2
                      ENDDO
                  ENDIF
              ENDIF
          ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NTP.GT.0) THEN
          CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
        ELSE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Generate a fully analytic potential in subroutine POTGEN ***********
c* Potentials generated in cm-1 with potential asymptote at energy VLIM
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR & NPAR are integers for specifying potential types.
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass 
c      dependent Born-Oppenheimer breakdown corrections will be included
c** For all functions considered, well depth and equilibrium distance
c  are read as  DSCM (cm-1)  and REQ (Angstroms), respectively.
c* [Most read-in parameters are dimensionless (scaled by DSCM & REQ).]
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NPAR) potential.
c** If IPOTL=2  generate an MLJ(NPAR) potential [JCP 112, 3949 (2000)]
c      If MPAR.ge.0  exponent parameter is polynomial of order (NVARB-1)
c           in z=(R-Re)/(R+Re), with the NVARB coefficients  PARM(j)
c      If MPAR < 0  exponent polynomial in  z  has order (NVARB-4) with 
c           coefficients PARM(i) (i= 1,NVARB-3), & includes a switching
c           function with exponent coefficient  ALPHA= PARM(NVARB)  and  
c           RSW= PARM(NVARB-1),  defined to yield limiting inverse-power
c           potential coefficient  Cn= PARM(NVARB-2).
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c       with exponent factor "beta" defined as a power series of order
c       (NVARB-1) in  z=(R-Re)/(R+Re)  with NVARB coefficients PARM(i).
c       Set  NVARB= 1 for conventional "simple" Morse potential.
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter 
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  use Seto's modification of Surkus' GPEF expansion in
c       z = [R^NPAR - Re^NPAR]/[a*R^NPAR + b*Re^NPAR] where 
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NPAR=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NPAR=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NPAR=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus NPAR < 0 case:  z(NPAR,a,b)= z(|NPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=5  generate generalized HFD(NPAR,6,8,10,12,14) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c       PARM(6-11)  are the reduced C_NPAR, C_6, C_8, C_10, C_12 and C14
c       parameters (NPAR < 6), while  AREP and  beta  are defined
c       by having the potential minimum at  x=1.  For NVARB < 11, higher
c       C_m coefficients automatically zero;  necessarily  NVARB.ge.7 .
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** IBOB > 0, add atomic-mass-dependent Born-Openheimer breakdown 
c  correction functions to rotationless and/or centrifugal potential(s).
c  Both expressed as power series in  z= (R-Re)/(R+Re) starting with the
c  constant term, using the mass shift convention of Le Roy [J.Mol.Spec. 
c  194, 189 (1999)].  Adiabatic B-O-B potential correction fx. defined 
c  by polynomials of order NC1 with (NC1+1) coefficients {CA1(i)} for 
c  atom-1  and order NC2 with (NC2+1) coefficients {CA2(i)} for atom-2, 
c  while centrifugal correction fx. defined polynomial of order NG1 with
c  (NG1+1) coefficients {GA1(i)} for atom-1 and order NG2 with (NG2+1)
c  coefficients {GA2(i)} for atom-2.
c** Input parameters IANi & IMNi are the atomic & mass number of atom-i
c  (i=1,2), while integers RMN1 & RMN2 read here are the mass numbers of
c  the reference isotopes defining the B-O-B correction functions.
c** NC1 & NC2 are orders of polynomials DELTA(V,atom-i) defining 
c  'adiabatic' corrections to the rotationless potential for atoms 1 & 2
c  DELTA(V)= (1-M1ref/M1)*DELTA(V,atom-1) + (1-M2ref/M2)*DELTA(V,atom-2)
c** NG1 & NG2 are orders of polynomials q1(z) & q2(z) defining B-O-B
c   correction to the centrifugal potential:
c      V(centrifugal)= [1 + (M1ref/M1)*q1(z) + (M2ref/M2)*q2(z)]/R**2
c ... to omit a particular correction set associated NCi or NGi .lt.0
c** RX > 0.0  invokes Coxon's (older) expansions in (R-Re) for potential
c     correction and in  [(R-Rx)**j - (Re-Rx)**j] for centrifugal corrn.
c ... OTHERWISE (to use Le Roy B-O-B formalism) set  RX.le.0.d0 !!
c-----------------------------------------------------------------------
c** Read inside subroutine POTGEN
c         IF(LNPT.GT.0) THEN
c             READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
c             IF(NVARB.GT.0) READ(5,*) (PARM(I), I=1,NVARB)
c             IF(IBOB.GT.0) THEN
c                 READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c                 IF(NC1.GE.0) READ(5,*) (CA1(I), I=0,NC1)
c                 IF(NC2.GE.0) READ(5,*) (CA2(I), I=0,NC2)
c                 IF(NG1.GE.0) READ(5,*) (GA1(I), I=0,NG1)
c                 IF(NG2.GE.0) READ(5,*) (GA2(I), I=0,NG2)
c             ENDIF
c         ENDIF
c-----------------------------------------------------------------------
          NCN= 99
          CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,
     1                                                        NCN,CNN)
        ENDIF
      IF(LPPOT.NE.0) THEN
c** If desired, on the first pass (i.e. if LNPT > 0) print the potential
          RH= RR(2)-RR(1)
          INPTS= IABS(LPPOT)
          IF(LPPOT.LT.0) THEN
c** Option to write resulting function compactly to channel-8. 
              RMIN= RR(1)
              NLIN= NPP/INPTS+ 1
              WRITE(8,800) NLIN,VLIM
              WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
            ELSE
c** Option to print potential & its 1-st three derivatives, the latter
c  calculated by differences, assuming equally spaced RR(I) values.
              WRITE(6,620)
              NPRS= MAX(1,(NPRS- 15*INPTS))
              NPRF= MIN(NPP,(NPRF+ 15*INPTS))
              NLIN= (NPRF-NPRS+1)/(2*INPTS)+1
              RH= INPTS*RH
              DO  I= 1,NLIN
                  LWR= NPRS+INPTS*(I-1)
                  DO  J= 1,2
                      JWR= LWR+(J-1)*NLIN*INPTS
                      IF(JWR.LE.NPP) THEN
                          RWR(J)= RR(JWR)
                          VWR(J)= VV(JWR)
                          D1V(J)= (VWR(J)-VWRB(J))/RH
                          VWRB(J)= VWR(J)
                          D2V(J)= (D1V(J)-D1VB(J))/RH
                          D1VB(J)= D1V(J)
                        ELSE
                          RWR(J)= 0.d0
                          VWR(J)= 0.d0
                        ENDIF
                      IF(I.LE.2) THEN
                          D2V(J)= 0.d0
                          IF(I.EQ.1) D1V(J)= 0.d0
                          ENDIF
                      ENDDO
                  WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
                  ENDDO
            ENDIF
          ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' Absolute energy at asymptote:   Y(lim)=',F12.4,'(cm-1)')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'
     1 ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  606 FORMAT(' Perform cubic spline interpolation over the',I5,' input p
     1oints' )
  608 FORMAT(' Interpolation performed over modified input array:   Y(I)
     1 * R(I)**2')
  610 FORMAT(/' Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  612 FORMAT(/' Scale input points:  (distance)*',1PD16.9,'   &  (energy
     1)*',D16.9/13x,'to get required internal units  [Angstroms & cm-1 f
     2or potentials]'/3('      R(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  618 FORMAT(3x,24('---')/' To make above input points consistent with',
     1  '  Y(lim),  add  Y(shift)=',F12.4)
  620 FORMAT(/'  Function and first 2 derivatives by differences'/
     1  2('     R       Y(R)     d1Y/dR1    d2Y/dR2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
  624 FORMAT(1x,39('--'))
  800 FORMAT(I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,3(F12.8,F14.6)))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,
     1                                              NCN,CNN,NPRS,NPRF)
c** GENINT produces a smooth function YY(i) at the NPP input distances
c  XX(i) by performing numerical interpolation over the range of the 
c  NTP input function values YI(j) at the distances XI(j), and using
c  analytic functions to extrapolate beyond their range to with an
c  exponential at short range and a form specified by ILR, NCN & CNN
c** ILR specifies how to extrapolate beyond largest given turning pts
c   If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c   If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c   If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c* If(ILR.ge.2) fit last turning points to:  VLIM - sum(of ILR
c  inverse-power terms beginning with  1/R**NCN). *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  ((cm-1)(Angstroms)**'NCN').
c  If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c  If ILR > 3 : this factor is  1/R .
c=== Calls subroutines PLYINTRP, SPLINT & SPLINE ==================
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      INTEGER  I,J,IFXCN,IDER,IR2,ILR,ISR,LNPT,MBEG,MFIN,MINNER,
     1  NN,NPP,NUSE,NUST,NORD,NCN,NCN2,NCN4,NPRS,NPRF,NTP,
     2  IMX1,NMX,JR2,JMAX,MI(10),MF(10)
      REAL*8  ASR,BSR,CSR,ALR,BLR,CLR,DCSR,ADCSR,PDCSR,VRAT,
     1  DX1,DX2,DX3,EX1,EX2,EX3,CNN,VLIM,X1,X2,X3,Y1,Y2,Y3,
     1  XX(NPP),YY(NPP),XI(NTP),YI(NTP),XJ(20),YJ(20),DUMM(20)
c
      SAVE ASR,BSR,CSR,ISR,ALR,BLR,CLR,IMX1,NMX,JR2,JMAX
c
      NUST= NUSE/2
      IF(NUSE.LE.0) NUST= 2
      IDER= 0
      NPRS= 1
      NPRF= NPP
c** Determine if/where need to begin extrapolation beyond input data
c  XX(MI(J))  is the 1-st mesh point past turning point  XI(J) .
c  XX(MF(J))  is the last mesh point before turning point  XI(NTP+1-J)
      DO 6 J = 1,NUST
          MI(J)= 1
          MF(J)= 0
          DO  I= 1,NPP
              IF(XX(I).LE.XI(J)) MI(J)= I+ 1
              IF(XX(I).GE.XI(NTP+1-J)) GO TO 6
              MF(J)= I
              ENDDO
    6     CONTINUE
      IF(NUST.LT.2) THEN
          MFIN= MI(1)-1
        ELSE
          MFIN= MI(2)-1
        ENDIF
      IF(LNPT.GT.0) THEN
c-----------------------------------------------------------------------
c** For a new case determine analytic functions for extrapolating beyond
c  the range of the input points (if necessary) on this or later calls.
c** Try to fit three innermost turning points to  V(R)=A+B*DEXP(-C*R).
c** If unsatisfactory, extrapolate inward with inverse power function
          IF(IR2.LE.0) THEN
              DO  I= 1,4
                  YJ(I)= YI(I)
                  ENDDO
            ELSE
              DO  I= 1,4
                  YJ(I)= YI(I)/XI(I)**2
                  ENDDO
            ENDIF
          X1= XI(1)
          X2= XI(2)
          X3= XI(3)
          Y1= YJ(1)
          Y2= YJ(2)
          Y3= YJ(3)
          IF((Y1-Y2)*(Y2-Y3).LE.0.d0) THEN
c** If 3 innermost points not monotonic, use A+B/X inward extrapoln.
              ISR= 0
              WRITE(6,600)
            ELSE
c** Use cubic through innermost points to get initial trial exponent
c  from ratio of derivatives,  Y''/Y'
              IDER= 2
              ISR= 4
              CALL PLYINTRP(XI,YJ,ISR,X2,XJ,ISR,IDER)
              CSR= XJ(3)/XJ(2)
              DCSR= DABS(CSR*X2)
              IF(DCSR.GT.1.5D+2) THEN
c** If exponential causes overflows, use inverse power inward extrapoln.
                  ISR= 0
                  WRITE(6,602) CSR
                  GO TO 20
                  ENDIF
c** Prepare parameters for inward exponential extrapolation
              VRAT= (Y3- Y2)/(Y1- Y2)
              DX1= X1- X2
              DX3= X3- X2
              EX2= 1.D0
              ADCSR= 1.d99
c** Now iterate (with actual point) to get exact exponent coefficient 
              DO  J= 1,15
                  PDCSR= ADCSR
                  EX1= DEXP( CSR*DX1)
                  EX3= DEXP( CSR*DX3)
                  DCSR= (VRAT- (EX3- EX2)/(EX1- EX2)) /
     1   ((X3*EX3- X2 - (X1*EX1- X2)*(EX3-EX2)/(EX1- EX2))/(EX1- EX2))
                  ADCSR= ABS(DCSR)
                  IF((ADCSR.GT.PDCSR).AND.(ADCSR.LT.1.d-8)) GO TO 12
                  IF(ADCSR.LT.1.d-12) GO TO 12
                  CSR= CSR+ DCSR 
                  ENDDO
              WRITE(6,604) DCSR
   12         BSR= (Y1-Y2)/(EX1-EX2)
              ASR= Y2-BSR*EX2
              BSR= BSR*DEXP(-CSR*X2)
              WRITE(6,606) X2,ASR,BSR,CSR
            ENDIF
   20     IF(ISR.LE.0) THEN
              IF((X1*X2).LE.0.d0) THEN
c** If 1'st two mesh points of opposite sign, extrapolate linearly
                  ISR= -1
                  ASR= Y2
                  BSR= (Y2- Y1)/(X2- X1)
                  CSR= X2
                  WRITE(6,608) X2,ASR,BSR,CSR
                ELSE
c** For inward extrapolation as inverse power through 1'st two points ..
                  BSR= (Y1-Y2)* X1*X2/(X2- X1)
                  ASR= Y1-BSR/X1
                  CSR= X2
                  WRITE(6,610) X2,ASR,BSR
                ENDIF
              ENDIF
          ENDIF
  600 FORMAT('  ** CAUTION ** Exponential inward extrapolation fails sin
     1ce first 3 points not monotonic, ... so ...')
  602 FORMAT(' *** CAUTION ** inward extrapolation exponent coefficient
     1   C=',D12.4/10x,'could cause overflows, ... so ...')
  604 FORMAT(' *** CAUTION ** after 15 tries inward extrap. exponent coe
     1fft change is',1PD9.1)
  606 FORMAT(' Extrapolate to   X .le.',F7.4,'  with'/'   Y=',F13.3,
     1  SP,1PD15.6,' * exp(',SS,D13.6,'*R)')
  608 FORMAT(' Extrapolate to   X .le.',F8.4,'   with'/'   Y=',F13.3,
     1  SP,1PD16.7,' * [R - (',SS,F8.4,')]')
  610 FORMAT(' Extrapolate to   X .le.',F8.4,'   with   Y=',F12.3,
     1  SP,1PD15.6,')/R**1')
c
      IF(MFIN.GT.0) THEN
c** If needed, calculate function in inner extrapolation region
          IF(ISR.GT.0) THEN
c ... either as an exponential
              DO  I= 1,MFIN
                  EX1= CSR*XX(I)
                  IF(DABS(EX1).GT.1.D+2) EX1= 1.D+2*DSIGN(1.d0,EX1)
                  YY(I)= ASR+BSR*DEXP(EX1)
                  ENDDO
            ELSEIF(ISR.EQ.0) THEN
c ... or if that fails, as an inverse power
              DO  I= 1,MFIN
                  YY(I)= ASR+BSR/XX(I)
                  ENDDO
            ELSEIF(ISR.LT.0) THEN
c ... or if X changes sign, extrapolate inward linearly
              DO  I= 1,MFIN
                  YY(I)= ASR+ BSR*(XX(I)- CSR)
                  ENDDO
            ENDIF
          ENDIF
c** End of inward extrapolation procedure
c-----------------------------------------------------------------------
      MINNER= MFIN
      IF(NUST.GT.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near inner end of range
          DO  J= 3,NUST
              NORD= 2*(J-1)
              MBEG= MI(J-1)
              MFIN= MI(J)-1
              IF(MFIN.GE.MBEG) THEN
                  DO  I=  MBEG,MFIN
                      CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                      YY(I)= DUMM(1)
                      ENDDO
                  ENDIF
              ENDDO
          ENDIF
c** Main interpolation step begins here
c=======================================================================
      MBEG= MI(NUST)
      MFIN= MF(NUST)
      IF(MFIN.GE.MBEG) THEN
          IF(NUSE.LE.0) THEN
c** Either ... use cubic spline for main interpolation step
              CALL SPLINT(LNPT,NTP,XI,YI,MBEG,MFIN,XX,YY)
            ELSE
c ... or use piecewise polynomials for main interpolation step
              DO  I= MBEG,MFIN
                  CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NUSE,IDER)
                  YY(I)= DUMM(1)
                  ENDDO
            ENDIF
          ENDIF
      IF(MFIN.LT.NPP) THEN
          IF(NUST.LE.2) THEN
c** If(NUSE.gt.5) minimize spurious behaviour by interpolating with
c  order less than NUSE on intervals near outer end of range
              MBEG= MF(NUST)+1
            ELSE
              NN= NUST-2
              DO  J= 1,NN
                  NORD= 2*(NUST-J)
                  MBEG= MF(NUST-J+1)+1
                  MFIN= MF(NUST-J)
                  IF(MFIN.GE.MBEG) THEN
                      DO  I= MBEG,MFIN
                          CALL PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
                          YY(I)= DUMM(1)
                          ENDDO
                      END IF
                  ENDDO
            ENDIF
          ENDIF
      MBEG= MFIN+1
      IF((MFIN.GT.MINNER).AND.(IR2.GT.0)) THEN
c** In (IR2.gt.0) option, now remove X**2 from the interpolated function
          DO  I= MINNER+1,MFIN
              YY(I)= YY(I)/XX(I)**2
              ENDDO
          ENDIF
c** Print test of smoothness at join with analytic inward extrapolation
c     IF(LNPT.GT.0) THEN
c         MST= MAX0(MINNER-4,1)
c         NPRS= MST
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         IF(MFN.GT.MFIN) MFN= MFIN
c         IF(MINNER.GT.0) WRITE(6,611) X2,((XX(I),YY(I),I= J,MFN,3),
c    1        J= MST,MST+2)
c 611 FORMAT('     Verify smoothness of inner join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
c         ENDIF
c-----------------------------------------------------------------------
c** To extrapolate potential beyond range of given turning points ...
      IF(LNPT.GT.0) THEN
c** On first entry, calculate things needed for extrapolation constants
          Y1= YI(NTP)
          Y2= YI(NTP-1)
          Y3= YI(NTP-2)
          X1= XI(NTP)
          X2= XI(NTP-1)
          X3= XI(NTP-2)
          IF(IR2.GT.0) THEN
              Y1= Y1/X1**2
              Y2= Y2/X2**2
              Y3= Y3/X3**2
              ENDIF
          ENDIF
c** Check inverse-power tail power ...
      IF(NCN.LE.0) NCN= 6
      IF(ILR.LT.0) THEN
          IF(LNPT.GT.0) THEN
C** For  ILR.lt.0  use  Y = VLIM - ALR * exp[-CLR*(X - BLR)**2]
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              BLR= (X1+X2 - (X2+X3)*EX1/EX2)/(2.d0- 2.d0*EX1/EX2)
              CLR= -EX1/(X1+X2-2.d0*BLR)
              ALR= (VLIM-Y1)*DEXP(CLR*(X1-BLR)**2)
              WRITE(6,614) X2,VLIM,ALR,CLR,BLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*DEXP(-CLR*(XX(I) - BLR)**2)
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.0) THEN
c** For ILR.le.0  use  Y = VLIM - ALR * X**p * exp(-CLR*X)
          IF(LNPT.GT.0) THEN
              EX1= DLOG((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
              EX2= DLOG((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
              DX1= DLOG(X1/X2)/(X1-X2)
              DX2= DLOG(X2/X3)/(X2-X3)
              BLR= (EX1-EX2)/(DX1-DX2)
              CLR= BLR*DX1- EX1
              ALR= (VLIM-Y1)* DEXP(CLR*X1)/X1**BLR 
              WRITE(6,616) X2,VLIM,ALR,BLR,CLR
              IF(CLR.LT.0.d0) THEN
c ... but replace it by an inverse power of exponent constant negative
                  WRITE(6,612)
                  ILR= 1
                  GO TO 50
                  ENDIF
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM- ALR*XX(I)**BLR *DEXP(-CLR*XX(I))
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
   50 IF(ILR.EQ.1) THEN
c** For  ILR=1 ,  use     Y = VLIM + ALR/X**BLR
          IF(LNPT.GT.0) THEN
              BLR= DLOG((VLIM-Y2)/(VLIM-Y1))/DLOG(X1/X2)
              ALR= (Y1- VLIM)*X1**BLR
              NCN= BLR
              WRITE(6,618) X2,VLIM,ALR,BLR,NCN
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM+ ALR/XX(I)**BLR
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
c** Set constants for long-range extrapolation
      IFXCN= 0
      IF((CNN.GT.0.d0).OR.(CNN.LT.0.d0)) IFXCN= 1
      NCN2= NCN+2
      IF(ILR.EQ.2) THEN
c** For ILR=2 ,  use   Y = VLIM - CNN/X**NCN - BSR/X**(NCN+2)
c*  If CNN held fixed need ILR > 2  to prevent discontinuity
          IF(LNPT.GT.0) THEN
              IF(IFXCN.LE.0) THEN
                  CNN= ((VLIM-Y1)*X1**NCN2 -
     1                 (VLIM-Y2)*X2**NCN2)/(X1**2-X2**2)
                  ENDIF
              ALR= CNN
              BLR= (VLIM-Y1)*X1**NCN2 - CNN*X1**2
              WRITE(6,620) X2,VLIM,CNN,NCN,BLR,NCN2
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  YY(I)= VLIM-(ALR+BLR/XX(I)**2)/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.EQ.3) THEN
c** For ILR=3 , use   Y = VLIM - (CN + CN2/X**2 + CN4/X**4)/X**NCN
          IF(LNPT.GT.0) THEN
              NCN4= NCN+4
              IF(IFXCN.GT.0) THEN
                  ALR= CNN
                  BLR= (((VLIM-Y1)*X1**NCN-ALR)*X1**4-((VLIM-Y2)
     1                     *X2**NCN-ALR)*X2**4)/(X1**2-X2**2)
                  CLR= ((VLIM-Y1)*X1**NCN-ALR-BLR/X1**2)*X1**4
                ELSE
                  EX1= X1**2
                  EX2= X2**2
                  EX3= X3**2
                  DX1= (VLIM-Y1)*X1**NCN4
                  DX2= (VLIM-Y2)*X2**NCN4
                  DX3= (VLIM-Y3)*X3**NCN4
                  BLR= (DX1-DX2)/(EX1-EX2)
                  ALR= (BLR-(DX2-DX3)/(EX2-EX3))/(EX1-EX3)
                  BLR= BLR-ALR*(EX1+EX2)
                  CLR= DX1-(ALR*EX1+BLR)*EX1
                ENDIF
              WRITE(6,622) X2,VLIM,ALR,NCN,BLR,NCN2,CLR,NCN4
              ENDIF
          IF(MBEG.LE.NPP) THEN
              DO  I= MBEG,NPP
                  EX2= 1.d0/XX(I)**2
                  YY(I)= VLIM-(ALR+EX2*(BLR+EX2*CLR))/XX(I)**NCN
                  ENDDO
              ENDIF
          GO TO 90
          ENDIF
      IF(ILR.GE.4) THEN
c** For ILR.ge.4,   Y = VLIM-SUM(BB(K)/X**K) , (K=NCN,NMX=NCN+ILR-1)
          IF(LNPT.GT.0) THEN
              IF(NCN.LE.0) NCN= 1
              IMX1= ILR-1
              NMX= NCN+IMX1
              JR2= 0
              IF(IR2.GT.0) JR2= 2
              IDER= 0
              JMAX= ILR
              IF(IFXCN.GT.0) JMAX= IMX1
              WRITE(6,624) X2,ILR,NCN,VLIM
              IF(IFXCN.GT.0) WRITE(6,626) NCN,CNN
              ENDIF
c** Actually extrapolate with polynomial fitted to the last JMAX
c  values of  (VLIM - YI(I))*XI(I)**NMX  , & then convert back to  YY(I).
          IF(MBEG.LE.NPP) THEN
              J= NTP- JMAX
              DO  I= 1,JMAX
                  J= J+1
                  XJ(I)= XI(J)
                  YJ(I)= (VLIM-YI(J)/XI(J)**JR2)*XI(J)**NMX
                  IF(IFXCN.GT.0) YJ(I)= YJ(I)- CNN*XI(J)**IMX1
                  ENDDO
              DO  I= MBEG,NPP
                  CALL PLYINTRP(XJ,YJ,JMAX,XX(I),DUMM,JMAX,IDER)
                  YY(I)= DUMM(1)
                  IF(IFXCN.GT.0) YY(I)= YY(I)+ CNN*XX(I)**IMX1
                  YY(I)= VLIM-YY(I)/XX(I)**NMX
                  ENDDO
              ENDIF
          ENDIF
c** Finished extrapolation section.
   90 CONTINUE
c** Test smoothness at outer join to analytic extrapolation function
c     IF((LNPT.GT.0).AND.(MBEG.LE.NPP)) THEN
c         MST= MBEG-5
c         IF(MST.LT.1) MST= 1
c         MFN= MST+8
c         IF(MFN.GT.NPP) MFN= NPP
c         WRITE(6,627) X2,((XX(I),YY(I),I= J,MFN,3),J= MST,MST+2)
c         NPRF= MFN
c         ENDIF
c 627 FORMAT('     Verify smoothness of outer join at   X=',F9.5/
c    1  (3X,3(F10.5,G15.7)))
      RETURN
  612 FORMAT('  *** BUT *** since exponent has positive coefficient, swi
     1tch form ...')
  614 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1  F12.4,' - (',1PD13.6,') * exp{-',0PF10.6,' * (R -',F9.6,')**2}')
  616 FORMAT(' Function for  X .GE.',F8.4,'   generated as'/'   Y=',
     1 F12.4,' - (',1PD13.6,') * R**',0PF10.6,'  * exp{-(',F11.6,'*R)}')
  618 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,SP,1PD15.6,'/X**(',SS,D13.6,')] ,  yielding   NCN=',I3)
  620 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/'   Y=',
     1  F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',SS,I1,']')
  622 FORMAT(' Extrapolate to  X .GE.',F8.4,'  using'/
     1  '   Y=',F12.4,' - [',1PD13.6,'/X**',I1,SP,D14.6,'/X**',
     2  SS,I1,SP,D14.6,'/X**',SS,I2,']')
  624 FORMAT(' Function for  X .GE.',F7.3,'  generated by',I3,
     1 '-point inverse-power interpolation'/'   with leading term  1/R**
     2',I1,'  relative to dissociation limit   YLIM=',F11.3)
  626 FORMAT('   and (dimensionless) leading coefficient fixed as   C',
     1  I1,'=',G15.8)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PLYINTRP(XI,YI,NPT,RR,C,NCFT,IDER)
c* From the NPT known mesh points (XI,YI) ,given in order of increasing
c  or decreasing XI(I), select the NCFT points (XJ,YJ) surrounding the 
c  given point RR, and by fitting an (NCFT-1)-th degree polynomial through
c  them, interpolate to find the function CC(1) and its first IDER 
c  derivatives (CC(I+1),I=1,IDER) evaluated at RR.
c* Adapted by  R.J. Le Roy  from algorithm #416,Comm.A.C.M.;  27/02/1988
c=======================================================================
      INTEGER  I,J,K,I1,I2,IFC,IM,IDER,J1,NH,NPT,NCFT
      REAL*8  RR,XX,XI(NPT),YI(NPT),C(NCFT),XJ(20),YJ(20)
c
      IF((NCFT.GT.20).OR.(NCFT.GT.NPT)) GO TO 101
      NH= NCFT/2
c** First locate the known mesh points (XJ,YJ) bracketing RR
      I1= 1
      I2= NCFT
      IF(NCFT.NE.NPT) THEN
          IF(XI(NPT).LE.XI(1)) THEN
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).LT.RR) GO TO 20
                  ENDDO
            ELSE
              DO  I= 1,NPT
                  IM= I
                  IF(XI(I).GT.RR) GO TO 20
                  ENDDO
            ENDIF
   20     I1= IM-NH
          IF(I1.LE.0) I1= 1
          I2= I1+NCFT-1
          IF(I2.GT.NPT) THEN
              I2= NPT
              I1= I2-NCFT+1
              ENDIF
          ENDIF
      J= 0
      DO  I= I1,I2
          J= J+1
          XJ(J)= XI(I)-RR
          YJ(J)= YI(I)
          ENDDO
c** Now determine polynomial coefficients C(I).
      DO  I= 2,NCFT
          I1= I-1
          K= I1+1
          DO  J= 1,I1
              K= K-1
              YJ(K)= (YJ(K+1)-YJ(K))/(XJ(I)-XJ(K))
              ENDDO
          ENDDO
      C(1)= YJ(1)
      DO  I= 2,NCFT
          XX= XJ(I)
          C(I)= C(I-1)
          IF(I.NE.2) THEN
              I1= I-1
              K= I1+1
              DO  J= 2,I1
                  K= K-1
                  C(K)= -XX*C(K)+C(K-1)
                  ENDDO
              ENDIF
          C(1)= YJ(I)-XX*C(1)
          ENDDO
c** Finally, convert polynomial coefficients to derivatives at RR.
      IFC= 1
      IF(IDER.GE.NCFT) IDER= NCFT-1
      IF(IDER.LE.1) GO TO 99
      DO  I= 2,IDER
          J= I+1
          IFC= IFC*I
          C(J)= C(J)*IFC
          ENDDO
      IF(J.LT.NCFT) THEN
          J1= J+1
          DO  I= J1,NCFT
              C(I)= 0.D+0
              ENDDO
          ENDIF
   99 RETURN
  101 WRITE(6,601) NCFT,NCFT,NPT
      STOP
  601 FORMAT(/' *** Dimensioning ERROR in PLYINTRP :  either   (NCFT=',
     1  I2,' .GT. 20)   or   (NCFT=',I2,' .GT. NPT=',I3,')')
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c**********************************************************************
      SUBROUTINE SPLINT(LNPT,NTP,R1,V1,MBEG,MEND,XX,YY)
c** Subroutine to generate (if LNPT.ge.0) 4*NTP coefficients CSP(J)
c  of a cubic spline passing through the NTP points (R1(J),V1(J))
c  and to then calculate values of the resulting function YY(I) at the
c  entering abscissae values XX(I) for  I=MBEG to MEND.
c** If LNPT < 0 , generate function values at the given XX(I) using
c  the coefficients CSP(J) obtained and SAVEd on a preceding call.
c** Assumes both R1(J) & XX(I) are monotonic increasing.
c+++++ Calls only subroutine SPLINE +++++++++++++++++++++++++++++++++++
c======================================================================
      INTEGER MAXSP
      PARAMETER (MAXSP=6400)
      INTEGER  I,IER,I1ST,IDER,JK,K,KK,LNPT,N2,N3,NIPT,NTP,MBEG,MEND
      REAL*8 EPS,R2,RI,RRR,TTMP,R1(NTP),V1(NTP),CSP(MAXSP),
     1  YY(MEND),XX(MEND)
      SAVE CSP
c
      IF(4*NTP.GT.MAXSP) THEN
          WRITE(6,602) MAXSP,NTP
          STOP
          ENDIF
      EPS= 1.D-6*(R1(2)-R1(1))
      N2= 2*NTP
      N3= 3*NTP
      IF(LNPT.GT.0) THEN
c** On first pass for a given data set, generate spline function
c  coefficients in subroutine SPLINE
c** Start by using a cubic polynomial at each end of the range to get
c  the first derivative at each end for use in defining the spline.
          IDER= 1
          NIPT= 4
          I1ST= NTP-3
          CALL PLYINTRP(R1(I1ST),V1(I1ST),NIPT,R1(NTP),CSP,NIPT,IDER)
          TTMP= CSP(2)
          CALL PLYINTRP(R1,V1,NIPT,R1(1),CSP,NIPT,IDER)
          CSP(1)= CSP(2)
          CSP(2)= TTMP
c** Now call routine to actually generate spline coefficients
          CALL SPLINE(R1,V1,NTP,3,CSP,MAXSP,IER)
          IF(IER .NE. 0) THEN
              WRITE(6,604)
              STOP
              ENDIF
          ENDIF
      IF(MEND.LT.MBEG) GO TO 99
c** Now, use spline to generate function at desired points XX(I)
      DO  I= MBEG,MEND
          RI= XX(I)
          RRR= RI-EPS
          KK= 1
c** For a monotonic increasing distance array XX(I),  this statement 
c  speeds up the search for which set of cubic coefficients to use.
          IF(I.GT.MBEG) THEN
              IF(XX(I).GT.XX(I-1)) KK= JK
              ENDIF
          DO  K= KK,NTP
              JK= K
              IF(R1(K).GE.RRR) GO TO 64
              ENDDO
   64     CONTINUE
          JK= JK-1
          IF(JK.LT.1) JK= 1
          R2= RI-R1(JK)
          YY(I)= CSP(JK)+R2*(CSP(NTP+JK)+R2*(CSP(N2+JK)+R2*CSP(N3+JK)))
          ENDDO
   99 RETURN
  602 FORMAT(' *** ERROR in SPLINT ***  Array dimension  MAXSP=',I4,
     1  ' cannot contain spline coefficients for  NTP=',I4)
  604 FORMAT(' *** ERROR in generating spline coefficients in SPLINE')
      END
c**********************************************************************
      SUBROUTINE SPLINE(X,Y,N,IOPT,C,N4,IER)
c** Subroutine for generating cubic spline coefficients
c  C(J), (J=1,N4=4*N) through the N points X(I), Y(I).
c** C(I+M*N), M=0-3  are the coefficients of order  0-3  of cubic
c  polynomial expanded about X(I) so as to describe the interval:
c             -  X(I) to X(I+1)  , if  X(I)  in increasing order
c             -  X(I-1) to X(I)  , if  X(I)  in decreasing order.
c** IOPT indicates boundary conditions used in creating the  spline .
c*  If (IOPT=0)  second derivatives = zero at both ends of range.
c*  If (IOPT=1)  1st derivative at first point X(1) fixed at C(1),
c                and 2nd derivative at X(N) = zero.
c*  If (IOPT=2)  1st derivative at last point X(N) fixed at C(2),
c                and 2nd derivative at X(1) = zero.
c*  If (IOPT=3)  constrain first derivatives at end points to have
c                (read in) values  C(1)  at  X(1)  &  C(2)  at  X(N)
c** IER is the error flag.  IER=0  on return if routine successful.
c-----------------------------------------------------------------------
      INTEGER I,II,IER,IOH,IOL,IOPT,J,J1,J2,J3,NER,N,N4,JMP
      REAL*8  A,H,R,DY2,DYA,DYB,XB,XC,YA,YB, X(N),Y(N),C(N4)
c
      JMP= 1
      NER= 1000
      IF(N.LE.1) GO TO 250
c** Initialization
      XC= X(1)
      YB= Y(1)
      H= 0.D0
      A= 0.D0
      R= 0.D0
      DYB= 0.D0
      NER= 2000
c
c  IOL=0 - given derivative at firstpoint
c  IOH=0 - given derivative at last point
c
      IOL= IOPT-1
      IOH= IOPT-2
      IF(IOH.EQ.1) THEN
          IOL= 0
          IOH= 0
          ENDIF
      DY2= C(2)
c
c  Form the system of linear equations
c  and eliminate subsequentially
c
      J= 1
      DO  I= 1,N
          J2= N+I
          J3= J2+N
          A= H*(2.D0-A)
          DYA= DYB+H*R
          IF(I.GE.N) THEN
c
c  set derivative dy2 at last point
c
              DYB= DY2
              H= 0.D0
              IF(IOH.EQ.0) GOTO 200
              DYB= DYA
              GOTO 220
              ENDIF
          J= J+JMP
          XB= XC
          XC= X(J)
          H= XC-XB
c
c  II= 0 - increasing abscissae
c  II= 1 - decreasing abscissae
c
          II= 0
          IF(H.LT.0) II= 1
          IF(H.EQ.0) GO TO 250
          YA= YB
          YB= Y(J)
          DYB= (YB-YA)/H
          IF(I.LE.1) THEN
              J1= II
              IF(IOL.NE.0) GO TO 220
              DYA= C(1)
              ENDIF
200       IF(J1.NE.II) GO TO 250
          A= 1.D0/(H+H+A)
220       R= A*(DYB-DYA)
          C(J3)= R
          A= H*A
          C(J2)= A
          C(I)= DYB
          ENDDO
c
c  back substitution of the system of linear equations
c     and computation of the other coefficients
c
      A= 1.D0
      J1= J3+N+II-II*N
      I= N
      DO  IOL= 1,N
          XB= X(J)
          H= XC-XB
          XC= XB
          A= A+H
          YB= R
          R= C(J3)-R*C(J2)
          YA= R+R
          C(J3)= YA+R
          C(J2)= C(I)-H*(YA+YB)
          C(J1)= (YB-R)/A
          C(I)= Y(J)
          A= 0.D0
          J= J-JMP
          I= I-1
          J2= J2-1
          J3= J3-1
          J1= J3+N+II
          ENDDO
      IER= 0
      RETURN
250   IER= NER
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,XO,RM2,VV,
     1                                                        NCN,CNN)
c** Generate analytic potential  VV(i)  as specified by the choice
c  of parameter IPOTL (see comments in PREPOT (& in main program))
c** All potentials generated in units cm-1 with absolute asymptote at
c  (input) energy VLIM for distance array  X0(i) Angstroms.
c** Return with NCN equal to power of asymptotically dominant inverse
c  power term in long range part of potential
c** Born-Oppenheimer correction functions in IPOTL=3 option may have up
c  to NBOB+1 terms.
c-----------------------------------------------------------------------
      INTEGER NBOB
      PARAMETER (NBOB=20)
      INTEGER  I,J,M,IBOB,IAN1,IAN2,IMN1,IMN2,RMN1,RMN2,IORD,IPOTL,
     1  NC1,NC2,NG1,NG2,NCMAX,NPAR,MPAR,NVARB,NPP,LNPT,NCN,GNS,GEL
      CHARACTER*2 NAME1,NAME2
      REAL*8  A0,A1,A2,A3,ALFA,BETA,BINF,B1,B2,CSAV,
     1  ABUND,CNN,DSCM,DX,DX1,FCT,
     2  FC1,FC2,MASS1,MASS2,RMASS1,RMASS2,RC6,RC8,RC10,RC12,RC14,
     3  RCNPAR,RD,RDIF,REQ,RPOW,RX,SC1,SC2,SG1,SG2,VLIM,
     4  XDF,X1,XM,XN,XM2C,XP1,ZZ,ZP, CA1(0:NBOB),CA2(0:NBOB),
     5  GA1(0:NBOB),GA2(0:NBOB),PARM(20),XO(NPP),VV(NPP),RM2(NPP)
c
      SAVE IBOB,IPOTL,IORD,MPAR,NPAR,NVARB,DSCM,REQ,PARM,CA1,CA2,GA1,
     1     GA2,RX,CSAV
c
      IF(LNPT.GT.0) THEN
c** Parameter definitions listed preceeding CALL in subroutine PREPOT
c-----------------------------------------------------------------------
          READ(5,*) IPOTL, MPAR, NPAR, NVARB, IBOB, DSCM, REQ
          IF(IPOTL.EQ.1) NVARB= 0
          IF((IPOTL.EQ.3).AND.(NPAR.EQ.-1)) NVARB= 2
          IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
          IF(IBOB.GT.0) THEN
              READ(5,*) RMN1, RMN2, NC1, NC2, NG1, NG2, RX
c-----------------------------------------------------------------------
              NCMAX= MAX0(NC1,NC2,NG1,NG2)
              IF(NCMAX.LT.0) THEN
                  IBOB= 0
                ELSE
c** If appropriate, read parameters & prepare to add mass-dep. BOB corrn
                  CALL MASSES(IAN1,IMN1,NAME1,GEL,GNS,MASS1,ABUND)
                  CALL MASSES(IAN1,RMN1,NAME1,GEL,GNS,RMASS1,ABUND)
                  CALL MASSES(IAN2,IMN2,NAME2,GEL,GNS,MASS2,ABUND)
                  CALL MASSES(IAN2,RMN2,NAME2,GEL,GNS,RMASS2,ABUND)
                  WRITE(6,628)
c  For simplicity, first zero out all correction function coefficients
                  DO  I=0,NCMAX
                      CA1(I)= 0.d0
                      CA2(I)= 0.d0
                      GA1(I)= 0.d0
                      GA2(I)= 0.d0
                      ENDDO
                  FC1= 0.d0
                  FC2= 0.d0
c=======================================================================
c** Read actual B-O-B polynomial expansion coefficients
c=======================================================================
                  IF(NC1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA1(I), I=0,NC1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 1,MASS1,NC1,NAME1,RMN1,NAME1,
     1                                    IMN1,NC1+1,(CA1(I),I= 0,NC1)
                          FC1= 1.d0 - RMASS1/MASS1
                        ELSE
                          WRITE(6,632) 1,MASS1,NC1,NAME1,IMN1,NC1+1,
     1                                               (CA1(I),I= 0,NC1)
                          FC1= 1.d0/MASS1
                        ENDIF
                      ENDIF
                  IF(NC2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (CA2(I), I=0,NC2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,630) 2,MASS2,NC2,NAME2,RMN2,NAME2,
     1                                    IMN2,NC2+1,(CA2(I),I= 0,NC2)
                          FC2= 1.d0 - RMASS2/MASS2
                        ELSE
                          WRITE(6,632) 2,MASS2,NC2,NAME2,IMN2,NC2+1,
     1                                                (CA2(I),I=0,NC2)
                          FC2= 1.d0/MASS2
                        ENDIF
                      ENDIF
                  IF(NG1.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA1(I), I=0,NG1)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 1,MASS1,NG1,NAME1,RMN1,NAME1,
     1                                    IMN1,NG1+1,(GA1(I),I= 0,NG1)
                        ELSE
                          WRITE(6,636) 1,MASS1,NG1,NAME1,IMN1,RX,
     1                                         NG1+1,(GA1(I),I= 0,NG1)
                        ENDIF                     
                      ENDIF
                  IF(NG2.GE.0) THEN
c-----------------------------------------------------------------------
                      READ(5,*) (GA2(I), I=0,NG2)
c-----------------------------------------------------------------------
                      IF(RX.LE.0.d0) THEN
                          WRITE(6,634) 2,MASS2,NG2,NAME2,RMN2,NAME2,
     1                                    IMN2,NG2+1,(GA2(I),I= 0,NG2)
                        ELSE
                          WRITE(6,636) 2,MASS2,NG2,NAME2,IMN2,RX,
     1                                         NG2+1,(GA2(I),I= 0,NG2)
                        ENDIF
                      ENDIF
                  DO  I=0,NCMAX
                      CA1(I)= CA1(I)*FC1
                      CA2(I)= CA2(I)*FC2
                      IF(RX.LE.0.d0) THEN
                          GA1(I)= GA1(I)*(1.d0-FC1)
                          GA2(I)= GA2(I)*(1.d0-FC2)
                        ELSE
                          GA1(I)= GA1(I)*FC1
                          GA2(I)= GA2(I)*FC2
                        ENDIF
                      ENDDO
                ENDIF
              ENDIF
          ENDIF
      IF(IPOTL.EQ.1) THEN 
c=======================================================================
c** Generate a  Lennard-Jones(MPAR,NPAR)  potential here.
c=======================================================================
          XM= MPAR
          XN= NPAR
          XDF= DSCM/(XM-XN)
          IF(LNPT.GE.0) WRITE(6,600) MPAR,NPAR,DSCM,REQ
          NCN= NPAR
          CNN= XM*XDF*REQ**NPAR
          DO  I= 1,NPP
              VV(I)= (XN*(REQ/XO(I))**MPAR - XM*(REQ/XO(I))**NPAR)*XDF
     1                  +VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.2) THEN
c=======================================================================
c** Generate an MLJ potential [as per JCP 112, 3949 (2000)] here ...
c=======================================================================
          NCN= NPAR
          IORD= NVARB-1
          IF(MPAR.LT.0) THEN
c  If appropriate, prepare to calculate switching function
              IORD= IORD- 3
              CNN= PARM(NVARB-2)
              BINF= DLOG(2.d0*DSCM*REQ**NPAR/CNN)
              ALFA= PARM(NVARB- 1)
              RX= PARM(NVARB)
            ELSE
c  Generate limiting Cn value for non-switching case from BINF
              BINF= 0.D0
              DO  I= 1,NVARB
                  BINF= BINF+ PARM(I)
                  ENDDO
              CNN= 2.d0*DSCM*REQ**NPAR *DEXP(-BINF)
            ENDIF
          IF(LNPT.GT.0) THEN
              WRITE(6,602) DSCM,REQ,IORD,(PARM(J),J= 1,IORD+1)
              IF(MPAR.LT.0) WRITE(6,604) NPAR,NPAR,NPAR,PARM(NVARB-2),
     1                                       PARM(NVARB-1),PARM(NVARB)
              ENDIF
          NCN= NPAR
c  Loop over distance array XO(I)
          DO  I= 1,NPP
              ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
              BETA= 0.d0
              DO  J= IORD,0,-1
                  BETA= BETA*ZZ+ PARM(J+1)
                  ENDDO
c  Calculate and apply switching function to MLJ exponent coefficient
              IF(MPAR.LT.0) BETA= BINF+ (BETA- BINF)/
     1                                 (1.d0+ DEXP(ALFA*(XO(I) - RX)))
              VV(I)= DSCM*(1.d0 - (REQ/XO(I))**NPAR *DEXP(-BETA*ZZ))**2
     1                                                    - DSCM+ VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.3) THEN
c=======================================================================
c** Generate a simple Morse, or Extended (EMO) Morse potential, or as
c  special cases, Coxon's GMO or Wei Hua's generalized Morse
c=======================================================================
          BETA= PARM(1)
          NCN= 99
          IF(LNPT.GE.0) THEN
              IF(MPAR.EQ.-1) THEN
c** Option to generate Wei Hua's extended 4-parameter Morse-type potl.
                  CSAV= PARM(2)
                  WRITE(6,605) DSCM,REQ,CSAV,BETA
                ELSE 
                  IF(NVARB.LE.1) WRITE(6,606) DSCM,REQ,BETA
                  IF(NVARB.GT.1) THEN
                      IF(MPAR.GE.0) WRITE(6,608) DSCM,REQ,NVARB-1,
     1                                      NVARB,(PARM(i),i= 1,NVARB)
                      IF(MPAR.EQ.-2) WRITE(6,610) DSCM,REQ,NVARB-1,
     1                                            (PARM(i),i= 1,NVARB)
                      ENDIF
                ENDIF
              ENDIF
c  Loop over distance array XO(I)
          DO  I= 1,NPP
c ... for Wei Hua's extended Morse function ...
              IF(MPAR.EQ.-1) THEN
                  VV(I)= DSCM*((1.d0 - DEXP(-BETA*(XO(I)-REQ)))/(1.d0 
     1                - CSAV*DEXP(-BETA*(XO(I)-REQ))))**2 - DSCM+ VLIM
                ELSE 
                  IF(NVARB.GT.1) THEN
                      ZZ= (XO(I)- REQ)/(XO(I)+ REQ)
c ... for Coxon-Hajigeorgiou "GMO" potential
                      IF(MPAR.EQ.-2) ZZ= (XO(I)- REQ)
                      BETA= 0.d0
                      DO  J= NVARB,1,-1             
                          BETA= BETA*ZZ+ PARM(J)
                          ENDDO
                      ENDIF
                  VV(I)= DSCM*(1.d0 - DEXP(-BETA*(XO(I)-REQ)))**2 
     1                                                    - DSCM+ VLIM
                ENDIF
              ENDDO
          ENDIF
      IF(IPOTL.EQ.4) THEN
c=======================================================================
c** Generate Seto-modified form of Surkus' GPEF function which includes
c  Dunham, SPF and OT forms as special cases.
c=======================================================================
          A0= DSCM
          IORD= NVARB-2
          X1= 1.d0
          FCT= PARM(NVARB-1)
          IF(DABS(FCT).GT.0.d0) THEN
              FCT= 1.d0/PARM(NVARB-1)
              DO  J=1,IORD
                  X1= X1+ PARM(J)*FCT**J
                  ENDDO
              DSCM= DSCM*X1
              ENDIF
          IF(NPAR.EQ.1) THEN
c  Cases with power =1 (including Dunham, SPF & O-T expansions).
              IF(DABS(PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Dunham expansion ...
                  WRITE(6,612) PARM(NVARB),REQ,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= -99
                  CNN= 0.d0
                  ENDIF
              IF(DABS(PARM(NVARB)).LE.0.d0) THEN
c ... print for Simons-Parr-Finlan expansion ...
                  WRITE(6,614) PARM(NVARB-1),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              IF(DABS(PARM(NVARB)-PARM(NVARB-1)).LE.0.d0) THEN
c ... print for Ogilvie-Tipping expansion ...
                  WRITE(6,616) PARM(NVARB),REQ,DSCM,A0,NVARB-2,
     1                                          (PARM(I),I= 1,NVARB-2)
                  NCN= 1
                  ENDIF
              ENDIF
          IF((NPAR.NE.1).OR.((DABS(PARM(NVARB)-PARM(NVARB-1)).GT.0.d0)
     1           .AND.(DABS(PARM(NVARB)*PARM(NVARB-1)).GT.0.d0))) THEN
c ... print for general GPEF expansion variable ...
              IF(NPAR.LT.0) THEN
c ... for negative NPAR, convert to equivalent positive NPAR case
                  NPAR= -NPAR
                  A1= PARM(NVARB)
                  PARM(NVARB)= -PARM(NVARB-1)
                  PARM(NVARB-1)= -A1
                  ENDIF
              WRITE(6,618) NPAR,NPAR,PARM(NVARB-1),NPAR,PARM(NVARB),
     1                 NPAR,REQ,DSCM,A0,NVARB-2,(PARM(I),I= 1,NVARB-2)
              NCN= NPAR
              ENDIF
          DO  I= 1, NPP
              ZZ= (XO(I)**NPAR - REQ**NPAR)/(PARM(NVARB-1)*XO(I)**NPAR
     1                                       + PARM(NVARB)*REQ**NPAR)
              A1= 1.d0
              ZP= 1.d0
              DO  J=1, NVARB-2
                  ZP= ZP*ZZ
                  A1= A1+ PARM(J)*ZP
                  ENDDO
              VV(I)= A0*ZZ*ZZ*A1 + VLIM
              ENDDO
          ENDIF
      IF(IPOTL.EQ.5) THEN
c=======================================================================
c** For generalized  H.F.D.(NPAR,6,8,10,12,14)  potential with reduced 
c  form   VBAR = ALFA*x**PARM(5) * exp[-BETR*x - PARM(4)*x**2] - D(x)*
c       [PARM(6)/x**NPAR + PARM(7)/x**6 + PARM(8)/x**8 + PARM(9)/x**10 
c       + PARM(10)/X**12 + PARM(11)/X**14]   where   x=r/R_e ,  
c  VBAR= V/epsilon   and   D(x)= exp[-PARM(1)*(PARM(2)/x - 1)**PARM(3)]
c  for  x < PARM(2)
c=======================================================================
          A1= PARM(1)
          A2= PARM(2)
          A3= PARM(3)
          B2= PARM(4)
          RC8= 0.d0
          RC10= 0.d0
          RC12= 0.d0
          RC14= 0.d0
          RCNPAR= PARM(6)
          NCN= 6
          IF(RCNPAR.GT.0.d0) NCN= NPAR
          RC6= PARM(7)
          IF(NVARB.ge.8)  RC8= PARM(8)
          IF(NVARB.ge.9)  RC10= PARM(9)
          IF(NVARB.ge.10) RC12= PARM(10)
          IF(NVARB.ge.11) RC14= PARM(11)
          DX= 1.d0
          DX1= 0.d0
          IF(A2.GT.1.d0) THEN
              DX= DEXP(-A1*(A2- 1.d0)**A3)
              DX1= A1*A2*A3*DX*(A2- 1.d0)**(A3- 1.d0)
              ENDIF
          ALFA= -1.D0+ (RCNPAR+ RC6+ RC8+ RC10+ RC12+ RC14)*DX
          IF(ALFA.LE.0.d0) THEN
              WRITE(6,622) RCNPAR,RC6,RC8,RC10,RC12,RC14,ALFA
              STOP
              ENDIF
          B1= ((NPAR*RCNPAR+6.D0*RC6+8.D0*RC8+10.D0*RC10+12.d0*RC12+
     1      14.d0*RC14)*DX - (RCNPAR+RC6+RC8+RC10+RC12+RC14)*DX1)/ALFA
     2      + PARM(5) - 2.D0*B2
          ALFA= ALFA*DEXP(B1+B2)
          IF(LNPT.GE.0) WRITE(6,624) NPAR,PARM(5),B1,B2,ALFA*DSCM,
     1                 RCNPAR,RC6,RC8,RC10,RC12,RC14,A1,A2,A3,DSCM,REQ
          DO  I= 1,NPP
              X1= XO(I)/REQ
              XP1= 0.0D0
              IF((B1*X1+ B2*X1**2).LT.170.D0) XP1= DEXP(-X1*(B1+ B2*X1))
              XP1= XP1*X1**PARM(5)
              FC1= 1.D0
              IF(A2.GT.X1) FC1= DEXP(-A1*(A2/X1- 1.d0)**A3)
              XM2C= (REQ/XO(I))**2
              VV(I)= DSCM*(ALFA*XP1- FC1*(((((RC14*XM2C+RC12)*XM2C+RC10)
     1         *XM2C+ RC8)*XM2C+ RC6)*XM2C**3 + RCNPAR/X1**NPAR)) + VLIM
              ENDDO
          ENDIF
      IF(IBOB.GT.0) THEN
c=======================================================================
c** If appropriate, generate Born-Oppenheimer breakdown correction 
c  functions to rotationless and/or centrifugal potential(s).
c [Special "Coxon" option: if  RX > 0.0, expand as per older Coxon work]
c=======================================================================
          IF(RX.GE.0.D0) RDIF= REQ-RX
          DO  I=1,NPP
              IF(RX.LE.0.d0) THEN
                  ZZ= (XO(I)-REQ)/(XO(I)+REQ)
                ELSE
                  ZZ= XO(I)- REQ
                  RD= XO(I)- RX
                ENDIF
              SC1= 0.d0
              SC2= 0.d0
              SG1= 0.d0
              SG2= 0.d0
              RPOW= 1.d0
              DO  J= 0,NCMAX
                  SC1= SC1+ RPOW*CA1(J)
                  SC2= SC2+ RPOW*CA2(J)
                  IF(RX.LE.0.d0) THEN
                      SG1= SG1+ RPOW*GA1(J)
                      SG2= SG2+ RPOW*GA2(J)
                    ELSE 
                      M= J-1
                      SG1= SG1+ (RD**J -RDIF**J)*GA1(J)
                      SG2= SG2+ (RD**J -RDIF**J)*GA2(J)
                    ENDIF
                  RPOW= RPOW*ZZ
                  ENDDO
              RM2(I)= (1.d0+ SG1+ SG2)/XO(i)**2
              VV(I)= VV(I) + SC1 + SC2
              ENDDO
          ENDIF
      RETURN
  600 FORMAT(/' Lennard-Jones(',I2,',',I2,') potential with   De=',
     1  F10.3,'(cm-1)   Re =',F10.6,'(A)')
  602 FORMAT(/' Use an MLJ potential with   De =',F10.3,
     1  '(cm-1)    Re =',F12.8,'(A)'/3x,'with exponent parameter BETA an
     2 order-',i2,' polynomial in  z=(R-Re)/(R+Re)  with'/
     3  '   coefficients:',1PD16.8,3D16.8:/(5D16.8:))
  604 FORMAT(' & exponent switching function yielding limiting C',i1,
     1 '/R^',i1,' with   C_',i1,'=',1PD13.6/10x,'defined by   ALPHA_s=',
     2  0Pf9.6,'   R_s=',f10.6)
  605 FORMAT(/' Potential is a Hua-Wei 4-parameter Morse type function w
     1ith   De =',F11.4/11x,'Re =',F12.9,'   C=',f7.4,'   &   beta=',
     1  F13.10,' [1/Angstroms]')
  606 FORMAT(/' Potential is a simple Morse function with   De =',F11.4,
     1  '    Re =',F12.9/39x,'and   beta =',F13.10,' [1/Angstroms]')
  608 FORMAT(/' Potential is Extended Morse Oscillator with   De=',
     1  F11.4,'    Re=',F12.9/5x,'Exponent factor "beta" is order-',i2,
     2 ' power series in  z=(R-Re)/(R+Re)  with'/5x,I2,' coefficients:',
     3  1x,1PD18.9,2D18.9:/(4X,4D18.9:))
  610 FORMAT(/' Potential is Generalized Morse Oscillator with   De=',
     1 F10.3,'   Re=',F11.8/4x,'Exponent factor "beta" is',i3,' order po
     2wer series in (R-Re) with coefficients:'/4x,1PD18.9,3D18.9:/
     3 (4X,4D18.9:))
  612 FORMAT(/' Potential is a Dunham expansion in  (R-Re)/(',f5.2,
     1  ' * Re)  with   Re=',f12.9/'  a0=',1PD16.9,'   and',i3,
     2  '  a_i coefficients:'/(5D16.8))
  614 FORMAT(/' Potential is an SPF expansion in  (R-Re)/(',F5.2,
     1  '* R)  with   Re=',f12.9/5x,'De=',g18.10,'   b0=',
     2  1PD16.9,'   and',i3,'  b_i  coefficients:'/(5D16.8))
  616 FORMAT(/' Potential is an O-T expansion in  (R-Re)/[',f5.2,
     1  '*(R+Re)]  with   Re=',f12.9/5x,'De=',G18.10,
     2  '   c0=',1PD16.9,'   and',i3,'  c_i coefficients:'/(5D16.8))
  618 FORMAT(/' Potential is a general GPEF expansion in  (R^',i1,
     1  ' - Re^',i1,')/(',SP,F5.2,'*R^',SS,i1,SP,F6.2,'*Re^',SS,i1,')'/
     2  5x,'with   Re=',f12.9,'   De=',g18.10,'   g0=',1PD16.9/
     2  5x,'and',i3,'  g_i coefficients:  ',3D16.8/(5D16.8:))
  622 FORMAT(/' *** ERROR in generating HFD potential *** C',i1,
     1 ', C6, C8, C10, C12, C14  =',6G15.7/10X,'yield    ALFA =',G15.7)
  624 FORMAT(/' Potential is Generalized HFD(',i1,',6,8,10,12,14) with',
     1 '  gamma=',f9.6/'    beta1=',f12.8,'    beta2=',f9.6,'    A=',
     2 1PD16.9/"    reduced {Cn's}:",3D14.6/19x,3D14.6/'    Damping func
     3tion  D(R)= exp[ -',0Pf6.4,'*(',f7.4,'/X -1.0)**',f5.2,']' /
     4 '  & DSCM=',f10.4,'[cm-1]   Re=',f9.6,'[Angst.]')
  628 FORMAT(' ')
  630 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1  I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-R
     2e)/(R+Re)}] * [1- MASS(',A2,i3,')/MASS(',A2,I3,')]'/5x,'with',i3,
     3  ' coefficients:',3G18.10:/(8x,4G18.10:))
  632 FORMAT(' B-O-B correction to rotationless potential for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in (R-Re)
     2]/[MASS(',A2,I3,')]   with',i3,' coefficients:'/(5x,4G18.10:))
  634 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Re
     2/(R+Re)}] * [MASS(',A2,I3,')/MASS(',A2,I3,')]'/5x,'with',i3,' coef
     3ficients:',3G18.10:/(8x,4G18.10:))
  636 FORMAT(' B-O-Breakdown correction to centrifugal term for atom-',
     1 I1,'  of mass ',f14.10/5x,'is  [order-',I2,' polynomial in {(R-Rx
     2)**i - (Re-Rx)**i}] / [MASS(',A2,I3,')]'/5x,'with  Rx=',f6.3,
     3 '  &',i3,' coefficients: ',1PD18.9,D18.9:/(5x,4D18.9:))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

