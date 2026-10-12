c***********************************************************************
      SUBROUTINE FUNUNC(NSTATES,OSEL,NPPS,PU,CM,FU)
c***********************************************************************
c** This subroutine will calculate the radially dependent uncertainties
c   for V(r), B(r), UA(r), UB(r), qA(r), qB(r) and lambdaR(r) and print 
c   them out in tecplot format.
c** This subroutine calls the subroutine MMCALC.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c++       COPYRIGHT 1997  by  Jenning Seto       (ver of 08/10/1998) +++
c         Modified by Yiye Huang
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the authors.    +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c    OSEL
c    NSTATES
c    NPPS(s)
c    PU(n)
c    CM(n,n)
c
c** On exit:
c    FU
c
      INCLUDE 'arrsizes.h'

      INTEGER NSTATES, OSEL(8,NSTATEMX), NPPS(NSTATEMX),
     1        ISTATE, IVAL, IPV, I, PI(NPARMX), FSEL

      REAL*8 FU, CM(NPARMX,NPARMX), PU(NPARMX), RDVAL, PT(NPARMX)
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     4  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR

      REAL*8 Rd,Vd,Vsr(NPNTMX,NSTATEMX),bM,Bsr(NPNTMX,NSTATEMX)
      INTEGER nPointSR(NSTATEMX)
      COMMON /VsrBLK/Vsr,Bsr,nPointSR

      REAL*8 Plr(NPNTMX,NSTATEMX),Blr(NPNTMX,NSTATEMX)
      INTEGER nPointLR(NSTATEMX)
      COMMON /PlrBLK/Plr,Blr,nPointLR

      REAL*8 RMAX(NSTATEMX)
      COMMON /RMAXBLK/RMAX
c
c** Born-Oppenheimer Breakdown Functions.
c
      REAL*8 UAR(NPNTMX,NSTATEMX),UBR(NPNTMX,NSTATEMX),
     1  QAR(NPNTMX,NSTATEMX),QBR(NPNTMX,NSTATEMX),
     2  ZMUA(NISTPMX),ZMUB(NISTPMX),ZMQA(NISTPMX),ZMQB(NISTPMX)

      COMMON /BFBLK/UAR,UBR,QAR,QBR,ZMUA,ZMUB,ZMQA,ZMQB
**
** Lamda doubling function
**
      INTEGER nLamdCoef(NSTATEMX),ifLamdCoef(0:NRNAMX,NSTATEMX)
      REAL*8 lamdCoef(0:NRNAMX,NSTATEMX),lambdaR(NPNTMX,NSTATEMX)

      COMMON /lambdaBLK/nLamdCoef,ifLamdCoef,lamdCoef,lambdaR
c
c** Partial derivative arrays.
c
      REAL*8 DVBODP(NPNTMX,NPOTMX),DUADP(NPNTMX,NPOTMX),
     1  DUBDP(NPNTMX,NPOTMX),DQADP(NPNTMX,NPOTMX),DQBDP(NPNTMX,NPOTMX)
     2  ,dLdP(NPNTMX,NPOTMX),DBDP(NPNTMX,NPOTMX)

      COMMON /PDVBLK/DVBODP,DUADP,DUBDP,DQADP,DQBDP,DBDP
      COMMON /PdVLBLK/dLdP
      integer il
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DO ISTATE=1,NSTATES
        DO I = 1,NPARMX
          PI(I) = 0
          PT(I) = 0.0d0
        END DO
        IPV = 0
        DO I = 1,ISTATE-1
          IPV = IPV + NPPS(I)
        END DO
        DO I = IPV+1,IPV+NPPS(ISTATE)
          PI(I) = 1
        END DO
        DO FSEL = 1,7
          IF (OSEL(FSEL,ISTATE).GT.0) THEN
            IF (FSEL.EQ.1) THEN
              WRITE(10,900) 'V(r) ', 'V(r) ', ISTATE
              DO IVAL = 1,nPointSR(ISTATE),OSEL(1,ISTATE)
                RDVAL = RH(ISTATE)*DBLE(IVAL-1) 
                WRITE(10,909) RDVAL,Vsr(IVAL,ISTATE)
              END DO   
            ELSEIF (FSEL.EQ.2) THEN
              WRITE(11,900) 'B(r) ', 'B(r) ', ISTATE
              DO IVAL = 1,nPointSR(ISTATE),OSEL(2,ISTATE)
                RDVAL = RH(ISTATE)*DBLE(IVAL-1)
                WRITE(11,909) RDVAL,Bsr(IVAL,ISTATE)
              END DO
            ELSEIF (FSEL.EQ.3) THEN
              WRITE(12,900) 'UA(r)', 'UA(r)', ISTATE
            ELSEIF (FSEL.EQ.4) THEN
              WRITE(13,900) 'UB(r)', 'UB(r)', ISTATE
            ELSEIF (FSEL.EQ.5) THEN
              WRITE(14,900) 'qA(r)', 'qA(r)', ISTATE
            ELSEIF (FSEL.EQ.6) THEN
              WRITE(15,900) 'qB(r)', 'qB(r)', ISTATE
            ELSEIF (FSEL.EQ.7) THEN
              WRITE(16,902) 'lambda(r)', 'lambda(r)', ISTATE
            END IF
            DO IVAL = 1,NDATPT(ISTATE),OSEL(FSEL,ISTATE)
              RDVAL = RMIN(ISTATE) + RH(ISTATE) * DBLE(IVAL-1)
              DO I = IPV+1,IPV+NPPS(ISTATE)
                IF (FSEL.EQ.1) THEN
                  PT(I) = PU(I) * DVBODP(IVAL,I)
                ELSEIF (FSEL.EQ.2) THEN
                  PT(I) = PU(I) * DBDP(IVAL,I)
                ELSEIF (FSEL.EQ.3) THEN
                  PT(I) = PU(I) * DUADP(IVAL,I)
                ELSEIF (FSEL.EQ.4) THEN
                  PT(I) = PU(I) * DUBDP(IVAL,I)
                ELSEIF (FSEL.EQ.5) THEN
                  PT(I) = PU(I) * DQADP(IVAL,I)*RDVAL**2
                ELSEIF (FSEL.EQ.6) THEN
                  PT(I) = PU(I) * DQBDP(IVAL,I)*RDVAL**2
                ELSEIF (FSEL.EQ.7) THEN
                  PT(I) = PU(I) * dLdP(IVAL,I)*RDVAL**2
                END IF
              END DO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now call subroutine MMCALC to deterimine the matrix multiplication
c   of the correlation matrix CM with the vector PT.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
              CALL MMCALC(IPV+NPPS(ISTATE),NPARMX,PI,PT,CM,FU)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Now write out the function and the uncertainties.
c
              IF (FSEL.EQ.1) THEN
                WRITE(10,910) RDVAL,VPOT(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.2) THEN
                WRITE(11,910) RDVAL,BMORSE(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.3) THEN
                WRITE(12,910) RDVAL,UAR(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.4) THEN
                WRITE(13,910) RDVAL,UBR(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.5) THEN
                WRITE(14,910) RDVAL,QAR(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.6) THEN
                WRITE(15,910) RDVAL,QBR(IVAL,ISTATE),FU
              ELSEIF (FSEL.EQ.7) THEN
                WRITE(16,910) RDVAL,lambdaR(IVAL,ISTATE),FU
              END IF
            END DO
            IF (FSEL.EQ.1) THEN
              DO IVAL = 1,nPointLR(ISTATE)
                RDVAL=RMAX(ISTATE)+RH(ISTATE)*DBLE(IVAL*OSEL(1,ISTATE))
                WRITE(10,909) RDVAL,Plr(IVAL,ISTATE)
              END DO
            ELSEIF (FSEL.EQ.2) THEN
              DO IVAL = 1,nPointLR(ISTATE)
                RDVAL=RMAX(ISTATE)+RH(ISTATE)*DBLE(IVAL*OSEL(1,ISTATE))
                WRITE(11,909) RDVAL,Blr(IVAL,ISTATE)
              END DO
            END IF
          END IF 
        END DO
      END DO
      RETURN
c-----------------------------------------------------------------------
  900 FORMAT(/'variables = "r", "',A5,'", "Uncertainty"'/
     1'zone T = " Plot of ',A5,' for state',I3,'"',/)
  902 FORMAT(/'variables = "r", "',A9,'", "Uncertainty"'/
     1'zone T = " Plot of ',A9,' for state',I3,'"',/)
  909 FORMAT(5X,G18.8,5X,G18.8,10X,'0.0')
  910 FORMAT(5X,G18.8,5X,G18.8,5X,G18.8)
      END
c***********************************************************************
      SUBROUTINE MMCALC(NPTOT,NPMAX,PI,PT,CM,FU)
c***********************************************************************
c** This subroutine will calculate the radial uncertainties for a 
c   statistically calculated function given the correlation matrix and
c   partial derivatives and parameter uncertaintes.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c                  COPYRIGHT 1998  by  Jenning Seto                    +
c   Dept. of Chemistry, Univ. of Waterloo, Waterloo, Ontario, Canada   +
c    This software may not be sold or any other commercial use made    +
c      of it without the express written permission of the author.     +
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** This calculation involves a straightforward matrix multiplication of
c   the vector (PT) which is formed by the partial derivatives of the
c   function in question with the parameters used in the fit (dF/dPk)
c   and the uncertainties associated with those parameters (uPk) with
c   the correlation matrix (CM): 
c
c                          PT(k) = (dF/dPk) * (uPk)
c
c                            FU = (PT)T * CM * PT
c
c** On entry:
c    NPTOT   is the actual number of elements used in the parameter
c            indicies vector (PI), the parameter uncertainty vector
c            (PT), and the correlation matrix (CM).
c    NPMAX   is the maximum array sizes for PI, PT and CM.
c    PI(n)   are indicies for which the elements in PT are not zero.
c    PT(n)   is the vector of the partial derivatives of the function
c            w.r.t. the parameters defining that function times the
c            uncertainties in those parameters.
c    CM(n,n) is the correlation matrix obtained by normalizing
c            variance/covariance matrix:
c                      CM(i,j) = VM(i,j)/SQRT[VM(i,i)*VM(j,j)]
c
c** On Exit:   
c    FU      is the uncertainty calculated.
c
c** Externally dimension:  PI and PT to have size NPMAX
c             CM   as a square matrix with column & row length  NPMAX
c***********************************************************************
      INTEGER NPTOT, NPMAX, PI(NPMAX), I, J

      REAL*8  PT(NPMAX), CM(NPMAX,NPMAX), TVAL, FU
c=======================================================================
c** Initialize the variables.
c
      FU = 0.0d0
c
c** Begin by determining which parts of the arrays will be used.
c
      DO I = 1,NPTOT
        IF (PI(I).EQ.1) THEN
          TVAL = 0.0d0
c
c** Multiply the vector PT with the correlation matrix (CM).
c
          DO J = 1,NPTOT
            IF (PI(J).EQ.1) TVAL = TVAL + PT(J)*CM(I,J)
          END DO
c
c** Now to multiply the vector (PT*CM) with the vector PT.
c
          FU = FU + TVAL * PT(I)
        END IF
      END DO
c
c** Complete calculation of the uncertainty of the function.
c
      FU = DSQRT(FU)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
