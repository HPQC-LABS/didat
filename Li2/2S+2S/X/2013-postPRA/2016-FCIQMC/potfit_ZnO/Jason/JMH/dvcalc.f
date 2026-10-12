c***********************************************************************
      SUBROUTINE DVCALC(NSTATES,OSEL)
c***********************************************************************
c** This subroutine will calculate the first five partial derivatives of
c   the potential function V(R) for ISTATEth electronic state and print
c   them out to channel-31 to 35.
c
      INCLUDE 'arrsizes.h'

      INTEGER NSTATES, OSEL(8,NSTATEMX), ISTATE, ESEL, I, J

      REAL*8 DEA, RDVAL, FAR, BOS(0:NBETAMX), GVAL, BTEMP,
     1  PARTV(NPNTMX,0:5), DMVDRM
c
c** Born-Oppenheimer Potential Variables.
c
      INTEGER PSEL(NSTATEMX),IOMEG(NSTATEMX),NBETA(NSTATEMX),
     1  IFDE(NSTATEMX),IFRE(NSTATEMX),IFBETA(0:NBETAMX,NSTATEMX),
     2  NDATPT(NSTATEMX),NCN(NSTATEMX),KGPEF(NSTATEMX),IRSR(NSTATEMX)

      REAL*8 DE(NSTATEMX),RE(NSTATEMX),BETA(0:NBETAMX,NSTATEMX),
     1  BMORSE(NPNTMX,NSTATEMX),RMIN(NSTATEMX),
     2  RH(NSTATEMX),VPOT(NPNTMX,NSTATEMX),VLIM(NSTATEMX),
     3  EPS(NSTATEMX),
     4  AGPEF(NSTATEMX),BGPEF(NSTATEMX),
     5  CNVAL(NSTATEMX),ALPHA(NSTATEMX),RX(NSTATEMX)

      COMMON /POTBLK/DE,RE,BETA,BMORSE,RMIN,RH,VPOT,VLIM,
     1  EPS,AGPEF,BGPEF,CNVAL,ALPHA,RX,PSEL,IOMEG,NBETA,
     2  IFDE,IFRE,IFBETA,NDATPT,NCN,KGPEF,IRSR
c
c** Initializing variables.
c
      DO ISTATE = 1,NSTATES
        IF (OSEL(8,ISTATE).GT.0) THEN
          DO I = 0,NBETA(ISTATE)
            BOS(I) = BETA(I,ISTATE)
          END DO
c
c** Normalizing the Dissociation Energy for the MMO.
c
          IF (PSEL(ISTATE).EQ.3) THEN
            BTEMP = 0.0d0
            DO I = 0,NBETA(ISTATE)
              BTEMP = BTEMP + BETA(I,ISTATE)
            END DO
            DEA = DE(ISTATE)/DEXP(BTEMP)
          ELSE
            DEA = DE(ISTATE)
          END IF
c
c** Now begin loop over the radial mesh.
c
          DO I = 1,NDATPT(ISTATE),OSEL(8,ISTATE)
            RDVAL = RMIN(ISTATE) + RH(ISTATE) * DBLE(I-1)
            IF (PSEL(ISTATE).EQ.2) THEN
              ESEL = 1
              GVAL = (RDVAL - RE(ISTATE))
            ELSE IF ((PSEL(ISTATE).EQ.3).OR.(PSEL(ISTATE).EQ.4)) THEN
              ESEL = 2
              GVAL = (RDVAL-RE(ISTATE))/(RDVAL+RE(ISTATE))
            ELSE
c
c** If the potential is not the EMO, MMO, or MLJ, then print a warning
c   and return to calling program.
c
              WRITE(6,600) PSEL(ISTATE)
              RETURN
            END IF
c
c** Now calculate the general Gruebele function FAR.
c
            FAR = DEXP(-GVAL*BMORSE(I,ISTATE))
            IF (PSEL(ISTATE).EQ.4) THEN
              FAR = FAR*(RE(ISTATE)/RDVAL)**NCN(ISTATE)
            END IF
c
c** Now to calculate the partial derivatives.
c
            DO J = 1,5
              PARTV(I,J) = DMVDRM(J,RE(ISTATE),RDVAL,
     1        NBETA(ISTATE),BOS,NCN(ISTATE),ESEL,DEA,FAR,VLIM(ISTATE))
            END DO
          END DO
c
c** Now to print out the partial derivatives to channel-31 to 35.
c
          DO J = 1,5
            WRITE(30+J,3000) J, J, ISTATE
            DO I = 1,NDATPT(ISTATE),OSEL(8,ISTATE)
              WRITE(30+J,3010) RDVAL, PARTV(I,J)
            END DO
          END DO
        END IF
      END DO
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** DVCALC ERROR ***',/4x,'Cannot calculate the partial
     1derivatives for potentials of type',I3,'.')
 3000 FORMAT('variables = "r", "PV(',I1,')"',/,'zone T = " Partial deriv
     1ative (',I1,') of state ',I1,' of V(R)."'/)
 3010 FORMAT(5X,G18.8,5X,G18.8)
      END
c***********************************************************************
      REAL*8 FUNCTION DMVDRM(M,RE,R,NBETA,BETA,NCN,ESEL,DE,FAR,VLIM)
c***********************************************************************
c** This function will calculate the Mth partial derivative of V(R).
c
c** On entry:
c    M       is the order of the derivative.
c    RE      is the eqilibrium bond length.
c    R       is the intermolecular distance.
c    NBETA   is the number of beta parameters in the BETA expansion.
c    BETA(n) is the array of beta parameters.
c    NCN     is the power of the radial term.
c    ESEL    is a flag denoting which ALPHA expansion is to be used:
c            = 1 : use the radial expansion - (R-RE)
c            = 2 : use the Ogilvie-Tipping expansion - Z=(R-RE)/(R+RE)
c    DE      is the normalized dissociation energy of the potential.
c    FAR     is the exponent term.
c    VLIM    is the long-range potential limit.
c
c** On exit:
c    DMVDRM  is the Mth derivative of V(R) with respect to R.
c
      INTEGER M, NBETA, NCN, ESEL

      REAL*8 RE, R, BETA(0:NBETA), DE, FAR, VLIM, DMADRM, TOTAL

      IF (M.EQ.0) THEN
        DMVDRM = DE*(1.0d0 - FAR)**2-DE+VLIM
      ELSE IF (M.EQ.1) THEN
        DMVDRM = 2.0d0*DE*FAR*(1.0d0 - FAR)
     1           *DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)
      ELSE IF (M.EQ.2) THEN
        TOTAL = (2.0d0*FAR-1.0d0)*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**2
        TOTAL = TOTAL + (1.0d0-FAR)*DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)
        DMVDRM = 2.0d0*DE*FAR*TOTAL
      ELSE IF (M.EQ.3) THEN
        TOTAL = (1.0d0-4.0d0*FAR)*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**3
        TOTAL = TOTAL+(6.0d0*FAR-3.0d0)
     1          *DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)
     2          *DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)
        TOTAL = TOTAL+(1.0d0-FAR)*DMADRM(2,RE,R,NBETA,BETA,NCN,ESEL)
        DMVDRM = 2.0d0*DE*FAR*TOTAL
      ELSE IF (M.EQ.4) THEN
        TOTAL = (8.0d0*FAR-1.0d0)*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**4
        TOTAL = TOTAL+(6.0d0-24.0d0*FAR)
     1          *DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**2
     2          *DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)
        TOTAL = TOTAL+(2.0d0*FAR-1.0d0)
     1          *(3.0d0*DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)**2
     2          +4.0d0*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)
     3          *DMADRM(2,RE,R,NBETA,BETA,NCN,ESEL))
        TOTAL = TOTAL+(1.0d0-FAR)*DMADRM(3,RE,R,NBETA,BETA,NCN,ESEL)
        DMVDRM = 2.0d0*DE*FAR*TOTAL
      ELSE IF (M.EQ.5) THEN
        TOTAL = (1.0d0-16.0d0*FAR)*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**5
        TOTAL = TOTAL+(80.0d0*FAR-10.0d0)
     1          *DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**3
     2          *DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)
        TOTAL = TOTAL+(1.0d0-4.0d0*FAR)
     1          *(15.0d0*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)
     2          *DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)**2
     3          +10.0d0*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)**2
     4          *DMADRM(2,RE,R,NBETA,BETA,NCN,ESEL))
        TOTAL = TOTAL+(2.0d0*FAR-1.0d0)
     1          *(10.0d0*DMADRM(1,RE,R,NBETA,BETA,NCN,ESEL)
     2          *DMADRM(2,RE,R,NBETA,BETA,NCN,ESEL)
     3          +5.0d0*DMADRM(0,RE,R,NBETA,BETA,NCN,ESEL)
     4          *DMADRM(3,RE,R,NBETA,BETA,NCN,ESEL))
        TOTAL = TOTAL+(1.0d0-FAR)*DMADRM(4,RE,R,NBETA,BETA,NCN,ESEL)
        DMVDRM = 2.0d0*DE*FAR*TOTAL
      ELSE
c
c** If a value of M is chosen that is beyond the range of this
c   function, write out a warning message and return a value of zero.
c
        WRITE(6,600) M
        DMVDRM = 0.0d0
      END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** DMVDRM ERROR ***',/4x,'Unable to calculate the',I3,'
     1th derivative of the potential V(R).')
      END
c***********************************************************************
      REAL*8 FUNCTION DMADRM(M,RE,R,NBETA,BETA,NCN,ESEL)
c***********************************************************************
c** This function will calculate the Mth partial derivative of ALPHA.
c
c** On entry:
c    M       is the order of the derivative.
c    RE      is the eqilibrium bond length.
c    R       is the intermolecular distance.
c    NBETA   is the number of beta parameters in the BETA expansion.
c    BETA(n) is the array of beta parameters.
c    NCN     is the power of the radial term.
c    ESEL    is a flag denoting which ALPHA expansion is to be used:
c            = 1 : use the radial expansion - (R-RE)
c            = 2 : use the Ogilvie-Tipping expansion - Z=(R-RE)/(R+RE)
c
c** On exit:
c    DMADRM  is the Mth derivative of ALPHA with respect to R.
c
      INTEGER M, NBETA, NCN, ESEL, I

      REAL*8 RE, R, BETA(0:NBETA), DMZDRM, DMBDRM, MCL, TOTAL
c
c** Calculating the Mth derivative of ALPHA for the (R-RE) expansion.
c
      IF (ESEL.EQ.1) THEN
        DMADRM = DBLE(M+1)*DMBDRM(M,RE,R,NBETA,BETA)
     1           + (R-RE)*DMBDRM(M+1,RE,R,NBETA,BETA)
c
c** Calculating the Mth derivative of ALPHA for the Z expansion.
c
      ELSE IF (ESEL.EQ.2) THEN
        TOTAL = DBLE(NCN)/R
        DO I = 1,M
          TOTAL = TOTAL*(-DBLE(I)/R)
        END DO
        DO I = 0,M+1
          TOTAL = TOTAL + MCL(M+1,I)*DMZDRM(M+1-I,RE,R)
     1                      *DMBDRM(I,RE,R,NBETA,BETA) 
        END DO
        DMADRM = TOTAL
      ELSE
c
c** If a value of ESEL is chosen that is beyond the range of this
c   function, write out a warning message and return a value of zero.
c
        WRITE(6,600) ESEL
        DMADRM = 0.0d0
      END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** DMADRM ERROR ***',/4x,'The value of ESEL=(',I3,') is
     1beyond the range of this function.')
      END
c***********************************************************************
      REAL*8 FUNCTION MCL(M,L)
c***********************************************************************
c** This function will calculate the M choose L.
c
c** On entry:
c    M   is the top value.
c    L   is the bottom value.
c
c** On exit:
c    MCL is the value of M choose L.
c
      INTEGER M, L, I

      REAL*8 TOTAL

      IF ((M.GE.0).AND.(L.GE.0).AND.(M.GE.L)) THEN
        TOTAL = 1.0d0
        DO I = M-L+1,M
          TOTAL = TOTAL*DBLE(I)
        END DO
        DO I = 1,L
          TOTAL = TOTAL/DBLE(I)
        END DO
        MCL = TOTAL
      ELSE
c
c** If the chosen values of M and L are beyond the range of this
c   function, write out a warning message and return a value of zero.
c
        WRITE(6,600) M, L
        MCL = 0.0d0
      END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** MCL ERROR ***',/4x,'The value of M=(',I3,') and/or L
     1=(',I3,') are beyond the range of this function.')
      END
c***********************************************************************
      REAL*8 FUNCTION DMBDRM(M,RE,R,NBETA,BETA)
c***********************************************************************
c** This function will calculate the Mth partial derivative of BETA(Z)
c
c                     DMBDRM = (d^M/dR^M)BETA(Z)
c
c   where BETA(Z) is a polynomial expansion in terms of Z=(R-RE)/(R+RE):
c
c           BETA(Z) = BETA(0) + BETA(1)*Z + BETA(2)*Z^2 + ...
c
c** On entry:
c    M       is the order of the derivative.
c    RE      is the eqilibrium bond length.
c    R       is the intermolecular distance.
c    NBETA   is the number of beta parameters in the expansion.
c    BETA(n) is the array of beta parameters.
c
c** On exit:
c    DMBDRM  is the Mth derivative of BETA(Z) with respect to R.
c
      INTEGER M, NBETA, I

      REAL*8 RE, R, BETA(0:NBETA), DMZDRM, ZTEMP, TOTAL(5)

      ZTEMP = (R-RE)/(R+RE)
      DO I = 1,5
        TOTAL(I) = 0.0d0
      END DO
c
c** Now calculate the temporary sums.
c
      DO I = NBETA,1,-1
        TOTAL(1) = TOTAL(1)*ZTEMP+BETA(I)*DBLE(I)
      END DO
      DO I = NBETA,2,-1
        TOTAL(2) = TOTAL(2)*ZTEMP+BETA(I)*DBLE(I*(I-1))
      END DO
      DO I = NBETA,3,-1
        TOTAL(3) = TOTAL(3)*ZTEMP+BETA(I)*DBLE(I*(I-1)*(I-2))
      END DO
      DO I = NBETA,4,-1
        TOTAL(4) = TOTAL(4)*ZTEMP+BETA(I)*DBLE(I*(I-1)*(I-2)*(I-3))
      END DO
      DO I = NBETA,5,-1
        TOTAL(5) =TOTAL(5)*ZTEMP+BETA(I)*DBLE(I*(I-1)*(I-2)*(I-3)*(I-4))
      END DO
c
c** Now calculate the actual contributions from each sum.
c
      IF (M.EQ.0) THEN
        DO I = NBETA,0,-1
          TOTAL(1) = TOTAL(1)*ZTEMP+BETA(I)
        END DO
        DMBDRM = TOTAL(1)
      ELSE IF (M.EQ.1) THEN
        TOTAL(1) = TOTAL(1) * DMZDRM(1,RE,R)
        DMBDRM = TOTAL(1)
      ELSE IF (M.EQ.2) THEN
        TOTAL(1) = TOTAL(1) * DMZDRM(2,RE,R)
        TOTAL(2) = TOTAL(2) * DMZDRM(1,RE,R)**2
        DMBDRM = TOTAL(1) + TOTAL(2)
      ELSE IF (M.EQ.3) THEN
        TOTAL(1) = TOTAL(1) * DMZDRM(3,RE,R)
        TOTAL(2) = TOTAL(2) * 3.0d0*DMZDRM(2,RE,R)*DMZDRM(1,RE,R)
        TOTAL(3) = TOTAL(3) * DMZDRM(1,RE,R)**3
        DMBDRM = TOTAL(1) + TOTAL(2) + TOTAL(3)
      ELSE IF (M.EQ.4) THEN
        TOTAL(1) = TOTAL(1) * DMZDRM(4,RE,R)
        TOTAL(2) = TOTAL(2) * (4.0d0*DMZDRM(3,RE,R)*DMZDRM(1,RE,R)
     1                        + 3.0d0*DMZDRM(2,RE,R)**2)
        TOTAL(3) = TOTAL(3) * 6.0d0*DMZDRM(2,RE,R)*DMZDRM(1,RE,R)**2
        TOTAL(4) = TOTAL(4) * DMZDRM(1,RE,R)**4
        DMBDRM = TOTAL(1) + TOTAL(2) + TOTAL(3) + TOTAL(4)
      ELSE IF (M.EQ.5) THEN
        TOTAL(1) = TOTAL(1) * DMZDRM(5,RE,R)
        TOTAL(2) = TOTAL(2) * (5.0d0*DMZDRM(4,RE,R)*DMZDRM(1,RE,R)
     1                        + 10.0d0*DMZDRM(3,RE,R)*DMZDRM(2,RE,R))
        TOTAL(3) = TOTAL(3) * (10.0d0*DMZDRM(3,RE,R)*DMZDRM(1,RE,R)**2
     1                        + 15.0d0*DMZDRM(2,RE,R)**2*DMZDRM(1,RE,R))
        TOTAL(4) = TOTAL(4) * 10.0d0*DMZDRM(2,RE,R)*DMZDRM(1,RE,R)**3
        TOTAL(5) = TOTAL(5) * DMZDRM(1,RE,R)**5
        DMBDRM = TOTAL(1) + TOTAL(2) + TOTAL(3) + TOTAL(4) + TOTAL(5)
      ELSE
c
c** If a value of M is chosen that is beyond the range of this
c   function, write out a warning message and return a value of zero.
c
        WRITE(6,600) M
        DMBDRM = 0.0d0
      END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** DMBDRM ERROR ***',/4x,'Unable to calculate the',I3,'
     1th derivative of the BETA(Z) function.')
      END 
c***********************************************************************
      REAL*8 FUNCTION DMZDRM(M,RE,R)
c***********************************************************************
c** This function will calculate the Mth partial derivative d^MZ/dR^M.
c
c** On entry:
c    M      is the order of the derivative.
c    RE     is the eqilibrium bond length.
c    R      is the intermolecular distance.
c
c** On exit:
c    DMZDRM is the Mth derivative of Z=(R-RE)/(R+RE) with respect to R.
c
      INTEGER M, I

      REAL*8 RE, R, TOTAL

      TOTAL = -2.0d0*RE/(R+RE)
      DO I = 1,M
        TOTAL = -TOTAL*DBLE(M)/(R+RE)
      END DO
      IF (M.EQ.0) THEN
        DMZDRM = (R-RE)/(R+RE)
      ELSE IF (M.GT.0) THEN
        DMZDRM = TOTAL
      ELSE
c
c** If a value of M is chosen that is beyond the range of this
c   function, write out a warning message and return a value of zero.
c
        WRITE(6,600) M
        DMZDRM = 0.0d0
      END IF
      RETURN
c-----------------------------------------------------------------------
  600 FORMAT('  *** DMZDRM ERROR ***',/4x,'Unable to calculate the',I3,'
     1th derivative of Z.')
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
