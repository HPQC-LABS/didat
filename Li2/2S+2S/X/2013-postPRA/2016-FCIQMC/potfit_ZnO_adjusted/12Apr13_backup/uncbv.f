c***********************************************************************
      SUBROUTINE UNCBV(NPTOT,PV,PU,CM)
c***********************************************************************
c**  Subroutine to compute the uncertainties in calcuated Bv values, 
c  1) by the finite difference R(0) approximation
c  2) from the LIDE expression
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** On entry:
c ... OSEL     print values at every OSEL'th mesh point
c ... ISTATE   electronic state counter
c ... PU(n)    uncertainties in the TOTPOTPAR parameters of the model
c ... CM(n,n)  (symmetric) correlation matrix from the fit
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKBOB.h'
      INCLUDE 'BLKPARAM.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDATA.h'
c-----------------------------------------------------------------------
      INTEGER IISTP, IV, I, NPTOT,NB1,COUNT1,ISTATE
      REAL*8 PV(NPARMX), PU(NPARMX), PD(NPARMX), CM(NPARMX,NPARMX), 
     1  PT(NPARMX), Bunc, CALC
c** First - set up fake R(0) datum
      COUNT1= COUNTOT+1
      NB1= NBANDTOT+1
      IB(COUNT1)= NB1
      JPP(COUNT1)= 0
      JP(COUNT1)= 1
      FREQ(COUNT1)= 0.d0
      UFREQ(COUNT1)= 0.d0
      DO ISTATE= 1, NSTATES
          IEP(NB1)= ISTATE
          IEPP(NB1)= ISTATE
          DO  IISTP= 1,NISTP
              ISTP(NB1)= IISTP
              WRITE(7,600) IISTP, SLABL(ISTATE)
c** for each vibrational level of each isotopologue ....
              DO  IV= VMIN(ISTATE,IISTP), VMAX(ISTATE,IISTP)
                  VP(NB1)= IV
                  VPP(NB1)= IV
c** First get partial derivatives from Bv as half of an R(0) energy
                  CALL DYIDPJ(COUNT1,COUNT1,NPTOT,FREQ(COUNT1),CALC,PV,
     1                                                             PD)
c ... set up column vector for uncertainty calculation
                  DO  I= 1,NPTOT
                      PT(I)= 0.5*PD(I)*PU(I)
                      ENDDO
                  CALL MMCALC(1,NPTOT,PT,CM,Bunc)
c   
                  WRITE(7,602)  IV,CALC/2,Bunc
                  ENDDO
              ENDDO
          ENDDO
      RETURN
  600 FORMAT(/'   Predict Bv uncertainties for isotopologue-',I1,' in st
     1ate ',A2)
  602 FORMAT(' For  v=',I3,'  Bv approx R(0)/2 =', F11.8,'  unc(Bv)=',
     1                                                        1PD12.6)
      END
c234567890 234567890 234567890 234567890 234567890 234567890 234567890
