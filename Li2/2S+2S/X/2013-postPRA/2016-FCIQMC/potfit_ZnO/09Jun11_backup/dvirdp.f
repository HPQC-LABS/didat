c***********************************************************************
      SUBROUTINE DVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
       INTEGER check(0:NPARMX+3),j,counter,nParams,ISTATE,ISTART,
     1 IDAT,i,m,k,kk,jj,ii,f
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,BTOT(3),
     1 dBTOTdp(NPARMX),CONST(3),jump,a,b,YVAL(8),RVAL(8),
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,BETADISTa,dVdR(8),d2VdR(8),
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,etime,
     5 elapsed(2),timestart,timestop
c
       REAL*8, PARAMETER :: k_boltz=1.3806504D-23, NA=6.02214179D23,
     1  h = 6.62606896D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538782D-27,
     3  Cu= 16.8576292D0
c* Common block data for quadrature weights and points
      data f/0/
      data x/0.960289856497536d0, 0.796666477413627d0,
     1       0.525532409916329d0, 0.183434642495650d0/,
     2     w/0.101228536290376d0, 0.222381034453374d0,
     3       0.313706645877887d0, 0.362683783378362d0/
      save f
c***********************************************************************
c.. runs a loop to set all quadrature weights and points
      DO j= 1,4
          XG(j)= -x(j)
          WG(j)= w(j)
          XG(9-j)= x(j)
          WG(9-j)= w(j)
          ENDDO

c.. initializes the array dBVIRdP to zero.
      DO J= 1,HPARMX
          dBVIRdP(J)= 0.d0
          ENDDO

c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)

c.. inverts BTemp
      BTempInv = 1.d0/BTemp

c.. initializes the total correction term and its derivative to zero
      BTOT = 0.d0
      dBTOTdP = 0.d0

c.. sets the constants for the intergration of each correction term
      Const(1)= -2.d0*pi*0.602214179d0
      Const(2)= Cu*pi*BTempInv**3/(6.d0*ZMU)
      Const(3)= -(pi/6.d0)*(Cu*BTempInv**2/ZMU)**2

c.. sets the number of paramters used for the integration and initializes
c.  the integrals to zero
      nParams = 3 + NCMM(ISTATE) + NLpow(ISTATE)
      DO ii=0,nParams+3
          check(ii) = 1
          ENDDO
      DO J= 1,nParams+3
          INTEGRALS(J)= 0.d0
          ENDDO

c.. this first loop is here to repeat the calculations until the values
c.  converge to within error, or until they are repeated 7 times without
c.  convergence
      DO  kk= 0,7
          IF(check(0).EQ.-1) EXIT
          counter= 2.d0**kk

c.. futher subdivides the interval for gaussian integration formula
          jump= 2.d0/DBLE(counter)
          b= -1.d0
c. the first loop is over each subinterval within (-1,1)
          DO m= 1,counter
              a= b
              b= a+ jump
              IF(m.EQ.counter) b= 1.d0
              DO i= 1,8
c.. sets the y value for the gaussian formula
                  YVAL(i)= 0.5d0*((b - a)*XG(i) + (b + a))
c... maping  y<->r for a particular choice of dimensionless variable ...
c   where 'real' range is  0.01*Re  to infty - with  'p=2'
                  RVAL(i)= Re(ISTATE)*DSQRT((1.0001d0 
     1                           + 0.9999d0*YVAL(i))/(1.d0 - YVAL(i)))
                  VDIST(i)= 0.d0
                  ENDDO

c.. calls VGENp to find the neccesary derivatives
              CALL VGENp(ISTATE,RVAL,VDIST,BETADISTa,IDAT,dVdR,d2VdR)

c.. the next loop is over each Gaussian point within each subinterval
              DO i= 1,8
c.. some of the terms that will be used in later calculations of the
c.  virial coefficients are constructed here
                  EXP_TERM= DEXP(-VDIST(i)*BTempINV)
                  RINV= 1.d0/RVAL(i)
                  int_fact= (Re(ISTATE)/(1.d0 - YVAL(i)))**2
     1                                          *RINV*(b - a)*0.5d0
                  Vsq= VDIST(i)**2
                  VPsq= dVdR(i)**2
                  Rsq= RVAL(i)**2
                  XTEMP= d2VdR(i)**2/10.d0 + VPsq*RINV**2/5.d0
     1      + dVdR(i)**3*BTempInv*RINV/9.d0 - (VPsq*BTempInv)**2/72.d0
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c*     Finally the relevant partial derivatives of the virial 
c*     coefficients are calculated
                  DO J = 1,nParams
                      jj = ISTART + J
                      IF(check(J).eq.1) THEN
c* the derivative of the classical expression
                          dBcdP(J)= -BTempInv*EXP_TERM
     1                                       *DVtot(jj,i)*Rsq*int_fact
c* and of the first quantum correction 
                          XTEMP1 = -VPsq*BTempInv*DVtot(jj,i)
     1                                         + 2*dVdR(i)*dVpdP(jj,i)
                          dBq1dP(J)= EXP_TERM*XTEMP1*Rsq*int_fact
  
c..  As the final step these terms all added together in a weighted sum
  
                          INTEGRALS(J) = INTEGRALS(J) + Const(1)
     1                      *dBcdP(J)*WG(I) + Const(2)*dBq1dP(J)*WG(I)
                          ENDIF
                      ENDDO
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c. now the integrands are evaluated at each particular Gaussian point
c. and summed together with the proper weighting
                  IF(check(nParams+1).EQ.1) THEN
                      Bclass= (EXP_TERM - 1.d0)*Rsq*int_fact
                      INTEGRALS(nParams+1)= INTEGRALS(NParams+1) 
     1                                                 + Bclass*WG(i)
                      ENDIF
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
                  ENDDO
              ENDDO
c. this section checks each integral to see if it has converged,
c. assigning a check value of -1 to those that have, which will stop
c. further individual calculations for them
          IF(counter.GE.2) THEN
              error= 0.000001
              DO i= 1,nParams
                  IF(integrals(i).NE.0.d0) THEN
                      print*,i,temp(idat), integrals(i), check(i)
                      ENDIF
                  IF(check(i).EQ.1) THEN
                      IF(abs(INTEGRALS(i) - dBTOTdp(i)).LE.(abs(error
     1                                  *INTEGRALS(i)))) THEN
                         check(i)= -1
                         ENDIF
                      dBTOTdp(i)= INTEGRALS(i)
                      INTEGRALS(i)= 0.d0
                      ENDIF
                  ENDDO
              DO i= 1,3
                  IF(check(nParams+i).EQ.1) THEN
                     IF(abs(INTEGRALS(nParams+i)-BTOT(i)).LE.(abs(error
     1                                   *INTEGRALS(nParams+i)))) THEN
                         check(nParams+i)= -1
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
          ENDDO
      BVIR= Const(1)*BTOT(1) + Const(2)*BTOT(2) + Const(3)*BTOT(3)
      DO J= 1,nParams
          dBVIRdP(J)= dBTOTdp(J)
          ENDDO
      END 
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
