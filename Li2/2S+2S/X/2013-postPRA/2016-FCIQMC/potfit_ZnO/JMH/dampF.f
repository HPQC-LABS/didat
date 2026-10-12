************************************************************************
      SUBROUTINE dampF(r,beta,n,Dn,Dnp,Dnpp,form)
*** Subroutine to generate values 'Dn' of the chosen potential damping
*   function used
*   form = 1 :  Tang-Toennies damping function (J.C.P.80,3726(1984))
*   form = 2 :  damped-dispersion damping function (M.P.,52,763(1984))
*   form = 3 :  damped-dispersion damping function (M.P.,52,763(1984))
*               (with correction on unit)
*   Modified on Wed Jul 18 14:28:49 EDT 2001
************************************************************************

      INTEGER i,n,form
      REAL*8 Dn,beta,r,Dnp,Dnpp,power,Coef,expPower
      REAL*8 betaR,eBetaR,sum(n),tempSum

      IF (form .EQ. 1) THEN
        betaR = r*beta
        eBetaR = DEXP(-betaR)
        DO i = 1,n
          sum(i) = 1.0d0
        END DO
        tempSum = 1.0d0
        
        DO i = 1,n
          tempSum = tempSum*betaR/i
          sum(i) = sum(i)+tempSum
        ENDDO
        Dn = 1.0d0-eBetaR*sum(n)
        Dnp = beta*eBetaR*(sum(n) - sum(n-1))
        Dnpp = beta**2*eBetaR*(-sum(n) + 2*sum(n-1) - sum(n-2))
      ELSE IF (form .EQ. 2) THEN
        power = -2.1d0*beta*r/n-0.109d0*(beta*r)**2/SQRT(n*1.0d0)
        expPower = exp(power)
        Dn = (1.0d0 - expPower)**n
        Coef = 2.1d0*beta/n + 0.218d0*beta**2*r/SQRT(n*1.0d0)
        Dnp = n*Coef*expPower*(1.0d0 - expPower)**(n-1)
        Dnpp = expPower*(1.0d0 - expPower)**(n-1)*
     &    (0.218d0*SQRT(n*1.0d0)*beta**2-n*Coef**2
     &    *(1.0d0-n*expPower)/(1.0d0 - expPower))
      ELSE IF (form .EQ. 3 .OR. form .EQ. 4) THEN
        power = -2.1d0/0.5291772083d0*beta*r/n-0.109d0/0.5291772083d0**2
     &          *(beta*r)**2/SQRT(n*1.0d0)
        expPower = exp(power)
        Dn = (1.0d0 - expPower)**n
        Coef = 2.1d0/0.5291772083d0*beta/n + 0.218d0/0.5291772083d0**2
     &         *beta**2*r/SQRT(n*1.0d0)
        Dnp = n*Coef*expPower*(1.0d0 - expPower)**(n-1)
        Dnpp = expPower*(1.0d0 - expPower)**(n-1)*
     &    (0.218d0/0.5291772083d0**2*SQRT(n*1.0d0)*beta**2-n*Coef**2
     &    *(1.0d0-n*expPower)/(1.0d0 - expPower))
      END IF

      RETURN
      END
