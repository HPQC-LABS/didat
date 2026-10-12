c***********************************************************************
      SUBROUTINE QDVIRDP(IDAT,ISTATE,ZMU,BVIR,dBVIRdP)
c=======================================================================
c  This subroutine calculates the second virial coefficient BVIR and its
c  partial derivatives with respect to the various parameters. It is 
c  used when virial data has been input into the program.
c=======================================================================
      INCLUDE 'arrsizes.h'
      INCLUDE 'BLKPOT.h'
      INCLUDE 'BLKBOBRF.h'
      INCLUDE 'BLKCOUNT.h'
      INCLUDE 'BLKDVDP.h'
      INCLUDE 'BLKISOT.h'
      INCLUDE 'BLKDATA.h'
c* Define types for local variables
       INTEGER check(0:NPARMX+3),j,counter,nParams,ISTATE,ISTART,
     1 IDAT,i,m,k,kk,jj,ii,f,BANDN,L,IISTP,NODX,N,IWR,LIFE,NARRAY
       REAL*8  XG(8),WG(8),x(4),w(4),BTemp,BTempInv,BTOT(3),
     1 dBTOTdp(NPARMX),CONST(3),jump,a,b,YVAL(8),RVAL(8),
     2 VDIST(8),EXP_TERM,int_fact,Vsq,VPsq,Rsq,RINV,XTEMP,Bclass,Bq1,
     3 Bq2,INTEGRALS(NPARMX+3),error,BVIR,BETADISTa,dVdR(8),d2VdR(8),
     4 XTEMP1,dBcdP(NPARMX),dBq1dP(NPARMX),dBVIRdP(NPARMX),ZMU,etime,
     5 elapsed(2),timestart,timestop,muFCT,BFCT,V(NPNTMX),RM2(NPNTMX),
     6 SUML,ZE,ADLT,STD,BZ,PHTST,TDTST
c
       REAL*8, PARAMETER :: k_boltz=1.3806504D-23, NA=6.02214179D23,
     1  h = 6.62606896D-34, c = 2.99792458D10, k_cm = k_boltz/(h*c),
     2  h_cm = 1.d0/c, pi = 3.14159265358979323846, amu=1.660538782D-27,
     3  Cu= 16.8576292D0

      IWR= 1
      PHTST= 1.d-5
      LIFE= 0
      BZ= ZMU/Cu
      NARRAY= NDATPT(ISTATE)

c    integral over energy  enclosing
c    sum over  l

c.. takes the current temperature from the virial data and changes it to E cm
      BTemp = k_cm*Temp(IDAT)
c.. inverts BTemp
      BTempInv = 1.d0/BTemp

      ZE= BTEMP

c.. initializes the total correction term and its derivative to zero
      BTOT = 0.d0
      dBTOTdP = 0.d0


c?????
c
c** Generating potential function including centrifugal, adiabatic BOB,
c   and rotational non-adiabatic BOB, and if appropriate, Lambda 
c   doubling or 2\Sigma doubling radial functions.
cc
c==> Model-A <==
c** First - for Li2(A), add BOB centrifugal shift  2*B(r)
cc    IF((NCMM(ISTATE).GE.3).AND.(MMLR(2,ISTATE).LE.0)
cc   1                      .AND.(IOMEG(ISTATE).NE.-2))
cc   2              JFCT= JFCT + 2.d0*ZMASS(3,1)/ZMASS(3,IISTP) - 2.d0
c==> Model-A <==
cc



c** To calculate the values for <vj|dV/dP(k)|vj>
c   we must first determine the wave equation for the vj state.
      muFCT= 16.85762920d0/ZMU
      BFCT= RH(ISTATE)*RH(ISTATE)/muFCT
      bandN= IB(IDAT)

      DO I= 1,NDATPT(ISTATE)
          V(I)= BFCT*VPOT(I,ISTATE)
cc        V(I)= BFCT*(VPOT(I,ISTATE) + ZMUA(IISTP,ISTATE)*UAR(I,ISTATE)
cc   1                              + ZMUB(IISTP,ISTATE)*UBR(I,ISTATE))
          RM2(I)= 1.0d0/RD(I,ISTATE)**2
cc        RM2(I)= (1.0d0 + ZMTA(IISTP,ISTATE)*TAR(I,ISTATE)
cc   1             + ZMTB(IISTP,ISTATE)*TBR(I,ISTATE))/RD(I,ISTATE)**2
          ENDDO



      SUML=0.0D0
      DO  L= 0,1000
          CALL PHSHTD(L,ZE,ADLT,STD,NODX,NARRAY,V,RM2,BZ,RH,RMIN,PHTST,
     1               TDTST,IWR,LIFE)
          SUML=SUML+(2*L+1)*ADLT
          WRITE(6,600) SUML 
          ENDDO
  600 FORMAT(1X,'Virial Sum',D10.4)


      BVIR= Const(1)*BTOT(1) + Const(2)*BTOT(2) + Const(3)*BTOT(3)
      DO J= 1,nParams
          dBVIRdP(J)= dBTOTdp(J)
          ENDDO
      RETURN
      END 

c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
