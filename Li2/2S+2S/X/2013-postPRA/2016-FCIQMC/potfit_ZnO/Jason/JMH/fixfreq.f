***********************************************************************
      SUBROUTINE DIFFSTATS(NSTATES,PROMC)
c** Subroutine to summarise dimensionless standard errors on a band-by-
c  band basis, and (if desired) print [obs.-calc.] values to channel-8.
c-----------------------------------------------------------------------
c                 Version of 6 April 2001
c-----------------------------------------------------------------------
      INCLUDE 'arrsizes.h'
c
      REAL*8 SSQTOT,AVE,AVETOT,RMSR,RMSTOT
      INTEGER I,IBB,ISOT,ISTATE,ISTATEE,J,NSTATES,PROMC
c
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)
      INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)
      CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX)
      CHARACTER*3 NEF(-1:1)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
c
c** Type statements & common blocks for characterizing transitions
c
      REAL*8  AVEUFREQ(NBANDMX),MAXUFREQ(NBANDMX)
      INTEGER NTRANSFS(NISTPMX,NSTATEMX),
     1  NTRANSVIS(NISTPMX,NSTATEMX,NSTATEMX),
     1  NBANDEL(NISTPMX,NSTATEMX,NSTATEMX),
     2  NTRANSIR(NISTPMX,NSTATEMX),NTRANSMW(NISTPMX,NSTATEMX),
     3  NBANDFS(NISTPMX,NSTATEMX),NBANDVIS(NISTPMX,NSTATEMX),
     4  NBANDIR(NISTPMX,NSTATEMX),NBANDMW(NISTPMX,NSTATEMX),
     5  NBVPP(NISTPMX,NSTATEMX),NWIDTH(NISTPMX,NSTATEMX),
     6  NVD(NISTPMX,NSTATEMX),
     6  NBANDS(NISTPMX),YPR(NISTPMX,NSTATEMX,7,6,NBANDMX)
c
      COMMON /TYPEBLK/AVEUFREQ,MAXUFREQ,NTRANSFS,NTRANSVIS,NTRANSIR,
     1  NTRANSMW,NBANDFS,NBANDEL,NBANDVIS,NBANDIR,NBANDMW,NBVPP,NWIDTH,
     2  NVD,NBANDS,YPR
      DATA  NEF/'  f','   ','  e'/
c========================================================================
      ISOT= 1
      SSQTOT= 0.d0
c** Summarize data discrepancies for one isotopomer at a time.
   10 WRITE(6,602) NBANDS(ISOT),(NAME(I),MN(I,ISOT),I= 1,2)
c
c** Loop over bands for each (lower) electronic state, in turm
      DO 50 ISTATE= 1,NSTATES
      IF(NTRANSMW(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Micowave data
          WRITE(6,604) NTRANSMW(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          WRITE(6,605)
          IF(ABS(PROMC).GT.1) WRITE(8,604) NTRANSMW(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDMW(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDMW(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,4,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,SSQTOT,
     1                                              DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,4,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,4,3,I)*AVE
              WRITE(6,606)YPR(ISOT,ISTATE,4,2,I),YPR(ISOT,ISTATE,4,1,I),
     1                  YPR(ISOT,ISTATE,4,3,I),YPR(ISOT,ISTATE,4,5,I),
     2                  YPR(ISOT,ISTATE,4,6,I),
     3                  AVEUFREQ(YPR(ISOT,ISTATE,4,4,I)),
     4                  MAXUFREQ(YPR(ISOT,ISTATE,4,4,I)),AVE,RMSR
              IF(ABS(PROMC).GT.1) THEN
                  WRITE(8,605)
                  WRITE(8,606) YPR(ISOT,ISTATE,4,2,I),
     1                         YPR(ISOT,ISTATE,4,1,I),
     2                  YPR(ISOT,ISTATE,4,3,I),YPR(ISOT,ISTATE,4,5,I),
     3                  YPR(ISOT,ISTATE,4,6,I),
     4                  AVEUFREQ(YPR(ISOT,ISTATE,4,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,4,4,I)),AVE,RMSR
                  CALL PBNDERR(IBB,NEF)
                  ENDIF
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSMW(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSMW(ISOT,ISTATE)
          WRITE(6,630) NTRANSMW(ISOT,ISTATE),AVETOT,RMSTOT
	    ENDIF
c
      IF(NTRANSIR(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Infrared data
          WRITE(6,608) NTRANSIR(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                         MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          WRITE(6,605)
          IF(ABS(PROMC).GT.1) WRITE(8,608) NTRANSIR(ISOT,ISTATE),
     1  SLABL(ISTATE),(NAME(I),MN(I,ISOT),I= 1,2),NBANDIR(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDIR(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,3,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,SSQTOT,
     1                                             DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,3,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,3,3,I)*AVE
              WRITE(6,606)YPR(ISOT,ISTATE,3,2,I),YPR(ISOT,ISTATE,3,1,I),
     1                  YPR(ISOT,ISTATE,3,3,I),YPR(ISOT,ISTATE,3,5,I),
     2                  YPR(ISOT,ISTATE,3,6,I),
     3                  AVEUFREQ(YPR(ISOT,ISTATE,3,4,I)),
     4                  MAXUFREQ(YPR(ISOT,ISTATE,3,4,I)),AVE,RMSR
              IF(ABS(PROMC).GT.1) THEN
                  WRITE(8,605)
                  WRITE(8,606) YPR(ISOT,ISTATE,3,2,I),
     1                         YPR(ISOT,ISTATE,3,1,I),
     2                  YPR(ISOT,ISTATE,3,3,I),YPR(ISOT,ISTATE,3,5,I),
     3                  YPR(ISOT,ISTATE,3,6,I),
     4                  AVEUFREQ(YPR(ISOT,ISTATE,3,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,3,4,I)),AVE,RMSR
                  CALL PBNDERR(IBB,NEF)
                  ENDIF
              ENDDO
          RMSTOT= DSQRT(RMSTOT/NTRANSIR(ISOT,ISTATE))
          AVETOT= AVETOT/NTRANSIR(ISOT,ISTATE)
          WRITE(6,630) NTRANSIR(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
c
c** Book-keeping for Electronic vibrational band data
      DO  ISTATEE= 1,NSTATES
          IF((ISTATEE.NE.ISTATE).AND.
     1                 (NTRANSVIS(ISOT,ISTATEE,ISTATE).GT.0)) THEN
c ... for ISTATEE{upper}-ISTATE{lower} electronic vibrational bands
              WRITE(6,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     1                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              WRITE(6,605)
              IF(ABS(PROMC).GT.1) 
     1                    WRITE(8,610) NTRANSVIS(ISOT,ISTATEE,ISTATE),
     2                      (NAME(I),MN(I,ISOT),I=1,2),SLABL(ISTATEE),
     2                      SLABL(ISTATE),NBANDEL(ISOT,ISTATEE,ISTATE)
              RMSTOT= 0.d0
              AVETOT= 0.d0
              DO  I= 1,NBANDVIS(ISOT,ISTATE)
                  IBB= YPR(ISOT,ISTATE,2,4,I)
                  IF(IEP(IBB).EQ.ISTATEE) THEN

                      CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,
     1                                             SSQTOT,DFREQ,UFREQ)
                      RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,2,3,I)*RMSR**2
                      AVETOT= AVETOT+ YPR(ISOT,ISTATE,2,3,I)*AVE
                      WRITE(6,606) YPR(ISOT,ISTATE,2,2,I),
     1                             YPR(ISOT,ISTATE,2,1,I),
     2                  YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3                  YPR(ISOT,ISTATE,2,6,I),
     4                  AVEUFREQ(YPR(ISOT,ISTATE,2,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,2,4,I)),AVE,RMSR
                      IF(ABS(PROMC).GT.1) THEN
                          WRITE(8,605)
                          WRITE(8,606) YPR(ISOT,ISTATE,2,2,I),
     1                                 YPR(ISOT,ISTATE,2,1,I),
     2            YPR(ISOT,ISTATE,2,3,I),YPR(ISOT,ISTATE,2,5,I),
     3        YPR(ISOT,ISTATE,2,6,I),AVEUFREQ(YPR(ISOT,ISTATE,2,4,I)),
     5        MAXUFREQ(YPR(ISOT,ISTATE,2,4,I)),AVE,RMSR
                          CALL PBNDERR(IBB,NEF)
                          ENDIF
                      ENDIF
                  ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSVIS(ISOT,ISTATEE,ISTATE))
              AVETOT= AVETOT/NTRANSVIS(ISOT,ISTATEE,ISTATE)
              WRITE(6,630) NTRANSVIS(ISOT,ISTATEE,ISTATE),AVETOT,RMSTOT
              ENDIF
          ENDDO
      IF(NTRANSFS(ISOT,ISTATE).GT.0)THEN
c** Book-keeping for Fluorescence data
          WRITE(6,612) NTRANSFS(ISOT,ISTATE),SLABL(ISTATE),(NAME(I),
     1                          MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          WRITE(6,617)
          IF(ABS(PROMC).GT.1) WRITE(8,612) NTRANSFS(ISOT,ISTATE),
     1   SLABL(ISTATE),(NAME(I),MN(I,ISOT),I=1,2),NBANDFS(ISOT,ISTATE)
          RMSTOT= 0.d0
          AVETOT= 0.d0
          DO  I= 1,NBANDFS(ISOT,ISTATE)
              IBB= YPR(ISOT,ISTATE,1,4,I)
              CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,SSQTOT,
     1                                              DFREQ,UFREQ)
              RMSTOT= RMSTOT+ YPR(ISOT,ISTATE,1,3,I)*RMSR**2
              AVETOT= AVETOT+ YPR(ISOT,ISTATE,1,3,I)*AVE
              WRITE(6,614) YPR(ISOT,ISTATE,1,1,I),
     1                     YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(YPR(ISOT,ISTATE,1,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,1,4,I)),AVE,RMSR
              IF(ABS(PROMC).GT.1) THEN
                  WRITE(8,617)
                  WRITE(8,614) YPR(ISOT,ISTATE,1,1,I),
     1                  YPR(ISOT,ISTATE,1,2,I),NEF(EFP(IFIRST(IBB))),
     2                  YPR(ISOT,ISTATE,1,3,I),YPR(ISOT,ISTATE,1,5,I),
     3                  YPR(ISOT,ISTATE,1,6,I),
     4                  AVEUFREQ(YPR(ISOT,ISTATE,1,4,I)),
     5                  MAXUFREQ(YPR(ISOT,ISTATE,1,4,I)),AVE,RMSR
                  CALL PBNDERR(IBB,NEF)
                  ENDIF
              ENDDO
              RMSTOT= DSQRT(RMSTOT/NTRANSFS(ISOT,ISTATE))
              AVETOT= AVETOT/NTRANSFS(ISOT,ISTATE)
              WRITE(6,630) NTRANSFS(ISOT,ISTATE),AVETOT,RMSTOT
          ENDIF
      IF(NBVPP(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Bv  data
          IBB= YPR(ISOT,ISTATE,5,4,1)
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,SSQTOT,DFREQ,
     1                                                          UFREQ)

          WRITE(6,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(YPR(ISOT,ISTATE,5,4,1)),
     4              MAXUFREQ(YPR(ISOT,ISTATE,5,4,1)),AVE,RMSR
          IF(ABS(PROMC).GT.0) THEN
              WRITE(8,616) NBVPP(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,5,3,1),
     2              YPR(ISOT,ISTATE,5,5,1),YPR(ISOT,ISTATE,5,6,1),
     3              AVEUFREQ(YPR(ISOT,ISTATE,5,4,1)),
     4              MAXUFREQ(YPR(ISOT,ISTATE,5,4,1)),AVE,RMSR
              ENDIF
          IBB= YPR(ISOT,ISTATE,5,4,1)
          DO  J= IFIRST(YPR(ISOT,ISTATE,5,4,1)),
     1                                   ILAST(YPR(ISOT,ISTATE,5,4,1))
              WRITE(6,618) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              IF(ABS(PROMC).GT.0)
     1            WRITE(8,618) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     2                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
c
      IF(NWIDTH(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for  Tunneling Width  data
          IBB= YPR(ISOT,ISTATE,6,4,1)    
          CALL BNDERR(IFIRST(IBB),ILAST(IBB),AVE,RMSR,SSQTOT,DFREQ,
     1                                                          UFREQ)
          WRITE(6,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(YPR(ISOT,ISTATE,6,4,1)),
     4              MAXUFREQ(YPR(ISOT,ISTATE,6,4,1)),AVE,RMSR
          IF(ABS(PROMC).GT.0) THEN
              WRITE(8,620) NWIDTH(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2),YPR(ISOT,ISTATE,6,3,1),
     2              YPR(ISOT,ISTATE,6,5,1),YPR(ISOT,ISTATE,6,6,1),
     3              AVEUFREQ(YPR(ISOT,ISTATE,6,4,1)),
     4              MAXUFREQ(YPR(ISOT,ISTATE,6,4,1)),AVE,RMSR
              ENDIF
          IBB= YPR(ISOT,ISTATE,6,4,1)    
          DO  J= IFIRST(YPR(ISOT,ISTATE,6,4,1)),
     1                                   ILAST(YPR(ISOT,ISTATE,6,4,1))
              WRITE(6,622) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              IF(ABS(PROMC).GT.0) 
     1            WRITE(8,622) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     2                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDDO
          ENDIF
      IF(NVD(ISOT,ISTATE).GT.0) THEN
c** Book-keeping for v_D data
          WRITE(6,621) NVD(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2)
          IF(ABS(PROMC).GT.0) THEN
              WRITE(8,621) NVD(ISOT,ISTATE),SLABL(ISTATE),
     1              (NAME(I),MN(I,ISOT),I=1,2)
              ENDIF
          J=IFIRST(YPR(ISOT,ISTATE,7,1,1))
          WRITE(6,622) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
          IF(ABS(PROMC).GT.0) THEN
              WRITE(8,622) JPP(J),NEF(EFPP(J)),FREQ(J),UFREQ(J),
     1                                      DFREQ(J),DFREQ(J)/UFREQ(J)
              ENDIF
          ENDIF
c** End of loop over the various (lower) electronic states
   50 CONTINUE
c=======================================================================
      IF(ISOT.LT.NISTP) THEN
c** If NISTP > 1, return to print data summaries for other isotopomers
          ISOT= ISOT+1
          GO TO 10
          ENDIF 
      RMSR= DSQRT(SSQTOT/COUNTOT)
      WRITE(6,624) COUNTOT,RMSR
      RETURN
  602 FORMAT(/1x,21('===')/'  *** Discrepancies for',I5,' bands/series o
     1f ',A2,'(',I3,')-',A2,'(',I3,') ***'/1x,21('==='))
  604 FORMAT(/1x,21('===')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') MW transitions in',i4,' vib. levels')
  605 FORMAT(1x,16('==='),'== Avge. ========'/"   v' ",
     2  ' v" #data  J"min  J"max  Av.Unc.  Max.Unc.   Err/Unc    RMSR'/
     1  1x,13('-----'))
  606 FORMAT(2I4,I6,3x,I4,3x,I4,1x,1P2D9.1,0PF11.5,F8.3)
  608 FORMAT(/1x,63('=')/I5,' State ',A2,1x,A2,'(',I3,')-',A2,'(',I3,
     1 ') InfraRed transitions in',I4,' bands')
  610 FORMAT(/1x,35('==')/I6,1x,A2,'(',I3,')-',A2,'(',i3,')  {State ',
     1  A2,'}--{State ',A2,'} Transitions in',i4,' bands')
  612 FORMAT(/1x,76('=')/I5,' Fluorescence transitions into State ',A2,
     1 2x,A2,'(',I3,')-',A2,'(',I3,')  in',i5,' series')
  617 FORMAT(1x,52('='),'= Avge. ',18('=')/"   v'  j' p' ",
     2  '#data  v"min  v"max','  AvgeUnc  Max.Unc.   Err/Unc    RMSR'/
     3  1x,40('--'))
  614 FORMAT(2I4,A3,I6,2I7,1x,1P2D10.1,0PF11.5,F8.3)
  616 FORMAT(/1x,66('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Bv values treated as independent data'/1x,20('=='),
     2 '  Avge.  ',17('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc   RMSR'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,30('==')/'    v  p',8x,'Bv',7x,'u(Bv)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  618 FORMAT(I5,A3,2x,F12.8,1PD9.1,0PF13.8,F12.4)
  620 FORMAT(/1x,73('=')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Tunneling Widths treated as independent data'/1x,20('=='),
     2 '  Avge.  ',24('=')/' #data   v"min  v"max  AvgeUnc  Max.Unc.  Er
     3r/Unc   RMSR'/1x,55('-')/I5,2I7,2x,1P2D9.1,0PF9.3,F8.3/
     4 1x,59('=')/'   v  p     Width',7x,'u(Width)  [calc-obs]  [calc-ob
     5s]/unc'/1x,59('-'))
  621 FORMAT(/1x,30('--')/1x,I3,' State ',A2,1x,A2,'(',I3,')-',A2,'(',
     1 I3,') Quantum number at dissociation'/1x,30('--')/
     4 '    J  p',8x,'vD',7x,'u(vD)',4x,
     5 '[calc-obs]  [calc-obs]/unc',/1x,30('--'))
  622 FORMAT(I4,A3,1PD14.6,D10.1,2D13.2)
  624 FORMAT(1x,29('==')/' For overall fit to',i6,' data,  RMS(residuals
     1)=',G11.4/1x,29('==')) 
  630 FORMAT(1x,7('--'),' For these',i6,' lines, overall:',F11.5,F8.3)
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE BNDERR(FIRST,LAST,AVEDD,RMSDD,SSQTOT,DFREQ,UFREQ)
c** Calculate the average (AVEDD) & the root mean square dimensionless 
c  deviation (RSMDD) for the band running from datum # FIRST to LAST.
      INCLUDE 'arrsizes.h'
c
      REAL*8 DFREQ(NDATAMX),UFREQ(NDATAMX),AVEDD,RMSDD,SSQTOT
      INTEGER FIRST,LAST,NDAT,I
c
      AVEDD= 0.d0
      RMSDD= 0.d0
      DO 10 I= FIRST,LAST
          AVEDD= AVEDD+ DFREQ(I)/UFREQ(I)
   10     RMSDD= RMSDD+ (DFREQ(I)/UFREQ(I))**2
      SSQTOT= SSQTOT+ RMSDD
      NDAT= LAST-FIRST+1
      AVEDD= AVEDD/NDAT
      RMSDD= DSQRT(RMSDD/NDAT)
      RETURN
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

c***********************************************************************
      SUBROUTINE PBNDERR(IBB,NEF)
c** Print to channel-8 a listing of the [obs.-calc.] values for the band
c  running from datum # FIRST to LAST.                           
      INCLUDE 'arrsizes.h'             
c
cc    INTEGER IBB,I,J
      INTEGER IBB,I
c
c** Type statements & common block for data
c
      REAL*8  FREQ(NDATAMX),UFREQ(NDATAMX),DFREQ(NDATAMX),
     1  ZMASS(3,NISTPMX),RSQMU(NISTPMX),RSQMUP(0:NDUNMX,NISTPMX),
     2  RMUP(0:8,NISTPMX)
      INTEGER  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN(2),MN(2,NISTPMX),
     1  IB(NDATAMX),JP(NDATAMX),JPP(NDATAMX),VP(NBANDMX),VPP(NBANDMX),
     2  EFP(NDATAMX),EFPP(NDATAMX),FSBAND(NBANDMX),NFS(NBANDMX),
     3  IEP(NBANDMX),IEPP(NBANDMX),ISTP(NBANDMX),IFIRST(NBANDMX),
     4  ILAST(NBANDMX)
      CHARACTER*3 NEF(-1:1)
      CHARACTER*2 NAME(2),SLABL(-3:NSTATEMX)
      COMMON /DATABLK/FREQ,UFREQ,DFREQ,ZMASS,RSQMU,RSQMUP,RMUP,
     1  COUNTOT,NISTP,NFSTOT,NBANDTOT,AN,MN,IB,JP,JPP,EFP,EFPP,VP,VPP,
     2  FSBAND,NFS,IEP,IEPP,ISTP,IFIRST,ILAST, NAME,SLABL
c----------------------------------------------------------------------- 
      WRITE(8,600)
cc    WRITE(10,*)
      DO  I= IFIRST(IBB),ILAST(IBB)
          IF(IEP(IBB).GT.0) WRITE(8,602) VP(IBB),JP(I),NEF(EFP(I)),
     1         VPP(IBB),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                               DFREQ(I)/UFREQ(I)
cc        IF(IEP(IBB).GT.0) WRITE(10,622) VP(IBB),JP(I),NEF(EFP(I)),
cc   1         VPP(IBB),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
cc   2                                         (MN(J,ISTP(IBB)),J=1,2)
          IF(IEP(IBB).EQ.0) WRITE(8,602) VP(IBB),VPP(IBB),NEF(EFP(I)),
     1            JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
     2                                               DFREQ(I)/UFREQ(I)
cc        IF(IEP(IBB).EQ.0) WRITE(10,622) VP(IBB),VPP(IBB),NEF(EFP(I)),
cc   1            JP(I),JPP(I),NEF(EFPP(I)),FREQ(I),UFREQ(I),DFREQ(I),
cc   2                                               DFREQ(I)/UFREQ(I)
          ENDDO
      WRITE(8,604)
      RETURN                       
  600 FORMAT(1x,59('='),'  calc-obs'/ "   v'  j' p'",
     1  '  v"  j" p"    FREQ(obs)     u(FREQ)    calc-obs  /u(FREQ)'/
     2  1x,69('-'))
  602 FORMAT(2(2I4,A3),f14.6,2f12.6,f10.4)
cc622 FORMAT(2(i3,' &',i3,' &',A3,' &'),f12.5,'  &',f9.5,' &',f9.5,
cc   1  ' &',i5,' &',i4)
  604 FORMAT(1x,69('-'))
      END   
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12

