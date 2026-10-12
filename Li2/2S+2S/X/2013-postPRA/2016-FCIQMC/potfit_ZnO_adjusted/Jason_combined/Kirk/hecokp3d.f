c     Evaluation routine for 3-dimensional version of Kirk Peterson's
c     He-CO ab initio surface.  This is a simple driver for the kpcoef
c     routine, which does almost all of the work.  For the 2d
c     vibrationally averaged version, see hecokp2d.f.

c     Arguments: 

c     RCO: CO internuclear distance, in bohr.

c     R: distance, in bohr, between He atom and CO center of
c     mass.  

C     CT: cosine of the angle between the vector pointing from CO center
C     of mass to He and the vector pointing from C to O.  CT=1
C     corresponds to the He-O-C linear arrangement.

c     Return value: He-CO interaction energy, in cm-1.  This value
c     does not include the CO monomer vibrational energy.

c     Author: G. C. McBane, 21 June 2002
c
      program driver
      implicit none
      double precision RCO, R, CT, result, hecokp3d

      RCO=2.1322d0
      R=6.418d0
      CT=0.328372127d0

      result=hecokp3d(RCO,R,CT)
      print *, RCO, R, acosd(CT), result

      end

c
      double precision function hecokp3d(RCO, R, CT)
      implicit none
      double precision RCO, R, CT

c     Common block keep2 communicates with kpcoef routine,
c     sets pivot for Taylor expansion in CO internuclear distance.

      common/keep2/r1e    
      double precision r1e
      data r1e /2.1322d0/
c     save /keep2/

c     set number of terms in polynomial expansion in r1 
c     (CO internuclear distance).  4 gives cubic polynomial.
c     This value may be decreased safely; increasing it might
c     require changes in array sizes in kpcoef.f and/or polyfit.fi.

       integer numterms
       parameter (numterms=4)

       integer ncall, i
       data ncall /1/
       save ncall

       double precision vr1coef(numterms), x

       if (ncall.eq.1) then
c     initialize
          call kpcoef(R, CT, vr1coef, numterms)
          ncall = 2
       end if

c     now evaluate potential coefficients and evaluate polynomial
c     in RCO

       call kpcoef(R, CT, vr1coef, -1)

       x = RCO-r1e
       hecokp3d = vr1coef(numterms)
       do i = numterms-1, 1, -1
          hecokp3d = x*hecokp3d+vr1coef(i)
       end do

       return
       end
