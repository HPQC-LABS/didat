c     Main routines for He-CO potential of Kirk Peterson.  For driver
c     functions, see hecokp2d.f and hecokp3d.f.  

c     The interpolation uses direct-product RKHS interpolation in R and
c     theta. (See Ho and Rabitz, JCP 104, 2584 (1996).) The RKHS uses
c     the distance-like kernel q[2,5] in R, and the angle-like kernel
c     for theta.  At each R, theta pair, it constructs a polynomial
c     least-squares fit (order selected by the calling routine) of
c     energies for CO internuclear distance.  The routine kpcoef
c     then returns the fit coefficients of that polynomial.

c     This routine reads in a data file of ab initio points, which
c     must be named kpdata.pot, and must follow a strict format.
c     It also needs the include file polyfit.fi which defines
c     a common block.   BLAS and LAPACK routines are used in several
c     places.

c     Authors: K. A. Peterson and G. C. McBane.
c     Updated 21 June 2002 GCM.

c     3/3/01 GCM
c     modified version of evalrkhs for Peterson He-CO potential.
c     This version returns the polynomial fit coefficients for the
c     interaction potential expanded as a polynomial in (r1-r1e), where
c     r1 is the C-O internuclear distance.  The order of polynomial expansion
c     is set on the initialization call as the value of the last
c     argument nnp_req.  A value of 1 returns the interaction potential
c     at r1=r1e, value of 2 returns a linear expansion, etc.  All work
c     so far has used nnp=4 to give a cubic expansion.

c     For an evaluation call, set nnp_req .le. 0 and supply x2(intermolecular
c     distance coordinate R, in bohr) and x3 (cos(theta)).  Fit coefficients
c     will be returned in vr1coef (units cm-1) such that 
c      V_{int}(r1, R, theta) = vr1coef(1) + vr1coef(2)*(r1-r1e) + 
c     vr1coef(3)*(r1-r1e)**2 + ... vr1coef(nnp_req)*(r1-r1e)**(nnp_req-1)
c     Calling program must ensure double precision array vr1coef has length
c     at least nnp_req; if this requirement is not met, the warranty is void.

c     Value of r1e (CO equilibrium interatomic distance, in bohr) must be set in 
c     common block /keep2/ before the initial call.



      subroutine kpcoef(x2,x3,vr1coef,nnp_req)
      implicit double precision  (a-h,o-z)
      double precision x2, x3, vr1coef(*)
      integer nnp_req

      PARAMETER(NN=50,maxpts=2000,nnn=nn*nn)
      common/rkhs/x1a(NN),x2a(NN),x3a(NN),m1,m2,m3
      common/indexr/ind2(maxpts),ind3(maxpts),n2
      common/keep/alp2,ya(nn,nnn),y2a(nn,nnn)
      common/keep2/r1e
      dimension ytmp(nnn),y2tmp(nnn)
      character*128 label

c   x2: r2
c   x3: cos(theta)
c
      if(nnp_req .gt. 0) then     !initialization call

         nin=15
c     
C**** generate 2-d RKHS coefficients from fort(nin) at first time*********** 
c...  Data is sorted on column 1, ascending order
c     
         open(unit=nin,file='kpdata.pot',status='OLD',form='FORMATTED')
c
c     open(unit=nin,file='cbs-corr.pot',status='OLD',form='FORMATTED')
c     open(unit=nin,file='cbs.pot',status='OLD',form='FORMATTED')
c     open(unit=nin,file='davdz.pot',status='OLD',form='FORMATTED')
c     open(unit=nin,file='davtz.pot',status='OLD',form='FORMATTED')
c     open(unit=nin,file='davqz.pot',status='OLD',form='FORMATTED')
c     
         read(nin,1000) label
         write(*,1000) label
 1000    format(A120)
         write(*, 1005) nnp_req, r1e
 1005    format(1x,'kpcoef: setting up for ', I2, '-term polynomial ',
     1        'fit of', /,'He-CO interaction potential, with r1e = ', 
     2        F9.5, ' bohr.')

         read (nin,*) m1,m2,m3
         n2=m2*m3
         np=0
         do 13 i=1,m1
            inn=0
            do 13 j=1,m2
               do 13 k=1,m3
                  inn=inn+1
                  np=np+1
                  read (nin,*) x1a(i),x2a(j),x3a(k),ya(i,inn)
                  ind2(np)=j
                  ind3(np)=k
 13      continue

         close(nin)
c     
         do 22 k=1,m3
            x3a(k)=(1.d0-cosd(x3a(k)))/2.d0 ! scale angle to interval [0,1]: angle kernel
 22      continue
c     
         do i=1,m1
            do j=1,n2
               ytmp(j)=ya(i,j)
            enddo
            call fit_rkhs2(i,ytmp,y2tmp)
            do j=1,n2
               y2a(i,j)=y2tmp(j)
            enddo
         enddo

c     set up polynomial fitting matrix for vibrational averaging
c     (results stored in polyfit common block for use by rkhs2coef)
         call gen_polyfit(x1a, m1, nnp_req)

      else  !  evaluation call
         
                 x3i=(1.d0-(x3))/2.d0   ! scale angle to interval [0,1]
         call rkhs2coef(x2,x3i,y2a, vr1coef)
      endif
      return
      end 
c

c  distance-like kernel q[2,5](x,y)
c  q[2,5](x,y) = (2/21)*xg**(-6)*(1-(3/4)*(xl/xg))
c  where xg is the greater of x and y and xl is the lesser
c  definition given in Ho and Rabitz, JCP 104, 2584 (1996), eqn. (17).
c  G. McBane 15 June 2002
c  
c  obscure coding intended for speed, though it might not be necessary
       double precision function rkhsdk(x,y)
       implicit none
       double precision x, y
       double precision oob7, a, b
       double precision to21, moo14
       parameter(to21 = 2.0d0/21.0d0, moo14 = -1.0d0/14.0d0)

       if (x .lt. y) then
          a = x
          b = y
       else
          a = y
          b = x
       end if
       oob7 = b*b*b
       oob7 = oob7*oob7
       oob7 = oob7*b
       oob7 = 1.0d0/oob7

       rkhsdk = oob7*(b*to21 + a*moo14)
       return
       end
c

c$$$c  here is Kirk's new rkhsdk, unmolested by gcm
c$$$        double precision function rkhsdk(x,xp)
c$$$c.. Evaluate the distance-like kernel of order nk with weight function
c$$$c.. x**(-mk), given (x,x')
c$$$c.. Currently, mk is limited to 6 and nk to 2.
c$$$c -----------------------------------------------------------------
c$$$        implicit real*8 (a-h, o-z)
c$$$        dimension F(0:9)
c$$$c.. Store a few factorials as floating point numbers.
c$$$        data (F(i),i=0,9)/1.0d0, 1.0d0, 2.0d0, 6.0d0, 24.0d0, 120.0d0,
c$$$     >                    720.0d0, 5040.0d0, 40320.0d0, 362880.0d0/
c$$$c
c$$$
c$$$c  use q[2,5] kernel to give long range behavior of C6*R**-6 + C7*R**-7
c$$$c  for He-CO
c$$$
c$$$        nk=2
c$$$        mk=5
c$$$
c$$$        if (x .le. xp) then
c$$$            z = x/xp
c$$$            A = 1.0d0/(xp**(mk+1))
c$$$        else
c$$$            z = xp/x
c$$$            A = 1.0d0/(x**(mk+1))
c$$$        endif
c$$$
c$$$c.. Beta function:
c$$$        B = F(nk-1)*F(mk)/F(nk+mk)
c$$$
c$$$c.. The special case of nk = 2.
c$$$        C = float(mk+1)/float(mk+3)
c$$$        rkhsdk = (nk**2)*A*B*(1.0d0 - C*z)
c$$$c
c$$$        return
c$$$        end
c$$$

c
c$$$       function rkhsak(x,y)
c$$$       implicit double precision (a-h,o-z)
c$$$       if((x.eq.0.d0).or.(y.eq.0.d0)) then
c$$$        rkhsak=1.d0
c$$$        return
c$$$       endif
c$$$       if(x.lt.y) then
c$$$          rkhsak=1.d0+x*y+2.d0*x*x*y*(1.d0-x/(3*y))
c$$$       else
c$$$          rkhsak=1.d0+x*y+2.d0*x*y*y*(1.d0-y/(3*x))
c$$$       endif
c$$$       return
c$$$       end

c  rewritten for speed

       double precision function rkhsak(x,y)
       implicit none
       double precision x, y
       double precision a, b
       double precision moo3
       parameter(moo3 = -1.0d0/3.0d0)
       
       if((x.eq.0.d0).or.(y.eq.0.d0)) then
          rkhsak=1.d0
          return  
       endif

       if (x .lt. y) then
          a = x
          b = y
       else
          a = y
          b = x
       end if

       rkhsak = 1.0d0 + a*(b+2.0d0*a*(b+moo3*a))
       return
       end

c  this version substantially modified
c  G. C. McBane 5/6 September 2000

       subroutine rkhs2coef(x2,x3,coef, vr1coef)

       implicit double precision (a-h,o-z)
       double precision vr1coef(*)
       parameter(nn=50,maxpts=2000,nnn=nn*nn)
       common/rkhs/r1(nn),r2(nn),theta(nn),m1,m2,m3
       common/indexr/ind2(maxpts),ind3(maxpts),n2
       common/keep2/r1e

c  Polyfit common block appears only here and in gen_polyfit
       include 'polyfit.fi'

       dimension coef(nn,nnn),e(nn)
       double precision kern_dp(nnn), kern_r2(nn), kern_theta(nn)
       integer index

c     get rkhs kernel values for r2 and theta (we have a rectangular
c     grid of ab points, so the positions of available (r2, theta)
c     points are the same no matter what r1 is.) There is still some
c     redundancy here; since x2 and x3 are the same on every call within
c     the loops, some work is being duplicated in the kernel functions
c     on calls after the first.  And since the r2 and theta vectors
c     never change, functions that involve only them are being
c     needlessly recalculated also. But I think this inefficiency is
c     minor.

       do j = 1, m2
          kern_r2(j) = rkhsdk(x2, r2(j))
       end do
       do j = 1, m3
          kern_theta(j) = rkhsak(x3, theta(j)) 
       end do

c     form direct product vector, arranged as it is in
c     the setup code (loop over R and theta with theta
c     innermost.) 
       index = 0
       do j = 1, m2
          do i = 1, m3
             index = index+1
             kern_dp(index)  = kern_r2(j)*kern_theta(i)
          end do
       end do

c     now generate e for different r1 at these values of R, theta.  
c     All the potential information is contained in coef matrix.

       call dgemv('N', m1, n2, 1.0d0, coef, nn, kern_dp, 1, 0.0d0, e, 1)

c     determine polynomial fit coefficients.  fitmat was set up
c     in gen_polyfit by SVD.
       
       call dgemv('N', nnp, m1, 1.0d0, fitmat, maxnnp, e, 1, 0.0d0, 
     1      vr1coef, 1)

       return
       end
c
c
c  fit to 2-d direct product of distance & angle kernels
       subroutine fit_rkhs2(n,y3,coef)
       implicit double precision (a-h,o-z)
       parameter(nn=50,maxpts=2000,nnn=nn*nn)
       parameter(lwork=nnn*20)
       dimension work(lwork),s(nnn)
       common/rkhs/x1a(NN),x2a(NN),x3a(NN),m1,m2,m3
       common/indexr/ind2(maxpts),ind3(maxpts),n2
       dimension coef(nnn),a(nnn,nnn),b(nnn),y3(nnn)
       data rcond/1.d-12/
c
       ii=(n-1)*n2
       do k=1,n2
          do l=1,n2
            a(k,l) =
     >               rkhsdk(x2a(ind2(ii+k)),x2a(ind2(ii+l)))*
     >               rkhsak(x3a(ind3(ii+k)),x3a(ind3(ii+l)))
          enddo
          b(k)=y3(k)
       enddo
c  do SVD fit and return coefficients
       call dgelss(n2,n2,1,a,nnn,b,nnn,s,rcond,irank,work,
     >              lwork,info)
c       print *, 'info= ',info,', nparms= ',n2,', rank= ',irank
c       write(6,10) rcond,s(1)/s(n2)
10     format('rcond= ',e9.1,', s(1)/s(np)= ',e9.1/)
c      write(35,'(e9.3)') (s(i),i=1,n2)
       do j=1,n2
        coef(j)=b(j)
       enddo
c
       return
       end
c

c  generate fitting matrix for polynomial fit in r1
c  by singular value decomposition.  Uses recommended sequence of
c  LAPACK calls for overdetermined system.
c  Begun 5 September 2000 G. C. McBane

       subroutine gen_polyfit(r, np, nnp_req)
       implicit none

c     r : vector of ab initio positions. 
c     np : length of r
c     nnp_req: requested order of polynomial fit in r1

       integer np, nnp_req
       double precision r(np)

       integer nn
       parameter(nn=50)

c     Polyfit common block is used for communication between this routine
c     and rkhs2coef
       include 'polyfit.fi'

c     fitmat : fitting matrix (effective inverse of design matrix
c     A).  Given a vector V of potential values at the np fitting
c     points, the coefficients of the best-fit polynomial through those
c     points will be fitmat*V.
c     nnp_req: the requested number of coefficients of the polynomial expansion.  
c     maxnnp: maximum number of coefficients for which storage
c     is available.
      
       common/keep2/r1e
       double precision r1e
       integer i, k, l
       integer  lwork

       parameter(lwork=nn*20)
       double precision work(lwork)
       double precision a(nn,maxnnp)
       double precision  d(maxnnp), e(maxnnp)
       double precision tauq(maxnnp), taup(maxnnp)
       double precision q(nn, maxnnp), pt(nn, maxnnp)
       integer info

       double precision diff, rcond, maxs
       data rcond/1.d-12/

       if (nnp_req .gt. maxnnp) then
          write(*, 1000) nnp_req, maxnnp
          stop
       else
          nnp = nnp_req
       end if
 1000  format(1x, 'nnp_req (',I2, ') too large in gen_polyfit.',
     1       '  Bailing out.', /, 
     2      1x, 'Change max_nnp in polyfit.fi and recompile.' )

c  form the (Vandermonde) design matrix A.
       do k = 1, np
          a(k,1) = 1.0d0
          diff = r(k) - r1e
          do l = 2, nnp
             a(k, l) = a(k, l-1)*diff
          end do
       end do

c     reduce A to band diagonal form

       call dgebrd(np, nnp, a, nn, d, e, tauq, taup,
     1      work, lwork, info)

c  copy a into q and pt (only one copy really necessary)
      do i = 1, nnp
         call dcopy(np, a(1,i), 1, q(1,i), 1)
         call dcopy(np, a(1,i), 1, pt(1,i), 1)
      end do 

c     generate first nnp columns of Q matrix 
      call dorgbr('Q', np, nnp, nnp, q, nn, tauq, work, lwork, info)
c     generate PT matrix (square, nnp x nnp)
      call dorgbr('P', nnp, nnp, np, pt, nn, taup, work, lwork, info)

c     perform SVD of band diagonal matrix. a is here as a dummy argument;
c     nothing will be done with it.  U and V^T matrices from SVD 
c     are returned in q and pt respectively.

      call dbdsqr('U', nnp, nnp, np, 0, d, e, pt, nn, q, nn, a, nn, 
     1     work, info)


c     scale columns of U (stored in q) by reciprocal singular values,
c     unless the singular values are too small, in which case we zero
c     those columns

      maxs = d(1)
      do i = 1, nnp
         if (d(i)/maxs .ge. rcond) then
            call dscal(np, 1.0d0/d(i), q(1,i), 1)
         else
            call dscal(np, 0.0d0, q(1,i), 1)
            write(*, 1010) i, d(i)/maxs
         end if
      end do

 1010 format('Note; polynomial fitting matrix is rank-deficient.  ',/,
     1     's(i)/smax for column ', I2, 'is ', G10.2)

c     now V*U^T gets us the effective inverse matrix

      call dgemm('T', 'T', nnp, np, nnp, 1.0d0, pt, nn,
     1     q, nn, 0.0d0, fitmat, maxnnp)



      return
      end





