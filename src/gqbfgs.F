ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c note: all real variables in this routine are double precision
c so make sure that any variable transfered into the routine is
c double precision (64 bit)
      subroutine gqbfgs(nvar,x, iswtch,maxf,iter,iprint,tolx,tolg,
     1           rfn,sprec,extbnd,objf,ihess,hess,g,ierr,nfcall,srhvec,
     2           solvec,grdvec,scrvec,histry,irhess,iwhess,funct)
c
      use defaults
      Implicit None
c      implicit real(kind=RDbl) (a-h,o-z)
c
      Real(Kind=RDbl) :: x, hess, g, srhvec, solvec, grdvec, scrvec
      Real(Kind=RDbl) :: temp, objf, sqgrad, euclid, tobjf, z
      Real(Kind=RDbl) :: dirdev, redfcn, rfn, alpha, sprec, extbnd
      Real(Kind=RDbl) :: t, tolg, tolx, histry
      Integer         :: nvar, ihess, nfcall, iter, ierr, iii, irhess
      Integer         :: iprint, k, j, ii, i, ij, iconv, iswtch
      Integer         :: iexit, iwhess, maxf
      dimension  x(nvar),hess(ihess),g(nvar),srhvec(nvar),
     1           solvec(nvar),grdvec(nvar),scrvec(nvar)
c
      external funct
c
c     this subroutine minimizes a function of several variables
c     objf(x(1),x(2),...,x(nvar))  using a quasi-newton method
c     optionally employing the  dfp  and  bfgs  updating formulas.
c     the user must provide a subroutine to calculate objf(x) and
c     its gradient (first partial derivative) vector g(x).
c
c     **************************************************************
c     this routine invokes the package modules  search  and  uphess
c     and the user supplied subroutine  funct.
c     **************************************************************
c
c     on input:
c
c        nvar  is the number of variables. it is also the dimension of
c           the vectors  x, g, srhvec, solvec, grdvec  and  scrvec.
c
c        x  contains an estimate of the solution vector
c           (x(1),x(2),...,x(nvar)).
c
c        iswtch  is a parameter set equal to  k  which selects the
c           formula used to update the approximation to the hessian
c           inverse. for
c
c           k = 1 - the  dfp  update,
c           k = 2 - the  bfgs  update.
c
c           the  bfgs  update is recommended.
c
c        maxf  is the limit on the number of calls to the function
c           evaluation routine  funct.
c
c        funct  is a user supplied subroutine to evaluate  objf(x)
c           and the components of the gradient  g(x)  at the estimate
c           x(i), i = 1,2,...,nvar. Declared external in calling
c           routine.
c
c        tolx,tolg are the accuracies required in the solution, i.e. a
c           normal return from the routine occurs if the difference
c           between the components of two successive estimates of the
c           solution are not greater than  max(tolx*abs(x(i)),tolx)
c           for all  i, and the l2 norm of the gradient is not greater
c           than  tolg.
c
c        rfn  is an estimate of the expected reduction in  objf(x).
c           this estimate is used only on the first iteration so an
c           order of magnitude estimate will suffice. the information
c           can be provided in the following ways depending upon the
c           value of  rfn. for
c
c           rfn .gt. 0.0 - the setting of  rfn  itself will be taken
c                          as the expected reduction in  objf(x),
c
c           rfn = 0.0    - it is assumed that an estimate of the
c                          minimum value of  objf(x)  has been set
c                          in the argument  objf, and the expected
c                          reduction in  objf(x)  will be computed
c                          as (initial function value) minus objf,
c
c           rfn .lt. 0.0 - a multiple  abs(rfn)  of the modulus of
c                          the initial function value will be taken
c                          as the expected reduction.
c
c        sprec  is the accuracy required in the linear search technique
c           invoked by  gqbfgs, i.e. a point  xm  is accepted as the
c           minimum along the search direction if the ratio of the
c           directional derivative at  xm  over the directional
c           derivative at the initial point is not greater than  sprec.
c
c           the setting  0.100  is recommended.
c
c        extbnd  is the upper bound on the multiplicative increase in
c           the search scalar during the extrapolation phase of the
c           linear search technique.
c
c           the setting  2.000  is recommended.
c
c        objf  contains an estimate of the minimum value of  objf(x)
c           if  rfn = 0.0. otherwise it is only an output parameter.
c
c        ihess  is a parameter set equal to the dimension of  hess
c           which is at least  nvar*(nvar+1)/2.
c
c        histry  is a dummy parameter.
c
c     on output
c
c
c        x  contains the best available estimate of the solution vector.
c
c        objf  contains the function value at  x.
c
c
c        hess  is an array of dimension  ihess  which contains the
c           upper triangle of the most recent approximation to the
c           hessian inverse stored row-wise.
c
c        g  contains the components of the gradient at  x.
c
c        ierr  is a parameter set equal to  k  which gives the following
c           termination indications
c
c           normal termination,
c             k = 0,
c           intermediate termination,
c             k = -n(n any integer) - user termination,
c             k = 1 - failure to converge in  maxf  calls of  funct,
c             k = 2 - linear search technique indicates that it is
c                     likely that no minimum exists,
c
c        nfcall  is the number of calls to  funct.
c
c        srhvec  contains the current search direction vector.
c
c        solvec  contains the current solution vector.
c
c        grdvec  contains the current gradient vector.
c
c        scrvec  is a scratch vector.
c
c     written by k. e. hillstrom, march, 1976.
c
c
c     initialize the following parameters
c
c       ierr   - the termination indicator
c       nfcall - the number of calls to  funct
c       redfcn - the initial predicted reduction in  objf
c       iter   - the current iteration number
c
      ierr = 0
      iter=0
      nfcall = 1
      temp = objf
c
c     *****************************************************************
      call funct(nvar, x, objf, g )
c     *****************************************************************
c
      sqgrad=0.0
      do 220 iii=1,nvar
  220     sqgrad=sqgrad+g(iii)**2
c
      if (ierr .lt. 0) go to 410
      redfcn = rfn
      if (rfn .eq. 0.0) redfcn = objf - temp
      if (rfn .lt. 0.0) redfcn = abs(redfcn * objf)
      if (redfcn .le. 0.0) redfcn = 1.0
c
c    read initial estimate of hessian inverse, if desired
      if(irhess.eq.0) go to 200
      read(11,202) (hess(ii),ii=1,ihess)
  202   format(5e16.9)
      go to 300
c
c     begin the quasi-newton process by initializing the approximation
c       to the hessian inverse to unity
c
  200 if(iprint.ne.0) write(10,201)
  201 format(//2x,'>>>> quasi-newton procedure started, with search',
     1      ' direction set to -g')
      k = 1 + nvar * (nvar+1) / 2
c
      do 210 i = 1, nvar
c
         do 205 j = 1, i
            k = k - 1
            hess(k) = 0.0
  205    continue
c
         hess(k) = 1.0
  210 continue
c
  300 if(iprint.eq.0) go to 301
      euclid=sqrt(sqgrad)
      if(iprint.ge.20.and.mod(iter,10).ne.0) go to 308
      write(6,307) iter,nfcall,objf,euclid
  307 format(1x,'iteration:',i5,26x,'function evaluation:',i6,
     1      /1x,'objective function:',e20.13,2x,'gradient norme:',
     2       e20.13)
  308 if(mod(iter,iprint).ne.0) go to 301
      write(10,302) iter,nfcall
  302 format(////2x,'iteration no ',i5//2x,'number of function and ',
     1       'gradient evaluations = ',i5//2x,'parameter values')
      write(10,303) (j,x(j),j=1,nvar)
  303 format(/3(2x,'x(',i4,') = ',e16.8))
      write(10,304) objf
  304 format(//2x,'function value objf = ',e16.8///2x,'gradient')
      write(10,306) (j,g(j),j=1,nvar)
  306 format(/3(2x,'g(',i4,') = ',e16.8))
  301 iter=iter+1
c
c     begin an iteration by saving the current best estimate of the
c       function and the solution and gradient vectors.
c
      do 310 i = 1, nvar
         solvec(i) = x(i)
         grdvec(i) = g(i)
  310 continue
c
      tobjf = objf
c
c     calculate the search direction vector in  srhvec  and the
c       directional derivative in  dirdev
c
      do 340 i = 1, nvar
         ij = i
         z = 0.0
c
         do 330 j = 1, nvar
            z = z - g(j) * hess(ij)
            if (j .ge. i) go to 325
            ij = ij + nvar - j
            go to 330
  325       ij = ij + 1
  330    continue
c
         srhvec(i) = z
  340 continue
c
      dirdev = 0.0
c
      do 350 i = 1, nvar
         dirdev = dirdev + srhvec(i) * g(i)
  350 continue
c
c     if the directional derivative  dirdev  is .gt. 0, there is no
c       guarantee that a search in the  w  direction  will result in a
c       smaller  objf. therefore, the  quasi-newton  process is
c       restarted at the current estimate of the solution with  srhvec
c       set to  -g.
c
      if (dirdev .gt. 0.0) go to 200
      if (dirdev .eq. 0.0) go to 500
c
c     compute the initial search scalar  alpha  and conduct the
c       linear search by means of a call to  search
c
      alpha = -2.0 * redfcn / dirdev
      if (alpha .gt. 1.0) alpha = 1.0
      redfcn = objf
c
c     **************************************************************
      call search(nvar,x, g,srhvec,objf,alpha,dirdev,sprec,
     1            extbnd,nfcall,scrvec,ierr,funct,maxf)
c     **************************************************************
c
c     test for abnormal termination
c
      if (ierr .lt. 0) go to 500
      if (nfcall .ge. maxf) go to 400
      if ((alpha .lt. 1.0e-20) .or.
     1    (alpha .gt. 1.0e20)) go to 410
c
c     test for convergence
c
      sqgrad = 0.0
      iconv = 0
c
      do 360 i = 1, nvar
         temp = alpha * srhvec(i)
         sqgrad = sqgrad + g(i) * g(i)
         t = tolx * abs(x(i))
         if (t .le. tolx) t = tolx
         if (abs(temp) .gt. t) iconv = 1
  360 continue
c
      if (sqgrad .gt. tolg*tolg) iconv = 1
      if (sqgrad .eq. 0.0) iconv = 0
      if (iconv .eq. 0) go to 500
c
c     the linear search technique has located a minimum. call  uphess
c       to update the approximation to the hessian inverse using the
c       dfp  or  bfgs  updating formulas
c
c     ***************************************************************
      call uphess(nvar,x,g,ihess,hess,solvec,grdvec,scrvec,iswtch,iexit)
c     ***************************************************************
c
c     if the update is not successful the quasi-newton process is
c       restarted with a descent step at the current estimate of the
c       solution
c
      redfcn = redfcn - objf
      if (iexit .ne. 0) go to 200
c
c     now start a new iteration
c
c  write the current estimate of the hessian inverse, if desired.
      if(iwhess.eq.0) go to 300
      write(12,202) (hess(ii),ii=1,ihess)
c
      go to 300
c
c     error return because there have been at least  maxf  calls of
c       funct
c
  400 ierr = 1
      go to 450
c
c     error return because linear search technique indicates that it is
c       likely that no minimum exists
c
  410 ierr = 2
c
  450 do 455 i = 1, nvar
         x(i) = solvec(i)
         g(i) = grdvec(i)
  455 continue
c
      objf = tobjf
  500 return
c
      end

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine uphess(n,x,g,ih,h,solvec,grdvec,scrvec,iswtch,iexit)
c
      Use defaults, Only: RDbl

      Implicit None
c      implicit real*8 (a-h,o-z)
c
      Real(Kind=RDbl) :: x, g, h, solvec, grdvec, scrvec, z, alpha
      Real(Kind=RDbl) :: wt
      Integer :: n, ih, iexit, i, k, j, iswtch
      dimension  x(n),g(n),h(ih),solvec(n),
     1           grdvec(n),scrvec(n)
c
c     this subroutine updates an approximation to the hessian inverse
c       using the  dfp  or  bfgs  formula
c
c     on input
c
c        n  is the dimension of the vectors  x, g, solvec, grdvec  and
c           scrvec.
c
c        x  contains an estimate of the solution vector.
c
c        g  contains the components of the gradient corresponding to
c           the  x  vector.
c
c        ih  is a parameter set equal to the dimension of  h  which is
c           at least  n*(n+1)/2.
c
c        h  is an array of dimension  ih  which contains the upper
c           triangle of an approximation to the hessian inverse stored
c           by rows.
c
c        solvec  contains the current solution vector.
c
c        grdvec  contains the current gradient vector.
c
c        iswtch  is a parameter set equal to  k  which selects the
c           updating formula. for
c
c           k = 1 - the  dfp  formula is used,
c           k = 2 - the  bfgs  formula is used.
c
c     on output
c
c        iexit  is a parameter set equal to  k  which indicates the
c           following. for
c
c           k = 0 - the update was successful,
c           k = 1 - the update failed due to zero divisors.
c
c        h  contains the updated approximation to the hessian inverse
c           if  iexit = 0.
c
c        scrvec  is a scratch vector.
c
c     written by k. e. hillstrom, march, 1976.
c
c
c     initialize the exit indicator  iexit
c
      iexit = 0
c
c     calculate the solution and gradient difference vectors from two
c       consecutive iterations. from this section on
c
c       solvec - contains  delta, the solution difference vector
c       grdvec - contains  gamma, the gradient difference vector
c
  100 do 110 i = 1, n
         solvec(i) = x(i) - solvec(i)
         grdvec(i) = g(i) - grdvec(i)
  110 continue
c
c     calculate  z = (gamma transpose) * delta and alpha =
c       (gamma transpose) * (hessian inverse) * gamma  occuring as
c       denominators in the  dfp  formula. from this section on
c
c       h      - contains the approximation to the hessian inverse
c       scrvec - contains the successive elements of  (gamma transpose)
c                * (hessian inverse)
c
      z = 0.0
      alpha = 0.0
c
      do 130 i = 1, n
         wt = grdvec(i)
         z = z + wt * solvec(i)
         k = i
         wt = 0.0
c
         do 120 j = 1, n
            wt = wt + grdvec(j) * h(k)
            if (j .ge. i) go to 115
            k = k + n - j
            go to 120
  115       k = k + 1
  120    continue
c
         alpha = alpha + wt * grdvec(i)
         scrvec(i) = wt
  130 continue
c
c     error exit if the  dfp  or  bfgs  formula breaks down due to zero
c       divisors  z  and/or  alpha
c
      if ((z .eq. 0.0) .or.
     1    (alpha .eq. 0.0 .and. iswtch .eq. 1)) go to 200
c
c     update the approximation to the hessian inverse using the  dfp
c       or  bfgs  updating formula
c
      k = 1
c
      do 160 i = 1, n
c
         do 150 j = i, n
            if (iswtch .eq. 1) go to 135
            h(k) = h(k) - (solvec(i) * scrvec(j) + scrvec(i) *
     1             solvec(j)) / z + (1.0 + alpha / z) * (solvec(i) *
     2             solvec(j) / z)
            go to 140
  135       h(k) = h(k) + solvec(i) * solvec(j) / z - scrvec(i) *
     1             scrvec(j) / alpha
  140       k = k + 1
  150    continue
c
  160 continue
c
      go to 300
c
c     error return due to zero divisors in the updating formula
c
  200 iexit = 1
  300 return
c
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine search(n,x,g,s,f,alpha,dirdev,sprec,extbnd,
     1                  nfcall,w,ierr,funct, maxf)
c
      Use defaults, Only: RDbl
      Implicit None
c      implicit real*8 (a-h,o-z)
c
      Real(Kind=RDbl) :: x, g, s, w, tot, cdirev, dirdev, pdirev
      Real(Kind=RDbl) :: alpha, ftest, f, sprec, temp, extbnd, wt
      Real(Kind=RDbl) :: ww
      Integer :: n, i, nfcall, maxf, ierr
      dimension x(n),g(n),s(n),w(n)
c
c       this subroutine obtains an estimate of the solution
c     xm = x0 + alpha * s  which minimizes  f  by means of a linear
c     search in the  s  direction.
c
c     on input
c
c        n  is the number of variables. it is also the dimension of the
c           vectors  x, g  and  s.
c
c        x  contains an estimate of the solution vector
c           (x0(1),x0(2),...,x0(n)).
c
c        s  contains the search direction vector.
c
c        f  contains the objective function  f(x).
c
c        funct  is a user supplied routine to evaluate  f(x)  and the
c           components of the gradient  g(x)  at the estimate in  x.
c
c        alpha  is  the initial step scalar.
c
c        dirdev  is the directional derivative at x.
c
c        sprec  is the accuracy required in the search, i.e. a point  xm
c           is accepted as the minimum along direction  s  if the ratio
c           of the directional derivative at  xm  over the directional
c           derivative at  x0  is not greater than  sprec.
c
c        extbnd  is the upper bound on the multiplicative increase in
c           alpha  during extrapolation.
c
c        nfcall  is the number of calls to the function evaluation
c           subroutine  fcn.
c
c        w  is a scratch vector.
c
c     on output
c
c        x  contains the estimate of the minimum
c           (xm(1),xm(2),...,xm(n))
c
c     ierr  is a parameter set to a negative integer if the user
c           wishes to force an exit from  search. otherwise it is
c           unaltered.
c
c        g  contains the components of the gradient at  x.
c
c        f  contains the function value  f(x).
c
c        alpha  is the final step scalar.
c
c        dirdev  is the directional derivative at  x.
c
c        nfcall  is the number of calls to the function evaluation
c           subroutine  funct.
c
c
c     initialize the following parameters and indicators
c
c       tot       -  the sum of the extrapolation steps
c       cdirev    -  the current directional derivative
c       pdirev    -  the previous directional derivative
c       ierr      -  the error indicator
c
      tot = 0.0
      cdirev = dirdev
      pdirev = dirdev
c
c     test whether  alpha  is too small
c
  105 if (alpha .le. 1.0e-20) go to 150
c
c     begin the linear search by incrementing the solution vector  x
c       and calculating the function and gradient at the incremented  x.
c
      do 108 i = 1, n
         w(i) = x(i)
         x(i) = x(i) + alpha * s(i)
  108 continue
c
c     ***************************************************************
      call funct(n,x, ftest,g)
c     ***************************************************************

c
      nfcall = nfcall + 1
      if (maxf.lt.nfcall) go to 160
      if (ierr .lt. 0) go to 150
c
c     compute the directional derivative  dirdev  at  x + alpha * s
c
      dirdev = 0.0
c
      do 110 i = 1, n
         dirdev = dirdev + g(i) * s(i)
  110 continue
c
c     test whether  f(x + alpha * s)  is less than  f(x).
c
      if (ftest .ge. f) go to 120
c
c     if (dirdev / pdirev) is less than the search precision  sprec,
c       alpha  is accepted. otherwise  alpha  is modified
c
      if (abs(dirdev / pdirev) .le. sprec) go to 140
c
c     alpha  is modified. test whether  alpha  is to be revised by
c       extrapolation or interpolation
c
      if (dirdev .gt. 0.0) go to 120
c
c     alpha  is revised using an extrapolation formula and a new step
c       is taken if the sum of the steps already made is not too large.
c       the input parameter  extbnd  limits the multiplicative change
c       in  alpha
c
      tot = tot + alpha
      if (tot .gt. 1.0e10) go to 145
      temp = extbnd
      if (cdirev .lt. dirdev) temp = dirdev / (cdirev - dirdev)
      if (temp .gt. extbnd) temp = extbnd
      f = ftest
      cdirev = dirdev
      alpha = alpha * temp
      go to 105
c
c     x  is reset to the current estimate, alpha  is revised using the
c       cubic interpolation formula and a new step is taken if the
c       convergence criteria have not been satisfied.
c
  120 do 130 i = 1, n
         x(i) = w(i)
  130 continue
c
      temp = 3.0 * (f - ftest) / alpha + dirdev + cdirev
      wt = abs(temp)
      if(wt.lt.abs(dirdev)) wt=abs(dirdev)
      if(wt.lt.abs(cdirev)) wt=abs(cdirev)
      ww = temp / wt
      ww = ww * ww - cdirev / wt * dirdev / wt
      if (ww .lt. 0.0) ww = 0.0
      ww = sqrt(ww) * wt
      temp = 1.0 - (dirdev + ww - temp) / (2.0 * ww + dirdev -
     1       cdirev)
      alpha = alpha * temp
      go to 105
c
c     alpha is accepted
c
  140 f = ftest
  145 alpha = tot + alpha
  150 return
c
  160 do 170 i = 1, n
         x(i) = w(i)
  170 continue

      return
c
      end
