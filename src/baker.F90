!---------------------------------------------------------------------
! This module is an implementation of the classic Baker's algorithm
! Ref: Journal of Computation Chemistry, vol. 7, 385--395, 1986.
!---------------------------------------------------------------------
Module baker
  Use defaults, Only: dashedline, RDbl, lstrLen, strLen
  Use utils, Only: filesrchstr, isfileopen, getMachPrec
  Use file, Only: file_getunit
  Use config, Only: AtMolCoords
  Use simcell, Only: Simcell_Params
  Use bakersubs, Only: bakersubs_calcgrad, bakersubs_calchess

  Implicit None
  Save

  Private
  Public :: Baker_Params, baker_isSaddleSearch, baker_display, baker_baker, &
      baker_init

  Type BAKER_Params
    Integer        :: isaddle ! Find the ts(1) or minima(0)
    Integer        :: ndf     ! Degrees of freedom
    Integer        :: maxstep ! Maximum no. of iterations
    Integer        :: nchess  ! Calc. Hessian every nchess steps
    Integer        :: dof     ! Degrees of freedom
    Integer        :: lastiterunit ! File for dumping for debug purposes
    Real(kind=RDbl):: hmax
    Real(kind=RDbl):: toldf, tolgrad
  End Type BAKER_Params
  
  Character(len=strLen)    :: default_baker_tag = &
      "Baker Algorithm Parameters"
  
Contains
  !------------------------------------------------------------------
  ! Initialize the parameters for Baker's Algorithm from the control
  ! file "ctrl_file" and store them in the "bakersparams" object
  !------------------------------------------------------------------
  Subroutine baker_init(bakerparams, ctrl_filename, opt_bakertag, optlineno)
    Type(BAKER_Params), Intent(out)  :: bakerparams
    Character(*), Intent(in)         :: ctrl_filename
    Character(*), Intent(in),Optional:: opt_bakertag
    Integer, Intent(in), Optional    :: optlineno

    Character(len=strLen)      :: tag
    Character(len=lstrLen)     :: line
    Integer                    :: lineno, unitno
    Real(kind=RDbl)            :: prec
    
    If (Present(opt_bakertag)) Then
      tag = opt_bakertag
    Else
      tag = default_baker_tag
    End If
    
    !** Open the ctrl_file if it is not opened
    unitno = isfileopen(ctrl_filename)
    If (unitno < 0) Then
      unitno = file_getunit(ctrl_filename)
      Open(file=ctrl_filename, unit=unitno)
    Endif
    
    !** Find the MEP section
    lineno = filesrchstr(unitno, tag, line)
    If (lineno == 0) Then
      Write(0,'(1x,2a,i4,3a)') __FILE__," : ",__LINE__, &
          " Could not find the tag ", tag, " in the control file"
      Stop
    Endif

    !** If the optional line number is present it implies that the init
    !** routine was called from within another section and we want to make
    !** sure that the tag it finds is different from the tag in the calling
    !** section
    If (Present(optlineno)) Then
      If (lineno > optlineno) Then
        Write(0,'(1x,2a,i4, 3a)') __FILE__," : ",__LINE__, &
            " Could not find the section '", Trim(tag), "'"
        Write(0,'(1x,2a)') "Either it is not defined in the control file or", &
            " it is defined AFTER the section referring to it"
        Stop
      End If
    End If
    
    !** Read the rest of the stuff
    Read(unitno,*) bakerparams%isaddle
    Read(unitno,*) bakerparams%maxstep
    Read(unitno,*) bakerparams%toldf
    Read(unitno,*) bakerparams%tolgrad
    Read(unitno,*) bakerparams%hmax
    Read(unitno,*) bakerparams%nchess

    !** Make sure the input make sense.  Particularly, check for precision
    ! Since the tolgrad is squared in the code we should check its square
    ! against the value of 'prec'
    prec = getMachPrec(bakerparams%tolgrad)
    If (bakerparams%tolgrad**2 < prec) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " 'tolgrad' has a precision lower than machine precision"
      Stop
    End If

    prec = getMachPrec(bakerparams%toldf)
    If (bakerparams%toldf < prec) Then
      Write(0,'(1x,2a,i4)') __FILE__," : ",__LINE__, &
          " 'toldf' has a precision lower than machine precision"
      Stop
    End If
  End Subroutine baker_init

  !--------------------------------------------------------
  ! Returns true if the parameters in "bakerparams" are set
  ! up to do a saddle point search
  !--------------------------------------------------------
  Logical Function baker_isSaddleSearch(bakerparams)
    Type(BAKER_Params), Intent(in)  :: bakerparams

    If (bakerparams%isaddle == 0) Then
      baker_isSaddleSearch = .False.
    Else
      baker_isSaddleSearch = .True.
    End If
  End Function baker_isSaddleSearch

  !----------------------------------------------------------------------
  !Baker's Algorithm subroutine.
  !Received (by mlg) from J. Scott Shaffer   Thu Apr 22 13:53:23 PDT 1993
  !Then passed on (unmodified) to me (Randy Snurr)
  !From Randy to Amit on  Tue Mar 14 14:31:52 CST 2000
  !----------------------------------------------------------------------
  !BAKER: This subroutine searches for stationary points (either minima
  !or first-order saddle points) on multi-dimensional potential energy
  !surfaces, using the algorithm described by 
  !Baker, J. Comp. Chem. 7, 385 (1986),
  !based on the rational function optimization method developed by
  !
  !Banerjee et al., J. Phys. Chem. 89, 52 (1985).
  !
  !The routine converges either to a minimum or a first-order saddle
  !point according to the flag ISADDLE, which will be zero for a minimum
  !search, or one for a transition-state search.
  !The user must supply the names of two potential function routines:
  !GRADIENT_SUB, which calculates the energy and gradient; and
  !HESSIAN_SUB, which calculates the energy, gradient, and Hessian.
  !
  !INPUT
  !    isaddle..........This flag signals the type of stationary point to
  !                     be located:
  !
  !                        isaddle = 0, search for a minimum
  !                        isaddle = 1, search for a first-order saddle
  !
  !    ndf..............Number of degrees of freedom in the potential
  !                     energy surface
  !    df...............Starting point for the degrees of freedom
  !    toldf,
  !    tolgrad..........Tolerances on the step size and norm of the
  !                     gradient at the saddle point
  !    hmax.............Maximum step magnitude allowed at each iteration
  !    maxstep..........Maximum number of steps allowed
  !    nchess...........The Hessian will be recalculated (as opposed to
  !                     just updated approximately) every nchess steps
  !    fv1..............Workspace for the matrix diagonalization
  !    
  !OUTPUT
  !    df.............Stationary point on the potential energy surface
  !    objf,
  !    grad,
  !    eigval...........Objective function, gradient vector, and
  !                     eigenvalues of the Hessian matrix at the
  !                     stationary point 
  !    eigvec..........The kth column of this matrix contains the kth
  !                     eigenvector at the stationary point
  !    hess.............The Hessian matrix at the stationary point
  !    ierr.............Error condition: 
  !
  !                        ierr = 0, successful location of stationary
  !                                  point 
  !                        ierr = 1, maximum number of steps exceeded
  !                                  without locating the stationary
  !                                  point 
  !                        ierr = 2, error in Hessian diagonalization
  !                        ierr = 3, error in Newton's method solution
  !                                  for eigenvalue shift parameter
  !                        ierr = 4, error in energy calculation
  !----------------------------------------------------------------------
  Subroutine baker_baker(bakerparams, sorbates, simcell, mtypes, &
      df, objf, grad, eigval, eigvec, hess, nstep, ierr)
    Type(BAKER_Params), Intent(in)  :: bakerparams
    Type(AtMolCoords), Dimension(:), Intent(inout) :: sorbates
    Type(Simcell_Params), Intent(inout) :: simcell
    Integer, Intent(in), Dimension(:)   :: mtypes
    Real(kind=RDbl)                 :: objf
    Real(kind=RDbl), Dimension(:)   :: df, grad, eigval
    Real(kind=RDbl), Dimension(:,:) :: eigvec, hess
    Integer, Intent(out)            :: ierr
    
    ! -- LOCAL ------------
    Integer          :: ndf, maxstep, nchess, isaddle
    Real(kind=RDbl)  :: toldf, tolgrad, zero, small, tolnr, hmax
    Real(kind=RDbl), Dimension(bakerparams%dof) :: fv1
    Integer          :: itnrmax, lastiterunit
    Parameter        (itnrmax = 30)
    Parameter        (zero = 0.0d0)
    Parameter        (small = 0.0001d0)
    Parameter        (tolnr = 1.0d-08)
    
    Integer         :: idf, ierrdiag, itnr, jdf, nstep, i, j, err
    Real(kind=RDbl) :: derivlam, flam, gnormsqr, hcompmax, hmagsqr
    Real(kind=RDbl) :: hTh, lambdan, lambdap, scaleh, scalevec, sum1, sum2,VTh 
    Real(kind=RDbl), Dimension(bakerparams%dof) :: F, gradold, h, V
    Logical         :: conv, convnr
    
    !** Get some sundry information and make sure we don't have any
    !** absurd values
    isaddle = bakerparams%isaddle
    maxstep = bakerparams%maxstep
    toldf   = bakerparams%toldf
    tolgrad = bakerparams%tolgrad
    hmax    = bakerparams%hmax
    nchess  = bakerparams%nchess
    ndf     = bakerparams%dof
    lastiterunit = bakerparams%lastiterunit

    If (hmax < 1.0e-6_RDbl) Then
      Write(0,'(1x,2a,i4,a)') __FILE__," : ",__LINE__, &
          " The value of 'hmax' is too small"
      Stop
    End If

    !==========
    !Initialize
    !==========
    nstep = 0
    conv  = .False.
    ierr  = 1
    
    Call bakersubs_calcGrad(sorbates, simcell, mtypes, df, objf, grad, err) 
    If (err /= 0) Then
      ierr = 4
      Return
    End If
    Call bakersubs_calcHess(sorbates, simcell, mtypes, df, hess, err)
    If (err /= 0) Then
      ierr = 4
      Return
    End If

!!$    Write(lastiterunit,'(a,6f10.6)') 'grad: ', grad
!!$    Do i=1, ndf
!!$      Do j=1, ndf
!!$        Write(lastiterunit,'(f10.4)', Advance='No') hess(i,j)
!!$      End Do
!!$      Write(lastiterunit,*)
!!$    End Do
    !=================================================
    !Iterate to convergence or maximum number of steps
    !=================================================
    Do 
      If ((nstep > maxstep) .Or. (conv)) Exit
      
      nstep = nstep + 1
      
      !-------------------------------
      !Diagonalize the Hessian
      !-------------------------------
      !** eigval has the eigenvalues in ASCENDING order
      Call tred2 (ndf, ndf, hess, eigval, fv1, eigvec)
      Call tql2 (ndf, ndf, eigval, fv1, eigvec, ierrdiag)
!!$      Write(lastiterunit,'(a,6f8.3)') 'eigenvalues :', eigval
!!$      Do i=1, ndf
!!$        Write(lastiterunit,'(6f8.3)') eigvec(i, 1:6)
!!$      End Do
      
      !<< Bail out if diagonalization fails >>
      If (ierrdiag == 1) Then
        conv = .False.
        ierr = 2
        Return
      Endif

      !------------------------------------------------------------
      !Transform the gradient into the local Hessian modes
      !------------------------------------------------------------
      !** Convert the gradient vector into "eigenvector" coordinate
      !** system
      Do jdf = 1, ndf
        F(jdf) = 0.0_RDbl
        
        Do idf = 1, ndf
          F(jdf) = F(jdf) + eigvec(idf,jdf)*grad(idf)
        Enddo
      Enddo
      
      !=====================================================================
      !This is where the minimum search differs from the saddle point search
      !=====================================================================
      If (isaddle == 0) Then
        
        !--------------
        !Minimum search
        !--------------
        !<< CALCULATE THE SHIFT PARAMETER LAMBDA_N >>
        
        ! This must be found iteratively; Newton's method is coded below.  We
        ! can be confident that Newton's method will converge, because the root
        ! we seek is bounded from above by zero and the first eigenvalue.
        
        !<< Initial guess >>
        If (eigval(1) .Ge. small) Then
          lambdan = zero
        Else
          lambdan = eigval(1) - small
        Endif
        
        !<< Start Newton's method iterations >>
        itnr = 1
        convnr = .False.
        Do 
          If ((itnr > itnrmax) .Or. (convnr)) Exit
          sum1 = zero
          sum2 = zero
          
          Do idf = 1, ndf
            sum1 = sum1 + F(idf)*F(idf)/(lambdan - eigval(idf))
            sum2 = sum2 + F(idf)*F(idf)/((lambdan - eigval(idf))**2)
          Enddo
          
          flam = lambdan - sum1
          derivlam = 1.0d0 + sum2
          
          !<< Check for convergence >>
          If (flam .Le. tolnr) Then
            convnr = .True.
          Else
            lambdan = lambdan - flam/derivlam
          Endif
          
          itnr = itnr + 1
        End Do
        
        If (.Not. convnr) Then
          ierr = 3
          Return
        Endif
        
        !<< CONSTRUCT THE STEP VECTOR >>
        Do idf = 1, ndf
          h(idf) = 0.0d0
        Enddo
        
        Do jdf = 1, ndf
          scalevec = -F(jdf)/(eigval(jdf) - lambdan)
          Do idf = 1, ndf
            h(idf) = h(idf) + scalevec*eigvec(idf,jdf)
          Enddo
        Enddo
      Else
        ! -------------------
        ! Saddle-point search
        ! -------------------
        !<< CALCULATE THE SHIFT PARAMETER LAMBDA_P >>
        lambdap = 0.5d0*(eigval(1) + &
            Sqrt(eigval(1)*eigval(1) + 4.d0*F(1)*F(1))) 
        
        !<< CALCULATE THE SHIFT PARAMETER LAMBDA_N >>
        ! This must be found iteratively; Newton's method is coded below.  We
        ! can be confident that Newton's method will converge, because the root
        ! we seek must be less than zero, and must be bracketed by the first
        ! and second eigenvalues.
        
        !<< Initial guess >>
        ! The initial guess is zero, unless the second eigenvector is negative
        ! (or close enough to zero to cause numerical problems), in which case
        ! the initial guess is halfway between the first and second
        ! eigenvalues... 
      
        If (eigval(2) >= small) Then
          lambdan = zero
        Else
          lambdan = 0.5d0*(eigval(1) + eigval(2))
        End If
        
        !<< Start Newton's method iterations >>
        itnr = 1
        convnr = .False.

        Do 
          If ((itnr > itnrmax) .Or. (convnr)) Exit
          sum1 = zero
          sum2 = zero
          
          Do idf = 1, ndf
            sum1 = sum1 + F(idf)*F(idf)/(lambdan - eigval(idf))
            sum2 = sum2 + F(idf)*F(idf)/ ((lambdan - eigval(idf))**2)
          Enddo
          
          flam = lambdan - sum1
          derivlam = 1.0d0 + sum2
          
          !<< Check for convergence >>
          If (flam <= tolnr) Then
            convnr = .True.
          Else
            lambdan = lambdan - flam/derivlam
          Endif
          
          itnr = itnr + 1
        End Do
        
        If (.Not. convnr) Then
          ierr = 3
          Return
        End If
        
        !<< CONSTRUCT THE STEP VECTOR >>
        Write(lastiterunit,*) __LINE__, eigval(1), lambdap
        scalevec = -F(1)/(eigval(1) - lambdap)

        Do idf = 1, ndf
          h(idf) = scalevec*eigvec(idf,1)
        Enddo
        
        Do jdf = 2, ndf
          scalevec = -F(jdf)/(eigval(jdf) - lambdan)
          Do idf = 1, ndf
            h(idf) = h(idf) + scalevec*eigvec(idf,jdf)
          Enddo
        End Do
      End If     ! Do a minima search or a saddle-pt search
      
      !-------------------------------
      !Check the magnitude of the step
      !-------------------------------
      hmagsqr = 0.0_RDbl
      Do idf = 1, ndf
        hmagsqr = hmagsqr + h(idf)**2
      End Do
      
      If (hmagsqr > hmax*hmax) Then
        scaleh = hmax/Sqrt(hmagsqr)

        Do idf = 1, ndf
          h(idf) = scaleh*h(idf)
        End Do
      End If
      
      !-------------
      !Take the step
      !-------------
      Do idf = 1, ndf
        df(idf) = df(idf) + h(idf)
      End Do
      
      !---------------------
      !Check for convergence
      !---------------------
      !<< Gradient >>
      gnormsqr = zero
      
      Do idf = 1, ndf
        gnormsqr = gnormsqr + grad(idf)**2
      Enddo
      Write(lastiterunit,'(i5,a,i7,f11.5)') &
          __LINE__, "step, gnormsqr : ", nstep, gnormsqr
      
      !<< Degrees of freedom >>
      hcompmax = h(1)
      Do idf = 2, ndf
        hcompmax = dmax1(hcompmax,dabs(h(idf)))
      End Do

      If ((gnormsqr <= tolgrad*tolgrad) .And. (hcompmax <= toldf)) Then
        conv = .True.
        ierr = 0
      Endif
      
      !---------------------------
      !UPDATE GRADIENT AND HESSIAN
      !---------------------------
      If (.Not. conv) Then
        
        Do idf = 1, ndf
          gradold(idf) = grad(idf)
        End Do
        
        Call bakersubs_calcGrad(sorbates, simcell, mtypes, df, objf, grad, err)
        Write(lastiterunit,'(i5, a)', Advance='No') __LINE__, ' grad: ' 
        Do j=1, ndf
          Write(lastiterunit,'(f11.5)', Advance='No') grad(j)
        End Do
        Write(lastiterunit,*)

        If (err /= 0) Then
          ierr = 4
          Return
        End If
        
        If (Mod(nstep,nchess) .Eq. 0) Then
          ! << Exact calculation of Hessian >>
          Call bakersubs_calcHess(sorbates, simcell, mtypes, df, hess, err)
          If (err /= 0) Then
            ierr = 4
            Return
          End If

        Else
          !<< Numerical updating scheme from Powell >>
          !   The vector V is defined by Baker...
          Do idf = 1, ndf
            V(idf) = grad(idf) - gradold(idf)
            
            Do jdf = 1, ndf
              V(idf) = V(idf) - hess(idf,jdf)*h(jdf)
            Enddo
          Enddo
          
          hTh = zero
          VTh = zero
          Do idf = 1, ndf
            hTh = hTh + h(idf)*h(idf)
            VTh = VTh + V(idf)*h(idf)
          Enddo
          
          Do jdf = 1, ndf
            Do idf = 1, ndf
              hess(idf,jdf) = hess(idf,jdf)+ (V(idf)*h(jdf) + &
                  V(jdf)*h(idf) - VTh*h(idf)*h(jdf)/hTh) / hTh
            Enddo
          Enddo
        Endif   ! Exact calculation of Hessian or Update
      Endif     ! If not converged
    Enddo

    !======
    !RETURN
    !======
    Return
  End Subroutine baker_baker


  !-------------------------------------------------------------
  ! Display the parameters for the baker's algorithm
  !-------------------------------------------------------------
  Subroutine baker_display(bakerparams, nspc, optunitno)
    Type(BAKER_Params), Intent(in) :: bakerparams
    Integer, Intent(in)            :: nspc
    Integer, Intent(in), Optional  :: optunitno
    
    Character(len=nspc)  :: spc
    Integer              :: unitno
    
    spc = " "
    If (Present(optunitno)) Then
      unitno = optunitno
    Else
      unitno = 6
    End If
    
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(2a)') spc, "Parameters for Baker's Algorithm"
    Write(unitno, '(2a)') spc, dashedline
    Write(unitno, '(a,a,i5)') spc, "Search Option (1=TS, 0=Min): ", &
        bakerparams%isaddle
    Write(unitno, '(a,a,i8)') spc, "Maximum Steps              : ", &
        bakerparams%maxstep
    Write(unitno, '(a,a,i5)') spc, "Degrees of Freedom         : ", &
        bakerparams%dof
    Write(unitno, '(a,a,e9.3)') spc, "Step Size Tolerance        : ", &
        bakerparams%toldf
    Write(unitno, '(a,a,e9.3)') spc, "Norm of Gradient Tolerance : ", &
        bakerparams%tolgrad
    Write(unitno, '(a,a,e9.3)') spc, "Max. Step Magnitude        : ", &
        bakerparams%tolgrad
    Write(unitno, '(a,a,i5)') spc, "Steps Between Hessian Calc.: ", &
        bakerparams%nchess
  End Subroutine baker_display

End Module baker
