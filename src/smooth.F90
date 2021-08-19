!------------------------------------------------------------------------------
! This module handles the smoothing of functions.  Currently it is used 
! to smooth the pairwise coulombic interactions that are evaluated using
! a cut-off.  It has facilities for using two different methods:
! 
! Function Multiplication
!   Uses a cubic function to smooth the given potential.  The function is
!   such that it is one at r1 and zero at rc.  It's derivatives are also
!   constrained to be zero at these points.  The result is that
!   multiplication by this function effectively smooths the potential in
!   the same way that a spline would, but without knowing the function
!   value or its derivative, f(r) = a + br + cr^2 + dr^3
!
! Spline
!   This method uses a cubic spline to smooth the potential and force to
!   zero at the cutoff while maintaining continuity of the potential and
!   force at the point when the spline kicks in, u(r) = a + br + cr^2 + dr^3.
!   This does not require the current value of the function at r because it 
!   is simply a fit.  However, it must have the value of the function and 
!   its derivative at the smoothing radius.
!------------------------------------------------------------------------------

Module smooth

  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only: real2str, ToUpper, split, allocerrdisplay, deallocerrdisplay
  Use file, Only: file_open
  Use store, Only: Store_Level_Pair,store_scaleall
  Use storebase, Only: EnergyPlus,storebase_scalarmult

  Implicit None
  Save

  Private
  Public :: Smooth_Method,CubicSmooth,CubicSpline,smooth_init,smooth_perform, &
      smooth_perform1,smooth_do,smooth_check,smooth_display,smooth_clean

  Type Smooth_Method
    Type(CubicSmooth), Pointer   :: fnmult
    Type(CubicSpline), Pointer   :: spline
  End Type Smooth_Method

  Type CubicSmooth
    Real(kind=RDbl)        :: smoothrad,cutrad,delta
    Real(kind=RDbl)        :: a,b,c,d
  End Type CubicSmooth

  Type CubicSpline
    Real(kind=RDbl)        :: smoothrad,cutrad,delta
    Real(kind=RDbl)        :: pot,deriv
  End Type CubicSpline

Contains
  !----------------------------------------------------------------------------
  ! Initialize the smoothing method.  The potential and derivation at the
  ! smoothing radius are only necessary if the spline method will be used.
  ! Requires:  params -- Smoothing parameters to initialize
  !            line -- string containing initialization information
  !            smoothrad -- radius to begin smoothing
  !            cutrad -- radius at which function should be zero
  !            pot -- function value at smooth radius
  !            deriv -- function derivative at smooth radius
  !----------------------------------------------------------------------------
  Subroutine smooth_init(params,line,smoothrad,cutrad,pot,deriv)
    Type(Smooth_Method), Intent(InOut)        :: params
    Character(*), Intent(In)                  :: line
    Real(kind=RDbl), Intent(In)               :: smoothrad,cutrad
    Real(kind=RDbl), Intent(In), Optional     :: pot,deriv

    Integer                                   :: error,nfields
    Character(len=strLen), Dimension(20)      :: fields

    !** Nullify the pointer set
    Nullify(params%fnmult)
    Nullify(params%spline)

    !** Select the smoothing method to initialize
    nfields = split(line,fields)
    Select Case(ToUpper(fields(1)))
    Case ('FNMULT')
      Allocate(params%fnmult,stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)   

      Call smooth_initfn(params%fnmult,smoothrad,cutrad)

    Case ('SPLINE')
      Allocate(params%spline,stat=error)    
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)   
      params%spline%smoothrad = smoothrad
      params%spline%cutrad = cutrad
      params%spline%delta = cutrad - smoothrad
      params%spline%pot = pot
      params%spline%deriv = deriv

      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' WARNING: Spline smoothing has not been checked'

    Case Default
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          ' Unable to identify smoothing specification string ',Trim(fields(1))
      Stop
    End Select

    !** Perform a couple single checks
!    Call smooth_check(params,'smooth_check.txt')

  End Subroutine smooth_init

  !----------------------------------------------------------------------------
  ! Initialize the function multiplication smoothing method
  ! Requires:  params -- Smoothing parameters to initialize
  !            smoothrad -- radius to begin smoothing
  !            cutrad -- radius at which function should be zero
  !----------------------------------------------------------------------------
  Subroutine smooth_initfn(params,smoothrad,cutrad)
    Type(CubicSmooth), Intent(InOut)       :: params
    Real(kind=RDbl), Intent(In)            :: smoothrad,cutrad

    Real(kind=RDbl)             :: delta,delta3,delta3inv

    params%smoothrad = smoothrad
    params%cutrad = cutrad
    params%delta = cutrad - smoothrad

    delta = smoothrad - cutrad
    delta3 = delta*delta*delta
    delta3inv = 1.0e0/delta3
    
    params%a = (3.0e0 * smoothrad - cutrad)*(cutrad**2)*delta3inv
    params%b = -6.0e0*smoothrad*cutrad*delta3inv
    params%c = 3.0e0*(smoothrad + cutrad)*delta3inv
    params%d = -2.0e0*delta3inv

  End Subroutine smooth_initfn

  !----------------------------------------------------------------------------
  ! Perform the smoothing method on a store level pair structure
  ! Requires:  params -- Smoothing parameters to initialize
  !            radius -- radius at which to do evaluation
  !            lp -- the Level Pair structure
  !----------------------------------------------------------------------------
  Subroutine smooth_perform(params,radius,lp)
    Type(Smooth_Method), Intent(In)        :: params
    Real(kind=RDbl), Intent(In)            :: radius
    Type(Store_Level_Pair), Intent(InOut)  :: lp

    Real(kind=RDbl)            :: factor

    !** Quick check
    If (.Not. Associated(params%fnmult)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Not setup to use any method other than function multiplication'
      Stop
    End If

    !** get the multiplication factor
    factor = smooth_do(params,radius,1.0_RDbl)

    !** Scale the interactions
    Call store_scaleall(lp,factor)

  End Subroutine smooth_perform

  !----------------------------------------------------------------------------
  ! Perform the smoothing method on just an EnergyPlus structure
  ! Requires:  params -- Smoothing parameters to initialize
  !            radius -- radius at which to do evaluation
  !            results -- EnergyPlus structure
  !----------------------------------------------------------------------------
  Subroutine smooth_perform1(params,radius,results)
    Type(Smooth_Method), Intent(In)   :: params
    Real(kind=RDbl), Intent(In)       :: radius
    Type(EnergyPlus), Intent(InOut)   :: results

    Real(kind=RDbl)            :: factor

    !** Quick check
    If (.Not. Associated(params%fnmult)) Then
      Write(0,'(2a,i4,a)') __FILE__,": ",__LINE__, &
          ' Not setup to use any method other than function multiplication'
      Stop
    End If

    !** get the multiplication factor
    factor = smooth_do(params,radius,1.0_RDbl)

    !** Scale the interactions
    Call storebase_scalarmult(results,factor)

  End Subroutine smooth_perform1

  !----------------------------------------------------------------------------
  ! Perform the smoothing method.  The function value need only be passed
  ! if the function multiplication method is initialized.  Returns the new 
  ! value of the function.
  ! Requires:  params -- Smoothing parameters to initialize
  !            radius -- radius at which to do evaluation
  !            invalue -- function value at radius
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function smooth_do(params,radius,invalue)
    Type(Smooth_Method), Intent(In)        :: params
    Real(kind=RDbl), Intent(In)            :: radius
    Real(kind=RDbl), Intent(In), Optional  :: invalue

    If (Associated(params%fnmult)) Then
      smooth_do = smooth_fnmult(params%fnmult,radius,invalue)

    Else If (Associated(params%spline)) Then
      smooth_do = smooth_spline(params%spline,radius)

    Else
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' No recognized pointers associated for smoothing method'
      Stop
    End If

  End Function smooth_do

  !----------------------------------------------------------------------------
  ! Perform the cubic function multiplication method
  ! Requires:  params -- Smoothing parameters to initialize
  !            r -- radius at which to do evaluation
  !            invalue -- function value at radius
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function smooth_fnmult(params,r,invalue)
    Type(CubicSmooth), Intent(In)  :: params
    Real(kind=RDbl), Intent(In)    :: r,invalue

    Real(kind=RDbl)        :: r2

    r2 = r*r
    smooth_fnmult = invalue*(params%a + params%b*r + params%c*r2 + params%d*r*r2)

  End Function smooth_fnmult

  !----------------------------------------------------------------------------
  ! Perform the spline fit method
  ! Requires:  params -- Smoothing parameters to initialize
  !            r -- radius at which to do evaluation
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function smooth_spline(params,r)
    Type(CubicSpline), Intent(In)  :: params
    Real(kind=RDbl), Intent(In)    :: r

    Real(kind=RDbl)             :: r2,rc2,rc,r1,ur1,u1r1
    Real(kind=RDbl)             :: delta,delta3,delta3inv
    Real(kind=RDbl)             :: a,b,c,d
    Real(kind=RDbl)             :: r1ur1,rcur1,rcu1r1,r1u1r1
    Real(kind=RDbl)             :: r1rcu1r1,rc2u1r1,r12u1r1

    rc = params%cutrad
    r1 = params%smoothrad
    ur1 = params%pot
    u1r1 = params%deriv

    delta = r1 - rc
    delta3 = delta*delta*delta
    delta3inv = 1.0e0/delta3
    
    r2 = r*r
    rc2 = rc*rc
    r1ur1 = r1*ur1
    rcur1 = rc*ur1
    rcu1r1 = rc*u1r1
    rc2u1r1 = rcu1r1*rc
    r1u1r1 = r1*u1r1
    r12u1r1 = r1u1r1*r1
    r1rcu1r1 = r1u1r1*rc

    !** Calculate the spline coefficients
    a = rc2 * (-r12u1r1 + r1rcu1r1 + 3*r1ur1 - rcur1) * delta3inv
    b = -rc * (-2.0e0*r12u1r1 + r1rcu1r1 + rc2u1r1 + 6.0e0*r1ur1) * delta3inv
    c = (-r12u1r1 - r1rcu1r1 + 2.0e0*rc2u1r1 + 3.0e0*r1ur1 + &
        3.0e0*rcur1)*delta3inv
    d = (r1u1r1 - rcu1r1 - 2.0e0*ur1) * delta3inv

    !** Calculate the value of the function at the radius
    smooth_spline = a + b*r + c*r2 + d*r*r2

  End Function smooth_spline

  !----------------------------------------------------------------------------
  ! Check the smoothing method.  Will dump a set of data points to a file
  ! if a name is supplied.
  ! Requires:  params -- Smoothing method pointer set
  !            filename -- optional filename
  !----------------------------------------------------------------------------
  Subroutine smooth_check(params,filename)
    Type(Smooth_Method), Intent(In)     :: params
    Character(*), Intent(In), Optional  :: filename

    Integer                     :: i,npts,unit
    Real(kind=RDbl)             :: r,r0,fn

    !** Get function value at the cut-off radius
    If (Associated(params%fnmult)) Then
      fn = smooth_do(params,params%fnmult%cutrad,1.0_Rdbl) 
    Else If (Associated(params%spline)) Then
      fn = smooth_do(params,params%spline%cutrad) 
    Else
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' No recognized pointers associated for smoothing method'
      Stop
    End If

    !** Make sure the result is zero at the end point
    If (Abs(fn) > 1.0e-5_RDbl) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' Smoothing function does not satisfy zero at r_cut criteria'
      Stop
    End If

    !** Make sure the first derivative of result is zero at the end point
    If (Associated(params%fnmult)) Then
      fn = (params%fnmult%b + 2*params%fnmult%c*params%fnmult%cutrad + &
          3*params%fnmult%d*(params%fnmult%cutrad**2))
    End If
    If (Abs(fn) > 1.0e-5_RDbl) Then
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' Smoothing function does not have zero derivative at r_cut'
      Stop
    End If

    !** Generate a datafile with function values between r_smooth and r_cut
    npts = 100
    If (Present(filename)) Then
      unit = file_open(filename)

      If (Associated(params%fnmult)) Then
        Write(unit,'(a)') '# blah blah'
        r0 = params%fnmult%smoothrad
        Do i = 1,npts
          r = r0 + (i-1)*(params%fnmult%delta/npts)
          Write(unit,'(f8.3,e14.6)') r,smooth_do(params,r,1.0_Rdbl) 
        End Do
  
      Else If (Associated(params%spline)) Then
        Write(unit,'(a)') '# blah blah'
        r0 = params%spline%smoothrad
        Do i = 1,npts
          r = r0 + (i-1)*(params%spline%delta/npts)
          Write(unit,'(f8.3,e14.6)') r,smooth_do(params,r) 
        End Do
  
      Else
        Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
            ' No recognized pointers associated for smoothing method'
        Stop
      End If
    End If

  End Subroutine smooth_check

  !----------------------------------------------------------------------------
  ! Displays the smoothing method parameters
  ! Requires:  params -- Smoothing parameters to initialize
  !            indent -- no. of spaces from the left margin
  !            unitno -- optional display unit number
  !----------------------------------------------------------------------------
  Subroutine smooth_display(params,indent,unit)
    Type(Smooth_Method), Intent(In)  :: params
    Integer, Intent(In)              :: indent
    Integer, Intent(In)              :: unit

    Character(len=indent)      :: blank

    blank = Repeat(' ',indent)

    If (Associated(params%fnmult)) Then
      Write(unit,'(2a)') blank,'Using cubic function multiplication method'
      Write(unit,'(2a,f8.3)') blank,'smoothing radius = ',params%fnmult%smoothrad
    
    Else If (Associated(params%spline)) Then
      Write(unit,'(2a)') blank,'Using cubic spline method'
      Write(unit,'(2a,f8.3)') blank,'smoothing radius = ',params%spline%smoothrad
    
    Else
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' No recognized pointers associated for smoothing method'
      Stop
    End If

  End Subroutine smooth_display

  !----------------------------------------------------------------------------
  ! Cleans the smoothing method pointer set
  ! Requires:  params -- Smoothing parameters to initialize
  !----------------------------------------------------------------------------
  Subroutine smooth_clean(params)
    Type(Smooth_Method), Intent(InOut)        :: params

    Integer          :: error

    If (Associated(params%fnmult)) Then
      Deallocate(params%fnmult, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)    

    Else If (Associated(params%spline)) Then
      Deallocate(params%spline, STAT=error)
      If (error /= 0) Call deallocErrDisplay(__FILE__,__LINE__)        

    Else
      Write(0,'(2a,i4,a)') __FILE__,":",__LINE__, &
          ' No recognized pointers associated for smoothing method'
      Stop
    End If

  End Subroutine smooth_clean

End Module smooth
