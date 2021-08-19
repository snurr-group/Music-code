!----------------------------------------------------------------------------
! contains subroutines used for calculating binary adsn isotherm from ias
!-----------------------------------------------------------------------------
Module iasbinary

  Use defaults, Only: strLen, lstrLen, RDbl
  Use isotherm, Only: Isotherm_Info

  Implicit None
  Save

  Private
  Public :: iasbinary_interpolate , iasbinary_lookup, iasbinary_solveGivenSigP

Contains




  !---------------------------------------------------------------------
  ! S U B R O U T I N E S
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! given tsingle component isotherms and a particular surface 
  ! pressure(sigma) and total pressure(P) calculates the molefractions 
  ! and loading using IAS theory (Myers, AIChe 1965). Note the order of 
  ! isotherm objects
  !---------------------------------------------------------------------
  Subroutine iasbinary_SolveGivenSigP(sigma, P, topIso, botIso, x1, x2, y1, &
      y2, n_tot)
    Real(kind=RDbl), Intent(in)  :: sigma,P
    Real(kind=RDbl), Intent(out) :: x1, x2, y1, y2, n_tot
    Type(Isotherm_Info), Intent(in) :: topIso, botIso

    Real(kind=RDbl) :: p10, p20, pi1, pi2, n10, n20
    Integer :: length
    Logical :: outOfRange
    ! get the individual single component pressure at the sigma
    p10=iasbinary_interpolate(sigma, topIso%sigma, topiso%P,outOfRange)
    n10=iasbinary_interpolate(p10, topiso%P,topIso%N, outOfRange)
    If (outOfRange) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif
    p20=iasbinary_interpolate(sigma, botIso%sigma, botiso%P,outOfRange)
    If (outOfRange) Then
      ! this could be because of flat sigma vs P curve
      ! so let extrapolate
      length=Size(botIso%sigma,1)
      p20=botIso%P(length) * &
          Exp(( sigma- botIso%sigma(length) ) / botIso%maxUCLoading )
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "extrpolated p20, maybe because of flat curve ", p20
      n20=botIso%maxUCLoading
    Else
      n20=iasbinary_interpolate(p20, botiso%P,botIso%N, outOfRange)
    Endif

    x1=(p20-P)/(p20-p10)
    x2=(P-p10)/(p20-p10)
    pi1=(P/p10 )*sigma     ! point F in fig.1
    pi2=(P/p20 )*sigma     ! point E in fig.1
    y1=(sigma-pi2)/(pi1-pi2)
    y2=(pi1-sigma)/(pi1-pi2)
    If ((x1<-0.00001).Or.(x1>1.00001).Or.(x2<-0.00001).Or.&
        (x2>1.00001).Or.(Abs(x1+x2-1.00000000)>1.00e-5)) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif


    If ((y1<-0.00001).Or.(y1>1.00001).Or.(y2<-0.00001).Or.&
        (y2>1.00001).Or.(Abs(y1+y2-1.00000000)>1.00e-5)) Then
      Write(*,*) y1, y2
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      stop
    Endif
    n_tot= 1/(  x1/n10 + x2/n20 )

  End Subroutine iasbinary_SolveGivenSigP

  !---------------------------------------------------------
  ! given two ordered arrays xarr and yarr calculates the iasbinary_interpolated value
  ! of y at given x
  !---------------------------------------------------------
  Real(kind=RDbl) Function iasbinary_interpolate (x , xarr, yarr, outOfRange)
    Real(kind=RDbl), Intent(in) :: x
    Real(kind=RDbl), Dimension(:), Intent(in) :: xarr, yarr
    Logical, Intent(out)::outOfRange
    Real(kind=RDbl) :: y, dy_dx, deltax
    Integer :: high, low
    outOfRange=.True.
    If ((x<xarr(1)).Or.(x>xarr(Size(xarr)))) Then
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write (*,*) "cant inter polate; out of range "
      Return
    Endif

    high=iasbinary_lookup(x,xarr)  ! the interval
    low=high-1           !
    dy_dx=(yarr(high)-yarr(low))/(xarr(high)-xarr(low))
    deltax=x-xarr(low)
    y=yarr(low)+deltax* dy_dx
    iasbinary_interpolate=y
    outOfRange=.False.
  End Function iasbinary_interpolate



  !-----------------------------------------------
  ! given a value of x finds the index of the smallest xarr member 
  ! with value higher  than x, assumes xarr is sorted
  ! Heavily hacked dont reuse anywhere else copy from; bmap code instead
  !-----------------------------------------------
  Integer Function iasbinary_lookup(x, xarr)
    Real(kind=RDbl), Dimension(:), Intent(in) :: xarr
    Real(kind=RDbl), Intent(in) :: x

    Integer   :: low, high, mid

    Real(kind=RDbl) :: mid_x
    low = 1
    high = size(xarr)

    !** Do a binary search to find the first index i so that xarr(i) > x
    Do
!!$$$$      If (low >= high) Exit ! Wrong ?
      If (low > high) Exit
      mid = (low + high)/2
      mid_x = xarr(mid) 

      If (x < mid_x) Then
        high = mid - 1
      Else If (x > mid_x ) Then
        low = mid + 1
        mid = low
      Else If (x == mid_x) Then
        mid=mid+1
        Exit
      Endif
    Enddo

    !** low is equal to high
    iasbinary_lookup=mid
    Return

  End Function iasbinary_lookup
End Module iasbinary
