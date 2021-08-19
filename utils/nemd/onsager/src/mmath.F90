!-------------------------------------------------------------------------
! This module contains general math subroutines and functions.
! written in music style
! examples : least square fiting
!-------------------------------------------------------------------------

Module mmath

  Use defaults, Only: strLen, lstrLen, RDbl, COMMENT_CHAR, twopi, xlstrlen, &
      one, zero, dbgflag

  Implicit None
  Save

  Private
  Public :: mmath_lsFit, mmath_simpson

Contains

  !----------------------------------------------------------
  ! Curve Fitting subroutines
  !----------------------------------------------------------

  !------------------------------------------
  ! Gives the slope and intercept,given x and y arrays,and number of elements
  ! see notes on this in documentation
  ! percentage error is returned in err
  ! fits to line y=mx+c
  !------------------------------------------
  Subroutine mmath_LSFit(x,y,N,m,c,Err)
    Real(kind=RDbl), Dimension(:),Intent(in)::x,y
    Real(kind=RDbl),Intent(out)::m,c,Err
    Integer,Intent(in):: N
    
    Integer::i
    Real(kind=RDbl)::xsum,ysum,xysum,x2sum,y2sum
    xsum=0.0_RDbl
    ysum=0.0_RDbl
    xysum=0.0_RDbl
    x2sum=0.0_RDbl
    y2sum=0.0_RDbl

    Do i=1,N
      xsum=xsum+x(i)
      ysum=ysum+y(i)
      xysum=xysum+x(i)*y(i)
      x2sum=x2sum+x(i)*x(i)
      y2sum=y2sum+y(i)*y(i)
    End Do

    If (ysum == 0.0_RDbl) Then
      m = 0.0_Rdbl
      c = 0.0_RDbl
      Err = 0.0_RDbl
      Return
    End If

    m=(xsum*ysum-N*xysum)/(xsum*xsum-N*x2sum)
    c=(xysum*xsum-x2sum*ysum)/(xsum*xsum-N*x2sum)
    Err=(y2sum+m*m*x2sum+N*c*c-2*m*xysum-2*c*ysum+2*m*c*xsum)/N


    !** To make sure error is +ve
    !** , truncation errors can make it negative
    Err= abs(Err)

    Err=100*Sqrt(Err/y2sum)/N

  End Subroutine  mmath_lsFit  

  !----------------------------------------------------------------------
  ! does simpson integration using quadratic interpolation
  !  -> calculates \integ y dx
  !  -> h =dx, n =size(y)
  !----------------------------------------------------------------------
  Real(kind=RDbl) Function mmath_simpson(y,h,n)
    Real(kind=RDbl),Dimension(:) , Intent(in) :: y
    Real(kind=RDbl), Intent(in) :: h
    Integer, Intent(in) :: n

    Real(kind=RDbl) :: odd_sum, even_sum, sum, integ
    Integer :: n_even, n_odd, i
    Logical :: even
    If (n<3) Then
      Write(*,*) "simpson integration needs minimum 3 y values"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
    !    figure out we have odd or even number of time-points
    even=.false.
    If (Mod(n,2)==0) even=.True.

    ! points are 1,2,3.................n
    ! odd points are 1,3,......
    ! evenpoitns are 2,4,.......

    If (even) Then
      n_odd = n/2
    Else
      n_odd = (n+1)/2
    Endif
    n_even=n_odd-1

    ! find odd sum
    odd_sum=zero
    odd_sum=odd_sum+y(1)+y(n_odd*2-1)
    Do i=2,n_odd-1
      odd_sum=odd_sum+2*y(i*2-1)
    End Do

    ! find even sum
    even_sum=zero
    Do i=1,n_even
      even_sum=even_sum+4*y(i*2)
    End Do

    sum=even_sum+odd_sum

    integ=sum*h/3
    If (even) integ=integ+(y(n)+y(n-1))*(h/2) ! add the last term by 
                                              ! trapezoidal rule 
    mmath_simpson=integ
  End Function mmath_simpson
End Module mmath

