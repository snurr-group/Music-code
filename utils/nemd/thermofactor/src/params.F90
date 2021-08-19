!-----------------------------------------------------------------------------
! deals with an isotherm as a parameter set of Jensen and Seaton
! Ref: Heuchel,  snurr and Buss Langmuir , vol.13, No.25, 1997, 
! pages: 6795-6804, Eqn-12
! Eqn is :   N_p=(K*p) / ( one+ (((K/Lmax)*p/(1+b*P))**c) )**(one/c)
!
!-------------------------------------------------------------------
Module params

  Use defaults, Only: strLen, lstrLen, RDbl, zero, one
  Use file, Only : file_open
  Use random, only : rranf
  Use utils, Only   : stripcmnt, split, allocErrDisplay, toint, toreal

  Implicit None
  Save

  Private
  Public :: ParamsSet, params_init, params_guess, params_solve, &
      params_eval, params_copy, params_interpol

  Type ParamsSet
    Real(kind=RDbl)      :: y
    Real(kind=RDbl)      :: K, Lmax, c, b
  End Type ParamsSet

  Character(len=lstrLen)       :: description 

Contains

  !----------------------------------------------------------------------------
  ! Reads the params and initializes
  ! here y reference to the y of the compound used in binary-gcmc
  !----------------------------------------------------------------------------
  Subroutine params_init(iso, unitno)
    Type(ParamsSet), Intent(inout) :: iso
    Integer, Intent(In)            :: unitno
    Integer                :: error, i, nfields, npts
    Character(len=strlen) :: str
    Read(unitno,*) str, iso%y, iso%K, iso%Lmax, iso%c, iso%b 
    ! convert from % to mole frac
    iso%y=iso%y/100
  End Subroutine params_init

  !----------------------------------------------------------------------------
  ! copies params
  !----------------------------------------------------------------------------
  Subroutine params_copy(ttarget, ssource)
    Type(ParamsSet), Intent(inout) :: ttarget
    Type(ParamsSet), Intent(inout) :: ssource
    ttarget%y=ssource%y 
    ttarget%K=ssource%K 
    ttarget%Lmax=ssource%Lmax 
    ttarget%c=ssource%c 
    ttarget%b=ssource%b 
  End Subroutine params_copy
!!$
!!$  !----------------------------------------------------------------------------
!!$  ! copies params
!!$  !----------------------------------------------------------------------------
!!$  Subroutine params_copy(ttarget, ssource)
!!$    Type(ParamsSet), Intent(inout) :: ttarget
!!$    Type(ParamsSet), Intent(inout) :: ssource
!!$    ttarget%y=ssource%y 
!!$    ttarget%K=ssource%K 
!!$    ttarget%Lmax=ssource%Lmax 
!!$    ttarget%c=ssource%c 
!!$    ttarget%b=ssource%b 
!!$  End Subroutine params_copy
!!$
  !----------------------------------------------------------------------------
  ! interploates params
  !----------------------------------------------------------------------------
  Subroutine params_interpol(iso, p0, p1, y)
    Type(ParamsSet), Intent(inout) :: iso
    Type(ParamsSet), Intent(in) :: p0, p1
    Real(kind=Rdbl),Intent(in) :: y
    Real(kind=Rdbl) :: y0, y1, fac
    y0=p0%y
    y1=p1%y
    fac=(y-y0)/(y1-y0)
    iso%y   = y
    iso%K   = p0%K + fac*(p1%K-p0%K)
    iso%Lmax= p0%Lmax + fac*(p1%Lmax-p0%Lmax)
    iso%c   = p0%c+fac*(p1%c-p0%c)
    iso%b   = p0%b+fac*(p1%b-p0%b)
  End Subroutine params_interpol

  !----------------------------------------------------------------------------
  ! Given P ealuates N, just calls function
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function params_eval(isoparam, P)
    Type(ParamsSet), Intent(in) :: isoparam
    Real(kind=RDbl) , Intent(in) :: P
    params_eval=params_function(isoparam, P)
  End Function params_eval

  !----------------------------------------------------------------------------
  ! Guesses a soln based on approximate solution
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function params_guess(isoparam, N)
    Type(ParamsSet), Intent(in) :: isoparam
    Real(kind=RDbl) , Intent(in) :: N
    Real(kind=RDbl):: P
    ! langmuir
    P=(N/isoparam%K)/(1-N/isoparam%Lmax)

    ! linear
    If (P<zero) Then
      P=(N/isoparam%K)
    Endif

    params_guess=P
  End Function params_guess

  !----------------------------------------------------------------------------
  ! Solves for the correct P using bisection, between P0 and P1
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function params_solve(isoparam, N, P0, P1)
    Type(ParamsSet), Intent(in) :: isoparam
    Real(kind=RDbl) , Intent(in) :: P0, P1, N
    Real(kind=RDbl) , Parameter :: N_toler=0.05, P_toler=1.0
    Integer, Parameter :: MAX_ITER=100
    Real(kind=RDbl) :: lowP, highP, newP, newN, deltaN, deltaP
    Integer :: i, count
    Logical :: found
    count=zero
    found=.false.
    Do i=1,MAX_ITER
      count=count+1
      newP=(lowP+highP)/2
      newN=params_eval(isoparam, newP)
      deltaN=Abs(newN-N)
      deltaP=Abs(newP-lowP)
      If((deltaP<P_toler).and.(deltaN<N_toler)) Then
        found=.True.
        Exit
      Endif
      If (newN<N) Then
        lowP=newP
      Else
        highP=newP
      Endif
    End Do
    If(found) Then
      params_solve=newP
      Write(*,*)"found solution in ", count, " trials"
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif
  End Function  params_solve

  !----------------------------------------------------------------------------
  ! This is the fitting function Toth Eqn, Heuchel , snurr and buss eqn-12
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function params_function(isoparam, p)
    Type(ParamsSet), Intent(in) :: isoparam
    Real(kind=RDbl) , Intent(in) :: p
    Real(kind=RDbl) :: K, Lmax, c,b
    K=isoparam%K
    Lmax=isoparam%Lmax
    c=isoparam%c
    b=isoparam%b
    If (((K/Lmax)*p/(1+b*P))<zero) then
      Write(*,*) K, Lmax, b, c
      Write(*,*) b,P,b*P
 
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    endif
    params_function=(K*p) / ( one+ (((K/Lmax)*p/(1+b*P))**c) )**(one/c)
  End Function params_function


End Module params
