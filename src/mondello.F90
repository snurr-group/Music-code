!---------------------------------------------------------------------------
! Model for calculating sp3-inversion interaction at a tertiary carbon atom
! Reference: M. Mondello and G.S. Grest JCP vol. 103 p. 7156  (year?!?)
! (see notes for understanding the formulas used here)
! Definition of atom numbers:
!   Atom 1,2,3 refers to the three atoms attached to the tertiary Carbon
!              in anti-clockwise order
!   Atom 4 refers to the tertiary carbon
!   Atom 1 is the one that we want to prevent from sp3-inversion
!
! For example in case of 1-3 dimehtyl CycloPropane:
!   atom 1 is the CH3        1
!   atom 2 is the CH2        |
!   atom 3 is the CH         4 
!   atom 4 is the CH        / \
!     bonded to atom1      /   \
!                         2 --- 3
!                              /
!                             5
!---------------------------------------------------------------------------

Module mondello

  Use defaults, Only: strLen, RDbl, lstrLen, one, zero, degTorad, pi
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use utils, Only: toreal, deallocerrdisplay

  Implicit None
  Save

  Private
  Public :: MondelloModel, mondello_idstring, mondello_init, mondello_display, &
      mondello_getinteraction

  Character(len=strLen), Parameter  :: mondello_idstring = 'MONDELLO'

  Type MondelloModel
    Real(kind=RDbl)          :: k,theta_0
  End Type MondelloModel

Contains
  
  !-----------------------------------------------------------------------
  ! Initializes the Mondello parameters.  
  ! Requires:  mparams -- Mondello data structure
  !            params -- array of 2 strings containing parameters
  !            acoords -- coordinates for the four atoms, might be needed   
  ! Parameters are:  k, equilibrium angle.
  ! CALC option for equilibrium angle is supported
  !-----------------------------------------------------------------------
  Subroutine mondello_init(mparams,params,acoords)
    Type(MondelloModel), Intent(Out)        :: mparams
    Character(*), Dimension(:), Intent(In)  :: params    
    Type(VecType), Dimension(4), Intent(In) :: acoords

    mparams%k = toreal(params(1))

    !** Set the equilibrium angle, calculate if necessary
    If (Index(params(2),'CALC') /= 0) Then
      mparams%theta_0 = mondello_getangle(acoords)
    Else
      mparams%theta_0 = (toreal(params(2)))*degTorad
    End If

!    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
!    Write(*,*) mparams%k,mparams%theta_0,params(1),params(2)

  End Subroutine Mondello_init
    
  !-----------------------------------------------------------------------
  ! r(1:4) contains the 4 position of atom centers comprising the 
  ! torsion list.  It returns u (potential energy) and vf(force vectors)
  ! Requires:  mparams -- Mondello data structure
  !            r -- position vectors for the four atoms
  !            u -- returned potential energy
  !            vf -- returned force vectors
  !-----------------------------------------------------------------------
  Subroutine mondello_getinteraction(params,r,ifvec,u,vf)
    Type(MondelloModel), Pointer                       :: params
    Type(VecType), Intent(In), Dimension(:)            :: r
    Logical, Intent(In) :: ifvec
    Real(kind=RDbl), Intent(Out)                       :: u       
    Type(VecType), Intent(Out), Dimension(:), Optional :: vf

    Real(kind=RDbl)            :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)            :: T1,T2,T3,T1T2,T1T3,T2T3_invroot
    Real(kind=RDbl)            :: cos_theta,sin_theta,ffac,theta
    Type(VecType)              :: a,b,c
    Type(VecType)              :: dcos_da,dcos_db,dcos_dc
    Type(VecType),Dimension(4) :: dcos_dr

    !** Determine whether coordinates or separation vectors have been passed
    If (ifvec) Then
      !** The coordinates passed are actually the separation vectors
      a = r(1)
      b = r(2)
      c = r(3)
    Else
      !** The cooridnates have been passed, calculate the sep vecs
     a=r(2)-r(1)
     b=r(3)-r(2)
     c=r(4)-r(3)
    End If    

    !**calculate the dot products 
    a_a=a*a
    a_b=a*b
    a_c=a*c
    b_b=b*b
    b_c=b*c
    c_c=c*c

    !**calculate T1,T2,T3, term that will have to used later
    T1= a_b*b_c - a_c*b_b
    T2= a_a*b_b - a_b*a_b
    T3= b_b*c_c - b_c*b_c
        
    T2T3_invroot = one/Sqrt(T2*T3)
    T1T2=T1/T2
    T1T3=T1/T3
        
    cos_theta=(T1*T2T3_invroot)
    If (cos_theta>0.99999999) Then
      theta=zero
    Else
      theta= dacos(cos_theta)
    End If

    !**Calculate potential
    u = params%k*(theta-params%theta_0)*(theta-params%theta_0)/2.0_RDbl

    If (Present(vf)) Then
      !** dcos /da, dcos /db, dcos /dc 
      dcos_da = (  (b*( b_c+a_b*T1T2 )) - (a*(T1T2*b_b)) - (c*b_b) ) * &
          T2T3_invroot
      dcos_db=(  (a*(b_c+a_b*T1T2)) + (c*(a_b+T1T3*b_c)) - &
          (b*(2*a_c+T1T2*a_a+T1T3*c_c))   ) * T2T3_invroot
      dcos_dc=( (b*(a_b+b_c*T1T3)) - (a*b_b) - (c*(T1T3*b_b)) ) * T2T3_invroot

      !** The derivative of cos with respect to position vector of 
      !** particles 1..4 : dcos /dr1, dcos /dr2 , dcos /dr3 and dcos /dr4
      dcos_dr(1)=dcos_da*(-one)
      dcos_dr(2)=dcos_da-dcos_db
      dcos_dr(3)=dcos_db-dcos_dc
      dcos_dr(4)=dcos_dc

      !** the pre-factor for forces = k(theta-theta_0)/sin_theta
      If (cos_theta==one.Or.theta==zero) Then
        ffac=zero
      Else
        sin_theta=Sqrt(1.0_RDbl-cos_theta*cos_theta)
        ffac=params%k * ( theta - params%theta_0 ) / sin_theta
      End If

      !** vector forces on particles 1..4
      vf(1)=dcos_dr(1)*(ffac)
      vf(2)=dcos_dr(2)*(ffac)
      vf(3)=dcos_dr(3)*(ffac)
      vf(4)=dcos_dr(4)*(ffac)
    End If

  End Subroutine mondello_getinteraction

  !-----------------------------------------------------------------------
  ! Just calculates the important angle.  Used for initialization only.
  ! Requires:  r -- position vectors for the four atoms
  !-----------------------------------------------------------------------
  Real(kind=RDbl) Function mondello_getangle(r)
    Type(VecType), Dimension(4), Intent(In)  :: r

    Real(kind=RDbl)            :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)            :: T1,T2,T3,T2T3_invroot
    Real(kind=RDbl)            :: cos_theta
    Type(VecType)              :: a,b,c

    !** Calculate a b c, these represent the vectors along the bonds
    a = r(2)-r(1)
    b = r(3)-r(2)
    c = r(4)-r(3)
    
    !** Calculate the dot products 
    a_a = a*a
    a_b = a*b
    a_c = a*c
    b_b = b*b
    b_c = b*c
    c_c = c*c

    !** Calculate T1,T2,T3, term that will have to used later
    T1 = a_b*b_c - a_c*b_b
    T2 = a_a*b_b - a_b*a_b
    T3 = b_b*c_c - b_c*b_c
        
    T2T3_invroot = one/Sqrt(T2*T3)
        
    cos_theta = (T1*T2T3_invroot)
    If (cos_theta > 0.99999999) Then
      mondello_getangle = zero
    Else
      mondello_getangle = dacos(cos_theta)
    End If

  End Function mondello_getangle

  !------------------------------------------------------------------
  ! Displays the contents of the mparams
  ! Requires:  mparams -- Mondello data structure
  !------------------------------------------------------------------
  Character(lstrLen) Function mondello_display(mparams)
    Type(MondelloModel) :: mparams
    
    Write(mondello_display,'(1x,a,2f10.4)')  "k, theta_0", &
        mparams%k, mparams%theta_0

  End Function Mondello_display
  
  !-----------------------------------------------------
  ! Cleanup the Mondello parameters
  ! Requires:  mparams -- Mondello data structure
  !-----------------------------------------------------
  Subroutine mondello_cleanup(mparams)
    Type(MondelloModel), Pointer :: mparams

    Integer :: error
    
    Deallocate(mparams,STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
    
  End Subroutine mondello_cleanup
  
End Module mondello


