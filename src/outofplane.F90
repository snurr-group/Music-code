!-----------------------------------------------------------------------
! Module for calculating out-of-plane potentials
! (see Shaji's notes for understanding the formulas used here)
!   Potential energy = k*D*D
!   Force = -2k*D *  (d_D/d_r)
!  where D is the difference from first atom
!  to the plane formed by the other 3 atoms
!-----------------------------------------------------------------------

Module outofplane

  Use defaults, Only: strLen, RDbl, lstrLen, one, zero
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag, vector_crossprod
  Use utils, Only: toreal, deallocerrdisplay

  Implicit None
  Save

  Private
  Public :: outofplaneModel, outofplane_idstring, outofplane_init, &
      outofplane_display, outofplane_getinteraction, outofplane_cleanup

  Character(len=strLen), Parameter  :: outofplane_idstring = 'OUTOFPLANE'

  !** Potential is kD^2 or 1/2 k D^2, coded here as kD^2
  Type OutofplaneModel
    Real(kind=RDbl)   ::  k  
  End Type OutofplaneModel

Contains
  
  !-----------------------------------------------------------------
  ! Initializes the Out-of-Plane parameters.  
  ! Requires:  oopParams -- Out-of-Plane data structure
  !            params -- array of strings containing parameters
  ! Parameters are:  k
  !-----------------------------------------------------------------
  Subroutine outofplane_init(oopParams,params)
    Type(OutofplaneModel) :: oopParams
    Character(*), Dimension(:) :: params

    oopParams%k = toreal(params(1))
    
  End Subroutine outofplane_init
    
  !-----------------------------------------------------------------
  ! Get the Out-of-plane interactions
  ! Requires:  oopParams -- Out-of-Plane data structure
  !            r -- position vectors for the four atoms
  !            u -- returned potential energy
  !            vf -- returned force vectors
  !-----------------------------------------------------------------
  Subroutine outofplane_getinteraction(params,r,ifvec,u,vf)
    Type(OutofplaneModel), Intent(In)                  :: params
    Type(VecType), Intent(IN), Dimension(:) :: r ! positions of atoms, or vectors
    Logical, Intent(In) :: ifvec
     Real(kind=RDbl), Intent(Out)                       :: u
    Type(VecType), Intent(Out), Dimension(:), Optional :: vf

    Integer                :: i
    Real(kind=RDbl)        :: a_a,a_b,b_b
    Type(VecType)          :: aCROSb,bCROSc,cCROSa,a,b,c
    Real(kind=RDbl)        :: T1,T2,D,T1_T2,A2,A1
    
    !** Determine whether coordinates or separation vectors have been passed
    If (ifvec) Then
      !** The coordinates passed are actually the separation vectors
      a = r(2)
      b = r(2)+r(3)
      c%comp(1) = -1.0*r(1)%comp(1)
      c%comp(2) = -1.0*r(1)%comp(2)
      c%comp(3) = -1.0*r(1)%comp(3)
      !** c=-r(1)
       Else
      !** The cooridnates have been passed, calculate the sep vecs
    a=r(3)-r(2)
    b=r(4)-r(2)
    c=r(1)-r(2)
    End If   

    !**calculate the required  dot products  between the vectors

    a_a=a*a
    a_b=a*b
    b_b=b*b
    
    !** Calculate the required cross products
    aCROSb = vector_crossprod(a,b)
    bCROSc = vector_crossprod(b,c)
    cCROSa = vector_crossprod(c,a)

    !** Calculate all those prefactors here
    T1 = aCROSb * c       ! a dot product
    T2 = aCROSb * aCROSb
    D = T1/sqrt(T2)   ! D = distance we are interested in
                      ! But we dont need to keep trak of it, so can be 
                      ! removed after debugging

!    Write(*,*) '------OUTOFPLANE-----------'   
!    Write(*,*) "Distance ",  D
!    Write(*,*) "Energy   ", params%k*D*D
    T1_T2=T1/T2

    !** energy
    A2 = 2 * params%k *  T1_T2
    A1 = A2 * T1_T2 

    !** Lets get into it
    u=  (A2 * T1)/2
!    Write(*,*) "U",u

    !** Calculate forces on atoms if force vectors are present
    If (Present(vf)) Then
      vf(1) = aCROSb *(-A2)
      vf(3) =  a *(A1*b_b) - ( b*(A1*a_b) + bCROSc*(A2) )
      vf(4) =  b *(A1*a_a) - ( a*(A1*a_b) + cCROSa*(A2) )
      vf(2) =  ( (vf(4) + vf(3)) + vf(1) )*(-one)
    End If

  End Subroutine outofplane_getinteraction

  !---------------------------------------------------------------
  ! Displays the contents of the Out-of-Plane data structure
  ! Requires:  oopParams -- Out-of-Plane data structure
  !---------------------------------------------------------------
  Character(lstrLen) Function outofplane_display(oopParams)
    Type(OutofplaneModel) :: oopParams
    
    Write(outofplane_display,'(1x,f10.4)') oopParams%k
    
  End Function outofplane_display
  
  !-----------------------------------------------------
  ! Cleanup the Out-of-Plane data structure
  ! Requires:  oopParams -- Out-of-Plane data structure
  !-----------------------------------------------------
  Subroutine outofplane_cleanup(oopParams)
    Type(OutofplaneModel), Pointer :: oopParams

    Integer :: error
    
    Deallocate(oopParams,STAT=error)
    If (error/=0) Call deallocErrDisplay(__FILE__,__LINE__)    
    
  End Subroutine outofplane_cleanup
  
End Module outofplane



















