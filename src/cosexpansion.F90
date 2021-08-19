!------------------------------------------------------------------------
! Model for calculating torional interactions with flexible bonds
! see notes for understanding the formulas used here
!
! potential energy = SUM(k=0,5) (V_k * Cos(phi)^k)
!   Where phi is the angle between the two plane normals
!         i+1       i-1
!          o         o       plane1: defined by bonds bi-1, bi
!         / \       /        plane2: defined by bonds bi, bi+1
!   bi+1 /   \ bi  / bi-1     
!       /     \   /          
!      /       \ /           
! i+2 o       i o     
!
! In the 'trans' conformation above, the angle is zero because
! the two plane normals are parallel.  In the 'gauche' conformation,
! it is Pi.
!------------------------------------------------------------------------
Module cosexpansion

  Use defaults, Only: strLen, RDbl, lstrLen, one, zero
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use utils, Only: toreal

  Implicit None
  Save

  Private
  Public :: CosExpansionModel, cosexpansion_idstring, cosexpansion_init, &
      cosexpansion_display, cosexpansion_getinteraction, cosexpansion_cleanup

  Integer, Parameter                :: MAX_PARAMS = 6
  Character(len=strLen), Parameter  :: cosexpansion_idstring = 'COSEXPANSION'
  
  Type CosExpansionModel
    Real(kind=RDbl), Dimension(0:MAX_PARAMS-1) :: cn
  End Type CosExpansionModel

  Interface cosexpansion_init
    module procedure cosexpansion_initparams
  End Interface
  
  Interface getinteraction
    Module Procedure cosexpansion_getinteraction
  End Interface
  
Contains
  
  !--------------------------------------------------------------
  ! Initializes the c arrary and fills it with data
  ! Requires:  params -- Cos Expansion torsion parameters 
  !--------------------------------------------------------------
  Subroutine cosexpansion_initparams(cparams,params)
    Type(CosExpansionModel) :: cparams
    Character(*), Dimension(:) :: params
    
    Integer :: i
    
    Do i = 1, MAX_PARAMS
      cparams%cn(i-1) = toreal(params(i))
    End Do
    
  End Subroutine Cosexpansion_initparams
  
  !-------------------------------------------------------------
  ! Displays the contents of the cparams
  ! Requires:  params -- Cos Expansion torsion parameters 
  !-------------------------------------------------------------
  Function cosexpansion_display(cparams)
    Character(lstrLen)          :: cosexpansion_display
    Type(CosExpansionModel)     :: cparams

    Integer :: i
    
    Write(cosexpansion_display,'(1x,6f7.3)')  &
        (cparams%cn(i), i=Lbound(cparams%cn,1),Ubound(cparams%cn,1))

  End Function Cosexpansion_display
  
  !--------------------------------------------------------------
  ! Evaluation the cosexpansion interactions
  ! Requires:  params -- Cos Expansion torsion parameters 
  !            r -- position vectors for atoms
  !            u -- potential energy
  !            vf -- force on each atom (optional)
  !--------------------------------------------------------------
  Subroutine cosexpansion_getinteraction(params,r,ifvec,u,vf)
    Type(CosExpansionModel), Intent(In)                :: params
    Type(VecType), Intent(IN), Dimension(:)            :: r
    Logical, Intent(In) :: ifvec
    Real(kind=RDbl), Intent(OUT)                       :: u
    Type(VecType), Intent(OUT), Dimension(:), Optional :: vf

    Integer                :: i
    Real(kind=RDbl)        :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)        :: T1,T2,T3,T1T2,T1T3,T2T3_invroot
    Real(kind=RDbl)        :: cos,costerm,sum,ffac
    Type(VecType),Dimension(4)        :: dcos_dr
    Type(VecType)   :: a,b,c
    Type(VecType)   ::dcos_da,dcos_db,dcos_dc

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
    T1 = a_b*b_c - a_c*b_b
    T2 = a_a*b_b - a_b*a_b
    T3 = b_b*c_c - b_c*b_c
    If (Abs(T2*T3)<1.0e-11) Then
      ! rare occurence of colinear placement of 4 atoms?
      Write(*,*) "Problem in cosexpansion, (torsion)"
      Write(*,*) "rare occurence of colinear placement of 4 atoms?"
      u=zero
      If (Present(vf)) Then
        Do i=1,4
          vf(i) = zero
        End Do
      Endif
      Return
    Endif
    T2T3_invroot = one/Sqrt(T2*T3)
    T1T2 = T1/T2
    T1T3 = T1/T3

    !** Note the sign, When theta is zero the potential enrgy is minimum(zero)
    !** theta=0, when atom.1 and atom.4 are farthest from each other.
    !** cos(theta) is defined so that the above condn on theta is satisfied
    cos = -(T1*T2T3_invroot)

    !** Calculate potential
    costerm = one
    sum = zero
    Do i = 0,MAX_PARAMS-1
      sum = sum + params%cn(i)*costerm
      costerm = costerm*cos
    End Do

    !** torsional potential
    u = sum

    !** Calculate forces if requested
    If (Present(vf)) Then
      !** dcos /da, dcos /db, dcos /dc
      dcos_da=( (a*(T1T2*b_b)) + (c*b_b)- (b*(b_c+a_b*T1T2)) )* T2T3_invroot
      dcos_db=( (b*(2*a_c+T1T2*a_a+T1T3*c_c)) - (a*(b_c+a_b*T1T2)) - &
          (c*(a_b+T1T3*b_c)) ) * T2T3_invroot
      dcos_dc=( (a*b_b) +(c*(T1T3*b_b)) -(b*(a_b+b_c*T1T3)) ) * T2T3_invroot

      !** The derivative of cos with respect to position vector of 
      !** particles 1..4 : dcos /dr1, dcos /dr2 , dcos /dr3 and dcos /dr4
      dcos_dr(1)=dcos_da*(-one)
      dcos_dr(2)=dcos_da-dcos_db
      dcos_dr(3)=dcos_db-dcos_dc
      dcos_dr(4)=dcos_dc
      
      !the pre-factor for forces
      costerm=one
      sum=zero
      Do i=1,MAX_PARAMS-1
        sum=sum+dble(i)*params%cn(i)*costerm
        costerm=costerm*cos
      End Do

      !** minus sign arising from   F= -Dv/Dr
      ffac = - sum

      !** vector forces on particles 1..4
      vf(1)=dcos_dr(1)*(ffac)
      vf(2)=dcos_dr(2)*(ffac)
      vf(3)=dcos_dr(3)*(ffac)
      vf(4)=dcos_dr(4)*(ffac)

      !**ah..... those were the forces

    Endif

  End Subroutine cosexpansion_getinteraction
  
  !-------------------------------------------------------------------
  ! cleanup the c parameters
  ! Requires:  params -- Cos Expansion torsion parameters to clean
  !-------------------------------------------------------------------
  Subroutine cosexpansion_cleanup(cparams)
    Type(CosExpansionModel)        :: cparams

    !** nothing to do here

  End Subroutine cosexpansion_cleanup
  
End Module cosexpansion













