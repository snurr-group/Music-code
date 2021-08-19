!-----------------------------------------------
! Model for calculating torional interactions with flexible bonds
! see notes for understanding the formulas used here
! DREIDING form
!-----------------------------------------------
Module cosexpansiondr

  Use defaults, Only: strLen, RDbl, lstrLen, one, zero, pi
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag
  Use utils, Only: toreal

  Implicit None
  Save

  Private
  Public :: CosExpansionModelDr, cosexpansiondr_init, cosexpansiondr_display, &
      cosexpansiondr_getinteraction

  Integer, Parameter :: MAX_PARAMS = 6


  Type CosExpansionModelDr
    Real(kind=RDbl), Dimension(0:MAX_PARAMS-1) :: cn
  End Type CosExpansionModelDr

  Interface cosexpansiondr_init
    module procedure cosexpansiondr_initparams
  End Interface

  Interface cosexpansiondr_force
    module procedure coexpansiondr_force_perturb
  End Interface

  Interface getinteraction
    Module Procedure cosexpansiondr_getinteraction
  End Interface

Contains

  !-----------------------------------------------------
  ! Initializes the c arrary and fills it with data
  !-----------------------------------------------------
  Subroutine cosexpansiondr_initparams(cparams,params)

    Type(CosExpansionModelDr) :: cparams
    Character(*), Dimension(:) :: params

    Integer :: i

    Do i = 1, MAX_PARAMS
      cparams%cn(i-1) = toreal(params(i))
      !       print*, 'Params: ', i, cparams%cn(i-1)
    End Do

    ! cparams%cn(0)<=>Vjk [kcal/mol]
    ! cparams%cn(1)<=> n (multiplicity)
    ! cparams%cn(2)<=> phi0 (initial angle)



  End Subroutine Cosexpansiondr_initparams

  !-----------------------------------------------------
  ! Displays the contents of the cparams
  !-----------------------------------------------------
  Function cosexpansiondr_display(cparams)
    Character(lstrLen)          :: cosexpansiondr_display
    Type(CosExpansionModelDr)     :: cparams

    Integer :: i

    Write(cosexpansiondr_display,'(1x,6f7.3)')  &
        (cparams%cn(i), i=Lbound(cparams%cn,1),Ubound(cparams%cn,1))

  End Function Cosexpansiondr_display

  !-----------------------------------------------------
  ! params contains the co-effs of the cos expansion
  ! r(1:4) contains the 4 position of atom centers comprising the torsion list
  ! It returns u (potential energy) and vf(force vectors)
  !-----------------------------------------------------
  Subroutine cosexpansiondr_getinteraction(params,r,ifvec,u,vf)
    Type(CosExpansionModelDr), Pointer :: params
    Type(VecType), Intent(IN), Dimension(:) :: r
    Logical, Intent(In) :: ifvec
    Type(VecType), Intent(OUT), Dimension(:), Optional :: vf
    Real(kind=RDbl), Intent(OUT) :: u
    Integer                :: i
    Real(kind=RDbl)        :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)        :: T1,T2,T3,T1T2,T1T3,T2T3_invroot
    Real(kind=RDbl)        :: cosa,costerm,sum,ffac
    Type(VecType),Dimension(4)        :: dcos_dr
    Type(VecType)   :: a,b,c
    Type(VecType)   ::dcos_da,dcos_db,dcos_dc
    Real(kind=RDbl) :: phi

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
    !**Note the sign, When theta is zero the potential enrgy is minimum(zero)
    !** theta=0, when atom.1 and atom.4 are farthest from each other.
    !**cosa(theta) is defined so that the above condn on theta is satisfied
    cosa=(T1*T2T3_invroot)
    phi=acos(real(cosa))


    !**Calculate potential
    costerm=one
    sum=zero
    sum=0.5*params%cn(0)*(1-cos(params%cn(1)*(phi-pi*params%cn(2)/180.0)))
    !** torsional potential
    u=sum

    !       print*, '----------COSDREIDINGDR--------'
    !       print*, u, phi*180.0/pi

    If (Present(vf)) Then

      ! here is analytical form

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
      If(cosa==one.or.dsin(phi)==zero) then
        sum=zero
      else
        sum=0.5_RDbl*params%cn(0)*params%cn(1)*dsin(params%cn(1)*(phi-pi*params%cn(2)/180.0_RDbl))&
            /dsin(phi)
      End IF

      !** minus sign arising from   F= -Dv/Dr
      ffac = - sum

      !** vector forces on particles 1..4
      vf(1)=dcos_dr(1)*(ffac)
      vf(2)=dcos_dr(2)*(ffac)
      vf(3)=dcos_dr(3)*(ffac)
      vf(4)=dcos_dr(4)*(ffac)

      ! Numerical call
      !   Do i=1, 4
      !    Do j=1,3
      !    call coexpansiondr_force_perturb(params,i,j,r,force)
      !        vf(i)%comp(j)=force
      !       print*, i,j, force, vf(i)%comp(j)
      !    End Do
      !   End Do

    Endif
  End Subroutine cosexpansiondr_getinteraction

  !-----------------------------------------------------
  ! cleanup the c parameters
  !-----------------------------------------------------
  Subroutine cosexpansiondr_cleanup(cparams)

    Type(CosExpansionModelDr), Pointer :: cparams
    Integer :: error

    Deallocate(cparams,STAT=error)
    If (error /= 0) Then
      Write(0,'(2a,i4,2a)') __FILE__,": ",__LINE__, &
          " : Error deallocating torsional parameters"
      Stop
    End If

  End Subroutine cosexpansiondr_cleanup



  !-------------------------------------------------------
  ! Here we use a numerical apporach to calculate 
  ! force given a perturbation of coordinate coordp of atom atomp
  !-------------------------------------------------------
  Subroutine coexpansiondr_force_perturb(params,atomp,coordp,r,force)
    Type(CosExpansionModelDr), Pointer :: params
    Type(VecType), Intent(IN), Dimension(:) :: r
    Real(kind=RDbl), Intent(OUT) :: force
    Integer                :: atomp, coordp
    Real(kind=RDbl)        :: a_a,a_b,a_c,b_b,b_c,c_c
    Real(kind=RDbl)        :: T1,T2,T3,T1T2,T1T3,T2T3_invroot
    Real(kind=RDbl)        :: cosa,sum1,sum2,dsum,delta
    Type(VecType)   :: a,b,c
    Type(VecType), Dimension(1:4) :: rtemp
    Real(kind=RDbl) :: phi

    ! Here atomp is the number of perturbed atom (1:4) and coordp 
    ! is the number of perturbed coordinate (1:3)

    delta=0.00001_RDbl

    rtemp=r
    rtemp(atomp)%comp(coordp)=rtemp(atomp)%comp(coordp)+delta


    a=rtemp(2)-rtemp(1)
    b=rtemp(3)-rtemp(2)
    c=rtemp(4)-rtemp(3)

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
    !**Note the sign, When theta is zero the potential enrgy is minimum(zero)
    !** theta=0, when atom.1 and atom.4 are farthest from each other.
    !**cosa(theta) is defined so that the above condn on theta is satisfied
    cosa=(T1*T2T3_invroot)

    phi=acos(real(cosa))

    !**Calculate potential FIRST POINT

    sum1=zero
    sum1=0.5_RDbl*params%cn(0)*(1.0_RDbl-dcos(params%cn(1)*(phi-pi*params%cn(2)/180.0_RDbl)))


    rtemp(atomp)%comp(coordp)=rtemp(atomp)%comp(coordp)-2.0_RDbl*delta

    a=rtemp(2)-rtemp(1)
    b=rtemp(3)-rtemp(2)
    c=rtemp(4)-rtemp(3)

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
    !**Note the sign, When theta is zero the potential enrgy is minimum(zero)
    !** theta=0, when atom.1 and atom.4 are farthest from each other.
    !**cosa(theta) is defined so that the above condn on theta is satisfied
    cosa=(T1*T2T3_invroot)

    phi=acos(real(cosa))

    !**Calculate potential SECOND POINT
    sum2=zero
    sum2=0.5_RDbl*params%cn(0)*(1.0_RDbl-dcos(params%cn(1)*(phi-pi*params%cn(2)/180.0_RDbl)))

    ! **Calculate DIFFERENCE in energy over 2 deltas
    dsum=sum1-sum2
    force=-dsum/(2.0_RDbl*delta)

  End Subroutine coexpansiondr_force_perturb

End Module cosexpansiondr

