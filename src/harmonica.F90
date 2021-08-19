!-----------------------------------------------------------------
! Harmonic angle module
! Calculates the potential and forces associated with bond bending
! for harmonic potential of the bond angle
! It works with flexible bonds too!!!!
!-----------------------------------------------------------------

Module harmonica

  Use utils, Only: toreal
  Use defaults, Only: RDbl, strLen, degTorad, zero
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), mag

  Implicit None
  Save

  Private
  Public :: HarmonicAModel, harmonica_init, harmonica_getinteraction, &
      harmonica_display, harmonica_getparams,HarmonicaString

  Integer, Parameter :: NO_PARAMS = 2
  Type HarmonicAModel
    Real(kind=RDbl)              :: ktheta
    Real(kind=RDbl)              :: thetaeq
  End Type HarmonicAModel

  Character(len=strLen) :: HarmonicaString="HARMONICA"

Contains
  
  !------------------------------------------------------------------------
  ! Initialize the parameters; assign to variables
  ! Requires: hparams -- type to be initialized
  !           params -- incoming strings to be used
  !           angle -- optional angle (in radians!)
  ! Uses the passed angle if it's present
  !------------------------------------------------------------------------
  Subroutine harmonica_init(hparams,params,angle)
    Type(HarmonicAModel), Pointer :: hparams
    Character(*), Dimension(:)    :: params
    Real(Kind=RDbl), Optional     :: angle
    
    hparams%ktheta = toreal(params(1))
    If (Present(angle)) Then
      hparams%thetaeq = angle
    Else 
      hparams%thetaeq = toreal(params(2))*degTorad
    End If

  End Subroutine harmonica_init
  
  !----------------------------------------------------
  ! Display the contents of params (angle in degrees)
  !----------------------------------------------------
  Function harmonica_display(hparams)
    Character(len=strLen)         :: harmonica_display
    Type(HarmonicAModel), Pointer :: hparams
    
    Write(harmonica_display,'(1x,e10.4,1x,f8.3)') hparams%ktheta, &
        hparams%thetaeq/degTorad
    
  End Function harmonica_display

  !----------------------------------------------------
  ! Returns the params as a 2D array (angle in radians )
  !----------------------------------------------------
  Subroutine harmonica_getparams(hparams,paramsArr)
    Type(HarmonicAModel), Pointer               :: hparams
    Real(kind=RDbl) , Dimension(2), Intent(out) ::paramsArr 
    paramsArr(1)=hparams%ktheta
    paramsArr(2)=hparams%thetaeq
  End Subroutine harmonica_getparams
  
  !--------------------------------------------------------
  ! Calculates the interaction energy and optionally forces
  !--------------------------------------------------------
  Subroutine harmonica_getinteraction(params,coords,sepvec,u,vf)

    Type(HarmonicAModel), Pointer :: params
    Type(VecType), Dimension(:), Intent(IN) :: coords
    Logical, Intent(In) :: sepvec
    Real(kind=RDbl), Intent(OUT) :: u
    Type(VecType), Dimension(:), Intent(OUT), Optional :: vf

    Type(VecType), Dimension(2) :: vr
    Real(kind=RDbl) :: r1, r12, r2, r22
    Real(kind=RDbl) :: theta, atheta, costheta,sintheta, adotb
    Real(kind=RDbl) :: rad_0,k

    k = params%ktheta
    rad_0 = params%thetaeq

    ! The case when angle=180 degrees leads to some irregularities. I
    ! think I have fixed these problems. So Lev's fix with DREIDING form
    ! might not be required here. So I have commented out them. If the
    ! DREIDING form is required then it can be put in a new module 
    ! - Shaji. 24/Aug/2002

!!$
!!$    If(rad_0/=180.00*degTorad) then

    !** Determine whether coordinates or separation vectors have been passed
    If (sepvec) Then
      !** The coordinates passed are actually the separation vectors
      vr(1) = coords(1)
      vr(2) = coords(2)
    Else
      !** The cooridnates have been passed, calculate the sep vecs
      vr(1) = coords(1)-coords(2) !vector towards atom.1
      vr(2) = coords(3)-coords(2) !vector towards atom.3
    End If

    !** dot product with itself to find norm of a vector
    r12 = vr(1)*vr(1)
    r22 = vr(2)*vr(2)

    r1=Sqrt(r12)
    r2=Sqrt(r22)

    !adotb= dotproduct of a  and b
    !theta= angle between a and b
    adotb= vr(2)*vr(1)
    costheta = adotb/(r1*r2)
    theta = dacos(costheta)

    !potential
    u = 0.5_RDbl*k*(theta-rad_0)*(theta-rad_0)

    !** Calculate force vector
    If (Present(vf)) Then
      sintheta= Sqrt(1.0_RDbl-costheta*costheta)

      If (Abs(sintheta)<1.00e-10) Then
        ! this is the case when angle=180, or 0. Force should be zero there
        atheta= zero
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Write(*,*) "There is a problem with angle bending potential"
        Write(*,*) "angle too close to 180 or 0"
        Write(*,'(a,f14.9)') "Sin(\theta) = ", sintheta
      Else
        atheta = k*(theta-rad_0)/(r1*r2*sintheta)
      Endif

      !** calc forces
      !** Refer to the documentation for notes on these equations
      !** 
      vf(1) = (vr(2)-vr(1)*(adotb/r12))*atheta
      vf(3) = (vr(1)-vr(2)*(adotb/r22))*atheta

      !newton's third law
      vf(2) = (vf(1)+vf(3))*(-1.0_RDbl)

    End If
!!$    
!!$    Else
!!$    !** If equilibrium angle is 0 we implement DREIDING way 
!!$    !** of handling things: 
!!$    !** Functional for of the potential is now : Energy=K*(1+cos(angle))
!!$    !** J. Phys Chem 94(26) 1990, Mayo et.al.
!!$
!!$  
!!$    
!!$    !** Determine whether coordinates or separation vectors have been passed
!!$    If (sepvec) Then
!!$      !** The coordinates passed are actually the separation vectors
!!$      vr(1) = coords(1)
!!$      vr(2) = coords(2)
!!$    Else
!!$      !** The cooridnates have been passed, calculate the sep vecs
!!$      vr(1) = coords(1)-coords(2) !vector towards atom.1
!!$      vr(2) = coords(3)-coords(2) !vector towards atom.3
!!$    End If
!!$
!!$    !** dot product with itself to find norm of a vector
!!$    r12 = vr(1)*vr(1)
!!$    r22 = vr(2)*vr(2)
!!$
!!$    r1=sqrt(r12)
!!$    r2=sqrt(r22)
!!$
!!$    !adotb= dotproduct of a  and b
!!$    !theta= angle between a and b
!!$    adotb= vr(2)*vr(1)
!!$    costheta = adotb/(r1*r2)
!!$    theta = dacos(costheta)
!!$    
!!$    !potential
!!$    u = k*(1.0_RDbl+costheta)
!!$    
!!$    !** Calculate force vector
!!$    If (Present(vf)) Then
!!$      atheta = -k
!!$      
!!$      !** calc forces
!!$      !** Refer to the documentation for notes on these equations
!!$      !** 
!!$      vf(1) = (vr(2)-vr(1)*(adotb/r12))*atheta
!!$      vf(3) = (vr(1)-vr(2)*(adotb/r22))*atheta
!!$      
!!$      !newton's third law??, I'm not sure
!!$      vf(2) = (vf(1)+vf(3))*(-1.0_RDbl)
!!$      
!!$      !       vf(2) = (  ( vr(1) + vr(2) ) * (-1.0_RDbl)    + &
!!$      !           ( vr(2)/r22 + vr(1)/r12 ) *adotb)*atheta
!!$      
!!$    End If
!!$      
!!$    End If
    
  End Subroutine harmonica_getinteraction
  
  
End Module harmonica







