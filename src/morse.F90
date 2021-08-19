!-------------------------------------------------------------------------------
! Morse potential data structure, calculation routine and misc routines
! The form of the potential is:
!   U(r) = D*(1 - exp{-beta*(r - req)})^2
!-------------------------------------------------------------------------------

Module morse

  Use defaults, Only: RDbl, strLen, one
  Use vector, Only: VecType, mag, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use utils, Only: toreal

  Implicit None
  Save

  Private
  Public :: MorseParams, morse_idstring, morse_init, morse_getinteraction, &
      morse_display

  Integer :: NO_PARAMS = 3
  !** Model parameters
  !** Morse Potential: U(r) = D*(1 - exp{-beta*(r - req)})^2
  
  Type MorseParams
    Real(kind=RDbl)         :: D
    Real(kind=RDbl)         :: beta
    Real(kind=RDbl)         :: req
  End Type MorseParams

  Character(len=strLen), Parameter         :: morse_idstring = 'MORSE'

Contains
  
  !---------------------------------------------------------------------
  ! Initialization subroutine
  ! Requires: moparams -- type to be initialized
  !           params -- incoming strings to be used
  !           blength -- optional bond length (in Angstroms)
  ! Uses the passed bond length if it's present
  !---------------------------------------------------------------------
  Subroutine morse_init(moparams,params,blength)
    Type(MorseParams), Intent(OUT)       :: moparams
    Character(*), Dimension(:)           :: params
    Real(Kind=RDbl), Optional            :: blength

    Integer              :: i

    If (Size(params) /= NO_PARAMS) Then
      Write(0,'(1x,2a,i4,a,i3)') __FILE__,": ",__LINE__, & 
          " : Incorrect number of parameters for morse pot ",Size(params)
      Write(0,'(1x,5a)') "params line is : ", &
          (Trim(params(i))//'  ',i=1,Size(params))
      Stop
    End If
    
    moparams%D = toreal(params(1))
    moparams%beta = toreal(params(2))

    If (Present(blength)) Then
      moparams%req = blength
    Else 
      moparams%req = toreal(params(3))
    End If
    
  End Subroutine Morse_init
  
  
  !----------------------------------------------------
  ! Potential AND force calculation
  !----------------------------------------------------
  Subroutine morse_getinteraction(obj,coords,sepvec,u,vf)
    
    !** Morse Potential: U(r) = D*(1 - exp{-beta*(r - req)})^2
    !** Morse Force    : F(r) = -2D*(1-exp{-beta*(r-req)})*
    !**                             (-(-beta)*exp{-beta*(r-req)})
    !**                       = -2D*beta*exp{-beta*(r-req)}*
    !**                             (1-exp{-beta*(r-req)})
    
    !** Passed parameters
    Type(MorseParams), Intent(IN)                      :: obj
    Type(VecType), Intent(IN), Dimension(:)            :: coords
    Logical, Intent(In)                                :: sepvec
    Type(VecType), Optional, Intent(OUT), Dimension(:) :: vf ! force vector
    Real(kind=RDbl), Intent(OUT)                       :: u  ! pot energy
    
    !** Local variables
    Real(kind=RDbl)   :: r  ! magnitude(distance)
    Real(kind=RDbl)   :: D  ,beta, req! morse params
    Type(VecType)     :: vr,vru ! distance vector,unit vector in same dirn.
    Real(kind=RDbl)   :: f  ! magnitude(force)
    
    !** Check to see if coordinates or sepvecs were passed
    If (sepvec) Then
      !** It's the separation vector
      vr = coords(1)
    Else
      !** Get the distance vector and its magnitude
      vr = coords(2)-coords(1)
    End If

    r = mag(vr)
    vru=vr/r
    
    D=obj%D
    beta=obj%beta
    req=obj%req
    
    !** Calculate the potential and the force
    u = D*(1-Exp(-beta*(r-req)))**2
    f = -2*D*beta*Exp(-beta*(r-req))*(1-Exp(-beta*(r-req)))
    
    !** Convert magnitude to a vector
    If (Present(vf)) Then
      vf(2) = vru*f
      vf(1) = vf(2)*(-one)
    End If
    
  End Subroutine morse_getinteraction
  
  !----------------------------------------------------
  ! Display morse potential params
  !----------------------------------------------------
!!$  Character(len=strLen) Function morse_display(moparams)
  Function morse_display(moparams)
    Character(len=strLen) :: morse_display
    Type(MorseParams), Intent(IN) :: moparams
    Write(morse_display,'(3(1x,e10.4))') moparams%D, moparams%beta, &
        moparams%req
  End Function Morse_display
  
End Module Morse





