!-------------------------------------------------------------------------------
! Harmonic potential data structure, calculation routine and misc routines
! The form of the potential is:
!    Harmonic Potential: U(r) = k/2*(r - req)^2
!    Harmonic Force    : F(r) = -k/2*(2r - 2req)
!-------------------------------------------------------------------------------
Module harmonic

  Use defaults, Only: RDbl, strLen
  Use utils, Only: toreal
  Use vector, Only: VecType, mag, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use storebase, Only: EnergyPlus,storebase_inc,storebase_nderivs

  Implicit None 
  Save

  Private
  Public :: HarmonicParams, harmonic_idstring, harmonic_init, &
      harmonic_simpleinit, harmonic_getinteraction, harmonic_snglint, &
      harmonic_display

  Integer :: NO_PARAMS = 2

  Type HarmonicParams
    Real(kind=RDbl)      :: k
    Real(kind=RDbl)      :: req
  End Type HarmonicParams

  Character(len=strLen), Parameter         :: harmonic_idstring = 'HARMONIC'

Contains

  !----------------------------------------------------------------------
  ! Initialize the potential parameters
  ! Requires: harparams -- type to be initialized
  !           params -- incoming strings to be used
  !           blength -- optional bond length (in Angstroms)
  ! Uses the passed bond length if it's present
  !----------------------------------------------------------------------
  Subroutine harmonic_init(harparams,params,blength)
    Type(HarmonicParams), Pointer          :: harparams
    Character(*), Dimension(:), Intent(In) :: params
    Real(Kind=RDbl), Optional              :: blength

    Integer              :: i

    If (Size(params) < NO_PARAMS) Then
      Write(0,'(1x,2a,i4,a)') __FILE__,": ",__LINE__, & 
          " : Incorrect number of parameters for harmonic pot"
      Write(0,'(1x,5a)') "params line is : ", &
          (Trim(params(i))//'  ',i=1,Size(params))
      Stop
    End If

    harparams%k = toreal(params(1))

    If (Present(blength)) Then
      harparams%req = blength
    Else 
      harparams%req = toreal(params(2))
    End If

  End Subroutine harmonic_init

  !----------------------------------------------------------------------
  ! Initialize the potential parameters from given numbers
  ! Requires:  harparams -- type to be initialized
  !            params -- array containing spring const and equilibrium len
  !----------------------------------------------------------------------
  Subroutine harmonic_simpleinit(harparams,params)
    Type(HarmonicParams), Intent(Out)         :: harparams
    Real(kind=RDbl), Dimension(2), Intent(In) :: params

    harparams%k = params(1)
    harparams%req = params(2)

  End Subroutine harmonic_simpleinit

  !---------------------------------------------------------------------------
  ! An interface for _getinteraction that takes the EnergyPlus structure. 
  ! This function has the same call structure as the pairwise interaction
  ! _snglint routines.  'results2' must be present if forces are requested.
  ! Requires:  params -- harmonic parameters structure
  !            sepvec -- atom-atom separation vector
  !            results -- resultant energies and derivatives
  !            results2 -- resultant energies and opposite forces
  !---------------------------------------------------------------------------
  Logical Function harmonic_snglint(params,sepvec,results,results2)
    Type(HarmonicParams), Intent(In)          :: params
    Type(VecType), Intent(In)                 :: sepvec
    Type(EnergyPlus), Intent(InOut)           :: results
    Type(EnergyPlus), Intent(InOut), Optional :: results2

    Real(kind=RDbl)                :: nrg
    Type(VecType), Dimension(2)    :: force

    harmonic_snglint = .True.

    If (storebase_nderivs(results) > 0) Then
      Call harmonic_getinteraction(params,(/sepvec/),.True.,nrg,force)
      Call storebase_inc(results,nrg)
      Call storebase_inc(results,force(1))
      Call storebase_inc(results2,force(2))
    Else
      Call harmonic_getinteraction(params,(/sepvec/),.True.,nrg)
      Call storebase_inc(results,nrg)
    End If

  End Function harmonic_snglint

  !---------------------------------------------------------------------------
  ! Calculate the potential AND the optional force 
  ! Requires: obj -- harmonic parameters structure
  !           coords -- either two atomic coordinates or a separation vector
  !           sepvec -- flag indicating that separation vectors are passed
  !           u -- potential 
  !           vf -- optional forces on the two atoms
  ! NOTE: uses expansions   U(r) = k/2*(r^2 - 2r*req + req^2)
  !                         F(r) = -k*(r-req)
  !---------------------------------------------------------------------------
  Subroutine harmonic_getinteraction(obj,coords,sepvec,u,vf)
    Type(HarmonicParams), Intent(In)                   :: obj
    Type(VecType), Dimension(:), Intent(In)            :: coords
    Logical, Intent(In)                                :: sepvec
    Real(kind=RDbl), Intent(Out)                       :: u
    Type(VecType), Dimension(:), Intent(Out), Optional :: vf 

    Real(kind=RDbl)                :: r     ! magnitude(distance)
    Real(kind=RDbl)                :: req,k ! harmonic params
    Real(kind=RDbl)                :: f     ! magnitude(force)
    Type(VecType)                  :: vr    ! distance vector
    Type(VecType)                  :: vru   ! unit 
    
    !** Zero the force and pot
    u = 0.0_RDbl
    If (Present(vf)) vf = VecType(0.0_RDbl)

    k = obj%k
    req = obj%req

    !** Determine whether coordinates or separation vectors have been passed
    If (sepvec) Then
      vr = coords(1)
    Else
      vr = coords(2) - coords(1)     ! direction toward (2)
    End If

    r = mag(vr)

    !** If r is exactly zero, we must set vru to 0 
    !** avoid a divide by zero error.
    If (r == 0.0_RDbl) Then
      vru = 0.0_RDbl
    Else
      vru = vr/r
    End If

    u = 0.5*k*(r-req)*(r-req)
    f = -k*(r-req)

    If (Present(vf)) Then
      vf(2) = vru*f
      vf(1) = vru*f*(-1.0d0)
    End If

  End Subroutine harmonic_getinteraction

  !---------------------------------------------------------
  ! Displays the values of the parameters
  ! Requires: harparams -- harmonic parameters structure
  !---------------------------------------------------------
  Function harmonic_display(harparams)
    Character(len=strLen) :: harmonic_display
    Type(HarmonicParams), Intent(IN) :: harparams

    Write(harmonic_display,'(1x,e10.4,1x,e10.4)') harparams%k,harparams%req

  End function Harmonic_display

End Module harmonic










