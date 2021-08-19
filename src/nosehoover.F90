!------------------------------------------------------------------------------
! Stores the information about the Nose-Hoover thermostat parameters
!------------------------------------------------------------------------------
Module nosehoover

  Use defaults, Only: RDbl, strLen, lstrLen, STATS_BLOCKSIZE, one, zero
  Use utils, Only: stripcmnt, split, toint, allocErrDisplay
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use molecules, Only: molecules_getnsorbs
  Use stats, Only: Statistics, stats_init

  Implicit None

  Private
  Public :: NHCoords, NHSorbates, NoseHooverParams, nosehoover_init, &
      nosehoover_display, nosehoover_getcoord, nosehoover_sampleCF
  
  Type NHCoords
    Real(kind=RDbl) :: v
    Real(kind=RDbl) :: r
  End Type NHCoords

  Type NHSorbates
    Type(NHCoords) :: s    ! Extended coordinate
    Type(NHCoords) :: tf   ! Thermal friction coord
  End Type NHSorbates

  Type NoseHooverParams
    Real(kind=RDbl)  :: Q      ! Thermal inertia, one per system
    Type(NHSorbates) :: coords ! This makes it look like the coords in config
    Type(NHSorbates),Dimension(:), Pointer :: spc ! for species_wise thermostat
    Type(Statistics) :: Hnose  ! Nose-Hoover hamiltonian
    Logical :: species_wise
  End Type NoseHooverParams

  Character(len=strLen), Parameter :: nhUnits = "ps^2 kcal/mol"

Contains

  !----------------------------------------------------------------------------
  ! Initializes the parameters for Nose-Hoover
  !----------------------------------------------------------------------------
  Subroutine nosehoover_init(nh,unitno,opt_spc_flag)
    Type(NoseHooverParams), Intent(InOut) :: nh
    Integer, Intent(In) :: unitno
    Logical, Intent(in), Optional :: opt_spc_flag

    Character(len=strLen), Dimension(strLen) :: params
    Character(len=lstrLen) :: text
    Integer :: i, nmols, error

    nh%species_wise=.False.
    If (Present(opt_spc_flag)) Then
      If (opt_spc_flag) nh%species_wise=.True.
    Endif

    !** Read the parameters
    Read(unitno,'(a)') text
    text = stripcmnt(text)
    i = split(text,params)
    nh%Q = toint(params(1))

    Call stats_init(nh%Hnose,"H(Nose)",STATS_BLOCKSIZE,.False.,'f10.3')

    If (nh%species_wise) Then
      nmols=molecules_getnsorbs()
      Allocate(nh%spc(nmols),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
      Do i=1,nmols
        nh%spc(i)%s%r  = one
        nh%spc(i)%s%v  = zero
        nh%spc(i)%tf%r = zero
        nh%spc(i)%tf%v = zero
      Enddo
    Else
      nh%coords%s%r  = one
      nh%coords%s%v  = zero
      nh%coords%tf%r = zero
      nh%coords%tf%v = zero
    Endif
  End Subroutine nosehoover_init

  !----------------------------------------------------------------------------
  ! Writes a sample of the required control file information to unit unitno
  !----------------------------------------------------------------------------
  Subroutine nosehoover_sampleCF(unitno)
    Integer, Intent(In) :: unitno

    Write(unitno,'(a,t30,3a)') 'Real','# Q(Nose-Hoover), (',Trim(nhUnits),')' 
  End Subroutine nosehoover_sampleCF

  !----------------------------------------------------------------------------
  ! Returns the value of Q, the thermal intertia
  !----------------------------------------------------------------------------
  Real(kind=RDbl) Function nosehoover_Q(nh)
    Type(NoseHooverParams), Intent(In) :: nh
    
    nosehoover_q = nh%Q
  End Function nosehoover_Q

  !----------------------------------------------------------------------------
  ! Returns a pointer to the s (extended Nose-Hoover) coordinates in spntr
  !----------------------------------------------------------------------------
  Subroutine nosehoover_getcoord(nh,spntr)
    Type(NoseHooverParams), Intent(In), Target :: nh
    Type(NHSorbates), Pointer   :: spntr
      spntr => nh%coords
  End Subroutine nosehoover_getcoord

  !----------------------------------------------------------------------------
  ! Returns a string with the contents of Nose-Hoover variable
  !----------------------------------------------------------------------------
  Subroutine nosehoover_display(nh,dUnit)
    Type(NoseHooverParams), Intent(In) :: nh
    Integer, Intent(In), Optional :: dUnit

    Write(dUnit,'(a)') "Nose-Hoover Thermostat Parameters"
    Write(dUnit,'(a,f10.3,1x,a)') "Q(Nose) (Thermal Inertia) : ", &
        nh%Q,Trim(nhUnits)

  End Subroutine nosehoover_display

End Module nosehoover

