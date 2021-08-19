!------------------------------------------------------------------------------
! Stores the information about the Gaussian Constraint thermostat parameters
!------------------------------------------------------------------------------
Module gaussts

  Use defaults, Only: RDbl, strLen, lstrLen, STATS_BLOCKSIZE, one, zero
  Use utils, Only: stripcmnt, split, toint, toupper, allocErrDisplay
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/)
  Use molecules, Only: molecules_getnsorbs
  Use stats, Only: Statistics, stats_init

  Implicit None

  Private
  Public :: GausstsParams, gaussts_init, gaussts_display

  Type GausstsParams
    Real(kind=RDbl) :: Tk
    Logical :: absolute         
    ! thermostat wrt absolute or streaming velocities 
  End Type GausstsParams

Contains

  !----------------------------------------------------------------------------
  ! Initializes the parameters for Nose-Hoover
  !----------------------------------------------------------------------------
  Subroutine gaussts_init(gauss,unitno)
    Type(GausstsParams), Intent(InOut) :: gauss
    Integer, Intent(in) :: unitno
    Character(len=strLen) :: text

    Read(unitno,'(a)') text
    text = stripcmnt(text)
    Select Case(Trim(toupper(text)))
    Case("ABSOLUTE")
      gauss%absolute=.True.
    Case("STREAMING")
      gauss%absolute=.False.
    Case default
      Write(*,*) "could not undrstand line -- "//&
          Trim(text)//" -- forgauss thermostat"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End Select
  End Subroutine gaussts_init

  !----------------------------------------------------------------------------
  ! Returns a string with the contents of Nose-Hoover variable
  !----------------------------------------------------------------------------
  Subroutine gaussts_display(gauss,dUnit)
    Type(GausstsParams), Intent(In) :: gauss
    Integer, Intent(In), Optional :: dUnit
    Write(dUnit,'(a)') "Gauss Thermostat Parameters"
  End Subroutine gaussts_display

End Module gaussts

