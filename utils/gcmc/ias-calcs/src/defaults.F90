!------------------------------------------------------------------------------
! This module contains upper limits on arrays difficult to allocate on the fly.
! It also contains the filename tag definitions
!------------------------------------------------------------------------------
Module defaults
  Implicit None
  Save

  Public 

  !-------------------------------------------------
  ! Input file comment character
  !-------------------------------------------------
  Character(len=1), Parameter :: COMMENT_CHAR = "#"

  !--------------------------------------------
  ! Type Definitions
  !--------------------------------------------
  Integer, Parameter   :: RDbl = Selected_real_kind(10, 50)
  Integer, Parameter   :: RSgl = Selected_real_kind(5,10)
  Integer, Parameter   :: IDbl = Selected_int_kind(10)
  Integer, Parameter   :: ISgl = Selected_int_kind(4)
  Integer, Parameter   :: strLen = 48 ! Default character string length
  Integer, Parameter   :: lstrLen = 2*strLen ! Long string length
  Integer, Parameter   :: xlstrLen = 3*strLen ! Extra long string length
  Real(kind=RDbl), Parameter  :: zeroTolerance = 1.0e-6_Rdbl ! For comparing
                                                             ! Reals

  !--------------------------------------------
  ! Mathematical Constants
  !--------------------------------------------
  Real(kind=RDbl), Parameter  :: pi    = 3.1415927_RDbl
  Real(kind=RDbl), Parameter  :: twopi = 2.0_RDbl*pi
  Real(kind=RDbl), Parameter  :: degTorad = pi/180.0_RDbl
  Real(kind=RDbl), Parameter  :: radTodeg = 1.0_RDbl/degTorad
  Real(kind=RDbl), Parameter  :: zero = 0.0_RDbl
  Real(kind=Rdbl), Parameter  :: one = 1.0_RDbl

  !----------------------------------------------------------------------------
  ! Physical Constants
  ! Note that the internal units are as follows:
  ! length = Ang
  ! time   = ps
  ! mass   = amu
  ! temperature = K
  ! velocity = Ang/ps
  ! accelerations = Ang/ps^2
  ! energy = kcal (when calculated)
  ! energy = kJ (when stored and reported)
  !----------------------------------------------------------------------------
  Character(len=strLen), Parameter :: lengthUnit = "Ang"
  Character(len=strLen), Parameter :: timeUnit = "ps"
  Character(len=strLen), Parameter :: massUnit = "amu"
  Character(len=strLen), Parameter :: TUnit = "K"
  Character(len=strLen), Parameter :: nrgUnit = "kJ/mol"
  Character(len=strLen) :: velUnit = Trim(lengthUnit)//"/"//Trim(timeUnit)
 
  Real(kind=RDbl), Parameter  :: Nav = 6.0221367E+23_RDbl   ! Avogadro's Number
  Real(kind=RDbl), Parameter  :: bohrToAng  = 0.529177_RDbl
   ! HartreesTokJ/mol
  Real(kind=RDbl), Parameter  :: HTokJmol = 4.3598E-21_RDbl*Nav 
  Real(kind=RDbl), Parameter  :: calToj = 4.184_RDbl  
  Real(kind=RDbl), Parameter  :: eVTokcalmol = 23.0605_RDbl  

  Real(kind=RDbl), Parameter :: Rgas    = 8.3144_Rdbl ! J/mol/K
  Real(kind=RDbl), Parameter :: hplanck = 6.625E-34_RDbl   ! J.s
  Real(kind=RDbl), Parameter :: scalepe = 4.184_RDbl ! (k)cal -> (k)J Unit conv
  Real(kind=RDbl), Parameter :: scalef =418.4_RDbl ! kcal/mol/Ang->amu*Ang/ps^2
  Real(kind=RDbl), Parameter :: scaleke = 0.01_RDbl   !(Ang/ps)^2*amu -> KJ/mol
  Real(kind=RDbl), Parameter :: jmolec_kb    = Rgas/Nav       ! J/molecule/K
  Real(kind=RDbl), Parameter :: kjmole_kb    = 8.314E-3_RDbl ! R, KJ/mole/K
  Real(kind=RDbl), Parameter :: jmole_kb     = 8.314_RDbl         ! J/mole/K
  Real(kind=RDbl), Parameter :: kcalmole_kb  = 1.9865E-3_RDbl ! kcal/mole/K
  Real(kind=RDbl), Parameter :: calmole_kb   = 1.9865_RDbl         ! cal/mole/K

  !--------------------------------------------
  ! Maxima
  !--------------------------------------------

  !maximum to be used as the argument of exponential functions
  Real(kind=RDbl), Parameter :: default_MAX_EXP =30.000_RDbl

  !--------------------------------------------
  ! Filename "Tag" definitions
  !--------------------------------------------
  Character(len=strLen)   :: d_ctrl_file = 'general control file'

  !--------------------------------------------
  ! Miscellaneous strings
  !--------------------------------------------
  Character(len=60)       :: dashedline  = &
      '------------------------------------------------------------'
  Character(len=60)       :: dashedline2  = &
      '  ----------------------------------------------------------'

  !file stuff
  Integer, Parameter:: FIRST_UNIT=12
  Integer, Parameter::MAX_FILES = 100

End Module defaults


