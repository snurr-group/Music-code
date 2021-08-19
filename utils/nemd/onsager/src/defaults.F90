!------------------------------------------------------------------------------
! This module contains upper limits on arrays difficult to allocate on the fly.
! It also contains the filename tag definitions
!------------------------------------------------------------------------------
Module defaults
  Implicit None
  Save

  Public 

  Integer :: pcalls

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

  !--------------------------------------------
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
  !--------------------------------------------
  Character(len=strLen), Parameter :: lengthUnit = "Ang"
  Character(len=strLen), Parameter :: timeUnit = "ps"
  Character(len=strLen), Parameter :: massUnit = "amu"
  Character(len=strLen), Parameter :: TUnit = "K"
  Character(len=strLen), Parameter :: nrgUnit = "kJ/mol"
  Character(len=strLen) :: velUnit = Trim(lengthUnit)//"/"//Trim(timeUnit)

  Real(kind=RDbl), Parameter  :: Echarge = 1.60E-19 
  ! electron chrage in coulombs

  Real(kind=RDbl), Parameter  :: Nav = 6.0221367E+23_RDbl   ! Avogadro's Number
  Real(kind=RDbl), Parameter  :: calToj = 4.184_RDbl  
  Real(kind=RDbl), Parameter  :: bohrToAng  = 0.529177_RDbl
  Real(kind=RDbl), Parameter  :: HTokJmol = 4.3598E-21_RDbl*Nav !Hartrees->kJ/mol
  Real(kind=RDbl), Parameter  :: HTokcalmol = HTokJmol/calToj   !H->kcal/mol
  Real(kind=RDbl), Parameter  :: eVTokcalmol = 23.0605_RDbl  

  Real(kind=RDbl), Parameter :: Rgas    = 8.314_Rdbl ! J/mol/K
  Real(kind=RDbl), Parameter :: hplanck = 6.625E-34_RDbl   ! J.s
  Real(kind=RDbl), Parameter :: scalepe = 4.184_RDbl ! (k)cal -> (k)J Unit conv
  Real(kind=RDbl), Parameter :: scalef  = 418.4_RDbl 
  ! Actual conversion : kcal -> amu. ang^2 / ps^2      (energy)
  !  OR                 kcal/ang -> amu. ang/ps^2      (force)
  Real(kind=RDbl), Parameter :: scaleke = 0.01_RDbl   !(Ang/ps)^2*amu -> KJ/mol
  Real(kind=RDbl), Parameter :: jmolec_kb    = Rgas/Nav       ! J/molecule/K
  Real(kind=RDbl), Parameter :: kjmole_kb    = 8.314E-3_RDbl ! R, KJ/mole/K
  Real(kind=RDbl), Parameter :: jmole_kb     = 8.314_RDbl         ! J/mole/K
  Real(kind=RDbl), Parameter :: kcalmole_kb  =  kjmole_kb/calToj! kcal/mole/K
  Real(kind=RDbl), Parameter :: calmole_kb   = jmole_kb/calToj         ! cal/mole/K

  !--------------------------------------------
  ! Maxima
  !--------------------------------------------
  !Max. number of maps. of each kind i.e. you can 
  !have MAX_MAPS pmaps, MAX_MAPS bias maps etc.
  Integer, Parameter   :: MAX_MAPS = 10

  !Max. number of sites that can be excluded
  Integer, Parameter   :: MAX_EXCLSITES = 10
  Integer, Parameter   :: MAX_SITES = 10


  !Max. number of atom types
  Integer, Parameter   :: MAX_ATOMTYPES = 20

  !Max. number of zeolite atom TYPES
  Integer, Parameter   :: MAX_ZEOATOMTYPES = 15

  !Max. number of sorbate types
  Integer, Parameter   :: MAX_SORBS = 10

  !Max. number of molecules of each sorbate type
  Integer, Parameter   :: MAX_MOLECULES = 1  

  ! Initialize the arrays to this size if no. of molecules
  ! is not specified
  Integer, Parameter   :: INITIAL_SIZE = 800

  ! The maximum number of atoms in a molecule -- AVOID USING THIS
  Integer, Parameter   :: MAX_ATOMS    = 40

  ! The maximum degrees of freedom in a molecule
  Integer, Parameter   :: MAX_DOF      = 20

  !Amount by which to increment the maximum no. of molecules
  !that can be stored in COORDINATE structure once the max. is exceeded
  Integer, Parameter   :: INCR_SIZE = 100

  !Amount by which to increment the maximum no. of molecules
  !that can be stored in INTERACTIONS structure once the max. is exceeded
  Integer, Parameter   :: INCR_STORAGE_SIZE = 2

  !The first file unit number used
  Integer, Parameter   :: FIRST_UNIT = 12

  !Max. number of files
  Integer, Parameter   :: MAX_FILES = 1200

  !Max. number of datafiles or configfiles
  Integer, Parameter   :: MAX_DAT_FILES = 30

  !maximum to be used as the argument of exponential functions
  Real(kind=RDbl), Parameter :: default_MAX_EXP =30.000_RDbl

  !The default size of the blocksize for keeping statistics
  Integer, Parameter   :: STATS_BLOCKSIZE = 10000

  !The maximum no. of lines possible in the control file
  Integer, Parameter   :: MAX_CTRL_LINES = 200

  !default cut-offs, can be used if nothing else specified
  Real(kind=RDbl), Parameter  :: LOW_DIST_CUTOFF = 1.0e-3_RDbl
  Real(kind=RDbl), Parameter  :: HIGH_DIST_CUTOFF = 13.0_RDbl

  !-----------------------------------------------------------
  ! Number of trial moves used in inserting branched molecules 
  !-----------------------------------------------------------
  !Number of trial moves in building up the angle libraries 
!  Integer, Parameter   :: lib_ntrials = 10000 
  !Number of trial moves in sampling the bond angles  
!  Integer, Parameter   :: ba_ntrials = 100
  !Number of trial moves in sampling the dihedral angles  
!  Integer, Parameter   :: ta_ntrials = 6  
  !Number of trial moves in computing the corrections added
  ! to the chemical potential due to intramolecular 
  ! interactions 
!  Integer, Parameter   :: correc_ntrials = 100000


  !--------------------------------------------
  ! Filename "Tag" definitions
  !--------------------------------------------
  Character(len=strLen)   :: d_ctrl_file = 'general control file'
  Character(len=strLen)   :: d_res_file  = 'restart file'
  Character(len=strLen)   :: d_con_file  = 'configuration file'
  Character(len=strLen)   :: d_dat_file  = ' data(config) file'
  Character(len=strLen)   :: d_ss_file   = 'sorbate pair file'
  Character(len=strLen)   :: d_aa_file   = 'atom pair file'
  Character(len=strLen)   :: d_mol_file  = 'molecule definition file'
  Character(len=strLen)   :: d_uc_file   = 'zeolite definition file'

  !--------------------------------------------
  ! Miscellaneous strings
  !--------------------------------------------
  Character(len=20)       :: shortdashedline  = &
      '--------------------'
  Character(len=60)       :: dashedline  = &
      '------------------------------------------------------------'
  Character(len=60)       :: dashedline2  = &
      '  ----------------------------------------------------------'

  !------------------------------------------------------------------
  ! Integers corresponding to different routines.  Used for timing
  ! by wtime
  !------------------------------------------------------------------
  Integer, Parameter :: MAX_ROUTINES        = 20
  Integer, Parameter :: ROUT_MAIN           = 1
  Integer, Parameter :: ROUT_HERMITE        = 2    
  Integer, Parameter :: ROUT_GCMODELSINSERT = 3
  Integer, Parameter :: ROUT_MSDRIVER       = 4
  Integer, Parameter :: ROUT_MSDRIVER_MAP   = 5
  Integer, Parameter :: ROUT_MSDRIVER_BASIC = 6
  Integer, Parameter :: ROUT_GCMC           = 7

  !------------------------------------------------------------------
  ! Integers corresponding to index of specific intramolecular 
  ! interaction, used in calls to intramolecular routines
  ! NOTE: these should be put in intramolecular if the energy
  !       accumulation is ever moved out of molecules
  !------------------------------------------------------------------
  Integer, Parameter :: STRETCH_INDEX        = 1
  Integer, Parameter :: BENDING_INDEX        = 2
  Integer, Parameter :: TORSION_INDEX        = 3
  Integer, Parameter :: CONSTRAINT_INDEX     = 4
  Integer, Parameter :: INTRAPAIR_INDEX      = 5
  Integer, Parameter :: INTRACOUL_INDEX      = 6
  Integer, Parameter :: EXTERNAL_INDEX       = 7
  Integer, Parameter :: TOTAL_INDEX          = 8
  Integer, Parameter :: NO_OF_INTRA_POTS     = 8
  Integer, Parameter :: INTRA_NRG_BLOCKSIZE = 1000

  !------------------------------------------------------------------
  ! Just a flag that can make debugging easier
  !------------------------------------------------------------------
  Logical            :: dbgflag = .False.


  !----------------------------------------------------------------
  ! size of some typicall arrays used
  !________________________________________________________________
  Integer, Parameter :: SIZE_OF_HOT=8  ! hot arry contains: 
                                       ! V, Vx, Vy, Vz, Vxy, Vxz, Vyz, Vxyz, 
                                       ! Where V=nrg, Vx= dV/dX etc....

End Module defaults




