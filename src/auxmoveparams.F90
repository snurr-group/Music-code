!-------------------------------------------------------
! Usually the moves modules needs objects like : sorbates, forcefield,
! movetypes, simcell. Any other auxiliary object that needs to passed to moves 
! is passed through this module 
! currently the only auxiliary object we have is cavitylist
!-------------------------------------------------------

Module auxmoveparams

  Use cavitylist, Only: Cavity_Params,cavitylist_init
  Use config, Only: AtMolCoords
  Use defaults, Only: strLen, lstrLen, RDbl, RSgl, MAX_MAPS, dashedline,&
      one, zero, d_ctrl_file, MAX_ATOMS
  Use file, Only: file_open
  Use simcell, Only: SimCell_Params
  Use utils, Only: split, getpath, allocErrDisplay, filesrchstr, int2str, &
      stripcmnt

  Implicit None
  Save

  Private
  Public :: AuxMoveObjects,  auxmoveparams_init


  Type AuxMoveObjects
    Type(Cavity_Params),Pointer :: cavity
    Logical :: cavbiasON
  End Type AuxMoveObjects

Contains
  !---------------------------------------------------------------
  ! initialize auxmv object
  ! requires :
  ! auxmv - 
  ! filename - main control file
  !---------------------------------------------------------------
  Subroutine auxmoveparams_init(auxmv, species, scell, filename )
    Type(AuxMoveObjects),Pointer :: auxmv
    Character(*), Intent(In)                       :: filename
    Type(AtMolCoords), Dimension(:), Intent(in)     :: species
    Type(SimCell_Params),Intent(in)  :: scell
    Logical :: found_cav
    Allocate(auxmv)
    
    auxmv%cavbiasON=.False.
    Call cavitylist_init(auxmv%cavity, species, scell, filename, found_cav )
    auxmv%cavbiasON= found_cav

  End Subroutine auxmoveparams_init

End Module auxmoveparams
