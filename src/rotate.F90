!------------------------------------------------------------------------
! This module acts as a common interface for rotation movements in 
! different coordinate systems.
!
! NOTE : As of now only rigid-random-rotate works.  What's the point of
! this module?  We only have one type of rotation move.
!------------------------------------------------------------------------
Module rotate
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen, one
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use rigidmoves, Only: RigidRRot_Params, &
      rigidmoves_RRotInit,rigidmoves_randomrotate,rigidmoves_adjustdeltarot,&
      rigidmoves_displayparams,rigidmoves_dyninfo, rigidmoves_postadjust
  Use brmoves, Only: BranchedRRot_Params, brmoves_RRotInit, &
      brmoves_randomrotate, brmoves_adjustdeltarot, brmoves_displayparams, &
      brmoves_postadjust
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      checkandstop, allocErrDisplay

  Implicit None
  Save

  Private
  Public :: Rotate_Params, RRotate_Params, &
      rotate_move , rotate_display, rotate_init, &
      rotate_adjustdelta, rotate_dyninfo, rotate_idstrings, &
      rotate_postadjust, rotate_setmaxnrg

  !***  Main Rotation Parameter Object
  Type Rotate_Params
    Type(RRotate_Params), Pointer     :: random 
  End Type Rotate_Params
  
  !** Random rotate params
  Type RRotate_Params
    Type(RigidRRot_Params), Pointer :: rigid
    Type(BranchedRRot_Params), Pointer  :: branched
  End Type RRotate_Params

  Character(len=strLen), Dimension(1), Parameter :: rotate_idstrings = &
      (/'RROTATE'/)
  
Contains


  !---------------------------------------------------------------------------
  ! Initializes the Rotate type move
  ! Requires:  params -- generic rotate object to be inited 
  !            movetype -- string identifying move type
  !            modeltype -- string identifying coordinate system model
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- other auxiliary stuff (ex. cavity bias)
  !---------------------------------------------------------------------------
  Subroutine rotate_init(params,movetype,modeltype,filename,spc, &
      species,simcell,auxmv)
    Type(Rotate_Params), Intent(Out)            :: params
    Character(*), Intent(In)                    :: movetype,modeltype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: simcell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    
    Integer      :: error
    Logical      :: nomatch
    
    !** Nullify all the pointer instances in params
    Call rotate_nullify(params)
    nomatch = .True.

    !** Call appropriate module/routine based on movetype and coordinate model
    Select Case (Trim(ToUpper(movetype)))
    Case (rotate_idstrings(1))   
      !** ------- Rotation by  a random amount -------- 

      Allocate(params%random, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random params")
      Nullify(params%random%rigid)
      Nullify(params%random%branched)

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%random%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random rigid")
        Call rigidmoves_RRotInit(params%random%rigid, simcell, auxmv, &
            spc, filename)
        nomatch = .False.

      Case('BRANCHED')
        Allocate(params%random%branched, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random branched")
        Call brmoves_RRotInit(params%random%branched, spc, auxmv, filename)
        nomatch = .False.

      End Select

      
    End Select

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get move-parameters for modeltype "//Trim(modeltype) )

  End Subroutine rotate_init

  !----------------------------------------------------------------------
  ! Set the maximum allowed energy
  ! Requires:  params -- generic translate object 
  !            max_nrg -- maximum energy (in kcal/mol)
  !----------------------------------------------------------------------
  Subroutine rotate_setmaxnrg(params,max_nrg)
    Type(Rotate_Params), Intent(InOut)   :: params
    Real(kind=RDbl), Intent(In)             :: max_nrg
    
    If (Associated(params%random)) Then
      If (Associated(params%random%rigid)) Then
        params%random%rigid%max_nrg = max_nrg
      Else
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' Unable to recognized sub-initialized pointer'
        Stop
      End If

    Else
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Unable to recognized initialized pointer'
      Stop
    End If

  End Subroutine rotate_setmaxnrg

  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic rotate pointer set
  !---------------------------------------------------------------------------
  Subroutine rotate_nullify(params)
    Type(Rotate_Params), Intent(InOut)          :: params

    Nullify(params%random)

  End Subroutine rotate_nullify

  !----------------------------------------------------------------------
  ! rotates a molecule and gets the  energies 
  ! Requires:  params  -- parameters for this rotate-move
  !            subset -- system subset to which to apply move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias used during rotateion of this molecule
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Logical Function rotate_move(params,subset,species,simcell, &
      biasfactor,subint)
    Type(Rotate_Params), Intent(Inout)              :: params
    Integer, Dimension(:), Intent(In)               :: subset
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Real(kind=RDbl), Intent(Out)                    :: biasfactor
    Type(Subset_Interactions), Intent(InOut)        :: subint

    Integer     :: natoms,atomno,spc,mol
    Logical     :: nomatch,nrg_calc_flag

    rotate_move = .False.
    nomatch = .True.
    spc = subset(1)
    mol = subset(2)

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%random)) Then   
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_randomrotate(params%random%rigid,species,spc,mol,&
            simcell,nrg_calc_flag, subint)
        biasfactor=one
        rotate_move = nrg_calc_flag
        nomatch = .False.
      Else If (Associated(params%random%branched)) Then
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End If

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not find appropriate move-parameters or modeltype ")

  End Function rotate_move

  !----------------------------------------------------------------------
  ! Displays parameters for rotations
  !----------------------------------------------------------------------
  Subroutine rotate_display(params, indent, unitno)
    Type(Rotate_Params), Intent(in)              :: params
    Integer, Intent(in)                             :: indent, unitno
    
    Logical :: nomatch
    nomatch = .True.

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%random)) Then    
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_displayparams(params%random%rigid,indent,unitno)
        nomatch = .False.
      Else If(Associated(params%random%branched)) Then
        Call rigidmoves_displayparams(params%random%rigid,indent,unitno)
        nomatch = .False.
      End If
    End If
    
    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get appropriate move-parameters or modeltype ")

  End Subroutine rotate_display

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Rotation parameters
  ! Requires:  params -- rotation parameters
  !----------------------------------------------------------------------
  Function rotate_dyninfo(params)
    Character(len=lstrLen)              :: rotate_dyninfo
    Type(Rotate_Params), Intent(In)     :: params

    rotate_dyninfo  = Repeat(' ',strLen)
    If (Associated(params%random))  Then 
      If (Associated(params%random%rigid))  &
          rotate_dyninfo  = rigidmoves_dyninfo(params%random%rigid)
    End If

  End Function rotate_dyninfo

  !----------------------------------------------------------------------
  ! adjusts parameters for Random rotations
  !----------------------------------------------------------------------
  Subroutine rotate_adjustdelta(params,ratio)
    Type(Rotate_Params), Intent(inout)           :: params
    Real(kind=RDbl)                              :: ratio

    If (Associated(params%random))  Then 
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_adjustdeltarot(params%random%rigid,ratio)
      Else If (Associated(params%random%branched))  Then 
        Call brmoves_adjustdeltarot(params%random%branched,ratio)
      End If
    End If

  End Subroutine rotate_adjustdelta

  !----------------------------------------------------------------------
  ! does adjust ment s for cavity bias after a move is accepted 
  ! Its coded only for rigid moveas as of now
  !----------------------------------------------------------------------
  Subroutine rotate_postadjust(params, species, sorb, mol )
    Type(Rotate_Params), Intent(In)                 :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species    
    Integer, Intent(In)                             :: sorb, mol

    Logical :: do_adjust

    If (Associated(params%random)) Then
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_postadjust(params%random%rigid, species, sorb, mol) 
      elseif (Associated(params%random%rigid)) Then
        Call brmoves_postadjust(params%random%branched, species, sorb, mol) 
      endif
    End If

  End Subroutine rotate_postadjust
  
End Module rotate








