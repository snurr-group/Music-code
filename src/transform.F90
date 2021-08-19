!------------------------------------------------------------------------
! This module acts as a common interface for transformation movements in 
! different coordinate systems.  Transformation movements are any moves
! that also change the internal configuration of a molecule.
!
! Currently available transformation moves:
!   cut and regrow - deletion of molecule part with biased regrowth
!   eigenmode perturb - perturbation along molecular non-rot/trans eigenmodes
!   internal rotation perturb - rotates internal atom groups
!
! Note that the coordinate systems in which these move types are 
! available is restricted.
!------------------------------------------------------------------------

Module transform
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only:  filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay
  Use random, Only: rranf
  Use molecules, Only: molecules_getgcmodeltype
  Use config, Only: AtMolCoords, config_getnatoms
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use brmoves, Only: BranchedRegrow_Params, brmoves_reGrowInit, &
      brmoves_reGrow,brmoves_displayparams, brmoves_postadjust
  Use emperturb, Only: EMPerturb_Params,emperturb_init,emperturb_move, &
      emperturb_display,emperturb_clean
  Use irotperturb, Only: IRotPerturb_Params,irotperturb_init,irotperturb_move, &
      irotperturb_dyninfo,irotperturb_adjust, &
      irotperturb_display,irotperturb_clean

  Implicit None
  Save

  Private
  Public :: Transform_Moves,transform_idstrings,transform_init, &
      transform_move,transform_setmaxnrg,transform_dyninfo, &
      transform_adjustdelta,transform_display,transform_clean, &
      Regrow_Params,EMode_Pointer_Set,IRot_Pointer_Set, transform_postadjust

  Character(len=strLen), Dimension(3), Parameter :: transform_idstrings = &
      (/'CUT_REGROW  ','EMODEPERTURB','IROTPERTURB '/)

  !** Collection of transformation move pointer sets
  Type Transform_Moves
    Type(ReGrow_Params), Pointer          :: regrow
    Type(EMode_Pointer_Set), Pointer      :: emode
    Type(IRot_Pointer_Set), Pointer       :: irot
  End Type Transform_Moves
  
  !** Regrowth transform params
  Type Regrow_Params
    Type(BranchedReGrow_Params), Pointer  :: branched
  End Type Regrow_Params

  !** Eigenmode perturbation parameters
  Type EMode_Pointer_Set
    Type(EMPerturb_Params), Pointer    :: rigid
  End Type EMode_Pointer_Set

  !** Internal rotation perturbation parameters
  Type IRot_Pointer_Set
    Type(IRotPerturb_Params), Pointer  :: rigid
  End Type IRot_Pointer_Set
  
Contains
  !---------------------------------------------------------------------------
  ! Initializes the transform type move
  ! Requires:  params -- generic transform pointer set
  !            movetype -- string identifying move type
  !            modeltype -- string identifying coordinate system model
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- other auxiliary stuff (ex. cavity bias)
  !---------------------------------------------------------------------------
  Subroutine transform_init(params,movetype,modeltype,filename,spc, &
      species,simcell,auxmv)
    Type(Transform_Moves), Intent(Out)          :: params
    Character(*), Intent(In)                    :: movetype,modeltype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: simcell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    
    Integer      :: error
    Logical      :: nomatch
    
    !** Nullify all the pointer instances in params
    Call transform_nullify(params)

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    Select Case (Trim(ToUpper(movetype)))
    Case (transform_idstrings(1))   !** Regrow move

      Allocate(params%regrow, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"regrow params")

      !** Identify coordinate model
      Select Case (Trim(toupper(modeltype)))
      Case('BRANCHED')
        Allocate(params%regrow%branched, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"branched regrow")
        Call brmoves_RegrowInit(params%regrow%branched,simcell, &
            auxmv, spc,species,filename)
        nomatch = .False.
      End Select

    Case (transform_idstrings(2))   !** Eigenmode perturbation move

      Allocate(params%emode, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"emode params")

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%emode%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"emode rigid")
        Call emperturb_init(params%emode%rigid,spc,species,filename)
        nomatch = .False.
      End Select

    Case (transform_idstrings(3))   !** Internal rotation move

      Allocate(params%irot, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"irot params")

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%irot%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"irot rigid")
        Call irotperturb_init(params%irot%rigid,spc,filename)
        nomatch = .False.
      End Select

    End Select

    !** Give error feedback if movetype and coordinate model couldn't be matched
    If (nomatch) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not get move-parameters for modeltype ", Trim(modeltype)
      Stop
    End If 

  End Subroutine transform_init

  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic transform pointer set
  !---------------------------------------------------------------------------
  Subroutine transform_nullify(params)
    Type(Transform_Moves), Intent(InOut)          :: params

    Nullify(params%regrow)
    Nullify(params%emode)
    Nullify(params%irot)

  End Subroutine transform_nullify

  !-----------------------------------------------------------------------------
  ! Does the transformation based on the movetype and coordinate system
  ! type.  Returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !            subint -- subset interactions (molecule-system)
  !-----------------------------------------------------------------------------
  Logical Function transform_move(params,subset,species,simcell, &
      biasfactor,subint)
    Type(Transform_Moves), Intent(In)               :: params
    Integer, Dimension(:), Intent(In)               :: subset
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Real(kind=RDbl), Intent(Out)                    :: biasfactor
    Type(Subset_Interactions), Intent(InOut)        :: subint

    Integer     :: natoms,atomno,spc,mol
    Logical     :: nomatch,succ_flag

    transform_move = .False.

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%regrow)) Then   !** Regrow move
      spc = subset(1)
      mol = subset(2)

      !** regrow will change the energy of all atoms above atomno, so oldnrgs 
      !** and newnrgs will return Only energy of these atoms. STILL VALID??
      If (Associated(params%regrow%branched)) Then

        !** the chain will be regrown by selecting new positions for 
        !** atoms - 'atomno+1' to 'natoms'
        natoms = config_getnatoms(species,spc)
        atomno = Int(rranf()*(natoms-1)) + 1

        Call brmoves_regrow(species,spc,mol,atomno,params%regrow%branched, &
            biasfactor,succ_flag,subint)
        transform_move = succ_flag
        nomatch = .False.

      End If

    Else If (Associated(params%emode)) Then   !** Eigenmode perturbation move

      If (Associated(params%emode%rigid)) Then
        Call emperturb_move(params%emode%rigid,subset(1),subset(2),species, &
            simcell,biasfactor,succ_flag,subint)
        transform_move = succ_flag
        nomatch = .False.
      End If

    Else If (Associated(params%irot)) Then   !** Internal rotation move

      If (Associated(params%irot%rigid)) Then
        Call irotperturb_move(params%irot%rigid,subset(1),subset(2),species, &
            simcell,biasfactor,succ_flag,subint)
        transform_move = succ_flag
        nomatch = .False.
      End If

    End If

    !** Give error feedback if movetype and coordinate model couldn't be matched
    If (nomatch) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop
    End If 

  End Function transform_move

  !----------------------------------------------------------------------
  ! For setting the maximum allowed energy during an insert-move
  ! Useful for cbgcmc, currently only valid for branchedinsert
  ! Requires:  params -- transform move parameters
  !            max_nrg -- maximum energy
  !----------------------------------------------------------------------
  Subroutine transform_setmaxnrg(params,max_nrg)
    Type(Transform_Moves), Intent(InOut)   :: params
    Real(kind=RDbl), Intent(In)            :: max_nrg

    Logical      :: nomatch

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%regrow)) Then   !** Regrow move

      If (Associated(params%regrow%branched)) Then
        params%regrow%branched%max_nrg = max_nrg  !HACK, reaches in
        nomatch = .False.
      End If

    Else If (Associated(params%emode)) Then   !** Eigenmode perturbation move
      !** nothing yet
      nomatch = .False.
    Else If (Associated(params%irot)) Then    !** Internal rotation move
      !** nothing yet
      nomatch = .False.
    End If

    !** Give error feedback if movetype and coordinate model couldn't be matched
    If (nomatch) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop
    End If 

  End Subroutine transform_setmaxnrg

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Transform parameters
  ! Requires:  params -- transform move parameters
  !----------------------------------------------------------------------
  Function transform_dyninfo(params)
    Character(len=lstrLen)                :: transform_dyninfo
    Type(Transform_Moves), Intent(In)     :: params

    transform_dyninfo = ''

    !** Select information routine based on associated pointers
    If (Associated(params%irot)) Then 
      If (Associated(params%irot%rigid))  &
          transform_dyninfo = irotperturb_dyninfo(params%irot%rigid)
    End If

  End Function transform_dyninfo

  !----------------------------------------------------------------------
  ! Adjusts parameters for Transformation move types.  Not all the move
  ! types support or need adjustment.
  ! Requires:  params -- transform move parameters
  !            ratio -- current acceptance ratio
  !----------------------------------------------------------------------
  Subroutine transform_adjustdelta(params,ratio)
    Type(Transform_Moves), Intent(InOut)     :: params
    Real(kind=RDbl)                          :: ratio

    !** Select adjustment routine based on associated pointers
    If (Associated(params%irot)) Then 
      If (Associated(params%irot%rigid))  &
        Call irotperturb_adjust(params%irot%rigid,ratio)
    End If

  End Subroutine transform_adjustdelta

  !----------------------------------------------------------------------
  ! Displays the transform move parameters
  ! Requires:  params -- transform move parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine transform_display(params,indent,unit)
    Type(Transform_Moves), Intent(In)   :: params
    Integer, Intent(In)                 :: indent,unit

    Logical      :: nomatch

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%regrow)) Then   !** Regrow move

      If (Associated(params%regrow%branched)) Then
        Call brmoves_displayparams(params%regrow%branched,unit)
        nomatch = .False.
      End If

    Else If (Associated(params%emode)) Then   !** Eigenmode perturbation move

      If (Associated(params%emode%rigid)) Then
        Call emperturb_display(params%emode%rigid,indent,unit)
        nomatch = .False.
      End If

    Else If (Associated(params%irot)) Then    !** Internal rotation move

      If (Associated(params%irot%rigid)) Then
        Call irotperturb_display(params%irot%rigid,indent,unit)
        nomatch = .False.
      End If

    End If

    !** Give error feedback if movetype and coordinate model couldn't be matched
    If (nomatch) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop
    End If 

  End Subroutine transform_display

  !----------------------------------------------------------------------
  ! Clean the transform move parameters
  ! Requires:  params -- transform move parameters
  !----------------------------------------------------------------------
  Subroutine transform_clean(params)
    Type(Transform_Moves), Intent(InOut)   :: params

    Logical      :: nomatch

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%regrow)) Then   !** Regrow move
      !** nothing yet
      nomatch = .False.

    Else If (Associated(params%emode)) Then   !** Eigenmode perturbation move

      If (Associated(params%emode%rigid)) Then
        Call emperturb_clean(params%emode%rigid)
        nomatch = .False.
      End If

    Else If (Associated(params%irot)) Then    !** Internal rotation move

      If (Associated(params%irot%rigid)) Then
        Call irotperturb_clean(params%irot%rigid)
        nomatch = .False.
      End If

    End If

    !** Give error feedback if movetype and coordinate model couldn't be matched
    If (nomatch) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop
    End If 

  End Subroutine transform_clean

  !----------------------------------------------------------------------
  ! does adjust ment s for cavity bias after a move is accepted 
  !----------------------------------------------------------------------
  Subroutine transform_postadjust(params, species, sorb, mol )
    Type(Transform_Moves), Intent(In)              :: params
    Type(AtMolCoords), Dimension(:), Intent(In)  :: species    
    Integer, Intent(In)                             :: sorb, mol

    If (Associated(params%regrow)) Then
      If (Associated(params%regrow%branched)) Then
        Call brmoves_postadjust(params%regrow%branched, species, sorb, mol)
      Endif
    End If

  End Subroutine transform_postadjust
  
    
End Module transform
