!------------------------------------------------------------------------
! This module acts as a common interface for translation movements in 
! different coordinate systems.
!
! At the moment, this routine handles two "translation" move types:
! 1) random translation of molecule by a randomly generated vector
! 2) reinsertion of molecule in new position, with new orientation
!
! Note that the "reinsertion" move is best handled here because it
! involves translation and can substitute for a translation move. 
! It could also substitute for a rotation move, but in the simpliest
! uniatomic case we need only a translation move.
!------------------------------------------------------------------------
Module translate
  Use auxmoveparams, Only: AuxMoveObjects
  Use brmoves, Only: BranchedFTrans_Params,BranchedRTRans_Params, & 
      brmoves_rtransinit, brmoves_ftransinit, brmoves_randomtranslate, &
      brmoves_adjustdeltatrans , brmoves_displayparams, brmoves_postadjust
  Use config, Only: AtMolCoords
  Use defaults, Only: RDbl, strLen, lstrLen, kjmole_kb, scaleke, one
  Use file, Only : file_getunit
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use molecules, Only: molecules_getgcmodeltype
  Use rigidmoves, Only: RigidRTrans_Params,RigidReinsert_Params, &
      rigidmoves_RTransInit,rigidmoves_randomtranslate,rigidmoves_reinsert, &
      rigidmoves_adjustdeltatrans, rigidmoves_postadjust, &
      rigidmoves_displayparams,rigidmoves_dyninfo,rigidmoves_ReinsertInit
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop
  Use vector, Only: VecType   

  Implicit None
  Save

  Private
  Public :: Translate_Params, RTranslate_Params, &
      ReInsert_Pset, translate_move, translate_display, translate_init, &
      translate_adjustdelta, translate_dyninfo, translate_idstrings, &
      translate_postadjust, translate_setmaxnrg

  !** Main Translation Parameter Object
  Type Translate_Params
    Type(RTranslate_Params), Pointer      :: random 
    Type(Reinsert_Pset), Pointer          :: reins 
  End Type Translate_Params
  
  !** Random translate params
  Type RTranslate_Params
    Type(RigidRTrans_Params), Pointer     :: rigid
    Type(BranchedRTRans_Params), Pointer  :: branched
  End Type RTranslate_Params

  !** Reinsertion params
  Type Reinsert_Pset
    Type(RigidReinsert_Params), Pointer :: rigid
  End Type Reinsert_Pset

  Character(len=strLen), Dimension(2), Parameter :: translate_idstrings = &
      (/'RTRANSLATE', 'REINSERT  '/)
  
Contains

  !---------------------------------------------------------------------------
  ! Initializes the Translate type move
  ! Requires:  params -- generic translate object to be inited 
  !            movetype -- string identifying move type
  !            modeltype -- string identifying coordinate system model
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- other auxiliary stuff (ex. cavity bias)
  !---------------------------------------------------------------------------
  Subroutine translate_init(params,movetype,modeltype,filename,spc, &
      species,simcell,auxmv)
    Type(Translate_Params), Intent(Out)            :: params
    Character(*), Intent(In)                    :: movetype,modeltype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: simcell
     Type(AuxMoveObjects),Pointer           :: auxmv    
    Integer      :: error
    Logical      :: nomatch
    
    !** Nullify all the pointer instances in params
    Call translate_nullify(params)
    nomatch = .True.


    !** Call appropriate module/routine based on movetype and coordinate model
    Select Case (Trim(ToUpper(movetype)))
    Case (translate_idstrings(1))   
      !** ------- Translation by  a random amount -------- 

      Allocate(params%random, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random params")
      Nullify(params%random%rigid)
      Nullify(params%random%branched)

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%random%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random rigid")
        Call rigidmoves_RTransInit(params%random%rigid, simcell, auxmv, &
            spc, filename)
        nomatch = .False.

      Case('BRANCHED')
        Allocate(params%random%branched, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random branched")
        Call brmoves_RTransInit(params%random%branched, simcell, auxmv, &
            spc, filename)
        nomatch = .False.

      End Select

    Case (translate_idstrings(2))   
      !** ----- params for Re-Translateing molecules --------

      Allocate(params%reins, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"reinsert")

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%reins%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Call rigidmoves_ReInsertInit(params%reins%rigid,spc,filename)
        nomatch = .False.        
      End Select

    End Select

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get move-parameters for modeltype "//Trim(modeltype) )

  End Subroutine translate_init

  !----------------------------------------------------------------------
  ! Set the maximum allowed energy
  ! Requires:  params -- generic translate object 
  !            max_nrg -- maximum energy (in kcal/mol)
  !----------------------------------------------------------------------
  Subroutine translate_setmaxnrg(params,max_nrg)
    Type(Translate_Params), Intent(InOut)   :: params
    Real(kind=RDbl), Intent(In)             :: max_nrg
    
    If (Associated(params%random)) Then
      If (Associated(params%random%rigid)) Then
        params%random%rigid%max_nrg = max_nrg
      Else If (Associated(params%random%branched)) Then
        params%random%branched%max_nrg = max_nrg
      Else
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' Unable to recognized sub-initialized pointer'
        Stop
      End If

    Else If (Associated(params%reins)) Then
      If (Associated(params%reins%rigid)) Then
        params%reins%rigid%max_nrg = max_nrg
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

  End Subroutine translate_setmaxnrg

  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic translate pointer set
  !---------------------------------------------------------------------------
  Subroutine translate_nullify(params)
    Type(Translate_Params), Intent(InOut)          :: params

    Nullify(params%random)
    Nullify(params%reins)

  End Subroutine translate_nullify

  !----------------------------------------------------------------------
  ! translates a molecule and gets the  energies 
  ! Requires:  params  -- parameters for this translate-move
  !            subset -- system subset to which to apply move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias used during translateion of this molecule
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Logical Function translate_move(params,subset,species,simcell, &
      biasfactor,subint)
    Type(Translate_Params), Intent(InOut)           :: params
    Integer, Dimension(:), Intent(In)               :: subset
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Real(kind=RDbl), Intent(Out)                    :: biasfactor
    Type(Subset_Interactions), Intent(InOut)        :: subint

    Integer     :: natoms,atomno,spc,mol
    Logical     :: nomatch,nrg_calc_flag

    translate_move = .False.
    nomatch = .True.
    spc = subset(1)
    mol = subset(2)

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%random)) Then   
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_randomtranslate(params%random%rigid,species,spc,mol,&
            simcell,nrg_calc_flag, subint)
        biasfactor=one
        translate_move = nrg_calc_flag
        nomatch = .False.

      Else If (Associated(params%random%branched)) Then
        Call brmoves_randomtranslate(params%random%branched,species,spc,mol, &
            translate_move,subint)
        biasfactor = one
        nomatch = .False.
      End If

    Else If (Associated(params%reins)) Then
      Call rigidmoves_reinsert(params%reins%rigid,species,spc,mol,simcell, &
          biasfactor,translate_move,subint)
      nomatch = .False.

    End If

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not find appropriate move-parameters or modeltype ")

  End Function translate_move
  
  !----------------------------------------------------------------------
  ! Displays parameters for translations
  !----------------------------------------------------------------------
  Subroutine translate_display(params, indent, unitno)
    Type(Translate_Params), Intent(in)              :: params
    Integer, Intent(in)                             :: indent, unitno
    
    Logical :: nomatch
    nomatch = .True.

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%random)) Then    
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_displayparams(params%random%rigid,indent,unitno)
        nomatch = .False.
      Else If(Associated(params%random%branched)) Then
        Call brmoves_displayparams(params%random%branched,indent,unitno)
        nomatch = .False.
      End If

    Else If (Associated(params%reins)) Then   
      If (Associated(params%reins%rigid)) Then
        Call rigidmoves_displayparams(params%reins%rigid,indent,unitno)
        nomatch = .False.
      End If
    End If
    
    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get appropriate move-parameters or modeltype ")

  End Subroutine translate_display

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires:  params -- translation parameters
  !----------------------------------------------------------------------
  Function translate_dyninfo(params)
    Character(len=lstrLen)                 :: translate_dyninfo
    Type(Translate_Params), Intent(In)     :: params

    translate_dyninfo  = Repeat(' ',strLen)
    If (Associated(params%random))  Then 
      If (Associated(params%random%rigid))  &
          translate_dyninfo  = rigidmoves_dyninfo(params%random%rigid)
    End If

  End Function translate_dyninfo

  !----------------------------------------------------------------------
  ! adjusts parameters for Random translations
  !----------------------------------------------------------------------
  Subroutine translate_adjustdelta(params,ratio)
    Type(Translate_Params), Intent(inout)           :: params
    Real(kind=RDbl)                                  :: ratio
    If (Associated(params%random))  Then 
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_adjustdeltatrans(params%random%rigid,ratio)
      Else If (Associated(params%random%branched))  Then 
        Call brmoves_adjustdeltatrans(params%random%branched,ratio)
      Endif
    Endif
  End Subroutine translate_adjustdelta


  !----------------------------------------------------------------------
  ! does adjust ment s for cavity bias after a move is accepted 
  !----------------------------------------------------------------------
  Subroutine translate_postadjust(params, species, sorb, mol )
    Type(Translate_Params), Intent(In)              :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species    
    Integer, Intent(In)                             :: sorb, mol

    Logical :: do_adjust

    If (Associated(params%random)) Then
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_postadjust(params%random%rigid, species, sorb, mol)
      ELseif (Associated(params%random%branched)) Then
        Call brmoves_postadjust(params%random%branched, species, sorb, mol)
      Endif

    Else If (Associated(params%reins)) Then

    End If

  End Subroutine translate_postadjust
  
End Module translate
