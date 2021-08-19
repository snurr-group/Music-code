!----------------------------------------------------------------------------
! This module acts as a common interface for deletion movements in 
! different coordinate systems
!----------------------------------------------------------------------------

Module delete
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen, TOTAL_INDEX
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, toupper, &
      checkandstop, allocErrDisplay
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use brmoves, Only: BranchedCBDel_Params,BranchedRDel_Params, &
      brmoves_cbdelete, brmoves_rdelete, brmoves_CBDelInit, &
      brmoves_displayparams, brmoves_postadjust
  Use rigidmoves, Only: RigidRDel_Params, RigidBdel_Params, &
      rigidmoves_RDelInit, rigidmoves_BDelInit, rigidmoves_rdelete, &
      rigidmoves_bdelete, rigidmoves_displayparams, rigidmoves_postadjust

  Implicit None
  Save

  Private
  Public :: Delete_Params, RDelete_Params, BDelete_Params, CBDelete_Params, &
      delete_move, delete_init,delete_display, &
      delete_setmaxnrg, delete_postadjust, delete_idstrings, delete_checkIntra

  !***  Main Deletion Parameter Object
  !*** 
  Type Delete_Params
    Type(BDelete_Params), Pointer     :: biased 
    Type(CBDelete_Params), Pointer    :: cbiased
    Type(RDelete_Params), Pointer     :: random
  End Type Delete_Params

  Type BDelete_Params
    Type(RigidBDel_Params), Pointer     :: rigid
  End Type BDelete_Params

  Type CBDelete_Params
    Type(BranchedCBDel_Params), Pointer      :: branched
  End Type CBDelete_Params

  Type RDelete_Params
    Type(RigidRDel_Params), Pointer     :: rigid
    Type(BranchedRDel_Params), Pointer      :: branched
  End Type RDelete_Params

  Character(len=strLen), Dimension(3), Parameter :: delete_idstrings = &
      (/'CBDELETE','BDELETE ', 'RDELETE '/)

Contains


  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic delete pointer set
  !---------------------------------------------------------------------------
  Subroutine delete_nullify(params)
    Type(Delete_Params), Intent(InOut)          :: params

    Nullify(params%biased)
    Nullify(params%random)
    Nullify(params%cbiased)

  End Subroutine delete_nullify


  !---------------------------------------------------------------------------
  ! Initializes the delete type move
  ! Requires:  params -- generic delete object to be inited 
  !            movetype -- string identifying move type
  !            modeltype -- string identifying coordinate system model
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- contains info like cavity bias etc.....
  !---------------------------------------------------------------------------
  Subroutine delete_init(params,movetype,modeltype,filename,spc, &
      species,simcell,auxmv)
    Type(Delete_Params), Intent(Out)             :: params
    Character(*), Intent(In)                    :: movetype,modeltype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: simcell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    
    Integer      :: error
    Logical      :: nomatch
    
    !** Nullify all the pointer instances in params
    Call delete_nullify(params)

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.


    Select Case (Trim(ToUpper(movetype)))


    Case (delete_idstrings(1))   
      !** -------- config-baised delete move -------

      Allocate(params%cbiased, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbiased params")

      !** Identify coordinate model
      Select Case (Trim(toupper(modeltype)))
      Case('BRANCHED')
        Allocate(params%cbiased%branched, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"branched regrow")
        Call brmoves_CBDelInit(params%cbiased%branched,simcell, auxmv, &
            spc,species,filename)
        nomatch = .False.
      End Select



    Case (delete_idstrings(2))   
      !** ------- Biased Deletion -------- 

      Allocate(params%biased, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"biased params")

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%biased%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"biased rigid")
        Call rigidmoves_BDelInit(params%biased%rigid, simcell, spc, filename)
        nomatch = .False.
      End Select



    Case (delete_idstrings(3))   
      !** ----- Random unbiased deleteions --------

      Allocate(params%random, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random params")
      
      !** Nullify all the pointer instances in DelParams
      Nullify(params%random%rigid)
      Nullify(params%random%branched)

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%random%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Call rigidmoves_RDelInit(params%random%rigid,simcell,spc)
        nomatch = .False.        
      Case Default
        Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
            " Could not get move-parameters for modeltype ", Trim(modeltype)
        Stop
      End Select



    End Select

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get move-parameters for modeltype "//Trim(modeltype) )

  End Subroutine delete_init


  !----------------------------------------------------------------------
  ! gets the bias factor, and energies for deletion of one molecule 
  ! of the  given sorb at given molec number
  ! Requires:  params  -- parameters for this delete-move
  !            subset -- system subset to which to apply move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias used during insertion of this molecule
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Logical Function delete_move(params,subset,species,simcell, &
      biasfactor,subint)
    Type(Delete_Params), Intent(In)               :: params
    Integer, Dimension(:), Intent(In)               :: subset
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Real(kind=RDbl), Intent(Out)                    :: biasfactor
    Type(Subset_Interactions), Intent(InOut)        :: subint

    Integer     :: natoms,atomno,spc,mol
    Logical     :: nomatch,nrg_calc_flag

    delete_move = .False.
    nomatch = .True.
    spc = subset(1)
    mol = subset(2)

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%biased)) Then   
      !** Biased Deletion
      Call rigidmoves_bdelete(params%biased%rigid, species,spc,mol, &
          biasfactor,nrg_calc_flag,subint)
      delete_move = nrg_calc_flag
      nomatch = .False.

    Else If (Associated(params%random)) Then   
      !** Unbiased deletion move
      If (Associated(params%random%rigid)) Then
        Call rigidmoves_rdelete(params%random%rigid, species, spc, mol, &
            biasfactor, nrg_calc_flag,subint)
        delete_move = nrg_calc_flag
        nomatch = .False.
      End If

    Else If (Associated(params%cbiased)) Then   
      !** Configurational bised deletion move
      Call brmoves_cbdelete(species,spc,mol,params%cbiased%branched,&
          biasfactor,subint,nrg_calc_flag)
      delete_move = nrg_calc_flag
      nomatch = .False.
    End If

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not find appropriate move-parameters or modeltype ")

  End Function delete_move


  !----------------------------------------------------------------------
  ! checks whether intra should be used for acceptance
  !----------------------------------------------------------------------
  Subroutine delete_checkIntra(params,use_intra)
    Type(Delete_Params), Intent(in)  :: params
    Logical, Intent(inout) :: use_intra
    
    If (Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) Then
        use_intra=.True.
      Endif
    Endif
  End Subroutine delete_checkIntra

  !----------------------------------------------------------------------
  ! Does some post acceptance move for cavity bias etc
  ! Requires:  params -- delete move parameters
  !            species, spc, molec
  !----------------------------------------------------------------------
  Subroutine delete_postadjust(params, species, sorb, molec )
    Type(Delete_Params), Intent(In)   :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species    
    Integer, Intent(In)                 :: sorb, molec

    If (Associated(params%biased)) Then
      Call rigidmoves_postadjust(params%biased%rigid,&
          species, sorb, molec)
    ELseif (Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) &       
          Call brmoves_postadjust(params%cbiased%branched,species,sorb,molec)
    Endif

  End Subroutine delete_postadjust

  !----------------------------------------------------------------------
  ! Displays the delete move parameters
  ! Requires:  params -- delete move parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine delete_display(params,indent,unit)
    Type(Delete_Params), Intent(In)   :: params
    Integer, Intent(In)                 :: indent,unit

    Logical      :: nomatch

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%cbiased)) Then   

      If (Associated(params%cbiased%branched)) Then
        Call brmoves_displayparams(params%cbiased%branched,unit)
        nomatch = .False.
      End If

    Else If (Associated(params%biased)) Then   

      If (Associated(params%biased%rigid)) Then
        Call rigidmoves_displayparams(params%biased%rigid,indent,unit)
        nomatch = .False.
      End If

    Else If (Associated(params%random)) Then    !** Internal rotation move

      If (Associated(params%random%rigid)) Then
!!$        Call rigidmoves_displayparams(params%random%rigid,indent,unit)
        nomatch = .False.
      End If

    End If
    
    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get appropriate move-parameters or modeltype ")

  End Subroutine delete_display

  !----------------------------------------------------------------------
  ! for setting the maximum alllowed energy during an insert-move
  ! we need it for calculating the bias of insertion used
  ! Useful for cbgcmc, currently only valid for brancheddelete 
  !----------------------------------------------------------------------
  Subroutine delete_setmaxnrg(params,max_nrg)
    Type(Delete_Params), Intent(inout)  :: params
    Real(kind=RDbl), Intent(in) :: max_nrg
    If (Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) Then
        params%cbiased%branched%max_nrg=max_nrg
      Endif
    Endif
  End Subroutine delete_setmaxnrg
  
End Module delete
