!-------------------------------------------------------------------------
! This module acts as a common interface for insertion movements in 
! different coordinate systems
!-------------------------------------------------------------------------

Module insert
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen, zero, TOTAL_INDEX
  Use config, Only: AtMolCoords
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use rigidmoves, Only: RigidBIns_Params,RigidRIns_Params, &
      RigidLIns_Params, rigidmoves_RInsInit, &
      rigidmoves_BInsInit, rigidmoves_LInsInit, &
      rigidmoves_rinsert, rigidmoves_binsert, rigidmoves_linsert, &
      rigidmoves_displayparams, rigidmoves_postadjust
  Use brmoves, Only : BranchedCBIns_Params, &     
      BranchedRIns_Params, brmoves_CBInsinit, brmoves_cbinsert, &
      brmoves_rinsert , brmoves_displayparams, brmoves_postadjust
  Use utils, Only: isfileopen, filesrchstr, stripcmnt, split, &
      toint, toupper, allocErrDisplay, checkandstop

  Implicit None
  Save

  Private 
  Public :: Insert_Params, BInsert_Params, CBInsert_Params, RInsert_Params, &
      LInsert_Params, insert_init, insert_move, insert_display, &
      insert_setmaxnrg, insert_postadjust, insert_idstrings, insert_checkIntra 

  !***  Main Insertion Parameter Object
  !*** 
  Type Insert_Params
    Type(RInsert_Params), Pointer     :: random
    Type(BInsert_Params), Pointer     :: biased 
    Type(LInsert_Params), Pointer     :: lib
    Type(CBInsert_Params), Pointer    :: cbiased
  End Type Insert_Params
  
  !** Random insertion params
  Type RInsert_Params
    Type(RigidRIns_Params), Pointer            :: rigid
    Type(BranchedRIns_Params), Pointer             :: branched
  End Type RInsert_Params

  !** Biased insertion params
  Type BInsert_Params
    Type(RigidBIns_Params), Pointer            :: rigid
  End Type BInsert_Params
  

  !** Biased insertion params for hybrid-gcmc, Using libraries
  Type LInsert_Params
    Type(RigidLIns_Params), Pointer            :: rigid
  End Type LInsert_Params

  !** Biased insertion params
  Type CBInsert_Params
    Type(BranchedCBIns_Params), Pointer             :: branched
  End Type CBInsert_Params
  

  Character(len=strLen), Dimension(4), Parameter :: insert_idstrings = &
      (/'RINSERT ', 'BINSERT ',  'LINSERT ', 'CBINSERT' /)

Contains


  !---------------------------------------------------------------------------
  ! Initializes the Insert type move
  ! Requires:  params -- generic insert object to be inited 
  !            movetype -- string identifying move type
  !            modeltype -- string identifying coordinate system model
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- contains info like cavity bias etc.....
  !---------------------------------------------------------------------------
  Subroutine insert_init(params,movetype,modeltype,filename,spc, &
      species,simcell,auxmv)
    Type(Insert_Params), Intent(Out)            :: params
    Character(*), Intent(In)                    :: movetype,modeltype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: simcell
    Type(AuxMoveObjects),Pointer           :: auxmv 
    
    Integer      :: error
    Logical      :: nomatch
    
    !** Nullify all the pointer instances in params
    Call insert_nullify(params)
    nomatch = .True.


    !** Call appropriate module/routine based on movetype and coordinate model
    Select Case (Trim(ToUpper(movetype)))


    Case (insert_idstrings(1))   
      !** -------- Random Insertion move ---------

      Allocate(params%random, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random params")

      Nullify(params%random%branched)

      !** Identify coordinate model
      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%random%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"random rigid")
        Call rigidmoves_RInsInit(params%random%rigid, simcell, auxmv, spc)
        nomatch = .False.
      End Select



    Case (insert_idstrings(2))   
      !** ------- Biased Insertion -------- 

      Allocate(params%biased, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"biased params")

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%biased%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"biased rigid")
        Call rigidmoves_BInsInit(params%biased%rigid, simcell, auxmv, spc, &
          filename)
        nomatch = .False.
      End Select


    Case (insert_idstrings(3))   
      !** ----- Insertion from library --------

      Allocate(params%lib, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"library ")
      
      !** Nullify all the pointer instances in DelParams
      Nullify(params%lib%rigid)

      Select Case (Trim(toupper(modeltype)))
      Case('RIGID')
        Allocate(params%lib%rigid, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
        Call rigidmoves_LInsInit(params%lib%rigid, simcell, auxmv, &
            spc, filename)
        nomatch = .False.        
      End Select



    Case (insert_idstrings(4))   
      !** ----- Config biased Insertions --------
      Allocate(params%cbiased, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"cbiased params")

      !** Identify coordinate model
      Select Case (Trim(toupper(modeltype)))
      Case('BRANCHED')
        Allocate(params%cbiased%branched, STAT=error)
        If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"branched insert")
        Call brmoves_CBInsInit(params%cbiased%branched,simcell, auxmv, &
            spc,species,filename)
        nomatch = .False.
      End Select



    End Select

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get move-parameters for modeltype "//Trim(modeltype) )

  End Subroutine insert_init


  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic insert pointer set
  !---------------------------------------------------------------------------
  Subroutine insert_nullify(params)
    Type(Insert_Params), Intent(InOut)          :: params

    Nullify(params%random)
    Nullify(params%biased)
    Nullify(params%lib)
    Nullify(params%cbiased)

  End Subroutine insert_nullify


  !----------------------------------------------------------------------
  ! gets the bias factor, and energies for the insertion of one molecule 
  ! of the  given sorb at given molec number
  ! Requires:  params  -- parameters for this insert-move
  !            subset -- system subset to which to apply move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias used during insertion of this molecule
  !            subint -- subset interactions (molecule-system)
  !---------------------------------------------------------------------------
  Logical Function insert_move(params,subset,species,simcell, &
      biasfactor,subint)
    Type(Insert_Params), Intent(Inout)               :: params
    Integer, Dimension(:), Intent(In)               :: subset
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species
    Type(SimCell_Params), Intent(In)                :: simcell
    Real(kind=RDbl), Intent(Out)                    :: biasfactor
    Type(Subset_Interactions), Intent(InOut)        :: subint

    Integer     :: natoms,atomno,spc,mol
    Logical     :: nomatch,nrg_calc_flag

    insert_move = .False.
    nomatch = .True.
    spc = subset(1)
    mol = subset(2)

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%random)) Then   
      !** Unbiased Insertion move
      If (Associated(params%random%rigid)) Then
      Call rigidmoves_rinsert(params%random%rigid, species, spc, mol, &
          biasfactor, nrg_calc_flag,subint)
        insert_move = nrg_calc_flag
        nomatch = .False.
      End If

    Else If (Associated(params%biased)) Then   

      !** Biased Insertion
      Call rigidmoves_binsert(params%biased%rigid, species, spc, mol, &
          biasfactor, nrg_calc_flag,subint)
      insert_move = nrg_calc_flag
      nomatch = .False.

    Else If (Associated(params%lib)) Then
      !** Librray insertion move 
      Call rigidmoves_linsert(params%lib%rigid, species, spc, mol, &
          biasfactor, nrg_calc_flag,subint)
      insert_move = nrg_calc_flag
      nomatch = .False.

    Else If (Associated(params%cbiased)) Then
      !** Configurational biased Insertion move
      Call brmoves_cbinsert(species,spc,mol,params%cbiased%branched,&
          biasfactor,nrg_calc_flag, subint)

      insert_move = nrg_calc_flag
      nomatch = .False.
    End If

    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not find appropriate move-parameters or modeltype ")

  End Function insert_move



  !----------------------------------------------------------------------
  ! for setting the maximum alllowed energy during an insert-move
  ! Useful for cbgcmc, currently only valid for branchedinsert
  !----------------------------------------------------------------------
  Subroutine insert_setmaxnrg(params,max_nrg)
    Type(Insert_Params), Intent(inout)  :: params
    Real(kind=RDbl), Intent(in) :: max_nrg
    
    If (Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) Then
        params%cbiased%branched%max_nrg=max_nrg
      Else
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Write(*,*) "Movetype or Coord type does not match this routine"
      stop
      
    Endif
  End Subroutine insert_setmaxnrg


  !----------------------------------------------------------------------
  ! checks whether intra should be used for acceptance
  !----------------------------------------------------------------------
  Subroutine insert_checkIntra(params,use_intra)
    Type(Insert_Params), Intent(in)  :: params
    Logical, Intent(inout) :: use_intra
    
    If (Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) Then
        use_intra=.True.
      Endif
    Endif
  End Subroutine insert_checkIntra


  !----------------------------------------------------------------------
  ! Does some post acceptance move for cavity bias etc
  ! Requires:  params -- insert move parameters
  !            species, spc, molec
  !----------------------------------------------------------------------
  Subroutine insert_postadjust(params, species, sorb, molec )
    Type(Insert_Params), Intent(In)   :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)  :: species    
    Integer, Intent(In)                 :: sorb, molec

    If (Associated(params%biased)) Then
      If (Associated(params%biased%rigid)) &
          Call rigidmoves_postadjust(params%biased%rigid, species, sorb,molec)
    Elseif(Associated(params%lib)) Then
      If (Associated(params%lib%rigid)) &
          Call rigidmoves_postadjust(params%lib%rigid, species, sorb,molec) 
    Elseif(Associated(params%cbiased)) Then
      If (Associated(params%cbiased%branched)) &       
          Call brmoves_postadjust(params%cbiased%branched,species,sorb,molec) 
    Endif

  End Subroutine insert_postadjust



  !----------------------------------------------------------------------
  ! Displays the insert move parameters
  ! Requires:  params -- insert move parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine insert_display(params,indent,unit)
    Type(Insert_Params), Intent(In)     :: params
    Integer, Intent(In)                 :: indent,unit

    Logical      :: nomatch

    !** Call appropriate module/routine based on movetype and coordinate model
    nomatch = .True.
    If (Associated(params%random)) Then    
      If (Associated(params%random%rigid)) Then
!!$        Call rigidmoves_displayparams(params%random%rigid,indent,unit)
        nomatch = .False.
      End If
    Else If (Associated(params%biased)) Then   
      
      If (Associated(params%biased%rigid)) Then
        Call rigidmoves_displayparams(params%biased%rigid,indent,unit)
        nomatch = .False.
      End If
    Else If (Associated(params%lib)) Then   
      If (Associated(params%lib%rigid)) Then
        Call rigidmoves_displayparams(params%lib%rigid,indent,unit)
        nomatch = .False.
      End If
    Else If (Associated(params%cbiased)) Then   
      If (Associated(params%cbiased%branched)) Then
        Call brmoves_displayparams(params%cbiased%branched,unit)
        nomatch = .False.
      End If
    End If
    
    !** Give error feedback if movetype and coord model couldn't be matched
    Call checkandstop(.Not.nomatch, __FILE__, __LINE__, &
        " Could not get appropriate move-parameters or modeltype ")

  End Subroutine insert_display
  
End Module insert
