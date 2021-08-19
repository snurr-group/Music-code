!------------------------------------------------------------------------
! This module acts as a common interface for Identity change moves
!------------------------------------------------------------------------

Module idchange
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only:  filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay
  Use random, Only: rranf
  Use molecules, Only: molecules_getgcmodeltype
  Use config, Only: AtMolCoords, config_getnatoms
  Use simcell, Only: SimCell_Params
  Use subinteract, Only: Subset_Interactions
  Use idflip, Only: IDflip_params, idflip_idstring, idflip_init, idflip_move, &
      idflip_spcs, idflip_restore, idflip_dyninfo, idflip_display, &
      idflip_clean, idflip_postadjust
  Use idflipregrow, Only: IDflipRegrow_params, idflipRegrow_idstring, &
      idflipregrow_init, idflipregrow_move, &
      idflipregrow_spcs, idflipregrow_restore, &
      idflipregrow_display, idflipregrow_clean, idflipregrow_postadjust
  Use dimhop, Only: Expansion_Params, Contraction_Params, &
      dimhop_expand_idstring, dimhop_contract_idstring, dimhop_initexpand, &
      dimhop_initcontract, dimhop_expand, dimhop_contract, &
      dimhop_unexpand, dimhop_uncontract, dimhop_spcs, &
      dimhop_dispexpand, dimhop_dispcontract, dimhop_cleanexpand, &
      dimhop_cleancontract
  Use idlibrary, Only: IDLibrary_Params,idlibrary_init,idlibrary_idstring, &
      idlibrary_spcs,idlibrary_move,idlibrary_restore,idlibrary_display, &
      idlibrary_clean

  Implicit None
  Save

  Private
  Public :: IDchange_Moves,idchange_init,idchange_idstrings, &
      idchange_move,idchange_restore,idchange_spcs,idchange_dyninfo, &
      idchange_display,idchange_clean, idchange_postadjust

  !** Collection of ID change move pointer sets
  Type IDchange_Moves
    Type(IDflip_Params), Pointer          :: idflip
    Type(IDflipReGrow_Params), Pointer    :: idflipregrow
    Type(Expansion_Params), Pointer       :: expand
    Type(Contraction_Params), Pointer     :: contract
    Type(IDLibrary_Params), Pointer       :: idlibrary
  End Type IDchange_Moves

  Character(len=strLen), Dimension(5), Parameter :: idchange_idstrings = &
      (/idflip_idstring,dimhop_expand_idstring,dimhop_contract_idstring, &
      idlibrary_idstring,idflipregrow_idstring/)

Contains
  !---------------------------------------------------------------------------
  ! Initializes the identity change type move
  ! Requires:  params -- generic idchange pointer set
  !            movetype -- string identifying move type
  !            filename -- filename where initialization data can be found
  !            spc -- species number
  !            species -- species data structure
  !            simcell -- the simulation cell information
  !            auxmv -- other auxiliary stuff (ex. cavity bias)
  !---------------------------------------------------------------------------
  Subroutine idchange_init(params,movetype,filename,spc,species,scell,auxmv)
    Type(IDchange_Moves), Intent(Out)           :: params
    Character(*), Intent(In)                    :: movetype,filename
    Integer, Intent(In)                         :: spc   
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: scell
    Type(AuxMoveObjects),Pointer           :: auxmv 

    Integer      :: error

    !** Nullify all the pointer instances in params
    Call idchange_nullify(params)

    !** Call appropriate module/routine based on movetype and coordinate model
    Select Case (Trim(ToUpper(movetype)))
    Case (idflip_idstring)   !** Identity flip move

      Allocate(params%idflip, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"idflip")
      Call idflip_init(params%idflip,spc,auxmv,filename)

    Case (idflipregrow_idstring)   !** Identity flip move (cb)

      Allocate(params%idflipregrow, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"idflip")
      Call idflipregrow_init(params%idflipregrow, scell, auxmv, &
          spc, species, filename)

    Case (dimhop_expand_idstring)   !** Expansion move

      Allocate(params%expand, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"expand")
      Call dimhop_initexpand(params%expand,spc,filename)

    Case (dimhop_contract_idstring)   !** Contraction move

      Allocate(params%contract, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"contract")
      Call dimhop_initcontract(params%contract,spc,filename)

    Case (idlibrary_idstring)   !** Contraction move

      Allocate(params%idlibrary, STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__,"idlibrary")
      Call idlibrary_init(params%idlibrary,spc,filename)

    Case Default
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not get move-parameters for movetype ", Trim(movetype)
      Stop
    End Select

  End Subroutine idchange_init

  !---------------------------------------------------------------------------
  ! Nullify the main pointer set
  ! Requires:  params -- generic idchange pointer set
  !---------------------------------------------------------------------------
  Subroutine idchange_nullify(params)
    Type(IDchange_Moves), Intent(InOut)          :: params

    Nullify(params%idflip)
    Nullify(params%idflipregrow)
    Nullify(params%expand)
    Nullify(params%contract)
    Nullify(params%idlibrary)

  End Subroutine idchange_nullify

  !-----------------------------------------------------------------------------
  ! Does the identity change move based on the movetype and coordinate system
  ! type.  Returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !            opt_spc2 -- type into which the id of subset should be changed
  !-----------------------------------------------------------------------------
  Logical Function idchange_move(params,subset,subints,species,simcell, &
      biasfactor)
    Type(IDchange_Moves), Intent(InOut)                    :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor
    Integer     :: natoms,atomno,spc2
    Logical     :: success

    idchange_move = .False.

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%idflip)) Then   !** ID flip
      idchange_move = idflip_move(params%idflip,subset,subints,species, &
          simcell,biasfactor)

    Else If (Associated(params%idflipRegrow)) Then   !** ID flip
      idchange_move = idflipregrow_move(params%idflipregrow, subset, &
          subints,species, simcell,biasfactor)

    Else If (Associated(params%expand)) Then   !** Expansion
      idchange_move = dimhop_expand(params%expand,subset,subints,species, &
          simcell,biasfactor)

    Else If (Associated(params%contract)) Then   !** Contraction
      idchange_move = dimhop_contract(params%contract,subset,subints,species, &
          simcell,biasfactor)

    Else If (Associated(params%idlibrary)) Then   !** ID library move
      idchange_move = idlibrary_move(params%idlibrary,subset,subints,species, &
          simcell,biasfactor)
    End If

  End Function idchange_move

  !-----------------------------------------------------------------------------
  ! Restores the system to its state before the identity change move
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-----------------------------------------------------------------------------
  Subroutine idchange_restore(params,subset,subints,species,simcell)
    Type(IDchange_Moves), Intent(In)                       :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Integer     :: natoms,atomno,spc,mol
    Logical     :: success

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%idflip)) Then   !** ID flip
      Call idflip_restore(params%idflip,subset,subints,species,simcell)

    Elseif (Associated(params%idflipregrow)) Then   !** ID flip
      Call idflipregrow_restore(params%idflipregrow, subset, subints, &
          species,simcell)

    Else If (Associated(params%expand)) Then   !** Expansion
      Call dimhop_unexpand(params%expand,subset,subints,species,simcell)

    Else If (Associated(params%contract)) Then   !** Contraction
      Call dimhop_uncontract(params%contract,subset,subints,species,simcell)

    Else If (Associated(params%idlibrary)) Then   !** ID library
      Call idlibrary_restore(params%idlibrary,subset,subints,species,simcell)
    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine idchange_restore



  !---------------------------------------------------------------------------
  ! adjusts cavity bias info after a move 
  ! Requires:  params  -- parameters for this move
  !            species -- species data structure
  !            simcell -- simulation cell information
  !--------------------------------------------------------------------------
  Subroutine idchange_postadjust(params,species)
    Type(IDchange_Moves), Intent(InOut)                    :: params
    Type(AtMolCoords), Dimension(:), Intent(In)         :: species

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%idflip)) Then   !** ID flip
      Call idflip_postadjust(params%idflip,species)

    Elseif (Associated(params%idflipregrow)) Then   !** ID flip
      Call idflipregrow_postadjust(params%idflipregrow, species)
    Else
      Return
    End If

  End Subroutine idchange_postadjust



  !-------------------------------------------------------------------------
  ! Get the species number identity of the new molecule created
  ! Requires:  params  -- parameters for this move
  !            spcs -- returned species numbers
  !-------------------------------------------------------------------------
  Subroutine idchange_spcs(params,spcs)
    Type(IDchange_Moves), Intent(In)         :: params
    Integer, Dimension(:), Intent(Out)       :: spcs

    !** Default, zero the returned array
    spcs = 0

    If (Associated(params%idflip)) Then   !** ID flip
      Call idflip_spcs(params%idflip,spcs(1:2))

    Else If (Associated(params%idflipregrow)) Then   !** ID flip regrow
      Call idflipregrow_spcs(params%idflipregrow,spcs(1:2))

    Else If (Associated(params%expand)) Then   !** Expansion
      Call dimhop_spcs(params%expand,spcs)
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'species changes after expand: ',spcs

    Else If (Associated(params%contract)) Then   !** Contraction
      Call dimhop_spcs(params%contract,spcs)
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'species changes after contract: ',spcs

    Else If (Associated(params%idlibrary)) Then   !** ID library move
      Call idlibrary_spcs(params%idlibrary,spcs)
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Write(*,*) 'ID library changes to species: ',spcs

    Else
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    End If

  End Subroutine idchange_spcs

  !----------------------------------------------------------------------
  ! Single line display information for dynamic IDchange parameters
  ! Requires:  params -- swap move parameters
  !----------------------------------------------------------------------
  Function idchange_dyninfo(params)
    Character(len=lstrLen)             :: idchange_dyninfo
    Type(IDchange_Moves), Intent(In)   :: params

    idchange_dyninfo = ''

    !** Select information routine based on associated pointers
    If (Associated(params%idflip)) Then 
      idchange_dyninfo  = idflip_dyninfo(params%idflip)
    Else If (Associated(params%expand)) Then   !** Expansion

    Else If (Associated(params%contract)) Then   !** Contraction

    Else If (Associated(params%contract)) Then   !** ID library move

    End If

  End Function idchange_dyninfo

  !----------------------------------------------------------------------
  ! Displays the ID change move parameters
  ! Requires:  params -- swap move parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine idchange_display(params,indent,unit)
    Type(IDchange_Moves), Intent(In)    :: params
    Integer, Intent(In)                 :: indent,unit

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%idflip)) Then   !** ID flip
      Call idflip_display(params%idflip,indent,unit)

    ElseIf (Associated(params%idflipregrow)) Then   !** ID flip
      Call idflipregrow_display(params%idflipregrow,indent,unit)

    Else If (Associated(params%expand)) Then   !** Expansion
      Call dimhop_dispexpand(params%expand,indent,unit)

    Else If (Associated(params%contract)) Then   !** Contraction
      Call dimhop_dispcontract(params%contract,indent,unit)

    Else If (Associated(params%idlibrary)) Then   !** ID library move
      Call idlibrary_display(params%idlibrary,indent,unit)

    Else
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop
    End If

  End Subroutine idchange_display

  !----------------------------------------------------------------------
  ! Clean the ID change move parameters
  ! Requires:  params -- idchange move parameters
  !----------------------------------------------------------------------
  Subroutine idchange_clean(params)
    Type(IDchange_Moves), Intent(InOut)   :: params

    !** Call appropriate module/routine based on movetype and coordinate model
    If (Associated(params%idflip)) Then   !** Identity swap 
      Call idflip_clean(params%idflip)

    Else If (Associated(params%expand)) Then   !** Expansion
      Call dimhop_cleanexpand(params%expand)

    Else If (Associated(params%contract)) Then   !** Contraction
      Call dimhop_cleancontract(params%contract)

    Else If (Associated(params%idlibrary)) Then   !** ID library move
      Call idlibrary_clean(params%idlibrary)

    Else
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Could not find matching movetype and coordinate model pointer"
      Stop

    End If

  End Subroutine idchange_clean

End Module idchange
