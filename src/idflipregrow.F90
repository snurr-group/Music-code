!------------------------------------------------------------------------
! This module handles the identity flip move type in branched coords
! [There is another movetype idflip]
!------------------------------------------------------------------------

Module idflipregrow
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop, int2str, cleanstring
  Use file, Only: file_open
  Use random, Only: rranf
  Use molecules, Only: molecules_getgcmodeltype

  ! be vary careful, do not use cavity bias from this module, 
  ! Only postadjust should be done from here -shaji
  Use cavitylist, Only: Cavity_Params, &
      cavitylist_updatemolecinfo, cavitylist_checkandIncrMemory

  Use config, Only: AtMolCoords, config_getnmoles, config_checkandincr, &
      config_setnmoles, config_clean, config_copymolec, config_delmol, &
      config_allocfields, config_idrebuild, config_getnatoms, &
      config_subset2xyz, config_isfixed
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use subinteract, Only: Subset_Interactions, subinteract_changenmoles, &
      subinteract_int, subinteract_newnrg
  Use brmoves, Only: branchedidflip_params, brmoves_idflip, &
      brmoves_displayparams, brmoves_idflipinit
  Implicit None
  Save

  Private
  Public :: Idflipregrow_Params, idflipregrow_init, &
      idflipregrow_idstring, idflipregrow_move, idflipregrow_restore, &
      idflipregrow_display, idflipregrow_clean, idflipregrow_spcs, &
      idflipregrow_postadjust

  !** Simple Identity Flip parameters
  Type Idflipregrow_Params
    Integer                          :: spc     !** species to operate on
    Integer                          :: ntospc
    Integer, Dimension(3)            :: osubset !** old subset (spc,mol,0)
    Integer, Dimension(3)            :: nsubset !** new subset (spc,mol,0)
    Real(kind=RDbl)                  :: max_nrg
    Integer, Dimension(:), Pointer   :: tospc
    Type(AtMolCoords)                :: molcopy !** storage for one molecule
    Type(branchedidflip_Params)      :: brmove
    Type(Cavity_Params),Pointer      :: cavity
  End Type Idflipregrow_Params

  Character(len=strLen), Parameter   :: idflipregrow_idstring = 'IDFLIPREGROW'
  Real(kind=RDbl), Parameter         :: VERY_LRG_NRG = 1.0e10_RDbl

Contains
  !---------------------------------------------------------------------------
  ! Initializes the identity flip type move
  ! Requires:  params -- identity flip parameters
  !            spc -- species number
  !            auxmv -- cavitylist etc.
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine idflipregrow_init(params,scell,auxmv,spc,species,filename)
    Type(Idflipregrow_Params), Intent(Out)            :: params
    Integer, Intent(In)                         :: spc   
    Character(*), Intent(In)                    :: filename
    Type(AtMolCoords), Dimension(:), Intent(In) :: species 
    Type(SimCell_Params), Intent(In)            :: scell
    Type(AuxMoveObjects),Pointer           :: auxmv

    Integer                              :: error,unit,i, newspc
    Character(len=strLen)                   :: line     
    Character(len=strLen), Dimension(20)    :: fields

    !** Initialize the cavity bais pointer
    Nullify(params%cavity)
    If (auxmv%cavBiasON)     params%cavity=>auxmv%cavity
    Write(*,*) auxmv%cavBiasON
    Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__

    !** Defaults
    params%max_nrg = VERY_LRG_NRG
    params%osubset = (/0,0,0/)
    params%nsubset = (/0,0,0/)

    !** Set the species to operate on
    params%spc = spc
    !** as of now it works only for branched coords , put a check here
    !** IMPORTANT : if this check is removed then an equivalent check 
    !** should be done somewhere else
    If ( Trim(Toupper(molecules_getgcmodeltype(spc)))/='BRANCHED') Then
      Write(*,*) molecules_getgcmodeltype(spc), 'BRANCHED'
      Write(*,*) " ALL IDFLIP REGROW SHOULD BE BETWEEN N-ALKANES"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    !** Open the control file 
    unit = file_open(filename)

    !** Read the line containing destination species numbers
    Read(unit,'(a)') line
    line = cleanstring(line)

    !** Place the destination species numbers in a list
    params%ntospc = split(line,fields,',')
    Allocate(params%tospc(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,params%ntospc
      newspc=toint(fields(i))
      params%tospc(i) = newspc
      ! see comment above about this check, IMPORTANT
      If ( Trim(Toupper(molecules_getgcmodeltype(newspc)))/='BRANCHED') Then
        Write(*,*) " ALL IDFLIP REGROW SHOULD BE BETWEEN N-ALKANES"
        Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
        Stop
      Endif
    End Do

    !** Allocate space for storage of one molecule of this species
    Call config_allocfields(params%molcopy,params%spc,1,.True.)
    Call brmoves_idflipInit(params%brmove,scell,auxmv,spc,species,filename)

  End Subroutine idflipregrow_init


  !--------------------------------------------------------------------------
  ! Removes old molecule from cavitylist, and add new molecule 
  ! Requires:  params  -- parameters for this move
  !            species -- species data structure
  !--------------------------------------------------------------------------
  Subroutine idflipregrow_postadjust(params,species)
    Type(IDFlipReGrow_Params), Intent(InOut)            :: params
    Type(AtMolCoords), Dimension(:), Intent(In)         :: species

    Integer     :: spc1, spc2, mol1, mol2
    Logical     :: succ_flag


    If (.Not.Associated(params%cavity)) Return

    ! old molecule - remove from cavity just like regular delete
    spc1=params%osubset(1)
    mol1=params%osubset(2)
    Call cavitylist_updateMolecInfo(params%cavity, species(spc1), &
        spc1, mol1, 110 )

    ! new molecule - add to cavity , just like regular insert
    spc2=params%nsubset(1)
    mol2=params%nsubset(2)

    !** Also check whether memory in cubelist is alright
    !** has to be done before updatemolecinfo
    Call cavitylist_checkandIncrMemory(params%cavity, species(spc2), spc2 )
    Call cavitylist_updateMolecInfo(params%cavity,species(spc2), spc2, &
        mol2, 101 )

  End Subroutine idflipregrow_postadjust


  !------------------------------------------------------------------------
  ! Get the species number identities of the molecules changed.  The output
  ! follows the convention in subinteract_update, species numbers are 
  ! listed as positive if they are inserted, negative if they are deleted
  ! and as both if they are perturbed.
  ! Requires:  params  -- parameters for this move
  !            spcs  -- returned, changed species numbers
  !------------------------------------------------------------------------
  Subroutine idflipregrow_spcs(params,spcs)
    Type(IDFlipRegrow_Params), Intent(In)     :: params
    Integer, Dimension(:), Intent(Out)        :: spcs

    spcs(1) = -params%osubset(1)
    spcs(2) = params%nsubset(1)

  End Subroutine idflipregrow_spcs


  !-------------------------------------------------------------------------
  ! Does the flip move based on the movetype and coordinate system
  ! type.  Returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply identity flip
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-------------------------------------------------------------------------
  Logical Function idflipregrow_move(params,subset,subints,species,simcell, &
      biasfactor)
    Type(Idflipregrow_Params), Intent(InOut)                     :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Integer     :: nmoles2, spc1,mol1,spc2,mol2,maxmoves
    Logical     :: succ_flag,nrg_calc_flag,skip_intra_flag,fast

    !** Set useful values and defaults
    params%osubset = subset
    idflipregrow_move = .True.
    skip_intra_flag = .False.
    fast = .True.  !** hard-coded HACK
    spc1=subset(1)
    mol1=subset(2)

    !** Copy the specified molecule's configuration
    Call config_copymolec(species(spc1),mol1,params%molcopy,1)

    !** Evaluate the energy of the removed molecule, before deleting
    succ_flag = subinteract_int(subints(spc1),species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),subset)

    Call checkandstop(succ_flag,__FILE__,__LINE__, &
        'Problem with nrg calculation of an existing structure')

    !** Get the new species and molecule number 
    spc2 = params%tospc(Int(rranf()*params%ntospc) + 1)

    If ((spc1==spc2).Or.(config_isfixed(species(spc2)))) Then
      Write(*,*) "trying to change spc1 to spc2: ", spc1, spc2
      Write(*,*) "check ctrl file for idflip sorbs specification"
      Write(*,'(1x,2a,i4)') __FILE__," : ",__LINE__
      Stop
    Endif

    nmoles2 = config_getnmoles(species,spc2) + 1
    mol2 = nmoles2
    params%nsubset = (/spc2,mol2,0/)

    Call brmoves_idflip(params%brmove, species, spc1, mol1, &
        spc2, mol2, biasfactor, nrg_calc_flag, subints(spc2))
    idflipregrow_move = nrg_calc_flag

  End Function idflipregrow_move

  !--------------------------------------------------------------------------
  ! Restores the system to its state before the flip move
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !-----------------------------------------------------------------------------
  Subroutine idflipregrow_restore(params,subset,subints,species,simcell)
    Type(Idflipregrow_Params), Intent(In)                        :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Integer     :: nmoles1,nmoles2
    Logical     :: succ_flag

    nmoles1 = config_getnmoles(species,params%osubset(1))
    nmoles2 = config_getnmoles(species,params%nsubset(1))

    !** Remove the new molecule from the new species list
    Call config_setnmoles(species,params%nsubset(1),nmoles2-1)
    succ_flag = subinteract_changenmoles(subints(params%nsubset(1)), &
        params%nsubset(1),nmoles2-1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Copy the molecule in the previous old molecule position 
    !** back to the end of the old species list
    Call config_copymolec(species(params%osubset(1)),params%osubset(2), &
        species(params%osubset(1)),nmoles1+1)
    Call config_setnmoles(species,params%osubset(1),nmoles1+1)

    !** Put the removed old molecule back into its old position
    Call config_copymolec(params%molcopy,1,species(params%osubset(1)), &
        params%osubset(2))

  End Subroutine idflipregrow_restore

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information 
  ! Requires:  unitno
  !----------------------------------------------------------------------------
  Subroutine idflipregrow_sampleCF(unitno)
    Integer, Intent(In)           :: unitno

    Write(unitno,'(a,t30,a)') 'IDFLIPREGROW','# type of move'
    Write(unitno,'(a,t30,a)') 'Integer List', &
        '# species numbers into which to transform'

  End Subroutine idflipregrow_sampleCF

  !----------------------------------------------------------------------
  ! Displays the flip move parameters
  ! Requires:  params -- identity flip parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine idflipregrow_display(params,indent,unit)
    Type(Idflipregrow_Params), Intent(In)     :: params
    Integer, Intent(In)                 :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Identity Flip Regrow Move Parameters :"
    string1 = int2str(params%spc)
    Write(unit,'(4a,10i2)') blank,'Species ',Trim(string1),' -> ', &
        (params%tospc(i),i=1,params%ntospc)
    Write(unit,'(2a)') blank,'Regrow move information:'
    Call brmoves_displayparams(params%brmove,indent+2,unit)

  End Subroutine idflipregrow_display

  !----------------------------------------------------------------------
  ! Clean the flip move parameters
  ! Requires:  params -- identity flip parameters
  !----------------------------------------------------------------------
  Subroutine idflipregrow_clean(params)
    Type(Idflipregrow_Params), Intent(InOut)   :: params

    Integer            :: error

    Deallocate(params%tospc, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Call config_clean(params%molcopy)

  End Subroutine idflipregrow_clean

End Module idflipregrow
