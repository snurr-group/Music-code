!------------------------------------------------------------------------
! This module handles subspace jumps involving species identity changes.
! The expansion move takes a molecule of one species and splits it apart
! to form a pair of molecules of different species.  The contraction move
! takes a pair of molecules and reforms them into another species.  This
! functionality is intended to provide a way to hop from a reaction 
! dividing surface represented by one species to the full space 
! represented by a pair of reactants.
!------------------------------------------------------------------------

Module dimhop

  Use defaults, Only: RDbl, strLen, lstrLen, zero
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop, int2str, swap, findint, toreal
  Use file, Only: file_open
  Use random, Only: rranf
  Use distrib, Only: Distribution_Params, distrib_pick, distrib_prob, &
      distrib_init, distrib_display, distrib_dump, distrib_clean
  Use vector, Only: VecType, Assignment(=), Operator(+), Operator(-), &
      Operator(*), Operator(/), vector_getunitvec, vector_getnorm, mag, &
      vector_display, swap
  Use matrixops, Only: matrixops_ludcmp, matrixops_lubksb
  Use molecules, Only: molecules_getnatoms,molecules_getatomcoords, &
      molecules_getmaxnatoms,molecules_getnsorbs
  Use config, Only: AtMolCoords, config_getnmoles, config_checkandincr, &
      config_setnmoles, config_clean, config_copymolec, config_delmol, &
      config_allocfields, config_idrebuild, config_getnatoms, config_applypbcs, &
      config_getrp, config_subset2xyz, config_setfromrp, config_getfixed
  Use molmap, Only: MolecMolec_Map,molmap_init,molmap_map,molmap_mapatm, &
      molmap_findbrokenbonds,molmap_brokenbond,molmap_display,molmap_clean, &
      molmap_whichspc,molmap_doublemap,molmap_destspcs
  Use atomadjust, Only: AtomAdjust_Params, atomadjust_init, atomadjust_move, &
      atomadjust_clean, atomadjust_display
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use subinteract, Only: Subset_Interactions, subinteract_changenmoles, &
      subinteract_int
  Use rigidmoves, Only: RigidBiasRT_Params, rigidmoves_biasrtInit,  &
      rigidmoves_biasrt, rigidmoves_adjustBiasRT, rigidmoves_displayparams, &
      rigidmoves_dyninfo

  Implicit None
  Save

  Private
  Public :: Expansion_Params, Contraction_Params, dimhop_expand_idstring, &
      dimhop_contract_idstring, dimhop_initexpand, dimhop_initcontract, &
      dimhop_expand, dimhop_contract, dimhop_unexpand, dimhop_uncontract, &
      dimhop_spcs, dimhop_dyninfo, dimhop_dispexpand, &
      dimhop_dispcontract, dimhop_cleanexpand, &
      dimhop_cleancontract

  !** Expansion to two separated species 
  Type Expansion_Params
    Integer                          :: spc     !** species to operate on
    Integer                          :: dspc1,dspc2
    Integer, Dimension(3)            :: osubset !** old subset (spc,mol,0)
    Integer, Dimension(3)            :: nsubset1,nsubset2 !** new subsets
    Real(kind=RDbl)                  :: max_nrg,gaussctr,sigma,binsize,halfwidth
    Type(VecType)                    :: sepvec
    Type(Distribution_Params)        :: distrib
    Type(AtMolCoords)                :: molcopy !** storage for one molecule
    Type(MolecMolec_Map)             :: mmmap   !** molecule map for spc
    Type(AtomAdjust_Params), Dimension(2) :: adjust
  End Type Expansion_Params

  !** Expansion to two separated species 
  Type Contraction_Params
    Integer                          :: ospc1,ospc2  !** species to operate on
    Integer                          :: ntospc
    Integer, Dimension(:), Pointer   :: tospc
    Integer, Dimension(3)            :: nsubset !** new subset (spc,mol,0)
    Integer, Dimension(3)            :: osubset1,osubset2 !** old subsets 
    Real(kind=RDbl)                  :: max_nrg,gaussctr,sigma,binsize,halfwidth
    Type(VecType)                    :: sepvec
    Type(Distribution_Params)        :: distrib
    Type(AtMolCoords), Dimension(2)  :: molcopy !** storage for old molecules
    Real(kind=RDbl), Dimension(:), Pointer         :: sepdist
    Type(MolecMolec_Map), Dimension(:), Pointer    :: mmmap
    Type(AtomAdjust_Params), Dimension(:), Pointer :: adjust
  End Type Contraction_Params

  Character(len=strLen), Parameter    :: dimhop_expand_idstring = 'IDEXPAND'
  Character(len=strLen), Parameter    :: dimhop_contract_idstring = 'IDCONTRACT'
  Character(len=strLen), Dimension(2) :: dimhop_idstrings = &
      (/dimhop_expand_idstring,dimhop_contract_idstring/)
  Real(kind=RDbl), Parameter         :: VERY_LRG_NRG = 1.0e10_RDbl

  Interface dimhop_spcs
    Module Procedure dimhop_expand_spcs
    Module Procedure dimhop_contract_spcs
  End Interface
  
Contains
  !---------------------------------------------------------------------------
  ! Initializes the EXPANSION dimension hop move type
  ! Requires:  params -- identity flip parameters
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine dimhop_initexpand(params,spc,filename)
    Type(Expansion_Params), Intent(Out)         :: params
    Integer, Intent(In)                         :: spc   
    Character(*), Intent(In)                    :: filename
    
    Integer                              :: error,unit,n
    Character(len=255)                   :: line     
    Character(len=strLen), Dimension(20) :: fields

    !** Defaults
    params%max_nrg = VERY_LRG_NRG
    params%osubset = (/0,0,0/)
    params%nsubset1 = (/0,0,0/)
    params%nsubset2 = (/0,0,0/)
    params%sepvec = 0.0_RDbl

    !** Set the species to operate on
    params%spc = spc

    !** Open the control file 
    unit = file_open(filename)

    !** Read the line containing separation distance
    Read(unit,'(a)') line
    line = stripcmnt(line)
    n = split(line,fields,',')
    params%gaussctr = toreal(fields(1),'separation distance for gaussian center')
    params%sigma = toreal(fields(2),'separation distance for gaussian center')

    !** Initialize the Gaussian distribution
    params%binsize = 0.02_RDbl
    params%halfwidth = params%sigma*2.0_RDbl
    Call distrib_init(params%distrib,'GAUSSIAN', &
        (/params%sigma,params%gaussctr,params%halfwidth,params%binsize/))

    !** Read the line containing destination species numbers
    Read(unit,'(a)') line
    line = stripcmnt(line)
    n = split(line,fields,',')
    If (n /= 2) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' number of resulting species must be two for expansion move'
      Stop
    End If

    !** Place the destination species numbers in a list
    params%dspc1 = toint(fields(1),'destination species #1')
    params%dspc2 = toint(fields(2),'destination species #2')

    !** Allocate space for storage of one molecule of the original species
    Call config_allocfields(params%molcopy,params%spc,1,.True.)

    !** Get the molecule maps for the origin species
    Call molmap_init(params%mmmap,params%spc)
    Call molmap_findbrokenbonds(params%mmmap)

    !** Get the atomic adjustment parameters for the destination species
    Call atomadjust_init(params%adjust(1),params%dspc1)
    Call atomadjust_init(params%adjust(2),params%dspc2)
    
  End Subroutine dimhop_initexpand

  !---------------------------------------------------------------------------
  ! Initializes the dimension hop move type
  ! Requires:  params -- identity flip parameters
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine dimhop_initcontract(params,spc,filename)
    Type(Contraction_Params), Intent(Out)       :: params
    Integer, Intent(In)                         :: spc   
    Character(*), Intent(In)                    :: filename
    
    Integer                              :: error,unit,i,dspc,nbroken,n
    Type(VecType)                        :: avec1,avec2,refvec
    Character(len=255)                   :: line     
    Integer, Dimension(2)                :: bondpair,spcpair
    Character(len=strLen), Dimension(20) :: fields

    !** Defaults
    params%max_nrg = VERY_LRG_NRG
    params%nsubset = (/0,0,0/)
    params%osubset1 = (/0,0,0/)
    params%osubset2 = (/0,0,0/)
    params%sepvec = 0.0_RDbl

    !** Set the species to operate on
    params%ospc1 = spc

    !** Open the control file 
    unit = file_open(filename)

    !** Read the cooperating species
    Read(unit,*) params%ospc2

    !** Read the line containing separation distance
    Read(unit,'(a)') line
    line = stripcmnt(line)
    n = split(line,fields,',')
    params%gaussctr = toreal(fields(1),'separation distance for gaussian center')
    params%sigma = toreal(fields(2),'separation distance for gaussian center')

    !** Initialize the Gaussian distribution
    params%binsize = 0.02_RDbl
    params%halfwidth = params%sigma*2.0_RDbl
    Call distrib_init(params%distrib,'GAUSSIAN', &
        (/params%sigma,params%gaussctr,params%halfwidth,params%binsize/))

    !** Read the line containing destination species numbers
    Read(unit,'(a)') line
    line = stripcmnt(line)

    !** Place the destination species numbers in a list
    params%ntospc = split(line,fields,',')
    Allocate(params%tospc(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,params%ntospc
      params%tospc(i) = toint(fields(i))
    End Do

    !** Allocate space for storage of one molecule of each origin species
    Call config_allocfields(params%molcopy(1),params%ospc1,1,.True.)
    Call config_allocfields(params%molcopy(2),params%ospc2,1,.True.)

    !** Allocate space for a molecule map of each 'tospc'
    Allocate(params%mmmap(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    
    !** Get the molecule maps for these species
    Do i = 1,params%ntospc
      Call molmap_init(params%mmmap(i),params%tospc(i))
      Call molmap_findbrokenbonds(params%mmmap(i))
    End Do

    !** Get the atomic adjustment parameters for the destination species
    Allocate(params%adjust(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,params%ntospc
      Call atomadjust_init(params%adjust(i),params%tospc(i))
    End Do

    !** Calculate the reference separation distances in the molecule structure 
    Allocate(params%sepdist(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    Do i = 1,params%ntospc
      dspc = params%tospc(i)
      Call molmap_brokenbond(params%mmmap(i),1,spcpair,bondpair,nbroken)
      If (nbroken > 1) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' Number of broken bonds more than one, expansion method undefined'
        Stop
      End If
      avec1 = molecules_getatomcoords(dspc,bondpair(1))
      avec2 = molecules_getatomcoords(dspc,bondpair(2))
      refvec = avec2 - avec1
      params%sepdist(i) = vector_getnorm(refvec)
    End Do

  End Subroutine dimhop_initcontract

  !-----------------------------------------------------------------------------
  ! Performs the expansion move and returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply identity flip
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-----------------------------------------------------------------------------
  Logical Function dimhop_expand(params,subset,subints,species,simcell, &
      biasfactor)
    Type(Expansion_Params), Intent(InOut)                  :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Integer           :: onmoles
    Integer           :: spc1,mol1,nmoles1
    Integer           :: spc2,mol2,nmoles2
    Logical           :: succ_flag,skip_intra_flag,fast
    Real(kind=RDbl)   :: newdist,prob,cutoff

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call config_subset2xyz(species,'start.xyz','pre-expand',.False.,subset)

    !** Set useful values and defaults
    biasfactor = 1.0_RDbl
    params%osubset = subset
    dimhop_expand = .True.
    skip_intra_flag = .False.
    fast = .True.  !** hard-coded HACK

    !** Get new species #1 with molecule number and make space for it
    spc1 = params%dspc1
    nmoles1 = config_getnmoles(species,spc1) + 1
    mol1 = nmoles1
    params%nsubset1 = (/spc1,mol1,0/)
    Call config_setnmoles(species,spc1,nmoles1)
    Call config_checkandincr(species,spc1,nmoles1)
    succ_flag = subinteract_changenmoles(subints(spc1),spc1,nmoles1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Create new molecule of species #1 by mapping from original species
    Call molmap_map(params%mmmap,species,subset(1),subset(2), &
        spc1,species(spc1)%coords(:,mol1)%rp)

    !** Get new species #2 with molecule number and make space for it
    spc2 = params%dspc2
    nmoles2 = config_getnmoles(species,spc2) + 1
    mol2 = nmoles2
    params%nsubset2 = (/spc2,mol2,0/)
    Call config_setnmoles(species,spc2,nmoles2)
    Call config_checkandincr(species,spc2,nmoles2)
    succ_flag = subinteract_changenmoles(subints(spc2),spc2,nmoles2)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Create new molecule of species #2 by mapping from original species
    Call molmap_map(params%mmmap,species,subset(1),subset(2), &
        spc2,species(spc2)%coords(:,mol2)%rp)

    !** Make adjustments on the two molecules 
    Call atomadjust_move(params%adjust(1),species(spc1)%coords(:,mol1)%rp)
    Call atomadjust_move(params%adjust(2),species(spc2)%coords(:,mol2)%rp)

    !** Separate the two new molecules
    newdist = dimhop_separate(params,species,spc1,mol1,spc2,mol2)

    !** Get the probability of picking the new separation distance
    prob = distrib_prob(params%distrib,newdist)/params%binsize
    biasfactor = 1.0_RDbl  !** HaCK

    !** Check atom-atom distance overlap criteria, reject if necessary
    cutoff = 1.0_RDbl
    If (dimhop_chkovlp(species,cutoff,params%mmmap,(/spc1,mol1,0/), &
       (/spc2,mol2,0/))) Then
      dimhop_expand = .False.
    End If

    !** Necessary stuff if the move might still be accepted
    If (dimhop_expand) Then
      !** Apply periodic boundary conditions to new molecules and reset gcoords
      Call config_applypbcs(species,simcell,(/spc1,mol1,0/))
      Call config_applypbcs(species,simcell,(/spc2,mol2,0/))
      Call config_setfromrp(species,spc1,mol1)
      Call config_setfromrp(species,spc2,mol2)
  
      !** Evaluate the energy of the removed molecule
      succ_flag = subinteract_int(subints(subset(1)),species,simcell,fast, &
          .False.,skip_intra_flag,(/params%max_nrg/),subset)
      Call checkandstop(succ_flag,__FILE__,__LINE__, &
          'Problem with nrg calculation of an existing structure')
    End If

    Call config_subset2xyz(species,'end.xyz','post-expand',.False., &
        (/spc1,mol1,0/),(/spc2,mol2,0/))
  
    !** Copy the specified molecule's configuration enabling restoration
    Call config_copymolec(species(subset(1)),subset(2),params%molcopy,1)

    !** Remove the specified molecule and its space
    onmoles = config_getnmoles(species,subset(1))
    Call config_delmol(species,subset(1),subset(2))

    !** Exit if we already know this won't be accepted
    If (.Not. dimhop_expand) Return

    !** Get the new energies
    dimhop_expand = subinteract_int(subints(spc1),species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc1,mol1,0/))
    If (.Not. dimhop_expand) Return
    dimhop_expand = subinteract_int(subints(spc2),species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc2,mol2,0/))

  End Function dimhop_expand

  !-----------------------------------------------------------------------------
  ! Performs the expansion move and returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply identity flip
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-----------------------------------------------------------------------------
  Logical Function dimhop_contract(params,subset,subints,species,simcell, &
      biasfactor)
    Type(Contraction_Params), Intent(InOut)                :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Integer           :: nmoles1,nmoles2,mol2,spc2
    Integer           :: nmoles,spc,mol,n,mapno,betterspc
    Logical           :: succ_flag,skip_intra_flag,fast
    Real(kind=RDbl)   :: olddist,prob,cutoff

    !** Set useful values and defaults
    biasfactor = 1.0_RDbl
    params%osubset1 = subset
    dimhop_contract = .True.
    skip_intra_flag = .False.
    fast = .True.  !** hard-coded HACK

    !** Pick a molecule of the second origin species randomly
    spc2 = params%ospc2
    nmoles2 = config_getnmoles(species,spc2)
    mol2 = Int(rranf()*nmoles2) + 1
    params%osubset2 = (/spc2,mol2,0/)

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call config_subset2xyz(species,'start2.xyz','pre-contract',.False., &
        subset,(/spc2,mol2,0/))

    !** Make space for a new molecule of the destination species
    n = Int(rranf()*params%ntospc) + 1
    spc = params%tospc(n)
    nmoles = config_getnmoles(species,spc) + 1
    mol = nmoles
    params%nsubset = (/spc,mol,0/)
    Call config_setnmoles(species,spc,nmoles)
    Call config_checkandincr(species,spc,nmoles)
    succ_flag = subinteract_changenmoles(subints(spc),spc,nmoles)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Get the map number for the new species
    mapno = findint(params%tospc,spc)

    !** Create new molecule by mapping from original species
    Call molmap_map(params%mmmap(mapno),species,subset(1),subset(2), &
        spc,species(spc)%coords(:,mol)%rp)
    Call molmap_map(params%mmmap(mapno),species,spc2,mol2, &
        spc,species(spc)%coords(:,mol)%rp)

    !** Make adjustments on the new molecule
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call atomadjust_move(params%adjust(mapno),species(spc)%coords(:,mol)%rp)

    !** Contract the pieces of the new molecule together and apply PBCs
    olddist = dimhop_together(params,species,spc,mol)

    !** Get the probability of picking the old separation distance
    prob = distrib_prob(params%distrib,olddist)/params%binsize
    biasfactor = 1.0_RDbl  !** HaCK

    !** Make sure we're not outside distribution range, if so, reject move
    If (prob == zero) dimhop_contract = .False.

    !** Check atom-atom distance overlap criteria, reject if necessary
    cutoff = Min(1.5_RDbl,params%sepdist(mapno))
    If (dimhop_chkovlp(species,cutoff,params%mmmap(mapno),(/spc,mol,0/))) Then
      dimhop_contract = .False.
    End If

    !** Check destination species similarity rejection criteria (bond distances)
    If (dimhop_contract) Then
      betterspc = dimhop_chkbonds(params,species,spc,mol)
      If (betterspc /= 0) Then
        Write(0,'(1x,2a,i4,a,i4,a)') __FILE__," : ",__LINE__, &
            ' species possibility ',betterspc,' has bond length below reference'
        dimhop_contract = .False.
      End If
    End If

    !** Necessary stuff if the move might still be accepted
    If (dimhop_contract) Then
      !** Apply PBCs and reset the generalized coordinates
      Call config_applypbcs(species,simcell,(/spc,mol,0/))
      Call config_setfromrp(species,spc,mol)
  
      !** Evaluate the energy of the molecules to be removed
      succ_flag = subinteract_int(subints(subset(1)),species,simcell,fast, &
          .True.,skip_intra_flag,(/params%max_nrg/),subset)
      Call checkandstop(succ_flag,__FILE__,__LINE__, &
          'Problem with nrg calculation of an existing structure')
      succ_flag = subinteract_int(subints(spc2),species,simcell,fast, &
          .True.,skip_intra_flag,(/params%max_nrg/),(/spc2,mol2,0/))
      Call checkandstop(succ_flag,__FILE__,__LINE__, &
          'Problem with nrg calculation of an existing structure')
    End If

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Call config_subset2xyz(species,'end2.xyz','post-contract', &
        .False.,(/spc,mol,0/))

    !** Copy the origin molecules enabling restoration
    Call config_copymolec(species(params%osubset1(1)), &
        params%osubset1(2),params%molcopy(1),1)
    Call config_copymolec(species(params%osubset2(1)), &
        params%osubset2(2),params%molcopy(2),1)

    !** Remove the origin molecules from the configuration
    Call config_delmol(species,subset(1),subset(2))
    Call config_delmol(species,spc2,mol2)

    !** Exit if we already know this won't be accepted
    If (.Not. dimhop_contract) Return

    !** Get the energy of the new molecule
    dimhop_contract = subinteract_int(subints(subset(1)),species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),subset)

  End Function dimhop_contract

  !-----------------------------------------------------------------------------
  ! Separates the two new molecules during an expansion move.  Operates only
  ! on the parent coordinates.
  ! Requires:  params  -- parameters for this move
  !            species -- species data structure
  !            spc1 -- first species 
  !            mol1 -- first molecule
  !            spc2 -- second species 
  !            mol2 -- second molecule
  !-----------------------------------------------------------------------------
  Real(kind=RDbl) Function dimhop_separate(params,species,spc1,mol1,spc2,mol2)
    Type(Expansion_Params), Intent(InOut)            :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Integer, Intent(In)                              :: spc1,mol1,spc2,mol2

    Integer                 :: nbroken,atm,atm1,atm2
    Real(kind=RDbl)         :: dist
    Integer, Dimension(2)   :: bondpair,spcpair
    Type(VecType)           :: bondvec,unitvec,dispvec

    !** Get the unit vector along the broken bond
    !** The bond vector points in the direction of the second new species
    Call molmap_brokenbond(params%mmmap,1,spcpair,bondpair,nbroken)
    If (nbroken > 1) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Number of broken bonds more than one, expansion method undefined'
      Stop
    End If
    If (spcpair(1) /= spc1) Then
      Call swap(spcpair(1),spcpair(2))
      Call swap(bondpair(1),bondpair(2))
    End If
    atm1 = molmap_mapatm(params%mmmap,bondpair(1))
    atm2 = molmap_mapatm(params%mmmap,bondpair(2))
    bondvec = config_getrp(species,(/spc2,mol2,atm2/)) - &
        config_getrp(species,(/spc1,mol1,atm1/))
    unitvec = vector_getunitvec(bondvec)

    !** Draw the separation distance from a distribution
    dist = distrib_pick(params%distrib)
    dimhop_separate = dist

    !** Move the first molecule 
    dispvec = (-0.5_RDbl*dist)*unitvec
    Do atm = 1,molecules_getnatoms(spc1)
      species(spc1)%coords(atm,mol1)%rp = &
          species(spc1)%coords(atm,mol1)%rp + dispvec
    End Do

    !** Move the second molecule 
    dispvec = 0.5_RDbl*dist*unitvec
    Do atm = 1,molecules_getnatoms(spc2)
      species(spc2)%coords(atm,mol2)%rp = &
          species(spc2)%coords(atm,mol2)%rp + dispvec
    End Do
    
    !** Save the unit displacement vector
    params%sepvec = unitvec

  End Function dimhop_separate

  !-----------------------------------------------------------------------------
  ! Brings two pieces of a single molecule together after a contraction move.
  ! Returns the original separation distance.
  ! Requires:  params  -- parameters for this move
  !            species -- species data structure
  !            spc -- first species 
  !            mol -- first molecule
  !-----------------------------------------------------------------------------
  Real(kind=RDbl) Function dimhop_together(params,species,spc,mol)
    Type(Contraction_Params), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Integer, Intent(In)                              :: spc,mol

    Integer                 :: nbroken,atm,mapno,ospc,dspc
    Real(kind=RDbl)         :: dist
    Integer, Dimension(2)   :: bondpair,spcpair
    Type(VecType)           :: dispvec1,dispvec2,avec1,avec2
    Type(VecType)           :: bondvec,unitvec

    !** Get the map number of the species in the contraction parameters
    mapno = findint(params%tospc,spc)

    !** Get the unit vector along the broken bond
    !** The bond vector points in the direction of the second new species
    Call molmap_brokenbond(params%mmmap(mapno),1,spcpair,bondpair,nbroken)
    If (nbroken > 1) Then
      Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
          ' Number of broken bonds more than one, expansion method undefined'
      Stop
    End If
    bondvec = config_getrp(species,(/spc,mol,bondpair(2)/)) - &
        config_getrp(species,(/spc,mol,bondpair(1)/))
    unitvec = vector_getunitvec(bondvec)

    !** Calculate the distance the pieces need to be moved together
    dimhop_together = vector_getnorm(bondvec)
    dist = dimhop_together - params%sepdist(mapno)

    !** Calculate the displacement vector for the two pieces
    dispvec1 = (0.5_RDbl*dist)*unitvec
    dispvec2 = (-0.5_RDbl*dist)*unitvec
    If (spcpair(1) /= params%ospc1) Call swap(dispvec1,dispvec2)

    !** Move the two pieces
    Do atm = 1,molecules_getnatoms(spc)
      ospc = molmap_whichspc(params%mmmap(mapno),atm)
      If (ospc == params%ospc1) Then
        species(spc)%coords(atm,mol)%rp = &
            species(spc)%coords(atm,mol)%rp + dispvec1
      Else If (ospc == params%ospc2) Then
        species(spc)%coords(atm,mol)%rp = &
            species(spc)%coords(atm,mol)%rp + dispvec2
      Else
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Stop
      End If
    End Do

    !** Calculate the corresponding distance in the molecule structure 
    avec1 = molecules_getatomcoords(spc,bondpair(1))
    avec2 = molecules_getatomcoords(spc,bondpair(2))
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Write(*,*) 'CONTRACTED TO: ',vector_getnorm(avec2 - avec1)

    If ((vector_getnorm(avec2 - avec1) - &
        vector_getnorm(bondvec)) > 0.001_RDbl) Then
      Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
      Stop
    End If
    
    !** Save the unit displacement vector
    params%sepvec = unitvec

  End Function dimhop_together

  !-----------------------------------------------------------------------------
  ! Check if the chosen destination species in a contraction move has broken
  ! bond length to any of the other possible species that are within the 
  ! critical bond lengths of those species.  This would mean that the new 
  ! configuration is more like another species than the intended destination
  ! species and would be a good criteria to reject the move before the 
  ! forcefield calculation.  Returns zero if all is ok, otherwise it returns
  ! the first destination species number in which the critical bond length
  ! is less than in the reference structure.
  ! Requires:  params  -- parameters for the contraction move
  !            species -- species data structure
  !            spc -- species number of freshly created molecule
  !            mol -- molecule number
  !-----------------------------------------------------------------------------
  Integer Function dimhop_chkbonds(params,species,spc,mol)
    Type(Contraction_Params), Intent(InOut)          :: params
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Integer, Intent(In)                              :: spc,mol

    Integer                 :: nbroken,max,nfixed,error,natoms
    Integer                 :: atm,ospc,dspc,mapno,refmapno
    Logical                 :: firsttime = .True.
    Real(kind=RDbl)         :: dist
    Type(VecType)           :: bondvec
    Integer, Dimension(2)   :: bondpair,spcpair,spclist
    Type(VecType), Dimension(:), Allocatable, Save  :: tempcoords

    !** Default 
    dimhop_chkbonds = 0
    refmapno = findint(params%tospc,spc)

    !** Allocate temporary space for coordinates
    If (firsttime) Then
      firsttime = .False.
      nfixed = config_getfixed(species,spclist)
      max = molecules_getmaxnatoms(spclist(1:nfixed))
      Allocate(tempcoords(max),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If 

    !** Check the bond distances to the other destination species 
    !** Return a .False. flag if bond length is less than that for another
    If (params%ntospc > 1) Then
      Do mapno = 1,params%ntospc
        dspc = params%tospc(mapno)
        If (dspc == spc) Cycle

        natoms = molecules_getnatoms(dspc)
        Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
        Write(*,*) params%ospc1,dspc,natoms
        Write(*,*) params%tospc

        !** Create a temporary new molecule by mapping from original species
        Call molmap_doublemap(params%mmmap(refmapno),params%mmmap(mapno), &
            species,spc,mol,dspc,tempcoords(1:natoms))

        !** Get the atoms that form the broken bond for this destination species
        Call molmap_brokenbond(params%mmmap(mapno),1,spcpair,bondpair,nbroken)
        If (nbroken > 1) Then
          Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
              ' Number of broken bonds more than one, expansion method undefined'
          Stop
        End If

        !** Get the broken bond vector for this destination species
        bondvec = tempcoords(bondpair(2)) - tempcoords(bondpair(1))
        dist = vector_getnorm(bondvec)
        If (dist < params%sepdist(mapno)) dimhop_chkbonds = dspc

      End Do
    End If

  End Function dimhop_chkbonds

  !--------------------------------------------------------------------------
  ! Check if the new species created during the expansion or contraction
  ! moves has atoms that would overlap with itself or it's co-fragment.  
  ! After an expansion move, pass both subset1,subset2, but after a 
  ! contraction move, pass only the new molecule info in subset1.  Basically
  ! the routine identifies two fragments and checks atom-atom distances 
  ! between the fragments.
  ! Requires:  species -- species data structure
  !            cutoff -- atom-atom distance under which they are overlapping
  !            mmmap  -- molecule map to split subset1 if necessary
  !            subset1 -- species and molecule number of 1st frag or new mol
  !            subset2 -- optional info for 2nd fragment (expansion only)
  !--------------------------------------------------------------------------
  Logical Function dimhop_chkovlp(species,cutoff,mmmap,subset1,subset2)
    Type(AtMolCoords), Dimension(:), Intent(InOut)   :: species
    Real(kind=RDbl), Intent(In)                      :: cutoff
    Type(MolecMolec_Map), Intent(In)                 :: mmmap
    Integer, Dimension(3), Intent(In)                :: subset1
    Integer, Dimension(3), Intent(In), Optional      :: subset2

    Integer                 :: max,error,nfixed,n
    Integer                 :: natoms1,natoms2,atm1,atm2
    Logical                 :: firsttime = .True.
    Real(kind=RDbl)         :: dist
    Type(VecType)           :: bondvec
    Integer, Dimension(2)   :: spclist
    Type(VecType), Dimension(:,:), Allocatable, Save  :: tempcoords

    !** Default 
    dimhop_chkovlp = .False.

    !** Allocate temporary space for coordinates
    If (firsttime) Then
      firsttime = .False.
      nfixed = config_getfixed(species,spclist)
      max = molecules_getmaxnatoms(spclist(1:nfixed))
      Allocate(tempcoords(2,max),STAT=error)
      If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    End If 

    !** Copy or create the two fragments
    If (Present(subset2)) Then  !** just copy into temporary arrays
      natoms1 = molecules_getnatoms(subset1(1))
      natoms2 = molecules_getnatoms(subset2(1))
      tempcoords(1,1:natoms1) = &
          species(subset1(1))%coords(1:natoms1,subset1(2))%rp
      tempcoords(1,1:natoms2) = &
          species(subset2(1))%coords(1:natoms2,subset2(2))%rp

    Else  !** Use molmap to split subset1 into fragments
      n = molmap_destspcs(mmmap,spclist)
      If (n > 2) Then
        Write(0,'(1x,2a,i4,2a)') __FILE__," : ",__LINE__, &
            ' Molecule map with more than two maps passed'
        Stop
      End If
      natoms1 = molecules_getnatoms(spclist(1))
      natoms2 = molecules_getnatoms(spclist(2))
      Call molmap_map(mmmap,species,subset1(1),subset1(2), &
          spclist(1),tempcoords(1,1:natoms1))
      Call molmap_map(mmmap,species,subset1(1),subset1(2), &
          spclist(2),tempcoords(2,1:natoms2))
    End If

    !** Loop through the atom pairs in separate fragments and check distances
    Do atm1 = 1,natoms1
      Do atm2 = 1,natoms2
        dist = vector_getnorm(tempcoords(1,atm1) - tempcoords(2,atm2))
        If (dist < cutoff) Then
          dimhop_chkovlp = .True.
          Return
        End If
      End Do
    End Do

  End Function dimhop_chkovlp

  !-----------------------------------------------------------------------------
  ! Restores the system to its state before the expansion move
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !-----------------------------------------------------------------------------
  Subroutine dimhop_unexpand(params,subset,subints,species,simcell)
    Type(Expansion_Params), Intent(In)                     :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Integer     :: onmoles,nmoles1,nmoles2
    Logical     :: succ_flag

    !** Get the number of molecules of each important type
    onmoles = config_getnmoles(species,params%osubset(1))
    nmoles1 = config_getnmoles(species,params%nsubset1(1))
    nmoles2 = config_getnmoles(species,params%nsubset2(1))

    !** Remove the two new molecules from the new species list
    Call config_setnmoles(species,params%nsubset1(1),nmoles1-1)
    Call config_setnmoles(species,params%nsubset2(1),nmoles2-1)
    succ_flag = subinteract_changenmoles(subints(params%nsubset1(1)), &
        params%nsubset1(1),nmoles1-1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')
    succ_flag = subinteract_changenmoles(subints(params%nsubset1(1)), &
        params%nsubset2(1),nmoles2-1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Copy the molecule in the previous old molecule position 
    !** back to the end of the old species list
    Call config_copymolec(species(params%osubset(1)),params%osubset(2), &
        species(params%osubset(1)),onmoles+1)
    Call config_setnmoles(species,params%osubset(1),onmoles+1)

    !** Put the removed old molecule back into its old position
    Call config_copymolec(params%molcopy,1,species(params%osubset(1)), &
        params%osubset(2))

  End Subroutine dimhop_unexpand

  !-----------------------------------------------------------------------------
  ! Restores the system to its state before the expansion move
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !-----------------------------------------------------------------------------
  Subroutine dimhop_uncontract(params,subset,subints,species,simcell)
    Type(Contraction_Params), Intent(In)                   :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell

    Integer     :: onmoles1,onmoles2,nmoles
    Logical     :: succ_flag

    !** Get the number of molecules of each important type
    onmoles1 = config_getnmoles(species,params%osubset1(1))
    onmoles2 = config_getnmoles(species,params%osubset2(1))
    nmoles = config_getnmoles(species,params%nsubset(1))

    !** Remove the new molecule from the new species list
    Call config_setnmoles(species,params%nsubset(1),nmoles-1)
    succ_flag = subinteract_changenmoles(subints(params%nsubset(1)), &
        params%nsubset(1),nmoles-1)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Copy the current molecules in the previous old molecule positions 
    !** back to the end of the old species list
    Call config_copymolec(species(params%osubset1(1)),params%osubset1(2), &
        species(params%osubset1(1)),onmoles1+1)
    Call config_setnmoles(species,params%osubset1(1),onmoles1+1)
    Call config_copymolec(species(params%osubset2(1)),params%osubset2(2), &
        species(params%osubset2(1)),onmoles2+1)
    Call config_setnmoles(species,params%osubset2(1),onmoles2+1)

    !** Put the removed old molecules back into their old positions
    Call config_copymolec(params%molcopy(1),1,species(params%osubset1(1)), &
        params%osubset1(2))
    Call config_copymolec(params%molcopy(2),1,species(params%osubset2(1)), &
        params%osubset2(2))

  End Subroutine dimhop_uncontract

  !-----------------------------------------------------------------------------
  ! Get the species number identity of the new molecule created.  The output
  ! follows the convention in subinteract_update, species numbers are 
  ! listed as positive if they are inserted, negative if they are deleted
  ! and as both if they are perturbed.
  ! Requires:  params  -- parameters for this move
  !            spcs  -- returned, changed species numbers
  !-----------------------------------------------------------------------------
  Subroutine dimhop_expand_spcs(params,spcs)
    Type(Expansion_Params), Intent(In)     :: params
    Integer, Dimension(:), Intent(Out)     :: spcs

    spcs(1) = -params%osubset(1)
    spcs(2) = params%nsubset1(1)
    spcs(3) = params%nsubset2(1)

  End Subroutine dimhop_expand_spcs

  !-----------------------------------------------------------------------------
  ! Get the species number identity of the new molecule created.  The output
  ! follows the convention in subinteract_update, species numbers are 
  ! listed as positive if they are inserted, negative if they are deleted
  ! and as both if they are perturbed.
  ! Requires:  params  -- parameters for this move
  !            spcs  -- returned, changed species numbers
  !-----------------------------------------------------------------------------
  Subroutine dimhop_contract_spcs(params,spcs)
    Type(Contraction_Params), Intent(In)     :: params
    Integer, Dimension(:), Intent(Out)       :: spcs

    spcs(1) = -params%osubset1(1)
    spcs(2) = -params%osubset2(1)
    spcs(3) = params%nsubset(1)

  End Subroutine dimhop_contract_spcs

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information 
  ! Requires:  unitno -- unit number to dump into
  !----------------------------------------------------------------------------
  Subroutine dimhop_sampleCF(unitno)
    Integer, Intent(In)           :: unitno
    
    Write(unitno,'(a,t30,a)') 'IDEXPAND','# type of move'
    Write(unitno,'(a,t30,a)') 'Integer List', &
        '# species numbers into which to transform'
    Write(unitno,*) 
    Write(unitno,'(a,t30,a)') 'IDCONTRACT','# type of move'
    Write(unitno,'(a,t30,a)') 'Integer   ', &
        '# species number to contract with'
    Write(unitno,'(a,t30,a)') 'Integer List', &
        '# species numbers into which to transform'

  End Subroutine dimhop_sampleCF

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires: params -- random translation parameters
  !----------------------------------------------------------------------
  Function dimhop_dyninfo(params)
    Character(len=strLen)                  :: dimhop_dyninfo
    Type(Expansion_Params), Intent(In)     :: params

    dimhop_dyninfo = 'nothing to report'

  End Function dimhop_dyninfo

  !----------------------------------------------------------------------
  ! Displays the flip move parameters
  ! Requires:  params -- identity flip parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine dimhop_dispexpand(params,indent,unit)
    Type(Expansion_Params), Intent(In)     :: params
    Integer, Intent(In)                    :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2,string3

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Expansion Move Parameters :"
    string1 = int2str(params%spc)
    string2 = int2str(params%dspc1)
    string3 = int2str(params%dspc2)
    Write(unit,'(7a)') blank,'Species ',Trim(string1),' -> ', &
        Trim(string2),',',Trim(string3)
    Call distrib_display(params%distrib,indent+1,6)
    Write(unit,'(2a,i2)') blank,'Atom adjustments for species: ',params%dspc1
    Call atomadjust_display(params%adjust(1),indent+2,unit)
    Write(unit,'(2a,i2)') blank,'Atom adjustments for species: ',params%dspc2
    Call atomadjust_display(params%adjust(2),indent+2,unit)

  End Subroutine dimhop_dispexpand

  !----------------------------------------------------------------------
  ! Displays the flip move parameters
  ! Requires:  params -- identity flip parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine dimhop_dispcontract(params,indent,unit)
    Type(Contraction_Params), Intent(In)   :: params
    Integer, Intent(In)                    :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Contraction Move Parameters :"
    string1 = int2str(params%ospc1)
    string2 = int2str(params%ospc2)
    Write(unit,'(6a,10i2)') blank,'Species ',Trim(string1),',', &
        Trim(string2),' -> ',(params%tospc(i),i=1,params%ntospc)
    Call distrib_display(params%distrib,indent+1,6)

    Do i = 1,params%ntospc
      Write(unit,'(2a,i2)') blank,'Atom adjustments for species: ', &
          params%tospc(i)
      Call atomadjust_display(params%adjust(i),indent+2,unit)
    End Do

  End Subroutine dimhop_dispcontract

  !----------------------------------------------------------------------
  ! Clean the Expansion move parameters
  ! Requires:  params -- identity flip parameters
  !----------------------------------------------------------------------
  Subroutine dimhop_cleanexpand(params)
    Type(Expansion_Params), Intent(InOut)   :: params

    Call config_clean(params%molcopy)

    !** Deallocate the molecule map
    Call molmap_clean(params%mmmap)

    Call atomadjust_clean(params%adjust(1))
    Call atomadjust_clean(params%adjust(2))

  End Subroutine dimhop_cleanexpand

  !----------------------------------------------------------------------
  ! Clean the Contraction move parameters
  ! Requires:  params -- identity flip parameters
  !----------------------------------------------------------------------
  Subroutine dimhop_cleancontract(params)
    Type(Contraction_Params), Intent(InOut)   :: params

    Integer            :: error,i

    Deallocate(params%tospc, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Call config_clean(params%molcopy(1))
    Call config_clean(params%molcopy(2))

    !** Deallocate the molecule maps and adjustment parameters
    Do i = 1,params%ntospc
      Call molmap_clean(params%mmmap(i))
      Call atomadjust_clean(params%adjust(i))
    End Do
    Deallocate(params%mmmap, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

  End Subroutine dimhop_cleancontract
    
End Module dimhop
