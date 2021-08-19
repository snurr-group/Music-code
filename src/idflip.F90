!------------------------------------------------------------------------
! This module handles the identity flip move type.  This is just the 
! substitution of one molecule for another.
! [There is another movetype idflipregrow which works for branchedcoords,
!  That also has skeletal code for a route to rigidcoords moves, but not 
!  functional as of now, The rigidcoords moves from here uses some sort of 
!  postmove which I am not sure works -shaji , dec,16,2002]
!------------------------------------------------------------------------

Module idflip
  Use auxmoveparams, Only: AuxMoveObjects
  Use defaults, Only: RDbl, strLen, lstrLen, one
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop, int2str
  Use file, Only: file_open
  Use random, Only: rranf

  ! be vary careful, do not use cavity bias from this module, 
  ! Only postadjust should be done from here -shaji
  Use cavitylist, Only: Cavity_Params, &
      cavitylist_updatemolecinfo, cavitylist_checkandIncrMemory

  Use molecules, Only: molecules_getgcmodeltype
  Use config, Only: AtMolCoords, config_getnmoles, config_checkandincr, &
      config_setnmoles, config_clean, config_copymolec, config_delmol, &
      config_allocfields, config_idrebuild, config_getnatoms, &
      config_subset2xyz
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use subinteract, Only: Subset_Interactions, subinteract_changenmoles, &
      subinteract_int, subinteract_newnrg
  ! could be useful for idflip with biasrt
!!$  Use rigidmoves, Only: RigidBiasRT_Params, rigidmoves_biasrtInit,  &
!!$      rigidmoves_biasrt, rigidmoves_rtretry, rigidmoves_adjustBiasRT, &
!!$      rigidmoves_displayparams, rigidmoves_dyninfo

  Implicit None
  Save

  Private
  Public :: IDflip_Params,idflip_init,idflip_spcs,idflip_idstring, &
      idflip_move,idflip_restore,idflip_dyninfo,idflip_display,idflip_clean, &
      idflip_postadjust

  !** Simple Identity Flip parameters
  Type IDflip_Params
    Integer                          :: spc     !** species to operate on
    Integer                          :: ntospc
    Integer, Dimension(3)            :: osubset !** old subset (spc,mol,0)
    Integer, Dimension(3)            :: nsubset !** new subset (spc,mol,0)
    Real(kind=RDbl)                  :: max_nrg
    Integer, Dimension(:), Pointer   :: tospc
    Type(AtMolCoords)                :: molcopy !** storage for one molecule
    Type(Cavity_Params),Pointer      :: cavity
  End Type IDflip_Params

  Character(len=strLen), Parameter   :: idflip_idstring = 'IDFLIP'
  Real(kind=RDbl), Parameter         :: VERY_LRG_NRG = 1.0e10_RDbl

Contains
  !---------------------------------------------------------------------------
  ! Initializes the identity flip type move
  ! Requires:  params -- identity flip parameters
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine idflip_init(params,spc,auxmv,filename)
    Type(IDFlip_Params), Intent(Out)            :: params
    Integer, Intent(In)                         :: spc   
    Type(AuxMoveObjects),Pointer           :: auxmv
    Character(*), Intent(In)                    :: filename

    Integer                              :: error,unit,i
    Character(len=255)                   :: line     
    Character(len=strLen), Dimension(20) :: fields

    !** Defaults
    params%max_nrg = VERY_LRG_NRG
    params%osubset = (/0,0,0/)
    params%nsubset = (/0,0,0/)

    !** Initialize the cavity bais pointer
    Nullify(params%cavity)
    If (auxmv%cavBiasON)     params%cavity=>auxmv%cavity

    !** Set the species to operate on
    params%spc = spc

    !** Open the control file 
    unit = file_open(filename)

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

    !** Allocate space for storage of one molecule of this species
    Call config_allocfields(params%molcopy,params%spc,1,.True.)
!!$
!!$    !** Initialize the biased rotation/translation post move type
!!$    Call rigidmoves_biasrtInit(params%postmove,'NULL',(/20.0,2.0,0.0,2.0,0.0/))
!!$    Call rigidmoves_biasrtInit(params%postmove,filename)

  End Subroutine idflip_init

  !---------------------------------------------------------------------------
  ! Does the flip move based on the movetype and coordinate system
  ! type.  Returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply identity flip
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-----------------------------------------------------------------------------
  Logical Function idflip_move(params,subset,subints,species,simcell, &
      biasfactor)
    Type(IDFlip_Params), Intent(InOut)                     :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Integer     :: natoms2,nmoles1,nmoles2,spc2,mol2,maxmoves
    Logical     :: succ_flag,skip_intra_flag,fast

    !** Set useful values and defaults
    params%osubset = subset
    idflip_move = .True.
    skip_intra_flag = .False.
    fast = .True.  !** hard-coded HACK
    biasfactor = one

    !** Copy the specified molecule's configuration
    Call config_copymolec(species(subset(1)),subset(2),params%molcopy,1)

    !** Evaluate the energy of the removed molecule
    succ_flag = subinteract_int(subints(subset(1)),species,simcell,fast, &
        .False.,skip_intra_flag,(/params%max_nrg/),subset)
    Call checkandstop(succ_flag,__FILE__,__LINE__, &
        'Problem with nrg calculation of an existing structure')

    !** Remove the specified molecule and its space
    nmoles1 = config_getnmoles(species,subset(1))
    Call config_delmol(species,subset(1),subset(2))

    !** Get the new species and molecule number 
    spc2 = params%tospc(Int(rranf()*params%ntospc) + 1)

    nmoles2 = config_getnmoles(species,spc2) + 1
    natoms2 = config_getnatoms(species,spc2)
    mol2 = nmoles2
    params%nsubset = (/spc2,mol2,0/)

    !** Create space for a new molecules of spc2 
    Call config_setnmoles(species,spc2,nmoles2)
    Call config_checkandincr(species,spc2,nmoles2)
    succ_flag = subinteract_changenmoles(subints(spc2),spc2,nmoles2)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Create the new molecule using old molecule's generalized coordinates
    !    Call config_idrebuild(species(subset(1)),subset(2),species(spc2),mol2)
    Call config_idrebuild(params%molcopy,1,species(spc2),mol2)

    !** Apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(simcell, &
        species(spc2)%coords(1:natoms2,mol2)%rp, &
        species(spc2)%coords(1:natoms2,mol2)%r, &
        species(spc2)%coords(1:natoms2,mol2)%cr)

    !** Get the new energies
    idflip_move = subinteract_int(subints(spc2),species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc2,mol2,0/))
!!$ -Not microscopically reversible, bias for revers move is never estimated
!!$    !** If the energy isn't calculable, make random moves on the
!!$    !** molecule until it is (or maxmoves reached).  I believe this
!!$    !** is legal as long as they're random.  The new configuration 
!!$    !** could have been picked randomly, this is effectively doing that
!!$    If (.Not. idflip_move) Then
!!$      maxmoves = 100
!!$      Call rigidmoves_rtretry(params%postmove,maxmoves,species,spc2,mol2, &
!!$          simcell,idflip_move,subints(spc2))
!!$    End If
!!$
!!$    !** Make a biased rotation/translation move on the new molecule
!!$    Call rigidmoves_biasrt(params%postmove,species,spc2,mol2,simcell, &
!!$        biasfactor,idflip_move,subints(spc2))

  End Function idflip_move

  !--------------------------------------------------------------------------
  ! Restores the system to its state before the flip
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !-----------------------------------------------------------------------------
  Subroutine idflip_restore(params,subset,subints,species,simcell)
    Type(IDFlip_Params), Intent(In)                        :: params
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

  End Subroutine idflip_restore



  !--------------------------------------------------------------------------
  ! Removes old molecule from cavitylist, and add new molecule 
  ! Requires:  params  -- parameters for this move
  !            species -- species data structure
  !--------------------------------------------------------------------------
  Subroutine idflip_postadjust(params,species)
    Type(IDFlip_Params), Intent(InOut)                        :: params
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

  End Subroutine idflip_postadjust


  !------------------------------------------------------------------------
  ! Get the species number identities of the molecules changed.  The output
  ! follows the convention in subinteract_update, species numbers are 
  ! listed as positive if they are inserted, negative if they are deleted
  ! and as both if they are perturbed.
  ! Requires:  params  -- parameters for this move
  !            spcs  -- returned, changed species numbers
  !------------------------------------------------------------------------
  Subroutine idflip_spcs(params,spcs)
    Type(IDFlip_Params), Intent(In)     :: params
    Integer, Dimension(:), Intent(Out)  :: spcs

    spcs(1) = -params%osubset(1)
    spcs(2) = params%nsubset(1)

  End Subroutine idflip_spcs

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information 
  ! Requires:  unitno
  !----------------------------------------------------------------------------
  Subroutine idflip_sampleCF(unitno)
    Integer, Intent(In)           :: unitno

    Write(unitno,'(a,t30,a)') 'IDFLIP','# type of move'
    Write(unitno,'(a,t30,a)') 'Integer List', &
        '# species numbers into which to transform'

  End Subroutine idflip_sampleCF

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires: params -- random translation parameters
  !----------------------------------------------------------------------
  Function idflip_dyninfo(params)
    Character(len=strLen)               :: idflip_dyninfo
    Type(IDFlip_Params), Intent(In)     :: params

    idflip_dyninfo = " " ! rigidmoves_dyninfo(params%postmove)

  End Function idflip_dyninfo

  !----------------------------------------------------------------------
  ! Displays the flip move parameters
  ! Requires:  params -- identity flip parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine idflip_display(params,indent,unit)
    Type(IDFlip_Params), Intent(In)     :: params
    Integer, Intent(In)                 :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Identity Flip Move Parameters :"
    string1 = int2str(params%spc)
    Write(unit,'(4a,10i2)') blank,'Species ',Trim(string1),' -> ', &
        (params%tospc(i),i=1,params%ntospc)
!!$    Write(unit,'(2a)') blank,'Post move information:'
!!$    Call rigidmoves_displayparams(params%postmove,indent+2,unit)

  End Subroutine idflip_display

  !----------------------------------------------------------------------
  ! Clean the flip move parameters
  ! Requires:  params -- identity flip parameters
  !----------------------------------------------------------------------
  Subroutine idflip_clean(params)
    Type(IDFlip_Params), Intent(InOut)   :: params

    Integer            :: error

    Deallocate(params%tospc, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Call config_clean(params%molcopy)

  End Subroutine idflip_clean

End Module idflip
