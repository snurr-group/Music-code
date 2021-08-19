!------------------------------------------------------------------------
! This module handles the identity library move type.  A library of 
! configurations is built out of libraries for each possibility species.
! This library is randomly distributed.  To make a move, a configuration
! (of a random species) is picked from the library according to the 
! cannonical ensemble weighting.
!------------------------------------------------------------------------

Module idlibrary

  Use defaults, Only: RDbl, strLen, lstrLen
  Use utils, Only: filesrchstr, stripcmnt, split, toint, toupper, &
      allocErrDisplay, checkandstop, int2str
  Use file, Only: file_open
  Use random, Only: rranf
  Use molecules, Only: molecules_getgcmodeltype,molecules_name
  Use config, Only: AtMolCoords, config_getnmoles, config_checkandincr, &
      config_setnmoles, config_clean, config_copymolec, config_delmol, &
      config_allocfields, config_idrebuild, config_getnatoms, &
      config_subset2xyz
  Use conflib, Only: ConfLibrary, conflib_basicInit, conflib_dataInit, &
      conflib_getrandomconfig
  Use simcell, Only: SimCell_Params, simcell_pbc
  Use subinteract, Only: Subset_Interactions, subinteract_changenmoles, &
      subinteract_int, subinteract_newnrg

  Implicit None
  Save

  Private
  Public :: IDLibrary_Params,idlibrary_init,idlibrary_spcs, &
      idlibrary_idstring,idlibrary_move,idlibrary_restore, &
      idlibrary_dyninfo,idlibrary_display,idlibrary_clean

  !** Identity Library move parameters
  Type IDLibrary_Params
    Integer                        :: spc       !** species to operate on
    Integer                        :: ntospc  
    Integer, Dimension(3)          :: osubset   !** old subset (spc,mol,0)
    Integer, Dimension(3)          :: nsubset   !** new subset (spc,mol,0)
    Real(kind=RDbl)                :: max_nrg  
    Type(AtMolCoords)              :: molcopy   !** storage for one molecule
    Type(ConfLibrary)              :: library 
    Integer, Dimension(:), Pointer                 :: tospc
    character(len=lstrLen), Dimension(:), Pointer  :: libfile
  End Type IDLibrary_Params

  Character(len=strLen), Parameter   :: idlibrary_idstring = 'IDLIBRARY'
  Real(kind=RDbl), Parameter         :: VERY_LRG_NRG = 1.0e10_RDbl
  
Contains
  !---------------------------------------------------------------------------
  ! Initializes the identity library type move
  ! Requires:  params -- identity library parameters
  !            spc -- species number
  !            filename -- filename where initialization data can be found
  !---------------------------------------------------------------------------
  Subroutine idlibrary_init(params,spc,filename)
    Type(IDLibrary_Params), Intent(Out)         :: params
    Integer, Intent(In)                         :: spc   
    Character(*), Intent(In)                    :: filename
    
    Integer                              :: error,unit,i,n
    Character(len=255)                   :: line     
    Character(len=strLen), Dimension(20) :: fields
    Character(len=strLen), Dimension(30) :: names

    !** Defaults
    params%max_nrg = VERY_LRG_NRG
    params%osubset = (/0,0,0/)
    params%nsubset = (/0,0,0/)

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

    !** Read the line for and process the configuration filenames
    Read(unit,'(a)') line
    line = stripcmnt(line)
    Allocate(params%libfile(params%ntospc),STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)
    n = split(line,fields,',')
    If (n /= params%ntospc) Then
      Write(0,'(1x,2a,i4, 2a)') __FILE__," : ",__LINE__, &
          " Number of destination species does not match number of config files"
      Stop
    End If
    Do i = 1,params%ntospc
      params%libfile(i) = fields(i)
    End Do

    !** Allocate space for storage of one molecule of this species
    Call config_allocfields(params%molcopy,params%spc,1,.True.)

    !** Initialize the library basics
    Do i = 1,params%ntospc
      names(i) = molecules_name(params%tospc(i))
    End Do
    Call conflib_basicInit(params%library,names(1:params%ntospc), &
        .False.,.False.)

    !** Initialize the library from the configuration files
    Do i = 1,params%ntospc
      Call conflib_dataInit(params%library,names(i),params%libfile(i))
      params%libfile(i) = fields(i)
    End Do

  End Subroutine idlibrary_init

  !-----------------------------------------------------------------------------
  ! Does the library move based on the movetype and coordinate system
  ! type.  Returns True if the move was made.  
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply identity library
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !            biasfactor -- bias ratio of final structure to initial structure
  !-----------------------------------------------------------------------------
  Logical Function idlibrary_move(params,subset,subints,species,simcell, &
      biasfactor)
    Type(IDLibrary_Params), Intent(InOut)                  :: params
    Integer, Dimension(:), Intent(In)                      :: subset
    Type(Subset_Interactions), Dimension(:), Intent(InOut) :: subints
    Type(AtMolCoords), Dimension(:), Intent(InOut)         :: species
    Type(SimCell_Params), Intent(In)                       :: simcell
    Real(kind=RDbl), Intent(Out)                           :: biasfactor

    Integer     :: natoms2,nmoles1,nmoles2,spc2,mol2,maxmoves
    Logical     :: succ_flag,skip_intra_flag,fast

    !** Set useful values and defaults
    params%osubset = subset
    idlibrary_move = .True.
    skip_intra_flag = .False.
    fast = .True.  !** hard-coded HACK

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

    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT
    Stop

    !** Get the new species and molecule number 
!    spc2 = ?
    nmoles2 = config_getnmoles(species,spc2) + 1
    natoms2 = config_getnatoms(species,spc2)
    mol2 = nmoles2
    params%nsubset = (/spc2,mol2,0/)

    !** Create space for a new molecules of spc2 
    Call config_setnmoles(species,spc2,nmoles2)
    Call config_checkandincr(species,spc2,nmoles2)
    succ_flag = subinteract_changenmoles(subints(spc2),spc2,nmoles2)
    Call checkandstop(succ_flag,__FILE__,__LINE__,'changenmoles failed')

    !** Copy the new molecule configuration into the species storage
    Write(*,'(2a,i4)') __FILE__," : ",__LINE__  !**DEBUG_INSERT

    !** Apply periodic boundary conditions before forcefield evaluation
    Call simcell_pbc(simcell, &
        species(spc2)%coords(1:natoms2,mol2)%rp, &
        species(spc2)%coords(1:natoms2,mol2)%r, &
        species(spc2)%coords(1:natoms2,mol2)%cr)

    !** Get the new energies
    idlibrary_move = subinteract_int(subints(spc2),species,simcell,fast, &
        .True.,skip_intra_flag,(/params%max_nrg/),(/spc2,mol2,0/))

  End Function idlibrary_move

  !-----------------------------------------------------------------------------
  ! Restores the system to its state before the library move
  ! Requires:  params  -- parameters for this move
  !            subset -- system subset to which to apply move
  !            subints -- subset interactions for each species
  !            species -- species data structure
  !            simcell -- simulation cell information
  !-----------------------------------------------------------------------------
  Subroutine idlibrary_restore(params,subset,subints,species,simcell)
    Type(IDLibrary_Params), Intent(In)                     :: params
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

  End Subroutine idlibrary_restore

  !-----------------------------------------------------------------------------
  ! Get the species number identity of the new molecule created.  The output
  ! follows the convention in subinteract_update, species numbers are 
  ! listed as positive if they are inserted, negative if they are deleted
  ! and as both if they are perturbed.
  ! Requires:  params  -- parameters for this move
  !            spcs  -- returned, changed species numbers
  !-----------------------------------------------------------------------------
  Subroutine idlibrary_spcs(params,spcs)
    Type(IDLibrary_Params), Intent(In)     :: params
    Integer, Dimension(:), Intent(Out)     :: spcs

    spcs(1) = -params%osubset(1)
    spcs(2) = params%nsubset(1)

  End Subroutine idlibrary_spcs

  !----------------------------------------------------------------------------
  ! Writes a sample section of the control file information 
  ! Requires:  unitno
  !----------------------------------------------------------------------------
  Subroutine idlibrary_sampleCF(unitno)
    Integer, Intent(In)           :: unitno
    
    Write(unitno,'(a,t30,a)') 'IDLIBRARY','# type of move'
    Write(unitno,'(a,t30,a)') 'Integer List', &
        '# species numbers into which to transform'

  End Subroutine idlibrary_sampleCF

  !----------------------------------------------------------------------
  ! Single line display information for dynamic Translation parameters
  ! Requires: params -- random translation parameters
  !----------------------------------------------------------------------
  Function idlibrary_dyninfo(params)
    Character(len=strLen)               :: idlibrary_dyninfo
    Type(IDLibrary_Params), Intent(In)  :: params

    idlibrary_dyninfo = ''

  End Function idlibrary_dyninfo

  !----------------------------------------------------------------------
  ! Displays the library move parameters
  ! Requires:  params -- identity library parameters
  !            indent -- indentation from left margin
  !            unit -- unit to dump into
  !----------------------------------------------------------------------
  Subroutine idlibrary_display(params,indent,unit)
    Type(IDLibrary_Params), Intent(In)  :: params
    Integer, Intent(In)                 :: indent,unit

    Integer                     :: i
    Character(len=indent)       :: blank
    Character(len=strLen)       :: string1,string2

    blank = Repeat(' ',indent)

    Write(unit,'(2a)') blank,"Identity Library Move Parameters :"
    string1 = int2str(params%spc)
    Write(unit,'(4a,10i2)') blank,'Species ',Trim(string1),' -> ', &
        (params%tospc(i),i=1,params%ntospc)

  End Subroutine idlibrary_display

  !----------------------------------------------------------------------
  ! Clean the library move parameters
  ! Requires:  params -- identity library parameters
  !----------------------------------------------------------------------
  Subroutine idlibrary_clean(params)
    Type(IDLibrary_Params), Intent(InOut)   :: params

    Integer            :: error

    Deallocate(params%tospc, STAT=error)
    If (error/=0) Call allocErrDisplay(__FILE__,__LINE__)

    Call config_clean(params%molcopy)

  End Subroutine idlibrary_clean
    
End Module idlibrary
